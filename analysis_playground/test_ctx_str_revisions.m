%% Playground for assorted paper revision code tests

%% Load in data

% New datasets

% data_fn = 'trial_activity_choiceworld_muafilt';
% data_fn = 'trial_activity_choiceworld_3depth_muafilt';

% data_fn = 'trial_activity_AP_choiceWorldStimPassive_trained_muafilt';
% data_fn = 'trial_activity_AP_choiceWorldStimPassive_naive_muafilt';

% data_fn = 'trial_activity_AP_choiceWorldStimPassive_trained_3depth_muafilt';
% data_fn = 'trial_activity_AP_choiceWorldStimPassive_naive_3depth_muafilt';

% data_fn = 'trial_activity_AP_lcrGratingPassive_ctxstrephys_str';
% data_fn = 'trial_activity_AP_lcrGratingPassive_ctxstrephys_ctx';

% data_fn = 'trial_activity_AP_lcrGratingPassive_pre_muscimol';

AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_recording;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));




%% Testing: "baseline" cortical prediction: remove cortical regions

% Set example experiment to use
animal = 'AP025';
day = '2017-10-04';
experiment = 1;
AP_load_experiment;


% Align U, recast V's excluding high kernel weight regions
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
Udf_aligned = AP_align_widefield(Udf,animal,day);

kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
load(kernel_roi_fn);

zero_region_idx = 1;
wf_region_zero = kernel_roi.bw(:,:,zero_region_idx) + ...
    AP_reflect_widefield(kernel_roi.bw(:,:,zero_region_idx));
Udf_aligned_regionzero = Udf_aligned.*~wf_region_zero;
fVdf_regionzero_altU = ChangeU(Udf_aligned,fVdf,Udf_aligned_regionzero);

% (those U's aren't orthonormal, recast back to original Udf_aligned)
fVdf_regionzero = ChangeU(Udf_aligned_regionzero,fVdf_regionzero_altU,Udf_aligned);



fVdf_deconv = AP_deconv_wf(fVdf_regionzero);






% Parameters for regression
regression_params.use_svs = 1:100;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

% Get time points to bin
sample_rate = framerate*regression_params.upsample_factor;
time_bins = frame_t(find(frame_t > ...
    regression_params.skip_seconds,1)):1/sample_rate: ...
    frame_t(find(frame_t-frame_t(end) < ...
    -regression_params.skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Resample deconvolved fluorescence
fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';

% Bin spikes to match widefield frames
binned_spikes = nan(n_depths,length(time_bin_centers));
for curr_depth = 1:n_depths
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    
    % (skip if no spikes at this depth)
    if isempty(curr_spike_times)
        continue
    end
    
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
end
binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);

% Load lambda from previously estimated and saved
lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
load(lambda_fn);
curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
if any(curr_animal_idx)
    curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
    if any(curr_day_idx)
        lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
    end
end

kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
    round(regression_params.kernel_t(2)*sample_rate);

% Regress cortex to striatum (region excluded)
[ctx_str_k,ctxpred_spikes_std,explained_var] = ...
    AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
    binned_spikes_std,kernel_frames,lambda, ...
    regression_params.zs,regression_params.cvfold, ...
    true,regression_params.use_constant);


% Reshape kernel and convert to pixel space
k_px = zeros(size(Udf_aligned_regionzero,1),size(Udf_aligned_regionzero,2), ...
    size(ctx_str_k{1},2),size(ctx_str_k{1},3));
for curr_spikes = 1:size(ctx_str_k{1},3)
    k_px(:,:,:,curr_spikes) = ...
        svdFrameReconstruct(Udf_aligned_regionzero(:,:,regression_params.use_svs), ...
        ctx_str_k{1}(:,:,curr_spikes));
end





%% Predict striatum WITHOUT input ctx regions (QUICK DIRTY VERSION)
% I doubt this is even usable because it's not comparable to doing it on
% each experiment basis so this whole thing is just a goddamn proof of
% concept



% (load in a dataset first)

% Load kernel templates
n_aligned_depths = 4;
kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template_' num2str(n_aligned_depths) '_depths.mat'];
load(kernel_template_fn);

% Get regions to zero based on kernel template weight
frac_max_weight = 0.1; % zero pixels > max weight * this
min_px = 1000; % get rid of small islands
dilate_px = 30; % pixels to dilate the zeroed regions
wf_region_zero_depths = false(size(kernel_template));
for curr_depth = 1:n_depths
    wf_region_zero_depths(:,:,curr_depth) = ...
        imdilate(bwareaopen(kernel_template(:,:,curr_depth) > ...
        max(reshape(kernel_template(:,:,curr_depth), ...
        [],1),[],1)*frac_max_weight,min_px),ones(dilate_px+1));
end

% Regress kernel ROI activity to striatum domain activity (per recording)
regression_params.use_svs = 1:100;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;
lambda = 10;
kernel_frames = floor(regression_params.kernel_t(1)*sample_rate): ...
    ceil(regression_params.kernel_t(2)*sample_rate);

mua_allcat_exp = mat2cell(mua_allcat,trials_recording,length(t),n_depths);
fluor_allcat_deconv_exp = mat2cell(fluor_allcat_deconv,use_split,length(t),n_vs);

mua_ctxtrialpred_exp = cellfun(@(x) nan(size(x)),mua_allcat_exp,'uni',false);
mua_ctxtrialpred_k = nan(length(regression_params.use_svs),length(kernel_frames),n_depths,length(use_split));
mua_ctxtrialpred_regionzero_exp = cellfun(@(x) nan(size(x)),mua_allcat_exp,'uni',false);
mua_ctxtrialpred_regionzero_k = nan(length(regression_params.use_svs),length(kernel_frames),n_depths,length(use_split));

for curr_exp = 1:length(trials_recording)
    for curr_depth = 1:n_depths
        
        curr_mua = reshape(mua_allcat_exp{curr_exp}(:,:,curr_depth)',1,[]);
        curr_fluor = reshape(permute( ...
            fluor_allcat_deconv_exp{curr_exp},[2,1,3]),[],n_vs)';
        
        % Zero out only current depth region
        U_master_regionzero = U_master.*~wf_region_zero_depths(:,:,curr_depth);
        fVdf_regionzero_altU = ChangeU(U_master(:,:,1:n_vs),curr_fluor,U_master_regionzero(:,:,1:n_vs));
        % (those U's aren't orthonormal, recast back to original Udf_aligned)
        fVdf_regionzero = ChangeU(U_master_regionzero(:,:,1:n_vs),fVdf_regionzero_altU,U_master(:,:,1:n_vs));
        curr_fluor_regionzero = fVdf_regionzero;
        
        % Skip if no data
        if all(isnan(curr_mua))
            continue
        end
        
        % Set discontinuities in trial data
        trial_discontinuities = false(size(mua_allcat_exp{curr_exp}(:,:,curr_depth)));
        trial_discontinuities(:,1) = true;
        trial_discontinuities = reshape(trial_discontinuities',[],1)';
        
        % Full cortex regression
        [k_fluor,curr_mua_fluorpred,explained_var] = ...
            AP_regresskernel(curr_fluor(regression_params.use_svs,:), ...
            curr_mua,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant,trial_discontinuities);
        
        mua_ctxtrialpred_k(:,:,curr_depth,curr_exp) = k_fluor;
        mua_ctxtrialpred_exp{curr_exp}(:,:,curr_depth) = ...
            reshape(curr_mua_fluorpred,length(t),[])';
        
        % Region-zeroed cortex regression
        [k_fluorregionzero,curr_mua_fluorregionzeropred,explained_var] = ...
            AP_regresskernel(curr_fluor_regionzero(regression_params.use_svs,:), ...
            curr_mua,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant,trial_discontinuities);
        
        mua_ctxtrialpred_regionzero_k(:,:,curr_depth,curr_exp) = k_fluorregionzero;
        mua_ctxtrialpred_regionzero_exp{curr_exp}(:,:,curr_depth) = ...
            reshape(curr_mua_fluorregionzeropred,length(t),[])';
        
    end
    AP_print_progress_fraction(curr_exp,length(trials_recording));
end

% Get average cortex->striatum kernel
ctx_str_k_mean = nanmean(cell2mat(permute(vertcat(ctx_str_k_all{:}),[2,3,4,5,1])),5);
ctx_str_k_mean_px = cell2mat(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:100), ...
    ctx_str_k_mean(:,:,x)),permute(1:n_depths,[1,3,4,2]),'uni',false));

mua_ctxtrialpred_k_mean = nanmean(mua_ctxtrialpred_k,4);
mua_ctxtrialpred_k_mean_px = cell2mat(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:100), ...
    mua_ctxtrialpred_k_mean(:,:,x)),permute(1:n_depths,[1,3,4,2]),'uni',false));

mua_ctxtrialpred_regionzero_k_mean = nanmean(mua_ctxtrialpred_regionzero_k,4);
mua_ctxtrialpred_regionzero_k_mean_px = cell2mat(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:100), ...
    mua_ctxtrialpred_regionzero_k_mean(:,:,x)),permute(1:n_depths,[1,3,4,2]),'uni',false));

AP_image_scroll([ctx_str_k_mean_px,mua_ctxtrialpred_k_mean_px,mua_ctxtrialpred_regionzero_k_mean_px]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(gca,brewermap([],'*RdBu'));
axis image;

% Plot stim-aligned
plot_trials = trial_stim_allcat == 1;
plot_depth = 4;
figure; hold on;
plot(nanmean(mua_allcat(plot_trials,:,plot_depth)));
plot(nanmean(mua_ctxpred_allcat(plot_trials,:,plot_depth)));
plot(nanmean(AP_index_ans(vertcat(mua_ctxtrialpred_exp{:}),plot_trials,[],plot_depth)));
plot(nanmean(AP_index_ans(vertcat(mua_ctxtrialpred_regionzero_exp{:}),plot_trials,[],plot_depth)));
legend({'MUA','MUA ctxpred','MUA ctxtrialpred','MUA ctxtrialpred regionzero'});



% Get R^2 for task, cortex full, and cortex ROI predictions
mua_ctxtrialpred_allcat = cell2mat(mua_ctxtrialpred_exp);
mua_ctxtrialpred_regionzero_allcat = cell2mat(mua_ctxtrialpred_regionzero_exp);

ctxpred_r2 = nan(max(split_idx),n_depths);
ctxtrialpred_r2 = nan(max(split_idx),n_depths);
ctxtrialpred_regionzero_r2 = nan(max(split_idx),n_depths);

for curr_exp = 1:max(split_idx)
    
    curr_data = reshape(permute(mua_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_ctxpred_data = reshape(permute(mua_ctxpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_ctxtrialpred_data = reshape(permute(mua_ctxtrialpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_ctxtrialpred_regionzero_data = reshape(permute(mua_ctxtrialpred_regionzero_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_ctxpred_data) ...
        | isnan(curr_ctxtrialpred_data) | isnan(curr_ctxtrialpred_regionzero_data);
    curr_data(nan_samples) = NaN;
    curr_ctxpred_data(nan_samples) = NaN;
    curr_ctxtrialpred_data(nan_samples) = NaN;
    curr_ctxtrialpred_regionzero_data(nan_samples) = NaN;
    
    ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    ctxtrialpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxtrialpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    ctxtrialpred_regionzero_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxtrialpred_regionzero_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    
end
figure; hold on;
errorbar(nanmean(ctxpred_r2,1),AP_sem(ctxpred_r2,1),'linewidth',2,'CapSize',0);
errorbar(nanmean(ctxtrialpred_r2,1),AP_sem(ctxpred_r2,1),'linewidth',2,'CapSize',0);
errorbar(nanmean(ctxtrialpred_regionzero_r2,1),AP_sem(ctxpred_r2,1),'linewidth',2,'CapSize',0);

xlabel('Striatum depth');
ylabel('Explained variance');
legend({'Cortex (full)','Cortex (trial)','Cortex (trial,region-zero)'})

%% Predict striatum with zeroed-out cortical region

% (load in a dataset first)

% Choose depths to run
plot_depth = 1;

% Regress kernel ROI activity to striatum domain activity (per recording)
regression_params.use_svs = 1:200;
regression_params.kernel_t = [0,0];
regression_params.zs = [false,false];
regression_params.cvfold = 2;
regression_params.use_constant = false;
lambda = 50;
kernel_frames = floor(regression_params.kernel_t(1)*sample_rate): ...
    ceil(regression_params.kernel_t(2)*sample_rate);

% use_t = t < 0;
% use_t = t > 0 & t < 0.2;
% use_t = t > 0.5 & t < 1;
% use_t = t > 0.5;
% use_t = t > 1.5;
use_t = true(size(t));

mua_allcat_exp = mat2cell(mua_allcat,trials_recording,length(t),n_depths);
fluor_allcat_deconv_exp = mat2cell(fluor_allcat_deconv,use_split,length(t),n_vs);

mua_ctxtrialpred_exp = cellfun(@(x) nan(size(x)),mua_allcat_exp,'uni',false);
mua_ctxtrialpred_k = nan(length(regression_params.use_svs),length(kernel_frames),n_depths,length(use_split));
mua_ctxtrialpred_regionzero_exp = cellfun(@(x) nan(size(x)),mua_allcat_exp,'uni',false);
mua_ctxtrialpred_regionzero_k = nan(length(regression_params.use_svs),length(kernel_frames),n_depths,length(use_split));

for curr_exp = 1:length(trials_recording)
    for curr_depth = plot_depth
        
        curr_mua = reshape(mua_allcat_exp{curr_exp}(:,use_t,curr_depth)',1,[]);
        
        curr_fluor = reshape(permute( ...
            fluor_allcat_deconv_exp{curr_exp}(:,use_t,:),[2,1,3]),[],n_vs)';
        curr_fluor_full = reshape(permute( ...
            fluor_allcat_deconv_exp{curr_exp}(:,:,:),[2,1,3]),[],n_vs)';
        
        % Zero out region
        ctx_zero = false(size(U_master(:,:,1)));
        
        % (zero anterior)
%         ctx_zero(1:220,:) = true;
        % (zero posterior)
        %         ctx_zero(220:end,:) = true;
        % (zero left)
        %         ctx_zero(:,1:212) = true;
        % (zero right)
        %         ctx_zero(:,212:end) = true;
        % (zero bottom left)
                ctx_zero(220:end,1:212) = true;
        % (zero top left)
        %         ctx_zero(1:220,1:212) = true;
        
        U_master_regionzero = U_master.*~ctx_zero;
        fVdf_regionzero_altU = ChangeU(U_master(:,:,1:n_vs),curr_fluor,U_master_regionzero(:,:,1:n_vs));
        % (those U's aren't orthonormal, recast back to original Udf_aligned)
        fVdf_regionzero = ChangeU(U_master_regionzero(:,:,1:n_vs),fVdf_regionzero_altU,U_master(:,:,1:n_vs));
        curr_fluor_regionzero = fVdf_regionzero;
        
        % Skip if no data
        if all(isnan(curr_mua))
            continue
        end
        
        % Set discontinuities in trial data
        trial_discontinuities = false(size(mua_allcat_exp{curr_exp}(:,use_t,curr_depth)));
        trial_discontinuities(:,1) = true;
        trial_discontinuities = reshape(trial_discontinuities',[],1)';
        
        % Full cortex regression
        [k_fluor,curr_mua_fluorpred,explained_var_trial] = ...
            AP_regresskernel(curr_fluor(regression_params.use_svs,:), ...
            curr_mua,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant,trial_discontinuities);
        
        mua_ctxtrialpred_k(:,:,curr_depth,curr_exp) = k_fluor;
        %         mua_ctxtrialpred_exp{curr_exp}(:,use_t,curr_depth) = ...
        %             reshape(curr_mua_fluorpred,sum(use_t),[])';
        % (apply kernel to full time)
        mua_ctxtrialpred_exp{curr_exp}(:,:,curr_depth) = ...
            sum(cell2mat(arrayfun(@(x) ...
            convn(fluor_allcat_deconv_exp{curr_exp}(:,:,regression_params.use_svs(x)), ...
            k_fluor(x,:)','same'),permute(1:length(regression_params.use_svs),[1,3,2]),'uni',false)),3);
        
        % Region-zeroed cortex regression
        [k_fluorregionzero,curr_mua_fluorregionzeropred,explained_var_trial_regionzero] = ...
            AP_regresskernel(curr_fluor_regionzero(regression_params.use_svs,:), ...
            curr_mua,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant,trial_discontinuities);
        
        mua_ctxtrialpred_regionzero_k(:,:,curr_depth,curr_exp) = k_fluorregionzero;
        %         mua_ctxtrialpred_regionzero_exp{curr_exp}(:,use_t,curr_depth) = ...
        %             reshape(curr_mua_fluorregionzeropred,sum(use_t),[])';
        % (apply kernel to full time)
        mua_ctxtrialpred_regionzero_exp{curr_exp}(:,:,curr_depth) = ...
            sum(cell2mat(arrayfun(@(x) ...
            convn(fluor_allcat_deconv_exp{curr_exp}(:,:,regression_params.use_svs(x)), ...
            k_fluorregionzero(x,:)','same'),permute(1:length(regression_params.use_svs),[1,3,2]),'uni',false)),3);
        
    end
    AP_print_progress_fraction(curr_exp,length(trials_recording));
end

% Get average cortex->striatum kernel
mua_ctxtrialpred_k_mean = nanmean(mua_ctxtrialpred_k,4);
mua_ctxtrialpred_k_mean_px = cell2mat(arrayfun(@(x) ...
    svdFrameReconstruct(U_master(:,:,regression_params.use_svs), ...
    mua_ctxtrialpred_k_mean(:,:,x)),permute(1:n_depths,[1,3,4,2]),'uni',false));

mua_ctxtrialpred_regionzero_k_mean = nanmean(mua_ctxtrialpred_regionzero_k,4);
mua_ctxtrialpred_regionzero_k_mean_px = cell2mat(arrayfun(@(x) ...
    svdFrameReconstruct(U_master(:,:,regression_params.use_svs), ...
    mua_ctxtrialpred_regionzero_k_mean(:,:,x)),permute(1:n_depths,[1,3,4,2]),'uni',false));

% AP_image_scroll([ctx_str_k_mean_px,mua_ctxtrialpred_k_mean_px,mua_ctxtrialpred_regionzero_k_mean_px]);
AP_image_scroll([mua_ctxtrialpred_k_mean_px,mua_ctxtrialpred_regionzero_k_mean_px]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(gca,brewermap([],'*RdBu'));
axis image;

% Plot stim-aligned
figure; hold on;
for curr_depth = 1:length(plot_depth)
    subplot(length(plot_depth),2,(curr_depth-1)*2+1); hold on;
    plot_trials = trial_stim_allcat == 1;
    plot(t,nanmean(mua_allcat(plot_trials,:,plot_depth(curr_depth))));
    plot(t,nanmean(mua_ctxpred_allcat(plot_trials,:,plot_depth(curr_depth))));
    plot(t,nanmean(AP_index_ans(vertcat(mua_ctxtrialpred_exp{:}),plot_trials,[],plot_depth(curr_depth))));
    plot(t,nanmean(AP_index_ans(vertcat(mua_ctxtrialpred_regionzero_exp{:}),plot_trials,[],plot_depth(curr_depth))));
    title('Contra');
    
    subplot(length(plot_depth),2,(curr_depth-1)*2+2); hold on;
    plot_trials = trial_stim_allcat == -1;
    plot(t,nanmean(mua_allcat(plot_trials,:,plot_depth(curr_depth))));
    plot(t,nanmean(mua_ctxpred_allcat(plot_trials,:,plot_depth(curr_depth))));
    plot(t,nanmean(AP_index_ans(vertcat(mua_ctxtrialpred_exp{:}),plot_trials,[],plot_depth(curr_depth))));
    plot(t,nanmean(AP_index_ans(vertcat(mua_ctxtrialpred_regionzero_exp{:}),plot_trials,[],plot_depth(curr_depth))));
    title('Ipsi');
    legend({'MUA','MUA ctxpred','MUA ctxtrialpred','MUA ctxtrialpred regionzero'});
end
linkaxes(get(gcf,'Children'));

% Get R^2 for task, cortex full, and cortex ROI predictions
mua_ctxtrialpred_allcat = cell2mat(mua_ctxtrialpred_exp);
mua_ctxtrialpred_regionzero_allcat = cell2mat(mua_ctxtrialpred_regionzero_exp);

ctxpred_r2 = nan(max(split_idx),n_depths);
ctxtrialpred_r2 = nan(max(split_idx),n_depths);
ctxtrialpred_regionzero_r2 = nan(max(split_idx),n_depths);

for curr_exp = 1:max(split_idx)
    
    curr_data = reshape(permute(mua_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_ctxpred_data = reshape(permute(mua_ctxpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_ctxtrialpred_data = reshape(permute(mua_ctxtrialpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_ctxtrialpred_regionzero_data = reshape(permute(mua_ctxtrialpred_regionzero_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_ctxpred_data) ...
        | isnan(curr_ctxtrialpred_data) | isnan(curr_ctxtrialpred_regionzero_data);
    curr_data(nan_samples) = NaN;
    curr_ctxpred_data(nan_samples) = NaN;
    curr_ctxtrialpred_data(nan_samples) = NaN;
    curr_ctxtrialpred_regionzero_data(nan_samples) = NaN;
    
    ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    ctxtrialpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxtrialpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    ctxtrialpred_regionzero_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxtrialpred_regionzero_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    
end
figure; hold on;
errorbar(nanmean(ctxpred_r2,1),AP_sem(ctxpred_r2,1),'.-','linewidth',2,'CapSize',0,'MarkerSize',50);
errorbar(nanmean(ctxtrialpred_r2,1),AP_sem(ctxpred_r2,1),'.-','linewidth',2,'CapSize',0,'MarkerSize',50);
errorbar(nanmean(ctxtrialpred_regionzero_r2,1),AP_sem(ctxpred_r2,1),'.-','linewidth',2,'CapSize',0,'MarkerSize',50);

xlabel('Striatum depth');
ylabel('Explained variance');
legend({'Cortex (full)','Cortex (trial)','Cortex (trial,region-zero)'})



%% Apply kernel (from above) to whole dataset

plot_depth = 1;

use_k = asdf; %nanmean(mua_ctxtrialpred_k(:,:,plot_depth,:),4);

mua_kpred = ...
    sum(cell2mat(arrayfun(@(x) ...
    convn(fluor_allcat_deconv(:,:,x), ...
    use_k(x,:)','same'),permute(1:size(use_k,1),[1,3,2]),'uni',false)),3);


% Plot stim-aligned
figure; hold on;
subplot(1,2,1); hold on;
plot_trials = trial_stim_allcat == 1;
plot(t,nanmean(mua_allcat(plot_trials,:,plot_depth)));
plot(t,nanmean(mua_kpred(plot_trials,:)));
title('Contra');

subplot(1,2,2); hold on;
plot_trials = trial_stim_allcat == -1;
plot(t,nanmean(mua_allcat(plot_trials,:,plot_depth)));
plot(t,nanmean(mua_kpred(plot_trials,:)));
title('Ipsi');
legend({'MUA','MUA saved kernel'});

linkaxes(get(gcf,'Children'));











%% Predict striatum EACH TIME POINT (for passive)

% (load in a dataset first)
% (just do one depth at a time)
plot_depth = 1;

% Regress kernel ROI activity to striatum domain activity (per recording)
regression_params.use_svs = 1:100;
regression_params.kernel_t = [0,0];
regression_params.zs = [false,false];
regression_params.cvfold = 2;
regression_params.use_constant = false;
lambda = 20;
kernel_frames = floor(regression_params.kernel_t(1)*sample_rate): ...
    ceil(regression_params.kernel_t(2)*sample_rate);

mua_allcat_exp = mat2cell(mua_allcat,trials_recording,length(t),n_depths);
fluor_allcat_deconv_exp = mat2cell(fluor_allcat_deconv,use_split,length(t),n_vs);

mua_ctxtrialpred_exp = cellfun(@(x) nan([size(x),length(t)]),mua_allcat_exp,'uni',false);
mua_ctxtrialpred_k = nan(length(regression_params.use_svs), ...
    length(kernel_frames),n_depths,length(t),length(use_split));

for curr_exp = 1:length(trials_recording)
    for curr_depth = plot_depth
        for curr_t = 1:length(t)
            
            use_t = false(size(t));
            use_t(curr_t) = true;
            
            curr_mua = reshape(mua_allcat_exp{curr_exp}(:,use_t,curr_depth)',1,[]);
            curr_fluor = reshape(permute( ...
                fluor_allcat_deconv_exp{curr_exp}(:,use_t,:),[2,1,3]),[],n_vs)';
            
            % Skip if no data
            if all(isnan(curr_mua))
                continue
            end
            
            % Set discontinuities in trial data
            trial_discontinuities = false(size(mua_allcat_exp{curr_exp}(:,use_t,curr_depth)));
            trial_discontinuities(:,1) = true;
            trial_discontinuities = reshape(trial_discontinuities',[],1)';
            
            % Full cortex regression
            [k_fluor,curr_mua_fluorpred,explained_var_trial] = ...
                AP_regresskernel(curr_fluor(regression_params.use_svs,:), ...
                curr_mua,kernel_frames,lambda, ...
                regression_params.zs,regression_params.cvfold, ...
                false,regression_params.use_constant,trial_discontinuities);
            
            mua_ctxtrialpred_k(:,:,curr_depth,curr_t,curr_exp) = k_fluor;
            %         mua_ctxtrialpred_exp{curr_exp}(:,use_t,curr_depth) = ...
            %             reshape(curr_mua_fluorpred,sum(use_t),[])';
            % (apply kernel to full time)
            mua_ctxtrialpred_exp{curr_exp}(:,:,curr_depth,curr_t) = ...
                sum(cell2mat(arrayfun(@(x) ...
                convn(fluor_allcat_deconv_exp{curr_exp}(:,:,regression_params.use_svs(x)), ...
                k_fluor(x,:)','same'),permute(1:length(regression_params.use_svs),[1,3,2]),'uni',false)),3);
            
        end
    end
    AP_print_progress_fraction(curr_exp,length(trials_recording));
end

% Get average cortex->striatum kernel
ctx_str_k_mean = nanmean(cell2mat(permute(vertcat(ctx_str_k_all{:}),[2,3,4,5,1])),5);
ctx_str_k_mean_px = cell2mat(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:100), ...
    ctx_str_k_mean(:,:,x)),permute(1:n_depths,[1,3,4,2]),'uni',false));

mua_ctxtrialpred_k_mean = nanmean(mua_ctxtrialpred_k,5);
mua_ctxtrialpred_k_mean_px = cell2mat(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:100), ...
    reshape(mua_ctxtrialpred_k_mean(:,:,x,:),100,[])),permute(1:n_depths,[1,3,4,2]),'uni',false));

AP_image_scroll(mua_ctxtrialpred_k_mean_px);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(gca,brewermap([],'*RdBu'));
axis image;

% Plot stim-aligned
figure; hold on;
p1 = subplot(2,2,1); hold on;
set(gca,'ColorOrder',[copper(length(t));1,0,0;0,0.7,0]);
plot_trials = trial_stim_allcat == 1;
plot(t,squeeze(nanmean(AP_index_ans(vertcat(mua_ctxtrialpred_exp{:}),plot_trials,[],plot_depth,[]))));
plot(t,nanmean(mua_allcat(plot_trials,:,plot_depth)),'linewidth',2);
plot(t,nanmean(mua_ctxpred_allcat(plot_trials,:,plot_depth)),'linewidth',2);
title('Contra');

subplot(2,2,3);
imagesc(squeeze(nanmean(mua_allcat(plot_trials,:,plot_depth) - ...
    AP_index_ans(vertcat(mua_ctxtrialpred_exp{:}),plot_trials,[],plot_depth,[]))));
caxis([-0.5,0.5])
colormap(gca,brewermap([],'*RdBu'));
xlabel('Time trained');
ylabel('Time tested');

p2 = subplot(2,2,2); hold on;
set(gca,'ColorOrder',[copper(length(t));1,0,0;0,0.7,0]);
plot_trials = trial_stim_allcat == -1;
plot(t,squeeze(nanmean(AP_index_ans(vertcat(mua_ctxtrialpred_exp{:}),plot_trials,[],plot_depth,[]))));
plot(t,nanmean(mua_allcat(plot_trials,:,plot_depth)),'linewidth',2);
plot(t,nanmean(mua_ctxpred_allcat(plot_trials,:,plot_depth)),'linewidth',2);
title('Contra');
legend({'MUA','MUA ctxpred','MUA ctxtrialpred','MUA ctxtrialpred regionzero'});

subplot(2,2,4);
imagesc(squeeze(nanmean(mua_allcat(plot_trials,:,plot_depth) - ...
    AP_index_ans(vertcat(mua_ctxtrialpred_exp{:}),plot_trials,[],plot_depth,[]))));
caxis([-0.5,0.5])
colormap(gca,brewermap([],'*RdBu'));
xlabel('Time trained');
ylabel('Time tested');

linkaxes([p1,p2]);

%% Cortex > striatum kernels by depth


protocol = 'vanillaChoiceworld';

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
k_fn = [data_path filesep 'ctx_str_kernels_' protocol ,'_15strdepths'];
load(k_fn);

framerate = 35;
upsample_factor = 1;
sample_rate = framerate*upsample_factor;
kernel_t = [-0.1,0.1];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
t = kernel_frames/sample_rate;

% Concatenate explained variance
expl_var_animal = cell2mat(cellfun(@(x) nanmean(horzcat(x{:}),2),ctx_str_expl_var','uni',false));
figure('Name',protocol);
errorbar(nanmean(expl_var_animal,2),AP_sem(expl_var_animal,2),'k','linewidth',2);
xlabel('Striatal depth');
ylabel('Fraction explained variance');

% Concatenate and mean
% (kernel is -:+ fluorescence lag, flip to be spike-oriented)
k_px_timeflipped = cellfun(@(x) cellfun(@(x) x(:,:,end:-1:1,:),x,'uni',false),ctx_str_kernel,'uni',false);
k_px_animal = cellfun(@(x) nanmean(cat(5,x{:}),5),k_px_timeflipped,'uni',false);
k_px = nanmean(double(cat(5,k_px_animal{:})),5);

% Get center-of-mass maps
k_px_positive = k_px;
k_px_positive(k_px_positive < 0) = 0;
k_px_com = sum(k_px_positive.*permute(1:n_aligned_depths,[1,3,4,2]),4)./sum(k_px_positive,4);
k_px_com_colored = nan(size(k_px_com,1),size(k_px_com,2),3,size(k_px_com,3));

use_colormap = min(jet(255)-0.2,1);
for curr_frame = 1:size(k_px_com,3)
    k_px_com_colored(:,:,:,curr_frame) = ...
        ind2rgb(round(mat2gray(k_px_com(:,:,curr_frame),...
        [1,n_aligned_depths])*size(use_colormap,1)),use_colormap);
end

% Plot center kernel frames independently at t = 0
figure('Name',protocol);
plot_frame = kernel_frames == 0;
for curr_depth = 1:n_aligned_depths
    subplot(n_aligned_depths,1,curr_depth);
    imagesc(k_px(:,:,plot_frame,curr_depth));
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    axis image off;
    colormap(brewermap([],'*RdBu'));
    caxis([-0.01,0.01]);
end

% Plot center-of-mass color at select time points
plot_t = [-0.05:0.025:0.05];

k_px_com_colored_t = ...
    permute(reshape(interp1(t,permute(reshape(k_px_com_colored,[],3,length(t)), ...
    [3,1,2]),plot_t),length(plot_t),size(k_px_com_colored,1), ...
    size(k_px_com_colored,2),3),[2,3,4,1]);

k_px_max = squeeze(max(k_px,[],4));
k_px_max_t = ...
    permute(reshape(interp1(t,reshape(k_px_max,[],length(t))', ...
    plot_t),length(plot_t),size(k_px_max,1), ...
    size(k_px_max,2)),[2,3,1]);

weight_max = 0.005;
figure('Name',protocol);
for t_idx = 1:length(plot_t)
    subplot(1,length(plot_t),t_idx);
    p = image(k_px_com_colored_t(:,:,:,t_idx));
    set(p,'AlphaData', ...
        mat2gray(k_px_max_t(:,:,t_idx),[0,weight_max]));
    axis image off;
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title([num2str(plot_t(t_idx)),' s']);
end

% Plot movie of kernels
AP_image_scroll(reshape(permute(k_px,[1,4,2,3]),size(k_px,1)*size(k_px,4),size(k_px,2),length(t)),t);
colormap(brewermap([],'*RdBu'));
caxis([-max(caxis),max(caxis)]);
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5],[],[size(k_px,1),size(k_px,2),size(k_px,4),1]);
axis image off




% Plot STA for each depth    
ctx_str_kernel_cat = cell2mat(permute(horzcat(ctx_str_kernel{:}),[1,3,4,5,2]));
ctx_str_sta_cat = cell2mat(permute(horzcat(ctx_str_sta{:}),[1,3,4,2]));

ctx_str_kernel_cat_mean = squeeze(nanmean(ctx_str_kernel_cat(:,:,median(1:size(ctx_str_kernel_cat,3)),:,:),5));
ctx_str_sta_cat_mean = nanmean(ctx_str_sta_cat,4);

ctx_str_kernel_cat_mean_norm = ctx_str_kernel_cat_mean./max(max(ctx_str_kernel_cat_mean,[],1),[],2);
ctx_str_sta_cat_mean_norm = ctx_str_sta_cat_mean./max(max(ctx_str_sta_cat_mean,[],1),[],2);

figure;imagesc( ...
    [reshape(ctx_str_sta_cat_mean_norm,size(ctx_str_sta_cat_mean_norm,1),[]); ...
    reshape(ctx_str_kernel_cat_mean_norm,size(ctx_str_kernel_cat_mean_norm,1),[])]);
caxis([-1,1]);
colormap(brewermap([],'*RdBu'));


% Plot center-of-mass by color
k_px_com = sum(ctx_str_kernel_cat_mean_norm.* ...
    permute(1:size(ctx_str_kernel_cat_mean_norm,3),[1,3,2]),3)./ ...
    sum(ctx_str_kernel_cat_mean_norm,3);

use_colormap = min(jet(255)-0.2,1);
k_px_com_colored = ...
        ind2rgb(round(mat2gray(k_px_com,...
        [1,size(ctx_str_kernel_cat_mean_norm,3)])* ...
        size(use_colormap,1)),use_colormap);

weight_max = 1;
ctx_str_kernel_cat_mean_norm_max = max(ctx_str_kernel_cat_mean_norm,[],3);
figure;
image(k_px_com_colored,'AlphaData',mat2gray(ctx_str_kernel_cat_mean_norm_max,[0,weight_max]));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
axis image off





%% Predict striatum from cortex ephys
% This is garbage and I think explained variance is a garbage measure

% % Load data
% data_fn = 'trial_activity_vanillaChoiceworld_ctxstrephys_ctx';
% AP_load_concat_normalize_ctx_str;
% ctx_exp = vertcat(mua_all{:});
% 
% data_fn = 'trial_activity_vanillaChoiceworld_ctxstrephys_str';
% AP_load_concat_normalize_ctx_str;
% str_exp = vertcat(mua_all{:});

% Load data
data_fn = 'trial_activity_AP_lcrGratingPassive_ctxstrephys_ctx';
AP_load_concat_normalize_ctx_str;
ctx_exp = vertcat(mua_all{:});

data_fn = 'trial_activity_AP_lcrGratingPassive_ctxstrephys_str';
AP_load_concat_normalize_ctx_str;
str_exp = vertcat(mua_all{:});



% Choose depths to run
plot_depth = 1:3;

% Regress kernel ROI activity to striatum domain activity (per recording)
regression_params.kernel_t = [-0.1,0.1];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;
lambda = 0;
kernel_frames = floor(regression_params.kernel_t(1)*sample_rate): ...
    ceil(regression_params.kernel_t(2)*sample_rate);

use_t = true(size(t));

mua_ctxtrialpred_exp = cellfun(@(x) nan(size(x)),str_exp,'uni',false);
mua_ctxtrialpred_k = nan(1,length(kernel_frames),n_depths,length(str_exp));

for curr_exp = 1:length(str_exp)
    for curr_depth = plot_depth
        
        % (sum all cortical depths)
        curr_ctx = reshape(permute(sum(ctx_exp{curr_exp}(:,use_t,:),3),[2,1,3]),[],1)';
        curr_str = reshape(permute(str_exp{curr_exp}(:,use_t,:),[2,1,3]),[],n_depths)';
        
        % Skip if no data
        if any(all(isnan(curr_str),2))
            continue
        end
        
        % Set discontinuities in trial data
        trial_discontinuities = false(size(str_exp{curr_exp}(:,use_t,curr_depth)));
        trial_discontinuities(:,1) = true;
        trial_discontinuities = reshape(trial_discontinuities',[],1)';
        
        % Regress current domains from other domains
        [k_mua_alt,curr_mua_muaaltpred,explained_var_trial] = ...
            AP_regresskernel(curr_ctx, ...
            curr_str(curr_depth,:),kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant,trial_discontinuities);
        
        mua_ctxtrialpred_k(:,:,curr_depth,curr_exp) = k_mua_alt;
        mua_ctxtrialpred_exp{curr_exp}(:,use_t,curr_depth) = ...
            reshape(curr_mua_muaaltpred,sum(use_t),[])';
      
    end
    AP_print_progress_fraction(curr_exp,length(str_exp));
end


% Get R^2 for task, cortex full, and cortex ROI predictions
mua_ctxpred_exp = vertcat(mua_ctxpred_all{:});

ctxpred_r2 = nan(length(str_exp),n_depths);
ctxtrialpred_r2 = nan(length(str_exp),n_depths);

for curr_exp = 1:length(str_exp)
    
    curr_data = reshape(permute(str_exp{curr_exp},[2,1,3]),[],n_depths);
    curr_ctxpred_data = reshape(permute(mua_ctxpred_exp{curr_exp},[2,1,3]),[],n_depths);
    curr_ctxtrialpred_data = reshape(permute(mua_ctxtrialpred_exp{curr_exp},[2,1,3]),[],n_depths);
    
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_ctxpred_data) ...
        | isnan(curr_ctxtrialpred_data);
    curr_data(nan_samples) = NaN;
    curr_ctxpred_data(nan_samples) = NaN;
    curr_ctxtrialpred_data(nan_samples) = NaN;
    
    ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    ctxtrialpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxtrialpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    
end

figure;
subplot(1,2,1); hold on
errorbar(nanmean(ctxpred_r2,1),AP_sem(ctxpred_r2,1),'.-','linewidth',2,'CapSize',0,'MarkerSize',30);
errorbar(nanmean(ctxtrialpred_r2,1),AP_sem(ctxtrialpred_r2,1),'.-','linewidth',2,'CapSize',0,'MarkerSize',30);
legend({'Ctx','MUA alt'});
xlabel('Striatum depth');
ylabel('Explained variance');

mua_ctxtrialpred_k_mean = nanmean(mua_ctxtrialpred_k,4);
mua_col = lines(n_depths);
for curr_depth = 1:n_depths
    subplot(n_depths,2,curr_depth*2); hold on
    set(gca,'ColorOrder',mua_col(setdiff(1:n_depths,curr_depth),:));
    plot(kernel_frames,mua_ctxtrialpred_k_mean(:,:,curr_depth)','linewidth',2);
    ylabel('Weight');
    title(['Str ' num2str(curr_depth)]);
end



%% ~~~~~~~~~ THESE ARE USED - UPDATE AND TRANSFER

%% New corticostriatal kernel figure

% Load data
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\str_ctxpred.mat')

% Load Master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Get time
framerate = 35;
upsample_factor = 1;
sample_rate = framerate*upsample_factor;
kernel_t = [-0.1,0.1];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
t = kernel_frames/sample_rate;

% Concatenate kernels and convert to pixels
% (and flip in time so it's fluorescence lead:lag spikes)
ctx_str_k_animal = [str_ctxpred.ctx_str_k]';
ctx_str_k_px_animal = cellfun(@(x) cellfun(@(x) ...
    flip(AP_svdFrameReconstruct(U_master(:,:,1:100),x),3),x,'uni',false), ...
    ctx_str_k_animal,'uni',false);


% Get mean kernels and plot
ctx_str_k_px_cat = vertcat(ctx_str_k_px_animal{:}); 

ctx_str_k_px_task_mean = nanmean(cat(5,ctx_str_k_px_cat{:,1}),5);
ctx_str_k_px_notask_mean = nanmean(cat(5,ctx_str_k_px_cat{:,2}),5);

AP_image_scroll([ctx_str_k_px_task_mean,ctx_str_k_px_notask_mean]);
axis image;
colormap(brewermap([],'*RdBu'));
caxis([-max(abs(caxis)),max(abs(caxis))]);
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5],[],[size(U_master,1),size(U_master,2),1,2]);

figure;
for curr_depth = 1:n_depths
    subplot(n_depths,2,(curr_depth-1)*2+1);
    imagesc(ctx_str_k_px_task_mean(:,:,t == 0,curr_depth));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'*RdBu'));
    axis image off;
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title('Task');
    
    subplot(n_depths,2,(curr_depth-1)*2+2);
    imagesc(ctx_str_k_px_notask_mean(:,:,t == 0,curr_depth));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'*RdBu'));
    axis image off;
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title('Passive');
end

figure;
for curr_depth = 1:n_depths
    subplot(n_depths,1,curr_depth);
    imagesc(reshape(ctx_str_k_px_task_mean(:,:,:,curr_depth), ...
        size(ctx_str_k_px_task_mean,1),[]));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'*RdBu'));
    axis image off;
    title('Task');
end

figure;
for curr_depth = 1:n_depths
    subplot(n_depths,1,curr_depth);
    imagesc(reshape(ctx_str_k_px_notask_mean(:,:,:,curr_depth), ...
        size(ctx_str_k_px_notask_mean,1),[]));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'*RdBu'));
    axis image off;
    title('Passive');
end


% Get within- across-condition kernel correlation
n_depths = size(ctx_str_k_px_task_mean,4);

task_notask_k_corr = nan(4,n_depths,length(ctx_str_k_px_animal));
for curr_animal = 1:length(ctx_str_k_px_animal)
    
    curr_px = cellfun(@(x) reshape(x,[],n_depths), ...
        ctx_str_k_px_animal{curr_animal},'uni',false);
    
    curr_px_task = cat(3,curr_px{:,1});
    curr_px_notask = cat(3,curr_px{:,2});
     
    % Correlate kernel task/notask within domain
    task_notask_k_corr(1,:,curr_animal) = ...
        nanmean(cell2mat(cellfun(@(x,y) diag(corr(x,y))', ...
        curr_px(:,1),curr_px(:,2),'uni',false)));
    
    % Correlate kernel task across domains
    task_notask_k_corr(2,:,curr_animal) = ...
        nanmean(cell2mat(cellfun(@(x,y) ...
        nansum(tril(corr(x),-1)+triu(corr(x),1),1)./(n_depths-1)', ...
        curr_px(:,1),'uni',false)));
       
    % Correlate kernel within task/notask across days within domain
    task_notask_k_corr(3,:,curr_animal) = arrayfun(@(depth) ...
        nanmean(AP_itril(corr(permute(curr_px_task(:,depth,:),[1,3,2])),-1)),1:n_depths);
    task_notask_k_corr(4,:,curr_animal) = arrayfun(@(depth) ...
        nanmean(AP_itril(corr(permute(curr_px_notask(:,depth,:),[1,3,2])),-1)),1:n_depths);  

end

% Get mean across domains
task_notask_k_corr_strmean = squeeze(nanmean(task_notask_k_corr,2));

% Plot mean and split by domains
figure; 

subplot(2,1,1);hold on; set(gca,'ColorOrder',copper(n_depths));
plot(task_notask_k_corr_strmean,'color',[0.5,0.5,0.5]);
errorbar(nanmean(task_notask_k_corr_strmean,2), ...
    AP_sem(task_notask_k_corr_strmean,2),'k','linewidth',2);
set(gca,'XTick',1:4,'XTickLabelRotation',20,'XTickLabel', ...
    {'Task-no task within day','Task within day across domains','Task across days','No task across days'})
ylabel('Spatiotemporal correlation');
xlim([0.5,4.5]);

subplot(2,1,2);hold on; set(gca,'ColorOrder',copper(n_depths));
errorbar(nanmean(task_notask_k_corr,3), ...
    AP_sem(task_notask_k_corr,3),'linewidth',2)
set(gca,'XTick',1:4,'XTickLabelRotation',20,'XTickLabel', ...
    {'Task-no task within day','Task within day across domains','Task across days','No task across days'})
ylabel('Spatiotemporal correlation');
xlim([0.5,4.5]);
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false))

% (within task-passive v task-task domains statistics)
disp('Task/passive vs task/task cross-domain:')
curr_p = signrank(squeeze(task_notask_k_corr_strmean(1,:)), ...
    squeeze(task_notask_k_corr_strmean(2,:)));
disp(['All str ' num2str(curr_depth) ' p = ' num2str(curr_p)]);
for curr_depth = 1:n_depths  
    curr_p = signrank(squeeze(task_notask_k_corr(1,curr_depth,:)), ...
        squeeze(task_notask_k_corr(2,curr_depth,:)));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end

% (within vs across statistics)
disp('Task/passive-within vs task-across:')
curr_p = signrank(squeeze(task_notask_k_corr_strmean(1,:)), ...
    squeeze(task_notask_k_corr_strmean(3,:)));
disp(['All str ' num2str(curr_depth) ' p = ' num2str(curr_p)]);
for curr_depth = 1:n_depths
    curr_p = signrank(squeeze(task_notask_k_corr(1,curr_depth,:)), ...
        squeeze(task_notask_k_corr(3,curr_depth,:)));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end

% (cross task vs no task statistics)
disp('Task-across vs passive-across');
curr_p = signrank(squeeze(task_notask_k_corr_strmean(3,:)), ...
    squeeze(task_notask_k_corr_strmean(4,:)));
disp(['All str ' num2str(curr_depth) ' p = ' num2str(curr_p)]);
for curr_depth = 1:n_depths
    curr_p = signrank(squeeze(task_notask_k_corr(3,curr_depth,:)), ...
        squeeze(task_notask_k_corr(4,curr_depth,:)));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end


% Total weight sum in time
ctx_str_k_px_sum = cellfun(@(x) ...
    permute(sum(abs(reshape(x,[],size(x,3),size(x,4))),1),[2,3,1]), ...
    ctx_str_k_px_cat,'uni',false);

figure;

task_sum_norm_tdiff = nan(n_depths,size(ctx_str_k_px_sum,1));
for curr_depth = 1:n_depths
    subplot(n_depths,1,curr_depth); hold on;
    
    % Min/max-normalized sum of weights    
    curr_task_sum_norm = cell2mat(cellfun(@(x) ...
        mat2gray(x(:,curr_depth)), ...
        ctx_str_k_px_sum(:,1)','uni',false));
    curr_notask_sum_norm = cell2mat(cellfun(@(x) ...
        mat2gray(x(:,curr_depth)), ...
        ctx_str_k_px_sum(:,2)','uni',false));
    
    % Get weight difference t < 0 and t > 0
    task_sum_norm_tdiff(curr_depth,:) = ...
        nanmean(curr_task_sum_norm(t < 0,:) - ...
        flipud(curr_task_sum_norm(t > 0,:)),1);
    
    AP_errorfill(t,nanmean(curr_task_sum_norm,2),AP_sem(curr_task_sum_norm,2),'k');
    AP_errorfill(t,nanmean(curr_notask_sum_norm,2),AP_sem(curr_notask_sum_norm,2),'r');
    xlabel('Ctx lead:lag str (s)');
    ylabel('Sum(abs(W))');
    line([0,0],ylim,'color','k');
    
end

% (cross task vs no task statistics)
disp('Weight t<0 vs t>0:')
curr_p = signrank(nanmean(task_sum_norm_tdiff,1));
disp(['All str ' num2str(curr_depth) ' p = ' num2str(curr_p)]);

for curr_depth = 1:n_depths
    curr_p = signrank(task_sum_norm_tdiff(curr_depth,:));
    disp(['Str ' num2str(curr_depth) ' weight t <0 vs > 0 = ' num2str(curr_p)]); 
end



%%%% MAYBE DON'T DO THE BELOW
k_px = ctx_str_k_px_task_mean;
k_px = k_px./max(max(k_px,[],1),[],2);


% Plot center-of-mass color at select time points
k_px_com = sum(k_px.*permute(1:n_aligned_depths,[1,3,4,2]),4)./sum(k_px,4);
k_px_com_colored = nan(size(k_px_com,1),size(k_px_com,2),3,size(k_px_com,3));

use_colormap = min(jet(255)-0.2,1);
for curr_frame = 1:size(k_px_com,3)
    k_px_com_colored(:,:,:,curr_frame) = ...
        ind2rgb(round(mat2gray(k_px_com(:,:,curr_frame),...
        [1,n_aligned_depths])*size(use_colormap,1)),use_colormap);
end

plot_t = [-0.05:0.025:0.05];
k_px_com_colored_t = ...
    permute(reshape(interp1(t,permute(reshape(k_px_com_colored,[],3,length(t)), ...
    [3,1,2]),plot_t),length(plot_t),size(k_px_com_colored,1), ...
    size(k_px_com_colored,2),3),[2,3,4,1]);

k_px_max = squeeze(max(k_px_norm,[],4));
k_px_max_t = ...
    permute(reshape(interp1(t,reshape(k_px_max,[],length(t))', ...
    plot_t),length(plot_t),size(k_px_max,1), ...
    size(k_px_max,2)),[2,3,1]);

weight_max = 1;
figure;
for t_idx = 1:length(plot_t)
    subplot(1,length(plot_t),t_idx);
    p = image(k_px_com_colored_t(:,:,:,t_idx));
    set(p,'AlphaData', ...
        mat2gray(k_px_max_t(:,:,t_idx),[0,weight_max]));
    axis image off;
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title([num2str(plot_t(t_idx)),' s']);
end


%% Predict striatum with mutiple zeroed-out cortical regions

% (load in a dataset first)

% Choose depths to run
plot_depth = 1:n_depths;

% Regress kernel ROI activity to striatum domain activity (per recording)
regression_params.use_svs = 1:100;
regression_params.kernel_t = [-0.1,0.1];
regression_params.zs = [false,false];
regression_params.cvfold = 2;
regression_params.use_constant = true;
lambda = 50;
kernel_frames = floor(regression_params.kernel_t(1)*sample_rate): ...
    ceil(regression_params.kernel_t(2)*sample_rate);

% use_t = t < 0;
% use_t = t > 0.05 & t < 0.1;
% use_t = t > 0.5 & t < 1;
% use_t = t > 0.5;
use_t = true(size(t));

% Set regions to zero
% (zero all EXCEPT quadrants)
ctx_zero = true(size(U_master,1),size(U_master,2),6);
ctx_zero(1:260,1:220,1) = false;
ctx_zero(1:260,220:end,2) = false;
ctx_zero(260:end,1:220,3) = false;
ctx_zero(260:end,220:end,4) = false;
ctx_zero(:,220:end,5) = false;
ctx_zero(:,1:220,6) = false;

% (use raw data for trial regression)
mua_exp = vertcat(mua_all{:});
fluor_exp = vertcat(fluor_all{:});
fluor_deconv_exp = cellfun(@AP_deconv_wf,fluor_exp,'uni',false);

mua_ctxtrialpred_exp = cellfun(@(x) nan(size(x)),mua_exp,'uni',false);
mua_ctxtrialpred_k = nan(length(regression_params.use_svs),length(kernel_frames),n_depths,length(use_split));
mua_ctxtrialpred_regionzero_exp = cellfun(@(x) nan([size(x),size(ctx_zero,3)]),mua_exp,'uni',false);
mua_ctxtrialpred_regionzero_k = nan(length(regression_params.use_svs),length(kernel_frames),n_depths,length(use_split));

for curr_exp = 1:length(trials_recording)
    for curr_depth = plot_depth
        
        curr_mua = reshape(mua_exp{curr_exp}(:,use_t,curr_depth)',1,[]);
        curr_mua_std = curr_mua./nanstd(curr_mua);
        
        curr_fluor = reshape(permute( ...
            fluor_deconv_exp{curr_exp}(:,use_t,:),[2,1,3]),[],n_vs)';
        curr_fluor_full = reshape(permute( ...
            fluor_deconv_exp{curr_exp}(:,:,:),[2,1,3]),[],n_vs)';
        
        curr_fluor_regionzero = nan(size(curr_fluor,1),size(curr_fluor,2),size(ctx_zero,3));
        for curr_zero = 1:size(ctx_zero,3)
            U_master_regionzero = U_master.*~ctx_zero(:,:,curr_zero);
            fVdf_regionzero_altU = ChangeU(U_master(:,:,1:n_vs),curr_fluor,U_master_regionzero(:,:,1:n_vs));
            % (those U's aren't orthonormal, recast back to original Udf_aligned)
            fVdf_regionzero = ChangeU(U_master_regionzero(:,:,1:n_vs),fVdf_regionzero_altU,U_master(:,:,1:n_vs));
            curr_fluor_regionzero(:,:,curr_zero) = fVdf_regionzero;
        end
        
        % Skip if no data
        if all(isnan(curr_mua))
            continue
        end
        
        % Set discontinuities in trial data
        trial_discontinuities = false(size(mua_exp{curr_exp}(:,use_t,curr_depth)));
        trial_discontinuities(:,1) = true;
        trial_discontinuities = reshape(trial_discontinuities',[],1)';
        
        % Full cortex regression
        [ctx_str_k,curr_mua_fluorpred_std,explained_var_trial] = ...
            AP_regresskernel(curr_fluor(regression_params.use_svs,:), ...
            curr_mua_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant,trial_discontinuities);
                
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        ctxpred_spikes = (curr_mua_fluorpred_std - squeeze(ctx_str_k{end})).* ...
            nanstd(curr_mua,[],2) + ...
            nanstd(curr_mua,[],2).*squeeze(ctx_str_k{end});
        
        mua_ctxtrialpred_k(:,:,curr_depth,curr_exp) = ctx_str_k{1};
        mua_ctxtrialpred_exp{curr_exp}(:,use_t,curr_depth) = ...
            reshape(ctxpred_spikes,sum(use_t),[])';
        %         % (apply kernel to full time)
        %         mua_ctxtrialpred_exp{curr_exp}(:,:,curr_depth) = ...
        %             sum(cell2mat(arrayfun(@(x) ...
        %             convn(fluor_allcat_deconv_exp{curr_exp}(:,:,regression_params.use_svs(x)), ...
        %             k_fluor(x,:)','same'),permute(1:length(regression_params.use_svs),[1,3,2]),'uni',false)),3);
        
        % Region-zeroed cortex regression
        for curr_zero = 1:size(ctx_zero,3)
            [ctx_str_k_regionzero,curr_mua_fluorregionzeropred_std,explained_var_trial_regionzero] = ...
                AP_regresskernel(curr_fluor_regionzero(regression_params.use_svs,:,curr_zero), ...
                curr_mua_std,kernel_frames,lambda, ...
                regression_params.zs,regression_params.cvfold, ...
                true,regression_params.use_constant,trial_discontinuities);
            
            % Re-scale the prediction (subtract offset, multiply, add scaled offset)
            ctxpred_spikes_regionzero = (curr_mua_fluorregionzeropred_std - squeeze(ctx_str_k_regionzero{end})).* ...
                nanstd(curr_mua,[],2) + ...
                nanstd(curr_mua,[],2).*squeeze(ctx_str_k_regionzero{end});
                        
            mua_ctxtrialpred_regionzero_k(:,:,curr_depth,curr_exp,curr_zero) = ctx_str_k_regionzero{1};
            mua_ctxtrialpred_regionzero_exp{curr_exp}(:,use_t,curr_depth,curr_zero) = ...
                reshape(ctxpred_spikes_regionzero,sum(use_t),[])';
            %         % (apply kernel to full time)
            %         mua_ctxtrialpred_regionzero_exp{curr_exp}(:,:,curr_depth) = ...
            %             sum(cell2mat(arrayfun(@(x) ...
            %             convn(fluor_allcat_deconv_exp{curr_exp}(:,:,regression_params.use_svs(x)), ...
            %             k_fluorregionzero(x,:)','same'),permute(1:length(regression_params.use_svs),[1,3,2]),'uni',false)),3);
        end
        
    end
    AP_print_progress_fraction(curr_exp,length(trials_recording));
end

% Get average cortex->striatum kernel
ctx_str_k_mean = nanmean(cell2mat(permute(vertcat(ctx_str_k_all{:}),[2,3,4,5,1])),5);
ctx_str_k_mean_px = cell2mat(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:100), ...
    ctx_str_k_mean(:,:,x)),permute(1:n_depths,[1,3,4,2]),'uni',false));

mua_ctxtrialpred_k_mean = nanmean(mua_ctxtrialpred_k,4);
mua_ctxtrialpred_k_mean_px = cell2mat(arrayfun(@(x) ...
    svdFrameReconstruct(U_master(:,:,regression_params.use_svs), ...
    mua_ctxtrialpred_k_mean(:,:,x)),permute(1:n_depths,[1,3,4,2]),'uni',false));

mua_ctxtrialpred_regionzero_k_mean = nanmean(mua_ctxtrialpred_regionzero_k,4);
mua_ctxtrialpred_regionzero_k_mean_px = cell2mat(arrayfun(@(x) ...
    svdFrameReconstruct(U_master(:,:,regression_params.use_svs), ...
    reshape(mua_ctxtrialpred_regionzero_k_mean(:,:,x,:),length(regression_params.use_svs),[])), ...
    permute(1:n_depths,[1,3,4,2]),'uni',false));

AP_image_scroll(cat(3,mua_ctxtrialpred_k_mean_px,mua_ctxtrialpred_regionzero_k_mean_px));
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(gca,brewermap([],'*RdBu'));
axis image;

% Plot stim-aligned
% (doesn't make sense to do: would have to normalize the trial pred)
figure; hold on;
for curr_depth = 1:length(plot_depth)
    subplot(length(plot_depth),2,(curr_depth-1)*2+1); hold on;
    plot_trials = trial_stim_allcat == 1;
    plot(t,nanmean(mua_allcat(plot_trials,:,plot_depth(curr_depth))));
    plot(t,nanmean(mua_ctxpred_allcat(plot_trials,:,plot_depth(curr_depth))));
    plot(t,nanmean(AP_index_ans(vertcat(mua_ctxtrialpred_exp{:}),plot_trials,[],plot_depth(curr_depth))));
    plot(t,squeeze(nanmean(AP_index_ans(vertcat(mua_ctxtrialpred_regionzero_exp{:}),plot_trials,[],plot_depth(curr_depth),[]))));
    title('Contra');
    
    subplot(length(plot_depth),2,(curr_depth-1)*2+2); hold on;
    plot_trials = trial_stim_allcat == -1;
    plot(t,nanmean(mua_allcat(plot_trials,:,plot_depth(curr_depth))));
    plot(t,nanmean(mua_ctxpred_allcat(plot_trials,:,plot_depth(curr_depth))));
    plot(t,nanmean(AP_index_ans(vertcat(mua_ctxtrialpred_exp{:}),plot_trials,[],plot_depth(curr_depth))));
    plot(t,squeeze(nanmean(AP_index_ans(vertcat(mua_ctxtrialpred_regionzero_exp{:}),plot_trials,[],plot_depth(curr_depth),[]))));
    title('Ipsi');
    legend({'MUA','MUA ctxpred','MUA ctxtrialpred','MUA ctxtrialpred regionzero'});
end
linkaxes(get(gcf,'Children'));

% Get R^2 for task, cortex full, and cortex ROI predictions
mua_ctxpred_exp = vertcat(mua_ctxpred_all{:});
mua_ctxtrialpred_allcat = cell2mat(mua_ctxtrialpred_exp);
mua_ctxtrialpred_regionzero_allcat = cell2mat(mua_ctxtrialpred_regionzero_exp);

ctxpred_r2 = nan(max(split_idx),n_depths);
ctxtrialpred_r2 = nan(max(split_idx),n_depths);
ctxtrialpred_regionzero_r2 = nan(max(split_idx),n_depths,size(ctx_zero,3));

for curr_exp = 1:max(split_idx)
    
    curr_data = reshape(permute(mua_exp{curr_exp},[2,1,3]),[],n_depths);
    curr_ctxpred_data = reshape(permute(mua_ctxpred_exp{curr_exp},[2,1,3]),[],n_depths);
    curr_ctxtrialpred_data = reshape(permute(mua_ctxtrialpred_exp{curr_exp},[2,1,3]),[],n_depths);
    curr_ctxtrialpred_regionzero_data = ...
        reshape(permute(mua_ctxtrialpred_regionzero_exp{curr_exp},[2,1,3,4]),[],n_depths,size(ctx_zero,3));
    
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_ctxpred_data) ...
        | isnan(curr_ctxtrialpred_data) | any(isnan(curr_ctxtrialpred_regionzero_data),3);
    curr_data(nan_samples) = NaN;
    curr_ctxpred_data(nan_samples) = NaN;
    curr_ctxtrialpred_data(nan_samples) = NaN;
    curr_ctxtrialpred_regionzero_data(repmat(nan_samples,1,1,size(ctx_zero,3))) = NaN;
    
    ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    ctxtrialpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxtrialpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    ctxtrialpred_regionzero_r2(curr_exp,:,:) = 1 - (nansum((curr_data-curr_ctxtrialpred_regionzero_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    
end

figure;
subplot(2,2,1); hold on;
errorbar(nanmean(ctxpred_r2,1),AP_sem(ctxpred_r2,1),'.-','linewidth',2,'CapSize',0,'MarkerSize',30);
errorbar(nanmean(ctxtrialpred_r2,1),AP_sem(ctxtrialpred_r2,1),'.-','linewidth',2,'CapSize',0,'MarkerSize',30);
legend({'Ctx full','Ctx trial'});
xlabel('Striatum depth');
ylabel('Explained variance');
subplot(2,2,2);
errorbar(permute(nanmean(ctxtrialpred_regionzero_r2,1),[3,2,1]), ...
    permute(AP_sem(ctxtrialpred_regionzero_r2,1),[3,2,1]),'.-','linewidth',2,'CapSize',0,'MarkerSize',30);
xlabel('Region zeroed');
ylabel('Explained variance');
legend({'Str 1','Str 2','Str 3','Str 4'})
linkaxes(get(gcf,'Children'),'y');
for i = 1:size(ctx_zero,3)
   subplot(2,size(ctx_zero,3),size(ctx_zero,3)+i);
   imagesc(~ctx_zero(:,:,i));
   AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
   axis image off;
   colormap(gray);
end


%% Goodness-of-fit update


% Regress kernel ROI activity to striatum domain activity (per recording)
regression_params.kernel_t = [0,0];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;
lambda = 0;
kernel_frames = floor(regression_params.kernel_t(1)*sample_rate): ...
    ceil(regression_params.kernel_t(2)*sample_rate);

mua_ctxroi_k = arrayfun(@(x) nan(length(kernel_frames),n_depths),1:length(use_split),'uni',false);
mua_ctxroipred_allcat = nan(size(mua_allcat));
for curr_exp = 1:length(trials_recording)
    for curr_depth = 1:n_depths
        
        curr_mua = reshape(mua_allcat(split_idx == curr_exp,:,curr_depth)',[],1)';
        curr_fluor_kernelroi = reshape(fluor_kernelroi_deconv(split_idx == curr_exp,:,curr_depth)',[],1)';
        
        % Skip if no data
        if all(isnan(curr_mua))
            continue
        end
        
        % Set discontinuities in trial data
        trial_discontinuities = false(size(mua_allcat(split_idx == curr_exp,:,curr_depth)));
        trial_discontinuities(:,1) = true;
        trial_discontinuities = reshape(trial_discontinuities',[],1)';
        
        % Do regression
        [k,curr_mua_kernelroipred,explained_var] = ...
            AP_regresskernel(curr_fluor_kernelroi, ...
            curr_mua,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant,trial_discontinuities);
              
        mua_ctxroi_k{curr_exp}(:,curr_depth) = k;
        mua_ctxroipred_allcat(split_idx == curr_exp,:,curr_depth) = ...
            reshape(curr_mua_kernelroipred,length(t),[])';

    end
end






% Get R^2 for task, cortex full, and cortex ROI predictions
taskpred_r2 = nan(max(split_idx),n_depths);
ctxpred_r2 = nan(max(split_idx),n_depths);
ctxroipred_r2 = nan(max(split_idx),n_depths);
for curr_exp = 1:max(split_idx)
       
    curr_data = reshape(permute(mua_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_taskpred_data = reshape(permute(mua_taskpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_ctxpred_data = reshape(permute(mua_ctxpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_ctxroipred_data = reshape(permute(mua_ctxroipred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_taskpred_data) | isnan(curr_ctxpred_data) | isnan(curr_ctxroipred_data);
    curr_data(nan_samples) = NaN;
    curr_taskpred_data(nan_samples) = NaN;
    curr_ctxpred_data(nan_samples) = NaN;
    curr_ctxroipred_data(nan_samples) = NaN;

    taskpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_taskpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    ctxroipred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxroipred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
end
figure; hold on;
errorbar(nanmean(taskpred_r2,1),AP_sem(taskpred_r2,1),'b','linewidth',2,'CapSize',0);
errorbar(nanmean(ctxpred_r2,1),AP_sem(ctxpred_r2,1),'color',[0,0.6,0],'linewidth',2,'CapSize',0);
errorbar(nanmean(ctxroipred_r2,1),AP_sem(ctxroipred_r2,1),'color',[1,0.5,0],'linewidth',2,'CapSize',0);
xlabel('Striatum depth');
ylabel('Explained variance');
legend({'Task','Cortex (Full)','Cortex (ROI)'});

% Get significance between cortex kernel and ROI
ctx_kernel_roi_p = nan(n_depths,1);
for curr_depth = 1:n_depths
   ctx_kernel_roi_p(curr_depth) = signrank(ctxroipred_r2(:,curr_depth), ...
       ctxpred_r2(:,curr_depth));
   disp(['Str ' num2str(curr_depth) ' kernel vs ROI: p = ' ...
       num2str(ctx_kernel_roi_p(curr_depth))]);
end



% Get R^2 for task in cortex ROI
taskpred_r2 = nan(max(split_idx),n_depths);
for curr_exp = 1:max(split_idx)
       
    curr_data = reshape(permute(fluor_kernelroi_deconv(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_taskpred_data = reshape(permute(fluor_kernelroi_taskpred(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_taskpred_data);
    curr_data(nan_samples) = NaN;
    curr_taskpred_data(nan_samples) = NaN;

    taskpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_taskpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));

end
figure; hold on;
errorbar(nanmean(taskpred_r2,1),AP_sem(taskpred_r2,1),'b','linewidth',2,'CapSize',0);
xlabel('Cortex ROI');
ylabel('Task explained variance');


% Cortex explained variance

% (spatial explained variance in pixels)
px_taskpred_r2 = nan(size(U_master,1),size(U_master,2),max(split_idx));
for curr_exp = 1:max(split_idx)  
    px_taskpred_r2(:,:,curr_exp) = AP_spatial_explained_var(U_master(:,:,1:n_vs), ...
        reshape(permute(fluor_allcat_deconv(split_idx == curr_exp,:,:),[2,1,3]),[],n_vs)', ...
        reshape(permute(fluor_taskpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_vs)',10);
    AP_print_progress_fraction(curr_exp,max(split_idx));
end

figure;imagesc(nanmedian(px_taskpred_r2,3));
axis image off; 
colormap(brewermap([],'Reds'));
caxis([0,1]); 
c = colorbar;
ylabel(c,'Task R^2');
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);


%% Deconv explained variance in cortex ephys


data_fn = 'trial_activity_AP_lcrGratingPassive_ctxstrephys_ctx';
% data_fn = 'trial_activity_vanillaChoiceworld_ctxstrephys_ctx';

AP_load_concat_normalize_ctx_str;


if exist('mua_taskpred_all','var') 
    mua_exp = vertcat(mua_all{:});
    mua_ctxpred_exp = vertcat(mua_ctxpred_all{:});
    mua_taskpred_exp = vertcat(mua_taskpred_all{:});
    
    % MUA explained variance (widefield + task)
    taskpred_r2 = nan(length(mua_exp),n_depths);
    ctxpred_r2 = nan(length(mua_exp),n_depths);
    ctxroipred_r2 = nan(length(mua_exp),n_depths);
    for curr_exp = 1:length(mua_exp)
        
        curr_data = reshape(permute(mua_exp{curr_exp},[2,1,3]),[],n_depths);
        curr_data_baselinesub = reshape(permute(mua_exp{curr_exp},[2,1,3]),[],n_depths) - ...
            (nanmean(reshape(mua_exp{curr_exp}(:,t < 0,:),[],size(mua_exp{curr_exp},3)),1));
        curr_taskpred_data = reshape(permute(mua_taskpred_exp{curr_exp},[2,1,3]),[],n_depths);
        curr_ctxpred_data = reshape(permute(mua_ctxpred_exp{curr_exp},[2,1,3]),[],n_depths);
        
        % Set common NaNs
        nan_samples = isnan(curr_data) | isnan(curr_data_baselinesub) | ...
            isnan(curr_taskpred_data) | isnan(curr_ctxpred_data);
        curr_data(nan_samples) = NaN;
        curr_data_baselinesub(nan_samples) = NaN;
        curr_taskpred_data(nan_samples) = NaN;
        curr_ctxpred_data(nan_samples) = NaN;
        
        % (task regressed from average baseline-subtracted data)
        taskpred_r2(curr_exp,:) = 1 - (nansum((curr_data_baselinesub-curr_taskpred_data).^2,1)./ ...
            nansum((curr_data_baselinesub-nanmean(curr_data_baselinesub,1)).^2,1));
        % (cortex regressed from raw data)
        ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
            nansum((curr_data-nanmean(curr_data,1)).^2,1));
        
    end
    figure; hold on;
    errorbar(nanmean(taskpred_r2,1),AP_sem(taskpred_r2,1),'b','linewidth',2,'CapSize',0);
    errorbar(nanmean(ctxpred_r2,1),AP_sem(ctxpred_r2,1),'color',[0,0.6,0],'linewidth',2,'CapSize',0);
    xlabel('Cortex depth');
    ylabel('Task explained variance');  
end




mua_exp = vertcat(mua_all{:});
mua_ctxpred_exp = vertcat(mua_ctxpred_all{:});


% MUA explained variance (widefield only)
ctxpred_r2 = nan(length(mua_exp),n_depths);
for curr_exp = 1:length(mua_exp)
       
    curr_data = reshape(permute(mua_exp{curr_exp},[2,1,3]),[],n_depths);
    curr_ctxpred_data = reshape(permute(mua_ctxpred_exp{curr_exp},[2,1,3]),[],n_depths);
    
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_ctxpred_data);
    curr_data(nan_samples) = NaN;
    curr_ctxpred_data(nan_samples) = NaN;
   
    % (cortex regressed from raw data)
    ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    
end
figure; hold on;
errorbar(nanmean(ctxpred_r2,1),AP_sem(ctxpred_r2,1),'color',[0,0.6,0],'linewidth',2,'CapSize',0);
xlabel('Cortex depth');
ylabel('Explained variance');


%% Predict striatum from other domains
% (load in a dataset first)

% Choose depths to run
plot_depth = 1:n_depths;

% Regress kernel ROI activity to striatum domain activity (per recording)
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 2;
regression_params.use_constant = false;
lambda = 0;
kernel_frames = floor(regression_params.kernel_t(1)*sample_rate): ...
    ceil(regression_params.kernel_t(2)*sample_rate);

% use_t = t < 0;
% use_t = t > 0.05 & t < 0.1;
% use_t = t > 0.5 & t < 1;
% use_t = t > 0.5;
use_t = true(size(t));

mua_allcat_exp = mat2cell(mua_allcat,trials_recording,length(t),n_depths);

mua_ctxtrialpred_exp = cellfun(@(x) nan(size(x)),mua_allcat_exp,'uni',false);
mua_ctxtrialpred_k = nan(n_depths-1,length(kernel_frames),n_depths,length(use_split));

for curr_exp = 1:length(trials_recording)
    for curr_depth = plot_depth
        
        curr_mua = reshape(permute(mua_allcat_exp{curr_exp}(:,use_t,:),[2,1,3]),[],n_depths)';
        
        % Skip if no data
        if any(all(isnan(curr_mua),2))
            continue
        end
        
        % Set discontinuities in trial data
        trial_discontinuities = false(size(mua_allcat_exp{curr_exp}(:,use_t,curr_depth)));
        trial_discontinuities(:,1) = true;
        trial_discontinuities = reshape(trial_discontinuities',[],1)';
        
        % Regress current domains from other domains
        [k_mua_alt,curr_mua_muaaltpred,explained_var_trial] = ...
            AP_regresskernel(curr_mua(setdiff(1:n_depths,curr_depth),:), ...
            curr_mua(curr_depth,:),kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant,trial_discontinuities);
        
        mua_ctxtrialpred_k(:,:,curr_depth,curr_exp) = k_mua_alt;
        mua_ctxtrialpred_exp{curr_exp}(:,use_t,curr_depth) = ...
            reshape(curr_mua_muaaltpred,sum(use_t),[])';
      
    end
    AP_print_progress_fraction(curr_exp,length(trials_recording));
end


% Get R^2 for task, cortex full, and cortex ROI predictions
mua_ctxtrialpred_allcat = cell2mat(mua_ctxtrialpred_exp);

ctxpred_r2 = nan(max(split_idx),n_depths);
ctxtrialpred_r2 = nan(max(split_idx),n_depths);

for curr_exp = 1:max(split_idx)
    
    curr_data = reshape(permute(mua_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_ctxpred_data = reshape(permute(mua_ctxpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_ctxtrialpred_data = reshape(permute(mua_ctxtrialpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_ctxpred_data) ...
        | isnan(curr_ctxtrialpred_data);
    curr_data(nan_samples) = NaN;
    curr_ctxpred_data(nan_samples) = NaN;
    curr_ctxtrialpred_data(nan_samples) = NaN;
    
    ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    ctxtrialpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxtrialpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    
end

figure;
subplot(1,2,1); hold on
errorbar(nanmean(ctxpred_r2,1),AP_sem(ctxpred_r2,1),'.-','linewidth',2,'CapSize',0,'MarkerSize',30);
errorbar(nanmean(ctxtrialpred_r2,1),AP_sem(ctxtrialpred_r2,1),'.-','linewidth',2,'CapSize',0,'MarkerSize',30);
legend({'Ctx','MUA alt'});
xlabel('Striatum depth');
ylabel('Explained variance');

mua_ctxtrialpred_k_mean = nanmean(mua_ctxtrialpred_k,4);
mua_col = lines(n_depths);
for curr_depth = 1:n_depths
    subplot(n_depths,2,curr_depth*2); hold on
    set(gca,'ColorOrder',mua_col(setdiff(1:n_depths,curr_depth),:));
    plot(kernel_frames,mua_ctxtrialpred_k_mean(:,:,curr_depth)','linewidth',2);
    ylabel('Weight');
    title(['Str ' num2str(curr_depth)]);
end



%% Pre/post learning passive

% data_fns = { ...
%     'trial_activity_AP_choiceWorldStimPassive_naive', ...
%     'trial_activity_AP_choiceWorldStimPassive_trained'};

data_fns = { ...
    'trial_activity_AP_choiceWorldStimPassive_naive', ...
    {'trial_activity_AP_choiceWorldStimPassive_trained', ...
    'trial_activity_AP_lcrGratingPassive_ctxstrephys_str', ...
    'trial_activity_AP_lcrGratingPassive_pre_muscimol'}};

mua_prepost_norm = cell(2,1);
fluor_kernelroi_prepost_norm = cell(2,1);

stimIDs = cell(2,1);
mua_training = cell(2,1);
mua_ctxpred_training = cell(2,1);
fluor_training = cell(2,1);
fluor_roi_training = cell(2,1);
fluor_kernelroi_training = cell(2,1);

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;

    % Split by experiment
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    
    % Get trials with movement during stim to exclude
    wheel_thresh = 0.025;
    quiescent_trials = ~any(abs(wheel_allcat(:,t >= -0.5 & t <= 0.5)) > wheel_thresh,2);
    quiescent_trials_exp = mat2cell(quiescent_trials,trials_recording,1);
    
    % Get stim and activity by experiment
    trial_stim_allcat_exp = mat2cell(trial_stim_allcat,trials_recording,1);  
    fluor_deconv_exp = mat2cell(fluor_allcat_deconv,trials_recording,length(t),n_vs);
    fluor_roi_deconv_exp = mat2cell(fluor_roi_deconv,trials_recording,length(t),numel(wf_roi));
    fluor_kernelroi_deconv_exp = mat2cell(fluor_kernelroi_deconv,trials_recording,length(t),n_depths);
    mua_allcat_exp = mat2cell(mua_allcat,trials_recording,length(t),n_depths);
    mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,trials_recording,length(t),n_depths);
      
    warning('Pick normalization')    
%     % (TO RENORMALIZE: MUA = day baseline, fluor = nothing)  
%     mua_prepost_norm{curr_data} = vertcat(mua_norm{:});
%     mua_allcat_exp = cellfun(@(act,curr_norm,use_norm) ...
%         (act.*curr_norm)./use_norm,mua_allcat_exp, ...
%         vertcat(mua_norm{:}),vertcat(mua_day_baseline{:}),'uni',false);
%     mua_ctxpred_allcat_exp = cellfun(@(act,curr_norm,use_norm) ...
%         (act.*curr_norm)./use_norm,mua_ctxpred_allcat_exp, ...
%         vertcat(mua_norm{:}),vertcat(mua_day_baseline{:}),'uni',false);
%     
%     fluor_kernelroi_prepost_norm{curr_data} = fluor_kernelroi_norm;
%     fluor_kernelroi_deconv_exp = cellfun(@(act,curr_norm) ...
%         (act.*curr_norm),fluor_kernelroi_deconv_exp, ...
%         fluor_kernelroi_norm,'uni',false);
    
    % Exclude trials with fluorescence spikes
    % (this is a dirty way to do this but don't have a better alt)
    fluor_spike_thresh = 100;
    fluor_spike_trial = cellfun(@(x) any(any(x > fluor_spike_thresh,2),3), ...
        fluor_kernelroi_deconv_exp,'uni',false);
    
    % Grab data
    stimIDs{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        trial_stim_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    mua_training{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        mua_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    mua_ctxpred_training{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        mua_ctxpred_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_training{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_roi_training{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_roi_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_kernelroi_training{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_kernelroi_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
end

% Plot average fluorescence
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');

use_stim = 1;

fluor_pretrain_mean = ...
    svdFrameReconstruct(U_master(:,:,1:200), ...
    permute(nanmean(cell2mat(cellfun(@(x,stim) ...
    nanmean(x(stim == use_stim,:,:),1),fluor_training{1},stimIDs{1},'uni',false)),1),[3,2,1]));
fluor_posttrain_mean = ...
    svdFrameReconstruct(U_master(:,:,1:200), ...
    permute(nanmean(cell2mat(cellfun(@(x,stim) ...
    nanmean(x(stim == use_stim,:,:),1),fluor_training{2},stimIDs{2},'uni',false)),1),[3,2,1]));

AP_image_scroll([fluor_pretrain_mean,fluor_posttrain_mean,fluor_posttrain_mean-fluor_pretrain_mean]);
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5],[],[size(U_master,1),size(U_master,2),1,3]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
axis image;
colormap(brewermap([],'*RdBu'));


% Get average activity in relevant stim period and bin distribution
stim_avg_t = [0,0.2];
stim_avg_t_idx = t >= stim_avg_t(1) & t <= stim_avg_t(2);

stim_contra_act = cellfun(@(str,ctx,stim) cellfun(@(str,ctx,stim) ...
    [nanmean(str(stim == 1,stim_avg_t_idx,:),2), ...
    nanmean(ctx(stim == 1,stim_avg_t_idx,:),2)], ...
    str,ctx,stim,'uni',false), ...
    mua_training, fluor_kernelroi_training, stimIDs,'uni',false);

stim_ipsi_act = cellfun(@(str,ctx,stim) cellfun(@(str,ctx,stim) ...
    [nanmean(str(stim == -1,stim_avg_t_idx,:),2), ...
    nanmean(ctx(stim == -1,stim_avg_t_idx,:),2)], ...
    str,ctx,stim,'uni',false), ...
    mua_training, fluor_kernelroi_training, stimIDs,'uni',false);

mua_bins = linspace(-1,2,100);
mua_bin_centers = mua_bins(1:end-1) + diff(mua_bins)./2;
fluor_bins = linspace(-1,2,100);
fluor_bin_centers = fluor_bins(1:end-1) + diff(fluor_bins)./2;

stim_contra_dist = cellfun(@(x) ...
    cell2mat(permute(cellfun(@(x) cell2mat(arrayfun(@(str) ...
    histcounts2(x(:,1,str),x(:,2,str),mua_bins,fluor_bins), ...
    permute(1:n_depths,[1,3,2]),'uni',false)), ...
    x,'uni',false),[2,3,4,1])),stim_contra_act,'uni',false);

stim_ipsi_dist = cellfun(@(x) ...
    cell2mat(permute(cellfun(@(x) cell2mat(arrayfun(@(str) ...
    histcounts2(x(:,1,str),x(:,2,str),mua_bins,fluor_bins), ...
    permute(1:n_depths,[1,3,2]),'uni',false)), ...
    x,'uni',false),[2,3,4,1])),stim_ipsi_act,'uni',false);


% Plot average cortex + striatum responses and 2D distributions
use_stim = 1;

mua_mean = cellfun(@(act,stim) cell2mat(cellfun(@(act,stim) ...
    nanmean(act(stim == use_stim,:,:),1),act,stim,'uni',false)), ...
    mua_training,stimIDs,'uni',false);

fluor_kernelroi_mean = cellfun(@(act,stim) cell2mat(cellfun(@(act,stim) ...
    nanmean(act(stim == use_stim,:,:),1),act,stim,'uni',false)), ...
    fluor_kernelroi_training,stimIDs,'uni',false);

figure;
p = gobjects(n_depths,5);
for curr_str = 1:n_depths
    
    p(curr_str,1) = subplot(n_depths,5,(curr_str-1)*5+1);
    AP_errorfill(t,nanmean(fluor_kernelroi_mean{1}(:,:,curr_str),1)', ...
        AP_sem(fluor_kernelroi_mean{1}(:,:,curr_str),1)','k');
    AP_errorfill(t,nanmean(fluor_kernelroi_mean{2}(:,:,curr_str),1)', ...
        AP_sem(fluor_kernelroi_mean{2}(:,:,curr_str),1),'r');
    ylabel('Cortex (std)');
    line(repmat(stim_avg_t(1),2,1),ylim);
    line(repmat(stim_avg_t(2),2,1),ylim);
    
    p(curr_str,2) = subplot(n_depths,5,(curr_str-1)*5+2);
    AP_errorfill(t,nanmean(mua_mean{1}(:,:,curr_str),1)', ...
        AP_sem(mua_mean{1}(:,:,curr_str),1)','k');
    AP_errorfill(t,nanmean(mua_mean{2}(:,:,curr_str),1)', ...
        AP_sem(mua_mean{2}(:,:,curr_str),1)','r');
    xlabel('Time from stim (s)');
    ylabel('Striatum (std)');
    title(['Str ' num2str(curr_str)]);
    line(repmat(stim_avg_t(1),2,1),ylim);
    line(repmat(stim_avg_t(2),2,1),ylim);
    
    p(curr_str,3) = subplot(n_depths,5,(curr_str-1)*5+3);
    curr_contra_dist = imgaussfilt(squeeze(sum(stim_contra_dist{1}(:,:,curr_str,:),4)),3);
    curr_ipsi_dist = imgaussfilt(squeeze(sum(stim_ipsi_dist{1}(:,:,curr_str,:),4)),3);    
    imagesc(fluor_bin_centers,mua_bin_centers, ...
        curr_contra_dist./sum(curr_contra_dist(:)) - ...
        curr_ipsi_dist./sum(curr_ipsi_dist(:)));
    colormap(brewermap([],'*RdBu'));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    set(gca,'YDir','normal');
    
    p(curr_str,4) = subplot(n_depths,5,(curr_str-1)*5+4);
    curr_contra_dist = imgaussfilt(squeeze(sum(stim_contra_dist{2}(:,:,curr_str,:),4)),3);
    curr_ipsi_dist = imgaussfilt(squeeze(sum(stim_ipsi_dist{2}(:,:,curr_str,:),4)),3);    
    imagesc(fluor_bin_centers,mua_bin_centers, ...
        curr_contra_dist./sum(curr_contra_dist(:)) - ...
        curr_ipsi_dist./sum(curr_ipsi_dist(:)));
    colormap(brewermap([],'*RdBu'));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    set(gca,'YDir','normal');
    
    p(curr_str,5) = subplot(n_depths,5,(curr_str-1)*5+5);  hold on;
    
    curr_contra_untrained = cell2mat(cellfun(@(x) nanmean(x(:,:,curr_str,1),1),stim_contra_act{1},'uni',false));
    curr_ipsi_untrained = cell2mat(cellfun(@(x) nanmean(x(:,:,curr_str,1),1),stim_ipsi_act{1},'uni',false));
    curr_stim_untrained = permute(cat(3,curr_contra_untrained,curr_ipsi_untrained),[3,2,1]);
    
    curr_contra_trained = cell2mat(cellfun(@(x) nanmean(x(:,:,curr_str,1),1),stim_contra_act{2},'uni',false));
    curr_ipsi_trained = cell2mat(cellfun(@(x) nanmean(x(:,:,curr_str,1),1),stim_ipsi_act{2},'uni',false));
    curr_stim_trained = permute(cat(3,curr_contra_trained,curr_ipsi_trained),[3,2,1]);
        
    plot( ...
        cell2mat(cellfun(@(x) x(:,2,curr_str),stim_contra_act{1},'uni',false)), ...
        cell2mat(cellfun(@(x) x(:,1,curr_str),stim_contra_act{1},'uni',false)), ...
        '.','color',[0.5,0.5,0.5]);
    plot( ...
        cell2mat(cellfun(@(x) x(:,2,curr_str),stim_contra_act{2},'uni',false)), ...
        cell2mat(cellfun(@(x) x(:,1,curr_str),stim_contra_act{2},'uni',false)), ...
        '.','color',[1,0.5,0.5]);
    
    errorbar(nanmean(curr_stim_untrained(:,2,:),3),nanmean(curr_stim_untrained(:,1,:),3), ...
        AP_sem(curr_stim_untrained(:,1,:),3),AP_sem(curr_stim_untrained(:,1,:),3), ...
        AP_sem(curr_stim_untrained(:,2,:),3),AP_sem(curr_stim_untrained(:,2,:),3), ...
        'color','k','linewidth',2,'linestyle','none');
    scatter(nanmean(curr_stim_untrained(:,2,:),3),nanmean(curr_stim_untrained(:,1,:),3), ...
        100,[1,0,0;0,0,1],'filled','MarkerEdgeColor',[0.5,0.5,0.5]);
    
    errorbar(nanmean(curr_stim_trained(:,2,:),3),nanmean(curr_stim_trained(:,1,:),3), ...
        AP_sem(curr_stim_trained(:,1,:),3),AP_sem(curr_stim_trained(:,1,:),3), ...
        AP_sem(curr_stim_trained(:,2,:),3),AP_sem(curr_stim_trained(:,2,:),3), ...
        'color','k','linewidth',2,'linestyle','none');
    scatter(nanmean(curr_stim_trained(:,2,:),3),nanmean(curr_stim_trained(:,1,:),3), ...
        100,[1,0,0;0,0,1],'filled','MarkerEdgeColor',[0,0,0]);
    
    
end
linkaxes(p(:,1:2),'xy');














