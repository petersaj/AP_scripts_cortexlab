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
data_fn = 'trial_activity_AP_lcrGratingPassive_ctxstrephys_ctx';

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
ctx_str_k_mean = nanmean(cell2mat(permute(vertcat(ctx_str_k_all{:}),[2,3,4,5,1])),5);
ctx_str_k_mean_px = cell2mat(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:200), ...
    ctx_str_k_mean(:,:,x)),permute(1:n_depths,[1,3,4,2]),'uni',false));

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

%% Predict striatum with mutiple zeroed-out cortical regions

error('NEED TO FIX THIS: use non-normalized/baselined data like in fig');

% (load in a dataset first)

% Choose depths to run
plot_depth = 1:n_depths;

% Regress kernel ROI activity to striatum domain activity (per recording)
regression_params.use_svs = 1:200;
regression_params.kernel_t = [-0.1,0.1];
regression_params.zs = [false,false];
regression_params.cvfold = 2;
regression_params.use_constant = false;
lambda = 50;
kernel_frames = floor(regression_params.kernel_t(1)*sample_rate): ...
    ceil(regression_params.kernel_t(2)*sample_rate);

% use_t = t < 0;
% use_t = t > 0.05 & t < 0.1;
% use_t = t > 0.5 & t < 1;
% use_t = t > 0.5;
use_t = true(size(t));

% % Set regions to zero
% ctx_zero = false(size(U_master,1),size(U_master,2));
% % (zero anterior)
% ctx_zero(1:220,:,1) = true;
% % (zero posterior)
% ctx_zero(220:end,:,2) = true;
% % (zero left)
% ctx_zero(:,1:212,3) = true;
% % (zero right)
% ctx_zero(:,212:end,4) = true;
% % (zero bottom left)
% ctx_zero(220:end,1:212,5) = true;
% % (zero top left)
% ctx_zero(1:220,1:212,6) = true;

% (zero all EXCEPT quadrants)
ctx_zero = true(size(U_master,1),size(U_master,2),6);
ctx_zero(1:260,1:220,1) = false;
ctx_zero(1:260,220:end,2) = false;
ctx_zero(260:end,1:220,3) = false;
ctx_zero(260:end,220:end,4) = false;
ctx_zero(:,220:end,5) = false;
ctx_zero(:,1:220,6) = false;

mua_allcat_exp = mat2cell(mua_allcat,trials_recording,length(t),n_depths);
fluor_allcat_deconv_exp = mat2cell(fluor_allcat_deconv,use_split,length(t),n_vs);

mua_ctxtrialpred_exp = cellfun(@(x) nan(size(x)),mua_allcat_exp,'uni',false);
mua_ctxtrialpred_k = nan(length(regression_params.use_svs),length(kernel_frames),n_depths,length(use_split));
mua_ctxtrialpred_regionzero_exp = cellfun(@(x) nan([size(x),size(ctx_zero,3)]),mua_allcat_exp,'uni',false);
mua_ctxtrialpred_regionzero_k = nan(length(regression_params.use_svs),length(kernel_frames),n_depths,length(use_split));

for curr_exp = 1:length(trials_recording)
    for curr_depth = plot_depth
        
        curr_mua = reshape(mua_allcat_exp{curr_exp}(:,use_t,curr_depth)',1,[]);
        curr_mua_std = curr_mua./nanstd(curr_mua);
        
        curr_fluor = reshape(permute( ...
            fluor_allcat_deconv_exp{curr_exp}(:,use_t,:),[2,1,3]),[],n_vs)';
        curr_fluor_full = reshape(permute( ...
            fluor_allcat_deconv_exp{curr_exp}(:,:,:),[2,1,3]),[],n_vs)';
        
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
        trial_discontinuities = false(size(mua_allcat_exp{curr_exp}(:,use_t,curr_depth)));
        trial_discontinuities(:,1) = true;
        trial_discontinuities = reshape(trial_discontinuities',[],1)';
        
        % Full cortex regression
        [k_fluor,curr_mua_fluorpred_std,explained_var_trial] = ...
            AP_regresskernel(curr_fluor(regression_params.use_svs,:), ...
            curr_mua_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant,trial_discontinuities);
        curr_mua_fluorpred = curr_mua_fluorpred_std*nanstd(curr_mua);
        
        mua_ctxtrialpred_k(:,:,curr_depth,curr_exp) = k_fluor;
        mua_ctxtrialpred_exp{curr_exp}(:,use_t,curr_depth) = ...
            reshape(curr_mua_fluorpred,sum(use_t),[])';
        %         % (apply kernel to full time)
        %         mua_ctxtrialpred_exp{curr_exp}(:,:,curr_depth) = ...
        %             sum(cell2mat(arrayfun(@(x) ...
        %             convn(fluor_allcat_deconv_exp{curr_exp}(:,:,regression_params.use_svs(x)), ...
        %             k_fluor(x,:)','same'),permute(1:length(regression_params.use_svs),[1,3,2]),'uni',false)),3);
        
        % Region-zeroed cortex regression
        for curr_zero = 1:size(ctx_zero,3)
            [k_fluorregionzero,curr_mua_fluorregionzeropred_std,explained_var_trial_regionzero] = ...
                AP_regresskernel(curr_fluor_regionzero(regression_params.use_svs,:,curr_zero), ...
                curr_mua_std,kernel_frames,lambda, ...
                regression_params.zs,regression_params.cvfold, ...
                false,regression_params.use_constant,trial_discontinuities);
            curr_mua_fluorregionzeropred = curr_mua_fluorregionzeropred_std*nanstd(curr_mua);
            
            mua_ctxtrialpred_regionzero_k(:,:,curr_depth,curr_exp,curr_zero) = k_fluorregionzero;
            mua_ctxtrialpred_regionzero_exp{curr_exp}(:,use_t,curr_depth,curr_zero) = ...
                reshape(curr_mua_fluorregionzeropred,sum(use_t),[])';
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
ctx_str_k_mean_px = cell2mat(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:200), ...
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
mua_ctxtrialpred_allcat = cell2mat(mua_ctxtrialpred_exp);
mua_ctxtrialpred_regionzero_allcat = cell2mat(mua_ctxtrialpred_regionzero_exp);

ctxpred_r2 = nan(max(split_idx),n_depths);
ctxtrialpred_r2 = nan(max(split_idx),n_depths);
ctxtrialpred_regionzero_r2 = nan(max(split_idx),n_depths,size(ctx_zero,3));

for curr_exp = 1:max(split_idx)
    
    curr_data = reshape(permute(mua_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_ctxpred_data = reshape(permute(mua_ctxpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_ctxtrialpred_data = reshape(permute(mua_ctxtrialpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_ctxtrialpred_regionzero_data = reshape(permute(mua_ctxtrialpred_regionzero_allcat(split_idx == curr_exp,:,:,:),[2,1,3,4]),[],n_depths,size(ctx_zero,3));
    
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






%% Goodness-of-fit retry


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



%% Pre/post learning passive
% (currently using copied from muscimol so labelled as that atm)

data_fns = { ...
    'trial_activity_AP_choiceWorldStimPassive_naive', ...
    'trial_activity_AP_choiceWorldStimPassive_trained'};

mua_prepost_norm = cell(2,1);
fluor_kernelroi_prepost_norm = cell(2,1);

stimIDs = cell(2,1);
mua_muscimol = cell(2,1);
mua_ctxpred_muscimol = cell(2,1);
fluor_muscimol = cell(2,1);
fluor_roi_muscimol = cell(2,1);
fluor_kernelroi_muscimol = cell(2,1);

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
            
    % (RENORMALIZE: MUA = day baseline, fluor = nothing)  
    mua_prepost_norm{curr_data} = vertcat(mua_norm{:});
    mua_allcat_exp = cellfun(@(act,curr_norm,use_norm) ...
        (act.*curr_norm)./use_norm,mua_allcat_exp, ...
        vertcat(mua_norm{:}),vertcat(mua_day_baseline{:}),'uni',false);
    mua_ctxpred_allcat_exp = cellfun(@(act,curr_norm,use_norm) ...
        (act.*curr_norm)./use_norm,mua_ctxpred_allcat_exp, ...
        vertcat(mua_norm{:}),vertcat(mua_day_baseline{:}),'uni',false);
    
    fluor_kernelroi_prepost_norm{curr_data} = fluor_kernelroi_norm;
    fluor_kernelroi_deconv_exp = cellfun(@(act,curr_norm) ...
        (act.*curr_norm),fluor_kernelroi_deconv_exp, ...
        fluor_kernelroi_norm,'uni',false);
    
    % Exclude trials with fluorescence spikes
    % (this is a dirty way to do this but don't have a better alt)
    fluor_spike_thresh = 100;
    fluor_spike_trial = cellfun(@(x) any(any(x > fluor_spike_thresh,2),3), ...
        fluor_kernelroi_deconv_exp,'uni',false);
    
    % Grab data
    stimIDs{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        trial_stim_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    mua_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        mua_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    mua_ctxpred_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        mua_ctxpred_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_roi_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_roi_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_kernelroi_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_kernelroi_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
end


%%%%%%%%%%%%% TESTING

figure; hold on;

bar(reshape([nanmean(a(:,:,1),1);nanmean(b(:,:,1),1)],[],1));
errorbar(reshape([nanmean(a(:,:,1),1);nanmean(b(:,:,1),1)],[],1), ...
    reshape([AP_sem(a(:,:,1),1);AP_sem(b(:,:,1),1)],[],1),'.k','linewidth',3);
set(gca,'XTick',1:4,'XTickLabel',{'Str naive','Str trained','Ctx naive','Ctx trained'});

figure; hold on;
bar(reshape([nanmean(b(:,:,1),1)./nanmean(a(:,:,1),1)],[],1));


%%%%%%%%%%%%%%%%%%% CURRENTLY HERE DO THIS BY EXPT

stim_avg_t = [0,0.2];
stim_avg_t_idx = t >= stim_avg_t(1) & t <= stim_avg_t(2);
t_smooth_n = sum(t >= stim_avg_t(1) & t <= stim_avg_t(2));
t_smooth = ones(t_smooth_n,1)./t_smooth_n;

curr_data = 2;
stim_1_act = cellfun(@(str_act,ctx_act,stim) ...
    [nanmean(str_act(stim == 1,stim_avg_t_idx,:),2), ...
    nanmean(ctx_act(stim == 1,stim_avg_t_idx,:),2)], ...
    vertcat(mua_muscimol{curr_data}), ...
    vertcat(fluor_kernelroi_muscimol{curr_data}), ...
    vertcat(stimIDs{curr_data}),'uni',false);
stim_2_act = cellfun(@(str_act,ctx_act,stim) ...
    [nanmean(str_act(stim == -1,stim_avg_t_idx,:),2), ...
    nanmean(ctx_act(stim == -1,stim_avg_t_idx,:),2)], ...
    vertcat(mua_muscimol{curr_data}), ...
    vertcat(fluor_kernelroi_muscimol{curr_data}), ...
    vertcat(stimIDs{curr_data}),'uni',false);



mua_bins = linspace(-2,2,100);
fluor_bins = linspace(-0.1,0.1,100);

act_dist_1 = cellfun(@(x) ...
    histcounts2(x(:,1,1),x(:,2,1),mua_bins,fluor_bins), ...
    stim_1_act,'uni',false);

act_dist_2 = cellfun(@(x) ...
    histcounts2(x(:,1,1),x(:,2,1),mua_bins,fluor_bins), ...
    stim_2_act,'uni',false);



a = cat(3,act_dist_1{:});
b = cat(3,act_dist_2{:});




use_bin_centers = conv(use_bins,[0.5,0.5],'valid');
act_dist_1 = histcounts2(curr_act_pred_avg_cat(timeavg_trial_conditions{1}(:,1)), ...
    curr_act_avg_cat(timeavg_trial_conditions{1}(:,1)),use_bins,use_bins);
act_dist_2 = histcounts2(curr_act_pred_avg_cat(timeavg_trial_conditions{1}(:,2)), ...
    curr_act_avg_cat(timeavg_trial_conditions{1}(:,2)),use_bins,use_bins);
act_dist_iti = histcounts2( ...
    conv(reshape(curr_act_pred_allcat(:,t < 0 | t > 1.5,plot_areas),[],1),t_smooth,'same'), ...
    conv(reshape(curr_act_allcat(:,t < 0 | t > 1.5,plot_areas),[],1),t_smooth,'same'), ...
    use_bins,use_bins);

subplot(2,2,2); hold on;
imagesc(use_bin_centers,use_bin_centers,imgaussfilt(act_dist_1',3)-imgaussfilt(act_dist_2',3));
line(xlim,xlim)
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(gca,brewermap([],'*RdBu'));
title('Stim act distr difference');

subplot(2,2,3); hold on;
imagesc(use_bin_centers,use_bin_centers,act_dist_1');
line(xlim,xlim)
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(gca,brewermap([],'*RdBu'));
title('Contra distribution');

subplot(2,2,4); hold on;
imagesc(use_bin_centers,use_bin_centers,act_dist_2');
line(xlim,xlim)
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(gca,brewermap([],'*RdBu'));
title('Ipsi distribution');


%%%%%%%%%%%%%%%%%%%%%





% Plot average fluorescence
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');

use_stim = 1;

fluor_premuscimol_mean = ...
    svdFrameReconstruct(U_master(:,:,1:200), ...
    permute(nanmean(cell2mat(cellfun(@(x,stim) ...
    nanmean(x(stim == use_stim,:,:),1),fluor_muscimol{1},stimIDs{1},'uni',false)),1),[3,2,1]));
fluor_postmuscimol_mean = ...
    svdFrameReconstruct(U_master(:,:,1:200), ...
    permute(nanmean(cell2mat(cellfun(@(x,stim) ...
    nanmean(x(stim == use_stim,:,:),1),fluor_muscimol{2},stimIDs{2},'uni',false)),1),[3,2,1]));

AP_image_scroll([fluor_premuscimol_mean,fluor_postmuscimol_mean]);
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5],[],[size(U_master,1),size(U_master,2),1,2]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
axis image;
colormap(brewermap([],'*RdBu'));


figure;
t_stim = t >= 0.05 & t <= 0.15;

subplot(1,2,1)
imagesc(nanmean(fluor_premuscimol_mean(:,:,t_stim),3));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
c = caxis;
axis image off;
colormap(brewermap([],'*RdBu'));
title('Pre-muscimol');

subplot(1,2,2)
imagesc(nanmean(fluor_postmuscimol_mean(:,:,t_stim),3));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
caxis(c);
axis image off;
colormap(brewermap([],'*RdBu'));
title('Post-muscimol');

% Get pre/post stim response
use_stim = 1;

mua_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_muscimol{1},stimIDs{1},'uni',false));
mua_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_muscimol{2},stimIDs{2},'uni',false));

mua_ctxpred_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_ctxpred_muscimol{1},stimIDs{1},'uni',false));
mua_ctxpred_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_ctxpred_muscimol{2},stimIDs{2},'uni',false));

fluor_roi_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_roi_muscimol{1},stimIDs{1},'uni',false));
fluor_roi_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_roi_muscimol{2},stimIDs{2},'uni',false));

fluor_kernelroi_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_kernelroi_muscimol{1},stimIDs{1},'uni',false));
fluor_kernelroi_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_kernelroi_muscimol{2},stimIDs{2},'uni',false));

% Plot all cortex + striatum responses
figure;
p = gobjects(n_depths,2);
for curr_str = 1:n_depths
    
    p(curr_str,1) = subplot(n_depths,2,curr_str*2-1);
    AP_errorfill(t,nanmean(fluor_kernelroi_premuscimol_mean(:,:,curr_str),1)', ...
        AP_sem(fluor_kernelroi_premuscimol_mean(:,:,curr_str),1)','k');
    AP_errorfill(t,nanmean(fluor_kernelroi_postmuscimol_mean(:,:,curr_str),1)', ...
        AP_sem(fluor_kernelroi_postmuscimol_mean(:,:,curr_str),1),'r');
    ylabel('Cortex ROI');
    
    p(curr_str,2) = subplot(n_depths,2,curr_str*2);
    AP_errorfill(t,nanmean(mua_premuscimol_mean(:,:,curr_str),1)', ...
        AP_sem(mua_premuscimol_mean(:,:,curr_str),1)','k');
    AP_errorfill(t,nanmean(mua_postmuscimol_mean(:,:,curr_str),1)', ...
        AP_sem(mua_postmuscimol_mean(:,:,curr_str),1)','r');
    xlabel('Time from stim (s)');
    ylabel('Spikes (std)');
    title(['Str ' num2str(curr_str)]);
       
end
linkaxes(p(:,1),'xy');
linkaxes(p(:,2),'xy');




%% ~~~~~~~~~~ Trial activity tests

%% CTX EPHYS: passive trial activity WITH NEW MUA FILTER
clear all

animals = {'AP043','AP060','AP061'};
protocol = 'AP_lcrGratingPassive';
recording_site = {'ctxstrephys_str_muafilt','ctxstrephys_ctx_muafilt'};

% Load ephys alignment
vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
load(vis_ctx_ephys_fn);

for curr_site = 1:2
    
    % Initialize save variable
    trial_data_all = struct;
    
    for curr_animal = 1:length(animals)
        
        animal = animals{curr_animal};
        
        experiments = AP_find_experiments(animal,protocol);
        experiments = experiments([experiments.imaging] & [experiments.ephys]);
        
        disp(['Loading ' animal]);
        
        for curr_day = 1:length(experiments)
            
            day = experiments(curr_day).day;
            experiment = experiments(curr_day).experiment(1);
            
            % Load experiment
            site = curr_site;
            switch site
                case 1
                    % Site 1 = striautm
                    str_align = 'kernel';
                case 2
                    % Site 2 = cortex
                    str_align = 'none';
            end
            AP_load_experiment;
            
            % If cortical experiment - depth-align cortical MUA
            if site == 2
                % Get depth group boundaries
                % (these are just manual based on the CSD)
                n_aligned_depths = 2;
                ctx_depth_edges = linspace(0,1200,n_aligned_depths+1);
                
                % Depth-align templates
                curr_animal_idx = strcmp(animal,{vis_ctx_ephys.animal});
                curr_day_idx = strcmp(day,vis_ctx_ephys(curr_animal_idx).day);
                curr_csd_depth = vis_ctx_ephys(curr_animal_idx).stim_csd_depth{curr_day_idx};
                curr_csd_depth_aligned = vis_ctx_ephys(curr_animal_idx).stim_csd_depth_aligned{curr_day_idx};
                
                template_depths_aligned = interp1(curr_csd_depth,curr_csd_depth_aligned,template_depths);
                spike_depths_aligned = template_depths_aligned(spike_templates);
                
                % Find cortex end by largest gap between templates
                sorted_template_depths = sort([template_depths_aligned]);
                [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
                ctx_end = sorted_template_depths(max_gap_idx)+1;
                ctx_depth = [sorted_template_depths(1),ctx_end];
                ctx_units = template_depths_aligned <= ctx_depth(2);
                
                % Assign templates to depth groups
                % (NOTE: variable still called 'str' to make it easier to
                % fit into other code)
                ctx_spike_depths = spike_depths_aligned;
                ctx_spike_depths(spike_depths_aligned < ctx_depth(1) | spike_depths_aligned > ctx_depth(2)) = NaN;
                aligned_str_depth_group = discretize(ctx_spike_depths,ctx_depth_edges);
            end
            
            % Pull out trial data
            filter_mua = true; % Filter MUA to match widefield
            AP_ctx_str_grab_trial_data;
            
            % Store trial data into master structure
            trial_data_fieldnames = fieldnames(trial_data);
            for curr_trial_data_field = trial_data_fieldnames'
                trial_data_all.(cell2mat(curr_trial_data_field)){curr_animal,1}{curr_day,1} = ...
                    trial_data.(cell2mat(curr_trial_data_field));
            end
            
            % Store general info
            trial_data_all.animals = animals;
            trial_data_all.t = t;
            
            % Clear for next loop
            clearvars -except ...
                vis_ctx_ephys ...
                recording_site curr_site ...
                animals curr_animal animal protocol experiments curr_day ...
                trial_data_all
            
            AP_print_progress_fraction(curr_day,length(experiments));
        end
    end
    
    clearvars -except ...
        vis_ctx_ephys ...
        recording_site curr_site animals protocol ...
        trial_data_all
    
    disp('Finished loading all')
    
    save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
    save_fn = ['trial_activity_' protocol '_' recording_site{curr_site}];
    save([save_path filesep save_fn],'-v7.3');
    disp(['Saved ' save_fn]);
    
end



%% Choiceworld trial activity (striatum domain) - MUA filt + 3 depth

clear all
disp('Choiceworld trial activity (striatum domain)')

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        n_aligned_depths = 3;
        str_align = 'kernel';
        AP_load_experiment;
        
        % Pull out trial data
        filter_mua = true; % Filter MUA to match widefield
        AP_ctx_str_grab_trial_data;
        
        % Store trial data into master structure
        trial_data_fieldnames = fieldnames(trial_data);
        for curr_trial_data_field = trial_data_fieldnames'
            trial_data_all.(cell2mat(curr_trial_data_field)){curr_animal,1}{curr_day,1} = ...
                trial_data.(cell2mat(curr_trial_data_field));
        end
        
        % Store general info
        trial_data_all.animals = animals;
        trial_data_all.t = t;
        trial_data_all.task_regressor_labels = task_regressor_labels;
        trial_data_all.task_regressor_sample_shifts = task_regressor_sample_shifts;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars -except animals curr_animal animal protocol experiments curr_day ...
            trial_data_all
        
    end
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_filename = [save_path filesep 'trial_activity_choiceworld_3depth_muafilt'];
save(save_filename,'-v7.3');
disp(['Saved ' save_filename]);

%% Passive trial activity (3 depths + mua filt)

clear all

trained_animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
naive_animals = {'AP032','AP033','AP034','AP035','AP036'};
protocols = {'AP_choiceWorldStimPassive'};

for animal_group = {'trained','naive'}
    animal_group = cell2mat(animal_group);
    
    switch animal_group
        case 'trained'
            animals = trained_animals;
        case 'naive'
            animals = naive_animals;
    end
    
    for protocol = protocols
        protocol = cell2mat(protocol);
        
        disp(['Passive trial activity: ' animal_group ' ' protocol])
        
        % Initialize save variable
        trial_data_all = struct;
        
        for curr_animal = 1:length(animals)
            
            animal = animals{curr_animal};
            
            if strcmp(animal_group,'trained')
                % Only use days with choiceworld (sometimes recorded in cortex, no bhv)
                bhv_protocol = 'vanillaChoiceworld';
                behavior_experiments = AP_find_experiments(animal,bhv_protocol);
                passive_experiments = AP_find_experiments(animal,protocol);
                behavior_day = ismember({passive_experiments.day},{behavior_experiments.day});
                experiments = passive_experiments([passive_experiments.imaging] & [passive_experiments.ephys] & behavior_day);
            elseif strcmp(animal_group,'naive')
                experiments = AP_find_experiments(animal,protocol);
                experiments = experiments([experiments.imaging] & [experiments.ephys]);
            end
            
            disp(['Loading ' animal]);
            
            for curr_day = 1:length(experiments)
                
                day = experiments(curr_day).day;
                experiment = experiments(curr_day).experiment(end);
                
                % Load experiment
                n_aligned_depths = 3;
                str_align = 'kernel';
                AP_load_experiment;
                
                % Pull out trial data
                filter_mua = true; % Filter MUA to match widefield
                AP_ctx_str_grab_trial_data;
                
                % Store trial data into master structure
                trial_data_fieldnames = fieldnames(trial_data);
                for curr_trial_data_field = trial_data_fieldnames'
                    trial_data_all.(cell2mat(curr_trial_data_field)){curr_animal,1}{curr_day,1} = ...
                        trial_data.(cell2mat(curr_trial_data_field));
                end
                
                % Store general info
                trial_data_all.animals = animals;
                trial_data_all.t = t;
                
                AP_print_progress_fraction(curr_day,length(experiments));
            end
            
            % Clear for next loop
            clearvars -except ...
                trained_animals naive_animals protocols ...
                protocol animal_group animals curr_animal ...
                animals curr_animal animal protocol experiments curr_day ...
                trial_data_all
            
        end
        
        clearvars -except ...
            trained_animals naive_animals protocols ...
            protocol animal_group animals...
            trial_data_all
        disp('Finished loading all')
        
        save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
        save_fn = ['trial_activity_' protocol '_' animal_group '_3depth_muafilt'];
        save([save_path filesep save_fn],'-v7.3');
        disp(['Saved ' save_fn]);
        
    end
end


%% Striatum cortical kernels (3 depth, MUA filt, more SVs, less time)

disp('Cortex -> striatum regression maps across protocols');

n_aligned_depths = 3;

% Parameters for regression
regression_params.use_svs = 1:200;
regression_params.skip_seconds = 30;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.1,0.1];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

protocols = {'vanillaChoiceworld'};

for protocol = protocols
    protocol = cell2mat(protocol);
    
    animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
    
    init_array = cellfun(@(x) cell(0,0),animals','uni',false);
    ctx_str_kernel = init_array;
    ctx_str_expl_var = init_array;
    
    for curr_animal = 1:length(animals)
        
        animal = animals{curr_animal};
        disp(['Loading ' animal ' ' protocol]);
        
        % Only use days with choiceworld (sometimes recorded in cortex, no bhv)
        behavior_protocol = 'vanillaChoiceworld';
        behavior_experiments = AP_find_experiments(animal,behavior_protocol);
        
        curr_experiments = AP_find_experiments(animal,protocol);
        
        behavior_day = ismember({curr_experiments.day},{behavior_experiments.day});
        
        experiments = curr_experiments([curr_experiments.imaging] & [curr_experiments.ephys] & behavior_day);
        
        % Skip if this animal doesn't have this experiment
        if isempty(experiments)
            continue
        end
        
        disp(animal);
        
        load_parts.cam = false;
        load_parts.imaging = true;
        load_parts.ephys = true;
        
        for curr_day = 1:length(experiments)
            
            day = experiments(curr_day).day;
            experiment = experiments(curr_day).experiment(end);
            
            % Load data and align striatum by depth
            str_align = 'kernel';
            AP_load_experiment;
            
            %%% Load lambda from previously estimated and saved
            lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
            load(lambda_fn);
            curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
            if any(curr_animal_idx)
                curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
                if any(curr_day_idx)
                    lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
                end
            end
            lambda = 20;
            
            %%% Prepare data for regression
            
            % Get time points to bin
            sample_rate = framerate*regression_params.upsample_factor;
            time_bins = frame_t(find(frame_t > ...
                regression_params.skip_seconds,1)):1/sample_rate: ...
                frame_t(find(frame_t-frame_t(end) < ...
                -regression_params.skip_seconds,1,'last'));
            time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
            
            % Deconvolve fluorescence
            fVdf_deconv = AP_deconv_wf(fVdf);
            fVdf_deconv(isnan(fVdf_deconv)) = 0;
            fVdf_deconv_resample = interp1(frame_t,fVdf_deconv(regression_params.use_svs,:)',time_bin_centers)';
            
            % Get striatum depth group by across-experiment alignment
            n_depths = n_aligned_depths;
            depth_group = aligned_str_depth_group;
            
            binned_spikes = zeros(n_depths,length(time_bin_centers));
            for curr_depth = 1:n_depths
                curr_spike_times = spike_times_timeline(depth_group == curr_depth);
                binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
            end
            
            binned_spikes_filt = AP_deconv_wf(binned_spikes,true);
            binned_spikes_std = binned_spikes_filt./nanstd(binned_spikes_filt,[],2);
            
            %%% Regress MUA from cortex
            kernel_frames = floor(regression_params.kernel_t(1)*sample_rate): ...
                ceil(regression_params.kernel_t(2)*sample_rate);
            
            [k,ctxpred_spikes_std,explained_var] = ...
                AP_regresskernel(fVdf_deconv_resample, ...
                binned_spikes_std,kernel_frames,lambda, ...
                regression_params.zs,regression_params.cvfold, ...
                false,regression_params.use_constant);
            
            % Reshape kernel and convert to pixel space
            aUdf = AP_align_widefield(Udf,animal,day);
            k_px = zeros(size(aUdf,1),size(aUdf,2),size(k,2),size(k,3),'single');
            for curr_spikes = 1:size(k,3)
                k_px(:,:,:,curr_spikes) = svdFrameReconstruct(aUdf(:,:,regression_params.use_svs),k(:,:,curr_spikes));
            end
            
            % Store
            ctx_str_kernel{curr_animal}{curr_day} = k_px;
            ctx_str_expl_var{curr_animal}{curr_day} = explained_var.total;
            
            AP_print_progress_fraction(curr_day,length(experiments));
            clearvars -except regression_params n_aligned_depths ...
                animals animal curr_animal protocol ...
                experiments curr_day animal load_parts ...
                ctx_str_kernel ctx_str_expl_var
            
        end
        
        disp(['Finished ' animal]);
        
    end
    
    % Save
    save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data'];
    save([save_path filesep 'ctx_str_kernels_' protocol '_TESTING'],'-v7.3');
    warning('saving -v7.3');
    disp(['Finished ' protocol]);
    
end









