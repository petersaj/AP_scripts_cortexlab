%% Playground for assorted paper revision code tests

%% Load in data

% New datasets

data_fn = 'trial_activity_choiceworld_filtmua';
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

%% Predict striatum with zeroed-out cortical regions (for passive)

% (load in a dataset first)

% (just plot one depth)
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

% use_t = t < 0;
% use_t = t > 0.05 & t < 0.1;
% use_t = t > 0.5 & t < 1;
% use_t = t > 0.5;
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
ctx_str_k_mean_px = cell2mat(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:100), ...
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
subplot(1,2,1); hold on;
plot_trials = trial_stim_allcat == 1;
plot(t,nanmean(mua_allcat(plot_trials,:,plot_depth)));
plot(t,nanmean(mua_ctxpred_allcat(plot_trials,:,plot_depth)));
plot(t,nanmean(AP_index_ans(vertcat(mua_ctxtrialpred_exp{:}),plot_trials,[],plot_depth)));
plot(t,nanmean(AP_index_ans(vertcat(mua_ctxtrialpred_regionzero_exp{:}),plot_trials,[],plot_depth)));
title('Contra');

subplot(1,2,2); hold on;
plot_trials = trial_stim_allcat == -1;
plot(t,nanmean(mua_allcat(plot_trials,:,plot_depth)));
plot(t,nanmean(mua_ctxpred_allcat(plot_trials,:,plot_depth)));
plot(t,nanmean(AP_index_ans(vertcat(mua_ctxtrialpred_exp{:}),plot_trials,[],plot_depth)));
plot(t,nanmean(AP_index_ans(vertcat(mua_ctxtrialpred_regionzero_exp{:}),plot_trials,[],plot_depth)));
title('Ipsi');
legend({'MUA','MUA ctxpred','MUA ctxtrialpred','MUA ctxtrialpred regionzero'});

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
errorbar(nanmean(ctxpred_r2,1),AP_sem(ctxpred_r2,1),'.','linewidth',2,'CapSize',0,'MarkerSize',50);
errorbar(nanmean(ctxtrialpred_r2,1),AP_sem(ctxpred_r2,1),'.','linewidth',2,'CapSize',0,'MarkerSize',50);
errorbar(nanmean(ctxtrialpred_regionzero_r2,1),AP_sem(ctxpred_r2,1),'.','linewidth',2,'CapSize',0,'MarkerSize',50);

xlabel('Striatum depth');
ylabel('Explained variance');
legend({'Cortex (full)','Cortex (trial)','Cortex (trial,region-zero)'})

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
lambda = 10;
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







