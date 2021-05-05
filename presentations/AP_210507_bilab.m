%% For bi-lab meeting: 2021-05-07

%% Plot average stim response

plot_trials = move_t < 0.5 & trial_stim_allcat == 1 & trial_choice_allcat == -1;

% Make move-aligned fluorescence
fluor_allcat_deconv_move = fluor_allcat_deconv;
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
for i = 1:size(fluor_allcat_deconv,1)
    fluor_allcat_deconv_move(i,:,:,:) = circshift(fluor_allcat_deconv_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
end

avg_px = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    permute(nanmean(fluor_allcat_deconv(plot_trials,:,:),1),[3,2,1]));

avg_px_move = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    permute(nanmean(fluor_allcat_deconv_move(plot_trials,:,:),1),[3,2,1]));

% Plot average around times
px_stim_avg = nanmean(avg_px(:,:,t > 0 & t < 0.1),3);
px_move_avg = nanmean(avg_px_move(:,:,t > 0 & t < 0.1),3);

figure;
subplot(1,2,1);
imagesc(px_stim_avg);
caxis([-max(abs(caxis)),max(abs(caxis))]);
axis image off;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

subplot(1,2,2);
imagesc(px_move_avg);
caxis([-max(abs(caxis)),max(abs(caxis))]);
axis image off;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

colormap(brewermap([],'PrGn'));


%% Cortical task kernels

% Get task>cortex parameters
n_regressors = length(task_regressor_labels);
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

% Concatenate and average task>cortex kernels
fluor_taskpred_k_cat = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(k) k{regressor}, ...
    vertcat(fluor_taskpred_k_all{:}),'uni',false),[2,3,4,1])), ...
    1:n_regressors,'uni',false);

fluor_taskpred_k_avg = cellfun(@(x) nanmean(x,4), fluor_taskpred_k_cat,'uni',false);
fluor_taskpred_k_avg_px = cellfun(@(x) ...
    AP_svdFrameReconstruct(U_master(:,:,1:n_vs),permute(x,[3,2,1])), ...
    fluor_taskpred_k_avg,'uni',false);

% Plot kernel max across time
max_subregressors = max(cellfun(@(x) size(x,1),fluor_taskpred_k_avg));

figure;
for curr_regressor = 1:length(task_regressor_labels)

    c_max = max(fluor_taskpred_k_avg_px{curr_regressor}(:));  
    for curr_subregressor = 1:size(fluor_taskpred_k_avg_px{curr_regressor},4)
        subplot(length(task_regressor_labels),max_subregressors, ...
            (curr_regressor-1)*max_subregressors+curr_subregressor);
        imagesc(max(fluor_taskpred_k_avg_px{curr_regressor}(:,:,:,curr_subregressor),[],3));
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        axis image off;
        colormap(brewermap([],'Greys'));
        caxis([0,c_max]);
        title([task_regressor_labels{curr_regressor}]);       
    end
    
end

px_stim_kernel_max = squeeze(max(fluor_taskpred_k_avg_px{1}(:,:,:,end),[],3));
px_move_kernel_max = squeeze(max(fluor_taskpred_k_avg_px{2}(:,:,:,1),[],3));

figure;

subplot(1,2,1);
imagesc(px_stim_kernel_max);
caxis([-max(abs(caxis)),max(abs(caxis))]);
axis image off;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

subplot(1,2,2);
imagesc(px_move_kernel_max);
caxis([-max(abs(caxis)),max(abs(caxis))]);
axis image off;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

colormap(brewermap([],'PrGn'));

%% Trained vs naive

data_fns = { ...
    'trial_activity_AP_choiceWorldStimPassive_naive', ...
    {'trial_activity_AP_choiceWorldStimPassive_trained', ...
    'trial_activity_AP_lcrGratingPassive_ctxstrephys_str', ...
    'trial_activity_AP_lcrGratingPassive_pre_muscimol'}};

% % (leave out original dataset: 1s instead of 0.5s stim)
% data_fns = { ...
%     'trial_activity_AP_choiceWorldStimPassive_naive', ...
%     {'trial_activity_AP_lcrGratingPassive_ctxstrephys_str', ...
%     'trial_activity_AP_lcrGratingPassive_pre_muscimol'}};

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

fluor_avg = cellfun(@(fluor,stim) nanmean(cell2mat(permute(cellfun(@(fluor,stim) ...
    permute(nanmean(fluor(stim == 1,:,:),1),[3,2,1]),fluor,stim,'uni',false), ...
    [2,3,1])),3),fluor_training,stimIDs,'uni',false);

fluor_avg_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    cat(3,fluor_avg{:}));

fluor_stim_px = squeeze(nanmean(fluor_avg_px(:,:,t > 0 & t < 0.1,:),3));


figure;

subplot(1,2,1);
imagesc(fluor_stim_px(:,:,1));
caxis([-max(abs(caxis)),max(abs(caxis))]);
axis image off;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

subplot(1,2,2);
imagesc(fluor_stim_px(:,:,2));
caxis([-max(abs(caxis)),max(abs(caxis))]);
axis image off;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

colormap(brewermap([],'PrGn'));


















