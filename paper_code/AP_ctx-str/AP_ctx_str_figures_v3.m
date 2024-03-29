% Generate figures for ctx-str paper 
% (trials data is prepared in AP_ctx_str_trial_preprocessing)
%
% Submission 1
% (note: changed slightly after submission when working on new figs)

%% Load in task data

% Load data

% (task)
data_fn = 'trial_activity_choiceworld'; % Primary dataset
% data_fn = 'trial_activity_choiceworld_15strdepth'; % Depth-aligned striatum
exclude_data = true;

% (passive)
% data_fn = 'trial_activity_AP_choiceWorldStimPassive_trained';
% data_fn = 'trial_activity_AP_choiceWorldStimPassive_naive';
% data_fn = 'trial_activity_stimKalatsky_naive';
% data_fn = 'trial_activity_stimKalatsky_trained';
% exclude_data = false;

% (unused at the moment)
% data_fn = 'trial_activity_choiceworld_wfonly'; % Widefield-only days (no craniotomy, so cleaner)
% exclude_data = true;

AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_recording;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));


%% Fig 1a, S1: Psychometric and reaction time

% Plot psychometric
stim_conditions = unique(trial_stim_allcat);
[~,stim_idx] = ismember(trial_stim_allcat,stim_conditions,'rows');

trial_stim_idx_allcat_exp = mat2cell(stim_idx,use_split,1);
trial_choice_allcat_exp = mat2cell(trial_choice_allcat,use_split,1);

frac_orient_right = cell2mat(cellfun(@(stim,choice) ...
    accumarray(stim,choice == -1,[length(stim_conditions),1],@nanmean,NaN), ...
    trial_stim_idx_allcat_exp,trial_choice_allcat_exp,'uni',false)');

figure; hold on; axis square;
plot(stim_conditions,frac_orient_right,'color',[0.5,0.5,0.5]);
plot(stim_conditions,nanmean(frac_orient_right,2),'k','linewidth',3);
line([0,0],ylim,'color','k','linestyle','--');
line(xlim,[0.5,0.5],'color','k','linestyle','--');
xlabel('Stimulus side and contrast');
ylabel('Fraction orient right');

% Plot reaction time by trial percentile within session
n_trial_prctile = 4;
trial_bin = arrayfun(@(x) min(floor(linspace(1,n_trial_prctile+1,x)), ...
    n_trial_prctile)',trials_recording,'uni',false);

move_t_bins = -0.2:1/sample_rate:1;
move_t_bin_centers = move_t_bins(1:end-1) + diff(move_t_bins)./2;
move_t_bin = mat2cell(discretize(move_t,move_t_bins),trials_recording,1);

move_t_hist = cell2mat(permute(cellfun(@(trial_bin,move_t_bin) ...
    accumarray([trial_bin(~isnan(move_t_bin)),move_t_bin(~isnan(move_t_bin))], ...
    1/sum(~isnan(move_t_bin)),[n_trial_prctile,length(move_t_bins)-1],@nansum,0), ...
    trial_bin,move_t_bin,'uni',false),[2,3,1]));

figure;
for curr_trial_prctile = 1:n_trial_prctile
    subplot(n_trial_prctile,1,curr_trial_prctile); hold on;
    plot(move_t_bin_centers, ...
        squeeze(move_t_hist(curr_trial_prctile,:,:)), ...
        'color',[0.5,0.5,0.5]);
    plot(move_t_bin_centers, ...
        nanmean(squeeze(move_t_hist(curr_trial_prctile,:,:)),2), ...
        'color','k','linewidth',3);
    xlabel('Time from stim onset');
    ylabel('Fraction of reaction times');
    line([0,0],ylim,'color','k');
    title(['Trial percentile ' num2str(curr_trial_prctile)]);
end
linkaxes(get(gcf,'Children'),'xy');


%% Fig 1b: Example recording

animal = 'AP025'; 
day = '2017-10-01'; 
experiment = 1; 
str_align = 'none'; 
verbose = true; 
AP_load_experiment;

%%% Plot example data

% Align U's, deconvolve widefield
use_components = 1:200;
aUdf = AP_align_widefield(Udf,animal,day);
fVdf_deconv = AP_deconv_wf(fVdf);

% Set time to plot
plot_t = [134,152];

raster_fig = figure;

% (wheel velocity)
wheel_axes = subplot(6,1,6);
plot_wheel_idx = Timeline.rawDAQTimestamps >= plot_t(1) & ...
    Timeline.rawDAQTimestamps <= plot_t(2);
plot(wheel_axes,Timeline.rawDAQTimestamps(plot_wheel_idx), ...
    wheel_velocity(plot_wheel_idx),'k','linewidth',2);
ylabel('Wheel velocity');
axis off

% (stimuli)
stim_col = colormap_BlueWhiteRed(5);
[~,trial_contrast_idx] = ...
    ismember(trial_conditions(:,1).*trial_conditions(:,2),unique(contrasts'.*sides),'rows');
stim_lines = arrayfun(@(x) line(wheel_axes,repmat(stimOn_times(x),1,2),ylim(wheel_axes),'color', ...
    stim_col(trial_contrast_idx(x),:),'linewidth',2), ...
    find(stimOn_times >= plot_t(1) & stimOn_times <= plot_t(2)));

% (movement starts)
move_col = [0.6,0,0.6;0,0.6,0];
[~,trial_choice_idx] = ismember(trial_conditions(:,3),[-1;1],'rows');
move_lines = arrayfun(@(x) line(wheel_axes,repmat(wheel_move_time(x),1,2),ylim(wheel_axes),'color', ...
    move_col(trial_choice_idx(x),:),'linewidth',2), ...
    find(wheel_move_time >= plot_t(1) & wheel_move_time <= plot_t(2)));

% (go cues)
go_col = [0.8,0.8,0.2];
go_cue_times = signals_events.interactiveOnTimes(1:n_trials);
go_cue_lines = arrayfun(@(x) line(wheel_axes,repmat(go_cue_times(x),1,2),ylim(wheel_axes),'color', ...
    go_col,'linewidth',2), ...
    find(go_cue_times >= plot_t(1) & go_cue_times <= plot_t(2)));

% (outcomes)
outcome_col = [0,0,0.8;0.5,0.5,0.5];
reward_lines = arrayfun(@(x) line(wheel_axes,repmat(reward_t_timeline(x),1,2),ylim(wheel_axes),'color', ...
    outcome_col(1,:),'linewidth',2), ...
    find(reward_t_timeline >= plot_t(1) & reward_t_timeline <= plot_t(2)));
punish_times = signals_events.responseTimes(trial_outcome == -1);
punish_lines = arrayfun(@(x) line(wheel_axes,repmat(punish_times(x),1,2),ylim(wheel_axes),'color', ...
    outcome_col(2,:),'linewidth',2), ...
    find(punish_times >= plot_t(1) & punish_times <= plot_t(2)));

% (striatum raster)
raster_axes = subplot(6,1,3:5,'YDir','reverse'); hold on;
plot_spikes = spike_times_timeline >= plot_t(1) & ...
    spike_times_timeline <= plot_t(2) & ...
    spike_depths >= str_depth(1) & spike_depths <= str_depth(2);
plot(raster_axes,spike_times_timeline(plot_spikes),spike_depths(plot_spikes),'.k');
ylabel('Depth (\mum)');
xlabel('Time (s)');
depth_scale = 1000;
line(repmat(min(xlim),2,1),[min(ylim),min(ylim) + depth_scale],'color','k','linewidth',3);
axis off

% (fluorescence from select ROIs)
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
roi_trace = AP_svd_roi(aUdf(:,:,use_components),fVdf_deconv(use_components,:),[],[],cat(3,wf_roi.mask));

plot_rois = [1,9,10];
fluor_spacing = 0.04;
fluor_axes = subplot(6,1,1:2); hold on;
plot_fluor_idx = frame_t >= plot_t(1) & frame_t <= plot_t(2);
AP_stackplot(roi_trace(plot_rois,plot_fluor_idx)', ...
    frame_t(plot_fluor_idx),fluor_spacing,false,[0,0.7,0],{wf_roi(plot_rois).area});

y_scale = 0.02;
t_scale = 2;
line([min(xlim),min(xlim) + t_scale],repmat(min(ylim),2,1),'color','k','linewidth',3);
line(repmat(min(xlim),2,1),[min(ylim),min(ylim) + y_scale],'color','k','linewidth',3);
axis off

linkaxes([wheel_axes,raster_axes,fluor_axes],'x');

% % Write legend
% [~,unique_contrasts_h] = unique(trial_contrast_idx);
% [~,unique_move_h] = unique(trial_choice_idx(trial_choice_idx > 0));
% legend([stim_lines(unique_contrasts_h),move_lines(unique_move_h), ...
%     go_cue_lines(1),reward_lines(1),punish_lines(1)], ...
%     [cellfun(@(x) ['Stim ' num2str(x)],num2cell(unique(contrasts'.*sides)),'uni',false); ...
%     {'Move L';'Move R';'Go cue';'Reward';'Punish'}]);

% Plot ROIs
figure; hold on
set(gca,'YDir','reverse');
AP_reference_outline('ccf_aligned','k');
for curr_roi = plot_rois
    curr_roi_boundary = cell2mat(bwboundaries(wf_roi(curr_roi).mask()));
    patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),[0,0.8,0]);   
    text(nanmean(curr_roi_boundary(:,2)),nanmean(curr_roi_boundary(:,1)), ...
        wf_roi(curr_roi).area,'FontSize',12,'HorizontalAlignment','center')
end
axis image off;

%% Fig 1c: Average stim-aligned cortex

% Get average stim-aligned fluorescence 
plot_trials = move_t < 0.5 & trial_stim_allcat == 1 & trial_choice_allcat == -1;
plot_trials_exp = mat2cell(plot_trials,use_split,1);

fluor_allcat_deconv_exp = mat2cell(fluor_allcat_deconv,use_split,length(t),n_vs);

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_stim_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set windows to average activity
use_align_labels = {'Stim','Move onset','Outcome'};
use_align = {stim_align,move_align,outcome_align};
plot_t = [0.08,0,0.08];

figure;
for curr_align = 1:length(use_align)
    
    % (re-align activity)
    curr_ctx_act = cellfun(@(act,trials,shift) cell2mat(arrayfun(@(trial) ...
        circshift(act(trial,:,:),shift(trial),2), ...
        find(trials),'uni',false)), ...
        fluor_allcat_deconv_exp,plot_trials_exp, ...
        mat2cell(use_align{curr_align},use_split,1),'uni',false);
    
    curr_ctx_act_mean = ...
        permute(nanmean(cell2mat(cellfun(@(x) nanmean(x,1), ...
        curr_ctx_act,'uni',false)),1),[3,2,1]);
    
    curr_ctx_act_mean_t = interp1(t,curr_ctx_act_mean',plot_t(curr_align))';
    curr_ctx_act_mean_t_px = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
        curr_ctx_act_mean_t);
    
    subplot(length(use_align),1,curr_align);
    imagesc(curr_ctx_act_mean_t_px);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    axis image off;
    colormap(brewermap([],'PRGn'));
    caxis([-0.02,0.02]);
    title([use_align_labels{curr_align} ': ' num2str(plot_t(curr_align)) ' sec']);
    
end


%% Fig 1c,d, S2: Task > cortex kernels

% Get task>cortex parameters
n_regressors = length(task_regressor_labels);
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

% Get average task > cortex kernels (V's and ROIs)
regressor_v = cell(n_regressors,1);
regressor_roi = cell(n_regressors,1);
for curr_regressor = 1:n_regressors
    
    curr_k = cell2mat(cellfun(@(x) x{curr_regressor}, ...
        permute(vertcat(fluor_taskpred_k_all{:}),[2,3,4,1]),'uni',false));
    
    curr_k_roi = nan(size(curr_k,1),size(curr_k,2),n_rois,size(curr_k,4));
    for curr_subregressor = 1:size(curr_k,1)
        for curr_exp = 1:size(curr_k,4)
            curr_k_roi(curr_subregressor,:,:,curr_exp) = ...
                permute(AP_svd_roi(U_master(:,:,1:n_vs), ...
                permute(curr_k(curr_subregressor,:,:,curr_exp),[3,2,1]), ...
                [],[],cat(3,wf_roi.mask)),[3,2,1]);        
        end
    end
    
    curr_k_v = nanmean(curr_k,4);
    
    regressor_v{curr_regressor} = curr_k_v;
    regressor_roi{curr_regressor} = curr_k_roi;
    
    AP_print_progress_fraction(curr_regressor,n_regressors);
end

% Plot example frame for each kernel group
plot_subregressor = [10,1,1,1];
plot_t = [0.08,0,0.05,0.08];

figure;
for curr_regressor = 1:length(task_regressor_labels)
    subplot(length(task_regressor_labels),1,curr_regressor);
    
    curr_k_v_t = interp1(task_regressor_t_shifts{curr_regressor}, ...
        permute(regressor_v{curr_regressor}(plot_subregressor(curr_regressor),:,:),[2,3,1]), ...
        plot_t(curr_regressor))';
    curr_k_px_t = svdFrameReconstruct(U_master(:,:,1:n_vs),curr_k_v_t);
    
    imagesc(curr_k_px_t);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    axis image off;
    colormap(brewermap([],'Greys'));
%     caxis([-0.02,0.02]);
    caxis([0,max(curr_k_px_t(:))]);
    title([task_regressor_labels{curr_regressor} ': ' num2str(plot_t(curr_regressor)) ' sec']);
    
end

% Get ROI traces for each subregressor
plot_rois = [1,9,10];

stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = 0.5;
move_col = [0.6,0,0.6;0.8,0.5,0];
go_col = [0.5,0.5,0.5;0.8,0.8,0.2];
outcome_col = [0,0,0.8;0.8,0,0];
task_regressor_cols = {stim_col,move_col,go_col,outcome_col};
figure;
p = nan(length(plot_rois),n_regressors);
for curr_roi_idx = 1:length(plot_rois)
    curr_roi = plot_rois(curr_roi_idx);    
    for curr_regressor = 1:n_regressors
        p(curr_roi_idx,curr_regressor) = ...
            subplot(length(plot_rois),n_regressors, ...
            sub2ind([n_regressors,length(plot_rois)],curr_regressor,curr_roi_idx)); hold on;
        
        for curr_subregressor = 1:size(regressor_roi{curr_regressor},1)
           AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
            nanmean(regressor_roi{curr_regressor}(curr_subregressor,:,curr_roi,:),4), ...
            AP_sem(regressor_roi{curr_regressor}(curr_subregressor,:,curr_roi,:),4), ...
            task_regressor_cols{curr_regressor}(curr_subregressor,:), ...
            0.5,true);
        end
        ylabel([wf_roi(curr_roi).area ' \DeltaF/F_0']);
        xlabel('Time (s)');
        line([0,0],ylim,'color','k');
    end
end
linkaxes(get(gcf,'Children'),'x')

% Plot ROI average trace and task fit
plot_rois = [1,9,10];
plot_trials = move_t < 0.5 & trial_stim_allcat > 0 & trial_choice_allcat == -1;

figure;
for curr_roi_idx = 1:length(plot_rois)
    subplot(length(plot_rois),1,curr_roi_idx);
    
    curr_roi = plot_rois(curr_roi_idx);
    
    AP_errorfill(t,nanmean(fluor_roi_deconv(plot_trials,:,curr_roi),1), ...
        AP_sem(fluor_roi_deconv(plot_trials,:,curr_roi),1),'k',0.5,true);
    AP_errorfill(t,nanmean(fluor_roi_taskpred(plot_trials,:,curr_roi),1), ...
        AP_sem(fluor_roi_taskpred(plot_trials,:,curr_roi),1),'b',0.5,true);  
    
end

% Plot task>cortex ROI regression examples
plot_rois = [1,9,10];
plot_prctiles = [25,50,75];

figure;
for curr_roi_idx = 1:length(plot_rois)
    
    curr_roi = plot_rois(curr_roi_idx);
    
    % Set current data (pad trials with NaNs for spacing)
    n_pad = 10;
    curr_data = padarray(fluor_roi_deconv(:,:,curr_roi),[0,n_pad],NaN,'post');
    curr_pred_data = padarray(fluor_roi_taskpred(:,:,curr_roi),[0,n_pad],NaN,'post');
    nan_samples = isnan(curr_data) | isnan(curr_pred_data);

    % Smooth
    n_smooth = 1;
    smooth_filt = ones(1,n_smooth)/n_smooth;
    curr_data = conv2(curr_data,smooth_filt,'same');
    curr_pred_data = conv2(curr_pred_data,smooth_filt,'same');
    
    % Set common NaNs for R^2    
    curr_data_nonan = curr_data; 
    curr_data_nonan(nan_samples) = NaN;
    
    curr_pred_data_nonan = curr_pred_data; 
    curr_pred_data(nan_samples) = NaN; 
    
    % Get squared error for each trial
    trial_r2 = 1 - (nansum((curr_data_nonan-curr_pred_data_nonan).^2,2)./ ...
        nansum((curr_data_nonan-nanmean(curr_data_nonan,2)).^2,2));
    
    trial_r2_nonan_idx = find(~isnan(trial_r2));
    [~,trial_r2_rank] = sort(trial_r2(trial_r2_nonan_idx));
    
    plot_prctile_trials = round(prctile(1:length(trial_r2_nonan_idx),plot_prctiles));
    plot_trials = trial_r2_nonan_idx(trial_r2_rank(plot_prctile_trials));
    
    
    curr_t = (1:length(reshape(curr_data(plot_trials,:)',[],1)))/sample_rate;
    
    subplot(length(plot_rois),1,curr_roi_idx); hold on;
    plot(curr_t,reshape(curr_data(plot_trials,:)',[],1),'color',[0,0.6,0],'linewidth',2);
    plot(curr_t,reshape(curr_pred_data(plot_trials,:)',[],1),'b','linewidth',2);
    title([wf_roi(curr_roi,1).area ', Percentiles: ' num2str(plot_prctiles)]);
    axis off
    
end
linkaxes(get(gcf,'Children'),'xy');
y_scale = 10e-3;
t_scale = 0.5;
line([min(xlim),min(xlim)+t_scale],repmat(min(ylim),2,1),'linewidth',3,'color','k');
line(repmat(min(xlim),2,1),[min(ylim),min(ylim)+y_scale],'linewidth',3,'color','k');



%% Fig 1e: Striatal unit PSTH/raster examples

animal = 'AP025'; 
day = '2017-10-01'; 
experiment = 1; 
str_align = 'kernel'; 
load_parts.ephys = true;
verbose = true; 
AP_load_experiment;

outcome_time = signals_events.responseTimes';

% AP_cellraster({stimOn_times,wheel_move_time,outcome_time}, ...
%     {trial_conditions(1:n_trials,1).*trial_conditions(1:n_trials,2), ...
%     trial_choice(1:n_trials),trial_outcome(1:n_trials)});

plot_trials = ...
    (trial_conditions(1:n_trials,1).*trial_conditions(1:n_trials,2)) > 0 & ...
    trial_choice(1:n_trials) == -1 & ...
    trial_outcome(1:n_trials == 1);

AP_cellraster({stimOn_times(plot_trials), ...
    wheel_move_time(plot_trials),outcome_time(plot_trials)});

% Examples used: 
% stim (depth 1) = 545
% move (depth 2) = 328
% move (depth 3) = 258
% reward (depth 4) = 79 (alt: 467 late, 109/96 sharp)

%% Fig 1f: Striatum multiunit end-of-striatum depth aligned

% Plot average stimulus-aligned activity in striatum
plot_trials = move_t < 0.5 & trial_stim_allcat > 0 & trial_choice_allcat == -1;
plot_trials_exp = mat2cell(plot_trials,use_split,1);

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_stim_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set windows to average activity
use_align_labels = {'Stim','Move onset','Outcome'};
use_align = {stim_align,move_align,outcome_align};

figure;
p = gobjects(n_depths,1);
align_col = [1,0,0;0.8,0,0.8;0,0,0.8];
for curr_align = 1:length(use_align)
        
    % (re-align and split activity)
    curr_str_act = mat2cell(...
        cell2mat(arrayfun(@(trial) circshift( ...
        mua_allcat(trial,:,:), ...
        use_align{curr_align}(trial),2),transpose(1:size(mua_allcat,1)),'uni',false)), ...
        use_split,length(t),size(mua_allcat,3));
    
    curr_str_act_plottrial_mean = cell2mat(permute( ...
        cellfun(@(act,trials) permute(nanmean(act(trials,:,:),1),[3,2,1]), ...
        curr_str_act,plot_trials_exp,'uni',false),[2,3,1]));
        
    curr_t_offset = (nanmean(stim_align(cell2mat(plot_trials_exp)) - ...
        use_align{curr_align}(cell2mat(plot_trials_exp))))/sample_rate;
    
    for curr_depth = 1:n_depths
        p(curr_depth,1) = subplot(n_depths,1,curr_depth); hold on;
        
        AP_errorfill(t + curr_t_offset,nanmean(curr_str_act_plottrial_mean(curr_depth,:,:),3)', ...
            AP_sem(curr_str_act_plottrial_mean(curr_depth,:,:),3)',align_col(curr_align,:),0.5);
        xlabel(['Time from ' use_align_labels{curr_align}]);
        ylabel('Spikes (std)');
        line(repmat(curr_t_offset,2,1),ylim,'color','k');
    end
    
end
linkaxes(get(gcf,'Children'),'xy');
xlim([-0.1,1]);

y_scale = 1;
t_scale = 0.2;
line([min(xlim),min(xlim)+t_scale],repmat(min(ylim),2,1),'linewidth',3,'color','k');
line(repmat(min(xlim),2,1),[min(ylim),min(ylim)+y_scale],'linewidth',3,'color','k');




%% Fig 2a: Example cortex > striatum regression by depth

% Set example experiment to use
animal = 'AP025';
day = '2017-10-04';
experiment = 1;

% Parameters for regression
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 60;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

% Load full data
str_align = 'none';
verbose = true;
AP_load_experiment

% Get time points to query
sample_rate = framerate*regression_params.upsample_factor;
time_bins = frame_t(find(frame_t > ...
    regression_params.skip_seconds,1)):1/sample_rate: ...
    frame_t(find(frame_t-frame_t(end) < ...
    -regression_params.skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Deconvolve and resample V
fVdf_deconv = AP_deconv_wf(fVdf);
fVdf_deconv(isnan(fVdf_deconv)) = 0;
fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';

% Get striatum multiunit in 4 depths for example plot
n_depths = 4;
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
[depth_group_n,depth_group] = histc(spike_depths,depth_group_edges);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

binned_spikes = zeros(n_depths,length(time_bins)-1);
for curr_depth = 1:n_depths
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
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

% Regress fluorescence to spikes
kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
    round(regression_params.kernel_t(2)*sample_rate);

[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
    binned_spikes_std,kernel_frames,lambda, ...
    regression_params.zs,regression_params.cvfold, ...
    false,regression_params.use_constant);

Udf_aligned = single(AP_align_widefield(animal,day,Udf));
k_px = zeros(size(Udf_aligned,1),size(Udf_aligned,2),size(k,2),size(k,3),'single');
for curr_spikes = 1:size(k,3)
    k_px(:,:,:,curr_spikes) = ...
        svdFrameReconstruct(Udf_aligned(:,:,regression_params.use_svs),k(:,:,curr_spikes));
end

% NaN-out depths with no spikes
k_px(:,:,:,~any(binned_spikes,2)) = NaN;

% Keep kernel one frame (t == 0)
k_px_frame = squeeze(k_px(:,:,kernel_frames == 0,:));

% Define ROIs and get fluorescence traces
roi_circle_size = 20;
roi_x = [136,176,142,116]; % roi_x = [131,174,110,51];
roi_y = [299,91,79,87]; % roi_y = [297,96,71,144];
[x,y] = meshgrid(1:size(Udf_aligned,1),1:size(Udf_aligned,2));
roi_mask = cell2mat(arrayfun(@(roi) sqrt((x-roi_x(roi)).^2 + (y-roi_y(roi)).^2) <= ...
    roi_circle_size,permute(1:length(roi_x),[1,3,2]),'uni',false));
roi_trace = AP_svd_roi(Udf_aligned,fVdf_deconv_resample,[],[],roi_mask);

% Plot ROIs
figure;
for i = 1:n_depths
   subplot(n_depths,1,i,'YDir','reverse'); hold on;
   curr_roi_boundaries = cell2mat(bwboundaries(roi_mask(:,:,i)));
   plot(curr_roi_boundaries(:,2),curr_roi_boundaries(:,1),'color',[0,0.8,0],'linewidth',2);
   axis image off;
   AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
end

% Plot overlaid striatal activity and and ROI fluorescence
figure; hold on
AP_stackplot(binned_spikes',time_bin_centers,10,true,'k');
AP_stackplot(roi_trace',time_bin_centers,10,true,[0,0.7,0]);
xlim([185,205]);

% Plot kernel frames
figure;
for i = 1:n_depths
   subplot(n_depths,1,i);
   imagesc(k_px_frame(:,:,i)); hold on;
   axis image off;
   colormap(brewermap([],'*RdBu'));
   caxis([-0.012,0.012]); 
   AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
end



%% Fig 2c: Cortex > striatum projections (Allen connectivity database)

% Define the probe vector manually according to the targeted trajectory
probe_vector_ccf = [520,240,510;520,511,239];

%%% Get the average relative depth of each kernel template

% Load kernels by depths, get depth relative to maximum extent
kernel_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
kernel_fn = ['ephys_kernel_depth'];
load([kernel_path filesep kernel_fn])
total_depths = 1:max(cellfun(@(x) size(x,3),[ephys_kernel_depth.k_px]));
k_px_depths = cellfun(@(x) total_depths(end-size(x,3)+1:end),[ephys_kernel_depth.k_px],'uni',false);

% Load the kernel template matches
n_aligned_depths = 4;
kernel_match_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
kernel_match_fn = ['ephys_kernel_align_' num2str(n_aligned_depths) '_depths.mat'];
load([kernel_match_path filesep kernel_match_fn]);

% Concatenate all relative depths and kernel matches
k_depths = cell2mat(k_px_depths);
k_matches = cell2mat([ephys_kernel_align.kernel_match]')';
k_match_depths_relative = grpstats(k_depths,k_matches)./max(total_depths);

%%% Query allen at each point using targeted trajectory

% Load in the annotated Allen volume and names
allen_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
av = readNPY([allen_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean

% Get probe location per micron
probe_size = pdist2(probe_vector_ccf(1,:),probe_vector_ccf(2,:))*10;
probe_depths = ...
    round([linspace(probe_vector_ccf(1,1)',probe_vector_ccf(2,1)',probe_size); ...
    linspace(probe_vector_ccf(1,2)',probe_vector_ccf(2,2)',probe_size); ...
    linspace(probe_vector_ccf(1,3)',probe_vector_ccf(2,3)',probe_size)]');

% Eliminiate trajectory points that are off the atlas
eliminate_depths = ...
    probe_depths(:,1) < 1 | probe_depths(:,1) > size(av,1) | ...
    probe_depths(:,2) < 1 | probe_depths(:,2) > size(av,2) | ...
    probe_depths(:,3) < 1 | probe_depths(:,3) > size(av,3);
probe_depths(eliminate_depths,:) = [];

% Convert probe depths subscripts to indicies
probe_depths_ind = sub2ind(size(av),probe_depths(:,1),probe_depths(:,2),probe_depths(:,3));

% Get structures that the probe is in
probe_structures = av(probe_depths_ind);

% Get target relative depths through striatum
str_id = find(strcmp(st.safe_name,'Caudoputamen'));
probe_structures_str = probe_structures == str_id;
probe_str = [probe_depths(find(probe_structures_str,1,'first'),:); ...
    probe_depths(find(probe_structures_str,1,'last'),:)];
kernel_depth_ccf = interp1([0,1],probe_str,k_match_depths_relative);

%%%% Just use regular depths?
regular_centers_borders = linspace(0,1,n_aligned_depths*2+1);
kernel_depth_ccf = interp1([0,1],probe_str,regular_centers_borders(2:2:end));

% Plot brain to overlay probes
% (note the CCF is rotated to allow for dim 1 = x)
h = figure; ccf_axes = axes; hold on
slice_spacing = 10;
target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) > 1,[2,1,3]);
structure_patch = isosurface(target_volume,0);
structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
target_structure_color = [0.7,0.7,0.7];
brain_outline = patch('Vertices',structure_wire.vertices*slice_spacing, ...
    'Faces',structure_wire.faces, ...
    'FaceColor','none','EdgeColor',target_structure_color);

str_id = find(strcmp(st.safe_name,'Caudoputamen'));
target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) == str_id,[2,1,3]);
structure_patch = isosurface(target_volume,0);
structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
target_structure_color = [0.7,0,0.7];
striatum_outline = patch('Vertices',structure_wire.vertices*slice_spacing, ...
    'Faces',structure_wire.faces, ...
    'FaceColor','none','EdgeColor',target_structure_color);

axis image vis3d off;
view([-30,25]);
cameratoolbar(h,'SetCoordSys','y');
cameratoolbar(h,'SetMode','orbit');

scatter3(kernel_depth_ccf(:,1),kernel_depth_ccf(:,2),kernel_depth_ccf(:,3), ...
    100,copper(n_aligned_depths),'filled');

% Mirror all of the locations across the midline and get projections from
% both (because the cortical injections in the database aren't evenly
% distributed, so this is more comprehensive)
kernel_depth_um = round(kernel_depth_ccf*10);
bregma_um = allenCCFbregma*10;
ccf_midline = bregma_um(3);
hemisphere = sign(kernel_depth_um(1,3) - ccf_midline);
str_depths_mirror = [kernel_depth_um(:,1:2),ccf_midline - hemisphere*abs(kernel_depth_um(:,3)-ccf_midline)];

str_depths_query = [kernel_depth_um;str_depths_mirror];
injection_parameters = get_allen_projection(str_depths_query);
injection_coordinates = {injection_parameters.coordinates};

% Standardize injection coordinates by hemisphere (left = contra, right =
% ipsi)
injection_coordinates_standardized = injection_coordinates;
for curr_coord = 1:length(injection_coordinates)
    
    target_hemisphere = sign(ccf_midline - str_depths_query(curr_coord,3));
    injection_coords_ml_offset = abs(injection_coordinates{curr_coord}(:,3) - ccf_midline);
    injection_coordinates_hemisphere = sign(injection_coordinates{curr_coord}(:,3) - ccf_midline);
    
    injection_coords_ipsi = injection_coordinates_hemisphere == target_hemisphere;
    injection_coords_contra = injection_coordinates_hemisphere == -target_hemisphere;
    
    injection_coordinates_standardized{curr_coord}(injection_coords_ipsi,3) = ...
        ccf_midline + injection_coords_ml_offset(injection_coords_ipsi);
    injection_coordinates_standardized{curr_coord}(injection_coords_contra,3) = ...
        ccf_midline - injection_coords_ml_offset(injection_coords_contra);
    
end

% Combine side-standardized injection data across hemispheres
injection_coordinates_bilateral = arrayfun(@(x) ...
    vertcat(injection_coordinates_standardized{x}, ...
    injection_coordinates_standardized{n_aligned_depths+x}), ...
    1:n_aligned_depths,'uni',false);

% Convert points from CCF to widefield
alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment';
ccf_tform_fn = [alignment_path filesep 'ccf_tform.mat'];
load(ccf_tform_fn);

um2pixel = 20.6;
injection_coordinates_bilateral_wf = cellfun(@(x) ...
    [x(:,[3,1]).*(1/(um2pixel)),ones(size(x,1),1)]*ccf_tform.T, ...
    injection_coordinates_bilateral,'uni',false);

% Get projection density
projection_strength = {injection_parameters.density};

projection_strength_bilateral = arrayfun(@(x) ...
    vertcat(projection_strength{x}', ...
    projection_strength{n_aligned_depths+x}'), ...
    1:n_aligned_depths,'uni',false);

% Load kernel templates for overlay
n_aligned_depths = 4;
kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template_' num2str(n_aligned_depths) '_depths.mat'];
load(kernel_template_fn);

figure; 
colormap(brewermap([],'*RdBu'));
for curr_depth = 1:n_aligned_depths
    subplot(n_aligned_depths,1,curr_depth); hold on; axis image off;
    set(gca,'YDir','reverse');
    
%     imagesc(kernel_template(:,:,curr_depth));
%     caxis([-prctile(abs(kernel_template(:)),99),prctile(abs(kernel_template(:)),99)]);
    
    scatter(injection_coordinates_bilateral_wf{curr_depth}(:,1), ...
        injection_coordinates_bilateral_wf{curr_depth}(:,2), ...
        projection_strength_bilateral{curr_depth}*50 + 10, ...
        'k','filled');
    
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
end

% Color all injection sites and plot overlaid
jet_basic = jet(255);
dark_colors = max(jet_basic,[],2) ~= 1;
jet_alt = max(interp1(1:255,jet_basic,linspace(find(~dark_colors,1,'first'), ...
    find(~dark_colors,1,'last'),255)) - 0.2,0);

cmap = jet_alt;
% depth_color_idx = round(conv(linspace(1,size(cmap,1),n_aligned_depths+1),[0.5,0.5],'valid'));
depth_color_idx = round(linspace(1,size(cmap,1),n_aligned_depths));
plot_colors = cmap(depth_color_idx,:);

figure;
hold on; set(gca,'YDir','reverse'); axis image off;
AP_reference_outline('ccf_aligned','k');
for curr_depth = 1:n_aligned_depths
    % (plot points from both hemispheres)
    scatter(injection_coordinates_bilateral_wf{curr_depth}(:,1), ...
        injection_coordinates_bilateral_wf{curr_depth}(:,2), ...
        projection_strength_bilateral{curr_depth}*100 + 10, ...
        plot_colors(curr_depth,:),'filled');
end


%% Fig 2b,d: Average cortex > striatum domain regression kernels

protocols = {'vanillaChoiceworld','stimSparseNoiseUncorrAsync'};

for protocol = protocols 
    protocol = cell2mat(protocol);
    
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
    k_fn = [data_path filesep 'ctx_str_kernels_' protocol];
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
    AP_imscroll(reshape(permute(k_px,[1,4,2,3]),size(k_px,1)*size(k_px,4),size(k_px,2),length(t)),t);
    colormap(brewermap([],'*RdBu'));
    caxis([-max(caxis),max(caxis)]);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5],[],[size(k_px,1),size(k_px,2),size(k_px,4),1]);
    axis image off
    
    drawnow;
    
end


%% Fig 3, S6: Striatal domain activity and task regression

% Plot stim-aligned/sorted measured and predicted striatum activity
% (correct contra trials)
for curr_trial_set = 1:2
    switch curr_trial_set
        case 1
            plot_trials = move_t < Inf & trial_stim_allcat > 0 & trial_choice_allcat == -1;
            figure('Name','Correct contra trials'); 
        case 2
            plot_trials = move_t < Inf & trial_stim_allcat < 0 & trial_choice_allcat == 1;
            figure('Name','Correct ipsi trials'); 
    end
    
    p = gobjects(n_depths,4);
    colormap(brewermap([],'Greys'));
    for curr_depth = 1:n_depths
        
        % Get trials to plot, sort by reaction time
        curr_trials = plot_trials & ~all(isnan(mua_allcat(:,:,curr_depth)),2);
        curr_trials_idx = find(curr_trials);
        [~,rxn_sort_idx] = sort(move_t(curr_trials_idx));
        
        sorted_plot_trials = curr_trials_idx(rxn_sort_idx);
        
        curr_plot = mua_allcat(sorted_plot_trials,:,curr_depth);
        curr_taskpred_plot = mua_taskpred_allcat(sorted_plot_trials,:,curr_depth);
        curr_ctxpred_plot = mua_ctxpred_allcat(sorted_plot_trials,:,curr_depth);
        
        % Smooth and plot with stim/move/reward times
        % (as conv(nans-zeroed)./conv(non-nan) to ignore in nans in conv)
        smooth_filt = [50,1]; % (trials x frames)
        
        curr_plot_smooth = conv2(curr_plot,ones(smooth_filt),'same')./ ...
            conv2(~isnan(curr_plot),ones(smooth_filt),'same');
        
        curr_taskpred_plot_smooth = curr_taskpred_plot;
        curr_taskpred_plot_smooth(isnan(curr_taskpred_plot_smooth)) = 0;
        curr_taskpred_plot_smooth = conv2(curr_taskpred_plot_smooth,ones(smooth_filt),'same')./ ...
            conv2(~isnan(curr_taskpred_plot),ones(smooth_filt),'same');
        
        curr_ctxpred_plot_smooth = curr_ctxpred_plot;
        curr_ctxpred_plot_smooth(isnan(curr_ctxpred_plot_smooth)) = 0;
        curr_ctxpred_plot_smooth = conv2(curr_ctxpred_plot_smooth,ones(smooth_filt),'same')./ ...
            conv2(~isnan(curr_ctxpred_plot),ones(smooth_filt),'same');
        
        p(curr_depth,1) = subplot(n_depths,4,1+(curr_depth-1)*4,'YDir','reverse'); hold on
        imagesc(t,[],curr_plot_smooth);
        plot(zeros(size(curr_trials_idx)),1:length(curr_trials_idx),'MarkerSize',1,'color','r');
        plot(move_t(sorted_plot_trials),1:length(curr_trials_idx),'MarkerSize',1,'color',[0.8,0,0.8]);
%         plot(outcome_t(sorted_plot_trials),1:length(curr_trials_idx),'.','MarkerSize',1,'color','b');
        axis tight;
        xlim([-0.2,1]);
        xlabel('Time from stim');
        ylabel('Trials (rxn sorted)');
        title('Measured');
        
        p(curr_depth,2) = subplot(n_depths,4,2+(curr_depth-1)*4,'YDir','reverse'); hold on
        imagesc(t,[],curr_taskpred_plot_smooth);
        plot(zeros(size(curr_trials_idx)),1:length(curr_trials_idx),'MarkerSize',1,'color','r');
        plot(move_t(sorted_plot_trials),1:length(curr_trials_idx),'MarkerSize',1,'color',[0.8,0,0.8]);
%         plot(outcome_t(sorted_plot_trials),1:length(curr_trials_idx),'.','MarkerSize',1,'color','b');
        axis tight;
        xlim([-0.2,1]);
        xlabel('Time from stim');
        ylabel('Trials (rxn sorted)');
        title('Task-predicted');
        
        p(curr_depth,3) = subplot(n_depths,4,3+(curr_depth-1)*4,'YDir','reverse'); hold on
        imagesc(t,[],curr_ctxpred_plot_smooth);
        plot(zeros(size(curr_trials_idx)),1:length(curr_trials_idx),'MarkerSize',1,'color','r');
        plot(move_t(sorted_plot_trials),1:length(curr_trials_idx),'MarkerSize',1,'color',[0.8,0,0.8]);
%         plot(outcome_t(sorted_plot_trials),1:length(curr_trials_idx),'.','MarkerSize',1,'color','b');
        axis tight;
        xlim([-0.2,1]);
        xlabel('Time from stim');
        ylabel('Trials (rxn sorted)');
        title('Cortex-predicted');
        
        % Split and average trials by animal
        curr_trials_exp = mat2cell(curr_trials,use_split,1);
        
        curr_mua_exp = mat2cell(mua_allcat(:,:,curr_depth),use_split,length(t));
        curr_mua_exp_mean = cell2mat(cellfun(@(data,trials) nanmean(data(trials,:),1), ...
            curr_mua_exp,curr_trials_exp,'uni',false));
        
        curr_taskpred_mua_exp = mat2cell(mua_taskpred_allcat(:,:,curr_depth),use_split,length(t));
        curr_taskpred_mua_exp_mean = cell2mat(cellfun(@(data,trials) nanmean(data(trials,:),1), ...
            curr_taskpred_mua_exp,curr_trials_exp,'uni',false));
        
        curr_ctxpred_mua_exp = mat2cell(mua_ctxpred_allcat(:,:,curr_depth),use_split,length(t));
        curr_ctxpred_mua_exp_mean = cell2mat(cellfun(@(data,trials) nanmean(data(trials,:),1), ...
            curr_ctxpred_mua_exp,curr_trials_exp,'uni',false));
        
        % Plot PSTH (measured, task-predicted, cortex-predicted);
        p(curr_depth,4) = subplot(n_depths,4,4+(curr_depth-1)*4); hold on
        p1 = AP_errorfill(t,nanmean(curr_mua_exp_mean,1)', ...
            AP_sem(curr_mua_exp_mean,1)','k',0.5);
        p2 = AP_errorfill(t,nanmean(curr_taskpred_mua_exp_mean,1)', ...
            AP_sem(curr_taskpred_mua_exp_mean,1)',[0,0,0.7],0.5);
        p3 = AP_errorfill(t,nanmean(curr_ctxpred_mua_exp_mean,1)', ...
            AP_sem(curr_ctxpred_mua_exp_mean,1)',[0,0.7,0],0.5);
        xlim([-0.2,1])
        line([0,0],ylim,'color','r');
        line(repmat(median(move_t(sorted_plot_trials)),1,2),ylim,'color',[0.8,0,0.8],'linestyle','--');
        line(repmat(median(outcome_t(sorted_plot_trials)),1,2),ylim,'color','b','linestyle','--');
        xlabel('Time from stim');
        ylabel('Spikes (std)');
        legend([p1,p2,p3],{'Measured','Task-predicted','Cortex-predicted'});
        
    end
    % Link the x-axes, set the c/y-axes same within a row
    linkaxes(p(:),'x');
    
    for curr_row = 1:size(p,1)
        curr_ylim = ylim(p(curr_row,4));
        caxis(p(curr_row,1),[0,curr_ylim(2)]);
        caxis(p(curr_row,2),[0,curr_ylim(2)]);
        caxis(p(curr_row,3),[0,curr_ylim(2)]);
    end
    
    trial_scale = 500;
    t_scale = 0.5;
    y_scale = 1;
    line(p(1,1),min(xlim(p(1,1))) + [0,t_scale],repmat(min(ylim(p(1,1))),2,1),'color','k','linewidth',3);
    line(p(1,4),min(xlim(p(1,4))) + [0,t_scale],repmat(min(ylim(p(1,4))),2,1),'color','k','linewidth',3);
    line(p(1,1),repmat(min(xlim(p(1,1))),2,1),min(ylim(p(1,1))) + [0,trial_scale],'color','k','linewidth',3);
    line(p(1,4),repmat(min(xlim(p(1,4))),2,1),min(ylim(p(1,4))) + [0,y_scale],'color','k','linewidth',3);
    
end

% Get task>striatum parameters
n_regressors = length(task_regressor_labels);

% Normalize task > striatum kernels across experiments with mua_norm
mua_taskpred_k_allcat_norm = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    vertcat(mua_taskpred_k_all{:}),vertcat(mua_norm{:}),'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

% Plot task>striatum kernels
stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = [];
move_col = [0.6,0,0.6;0,0.6,0];
go_col = [0.8,0.8,0.2;0.5,0.5,0.5];
outcome_col = [0,0,0.8;0.8,0,0];
task_regressor_cols = {stim_col,move_col,go_col,outcome_col};
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

figure;
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors  
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);  
        
        curr_kernels = mua_taskpred_k_allcat_norm{curr_regressor}(:,:,curr_depth,:);
        n_subregressors = size(mua_taskpred_k_allcat_norm{curr_regressor},1);
        col = task_regressor_cols{curr_regressor};
        for curr_subregressor = 1:n_subregressors
            AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
                nanmean(curr_kernels(curr_subregressor,:,:,:),4)', ...
                AP_sem(curr_kernels(curr_subregressor,:,:,:),4)', ...
                col(curr_subregressor,:),0.5);
        end
        
        xlabel('Time (s)');
        ylabel('Weight');
        title(task_regressor_labels{curr_regressor});
        line([0,0],ylim,'color','k');
        
    end
end
linkaxes(p);
y_scale = 1;
t_scale = 0.5;
line(min(xlim) + [0,t_scale],repmat(min(ylim),2,1),'color','k','linewidth',3);
line(repmat(min(xlim),2,1),min(ylim) + [0,y_scale],'color','k','linewidth',3);


%% Fig 4a,b: Striatum vs Cortex-predicted (and task response fix)

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(mua_allcat,1),1);
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set windows to average activity
timeavg_labels = {'Stim','Move early','Move late','Outcome'};
timeavg_t = {[0.05,0.15],[-0.05,0.05],[0.05,0.15],[0.05,0.15]};
timeavg_align = {stim_align,move_align,move_align,outcome_align};
timeavg_trial_conditions = ...
    {[trial_stim_allcat > 0,trial_stim_allcat == 0,trial_stim_allcat < 0], ...
    [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
    [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
    [trial_outcome_allcat == 1, trial_outcome_allcat == -1]};
timeavg_condition_colors = ...
    {[1,0,0;0,0,0;0,0,1], ...
    [0.6,0,0.6;0,0.6,0], ...
    [0.6,0,0.6;0,0.6,0], ...
    [0,0,0.7;0,0,0]};
timeavg_event_offset = {'Stim','Move onset','Move onset','Outcome'};

% Set activity percentiles and bins
act_prctile = [5,95];
n_act_bins = 5;

% Set areas and conditions
plot_areas = [1:n_depths];

% Loop across area pairs, plot binned predicted v measured activity
curr_act_allcat = mua_allcat;
curr_act_pred_allcat = mua_ctxpred_allcat;

% (old: to use all events)
% task_fix = mua_taskpred_allcat - mua_ctxpred_taskpred_allcat;
% curr_act_pred_fix_allcat = curr_act_pred_allcat + task_fix;

% (to use individual events)
task_fix = (mua_taskpred_allcat - mua_taskpred_reduced_allcat) - ...
    (mua_ctxpred_taskpred_allcat - mua_ctxpred_taskpred_reduced_allcat);
curr_act_pred_fix_allcat = curr_act_pred_allcat + task_fix;

measured_v_pred_fig = figure('color','w');
for curr_area_idx = 1:length(plot_areas)   
    
    plot_area = plot_areas(curr_area_idx);  
    
    for curr_timeavg = 1:length(timeavg_labels)
        
        trial_conditions = timeavg_trial_conditions{curr_timeavg};
        curr_event_offset_idx = find(strcmp(timeavg_event_offset{curr_timeavg},task_regressor_labels));
        curr_act_pred_fix_allcat_event = ...
            curr_act_pred_fix_allcat(:,:,:,curr_event_offset_idx);
        
        % (re-align and split activity)
        act_title = timeavg_labels{curr_timeavg};
        
        curr_act = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        curr_act_pred = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_pred_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        curr_act_predfix = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_pred_fix_allcat_event(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_fix_allcat_event,1)),'uni',false)), ...
            use_split,length(t));
        
        trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
        
        % (get average activity within window)
        curr_event_t = t >= timeavg_t{curr_timeavg}(1) & t <= timeavg_t{curr_timeavg}(2);
        curr_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act,'uni',false);
        curr_act_pred_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_pred,'uni',false);
        curr_act_predfix_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_predfix,'uni',false);
        
        % (bin predicted data across percentile range)
        pred_bin_edges = prctile(cell2mat(curr_act_pred_avg),linspace(act_prctile(1),act_prctile(2),n_act_bins+1));
        pred_bin_centers = pred_bin_edges(1:end-1) + diff(pred_bin_edges)./2;        
        pred_trial_bins = cellfun(@(x) discretize(x,pred_bin_edges),curr_act_pred_avg,'uni',false);         
        
        % (get activity binned by predicted)        
        pred_use_trials = cellfun(@(act,act_pred,trial_bins) ...
            squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
            squeeze(~any(isnan(act_pred(:,curr_event_t)),2)) & ...
            ~isnan(trial_bins), ...
            curr_act,curr_act_pred,pred_trial_bins,'uni',false);
        
        act_predbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_predbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        % (get "fixed" predicted activity binned by predicted)
        act_predfix_predbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_predfix_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        % Get condition difference significance from shuffle
        n_shuff = 1000;
        trial_conditions_shuff = cellfun(@(x) AP_shake(x,1), ...
            repmat(trial_conditions_exp,1,n_shuff),'uni',false);
        
        act_predbinmean_condshuff = cell2mat(arrayfun(@(curr_shuff) ...
            cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,pred_trial_bins,trial_conditions_shuff(:,curr_shuff), ...
            pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false)), ...
            permute(1:n_shuff,[1,3,4,2]),'uni',false));
        
        act_predfix_predbinmean_condshuff = cell2mat(arrayfun(@(curr_shuff) ...
            cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_predfix_avg,pred_trial_bins,trial_conditions_shuff(:,curr_shuff), ...
            pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false)), ...
            permute(1:n_shuff,[1,3,4,2]),'uni',false));
        
        % (measured data: null = no difference between conditions)
        cond_combos = nchoosek(1:size(trial_conditions,2),2);
        cond_sig_diff = false(n_act_bins,size(cond_combos,1));
        predfix_cond_sig_diff = false(n_act_bins,size(cond_combos,1));
        for curr_cond_combo = 1:size(cond_combos,1)
            curr_combo_diff = nanmean(act_predbinmean(:,:,cond_combos(curr_cond_combo,1)) - ...
                 act_predbinmean(:,:,cond_combos(curr_cond_combo,2)),2);
            shuff_prctile = squeeze(prctile(nanmean(act_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,1),:) - ...
                act_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,2),:),2),[2.5,97.5],4));
            cond_sig_diff(:,curr_cond_combo) = ...
                curr_combo_diff < shuff_prctile(:,1) | ...
                curr_combo_diff > shuff_prctile(:,2);
            
            predfix_curr_combo_diff = nanmean(act_predfix_predbinmean(:,:,cond_combos(curr_cond_combo,1)) - ...
                 act_predfix_predbinmean(:,:,cond_combos(curr_cond_combo,2)),2);
            predfix_shuff_prctile = squeeze(prctile(nanmean(act_predfix_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,1),:) - ...
                act_predfix_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,2),:),2),[2.5,97.5],4));
            predfix_cond_sig_diff(:,curr_cond_combo) = ...
                predfix_curr_combo_diff < predfix_shuff_prctile(:,1) | ...
                predfix_curr_combo_diff > predfix_shuff_prctile(:,2);
        end 
      
        % Plot binned measured, predicted, and error (by predicted bins)
        measured_pred_fig = figure('color','w','Name', ...
            ['Str ' num2str(plot_area)' ', ' timeavg_labels{curr_timeavg}]);
        n_col_bins = n_act_bins + 2;
        
        [binned_act_pred_t,binned_act_pred_grp] = grpstats(cell2mat(curr_act_pred), ...
            [cell2mat(pred_trial_bins),trial_conditions],{'nanmean','gname'});
        binned_act_pred_grp = cellfun(@str2num,binned_act_pred_grp);        
        
        [binned_act_t,binned_act_grp] = grpstats(cell2mat(curr_act), ...
            [cell2mat(pred_trial_bins),trial_conditions],{'nanmean','gname'});
        binned_act_grp = cellfun(@str2num,binned_act_grp);
        
        % (plot predicted data timecourse)
        for curr_cond = 1:size(trial_conditions,2)           
            subplot(2,size(trial_conditions,2), ...
                sub2ind(fliplr([2,size(trial_conditions,2)]),curr_cond,1)); hold on;
                     
            curr_mean = nanmean(cell2mat(cellfun(@(act,use_trials,cond) ...
                act(use_trials & cond(:,curr_cond),:), ...
                curr_act_pred,pred_use_trials,trial_conditions_exp,'uni',false)),1);
            curr_std = nanstd(cell2mat(cellfun(@(act,use_trials,cond) ...
                act(use_trials & cond(:,curr_cond),:), ...
                curr_act_pred,pred_use_trials,trial_conditions_exp,'uni',false)),[],1);
            
            AP_errorfill(t,curr_mean',curr_std','k',0.5,true);
            
            set(gca,'ColorOrder',brewermap(n_col_bins,'*Greens'));
            plot(t,binned_act_pred_t(binned_act_pred_grp(:,curr_cond + 1) == 1,:)','linewidth',2);
            xlabel('Time'); ylabel('Predicted data'); 
            title(['Condition ' num2str(curr_cond)]);
            line(repmat(timeavg_t{curr_timeavg}(1),2,1),ylim,'color','k');
            line(repmat(timeavg_t{curr_timeavg}(2),2,1),ylim,'color','k');            
        end
 
        % (plot measured data timecourse)        
        for curr_cond = 1:size(trial_conditions,2)           
            subplot(2,size(trial_conditions,2), ...
                sub2ind(fliplr([2,size(trial_conditions,2)]),curr_cond,2)); hold on;
            set(gca,'ColorOrder',[brewermap(n_col_bins,'*Greys')]);
            plot(t,binned_act_t(binned_act_grp(:,curr_cond + 1) == 1,:)','linewidth',2);
            xlabel('Time'); ylabel('Measured data'); 
            title(['Condition ' num2str(curr_cond)]);
            line(repmat(timeavg_t{curr_timeavg}(1),2,1),ylim,'color','k');
            line(repmat(timeavg_t{curr_timeavg}(2),2,1),ylim,'color','k');            
        end        
        
        linkaxes(get(measured_pred_fig,'Children'),'xy');       
        
        % Plot measured v binned predicted
        figure(measured_v_pred_fig)
        
        % (measured vs binned predicted)
        subplot(length(plot_areas),length(timeavg_labels), ...
            sub2ind(fliplr([length(plot_areas),length(timeavg_labels)]),curr_timeavg,curr_area_idx));
        hold on;
        set(gca,'ColorOrder',timeavg_condition_colors{curr_timeavg});

        fill_cols = min(timeavg_condition_colors{curr_timeavg} + 0.5,1);
        for curr_cond = 1:size(trial_conditions,2)
            AP_errorfill( ...
                squeeze(nanmean(act_pred_predbinmean(:,:,curr_cond),2)), ...
                squeeze(nanmean(act_predfix_predbinmean(:,:,curr_cond),2)), ...
                squeeze(AP_sem(act_predfix_predbinmean(:,:,curr_cond),2)), ...
                fill_cols(curr_cond,:),1,false);
        end
        
        errorbar( ...
            squeeze(nanmean(act_pred_predbinmean,2)), ...
            squeeze(nanmean(act_predbinmean,2)), ...
            squeeze(AP_sem(act_predbinmean,2)), ...
            'linestyle','none','linewidth',2,'CapSize',10);
        
        ctx_col = brewermap(n_col_bins,'*Greens');        
        scatter( ...
            reshape(squeeze(nanmean(act_pred_predbinmean,2)),[],1), ...
            reshape(squeeze(nanmean(act_predbinmean,2)),[],1), ...
            60,repmat(ctx_col(1:n_act_bins,:),size(trial_conditions,2),1),'filled');
        
        xlabel(['Predicted (' num2str(plot_area) ')']);
        ylabel(['Measured (' num2str(plot_area) ')'])
        axis square tight;
        xlim(xlim + [-0.1,0.1]);
        title(act_title);
        
        % (plot significant measured condition differences)
        % (* and o = significant in both measured and "fixed" predicted)
        curr_ylim = max(ylim);
        for curr_cond_combo = 1:size(cond_combos,1)
            % (plot * for measured condition differences)
            sig_y = curr_ylim + 0.1*curr_cond_combo;
            sig_x = pred_bin_centers;
            if any(cond_sig_diff(:,curr_cond_combo))
                plot(sig_x(cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                    '*','MarkerSize',10,'color', ...
                    timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,1),:));
                plot(sig_x(cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                    '*','MarkerSize',5,'color', ...
                    timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,2),:));
            end
            % (plot o for [predicted condition differences)
            if any(predfix_cond_sig_diff(:,curr_cond_combo))
                plot(sig_x(predfix_cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                    'o','MarkerSize',15,'color', ...
                    timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,1),:));
                plot(sig_x(predfix_cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                    'o','MarkerSize',10,'color', ...
                    timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,2),:));
            end
        end
                   
        drawnow;
    end
end


%% Fig 4c: Striatum vs Cortex-predicted (passive choiceworld)
% (do separately for trained and naive)

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(mua_allcat,1),1);

% Get trials with movement during stim to exclude
wheel_thresh = 0.025;
quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > wheel_thresh,2);

% Set windows to average activity
timeavg_labels = {'Stim'};
timeavg_t = {[0.05,0.15]};
timeavg_align = {stim_align};
timeavg_trial_conditions = ...
    {[trial_stim_allcat > 0 & quiescent_trials, ...
    trial_stim_allcat < 0 & quiescent_trials]};
timeavg_condition_colors = ...
    {[1,0,0;0,0,1]};

% Set activity percentiles and bins
act_prctile = [5,95];
n_act_bins = 5;

% Set areas and conditions
plot_areas = [1];

% Loop across area pairs, plot binned predicted v measured activity
curr_act_allcat = mua_allcat;

% (ctx-predicted)
curr_act_pred_allcat = mua_ctxpred_allcat;

% (fix by average stim response within experiment)
mua_allcat_exp = mat2cell(mua_allcat,use_split,length(t),n_depths);
mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,use_split,length(t),n_depths);
trial_stim_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
trial_conditions_exp = mat2cell(timeavg_trial_conditions{1},use_split, ...
    size(timeavg_trial_conditions{1},2));

str_stim_exp = nan(length(unique(trial_stim_allcat)),length(t),n_depths,length(use_split));
ctx_stim_exp = nan(length(unique(trial_stim_allcat)),length(t),n_depths,length(use_split));
curr_act_pred_fix_allcat_exp = mua_ctxpred_allcat_exp;
for curr_exp = 1:length(use_split)
    
    curr_stim = unique(trial_stim_allcat_exp{curr_exp});
    for curr_stim_idx = 1:length(curr_stim)
        curr_trials = any(trial_conditions_exp{curr_exp},2) & ...
            trial_stim_allcat_exp{curr_exp} == curr_stim(curr_stim_idx);
        curr_act_stim_avg = nanmean(mua_allcat_exp{curr_exp}(curr_trials,:,:),1);
        curr_act_pred_stim_avg = nanmean(mua_ctxpred_allcat_exp{curr_exp}(curr_trials,:,:),1);
        
        str_stim_exp(curr_stim_idx,:,:,curr_exp) = curr_act_stim_avg;
        ctx_stim_exp(curr_stim_idx,:,:,curr_exp) = curr_act_pred_stim_avg;
        
        curr_stim_fix = curr_act_stim_avg - curr_act_pred_stim_avg;        
        curr_act_pred_fix_allcat_exp{curr_exp}(curr_trials,:,:) = ...
            curr_act_pred_fix_allcat_exp{curr_exp}(curr_trials,:,:) + curr_stim_fix;       
    end
end
curr_act_pred_fix_allcat = cell2mat(curr_act_pred_fix_allcat_exp);

% Plot stim response
stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = 0.5;
stim_col_contrastside = unique([0,0.06,0.125,0.25,0.5,1].*[-1;1]);
[~,used_stim_idx] = ismember(unique(trial_stim_allcat), ...
    stim_col_contrastside,'rows');

str_stim_avg = nanmean(str_stim_exp,4);
figure('color','w');
for curr_area_idx = 1:length(plot_areas)
    curr_area = plot_areas(curr_area_idx);
    subplot(length(plot_areas),2,length(plot_areas)*(curr_area_idx-1)+1); 
    hold on;
    for curr_stim_idx = 1:size(str_stim_avg,1)
        AP_errorfill(t, ...
            nanmean(str_stim_exp(curr_stim_idx,:,curr_area,:),4)', ...
            AP_sem(str_stim_exp(curr_stim_idx,:,curr_area,:),4)', ...
            stim_col(used_stim_idx(curr_stim_idx),:),0.5,true);
    end
    line([0,0],ylim,'color','k');
    ylabel(['Measured (' num2str(curr_area) ')']);
    xlabel('Time from stim');
    
    subplot(length(plot_areas),2,length(plot_areas)*(curr_area_idx-1)+2); 
    hold on;
    for curr_stim_idx = 1:size(str_stim_avg,1)
        AP_errorfill(t, ...
            nanmean(ctx_stim_exp(curr_stim_idx,:,curr_area,:),4)', ...
            AP_sem(ctx_stim_exp(curr_stim_idx,:,curr_area,:),4)', ...
            stim_col(used_stim_idx(curr_stim_idx),:),0.5,true);
    end
    line([0,0],ylim,'color','k');
    ylabel(['Predicted (' num2str(curr_area) ')']);
    xlabel('Time from stim');
end
linkaxes(get(gcf,'Children'),'xy');

measured_v_pred_fig = figure('color','w');
for curr_area_idx = 1:length(plot_areas)   
    
    plot_area = plot_areas(curr_area_idx);  
    
    for curr_timeavg = 1:length(timeavg_labels)
        
        trial_conditions = timeavg_trial_conditions{curr_timeavg};
        
        % (re-align and split activity)
        act_title = timeavg_labels{curr_timeavg};
        
        curr_act = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        curr_act_pred = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_pred_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        curr_act_predfix = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_pred_fix_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_fix_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
        
        % (get average activity within window)
        curr_event_t = t >= timeavg_t{curr_timeavg}(1) & t <= timeavg_t{curr_timeavg}(2);
        curr_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act,'uni',false);
        curr_act_pred_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_pred,'uni',false);
        curr_act_predfix_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_predfix,'uni',false);
        
        % (bin predicted data across percentile range)
        pred_bin_edges = prctile(cell2mat(curr_act_pred_avg),linspace(act_prctile(1),act_prctile(2),n_act_bins+1));
        pred_bin_centers = pred_bin_edges(1:end-1) + diff(pred_bin_edges)./2;        
        pred_trial_bins = cellfun(@(x) discretize(x,pred_bin_edges),curr_act_pred_avg,'uni',false);         
        
        % (get activity binned by predicted)        
        pred_use_trials = cellfun(@(act,act_pred,trial_bins) ...
            squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
            squeeze(~any(isnan(act_pred(:,curr_event_t)),2)) & ...
            ~isnan(trial_bins), ...
            curr_act,curr_act_pred,pred_trial_bins,'uni',false);
        
        act_predbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_predbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        % (get "fixed" predicted activity binned by predicted)
        act_predfix_predbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_predfix_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        % Get condition difference significance from shuffle
        n_shuff = 1000;
        trial_conditions_shuff = cellfun(@(x) AP_shake(x,1), ...
            repmat(trial_conditions_exp,1,n_shuff),'uni',false);
        
        act_predbinmean_condshuff = cell2mat(arrayfun(@(curr_shuff) ...
            cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,pred_trial_bins,trial_conditions_shuff(:,curr_shuff), ...
            pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false)), ...
            permute(1:n_shuff,[1,3,4,2]),'uni',false));
        
        act_predfix_predbinmean_condshuff = cell2mat(arrayfun(@(curr_shuff) ...
            cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_predfix_avg,pred_trial_bins,trial_conditions_shuff(:,curr_shuff), ...
            pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false)), ...
            permute(1:n_shuff,[1,3,4,2]),'uni',false));
        
        % (measured data: null = no difference between conditions)
        cond_combos = nchoosek(1:size(trial_conditions,2),2);
        cond_sig_diff = false(n_act_bins,size(cond_combos,1));
        predfix_cond_sig_diff = false(n_act_bins,size(cond_combos,1));
        for curr_cond_combo = 1:size(cond_combos,1)
            curr_combo_diff = nanmean(act_predbinmean(:,:,cond_combos(curr_cond_combo,1)) - ...
                 act_predbinmean(:,:,cond_combos(curr_cond_combo,2)),2);
            shuff_prctile = squeeze(prctile(nanmean(act_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,1),:) - ...
                act_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,2),:),2),[2.5,97.5],4));
            cond_sig_diff(:,curr_cond_combo) = ...
                curr_combo_diff < shuff_prctile(:,1) | ...
                curr_combo_diff > shuff_prctile(:,2);
            
            predfix_curr_combo_diff = nanmean(act_predfix_predbinmean(:,:,cond_combos(curr_cond_combo,1)) - ...
                 act_predfix_predbinmean(:,:,cond_combos(curr_cond_combo,2)),2);
            predfix_shuff_prctile = squeeze(prctile(nanmean(act_predfix_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,1),:) - ...
                act_predfix_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,2),:),2),[2.5,97.5],4));
            predfix_cond_sig_diff(:,curr_cond_combo) = ...
                predfix_curr_combo_diff < predfix_shuff_prctile(:,1) | ...
                predfix_curr_combo_diff > predfix_shuff_prctile(:,2);
        end 
        
        % Plot measured v predicted in bins
        figure(measured_v_pred_fig)
        
        % (measured vs binned predicted)
        subplot(length(plot_areas),length(timeavg_labels), ...
            sub2ind(fliplr([length(plot_areas),length(timeavg_labels)]),curr_timeavg,curr_area_idx));
        hold on;
        set(gca,'ColorOrder',timeavg_condition_colors{curr_timeavg});

        fill_cols = min(timeavg_condition_colors{curr_timeavg} + 0.5,1);
%         for curr_cond = 1:size(trial_conditions,2)
%             AP_errorfill( ...
%                 squeeze(nanmean(act_pred_predbinmean(:,:,curr_cond),2)), ...
%                 squeeze(nanmean(act_predfix_predbinmean(:,:,curr_cond),2)), ...
%                 squeeze(AP_sem(act_predfix_predbinmean(:,:,curr_cond),2)), ...
%                 fill_cols(curr_cond,:),1,false);
%         end
        
        errorbar( ...
            squeeze(nanmean(act_pred_predbinmean,2)), ...
            squeeze(nanmean(act_predbinmean,2)), ...
            squeeze(AP_sem(act_predbinmean,2)), ...
            'linestyle','-','linewidth',2,'CapSize',10);
        
        n_col_bins = n_act_bins + 2;
        ctx_col = brewermap(n_col_bins,'*Greens');        
        scatter( ...
            reshape(squeeze(nanmean(act_pred_predbinmean,2)),[],1), ...
            reshape(squeeze(nanmean(act_predbinmean,2)),[],1), ...
            60,repmat(ctx_col(1:n_act_bins,:),size(trial_conditions,2),1),'filled');
        
        xlabel(['Predicted (' num2str(plot_area) ')']);
        ylabel(['Measured (' num2str(plot_area) ')'])
        axis square;
        title(act_title);
        
        % (plot significant measured condition differences)
        % (* and o = significant in both measured and "fixed" predicted)
        curr_ylim = max(ylim);
        for curr_cond_combo = 1:size(cond_combos,1)
            % (plot * for measured condition differences)
            sig_y = curr_ylim + 0.1*curr_cond_combo;
            sig_x = pred_bin_centers;
            if any(cond_sig_diff(:,curr_cond_combo))
                plot(sig_x(cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                    '*','MarkerSize',10,'color', ...
                    timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,1),:));
                plot(sig_x(cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                    '*','MarkerSize',5,'color', ...
                    timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,2),:));
            end
            % (plot o for [predicted condition differences)
            if any(predfix_cond_sig_diff(:,curr_cond_combo))
                plot(sig_x(predfix_cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                    'o','MarkerSize',15,'color', ...
                    timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,1),:));
                plot(sig_x(predfix_cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                    'o','MarkerSize',10,'color', ...
                    timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,2),:));
            end
        end
        
        drawnow;
        
    end
end

linkaxes(get(measured_v_pred_fig,'Children'),'xy');



%% Fig S3a: Widefield hemodynamic correction and deconvolution

% Load and plot deconvolution kernel
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\gcamp_kernel\gcamp6s_kernel.mat');
gcamp6s_kernel_cat = vertcat(gcamp6s_kernel.regression{:});
gcamp6s_kernel_cat_norm = gcamp6s_kernel_cat./max(abs(gcamp6s_kernel_cat),[],2);
gcamp6s_kernel_mean = nanmean(gcamp6s_kernel_cat_norm,1);

figure; hold on;
plot(-gcamp6s_kernel.regression_t,gcamp6s_kernel_cat_norm','color',[0.5,0.5,0.5]);
plot(-gcamp6s_kernel.regression_t,gcamp6s_kernel_mean,'color','k','linewidth',2);
line([0,0],ylim,'color','k','linestyle','--');
line(xlim,[0,0],'color','k','linestyle','--');
xlabel('Time from spike');
ylabel('Max-normalized weight');

% Load sample cortical recording
animal = 'AP026';
day = '2017-12-09';
experiment = 2;
verbose = true;
kilosort_version = 1;
AP_load_experiment;

% Get raw, hemo-corrected and deconvolved fluorescence in craniotomy

[roi_trace_raw,roi_mask] = AP_svd_roi(Un,fVn,avg_im);

roi_trace_hemo = AP_svd_roi(Udf,fVdf,avg_im,[],roi_mask);

fVdf_deconv = AP_deconv_wf(fVdf);
roi_trace_deconv = AP_svd_roi(Udf,fVdf_deconv,avg_im,[],roi_mask);

% Get cortical multiunit binned by imaged frame
ctx_depth = [0,1600];

% Get event-aligned activity
raster_window = [-0.5,5];
raster_sample_rate = 1/framerate;
t = raster_window(1):raster_sample_rate:raster_window(2);
t_peri_event = bsxfun(@plus,stimOn_times,t);
t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];

ctx_spikes = spike_times_timeline(spike_depths >= ctx_depth(1) & spike_depths <= ctx_depth(2));
stim_aligned_mua = cell2mat(arrayfun(@(x) ...
    histcounts(ctx_spikes,t_peri_event_bins(x,:)), ...
    [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;

stim_aligned_fluor_raw = interp1(frame_t,roi_trace_raw',t_peri_event);
stim_aligned_fluor_hemo = interp1(frame_t,roi_trace_hemo',t_peri_event);
stim_aligned_fluor_deconv = interp1(frame_t,roi_trace_deconv',t_peri_event);

% Subtract baseline, average across one stimulus
baseline_t = t < 0;
stim_aligned_mua_baselined = stim_aligned_mua - nanmean(stim_aligned_mua(:,baseline_t,:),2);
stim_aligned_fluor_raw_baselined = stim_aligned_fluor_raw - nanmean(stim_aligned_fluor_raw(:,baseline_t,:),2);
stim_aligned_fluor_hemo_baselined = stim_aligned_fluor_hemo - nanmean(stim_aligned_fluor_hemo(:,baseline_t,:),2);
stim_aligned_fluor_deconv_baselined = stim_aligned_fluor_deconv - nanmean(stim_aligned_fluor_deconv(:,baseline_t,:),2);

plot_stim = 2;
stim_mua = nanmean(stim_aligned_mua_baselined(stimIDs == plot_stim,:),1);
stim_fluor_raw = nanmean(stim_aligned_fluor_raw_baselined(stimIDs == plot_stim,:),1);
stim_fluor_hemo = nanmean(stim_aligned_fluor_hemo_baselined(stimIDs == plot_stim,:),1);
stim_fluor_deconv = nanmean(stim_aligned_fluor_deconv_baselined(stimIDs == plot_stim,:),1);

figure; hold on
plot(t,stim_mua./max(stim_mua),'k','linewidth',2);
plot(t,stim_fluor_raw./max(stim_fluor_raw),'g','linewidth',2);
plot(t,stim_fluor_hemo./max(stim_fluor_hemo),'r','linewidth',2);
plot(t,stim_fluor_deconv./max(stim_fluor_deconv),'b','linewidth',2);
line([0,0],ylim,'color','k','linestyle','--');
xlabel('Time from stimulus');
ylabel('Max-normalized activity');
legend({'Multiunit','Fluorescence (Raw)','Fluorescence (Hemo-corrected)','Fluorescence(Deconvolved)'});


%% Fig S3b: Widefield alignment

% Show average images across days for one animal
animal = 'AP025';

protocol = 'vanillaChoiceworld';
experiments = AP_find_experiments(animal,protocol);
experiments = experiments([experiments.imaging] & [experiments.ephys]);

figure;
for curr_day = 1:length(experiments)
    
    day = experiments(curr_day).day;
    
    [img_path,img_exists] = AP_cortexlab_filename(animal,day,[],'imaging');
    avg_im = readNPY([img_path filesep 'meanImage_purple.npy']);
        
    subplot(1,length(experiments),curr_day);
    imagesc(avg_im);
    axis image off;
    colormap(gray);
    caxis([0,40000]);
    title([animal ' Day ' num2str(curr_day)]);
    
end

% Show retinotopy for all animals
% (grab retinotopy from figs)
retinotopy_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\retinotopy';
retinotopy_dir = dir(retinotopy_path);

animal_retinotopy_idx = cellfun(@(x) ~isempty(x), regexp({retinotopy_dir.name},'AP\d*_retinotopy'));
animals_tokens = cellfun(@(x) regexp({x},'(AP\d*)_retinotopy','tokens'),{retinotopy_dir.name});
animals = cellfun(@(x) cell2mat(x{:}),animals_tokens(cellfun(@(x) ~isempty(x),animals_tokens)),'uni',false);

retinotopy_unaligned = cell(size(animals));
for curr_animal = 1:length(animals)
    h = open([retinotopy_path filesep animals{curr_animal} '_retinotopy.fig']);
    retinotopy_unaligned{curr_animal} = get(get(subplot(1,2,1),'Children'),'CData');
    close(h);
end

figure;
for curr_animal = 1:6
   subplot(1,6,curr_animal);
   imagesc(retinotopy_unaligned{curr_animal});
   axis image off
   colormap(brewermap([],'*RdBu'));
   caxis([-1,1]);
   title(animals{curr_animal});
end

% Show colored CCF and master retinotopy with aligned CCF

% (these steps are normally done in AP_vfs_ccf_align)

% Load CCF av and st
ccf_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
av = readNPY([ccf_path filesep 'annotation_volume_10um_by_index.npy']); 
st = loadStructureTree([ccf_path filesep 'structure_tree_safe_2017.csv']);

% Get first brain pixel from top-down, get annotation at that point
[~,top_down_depth] = max(av>1, [], 2);
top_down_depth = squeeze(top_down_depth);

[xx,yy] = meshgrid(1:size(top_down_depth,2), 1:size(top_down_depth,1));
top_down_annotation = reshape(av(sub2ind(size(av),yy(:),top_down_depth(:),xx(:))), size(av,1), size(av,3));

% Get all labelled areas
used_areas = unique(top_down_annotation(:));

% Restrict to only cortical areas
structure_id_path = cellfun(@(x) textscan(x(2:end),'%d', 'delimiter',{'/'}),st.structure_id_path);

ctx_path = [997,8,567,688,695,315];
ctx_idx = find(cellfun(@(id) length(id) > length(ctx_path) & ...
    all(id(min(length(id),length(ctx_path))) == ctx_path(min(length(id),length(ctx_path)))),structure_id_path));

plot_areas = intersect(used_areas,ctx_idx);

bregma = allenCCFbregma;

% Get outlines of all areas
top_down_cortical_area_boundaries = cell(size(plot_areas));
for curr_area_idx = 1:length(plot_areas)
    top_down_cortical_area_boundaries{curr_area_idx} = bwboundaries(top_down_annotation == plot_areas(curr_area_idx));
end

% Color CCF by VFS
a_idx = find(cellfun(@(name) strcmp(name,'Anterior area layer 1'),st.safe_name(used_areas)));
al_idx = find(cellfun(@(name) strcmp(name,'Anterolateral visual area layer 1'),st.safe_name(used_areas)));
am_idx = find(cellfun(@(name) strcmp(name,'Anteromedial visual area layer 1'),st.safe_name(used_areas)));
lm_idx = find(cellfun(@(name) strcmp(name,'Lateral visual area layer 1'),st.safe_name(used_areas)));
v1_idx = find(cellfun(@(name) strcmp(name,'Primary visual area layer 1'),st.safe_name(used_areas)));
p_idx = find(cellfun(@(name) strcmp(name,'Posterolateral visual area layer 1'),st.safe_name(used_areas)));
pm_idx = find(cellfun(@(name) strcmp(name,'posteromedial visual area layer 1'),st.safe_name(used_areas)));
li_idx = find(cellfun(@(name) strcmp(name,'Laterointermediate area layer 1'),st.safe_name(used_areas)));
rl_idx = find(cellfun(@(name) strcmp(name,'Rostrolateral area layer 1'),st.safe_name(used_areas)));

ccf_vfs = zeros(size(top_down_annotation));
ccf_vfs(ismember(top_down_annotation,used_areas([v1_idx,am_idx,al_idx,li_idx]))) = 1;
ccf_vfs(ismember(top_down_annotation,used_areas([a_idx,p_idx,pm_idx,rl_idx,lm_idx]))) = -1;

% Load master retinotopy
combined_retinotopy_filename = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\retinotopy\combined_retinotopy.fig';
h = open(combined_retinotopy_filename);
imaged_vfs = get(get(gca,'Children'),'CData');
close(h);

% Load widefield correlation borders
combined_borders_filename = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_borders\combined_wf_borders.fig';
h = open(combined_borders_filename);
imaged_corr_borders = get(get(gca,'Children'),'CData');
close(h);


% Plot CCF VFS, imaged VFS, and imaged correlation borders
figure;

subplot(1,3,1,'YDir','reverse'); hold on;
imagesc(ccf_vfs);
cellfun(@(area) plot(area(:,2),area(:,1),'color',[0.5,0.5,0.5]), ...
    vertcat(top_down_cortical_area_boundaries{:}),'uni',false);
axis image off;
colormap(brewermap([],'*RdBu'));
caxis([-1,1]);    
title('CCF');

subplot(1,3,2);
imagesc(imaged_vfs);
axis image off;
colormap(brewermap([],'*RdBu'));
caxis([-1,1]);
title('Combined');
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

subplot(1,3,3);
imagesc(imaged_corr_borders);
axis image off;
colormap(brewermap([],'*RdBu'));
caxis([-max(caxis),max(caxis)]);
title('Correlation borders');
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);



%% Fig S4: Striatum border estimation (histology and electrophysiology)

animal = 'AP032';

% Get probe areas estimated from histology
load(['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\histology\' animal '\processed\probe_areas']);

% Get and plot striatum boundaries
ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
ephys_align_fn = ['ephys_depth_align'];
load([ephys_align_path filesep ephys_align_fn])

% Parameters from batch
n_corr_groups = 40;
depth_group_edges = linspace(0,3820,n_corr_groups+1);
depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

curr_animal = find(strcmp({ephys_depth_align.animal},animal));

% Plot first day and histology-estimated boundaries
str_id = 574;
str_start_histology = probe_areas.y(find(probe_areas.av == str_id,1,'first'));
str_end_histology = probe_areas.y(find(probe_areas.av == str_id,1,'last'));

curr_day = 1;
mua_corr = ephys_depth_align(curr_animal).mua_corr{curr_day};
template_depths = ephys_depth_align(curr_animal).template_depths{curr_day};
str_depth = ephys_depth_align(curr_animal).str_depth(curr_day,:);

figure;
imagesc(depth_group_centers,depth_group_centers,mua_corr);
caxis([0,1]); colormap(hot); axis square;
line(xlim,[str_start_histology,str_start_histology],'color','w')
line([str_start_histology,str_start_histology],ylim,'color','w')
line(xlim,[str_end_histology,str_end_histology],'color','w')
line([str_end_histology,str_end_histology],ylim,'color','w')
xlabel('Depth (\mum)');
ylabel('Depth (\mum)');
title(['Day ' num2str(curr_day) ' (histology-approximated)']);
    
% Plot all days and ephys-estimated boundaries
n_days = length(ephys_depth_align(curr_animal).mua_corr);
figure;
for curr_day = 1:n_days
    mua_corr = ephys_depth_align(curr_animal).mua_corr{curr_day};
    template_depths = ephys_depth_align(curr_animal).template_depths{curr_day};
    str_depth = ephys_depth_align(curr_animal).str_depth(curr_day,:);
    
    subplot(1,n_days,curr_day);
    imagesc(depth_group_centers,depth_group_centers,mua_corr);
    caxis([0,1]); colormap(hot); axis square;
    line(xlim,[str_depth(1),str_depth(1)],'color','w');
    line([str_depth(1),str_depth(1)],ylim,'color','w');
    line(xlim,[str_depth(2),str_depth(2)],'color','w');
    line([str_depth(2),str_depth(2)],ylim,'color','w');
    xlabel('Depth (\mum)');
    ylabel('Depth (\mum)');
    title(['Day ' num2str(curr_day)]);
end
set(gcf,'Name',ephys_depth_align(curr_animal).animal);


%% Fig S6: Corticostriatal kernel alignment

% Load kernels by depths
kernel_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
kernel_fn = ['ephys_kernel_depth'];
load([kernel_path filesep kernel_fn])

plot_animals = [4,5,6];
plot_day = 1;

figure;
for curr_animal = 1:length(plot_animals)
    subplot(1,length(plot_animals),curr_animal);
    
    curr_k = ephys_kernel_depth(plot_animals(curr_animal)).k_px{plot_day};
    curr_k_norm = curr_k./nanstd(reshape(curr_k,[],1,size(curr_k,3)),[],1);
    
    
    max(max(abs(curr_k),[],1),[],2);
    
    imagesc(reshape(permute(curr_k_norm,[1,3,2]),[],size(curr_k,2)));
    axis image off;
    caxis([-5,5]);
    colormap(brewermap([],'*RdBu'));
    title(ephys_kernel_depth(plot_animals(curr_animal)).animal);
    
end


%% Fig S8: Striatum task v cortex regression

% Get task>striatum parameters
n_regressors = length(task_regressor_labels);

% Normalize task > striatum kernels across experiments with mua_norm
mua_taskpred_k_allcat_norm = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    vertcat(mua_taskpred_k_all{:}),vertcat(mua_norm{:}),'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

mua_ctxpred_taskpred_k_allcat_norm = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    vertcat(mua_ctxpred_taskpred_k_all{:}),vertcat(mua_norm{:}),'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

% Get sum for each kernel
mua_taskpred_k_allcat_norm_sum = cellfun(@(x) squeeze(sum(x,2)), ...
    mua_taskpred_k_allcat_norm,'uni',false);

mua_ctxpred_taskpred_k_allcat_norm_sum = cellfun(@(x) squeeze(sum(x,2)), ...
    mua_ctxpred_taskpred_k_allcat_norm,'uni',false);

% Plot sum task>striatum and task>ctx-str kernel by condition
figure;
for curr_regressor = 1:n_regressors
    if curr_regressor == 1
        x = unique([0.06,0.125,0.25,0.5,1].*[-1;1]);
    else
        x = 1:size(mua_taskpred_k_allcat_norm_sum{curr_regressor},1);
    end
      
    for curr_depth = 1:n_depths
        subplot(n_depths,n_regressors, ...
            sub2ind([n_regressors,n_depths],curr_regressor,curr_depth));
        hold on
        
        errorbar(x,nanmean(mua_taskpred_k_allcat_norm_sum{curr_regressor}(:,curr_depth,:),3), ...
            AP_sem(mua_taskpred_k_allcat_norm_sum{curr_regressor}(:,curr_depth,:),3), ...
            'color','b','linewidth',2);
        
        errorbar(x,nanmean(mua_ctxpred_taskpred_k_allcat_norm_sum{curr_regressor}(:,curr_depth,:),3), ...
            AP_sem(mua_ctxpred_taskpred_k_allcat_norm_sum{curr_regressor}(:,curr_depth,:),3), ...
            'color',[0,0.7,0],'linewidth',2);
        
        title(task_regressor_labels{curr_regressor});
        ylabel('Sum weight');
        xlabel('Condition');
        
    end
end
ax_handles = reshape(flipud(get(gcf,'Children')),n_depths,curr_regressor);
for curr_depth = 1:n_depths
    linkaxes(ax_handles(curr_depth,:),'y');
end
legend(ax_handles(1,1),{'Task>Str','Task>Ctx-Str'});

% Plot task>striatum regression examples
figure;
plot_prctiles = linspace(0,100,5);
for curr_depth = 1:n_depths   
    
    % Set current data (pad trials with NaNs for spacing)
    n_pad = 10;
    curr_data = padarray(mua_allcat(:,:,curr_depth),[0,n_pad],NaN,'post');
    curr_taskpred_data = padarray(mua_taskpred_allcat(:,:,curr_depth),[0,n_pad],NaN,'post');
    curr_ctxpred_data = padarray(mua_ctxpred_allcat(:,:,curr_depth),[0,n_pad],NaN,'post');
    nan_samples = isnan(curr_data) | isnan(curr_taskpred_data) | isnan(curr_ctxpred_data);

    % Smooth
    n_smooth = 3;
    smooth_filt = ones(1,n_smooth)/n_smooth;
    curr_data = conv2(curr_data,smooth_filt,'same');
    curr_taskpred_data = conv2(curr_taskpred_data,smooth_filt,'same');
    curr_ctxpred_data = conv2(curr_ctxpred_data,smooth_filt,'same');
    
    % Set common NaNs for R^2    
    curr_data_nonan = curr_data; 
    curr_data_nonan(nan_samples) = NaN;
    
    curr_taskpred_data_nonan = curr_taskpred_data; 
    curr_taskpred_data(nan_samples) = NaN; 
    
    curr_ctxpred_data_nonan = curr_ctxpred_data; 
    curr_ctxpred_data(nan_samples) = NaN; 
    
    % Get squared error for each trial (for task prediction)
    trial_r2 = 1 - (nansum((curr_data_nonan-curr_taskpred_data_nonan).^2,2)./ ...
        nansum((curr_data_nonan-nanmean(curr_data_nonan,2)).^2,2));
    
    trial_r2_nonan_idx = find(~isnan(trial_r2));
    [~,trial_r2_rank] = sort(trial_r2(trial_r2_nonan_idx));
    
    plot_prctile_trials = round(prctile(1:length(trial_r2_rank),plot_prctiles));
    plot_trials = trial_r2_nonan_idx(trial_r2_rank(plot_prctile_trials));
    
    curr_t = (1:length(reshape(curr_data(plot_trials,:)',[],1)))/sample_rate;
    
    subplot(n_depths,1,curr_depth); hold on;
    plot(curr_t,reshape(curr_data(plot_trials,:)',[],1),'k','linewidth',2);
    plot(curr_t,reshape(curr_taskpred_data(plot_trials,:)',[],1),'b','linewidth',2);
    plot(curr_t,reshape(curr_ctxpred_data(plot_trials,:)',[],1),'color',[0,0.7,0],'linewidth',2);
    
    xlabel('Time (s)');
    ylabel('Spikes (std)');
    title(['R^2 percentiles plotted: ' num2str(plot_prctiles)]);
    legend({'Measured','Task-predicted','Cortex-predicted'});

end
linkaxes(get(gcf,'Children'),'xy');
y_scale = 2;
t_scale = 1;
line([min(xlim),min(xlim)+t_scale],repmat(min(ylim),2,1),'linewidth',3,'color','k');
line(repmat(min(xlim),2,1),[min(ylim),min(ylim)+y_scale],'linewidth',3,'color','k');


% Get R^2 for task regression 
taskpred_r2 = nan(max(split_idx),n_depths);
ctxpred_r2 = nan(max(split_idx),n_depths);
for curr_exp = 1:max(split_idx)
    curr_data = reshape(permute(mua_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_taskpred_data = reshape(permute(mua_taskpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_ctxpred_data = reshape(permute(mua_ctxpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_taskpred_data) | isnan(curr_ctxpred_data);
    curr_data(nan_samples) = NaN;
    curr_taskpred_data(nan_samples) = NaN;
    curr_ctxpred_data(nan_samples) = NaN;

    taskpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_taskpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
end
figure; hold on;
errorbar(nanmean(taskpred_r2,1),AP_sem(taskpred_r2,1),'b','linewidth',2);
errorbar(nanmean(ctxpred_r2,1),AP_sem(ctxpred_r2,1),'color',[0,0.7,0],'linewidth',2);
xlabel('Striatum depth');
ylabel('Task explained variance');
legend({'Task','Cortex'});


%% Fig S8a (addition): comparison of cortex kernel vs ROI

% Load kernel templates
n_aligned_depths = 3;
kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template_' num2str(n_aligned_depths) '_depths.mat'];
load(kernel_template_fn);

% Plot the kernels and ROIs
figure; colormap(gray);
for i = 1:n_aligned_depths
    p1 = subplot(n_aligned_depths,2,(i-1)*2+1);
    imagesc(kernel_template(:,:,i));
    caxis([-max(abs(caxis)),max(abs(caxis))])
    colormap(p1,brewermap([],'*RdBu'));
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    axis image off;
    
    p2 = subplot(n_aligned_depths,2,(i-1)*2+2);
    imagesc(kernel_roi.bw(:,:,i));
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    colormap(p2,brewermap([],'Greys'));
    axis image off;
end

% Get fluorescence within kernel ROIs
fluor_kernelroi_deconv = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],kernel_roi.bw), ...
    size(kernel_roi.bw,3),[],size(fluor_allcat_deconv,1)),[3,2,1]);

% Regress kernel ROI activity to striatum domain activity (per recording)
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;
lambda = 0;
kernel_frames = floor(regression_params.kernel_t(1)*sample_rate): ...
    ceil(regression_params.kernel_t(2)*sample_rate);

fluor_kernelroi_deconv_exp = mat2cell(fluor_kernelroi_deconv,trials_recording,length(t),n_depths);
mua_allcat_exp = mat2cell(mua_allcat,trials_recording,length(t),n_depths);

mua_ctxroipred_exp = cellfun(@(x) nan(size(x)),mua_allcat_exp,'uni',false);
for curr_exp = 1:length(trials_recording)
    for curr_depth = 1:n_depths
        
        curr_mua = reshape(mua_allcat_exp{curr_exp}(:,:,curr_depth)',[],1)';
        curr_fluor_kernelroi = reshape(fluor_kernelroi_deconv_exp{curr_exp}(:,:,curr_depth)',[],1)';
        
        % Skip if no data
        if all(isnan(curr_mua))
            continue
        end
        
        % Set discontinuities in trial data
        trial_discontinuities = false(size(mua_allcat_exp{curr_exp}(:,:,curr_depth)));
        trial_discontinuities(:,1) = true;
        trial_discontinuities = reshape(trial_discontinuities',[],1)';
        
        % Do regression
        [k,curr_mua_kernelroipred,explained_var] = ...
            AP_regresskernel(curr_fluor_kernelroi, ...
            curr_mua,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant,trial_discontinuities);
              
        mua_ctxroipred_exp{curr_exp}(:,:,curr_depth) = ...
            reshape(curr_mua_kernelroipred,length(t),[])';

    end
    AP_print_progress_fraction(curr_exp,length(trials_recording));
end
mua_ctxroipred_allcat = cell2mat(mua_ctxroipred_exp);

% Plot task>striatum regression examples
figure;
plot_prctiles = [25,50,75];
for curr_depth = 1:n_depths   
    
    % Set current data (pad trials with NaNs for spacing)
    n_pad = 10;
    curr_data = padarray(mua_allcat(:,:,curr_depth),[0,n_pad],NaN,'post');
    curr_taskpred_data = padarray(mua_taskpred_allcat(:,:,curr_depth),[0,n_pad],NaN,'post');
    curr_ctxpred_data = padarray(mua_ctxpred_allcat(:,:,curr_depth),[0,n_pad],NaN,'post');
    curr_ctxroipred_data = padarray(mua_ctxroipred_allcat(:,:,curr_depth),[0,n_pad],NaN,'post');
    nan_samples = isnan(curr_data) | isnan(curr_taskpred_data) | isnan(curr_ctxpred_data) | isnan(curr_ctxroipred_data);

    % Smooth
    n_smooth = 3;
    smooth_filt = ones(1,n_smooth)/n_smooth;
    curr_data = conv2(curr_data,smooth_filt,'same');
    curr_taskpred_data = conv2(curr_taskpred_data,smooth_filt,'same');
    curr_ctxpred_data = conv2(curr_ctxpred_data,smooth_filt,'same');
    curr_ctxroipred_data = conv2(curr_ctxroipred_data,smooth_filt,'same');
    
    % Set common NaNs for R^2    
    curr_data_nonan = curr_data; 
    curr_data_nonan(nan_samples) = NaN;
    
    curr_taskpred_data_nonan = curr_taskpred_data; 
    curr_taskpred_data(nan_samples) = NaN; 
    
    curr_ctxpred_data_nonan = curr_ctxpred_data; 
    curr_ctxpred_data(nan_samples) = NaN; 
    
    curr_ctxroipred_data_nonan = curr_ctxroipred_data; 
    curr_ctxroipred_data(nan_samples) = NaN; 
    
    % Get squared error for each trial (for task prediction)
    trial_r2 = 1 - (nansum((curr_data_nonan-curr_taskpred_data_nonan).^2,2)./ ...
        nansum((curr_data_nonan-nanmean(curr_data_nonan,2)).^2,2));
    
    trial_r2_nonan_idx = find(~isnan(trial_r2));
    [~,trial_r2_rank] = sort(trial_r2(trial_r2_nonan_idx));
    
    plot_prctile_trials = round(prctile(1:length(trial_r2_rank),plot_prctiles));
    plot_trials = trial_r2_nonan_idx(trial_r2_rank(plot_prctile_trials));
    
    curr_t = (1:length(reshape(curr_data(plot_trials,:)',[],1)))/sample_rate;
    
    subplot(n_depths,1,curr_depth); hold on;
    plot(curr_t,reshape(curr_data(plot_trials,:)',[],1),'k','linewidth',2);
    plot(curr_t,reshape(curr_taskpred_data(plot_trials,:)',[],1),'b','linewidth',2);
    plot(curr_t,reshape(curr_ctxpred_data(plot_trials,:)',[],1),'color',[0,0.7,0],'linewidth',2);
    plot(curr_t,reshape(curr_ctxroipred_data(plot_trials,:)',[],1),'color',[1,0.5,0],'linewidth',2);
    
    xlabel('Time (s)');
    ylabel('Spikes (std)');
    title(['R^2 percentiles plotted: ' num2str(plot_prctiles)]);
    legend({'Measured','Task-predicted','Cortex-predicted (full)','Cortex-predicted (ROI)'});

end
linkaxes(get(gcf,'Children'),'xy');
y_scale = 2;
t_scale = 1;
line([min(xlim),min(xlim)+t_scale],repmat(min(ylim),2,1),'linewidth',3,'color','k');
line(repmat(min(xlim),2,1),[min(ylim),min(ylim)+y_scale],'linewidth',3,'color','k');

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
ylabel('Task explained variance');
legend({'Task','Cortex (Full)','Cortex (ROI)'});

% Get significance between cortex kernel and ROI
ctx_kernel_roi_p = nan(n_depths,1);
for curr_depth = 1:n_depths
   ctx_kernel_roi_p(curr_depth) = signrank(ctxroipred_r2(:,curr_depth), ...
       ctxpred_r2(:,curr_depth));
   disp(['Str ' num2str(curr_depth) ' kernel vs ROI: p = ' ...
       num2str(ctx_kernel_roi_p(curr_depth))]);
end


%% Movie 1: Average widefield 

% Get average stim-aligned fluorescence on trial subsets
plot_trials = { ...
    move_t < 0.5 & trial_stim_allcat > 0 & trial_choice_allcat == -1, ...
    move_t < 0.5 & trial_stim_allcat < 0 & trial_choice_allcat == 1, ...
    move_t < 0.5 & trial_stim_allcat == 0 & trial_outcome_allcat == 1, ...
    move_t < 0.5 & trial_stim_allcat > 0 & trial_choice_allcat == 1, ...
    move_t < 0.5 & trial_stim_allcat < 0 & trial_choice_allcat == -1, ...
    move_t < 0.5 & trial_stim_allcat == 0 & trial_outcome_allcat == -1};

plot_trials_exp = cellfun(@(x) mat2cell(x,use_split,1),plot_trials,'uni',false);
fluor_allcat_deconv_exp = mat2cell(fluor_allcat_deconv,use_split,length(t),n_vs);
movie_annotation = {'Correct right stim','Correct left stim','Rewarded 0%', ...
    'Incorrect right stim','Uncorrect left stim','Unrewarded 0%'};

fluor_allcat_deconv_mean_px = cellfun(@(trials) nanmean(cell2mat( ...
    permute(cellfun(@(trials,fluor) ...
    svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    permute(nanmean(fluor(trials,:,:),1),[3,2,1])), ...
    trials,fluor_allcat_deconv_exp,'uni',false),[2,3,4,1])),4), ...
    plot_trials_exp,'uni',false);
    
% Make movie of average fluorescence for all trial groups
movie_rate = sample_rate/5;
color_map = brewermap([],'Greens');
color_axis = [0,0.02];
figure_position = [24,408,1868,399];
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\figs\movies';
savefile = [save_path filesep 'movie_1_avg_fluor'];
t_annotation = cellfun(@(x) sprintf('Time from stimulus: %0.2f sec',x),num2cell(t),'uni',false);
AP_movie2avi(cat(4,fluor_allcat_deconv_mean_px{:}), ...
    movie_rate,color_map,color_axis,figure_position,savefile,t_annotation,movie_annotation);


%% Movie 2-5: Task->cortex kernels

% Get task>cortex parameters
n_regressors = length(task_regressor_labels);
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

% Get average task > cortex kernels (V's and ROIs)
regressor_v = cell(n_regressors,1);
for curr_regressor = 1:n_regressors
    
    curr_k = cell2mat(cellfun(@(x) x{curr_regressor}, ...
        permute(vertcat(fluor_taskpred_k_all{:}),[2,3,4,1]),'uni',false));
    
    curr_k_roi = nan(size(curr_k,1),size(curr_k,2),n_rois,size(curr_k,4));
    for curr_subregressor = 1:size(curr_k,1)
        for curr_exp = 1:size(curr_k,4)
            curr_k_roi(curr_subregressor,:,:,curr_exp) = ...
                permute(AP_svd_roi(U_master(:,:,1:n_vs), ...
                permute(curr_k(curr_subregressor,:,:,curr_exp),[3,2,1]), ...
                [],[],cat(3,wf_roi.mask)),[3,2,1]);        
        end
    end
    
    curr_k_v = nanmean(curr_k,4);
    
    regressor_v{curr_regressor} = curr_k_v;
    
    AP_print_progress_fraction(curr_regressor,n_regressors);
end

% Get regressor pixels
regressor_px = cellfun(@(v) cell2mat(arrayfun(@(subregressor) ...
    svdFrameReconstruct(U_master(:,:,1:n_vs),permute(v(subregressor,:,:),[3,2,1])), ...
    permute(1:size(v,1),[1,3,4,2]),'uni',false)),regressor_v,'uni',false);

% Make movies of regressor groups
movie_rate = sample_rate/5;
% color_map = brewermap([],'PRGn');
% color_axis = [-0.015,0.015];
color_map = brewermap([],'Greys');
figure_position = [24,408,1868,399];
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\figs\movies';

for curr_regressor = 1:length(regressor_px)
    movie_num = 2 + (curr_regressor-1);
    savefile = [save_path filesep 'movie_' num2str(movie_num) '_ctx_' task_regressor_labels{curr_regressor} '_kernel'];
    t_annotation = cellfun(@(x) sprintf('Time from event: %0.2f sec',x), ...
        num2cell(task_regressor_t_shifts{curr_regressor}),'uni',false);
    switch curr_regressor
        case 1
            movie_annotation = {'Left 100% contrast','Left 50% contrast','Left 25% contrast','Left 12.5% contrast','Left 6% contrast',...
                'Right 100% contrast','Right 50% contrast','Right 25% contrast','Right 12.5% contrast','Right 6% contrast'};
            curr_im = regressor_px{curr_regressor}(:,:,:,[1:5,10:-1:6]);
        case 2
            movie_annotation = {'Orient right','Orient left'};
            curr_im = regressor_px{curr_regressor};
        case 3
            movie_annotation = {'Go cue (already moving)','Go cue (not moving)'};
            curr_im = regressor_px{curr_regressor};
        case 4
            movie_annotation = {'Rewarded','Not rewarded'};
            curr_im = regressor_px{curr_regressor};
    end
    
    color_axis = [0,max(curr_im(:))];
    
    AP_movie2avi(curr_im,movie_rate,color_map,color_axis, ...
        figure_position,savefile,t_annotation,movie_annotation);
    disp(['Saved ' savefile]);
end


%% Movie 6: Cortex->striatum kernels

protocols = {'vanillaChoiceworld'};

for protocol = protocols 
    protocol = cell2mat(protocol);
    
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
    k_fn = [data_path filesep 'ctx_str_kernels_' protocol];
    load(k_fn);
    
    framerate = 35;
    upsample_factor = 1;
    sample_rate = framerate*upsample_factor;
    kernel_t = [-0.5,0.5];
    kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
    t = kernel_frames/sample_rate;
    
    % Concatenate and mean
    % (kernel is -:+ fluorescence lag, flip to be spike-oriented)
    k_px_timeflipped = cellfun(@(x) cellfun(@(x) x(:,:,end:-1:1,:),x,'uni',false),ctx_str_kernel,'uni',false);
    k_px_animal = cellfun(@(x) nanmean(cat(5,x{:}),5),k_px_timeflipped,'uni',false);
    k_px = nanmean(double(cat(5,k_px_animal{:})),5);
    
    % Normalize to max weight
    k_px_norm = k_px./max(max(max(k_px,[],1),[],2),[],3);
    
    % Make movie of ctx-str kernels by depth
    movie_rate = sample_rate/5;
    color_map = brewermap([],'*RdBu');
    color_axis = [-1,1];
    figure_position = [24,408,1868,399];
    save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\figs\movies';
    savefile = [save_path filesep 'movie_6_ctx_str_kernel'];
    t_annotation = cellfun(@(x) sprintf('Lag from multiunit: %0.2f sec',x), ...
        num2cell(t),'uni',false);    
    movie_annotation = {'Medial domain','Centromedial domain','Centrolateral domain','Lateral domain'};
    
    AP_movie2avi(k_px_norm,movie_rate,color_map,color_axis,figure_position,savefile,t_annotation,movie_annotation);
  
%     % (Colored: not done at the moment: requires a colored AP_movie2avi)
%     % Get center-of-mass maps
%     k_px_positive = k_px;
%     k_px_positive(k_px_positive < 0) = 0;
%     k_px_com = sum(k_px_positive.*permute(1:n_aligned_depths,[1,3,4,2]),4)./sum(k_px_positive,4);
%     k_px_com_colored = nan(size(k_px_com,1),size(k_px_com,2),3,size(k_px_com,3));
%     
%     use_colormap = min(jet(255)-0.2,1);
%     for curr_frame = 1:size(k_px_com,3)
%         k_px_com_colored(:,:,:,curr_frame) = ...
%             ind2rgb(round(mat2gray(k_px_com(:,:,curr_frame),...
%             [1,n_aligned_depths])*size(use_colormap,1)),use_colormap);
%     end
%            
%     % Whiten colored kernel by transparency
%     k_px_com_colored_transparent = k_px_com_colored + ...
%         ((1-(permute(max(abs(k_px),[],4)./ ...
%         prctile(abs(k_px(:)),100),[1,2,4,3]))) .* ...
%         (ones(1,1,3,length(t)) - k_px_com_colored));
%     
%     % Make movie of colored kernel
%     movie_rate = sample_rate/3;
%     color_map = brewermap([],'PRGn');
%     color_axis = [-0.015,0.015];
%     figure_position = [651,371,574,450];
%     save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\figs\movies';
%             
%     savefile = [save_path filesep '_ctx_str_colored'];
%     annotation_text = cellfun(@(x) sprintf('Time from event: %0.2f sec',x), ...
%         num2cell(task_regressor_t_shifts{curr_regressor}),'uni',false);    
% %     AP_movie2avi(k_px_com_colored_transparent,movie_rate,color_map,color_axis,figure_position,savefile,annotation_text,n_subregressors)
       
end





%% ~~~~~~~~~~~~~ UNUSED ~~~~~~~~~~~

%% Spatial map of explained variance in cortex 

use_t = true(size(t));
spatial_downsample = 10;

correct_trials = trial_outcome_allcat == 1;
incorrect_trials = trial_outcome_allcat == -1;

visual_trials = trial_stim_allcat ~= 0;
zero_trials = trial_stim_allcat == 0;

ctx_expl_var_correct = AP_spatial_explained_var(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv(correct_trials,use_t,:),[2,1,3]),[],n_vs)', ...
    reshape(permute(fluor_taskpred_allcat(correct_trials,use_t,:),[2,1,3]),[],n_vs)',spatial_downsample);
ctx_expl_var_incorrect = AP_spatial_explained_var(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv(incorrect_trials,use_t,:),[2,1,3]),[],n_vs)', ...
    reshape(permute(fluor_taskpred_allcat(incorrect_trials,use_t,:),[2,1,3]),[],n_vs)',spatial_downsample);

ctx_expl_var_visual = AP_spatial_explained_var(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv(visual_trials,use_t,:),[2,1,3]),[],n_vs)', ...
    reshape(permute(fluor_taskpred_allcat(visual_trials,use_t,:),[2,1,3]),[],n_vs)',spatial_downsample);
ctx_expl_var_zero = AP_spatial_explained_var(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv(zero_trials,use_t,:),[2,1,3]),[],n_vs)', ...
    reshape(permute(fluor_taskpred_allcat(zero_trials,use_t,:),[2,1,3]),[],n_vs)',spatial_downsample);

figure;
subplot(2,3,1);
imagesc(ctx_expl_var_correct);
axis image off; 
colormap(brewermap([],'*RdBu'))
caxis([-1,1]); c = colorbar; ylabel(c,'Variance explained');
AP_reference_outline('ccf_aligned','k');
title('Correct trials');

subplot(2,3,2);
imagesc(ctx_expl_var_incorrect);
axis image off; 
colormap(brewermap([],'*RdBu'))
caxis([-1,1]); c = colorbar; ylabel(c,'Variance explained');
AP_reference_outline('ccf_aligned','k');
title('Incorrect trials');

subplot(2,3,3);
imagesc(ctx_expl_var_correct - ctx_expl_var_incorrect);
axis image off; 
colormap(brewermap([],'*RdBu'))
caxis([-1,1]); c = colorbar; ylabel(c,'Variance explained');
AP_reference_outline('ccf_aligned','k');
title('Correct - incorrect');

subplot(2,3,4);
imagesc(ctx_expl_var_visual);
axis image off; 
colormap(brewermap([],'*RdBu'))
caxis([-1,1]); c = colorbar; ylabel(c,'Variance explained');
AP_reference_outline('ccf_aligned','k');
title('Visual trials');

subplot(2,3,5);
imagesc(ctx_expl_var_zero);
axis image off; 
colormap(brewermap([],'*RdBu'))
caxis([-1,1]); c = colorbar; ylabel(c,'Variance explained');
AP_reference_outline('ccf_aligned','k');
title('Zero trials');

subplot(2,3,6);
imagesc(ctx_expl_var_visual - ctx_expl_var_zero);
axis image off; 
colormap(brewermap([],'*RdBu'))
caxis([-1,1]); c = colorbar; ylabel(c,'Variance explained');
AP_reference_outline('ccf_aligned','k');
title('Visual - zero');


%% Compare kernels with domain maps

%%% Task>cortex regression results

% Get task>cortex parameters
n_regressors = 4;
t_shifts = {[0,0.5]; ... % stim
    [-0.5,1]; ... % move
    [-0.1,0.5]; ... % go cue
    [-0.5,1]}; % outcome
regressor_labels = {'Stim','Move onset','Go cue','Outcome'};

sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
    round(x(2)*(sample_rate)),t_shifts,'uni',false);
t_shifts = cellfun(@(x) x/sample_rate,sample_shifts,'uni',false);

% Get average task>striatum kernels
regressor_px = cell(n_regressors,1);
for curr_regressor = 1:n_regressors
    curr_k_cell = cellfun(@(x) x(curr_regressor,:),vertcat(fluor_taskpred_k_all{:}),'uni',false);
    curr_k_cell = vertcat(curr_k_cell{:});
    curr_k = permute(cell2mat(permute(arrayfun(@(x) ...
        nanmean(cat(3,curr_k_cell{:,x}),3),1:n_vs,'uni',false),[1,3,2])),[3,2,1]);
    curr_k_px = cell2mat(permute(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:size(curr_k,1)), ...
        curr_k(:,:,x)),1:size(curr_k,3),'uni',false),[1,3,4,2]));
    AP_imscroll(curr_k_px,t_shifts{curr_regressor});
    axis image;
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'*RdBu'));
    AP_reference_outline('ccf_aligned','k');
    set(gcf,'Name',regressor_labels{curr_regressor});
    
    regressor_px{curr_regressor} = curr_k_px;
end

% Plot max for each kernel across time
regressor_t_max = cellfun(@(x) squeeze(max(x,[],3)),regressor_px,'uni',false);

figure;
max_subregressors = max(cellfun(@(x) size(x,3),regressor_t_max));
max_c = max(abs(cell2mat(cellfun(@(x) x(:),regressor_t_max,'uni',false))));
for curr_regressor = 1:n_regressors
    for curr_subregressor = 1:size(regressor_t_max{curr_regressor},3)
        subplot(n_regressors,max_subregressors, ...
            curr_subregressor+(curr_regressor-1)*max_subregressors);
        imagesc(regressor_t_max{curr_regressor}(:,:,curr_subregressor));
        AP_reference_outline('ccf_aligned','k');
        axis image off; 
        colormap(brewermap([],'PRGn'));
        caxis([-max_c,max_c]);
    end
end





% Get average cortex->striatum kernel
ctx_str_k_mean = nanmean(cell2mat(permute(vertcat(ctx_str_k_all{:}),[2,3,4,5,1])),5);
ctx_str_k_mean_px = cell2mat(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:100), ...
    ctx_str_k_mean(:,:,x)),permute(1:n_depths,[1,3,4,2]),'uni',false));

ctx_str_kernel_frames_t = [-0.5,0.5];
ctx_str_kernel_frames = round(ctx_str_kernel_frames_t(1)*sample_rate): ...
    round(ctx_str_kernel_frames_t(2)*sample_rate);
ctx_str_kernel_t = ctx_str_kernel_frames./sample_rate;

% (use t = 0)
ctx_str_kernel_t0 = squeeze(ctx_str_k_mean_px(:,:,ctx_str_kernel_t == 0,:));

% Get correlation between task->cortex and cortex->striatum regressor
ctx_task_str_corr = cellfun(@(k) cell2mat(arrayfun(@(x) 1 - pdist2(reshape(k(:,:,x),[],1)', ...
    reshape(ctx_str_kernel_t0,[],n_depths)','correlation'), ...
    transpose(1:size(k,3)),'uni',false)),regressor_t_max,'uni',false); 

figure;
for curr_regressor = 1:n_regressors
    if curr_regressor == 1
        col = colormap_BlueWhiteRed(size(ctx_task_str_corr{curr_regressor},1)/2);
        col(median(1:size(col,1)),:) = [];
    else
        col = lines(size(ctx_task_str_corr{curr_regressor}));
    end
    
    subplot(n_regressors,1,curr_regressor); hold on;
    set(gca,'ColorOrder',col);
    bar(ctx_task_str_corr{curr_regressor}');    
    set(gca,'XTick',1:4);
    xlabel('Striatum domain')
    ylabel('Kernel correlation')
    title(regressor_labels{curr_regressor});
end

%% Example data with task regression overlay

animal = 'AP025'; 
day = '2017-10-01'; 
experiment = 1; 
str_align = 'none'; 
verbose = true; 
AP_load_experiment;

% Align U deconvolve V, set components to use
use_components = 1:200;
aUdf = AP_align_widefield(animal,day,Udf);
fVdf_deconv = AP_deconv_wf(fVdf);

% Regress task events to fluorescence
% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

% Get event-aligned activity
raster_window = [-0.5,2];
upsample_factor = 1;
raster_sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):raster_sample_rate:raster_window(2);

% Get align times
use_align = stimOn_times;
use_align(isnan(use_align)) = 0;

t_peri_event = bsxfun(@plus,use_align,t);
t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];

%%% Trial-align wheel velocity
event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
    wheel_velocity,t_peri_event);

%%% Trial-align facecam movement
event_aligned_movement = interp1(facecam_t(~isnan(facecam_t)), ...
    frame_movement(~isnan(facecam_t)),t_peri_event);

%%% Trial-align outcome (reward page 1, punish page 2)
% (note incorrect outcome imprecise from signals, but looks good)
event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);

event_aligned_outcome(trial_outcome == 1,:,1) = ...
    (cell2mat(arrayfun(@(x) ...
    histcounts(reward_t_timeline,t_bins(x,:)), ...
    find(trial_outcome == 1),'uni',false))) > 0;

event_aligned_outcome(trial_outcome == -1,:,2) = ...
    (cell2mat(arrayfun(@(x) ...
    histcounts(signals_events.responseTimes,t_bins(x,:)), ...
    find(trial_outcome == -1),'uni',false))) > 0;

% Pick trials to keep
use_trials = ...
    trial_outcome ~= 0 & ...
    ~signals_events.repeatTrialValues(1:n_trials)' & ...
    stim_to_feedback < 1.5;

% Get behavioural data
D = struct;
D.stimulus = zeros(sum(use_trials),2);

L_trials = signals_events.trialSideValues(1:n_trials) == -1;
R_trials = signals_events.trialSideValues(1:n_trials) == 1;

D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials');
D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials');

D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
D.repeatNum = ones(sum(use_trials),1);

D.outcome = reshape(trial_outcome(use_trials),[],1);

%%% Regress task to cortex/striatum/cortex-predicted striatum

% Get time points to bin
sample_rate = framerate*regression_params.upsample_factor;
time_bins = frame_t(find(frame_t > ...
    regression_params.skip_seconds,1)):1/sample_rate: ...
    frame_t(find(frame_t-frame_t(end) < ...
    -regression_params.skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Get reaction time for building regressors
[move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
move_idx(~move_trial) = NaN;
move_t = nan(size(move_idx));
move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';

% Stim regressors
unique_stim = unique(contrasts(contrasts > 0).*sides');
stim_contrastsides = ...
    signals_events.trialSideValues(1:length(stimOn_times))'.* ...
    signals_events.trialContrastValues(1:length(stimOn_times))';

stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
for curr_stim = 1:length(unique_stim)
    curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
    stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
end

% Stim move regressors (one for each stim when it starts to move)
stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
for curr_stim = 1:length(unique_stim)
    
    % (find the first photodiode flip after the stim azimuth has
    % moved past a threshold)
    
    curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
        stim_contrastsides == unique_stim(curr_stim));
    
    azimuth_move_threshold = 5; % degrees to consider stim moved
    stim_move_times_signals = ...
        signals_events.stimAzimuthTimes( ...
        abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
    curr_stim_move_times_signals = arrayfun(@(x) ...
        stim_move_times_signals(find(stim_move_times_signals > ...
        curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
    
    curr_stim_move_times_photodiode = arrayfun(@(x) ...
        photodiode_flip_times(find(photodiode_flip_times > ...
        curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));
    
    stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
    
end

% Stim center regressors (one for each stim when it's stopped during reward)
unique_contrasts = unique(contrasts(contrasts > 0));

stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
for curr_contrast = 1:length(unique_contrasts)
    
    % (find the last photodiode flip before the reward)
    curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ...
        abs(stim_contrastsides) == unique_contrasts(curr_contrast));
    
    curr_reward_times = arrayfun(@(x) ...
        reward_t_timeline(find(reward_t_timeline > ...
        curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
    
    curr_prereward_photodiode_times = arrayfun(@(x) ...
        photodiode_flip_times(find(photodiode_flip_times < ...
        curr_reward_times(x),1,'last')),1:length(curr_reward_times));
    
    stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
    
end

% Move onset regressors (L/R)
move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
    find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == -1));
move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
    find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == 1));

move_onset_regressors = zeros(2,length(time_bin_centers));
move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);

% Move onset x stim regressors (one for each contrast/side)
move_onset_stim_time_absolute = arrayfun(@(curr_stim) ...
    arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
    find(~isnan(move_idx) & stim_contrastsides == unique_stim(curr_stim))), ...
    1:length(unique_stim),'uni',false);

move_onset_stim_regressors = zeros(10,length(time_bin_centers));
for curr_stim = 1:length(unique_stim)
    move_onset_stim_regressors(curr_stim,:) = ...
        histcounts(move_onset_stim_time_absolute{curr_stim},time_bins);
end

% Move ongoing regressors (L/R choice for duration of movement)
wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);

move_stopped_t = 0.5;
move_stopped_samples = round(sample_rate*move_stopped_t);
wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.02), ...
    ones(1,move_stopped_samples),'full') > 0;
wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);

move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
    find(time_bin_centers > x): ...
    find(time_bin_centers > x & ~wheel_moving,1), ...
    move_time_L_absolute','uni',false));
move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
    find(time_bin_centers > x): ...
    find(time_bin_centers > x & ~wheel_moving,1), ...
    move_time_R_absolute','uni',false));

move_ongoing_regressors = zeros(2,length(time_bin_centers));
move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
move_ongoing_regressors(2,move_ongoing_R_samples) = 1;

% Go cue regressors - separate for early/late move
% (using signals timing - not precise but looks good)
if length(signals_events.interactiveOnTimes) ~= length(move_t)
    error('Different number of interactive ons and move times')
end

go_cue_regressors = zeros(2,length(time_bin_centers));
go_cue_regressors(1,:) = histcounts( ...
    signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
go_cue_regressors(2,:) = histcounts( ...
    signals_events.interactiveOnTimes(move_t > 0.5),time_bins);

% Outcome regressors
% (using signals timing - not precise but looks good)
outcome_regressors = zeros(2,length(time_bin_centers));

outcome_regressors(1,:) = histcounts( ...
    reward_t_timeline,time_bins);
outcome_regressors(2,:) = histcounts( ...
    signals_events.responseTimes(trial_outcome == -1),time_bins);

% Concatenate selected regressors, set parameters
regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
regressor_labels = {'Stim','Move onset','Go cue','Outcome'};

t_shifts = {[0,0.5]; ... % stim
    [-0.5,1]; ... % move
    [-0.1,0.5]; ... % go cue
    [-0.5,1]}; % outcome

sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
    round(x(2)*(sample_rate)),t_shifts,'uni',false);
lambda = 0;
zs = [false,false];
cvfold = 5;
use_constant = false;
return_constant = false;

% Regress task -> fluor
event_aligned_V_deconv = ...
    interp1(frame_t,fVdf_deconv(use_components,:)',t_peri_event);

baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,:),[],size(event_aligned_V_deconv,3)))';
activity = interp1(frame_t,fVdf_deconv(use_components,:)',time_bin_centers)' - baseline;

[~,fluor_taskpred_short,~,~] = ...
    AP_regresskernel(regressors,activity,sample_shifts, ...
    lambda,zs,cvfold,return_constant,use_constant);

% (interpolate the task-predicted fluorescence and add back baseline)
fluor_taskpred = interp1(time_bin_centers,fluor_taskpred_short',frame_t)' + baseline;

%%% Plot example data

% Set time to plot
plot_t = [134,152];

raster_fig = figure;

% (wheel velocity)
wheel_axes = subplot(6,1,1);
plot_wheel_idx = Timeline.rawDAQTimestamps >= plot_t(1) & ...
    Timeline.rawDAQTimestamps <= plot_t(2);
plot(wheel_axes,Timeline.rawDAQTimestamps(plot_wheel_idx), ...
    wheel_velocity(plot_wheel_idx),'k','linewidth',2);
ylabel('Wheel velocity');

% (stimuli)
stim_col = colormap_BlueWhiteRed(5);
[~,trial_contrast_idx] = ...
    ismember(trial_conditions(:,1).*trial_conditions(:,2),unique(contrasts'.*sides),'rows');
stim_lines = arrayfun(@(x) line(wheel_axes,repmat(stimOn_times(x),1,2),ylim(wheel_axes),'color', ...
    stim_col(trial_contrast_idx(x),:),'linewidth',2), ...
    find(stimOn_times >= plot_t(1) & stimOn_times <= plot_t(2)));

% (movement starts)
move_col = [0.6,0,0.6;0,0.6,0];
[~,trial_choice_idx] = ismember(trial_conditions(:,3),[-1;1],'rows');
move_lines = arrayfun(@(x) line(wheel_axes,repmat(wheel_move_time(x),1,2),ylim(wheel_axes),'color', ...
    move_col(trial_choice_idx(x),:),'linewidth',2), ...
    find(wheel_move_time >= plot_t(1) & wheel_move_time <= plot_t(2)));

% (go cues)
go_col = [0.8,0.8,0.2];
go_cue_times = signals_events.interactiveOnTimes(1:n_trials);
go_cue_lines = arrayfun(@(x) line(wheel_axes,repmat(go_cue_times(x),1,2),ylim(wheel_axes),'color', ...
    go_col,'linewidth',2,'linestyle','--'), ...
    find(go_cue_times >= plot_t(1) & go_cue_times <= plot_t(2)));

% (outcomes)
outcome_col = [0,0,0.8;0.5,0.5,0.5];
reward_lines = arrayfun(@(x) line(wheel_axes,repmat(reward_t_timeline(x),1,2),ylim(wheel_axes),'color', ...
    outcome_col(1,:),'linewidth',2,'linestyle','--'), ...
    find(reward_t_timeline >= plot_t(1) & reward_t_timeline <= plot_t(2)));
punish_times = signals_events.responseTimes(trial_outcome == -1);
punish_lines = arrayfun(@(x) line(wheel_axes,repmat(punish_times(x),1,2),ylim(wheel_axes),'color', ...
    outcome_col(2,:),'linewidth',2,'linestyle','--'), ...
    find(punish_times >= plot_t(1) & punish_times <= plot_t(2)));

% (striatum raster)
raster_axes = subplot(6,1,2:4,'YDir','reverse'); hold on;
plot_spikes = spike_times_timeline >= plot_t(1) & ...
    spike_times_timeline <= plot_t(2) & ...
    spike_depths >= str_depth(1) & spike_depths <= str_depth(2);
plot(raster_axes,spike_times_timeline(plot_spikes),spike_depths(plot_spikes),'.k');
ylabel('Depth (\mum)');
xlabel('Time (s)')

% (fluorescence from select ROIs)
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
roi_trace = AP_svd_roi(aUdf(:,:,use_components),fVdf_deconv(use_components,:),[],[],cat(3,wf_roi.mask));
roi_trace_taskpred = AP_svd_roi(aUdf(:,:,use_components),fluor_taskpred,[],[],cat(3,wf_roi.mask));

plot_rois = [1,7,9];
fluor_spacing = 70;
fluor_axes = subplot(6,1,5:6); hold on;
plot_fluor_idx = frame_t >= plot_t(1) & frame_t <= plot_t(2);
AP_stackplot(roi_trace(plot_rois,plot_fluor_idx)', ...
    frame_t(plot_fluor_idx),fluor_spacing,false,[0,0.7,0],{wf_roi(plot_rois).area});
AP_stackplot(roi_trace_taskpred(plot_rois,plot_fluor_idx)', ...
    frame_t(plot_fluor_idx),fluor_spacing,false,'b',{wf_roi(plot_rois).area});


linkaxes([wheel_axes,raster_axes,fluor_axes],'x');


% % Write legend
% [~,unique_contrasts_h] = unique(trial_contrast_idx);
% [~,unique_move_h] = unique(trial_choice_idx(trial_choice_idx > 0));
% legend([stim_lines(unique_contrasts_h),move_lines(unique_move_h), ...
%     go_cue_lines(1),reward_lines(1),punish_lines(1)], ...
%     [cellfun(@(x) ['Stim ' num2str(x)],num2cell(unique(contrasts'.*sides)),'uni',false); ...
%     {'Move L';'Move R';'Go cue';'Reward';'Punish'}]);


% Plot fluorescence at regular intervals within time range

% plot_frames_idx = [4721,4833,4924,5008]; % (across 4 trials)
plot_frames_idx = [4827,4831,4834,4838]; % (within one trial)

plot_frames = svdFrameReconstruct(aUdf(:,:,use_components),fVdf_deconv(use_components,plot_frames_idx));
plot_frames_taskpred = svdFrameReconstruct(aUdf(:,:,use_components),fluor_taskpred(:,plot_frames_idx));

wf_fig = figure;
for curr_frame = 1:length(plot_frames_idx)
    subplot(2,length(plot_frames_idx),curr_frame);
    imagesc(plot_frames(:,:,curr_frame));
    colormap(brewermap([],'PrGn'));
    caxis([-0.03,0.03]);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title(sprintf('Deconv - %0.2fs',frame_t(plot_frames_idx(curr_frame))));
    axis image off;
    
    subplot(2,length(plot_frames_idx),length(plot_frames_idx) + curr_frame);
    imagesc(plot_frames_taskpred(:,:,curr_frame));
    colormap(brewermap([],'PrGn'));
    caxis([-0.03,0.03]);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title(sprintf('Task predicted - %0.2fs',frame_t(plot_frames_idx(curr_frame))));
    axis image off;
end

% (draw lines on the ROI plot where the frames were taken from)
frame_lines = arrayfun(@(x) line(fluor_axes,repmat(frame_t(x),1,2), ...
    ylim(fluor_axes),'color','k','linewidth',2),plot_frames_idx);


% Plot ROIs
figure; hold on
set(gca,'YDir','reverse');
AP_reference_outline('ccf_aligned','k');
for curr_roi = plot_rois
    curr_roi_boundary = cell2mat(bwboundaries(wf_roi(curr_roi).mask()));
    patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),[0,0.8,0]);   
    text(nanmean(curr_roi_boundary(:,2)),nanmean(curr_roi_boundary(:,1)), ...
        wf_roi(curr_roi).area,'FontSize',12,'HorizontalAlignment','center')
end
axis image off;


%% Striatal unit regression examples

% (animal 2 day 3, AP025 2017-10-01, run unit regression first)

smooth_size = 11;
gw = gausswin(smooth_size,5)';
smWin = gw./sum(gw);
% smWin = 1;

binned_spikes_smoothed = conv2(binned_spikes,smWin,'same');
binned_spikes_taskpred_smoothed = conv2(binned_spikes_taskpred,smWin,'same');

figure;
example_units = [237,150,170,66];
p = nan(length(example_units),1);
for curr_unit_idx = 1:length(example_units)
    curr_unit = example_units(curr_unit_idx);
    
    p(curr_unit_idx) = subplot(length(example_units),1,curr_unit_idx);
    hold on;
    plot(binned_spikes_smoothed(curr_unit,:),'k');
    plot(binned_spikes_taskpred_smoothed(curr_unit,:),'r');
    title(curr_unit);
    xlabel('Time');
    ylabel('Spikes');
end

linkaxes(p,'x');

x_bounds = [10877,12010];
xlim(p(1),x_bounds);

%% Task > striatum units explained variance histogram

% Load the unit kernel results
unit_kernel_dir = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
unit_kernel_fn = 'unit_kernel_all.mat';
load([unit_kernel_dir filesep unit_kernel_fn]);
task_regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
task_regressor_cols = [1,0,0;1,0,1;0.8,0.7,0.2;0,0,1];

% Get estimation of end of striatum for each recording
ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
ephys_align_fn = ['ephys_depth_align.mat'];
load([ephys_align_path filesep ephys_align_fn]);

% Concatenate all units and re-depth by end of striatum

% (get estimation of end of striatum for each recording)
ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
ephys_align_fn = ['ephys_depth_align.mat'];
load([ephys_align_path filesep ephys_align_fn]);
str_depth_cat = vertcat(ephys_depth_align(1:6).str_depth);

animaldays = reshape(arrayfun(@(x) ~isempty(unit_kernel_all(x).template_depths), ...
    1:numel(unit_kernel_all)),size(unit_kernel_all));

str_units = cellfun(@(unit_depths,str_depths) ...
    unit_depths >= str_depths(1) & unit_depths <= str_depths(2), ...
    AP_index_ans(reshape({unit_kernel_all.template_depths}, ...
    size(unit_kernel_all))',animaldays'), ...
    mat2cell(str_depth_cat,ones(size(str_depth_cat,1),1),2),'uni',false);

template_depths_cat = cellfun(@(unit_depths,str_depths,str_units) ...
    unit_depths(str_units) - str_depths(2), ...
    AP_index_ans(reshape({unit_kernel_all.template_depths}, ...
    size(unit_kernel_all))',animaldays'), ...
    mat2cell(str_depth_cat,ones(size(str_depth_cat,1),1),2), ...
    str_units,'uni',false);
unit_expl_var_total_cat = cellfun(@(x,str_units) x(str_units), ...
    AP_index_ans(reshape({unit_kernel_all.unit_expl_var_total}, ...
    size(unit_kernel_all))',animaldays'),str_units,'uni',false);
unit_expl_var_partial_cat = cellfun(@(x,str_units) x(str_units,:,:), ...
    AP_index_ans(reshape({unit_kernel_all.unit_expl_var_partial}, ...
    size(unit_kernel_all))',animaldays'),str_units,'uni',false);
spike_rate_cat = cellfun(@(x,str_units) x(str_units), ...
    AP_index_ans(reshape({unit_kernel_all.spike_rate}, ...
    size(unit_kernel_all))',animaldays'),str_units,'uni',false);

% Plot all units explained by kernels and expl var histogram by depth
use_partial = 2;
rate_cutoff = 0.1;
use_units_rate = cell2mat(spike_rate_cat) > rate_cutoff;

template_depths_allcat = cell2mat(template_depths_cat);

depth_bin_um = 400;
depth_bins = linspace(min(template_depths_allcat), ...
    max(template_depths_allcat),round(range(template_depths_allcat)/depth_bin_um));
depth_bin_centers = depth_bins(1:end-1) + diff(depth_bins)./2;
unit_depth_bins = discretize(template_depths_allcat,depth_bins);

unit_expl_var_partial_allcat = cell2mat(cellfun(@(x) x(:,:,use_partial), ...
    unit_expl_var_partial_cat,'uni',false));

figure; 
for curr_regressor = 1:length(task_regressor_labels)
    subplot(1,length(task_regressor_labels)+1,curr_regressor);
    hold on; set(gca,'YDir','reverse');
    plot_units = use_units_rate & unit_expl_var_partial_allcat(:,curr_regressor) > 0;
    scatter(unit_expl_var_partial_allcat(plot_units,curr_regressor), ...
        template_depths_allcat(plot_units),10,task_regressor_cols(curr_regressor,:),'filled');
    xlabel('Expl var');
    ylabel('Depth (\mum)');
    
    subplot(1,length(task_regressor_labels)+1,length(task_regressor_labels)+1);   
    hold on; set(gca,'YDir','reverse');
    curr_expl_var_norm = unit_expl_var_partial_allcat(:,curr_regressor);
    curr_expl_var_norm(curr_expl_var_norm < 0) = 0;
    curr_expl_var_norm_depth = accumarray(unit_depth_bins(plot_units), ...
        curr_expl_var_norm(plot_units,:), ...
        [length(depth_bin_centers),1],@nanmean,NaN);   
    plot(curr_expl_var_norm_depth,depth_bin_centers, ...
        'linewidth',2,'color',task_regressor_cols(curr_regressor,:));
    xlabel('Mean expl var (> 0)');
    ylabel('Depth (\mum)');
end

% (draw striatum starts)
for curr_plot = 1:length(task_regressor_labels)+1
    subplot(1,length(task_regressor_labels)+1,curr_plot)
    str_start_line = nan(size(str_depth_cat,1),1);
    for i = 1:size(str_depth_cat,1)
        str_start_line(i) = line(xlim,repmat(-diff(str_depth_cat(i,:)),1,2),'color',[0.8,0.8,0.8]);
    end
    set(gca,'Children',circshift(get(gca,'Children'),-size(str_depth_cat,1)));
end

%% TESTING supplemental task>ctx-str

% Get task>striatum parameters
n_regressors = length(task_regressor_labels);

% Normalize task > striatum kernels across experiments with mua_norm
mua_taskpred_k_all_norm = cellfun(@(kernel_animal,mua_norm_animal) ...
    cellfun(@(kernel_exp,mua_norm_exp) ...
    cellfun(@(kernel_regressor) ...
    kernel_regressor./(mua_norm_exp/sample_rate), ...
    kernel_exp,'uni',false),kernel_animal,mua_norm_animal,'uni',false), ...
    mua_taskpred_k_all,mua_norm,'uni',false);

mua_ctxpred_taskpred_k_all_norm = cellfun(@(kernel_animal,mua_norm_animal) ...
    cellfun(@(kernel_exp,mua_norm_exp) ...
    cellfun(@(kernel_regressor) ...
    kernel_regressor./(mua_norm_exp/sample_rate), ...
    kernel_exp,'uni',false),kernel_animal,mua_norm_animal,'uni',false), ...
    mua_ctxpred_taskpred_k_all,mua_norm,'uni',false);

% Average and concatenate task>striatum kernels within animals
task_str_k_animal = cell(n_regressors,1);
task_ctxpred_str_k_animal = cell(n_regressors,1);
for curr_animal = 1:length(mua_taskpred_k_all_norm)
    if isempty(mua_taskpred_k_all_norm{curr_animal})
        continue
    end
    
    curr_k = cat(2,mua_taskpred_k_all_norm{curr_animal}{:});
    curr_ctxpred_k = cat(2,mua_ctxpred_taskpred_k_all_norm{curr_animal}{:});
    for curr_regressor = 1:n_regressors
        curr_k_mean = nanmean(cat(4,curr_k{curr_regressor,:}),4);        
        task_str_k_animal{curr_regressor} = cat(4, ...
            task_str_k_animal{curr_regressor},curr_k_mean);
        
        curr_ctxpred_k_mean = nanmean(cat(4,curr_ctxpred_k{curr_regressor,:}),4);        
        task_ctxpred_str_k_animal{curr_regressor} = cat(4, ...
            task_ctxpred_str_k_animal{curr_regressor},curr_ctxpred_k_mean);
    end
end

% Plot task>striatum kernels
stim_col = colormap_BlueWhiteRed(5);
move_col = [0.6,0,0.6;0,0.6,0];
go_col = [0.8,0.8,0.2;0.5,0.5,0.5];
outcome_col = [0,0,0.8;0.8,0,0];
task_regressor_cols = {stim_col,move_col,go_col,outcome_col};
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

figure;
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors  
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);  
        
        curr_kernels = permute(task_ctxpred_str_k_animal{curr_regressor}(:,:,curr_depth,:),[1,2,4,3]);
        n_subregressors = size(task_ctxpred_str_k_animal{curr_regressor},1);
        col = task_regressor_cols{curr_regressor};
        for curr_subregressor = 1:n_subregressors
            AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
                nanmean(curr_kernels(curr_subregressor,:,:),3), ...
                AP_sem(curr_kernels(curr_subregressor,:,:),3), ...
                col(curr_subregressor,:),0.5);
        end
        
        xlabel('Time (s)');
        ylabel('Weight');
        title(task_regressor_labels{curr_regressor});
        line([0,0],ylim,'color','k');
        
    end
end
linkaxes(p);

%% Fig 4x: average measured/task-/cortex-predicted activity by condition

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_stim_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set windows to average activity
align_times = {stim_align,move_align,outcome_align};
align_labels = {'Stim','Move onset','Outcome'};
align_reduction = [1,2,4];
% align_trials = { ...
%     [move_t < 0.5 & trial_stim_allcat > 0], ...
%     [move_t < 0.5 & trial_choice_allcat == -1], ...
%     [move_t < 0.5 & trial_outcome_allcat == 1]};
align_trials = { ...
    [move_t < 0.5 & trial_stim_allcat < 0], ...
    [move_t < 0.5 & trial_choice_allcat == 1], ...
    [move_t < 0.5 & trial_outcome_allcat == -1]};

figure;
for curr_align = 1:length(align_times)      
    
    % Re-align and split activity
    curr_act = mat2cell(...
        cell2mat(arrayfun(@(trial) circshift( ...
        mua_allcat(trial,:,:), ...
        align_times{curr_align}(trial),2),transpose(1:size(mua_allcat,1)),'uni',false)), ...
        use_split,length(t),n_depths);
    
    curr_act_taskpred = mat2cell(...
        cell2mat(arrayfun(@(trial) circshift( ...
        mua_taskpred_allcat(trial,:,:), ...
        align_times{curr_align}(trial),2),transpose(1:size(mua_allcat,1)),'uni',false)), ...
        use_split,length(t),n_depths);
    
    curr_act_ctxpred = mat2cell(...
        cell2mat(arrayfun(@(trial) circshift( ...
        mua_ctxpred_allcat(trial,:,:), ...
        align_times{curr_align}(trial),2),transpose(1:size(mua_allcat,1)),'uni',false)), ...
        use_split,length(t),n_depths);
    
    % Re-align and split activity (reduced)
    curr_act_reduced = mat2cell(...
        cell2mat(arrayfun(@(trial) circshift( ...
        mua_allcat(trial,:,:) - mua_taskpred_reduced_allcat(trial,:,:,align_reduction(curr_align)), ...
        align_times{curr_align}(trial),2),transpose(1:size(mua_allcat,1)),'uni',false)), ...
        use_split,length(t),n_depths);
    
    curr_act_taskpred_reduced = mat2cell(...
        cell2mat(arrayfun(@(trial) circshift( ...
        mua_taskpred_allcat(trial,:,:) - mua_taskpred_reduced_allcat(trial,:,:,align_reduction(curr_align)), ...
        align_times{curr_align}(trial),2),transpose(1:size(mua_allcat,1)),'uni',false)), ...
        use_split,length(t),n_depths);
    
    curr_act_ctxpred_reduced = mat2cell(...
        cell2mat(arrayfun(@(trial) circshift( ...
        mua_ctxpred_allcat(trial,:,:) - mua_ctxpred_taskpred_reduced_allcat(trial,:,:,align_reduction(curr_align)), ...
        align_times{curr_align}(trial),2),transpose(1:size(mua_allcat,1)),'uni',false)), ...
        use_split,length(t),n_depths);
  
    % Get average for chosen trials
    curr_trials = mat2cell(align_trials{curr_align},use_split);
    
    curr_act_trialmean = cell2mat(cellfun(@(act,trials) ...
        nanmean(act(trials,:,:),1),curr_act,curr_trials,'uni',false));
    curr_act_taskpred_trialmean = cell2mat(cellfun(@(act,trials) ...
        nanmean(act(trials,:,:),1),curr_act_taskpred,curr_trials,'uni',false));
    curr_act_ctxpred_trialmean = cell2mat(cellfun(@(act,trials) ...
        nanmean(act(trials,:,:),1),curr_act_ctxpred,curr_trials,'uni',false));
    
    curr_act_reduced_trialmean = cell2mat(cellfun(@(act,trials) ...
        nanmean(act(trials,:,:),1),curr_act_reduced,curr_trials,'uni',false));
    curr_act_taskpred_reduced_trialmean = cell2mat(cellfun(@(act,trials) ...
        nanmean(act(trials,:,:),1),curr_act_taskpred_reduced,curr_trials,'uni',false));
    curr_act_ctxpred_reduced_trialmean = cell2mat(cellfun(@(act,trials) ...
        nanmean(act(trials,:,:),1),curr_act_ctxpred_reduced,curr_trials,'uni',false));
    
    % Plot
    for curr_depth = 1:n_depths
        subplot(n_depths,length(align_times), ...
            sub2ind([length(align_times),n_depths],curr_align,curr_depth));
        hold on
               
        plot(t,nanmean(curr_act_trialmean(:,:,curr_depth),1),'--k','linewidth',1);
        plot(t,nanmean(curr_act_reduced_trialmean(:,:,curr_depth),1),'k','linewidth',2);
        AP_errorfill(t,nanmean(curr_act_taskpred_reduced_trialmean(:,:,curr_depth),1), ...
            AP_sem(curr_act_taskpred_reduced_trialmean(:,:,curr_depth),1),'b',0.5,false);
        AP_errorfill(t,nanmean(curr_act_ctxpred_reduced_trialmean(:,:,curr_depth),1), ...
            AP_sem(curr_act_ctxpred_reduced_trialmean(:,:,curr_depth),1),[0,0.7,0],0.5,false);
        
        xlabel(['Time from ' align_labels{curr_align}]);
        ylabel(['Str ' num2str(curr_depth)])
        line([0,0],ylim,'color','k');
        if curr_align == 1 && curr_depth == 1
            legend({'Measured','Measured (task-reduced)','Task-pred','Ctx-pred'});
        end
    end
    
end

linkaxes(get(gcf,'Children'),'y');


%% Fig 4c: Plot kernel ROIs

% Plot kernel ROIs
figure;
for curr_roi = 1:n_depths
    subplot(n_depths,1,curr_roi,'YDir','reverse');
    
    curr_roi_boundary = bwboundaries(kernel_roi.bw(:,:,curr_roi));
    for curr_boundary = 1:length(curr_roi_boundary)
        patch(curr_roi_boundary{curr_boundary}(:,2), ...
            curr_roi_boundary{curr_boundary}(:,1),[0,0.8,0]);
    end
    
    AP_reference_outline('ccf_aligned','k');
    axis image off;
    title(['ROI for Str ' num2str(curr_roi)])
end


%% Fig 4d,e,f: Striatum v Cortex by condition

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_stim_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set windows to average activity
timeavg_labels = {'Stim','Move onset','Outcome'};
timeavg_t = {[0.05,0.15],[-0.05,0.05],[0.05,0.15]};
timeavg_align = {stim_align,move_align,outcome_align};
timeavg_trial_conditions = ...
    {[sign(trial_stim_allcat) == 1,sign(trial_stim_allcat) == -1], ...
    [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
    [trial_outcome_allcat == 1, trial_outcome_allcat == -1]};
timeavg_task_reduction = [1,2,4];

% Set activity percentiles and bins
act_prctile = [10,90];
n_act_bins = 5;

% Set areas and conditions
plot_str_ctx = {[1,1],[2,2],[3,3],[4,4]};

% Loop across area pairs, plot binned predicted v measured activity
str_v_ctx_fig = figure('color','w');
for curr_str_ctx = 1:length(plot_str_ctx)
    
    plot_str = plot_str_ctx{curr_str_ctx}(1);
    plot_ctx = plot_str_ctx{curr_str_ctx}(2);

    for curr_mua = 1:2
        
        % (set striatum activity to use)
        switch curr_mua
            case 1
                curr_str_act_allcat = single(mua_allcat);
                curr_str_act_taskpred_reduced_allcat = single(mua_taskpred_reduced_allcat);
                mua_label = 'Measured';
            case 2
                curr_str_act_allcat = single(mua_ctxpred_allcat);
                curr_str_act_taskpred_reduced_allcat = single(mua_ctxpred_taskpred_reduced_allcat);
                mua_label = 'Ctx-pred';
        end
        
        for curr_timeavg = 1:length(timeavg_labels)
            
            trial_conditions = timeavg_trial_conditions{curr_timeavg};
            curr_task_reduction = timeavg_task_reduction(curr_timeavg);
            
            % (set cortex activity to use)
%             curr_ctx_act_allcat = fluor_roi_deconv;
            curr_ctx_act_allcat = fluor_kernelroi_deconv;   
            curr_ctx_act_taskpred_reduced_allcat = fluor_kernelroi_taskpred_reduced;
            
            % (re-align and split activity and conditions)
            curr_ctx_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_ctx_act_allcat(trial,:,:) - curr_ctx_act_taskpred_reduced_allcat(trial,:,:,curr_task_reduction), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_ctx_act_allcat,1)),'uni',false)), ...
                use_split,length(t),size(curr_ctx_act_allcat,3));
            
            curr_str_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_str_act_allcat(trial,:,:) - curr_str_act_taskpred_reduced_allcat(trial,:,:,curr_task_reduction), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_str_act_allcat,1)),'uni',false)), ...
                use_split,length(t),size(curr_str_act_allcat,3));
            
            trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
            
            % (get average activity within window)
            curr_event_t = t > timeavg_t{curr_timeavg}(1) & t < timeavg_t{curr_timeavg}(2);
            curr_ctx_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_ctx_act,'uni',false);
            curr_str_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_str_act,'uni',false);
            
            % (bin measured data across percentile range)
            bin_range = prctile(cell2mat(cellfun(@(x) x(:,plot_ctx),curr_ctx_act_avg,'uni',false)),act_prctile);
            bin_edges = linspace(bin_range(1),bin_range(2),n_act_bins+1);
            bin_centers = bin_edges(1:end-1) + diff(bin_edges)./2;
            
            trial_bins = cellfun(@(x) discretize(x,bin_edges),curr_ctx_act_avg,'uni',false);
            total_bins = cellfun(@(x) discretize(x,[bin_edges(1),bin_edges(end)]),curr_ctx_act_avg,'uni',false);
            
            % (get trials to use: no NaNs in either time series or bins)
            nonan_trials = cellfun(@(ctx_act,str_act,trial_bins) ...
                squeeze(~any(isnan(ctx_act(:,curr_event_t,plot_ctx)),2)) & ...
                squeeze(~any(isnan(str_act(:,curr_event_t,plot_str)),2)) & ...
                ~isnan(trial_bins), ...
                curr_ctx_act,curr_str_act,trial_bins,'uni',false);
            
            % (get average binned activity for measured/predicted by condition)
            ctx_act_binmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_ctx_act_avg,trial_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            str_act_binmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_str_act_avg,trial_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            % (get the average total activity)
            ctx_act_totalmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [1,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_ctx_act_avg,total_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            str_act_totalmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [1,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_str_act_avg,total_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            % Plot binned predicted v measured           
            subplot(length(plot_str_ctx),length(timeavg_labels), ...
                sub2ind(fliplr([length(plot_str_ctx),length(timeavg_labels)]),curr_timeavg,curr_str_ctx));
            hold on;
            col = lines(size(trial_conditions,2));
            switch curr_mua
                case 1
                    errorbar( ...
                        squeeze(nanmean(ctx_act_binmean(:,plot_ctx,:,:),3)), ...
                        squeeze(nanmean(str_act_binmean(:,plot_str,:,:),3)), ...
                        squeeze(AP_sem(str_act_binmean(:,plot_str,:,:),3)), ...
                        'linewidth',2,'CapSize',0);
%                     errorbar( ...
%                         squeeze(nanmean(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
%                         squeeze(nanmean(str_act_totalmean(:,plot_str,:,:),3)), ...
%                         squeeze(AP_sem(str_act_totalmean(:,plot_str,:,:),3)), ...
%                         squeeze(AP_sem(str_act_totalmean(:,plot_str,:,:),3)), ...
%                         squeeze(AP_sem(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
%                         squeeze(AP_sem(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
%                         '.','linewidth',3,'CapSize',0);
                case 2
                    for curr_cond = 1:size(trial_conditions,2)
                        AP_errorfill( ...
                            squeeze(nanmean(ctx_act_binmean(:,plot_ctx,:,curr_cond),3)), ...
                            squeeze(nanmean(str_act_binmean(:,plot_str,:,curr_cond),3)), ...
                            squeeze(AP_sem(str_act_binmean(:,plot_str,:,curr_cond),3)), ...
                            col(curr_cond,:),0.5,false);
                    end
            end
            xlabel(['Ctx (' num2str(plot_ctx) ')']);
            ylabel([' Str (' num2str(plot_str) ')'])
            title([timeavg_labels{curr_timeavg} '(' task_regressor_labels{curr_task_reduction} '-reduced)']);
            
        end
        
    end
end

% Link axes of all plots
linkaxes(get(str_v_ctx_fig,'Children'));


%% Fig 4g: Cortex-predicted striatum error

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_stim_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set windows to average activity
timeavg_labels = {'Stim','Move onset','Outcome'};
timeavg_t = {[0.05,0.15],[-0.05,0.05],[0.05,0.15]};
timeavg_align = {stim_align,move_align,outcome_align};
timeavg_trial_conditions = ...
    {[sign(trial_stim_allcat) == 1,sign(trial_stim_allcat) == -1], ...
    [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
    [trial_outcome_allcat == 1, trial_outcome_allcat == -1]};
timeavg_task_reduction = [1,2,4];

% timeavg_labels = {'Pre-stim','Stim'};
% timeavg_t = {[-0.2,-0.1],[0.05,0.15]};
% timeavg_align = {stim_align,stim_align};

% Set activity percentiles and bins
act_prctile = [10,90];
n_act_bins = 5;

% Set areas and conditions
plot_areas = [1,2,3,4];

% Loop across area pairs, plot binned predicted v measured activity
curr_act_allcat = mua_allcat;
curr_act_taskpred_reduced_allcat = mua_taskpred_reduced_allcat;
curr_act_pred_allcat = mua_ctxpred_allcat;
curr_act_pred_taskpred_reduced_allcat = mua_ctxpred_taskpred_reduced_allcat;

measured_v_pred_fig = figure('color','w');
for curr_area_idx = 1:length(plot_areas)
    
    plot_area = plot_areas(curr_area_idx);
    
    % Set up the summary values
    curr_act_pred_diff = nan(2,length(use_split),length(timeavg_labels));
    curr_act_pred_rank_condition_diff = nan(length(use_split),length(timeavg_labels));
    n_shuff = 1000;
    curr_act_pred_rank_condition_diff_shuff = nan(length(use_split),n_shuff,length(timeavg_labels));
    
    for curr_timeavg = 1:length(timeavg_labels)
        
        trial_conditions = timeavg_trial_conditions{curr_timeavg};
        curr_task_reduction = timeavg_task_reduction(curr_timeavg);
        
        % (re-align and split activity)
        curr_act = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_allcat(trial,:,plot_area) - curr_act_taskpred_reduced_allcat(trial,:,plot_area,curr_task_reduction), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        curr_act_pred = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_pred_allcat(trial,:,plot_area) - curr_act_pred_taskpred_reduced_allcat(trial,:,plot_area,curr_task_reduction), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
        
        % (get average activity within window)
        curr_event_t = t > timeavg_t{curr_timeavg}(1) & t < timeavg_t{curr_timeavg}(2);
        curr_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act,'uni',false);
        curr_act_pred_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_pred,'uni',false);

        % (bin predicted data across percentile range)
        bin_range = prctile(cell2mat(curr_act_pred_avg),act_prctile);
        bin_edges = linspace(bin_range(1),bin_range(2),n_act_bins+1);
        bin_centers = bin_edges(1:end-1) + diff(bin_edges)./2;
        
        trial_bins = cellfun(@(x) discretize(x,bin_edges),curr_act_pred_avg,'uni',false);
        total_bins = cellfun(@(x) discretize(x,[bin_edges(1),bin_edges(end)]),curr_act_pred_avg,'uni',false);
        
        % (get trials to use: no NaNs in either time series or bins)
        use_trials = cellfun(@(act,act_taskpred,trial_bins) ...
            squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
            squeeze(~any(isnan(act_taskpred(:,curr_event_t)),2)) & ...
            ~isnan(trial_bins), ...
            curr_act,curr_act_pred,trial_bins,'uni',false);
        
        % (get average binned activity)
        act_binmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,trial_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_binmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,trial_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));      

        % (get the average total activity)       
        act_totalmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [1,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,total_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_totalmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [1,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,total_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        % Get average act-pred difference (ALL TRIALS)
        curr_act_pred_diff(:,:,curr_timeavg) = ...
            cell2mat(arrayfun(@(cond) cellfun(@(act,pred,trial_cond,use_trials) ...
            nanmean(act(trial_cond(:,cond)) - pred(trial_cond(:,cond))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp,use_trials), ...
            1:size(trial_conditions,2),'uni',false))';
        
        % Get condition  difference (assume 2 conditions)
        curr_act_avg_rank = cellfun(@(x) tiedrank(x)./max(tiedrank(x)),curr_act_avg,'uni',false);
        curr_act_pred_avg_rank = cellfun(@(x) tiedrank(x)./max(tiedrank(x)),curr_act_pred_avg,'uni',false);        
        curr_act_pred_condition_diff(:,curr_timeavg) = ...
            cellfun(@(act,pred,trial_cond,use_trials) ...       
            nanmean(act(trial_cond(:,1)) - pred(trial_cond(:,1))) - ...
            nanmean(act(trial_cond(:,2)) - pred(trial_cond(:,2))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp,use_trials);
        
        % Get condition-shuffled act-pred difference
        trial_conditions_exp_shuff = cellfun(@(x) AP_shake(repmat(x,1,1,n_shuff),1), ...
            trial_conditions_exp,'uni',false);   
                
        curr_act_pred_condition_diff_shuff(:,:,curr_timeavg) = ...
            cell2mat(arrayfun(@(shuff) ...
            cellfun(@(act,pred,trial_cond,use_trials) ...
            nanmean(act(trial_cond(:,1,shuff)) - pred(trial_cond(:,1,shuff))) - ...
            nanmean(act(trial_cond(:,2,shuff)) - pred(trial_cond(:,2,shuff))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp_shuff,use_trials), ...
            1:n_shuff,'uni',false));
        
        % Plot binned predicted v measured
        figure(measured_v_pred_fig);
        subplot(length(plot_areas),length(timeavg_labels)+2, ...
            sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),curr_timeavg,curr_area_idx));
        hold on;
        
        errorbar( ...
            squeeze(nanmean(act_pred_binmean,2)), ...
            squeeze(nanmean(act_binmean,2)), ...
            squeeze(AP_sem(act_binmean,2)), ...
            'linewidth',2,'CapSize',0);
        errorbar( ...
            squeeze(nanmean(act_pred_totalmean,2)), ...
            squeeze(nanmean(act_totalmean,2)), ...
            squeeze(AP_sem(act_totalmean,2)), ...
            squeeze(AP_sem(act_totalmean,2)), ...
            squeeze(AP_sem(act_pred_totalmean,2)), ...
            squeeze(AP_sem(act_pred_totalmean,2)), ...
            '.','linewidth',3,'CapSize',0);
        xlabel(['Predicted (' num2str(plot_area) ')']);
        ylabel(['Measured (' num2str(plot_area) ')'])
        ylim(xlim); axis square;
        title([timeavg_labels{curr_timeavg} ' (' task_regressor_labels{curr_task_reduction} '-reduced)']);
        
    end
    
    % Plot measured - predicted
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+1,curr_area_idx));
     hold on;
    errorbar( ...
        permute(nanmean(curr_act_pred_diff,2),[3,1,2]), ...
        permute(AP_sem(curr_act_pred_diff,2),[3,1,2]),'linewidth',2,'CapSize',0);
    ylabel('Meas - Pred');
    set(gca,'XTick',1:4,'XTickLabels',timeavg_labels,'XTickLabelRotation',45)
    xlim([0.5,length(timeavg_labels)+0.5]);
    
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+2,curr_area_idx));
     hold on;
    errorbar( ...
        nanmean(curr_act_pred_condition_diff,1), ...
        AP_sem(curr_act_pred_condition_diff,1),'linewidth',2,'CapSize',0);
    ylabel('Condition difference');
    set(gca,'XTick',1:4,'XTickLabels',timeavg_labels,'XTickLabelRotation',45)
    xlim([0.5,length(timeavg_labels)+0.5]);
    
    % Get and plot significance
    real_diff = nanmean(curr_act_pred_condition_diff,1); 
    alpha = [5/2/3,100-(5/2/3)];
    shuff_diff_ci = permute(nanmean(prctile(curr_act_pred_condition_diff_shuff,alpha,2),1),[2,3,1]);
    sig_diff = real_diff < shuff_diff_ci(1,:) | real_diff > shuff_diff_ci(2,:);
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+1,curr_area_idx));
    plot(find(sig_diff),repmat(max(ylim),sum(sig_diff),1),'*','color','k','MarkerSize',5);
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+2,curr_area_idx));
    plot(shuff_diff_ci','r');
    plot(find(sig_diff),repmat(max(ylim),sum(sig_diff),1),'*','color','k','MarkerSize',5);
    
end

% Link axes of all predicted v measured and all explained var
all_axes = get(measured_v_pred_fig,'Children');
meas_pred_axes = all_axes(setdiff(1:length(all_axes), ...
    [1:length(timeavg_labels)+2:length(all_axes), ...
    2:length(timeavg_labels)+2:length(all_axes)]));
linkaxes(meas_pred_axes)
linkaxes(all_axes(1:length(timeavg_labels)+2:length(all_axes)));
linkaxes(all_axes(2:length(timeavg_labels)+2:length(all_axes)));

for curr_axes = 1:length(meas_pred_axes)
   line(meas_pred_axes(curr_axes),xlim(meas_pred_axes(curr_axes)), ...
       xlim(meas_pred_axes(curr_axes)),'color','k'); 
end
for curr_axes = 1:length(timeavg_labels)+2:length(all_axes)
   line(all_axes(curr_axes),xlim(all_axes(curr_axes)),[0,0],'color','k');
end
for curr_axes = 2:length(timeavg_labels)+2:length(all_axes)
   line(all_axes(curr_axes),xlim(all_axes(curr_axes)),[0,0],'color','k');
end


%% Fig 4h: Cortex-predicted striatum error (passive choiceworld stim: trained/naive)
% run separately for trained/naive then combine posthoc with copyobj

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(mua_allcat,1),1);

% Set windows to average activity
timeavg_labels = {'Stim'};
timeavg_t = {[0.05,0.15]};
timeavg_align = {stim_align};
timeavg_trial_conditions = ...
    {[sign(trial_stim_allcat) == 1,sign(trial_stim_allcat) == -1]};

% Set activity percentiles and bins
act_prctile = [10,90];
n_act_bins = 5;

% Set areas and conditions
plot_areas = [1,2,3,4];

% Loop across area pairs, plot binned predicted v measured activity
curr_act_allcat = mua_allcat;
curr_act_pred_allcat = mua_ctxpred_allcat;

measured_v_pred_fig = figure('Name',data_fn,'color','w');
for curr_area_idx = 1:length(plot_areas)
    
    plot_area = plot_areas(curr_area_idx);
    
    % Set up the summary values
    curr_act_pred_diff = nan(size(timeavg_trial_conditions{1},2),length(use_split),length(timeavg_labels));
    curr_act_pred_rank_condition_diff = nan(length(use_split),length(timeavg_labels));
    n_shuff = 1000;
    curr_act_pred_rank_condition_diff_shuff = nan(length(use_split),n_shuff,length(timeavg_labels));
    
    for curr_timeavg = 1:length(timeavg_labels)
        
        trial_conditions = timeavg_trial_conditions{curr_timeavg};
        
        % (re-align and split activity)
        curr_act = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        curr_act_pred = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_pred_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
        
        % (get average activity within window)
        curr_event_t = t > timeavg_t{curr_timeavg}(1) & t < timeavg_t{curr_timeavg}(2);
        curr_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act,'uni',false);
        curr_act_pred_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_pred,'uni',false);

        % (bin predicted data across percentile range)
        bin_range = prctile(cell2mat(curr_act_pred_avg),act_prctile);
        bin_edges = linspace(bin_range(1),bin_range(2),n_act_bins+1);
        bin_centers = bin_edges(1:end-1) + diff(bin_edges)./2;
        
        trial_bins = cellfun(@(x) discretize(x,bin_edges),curr_act_pred_avg,'uni',false);
        total_bins = cellfun(@(x) discretize(x,[bin_edges(1),bin_edges(end)]),curr_act_pred_avg,'uni',false);
        
        % (get trials to use: no NaNs in either time series or bins)
        use_trials = cellfun(@(act,act_taskpred,trial_bins) ...
            squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
            squeeze(~any(isnan(act_taskpred(:,curr_event_t)),2)) & ...
            ~isnan(trial_bins), ...
            curr_act,curr_act_pred,trial_bins,'uni',false);
        
        % (get average binned activity)
        act_binmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,trial_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_binmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,trial_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));      

        % (get the average total activity)       
        act_totalmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [1,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,total_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_totalmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [1,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,total_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        % Get average act-pred difference (ALL TRIALS)
        curr_act_pred_diff(:,:,curr_timeavg) = ...
            cell2mat(arrayfun(@(cond) cellfun(@(act,pred,trial_cond,use_trials) ...
            nanmean(act(trial_cond(:,cond)) - pred(trial_cond(:,cond))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp,use_trials), ...
            1:size(trial_conditions,2),'uni',false))';
        
        % Get condition  difference (assume 2 conditions)
        curr_act_avg_rank = cellfun(@(x) tiedrank(x)./max(tiedrank(x)),curr_act_avg,'uni',false);
        curr_act_pred_avg_rank = cellfun(@(x) tiedrank(x)./max(tiedrank(x)),curr_act_pred_avg,'uni',false);        
        curr_act_pred_condition_diff(:,curr_timeavg) = ...
            cellfun(@(act,pred,trial_cond,use_trials) ...       
            nanmean(act(trial_cond(:,1)) - pred(trial_cond(:,1))) - ...
            nanmean(act(trial_cond(:,2)) - pred(trial_cond(:,2))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp,use_trials);
        
        % Get condition-shuffled act-pred difference
        trial_conditions_exp_shuff = cellfun(@(x) AP_shake(repmat(x,1,1,n_shuff),1), ...
            trial_conditions_exp,'uni',false);   
                
        curr_act_pred_condition_diff_shuff(:,:,curr_timeavg) = ...
            cell2mat(arrayfun(@(shuff) ...
            cellfun(@(act,pred,trial_cond,use_trials) ...
            nanmean(act(trial_cond(:,1,shuff)) - pred(trial_cond(:,1,shuff))) - ...
            nanmean(act(trial_cond(:,2,shuff)) - pred(trial_cond(:,2,shuff))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp_shuff,use_trials), ...
            1:n_shuff,'uni',false));
        
        % Plot binned predicted v measured
        figure(measured_v_pred_fig);
        subplot(length(plot_areas),length(timeavg_labels)+2, ...
            sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),curr_timeavg,curr_area_idx));
        hold on;
        
        errorbar( ...
            squeeze(nanmean(act_pred_binmean,2)), ...
            squeeze(nanmean(act_binmean,2)), ...
            squeeze(AP_sem(act_binmean,2)), ...
            'linewidth',2,'CapSize',0);
        errorbar( ...
            squeeze(nanmean(act_pred_totalmean,2)), ...
            squeeze(nanmean(act_totalmean,2)), ...
            squeeze(AP_sem(act_totalmean,2)), ...
            squeeze(AP_sem(act_totalmean,2)), ...
            squeeze(AP_sem(act_pred_totalmean,2)), ...
            squeeze(AP_sem(act_pred_totalmean,2)), ...
            '.','linewidth',3,'CapSize',0);
        xlabel(['Predicted (' num2str(plot_area) ')']);
        ylabel(['Measured (' num2str(plot_area) ')'])
        ylim(xlim); axis square;
        title([timeavg_labels{curr_timeavg}]);
        
    end
    
    % Plot measured - predicted
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+1,curr_area_idx));
     hold on;
    errorbar( ...
        permute(nanmean(curr_act_pred_diff,2),[3,1,2]), ...
        permute(AP_sem(curr_act_pred_diff,2),[3,1,2]),'linewidth',2,'CapSize',0);
    ylabel('Meas - Pred');
    xlabel('Stim');
    set(gca,'XTick',1:2,'XTickLabel',{'Contra','Ipsi'});
    xlim([0.5,size(timeavg_trial_conditions{1},2)+0.5]);
    
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+2,curr_area_idx));
     hold on;
    errorbar( ...
        nanmean(curr_act_pred_condition_diff,1), ...
        AP_sem(curr_act_pred_condition_diff,1),'linewidth',2,'CapSize',0);
    ylabel('Condition difference');
    xlabel('Stim');
    set(gca,'XTick',1:2,'XTickLabel',{'Contra','Ipsi'});
    xlim([0.5,size(timeavg_trial_conditions{1},2)+0.5]);
    
    % Get and plot significance
    real_diff = nanmean(curr_act_pred_condition_diff,1); 
    alpha = [5/2/3,100-(5/2/3)];
    shuff_diff_ci = permute(nanmean(prctile(curr_act_pred_condition_diff_shuff,alpha,2),1),[2,3,1]);
    sig_diff = real_diff < shuff_diff_ci(1,:) | real_diff > shuff_diff_ci(2,:);
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+1,curr_area_idx));
    plot(find(sig_diff),repmat(max(ylim),sum(sig_diff),1),'*','color','k','MarkerSize',5);
    
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+2,curr_area_idx));
    plot(shuff_diff_ci','r');
    plot(find(sig_diff),repmat(max(ylim),sum(sig_diff),1),'*','color','k','MarkerSize',5);
    
end

% Link axes of all predicted v measured and all explained var
all_axes = get(measured_v_pred_fig,'Children');
meas_pred_axes = all_axes(setdiff(1:length(all_axes), ...
    [1:length(timeavg_labels)+2:length(all_axes), ...
    2:length(timeavg_labels)+2:length(all_axes)]));
linkaxes(meas_pred_axes)
linkaxes(all_axes(1:length(timeavg_labels)+2:length(all_axes)));
linkaxes(all_axes(2:length(timeavg_labels)+2:length(all_axes)));

for curr_axes = 1:length(meas_pred_axes)
   line(meas_pred_axes(curr_axes),xlim(meas_pred_axes(curr_axes)), ...
       xlim(meas_pred_axes(curr_axes)),'color','k'); 
end
for curr_axes = 1:length(timeavg_labels)+2:length(all_axes)
   line(all_axes(curr_axes),xlim(all_axes(curr_axes)),[0,0],'color','k');
end
for curr_axes = 2:length(timeavg_labels)+2:length(all_axes)
   line(all_axes(curr_axes),xlim(all_axes(curr_axes)),[0,0],'color','k');
end


%% Striatum 2 stim/move activity interaction

% Depth to plot
plot_depth = 2;

% Set number of shuffles for significance testing
n_shuff = 1000; 

%%% STIM ACTIVITY

% Set time to average activity
use_stim_t = t > 0 & t < 0.2;

% Get stim bins for each trial
stims = unique([0.06,0.125,0.25,0.5,1].*[-1;1]);
stim_bins_idx = discretize(trial_stim_allcat,stims);

% Set contexts for activity
stim_contexts = [trial_choice_allcat == -1, ...
    trial_choice_allcat == 1];

% Get stim-isolated activity by trial
stim_activity = mua_allcat - mua_taskpred_reduced_allcat(:,:,:,1);
stim_trial_activity = nanmean(stim_activity(:,use_stim_t,:),2);

% Split activity by animal, stim, and correct/incorrect
stim_bins = stims;
stim_trial_activity_split = cell(max(split_idx),length(stim_bins),size(stim_contexts,2),n_depths);
for curr_exp = 1:max(split_idx)
    for curr_bin = 1:length(stim_bins)
        for curr_context = 1:size(stim_contexts,2)               
            for curr_depth = 1:n_depths
                % Get activity of all group trials (exclude NaNs)
                curr_trials = ...
                    split_idx == curr_exp & ...
                    trial_stim_allcat == stim_bins(curr_bin) & ...
                    stim_contexts(:,curr_context);             
                curr_activity = stim_trial_activity(curr_trials,:,curr_depth);
                stim_trial_activity_split{curr_exp,curr_bin,curr_context,curr_depth} = ...
                    curr_activity(~isnan(curr_activity));                
            end
            
        end
    end
end

stim_trial_activity_split_mean = cellfun(@nanmean,stim_trial_activity_split);

% Plot stim activity by movement direction
figure('Name',['Str ' num2str(plot_depth)]);
subplot(1,2,1); hold on;
line_col = [0.6,0,0.6;0,0.6,0];

dot_col = colormap_BlueWhiteRed(5);
dot_col(6,:) = [];

p = nan(size(stim_contexts,2),1);
for curr_context = 1:size(stim_contexts,2)
    p(curr_context) = ...
        errorbar(stims,nanmean(stim_trial_activity_split_mean(:,:,curr_context,plot_depth),1), ...
        AP_sem(stim_trial_activity_split_mean(:,:,curr_context,plot_depth),1),'color',line_col(curr_context,:),'linewidth',3);
    scatter(stims,nanmean(stim_trial_activity_split_mean(:,:,curr_context,plot_depth),1),80,dot_col, ...
        'Filled','MarkerEdgeColor',line_col(curr_context,:),'linewidth',3);
end
xlabel('Stimulus');
ylabel('Activity');
legend(p,{'Move left','Move right'});

% Get condition difference and compare to shuffle
subplot(1,2,2); hold on;

stim_condition_diff = nanmean(stim_trial_activity_split_mean(:,:,1,plot_depth) - ...
    stim_trial_activity_split_mean(:,:,2,plot_depth),1);

stim_condition_shuff_diff = nan(max(split_idx),length(stim_bins),n_shuff);
for curr_exp = 1:max(split_idx)
    for curr_bin = 1:length(stim_bins)
        curr_act = stim_trial_activity_split(curr_exp,curr_bin,:,plot_depth);
        % Shuffle, split, difference
        curr_act_shuff = mat2cell(AP_shake(repmat(vertcat( ...
            curr_act{:}),1,n_shuff),1),cellfun(@length,curr_act),n_shuff);
        stim_condition_shuff_diff(curr_exp,curr_bin,:) = ...
            nanmean(curr_act_shuff{1},1) - nanmean(curr_act_shuff{2},1);            
    end
end
stim_condition_diff_ci = squeeze(prctile(nanmean(stim_condition_shuff_diff,1),[2.5,97.5],3));

plot(stims,stim_condition_diff,'k','linewidth',2);
plot(stims,stim_condition_diff_ci,'k','linewidth',2,'linestyle','--');
xlabel('Stimulus');
ylabel('Condition difference');


%%% MOVEMENT ACTIVITY

% Set time to average activity
use_move_t = t > 0 & t < 0.2;

% Get velocity bins for each trial
max_vel = AP_signed_max(wheel_velocity_allcat_move(:,t > 0 & t < 0.2),2);

% Normalize velocity for each day
max_vel_exp = mat2cell(max_vel,trials_recording,1);
max_vel_norm = cell2mat(cellfun(@(x) x./nanmean(abs(x)),max_vel_exp,'uni',false));

n_vel_bins = 4;
vel_edges = prctile(abs(max_vel_norm),linspace(0,100,n_vel_bins+1));
vel_edges = linspace(min(abs(max_vel_norm)),max(abs(max_vel_norm)),n_vel_bins+1);

vel_edges = sort([vel_edges,-vel_edges]);
vel_centers = vel_edges(1:end-1) + diff(vel_edges)/2;

vel_bins = discretize(max_vel_norm,vel_edges);

% Set contexts for activity
move_contexts = [trial_stim_allcat > 0, ...
    trial_stim_allcat < 0, ...
    trial_stim_allcat == 0];

% Get move-isolated activity by trial
move_activity = mua_allcat_move - mua_taskpred_reduced_allcat_move(:,:,:,2);
move_trial_activity = nanmean(move_activity(:,use_move_t,:),2);

% Split activity by animal, velocity, and stim side/zero contrast
move_trial_activity_split = cell(max(split_idx),max(vel_bins),size(move_contexts,2),n_depths);
for curr_exp = 1:max(split_idx)
    for curr_bin = 1:max(vel_bins)
        for curr_context = 1:size(move_contexts,2)               
            for curr_depth = 1:n_depths
                % Get activity of all group trials (exclude NaNs)
                curr_trials = ...
                    split_idx == curr_exp & ...
                    vel_bins == curr_bin & ...
                    move_contexts(:,curr_context);              
                curr_activity = move_trial_activity(curr_trials,:,curr_depth);
                move_trial_activity_split{curr_exp,curr_bin,curr_context,curr_depth} = ...
                    curr_activity(~isnan(curr_activity));                
            end
            
        end
    end
end

move_trial_activity_split_mean = cellfun(@nanmean,move_trial_activity_split);

% Plot move activity by stim side
figure('Name',['Str ' num2str(plot_depth)]);
subplot(1,2,1); hold on;
line_col = [0.7,0,0;0,0,0.7;0.5,0.5,0.5];

dot_col = [linspace(0.2,1,n_vel_bins)',0.1*ones(n_vel_bins,1),linspace(0.2,1,n_vel_bins)'; ...
    1,0,0;
    0.1*ones(n_vel_bins,1),linspace(1,0.2,n_vel_bins)',0.1*ones(n_vel_bins,1)];

p = nan(size(move_contexts,2),1);
for curr_context = 1:size(move_contexts,2)
    p(curr_context) = ...
        errorbar(vel_centers,nanmean(move_trial_activity_split_mean(:,:,curr_context,plot_depth),1), ...
        AP_sem(move_trial_activity_split_mean(:,:,curr_context,plot_depth),1),'color',line_col(curr_context,:),'linewidth',3);
    scatter(vel_centers,nanmean(move_trial_activity_split_mean(:,:,curr_context,plot_depth),1),80,dot_col, ...
        'Filled','MarkerEdgeColor',line_col(curr_context,:),'linewidth',3);
end
xlabel('Velocity');
ylabel('Activity');
legend(p,{'Stim right','Stim left','No stim'});
 
% Get condition difference and compare to shuffle
subplot(1,2,2); hold on;

move_condition_diff_1 = nanmean(move_trial_activity_split_mean(:,:,1,plot_depth) - ...
    move_trial_activity_split_mean(:,:,2,plot_depth),1);
move_condition_diff_2 = nanmean(move_trial_activity_split_mean(:,:,2,plot_depth) - ...
    move_trial_activity_split_mean(:,:,3,plot_depth),1);
move_condition_diff_3 = nanmean(move_trial_activity_split_mean(:,:,1,plot_depth) - ...
    move_trial_activity_split_mean(:,:,3,plot_depth),1);

move_condition_shuff_diff_1 = nan(max(split_idx),max(vel_bins),n_shuff);
move_condition_shuff_diff_2 = nan(max(split_idx),max(vel_bins),n_shuff);
move_condition_shuff_diff_3 = nan(max(split_idx),max(vel_bins),n_shuff);
for curr_exp = 1:max(split_idx)
    for curr_bin = 1:max(vel_bins)
        curr_act = move_trial_activity_split(curr_exp,curr_bin,:,plot_depth);
        % Shuffle, split, difference
        curr_act_shuff = mat2cell(AP_shake(repmat(vertcat( ...
            curr_act{:}),1,n_shuff),1),cellfun(@length,curr_act),n_shuff);
        move_condition_shuff_diff_1(curr_exp,curr_bin,:) = ...
            nanmean(curr_act_shuff{1},1) - nanmean(curr_act_shuff{2},1);            
        move_condition_shuff_diff_2(curr_exp,curr_bin,:) = ...
            nanmean(curr_act_shuff{2},1) - nanmean(curr_act_shuff{3},1);            
        move_condition_shuff_diff_3(curr_exp,curr_bin,:) = ...
            nanmean(curr_act_shuff{1},1) - nanmean(curr_act_shuff{3},1);            
    end
end
move_condition_diff_1_ci = squeeze(prctile(nanmean(move_condition_shuff_diff_1,1),[2.5,97.5],3));
move_condition_diff_2_ci = squeeze(prctile(nanmean(move_condition_shuff_diff_2,1),[2.5,97.5],3));
move_condition_diff_3_ci = squeeze(prctile(nanmean(move_condition_shuff_diff_3,1),[2.5,97.5],3));

col = lines(3);
plot(vel_centers,move_condition_diff_1,'color',col(1,:),'linewidth',2);
plot(vel_centers,move_condition_diff_1_ci,'color',col(1,:),'linewidth',2,'linestyle','--');

plot(vel_centers,move_condition_diff_2,'color',col(2,:),'linewidth',2);
plot(vel_centers,move_condition_diff_2_ci,'color',col(2,:),'linewidth',2,'linestyle','--');

plot(vel_centers,move_condition_diff_3,'color',col(3,:),'linewidth',2);
plot(vel_centers,move_condition_diff_3_ci,'color',col(3,:),'linewidth',2,'linestyle','--');
xlabel('Velocity');
ylabel('Condition difference');























