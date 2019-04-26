% Generate figures for ctx-str paper

% Anything that takes a lot of time is done in
% AP_ctx_str_trial_preprocessing and saved for plotting here

% The original scripts here were in test_wf_ephys_choiceworld_analysis

% (this is addendum to AP_ctx_str_figures_v2, make comprehensive later)

%% Load in task data

% Load data
data_fn = 'trial_activity_choiceworld';
exclude_data = true;
AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(mua_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_animal;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));


%% Figure 1f: Cortical task-explained variance

use_t = true(size(t));
spatial_downsample = 10;

correct_trials = trial_side_allcat == -trial_choice_allcat;
incorrect_trials = trial_side_allcat == trial_choice_allcat;

visual_trials = trial_contrast_allcat > 0;
zero_trials = trial_contrast_allcat == 0;

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



%% Figure 1b: Example recordings

animal = 'AP025'; 
day = '2017-10-01'; 
experiment = 1; 
str_align = 'none'; 
verbose = true; 
AP_load_experiment;

% Load task->unit regression
unit_kernel_dir = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
unit_kernel_fn = 'unit_kernel_all_triaged.mat';
load([unit_kernel_dir filesep unit_kernel_fn]);
regressor_labels = {'Stim','Move','Go cue','Outcome'};
regressor_cols = [01,0,0;1,0,1;0.8,0.7,0.2;0,0,1];

curr_animal = 2;
curr_day = 3;
use_partial = 2;
[~,max_regressor_idx] = max(unit_kernel_all(curr_animal,curr_day).unit_expl_var_partial(:,:,use_partial),[],2);

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
    colormap(crameri('cork'));
    caxis([-0.03,0.03]);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title(sprintf('Deconv - %0.2fs',frame_t(plot_frames_idx(curr_frame))));
    axis image off;
    
    subplot(2,length(plot_frames_idx),length(plot_frames_idx) + curr_frame);
    imagesc(plot_frames_taskpred(:,:,curr_frame));
    colormap(crameri('cork'));
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






%% Figure 1g-h: Task -> striatum unit regression (example and summary histogram)

% Load the unit kernel results
unit_kernel_dir = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
unit_kernel_fn = 'unit_kernel_all_triaged.mat';
load([unit_kernel_dir filesep unit_kernel_fn]);
regressor_labels = {'Stim','Move','Go cue','Outcome'};
regressor_cols = [01,0,0;1,0,1;0.8,0.7,0.2;0,0,1];

% Get estimation of end of striatum for each recording
ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
ephys_align_fn = ['ephys_depth_align.mat'];
load([ephys_align_path filesep ephys_align_fn]);

% Plot example recording
curr_animal = 2;
curr_day = 3;

curr_str_depths = ephys_depth_align(curr_animal).str_depth(curr_day,:);

use_partial = 2;
[~,max_regressor_idx] = max(unit_kernel_all(curr_animal,curr_day).unit_expl_var_partial(:,:,use_partial),[],2);

figure;

subplot(1,length(regressor_labels)+2,1,'YDir','reverse'); hold on;
curr_expl_var = unit_kernel_all(curr_animal,curr_day).unit_expl_var_total;
scatter3(log10(unit_kernel_all(curr_animal,curr_day).spike_rate), ...
    unit_kernel_all(curr_animal,curr_day).template_depths, ...
    1:length(unit_kernel_all(curr_animal,curr_day).template_depths),20, ...
    regressor_cols(max_regressor_idx,:),'filled')
line(xlim,repmat(curr_str_depths(1),1,2),'color','k')
line(xlim,repmat(curr_str_depths(2),1,2),'color','k')
xlabel('Normalized n spikes');
ylabel('Depth (\mum)');
title({curr_animal,curr_day,'Best regressor'})

subplot(1,length(regressor_labels)+2,2,'YDir','reverse'); hold on;
curr_expl_var = unit_kernel_all(curr_animal,curr_day).unit_expl_var_total;
curr_expl_var_dotsize = 100*rescale(curr_expl_var.*(curr_expl_var > 0),0,1) + 1;
scatter3(log10(unit_kernel_all(curr_animal,curr_day).spike_rate), ...
    unit_kernel_all(curr_animal,curr_day).template_depths, ...
    1:length(unit_kernel_all(curr_animal,curr_day).template_depths),curr_expl_var_dotsize, ...
    regressor_cols(max_regressor_idx,:),'filled')
line(xlim,repmat(curr_str_depths(1),1,2),'color','k')
line(xlim,repmat(curr_str_depths(2),1,2),'color','k')
xlabel('Normalized n spikes');
ylabel('Depth (\mum)');
title('Total expl var');

% (plot units size-scaled by explained variance)
for curr_regressor = 1:length(regressor_labels)
    subplot(1,length(regressor_labels)+2,curr_regressor+2,'YDir','reverse');
    hold on;
    curr_expl_var = unit_kernel_all(curr_animal,curr_day).unit_expl_var_partial(:,curr_regressor,use_partial);
    curr_expl_var_dotsize = 100*rescale(curr_expl_var.*(curr_expl_var > 0),0,1) + 1;
    
    scatter3(log10(unit_kernel_all(curr_animal,curr_day).spike_rate), ...
        unit_kernel_all(curr_animal,curr_day).template_depths, ...
        1:length(unit_kernel_all(curr_animal,curr_day).template_depths), ...
        curr_expl_var_dotsize, ...
        regressor_cols(curr_regressor,:),'filled')
    line(xlim,repmat(curr_str_depths(1),1,2),'color','k')
    line(xlim,repmat(curr_str_depths(2),1,2),'color','k')
    title(regressor_labels{curr_regressor});
    xlabel('Normalized n spikes');
    ylabel('Depth (\mum)');
end



% Concatenate all units (normalized distance)

% Get estimation of end of striatum for each recording
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

use_partial = 2;

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

[~,max_regressor_idx_cat] = cellfun(@(x) max(x(:,:,use_partial),[],2),unit_expl_var_partial_cat,'uni',false);

h = figure;

subplot(1,length(regressor_labels)+3,1,'YDir','reverse'); hold on;
scatter(log10(cell2mat(spike_rate_cat)), ...
    cell2mat(template_depths_cat),20, ...
    regressor_cols(cell2mat(max_regressor_idx_cat),:),'filled')
xlabel('log10(spike rate)');
ylabel('Depth from str end (\mum)');
    title({'All units','Best regressor'})
line(xlim,[0,0],'color','k','linewidth',2);

subplot(1,length(regressor_labels)+3,2,'YDir','reverse'); hold on;
curr_expl_var = cell2mat(unit_expl_var_total_cat);
curr_expl_var_dotsize = 100*rescale(curr_expl_var.*(curr_expl_var > 0),0,1) + 1;
curr_expl_var_dotsize(isnan(curr_expl_var) | curr_expl_var < -1 | curr_expl_var > 1) = NaN;
scatter(log10(cell2mat(spike_rate_cat)), ...
    cell2mat(template_depths_cat),curr_expl_var_dotsize, ...
    regressor_cols(cell2mat(max_regressor_idx_cat),:),'filled')
xlabel('log10(spike rate)');
ylabel('Depth from str end (\mum)');
    title({'All regressors'})
line(xlim,[0,0],'color','k','linewidth',2);

% (plot units size-scaled by explained variance)
for curr_regressor = 1:length(regressor_labels)
    
    p1 = subplot(1,length(regressor_labels)+3,curr_regressor+2,'YDir','reverse');
    hold on;
    curr_expl_var = cellfun(@(x) x(:,curr_regressor,use_partial),unit_expl_var_partial_cat,'uni',false);
  
    % Normalize across experiments
    curr_expl_var_cat = cell2mat(curr_expl_var);
    curr_expl_var_norm = curr_expl_var_cat./ ...
        max(curr_expl_var_cat(curr_expl_var_cat < 1));
    curr_expl_var_dotsize = 100*rescale(curr_expl_var_norm.*(curr_expl_var_norm > 0),0,1) + 1;
    norm_type = 'concat normalized';    
    
    scatter(log10(cell2mat(spike_rate_cat)), ...
        cell2mat(template_depths_cat), ...
        curr_expl_var_dotsize, ...
        regressor_cols(curr_regressor,:),'filled');
    title({regressor_labels{curr_regressor},norm_type});
    xlabel('log10(spike rate)');
    ylabel('Depth from str end (\mum)');
    line(xlim,[0,0],'color','k','linewidth',2);
    
    p2 = subplot(1,length(regressor_labels)+3,length(regressor_labels)+3,'YDir','reverse');
    hold on;
    depth_bins = linspace(min(cell2mat(template_depths_cat)), ...
        max(cell2mat(template_depths_cat)),round(range(cell2mat(template_depths_cat))/400));
    depth_bin_centers = depth_bins(1:end-1) + diff(depth_bins)./2;
    curr_depth_groups = discretize(cell2mat(template_depths_cat),depth_bins);
    
    rate_cutoff = -0.5;
    use_units = log10(cell2mat(spike_rate_cat)) > rate_cutoff;
    line(p1,repmat(rate_cutoff,1,2),ylim(p1),'color','k');
    
    curr_expl_var_norm_depth = accumarray(curr_depth_groups(use_units), ...
        curr_expl_var_norm(use_units),[length(depth_bin_centers),1],@nanmean);
    plot(rescale(curr_expl_var_norm_depth,0,1),depth_bin_centers,'color',regressor_cols(curr_regressor,:),'linewidth',2);
    title({'Binned explained var',norm_type});
    xlabel('Normalized explained var');
    ylabel('Depth from str end (\mum)')
    line(xlim,[0,0],'color','k','linewidth',2);
    
end

% (draw striatum starts)
subplot(1,length(regressor_labels)+3,length(regressor_labels)+3)
str_start_line = nan(size(str_depth_cat,1),1);
for i = 1:size(str_depth_cat,1)
    str_start_line(i) = line(xlim,repmat(-diff(str_depth_cat(i,:)),1,2),'color',[0.8,0.8,0.8]);
end
set(gca,'Children',circshift(get(gca,'Children'),-size(str_depth_cat,1)));

linkaxes(get(h,'Children'),'y');


%% Figure 1: Striatal unit regression examples

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


%% ~~~~~~~~~~~~~ UNUSED ~~~~~~~~~~~



%% Compare kernels with domain maps

% Load data
data_fn = 'trial_activity_choiceworld';
exclude_data = true;
AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(mua_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_animal;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));



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
    AP_image_scroll(curr_k_px,t_shifts{curr_regressor});
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
        colormap(brewermap([],'Greens'));
        caxis([0,max_c]);
    end
end





% Get average cortex->striatum kernel
ctx_str_k_mean = nanmean(cell2mat(permute(vertcat(ctx_str_k_all{:}),[2,3,4,5,1])),5);
ctx_str_k_mean_px = cell2mat(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:50), ...
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









