%% Test analysis for recording across learning

%% Passive stim during choiceworld
% AP040/AP041 so far: plot average passive response over days
% (make sure to take out traces with wheel movement, but really should
% probably get movement based on camera with AP_mouse_movie_movement)


% get within stim window: 
% frame_movement, wheel velocity


animal = 'AP040';

protocol = 'AP_lcrGratingPassive';
experiments = AP_find_experiments(animal,protocol);

im_stim_all = cell(size(experiments));
for curr_day = 1:length(experiments)
    
    % Load data
    day = experiments(curr_day).day;
    experiment = experiments(curr_day).experiment(end);
    load_parts.imaging = true;
    AP_load_experiment;
    
    fVdf_deconv = AP_deconv_wf(fVdf);
    
    % Get wheel movements during stim
    wheel_window = [0,0.5];
    wheel_window_t = wheel_window(1):1/framerate:wheel_window(2);
    wheel_window_t_peri_event = bsxfun(@plus,stimOn_times,wheel_window_t);
    event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
        wheel_velocity,wheel_window_t_peri_event);
    wheel_thresh = 0.025;
    quiescent_trials = ~any(abs(event_aligned_wheel) > wheel_thresh,2);
    
    % Set options
    surround_window = [-0.5,1];
    baseline_window = [-0.1,0];
    
    surround_samplerate = 1/(framerate*1);
    surround_time = surround_window(1):surround_samplerate:surround_window(2);
    baseline_surround_time = baseline_window(1):surround_samplerate:baseline_window(2);
    
    % Average (time course) responses
    use_vs = 1:size(U,3);
    
    conditions = unique(stimIDs);
    im_stim = nan(size(U,1),size(U,2),length(surround_time),length(conditions));
    for curr_condition_idx = 1:length(conditions)
        curr_condition = conditions(curr_condition_idx);
        
        use_stims = find(stimIDs == curr_condition & quiescent_trials);
        use_stimOn_times = stimOn_times(use_stims);
        
        stim_surround_times = bsxfun(@plus, use_stimOn_times(:), surround_time);
        stim_baseline_surround_times = bsxfun(@plus, use_stimOn_times(:), baseline_surround_time);
        
        peri_stim_v = permute(interp1(frame_t,fVdf',stim_surround_times),[3,2,1]);
        baseline_v = permute(nanmean(interp1(frame_t,fVdf',stim_baseline_surround_times),2),[3,2,1]);
        
        stim_v_mean = nanmean(peri_stim_v - baseline_v,3);
        
        im_stim(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf(:,:,use_vs), ...
            stim_v_mean(use_vs,:));
    end
    
    im_stim_aligned = AP_align_widefield(im_stim,animal,day);
    
    im_stim_all{curr_day} = im_stim_aligned;

    AP_print_progress_fraction(curr_day,length(experiments));
end

im_stim_cat = cell2mat(permute(cellfun(@(im) ...
    reshape(permute(im,[1,2,4,3]),size(im,1),[],size(im,3)), ...
    im_stim_all,'uni',false),[2,3,4,1]));

AP_image_scroll(im_stim_cat);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));

t_stim = 20:22;
im_stim_avg = squeeze(max(im_stim_cat(:,:,t_stim,:),[],3));
AP_image_scroll(im_stim_avg);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));




%% Get task -> widefield kernel

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

% Prepare fluorescence

% NOT YET ANIMAL ALIGNED
% % Convert U to master U
% load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
% Udf_aligned = AP_align_widefield(Udf,animal,day);
% fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);

U_master = Udf;
Udf_aligned = Udf;
fVdf_recast = fVdf;

% Set components to keep
use_components = 1:200;

% Get event-aligned activity
raster_window = [-0.5,2];
upsample_factor = 1;
raster_sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):raster_sample_rate:raster_window(2);

% Get align times
use_align = stimOn_times;
use_align(isnan(use_align)) = 0;

t_peri_event = bsxfun(@plus,use_align,t);
t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];

%%% Trial-align cortex
event_aligned_V = ...
    interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);

% Get time points to bin
sample_rate = framerate*regression_params.upsample_factor;
time_bins = frame_t(find(frame_t > ...
    regression_params.skip_seconds,1)):1/sample_rate: ...
    frame_t(find(frame_t-frame_t(end) < ...
    -regression_params.skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Deconvolve fluoresence
fVdf_deconv = AP_deconv_wf(fVdf);
fVdf_deconv(isnan(fVdf_deconv)) = 0;
fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';

%%% Trial-align wheel velocity
event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
    wheel_velocity,t_peri_event);

%%% Regress task to cortex/striatum/cortex-predicted striatum

% Get reaction time for building regressors
% (wheel thresh is the max in the quiescent period)
wheel_thresh = max(reshape(abs(event_aligned_wheel(:,t < -0.2)),[],1));
[move_trial,move_idx] = max(abs(event_aligned_wheel) > wheel_thresh,[],2);
move_idx(~move_trial) = NaN;
move_t = nan(size(move_idx));
move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';

% Build regressors (only a subset of these are used)

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

move_onset_stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
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
% (for go cue only on late move trials)
%         go_cue_regressors = histcounts( ...
%             signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
% (for go cue with early/late move trials)
go_cue_regressors = zeros(1,length(time_bin_centers));
go_cue_regressors(1,:) = histcounts( ...
    signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
go_cue_regressors(2,:) = histcounts( ...
    signals_events.interactiveOnTimes(move_t > 0.5),time_bins);

% Outcome regressors
% (using signals timing - not precise but looks good)
% (regressors for hit only)
%         outcome_regressors = histcounts(reward_t_timeline,time_bins);
% (regressors for both hit and miss)
outcome_regressors = zeros(2,length(time_bin_centers));
outcome_regressors(1,:) = histcounts( ...
    reward_t_timeline,time_bins);
outcome_regressors(2,:) = histcounts( ...
    signals_events.responseTimes(trial_outcome == -1),time_bins);

% Concatenate selected regressors, set parameters

task_regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
task_regressor_labels = {'Stim','Move onset','Go cue','Outcome'};

task_t_shifts = { ...
    [0,0.5]; ... % stim
    [-0.5,1]; ... % move
    [0,0.5]; ... % go cue
    [0,0.5]}; % outcome

task_regressor_sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
    round(x(2)*(sample_rate)),task_t_shifts,'uni',false);

lambda = 0;
zs = [false,false];
cvfold = 5;
use_constant = false;
return_constant = false;

% Regression task -> (master U, deconvolved) fluor
event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);

baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,:),[],size(event_aligned_V_deconv,3)))';
activity = single(fVdf_deconv_resample_recast(use_components,:))-baseline;

[fluor_taskpred_k,fluor_taskpred_long,fluor_taskpred_expl_var,fluor_taskpred_reduced_long] = ...
    AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
    lambda,zs,cvfold,return_constant,use_constant);


% Get regressor pixels
regressor_px = cellfun(@(v) cell2mat(arrayfun(@(subregressor) ...
    svdFrameReconstruct(U_master(:,:,use_components),permute(v(subregressor,:,:),[3,2,1])), ...
    permute(1:size(v,1),[1,3,4,2]),'uni',false)),fluor_taskpred_k,'uni',false);


AP_image_scroll(regressor_px{1});
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));













