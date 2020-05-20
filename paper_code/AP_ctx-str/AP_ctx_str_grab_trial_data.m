%% ~~~ AP_ctx_str_grab_trial_data: pull out data organized by trials ~~~

%% Set flags

% filter_mua: convolve MUA with filter to match widefield fluorescence
if ~exist('filter_mua','var')
    filter_mua = false;
end

%% Get if task if passive dataset

% Task dataset if signals (expDef) and contains 'vanillaChoiceworld'
task_dataset = exist('expDef','var') && contains(expDef,'vanillaChoiceworld');


%% Set parameters for cortical fluoresence

if verbose; disp('Re-casting and deconvolving fluorescence...'); end;

% Convert U to master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
Udf_aligned = AP_align_widefield(Udf,animal,day);
fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);

% Set components to keep
use_components = 1:200;

% Deconvolve fluoresence (for use in regression)
fVdf_deconv = AP_deconv_wf(fVdf);
fVdf_deconv(isnan(fVdf_deconv)) = 0;


%% Set parameters for striatal multiunit

n_depths = n_aligned_depths;
depth_group = aligned_str_depth_group;


%% Set parameters for trials

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

% Pick trials to keep
if task_dataset
    use_trials = ...
        trial_outcome ~= 0 & ...
        ~signals_events.repeatTrialValues(1:n_trials)' & ...
        stim_to_feedback < 1.5;    
else
    use_trials = true(size(stimIDs));
end

%% Get trial info

trial_info = struct;

if task_dataset
    
    stim_contrastside = signals_events.trialSideValues(1:n_trials)'.* ...
        signals_events.trialContrastValues(1:n_trials)';
    trial_info.stimulus = stim_contrastside(use_trials);
    
    trial_info.response = 3-(abs((trial_choice(use_trials)+1)/2)+1);
    trial_info.repeatNum = ones(sum(use_trials),1);
    
    trial_info.outcome = reshape(trial_outcome(use_trials),[],1);
    
    trial_info.stim_to_move = stim_to_move(use_trials);
    
else % If passive dataset
    
    trial_info.stimulus = stimIDs;
    
end

%% Trial-align data

if verbose; disp('Trial-aligning data...'); end;

% Cortical fluorescence
event_aligned_V = ...
    interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);

% Striatal multiunit
event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
for curr_depth = 1:n_depths
    curr_spikes = spike_times_timeline(depth_group == curr_depth);
    
    % (skip if no spikes at this depth)
    if isempty(curr_spikes)
        continue
    end
    
    event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
        histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
        [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
end

% (filter MUA if selected)
if filter_mua
    event_aligned_mua = AP_deconv_wf(event_aligned_mua,true);
end

% Wheel velocity
event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
    wheel_velocity,t_peri_event);

% Facecam movement
if exist('frame_movement','var')
    event_aligned_movement = interp1(facecam_t(~isnan(facecam_t)), ...
        frame_movement(~isnan(facecam_t)),t_peri_event);
end

if task_dataset
    
% Outcome (reward page 1, punish page 2)
% (note incorrect outcome imprecise from signals, but looks good)
event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);

event_aligned_outcome(trial_outcome == 1,:,1) = ...
    (cell2mat(arrayfun(@(x) ...
    histcounts(reward_t_timeline,t_peri_event_bins(x,:)), ...
    find(trial_outcome == 1),'uni',false))) > 0;

event_aligned_outcome(trial_outcome == -1,:,2) = ...
    (cell2mat(arrayfun(@(x) ...
    histcounts(signals_events.responseTimes,t_peri_event_bins(x,:)), ...
    find(trial_outcome == -1),'uni',false))) > 0;

end

%% Regress cortex to striatum and trial-align

if verbose; disp('Regressing cortex to striatum...'); end;

% Parameters for regression
regression_params.use_svs = 1:200;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.1,0.1];
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

% (filter MUA if selected)
if filter_mua
    binned_spikes = AP_deconv_wf(binned_spikes,true);
    binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
end

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

% Regress cortex to striatum
[ctx_str_k,ctxpred_spikes_std,explained_var] = ...
    AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
    binned_spikes_std,kernel_frames,lambda, ...
    regression_params.zs,regression_params.cvfold, ...
    true,regression_params.use_constant);

% Recast the k's into the master U
ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
    reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
    size(ctx_str_k{1}));

% Re-scale the prediction (subtract offset, multiply, add scaled offset)
ctxpred_spikes = (ctxpred_spikes_std - squeeze(ctx_str_k{end})).* ...
    nanstd(binned_spikes,[],2) + ...
    nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});

event_aligned_mua_ctxpred = ...
    interp1(time_bin_centers,ctxpred_spikes',t_peri_event)./raster_sample_rate;


%% Regress cortex to wheel velocity and speed and trial-align

if verbose; disp('Regressing cortex to wheel...'); end;

wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);

[ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
    AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
    wheel_velspeed_resample_std,kernel_frames,lambda, ...
    regression_params.zs,regression_params.cvfold, ...
    false,false);

predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
    std(wheel_velocity_resample);

% Recast the k's into the master U
ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
    reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
    size(ctx_wheel_k));

event_aligned_wheel_ctxpred = ...
    interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);


%% Regress task to cortex/striatum/cortex-predicted striatum

if task_dataset
    
    if verbose; disp('Regressing task to neural data...'); end;
    
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
    move_onset_regressors = zeros(2,length(time_bin_centers));
    move_onset_regressors(1,:) = histcounts(wheel_move_time(trial_choice == -1),time_bins);
    move_onset_regressors(2,:) = histcounts(wheel_move_time(trial_choice == 1),time_bins);
    
    % Go cue regressors - separate for early/late move
    % (using signals timing - not precise but looks good)
    % (for go cue only on late move trials)
    %         go_cue_regressors = histcounts( ...
    %             signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
    % (for go cue with early/late move trials)
    go_cue_regressors = zeros(1,length(time_bin_centers));
    go_cue_regressors(1,:) = histcounts( ...
        signals_events.interactiveOnTimes(stim_to_move <= 0.5),time_bins);
    go_cue_regressors(2,:) = histcounts( ...
        signals_events.interactiveOnTimes(stim_to_move > 0.5),time_bins);
    
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
    
    % Regression task -> MUA
    baseline = nanmean(reshape(event_aligned_mua(:,t < 0,:),[], ...
        size(event_aligned_mua,3))*raster_sample_rate,1)';
    activity = single(binned_spikes) - baseline;
    
    [mua_taskpred_k,mua_taskpred_long,mua_taskpred_expl_var,mua_taskpred_reduced_long] = ...
        AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
        lambda,zs,cvfold,return_constant,use_constant);
    
    mua_taskpred = ...
        interp1(time_bin_centers,mua_taskpred_long',t_peri_event)./raster_sample_rate;
    
    mua_taskpred_reduced = cell2mat(arrayfun(@(x) ...
        interp1(time_bin_centers,mua_taskpred_reduced_long(:,:,x)', ...
        t_peri_event)./raster_sample_rate,permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
    
    % Regression task -> MUA-ctxpred
    baseline = nanmean(reshape(event_aligned_mua_ctxpred(:,t < 0,:),[], ...
        size(event_aligned_mua_ctxpred,3))*raster_sample_rate,1)';
    activity = single(ctxpred_spikes) - baseline;
    
    [mua_ctxpred_taskpred_k,mua_ctxpred_taskpred_long,mua_ctxpred_taskpred_expl_var,mua_ctxpred_taskpred_reduced_long] = ...
        AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
        lambda,zs,cvfold,return_constant,use_constant);
    
    mua_ctxpred_taskpred = ...
        interp1(time_bin_centers,mua_ctxpred_taskpred_long',t_peri_event)./raster_sample_rate;
    
    mua_ctxpred_taskpred_reduced = cell2mat(arrayfun(@(x) ...
        interp1(time_bin_centers,mua_ctxpred_taskpred_reduced_long(:,:,x)', ...
        t_peri_event)./raster_sample_rate,permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
    
    % Regression task -> (master U, deconvolved) fluor
    event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
    fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);
    
    baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,:),[],size(event_aligned_V_deconv,3)))';
    activity = single(fVdf_deconv_resample_recast(use_components,:))-baseline;
    
    [fluor_taskpred_k,fluor_taskpred_long,fluor_taskpred_expl_var,fluor_taskpred_reduced_long] = ...
        AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
        lambda,zs,cvfold,return_constant,use_constant);
    
    fluor_taskpred = ...
        interp1(time_bin_centers,fluor_taskpred_long',t_peri_event);
    
    fluor_taskpred_reduced = cell2mat(arrayfun(@(x) ...
        interp1(time_bin_centers,fluor_taskpred_reduced_long(:,:,x)', ...
        t_peri_event),permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
    
end

%% Store everything into structure

trial_data = struct;

trial_data.trial_info_all = trial_info;

trial_data.fluor_all = event_aligned_V(use_trials,:,:,:);
trial_data.mua_all = event_aligned_mua(use_trials,:,:,:);

trial_data.ctx_str_k_all = ctx_str_k_recast;
trial_data.mua_ctxpred_all = event_aligned_mua_ctxpred(use_trials,:,:,:);

trial_data.wheel_all = event_aligned_wheel(use_trials,:,:);
if exist('event_aligned_movement','var')
    trial_data.movement_all = event_aligned_movement(use_trials,:,:);
end

trial_data.ctx_wheel_k_all = ctx_wheel_k_recast;
trial_data.wheel_ctxpred_all = event_aligned_wheel_ctxpred(use_trials,:,:);

if task_dataset
    
    trial_data.outcome_all = event_aligned_outcome(use_trials,:,:);
    
    trial_data.mua_taskpred_k_all = mua_taskpred_k;
    trial_data.mua_taskpred_all = mua_taskpred(use_trials,:,:,:);
    trial_data.mua_taskpred_reduced_all = mua_taskpred_reduced(use_trials,:,:,:);
    trial_data.mua_taskpred_expl_var_total_all = mua_taskpred_expl_var.total;
    trial_data.mua_taskpred_expl_var_partial_all = mua_taskpred_expl_var.partial;
    
    trial_data.mua_ctxpred_taskpred_k_all = mua_ctxpred_taskpred_k;
    trial_data.mua_ctxpred_taskpred_all = mua_ctxpred_taskpred(use_trials,:,:,:);
    trial_data.mua_ctxpred_taskpred_reduced_all = mua_ctxpred_taskpred_reduced(use_trials,:,:,:);
    trial_data.mua_ctxpred_taskpred_expl_var_total_all = mua_ctxpred_taskpred_expl_var.total;
    trial_data.mua_ctxpred_taskpred_expl_var_partial_all = mua_ctxpred_taskpred_expl_var.partial;
    
    trial_data.fluor_taskpred_k_all = fluor_taskpred_k;
    trial_data.fluor_taskpred_all = fluor_taskpred(use_trials,:,:,:);
    trial_data.fluor_taskpred_reduced_all = fluor_taskpred_reduced(use_trials,:,:,:);
    trial_data.fluor_taskpred_expl_var_total_all = fluor_taskpred_expl_var.total;
    trial_data.fluor_taskpred_expl_var_partial_all = fluor_taskpred_expl_var.partial;
    
end

if verbose; disp('Done getting trial data.'); end;





