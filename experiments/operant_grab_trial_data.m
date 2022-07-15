% operant_grab_trial_data
% pull out data organized by trials

%% Set flags

do_regression = false; % manually set

% Task dataset if signals (expDef) and includes task expDef
task_dataset = exist('expDef','var') && contains(expDef,'stimWheel');

% Data modalities
imaging_exists = experiments(curr_day).imaging;
ephys_exists = experiments(curr_day).ephys;


%% Set parameters for cortical fluoresence

if imaging_exists
    
if verbose; disp('Re-casting and deconvolving fluorescence...'); end;

% Deconvolve fluoresence (for use in regression)
fVdf_deconv = AP_deconv_wf(fVdf);
fVdf_deconv(isnan(fVdf_deconv)) = 0;

% Convert U to master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
Udf_aligned = AP_align_widefield(Udf,animal,day);
fVdf_deconv_recast = ChangeU(Udf_aligned,fVdf_deconv,U_master);

% Set components to keep
use_components = 1:200;

end

%% Get trial info

% Pick trials to keep (currently: all of them)
if task_dataset
    use_trials = true(n_trials,1);
else
    use_trials = true(size(stimIDs));
end

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

%% Get whisker movement
% (not general enough to go in AP_load_experiment)

if verbose; disp('Getting whisker movements...'); end;

% Get aligned whisker mask for animal/day
facecam_processing_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\facecam_processing';
facecam_align_fn = fullfile(facecam_processing_path,'facecam_align.mat');
load(facecam_align_fn);

facecam_align_animalidx = find(strcmp(animal,{facecam_align.animal}));
facecam_align_dayidx = find(strcmp(day,[facecam_align(facecam_align_animalidx).day]));

try
    whisker_mask = facecam_align(facecam_align_animalidx).whisker_mask{facecam_align_dayidx};

    % Get whisker movement
    if ~isempty(whisker_mask)
        vr = VideoReader(facecam_fn);
        n_frames = vr.NumFrames;

        whisker_px = nan(n_frames,sum(whisker_mask(:)));
        for curr_frame = 1:n_frames
            curr_im = read(vr,curr_frame);
            whisker_px(curr_frame,:) = curr_im(whisker_mask);
        end

        whisker_move = [NaN;sum(abs(diff(whisker_px,[],1)),2)];
    end
end

%% Set parameters for trials

% Get event-aligned activity
raster_window = [-0.5,1];

% (either imaging or ephys, set separate sample rates)
if imaging_exists
    raster_sample_rate = 50;
elseif ephys_exists
    raster_sample_rate = 50;
end

raster_sample_time = 1/raster_sample_rate;
t = raster_window(1):raster_sample_time:raster_window(2);


%% Set align times (stimulus)

stim_align = stimOn_times;
t_peri_stim = bsxfun(@plus,stim_align,t);
t_peri_stim_bins = [t_peri_stim-raster_sample_time/2,t_peri_stim(:,end)+raster_sample_time/2];

%% Set align times (movement - first post-stim)

if task_dataset
    move_align = wheel_starts(wheel_move_stim_idx);
    t_peri_move = bsxfun(@plus,move_align,t);
    t_peri_move_bins = [t_peri_move-raster_sample_time/2,t_peri_move(:,end)+raster_sample_time/2];
end

%% Set align times (delay-period movement)
% (movements without stim on screen that hit reward threshold)

if task_dataset

    if verbose; disp('Getting rewardable movements without stimuli...'); end;

    % Grab all movements separately (in degrees from move start)
    wheel_moves_deg = arrayfun(@(x) wheel_position_deg( ...
        Timeline.rawDAQTimestamps >= wheel_starts(x) & ...
        Timeline.rawDAQTimestamps <= wheel_stops(x)) - ...
        wheel_position_deg(find(Timeline.rawDAQTimestamps >= wheel_starts(x),1)), ...
        1:length(wheel_starts),'uni',false);

    % Find movements that hit reward limit (and don't hit punish limit)
    deg_reward = -90;
    deg_punish = 90;

    wheel_moves_deg_rewardlimit = find(cellfun(@(x) ...
        any(x <= deg_reward) && ...
        ~any(x >= deg_punish),wheel_moves_deg));

    wheel_move_nostim_rewardable_idx = ...
        intersect(wheel_move_nostim_idx,wheel_moves_deg_rewardlimit);

    % Less strict: find movements that are just net leftwards
    wheel_move_leftward_idx = find(cellfun(@sum,wheel_moves_deg) < 0);
    wheel_move_nostim_leftward_idx = ...
        intersect(wheel_move_nostim_idx,wheel_move_leftward_idx);

    % Get move no-stim align times
    move_nostim_rewardable_align = wheel_starts(wheel_move_nostim_rewardable_idx);
    t_peri_move_nostim_rewardable = bsxfun(@plus,move_nostim_rewardable_align,t);
    t_peri_move_nostim_rewardable_bins = [t_peri_move_nostim_rewardable-raster_sample_time/2,t_peri_move_nostim_rewardable(:,end)+raster_sample_time/2];

    move_nostim_leftward_align = wheel_starts(wheel_move_nostim_leftward_idx);
    t_peri_move_nostim_leftward = bsxfun(@plus,move_nostim_leftward_align,t);
    t_peri_move_nostim_leftward_bins = [t_peri_move_nostim_leftward-raster_sample_time/2,t_peri_move_nostim_leftward(:,end)+raster_sample_time/2];

end

%% Grab aligned data

if verbose; disp('Trial-aligning data...'); end;

% Cortical fluorescence
if imaging_exists
    stim_aligned_V_deconv = ...
        interp1(frame_t,fVdf_deconv_recast(use_components,:)',t_peri_stim,'previous');

    if task_dataset
        move_aligned_V_deconv = ...
            interp1(frame_t,fVdf_deconv_recast(use_components,:)',t_peri_move,'previous');

        move_nostim_rewardable_aligned_V_deconv = ...
            interp1(frame_t,fVdf_deconv_recast(use_components,:)',t_peri_move_nostim_rewardable,'previous');

        move_nostim_leftward_aligned_V_deconv = ...
            interp1(frame_t,fVdf_deconv_recast(use_components,:)',t_peri_move_nostim_leftward,'previous');
    end
end

% Cortical electrophysiology (multiunit by area)
% (split by wide/narrow waveform)
if ephys_exists
    
    % (to just pool all spikes)
    use_celltypes = {true(size(unique(spike_templates)))};
    % (to pool by cell type)
%     use_celltypes = {wide,narrow}; 

    stim_aligned_mua_area = nan(size(t_peri_stim_bins,1), ...
        length(t),length(probe_areas),length(use_celltypes));

    if task_dataset
        move_nostim_rewardable_aligned_mua_area = ...
            nan(size(t_peri_move_nostim_rewardable_bins,1), ...
            length(t),length(probe_areas),length(use_celltypes));
    end

    for curr_area = 1:length(probe_areas)
        for curr_celltype = 1:length(use_celltypes) % (wide, narrow)

            curr_celltype_idx = find(use_celltypes{curr_celltype});

            curr_spikes = spike_times_timeline(...
                spike_depths >= probe_area_boundaries{curr_area}(1) & ...
                spike_depths <= probe_area_boundaries{curr_area}(2) & ...
                ismember(spike_templates,curr_celltype_idx));

            if isempty(curr_spikes)
                continue
            end

            stim_aligned_mua_area(:,:,curr_area,curr_celltype) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_stim_bins(x,:)), ...
                [1:size(t_peri_stim_bins,1)]','uni',false))*raster_sample_rate;

            if task_dataset
                move_nostim_rewardable_aligned_mua_area(:,:,curr_area,curr_celltype) = cell2mat(arrayfun(@(x) ...
                    histcounts(curr_spikes,t_peri_move_nostim_rewardable_bins(x,:)), ...
                    [1:size(t_peri_move_nostim_rewardable_bins,1)]','uni',false))*raster_sample_rate;
            end
        end
    end
   
end

% Wheel velocity
stim_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
    wheel_velocity,t_peri_stim,'previous');

% Pupil diameter
if exist('pupil_diameter','var')
    stim_aligned_pupil= interp1(eyecam_t(~isnan(eyecam_t)), ...
        pupil_diameter(~isnan(eyecam_t)),t_peri_stim,'previous');
end

% Whisker movement
if exist('whisker_move','var')
    stim_aligned_whisker = interp1(facecam_t(~isnan(facecam_t)), ...
        whisker_move(~isnan(facecam_t)),t_peri_stim,'previous');
end

if task_dataset

    % Outcome (reward page 1, punish page 2)
    % (note incorrect outcome imprecise from signals, but looks good)
    stim_aligned_outcome = zeros(size(t_peri_stim,1),size(t_peri_stim,2),2);

    stim_aligned_outcome(trial_outcome == 1,:,1) = ...
        (cell2mat(arrayfun(@(x) ...
        histcounts(reward_t_timeline,t_peri_stim_bins(x,:)), ...
        find(trial_outcome == 1),'uni',false))) > 0;

    stim_aligned_outcome(trial_outcome == -1,:,2) = ...
        (cell2mat(arrayfun(@(x) ...
        histcounts(signals_events.responseTimes,t_peri_stim_bins(x,:)), ...
        find(trial_outcome == -1),'uni',false))) > 0;

end


%% Task regression: set time points and grab continuous data

if task_dataset && do_regression
    
    % Parameters for regression
    regression_params.use_svs = 1:100;
    regression_params.skip_seconds = 20;
    regression_params.sample_rate = raster_sample_rate;
    regression_params.kernel_t = [-0.1,0.1];
    regression_params.zs = [false,false];
    regression_params.cvfold = 5;
    regression_params.use_constant = true;
    
    if imaging_exists
        
        % Continuous time bins based on frame start/end
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/regression_params.sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Resample deconvolved fluorescence
        fVdf_deconv_recast_resample = interp1(frame_t,fVdf_deconv_recast',time_bin_centers,'previous')';
        
    elseif ~imaging_exists && ephys_exists
        
        % Continuous time bins based on experiment start/end
        time_bins = (acqLive_timeline(1) + regression_params.skip_seconds): ...
            1/regression_params.sample_rate: ...
            (acqLive_timeline(end) - regression_params.skip_seconds);
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Bin continuous muliunit by area
        binned_spikes = nan(length(probe_areas),length(time_bin_centers));
        for curr_area = 1:length(probe_areas)
            curr_spikes = spike_times_timeline(...
                spike_depths >= probe_area_boundaries{curr_area}(1) & ...
                spike_depths <= probe_area_boundaries{curr_area}(2));
            
            if isempty(curr_spike_times)
                continue
            end
            
            binned_spikes(curr_area,:) = ...
                histcounts(curr_spikes,time_bins)*raster_sample_rate;
        end
                
    end
    
end

%% Task regression: execute

if task_dataset && do_regression
    
    if verbose; disp('Regressing task to neural data...'); end;
    
    % Build regressors (only a subset of these are used)
    
    % Stim regressors
    stim_contrastsides = ...
        signals_events.trialSideValues(1:length(stimOn_times))'.* ...
        signals_events.trialContrastValues(1:length(stimOn_times))';
    % (only 1 stim in AP_stimWheelRight/Left)
    unique_stim = unique(stim_contrastsides);
    
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
    
    % Move onset regressors (L/R)
%     move_onset_regressors = zeros(2,length(time_bin_centers));
%     move_onset_regressors(1,:) = histcounts(wheel_move_time(trial_choice == -1),time_bins);
%     move_onset_regressors(2,:) = histcounts(wheel_move_time(trial_choice == 1),time_bins);    
    
    % (trying new thing: all rewarded/unrewarded move onsets)
    wheel_move_resample = interp1(Timeline.rawDAQTimestamps,+wheel_move,time_bin_centers,'previous');
    
    move_onset_resample_t = time_bin_centers([false,diff(wheel_move_resample) == 1]);
    move_offset_resample_t = time_bin_centers([false,diff(wheel_move_resample) == -1]);   
    
    % (rewarded movements: which comes first, a reward or move offset)
    use_move_offsets = ~ismember(move_offset_resample_t,reward_t_timeline) & ...
        ~ismember(move_offset_resample_t,time_bin_centers(end));
    rewarded_movements = logical(interp1( ...
        [reward_t_timeline,move_offset_resample_t(use_move_offsets),time_bin_centers(end)], ...
        [ones(size(reward_t_timeline)),zeros(1,sum(use_move_offsets)),0], ...
        move_onset_resample_t,'next','extrap'));
    
    move_onset_regressors = +[... 
        ismember(time_bin_centers,move_onset_resample_t(rewarded_movements)); ...
        ismember(time_bin_centers,move_onset_resample_t(~rewarded_movements))];   
    
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
    
    % Trial offset regressors
    % (separate constant offset "baseline" regressor for each trial)
    % (starts 0.5s before stim which is min quiescence time)
    % (--> didn't use this, made no obvious difference)
    trial_prestim_time = -0.5;
    trial_regressor_start = (time_bin_centers - stimOn_times) >= trial_prestim_time;
    trial_regressor_end = padarray(time_bin_centers - stimOn_times(2:end) < ...
        trial_prestim_time,[1,0],true,'post');    
    trial_offset_regressors = +(trial_regressor_start & trial_regressor_end);
    
    % Concatenate selected regressors, set parameters
    task_regressors = {stim_regressors;move_onset_regressors;outcome_regressors};
    task_regressor_labels = {'Stim','Move onset','Outcome'};
    
    task_t_shifts = { ...
        [0,1]; ... % stim
        [-0.5,0.5]; ... % move
        [0,1]}; % outcome
    
    task_regressor_sample_shifts = cellfun(@(x) round(x(1)*(regression_params.sample_rate)): ...
        round(x(2)*(regression_params.sample_rate)),task_t_shifts,'uni',false);
    lambda = 0;
    zs = [false,false];
    cvfold = 5;
    use_constant = false;
    return_constant = false;
   
    % Regress task to fluor (master U, deconvolved)
    if imaging_exists
        baseline = nanmean(reshape(stim_aligned_V_deconv(:,t < 0,:),[],size(stim_aligned_V_deconv,3)))';
        activity = single(fVdf_deconv_recast_resample(use_components,:))-baseline;
        
        [fluor_taskpred_k,fluor_taskpred_long,fluor_taskpred_expl_var,fluor_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        fluor_taskpred = ...
            interp1(time_bin_centers,fluor_taskpred_long',t_peri_stim,'previous');
        
        fluor_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,fluor_taskpred_reduced_long(:,:,x)', ...
            t_peri_stim,'previous'),permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
    end
    
    % Regress task to multiunit
    if ephys_exists
        
        baseline = nanmean(reshape(stim_aligned_mua_area(:,t < 0,:),[], ...
            size(stim_aligned_mua_area,3)),1)';
        activity = single(binned_spikes) - baseline;
        
        [mua_taskpred_k,mua_taskpred_long,mua_taskpred_expl_var,mua_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        mua_taskpred = ...
            interp1(time_bin_centers,mua_taskpred_long',t_peri_stim);
        
        mua_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,mua_taskpred_reduced_long(:,:,x)', ...
            t_peri_stim),permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
    end
    
end

%% Store everything into structure

trial_data = struct;

% Task/behavioral data
trial_data.trial_info_all = trial_info;
trial_data.wheel_all = stim_aligned_wheel(use_trials,:);
if exist('stim_aligned_whisker','var')
    trial_data.whisker_move_all = stim_aligned_whisker(use_trials,:);
end
if exist('stim_aligned_pupil','var')
    trial_data.pupil_diameter_all = stim_aligned_pupil(use_trials,:);
end

% Aligned neural data
if imaging_exists
    trial_data.fluor_all = stim_aligned_V_deconv(use_trials,:,:,:);
    if task_dataset
        trial_data.fluor_move_all = move_aligned_V_deconv(use_trials,:,:,:);

        % (move no-stim: just average)
        trial_data.fluor_move_nostim_rewardable_all = nanmean(move_nostim_rewardable_aligned_V_deconv,1);
        trial_data.fluor_move_nostim_leftward_all = nanmean(move_nostim_leftward_aligned_V_deconv,1);
    end
end
if ephys_exists
    trial_data.probe_areas_all = probe_areas; 
    trial_data.mua_area_all = stim_aligned_mua_area(use_trials,:,:,:);
    if task_dataset
        trial_data.mua_area_move_nostim_rewardable_all = nanmean(move_nostim_rewardable_aligned_mua_area,1);
    end
end

% Task information
if task_dataset
    trial_data.outcome_all = stim_aligned_outcome(use_trials,:,:);
end

% Task regression
if task_dataset && do_regression
    if imaging_exists
        trial_data.fluor_taskpred_k_all = fluor_taskpred_k;
        trial_data.fluor_taskpred_all = fluor_taskpred(use_trials,:,:,:);
        trial_data.fluor_taskpred_reduced_all = fluor_taskpred_reduced(use_trials,:,:,:);
        trial_data.fluor_taskpred_expl_var_total_all = fluor_taskpred_expl_var.total;
        trial_data.fluor_taskpred_expl_var_partial_all = fluor_taskpred_expl_var.partial;
    end

    if ephys_exists
        trial_data.mua_taskpred_k_all = mua_taskpred_k;
        trial_data.mua_taskpred_all = mua_taskpred(use_trials,:,:,:);
        trial_data.mua_taskpred_reduced_all = mua_taskpred_reduced(use_trials,:,:,:);
        trial_data.mua_taskpred_expl_var_total_all = mua_taskpred_expl_var.total;
        trial_data.mua_taskpred_expl_var_partial_all = mua_taskpred_expl_var.partial;
    end
end
    

if verbose; disp('Done getting trial data.'); end;





