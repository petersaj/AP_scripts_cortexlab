% AP_load_experiment(animal,day,experiment,site)
%
% Loads data from experiments


%% Display progress or not
if ~exist('verbose','var')
    verbose = false;
end

%% Define what to load

% Site is optional
if ~exist('site','var')
    site = [];
end

% If nothing specified, load everything
if ~exist('load_parts','var')
    load_parts.cam = true;
    load_parts.imaging = true;
    load_parts.ephys = true;
else
    % If only some things specified, don't load others
    if ~isfield(load_parts,'cam')
        load_parts.cam = false;
    end
    if ~isfield(load_parts,'imaging')
        load_parts.imaging = false;
    end
    if ~isfield(load_parts,'ephys')
        load_parts.ephys = false;
    end
end

%% Load timeline and associated inputs

[timeline_filename,timeline_exists] = AP_cortexlab_filename(animal,day,experiment,'timeline');
if ~timeline_exists
    error([animal ' ' day ' ' num2str(experiment) ': no timeline']);
end

if timeline_exists
    if verbose; disp('Loading timeline...'); end
    
    load(timeline_filename);
       
    % Get camera times
    cam_name = 'pcoExposure';
    timeline_cam_idx = strcmp({Timeline.hw.inputs.name}, cam_name);
    
    cam_expose_starts = Timeline.rawDAQTimestamps( ...
        find(Timeline.rawDAQData(1:end-1,timeline_cam_idx) <= 2 & ...
        Timeline.rawDAQData(2:end,timeline_cam_idx) > 2) + 1);
    cam_expose_stops = Timeline.rawDAQTimestamps( ...
        find(Timeline.rawDAQData(1:end-1,timeline_cam_idx) >= 2 & ...
        Timeline.rawDAQData(2:end,timeline_cam_idx) < 2) + 1);
    
    cam_time = cam_expose_starts;
    cam_expose_times = cam_expose_stops - cam_expose_starts;
    
    % Get acqLive signal
    acqLive_name = 'acqLive';
    acqLive_idx = strcmp({Timeline.hw.inputs.name}, acqLive_name);
    thresh = max(Timeline.rawDAQData(:,acqLive_idx))/2;
    acqLive_trace = Timeline.rawDAQData(:,acqLive_idx) > thresh;
    acqLive_timeline = Timeline.rawDAQTimestamps( ...
        [find(acqLive_trace,1),find(acqLive_trace,1,'last')+1]);
    
    % Get wheel position
    rotaryEncoder_idx = strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder');
    % (this is a very strange hack to overcome a problem in the rotary
    % encoder that's known in the lab and was put on the wiki)
    wheel_position = Timeline.rawDAQData(:,rotaryEncoder_idx);
    wheel_position(wheel_position > 2^31) = wheel_position(wheel_position > 2^31) - 2^32;
    
    % Get wheel velocity by smoothing the wheel trace and taking deriv
    wheel_smooth_t = 0.05; % seconds
    wheel_smooth_samples = wheel_smooth_t/Timeline.hw.samplingInterval;
    wheel_velocity = interp1(conv(Timeline.rawDAQTimestamps,[1,1]/2,'valid'), ...
        diff(smooth(wheel_position,wheel_smooth_samples)),Timeline.rawDAQTimestamps)';
       
    % Get whether stim was flickering
    stimScreen_idx = strcmp({Timeline.hw.inputs.name}, 'stimScreen');
    if any(stimScreen_idx)
        stimScreen_flicker = max(Timeline.rawDAQData(:,stimScreen_idx)) - ...
            min(Timeline.rawDAQData(:,stimScreen_idx)) > 2;
    end
    
    % Get photodiode flips (compensate for screen flicker)
    photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
    % (define stim screen on from photodiode - sometimes sample-length
    % offset maybe because of backlight onset delay)
    stimScreen_on = Timeline.rawDAQData(:,photodiode_idx) > 0.15;
    stimScreen_on_t = Timeline.rawDAQTimestamps(stimScreen_on);
    photodiode_thresh = 2; % old: max(Timeline.rawDAQData(:,photodiode_idx))/2
    photodiode_trace = Timeline.rawDAQData(stimScreen_on,photodiode_idx) > photodiode_thresh;
    % (medfilt because photodiode can be intermediate value when backlight
    % coming on)
    photodiode_trace_medfilt = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
        photodiode_idx),3) > photodiode_thresh;
    photodiode_flip = find((~photodiode_trace_medfilt(1:end-1) & photodiode_trace_medfilt(2:end)) | ...
        (photodiode_trace_medfilt(1:end-1) & ~photodiode_trace_medfilt(2:end)))+1;
    photodiode_flip_times = stimScreen_on_t(photodiode_flip)';
    
    % Get flipper signal (this was added late, might not be present)
    flipper_name = 'flipper';
    flipper_idx = strcmp({Timeline.hw.inputs.name}, flipper_name);
    flipper_thresh = 2; % TTL threshold
    flipper_trace = Timeline.rawDAQData(:,flipper_idx) > flipper_thresh;
    flipper_flip = find((~flipper_trace(1:end-1) & flipper_trace(2:end)) | ...
        (flipper_trace(1:end-1) & ~flipper_trace(2:end)))+1;
    flipper_flip_times_timeline = Timeline.rawDAQTimestamps(flipper_flip)';
    
end

%% Load mpep protocol

[protocol_filename,protocol_exists] = AP_cortexlab_filename(animal,day,experiment,'protocol');

if protocol_exists
    
    if verbose; disp('Loading mpep protocol...'); end
    
    load(protocol_filename);
    
    % Load in hardware info
    hwinfo_filename = AP_cortexlab_filename(animal,day,experiment,'hardware');
    load(hwinfo_filename);
        
    % Stim times should just be odd (on) and even (off)
    if mod(length(photodiode_flip_times),2) == 0
        photodiode_onsets = photodiode_flip_times(1:2:end);
        photodiode_offsets = photodiode_flip_times(2:2:end);
    else
        error('Odd number of photodiode flips')
    end
    
    % Get flicker/steady photodiode mode
    photodiode_type = lower(myScreenInfo.SyncSquare.Type);
    
    % Get stim on times
    if strcmp(Protocol.xfile,'stimSparseNoiseUncorrAsync.x')
        % Sparse noise
        % (at the moment this is in a sparse noise-specific script)       
    else
        % Anything else
        if length(photodiode_onsets) == numel(Protocol.seqnums)
            % If photodiode onsets matches number of stimuli, use those
            stimOn_times = photodiode_onsets;
        else
            % Get specific stim onsets by time between last offset and new onset
            % (occasionally there a bad frame so flip but not new stim)
            refresh_rate_cutoff = 1/5;
            stimOn_times = photodiode_onsets( ...
                [1;find(photodiode_onsets(2:end) - photodiode_offsets(1:end-1) > refresh_rate_cutoff) + 1]);
            
            if length(stimOn_times) ~= numel(Protocol.seqnums)
                error('MPEP/Photodiode error: photodiode doesn''t match stim')
            end
        end
        stimIDs = zeros(size(stimOn_times));
        for q = 1:size(Protocol.seqnums,1)
            stimIDs(Protocol.seqnums(q,:)) = q;
        end
    end
    
end

%% Load task/behavior

% Load the block
[block_filename, block_exists] = AP_cortexlab_filename(animal,day,experiment,'block');

if block_exists
    
    if verbose; disp('Loading block file...'); end
    
    load(block_filename);
    
    signals_events = block.events;
    
    % If reward information exists, use that to align signals/timeline
    % (bad now because manual reward possible - use flipper in future)
    if exist('Timeline','var') && isfield(block.outputs,'rewardTimes')
        reward_t_block = block.outputs.rewardTimes(block.outputs.rewardValues > 0);
        
        timeline_reward_idx = strcmp({Timeline.hw.inputs.name}, 'rewardEcho');
        reward_thresh = max(Timeline.rawDAQData(:,timeline_reward_idx))/2;
        reward_trace = Timeline.rawDAQData(:,timeline_reward_idx) > reward_thresh;
        reward_t_timeline = Timeline.rawDAQTimestamps(find(reward_trace(2:end) & ~reward_trace(1:end-1))+1);
        
        % If there's a different number of block and timeline rewards (aka
        % manual rewards were given), try to fix this
        if length(reward_t_block) ~= length(reward_t_timeline)
            % (this is really inelegant but I think works - find the most
            % common offset between block/timeline rewards)
            reward_t_offset = bsxfun(@minus,reward_t_block',reward_t_timeline);
            blunt_reward_offset = mode(round(reward_t_offset(:)*10))/10;
            reward_t_offset_shift = reward_t_offset - blunt_reward_offset;
            t_offset_tolerance = 0.1;
            reward_t_offset_binary = abs(reward_t_offset_shift) < t_offset_tolerance;
            if all(sum(reward_t_offset_binary,2) == 1)
                % one timeline reward for each block reward, you're good
                % (eliminate the timeline rewards with no match)
                manual_timeline_rewards = sum(reward_t_offset_binary,1) == 0;
                reward_t_timeline(manual_timeline_rewards) = [];
                warning('Manual rewards included - removed successfully');
            else
                % otherwise, you're in trouble
                error('Manual rewards included - couldn''t match to block');
            end
        end
        
        % Go through all block events and convert to timeline time
        % (uses reward as reference)
        block_fieldnames = fieldnames(block.events);
        block_values_idx = cellfun(@(x) ~isempty(x),strfind(block_fieldnames,'Values'));
        block_times_idx = cellfun(@(x) ~isempty(x),strfind(block_fieldnames,'Times'));
        for curr_times = find(block_times_idx)'
            if isempty(signals_events.(block_fieldnames{curr_times}))
                % skip if empty
                continue
            end
            signals_events.(block_fieldnames{curr_times}) = ...
                interp1(reward_t_block,reward_t_timeline,block.events.(block_fieldnames{curr_times}),'linear','extrap');
        end
    end   
    
    % SPECIFIC TO PROTOCOL
    [~,expDef] = fileparts(block.expDef);
    if strcmp(expDef,'vanillaChoiceworld')       
        
        % Hit/miss recorded for last trial, circshift to align
        signals_events.hitValues = circshift(signals_events.hitValues,[0,-1]);
        signals_events.missValues = circshift(signals_events.missValues,[0,-1]);
        
        % Get number of completed trials (if uncompleted last trial)
        n_trials = length(signals_events.endTrialTimes);
        
        % Get stim on times by closest photodiode flip
        [~,closest_stimOn_photodiode] = ...
            arrayfun(@(x) min(abs(signals_events.stimOnTimes(x) - ...
            photodiode_flip_times)), ...
            1:n_trials);
        stimOn_times = photodiode_flip_times(closest_stimOn_photodiode);
        
        % Get time from stim on to first wheel movement     
        surround_time = [-0.5,2];
        surround_samples = surround_time/Timeline.hw.samplingInterval;
     
        surround_time = surround_time(1):Timeline.hw.samplingInterval:surround_time(2);
        pull_times = bsxfun(@plus,stimOn_times,surround_time);
        
        stim_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,pull_times);
        
        thresh_displacement = 0.02;
        [~,wheel_move_sample] = max(abs(stim_aligned_wheel) > thresh_displacement,[],2);
        wheel_move_time = arrayfun(@(x) pull_times(x,wheel_move_sample(x)),1:size(pull_times,1))';
        wheel_move_time(wheel_move_sample == 1) = NaN;
        
        % Get conditions for all trials      
        
        % (trial_timing)
        stim_to_move = padarray(wheel_move_time - stimOn_times,[n_trials-length(stimOn_times),0],NaN,'post');
        stim_to_feedback = signals_events.responseTimes(1:n_trials)' - stimOn_times(1:n_trials);
        
        % (early vs late move)
        trial_timing = 1 + (stim_to_move > 0.5);
        
        % (choice and outcome)
        go_left = (signals_events.trialSideValues == 1 & signals_events.hitValues == 1) | ...
            (signals_events.trialSideValues == -1 & signals_events.missValues == 1);
        go_right = (signals_events.trialSideValues == -1 & signals_events.hitValues == 1) | ...
            (signals_events.trialSideValues == 1 & signals_events.missValues == 1);
        trial_choice = go_right' - go_left';
        trial_outcome = signals_events.hitValues(1:n_trials)'-signals_events.missValues(1:n_trials)';
        
        % (trial conditions: [contrast,side,choice,timing])
        contrasts = [0,0.06,0.125,0.25,0.5,1];
        sides = [-1,1];
        choices = [-1,1];
        timings = [1,2];
        
        conditions = combvec(contrasts,sides,choices,timings)';
        n_conditions = size(conditions,1);
        
        trial_conditions = ...
            [signals_events.trialContrastValues(1:n_trials)', signals_events.trialSideValues(1:n_trials)', ...
            trial_choice(1:n_trials), trial_timing(1:n_trials)];
        [~,trial_id] = ismember(trial_conditions,conditions,'rows');
        
    elseif strcmp(expDef,'AP_visAudioPassive')
        %         min_stim_downtime = 0.5; % minimum time between pd flips to get stim
        %         stimOn_times_pd = photodiode_flip_times([true;diff(photodiode_flip_times) > min_stim_downtime]);
        %         stimOff_times_pd = photodiode_flip_times([diff(photodiode_flip_times) > min_stim_downtime;true]);
        %         warning('visAudioPassive: THIS IS TEMPORARY BECAUSE NO BUFFER TIME')
        %
        %         stimOn_times = nan(size(signals_events.visualOnsetTimes));
        %         stimOn_times(end-(length(stimOn_times_pd)-1):end) = stimOn_times_pd;
        %
        %         stimOff_times = nan(size(signals_events.visualOnsetTimes));
        %         stimOff_times(end-(length(stimOff_times_pd)-1):end) = stimOff_times_pd;
        %
        %
        %         % sanity check
        %         if length(signals_events.visualOnsetValues) ~= length(stimOn_times)
        %             error('Different number of signals/timeline stim ons')
        %         end
        error('AP_visAudioPassive isn''t reliable yet')
        
    elseif strcmp(expDef,'AP_choiceWorldStimPassive')
        % This is kind of a dumb hack to get the stimOn times, maybe not
        % permanent unless it works fine: get stim times by checking for
        % close to the median photodiode flip difference
        block_stim_iti = mean(diff(block.stimWindowUpdateTimes));
        
        photodiode_flip_diff = diff(stimScreen_on_t(photodiode_flip));
        median_photodiode_flip_diff = mode(round(photodiode_flip_diff*10)/10);
        
        stimOn_idx = find(abs(photodiode_flip_diff-median_photodiode_flip_diff) < 0.1);
        
        stimOn_times = stimScreen_on_t(photodiode_flip(stimOn_idx))';
        
        % assume the times correspond to the last n values (this is because
        % sometimes if the buffer time wasn't enough, the first stimuli
        % weren't shown or weren't shown completely)
        [conditions,conditions_idx,stimIDs] = unique(signals_events.visualParamsValues(:, ...
            signals_events.visualOnsetValues(end-length(stimOn_times)+1:end))','rows');
        conditions_params = signals_events.visualParamsValues(:,conditions_idx);
        
    elseif strcmp(expDef,'DS_choiceWorldStimPassive')
        % get stim times - first stim photodiode is messed up so throw it out
        stimOn_times = photodiode_flip_times(2:2:end);
        
        % sanity check: times between stim on times in signals
        signals_photodiode_iti_diff = diff(signals_events.stimOnTimes(2:end)) - diff(stimOn_times)';
        if any(signals_photodiode_iti_diff > 0.1)
            error('mismatching signals/photodiode stim ITIs')
        end
        
        % Get stim ID and conditions
        contrasts = unique(signals_events.stimContrastValues);
        azimuths = unique(signals_events.stimAzimuthValues);
        
        conditions = combvec(contrasts,azimuths)';
        n_conditions = size(conditions,1);
        
        trial_conditions = ...
            [signals_events.stimContrastValues(2:end); signals_events.stimAzimuthValues(2:end)]';
        [~,stimIDs] = ismember(trial_conditions,conditions,'rows');
        
    elseif strcmp(expDef,'AP_localize_choiceWorldStimPassive')
        % get stim times - first stim photodiode is messed up so throw it out
        stimOn_times = photodiode_flip_times(2:2:end);
        
        % sanity check: times between stim on times in signals
        signals_photodiode_iti_diff = diff(signals_events.stimOnTimes(2:end)) - diff(stimOn_times)';
        if any(signals_photodiode_iti_diff > 0.1)
            error('mismatching signals/photodiode stim ITIs')
        end
        
        % Get stim ID and conditions
        azimuths = unique(signals_events.stimAzimuthValues);
        altitudes = unique(signals_events.stimAltitudeValues);

        trial_conditions = reshape(signals_events.visualParamsValues,2,[])';
        
        conditions = unique(trial_conditions,'rows');
        n_conditions = size(conditions,1);
        
        [~,stimIDs] = ismember(trial_conditions,conditions,'rows');
        
        % Get rid of the first one for now
        trial_conditions = trial_conditions(2:end);
        stimIDs = stimIDs(2:end);
        
    elseif strcmp(expDef,'sparseNoiseAsync_NS2')
        if length(photodiode_flip_times) - length(block.stimWindowRenderTimes) ~= -1
            error('Maybe skipped frames?');
        end
        
        stimOn_times = photodiode_flip_times;
        
    else
        warning(['Signals protocol with no analysis script:' expDef]);
    end
    
    
    
end


%% Load face/eyecam and processing

% Don't load if no timeline
if exist('Timeline','var') && load_parts.cam
    
    % Get cam sync from timeline
    camSync_idx = strcmp({Timeline.hw.inputs.name}, 'camSync');
    camSync_thresh = max(Timeline.rawDAQData(:,camSync_idx))/2;
    camSync = Timeline.rawDAQData(:,camSync_idx) > camSync_thresh;
    camSync_flip = find((camSync(1:end-1) ~= camSync(2:end)))+1;
    if length(camSync_flip) ~= 4
        error('camSync flip number ~= 4')
    end
    
    % EYECAM
    [eyecam_dir,eyecam_exists] = AP_cortexlab_filename(animal,day,experiment,'eyecam');
    
    if eyecam_exists
        if verbose; disp('Loading eyecam...'); end
        
        % Load camera processed data
        [eyecam_processed_filename,eyecam_processed_exists] = AP_cortexlab_filename(animal,day,experiment,'eyecam_processed');
        if eyecam_processed_exists
            eyecam = load(eyecam_processed_filename);
        end
        
        % Get camera times
        eyecam_fn = AP_cortexlab_filename(animal,day,experiment,'eyecam');
        eyecam_dir = fileparts(eyecam_fn);
        eyecam_t_savefile = [eyecam_dir filesep 'eyecam_t.mat'];
        
        if exist(eyecam_fn,'file') && ~exist(eyecam_t_savefile,'file')
            % Get facecam strobes
            eyeCamStrobe_idx = strcmp({Timeline.hw.inputs.name}, 'eyeCameraStrobe') | ...
                strcmp({Timeline.hw.inputs.name}, 'eyeCamStrobe');
            eyeCamStrobe_thresh = max(Timeline.rawDAQData(:,eyeCamStrobe_idx))/5;
            eyeCamStrobe = Timeline.rawDAQData(:,eyeCamStrobe_idx) > eyeCamStrobe_thresh;
            eyeCamStrobe_up = find((~eyeCamStrobe(1:end-1) & eyeCamStrobe(2:end)))+1;
            eyeCamStrobe_up_t = Timeline.rawDAQTimestamps(eyeCamStrobe_up);
            
            % Get sync times for cameras (or load if already done)
            [eyecam_sync_frames,n_eyecam_frames] = AP_get_cam_sync_frames(eyecam_fn);
            
            if ~isempty(eyecam_sync_frames)
                % Get the closest cam strobe to sync start, find offset and frame idx
                [~,eyecam_strobe_sync] = min(abs(camSync_flip(1) - eyeCamStrobe_up));
                eyecam_frame_offset = eyecam_sync_frames(1) - eyecam_strobe_sync;
                eyecam_frame_idx = [1:length(eyeCamStrobe_up)] + eyecam_frame_offset;
                
                % Check that the number of frames between synchs matches
                % video and timeline
                n_eyecam_frames_syncd_movie = diff(eyecam_sync_frames) + 1;
                [~,eyecam_strobe_sync_end] = min(abs(camSync_flip(3) - eyeCamStrobe_up));
                n_eyecam_frames_syncd_timeline = eyecam_strobe_sync_end - eyecam_strobe_sync;
                if abs(n_eyecam_frames_syncd_movie - n_eyecam_frames_syncd_timeline) > 1
                   error('Eyecam: different n frames video vs timeline'); 
                end
                
                % Get times of cam frames in timeline
                eyecam_t = nan(n_eyecam_frames,1);
                eyecam_t(eyecam_frame_idx(eyecam_frame_idx > 0)) = eyeCamStrobe_up_t(eyecam_frame_idx > 0);
                
                save(eyecam_t_savefile,'eyecam_t');
            end
        elseif exist(eyecam_fn,'file') && exist(eyecam_t_savefile,'file')
            load(eyecam_t_savefile);
        end
        
    end
    
    % FACECAM
    [facecam_dir,facecam_exists] = AP_cortexlab_filename(animal,day,experiment,'facecam');
    
    if facecam_exists
        if verbose; disp('Loading facecam...'); end
        
        % Get camera times
        facecam_fn = AP_cortexlab_filename(animal,day,experiment,'facecam');
        facecam_dir = fileparts(facecam_fn);
        facecam_t_savefile = [facecam_dir filesep 'facecam_t.mat'];
        
        if exist(facecam_fn,'file') && ~exist(facecam_t_savefile,'file')
            % Get facecam strobes
            faceCamStrobe_idx = strcmp({Timeline.hw.inputs.name}, 'faceCamStrobe');
            faceCamStrobe_thresh = max(Timeline.rawDAQData(:,faceCamStrobe_idx))/5;
            faceCamStrobe = Timeline.rawDAQData(:,faceCamStrobe_idx) > faceCamStrobe_thresh;
            faceCamStrobe_up = find((~faceCamStrobe(1:end-1) & faceCamStrobe(2:end)))+1;
            faceCamStrobe_up_t = Timeline.rawDAQTimestamps(faceCamStrobe_up);
            
            % Get sync times for cameras (or load if already done)
            [facecam_sync_frames,n_facecam_frames] = AP_get_cam_sync_frames(facecam_fn);
            
            if ~isempty(facecam_sync_frames)
                % Get the closest cam strobe to sync start, find offset and frame idx
                [~,facecam_strobe_sync] = min(abs(camSync_flip(1) - faceCamStrobe_up));
                facecam_frame_offset = facecam_sync_frames(1) - facecam_strobe_sync;
                facecam_frame_idx = [1:length(faceCamStrobe_up)] + facecam_frame_offset;
                
                % Check that the number of frames between syncs matches
                % video and timeline
                n_facecam_frames_syncd_movie = diff(facecam_sync_frames) + 1;
                [~,facecam_strobe_sync_end] = min(abs(camSync_flip(3) - faceCamStrobe_up));
                n_facecam_frames_syncd_timeline = facecam_strobe_sync_end - facecam_strobe_sync;
                if abs(n_facecam_frames_syncd_movie - n_facecam_frames_syncd_timeline) > 1
                   error('Facecam: different n frames video vs timeline'); 
                end
                
                % Get times of cam frames in timeline
                facecam_t = nan(n_facecam_frames,1);
                facecam_t(facecam_frame_idx) = faceCamStrobe_up_t;
                
                save(facecam_t_savefile,'facecam_t');
            end
        elseif exist(facecam_fn,'file') && exist(facecam_t_savefile,'file')
            load(facecam_t_savefile);
        end
        
        % (old/unused: etGUI and facemap)
        [facecam_processed_filename,facecam_processed_exists] = AP_cortexlab_filename(animal,day,experiment,'facecam_processed');
        if facecam_processed_exists
            facecam = load(facecam_processed_filename);
        end
        
        % (output from AP_mouse_movie_movement)
        [facecam_movement_filename,facecam_movement_exists] = AP_cortexlab_filename(animal,day,experiment,'facecam_movement');
        if facecam_movement_exists
            load(facecam_movement_filename);
        end
        
    end
    
end

%% Load imaging data

[data_path,data_path_exists] = AP_cortexlab_filename(animal,day,experiment,'imaging',site);
experiment_path = [data_path filesep num2str(experiment)];

% (check for specific imaging file since data path is just root)
spatialComponents_fns = dir([data_path filesep 'svdSpatialComponents*']);
imaging_exists = ~isempty(spatialComponents_fns);

if imaging_exists && load_parts.imaging
    if verbose; disp('Loading imaging data...'); end
    
    % Get the imaging file locations
    spatialComponents_dir = dir([data_path filesep 'svdSpatialComponents*']);   
    meanImage_dir = dir([data_path filesep 'meanImage*']);
    
    cam_color_n = length(spatialComponents_dir);
    cam_color_signal = 'blue';
    cam_color_hemo = 'purple';
    
    if cam_color_n == 1
        
        U = readUfromNPY([data_path filesep spatialComponents_dir.name]);
        V = readVfromNPY([experiment_path filesep strrep(spatialComponents_dir.name,'Spatial','Temporal')]);
        frame_t = cam_time;
        
        % Get average framerate (nearest interger to prevent minor errors)
        framerate = round(1./nanmean(diff(frame_t)));
        
        % Detrend and high-pass filter
        highpassCutoff = 0.01; % Hz
        [b100s, a100s] = butter(2, highpassCutoff/(framerate/2), 'high');
        dV = detrend(V', 'linear')';
        fV = single(filter(b100s,a100s,double(dV)')');
        
        avg_im = readNPY([data_path filesep meanImage_dir.name]);
        
    elseif cam_color_n == 2
        
        % Load in all things as neural (n) or hemodynamic (h)        
        Un = readUfromNPY([data_path filesep 'svdSpatialComponents_' cam_color_signal '.npy']);
        Vn = readVfromNPY([experiment_path filesep 'svdTemporalComponents_' cam_color_signal '.npy']);
        dataSummary_n = load([data_path filesep 'dataSummary_' cam_color_signal '.mat']);
        avg_im_n = readNPY([data_path filesep 'meanImage_' cam_color_signal '.npy']);
        
        Uh = readUfromNPY([data_path filesep 'svdSpatialComponents_' cam_color_hemo '.npy']);
        Vh = readVfromNPY([experiment_path filesep 'svdTemporalComponents_' cam_color_hemo '.npy']);
        dataSummary_h = load([data_path filesep 'dataSummary_' cam_color_signal '.mat']);
        avg_im_h = readNPY([data_path filesep 'meanImage_' cam_color_hemo '.npy']);
        
        % Get frame timestamps (assume odd = blue, even = purple for now)
        tn = cam_time(1:2:end);
        th = cam_time(2:2:end);
        
        % Get average framerate (nearest interger to prevent minor errors)
        framerate = round(1./nanmean(diff(tn)));
        
        % Correct hemodynamic signal in blue from green
        % First need to shift alternating signals to be temporally aligned
        % (shifts neural to hemo)
        if verbose; disp('Correcting hemodynamics...'); end
        
        % Check if number of timeline frames matches imaged frames
        cam_tl_imaged_diff = length(cam_time) - (size(Vn,2) + size(Vh,2));
        if cam_tl_imaged_diff ~= 0
            warning(sprintf( ...
                '\n %s %s: %d timeline-imaged frames, assuming dropped at end', ...
                animal,day,cam_tl_imaged_diff));
        end
        
        % Eliminate odd frames out (unpaired colors)
        min_frames = min(size(Vn,2),size(Vh,2));
        Vn = Vn(:,1:min_frames);
        tn = tn(1:min_frames);
        
        Vh = Vh(:,1:min_frames);
        th = th(1:min_frames);
   
        % This was to get rid of bad exposures: not sure I want this though
%         cam_expose_time_reshape = ...
%             reshape(cam_expose_times(1:end-mod(length(cam_expose_times),2)),2,[]);
%         bad_cam_expose = any(cam_expose_time_reshape > ...
%             median(cam_expose_time_reshape(:))*2,1);
%         
%         if any(bad_cam_expose)
%             warning(['Bad cam expose time: ' num2str(find(bad_cam_expose)) '/' num2str(min_frames)]);
%             Vn = Vn(:,~bad_cam_expose);
%             tn = tn(~bad_cam_expose);
%             
%             Vh = Vh(:,~bad_cam_expose);
%             th = th(~bad_cam_expose);
%         end
        
        Vn_th = SubSampleShift(Vn,1,2);
        
        Vh_Un = ChangeU(Uh,Vh,Un);
        
        hemo_tform_fn = [experiment_path filesep 'hemo_tform.mat'];
        if exist(hemo_tform_fn,'file')
            % If the hemo tform matrix has been computed, load and fix
            if verbose; disp('Using old hemo tform...'); end
            load(hemo_tform_fn)
            zVh_Un = bsxfun(@minus, Vh_Un, nanmean(Vh_Un,2));
            Vn_hemo = transpose(Vn_th' - zVh_Un'*hemo_tform');
        else
            % If no p hemo tform matrix, compute and save
            if verbose; disp('Computing hemo tform...'); end
            %hemo_freq = [0.1,1];
            hemo_freq = [7,13];
            [Vn_hemo,hemo_tform] = HemoCorrectLocal(Un,Vn_th,Vh_Un,framerate,hemo_freq,3);
            save(hemo_tform_fn,'hemo_tform');
            % Close the figures (hacky - but function isn't mine)
            close(gcf)
            close(gcf)
        end
        
        if verbose; disp('Filtering...'); end
        % Don't bother filtering heartbeat, just detrend and highpass
        % fVn_hemo = detrendAndFilt(Vn_hemo, framerate);
        highpassCutoff = 0.01; % Hz
        [b100s, a100s] = butter(2, highpassCutoff/(framerate/2), 'high');
        
        dVn_hemo = detrend(Vn_hemo', 'linear')';
        
        % non-zero-lag filter, but causal (only moves forwards in time)
        fVn_hemo = filter(b100s,a100s,dVn_hemo,[],2);        
        % non-causal but zero-lag filter: changed because can introduce
        % artifacts with single wonky points, also big changes propogate
        % backwards in time which potentially gives bad causality
        %fVn_hemo = single(filtfilt(b100s,a100s,double(dVn_hemo)')');
        
        % Do this for the colors individually, in case they're used
        dVn = detrend(Vn', 'linear')';
        fVn = single(filter(b100s,a100s,double(dVn)')');
        
        dVh = detrend(Vh', 'linear')';
        fVh = single(filter(b100s,a100s,double(dVh)')');
        
        % set final U/V to use
        fV = fVn_hemo;
        U = Un;
        avg_im = avg_im_n;
        frame_t = th; % shifted to use hemo color times
        
    end
    if verbose; disp('Done.'); end
    
    % Make dF/F
    [Udf,fVdf] = dffFromSVD(U,fV,avg_im);
    % zero out NaNs in the Udfs (from saturated pixels?)
    Udf(isnan(Udf)) = 0;
end

%% Load ephys data (single long recording)

% Pick kilosort version (2 by default, 1 old if selected)
if ~exist('kilosort_version','var') || kilosort_version == 2
    [ephys_path,ephys_exists] = AP_cortexlab_filename(animal,day,experiment,'ephys',site);
elseif exist('kilosort_version','var') && kilosort_version == 1
    [ephys_path,ephys_exists] = AP_cortexlab_filename(animal,day,experiment,'ephys_ks1',site);
end

if ephys_exists && load_parts.ephys
    
    if verbose; disp('Loading ephys...'); end
    
    % These are the digital channels going into the FPGA
    photodiode_sync_idx = 1;
    acqLive_sync_idx = 2;
    pcoExposure_sync_idx = 3;
    flipper_sync_idx = 4;

    load_lfp = false;
    
    % Load phy sorting if it exists 
    % (old = cluster_groups.csv, new = cluster_group.tsv because fuck me)
    cluster_filepattern = [ephys_path filesep 'cluster_group*'];
    cluster_filedir = dir(cluster_filepattern);
    if ~isempty(cluster_filedir)
        cluster_filename = [ephys_path filesep cluster_filedir.name];
        fid = fopen(cluster_filename);
        cluster_groups = textscan(fid,'%d%s','HeaderLines',1);
        fclose(fid);
    end
    
    % Load sync/photodiode
    load(([ephys_path filesep 'sync.mat']));
    
    % Read header information
    header_path = [ephys_path filesep 'dat_params.txt'];
    header_fid = fopen(header_path);
    header_info = textscan(header_fid,'%s %s', 'delimiter',{' = '});
    fclose(header_fid);
    
    header = struct;
    for i = 1:length(header_info{1})
        header.(header_info{1}{i}) = header_info{2}{i};
    end
    
    % Load spike data
    if isfield(header,'sample_rate')
        ephys_sample_rate = str2num(header.sample_rate);
    elseif isfield(header,'ap_sample_rate')
        ephys_sample_rate = str2num(header.ap_sample_rate);
    end
    spike_times = double(readNPY([ephys_path filesep 'spike_times.npy']))./ephys_sample_rate;
    spike_templates_0idx = readNPY([ephys_path filesep 'spike_templates.npy']);
    templates_whitened = readNPY([ephys_path filesep 'templates.npy']);
    channel_positions = readNPY([ephys_path filesep 'channel_positions.npy']);
    channel_map = readNPY([ephys_path filesep 'channel_map.npy']);
    winv = readNPY([ephys_path filesep 'whitening_mat_inv.npy']);
    template_amplitudes = readNPY([ephys_path filesep 'amplitudes.npy']);

    % Default channel map/positions are from end: make from surface
    % (hardcode this: kilosort2 drops channels)
    max_depth = 3840;
    channel_positions(:,2) = max_depth - channel_positions(:,2);
    
    % Unwhiten templates
    templates = zeros(size(templates_whitened));
    for t = 1:size(templates_whitened,1)
        templates(t,:,:) = squeeze(templates_whitened(t,:,:))*winv;
    end
    
    % Get the waveform of all templates (channel with largest amplitude)
    [~,max_site] = max(max(abs(templates),[],2),[],3);
    templates_max = nan(size(templates,1),size(templates,2));
    for curr_template = 1:size(templates,1)
        templates_max(curr_template,:) = ...
            templates(curr_template,:,max_site(curr_template));
    end
    waveforms = templates_max;
    
    % Get depth of each template
    % (get min-max range for each channel)
    template_chan_amp = squeeze(range(templates,2));
    % (zero-out low amplitude channels)
    template_chan_amp_thresh = max(template_chan_amp,[],2)*0.3;
    template_chan_amp_overthresh = template_chan_amp.*(template_chan_amp >= template_chan_amp_thresh);
    % (get center-of-mass on thresholded channel amplitudes)
    template_depths = sum(template_chan_amp_overthresh.*channel_positions(:,2)',2)./sum(template_chan_amp_overthresh,2);
    
    % Get the depth of each spike (templates are zero-indexed)
    spike_depths = template_depths(spike_templates_0idx+1);
            
    % Get trough-to-peak time for each template
    templates_max_signfix = bsxfun(@times,templates_max, ...
        sign(abs(min(templates_max,[],2)) - abs(max(templates_max,[],2))));
    
    [~,waveform_trough] = min(templates_max,[],2);
    [~,waveform_peak_rel] = arrayfun(@(x) ...
        max(templates_max(x,waveform_trough(x):end),[],2), ...
        transpose(1:size(templates_max,1)));
    waveform_peak = waveform_peak_rel + waveform_trough;
    
    templateDuration = waveform_peak - waveform_trough;
    templateDuration_us = (templateDuration/ephys_sample_rate)*1e6;
    
    % Load LFP
    n_channels = str2num(header.n_channels);
    %lfp_filename = [ephys_path filesep 'lfp.dat']; (this is old)
    [data_path,data_path_exists] = AP_cortexlab_filename(animal,day,experiment,'ephys_dir',site);
    lfp_dir = dir([data_path 'experiment*-1_0.dat']);
    lfp_filename = [data_path lfp_dir.name];
    if load_lfp && exist(lfp_filename,'file')
        lfp_sample_rate = str2num(header.lfp_sample_rate);
        lfp_cutoff = str2num(header.filter_cutoff);
        
        fid = fopen(lfp_filename);
        % define where/how much of LFP to load
        lfp_skip_minutes = 10; % move to N minutes after recording start
        lfp_load_start = (lfp_sample_rate*60*lfp_skip_minutes*n_channels);
        lfp_load_samples = 1e6;
        % load LFP
        fseek(fid,lfp_load_start,'bof');
        lfp_all = fread(fid,[n_channels,lfp_load_samples],'int16'); % pull snippet
        fclose(fid);
        % eliminate non-connected channels
        lfp = lfp_all(channel_map+1,:);
        clear lfp_all;
        
        lfp_t = [(lfp_load_start/n_channels):(lfp_load_start/n_channels)+lfp_load_samples-1]/lfp_sample_rate;
    end
    
    % Get sync points for alignment
 
    % Get experiment index by finding numbered folders
    protocols = AP_list_experiments(animal,day);
    experiment_idx = experiment == [protocols.experiment];
 
    if exist('flipper_flip_times_timeline','var') && length(sync) >= flipper_sync_idx
        % (if flipper, use that)
        % (at least one experiment the acqLive connection to ephys was bad
        % so it was delayed - ideally check consistency since it's
        % redundant)
        bad_flipper = false; 
        
        % Get flipper experiment differences by long delays
        flip_diff_thresh = 1; % time between flips to define experiment gap (s)
        flipper_expt_idx = [1;find(diff(sync(flipper_sync_idx).timestamps) > ...
            flip_diff_thresh)+1;length(sync(flipper_sync_idx).timestamps)+1];
        
        flipper_flip_times_ephys = sync(flipper_sync_idx).timestamps( ...
            flipper_expt_idx(find(experiment_idx)):flipper_expt_idx(find(experiment_idx)+1)-1);
        
        % Check that number of flipper flips in timeline matches ephys
        if length(flipper_flip_times_ephys) ~= length(flipper_flip_times_timeline)
            warning([animal ' ' day ':Flipper flip times different in timeline/ephys'])
            bad_flipper = true;
        end
        
%         % Plot stim aligned by flipper times for sanity check
%         figure; hold on;
%         stim_ephys_timestamps = interp1(flipper_flip_times_ephys,flipper_flip_times_timeline,sync(1).timestamps,'linear','extrap');
%         plot(stim_ephys_timestamps,sync(1).values*5,'.g')
%         plot(stimOn_times,5,'ob');
        
        sync_timeline = flipper_flip_times_timeline;
        sync_ephys = flipper_flip_times_ephys;
    else
        bad_flipper = true;
    end
    
    if bad_flipper
        % (if no flipper or flipper problem, use acqLive)
        
        % Get acqLive times for current experiment
        experiment_ephys_starts = sync(acqLive_sync_idx).timestamps(sync(acqLive_sync_idx).values == 1);
        experiment_ephys_stops = sync(acqLive_sync_idx).timestamps(sync(acqLive_sync_idx).values == 0);
        acqlive_ephys_currexpt = [experiment_ephys_starts(experiment_idx), ...
            experiment_ephys_stops(experiment_idx)];
        
        sync_timeline = acqLive_timeline;
        sync_ephys = acqlive_ephys_currexpt;
        
        % Check that the experiment time is the same within threshold
        % (it should be almost exactly the same)
        if abs(diff(acqLive_timeline) - diff(acqlive_ephys_currexpt)) > 1
            error([animal ' ' day ': acqLive duration different in timeline and ephys']);
        end                     
    end
    
    % Get the spike/lfp times in timeline time (accounts for clock drifts)
    spike_times_timeline = interp1(sync_ephys,sync_timeline,spike_times,'linear','extrap');
    if load_lfp && exist(lfp_filename,'file')
        lfp_t_timeline = interp1(sync_ephys,sync_timeline,lfp_t,'linear','extrap');
    end
    
    % Get "good" templates from labels  
    if exist('cluster_groups','var')
        % If there's a manual classification
        if verbose; disp('Keeping manually labelled good units...'); end
        
        % Check that all used spike templates have a label
        spike_templates_0idx_unique = unique(spike_templates_0idx);
        if ~all(ismember(spike_templates_0idx_unique,uint32(cluster_groups{1}))) || ...
                ~all(ismember(cluster_groups{2},{'good','mua','noise'}))
            warning([animal ' ' day ': not all templates labeled]']);
        end
        
        % Define good units from labels
        good_templates_idx = uint32(cluster_groups{1}( ...
            strcmp(cluster_groups{2},'good') | strcmp(cluster_groups{2},'mua')));
        good_templates = ismember(0:size(templates,1)-1,good_templates_idx);
        
    elseif exist([ephys_path filesep 'cluster_AP_triage.tsv'],'file')
        % If no manual but AP_triage clusters are available
        if verbose; disp('Keeping AP_triage good units...'); end

        % Load triage labels
        triage_label_filename = [ephys_path filesep 'cluster_AP_triage.tsv'];        
        
        fid = fopen(triage_label_filename);
        triage_labels = textscan(fid,'%d%s','HeaderLines',1);
        fclose(fid);
        
        triage_good_templates = strcmp(triage_labels{2},'good');      
      
        good_templates = ...
            triage_good_templates;
        good_templates_idx = find(good_templates)-1;
        
    else
        % If no cluster groups at all, keep all
        warning([animal ' ' day ' - no cluster groups']);      
        if verbose; disp('No manual labeling, keeping all and re-indexing'); end
        good_templates_idx = unique(spike_templates_0idx);
        good_templates = ismember(0:size(templates,1)-1,good_templates_idx);
    end
    
    % Throw out all non-good template data
    templates = templates(good_templates,:,:);
    template_depths = template_depths(good_templates);
    waveforms = waveforms(good_templates,:);
    templateDuration = templateDuration(good_templates);
    templateDuration_us = templateDuration_us(good_templates);
    
    % Throw out all non-good spike data
    good_spike_idx = ismember(spike_templates_0idx,good_templates_idx);
    spike_times = spike_times(good_spike_idx);
    spike_templates_0idx = spike_templates_0idx(good_spike_idx);
    template_amplitudes = template_amplitudes(good_spike_idx);
    spike_depths = spike_depths(good_spike_idx);
    spike_times_timeline = spike_times_timeline(good_spike_idx);
    
    % Rename the spike templates according to the remaining templates
    % (and make 1-indexed from 0-indexed)
    new_spike_idx = nan(max(spike_templates_0idx)+1,1);
    new_spike_idx(good_templates_idx+1) = 1:length(good_templates_idx);
    spike_templates = new_spike_idx(spike_templates_0idx+1);   
    
end

%% Regress light artifact out of LFP (unused at the moment)

% if ephys_exists && load_parts.ephys
%     
%     if verbose; disp('Cleaning LFP...'); end;
%     
%     % LFP correlation
%     % (fix the light artifact on each channel with regression)
%     light_timeline = interp1(acqlive_ephys_currexpt,acqLive_timeline,sync(3).timestamps,'linear','extrap');
%     light_on = light_timeline(sync(3).values == 1);
%     light_off = light_timeline(sync(3).values == 0);
%     
%     blue_on = light_on(1:2:end);
%     blue_off = light_off(1:2:end);
%     violet_on = light_on(2:2:end);
%     violet_off = light_off(2:2:end);
%     
%     lfp_t_bins = [lfp_t_timeline-0.5/lfp_sample_rate,lfp_t_timeline(end)+0.5/lfp_sample_rate];
%     blue_on_vector = histcounts(blue_on,lfp_t_bins);
%     blue_off_vector = histcounts(blue_off,lfp_t_bins);
%     violet_on_vector = histcounts(violet_on,lfp_t_bins);
%     violet_off_vector = histcounts(violet_off,lfp_t_bins);
%     
%     light_vectors = [blue_on_vector;blue_off_vector;violet_on_vector;violet_off_vector];
%     
%     t_shift = round((1/35)*lfp_sample_rate*1.5);
%     t_shifts = [-t_shift:t_shift];
%     lambda = 0;
%     zs = [false,false];
%     cvfold = 1;
%     
%     % (in chunks: necessary memory-wise, also allows changing light)
%     n_chunks = 10;
%     lfp_t_chunk = round(linspace(1,size(lfp,2),n_chunks+1));
%     
%     lfp_lightfix = nan(size(lfp));
%     for curr_chunk = 1:n_chunks
%         curr_chunk_t = lfp_t_chunk(curr_chunk):lfp_t_chunk(curr_chunk+1);
%         [light_k,artifact_lfp] = AP_regresskernel(light_vectors(:,curr_chunk_t),lfp(:,curr_chunk_t),t_shifts,lambda,zs,cvfold);
%         
%         lfp_lightfix(:,curr_chunk_t) = lfp(:,curr_chunk_t)-artifact_lfp;
%         %             AP_print_progress_fraction(curr_chunk,n_chunks);
%     end
%     lfp_lightfix(isnan(lfp_lightfix)) = 0;
%     
%     % (group channels by depth)
%     channel_depth_grp = discretize(channel_positions(:,2),depth_group_edges);
%     lfp_depth_median = grpstats(lfp_lightfix,channel_depth_grp,'median');
%     
%     % (low-pass filter: sometimes bunch of junk at high freq?)
%     freqCutoff = 300; % Hz
%     [b100s, a100s] = butter(2,freqCutoff/(lfp_sample_rate/2),'low');
%     lfp_depth_median_filt = single(filtfilt(b100s,a100s,double(lfp_depth_median)')');
%     
% end

%% Estimate striatal boundaries on probe

if ephys_exists && load_parts.ephys
    if verbose; disp('Estimating striatum boundaries on probe...'); end
    
    % str_align = alignment method ('none', 'depth', or 'kernel')
    
    % requires n_aligned_depths for alignment, default 4
    if ~exist('n_aligned_depths','var')
        n_aligned_depths = 4;
    end
    
    % if no alignment specified, default kernel
    if ~exist('str_align','var')
        str_align = 'kernel';
    end
    
    AP_align_striatum_ephys
    
    if verbose
       figure;
       imagesc(depth_group_centers,depth_group_centers,mua_corr);
       axis tight equal;
       colormap(hot)
       line([str_depth(1),str_depth(1)],ylim,'color','b','linewidth',2);
       line([str_depth(2),str_depth(2)],ylim,'color','b','linewidth',2);
       line(xlim,[str_depth(1),str_depth(1)],'color','b','linewidth',2);
       line(xlim,[str_depth(2),str_depth(2)],'color','b','linewidth',2);
       xlabel('Probe depth (\mum)');
       ylabel('Probe depth (\mum)');
       title('MUA correlation: striatum location');
    end
    
end

%% Classify spikes

if ephys_exists && load_parts.ephys
    if verbose; disp('Classifying spikes...'); end
    
    str_templates = template_depths >= str_depth(1) & template_depths <= str_depth(2);
    non_str_templates = ~str_templates;
    
    % Define the window to look for spiking statistics in (spikes go in and
    % out, so take the bin with the largest firing rate for each cell and work
    % with that one)
    % spiking_stat_window = 60*5; % seconds
    % spiking_stat_bins = min(spike_times_timeline):spiking_stat_window: ...
    %     max(spike_times_timeline);
    
    % % (for whole session)
    spiking_stat_window = max(spike_times_timeline)-min(spike_times_timeline);
    spiking_stat_bins = [min(spike_times_timeline),max(spike_times_timeline)];
    
    % Get firing rate across the session
    bin_spikes = nan(size(templates,1), ...
        length(spiking_stat_bins)-1);
    for curr_template = unique(spike_templates)'
        bin_spikes(curr_template,:) = ...
            histcounts(spike_times_timeline(spike_templates == curr_template), ...
            spiking_stat_bins);
    end
    min_spikes = 10;
    use_spiking_stat_bins = bsxfun(@ge,bin_spikes,prctile(bin_spikes,80,2)) & bin_spikes > min_spikes;
    spike_rate = sum(bin_spikes.*use_spiking_stat_bins,2)./ ...
        (sum(use_spiking_stat_bins,2)*spiking_stat_window);
    
    % Get proportion of ISI > 2s (Yamin/Cohen 2013) and CV2 (Stalnaker/Schoenbaum 2016)
    prop_long_isi = nan(size(templates,1),1);
    cv2 = nan(size(templates,1),1);
    for curr_template = unique(spike_templates)'
        
        long_isi_total = 0;
        isi_ratios = [];
        for curr_bin = find(use_spiking_stat_bins(curr_template,:))
            curr_spike_times = spike_times_timeline( ...
                spike_times_timeline > spiking_stat_bins(curr_bin) & ...
                spike_times_timeline < spiking_stat_bins(curr_bin+1) & ...
                spike_templates == curr_template);
            curr_isi = diff(curr_spike_times);
            
            long_isi_total = long_isi_total + sum(curr_isi(curr_isi > 2));
            
            isi_ratios = [isi_ratios;(2*abs(curr_isi(2:end) - curr_isi(1:end-1)))./ ...
                (curr_isi(2:end) + curr_isi(1:end-1))];
        end
        
        prop_long_isi(curr_template) = long_isi_total/ ...
            (sum(use_spiking_stat_bins(curr_template,:))*spiking_stat_window);
        cv2(curr_template) = nanmean(isi_ratios);
        
    end
    
    % Cortical classification (like Bartho JNeurophys 2004)
    waveform_duration_cutoff = 400;
    narrow = non_str_templates & templateDuration_us <= waveform_duration_cutoff;
    wide = non_str_templates & templateDuration_us > waveform_duration_cutoff;
    
    % Striatum classification
    prop_long_isi_cutoff = 0.35;
    cv2_cutoff = 0.8;
    
    msn = str_templates & ...
        templateDuration_us > waveform_duration_cutoff & ...
        prop_long_isi >= prop_long_isi_cutoff;
    
    fsi = str_templates & ...
        templateDuration_us <= waveform_duration_cutoff & ...
        prop_long_isi < prop_long_isi_cutoff;
    
    tan = str_templates & ...
        templateDuration_us > waveform_duration_cutoff & ...
        prop_long_isi < prop_long_isi_cutoff;
    
    uin = str_templates & ~msn & ~fsi & ~tan;
    
    waveform_t = 1e3*((0:size(templates,2)-1)/ephys_sample_rate);
    
    if verbose
        
        % Plot the waveforms and spike statistics
        figure;
        
        if any(non_str_templates)
            subplot(2,2,1); hold on;
            p = plot(waveform_t,waveforms(non_str_templates,:)');
            set(p(wide(non_str_templates)),'color','k')
            set(p(narrow(non_str_templates)),'color','r')
            xlabel('Time (ms)')
            title('Not striatum');
            legend([p(find(wide(non_str_templates),1)),p(find(narrow(non_str_templates),1))],{'Wide','Narrow'})
        end
        
        subplot(2,2,2); hold on;
        p = plot(waveform_t,waveforms(str_templates,:)');
        set(p(msn(str_templates)),'color','m')
        set(p(fsi(str_templates)),'color','b')
        set(p(tan(str_templates)),'color','g')
        set(p(uin(str_templates)),'color','c')
        xlabel('Time (ms)')
        title('Striatum');
        legend([p(find(msn(str_templates),1)),p(find(fsi(str_templates),1)), ...
            p(find(tan(str_templates),1)),p(find(uin(str_templates),1))],{'MSN','FSI','TAN','UIN'});
        
        subplot(2,2,3); hold on;
        
        stem3( ...
            templateDuration_us(wide)/1000, ...
            prop_long_isi(wide), ...
            spike_rate(wide),'k');
        
        stem3( ...
            templateDuration_us(narrow)/1000, ...
            prop_long_isi(narrow), ...
            spike_rate(narrow),'r');
        
        xlabel('waveform duration (ms)')
        ylabel('frac long ISI')
        zlabel('spike rate')
        
        set(gca,'YDir','reverse')
        set(gca,'XDir','reverse')
        view(3);
        grid on;
        axis vis3d;
        
        subplot(2,2,4); hold on;
        stem3( ...
            templateDuration_us(msn)/1000, ...
            prop_long_isi(msn), ...
            spike_rate(msn),'m');
        
        stem3( ...
            templateDuration_us(fsi)/1000, ...
            prop_long_isi(fsi), ...
            spike_rate(fsi),'b');
        
        stem3( ...
            templateDuration_us(tan)/1000, ...
            prop_long_isi(tan), ...
            spike_rate(tan),'g');
        
        stem3( ...
            templateDuration_us(uin)/1000, ...
            prop_long_isi(uin), ...
            spike_rate(uin),'c');
        
        xlabel('waveform duration (ms)')
        ylabel('frac long ISI')
        zlabel('spike rate')
        
        set(gca,'YDir','reverse')
        set(gca,'XDir','reverse')
        view(3);
        grid on;
        axis vis3d;
        
        % Plot cell type by depth
        celltype_labels = {'Wide','Narrow','MSN','FSI','TAN','UIN'};
        celltypes = wide.*1 + narrow.*2 + msn.*3 + fsi.*4 + tan.*5 + uin.*6;
        use_colors = ...
            [0,0,0;
            1,0,0;
            1,0,1;
            0,0,1;
            0,1,0;
            0,1,1];
        
        plot_celltypes = any([wide,narrow,msn,fsi,tan,uin],1);
        
        figure('Position',[94,122,230,820]);
        scatter(rand(size(template_depths))-0.5,template_depths,10,use_colors(celltypes,:),'filled');
        xlim([-1,1])
        set(gca,'XTick',[]);
        set(gca,'YDir','reverse');
        ylabel('Depth (\mum)');
        legend(celltype_labels(plot_celltypes));
        ylim([0,max(channel_positions(:,2))])
        
        drawnow;
        
    end
    
end


%% Finished
if verbose; disp('Finished loading experiment.'); end








