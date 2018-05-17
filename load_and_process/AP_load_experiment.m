% AP_load_experiment(animal,day,experiment,site)
%
% Loads data from experiments
% assumes kilotrode, among other things
%
% Not a function at the moment because nothing is packaged

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
    if ~isfield(load_parts,'cam');
        load_parts.cam = false;
    end
    if ~isfield(load_parts,'imaging');
        load_parts.imaging = false;
    end
    if ~isfield(load_parts,'ephys');
        load_parts.ephys = false;
    end
end

%% Load timeline

[timeline_filename,timeline_exists] = AP_cortexlab_filename(animal,day,experiment,'timeline');

if timeline_exists
    if verbose; disp('Loading timeline...'); end;
    
    load(timeline_filename);
    
    % Set rig-specific timeline names
    cam_name = 'pcoExposure';
    acqLive_name = 'acqLive';
    
    % Get camera times
    timeline_cam_idx = strcmp({Timeline.hw.inputs.name}, cam_name);
    cam_samples = find(Timeline.rawDAQData(1:end-1,timeline_cam_idx) <= 2 & ...
        Timeline.rawDAQData(2:end,timeline_cam_idx) > 2) + 1;
    cam_time = Timeline.rawDAQTimestamps(cam_samples);
    
    % Get acqLive signal
    acqLive_idx = strcmp({Timeline.hw.inputs.name}, acqLive_name);
    thresh = max(Timeline.rawDAQData(:,acqLive_idx))/2;
    acqLive_trace = Timeline.rawDAQData(:,acqLive_idx) > thresh;
    acqLive_timeline = Timeline.rawDAQTimestamps( ...
        [find(acqLive_trace,1),find(acqLive_trace,1,'last')+1]);
end

%% Load mpep protocol

[protocol_filename,protocol_exists] = AP_cortexlab_filename(animal,day,experiment,'protocol');

if protocol_exists
    
    if verbose; disp('Loading mpep protocol...'); end;
    
    load(protocol_filename);
    
    % Load in hardware info
    hwinfo_filename = AP_cortexlab_filename(animal,day,experiment,'hardware');
    load(hwinfo_filename);
    
    % Get flicker or steady photodiode
    photodiode_type = myScreenInfo.SyncSquare.Type;
    
    % Get stimulus onsets and parameters
    photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
    
    % Get stim screen signal index
    stimScreen_idx = strcmp({Timeline.hw.inputs.name}, 'stimScreen');
    if any(stimScreen_idx)
        stimScreen_flicker = max(Timeline.rawDAQData(:,stimScreen_idx)) - ...
            min(Timeline.rawDAQData(:,stimScreen_idx)) > 2;
        stimScreen_thresh = max(Timeline.rawDAQData(:,stimScreen_idx))/2;
        stimScreen_on = Timeline.rawDAQData(:,stimScreen_idx) > stimScreen_thresh;
    end
    
    switch lower(photodiode_type)
        case 'flicker'
            warning('if flickering photodiode and steady screen, write diff')
            %         % get differential of photodiode
            %         photodiode_diff = [diff(Timeline.rawDAQData(:,photodiode_idx));0];
            %
            %         stimScreen_idx = strcmp({Timeline.hw.inputs.name}, 'stimScreen');
            %         if any(Timeline.rawDAQData(:,stimScreen_idx) < 1);
            %             stimScreen_thresh = max(Timeline.rawDAQData(:,stimScreen_idx))/2;
            %             stimScreen_on = Timeline.rawDAQData(:,stimScreen_idx) > stimScreen_thresh;
            %             photodiode_diff(~stimScreen_on) = NaN;
            %             photodiode_
            %
            %         end
            
            % This is the same as below... can probably just use
            stimScreen_on = Timeline.rawDAQData(:,photodiode_idx) > 0.15;
            stimScreen_on_t = Timeline.rawDAQTimestamps(stimScreen_on);
            photodiode_thresh = max(Timeline.rawDAQData(:,photodiode_idx))/2;
            % median filter because of weird effect where
            % photodiode dims instead of off for one sample
            % while backlight is turning off
            photodiode_trace = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
                photodiode_idx),5) > photodiode_thresh;
            photodiode_flip = find((~photodiode_trace(1:end-1) & photodiode_trace(2:end)) | ...
                (photodiode_trace(1:end-1) & ~photodiode_trace(2:end)))+1;
            
            photodiode = struct('timestamps',[],'values',[]);
            photodiode.timestamps = stimScreen_on_t(photodiode_flip)';
            photodiode.values = photodiode_trace(photodiode_flip);
            
        case 'steady'
            
            % Take into account if the screen flickers
            
            % have to redefine periods of screen on, because
            % sometimes there's a sample or so difference
            stimScreen_on = Timeline.rawDAQData(:,photodiode_idx) > 0.3;
            stimScreen_on_t = Timeline.rawDAQTimestamps(stimScreen_on);
            photodiode_thresh = (max(Timeline.rawDAQData(:,photodiode_idx)) ...
                - min(Timeline.rawDAQData(:,photodiode_idx)))/2 + ...
                min(Timeline.rawDAQData(:,photodiode_idx));
            % median filter because of weird effect where
            % photodiode dims instead of off for one sample
            % while backlight is turning off
            photodiode_trace = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
                photodiode_idx),10) > photodiode_thresh;
            photodiode_flip = find((~photodiode_trace(1:end-1) & photodiode_trace(2:end)) | ...
                (photodiode_trace(1:end-1) & ~photodiode_trace(2:end)))+1;
            
            photodiode = struct('timestamps',[],'values',[]);
            photodiode.timestamps = stimScreen_on_t(photodiode_flip)';
            photodiode.values = photodiode_trace(photodiode_flip);
    end
    
    photodiode_offsets = photodiode.timestamps(photodiode.values == 0);
    photodiode_onsets = photodiode.timestamps(photodiode.values == 1);
    
    % Get specific stim onsets by time between last offset and new onset
    % (occasionally there a bad frame so flip but not new stim)
    refresh_rate_cutoff = 1/5;
    stimOn_times = photodiode_onsets( ...
        [1;find(photodiode_onsets(2:end) - photodiode_offsets(1:end-1) > refresh_rate_cutoff) + 1]);
    
    if length(stimOn_times) ~= numel(Protocol.seqnums)
        error('MPEP/Photodiode error: photodiode doesn''t match stim')
    end
    
    stimIDs = zeros(size(stimOn_times));
    for q = 1:size(Protocol.seqnums,1)
        stimIDs(Protocol.seqnums(q,:)) = q;
    end
    
end

%% Load task/behavior

% Load the block
[block_filename, block_exists] = AP_cortexlab_filename(animal,day,experiment,'block');

if block_exists
    
    if verbose; disp('Loading block file...'); end;
    
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
            if isempty(signals_events.(block_fieldnames{curr_times}));
                % skip if empty
                continue
            end
            signals_events.(block_fieldnames{curr_times}) = ...
                AP_clock_fix(block.events.(block_fieldnames{curr_times}),reward_t_block,reward_t_timeline);
        end
    end
    
    % Get photodiode flips
    % (get stim screen flickering, if that happens)
    stimScreen_idx = strcmp({Timeline.hw.inputs.name}, 'stimScreen');
    if any(stimScreen_idx)
        stimScreen_flicker = max(Timeline.rawDAQData(:,stimScreen_idx)) - ...
            min(Timeline.rawDAQData(:,stimScreen_idx)) > 2;
        stimScreen_thresh = max(Timeline.rawDAQData(:,stimScreen_idx))/2;
        stimScreen_on = Timeline.rawDAQData(:,stimScreen_idx) > stimScreen_thresh;
        stimScreen_on_t = Timeline.rawDAQTimestamps(stimScreen_on);
    end
    % median filter because of weird effect where
    % photodiode dims instead of off for one sample
    % while backlight is turning off
    photodiode_name = 'photoDiode';
    photodiode_idx = strcmp({Timeline.hw.inputs.name}, photodiode_name);
    photodiode_trace = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
        photodiode_idx),10) > 2;
    photodiode_flip = find((~photodiode_trace(1:end-1) & photodiode_trace(2:end)) | ...
        (photodiode_trace(1:end-1) & ~photodiode_trace(2:end)))+1;
    photodiode_flip_times = stimScreen_on_t(photodiode_flip)';
    
    % SPECIFIC TO PROTOCOL
    [~,expDef] = fileparts(block.expDef);
    if strcmp(expDef,'vanillaChoiceworld');
        
        % dumb signals thing, fix
        signals_events.hitValues = circshift(signals_events.hitValues,[0,-1]);
        signals_events.missValues = circshift(signals_events.missValues,[0,-1]);
        
        % Get stim on times by closest photodiode flip
        [~,closest_stimOn_photodiode] = ...
            arrayfun(@(x) min(abs(signals_events.stimOnTimes(x) - ...
            photodiode_flip_times)), ...
            1:length(signals_events.stimOnTimes));
        stimOn_times = photodiode_flip_times(closest_stimOn_photodiode);
        
        % Get time from stim on to first wheel movement
        rotaryEncoder_idx = strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder');
        % (this is a very strange hack to overcome a problem in the rotary
        % encoder that's known in the lab and was put on the wiki)
        wheel_position = Timeline.rawDAQData(:,rotaryEncoder_idx);
        wheel_position(wheel_position > 2^31) = wheel_position(wheel_position > 2^31) - 2^32;
        
        surround_time = [-0.5,2];
        surround_samples = surround_time/Timeline.hw.samplingInterval;
        
        % (wheel velocity by smoothing the wheel trace and taking dt/t)
        wheel_smooth_t = 0.05; % seconds
        wheel_smooth_samples = wheel_smooth_t/Timeline.hw.samplingInterval;
        wheel_velocity = diff(smooth(wheel_position,wheel_smooth_samples));
        
        surround_time = surround_time(1):Timeline.hw.samplingInterval:surround_time(2);
        pull_times = bsxfun(@plus,stimOn_times,surround_time);
        
        stim_aligned_wheel_raw = interp1(Timeline.rawDAQTimestamps, ...
            wheel_position,pull_times);
        stim_aligned_wheel = bsxfun(@minus,stim_aligned_wheel_raw, ...
            nanmedian(stim_aligned_wheel_raw(:,surround_time < 0),2));
        
        thresh_displacement = 2;
        [~,wheel_move_sample] = max(abs(stim_aligned_wheel) > thresh_displacement,[],2);
        wheel_move_time = arrayfun(@(x) pull_times(x,wheel_move_sample(x)),1:size(pull_times,1));
        wheel_move_time(wheel_move_sample == 1) = NaN;
        
        % Get conditions for all trials
        % (trial timing)
        n_trials = length(block.paramsValues);
        trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
        stim_to_move = padarray(wheel_move_time - stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
        stim_to_feedback = padarray(signals_events.responseTimes, ...
            [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
            padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
        
        % (early vs late move)
        trial_timing = 1 + (stim_to_move > 0.5);
        
        % (left vs right choice)
        go_left = (signals_events.trialSideValues == 1 & signals_events.hitValues == 1) | ...
            (signals_events.trialSideValues == -1 & signals_events.missValues == 1);
        go_right = (signals_events.trialSideValues == -1 & signals_events.hitValues == 1) | ...
            (signals_events.trialSideValues == 1 & signals_events.missValues == 1);
        trial_choice = go_right - go_left;
        
        % (trial conditions: [contrast,side,choice,timing])
        contrasts = [0,0.06,0.125,0.25,0.5,1];
        sides = [-1,1];
        choices = [-1,1];
        timings = [1,2];
        
        conditions = combvec(contrasts,sides,choices,timings)';
        n_conditions = size(conditions,1);
        
        trial_conditions = ...
            [signals_events.trialContrastValues(1:n_trials); signals_events.trialSideValues(1:n_trials); ...
            trial_choice(1:n_trials); trial_timing(1:n_trials)]';
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
        % permanent unless it works fine: get stim times by photodiode
        % flips that are separated by 1s
        photodiode_flip_diff = diff(stimScreen_on_t(photodiode_flip));
        stimOn_idx = find(photodiode_flip_diff > 0.9 & photodiode_flip_diff < 1.1);
        
        stimOn_times = stimScreen_on_t(photodiode_flip(stimOn_idx));
        
        % assume the times correspond to the last n values (this is because
        % sometimes if the buffer time wasn't enough, the first stimuli
        % weren't shown or weren't shown completely)
        [conditions,conditions_idx,stimIDs] = unique(signals_events.visualParamsValues(:, ...
            signals_events.visualOnsetValues(end-length(stimOn_times)+1:end))','rows');
        conditions_params = signals_events.visualParamsValues(:,conditions_idx);
        
    elseif strcmp(expDef,'DS_choiceWorldStimPassive')
        % photodiode turns gray? change threshold
        photodiode_idx = strcmp({Timeline.hw.inputs.name}, photodiode_name);
        photodiode_trace = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
            photodiode_idx),10) > 6;
        photodiode_flip = find((~photodiode_trace(1:end-1) & photodiode_trace(2:end)) | ...
            (photodiode_trace(1:end-1) & ~photodiode_trace(2:end)))+1;
        photodiode_flip_times = stimScreen_on_t(photodiode_flip)';
        
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
        % photodiode turns gray? change threshold
        photodiode_idx = strcmp({Timeline.hw.inputs.name}, photodiode_name);
        photodiode_trace = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
            photodiode_idx),10) > 6;
        photodiode_flip = find((~photodiode_trace(1:end-1) & photodiode_trace(2:end)) | ...
            (photodiode_trace(1:end-1) & ~photodiode_trace(2:end)))+1;
        photodiode_flip_times = stimScreen_on_t(photodiode_flip)';
        
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
        
    else
        error('Signals protocol with no analysis script')
    end
    
    
    
end


%% Load face/eyecam processing (with eyeGUI)

% Don't load if no timeline
if exist('Timeline','var') && load_parts.cam
    
    % Get cam sync from timeline
    camSync_idx = strcmp({Timeline.hw.inputs.name}, 'camSync');
    camSync_thresh = max(Timeline.rawDAQData(:,camSync_idx))/2;
    camSync = Timeline.rawDAQData(:,camSync_idx) > camSync_thresh;
    camSync_up = find((~camSync(1:end-1) & camSync(2:end)))+1;
    
    % EYECAM
    [eyecam_dir,eyecam_exists] = AP_cortexlab_filename(animal,day,experiment,'eyecam');
    
    if eyecam_exists
        if verbose; disp('Loading eyecam...'); end;
        
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
            eyeCamStrobe_idx = strcmp({Timeline.hw.inputs.name}, 'eyeCameraStrobe');
            eyeCamStrobe_thresh = max(Timeline.rawDAQData(:,eyeCamStrobe_idx))/2;
            eyeCamStrobe = Timeline.rawDAQData(:,eyeCamStrobe_idx) > eyeCamStrobe_thresh;
            eyeCamStrobe_up = find((~eyeCamStrobe(1:end-1) & eyeCamStrobe(2:end)))+1;
            eyeCamStrobe_up_t = Timeline.rawDAQTimestamps(eyeCamStrobe_up);
            
            % Get sync times for cameras (or load if already done)
            [eyecam_sync_frames,n_eyecam_frames] = AP_get_cam_sync_frames(eyecam_fn);
            
            if ~isempty(eyecam_sync_frames)
                % Get the closest facecam strobe to sync start, find offset and frame idx
                [~,eyecam_strobe_sync] = min(abs(camSync_up(1) - eyeCamStrobe_up));
                eyecam_frame_offset = eyecam_sync_frames(1) - eyecam_strobe_sync;
                eyecam_frame_idx = [1:length(eyeCamStrobe_up)] + eyecam_frame_offset;
                
                % Get times of facecam frames in timeline
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
        if verbose; disp('Loading facecam...'); end;
        
        [facecam_processed_filename,facecam_processed_exists] = AP_cortexlab_filename(animal,day,experiment,'facecam_processed');
        if facecam_processed_exists
            facecam = load(facecam_processed_filename);
        end
        
        % Get camera times
        facecam_fn = AP_cortexlab_filename(animal,day,experiment,'facecam');
        facecam_dir = fileparts(facecam_fn);
        facecam_t_savefile = [facecam_dir filesep 'facecam_t.mat'];
        
        if exist(facecam_fn,'file') && ~exist(facecam_t_savefile,'file')
            % Get facecam strobes
            faceCamStrobe_idx = strcmp({Timeline.hw.inputs.name}, 'faceCamStrobe');
            faceCamStrobe_thresh = max(Timeline.rawDAQData(:,faceCamStrobe_idx))/2;
            faceCamStrobe = Timeline.rawDAQData(:,faceCamStrobe_idx) > faceCamStrobe_thresh;
            faceCamStrobe_up = find((~faceCamStrobe(1:end-1) & faceCamStrobe(2:end)))+1;
            faceCamStrobe_up_t = Timeline.rawDAQTimestamps(faceCamStrobe_up);
            
            % Get sync times for cameras (or load if already done)
            [facecam_sync_frames,n_facecam_frames] = AP_get_cam_sync_frames(facecam_fn);
            
            if ~isempty(facecam_sync_frames)
                % Get the closest facecam strobe to sync start, find offset and frame idx
                [~,facecam_strobe_sync] = min(abs(camSync_up(1) - faceCamStrobe_up));
                facecam_frame_offset = facecam_sync_frames(1) - facecam_strobe_sync;
                facecam_frame_idx = [1:length(faceCamStrobe_up)] + facecam_frame_offset;
                
                % Get times of facecam frames in timeline
                facecam_t = nan(n_facecam_frames,1);
                facecam_t(facecam_frame_idx) = faceCamStrobe_up_t;
                
                save(facecam_t_savefile,'facecam_t');
            end
        elseif exist(facecam_fn,'file') && exist(facecam_t_savefile,'file')
            load(facecam_t_savefile);
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
    if verbose; disp('Loading imaging data...'); end;
    
    % Get the imaging file locations
    spatialComponents_dir = dir([data_path filesep 'svdSpatialComponents*']);
    temporalComponents_dir = dir([experiment_path filesep 'svdTemporalComponents*']);
    imaging_timestamps_idx = cellfun(@any,strfind({temporalComponents_dir.name},'timestamps'));
    
    meanImage_dir = dir([data_path filesep 'meanImage*']);
    
    cam_color_n = length(spatialComponents_dir);
    cam_color_signal = 'blue';
    cam_color_hemo = 'purple';
    
    if cam_color_n == 1
        
        U = readUfromNPY([data_path filesep spatialComponents_dir.name]);
        V = readVfromNPY([experiment_path filesep temporalComponents_dir(~imaging_timestamps_idx).name]);
        frame_t = readNPY([experiment_path filesep temporalComponents_dir(imaging_timestamps_idx).name]);
        
        framerate = 1./nanmedian(diff(frame_t));
        
        % Detrend and high-pass filter
        highpassCutoff = 0.01; % Hz
        [b100s, a100s] = butter(2, highpassCutoff/(framerate/2), 'high');
        dV = detrend(V', 'linear')';
        fV = single(filtfilt(b100s,a100s,double(dV)')');
        
        avg_im = readNPY([data_path filesep meanImage_dir.name]);
        
    elseif cam_color_n == 2
        
        % Load in all things as neural (n) or hemodynamic (h)
        
        tn = readNPY([experiment_path filesep 'svdTemporalComponents_' cam_color_signal '.timestamps.npy']);
        Un = readUfromNPY([data_path filesep 'svdSpatialComponents_' cam_color_signal '.npy']);
        Vn = readVfromNPY([experiment_path filesep 'svdTemporalComponents_' cam_color_signal '.npy']);
        dataSummary_n = load([data_path filesep 'dataSummary_' cam_color_signal '.mat']);
        avg_im_n = readNPY([data_path filesep 'meanImage_' cam_color_signal '.npy']);
        
        th = readNPY([experiment_path filesep 'svdTemporalComponents_' cam_color_hemo '.timestamps.npy']);
        Uh = readUfromNPY([data_path filesep 'svdSpatialComponents_' cam_color_hemo '.npy']);
        Vh = readVfromNPY([experiment_path filesep 'svdTemporalComponents_' cam_color_hemo '.npy']);
        dataSummary_h = load([data_path filesep 'dataSummary_' cam_color_signal '.mat']);
        avg_im_h = readNPY([data_path filesep 'meanImage_' cam_color_hemo '.npy']);
        
        framerate = 1./nanmedian(diff(tn));
        
        % Correct hemodynamic signal in blue from green
        % First need to shift alternating signals to be temporally aligned
        % (shifts neural to hemo)
        % Eliminate odd frames out
        if verbose; disp('Correcting hemodynamics...'); end
        
        min_frames = min(size(Vn,2),size(Vh,2));
        Vn = Vn(:,1:min_frames);
        Vh = Vh(:,1:min_frames);
        
        Vn_th = SubSampleShift(Vn,1,2);
        
        Vh_Un = ChangeU(Uh,Vh,Un);
        
        hemo_tform_fn = [experiment_path filesep 'hemo_tform.mat'];
        if exist(hemo_tform_fn,'file')
            % If the hemo tform matrix has been computed, load and fix
            if verbose; disp('Using old hemo tform...'); end;
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
        
        if verbose; disp('Filtering...'); end;
        % Don't bother filtering heartbeat, just detrend and highpass
        % fVn_hemo = detrendAndFilt(Vn_hemo, framerate);
        highpassCutoff = 0.01; % Hz
        [b100s, a100s] = butter(2, highpassCutoff/(framerate/2), 'high');
        
        dVn_hemo = detrend(Vn_hemo', 'linear')';
        % wasn't zero-lag filtered before? why not?
        %fVn_hemo = filter(b100s,a100s,dVn_hemo,[],2);
        fVn_hemo = single(filtfilt(b100s,a100s,double(dVn_hemo)')');
        
        % Do this for the colors individually, in case they're used
        dVn = detrend(Vn', 'linear')';
        fVn = single(filtfilt(b100s,a100s,double(dVn)')');
        
        dVh = detrend(Vh', 'linear')';
        fVh = single(filtfilt(b100s,a100s,double(dVh)')');
        
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
[ephys_path,ephys_exists] = AP_cortexlab_filename(animal,day,experiment,'ephys',site);

if ephys_exists && load_parts.ephys
    
    if verbose; disp('Loading ephys...'); end;
    
    acqLive_channel = 2;
    load_lfp = true;
    
    % Load clusters, if they exist
    cluster_filename = [ephys_path filesep 'cluster_groups.csv'];
    if exist(cluster_filename,'file')
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
    spike_templates = readNPY([ephys_path filesep 'spike_templates.npy']);
    templates = readNPY([ephys_path filesep 'templates.npy']);
    channel_positions = readNPY([ephys_path filesep 'channel_positions.npy']);
    channel_map = readNPY([ephys_path filesep 'channel_map.npy']);
    winv = readNPY([ephys_path filesep 'whitening_mat_inv.npy']);
    template_amplitudes = readNPY([ephys_path filesep 'amplitudes.npy']);
    
    % Flip channel map and positions if banks are reversed
    % (this was only for phase 2, so setting false by default)
    flipped_banks = false;
    if flipped_banks
        channel_map = [channel_map(61:end);channel_map(1:60)];
        channel_positions = [channel_positions(61:end,:);channel_positions(1:60,:)];
    end
    
    % Default channel map/positions are from end: make from surface
    channel_positions(:,2) = max(channel_positions(:,2)) - channel_positions(:,2);
    
    % Load LFP
    n_channels = str2num(header.n_channels);
    %lfp_filename = [ephys_path filesep 'lfp.dat']; (this is old)
    [data_path,data_path_exists] = AP_cortexlab_filename(animal,day,experiment,'ephysraw',site);
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
    
    % Get acqLive times for current experiment
    experiment_ephys_starts = sync(acqLive_channel).timestamps(sync(acqLive_channel).values == 1);
    experiment_ephys_stops = sync(acqLive_channel).timestamps(sync(acqLive_channel).values == 0);
    
    % (get folders with only a number - those're the experiment folders)
    experiments_dir = dir(AP_cortexlab_filename(animal,day,experiment,'expInfo'));
    experiments_num_idx = cellfun(@(x) ~isempty(x), regexp({experiments_dir.name},'^\d*$'));
    experiment_num = experiment == cellfun(@str2num,{experiments_dir(experiments_num_idx).name});
    acqlive_ephys_currexpt = [experiment_ephys_starts(experiment_num), ...
        experiment_ephys_stops(experiment_num)];
    
    % Get the spike/lfp times in timeline time (accounts for clock drifts)
    spike_times_timeline = AP_clock_fix(spike_times,acqlive_ephys_currexpt,acqLive_timeline);
    if load_lfp && exist(lfp_filename,'file')
        lfp_t_timeline = AP_clock_fix(lfp_t,acqlive_ephys_currexpt,acqLive_timeline);
    end
    
    % Get the depths of each template
    % (by COM - this used to not work but now looks ok)
    [spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = ...
        templatePositionsAmplitudes(templates,winv,channel_positions(:,2),spike_templates,template_amplitudes);
    %     % (by max waveform channel)
    %     template_abs = permute(max(abs(templates),[],2),[3,1,2]);
    %     [~,max_channel_idx] =  max(template_abs,[],1);
    %     templateDepths = channel_positions(max_channel_idx,2);
    %     % Get each spike's depth
    %     spikeDepths = templateDepths(spike_templates+1);
    
    % Get the waveform duration of all templates (channel with largest amp)
    [~,max_site] = max(max(abs(templates),[],2),[],3);
    templates_max = nan(size(templates,1),size(templates,2));
    for curr_template = 1:size(templates,1)
        templates_max(curr_template,:) = ...
            templates(curr_template,:,max_site(curr_template));
    end
    waveforms = templates_max;
    
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
    
    % Eliminate spikes that were classified as not "good"
    if exist('cluster_groups','var')
        
        if verbose; disp('Removing non-good templates'); end;
        
        good_templates_idx = uint32(cluster_groups{1}(strcmp(cluster_groups{2},'good')));
        good_templates = ismember(0:size(templates,1)-1,good_templates_idx);
        
        % Throw out all non-good template data
        templates = templates(good_templates,:,:);
        templateDepths = templateDepths(good_templates);
        waveforms = waveforms(good_templates,:);
        templateDuration = templateDuration(good_templates);
        templateDuration_us = templateDuration_us(good_templates);
        
        % Throw out all non-good spike data
        good_spike_idx = ismember(spike_templates,good_templates_idx);
        spike_times = spike_times(good_spike_idx);
        spike_templates = spike_templates(good_spike_idx);
        template_amplitudes = template_amplitudes(good_spike_idx);
        spikeDepths = spikeDepths(good_spike_idx);
        spike_times_timeline = spike_times_timeline(good_spike_idx);
        
        % Re-name the spike templates according to the remaining templates
        % (and make 1-indexed from 0-indexed)
        new_spike_idx = nan(max(spike_templates)+1,1);
        new_spike_idx(good_templates_idx+1) = 1:length(good_templates_idx);
        spike_templates = new_spike_idx(spike_templates+1);
        
    elseif ~exist('cluster_groups','var')
        if verbose; disp('Clusters not yet sorted'); end;
    end
    
end

%% Estimate striatal boundaries on probe

if ephys_exists && load_parts.ephys
    if verbose; disp('Estimating striatum boundaries on probe...'); end;
    
    %%% Get correlation of MUA and LFP
    n_corr_groups = 40;
    depth_group_edges = linspace(0,max(channel_positions(:,2)),n_corr_groups+1);
    depth_group = discretize(templateDepths,depth_group_edges);
    depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);
    unique_depths = 1:length(depth_group_edges)-1;
    
    spike_binning = 0.01; % seconds
    corr_edges = spike_times_timeline(1):spike_binning:spike_times_timeline(end);
    corr_centers = corr_edges(1:end-1) + diff(corr_edges);
    
    binned_spikes_depth = zeros(length(unique_depths),length(corr_edges)-1);
    for curr_depth = 1:length(unique_depths);
        binned_spikes_depth(curr_depth,:) = histcounts(spike_times_timeline( ...
            ismember(spike_templates,find(depth_group == unique_depths(curr_depth)))), ...
            corr_edges);
    end
    
    mua_corr = corrcoef(binned_spikes_depth');
    
    %%% Estimate start and end depths of striatum

    % end of striatum: biggest (smoothed) drop in MUA correlation near end
    groups_back = 15;
    mua_corr_end = medfilt2(mua_corr(end-groups_back+1:end,end-groups_back+1:end),[3,3]);
    mua_corr_end(triu(true(length(mua_corr_end)),0)) = nan;
    median_corr = medfilt1(nanmedian(mua_corr_end,2),3);
    [x,max_corr_drop] = min(diff(median_corr));
    str_end = depth_group_centers(end-groups_back+max_corr_drop);
    
    % start of striatum: look for ventricle (dropoff in templates)
    
    % (by template density)
    %     n_template_bins = 40;
    %     size_template_bins = max(channel_positions(:,2))/n_template_bins;
    %     template_density_bins = linspace(0,max(channel_positions(:,2)),n_template_bins);
    %     template_density = histcounts(templateDepths,template_density_bins);
    %
    %     str_end_bin = floor(str_end/size_template_bins);
    %
    %     n_bins_check = 3;
    %     bins_conv = ones(1,n_bins_check)/n_bins_check;
    %     template_gaps = conv(+(fliplr(template_density(1:str_end_bin)) < 2),bins_conv);
    %
    %     sorted_template_depths = sort([0;templateDepths]);
    %
    %     if any(template_gaps)
    %         str_gap_stop = length(template_gaps) - n_bins_check - find(template_gaps(n_bins_check:end),1);
    %         str_start = sorted_template_depths(find(sorted_template_depths > template_density_bins(str_gap_stop),1)) - 1;
    %     else
    %         str_start = sorted_template_depths(2);
    %     end
    
    % (by biggest gap)
    min_gap = 200;
    sorted_template_depths = sort([0;templateDepths]);
    [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
    if max_gap > min_gap
        str_start = sorted_template_depths(max_gap_idx+1)-1;
    else
        str_start = sorted_template_depths(2);
    end
    
    str_depth = [str_start,str_end];
    
    %%% Align striatal recordings using saved alignment
    n_aligned_depths = 6;
    ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
    ephys_align_fn = ['ephys_depth_align.mat'];
    if exist([ephys_align_path filesep ephys_align_fn],'file')
        load([ephys_align_path filesep ephys_align_fn]);
        % Only use if this animal was included in the alignment
        curr_animal_idx = strcmp(animal,{ephys_depth_align.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ephys_depth_align(curr_animal_idx).day);
            if any(curr_day_idx)
                if verbose; disp('Aligning striatum by saved depths...'); end;                
                
                % Get the maximum striatal length, split evenly
                all_str_lengths = diff(vertcat(ephys_depth_align(:).str_depth),[],2);
                max_str_length = max(all_str_lengths);
                global_str_depth_edges = linspace(0,max_str_length,n_aligned_depths+1);
                
                % Set current striatum depth edges by global length
                % assuming that the end is actually the end
                str_depth_edges = sort(str_depth(2) - global_str_depth_edges);
                
                % Get spike depths, setting all outside the striatum to NaN
                str_spikeDepths = spikeDepths;
                str_spikeDepths(spikeDepths < str_depth(1) | spikeDepths > str_depth(2)) = NaN;
                
                % Group the striatal spike depths into the um-standardized bins
                aligned_str_depth_group = discretize(str_spikeDepths,str_depth_edges);                
                
            end
        end
    end
        
    %%% Align by cortex-striatum kernel template
    % (this overrides the previous alignments)
    ephys_kernel_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
    ephys_kernel_align_fn = 'ephys_kernel_align.mat';
    if exist([ephys_kernel_align_path filesep ephys_kernel_align_fn],'file')
        load([ephys_kernel_align_path filesep ephys_kernel_align_fn]);
        % If alignment exists for this dataset, align
        curr_animal_idx = strcmp(animal,{ephys_kernel_align.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ephys_kernel_align(curr_animal_idx).days);
            if any(curr_day_idx)
                if verbose; disp('Aligning striatum by kernel alignment...'); end;
                % (use previously saved depth groups)
                aligned_str_depth_group = ephys_kernel_align(curr_animal_idx).aligned_str_depth_group{curr_day_idx};
                n_aligned_depths = ephys_kernel_align(curr_animal_idx).n_aligned_depths(curr_day_idx);
            end
        end
    end
    
%     % (to plot the aligned groups)    
%     [~,idx,~] = unique(spike_templates);
%     template_aligned_depth = aligned_str_depth_group(idx)+1;
%     template_aligned_depth(isnan(template_aligned_depth)) = 1;
%     
%     col = [0,0,1;copper(n_aligned_depths)];
%     
%     figure;plotSpread(templateDepths,'distributionIdx', ...
%         template_aligned_depth,'distributionColors',col(unique(template_aligned_depth),:));
%     set(gca,'YDir','reverse');
%     line(xlim,[str_depth(1),str_depth(1)]);
%     line(xlim,[str_depth(2),str_depth(2)]);
    
end

%% Classify spikes

if ephys_exists && load_parts.ephys && exist('cluster_groups','var')
    if verbose; disp('Classifying spikes...'); end;
    
    str_templates = templateDepths >= str_depth(1) & templateDepths <= str_depth(2);
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
    bin_spikes = nan(max(spike_templates), ...
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
    prop_long_isi = nan(max(spike_templates),1);
    cv2 = nan(max(spike_templates),1);
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
    
end


%% Finished
if verbose; disp('Finished loading experiment.'); end








