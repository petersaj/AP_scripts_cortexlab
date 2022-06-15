% Loads data from experiments
%
% Settings:
% (what to load)
% load_parts.(cam/imaging/ephys)
%
% (ephys)
% site (if multiple probes)
% recording (if multiple on/off recordings within probe)

% Turn warnings off
warning on

%% Display progress or not
if ~exist('verbose','var')
    verbose = false;
end

%% Define what to load

% Site (multiple probes) is optional
if ~exist('site','var')
    site = [];
end

% Recording (multiple recordings on one probe) is optional
if ~exist('recording','var')
    recording = [];
end

% If nothing specified, load everything (but not LFP)
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
    
    % Get wheel position and velocity
    rotaryEncoder_idx = strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder');
    % (this is a very strange hack to overcome a problem in the rotary
    % encoder that's known in the lab and was put on the wiki)
    wheel_position = Timeline.rawDAQData(:,rotaryEncoder_idx);
    wheel_position(wheel_position > 2^31) = wheel_position(wheel_position > 2^31) - 2^32;    
    [wheel_velocity,wheel_move] = AP_parse_wheel(wheel_position,Timeline.hw.daqSampleRate);
    
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
    stimScreen_on = Timeline.rawDAQData(:,photodiode_idx) > 0.2;
    stimScreen_on_t = Timeline.rawDAQTimestamps(stimScreen_on);
    photodiode_thresh = 2; % old: max(Timeline.rawDAQData(:,photodiode_idx))/2
    photodiode_trace = Timeline.rawDAQData(stimScreen_on,photodiode_idx) > photodiode_thresh;
    % (medfilt because photodiode can be intermediate value when backlight
    % coming on)
    % (OLD: this worked fine until photodiode bug: sometimes gray)
%     photodiode_trace_medfilt = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
%         photodiode_idx),3) > photodiode_thresh;
%     photodiode_flip = find((~photodiode_trace_medfilt(1:end-1) & photodiode_trace_medfilt(2:end)) | ...
%         (photodiode_trace_medfilt(1:end-1) & ~photodiode_trace_medfilt(2:end)))+1;
%     photodiode_flip_times = stimScreen_on_t(photodiode_flip)';
    % (NEW: accomodating photodiode bug flipping sometimes to gray)
    photodiode_trace_medfilt = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
        photodiode_idx),3);
    photodiode_diff_thresh = range(Timeline.rawDAQData(:,photodiode_idx))*0.2;
    photodiode_diff_t = 50; % time (in ms) to get delayed differential
    photodiode_diff_samples = round(Timeline.hw.daqSampleRate/1000*photodiode_diff_t);
    photodiode_diff_filt = [1,zeros(1,photodiode_diff_samples),-1];
    photodiode_trace_diff = abs(conv(photodiode_trace_medfilt,photodiode_diff_filt,'valid')) > ...
        photodiode_diff_thresh;
    photodiode_flip = find(~photodiode_trace_diff(1:end-1) & ...
        photodiode_trace_diff(2:end))+ photodiode_diff_samples + 1;
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
                  
        if length(reward_t_block) == length(reward_t_timeline) && ...
                ~isempty(reward_t_block)
            % If same number of timeline/block rewards, use those for sync
            timeline2block = reward_t_timeline;
            block2timeline = reward_t_block;
            
        elseif length(reward_t_block) ~= length(reward_t_timeline) && ...
                ~isempty(reward_t_block)
            % If there's a different number of block and timeline rewards (aka
            % manual rewards were given), try to fix this            
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
                
                timeline2block = reward_t_timeline;
                block2timeline = reward_t_block;
            else
                % otherwise, you're in trouble
                error('Manual rewards included - couldn''t match to block');
            end
            
        elseif isempty(reward_t_block)
            % If no rewards: probably much less robust, but use stim on
            % times and photodiode flips
            warning('No rewards, aligning block with estimated stimOn times');
            stim_t_offset = block.events.stimOnTimes' - photodiode_flip_times';
            blunt_stim_offset = mode(round(stim_t_offset(:)*10))/10;
            stim_t_offset_shift = stim_t_offset - blunt_stim_offset;
            t_offset_tolerance = 0.05;
            stim_t_offset_binary = abs(stim_t_offset_shift) < t_offset_tolerance;
            stim_tl_block_matched = sum(stim_t_offset_binary,2) == 1;           
            if nanmean(stim_tl_block_matched) > 0.9
                % (if 90% stim matched, use this)
                [~,estimated_stimOn_idx] = max(stim_t_offset_binary,[],2);
                
                timeline2block = photodiode_flip_times(estimated_stimOn_idx(stim_tl_block_matched));
                block2timeline = block.events.stimOnTimes(stim_tl_block_matched);
            else
                % (if >10% unmatched, don't use)
                error('Attempted stim on alignment - bad match');
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
                interp1(block2timeline,timeline2block,block.events.(block_fieldnames{curr_times}),'linear','extrap');
        end
    end
    
    % Pull signals specific to protocol
    AP_pull_signals;

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
        % (from facemap)
        [eyecam_processed_filename,eyecam_processed_exists] = AP_cortexlab_filename(animal,day,experiment,'eyecam_processed');
        if eyecam_processed_exists
            eyecam = load(eyecam_processed_filename);
        end
        % (from deeplabcut)
        [eyecam_dlc_filename,eyecam_dlc_exists] = AP_cortexlab_filename(animal,day,experiment,'eyecam_dlc');
        if eyecam_dlc_exists
            eyecam_dlc = AP_load_dlc(eyecam_dlc_filename);
            % (pupil position = average half-way between L/R and T/B)
            pupil_position = ...
                ([eyecam_dlc.pupil_right.x + eyecam_dlc.pupil_left.x, ...
                eyecam_dlc.pupil_right.y + eyecam_dlc.pupil_left.y]./2 + ...
                [eyecam_dlc.pupil_top.x + eyecam_dlc.pupil_bot.x, ...
                eyecam_dlc.pupil_top.y + eyecam_dlc.pupil_bot.y]./2)./2;
            
            % (pupil diameter = average distance between L/R and T/B)
            pupil_diameter = ...
                (vecnorm([eyecam_dlc.pupil_right.x - eyecam_dlc.pupil_left.x, ...
                eyecam_dlc.pupil_right.y - eyecam_dlc.pupil_left.y]') + ...
                vecnorm([eyecam_dlc.pupil_top.x - eyecam_dlc.pupil_bot.x, ...
                eyecam_dlc.pupil_top.y - eyecam_dlc.pupil_bot.y]'))/2;
            
            % (blinks = uncertain pupil points)
            eyecam_dlc_params.surroundingBlinks = 5;
            eyecam_dlc_params.thresLikelihood = 0.8;
            pupil_likelihood = [ ...
                eyecam_dlc.pupil_bot.likelihood, ...
                eyecam_dlc.pupil_left.likelihood, ...
                eyecam_dlc.pupil_right.likelihood, ...
                eyecam_dlc.pupil_top.likelihood];
            uncertain_pupil_points = conv(sum(pupil_likelihood > ...
                eyecam_dlc_params.thresLikelihood, 2) <= 2, ...
                ones(2*eyecam_dlc_params.surroundingBlinks+1,1), 'same') > 0;
            
            % (NaN-out blinks)
            pupil_position(uncertain_pupil_points,:) = NaN;
            pupil_diameter(uncertain_pupil_points) = NaN;
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
                if abs(n_eyecam_frames_syncd_movie - n_eyecam_frames_syncd_timeline) > 2
                    warning('[%s %s %d] Eyecam: %d frames video vs timeline', ...
                        animal,day,experiment,n_eyecam_frames_syncd_movie - n_eyecam_frames_syncd_timeline);
                end
                
                % Get times of cam frames in timeline
                eyecam_t = nan(n_eyecam_frames,1);
                eyecam_t(eyecam_frame_idx(eyecam_frame_idx > 0)) = eyeCamStrobe_up_t(eyecam_frame_idx > 0);
                
                try
                    save(eyecam_t_savefile,'eyecam_t');
                catch me
                    warning (['Can''t save: ' eyecam_t_savefile])
                end
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
                if abs(n_facecam_frames_syncd_movie - n_facecam_frames_syncd_timeline) > 2
                    warning('[%s %s %d] Facecam: %d frames video vs timeline', ...
                        animal,day,experiment,n_facecam_frames_syncd_movie - n_facecam_frames_syncd_timeline);
                end
                
                % Get times of cam frames in timeline
                facecam_t = nan(n_facecam_frames,1);
                facecam_t(facecam_frame_idx(facecam_frame_idx > 0)) = faceCamStrobe_up_t(facecam_frame_idx > 0);
                
                try
                    save(facecam_t_savefile,'facecam_t');
                catch me
                    warning (['Can''t save: ' facecam_t_savefile])
                end
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

[data_path,data_path_exists] = AP_cortexlab_filename(animal,day,experiment,'imaging');
experiment_path = [data_path filesep num2str(experiment)];

% (check for specific imaging file since data path is just root)
spatialComponents_fns = dir([data_path filesep 'svdSpatialComponents*']);
imaging_exists = ~isempty(spatialComponents_fns);

if imaging_exists && load_parts.imaging
    if verbose; disp('Loading imaging data...'); end
    
    % Get the imaging file locations
    spatialComponents_dir = dir([data_path filesep 'svdSpatialComponents*']);
    meanImage_dir = dir([data_path filesep 'meanImage*']);
    
    % Get colors
    cam_color_n = length(spatialComponents_dir);
    
    cam_color_names = regexp([spatialComponents_dir.name],'svdSpatialComponents_(\w*).npy','tokens');
    
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
        avg_im_n = readNPY([data_path filesep 'meanImage_' cam_color_signal '.npy']);
        
        Uh = readUfromNPY([data_path filesep 'svdSpatialComponents_' cam_color_hemo '.npy']);
        Vh = readVfromNPY([experiment_path filesep 'svdTemporalComponents_' cam_color_hemo '.npy']);
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
        
        % Shift neural frames to hemo timestamps
        Vn_th = SubSampleShift(Vn,1,2);
        
        % Re-cast hemo V's into neural U's
        Vh_Un = ChangeU(Uh,Vh,Un);
        
        % Make/load hemo tform
        hemo_tform_fn = [experiment_path filesep 'hemo_tform.mat'];
        if exist(hemo_tform_fn,'file')
            % If the hemo tform matrix has been computed, load and fix
            if verbose; disp('Using old hemo tform...'); end
            load(hemo_tform_fn)
            zVh_Un = bsxfun(@minus, Vh_Un, nanmean(Vh_Un,2));
            Vn_hemo = transpose(Vn_th' - zVh_Un'*hemo_tform');
        else
            % If no hemo tform matrix, compute and save
            if verbose; disp('Computing hemo tform...'); end
            hemo_freq = [5,15];
            skip_seconds = 20; % the beginning and end can be wonky
            skip_frames = 1 + round(skip_seconds*framerate);
            [~,hemo_tform] = HemoCorrectLocal(Un, ...
                Vn_th(:,skip_frames:end-skip_frames-1), ...
                Vh_Un(:,skip_frames:end-skip_frames-1), ...
                framerate,hemo_freq,3);
            
            zVh_Un = bsxfun(@minus, Vh_Un, nanmean(Vh_Un,2));
            Vn_hemo = transpose(Vn_th' - zVh_Un'*hemo_tform');
            
            try
                save(hemo_tform_fn,'hemo_tform');
            catch me
                warning(['Can''t save ' hemo_tform_fn]);
            end
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

% Pick quality control
if ~exist('ephys_quality_control','var')
    ephys_quality_control = true;
end

% Pick kilosort version (2 by default, 1 old if selected)
if ~exist('kilosort_version','var') || kilosort_version == 2
    [ephys_path,ephys_exists] = AP_cortexlab_filename(animal,day,experiment,'ephys',site,recording);
elseif exist('kilosort_version','var') && kilosort_version == 1
    [ephys_path,ephys_exists] = AP_cortexlab_filename(animal,day,experiment,'ephys_ks1',site,recording);
end

if ephys_exists && load_parts.ephys
    
    if verbose; disp('Loading ephys...'); end
    
    % These are the digital channels going into the FPGA
    photodiode_sync_idx = 1;
    acqLive_sync_idx = 2;
    led_sync_idx = 3;
    flipper_sync_idx = 4;
    
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
    template_chan_amp_thresh = max(template_chan_amp,[],2)*0.5;
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
    
    % Get sync points for alignment
    
    % Get experiment index by finding numbered folders
    protocols_list = AP_list_experiments(animal,day);
    experiment_idx = experiment == [protocols_list.experiment];
    
    if exist('flipper_flip_times_timeline','var') && length(sync) >= flipper_sync_idx
        % (if flipper, use that)
        % (at least one experiment the acqLive connection to ephys was bad
        % so it was delayed - ideally check consistency since it's
        % redundant)
        bad_flipper = false;
        
        % Get flipper experiment differences by long delays
        % (note: this is absolute difference, if recording stopped and
        % started then the clock starts over again, although I thought it
        % wasn't supposed to when I grab the concatenated sync, so
        % something might be wrong)
        flip_diff_thresh = 1; % time between flips to define experiment gap (s)
        flipper_expt_idx = [1;find(abs(diff(sync(flipper_sync_idx).timestamps)) > ...
            flip_diff_thresh)+1;length(sync(flipper_sync_idx).timestamps)+1];
        
        flipper_flip_times_ephys = sync(flipper_sync_idx).timestamps( ...
            flipper_expt_idx(find(experiment_idx)):flipper_expt_idx(find(experiment_idx)+1)-1);
        
        % Pick flipper times to use for alignment
        if length(flipper_flip_times_ephys) == length(flipper_flip_times_timeline)
            % If same number of flips in ephys/timeline, use all
            sync_timeline = flipper_flip_times_timeline;
            sync_ephys = flipper_flip_times_ephys;
        elseif length(flipper_flip_times_ephys) ~= length(flipper_flip_times_timeline)
            % If different number of flips in ephys/timeline, best
            % contiguous set via xcorr of diff
            warning([animal ' ' day ':Flipper flip times different in timeline/ephys'])
            warning(['The fix for this is probably not robust: always check'])
            [flipper_xcorr,flipper_lags] = ...
                xcorr(diff(flipper_flip_times_timeline),diff(flipper_flip_times_ephys));
            [~,flipper_lag_idx] = max(flipper_xcorr);
            flipper_lag = flipper_lags(flipper_lag_idx);
            % (at the moment, assuming only dropped from ephys)
            sync_ephys = flipper_flip_times_ephys;
            sync_timeline = flipper_flip_times_timeline(flipper_lag+1: ...
                flipper_lag+1:flipper_lag+length(flipper_flip_times_ephys));
        end
        
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
    
    % Get spike times in timeline time
    spike_times_timeline = interp1(sync_ephys,sync_timeline,spike_times,'linear','extrap');
    
    % Get "good" templates from labels if quality control selected and
    % manual labels exist
    if ephys_quality_control && exist('cluster_groups','var')
        % If there's a manual classification
        if verbose; disp('Keeping manually labelled good units...'); end
        
        % Check that all used spike templates have a label
        spike_templates_0idx_unique = unique(spike_templates_0idx);
        if ~all(ismember(spike_templates_0idx_unique,uint32(cluster_groups{1}))) || ...
                ~all(ismember(cluster_groups{2},{'good','mua','noise'}))
            warning([animal ' ' day ': not all templates labeled']);
        end
        
        % Define good units from labels
        good_templates_idx = uint32(cluster_groups{1}( ...
            strcmp(cluster_groups{2},'good') | strcmp(cluster_groups{2},'mua')));
        good_templates = ismember(0:size(templates,1)-1,good_templates_idx);        
    else
        % If no cluster groups at all, keep all
        warning([animal ' ' day ' - no ephys quality control']);
        if verbose; disp('No ephys quality control, keeping all and re-indexing'); end
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

%% Load LFP
% (either single channel full or all channel snippet)

if ephys_exists && load_parts.ephys && exist('lfp_channel','var')
    
    % Get LFP file info
    n_channels = str2num(header.n_channels);
    [data_path,data_path_exists] = AP_cortexlab_filename(animal,day,experiment,'ephys_dir',site);
    
    % (get all recordings within site - assume concat at this point)
    lfp_recordings = dir([data_path filesep 'experiment*']);
    lfp_filenames = cellfun(@(x) ...
        [data_path filesep x filesep 'recording1\continuous\Neuropix-3a-100.1\continuous.dat'], ...
        {lfp_recordings.name},'uni',false);
    
    % Get LFP properties
    % (NOTE: LFP channel map is different from kilosort channel map because
    % kilosort2 drops channels without spikes)
    channel_map_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\kilosort_channelmaps\forPRBimecP3opt3.mat';
    channel_map_full = load(channel_map_fn);
    max_depth = 3840;
    lfp_channel_positions = max_depth - channel_map_full.ycoords;
    lfp_sample_rate = str2num(header.lfp_sample_rate);
    lfp_cutoff = str2num(header.filter_cutoff);
    
    % Memory map LFP
    n_bytes = 2; % LFP = int16 = 2 bytes
    n_lfp_samples = nan(size(lfp_filenames));
    lfp_memmap = cell(size(lfp_filenames));
    for curr_lfp_filename = 1:length(lfp_filenames)
        lfp_fileinfo = dir(lfp_filenames{curr_lfp_filename});
        n_lfp_samples(curr_lfp_filename) = lfp_fileinfo.bytes/n_bytes/n_channels;
        lfp_memmap{curr_lfp_filename} = ...
            memmapfile(lfp_filenames{curr_lfp_filename}, ...
            'Format',{'int16',[n_channels n_lfp_samples(curr_lfp_filename)],'lfp'});
    end
            
    if isnumeric(lfp_channel)
        
        % Load LFP of whole current experiment from one channel
        if verbose; disp(['Loading LFP (channel ' num2str(lfp_channel) ')...']); end;
        
        % Load single LFP channel within experiment bounds
        % (treat as concatenated if multiple files)
        lfp_load_start = round((lfp_sample_rate*sync_ephys(1)));
        lfp_load_stop = round((lfp_sample_rate*sync_ephys(end)));
        
        lfp_concat_length = [0,cumsum(n_lfp_samples)];
        lfp_load_file = find(lfp_load_start < lfp_concat_length,1)-1;
        lfp_load_start_rel = lfp_load_start - lfp_concat_length(lfp_load_file);
        lfp_load_stop_rel = lfp_load_stop - lfp_concat_length(lfp_load_file);
        
        lfp = double(lfp_memmap{lfp_load_file}.Data.lfp(lfp_channel,lfp_load_start_rel:lfp_load_stop_rel));
        
        % Get LFP times and convert to timeline time
        lfp_load_start_t = lfp_load_start/lfp_sample_rate;
        lfp_t = [0:size(lfp,2)-1]/lfp_sample_rate + lfp_load_start_t;
        lfp_t_timeline = interp1(sync_ephys,sync_timeline,lfp_t,'linear','extrap');
        
        %%% Remove light artifact
        if verbose; disp('Cleaning LFP...'); end;
        
        % Get light times (assume blue/violet alternate)
        light_t_timeline = interp1(sync_ephys,sync_timeline,sync(led_sync_idx).timestamps,'linear','extrap');
        use_light_times = light_t_timeline >= lfp_t_timeline(1) & light_t_timeline <= lfp_t_timeline(end);
        light_on = light_t_timeline(sync(led_sync_idx).values == 1 & use_light_times);
        light_off = light_t_timeline(sync(led_sync_idx).values == 0 & use_light_times);
        
        % (cut uncoupled off/on from start/end)
        if light_off(1) < light_on(1)
            light_off(1) = [];
        end
        if light_on(end) > light_off(end)
            light_on(end) = [];
        end
        
        blue_on = light_on(1:2:end);
        blue_off = light_off(1:2:end);
        violet_on = light_on(2:2:end);
        violet_off = light_off(2:2:end);
        
        light_on_mean = mean(light_off - light_on);
        light_off_mean = mean(light_on(2:end) - light_off(1:end-1));
        light_surround_t = [-(light_off_mean/2):1/lfp_sample_rate:(light_on_mean+(light_off_mean/2))];
        
        % Pull out LFP around light on
        use_blue_on = blue_on >= lfp_t_timeline(1) & blue_on <= lfp_t_timeline(end);
        blue_on_pull_t = blue_on(use_blue_on) + light_surround_t;
        blue_on_lfp = interp1(lfp_t_timeline,lfp',blue_on_pull_t);
        
        use_violet_on = violet_on >= lfp_t_timeline(1) & violet_on <= lfp_t_timeline(end);
        violet_on_pull_t = violet_on(use_violet_on) + light_surround_t;
        violet_on_lfp = interp1(lfp_t_timeline,lfp',violet_on_pull_t);
        
        % Subtract baseline
        baseline_t = find(light_surround_t < 0,1,'last');
        blue_on_lfp_baselinesub = blue_on_lfp - blue_on_lfp(:,baseline_t,:);
        violet_on_lfp_baselinesub = violet_on_lfp - violet_on_lfp(:,baseline_t,:);
        
        % Get rolling median (allow light artifact to change slightly)
        n_light = 500;
        blue_on_lfp_baselinesub_med = movmedian(blue_on_lfp_baselinesub,n_light,1);
        violet_on_lfp_baselinesub_med = movmedian(violet_on_lfp_baselinesub,n_light,1);
        
        % Interpolate out the artifact to remove
        n_lfp_channels = size(lfp,1);
        blue_light_remove = interp1( ...
            reshape(permute(blue_on_pull_t,[2,1]),[],1), ...
            reshape(permute(blue_on_lfp_baselinesub_med,[2,1,3]),[],n_lfp_channels), ...
            reshape(lfp_t_timeline,[],1))';
        violet_light_remove = interp1( ...
            reshape(permute(violet_on_pull_t,[2,1]),[],1), ...
            reshape(permute(violet_on_lfp_baselinesub_med,[2,1,3]),[],n_lfp_channels), ...
            reshape(lfp_t_timeline,[],1))';
        
        % Zero-out any NaNs (e.g. remove nothing)
        blue_light_remove(isnan(blue_light_remove)) = 0;
        violet_light_remove(isnan(violet_light_remove)) = 0;
        
        % Remove the artifact
        lfp_lightfix = lfp - (blue_light_remove + violet_light_remove);
        
        % NOT DOING THIS: IS THIS NECESSARY? TOO MUCH MEMORY
        %     % (low-pass filter: sometimes bunch of junk at high freq?)
        %     freqCutoff = 300; % Hz
        %     [b100s, a100s] = butter(2,freqCutoff/(lfp_sample_rate/2),'low');
        %     lfp_lightfix = single(filtfilt(b100s,a100s,double(lfp_lightfix)')');
        
    elseif strcmp(lfp_channel,'all')
        
        % Load short LFP segment (from start = no light) from all channels
        if verbose; disp('Loading LFP (all channels snippet)...'); end;
        
        % Choose snippet of recording time before first experiment (no light)
        t_load = 10; % time to load (in seconds)
        t_load_pre_exp = 1; % time before first experiment to load up to
        experiment_ephys_starts = sync(acqLive_sync_idx).timestamps(sync(acqLive_sync_idx).values == 1);
        lfp_load_start = round((lfp_sample_rate*(experiment_ephys_starts(1)-t_load_pre_exp-t_load)));
        lfp_load_stop = round((lfp_sample_rate*(experiment_ephys_starts(1)-t_load_pre_exp)));
        
        % Load all LFP channels in snippet (before recording = first file)
        lfp = lfp_memmap{1}.Data.lfp(:,lfp_load_start:lfp_load_stop);
        
        % Sort LFP so it goes from surface to depth
        [~,lfp_sort_idx] = sort(lfp_channel_positions);
        lfp_channel_positions = lfp_channel_positions(lfp_sort_idx);
        lfp = lfp(lfp_sort_idx,:);
        
        % Get LFP times and convert to timeline time
        lfp_load_start_t = lfp_load_start/lfp_sample_rate;
        lfp_t = [0:size(lfp,2)-1]/lfp_sample_rate + lfp_load_start_t;
        
        % Get power spectrum of LFP
        window_length = 2; % in seconds
        window_overlap = 1; % in seconds
        window_length_samples = round(window_length/(1/lfp_sample_rate));
        window_overlap_samples = round(window_overlap/(1/lfp_sample_rate));
        [lfp_power,lfp_power_freq] = pwelch(zscore(double(lfp),[],2)', ...
            window_length_samples,window_overlap_samples,[],lfp_sample_rate);
        
        if verbose
            figure;        
            imagesc(lfp_power_freq,lfp_channel_positions,log10(lfp_power'));
            xlabel('Frequency');
            ylabel('Depth (\mum)');
            c = colorbar;
            ylabel(c,'Log_{10} power');
            xlim([0,100]);
            colormap(hot);
            title('LFP power');            
        end        
    end    
end


%% Estimate area boundaries on probe
% (from old experiments, expects striatum alignment by default unless an
% 'ephys_align' variable is set)

if ephys_exists && load_parts.ephys && ~exist('ephys_align','var')
    if verbose; disp('Estimating striatum boundaries on probe...'); end
    
    % str_align = alignment method ('none', 'depth', or 'kernel')
    
    % requires n_aligned_depths for alignment, set default
    if ~exist('n_aligned_depths','var')
        n_aligned_depths = 3;
    end
    
    % if no alignment specified, default kernel
    if ~exist('str_align','var')
        str_align = 'kernel';
    end
    
    [str_depth,aligned_str_depth_group] = AP_align_striatum_ephys;
    
elseif ephys_exists && load_parts.ephys && exist('ephys_align','var') && ...
        strcmp(ephys_align,'cortex')    
    if verbose; disp('Estimating cortex boundaries on probe...'); end
    
    % Cortex alignment: assume top of probe out of brain and very
    % correlated, look for big drop in correlation from top-down
    lfp_corr = corrcoef(double(transpose(lfp-nanmedian(lfp,1))));
    lfp_corr_diag = lfp_corr;
    lfp_corr_diag(triu(true(size(lfp_corr)),0)) = NaN;
    lfp_corr_from_top = nanmean(lfp_corr_diag,2)';
    
    n_lfp_medfilt = 10;
    lfp_corr_from_top_medfilt = medfilt1(lfp_corr_from_top,n_lfp_medfilt);
    lfp_corr_from_top_medfilt_thresh = ...
        (max(lfp_corr_from_top_medfilt) - ...
        range(lfp_corr_from_top_medfilt)*0.2);
    ctx_start = lfp_channel_positions( ...
        find(lfp_corr_from_top_medfilt > lfp_corr_from_top_medfilt_thresh,...
        1,'last'));

    if verbose
        figure;
        imagesc(lfp_channel_positions,lfp_channel_positions,lfp_corr_diag);
        axis image
        colormap(brewermap([],'*RdBu'));
        caxis([-1,1])
        xlabel('Depth (\mum)');
        ylabel('Depth (\mum)');
        c = colorbar;
        ylabel(c,'Med. sub. correlation');
        line(xlim,[ctx_start,ctx_start],'color','k','linewidth',2);
        line([ctx_start,ctx_start],ylim,'color','k','linewidth',2);
        title(sprintf('%s %s: LFP correlation and cortex start',animal,day));
    end
            
    % (give leeway for cortex to start relative to first unit - if within
    % this leeway then back up cortex start, if outside then keep but throw
    % a warning)
    ctx_lfp_spike_leeway = 100; % um leeway for lfp/unit match
    ctx_lfp_spike_diff = ctx_start-min(template_depths);
    if ctx_lfp_spike_diff > ctx_lfp_spike_leeway
        warning('%s %s: LFP-estimated cortex start is after first unit %.0f um', ...
            animal,day,ctx_lfp_spike_diff);
    elseif ctx_lfp_spike_diff > 0 && ctx_lfp_spike_diff < ctx_lfp_spike_leeway
        ctx_start = min(template_depths)-1; % (back up to 1um before first unit)
    end  
    
    %%% If histology is aligned, get areas by depth
    [probe_ccf_fn,probe_ccf_fn_exists] = AP_cortexlab_filename(animal,[],[],'probe_ccf');
    if ~probe_ccf_fn_exists
        if verbose; disp('No histology alignment for probe...'); end
    elseif probe_ccf_fn_exists
        if verbose; disp('Estimating histology-aligned cortical areas on probe...'); end
        
        % (load the CCF structure tree)
        allen_atlas_path = fileparts(which('template_volume_10um.npy'));
        st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);
        
        % Load probe CCF alignment
        load(probe_ccf_fn);
        
        % Get area names and borders across probe
        [~,dv_sort_idx] = sort(probe_ccf.trajectory_coords(:,2));
        dv_voxel2um = 10*0.945; % CCF DV estimated scaling
        probe_trajectory_depths = ...
            pdist2(probe_ccf.trajectory_coords, ...
            probe_ccf.trajectory_coords((dv_sort_idx == 1),:))*dv_voxel2um;
        
        probe_depths = probe_trajectory_depths + ctx_start;
        
        % Get recorded areas and boundaries (ignore layer distinctions)
        probe_areas = unique(regexprep( ...
            st(probe_ccf.trajectory_areas(probe_depths > 0 & ...
            probe_depths < max(channel_positions(:,2))),:).safe_name, ...
            ' layer .*',''));
        
        probe_area_boundaries = cellfun(@(area) ...
            prctile(probe_depths(contains( ...
            st(probe_ccf.trajectory_areas,:).safe_name,area)),[0,100]), ...
            probe_areas,'uni',false);
    end
    
end

%% Classify spikes

if ephys_exists && load_parts.ephys
    if verbose; disp('Classifying spikes...'); end
    
    if exist('str_depth','var') && ~isempty(str_depth)
        % If striatum depths defined - split striatal/nonstriatal cells
        str_templates = template_depths >= str_depth(1) & template_depths <= str_depth(2);
    else
        % If not - just assume everything cortical
        str_templates = false(size(template_depths));
    end
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
        
        if any(str_templates)
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
        end
        
        % Plot depth vs. firing rate colored by cell type
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
        
        norm_spike_n = mat2gray(log10(accumarray(spike_templates,1)+1));
        
        figure('Position',[94,122,230,820]);
        gscatter(norm_spike_n,template_depths,celltypes,use_colors,[],10);
        xlim([0,1])
        set(gca,'YDir','reverse');
        xlabel('Norm log_{10} spike rate');
        ylabel('Depth (\mum)');
        legend(celltype_labels(plot_celltypes),'location','NW');
        ylim([0,max(channel_positions(:,2))])
        
        drawnow;
        
    end
    
end


%% Finished
if verbose; disp('Finished loading experiment.'); end








