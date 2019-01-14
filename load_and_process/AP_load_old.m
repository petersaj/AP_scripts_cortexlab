% AP_load_old(animal,day,experiment,site)
%
% Loads old data (Phase 2 neuropixels?)

%% Old MPEP mouse name
% Old MPEP required a bad naming convention, check for files using this
mpep_animal = ['M111111_' animal];

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

%% Load timeline

[timeline_filename,timeline_exists] = AP_cortexlab_filename(animal,day,experiment,'timeline');

% (check old MPEP name if nonexistant)
if ~timeline_exists
    [timeline_filename,timeline_exists] = AP_cortexlab_filename(mpep_animal,day,experiment,'timeline');
end

if timeline_exists
    if verbose; disp('Loading timeline...'); end;
    
    load(timeline_filename);
    
    % Set rig-specific timeline names
    cam_name = 'pcoExposure';
    acqLive_name = 'acqLive';
    
    % Set rig-specific timeline names
    if any(strcmp('pcoExposure',{Timeline.hw.inputs.name}))
        rig = 'kilotrode';
    else
        rig = 'bigrig';
    end
    
    switch rig
        case 'bigrig'
            cam_name = 'cam2';
            acqLive_name = 'acqLiveEcho';
        case 'kilotrode'
            cam_name = 'pcoExposure';
            acqLive_name = 'acqLive';
    end
    
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
    
    if verbose; disp('Loading mpep protocol...'); end
    
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
    
    % If more offsets than offsets, clear the first one
    if photodiode.values(1) == 0
        photodiode_offsets(1) = [];
    end
    
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
    
    % SPECIFIC
    [~,expDef] = fileparts(block.expDef);
    if strcmp(expDef,'vanillaChoiceworld');
        % dumb signals thing, fix
        signals_events.hitValues = circshift(signals_events.hitValues,[0,-1]);
        signals_events.missValues = circshift(signals_events.missValues,[0,-1]);
        
        % Get stim times by closest photodiode flip
        [~,closest_stimOn_photodiode] = ...
            arrayfun(@(x) min(abs(signals_events.stimOnTimes(x) - ...
            photodiode_flip_times)), ...
            1:length(signals_events.stimOnTimes));
        stimOn_times = photodiode_flip_times(closest_stimOn_photodiode);
        
    elseif strcmp(expDef,'AP_visAudioPassive')
        min_stim_downtime = 0.5; % minimum time between pd flips to get stim
        stimOn_times_pd = photodiode_flip_times([true;diff(photodiode_flip_times) > min_stim_downtime]);
        stimOff_times_pd = photodiode_flip_times([diff(photodiode_flip_times) > min_stim_downtime;true]);
        warning('visAudioPassive: THIS IS TEMPORARY BECAUSE NO BUFFER TIME')
        
        stimOn_times = nan(size(signals_events.visualOnsetTimes));
        stimOn_times(end-(length(stimOn_times_pd)-1):end) = stimOn_times_pd;
        
        stimOff_times = nan(size(signals_events.visualOnsetTimes));
        stimOff_times(end-(length(stimOff_times_pd)-1):end) = stimOff_times_pd;
        
        % sanity check
        if length(signals_events.visualOnsetValues) ~= length(stimOn_times)
            error('Different number of signals/timeline stim ons')
        end
        
    else
        % Specialized: get stim on/off times BY ASSUMING MINIMUM STIM DOWNTIME
        min_stim_downtime = 1; % minimum time between pd flips to get stim
        stimOn_times = photodiode_flip_times([true;diff(photodiode_flip_times) > min_stim_downtime]);
        stimOff_times = photodiode_flip_times([diff(photodiode_flip_times) > min_stim_downtime;true]);
        % sanity check
        if length(signals_events.stimOnTimes) ~= length(stimOn_times)
            error('Different number of signals/timeline stim ons')
        end
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
        if verbose; disp('Loading facecam...'); end
        
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

% (check old MPEP name if nonexistant)
if ~data_path_exists
    [data_path,data_path_exists] = AP_cortexlab_filename(mpep_animal,day,experiment,'imaging',site);
end

if data_path_exists && load_parts.imaging
    if verbose; disp('Loading imaging data...'); end
    
    % Get the imaged colors
    spatialComponents_fns = dir([data_path filesep 'svdSpatialComponents*']);
    spatialComponents_names = {spatialComponents_fns.name};
    
    cam_color_n = length(spatialComponents_names);
    cam_color_signal = 'blue';
    cam_color_hemo = 'purple';
    
    if cam_color_n == 1
        
        experiment_path = [data_path filesep num2str(experiment)];
        
        frame_t_dir = dir([experiment_path filesep 'svdTemporalComponents_*.timestamps.npy']);
        U_dir = dir([data_path filesep 'svdSpatialComponents_*.npy']);
        V_dir = dir([experiment_path filesep 'svdTemporalComponents_*.npy']);
        
        frame_t = readNPY([frame_t_dir.folder filesep frame_t_dir.name]);
        U = readUfromNPY([U_dir.folder filesep U_dir.name]);
        V = readVfromNPY([V_dir(1).folder filesep V_dir(1).name]);
        
        framerate = 1./nanmedian(diff(frame_t));
        
        % Detrend and high-pass filter
        highpassCutoff = 0.01; % Hz
        [b100s, a100s] = butter(2, highpassCutoff/(framerate/2), 'high');
        dV = detrend(V', 'linear')';
        fV = single(filtfilt(b100s,a100s,double(dV)')');
        
        avg_im_dir = dir([data_path filesep 'meanImage_*.npy']);
        avg_im = readNPY([avg_im_dir.folder filesep avg_im_dir.name]);
        
    elseif cam_color_n == 2
        
        % Load in all things as neural (n) or hemodynamic (h)
        experiment_path = [data_path filesep num2str(experiment)];
        
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
        
        hemo_tform_fn = [data_path filesep 'hemo_tform.mat'];
        if exist(hemo_tform_fn,'file')
            % If the hemo tform matrix has been computed, load and fix
            if verbose; disp('Using old hemo tform...'); end
            load(hemo_tform_fn)
            zVh_Un = bsxfun(@minus, Vh_Un, mean(Vh_Un));
            Vn_hemo = transpose(Vn_th' - zVh_Un'*hemo_tform');
        else
            % If no p hemo tform matrix, compute and save
            if verbose; disp('Computing hemo tform...'); end
            %hemo_freq = [0.2,3];
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
    
    if verbose; disp('Loading ephys...'); end
    
    acqLive_channel = 1;
    load_lfp = false;
    
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
    % (true only for phase 2)
    flipped_banks = true;
    if flipped_banks
        channel_map = [channel_map(61:end);channel_map(1:60)];
        channel_positions = [channel_positions(61:end,:);channel_positions(1:60,:)];
    end
    
    % Default channel map/positions are from end: make from surface
    channel_positions(:,2) = max(channel_positions(:,2)) - channel_positions(:,2);
    
    % Load LFP (random snippet: just for correlation)
    n_channels = str2num(header.n_channels);
    % (this is where the old LFP is)
    lfp_filename = [ephys_path 'lfp.dat']; 
    % (this is where the new LFP is)
%     [data_path,data_path_exists] = AP_cortexlab_filename(animal,day,experiment,'ephysraw',site);
%     lfp_dir = dir([data_path 'experiment*_100-1_0.dat']);
%     lfp_filename = [data_path lfp_dir.name];
    if load_lfp && exist(lfp_filename,'file')
        
        if isfield(header,'lfp_sample_rate')
            lfp_sample_rate = str2num(header.lfp_sample_rate);
            lfp_cutoff = str2num(header.filter_cutoff);
            lfp_downsamp = (lfp_sample_rate/lfp_cutoff)/2;
        elseif isfield(header,'sample_rate')
            lfp_sample_rate = str2num(header.sample_rate);
            lfp_cutoff = str2num(header.lfp_cutoff);
            lfp_downsamp = (lfp_sample_rate/lfp_cutoff)/2;
        end       
        
        fid = fopen(lfp_filename);
        % (to load all LFP)
        lfp_all = fread(fid,[n_channels,Inf],'int16');
        % (to only load a snippet of LFP)
%         fseek(fid,(lfp_sample_rate*60*10*n_channels),'bof'); % move to 10 minutes after recording start
%         lfp_all = fread(fid,[n_channels,1e6],'int16'); % pull snippet
        fclose(fid);
        % eliminate non-connected channels
        lfp = lfp_all(channel_map+1,:);
        clear lfp_all;
        % get time of LFP sample points (NOTE: this is messy, based off of sample
        % rate and knowing what kwik2dat does, not sure how accurate)
        lfp_t = ([1:size(lfp,2)]*lfp_downsamp)/lfp_sample_rate;
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
    spike_times_timeline = interp1(acqlive_ephys_currexpt,acqLive_timeline,spike_times,'linear','extrap');
    if load_lfp && exist(lfp_filename,'file')
        lfp_t_timeline = interp1(acqlive_ephys_currexpt,acqLive_timeline,lfp_t,'linear','extrap');
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
        
        if verbose; disp('Removing non-good templates'); end
        
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
        if verbose; disp('Clusters not yet sorted'); end
    end
    
end

%% Finished
if verbose; disp('Finished loading experiment.'); end








