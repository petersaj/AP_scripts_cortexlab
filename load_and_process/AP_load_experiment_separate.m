%% Define experiment

animal = 'AP017';
day = '2017-06-07';
experiment = '1';
rig = 'kilotrode'; % kilotrode or bigrig

% Get the number of colors imaged during the experiment
data_path = AP_cortexlab_filename(animal,day,experiment,'datapath');

spatialComponents_fns = dir([data_path filesep 'svdSpatialComponents*']);
spatialComponents_names = {spatialComponents_fns.name};

cam_color_n = length(spatialComponents_names);
cam_color_signal = 'blue';
cam_color_hemo = 'purple';

%% Load timeline

% Set rig-specific timeline names
switch rig
    case 'bigrig'
        cam_name = 'cam2';
        acqLive_name = 'acqLiveEcho';
    case 'kilotrode'
        cam_name = 'pcoExposure';
        acqLive_name = 'acqLive';
end

% Load timeline
timeline_filename = AP_cortexlab_filename(animal,day,experiment,'timeline');
load(timeline_filename);
timeline_sample_rate = Timeline.hw.daqSampleRate;

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

% Load in protocol and get photodiode, if there was an mpep experiment
if exist(fileparts(AP_cortexlab_filename(animal,day,experiment,'protocol')),'dir')
    
    % Load in hardware info
    hwinfo_filename = AP_cortexlab_filename(animal,day,experiment,'hardware');
    load(hwinfo_filename);
    
    % Get flicker or steady photodiode
    photodiode_type = myScreenInfo.SyncSquare.Type;
    
    try
        protocol_filename = AP_cortexlab_filename(animal,day,experiment,'protocol','8digit');
        load(protocol_filename);
    catch me
        protocol_filename = AP_cortexlab_filename(animal,day,experiment,'protocol','day_dash');
        load(protocol_filename);
    end
    
    % Get stimulus onsets and parameters
    photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
    
    % If the stim screen is flickered, interpolate photodiode when off
    disp(['Getting mpep data: ' rig ', photodiode: ' photodiode_type]);
    switch rig        
        case 'kilotrode'
            
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
            
        case 'bigrig'
            switch lower(photodiode_type)
                case 'flicker'
                    photodiode_thresh = max(Timeline.rawDAQData(:,photodiode_idx))/2;
                    photodiode_trace = Timeline.rawDAQData(:,photodiode_idx) > photodiode_thresh;
                    photodiode_flip = find((~photodiode_trace(1:end-1) & photodiode_trace(2:end)) | ...
                        (photodiode_trace(1:end-1) & ~photodiode_trace(2:end)))+1;
                    
                    photodiode = struct('timestamps',[],'values',[]);
                    photodiode.timestamps = Timeline.rawDAQTimestamps(photodiode_flip)';
                    photodiode.values = photodiode_trace(photodiode_flip);
                case 'steady'
                    error('write this');
            end
            
    end
     
    photodiode_onsets = photodiode.timestamps(photodiode.values == 1);
    
    refresh_rate_cutoff = 1/5;
    stim_onsets = photodiode_onsets( ...
        [1;find(diff(photodiode_onsets) > refresh_rate_cutoff) + 1]);
    
    stimIDs = zeros(size(stim_onsets));
    for q = 1:size(Protocol.seqnums,1)
        stimIDs(Protocol.seqnums(q,:)) = q;
    end
    
end

%% Load task/behavior

% Load the block
[block_filename, block_exists] = AP_cortexlab_filename(animal,day,experiment,'block');
load(block_filename);

% Get reward times in block and timeline
reward_t_block = block.outputs.rewardTimes;

timeline_reward_idx = strcmp({Timeline.hw.inputs.name}, 'rewardEcho');
reward_thresh = max(Timeline.rawDAQData(:,timeline_reward_idx))/2;
reward_trace = Timeline.rawDAQData(:,timeline_reward_idx) > reward_thresh;
reward_t_timeline = Timeline.rawDAQTimestamps(find(reward_trace(2:end) & ~reward_trace(1:end-1))+1);

% Go through all block events and convert to timeline time using the reward
% as the reference event
block_fieldnames = fieldnames(block.events);
block_values_idx = cellfun(@(x) ~isempty(x),strfind(block_fieldnames,'Values'));
block_times_idx = cellfun(@(x) ~isempty(x),strfind(block_fieldnames,'Times'));

choiceworld = block.events;
for curr_times = find(block_times_idx)'
    choiceworld.(block_fieldnames{curr_times}) = ...
        AP_clock_fix(block.events.(block_fieldnames{curr_times}),reward_t_block,reward_t_timeline);
end


%% Get wheel velocity and licking
% this is super preliminary

rotaryEncoder_idx = strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder');

wheel_interp = interp1(Timeline.rawDAQTimestamps,Timeline.rawDAQData(:,rotaryEncoder_idx),frame_t);
% subtract median filtered because of these crazy jumps sometimes??
wheel_interp_medfilt = medfilt1(wheel_interp,100);
wheel_interp_medfiltsub = wheel_interp - wheel_interp_medfilt;
% remove ridiculous outliers
wheel_interp_medfiltsub(abs(wheel_interp_medfiltsub) > 1e5) = 0;

wheel_velocity = [0;diff(smooth(wheel_interp_medfiltsub,10))];
wheel_speed = abs(hilbert(wheel_velocity))';

lickPiezo_idx = strcmp({Timeline.hw.inputs.name}, 'piezoLickDetector');
lickPiezo_interp = interp1(Timeline.rawDAQTimestamps,Timeline.rawDAQData(:,lickPiezo_idx),frame_t);
licking_trace = abs(hilbert(lickPiezo_interp));


%% Load eye tracking data (old, with etGUI)

% Load in processed eye tracking data
eyetracking_processed_filename = AP_cortexlab_filename(animal,day,experiment,'eyecam_processed');
eyetracking = load(eyetracking_processed_filename);

% Get frame times for eyetracking data
eyetracking_t = eyetracker_getFrameTimes_AP(animal,day,experiment);


%% Load face/eyecam processing (with eyeGUI)

% Get folder with camera data
cam_dir = fileparts(AP_cortexlab_filename(animal,day,experiment,'eyecam'));

% Get cam sync
camSync_idx = strcmp({Timeline.hw.inputs.name}, 'camSync');
camSync_thresh = max(Timeline.rawDAQData(:,camSync_idx))/2;
camSync = Timeline.rawDAQData(:,camSync_idx) > camSync_thresh;
camSync_up = find((~camSync(1:end-1) & camSync(2:end)))+1;

% EYECAM
% Load camera processed data
eyecam_processed_filename = AP_cortexlab_filename(animal,day,experiment,'eyecam_processed');
if exist(eyecam_processed_filename,'file')
    eyecam = load(eyecam_processed_filename);
end

% Get camera times
eyecam_fn = AP_cortexlab_filename(animal,day,experiment,'eyecam');
eyecam_t_savefile = [cam_dir filesep 'eyecam_t.mat'];

if exist(eyecam_fn,'file') && ~exist(eyecam_t_savefile,'file')   
    % Get facecam strobes
    eyeCamStrobe_idx = strcmp({Timeline.hw.inputs.name}, 'eyeCameraStrobe');
    eyeCamStrobe_thresh = max(Timeline.rawDAQData(:,eyeCamStrobe_idx))/2;
    eyeCamStrobe = Timeline.rawDAQData(:,eyeCamStrobe_idx) > eyeCamStrobe_thresh;
    eyeCamStrobe_up = find((~eyeCamStrobe(1:end-1) & eyeCamStrobe(2:end)))+1;
    eyeCamStrobe_up_t = Timeline.rawDAQTimestamps(eyeCamStrobe_up);
    
    % Get sync times for cameras (or load if already done)
    [eyecam_sync_frames,n_eyecam_frames] = AP_get_cam_sync_frames(eyecam_fn);
    
    % Get the closest facecam strobe to sync start, find offset and frame idx
    [~,eyecam_strobe_sync] = min(abs(camSync_up(1) - eyeCamStrobe_up));
    eyecam_frame_offset = eyecam_sync_frames(1) - eyecam_strobe_sync;
    eyecam_frame_idx = [1:length(eyeCamStrobe_up)] + eyecam_frame_offset;
    
    % Get times of facecam frames in timeline
    eyecam_t = nan(n_eyecam_frames,1);
    eyecam_t(eyecam_frame_idx) = eyeCamStrobe_up_t;
            
    save(eyecam_t_savefile,'eyecam_t');   
elseif exist(eyecam_fn,'file') && exist(eyecam_t_savefile,'file')   
    load(eyecam_t_savefile);
end

% FACECAM
% Load camera processed data
facecam_processed_filename = AP_cortexlab_filename(animal,day,experiment,'facecam_processed');
if exist(facecam_processed_filename,'file')
    facecam = load(facecam_processed_filename);
end

% Get camera times
facecam_fn = AP_cortexlab_filename(animal,day,experiment,'facecam');
facecam_t_savefile = [cam_dir filesep 'facecam_t.mat'];

if exist(facecam_fn,'file') && ~exist(facecam_t_savefile,'file')      
    % Get facecam strobes
    faceCamStrobe_idx = strcmp({Timeline.hw.inputs.name}, 'faceCamStrobe');
    faceCamStrobe_thresh = max(Timeline.rawDAQData(:,faceCamStrobe_idx))/2;
    faceCamStrobe = Timeline.rawDAQData(:,faceCamStrobe_idx) > faceCamStrobe_thresh;
    faceCamStrobe_up = find((~faceCamStrobe(1:end-1) & faceCamStrobe(2:end)))+1;
    faceCamStrobe_up_t = Timeline.rawDAQTimestamps(faceCamStrobe_up);
    
    % Get sync times for cameras (or load if already done)
    [facecam_sync_frames,n_facecam_frames] = AP_get_cam_sync_frames(facecam_fn);
    
    % Get the closest facecam strobe to sync start, find offset and frame idx
    [~,facecam_strobe_sync] = min(abs(camSync_up(1) - faceCamStrobe_up));
    facecam_frame_offset = facecam_sync_frames(1) - facecam_strobe_sync;
    facecam_frame_idx = [1:length(faceCamStrobe_up)] + facecam_frame_offset;
    
    % Get times of facecam frames in timeline
    facecam_t = nan(n_facecam_frames,1);
    facecam_t(facecam_frame_idx) = faceCamStrobe_up_t;
            
    save(facecam_t_savefile,'facecam_t');   
elseif exist(facecam_fn,'file') && exist(facecam_t_savefile,'file')    
    load(facecam_t_savefile);    
end

%% Load imaging data

local_dir = false;

if cam_color_n == 1
    
    if local_dir
        data_path = ['D:\data\' animal filesep day];
    else
        data_path = ['\\zserver.cortexlab.net\Data\Subjects\' animal filesep day];
    end
    
    experiment_path = [data_path filesep num2str(experiment)];
    
    disp('Loading imaging data...')
    frame_t = readNPY([experiment_path filesep 'svdTemporalComponents_blue.timestamps.npy']);
    U = readUfromNPY([data_path filesep 'svdSpatialComponents_blue.npy']);
    V = readVfromNPY([experiment_path filesep 'svdTemporalComponents_blue.npy']);
    
    disp('Done loading')
    
    framerate = 1./nanmedian(diff(frame_t));
    
    % Detrend and high-pass filter
    highpassCutoff = 0.01; % Hz
    [b100s, a100s] = butter(2, highpassCutoff/(framerate/2), 'high');
    dV = detrend(V', 'linear')';
    fV = single(filtfilt(b100s,a100s,double(dV)')');
        
    avg_im = readNPY([data_path filesep 'meanImage_blue.npy']);
    
elseif cam_color_n == 2
    
    % Load in all things as neural (n) or hemodynamic (h)
    
    if local_dir
        data_path = ['D:\data\' animal filesep day];
    else
        data_path = ['\\zserver.cortexlab.net\Data\Subjects\' animal filesep day];
    end
    
    experiment_path = [data_path filesep experiment];
    
    disp('Loading imaging data...')
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
           
    disp('Done loading')
    
    framerate = 1./nanmedian(diff(tn));

    % Correct hemodynamic signal in blue from green
    % First need to shift alternating signals to be temporally aligned
    % (shifts neural to hemo)
    % Eliminate odd frames out
    min_frames = min(size(Vn,2),size(Vh,2));
    Vn = Vn(:,1:min_frames);
    Vh = Vh(:,1:min_frames);
    
    Vn_th = SubSampleShift(Vn,1,2);
    
    Vh_Un = ChangeU(Uh,Vh,Un);
    
    %hemo_freq = [0.2,3];
    hemo_freq = [7,13];
    disp('Correcting hemodynamics...')
    Vn_hemo = HemoCorrectLocal(Un,Vn_th,Vh_Un,framerate,hemo_freq,3);
    
    disp('Filtering...')
    % Don't bother filtering heartbeat, just detrend and highpass   
    % fVn_hemo = detrendAndFilt(Vn_hemo, framerate);    
    highpassCutoff = 0.01; % Hz
    [b100s, a100s] = butter(2, highpassCutoff/(framerate/2), 'high');
    
    dVn_hemo = detrend(Vn_hemo', 'linear')';
    % wasn't zero-lag filtered before? why not?
    %fVn_hemo = filter(b100s,a100s,dVn_hemo,[],2);
    fVn_hemo = single(filtfilt(b100s,a100s,double(dVn_hemo)')');
    
    % set final U/V to use
    fV = fVn_hemo;
    U = Un;
    avg_im = avg_im_n;
    frame_t = th; % shifted to use hemo color times
        
end
disp('Done.')

% Make dF/F
[Udf,fVdf] = dffFromSVD(U,fV,avg_im);


%% Load ephys data (separated by experiment)

% This is just for a few experiments, or true with combined bank
flipped_banks = true;
load_lfp = false;

data_path = ['\\basket.cortexlab.net\data\ajpeters\' animal filesep day filesep 'ephys' filesep num2str(experiment)];

% Load clusters, if they exist
cluster_filename = [data_path filesep 'cluster_groups.csv'];
if exist(cluster_filename,'file')
    fid = fopen(cluster_filename);
    cluster_groups = textscan(fid,'%d%s','HeaderLines',1);
    fclose(fid);
end

% Load sync/photodiode
load(([data_path filesep 'sync.mat']));

% Read header information
header_path = [data_path filesep 'dat_params.txt'];
header_fid = fopen(header_path);
header_info = textscan(header_fid,'%s %s', 'delimiter',{' = '});
fclose(header_fid);

header = struct;
for i = 1:length(header_info{1})
    header.(header_info{1}{i}) = header_info{2}{i};
end

% Load spike data 
ephys_sample_rate = str2num(header.sample_rate);
spike_times = double(readNPY([data_path filesep 'spike_times.npy']))./ephys_sample_rate;
spike_templates = readNPY([data_path filesep 'spike_templates.npy']);
templates = readNPY([data_path filesep 'templates.npy']);
channel_positions = readNPY([data_path filesep 'channel_positions.npy']);
channel_map = readNPY([data_path filesep 'channel_map.npy']);
winv = readNPY([data_path filesep 'whitening_mat_inv.npy']);
template_amplitudes = readNPY([data_path filesep 'amplitudes.npy']);

% Flip channel map and positions if banks are reversed
if flipped_banks
    channel_map = [channel_map(61:end);channel_map(1:60)];
    channel_positions = [channel_positions(61:end,:);channel_positions(1:60,:)];
end

% Load LFP
if load_lfp
    n_channels = str2num(header.n_channels);
    lfp_filename = [data_path filesep 'lfp.dat'];
    fid = fopen(lfp_filename);
    lfp_all = fread(fid,[n_channels,inf],'int16');
    fclose(fid);
    % eliminate non-connected channels and sort by position
    lfp = lfp_all(channel_map+1,:);
    % get time of LFP sample points (NOTE: this is messy, based off of sample
    % rate and knowing what kwik2dat does, not sure how accurate)
    sample_rate = str2num(header.sample_rate);
    lfp_cutoff = str2num(header.lfp_cutoff);
    lfp_downsamp = (sample_rate/lfp_cutoff)/2;
    lfp_t = ([1:size(lfp,2)]*lfp_downsamp)/sample_rate;
end

% Get stim onset times
% Check that sync matches photodiode number
% Get stim onset times
% Check that sync matches photodiode number
if length(sync.timestamps) == length(photodiode.timestamps)
    refresh_rate_cutoff = 1/10;
    stim_onset_idx_ephys = [1,find(diff(sync.timestamps) > refresh_rate_cutoff) + 1];
    stim_onset_idx_timeline = [1;find(diff(photodiode.timestamps) > refresh_rate_cutoff) + 1];
else
    warning(['Ephys vs. Timeline photodiode = ' num2str(length(sync.timestamps) - length(photodiode.timestamps))]);
    
    refresh_rate_cutoff = 1/10;
    stim_onset_idx_ephys = [1,find(diff(sync.timestamps) > refresh_rate_cutoff) + 1];
    stim_onset_idx_timeline = [1;find(diff(photodiode.timestamps) > refresh_rate_cutoff) + 1];
    if length(stim_onset_idx_ephys) == length(stim_onset_idx_timeline)
        stim_onset_idx = stim_onset_idx_ephys;
        warning('But same number of stimuli');
    else
        error('Ephys vs. Timeline photodiode: different number of stimuli')
    end
end

% Get the spike/lfp times in timeline time (accounts for clock drifts)
spike_times_timeline = AP_clock_fix(spike_times,sync(1).timestamps,acqLive_timeline);
if load_lfp
    lfp_t_timeline = AP_clock_fix(lfp_t,sync(1).timestamps,acqLive_timeline);
end

% old, in case above breaks
%spike_times_timeline = AP_clock_fix(spike_times,sync.timestamps(stim_onset_idx_ephys),photodiode.timestamps(stim_onset_idx_timeline));
%lfp_t_timeline = AP_clock_fix(lfp_t,sync.timestamps(stim_onset_idx_ephys),photodiode.timestamps(stim_onset_idx_timeline));



% (by max waveform channel)
template_abs = permute(max(abs(templates),[],2),[3,1,2]);
[~,max_channel_idx] =  max(template_abs,[],1);
templateDepths = channel_positions(max_channel_idx,2);

% Get each spike's depth
spikeDepths = templateDepths(spike_templates+1);

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

disp('Done');

%% Load ephys data (single long recording)

% This is just for a few experiments where flipped or not corrected in ks
flipped_banks = false;
acqLive_channel = 2;
load_lfp = false;

disp('Loading ephys');

data_path = ['\\basket.cortexlab.net\data\ajpeters\' animal filesep day filesep 'ephys'];

% Load clusters, if they exist
cluster_filename = [data_path filesep 'cluster_groups.csv'];
if exist(cluster_filename,'file')
    fid = fopen(cluster_filename);
    cluster_groups = textscan(fid,'%d%s','HeaderLines',1);
    fclose(fid);
end

% Load sync/photodiode
load(([data_path filesep 'sync.mat']));

% Read header information
header_path = [data_path filesep 'dat_params.txt'];
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
spike_times = double(readNPY([data_path filesep 'spike_times.npy']))./ephys_sample_rate;
spike_templates = readNPY([data_path filesep 'spike_templates.npy']);
templates = readNPY([data_path filesep 'templates.npy']);
channel_positions = readNPY([data_path filesep 'channel_positions.npy']);
channel_map = readNPY([data_path filesep 'channel_map.npy']);
winv = readNPY([data_path filesep 'whitening_mat_inv.npy']);
template_amplitudes = readNPY([data_path filesep 'amplitudes.npy']);

% Flip channel map and positions if banks are reversed
if flipped_banks
    channel_map = [channel_map(61:end);channel_map(1:60)];
    channel_positions = [channel_positions(61:end,:);channel_positions(1:60,:)];
end

% Default channel map/positions are from end: make from surface
channel_positions(:,2) = max(channel_positions(:,2)) - channel_positions(:,2);

% Load LFP
n_channels = str2num(header.n_channels);
lfp_filename = [data_path filesep 'lfp.dat'];
if load_lfp && exist(lfp_filename,'file')
    fid = fopen(lfp_filename);
    lfp_all = fread(fid,[n_channels,inf],'int16');
    fclose(fid);
    % eliminate non-connected channels and sort by position (surface to deep)
    lfp = lfp_all(flipud(channel_map)+1,:);
    % get time of LFP sample points (NOTE: this is messy, based off of sample
    % rate and knowing what kwik2dat does, not sure how accurate)
    sample_rate = str2num(header.sample_rate);
    lfp_cutoff = str2num(header.lfp_cutoff);
    lfp_downsamp = (sample_rate/lfp_cutoff)/2;
    lfp_t = ([1:size(lfp,2)]*lfp_downsamp)/sample_rate;
end

% Get acqLive times for current experiment
experiment_ephys_starts = sync(acqLive_channel).timestamps(sync(acqLive_channel).values == 1);
experiment_ephys_stops = sync(acqLive_channel).timestamps(sync(acqLive_channel).values == 0);

experiments_dir = dir(fileparts(AP_cortexlab_filename(animal,day,[],'timeline')));
experiment_num = strmatch(experiment,{experiments_dir.name})-2;
acqlive_ephys_currexpt = [experiment_ephys_starts(experiment_num), ...
    experiment_ephys_stops(experiment_num)];

% Get stim onset times
% Check that sync matches photodiode number (assume 1 = acqLive,2 = pd)
% I DON'T THINK THIS IS NECESSARY
% if exist('photodiode','var')
%     
%     photodiode_currexpt = sync(2).timestamps(sync(2).timestamps > ...
%     acqlive_ephys_currextp(1) & sync(2).timestamps < ...
%     acqlive_ephys_currextp(2));
% 
%     if length(photodiode_currexpt) == length(photodiode.timestamps)
%         refresh_rate_cutoff = 1/10;
%         stim_onset_idx_ephys = [1,find(diff(photodiode_currexpt) > refresh_rate_cutoff) + 1];
%         stim_onset_idx_timeline = [1;find(diff(photodiode.timestamps) > refresh_rate_cutoff) + 1];
%     else
%         warning(['Ephys vs. Timeline photodiode = ' num2str(length(sync.timestamps) - length(photodiode.timestamps))]);
%         
%         refresh_rate_cutoff = 1/10;
%         stim_onset_idx_ephys = [1,find(diff(photodiode_currexpt) > refresh_rate_cutoff) + 1];
%         stim_onset_idx_timeline = [1;find(diff(photodiode.timestamps) > refresh_rate_cutoff) + 1];
%         if length(stim_onset_idx_ephys) == length(stim_onset_idx_timeline)
%             stim_onset_idx = stim_onset_idx_ephys;
%             warning('But same number of stimuli');
%         else
%             error('Ephys vs. Timeline photodiode: different number of stimuli')
%         end
%     end
% end

% Get the spike/lfp times in timeline time (accounts for clock drifts)
spike_times_timeline = AP_clock_fix(spike_times,acqlive_ephys_currexpt,acqLive_timeline);
if load_lfp && exist(lfp_filename,'file')
    lfp_t_timeline = AP_clock_fix(lfp_t,acqlive_ephys_currexpt,acqLive_timeline);
end

% Get the depths of each template 
% (by COM: this gives totally wonky answers because of artifacts maybe?)
%[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = ...
%    templatePositionsAmplitudes(templates,winv,channel_positions(:,2),spike_templates,template_amplitudes);

% (by max waveform channel)
template_abs = permute(max(abs(templates),[],2),[3,1,2]);
[~,max_channel_idx] =  max(template_abs,[],1);
templateDepths = channel_positions(max_channel_idx,2);

% Get each spike's depth
spikeDepths = templateDepths(spike_templates+1);

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

disp('Done');

%% Eliminate spikes that were classified as not "good"

if exist('cluster_groups','var') && ~exist('good_templates','var')
    
    disp('Removing non-good templates')
    
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
    
elseif exist('cluster_groups','var') && exist('good_templates','var')
    disp('Good templates already identified, skipping')
elseif ~exist('cluster_groups','var')
    disp('Clusters not yet sorted');
end













