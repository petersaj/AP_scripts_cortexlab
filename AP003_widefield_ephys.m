% Dual widefield/electrophysiology experiments in AP003 cortex

% This gives you these variables:
%
% U = imaging U from SVD
% fV = detrended V from SVD
% frame_t = the time of each frame relative to timeline
%
% spike_templates = the template that each spike belongs to
% spike_times_timeline = times of all spikes relative to timeline


%% Define experiment

% Pick here: 
% 2016-06-02 = M2/cingulate
% 2016-06-08 = Motor/Retrosplenial
day = '2016-06-02'; % 

% Doesn't change
animal = 'AP003';
experiment = '1';
rig = 'bigrig';
cam_color_signal = 'blue';
mpep_animal = ['M111111_' animal];


%% Load experiment info

cam_name = 'cam2';
acqLive_name = 'acqLiveEcho';

% Load timeline
timeline_filename = get_cortexlab_filename(mpep_animal,day,experiment,'timeline');
load(timeline_filename);
timeline_sample_rate = Timeline.hw.daqSampleRate;

% Get camera times
timeline_cam_idx = strcmp({Timeline.hw.inputs.name}, cam_name);
cam_samples = find(Timeline.rawDAQData(1:end-1,timeline_cam_idx) <= 2 & ...
    Timeline.rawDAQData(2:end,timeline_cam_idx) > 2);
cam_time = cam_samples./timeline_sample_rate;

% Load in hardware
hwinfo_filename = get_cortexlab_filename(mpep_animal,day,experiment,'hardware');
load(hwinfo_filename);

% Load in protocol
try
    protocol_filename = get_cortexlab_filename(mpep_animal,day,experiment,'protocol','8digit');
    load(protocol_filename);
catch me
    protocol_filename = get_cortexlab_filename(mpep_animal,day,experiment,'protocol','day_dash');
    load(protocol_filename);
end

% Get acqLive signal
acqLive_idx = strcmp({Timeline.hw.inputs.name}, acqLive_name);
thresh = max(Timeline.rawDAQData(:,acqLive_idx))/2;
acqLive_trace = Timeline.rawDAQData(:,acqLive_idx) > thresh;
acqLive_timeline = Timeline.rawDAQTimestamps( ...
    [find(acqLive_trace,1),find(acqLive_trace,1,'last')+1]);

% Get stimulus onsets and parameters
photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
thresh = max(Timeline.rawDAQData(:,photodiode_idx))/2;
photodiode_trace = Timeline.rawDAQData(:,photodiode_idx) > thresh;
photodiode_flip = find((~photodiode_trace(1:end-1) & photodiode_trace(2:end)) | ...
    (photodiode_trace(1:end-1) & ~photodiode_trace(2:end)))+1;

photodiode = struct('timestamps',[],'values',[]);
photodiode.timestamps = Timeline.rawDAQTimestamps(photodiode_flip)';
photodiode.values = photodiode_trace(photodiode_flip);

photodiode_onsets = photodiode.timestamps(photodiode.values == 1);

refresh_rate_cutoff = 1/5;
stim_onsets = photodiode_onsets( ...
    [1;find(diff(photodiode_onsets) > refresh_rate_cutoff) + 1]);

stimIDs = zeros(size(stim_onsets));
for q = 1:size(Protocol.seqnums,1)
    stimIDs(Protocol.seqnums(q,:)) = q;
end


%% Load imaging data

data_path = ['\\zserver.cortexlab.net\Data\Subjects\' mpep_animal filesep day];
experiment_path = [data_path filesep num2str(experiment)];

frame_t = readNPY([experiment_path filesep 'svdTemporalComponents_cam2.timestamps.npy']);
U = readUfromNPY([data_path filesep 'svdSpatialComponents_cam2.npy']);
V = readVfromNPY([experiment_path filesep 'svdTemporalComponents_cam2.npy']);

framerate = 1./nanmedian(diff(frame_t));

% Detrend data, remove heartbeat if that's within frequency range
if framerate > 28
    fV = detrendAndFilt(V, framerate);
else
    highpassCutoff = 0.01; % Hz
    [b100s, a100s] = butter(2, highpassCutoff/(framerate/2), 'high');
    
    dV = detrend(V', 'linear')';
    fV = filter(b100s,a100s,dV,[],2);
end

avg_im = readNPY([data_path filesep 'meanImage_cam2.npy']);


%% Load ephys data

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

% Get the spike times in timeline time
spike_times_timeline = AP_clock_fix(spike_times,sync.timestamps(stim_onset_idx_ephys),photodiode.timestamps(stim_onset_idx_timeline));


%% Eliminate spikes that were classified as not "good"

good_templates = uint32(cluster_groups{1}(strcmp(cluster_groups{2},'good')));
good_spike_idx = ismember(spike_templates,good_templates);

spike_times = spike_times(good_spike_idx);
spike_templates = spike_templates(good_spike_idx);
templates = templates(good_spike_idx);
template_amplitudes = template_amplitudes(good_spike_idx);
spike_times_timeline = spike_times_timeline(good_spike_idx);






