%% Hemodynamic test: WT animal with gratings, green illumination only


%% Load experiment info

animal = 'M111111_AP004';
day = '2016-06-21';
experiment = '11';

% Load timeline
timeline_filename = get_cortexlab_filename(animal,day,experiment,'timeline');
load(timeline_filename);
timeline_sample_rate = Timeline.hw.daqSampleRate;

% Load in protocol
protocol_filename = get_cortexlab_filename(animal,day,experiment,'protocol');
load(protocol_filename);

% Get stimulus onsets and parameters

photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
thresh = max(Timeline.rawDAQData(:,photodiode_idx))/2;
photodiode_onsets = Timeline.rawDAQTimestamps((Timeline.rawDAQData(1:end-1,photodiode_idx) <= thresh) & ...
    (Timeline.rawDAQData(2:end,photodiode_idx) > thresh));

refresh_rate_cutoff = 1/10;
stim_onsets = photodiode_onsets( ...
    [1,find(diff(photodiode_onsets) > refresh_rate_cutoff) + 1]);

stimIDs = zeros(size(stim_onsets));
for q = 1:size(Protocol.seqnums,1)
    stimIDs(Protocol.seqnums(q,:)) = q;
end




%% Load data

data_path = ['\\zserver.cortexlab.net\Data\Subjects\' animal filesep day];
experiment_path = [data_path filesep experiment];

Fs = 35;

U = readUfromNPY([data_path filesep 'svdSpatialComponents_cam2.npy']);
V = readVfromNPY([experiment_path filesep 'svdTemporalComponents_cam2.npy']);
fV = detrendAndFilt(V, Fs);
t = readNPY([experiment_path filesep 'svdTemporalComponents_cam2.timestamps.npy']);

avg_im = readNPY([data_path filesep 'meanImage_cam2.npy']);

pixelTuningCurveViewerSVD(U,fV,t,stim_onsets,stimIDs,[-1 3])















