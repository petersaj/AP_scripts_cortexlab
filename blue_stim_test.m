%% Load stim times for blue stim test

animal = 'M111111_AP004';
day = '2016-06-21';
experiment = '11';

% Find filenames for behavior/input
timeline_filename = get_cortexlab_filename(animal,day,experiment,'timeline');

% Load timeline
load(timeline_filename);
timeline_sample_rate = Timeline.hw.daqSampleRate;

% Load in protocol
protocol_filename = get_cortexlab_filename(animal,day,experiment,'protocol');
load(protocol_filename);

% Get frames in timeline
frame_idx = strcmp({Timeline.hw.inputs.name}, 'cam2');
thresh = 1;
frame_onsets = find((Timeline.rawDAQData(1:end-1,frame_idx) <= thresh) & ...
    (Timeline.rawDAQData(2:end,frame_idx) > thresh))/timeline_sample_rate;
frame_offsets= find((Timeline.rawDAQData(1:end-1,frame_idx) > thresh) & ...
    (Timeline.rawDAQData(2:end,frame_idx) <= thresh))/timeline_sample_rate;

% (note: first frame light flip is 1 sample late, so use frame start + n)
light_frame_leeway = 3;

arduinoBlue_idx = strcmp({Timeline.hw.inputs.name}, 'arduinoBlue');
frames_b = Timeline.rawDAQData((round(frame_onsets*timeline_sample_rate)) + ...
    light_frame_leeway, ...
    arduinoBlue_idx) > thresh;

arduinoGreen_idx = strcmp({Timeline.hw.inputs.name}, 'arduinoGreen');
frames_g = Timeline.rawDAQData((round(frame_onsets*timeline_sample_rate)) + ...
    light_frame_leeway, ...
    arduinoGreen_idx) > thresh;

frame_onsets_b = frame_onsets(frames_b);
frame_offsets_b = frame_offsets(frames_b);

frame_onsets_g = frame_onsets(frames_g);
frame_onsets_g = frame_offsets(frames_g);

% Get stimulus onsets and parameters
photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
% (the only way to find stim onsets is to look for gaps in photodiode,
% which seems messy, but for now just look for gaps that are > the refresh
% rate of the monitor, which looks to be about 30 Hz)
% (because time between stimuli not recorded?!?!)
photodiode_onsets = find((Timeline.rawDAQData(1:end-1,photodiode_idx) <= thresh) & ...
    (Timeline.rawDAQData(2:end,photodiode_idx) > thresh))/timeline_sample_rate;

refresh_rate_cutoff = 1/10;
stim_onsets = photodiode_onsets( ...
    [1;find(diff(photodiode_onsets) > refresh_rate_cutoff) + 1]);

stimIDs = zeros(size(stim_onsets));
for q = 1:size(Protocol.seqnums,1)
    stimIDs(Protocol.seqnums(q,:)) = q;
end



%% Load in SVD data
svd_path = 'C:\Users\Andrew\Documents\CarandiniHarrisLab\data\bluestim_noshield_test';

nSV = 2000;
Fs = 35;

U = readUfromNPY([svd_path filesep 'U.npy'], nSV);
V = readVfromNPY([svd_path filesep 'V.npy'], nSV);
fV = detrendAndFilt(V, Fs);
t = readNPY([svd_path filesep 't.npy']);

% Get frames corresponding to time around stim (hardcoded now, but can get
% stim time)

framerate = round(1./nanmedian(diff(frame_onsets_b)));
frames_back = framerate*1;
frames_forward = framerate*1;
total_frames = frames_back + frames_forward + 1;

[~,stim_frames] = arrayfun(@(x) min(abs(frame_onsets_b - ...
    stim_onsets(x))),1:length(stim_onsets));

frames_pull = repmat(stim_frames,total_frames,1) + ...
    repmat(transpose(-frames_back:frames_forward),1,length(stim_frames));

fV_stim = permute(cell2mat(arrayfun(@(x) reshape(fV(x,frames_pull),size(frames_pull,1), ...
    size(frames_pull,2)),permute(1:size(V,1),[1,3,2]),'uni',false)),[1,3,2]);

% Average activity for each type of stimulus
n_stim = Protocol.npfilestimuli;
n_repeats = Protocol.nrepeats;
stim_sequence = Protocol.seqnums;
fV_stim_grp = nan(size(fV_stim,1),size(fV_stim,2),n_stim);
for curr_stim = 1:n_stim
    fV_stim_grp(:,:,curr_stim) = ...
        nanmean(fV_stim(:,:,stim_sequence(curr_stim,:)),3);
end

% Recover spatial activity from SVs
a = permute(reshape(fV_stim_grp(:,:,1)*reshape(U,[],size(U,3))', ...
    total_frames,size(U,1),size(U,2)),[2,3,1]);



%% Widefield GUI tests

pixelTuningCurveViewerSVD(U,fV,t,stim_onsets,stimIDs,[-1 5])

pixelCorrelationViewerSVD(U,V);

movieWithTracesSVD(U,fV,t,[],[]);


%% Load in new pipeline named SVD data
data_path = ['\\zserver.cortexlab.net\Data\Subjects\' mpep_animal filesep day];
experiment_path = [data_path filesep experiment];
Fs = 35;

Ub = readUfromNPY([data_path filesep 'svdSpatialComponents_blue.npy']);
Vb = readVfromNPY([experiment_path filesep 'svdTemporalComponents_blue.npy']);
fVb = detrendAndFilt(Vb, Fs);
tb = readNPY([experiment_path filesep 'svdTemporalComponents_blue.timestamps.npy']);

Ug = readUfromNPY([data_path filesep 'svdSpatialComponents_green.npy']);
Vg = readVfromNPY([experiment_path filesep 'svdTemporalComponents_green.npy']);
fVg = detrendAndFilt(Vg, Fs);
tg = readNPY([experiment_path filesep 'svdTemporalComponents_green.timestamps.npy']);

% Correct hemodynamic signal in blue from green
% First need to shift alternating signals to be temporally aligned
% NOTE! this assumes blue comes first (which is probably the new standard)
% Eliminate odd frames out
min_frames = min(size(fVb,2),size(fVg,2));
fVb = fVb(:,1:min_frames);
fVg = fVg(:,1:min_frames);

fVb_tg = SubSampleShift(fVb,1,2);

fVg_Ub = ChangeU(Ug,fVg,Ub);

hemo_freq = [0.2,3];
fVbh = HemoCorrectLocal(Ub,fVb_tg,fVg_Ub,35,hemo_freq,3);

pixelTuningCurveViewerSVD(Ub,fVbh,tg(1:min_frames),stim_onsets,stimIDs,[-1 5])




















