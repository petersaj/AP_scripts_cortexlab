%% No illumination test: GCaMP animal with gratings


%% Load experiment info

animal = 'M160302_NS1HEY';
day = '2016-01-01';
experiment = '2';

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


% Load data

data_path = ['\\zserver.cortexlab.net\Data\Subjects\' animal filesep day];
experiment_path = [data_path filesep experiment];

Fs = 35;

U = readUfromNPY([data_path filesep 'svdSpatialComponents_all.npy']);
V = readVfromNPY([experiment_path filesep 'svdTemporalComponents_all.npy']);
t = readNPY([experiment_path filesep 'svdTemporalComponents_all.timestamps.npy']);

avg_im = readNPY([data_path filesep 'meanImage_all.npy']);

%pixelTuningCurveViewerSVD(U,V,t,stim_onsets,stimIDs,[-1 3])


% Get average image during stim
stim_time = 2;
baseline_time = 0.5;
avg_stim_im = nan(size(U,1),size(U,2),max(stimIDs));
avg_baseline_im = nan(size(U,1),size(U,2),max(stimIDs));
for curr_stim = unique(stimIDs)
   curr_onsets = stim_onsets(stimIDs == curr_stim); 
   
   stim_t = arrayfun(@(x) find(t >= curr_onsets(x) & t <= curr_onsets(x) + ...
       stim_time),1:length(curr_onsets),'uni',false);
   baseline_t = arrayfun(@(x) find(t >= curr_onsets(x)-baseline_time & t < ...
       curr_onsets(x)),1:length(curr_onsets),'uni',false);
   
   avg_stim_im(:,:,curr_stim) = ...
       nanmean(svdFrameReconstruct(U,V(:,horzcat(stim_t{:}))),3);    
   avg_baseline_im(:,:,curr_stim) = ...
       nanmean(svdFrameReconstruct(U,V(:,horzcat(baseline_t{:}))),3);  
end

avg_stim_im_baselinesub = avg_stim_im - avg_baseline_im;

figure;
for curr_stim = 1:12;
   subplot(2,6,curr_stim);
   imagesc(avg_stim_im_baselinesub(:,:,curr_stim));
   colormap(gray);
   axis off
   title(['Stim ' num2str(curr_stim)]);
end








