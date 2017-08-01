%% Preprocess and kilosort data (IMEC Phase 3)

error('Use AP_preprocess_phase3 instead')

animal = 'AP018';
day = '2017-05-23';

data_path =  ...
    ['\\zserver.cortexlab.net\Data\Subjects\' animal filesep day '\ephys'];
save_path =  ...
    ['\\basket.cortexlab.net\data\ajpeters\' animal filesep day '\ephys'];

if ~exist(save_path,'dir')
    mkdir(save_path)
end

% Hardcode the filenames here... are these always like this?
ap_data_filename = [data_path filesep 'experiment1_100-0_0.dat'];
lfp_data_filename = [data_path filesep 'experiment1_100-1_0.dat'];
sync_filename = [data_path filesep 'experiment1_all_channels_0.events'];
messages_filename = [data_path filesep 'experiment1_messages_0.events'];
settings_filename = [data_path filesep 'settings.xml'];

%%%%%%%%%%%%%%%%%%%%%%%%% MAKE THIS WORK WHEN DATA FOR REAL
% Get index of electrophysiology channels in recordings
ephys_settings = xml2struct(settings_filename);

% Get sample rate, gain, cutoff frequency

apGain = ephys_settings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.EDITOR.NEUROPIXELS.Attributes.apGainValue;
lfpGain = ephys_settings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.EDITOR.NEUROPIXELS.Attributes.lfpGainValue;
filterCut = ephys_settings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.EDITOR.NEUROPIXELS.Attributes.filterCut;

%%%%% MAKE THESE PARAMETERS THE RELEVANT ONES
% Nick must've done this already...
% relationship between gain and 0.195x for int16 to uV?

params = {'raw_path',['''' data_path '''']; ...
    'n_channels',num2str(384); ...
    'sample_rate',num2str(sample_rate); ... % this should be 30000 AP, 2500 LFP
    'gain',num2str(ch_gain); ...
    'lfp_cutoff',num2str(lfp_cutoff)};

param_filename = [save_path filesep 'dat_params.txt'];

formatSpec = '%s = %s\r\n';
fid = fopen(param_filename,'w');
for curr_param = 1:size(params,1)
    fprintf(fid,formatSpec,params{curr_param,:});
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%

% Pull out and save events (saved in Open Ephys format)
% Need to fix timing
% From the thread about it:
% 2) It seems that event times require this "get_experiment_start_time" correction, but I cannot find where that start time is recorded. Do you know where it is?
% 
% This is in the messages_0.events files. In python I read it in with:
% 
% int(open(os.path.join(self.datapath,'experiment1_messages_0.events'),'r').readlines()[0].split(' ')[0])
error('FIX THE TIMING STUFF?')

% Get experiment start time (these messages are saved in a super dumb way
% that are impossible to parse generally, so this is messy)
messages_id = fopen(messages_filename);
messages_text = textscan(messages_id,'%*d %s %s', 'delimiter',{': '});
fclose(messages_id);

start_time_idx = strcmp(messages_text{1},'start time');
start_time = str2num(messages_text{2}{start_time_idx}(1:strfind(messages_text{2}{start_time_idx},'@')-1));
start_time_freq = str2num(messages_text{2}{start_time_idx}(strfind(messages_text{2}{start_time_idx},'@')+1: ...
    strfind(messages_text{2}{start_time_idx},'Hz')-1));
start_time_sec = start_time/start_time_freq;
% TO DO: make sure this is real (the number looks potentially resonable,
% for the test it gave a 36 second delay between ephys recording and
% acqLive and a 5 minute delay to hit record...)

[sync_data, sync_timestamps, sync_info] = load_open_ephys_data_faster(sync_filename);
sync_channels = unique(sync_data);
sync = struct('timestamps',cell(size(sync_channels)),'values',cell(size(sync_channels)));
for curr_sync = 1:length(sync_channels)
    sync_events = sync_data == (sync_channels(curr_sync));
    sync(curr_sync).timestamps = sync_timestamps(sync_events);
    sync(curr_sync).values = logical(sync_info.eventId(sync_events));
end
sync_save_filename = [save_path filesep 'sync.mat'];
save(sync_save_filename,'sync');

% Copy files to local drive to speed up loading
disp('Copying data to local drive...')
temp_path = 'C:\data_temp\kilosort';
ap_temp_filename = [temp_path filesep animal '_' day  '_' 'ephys_apband.dat'];
if ~exist(temp_path,'dir')
    mkdir(temp_path)
end
copyfile(ap_data_filename,ap_temp_filename);
disp('Done');

% Subtract common median across AP-band channels (hardcode channels?)
ops.NchanTOT = 384;
medianTrace = applyCARtoDat(ap_temp_filename, ops.NchanTOT);
ap_temp_car_filename = [ap_temp_filename(1:end-4) '_CAR.dat'];

% Get rid of the original non-CAR (usually not enough disk space)
delete(ap_temp_filename);

% Run kilosort on CAR data
sample_rate = 30000; % in the future, get this from somewhere?
AP_run_kilosort(ap_temp_car_filename,sample_rate)







%% Convert and kilosort data (IMEC Phase 2)

animal = 'AP009';
day = '2016-11-23';
combined_bank = true;
subtract_light = true;

kwik_path =  ...
    ['\\zserver.cortexlab.net\Data\Subjects\' animal filesep day '\ephys'];

save_path =  ...
    ['\\basket.cortexlab.net\data\ajpeters\' animal filesep day '\ephys'];

% Convert from Kwik format to raw
sync_channel = [1,2,3,4];
sync_input = 'adc';
local_copy = true;
kwik2dat(kwik_path,save_path,sync_channel,sync_input,local_copy,subtract_light);

% Run Kilosort
% (get header info to get sample rate)
header_path = [save_path filesep 'dat_params.txt'];
header_fid = fopen(header_path);
header_info = textscan(header_fid,'%s %s', 'delimiter',{' = '});
fclose(header_fid);

header = struct;
for i = 1:length(header_info{1})
    header.(header_info{1}{i}) = header_info{2}{i};
end

ephys_sample_rate = str2num(header.sample_rate);

input_board = 'oe';
data_filename = [save_path filesep 'spikes.dat'];
AP_run_kilosort_old(data_filename,input_board,ephys_sample_rate,combined_bank)


%% Cut data and do kilosort (if bad section of recording)
% I know this is probably a dumb way to do it: but at the moment, just load
% the data, cut it, save a cut version, then do kilosort on that

animal = 'AP009';
day = '2016-11-03';
combined_bank = true;


data_path = ['\\basket.cortexlab.net\data\ajpeters\' animal filesep day filesep 'ephys'];

% Read header information
header_path = [data_path filesep 'dat_params.txt'];
header_fid = fopen(header_path);
header_info = textscan(header_fid,'%s %s', 'delimiter',{' = '});
fclose(header_fid);

header = struct;
for i = 1:length(header_info{1})
    header.(header_info{1}{i}) = header_info{2}{i};
end

% Load LFP
disp('Loading LFP...')
n_channels = str2num(header.n_channels);
lfp_filename = [data_path filesep 'lfp.dat'];
fid = fopen(lfp_filename);
lfp_all = fread(fid,[n_channels,inf],'int16');
fclose(fid);

lfp_fig = figure;
imagesc(lfp_all);
colormap(gray)
title('Select last timepoint to use');
[end_time,~] = ginput(1);

close(lfp_fig);

sample_rate = str2num(header.sample_rate);
lfp_cutoff = str2num(header.lfp_cutoff);
lfp_downsamp = (sample_rate/lfp_cutoff)/2;

end_sample = end_time*lfp_downsamp;

clear lfp_all;

% Bring spikes to local
disp('Copying data to local SSD...');
save_path =  ...
    ['\\basket.cortexlab.net\data\ajpeters\' animal filesep day '\ephys'];
data_filename = [save_path filesep 'spikes.dat'];
[data_path,data_file,data_ext] = fileparts(data_filename);
local_path = 'C:\Users\Andrew\Documents\CarandiniHarrisLab\data\kilosort_temp';
local_file = [data_file '_cut' data_ext];
local_data_filename = [local_path filesep local_file];
if ~exist(local_data_filename)
    copyfile(data_filename,local_data_filename);
else
    disp('Already copied');
end
disp('Done');

% Load and save spikes, cutting off everything past selected point
disp('Loading data until cut point, re-saving (locally)...')
fid = fopen(local_data_filename, 'r');
spikes = fread(fid,[str2num(header.n_channels),end_sample],'*int16');
fclose(fid);

spikes_fid = fopen(local_data_filename, 'w');
fwrite(spikes_fid,spikes,'int16');
fclose(spikes_fid);

clear spikes;

% Run Kilosort
title('Running kilosort...')
data_filename = [save_path filesep 'spikes_cut.dat'];
ephys_sample_rate = str2num(header.sample_rate);
input_board = 'oe';
AP_run_kilosort(data_filename,input_board,ephys_sample_rate,combined_bank);

disp('Done')

%% Load in data from one day

animal = '65';
day = '20151102';

processed_data_path = '\\basket.cortexlab.net\data\ajpeters';
raw_data_filename = ['\\zserver.cortexlab.net\Data\multichanspikes\' ...
    animal filesep day filesep day '_1.dat'];

% Get ephys sample rate from header
ephys_header_filename = [raw_data_filename(1:end-4) '.meta'];
ephys_header_id = fopen(ephys_header_filename);
ephys_header = textscan(ephys_header_id,'%s %s', 'delimiter',{' ='});
fclose(ephys_header_id);

sample_rate_idx = strcmp(ephys_header{1},'sRateHz');
ephys_sample_rate = str2num(ephys_header{2}{sample_rate_idx});

% Find filenames for behavior/input
day_dash = datestr(datenum(day,'yyyymmdd'),'yyyy-mm-dd');
timeline_filename = get_cortexlab_filename(animal,day_dash,1,'timeline');
parameters_filename = get_cortexlab_filename(animal,day_dash,1,'parameters');
block_filename = get_cortexlab_filename(animal,day_dash,1,'block');

% Load behavior/input
load(timeline_filename);
load(parameters_filename);
load(block_filename);

% Load ephys
ephys_path = [processed_data_path filesep animal filesep day];
spike_times = readNPY([ephys_path filesep 'spike_times.npy']);
spike_templates = readNPY([ephys_path filesep 'spike_templates.npy']);
templates = readNPY([ephys_path filesep 'templates.npy']);
channel_positions = readNPY([ephys_path filesep 'channel_positions.npy']);

spike_times_sec = double(spike_times) / ephys_sample_rate;

sorted = exist([ephys_path filesep 'spike_clusters.npy'],'file');
if sorted
    
    % Load clusters if sorted in Phy
    spike_clusters = readNPY([ephys_path filesep 'spike_clusters.npy']);
    
    % Load cluster labels
    cluster_groups_filename = [ephys_path filesep 'cluster_groups.csv'];    
    cluster_groups_id = fopen(cluster_groups_filename);
    cluster_groups_read = textscan(cluster_groups_id,'%f %s', 'delimiter','\t','HeaderLines',1);   
    
    % 1-index clusters
    cluster_groups_1idx = cluster_groups_read{1} + 1;
    
    % Make cluster labels in indicies
    cluster_groups = struct( ...
        'good',false(max(cluster_groups_1idx),1), ...
        'mua',false(max(cluster_groups_1idx),1), ...
        'noise',false(max(cluster_groups_1idx),1), ...
        'unsorted',false(max(cluster_groups_1idx),1));
    
    cluster_groups.good(cluster_groups_1idx( ...
        strcmp(cluster_groups_read{2},'good'))) = true;
    cluster_groups.mua(cluster_groups_1idx( ...
        strcmp(cluster_groups_read{2},'mua'))) = true;
    cluster_groups.noise(cluster_groups_1idx( ...
        strcmp(cluster_groups_read{2},'noise'))) = true;
    cluster_groups.unsorted(cluster_groups_1idx( ...
        strcmp(cluster_groups_read{2},'unsorted'))) = true;
else
    
    % If no clusters, use template indicies
    spike_clusters = spike_templates;
    
    % Make all clusters good
    % Make cluster labels in indicies
    cluster_groups = struct( ...
        'good',true(max(spike_clusters),1), ...
        'mua',false(max(spike_clusters),1), ...
        'noise',false(max(spike_clusters),1), ...
        'unsorted',false(max(spike_clusters),1));
    
end

% Templates/clusters are 0-indexed, fix that
spike_templates_1idx = spike_templates + 1;
spike_clusters_1idx = spike_clusters + 1;

% FOR NOW: simplify spike times by grouping only used clusters
cluster_spike_times = cell(max(spike_clusters_1idx),1);
for curr_cluster = 1:max(spike_clusters_1idx)'
    cluster_spike_times{curr_cluster} = ...
        double(spike_times(spike_clusters_1idx == curr_cluster))/ephys_sample_rate;
end

% Get sync signal from ephys
ephys_sync_samples = AP_sync_from_channel(raw_data_filename,processed_data_path);
% four pulses: only first two are valid??
ephys_sync_samples = ephys_sync_samples(1:2);

% Sync signal from timeline 
timeline_sync_idx = strcmp({Timeline.hw.inputs.name}, 'camSync');
timeline_sync_samples = find(Timeline.rawDAQData(1:end-1,timeline_sync_idx) <= 2 & ...
    Timeline.rawDAQData(2:end,timeline_sync_idx) > 2);
timeline_sample_rate = Timeline.hw.daqSampleRate;

% Get conversion from timeline to ephys
timeline_ephys_scale = diff(ephys_sync_samples)/diff(timeline_sync_samples);
timeline_ephys_offset = ephys_sync_samples(1) - timeline_sync_samples(1)*timeline_ephys_scale;

% Photodiode signal from timeline
photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
photodiode = Timeline.rawDAQData(:,photodiode_idx);
photodiode_flip = find((photodiode(1:end-1) <= 1.5 & photodiode(2:end) > 1.5) | ...
    (photodiode(1:end-1) >= 1.5 & photodiode(2:end) < 1.5));
photodiode_sec = photodiode_flip/timeline_sample_rate;

% Get stimulus presentation times relative to ephys
n_trials = length(block.trial);
stim_times_block = [block.trial.stimulusCueStartedTime];
stim_samples_photodiode = photodiode_flip(arrayfun(@(x) find(photodiode_sec > ...
    stim_times_block(x),1),1:n_trials));
stim_samples_ephys = stim_samples_photodiode*timeline_ephys_scale + timeline_ephys_offset;
stim_times_ephys = stim_samples_ephys/ephys_sample_rate;

% Get stimulus conditions
trial_condition = cell2mat(cellfun(@(x) x.visCueContrast,{block.trial(:).condition},'uni',false));

% Get most used template for each cluster
cluster_template = nan(max(spike_clusters_1idx),1);
for curr_clust = unique(spike_clusters_1idx)';
    cluster_template(curr_clust) = ...
        mode(spike_templates_1idx(spike_clusters_1idx == curr_clust));
end
used_clusters = ~isnan(cluster_template);

% Location (center of mass) of each template
template_abs = permute(max(abs(templates),[],2),[3,1,2]);

template_com_x = (template_abs'*channel_positions(:,1))./sum(template_abs,1)';
template_com_y = (template_abs'*channel_positions(:,2))./sum(template_abs,1)';

% Location (center of mass) of each cluster
cluster_com_x = nan(max(spike_clusters),1);
cluster_com_y = nan(max(spike_clusters),1);

cluster_com_x(used_clusters) = template_com_x(cluster_template(used_clusters));
cluster_com_y(used_clusters) = template_com_y(cluster_template(used_clusters));

% Get the largest amplitude waveform for each template
[~,max_site] = max(max(abs(templates),[],2),[],3);
templates_max = nan(size(templates,1),size(templates,2));
for curr_template = 1:size(templates,1)
    templates_max(curr_template,:) = ...
        templates(curr_template,:,max_site(curr_template));
end
used_templates = any(templates_max,2);

% Assign waveform to each cluster
cluster_waveform = nan(max(spike_clusters),size(templates_max,2));
cluster_waveform(used_clusters,:) = templates_max(cluster_template(used_clusters),:);

% Get FWHM of each waveform
waveform_halfmax = max(abs(templates_max),[],2)./2;

good_waveforms = false(size(templates_max,1),1);
waveform_halfmax_start = nan(size(templates_max,1),1);
waveform_halfmax_stop = nan(size(templates_max,1),1);

waveform_halfmax_start(used_templates) = arrayfun(@(x) find(abs(templates_max(x,:)) > ...
    waveform_halfmax(x),1),find(used_templates));

good_waveforms(used_templates) = cellfun(@(x) ~isempty(x),arrayfun(@(x) ...
    find(abs(templates_max(x,waveform_halfmax_start(x):end)) < ...
    waveform_halfmax(x),1),find(used_templates),'uni',false));

waveform_halfmax_stop(used_templates & good_waveforms) = ...
    arrayfun(@(x) ...
    find(abs(templates_max(x,waveform_halfmax_start(x):end)) < ...
    waveform_halfmax(x),1),find(used_templates & good_waveforms)) + ...
    waveform_halfmax_start(used_templates & good_waveforms);

template_fwhm = waveform_halfmax_stop - waveform_halfmax_start;

cluster_fwhm = nan(max(spike_clusters),1);
cluster_fwhm(used_clusters) = template_fwhm(cluster_template(used_clusters));


% Store relevant spike data in master structure
spikes = struct;

spikes.sorted = sorted;

spikes.sample_rate = ephys_sample_rate;
spikes.spike_samples = spike_times;
spikes.spike_times = spike_times_sec;
spikes.spike_templates = spike_templates_1idx;
spikes.spike_clusters = spike_clusters_1idx;
spikes.channel_positions = channel_positions;

spikes.template.waveform = templates_max;
spikes.template.location = [template_com_x,template_com_y];
spikes.template.fwhm = template_fwhm; 

spikes.cluster.spike_times = cluster_spike_times;
spikes.cluster.groups = cluster_groups;
spikes.cluster.template = cluster_template;
spikes.cluster.waveform = cluster_waveform;
spikes.cluster.location = [cluster_com_x,cluster_com_y];
spikes.cluster.fwhm = cluster_fwhm;

% Store experiment data in master structure
xpr = struct;

xpr.animal = animal;
xpr.day = day;
xpr.block = block;
xpr.timeline = Timeline;

xpr.timeline.sample_rate = timeline_sample_rate;
xpr.timeline.photodiode_time = photodiode_sec;

xpr.bhv.condition = trial_condition;
xpr.bhv.stim_onset_ephys = stim_times_ephys;

xpr.timeline2ephys = [timeline_ephys_scale,timeline_ephys_offset];


%% PSTH stuff

trial_condition_rl = (trial_condition(1,:) > trial_condition(2,:)) - ...
    (trial_condition(2,:) > trial_condition(1,:));

% Manually scroll through PSTH
spike_times_sec = double(spike_times) / ephys_sample_rate;
raster_window = [-2,2];
psthViewer(spike_times_sec,spike_clusters_1idx, ...
    stim_times_ephys,raster_window,trial_condition_rl);

% Get PSTH for all clusters
raster_window = [-2,2];
psth_bin_size = 0.001;
n_bins = diff(raster_window)/psth_bin_size + 1;

cluster_psth = nan(length(cluster_spike_times),n_bins);
for curr_cluster = unique(spike_clusters_1idx)';
    
    [psth,bins,rasterX,rasterY,spikeCounts] = ...
        psthRasterAndCounts(cluster_spike_times{curr_cluster}, ...
        stim_times_ephys, raster_window, psth_bin_size);
    
    cluster_psth(curr_cluster,:) = psth;
    
end

% The first and last bins are weird, eliminate
cluster_psth = cluster_psth(:,2:end-1);
bins = bins(2:end-1);

% Gaussian smooth psth
smooth_size = 50;
gw = gausswin(round(smooth_size*6),3)';
smWin = gw./sum(gw);
cluster_psth_smooth = conv2(cluster_psth, smWin, 'same');

cluster_psth_smooth_zscore = zscore(cluster_psth_smooth,[],2);

figure;plot(bins,nanmean(cluster_psth_smooth_zscore),'k','linewidth',2);



%% Get average normalized PSTH for all days

animal = '65';

days = {'20151028','20151029','20151030','20151031','20151101','20151102','20151103'};
processed_data_path = '\\basket.cortexlab.net\data\ajpeters';

 raster_window = [-2,2];
 psth_bin_size = 0.001;
 n_bins = (diff(raster_window)/psth_bin_size + 1) - 2;
 smoothed_norm_mean_psth = nan(length(days),n_bins);

for curr_day = 1:length(days)
    
    day = days{curr_day};
    
    raw_data_filename = ['\\zserver.cortexlab.net\Data\multichanspikes\' ...
        animal filesep day filesep day '_1.dat'];
    
    % Get ephys sample rate from header
    ephys_header_filename = [raw_data_filename(1:end-4) '.meta'];
    ephys_header_id = fopen(ephys_header_filename);
    ephys_header = textscan(ephys_header_id,'%s %s', 'delimiter',{' ='});
    fclose(ephys_header_id);
    
    sample_rate_idx = strcmp(ephys_header{1},'sRateHz');
    ephys_sample_rate = str2num(ephys_header{2}{sample_rate_idx});
    
    % Find filenames for behavior/input
    day_dash = datestr(datenum(day,'yyyymmdd'),'yyyy-mm-dd');
    timeline_filename = get_cortexlab_filename(animal,day_dash,1,'timeline');
    parameters_filename = get_cortexlab_filename(animal,day_dash,1,'parameters');
    block_filename = get_cortexlab_filename(animal,day_dash,1,'block');
    
    % Load behavior/input
    load(timeline_filename);
    load(parameters_filename);
    load(block_filename);
    
    % Load ephys
    ephys_path = ['\\basket.cortexlab.net\data\ajpeters\' animal filesep day];
    spike_times = readNPY([ephys_path filesep 'spike_times.npy']);
    spike_templates = readNPY([ephys_path filesep 'spike_templates.npy']);
    templates = readNPY([ephys_path filesep 'templates.npy']);
    
    sorted = exist([ephys_path filesep 'spike_clusters.npy'],'file');
    % Load clusters if manually sorted already
    if sorted
        spike_clusters = readNPY([ephys_path filesep 'spike_clusters.npy']);
        % If not, use template indicies as clusters
    else
        spike_clusters = spike_templates;
    end
    
    % Clusters are 0-indexed, fix that
    spike_clusters_1idx = spike_clusters + 1;
    
    % FOR NOW: simplify spike times by grouping only used clusters
    cluster_spike_times = cell(max(spike_clusters_1idx),1);
    for curr_cluster = 1:max(spike_clusters_1idx)'
        cluster_spike_times{curr_cluster} = ...
            double(spike_times(spike_clusters_1idx == curr_cluster))/ephys_sample_rate;
    end
    
    % Get sync signal from ephys
    ephys_sync_samples = AP_sync_from_channel(raw_data_filename,processed_data_path);
    % four pulses: only first two are valid??
    ephys_sync_samples(3:4) = [];
    
    % Sync signal from timeline
    timeline_sync_idx = strcmp({Timeline.hw.inputs.name}, 'camSync');
    timeline_sync_samples = find(Timeline.rawDAQData(1:end-1,timeline_sync_idx) <= 2 & ...
        Timeline.rawDAQData(2:end,timeline_sync_idx) > 2);
    timeline_sample_rate = Timeline.hw.daqSampleRate;
    
    % Get conversion from timeline to ephys
    timeline_ephys_scale = diff(ephys_sync_samples)/diff(timeline_sync_samples);
    timeline_ephys_offset = ephys_sync_samples(1) - timeline_sync_samples(1)*timeline_ephys_scale;
    
    % Photodiode signal from timeline
    photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
    photodiode = Timeline.rawDAQData(:,photodiode_idx);
    photodiode_flip = find((photodiode(1:end-1) <= 1.5 & photodiode(2:end) > 1.5) | ...
        (photodiode(1:end-1) >= 1.5 & photodiode(2:end) < 1.5))+1;
    photodiode_sec = photodiode_flip/timeline_sample_rate;
    
    % Get stimulus presentation times relative to ephys
    n_trials = length(block.trial);
    stim_times_block = [block.trial.stimulusCueStartedTime];
    stim_samples_photodiode = photodiode_flip(arrayfun(@(x) find(photodiode_sec > ...
        stim_times_block(x),1),1:n_trials));
    stim_samples_ephys = stim_samples_photodiode*timeline_ephys_scale + timeline_ephys_offset;
    stim_times_ephys = stim_samples_ephys/ephys_sample_rate;
    
    
    % Get PSTH for all clusters
    raster_window = [-2,2];
    psth_bin_size = 0.001;
    n_bins = diff(raster_window)/psth_bin_size + 1;
    
    cluster_psth = nan(length(cluster_spike_times),n_bins);
    for curr_cluster = unique(spike_clusters_1idx)';
        
        [psth,bins,rasterX,rasterY,spikeCounts] = ...
            psthRasterAndCounts(cluster_spike_times{curr_cluster}, ...
            stim_times_ephys, raster_window, psth_bin_size);
        
        cluster_psth(curr_cluster,:) = psth;
        
    end
    
    % The first and last bins are weird, eliminate
    cluster_psth = cluster_psth(:,2:end-1);
    bins = bins(2:end-1);
    
    % Gaussian smooth psth
    smooth_size = 15;
    gw = gausswin(round(smooth_size*6),3)';
    smWin = gw./sum(gw);
    cluster_psth_smooth = conv2(cluster_psth, smWin, 'same');
    
    cluster_psth_smooth_zscore = zscore(cluster_psth_smooth,[],2);
        
    smoothed_norm_mean_psth(curr_day,:) = nanmean(cluster_psth_smooth_zscore,1);
    
    disp(['Done ' day]);
    
end

plot_space = 2;
smoothed_norm_mean_psth_spaced = bsxfun(@plus, ...
    smoothed_norm_mean_psth,plot_space*transpose(1:length(days)));
figure;
plot(bins,smoothed_norm_mean_psth_spaced','k','linewidth',2);
line([0,0],ylim,'color','r','linewidth',1,'linestyle','--');

xlabel('Time from stimulus onset (s)');
ylabel('Mean normalized PSTH');
set(gca,'YTick',plot_space:plot_space:plot_space*length(days));
set(gca,'YTickLabel',{'V1 (loc 1)','RL','AL','V1 (loc 2)', 'PM','V1 (loc 1,RF match)','Cg'});



%% Load all data !!! USING STRUCTURES FROM AP_load_ephys !!!

animal = 'AP002';
days = {'2017-07-13'};

curr_day = 1;
[spikes,xpr] = AP_load_ephys(animal,days{curr_day});


%% Tuning curve for spikes

% Look at number of templates at each depth to estimate cortex boundaries
h = figure('Position',[94,122,230,820]);
plotSpread(spikes.template.depth,'distributionColors','k')
axis off
title('Mark cortex top')
[~,cortex_top] = ginput(1);
title('Mark cortex bottom');
[~,cortex_bottom] = ginput(1);
close(h);
drawnow;
cortex_depth = cortex_top - cortex_bottom;

% Make depths relative to cortex top
cluster_depth_cortex = cortex_top - spikes.cluster.depth;

% Group spikes by evenly spaced location in cortex
% (0 = out of cortex, n_depth_groups+1 = subcortex)
n_depth_groups = 4;
depth_group_edges = [linspace(0,cortex_depth,n_depth_groups+1),Inf];
[depth_group_n,depth_group] = histc(cluster_depth_cortex,depth_group_edges);
depth_groups_used = unique(depth_group);
depth_group_centers = grpstats(cluster_depth_cortex,depth_group);

% Get PSTH for all clusters
raster_window = [-2,2];
spike_count_window = [0,2];
psth_bin_size = 0.001;
n_bins = diff(raster_window)/psth_bin_size + 1;

n_clusters = length(spikes.cluster.spike_times);
n_trials = size(xpr.bhv.condition,1);

cluster_spike_count = nan(n_clusters,n_trials);
cluster_spike_count_correct = nan(n_clusters,sum(xpr.bhv.correct));
cluster_spike_count_incorrect = nan(n_clusters,sum(~xpr.bhv.correct));
for curr_cluster = 1:n_clusters;
    
    % Spike counts for all trials
    [~,~,~,~,cluster_spike_count(curr_cluster,:)] = ...
        psthRasterAndCounts(spikes.cluster.spike_times{curr_cluster}, ...
        xpr.bhv.stim_onset_ephys, spike_count_window, psth_bin_size);
    
    % Spike counts for correct trials
    [~,~,~,~,cluster_spike_count_correct(curr_cluster,:)] = ...
        psthRasterAndCounts(spikes.cluster.spike_times{curr_cluster}, ...
        xpr.bhv.stim_onset_ephys(xpr.bhv.correct), spike_count_window, psth_bin_size);
    
    % Spike counts for incorrect trials
    [~,~,~,~,cluster_spike_count_incorrect(curr_cluster,:)] = ...
        psthRasterAndCounts(spikes.cluster.spike_times{curr_cluster}, ...
        xpr.bhv.stim_onset_ephys(~xpr.bhv.correct), spike_count_window, psth_bin_size);
    
end

% Get spike counts by condition
[count_mean,count_sem,condition] = grpstats(cluster_spike_count',xpr.bhv.condition,{'mean','sem','gname'});
count_mean_zscore = zscore(count_mean,[],1);

[count_correct_mean,count_correct_sem,correct_condition] = ...
    grpstats(cluster_spike_count_correct',xpr.bhv.condition(xpr.bhv.correct,:),{'mean','sem','gname'});
count_correct_mean_zscore = zscore(count_correct_mean,[],1);

[count_incorrect_mean,count_incorrect_sem,incorrect_condition] = ...
    grpstats(cluster_spike_count_incorrect',xpr.bhv.condition(~xpr.bhv.correct,:),{'mean','sem','gname'});
count_incorrect_mean_zscore = zscore(count_incorrect_mean,[],1);

condition = cellfun(@(x) str2num(x),condition)*100;
condition_r_idx = find(condition(:,2) == 0);
condition_l_idx = find(condition(:,1) == 0);

% Get spike counts by depth
count_mean_depth = grpstats(count_mean_zscore',depth_group);
count_correct_mean_depth = grpstats(count_correct_mean_zscore',depth_group);
count_incorrect_mean_depth = grpstats(count_incorrect_mean_zscore',depth_group);

% Plot contrast response functions
% (in cortex)
cortex_groups = depth_groups_used > 0 & depth_groups_used <= n_depth_groups+1;
subcortex_group = depth_groups_used == n_depth_groups+1;

figure('Name',['Animal ' xpr.animal ', day ' xpr.day]);
use_conditions = condition_l_idx;

subplot(1,3,1); hold on;
col = repmat(linspace(0,0.8,sum(cortex_groups))',1,3);
set(gca,'ColorOrder',col);
plot(count_mean_depth(cortex_groups,use_conditions)','linewidth',2);
plot(count_mean_depth(subcortex_group,use_conditions),'r','linewidth',2);
set(gca,'XTick',1:length(use_conditions));
set(gca,'XTickLabel',cellfun(@num2str, ...
    mat2cell(condition(use_conditions,:),ones(size(use_conditions)),2),'uni',false))
ylim([-0.5,1]);
ylabel('Z-scored spike count');
xlabel('Contrasts')
title('All trials');

subplot(1,3,2); hold on;
col = repmat(linspace(0,0.8,sum(cortex_groups))',1,3);
set(gca,'ColorOrder',col);
plot(count_correct_mean_depth(cortex_groups,use_conditions)','linewidth',2);
plot(count_correct_mean_depth(subcortex_group,use_conditions),'r','linewidth',2);
set(gca,'XTick',1:length(use_conditions));
set(gca,'XTickLabel',cellfun(@num2str, ...
    mat2cell(condition(use_conditions,:),ones(size(use_conditions)),2),'uni',false));
ylim([-0.5,1]);
ylabel('Z-scored spike count');
xlabel('Contrasts')
title('Correct trials');

subplot(1,3,3); hold on;
col = repmat(linspace(0,0.8,sum(cortex_groups))',1,3);
set(gca,'ColorOrder',col);
plot(count_incorrect_mean_depth(cortex_groups,use_conditions)','linewidth',2);
plot(count_incorrect_mean_depth(subcortex_group,use_conditions),'r','linewidth',2);
set(gca,'XTick',1:length(use_conditions));
set(gca,'XTickLabel',cellfun(@num2str, ...
    mat2cell(condition(use_conditions,:),ones(size(use_conditions)),2),'uni',false));
ylim([-0.5,1]);
ylabel('Z-scored spike count');
xlabel('Contrasts')
title('Incorrect trials');

legend([cellfun(@(x) ['Depth: ' num2str(round(x*100)/100)], ...
    num2cell(depth_group_centers(cortex_groups)),'uni',false);'Subcortex']);


%% Tuning curve for mua based on depth

% Look at number of templates at each depth to estimate cortex boundaries
h = figure('Position',[94,122,230,820]);
plotSpread(spikes.template.depth,'distributionColors','k')
title('Mark cortex top')
[~,cortex_top] = ginput(1);
title('Mark cortex bottom');
[~,cortex_bottom] = ginput(1);
close(h);
drawnow;
cortex_depth = cortex_top - cortex_bottom;

% Make depths relative to cortex top
cluster_depth_cortex = cortex_top - spikes.cluster.depth;

% Group spikes by evenly spaced location in cortex
% (0 = out of cortex, n_depth_groups+1 = subcortex)
n_depth_groups = 4;
depth_group_edges = [linspace(0,cortex_depth,n_depth_groups+1),Inf];
[depth_group_n,depth_group] = histc(cluster_depth_cortex,depth_group_edges);
depth_groups_used = unique(depth_group);
depth_group_centers = grpstats(cluster_depth_cortex,depth_group);

% Create MUA times grouped according to depth
mua_depth = nan(size(spikes.spike_clusters));
for curr_cluster = unique(spikes.spike_clusters)'
    mua_depth(spikes.spike_clusters == curr_cluster) = ...
        depth_group(curr_cluster);
end

mua_times = cell(max(depth_groups_used),1);
for curr_depth = 1:max(depth_groups_used)
   mua_times{curr_depth} = spikes.spike_times(mua_depth == curr_depth);
end

% Get PSTH for all clusters
raster_window = [-2,2];
spike_count_window = [0,2];
psth_bin_size = 0.001;
n_bins = diff(raster_window)/psth_bin_size + 1;

n_depths = length(mua_times);
n_trials = size(xpr.bhv.condition,1);

depth_spike_count = nan(n_depths,n_trials);
depth_spike_count_correct = nan(n_depths,sum(xpr.bhv.correct));
depth_spike_count_incorrect = nan(n_depths,sum(~xpr.bhv.correct));
for curr_depth = 1:n_depths;
        
    % Spike counts for all trials
    [~,~,~,~,depth_spike_count(curr_depth,:)] = ...
        psthRasterAndCounts(mua_times{curr_depth}, ...
        xpr.bhv.stim_onset_ephys, spike_count_window, psth_bin_size);
    
    % Spike counts for correct trials
    [~,~,~,~,depth_spike_count_correct(curr_depth,:)] = ...
        psthRasterAndCounts(mua_times{curr_depth}, ...
        xpr.bhv.stim_onset_ephys(xpr.bhv.correct), spike_count_window, psth_bin_size);
    
    % Spike counts for incorrect trials
    [~,~,~,~,depth_spike_count_incorrect(curr_depth,:)] = ...
        psthRasterAndCounts(mua_times{curr_depth}, ...
        xpr.bhv.stim_onset_ephys(~xpr.bhv.correct), spike_count_window, psth_bin_size);
    
end

% Get spike counts by condition
[count_mean,count_sem,condition] = grpstats(depth_spike_count',xpr.bhv.condition,{'mean','sem','gname'});
count_mean_zscore = zscore(count_mean,[],1);

[count_correct_mean,count_correct_sem,correct_condition] = ...
    grpstats(depth_spike_count_correct',xpr.bhv.condition(xpr.bhv.correct,:),{'mean','sem','gname'});
count_correct_mean_zscore = zscore(count_correct_mean,[],1);

[count_incorrect_mean,count_incorrect_sem,incorrect_condition] = ...
    grpstats(depth_spike_count_incorrect',xpr.bhv.condition(~xpr.bhv.correct,:),{'mean','sem','gname'});
count_incorrect_mean_zscore = zscore(count_incorrect_mean,[],1);

condition = cellfun(@(x) str2num(x),condition)*100;
condition_r_idx = find(condition(:,2) == 0);
condition_l_idx = find(condition(:,1) == 0);

% Plot contrast response functions
% (in cortex)
cortex_groups = 1:n_depth_groups;
subcortex_group = n_depth_groups+1;

figure('Name',['Animal ' xpr.animal ', day ' xpr.day]);
use_conditions = condition_l_idx;

subplot(1,3,1); hold on;
col = repmat(linspace(0,0.8,length(cortex_groups))',1,3);
set(gca,'ColorOrder',col);
plot(count_mean_zscore(use_conditions,cortex_groups),'linewidth',2);
plot(count_mean_zscore(use_conditions,subcortex_group),'r','linewidth',2);
set(gca,'XTick',1:length(use_conditions));
set(gca,'XTickLabel',cellfun(@num2str, ...
    mat2cell(condition(use_conditions,:),ones(size(use_conditions)),2),'uni',false))
ylabel('Z-scored spike count');
xlabel('Contrasts')
title('All trials');

subplot(1,3,2); hold on;
col = repmat(linspace(0,0.8,length(cortex_groups))',1,3);
set(gca,'ColorOrder',col);
plot(count_correct_mean_zscore(use_conditions,cortex_groups),'linewidth',2);
plot(count_correct_mean_zscore(use_conditions,subcortex_group),'r','linewidth',2);
set(gca,'XTick',1:length(use_conditions));
set(gca,'XTickLabel',cellfun(@num2str, ...
    mat2cell(condition(use_conditions,:),ones(size(use_conditions)),2),'uni',false));
ylabel('Z-scored spike count');
xlabel('Contrasts')
title('Correct trials');

subplot(1,3,3); hold on;
col = repmat(linspace(0,0.8,length(cortex_groups))',1,3);
set(gca,'ColorOrder',col);
plot(count_incorrect_mean_zscore(use_conditions,cortex_groups),'linewidth',2);
plot(count_incorrect_mean_zscore(use_conditions,subcortex_group),'r','linewidth',2);
set(gca,'XTick',1:length(use_conditions));
set(gca,'XTickLabel',cellfun(@num2str, ...
    mat2cell(condition(use_conditions,:),ones(size(use_conditions)),2),'uni',false));
ylabel('Z-scored spike count');
xlabel('Contrasts')
title('Incorrect trials');

legend([cellfun(@(x) ['Depth: ' num2str(round(x))], ...
    num2cell(depth_group_centers(depth_groups_used > 0 & ...
    depth_groups_used < n_depth_groups+1)),'uni',false);'Subcortex']);


%% Raster plot by depth (OLD LOAD STRUCTURE)

animal = '65';
days = {'20151028','20151029','20151030','20151031','20151101','20151102','20151103'};

for curr_day = 1:length(days)
    clearvars -except animal days curr_day
    
    [spikes,xpr] = AP_load_ephys(animal,days{curr_day});
    
%     % Look at number of templates at each depth to estimate cortex boundaries
%     h = figure('Position',[94,122,230,820]);
%     plotSpread(spikes.template.depth,'distributionColors','k')
%     title('Mark cortex top')
%     [~,cortex_top] = ginput(1);
%     title('Mark cortex bottom');
%     [~,cortex_bottom] = ginput(1);
%     close(h);
%     drawnow;
    % just use whole depth for now
    cortex_top = max(spikes.template.depth)+1;
    cortex_bottom = min(spikes.template.depth)-1;
    cortex_depth = cortex_top - cortex_bottom;
    
    % Make depths relative to cortex top
    template_depth_cortex = cortex_top - spikes.template.depth;
    cluster_depth_cortex = cortex_top - spikes.cluster.depth;
    
    % Group spikes by evenly spaced location in cortex
    % (0 = out of cortex, n_depth_groups+1 = subcortex)
    n_depth_groups = 15;
    depth_group_edges = [linspace(0,cortex_depth,n_depth_groups+1),Inf];
    [depth_group_n,depth_group] = histc(cluster_depth_cortex,depth_group_edges);
    depth_groups_used = unique(depth_group);
    depth_group_centers = grpstats(cluster_depth_cortex,depth_group);
    
    % Create MUA times grouped according to depth
    mua_depth = nan(size(spikes.spike_clusters));
    for curr_cluster = unique(spikes.spike_clusters)'
        mua_depth(spikes.spike_clusters == curr_cluster) = ...
            depth_group(curr_cluster);
    end
    
    mua_times = cell(max(depth_groups_used),1);
    for curr_depth = 1:max(depth_groups_used)
        mua_times{curr_depth} = spikes.spike_times(mua_depth == curr_depth);
    end
    
    % Get PSTH for all clusters
    raster_window = [-2,2];
    spike_count_window = [0,2];
    psth_bin_size = 0.001;
    n_bins = diff(raster_window)/psth_bin_size + 1;
    
    n_depths = length(mua_times);
    n_trials = size(xpr.bhv.condition,1);
    
    depth_psth = nan(n_depths,n_bins);
    depth_spike_count = nan(n_depths,n_trials);
    depth_spike_count_correct = nan(n_depths,sum(xpr.bhv.correct));
    depth_spike_count_incorrect = nan(n_depths,sum(~xpr.bhv.correct));
    for curr_depth = 1:n_depths;
        
        % Raster for all trials
        [depth_psth(curr_depth,:),bins] = ...
            psthRasterAndCounts(mua_times{curr_depth}, ...
            xpr.bhv.stim_onset_ephys, raster_window, psth_bin_size);
        
        % Spike counts for all trials
        [~,~,~,~,depth_spike_count(curr_depth,:)] = ...
            psthRasterAndCounts(mua_times{curr_depth}, ...
            xpr.bhv.stim_onset_ephys, spike_count_window, psth_bin_size);
        
        % Spike counts for correct trials
        [~,~,~,~,depth_spike_count_correct(curr_depth,:)] = ...
            psthRasterAndCounts(mua_times{curr_depth}, ...
            xpr.bhv.stim_onset_ephys(xpr.bhv.correct), spike_count_window, psth_bin_size);
        
        % Spike counts for incorrect trials
        [~,~,~,~,depth_spike_count_incorrect(curr_depth,:)] = ...
            psthRasterAndCounts(mua_times{curr_depth}, ...
            xpr.bhv.stim_onset_ephys(~xpr.bhv.correct), spike_count_window, psth_bin_size);
        
    end
    
    % Gaussian smooth psth
    smooth_size = 15;
    gw = gausswin(round(smooth_size*6),3)';
    smWin = gw./sum(gw);
    depth_psth_smooth = conv2(depth_psth, smWin, 'same');
    
    depth_psth_zscore = zscore(depth_psth_smooth,[],2);
    
    % Plot grouping and PSTH by depth
    figure('Name',['Animal ' xpr.animal ', day ' xpr.day]);
    
    subplot(1,4,1);
    plotSpread(template_depth_cortex,'distributionColors','k')
    set(gca,'XTick',[]);
    ylabel('Template depth (\mum)');
    set(gca,'YDir','Reverse');
    for curr_edge = 2:length(depth_group_edges-1)
        line(xlim,repmat(depth_group_edges(curr_edge),2,1),'color','r')
    end
    
    subplot(1,4,2:4);
    plot_spacing = 5;
    depth_psth_plot = bsxfun(@plus,depth_psth_zscore, ...
        plot_spacing*transpose(size(depth_psth,1):-1:1))';
    plot(bins,depth_psth_plot,'k','linewidth',2);
    set(gca,'YTick',plot_spacing:plot_spacing:plot_spacing*n_depths);
    set(gca,'YTickLabel',flipud(round(depth_group_centers(depth_groups_used > 0))));
    ylabel('MUA depth (\mum)')
    xlabel('Time from stimulus onset (s)')
    title('Z-scored PSTH')
    
    
end

%% Raster plot by depth

align_times = stim_onsets(ismember(stimIDs,[90]));

% Group by depth
n_depth_groups = 6;
depth_group_edges = linspace(1000,max(spikeDepths),n_depth_groups+1);
depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
depth_group_edges(end) = Inf;
depth_group = discretize(spikeDepths,depth_group_edges);
depth_groups_used = unique(depth_group);

% Create MUA times grouped according to depth
mua_times = cell(n_depth_groups,1);
for curr_depth = 1:n_depth_groups
    mua_times{curr_depth} = spike_times_timeline(depth_group == curr_depth);
end

% PSTHs
raster_window = [-0.5,2.5];
psth_bin_size = 0.001;

depth_psth = nan(n_depth_groups,diff(raster_window)/psth_bin_size);

for curr_depth = 1:n_depth_groups
    [depth_psth(curr_depth,:),bins,rasterX,rasterY,spikeCounts] = psthAndBA( ...
        mua_times{curr_depth}, ...
        align_times, ...
        raster_window, psth_bin_size);
end

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
psth_smooth = conv2(depth_psth, smWin, 'same');
trace_spacing = max(psth_smooth(:));
figure; AP_stackplot(psth_smooth(:,20:end-20)',bins(20:end-20), ...
    10,true,'k',depth_group_centers);
line([0,0],ylim,'linestyle','--','color','k');
ylabel('Depth (\mum)');
xlabel('Time from stim onset (s)')
title('Population raster by depth');

%% Stim-triggered LFP by depth

% Group by depth
n_depth_groups = 20;
depth_group_edges = linspace(0,max(templateDepths),n_depth_groups+1);
depth_group_edges(end) = Inf;
depth_group_centers = depth_group_edges(1:end-1) + diff(depth_group_edges)./2;

% Get mean LFP by depth
channel_depth_grp = discretize(sort(channel_positions(:,2),'ascend'),depth_group_edges);
lfp_depth_mean = grpstats(lfp,channel_depth_grp);

% Stimulus-triggered LFP by stim
n_stim = length(unique(stimIDs));
align_times = stim_onsets;

lfp_window = [-1,5];
t_space = 0.001;
lfp_stim_mean = nan(n_stim,1+diff(lfp_window)/t_space,size(lfp_depth_mean,1));
lfp_stim_sem = nan(n_stim,1+diff(lfp_window)/t_space,size(lfp_depth_mean,1));

for curr_stim = 1:n_stim
    
    use_align = align_times(stimIDs == curr_stim);
    use_align_t = bsxfun(@plus,use_align,lfp_window(1):t_space:lfp_window(2));
    
    curr_lfp_stim = interp1(lfp_t_timeline,lfp_depth_mean',use_align_t);
    
    lfp_stim_mean(curr_stim,:,:) = nanmean(curr_lfp_stim,1);
    lfp_stim_sem(curr_stim,:,:) = nanstd(curr_lfp_stim,[],1)./sqrt(sum(~isnan(curr_lfp_stim),1));

end

plot_t = lfp_window(1):t_space:lfp_window(2);

% Plot one stim across depths
plot_stim = 5;
trace_spacing = 5000;
plot_lfp = squeeze(lfp_stim_mean(plot_stim,:,:));
yvals = trace_spacing*[1:size(plot_lfp,2)];
figure;AP_stackplot(plot_lfp,plot_t,trace_spacing,[],copper(size(plot_lfp,2)))
set(gca,'YTick',yvals);
set(gca,'YTickLabel',sort(depth_group_centers,'descend'));
ylabel('Depth (\mum)');
title('Stimulus-triggered LFP');
xlabel('Time from stim onset (s)');

% Plot one depth across stims
plot_depth = 1;
trace_spacing = 2000;
plot_lfp = squeeze(lfp_stim_mean(:,:,plot_depth))';
yvals = 1500*[1:size(plot_lfp,2)];
figure;AP_stackplot(plot_lfp,plot_t,trace_spacing,[],copper(size(plot_lfp,2)))
set(gca,'YTick',yvals);
set(gca,'YTickLabel',n_stim:-1:1);
ylabel('Stimulus)');
title('Stimulus-triggered LFP');
xlabel('Time from stim onset (s)');



%% Raster GUI with new load structures

trial_condition_rl = (xpr.bhv.condition(1,:) > xpr.bhv.condition(2,:)) - ...
    (xpr.bhv.condition(2,:) > xpr.bhv.condition(1,:));
spike_times_sec = double(spikes.spike_times);
raster_window = [-2,2];
psthViewer(spike_times_sec,spikes.spike_clusters, ...
    xpr.bhv.stim_onset_ephys,raster_window,trial_condition_rl);

%% Raster plot by depth, selected trials

clearvars -except animal days curr_day

[spikes,xpr] = AP_load_ephys(animal,days{curr_day});

% Look at number of templates at each depth to estimate cortex boundaries
h = figure('Position',[94,122,230,820]);
plotSpread(spikes.template.depth,'distributionColors','k')
title('Mark cortex top')
[~,cortex_top] = ginput(1);
title('Mark cortex bottom');
[~,cortex_bottom] = ginput(1);
close(h);
drawnow;
cortex_depth = cortex_top - cortex_bottom;

% Make depths relative to cortex top
template_depth_cortex = cortex_top - spikes.template.depth;
cluster_depth_cortex = cortex_top - spikes.cluster.depth;

% Group spikes by evenly spaced location in cortex
% (0 = out of cortex, n_depth_groups+1 = subcortex)
n_depth_groups = 6;
depth_group_edges = [linspace(0,cortex_depth,n_depth_groups+1),Inf];
[depth_group_n,depth_group] = histc(cluster_depth_cortex,depth_group_edges);
depth_groups_used = unique(depth_group);
depth_group_centers = grpstats(cluster_depth_cortex,depth_group);

% Create MUA times grouped according to depth
mua_depth = nan(size(spikes.spike_clusters));
for curr_cluster = unique(spikes.spike_clusters)'
    mua_depth(spikes.spike_clusters == curr_cluster) = ...
        depth_group(curr_cluster);
end

mua_times = cell(max(depth_groups_used),1);
for curr_depth = 1:max(depth_groups_used)
    mua_times{curr_depth} = spikes.spike_times(mua_depth == curr_depth);
end

% Get PSTH for all clusters
raster_window = [-2,2];
spike_count_window = [0,2];
psth_bin_size = 0.001;
n_bins = diff(raster_window)/psth_bin_size + 1;

% Trials to use
completed_trials = xpr.block.numCompletedTrials;
use_trials = (xpr.bhv.condition(1:completed_trials,1) == 0 & ...
    xpr.bhv.condition(1:completed_trials,2) == 0.5);

n_depths = length(mua_times);
n_trials = sum(use_trials);

%%%% WORKING ON THIS HERE: psth by correct
depth_psth = nan(n_depths,n_bins,2);

for curr_depth = 1:n_depths;
    
    % Raster for all trials
    [depth_psth(curr_depth,:),bins] = ...
        psthRasterAndCounts(mua_times{curr_depth}, ...
        xpr.bhv.stim_onset_ephys, raster_window, psth_bin_size);
    
    % Spike counts for all trials
    [~,~,~,~,depth_spike_count(curr_depth,:)] = ...
        psthRasterAndCounts(mua_times{curr_depth}, ...
        xpr.bhv.stim_onset_ephys, spike_count_window, psth_bin_size);
    
    % Spike counts for correct trials
    [~,~,~,~,depth_spike_count_correct(curr_depth,:)] = ...
        psthRasterAndCounts(mua_times{curr_depth}, ...
        xpr.bhv.stim_onset_ephys(xpr.bhv.correct), spike_count_window, psth_bin_size);
    
    % Spike counts for incorrect trials
    [~,~,~,~,depth_spike_count_incorrect(curr_depth,:)] = ...
        psthRasterAndCounts(mua_times{curr_depth}, ...
        xpr.bhv.stim_onset_ephys(~xpr.bhv.correct), spike_count_window, psth_bin_size);
    
end

% Gaussian smooth psth
smooth_size = 15;
gw = gausswin(round(smooth_size*6),3)';
smWin = gw./sum(gw);
depth_psth_smooth = conv2(depth_psth, smWin, 'same');

depth_psth_zscore = zscore(depth_psth_smooth,[],2);

% Plot grouping and PSTH by depth
figure('Name',['Animal ' xpr.animal ', day ' xpr.day]);

subplot(1,4,1);
plotSpread(template_depth_cortex,'distributionColors','k')
set(gca,'XTick',[]);
ylabel('Template depth (\mum)');
set(gca,'YDir','Reverse');
for curr_edge = 2:length(depth_group_edges-1)
    line(xlim,repmat(depth_group_edges(curr_edge),2,1),'color','r')
end

subplot(1,4,2:4);
plot_spacing = 5;
depth_psth_plot = bsxfun(@plus,depth_psth_zscore, ...
    plot_spacing*transpose(size(depth_psth,1):-1:1))';
plot(bins,depth_psth_plot,'k','linewidth',2);
set(gca,'YTick',plot_spacing:plot_spacing:plot_spacing*n_depths);
set(gca,'YTickLabel',flipud(round(depth_group_centers(depth_groups_used > 0))));
ylabel('MUA depth (\mum)')
xlabel('Time from stimulus onset (s)')
title('Z-scored PSTH')



%% For demonstration: plot narrow vs. broad spiking waveforms

[spikes,xpr] = AP_load_ephys('65',20151102);

figure;

subplot(1,2,1); hold on;
plot(spikes.cluster.duration(spikes.cluster.group.good),spikes.cluster.fwhm(spikes.cluster.group.good),'.k')
line([11,11],ylim,'color','r');
xlabel('Duration (arb)')
ylabel('FWHM')

subplot(1,2,2); hold on;
plot(spikes.cluster.waveform(spikes.cluster.broad & spikes.cluster.group.good,:)','k')
plot(spikes.cluster.waveform(spikes.cluster.narrow & spikes.cluster.group.good,:)','r')
xlim([20,80]);
axis off;

%% CB1 test recordings: xypos

% days = 160510,2016-05-23

animal = 'CB1';
day = '2016-05-23';
experiment = 2;

data_path = ['\\basket.cortexlab.net\data\ajpeters\' animal filesep day filesep num2str(experiment)];

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

% (to load in raw data)
% ephys_fid = fopen([data_dir filesep 'ephys.dat'],'r');
% ephys = fread(ephys_fid,'int16');
% fclose(ephys_fid);

timeline_filename = get_cortexlab_filename(['M111111_' animal],day,experiment,'timeline');
load(timeline_filename);
timeline_sample_rate = Timeline.hw.daqSampleRate;
photodiode_idx = strcmp({Timeline.hw.inputs.name},'photoDiode');

thresh = max(Timeline.rawDAQData(:,photodiode_idx))/2;
photodiode_trace = Timeline.rawDAQData(:,photodiode_idx) > thresh;
photodiode_flip = find((~photodiode_trace(1:end-1) & photodiode_trace(2:end)) | ...
    (photodiode_trace(1:end-1) & ~photodiode_trace(2:end)));

photodiode = struct('timestamps',[],'values',[]);
photodiode.timestamps = photodiode_flip./timeline_sample_rate;
photodiode.values = photodiode_trace(photodiode_flip);

% Get stim onset times
% Check that sync matches photodiode number
if length(sync.timestamps) == length(photodiode.timestamps)
    refresh_rate_cutoff = 1/10;
    stim_onset_idx = [1,find(diff(sync.timestamps) > refresh_rate_cutoff) + 1];
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

% Load mpep protocol
protocol_filename = get_cortexlab_filename(['M111111_' animal],day,experiment,'protocol');
load(protocol_filename);

% Get stimulus IDs
stimIDs = zeros(size(stim_onset_idx));
for q = 1:size(Protocol.seqnums,1)
    stimIDs(Protocol.seqnums(q,:)) = q;
end

% Manually go through individual PSTHs
raster_window = [-5,5];
psthViewer(spike_times,spike_templates, ...
    sync.timestamps(stim_onset_idx),raster_window,stimIDs);

% Plot PSTH of whole population per stimulus
psth_bin_size = 0.001;
psth = nan(max(stimIDs),diff(raster_window)/psth_bin_size + 1);
for curr_stim = 1:max(stimIDs);
    [psth(curr_stim,:),bins,rasterX,rasterY,spikeCounts] = ...
        psthRasterAndCounts(spike_times, ...
        sync.timestamps(stim_onset_idx(stimIDs == curr_stim))', ...
        raster_window, psth_bin_size);
end

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
psth_smooth = conv2(psth, smWin, 'same');
psth_plot = bsxfun(@plus,mat2gray(psth_smooth),transpose(1:size(psth,1)));
figure; hold on;
plot(bins,psth_plot','k','linewidth',2);
line([0,0],ylim,'linestyle','--','color','k');
ylabel('Stim ID');
xlabel('Time from stim onset')

% Plot raster by depth
n_depth_groups = 15;
depth_group_edges = linspace(0,max(templateDepths),n_depth_groups+1);
depth_group_edges(end) = Inf;
[depth_group_n,depth_group] = histc(spikeDepths,depth_group_edges);
depth_groups_used = unique(depth_group);
depth_group_centers = grpstats(spikeDepths,depth_group);

raster_window = [-2,3];
psth_bin_size = 0.001;
psth = nan(n_depth_groups,diff(raster_window)/psth_bin_size + 1);
for curr_depth = 1:n_depth_groups;
    [psth(curr_depth,:),bins,rasterX,rasterY,spikeCounts] = ...
        psthRasterAndCounts(spike_times(depth_group == curr_depth), ...
        sync.timestamps(stim_onset_idx)', ...
        raster_window, psth_bin_size);
end

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
psth_smooth = conv2(psth, smWin, 'same');
psth_plot = bsxfun(@plus,mat2gray(psth_smooth),transpose(1:size(psth,1)));
figure; hold on;
plot(bins,psth_plot','k','linewidth',2);
line([0,0],ylim,'linestyle','--','color','k');
ylabel('Depth group');
xlabel('Time from stim onset')

% Load LFP
lfp_filename = [data_path filesep 'lfp.dat'];
fid = fopen(lfp_filename);
lfp = fread(fid,[128,inf],'int16');
fclose(fid);

% Sort LFP by depth, keep depth
[lfp_depth,depth_sort_idx] = sort(channel_positions(:,2));
lfp = lfp(channel_map(depth_sort_idx)+1,:);

lfp_sample_rate = 1000;
% LFP power spectrum
NFFT = 2^nextpow2(size(lfp,2));
hippocampus_lfp = nanmean(lfp(1:37,:),1);
cortex_lfp = nanmean(lfp(77:end,:),1);
[Pxx,F] = pwelch([hippocampus_lfp',cortex_lfp'],10000,[],NFFT,lfp_sample_rate);
figure;loglog(F,Pxx);
legend({'Hippocampus','Cortex'})

% Stimulus-triggered LFP
lfp_t = (1:size(lfp,2))/lfp_sample_rate;
lfp_stim_window_sec = [-1,3];
lfp_stim_window_samples = round(lfp_stim_window_sec*lfp_sample_rate);
lfp_stim_window_plot = linspace(lfp_stim_window_sec(1),lfp_stim_window_sec(2),diff(lfp_stim_window_samples));
stim_onset_t = sync.timestamps(stim_onset_idx);
stim_onset_lfp = round(stim_onset_t*lfp_sample_rate);
lfp_stim_surround_x = arrayfun(@(x) stim_onset_lfp(x)+lfp_stim_window_samples(1): ...
    stim_onset_lfp(x)+lfp_stim_window_samples(2),1:length(stim_onset_lfp),'uni',false);
lfp_stim_surround = permute(nanmean(reshape(lfp(:,horzcat(lfp_stim_surround_x{:}))', ...
    [],length(stim_onset_t),size(lfp,1)),2),[3,1,2]);

figure;imagesc(lfp_stim_surround);colormap(gray);
set(gca,'YDir','normal');
xtick = get(gca,'XTick');
set(gca,'XTickLabel',round(lfp_stim_window_plot(xtick)*10)/10);
line(repmat(find(lfp_stim_window_plot > 0,1),1,2),ylim,'color','r','linewidth',2)
ylabel('Channel (sorted by depth)')
xlabel('Time from stim onset (sec)')

% Stimulus-triggered LFP by stim
for curr_stim = 1:max(stimIDs)
    lfp_t = (1:size(lfp,2))/lfp_sample_rate;
    lfp_stim_window_sec = [-1,3];
    lfp_stim_window_samples = round(lfp_stim_window_sec*lfp_sample_rate);
    lfp_stim_window_plot = linspace(lfp_stim_window_sec(1),lfp_stim_window_sec(2),diff(lfp_stim_window_samples));
    stim_onset_t = sync.timestamps(stim_onset_idx(stimIDs == curr_stim));
    stim_onset_lfp = round(stim_onset_t*lfp_sample_rate);
    lfp_stim_surround_x = arrayfun(@(x) stim_onset_lfp(x)+lfp_stim_window_samples(1): ...
        stim_onset_lfp(x)+lfp_stim_window_samples(2),1:length(stim_onset_lfp),'uni',false);
    lfp_stim_surround = permute(nanmean(reshape(lfp(:,horzcat(lfp_stim_surround_x{:}))', ...
        [],length(stim_onset_t),size(lfp,1)),2),[3,1,2]);
    
    hipp = 1:37;
    ctx = 77:120;
    
    figure; hold on
    x = linspace(min(lfp_stim_window_sec),max(lfp_stim_window_sec),size(lfp_stim_surround,2));
    plot(x,nanmean(lfp_stim_surround(hipp,:),1),'k');
    plot(x,nanmean(lfp_stim_surround(ctx,:),1),'r');
    ylabel('Voltage (arb)')
    xlabel('Time from stim onset (sec)')
    
    title(['Stim ' num2str(curr_stim)]);
end



%% CB1 test recordings: sparse noise

animal = 'CB1';
day = '2016-05-23';
experiment = 4;

data_path = ['\\basket.cortexlab.net\data\ajpeters\' animal filesep day filesep num2str(experiment)];

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

% (to load in raw data)
% ephys_fid = fopen([data_dir filesep 'ephys.dat'],'r');
% ephys = fread(ephys_fid,'int16');
% fclose(ephys_fid);

timeline_filename = get_cortexlab_filename(['M111111_' animal],day,experiment,'timeline');
load(timeline_filename);
timeline_sample_rate = Timeline.hw.daqSampleRate;
photodiode_idx = strcmp({Timeline.hw.inputs.name},'photoDiode');

thresh = max(Timeline.rawDAQData(:,photodiode_idx))/2;
photodiode_trace = Timeline.rawDAQData(:,photodiode_idx) > thresh;
photodiode_flip = find((~photodiode_trace(1:end-1) & photodiode_trace(2:end)) | ...
    (photodiode_trace(1:end-1) & ~photodiode_trace(2:end)));

photodiode = struct('timestamps',[],'values',[]);
photodiode.timestamps = photodiode_flip./timeline_sample_rate;
photodiode.values = photodiode_trace(photodiode_flip);

% Check that sync matches photodiode number
if length(sync.timestamps) ~= length(photodiode.timestamps)
    error('Different number of timeline vs ephys pulses');
end

% Load mpep protocol and generate stimuli
protocol_filename = get_cortexlab_filename(['M111111_' animal],day,experiment,'protocol');
load(protocol_filename);
hwinfo_filename = get_cortexlab_filename(['M111111_' animal],day,experiment,'hardware');
load(hwinfo_filename);

myScreenInfo.windowPtr = NaN; % so we can call the stimulus generation and it won't try to display anything
stimNum = 1;
ss = eval([Protocol.xfile(1:end-2) '(myScreenInfo, Protocol.pars(:,stimNum));']);

% convert stimulus info to array
imTextSeq = ss.ImageTextures(ss.ImageSequence(1:end-1)); % excluding the last ImageSequence here as a hack to make the right number of photodiode events (?)
q = cat(3,imTextSeq{:});
stimArray = repmat(q, [1 1 Protocol.nrepeats]); clear q; % in my case I had three repeats of the same stimulus. Otherwise you want to generate a differen "ss" object for each stimulus, above

nX = size(stimArray,1);
nY = size(stimArray,2);
stimArrayZeroPad = cat(3,zeros(size(stimArray,1), size(stimArray,2),1), stimArray);
for x = 1:nX
    for y = 1:nY
        stimEventTimes{x,y,1} = photodiodeFlips(stimArrayZeroPad(x,y,1:end-1)==0 & ...
            stimArrayZeroPad(x,y,2:end)==1); % going from grey to white
        stimEventTimes{x,y,2} = photodiodeFlips(stimArrayZeroPad(x,y,1:end-1)==0 & ...
            stimArrayZeroPad(x,y,2:end)==-1); % going from grey to black
    end
end


%% Raster aligned to stimuli

% use_spikes_idx = ismember(spike_templates,find(templateDepths >= 0 & templateDepths <= 1500));
use_spikes_idx = ismember(spike_templates,find(templateDepths > 3500 & templateDepths < 4000)) & ...
   (ismember(spike_templates,find(msn)));

% use_spikes_idx = true(size(spike_times_timeline));

use_spikes = spike_times_timeline(use_spikes_idx);
use_spike_templates = spike_templates(use_spikes_idx);

align_times = stim_onsets(ismember(stimIDs,[0]));

% PSTHs
raster_window = [-0.5,3];
psthViewer(use_spikes,use_spike_templates, ...
    stim_onsets,raster_window,stimIDs);

psth_bin_size = 0.001;
[psth,bins,rasterX,rasterY,spikeCounts] = psthAndBA( ...
    use_spikes,align_times, ...
    raster_window, psth_bin_size);

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
psth_smooth = conv2(psth,smWin,'same');
figure; hold on;
plot(bins(20:end-20),psth_smooth(20:end-20)','k','linewidth',2);
line([0,0],ylim,'linestyle','--','color','k');
ylabel('Population spikes');
xlabel('Time from stim onset')

% PSTH by stim condition
stim_psth = nan(length(unique(stimIDs)),length(psth_smooth));
unique_stims = unique(stimIDs);
for curr_stim_idx = 1:length(unique_stims);   
    
    [psth,bins,rasterX,rasterY,spikeCounts] = psthAndBA( ...
        use_spikes,stim_onsets(stimIDs == unique_stims(curr_stim_idx)), ...
        raster_window, psth_bin_size);
    
    stim_psth(curr_stim_idx,:) = psth;
    
end
stim_psth_smooth = conv2(stim_psth,smWin,'same');
figure; hold on;
trace_spacing = max(stim_psth_smooth(:));
AP_stackplot(stim_psth_smooth(:,20:end-20)',bins(20:end-20),trace_spacing,false,[],unique(stimIDs))
xlabel('Time from stim onset')
ylabel('Population spikes (by stim)');
line([0,0],ylim,'linestyle','--','color','k');


%% Sparse noise receptive fields

% Generate stimuli from protocol and hardware
myScreenInfo.windowPtr = NaN; % so we can call the stimulus generation and it won't try to display anything
stimNum = 1;
ss = eval([Protocol.xfile(1:end-2) '(myScreenInfo, Protocol.pars(:,stimNum));']);
stim_screen = cat(3,ss.ImageTextures{:});
ny = size(stim_screen,1);
nx = size(stim_screen,2);

switch lower(photodiode_type)
    case 'flicker'
        % Check for case of mismatch between photodiode and stimuli:
        % odd number of stimuli, but one extra photodiode flip to come back down
        if mod(size(stim_screen,3),2) == 1 && ...
                length(photodiode.timestamps) == size(stim_screen,3) + 1;
            photodiode.timestamps(end) = [];
            photodiode.values(end) = [];
            warning('Odd number of stimuli, removed last photodiode');
        end
        
        % If there's still a mismatch, break
        if size(stim_screen,3) ~= length(photodiode.timestamps);
            warning([num2str(size(stim_screen,3)) ' stimuli, ', ...
                num2str(length(photodiode.timestamps)) ' photodiode pulses']);
            
            % Try to estimate which stim were missed by time difference
            photodiode_diff = diff(photodiode.timestamps);
            max_regular_diff_time = prctile(diff(photodiode.timestamps),99);
            skip_cutoff = max_regular_diff_time*2;
            photodiode_skip = find(photodiode_diff > skip_cutoff);
            est_n_pulse_skip = ceil(photodiode_diff(photodiode_skip)/max_regular_diff_time)-1;
            stim_skip = cell2mat(arrayfun(@(x) photodiode_skip(x):photodiode_skip(x)+est_n_pulse_skip(x)-1, ...
                1:length(photodiode_skip),'uni',false));
            
            if isempty(est_n_pulse_skip) || length(photodiode.timestamps) + sum(est_n_pulse_skip) ~= size(stim_screen,3)
                error('Can''t match photodiode events to stimuli')
            end
        end
        
        stim_times = photodiode.timestamps;
        
    case 'steady'
        % If the photodiode is on steady: extrapolate the stim times
        if length(photodiode.timestamps) ~= 2
            error('Steady photodiode, but not 2 flips')
        end
        stim_duration = diff(photodiode.timestamps)/size(stim_screen,3);
        stim_times = linspace(photodiode.timestamps(1), ...
            photodiode.timestamps(2)-stim_duration,size(stim_screen,3))';
        
end

% Get stim times vector (x,y)
nY = size(stim_screen,1);
nX = size(stim_screen,2);
stim_times_grid = cell(nY,nX);
for x = 1:nX
    for y = 1:nY      
        align_stims = (stim_screen(y,x,2:end)~= 0) & ...
            (diff(stim_screen(y,x,:),[],3) ~= 0);
        align_times = stim_times(find(align_stims)+1);
        stim_times_grid{y,x} = align_times;
    end
end

% Vectorize stim times by y,x position
[stim_x,stim_y] = meshgrid(1:nX,1:nY);
stim_positions = cellfun(@(x,y,times) repmat([y,x],length(times),1), ...
    num2cell(stim_x),num2cell(stim_y),stim_times_grid,'uni',false);

params = struct;
params.makePlots = false;
params.useSVD = false;
rf_map = nan(nY,nX,max(spike_templates));
for curr_template = unique(spike_templates)'
    [rf_map(:,:,curr_template+1),stats] = sparseNoiseRF(spike_times_timeline(spike_templates == curr_template), ...
        vertcat(stim_times_grid{:}),vertcat(stim_positions{:}),params);
end

disp('done')

gauss_sigma = 1;
%rf_map_smooth = imgaussfilt(rf_map,gauss_sigma);
gauss_filt = fspecial('gaussian',[nY,nX],gauss_sigma);
rf_map_smooth = imfilter(rf_map,gauss_filt);

% % Get significant RF by shuffle
% n_shuff = 1000;
% rf_map_reshape = reshape(rf_map,nY*nX,[]);
% max_resp = nan(n_shuff,size(rf_map,3));
% for i = 1:n_shuff
%     warning off
%     %max_resp(i,:) = max(reshape(imgaussfilt(reshape(shake(rf_map_reshape, ...
%     %    1),nY,nX,[]),gauss_sigma),nY*nX,[]),[],1);
%     max_resp(i,:) = max(reshape(imfilter(reshape(shake(rf_map_reshape, ...
%         1),nY,nX,[]),gauss_filt),nY*nX,[]),[],1);
%     warning on
%     disp(i);
% end
% 
% cutoff_resp = prctile(max_resp,95);
% sig_rf = bsxfun(@gt,rf_map_smooth,permute(cutoff_resp,[1,3,2]));
% 



% Get stim-triggered average for each stimulus
use_spikes_idx = ismember(spike_templates,find(templateDepths >= 1300 & templateDepths <= 2500));
use_spikes = spike_times_timeline(use_spikes_idx);
use_spike_templates = spike_templates(use_spikes_idx);

% Get stim times vector (x,y)
stim_aligned_avg = cell(nY,nX);
raster_window = [-0.2,0.2];
smooth_size = 2;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
for x = 1:nX
    for y = 1:nY      
        
        psth_bin_size = 0.001;
        [psth,bins,rasterX,rasterY,spikeCounts] = psthAndBA( ...
            use_spikes,stim_times_grid{y,x}, ...
            raster_window, psth_bin_size);
        
        psth_smooth = conv2(psth,smWin,'same');
        stim_aligned_avg{y,x} = psth_smooth;
        
    end
end
stim_aligned_avg_cat = cell2mat(cellfun(@(x) permute(x,[1,3,2]),stim_aligned_avg,'uni',false));
AP_image_scroll(stim_aligned_avg_cat,raster_window(1):psth_bin_size:raster_window(2));
axis equal




%% Classify cell type

% Define cortical and striatal cells
str_depth = [0,Inf];

str_templates = templateDepths >= str_depth(1) & templateDepths <= str_depth(2);
non_str_templates = ~str_templates;

% Get firing rate
spike_rate = nan(max(spike_templates),1);
for curr_template = unique(spike_templates)'
    spike_rate(curr_template) = ...
        sum(spike_templates == curr_template)./ ...
        (max(spike_times_timeline) - min(spike_times_timeline));
end

% Get proportion of ISI > 2s
prop_long_isi = nan(max(spike_templates),1);
for curr_template = unique(spike_templates)'
    curr_spike_times = spike_times_timeline(spike_templates == curr_template);
    curr_isi = diff(curr_spike_times);
    
    prop_long_isi(curr_template) = sum(curr_isi(curr_isi > 2))./ ...
        (max(spike_times_timeline) - min(spike_times_timeline));
end

waveform_duration_cutoff = 400;

% Cortical classification (like Bartho JNeurophys 2004)
narrow = non_str_templates & templateDuration_us <= waveform_duration_cutoff;
wide = non_str_templates & templateDuration_us > waveform_duration_cutoff;

% Striatum classification (like Yamin/Cohen 2013)
long_isi_cutoff = 0.35;

msn = str_templates & ...
    templateDuration_us >= waveform_duration_cutoff & ...
    prop_long_isi >= long_isi_cutoff;

fsi = str_templates & ...
    templateDuration_us <= waveform_duration_cutoff & ...
    prop_long_isi <= long_isi_cutoff;

tan = str_templates & ...
    templateDuration_us > waveform_duration_cutoff & ...
    prop_long_isi < long_isi_cutoff;

uin = str_templates & ...
    templateDuration_us < waveform_duration_cutoff & ...
    prop_long_isi > long_isi_cutoff;

waveform_t = 1e3*((0:size(templates,2)-1)/ephys_sample_rate);

% Plot the waveforms and spike statistics
figure;

if any(wide) || any(narrow)
    subplot(2,2,1); hold on;
    p1 = plot(waveform_t,waveforms(wide,:)','k');
    p2 = plot(waveform_t,waveforms(narrow,:)','r');
    xlabel('Time (ms)')
    title('Not striatum');
    legend([p1(1),p2(1)],{'Wide','Narrow'})
end

subplot(2,2,2); hold on;
p1 = plot(waveform_t,waveforms(msn,:)','m');
p2 = plot(waveform_t,waveforms(fsi,:)','b');
p3 = plot(waveform_t,waveforms(tan,:)','g');
p4 = plot(waveform_t,waveforms(uin,:)','c');
xlabel('Time (ms)')
title('Striatum');
legend([p1(1),p2(1),p3(1),p4(1)],{'MSN','TAN','FSI','UIN'});

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
    templateDuration_us(tan)/1000, ...
    prop_long_isi(tan), ...
    spike_rate(tan),'b');

stem3( ...
    templateDuration_us(fsi)/1000, ...
    prop_long_isi(fsi), ...
    spike_rate(fsi),'g');

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
celltype_labels = {'Wide','Narrow','MSN','TAN','FSI','UIN'};
celltypes = wide.*1 + narrow.*2 + msn.*3 + tan.*4 + fsi.*5 + uin.*6;
use_colors = {'k','r','m','b','g','c'};

plot_celltypes = any([wide,narrow,msn,tan,fsi,uin],1);

figure; plotSpread(templateDepths,'categoryIdx', ...
    celltypes,'categoryColors',use_colors(plot_celltypes));
set(gca,'XTick',[]);
set(gca,'YDir','reverse');
ylabel('Depth (\mum)');
legend(celltype_labels(plot_celltypes));

%% MUA/LFP correlation by depth 

n_depth_groups = 50;
depth_group_edges = linspace(0,max(channel_positions(:,2)),n_depth_groups+1);
depth_group = discretize(templateDepths,depth_group_edges);
depth_group_centers = depth_group_edges(1:end-1)+diff(depth_group_edges);
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

% LFP (subtract the median across channels)
channel_depth_grp = discretize(channel_positions(:,2),depth_group_edges);
lfp_depth_mean = grpstats(lfp,channel_depth_grp);
lfp_depth_mean_mediansub = bsxfun(@minus,lfp_depth_mean,nanmedian(lfp_depth_mean,1));
lfp_corr = corrcoef(lfp_depth_mean_mediansub');

% Plot
figure; colormap(parula);

subplot(1,2,1);
imagesc(depth_group_centers,depth_group_centers,mua_corr);
xlabel('Depth (\mum)')
ylabel('Depth (\mum)')
title(['MUA correlation, ' num2str(spike_binning) 's bins'])
caxis([min(AP_itril(mua_corr,-1)),max(AP_itril(mua_corr,-1))]);
colorbar
axis square;

subplot(1,2,2);
imagesc(depth_group_centers,depth_group_centers,lfp_corr);
xlabel('Depth (\mum)')
ylabel('Depth (\mum)')
title('LFP correlation')
caxis([min(AP_itril(lfp_corr,-1)),max(AP_itril(lfp_corr,-1))]);
colorbar
axis square;

%% Spectral analysis
% Power spectrum
use_trace = roi_trace(1:end-1);
use_t = frame_t;

Fs = 1./median(diff(use_t));
L = length(use_trace);
NFFT = 2^nextpow2(L);
[P,F] = pwelch(double(use_trace)',[],[],NFFT,Fs);
Pc = smooth(P,50); 
figure;plot(F,log10(Pc),'k')
xlabel('Frequency');
ylabel('Log Power');

% Notch filter
freqs = [3 5];
for ff = 1:length(freqs)
    f = fdesign.notch('N,F0,Q',2,freqs(ff),10,Fs);
    h = design(f);
    % hfvt= fvtool(h,'Color','white');
    use_trace_filt = filter(h, use_trace')';
end

p = bandpower(use_trace,Fs,[3,5]);

% Spectrogram
spect_overlap = 50;
window_length = 5; % in seconds
window_length_samples = window_length/(1/Fs);
figure;spectrogram(use_trace,window_length_samples, ...
    round(spect_overlap/100*window_length_samples),[],Fs,'yaxis')
colormap(hot)

% Band power over time
spect_overlap = 80;
window_length = 3; % in seconds
window_length_samples = window_length/(1/Fs);

[s,f,t] = spectrogram(use_trace,window_length_samples, ...
    round(spect_overlap/100*window_length_samples),[],Fs);

N = window_length_samples; % window length
df = Fs/N; % frequency increment
s_squared = (s/Fs).*conj(s/Fs);  % Fs is used to normalize the FFT amplitudes
power_0_2 = 2*sum(s_squared( f >= 0.1 & f <= 2,:))*df; 
power_3_6 = 2*sum(s_squared( f >= 3 & f <= 6,:))*df; 
power_10_14 = 2*sum(s_squared( f >= 10 & f <= 14,:))*df; 

%% Spike --> spike regression by depth

spike_binning = 0.01; % seconds
corr_edges = spike_times_timeline(1):spike_binning:spike_times_timeline(end);

ctx_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 0 & templateDepths < 2400)));
binned_spikes_ctx = single(histcounts(ctx_spikes,corr_edges));

use_templates = good_templates(templateDepths(good_templates) > 2400);
binned_spikes_str = zeros(length(use_templates),length(corr_edges)-1,'single');
for curr_template_idx = 1:length(use_templates)    
    curr_template = use_templates(curr_template_idx);   
    use_spikes = spike_times_timeline(spike_templates == curr_template);    
    curr_binned_spikes = histcounts(use_spikes,corr_edges);      
    binned_spikes_str(curr_template_idx,:) = curr_binned_spikes;    
end

kernel_frames = -10:1:5;
lambda = 0;

[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(binned_spikes_str,binned_spikes_ctx,kernel_frames,lambda,true);
k = reshape(k,size(binned_spikes_str,1),[]);

figure;plot(kernel_frames*spike_binning,k);
line([0,0],ylim,'color','k')

k_max = max(k,[],2);
use_template_depths = templateDepths(use_templates+1);
figure;plot(use_template_depths,k_max,'.k');
ylabel('Maximum weight over time');
xlabel('Template depth')

% Group multiunit by depth
n_depth_groups = 10;
depth_group_edges = linspace(2000,double(max(channel_positions(:,2))),n_depth_groups+1);
depth_group_edges_use = depth_group_edges;
depth_group_edges_use(end) = Inf;

[depth_group_n,depth_group] = histc(spikeDepths,depth_group_edges_use);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges_use(1:end-1)+diff(depth_group_edges_use);

binned_spikes_str = zeros(n_depth_groups,length(corr_edges)-1,'single');
for curr_depth = 1:length(depth_group_edges_use)-1   
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    % Discretize spikes into frames and count spikes per frame
    binned_spikes_str(curr_depth,:) = histcounts(curr_spike_times,corr_edges);  
end

kernel_frames = -10:1:5;
lambda = 1000;

[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(binned_spikes_ctx,binned_spikes_str,kernel_frames,lambda,false);
k = reshape(k,size(binned_spikes_str,1),[]);

k_max = max(k,[],2);
figure;plot(depth_group_centers,k_max,'k','linewidth',2);
ylabel('Max weight over time');
xlabel('MUA depth');

%% Spike --> spike regression by cell type

% Bin all spikes
spike_binning = 0.01; % seconds
corr_edges = spike_times_timeline(1):spike_binning:spike_times_timeline(end);

binned_spikes = zeros(size(templates,1),length(corr_edges)-1,'single');
for curr_template = 1:size(templates,1)    
    use_spikes = spike_times_timeline(spike_templates == curr_template);    
    curr_binned_spikes = histcounts(use_spikes,corr_edges);      
    binned_spikes(curr_template,:) = curr_binned_spikes;    
end

kernel_time = [-0.5,0.5];
kernel_timepoints = kernel_time(1)/spike_binning:kernel_time(2)/spike_binning;
lambda = 0;

use_depth_templates = templateDepths > 0 & templateDepths < 2000;

msn_spikes = sum(binned_spikes(msn & use_depth_templates,:),1);
tan_spikes = sum(binned_spikes(tan & use_depth_templates,:),1);
fsi_spikes = sum(binned_spikes(fsi & use_depth_templates,:),1);
uin_spikes = sum(binned_spikes(uin & use_depth_templates,:),1);

[k,predicted_spikes,explained_var] = ...
    AP_regresskernel([tan_spikes;fsi_spikes;uin_spikes],msn_spikes, ...
    kernel_timepoints,lambda,true);

k = reshape(k,3,length(kernel_timepoints),[]);

figure;
subplot(1,2,1);
plot(kernel_timepoints*spike_binning,k','linewidth',2);
legend({'TAN','FSI','UIN'})

plot_frames = 100;
soft_reg_factor = 1e9;
x_autonorm = ifft((fft(msn_spikes).*conj(fft(fsi_spikes)))./(soft_reg_factor+fft(fsi_spikes).*conj(fft(fsi_spikes))));

t_shift = [frame_t(end-plot_frames+1:end)-frame_t(end)-1/framerate,frame_t(1:plot_frames)-frame_t(1)];

subplot(1,2,2);
plot(t_shift,[x_autonorm(end-plot_frames+1:end),x_autonorm(1:plot_frames)],'k','linewidth',2);
xlabel('Time (s)');
ylabel('Impluse response');



%% Raster plots to choiceworld events

stimIDs = signals_events.trialSideValues.*signals_events.trialContrastValues;

use_spikes_idx = ismember(spike_templates,find(templateDepths >= 1300 & templateDepths <= 2500));
use_spikes = spike_times_timeline(use_spikes_idx);
use_spike_templates = spike_templates(use_spikes_idx);

raster_window = [-2,4];
psth_bin_size = 0.001;

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

% PSTH by stim condition
stim_psth_hit = nan(length(unique(stimIDs)),diff(raster_window)/psth_bin_size);
stim_psth_miss = nan(length(unique(stimIDs)),diff(raster_window)/psth_bin_size);

unique_stims = unique(stimIDs);
for curr_stim_idx = 1:length(unique_stims);   
    
    curr_trials_hit = (stimIDs == unique_stims(curr_stim_idx)) & ...
        signals_events.hitValues == 1;
    if sum(curr_trials_hit) > 0
        [psth_hit,bins] = psthAndBA( ...
            use_spikes,signals_events.stimOnTimes(curr_trials_hit), ...
            raster_window, psth_bin_size);
        stim_psth_hit(curr_stim_idx,:) = psth_hit;
    end
    
    curr_trials_miss = (stimIDs == unique_stims(curr_stim_idx)) & ...
        signals_events.hitValues == 0;
    if sum(curr_trials_miss) > 0
        [psth_miss,bins] = psthAndBA( ...
            use_spikes,signals_events.stimOnTimes(curr_trials_miss), ...
            raster_window, psth_bin_size);
        stim_psth_miss(curr_stim_idx,:) = psth_miss;
    end
    
end
stim_psth_hit_smooth = conv2(stim_psth_hit,smWin,'same');
stim_psth_miss_smooth = conv2(stim_psth_miss,smWin,'same');

figure; hold on;
trace_spacing = max([stim_psth_hit_smooth(:);stim_psth_miss_smooth(:)]);
AP_stackplot(stim_psth_hit_smooth(:,20:end-20)',bins(20:end-20),trace_spacing,false,[],unique(stimIDs))
AP_stackplot(stim_psth_miss_smooth(:,20:end-20)',bins(20:end-20),trace_spacing,false,[],unique(stimIDs))
xlabel('Time from stim onset')
ylabel('Population spikes (by stim)');
line([0,0],ylim,'linestyle','--','color','k');







