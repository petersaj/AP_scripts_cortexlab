function [spikes,xpr] = AP_load_ephys(animal,day)
% [spikes,xpr] = AP_load_ephys(animal,day);
%
% INPUTS: 
% animal = animal string
% day = 8-digit day
%
% OUTPUTS: 
% spikes = structure containing all electrophysiology data
% xpr = structure containing other experiment data
%
% NOTE: this makes a lot of assumptions about file locations at the moment

error('The alignment for this is probably off: don''t use it');

%% Set paths, load data

% If day entered as number, convert to string
if isnumeric(day)
    day = num2str(day);
end

% Set paths for raw and processed data
raw_data_filename = ['\\zserver.cortexlab.net\Data\multichanspikes\' ...
    animal filesep day filesep day '_1.dat'];
processed_data_path = '\\basket.cortexlab.net\data\ajpeters';

if ~exist(raw_data_filename)
    error('No raw data found for animal and day')
end

% Get ephys sample rate from header
ephys_header_filename = [raw_data_filename(1:end-4) '.meta'];
ephys_header_id = fopen(ephys_header_filename);
ephys_header = textscan(ephys_header_id,'%s %s', 'delimiter',{' ='});
fclose(ephys_header_id);

sample_rate_idx = strcmp(ephys_header{1},'sRateHz');
ephys_sample_rate = str2num(ephys_header{2}{sample_rate_idx});

% Find filenames for behavior/input
timeline_filename = get_cortexlab_filename(animal,day,1,'timeline','dash');
parameters_filename = get_cortexlab_filename(animal,day,1,'parameters','dash');
block_filename = get_cortexlab_filename(animal,day,1,'block','dash');

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
    fclose(cluster_groups_id);
    
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
    cluster_groups_1idx = 1:(max(spike_clusters)+1);
    cluster_groups = struct( ...
        'good',true(max(cluster_groups_1idx),1), ...
        'mua',false(max(cluster_groups_1idx),1), ...
        'noise',false(max(cluster_groups_1idx),1), ...
        'unsorted',false(max(cluster_groups_1idx),1));
        
end

% Templates/clusters are 0-indexed, fix that
spike_templates_1idx = spike_templates + 1;
spike_clusters_1idx = spike_clusters + 1;

% Group spike times by cluster for easy access
cluster_spike_times = cell(max(spike_clusters_1idx),1);
for curr_cluster = 1:max(spike_clusters_1idx)
    cluster_spike_times{curr_cluster} = ...
        double(spike_times(spike_clusters_1idx == curr_cluster))/ephys_sample_rate;
end

%% Get time conversions and trial information

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
    (photodiode(1:end-1) >= 1.5 & photodiode(2:end) < 1.5))+1;
photodiode_sec = photodiode_flip/timeline_sample_rate;

% Get stimulus presentation times relative to ephys
n_trials = length(block.trial);
stim_times_block = [block.trial.stimulusCueStartedTime];
stim_samples_photodiode = photodiode_flip(arrayfun(@(x) find(photodiode_sec >= ...
    stim_times_block(x),1),1:n_trials));
stim_samples_ephys = stim_samples_photodiode*timeline_ephys_scale + timeline_ephys_offset;
stim_times_ephys = stim_samples_ephys/ephys_sample_rate;

% Get stimulus conditions
trial_condition = cell2mat(cellfun(@(x) x.visCueContrast,{block.trial(:).condition},'uni',false))';

%% Get cluster information

% Get most used template for each cluster
cluster_template = nan(max(spike_clusters_1idx),1);
for curr_clust = unique(spike_clusters_1idx)';
    cluster_template(curr_clust) = ...
        mode(spike_templates_1idx(spike_clusters_1idx == curr_clust));
end
used_clusters = ~isnan(cluster_template);

% Location (center of mass) of each template (from bottom of probe??)
template_abs = permute(max(abs(templates),[],2),[3,1,2]);

template_com_x = (template_abs'*channel_positions(:,1))./sum(template_abs,1)';
%template_com_y = (template_abs'*channel_positions(:,2))./sum(template_abs,1)';
[~,max_channel_idx] =  max(template_abs,[],1); % COM gets too much crap
template_com_y = channel_positions(max_channel_idx,2);

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

% Get duration (trough-to-peak, Bartho JNeurophys 2004)
templates_max_signfix = bsxfun(@times,templates_max, ...
    sign(abs(min(templates_max,[],2)) - abs(max(templates_max,[],2))));

[~,waveform_trough] = min(templates_max,[],2);
[~,waveform_peak_rel] = arrayfun(@(x) ...
    max(templates_max(x,waveform_trough(x):end),[],2), ...
    transpose(1:size(templates_max,1)));
waveform_peak = waveform_peak_rel + waveform_trough;

template_duration = waveform_peak - waveform_trough;

cluster_duration = nan(max(spike_clusters),1);
cluster_duration(used_clusters) = template_duration(cluster_template(used_clusters));

% Group waveforms by narrow/broad spike duration
waveform_narrow = cluster_duration <= 11;
waveform_broad = cluster_duration > 11;


%% Store relevant spike data in master structure
spikes = struct;

spikes.sorted = sorted;

spikes.sample_rate = ephys_sample_rate;
spikes.spike_samples = spike_times;
spikes.spike_times = spike_times_sec;
spikes.spike_templates = spike_templates_1idx;
spikes.spike_clusters = spike_clusters_1idx;
spikes.channel_positions = channel_positions;

spikes.template.waveform = templates_max;
spikes.template.depth = template_com_y;
spikes.template.fwhm = template_fwhm; 
spikes.template.duration = template_duration;

spikes.cluster.spike_times = cluster_spike_times;
spikes.cluster.group = cluster_groups;
spikes.cluster.template = cluster_template;
spikes.cluster.waveform = cluster_waveform;
spikes.cluster.depth = cluster_com_y;
spikes.cluster.fwhm = cluster_fwhm;
spikes.cluster.duration = cluster_duration;

spikes.cluster.narrow = waveform_narrow;
spikes.cluster.broad = waveform_broad;

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
xpr.bhv.correct = [block.trial(:).feedbackType] == 1;

xpr.timeline2ephys = [timeline_ephys_scale,timeline_ephys_offset];





