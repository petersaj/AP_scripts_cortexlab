%% Define animal and day

animal = '80';
day = '20160527';

% experiment is always the last one of the day
day_dash = [day(1:4) '-' day(5:6) '-' day(7:8)];
expInfo_path = ['\\zserver.cortexlab.net\Data\expInfo\80\' day_dash];
expInfo_dir = dir(expInfo_path);
expInfo_expts = sort({expInfo_dir.name});
experiment = str2num(expInfo_expts{end});

%% Load experiment info

% Find filenames for behavior/input
timeline_filename = get_cortexlab_filename(animal,day,experiment,'timeline','dash');
parameters_filename = get_cortexlab_filename(animal,day,experiment,'parameters','dash');
block_filename = get_cortexlab_filename(animal,day,experiment,'block','dash');

% Load behavior/input
load(timeline_filename);
load(parameters_filename);
load(block_filename);

% Get acquisition live in timeline (used for sync)
acqLive_idx = strcmp({Timeline.hw.inputs.name}, 'acqLiveEcho');
thresh = max(Timeline.rawDAQData(:,acqLive_idx))/2;
acqLive_timeline = Timeline.rawDAQTimestamps( ...
    ((Timeline.rawDAQData(1:end-1,acqLive_idx) < thresh) & ...
    (Timeline.rawDAQData(2:end,acqLive_idx) > thresh)) | ...
    ((Timeline.rawDAQData(1:end-1,acqLive_idx) > thresh) & ...
    (Timeline.rawDAQData(2:end,acqLive_idx) < thresh)));

% Get photodiode signal from timeline
photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
thresh = max(Timeline.rawDAQData(:,photodiode_idx))/2;
photodiode_timeline = Timeline.rawDAQTimestamps( ...
    ((Timeline.rawDAQData(1:end-1,photodiode_idx) < thresh) & ...
    (Timeline.rawDAQData(2:end,photodiode_idx) > thresh)) | ...
    ((Timeline.rawDAQData(1:end-1,photodiode_idx) > thresh) & ...
    (Timeline.rawDAQData(2:end,photodiode_idx) < thresh)));

% Get stimulus presentation times (leave out the last trial - incomplete)
n_trials = length(block.trial)-1;
stim_times_block = [block.trial.stimulusCueStartedTime];
stim_times_timeline = photodiode_timeline(arrayfun(@(x) find(photodiode_timeline >= ...
    stim_times_block(x),1),1:n_trials));

% Get stimulus conditions
trial_condition = cell2mat(cellfun(@(x) x.visCueContrast,{block.trial(1:n_trials).condition},'uni',false))';
trial_correct = [block.trial(1:n_trials).feedbackType] == 1;

% Get percent correct by condition
all_conditions = diff(trial_condition,[],2);
[plot_conditions,condition_correct] =  grpstats(trial_correct,all_conditions,{'gname','mean'});
plot_conditions = cellfun(@(x) str2num(x),plot_conditions);


%% Load electrophysiology

flipped_banks = true; % plugged the banks in in reverse 

data_path = ['\\basket.cortexlab.net\data\ajpeters\' animal filesep day];

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

% Get the spike times in timeline time
spike_times_timeline = AP_clock_fix(spike_times,sync.timestamps,acqLive_timeline);

% Load clusters, if they exist
cluster_filename = [data_path filesep 'cluster_groups.csv'];
if exist(cluster_filename,'file')
    fid = fopen(cluster_filename);
    cluster_groups = textscan(fid,'%d%s','HeaderLines',1);
    fclose(fid);
end

% Eliminate spikes that were classified as not "good"

% Saftey check: if this variable exists, don't do it again
if exist('cluster_groups','var') && ~exist('good_templates','var')
    disp('Removing non-good templates')
    good_templates = uint32(cluster_groups{1}(strcmp(cluster_groups{2},'good')));
    good_spike_idx = ismember(spike_templates,good_templates);
    
    spike_times = spike_times(good_spike_idx);
    spike_templates = spike_templates(good_spike_idx);
    template_amplitudes = template_amplitudes(good_spike_idx);
    spike_times_timeline = spike_times_timeline(good_spike_idx);
else
    disp('Good templates already identified')
end

% Flip channel map and positions if banks are reversed
if flipped_banks
    channel_map = [channel_map(61:end);channel_map(1:60)];
    channel_positions = [channel_positions(61:end,:);channel_positions(1:60,:)];
    templates = cat(3,templates(:,:,61:end),templates(:,:,1:60));
end


% Default channel map/positions are from end: make from surface
channel_positions(:,2) = max(channel_positions(:,2)) - channel_positions(:,2);

% Load LFP
n_channels = str2num(header.n_channels);
lfp_filename = [data_path filesep 'lfp.dat'];
fid = fopen(lfp_filename);
lfp_all = fread(fid,[n_channels,inf],'int16');
fclose(fid);
% eliminate non-connected channels and sort by position
lfp = lfp_all(flipud(channel_map+1),:);
% get time of LFP sample points (NOTE: this is messy, based off of sample
% rate and knowing what kwik2dat does, not sure how accurate)
sample_rate = str2num(header.sample_rate);
lfp_cutoff = str2num(header.lfp_cutoff);
lfp_downsamp = (sample_rate/lfp_cutoff)/2;
lfp_t = ([1:size(lfp,2)]*lfp_downsamp)/sample_rate;
lfp_t_timeline = AP_clock_fix(lfp_t,sync.timestamps,acqLive_timeline);

% Get the depths of each template 
% (by COM: this gives totally wonky answers because of artifacts maybe?)
%[spikeAmps, spike_depths, template_depths, tempAmps, tempsUnW, templateDuration, waveforms] = ...
%    templatePositionsAmplitudes(templates,winv,channel_positions(:,2),spike_templates,template_amplitudes);

% (by max waveform channel)
template_abs = permute(max(abs(templates),[],2),[3,1,2]);
[~,max_channel_idx] =  max(template_abs,[],1);
template_depths = channel_positions(max_channel_idx,2);


% Get each spike's depth
spike_depths = template_depths(spike_templates+1);

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

%% MAKE PLOTS

%% PSTH - superficial

start_depth = 0;
end_depth = 400;

use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > start_depth & template_depths < end_depth)-1));
use_spike_templates = spike_templates(ismember(spike_templates,find(template_depths > start_depth & template_depths < end_depth)-1));

align_times = stim_times_timeline;

all_conditions = diff(trial_condition,[],2);
unique_conditions = unique(all_conditions);

raster_window = [-0.5,0.8];

% Plot PSTH of whole population per stimulus
psth_bin_size = 0.001;
all_bins = raster_window(1):psth_bin_size:raster_window(2);
bins = diff(all_bins)/2+all_bins(1:end-1);

psth = nan(length(unique_conditions),diff(raster_window)/psth_bin_size);
psth_sem = nan(length(unique_conditions),diff(raster_window)/psth_bin_size);
for curr_stim = 1:length(unique_conditions);

    use_align = align_times((all_conditions == unique_conditions(curr_stim)))';
    
    curr_psth = nan(length(use_align),diff(raster_window)/psth_bin_size);
    for curr_trial = 1:length(use_align)
        curr_psth(curr_trial,:) = histcounts(use_spikes, ...
            [use_align(curr_trial) + raster_window(1):psth_bin_size:...
            use_align(curr_trial) + raster_window(2)]);
    end
        
    psth(curr_stim,:) = nanmean(curr_psth,1);
    psth_sem(curr_stim,:) = nanstd(curr_psth,[],1)./sqrt(sum(~isnan(curr_psth),1));

end

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw.*psth_bin_size);
psth_smooth = conv2(psth, smWin, 'same');
plot_spacing = 500;
psth_plot = bsxfun(@plus,psth_smooth,transpose(1:size(psth_smooth,1))*plot_spacing);
psth_plot(all(isnan(psth_smooth),2),:) = NaN;

psth_sem_smooth = conv2(psth_sem, smWin, 'same');

figure; hold on;

% Draw shaded error bars
for curr_stim = 1:size(psth_plot,1)
    fill([bins,fliplr(bins)], ...
        [psth_plot(curr_stim,:) + psth_sem_smooth(curr_stim,:), ...
        fliplr(psth_plot(curr_stim,:) - psth_sem_smooth(curr_stim,:))],[0.5,0.5,0.5])
end

% Plot the means
plot(bins,psth_plot','k','linewidth',2);
line([0,0],ylim,'linestyle','--','color','k');
ylabel('Contrast difference');
xlabel('Time from stimulus onset')
set(gca,'YTick',transpose(1:size(psth_smooth,1))*plot_spacing)
set(gca,'YTickLabel',unique_conditions);

title(['MUA, ' num2str(start_depth) '-' num2str(end_depth) '\mum'])

%% PSTH - deep

start_depth = 400;
end_depth = 800;

use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > start_depth & template_depths < end_depth)-1));
use_spike_templates = spike_templates(ismember(spike_templates,find(template_depths > start_depth & template_depths < end_depth)-1));

align_times = stim_times_timeline;

all_conditions = diff(trial_condition,[],2);
unique_conditions = unique(all_conditions);

raster_window = [-0.5,0.8];

% Plot PSTH of whole population per stimulus
psth_bin_size = 0.001;
all_bins = raster_window(1):psth_bin_size:raster_window(2);
bins = diff(all_bins)/2+all_bins(1:end-1);

psth = nan(length(unique_conditions),diff(raster_window)/psth_bin_size);
psth_sem = nan(length(unique_conditions),diff(raster_window)/psth_bin_size);
for curr_stim = 1:length(unique_conditions);
     
    use_align = align_times((all_conditions == unique_conditions(curr_stim)))';
    
    curr_psth = nan(length(use_align),diff(raster_window)/psth_bin_size);
    for curr_trial = 1:length(use_align)
        curr_psth(curr_trial,:) = histcounts(use_spikes, ...
            [use_align(curr_trial) + raster_window(1):psth_bin_size:...
            use_align(curr_trial) + raster_window(2)]);
    end
        
    psth(curr_stim,:) = nanmean(curr_psth,1);
    psth_sem(curr_stim,:) = nanstd(curr_psth,[],1)./sqrt(sum(~isnan(curr_psth),1));

end

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw.*psth_bin_size);
psth_smooth = conv2(psth, smWin, 'same');
plot_spacing = 500;
psth_plot = bsxfun(@plus,psth_smooth,transpose(1:size(psth_smooth,1))*plot_spacing);
psth_plot(all(isnan(psth_smooth),2),:) = NaN;

psth_sem_smooth = conv2(psth_sem, smWin, 'same');

figure; hold on;

% Draw shaded error bars
for curr_stim = 1:size(psth_plot,1)
    fill([bins,fliplr(bins)], ...
        [psth_plot(curr_stim,:) + psth_sem_smooth(curr_stim,:), ...
        fliplr(psth_plot(curr_stim,:) - psth_sem_smooth(curr_stim,:))],[0.5,0.5,0.5])
end

% Plot the means
plot(bins,psth_plot','k','linewidth',2);
line([0,0],ylim,'linestyle','--','color','k');
ylabel('Contrast difference');
xlabel('Time from stimulus onset')
set(gca,'YTick',transpose(1:size(psth_smooth,1))*plot_spacing)
set(gca,'YTickLabel',unique_conditions);

title(['MUA, ' num2str(start_depth) '-' num2str(end_depth) '\mum'])

%% Stim-triggered LFP (superficial and deep)

% group LFP by depth
lfp_use_depths = [0 400 800];
channel_depth_grp = discretize(sort(channel_positions(:,2),'ascend'),lfp_use_depths);
lfp_depth_mean = grpstats(lfp,channel_depth_grp);

all_conditions = diff(trial_condition,[],2);
unique_conditions = unique(all_conditions);

align_times = stim_times_timeline;

% Stimulus-triggered LFP by stim
lfp_window = [-0.5,0.8];
t_space = 0.001;
lfp_stim_mean = nan(length(unique_conditions),1+diff(lfp_window)/t_space,size(lfp_depth_mean,1));
lfp_stim_sem = nan(length(unique_conditions),1+diff(lfp_window)/t_space,size(lfp_depth_mean,1));

for curr_stim = 1:length(unique_conditions);
       
    use_align = align_times((all_conditions == unique_conditions(curr_stim)))';
    use_align_t = bsxfun(@plus,use_align,lfp_window(1):t_space:lfp_window(2));

    curr_lfp_stim = interp1(lfp_t_timeline,lfp_depth_mean',use_align_t);
    
    lfp_stim_mean(curr_stim,:,:) = nanmedian(curr_lfp_stim,1);
    lfp_stim_sem(curr_stim,:,:) = nanstd(curr_lfp_stim,[],1)./sqrt(sum(~isnan(curr_lfp_stim),1));
    
end

plot_spacing = 2000;
lfp_stim_mean_plot = bsxfun(@plus,lfp_stim_mean,transpose(1:size(lfp_stim_mean,1))*plot_spacing);

t = lfp_window(1):t_space:lfp_window(2);

% Plot superficial
figure; hold on;

plot_lfp = 1;
% Shaded error bars
for curr_stim = 1:size(lfp_stim_mean,1)
    fill([t,fliplr(t)], ...
        [lfp_stim_mean_plot(curr_stim,:,plot_lfp) + lfp_stim_sem(curr_stim,:,plot_lfp), ...
        fliplr(lfp_stim_mean_plot(curr_stim,:,plot_lfp) - lfp_stim_sem(curr_stim,:,plot_lfp))],[0.5,0.5,0.5])
end
% Means
plot(t,lfp_stim_mean_plot(:,:,plot_lfp)','k','linewidth',2)
line([0,0],ylim,'linestyle','--','color','k');
ylabel('Contrast difference');
xlabel('Time from stimulus onset')
set(gca,'YTick',transpose(1:length(unique_conditions))*plot_spacing)
set(gca,'YTickLabel',unique_conditions);

title(['LFP, ' num2str(lfp_use_depths(1)) '-' num2str(lfp_use_depths(2)) '\mum'])

% Plot deep
figure; hold on;

plot_lfp = 2;
% Shaded error bars
for curr_stim = 1:size(lfp_stim_mean,1)
    fill([t,fliplr(t)], ...
        [lfp_stim_mean_plot(curr_stim,:,plot_lfp) + lfp_stim_sem(curr_stim,:,plot_lfp), ...
        fliplr(lfp_stim_mean_plot(curr_stim,:,plot_lfp) - lfp_stim_sem(curr_stim,:,plot_lfp))],[0.5,0.5,0.5])
end
% Means
plot(t,lfp_stim_mean_plot(:,:,plot_lfp)','k','linewidth',2)
line([0,0],ylim,'linestyle','--','color','k');
ylabel('Contrast difference');
xlabel('Time from stimulus onset')
set(gca,'YTick',transpose(1:length(unique_conditions))*plot_spacing)
set(gca,'YTickLabel',unique_conditions);

title(['LFP, ' num2str(lfp_use_depths(2)) '-' num2str(lfp_use_depths(3)) '\mum'])
