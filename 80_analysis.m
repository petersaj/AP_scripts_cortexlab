%% Define animal and day

animal = '80';
day = '20160602';

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

%% PSTHs

align_times = [block.trial.interactiveStartedTime];
%align_times = [block.trial.feedbackStartedTime];
%align_times = cellfun(@(x) x(1),{block.trial.interactiveMovementTime});
%align_times = stim_times_timeline;

all_conditions = diff(trial_condition,[],2);
unique_conditions = unique(all_conditions);

raster_window = [-1,2];
psthViewer(spike_times_timeline,spike_templates, ...
    align_times,raster_window,all_conditions);

% Plot PSTH of whole population per stimulus (correct vs. incorrect)
psth_bin_size = 0.001;
psth_correct = nan(length(unique_conditions),diff(raster_window)/psth_bin_size);
psth_incorrect = nan(length(unique_conditions),diff(raster_window)/psth_bin_size);
for curr_stim = 1:length(unique_conditions);
%     [psth(curr_stim,:),bins,rasterX,rasterY,spikeCounts] = ...
%         psthRasterAndCounts(spike_times_timeline, ...
%         stim_times_timeline(plot_conditions == unique_conditions(curr_stim))', ...
%         raster_window, psth_bin_size);
%     
    use_align = align_times((all_conditions == unique_conditions(curr_stim)) & trial_correct')';
    if ~isempty(use_align)
    [psth_correct(curr_stim,:),bins,rasterX,rasterY,spikeCounts] = psthAndBA( ...
        spike_times_timeline,use_align, ...
        raster_window, psth_bin_size);
    end
    
    use_align = align_times((all_conditions == unique_conditions(curr_stim)) & ~trial_correct')';
    if ~isempty(use_align)
    [psth_incorrect(curr_stim,:),bins,rasterX,rasterY,spikeCounts] = psthAndBA( ...
        spike_times_timeline,use_align, ...
        raster_window, psth_bin_size);
    end

end

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
psth_correct_smooth = conv2(psth_correct, smWin, 'same');
psth_incorrect_smooth = conv2(psth_incorrect, smWin, 'same');
plot_spacing = 1;
plot_min = min([psth_correct_smooth(:);psth_incorrect_smooth(:)]);
plot_max = max([psth_correct_smooth(:);psth_incorrect_smooth(:)]);
psth_correct_plot = bsxfun(@plus,mat2gray(psth_correct_smooth,[plot_min,plot_max]),transpose(1:size(psth_correct_smooth,1))*plot_spacing);
psth_incorrect_plot = bsxfun(@plus,mat2gray(psth_incorrect_smooth,[plot_min,plot_max]),transpose(1:size(psth_incorrect_smooth,1))*plot_spacing);

psth_correct_plot(all(isnan(psth_correct_smooth),2),:) = NaN;
psth_incorrect_plot(all(isnan(psth_incorrect_smooth),2),:) = NaN;

figure; hold on;
plot(bins,psth_correct_plot','k','linewidth',2);
plot(bins,psth_incorrect_plot','r','linewidth',2);
line([0,0],ylim,'linestyle','--','color','k');
ylabel('Contrast difference');
xlabel('Time from interactive onset')
legend({'Correct','Incorrect'});
set(gca,'YTickLabel',unique_conditions);
title([day ': ' num2str(round(nanmean(trial_correct)*100)) '%']);


%% Wheel movement

%align_times = [block.trial.interactiveStartedTime];
%align_times = [block.trial.feedbackStartedTime];
align_times = cellfun(@(x) x(1),{block.trial.interactiveMovementTime});

surround_time = [-1,2];
surround_samples = surround_time*Timeline.hw.daqSampleRate;

% rotary encoder goes backwards in time?? sort it I guess...
[~,rotary_sort_idx] = sort(block.inputSensorPositionTimes);

wheel_times = block.inputSensorPositionTimes(rotary_sort_idx);
wheel_positions = block.inputSensorPositions(rotary_sort_idx);

wheel_align_times = cellfun(@(x) ...
    wheel_times(wheel_times > x + surround_time(1) & ...
    wheel_times < x + surround_time(2)) - x,num2cell(align_times),'uni',false);
wheel_align_position = cellfun(@(x) ...
    wheel_positions(wheel_times > x + surround_time(1) & ...
    wheel_times < x + surround_time(2)),num2cell(align_times),'uni',false);

use_trials = trial_correct & (diff(trial_condition,[],2) == 0.12)';

figure; hold on
for trial = find(use_trials)
    plot(wheel_align_times{trial},wheel_align_position{trial} - wheel_align_position{trial}(1),'k')
end

%% Get PSTHs by spike width (only for one contrast)

% Get waveform information
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = ...
    templatePositionsAmplitudes(templates,winv,channel_positions(:,2),spike_templates,template_amplitudes);

% Group by waveform duration
duration_cutoff = round(0.5*30);
narrow_templates = find(templateDuration <= duration_cutoff);
narrow_spikes = ismember(spike_templates,narrow_templates);

%align_times = [block.trial.interactiveStartedTime];
%align_times = [block.trial.feedbackStartedTime];
align_times = cellfun(@(x) x(end),{block.trial.interactiveMovementTime});
%align_times = stim_times_timeline;

all_conditions = diff(trial_condition,[],2);
unique_conditions = unique(all_conditions);

raster_window = [-1,2];

% Plot PSTH of whole population per stimulus (correct vs. incorrect)
psth_bin_size = 0.001;
curr_stim = length(unique_conditions);

use_align = align_times((all_conditions == unique_conditions(curr_stim)) & trial_correct')';

[psth_wide,bins,rasterX,rasterY,spikeCounts] = psthAndBA( ...
    spike_times_timeline(~narrow_spikes),use_align, ...
    raster_window, psth_bin_size);

[psth_narrow,bins,rasterX,rasterY,spikeCounts] = psthAndBA( ...
    spike_times_timeline(narrow_spikes),use_align, ...
    raster_window, psth_bin_size);

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
psth_wide_smooth = conv2(psth_wide, smWin, 'same');
psth_narrow_smooth = conv2(psth_narrow, smWin, 'same');

figure; hold on;
plot(bins,zscore(psth_wide_smooth),'k','linewidth',2);
plot(bins,zscore(psth_narrow_smooth),'r','linewidth',2);
line([0,0],ylim,'linestyle','--','color','k');
ylabel('Normalized PSTH');
xlabel('Time from interactive onset')
legend({'Wide waveforms','Narrow waveforms'});
set(gca,'YTickLabel',unique_conditions);
title([day ': ' num2str(round(nanmean(trial_correct)*100)) '%']);

%% Population spikes vs. cumulative wheel position

start_times = [block.trial(1:n_trials).interactiveStartedTime];
stop_times  = [block.trial(1:n_trials).feedbackStartedTime];

% rotary encoder goes backwards in time?? sort it I guess...
[~,rotary_sort_idx] = sort(block.inputSensorPositionTimes);

wheel_times = block.inputSensorPositionTimes(rotary_sort_idx);
wheel_positions = block.inputSensorPositions(rotary_sort_idx);

wheel_times_trial = cell(length(start_times),1);
wheel_positions_trial_rel = cell(length(start_times),1);
for curr_trial = 1:length(start_times);
    
    % Get Timeline times within this trial
    curr_times = Timeline.rawDAQTimestamps( ...
        Timeline.rawDAQTimestamps > start_times(curr_trial) & ...
        Timeline.rawDAQTimestamps < stop_times(curr_trial));
    
    % Get wheel times and position within this trial
    curr_wheel_idx = wheel_times > start_times(curr_trial) & ...
        wheel_times < stop_times(curr_trial);
    curr_wheel_times = wheel_times(curr_wheel_idx);
    curr_wheel_positions = wheel_positions(curr_wheel_idx);
    
    % Interpolate wheel across times
    curr_wheel_positions_interp = interp1(curr_wheel_times, ...
        curr_wheel_positions,curr_times);
    
    % Set beginning/end NaNs to closest values
    start_val = curr_wheel_positions_interp(find(~isnan( ...
        curr_wheel_positions_interp),1));
    end_val = curr_wheel_positions_interp(find(~isnan( ...
        curr_wheel_positions_interp),1,'last'));
    
    curr_wheel_positions_interp(1:find(~isnan( ...
        curr_wheel_positions_interp),1)) = start_val;
    curr_wheel_positions_interp(end:-1:find(~isnan( ...
        curr_wheel_positions_interp),1,'last')) = end_val;
    
    % Store
    wheel_times_trial{curr_trial} = curr_times;
    wheel_positions_trial_rel{curr_trial} = curr_wheel_positions_interp;
    
    disp(curr_trial);
end

% Fix the wheel positions based on contrast and correctness
trial_contrast = diff(trial_condition,[],2);
out_center_distance = median(cellfun(@(x) abs(x(end) - x(1)), ...
    wheel_positions_trial_rel(trial_contrast ~= 0 & trial_correct')));
wheel_positions_trial = cellfun(@(wheel,contrast) (wheel-wheel(1)) + ...
    out_center_distance*sign(contrast),wheel_positions_trial_rel, ...
    num2cell(trial_contrast),'uni',false);
  
% Bin spike times according to wheel position
n_bins = 20;
wheel_bin_edges= [linspace(min(horzcat(wheel_positions_trial{:})), ...
    max(horzcat(wheel_positions_trial{:})),n_bins)];
wheel_bin_centers = wheel_bin_edges(1:end-1) + diff(wheel_bin_edges)/2;
wheel_bins = cellfun(@(x) discretize(x,wheel_bin_edges),wheel_positions_trial,'uni',false);

wheel_bin_spike_count = zeros(n_bins-1,1);
wheel_bin_n = zeros(n_bins-1,1);
for curr_trial = find(trial_contrast == -0.5 & trial_correct')'
    
    curr_spike_count = histcounts(spike_times_timeline,wheel_times_trial{curr_trial});
    
    for curr_bin = 1:n_bins-1
        wheel_bin_spike_count(curr_bin) = ...
            wheel_bin_spike_count(curr_bin) + ...
            sum(curr_spike_count(wheel_bins{curr_trial}(1:end-1) == curr_bin));
        
        wheel_bin_n(curr_bin) = ...
            wheel_bin_n(curr_bin) + ...
            sum(wheel_bins{curr_trial}(1:end-1) == curr_bin);
    end
    
end

wheel_bin_spike_mean = wheel_bin_spike_count./wheel_bin_n;
figure;plot(wheel_bin_centers,wheel_bin_spike_mean,'k');












