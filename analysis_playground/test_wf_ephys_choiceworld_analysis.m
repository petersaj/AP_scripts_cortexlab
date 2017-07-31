%% Raster plots to choiceworld conditions

stimIDs = signals_events.trialSideValues.*signals_events.trialContrastValues;

use_spikes_idx = ismember(spike_templates,find(templateDepths >= 1500 & templateDepths <= 2000));
use_spikes = spike_times_timeline(use_spikes_idx);

raster_window = [-0.5,2];
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
AP_stackplot(stim_psth_hit_smooth(:,20:end-20)',bins(20:end-20),trace_spacing,false,'k',unique(stimIDs));
AP_stackplot(stim_psth_miss_smooth(:,20:end-20)',bins(20:end-20),trace_spacing,false,'r',unique(stimIDs));
xlabel('Time from stim onset')
ylabel('Population spikes (by stim)');
line([0,0],ylim,'linestyle','--','color','k');


%% Raster plots for left vs. right stim, choose left vs. right stim (by depth)

% Group multiunit by depth
n_depth_groups = 8;
depth_group_edges = linspace(1300,max(channel_positions(:,2)),n_depth_groups+1);
%depth_group_edges = [0 1300];
depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);

depth_group = discretize(spikeDepths,depth_group_edges);


raster_window = [-0.5,2];
psth_bin_size = 0.001;
smooth_size = 100;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

psth_right_hit = nan(n_depth_groups,diff(raster_window)/psth_bin_size);
psth_right_miss = nan(n_depth_groups,diff(raster_window)/psth_bin_size);
psth_left_hit = nan(n_depth_groups,diff(raster_window)/psth_bin_size);
psth_left_miss = nan(n_depth_groups,diff(raster_window)/psth_bin_size);

for curr_depth = 1:n_depth_groups
    
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    
    use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 1;
    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,signals_events.stimOnTimes(use_trials),raster_window,psth_bin_size);
    psth_smooth = conv2(psth,smWin,'same');
    psth_right_hit(curr_depth,:) = psth_smooth;
    
    use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 0;
    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,signals_events.stimOnTimes(use_trials),raster_window,psth_bin_size);
    psth_smooth = conv2(psth,smWin,'same');
    psth_right_miss(curr_depth,:) = psth_smooth;   
        
    use_trials = signals_events.trialSideValues == -1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 1;
    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,signals_events.stimOnTimes(use_trials),raster_window,psth_bin_size);
    psth_smooth = conv2(psth,smWin,'same');
    psth_left_hit(curr_depth,:) = psth_smooth;
    
    use_trials = signals_events.trialSideValues == -1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 0;
    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,signals_events.stimOnTimes(use_trials),raster_window,psth_bin_size);
    psth_smooth = conv2(psth,smWin,'same');
    psth_left_miss(curr_depth,:) = psth_smooth;
    
end

figure; hold on;
trace_spacing = max([max(psth_right_hit(:)),max(psth_right_miss(:)),max(psth_left_hit(:)),max(psth_left_hit(:))]);
p_rh = AP_stackplot(psth_right_hit',bins,trace_spacing,[],'k',depth_group_centers);
p_rm = AP_stackplot(psth_right_miss',bins,trace_spacing,[],'r',depth_group_centers);
p_lh = AP_stackplot(psth_left_hit',bins,trace_spacing,[],'b',depth_group_centers);
p_lm = AP_stackplot(psth_left_miss',bins,trace_spacing,[],'m',depth_group_centers);

line([0,0],ylim,'color','k','linestyle','--');

legend([p_rh(1),p_rm(1),p_lh(1),p_lm(1)],{'Stim right hit','Stim right miss','Stim left hit','Stim left miss'})
xlabel('Time from stim');
ylabel('Depth (\mum)');


%% Raster plots for left vs. right stim, choose left vs. right stim (errorbars)

stimIDs = signals_events.trialSideValues.*signals_events.trialContrastValues;

use_spikes_idx = ismember(spike_templates,find(templateDepths >= 1300 & templateDepths <= 2500));
use_spikes = spike_times_timeline(use_spikes_idx);

raster_window = [-0.5,2];
psth_bin_size = 0.02;

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

% PSTH by stim/move left/right
figure;

subplot(1,3,1); hold on;

use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 1;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(use_spikes,signals_events.stimOnTimes(use_trials),raster_window,psth_bin_size);
AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'k');

use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 0;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(use_spikes,signals_events.stimOnTimes(use_trials),raster_window,psth_bin_size);
AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'r');

use_trials = signals_events.trialSideValues == -1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 1;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(use_spikes,signals_events.stimOnTimes(use_trials),raster_window,psth_bin_size);
AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'b');

use_trials = signals_events.trialSideValues == -1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 0;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(use_spikes,signals_events.stimOnTimes(use_trials),raster_window,psth_bin_size);
AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'m');

legend({'Right hit','Right miss','Left hit','Left miss'})
xlim([bins(1),bins(end)])
xlabel('Time from stim');
ylabel('Spikes');

subplot(1,3,2); hold on;

use_trials = (signals_events.trialSideValues == 1 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == -1 & signals_events.hitValues == 0);
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(use_spikes,signals_events.stimOnTimes(use_trials),raster_window,psth_bin_size);
AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'k');

use_trials = (signals_events.trialSideValues == -1 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == 1 & signals_events.hitValues == 0);
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(use_spikes,signals_events.stimOnTimes(use_trials),raster_window,psth_bin_size);
AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'r');

legend({'Move left','Move right'});
xlim([bins(1),bins(end)])
xlabel('Time from stim');
ylabel('Spikes');


subplot(1,3,3); hold on;

use_trials = (signals_events.trialSideValues == 1 & signals_events.trialContrastValues == 0 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == -1 & signals_events.trialContrastValues == 0 & signals_events.hitValues == 0);
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(use_spikes,signals_events.stimOnTimes(use_trials),raster_window,psth_bin_size);
AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'k');

use_trials = (signals_events.trialSideValues == -1 & signals_events.trialContrastValues == 0 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == 1 & signals_events.trialContrastValues == 0 & signals_events.hitValues == 0);
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(use_spikes,signals_events.stimOnTimes(use_trials),raster_window,psth_bin_size);
AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'r');

legend({'Move left 0 contrast','Move right 0 contrast'});
xlim([bins(1),bins(end)])
xlabel('Time from stim');
ylabel('Spikes');




