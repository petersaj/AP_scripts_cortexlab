%% PSTH to choiceworld conditions

stimIDs = signals_events.trialSideValues.*signals_events.trialContrastValues;
% stimIDs = discretize(stimIDs,[-Inf,-0.125,-0.01,0.01,0.25,Inf],[-2,-1,0,1,2]);

use_spikes_idx = ismember(spike_templates,find(templateDepths >= str_depth(1) & templateDepths <= str_depth(1)+100));
% use_spikes_idx = ismember(spike_templates,find(templateDepths >= 500 & templateDepths <= 800));
% use_spikes_idx = ismember(spike_templates,intersect(find(templateDepths >= 500 & templateDepths <= 1500),find(msn)));

use_spikes = spike_times_timeline(use_spikes_idx);

raster_window = [-0.5,2.5];
psth_bin_size = 0.001;

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

% PSTH by stim condition
stim_psth_hit = nan(length(unique(stimIDs)),diff(raster_window)/psth_bin_size);
stim_psth_miss = nan(length(unique(stimIDs)),diff(raster_window)/psth_bin_size);

unique_stims = unique(stimIDs);
for curr_stim_idx = 1:length(unique_stims);
    
    min_trials = 5;
    
    curr_trials_hit = (stimIDs == unique_stims(curr_stim_idx)) & ...
        signals_events.hitValues == 1;
    if sum(curr_trials_hit) > min_trials
        [psth_hit,bins] = psthAndBA( ...
            use_spikes,stimOn_times(curr_trials_hit), ...
            raster_window, psth_bin_size);
        stim_psth_hit(curr_stim_idx,:) = psth_hit;
    end
    
    curr_trials_miss = (stimIDs == unique_stims(curr_stim_idx)) & ...
        signals_events.hitValues == 0;
    if sum(curr_trials_miss) > min_trials
        [psth_miss,bins] = psthAndBA( ...
            use_spikes,stimOn_times(curr_trials_miss), ...
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

median_move_time = nanmedian(wheel_move_time' - stimOn_times);
median_reward_time = nanmedian(reward_t_timeline - stimOn_times(signals_events.hitValues == 1)');

line([0,0],ylim,'color','k','linestyle','--');
line([median_move_time,median_move_time],ylim,'color','r','linestyle','--');
line([median_reward_time,median_reward_time],ylim,'color','b','linestyle','--');

% Plot the hit PSTHs on top of each other
figure; hold on
set(gca,'ColorOrder',copper(length(unique(stimIDs))));
plot(bins(20:end-20),stim_psth_hit_smooth(:,20:end-20)','linewidth',2)

%% Plot contrast response curves response by depth

stimIDs = signals_events.trialSideValues.*signals_events.trialContrastValues;

% Group multiunit by depth
n_depth_groups = 6;
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);

depth_group = discretize(spikeDepths,depth_group_edges);

raster_window = [-0.2,1];
psth_bin_size = 0.001;
smooth_size = 100;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

% Contrast response function
unique_stims = unique(stimIDs);
contrast_response = nan(n_depth_groups,length(unique_stims));
for curr_depth = 1:n_depth_groups
    
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
%     curr_spike_times = spike_times_timeline((depth_group == curr_depth) & ...
%         ismember(spike_templates,find(msn)));
        
    stim_psth_hit = nan(length(unique(stimIDs)),diff(raster_window)/psth_bin_size);
    for curr_stim_idx = 1:length(unique_stims);
        
        min_trials = 5;
        
        curr_trials_hit = (stimIDs == unique_stims(curr_stim_idx)) & ...
            signals_events.hitValues == 1;
        
        if sum(curr_trials_hit) > min_trials
            [psth_hit,bins] = psthAndBA( ...
                curr_spike_times,stimOn_times(curr_trials_hit), ...
                raster_window, psth_bin_size);
            stim_psth_hit(curr_stim_idx,:) = psth_hit;
        end
        
    end
    
    stim_psth_hit_smooth = conv2(stim_psth_hit,smWin,'same');
    stim_psth_hit_smooth_max = max(stim_psth_hit_smooth,[],2);
    contrast_response(curr_depth,:) = stim_psth_hit_smooth_max;   
    
end

contrast_response_norm = bsxfun(@rdivide,bsxfun(@minus,contrast_response,min(contrast_response,[],2)),min(contrast_response,[],2));
figure; hold on;
set(gca,'ColorOrder',copper(n_depth_groups));
plot(unique_stims,contrast_response_norm','linewidth',2)
xlabel('Contrast');
ylabel('Normalized response');
legend(cellfun(@(x) ['Depth ' num2str(x)],num2cell(1:n_depth_groups),'uni',false));

%% PSTH for left vs. right stim, choose left vs. right stim (by depth)

% Group multiunit by depth
n_depth_groups = 6;
%depth_group_edges = linspace(0,max(channel_positions(:,2)),n_depth_groups+1);
%depth_group_edges = linspace(0,4000,n_depth_groups+1);
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
%depth_group_edges = [0 1300];
depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);

depth_group = discretize(spikeDepths,depth_group_edges);

raster_window = [-0.5,2.5];
psth_bin_size = 0.001;
smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

psth_right_hit = nan(n_depth_groups,diff(raster_window)/psth_bin_size);
psth_right_miss = nan(n_depth_groups,diff(raster_window)/psth_bin_size);
psth_left_hit = nan(n_depth_groups,diff(raster_window)/psth_bin_size);
psth_left_miss = nan(n_depth_groups,diff(raster_window)/psth_bin_size);

for curr_depth = 1:n_depth_groups
    
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
%     curr_spike_times = spike_times_timeline((depth_group == curr_depth) & ...
%         ismember(spike_templates,find(msn)));
    
    use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 1;
    if sum(use_trials) > 0
        [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,stimOn_times(use_trials),raster_window,psth_bin_size);
        psth_smooth = conv2(psth,smWin,'same');
        psth_right_hit(curr_depth,:) = psth_smooth;
    end
    
    use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 0;
    if sum(use_trials) > 0
        [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,stimOn_times(use_trials),raster_window,psth_bin_size);
        psth_smooth = conv2(psth,smWin,'same');
        psth_right_miss(curr_depth,:) = psth_smooth;
    end
    
    use_trials = signals_events.trialSideValues == -1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 1;
    if sum(use_trials) > 0
        [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,stimOn_times(use_trials),raster_window,psth_bin_size);
        psth_smooth = conv2(psth,smWin,'same');
        psth_left_hit(curr_depth,:) = psth_smooth;
    end
    
    use_trials = signals_events.trialSideValues == -1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 0;
    if sum(use_trials) > 0
        [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,stimOn_times(use_trials),raster_window,psth_bin_size);
        psth_smooth = conv2(psth,smWin,'same');
        psth_left_miss(curr_depth,:) = psth_smooth;
    end
    
end

figure; hold on;
trace_spacing = max([max(psth_right_hit(:)),max(psth_right_miss(:)),max(psth_left_hit(:)),max(psth_left_hit(:))]);
zs = true;
if zs
    trace_spacing = 5;
end
p_rh = AP_stackplot(psth_right_hit',bins,trace_spacing,zs,'k',depth_group_centers);
p_rm = AP_stackplot(psth_right_miss',bins,trace_spacing,zs,'r',depth_group_centers);
p_lh = AP_stackplot(psth_left_hit',bins,trace_spacing,zs,'b',depth_group_centers);
p_lm = AP_stackplot(psth_left_miss',bins,trace_spacing,zs,'m',depth_group_centers);

median_move_time = nanmedian(wheel_move_time' - stimOn_times);
median_reward_time = nanmedian(reward_t_timeline - stimOn_times(signals_events.hitValues == 1)');

line([0,0],ylim,'color','k','linestyle','--');
line([median_move_time,median_move_time],ylim,'color','r','linestyle','--');
line([median_reward_time,median_reward_time],ylim,'color','b','linestyle','--');

legend([p_rh(1),p_rm(1),p_lh(1),p_lm(1)],{'Stim right hit','Stim right miss','Stim left hit','Stim left miss'})
xlabel('Time from stim');
ylabel('Depth (\mum)');


%% PSTH for left vs. right stim, choose left vs. right stim (errorbars)

stimIDs = signals_events.trialSideValues.*signals_events.trialContrastValues;

% use_spikes_idx = ismember(spike_templates,find(templateDepths >= 1000 & templateDepths <= 2000));
use_spikes_idx =ismember(spike_templates,find(templateDepths > 1000 & templateDepths < 2000)) &...
    ismember(spike_templates,find(msn));
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
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(use_spikes,stimOn_times(use_trials),raster_window,psth_bin_size);
AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'k');

use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 0;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(use_spikes,stimOn_times(use_trials),raster_window,psth_bin_size);
AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'r');

use_trials = signals_events.trialSideValues == -1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 1;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(use_spikes,stimOn_times(use_trials),raster_window,psth_bin_size);
AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'b');

use_trials = signals_events.trialSideValues == -1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 0;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(use_spikes,stimOn_times(use_trials),raster_window,psth_bin_size);
AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'m');

legend({'Right hit','Right miss','Left hit','Left miss'})
xlim([bins(1),bins(end)])
xlabel('Time from stim');
ylabel('Spikes');

subplot(1,3,2); hold on;

use_trials = (signals_events.trialSideValues == 1 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == -1 & signals_events.hitValues == 0);
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(use_spikes,wheel_move_time(use_trials),raster_window,psth_bin_size);
p1 = AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'k');

use_trials = (signals_events.trialSideValues == -1 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == 1 & signals_events.hitValues == 0);
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(use_spikes,wheel_move_time(use_trials),raster_window,psth_bin_size);
p2 = AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'r');

legend([p1(1),p2(1)],{'Move left','Move right'});
xlim([bins(1),bins(end)])
xlabel('Time from movement');
ylabel('Spikes');


subplot(1,3,3); hold on;

use_trials = (signals_events.trialSideValues == 1 & signals_events.trialContrastValues == 0 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == -1 & signals_events.trialContrastValues == 0 & signals_events.hitValues == 0);
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(use_spikes,wheel_move_time(use_trials),raster_window,psth_bin_size);
AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'k');

use_trials = (signals_events.trialSideValues == -1 & signals_events.trialContrastValues == 0 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == 1 & signals_events.trialContrastValues == 0 & signals_events.hitValues == 0);
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(use_spikes,wheel_move_time(use_trials),raster_window,psth_bin_size);
AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'r');

legend({'Move left 0 contrast','Move right 0 contrast'});
xlim([bins(1),bins(end)])
xlabel('Time from movement');
ylabel('Spikes');


%% PSTH population choose left vs. right stim (move-aligned)

stimIDs = signals_events.trialSideValues.*signals_events.trialContrastValues;

use_spikes_idx = ismember(spike_templates,find(templateDepths >= depth_group_edges(3) & templateDepths <= depth_group_edges(4)));
%use_spikes_idx = ismember(spike_templates,find(templateDepths >= 3000 & templateDepths <= 4000));
% use_spikes_idx = ismember(spike_templates,intersect(find(templateDepths >= 1000 & templateDepths <= 2000),find(msn)));

use_spikes = spike_times_timeline(use_spikes_idx);

raster_window = [-1 2];
psth_bin_size = 0.02;

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

% Plot
figure; hold on;

use_trials = (signals_events.trialSideValues == 1 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == -1 & signals_events.hitValues == 0) & ...
    ~signals_events.repeatTrialValues;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(use_spikes,wheel_move_time(use_trials),raster_window,psth_bin_size);
p1 = AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'k');

use_trials = (signals_events.trialSideValues == -1 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == 1 & signals_events.hitValues == 0) & ...
    ~signals_events.repeatTrialValues;
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(use_spikes,wheel_move_time(use_trials),raster_window,psth_bin_size);
p2 = AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'r');

legend([p1(1),p2(1)],{'Move left','Move right'});
xlim([bins(1),bins(end)])
xlabel('Time from movement');
ylabel('Spikes');

line([0,0],ylim,'linestyle','--','color','k');


%% PSTH templates choose left vs. right stim (move-aligned)

stimIDs = signals_events.trialSideValues.*signals_events.trialContrastValues;

use_spikes_idx = ismember(spike_templates,find(templateDepths >= 500 & templateDepths <= 3500));
% use_spikes_idx = ismember(spike_templates,intersect(find(templateDepths >= 2000 & templateDepths <= 3000),find(msn)));
use_spikes = spike_times_timeline(use_spikes_idx);

use_templates = unique(spike_templates(use_spikes_idx));

% Plot
go_left_trials = (signals_events.trialSideValues == 1 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == -1 & signals_events.hitValues == 0) & ...
    ~signals_events.repeatTrialValues;

go_right_trials = (signals_events.trialSideValues == -1 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == 1 & signals_events.hitValues == 0) & ...
    ~signals_events.repeatTrialValues;

trial_choice = go_left_trials + 2.*go_right_trials;

raster_window = [-2,2];
psth_bin_size = 0.001;

template_psth_left = nan(length(use_templates),diff(raster_window/psth_bin_size));
for curr_template_idx = 1:length(use_templates)
    curr_template = use_templates(curr_template_idx);
    curr_use_spikes = spike_times_timeline(spike_templates == curr_template);
    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = ...
        psthAndBA(curr_use_spikes,wheel_move_time(go_left_trials),raster_window,psth_bin_size);
    template_psth_left(curr_template_idx,:) = psth;
end

template_psth_right = nan(length(use_templates),diff(raster_window/psth_bin_size));
for curr_template_idx = 1:length(use_templates)
    curr_template = use_templates(curr_template_idx);
    curr_use_spikes = spike_times_timeline(spike_templates == curr_template);
    [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = ...
        psthAndBA(curr_use_spikes,wheel_move_time(go_right_trials),raster_window,psth_bin_size);
    template_psth_right(curr_template_idx,:) = psth;
end

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

template_psth_left_smooth = conv2(template_psth_left,smWin,'same');
template_psth_right_smooth = conv2(template_psth_right,smWin,'same');

figure; hold on;
plot(nanmean(zscore(template_psth_left_smooth,[],2)),'k');
plot(nanmean(zscore(template_psth_right_smooth,[],2)),'r');

% Get left/right differences for correct non-repeat trials, plot raster and
% sort by difference
sort_time = raster_window(1):psth_bin_size:raster_window(2) > 0 & ...
    raster_window(1):psth_bin_size:raster_window(2) < 0.5;

l_r_diff = (sum(template_psth_left_smooth(:,sort_time),2) - sum(template_psth_right_smooth(:,sort_time),2))./ ...
    (sum(template_psth_left_smooth(:,sort_time),2) + sum(template_psth_right_smooth(:,sort_time),2));

[~,sort_idx] = sort(l_r_diff);
sort_templates = nan(size(templates,1),1);
sort_templates(use_templates(sort_idx)) = 1:length(use_templates);
template_sort = sort_templates(spike_templates(use_spikes_idx));

raster_window = [-1,1];
psthViewer(use_spikes,template_sort, ...
    [wheel_move_time(go_right_trials),wheel_move_time(go_left_trials)], ...
    raster_window,[ones(1,sum(go_right_trials)),2*ones(1,sum(go_left_trials))]);

sorted_diff = template_psth_left_smooth(sort_idx,:) - template_psth_right_smooth(sort_idx,:);
figure;imagesc(sorted_diff);
colormap(colormap_BlueWhiteRed);
caxis([-max(abs(sorted_diff(:))),max(abs(sorted_diff(:)))]);

%% Get average fluorescence to Signals event

% Define the window to get an aligned response to
surround_window = [-0.5,2];

% Define the times to align to
use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues >= 0.5 & signals_events.hitValues == 1;
align_times = stimOn_times(use_trials(1:length(stimOn_times)));

% Get the surround time
framerate = 1./nanmedian(diff(frame_t));
surround_samplerate = 1/(framerate*1);
t_surround = surround_window(1):surround_samplerate:surround_window(2);

% Don't use times that fall outside of imaging
align_times(align_times + t_surround(1) < frame_t(2) | ...
    align_times + t_surround(2) > frame_t(end)) = [];

% Use closest frames to times
align_surround_times = bsxfun(@plus, align_times, t_surround);
frame_edges = [frame_t,frame_t(end)+1/framerate];
align_frames = discretize(align_surround_times,frame_edges);

% If any aligned V's are NaNs (when does this happen?), don't use
align_frames(any(isnan(align_frames),2),:) = [];

aligned_V = reshape(fVdf(:,align_frames'), ...
    size(fV,1),size(align_frames,2),size(align_frames,1));

mean_aligned_V = nanmean(aligned_V,3);

% Get and plot the average fluorescence around event
mean_aligned_px = svdFrameReconstruct(Udf,mean_aligned_V);



a = diff(imgaussfilt(mean_aligned_px,2),[],3);
a(a < 0) = 0;
AP_image_scroll(a,t_surround);
axis image;

%% Align fluorescence and MUA to task event across trials

depth_edges = [depth_group_edges(4),depth_group_edges(5)];
sample_rate_factor = 3;

% Define times to align
% SIGNALS - CHOICEWORLD

n_trials = length(block.paramsValues);
trial_conditions = signals_events.trialSideValues(1:n_trials).*signals_events.trialContrastValues(1:n_trials);
trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
stim_to_move = padarray(wheel_move_time - stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
stim_to_feedback = padarray(signals_events.responseTimes, ...
    [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
    padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');

use_trials = ...
    trial_outcome ~= 0 & ...
    ~signals_events.repeatTrialValues(1:n_trials) & ...
    stim_to_move < 0.5 & ...
    stim_to_feedback < 1.5;

use_trials_idx = find(use_trials);
align_times = reshape(wheel_move_time(use_trials),[],1);
hit_trials = signals_events.hitValues(use_trials) == 1;
stim_conditions = signals_events.trialContrastValues(use_trials);

% MPEP PASSIVE
% align_times = reshape(stimOn_times,[],1);
% hit_trials = rand(size(align_times)) > 0.5;
% stim_conditions = stimIDs;

interval_surround = [-0.5,1.5];
sample_rate = framerate*sample_rate_factor;
t_surround = interval_surround(1):1/sample_rate:interval_surround(2);
t_peri_event = bsxfun(@plus,align_times,t_surround);

% Draw ROI and align fluorescence
Udf_aligned = AP_align_widefield(animal,day,Udf);
[roi_trace,roi_mask] = AP_svd_roi(Udf_aligned,fVdf,'master'); % weight_im, retinotopic_map, response_im
event_aligned_f = interp1(frame_t,roi_trace,t_peri_event);
event_aligned_df = interp1(conv(frame_t,[1,1]/2,'valid'),diff(roi_trace),t_peri_event);
event_aligned_df(event_aligned_df < 0) = 0;

% Pull out MUA at a given depth
use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > depth_edges(1) & templateDepths < depth_edges(2))));
% use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 500 & templateDepths < 1500)) &...
%     ismember(spike_templates,find(msn)));
t_peri_event_bins = [t_peri_event - 1/(sample_rate*2), ...
    t_peri_event(:,end) + 1/(sample_rate*2)];
event_aligned_spikes = cell2mat(arrayfun(@(x) ...
    histcounts(use_spikes,t_peri_event_bins(x,:)),[1:length(align_times)]','uni',false));

figure; colormap(gray);
subplot(2,2,1);
imagesc(t_surround,1:length(align_times),event_aligned_df)
line([0,0],ylim,'color','r');
ylabel('Event number')
xlabel('Time (s)')
title('\DeltaF')

subplot(2,2,2);
imagesc(t_surround,1:length(align_times),event_aligned_spikes)
line([0,0],ylim,'color','r');
ylabel('Event number')
xlabel('Time (s)')
title('Spikes')

subplot(2,2,3:4); hold on;
plot(t_surround,mat2gray(nanmean(event_aligned_df,1)),'k','linewidth',2)
plot(t_surround,mat2gray(nanmean(event_aligned_spikes,1)),'b','linewidth',2)
line([0,0],ylim,'color','r');
ylabel('Normalized units')
xlabel('Time from event')

% Plot the average response
t_avg = t_surround > 0.05 & t_surround < 0.15;
plot(t_surround,t_avg,'m')

figure;

% (trial-trial)
ax1 = subplot(1,2,1); hold on;

df_peak_hit = nanmean(event_aligned_df(hit_trials,t_avg),2);
spikes_peak_hit = nanmean(event_aligned_spikes(hit_trials,t_avg),2);

df_peak_miss = nanmean(event_aligned_df(~hit_trials,t_avg),2);
spikes_peak_miss = nanmean(event_aligned_spikes(~hit_trials,t_avg),2);

fit_hit = robustfit(df_peak_hit,spikes_peak_hit);
fit_miss = robustfit(df_peak_miss,spikes_peak_miss);

[used_conditions,~,revalued_conditions] = unique(stim_conditions);
col = copper(length(used_conditions));
scatter(df_peak_hit,spikes_peak_hit,50,col(revalued_conditions(hit_trials),:),'filled','MarkerEdgeColor','b');
scatter(df_peak_miss,spikes_peak_miss,50,col(revalued_conditions(~hit_trials),:),'filled','MarkerEdgeColor','r');
axis square tight;

line(xlim,xlim*fit_hit(2)+fit_hit(1)','color','b')
line(xlim,xlim*fit_miss(2)+fit_miss(1),'color','r')
xlabel('\DeltaF');
ylabel('Spikes');
title('Trial');

% (contrast mean)
ax2 = subplot(1,2,2); hold on;

df_peak_hit_contrastmean = grpstats(df_peak_hit,revalued_conditions(hit_trials),'mean');
spikes_peak_hit_contrastmean = grpstats(spikes_peak_hit,revalued_conditions(hit_trials),'mean');

df_peak_miss_contrastmean = grpstats(df_peak_miss,revalued_conditions(~hit_trials),'mean');
spikes_peak_miss_contrastmean = grpstats(spikes_peak_miss,revalued_conditions(~hit_trials),'mean');

fit_hit_contrastmean = robustfit(df_peak_hit_contrastmean,spikes_peak_hit_contrastmean);
fit_miss_contrastmean = robustfit(df_peak_miss_contrastmean,spikes_peak_miss_contrastmean);

col = copper(length(used_conditions));
hit_contrasts = unique(revalued_conditions(hit_trials));
miss_contrasts = unique(revalued_conditions(~hit_trials));
scatter(df_peak_hit_contrastmean,spikes_peak_hit_contrastmean,100,col(hit_contrasts,:),'filled','MarkerEdgeColor','b');
scatter(df_peak_miss_contrastmean,spikes_peak_miss_contrastmean,100,col(miss_contrasts,:),'filled','MarkerEdgeColor','r');
axis square tight;

line(xlim,xlim*fit_hit_contrastmean(2)+fit_hit_contrastmean(1)','color','b')
line(xlim,xlim*fit_miss_contrastmean(2)+fit_miss_contrastmean(1),'color','r')
xlabel('\DeltaF');
ylabel('Spikes');
title('Contrast mean');

% (total session, put on both plots)
skip_seconds = 60;
use_frames = frame_t > skip_seconds & abs(frame_t-frame_t(end)) > skip_seconds;
roi_trace_df = diff(roi_trace);
frame_t_df = conv(frame_t,[1,1]/2,'valid');
frame_t_df_resample = frame_t_df(1):1/sample_rate:frame_t_df(end);
roi_trace_df_resample = interp1(frame_t_df,roi_trace_df,frame_t_df_resample);
roi_trace_df_resample(roi_trace_df < 0) = 0;
spike_bins = [frame_t_df_resample-1/sample_rate,frame_t_df_resample(end)+1/sample_rate];
binned_spikes = histcounts(use_spikes,spike_bins);

fit_total = robustfit(roi_trace_df_resample(use_frames)',binned_spikes(use_frames)');
line(ax1,xlim(ax1),xlim(ax1)*fit_total(2)+fit_total(1)','color','k','linestyle','--')
line(ax2,xlim(ax2),xlim(ax2)*fit_total(2)+fit_total(1)','color','k','linestyle','--')

legend(ax2,{'Correct','Incorrect','Correct fit','Incorrect fit','Session fit'});

% Get significant correlations from df to spikes
corr_real = (zscore(event_aligned_df,1)'*zscore(event_aligned_spikes,1))./(size(event_aligned_df,1)-1);
n_shuff = 1000;
corr_shuff = nan(size(corr_real,1),size(corr_real,2),n_shuff);
warning off
for curr_shuff = 1:n_shuff
    corr_shuff(:,:,curr_shuff) = (...
        zscore(shake(event_aligned_df,1),1)'* ...
        zscore(event_aligned_spikes,1))./(size(event_aligned_df,1)-1);
end
warning on
corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) - (corr_real < corr_shuff_cutoff(:,:,1));

figure;
subplot(4,4,[1,2,3,5,6,7,9,10,11])
imagesc(t_surround,t_surround,corr_sig);
colormap(gray);
ylabel('\DeltaF');
xlabel('Spikes');

if exist('reward_t_timeline','var')
    median_stim_time_rel = nanmedian(stimOn_times(use_trials) - align_times);
    median_move_time_rel = nanmedian(wheel_move_time(use_trials)' - align_times);
    reward_trial_idx = find(signals_events.hitValues == 1);
    median_reward_time_rel = nanmedian( ...
        reward_t_timeline(ismember(reward_trial_idx,use_trials_idx(hit_trials))) - align_times(hit_trials)');
end

line(xlim,ylim,'color','r');
line([0,0],ylim,'color','m');
line(xlim,[0,0],'color','m');
if exist('reward_t_timeline','var')
    line(repmat(median_stim_time_rel,1,2),ylim,'color','g','linestyle','--');
    line(xlim,repmat(median_stim_time_rel,1,2),'color','g','linestyle','--');
    line(repmat(median_move_time_rel,1,2),ylim,'color','r','linestyle','--');
    line(xlim,repmat(median_move_time_rel,1,2),'color','r','linestyle','--');
    line(repmat(median_reward_time_rel,1,2),ylim,'color','b','linestyle','--');
    line(xlim,repmat(median_reward_time_rel,1,2),'color','b','linestyle','--');
end

subplot(4,4,[4,8,12]);
plot(nanmean(event_aligned_df,1),t_surround,'k','linewidth',2);
line(xlim,[0,0],'color','m');
if exist('reward_t_timeline','var')
    line(xlim,repmat(median_stim_time_rel,1,2),'color','g','linestyle','--');
    line(xlim,repmat(median_move_time_rel,1,2),'color','r','linestyle','--');
    line(xlim,repmat(median_reward_time_rel,1,2),'color','b','linestyle','--');
end
axis tight off; set(gca,'YDir','reverse');

subplot(4,4,[13,14,15]);
plot(t_surround,nanmean(event_aligned_spikes,1),'k','linewidth',2);
line([0,0],ylim,'color','m');
if exist('reward_t_timeline','var')
line(repmat(median_stim_time_rel,1,2),ylim,'color','g','linestyle','--');
    line(repmat(median_move_time_rel,1,2),ylim,'color','r','linestyle','--');
    line(repmat(median_reward_time_rel,1,2),ylim,'color','b','linestyle','--');
end
axis tight off;

% Get time-varying fit from fluorescence to spikes 
smooth_size = 10;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
event_aligned_df_smooth = conv2(event_aligned_df,smWin,'same');
event_aligned_spikes_smooth = conv2(event_aligned_spikes,smWin,'same');

f = figure('Position',[464,558,776,342]);
ax1 = axes; hold on; axis square
subplot(1,2,1,ax1);
ylabel('Normalized units');
xlabel('Time from event');
plot(t_surround,mat2gray(nanmean(event_aligned_df_smooth,1)),'k','linewidth',2);
plot(t_surround,mat2gray(nanmean(event_aligned_spikes_smooth,1)),'b','linewidth',2);
t_line = line([0,0],ylim,'color','k');

ax2 = axes; hold on; axis square
subplot(1,2,2,ax2);
xlabel('\DeltaF');
ylabel('Spikes');
xlim([0,max(event_aligned_df_smooth(:))]);
ylim([0,max(event_aligned_spikes_smooth(:))])

hit_plot = plot(event_aligned_df_smooth(hit_trials,1),event_aligned_spikes_smooth(hit_trials,1),'.b');
miss_plot = plot(event_aligned_df_smooth(~hit_trials,1),event_aligned_spikes_smooth(~hit_trials,1),'.r');

hit_fit_line = line(xlim,ylim,'color','b');
miss_fit_line = line(xlim,ylim,'color','r');

hit_fit_t = nan(size(event_aligned_df_smooth,2),2);
miss_fit_t = nan(size(event_aligned_df_smooth,2),2);
plot_fit = diag(corr_sig);
warning off;
for i = 1:size(event_aligned_df_smooth,2)
    set(hit_plot,'XData',event_aligned_df_smooth(hit_trials,i),'YData',event_aligned_spikes_smooth(hit_trials,i));
    set(miss_plot,'XData',event_aligned_df_smooth(~hit_trials,i),'YData',event_aligned_spikes_smooth(~hit_trials,i));
    
    curr_fit_hit = robustfit(event_aligned_df_smooth(hit_trials,i),event_aligned_spikes_smooth(hit_trials,i));
    curr_fit_miss = robustfit(event_aligned_df_smooth(~hit_trials,i),event_aligned_spikes_smooth(~hit_trials,i));
    
    if plot_fit(i)
        set(hit_fit_line,'XData',xlim,'YData',xlim*curr_fit_hit(2)+curr_fit_hit(1));
        set(miss_fit_line,'XData',xlim,'YData',xlim*curr_fit_miss(2)+curr_fit_miss(1));
    else
        set(hit_fit_line,'XData',[],'YData',[]);
        set(miss_fit_line,'XData',[],'YData',[]);
    end
    
    set(t_line,'XData',repmat(t_surround(i),1,2));
    
    hit_fit_t(i,:) = curr_fit_hit;
    miss_fit_t(i,:) = curr_fit_miss;
   
    frames(i) = getframe(f);
end
warning on;
close(f);

% % (write this movie to a file)
% fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\presentations\171114_sfn\figs\aligned_fluor_spikes_regression.mov';
% writerObj = VideoWriter(fn);
% writerObj.FrameRate = samplerate/4;
% open(writerObj);
% writeVideo(writerObj,frames);
% close(writerObj);

hit_fit_t(~plot_fit,:) = NaN;
miss_fit_t(~plot_fit,:) = NaN;

figure;
subplot(2,1,1); hold on;
plot(t_surround,hit_fit_t(:,1),'b','linewidth',2);
plot(t_surround,miss_fit_t(:,1),'r','linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Fit intercept');
xlabel('Time from event');
legend({'Hit','Miss'})
subplot(2,1,2); hold on;
plot(t_surround,hit_fit_t(:,2),'b','linewidth',2);
plot(t_surround,miss_fit_t(:,2),'r','linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Fit slope');
xlabel('Time from event');

%% Compare predicted to actual spikes by condition

depth_edges = [str_depth(1),str_depth(1)+1000];
sample_rate_factor = 1;

% SIGNALS - CHOICEWORLD
use_trials = signals_events.hitValues == 1 & signals_events.repeatTrialValues == 0;
align_times = reshape(stimOn_times(use_trials),[],1);
% align_times = reshape(signals_events.responseTimes(use_trials),[],1);
% align_times = reshape(reward_t_timeline,[],1);
% align_times = reshape(wheel_move_time(hit_trials),[],1);
stim_conditions = signals_events.trialContrastValues(use_trials).*signals_events.trialSideValues(use_trials);

% MPEP PASSIVE
% align_times = reshape(stimOn_times,[],1);
% stim_conditions = stimIDs;

interval_surround = [-0.5,1.5];
sample_rate = (1/median(diff(frame_t)))*sample_rate_factor;

t_surround = interval_surround(1):1/sample_rate:interval_surround(2);
t_peri_event = bsxfun(@plus,align_times,t_surround);

% Skip the first/last n seconds for prediction
skip_seconds = 60*1;

time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > depth_edges(1) & templateDepths < depth_edges(2))));
% use_spikes = spike_times_timeline(ismember(spike_templates, ...
%     find(templateDepths > depth_edges(1) & templateDepths < depth_edges(2))) & ...
%     ismember(spike_templates,find(msn)));

binned_spikes = histcounts(use_spikes,time_bins);

use_svs = 1:50;
kernel_frames = -35:17;
lambda = 2e5;
zs = [false,true];
cvfold = 5;

fVdf_resample = interp1(frame_t,fVdf(use_svs,:)',time_bin_centers)';
dfVdf_resample = interp1(conv(frame_t,[1,1]/2,'valid'),diff(fVdf(use_svs,:),[],2)',time_bin_centers)';

% [k,predicted_spikes,explained_var] = ...
%     AP_regresskernel(fVdf_resample, ...
%     binned_spikes,kernel_frames,lambda,zs,cvfold);
[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(dfVdf_resample, ...
    binned_spikes,kernel_frames,lambda,zs,cvfold);

binned_spikes_std = std(binned_spikes);
binned_spikes_mean = mean(binned_spikes);
predicted_spikes_reranged = predicted_spikes*binned_spikes_std+binned_spikes_mean;



% find predicted -> real nonlinearity
[grp_binned_spikes,grp_predicted_spikes]= grpstats(predicted_spikes_reranged,binned_spikes,{'gname','median'});
grp_binned_spikes = cellfun(@str2num,grp_binned_spikes);

pred_offset = grp_predicted_spikes(grp_binned_spikes == 0);
pred_exp = log(grp_predicted_spikes(grp_binned_spikes ~= 0))\ ...
    log(grp_binned_spikes(grp_binned_spikes ~= 0));

% figure;
% plot((grp_predicted_spikes-pred_offset).^pred_exp,grp_binned_spikes,'.k')
% xlim([0,100]);ylim([0,100]);
% line(xlim,ylim,'color','r');

predicted_spikes_rectified = (predicted_spikes_reranged-pred_offset);
predicted_spikes_rectified(predicted_spikes_rectified < 0) = 0;
predicted_spikes_nlin = predicted_spikes_rectified.^pred_exp;

% nonlinear-fixed explained variance
sse_signals = sum(binned_spikes.^2,2);
sse_total_residual = sum(bsxfun(@minus,predicted_spikes_nlin,binned_spikes).^2,2);
explained_var_nlin = (sse_signals - sse_total_residual)./sse_signals;



spikes_real_aligned = interp1(time_bin_centers,binned_spikes,t_peri_event);
spikes_pred_aligned = interp1(time_bin_centers,predicted_spikes_nlin,t_peri_event);

spikes_real_aligned_mean = grpstats(spikes_real_aligned,stim_conditions);
spikes_pred_aligned_mean = grpstats(spikes_pred_aligned,stim_conditions);

% Plot all responses
figure; hold on
p1 = AP_stackplot(spikes_real_aligned_mean',t_surround,20,false,'k',unique(stim_conditions));
p2 = AP_stackplot(spikes_pred_aligned_mean',t_surround,20,false,'r');
ylabel('Stim');
xlabel('Time from event onset');
legend([p1(1),p2(1)],{'Real','Predicted'});
line([0,0],ylim,'color','k');

% Plot responses by condition and error
response_t = [0,0.3];
response_t_use = t_surround >= response_t(1) & t_surround <= response_t(2);
spikes_real_response = nanmean(spikes_real_aligned(:,response_t_use),2);
spikes_pred_response = nanmean(spikes_pred_aligned(:,response_t_use),2);

figure;
pred_error = spikes_real_response - spikes_pred_response;
[pred_error_mean,pred_error_sem] = grpstats(pred_error,stim_conditions,{'mean','sem'});
errorbar(unique(stim_conditions),pred_error_mean,pred_error_sem,'k','linewidth',2)
line(xlim,[0,0],'color','r')
ylabel('Predicted spikes error');
xlabel('Stim condition');

%% Align fluorescence, MUA, and predicted MUA, compare across conditions

depth_edges = [str_depth(1),str_depth(1)+200];
sample_rate_factor = 1;

% Define times to align
% CHOICEWORLD
% use_trials = signals_events.trialSideValues == 1 & ~isnan(wheel_move_time);
use_trials = ~isnan(wheel_move_time);
% align_times = reshape(stimOn_times(use_trials(1:length(stimOn_times))),[],1);
% stim_conditions = signals_events.trialContrastValues(use_trials);
align_times = reshape(stimOn_times(use_trials(1:length(stimOn_times))),[],1);
stim_conditions = signals_events.trialContrastValues(use_trials).*signals_events.trialSideValues(use_trials);

% % MPEP PASSIVE
% use_trials = true(size(stimIDs));
% align_times = reshape(stimOn_times,[],1);
% stim_conditions = stimIDs;

[used_conditions,~,revalued_conditions] = unique(stim_conditions);

interval_surround = [-0.5,2.5];
sample_rate = framerate*sample_rate_factor;
t_surround = interval_surround(1):1/sample_rate:interval_surround(2);
t_peri_event = bsxfun(@plus,align_times,t_surround);

% Pull out MUA at a given depth

use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > depth_edges(1) & templateDepths < depth_edges(2))));
% use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > depth_edges(1) & templateDepths < depth_edges(2))) &...
%     ismember(spike_templates,find(msn)));

skip_seconds = 60*1;

time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > depth_edges(1) & templateDepths < depth_edges(2))));
binned_spikes = histcounts(use_spikes,time_bins);

event_aligned_spikes = interp1(time_bin_centers,binned_spikes,t_peri_event);

% Define spikes time of interest
[~,sort_idx] = sort(revalued_conditions);
h = figure;
imagesc(t_surround,1:length(align_times),event_aligned_spikes(sort_idx,:));
colormap(flipud(gray));
title('Mark range start')
[t_start,~] = ginput(1);
title('Mark range end');
[t_stop,~] = ginput(1);
close(h);drawnow;

t_avg = t_surround >= t_start & t_surround <= t_stop;

% Predict and align spikes from fluorescence
use_svs = 1:50;
kernel_frames = -35:17;
lambda = 2e5;
zs = [false,true];
cvfold = 5;

dfVdf_resample = interp1(conv(frame_t,[1,1]/2,'valid'),diff(fVdf(use_svs,:),[],2)',time_bin_centers)';
[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(dfVdf_resample, ...
    binned_spikes,kernel_frames,lambda,zs,cvfold);

binned_spikes_std = std(binned_spikes);
binned_spikes_mean = mean(binned_spikes);
predicted_spikes_reranged = predicted_spikes*binned_spikes_std+binned_spikes_mean;

event_aligned_spikes_predicted = interp1(time_bin_centers,predicted_spikes_reranged,t_peri_event);

% Divide and plot everything by contrast/side/movement direction/hit-miss
spikes_peak = nanmean(event_aligned_spikes(:,t_avg),2);
spikes_predicted_peak = nanmean(event_aligned_spikes_predicted(:,t_avg),2);

[spikes_grp_mean,spikes_grp_sem,spikes_grp] = ...
    grpstats(spikes_peak,stim_conditions,{'mean','sem','gname'});
spikes_grp = cellfun(@str2num,spikes_grp);

[spikes_predicted_grp_mean,spikes_predicted_grp_sem,spikes_predicted_grp] = ...
    grpstats(spikes_predicted_peak,stim_conditions,{'mean','sem','gname'});
spikes_predicted_grp = cellfun(@str2num,spikes_predicted_grp);

figure;

subplot(2,2,1); hold on;
plot(t_surround,nanmean(event_aligned_spikes_predicted,1),'r','linewidth',2)
plot(t_surround,nanmean(event_aligned_spikes,1),'k','linewidth',2)
patch([t_start,t_stop,t_stop,t_start],[repmat(min(ylim),1,2),repmat(max(ylim),1,2)],'y','EdgeColor','none');
set(gca,'Children',flipud(get(gca,'Children')))
ylabel('Spikes');
xlabel('Time from event');
legend({'Time used','Real','Predicted'})
axis tight;
line([0,0],ylim,'color','k');

subplot(2,2,2); hold on;
if exist('signals_events','var')
    hit_trials = signals_events.hitValues(use_trials) == 1;
    miss_trials = signals_events.missValues(use_trials) == 1;
    ipsi_stim_trials = signals_events.trialSideValues(use_trials) == -1;
    contra_stim_trials = signals_events.trialSideValues(use_trials) == 1;
    
    col = copper(length(used_conditions));
    p1 = scatter(spikes_peak(hit_trials & contra_stim_trials),spikes_predicted_peak(hit_trials & contra_stim_trials),70, ...
        col(revalued_conditions(hit_trials & contra_stim_trials),:),'filled','MarkerEdgeColor','b','LineWidth',2);
    p2 = scatter(spikes_peak(miss_trials & contra_stim_trials),spikes_predicted_peak(miss_trials & contra_stim_trials),70, ...
        col(revalued_conditions(miss_trials & contra_stim_trials),:),'filled','MarkerEdgeColor','r','LineWidth',2);
    legend([p1(1),p2(1)],{'Correct','Incorrect'});
else
    col = copper(length(used_conditions));
    p1 = scatter(spikes_peak(contra_stim_trials),spikes_predicted_peak(contra_stim_trials),70, ...
        col(revalued_conditions(contra_stim_trials),:),'filled','MarkerEdgeColor','k');
end
line([0,max([spikes_peak;spikes_predicted_peak])],[0,max([spikes_peak;spikes_predicted_peak])],'color','k');
axis square tight;

xlabel('Spikes');
ylabel('Predicted spikes');
title('Trials')

% Plot average responses and predictions by contrast
subplot(2,2,3); hold on;
errorbar(spikes_grp,spikes_grp_mean,spikes_grp_sem,'k','linewidth',2);
errorbar(spikes_predicted_grp,spikes_predicted_grp_mean,spikes_predicted_grp_sem,'r','linewidth',2);
axis tight;
ylabel('Spikes');
xlabel('Contrast');

% Fit by condition and correctness
subplot(2,2,4);

spike_lim = [0,max([spikes_peak;spikes_predicted_peak])];
line(spike_lim,spike_lim);

fit_conditions = nan(length(used_conditions),2);
for curr_condition = 1:length(used_conditions)
    fit_conditions(curr_condition,:) = ...
        robustfit(spikes_peak(revalued_conditions == curr_condition),spikes_predicted_peak(revalued_conditions == curr_condition));   
    line(spike_lim,spike_lim*fit_conditions(curr_condition,2)+ ...
        fit_conditions(curr_condition,1)','color',col(curr_condition,:),'linewidth',2)
end
if exist('signals_events','var')
    fit_hit = robustfit(spikes_peak(hit_trials),spikes_predicted_peak(hit_trials));
    fit_miss = robustfit(spikes_peak(miss_trials),spikes_predicted_peak(miss_trials));
    
    line(spike_lim,spike_lim*fit_hit(2)+ ...
        fit_hit(1)','color','b','linewidth',2)
    line(spike_lim,spike_lim*fit_miss(2)+ ...
        fit_miss(1)','color','r','linewidth',2)
end

axis square tight;
xlabel('Spikes');
ylabel('Predicted spikes');
title('Fits')


%% THIS WAS FROM ABOVE: TO PLOT FLUOR, BUT NOT DIRECTLY RELATED, SO MOVED

%%%%

% Draw ROI and align fluorescence
[roi_trace,roi_mask] = AP_svd_roi(Udf,fVdf,weight_im,retinotopic_map); % weight_im, retinotopic_map, response_im
event_aligned_f = interp1(frame_t,roi_trace,t_peri_event);
event_aligned_df = interp1(conv(frame_t,[1,1]/2,'valid'),diff(roi_trace),t_peri_event);
event_aligned_df(event_aligned_df < 0) = 0;

%%%%

% % To use ROI fluorescence to predict spikes
% roi_trace_df_resample = interp1(conv(frame_t,[1,1]/2,'valid'),diff(roi_trace)',time_bin_centers);
% roi_trace_df_resample(roi_trace_df_resample < 0) = 0;
% [k,predicted_spikes,explained_var] = ...
%     AP_regresskernel(roi_trace_df_resample, ...
%     binned_spikes,kernel_frames,lambda,zs,cvfold);

%%%%

df_peak = nanmean(event_aligned_df(:,t_avg),2);

[fluor_grp_mean,fluor_grp_sem,fluor_grp] = ...
    grpstats(df_peak,stim_conditions,{'mean','sem','gname'});
fluor_grp = cellfun(@str2num,fluor_grp);

%%%%

figure;
subplot(2,2,1); hold on; axis tight;
errorbar(fluor_grp,fluor_grp_mean,fluor_grp_sem,'k','linewidth',2)
xlabel('Contrast');
ylabel('Fluorescence');

subplot(2,2,2); hold on; axis tight;
errorbar(spikes_grp,spikes_grp_mean,spikes_grp_sem,'k','linewidth',2)
errorbar(spikes_predicted_grp,spikes_predicted_grp_mean,spikes_predicted_grp_sem,'r','linewidth',2)
xlabel('Contrast');
ylabel('Spikes');
legend({'Real','Predicted'});


%% !!!!!!!! BATCH PROCESSED ANALYSIS !!!!!!!!!

%% Load batch widefield passive

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% protocol = 'stimKalatsky';
protocol = 'AP_choiceWorldStimPassive';

surround_window = [-0.5,5];
framerate = 35;
surround_samplerate = 1/(framerate*1);
t_surround = surround_window(1):surround_samplerate:surround_window(2);

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\passive';
im = cell(size(animals));
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    fn = [data_path filesep animal '_' protocol '.mat'];
    if ~exist(fn,'file');
        continue
    end
    load(fn);
    im{curr_animal} = im_aligned_average;
    AP_print_progress_fraction(curr_animal,length(animals));
end

% Get ddf
ddf_im = im;
for curr_animal = 1:length(im)
    if isempty(ddf_im{curr_animal})
        continue
    end
    curr_im = ddf_im{curr_animal};
    curr_im(isnan(curr_im)) = 0;
    curr_im = diff(curr_im,[],3);
    curr_im(curr_im < 0) = 0;
    ddf_im{curr_animal} = curr_im;
end

im = nanmean(cat(5,im{:}),5);
ddf_im = nanmean(cat(5,ddf_im{:}),5);




%% Make batch widefield choiceworld mean (stim aligned)

clear all
disp('Early move hit');
animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
im_stim_earlymove_hit_avg_combined = zeros(437,416,265,11);
n_conditions = [];
for curr_animal = 1:length(animals) 
    animal = animals{curr_animal};
    fn = [data_path filesep animal '_im_stim_earlymove_hit_avg.mat'];
    load(fn);
    
    im_stim_earlymove_hit_avg_combined = nansum(cat(5,im_stim_earlymove_hit_avg_combined,im_stim_earlymove_hit_avg),5); 
    n_conditions = sum(cat(5,n_conditions,any(any(any(im_stim_earlymove_hit_avg,1),2),3)),5);
    AP_print_progress_fraction(curr_animal,length(animals));
end
im_stim_earlymove_hit_avg_combined = bsxfun(@rdivide,im_stim_earlymove_hit_avg_combined,n_conditions);
save_fn = [data_path filesep 'im_stim_earlymove_hit_avg_combined.mat'];
save(save_fn,'im_stim_earlymove_hit_avg_combined','-v7.3');
disp('Done')

clear all
disp('Late move hit');
animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
im_stim_latemove_hit_avg_combined = zeros(437,416,194,11);
n_conditions = [];
for curr_animal = 1:length(animals) 
    animal = animals{curr_animal};
    fn = [data_path filesep animal '_im_stim_latemove_hit_avg.mat'];
    load(fn);
    
    im_stim_latemove_hit_avg_combined = nansum(cat(5,im_stim_latemove_hit_avg_combined,im_stim_latemove_hit_avg),5); 
    n_conditions = sum(cat(5,n_conditions,any(any(any(im_stim_latemove_hit_avg,1),2),3)),5);
    AP_print_progress_fraction(curr_animal,length(animals));
end
im_stim_latemove_hit_avg_combined = bsxfun(@rdivide,im_stim_latemove_hit_avg_combined,n_conditions);
save_fn = [data_path filesep 'im_stim_latemove_hit_avg_combined.mat'];
save(save_fn,'im_stim_latemove_hit_avg_combined','-v7.3');
disp('Done')

clear all
disp('Early move miss');
animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
im_stim_earlymove_miss_avg_combined = zeros(437,416,194,11);
n_conditions = [];
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    fn = [data_path filesep animal '_im_stim_earlymove_miss_avg.mat'];
    load(fn);
    
    im_stim_earlymove_miss_avg_combined = nansum(cat(5,im_stim_earlymove_miss_avg_combined,im_stim_earlymove_miss_avg),5); 
    n_conditions = sum(cat(5,n_conditions,any(any(any(im_stim_earlymove_miss_avg,1),2),3)),5);
    AP_print_progress_fraction(curr_animal,length(animals));
end
im_stim_earlymove_miss_avg_combined = bsxfun(@rdivide,im_stim_earlymove_miss_avg_combined,n_conditions);
save_fn = [data_path filesep 'im_stim_earlymove_miss_avg_combined.mat'];
save(save_fn,'im_stim_earlymove_miss_avg_combined','-v7.3');
disp('Done')

clear all
disp('Late move miss');
animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
im_stim_latemove_miss_avg_combined = zeros(437,416,194,11);
n_conditions = [];
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    fn = [data_path filesep animal '_im_stim_latemove_miss_avg.mat'];
    load(fn);
    
    im_stim_latemove_miss_avg_combined = nansum(cat(5,im_stim_latemove_miss_avg_combined,im_stim_latemove_miss_avg),5); 
    n_conditions = sum(cat(5,n_conditions,any(any(any(im_stim_latemove_miss_avg,1),2),3)),5);
    AP_print_progress_fraction(curr_animal,length(animals));
end
im_stim_latemove_miss_avg_combined = bsxfun(@rdivide,im_stim_latemove_miss_avg_combined,n_conditions);
save_fn = [data_path filesep 'im_stim_latemove_miss_avg_combined.mat'];
save(save_fn,'im_stim_latemove_miss_avg_combined','-v7.3');
disp('Done')

%% Make batch widefield choiceworld mean (move aligned)

clear all
disp('Early move hit');
animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
im_move_earlymove_hit_avg_combined = zeros(437,416,265,11);
n_conditions = [];
for curr_animal = 1:length(animals) 
    animal = animals{curr_animal};
    fn = [data_path filesep animal '_im_move_earlymove_hit_avg.mat'];
    load(fn);
    
    im_move_earlymove_hit_avg_combined = nansum(cat(5,im_move_earlymove_hit_avg_combined,im_move_earlymove_hit_avg),5); 
    n_conditions = sum(cat(5,n_conditions,any(any(any(im_move_earlymove_hit_avg,1),2),3)),5);
    AP_print_progress_fraction(curr_animal,length(animals));
end
im_move_earlymove_hit_avg_combined = bsxfun(@rdivide,im_move_earlymove_hit_avg_combined,n_conditions);
save_fn = [data_path filesep 'im_move_earlymove_hit_avg_combined.mat'];
save(save_fn,'im_move_earlymove_hit_avg_combined','-v7.3');
disp('Done')

clear all
disp('Late move hit');
animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
im_move_latemove_hit_avg_combined = zeros(437,416,124,11);
n_conditions = [];
for curr_animal = 1:length(animals) 
    animal = animals{curr_animal};
    fn = [data_path filesep animal '_im_move_latemove_hit_avg.mat'];
    load(fn);
    
    im_move_latemove_hit_avg_combined = nansum(cat(5,im_move_latemove_hit_avg_combined,im_move_latemove_hit_avg),5); 
    n_conditions = sum(cat(5,n_conditions,any(any(any(im_move_latemove_hit_avg,1),2),3)),5);
    AP_print_progress_fraction(curr_animal,length(animals));
end
im_move_latemove_hit_avg_combined = bsxfun(@rdivide,im_move_latemove_hit_avg_combined,n_conditions);
save_fn = [data_path filesep 'im_move_latemove_hit_avg_combined.mat'];
save(save_fn,'im_move_latemove_hit_avg_combined','-v7.3');
disp('Done')

clear all
disp('Early move miss');
animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
im_move_earlymove_miss_avg_combined = zeros(437,416,124,11);
n_conditions = [];
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    fn = [data_path filesep animal '_im_move_earlymove_miss_avg.mat'];
    load(fn);
    
    im_move_earlymove_miss_avg_combined = nansum(cat(5,im_move_earlymove_miss_avg_combined,im_move_earlymove_miss_avg),5); 
    n_conditions = sum(cat(5,n_conditions,any(any(any(im_move_earlymove_miss_avg,1),2),3)),5);
    AP_print_progress_fraction(curr_animal,length(animals));
end
im_move_earlymove_miss_avg_combined = bsxfun(@rdivide,im_move_earlymove_miss_avg_combined,n_conditions);
save_fn = [data_path filesep 'im_move_earlymove_miss_avg_combined.mat'];
save(save_fn,'im_move_earlymove_miss_avg_combined','-v7.3');
disp('Done')

clear all
disp('Late move miss');
animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
im_move_latemove_miss_avg_combined = zeros(437,416,124,11);
n_conditions = [];
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    fn = [data_path filesep animal '_im_move_latemove_miss_avg.mat'];
    load(fn);
    
    im_move_latemove_miss_avg_combined = nansum(cat(5,im_move_latemove_miss_avg_combined,im_move_latemove_miss_avg),5); 
    n_conditions = sum(cat(5,n_conditions,any(any(any(im_move_latemove_miss_avg,1),2),3)),5);
    AP_print_progress_fraction(curr_animal,length(animals));
end
im_move_latemove_miss_avg_combined = bsxfun(@rdivide,im_move_latemove_miss_avg_combined,n_conditions);
save_fn = [data_path filesep 'im_move_latemove_miss_avg_combined.mat'];
save(save_fn,'im_move_latemove_miss_avg_combined','-v7.3');
disp('Done')

%% Load and process widefield choiceworld mean (stim aligned)

surround_window = [-0.5,2];
upsample_rate = 3;

framerate = 35.2;
surround_samplerate = 1/(framerate*upsample_rate);
t_surround = surround_window(1):surround_samplerate:surround_window(2);
t_df = conv2(t_surround,[1,1]/2,'valid');
conditions = [-1,-0.5,-0.25,-0.125,-0.06,0,0.06,0.125,0.25,0.5,1];

early_hit_data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
fn = [early_hit_data_path filesep 'im_stim_earlymove_hit_avg_combined'];
load(fn);

late_hit_data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
fn = [late_hit_data_path filesep 'im_stim_latemove_hit_avg_combined'];
load(fn);

early_miss_data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
fn = [early_miss_data_path filesep 'im_stim_earlymove_miss_avg_combined'];
load(fn);

late_miss_data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
fn = [late_miss_data_path filesep 'im_stim_latemove_miss_avg_combined'];
load(fn);

% ddf_earlymove_hit = im_stim_earlymove_hit_avg_combined;
ddf_earlymove_hit = diff(im_stim_earlymove_hit_avg_combined,[],3);
ddf_earlymove_hit(ddf_earlymove_hit < 0) = 0;
clear im_stim_earlymove_hit_avg_combined

% ddf_latemove_hit = im_stim_latemove_hit_avg_combined;
ddf_latemove_hit = diff(im_stim_latemove_hit_avg_combined,[],3);
ddf_latemove_hit(ddf_latemove_hit < 0) = 0;
clear im_stim_latemove_hit_avg_combined

% ddf_earlymove_miss = im_stim_earlymove_miss_avg_combined;
ddf_earlymove_miss = diff(im_stim_earlymove_miss_avg_combined,[],3);
ddf_earlymove_miss(ddf_earlymove_miss < 0) = 0;
clear im_stim_earlymove_miss_avg_combined

% ddf_latemove_miss = im_stim_latemove_miss_avg_combined;
ddf_latemove_miss = diff(im_stim_latemove_miss_avg_combined,[],3);
ddf_latemove_miss(ddf_latemove_miss < 0) = 0;
clear im_stim_latemove_miss_avg_combined

% % Get conditions which have data in all cases
% % (not used at the moment - all cases have data? how??)
% use_data_grid = ...
%     +squeeze(~any(isnan(cat(5, ...
%     ddf_earlymove_hit(:,:,1,7:end), ...
%     ddf_earlymove_hit(:,:,1,5:-1:1), ...
%     ddf_earlymove_miss(:,:,1,7:end), ...
%     ddf_earlymove_miss(:,:,1,5:-1:1), ...
%     ddf_latemove_hit(:,:,1,7:end), ...
%     ddf_latemove_hit(:,:,1,5:-1:1), ...
%     ddf_latemove_miss(:,:,1,7:end), ...
%     ddf_latemove_miss(:,:,1,5:-1:1))),5));
% use_data_grid(~logical(use_data_grid)) = NaN;

% Plot contrast response functions across the cortex
r_earlymove = nan(size(ddf_earlymove_hit,1),size(ddf_earlymove_hit,2),size(ddf_earlymove_hit,3));
l_earlymove = nan(size(ddf_earlymove_hit,1),size(ddf_earlymove_hit,2),size(ddf_earlymove_hit,3));

r_latemove = nan(size(ddf_latemove_hit,1),size(ddf_latemove_hit,2),size(ddf_latemove_hit,3));
l_latemove = nan(size(ddf_latemove_hit,1),size(ddf_latemove_hit,2),size(ddf_latemove_hit,3));

for curr_frame = 1:size(ddf_earlymove_hit,3);
    curr_r = gpuArray(reshape(squeeze(ddf_earlymove_hit(:,:,curr_frame,6:end)),[],6));
    curr_r = bsxfun(@minus,curr_r,mean(curr_r,2));
    r_fit = gather(curr_r/[1:6]);
    r_earlymove(:,:,curr_frame) = reshape(r_fit,size(ddf_earlymove_hit,1),size(ddf_earlymove_hit,2));
    
    curr_l = gpuArray(reshape(squeeze(ddf_earlymove_hit(:,:,curr_frame,6:-1:1)),[],6));
    curr_l = bsxfun(@minus,curr_l,mean(curr_l,2));
    l_fit = gather(curr_l/[1:6]);
    l_earlymove(:,:,curr_frame) = reshape(l_fit,size(ddf_earlymove_hit,1),size(ddf_earlymove_hit,2));
    
%     curr_r = gpuArray(reshape(squeeze(ddf_latemove_hit(:,:,curr_frame,6:end)),[],6));
%     curr_r = bsxfun(@minus,curr_r,mean(curr_r,2));
%     r_fit = gather(curr_r/[1:6]);
%     r_latemove(:,:,curr_frame) = reshape(r_fit,size(ddf_latemove_hit,1),size(ddf_latemove_hit,2));
%     
%     curr_l = gpuArray(reshape(squeeze(ddf_latemove_hit(:,:,curr_frame,6:-1:1)),[],6));
%     curr_l = bsxfun(@minus,curr_l,mean(curr_l,2));
%     l_fit = gather(curr_l/[1:6]);
%     l_latemove(:,:,curr_frame) = reshape(l_fit,size(ddf_latemove_hit,1),size(ddf_latemove_hit,2));
    
    AP_print_progress_fraction(curr_frame,size(ddf_earlymove_hit,3));
end

AP_image_scroll([r_earlymove-l_earlymove],t_surround);
caxis([-max(abs(caxis)),max(abs(caxis))]);
axis image; colormap(colormap_BlueWhiteRed);

% Plot stim/decision difference
decision_earlymove_diff_r = ddf_earlymove_hit(:,:,:,8)-ddf_earlymove_miss(:,:,:,8);
stim_earlymove_diff_r = ddf_earlymove_hit(:,:,:,8)-ddf_earlymove_miss(:,:,:,4);
stim_decision_earlymove_diff_r = nanmean(abs(stim_earlymove_diff_r)-abs(decision_earlymove_diff_r),4);

decision_earlymove_diff_l = ddf_earlymove_hit(:,:,:,4)-ddf_earlymove_miss(:,:,:,4);
stim_earlymove_diff_l = ddf_earlymove_hit(:,:,:,4)-ddf_earlymove_miss(:,:,:,8);
stim_decision_earlymove_diff_l = nanmean(abs(stim_earlymove_diff_l)-abs(decision_earlymove_diff_l),4);

stim_decision_earlymove_mean = ...
    nanmean(cat(4,stim_decision_earlymove_diff_l,stim_decision_earlymove_diff_r),4);

AP_image_scroll(stim_decision_earlymove_mean,t_surround);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(colormap_RedWhiteBlue);

% Get timing of activity

% (time of max or com, but bad because often multiphasic)
t_use = t_surround < 0.5;
ddf_time_com = squeeze(bsxfun(@rdivide,sum(bsxfun(@times, ...
    ddf_earlymove_hit(:,:,t_use,:),permute(t_surround(t_use),[1,3,2])),3), ...
    sum(ddf_earlymove_hit(:,:,t_use,:),3)));

t_use = t_surround > 0.5 & t_surround < 1;
ddf_time_com = squeeze(bsxfun(@rdivide,sum(bsxfun(@times, ...
    ddf_latemove_hit(:,:,t_use,:),permute(t_surround(t_use),[1,3,2])),3), ...
    sum(ddf_latemove_hit(:,:,t_use,:),3)));

% (first time above baseline)
t_baseline = t_surround < 0;
t_use = t_surround > 0 & t_surround < 0.5;
t_use_vals = t_surround(t_use);
ddf_baseline = nanmean(ddf_earlymove_hit(:,:,t_baseline,:),3);
ddf_baseline_mean = nanmean(ddf_baseline,4);
[~,first_rise] = max(bsxfun(@gt,ddf_earlymove_hit(:,:,t_use,:),ddf_baseline_mean*10),[],3);
first_rise_t = squeeze(t_use_vals(first_rise));

% Get traces for all pre-drawn ROIs
% Load widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = length(wf_roi);

figure('Name','Early hit');
for curr_roi = 1:n_rois   
    
    curr_mask_l = wf_roi(curr_roi,1).mask;
    curr_mask_r = curr_mask_l-wf_roi(curr_roi,2).mask;
    
    curr_traces_l = squeeze(...
        sum(sum(bsxfun(@times,ddf_earlymove_hit,curr_mask_l),1),2)./sum(curr_mask_l(:) ~= 0));
    curr_traces_r = squeeze(...
        sum(sum(bsxfun(@times,ddf_earlymove_hit,curr_mask_r),1),2)./sum(curr_mask_r(:) ~= 0));
    
    subplot(2,n_rois,curr_roi); hold on;
    set(gca,'ColorOrder',colormap_BlueWhiteRed((length(conditions)-1)/2));
    plot(t_df,curr_traces_l','linewidth',2);
    plot(t_df,curr_traces_l(:,conditions == 0),'k','linewidth',1);
    axis tight;
    line([0,0],ylim,'color','k');
    xlabel('Time from stim onset (s)')
    ylabel('\Delta\DeltaF/F');
    title(wf_roi(curr_roi,1).area);
    
    subplot(2,n_rois,n_rois+curr_roi); hold on;
    set(gca,'ColorOrder',colormap_BlueWhiteRed((length(conditions)-1)/2));
    plot(t_df,curr_traces_r','linewidth',2);
    plot(t_df,curr_traces_r(:,conditions == 0),'k','linewidth',1);
    axis tight;
    line([0,0],ylim,'color','k');
    xlabel('Time from stim onset (s)')
    ylabel('\Delta\DeltaF/F');
    title([wf_roi(curr_roi,1).area '-' wf_roi(curr_roi,2).area]);
    
end


%% Load and process widefield choiceworld mean (move aligned)

surround_window = [-0.5,2];
upsample_rate = 3;

framerate = 35.2;
surround_samplerate = 1/(framerate*upsample_rate);
t_surround = surround_window(1):surround_samplerate:surround_window(2);
t_df = conv2(t_surround,[1,1]/2,'valid');
conditions = [-1,-0.5,-0.25,-0.125,-0.06,0,0.06,0.125,0.25,0.5,1];

early_hit_data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
fn = [early_hit_data_path filesep 'im_move_earlymove_hit_avg_combined'];
load(fn);

late_hit_data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
fn = [late_hit_data_path filesep 'im_move_latemove_hit_avg_combined'];
load(fn);

early_miss_data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
fn = [early_miss_data_path filesep 'im_move_earlymove_miss_avg_combined'];
load(fn);

late_miss_data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
fn = [late_miss_data_path filesep 'im_move_latemove_miss_avg_combined'];
load(fn);

% ddf_earlymove_hit = im_move_earlymove_hit_avg_combined;
ddf_earlymove_hit = diff(im_move_earlymove_hit_avg_combined,[],3);
ddf_earlymove_hit(ddf_earlymove_hit < 0) = 0;
clear im_move_earlymove_hit_avg_combined

% ddf_latemove_hit = im_move_latemove_hit_avg_combined;
ddf_latemove_hit = diff(im_move_latemove_hit_avg_combined,[],3);
ddf_latemove_hit(ddf_latemove_hit < 0) = 0;
clear im_move_latemove_hit_avg_combined

% ddf_earlymove_miss = im_move_earlymove_miss_avg_combined;
ddf_earlymove_miss = diff(im_move_earlymove_miss_avg_combined,[],3);
ddf_earlymove_miss(ddf_earlymove_miss < 0) = 0;
clear im_move_earlymove_miss_avg_combined

% ddf_latemove_miss = im_move_latemove_miss_avg_combined;
ddf_latemove_miss = diff(im_move_latemove_miss_avg_combined,[],3);
ddf_latemove_miss(ddf_latemove_miss < 0) = 0;
clear im_move_latemove_miss_avg_combined

% % Get conditions which have data in all cases
% % (not used at the moment - all cases have data? how??)
% use_data_grid = ...
%     +squeeze(~any(isnan(cat(5, ...
%     ddf_earlymove_hit(:,:,1,7:end), ...
%     ddf_earlymove_hit(:,:,1,5:-1:1), ...
%     ddf_earlymove_miss(:,:,1,7:end), ...
%     ddf_earlymove_miss(:,:,1,5:-1:1), ...
%     ddf_latemove_hit(:,:,1,7:end), ...
%     ddf_latemove_hit(:,:,1,5:-1:1), ...
%     ddf_latemove_miss(:,:,1,7:end), ...
%     ddf_latemove_miss(:,:,1,5:-1:1))),5));
% use_data_grid(~logical(use_data_grid)) = NaN;

% Plot contrast response functions across the cortex
r_earlymove = nan(size(ddf_earlymove_hit,1),size(ddf_earlymove_hit,2),size(ddf_earlymove_hit,3));
l_earlymove = nan(size(ddf_earlymove_hit,1),size(ddf_earlymove_hit,2),size(ddf_earlymove_hit,3));

r_latemove = nan(size(ddf_latemove_hit,1),size(ddf_latemove_hit,2),size(ddf_latemove_hit,3));
l_latemove = nan(size(ddf_latemove_hit,1),size(ddf_latemove_hit,2),size(ddf_latemove_hit,3));

for curr_frame = 1:size(ddf_earlymove_hit,3);
    curr_r = gpuArray(reshape(squeeze(ddf_earlymove_hit(:,:,curr_frame,6:end)),[],6));
    curr_r = bsxfun(@minus,curr_r,mean(curr_r,2));
    r_fit = gather(curr_r/[1:6]);
    r_earlymove(:,:,curr_frame) = reshape(r_fit,size(ddf_earlymove_hit,1),size(ddf_earlymove_hit,2));
    
    curr_l = gpuArray(reshape(squeeze(ddf_earlymove_hit(:,:,curr_frame,6:-1:1)),[],6));
    curr_l = bsxfun(@minus,curr_l,mean(curr_l,2));
    l_fit = gather(curr_l/[1:6]);
    l_earlymove(:,:,curr_frame) = reshape(l_fit,size(ddf_earlymove_hit,1),size(ddf_earlymove_hit,2));
    
%     curr_r = gpuArray(reshape(squeeze(ddf_latemove_hit(:,:,curr_frame,6:end)),[],6));
%     curr_r = bsxfun(@minus,curr_r,mean(curr_r,2));
%     r_fit = gather(curr_r/[1:6]);
%     r_latemove(:,:,curr_frame) = reshape(r_fit,size(ddf_latemove_hit,1),size(ddf_latemove_hit,2));
%     
%     curr_l = gpuArray(reshape(squeeze(ddf_latemove_hit(:,:,curr_frame,6:-1:1)),[],6));
%     curr_l = bsxfun(@minus,curr_l,mean(curr_l,2));
%     l_fit = gather(curr_l/[1:6]);
%     l_latemove(:,:,curr_frame) = reshape(l_fit,size(ddf_latemove_hit,1),size(ddf_latemove_hit,2));
    
    AP_print_progress_fraction(curr_frame,size(ddf_earlymove_hit,3));
end

AP_image_scroll(r_earlymove,t_df);
caxis([-max(abs(caxis)),max(abs(caxis))]);
axis image; colormap(colormap_BlueWhiteRed);

% Plot stim/decision difference
decision_earlymove_diff_r = ddf_earlymove_hit(:,:,:,8)-ddf_earlymove_miss(:,:,:,8);
stim_earlymove_diff_r = ddf_earlymove_hit(:,:,:,8)-ddf_earlymove_miss(:,:,:,4);
stim_decision_earlymove_diff_r = nanmean(abs(stim_earlymove_diff_r)-abs(decision_earlymove_diff_r),4);

decision_earlymove_diff_l = ddf_earlymove_hit(:,:,:,4)-ddf_earlymove_miss(:,:,:,4);
stim_earlymove_diff_l = ddf_earlymove_hit(:,:,:,4)-ddf_earlymove_miss(:,:,:,8);
stim_decision_earlymove_diff_l = nanmean(abs(stim_earlymove_diff_l)-abs(decision_earlymove_diff_l),4);

stim_decision_earlymove_mean = ...
    nanmean(cat(4,stim_decision_earlymove_diff_l,stim_decision_earlymove_diff_r),4);

AP_image_scroll(stim_decision_earlymove_mean,t_surround);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(colormap_RedWhiteBlue);

% Get timing of activity

% (time of max or com, but bad because often multiphasic)
t_use = t_surround < 0.5;
ddf_time_com = squeeze(bsxfun(@rdivide,sum(bsxfun(@times, ...
    ddf_earlymove_hit(:,:,t_use,:),permute(t_surround(t_use),[1,3,2])),3), ...
    sum(ddf_earlymove_hit(:,:,t_use,:),3)));

t_use = t_surround > 0.5 & t_surround < 1;
ddf_time_com = squeeze(bsxfun(@rdivide,sum(bsxfun(@times, ...
    ddf_latemove_hit(:,:,t_use,:),permute(t_surround(t_use),[1,3,2])),3), ...
    sum(ddf_latemove_hit(:,:,t_use,:),3)));

% (first time above baseline)
t_baseline = t_surround < 0;
t_use = t_surround > 0 & t_surround < 0.5;
t_use_vals = t_surround(t_use);
ddf_baseline = nanmean(ddf_earlymove_hit(:,:,t_baseline,:),3);
ddf_baseline_mean = nanmean(ddf_baseline,4);
[~,first_rise] = max(bsxfun(@gt,ddf_earlymove_hit(:,:,t_use,:),ddf_baseline_mean*10),[],3);
first_rise_t = squeeze(t_use_vals(first_rise));

% Get traces for all pre-drawn ROIs
% Load widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = length(wf_roi);

figure('Name','Early hit');
for curr_roi = 1:n_rois   
    curr_mask = wf_roi(curr_roi).mask;
    curr_traces = squeeze(...
        sum(sum(bsxfun(@times,ddf_earlymove_hit,curr_mask),1),2)./sum(curr_mask(:)));
    
    subplot(2,7,curr_roi); hold on;
    set(gca,'ColorOrder',colormap_BlueWhiteRed((length(conditions)-1)/2));
    plot(t_df,curr_traces','linewidth',2);
    plot(t_df,curr_traces(:,conditions == 0),'k','linewidth',1);
    line([0,0],ylim,'color','k');
    xlabel('Time from move onset (s)')
    ylabel('\Delta\DeltaF/F');
    title(wf_roi(curr_roi).area);
end

figure('Name','Late hit');
for curr_roi = 1:n_rois   
    curr_mask = wf_roi(curr_roi).mask;
    curr_traces = squeeze(...
        sum(sum(bsxfun(@times,ddf_latemove_hit,curr_mask),1),2)./sum(curr_mask(:)));
    
    subplot(2,4,curr_roi); hold on;
    set(gca,'ColorOrder',colormap_BlueWhiteRed((length(conditions)-1)/2));
    plot(t_df,curr_traces','linewidth',2);
    plot(t_df,curr_traces(:,conditions == 0),'k','linewidth',1);
    line([0,0],ylim,'color','k');
    xlabel('Time from move onset (s)')
    ylabel('\Delta\DeltaF/F');
    title(wf_roi(curr_roi).area);
end


%% Load and process striatal MUA during choiceworld

data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
mua_fn = [data_path filesep 'mua_choiceworld'];
load(mua_fn);

conditions = [-1,-0.5,-0.25,-0.125,-0.06,0,0.06,0.125,0.25,0.5,1];
n_conditions = length(conditions);
n_depths = size(batch_vars(1).mua_stim_earlymove_hit,1);

raster_window = [-0.5,3];
psth_bin_size = 0.001;
t = raster_window(1):psth_bin_size:raster_window(2);
t_bins = t(1:end-1) + diff(t);

% Loop through all MUA fields, normalize and combine
smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

t_baseline = t_bins < 0;

softnorm = 1;

% (get baseline from all stim-aligned MUAs)
mua_baseline = cell(size(batch_vars));
for curr_animal = 1:length(batch_vars)
   mua_baseline{curr_animal} = nanmean(nanmean(cat(5, ...
       batch_vars(curr_animal).mua_stim_earlymove_hit(:,t_baseline,:,:), ...
       batch_vars(curr_animal).mua_stim_latemove_hit(:,t_baseline,:,:), ...
       batch_vars(curr_animal).mua_stim_earlymove_miss(:,t_baseline,:,:), ...
       batch_vars(curr_animal).mua_stim_latemove_miss(:,t_baseline,:,:)),2),5);   
end

mua_fieldnames = fieldnames(batch_vars);
mua = struct;
for curr_mua = mua_fieldnames'
    curr_field = curr_mua{:};
    curr_mua_smoothed = cellfun(@(x) convn(x,smWin,'same'),{batch_vars(:).(curr_field)},'uni',false);    
    curr_mua_norm = cellfun(@(x,baseline) bsxfun(@rdivide,x,baseline+softnorm),curr_mua_smoothed,mua_baseline,'uni',false);
    curr_mua_mean = cellfun(@(x) nanmean(x,4),curr_mua_norm,'uni',false);
    curr_mua_combined = nanmean(cat(4,curr_mua_mean{:}),4);
    
    mua.(curr_field(5:end)) = curr_mua_combined;    
end

% (to plot all early move contrasts across all depths)
% for plot_depth = 1:n_depths
%     figure; hold on;
%     trace_spacing = 5;
%     p1 = AP_stackplot(squeeze(mua.stim_earlymove_hit(plot_depth,:,:)),t_bins,trace_spacing,false,'k',conditions);
%     p2 = AP_stackplot(squeeze(mua.stim_earlymove_miss(plot_depth,:,:)),t_bins,trace_spacing,false,'r',conditions);
%     xlabel('Time from stim onset');
%     ylabel(['MUA depth ' num2str(plot_depth)])
%     legend([p1(1),p2(1)],{'Hit','Miss'});
%     line([0,0],ylim,'color','k','linestyle','--');
% end

% (to plot all late move contrasts across all depths)
% for plot_depth = 1:n_depths
%     figure; hold on;
%     trace_spacing = 5;
%     p1 = AP_stackplot(squeeze(mua.stim_latemove_hit(plot_depth,:,:)),t_bins,trace_spacing,false,'k',conditions);
%     p2 = AP_stackplot(squeeze(mua.stim_latemove_miss(plot_depth,:,:)),t_bins,trace_spacing,false,'r',conditions);
%     xlabel('Time from stim onset');
%     ylabel(['MUA depth ' num2str(plot_depth)])
%     legend([p1(1),p2(1)],{'Hit','Miss'});
%     line([0,0],ylim,'color','k','linestyle','--');
% end

% Get conditions which have data in all cases
use_data_grid = ...
    +squeeze(~any(isnan(cat(4,mua.stim_earlymove_hit(:,:,7:end), ...
    mua.stim_earlymove_hit(:,:,5:-1:1), ...
    mua.stim_earlymove_miss(:,:,7:end), ...
    mua.stim_earlymove_miss(:,:,5:-1:1), ...
    mua.stim_latemove_hit(:,:,7:end), ...
    mua.stim_latemove_hit(:,:,5:-1:1), ...
    mua.stim_latemove_miss(:,:,7:end), ...
    mua.stim_latemove_miss(:,:,5:-1:1))),4));
use_data_grid(~logical(use_data_grid)) = NaN;
% (TO ONLY USE 12.5% CONDITION)
use_data_grid = bsxfun(@times,use_data_grid,permute([NaN,1,NaN,NaN,NaN],[1,3,2]));

% Get average activity within decision
move_right_earlymove_hit = nanmean(mua.stim_earlymove_hit(:,:,5:-1:1).*use_data_grid,3);
move_right_earlymove_miss = nanmean(mua.stim_earlymove_miss(:,:,7:end).*use_data_grid,3);

move_left_earlymove_hit = nanmean(mua.stim_earlymove_hit(:,:,7:end).*use_data_grid,3);
move_left_earlymove_miss = nanmean(mua.stim_earlymove_miss(:,:,5:-1:1).*use_data_grid,3);

move_right_latemove_hit = nanmean(mua.stim_latemove_hit(:,:,5:-1:1).*use_data_grid,3);
move_right_latemove_miss = nanmean(mua.stim_latemove_miss(:,:,7:end).*use_data_grid,3);

move_left_latemove_hit = nanmean(mua.stim_latemove_hit(:,:,7:end).*use_data_grid,3);
move_left_latemove_miss = nanmean(mua.stim_latemove_miss(:,:,5:-1:1).*use_data_grid,3);

% Plot overlay all conditions for all depth
curr_data = cell2mat(permute(arrayfun(@(x) ...
    mat2gray(squeeze(mua.move_latemove_hit(x,:,:))),1:n_depths,'uni',false),[1,3,2]));

figure; hold on;
col = colormap_BlueWhiteRed((length(conditions)-1)/2);
for curr_cond = 1:n_conditions
    AP_stackplot(squeeze(curr_data(:,curr_cond,:)),t_bins,1,false,col(curr_cond,:));
end

% (plot decision by depth)
figure; trace_spacing = 6;

subplot(1,2,1);hold on;
p1 = AP_stackplot(move_right_earlymove_hit',t_bins,trace_spacing,false,'b',1:n_depths);
p2 = AP_stackplot(move_right_earlymove_miss',t_bins,trace_spacing,false,[0.7,0.7,1],1:n_depths);
p3 = AP_stackplot(move_left_earlymove_hit',t_bins,trace_spacing,false,'r',1:n_depths);
p4 = AP_stackplot(move_left_earlymove_miss',t_bins,trace_spacing,false,[1,0.7,0.7],1:n_depths);
ylim([trace_spacing,(n_depths+2)*trace_spacing]);
line([0,0],ylim,'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
xlabel('Time from stim onset (s)');
ylabel('MUA depth');
legend([p1(1),p2(1),p3(1),p4(1)],{'Move right hit','Move right miss','Move left hit','Move left miss'});
title('Early move');

subplot(1,2,2);hold on;
p1 = AP_stackplot(move_right_latemove_hit',t_bins,trace_spacing,false,'b',1:n_depths);
p2 = AP_stackplot(move_right_latemove_miss',t_bins,trace_spacing,false,[0.7,0.7,1],1:n_depths);
p3 = AP_stackplot(move_left_latemove_hit',t_bins,trace_spacing,false,'r',1:n_depths);
p4 = AP_stackplot(move_left_latemove_miss',t_bins,trace_spacing,false,[1,0.7,0.7],1:n_depths);
ylim([trace_spacing,(n_depths+2)*trace_spacing]);
line([0,0],ylim,'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
xlabel('Time from stim onset (s)');
ylabel('MUA depth');
legend([p1(1),p2(1),p3(1),p4(1)],{'Move right hit','Move right miss','Move left hit','Move left miss'});
title('Late move');

% (plot decision difference by depth)

% this probably isn't the way to do it, will include lopsided conditions
% decision_earlymove_diff = move_left_earlymove_hit-move_right_earlymove_miss;
% stim_earlymove_diff = move_left_earlymove_hit-move_left_earlymove_miss;
% stim_decision_earlymove_diff = stim_earlymove_diff-decision_earlymove_diff;
% 
% decision_latemove_diff = move_left_latemove_hit-move_right_latemove_miss;
% stim_latemove_diff = move_left_latemove_hit-move_left_latemove_miss;
% stim_decision_latemove_diff = stim_latemove_diff-decision_latemove_diff;

% this way only uses conditions that are represented in all used cases
stim_earlymove_diff = (mua.stim_earlymove_hit(:,:,7:end)-mua.stim_earlymove_miss(:,:,5:-1:1)).*use_data_grid;
decision_earlymove_diff = (mua.stim_earlymove_hit(:,:,7:end)-mua.stim_earlymove_miss(:,:,7:end)).*use_data_grid;
stim_decision_earlymove_diff = nanmean(abs(stim_earlymove_diff)-abs(decision_earlymove_diff),3);

stim_latemove_diff = (mua.stim_latemove_hit(:,:,7:end)-mua.stim_latemove_miss(:,:,5:-1:1)).*use_data_grid;
decision_latemove_diff = (mua.stim_latemove_hit(:,:,7:end)-mua.stim_latemove_miss(:,:,7:end)).*use_data_grid;
stim_decision_latemove_diff = nanmean(abs(stim_latemove_diff)-abs(decision_latemove_diff),3);

figure; trace_spacing = 10;
subplot(1,2,1); hold on;
pos_plot = ones(n_depths,length(t_bins)); pos_plot(stim_decision_earlymove_diff < 0) = NaN;
neg_plot = ones(n_depths,length(t_bins)); neg_plot(stim_decision_earlymove_diff > 0) = NaN;
p1 = AP_stackplot(abs(nanmean(stim_earlymove_diff,3))',t_bins,trace_spacing,false,[0.5,0.5,0.5],1:n_depths);
p2 = AP_stackplot(-abs(nanmean(decision_earlymove_diff,3))',t_bins,trace_spacing,false,[0.5,0.5,0.5],1:n_depths);
p3 = AP_stackplot(stim_decision_earlymove_diff'.*pos_plot',t_bins,trace_spacing,false,'b',1:n_depths);
p4 = AP_stackplot(stim_decision_earlymove_diff'.*neg_plot',t_bins,trace_spacing,false,'r',1:n_depths);
ylim([trace_spacing/2,(n_depths+0.5)*trace_spacing]);
line([0,0],ylim,'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
legend([p1(1),p3(1),p4(1)],{'Independent','Stim','Movement'})
xlabel('Time from stim onset (s)');
ylabel('MUA depth');
title('Early move');

subplot(1,2,2); hold on;
pos_plot = ones(n_depths,length(t_bins)); pos_plot(stim_decision_latemove_diff < 0) = NaN;
neg_plot = ones(n_depths,length(t_bins)); neg_plot(stim_decision_latemove_diff > 0) = NaN;
p1 = AP_stackplot(abs(nanmean(stim_latemove_diff,3))',t_bins,trace_spacing,false,[0.5,0.5,0.5],1:n_depths);
p2 = AP_stackplot(-abs(nanmean(decision_latemove_diff,3))',t_bins,trace_spacing,false,[0.5,0.5,0.5],1:n_depths);
p3 = AP_stackplot(stim_decision_latemove_diff'.*pos_plot',t_bins,trace_spacing,false,'b',1:n_depths);
p4 = AP_stackplot(stim_decision_latemove_diff'.*neg_plot',t_bins,trace_spacing,false,'r',1:n_depths);
ylim([trace_spacing/2,(n_depths+0.5)*trace_spacing]);
line([0,0],ylim,'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
legend([p1(1),p3(1),p4(1)],{'Independent','Stim','Movement'})
xlabel('Time from stim onset (s)');
ylabel('MUA depth');
title('Late move');

% (plot all traces from a depth on top of each other - can see decision)
plot_depth = 1;
figure; 

a1 = subplot(1,2,1); hold on
set(gca,'ColorOrder',colormap_BlueWhiteRed((n_conditions-1)/2));
p1 = plot(t_bins,squeeze(mua.stim_earlymove_hit(plot_depth,:,:)),'linewidth',2);
p2 = plot(t_bins,squeeze(mua.stim_earlymove_miss(plot_depth,:,:)),'linewidth',2,'linestyle','--');
legend([p1(1),p2(2)],{'Hit','Miss'})
title('Early move');

a2 = subplot(1,2,2); hold on
set(gca,'ColorOrder',colormap_BlueWhiteRed((n_conditions-1)/2));
p1 = plot(t_bins,squeeze(mua.stim_latemove_hit(plot_depth,:,:)),'linewidth',2);
p2 = plot(t_bins,squeeze(mua.stim_latemove_miss(plot_depth,:,:)),'linewidth',2,'linestyle','--');
legend([p1(1),p2(2)],{'Hit','Miss'})
title('Late move');

linkaxes([a1,a2],'x');

% (to plot all depths for left/right contrast 1)
figure;
trace_spacing = 6;

subplot(1,2,1); hold on;
p1 = AP_stackplot(mua.stim_earlymove_hit(:,:,1)',t_bins,trace_spacing,false,'k',1:n_depths);
p2 = AP_stackplot(mua.stim_latemove_hit(:,:,1)',t_bins,trace_spacing,false,'r',1:n_depths);
ylim([trace_spacing,(n_depths+1)*trace_spacing])
line([0,0],ylim,'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
legend([p1(1),p2(1)],{'Early move','Late move'});
xlabel('Time from stim onset (s)')
ylabel('Striatum depth')
title('Left Contrast 1')

subplot(1,2,2); hold on;
p1 = AP_stackplot(mua.stim_earlymove_hit(:,:,11)',t_bins,trace_spacing,false,'k',1:n_depths);
p2 = AP_stackplot(mua.stim_latemove_hit(:,:,11)',t_bins,trace_spacing,false,'r',1:n_depths);
ylim([trace_spacing,(n_depths+1)*trace_spacing])
line([0,0],ylim,'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
legend([p1(1),p2(1)],{'Early move','Late move'});
xlabel('Time from stim onset (s)')
ylabel('Striatum depth')
title('Right Contrast 1')

% Get contrast tuning by time point for hit/miss
r_earlymove_hit = nan(6,length(t_bins));
l_earlymove_hit = nan(6,length(t_bins));

r_latemove_hit = nan(6,length(t_bins));
l_latemove_hit = nan(6,length(t_bins));
for curr_t = 1:length(t_bins)
    
    curr_r = squeeze(mua.stim_earlymove_hit(:,curr_t,6:end));
    curr_r = bsxfun(@minus,curr_r,mean(curr_r,2));
    r_fit = curr_r/[1:6];
    r_earlymove_hit(:,curr_t) = r_fit;
    
    curr_l = squeeze(mua.stim_earlymove_hit(:,curr_t,6:-1:1));
    curr_l = bsxfun(@minus,curr_l,mean(curr_l,2));
    l_fit = curr_l/[1:6];
    l_earlymove_hit(:,curr_t) = l_fit;
    
    curr_r = squeeze(mua.stim_latemove_hit(:,curr_t,6:end));
    curr_r = bsxfun(@minus,curr_r,mean(curr_r,2));
    r_fit = curr_r/[1:6];
    r_latemove_hit(:,curr_t) = r_fit;
    
    curr_l = squeeze(mua.stim_latemove_hit(:,curr_t,6:-1:1));
    curr_l = bsxfun(@minus,curr_l,mean(curr_l,2));
    l_fit = curr_l/[1:6];
    l_latemove_hit(:,curr_t) = l_fit;
    
end

trace_spacing = 0.3;

figure; 
subplot(1,3,1); hold on;
p1 = AP_stackplot(r_earlymove_hit',t_bins,trace_spacing,false,'k',1:n_depths);
p2 = AP_stackplot(l_earlymove_hit',t_bins,trace_spacing,false,'r',1:n_depths);
ylim([trace_spacing/2,n_depths*trace_spacing+trace_spacing/2])
ylabel('Contrast slope')
xlabel('Time from stim onset (s)');
line([0,0],ylim,'linestyle','--','color','k');
line([0.5,0.5],ylim,'color','k','linestyle','--');
legend([p1(1),p2(1)],{'Right stim','Left stim'});
title('Early move hit');

subplot(1,3,2); hold on;
p1 = AP_stackplot(r_latemove_hit',t_bins,trace_spacing,false,'k',1:n_depths);
p2 = AP_stackplot(l_latemove_hit',t_bins,trace_spacing,false,'r',1:n_depths);
ylim([trace_spacing/2,n_depths*trace_spacing+trace_spacing/2])
ylabel('Contrast slope')
xlabel('Time from stim onset (s)');
line([0,0],ylim,'linestyle','--','color','k');
line([0.5,0.5],ylim,'color','k','linestyle','--');
legend([p1(1),p2(1)],{'Right stim','Left stim'});
title('Late move hit');

subplot(1,3,3); hold on;
p1 = AP_stackplot(r_earlymove_hit',t_bins,trace_spacing,false,'k',1:n_depths);
p2 = AP_stackplot(r_latemove_hit',t_bins,trace_spacing,false,'r',1:n_depths);
ylim([trace_spacing/2,n_depths*trace_spacing+trace_spacing/2])
ylabel('Contrast slope')
xlabel('Time from stim onset (s)');
line([0,0],ylim,'linestyle','--','color','k');
line([0.5,0.5],ylim,'color','k','linestyle','--');
legend([p1(1),p2(1)],{'Early move','Late move'});
title('Right stimuli hit');

vis_t = t_bins > 0 & t_bins < 0.2;
vis_response_earlymove_hit = squeeze(max(mua.stim_earlymove_hit(:,vis_t,:),[],2));
vis_response_latemove_hit = squeeze(max(mua.stim_latemove_hit(:,vis_t,:),[],2));
vis_response_earlymove_miss = squeeze(max(mua.stim_earlymove_miss(:,vis_t,:),[],2));
vis_response_latemove_miss = squeeze(max(mua.stim_latemove_miss(:,vis_t,:),[],2));
figure; 
subplot(1,2,1); hold on
p1 = AP_stackplot(vis_response_earlymove_hit',conditions,2,false,'k',1:n_depths);
p2 = AP_stackplot(vis_response_earlymove_miss',conditions,2,false,'r',1:n_depths);
legend([p1(1),p2(1)],{'Hit','Miss'})
ylabel('Depth')
xlabel('Condition')
title('Early move')
axis tight
subplot(1,2,2); hold on
p1 = AP_stackplot(vis_response_latemove_hit',conditions,2,false,'k',1:n_depths);
p2 = AP_stackplot(vis_response_latemove_miss',conditions,2,false,'r',1:n_depths);
legend([p1(1),p2(1)],{'Hit','Miss'})
ylabel('Depth')
xlabel('Condition')
title('Late move')
axis tight


%% Load and process striatal MUA during passive

data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\passive'];

protocol = 'stimKalatsky';
% protocol = 'AP_choiceWorldStimPassive';

mua_fn = [data_path filesep 'mua_stim_' protocol];
load(mua_fn);

raster_window = [-0.5,5];
psth_bin_size = 0.001;
t = raster_window(1):psth_bin_size:raster_window(2);
t_bins = t(1:end-1) + diff(t);

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

t_baseline = t_bins < 0;

softnorm = 1;

use_animals = cellfun(@(x) ~isempty(x),{batch_vars(:).mua_stim});

mua_stim_smoothed = cellfun(@(x) convn(x,smWin,'same'),{batch_vars(use_animals).mua_stim},'uni',false);
mua_stim_norm = cellfun(@(x) bsxfun(@rdivide,x,nanmean(x(:,t_baseline,:,:),2)+softnorm),mua_stim_smoothed,'uni',false);
mua_stim_mean = cellfun(@(x) nanmean(x,4),mua_stim_norm,'uni',false);
mua_stim_combined = nanmean(cat(4,mua_stim_mean{:}),4);

trace_spacing = 2;
n_conditions = size(mua_stim_combined,3);
n_depths = size(mua_stim_combined,1);
col = copper(n_conditions);
figure; hold on;
p = gobjects(n_depths,n_conditions);
for curr_condition = 1:n_conditions
    p(:,curr_condition) = ...
        AP_stackplot(mua_stim_combined(:,:,curr_condition)',t_bins,trace_spacing,false,col(curr_condition,:));
end
axis tight;
ylabel('Spikes')
xlabel('Time from stim onset (s)');
line([0,0],ylim,'linestyle','--','color','k');
line([1,1],ylim,'linestyle','--','color','k');
line([2,2],ylim,'linestyle','--','color','k');
legend(p(1,:),cellfun(@num2str,num2cell(1:n_conditions),'uni',false));


%% Load and average wf ephys maps
error('different weights for different lambdas, need to normalize')

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_ephys';

protocol = 'vanillaChoiceworld';
% protocol = 'stimSparseNoiseUncorrAsync';
% protocol = 'stimKalatsky';
% protocol = 'AP_choiceWorldStimPassive';

map_fn = [data_path filesep 'wf_ephys_maps_' protocol];
load(map_fn);

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

r_px = cell(length(animals),1);
com = nan(437,416,length(animals));
weight = nan(437,416,length(animals));
explained_var = nan(6,length(animals));
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);       
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    days = {experiments.day};
    
    % Skip if this animal doesn't have this experiment
    if isempty(experiments)
        continue
    end
    
    r_px_aligned = AP_align_widefield(animal,days,batch_vars(curr_animal).r_px);
    r_px_com_aligned = AP_align_widefield(animal,days,batch_vars(curr_animal).r_px_com);
    r_px_weight_aligned = AP_align_widefield(animal,days,batch_vars(curr_animal).r_px_weight);
    explained_var_cat = horzcat(batch_vars(curr_animal).explained_var{:});
    
    r_px{curr_animal} = nanmean(r_px_aligned,5);
    com(:,:,curr_animal) = nanmean(r_px_com_aligned,3);
    weight(:,:,curr_animal) = nanmean(r_px_weight_aligned,3);
    explained_var(:,curr_animal) = nanmean(explained_var_cat,2);
    
    AP_print_progress_fraction(curr_animal,length(animals));
    
end

r_px_mean = nanmean(cell2mat(permute(cellfun(@(x) nanmean(x,5),r_px,'uni',false),[2,3,4,5,1])),5);
com_mean = nanmean(com,3);
weight_mean = nanmean(weight,3);

% Plot map
n_depth_groups = 6;
com_leeway = 1;
c_range = [1+com_leeway,n_depth_groups-com_leeway];
com_colored = ind2rgb(round(mat2gray(com_mean,c_range)*255),jet(255));

figure;
p = imagesc(com_colored);
axis image off
weight_norm = mat2gray(max(weight_mean,[],3),[0,double(prctile(reshape(max(weight_mean,[],3),[],1),95))]);
set(p,'AlphaData',weight_norm);

c = colorbar;
ylabel(c,'Depth (fraction)');
colormap(c,jet);
set(c,'YDir','reverse');
set(c,'YTick',linspace(0,1,diff(c_range)+1));
set(c,'YTickLabel',c_range(1):c_range(2));

AP_reference_outline('retinotopy','m');
AP_reference_outline('ccf_aligned','k');

title(protocol);

% Plot weight
figure;imagesc(weight_mean);
axis image off
colormap(hot); 
caxis([0,prctile(reshape(max(weight_mean,[],3),[],1),95)]);
AP_reference_outline('retinotopy','m');
AP_reference_outline('ccf_aligned','k');
title([protocol ' - weights']);

% Plot explained variance by depth
figure;
errorbar(nanmean(explained_var,2),nanstd(explained_var,[],2)./sqrt(sum(~isnan(explained_var),2)),'k','linewidth',2);
ylabel('Fraction explained variance');
xlabel('Depth');
axis tight
title(protocol);

%% Load and process cortex-predicted striatal MUA during choiceworld

data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
mua_fn = [data_path filesep 'mua_stim_choiceworld_pred'];
load(mua_fn);

% This smoothing window was got empirically (should be found optimally): 
% it find the MUA smoothing which the predicted spiking gives best
smooth_size = 20;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
% smWin = ones(1,10)/10;
% smWin = 1;

raster_window = [-0.5,3];
t_length = size(batch_vars(1).mua_stim_earlymove_hit,2);
t = linspace(raster_window(1),raster_window(2),t_length);

n_depths = size(batch_vars(1).mua_stim_earlymove_hit,1);
conditions = [-1,-0.5,-0.25,-0.125,-0.06,0,0.06,0.125,0.25,0.5,1];
n_conditions = length(conditions);

t_baseline = t < 0;

softnorm = 1;

mua_stim_earlymove_hit_norm = cellfun(@(x) ...
    bsxfun(@rdivide,x,nanmean(x(:,t_baseline,:,:),2)+softnorm),{batch_vars(:).mua_stim_earlymove_hit},'uni',false);
mua_stim_earlymove_hit_mean = cellfun(@(x) nanmean(x,4),mua_stim_earlymove_hit_norm,'uni',false);
mua_stim_earlymove_hit_combined = convn(nanmean(cat(4,mua_stim_earlymove_hit_mean{:}),4),smWin,'same');

mua_stim_earlymove_hit_pred_norm = cellfun(@(x) ...
    bsxfun(@rdivide,x,nanmean(x(:,t_baseline,:,:),2)+softnorm),{batch_vars(:).mua_stim_earlymove_hit_pred},'uni',false);
mua_stim_earlymove_hit_pred_mean = cellfun(@(x) nanmean(x,4),mua_stim_earlymove_hit_pred_norm,'uni',false);
mua_stim_earlymove_hit_pred_combined = nanmean(cat(4,mua_stim_earlymove_hit_pred_mean{:}),4);

mua_stim_latemove_hit_norm = cellfun(@(x) ...
    bsxfun(@rdivide,x,nanmean(x(:,t_baseline,:,:),2)+softnorm),{batch_vars(:).mua_stim_latemove_hit},'uni',false);
mua_stim_latemove_hit_mean = cellfun(@(x) nanmean(x,4),mua_stim_latemove_hit_norm,'uni',false);
mua_stim_latemove_hit_combined = convn(nanmean(cat(4,mua_stim_latemove_hit_mean{:}),4),smWin,'same');

mua_stim_latemove_hit_pred_norm = cellfun(@(x) ...
    bsxfun(@rdivide,x,nanmean(x(:,t_baseline,:,:),2)+softnorm),{batch_vars(:).mua_stim_latemove_hit_pred},'uni',false);
mua_stim_latemove_hit_pred_mean = cellfun(@(x) nanmean(x,4),mua_stim_latemove_hit_pred_norm,'uni',false);
mua_stim_latemove_hit_pred_combined = nanmean(cat(4,mua_stim_latemove_hit_pred_mean{:}),4);

mua_stim_earlymove_miss_norm = cellfun(@(x) ...
    bsxfun(@rdivide,x,nanmean(x(:,t_baseline,:,:),2)+softnorm),{batch_vars(:).mua_stim_earlymove_miss},'uni',false);
mua_stim_earlymove_miss_mean = cellfun(@(x) nanmean(x,4),mua_stim_earlymove_miss_norm,'uni',false);
mua_stim_earlymove_miss_combined = convn(nanmean(cat(4,mua_stim_earlymove_miss_mean{:}),4),smWin,'same');

mua_stim_earlymove_miss_pred_norm = cellfun(@(x) ...
    bsxfun(@rdivide,x,nanmean(x(:,t_baseline,:,:),2)+softnorm),{batch_vars(:).mua_stim_earlymove_miss_pred},'uni',false);
mua_stim_earlymove_miss_pred_mean = cellfun(@(x) nanmean(x,4),mua_stim_earlymove_miss_pred_norm,'uni',false);
mua_stim_earlymove_miss_pred_combined = nanmean(cat(4,mua_stim_earlymove_miss_pred_mean{:}),4);

mua_stim_latemove_miss_norm = cellfun(@(x) ...
    bsxfun(@rdivide,x,nanmean(x(:,t_baseline,:,:),2)+softnorm),{batch_vars(:).mua_stim_latemove_miss},'uni',false);
mua_stim_latemove_miss_mean = cellfun(@(x) nanmean(x,4),mua_stim_latemove_miss_norm,'uni',false);
mua_stim_latemove_miss_combined = convn(nanmean(cat(4,mua_stim_latemove_miss_mean{:}),4),smWin,'same');

mua_stim_latemove_miss_pred_norm = cellfun(@(x) ...
    bsxfun(@rdivide,x,nanmean(x(:,t_baseline,:,:),2)+softnorm),{batch_vars(:).mua_stim_latemove_miss_pred},'uni',false);
mua_stim_latemove_miss_pred_mean = cellfun(@(x) nanmean(x,4),mua_stim_latemove_miss_pred_norm,'uni',false);
mua_stim_latemove_miss_pred_combined = nanmean(cat(4,mua_stim_latemove_miss_pred_mean{:}),4);

% Get conditions which have data in all cases
use_data_grid = ...
    +squeeze(~any(isnan(cat(4,mua_stim_earlymove_hit_combined(:,:,7:end), ...
    mua_stim_earlymove_hit_combined(:,:,5:-1:1), ...
    mua_stim_earlymove_miss_combined(:,:,7:end), ...
    mua_stim_earlymove_miss_combined(:,:,5:-1:1), ...
    mua_stim_latemove_hit_combined(:,:,7:end), ...
    mua_stim_latemove_hit_combined(:,:,5:-1:1), ...
    mua_stim_latemove_miss_combined(:,:,7:end), ...
    mua_stim_latemove_miss_combined(:,:,5:-1:1))),4));
use_data_grid(~logical(use_data_grid)) = NaN;

% Plot depth by condition
plot_depth = 4;
trace_spacing = 1.5;

figure; 
subplot(1,4,1); hold on
p1 = AP_stackplot(squeeze(mua_stim_earlymove_hit_combined(plot_depth,:,:)),t,trace_spacing,false,'k',conditions);
p2 = AP_stackplot(squeeze(mua_stim_earlymove_hit_pred_combined(plot_depth,:,:)),t,trace_spacing,false,'r',conditions);
ylim([trace_spacing,(n_conditions+2)*trace_spacing]);
line([0,0],ylim,'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
xlabel('Time from stim onset');
ylabel(['Condition (MUA depth ' num2str(plot_depth) ')'])
legend([p1(1),p2(1)],{'Real','Predicted'});
title('Early move hit');

subplot(1,4,2); hold on
p1 = AP_stackplot(squeeze(mua_stim_earlymove_miss_combined(plot_depth,:,:)),t,trace_spacing,false,'k',conditions);
p2 = AP_stackplot(squeeze(mua_stim_earlymove_miss_pred_combined(plot_depth,:,:)),t,trace_spacing,false,'r',conditions);
ylim([trace_spacing,(n_conditions+2)*trace_spacing]);
line([0,0],ylim,'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
xlabel('Time from stim onset');
ylabel(['Condition (MUA depth ' num2str(plot_depth) ')'])
legend([p1(1),p2(1)],{'Real','Predicted'});
title('Early move miss');

subplot(1,4,3); hold on
p1 = AP_stackplot(squeeze(mua_stim_latemove_hit_combined(plot_depth,:,:)),t,trace_spacing,false,'k',conditions);
p2 = AP_stackplot(squeeze(mua_stim_latemove_hit_pred_combined(plot_depth,:,:)),t,trace_spacing,false,'r',conditions);
ylim([trace_spacing,(n_conditions+2)*trace_spacing]);
line([0,0],ylim,'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
xlabel('Time from stim onset');
ylabel(['Condition (MUA depth ' num2str(plot_depth) ')'])
legend([p1(1),p2(1)],{'Real','Predicted'});
title('Late move hit');

subplot(1,4,4); hold on
p1 = AP_stackplot(squeeze(mua_stim_latemove_miss_combined(plot_depth,:,:)),t,trace_spacing,false,'k',conditions);
p2 = AP_stackplot(squeeze(mua_stim_latemove_miss_pred_combined(plot_depth,:,:)),t,trace_spacing,false,'r',conditions);
ylim([trace_spacing,(n_conditions+2)*trace_spacing]);
line([0,0],ylim,'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
xlabel('Time from stim onset');
ylabel(['Condition (MUA depth ' num2str(plot_depth) ')'])
legend([p1(1),p2(1)],{'Real','Predicted'});
title('Late move miss');

% Plot all conditions for one depth (early)
figure; 

curr_data = cell2mat(permute(arrayfun(@(x) mat2gray([squeeze(mua_stim_earlymove_hit_combined(x,:,:)); ...
    squeeze(mua_stim_earlymove_hit_pred_combined(x,:,:))]),1:n_depths,'uni',false),[1,3,2]));

col = colormap_BlueWhiteRed((length(conditions)-1)/2);
a1 = subplot(1,3,1); hold on; 
title('Real'); xlabel('Time from stim onset'); ylabel('Str depth');
a2 = subplot(1,3,2); hold on; 
title('Predicted'); xlabel('Time from stim onset'); ylabel('Str depth');
a3 = subplot(1,3,3); hold on; 
title('Difference'); xlabel('Time from stim onset'); ylabel('Str depth');

for curr_cond = 1:n_conditions
    axes(a1);
    AP_stackplot(squeeze(curr_data(1:size(curr_data,1)/2,curr_cond,:)),t,1,false,col(curr_cond,:));
    axes(a2);
    AP_stackplot(squeeze(curr_data(size(curr_data,1)/2+1:end,curr_cond,:)),t,1,false,col(curr_cond,:));
    axes(a3);
    AP_stackplot(squeeze(curr_data(1:size(curr_data,1)/2,curr_cond,:)) - ...
        squeeze(curr_data(size(curr_data,1)/2+1:end,curr_cond,:)),t,1,false,col(curr_cond,:));
end

linkaxes([a1,a2,a3],'xy');

line(a1,[0,0],ylim,'color','k','linestyle','--');
line(a1,[0.5,0.5],ylim,'color','k','linestyle','--');
line(a2,[0,0],ylim,'color','k','linestyle','--');
line(a2,[0.5,0.5],ylim,'color','k','linestyle','--')
line(a3,[0,0],ylim,'color','k','linestyle','--');
line(a3,[0.5,0.5],ylim,'color','k','linestyle','--')


% Plot overlay all conditions for all depth (late)
figure; 

curr_data = cell2mat(permute(arrayfun(@(x) mat2gray([squeeze(mua_stim_latemove_hit_combined(x,:,:)); ...
    squeeze(mua_stim_latemove_hit_pred_combined(x,:,:))]),1:n_depths,'uni',false),[1,3,2]));

col = colormap_BlueWhiteRed((length(conditions)-1)/2);
a1 = subplot(1,3,1); hold on; 
title('Real'); xlabel('Time from stim onset'); ylabel('Str depth');
a2 = subplot(1,3,2); hold on; 
title('Predicted'); xlabel('Time from stim onset'); ylabel('Str depth');
a3 = subplot(1,3,3); hold on; 
title('Difference'); xlabel('Time from stim onset'); ylabel('Str depth');

for curr_cond = 1:n_conditions
    axes(a1);
    AP_stackplot(squeeze(curr_data(1:size(curr_data,1)/2,curr_cond,:)),t,1,false,col(curr_cond,:));
    axes(a2);
    AP_stackplot(squeeze(curr_data(size(curr_data,1)/2+1:end,curr_cond,:)),t,1,false,col(curr_cond,:));
    axes(a3);
    AP_stackplot(squeeze(curr_data(1:size(curr_data,1)/2,curr_cond,:)) - ...
        squeeze(curr_data(size(curr_data,1)/2+1:end,curr_cond,:)),t,1,false,col(curr_cond,:));
end

linkaxes([a1,a2,a3],'xy');

line(a1,[0,0],ylim,'color','k','linestyle','--');
line(a1,[0.5,0.5],ylim,'color','k','linestyle','--');
line(a2,[0,0],ylim,'color','k','linestyle','--');
line(a2,[0.5,0.5],ylim,'color','k','linestyle','--')
line(a3,[0,0],ylim,'color','k','linestyle','--');
line(a3,[0.5,0.5],ylim,'color','k','linestyle','--')

% Plot difference of prediction amplitude for stim/movement
pred_diff_hit = mua_stim_earlymove_hit_combined - mua_stim_earlymove_hit_pred_combined;
pred_diff_miss = mua_stim_earlymove_miss_combined - mua_stim_earlymove_miss_pred_combined;
stim_diff = pred_diff_hit(:,:,8) - pred_diff_miss(:,:,4);
move_diff = pred_diff_hit(:,:,8) - pred_diff_miss(:,:,8);
stim_move_diff = abs(stim_diff) - abs(move_diff);

figure; hold on; trace_spacing = 1;
pos_plot = ones(n_depths,length(t)); pos_plot(stim_move_diff < 0) = NaN;
neg_plot = ones(n_depths,length(t)); neg_plot(stim_move_diff > 0) = NaN;
p1 = AP_stackplot(stim_move_diff'.*pos_plot',t,trace_spacing,false,'b',1:n_depths);
p2 = AP_stackplot(stim_move_diff'.*neg_plot',t,trace_spacing,false,'r',1:n_depths);
xlim([-0.2,01.5])
line([0,0],ylim,'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
xlabel('Time from stim onset');
ylabel('MUA depth');
title('Prediction error'); 

% Plot condition by depth 
plot_condition = 11;
trace_spacing = 4;

figure; 
subplot(1,2,1); hold on
p1 = AP_stackplot(squeeze(mua_stim_earlymove_hit_combined(:,:,plot_condition))',t,trace_spacing,false,'k',1:n_depths);
p2 = AP_stackplot(squeeze(mua_stim_earlymove_hit_pred_combined(:,:,plot_condition))',t,trace_spacing,false,'r',1:n_depths);
ylim([trace_spacing,(n_depths+1)*trace_spacing]);
line([0,0],ylim,'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
xlabel('Time from stim onset');
ylabel(['MUA depth (condition ' num2str(plot_condition) ')'])
legend([p1(1),p2(1)],{'Real','Predicted'});
title('Early move');

subplot(1,2,2); hold on
p1 = AP_stackplot(squeeze(mua_stim_latemove_hit_combined(:,:,plot_condition))',t,trace_spacing,false,'k',1:n_depths);
p2 = AP_stackplot(squeeze(mua_stim_latemove_hit_pred_combined(:,:,plot_condition))',t,trace_spacing,false,'r',1:n_depths);
ylim([trace_spacing,(n_depths+1)*trace_spacing]);
line([0,0],ylim,'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
xlabel('Time from stim onset');
ylabel(['MUA depth (condition ' num2str(plot_condition) ')'])
legend([p1(1),p2(1)],{'Real','Predicted'});
title('Late move');

% Stim vs. movement tuning
stim_earlymove_diff = (mua_stim_earlymove_hit_combined(:,:,7:end)-mua_stim_earlymove_miss_combined(:,:,5:-1:1)).*use_data_grid;
decision_earlymove_diff = (mua_stim_earlymove_hit_combined(:,:,7:end)-mua_stim_earlymove_miss_combined(:,:,7:end)).*use_data_grid;
stim_decision_earlymove_diff = nanmean(abs(stim_earlymove_diff)-abs(decision_earlymove_diff),3);

stim_earlymove_pred_diff = (mua_stim_earlymove_hit_pred_combined(:,:,7:end)-mua_stim_earlymove_miss_pred_combined(:,:,5:-1:1)).*use_data_grid;
decision_earlymove_pred_diff = (mua_stim_earlymove_hit_pred_combined(:,:,7:end)-mua_stim_earlymove_miss_pred_combined(:,:,7:end)).*use_data_grid;
stim_decision_earlymove_pred_diff = nanmean(abs(stim_earlymove_pred_diff)-abs(decision_earlymove_pred_diff),3);

stim_latemove_diff = (mua_stim_latemove_hit_combined(:,:,7:end)-mua_stim_latemove_miss_combined(:,:,5:-1:1)).*use_data_grid;
decision_latemove_diff = (mua_stim_latemove_hit_combined(:,:,7:end)-mua_stim_latemove_miss_combined(:,:,7:end)).*use_data_grid;
stim_decision_latemove_diff = nanmean(abs(stim_latemove_diff)-abs(decision_latemove_diff),3);

stim_latemove_pred_diff = (mua_stim_latemove_hit_pred_combined(:,:,7:end)-mua_stim_latemove_miss_pred_combined(:,:,5:-1:1)).*use_data_grid;
decision_latemove_pred_diff = (mua_stim_latemove_hit_pred_combined(:,:,7:end)-mua_stim_latemove_miss_pred_combined(:,:,7:end)).*use_data_grid;
stim_decision_latemove_pred_diff = nanmean(abs(stim_latemove_pred_diff)-abs(decision_latemove_pred_diff),3);

figure; trace_spacing = 3;
subplot(1,4,1); hold on;
pos_plot = ones(n_depths,length(t)); pos_plot(stim_decision_earlymove_diff < 0) = NaN;
neg_plot = ones(n_depths,length(t)); neg_plot(stim_decision_earlymove_diff > 0) = NaN;
p1 = AP_stackplot(abs(nanmean(stim_earlymove_diff,3))',t,trace_spacing,false,[0.5,0.5,0.5],1:n_depths);
p2 = AP_stackplot(-abs(nanmean(decision_earlymove_diff,3))',t,trace_spacing,false,[0.5,0.5,0.5],1:n_depths);
p3 = AP_stackplot(stim_decision_earlymove_diff'.*pos_plot',t,trace_spacing,false,'b',1:n_depths);
p4 = AP_stackplot(stim_decision_earlymove_diff'.*neg_plot',t,trace_spacing,false,'r',1:n_depths);
ylim([trace_spacing/2,(n_depths+0.5)*trace_spacing]);
line([0,0],ylim,'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
legend([p1(1),p3(1),p4(1)],{'Independent','Stim','Movement'})
xlabel('Time from stim onset (s)');
ylabel('MUA depth');
title('Early move');

subplot(1,4,2); hold on;
pos_plot = ones(n_depths,length(t)); pos_plot(stim_decision_earlymove_pred_diff < 0) = NaN;
neg_plot = ones(n_depths,length(t)); neg_plot(stim_decision_earlymove_pred_diff > 0) = NaN;
p1 = AP_stackplot(abs(nanmean(stim_earlymove_pred_diff,3))',t,trace_spacing,false,[0.5,0.5,0.5],1:n_depths);
p2 = AP_stackplot(-abs(nanmean(decision_earlymove_pred_diff,3))',t,trace_spacing,false,[0.5,0.5,0.5],1:n_depths);
p3 = AP_stackplot(stim_decision_earlymove_pred_diff'.*pos_plot',t,trace_spacing,false,'b',1:n_depths);
p4 = AP_stackplot(stim_decision_earlymove_pred_diff'.*neg_plot',t,trace_spacing,false,'r',1:n_depths);
ylim([trace_spacing/2,(n_depths+0.5)*trace_spacing]);
line([0,0],ylim,'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
legend([p1(1),p3(1),p4(1)],{'Independent','Stim','Movement'})
xlabel('Time from stim onset (s)');
ylabel('MUA depth');
title('Early move predicted');

subplot(1,4,3); hold on;
pos_plot = ones(n_depths,length(t)); pos_plot(stim_decision_latemove_diff < 0) = NaN;
neg_plot = ones(n_depths,length(t)); neg_plot(stim_decision_latemove_diff > 0) = NaN;
p1 = AP_stackplot(abs(nanmean(stim_latemove_diff,3))',t,trace_spacing,false,[0.5,0.5,0.5],1:n_depths);
p2 = AP_stackplot(-abs(nanmean(decision_latemove_diff,3))',t,trace_spacing,false,[0.5,0.5,0.5],1:n_depths);
p3 = AP_stackplot(stim_decision_latemove_diff'.*pos_plot',t,trace_spacing,false,'b',1:n_depths);
p4 = AP_stackplot(stim_decision_latemove_diff'.*neg_plot',t,trace_spacing,false,'r',1:n_depths);
ylim([trace_spacing/2,(n_depths+0.5)*trace_spacing]);
line([0,0],ylim,'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
legend([p1(1),p3(1),p4(1)],{'Independent','Stim','Movement'})
xlabel('Time from stim onset (s)');
ylabel('MUA depth');
title('Late move');

subplot(1,4,4); hold on;
pos_plot = ones(n_depths,length(t)); pos_plot(stim_decision_latemove_pred_diff < 0) = NaN;
neg_plot = ones(n_depths,length(t)); neg_plot(stim_decision_latemove_pred_diff > 0) = NaN;
p1 = AP_stackplot(abs(nanmean(stim_latemove_pred_diff,3))',t,trace_spacing,false,[0.5,0.5,0.5],1:n_depths);
p2 = AP_stackplot(-abs(nanmean(decision_latemove_pred_diff,3))',t,trace_spacing,false,[0.5,0.5,0.5],1:n_depths);
p3 = AP_stackplot(stim_decision_latemove_pred_diff'.*pos_plot',t,trace_spacing,false,'b',1:n_depths);
p4 = AP_stackplot(stim_decision_latemove_pred_diff'.*neg_plot',t,trace_spacing,false,'r',1:n_depths);
ylim([trace_spacing/2,(n_depths+0.5)*trace_spacing]);
line([0,0],ylim,'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
legend([p1(1),p3(1),p4(1)],{'Independent','Stim','Movement'})
xlabel('Time from stim onset (s)');
ylabel('MUA depth');
title('Late move predicted');

%% Load and process MUA-fluor correlations

interval_surround = [-0.5,1.5];
t = linspace(interval_surround(1),interval_surround(2),212);
sample_rate = 1/median(diff(t));

plot_t = [-0.2,0.2];
t_use = t > -0.1 & t < 0;

% Load correlations
corr_use = 'corr_mua_fluor_stim_earlymove_conditionshuff';
fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\' corr_use];
load(fn);
n_depths = 6;

% Load widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
% wf_roi = wf_roi(:,1); % (if only ipsi ROIs)
n_rois = numel(wf_roi);
wf_areas = {wf_roi.area};

corr_mua_mua = nanmean(cell2mat(permute(arrayfun(@(x) ...
    nanmean(cell2mat(batch_vars(x).corr_mua_mua),3),1:length(batch_vars),'uni',false),[1,3,2])),3);
corr_fluor_fluor = nanmean(cell2mat(permute(arrayfun(@(x) ...
    nanmean(cell2mat(batch_vars(x).corr_fluor_fluor),3),1:length(batch_vars),'uni',false),[1,3,2])),3);
corr_fluor_mua = nanmean(cell2mat(permute(arrayfun(@(x) ...
    nanmean(cell2mat(batch_vars(x).corr_fluor_mua),3),1:length(batch_vars),'uni',false),[1,3,2])),3);

corr_mua_wheel = nanmean(cell2mat(permute(arrayfun(@(x) ...
    nanmean(cell2mat(batch_vars(x).corr_mua_wheel),3),1:length(batch_vars),'uni',false),[1,3,2])),3);
corr_fluor_wheel = nanmean(cell2mat(permute(arrayfun(@(x) ...
    nanmean(cell2mat(batch_vars(x).corr_fluor_wheel),3),1:length(batch_vars),'uni',false),[1,3,2])),3);

corr_mua_choice = nanmean(cell2mat(permute(arrayfun(@(x) ...
    nanmean(cell2mat(batch_vars(x).corr_mua_choice),3),1:length(batch_vars),'uni',false),[1,3,2])),3);
corr_fluor_choice = nanmean(cell2mat(permute(arrayfun(@(x) ...
    nanmean(cell2mat(batch_vars(x).corr_fluor_choice),3),1:length(batch_vars),'uni',false),[1,3,2])),3);

corr_mua_mua_split = mat2cell(corr_mua_mua,repmat(length(t),n_depths,1),repmat(length(t),1,n_depths));
corr_fluor_fluor_split = mat2cell(corr_fluor_fluor,repmat(length(t),n_rois,1),repmat(length(t),1,n_rois));
corr_fluor_mua_split = mat2cell(corr_fluor_mua,repmat(length(t),n_rois,1),repmat(length(t),1,n_depths));

corr_mua_wheel_split = mat2cell(corr_mua_wheel,repmat(length(t),n_depths,1),length(t));
corr_fluor_wheel_split = mat2cell(corr_fluor_wheel,repmat(length(t),n_rois,1),length(t));

corr_mua_choice_split = mat2cell(corr_mua_choice,repmat(length(t),n_depths,1),1);
corr_fluor_choice_split = mat2cell(corr_fluor_choice,repmat(length(t),n_rois,1),1);

% Plot grid of correlations

% MUA-MUA
figure; colormap(colormap_BlueWhiteRed);
ax = tight_subplot(n_depths,n_depths,[0.01,0.01]);
for curr_depth1 = 1:n_depths
    for curr_depth2 = 1:n_depths
        axes(ax((curr_depth1-1)*n_depths+curr_depth2));
        imagesc(t,t,corr_mua_mua_split{curr_depth1,curr_depth2})
        caxis([-0.05,0.05])
        title(['Str ' num2str(curr_depth1) '- Str' num2str(curr_depth2)]);
        set(gca,'FontSize',8);
        line([0,0],ylim,'color','k');
        line(xlim,[0,0],'color','k');
        line(xlim,ylim,'color','k');
        axis square off;
        xlim(plot_t)
        ylim(plot_t)
    end
end

% Fluor-Fluor
figure; colormap(colormap_BlueWhiteRed);
ax = tight_subplot(n_rois,n_rois,[0.01,0.01]);
for curr_roi1 = 1:n_rois
    for curr_roi2 = 1:n_rois
        axes(ax((curr_roi1-1)*n_rois+curr_roi2));
        imagesc(t,t,corr_fluor_fluor_split{curr_roi1,curr_roi2})
        caxis([-0.15,0.15])
        title([wf_areas{curr_roi1} '-' wf_areas{curr_roi2}]);
        set(gca,'FontSize',8);
        line([0,0],ylim,'color','k');
        line(xlim,[0,0],'color','k');
        line(xlim,ylim,'color','k');
        axis square off;
        xlim(plot_t)
        ylim(plot_t)
    end
end

% Fluor-MUA
figure; colormap(colormap_BlueWhiteRed);
ax = tight_subplot(n_rois,n_depths,[0.01,0.01]);
for curr_roi = 1:n_rois
    for curr_depth = 1:n_depths
        axes(ax((curr_roi-1)*n_depths+curr_depth));
        imagesc(t,t,corr_fluor_mua_split{curr_roi,curr_depth})
        caxis([-0.05,0.05])
        title([wf_areas{curr_roi} '- Str' num2str(curr_depth)]);
        set(gca,'FontSize',8);
        line([0,0],ylim,'color','k');
        line(xlim,[0,0],'color','k');
        line(xlim,ylim,'color','k');
        axis square off;
        xlim(plot_t)
        ylim(plot_t)
    end
end

% MUA/Fluor-wheel
figure; colormap(colormap_BlueWhiteRed);
ax = tight_subplot(2,max(n_depths,n_rois),[0.01,0.01]);
for curr_depth = 1:n_depths
    axes(ax(curr_depth));
    imagesc(t,t,corr_mua_wheel_split{curr_depth})
    caxis([-0.1,0.1])
    title(['Str' num2str(curr_depth) '-Wheel']);
    set(gca,'FontSize',8);
    line([0,0],ylim,'color','k');
    line(xlim,[0,0],'color','k');
    line(xlim,ylim,'color','k');
    axis square off;
    xlim(plot_t)
    ylim(plot_t)
end
for curr_roi = 1:n_rois
    axes(ax(n_depths+(max(n_rois,n_depths)-n_depths)+curr_roi));
    imagesc(t,t,corr_fluor_wheel_split{curr_roi})
    caxis([-0.2,0.2])
    title([wf_areas{curr_roi} '-Wheel']);
    set(gca,'FontSize',8);
    line([0,0],ylim,'color','k');
    line(xlim,[0,0],'color','k');
    line(xlim,ylim,'color','k');
    axis square off;
    xlim(plot_t)
    ylim(plot_t)
end

% MUA/Fluor-choice
% (to smooth mua)
smooth_size = 5;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

figure; 

subplot(1,2,1); hold on
set(gca,'ColorOrder',copper(n_depths));
plot(t,conv2(horzcat(corr_mua_choice_split{:}),smWin','same'),'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Correlation with decision')
xlabel('Time from movement onset (s)');
title('MUA');
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:6),'uni',false))

subplot(1,2,2); hold on
set(gca,'ColorOrder',[autumn(size(wf_roi,1));winter(size(wf_roi,1))]);
plot(t,conv2(horzcat(corr_fluor_choice_split{:}),smWin','same'),'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Correlation with decision')
xlabel('Time from movement onset (s)')
title('Fluor');
legend(wf_areas);

% Plot pre-movement predicted/predictor
corr_mua_mua_use = cellfun(@(x) x(t_use,t_use),corr_mua_mua_split,'uni',false);
corr_mua_mua_diag = cellfun(@(y) arrayfun(@(x) nanmean(diag(y,x)),-(length(y)-1):length(y)-1),corr_mua_mua_use,'uni',false);

corr_fluor_fluor_use = cellfun(@(x) x(t_use,t_use),corr_fluor_fluor_split,'uni',false);
corr_fluor_fluor_diag = cellfun(@(y) arrayfun(@(x) nanmean(diag(y,x)),-(length(y)-1):length(y)-1),corr_fluor_fluor_use,'uni',false);

corr_fluor_mua_use = cellfun(@(x) x(t_use,t_use),corr_fluor_mua_split,'uni',false);
corr_fluor_mua_diag = cellfun(@(y) arrayfun(@(x) nanmean(diag(y,x)),-(length(y)-1):length(y)-1),corr_fluor_mua_use,'uni',false);

figure;
n_cols = max(n_depths,n_rois);
lags = (-(sum(t_use)-1):(sum(t_use)-1))/sample_rate;
use_lags = lags > 0 & lags < 0.05;

mua_mua_pred = nan(n_depths);
for curr_depth = 1:n_depths    
    subplot(4,n_cols,curr_depth); hold on;
    set(gca,'ColorOrder',copper(n_depths));    
    curr_diag_mean = vertcat(corr_mua_mua_diag{curr_depth,:});
    curr_diag_mean = curr_diag_mean - fliplr(curr_diag_mean);
    curr_diag_mean(curr_depth,:) = NaN;    
    
    max_diag_data = curr_diag_mean(:,use_lags);
    [max_mag,max_mag_idx] = max(abs(max_diag_data),[],2);
    max_sign = arrayfun(@(x) sign(max_diag_data(x,max_mag_idx(x))),1:size(max_diag_data))';    
    mua_mua_pred(:,curr_depth) = max_mag.*max_sign;
    
    plot(lags,curr_diag_mean','linewidth',2)
    line([0,0],ylim,'color','k')
    xlabel('Time (s)')
    ylabel('Predicts - predicted');
    title(['Str ' num2str(curr_depth)]);
    axis tight;    
    xlim([0,max(lags)]);
end

fluor_fluor_pred = nan(n_rois);
for curr_roi = 1:n_rois    
    subplot(4,n_cols,n_cols+curr_roi); hold on;
    set(gca,'ColorOrder',copper(n_rois));    
    curr_diag_mean = vertcat(corr_fluor_fluor_diag{curr_roi,:});
    curr_diag_mean = curr_diag_mean - fliplr(curr_diag_mean);
    curr_diag_mean(curr_roi,:) = NaN;
    
    max_diag_data = curr_diag_mean(:,use_lags);
    [max_mag,max_mag_idx] = max(abs(max_diag_data),[],2);
    max_sign = arrayfun(@(x) sign(max_diag_data(x,max_mag_idx(x))),1:size(max_diag_data))';    
    fluor_fluor_pred(:,curr_roi) = max_mag.*max_sign;
    
    plot(lags,curr_diag_mean','linewidth',2)
    line([0,0],ylim,'color','k')
    xlabel('Time (s)')
    ylabel('Predicts - predicted');
    title(wf_roi(curr_roi).area);
    axis tight;
    xlim([0,max(lags)]);
end

fluor_mua_pred = nan(n_depths,n_rois);
for curr_roi = 1:n_rois    
    subplot(4,n_cols,n_cols*2+curr_roi); hold on;
    set(gca,'ColorOrder',copper(n_depths));    
    curr_diag_mean = vertcat(corr_fluor_mua_diag{curr_roi,:});
    curr_diag_mean = curr_diag_mean - fliplr(curr_diag_mean);
    
    max_diag_data = curr_diag_mean(:,use_lags);
    [max_mag,max_mag_idx] = max(abs(max_diag_data),[],2);
    max_sign = arrayfun(@(x) sign(max_diag_data(x,max_mag_idx(x))),1:size(max_diag_data))';    
    fluor_mua_pred(:,curr_roi) = max_mag.*max_sign;
    
    plot(lags,curr_diag_mean','linewidth',2)
    line([0,0],ylim,'color','k')
    xlabel('Time (s)')
    ylabel('Predicts - predicted');
    title(wf_roi(curr_roi).area);
    axis tight;    
    xlim([0,max(lags)]);
end

for curr_depth = 1:n_depths    
    subplot(4,n_cols,n_cols*3+curr_depth); hold on;
    set(gca,'ColorOrder',copper(n_rois));    
    curr_diag_mean = vertcat(corr_fluor_mua_diag{:,curr_depth});
    curr_diag_mean = curr_diag_mean - fliplr(curr_diag_mean);
    plot(lags,-curr_diag_mean','linewidth',2); % neg b/c orientation
    line([0,0],ylim,'color','k')
    xlabel('Time (s)')
    ylabel('Predicts - predicted');
    title(['Str ' num2str(curr_depth)]);
    axis tight;    
    xlim([0,max(lags)]);
end

% Plot the max/min correlations spatially
figure; 

% Fluor-Fluor
subplot(1,3,1); hold on; axis image off;

ctx_circle_radius = 50;
ctx_circle_radius_text = ctx_circle_radius + 10;
ctx_circle_degs = linspace(0,360,n_rois+1);
ctx_plot_center = ctx_circle_radius*...
    [sind(ctx_circle_degs(1:n_rois))',cosd(ctx_circle_degs(1:n_rois))'];
ctx_text_center = ctx_circle_radius_text*...
    [sind(ctx_circle_degs(1:n_rois))',cosd(ctx_circle_degs(1:n_rois))'];
scatter(ctx_plot_center(:,1),ctx_plot_center(:,2),200, ...
    [autumn(size(wf_roi,1));winter(size(wf_roi,1))],'filled','MarkerEdgeColor','k');
text(ctx_text_center(:,1),ctx_text_center(:,2),{wf_roi.area}, ...
    'HorizontalAlignment','center','VerticalAlignment','middle');

for curr_roi1 = 1:n_rois
    for curr_roi2 = curr_roi1+1:n_rois
        
        curr_weight = 0.01+abs(fluor_fluor_pred(curr_roi2,curr_roi1))*50;
        curr_sign = sign(fluor_fluor_pred(curr_roi2,curr_roi1));
        
        line([ctx_plot_center(curr_roi1,1),ctx_plot_center(curr_roi2,1)], ...
            [ctx_plot_center(curr_roi1,2),ctx_plot_center(curr_roi2,2)], ...
            'linewidth',curr_weight, ...
            'color',[curr_sign == 1,0,curr_sign == -1]);
                
    end    
end

% MUA-MUA
subplot(1,3,2); hold on; axis image off;

str_circle_radius = 50;
str_circle_radius_text = str_circle_radius + 10;
str_circle_degs = linspace(0,360,n_depths+1);
str_plot_center = str_circle_radius*...
    [sind(str_circle_degs(1:n_depths))',cosd(str_circle_degs(1:n_depths))'];
str_text_center = str_circle_radius_text*...
    [sind(str_circle_degs(1:n_depths))',cosd(str_circle_degs(1:n_depths))'];
scatter(str_plot_center(:,1),str_plot_center(:,2),200,copper(n_depths),'filled','MarkerEdgeColor','k');
text(str_text_center(:,1),str_text_center(:,2), ...
    arrayfun(@(x) ['Str ' num2str(x)],1:n_depths,'uni',false), ...
    'HorizontalAlignment','center','VerticalAlignment','middle');

for curr_depth1 = 1:n_depths
    for curr_depth2 = curr_depth1+1:n_depths
        
        curr_weight = 0.01+abs(mua_mua_pred(curr_depth2,curr_depth1))*120;
        curr_sign = sign(mua_mua_pred(curr_depth2,curr_depth1));
        
        line([str_plot_center(curr_depth1,1),str_plot_center(curr_depth2,1)], ...
            [str_plot_center(curr_depth1,2),str_plot_center(curr_depth2,2)], ...
            'linewidth',curr_weight, ...
            'color',[curr_sign == 1,0,curr_sign == -1]);
        
    end
end

% Fluor-MUA
subplot(1,3,3); hold on; axis off;
set(gca,'YDir','reverse');

% (to use wf ROI spatial centers, this gets too messy)

% AP_reference_outline('ccf_aligned',[0.7,0.7,0.7]);AP_reference_outline('retinotopy',[0.7,0.7,1]);
% wf_roi_center = nan(n_rois,2);
% for curr_roi = 1:n_rois
%    curr_roi_boundaries = bwboundaries(wf_roi(curr_roi).mask); 
%    plot(curr_roi_boundaries{1}(:,2),curr_roi_boundaries{1}(:,1),'color',[0,0.7,0],'linewidth',2);    
%    wf_roi_center(curr_roi,:) = [nanmean(curr_roi_boundaries{1}(:,2)),nanmean(curr_roi_boundaries{1}(:,1))];  
% end

% (to just line wf ROIs up)
ctx_plot_center = ...
    [[linspace(-50,-50,size(wf_roi,1))',linspace(0,50,size(wf_roi,1))']; ...
    [linspace(50,50,size(wf_roi,1))',linspace(0,50,size(wf_roi,1))']];
ctx_text_center = ctx_plot_center + ...
    [-15*ones(size(wf_roi,1),1),zeros(size(wf_roi,1),1); ...
    15*ones(size(wf_roi,1),1),zeros(size(wf_roi,1),1)];
scatter(ctx_plot_center(:,1),ctx_plot_center(:,2),200, ...
    [autumn(size(wf_roi,1));winter(size(wf_roi,1))],'filled','MarkerEdgeColor','k');
text(ctx_text_center(:,1),ctx_text_center(:,2),{wf_roi.area}, ...
    'HorizontalAlignment','center','VerticalAlignment','middle');

str_plot_center = [linspace(0,0,n_depths)',linspace(0,50,n_depths)'];
str_text_center = str_plot_center + [zeros(n_depths,1),-2*ones(n_depths,1)];
% str_plot_center = [linspace(215,215,n_depths)',linspace(0,400,n_depths)'];
% rectangle('Position',[min(str_plot_center(:,1))-20,min(str_plot_center(:,2))-20,...
%     range(str_plot_center(:,1))+40,range(str_plot_center(:,2))+40],'FaceColor','w');
scatter(str_plot_center(:,1),str_plot_center(:,2),100,copper(n_depths),'filled','MarkerEdgeColor','k');
text(str_text_center(:,1),str_text_center(:,2), ...
    arrayfun(@(x) ['Str ' num2str(x)],1:n_depths,'uni',false), ...
    'HorizontalAlignment','center','VerticalAlignment','middle');

for curr_roi = 1:n_rois
    for curr_depth = 1:n_depths
        
        curr_weight = 0.01+abs(fluor_mua_pred(curr_depth,curr_roi))*50;
        curr_sign = sign(fluor_mua_pred(curr_depth,curr_roi));
        
        line([ctx_plot_center(curr_roi,1),str_plot_center(curr_depth,1)], ...
            [ctx_plot_center(curr_roi,2),str_plot_center(curr_depth,2)], ...
            'linewidth',curr_weight, ...
            'color',[curr_sign == 1,0,curr_sign == -1]);
        
    end
end
%% Feed-forward/back (first run above)

max_lag = 0.08;
plot_t = [-0.5,1];

max_lag_samples = round(max_lag*sample_rate);
use_t = t >= plot_t(1) & t <= plot_t(2);

local_corr_fluor_fluor = cellfun(@(curr_data) cell2mat(arrayfun(@(x) padarray(diag(curr_data,x),[abs(x),0],NaN,'post'), ...
    -max_lag_samples:max_lag_samples,'uni',false)),corr_fluor_fluor_split,'uni',false);
local_corr_fluor_fluor_diff = cellfun(@(local_corr) nanmean(local_corr(:,max_lag_samples+2:end) - ...
    local_corr(:,max_lag_samples:-1:1),2),local_corr_fluor_fluor,'uni',false);

local_corr_mua_mua = cellfun(@(curr_data) cell2mat(arrayfun(@(x) padarray(diag(curr_data,x),[abs(x),0],NaN,'post'), ...
    -max_lag_samples:max_lag_samples,'uni',false)),corr_mua_mua_split,'uni',false);
local_corr_mua_mua_diff = cellfun(@(local_corr) nanmean(local_corr(:,max_lag_samples+2:end) - ...
    local_corr(:,max_lag_samples:-1:1),2),local_corr_mua_mua,'uni',false);

local_corr_fluor_mua = cellfun(@(curr_data) cell2mat(arrayfun(@(x) padarray(diag(curr_data,x),[abs(x),0],NaN,'post'), ...
    -max_lag_samples:max_lag_samples,'uni',false)),corr_fluor_mua_split,'uni',false);
local_corr_fluor_mua_diff = cellfun(@(local_corr) nanmean(local_corr(:,max_lag_samples+2:end) - ...
    local_corr(:,max_lag_samples:-1:1),2),local_corr_fluor_mua,'uni',false);

% Plot all choice correlations
smooth_size = 5;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

fig_h = figure('Position',[31,109,1838,828]);
subplot(2,1,1);hold on; 

set(gca,'ColorOrder',[copper(n_depths);autumn(size(wf_roi,1));winter(size(wf_roi,1))]);
plot(t,conv2(horzcat(corr_mua_choice_split{:}),smWin','same'),'linewidth',2);
plot(t,conv2(horzcat(corr_fluor_choice_split{:}),smWin','same'),'linewidth',2);

axis tight; xlim(plot_t);
line([0,0],ylim,'color','k');
ylabel('Correlation with decision')
xlabel('Time from stim onset (s)')

t_h = line([0,0],ylim,'color','k');

% Plot the max/min correlations spatially
initial_weight = 0.1;
initial_color = 'k';

fluor_fluor_h = gobjects(n_rois);
mua_mua_h = gobjects(n_depths);
fluor_mua_h = gobjects(n_rois,n_depths);

% Fluor-Fluor
subplot(2,3,4); hold on; axis image off;

ctx_circle_radius = 50;
ctx_circle_radius_text = ctx_circle_radius + 10;
ctx_circle_degs = linspace(0,360,n_rois+1);
ctx_plot_center = ctx_circle_radius*...
    [sind(ctx_circle_degs(1:n_rois))',cosd(ctx_circle_degs(1:n_rois))'];
ctx_text_center = ctx_circle_radius_text*...
    [sind(ctx_circle_degs(1:n_rois))',cosd(ctx_circle_degs(1:n_rois))'];
scatter(ctx_plot_center(:,1),ctx_plot_center(:,2),200, ...
    [autumn(size(wf_roi,1));winter(size(wf_roi,1))],'filled','MarkerEdgeColor','k');
text(ctx_text_center(:,1),ctx_text_center(:,2),{wf_roi.area}, ...
    'HorizontalAlignment','center','VerticalAlignment','middle');

for curr_roi1 = 1:n_rois
    for curr_roi2 = curr_roi1+1:n_rois               
        fluor_fluor_h(curr_roi1,curr_roi2) = ...
            line([ctx_plot_center(curr_roi1,1),ctx_plot_center(curr_roi2,1)], ...
            [ctx_plot_center(curr_roi1,2),ctx_plot_center(curr_roi2,2)], ...
            'linewidth',initial_weight, ...
            'color',initial_color);             
    end    
end

% MUA-MUA
subplot(2,3,5); hold on; axis image off;

str_circle_radius = 50;
str_circle_radius_text = str_circle_radius + 10;
str_circle_degs = linspace(0,360,n_depths+1);
str_plot_center = str_circle_radius*...
    [sind(str_circle_degs(1:n_depths))',cosd(str_circle_degs(1:n_depths))'];
str_text_center = str_circle_radius_text*...
    [sind(str_circle_degs(1:n_depths))',cosd(str_circle_degs(1:n_depths))'];
scatter(str_plot_center(:,1),str_plot_center(:,2),200,copper(n_depths),'filled','MarkerEdgeColor','k');
text(str_text_center(:,1),str_text_center(:,2), ...
    arrayfun(@(x) ['Str ' num2str(x)],1:n_depths,'uni',false), ...
    'HorizontalAlignment','center','VerticalAlignment','middle');

for curr_depth1 = 1:n_depths
    for curr_depth2 = curr_depth1+1:n_depths
        mua_mua_h(curr_depth1,curr_depth2) = ...
            line([str_plot_center(curr_depth1,1),str_plot_center(curr_depth2,1)], ...
            [str_plot_center(curr_depth1,2),str_plot_center(curr_depth2,2)], ...
            'linewidth',initial_weight, ...
            'color',initial_color);        
    end
end

% Fluor-MUA
subplot(2,3,6); hold on; axis off;
set(gca,'YDir','reverse');

ctx_plot_center = ...
    [[linspace(-50,-50,size(wf_roi,1))',linspace(0,50,size(wf_roi,1))']; ...
    [linspace(50,50,size(wf_roi,1))',linspace(0,50,size(wf_roi,1))']];
ctx_text_center = ctx_plot_center + ...
    [-15*ones(size(wf_roi,1),1),zeros(size(wf_roi,1),1); ...
    15*ones(size(wf_roi,1),1),zeros(size(wf_roi,1),1)];
scatter(ctx_plot_center(:,1),ctx_plot_center(:,2),200, ...
    [autumn(size(wf_roi,1));winter(size(wf_roi,1))],'filled','MarkerEdgeColor','k');
text(ctx_text_center(:,1),ctx_text_center(:,2),{wf_roi.area}, ...
    'HorizontalAlignment','center','VerticalAlignment','middle');

str_plot_center = [linspace(0,0,n_depths)',linspace(0,50,n_depths)'];
str_text_center = str_plot_center + [zeros(n_depths,1),-4*ones(n_depths,1)];
scatter(str_plot_center(:,1),str_plot_center(:,2),100,copper(n_depths),'filled','MarkerEdgeColor','k');
text(str_text_center(:,1),str_text_center(:,2), ...
    arrayfun(@(x) ['Str ' num2str(x)],1:n_depths,'uni',false), ...
    'HorizontalAlignment','center','VerticalAlignment','middle');

for curr_roi = 1:n_rois
    for curr_depth = 1:n_depths
        fluor_mua_h(curr_roi,curr_depth) = ...
        line([ctx_plot_center(curr_roi,1),str_plot_center(curr_depth,1)], ...
            [ctx_plot_center(curr_roi,2),str_plot_center(curr_depth,2)], ...
            'linewidth',initial_weight, ...
            'color',initial_color);       
    end
end

% Make movie updating lines
clear movie_frames
curr_frame = 1;
for curr_t = find(use_t)
    
    set(t_h,'XData',repmat(t(curr_t),2,1));
    
    % Fluor-fluor
    for curr_roi1 = 1:n_rois
        for curr_roi2 = curr_roi1+1:n_rois
            curr_corr = local_corr_fluor_fluor_diff{curr_roi1,curr_roi2}(curr_t);
            set(fluor_fluor_h(curr_roi1,curr_roi2), ...
                'linewidth',abs(curr_corr)*50, ...
                'color',[sign(curr_corr) == 1,0,sign(curr_corr) == -1]);  
        end
    end
    % MUA-MUA
    for curr_depth1 = 1:n_depths
        for curr_depth2 = curr_depth1+1:n_depths
            curr_corr = local_corr_mua_mua_diff{curr_depth1,curr_depth2}(curr_t);
            set(mua_mua_h(curr_depth1,curr_depth2), ...
                'linewidth',abs(curr_corr)*120, ...
                'color',[sign(curr_corr) == 1,0,sign(curr_corr) == -1]);  
        end
    end
    % Fluor-MUA
    for curr_roi = 1:n_rois
        for curr_depth = 1:n_depths
            curr_corr = local_corr_fluor_mua_diff{curr_roi,curr_depth}(curr_t);
            set(fluor_mua_h(curr_roi,curr_depth), ...
                'linewidth',abs(curr_corr)*50, ...
                'color',[sign(curr_corr) == 1,0,sign(curr_corr) == -1]);  
        end
    end
    drawnow
    
    warning off;
    movie_frames(curr_frame) = getframe(fig_h);
    warning on;
    curr_frame = curr_frame + 1;
end

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
save_filename = [save_path filesep corr_use '.avi'];

writerObj = VideoWriter(save_filename);
writerObj.FrameRate = sample_rate/5;
open(writerObj);
writeVideo(writerObj,movie_frames);
close(writerObj);

%% Max abs and time of stim/move-aligned correlations

interval_surround = [-0.5,1.5];
t = linspace(interval_surround(1),interval_surround(2),212);
sample_rate = 1/median(diff(t));
corr_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';

corr_fn = {'corr_mua_fluor_stim_earlymove_conditionshuff', ...
    'corr_mua_fluor_move_earlymove_conditionshuff'};
t_use = [{t > 0 & t < 0.2}, ...
    {t > -0.2 & t < 0}];

% Striatum
n_depths = 6;

% Widefield
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);
wf_areas = {wf_roi.area};

str_max = nan(n_depths,2);
str_time = nan(n_depths,2);
wf_max = nan(n_rois,2);
wf_time = nan(n_rois,2);

for curr_align = 1:2
    
    % Load correlations
    load([corr_path filesep corr_fn{curr_align}]);   
    
    % Average correlations across days/animals and split
    corr_mua_choice = nanmean(cell2mat(permute(arrayfun(@(x) ...
        nanmean(cell2mat(batch_vars(x).corr_mua_choice),3),1:length(batch_vars),'uni',false),[1,3,2])),3);
    corr_fluor_choice = nanmean(cell2mat(permute(arrayfun(@(x) ...
        nanmean(cell2mat(batch_vars(x).corr_fluor_choice),3),1:length(batch_vars),'uni',false),[1,3,2])),3);
    
    corr_mua_choice_split = mat2cell(corr_mua_choice,repmat(length(t),n_depths,1),1);
    corr_fluor_choice_split = mat2cell(corr_fluor_choice,repmat(length(t),n_rois,1),1);
    
    % Smooth correlation and get max/time
    smooth_size = 5;
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    
    corr_mua_choice_smooth = conv2(horzcat(corr_mua_choice_split{:}),smWin','same');
    corr_fluor_choice_smooth = conv2(horzcat(corr_fluor_choice_split{:}),smWin','same');
    
    [str_max(:,curr_align),str_time_idx] = max(abs(corr_mua_choice_smooth(t_use{curr_align},:)),[],1);
    [wf_max(:,curr_align),wf_time_idx] = max(abs(corr_fluor_choice_smooth(t_use{curr_align},:)),[],1);
    
    curr_t = t(t_use{curr_align});
    str_time(:,curr_align) = curr_t(str_time_idx);
    wf_time(:,curr_align) = curr_t(wf_time_idx);
    
end

figure; 

% Plot stim vs. move peak
subplot(4,2,1);
scatter(wf_max(:,1),wf_max(:,2),100, ...
    [autumn(size(wf_roi,1));winter(size(wf_roi,1))],'filled','MarkerEdgeColor','k');
xlim([min(wf_max(:)),max(wf_max(:))]);
ylim([min(wf_max(:)),max(wf_max(:))]);
axis square;
line(xlim,ylim,'color','k');
xlabel('Choice corr (stim-aligned)');
ylabel('Choice corr (move-aligned)')
title('Cortex');

subplot(4,2,2);
scatter(str_max(:,1),str_max(:,2),100, ...
    copper(n_depths),'filled','MarkerEdgeColor','k');
xlim([min(str_max(:)),max(str_max(:))]);
ylim([min(str_max(:)),max(str_max(:))]);
axis square;
line(xlim,ylim,'color','k');
xlabel('Choice corr (stim-aligned)');
ylabel('Choice corr (move-aligned)')
title('Striatum');

% Plot times for stim/move separately
subplot(4,2,3);
scatter(1:n_rois,wf_time(:,1),100, ...
    [autumn(size(wf_roi,1));winter(size(wf_roi,1))], ...
    'filled','MarkerEdgeColor','k');
set(gca,'XTick',1:n_rois,'XTickLabel',{wf_roi.area},'XTickLabelRotation',45);
xlabel('ROI')
ylabel('Peak time (stim-aligned)');
subplot(4,2,5);
scatter(1:n_rois,wf_time(:,2),100, ...
    [autumn(size(wf_roi,1));winter(size(wf_roi,1))], ...
    'filled','MarkerEdgeColor','k');
set(gca,'XTick',1:n_rois,'XTickLabel',{wf_roi.area},'XTickLabelRotation',45);
xlabel('ROI')
ylabel('Peak time (move-aligned)');

subplot(4,2,4);
scatter(1:n_depths,str_time(:,1),100, ...
    copper(n_depths),'filled','MarkerEdgeColor','k');
xlabel('Depth')
ylabel('Peak time (stim-aligned)');
subplot(4,2,6);
scatter(1:n_depths,str_time(:,2),100, ...
    copper(n_depths),'filled','MarkerEdgeColor','k');
xlabel('Depth')
ylabel('Peak time (move-aligned)');

% Plot some kind of weighted average for timing
est_move_time = 0.2;
wf_time_scaled = sum(bsxfun(@plus,wf_time,[zeros(n_rois,1),est_move_time*ones(n_rois,1)]) ...
    .*bsxfun(@rdivide,wf_max,sum(wf_max,2)),2);
str_time_scaled = sum(bsxfun(@plus,str_time,[zeros(n_depths,1),est_move_time*ones(n_depths,1)]) ...
    .*bsxfun(@rdivide,str_max,sum(str_max,2)),2);

wf_scale = mat2gray(max(wf_max,[],2))*200+20;
str_scale = mat2gray(max(str_max,[],2))*200+20;

time_scaled_cat = [wf_time_scaled;str_time_scaled];
scale_cat = [wf_scale;str_scale];
label_cat = [{wf_roi.area},arrayfun(@(x) ['Str ' num2str(x)],1:n_depths,'uni',false)];
cmap_cat = [autumn(size(wf_roi,1));winter(size(wf_roi,1));copper(6)];

[~,time_scaled_idx] = sort(time_scaled_cat);

subplot(4,1,4); 
scatter(1:n_rois+n_depths,time_scaled_cat(time_scaled_idx), ...
    scale_cat(time_scaled_idx), ...
    cmap_cat(time_scaled_idx,:), ...
    'filled','MarkerEdgeColor','k');
set(gca,'XTick',1:n_rois+n_depths,'XTickLabel', ...
    label_cat(time_scaled_idx),'XTickLabelRotation',45);
xlabel('Area')
ylabel('Weighted time');

























