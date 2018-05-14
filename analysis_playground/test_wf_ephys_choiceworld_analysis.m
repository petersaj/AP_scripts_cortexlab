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

n_depths = 6;
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
depth_group = discretize(spikeDepths,depth_group_edges);

use_spikes_idx = depth_group == 2;

% use_spikes_idx = ismember(spike_templates,find(templateDepths >= 500 & templateDepths <= 3500));
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

%% !!!!!!!! NEW INDIVIDUAL DAY ANALYSIS (some goes into batch) !!!!!!!!!

%% PSTH viewer

% Spikes in striatum
% use_spikes_idx = spikeDepths >= str_depth(1) & spikeDepths <= str_depth(2);
% use_spikes_idx = aligned_str_depth_group == 4;
use_spikes_idx = spike_templates == 30;

use_spikes = spike_times_timeline(use_spikes_idx);
use_templates = spike_templates(use_spikes_idx);

% Trial properties 
n_trials = length(block.paramsValues);
trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
stim_to_move = padarray(wheel_move_time - stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
stim_to_feedback = padarray(signals_events.responseTimes, ...
    [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
    padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');

go_left = (signals_events.trialSideValues == 1 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == -1 & signals_events.missValues == 1);
go_right = (signals_events.trialSideValues == -1 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == 1 & signals_events.missValues == 1);
trial_choice = go_right - go_left;

% Trials and groupings to use
use_trials = ...
    trial_outcome ~= 0 & ...
    ~signals_events.repeatTrialValues & ...
    stim_to_feedback < 1.5;

conditions = combvec([-1,1],[-1,1])';
trial_conditions = ...
    [signals_events.trialSideValues; trial_choice]';
[~,trial_id] = ismember(trial_conditions,conditions,'rows');

condition_counts = histcounts(trial_id(use_trials), ...
    'BinLimits',[1,size(conditions,1)],'BinMethod','integers')';

raster_window = [-0.5,1];

% Population PSTH
psthViewer(use_spikes,ones(size(use_spikes)), ...
    wheel_move_time(use_trials)',raster_window,trial_id(use_trials));
set(gcf,'Name','Population');

% Template PSTH
psthViewer(use_spikes,use_templates, ...
    wheel_move_time(use_trials)',raster_window,trial_id(use_trials));
set(gcf,'Name','Templates');

% PSTH viewier sorted by choice
use_templates_unique = unique(use_templates);

use_t = [-0.2,-0.02];
use_align = wheel_move_time(use_trials)';

t_peri_event = bsxfun(@plus,use_align,use_t);
d_prime_diff = nan(size(use_templates_unique));
for curr_template_idx = 1:length(use_templates_unique)
    
    curr_template = use_templates_unique(curr_template_idx);
    curr_spikes = spike_times_timeline(spike_templates == curr_template);
    
    curr_spikes_binned = cell2mat(arrayfun(@(x) ...
        histcounts(curr_spikes,t_peri_event(x,:)), ...
        [1:size(t_peri_event,1)]','uni',false))./diff(use_t);
    
    stim_r_trials = signals_events.trialSideValues(use_trials) == 1;
    go_l_trials = trial_choice(use_trials) == -1;
    
    d_prime_choice = ...
        abs(mean(curr_spikes_binned(go_l_trials)) - mean(curr_spikes_binned(~go_l_trials)))./ ...
        sqrt(std(curr_spikes_binned(go_l_trials)).^2 + std(curr_spikes_binned(~go_l_trials)).^2);
    
    d_prime_stim = ...
        abs(mean(curr_spikes_binned(stim_r_trials & go_l_trials)) - mean(curr_spikes_binned(stim_r_trials & ~go_l_trials)))./ ...
        sqrt(std(curr_spikes_binned(stim_r_trials & go_l_trials)).^2 + std(curr_spikes_binned(stim_r_trials & ~go_l_trials)).^2);
    
    d_prime_diff(curr_template_idx) = d_prime_choice-d_prime_stim;
    
end

d_prime_diff(isnan(d_prime_diff)) = -inf;
[~,sort_idx] = sort(d_prime_diff,'descend');

sort_templates = nan(size(templates,1),1);
sort_templates(use_templates_unique(sort_idx)) = use_templates_unique;
template_sort = sort_templates(use_templates);

raster_window = [-0.5,0.5];
psthViewer(use_spikes,template_sort, ...
    wheel_move_time(use_trials)',raster_window,trial_id(use_trials));


%% Ephys: rasters by template

% Pull out spikes within striatum
use_spikes_idx = spikeDepths >= str_depth(1) & spikeDepths <= str_depth(2);

use_spikes = spike_times_timeline(use_spikes_idx);
use_templates = spike_templates(use_spikes_idx);

use_templates_unique = unique(use_templates);

% Define trials to use
n_trials = length(block.paramsValues);
trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
stim_to_move = padarray(wheel_move_time - stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
stim_to_feedback = padarray(signals_events.responseTimes, ...
    [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
    padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');

go_left = (signals_events.trialSideValues == 1 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == -1 & signals_events.missValues == 1);
go_right = (signals_events.trialSideValues == -1 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == 1 & signals_events.missValues == 1);
trial_choice = go_right - go_left;

trial_timing = 1 + (stim_to_move > 0.5);

use_trials = ...
    trial_outcome ~= 0 & ...
    ~signals_events.repeatTrialValues(1:n_trials) & ...
    stim_to_feedback < 1.5;

% Get trial conditions
% [contrast,side,choice,timing]
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];

conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);

trial_conditions = ...
    [signals_events.trialContrastValues; signals_events.trialSideValues; ...
    trial_choice; trial_timing]';
[~,trial_id] = ismember(trial_conditions,conditions,'rows');

condition_counts = histcounts(trial_id(use_trials), ...
    'BinLimits',[1,n_conditions],'BinMethod','integers')';

% Get raster for all chosen templates and all conditions

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t = raster_window(1):psth_bin_size:raster_window(2);
t_bins = conv2(t,[1,1]/2,'valid');

% PSTH for all conditions
template_psth = nan(n_conditions,length(t)-1,length(use_templates_unique),2);
for curr_align = 1:2
    switch curr_align
        case 1
            use_align = stimOn_times;
        case 2
            use_align = wheel_move_time';
            use_align(isnan(use_align)) = 0;
    end
    t_peri_event = bsxfun(@plus,use_align,t);
    for curr_template_idx = 1:length(use_templates_unique)
        
        curr_template = use_templates_unique(curr_template_idx);
        curr_spikes = spike_times_timeline(spike_templates == curr_template);
        
        curr_spikes_binned = cell2mat(arrayfun(@(x) ...
            histcounts(curr_spikes,t_peri_event(x,:)), ...
            [1:size(t_peri_event,1)]','uni',false))./psth_bin_size;
        
        [curr_ids_str,curr_mean_psth] = ...
            grpstats(curr_spikes_binned(use_trials,:), ...
            trial_id(use_trials),{'gname','mean'});
        curr_ids = cellfun(@str2num,curr_ids_str);
        
        template_psth(curr_ids,:,curr_template_idx,curr_align) = curr_mean_psth;
    end
end

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
template_psth_smooth = convn(template_psth,smWin,'same');
AP_image_scroll(template_psth_smooth)
line(repmat(find(t_bins > 0,1),2,1),ylim,'color','r');

% Choice difference psth
% [contrast,side,choice,timing]
left_early = ismember(conditions(:,3:4),[-1,1],'rows');
right_early = ismember(conditions(:,3:4),[1,1],'rows');

move_diff = permute(nanmean(template_psth_smooth(right_early,:,:,:),1) - ...
    nanmean(template_psth_smooth(left_early,:,:,:),1),[3,2,4,1]);

[~,sort_idx] = sort(templateDepths(use_templates_unique));

AP_image_scroll(move_diff(sort_idx,:,:))
caxis([-abs(max(move_diff(:))),abs(max(move_diff(:)))]);
colormap(colormap_BlueWhiteRed);
line(repmat(find(t_bins > 0,1),2,1),ylim,'color','k');


%% Each depth: get left/right PCA groups (after above)
% the point of this was to get 2 basic groups of cells at each depth to
% get at more population-dynamics stuff

% Group striatum depths
n_depths = 6;
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
depth_group = discretize(templateDepths(use_templates_unique),depth_group_edges);

go_left = ismember(conditions(:,2:3),[1,-1],'rows');
go_right = ismember(conditions(:,2:3),[1,1],'rows');

curr_depth = 4;
curr_align = 2;
pca_data = zscore([squeeze(nanmean(template_psth_smooth(go_left,:,depth_group == curr_depth,curr_align),1)); ...
    squeeze(nanmean(template_psth_smooth(go_right,:,depth_group == curr_depth,curr_align),1))],[],1);

[coeff,score,latent] = pca(pca_data);



%% Ephys: rasters by depth

% % (to group striatum depths)
% n_depths = 6;
% depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
% depth_group = discretize(spikeDepths,depth_group_edges);

% (to use aligned striatum depths)
n_depths = n_aligned_depths;
depth_group = aligned_str_depth_group;

% % (for manual depth)
% depth_group_edges = [1500,2200];
% n_depths = length(depth_group_edges) - 1;
% [depth_group_n,depth_group] = histc(spikeDepths,depth_group_edges);

% Define trials to use
n_trials = length(block.paramsValues);
trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
stim_to_move = padarray(wheel_move_time - stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
stim_to_feedback = padarray(signals_events.responseTimes, ...
    [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
    padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');

go_left = (signals_events.trialSideValues == 1 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == -1 & signals_events.missValues == 1);
go_right = (signals_events.trialSideValues == -1 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == 1 & signals_events.missValues == 1);
trial_choice = go_right - go_left;

trial_timing = 1 + (stim_to_move > 0.5);

use_trials = ...
    trial_outcome ~= 0 & ...
    ~signals_events.repeatTrialValues(1:n_trials) & ...
    stim_to_feedback < 1.5;

% Get trial conditions
% [contrast,side,choice,timing]
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];

conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);

trial_conditions = ...
    [signals_events.trialContrastValues; signals_events.trialSideValues; ...
    trial_choice; trial_timing]';
[~,trial_id] = ismember(trial_conditions,conditions,'rows');

condition_counts = histcounts(trial_id(use_trials), ...
    'BinLimits',[1,n_conditions],'BinMethod','integers')';

% Get raster for all chosen templates and all conditions

% Set times for PSTH
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t = raster_window(1):psth_bin_size:raster_window(2);
t_bins = conv2(t,[1,1]/2,'valid');

% PSTH for all conditions
depth_psth = nan(n_conditions,length(t)-1,n_depths,2);
for curr_align = 1:2
    switch curr_align
        case 1
            use_align = stimOn_times;
        case 2
            use_align = wheel_move_time';
            use_align(isnan(use_align)) = 0;
    end
    t_peri_event = bsxfun(@plus,use_align,t);
    for curr_depth = 1:n_depths
        
        curr_spikes = spike_times_timeline(depth_group == curr_depth);
        
        curr_spikes_binned = cell2mat(arrayfun(@(x) ...
            histcounts(curr_spikes,t_peri_event(x,:)), ...
            [1:size(t_peri_event,1)]','uni',false))./psth_bin_size;
        
        [curr_ids_str,curr_mean_psth] = ...
            grpstats(curr_spikes_binned(use_trials,:), ...
            trial_id(use_trials),{'gname','mean'});
        curr_ids = cellfun(@str2num,curr_ids_str);
        
        depth_psth(curr_ids,:,curr_depth,curr_align) = curr_mean_psth;
    end
end

smooth_size = 100;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
depth_psth_smooth = convn(depth_psth,smWin,'same');
AP_image_scroll(depth_psth_smooth)
line(repmat(find(t_bins > 0,1),2,1),ylim,'color','r');

% Choice psth
% [contrast,side,choice,timing]
left_early = ismember(conditions(:,3:4),[-1,1],'rows');
right_early = ismember(conditions(:,3:4),[1,1],'rows');

left_early_psth = squeeze(nansum(bsxfun(@times,depth_psth_smooth(left_early,:,:,:), ...
    condition_counts(left_early)),1)./sum(condition_counts(left_early)));
right_early_psth = squeeze(nansum(bsxfun(@times,depth_psth_smooth(right_early,:,:,:), ...
    condition_counts(right_early)),1)./sum(condition_counts(right_early)));



%% Widefield: rasters by ROI

% Get traces for all pre-drawn ROIs
% Load widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

aUdf = single(AP_align_widefield(animal,day,Udf));
roi_trace = AP_svd_roi(aUdf,fVdf,[],[],cat(3,wf_roi.mask));

% (get ddf)
ddf = diff(roi_trace,[],2);
ddf(ddf < 0) = 0;

% (make L-R traces)
ddf(size(wf_roi,1)+1:end,:) = ...
    ddf(1:size(wf_roi,1),:) - ddf(size(wf_roi,1)+1:end,:);

% Define trials to use
n_trials = length(block.paramsValues);
trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
stim_to_move = padarray(wheel_move_time - stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
stim_to_feedback = padarray(signals_events.responseTimes, ...
    [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
    padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');

go_left = (signals_events.trialSideValues == 1 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == -1 & signals_events.missValues == 1);
go_right = (signals_events.trialSideValues == -1 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == 1 & signals_events.missValues == 1);
trial_choice = go_right - go_left;

trial_timing = 1 + (stim_to_move > 0.5);

use_trials = ...
    trial_outcome ~= 0 & ...
    ~signals_events.repeatTrialValues(1:n_trials) & ...
    stim_to_feedback < 1.5;

% Get trial conditions
% [contrast,side,choice,timing]
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];

conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);

trial_conditions = ...
    [signals_events.trialContrastValues; signals_events.trialSideValues; ...
    trial_choice; trial_timing]';
[~,trial_id] = ismember(trial_conditions,conditions,'rows');

condition_counts = histcounts(trial_id(use_trials), ...
    'BinLimits',[1,n_conditions],'BinMethod','integers')';

% Get event-aligned fluorescence
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):sample_rate:raster_window(2);

roi_psth = nan(n_conditions,length(t),n_rois,2);
for curr_align = 1:2
    switch curr_align
        case 1
            use_align = stimOn_times;
        case 2
            use_align = wheel_move_time';
            use_align(isnan(use_align)) = 0;
    end
    
    t_peri_event = bsxfun(@plus,use_align,t);
    event_aligned_ddf = interp1(conv2(frame_t,[1,1]/2,'valid'),ddf',t_peri_event);
    
    for curr_roi = 1:n_rois      
        
        [curr_ids,curr_mean_psth] = ...
            grpstats(event_aligned_ddf(use_trials,:,curr_roi), ...
            trial_id(use_trials),{'gname',@(x) mean(x,1)});
        curr_ids = cellfun(@str2num,curr_ids);
        
        roi_psth(curr_ids,:,curr_roi,curr_align) = curr_mean_psth;
        
    end
end

AP_image_scroll(roi_psth,{wf_roi.area})
line(repmat(find(t > 0,1),2,1),ylim,'color','r');
colormap(colormap_BlueWhiteRed);
caxis([-max(abs(caxis)),max(abs(caxis))]);

% Choice difference psth
% [contrast,side,choice,timing]
left_early = ismember(conditions(:,3:4),[-1,1],'rows');
right_early = ismember(conditions(:,3:4),[1,1],'rows');

move_diff = permute(squeeze(nanmean(roi_psth(right_early,:,:,:),1) - ...
    nanmean(roi_psth(left_early,:,:,:),1)),[2,1,3]);

AP_image_scroll(move_diff,{'Stim-aligned','Move-aligned'})
line(repmat(find(t > 0,1),2,1),ylim,'color','r');
axis on;
colormap(colormap_BlueWhiteRed);
caxis([-max(abs(caxis)),max(abs(caxis))]);
set(gca,'YTick',1:n_rois,'YTickLabel',{wf_roi.area})

% Plot difference across trials (stim/move-aligned)
early_condition = find(ismember(conditions(:,4),[1],'rows'));
curr_trials = ismember(trial_id,early_condition) & use_trials';

figure;
subplot(2,1,1);
use_align = stimOn_times;
use_align(isnan(use_align)) = 0;
t_peri_event = bsxfun(@plus,use_align,t);
event_aligned_ddf = interp1(conv2(frame_t,[1,1]/2,'valid'),ddf',t_peri_event);

use_t = t > 0.05 & t < 0.15;
ddf_mean = squeeze(nanmean(event_aligned_ddf(:,use_t,:),2));
ddf_group = reshape(bsxfun(@plus,(trial_choice(curr_trials)'+1)/4,1:n_rois),[],1);

plotSpread(reshape(ddf_mean(curr_trials,:),[],1),'distributionIdx', ...
    ddf_group,'distributionColors',repmat({'r','b'},1,n_rois));
set(gca,'XTick',(1:n_rois)+0.25,'XTickLabel',{wf_roi.area});
axis tight
line(xlim,[0,0],'color','k')
ylabel('\Delta\DeltaF/F');
legend({'Go left','Go right'});
title('Stim-aligned')

subplot(2,1,2);
use_align = wheel_move_time';
use_align(isnan(use_align)) = 0;
t_peri_event = bsxfun(@plus,use_align,t);
event_aligned_ddf = interp1(conv2(frame_t,[1,1]/2,'valid'),ddf',t_peri_event);

use_t = t > -0.2 & t < 0;
ddf_mean = squeeze(nanmean(event_aligned_ddf(:,use_t,:),2));
ddf_group = reshape(bsxfun(@plus,(trial_choice(curr_trials)'+1)/4,1:n_rois),[],1);

plotSpread(reshape(ddf_mean(curr_trials,:),[],1),'distributionIdx', ...
    ddf_group,'distributionColors',repmat({'r','b'},1,n_rois));
set(gca,'XTick',(1:n_rois)+0.25,'XTickLabel',{wf_roi.area});
axis tight
line(xlim,[0,0],'color','k')
ylabel('\Delta\DeltaF/F');
legend({'Go left','Go right'});
title('Move-aligned')

% % Stats
% 
% %%%% TESTING ANOVAN via all comparisons? 
% % this returns everything as rank deficient, nothing is modelable
% 
% use_align = wheel_move_time';
% use_t = t > -0.2 & t < 0;
% 
% use_align(isnan(use_align)) = 0;
% t_peri_event = bsxfun(@plus,use_align,t);
% event_aligned_ddf = interp1(conv2(frame_t,[1,1]/2,'valid'),ddf',t_peri_event);
% ddf_mean = squeeze(nanmean(event_aligned_ddf(:,use_t,:),2));
% 
% model = [1,1,0,0,0;1,0,1,0,0;1,0,0,1,0; ...
%     1,1,1,0,0;1,0,1,1,0;1,1,1,1,0];
% [p,tbl,stats] = anovan(reshape(ddf_mean(use_trials,:),[],1), ...
%     [reshape(repmat(1:14,sum(use_trials),1),[],1), ...
%     repmat(trial_conditions(use_trials,:),n_rois,1)], ...
%     model,1,{'ROI';'contrast';'side';'choice';'timing'});
% [results,~,h,gnames] = multcompare(stats,'Dimension',[1,4],'Display','off');
% comp_idx = reshape(1:length(gnames),[],2);
% sig_rois = arrayfun(@(x) results(results(:,1) == comp_idx(x,1) & ...
%     results(:,2) == comp_idx(x,2),6),1:size(comp_idx,1)) < 0.5;
% rois_sig = sum(bsxfun(@times,cat(3,wf_roi.mask),permute((sig_rois-0.5)*2,[1,3,2])),3);
% figure;imagesc(rois_sig);
% caxis([-1,1])
% colormap(colormap_BlueWhiteRed);
% AP_reference_outline('ccf_aligned','k');AP_reference_outline('retinotopy','m');
% title('ROIs significantly modulated by ?');
% axis image off;


% Plot activity by contrast/side for one ROI (for data club)
use_align = stimOn_times;
use_align(isnan(use_align)) = 0;
t_peri_event = bsxfun(@plus,use_align,t);
event_aligned_ddf = interp1(conv2(frame_t,[1,1]/2,'valid'),ddf',t_peri_event);

use_t = t > 0.05 & t < 0.15;
ddf_mean = squeeze(nanmean(event_aligned_ddf(:,use_t,:),2));
ddf_group = reshape(bsxfun(@plus,(trial_choice(curr_trials)'+1)/4,1:n_rois),[],1);

curr_trials = [trial_conditions(:,4) == 1]';

plot_roi = 10;
figure; hold on;
col = [0,0,1;1,0,0];
choice_idx = round((trial_choice+1)/2+1);
choice_col = col(choice_idx,:);
scatter(trial_conditions(curr_trials,1).*trial_conditions(curr_trials,2), ...
    ddf_mean(curr_trials,plot_roi),50,choice_col(curr_trials,:),'filled');

n_shuff = 1000;
warning off;
use_comb = combvec(contrasts,sides)';
data_idx = reshape(1:n_trials*length(t), ...
    n_trials,length(t));
shuff_idx = nan(n_trials,length(t),n_shuff);
for curr_condition = 1:size(use_comb,1)
    curr_shuff_trials = ismember(trial_conditions(:,1:2),use_comb(curr_condition,:),'rows');
    shuff_idx(curr_shuff_trials,:,:) = ...
        shake(repmat(data_idx(curr_shuff_trials,:,:),1,1,n_shuff),1);
end
warning on;

real_diff = mean(ddf_mean(curr_trials & trial_choice == 1,plot_roi)) - mean(ddf_mean(curr_trials & trial_choice == -1,plot_roi));
shuff_diff = nan(n_shuff,1);
warning off;
for curr_shuff = 1:n_shuff
    curr_shuff_choice = trial_choice(:,shuff_idx(:,1,curr_shuff));
    shuff_diff(curr_shuff) = ...
        mean(ddf_mean(curr_trials & curr_shuff_choice == 1,plot_roi)) - ...
        mean(ddf_mean(curr_trials & curr_shuff_choice == -1,plot_roi));
end
warning on;

corr_rank = tiedrank([real_diff;shuff_diff]);
corr_p = corr_rank(1,:)/(n_shuff+1);

figure;hist(shuff_diff,100);
line([real_diff,real_diff],ylim,'color','r','linewidth',2);
xlabel('Activity move right - activity move left')
ylabel('Frequency');


%% Get trial activity by contrast and choice

% Prepare fluorescence
% (load widefield ROIs)
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

aUdf = single(AP_align_widefield(animal,day,Udf));
roi_trace = AP_svd_roi(aUdf,fVdf,[],[],cat(3,wf_roi.mask));

% (get ddf)
ddf = diff(roi_trace,[],2);
ddf(ddf < 0) = 0;

% (make L-R traces)
ddf(size(wf_roi,1)+1:end,:) = ...
    ddf(1:size(wf_roi,1),:) - ddf(size(wf_roi,1)+1:end,:);

% Prepare MUA
% (group striatum depths)
n_depths = 6;
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
depth_group = discretize(spikeDepths,depth_group_edges);

% Get event-aligned activity
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):sample_rate:raster_window(2);

event_aligned_ddf = nan(length(stimOn_times),length(t),n_rois,2);
event_aligned_mua = nan(length(stimOn_times),length(t),n_depths,2);
for curr_align = 1:2
    switch curr_align
        case 1
            use_align = stimOn_times;
            use_align(isnan(use_align)) = 0;
        case 2
            use_align = wheel_move_time';
            use_align(isnan(use_align)) = 0;
    end
    
    t_peri_event = bsxfun(@plus,use_align,t);
    
    % Fluorescence
    event_aligned_ddf(:,:,:,curr_align) = ...
        interp1(conv2(frame_t,[1,1]/2,'valid'),ddf',t_peri_event);
    
    % MUA
    % (raster times)
    t_bins = [t_peri_event-sample_rate/2,t_peri_event(:,end)+sample_rate/2];
    for curr_depth = 1:n_depths
        curr_spikes = spike_times_timeline(depth_group == curr_depth);
        event_aligned_mua(:,:,curr_depth,curr_align) = cell2mat(arrayfun(@(x) ...
            histcounts(curr_spikes,t_bins(x,:)), ...
            [1:size(t_peri_event,1)]','uni',false))./sample_rate;
    end
end

% Smooth MUA
smooth_size = 20;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
event_aligned_mua = convn(event_aligned_mua,smWin,'same');

% Pick trials to use
use_trials = ...
    trial_outcome ~= 0 & ...
    ~signals_events.repeatTrialValues(1:n_trials) & ...
    stim_to_feedback < 1.5;

% Pick times to average across
use_t_stim = t > 0.05 & t < 0.15;
use_t_move = t > -0.15 & t < 0.02;

% Get average activity for both alignments
event_aligned_ddf_avg = nan(sum(use_trials),n_rois,2);
event_aligned_mua_avg = nan(sum(use_trials),n_depths,2);

for curr_align = 1:2
    switch curr_align
        case 1
            use_t = use_t_stim;
        case 2
            use_t = use_t_move;
    end
    event_aligned_ddf_avg(:,:,curr_align) = ...
        squeeze(nanmean(event_aligned_ddf(use_trials,use_t,:,curr_align),2));
    event_aligned_mua_avg(:,:,curr_align) = ...
        squeeze(nanmean(event_aligned_mua(use_trials,use_t,:,curr_align),2));
end

% Package all activity by trial ID (must be more elegant way but whatever)
fluor_trial_act = cell(n_rois,n_conditions,2);
for curr_roi = 1:n_rois
    for curr_align = 1:2
        for curr_condition = 1:n_conditions
            fluor_trial_act{curr_roi,curr_condition,curr_align} = ...
                event_aligned_ddf_avg(trial_id(use_trials) == curr_condition,curr_roi,curr_align);
        end
    end
end

mua_trial_act = cell(n_depths,n_conditions,2);
for curr_depth = 1:n_depths
    for curr_align = 1:2
        for curr_condition = 1:n_conditions
            mua_trial_act{curr_depth,curr_condition,curr_align} = ...
                event_aligned_mua_avg(trial_id(use_trials) == curr_condition,curr_depth,curr_align);
        end
    end
end



%% !!!!!!!! BATCH PROCESSED ANALYSIS !!!!!!!!!


%% Load batch widefield passive

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

protocol = 'stimKalatsky';
% protocol = 'AP_choiceWorldStimPassive';

surround_window = [-0.5,5];
framerate = 35.2;
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

% % (used for df/ddf comparison for data club)
% a = [mat2gray(im(:,:,1:end-1,3),prctile(reshape(im(:,:,:,3),[],1),[1,99])), ...
%     mat2gray(ddf_im(:,:,:,3),prctile(reshape(ddf_im(:,:,:,3),[],1),[1,99]))];
% fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\passive\f_v_ddf';
% 
% plot_t = t_surround > -0.2 & t_surround < 3;
% AP_movie2avi(a(:,:,plot_t),35/2,[0,1],p,fn,cellfun(@num2str,num2cell(t_surround(plot_t)),'uni',false));



%% Make batch widefield choiceworld mean (GENERAL)

clear all
animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

for trialtype_align = {'stim','move'};
    trialtype_align = cell2mat(trialtype_align);
    for trialtype_timing = {'earlymove','latemove'};
        trialtype_timing = cell2mat(trialtype_timing);
        for trialtype_success = {'hit','miss'};
            trialtype_success = cell2mat(trialtype_success);
            
            data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
            im_aligned_avg_combined = zeros(437,416,265,11);
            n_conditions = [];
            for curr_animal = 1:length(animals)
                animal = animals{curr_animal};
                fn = [data_path filesep animal '_im_' trialtype_align '_' trialtype_timing '_' trialtype_success];
                load(fn);
                
                im_aligned_avg_combined = nansum(cat(5,im_aligned_avg_combined,im_aligned_avg),5);
                n_conditions = sum(cat(5,n_conditions,any(any(any(im_aligned_avg_combined,1),2),3)),5);
                AP_print_progress_fraction(curr_animal,length(animals));
            end
            im_aligned_avg_combined = bsxfun(@rdivide,im_aligned_avg_combined,n_conditions);
            save_fn = [data_path filesep 'im_' trialtype_align '_' trialtype_timing '_' trialtype_success '_combined.mat'];
            save(save_fn,'im_aligned_avg_combined','-v7.3');
            
            clearvars -except animals protocol trialtype_align ...
                        trialtype_timing trialtype_success curr_animal  ...
                        experiments curr_day animal batch_vars load_parts
            
        end
    end
end

%% Plot widefield ROIs

wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

roi_cat = cat(3,wf_roi.mask);
roi_col = [autumn(size(wf_roi,1));winter(size(wf_roi,1))];

figure; hold on
set(gca,'YDir','reverse');
AP_reference_outline('ccf_aligned','k');AP_reference_outline('retinotopy','m');
for curr_roi = 1:n_rois
    curr_roi_boundary = cell2mat(bwboundaries(roi_cat(:,:,curr_roi)));
    patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),roi_col(curr_roi,:));
end
axis image off;

%% Load and process widefield choiceworld mean (GENERAL - FULL IMAGE)

surround_window = [-0.5,2];
upsample_rate = 3;

framerate = 35.2;
surround_samplerate = 1/(framerate*upsample_rate);
t_surround = surround_window(1):surround_samplerate:surround_window(2);
t_df = conv2(t_surround,[1,1]/2,'valid');
conditions = [-1,-0.5,-0.25,-0.125,-0.06,0,0.06,0.125,0.25,0.5,1];

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';

trialtype_align = 'stim';
trialtype_timing = 'latemove';
trialtype_success = 'hit';

trialtype = [trialtype_align ' ' trialtype_timing ' ' trialtype_success];

wf_fn = [data_path filesep 'im_' trialtype_align '_' trialtype_timing '_' trialtype_success '_combined.mat'];
load(wf_fn);

% Get ddf
% ddf = diff(im_aligned_avg_combined,[],3);
% ddf(ddf < 0) = 0;
ddf = im_aligned_avg_combined(:,:,1:end-1,:);

AP_image_scroll(ddf,t_df);
axis image;
AP_reference_outline('ccf_aligned','r');AP_reference_outline('retinotopy','b');

% Get traces for all pre-drawn ROIs
% Load widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = length(wf_roi);

figure('Name',trialtype);
for curr_roi = 1:n_rois   
    
    curr_mask_l = wf_roi(curr_roi,1).mask;
    curr_mask_r = curr_mask_l-wf_roi(curr_roi,2).mask;
    
    curr_traces_l = squeeze(...
        sum(sum(bsxfun(@times,ddf,curr_mask_l),1),2)./sum(curr_mask_l(:) ~= 0));
    curr_traces_r = squeeze(...
        sum(sum(bsxfun(@times,ddf,curr_mask_r),1),2)./sum(curr_mask_r(:) ~= 0));
    
    subplot(2,n_rois,curr_roi); hold on;
    set(gca,'ColorOrder',colormap_BlueWhiteRed((length(conditions)-1)/2));
    plot(t_df,curr_traces_l','linewidth',2);
    plot(t_df,curr_traces_l(:,conditions == 0),'k','linewidth',1);
    axis tight;
    line([0,0],ylim,'color','k');
    xlabel(['Time from ' trialtype_align])
    ylabel('\Delta\DeltaF/F');
    title(wf_roi(curr_roi,1).area);
    
    subplot(2,n_rois,n_rois+curr_roi); hold on;
    set(gca,'ColorOrder',colormap_BlueWhiteRed((length(conditions)-1)/2));
    plot(t_df,curr_traces_r','linewidth',2);
    plot(t_df,curr_traces_r(:,conditions == 0),'k','linewidth',1);
    axis tight;
    line([0,0],ylim,'color','k');
    xlabel(['Time from ' trialtype_align])
    ylabel('\Delta\DeltaF/F');
    title([wf_roi(curr_roi,1).area '-' wf_roi(curr_roi,2).area]);
    
end

% (this shouldn't be necessary anymore, reflecting wf is the same as
% reflecting the ROIs which is how ROIs are made now)

% % Reflect widefield, get L-R 
% ddf(isnan(ddf)) = 0;
% ddf_diff = ddf - AP_reflect_widefield(ddf);
% AP_image_scroll(ddf_diff,t_df)
% AP_reference_outline('ccf_aligned','k');AP_reference_outline('retinotopy','m');
% axis image;
% colormap(colormap_BlueWhiteRed);
% caxis([-0.003,0.003]);
% 
% figure('Name',trialtype);
% for curr_roi = 1:n_rois   
%     
%     curr_mask_l = wf_roi(curr_roi,1).mask;
%     
%     curr_traces_l = squeeze(...
%         sum(sum(bsxfun(@times,ddf_diff,curr_mask_l),1),2)./sum(curr_mask_l(:) ~= 0));
%     
%     subplot(1,n_rois,curr_roi); hold on;
%     set(gca,'ColorOrder',colormap_BlueWhiteRed((length(conditions)-1)/2));
%     plot(t_df,curr_traces_l','linewidth',2);
%     plot(t_df,curr_traces_l(:,conditions == 0),'k','linewidth',1);
%     axis tight;
%     line([0,0],ylim,'color','k');
%     xlabel(['Time from ' trialtype_align])
%     ylabel('\Delta\DeltaF/F');
%     title(wf_roi(curr_roi,1).area);
%     
% end

%% Load and process widefield choiceworld mean (GENERAL - ROIs)

use_fluor = 1; % 1 = ddf, 2 = df

% Get times and conditions
surround_window = [-0.5,3];
upsample_rate = 5;

framerate = 35.2;
surround_samplerate = 1/(framerate*upsample_rate);
t = surround_window(1):surround_samplerate:surround_window(2);

contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';

% Load data
data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
load([data_path filesep 'roi_choiceworld']);

% Load widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

% Combine data
roi_psth = cellfun(@(x) squeeze(x(:,:,:,:,use_fluor,:)),{batch_vars.roi_psth},'uni',false);
roi_psth_mean = nanmean(cell2mat(permute(cellfun(@(x) nanmean(x,5),roi_psth,'uni',false),[1,3,4,5,2])),5);

AP_image_scroll(roi_psth_mean,{wf_roi.area})
line(repmat(find(t > 0,1),2,1),ylim,'color','r');

% (make L-R from here on)
roi_psth_mean(:,:,size(wf_roi,1)+1:end,:) = roi_psth_mean(:,:,1:size(wf_roi,1),:) - ...
    roi_psth_mean(:,:,size(wf_roi,1)+1:end,:);
AP_image_scroll(roi_psth_mean(:,:,size(wf_roi,1)+1:end,:),cellfun(@(x) ['\Delta' x],{wf_roi.area},'uni',false));
line(repmat(find(t > 0,1),2,1),ylim,'color','k');
colormap(colormap_BlueWhiteRed)
caxis([-max(abs(caxis)),max(abs(caxis))]);

% Plot all select conditions from given ROI
plot_roi = 13;
plot_success = 1;

figure('Name',wf_roi(plot_roi).area); 

plot_conditions = conditions(:,2) == -plot_success*conditions(:,3) & conditions(:,4) == 1;
plot_color = colormap_BlueWhiteRed(6);
plot_color = [[0,0,0];plot_color(5:-1:1,:);[0,0,0];plot_color(end-5:end,:)];

subplot(2,2,1); hold on;
set(gca,'ColorOrder',plot_color);
plot(t,roi_psth_mean(plot_conditions,:,plot_roi,1)','linewidth',2);
line([0,0],ylim,'color','k');
xlabel('Time from stim')
title('Early move')

subplot(2,2,2); hold on;
set(gca,'ColorOrder',plot_color);
plot(t,roi_psth_mean(plot_conditions,:,plot_roi,2)','linewidth',2);
line([0,0],ylim,'color','k');
xlabel('Time from move')
title('Early move')

plot_conditions = conditions(:,2) == -plot_success*conditions(:,3) & conditions(:,4) == 2;

subplot(2,2,3); hold on;
set(gca,'ColorOrder',plot_color);
plot(t,roi_psth_mean(plot_conditions,:,plot_roi,1)','linewidth',2);
line([0,0],ylim,'color','k');
xlabel('Time from stim')
title('Late move')

subplot(2,2,4); hold on;
set(gca,'ColorOrder',plot_color);
plot(t,roi_psth_mean(plot_conditions,:,plot_roi,2)','linewidth',2);
line([0,0],ylim,'color','k');
xlabel('Time from move')
title('Late move')

% Plot move L - move R from given ROI
plot_roi = 13;

figure('Name',wf_roi(plot_roi).area); 

plot_conditions = conditions(:,4) == 1;

subplot(2,2,1); hold on;
set(gca,'ColorOrder',[autumn(6);winter(6)]);
plot(t,diff(reshape(roi_psth_mean(plot_conditions,:,plot_roi,1)', ...
    [],sum(plot_conditions)/2,2),[],3),'linewidth',2);
line([0,0],ylim,'color','k');
xlabel('Time from stim')
title('Early move')

subplot(2,2,2); hold on;
set(gca,'ColorOrder',[autumn(6);winter(6)]);
plot(t,diff(reshape(roi_psth_mean(plot_conditions,:,plot_roi,2)', ...
    [],sum(plot_conditions)/2,2),[],3),'linewidth',2);line([0,0],ylim,'color','k');
xlabel('Time from move')
title('Early move')

plot_conditions = ismember(conditions(:,[4]),[2],'rows');

subplot(2,2,3); hold on;
set(gca,'ColorOrder',[autumn(6);winter(6)]);
plot(t,diff(reshape(roi_psth_mean(plot_conditions,:,plot_roi,1)', ...
    [],sum(plot_conditions)/2,2),[],3),'linewidth',2);line([0,0],ylim,'color','k');
xlabel('Time from stim')
title('Late move')

subplot(2,2,4); hold on;
set(gca,'ColorOrder',[autumn(6);winter(6)]);
plot(t,diff(reshape(roi_psth_mean(plot_conditions,:,plot_roi,2)', ...
    [],sum(plot_conditions)/2,2),[],3),'linewidth',2);line([0,0],ylim,'color','k');
xlabel('Time from move')
title('Late move')

% Plot all ROIs for one condition
plot_contrast = 0;
plot_align = 1;
plot_timing = 1;

figure; 

plot_condition = ismember(conditions,[plot_contrast,-1,-1,plot_timing],'rows');
subplot(2,2,1); hold on;
set(gca,'ColorOrder',autumn(size(wf_roi,1)));
plot(t,squeeze(roi_psth_mean(plot_condition,:,1:size(wf_roi,1),plot_align)),'linewidth',2);
line([0,0],ylim,'color','k');
title(['Contrast ' num2str(plot_contrast) ', Stim L/Move L'])

plot_condition = ismember(conditions,[plot_contrast,1,-1,plot_timing],'rows');
subplot(2,2,2); hold on;
set(gca,'ColorOrder',autumn(size(wf_roi,1)));
plot(t,squeeze(roi_psth_mean(plot_condition,:,1:size(wf_roi,1),plot_align)),'linewidth',2);
line([0,0],ylim,'color','k');
title(['Contrast ' num2str(plot_contrast) ', Stim R/Move L'])

plot_condition = ismember(conditions,[plot_contrast,-1,1,plot_timing],'rows');
subplot(2,2,3); hold on;
set(gca,'ColorOrder',autumn(size(wf_roi,1)));
plot(t,squeeze(roi_psth_mean(plot_condition,:,1:size(wf_roi,1),plot_align)),'linewidth',2);
line([0,0],ylim,'color','k');
title(['Contrast ' num2str(plot_contrast) ', Stim L/Move R'])

plot_condition = ismember(conditions,[plot_contrast,1,1,plot_timing],'rows');
subplot(2,2,4); hold on;
set(gca,'ColorOrder',autumn(size(wf_roi,1)));
plot(t,squeeze(roi_psth_mean(plot_condition,:,1:size(wf_roi,1),plot_align)),'linewidth',2);
line([0,0],ylim,'color','k');
title(['Contrast ' num2str(plot_contrast) ', Stim R/Move R'])

% % Fit contrast response
% 
% % [depth,time,align,timing,param,success]
% activity_model_params = nan(n_rois,length(t),2,2,3,2);
% for curr_success = 1:2
%     switch curr_success
%         case 1
%             use_success = 1;
%         case 2
%             use_success = -1;
%     end
%     
%     for curr_timing = 1:2
%         
%         use_conditions = conditions(:,2) == -use_success*conditions(:,3) & conditions(:,4) == curr_timing;
%         contrast_sides = zeros(size(conditions,1),2);
%         contrast_sides(conditions(:,2) == -1,1) = conditions(conditions(:,2) == -1,1);
%         contrast_sides(conditions(:,2) == 1,2) = conditions(conditions(:,2) == 1,1);
%         use_contrast_sides = contrast_sides(use_conditions,:);
%         
%         %     activity_model = @(P) P(1) + P(2).*use_contrast_sides(:,1).^P(4) + P(3).*use_contrast_sides(:,2).^P(4);
%         activity_model = @(P) P(1) + P(2).*use_contrast_sides(:,1).^1 + P(3).*use_contrast_sides(:,2).^1;
%         
%         for curr_roi = 1:n_rois
%             for curr_align = 1:2
%                 for curr_t = 1:length(t)
%                     
%                     curr_act = roi_psth_mean(use_conditions,curr_t,curr_roi,curr_align);
%                     model_sse = @(P) sum((curr_act-(activity_model(P))).^2);
%                     
%                     start_params = [0,0,0.005];
%                     
%                     options = struct;
%                     options.MaxFunEvals = 1000;
%                     options.Display = 'off';
%                     
%                     params_fit = fminsearch(model_sse,start_params,options);
%                     activity_model_params(curr_roi,curr_t,curr_align,curr_timing,:,curr_success) = ...
%                         params_fit;
%                     
%                 end
%             end
%             disp(curr_roi);
%         end
%     end
% end
% 
% 
% medfilt_t = 1;
% spacing = 0.005;
% plot_param = 1;
% plot_align = 1;
% 
% figure; 
% s1 = subplot(1,3,1);hold on;
% p1 = AP_stackplot(medfilt1(activity_model_params(:,:,plot_align,1,plot_param,1)',medfilt_t),t,spacing,false,'k',{wf_roi.area});
% p2 = AP_stackplot(medfilt1(activity_model_params(:,:,plot_align,1,plot_param,2)',medfilt_t),t,spacing,false,'r',{wf_roi.area});
% line([0,0],ylim,'color','k');
% line([0.2,0.2],ylim,'color','k');
% xlabel('Time from align')
% ylabel(['Param ' num2str(plot_param)])
% title('Early move')
% legend([p1(1),p2(1)],{'Early hit','Early miss'});
% 
% s2 = subplot(1,3,2);hold on;
% p1 = AP_stackplot(medfilt1(activity_model_params(:,:,plot_align,2,plot_param,1)',medfilt_t),t,spacing,false,'k',{wf_roi.area});
% p2 = AP_stackplot(medfilt1(activity_model_params(:,:,plot_align,2,plot_param,2)',medfilt_t),t,spacing,false,'r',{wf_roi.area});
% line([0,0],ylim,'color','k');
% line([0.5,0.5],ylim,'color','k');
% line([0.7,0.7],ylim,'color','k');
% xlabel('Time from align')
% ylabel(['Param ' num2str(plot_param)])
% title('Late move')
% legend([p1(1),p2(1)],{'Late hit','Late miss'});
% 
% s3 = subplot(1,3,3);hold on;
% p1 = AP_stackplot(medfilt1(activity_model_params(:,:,plot_align,1,plot_param,1)',medfilt_t),t,spacing,false,'b',{wf_roi.area});
% p2 = AP_stackplot(medfilt1(activity_model_params(:,:,plot_align,2,plot_param,1)',medfilt_t),t,spacing,false,'m',{wf_roi.area});
% line([0,0],ylim,'color','k');
% line([0.2,0.2],ylim,'color','k');
% line([0.5,0.5],ylim,'color','k');
% line([0.7,0.7],ylim,'color','k');
% xlabel('Time from align')
% ylabel(['Param ' num2str(plot_param)])
% title('Stim-aligned')
% legend([p1(1),p2(1)],{'Early move hit','Late move hit'});
% 
% linkaxes([s1,s2,s3]);


% Fit contrast response plus choice

% [depth,time,align,timing,param,success]
activity_model_params = nan(n_rois,length(t),2,2,4);

for curr_timing = 1:2
    
    use_conditions = conditions(:,4) == curr_timing;
    contrast_sides = zeros(size(conditions,1),2);
    contrast_sides(conditions(:,2) == -1,1) = conditions(conditions(:,2) == -1,1);
    contrast_sides(conditions(:,2) == 1,2) = conditions(conditions(:,2) == 1,1);
    
    use_choices = conditions(use_conditions,3);
    use_contrast_sides = contrast_sides(use_conditions,:);
    
    activity_model = @(P) P(1) + P(2).*use_contrast_sides(:,1).^1 + P(3).*use_contrast_sides(:,2).^1 + P(4).*use_choices;
    
    for curr_roi = 1:n_rois
        for curr_align = 1:2
            for curr_t = 1:length(t)
                
                curr_act = roi_psth_mean(use_conditions,curr_t,curr_roi,curr_align);
                params_fit = [ones(size(curr_act)),use_contrast_sides,use_choices]\curr_act;
                activity_model_params(curr_roi,curr_t,curr_align,curr_timing,:) = ...
                    params_fit;
                
            end
        end
    end
end


medfilt_t = 1;
plot_param = 4;

figure; 

% L ROIs
spacing = max(reshape(activity_model_params(1:size(wf_roi,1),:,:,:,plot_param),[],1));

s1 = subplot(2,2,1);hold on;
p1 = AP_stackplot(medfilt1(activity_model_params(1:size(wf_roi,1),:,1,1,plot_param)',medfilt_t),t,spacing,false,'k',{wf_roi(:,1).area});
p2 = AP_stackplot(medfilt1(activity_model_params(1:size(wf_roi,1),:,1,2,plot_param)',medfilt_t),t,spacing,false,'r',{wf_roi(:,1).area});
line([0,0],ylim,'color','k');
line([0.5,0.5],ylim,'color','k');
xlabel('Time from stim')
ylabel(['Param ' num2str(plot_param)])
legend([p1(1),p2(1)],{'Early move','Late move'});

s2 = subplot(2,2,2);hold on;
p1 = AP_stackplot(medfilt1(activity_model_params(1:size(wf_roi,1),:,2,1,plot_param)',medfilt_t),t,spacing,false,'k',{wf_roi(:,1).area});
p2 = AP_stackplot(medfilt1(activity_model_params(1:size(wf_roi,1),:,2,2,plot_param)',medfilt_t),t,spacing,false,'r',{wf_roi(:,1).area});
line([0,0],ylim,'color','k');
xlabel('Time from move')
ylabel(['Param ' num2str(plot_param)])
legend([p1(1),p2(1)],{'Early move','Late move'});

% dROIs
spacing = max(reshape(activity_model_params(size(wf_roi,1)+1:end,:,:,:,plot_param),[],1));

s3 = subplot(2,2,3);hold on;
p1 = AP_stackplot(medfilt1(activity_model_params(size(wf_roi,1)+1:end,:,1,1,plot_param)',medfilt_t),t,spacing,false,'k',{wf_roi(:,2).area});
p2 = AP_stackplot(medfilt1(activity_model_params(size(wf_roi,1)+1:end,:,1,2,plot_param)',medfilt_t),t,spacing,false,'r',{wf_roi(:,2).area});
line([0,0],ylim,'color','k');
line([0.5,0.5],ylim,'color','k');
xlabel('Time from stim')
ylabel(['Param ' num2str(plot_param)])
legend([p1(1),p2(1)],{'Early move','Late move'});

s4 = subplot(2,2,4);hold on;
p1 = AP_stackplot(medfilt1(activity_model_params(size(wf_roi,1)+1:end,:,2,1,plot_param)',medfilt_t),t,spacing,false,'k',{wf_roi(:,2).area});
p2 = AP_stackplot(medfilt1(activity_model_params(size(wf_roi,1)+1:end,:,2,2,plot_param)',medfilt_t),t,spacing,false,'r',{wf_roi(:,2).area});
line([0,0],ylim,'color','k');
xlabel('Time from move')
ylabel(['Param ' num2str(plot_param)])
legend([p1(1),p2(1)],{'Early move','Late move'});

linkaxes([s1,s2,s3,s4],'x');


%%% Params from trial activity within day

% [depth,time,align,timing,param]
activity_model_params = cellfun(@(x) squeeze(x(:,:,:,:,:,use_fluor,:)),{batch_vars.activity_model_params},'uni',false);
activity_model_params_mean = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmedian(x,6),activity_model_params,'uni',false),[1,3,4,5,6,2])),6);

plot_param = 4;

spacing = max(abs(reshape(activity_model_params_mean(:,:,:,:,plot_param),[],1)))*1.5;

figure; 
s1 = subplot(1,2,1);hold on;
p1 = AP_stackplot(activity_model_params_mean(:,:,1,1,plot_param)',t,spacing,false,'k',{wf_roi.area});
p2 = AP_stackplot(activity_model_params_mean(:,:,1,2,plot_param)',t,spacing,false,'b',{wf_roi.area});
line([0,0],ylim,'color','k');
line([0.5,0.5],ylim,'color','k');
xlabel('Time from stim')
ylabel(['Param ' num2str(plot_param)])
legend([p1(1),p2(1)],{'Early','Late'});

s2 = subplot(1,2,2);hold on;
p1 = AP_stackplot(activity_model_params_mean(:,:,2,1,plot_param)',t,spacing,false,'k',{wf_roi.area});
p2 = AP_stackplot(activity_model_params_mean(:,:,2,2,plot_param)',t,spacing,false,'b',{wf_roi.area});
line([0,0],ylim,'color','k');
xlabel('Time from move')
ylabel(['Param ' num2str(plot_param)])
legend([p1(1),p2(1)],{'Early','Late'});

linkaxes([s1,s2]);


% Plot max b3 v max b4
use_t_stim = t > 0.05 & t < 0.15;
use_t_move = t > -0.15 & t < -0.02;

max_b3_early = max(max(abs(activity_model_params_mean(:,use_t_stim,1,1,3)),[],2), ...
    max(abs(activity_model_params_mean(:,use_t_move,2,1,3)),[],2));
max_b3_late = max(max(abs(activity_model_params_mean(:,use_t_stim,1,2,3)),[],2), ...
    max(abs(activity_model_params_mean(:,use_t_move,2,2,3)),[],2));

max_b4_early = max(max(abs(activity_model_params_mean(:,use_t_stim,1,1,4)),[],2), ...
    max(abs(activity_model_params_mean(:,use_t_move,2,1,4)),[],2));
max_b4_late = max(max(abs(activity_model_params_mean(:,use_t_stim,1,2,4)),[],2), ...
    max(abs(activity_model_params_mean(:,use_t_move,2,2,4)),[],2));

figure; 
plot_col = [autumn(n_rois/2);winter(n_rois/2)];

subplot(4,1,1); hold on;
scatter(1:n_rois,max_b3_early,100,plot_col,'filled','MarkerEdgeColor','k','linewidth',2);
scatter(1:n_rois,max_b3_late,100,plot_col,'filled','MarkerEdgeColor',[0.5,0.5,0.5],'linewidth',2);
set(gca,'XTick',1:n_rois,'XTickLabel',{wf_roi.area});
ylabel('\beta_3');

subplot(4,1,2); hold on;
scatter(1:n_rois,max_b4_early,100,plot_col,'filled','MarkerEdgeColor','k','linewidth',2);
scatter(1:n_rois,max_b4_late,100,plot_col,'filled','MarkerEdgeColor',[0.5,0.5,0.5],'linewidth',2);
set(gca,'XTick',1:n_rois,'XTickLabel',{wf_roi.area});
ylabel('\beta_4');

subplot(4,1,3); hold on;
scatter(1:n_rois,(max_b3_early-max_b4_early)./(max_b3_early+max_b4_early), ...
    100,plot_col,'filled','MarkerEdgeColor','k','linewidth',2);
scatter(1:n_rois,(max_b3_late-max_b4_late)./(max_b3_late+max_b4_late), ...
    100,plot_col,'filled','MarkerEdgeColor',[0.5,0.5,0.5],'linewidth',2);
set(gca,'XTick',1:n_rois,'XTickLabel',{wf_roi.area});
ylabel('(\beta_3-\beta_4)/(\beta_3+\beta_4)');

subplot(4,1,4); hold on;
scatter(max_b3_early,max_b4_early,100,plot_col,'filled','MarkerEdgeColor','k','linewidth',2);
scatter(max_b3_late,max_b4_late,100,plot_col,'filled','MarkerEdgeColor',[0.5,0.5,0.5],'linewidth',2);
xlabel('\beta_3');ylabel('\beta_4')
max_weight = max([max_b3_early;max_b3_late;max_b4_early;max_b4_late]);
line([0,max_weight],[0,max_weight],'color','k');
axis image





%% Load and process striatal MUA during choiceworld (OLD)

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


% Plot overlay all conditions for all depths
figure; 

subplot(1,2,1); hold on;
curr_data = cell2mat(permute(arrayfun(@(x) ...
    mat2gray(squeeze(mua.stim_earlymove_hit(x,:,:))),1:n_depths,'uni',false),[1,3,2]));
col = colormap_BlueWhiteRed((length(conditions)-1)/2);
col(conditions == 0,:) = 0;
for curr_cond = 1:n_conditions
    AP_stackplot(squeeze(curr_data(:,curr_cond,:)),t_bins,1,false,col(curr_cond,:));
end
line([0,0],ylim,'color','k');
xlabel('Time from stim onset');

subplot(1,2,2); hold on;
curr_data = cell2mat(permute(arrayfun(@(x) ...
    mat2gray(squeeze(mua.move_earlymove_hit(x,:,:))),1:n_depths,'uni',false),[1,3,2]));
col = colormap_BlueWhiteRed((length(conditions)-1)/2);
col(conditions == 0,:) = 0;
for curr_cond = 1:n_conditions
    AP_stackplot(squeeze(curr_data(:,curr_cond,:)),t_bins,1,false,col(curr_cond,:));
end
line([0,0],ylim,'color','k');
xlabel('Time from move onset');


%% Load and process striatal MUA during choiceworld (OLD 2)

data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
mua_fn = [data_path filesep 'mua_choiceworld'];
load(mua_fn);

conditions = [-1,-0.5,-0.25,-0.125,-0.06,0,0.06,0.125,0.25,0.5,1];
n_conditions = length(conditions);
n_depths = 6;

raster_window = [-0.5,3];
psth_bin_size = 0.001;
t = raster_window(1):psth_bin_size:raster_window(2);
t_bins = conv2(t,[1,1]/2,'valid');

% Loop through all MUA fields, normalize and combine
smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

t_baseline = t_bins < 0;

softnorm = 1;

% (get baseline from all stim-aligned MUAs)
mua_baseline = cell(size(batch_vars));
for curr_animal = 1:length(batch_vars)
   mua_baseline{curr_animal} = cellfun(@(x) ...
       nanmean(nanmean(x(:,t_baseline,:),2),3), ...
       {batch_vars(curr_animal).mua.stim_earlymove_hit},'uni',false);
end

mua_fieldnames = fieldnames(batch_vars(1).mua);
mua = struct;
for curr_field = mua_fieldnames'
    curr_field = curr_field{:};
    
    curr_mua_daycat = arrayfun(@(x) cat(4,batch_vars(x).mua(:).(curr_field)),1:length(batch_vars),'uni',false);
    curr_mua_daycat_norm = cellfun(@(data,baseline) ...
        bsxfun(@rdivide,data,cat(4,baseline{:})+softnorm), ...
        curr_mua_daycat,mua_baseline,'uni',false);
    
    curr_mua_smoothed = cellfun(@(x) convn(x,smWin,'same'),curr_mua_daycat_norm,'uni',false);    
    curr_mua_mean = cellfun(@(x) nanmean(x,4),curr_mua_smoothed,'uni',false);
    curr_mua_combined = nanmean(cat(4,curr_mua_mean{:}),4);
    
    mua.(curr_field) = curr_mua_combined;    
end


trialtype_align = 'move';
trialtype_timing = 'earlymove';
trialtype_success = 'miss';
trialtype = [trialtype_align '_' trialtype_timing '_' trialtype_success];

figure; hold on;
curr_data = cell2mat(permute(arrayfun(@(x) ...
    mat2gray(squeeze(mua.(trialtype)(x,:,:))),1:n_depths,'uni',false),[1,3,2]));
col = colormap_BlueWhiteRed((length(conditions)-1)/2);
col(conditions == 0,:) = 0;
for curr_cond = 1:n_conditions
    AP_stackplot(squeeze(curr_data(:,curr_cond,:)),t_bins,1,false,col(curr_cond,:));
end
line([0,0],ylim,'color','k');
if strcmp(trialtype_align,'stim');
    line([0.5,0.5],ylim,'color','k');
end
xlabel(['Time from ' trialtype_align]);
title([trialtype_timing ' ' trialtype_success])

%% Load and process striatal MUA during choiceworld (NEW)

data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
mua_fn = [data_path filesep 'mua_choiceworld'];
load(mua_fn);

n_depths = 6;

contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';

raster_window = [-0.5,3];
psth_bin_size = 0.001;
t_edges = raster_window(1):psth_bin_size:raster_window(2);
t = conv2(t_edges,[1,1]/2,'valid');

% Normalize, smooth, and combine PSTH
% (PSTH: [condition,time,depth,align,day])
smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

t_baseline = t < 0;

softnorm = 50;

depth_psth = {batch_vars.depth_psth};
mua_baseline = cellfun(@(x) nanmean(x(:,t_baseline,:,1,:),2),depth_psth,'uni',false);
depth_psth_norm = cellfun(@(mua,baseline) ...
    bsxfun(@rdivide,mua,baseline + softnorm),depth_psth,mua_baseline,'uni',false);
depth_psth_smoothed = cellfun(@(x) convn(x,smWin,'same'),depth_psth_norm,'uni',false);
depth_psth_mean = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmean(x,5),depth_psth_smoothed,'uni',false),[1,3,4,5,2])),5);

AP_image_scroll(depth_psth_mean, ...
    cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false));
line(repmat(find(t > 0,1),2,1),ylim,'color','r');

% Stack success plots 
plot_color = colormap_BlueWhiteRed(6);
plot_color = [plot_color(5:-1:1,:);plot_color(end-5:end,:)];
plot_success = 1;
plot_timing = 2;
plot_conditions = find(conditions(:,1) ~= 0 & conditions(:,2) == -plot_success*conditions(:,3) & conditions(:,4) == plot_timing);
figure; 
p1 = subplot(1,2,1); hold on;
for curr_condition_idx = 1:length(plot_conditions);
    curr_condition = plot_conditions(curr_condition_idx);
    AP_stackplot(squeeze(depth_psth_mean(curr_condition,:,:,1)),t,3,false,plot_color(curr_condition_idx,:));
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');

p2 = subplot(1,2,2); hold on;
for curr_condition_idx = 1:length(plot_conditions);
    curr_condition = plot_conditions(curr_condition_idx);
    AP_stackplot(squeeze(depth_psth_mean(curr_condition,:,:,2)),t,3,false,plot_color(curr_condition_idx,:));
end
line([0,0],ylim,'color','k');
xlabel('Time from move');

linkaxes([p1,p2]);

% Stack L/R plots 
plot_conditions = [find(conditions(:,3) == -1 & conditions(:,4) == 1); ...
    find(conditions(:,3) == 1 & conditions(:,4) == 1)];
plot_color = [repmat([0,0,0],12,1);repmat([0,0,1],12,1)];

figure; 
p1 = subplot(1,2,1); hold on;
for curr_condition_idx = 1:length(plot_conditions);
    curr_condition = plot_conditions(curr_condition_idx);
    AP_stackplot(squeeze(depth_psth_mean(curr_condition,:,:,1)),t,3,false,plot_color(curr_condition_idx,:));
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');

p2 = subplot(1,2,2); hold on;
for curr_condition_idx = 1:length(plot_conditions);
    curr_condition = plot_conditions(curr_condition_idx);
    AP_stackplot(squeeze(depth_psth_mean(curr_condition,:,:,2)),t,3,false,plot_color(curr_condition_idx,:));
end
line([0,0],ylim,'color','k');
xlabel('Time from move');

linkaxes([p1,p2]);

% Plot all hit/miss conditions for one depth
plot_str = 3;
plot_success = 1;
plot_color = colormap_BlueWhiteRed(6);
plot_color = [[0,0,0];plot_color(5:-1:1,:);[0,0,0];plot_color(end-5:end,:)];

figure('Name',['Str ' num2str(plot_str)]); 

subplot(2,2,1); hold on;
set(gca,'ColorOrder',plot_color);
plot_conditions = conditions(:,2) == -plot_success*conditions(:,3) & conditions(:,4) == 1;
plot(t,depth_psth_mean(plot_conditions,:,plot_str,1)','linewidth',2);
line([0,0],ylim,'color','k');
xlabel('Time from stim')
title('Early move')

subplot(2,2,2); hold on;
set(gca,'ColorOrder',plot_color);
plot_conditions = conditions(:,2) == -plot_success*conditions(:,3) & conditions(:,4) == 1;
plot(t,depth_psth_mean(plot_conditions,:,plot_str,2)','linewidth',2);
line([0,0],ylim,'color','k');
xlabel('Time from move')
title('Early move')

subplot(2,2,3); hold on;
set(gca,'ColorOrder',plot_color);
plot_conditions = conditions(:,2) == -plot_success*conditions(:,3) & conditions(:,4) == 2;
plot(t,depth_psth_mean(plot_conditions,:,plot_str,1)','linewidth',2);
line([0,0],ylim,'color','k');
xlabel('Time from stim')
title('Late move')

subplot(2,2,4); hold on;
set(gca,'ColorOrder',plot_color);
plot_conditions = conditions(:,2) == -plot_success*conditions(:,3) & conditions(:,4) == 2;
plot(t,depth_psth_mean(plot_conditions,:,plot_str,2)','linewidth',2);
line([0,0],ylim,'color','k');
xlabel('Time from move')
title('Late move')

% Plot all depths for one condition
plot_contrast = 0.125;
plot_align = 1;
plot_timing = 2;

figure; 

subplot(2,2,1); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot_condition = ismember(conditions,[plot_contrast,-1,-1,plot_timing],'rows');
plot(t,squeeze(depth_psth_mean(plot_condition,:,:,plot_align)),'linewidth',2);
line([0,0],ylim,'color','k');
title(['Contrast ' num2str(plot_contrast) ', Stim L/Move L'])

subplot(2,2,2); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot_condition = ismember(conditions,[plot_contrast,1,-1,plot_timing],'rows');
plot(t,squeeze(depth_psth_mean(plot_condition,:,:,plot_align)),'linewidth',2);
line([0,0],ylim,'color','k');
title(['Contrast ' num2str(plot_contrast) ', Stim R/Move L'])

subplot(2,2,3); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot_condition = ismember(conditions,[plot_contrast,-1,1,plot_timing],'rows');
plot(t,squeeze(depth_psth_mean(plot_condition,:,:,plot_align)),'linewidth',2);
line([0,0],ylim,'color','k');
title(['Contrast ' num2str(plot_contrast) ', Stim L/Move R'])

subplot(2,2,4); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot_condition = ismember(conditions,[plot_contrast,1,1,plot_timing],'rows');
plot(t,squeeze(depth_psth_mean(plot_condition,:,:,plot_align)),'linewidth',2);
line([0,0],ylim,'color','k');
title(['Contrast ' num2str(plot_contrast) ', Stim R/Move R'])

% Plot move L-R difference for all depths
plot_contrast = 0.5;
plot_align = 2;
plot_timing = 1;

figure; 

subplot(1,2,1); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot_L = ismember(conditions,[plot_contrast,-1,-1,plot_timing],'rows');
plot_R = ismember(conditions,[plot_contrast,-1,1,plot_timing],'rows');
plot(t, ...
    squeeze(depth_psth_mean(plot_L,:,:,plot_align)) - ...
    squeeze(depth_psth_mean(plot_R,:,:,plot_align)),'linewidth',2);
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');
ylabel('Move L - Move R');
title(['Contrast ' num2str(plot_contrast) ', Stim L'])

subplot(1,2,2); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot_L = ismember(conditions,[plot_contrast,1,-1,plot_timing],'rows');
plot_R = ismember(conditions,[plot_contrast,1,1,plot_timing],'rows');
plot(t, ...
    squeeze(depth_psth_mean(plot_L,:,:,plot_align)) - ...
    squeeze(depth_psth_mean(plot_R,:,:,plot_align)),'linewidth',2);line([0,0],ylim,'color','k');
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');
ylabel('Move L - Move R');
title(['Contrast ' num2str(plot_contrast) ', Stim R'])

% Plot move L-R difference for one depth/all conditions
plot_str = 2;
plot_align = 2;
plot_timing = 1;

figure; 

subplot(1,2,1); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot_L = ismember(conditions(:,2:4),[-1,-1,plot_timing],'rows');
plot_R = ismember(conditions(:,2:4),[-1,1,plot_timing],'rows');
plot(t, ...
    squeeze(depth_psth_mean(plot_L,:,plot_str,plot_align))' - ...
    squeeze(depth_psth_mean(plot_R,:,plot_str,plot_align))','linewidth',2);
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');
ylabel('Move L - Move R');
title(['Str ' num2str(plot_str), ',Stim L']);

subplot(1,2,2); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot_L = ismember(conditions(:,2:4),[1,-1,plot_timing],'rows');
plot_R = ismember(conditions(:,2:4),[1,1,plot_timing],'rows');
plot(t, ...
    squeeze(depth_psth_mean(plot_L,:,plot_str,plot_align))' - ...
    squeeze(depth_psth_mean(plot_R,:,plot_str,plot_align))','linewidth',2);
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');
ylabel('Move L - Move R');
title(['Str ' num2str(plot_str), ',Stim R']);
legend(cellfun(@num2str,num2cell(contrasts),'uni',false))


% % Fit contrast response
% 
% % [depth,time,align,timing,param,success]
% activity_model_params = nan(n_depths,length(t),2,2,3,2);
% for curr_success = 1:2
%     switch curr_success
%         case 1
%             use_success = 1;
%         case 2
%             use_success = -1;
%     end
%     
%     for curr_timing = 1:2
%         
%         use_conditions = conditions(:,2) == -use_success*conditions(:,3) & conditions(:,4) == curr_timing;
%         contrast_sides = zeros(size(conditions,1),2);
%         contrast_sides(conditions(:,2) == -1,1) = conditions(conditions(:,2) == -1,1);
%         contrast_sides(conditions(:,2) == 1,2) = conditions(conditions(:,2) == 1,1);
%         use_contrast_sides = contrast_sides(use_conditions,:);
%         
%         %     activity_model = @(P) P(1) + P(2).*use_contrast_sides(:,1).^P(4) + P(3).*use_contrast_sides(:,2).^P(4);
%         activity_model = @(P) P(1) + P(2).*use_contrast_sides(:,1).^0.3 + P(3).*use_contrast_sides(:,2).^0.3;
%         
%         for curr_depth = 1:n_depths
%             for curr_align = 1:2
%                 for curr_t = 1:length(t)
%                     
%                     curr_act = depth_psth_mean(use_conditions,curr_t,curr_depth,curr_align);
%                     model_sse = @(P) sum((curr_act-(activity_model(P))).^2);
%                     
%                     start_params = [0,0,1];
%                     
%                     options = struct;
%                     options.MaxFunEvals = 1000;
%                     options.Display = 'off';
%                     
%                     params_fit = fminsearch(model_sse,start_params,options);
%                     activity_model_params(curr_depth,curr_t,curr_align,curr_timing,:,curr_success) = ...
%                         params_fit;
%                     
%                     % Doesn't make sense since not cross-validated
%                     %                 sse_act = sum(curr_act.^2);
%                     %                 sse_residual = sum((curr_act-activity_model(params_fit)).^2);
%                     %                 explained_var(curr_depth,curr_t,curr_align,curr_timing) = ...
%                     %                     (sse_act - sse_residual)/sse_act;
%                     
%                 end
%             end
%             disp(curr_depth);
%         end
%     end
% end
% 
% medfilt_t = 1;
% spacing = 2;
% plot_param = 3;
% plot_align = 1;
% 
% figure; 
% s1 = subplot(1,3,1);hold on;
% p1 = AP_stackplot(medfilt1(activity_model_params(:,:,plot_align,1,plot_param,1)',medfilt_t),t,spacing,false,'k');
% p2 = AP_stackplot(medfilt1(activity_model_params(:,:,plot_align,1,plot_param,2)',medfilt_t),t,spacing,false,'r');
% line([0,0],ylim,'color','k');
% line([0.2,0.2],ylim,'color','k');
% xlabel('Time from align')
% ylabel(['Param ' num2str(plot_param)])
% title('Early move')
% legend([p1(1),p2(1)],{'Early hit','Early miss'});
% 
% s2 = subplot(1,3,2);hold on;
% p1 = AP_stackplot(medfilt1(activity_model_params(:,:,plot_align,2,plot_param,1)',medfilt_t),t,spacing,false,'k');
% p2 = AP_stackplot(medfilt1(activity_model_params(:,:,plot_align,2,plot_param,2)',medfilt_t),t,spacing,false,'r');
% line([0,0],ylim,'color','k');
% line([0.5,0.5],ylim,'color','k');
% line([0.7,0.7],ylim,'color','k');
% xlabel('Time from align')
% ylabel(['Param ' num2str(plot_param)])
% title('Late move')
% legend([p1(1),p2(1)],{'Late hit','Late miss'});
% 
% s3 = subplot(1,3,3);hold on;
% p1 = AP_stackplot(medfilt1(activity_model_params(:,:,plot_align,1,plot_param,1)',medfilt_t),t,spacing,false,'b');
% p2 = AP_stackplot(medfilt1(activity_model_params(:,:,plot_align,2,plot_param,1)',medfilt_t),t,spacing,false,'m');
% line([0,0],ylim,'color','k');
% line([0.2,0.2],ylim,'color','k');
% line([0.5,0.5],ylim,'color','k');
% line([0.7,0.7],ylim,'color','k');
% xlabel('Time from align')
% ylabel(['Param ' num2str(plot_param)])
% title('Stim-aligned')
% legend([p1(1),p2(1)],{'Early move hit','Late move hit'});
% 
% linkaxes([s1,s2,s3]);


% Fit contrast response plus choice
contrast_exp = 0.3;

% [depth,time,align,timing,param,success]
activity_model_params = nan(n_depths,length(t),2,2,4);

for curr_timing = 1:2
    
    use_conditions = conditions(:,4) == curr_timing;
    contrast_sides = zeros(size(conditions,1),2);
    contrast_sides(conditions(:,2) == -1,1) = conditions(conditions(:,2) == -1,1);
    contrast_sides(conditions(:,2) == 1,2) = conditions(conditions(:,2) == 1,1);
    
    use_choices = conditions(use_conditions,3);
    use_contrast_sides = contrast_sides(use_conditions,:);
    
    activity_model = @(P) P(1) + P(2).*use_contrast_sides(:,1).^contrast_exp + ...
        P(3).*use_contrast_sides(:,2).^contrast_exp + P(4).*use_choices;
    
    for curr_depth = 1:n_depths
        for curr_align = 1:2
            for curr_t = 1:length(t)
                
                curr_act = depth_psth_mean(use_conditions,curr_t,curr_depth,curr_align);
                
                % to fit non-linear
                %                 model_sse = @(P) sum((curr_act-(activity_model(P))).^2);
                %
                %                 start_params = [0,0,1,0];
                %
                %                 options = struct;
                %                 options.MaxFunEvals = 1000;
                %                 options.Display = 'off';
                %
                %                 params_fit = fminsearch(model_sse,start_params,options);
                %                 activity_model_params(curr_roi,curr_t,curr_align,curr_timing,:) = ...
                %                     params_fit;
                
                % to fit linear
                params_fit = [ones(size(curr_act)),use_contrast_sides.^contrast_exp,use_choices]\curr_act;
                activity_model_params(curr_depth,curr_t,curr_align,curr_timing,:) = ...
                    params_fit;
                
            end
        end
    end
end



plot_param = 3;

spacing = max(abs(reshape(activity_model_params(:,:,:,:,plot_param),[],1)));

figure; 
s1 = subplot(1,2,1);hold on;
p1 = AP_stackplot(activity_model_params(:,:,1,1,plot_param)',t,spacing,false,'k',1:6);
p2 = AP_stackplot(activity_model_params(:,:,1,2,plot_param)',t,spacing,false,'r',1:6);
line([0,0],ylim,'color','k');
line([0.5,0.5],ylim,'color','k');
xlabel('Time from align')
ylabel(['Param ' num2str(plot_param)])
title('Stim aligned')
legend([p1(1),p2(1)],{'Early move','Late move'});

s2 = subplot(1,2,2);hold on;
p1 = AP_stackplot(activity_model_params(:,:,2,1,plot_param)',t,spacing,false,'k',1:6);
p2 = AP_stackplot(activity_model_params(:,:,2,2,plot_param)',t,spacing,false,'r',1:6);
line([0,0],ylim,'color','k');
xlabel('Time from align')
ylabel(['Param ' num2str(plot_param)])
title('Move aligned')
legend([p1(1),p2(1)],{'Early move','Late move'});

linkaxes([s1,s2]);


% Plot max b3 v max b4
use_t_stim = t > 0.05 & t < 0.15;
use_t_move = t > -0.15 & t < -0.02;

max_b3_early = max(max(abs(activity_model_params(:,use_t_stim,1,1,3)),[],2), ...
    max(abs(activity_model_params(:,use_t_move,2,1,3)),[],2));
max_b3_late = max(max(abs(activity_model_params(:,use_t_stim,1,2,3)),[],2), ...
    max(abs(activity_model_params(:,use_t_move,2,2,3)),[],2));

max_b4_early = max(max(abs(activity_model_params(:,use_t_stim,1,1,4)),[],2), ...
    max(abs(activity_model_params(:,use_t_move,2,1,4)),[],2));
max_b4_late = max(max(abs(activity_model_params(:,use_t_stim,1,2,4)),[],2), ...
    max(abs(activity_model_params(:,use_t_move,2,2,4)),[],2));

figure; 
plot_col = copper(n_depths);

subplot(4,1,1); hold on;
scatter(1:n_depths,max_b3_early,100,plot_col,'filled','MarkerEdgeColor','k','linewidth',2);
scatter(1:n_depths,max_b3_late,100,plot_col,'filled','MarkerEdgeColor',[0.5,0.5,0.5],'linewidth',2);
set(gca,'XTick',1:n_depths);
xlabel('Striatal depth');
ylabel('\beta_3');

subplot(4,1,2); hold on;
scatter(1:n_depths,max_b4_early,100,plot_col,'filled','MarkerEdgeColor','k','linewidth',2);
scatter(1:n_depths,max_b4_late,100,plot_col,'filled','MarkerEdgeColor',[0.5,0.5,0.5],'linewidth',2);
set(gca,'XTick',1:n_depths);
xlabel('Striatal depth');
ylabel('\beta_4');

subplot(4,1,3); hold on;
scatter(1:n_depths,(max_b3_early-max_b4_early)./(max_b3_early+max_b4_early), ...
    100,plot_col,'filled','MarkerEdgeColor','k','linewidth',2);
scatter(1:n_depths,(max_b3_late-max_b4_late)./(max_b3_late+max_b4_late), ...
    100,plot_col,'filled','MarkerEdgeColor',[0.5,0.5,0.5],'linewidth',2);
set(gca,'XTick',1:n_depths);
xlabel('Striatal depth');
ylabel('(\beta_3-\beta_4)/(\beta_3+\beta_4)');

subplot(4,1,4); hold on;
scatter(max_b3_early,max_b4_early,100,plot_col,'filled','MarkerEdgeColor','k','linewidth',2);
scatter(max_b3_late,max_b4_late,100,plot_col,'filled','MarkerEdgeColor',[0.5,0.5,0.5],'linewidth',2);
xlabel('\beta_3');ylabel('\beta_4')
max_weight = max([max_b3_early;max_b3_late;max_b4_early;max_b4_late]);
line([0,max_weight],[0,max_weight],'color','k');
axis image





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

softnorm = 20;

use_animals = cellfun(@(x) ~isempty(x),{batch_vars(:).mua_stim});

mua_stim_smoothed = cellfun(@(x) convn(x,smWin,'same'),{batch_vars(use_animals).mua_stim},'uni',false);
mua_stim_norm = cellfun(@(x) bsxfun(@rdivide,x,nanmean(x(:,t_baseline,:,:),2)+softnorm),mua_stim_smoothed,'uni',false);
mua_stim_mean = cellfun(@(x) nanmean(x,4),mua_stim_norm,'uni',false);
mua_stim_combined = nanmean(cat(4,mua_stim_mean{:}),4);

trace_spacing = 0.8;
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

%% Align ephys recordings by passive correlations

% Define parameters
n_str_depths = 6;
save_alignment = false;

% Load data from passive visual experiments
data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\passive'];
data_fn = 'mua_passive';
load([data_path filesep data_fn]);

n_experiments = cellfun(@length,{batch_vars.mua_corr});

n_depths = cellfun(@(x) cellfun(@length,x),{batch_vars.mua_corr},'uni',false);
n_depths_cat = horzcat(n_depths{:});
max_grps = max(horzcat(n_depths{:}));

mua_pad = cellfun(@(x) cellfun(@(x) padarray(x,[max_grps-length(x),max_grps-length(x)],NaN,'pre'),x,'uni',false),{batch_vars.mua_corr},'uni',false);
lfp_pad = cellfun(@(x) cellfun(@(x) padarray(x,[max_grps-length(x),max_grps-length(x)],NaN,'pre'),x,'uni',false),{batch_vars.lfp_corr},'uni',false);
vis_modulation_pad = cellfun(@(x) cell2mat(cellfun(@(x) padarray(x,[max_grps-length(x),0],NaN,'pre'),x,'uni',false)),{batch_vars.vis_modulation},'uni',false);

mua_cat = cell2mat(cellfun(@(x) cat(3,x{:}),permute(mua_pad,[1,3,2]),'uni',false));
lfp_cat = cell2mat(cellfun(@(x) cat(3,x{:}),permute(lfp_pad,[1,3,2]),'uni',false));

vis_modulation_cat = horzcat(vis_modulation_pad{:});
vis_modulation_mean = nanmean(vis_modulation_cat,2);

% Align LFP to templates, then align visual modulation by that
template_size = max_grps*2;
template_lfp_groups = 2;
template_lfp_corr = zeros(template_size);
template_lfp_corr(eye(template_size) > 0) = 1;
template_lfp_edges = round(linspace(1,template_size,template_lfp_groups+1));
for curr_grp = 1:template_lfp_groups
    template_lfp_corr(template_lfp_edges(curr_grp):template_lfp_edges(curr_grp+1), ...
        template_lfp_edges(curr_grp):template_lfp_edges(curr_grp+1)) = 1;
end

lfp_bigpad = cellfun(@(x) cellfun(@(x) padarray(x, ...
    [template_size-length(x),template_size-length(x)],0,'post'), ...
    x,'uni',false),{batch_vars.lfp_corr},'uni',false);
lfp_bigcat = cell2mat(cellfun(@(x) cat(3,x{:}),permute(lfp_bigpad,[1,3,2]),'uni',false));

shift_grp = nan(size(lfp_bigcat,3),1);
aligned_lfp = nan(size(lfp_bigcat));
for curr_exp = 1:size(lfp_bigcat,3)
    n_shifts = template_size-n_depths_cat(curr_exp);
    template_dot = nan(n_shifts,1);
    for i = 1:n_shifts
        shifted_lfp = circshift(lfp_bigcat(:,:,curr_exp),[i,i]);
        template_dot(i) = transpose(shifted_lfp(:))*template_lfp_corr(:);
    end
    [~,shift_grp(curr_exp)] = max(template_dot);
    aligned_lfp(:,:,curr_exp) = circshift(lfp_bigcat(:,:,curr_exp),[shift_grp(curr_exp),shift_grp(curr_exp)]);
end

vis_modulation_bigpad = cellfun(@(x) cell2mat(cellfun(@(x) padarray(x, ...
    [template_size-length(x),0],NaN,'post'),x,'uni',false)), ...
    {batch_vars.vis_modulation},'uni',false);
vis_modulation_bigcat = horzcat(vis_modulation_bigpad{:});
aligned_vis_modulation = nan(size(vis_modulation_bigcat));
for curr_exp = 1:size(lfp_bigcat,3)
    aligned_vis_modulation(:,curr_exp) = ...
        circshift(vis_modulation_bigcat(:,curr_exp),[shift_grp(curr_exp),0]);
end

% From that alignment, divide the total extent into 6 groups
depth_bin_size = 100; % from the batch script
used_bins = any(~isnan(aligned_vis_modulation),2);
n_used_bins = sum(used_bins);
total_depth = depth_bin_size*n_used_bins;
str_depth_edges = linspace(0,total_depth,n_str_depths+1);

% Get the striatum "offset" for each experiment, package by animal
first_bin = find(used_bins,1);
[~,str_offset_bin] = max(~isnan(aligned_vis_modulation(first_bin:end,:)),[],1);
str_offset = mat2cell((str_offset_bin-1)*depth_bin_size,1,n_experiments);

% Plot the template and mean aligned LFP

% % Above didn't work - align by just assuming the end is the correct border
% 
% % Plot
% vis_modulation_cutoff = 0.5;
% max_vis = cellfun(@(x) cellfun(@max,x),{batch_vars.vis_modulation},'uni',false);
% 
% figure;
% curr_exp = 1;
% for curr_animal = 1:length(batch_vars)
%     for curr_day = 1:n_experiments(curr_animal)
%         subplot(length(batch_vars),max(n_experiments), ...
%             (curr_animal-1)*max(n_experiments)+curr_day);
%         imagesc(lfp_cat(:,:,curr_exp));
%         axis image off;
%         colormap(gray);
%         
%         if max_vis{curr_animal}(curr_day) > vis_modulation_cutoff
%             title(num2str(round(max_vis{curr_animal}(curr_day)*100)))
%         end
%         
%         curr_exp = curr_exp + 1;
%     end
% end

% Divide the total extent into N depths
depth_bin_size = 100; % from the batch script
total_depth = depth_bin_size*max_grps;

str_depth_edges = linspace(0,total_depth,n_str_depths+1);

figure;imagesc(1:size(vis_modulation_cat,2),(0:max_grps-1)*depth_bin_size,vis_modulation_cat)
for i = 1:length(str_depth_edges)
    line(xlim,repmat(str_depth_edges(i),2,1),'linewidth',2,'color','r')
end

% Get the striatum "offset" for each experiment, package by animal
str_offset_bin = max_grps-n_depths_cat+1;
str_offset = mat2cell((str_offset_bin-1)*depth_bin_size,1,n_experiments);

% Package in structure
n_animal = length(batch_vars);
ephys_align = struct('animal',cell(n_animal,1),'days',cell(n_animal,1), ...
    'str_depth_edges',cell(n_animal,1),'str_offset',cell(n_animal,1));
[ephys_align.animal] = deal(batch_vars.animal);
[ephys_align.days] = deal(batch_vars.day);
[ephys_align.str_depth_edges] = deal(str_depth_edges);
[ephys_align.str_offset] = deal(str_offset{:});

% To save
if save_alignment
    ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
    ephys_align_fn = ['ephys_str_align_' num2str(n_str_depths) '_depths'];
    save([ephys_align_path filesep ephys_align_fn],'ephys_align');
    disp('Saved ephys alignment');
end

%% Load and average wf -> ephys maps
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_ephys';

protocol = 'vanillaChoiceworld';
% protocol = 'stimSparseNoiseUncorrAsync';
% protocol = 'stimKalatsky'; 
% protocol = 'AP_choiceWorldStimPassive';

map_fn = [data_path filesep 'wf_ephys_maps_' protocol];
load(map_fn);

n_depths = size(batch_vars(1).explained_var{1},1);

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% (scale r_px's because different lambdas give different weights)
% (do in a loop because memory can't handle a cellfun??)
r_px = nan(437,416,43,n_depths,length(animals));
for curr_animal = 1:length(animals)
    curr_animal_r_px = nan(437,416,43,n_depths,length(days));
    for curr_day = 1:length(batch_vars(curr_animal).r_px)        
        
        curr_r_px = batch_vars(curr_animal).r_px{curr_day};
        curr_scaled_r_px = mat2gray(bsxfun(@rdivide, ...
            curr_r_px,permute(prctile(abs( ...
            reshape(curr_r_px,[],size(curr_r_px,5))),95)',[2,3,4,5,1])),[-1,1])*2-1;     
        
        % Set any NaN explained (no MUA data probably) to NaN
        curr_scaled_r_px(:,:,:,isnan(batch_vars(curr_animal).explained_var{curr_day})) = NaN;
        
        curr_animal_r_px(:,:,:,:,curr_day) = curr_scaled_r_px;
        
    end
    r_px(:,:,:,:,curr_animal) = nanmean(curr_animal_r_px,5);
    disp(curr_animal);
end

com = cell2mat(shiftdim(cellfun(@(x) nanmean(cat(3,x{:}),3),{batch_vars.r_px_com},'uni',false),-1));
weight = cell2mat(shiftdim(cellfun(@(x) nanmean(cat(3,x{:}),3),{batch_vars.r_px_weight},'uni',false),-1));
explained_var = cell2mat(cellfun(@(x) nanmean(cat(2,x{:}),2),{batch_vars.explained_var},'uni',false));

r_px_mean = nanmean(r_px,5);
% flip r_px to go forward in time
r_px_mean = r_px_mean(:,:,end:-1:1,:);
com_mean = nanmean(com,3);
weight_mean = nanmean(weight,3);

% Load widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);

% Plot the kernel matches
if isfield(batch_vars,'kernel_match')
    max_n_kernel_match = max(cell2mat(cellfun(@(x) ...
        cellfun(@length,x),{batch_vars.kernel_match},'uni',false)));
    
    kernel_match_raw_all = cellfun(@(x) cell2mat(cellfun(@(x) ...
        padarray(x,max_n_kernel_match-length(x),NaN,'pre'),x,'uni',false)), ...
        {batch_vars.kernel_match_raw},'uni',false);
    
    kernel_match_all = cellfun(@(x) cell2mat(cellfun(@(x) ...
        padarray(x,max_n_kernel_match-length(x),NaN,'pre'),x,'uni',false)), ...
        {batch_vars.kernel_match},'uni',false);
    
    figure;
    for curr_animal = 1:length(animals)
        subplot(2,length(animals),curr_animal); hold on;
        set(gca,'ColorOrder',copper(size(kernel_match_raw_all{curr_animal},2)));
        plot(kernel_match_raw_all{curr_animal},'linewidth',2);
        
        subplot(2,length(animals),curr_animal+length(animals)); hold on;
        set(gca,'ColorOrder',copper(size(kernel_match_all{curr_animal},2)));
        plot(kernel_match_all{curr_animal},'linewidth',2);
    end
end

% Plot weights over time (only use ipsi ROIs)
t = linspace(-0.3,0.3,size(r_px_mean,3));
AP_image_scroll(r_px_mean,t)  
axis image;
caxis([-1,1]);
colormap(colormap_BlueWhiteRed);
AP_reference_outline('ccf_aligned','k');AP_reference_outline('retinotopy','m');

r_px_max = squeeze(max(r_px_mean,[],3));

wf_roi_masks = cat(3,wf_roi(:,1).mask);
wf_roi_mask_mult = bsxfun(@rdivide,wf_roi_masks, ...
    permute(sum(reshape(wf_roi_masks,[],size(wf_roi,1)),1),[1,3,2]));

roi_traces = permute(reshape(reshape(r_px_mean,[],length(t)*n_depths)'* ...
    reshape(wf_roi_mask_mult,[],size(wf_roi,1)), ...
    length(t),n_depths,size(wf_roi,1)),[1,3,2]);

figure;
for curr_depth = 1:n_depths
    subplot(2,n_depths,curr_depth)
    imagesc(r_px_max(:,:,curr_depth));
    axis image off;
    caxis([-1,1]);
    colormap(colormap_BlueWhiteRed);
    AP_reference_outline('ccf_aligned','k');AP_reference_outline('retinotopy','m');
    title(['Str ' num2str(curr_depth)]);
    
    subplot(2,n_depths,n_depths+curr_depth); hold on;
    set(gca,'ColorOrder',autumn(size(wf_roi,1)));
    plot(t,roi_traces(:,:,curr_depth),'linewidth',2);
    axis tight;
    line([0,0],ylim,'color','k');
    xlabel('Time from spike');
    ylabel('Normalized fluor weight');
end

% Plot time asymmetry in maps
r_px_asym = r_px_mean(:,:,t>=0 & t<=0.2,:) - r_px_mean(:,:,fliplr(find(t<=0 & t >=-0.2)),:);
figure;
for curr_depth = 1:n_depths
    subplot(1,n_depths,curr_depth)
    imagesc(r_px_asym(:,:,3,curr_depth))
    axis image off;
    caxis([-0.5,0.5]);
    colormap(colormap_BlueWhiteRed);
    AP_reference_outline('ccf_aligned','k');AP_reference_outline('retinotopy','m');
    title(['Str ' num2str(curr_depth)]);   
end

% Plot map (from mean kernel)
use_t = t > -0.1 & t < 0.1;
r_px_max = squeeze(max(r_px_mean(:,:,use_t,:),[],3)).^3;
for i = 1:n_depths
    r_px_max(:,:,i) = medfilt2(r_px_max(:,:,i),[10,10]);
end
r_px_max_norm = bsxfun(@rdivide,r_px_max, ...
    permute(max(reshape(r_px_max,[],n_depths),[],1),[1,3,2]));
r_px_max_norm(isnan(r_px_max_norm)) = 0;
r_px_com = sum(bsxfun(@times,r_px_max_norm,permute(1:n_depths,[1,3,2])),3)./sum(r_px_max_norm,3);
com_colored = ind2rgb(round(mat2gray(r_px_com,[1,n_depths])*255),jet(255));

figure;
p = imagesc(com_colored);
axis image off
weight_norm = mat2gray(max(r_px_max,[],3),[0,double(prctile(reshape(max(r_px_max,[],3),[],1),90))]);
set(p,'AlphaData',weight_norm);

c = colorbar;
ylabel(c,'Depth (fraction)');
colormap(c,jet);
set(c,'YDir','reverse');
set(c,'YTick',linspace(0,1,n_depths));
set(c,'YTickLabel',1:n_depths);
 
AP_reference_outline('retinotopy','m');
AP_reference_outline('ccf_aligned','k');

% %%%% SAVE TEMPLATE KERNELS
% kernel_template_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template';
% kernel_template = r_px_max_norm;
% save(kernel_template_fn,'kernel_template');
% disp('Saved kernel template');
% %%%%

% Plot map (from mean CoM)
com_leeway =1;
c_range = [1+com_leeway,n_depths-com_leeway];
com_colored = ind2rgb(round(mat2gray(com_mean,c_range)*255),jet(255));

figure;
p = imagesc(com_colored);
axis image off
weight_norm = mat2gray(max(weight_mean,[],3),[0,double(prctile(reshape(max(weight_mean,[],3),[],1),98))]);
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

%% Plot kernel matches

% Load the kernel template matches
kernel_match_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
kernel_match_fn = 'ephys_kernel_align';
load([kernel_match_path filesep kernel_match_fn]);

% Plot the kernel matches
max_n_kernel_match = max(cell2mat(cellfun(@(x) ...
    cellfun(@length,x),{ephys_kernel_align.kernel_match},'uni',false)));

kernel_match_raw_all = cellfun(@(x) cell2mat(cellfun(@(x) ...
    padarray(x,max_n_kernel_match-length(x),NaN,'pre'),x,'uni',false)), ...
    {ephys_kernel_align.kernel_match_raw},'uni',false);

kernel_match_all = cellfun(@(x) cell2mat(cellfun(@(x) ...
    padarray(x,max_n_kernel_match-length(x),NaN,'pre'),x,'uni',false)), ...
    {ephys_kernel_align.kernel_match},'uni',false);

figure;
for curr_animal = 1:length(ephys_kernel_align)
    subplot(2,length(ephys_kernel_align),curr_animal); hold on;
    set(gca,'ColorOrder',copper(size(kernel_match_raw_all{curr_animal},2)));
    plot(kernel_match_raw_all{curr_animal},'linewidth',2);
    
    subplot(2,length(ephys_kernel_align),curr_animal+length(ephys_kernel_align)); hold on;
    set(gca,'ColorOrder',copper(size(kernel_match_all{curr_animal},2)));
    plot(kernel_match_all{curr_animal},'linewidth',2);
end

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

softnorm = 20;

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

%% Load and process cortex-predicted striatal MUA during choiceworld (NEW)

data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
mua_fn = [data_path filesep 'mua_choiceworld_predicted_str-ctx'];
load(mua_fn);

n_depths = 6;

contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';

sample_rate = 35.2*2;
interval_surround = [-0.5,3];
t = interval_surround(1):1/sample_rate:interval_surround(2);

% Normalize, smooth, and combine PSTH
% (PSTH: [condition,time,depth,align,day])
t_baseline = t < 0;
softnorm = 500;

depth_psth = {batch_vars.depth_psth};
mua_baseline = cellfun(@(x) nanmean(x(:,t_baseline,:,1,:),2),depth_psth,'uni',false);
depth_psth_norm = cellfun(@(mua,baseline) ...
    bsxfun(@rdivide,mua,baseline + softnorm),depth_psth,mua_baseline,'uni',false);
depth_psth_mean = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmean(x,5),depth_psth_norm,'uni',false),[1,3,4,5,2])),5);

predicted_depth_psth = {batch_vars.predicted_depth_psth};
predicted_mua_baseline = cellfun(@(x) nanmean(x(:,t_baseline,:,1,:),2),predicted_depth_psth,'uni',false);
predicted_depth_psth_norm = cellfun(@(mua,baseline) ...
    bsxfun(@rdivide,mua,baseline + softnorm),predicted_depth_psth,predicted_mua_baseline,'uni',false);
predicted_depth_psth_mean = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmean(x,5),predicted_depth_psth_norm,'uni',false),[1,3,4,5,2])),5);

% Smooth the real psth to match the inherently smoothed predicted
% (by finding the best match)
start_params = [9,3];
options = struct;
options.MaxFunEvals = 1000;
options.Display = 'off';

model_sse = @(P) ...
    sum((predicted_depth_psth_mean(:) - ...
    reshape(convn(depth_psth_mean, ...
    (gausswin(P(1),P(2))'./sum(gausswin(P(1),P(2))')),'same'),[],1)).^2);

warning off
smooth_params = round(fminsearch(model_sse,start_params,options));
warning on

gw = gausswin(smooth_params(1),smooth_params(2))';
smWin = gw./sum(gw);

depth_psth_mean = convn(depth_psth_mean,smWin,'same');

% Plot all hit/miss conditions for one depth
plot_str = 3;
plot_success = -1;
plot_align = 2;
plot_timing = 1;
plot_color = colormap_BlueWhiteRed(6);
plot_color = [[0,0,0];plot_color(5:-1:1,:);[0,0,0];plot_color(end-5:end,:)];

figure; 
p1 = subplot(1,3,1); hold on;
set(gca,'ColorOrder',plot_color);
plot_conditions = conditions(:,2) == -plot_success*conditions(:,3) & conditions(:,4) == plot_timing;
plot(t,depth_psth_mean(plot_conditions,:,plot_str,plot_align)','linewidth',2);
line([0,0],ylim,'color','k');
ylabel(['Str ' num2str(plot_str)]);
xlabel('Time from stim')
title('Real')

p2 = subplot(1,3,2); hold on;
set(gca,'ColorOrder',plot_color);
plot_conditions = conditions(:,2) == -plot_success*conditions(:,3) & conditions(:,4) == plot_timing;
plot(t,predicted_depth_psth_mean(plot_conditions,:,plot_str,plot_align)','linewidth',2);
line([0,0],ylim,'color','k');
ylabel(['Str ' num2str(plot_str)]);
xlabel('Time from stim')
title('Predicted')

p3 = subplot(1,3,3); hold on;
set(gca,'ColorOrder',plot_color);
plot_conditions = conditions(:,2) == -plot_success*conditions(:,3) & conditions(:,4) == plot_timing;
plot(t,depth_psth_mean(plot_conditions,:,plot_str,plot_align)' - ...
    predicted_depth_psth_mean(plot_conditions,:,plot_str,plot_align)','linewidth',2);
line([0,0],ylim,'color','k');
ylabel(['Str ' num2str(plot_str)]);
xlabel('Time from stim')
title('Real - Predicted')

linkaxes([p1,p2,p3],'x');
linkaxes([p1,p2],'y');


% % Fit contrast response
% % [depth,time,align,timing,param,success]
% activity_model_params = nan(n_depths,length(t),2,2,3,2);
% predicted_activity_model_params = nan(n_depths,length(t),2,2,3,2);
% for curr_success = 1:2
%     switch curr_success
%         case 1
%             use_success = 1;
%         case 2
%             use_success = -1;
%     end
%     
%     for curr_timing = 1:2
%         
%         use_conditions = conditions(:,2) == -use_success*conditions(:,3) & conditions(:,4) == curr_timing;
%         contrast_sides = zeros(size(conditions,1),2);
%         contrast_sides(conditions(:,2) == -1,1) = conditions(conditions(:,2) == -1,1);
%         contrast_sides(conditions(:,2) == 1,2) = conditions(conditions(:,2) == 1,1);
%         use_contrast_sides = contrast_sides(use_conditions,:);
%         
%         %     activity_model = @(P) P(1) + P(2).*use_contrast_sides(:,1).^P(4) + P(3).*use_contrast_sides(:,2).^P(4);
%         activity_model = @(P) P(1) + P(2).*use_contrast_sides(:,1).^0.3 + P(3).*use_contrast_sides(:,2).^0.3;
%         
%         for curr_depth = 1:n_depths
%             for curr_align = 1:2
%                 for curr_t = 1:length(t)
%                     
%                     % Real
%                     curr_act = depth_psth_mean(use_conditions,curr_t,curr_depth,curr_align);
%                     model_sse = @(P) sum((curr_act-(activity_model(P))).^2);
%                     
%                     start_params = [0,1,0.3];
%                     
%                     options = struct;
%                     options.MaxFunEvals = 1000;
%                     options.Display = 'off';
%                     
%                     params_fit = fminsearch(model_sse,start_params,options);
%                     activity_model_params(curr_depth,curr_t,curr_align,curr_timing,:,curr_success) = ...
%                         params_fit;
%                     
%                     % Predicted
%                     curr_act = predicted_depth_psth_mean(use_conditions,curr_t,curr_depth,curr_align);
%                     model_sse = @(P) sum((curr_act-(activity_model(P))).^2);
%                     
%                     start_params = [0,1,0.3];
%                     
%                     options = struct;
%                     options.MaxFunEvals = 1000;
%                     options.Display = 'off';
%                     
%                     params_fit = fminsearch(model_sse,start_params,options);
%                     predicted_activity_model_params(curr_depth,curr_t,curr_align,curr_timing,:,curr_success) = ...
%                         params_fit;
%                     
%                 end
%             end
%             disp(curr_depth);
%         end
%     end
% end
% 
% spacing = 0.3;
% plot_param = 1;
% 
% figure; 
% s1 = subplot(1,2,1);hold on;
% p1 = AP_stackplot(activity_model_params(:,:,1,1,plot_param,1)',t,spacing,false,'k');
% p2 = AP_stackplot(predicted_activity_model_params(:,:,1,1,plot_param,1)',t,spacing,false,'r');
% p3 = AP_stackplot(activity_model_params(:,:,1,1,plot_param,2)',t,spacing,false,'b');
% p4 = AP_stackplot(predicted_activity_model_params(:,:,1,1,plot_param,2)',t,spacing,false,'m');
% line([0,0],ylim,'color','k');
% line([0.2,0.2],ylim,'color','k');
% xlabel('Time from stim onset')
% ylabel(['Param ' num2str(plot_param)])
% title('Early move')
% legend([p1(1),p2(1),p3(1),p4(1)],{'Real hit','Predicted hit','Real miss','Predicted miss'});
% 
% s2 = subplot(1,2,2);hold on;
% p1 = AP_stackplot(activity_model_params(:,:,1,2,plot_param,1)',t,spacing,false,'k');
% p2 = AP_stackplot(predicted_activity_model_params(:,:,1,2,plot_param,1)',t,spacing,false,'r');
% p3 = AP_stackplot(activity_model_params(:,:,1,2,plot_param,2)',t,spacing,false,'b');
% p4 = AP_stackplot(predicted_activity_model_params(:,:,1,2,plot_param,2)',t,spacing,false,'m');
% line([0,0],ylim,'color','k');
% line([0.5,0.5],ylim,'color','k');
% line([0.7,0.7],ylim,'color','k');
% xlabel('Time from stim onset')
% ylabel(['Param ' num2str(plot_param)])
% title('Late move')
% legend([p1(1),p2(1),p3(1),p4(1)],{'Real hit','Predicted hit','Real miss','Predicted miss'});
% 
% linkaxes([s1,s2,s3,s4]);

% Fit contrast response plus choice
contrast_exp = 0.3;

% [depth,time,align,timing,param,success]
activity_model_params = nan(n_depths,length(t),2,2,4);
predicted_activity_model_params = nan(n_depths,length(t),2,2,4);

for curr_timing = 1:2
    
    use_conditions = conditions(:,4) == curr_timing;
    contrast_sides = zeros(size(conditions,1),2);
    contrast_sides(conditions(:,2) == -1,1) = conditions(conditions(:,2) == -1,1);
    contrast_sides(conditions(:,2) == 1,2) = conditions(conditions(:,2) == 1,1);
    
    use_choices = conditions(use_conditions,3);
    use_contrast_sides = contrast_sides(use_conditions,:);
    
    for curr_depth = 1:n_depths
        for curr_align = 1:2
            for curr_t = 1:length(t)
                
                curr_act = depth_psth_mean(use_conditions,curr_t,curr_depth,curr_align);
                params_fit = [ones(size(curr_act)),use_contrast_sides.^contrast_exp,use_choices]\curr_act;
                activity_model_params(curr_depth,curr_t,curr_align,curr_timing,:) = ...
                    params_fit;
                
                curr_act = predicted_depth_psth_mean(use_conditions,curr_t,curr_depth,curr_align);
                params_fit = [ones(size(curr_act)),use_contrast_sides.^contrast_exp,use_choices]\curr_act;
                predicted_activity_model_params(curr_depth,curr_t,curr_align,curr_timing,:) = ...
                    params_fit;
                
            end
        end
    end
end

plot_param = 4;
spacing = max(abs(reshape(activity_model_params(:,:,:,:,plot_param),[],1)))*1.5;

figure; 
s1 = subplot(1,2,1);hold on;
p1 = AP_stackplot(activity_model_params(:,:,1,1,plot_param)',t,spacing,false,'k');
p2 = AP_stackplot(predicted_activity_model_params(:,:,1,1,plot_param)',t,spacing,false,'r');
p3 = AP_stackplot(activity_model_params(:,:,1,2,plot_param)',t,spacing,false,'b');
p4 = AP_stackplot(predicted_activity_model_params(:,:,1,2,plot_param)',t,spacing,false,'m');
line([0,0],ylim,'color','k');
line([0.5,0.5],ylim,'color','k');
xlabel('Time from stim')
ylabel(['Param ' num2str(plot_param)])
legend([p1(1),p2(1),p3(1),p4(1)],{'Early real','Early predicted','Late real','Late predicted'});

s2 = subplot(1,2,2);hold on;
p1 = AP_stackplot(activity_model_params(:,:,2,1,plot_param)',t,spacing,false,'k');
p2 = AP_stackplot(predicted_activity_model_params(:,:,2,1,plot_param)',t,spacing,false,'r');
p3 = AP_stackplot(activity_model_params(:,:,2,2,plot_param)',t,spacing,false,'b');
p4 = AP_stackplot(predicted_activity_model_params(:,:,2,2,plot_param)',t,spacing,false,'m');
line([0,0],ylim,'color','k');
xlabel('Time from move')
ylabel(['Param ' num2str(plot_param)])
legend([p1(1),p2(1),p3(1),p4(1)],{'Early real','Early predicted','Late real','Late predicted'});

linkaxes([s1,s2]);


%%% Params from trial activity within day

% [depth,time,align,timing,param]
activity_model_params = {batch_vars.activity_model_params};
activity_model_params_mean = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmedian(x,6),activity_model_params,'uni',false),[1,3,4,5,6,2])),6);
activity_model_params_mean_smoothed = convn(activity_model_params_mean,smWin,'same');

predicted_activity_model_params = {batch_vars.predicted_activity_model_params};
predicted_activity_model_params_mean = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmean(x,6),predicted_activity_model_params,'uni',false),[1,3,4,5,6,2])),6);

plot_param = 3;

spacing = max(abs(reshape(activity_model_params_mean_smoothed(:,:,:,:,plot_param),[],1)))*1.5;

figure; 
s1 = subplot(1,2,1);hold on;
p1 = AP_stackplot(activity_model_params_mean_smoothed(:,:,1,1,plot_param)',t,spacing,false,'k');
p2 = AP_stackplot(predicted_activity_model_params_mean(:,:,1,1,plot_param)',t,spacing,false,'r');
p3 = AP_stackplot(activity_model_params_mean_smoothed(:,:,1,2,plot_param)',t,spacing,false,'b');
p4 = AP_stackplot(predicted_activity_model_params_mean(:,:,1,2,plot_param)',t,spacing,false,'m');
line([0,0],ylim,'color','k');
line([0.5,0.5],ylim,'color','k');
xlabel('Time from stim')
ylabel(['Param ' num2str(plot_param)])
legend([p1(1),p2(1),p3(1),p4(1)],{'Early real','Early predicted','Late real','Late predicted'});

s2 = subplot(1,2,2);hold on;
p1 = AP_stackplot(activity_model_params_mean_smoothed(:,:,2,1,plot_param)',t,spacing,false,'k');
p2 = AP_stackplot(predicted_activity_model_params_mean(:,:,2,1,plot_param)',t,spacing,false,'r');
p3 = AP_stackplot(activity_model_params_mean_smoothed(:,:,2,2,plot_param)',t,spacing,false,'b');
p4 = AP_stackplot(predicted_activity_model_params_mean(:,:,2,2,plot_param)',t,spacing,false,'m');
line([0,0],ylim,'color','k');
xlabel('Time from move')
ylabel(['Param ' num2str(plot_param)])
legend([p1(1),p2(1),p3(1),p4(1)],{'Early real','Early predicted','Late real','Late predicted'});

linkaxes([s1,s2]);

% Plot unsmoothed
plot_param = 4;
spacing = max(abs(reshape(activity_model_params_mean(:,:,:,:,plot_param),[],1)))*1.5;

figure; 
s1 = subplot(1,2,1);hold on;
p1 = AP_stackplot(activity_model_params_mean(:,:,1,1,plot_param)',t,spacing,false,'k');
p2 = AP_stackplot(activity_model_params_mean(:,:,1,2,plot_param)',t,spacing,false,'b');
line([0,0],ylim,'color','k');
line([0.5,0.5],ylim,'color','k');
xlabel('Time from stim')
ylabel(['Param ' num2str(plot_param)])
legend([p1(1),p2(1)],{'Early real','Late real'});

s2 = subplot(1,2,2);hold on;
p1 = AP_stackplot(activity_model_params_mean(:,:,2,1,plot_param)',t,spacing,false,'k');
p2 = AP_stackplot(activity_model_params_mean(:,:,2,2,plot_param)',t,spacing,false,'b');
line([0,0],ylim,'color','k');
xlabel('Time from move')
ylabel(['Param ' num2str(plot_param)])
legend([p1(1),p2(1)],{'Early real','Late real'});

linkaxes([s1,s2]);


% Plot max b3 v max b4
use_t_stim = t > 0.05 & t < 0.15;
use_t_move = t > -0.15 & t < -0.02;

max_b3_early = max(max(abs(activity_model_params_mean(:,use_t_stim,1,1,3)),[],2), ...
    max(abs(activity_model_params_mean(:,use_t_move,2,1,3)),[],2));
max_b3_late = max(max(abs(activity_model_params_mean(:,use_t_stim,1,2,3)),[],2), ...
    max(abs(activity_model_params_mean(:,use_t_move,2,2,3)),[],2));

max_b4_early = max(max(abs(activity_model_params_mean(:,use_t_stim,1,1,4)),[],2), ...
    max(abs(activity_model_params_mean(:,use_t_move,2,1,4)),[],2));
max_b4_late = max(max(abs(activity_model_params_mean(:,use_t_stim,1,2,4)),[],2), ...
    max(abs(activity_model_params_mean(:,use_t_move,2,2,4)),[],2));

figure; 
plot_col = copper(n_depths);

subplot(4,1,1); hold on;
scatter(1:n_depths,max_b3_early,100,plot_col,'filled','MarkerEdgeColor','k','linewidth',2);
scatter(1:n_depths,max_b3_late,100,plot_col,'filled','MarkerEdgeColor',[0.5,0.5,0.5],'linewidth',2);
set(gca,'XTick',1:n_depths);
xlabel('Striatal depth');
ylabel('\beta_3');

subplot(4,1,2); hold on;
scatter(1:n_depths,max_b4_early,100,plot_col,'filled','MarkerEdgeColor','k','linewidth',2);
scatter(1:n_depths,max_b4_late,100,plot_col,'filled','MarkerEdgeColor',[0.5,0.5,0.5],'linewidth',2);
set(gca,'XTick',1:n_depths);
xlabel('Striatal depth');
ylabel('\beta_4');

subplot(4,1,3); hold on;
scatter(1:n_depths,(max_b3_early-max_b4_early)./(max_b3_early+max_b4_early), ...
    100,plot_col,'filled','MarkerEdgeColor','k','linewidth',2);
scatter(1:n_depths,(max_b3_late-max_b4_late)./(max_b3_late+max_b4_late), ...
    100,plot_col,'filled','MarkerEdgeColor',[0.5,0.5,0.5],'linewidth',2);
set(gca,'XTick',1:n_depths);
xlabel('Striatal depth');
ylabel('(\beta_3-\beta_4)/(\beta_3+\beta_4)');

subplot(4,1,4); hold on;
scatter(max_b3_early,max_b4_early,100,plot_col,'filled','MarkerEdgeColor','k','linewidth',2);
scatter(max_b3_late,max_b4_late,100,plot_col,'filled','MarkerEdgeColor',[0.5,0.5,0.5],'linewidth',2);
xlabel('\beta_3');ylabel('\beta_4')
max_weight = max([max_b3_early;max_b3_late;max_b4_early;max_b4_late]);
line([0,max_weight],[0,max_weight],'color','k');
axis image



%% Load and process MUA-fluor correlations

interval_surround = [-0.5,1.5];
t = linspace(interval_surround(1),interval_surround(2),212);
sample_rate = 1/median(diff(t));

plot_t = [-0.2,0.7];
t_use = t > 0.5 & t < 0.6;

% Load correlations
trialtype_align = 'move';
trialtype_timing = 'earlymove';

corr_use = ['corr_mua_fluor_' trialtype_align '_' trialtype_timing '_conditionshuff'];
fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\' corr_use];
load(fn);
n_depths = 6;

% Load widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
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
% gw = gausswin(smooth_size,3)';
% smWin = gw./sum(gw);
smWin = ones(1,smooth_size)/smooth_size;

figure; 

subplot(1,2,1); hold on
set(gca,'ColorOrder',copper(n_depths));
plot(t,conv2(horzcat(corr_mua_choice_split{:}),smWin','same'),'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Correlation with decision')
xlabel(['Time from ' trialtype_align ' onset (s)']);
title('MUA');
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:6),'uni',false))

subplot(1,2,2); hold on
set(gca,'ColorOrder',[autumn(size(wf_roi,1));winter(size(wf_roi,1))]);
plot(t,conv2(horzcat(corr_fluor_choice_split{:}),smWin','same'),'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Correlation with decision')
xlabel(['Time from ' trialtype_align ' onset (s)']);
title('Fluor');
legend(wf_areas);

% Plot predicted/predictor within window and lag
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

%% Choice correlation: plot individual animals (SUBSET OF ABOVE)

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

choice_animals = {batch_vars(:).corr_mua_choice};
choice_day = cell2mat(permute(cellfun(@(x) cell2mat(permute(x,[2,1,3])),choice_animals,'uni',false),[1,3,2]));
animalmean = cellfun(@(x) nanmean(cell2mat(permute(x,[2,1,3])),3),choice_animals,'uni',false);
animalmean = cat(3,animalmean{:});

smooth_size = 5;
smWin = ones(1,smooth_size)/smooth_size;

choice_day_smooth = convn(choice_day,smWin','same');
choice_animalmean_smooth = convn(animalmean,smWin','same');

plot_area = 2;

figure;
subplot(1,2,1); hold on
curr_plot = squeeze(choice_day_smooth(:,plot_area,:));
plot(t,curr_plot)
plot(t,nanmean(curr_plot,2),'k','linewidth',2)
axis tight;
line([0,0],ylim,'color','k');
xlabel(['Time from ' trialtype_align ' onset (s)']);
legend({'Animal-day'})

subplot(1,2,2); hold on;
curr_plot = squeeze(choice_animalmean_smooth(:,plot_area,:));
plot(t,curr_plot)
plot(t,nanmedian(curr_plot,2),'k','linewidth',2)
axis tight;
line([0,0],ylim,'color','k');
xlabel(['Time from ' trialtype_align ' onset (s)']);
legend(animals)


%% Feed-forward/back (first run above)

max_lag = 0.08;
plot_t = [-0.5,1.5];

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
            if ~isnan(curr_corr)
                set(fluor_fluor_h(curr_roi1,curr_roi2), ...
                    'linewidth',abs(curr_corr)*50, ...
                    'color',[sign(curr_corr) == 1,0,sign(curr_corr) == -1]);
            end
        end
    end
    % MUA-MUA
    for curr_depth1 = 1:n_depths
        for curr_depth2 = curr_depth1+1:n_depths
            curr_corr = local_corr_mua_mua_diff{curr_depth1,curr_depth2}(curr_t);
            if ~isnan(curr_corr)
                set(mua_mua_h(curr_depth1,curr_depth2), ...
                    'linewidth',abs(curr_corr)*120, ...
                    'color',[sign(curr_corr) == 1,0,sign(curr_corr) == -1]);
            end
        end
    end
    % Fluor-MUA
    for curr_roi = 1:n_rois
        for curr_depth = 1:n_depths
            curr_corr = local_corr_fluor_mua_diff{curr_roi,curr_depth}(curr_t);
            if ~isnan(curr_corr)
                set(fluor_mua_h(curr_roi,curr_depth), ...
                    'linewidth',abs(curr_corr)*50, ...
                    'color',[sign(curr_corr) == 1,0,sign(curr_corr) == -1]);
            end
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

close(fig_h);

% % this was doing PCA on correlation patterns, makes the labels
% str_labels = cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false);
% fluor_mua_labels = cellfun(@(fluor,mua) [fluor '-' mua], ...
%     repmat({wf_roi.area}',1,n_depths),repmat(str_labels,n_rois,1),'uni',false);
% 
% fluor_fluor_labels = cellfun(@(mua1,mua2) [mua1 '-' mua2], ...
%     repmat({wf_roi.area}',1,n_rois),repmat({wf_roi.area}',1,n_rois)','uni',false);
% 
% mua_mua_labels = cellfun(@(mua1,mua2) [mua1 '-' mua2], ...
%     repmat(str_labels,n_depths,1),repmat(str_labels,n_depths,1)','uni',false);




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

%% Use scales from above to combine feedforward/back plots (first run above)

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

local_corr_fluor_fluor_diff = cell(n_rois,n_rois,2);
local_corr_mua_mua_diff = cell(n_depths,n_depths,2);
local_corr_fluor_mua_diff = cell(n_rois,n_depths,2);

corr_mua_choice_cat = cell(2,1);
corr_fluor_choice_cat = cell(2,1);

for curr_align = 1:2
    
    % Load correlations
    load([corr_path filesep corr_fn{curr_align}]);
    
    % Create average area correlations
    corr_mua_mua = nanmean(cell2mat(permute(arrayfun(@(x) ...
        nanmean(cell2mat(batch_vars(x).corr_mua_mua),3),1:length(batch_vars),'uni',false),[1,3,2])),3);
    corr_fluor_fluor = nanmean(cell2mat(permute(arrayfun(@(x) ...
        nanmean(cell2mat(batch_vars(x).corr_fluor_fluor),3),1:length(batch_vars),'uni',false),[1,3,2])),3);
    corr_fluor_mua = nanmean(cell2mat(permute(arrayfun(@(x) ...
        nanmean(cell2mat(batch_vars(x).corr_fluor_mua),3),1:length(batch_vars),'uni',false),[1,3,2])),3);
 
    corr_mua_mua_split = mat2cell(corr_mua_mua,repmat(length(t),n_depths,1),repmat(length(t),1,n_depths));
    corr_fluor_fluor_split = mat2cell(corr_fluor_fluor,repmat(length(t),n_rois,1),repmat(length(t),1,n_rois));
    corr_fluor_mua_split = mat2cell(corr_fluor_mua,repmat(length(t),n_rois,1),repmat(length(t),1,n_depths));
  
    % Compute running feedforward/back
    max_lag = 0.08;
    max_lag_samples = round(max_lag*sample_rate);
    
    local_corr_fluor_fluor = cellfun(@(curr_data) cell2mat(arrayfun(@(x) padarray(diag(curr_data,x),[abs(x),0],NaN,'post'), ...
        -max_lag_samples:max_lag_samples,'uni',false)),corr_fluor_fluor_split,'uni',false);
    local_corr_fluor_fluor_diff(:,:,curr_align) = ...
        cellfun(@(local_corr) nanmean(local_corr(:,max_lag_samples+2:end) - ...
        local_corr(:,max_lag_samples:-1:1),2),local_corr_fluor_fluor,'uni',false);
    
    local_corr_mua_mua = cellfun(@(curr_data) cell2mat(arrayfun(@(x) padarray(diag(curr_data,x),[abs(x),0],NaN,'post'), ...
        -max_lag_samples:max_lag_samples,'uni',false)),corr_mua_mua_split,'uni',false);
    local_corr_mua_mua_diff(:,:,curr_align) = ...
        cellfun(@(local_corr) nanmean(local_corr(:,max_lag_samples+2:end) - ...
        local_corr(:,max_lag_samples:-1:1),2),local_corr_mua_mua,'uni',false);
    
    local_corr_fluor_mua = cellfun(@(curr_data) cell2mat(arrayfun(@(x) padarray(diag(curr_data,x),[abs(x),0],NaN,'post'), ...
        -max_lag_samples:max_lag_samples,'uni',false)),corr_fluor_mua_split,'uni',false);
    local_corr_fluor_mua_diff(:,:,curr_align) = ...
        cellfun(@(local_corr) nanmean(local_corr(:,max_lag_samples+2:end) - ...
        local_corr(:,max_lag_samples:-1:1),2),local_corr_fluor_mua,'uni',false);
    
    % Compute average choice correlations
    corr_mua_choice = nanmean(cell2mat(permute(arrayfun(@(x) ...
        nanmean(cell2mat(batch_vars(x).corr_mua_choice),3),1:length(batch_vars),'uni',false),[1,3,2])),3);
    corr_fluor_choice = nanmean(cell2mat(permute(arrayfun(@(x) ...
        nanmean(cell2mat(batch_vars(x).corr_fluor_choice),3),1:length(batch_vars),'uni',false),[1,3,2])),3);
    
    corr_mua_choice_split = mat2cell(corr_mua_choice,repmat(length(t),n_depths,1),1);
    corr_fluor_choice_split = mat2cell(corr_fluor_choice,repmat(length(t),n_rois,1),1);
        
    corr_mua_choice_cat{curr_align} = horzcat(corr_mua_choice_split{:});
    corr_fluor_choice_cat{curr_align} = horzcat(corr_fluor_choice_split{:});
        
end

est_move_time = 0.2;
t_plot = -0.2:1/sample_rate:0.4;

wf_ratio = bsxfun(@rdivide,wf_max,sum(wf_max,2));
str_ratio = bsxfun(@rdivide,str_max,sum(str_max,2));

% Scale and combine area correlations
fluor_fluor(:,:,1) = cellfun(@(x) interp1(t,x,t_plot),local_corr_fluor_fluor_diff(:,:,1),'uni',false);
fluor_fluor(:,:,2) = cellfun(@(x) interp1(t+est_move_time,x,t_plot),local_corr_fluor_fluor_diff(:,:,2),'uni',false);
[roi2,roi1] = meshgrid(1:n_rois);
fluor_fluor_scaled = cellfun(@(data1,data2,roi1,roi2) ...
    data1*wf_ratio(roi1,1)+data2*wf_ratio(roi2,2), ...
    fluor_fluor(:,:,1),fluor_fluor(:,:,2),num2cell(roi1),num2cell(roi2),'uni',false);

mua_mua(:,:,1) = cellfun(@(x) interp1(t,x,t_plot),local_corr_mua_mua_diff(:,:,1),'uni',false);
mua_mua(:,:,2) = cellfun(@(x) interp1(t+est_move_time,x,t_plot),local_corr_mua_mua_diff(:,:,2),'uni',false);
[depth2,depth1] = meshgrid(1:n_depths);
fluor_fluor_scaled = cellfun(@(data1,data2,depth1,depth2) ...
    data1*str_ratio(depth1,1)+data2*str_ratio(depth2,2), ...
    mua_mua(:,:,1),mua_mua(:,:,2),num2cell(depth1),num2cell(depth2),'uni',false);

fluor_mua(:,:,1) = cellfun(@(x) interp1(t,x,t_plot),local_corr_fluor_mua_diff(:,:,1),'uni',false);
fluor_mua(:,:,2) = cellfun(@(x) interp1(t+est_move_time,x,t_plot),local_corr_fluor_mua_diff(:,:,2),'uni',false);
[depth,roi] = meshgrid(1:n_depths,1:n_rois);
fluor_fluor_scaled = cellfun(@(data1,data2,roi,depth) ...
    data1*mean([wf_ratio(roi,1),str_ratio(depth,1)]) + ...
    data2*mean([wf_ratio(roi,2),str_ratio(depth,2)]), ...
    fluor_mua(:,:,1),fluor_mua(:,:,2),num2cell(roi),num2cell(depth),'uni',false);

% Scale and combine choice correlations
fluor_choice{1} = interp1(t,corr_fluor_choice_cat{1},t_plot);
fluor_choice{2} = interp1(t+est_move_time,corr_fluor_choice_cat{2},t_plot);
fluor_choice_scaled = sum(bsxfun(@times,cat(3,fluor_choice{:}),permute(wf_ratio,[3,1,2])),3);

mua_choice{1} = interp1(t,corr_mua_choice_cat{1},t_plot);
mua_choice{2} = interp1(t+est_move_time,corr_mua_choice_cat{2},t_plot);
mua_choice_scaled = sum(bsxfun(@times,cat(3,mua_choice{:}),permute(str_ratio,[3,1,2])),3);

% Plot scaled choice correlations
fig_h = figure('Position',[31,109,1838,828]);

smooth_size = 5;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

fig_h = figure('Position',[31,109,1838,828]);
subplot(2,1,1);hold on; 

set(gca,'ColorOrder',[copper(n_depths);autumn(size(wf_roi,1));winter(size(wf_roi,1))]);
plot(t_plot,conv2(mua_choice_scaled,smWin','same'),'linewidth',2);
plot(t_plot,conv2(fluor_choice_scaled,smWin','same'),'linewidth',2);

axis tight;
line([0,0],ylim,'color','k');
ylabel('Correlation with decision')
xlabel('Time from stim onset (s)')

t_h = line([0,0],ylim,'color','k');
t_mov = line(repmat(est_move_time,2,1),ylim,'color','k');

% Plot scaled area correlations
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
for curr_t = 1:length(t_plot)
    
    set(t_h,'XData',repmat(t_plot(curr_t),2,1));
    
    % Fluor-fluor
    for curr_roi1 = 1:n_rois
        for curr_roi2 = curr_roi1+1:n_rois
            curr_corr = fluor_fluor{curr_roi1,curr_roi2}(curr_t);
            set(fluor_fluor_h(curr_roi1,curr_roi2), ...
                'linewidth',abs(curr_corr)*50, ...
                'color',[sign(curr_corr) == 1,0,sign(curr_corr) == -1]);  
        end
    end
    % MUA-MUA
    for curr_depth1 = 1:n_depths
        for curr_depth2 = curr_depth1+1:n_depths
            curr_corr = mua_mua{curr_depth1,curr_depth2}(curr_t);
            set(mua_mua_h(curr_depth1,curr_depth2), ...
                'linewidth',abs(curr_corr)*120, ...
                'color',[sign(curr_corr) == 1,0,sign(curr_corr) == -1]);  
        end
    end
    % Fluor-MUA
    for curr_roi = 1:n_rois
        for curr_depth = 1:n_depths
            curr_corr = fluor_mua{curr_roi,curr_depth}(curr_t);
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
save_filename = [save_path filesep 'corr_combined.avi'];

writerObj = VideoWriter(save_filename);
writerObj.FrameRate = sample_rate/5;
open(writerObj);
writeVideo(writerObj,movie_frames);
close(writerObj);

%% Choice correlation: p-value of L/R difference

% Load batch analysis
fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\choice_p';
load(fn);

% Load widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_depths = 6;

interval_surround = [-0.5,1.5];
t = linspace(interval_surround(1),interval_surround(2),212);
sample_rate = 1/median(diff(t));

mua_choice_p = cellfun(@(x) cat(5,x{:}) < 0.05,{batch_vars(:).mua_choice_p},'uni',false);
mua_choice_p_mean = cellfun(@(x) nanmean(x,5),mua_choice_p,'uni',false);
mua_choice_p_mean = nanmean(cat(5,mua_choice_p_mean{:}),5);

fluor_choice_p = cellfun(@(x) cat(5,x{:}) < 0.05,{batch_vars(:).fluor_choice_p},'uni',false);
fluor_choice_p_mean = cellfun(@(x) nanmean(x,5),fluor_choice_p,'uni',false);
fluor_choice_p_mean = nanmean(cat(5,fluor_choice_p_mean{:}),5);

smooth_size = 5;
smWin = ones(1,smooth_size)/smooth_size;

figure('Name','Early move');

subplot(2,2,1); hold on
set(gca,'ColorOrder',copper(n_depths));
plot(t,conv2(mua_choice_p_mean(:,:,1,1)',smWin','same'),'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Fraction significant decision difference')
xlabel(['Time from stim onset (s)']);
title('MUA');
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:6),'uni',false))

subplot(2,2,2); hold on;
set(gca,'ColorOrder',[autumn(size(wf_roi,1));winter(size(wf_roi,1))]);
plot(t,conv2(fluor_choice_p_mean(:,:,1,1)',smWin','same'),'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Fraction significant decision difference')
xlabel(['Time from stim onset (s)']);
title('Fluor');
legend({wf_roi.area});

subplot(2,2,3); hold on
set(gca,'ColorOrder',copper(n_depths));
plot(t,conv2(mua_choice_p_mean(:,:,2,1)',smWin','same'),'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Fraction significant decision difference')
xlabel(['Time from move onset (s)']);
title('MUA');
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:6),'uni',false))

subplot(2,2,4); hold on;
set(gca,'ColorOrder',[autumn(size(wf_roi,1));winter(size(wf_roi,1))]);
plot(t,conv2(fluor_choice_p_mean(:,:,2,1)',smWin','same'),'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Fraction significant decision difference')
xlabel(['Time from move onset (s)']);
title('Fluor');
legend({wf_roi.area});

figure('Name','Late move');

subplot(2,2,1); hold on
set(gca,'ColorOrder',copper(n_depths));
plot(t,conv2(mua_choice_p_mean(:,:,1,2)',smWin','same'),'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Fraction significant decision difference')
xlabel(['Time from stim onset (s)']);
title('MUA');
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:6),'uni',false))

subplot(2,2,2); hold on;
set(gca,'ColorOrder',[autumn(size(wf_roi,1));winter(size(wf_roi,1))]);
plot(t,conv2(fluor_choice_p_mean(:,:,1,2)',smWin','same'),'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Fraction significant decision difference')
xlabel(['Time from stim onset (s)']);
title('Fluor');
legend({wf_roi.area});

subplot(2,2,3); hold on
set(gca,'ColorOrder',copper(n_depths));
plot(t,conv2(mua_choice_p_mean(:,:,2,2)',smWin','same'),'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Fraction significant decision difference')
xlabel(['Time from move onset (s)']);
title('MUA');
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:6),'uni',false))

subplot(2,2,4); hold on;
set(gca,'ColorOrder',[autumn(size(wf_roi,1));winter(size(wf_roi,1))]);
plot(t,conv2(fluor_choice_p_mean(:,:,2,2)',smWin','same'),'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Fraction significant decision difference')
xlabel(['Time from move onset (s)']);
title('Fluor');
legend({wf_roi.area});


%% Load both aligned fluor and MUA, plot against each other

contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';

% Get fluorescence
surround_window = [-0.5,3];
upsample_rate = 5;

framerate = 35.2;
surround_samplerate = 1/(framerate*upsample_rate);
t_fluor = surround_window(1):surround_samplerate:surround_window(2);

data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
load([data_path filesep 'roi_choiceworld']);

wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

roi_psth = {batch_vars.roi_psth};
roi_psth_mean = nanmean(cell2mat(permute(cellfun(@(x) nanmean(x,5),roi_psth,'uni',false),[1,3,4,5,2])),5);
roi_psth_mean(:,:,size(wf_roi,1)+1:end,:) = roi_psth_mean(:,:,1:size(wf_roi,1),:) - ...
    roi_psth_mean(:,:,size(wf_roi,1)+1:end,:);

% Get MUA
n_depths = 6;

raster_window = [-0.5,3];
psth_bin_size = 0.001;
t_edges = raster_window(1):psth_bin_size:raster_window(2);
t_mua = conv2(t_edges,[1,1]/2,'valid');

data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
mua_fn = [data_path filesep 'mua_choiceworld'];
load(mua_fn);

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

t_baseline = t_mua < 0;

softnorm = 20;

depth_psth = {batch_vars.depth_psth};
mua_baseline = cellfun(@(x) nanmean(x(:,t_baseline,:,1,:),2),depth_psth,'uni',false);
depth_psth_norm = cellfun(@(mua,baseline) ...
    bsxfun(@rdivide,mua,baseline + softnorm),depth_psth,mua_baseline,'uni',false);
depth_psth_smoothed = cellfun(@(x) convn(x,smWin,'same'),depth_psth_norm,'uni',false);
depth_psth_mean_t_mua = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmean(x,5),depth_psth_smoothed,'uni',false),[1,3,4,5,2])),5);

% Interpolate MUA to same time as fluor
t = t_fluor;
depth_psth_mean = permute(interp1(t_mua,permute(depth_psth_mean_t_mua,[2,1,3,4]),t),[2,1,3,4]);

% Plot all areas for one condition
plot_align = 2;
plot_condition = ismember(conditions,[1,1,-1,1],'rows');

plot_fluor = permute(roi_psth_mean(plot_condition,:,:,plot_align),[2,3,1]);
plot_fluor = bsxfun(@minus,plot_fluor,nanmean(plot_fluor(t < -0.2,:),1));
plot_fluor = bsxfun(@rdivide,plot_fluor,max(plot_fluor,[],1));

plot_mua = permute(depth_psth_mean(plot_condition,:,:,plot_align),[2,3,1]);
plot_mua = bsxfun(@minus,plot_mua,nanmean(plot_mua(t < -0.2,:),1));
plot_mua = bsxfun(@rdivide,plot_mua,max(plot_mua,[],1));

figure; 
p1 = subplot(1,3,1);hold on;
set(gca,'ColorOrder',[autumn(size(wf_roi,1))]);
plot(t,[plot_fluor(:,1:size(wf_roi,1))],'linewidth',2);
line([0,0],ylim,'color','k');

p2 = subplot(1,3,2);hold on;
set(gca,'ColorOrder',[winter(size(wf_roi,1))]);
plot(t,[plot_fluor(:,size(wf_roi,1)+1:end)],'linewidth',2);
line([0,0],ylim,'color','k');

p3 = subplot(1,3,3);hold on;
set(gca,'ColorOrder',[copper(n_depths)]);
plot(t,[plot_mua],'linewidth',2);
line([0,0],ylim,'color','k');

linkaxes([p1,p2,p3]);

% Plot
plot_success = 1;
plot_align = 2;
plot_timing = 1;
plot_t = t > -0.2 & t < 0;

plot_conditions = conditions(:,2) == -sign(plot_success)*conditions(:,3) & ...
    conditions(:,4) == plot_timing;

figure;
ax = tight_subplot(n_rois,n_depths,[0.01,0.01]);
for curr_roi = 1:n_rois
    for curr_depth = 1:n_depths
        axes(ax((curr_roi-1)*n_depths+curr_depth)); hold on;      
        set(gca,'ColorOrder',[autumn(length(contrasts));winter(length(contrasts))]);
        plot(roi_psth_mean(plot_conditions,plot_t,curr_roi,plot_align)', ...
            depth_psth_mean(plot_conditions,plot_t,curr_depth,plot_align)','linewidth',2);
        axis tight square off;
        drawnow;
    end
end

curr_roi = 10;
curr_depth = 3;
plot_t = t > -0.2 & t < 0;

figure; hold on;
set(gca,'ColorOrder',[autumn(length(contrasts));winter(length(contrasts))]);
plot(roi_psth_mean(plot_conditions,plot_t,curr_roi,plot_align)', ...
    depth_psth_mean(plot_conditions,plot_t,curr_depth,plot_align)','linewidth',2);
axis tight square;
xlabel(wf_roi(curr_roi).area);
ylabel(['Str ' num2str(curr_depth)]);

figure; 
p1 = subplot(3,1,1);hold on;
set(gca,'ColorOrder',[autumn(length(contrasts));winter(length(contrasts))]);
plot(t,roi_psth_mean(plot_conditions,:,curr_roi,plot_align)','linewidth',2);
title(wf_roi(curr_roi).area)

p2 = subplot(3,1,2);hold on;
set(gca,'ColorOrder',[autumn(length(contrasts));winter(length(contrasts))]);
plot(t,depth_psth_mean(plot_conditions,:,curr_depth,plot_align)','linewidth',2);
title(['Str ' num2str(curr_depth)]);

p3 = subplot(3,1,3);hold on;
set(gca,'ColorOrder',[autumn(length(contrasts));winter(length(contrasts))]);
plot(t, ...
    bsxfun(@rdivide,roi_psth_mean(plot_conditions,:,curr_roi,plot_align)', ...
    max(roi_psth_mean(plot_conditions,:,curr_roi,plot_align)',[],1)) - ...
    bsxfun(@rdivide,depth_psth_mean(plot_conditions,:,curr_depth,plot_align)', ...
    max(depth_psth_mean(plot_conditions,:,curr_depth,plot_align)',[],1)),'linewidth',2);
title([wf_roi(curr_roi).area ' - Str ' num2str(curr_depth)])





a = bsxfun(@rdivide,roi_psth_mean(plot_conditions,:,2,plot_align)', ...
    max(roi_psth_mean(plot_conditions,:,2,plot_align)',[],1)) - ...
    bsxfun(@rdivide,depth_psth_mean(plot_conditions,:,curr_depth,plot_align)', ...
    max(depth_psth_mean(plot_conditions,:,curr_depth,plot_align)',[],1));
b = bsxfun(@rdivide,roi_psth_mean(plot_conditions,:,10,plot_align)', ...
    max(roi_psth_mean(plot_conditions,:,10,plot_align)',[],1)) - ...
    bsxfun(@rdivide,depth_psth_mean(plot_conditions,:,curr_depth,plot_align)', ...
    max(depth_psth_mean(plot_conditions,:,curr_depth,plot_align)',[],1));
figure; hold on;
set(gca,'ColorOrder',[autumn(length(contrasts));winter(length(contrasts))]);
plot(t,a-b,'linewidth',2);




%% Load/process Peter model

fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\loglik_increase_early';
% fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\loglik_increase_late';
load(fn);

% Widefield ROs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

% Get MUA
n_depths = 6;

% Get time
raster_window = [-0.5,1];
upsample_factor = 3;
framerate = 35.2;
sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):sample_rate:raster_window(2);

% Combine and plot data
loglik_increase_emp = {batch_vars.loglik_increase_empirical};

loglik_increase_fluor = {batch_vars.loglik_increase_fluor};
loglik_increase_fluor_rel = cellfun(@(ll,emp) ...
    bsxfun(@rdivide,ll,permute(emp,[1,3,4,2])),loglik_increase_fluor,loglik_increase_emp,'uni',false);
loglik_increase_fluor_animal_mean = cellfun(@(x) nanmedian(x,4),loglik_increase_fluor,'uni',false);
loglik_increase_fluor_mean = nanmean(cat(4,loglik_increase_fluor_animal_mean{:}),4);

loglik_increase_mua = {batch_vars.loglik_increase_mua};
loglik_increase_mua_rel = cellfun(@(ll,emp) ...
    bsxfun(@rdivide,ll,permute(emp,[1,3,4,2])),loglik_increase_mua,loglik_increase_emp,'uni',false);
loglik_increase_mua_animal_mean = cellfun(@(x) nanmedian(x,4),loglik_increase_mua,'uni',false);
loglik_increase_mua_mean = nanmean(cat(4,loglik_increase_mua_animal_mean{:}),4);

figure; 

subplot(2,2,1); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot(t,loglik_increase_mua_mean(:,:,1),'linewidth',2)
ylabel('Relative log likelihood');
xlabel('Time from stim onset');
line([0,0],ylim,'color','k');

subplot(2,2,3); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot(t,loglik_increase_mua_mean(:,:,2),'linewidth',2)
ylabel('Relative log likelihood');
xlabel('Time from move onset');
line([0,0],ylim,'color','k');
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false))

subplot(2,2,2); hold on;
set(gca,'ColorOrder',[autumn(size(wf_roi,1));winter(size(wf_roi,1))]);
plot(t,loglik_increase_fluor_mean(:,:,1),'linewidth',2)
ylabel('Relative log likelihood');
xlabel('Time from stim onset');
line([0,0],ylim,'color','k');

subplot(2,2,4); hold on;
set(gca,'ColorOrder',[autumn(size(wf_roi,1));winter(size(wf_roi,1))]);
plot(t,loglik_increase_fluor_mean(:,:,2),'linewidth',2)
ylabel('Relative log likelihood');
xlabel('Time from move onset');
line([0,0],ylim,'color','k');
legend({wf_roi.area});


%% Load Peter mean time model

% fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\c_cn_modeling_ll';
% fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\c_cn_modeling_ll_hemisphere';
fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\c_cn_modeling_ll_earlymove';

load(fn);

% Widefield ROs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

% Get MUA
n_depths = 6;

% Get time
raster_window = [-0.5,1];
upsample_factor = 3;
framerate = 35.2;
sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):sample_rate:raster_window(2);

% Combine and plot data
loglik_increase_fluor = {batch_vars.loglik_increase_fluor};
loglik_increase_fluor_animal_mean = cellfun(@(x) nanmedian(x,3),loglik_increase_fluor,'uni',false);
loglik_increase_fluor_mean = nanmean(cat(3,loglik_increase_fluor_animal_mean{:}),3);

loglik_increase_mua = {batch_vars.loglik_increase_mua};
loglik_increase_mua_animal_mean = cellfun(@(x) nanmedian(x,3),loglik_increase_mua,'uni',false);
loglik_increase_mua_mean = nanmean(cat(3,loglik_increase_mua_animal_mean{:}),3);

figure;
subplot(2,1,1); hold on;
plot_col = [autumn(n_rois/2);winter(n_rois/2)];
scatter(1:n_rois,loglik_increase_fluor_mean(:,1),100,plot_col,'filled','MarkerEdgeColor','k','linewidth',2);
scatter(1:n_rois,loglik_increase_fluor_mean(:,2),100,plot_col,'filled','MarkerEdgeColor',[0.5,0.5,0.5],'linewidth',2);
ylabel('Relative log likelihood (bpt)');
line(xlim,[0,0],'color','k');
set(gca,'XTick',1:n_rois,'XTickLabel',{wf_roi.area});
legend({'Stim-aligned','Move-aligned'})

subplot(2,1,2); hold on;
plot_col = copper(n_depths);
scatter(1:n_depths,loglik_increase_mua_mean(:,1),100,plot_col,'filled','MarkerEdgeColor','k','linewidth',2);
scatter(1:n_depths,loglik_increase_mua_mean(:,2),100,plot_col,'filled','MarkerEdgeColor',[0.5,0.5,0.5],'linewidth',2);
ylabel('Relative log likelihood (bpt)');
line(xlim,[0,0],'color','k');
set(gca,'XTick',1:n_depths,'XTickLabel',cellfun(@(x) ...
    ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false));

% % (for plotting L/R params)
% fluor_params = {batch_vars.fluor_params};
% fluor_weight = cellfun(@(x) [abs(x(:,5,:,:)),abs(x(:,6,:,:)) ...
%     .*sign(x(:,5,:,:)).*sign(x(:,6,:,:))],fluor_params,'uni',false);
% fluor_weight_animal_mean = cellfun(@(x) nanmedian(x,4),fluor_weight,'uni',false);
% fluor_weight_mean = nanmean(cat(4,fluor_weight_animal_mean{:}),4);
% 
% figure;
% subplot(2,1,1);
% plot(squeeze(fluor_weight_mean(:,:,1)),'linewidth',2)
% legend({'L hemi','R hemi'})
% ylabel('Weight')
% set(gca,'XTick',1:n_rois/2,'XTickLabel',{wf_roi(:,1).area});
% title('Stim-aligned');
% 
% subplot(2,1,2);
% plot(squeeze(fluor_weight_mean(:,:,2)),'linewidth',2)
% legend({'L hemi','R hemi'})
% ylabel('Weight')
% set(gca,'XTick',1:n_rois/2,'XTickLabel',{wf_roi(:,1).area});
% title('Move-aligned')

%% Trial-by-trial activity time-averaged
error('Don''t use this - use the distribution plots below');

% Load trial activity
% [area, trial ID, alignment, day]
fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\trial_activity';
load(fn);

% Widefield ROs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

% Get MUA and normalize
n_depths = 6;
softnorm = 20;
mua_trial_act_norm = {batch_vars.mua_trial_act};
for curr_animal = 1:length(batch_vars)
    for curr_day = 1:size(mua_trial_act_norm{curr_animal},4);
        for curr_depth = 1:n_depths
            curr_median = nanmedian(vertcat(mua_trial_act_norm{curr_animal}{curr_depth,:,:,curr_day}));
            curr_std = nanstd(vertcat(mua_trial_act_norm{curr_animal}{curr_depth,:,:,curr_day}));
            mua_trial_act_norm{curr_animal}(curr_depth,:,:,curr_day) = ...
                cellfun(@(x) x/(curr_median + softnorm), ...
                mua_trial_act_norm{curr_animal}(curr_depth,:,:,curr_day),'uni',false);
        end
    end
end

% Get trial conditions
% [contrast,side,choice,timing]
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];

conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);



% Plot distribution of all trials 
plot_act = 'fluor';
plot_area = 1;
plot_conditions = conditions(:,4) == 1;
plot_align = 1;

switch plot_act
    case 'fluor'
        use_data = {batch_vars.fluor_trial_act};
        area_labels = {wf_roi.area};
    case 'mua'
        use_data =  mua_trial_act_norm;
        area_labels = cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false);
end

plot_contrasts = [unique(conditions(plot_conditions,1))*sides(1); ...
    unique(conditions(plot_conditions,1))*sides(2)];
[~,sort_idx] = sort(plot_contrasts);

curr_trial_act = cellfun(@(x) ...
    squeeze(x(plot_area,plot_conditions,plot_align,:)), ...
    use_data,'uni',false);
curr_trial_act_animalcat = horzcat(curr_trial_act{:});
curr_trial_act_allcat = reshape(arrayfun(@(y) ...
    vertcat(curr_trial_act_animalcat{y,:}), ...
    1:size(curr_trial_act_animalcat,1),'uni',false),[],2);
% (combine zeros)
zero_contrasts = find(plot_contrasts == 0);
zero_trials = {vertcat(curr_trial_act_allcat{zero_contrasts,1}), ...
    vertcat(curr_trial_act_allcat{zero_contrasts,2})};
curr_trial_act_allcat(zero_contrasts(1),:) = zero_trials;
curr_trial_act_allcat(zero_contrasts(2:end),:) = [];
curr_plot_contrasts = plot_contrasts(setdiff(1:length(plot_contrasts),zero_contrasts(2:end)));
[~,curr_sort_idx] = sort(curr_plot_contrasts);

curr_trial_act_allcat_median = cellfun(@nanmedian,curr_trial_act_allcat);

figure; hold on;
p1 = distributionPlot(curr_trial_act_allcat(curr_sort_idx,1),'distWidth',0.5, ...
    'xValues',(1:11)-0.25,'histOri','left','showMM',0,'color','r');
p2 = distributionPlot(curr_trial_act_allcat(curr_sort_idx,2),'distWidth',0.5, ...
    'xValues',(1:11)+0.25,'histOri','right','showMM',0,'color','b');
plot((1:11),curr_trial_act_allcat_median(curr_sort_idx,1),'color',[0.9,0.6,0.6],'linewidth',2)
plot((1:11),curr_trial_act_allcat_median(curr_sort_idx,2),'color',[0.6,0.6,0.9],'linewidth',2)

set(gca,'XTick',1:11,'XTickLabel',cellfun(@num2str,num2cell(curr_plot_contrasts(curr_sort_idx)),'uni',false));
xlabel('Contrast*Side');
ylabel(area_labels{plot_area});

plot_choices = unique(conditions(plot_conditions,3));
plot_timing = unique(conditions(plot_conditions,4));
condition_labels{1} = ['Choice ' num2str(plot_choices(1)), ', Timing ' num2str(plot_timing(1))];
condition_labels{2} = ['Choice ' num2str(plot_choices(end)), ', Timing ' num2str(plot_timing(end))];
legend([p1{1}(1),p2{1}(1)],condition_labels);

align_labels = {'Stim','Move','Beep'};
title([align_labels{plot_align} '-aligned']);






% Plot averaging across days then animals
curr_align = 2;

figure('Name','Within days, Within animals');
for curr_roi = 1:n_rois;    
    curr_trial_act = cellfun(@(x) ...
        squeeze(x(curr_roi,plot_conditions,curr_align,:)), ...
        {batch_vars.fluor_trial_act},'uni',false);
    
    curr_trial_act_animalmean = cellfun(@(x) cellfun(@nanmean,x),curr_trial_act,'uni',false);
    curr_trial_act_mean = reshape(nanmean(cell2mat(cellfun(@(x) nanmedian(x,2), ...
        curr_trial_act_animalmean,'uni',false)),2),[],2);
   
    subplot(3,8,curr_roi); hold on;   
    % (to plot just mean L and R)
    plot(plot_contrasts(sort_idx),curr_trial_act_mean(sort_idx,:));
    xlabel('Condition');
    ylabel(wf_roi(curr_roi).area);   
end
legend({'Move L','Move R'});

for curr_depth = 1:n_depths;    
     curr_trial_act = cellfun(@(x) ...
        squeeze(x(curr_depth,plot_conditions,curr_align,:)), ...
        mua_trial_act_norm,'uni',false);
    
    curr_trial_act_animalmean = cellfun(@(x) cellfun(@nanmean,x),curr_trial_act,'uni',false);
    curr_trial_act_mean = reshape(nanmean(cell2mat(cellfun(@(x) nanmedian(x,2), ...
        curr_trial_act_animalmean,'uni',false)),2),[],2);
   
    subplot(3,8,n_rois+curr_depth); hold on;
    % (to plot just mean L and R)
    plot(plot_contrasts(sort_idx),curr_trial_act_mean(sort_idx,:));
    xlabel('Condition');
    ylabel(['Str ' num2str(curr_depth)]);
end
legend({'Move L','Move R'});


% Plot combining trials within animals
figure('Name','Across days, within animals');
for curr_roi = 1:n_rois;    
    curr_trial_act = cellfun(@(x) ...
        squeeze(x(curr_roi,plot_conditions,curr_align,:)), ...
        {batch_vars.fluor_trial_act},'uni',false);
    curr_trial_act_cat = cellfun(@(x) reshape(arrayfun(@(y) vertcat(x{y,:}),1:size(x,1),'uni',false),[],2),curr_trial_act,'uni',false);
    curr_trial_act_mean = nanmean(cell2mat(permute(cellfun(@(x) cellfun(@nanmean,x),curr_trial_act_cat,'uni',false),[1,3,2])),3);
    curr_trial_act_std = nanmean(cell2mat(permute(cellfun(@(x) cellfun(@nanstd,x),curr_trial_act_cat,'uni',false),[1,3,2])),3);
    
    curr_trial_act_cat_choicecat = cellfun(@(x) arrayfun(@(y) ...
        vertcat(x{y,:}),1:size(x,1),'uni',false),curr_trial_act_cat,'uni',false);
    curr_trial_act_choicecat_mean = nanmean(cell2mat(permute(cellfun(@(x) ...
        cellfun(@nanmean,x),curr_trial_act_cat_choicecat,'uni',false),[1,3,2])),3)';

    subplot(3,8,curr_roi); hold on;
    % (to plot just mean L and R)
    plot(plot_contrasts(sort_idx),curr_trial_act_mean(sort_idx,:));
    % (to subtract mean of all trials combined)
%     plot(plot_contrasts(sort_idx),bsxfun(@minus, ...
%         curr_trial_act_mean(sort_idx,:),curr_trial_act_choicecat_mean(sort_idx,:)));
    % (to plot std)
%     errorbar(repmat(plot_contrasts(sort_idx)',1,2),...
%         bsxfun(@minus,curr_trial_act_mean(sort_idx,:), ...
%         curr_trial_act_choicecat_mean(sort_idx,:)), ...
%         curr_trial_act_std(sort_idx,:));

    xlabel('Condition');
    ylabel(wf_roi(curr_roi).area);   
end
legend({'Move L','Move R'});

for curr_depth = 1:n_depths;    
    curr_trial_act = cellfun(@(x) ...
        squeeze(x(curr_depth,plot_conditions,curr_align,:)), ...
        mua_trial_act_norm,'uni',false);
    curr_trial_act_cat = cellfun(@(x) arrayfun(@(y) vertcat(x{y,:}),1:size(x,1),'uni',false),curr_trial_act,'uni',false);
    curr_trial_act_mean = nanmean(cell2mat(permute(cellfun(@(x) reshape(cellfun(@nanmean,x),[],2),curr_trial_act_cat,'uni',false),[1,3,2])),3);
    
    subplot(3,8,n_rois+curr_depth); hold on;
    % (to plot just mean L and R)
        plot(plot_contrasts(sort_idx),curr_trial_act_mean(sort_idx,:));
    % (to subtract mean of all trials combined)
%     plot(plot_contrasts(sort_idx),bsxfun(@minus, ...
%         curr_trial_act_mean(sort_idx,:),curr_trial_act_choicecat_mean(sort_idx,:)));
    xlabel('Condition');
    ylabel(['Str ' num2str(curr_depth)]);
end
legend({'Move L','Move R'});


% Plot combining all trials across animals
figure('Name','Across days, across animals');
for curr_roi = 1:n_rois;
    curr_trial_act = cellfun(@(x) ...
        squeeze(x(curr_roi,plot_conditions,curr_align,:)), ...
        {batch_vars.fluor_trial_act},'uni',false);
    
    curr_trial_act_animalcat = horzcat(curr_trial_act{:});
    curr_trial_act_allcat = reshape(arrayfun(@(y) ...
        vertcat(curr_trial_act_animalcat{y,:}), ...
        1:size(curr_trial_act_animalcat,1),'uni',false),[],2);
    curr_trial_act_mean = cellfun(@nanmedian,curr_trial_act_allcat);
    
    subplot(3,8,curr_roi); hold on;
    plot(plot_contrasts(sort_idx),curr_trial_act_mean(sort_idx,:))
    xlabel('Condition');
    ylabel(wf_roi(curr_roi).area);
end
legend({'Move L','Move R'});

for curr_depth = 1:n_depths;
    curr_trial_act = cellfun(@(x) ...
        squeeze(x(curr_depth,plot_conditions,curr_align,:)), ...
        mua_trial_act_norm,'uni',false);
    
    curr_trial_act_animalcat = horzcat(curr_trial_act{:});
    curr_trial_act_allcat = reshape(arrayfun(@(y) ...
        vertcat(curr_trial_act_animalcat{y,:}), ...
        1:size(curr_trial_act_animalcat,1),'uni',false),[],2);
    curr_trial_act_mean = cellfun(@nanmedian,curr_trial_act_allcat);
    
    subplot(3,8,n_rois+curr_depth); hold on;
    plot(plot_contrasts(sort_idx),curr_trial_act_mean(sort_idx,:))
    xlabel('Condition');
    ylabel(['Str ' num2str(curr_depth)]);
end
legend({'Move L','Move R'});


%% Ctx (ROI) -> str regression, also tr-tr regression

% Load data
data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
mua_fn = [data_path filesep 'roi-mua_choiceworld_predicted'];
load(mua_fn);

% Widefield ROs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

n_depths = 6;

contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';

frame_rate = 35.2;
upsample_factor = 2;
sample_rate = frame_rate*upsample_factor;

interval_surround = [-0.5,3];
t = interval_surround(1):1/sample_rate:interval_surround(2);

% Smooth the real psth to match the inherently smoothed predicted
gw = gausswin(9,1)';
smWin = gw./sum(gw);

%%% Params from trial activity within day

% [depth,time,align,timing,param]
activity_model_params = {batch_vars.activity_model_params};
activity_model_params_mean = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmean(x,6),activity_model_params,'uni',false),[1,3,4,5,6,2])),6);
activity_model_params_mean_smoothed = convn(activity_model_params_mean,smWin,'same');

predicted_activity_model_params = {batch_vars.predicted_activity_model_params};
predicted_activity_model_params_mean = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmean(x,6),predicted_activity_model_params,'uni',false),[1,3,4,5,6,2])),6);

plot_param = 4;

spacing = max(abs(reshape(activity_model_params_mean_smoothed(:,:,:,:,plot_param),[],1)))*1.5;

figure; 
s1 = subplot(1,2,1);hold on;
p1 = AP_stackplot(activity_model_params_mean_smoothed(:,:,1,1,plot_param)',t,spacing,false,'k');
p2 = AP_stackplot(predicted_activity_model_params_mean(:,:,1,1,plot_param)',t,spacing,false,'r');
p3 = AP_stackplot(activity_model_params_mean_smoothed(:,:,1,2,plot_param)',t,spacing,false,'b');
p4 = AP_stackplot(predicted_activity_model_params_mean(:,:,1,2,plot_param)',t,spacing,false,'m');
line([0,0],ylim,'color','k');
line([0.5,0.5],ylim,'color','k');
xlabel('Time from stim')
ylabel(['Param ' num2str(plot_param)])
legend([p1(1),p2(1),p3(1),p4(1)],{'Early real','Early predicted','Late real','Late predicted'});

s2 = subplot(1,2,2);hold on;
p1 = AP_stackplot(activity_model_params_mean_smoothed(:,:,2,1,plot_param)',t,spacing,false,'k');
p2 = AP_stackplot(predicted_activity_model_params_mean(:,:,2,1,plot_param)',t,spacing,false,'r');
p3 = AP_stackplot(activity_model_params_mean_smoothed(:,:,2,2,plot_param)',t,spacing,false,'b');
p4 = AP_stackplot(predicted_activity_model_params_mean(:,:,2,2,plot_param)',t,spacing,false,'m');
line([0,0],ylim,'color','k');
xlabel('Time from move')
ylabel(['Param ' num2str(plot_param)])
legend([p1(1),p2(1),p3(1),p4(1)],{'Early real','Early predicted','Late real','Late predicted'});

linkaxes([s1,s2]);

% Plot unsmoothed
plot_param = 4;

spacing = max(abs(reshape(activity_model_params_mean(:,:,:,:,plot_param),[],1)))*1.5;

figure; 
s1 = subplot(1,2,1);hold on;
p1 = AP_stackplot(activity_model_params_mean(:,:,1,1,plot_param)',t,spacing,false,'k');
p2 = AP_stackplot(activity_model_params_mean(:,:,1,2,plot_param)',t,spacing,false,'b');
line([0,0],ylim,'color','k');
line([0.5,0.5],ylim,'color','k');
xlabel('Time from stim')
ylabel(['Param ' num2str(plot_param)])
legend([p1(1),p2(1)],{'Early real','Late real'});

s2 = subplot(1,2,2);hold on;
p1 = AP_stackplot(activity_model_params_mean(:,:,2,1,plot_param)',t,spacing,false,'k');
p2 = AP_stackplot(activity_model_params_mean(:,:,2,2,plot_param)',t,spacing,false,'b');
line([0,0],ylim,'color','k');
xlabel('Time from move')
ylabel(['Param ' num2str(plot_param)])
legend([p1(1),p2(1)],{'Early real','Late real'});

linkaxes([s1,s2]);


% Kernel for ROIs by depth
% [ROI, time, depth, alignment, timing];
roi_depth_kernel = {batch_vars.roi_depth_kernel};
roi_depth_kernel_mean = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmean(x,6),roi_depth_kernel,'uni',false),[1,3,4,5,6,2])),6);

%% Simulate data
% (this was initially to look at L/L-R differences, not sure if that made
% sense but at least it's useful for testing stuff)

contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];

n_trials = 10000;

% Make sides and contrasts
trial_sides = sides(randi(length(sides),n_trials,1));
trial_contrasts = contrasts(randi(length(contrasts),n_trials,1));

trial_contrasts_sides = zeros(n_trials,2);
trial_contrasts_sides(trial_sides == -1,1) = trial_contrasts(trial_sides == -1);
trial_contrasts_sides(trial_sides == 1,2) = trial_contrasts(trial_sides == 1);

% Generate activity
private_act_noise_coeff = 0.2;
shared_act_noise_coeff = 1;
choice_noise_coeff = 0.5;

shared_noise = shared_act_noise_coeff*randn(n_trials,1);
act_l = trial_contrasts_sides(:,2).^0.3 + private_act_noise_coeff*randn(n_trials,1);
act_r = trial_contrasts_sides(:,1).^0.3 + private_act_noise_coeff*randn(n_trials,1);

% % (unilateral: activity based on trial choice)
trial_choice = -sign((act_l - act_r) + choice_noise_coeff*randn(n_trials,1));

% (bilateral: trial choice based on L-R)
% act_choice_add = 0.5;
% act_choice_noise_coeff = 0.5;
% trial_choice = sign(-trial_contrasts_sides(:,1).^0.3+trial_contrasts_sides(:,2).^0.3 + choice_noise_coeff*randn(n_trials,1));
% act_l = act_l + act_choice_add*(trial_choice == -1) + act_choice_noise_coeff*randn(n_trials,1);
% act_r = act_r + act_choice_add*(trial_choice == 1) + act_choice_noise_coeff*randn(n_trials,1);

[psychometric,conditions] = grpstats(trial_choice == -1,trial_sides.*trial_contrasts,{'mean','gname'});
conditions = cellfun(@str2num,conditions);
figure;plot(conditions,psychometric,'color','k','linewidth',3);

% Regress L activity from task parameters
contrast_exp = 1;

regression_unilateral = [ones(n_trials,1), ...
    trial_contrasts_sides.^contrast_exp, ...
    trial_choice]\act_l;

regression_bilateral = [ones(n_trials,1), ...
    trial_contrasts_sides.^contrast_exp, ...
    trial_choice]\(act_l-act_r);

regression_choice_ratio = regression_bilateral(4)/regression_unilateral(4);
disp(['Regression ratio: ' num2str(regression_choice_ratio)]);       


% Fit Peter's model
cvfold = 10;

D = struct;
D.stimulus = trial_contrasts_sides;
D.response = ((trial_choice+1)/2)+1;
D.repeatNum = ones(n_trials,1);

% Fit all stim
use_model = 'AP_test_stim';
g_stim_all = GLM(D).setModel(use_model).fit;
behavParameterFit = g_stim_all.parameterFits;

% Fit without neural CV
use_model = 'AP_test_stim';

g_stim = GLM(D).setModel(use_model).fitCV(cvfold);
pL = g_stim.p_hat(:,1);    pR = g_stim.p_hat(:,2);
likelihood = pL.*(g_stim.data.response==1) + pR.*(g_stim.data.response==2);

stim_params = nanmean(g_stim.parameterFits,1);
loglik_bpt_stim = nanmean(log2(likelihood(likelihood ~= 0)));

% Fit with neural
use_model = 'AP_test_neur_stim';

D.neur = act_l; 
clear g_fluor
g_act = GLM(D).setModel(use_model).fitCV(cvfold);
pL = g_act.p_hat(:,1);
pR = g_act.p_hat(:,2);
likelihood = pL.*(g_act.data.response==1) + pR.*(g_act.data.response==2);
loglik_bpt_act_unilateral = nanmean(log2(likelihood(likelihood ~= 0)));

D.neur = act_l-act_r;
clear g_fluor
g_act = GLM(D).setModel(use_model).fitCV(cvfold);
pL = g_act.p_hat(:,1);
pR = g_act.p_hat(:,2);
likelihood = pL.*(g_act.data.response==1) + pR.*(g_act.data.response==2);
loglik_bpt_act_bilateral = nanmean(log2(likelihood(likelihood ~= 0)));

loglik_increase_unilateral = loglik_bpt_act_unilateral - loglik_bpt_stim;
loglik_increase_bilateral = loglik_bpt_act_bilateral - loglik_bpt_stim;

model_choice_ratio = loglik_increase_bilateral/loglik_increase_unilateral;
disp(['Model ratio: ' num2str(model_choice_ratio)]);       
        
%% Logistic regression from V's to choice

% Load
fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\logistic_regression_pV';
load(fn);

% Get time
framerate = 35.2;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):sample_rate:raster_window(2);

% Take the mean progressively and overwrite - file's too big
for curr_animal = 1:length(batch_vars)
    batch_vars(curr_animal).pV = nanmean(batch_vars(curr_animal).pV,5);
end

pV = nanmean(cat(5,batch_vars.pV),5);

AP_image_scroll(pV,t);
axis image;
caxis([-5e-11,5e-11]);
colormap(colormap_BlueWhiteRed);

%% Logistic regression: concatenate trials within animal, mean time

analysis_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
analysis_fn = ['activity_choice_logistic_regression_meantime'];
load([analysis_path filesep analysis_fn]);

% Widefield ROs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

% Get MUA
n_depths = 6;

% Get time
raster_window = [-0.5,1];
upsample_factor = 3;
framerate = 35.2;
sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):sample_rate:raster_window(2);

% Plot average data
loglik_increase_fluor = nanmean(batch_vars.loglik_increase_fluor,ndims(batch_vars.loglik_increase_fluor));
loglik_increase_mua = nanmean(batch_vars.loglik_increase_mua,ndims(batch_vars.loglik_increase_mua));

figure;
subplot(2,1,1);
plot(loglik_increase_fluor,'linewidth',2);
line(xlim,[0,0],'color','k');
set(gca,'XTick',1:n_rois,'XTickLabel',{wf_roi.area});
ylabel('Relative loglikelihood (bpt)');
legend({'Stim-aligned','Move-aligned'});

subplot(2,1,2);
plot(loglik_increase_mua,'linewidth',2);
line(xlim,[0,0],'color','k');
set(gca,'XTick',1:n_depths);
ylabel('Relative loglikelihood (bpt)');
legend({'Stim-aligned','Move-aligned'});


%% Load logistic regression on day-concatenated activity

% Load data
data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
data_fn = ['activity_sessioncat_logistic_regression_earlymove_kernel-str'];
% data_fn = ['activity_sessioncat_logistic_regression_earlymove'];
load([data_path filesep data_fn])

% Get time
framerate = 35.2;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):sample_rate:raster_window(2);

% Get widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

n_depths = size(loglik_increase_mua,2);

% Plot mean
yrange = [min(reshape(nanmedian(loglik_increase_fluor,4),[],1)), ...
    max(reshape(nanmedian(loglik_increase_fluor,4),[],1))];

figure; 
p1 = subplot(2,2,1); hold on;
set(gca,'ColorOrder',[autumn(n_rois/2);winter(n_rois/2)]);
plot(t,nanmean(loglik_increase_fluor(:,:,1,:),4))
line([0,0],yrange,'color','k');
xlabel('Time from stim')
ylabel('Relative loglikelihood (bpt)');

p2 = subplot(2,2,2); hold on;
set(gca,'ColorOrder',[autumn(n_rois/2);winter(n_rois/2)]);
plot(t,nanmean(loglik_increase_fluor(:,:,2,:),4))
line([0,0],yrange,'color','k');
xlabel('Time from move')
ylabel('Relative loglikelihood (bpt)');

p3 = subplot(2,2,3); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot(t,nanmean(loglik_increase_mua(:,:,1,:),4))
line([0,0],yrange,'color','k');
xlabel('Time from stim')
ylabel('Relative loglikelihood (bpt)');

p4 = subplot(2,2,4); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot(t,nanmean(loglik_increase_mua(:,:,2,:),4))
line([0,0],yrange,'color','k');
xlabel('Time from move')
ylabel('Relative loglikelihood (bpt)');

linkaxes([p1,p2,p3,p4])
ylim(yrange);

% Plot move-aligned all together
figure; hold on;
set(gca,'ColorOrder',[autumn(n_rois/2);winter(n_rois/2);copper(n_depths)]);
plot(t,nanmean(loglik_increase_fluor(:,:,2,:),4))
plot(t,nanmean(loglik_increase_mua(:,:,2,:),4))

% Stack plot across animals for each region
figure; hold on;
n_animals = size(loglik_increase_mua,4);
col = copper(n_animals);
p = gobjects(n_animals,1);
for curr_animal = 1:n_animals
    curr_p = AP_stackplot(squeeze(loglik_increase_mua(:,:,2,curr_animal)),t,0.15,false,col(curr_animal,:));
    p(curr_animal) = curr_p(1);
end
animals = {'AP024','AP025','AP027','AP028','AP029'};
legend(p,animals);

%% Load logistic regression (all within-modal data simultaneously)

% Load data
data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
data_fn = ['activity_sessioncat_logistic_regression_neucombined_earlymove_4str'];
load([data_path filesep data_fn])

% Get time
framerate = 35.2;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):sample_rate:raster_window(2);

% Get widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

n_depths = size(loglik_increase_mua,2);

% Get mean data
loglik_increase_fluor_mean = nanmean(loglik_increase_fluor,3);
loglik_increase_mua_mean = nanmean(loglik_increase_mua,3);

loglik_params_fluor_mean = nanmean(loglik_params_fluor,4);
loglik_params_mua_mean = nanmean(loglik_params_mua,4);

% Plot mean
figure; 
subplot(3,2,1); hold on;
plot(t,loglik_increase_fluor_mean(:,1),'color',[0,0.7,0],'linewidth',2);
plot(t,loglik_increase_mua_mean(:,1),'k','linewidth',2);
line([0,0],ylim,'color','k');
xlabel('Time from stim')
ylabel('Relative loglikelihood (bpt)');

subplot(3,2,2); hold on;
plot(t,loglik_increase_fluor_mean(:,2),'color',[0,0.7,0],'linewidth',2);
plot(t,loglik_increase_mua_mean(:,2),'k','linewidth',2);
line([0,0],ylim,'color','k');
xlabel('Time from move')
ylabel('Relative loglikelihood (bpt)');

subplot(3,2,3); hold on;
set(gca,'ColorOrder',[autumn(n_rois/2);winter(n_rois/2)]);
plot(t,loglik_params_fluor_mean(:,:,1));
line([0,0],ylim,'color','k');
xlabel('Time from stim')
ylabel('Weight');

subplot(3,2,4); hold on;
set(gca,'ColorOrder',[autumn(n_rois/2);winter(n_rois/2)]);
plot(t,loglik_params_fluor_mean(:,:,1));
line([0,0],ylim,'color','k');
xlabel('Time from stim')
ylabel('Weight');

subplot(3,2,5); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot(t,loglik_params_mua_mean(:,:,1));
line([0,0],ylim,'color','k');
xlabel('Time from stim')
ylabel('Weight');

subplot(3,2,6); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot(t,loglik_params_mua_mean(:,:,2));
line([0,0],ylim,'color','k');
xlabel('Time from move')
ylabel('Weight');

%% Plot mean day-concatenated activity

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = 'all_trial_activity_df_kernel-str_earlymove.mat';
% data_fn = 'all_trial_activity_df_earlymove.mat';

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35.2;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):sample_rate:raster_window(2);

% % Load use experiments and cut out bad ones
% bhv_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\bhv_processing\use_experiments';
% load(bhv_fn);
% D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
% fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
% mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
% 
% use_animals = cellfun(@(x) ~isempty(x),D_all);
% D_all = D_all(use_animals);
% fluor_all = fluor_all(use_animals);
% mua_all = mua_all(use_animals);

% Get widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

% MUA depths
n_depths = size(mua_all{1}{1},3);

% Set up conditions
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);

fluor_psth = nan(n_conditions,length(t),n_rois,2,length(D_all));
mua_psth = nan(n_conditions,length(t),n_depths,2,length(D_all));
for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    fluor_cat_hemidiff = cat(3,fluor_cat(:,:,1:n_rois/2,:),fluor_cat(:,:,1:n_rois/2,:) - ...
        fluor_cat(:,:,n_rois/2+1:end,:));
        
    fluor_cat_norm = bsxfun(@rdivide,fluor_cat,permute(std(reshape(permute(fluor_cat,[1,2,4,3]),[],n_rois),[],1),[1,3,2,4]));
    fluor_cat_hemidiff_norm = bsxfun(@rdivide,fluor_cat_hemidiff,permute(std(abs(reshape(permute(fluor_cat_hemidiff,[1,2,4,3]),[],n_rois)),[],1),[1,3,2,4]));     
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    smooth_size = 8;
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    
    mua_cat_raw_smoothed = convn(mua_cat_raw,smWin,'same');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,2),'uni',false));
    
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x,[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,2),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
          
    % Concatenate behavioural data
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
    D.response = cell2mat(cellfun(@(x) x.response,D_all{curr_animal},'uni',false));
    D.repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all{curr_animal},'uni',false));
    
    D.day = trial_day;
    
    % Get trial ID   
    trial_contrast = max(D.stimulus,[],2);
    [~,side_idx] = max(D.stimulus > 0,[],2);
    trial_side = (side_idx-1.5)*2;
    trial_choice = -(D.response-1.5)*2;
    
    trial_conditions = ...
        [trial_contrast, trial_side, ...
        trial_choice, ones(size(trial_day))];
    [~,trial_id] = ismember(trial_conditions,conditions,'rows');
    
    % Save aligned average by trial ID
    for curr_roi = 1:n_rois
        fluor_psth(unique(trial_id),:,curr_roi,1,curr_animal) = grpstats(fluor_cat_hemidiff_norm(:,:,curr_roi,1),trial_id);
        fluor_psth(unique(trial_id),:,curr_roi,2,curr_animal) = grpstats(fluor_cat_hemidiff_norm(:,:,curr_roi,2),trial_id);
    end    
    
    for curr_depth = 1:n_depths
        mua_psth(unique(trial_id),:,curr_depth,1,curr_animal) = grpstats(mua_cat_norm(:,:,curr_depth,1),trial_id);
        mua_psth(unique(trial_id),:,curr_depth,2,curr_animal) = grpstats(mua_cat_norm(:,:,curr_depth,2),trial_id);
    end    
 
end

fluor_psth_mean = nanmean(fluor_psth,5);
mua_psth_mean = nanmean(mua_psth,5);

% Set up plot conditions
plot_timing = 1; % always same, because slots filled in above

% (to plot only correct/incorrect)
plot_success = 1;
plot_conditions = find(conditions(:,1) ~= 0 & ...
    conditions(:,2) == -plot_success*conditions(:,3) & ...
    conditions(:,4) == plot_timing);
plot_color = colormap_BlueWhiteRed(5);
plot_color = [plot_color(end-4:end,:);plot_color(5:-1:1,:)];

% % (to plot one movement direction)
% plot_move = 1;
% plot_conditions = find(conditions(:,1) ~= 0 & ...
%     conditions(:,3) == plot_move & ...
%     conditions(:,4) == plot_timing);
% plot_color = colormap_BlueWhiteRed(6);
% plot_color = [plot_color(5:-1:1,:);plot_color(end-5:end,:)];

% % (to plot one side)
% plot_side = 1;
% plot_conditions = find(conditions(:,1) < Inf & ...
%     conditions(:,2) == plot_side & ...
%     conditions(:,4) == plot_timing);
% n_plot_contrasts = length(unique(conditions(plot_conditions,1)));
% plot_color = colormap_BlueWhiteRed(n_plot_contrasts);
% plot_color = [plot_color(n_plot_contrasts:-1:1,:);plot_color(end-n_plot_contrasts+1:end,:)];

% % (to plot by movement direction)
% plot_conditions = find(conditions(:,1) == 0.25 & ...
%     conditions(:,2) == -1 & ...
%     conditions(:,4) == plot_timing);
% plot_color = zeros(sum(plot_conditions),3);
% plot_color(conditions(plot_conditions,3) == -1,3) = 1;
% plot_color(conditions(plot_conditions,3) == 1,2) = 1;


% Stack MUA plots 
figure; 
p1 = subplot(1,2,1); hold on;
for curr_condition_idx = 1:length(plot_conditions);
    curr_condition = plot_conditions(curr_condition_idx);
    AP_stackplot(squeeze(mua_psth_mean(curr_condition,:,:,1)),t,1.5,false,plot_color(curr_condition_idx,:));
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');

p2 = subplot(1,2,2); hold on;
for curr_condition_idx = 1:length(plot_conditions);
    curr_condition = plot_conditions(curr_condition_idx);
    AP_stackplot(squeeze(mua_psth_mean(curr_condition,:,:,2)),t,1.5,false,plot_color(curr_condition_idx,:));
end
line([0,0],ylim,'color','k');
xlabel('Time from move');

linkaxes([p1,p2]);

% Stack fluor plots 
figure; 
p1 = subplot(1,2,1); hold on;
for curr_condition_idx = 1:length(plot_conditions);
    curr_condition = plot_conditions(curr_condition_idx);
    AP_stackplot(squeeze(fluor_psth_mean(curr_condition,:,:,1)),t,5,false,plot_color(curr_condition_idx,:),{wf_roi.area});
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');

p2 = subplot(1,2,2); hold on;
for curr_condition_idx = 1:length(plot_conditions);
    curr_condition = plot_conditions(curr_condition_idx);
    AP_stackplot(squeeze(fluor_psth_mean(curr_condition,:,:,2)),t,5,false,plot_color(curr_condition_idx,:),{wf_roi.area});
end
line([0,0],ylim,'color','k');
xlabel('Time from move');

linkaxes([p1,p2]);


%% Distribution plots of time-averaged day-concatenated activity

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = 'all_trial_activity_df_4str_earlymove.mat';
load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35.2;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):sample_rate:raster_window(2);

% % Load use experiments and cut out bad ones
% bhv_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\bhv_processing\use_experiments';
% load(bhv_fn);
% D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
% fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
% mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
% 
% use_animals = cellfun(@(x) ~isempty(x),D_all);
% D_all = D_all(use_animals);
% fluor_all = fluor_all(use_animals);
% mua_all = mua_all(use_animals);

% Get widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

% MUA depths
n_depths = size(mua_all{1}{1},3);

% Set up conditions
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);

sidecontrasts = unique(bsxfun(@times,contrasts',sides));

% Plot distribution of all trials 
plot_act = 'mua';
plot_area = 3;
plot_align = 2;
% plot_t = t >= 0.05 & t <= 0.15;
plot_t = t >= -0.05 & t <= 0;

figure;
trial_act_timeavg_split_all = cell(length(sidecontrasts),2);
for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    fluor_cat_hemidiff = cat(3,fluor_cat(:,:,1:n_rois/2,:),fluor_cat(:,:,1:n_rois/2,:) - ...
        fluor_cat(:,:,n_rois/2+1:end,:));
        
    fluor_cat_norm = bsxfun(@rdivide,fluor_cat,permute(std(reshape(permute(fluor_cat,[1,2,4,3]),[],n_rois),[],1),[1,3,2,4]));
    fluor_cat_hemidiff_norm = bsxfun(@rdivide,fluor_cat_hemidiff,permute(std(abs(reshape(permute(fluor_cat_hemidiff,[1,2,4,3]),[],n_rois)),[],1),[1,3,2,4]));     
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    smooth_size = 8;
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    
    mua_cat_raw_smoothed = convn(mua_cat_raw,smWin,'same');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,2),'uni',false));
    
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(std(reshape(permute(x,[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,2),'uni',false));
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std);
        
    % Concatenate behavioural data
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
    D.response = cell2mat(cellfun(@(x) x.response,D_all{curr_animal},'uni',false));
    D.repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all{curr_animal},'uni',false));
    
    D.day = trial_day;
    
    % Get trial ID   
    trial_contrast = max(D.stimulus,[],2);
    [~,side_idx] = max(D.stimulus > 0,[],2);
    trial_side = (side_idx-1.5)*2;
    trial_choice = -(D.response-1.5)*2;
    
    % Get time-averaged activity for each trial
    switch plot_act
        case 'fluor'
            use_data = fluor_cat_hemidiff_norm;
            area_labels = {wf_roi.area};
        case 'mua'
            use_data =  mua_cat_norm;
            area_labels = cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false);
    end

    trial_act_timeavg = squeeze(nanmean(use_data(:,plot_t,:,:),2));
    
    % Split activity by contrast*side/choice
    trial_act_timeavg_split = cell(length(sidecontrasts),2);
    trial_act_timeavg_split(:,1) = arrayfun(@(x) ...
        trial_act_timeavg(trial_side.*trial_contrast == sidecontrasts(x) & ...
        trial_choice == -1, ...
        plot_area, plot_align),1:length(sidecontrasts),'uni',false);
    trial_act_timeavg_split(:,2) = arrayfun(@(x) ...
        trial_act_timeavg(trial_side.*trial_contrast == sidecontrasts(x) & ...
        trial_choice == 1, ...
        plot_area, plot_align),1:length(sidecontrasts),'uni',false);
    trial_act_timeavg_split_median = cellfun(@nanmedian,trial_act_timeavg_split);
    
    trial_act_timeavg_split_all = cellfun(@(curr,all) [curr;all], ...
        trial_act_timeavg_split,trial_act_timeavg_split_all,'uni',false);
    
    % Plot average activity across time for reference
    subplot(2,length(D_all),curr_animal); hold on;
    condition_act = grpstats(use_data(:,:,plot_area,plot_align),trial_side.*trial_contrast);
    set(gca,'ColorOrder',colormap_BlueWhiteRed(length(contrasts)-1));
    plot(t,condition_act','linewidth',2);
    rectangle('Position',[t(find(plot_t,1)),min(condition_act(:)), ...
        range(t(plot_t)),range(condition_act(:))],'linewidth',2);
    axis tight;
    ylabel(area_labels{plot_area});
    title('Average activity');
    
    % Plot distribution plot
    subplot(2,length(D_all),length(D_all) + curr_animal); hold on;
    p1 = distributionPlot(trial_act_timeavg_split(:,1),'distWidth',0.5, ...
        'xValues',(1:11)-0.25,'histOri','left','showMM',0,'color',[0.6,0,0.6]);
    p2 = distributionPlot(trial_act_timeavg_split(:,2),'distWidth',0.5, ...
        'xValues',(1:11)+0.25,'histOri','right','showMM',0,'color',[0,0.6,0]);
    plot((1:11),trial_act_timeavg_split_median(:,1),'color',[0.9,0.6,0.9],'linewidth',2)
    plot((1:11),trial_act_timeavg_split_median(:,2),'color',[0.6,0.9,0.6],'linewidth',2)
    
    set(gca,'XTick',1:11,'XTickLabel',cellfun(@num2str,num2cell(sidecontrasts),'uni',false));
    xlabel('Contrast*Side');
    ylabel(area_labels{plot_area});
    legend([p1{1}(6),p2{1}(6)],{'Move left','Move right'});
    title('Trial activity');
    
end

% Plot trials from all animals concatenated
figure; hold on;
p1 = distributionPlot(trial_act_timeavg_split_all(:,1),'distWidth',0.5, ...
    'xValues',(1:11)-0.25,'histOri','left','showMM',0,'color',[0.6,0,0.6]);
p2 = distributionPlot(trial_act_timeavg_split_all(:,2),'distWidth',0.5, ...
    'xValues',(1:11)+0.25,'histOri','right','showMM',0,'color',[0,0.6,0]);

trial_act_timeavg_split_all_median = cellfun(@nanmedian,trial_act_timeavg_split_all);
plot((1:11),trial_act_timeavg_split_all_median(:,1),'color',[0.9,0.6,0.9],'linewidth',2)
plot((1:11),trial_act_timeavg_split_all_median(:,2),'color',[0.6,0.9,0.6],'linewidth',2)

set(gca,'XTick',1:11,'XTickLabel',cellfun(@num2str,num2cell(sidecontrasts),'uni',false));
xlabel('Contrast*Side');
ylabel(area_labels{plot_area});
legend([p1{1}(6),p2{1}(6)],{'Move left','Move right'});
title('Trial activity across animals');


%% Linear regression on day-concatenated activity (weights & expl var)

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = 'all_trial_activity_df_kernel-str_earlymove.mat';
% data_fn = 'all_trial_activity_df_earlymove.mat';
load([data_path filesep data_fn]);

% Get time
framerate = 35.2;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
bhv_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\bhv_processing\use_experiments';
load(bhv_fn);
D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);

% Get widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

n_depths = size(mua_all{1}{1},3);

% Set up variables for model parameters
use_contrasts = [0.06,0.125,0.25,0.5,1];
% use_contrasts = [0.06,0.125];
% use_contrasts = [];

n_regressors = length(use_contrasts)*2 + 1;
fluor_activity_model_params = nan(n_regressors,length(t),n_rois,2,length(D_all));
mua_activity_model_params = nan(n_regressors,length(t),n_depths,2,length(D_all));
mua_predicted_activity_model_params = nan(n_regressors,length(t),n_depths,2,length(D_all));

% Set up variables for explained variance
fluor_model_expl_var = nan(length(t),n_rois,2,2,length(D_all));
mua_model_expl_var = nan(length(t),n_depths,2,2,length(D_all));
mua_predicted_model_expl_var = nan(length(t),n_depths,2,2,length(D_all));

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate behavioural data and select trials to use
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
    D.response = cell2mat(cellfun(@(x) x.response,D_all{curr_animal},'uni',false));
    D.repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all{curr_animal},'uni',false));
    
    D.day = trial_day;
    
    max_contrast = max(D.stimulus,[],2);
    use_trials = ismember(max_contrast,[0,use_contrasts]);
    n_total_trials = length(max_contrast);
    
    % Concatenate fluorescence, get L-R, normalize
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    fluor_cat_hemidiff = cat(3,fluor_cat(:,:,1:n_rois/2,:),fluor_cat(:,:,1:n_rois/2,:) - ...
        fluor_cat(:,:,n_rois/2+1:end,:));
        
    fluor_cat_norm = bsxfun(@rdivide,fluor_cat,permute(std(reshape(permute(fluor_cat,[1,2,4,3]),[],n_rois),[],1),[1,3,2,4]));
    fluor_cat_hemidiff_norm = bsxfun(@rdivide,fluor_cat_hemidiff,permute(std(abs(reshape(permute(fluor_cat_hemidiff,[1,2,4,3]),[],n_rois)),[],1),[1,3,2,4]));     
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    smooth_size = 8;
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    
    mua_cat_raw_smoothed = convn(mua_cat_raw,smWin,'same');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,2),'uni',false));
    
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x,[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,2),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
 
    % Regress from fluor to MUA
    % (by depth: data-less days cause different trial numbers)
    kernel_t = [-0.2,0];
    kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
    zs = [false,false];
    cvfold = 1;
    lambda = 0;
    
    mua_nonan_trials = ~squeeze(any(any(isnan(mua_cat_norm),2),4));
    mua_cat_predicted = nan(size(mua_cat_norm));
    for curr_depth = 1:n_depths
        curr_valid_trials = use_trials & mua_nonan_trials(:,curr_depth);

        [~,predicted_spikes,~] = ...
            AP_regresskernel(reshape(permute( ...
            fluor_cat_norm(curr_valid_trials,:,:,:), ...
            [2,1,4,3]),[],n_rois)', ...
            reshape(permute(mua_cat_norm(curr_valid_trials,:,curr_depth,:), ...
            [2,1,4,3]),[],1)',kernel_frames,lambda,zs,cvfold);
        
        mua_cat_predicted(curr_valid_trials,:,curr_depth,:) = ...
            permute(reshape(predicted_spikes',length(t),sum(curr_valid_trials),2),[2,1,4,3]);
        
    end
            
    %%% Regression (separate regressors to get partial explained variance)
    
    % Set up stim regressors
    if ~isempty(use_contrasts)
        stim_regressors = zeros(n_total_trials,length(use_contrasts),2);
        [~,contrast_l_idx] = ismember(D.stimulus(:,1),use_contrasts','rows');
        [~,contrast_r_idx] = ismember(D.stimulus(:,2),use_contrasts','rows');
        for curr_trial = 1:n_total_trials
            if contrast_l_idx(curr_trial)
                stim_regressors(curr_trial,contrast_l_idx(curr_trial),1) = 1;
            end
            if contrast_r_idx(curr_trial)
                stim_regressors(curr_trial,contrast_r_idx(curr_trial),2) = 1;
            end
        end
        stim_regressors = reshape(stim_regressors,n_total_trials,[]);
    else
        stim_regressors = zeros(n_total_trials,0);
    end
    
    % Set up choice regressors
    choice_regressors = 2*(D.response-1.5);
    
    regressors = cellfun(@transpose, ...
        {stim_regressors, ...
        choice_regressors},'uni',false);
    
    t_shifts = {0,0};
    lambda = 0;
    zs = [false,false];
    cvfold = 10;
    
    % Regress fluorescence all together
    [fluor_params_fit,~,fluor_expl_var] = AP_regresskernel( ...
        cellfun(@(x) x(:,use_trials),regressors,'uni',false)', ...
        reshape(fluor_cat_hemidiff_norm(use_trials,:,:,:),sum(use_trials),[])',t_shifts,lambda,zs,cvfold);
    fluor_params_fit = reshape(fluor_params_fit,n_regressors,length(t),n_rois,2);
    fluor_expl_var = reshape(fluor_expl_var.reduced,length(t),n_rois,2,2);
    
    % Regress MUA by depth (because missing trials = uneven data)
    mua_params_fit = nan(n_regressors,length(t),n_depths,2);
    mua_expl_var = nan(length(t),n_depths,2,2);
    for curr_depth = 1:n_depths
        curr_valid_trials = use_trials & mua_nonan_trials(:,curr_depth);
        curr_regressors = cellfun(@(x) x(:,curr_valid_trials),regressors,'uni',false);
        try
            [curr_mua_params_fit,~,curr_mua_expl_var] = AP_regresskernel(curr_regressors', ...
                reshape(mua_cat_norm(curr_valid_trials, ...
                :,curr_depth,:),sum(curr_valid_trials),[])',t_shifts,lambda,zs,cvfold);
            
            mua_params_fit(:,:,curr_depth,:) = reshape(curr_mua_params_fit,n_regressors,length(t),1,2);
            mua_expl_var(:,curr_depth,:,:) = reshape(curr_mua_expl_var.reduced,[],2,2);            
        catch me
            continue
        end
    end
    
    mua_predicted_params_fit = nan(n_regressors,length(t),n_depths,2);
    mua_predicted_expl_var = nan(length(t),n_depths,2,2);
    for curr_depth = 1:n_depths
        curr_valid_trials =  use_trials & mua_nonan_trials(:,curr_depth);
        curr_regressors = cellfun(@(x) x(:,curr_valid_trials),regressors,'uni',false);
        try
            [curr_mua_predicted_params_fit,~,curr_mua_predicted_expl_var] = AP_regresskernel(curr_regressors', ...
                reshape(mua_cat_predicted(curr_valid_trials, ...
                :,curr_depth,:),sum(curr_valid_trials),[])',t_shifts,lambda,zs,cvfold);
            
            mua_predicted_params_fit(:,:,curr_depth,:) = reshape(curr_mua_predicted_params_fit,n_regressors,length(t),1,2);
            mua_predicted_expl_var(:,curr_depth,:,:) = reshape(curr_mua_predicted_expl_var.reduced,[],2,2);
        catch me
            continue
        end
    end
 
    % Package everything
    fluor_activity_model_params(:,:,:,:,curr_animal)= fluor_params_fit;
    mua_activity_model_params(:,:,:,:,curr_animal)= mua_params_fit;
    mua_predicted_activity_model_params(:,:,:,:,curr_animal)= ...
        reshape(mua_predicted_params_fit,n_regressors,length(t),n_depths,2);
    
    fluor_model_expl_var(:,:,:,:,curr_animal) = fluor_expl_var;
    mua_model_expl_var(:,:,:,:,curr_animal) = mua_expl_var;
    mua_predicted_model_expl_var(:,:,:,:,curr_animal) = mua_predicted_expl_var;
    
end

fluor_activity_model_params_mean = nanmean(fluor_activity_model_params,5);
mua_activity_model_params_mean = nanmean(mua_activity_model_params,5);
mua_predicted_activity_model_params_mean = nanmean(mua_predicted_activity_model_params,5);

fluor_model_expl_var_mean = nanmean(fluor_model_expl_var,5);
mua_model_expl_var_mean = nanmean(mua_model_expl_var,5);
mua_predicted_model_expl_var_mean = nanmean(mua_predicted_model_expl_var,5);

% Plot regressed stim weight
contrasts = [0.06,0.125,0.25,0.5,1];
use_contrasts_idx = ismember(use_contrasts,contrasts);
bwr_colors = colormap_BlueWhiteRed(length(contrasts));
l_colors = bwr_colors(length(contrasts):-1:1,:);
r_colors = bwr_colors(end-length(contrasts)+1:end,:);
plot_color = [l_colors(use_contrasts_idx,:); ...
    r_colors(use_contrasts_idx,:)];

% (Fluor)
figure;
spacing = max(fluor_activity_model_params_mean(:));

subplot(1,2,1); hold on;
for curr_stim = 1:length(use_contrasts)*2;
   AP_stackplot(squeeze(fluor_activity_model_params_mean(curr_stim,:,:,1)),t, ...
       spacing,false,plot_color(curr_stim,:),{wf_roi.area});
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');
ylabel(['\beta stim'])
title('Fluor');

subplot(1,2,2); hold on;
for curr_stim = 1:length(use_contrasts)*2;
   AP_stackplot(squeeze(fluor_activity_model_params_mean(curr_stim,:,:,2)),t, ...
       spacing,false,plot_color(curr_stim,:),{wf_roi.area});
end
line([0,0],ylim,'color','k');
xlabel('Time from move');
ylabel(['\beta stim'])
title('Fluor');

% (MUA)
figure;
spacing = max(mua_activity_model_params_mean(:));

subplot(1,2,1); hold on;
for curr_stim = 1:length(use_contrasts)*2;
   AP_stackplot(squeeze(mua_activity_model_params_mean(curr_stim,:,:,1)),t, ...
       spacing,false,plot_color(curr_stim,:),1:n_depths);
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');
ylabel(['\beta stim'])
title('MUA');

subplot(1,2,2); hold on;
for curr_stim = 1:length(use_contrasts)*2;
   AP_stackplot(squeeze(mua_activity_model_params_mean(curr_stim,:,:,2)),t, ...
       spacing,false,plot_color(curr_stim,:),1:n_depths);
end
line([0,0],ylim,'color','k');
xlabel('Time from move');
ylabel(['\beta stim'])
title('MUA');

% (MUA predicted)
figure;
spacing = max(mua_activity_model_params_mean(:));

subplot(1,2,1); hold on;
for curr_stim = 1:length(use_contrasts)*2;
   AP_stackplot(squeeze(mua_predicted_activity_model_params_mean(curr_stim,:,:,1)),t, ...
       spacing,false,plot_color(curr_stim,:),1:n_depths);
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');
ylabel(['\beta stim'])
title('MUA predicted');

subplot(1,2,2); hold on;
for curr_stim = 1:length(use_contrasts)*2;
   AP_stackplot(squeeze(mua_predicted_activity_model_params_mean(curr_stim,:,:,2)),t, ...
       spacing,false,plot_color(curr_stim,:),1:n_depths);
end
line([0,0],ylim,'color','k');
xlabel('Time from move');
ylabel(['\beta stim'])
title('MUA predicted');

% Plot regressed choice weight
figure;
curr_subplot = 1;
for curr_modality = 1:3
    for curr_alignment = 1:2
        switch curr_alignment
            case 1
                align_text = 'stim';
            case 2
                align_text = 'move';
        end
        subplot(3,2,curr_subplot); hold on;
        switch curr_modality
            case 1
                set(gca,'ColorOrder',[autumn(n_rois/2);winter(n_rois/2)]);
                plot(t,squeeze(fluor_activity_model_params_mean(end,:,:,curr_alignment)));
                xlabel(['Time from ' align_text]);
                ylabel(['\beta choice']);
                axis tight
                line([0,0],ylim,'color','k');
                line(xlim,[0,0],'color','k');
            case 2
                set(gca,'ColorOrder',copper(n_depths));
                plot(t,squeeze(mua_activity_model_params_mean(end,:,:,curr_alignment)));
                xlabel(['Time from ' align_text]);
                ylabel(['\beta choice']);
                axis tight
                line([0,0],ylim,'color','k');
                line(xlim,[0,0],'color','k');
            case 3
                set(gca,'ColorOrder',copper(n_depths));
                plot(t,squeeze(mua_predicted_activity_model_params_mean(end,:,:,curr_alignment)));
                xlabel(['Time from ' align_text]);
                ylabel(['\beta choice']);
                axis tight
                line([0,0],ylim,'color','k');
                line(xlim,[0,0],'color','k');
        end
        curr_subplot = curr_subplot + 1;
    end
end

% Plot explained variance
figure;
curr_subplot = 1;
for curr_alignment = 1:2
    switch curr_alignment
        case 1
            align_text = 'stim';
        case 2
            align_text = 'move';
    end
    for curr_modality = 1:3
        for curr_param = 1:2            
            switch curr_param
                case 1
                    param_text = 'stim';
                case 2
                    param_text = 'choice';
            end
            
            subplot(6,2,curr_subplot); hold on;
            switch curr_modality
                case 1
                    set(gca,'ColorOrder',[autumn(n_rois/2);winter(n_rois/2)]);
                    plot(t,squeeze(fluor_model_expl_var_mean(:,:,curr_alignment,curr_param)));
                    xlabel(['Time from ' align_text]);
                    ylabel(['Expl var: ' param_text]);
                    axis tight
                    line([0,0],ylim,'color','k');
                    line(xlim,[0,0],'color','k');
                case 2
                    set(gca,'ColorOrder',copper(n_depths));
                    plot(t,squeeze(mua_model_expl_var_mean(:,:,curr_alignment,curr_param)));
                    xlabel(['Time from ' align_text]);
                    ylabel(['Expl var: ' param_text]);
                    axis tight
                    line([0,0],ylim,'color','k');
                    line(xlim,[0,0],'color','k');
                case 3
                    set(gca,'ColorOrder',copper(n_depths));
                    plot(t,squeeze(mua_predicted_model_expl_var_mean(:,:,curr_alignment,curr_param)));                    
                    xlabel(['Time from ' align_text]);
                    ylabel(['Expl var: ' param_text]);
                    axis tight
                    line([0,0],ylim,'color','k');
                    line(xlim,[0,0],'color','k');
            end
            curr_subplot = curr_subplot + 1;
        end
    end
end















