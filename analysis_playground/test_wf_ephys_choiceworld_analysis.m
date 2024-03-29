%% PSTH to choiceworld conditions

stimIDs = signals_events.trialSideValues.*signals_events.trialContrastValues;
% stimIDs = discretize(stimIDs,[-Inf,-0.125,-0.01,0.01,0.25,Inf],[-2,-1,0,1,2]);

use_spikes_idx = ismember(spike_templates,find(template_depths >= str_depth(1) & template_depths <= str_depth(1)+100));
% use_spikes_idx = ismember(spike_templates,find(template_depths >= 500 & template_depths <= 800));
% use_spikes_idx = ismember(spike_templates,intersect(find(template_depths >= 500 & template_depths <= 1500),find(msn)));

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
for curr_stim_idx = 1:length(unique_stims)
    
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

depth_group = discretize(spike_depths,depth_group_edges);

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

depth_group = discretize(spike_depths,depth_group_edges);

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

% use_spikes_idx = ismember(spike_templates,find(template_depths >= 1000 & template_depths <= 2000));
use_spikes_idx =ismember(spike_templates,find(template_depths > 1000 & template_depths < 2000)) &...
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

use_spikes_idx = ismember(spike_templates,find(template_depths >= depth_group_edges(3) & template_depths <= depth_group_edges(4)));
%use_spikes_idx = ismember(spike_templates,find(template_depths >= 3000 & template_depths <= 4000));
% use_spikes_idx = ismember(spike_templates,intersect(find(template_depths >= 1000 & template_depths <= 2000),find(msn)));

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
depth_group = discretize(spike_depths,depth_group_edges);

use_spikes_idx = depth_group == 2;

% use_spikes_idx = ismember(spike_templates,find(template_depths >= 500 & template_depths <= 3500));
% use_spikes_idx = ismember(spike_templates,intersect(find(template_depths >= 2000 & template_depths <= 3000),find(msn)));

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
AP_imscroll(a,t_surround);
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
use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > depth_edges(1) & template_depths < depth_edges(2))));
% use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 500 & template_depths < 1500)) &...
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

use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > depth_edges(1) & template_depths < depth_edges(2))));
% use_spikes = spike_times_timeline(ismember(spike_templates, ...
%     find(template_depths > depth_edges(1) & template_depths < depth_edges(2))) & ...
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

use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > depth_edges(1) & template_depths < depth_edges(2))));
% use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > depth_edges(1) & template_depths < depth_edges(2))) &...
%     ismember(spike_templates,find(msn)));

skip_seconds = 60*1;

time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > depth_edges(1) & template_depths < depth_edges(2))));
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
use_spikes_idx = true(size(spike_depths));
% use_spikes_idx = spike_depths > 2000 & spike_depths < 3000;
% use_spikes_idx = spike_depths >= str_depth(1) & spike_depths <= str_depth(2);
% use_spikes_idx = aligned_str_depth_group == 1;
% use_spikes_idx = spike_templates == 30;

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
    stim_to_feedback < 1 & ...
    ~isnan(wheel_move_time);

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
    wheel_move_time(use_trials)',raster_window,trial_choice(use_trials));
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
set(gcf,'Name','Templates sorted by choice difference');


%% Ephys: rasters by template

% Pull out spikes within striatum
% use_spikes_idx = spike_depths >= str_depth(1) & spike_depths <= str_depth(2);
use_spikes_idx = aligned_str_depth_group == 1;

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
AP_imscroll(template_psth_smooth)
line(repmat(find(t_bins > 0,1),2,1),ylim,'color','r');

% Choice difference psth
% [contrast,side,choice,timing]
left_early = ismember(conditions(:,3:4),[-1,1],'rows');
right_early = ismember(conditions(:,3:4),[1,1],'rows');

move_diff = permute(nanmean(template_psth_smooth(right_early,:,:,:),1) - ...
    nanmean(template_psth_smooth(left_early,:,:,:),1),[3,2,4,1]);

[~,sort_idx] = sort(template_depths(use_templates_unique));

AP_imscroll(move_diff(sort_idx,:,:))
caxis([-abs(max(move_diff(:))),abs(max(move_diff(:)))]);
colormap(colormap_BlueWhiteRed);
line(repmat(find(t_bins > 0,1),2,1),ylim,'color','k');


%% Each depth: get left/right PCA groups (after above)
% the point of this was to get 2 basic groups of cells at each depth to
% get at more population-dynamics stuff

% Group striatum depths
n_depths = 6;
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
depth_group = discretize(template_depths(use_templates_unique),depth_group_edges);

go_left = ismember(conditions(:,2:3),[1,-1],'rows');
go_right = ismember(conditions(:,2:3),[1,1],'rows');

curr_depth = 4;
curr_align = 2;
pca_data = zscore([squeeze(nanmean(template_psth_smooth(go_left,:,depth_group == curr_depth,curr_align),1)); ...
    squeeze(nanmean(template_psth_smooth(go_right,:,depth_group == curr_depth,curr_align),1))],[],1);

[coeff,score,latent] = pca(pca_data);



%% Ephys: rasters by depth

% % (for equally separated depths)
n_depths = 4;
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
depth_group = discretize(spike_depths,depth_group_edges);

% (to use aligned striatum depths)
% n_depths = n_aligned_depths;
% depth_group = aligned_str_depth_group;

% % (for manual depth)
% depth_group_edges = [1500,2200];
% n_depths = length(depth_group_edges) - 1;
% [depth_group_n,depth_group] = histc(spike_depths,depth_group_edges);

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
AP_imscroll(depth_psth_smooth)
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

% (to no use ddf)
ddf = roi_trace(:,1:end-1);

% % (get ddf)
% ddf = diff(roi_trace,[],2);
% ddf(ddf < 0) = 0;

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

AP_imscroll(roi_psth,{wf_roi.area})
line(repmat(find(t > 0,1),2,1),ylim,'color','r');
colormap(colormap_BlueWhiteRed);
caxis([-max(abs(caxis)),max(abs(caxis))]);

% Choice difference psth
% [contrast,side,choice,timing]
left_early = ismember(conditions(:,3:4),[-1,1],'rows');
right_early = ismember(conditions(:,3:4),[1,1],'rows');

move_diff = permute(squeeze(nanmean(roi_psth(right_early,:,:,:),1) - ...
    nanmean(roi_psth(left_early,:,:,:),1)),[2,1,3]);

AP_imscroll(move_diff,{'Stim-aligned','Move-aligned'})
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
depth_group = discretize(spike_depths,depth_group_edges);

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

surround_window = [-0.5,3];
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

AP_imscroll(im);
axis image off;

% (used for df/ddf comparison for data club)
% plot_im = [mat2gray(im(:,:,1:end-1,3),prctile(reshape(im(:,:,:,3),[],1),[1,99])), ...
%     mat2gray(ddf_im(:,:,:,3),prctile(reshape(ddf_im(:,:,:,3),[],1),[1,99]))];
plot_im = [mat2gray(im(:,:,:,1),prctile(reshape(im,[],1),[1,99]))];
fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\passive\kalatsky_left';

plot_t = t_surround > -0.2 & t_surround < 3;
% AP_movie2avi(plot_im(:,:,plot_t),35*2,[0,1],p,fn,cellfun(@num2str,num2cell(t_surround(plot_t)),'uni',false));
AP_movie2avi(plot_im(:,:,plot_t),35*2,[0,1],p,fn);



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
AP_reference_outline('ccf_aligned','k');%AP_reference_outline('retinotopy','m');
for curr_roi = 1:n_rois
    curr_roi_boundary = cell2mat(bwboundaries(roi_cat(:,:,curr_roi)));
    patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),roi_col(curr_roi,:));
    
    text(nanmean(curr_roi_boundary(:,2)),nanmean(curr_roi_boundary(:,1)), ...
        wf_roi(curr_roi).area,'FontSize',12,'HorizontalAlignment','center')
end
axis image off;


%% Load and process widefield choiceworld mean (GENERAL - FULL IMAGE)

surround_window = [-0.5,2];
upsample_rate = 3;

framerate = 35;
surround_samplerate = 1/(framerate*upsample_rate);
t_surround = surround_window(1):surround_samplerate:surround_window(2);
t_df = conv2(t_surround,[1,1]/2,'valid');
conditions = [-1,-0.5,-0.25,-0.125,-0.06,0,0.06,0.125,0.25,0.5,1];

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';

trialtype_align = 'move';
trialtype_timing = 'earlymove';
trialtype_success = 'hit';

trialtype = [trialtype_align ' ' trialtype_timing ' ' trialtype_success];

wf_fn = [data_path filesep 'im_' trialtype_align '_' trialtype_timing '_' trialtype_success '_combined.mat'];
load(wf_fn);

% Get ddf
ddf = diff(im_aligned_avg_combined,[],3);
ddf(ddf < 0) = 0;

AP_imscroll(ddf,t_df);
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
% AP_imscroll(ddf_diff,t_df)
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

framerate = 35;
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

AP_imscroll(roi_psth_mean,{wf_roi.area})
line(repmat(find(t > 0,1),2,1),ylim,'color','r');

% (make L-R from here on)
roi_psth_mean(:,:,size(wf_roi,1)+1:end,:) = roi_psth_mean(:,:,1:size(wf_roi,1),:) - ...
    roi_psth_mean(:,:,size(wf_roi,1)+1:end,:);
AP_imscroll(roi_psth_mean(:,:,size(wf_roi,1)+1:end,:),cellfun(@(x) ['\Delta' x],{wf_roi.area},'uni',false));
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

AP_imscroll(depth_psth_mean, ...
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

n_depths = 4;
protocol = 'stimKalatsky';
% protocol = 'AP_choiceWorldStimPassive';

data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\passive'];

mua_fn = [data_path filesep 'mua_' protocol '_' num2str(n_depths) '_depths'];
load(mua_fn);

% Parameters defined in batch
raster_window = [-0.5,3];
psth_bin_size = 0.001;
t_edges = raster_window(1):psth_bin_size:raster_window(2);
t = t_edges(1:end-1) + diff(t_edges);

mua_all = cellfun(@transpose,{batch_vars.stim_aligned_mua}','uni',false);

mua_psth = nan(3,length(t),n_depths,length(batch_vars));
for curr_animal = 1:length(batch_vars)
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths),'uni',false));
    
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x,[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
          
    stimID_cat = vertcat(batch_vars(curr_animal).stimIDs{:});
    for curr_depth = 1:n_depths
        mua_psth(:,:,curr_depth,curr_animal) = grpstats(mua_cat_norm(:,:,curr_depth),stimID_cat);
    end
    
end

mua_psth_mean = nanmean(mua_psth,4);

trace_spacing = 0.8;
n_conditions = size(mua_psth_mean,1);
n_depths = size(mua_psth_mean,3);
col = copper(n_conditions);
figure; hold on;
p = gobjects(n_depths,n_conditions);
for curr_condition = 1:n_conditions
    p(:,curr_condition) = ...
        AP_stackplot(permute(mua_psth_mean(curr_condition,:,:),[2,3,1]), ...
        t,trace_spacing,false,col(curr_condition,:));
end
axis tight;
ylabel('Spikes')
xlabel('Time from stim onset (s)');
line([0,0],ylim,'linestyle','--','color','k');
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

n_aligned_depths = 4;

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_ephys';

protocol = 'vanillaChoiceworld';
% protocol = 'stimSparseNoiseUncorrAsync';
% protocol = 'stimKalatsky';
% protocol = 'AP_choiceWorldStimPassive';

map_fn = [data_path filesep 'wf_ephys_maps_' protocol '_' num2str(n_aligned_depths) '_depths_kernel'];
load(map_fn);

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
% animals = {'AP032','AP033','AP034','AP035','AP036'};

% (scale r_px's because different lambdas give different weights)
% (do in a loop because memory can't handle a cellfun??)
n_depths = n_aligned_depths;

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

% % Save kernels for extracting from fluorescence later (in pixel and V)
% % (turn kernel pixels into V from master U)
% load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');
% n_px = size(U_master,1)*size(U_master,2);
% kernel_V = reshape(reshape(U_master,n_px,[])\reshape(r_px_mean,n_px,[]), ...
%     size(U_master,3),size(r_px_mean,3),n_depths);
% 
% ctx_str_kernel.t = t;
% ctx_str_kernel.px = r_px_mean;
% ctx_str_kernel.V = kernel_V;
% 
% ctx_str_kernel_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\ctx_str_kernel';
% save(ctx_str_kernel_fn,'ctx_str_kernel');
% disp('Saved ctx-str kernels');


% Plot weights over time (only use ipsi ROIs)
t = linspace(-0.3,0.3,size(r_px_mean,3));
AP_imscroll(r_px_mean,t)  
axis image;
caxis([-1,1]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k');AP_reference_outline('retinotopy','m');

% Plot weights at t = 0
figure('Name','Weights at t = 0'); colormap(brewermap([],'*RdBu'));
for curr_depth = 1:n_depths
    subplot(1,n_depths,curr_depth);
    imagesc(r_px_mean(:,:,t == 0,curr_depth));
    AP_reference_outline('ccf_aligned','k');
    caxis([-1,1]);
    axis image off
end

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
use_t = t >= 0 & t <= 0;
r_px_max = squeeze(nanmean(r_px_mean(:,:,use_t,:),3));
r_px_max_norm = bsxfun(@rdivide,r_px_max, ...
    permute(max(reshape(r_px_max,[],n_depths),[],1),[1,3,2]));
r_px_max_norm(isnan(r_px_max_norm)) = 0;
r_px_max_norm(r_px_max_norm < 0) = 0;
r_px_com = sum(bsxfun(@times,r_px_max_norm,permute(1:n_depths,[1,3,2])),3)./sum(r_px_max_norm,3);
com_colored = ind2rgb(round(mat2gray(r_px_com,[1,n_depths])*255),jet(255));

figure;
p = imagesc(com_colored);
axis image off
weight_norm = mat2gray(max(r_px_max,[],3),[0,double(prctile(reshape(max(r_px_max,[],3),[],1),100))]);
set(p,'AlphaData',weight_norm);

c = colorbar;
ylabel(c,'Depth (fraction)');
colormap(c,jet);
set(c,'YDir','reverse');
set(c,'YTick',linspace(0,1,n_depths));
set(c,'YTickLabel',1:n_depths);
 
AP_reference_outline('retinotopy','m');
AP_reference_outline('ccf_aligned','k');

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

%% ~~~~~ PLOT ALIGNMENT PARAMETERS ~~~~~~

%% Plot all striatum boundaries

ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
ephys_align_fn = ['ephys_depth_align'];
load([ephys_align_path filesep ephys_align_fn])

% Parameters from batch
n_corr_groups = 40;
depth_group_edges = linspace(0,3820,n_corr_groups+1);
depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

for curr_animal = 1:length(ephys_depth_align)
    n_days = length(ephys_depth_align(curr_animal).mua_corr);
    figure; 
    for curr_day = 1:n_days
        mua_corr = ephys_depth_align(curr_animal).mua_corr{curr_day};
        template_depths = ephys_depth_align(curr_animal).template_depths{curr_day};
        str_depth = ephys_depth_align(curr_animal).str_depth(curr_day,:);
        
        subplot(2,n_days,curr_day);
        imagesc(depth_group_centers,depth_group_centers,mua_corr);
        caxis([0,1]); colormap(hot); axis square;
        line(xlim,[str_depth(1),str_depth(1)],'color','r')
        line([str_depth(1),str_depth(1)],ylim,'color','r')
        line(xlim,[str_depth(2),str_depth(2)],'color','r')
        line([str_depth(2),str_depth(2)],ylim,'color','r')
        
        subplot(2,n_days,n_days+curr_day);
        plotSpread(template_depths,'distributionColor','k');
        set(gca,'YDir','reverse');
        line(xlim,[str_depth(1),str_depth(1)]);
        line(xlim,[str_depth(2),str_depth(2)]);        
    end   
    set(gcf,'Name',ephys_depth_align(curr_animal).animal);
end

%% Plot lambda fits, save explained var use experiments

% Load lambda from previously estimated and saved
lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
load(lambda_fn);

max_days = max(cellfun(@length,{ctx_str_lambda.day}));

expl_vars = cellfun(@(x) cellfun(@max,x),{ctx_str_lambda.explained_var_lambdas},'uni',false);
max_expl_var = max(horzcat(expl_vars{:}));

figure;
for curr_animal = 1:length(ctx_str_lambda)
    n_days = length(ctx_str_lambda(curr_animal).day);
    for curr_day = 1:n_days       
        
        subplot(length(ctx_str_lambda),max_days,(curr_animal-1)*max_days+curr_day); hold on;
        plot(ctx_str_lambda(curr_animal).lambdas{curr_day}, ...
            ctx_str_lambda(curr_animal).explained_var_lambdas{curr_day},'.k','linewidth',2);
        plot(ctx_str_lambda(curr_animal).lambdas{curr_day}, ...
            smooth(ctx_str_lambda(curr_animal).explained_var_lambdas{curr_day},10),'linewidth',1);
        
%         ylim([0,max_expl_var]);

        line(repmat(ctx_str_lambda(curr_animal).best_lambda(curr_day),1,2),ylim,'color','r');
        
    end   
end

figure;

subplot(1,2,1); hold on; set(gca,'ColorOrder',jet(length({ctx_str_lambda.animal})));
for i = 1:length(expl_vars)
    plot(ctx_str_lambda(i).best_lambda, expl_vars{i},'.','MarkerSize',15);
end
xlabel('\lambda')
ylabel('Explained variance')
legend({ctx_str_lambda.animal})

subplot(1,2,2); hold on; set(gca,'ColorOrder',jet(length({ctx_str_lambda.animal})));
for i = 1:length(expl_vars)
    plot(ctx_str_lambda(i).best_lambda, sqrt(expl_vars{i}),'.','MarkerSize',15);
end
xlabel('\lambda')
ylabel('\surdExplained variance')
legend({ctx_str_lambda.animal})

lambda_cutoff = 18;
sqrt_expl_var_cutoff = 0.1;
axis tight;
line(repmat(lambda_cutoff,2,1),ylim,'linewidth',2);
line(xlim,repmat(sqrt_expl_var_cutoff,2,1),'linewidth',2);

% % Set cutoff for explained variance, eliminate other experiments
% use_experiments = cellfun(@(x) sqrt(x) > sqrt_expl_var_cutoff,expl_vars,'uni',false);
% 
% save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
% save_fn = 'expl_var_use_experiments';
% save([save_path filesep save_fn],'use_experiments');
% disp('Saved explained variance experiment cutoff');

%% Plot kernel matches

% Load and plot the kernel templates
n_aligned_depths = 3;
kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template_' num2str(n_aligned_depths) '_depths.mat'];
load(kernel_template_fn);
n_kernels = n_aligned_depths;
figure; 
for i = 1:n_kernels
    subplot(1,n_kernels,i);
    imagesc(kernel_template(:,:,i));
    axis image off;
    caxis([-prctile(abs(kernel_template(:)),99),prctile(abs(kernel_template(:)),99)]);
    colormap(brewermap([],'*RdBu'));
    AP_reference_outline('ccf_aligned','k');
end

% Load the kernel template matches
kernel_match_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
kernel_match_fn = ['ephys_kernel_align_' num2str(n_aligned_depths) '_depths.mat'];
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
    xlabel('Striatum depth')
    ylabel('Kernel match')
    title(ephys_kernel_align(curr_animal).animal);
    
    subplot(2,length(ephys_kernel_align),curr_animal+length(ephys_kernel_align)); hold on;
    set(gca,'ColorOrder',copper(size(kernel_match_all{curr_animal},2)));
    plot(kernel_match_all{curr_animal},'linewidth',2);
    xlabel('Striatum depth')
    ylabel('Kernel match (cleaned)')
    title(ephys_kernel_align(curr_animal).animal);
end

% Plot explained variance against kernel correlation

% Load lambda from previously estimated and saved
lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
load(lambda_fn);

expl_vars = cellfun(@(x) cellfun(@max,x),{ctx_str_lambda.explained_var_lambdas},'uni',false);

% Get maximum kernel correlation 
session_kernel_corr = cellfun(@(x) cellfun(@(x) mean(max(x,[],2)),x), ...
        {ephys_kernel_align.kernel_corr},'uni',false);

figure; hold on;
for i = 1:length(expl_vars)
    plot(session_kernel_corr{i}, expl_vars{i},'.','MarkerSize',15);
end
xlabel('Mean kernel correlation')
ylabel('Explained variance')
legend({ctx_str_lambda.animal})


% Plot correlations across kernel grouping transitions 

% Load kernels by depths
kernel_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
kernel_fn = ['ephys_kernel_depth'];
load([kernel_path filesep kernel_fn])

% Get within-recording k_px correlation
k_px_cat = [ephys_kernel_depth(:).k_px];
k_px_cat_corr = cellfun(@(x) corrcoef(reshape(x,[],size(x,3))),k_px_cat,'uni',false);

% Get all template matches
k_match_cat = [ephys_kernel_align(:).kernel_match];

% Pull out n to non-n transitions, align and average
max_depths = max(cellfun(@(x) size(x,3),[ephys_kernel_depth.k_px]));
k_px_cat_corr_pad = cellfun(@(x) padarray(x,2*max_depths-size(x),NaN,'post'), ...
    k_px_cat_corr,'uni',false);

k_corr_transition_mean = nan(2*max_depths,2*max_depths,n_aligned_depths,2);
for curr_template = 1:n_aligned_depths
    
    curr_leave_border = cellfun(@(x) find(x(1:end-1) == curr_template & ...
        x(2:end) ~= curr_template & ~isnan(x(2:end)),1),k_match_cat,'uni',false);
    curr_use_recordings = cellfun(@any,curr_leave_border);
    if any(curr_use_recordings)
        curr_corr_pad_aligned = cell2mat(permute(cellfun(@(corr,border) ...
            circshift(corr,[max_depths-border,max_depths-border]), ...
            k_px_cat_corr_pad(curr_use_recordings), ...
            curr_leave_border(curr_use_recordings),'uni',false),[1,3,2]));
        k_corr_transition_mean(:,:,curr_template,1) = nanmean(curr_corr_pad_aligned,3);
    end
    
    curr_enter_border = cellfun(@(x) find(x(1:end-1) ~= curr_template & ...
        x(2:end) == curr_template & ~isnan(x(1:end-1)),1),k_match_cat,'uni',false);
    curr_use_recordings = cellfun(@any,curr_enter_border);
    if any(curr_use_recordings)
        curr_corr_pad_aligned = cell2mat(permute(cellfun(@(corr,border) ...
            circshift(corr,[max_depths-border,max_depths-border]), ...
            k_px_cat_corr_pad(curr_use_recordings), ...
            curr_enter_border(curr_use_recordings),'uni',false),[1,3,2]));
        k_corr_transition_mean(:,:,curr_template,2) = nanmean(curr_corr_pad_aligned,3);
    end
    
end
figure;colormap(gray);
for curr_depth = 1:n_aligned_depths
   subplot(2,n_aligned_depths,curr_depth);
   imagesc(k_corr_transition_mean(:,:,curr_depth,1));
   line([max_depths+0.5,max_depths+0.5],ylim,'color','r');
   line(xlim,[max_depths+0.5,max_depths+0.5],'color','r');
   caxis([0,1]);
   axis square;
   title('Exiting');
   
   subplot(2,n_aligned_depths,n_aligned_depths+curr_depth);
   imagesc(k_corr_transition_mean(:,:,curr_depth,2));
   line([max_depths+0.5,max_depths+0.5],ylim,'color','r');
   line(xlim,[max_depths+0.5,max_depths+0.5],'color','r');
   caxis([0,1]);
   axis square;
   title('Entering');
end


%%%%%%%%% TO DO: plot the correlation across kernels?
%%%% also something to justify 4 groups
%%% also transitions between borders

% Get all depth kernels
k_px_cat = [ephys_kernel_depth(:).k_px];
k_px_cat = cat(3,k_px_cat{:});
k_px_cat_reshape = reshape(k_px_cat,[],size(k_px_cat,3));

% Concatenate kernels and matches
k_match_cat = cell2mat([ephys_kernel_align(:).kernel_match]');
k_px_cat = cell2mat(permute([ephys_kernel_depth(:).k_px],[1,3,2]));

% sort by depth instead?
total_depths = 1:max(cellfun(@(x) size(x,3),[ephys_kernel_depth.k_px]));
k_px_depths = cellfun(@(x) total_depths(end-size(x,3)+1:end),[ephys_kernel_depth.k_px],'uni',false);
k_px_depth_cat = horzcat(k_px_depths{:});



% % to look at inter/intra-kernel correlation?
% k_px_corr = corrcoef(k_px_cat_reshape(:,use_k_px(k_depth_sort_idx)));
% k_px_corr(logical(eye(size(k_px_corr)))) = NaN;
% n_kidx_depths = histcounts(kidx_depth);
% k_px_corr_split = mat2cell(k_px_corr,n_kidx_depths,n_kidx_depths);


%% Create ROIs from kernel templates

% Load kernel templates
n_aligned_depths = 4;
kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template_' num2str(n_aligned_depths) '_depths.mat'];
load(kernel_template_fn);

kernel_bw = false(size(kernel_template));

% % Set cutoff by percentile, get rid of small islands
% kernel_cutoff_prctile = 95;
% smallest_area = 2000;
% 
% for curr_kernel = 1:size(kernel_bw,3)
%     curr_kernel_cutoff = prctile(reshape(kernel_template(:,:,curr_kernel),[],1),kernel_cutoff_prctile);
%     curr_kernel_bw = kernel_template(:,:,curr_kernel) > curr_kernel_cutoff;
%     kernel_bw(:,:,curr_kernel) = bwareaopen(curr_kernel_bw,smallest_area);
% end

% Set cutoff by std, get rid of small islands
kernel_std_prctile = 2;
smallest_area = 2000;

for curr_kernel = 1:size(kernel_bw,3)
    curr_kernel_cutoff = std(abs(reshape(kernel_template(:,:,curr_kernel),[],1)));
    curr_kernel_bw = kernel_template(:,:,curr_kernel) > kernel_std_prctile*curr_kernel_cutoff;
    kernel_bw(:,:,curr_kernel) = bwareaopen(curr_kernel_bw,smallest_area);
end

% Only use the ipsilateral side
bregma = allenCCFbregma;
ccf_tform_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\ccf_tform'];
load(ccf_tform_fn);

um2pixel = 20.6;
bregma_resize = bregma*(10/um2pixel);
bregma_align = [bregma_resize([3,1]),1]*ccf_tform.T;
kernel_bw(:,round(bregma_align(1)):end,:) = false;

% Plot the kernels and ROIs
figure; colormap(gray);
for i = 1:n_aligned_depths
    subplot(2,n_aligned_depths,i);
    imagesc(kernel_template(:,:,i));
    AP_reference_outline('ccf_aligned','r');
    axis image off;
    
    subplot(2,n_aligned_depths,i+n_aligned_depths);
    imagesc(kernel_bw(:,:,i));
    AP_reference_outline('ccf_aligned','r');
    axis image off;
end

% Save kernel ROIs
kernel_roi = kernel_bw;
kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
save(kernel_roi_fn,'kernel_roi');
disp('Saved kernel ROIs');


%% ~~~~~ ~~~~~~~~~~~~~~~~ ~~~~~~


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
n_depths = 4;

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

plot_area = 3;

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
n_depths = 4;

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

framerate = 35;
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
framerate = 35;
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
framerate = 35;
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
framerate = 35;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):sample_rate:raster_window(2);

% Take the mean progressively and overwrite - file's too big
for curr_animal = 1:length(batch_vars)
    batch_vars(curr_animal).pV = nanmean(batch_vars(curr_animal).pV,5);
end

pV = nanmean(cat(5,batch_vars.pV),5);

AP_imscroll(pV,t);
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
framerate = 35;
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

n_aligned_depths = 4;

% Load data
data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
data_fn = ['activity_sessioncat_logistic_regression_earlymove_kernel-str_' num2str(n_aligned_depths) '_depths.mat'];
% data_fn = ['activity_sessioncat_logistic_regression_earlymove'];
load([data_path filesep data_fn])

% Get time
framerate = 35;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):sample_rate:raster_window(2);

% Get widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

n_depths = size(loglik_increase_mua,2);

% Set colors
fluor_colors = [autumn(n_rois/2);winter(n_rois/2)];
mua_colors = copper(n_depths);

% The offsets are messed up? fix
loglik_increase_fluor = bsxfun(@minus,loglik_increase_fluor,loglik_increase_fluor(10,:,:,:));
loglik_increase_mua = bsxfun(@minus,loglik_increase_mua,loglik_increase_mua(10,:,:,:));

% Plot mean
yrange = [min(reshape(nanmean(loglik_increase_fluor,4),[],1)), ...
    max(reshape(nanmean(loglik_increase_fluor,4),[],1))];

figure; 
p1 = subplot(2,2,1); hold on;
set(gca,'ColorOrder',fluor_colors);
plot(t,nanmean(loglik_increase_fluor(:,:,1,:),4))
line([0,0],yrange,'color','k');
xlabel('Time from stim')
ylabel('Relative loglikelihood (bpt)');

p2 = subplot(2,2,2); hold on;
set(gca,'ColorOrder',fluor_colors);
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
set(gca,'ColorOrder',[fluor_colors;mua_colors]);
plot(t,nanmean(loglik_increase_fluor(:,:,2,:),4))
plot(t,nanmean(loglik_increase_mua(:,:,2,:),4))

% Max predicatbility in cortex and striatum for each animal
% (max difference for all regions and mice)
n_animals = size(loglik_increase_fluor,4);

loglik_diff_fluor = squeeze(max(loglik_increase_fluor,[],1) - min(loglik_increase_fluor,[],1));
loglik_diff_mua = squeeze(max(loglik_increase_mua,[],1) - min(loglik_increase_mua,[],1));

figure; hold on;
scatter(squeeze(max(loglik_diff_fluor(10,2,:),[],1)), ...
    squeeze(max(loglik_diff_mua(3,2,:),[],1)), ...
    80,parula(n_animals),'filled','MarkerEdgeColor','k')

axis([0,max([ylim,xlim]),0,max([ylim,xlim])]);
axis image

line([0,max([ylim,xlim])],[0,max([ylim,xlim])],'color','k')

ylabel('Max fluor pred (\DeltaAM)');
xlabel('Max str pred (Depth 3)');

%% Load logistic regression (all within-modal data simultaneously)

% Load data
data_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
data_fn = ['activity_sessioncat_logistic_regression_neucombined_earlymove_kernel-str'];
load([data_path filesep data_fn])

% Get time
framerate = 35;
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

%% Plot mean day-concatenated activity (passive)

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_df_kernel-str_passive_' num2str(n_aligned_depths) '_depths.mat'];

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

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

% MUA depths
n_depths = size(mua_all{1}{1},3);

n_conditions = max(D_all{1}{1}.stimulus);

fluor_psth = nan(n_conditions,length(t),n_rois,2,length(D_all));
mua_psth = nan(n_conditions,length(t),n_depths,2,length(D_all));
mua_predicted_psth = nan(n_conditions,length(t),n_depths,2,length(D_all));

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
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths)), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths),'uni',false));
    
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x,[1,2,4,3]),[],n_depths),[]), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
          
    % Regress from fluor to MUA
    % (by depth: data-less days cause different trial numbers)
    kernel_t = [-0.1,0.1];
    kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
    zs = [false,false];
    cvfold = 1;
    lambda = 0;
    
    mua_nonan_trials = ~squeeze(any(any(isnan(mua_cat_norm),2),4));
    mua_cat_predicted = nan(size(mua_cat_norm));
    for curr_depth = 1:n_depths
        curr_valid_trials = mua_nonan_trials(:,curr_depth);

        [~,predicted_spikes,~] = ...
            AP_regresskernel(reshape(permute( ...
            fluor_cat_norm(curr_valid_trials,:,:,:), ...
            [2,1,4,3]),[],n_rois)', ...
            reshape(permute(mua_cat_norm(curr_valid_trials,:,curr_depth,:), ...
            [2,1,4,3]),[],1)',kernel_frames,lambda,zs,cvfold);
        
        mua_cat_predicted(curr_valid_trials,:,curr_depth,:) = ...
            permute(reshape(predicted_spikes',length(t),sum(curr_valid_trials)),[2,1,4,3]);
        
    end
    
    % Concatenate behavioural data
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
    
    D.day = trial_day;
    
    % Get trial ID     
    trial_conditions = D.stimulus;
    [~,trial_id] = ismember(trial_conditions,[1:max(D.stimulus)]','rows');
    
    % Save aligned average by trial ID
    for curr_roi = 1:n_rois
        fluor_psth(unique(trial_id),:,curr_roi,1,curr_animal) = grpstats(fluor_cat_hemidiff_norm(:,:,curr_roi),trial_id);
    end    
    
    for curr_depth = 1:n_depths
        mua_psth(unique(trial_id),:,curr_depth,1,curr_animal) = grpstats(mua_cat_norm(:,:,curr_depth),trial_id);
    end    
    
    for curr_depth = 1:n_depths
        mua_predicted_psth(unique(trial_id),:,curr_depth,1,curr_animal) = grpstats(mua_cat_predicted(:,:,curr_depth),trial_id);
    end    
 
end

fluor_psth_mean = nanmean(fluor_psth,5);
mua_psth_mean = nanmean(mua_psth,5);
mua_predicted_psth_mean = nanmean(mua_predicted_psth,5);

% Set up plot conditions
plot_conditions = 1:size(fluor_psth_mean,1);
plot_color = copper(length(plot_conditions));

% Stack MUA plots 
figure; hold on;
for curr_condition_idx = 1:length(plot_conditions)
    curr_condition = plot_conditions(curr_condition_idx);
    AP_stackplot(squeeze(mua_psth_mean(curr_condition,:,:,1)),t,1.5,false,plot_color(curr_condition_idx,:));
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');

% Stack MUA plots 
figure('Name','Predicted'); hold on;
for curr_condition_idx = 1:length(plot_conditions);
    curr_condition = plot_conditions(curr_condition_idx);
    AP_stackplot(squeeze(mua_predicted_psth_mean(curr_condition,:,:,1)),t,1.5,false,plot_color(curr_condition_idx,:));
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');

% Stack fluor plots 
figure; hold on;
for curr_condition_idx = 1:length(plot_conditions)
    curr_condition = plot_conditions(curr_condition_idx);
    AP_stackplot(squeeze(fluor_psth_mean(curr_condition,:,:,1)),t,5,false,plot_color(curr_condition_idx,:),{wf_roi.area});
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');


%% Plot mean day-concatenated activity

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_df_kernel-str_' num2str(n_aligned_depths) '_depths.mat'];
load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = framerate*upsample_factor;
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);

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
mua_predicted_psth = nan(n_conditions,length(t),n_depths,2,length(D_all));
wheel_surround = nan(n_conditions,length(t),2,length(D_all));

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
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
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
    kernel_t = [-0.1,0.1];
    kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
    zs = [false,false];
    cvfold = 1;
    lambda = 0;
    
    mua_nonan_trials = ~squeeze(any(any(isnan(mua_cat_norm),2),4));
    mua_cat_predicted = nan(size(mua_cat_norm));
    for curr_depth = 1:n_depths
        curr_valid_trials = mua_nonan_trials(:,curr_depth);

        [~,predicted_spikes,~] = ...
            AP_regresskernel(reshape(permute( ...
            fluor_cat_norm(curr_valid_trials,:,:,:), ...
            [2,1,4,3]),[],n_rois)', ...
            reshape(permute(mua_cat_norm(curr_valid_trials,:,curr_depth,:), ...
            [2,1,4,3]),[],1)',kernel_frames,lambda,zs,cvfold);
        
        mua_cat_predicted(curr_valid_trials,:,curr_depth,:) = ...
            permute(reshape(predicted_spikes',length(t),sum(curr_valid_trials),2),[2,1,4,3]);
        
    end
    
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
    wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
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
    grp_fun = @(x) nanmean(x,1);
    
    for curr_roi = 1:n_rois
        fluor_psth(unique(trial_id),:,curr_roi,1,curr_animal) = grpstats(fluor_cat_hemidiff_norm(:,:,curr_roi,1),trial_id,grp_fun);
        fluor_psth(unique(trial_id),:,curr_roi,2,curr_animal) = grpstats(fluor_cat_hemidiff_norm(:,:,curr_roi,2),trial_id,grp_fun);
    end    
    
    for curr_depth = 1:n_depths
        mua_psth(unique(trial_id),:,curr_depth,1,curr_animal) = grpstats(mua_cat_norm(:,:,curr_depth,1),trial_id,grp_fun);
        mua_psth(unique(trial_id),:,curr_depth,2,curr_animal) = grpstats(mua_cat_norm(:,:,curr_depth,2),trial_id,grp_fun);
    end    
    
    for curr_depth = 1:n_depths
        mua_predicted_psth(unique(trial_id),:,curr_depth,1,curr_animal) = grpstats(mua_cat_predicted(:,:,curr_depth,1),trial_id,grp_fun);
        mua_predicted_psth(unique(trial_id),:,curr_depth,2,curr_animal) = grpstats(mua_cat_predicted(:,:,curr_depth,2),trial_id,grp_fun);
    end    
    
    wheel_surround(unique(trial_id),:,1,curr_animal) = grpstats(wheel_cat_norm(:,:,1),trial_id,grp_fun);
    wheel_surround(unique(trial_id),:,2,curr_animal) = grpstats(wheel_cat_norm(:,:,2),trial_id,grp_fun);
   
end

fluor_psth_mean = nanmean(fluor_psth,5);
mua_psth_mean = nanmean(mua_psth,5);
mua_predicted_psth_mean = nanmean(mua_predicted_psth,5);
wheel_surround_mean = diff(nanmean(wheel_surround,4),[],2);

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
% plot_conditions = find(conditions(:,1) == 0.06 & ...
%     conditions(:,2) == 1 & ...
%     conditions(:,4) == plot_timing);
% plot_color = zeros(length(plot_conditions),3);
% plot_color(conditions(plot_conditions,3) == -1,1) = 1;
% plot_color(conditions(plot_conditions,3) == 1,3) = 1;

% % (to plot zero-contrast by movement)
% plot_conditions = find(conditions(:,1) == 0 & ...
%     conditions(:,4) == plot_timing);
% plot_color = zeros(length(plot_conditions),3);
% plot_color(conditions(plot_conditions,3) == -1,1) = 1;
% plot_color(conditions(plot_conditions,3) == 1,3) = 1;

% Stack MUA plots 
figure; 
p1 = subplot(1,2,1); hold on;
for curr_condition_idx = 1:length(plot_conditions)
    curr_condition = plot_conditions(curr_condition_idx);
    AP_stackplot(squeeze(mua_psth_mean(curr_condition,:,:,1)),t,2,false,plot_color(curr_condition_idx,:), ...
        cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false));
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');

p2 = subplot(1,2,2); hold on;
for curr_condition_idx = 1:length(plot_conditions)
    curr_condition = plot_conditions(curr_condition_idx);
    AP_stackplot(squeeze(mua_psth_mean(curr_condition,:,:,2)),t,2,false,plot_color(curr_condition_idx,:), ...
        cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false));
end
line([0,0],ylim,'color','k');
xlabel('Time from move');

linkaxes([p1,p2]);

% Stack MUA predicted plots 
figure('Name','Predicted'); 
p1 = subplot(1,2,1); hold on;
for curr_condition_idx = 1:length(plot_conditions)
    curr_condition = plot_conditions(curr_condition_idx);
    AP_stackplot(squeeze(mua_predicted_psth_mean(curr_condition,:,:,1)),t,1.5,false,plot_color(curr_condition_idx,:), ...
        cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false));
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');

p2 = subplot(1,2,2); hold on;
for curr_condition_idx = 1:length(plot_conditions)
    curr_condition = plot_conditions(curr_condition_idx);
    AP_stackplot(squeeze(mua_predicted_psth_mean(curr_condition,:,:,2)),t,1.5,false,plot_color(curr_condition_idx,:), ...
        cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false));
end
line([0,0],ylim,'color','k');
xlabel('Time from move');

linkaxes([p1,p2]);

% Stack fluor plots 
figure; 
p1 = subplot(1,2,1); hold on;
for curr_condition_idx = 1:length(plot_conditions)
    curr_condition = plot_conditions(curr_condition_idx);
    AP_stackplot(squeeze(fluor_psth_mean(curr_condition,:,:,1)),t,5,false,plot_color(curr_condition_idx,:),{wf_roi.area});
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');

p2 = subplot(1,2,2); hold on;
for curr_condition_idx = 1:length(plot_conditions)
    curr_condition = plot_conditions(curr_condition_idx);
    AP_stackplot(squeeze(fluor_psth_mean(curr_condition,:,:,2)),t,5,false,plot_color(curr_condition_idx,:),{wf_roi.area});
end
line([0,0],ylim,'color','k');
xlabel('Time from move');

linkaxes([p1,p2]);

% Stack fluor plots (unilateral)
figure; 
p1 = subplot(1,2,1); hold on;
for curr_condition_idx = 1:length(plot_conditions)
    curr_condition = plot_conditions(curr_condition_idx);
    AP_stackplot(squeeze(fluor_psth_mean(curr_condition,:,1:n_rois/2,1)),t,3,false,plot_color(curr_condition_idx,:),{wf_roi(:,1).area});
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');

p2 = subplot(1,2,2); hold on;
for curr_condition_idx = 1:length(plot_conditions)
    curr_condition = plot_conditions(curr_condition_idx);
    AP_stackplot(squeeze(fluor_psth_mean(curr_condition,:,1:n_rois/2,2)),t,3,false,plot_color(curr_condition_idx,:),{wf_roi(:,1).area});
end
line([0,0],ylim,'color','k');
xlabel('Time from move');

linkaxes([p1,p2]);

% Plot wheel
t_diff = conv(t,[1,1]/2,'valid');
figure;
p1 = subplot(1,2,1); hold on;
for curr_condition_idx = 1:length(plot_conditions)
    curr_condition = plot_conditions(curr_condition_idx);
    plot(t_diff,wheel_surround_mean(curr_condition,:,1),'color',plot_color(curr_condition_idx,:),'linewidth',2);
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');
ylabel('Wheel velocity');

p2 = subplot(1,2,2); hold on;
for curr_condition_idx = 1:length(plot_conditions)
    curr_condition = plot_conditions(curr_condition_idx);
    plot(t_diff,wheel_surround_mean(curr_condition,:,2),'color',plot_color(curr_condition_idx,:),'linewidth',2);
end
line([0,0],ylim,'color','k');
xlabel('Time from move');
ylabel('Wheel velocity');

linkaxes([p1,p2]);

%% Get all activity normalized and concatenated, sort by reaction time

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_df_kernel-str_' num2str(n_aligned_depths) '_depths.mat'];
load([data_path filesep data_fn]);

% Get time
framerate = 35;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = framerate*upsample_factor;
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);

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

% Concantenate D
D_cellcat = vertcat(D_all{:});
D_cellcat = struct2cell(vertcat(D_cellcat{:}));
D_allcat = cell2struct(arrayfun(@(x) vertcat(D_cellcat{x,:}),1:size(D_cellcat,1),'uni',false)',fieldnames(D_all{1}{1}));

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_rois,2);
fluor_unilateral_allcat = nan(0,length(t),n_rois,2);
mua_allcat = nan(0,length(t),n_depths,2);
wheel_allcat = nan(0,length(t),2);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

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
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
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
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
    wheel_cat_norm = wheel_cat;
%     % (to normalize wheel)
%     wheel_vel = diff(wheel_cat,[],2);
%     wheel_norm_factor = prctile(max(max(abs(wheel_vel),[],2),[],3),90);
%     wheel_vel_norm = (mat2gray(wheel_vel,[-wheel_norm_factor,wheel_norm_factor])-0.5)*2;
%     wheel_cat_norm = [zeros(size(wheel_cat,1),1,2),cumsum(wheel_vel_norm,2)];

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
    
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat_hemidiff_norm];
    fluor_unilateral_allcat = [fluor_unilateral_allcat;fluor_cat_norm];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat_norm];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];

end

% Concatenate unilateral 
fluor_unilateral_allcat = cat(1,fluor_unilateral_allcat(:,:,1:length(wf_roi),:), ...
    fluor_unilateral_allcat(:,:,length(wf_roi)+1:end,:));

% Sort trials by stim-to-move time
% (by position)
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 1,[],2);
% (by max speed)
% [~,move_idx] = max(abs(diff(wheel_allcat(:,:,1),[],2)),[],2);

move_t = t(move_idx);
[~,sort_idx] = sort(move_idx);

% Get maximum wheel velocity in chosen direction
wheel_velocity_allcat = [zeros(size(wheel_allcat,1),1,2),diff(wheel_allcat,[],2)];
[max_speed,max_vel_idx] = max(abs(wheel_velocity_allcat(:,:,1).* ...
    (bsxfun(@times,wheel_velocity_allcat(:,:,1),trial_choice_allcat) > 0)),[],2);
vel_t = t(max_vel_idx);

max_vel = max_speed.*trial_choice_allcat;

% Set trials to plot
correct_trials = trial_contrast_allcat > 0 & trial_side_allcat == -trial_choice_allcat;
incorrect_trials = trial_contrast_allcat > 0 & trial_side_allcat == trial_choice_allcat;
zero_l_trials = trial_contrast_allcat == 0 & trial_choice_allcat == -1;
zero_r_trials = trial_contrast_allcat == 0 & trial_choice_allcat == 1;

% Plot wheel
plot_wheel = abs(diff(wheel_allcat(:,:,1),[],2));
t_diff = conv(t,[1,1]/2,'valid');
figure; colormap(gray);

subplot(1,4,1);
imagesc(t_diff,[],plot_wheel(sort_idx(correct_trials(sort_idx)),:));
line([0,0],ylim,'color','r');
line([0.5,0.5],ylim,'color','r');
ylabel('Sorted trial');
xlabel('Time')
title('Correct');

subplot(1,4,2);
imagesc(t_diff,[],plot_wheel(sort_idx(incorrect_trials(sort_idx)),:));
line([0,0],ylim,'color','r');
line([0.5,0.5],ylim,'color','r');
ylabel('Sorted trial');
xlabel('Time')
title('Incorrect');

subplot(1,4,3);
imagesc(t_diff,[],plot_wheel(sort_idx(zero_l_trials(sort_idx)),:));
line([0,0],ylim,'color','r');
line([0.5,0.5],ylim,'color','r');
ylabel('Sorted trial');
xlabel('Time')
title('Zero move L');

subplot(1,4,4);
imagesc(t_diff,[],plot_wheel(sort_idx(zero_r_trials(sort_idx)),:));
line([0,0],ylim,'color','r');
line([0.5,0.5],ylim,'color','r');
ylabel('Sorted trial');
xlabel('Time')
title('Zero move R');

% Plot activity
plot_act = mua_allcat(:,:,3,1);

% % (if concatenating side - messy, just here for the moment)
% trial_contrast_allcat = [trial_contrast_allcat;trial_contrast_allcat];
% trial_side_allcat = [trial_side_allcat;-trial_side_allcat];
% trial_choice_allcat = [trial_choice_allcat;-trial_choice_allcat];
% move_t = repmat(move_t,1,2);
% max_vel = [max_vel;-max_vel];

cscale = 5;

figure; colormap(gray);

subplot(1,4,1);
imagesc(t,[],plot_act(sort_idx(correct_trials(sort_idx)),:));
caxis([0,cscale])
line([0,0],ylim,'color','r');
line([0.5,0.5],ylim,'color','r');
title('Correct');
subplot(1,4,2);
imagesc(t,[],plot_act(sort_idx(incorrect_trials(sort_idx)),:));
caxis([0,cscale])
line([0,0],ylim,'color','r');
line([0.5,0.5],ylim,'color','r');
title('Incorrect');
subplot(1,4,3);
imagesc(t,[],plot_act(sort_idx(zero_l_trials(sort_idx)),:));
caxis([0,cscale])
line([0,0],ylim,'color','r');
line([0.5,0.5],ylim,'color','r');
title('Zero move L');
subplot(1,4,4);
imagesc(t,[],plot_act(sort_idx(zero_r_trials(sort_idx)),:));
caxis([0,cscale])
line([0,0],ylim,'color','r');
line([0.5,0.5],ylim,'color','r');
title('Zero move R');

% Plot activity by stim-to-move time
use_t = t > -0.1 & t < 0.1;
avg_act = nanmean(plot_act(:,use_t),2);

rxn_times = linspace(0,0.5,6);
rxn_times_centers = rxn_times(1:end-1)+diff(rxn_times)/2;
move_t_discretize = discretize(move_t,rxn_times);

use_contrast = 0.06;
move_l_trials = trial_choice_allcat == -1 & trial_side_allcat == 1 & trial_contrast_allcat == use_contrast;
move_r_trials = trial_choice_allcat == 1 & trial_side_allcat == 1 & trial_contrast_allcat == use_contrast;

[move_t_grp_move_l,avg_act_med_move_l,avg_act_mad_move_l] = ...
    grpstats(avg_act(move_l_trials),move_t_discretize(move_l_trials),{'gname','median','mad'});
[move_t_grp_move_r,avg_act_med_move_r,avg_act_mad_move_r] = ...
    grpstats(avg_act(move_r_trials),move_t_discretize(move_r_trials),{'gname','median','mad'});

move_t_grp_move_l = cellfun(@str2num,move_t_grp_move_l);
move_t_grp_move_r = cellfun(@str2num,move_t_grp_move_r);

figure; hold on
AP_errorfill(move_t_grp_move_l,avg_act_med_move_l,avg_act_mad_move_l,'k');
AP_errorfill(move_t_grp_move_r,avg_act_med_move_r,avg_act_mad_move_r,'r');
xlabel('Time from stim to move');
ylabel('Average activity');
legend({'Correct','Incorrect'});

% Plot activity by contrast/correct, restrict reaction time
use_t = t > -0.005 & t < 0.005;
use_stim_to_move = [move_t > 0.2 & move_t < 0.3]';
avg_act = nanmean(plot_act(:,use_t),2);

[move_t_grp_move_l,avg_act_med_move_l] = grpstats(avg_act(use_stim_to_move & trial_choice_allcat == -1), ...
    trial_side_allcat(use_stim_to_move & trial_choice_allcat == -1).* ...
    trial_contrast_allcat(use_stim_to_move & trial_choice_allcat == -1),{'gname','median'});
[move_t_grp_move_r,avg_act_med_move_r] = grpstats(avg_act(use_stim_to_move & trial_choice_allcat == 1), ...
    trial_side_allcat(use_stim_to_move & trial_choice_allcat == 1).* ...
    trial_contrast_allcat(use_stim_to_move & trial_choice_allcat == 1),{'gname','median'});

move_t_grp_move_l = cellfun(@str2num,move_t_grp_move_l);
move_t_grp_move_r = cellfun(@str2num,move_t_grp_move_r);

figure; hold on;
plot(move_t_grp_move_l,avg_act_med_move_l,'b','linewidth',2);
plot(move_t_grp_move_r,avg_act_med_move_r,'r','linewidth',2);
xlabel('Contrast*Side');
ylabel('Activity');

% Plot average activity within condition restricted by reaction time
plot_act = mua_allcat(:,:,1,1);

rxn_times = linspace(0,0.5,6);
n_rxn_times = length(rxn_times)-1;

activity_split = arrayfun(@(x) ...
    plot_act( ... 
    move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
    trial_side_allcat == 1 & ...
    trial_contrast_allcat == 1 & ...
    trial_choice_allcat == -1,:), ...
    1:length(rxn_times)-1,'uni',false);

activity_split_mean = cell2mat(cellfun(@(x) nanmean(x,1),activity_split','uni',false));

n_boot = 1000;
activity_boot = cellfun(@ (act) ...
    prctile(bootstrp(n_boot,@(x) nanmean(x,1),act),[50,2.5,97.5]), ...
    activity_split,'uni',false);

plot_col = copper(n_rxn_times);
figure; hold on;
p = arrayfun(@(x) AP_errorfill(t,activity_boot{x}(1,:)', ...
    bsxfun(@minus,activity_boot{x}(2:3,:),activity_boot{x}(1,:))', ...
    plot_col(x,:),0.5),1:n_rxn_times,'uni',false);

for i = 1:length(rxn_times)-1
   line(repmat(rxn_times(i),2,1),ylim,'color',plot_col(i,:));
end

legend([p{:}],arrayfun(@(x) [num2str(rxn_times(x)) '-' num2str(rxn_times(x+1))], ...
    1:length(rxn_times)-1,'uni',false));

% Plot activity difference of conditions
plot_act = mua_allcat(:,:,1,1);

rxn_times = linspace(0,0.5,6);
n_rxn_times = length(rxn_times)-1;

activity_split_0 = arrayfun(@(x) ...
    plot_act( ... 
    move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
    trial_contrast_allcat == 0 & ...
    trial_choice_allcat == -1,:), ...
    1:length(rxn_times)-1,'uni',false);

activity_split_6 = arrayfun(@(x) ...
    plot_act( ... 
    move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
    trial_side_allcat == 1 & ...
    trial_contrast_allcat == 0.06 & ...
    trial_choice_allcat == -1,:), ...
    1:length(rxn_times)-1,'uni',false);

n_boot = 1000;
activity_6_0 = cellfun(@ (x0,x6) ...
    prctile(bootstrp(n_boot,@(x) nanmean(x,1),x6) - ...
    bootstrp(n_boot,@(x) nanmean(x,1),x0),[50,2.5,97.5]), ...
    activity_split_0,activity_split_6,'uni',false);

activity_6_0_sig = ...
    cell2mat(cellfun(@(x) x(2,:) > 0,activity_6_0','uni',false));

plot_col = copper(n_rxn_times);
figure; hold on;
arrayfun(@(x) AP_errorfill(t,activity_6_0{x}(1,:)', ...
    bsxfun(@minus,activity_6_0{x}(2:3,:),activity_6_0{x}(1,:))', ...
    plot_col(x,:),0.5),1:n_rxn_times,'uni',false);

plot_size = fliplr(linspace(10,(n_rxn_times*10),n_rxn_times));
max_val = max(reshape([activity_6_0{:}],[],1));
for curr_rxn = 1:n_rxn_times
    curr_sig = activity_6_0_sig(curr_rxn,:);
    plot(t(curr_sig),repmat(max_val,sum(curr_sig),1),'.', ...
        'color',plot_col(curr_rxn,:),'MarkerSize',plot_size(curr_rxn));
end

% Plot all areas in different reaction time bins
plot_align = 1;

use_condition = ...
    trial_contrast_allcat == 0 & ...
    trial_choice_allcat == -1;

rxn_times = linspace(0,0.5,2);
n_rxn_times = length(rxn_times)-1;

fluor_split = arrayfun(@(x) ...
    fluor_unilateral_allcat( ... 
    move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
    use_condition,:,:,plot_align), ...
    1:length(rxn_times)-1,'uni',false);
fluor_split_mean = cell2mat(cellfun(@(x) nanmean(x,1),fluor_split','uni',false));

mua_split = arrayfun(@(x) ...
    mua_allcat( ... 
    move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
    use_condition,:,:,plot_align), ...
    1:length(rxn_times)-1,'uni',false);
mua_split_mean = cell2mat(cellfun(@(x) nanmean(x,1),mua_split','uni',false));

figure; 
for i = 1:n_rxn_times
    subplot(n_rxn_times,2,i*2-1); hold on;
    set(gca,'ColorOrder',jet(8));
    plot(t,permute(fluor_split_mean(i,:,:),[2,3,1]),'linewidth',2);
    ylim([min(fluor_split_mean(:)),max(fluor_split_mean(:))]);
    line([0,0],ylim,'color','k');
    
    subplot(n_rxn_times,2,i*2); hold on;
    set(gca,'ColorOrder',copper(4));
    plot(t,permute(mua_split_mean(i,:,:),[2,3,1]),'linewidth',2);
    ylim([min(mua_split_mean(:)),max(mua_split_mean(:))]);
    line([0,0],ylim,'color','k');
end

% Plot grid of all widefield areas in different reaction time bins
plot_align = 1;

rxn_times = linspace(0,0.5,4);
n_rxn_times = length(rxn_times)-1;

activity_split = cell(length(contrasts),n_rxn_times,8);
for curr_roi = 1:8
    for curr_contrast = 1:length(contrasts)
        
        if contrasts(curr_contrast) == 0
            use_sides = true(size(trial_side_allcat));
        else
            use_sides = trial_side_allcat == 1;
        end
        
        activity_split(curr_contrast,:,curr_roi) = ...
            arrayfun(@(x) fluor_unilateral_allcat( ...
            move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
            use_sides & ...
            trial_contrast_allcat == contrasts(curr_contrast) & ...
            trial_choice_allcat == -1,:,curr_roi,plot_align), ...
            1:n_rxn_times,'uni',false);

    end
end

activity_split_mean = cellfun(@(x) nanmean(x,1),activity_split,'uni',false);

figure;
for curr_roi = 1:8
    for curr_contrast = 1:length(contrasts)
        subplot(8,length(contrasts),length(contrasts)*(curr_roi-1)+curr_contrast);
        hold on; set(gca,'ColorOrder',copper(n_rxn_times));
        curr_data = vertcat(activity_split_mean{curr_contrast,:,curr_roi});
        plot(t,curr_data','linewidth',2);
        ylim([0,5])
        if curr_contrast == 1
            ylabel(wf_roi(curr_roi).area)
        end
    end
end

figure('Name','Zero-subtracted');
for curr_roi = 1:8
    for curr_contrast = 1:length(contrasts)
        subplot(8,length(contrasts),length(contrasts)*(curr_roi-1)+curr_contrast);
        hold on; set(gca,'ColorOrder',copper(n_rxn_times));
        curr_zero = vertcat(activity_split_mean{1,:,curr_roi});
        curr_data = vertcat(activity_split_mean{curr_contrast,:,curr_roi});
        plot(t,curr_data'-curr_zero','linewidth',2);
        
        if curr_contrast == 1
            ylabel(wf_roi(curr_roi).area)
        end
    end
end

% Plot grid of all MUA areas in different reaction time bins
plot_align = 1;

rxn_times = linspace(0,0.5,4);
n_rxn_times = length(rxn_times)-1;

activity_split = cell(length(contrasts),n_rxn_times,n_depths);
for curr_roi = 1:n_depths
    for curr_contrast = 1:length(contrasts)
        
        if contrasts(curr_contrast) == 0
            use_sides = true(size(trial_side_allcat));
        else
            use_sides = trial_side_allcat == 1;
        end
        
        activity_split(curr_contrast,:,curr_roi) = ...
            arrayfun(@(x) mua_allcat( ...
            move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
            use_sides & ...
            trial_contrast_allcat == contrasts(curr_contrast) & ...
            trial_choice_allcat == -1,:,curr_roi,plot_align), ...
            1:length(rxn_times)-1,'uni',false);

    end
end

activity_split_mean = cellfun(@(x) nanmean(x,1),activity_split,'uni',false);

figure;
for curr_roi = 1:n_depths
    for curr_contrast = 1:length(contrasts)
        subplot(n_depths,length(contrasts),length(contrasts)*(curr_roi-1)+curr_contrast);
        hold on; set(gca,'ColorOrder',copper(n_rxn_times));
        curr_data = vertcat(activity_split_mean{curr_contrast,:,curr_roi});
        plot(t,curr_data','linewidth',2);
        
        if curr_contrast == 1
            ylabel(['Str ' num2str(curr_roi)])
        end
    end
end

figure('Name','Zero-subtracted');
for curr_roi = 1:n_depths
    for curr_contrast = 1:length(contrasts)
        subplot(n_depths,length(contrasts),length(contrasts)*(curr_roi-1)+curr_contrast);
        hold on; set(gca,'ColorOrder',copper(n_rxn_times));
        curr_zero = vertcat(activity_split_mean{1,:,curr_roi});
        curr_data = vertcat(activity_split_mean{curr_contrast,:,curr_roi});
        plot(t,curr_data'-curr_zero','linewidth',2);
        
        if curr_contrast == 1
            ylabel(['Str ' num2str(curr_roi)])
        end
    end
end

% Plot different trial types overlapping 

% vis_trials = trial_contrast_allcat == 0.06 & ...
%     trial_side_allcat == 1 & ...
%     trial_choice_allcat == -1;
% 
% nonvis_trials = (trial_contrast_allcat == 0 & ...
%     trial_choice_allcat == -1) | ...
%     (trial_contrast_allcat > 0 & ...
%     trial_side_allcat == -1 & ...
%     trial_choice_allcat == -1);
% 
% trial_types = [vis_trials,nonvis_trials];

% use_rxn_time = move_t' > 0.5 & move_t' < 1;
% 
% vis_l_hit_trials = trial_contrast_allcat > 0 & ...
%     trial_side_allcat == 1 & ...
%     trial_choice_allcat == -1 & ...
%     use_rxn_time;
% 
% vis_l_miss_trials = trial_contrast_allcat > 0 & ...
%     trial_side_allcat == -1 & ...
%     trial_choice_allcat == -1 & ...
%     use_rxn_time;
% 
% zl_trials = trial_contrast_allcat == 0 & ...
%     trial_choice_allcat == -1 & ...
%     use_rxn_time;
% 
% trial_types = [vis_l_hit_trials,vis_l_miss_trials,zl_trials];

% rl_trials = trial_contrast_allcat == 0.25 & ...
%     trial_side_allcat == 1 & ...
%     trial_choice_allcat == -1;
% 
% rr_trials = trial_contrast_allcat == 0.25 & ...
%     trial_side_allcat == 1 & ...
%     trial_choice_allcat == 1;
% 
% ll_trials = trial_contrast_allcat == 0.25 & ...
%     trial_side_allcat == -1 & ...
%     trial_choice_allcat == -1;
% 
% lr_trials = trial_contrast_allcat == 0.25 & ...
%     trial_side_allcat == -1 & ...
%     trial_choice_allcat == 1;
% 
% trial_types = [rl_trials,rr_trials,ll_trials,lr_trials];

% zl_trials = trial_contrast_allcat == 0 & ...
%     trial_choice_allcat == -1;
% 
% zr_trials = trial_contrast_allcat == 0 & ...
%     trial_choice_allcat == 1;
% 
% trial_types = [zl_trials,zr_trials];


use_rxn_time = move_t' > 0.1 & move_t' < 0.3;

rl_trials = trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == -1 & ...
    use_rxn_time;

zl_trials = trial_contrast_allcat == 0 & ...
    trial_choice_allcat == -1 & ...
    use_rxn_time;

zr_trials = trial_contrast_allcat == 0 & ...
    trial_choice_allcat == 1 & ...
    use_rxn_time;

trial_types = [rl_trials,zl_trials,zr_trials];

% rl_trials = trial_contrast_allcat == 0.25 & ...
%     trial_side_allcat == 1 & ...
%     trial_choice_allcat == -1;
% 
% rr_trials = trial_contrast_allcat == 0.25 & ...
%     trial_side_allcat == 1 & ...
%     trial_choice_allcat == 1;
% 
% zl_trials = trial_contrast_allcat == 0 & ...
%     trial_choice_allcat == -1;
% 
% trial_types = [rl_trials,rr_trials,zl_trials];

fluor_split = ...
    arrayfun(@(x) fluor_unilateral_allcat( ...
    trial_types(:,x),:,:,:), ...
    1:size(trial_types,2),'uni',false);
fluor_split_mean = permute(cell2mat(cellfun(@(x) nanmean(x,1),fluor_split','uni',false)),[2,1,3,4]);

mua_split = ...
    arrayfun(@(x) mua_allcat( ...
    trial_types(:,x),:,:,:), ...
    1:size(trial_types,2),'uni',false);
mua_split_mean = permute(cell2mat(cellfun(@(x) nanmean(x,1),mua_split','uni',false)),[2,1,3,4]);

figure;
for curr_roi = 1:8
    for curr_align = 1:2
        subplot(8,2,curr_align+(curr_roi-1)*2); hold on;
        set(gca,'ColorOrder',copper(size(trial_types,2)));
        plot(t,fluor_split_mean(:,:,curr_roi,curr_align),'linewidth',2)
        ylabel(wf_roi(curr_roi).area);
    end
end

figure;
for curr_roi = 1:4
    for curr_align = 1:2
        subplot(4,2,curr_align+(curr_roi-1)*2); hold on;
        set(gca,'ColorOrder',copper(size(trial_types,2)));
        plot(t,mua_split_mean(:,:,curr_roi,curr_align),'linewidth',2)
        ylabel(['Str ' num2str(curr_roi)]);
    end
end

% Plot activity amplitude by wheel velocity
use_rxn = move_t' > 0 & move_t' < 0.5;

use_trials = trial_contrast_allcat > 0 & ...
    trial_side_allcat == -trial_choice_allcat & use_rxn;

use_data = mua_allcat(:,:,1,2);

% (fit response as mean)
fit_response = mat2gray(nanmean(abs(use_data(use_rxn,:)),1));
% (fit response as first SVD component)
% [rU,~,~] = svd(abs(use_data(use_rxn & ~any(isnan(use_data),2),:))');
% fit_response = mat2gray(U(:,1).*sign(nanmean(U(:,1))))';

activity_amplitude = fit_response'\use_data';
activity_reconstruct = transpose(fit_response'*activity_amplitude);

use_vel = max(abs(max_vel));
% vel_amp_edges = linspace(-use_vel,use_vel,11);
vel_amp_edges = prctile(abs(max_vel),linspace(0,100,5));
vel_amp_edges = sort([vel_amp_edges,-vel_amp_edges]);
vel_amp_bins = discretize(max_vel(use_trials),vel_amp_edges);

vel_amp_centers = grpstats(max_vel(use_trials),vel_amp_bins,{'nanmean'});
[amplitude_binned_mean,amplitude_binned_sem] = ...
    grpstats(activity_amplitude(use_trials),vel_amp_bins,{'nanmean','sem'});

% vel_act_edges = linspace(-use_vel,use_vel,7);
vel_act_edges = prctile(abs(max_vel),linspace(0,100,5));
vel_act_edges = sort([vel_act_edges,-vel_act_edges]);
vel_act_bins = discretize(max_vel(use_trials),vel_act_edges);
vel_act_centers = grpstats(max_vel(use_trials),vel_act_bins,{'nanmean'});
activity_binned_mean = grpstats(use_data(use_trials,:),vel_act_bins,{'nanmean'});

figure;
subplot(1,3,1); hold on;
col = colormap_BlueWhiteRed(floor(size(activity_binned_mean,1)/2));
set(gca,'ColorOrder',col);
plot(t,activity_binned_mean','linewidth',2);
line([0,0],ylim,'color','k');
xlabel('Time from move onset');
legend(cellfun(@num2str,num2cell(vel_act_centers),'uni',false));
axis tight

subplot(1,3,2); hold on;
plot(t,nanmean(use_data(use_trials,:),1),'k','linewidth',2);
plot(t,nanmean(abs(use_data(use_trials,:)-activity_reconstruct(use_trials,:)),1),'r','linewidth',2)
line([0,0],ylim,'color','b');
xlabel('Time from move onset');
legend({'Mean activity','Mean abs residual'});
axis tight

subplot(1,3,3);
errorbar(vel_amp_centers,amplitude_binned_mean,amplitude_binned_sem,'k','linewidth',2);
line([0,0],ylim,'color','r');
xlabel('Wheel velocity');ylabel('Response amplitude');
axis tight

% ANOVAN between contrast and reaction time
use_trials = trial_contrast_allcat > 0 & trial_side_allcat == 1 & trial_choice_allcat == -1;

activity_cat = cat(3,fluor_allcat(:,:,1:n_rois/2,:),mua_allcat);
area_labels = [{wf_roi(:,1).area}, ...
    cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false)];

activity_p = nan(3,size(activity_cat,3));

for curr_act = 1:size(activity_cat,3)
    use_data = activity_cat(:,:,curr_act,2);
    fit_response = mat2gray(nanmean(use_data,1));
    activity_amplitude = fit_response'\use_data';
    
    [p,tbl,stats,terms] = anovan(activity_amplitude(use_trials), ...
        [trial_contrast_allcat(use_trials),move_t(use_trials)'], ...
        'continuous',1:2,'model','full','display','off');
    activity_p(:,curr_act) = p;
end
figure;imagesc(activity_p < 0.05);
colormap(gray);
set(gca,'YTick',1:3,'YTickLabel',{'Contrast','Response time','Interaction'});
set(gca,'XTick',1:size(activity_p,2),'XTickLabel',area_labels);
axis image;
title('ANOVAN p < 0.05');

% Correlation of fitted amplitudes across areas
use_rxn = move_t' > 0.1 & move_t' < 0.3;
template_trials = use_rxn;

use_trials = use_rxn;

activity_cat = cat(3,fluor_allcat(:,:,:,:),mua_allcat);
area_labels = [{wf_roi(:).area}, ...
    cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false)];

activity_amplitude = nan(size(activity_cat,1),size(activity_cat,3),size(activity_cat,4));

for curr_act = 1:size(activity_cat,3)
    use_data = activity_cat(:,:,curr_act,:);
    fit_response = squeeze(mat2gray(nanmean(abs(use_data(template_trials,:,:)),1)));
    for curr_align = 1:size(activity_cat,4)
        activity_amplitude(:,curr_act,curr_align) = ...
            fit_response(:,curr_align)\ ...
            use_data(:,:,:,curr_align)';   
    end    
end

nonan_trials = ~any(any(isnan(activity_amplitude),2),3);
amp_corr = corrcoef(reshape(activity_amplitude(nonan_trials & use_trials,:,:), ...
    sum(nonan_trials & use_trials),[]));
figure;imagesc(amp_corr);
axis image;
colormap(brewermap([],'*RdBu'));
caxis([-max(abs(amp_corr(:))),max(abs(amp_corr(:)))]);
set(gca,'YTick',1:length(amp_corr),'YTickLabel',area_labels)
set(gca,'XTick',1:length(amp_corr),'XTickLabel',area_labels,'XTickLabelRotation',90)


%% Inter-area correlations on concatenated data

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_df_kernel-str_' num2str(n_aligned_depths) '_depths.mat'];
load([data_path filesep data_fn]);

% Get time
framerate = 35;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = framerate*upsample_factor;
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);

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

% Concantenate D
D_cellcat = vertcat(D_all{:});
D_cellcat = struct2cell(vertcat(D_cellcat{:}));
D_allcat = cell2struct(arrayfun(@(x) vertcat(D_cellcat{x,:}),1:size(D_cellcat,1),'uni',false)',fieldnames(D_all{1}{1}));

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_rois,2);
fluor_unilateral_allcat = nan(0,length(t),n_rois,2);
mua_allcat = nan(0,length(t),n_depths,2);
wheel_allcat = nan(0,length(t),2);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

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
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
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
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
    wheel_cat_norm = wheel_cat;
%     % (to normalize wheel)
%     wheel_vel = diff(wheel_cat,[],2);
%     wheel_norm_factor = prctile(max(max(abs(wheel_vel),[],2),[],3),90);
%     wheel_vel_norm = (mat2gray(wheel_vel,[-wheel_norm_factor,wheel_norm_factor])-0.5)*2;
%     wheel_cat_norm = [zeros(size(wheel_cat,1),1,2),cumsum(wheel_vel_norm,2)];

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
    
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat_hemidiff_norm];
    fluor_unilateral_allcat = [fluor_unilateral_allcat;fluor_cat_norm];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat_norm];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];

end

% Sort trials by stim-to-move time
% (by position)
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 1,[],2);
% (by max speed)
% [~,move_idx] = max(abs(diff(wheel_allcat(:,:,1),[],2)),[],2);

move_t = t(move_idx);
[~,sort_idx] = sort(move_idx);

% Get maximum wheel velocity in chosen direction
wheel_velocity_allcat = [zeros(size(wheel_allcat,1),1,2),diff(wheel_allcat,[],2)];
[max_speed,max_vel_idx] = max(abs(wheel_velocity_allcat(:,:,1).* ...
    (bsxfun(@times,wheel_velocity_allcat(:,:,1),trial_choice_allcat) > 0)),[],2);
vel_t = t(max_vel_idx);

max_vel = max_speed.*trial_choice_allcat;

% Concatenate all activity
activity_cat = cat(3,fluor_unilateral_allcat(:,:,1:n_rois/2,:),mua_allcat);
area_labels = [{wf_roi(:,1).area}, ...
    cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false)];

% Correlate all pairs of areas
use_trials = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == -1 & ...
    move_t' > 0.2 & move_t' < 0.3;

use_align = 1;

corr_grid = cell(size(activity_cat,3));
for curr_area_1 = 1:size(activity_cat,3)
    for curr_area_2 = 1:curr_area_1-1
         
        nonan_trials = ~any(any(isnan(activity_cat(:,:,curr_area_1,:)),2),4) & ...
            ~any(any(isnan(activity_cat(:,:,curr_area_2,:)),2),4);
        
        compare_act_corr = 1-pdist2(activity_cat(nonan_trials & use_trials,:,curr_area_1,use_align)', ...
            activity_cat(nonan_trials & use_trials,:,curr_area_2,use_align)','correlation');
        
        corr_grid{curr_area_1,curr_area_2} = compare_act_corr;                
        
    end
    
    AP_print_progress_fraction(curr_area_1,size(activity_cat,3));
    
end

area_labels_grid = cellfun(@(x,y) [x '-' y], ...
    repmat(area_labels',1,length(area_labels)), ...
    repmat(area_labels,length(area_labels),1), ...
    'uni',false);
area_labels_grid_used = area_labels_grid(~cellfun(@isempty,corr_grid));

corr_grid_cat = cat(3,corr_grid{:});

AP_imscroll(corr_grid_cat,area_labels_grid_used);
axis image;
colormap(brewermap([],'*RdBu'));
caxis([-0.5,0.5])
line([54,54],ylim,'color','k');
line(xlim,[54,54],'color','k');
line(xlim,ylim,'color','k');

% Get zero-lag correlation
corr_diag = cell2mat(arrayfun(@(x) diag(corr_grid_cat(:,:,x)),1:size(corr_grid_cat,3),'uni',false))';
corr_diag_norm = bsxfun(@rdivide,corr_diag,max(abs(corr_diag),[],2));

% % Sort by max time
% [~,max_idx] = max(corr_diag,[],2);
% [~,sort_idx] = sort(max_idx);
% Sort by PC1
[coeff,score,latent] = pca(corr_diag_norm);
[~,sort_idx] = sort(score(:,1));

figure;
imagesc(t,[],corr_diag_norm(sort_idx,:));
caxis([-1,1]);
colormap(brewermap([],'*RdBu'));
set(gca,'YTick',1:size(corr_diag,1),'YTickLabel',area_labels_grid_used(sort_idx));
line([0,0],ylim,'color','k');
title('Zero-lag correlations')

% Get local-lag correlation
max_lag = 0.05;
max_lag_samples = round(max_lag*sample_rate);

corr_diag_local = squeeze(max(cell2mat(permute(arrayfun(@(y) cell2mat(arrayfun(@(x) padarray(diag(corr_grid_cat(:,:,y),x), ...
    [abs(x),0],NaN,'post'),-max_lag_samples:max_lag_samples,'uni',false)),1:size(corr_grid_cat,3),'uni',false),[1,3,2])),[],2))';
corr_diag_norm = bsxfun(@rdivide,corr_diag_local,max(abs(corr_diag_local),[],2));

% Sort by PC1
[coeff,score,latent] = pca(corr_diag_norm);
[~,sort_idx] = sort(score(:,1));

figure;
imagesc(t,[],corr_diag_norm(sort_idx,:));
caxis([-1,1]);
colormap(brewermap([],'*RdBu'));
set(gca,'YTick',1:size(corr_diag,1),'YTickLabel',area_labels_grid_used(sort_idx));
line([0,0],ylim,'color','k');
title('Local-lag correlations')

% Get time-local forward-back correlation
max_lag = 0.08;
max_lag_samples = round(max_lag*sample_rate);

corr_grid_local = cell2mat(permute(arrayfun(@(y) cell2mat(arrayfun(@(x) padarray(diag(corr_grid_cat(:,:,y),x), ...
    [abs(x),0],NaN,'post'),-max_lag_samples:max_lag_samples,'uni',false)),1:size(corr_grid_cat,3),'uni',false),[1,3,2]));

corr_grid_local_diff = squeeze(nanmean(corr_grid_local(:,max_lag_samples+2:end,:),2) - ...
    nanmean(corr_grid_local(:,max_lag_samples:-1:1,:),2))';
corr_grid_local_diff_norm = bsxfun(@rdivide,corr_grid_local_diff,max(abs(corr_grid_local_diff),[],2));

[~,sort_idx] = sort(max(corr_grid_local_diff,[],2));

figure;
imagesc(t,[],corr_grid_local_diff(sort_idx,:));
caxis([-max(abs(corr_grid_local_diff(:))),max(abs(corr_grid_local_diff(:)))]);
colormap(brewermap([],'*RdBu'));
set(gca,'YTick',1:size(corr_diag,1),'YTickLabel',area_labels_grid_used(sort_idx));
line([0,0],ylim,'color','k');
title('Forward-back correlations')


% Plot correlation vs activity mean (to see if one is just a function of
% the other)
mean_act = squeeze(nanmean(activity_cat(use_trials,:,:,use_align),1));
mean_act_geomean = sqrt(abs(bsxfun(@times,permute(mean_act,[2,3,1]),permute(mean_act,[3,2,1]))));

mean_act_geomean_pairs = ...
    cell2mat(arrayfun(@(x) AP_itril(mean_act_geomean(:,:,x),-1),1:length(t),'uni',false));



corr_mean_diff = bsxfun(@rdivide,corr_diag,max(abs(corr_diag),[],2)) - ...
    bsxfun(@rdivide,mean_act_geomean_pairs,max(abs(mean_act_geomean_pairs),[],2));
figure;imagesc(t,[],corr_mean_diff(sort_idx,:));
caxis([-1,1]);
colormap(brewermap([],'*RdBu'));
set(gca,'YTick',1:size(corr_diag,1),'YTickLabel',area_labels_grid_used(sort_idx));
line([0,0],ylim,'color','k');



compare_areas = [1,2];
curr_geo_mean = squeeze(mean_act_geomean(compare_areas(2),compare_areas(1),:));
curr_corr = diag(corr_grid{compare_areas(2),compare_areas(1)});

figure; 
subplot(1,2,1); hold on;
plot(t,mean_act(:,compare_areas(1))/max(mean_act(:,compare_areas(1))),'linewidth',2);
plot(t,mean_act(:,compare_areas(2))/max(mean_act(:,compare_areas(2))),'linewidth',2);
plot(t,curr_corr,'linewidth',2);
xlabel('Time');
ylabel('Max-normalized')
legend([area_labels(compare_areas),{'Corr'}]);

subplot(1,2,2);hold on
t0 = find(t > 0,1);
t1 = find(t > 0.1,1);
plot(curr_geo_mean,curr_corr,'k','linewidth',2)
plot(curr_geo_mean(t0),curr_corr(t0),'ok','MarkerSize',10)
plot(curr_geo_mean(t1),curr_corr(t1),'^k','MarkerSize',10)
xlabel('Geometric mean activity');
ylabel('Correlation');



% Plot ctx/str average correlations
ctx_str_binary = [zeros(size(activity_cat,3)-4,1);ones(4,1)];
ctx_str_used = AP_itril(bsxfun(@plus,ctx_str_binary,ctx_str_binary'),-1);

ctx_ctx_corr_diag = corr_diag_local(ctx_str_used == 0,:);
ctx_str_corr_diag = corr_diag_local(ctx_str_used == 1,:);
str_str_corr_diag = corr_diag_local(ctx_str_used == 2,:);

figure;
subplot(4,3,[1,4,7]);
[coeff,score,latent] = pca(ctx_ctx_corr_diag);
[~,sort_idx] = sort(score(:,1));
imagesc(t,[],ctx_ctx_corr_diag(sort_idx,:));
caxis([-1,1]);
colormap(brewermap([],'*RdBu'));
line([0,0],ylim,'color','k');

subplot(4,3,[2,5,8]);
[coeff,score,latent] = pca(ctx_str_corr_diag);
[~,sort_idx] = sort(score(:,1));
imagesc(t,[],ctx_str_corr_diag(sort_idx,:));
caxis([-1,1]);
colormap(brewermap([],'*RdBu'));
line([0,0],ylim,'color','k');

subplot(4,3,[3,6,9]);
[coeff,score,latent] = pca(str_str_corr_diag);
[~,sort_idx] = sort(score(:,1));
imagesc(t,[],str_str_corr_diag(sort_idx,:));
caxis([-1,1]);
colormap(brewermap([],'*RdBu'));
line([0,0],ylim,'color','k');

subplot(4,3,[10,11,12]); hold on
plot(t,nanmean(ctx_ctx_corr_diag,1),'linewidth',2);
plot(t,nanmean(ctx_str_corr_diag,1),'linewidth',2);
plot(t,nanmean(str_str_corr_diag,1),'linewidth',2);
line([0,0],ylim,'color','k');
legend({'Ctx-Ctx','Ctx-Str','Str-Str'})


% This doesn't really look different from above

% % Get sig correlations across all time points across modalities
% % (SHUFFLE ONLY WITHIN CONDITION!)
% 
% nonan_trials = ~any(any(any(isnan(activity_cat),2),4),3);
%    
% use_trials = move_t' > 0 & move_t' < 0.5 & nonan_trials;
% 
% trial_sidecontrast_allcat = trial_side_allcat.*trial_contrast_allcat;
% 
% n_shuff = 1000;
% warning off;
% use_conditions = unique(trial_sidecontrast_allcat(use_trials));
% data_idx = reshape(1:sum(use_trials)*length(t), ...
%     sum(use_trials),length(t));
% shuff_idx = nan(sum(use_trials),length(t),n_shuff);
% for curr_condition_idx = 1:length(use_conditions)
%     curr_condition = use_conditions(curr_condition_idx);
%     curr_trials = trial_sidecontrast_allcat(use_trials) == curr_condition;
%     shuff_idx(curr_trials,:,:) = ...
%         AP_shake(repmat(data_idx(curr_trials,:,:),1,1,n_shuff),1);
%     AP_print_progress_fraction(curr_condition_idx,length(use_conditions));
% end
% warning on ;
% 
% % correlations and shuffle
% use_align = 2;
% act_corr = cell(size(activity_cat,3),size(activity_cat,3));
% for curr_area_1 = 1:size(activity_cat,3)
%     for curr_area_2 = 1:curr_area_1-1
%               
%         curr_data1 = zscore(activity_cat(use_trials,:,curr_area_1,use_align),[],1);
%         curr_data2 = zscore(activity_cat(use_trials,:,curr_area_2,use_align),[],1);
%         curr_data2_shuff = curr_data2(shuff_idx);
%         
%         corr_real = (curr_data1'*curr_data2)./(sum(use_trials)-1);
%         
%         corr_shuff = ...
%             gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
%             gpuArray(curr_data2_shuff)))./(sum(use_trials)-1);
%         
%         corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
%         corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
%         
%         act_corr{curr_area_1,curr_area_2} = corr_real;
%         act_corr{curr_area_1,curr_area_2}(~corr_sig) = 0;
%         
%     end
%     AP_print_progress_fraction(curr_area_2,size(activity_cat,3));
% end
% 
% area_labels_grid = cellfun(@(x,y) [x '-' y], ...
%     repmat(area_labels',1,length(area_labels)), ...
%     repmat(area_labels,length(area_labels),1), ...
%     'uni',false);
% area_labels_grid_used = area_labels_grid(~cellfun(@isempty,act_corr));
% 
% a = cat(3,act_corr{:});
% 
% AP_imscroll(a,area_labels_grid_used);
% axis image;
% colormap(brewermap([],'*RdBu'));
% caxis([-0.5,0.5])
% line([54,54],ylim,'color','k');
% line(xlim,[54,54],'color','k');
% line(xlim,ylim,'color','k');


%% Get concatenated activity in V-space

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_Udf_kernel-str_4_depths'];
% data_fn = 'all_trial_activity_df_msn_earlymove.mat';

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);


% Get widefield ROIs
n_rois = 200;

% MUA depths
n_depths = size(mua_all{1}{1},3);

% Set up conditions
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_rois,2);
mua_allcat = nan(0,length(t),n_depths,2);
wheel_allcat = nan(0,length(t),2);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
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
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
%     wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
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
    
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];
   
end

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Settings to plot
use_rxn_time = move_t > 0.1 & move_t < 0.3;
plot_align = 2;
normalize_px = true;

% Get major trial types
vis_L_trials_hit = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == -1 & ...
    use_rxn_time;

vis_R_trials_hit = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == -1 & ...
    trial_choice_allcat == 1 & ...
    use_rxn_time;

vis_L_trials_miss = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == -1 & ...
    trial_choice_allcat == -1 & ...
    use_rxn_time;

vis_R_trials_miss = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == 1 & ...
    use_rxn_time;

zero_L_trials = ...
    trial_contrast_allcat == 0 & ...
    trial_choice_allcat == -1 & ...
    use_rxn_time;

zero_R_trials = ...
    trial_contrast_allcat == 0 & ...
    trial_choice_allcat == 1 & ...
    use_rxn_time;

trial_types = ...
    [vis_L_trials_hit, vis_R_trials_hit, ...
    vis_L_trials_miss, vis_R_trials_miss, ...
    zero_L_trials, zero_R_trials];

%%% Get and plot fluorescence

px_trial_types = nan(size(U_master,1),size(U_master,2), ...
    size(fluor_allcat,2)-1,size(trial_types,2));
for curr_trial_type = 1:size(trial_types,2)
    
%     % Straight fluorescence
%     curr_data = fluor_allcat(trial_types(:,curr_trial_type),:,:,plot_align);
%     curr_data = squeeze(nanmean(bsxfun(@minus,curr_data,nanmean(curr_data(:,t > -0.2 & t < 0,:),2)),1))';
%     curr_px = svdFrameReconstruct(U_master(:,:,1:n_rois),curr_data(:,1:end-1));
    
%     % Fluorescence derivative
%     curr_data = squeeze(nanmean(fluor_allcat(trial_types(:,curr_trial_type),:,:,plot_align),1))';
%     curr_px = diff(svdFrameReconstruct(U_master(:,:,1:n_rois),curr_data),[],3);
%     curr_px(curr_px < 0) = 0;
    
    % Smoothed Fluorescence derivative       
    smooth_factor = 10;
    curr_data = squeeze(nanmean(convn(fluor_allcat(trial_types(:,curr_trial_type),:,:,plot_align), ...
        ones(1,smooth_factor)/smooth_factor,'same'),1))';
    curr_px = diff(svdFrameReconstruct(U_master(:,:,1:n_rois),curr_data),[],3);
    curr_px(curr_px < 0) = 0;

    px_trial_types(:,:,:,curr_trial_type) = curr_px;

end

t_diff = conv(t,[1,1]/2,'valid');

% Flip and combine trial types
% 1) visual hit (contra), 2) visual miss (ipsi), 3) zero (L)
px_combined = cat(4, ...
    (px_trial_types(:,:,:,1) + AP_reflect_widefield(px_trial_types(:,:,:,2)))./2, ...
    (px_trial_types(:,:,:,3) + AP_reflect_widefield(px_trial_types(:,:,:,4)))./2, ...
    (px_trial_types(:,:,:,5) + AP_reflect_widefield(px_trial_types(:,:,:,6)))./2);

px_combined_hemidiff = px_combined - AP_reflect_widefield(px_combined);

if normalize_px
    % Normalize by dividing by std of each frame
%     px_combined = bsxfun(@rdivide,px_combined,nanstd(reshape( ...
%         px_combined,[],1,size(px_combined,3),size(px_combined,4)),[],1));
%     px_combined_hemidiff = bsxfun(@rdivide,px_combined_hemidiff,nanstd(reshape( ...
%         px_combined_hemidiff,[],1,size(px_combined_hemidiff,3),size(px_combined_hemidiff,4)),[],1));

    % Normalize by dividing by max of each frame
    px_combined = bsxfun(@rdivide,px_combined,max(abs(reshape( ...
         px_combined,[],1,size(px_combined,3),size(px_combined,4))),[],1));
    px_combined_hemidiff = bsxfun(@rdivide,px_combined_hemidiff,max(abs(reshape( ...
         px_combined_hemidiff,[],1,size(px_combined_hemidiff,3),size(px_combined_hemidiff,4))),[],1));
end

% Plot all concatenated
AP_imscroll([reshape(permute(px_combined,[1,4,2,3]), ...
    [],size(px_combined,2),size(px_combined,3)), ...
    reshape(permute(px_combined_hemidiff,[1,4,2,3]), ...
    [],size(px_combined_hemidiff,2),size(px_combined_hemidiff,3))],t_diff);
axis image;
caxis([-prctile(abs(px_combined(:)),99),prctile(abs(px_combined(:)),99)]);
colormap(colormap_BlueWhiteRed);
% colormap(brewermap([],'RdBu'))

%%% Get and plot MUA (1-visual hit, 2-visual miss, 3-zero, uni and diff)

mua_trial_types = cell2mat(arrayfun(@(x) permute(nanmean(mua_allcat( ...
    trial_types(:,x),:,:,plot_align),1),[3,2,1]), ...
    permute(1:size(trial_types,2),[1,3,2]),'uni',false));

% 1) visual hit, 2) visual miss, 3) zero
col = copper(3);
figure; 

p1 = subplot(1,2,1); hold on;
AP_stackplot(mua_trial_types(:,:,1)',t,2,false,col(1,:));
AP_stackplot(mua_trial_types(:,:,3)',t,2,false,col(2,:));
AP_stackplot(mua_trial_types(:,:,5)',t,2,false,col(3,:));
title('Move L');

p2 = subplot(1,2,2); hold on;
AP_stackplot(mua_trial_types(:,:,1)' - mua_trial_types(:,:,2)',t,2,false,col(1,:));
AP_stackplot(mua_trial_types(:,:,3)' - mua_trial_types(:,:,4)',t,2,false,col(2,:));
AP_stackplot(mua_trial_types(:,:,5)' - mua_trial_types(:,:,6)',t,2,false,col(3,:));
title('Move L - Move R');

% Plot overlayed so that ROIs can be plotted over
figure; 
subplot(1,3,1); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot(t,mua_trial_types(:,:,1)','linewidth',2);
title ('Vis hit L');

subplot(1,3,2); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot(t,mua_trial_types(:,:,3)','linewidth',2);
title ('Vis miss L');

subplot(1,3,3); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot(t,mua_trial_types(:,:,5)','linewidth',2);
title ('Zero L');


%% Get concatenated activity in V-space (passive)

% Load data
n_aligned_depths = 4;
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_Udf_kernel-str_passive_fullscreen_4_depths'];
% data_fn = ['all_trial_activity_Udf_kernel-str_passive_choiceworld_4_depths_naive'];

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% % Load use experiments and cut out bad ones
% exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
% exclude_fn{1} = 'bhv_use_experiments';
% % exclude_fn{2} = 'expl_var_use_experiments';
% use_experiments_all = {};
% for curr_exclude = 1:length(exclude_fn)
%     curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
%     use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
% end
% use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
%     1:size(use_experiments_all,2),'uni',false);
% 
% D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
% fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
% mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
% wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);

% Get widefield ROIs
n_vs = 200;

% MUA depths
n_depths = size(mua_all{1}{1},3);

% Set up conditions
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_vs,1);
mua_allcat = nan(0,length(t),n_depths,1);
wheel_allcat = nan(0,length(t),1);
reward_allcat = nan(0,length(t),1);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

D_cat = struct('stimulus',cell(1,1),'response',cell(1,1),'repeatNum',cell(1,1));

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    n_align = size(fluor_cat,4);
    
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    t_std = t > -0.5 & t < 0.5;
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x(:,t_std,:,1),[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
    
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
    
    % Concatenate stim
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
    D_cat.stimulus = vertcat(D_cat.stimulus,D.stimulus);
   
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    
end

%%% REMOVE MOVE TRIALS FOR PASSIVE
move_trials = any(abs(wheel_allcat(:,t >= -0.5 & t <= 2) > 2),2);
fluor_allcat(move_trials,:,:) = [];
mua_allcat(move_trials,:,:) = [];
wheel_allcat(move_trials,:,:) = [];
D_cat.stimulus(move_trials) = [];
%%%

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

%%% Get and plot fluorescence
px_stim = nan(size(U_master,1),size(U_master,2), ...
    size(fluor_allcat,2)-1,length(unique(D_cat.stimulus)));
for curr_stim = unique(D_cat.stimulus)'
    
    curr_trials = D_cat.stimulus == curr_stim;
    
    %             % Straight fluorescence
    %             curr_data = fluor_allcat(curr_trials,:,:,plot_align);
    %             curr_baseline = nanmean(fluor_allcat(curr_trials,t > -0.2 & t < 0,:,1),2);
    %             curr_data = squeeze(nanmean(bsxfun(@minus,curr_data,curr_baseline),1))';
    %             curr_px = svdFrameReconstruct(U_master(:,:,1:n_rois),curr_data(:,1:end-1));
    
    %         % Fluorescence derivative
    %         curr_data = squeeze(nanmean(fluor_allcat(curr_trials,:,:),1))';
    %         curr_px = diff(svdFrameReconstruct(U_master(:,:,1:n_rois),curr_data),[],3);
    %         curr_px(curr_px < 0) = 0;
    
    % Smoothed Fluorescence derivative
    smooth_factor = 3;
    curr_data = convn(fluor_allcat(curr_trials,:,:), ...
        ones(1,smooth_factor)/smooth_factor,'same');
    
    curr_data_mean = squeeze(nanmean(curr_data,1))';   
    curr_px = diff(svdFrameReconstruct(U_master(:,:,1:n_vs),curr_data_mean),[],3);
%     curr_px(curr_px < 0) = 0;
    
    px_stim(:,:,:,curr_stim) = curr_px;
end

t_diff = conv(t,[1,1]/2,'valid');

% Choose stim to plot
plot_stim = unique(D_cat.stimulus);
% plot_stim = [2,5,8,10];

% Plot all together
AP_imscroll(reshape(permute(px_stim(:,:,:,plot_stim),[1,2,4,3]), ...
    size(px_stim,1),[],size(px_stim,3)),t_diff);
axis image; caxis([-1.5e-3,1.5e-3]); 
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k',[], ...
    [size(px_stim,1),size(px_stim,2),1,length(plot_stim)]);

% Plot striatum
plot_cols = lines(length(plot_stim));
figure; hold on;
p = gobjects(n_depths,length(plot_stim));
for curr_plot_condition = 1:length(plot_stim)  
    curr_trials = D_cat.stimulus == plot_stim(curr_plot_condition);
    curr_data_mean = squeeze(nanmean(mua_allcat(curr_trials,:,:),1));   
    curr_col = plot_cols(curr_plot_condition,:);   
    p(:,curr_plot_condition) = AP_stackplot(curr_data_mean,t,1,false,curr_col,1:n_depths);    
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');
legend(p(1,:),cellfun(@num2str,num2cell(unique(D_cat.stimulus)),'uni',false));

% Get fluorescence in widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);
n_rois = numel(wf_roi);

roi_mask = cat(3,wf_roi.mask);

U_roi = bsxfun(@rdivide,transpose(reshape(U_master,[],size(U_master,3))'* ...
    reshape(roi_mask,[],size(roi_mask,3))), ...
    sum(reshape(roi_mask,[],size(roi_mask,3)),1)');

smooth_size = 1;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

fluor_roi = permute(reshape(transpose(U_roi(:,1:n_vs)* ...
    reshape(permute(fluor_allcat(:,:,:,1),[3,2,1,4]),n_vs,[])), ...
    size(fluor_allcat,2),size(fluor_allcat,1),n_rois),[2,1,3]);

fluor_roi_diff = diff(padarray(convn(fluor_roi,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both'),[],2);

% (zero negatives, normalize)
fluor_roi_diff(fluor_roi_diff < 0) = 0;
fluor_roi_diff = bsxfun(@rdivide,fluor_roi_diff,nanstd(reshape(fluor_roi_diff,[],1,size(wf_roi,1)),[],1));

t_diff =  conv(t,[1,1]/2,'valid');

figure; hold on;
plot_cols = lines(length(plot_stim));
for curr_plot_condition = 1:length(plot_stim)
    
    curr_trials = D_cat.stimulus == plot_stim(curr_plot_condition);
    curr_data = fluor_roi_diff(curr_trials,:,:);
   
    curr_data_mean = squeeze(nanmean(curr_data,1));

    curr_col = plot_cols(curr_plot_condition,:);   
    AP_stackplot(curr_data_mean,t_diff,3.5,false,curr_col,{wf_roi.area});
    
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');
title('Widefield ROI');


%% Get concatenated activity in V-space in reaction time bins

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_Udf_kernel-str_4_depths_long'];

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);
reward_all = cellfun(@(data,use_expts) data(use_expts),reward_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);
reward_all = reward_all(use_animals);

% Get widefield ROIs
n_vs = 200;

% MUA depths
n_depths = size(mua_all{1}{1},3);

% Set up conditions
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_vs,1);
mua_allcat = nan(0,length(t),n_depths,1);
wheel_allcat = nan(0,length(t),1);
reward_allcat = nan(0,length(t),1);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

D_cat = struct('stimulus',cell(1,1),'response',cell(1,1),'repeatNum',cell(1,1));

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    n_align = size(fluor_cat,4);
    
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x,[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
%     wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
    % Concatenate behavioural data
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
    D.response = cell2mat(cellfun(@(x) x.response,D_all{curr_animal},'uni',false));
    D.repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all{curr_animal},'uni',false));
    
    D.day = trial_day;
    
    D_cat.stimulus = vertcat(D_cat.stimulus,D.stimulus);
    D_cat.response = vertcat(D_cat.response,D.response);
    D_cat.repeatNum = vertcat(D_cat.repeatNum,D.repeatNum);
    
    % Get trial ID   
    trial_contrast = max(D.stimulus,[],2);
    [~,side_idx] = max(D.stimulus > 0,[],2);
    trial_side = (side_idx-1.5)*2;
    trial_choice = -(D.response-1.5)*2;
    
    trial_conditions = ...
        [trial_contrast, trial_side, ...
        trial_choice, ones(size(trial_day))];
    [~,trial_id] = ismember(trial_conditions,conditions,'rows');
    
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    reward_allcat = [reward_allcat;vertcat(reward_all{curr_animal}{:})];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];
   
end

% Get max velocity in chosen direction
t_leeway = 0.5;
leeway_samples = round(t_leeway*sample_rate);

wheel_velocity_allcat = [zeros(size(wheel_allcat,1),1,n_align),diff(wheel_allcat,[],2)];

max_vel = max(abs(wheel_velocity_allcat(:,:,n_align).* ...
    (bsxfun(@times,wheel_velocity_allcat(:,:,n_align),trial_choice_allcat) > 0)),[],2).* ...
    trial_choice_allcat;

[max_speed,max_vel_idx] = max(abs(wheel_velocity_allcat(:,:,1).* ...
    (bsxfun(@times,wheel_velocity_allcat(:,:,1),trial_choice_allcat) > 0)),[],2);
vel_t = t(max_vel_idx);

max_vel = max_speed.*trial_choice_allcat;

% % To deconvolve, but this just looked like a super-smoothed derivative
% load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\gcamp_kernel\gcamp6s_kernel.mat');
% 
% gcamp6s_kernel_resample = interp1(gcamp6s_kernel_t,gcamp6s_kernel,gcamp6s_kernel_t(1):1/sample_rate:gcamp6s_kernel_t(end));
% gcamp6s_kernel_long = gcamp6s_kernel_resample(1:size(fluor_allcat,2));
% 
% softnorm = 10;
% fluor_allcat_deconv = permute(real(ifft(fft(permute(fluor_allcat,[2,1,3,4])).* ...
%     conj(fft(gcamp6s_kernel_long'))./(abs(fft(gcamp6s_kernel_long')).^2 + softnorm))),[2,1,3,4]);

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Settings to plot
% rxn_time_bins = {[0,0.1],[0.1,0.2],[0.2,0.3],[0.3,0.4],[0.6,0.7]};
% rxn_time_bins = {[0,0.15],[0.15,0.4],[0.6,0.7]};
% rxn_time_bins = {[0.1,0.3],[0.6,0.7]};
rxn_time_bins = {[0,0.5]};

move_align = false;
normalize_px = true;

% Get major trial types
vis_L_trials_hit = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == -1;

vis_R_trials_hit = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == -1 & ...
    trial_choice_allcat == 1;

vis_L_trials_miss = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == -1 & ...
    trial_choice_allcat == -1;

vis_R_trials_miss = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == 1;

zero_L_trials = ...
    trial_contrast_allcat == 0 & ...
    trial_choice_allcat == -1;

zero_R_trials = ...
    trial_contrast_allcat == 0 & ...
    trial_choice_allcat == 1;

trial_types = ...
    [vis_L_trials_hit, vis_R_trials_hit, ...
    vis_L_trials_miss, vis_R_trials_miss, ...
    zero_L_trials, zero_R_trials];

%%% Get and plot fluorescence

px_trial_types = nan(size(U_master,1),size(U_master,2), ...
    size(fluor_allcat,2)-1,size(trial_types,2),length(rxn_time_bins));
for curr_rxn = 1:length(rxn_time_bins)
    for curr_trial_type = 1:size(trial_types,2)
        
        curr_trials = trial_types(:,curr_trial_type) & ...
            move_t > rxn_time_bins{curr_rxn}(1) & ...
            move_t < rxn_time_bins{curr_rxn}(2);
        
        %             % Straight fluorescence
        %             curr_data = fluor_allcat(curr_trials,:,:,plot_align);
        %             curr_baseline = nanmean(fluor_allcat(curr_trials,t > -0.2 & t < 0,:,1),2);
        %             curr_data = squeeze(nanmean(bsxfun(@minus,curr_data,curr_baseline),1))';
        %             curr_px = svdFrameReconstruct(U_master(:,:,1:n_rois),curr_data(:,1:end-1));
        
        %         % Fluorescence derivative
        %         curr_data = squeeze(nanmean(fluor_allcat(curr_trials,:,:),1))';
        %         curr_px = diff(svdFrameReconstruct(U_master(:,:,1:n_rois),curr_data),[],3);
        %         curr_px(curr_px < 0) = 0;
        
        % Smoothed Fluorescence derivative
        smooth_factor = 3;
        curr_data = convn(fluor_allcat(curr_trials,:,:), ...
            ones(1,smooth_factor)/smooth_factor,'same');
        
        if move_align
            % Re-align to movement onset
            t_leeway = 0.5;
            leeway_samples = round(t_leeway*sample_rate);
            curr_move_idx = move_idx(curr_trials);
            for i = 1:size(curr_data,1)
                curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
            end
        end
        
        curr_data_mean = squeeze(nanmean(curr_data,1))';
        
        curr_px = diff(svdFrameReconstruct(U_master(:,:,1:n_vs),curr_data_mean),[],3);
%         curr_px(curr_px < 0) = 0;
        
        px_trial_types(:,:,:,curr_trial_type,curr_rxn) = curr_px;
    end
    AP_print_progress_fraction(curr_rxn,length(rxn_time_bins));
end

t_diff = conv(t,[1,1]/2,'valid');

% Flip and combine trial types
% 1) visual hit (contra), 2) visual miss (ipsi), 3) zero (L)
px_combined = cat(4, ...
    (px_trial_types(:,:,:,1,:) + AP_reflect_widefield(px_trial_types(:,:,:,2,:)))./2, ...
    (px_trial_types(:,:,:,3,:) + AP_reflect_widefield(px_trial_types(:,:,:,4,:)))./2, ...
    (px_trial_types(:,:,:,5,:) + AP_reflect_widefield(px_trial_types(:,:,:,6,:)))./2);

px_combined_hemidiff = px_combined - AP_reflect_widefield(px_combined);

if normalize_px
    % Normalize by dividing by max of each frame
    px_dims = size(px_combined);
    px_combined = bsxfun(@rdivide,px_combined,max(abs(reshape( ...
         px_combined,[px_dims(1)*px_dims(2),1,px_dims(3:end)])),[],1));
     
    px_combined_hemidiff = bsxfun(@rdivide,px_combined_hemidiff,max(abs(reshape( ...
         px_combined_hemidiff,[px_dims(1)*px_dims(2),1,px_dims(3:end)])),[],1));
end

% Plot
plot_rxn = 1;
AP_imscroll(cat(4,px_combined(:,:,:,:,plot_rxn),px_combined_hemidiff(:,:,:,:,plot_rxn)),t_diff);
axis image; caxis([-1,1]); 
% colormap(colormap_BlueWhiteRed);
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k');

% Plot all concatenated
px_dims = size(px_combined);
AP_imscroll([reshape(permute(px_combined,[1,4,2,5,3]), ...
    [],size(px_combined,2)*size(px_combined,5),size(px_combined,3)), ...
    reshape(permute(px_combined_hemidiff,[1,4,2,5,3]), ...
    [],size(px_combined_hemidiff,2)*size(px_combined_hemidiff,5), ...
    size(px_combined_hemidiff,3))],t_diff);
axis image;
% caxis([-prctile(abs(px_combined(:)),99),prctile(abs(px_combined(:)),99)]);
caxis([-1,1]);
% colormap(colormap_BlueWhiteRed);
colormap(brewermap([],'*RdBu'))


%%% Get and plot MUA 
mua_trial_types = nan(n_depths, ...
    size(mua_allcat,2),size(trial_types,2),length(rxn_time_bins));
for curr_rxn = 1:length(rxn_time_bins)
    for curr_trial_type = 1:size(trial_types,2)
        
        curr_trials = trial_types(:,curr_trial_type) & ...
            move_t > rxn_time_bins{curr_rxn}(1) & ...
            move_t < rxn_time_bins{curr_rxn}(2);
        
        curr_data = mua_allcat(curr_trials,:,:);
        
        if move_align
            % Re-align to movement onset
            t_leeway = 0.5;
            leeway_samples = round(t_leeway*sample_rate);
            curr_move_idx = move_idx(curr_trials);
            for i = 1:size(curr_data,1)
                curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
            end
        end
        
        mua_trial_types(:,:,curr_trial_type,curr_rxn) = squeeze(nanmean(curr_data,1))';
               
    end
    AP_print_progress_fraction(curr_rxn,length(rxn_time_bins));
end

% mua_trial_types_pseudodiff = cat(3,mua_trial_types(:,:,[1,3,5],:), ...
%     mua_trial_types(:,:,[1,3,5],:) -  mua_trial_types(:,:,[2,4,6],:));
mua_trial_types_pseudodiff = cat(3,mua_trial_types(:,:,[1,3,5],:), ...
    mua_trial_types(:,:,[2,4,6],:));

figure;
for curr_rxn = 1:length(rxn_time_bins)
    for curr_trial_type = 1:size(trial_types,2)/2
        
        subplot(size(trial_types,2)/2,length(rxn_time_bins)*2, ...
            curr_rxn + length(rxn_time_bins)*2*(curr_trial_type-1)); 
        hold on; set(gca,'ColorOrder',copper(n_depths));
        plot(t, ...
            mua_trial_types_pseudodiff(:,:,curr_trial_type,curr_rxn)', ...
            'linewidth',2);
        ylim([min(mua_trial_types(:)),max(mua_trial_types(:))]);
        line([0,0],ylim,'color','k')
        
        subplot(size(trial_types,2)/2,length(rxn_time_bins)*2, ...
            curr_rxn + length(rxn_time_bins) + length(rxn_time_bins)*2*(curr_trial_type-1)); 
        hold on; set(gca,'ColorOrder',copper(n_depths));
        plot(t, ...
            mua_trial_types_pseudodiff(:,:,curr_trial_type+3,curr_rxn)', ...
            'linewidth',2);
       ylim([min(mua_trial_types(:)),max(mua_trial_types(:))]);
        line([0,0],ylim,'color','k')
        
    end
end

%%% Get and plot wheel 
wheel_trial_types = nan(...
    size(wheel_velocity_allcat,2),size(trial_types,2),length(rxn_time_bins));
for curr_rxn = 1:length(rxn_time_bins)
    for curr_trial_type = 1:size(trial_types,2)
        
        curr_trials = trial_types(:,curr_trial_type) & ...
            move_t > rxn_time_bins{curr_rxn}(1) & ...
            move_t < rxn_time_bins{curr_rxn}(2);
        
        curr_data = wheel_velocity_allcat(curr_trials,:);
        if move_align
            % Re-align to movement onset
            t_leeway = 0.5;
            leeway_samples = round(t_leeway*sample_rate);
            curr_move_idx = move_idx(curr_trials);
            for i = 1:size(curr_data,1)
                curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
            end
        end
        wheel_trial_types(:,curr_trial_type,curr_rxn) = squeeze(nanmean(curr_data,1))';
               
    end
    AP_print_progress_fraction(curr_rxn,length(rxn_time_bins));
end

figure;
for curr_rxn = 1:length(rxn_time_bins)
    subplot(1,length(rxn_time_bins),curr_rxn); hold on;
    plot(t,wheel_trial_types(:,:,curr_rxn),'linewidth',2);
    ylim([-max(abs(wheel_trial_types(:))),max(abs(wheel_trial_types(:)))]);
    line([0,0],ylim,'color','k')
end

% Flip and combine trial types
% 1) visual hit (contra), 2) visual miss (ipsi), 3) zero (L)
wheel_combined = cat(2, ...
    (wheel_trial_types(:,1,:) - wheel_trial_types(:,2,:))./2, ...
    (wheel_trial_types(:,3,:) - wheel_trial_types(:,4,:))./2, ...
    (wheel_trial_types(:,5,:) - wheel_trial_types(:,6,:))./2);

figure;
for curr_rxn = 1:length(rxn_time_bins)
    for curr_trial_type = 1:size(trial_types,2)/2
        
        subplot(size(trial_types,2)/2,length(rxn_time_bins), ...
            curr_rxn + length(rxn_time_bins)*(curr_trial_type-1)); 
        plot(t, ...
            wheel_combined(:,curr_trial_type,curr_rxn)', ...
            'k','linewidth',2);
        ylim([-max(abs(wheel_trial_types(:))),max(abs(wheel_trial_types(:)))]);
        line([0,0],ylim,'color','k')

    end
end




% Get pixels by any conditions


% Get conditions for each trial, plot selected
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];

% % contrast, side, choice
% plot_conditions = ...
%     [contrasts(2:end),contrasts(2:end); ...
%     -ones(1,5),ones(1,5); ...
%     ones(1,5),-ones(1,5)]';

% plot_conditions = ...
%     [1,0,1,1,0,1; ...
%     1,-1,-1,-1,-1,1; ...
%     -1,-1,-1,1,1,1]';

plot_conditions = ...
    [1,0,1; ...
    1,-1,-1; ...
    -1,-1,-1]';

% plot_conditions = ...
%     [0.06,0.06,0.06; ...
%     1,1,-1; ...
%     -1,1,-1]';

% plot_conditions = ...
%     [1; ...
%     1; ...
%     -1]';

use_rxn = move_t > 0 & move_t < 0.5;
    
[~,plot_id] = ismember( ...
    [trial_contrast_allcat > 0,trial_side_allcat,trial_choice_allcat], ...
    plot_conditions,'rows');

move_align = true;

px_trial_types = nan(size(U_master,1),size(U_master,2), ...
    size(fluor_allcat,2)-1,size(plot_conditions,1));
for curr_plot_condition = 1:size(plot_conditions,1)
    curr_trials = plot_id == curr_plot_condition & use_rxn;
    
%     % Straight fluorescence
%     curr_data = fluor_allcat(curr_trials,:,:);
%     curr_baseline = nanmean(fluor_allcat(curr_trials,t > -0.2 & t < 0,:,1),2);
%     curr_data = bsxfun(@minus,curr_data(:,2:end,:),curr_baseline);
    
    % Fluorescence derivative
    curr_data = diff(fluor_allcat(curr_trials,:,:),[],2);

    if move_align
        % Re-align to movement onset
        t_leeway = 0.5;
        leeway_samples = round(t_leeway*sample_rate);
        curr_move_idx = move_idx(curr_trials);
        for i = 1:size(curr_data,1)
            curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
        end
    end
    
    curr_data_mean = squeeze(nanmean(curr_data,1))';
    curr_px = svdFrameReconstruct(U_master(:,:,1:n_vs),curr_data_mean);
    %         curr_px(curr_px < 0) = 0;
    
    px_trial_types(:,:,:,curr_plot_condition) = curr_px;
    AP_print_progress_fraction(curr_plot_condition,size(plot_conditions,1));
end

AP_imscroll(px_trial_types,t_diff);
axis image; caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k');



%% Load and plot long V data / ctx->str predicted

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_Udf_kernel-str_4_depths_long_predicted'];

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
mua_predicted_all = cellfun(@(data,use_expts) data(use_expts),mua_predicted_all,use_experiments','uni',false);
mua_predicted_nlin_all = cellfun(@(data,use_expts) data(use_expts),mua_predicted_nlin_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);
reward_all = cellfun(@(data,use_expts) data(use_expts),reward_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
mua_predicted_all = mua_predicted_all(use_animals);
mua_predicted_nlin_all = mua_predicted_nlin_all(use_animals);
wheel_all = wheel_all(use_animals);
reward_all = reward_all(use_animals);

% Get widefield ROIs
n_rois = 200;

% MUA depths
n_depths = size(mua_all{1}{1},3);

% Set up conditions
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_rois,1);
mua_allcat = nan(0,length(t),n_depths,1);
mua_predicted_allcat = nan(0,length(t),n_depths,1);
mua_predicted_nlin_allcat = nan(0,length(t),n_depths,1);
wheel_allcat = nan(0,length(t),1);
reward_allcat = nan(0,length(t),1);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

D_cat = struct('stimulus',cell(1,1),'response',cell(1,1),'repeatNum',cell(1,1));

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    mua_predicted_cat_raw = cat(1,mua_predicted_all{curr_animal}{:});
    mua_predicted_nlin_cat_raw = cat(1,mua_predicted_nlin_all{curr_animal}{:});
    
    n_align = size(fluor_cat,4);
    
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    mua_predicted_cat_raw = bsxfun(@times,mua_predicted_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    mua_predicted_nlin_cat_raw = bsxfun(@times,mua_predicted_nlin_cat_raw,permute(filled_mua_trials,[1,3,2,4]));   
    
    % Smooth and normalize MUA
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');

    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x,[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
 
    % Smooth and normalize MUA predicted (with MUA values)
    mua_predicted_cat_raw_smoothed = padarray(convn(mua_predicted_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    mua_predicted_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_predicted_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm); 
    
    % Smooth and normalize MUA predicted nonlinear (with MUA values)
    mua_predicted_nlin_cat_raw_smoothed = padarray(convn(mua_predicted_nlin_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both'); 
    mua_predicted_nlin_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_predicted_nlin_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm); 
        
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
%     wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
    % Concatenate behavioural data
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
    D.response = cell2mat(cellfun(@(x) x.response,D_all{curr_animal},'uni',false));
    D.repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all{curr_animal},'uni',false));
    
    D.day = trial_day;
    
    D_cat.stimulus = vertcat(D_cat.stimulus,D.stimulus);
    D_cat.response = vertcat(D_cat.response,D.response);
    D_cat.repeatNum = vertcat(D_cat.repeatNum,D.repeatNum);
    
    % Get trial ID   
    trial_contrast = max(D.stimulus,[],2);
    [~,side_idx] = max(D.stimulus > 0,[],2);
    trial_side = (side_idx-1.5)*2;
    trial_choice = -(D.response-1.5)*2;
    
    trial_conditions = ...
        [trial_contrast, trial_side, ...
        trial_choice, ones(size(trial_day))];
    [~,trial_id] = ismember(trial_conditions,conditions,'rows');
    
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat];
    mua_allcat = [mua_allcat;mua_cat_norm];
    mua_predicted_allcat = [mua_predicted_allcat;mua_predicted_cat_norm];
    mua_predicted_nlin_allcat = [mua_predicted_nlin_allcat;mua_predicted_nlin_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    reward_allcat = [reward_allcat;vertcat(reward_all{curr_animal}{:})];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];
   
end

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Get maximum wheel velocity in chosen direction
wheel_velocity_allcat = [zeros(size(wheel_allcat,1),1,1),diff(wheel_allcat,[],2)];
[max_speed,max_vel_idx] = max(abs(wheel_velocity_allcat(:,:,1).* ...
    (bsxfun(@times,wheel_velocity_allcat(:,:,1),trial_choice_allcat) > 0)),[],2);
vel_t = t(max_vel_idx);

max_vel = max_speed.*trial_choice_allcat;

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Get conditions for each trial, plot selected
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];

contrast_side_col = colormap_BlueWhiteRed(5);
contrast_side_col(6,:) = 0;
contrast_side_val = unique(sort([-contrasts,contrasts]))';

% contrast, side, choice
plot_conditions = ...
    [contrasts,contrasts; ...
    -ones(1,6),-1,ones(1,5); ...
    ones(1,6),-ones(1,6)]';

use_rxn = move_t > 0 & move_t < 0.5;

[~,plot_id] = ismember( ...
    [trial_contrast_allcat,trial_side_allcat,trial_choice_allcat], ...
    plot_conditions,'rows');

% Plot striatum
figure; hold on;
for curr_plot = 1:3
    
    switch curr_plot
        case 1
            plot_data = mua_allcat;
            plot_title = 'Measured';
        case 2
            plot_data = mua_predicted_nlin_allcat;
            plot_title = 'Predicted';
        case 3
            plot_data = mua_allcat - mua_predicted_nlin_allcat;
            plot_title = 'Residual';
    end
    
    p_str(curr_plot) = subplot(1,3,curr_plot); hold on;
    for curr_plot_condition = 1:size(plot_conditions,1)
        
        curr_trials = plot_id == curr_plot_condition & use_rxn;
        curr_data = plot_data(curr_trials,:,:);
        
        % re-align to movement onset
        t_leeway = 0.5;
        leeway_samples = round(t_leeway*(1/mean(diff(t))));
        curr_move_idx = move_idx(curr_trials);
        for i = 1:size(curr_data,1)
            curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
        end
        
        curr_data_mean = squeeze(nanmean(curr_data,1));
        
        curr_contrast_side = max(trial_contrast_allcat(curr_trials))*sign(max(trial_side_allcat(curr_trials)));
        contrast_side_idx = find(curr_contrast_side == contrast_side_val);
        curr_col = contrast_side_col(contrast_side_idx,:);
        
        if curr_contrast_side == 0
            switch max(trial_choice_allcat(curr_trials))
                case -1
                    curr_col = 'm';
                case 1
                    curr_col = 'c';
            end
        end
        
        AP_stackplot(curr_data_mean,t,2,false,curr_col,1:n_depths);
        
    end
    line([0,0],ylim,'color','k');
    xlabel('Time from stim');
    title(plot_title);
end
linkaxes(p_str)




%% Plot ROI activity from V-space data

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_Udf_kernel-str_4_depths_long'];
% data_fn = 'all_trial_activity_df_msn_earlymove.mat';

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);
reward_all = cellfun(@(data,use_expts) data(use_expts),reward_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);
reward_all = reward_all(use_animals);

% Get widefield components
n_vs = 200;

% MUA depths
n_depths = size(mua_all{1}{1},3);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_vs,1);
mua_allcat = nan(0,length(t),n_depths,1);
wheel_allcat = nan(0,length(t),1);
reward_allcat = nan(0,length(t),1);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    n_align = size(fluor_cat,4);
    
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x,[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
%     wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
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
 
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    reward_allcat = [reward_allcat;vertcat(reward_all{curr_animal}{:})];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];
   
end

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Get maximum wheel velocity in chosen direction
wheel_velocity_allcat = [zeros(size(wheel_allcat,1),1,1),diff(wheel_allcat,[],2)];
[max_speed,max_vel_idx] = max(abs(wheel_velocity_allcat(:,:,1).* ...
    (bsxfun(@times,wheel_velocity_allcat(:,:,1),trial_choice_allcat) > 0)),[],2);
vel_t = t(max_vel_idx);

max_vel = max_speed.*trial_choice_allcat;

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Low-pass filter data
lowpassCutoff = 10; % Hz
[b100s, a100s] = butter(2, lowpassCutoff/((sample_rate)/2), 'low');
fluor_allcat = filter(b100s,a100s,fluor_allcat,[],2);
mua_allcat = filter(b100s,a100s,mua_allcat,[],2);

% Get conditions for each trial, plot selected
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];

contrast_side_col = colormap_BlueWhiteRed(5);
contrast_side_col(6,:) = 0;
contrast_side_val = unique(sort([-contrasts,contrasts]))';

% % contrast, side, choice
% plot_conditions = ...
%     [contrasts(2:end),contrasts(2:end); ...
%     -ones(1,5),ones(1,5); ...
%     ones(1,5),-ones(1,5)]';

% plot_conditions = ...
%     [1,0,1,1,0,1; ...
%     1,-1,-1,-1,-1,1; ...
%     -1,-1,-1,1,1,1]';

% plot_conditions = ...
%     [1,0,1,1; ...
%     -1,-1,1,-1; ...
%     -1,-1,-1,1]';

plot_conditions = ...
    [1; ...
    1,; ...
    -1]';

% plot_conditions = ...
%     [0.06,0.06,0.06; ...
%     1,1,-1; ...
%     -1,1,-1]';

use_rxn = move_t > 0.5 & move_t < 1;
    
[~,plot_id] = ismember( ...
    [trial_contrast_allcat,trial_side_allcat,trial_choice_allcat], ...
    plot_conditions,'rows');

move_align = false;

% Get fluorescence in widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);
n_rois = numel(wf_roi);

roi_mask = cat(3,wf_roi.mask);

U_roi = bsxfun(@rdivide,transpose(reshape(U_master,[],size(U_master,3))'* ...
    reshape(roi_mask,[],size(roi_mask,3))), ...
    sum(reshape(roi_mask,[],size(roi_mask,3)),1)');

smooth_size = 1;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

fluor_roi = permute(reshape(transpose(U_roi(:,1:n_vs)* ...
    reshape(permute(fluor_allcat(:,:,:,1),[3,2,1,4]),n_vs,[])), ...
    size(fluor_allcat,2),size(fluor_allcat,1),n_rois),[2,1,3]);

fluor_roi_diff = diff(padarray(convn(fluor_roi,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both'),[],2);

% % hemidiff
% fluor_roi_diff(fluor_roi_diff < 0) = 0;
% fluor_roi_diff = fluor_roi_diff(:,:,1:size(wf_roi,1)) - ...
%     fluor_roi_diff(:,:,size(wf_roi,1)+1:end);

% (zero negatives, normalize)
fluor_roi_diff(fluor_roi_diff < 0) = 0;
fluor_roi_diff = bsxfun(@rdivide,fluor_roi_diff,nanstd(reshape(fluor_roi_diff,[],1,size(wf_roi,1)),[],1));

t_diff =  conv(t,[1,1]/2,'valid');

figure; hold on;
for curr_plot_condition = 1:size(plot_conditions,1)
    
    curr_trials = plot_id == curr_plot_condition & use_rxn;    
    curr_data = fluor_roi_diff(curr_trials,:,:);
    
    if move_align
        % Re-align to movement onset
        t_leeway = 0.5;
        leeway_samples = round(t_leeway*sample_rate);
        curr_move_idx = move_idx(curr_trials);
        for i = 1:size(curr_data,1)
            curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
        end
    end
    
    curr_data_mean = squeeze(nanmean(curr_data,1));

    curr_contrast_side = max(trial_contrast_allcat(curr_trials))*sign(max(trial_side_allcat(curr_trials)));
    contrast_side_idx = find(curr_contrast_side == contrast_side_val);
    curr_col = contrast_side_col(contrast_side_idx,:);   
    
    curr_linewidth = mean(trial_choice_allcat(curr_trials))+2;
    
    if curr_contrast_side == 0
        switch max(trial_choice_allcat(curr_trials))
            case -1
                curr_col = 'm';
            case 1
                curr_col = 'c';
        end
    end
    
    AP_stackplot(curr_data_mean,t_diff,3.5,false,curr_col,{wf_roi.area});
    
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');
title('Widefield ROI');

% Get fluorescence in kernel ROIs
kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
load(kernel_roi_fn);
n_rois = size(kernel_roi,3);
roi_mask = kernel_roi;

fluor_roi = nan(size(fluor_allcat,1),size(fluor_allcat,2)-1,n_rois);
U_roi = bsxfun(@rdivide,transpose(reshape(U_master,[],size(U_master,3))'* ...
    reshape(roi_mask,[],size(roi_mask,3))), ...
    sum(reshape(roi_mask,[],size(roi_mask,3)),1)');

smooth_size = 1;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

fluor_roi = permute(reshape(transpose(U_roi(:,1:n_vs)* ...
    reshape(permute(fluor_allcat(:,:,:,1),[3,2,1,4]),n_vs,[])), ...
    size(fluor_allcat,2),size(fluor_allcat,1),n_rois),[2,1,3]);

fluor_roi_diff = diff(padarray(convn(fluor_roi,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both'),[],2);

% (zero negatives, normalize)
fluor_roi_diff(fluor_roi_diff < 0) = 0;
fluor_roi_diff = bsxfun(@rdivide,fluor_roi_diff,nanstd(reshape(fluor_roi_diff,[],1,n_rois),[],1));

t_diff =  conv(t,[1,1]/2,'valid');

figure; hold on;
for curr_plot_condition = 1:size(plot_conditions,1)
    
    curr_trials = plot_id == curr_plot_condition & use_rxn;
    curr_data = fluor_roi_diff(curr_trials,:,:);
    
    if move_align
        % Re-align to movement onset
        t_leeway = 0.5;
        leeway_samples = round(t_leeway*sample_rate);
        curr_move_idx = move_idx(curr_trials);
        for i = 1:size(curr_data,1)
            curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
        end
    end
    
    curr_data_mean = squeeze(nanmean(curr_data,1));
    
    curr_contrast_side = max(trial_contrast_allcat(curr_trials))*sign(max(trial_side_allcat(curr_trials)));
    contrast_side_idx = find(curr_contrast_side == contrast_side_val);
    curr_col = contrast_side_col(contrast_side_idx,:);   
        
    if curr_contrast_side == 0
        switch max(trial_choice_allcat(curr_trials))
            case -1
                curr_col = 'm';
            case 1
                curr_col = 'c';
        end
    end
    
    AP_stackplot(curr_data_mean,t_diff,5,false,curr_col,1:n_rois);
    
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');
title('Widefield striatum kernel ROI');

% Fluorescence with kernels
ctx_str_kernel_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\ctx_str_kernel';
load(ctx_str_kernel_fn)
% (interpolate to upsampled)
kernel_t_upsample = ctx_str_kernel.t(1):1/sample_rate:ctx_str_kernel.t(end);
kernel_upsample = permute(interp1(ctx_str_kernel.t,permute(ctx_str_kernel.V,[2,1,3]),kernel_t_upsample),[2,1,3]);

fluor_str_kernel_act = zeros(size(fluor_allcat,1),size(fluor_allcat,2)-1,n_depths,'single');
for curr_depth = 1:n_depths
    curr_conv = zeros(size(fluor_allcat,1),size(fluor_allcat,2)-1,size(fluor_allcat,3),'single');
    for curr_v = 1:n_vs
        curr_conv(:,:,curr_v) = ...
            conv2(diff(fluor_allcat(:,:,curr_v),[],2),permute(kernel_upsample(curr_v,:,curr_depth),[3,2,1]),'same');
    end
    fluor_str_kernel_act(:,:,curr_depth) = sum(curr_conv,3);
end

fluor_str_kernel_act = bsxfun(@rdivide,fluor_str_kernel_act, ...
    nanstd(reshape(fluor_str_kernel_act,[],1,n_depths),[],1));

figure; hold on;
for curr_plot_condition = 1:size(plot_conditions,1)
    
    curr_trials = plot_id == curr_plot_condition & use_rxn;
    curr_data = fluor_str_kernel_act(curr_trials,:,:);
    
    if move_align
        % Re-align to movement onset
        t_leeway = 0.5;
        leeway_samples = round(t_leeway*sample_rate);
        curr_move_idx = move_idx(curr_trials);
        for i = 1:size(curr_data,1)
            curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
        end
    end
    
    curr_data_mean = squeeze(nanmean(curr_data,1));
    
    curr_contrast_side = max(trial_contrast_allcat(curr_trials))*sign(max(trial_side_allcat(curr_trials)));
    contrast_side_idx = find(curr_contrast_side == contrast_side_val);
    curr_col = contrast_side_col(contrast_side_idx,:);
        
    if curr_contrast_side == 0
        switch max(trial_choice_allcat(curr_trials))
            case -1
                curr_col = 'm';
            case 1
                curr_col = 'c';
        end
    end
    
    AP_stackplot(curr_data_mean,t_diff,4,false,curr_col,1:n_depths);
    
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');
title('Widefield striatum kernels');


% Plot MUA
figure; hold on;
for curr_plot_condition = 1:size(plot_conditions,1)
    
    curr_trials = plot_id == curr_plot_condition & use_rxn;
    curr_data = mua_allcat(curr_trials,:,:);
    
    if move_align
        % Re-align to movement onset
        t_leeway = 0.5;
        leeway_samples = round(t_leeway*sample_rate);
        curr_move_idx = move_idx(curr_trials);
        for i = 1:size(curr_data,1)
            curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
        end
    end
    
    curr_data_mean = squeeze(nanmean(curr_data,1));
    
    curr_contrast_side = max(trial_contrast_allcat(curr_trials))*sign(max(trial_side_allcat(curr_trials)));
    contrast_side_idx = find(curr_contrast_side == contrast_side_val);
    curr_col = contrast_side_col(contrast_side_idx,:);
        
    if curr_contrast_side == 0
        switch max(trial_choice_allcat(curr_trials))
            case -1
                curr_col = 'm';
            case 1
                curr_col = 'c';
        end
    end
    
    AP_stackplot(curr_data_mean,t,2,false,curr_col,1:n_depths);
    
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');
title('Striatum');


%% Plot activity 3 areas together

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_Udf_kernel-str_4_depths_long'];
% data_fn = 'all_trial_activity_df_msn_earlymove.mat';

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);
reward_all = cellfun(@(data,use_expts) data(use_expts),reward_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);
reward_all = reward_all(use_animals);

% Get widefield components
n_vs = 200;

% MUA depths
n_depths = size(mua_all{1}{1},3);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_vs,1);
mua_allcat = nan(0,length(t),n_depths,1);
wheel_allcat = nan(0,length(t),1);
reward_allcat = nan(0,length(t),1);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    n_align = size(fluor_cat,4);
    
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x,[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
%     wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
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
 
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat];

    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    reward_allcat = [reward_allcat;vertcat(reward_all{curr_animal}{:})];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];
   
end

% Load the master U, get fluorescence in ROIs
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);
n_rois = numel(wf_roi);

roi_mask = cat(3,wf_roi.mask);

U_roi = bsxfun(@rdivide,transpose(reshape(U_master,[],size(U_master,3))'* ...
    reshape(roi_mask,[],size(roi_mask,3))), ...
    sum(reshape(roi_mask,[],size(roi_mask,3)),1)');

smooth_size = 1;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

fluor_roi = permute(reshape(transpose(U_roi(:,1:n_vs)* ...
    reshape(permute(fluor_allcat(:,:,:,1),[3,2,1,4]),n_vs,[])), ...
    size(fluor_allcat,2),size(fluor_allcat,1),n_rois),[2,1,3]);

fluor_roi_diff = [zeros(size(fluor_roi,1),1,size(fluor_roi,3)), ...
    diff(padarray(convn(fluor_roi,smWin,'valid'), ...
    [0,floor(size(smWin,2)/2)],'replicate','both'),[],2)];

% (zero negatives, normalize)
fluor_roi_diff(fluor_roi_diff < 0) = 0;
fluor_roi_diff = bsxfun(@rdivide,fluor_roi_diff,nanstd(reshape(fluor_roi_diff,[],1,size(wf_roi,1)),[],1));

% Concatenate unilateral 
% fluor_unilateral_allcat = cat(1,fluor_unilateral_allcat(:,:,1:length(wf_roi),:), ...
%     fluor_unilateral_allcat(:,:,length(wf_roi)+1:end,:));
% % (if concatenating side - messy, just here for the moment)
% trial_contrast_allcat = [trial_contrast_allcat;trial_contrast_allcat];
% trial_side_allcat = [trial_side_allcat;-trial_side_allcat];
% trial_choice_allcat = [trial_choice_allcat;-trial_choice_allcat];
% move_t = repmat(move_t,1,2);

% Low-pass filter fluorescence (where does this ~10 Hz crap come from?)
lowpassCutoff = 6; % Hz
[b100s, a100s] = butter(2, lowpassCutoff/((sample_rate)/2), 'low');
fluor_roi_diff_filt = filter(b100s,a100s,fluor_roi_diff,[],2);
mua_allcat_filt = filter(b100s,a100s,mua_allcat,[],2);

activity_cat = cat(3,fluor_roi_diff_filt,mua_allcat_filt);
area_labels = [{wf_roi(:,1).area}, ...
    cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false)];

% Sort trials by stim-to-move time
% (by position)
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
% (by max speed)
% [~,move_idx] = max(abs(diff(wheel_allcat(:,:,1),[],2)),[],2);

move_t = t(move_idx)';
[~,sort_idx] = sort(move_idx);

%%% 3 areas together
plot_areas = [9,10,11];

% Get conditions for each trial, plot selected
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];

contrast_side_col = colormap_BlueWhiteRed(5);
contrast_side_col(6,:) = 0;
contrast_side_val = unique(sort([-contrasts,contrasts]))';

% % contrast, side, choice
% plot_conditions = ...
%     [contrasts,contrasts; ...
%     -ones(1,7),ones(1,5); ...
%     ones(1,6),-ones(1,6)]';

plot_conditions = ...
    [1,0,1,0; ...
    1,-1,-1,-1; ...
    -1,-1,1,1]';

% plot_conditions = ...
%     [0.06,0.06,0.06; ...
%     1,1,-1; ...
%     -1,1,-1]';

use_rxn = move_t > 0 & move_t < 0.5;
    
[~,plot_id] = ismember( ...
    [trial_contrast_allcat > 0,trial_side_allcat,trial_choice_allcat], ...
    plot_conditions,'rows');

move_align = true;

figure; hold on;
for curr_plot_condition = 1:size(plot_conditions,1)
    
    curr_trials = plot_id == curr_plot_condition & use_rxn;    
    curr_data = activity_cat(curr_trials,:,plot_areas);
    
    if move_align
        % Re-align to movement onset
        t_leeway = 0.5;
        leeway_samples = round(t_leeway*sample_rate);
        curr_move_idx = move_idx(curr_trials);
        for i = 1:size(curr_data,1)
            curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
        end
    end
    
    curr_data_mean = squeeze(nanmean(curr_data,1));

    curr_contrast_side = max(trial_contrast_allcat(curr_trials))*sign(max(trial_side_allcat(curr_trials)));
    contrast_side_idx = find(curr_contrast_side == contrast_side_val);
    curr_col = contrast_side_col(contrast_side_idx,:);   
    
    curr_linewidth = mean(trial_choice_allcat(curr_trials))+2;
    
    if curr_contrast_side == 0
        switch max(trial_choice_allcat(curr_trials))
            case -1
                curr_col = 'm';
            case 1
                curr_col = 'c';
        end
    end
    
    plot3(curr_data_mean(:,1),curr_data_mean(:,2),curr_data_mean(:,3), ...
        'color',curr_col,'linewidth',2);
    plot_t = find(t > 0,1);
    plot3(curr_data_mean(plot_t,1),curr_data_mean(plot_t,2),curr_data_mean(plot_t,3), ...
        '.','color',curr_col,'MarkerSize',30);
 
end

view([45,45,45]);axis normal vis3d;

xlabel(area_labels{plot_areas(1)});
ylabel(area_labels{plot_areas(2)});
zlabel(area_labels{plot_areas(3)});



%% Time-averaged day-concatenated activity (move L/R) (animals cat'd)

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_Udf_kernel-str_4_depths_long'];

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);
reward_all = cellfun(@(data,use_expts) data(use_expts),reward_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);
reward_all = reward_all(use_animals);

% Get widefield components
n_vs = 200;

% MUA depths
n_depths = size(mua_all{1}{1},3);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_vs,1);
mua_allcat = nan(0,length(t),n_depths,1);
wheel_allcat = nan(0,length(t),1);
reward_allcat = nan(0,length(t),1);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    n_align = size(fluor_cat,4);
    
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x,[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
%     wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
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
 
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    reward_allcat = [reward_allcat;vertcat(reward_all{curr_animal}{:})];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];
   
end

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Get conditions for each trial, plot selected
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
sidecontrasts = unique(bsxfun(@times,contrasts',sides));

use_rxn = move_t > 0 & move_t < 0.5;

move_align = true;

% Get fluorescence in widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);
n_rois = numel(wf_roi);

roi_mask = cat(3,wf_roi.mask);

U_roi = bsxfun(@rdivide,transpose(reshape(U_master,[],size(U_master,3))'* ...
    reshape(roi_mask,[],size(roi_mask,3))), ...
    sum(reshape(roi_mask,[],size(roi_mask,3)),1)');

smooth_size = 1;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

fluor_roi = permute(reshape(transpose(U_roi(:,1:n_vs)* ...
    reshape(permute(fluor_allcat(:,:,:,1),[3,2,1,4]),n_vs,[])), ...
    size(fluor_allcat,2),size(fluor_allcat,1),n_rois),[2,1,3]);

% Low-pass filter fluorescence (where does this ~10 Hz crap come from?)
lowpassCutoff = 6; % Hz
[b100s, a100s] = butter(2, lowpassCutoff/((sample_rate)/2), 'low');
fluor_roi_filt = filter(b100s,a100s,fluor_roi,[],2);
mua_allcat_filt = filter(b100s,a100s,mua_allcat,[],2);

% Fluorescence derivative
fluor_roi_filt_diff = [zeros(size(fluor_roi_filt,1),1,size(fluor_roi_filt,3)), ...
    diff(padarray(convn(fluor_roi_filt,smWin,'valid'), ...
    [0,floor(size(smWin,2)/2)],'replicate','both'),[],2)];

% fluor_roi_filt_diff = fluor_roi_filt_diff(:,:,1:size(wf_roi,1)) - ...
%     fluor_roi_filt_diff(:,:,size(wf_roi,1)+1:end);

% (zero negatives, normalize)
fluor_roi_filt_diff(fluor_roi_filt_diff < 0) = 0;
fluor_roi_filt_diff = bsxfun(@rdivide,fluor_roi_filt_diff, ...
    nanstd(reshape(fluor_roi_filt_diff,[],1,size(wf_roi,1)),[],1));

activity_cat = cat(3,fluor_roi_filt_diff,mua_allcat_filt);
area_labels = [{wf_roi(:,1).area}, ...
    cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false)];

% re-align to movement onset
t_leeway = 0.5;
activity_cat_move = activity_cat;
leeway_samples = round(t_leeway*(sample_rate));
for i = 1:size(activity_cat_move,1)
    activity_cat_move(i,:,:) = circshift(activity_cat_move(i,:,:),-move_idx(i)+leeway_samples,2);
end

plot_area = 5;
plot_t = t > -0.05 & t < 0.05;
use_data = activity_cat_move(:,:,plot_area);

trial_act_timeavg = permute(nanmean(use_data(:,plot_t,:),2),[1,3,4,2]);

% Split activity by contrast*side/choice
trial_act_timeavg_split = cell(length(sidecontrasts),2);
trial_act_timeavg_split(:,1) = arrayfun(@(x) ...
    trial_act_timeavg(trial_side_allcat.*trial_contrast_allcat == sidecontrasts(x) & ...
    trial_choice_allcat == -1 & use_rxn),1:length(sidecontrasts),'uni',false);
trial_act_timeavg_split(:,2) = arrayfun(@(x) ...
    trial_act_timeavg(trial_side_allcat.*trial_contrast_allcat == sidecontrasts(x) & ...
    trial_choice_allcat == 1 & use_rxn),1:length(sidecontrasts),'uni',false);
trial_act_timeavg_split_median = cellfun(@nanmedian,trial_act_timeavg_split);

figure; 
% Plot trials from all animals concatenated
subplot(1,3,1); hold on;
p1 = distributionPlot(trial_act_timeavg_split(:,1),'distWidth',0.5, ...
    'xValues',(1:11)-0.25,'histOri','left','showMM',0,'color',[0.6,0,0.6]);
p2 = distributionPlot(trial_act_timeavg_split(:,2),'distWidth',0.5, ...
    'xValues',(1:11)+0.25,'histOri','right','showMM',0,'color',[0,0.6,0]);

trial_act_timeavg_split_all_median = cellfun(@nanmedian,trial_act_timeavg_split);
trial_act_timeavg_split_all_mad = cellfun(@(x) mad(x,1),trial_act_timeavg_split);
plot((1:11),trial_act_timeavg_split_all_median(:,1),'color',[0.9,0.6,0.9],'linewidth',2)
plot((1:11),trial_act_timeavg_split_all_median(:,2),'color',[0.6,0.9,0.6],'linewidth',2)

set(gca,'XTick',1:11,'XTickLabel',cellfun(@num2str,num2cell(sidecontrasts),'uni',false));
xlabel('Contrast*Side');
ylabel(area_labels{plot_area});
legend([p1{1}(6),p2{1}(6)],{'Move left','Move right'});

% Plot trials from both choices conconatenated
subplot(1,3,2); hold on;
trial_act_timeavg_combined = arrayfun(@(x) ...
    vertcat(trial_act_timeavg_split{x,:}), ...
    1:size(trial_act_timeavg_split),'uni',false)';
p1 = distributionPlot(trial_act_timeavg_combined,'distWidth',0.5, ...
    'xValues',1:11,'histOri','center','showMM',0,'color',[0.7,0.7,0.7]);

trial_act_timeavg_combined_all_median = cellfun(@nanmedian,trial_act_timeavg_combined);
plot(1:11,trial_act_timeavg_combined_all_median,'color','k','linewidth',2)

set(gca,'XTick',1:11,'XTickLabel',cellfun(@num2str,num2cell(sidecontrasts),'uni',false));
xlabel('Contrast*Side');
ylabel(area_labels{plot_area});

% Plot just medians by correctly scaled x-axis
subplot(1,3,3); hold on;

AP_errorfill(sidecontrasts,trial_act_timeavg_split_all_median(:,1), ...
    trial_act_timeavg_split_all_mad(:,1),[0.6,0,0.6],0.5);
AP_errorfill(sidecontrasts,trial_act_timeavg_split_all_median(:,2), ...
    trial_act_timeavg_split_all_mad(:,2),[0,0.6,0],0.5);

p1 = plot(sidecontrasts,trial_act_timeavg_split_all_median(:,1),'color',[0.6,0,0.6],'linewidth',2);
p2 = plot(sidecontrasts,trial_act_timeavg_split_all_median(:,2),'color',[0,0.6,0],'linewidth',2);
p3 = plot(sidecontrasts,trial_act_timeavg_combined_all_median,'color','k','linewidth',2);

xlabel('Contrast*Side');
ylabel(area_labels{plot_area});
legend([p1,p2,p3],{'Move left','Move right','All'})


%% Time-averaged day-concatenated activity (move L/R) (animals mean'd)

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_df_kernel-str_earlymove_' num2str(n_aligned_depths) '_depths.mat'];
% data_fn = 'all_trial_activity_df_msn_earlymove.mat';
load([data_path filesep data_fn]);

n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);


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
plot_act = 'fluor_cat';
plot_area = 2;
plot_align = 2;
% plot_t = t >= 0.08 & t <= 0.12; % stim response
plot_t = t >= -0.14 & t <= -0.08; % early premove
% plot_t = t >= -0.08 & t <= 0; % late premove
% plot_t = t >= 0.55 & t <= 0.6; % post-beep

figure;
trial_act_timeavg_split_all = cell(length(sidecontrasts),2,length(D_all));
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
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
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
    
    % Get time-averaged activity for each trial
    switch plot_act
        case 'fluor'
            use_data = fluor_cat_hemidiff_norm;
            area_labels = {wf_roi.area};
            area_labels(n_rois/2+1:end) = cellfun(@(x) ...
                ['\Delta' x(1:end-2)],area_labels(n_rois/2+1:end),'uni',false);
            
        case 'fluor_cat'            
            if plot_area > n_rois/2
                error('Concatenating fluor: only use L ROIs')
            end
            % concatenate L/R and reorder trial ID as ipsi/contra
            use_data = [fluor_cat_norm(:,:,1:n_rois/2,:); ...
                fluor_cat_norm(:,:,n_rois/2+1:end,:)];
            area_labels = cellfun(@(x) [x(1:end-2) '_{unilateral}'], ...
                {wf_roi(1:n_rois/2).area},'uni',false);
            
            trial_contrast = [trial_contrast;trial_contrast];
            trial_side = [trial_side;-trial_side];
            trial_choice = [trial_choice;-trial_choice];
            
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
    
    trial_act_timeavg_split_all(:,:,curr_animal) = trial_act_timeavg_split;
    
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

% Plot just medians by correctly scaled x-axis
trial_act_timeavg_split_all_median = nanmean(cellfun(@nanmedian,trial_act_timeavg_split_all),3);
trial_act_timeavg_split_all_std = nanstd(cellfun(@nanmedian,trial_act_timeavg_split_all),[],3)./ ...
    sqrt(sum(~isnan(cellfun(@nanmedian,trial_act_timeavg_split_all)),3));

figure; hold on;

AP_errorfill(sidecontrasts,trial_act_timeavg_split_all_median(:,1), ...
    trial_act_timeavg_split_all_std(:,1),[0.6,0,0.6],0.5);
AP_errorfill(sidecontrasts,trial_act_timeavg_split_all_median(:,2), ...
    trial_act_timeavg_split_all_std(:,2),[0,0.6,0],0.5);

p1 = plot(sidecontrasts,trial_act_timeavg_split_all_median(:,1),'color',[0.6,0,0.6],'linewidth',2)
p2 = plot(sidecontrasts,trial_act_timeavg_split_all_median(:,2),'color',[0,0.6,0],'linewidth',2)

xlabel('Contrast*Side');
ylabel(area_labels{plot_area});
legend([p1,p2],{'Move left','Move right','All'})


%% Plot move L-R activity difference of mean concatenated activity

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_df_kernel-str_earlymove_' num2str(n_aligned_depths) '_depths.mat'];
% data_fn = 'all_trial_activity_df_msn_earlymove.mat';
load([data_path filesep data_fn]);

n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);


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
trial_act_split_fluor_all = cell(length(sidecontrasts),2,n_rois,2);
trial_act_split_mua_all = cell(length(sidecontrasts),2,n_depths,2);
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
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
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
    
    % Get time-averaged activity for each trial
    % (fluor)
    for curr_align = 1:2
        for curr_roi = 1:n_rois
            % Split activity by contrast*side/choice
            trial_act_split_fluor = cell(length(sidecontrasts),2);
            trial_act_split_fluor(:,1) = arrayfun(@(x) ...
                fluor_cat_hemidiff_norm(trial_side.*trial_contrast == sidecontrasts(x) & ...
                trial_choice == -1,:,curr_roi,curr_align), ...
                1:length(sidecontrasts),'uni',false);
            trial_act_split_fluor(:,2) = arrayfun(@(x) ...
                fluor_cat_hemidiff_norm(trial_side.*trial_contrast == sidecontrasts(x) & ...
                trial_choice == 1,:,curr_roi,curr_align), ...
                1:length(sidecontrasts),'uni',false);
            
            trial_act_split_fluor_all(:,:,curr_roi,curr_align) = ...
                cellfun(@(curr,all) [curr;all], ...
                trial_act_split_fluor, ...
                trial_act_split_fluor_all(:,:,curr_roi,curr_align),'uni',false);
        end
    end
    
    % (mua)
    for curr_align = 1:2
        for curr_depth = 1:n_depths
            % Split activity by contrast*side/choice
            trial_act_split_mua = cell(length(sidecontrasts),2);
            trial_act_split_mua(:,1) = arrayfun(@(x) ...
                mua_cat_norm(trial_side.*trial_contrast == sidecontrasts(x) & ...
                trial_choice == -1,:,curr_depth,curr_align), ...
                1:length(sidecontrasts),'uni',false);
            trial_act_split_mua(:,2) = arrayfun(@(x) ...
                mua_cat_norm(trial_side.*trial_contrast == sidecontrasts(x) & ...
                trial_choice == 1,:,curr_depth,curr_align), ...
                1:length(sidecontrasts),'uni',false);
            
            trial_act_split_mua_all(:,:,curr_depth,curr_align) = ...
                cellfun(@(curr,all) [curr;all], ...
                trial_act_split_mua, ...
                trial_act_split_mua_all(:,:,curr_depth,curr_align),'uni',false);
        end
    end
        
          
end

% Get mean of all trials together
[d1,d2,d3,d4] = size(trial_act_split_fluor_all);
fluor_reshape = reshape(trial_act_split_fluor_all,d1,[]);
trial_act_cat_fluor_all_mean = reshape(arrayfun(@(x) nanmean(vertcat(fluor_reshape{:,x}),1), ...
    1:size(fluor_reshape,2),'uni',false),d2,d3,d4);
trial_act_cat_fluor_all_diff = squeeze(cellfun(@(x,y) y-x, ...
    trial_act_cat_fluor_all_mean(1,:,:), ...
    trial_act_cat_fluor_all_mean(2,:,:),'uni',false));

[d1,d2,d3,d4] = size(trial_act_split_mua_all);
mua_reshape = reshape(trial_act_split_mua_all,d1,[]);
trial_act_cat_mua_all_mean = reshape(arrayfun(@(x) nanmean(vertcat(mua_reshape{:,x}),1), ...
    1:size(mua_reshape,2),'uni',false),d2,d3,d4);
trial_act_cat_mua_all_diff = squeeze(cellfun(@(x,y) y-x, ...
    trial_act_cat_mua_all_mean(1,:,:), ...
    trial_act_cat_mua_all_mean(2,:,:),'uni',false));

% Get means of trials by condition
trial_act_split_fluor_all_mean = cellfun(@(x) ...
    nanmean(x,1),trial_act_split_fluor_all,'uni',false);
trial_act_split_fluor_all_mean_diff = squeeze(cellfun(@(move_l,move_r) ...
    move_l-move_r,trial_act_split_fluor_all_mean(:,1,:,:), ...
    trial_act_split_fluor_all_mean(:,2,:,:),'uni',false));

trial_act_split_mua_all_medn = cellfun(@(x) ...
    nanmean(x,1),trial_act_split_mua_all,'uni',false);
trial_act_split_mua_all_mean_diff = squeeze(cellfun(@(move_l,move_r) ...
    move_l-move_r,trial_act_split_mua_all_medn(:,1,:,:), ...
    trial_act_split_mua_all_medn(:,2,:,:),'uni',false));

% Plot by grand mean
figure('Name','Move R - Move L');
subplot(1,4,1);
AP_stackplot(vertcat(trial_act_cat_fluor_all_diff{:,1})',t,3,false,'k',{wf_roi.area});
line([0,0],ylim);
xlabel('Time from stim');

subplot(1,4,2);
AP_stackplot(vertcat(trial_act_cat_fluor_all_diff{:,2})',t,3,false,'k',{wf_roi.area});
line([0,0],ylim);
xlabel('Time from move');

subplot(1,4,3);
AP_stackplot(vertcat(trial_act_cat_mua_all_diff{:,1})',t,3,false,'r',1:4);
line([0,0],ylim);
xlabel('Time from stim');

subplot(1,4,4);
AP_stackplot(vertcat(trial_act_cat_mua_all_diff{:,2})',t,3,false,'r',1:4);
line([0,0],ylim);
xlabel('Time from move');

% Plot by condition
figure('Name','Move R - Move L: Fluor');
for curr_align = 1:2
    for curr_roi = 1:n_rois/2
        subplot(2,n_rois/2,n_rois/2*(curr_align-1)+curr_roi);
        imagesc(t,sidecontrasts,vertcat(trial_act_split_fluor_all_mean_diff{:,curr_roi,curr_align}))
        colormap(colormap_BlueWhiteRed);
        caxis([-max(abs(caxis)),max(abs(caxis))]);
        xlabel('Time (s)');
        ylabel('Contrast*Side');
        title([wf_roi(curr_roi).area ' Align ' num2str(curr_align)])
    end
end

figure('Name','Move R - Move L: dFluor');
for curr_align = 1:2
    for curr_roi = n_rois/2+1:n_rois
        subplot(2,n_rois/2,n_rois/2*(curr_align-1)+curr_roi-n_rois/2);
        imagesc(t,sidecontrasts,vertcat(trial_act_split_fluor_all_mean_diff{:,curr_roi,curr_align}))
        colormap(colormap_BlueWhiteRed);
        caxis([-max(abs(caxis)),max(abs(caxis))]);
        xlabel('Time (s)');
        ylabel('Contrast*Side');
        title([wf_roi(curr_roi).area ' Align ' num2str(curr_align)])
    end
end

figure('Name','Move R - Move L: MUA');
for curr_align = 1:2
    for curr_depth = 1:n_depths
        subplot(2,n_depths,n_depths*(curr_align-1)+curr_depth);
        imagesc(t,sidecontrasts,vertcat(trial_act_split_mua_all_mean_diff{:,curr_depth,curr_align}))
        colormap(colormap_BlueWhiteRed);
        caxis([-max(abs(caxis)),max(abs(caxis))]);
        xlabel('Time (s)');
        ylabel('Contrast*Side');
        title(['Str ' num2str(curr_depth) ' Align ' num2str(curr_align)]);
    end
end

%% Linear regression on day-concatenated activity (weights & expl var)

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_df_kernel-str_earlymove_' num2str(n_aligned_depths) '_depths.mat'];
% data_fn = 'all_trial_activity_df_earlymove.mat';
load([data_path filesep data_fn]);

% Get time
framerate = 35;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = framerate*upsample_factor;
t = raster_window(1):(1/sample_rate):raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);


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
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
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
    kernel_t = [-0.1,0.1];
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
   
    AP_print_progress_fraction(curr_animal,length(D_all))
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


%% Regress day-concatenated mua from cortex in time (timepoint specific)

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_df_kernel-str_earlymove_' num2str(n_aligned_depths) '_depths.mat'];
load([data_path filesep data_fn]);

% Get time
framerate = 35;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = framerate*upsample_factor;
t = raster_window(1):(1/sample_rate):raster_window(2);

% Set up conditions
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);


% Get widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

% Get number of MUA depths
n_depths = size(mua_all{1}{1},3);

interarea_kernel_all = cell(size(D_all));
interarea_expl_var_all = cell(size(D_all));
predicted_signal_all = nan(n_conditions,length(t),length(D_all));
for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate behavioural data and select trials to use
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
    D.response = cell2mat(cellfun(@(x) x.response,D_all{curr_animal},'uni',false));
    D.repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all{curr_animal},'uni',false));
    
    D.day = trial_day;
    
    % Get average predicted activity by condition
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
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
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
    
    % For each time point, regress fluorescence of one area from others   
    regressor_areas = 1:16;
    regress_area = 3;
    t_surround = -10:10;
    interarea_kernel = nan(length(regressor_areas),length(t),length(t_surround));
    interarea_expl_var = nan(1,length(t));
    predicted_signal = nan(size(D.stimulus,1),length(t));
    for curr_t = -min(t_surround)+1:length(t)-max(t_surround)
        
        curr_regressors = reshape(fluor_cat_norm( ...
            :,curr_t+t_surround,regressor_areas,2),[],length(regressor_areas)*length(t_surround));
        
        curr_signal = mua_cat_norm(:,curr_t,regress_area,2);   
        
        % Only use non-nan regressors and signals
        use_trials = ~all(isnan(curr_regressors),2) & all(~isnan(curr_signal),2);       
        lambda = 5e1;
        cv = 5;
        
        [curr_kernel,curr_predicted_signal,curr_expl_var] = AP_regresskernel( ...
            curr_regressors(use_trials,:)',curr_signal(use_trials,:)', ...
            [],lambda,[false,false],cv);
        
        interarea_kernel(:,curr_t,:) = permute(reshape(curr_kernel,[],length(regressor_areas)),[2,3,1]);
        interarea_expl_var(curr_t) = curr_expl_var.total;
        predicted_signal(use_trials,curr_t) = curr_predicted_signal;
        
    end 
    
    % Get the maximum amplitude kernel with sign
    interarea_kernel_max = max(interarea_kernel,[],3);
    interarea_kernel_min = min(interarea_kernel,[],3);
    interarea_kernel_pos = abs(interarea_kernel_max) > abs(interarea_kernel_min);
    interarea_kernel_magnitude = interarea_kernel_max.*interarea_kernel_pos + ...
        interarea_kernel_min.*~interarea_kernel_pos;
    
    interarea_kernel_mean = mean(interarea_kernel,3);
    
    interarea_kernel_all{curr_animal} = interarea_kernel_mean;
    interarea_expl_var_all{curr_animal} = interarea_expl_var;
  
    predicted_signal_all(unique(trial_id),:,curr_animal) = grpstats(predicted_signal,trial_id);
   
    AP_print_progress_fraction(curr_animal,length(D_all));
    
end

interarea_kernel_mean = nanmean(cat(3,interarea_kernel_all{:}),3);
interarea_expl_var_mean = nanmean(cat(3,interarea_expl_var_all{:}),3);
predicted_signal_all_mean = nanmean(predicted_signal_all,3);

figure;
subplot(3,1,1); hold on;
set(gca,'ColorOrder',[autumn(8);winter(8)]);
plot(t,interarea_kernel_mean')
axis tight
line([0,0],ylim,'color','k');
ylabel('Weight');

subplot(3,1,2);
plot(t,interarea_expl_var_mean,'k','linewidth',2);
axis tight
line([0,0],ylim,'color','k');
ylabel('Explained variance');

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
% plot_move = -1;
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
% plot_conditions = find(conditions(:,1) == 0.06 & ...
%     conditions(:,2) == 1 & ...
%     conditions(:,4) == plot_timing);
% plot_color = zeros(length(plot_conditions),3);
% plot_color(conditions(plot_conditions,3) == -1,1) = 1;
% plot_color(conditions(plot_conditions,3) == 1,3) = 1;

% % (to plot zero-contrast by movement)
% plot_conditions = find(conditions(:,1) == 0 & ...
%     conditions(:,4) == plot_timing);
% plot_color = zeros(length(plot_conditions),3);
% plot_color(conditions(plot_conditions,3) == -1,1) = 1;
% plot_color(conditions(plot_conditions,3) == 1,3) = 1;

% Plot predicted signal
subplot(3,1,3); 
hold on;
for curr_condition_idx = 1:length(plot_conditions);
    curr_condition = plot_conditions(curr_condition_idx);
    plot(t,squeeze(predicted_signal_all_mean(curr_condition,:)), ...
        'color',plot_color(curr_condition_idx,:),'linewidth',2);
end
axis tight;
line([0,0],ylim,'color','k');
xlabel('Time from stim');

% Plot weights by spatial ROI
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

roi_cat = cat(3,wf_roi.mask);

use_t_weights = t > -0.1 & t < 0;

interarea_kernel_mean = nanmean(cat(3,interarea_kernel_all{:}),3);
roi_weights = nanmean(interarea_kernel_mean(:,use_t_weights),2);

figure; hold on
set(gca,'YDir','reverse');
AP_reference_outline('ccf_aligned','k');
for curr_roi = 1:n_rois
    curr_roi_boundary = cell2mat(bwboundaries(roi_cat(:,:,curr_roi)));
    patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),roi_weights(curr_roi));
end
axis image off;
colormap(colormap_BlueWhiteRed);
caxis([-max(abs(caxis)),max(abs(caxis))]);

%% Regress day-concatenated cortex from cortex in time (timepoint specific)

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_df_kernel-str_earlymove_' num2str(n_aligned_depths) '_depths.mat'];
load([data_path filesep data_fn]);

% Get time
framerate = 35;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = framerate*upsample_factor;
t = raster_window(1):(1/sample_rate):raster_window(2);

% Set up conditions
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);


% Get widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

% Get number of MUA depths
n_depths = size(mua_all{1}{1},3);

interarea_kernel_all = cell(size(D_all));
interarea_expl_var_all = cell(size(D_all));
predicted_signal_all = nan(n_conditions,length(t),length(D_all));
for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate behavioural data and select trials to use
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
    D.response = cell2mat(cellfun(@(x) x.response,D_all{curr_animal},'uni',false));
    D.repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all{curr_animal},'uni',false));
    
    D.day = trial_day;
    
    % Get average predicted activity by condition
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
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
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
    
    % For each time point, regress fluorescence of one area from others   
    regress_area = 14;
    regressor_areas = setdiff(1:16,regress_area);    
    t_surround = -10:10;
    interarea_kernel = nan(length(regressor_areas),length(t),length(t_surround));
    interarea_expl_var = nan(1,length(t));
    predicted_signal = nan(size(D.stimulus,1),length(t));
    for curr_t = -min(t_surround)+1:length(t)-max(t_surround)
        
        curr_regressors = reshape(fluor_cat_hemidiff_norm( ...
            :,curr_t+t_surround,regressor_areas,2),[],length(regressor_areas)*length(t_surround));
        
        curr_signal = fluor_cat_hemidiff_norm(:,curr_t,regress_area,2);   
        
        % Only use non-nan regressors and signals
        use_trials = ~all(isnan(curr_regressors),2) & all(~isnan(curr_signal),2);       
        lambda = 5e1;
        cv = 5;
        
        [curr_kernel,curr_predicted_signal,curr_expl_var] = AP_regresskernel( ...
            curr_regressors(use_trials,:)',curr_signal(use_trials,:)', ...
            [],lambda,[false,false],cv);
        
        interarea_kernel(:,curr_t,:) = permute(reshape(curr_kernel,[],length(regressor_areas)),[2,3,1]);
        interarea_expl_var(curr_t) = curr_expl_var.total;
        predicted_signal(use_trials,curr_t) = curr_predicted_signal;
        
    end 
    
    % Get the maximum amplitude kernel with sign
    interarea_kernel_max = max(interarea_kernel,[],3);
    interarea_kernel_min = min(interarea_kernel,[],3);
    interarea_kernel_pos = abs(interarea_kernel_max) > abs(interarea_kernel_min);
    interarea_kernel_magnitude = interarea_kernel_max.*interarea_kernel_pos + ...
        interarea_kernel_min.*~interarea_kernel_pos;
    
    interarea_kernel_mean = mean(interarea_kernel,3);
    
    interarea_kernel_all{curr_animal} = interarea_kernel_mean;
    interarea_expl_var_all{curr_animal} = interarea_expl_var;
  
    predicted_signal_all(unique(trial_id),:,curr_animal) = grpstats(predicted_signal,trial_id);
   
    AP_print_progress_fraction(curr_animal,length(D_all));
    
end

interarea_kernel_mean = nanmean(cat(3,interarea_kernel_all{:}),3);
interarea_expl_var_mean = nanmean(cat(3,interarea_expl_var_all{:}),3);
predicted_signal_all_mean = nanmean(predicted_signal_all,3);

figure;
subplot(3,1,1); hold on;
set(gca,'ColorOrder',[autumn(8);winter(8)]);
plot(t,interarea_kernel_mean')
axis tight
line([0,0],ylim,'color','k');
ylabel('Weight');

subplot(3,1,2);
plot(t,interarea_expl_var_mean,'k','linewidth',2);
axis tight
line([0,0],ylim,'color','k');
ylabel('Explained variance');

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
% plot_move = -1;
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
% plot_conditions = find(conditions(:,1) == 0.06 & ...
%     conditions(:,2) == 1 & ...
%     conditions(:,4) == plot_timing);
% plot_color = zeros(length(plot_conditions),3);
% plot_color(conditions(plot_conditions,3) == -1,1) = 1;
% plot_color(conditions(plot_conditions,3) == 1,3) = 1;

% % (to plot zero-contrast by movement)
% plot_conditions = find(conditions(:,1) == 0 & ...
%     conditions(:,4) == plot_timing);
% plot_color = zeros(length(plot_conditions),3);
% plot_color(conditions(plot_conditions,3) == -1,1) = 1;
% plot_color(conditions(plot_conditions,3) == 1,3) = 1;

% Plot predicted signal
subplot(3,1,3); 
hold on;
for curr_condition_idx = 1:length(plot_conditions);
    curr_condition = plot_conditions(curr_condition_idx);
    plot(t,squeeze(predicted_signal_all_mean(curr_condition,:)), ...
        'color',plot_color(curr_condition_idx,:),'linewidth',2);
end
axis tight;
line([0,0],ylim,'color','k');
xlabel('Time from stim');

% Plot weights by spatial ROI
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

roi_cat = cat(3,wf_roi.mask);

use_t_weights = t > -0.1 & t < 0;

interarea_kernel_mean = nanmean(cat(3,interarea_kernel_all{:}),3);
roi_weights = nanmean(interarea_kernel_mean(:,use_t_weights),2);

figure; hold on
set(gca,'YDir','reverse');
AP_reference_outline('ccf_aligned','k');
for curr_roi_idx = 1:length(regressor_areas)
    curr_roi_boundary = cell2mat(bwboundaries(roi_cat(:,:,regressor_areas(curr_roi_idx))));
    patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),roi_weights(curr_roi_idx));
end
curr_roi_boundary = cell2mat(bwboundaries(roi_cat(:,:,regress_area)));
patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),[0,0,0]);
axis image off;
colormap(colormap_BlueWhiteRed);
caxis([-max(abs(caxis)),max(abs(caxis))]);


%% Load and process wheel regression 

data_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\wheel_regression';
load(data_fn);

% Set times
framerate = 35;
upsample_factor = 2;
sample_rate = framerate*upsample_factor;

kernel_samples = [-40:1:20];
kernel_t = kernel_samples/sample_rate;

% Wheel prediction kernels
str_wheel_kernel = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmean(cat(3,x{:}),3),{batch_vars(:).str_wheel_kernel},'uni',false),[1,3,2])),3);

ctx_wheel_kernel = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmean(cat(4,x{:}),4),{batch_vars(:).ctx_wheel_kernel},'uni',false),[1,3,4,2])),4);

figure;
hold on;
set(gca,'ColorOrder',copper(4));
plot(kernel_t,str_wheel_kernel','linewidth',2);
line([0,0],ylim,'color','k');
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:size(str_wheel_kernel,1)),'uni',false)')
xlabel('Time lag (s)');
ylabel('Weight');

AP_imscroll(ctx_wheel_kernel,kernel_t)
axis image; colormap(colormap_BlueWhiteRed);
caxis([-prctile(abs(ctx_wheel_kernel(:)),99.9),prctile(abs(ctx_wheel_kernel(:)),99.9)]);
AP_reference_outline('ccf_aligned','k');
AP_reference_outline('grid',[0.5,0.5,0.5]);

% Predicted v measured wheel velocity
n_bins = 19;
velocity_prctile = 95;

[~,~,wheel_velocity_bins] = cellfun(@(x) cellfun(@(x) ...
    histcounts(x/prctile(abs(x),velocity_prctile),n_bins),x,'uni',false), ...
    {batch_vars(:).wheel_velocity_resample},'uni',false);

str_predicted_wheel_binned = ...
    cellfun(@(x,y) cellfun(@(x,y) ...
    accumarray(x',(y/prctile(abs(y),velocity_prctile))',[],@median),x,y,'uni',false), ...
    wheel_velocity_bins,{batch_vars(:).str_predicted_wheel},'uni',false);
str_predicted_wheel_binned_mean = ...
    nanmean(cell2mat(cellfun(@(x) nanmean(horzcat(x{:}),2),str_predicted_wheel_binned,'uni',false)),2);

ctx_predicted_wheel_binned = ...
    cellfun(@(x,y) cellfun(@(x,y) ...
    accumarray(x',(y/prctile(abs(y),velocity_prctile))',[],@median),x,y,'uni',false), ...
    wheel_velocity_bins,{batch_vars(:).ctx_predicted_wheel},'uni',false);
ctx_predicted_wheel_binned_mean = ...
    nanmean(cell2mat(cellfun(@(x) nanmean(horzcat(x{:}),2),ctx_predicted_wheel_binned,'uni',false)),2);

figure; hold on;
bin_centers = conv(linspace(-1,1,n_bins+1),[1,1]/2,'valid');
plot(bin_centers,str_predicted_wheel_binned_mean,'r','linewidth',2)
plot(bin_centers,ctx_predicted_wheel_binned_mean,'b','linewidth',2)
line([-1,1],[-1,1],'color','k');
line([0,0],[-1,1],'color','k');
line([-1,1],[0,0],'color','k');
axis square;
xlabel('Measured wheel velocity (max norm)')
ylabel('Predicted wheel velocity (max norm)')
legend({'Striatum','Cortex'});


%% Load and process wheel regression (IMAGING ONLY DAYS)

data_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\wheel_regression_wfonly';
load(data_fn);

% Set times
framerate = 35;
upsample_factor = 2;
sample_rate = framerate*upsample_factor;

kernel_samples = [-40:1:20];
kernel_t = kernel_samples/sample_rate;

% Wheel prediction kernels
ctx_wheel_kernel = nanmean(cell2mat(permute(cellfun(@(x) ...
    nanmean(cat(4,x{:}),4),{batch_vars(:).ctx_wheel_kernel},'uni',false),[1,3,4,2])),4);

AP_imscroll(ctx_wheel_kernel,kernel_t)
axis image; colormap(colormap_BlueWhiteRed);
caxis([-prctile(abs(ctx_wheel_kernel(:)),99.9),prctile(abs(ctx_wheel_kernel(:)),99.9)]);
AP_reference_outline('ccf_aligned','k');
AP_reference_outline('grid',[0.5,0.5,0.5]);

% Kernel reflected and subtract (different - left) and add (same - right)
AP_imscroll([ctx_wheel_kernel - AP_reflect_widefield(ctx_wheel_kernel), ...
     ctx_wheel_kernel + AP_reflect_widefield(ctx_wheel_kernel)],...
     kernel_t);
axis image; colormap(colormap_BlueWhiteRed);
caxis([-prctile(abs(ctx_wheel_kernel(:)),99.9),prctile(abs(ctx_wheel_kernel(:)),99.9)]);
AP_reference_outline('ccf_aligned','k');
AP_reference_outline('grid',[0.5,0.5,0.5]);

% Predicted v measured wheel velocity
n_bins = 19;
velocity_prctile = 95;

[~,~,wheel_velocity_bins] = cellfun(@(x) cellfun(@(x) ...
    histcounts(x/prctile(abs(x),velocity_prctile),n_bins),x,'uni',false), ...
    {batch_vars(:).wheel_velocity_resample},'uni',false);

ctx_predicted_wheel_binned = ...
    cellfun(@(x,y) cellfun(@(x,y) ...
    accumarray(x',(y/prctile(abs(y),velocity_prctile))',[],@median),x,y,'uni',false), ...
    wheel_velocity_bins,{batch_vars(:).ctx_predicted_wheel},'uni',false);
ctx_predicted_wheel_binned_mean = ...
    nanmean(cell2mat(cellfun(@(x) nanmean(horzcat(x{:}),2),ctx_predicted_wheel_binned,'uni',false)),2);

figure; hold on;
bin_centers = conv(linspace(-1,1,n_bins+1),[1,1]/2,'valid');
plot(bin_centers,ctx_predicted_wheel_binned_mean,'b','linewidth',2)
line([-1,1],[-1,1],'color','k');
line([0,0],[-1,1],'color','k');
line([-1,1],[0,0],'color','k');
axis square;
xlabel('Measured wheel velocity (max norm)')
ylabel('Predicted wheel velocity (max norm)')


%% Regress from task events to full activity trace (stim + move choice)

n_aligned_depths = 4;
load_data = 'early';

% Load data (early/late/all)
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
early_data_fn = ['all_trial_activity_df_kernel-str_earlymove_' num2str(n_aligned_depths) '_depths.mat'];
late_data_fn = ['all_trial_activity_df_kernel-str_latemove_' num2str(n_aligned_depths) '_depths.mat'];

switch load_data
    case 'early'
        load([data_path filesep early_data_fn])
    case 'late'
        load([data_path filesep late_data_fn])
    case 'all'
        early_data = load([data_path filesep early_data_fn]);
        late_data = load([data_path filesep late_data_fn]);
        
        n_animals = length(early_data.D_all);
        
        D_all = cell(n_animals,1);
        fluor_all = cell(n_animals,1);
        mua_all = cell(n_animals,1);
        wheel_all = cell(n_animals,1);
        for curr_animal = 1:n_animals
            for curr_day = 1:length(early_data.D_all{curr_animal})
                D_all{curr_animal}{curr_day,1} = ...
                    cell2struct(cellfun(@vertcat,struct2cell(early_data.D_all{curr_animal}{curr_day}), ...
                    struct2cell(late_data.D_all{curr_animal}{curr_day}),'uni',0),fieldnames(early_data.D_all{curr_animal}{curr_day}),1);
                
                fluor_all{curr_animal}{curr_day,1} = ...
                    [early_data.fluor_all{curr_animal}{curr_day}; ...
                    late_data.fluor_all{curr_animal}{curr_day}];
                
                mua_all{curr_animal}{curr_day,1} = ...
                    [early_data.mua_all{curr_animal}{curr_day}; ...
                    late_data.mua_all{curr_animal}{curr_day}];
                
                wheel_all{curr_animal}{curr_day,1} = ...
                    [early_data.wheel_all{curr_animal}{curr_day}; ...
                    late_data.wheel_all{curr_animal}{curr_day}];
                
            end
        end
end

% Get time
framerate = 35;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = framerate*upsample_factor;
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);


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

% Concantenate D
D_cellcat = vertcat(D_all{:});
D_cellcat = struct2cell(vertcat(D_cellcat{:}));
D_allcat = cell2struct(arrayfun(@(x) vertcat(D_cellcat{x,:}),1:size(D_cellcat,1),'uni',false)',fieldnames(D_all{1}{1}));

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_rois,2);
mua_allcat = nan(0,length(t),n_depths,2);
wheel_allcat = nan(0,length(t),2);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

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
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
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
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
    wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
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
    
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat_hemidiff_norm];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];

end

% Downsample
downsample_factor = 3;
t_downsample = linspace(t(1),t(end),round(length(t)/downsample_factor));

wheel = interp1(t,wheel_allcat(:,:,1)',t_downsample)';

%%% Regression (separate regressors to get partial explained variance)

% Set up stim regressors
contrasts = unique(trial_contrast_allcat(trial_contrast_allcat > 0));
contrastsides = sort([-contrasts;contrasts]);

stim_regressors = zeros(size(wheel,1),size(wheel,2),length(contrastsides));
for curr_condition_idx = 1:length(contrastsides)
    stim_regressors(trial_contrast_allcat.*trial_side_allcat == ...
        contrastsides(curr_condition_idx),find(t_downsample > 0,1),curr_condition_idx) = 1;
end

% Choice regressors: one at time of stim, one at time of move
choice_stim_regressors = zeros(size(wheel,1),size(wheel,2),2);
choice_stim_regressors(trial_choice_allcat == -1,find(t_downsample > 0,1),1) = 1;
choice_stim_regressors(trial_choice_allcat == 1,find(t_downsample > 0,1),2) = 1;

[~,move_idx] = max(abs(wheel(:,:,1)) > 2,[],2);
choice_move_regressors = zeros(size(wheel,1),size(wheel,2),2);
for curr_trial = 1:size(choice_move_regressors,1)
    curr_choice = trial_choice_allcat(curr_trial);
    choice_move_regressors(curr_trial,move_idx(curr_trial),curr_choice/2+1.5) = 1;
end

% Wheel regressors - separate leftward/rightward velocity
wheel_vel_norm = [zeros(size(wheel,1),1),diff(wheel(:,:,1),[],2)];
wheel_vel_norm = wheel_vel_norm/prctile(abs(wheel_vel_norm(:)),95);
wheel_vel_norm(abs(wheel_vel_norm) > 1) = sign(wheel_vel_norm(abs(wheel_vel_norm) > 1));
wheel_regressors = abs(repmat(wheel_vel_norm,1,1,2).*cat(3,wheel_vel_norm < 0,wheel_vel_norm > 0)); 

regressors = {stim_regressors;choice_stim_regressors;choice_move_regressors;wheel_regressors};

% Do regression 
t_shifts = {[0:30];[0:20];[-40:20];[-20:20]};
lambda = 1e1;
zs = [false,false];
cvfold = 5;

% Set activity 
activity = mua_allcat(:,:,2,1);
activity = interp1(t,activity',t_downsample)';

activity_reshape = reshape(activity',[],1)';
regressors_reshape = cellfun(@(x) ...
    reshape(permute(x,[2,1,3]),[],size(x,3))',regressors,'uni',false);

[kernel_flat,activity_predicted_reshape,expl_var,activity_predicted_reduced_reshape] = AP_regresskernel( ...
    cellfun(@(x) x(:,~isnan(activity_reshape)),regressors_reshape,'uni',false), ...
    activity_reshape(~isnan(activity_reshape)),t_shifts,lambda,zs,cvfold);

activity_predicted = nan(size(activity))';
% (to use full model)
activity_predicted(~isnan(activity')) = activity_predicted_reshape;
% (to use reduced model)
% activity_predicted(~isnan(activity')) = activity_predicted_reduced_reshape(:,:,1);
activity_predicted = activity_predicted';

kernel = cellfun(@(x,t) reshape(x,[],length(t)), ...
    mat2cell(kernel_flat,cellfun(@(x) size(x,1),regressors_reshape).*cellfun(@length,t_shifts),1), ...
    t_shifts,'uni',false);

% Plot kernels
col = colormap_BlueWhiteRed(5);
col(6,:) = [];
figure; hold on;
set(gca,'ColorOrder',[col; ...
    0.8,0,0.8; ...
    0,0.8,0; ...
    0.5,0,0.5; ...
    0,0.5,0; ...
    0,0,0; ...
    0.5,0.5,0.5]);
cellfun(@(k,t) plot(t*(downsample_factor/sample_rate),k'),kernel,t_shifts,'uni',false);
xlabel('Lag time');
ylabel('Weight (\beta)');
legend([cellfun(@(x) ['Stim: ' num2str(x)],num2cell(contrastsides),'uni',false); ...
    {'Choice L (stim-aligned)';'Choice R (stim-aligned)'; ...
    'Choice L (move-aligned)';'Choice R (move-aligned)'; ...
    'Wheel velocity left';'Wheel velocity right'}]);

% Plot activity and predicted by condition
sidecontrasts = unique(bsxfun(@times,trial_contrast_allcat,trial_side_allcat));

activity_split = cell(length(sidecontrasts),2);
activity_split(:,1) = arrayfun(@(x) ...
    activity(trial_side_allcat.*trial_contrast_allcat == sidecontrasts(x) & ...
    trial_choice_allcat == -1,:), ...
    1:length(sidecontrasts),'uni',false);
activity_split(:,2) = arrayfun(@(x) ...
    activity(trial_side_allcat.*trial_contrast_allcat == sidecontrasts(x) & ...
    trial_choice_allcat == 1,:), ...
    1:length(sidecontrasts),'uni',false);
activity_split_mean = cellfun(@(x) nanmean(x,1),activity_split,'uni',false);

activity_predicted_split = cell(length(sidecontrasts),2);
activity_predicted_split(:,1) = arrayfun(@(x) ...
    activity_predicted(trial_side_allcat.*trial_contrast_allcat == sidecontrasts(x) & ...
    trial_choice_allcat == -1,:), ...
    1:length(sidecontrasts),'uni',false);
activity_predicted_split(:,2) = arrayfun(@(x) ...
    activity_predicted(trial_side_allcat.*trial_contrast_allcat == sidecontrasts(x) & ...
    trial_choice_allcat == 1,:), ...
    1:length(sidecontrasts),'uni',false);
activity_predicted_split_mean = cellfun(@(x) nanmean(x,1),activity_predicted_split,'uni',false);

figure; 
col = colormap_BlueWhiteRed(5);
col(6,:) = [0,0,0];

subplot(2,2,1); hold on
set(gca,'ColorOrder',col);
plot(t_downsample,vertcat(activity_split_mean{:,1}),'linewidth',2);
title('Move left');

subplot(2,2,2); hold on
set(gca,'ColorOrder',col);
plot(t_downsample,vertcat(activity_split_mean{:,2}),'linewidth',2);
title('Move right');

subplot(2,2,3); hold on
set(gca,'ColorOrder',col);
plot(t_downsample,vertcat(activity_predicted_split_mean{:,1}),'linewidth',2);
title('Move left (predicted');

subplot(2,2,4); hold on
set(gca,'ColorOrder',col);
plot(t_downsample,vertcat(activity_predicted_split_mean{:,2}),'linewidth',2);
title('Move right (predicted)');

%%% To plot predicted activity by reaction time
move_t = t_downsample(move_idx);

activity_predicted = nan(size(activity))';
% (to use full model)
activity_predicted(~isnan(activity')) = activity_predicted_reshape;
% (to use reduced model)
% activity_predicted(~isnan(activity')) = activity_predicted_reduced_reshape(:,:,4);
activity_predicted = activity_predicted';

% Plot average activity within condition restricted by reaction time
rxn_times = linspace(0,0.5,6);

sides = [-1,1];

activity_split = arrayfun(@(x) ...
    activity_predicted( ... 
    move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
    trial_side_allcat == 1 & ...
    trial_contrast_allcat == 1 & ...
    trial_choice_allcat == -1,:), ...
    1:length(rxn_times)-1,'uni',false);

activity_split_mean = cell2mat(cellfun(@(x) nanmean(x,1),activity_split','uni',false));

figure; hold on;
plot_col = copper(length(activity_split));
set(gca,'ColorOrder',plot_col);
plot(t_downsample,activity_split_mean,'linewidth',2);
legend(arrayfun(@(x) [num2str(rxn_times(x)) '-' num2str(rxn_times(x+1))], ...
    1:length(rxn_times)-1,'uni',false));

for i = 1:length(rxn_times)-1
   line(repmat(rxn_times(i),2,1),ylim,'color',plot_col(i,:));
end

%% Regress from task events to full activity trace (vis/non-vis choice)

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_df_kernel-str_' num2str(n_aligned_depths) '_depths.mat'];
load([data_path filesep data_fn]);

% Get time
framerate = 35;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = framerate*upsample_factor;
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);


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

% Concantenate D
D_cellcat = vertcat(D_all{:});
D_cellcat = struct2cell(vertcat(D_cellcat{:}));
D_allcat = cell2struct(arrayfun(@(x) vertcat(D_cellcat{x,:}),1:size(D_cellcat,1),'uni',false)',fieldnames(D_all{1}{1}));

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_rois,2);
fluor_unilateral_allcat = nan(0,length(t),n_rois,2);
mua_allcat = nan(0,length(t),n_depths,2);
wheel_allcat = nan(0,length(t),2);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

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
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
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
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
    wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
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
    
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat_hemidiff_norm];
    fluor_unilateral_allcat = [fluor_unilateral_allcat;fluor_cat_norm];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];

end

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 1,[],2);
move_t = t(move_idx)';

% Concatenate unilateral 
fluor_unilateral_allcat = cat(1,fluor_unilateral_allcat(:,:,1:length(wf_roi),:), ...
    fluor_unilateral_allcat(:,:,length(wf_roi)+1:end,:));

% Downsample
downsample_factor = 3;
t_downsample = linspace(t(1),t(end),round(length(t)/downsample_factor));

wheel = interp1(t,wheel_allcat(:,:,1)',t_downsample)';

%%% Regression (separate regressors to get partial explained variance)

% Set up stim regressors
contrasts = unique(trial_contrast_allcat(trial_contrast_allcat > 0));
contrastsides = sort([-contrasts;contrasts]);

stim_regressors = zeros(size(wheel,1),size(wheel,2),length(contrastsides));
for curr_condition_idx = 1:length(contrastsides)
    stim_regressors(trial_contrast_allcat.*trial_side_allcat == ...
        contrastsides(curr_condition_idx),find(t_downsample > 0,1),curr_condition_idx) = 1;
end

% Choice regressors: one at time of move (L/R and vis/non-vis)
choice_stim_regressors = zeros(size(wheel,1),size(wheel,2),2);
choice_stim_regressors(trial_choice_allcat == -1,find(t_downsample > 0,1),1) = 1;
choice_stim_regressors(trial_choice_allcat == 1,find(t_downsample > 0,1),2) = 1;

[~,move_idx] = max(abs(wheel(:,:,1)) > 2,[],2);
choice_move_regressors = zeros(size(wheel,1),size(wheel,2),4);
for curr_trial = 1:size(choice_move_regressors,1)
    
    if trial_choice_allcat(curr_trial) == 1 && ...
            (trial_contrast_allcat(curr_trial) == 0 || trial_side_allcat(curr_trial) == 1)
        choice_move_regressors(curr_trial,move_idx(curr_trial),1) = 1;
        
    elseif trial_choice_allcat(curr_trial) == 1 && ...
            trial_contrast_allcat(curr_trial) > 0 && trial_side_allcat(curr_trial) == -1
        choice_move_regressors(curr_trial,move_idx(curr_trial),2) = 1;
        
    elseif trial_choice_allcat(curr_trial) == -1 && ...
            (trial_contrast_allcat(curr_trial) == 0 || trial_side_allcat(curr_trial) == -1)
        choice_move_regressors(curr_trial,move_idx(curr_trial),3) = 1;
        
    elseif trial_choice_allcat(curr_trial) == -1 && ...
            trial_contrast_allcat(curr_trial) > 0 && trial_side_allcat(curr_trial) == 1
        choice_move_regressors(curr_trial,move_idx(curr_trial),4) = 1;
        
    end
end

% Wheel regressors - separate leftward/rightward velocity
wheel_vel_norm = [zeros(size(wheel,1),1),diff(wheel(:,:,1),[],2)];
wheel_vel_norm = wheel_vel_norm/prctile(abs(wheel_vel_norm(:)),95);
wheel_vel_norm(abs(wheel_vel_norm) > 1) = sign(wheel_vel_norm(abs(wheel_vel_norm) > 1));
wheel_regressors = abs(repmat(wheel_vel_norm,1,1,2).*cat(3,wheel_vel_norm < 0,wheel_vel_norm > 0)); 

regressors = {stim_regressors;choice_move_regressors;wheel_regressors};

% Set trials to use
use_trials = move_t > 0 & move_t < 0.5;

% Do regression on MUA
t_shifts = {[0:30];[-40:40];[-20:20]};
lambda = 0;
zs = [false,true];
cvfold = 5;

mua_kernel = cell(3,n_depths);
mua_expl_var = cell(3,n_depths);
for curr_depth = 1:n_depths
    
    activity = mua_allcat(use_trials,:,curr_depth,1);
    activity = interp1(t,activity',t_downsample)';
    
    activity_reshape = reshape(activity',[],1)';
    regressors_reshape = cellfun(@(x) ...
        reshape(permute(x(use_trials,:,:),[2,1,3]),[],size(x,3))',regressors,'uni',false);
    
    [kernel_flat,activity_predicted_reshape,expl_var,activity_predicted_reduced_reshape] = AP_regresskernel( ...
        cellfun(@(x) x(:,~isnan(activity_reshape)),regressors_reshape,'uni',false), ...
        activity_reshape(~isnan(activity_reshape)),t_shifts,lambda,zs,cvfold);
    
    activity_predicted = nan(size(activity))';
    % (to use full model)
    activity_predicted(~isnan(activity')) = activity_predicted_reshape;
    % (to use reduced model)
%     activity_predicted(~isnan(activity')) = activity_predicted_reduced_reshape(:,:,1);
    activity_predicted = activity_predicted';
    
    mua_kernel(:,curr_depth) = cellfun(@(x,t) reshape(x,[],length(t)), ...
        mat2cell(kernel_flat,cellfun(@(x) size(x,1),regressors_reshape).*cellfun(@length,t_shifts),1), ...
        t_shifts,'uni',false);
    
    AP_print_progress_fraction(curr_depth,n_depths);
    
end

figure;
for curr_regressor = 1:3
    for curr_depth = 1:n_depths
        subplot(3,n_depths,(curr_regressor-1)*n_depths+curr_depth)
        imagesc(mua_kernel{curr_regressor,curr_depth})
        colormap(hot);
    end
end

move_start_k = cat(3,mua_kernel{2,:});
wheel_k = cat(3,mua_kernel{3,:});

figure; 

p1 = subplot(1,5,1); hold on;
set(gca,'ColorOrder',copper(4));
plot(t_shifts{2}/sample_rate/downsample_factor, ...
    squeeze(move_start_k(4,:,:)),'linewidth',2)
line([0,0],ylim,'color','k');
ylabel('Weight')
title('Move L Vis')

p2 = subplot(1,5,2); hold on;
set(gca,'ColorOrder',copper(4));
plot(t_shifts{2}/sample_rate/downsample_factor, ...
    squeeze(nanmean(move_start_k([4],:,:) - ...
    move_start_k([3],:,:),1)),'linewidth',2)
line([0,0],ylim,'color','k');
ylabel('Weight')
title('Vis L - Non-vis L')

p3 = subplot(1,5,3); hold on;
set(gca,'ColorOrder',copper(4));
plot(t_shifts{2}/sample_rate/downsample_factor, ...
    squeeze(nanmean(move_start_k(3,:,:) - ...
    move_start_k(1,:,:),1)),'linewidth',2)
line([0,0],ylim,'color','k');
ylabel('Weight')
title('Non-vis L - Non-vis R')

p4 = subplot(1,5,4); hold on;
set(gca,'ColorOrder',copper(4));
plot(t_shifts{2}/sample_rate/downsample_factor, ...
    squeeze(nanmean(move_start_k(4,:,:) - ...
    move_start_k(2,:,:),1)),'linewidth',2)
line([0,0],ylim,'color','k');
ylabel('Weight')
title('Vis L - Vis R')

p5 = subplot(1,5,5); hold on;
set(gca,'ColorOrder',copper(4));
plot(t_shifts{3}/sample_rate/downsample_factor, ...
    squeeze(nanmean(wheel_k(1,:,:) - ...
    wheel_k(2,:,:),1)),'linewidth',2)
line([0,0],ylim,'color','k');
ylabel('Weight')
title('Wheel L - Wheel R')

linkaxes([p1,p2,p3,p4,p5],'y');

% Do regression on fluor (not combined unilateral at the moment)
t_shifts = {[0:30];[-40:40];[-20:20]};
lambda = 0;
zs = [false,true];
cvfold = 5;

fluor_kernel = cell(3,n_rois/2);
fluor_expl_var = cell(3,n_rois/2);
for curr_roi = 1:n_rois/2
    
    activity = fluor_allcat(use_trials,:,curr_roi,1);
    activity = interp1(t,activity',t_downsample)';
    
    activity_reshape = reshape(activity',[],1)';
    regressors_reshape = cellfun(@(x) ...
        reshape(permute(x(use_trials,:,:),[2,1,3]),[],size(x,3))',regressors,'uni',false);
    
    [kernel_flat,activity_predicted_reshape,expl_var,activity_predicted_reduced_reshape] = AP_regresskernel( ...
        cellfun(@(x) x(:,~isnan(activity_reshape)),regressors_reshape,'uni',false), ...
        activity_reshape(~isnan(activity_reshape)),t_shifts,lambda,zs,cvfold);
    
    activity_predicted = nan(size(activity))';
    % (to use full model)
    activity_predicted(~isnan(activity')) = activity_predicted_reshape;
    % (to use reduced model)
%     activity_predicted(~isnan(activity')) = activity_predicted_reduced_reshape(:,:,1);
    activity_predicted = activity_predicted';
    
    fluor_kernel(:,curr_roi) = cellfun(@(x,t) reshape(x,[],length(t)), ...
        mat2cell(kernel_flat,cellfun(@(x) size(x,1),regressors_reshape).*cellfun(@length,t_shifts),1), ...
        t_shifts,'uni',false);
    
    AP_print_progress_fraction(curr_roi,n_rois/2);
    
end

figure;
for curr_regressor = 1:3
    for curr_roi = 1:n_rois/2
        subplot(3,n_rois/2,(curr_regressor-1)*n_rois/2+curr_roi)
        imagesc(fluor_kernel{curr_regressor,curr_roi})
        colormap(hot);
        
        if curr_regressor == 1
            title(wf_roi(curr_roi).area)
        end
    end
end

move_start_k = cat(3,fluor_kernel{2,:});
wheel_k = cat(3,fluor_kernel{3,:});

figure; 
p1 = subplot(1,5,1); hold on;
set(gca,'ColorOrder',jet(n_rois/2));
plot(t_shifts{2}/sample_rate/downsample_factor, ...
    squeeze(move_start_k(3,:,:)),'linewidth',2)
line([0,0],ylim,'color','k');
ylabel('Weight')
title('Move L Vis')

p2 = subplot(1,5,2); hold on;
set(gca,'ColorOrder',jet(n_rois/2));
plot(t_shifts{2}/sample_rate/downsample_factor, ...
    squeeze(nanmean(move_start_k([4],:,:) - ...
    move_start_k([3],:,:),1)),'linewidth',2)
line([0,0],ylim,'color','k');
ylabel('Weight')
title('Vis L - Non-vis L')

p3 = subplot(1,5,3); hold on;
set(gca,'ColorOrder',jet(n_rois/2));
plot(t_shifts{2}/sample_rate/downsample_factor, ...
    squeeze(nanmean(move_start_k(3,:,:) - ...
    move_start_k(1,:,:),1)),'linewidth',2)
line([0,0],ylim,'color','k');
ylabel('Weight')
title('Non-vis L - Non-vis R')

p4 = subplot(1,5,4); hold on;
set(gca,'ColorOrder',jet(n_rois/2));
plot(t_shifts{2}/sample_rate/downsample_factor, ...
    squeeze(nanmean(move_start_k(4,:,:) - ...
    move_start_k(2,:,:),1)),'linewidth',2)
line([0,0],ylim,'color','k');
ylabel('Weight')
title('Vis L - Vis R')

p5 = subplot(1,5,5); hold on;
set(gca,'ColorOrder',jet(n_rois/2));
plot(t_shifts{3}/sample_rate/downsample_factor, ...
    squeeze(nanmean(wheel_k(1,:,:) - ...
    wheel_k(2,:,:),1)),'linewidth',2)
line([0,0],ylim,'color','k');
ylabel('Weight')
title('Wheel L - Wheel R')

linkaxes([p1,p2,p3,p4,p5],'y');

%% Regress from task events to full activity trace (fluor unilateral)

n_aligned_depths = 4;
load_data = 'late';

% Load data (early/late/all)
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
early_data_fn = ['all_trial_activity_df_kernel-str_earlymove_' num2str(n_aligned_depths) '_depths.mat'];
late_data_fn = ['all_trial_activity_df_kernel-str_latemove_' num2str(n_aligned_depths) '_depths.mat'];

switch load_data
    case 'early'
        load([data_path filesep early_data_fn])
    case 'late'
        load([data_path filesep late_data_fn])
    case 'all'
        early_data = load([data_path filesep early_data_fn]);
        late_data = load([data_path filesep late_data_fn]);
        
        n_animals = length(early_data.D_all);
        
        D_all = cell(n_animals,1);
        fluor_all = cell(n_animals,1);
        mua_all = cell(n_animals,1);
        wheel_all = cell(n_animals,1);
        for curr_animal = 1:n_animals
            for curr_day = 1:length(early_data.D_all{curr_animal})
                D_all{curr_animal}{curr_day,1} = ...
                    cell2struct(cellfun(@vertcat,struct2cell(early_data.D_all{curr_animal}{curr_day}), ...
                    struct2cell(late_data.D_all{curr_animal}{curr_day}),'uni',0),fieldnames(early_data.D_all{curr_animal}{curr_day}),1);
                
                fluor_all{curr_animal}{curr_day,1} = ...
                    [early_data.fluor_all{curr_animal}{curr_day}; ...
                    late_data.fluor_all{curr_animal}{curr_day}];
                
                mua_all{curr_animal}{curr_day,1} = ...
                    [early_data.mua_all{curr_animal}{curr_day}; ...
                    late_data.mua_all{curr_animal}{curr_day}];
                
                wheel_all{curr_animal}{curr_day,1} = ...
                    [early_data.wheel_all{curr_animal}{curr_day}; ...
                    late_data.wheel_all{curr_animal}{curr_day}];
                
            end
        end
end

% Get time
framerate = 35;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = framerate*upsample_factor;
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);


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

% Concantenate D
D_cellcat = vertcat(D_all{:});
D_cellcat = struct2cell(vertcat(D_cellcat{:}));
D_allcat = cell2struct(arrayfun(@(x) vertcat(D_cellcat{x,:}),1:size(D_cellcat,1),'uni',false)',fieldnames(D_all{1}{1}));

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_rois,2);
fluor_unilateral_allcat = nan(0,length(t),n_rois,2);
mua_allcat = nan(0,length(t),n_depths,2);
wheel_allcat = nan(0,length(t),2);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

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
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
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
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
    wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
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
    
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat_hemidiff_norm];
    fluor_unilateral_allcat = [fluor_unilateral_allcat;fluor_cat_norm];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];

end

% Downsample
downsample_factor = 3;
t_downsample = linspace(t(1),t(end),round(length(t)/downsample_factor));

wheel = interp1(t,wheel_allcat(:,:,1)',t_downsample)';

% Concatenate unilateral 
fluor_unilateral_allcat = cat(1,fluor_unilateral_allcat(:,:,1:length(wf_roi),:), ...
    fluor_unilateral_allcat(:,:,length(wf_roi)+1:end,:));

% % (if concatenating side - messy, just here for the moment)
trial_contrast_allcat = [trial_contrast_allcat;trial_contrast_allcat];
trial_side_allcat = [trial_side_allcat;-trial_side_allcat];
trial_choice_allcat = [trial_choice_allcat;-trial_choice_allcat];
wheel = [wheel;-wheel];
n_rois = n_rois/2;

%%% Regression (separate regressors to get partial explained variance)

% Set up stim regressors
contrasts = unique(trial_contrast_allcat(trial_contrast_allcat > 0));
contrastsides = sort([-contrasts;contrasts]);

stim_regressors = zeros(size(wheel,1),size(wheel,2),length(contrastsides));
for curr_condition_idx = 1:length(contrastsides)
    stim_regressors(trial_contrast_allcat.*trial_side_allcat == ...
        contrastsides(curr_condition_idx),find(t_downsample > 0,1),curr_condition_idx) = 1;
end

% Choice regressors: one at time of move (L/R and vis/non-vis)
choice_stim_regressors = zeros(size(wheel,1),size(wheel,2),2);
choice_stim_regressors(trial_choice_allcat == -1,find(t_downsample > 0,1),1) = 1;
choice_stim_regressors(trial_choice_allcat == 1,find(t_downsample > 0,1),2) = 1;

[~,move_idx] = max(abs(wheel(:,:,1)) > 2,[],2);
choice_move_regressors = zeros(size(wheel,1),size(wheel,2),4);
for curr_trial = 1:size(choice_move_regressors,1)
    
    if trial_choice_allcat(curr_trial) == 1 && ...
            (trial_contrast_allcat(curr_trial) == 0 || trial_side_allcat(curr_trial) == 1)
        choice_move_regressors(curr_trial,move_idx(curr_trial),1) = 1;
        
    elseif trial_choice_allcat(curr_trial) == 1 && ...
            trial_contrast_allcat(curr_trial) > 0 && trial_side_allcat(curr_trial) == -1
        choice_move_regressors(curr_trial,move_idx(curr_trial),2) = 1;
        
    elseif trial_choice_allcat(curr_trial) == -1 && ...
            (trial_contrast_allcat(curr_trial) == 0 || trial_side_allcat(curr_trial) == -1)
        choice_move_regressors(curr_trial,move_idx(curr_trial),3) = 1;
        
    elseif trial_choice_allcat(curr_trial) == -1 && ...
            trial_contrast_allcat(curr_trial) > 0 && trial_side_allcat(curr_trial) == 1
        choice_move_regressors(curr_trial,move_idx(curr_trial),4) = 1;
        
    end
end

% Wheel regressors - separate leftward/rightward velocity
wheel_vel_norm = [zeros(size(wheel,1),1),diff(wheel(:,:,1),[],2)];
wheel_vel_norm = wheel_vel_norm/prctile(abs(wheel_vel_norm(:)),95);
wheel_vel_norm(abs(wheel_vel_norm) > 1) = sign(wheel_vel_norm(abs(wheel_vel_norm) > 1));
wheel_regressors = abs(repmat(wheel_vel_norm,1,1,2).*cat(3,wheel_vel_norm < 0,wheel_vel_norm > 0)); 

regressors = {stim_regressors;choice_move_regressors;wheel_regressors};

% Do regression on fluor (not combined unilateral at the moment)
t_shifts = {[0:30];[-10:10];[-10:5]};
lambda = 1e2;
zs = [true,true];
cvfold = 5;

fluor_kernel = cell(3,n_rois);
fluor_expl_var = cell(3,n_rois);
for curr_roi = 1:n_rois
    
    activity = fluor_unilateral_allcat(:,:,curr_roi,1);
    activity = interp1(t,activity',t_downsample)';
    
    activity_reshape = reshape(activity',[],1)';
    regressors_reshape = cellfun(@(x) ...
        reshape(permute(x,[2,1,3]),[],size(x,3))',regressors,'uni',false);
    
    [kernel_flat,activity_predicted_reshape,expl_var,activity_predicted_reduced_reshape] = AP_regresskernel( ...
        cellfun(@(x) x(:,~isnan(activity_reshape)),regressors_reshape,'uni',false), ...
        activity_reshape(~isnan(activity_reshape)),t_shifts,lambda,zs,cvfold);
    
    activity_predicted = nan(size(activity))';
    % (to use full model)
    activity_predicted(~isnan(activity')) = activity_predicted_reshape;
    % (to use reduced model)
%     activity_predicted(~isnan(activity')) = activity_predicted_reduced_reshape(:,:,1);
    activity_predicted = activity_predicted';
    
    fluor_kernel(:,curr_roi) = cellfun(@(x,t) reshape(x,[],length(t)), ...
        mat2cell(kernel_flat,cellfun(@(x) size(x,1),regressors_reshape).*cellfun(@length,t_shifts),1), ...
        t_shifts,'uni',false);
    
    AP_print_progress_fraction(curr_roi,n_rois);
    
end

figure;
for curr_regressor = 1:3
    for curr_roi = 1:n_rois
        subplot(3,n_rois,(curr_regressor-1)*n_rois+curr_roi)
        imagesc(fluor_kernel{curr_regressor,curr_roi})
        colormap(hot);
    end
end

move_start_k = cat(3,fluor_kernel{2,:});
wheel_k = cat(3,fluor_kernel{3,:});

figure; 
subplot(1,3,1); hold on;
set(gca,'ColorOrder',jet(n_rois));
plot(t_shifts{2}/sample_rate/downsample_factor, ...
    squeeze(nanmean(move_start_k([2,4],:,:) - ...
    move_start_k([1,3],:,:),1)),'linewidth',2)
line([0,0],ylim,'color','k');
ylabel('Weight')
title('Vis - Non-vis')

subplot(1,3,2); hold on;
set(gca,'ColorOrder',jet(n_rois));
plot(t_shifts{2}/sample_rate/downsample_factor, ...
    squeeze(nanmean(move_start_k(4,:,:) - ...
    move_start_k(2,:,:),1)),'linewidth',2)
line([0,0],ylim,'color','k');
ylabel('Weight')
title('Vis L - Vis R')

subplot(1,3,3); hold on;
set(gca,'ColorOrder',jet(n_rois));
plot(t_shifts{3}/sample_rate/downsample_factor, ...
    squeeze(nanmean(wheel_k(1,:,:) - ...
    wheel_k(2,:,:),1)),'linewidth',2)
line([0,0],ylim,'color','k');
ylabel('Weight')
title('Wheel L - Wheel R')

%% Regress from task events to activity (V,long)
% this only uses stim, move onset, move ongoing, reward
% generates predicted activity to use in other analyses
%%%% TO DO HERE: low-pass 6hz?

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_Udf_kernel-str_4_depths_long'];

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);
reward_all = cellfun(@(data,use_expts) data(use_expts),reward_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);
reward_all = reward_all(use_animals);

% Get widefield ROIs
n_rois = 200;

% MUA depths
n_depths = size(mua_all{1}{1},3);

% Set up conditions
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_rois,1);
mua_allcat = nan(0,length(t),n_depths,1);
wheel_allcat = nan(0,length(t),1);
reward_allcat = nan(0,length(t),1);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

D_cat = struct('stimulus',cell(1,1),'response',cell(1,1),'repeatNum',cell(1,1));

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    n_align = size(fluor_cat,4);
    
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x,[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
%     wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
    % Concatenate behavioural data
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
    D.response = cell2mat(cellfun(@(x) x.response,D_all{curr_animal},'uni',false));
    D.repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all{curr_animal},'uni',false));
    
    D.day = trial_day;
    
    D_cat.stimulus = vertcat(D_cat.stimulus,D.stimulus);
    D_cat.response = vertcat(D_cat.response,D.response);
    D_cat.repeatNum = vertcat(D_cat.repeatNum,D.repeatNum);
    
    % Get trial ID   
    trial_contrast = max(D.stimulus,[],2);
    [~,side_idx] = max(D.stimulus > 0,[],2);
    trial_side = (side_idx-1.5)*2;
    trial_choice = -(D.response-1.5)*2;
    
    trial_conditions = ...
        [trial_contrast, trial_side, ...
        trial_choice, ones(size(trial_day))];
    [~,trial_id] = ismember(trial_conditions,conditions,'rows');
    
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    reward_allcat = [reward_allcat;vertcat(reward_all{curr_animal}{:})];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];
   
end

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Get maximum wheel velocity in chosen direction
wheel_velocity_allcat = [zeros(size(wheel_allcat,1),1,1),diff(wheel_allcat,[],2)];
[max_speed,max_vel_idx] = max(abs(wheel_velocity_allcat(:,:,1).* ...
    (bsxfun(@times,wheel_velocity_allcat(:,:,1),trial_choice_allcat) > 0)),[],2);
vel_t = t(max_vel_idx);

max_vel = max_speed.*trial_choice_allcat;

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Downsample (otherwise it's too much for regression) and d(smooth(fluor))
downsample_factor = 4;
t_downsample = linspace(t(1),t(end),round(length(t)/downsample_factor));

smooth_factor = 3;

t_diff =  conv(t,[1,1]/2,'valid');
t_downsample_diff = conv(t_downsample,[1,1]/2,'valid');

fluor_allcat_downsamp_diff = permute(interp1(t_diff,permute(diff(convn( ...
    fluor_allcat,ones(1,smooth_factor)/smooth_factor,'same'),[],2),[2,1,3,4]),t_downsample_diff),[2,1,3,4]);

mua_allcat_downsamp = permute(interp1(t,permute(mua_allcat,[2,1,3,4]),t_downsample_diff),[2,1,3,4]);

wheel = interp1(t,wheel_allcat(:,:,1)',t_downsample_diff)';

%%% Regression (separate regressors to get partial explained variance)

% Stim regressors
contrasts = unique(trial_contrast_allcat(trial_contrast_allcat > 0));
contrastsides = sort([-contrasts;contrasts]);

stim_regressors = zeros(size(wheel,1),size(wheel,2),length(contrastsides));
for curr_condition_idx = 1:length(contrastsides)
    stim_regressors(trial_contrast_allcat.*trial_side_allcat == ...
        contrastsides(curr_condition_idx),find(t_downsample_diff > 0,1),curr_condition_idx) = 1;
end

% Move onset regressors (L/R)
[~,move_idx] = max(abs(wheel(:,:,1)) > 2,[],2);
move_onset_regressors = zeros(size(wheel,1),size(wheel,2),2);
for curr_trial = 1:size(move_onset_regressors,1)
    
    % To use binary
    if trial_choice_allcat(curr_trial) == -1
        move_onset_regressors(curr_trial,move_idx(curr_trial),1) = 1;        
    elseif trial_choice_allcat(curr_trial) == 1
        move_onset_regressors(curr_trial,move_idx(curr_trial),2) = 1;         
    end
    
%     % To fold in maximum velocity in chosen direction
%     if trial_choice_allcat(curr_trial) == -1
%         move_onset_regressors(curr_trial,move_idx(curr_trial),1) = abs(max_vel(curr_trial));
%     elseif trial_choice_allcat(curr_trial) == 1
%         move_onset_regressors(curr_trial,move_idx(curr_trial),2) = abs(max_vel(curr_trial));
%     end    
    
end

% Move ongoing regressor - 1 if ongoing movement in that direction
wheel_vel = [zeros(size(wheel,1),1),diff(wheel(:,:,1),[],2)];
move_ongoing_regressors = zeros(size(wheel,1),size(wheel,2),2);
move_ongoing_regressors(:,:,1) = wheel_vel < -2;
move_ongoing_regressors(:,:,2) = wheel_vel > 2;

% Wheel regressors - separate leftward/rightward velocity
wheel_vel_norm = [zeros(size(wheel,1),1),diff(wheel(:,:,1),[],2)];
wheel_vel_norm = wheel_vel_norm/prctile(abs(wheel_vel_norm(:)),95);
wheel_vel_norm(abs(wheel_vel_norm) > 1) = sign(wheel_vel_norm(abs(wheel_vel_norm) > 1));
wheel_regressors = abs(repmat(wheel_vel_norm,1,1,2).*cat(3,wheel_vel_norm < 0,wheel_vel_norm > 0)); 

% Go cue regressors - separate for early/late move
go_cue_regressors = zeros(size(wheel,1),size(wheel,2));
go_cue_regressors(move_t <= 0.5,find(t_downsample_diff > 0.5,1),1) = 1;
go_cue_regressors(move_t > 0.5,find(t_downsample_diff > 0.5,1),2) = 1;

% Reward regressors
reward_allcat_downsamp = permute(interp1(t,permute(reward_allcat,[2,1,3,4]),t_downsample_diff,'nearest'),[2,1,3,4]);

reward_allcat_regressor = zeros(size(reward_allcat,1),length(t_downsample_diff));
for curr_trial = 1:size(reward_allcat,1)
   curr_reward = find(reward_allcat(curr_trial,:));
   for i = curr_reward
       curr_reward_t = t(i);
       reward_allcat_regressor(curr_trial,find(t_downsample_diff > curr_reward_t,1)) = 1;
   end
end

regressors = {stim_regressors;move_onset_regressors;move_ongoing_regressors;go_cue_regressors;reward_allcat_regressor};
regressor_labels = {'Stim','Move onset','Move ongoing','Go cue','Reward'};

% Set regression parameters
regress_trials = true(size(move_t));
t_shifts = {[-0.1,0.6];[-0.2,0.5];[0,0];[-0.05,0.3];[-0.1,0.5]};

sample_shifts = cellfun(@(x) round(x(1)*(sample_rate/downsample_factor)): ...
    round(x(2)*(sample_rate/downsample_factor)),t_shifts,'uni',false);
lambda = 0;
zs = [false,false];
cvfold = 5;
return_constant = true;


%%% Do regression on fluor
fluor_allcat_predicted = nan(size(fluor_allcat_downsamp_diff));

fluor_kernel = cell(length(regressors)+1,n_rois);
fluor_expl_var = cell(length(regressors),n_rois);
for curr_roi = 1:n_rois
    
    activity = fluor_allcat_downsamp_diff(regress_trials,:,curr_roi,1);
    
    activity_reshape = reshape(activity',[],1)';
    regressors_reshape = cellfun(@(x) ...
        reshape(permute(x(regress_trials,:,:),[2,1,3]),[],size(x,3))',regressors,'uni',false);
    
    [kernel_flat,activity_predicted_reshape,expl_var,activity_predicted_reduced_reshape] = AP_regresskernel( ...
        cellfun(@(x) x(:,~isnan(activity_reshape)),regressors_reshape,'uni',false), ...
        activity_reshape(~isnan(activity_reshape)),sample_shifts,lambda,zs,cvfold,return_constant);
    
    activity_predicted = nan(size(activity))';
    % (to use full model)
    activity_predicted(~isnan(activity')) = activity_predicted_reshape;
    % (to use reduced model)
%     activity_predicted(~isnan(activity')) = activity_predicted_reduced_reshape(:,:,1);
    activity_predicted = activity_predicted';    
    fluor_allcat_predicted(regress_trials,:,curr_roi,1) = activity_predicted;    
    
    fluor_kernel(:,curr_roi) = cellfun(@(x,t) reshape(x,[],length(t)), ...
        mat2cell(kernel_flat,[cellfun(@(x) size(x,1),regressors_reshape);1].*[cellfun(@length,sample_shifts);1],1), ...
        [sample_shifts;{1}],'uni',false);
    
    AP_print_progress_fraction(curr_roi,n_rois);
    
end

fluor_allcat_residual = fluor_allcat_downsamp_diff - fluor_allcat_predicted;

use_svs = 200;
% Plot regressors
for curr_regressor = 1:length(regressors)
    curr_k_v = permute(cell2mat(permute(fluor_kernel(curr_regressor,1:use_svs),[1,3,2])),[3,2,1]);
    curr_k_px = reshape(svdFrameReconstruct(U_master(:,:,1:use_svs), ...
        reshape(curr_k_v,use_svs,[])),size(U_master,1),size(U_master,2),[],size(curr_k_v,3));
    AP_imscroll(curr_k_px,sample_shifts{curr_regressor}/(sample_rate/downsample_factor));
    axis image
    caxis([-prctile(abs(curr_k_px(:)),100),prctile(abs(curr_k_px(:)),100)]);
    colormap(brewermap([],'*RdBu'))
    AP_reference_outline('ccf_aligned','k');
    set(gcf,'Name',regressor_labels{curr_regressor});
end
% Plot constant
curr_k_v = permute(cell2mat(permute(fluor_kernel(length(regressors)+1,1:use_svs),[1,3,2])),[3,2,1]);
curr_k_px = reshape(svdFrameReconstruct(U_master(:,:,1:use_svs), ...
    reshape(curr_k_v,use_svs,[])),size(U_master,1),size(U_master,2),[],size(curr_k_v,3));
figure;imagesc(curr_k_px);
axis image off
caxis([-prctile(abs(curr_k_px(:)),100),prctile(abs(curr_k_px(:)),100)]);
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k');
title('Constant');

% Settings to plot
% rxn_time_bins = {[0,0.1],[0.1,0.2],[0.2,0.3],[0.3,0.4]};
% rxn_time_bins = {[0,0.15],[0.15,0.4],[0.6,0.7]};
% rxn_time_bins = {[0,0.5],[0.5,1]};
rxn_time_bins = {[0.1,0.3],[0.6,0.7]};
% rxn_time_bins = {[0.1,0.3]};
% rxn_time_bins = {[0.6,0.7]};

normalize_px = true;

% Get major trial types
vis_L_trials_hit = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == -1;

vis_R_trials_hit = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == -1 & ...
    trial_choice_allcat == 1;

vis_L_trials_miss = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == -1 & ...
    trial_choice_allcat == -1;

vis_R_trials_miss = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == 1;

zero_L_trials = ...
    trial_contrast_allcat == 0 & ...
    trial_choice_allcat == -1;

zero_R_trials = ...
    trial_contrast_allcat == 0 & ...
    trial_choice_allcat == 1;

trial_types = ...
    [vis_L_trials_hit, vis_R_trials_hit, ...
    vis_L_trials_miss, vis_R_trials_miss, ...
    zero_L_trials, zero_R_trials];

%%% Get and plot fluorescence
px_trial_types = nan(size(U_master,1),size(U_master,2), ...
    length(t_downsample_diff),size(trial_types,2),length(rxn_time_bins));
for curr_rxn = 1:length(rxn_time_bins)
    for curr_trial_type = 1:size(trial_types,2)
        
        curr_trials = trial_types(:,curr_trial_type) & ...
            move_t > rxn_time_bins{curr_rxn}(1) & ...
            move_t < rxn_time_bins{curr_rxn}(2);
        
        curr_data = fluor_allcat_downsamp_diff(curr_trials,:,1:use_svs);

        % re-align to movement onset
        t_leeway = 0.5;
        leeway_samples = round(t_leeway*(sample_rate/downsample_factor));
        curr_move_idx = move_idx(curr_trials);
        for i = 1:size(curr_data,1)
            curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
        end
        
        curr_data_mean = squeeze(nanmean(curr_data,1))';
        
        curr_px = svdFrameReconstruct(U_master(:,:,1:use_svs),curr_data_mean);
        
        px_trial_types(:,:,:,curr_trial_type,curr_rxn) = curr_px;
        
    end
    AP_print_progress_fraction(curr_rxn,length(rxn_time_bins));
end

% Flip and combine trial types
% 1) visual hit (contra), 2) visual miss (ipsi), 3) zero (L)
px_combined = cat(4, ...
    (px_trial_types(:,:,:,1,:) + AP_reflect_widefield(px_trial_types(:,:,:,2,:)))./2, ...
    (px_trial_types(:,:,:,3,:) + AP_reflect_widefield(px_trial_types(:,:,:,4,:)))./2, ...
    (px_trial_types(:,:,:,5,:) + AP_reflect_widefield(px_trial_types(:,:,:,6,:)))./2);

px_combined_hemidiff = px_combined - AP_reflect_widefield(px_combined);

if normalize_px
    % Normalize by dividing by max of each frame
    px_dims = size(px_combined);
    px_combined = bsxfun(@rdivide,px_combined,max(abs(reshape( ...
         px_combined,[px_dims(1)*px_dims(2),1,px_dims(3:end)])),[],1));
     
    px_combined_hemidiff = bsxfun(@rdivide,px_combined_hemidiff,max(abs(reshape( ...
         px_combined_hemidiff,[px_dims(1)*px_dims(2),1,px_dims(3:end)])),[],1));
end

% Plot
plot_rxn = 1;
AP_imscroll(cat(4,px_combined(:,:,:,:,plot_rxn),px_combined_hemidiff(:,:,:,:,plot_rxn)),t_downsample_diff);
axis image; caxis([-1,1]); 
% colormap(colormap_BlueWhiteRed);
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k');

% Plot all concatenated
px_dims = size(px_combined);
AP_imscroll([reshape(permute(px_combined,[1,4,2,5,3]), ...
    [],size(px_combined,2)*size(px_combined,5),size(px_combined,3)), ...
    reshape(permute(px_combined_hemidiff,[1,4,2,5,3]), ...
    [],size(px_combined_hemidiff,2)*size(px_combined_hemidiff,5), ...
    size(px_combined_hemidiff,3))],t_downsample_diff);
axis image;
caxis([-1,1]);
colormap(brewermap([],'*RdBu'))


%%% Do regression on MUA

mua_allcat_predicted = nan(size(mua_allcat_downsamp));

mua_kernel = cell(length(regressors)+1,n_depths);
mua_expl_var = cell(length(regressors),n_depths);
for curr_depth = 1:n_depths
    
    activity = mua_allcat_downsamp(regress_trials,:,curr_depth,1);
    
    activity_reshape = reshape(activity',[],1)';
    regressors_reshape = cellfun(@(x) ...
        reshape(permute(x(regress_trials,:,:),[2,1,3]),[],size(x,3))',regressors,'uni',false);
    
    [kernel_flat,activity_predicted_reshape,expl_var,activity_predicted_reduced_reshape] = AP_regresskernel( ...
        cellfun(@(x) x(:,~isnan(activity_reshape)),regressors_reshape,'uni',false), ...
        activity_reshape(~isnan(activity_reshape)),sample_shifts,lambda,zs,cvfold,return_constant);
    
    activity_predicted = nan(size(activity))';
    % (to use full model)
    activity_predicted(~isnan(activity')) = activity_predicted_reshape;
    % (to use reduced model)
%     activity_predicted(~isnan(activity')) = activity_predicted_reduced_reshape(:,:,1);
    activity_predicted = activity_predicted';    
    mua_allcat_predicted(regress_trials,:,curr_depth,1) = activity_predicted;
    
    mua_kernel(:,curr_depth) = cellfun(@(x,t) reshape(x,[],length(t)), ...
        mat2cell(kernel_flat,[cellfun(@(x) size(x,1),regressors_reshape);1].*[cellfun(@length,sample_shifts);1],1), ...
        [sample_shifts;{1}],'uni',false);
    
    AP_print_progress_fraction(curr_depth,n_depths);
    
end

mua_allcat_residual = mua_allcat_downsamp - mua_allcat_predicted;

figure;
for curr_regressor = 1:length(regressors)
    for curr_depth = 1:n_depths
        subplot(length(regressors),n_depths,(curr_regressor-1)*n_depths+curr_depth)
        imagesc(sample_shifts{curr_regressor}/(sample_rate/downsample_factor),[], ...
            mua_kernel{curr_regressor,curr_depth})
        colormap(hot);
        title(regressor_labels{curr_regressor});
    end
end


mua_trial_types = nan(n_depths, ...
    length(t_downsample_diff),size(trial_types,2),length(rxn_time_bins));
for curr_rxn = 1:length(rxn_time_bins)
    for curr_trial_type = 1:size(trial_types,2)
        
        curr_trials = trial_types(:,curr_trial_type) & ...
            move_t > rxn_time_bins{curr_rxn}(1) & ...
            move_t < rxn_time_bins{curr_rxn}(2);
              
        curr_data = mua_allcat_residual(curr_trials,:,:);

        % re-align to movement onset
        t_leeway = 0.5;
        leeway_samples = round(t_leeway*(sample_rate/downsample_factor));
        curr_move_idx = move_idx(curr_trials);
        for i = 1:size(curr_data,1)
            curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
        end        
                
        mua_trial_types(:,:,curr_trial_type,curr_rxn) = squeeze(nanmean(curr_data,1))';
               
    end
    AP_print_progress_fraction(curr_rxn,length(rxn_time_bins));
end

for curr_rxn = 1:length(rxn_time_bins)
    figure;
    for curr_trial_type = 1:size(trial_types,2)
        subplot(size(trial_types,2)/2,2,curr_trial_type); hold on;
        set(gca,'ColorOrder',copper(4));
        plot(t_downsample_diff,mua_trial_types(:,:,curr_trial_type,curr_rxn)','linewidth',2);
    end
end

% Save regression results
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
save_fn = ['all_trial_activity_regressed'];
save([save_path filesep save_fn], ...
    'trial_contrast_allcat','trial_side_allcat','trial_choice_allcat', ...
    'downsample_factor','t_downsample_diff','wheel', ...
    'fluor_allcat_downsamp_diff','mua_allcat_downsamp', ...
    'regressors','regressor_labels','t_shifts','regress_trials', ...
    'fluor_allcat_predicted','fluor_kernel', ...
    'mua_allcat_predicted','mua_kernel','-v7.3');
disp('Saved regression');

%% Regression >> load/prepare regression results from above

% Load regression results
regression_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\all_trial_activity_regressed.mat';
load(regression_fn);

% Get number of V's/depths
n_vs = size(fluor_allcat_predicted,3);
n_depths = size(mua_allcat_predicted,3);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Get move onset index
[~,move_idx] = max(abs(wheel(:,:,1)) > 2,[],2);
move_t = t_downsample_diff(move_idx)';

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Set which components to use
use_svs = 200;

% Get time shifts in samples
sample_shifts = cellfun(@(x) round(x(1)*(sample_rate/downsample_factor)): ...
    round(x(2)*(sample_rate/downsample_factor)),t_shifts,'uni',false);

disp('Regression results loaded');

%% Regression >> spatial explained variance

use_t = t_downsample_diff > 0 & t_downsample_diff < 0.5;
use_trials = move_t > 0 & move_t < 0.5;

spatial_explained_var = AP_spatial_explained_var(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_downsamp_diff(use_trials,use_t,:),[2,1,3]),[],n_vs)', ...
    reshape(permute(fluor_allcat_predicted(use_trials,use_t,:),[2,1,3]),[],n_vs)',10);
figure;imagesc(spatial_explained_var);
axis image off; 
colormap(brewermap([],'Reds'));
caxis([0,1]); colorbar;
AP_reference_outline('ccf_aligned','k');

%% Regression >> plot kernels

% Plot fluorescence kernels
for curr_regressor = 1:length(regressors)
    curr_k_v = permute(cell2mat(permute(fluor_kernel(curr_regressor,1:use_svs),[1,3,2])),[3,2,1]);
    curr_k_px = reshape(svdFrameReconstruct(U_master(:,:,1:use_svs), ...
        reshape(curr_k_v,use_svs,[])),size(U_master,1),size(U_master,2),[],size(curr_k_v,3));
    AP_imscroll(curr_k_px,sample_shifts{curr_regressor}/(sample_rate/downsample_factor));
    axis image
    caxis([-prctile(abs(curr_k_px(:)),100),prctile(abs(curr_k_px(:)),100)]);
    colormap(brewermap([],'*RdBu'))
    AP_reference_outline('ccf_aligned','k');
    set(gcf,'Name',regressor_labels{curr_regressor});
end
% Plot fluorescence constant
curr_k_v = permute(cell2mat(permute(fluor_kernel(length(regressors)+1,1:use_svs),[1,3,2])),[3,2,1]);
curr_k_px = reshape(svdFrameReconstruct(U_master(:,:,1:use_svs), ...
    reshape(curr_k_v,use_svs,[])),size(U_master,1),size(U_master,2),[],size(curr_k_v,3));
figure;imagesc(curr_k_px);
axis image off
caxis([-prctile(abs(curr_k_px(:)),100),prctile(abs(curr_k_px(:)),100)]);
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k');
title('Constant');

% Plot fluorescence ROI kernels
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);
n_rois = numel(wf_roi);

roi_mask = cat(3,wf_roi.mask);

U_roi = bsxfun(@rdivide,transpose(reshape(U_master,[],size(U_master,3))'* ...
    reshape(roi_mask,[],size(roi_mask,3))), ...
    sum(reshape(roi_mask,[],size(roi_mask,3)),1)');

for curr_regressor = 1:length(regressors)
    curr_k_v = permute(cell2mat(permute(fluor_kernel(curr_regressor,1:use_svs),[1,3,2])),[3,2,1]);
    
    curr_k_roi = reshape(AP_svd_roi(U_master(:,:,1:n_vs), ...
        reshape(curr_k_v,n_vs,[]),[],[],roi_mask)', ...
        size(curr_k_v,2),size(curr_k_v,3),n_rois);
    
    if size(curr_k_roi,2) > 1
        curr_col = colormap_BlueWhiteRed(size(curr_k_roi,2)/2);
        curr_col(size(curr_k_roi,2)/2+1,:) = [];
    else
        curr_col = 'k'
    end
    
    figure; hold on;
    for curr_subk = 1:size(curr_k_roi,2)
        AP_stackplot(squeeze(curr_k_roi(:,curr_subk,:)), ...
            sample_shifts{curr_regressor}/(sample_rate/downsample_factor), ...
            1.5e-3,false,curr_col(curr_subk,:));
    end
    title(regressor_labels{curr_regressor});
end

% Plot MUA kernels
figure;
for curr_regressor = 1:length(regressors)
    
    curr_k_cat = permute(cat(3,mua_kernel{curr_regressor,:}),[2,1,3]);
    
    if size(curr_k_cat,2) > 1
        curr_col = colormap_BlueWhiteRed(size(curr_k_cat,2)/2);
        curr_col(size(curr_k_cat,2)/2+1,:) = [];
    else
        curr_col = 'k'
    end
    
    figure; hold on;
    for curr_subk = 1:size(curr_k_cat,2)
        AP_stackplot(squeeze(curr_k_cat(:,curr_subk,:)), ...
            sample_shifts{curr_regressor}/(sample_rate/downsample_factor), ...
            1,false,curr_col(curr_subk,:));
    end
    title(regressor_labels{curr_regressor});
end


%% Regression >> line plot results

% Get fluorescence in widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);
n_rois = numel(wf_roi);
roi_mask = cat(3,wf_roi.mask);

U_roi = bsxfun(@rdivide,transpose(reshape(U_master,[],size(U_master,3))'* ...
    reshape(roi_mask,[],size(roi_mask,3))), ...
    sum(reshape(roi_mask,[],size(roi_mask,3)),1)');

fluor_downsamp_roi = permute(reshape(transpose(U_roi(:,1:n_vs)* ...
    reshape(permute(fluor_allcat_downsamp_diff,[3,2,1,4]),n_vs,[])), ...
    size(fluor_allcat_downsamp_diff,2),size(fluor_allcat_downsamp_diff,1),n_rois),[2,1,3]);

fluor_predicted_roi = permute(reshape(transpose(U_roi(:,1:n_vs)* ...
    reshape(permute(fluor_allcat_predicted,[3,2,1,4]),n_vs,[])), ...
    size(fluor_allcat_downsamp_diff,2),size(fluor_allcat_downsamp_diff,1),n_rois),[2,1,3]);

% Get residuals
fluor_residual_roi = fluor_downsamp_roi - fluor_predicted_roi;
mua_allcat_residual = mua_allcat_downsamp - mua_allcat_predicted;

% Low-pass filter fluorescence (where does this ~10 Hz crap come from?)
lowpassCutoff = 6; % Hz
[b100s, a100s] = butter(2, lowpassCutoff/((sample_rate/downsample_factor)/2), 'low');
fluor_downsamp_roi = filter(b100s,a100s,fluor_downsamp_roi,[],2);
fluor_predicted_roi = filter(b100s,a100s,fluor_predicted_roi,[],2);
fluor_residual_roi = fluor_downsamp_roi - fluor_predicted_roi;

% Plot mean residuals
fluor_downsamp_roi_move = fluor_downsamp_roi;
fluor_residual_roi_move = fluor_residual_roi;
% re-align to movement onset
t_leeway = 0.5;
leeway_samples = round(t_leeway*(sample_rate/downsample_factor));
for i = 1:size(fluor_downsamp_roi_move,1)
    fluor_downsamp_roi_move(i,:,:) = circshift(fluor_downsamp_roi_move(i,:,:),-move_idx(i)+leeway_samples,2);
    fluor_residual_roi_move(i,:,:) = circshift(fluor_residual_roi_move(i,:,:),-move_idx(i)+leeway_samples,2);
end

mua_allcat_downsamp_move = mua_allcat_downsamp;
mua_allcat_residual_move = mua_allcat_residual;
% re-align to movement onset
t_leeway = 0.5;
leeway_samples = round(t_leeway*(sample_rate/downsample_factor));
for i = 1:size(mua_allcat_downsamp_move,1)
    mua_allcat_downsamp_move(i,:,:) = circshift(mua_allcat_downsamp_move(i,:,:),-move_idx(i)+leeway_samples,2);
    mua_allcat_residual_move(i,:,:) = circshift(mua_allcat_residual_move(i,:,:),-move_idx(i)+leeway_samples,2);
end

plot_residual_trials = regress_trials & move_t > 0.1 & move_t < 0.3;

fluor_downsamp_mean = squeeze(nanmean(fluor_downsamp_roi_move(plot_residual_trials,:,:),1));
fluor_residual_mean = squeeze(nanmean(fluor_residual_roi_move(plot_residual_trials,:,:),1));
fluor_residual_abs_mean = squeeze(nanmean(abs(fluor_residual_roi_move(plot_residual_trials,:,:)),1));
% fluor_residual_abs_mean = bsxfun(@minus,fluor_residual_abs_mean, ...
%     nanmean(fluor_residual_abs_mean(t_downsample_diff < -0.2,:),1));

mua_downsamp_mean = squeeze(nanmean(mua_allcat_downsamp_move(plot_residual_trials,:,:),1));
mua_residual_mean = squeeze(nanmean(mua_allcat_residual_move(plot_residual_trials,:,:),1));
mua_residual_abs_mean = squeeze(nanmean(abs(mua_allcat_residual_move(plot_residual_trials,:,:)),1));
% mua_residual_abs_mean = bsxfun(@minus,mua_residual_abs_mean, ...
%     nanmean(mua_residual_abs_mean(t_downsample_diff < -0.2,:),1));

fluor_sse_measured = squeeze(nansum(fluor_downsamp_roi_move(plot_residual_trials,:,:).^2,1));
fluor_sse_residual = squeeze(nansum(fluor_residual_roi_move(plot_residual_trials,:,:).^2,1));
fluor_expl_var = (fluor_sse_measured - fluor_sse_residual)./fluor_sse_measured;

mua_sse_measured = squeeze(nansum(mua_allcat_downsamp_move(plot_residual_trials,:,:).^2,1));
mua_sse_residual = squeeze(nansum(mua_allcat_residual_move(plot_residual_trials,:,:).^2,1));
mua_expl_var = (mua_sse_measured - mua_sse_residual)./mua_sse_measured;

figure;
subplot(4,2,1); hold on;
set(gca,'ColorOrder',jet(size(wf_roi,1)));
plot(t_downsample_diff,fluor_downsamp_mean,'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Mean activity')
xlabel('Time from move');
legend({wf_roi.area})

subplot(4,2,2); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot(t_downsample_diff,mua_downsamp_mean,'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Mean activity')
xlabel('Time from move');

subplot(4,2,3); hold on;
set(gca,'ColorOrder',jet(size(wf_roi,1)));
plot(t_downsample_diff,fluor_expl_var,'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Variance explained')
xlabel('Time from move');

subplot(4,2,4); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot(t_downsample_diff,mua_expl_var,'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Variance explained')
xlabel('Time from move');

subplot(4,2,5); hold on;
set(gca,'ColorOrder',jet(size(wf_roi,1)));
plot(t_downsample_diff,fluor_residual_mean,'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Mean residual')
xlabel('Time from move');

subplot(4,2,6); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot(t_downsample_diff,mua_residual_mean,'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Mean residual')
xlabel('Time from move');

subplot(4,2,7); hold on;
set(gca,'ColorOrder',jet(size(wf_roi,1)));
plot(t_downsample_diff,fluor_residual_abs_mean,'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Mean abs residual')
xlabel('Time from move');

subplot(4,2,8); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot(t_downsample_diff,mua_residual_abs_mean,'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Mean abs residual')
xlabel('Time from move');

% Plot total explained variance
figure;
subplot(2,1,1);
plot((sum(fluor_sse_measured,1)-sum(fluor_sse_residual,1))./ ...
    sum(fluor_sse_measured,1),'k','linewidth',2)
set(gca,'XTick',1:length(wf_roi),'XTickLabel',{wf_roi.area});
ylabel('Explained variance');

subplot(2,1,2);
plot((sum(mua_sse_measured,1)-sum(mua_sse_residual,1))./ ...
    sum(mua_sse_measured,1),'k','linewidth',2)
ylabel('Explained variance');

% Plot all trials measured and predicted
[~,sort_idx] = sort(move_idx);
AP_imscroll([fluor_downsamp_roi(sort_idx,:,:),fluor_predicted_roi(sort_idx,:,:),fluor_residual_roi(sort_idx,:,:)],{wf_roi.area})
AP_imscroll([mua_allcat_downsamp(sort_idx,:,:),mua_allcat_predicted(sort_idx,:,:),mua_allcat_residual(sort_idx,:,:)])

% Get conditions for each trial, plot selected
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];

contrast_side_col = colormap_BlueWhiteRed(5);
contrast_side_col(6,:) = 0;
contrast_side_val = unique(sort([-contrasts,contrasts]))';

% contrast, side, choice
% plot_conditions = ...
%     [contrasts,contrasts; ...
%     -ones(1,6),-1,ones(1,5); ...
%     ones(1,6),-ones(1,6)]';
% plot_conditions = ...
%     [contrasts(2:end),contrasts(2:end); ...
%     -ones(1,5),ones(1,5); ...
%     ones(1,5),-ones(1,5)]';
% plot_conditions = ...
%     [0,0; ...
%     -1,-1; ...
%     -1,1]';
plot_conditions = ...
    [1,0,1; ...
    1,-1,-1; ...
    -1,-1,-1]';

use_rxn = move_t > 0  & move_t < 0.5;
    
[~,plot_id] = ismember( ...
    [trial_contrast_allcat > 0,trial_side_allcat,trial_choice_allcat], ...
    plot_conditions,'rows');

% Plot cortex
figure;
for curr_plot = 1:3
    
    switch curr_plot
        case 1
            plot_data = fluor_downsamp_roi;
            plot_title = 'Measured';
        case 2
            plot_data = fluor_predicted_roi;
            plot_title = 'Predicted';
        case 3
            plot_data = fluor_residual_roi;
            plot_title = 'Residual';
    end
    
    p_ctx(curr_plot) = subplot(1,3,curr_plot); hold on;
    for curr_plot_condition = 1:size(plot_conditions,1)
        
        curr_trials = plot_id == curr_plot_condition & use_rxn;
        curr_data = plot_data(curr_trials,:,:);
        
        % re-align to movement onset
        t_leeway = 0.5;
        leeway_samples = round(t_leeway*(sample_rate/downsample_factor));
        curr_move_idx = move_idx(curr_trials);
        for i = 1:size(curr_data,1)
            curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
        end
        
        curr_data_mean = squeeze(nanmean(curr_data,1));
        
        curr_contrast_side = max(trial_contrast_allcat(curr_trials))*sign(max(trial_side_allcat(curr_trials)));
        contrast_side_idx = find(curr_contrast_side == contrast_side_val);
        curr_col = contrast_side_col(contrast_side_idx,:);
        
        if curr_contrast_side == 0
            switch max(trial_choice_allcat(curr_trials))
                case -1
                    curr_col = 'k';
                case 1
                    curr_col = 'c';
            end
        end
        
        AP_stackplot(curr_data_mean,t_downsample_diff,2e-3,false,curr_col,{wf_roi.area});
        
    end
    line([0,0],ylim,'color','k');
    xlabel('Time from stim');
    title(plot_title);   
end
linkaxes(p_ctx)

% Plot striatum
figure; hold on;
for curr_plot = 1:3
    
    switch curr_plot
        case 1
            plot_data = mua_allcat_downsamp;
            plot_title = 'Measured';
        case 2
            plot_data = mua_allcat_predicted;
            plot_title = 'Predicted';
        case 3
            plot_data = mua_allcat_residual;
            plot_title = 'Residual';
    end    
    
    p_str(curr_plot) = subplot(1,3,curr_plot); hold on;
    for curr_plot_condition = 1:size(plot_conditions,1)
        
        curr_trials = plot_id == curr_plot_condition & use_rxn;
        curr_data = plot_data(curr_trials,:,:);
        
        % re-align to movement onset
        t_leeway = 0.5;
        leeway_samples = round(t_leeway*(sample_rate/downsample_factor));
        curr_move_idx = move_idx(curr_trials);
        for i = 1:size(curr_data,1)
            curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
        end
        
        curr_data_mean = squeeze(nanmean(curr_data,1));
        
        curr_contrast_side = max(trial_contrast_allcat(curr_trials))*sign(max(trial_side_allcat(curr_trials)));
        contrast_side_idx = find(curr_contrast_side == contrast_side_val);
        curr_col = contrast_side_col(contrast_side_idx,:);
        
        if curr_contrast_side == 0
            switch max(trial_choice_allcat(curr_trials))
                case -1
                    curr_col = 'k';
                case 1
                    curr_col = 'c';
            end
        end
        
        AP_stackplot(curr_data_mean,t_downsample_diff,2,false,curr_col,1:n_depths);
        
    end
    line([0,0],ylim,'color','k');
    xlabel('Time from stim');
    title(plot_title);
end
linkaxes(p_str)

%% Regression >> movie results

t_shifts = {[-0.1,0.3];[-0.2,0.5];[0,0];[-0.05,0.3];[-0.1,0.5]};

% Settings to plot (only use one bin here)
rxn_time_bins = {[0.3,0.4]};

normalize_px = true;
move_align = true;

% Get major trial types
vis_L_trials_hit = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == -1;

vis_R_trials_hit = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == -1 & ...
    trial_choice_allcat == 1;

vis_L_trials_miss = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == -1 & ...
    trial_choice_allcat == -1;

vis_R_trials_miss = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == 1;

zero_L_trials = ...
    trial_contrast_allcat == 0 & ...
    trial_choice_allcat == -1;

zero_R_trials = ...
    trial_contrast_allcat == 0 & ...
    trial_choice_allcat == 1;

trial_types = ...
    [vis_L_trials_hit, vis_R_trials_hit, ...
    vis_L_trials_miss, vis_R_trials_miss, ...
    zero_L_trials, zero_R_trials];

%%% Get and plot fluorescence
px_trial_types = nan(size(U_master,1),size(U_master,2), ...
    length(t_downsample_diff),size(trial_types,2),length(rxn_time_bins));
px_predicted_trial_types = nan(size(U_master,1),size(U_master,2), ...
    length(t_downsample_diff),size(trial_types,2),length(rxn_time_bins));
for curr_rxn = 1:length(rxn_time_bins)
    for curr_trial_type = 1:size(trial_types,2)
        
        curr_trials = trial_types(:,curr_trial_type) & ...
            move_t > rxn_time_bins{curr_rxn}(1) & ...
            move_t < rxn_time_bins{curr_rxn}(2);
        
        curr_data = fluor_allcat_downsamp_diff(curr_trials,:,1:use_svs);
        curr_data_predicted = fluor_allcat_predicted(curr_trials,:,1:use_svs);
        
        if move_align
            % re-align to movement onset
            t_leeway = 0.5;
            leeway_samples = round(t_leeway*(sample_rate/downsample_factor));
            curr_move_idx = move_idx(curr_trials);
            for i = 1:size(curr_data,1)
                curr_data(i,:) = circshift(curr_data(i,:),-curr_move_idx(i)+leeway_samples,2);
                curr_data_predicted(i,:) = circshift(curr_data_predicted(i,:),-curr_move_idx(i)+leeway_samples,2);
            end
        end
        
        curr_data_mean = squeeze(nanmean(curr_data,1))';
        curr_data_predicted_mean = squeeze(nanmean(curr_data_predicted,1))';
        
        curr_px = svdFrameReconstruct(U_master(:,:,1:use_svs),curr_data_mean);
        curr_px_predicted = svdFrameReconstruct(U_master(:,:,1:use_svs),curr_data_predicted_mean);
        
        px_trial_types(:,:,:,curr_trial_type,curr_rxn) = curr_px;
        px_predicted_trial_types(:,:,:,curr_trial_type,curr_rxn) = curr_px_predicted;
        
    end
    AP_print_progress_fraction(curr_rxn,length(rxn_time_bins));
end

% Flip and combine trial types
% 1) visual hit (contra), 2) visual miss (ipsi), 3) zero (L)
px_combined = cat(4, ...
    (px_trial_types(:,:,:,1,:) + AP_reflect_widefield(px_trial_types(:,:,:,2,:)))./2, ...
    (px_trial_types(:,:,:,3,:) + AP_reflect_widefield(px_trial_types(:,:,:,4,:)))./2, ...
    (px_trial_types(:,:,:,5,:) + AP_reflect_widefield(px_trial_types(:,:,:,6,:)))./2);
px_combined_hemidiff = px_combined - AP_reflect_widefield(px_combined);

px_predicted_combined = cat(4, ...
    (px_predicted_trial_types(:,:,:,1,:) + AP_reflect_widefield(px_predicted_trial_types(:,:,:,2,:)))./2, ...
    (px_predicted_trial_types(:,:,:,3,:) + AP_reflect_widefield(px_predicted_trial_types(:,:,:,4,:)))./2, ...
    (px_predicted_trial_types(:,:,:,5,:) + AP_reflect_widefield(px_predicted_trial_types(:,:,:,6,:)))./2);
px_predicted_combined_hemidiff = px_predicted_combined - AP_reflect_widefield(px_predicted_combined);

if normalize_px
    % Normalize by dividing by max of each unpredicted frame
    px_dims = size(px_combined);
    px_norm = max(abs(reshape(px_combined, ...
        [px_dims(1)*px_dims(2),1,px_dims(3:end)])),[],1);
    px_hemidiff_norm = max(abs(reshape(px_combined_hemidiff, ...
        [px_dims(1)*px_dims(2),1,px_dims(3:end)])),[],1);
    
    px_combined = bsxfun(@rdivide,px_combined,px_norm);     
    px_combined_hemidiff = bsxfun(@rdivide,px_combined_hemidiff,px_hemidiff_norm);
    
    px_predicted_combined = bsxfun(@rdivide,px_predicted_combined,px_norm);     
    px_predicted_combined_hemidiff = bsxfun(@rdivide,px_predicted_combined_hemidiff,px_hemidiff_norm);
end

% Plot all concatenated
px_residual = px_combined - px_predicted_combined;
px_dims = size(px_combined);
AP_imscroll([reshape(permute(px_combined,[1,4,2,5,3]), ...
    [],size(px_combined,2)*size(px_combined,5),size(px_combined,3)), ...
    reshape(permute(px_predicted_combined,[1,4,2,5,3]), ...
    [],size(px_predicted_combined,2)*size(px_predicted_combined,5),size(px_predicted_combined,3)), ...
    reshape(permute(px_residual,[1,4,2,5,3]), ...
    [],size(px_residual,2)*size(px_residual,5),size(px_residual,3))],t_downsample_diff);
AP_reference_outline('ccf_aligned','k',[], ...
    [size(U_master,1),size(U_master,2),3,length(rxn_time_bins)*3]);
axis image;
caxis([-1,1]);
colormap(brewermap([],'*RdBu'))


%% Regression >> inter-area correlations by prediction/downsampled

% Get fluorescence in widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);
n_rois = numel(wf_roi);
roi_mask = cat(3,wf_roi.mask);

U_roi = bsxfun(@rdivide,transpose(reshape(U_master,[],size(U_master,3))'* ...
    reshape(roi_mask,[],size(roi_mask,3))), ...
    sum(reshape(roi_mask,[],size(roi_mask,3)),1)');

fluor_downsamp_roi = permute(reshape(transpose(U_roi(:,1:n_vs)* ...
    reshape(permute(fluor_allcat_downsamp_diff,[3,2,1,4]),n_vs,[])), ...
    size(fluor_allcat_downsamp_diff,2),size(fluor_allcat_downsamp_diff,1),n_rois),[2,1,3]);

fluor_predicted_roi = permute(reshape(transpose(U_roi(:,1:n_vs)* ...
    reshape(permute(fluor_allcat_predicted,[3,2,1,4]),n_vs,[])), ...
    size(fluor_allcat_downsamp_diff,2),size(fluor_allcat_downsamp_diff,1),n_rois),[2,1,3]);

% Get residuals
fluor_residual_roi = fluor_downsamp_roi - fluor_predicted_roi;
mua_allcat_residual = mua_allcat_downsamp - mua_allcat_predicted;

% Low-pass filter fluorescence (where does this ~10 Hz crap come from?)
lowpassCutoff = 6; % Hz
[b100s, a100s] = butter(2, lowpassCutoff/((sample_rate/downsample_factor)/2), 'low');

fluor_downsamp_roi_filt = filter(b100s,a100s,fluor_downsamp_roi,[],2);
fluor_predicted_roi_filt = filter(b100s,a100s,fluor_predicted_roi,[],2);
fluor_residual_roi_filt = fluor_downsamp_roi_filt - fluor_predicted_roi_filt;

mua_allcat_downsamp_filt = filter(b100s,a100s,mua_allcat_downsamp,[],2);
mua_allcat_predicted_filt = filter(b100s,a100s,mua_allcat_predicted,[],2);
mua_allcat_residual_filt = mua_allcat_downsamp_filt - mua_allcat_predicted_filt;

% re-align to movement onset
fluor_downsamp_roi_filt_move = fluor_downsamp_roi_filt;
fluor_predicted_roi_filt_move = fluor_predicted_roi_filt;
mua_allcat_downsamp_filt_move = mua_allcat_downsamp_filt;
mua_allcat_predicted_filt_move = mua_allcat_predicted_filt;
t_leeway = 0.5;
leeway_samples = round(t_leeway*(sample_rate/downsample_factor));
for i = 1:size(fluor_downsamp_roi_filt,1)
    fluor_downsamp_roi_filt_move(i,:,:) = circshift(fluor_downsamp_roi_filt_move(i,:,:),-move_idx(i)+leeway_samples,2);
    fluor_predicted_roi_filt_move(i,:,:) = circshift(fluor_predicted_roi_filt_move(i,:,:),-move_idx(i)+leeway_samples,2);
    mua_allcat_downsamp_filt_move(i,:,:) = circshift(mua_allcat_downsamp_filt_move(i,:,:),-move_idx(i)+leeway_samples,2);    
    mua_allcat_predicted_filt_move(i,:,:) = circshift(mua_allcat_predicted_filt_move(i,:,:),-move_idx(i)+leeway_samples,2);
end
fluor_residual_roi_filt_move = fluor_downsamp_roi_filt_move - fluor_predicted_roi_filt_move;
mua_allcat_residual_filt_move = mua_allcat_downsamp_filt_move - mua_allcat_predicted_filt_move;

% Concatenate all activity
%%% Real
% activity_cat = cat(3,fluor_downsamp_roi_filt_move,mua_allcat_downsamp_filt_move);
%%% Predicted
% activity_cat = cat(3,fluor_predicted_roi_filt_move,mua_allcat_predicted_filt_move);
%%% Predicted with independent noise
% activity_cat = cat(3,fluor_predicted_roi_filt + randn(size(fluor_predicted_roi))*0.1*nanstd(fluor_predicted_roi(:)), ...
%     mua_allcat_predicted_filt + randn(size(mua_allcat_predicted))*1.5*nanstd(mua_allcat_predicted(:)));
%%% Predicted with common trial noise
fluor_trial_noise = randn(size(fluor_predicted_roi(:,1,:)))*0.1*nanstd(fluor_predicted_roi(:));
mua_trial_noise = randn(size(mua_allcat_predicted(:,1,:)))*1.5*nanstd(mua_allcat_predicted(:));
activity_cat = cat(3, ...
    bsxfun(@plus,fluor_predicted_roi_filt_move,fluor_trial_noise), ...
    bsxfun(@plus,mua_allcat_predicted_filt_move,mua_trial_noise));

area_labels = [{wf_roi(:,1).area}, ...
    cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false)];

% Correlate all pairs of areas
use_trials = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == ...
    -trial_choice_allcat & ...
    move_t > 0 & move_t < 0.5;
% use_trials = ...
%     trial_contrast_allcat == 0 & ...
%     move_t > 0 & move_t < 0.5;
% use_trials = true(size(move_t));
% use_trials = ...
%     trial_contrast_allcat > 0;
% use_trials = ...
%     move_t > 0.5 & move_t < 1;

corr_grid = cell(size(activity_cat,3));
for curr_area_1 = 1:size(activity_cat,3)
    for curr_area_2 = 1:curr_area_1-1
         
        nonan_trials = ~any(any(isnan(activity_cat(:,:,curr_area_1,:)),2),4) & ...
            ~any(any(isnan(activity_cat(:,:,curr_area_2,:)),2),4);
        
        compare_act_corr = 1-pdist2(activity_cat(nonan_trials & use_trials,:,curr_area_1)', ...
            activity_cat(nonan_trials & use_trials,:,curr_area_2)','correlation');
        
        corr_grid{curr_area_1,curr_area_2} = compare_act_corr;                
        
    end
    
    AP_print_progress_fraction(curr_area_1,size(activity_cat,3));
    
end

area_labels_grid = cellfun(@(x,y) [x '-' y], ...
    repmat(area_labels',1,length(area_labels)), ...
    repmat(area_labels,length(area_labels),1), ...
    'uni',false);
area_labels_grid_used = area_labels_grid(~cellfun(@isempty,corr_grid));

corr_grid_cat = cat(3,corr_grid{:});

AP_imscroll(corr_grid_cat,area_labels_grid_used);
axis image;
colormap(brewermap([],'*RdBu'));
caxis([-0.5,0.5])
line(repmat(find(t_downsample_diff > 0,1),2,1),ylim,'color','k');
line(xlim,repmat(find(t_downsample_diff > 0,1),2,1),'color','k');
line(xlim,ylim,'color','k');

% Get zero-lag correlation
corr_diag = cell2mat(arrayfun(@(x) diag(corr_grid_cat(:,:,x)),1:size(corr_grid_cat,3),'uni',false))';

% Get local-lag correlation
max_lag = 0.05;
max_lag_samples = round(max_lag*(sample_rate/downsample_factor));

corr_diag_local = squeeze(max(cell2mat(permute(arrayfun(@(y) cell2mat(arrayfun(@(x) padarray(diag(corr_grid_cat(:,:,y),x), ...
    [abs(x),0],NaN,'post'),-max_lag_samples:max_lag_samples,'uni',false)),1:size(corr_grid_cat,3),'uni',false),[1,3,2])),[],2))';
corr_diag_norm = bsxfun(@rdivide,corr_diag_local,max(abs(corr_diag_local),[],2));

% Sort by PC1
[coeff,score,latent] = pca(corr_diag_norm);
[~,sort_idx] = sort(score(:,1));
sort_idx = 1:size(corr_diag_norm,1); % don't sort

figure;
imagesc(t_downsample_diff,[],corr_diag_local(sort_idx,:));
caxis([-1,1]);
colormap(brewermap([],'*RdBu'));
set(gca,'YTick',1:size(corr_diag_local,1),'YTickLabel',area_labels_grid_used(sort_idx));
line([0,0],ylim,'color','k');
title('Local-lag correlations')

% Get time-local forward-back correlation
max_lag = 0.08;
max_lag_samples = round(max_lag*(sample_rate/downsample_factor));

corr_grid_local = cell2mat(permute(arrayfun(@(y) cell2mat(arrayfun(@(x) padarray(diag(corr_grid_cat(:,:,y),x), ...
    [abs(x),0],NaN,'post'),-max_lag_samples:max_lag_samples,'uni',false)),1:size(corr_grid_cat,3),'uni',false),[1,3,2]));

corr_grid_local_diff = squeeze(nanmean(corr_grid_local(:,max_lag_samples+2:end,:),2) - ...
    nanmean(corr_grid_local(:,max_lag_samples:-1:1,:),2))';
corr_grid_local_diff_norm = bsxfun(@rdivide,corr_grid_local_diff,max(abs(corr_grid_local_diff),[],2));

[~,sort_idx] = sort(max(corr_grid_local_diff,[],2));

figure;
imagesc(t_downsample_diff,[],corr_grid_local_diff(sort_idx,:));
caxis([-max(abs(corr_grid_local_diff(:))),max(abs(corr_grid_local_diff(:)))]);
colormap(brewermap([],'*RdBu'));
set(gca,'YTick',1:size(corr_grid_local,1),'YTickLabel',area_labels_grid_used(sort_idx));
line([0,0],ylim,'color','k');
title('Forward-back correlations')


% Plot correlation vs activity mean (to see if one is just a function of
% the other)
mean_act = squeeze(nanmean(activity_cat(use_trials,:,:),1));
mean_act_geomean = sqrt(abs(bsxfun(@times,permute(mean_act,[2,3,1]),permute(mean_act,[3,2,1]))));

mean_act_geomean_pairs = ...
    cell2mat(arrayfun(@(x) AP_itril(mean_act_geomean(:,:,x),-1),1:length(t_downsample_diff),'uni',false));


compare_areas = [2,10];
curr_geo_mean = squeeze(mean_act_geomean(compare_areas(2),compare_areas(1),:));
curr_corr = diag(corr_grid{compare_areas(2),compare_areas(1)});

figure; 
subplot(1,2,1); hold on;
plot(t_downsample_diff,mean_act(:,compare_areas(1))/max(mean_act(:,compare_areas(1))),'linewidth',2);
plot(t_downsample_diff,mean_act(:,compare_areas(2))/max(mean_act(:,compare_areas(2))),'linewidth',2);
plot(t_downsample_diff,curr_corr,'linewidth',2);
xlabel('Time');
ylabel('Max-normalized')
legend([area_labels(compare_areas),{'Corr'}]);

subplot(1,2,2);hold on
t0 = find(t_downsample_diff > 0,1);
t1 = find(t_downsample_diff > 0.1,1);
plot(curr_geo_mean,curr_corr,'k','linewidth',2)
plot(curr_geo_mean(t0),curr_corr(t0),'ok','MarkerSize',10)
plot(curr_geo_mean(t1),curr_corr(t1),'^k','MarkerSize',10)
xlabel('Geometric mean activity');
ylabel('Correlation');



% Plot ctx/str average correlations
ctx_str_binary = [zeros(size(activity_cat,3)-4,1);ones(4,1)];
ctx_str_used = AP_itril(bsxfun(@plus,ctx_str_binary,ctx_str_binary'),-1);

ctx_ctx_corr_diag = corr_diag_local(ctx_str_used == 0,:);
ctx_str_corr_diag = corr_diag_local(ctx_str_used == 1,:);
str_str_corr_diag = corr_diag_local(ctx_str_used == 2,:);

figure;
subplot(4,3,[1,4,7]);
imagesc(t_downsample_diff,[],ctx_ctx_corr_diag);
caxis([-1,1]);
colormap(brewermap([],'*RdBu'));
line([0,0],ylim,'color','k');
title('Cortex-cortex correlations');
set(gca,'YTick',1:sum(ctx_str_used == 0),'YTickLabel', ...
    area_labels_grid_used(ctx_str_used == 0));

subplot(4,3,[2,5,8]);
imagesc(t_downsample_diff,[],ctx_str_corr_diag);
caxis([-1,1]);
colormap(brewermap([],'*RdBu'));
line([0,0],ylim,'color','k');
title('Cortex-striatum correlations');
set(gca,'YTick',1:sum(ctx_str_used == 1),'YTickLabel',...
    area_labels_grid_used(ctx_str_used == 1));

subplot(4,3,[3,6,9]);
imagesc(t_downsample_diff,[],str_str_corr_diag);
caxis([-1,1]);
colormap(brewermap([],'*RdBu'));
line([0,0],ylim,'color','k');
title('Striatum-striatum correlations');
set(gca,'YTick',1:sum(ctx_str_used == 2),'YTickLabel', ...
    area_labels_grid_used(ctx_str_used == 2));

subplot(4,3,[10,11,12]); hold on
plot(t_downsample_diff,nanmean(ctx_ctx_corr_diag,1),'linewidth',2);
plot(t_downsample_diff,nanmean(ctx_str_corr_diag,1),'linewidth',2);
plot(t_downsample_diff,nanmean(str_str_corr_diag,1),'linewidth',2);
line([0,0],ylim,'color','k');
legend({'Ctx-Ctx','Ctx-Str','Str-Str'})





%% Regress concatenated MUA from fluor, plot by reaction time

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_df_kernel-str_' num2str(n_aligned_depths) '_depths.mat'];
load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = framerate*upsample_factor;
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);


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

% Concantenate D
D_cellcat = vertcat(D_all{:});
D_cellcat = struct2cell(vertcat(D_cellcat{:}));
D_allcat = cell2struct(arrayfun(@(x) vertcat(D_cellcat{x,:}),1:size(D_cellcat,1),'uni',false)',fieldnames(D_all{1}{1}));

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_rois,2);
fluor_unilateral_allcat = nan(0,length(t),n_rois,2);
mua_allcat = nan(0,length(t),n_depths,2);
wheel_allcat = nan(0,length(t),2);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

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
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
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
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
    wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
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
    
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat_hemidiff_norm];
    fluor_unilateral_allcat = [fluor_unilateral_allcat;fluor_cat_norm];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];

end

% % Regress from fluor to MUA (all trials)
% % (by depth: data-less days cause different trial numbers)
% kernel_t = [-0.05,0.05];
% kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
% zs = [false,false];
% cvfold = 5;
% lambda = 0;
% 
% mua_nonan_trials = ~squeeze(any(any(isnan(mua_allcat),2),4));
% mua_allcat_predicted = nan(size(mua_allcat));
% for curr_depth = 1:n_depths
%     
%     curr_valid_trials = mua_nonan_trials(:,curr_depth);
%     
%     [~,predicted_spikes,~] = ...
%         AP_regresskernel(reshape(permute( ...
%         fluor_unilateral_allcat(curr_valid_trials,:,:,:), ...
%         [2,1,4,3]),[],n_rois)', ...
%         reshape(permute(mua_allcat(curr_valid_trials,:,curr_depth,:), ...
%         [2,1,4,3]),[],1)',kernel_frames,lambda,zs,cvfold);
%     
%     mua_allcat_predicted(curr_valid_trials,:,curr_depth,:) = ...
%         permute(reshape(predicted_spikes',length(t),sum(curr_valid_trials),2),[2,1,4,3]);
%     
% end

% Regress from fluor to MUA (kernel from non-visually guided trials)
kernel_t = [-0.1,0.1];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
zs = [false,false];
cvfold = 5;
lambda = 1e2;

% (define non-vis trials: 0 contrast and incorrect)
nonvis_trials = trial_contrast_allcat == 0 | ...
    (trial_contrast_allcat > 0 & trial_side_allcat == trial_choice_allcat);

mua_nonan_trials = ~squeeze(any(any(isnan(mua_allcat),2),4));
mua_allcat_predicted = nan(size(mua_allcat));
nonvis_k_reshape = nan(n_rois,length(kernel_frames),n_depths);
for curr_depth = 1:n_depths
    
    curr_nonan_trials = mua_nonan_trials(:,curr_depth);
    curr_kernel_trials = nonvis_trials & mua_nonan_trials(:,curr_depth);
    
    [nonvis_k,~,explained_var] = ...
        AP_regresskernel(reshape(permute( ...
        fluor_unilateral_allcat(curr_kernel_trials,:,:,:), ...
        [2,1,4,3]),[],n_rois)', ...
        reshape(permute(mua_allcat(curr_kernel_trials,:,curr_depth,:), ...
        [2,1,4,3]),[],1)',kernel_frames,lambda,zs,cvfold);
    
    nonvis_k_reshape(:,:,curr_depth) = reshape(nonvis_k,16,[]);
       
    % Create design matrix of all time-shifted regressors
    regressor_design = repmat(reshape(permute( ...
        fluor_unilateral_allcat(curr_nonan_trials,:,:,:), ...
        [2,1,4,3]),[],n_rois), ...
        [1,1,length(kernel_frames)]);
    
    % Temporally shift each page
    for curr_kernel_frame = 1:length(kernel_frames)
        regressor_design(:,:,curr_kernel_frame) = ...
            circshift(regressor_design(:,:,curr_kernel_frame), ...
            [kernel_frames(curr_kernel_frame),0,0]);
    end
    
    regressor_design = ...
        reshape(regressor_design,[],size(regressor_design,2)*size(regressor_design,3)); 
    
    predicted_spikes = regressor_design*nonvis_k;
    
    mua_allcat_predicted(curr_nonan_trials,:,curr_depth,:) = ...
        permute(reshape(predicted_spikes',length(t),sum(curr_nonan_trials),2),[2,1,4,3]);
    
end

% Get reaction time per trial
% (by position)
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
% (by max speed)
% [~,move_idx] = max(abs(diff(wheel_allcat(:,:,1),[],2)),[],2);

move_t = t(move_idx);

% Plot average activity within condition restricted by reaction time
rxn_times = linspace(0,0.5,6);
n_rxn_times = length(rxn_times)-1;
plot_align = 1;

activity_split = cell(length(contrasts),n_rxn_times,n_depths);
activity_predicted_split = cell(length(contrasts),n_rxn_times,n_depths);
for curr_depth = 1:n_depths
    for curr_contrast = 1:length(contrasts)
        
        if contrasts(curr_contrast) == 0
            use_sides = true(size(trial_side_allcat));
        else
            use_sides = trial_side_allcat == 1;
        end
        
        activity_split(curr_contrast,:,curr_depth) = ...
            arrayfun(@(x) mua_allcat( ...
            move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
            use_sides & ...
            trial_contrast_allcat == contrasts(curr_contrast) & ...
            trial_choice_allcat == -1,:,curr_depth,plot_align), ...
            1:length(rxn_times)-1,'uni',false);
        
        activity_predicted_split(curr_contrast,:,curr_depth) = ...
            arrayfun(@(x) mua_allcat_predicted( ...
            move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
            use_sides & ...
            trial_contrast_allcat == contrasts(curr_contrast) & ...
            trial_choice_allcat == -1,:,curr_depth,plot_align), ...
            1:length(rxn_times)-1,'uni',false);
        
    end
end

n_boot = 1000;

activity_boot = cellfun(@ (act) ...
    prctile(bootstrp(n_boot,@(x) nanmean(x,1),act),[50,2.5,97.5]), ...
    activity_split,'uni',false);

activity_predicted_boot = cellfun(@ (act) ...
    prctile(bootstrp(n_boot,@(x) nanmean(x,1),act),[50,2.5,97.5]), ...
    activity_predicted_split,'uni',false);

activity_residual_boot = cellfun(@ (act,act_pred) ...
    prctile(bootstrp(n_boot,@(x) nanmean(x,1),act-act_pred),[50,2.5,97.5]), ...
    activity_split,activity_predicted_split,'uni',false);

plot_mua = 1;
figure;
plot_col = copper(n_rxn_times);
for curr_contrast = 1:length(contrasts)
    
    subplot(length(contrasts),3,3*(curr_contrast-1)+1); hold on;
    p = arrayfun(@(x) AP_errorfill(t,activity_boot{curr_contrast,x,plot_mua}(1,:)', ...
        bsxfun(@minus,activity_boot{curr_contrast,x,plot_mua}(2:3,:), ...
        activity_boot{curr_contrast,x,plot_mua}(1,:))', ...
        plot_col(x,:),0.5),1:n_rxn_times,'uni',false);
    for i = 1:length(rxn_times)-1
        line(repmat(rxn_times(i),2,1),ylim,'color',plot_col(i,:));
    end
    ylabel('Measured');
    title(contrasts(curr_contrast));
    
    subplot(length(contrasts),3,3*(curr_contrast-1)+2); hold on;
    p = arrayfun(@(x) AP_errorfill(t,activity_predicted_boot{curr_contrast,x,plot_mua}(1,:)', ...
        bsxfun(@minus,activity_predicted_boot{curr_contrast,x,plot_mua}(2:3,:), ...
        activity_predicted_boot{curr_contrast,x,plot_mua}(1,:))', ...
        plot_col(x,:),0.5),1:n_rxn_times,'uni',false);
    for i = 1:length(rxn_times)-1
        line(repmat(rxn_times(i),2,1),ylim,'color',plot_col(i,:));
    end
    ylabel('Predicted');
    title(contrasts(curr_contrast));
    
    subplot(length(contrasts),3,3*(curr_contrast-1)+3); hold on;
    p = arrayfun(@(x) AP_errorfill(t,activity_residual_boot{curr_contrast,x,plot_mua}(1,:)', ...
        bsxfun(@minus,activity_residual_boot{curr_contrast,x,plot_mua}(2:3,:), ...
        activity_residual_boot{curr_contrast,x,plot_mua}(1,:))', ...
        plot_col(x,:),0.5),1:n_rxn_times,'uni',false);
    for i = 1:length(rxn_times)-1
        line(repmat(rxn_times(i),2,1),ylim,'color',plot_col(i,:));
    end
    ylabel('Residual');
    title(contrasts(curr_contrast));
    
end

% Plot residuals together across contrasts, compare to V1 latency
v1_split = cell(length(contrasts),n_rxn_times);
for curr_contrast = 1:length(contrasts)
    
    if contrasts(curr_contrast) == 0
        use_sides = true(size(trial_side_allcat));
    else
        use_sides = trial_side_allcat == 1;
    end
    
    v1_split(curr_contrast,:) = ...
        arrayfun(@(x) fluor_allcat( ...
        move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
        use_sides & ...
        trial_contrast_allcat == contrasts(curr_contrast) & ...
        trial_choice_allcat == -1,:,1,1), ...
        1:length(rxn_times)-1,'uni',false);   
    
end

v1_boot = cellfun(@ (act) ...
    prctile(bootstrp(n_boot,@(x) nanmean(x,1),act),[50]), ...
    v1_split,'uni',false);

use_t = t > 0 & t < 0.15;
plot_residuals = cellfun(@(x) x(1,:),activity_residual_boot,'uni',false);
[max_residuals,max_residuals_idx] = cellfun(@(x) max(x(:,use_t),[],2),plot_residuals);
[max_v1,max_v1_idx] = cellfun(@(x) max(x(:,use_t),[],2),v1_boot);

t_used = t(use_t);
max_residuals_t = t_used(max_residuals_idx);
max_v1_t = t_used(max_v1_idx)';

% % (to bootstrap max - doesn't really change it?)
max_residuals_idx_boot = cellfun(@ (act,act_pred) ...
    median(bootstrp(n_boot,@(x) find(nanmean(x(:,use_t),1) == ...
    max(nanmean(x(:,use_t),1),[],2),1),act-act_pred)), ...
    activity_split,activity_predicted_split);
max_residuals_t = t_used(round(max_residuals_idx_boot));

%(this part will break if more than one reaction time bin)
figure; 
subplot(n_depths+1,2,1); hold on;
set(gca,'ColorOrder',copper(length(contrasts)));
plot(t,bsxfun(@rdivide,cell2mat(v1_boot),max_v1),'linewidth',2);
title('V1 (peak-normalized)');

for curr_depth = 1:n_depths
    subplot(n_depths+1,2,3+(curr_depth-1)*2); hold on;
    set(gca,'ColorOrder',copper(length(contrasts)));
    plot(t,bsxfun(@rdivide,cell2mat(plot_residuals(:,:,curr_depth)), ...
        max_residuals(:,:,curr_depth)),'linewidth',2);
    title(['Str ' num2str(curr_depth) ' (peak-normalized)']);
end

subplot(1,2,2); hold on;
set(gca,'ColorOrder',[1,0,0;copper(4);copper(4)]);
plot(contrasts(2:end),max_v1_t(2:end),'linewidth',2);
plot(contrasts(2:end),squeeze(max_residuals_t(2:end,:,:)),'linewidth',2);
plot(contrasts(2:end),squeeze(bsxfun(@minus,max_residuals_t(2:end,:,:),max_v1_t(2:end))),'-.','linewidth',2);
xlabel('Contrast');
ylabel('Time to activity peak');

% Plot weights by spatial ROI
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

roi_cat = cat(3,wf_roi.mask);
use_kernel_t = kernel_frames >= -4 & kernel_frames <= 1;

figure;
for curr_depth = 1:n_depths   
    
    curr_k = nonvis_k_reshape(:,use_kernel_t,curr_depth);
    [~,max_ampl_weight_idx] = max(abs(curr_k),[],2);
    roi_weights = arrayfun(@(x) curr_k(x,max_ampl_weight_idx(x)),1:n_rois);
        
    subplot(1,n_depths,curr_depth); hold on
    set(gca,'YDir','reverse');
    AP_reference_outline('ccf_aligned','k');
    for curr_roi = 1:n_rois
        curr_roi_boundary = cell2mat(bwboundaries(roi_cat(:,:,curr_roi)));
        patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),roi_weights(curr_roi));
    end
    axis image off;
    colormap(colormap_BlueWhiteRed);
    caxis([-max(abs(caxis)),max(abs(caxis))]);
end



%% Regress move starts from all activity within condition
% if there's really a different kind of activity for spontaneous vs
% visually-guided, then should show up in kernel for different reaction
% time bins for a given contrast
%
% (this didn't really show anything obviously interesting)

n_aligned_depths = 4;
load_data = 'early';

% Load data (early/late/all)
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
early_data_fn = ['all_trial_activity_df_kernel-str_earlymove_' num2str(n_aligned_depths) '_depths.mat'];
late_data_fn = ['all_trial_activity_df_kernel-str_latemove_' num2str(n_aligned_depths) '_depths.mat'];

switch load_data
    case 'early'
        load([data_path filesep early_data_fn])
    case 'late'
        load([data_path filesep late_data_fn])
    case 'all'
        early_data = load([data_path filesep early_data_fn]);
        late_data = load([data_path filesep late_data_fn]);
        
        n_animals = length(early_data.D_all);
        
        D_all = cell(n_animals,1);
        fluor_all = cell(n_animals,1);
        mua_all = cell(n_animals,1);
        wheel_all = cell(n_animals,1);
        for curr_animal = 1:n_animals
            for curr_day = 1:length(early_data.D_all{curr_animal})
                D_all{curr_animal}{curr_day,1} = ...
                    cell2struct(cellfun(@vertcat,struct2cell(early_data.D_all{curr_animal}{curr_day}), ...
                    struct2cell(late_data.D_all{curr_animal}{curr_day}),'uni',0),fieldnames(early_data.D_all{curr_animal}{curr_day}),1);
                
                fluor_all{curr_animal}{curr_day,1} = ...
                    [early_data.fluor_all{curr_animal}{curr_day}; ...
                    late_data.fluor_all{curr_animal}{curr_day}];
                
                mua_all{curr_animal}{curr_day,1} = ...
                    [early_data.mua_all{curr_animal}{curr_day}; ...
                    late_data.mua_all{curr_animal}{curr_day}];
                
                wheel_all{curr_animal}{curr_day,1} = ...
                    [early_data.wheel_all{curr_animal}{curr_day}; ...
                    late_data.wheel_all{curr_animal}{curr_day}];
                
            end
        end
end

% Get time
framerate = 35;
raster_window = [-0.5,1];
upsample_factor = 3;
sample_rate = framerate*upsample_factor;
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);


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

% Concantenate D
D_cellcat = vertcat(D_all{:});
D_cellcat = struct2cell(vertcat(D_cellcat{:}));
D_allcat = cell2struct(arrayfun(@(x) vertcat(D_cellcat{x,:}),1:size(D_cellcat,1),'uni',false)',fieldnames(D_all{1}{1}));

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_rois,2);
fluor_unilateral_allcat = nan(0,length(t),n_rois,2);
mua_allcat = nan(0,length(t),n_depths,2);
wheel_allcat = nan(0,length(t),2);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

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
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
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
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
    wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
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
    
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat_hemidiff_norm];
    fluor_unilateral_allcat = [fluor_unilateral_allcat;fluor_cat_norm];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];

end

% Get reaction time per trial
% (by position)
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
% (by max speed)
% [~,move_idx] = max(abs(diff(wheel_allcat(:,:,1),[],2)),[],2);

move_t = t(move_idx);

move_starts = zeros(size(wheel_allcat,1),size(wheel_allcat,2),2);
for curr_trial = 1:size(move_starts,1)
    curr_choice = trial_choice_allcat(curr_trial);
    move_starts(curr_trial,move_idx(curr_trial),curr_choice/2+1.5) = 1;
end

% Plot average activity within condition restricted by reaction time
use_condition = ...
    trial_side_allcat == 1 & ...
    trial_contrast_allcat == 0.06 & ...
    trial_choice_allcat == -1;

rxn_times = linspace(0,0.5,6);

fluor_split = arrayfun(@(x) ...
    fluor_allcat( ... 
    move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
    use_condition,:,:,1), ...
    1:length(rxn_times)-1,'uni',false);

mua_split = arrayfun(@(x) ...
    mua_allcat( ... 
    move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
    use_condition,:,:,1), ...
    1:length(rxn_times)-1,'uni',false);

activity_split = cellfun(@(x,y) cat(3,x,y),fluor_split,mua_split,'uni',false);

move_starts_split = arrayfun(@(x) ...
    move_starts( ... 
    move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
    use_condition,:,:,1), ...
    1:length(rxn_times)-1,'uni',false);

% % Do regression together
% t_shifts = repmat({-30:30},2,1);
% lambda = 1e3;
% zs = [true,true];
% cvfold = 5;
% 
% kernel = cell(length(rxn_times)-1,1);
% for curr_rxn = 1:length(rxn_times)-1
%     
%     % this isn't great because it restricts to only when all mua was
%     % recorded
%     use_trials = all(any(fluor_split{curr_rxn},2),3) & ...
%         all(any(mua_split{curr_rxn},2),3);
%     
%     % Set activity
%     activity_regressors = {reshape(permute(fluor_split{curr_rxn}(use_trials,:,1:n_rois/2),[2,1,3]),[],n_rois/2)', ...
%         reshape(permute(mua_split{curr_rxn}(use_trials,:,:),[2,1,3]),[],n_depths)'}';
%     
%     move_regress = reshape(permute(move_starts_split{curr_rxn}(use_trials,:,1),[2,1,3]),[],1)';
%         
%     [kernel_flat,move_predicted,expl_var,move_predicted_reduced] = AP_regresskernel( ...
%         activity_regressors,move_regress,t_shifts,lambda,zs,cvfold);
%    
%     kernel{curr_rxn} = cellfun(@(x,t) reshape(x,[],length(t)), ...
%         mat2cell(kernel_flat,cellfun(@(x) size(x,1),activity_regressors).*cellfun(@length,t_shifts),1), ...
%         t_shifts,'uni',false);
%     
% end
% 
% figure;
% for curr_rxn = 1:length(rxn_times)-1
%     subplot(length(rxn_times)-1,1,curr_rxn); hold on;
%     set(gca,'ColorOrder',[jet(8);copper(4)]);
%     plot(vertcat(kernel{curr_rxn}{:})','linewidth',2);
% end

% Do regression separately 
t_shifts = -30:30;
lambda = 1e3;
zs = [true,true];
cvfold = 5;

kernel = cell(length(rxn_times)-1,1);
expl_var_all = cell(length(rxn_times)-1,1);
for curr_rxn = 1:length(rxn_times)-1
    for curr_area = 1:20
            
    use_trials = any(activity_split{curr_rxn}(:,:,curr_area),2);   
    
    % Set activity
    activity_regressors = reshape(permute(activity_split{curr_rxn}(use_trials,:,curr_area),[2,1,3]),[],1)';
    
    move_regress = reshape(permute(move_starts_split{curr_rxn}(use_trials,:,1),[2,1,3]),[],1)';
        
    [curr_kernel,move_predicted,expl_var,move_predicted_reduced] = AP_regresskernel( ...
        activity_regressors,move_regress,t_shifts,lambda,zs,cvfold);
   
    kernel{curr_rxn}(:,curr_area) = curr_kernel;
    expl_var_all{curr_rxn}(curr_area) = expl_var.total;
    
    end
end

figure;
for curr_rxn = 1:length(rxn_times)-1
    subplot(length(rxn_times)-1,1,curr_rxn); hold on;
    set(gca,'ColorOrder',[jet(8);copper(4)]);
    plot(kernel{curr_rxn},'linewidth',2);
end

figure; hold on
set(gca,'ColorOrder',copper(length(rxn_times)-1));
plot(vertcat(expl_var_all{:})','linewidth',2)
xlabel('Area');
set(gca,'XTick',1:20,'XTickLabel',[{wf_roi.area},num2cell(1:4)])
ylabel('Explained var');
legend(arrayfun(@(x) [num2str(rxn_times(x)) '-' num2str(rxn_times(x+1))], ...
    1:length(rxn_times)-1,'uni',false));


%%% One reaction time bin broken down by contrast
rxn_times = linspace(0,0.5,2);
n_rxn_times = length(rxn_times)-1;

fluor_split = cell(length(contrasts),n_rxn_times);
mua_split = cell(length(contrasts),n_rxn_times);
move_starts_split = cell(length(contrasts),n_rxn_times);
for curr_contrast = 1:length(contrasts)
    
    if contrasts(curr_contrast) == 0
        use_sides = true(size(trial_side_allcat));
    else
        use_sides = trial_side_allcat == 1;
    end
    
    fluor_split(curr_contrast,:) = ...
        arrayfun(@(x) fluor_allcat( ...
        move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
        use_sides & ...
        trial_contrast_allcat == contrasts(curr_contrast) & ...
        trial_choice_allcat == -1,:,:,1), ...
        1:length(rxn_times)-1,'uni',false);
    
    mua_split(curr_contrast,:) = ...
        arrayfun(@(x) mua_allcat( ...
        move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
        use_sides & ...
        trial_contrast_allcat == contrasts(curr_contrast) & ...
        trial_choice_allcat == -1,:,:,1), ...
        1:length(rxn_times)-1,'uni',false);
    
    move_starts_split(curr_contrast,:) = ...
        arrayfun(@(x) move_starts( ...
        move_t' > rxn_times(x) & move_t' < rxn_times(x+1) & ...
        use_sides & ...
        trial_contrast_allcat == contrasts(curr_contrast) & ...
        trial_choice_allcat == -1,:,:,1), ...
        1:length(rxn_times)-1,'uni',false);
    
end

activity_split = cellfun(@(x,y) cat(3,x,y),fluor_split,mua_split,'uni',false);

% Do regression separately 
t_shifts = -30:30;
lambda = 1e3;
zs = [true,true];
cvfold = 5;

kernel = cell(length(contrasts),n_rxn_times);
expl_var_all = cell(length(contrasts),n_rxn_times);
for curr_contrast = 1:length(contrasts)  
    for curr_rxn = 1:n_rxn_times
        for curr_area = 1:20
            
            use_trials = any(activity_split{curr_contrast,curr_rxn}(:,:,curr_area),2);
            
            % Set activity
            activity_regressors = reshape(permute(activity_split{curr_contrast,curr_rxn}(use_trials,:,curr_area),[2,1,3]),[],1)';
            
            move_regress = reshape(permute(move_starts_split{curr_contrast,curr_rxn}(use_trials,:,1),[2,1,3]),[],1)';
            
            [curr_kernel,move_predicted,expl_var,move_predicted_reduced] = AP_regresskernel( ...
                activity_regressors,move_regress,t_shifts,lambda,zs,cvfold);
            
            kernel{curr_contrast,curr_rxn}(:,curr_area) = curr_kernel;
            expl_var_all{curr_contrast,curr_rxn}(curr_area) = expl_var.total;
            
        end
    end
    disp(curr_contrast);
end

expl_var_grid = cell2mat(cellfun(@(x) permute(x,[1,3,2]),expl_var_all,'uni',false));

col = colormap_BlueWhiteRed(11);
col = [0,0,0;col(end-4:end,:)];
figure; hold on
set(gca,'ColorOrder',col);
plot(squeeze(expl_var_grid)','linewidth',2);
ylabel('Explained variance');
xlabel('Area');
set(gca,'XTick',1:20,'XTickLabel',[{wf_roi.area},cellfun(@(x) ...
    ['Str ' num2str(x)],num2cell(1:4),'uni',false)]);
legend(cellfun(@(x) ['Contrast ' num2str(x)],num2cell(contrasts),'uni',false));


% Regress move starts from activity on vis/non-vis

% (define vis/non-vis trials)
vis_trials = trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == -1;

nonvis_trials = (trial_contrast_allcat == 0 & ...
    trial_choice_allcat == -1) | ...
    (trial_contrast_allcat > 0 & ...
    trial_side_allcat == -1 & ...
    trial_choice_allcat == -1);

trial_types = [vis_trials,nonvis_trials];

%%% One reaction time bin broken down by contrast
fluor_split = ...
    arrayfun(@(x) fluor_allcat( ...
    trial_types(:,x),:,:,1), ...
    1:size(trial_types,2),'uni',false);

mua_split = ...
    arrayfun(@(x) mua_allcat( ...
    trial_types(:,x),:,:,1), ...
    1:size(trial_types,2),'uni',false);

move_starts_split = ...
    arrayfun(@(x) move_starts( ...
    trial_types(:,x),:,:,1), ...
    1:size(trial_types,2),'uni',false);

activity_split = cellfun(@(x,y) cat(3,x,y),fluor_split,mua_split,'uni',false);

% Do regression separately 
t_shifts = -30:30;
lambda = 1e2;
zs = [true,false];
cvfold = 5;

kernel = cell(size(trial_types,2),1);
expl_var_all = cell(size(trial_types,2),1);
for curr_trials = 1:size(trial_types,2)
        for curr_area = 1:20
            
            use_trials = any(activity_split{curr_trials}(:,:,curr_area),2);
            
            % Set activity
            activity_regressors = reshape(permute(activity_split{curr_trials}(use_trials,:,curr_area),[2,1,3]),[],1)';
            
            move_regress = reshape(permute(move_starts_split{curr_trials}(use_trials,:,1),[2,1,3]),[],1)';
            
            [curr_kernel,move_predicted,expl_var,move_predicted_reduced] = AP_regresskernel( ...
                activity_regressors,move_regress,t_shifts,lambda,zs,cvfold);
            
            kernel{curr_trials}(:,curr_area) = curr_kernel;
            expl_var_all{curr_trials}(curr_area) = expl_var.total;
            
            AP_print_progress_fraction(curr_area,20);
        end
end

expl_var_grid = cell2mat(cellfun(@(x) permute(x,[1,3,2]),expl_var_all,'uni',false));

col = [1,0,0;0,0,1];
figure; hold on
set(gca,'ColorOrder',col);
plot(squeeze(expl_var_grid)','linewidth',2);
ylabel('Explained variance');
xlabel('Area');
set(gca,'XTick',1:20,'XTickLabel',[{wf_roi.area},cellfun(@(x) ...
    ['Str ' num2str(x)],num2cell(1:4),'uni',false)]);
legend({'Vis','Non-vis'});

% plot kernels
cax = [-max(abs(reshape([kernel{:}],[],1))),max(abs(reshape([kernel{:}],[],1)))];
figure; 
subplot(1,3,1);
imagesc([],t_shifts/sample_rate,kernel{1});
xlabel('Area');
caxis(cax);
subplot(1,3,2);
imagesc([],t_shifts/sample_rate,kernel{2});
xlabel('Area');
caxis(cax);
subplot(1,3,3);
imagesc([],t_shifts/sample_rate,kernel{1}-kernel{2});
xlabel('Area');
caxis(cax);
colormap(colormap_BlueWhiteRed);

% get difference between kernels
kernel_diff = kernel{1} - kernel{2};
figure; hold on
set(gca,'ColorOrder',[jet(8);copper(4)]);
plot(t_shifts/sample_rate,kernel_diff(:,[1:8,17:20]),'linewidth',2);
line([0,0],ylim,'color','k');
legend([{wf_roi(:,1).area},cellfun(@(x) ...
    ['Str ' num2str(x)],num2cell(1:4),'uni',false)]);


%% Regress concatenated fluor -> mua kernel in V-space

% Load data
use_data = 'task';

n_aligned_depths = 4;
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
switch use_data
    case 'task'
        data_fn = ['all_trial_activity_Udf_kernel-str_4_depths_long'];
        load_task = true;
    case 'passive'
        data_fn = ['all_trial_activity_Udf_kernel-str_passive_fullscreen_4_depths'];
%         data_fn = ['all_trial_activity_Udf_kernel-str_passive_choiceworld_4_depths'];
        load_task = false;
end

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);
if load_task
    reward_all = cellfun(@(data,use_expts) data(use_expts),reward_all,use_experiments','uni',false);
end

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);
if load_task
    reward_all = reward_all(use_animals);
end

% Get widefield ROIs
n_rois = 200;

% MUA depths
n_depths = size(mua_all{1}{1},3);

% Set up conditions
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_rois,1);
mua_allcat = nan(0,length(t),n_depths,1);
wheel_allcat = nan(0,length(t),1);
reward_allcat = nan(0,length(t),1);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

D_cat = struct('stimulus',cell(1,1),'response',cell(1,1),'repeatNum',cell(1,1));

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    n_align = size(fluor_cat,4);
    
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    t_std = t > -0.5 & t < 0.5;
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x(:,t_std,:,1),[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
    
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
    %     wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
        
    if load_task            
        % Concatenate behavioural data
        D = struct;
        D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
        D.response = cell2mat(cellfun(@(x) x.response,D_all{curr_animal},'uni',false));
        D.repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all{curr_animal},'uni',false));
        
        D.day = trial_day;
        
        D_cat.stimulus = vertcat(D_cat.stimulus,D.stimulus);
        D_cat.response = vertcat(D_cat.response,D.response);
        D_cat.repeatNum = vertcat(D_cat.repeatNum,D.repeatNum);
        
        % Get trial ID
        trial_contrast = max(D.stimulus,[],2);
        [~,side_idx] = max(D.stimulus > 0,[],2);
        trial_side = (side_idx-1.5)*2;
        trial_choice = -(D.response-1.5)*2;
        
        trial_conditions = ...
            [trial_contrast, trial_side, ...
            trial_choice, ones(size(trial_day))];
        [~,trial_id] = ismember(trial_conditions,conditions,'rows');
    else
        % Concatenate stim
        D = struct;
        D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));     
        D_cat.stimulus = vertcat(D_cat.stimulus,D.stimulus);
    end
    
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    if load_task
        reward_allcat = [reward_allcat;vertcat(reward_all{curr_animal}{:})];        
        trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
        trial_side_allcat = [trial_side_allcat;trial_side];
        trial_choice_allcat = [trial_choice_allcat;trial_choice];
    end
    
end

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% %%% REMOVE MOVE TRIALS FOR PASSIVE
% move_trials = any(abs(wheel_allcat(:,t >= 0 & t <= 2) > 2),2);
% fluor_allcat(move_trials,:,:) = [];
% mua_allcat(move_trials,:,:) = [];
% wheel_allcat(move_trials,:,:) = [];
% D_cat.stimulus(move_trials) = [];
% %%%

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Downsample (otherwise it's too much for regression) and d(smooth(fluor))
downsample_factor = 4;
t_downsample = linspace(t(1),t(end),round(length(t)/downsample_factor));

t_diff =  conv(t,[1,1]/2,'valid');
t_downsample_diff = conv(t_downsample,[1,1]/2,'valid');

if load_task
    wheel = interp1(t,wheel_allcat(:,:,1)',t_downsample)';
    [~,move_idx] = max(abs(wheel(:,:,1)) > 2,[],2);
end

% (old - interp doesn't make sense b/c binned)
% mua_allcat_downsamp = permute(interp1(t,permute(mua_allcat,[2,1,3,4]),t_downsample),[2,1,3,4]);
% (new - sum)
mua_allcat_downsamp = permute(reshape(squeeze(sum(reshape(reshape(permute( ...
    mua_allcat(:,1:end-mod(size(mua_allcat,2),downsample_factor),:), ...
    [2,1,3]),[],n_depths),downsample_factor,[],n_depths),1)), ...
    [],size(mua_allcat,1),n_depths),[2,1,3]);

% % To use straight fluorescence (subtract t < 0)
% smooth_factor = 1;
% fluor_allcat_downsamp = permute(interp1(t,permute(convn( ...
%     fluor_allcat,ones(1,smooth_factor)/smooth_factor,'same'),[2,1,3,4]),t_downsample),[2,1,3,4]);
% fluor_allcat_downsamp = bsxfun(@minus,fluor_allcat_downsamp, ...
%     nanmedian(fluor_allcat_downsamp(:,t_downsample < 0,:),2));

% To use derivative
smooth_factor = 3;
fluor_allcat_downsamp = permute(interp1(t_diff,permute(diff(convn( ...
    fluor_allcat,ones(1,smooth_factor)/smooth_factor,'same'),[],2),[2,1,3,4]),t_downsample),[2,1,3,4]);

% (there's a nan in one trial??)
fluor_allcat_downsamp(isnan(fluor_allcat_downsamp)) = 0;

% Low-pass filter activity (where does this ~10 Hz crap come from?)
lowpassCutoff = 10; % Hz
[b100s, a100s] = butter(2, lowpassCutoff/((sample_rate/downsample_factor)/2), 'low');
fluor_allcat_downsamp_filt = filter(b100s,a100s,fluor_allcat_downsamp,[],2);
mua_allcat_downsamp_filt = filter(b100s,a100s,mua_allcat_downsamp,[],2);

% Regress from fluor to MUA (kernel from non-visually guided trials)
kernel_t = [-0.08,0.08];
kernel_frames = round(kernel_t(1)*sample_rate/downsample_factor):round(kernel_t(2)*sample_rate/downsample_factor);
zs = [false,false];
cvfold = 5;
lambda = 0;
return_constant = true;

% Set trials to do regression on (kernel then applied to all trials)
% (to use all - for cv)
kernel_trials = true(size(fluor_allcat,1),1);
% (to only use incorrect/zero trials)
% kernel_trials = trial_contrast_allcat == 0 | ...
%     trial_side_allcat == trial_choice_allcat;
% (to only use particular choice trials)
% kernel_trials = trial_choice_allcat == 1;
% (to only use timing)
% kernel_trials = move_t > 0.6 & move_t < 0.7;

mua_nonan_trials = ~squeeze(any(any(isnan(mua_allcat),2),4));
mua_allcat_predicted = nan(size(mua_allcat_downsamp_filt));
fluor_k_reshape = nan(n_rois,length(kernel_frames),n_depths);
for curr_depth = 1:n_depths
    
    curr_nonan_trials = mua_nonan_trials(:,curr_depth);
    curr_kernel_trials = kernel_trials & mua_nonan_trials(:,curr_depth);
    
    [fluor_k,curr_mua_predicted,explained_var] = ...
        AP_regresskernel(reshape(permute( ...
        fluor_allcat_downsamp_filt(curr_kernel_trials,:,:,:), ...
        [2,1,4,3]),[],n_rois)', ...
        reshape(permute(mua_allcat_downsamp_filt(curr_kernel_trials,:,curr_depth,:), ...
        [2,1,4,3]),[],1)',kernel_frames,lambda,zs,cvfold,return_constant);
    
    fluor_k_reshape(:,:,curr_depth) = reshape(fluor_k(1:end-1),n_rois,[]);
       
    % To use predicted straight from regression (cross-validated)
    mua_allcat_predicted(curr_kernel_trials,:,curr_depth,:) = ...
        permute(reshape(curr_mua_predicted',length(t_downsample),sum(curr_kernel_trials)),[2,1,4,3]);
    
%     % To apply kernel from regressed subset to all (CV NOT APPLICABLE)
%     % Create design matrix of all time-shifted regressors
%     regressor_design = repmat(reshape(permute( ...
%         fluor_allcat_downsamp_filt(curr_nonan_trials,:,:,:), ...
%         [2,1,4,3]),[],n_rois), ...
%         [1,1,length(kernel_frames)]);
%     
%     % Temporally shift each page
%     for curr_kernel_frame = 1:length(kernel_frames)
%         regressor_design(:,:,curr_kernel_frame) = ...
%             circshift(regressor_design(:,:,curr_kernel_frame), ...
%             [kernel_frames(curr_kernel_frame),0,0]);
%     end
%     
%     regressor_design = ...
%         [reshape(regressor_design,[],size(regressor_design,2)*size(regressor_design,3)) ...
%         ones(size(regressor_design,1),1)]; 
%     
%     predicted_spikes = regressor_design*fluor_k;
%     
%     mua_allcat_predicted(curr_nonan_trials,:,curr_depth,:) = ...
%         permute(reshape(predicted_spikes',length(t_downsample),sum(curr_nonan_trials)),[2,1,4,3]);
    
    AP_print_progress_fraction(curr_depth,n_depths);
    
end

fluor_k_px = reshape(svdFrameReconstruct(U_master(:,:,1:200),reshape(fluor_k_reshape,200,[])), ...
    size(U_master,1),size(U_master,2),length(kernel_frames),[]);

AP_imscroll(fluor_k_px);
caxis([-prctile(abs(fluor_k_px(:)),100),prctile(abs(fluor_k_px(:)),100)]);
colormap(brewermap([],'*RdBu'));
axis image;
AP_reference_outline('ccf_aligned','k');

fluor_k_px_max = squeeze(max(fluor_k_px(:,:,kernel_frames == 0,:),[],3));
figure;
for curr_depth = 1:n_depths
    subplot(1,n_depths,curr_depth)
    imagesc(fluor_k_px_max(:,:,curr_depth))
    AP_reference_outline('ccf_aligned','k');
    colormap(brewermap([],'*RdBu'));
    caxis([-prctile(abs(fluor_k_px_max(:)),99),prctile(abs(fluor_k_px_max(:)),99)]);
    axis image off;
end

% Apply empirical static nonlinearity
figure;
mua_allcat_predicted_nlin = mua_allcat_predicted;
for curr_depth = 1:n_depths       
    measured_data = reshape(mua_allcat_downsamp_filt(:,:,curr_depth),[],1);
    predicted_data = reshape(mua_allcat_predicted(:,:,curr_depth),[],1);
    
    n_bins = 500;
    activity_bounds = linspace(-1,6,n_bins+1);
    activity_bin_centers = conv2(activity_bounds,[1,1]/2,'valid');
  
    measured_bins = discretize(measured_data,activity_bounds);
    predicted_bins = discretize(predicted_data,activity_bounds);
    
    measured_data_binmean = accumarray(predicted_bins(~isnan(predicted_bins)), ...
        measured_data(~isnan(predicted_bins)),[n_bins,1],@median,nan);
    predicted_data_binmean = accumarray(predicted_bins(~isnan(predicted_bins)),...
        predicted_data(~isnan(predicted_bins)),[n_bins,1],@median,nan);
    
    % smooth out the measured data binmean
    measured_data_binmean_smooth = medfilt1(measured_data_binmean,10);
    
    predicted_data_nlin = nan(size(predicted_data));
    predicted_data_nlin(~isnan(predicted_bins)) = measured_data_binmean_smooth(predicted_bins(~isnan(predicted_bins)));
    
    predicted_data_nlin_bins = discretize(predicted_data_nlin,activity_bounds);
    predicted_data_nlin_binmean = accumarray( ...
        predicted_bins(~isnan(predicted_bins)), ...
        predicted_data_nlin(~isnan(predicted_bins)),[n_bins,1],@mean,nan);
    
    mua_allcat_predicted_nlin(:,:,curr_depth) = ...
        reshape(predicted_data_nlin, ...
        size(mua_allcat_predicted,1),size(mua_allcat_predicted,2));
    
    subplot(1,n_depths,curr_depth); hold on;
    plot(predicted_data,measured_data,'.')
    plot(predicted_data_binmean,measured_data_binmean,'linewidth',2);
    plot(predicted_data_binmean,measured_data_binmean_smooth,'linewidth',2);
    plot(predicted_data_nlin_binmean,measured_data_binmean,'linewidth',2);
    xlim([-2,5]);ylim([-2,5]);
    line(xlim,ylim,'color','k');
    xlabel('Predicted')
    ylabel('Measured')
    axis square;
end

% Get explained variance within time
mua_allcat_residual = mua_allcat_downsamp_filt - mua_allcat_predicted_nlin;

t_var = t_downsample > -inf & t_downsample < inf;
mua_sse_measured = squeeze(nansum(mua_allcat_downsamp_filt(:,t_var,:).^2,1));
mua_sse_residual = squeeze(nansum(mua_allcat_residual(:,t_var,:).^2,1));
mua_expl_var = (mua_sse_measured - mua_sse_residual)./mua_sse_measured;
figure;
plot((sum(mua_sse_measured,1)-sum(mua_sse_residual,1))./ ...
    sum(mua_sse_measured,1),'k','linewidth',2)
ylabel('Explained variance');
xlabel('Striatum depth');

if load_task
    
    % Get conditions for each trial, plot selected
    contrasts = [0,0.06,0.125,0.25,0.5,1];
    sides = [-1,1];
    choices = [-1,1];
    
    contrast_side_col = colormap_BlueWhiteRed(5);
    contrast_side_col(6,:) = 0;
    contrast_side_val = unique(sort([-contrasts,contrasts]))';
    
    % contrast, side, choice
%     plot_conditions = ...
%         [contrasts,contrasts; ...
%         -ones(1,6),-1,ones(1,5); ...
%         ones(1,6),-ones(1,6)]';
%     plot_conditions = ...
%         [0.125,1,0.125,1; ...
%         -1,-1,1,1; ...
%         1,1,-1,-1]';
%     % plot_conditions = ...
    %     [contrasts(2:end),contrasts(2:end); ...
    %     ones(1,10); ...
    %     ones(1,5),-ones(1,5)]';
    plot_conditions = ...
        [contrasts(2:end),contrasts(2:end); ...
        -ones(1,5),ones(1,5); ...
        ones(1,5),-ones(1,5)]';
    % plot_conditions = ...
    %     [contrasts(2:end),contrasts(2:end); ...
    %     -ones(1,5),ones(1,5); ...
    %     ones(1,5),ones(1,5)]';
    % plot_conditions = ...
    %     [0,0; ...
    %     -1,-1; ...
    %     -1,1]';
    
    use_rxn = move_t > 0 & move_t < 0.5;
    
    [~,plot_id] = ismember( ...
        [trial_contrast_allcat,trial_side_allcat,trial_choice_allcat], ...
        plot_conditions,'rows');    
    
    % Plot striatum
    figure; hold on;
    for curr_plot = 1:3
        
        switch curr_plot
            case 1
                plot_data = mua_allcat_downsamp_filt;
                plot_title = 'Measured';
            case 2
                plot_data = mua_allcat_predicted_nlin;
                plot_title = 'Predicted';
            case 3
                plot_data = mua_allcat_downsamp_filt - mua_allcat_predicted_nlin;
                plot_title = 'Residual';
        end
        
        p_str(curr_plot) = subplot(1,3,curr_plot); hold on;
        for curr_plot_condition = 1:size(plot_conditions,1)
            
            curr_trials = plot_id == curr_plot_condition & use_rxn;
            curr_data = plot_data(curr_trials,:,:);
            
            % re-align to movement onset
            t_leeway = 0.5;
            leeway_samples = round(t_leeway*(sample_rate/downsample_factor));
            curr_move_idx = move_idx(curr_trials);
            for i = 1:size(curr_data,1)
                curr_data(i,:) = circshift(curr_data(i,:),-curr_move_idx(i)+leeway_samples,2);
            end
            
            curr_data_mean = squeeze(nanmean(curr_data,1));
            
            curr_contrast_side = max(trial_contrast_allcat(curr_trials))*sign(max(trial_side_allcat(curr_trials)));
            contrast_side_idx = find(curr_contrast_side == contrast_side_val);
            curr_col = contrast_side_col(contrast_side_idx,:);
            
            if curr_contrast_side == 0
                switch max(trial_choice_allcat(curr_trials))
                    case -1
                        curr_col = 'm';
                    case 1
                        curr_col = 'c';
                end
            end
            
            AP_stackplot(curr_data_mean,t_downsample,7,false,curr_col,1:n_depths);
            
        end
        line([0,0],ylim,'color','k');
        xlabel('Time from stim');
        title(plot_title);
    end
    linkaxes(p_str)
    
else
    
    plot_cols = lines(max(D_cat.stimulus));
    
    % Plot striatum
    figure; hold on;
    for curr_plot = 1:3
        
        switch curr_plot
            case 1
                plot_data = mua_allcat_downsamp_filt;
                plot_title = 'Measured';
            case 2
                plot_data = mua_allcat_predicted_nlin;
                plot_title = 'Predicted';
            case 3
                plot_data = mua_allcat_downsamp_filt - mua_allcat_predicted_nlin;
                plot_title = 'Residual';
        end
        
        p_str(curr_plot) = subplot(1,3,curr_plot); hold on;
        for curr_plot_condition = 1:max(D_cat.stimulus)
            
            curr_trials = D_cat.stimulus == curr_plot_condition;
            curr_data = plot_data(curr_trials,:,:);
       
            curr_data_mean = squeeze(nanmean(curr_data,1));            

            curr_col = plot_cols(curr_plot_condition,:);

            AP_stackplot(curr_data_mean,t_downsample,2,false,curr_col,1:n_depths);
            
        end
        line([0,0],ylim,'color','k');
        xlabel('Time from stim');
        title(plot_title);
    end
    linkaxes(p_str)
    
end


% Test for significant differences in residual?
% re-align to movement onset
mua_allcat_downsamp_filt_move = mua_allcat_downsamp_filt;
mua_allcat_predicted_nlin_move = mua_allcat_predicted_nlin;
t_leeway = 0.5;
leeway_samples = round(t_leeway*(sample_rate/downsample_factor));
for i = 1:size(mua_allcat_downsamp,1)
    mua_allcat_downsamp_filt_move(i,:,:) = circshift(mua_allcat_downsamp_filt_move(i,:,:),-move_idx(i)+leeway_samples,2);
    mua_allcat_predicted_nlin_move(i,:,:) = circshift(mua_allcat_predicted_nlin_move(i,:,:),-move_idx(i)+leeway_samples,2);
end
mua_allcat_residual_move = mua_allcat_downsamp_filt_move - mua_allcat_predicted_nlin_move;

use_residual_t = t_downsample > -0.1 & t_downsample < 0.1;
mean_residual_t = squeeze(sum(mua_allcat_residual_move(:,use_residual_t,:),2));

grp_conditions = [ ...
    1,1,-1; ...
    1,-1,-1; ...
    0,-1,-1; ...
    1,1,1; ...
    1,-1,1; ...   
    0,-1,1];
[~,grp_id] = ismember( ...
    [trial_contrast_allcat > 0,trial_side_allcat,trial_choice_allcat], ...
    grp_conditions,'rows');

h = figure;
for curr_depth = 1:n_depths
    figure(h);
    subplot(n_depths,1,curr_depth);
    distributionPlot(mean_residual_t(:,curr_depth),'groups',grp_id);
    
    use_trials = move_t > 0 & move_t < 0.5;
    nonan = ~isnan(mean_residual_t(:,curr_depth));
    [p,tbl,stats] = anova1(mean_residual_t(nonan & use_trials,curr_depth),grp_id(nonan & use_trials));
    [results, means] = multcompare(stats);
end

figure; hold on;
use_trials = move_t > 0 & move_t < 0.5;
for curr_depth = 1:n_depths
    p(curr_depth) = subplot(n_depths,1,curr_depth); hold on;
    for curr_grp = 1:size(grp_conditions,1)
        %To plot residual
        plot(t_downsample,nanmean(mua_allcat_residual_move( ...
            grp_id == curr_grp & use_trials,:,curr_depth).^2,1),'linewidth',2);
        
%         % To plot explained variance
%         curr_meas = mua_allcat_downsamp_filt_move( ...
%             grp_id == curr_grp & use_trials,:,curr_depth);        
%         curr_res = mua_allcat_downsamp_filt_move( ...
%             grp_id == curr_grp & use_trials,:,curr_depth) - ...
%             mua_allcat_predicted_nlin_move( ...
%             grp_id == curr_grp & use_trials,:,curr_depth);      
%         expl_var = (nansum(curr_meas.^2,1) - nansum(curr_res.^2))./(nansum(curr_meas.^2));
%         plot(t_downsample,expl_var,'linewidth',2);
        
    end
    xlabel('Time from move');
    ylabel('Sum absolute error');
    line([0,0],ylim,'color','k');
    legend({'Stim R Move L','Stim L Move L','Stim 0 Move L', ...
        'Stim R Move R','Stim L Move R','Stim 0 Move R'});
end
linkaxes(p,'x');

%%%% TRYING HERE: fitting line to error
use_trials = find(trial_contrast_allcat > 0 & trial_side_allcat == 1 & ...
    trial_choice_allcat == -1 & move_t > 0 & move_t < 0.5);

use_measured = cell2mat(arrayfun(@(x) ...
    squeeze(mua_allcat_downsamp_filt(x,move_idx(x),:)),use_trials,'uni',false)')';
use_predicted = cell2mat(arrayfun(@(x) ...
    squeeze(mua_allcat_predicted_nlin(x,move_idx(x),:)),use_trials,'uni',false)')';

figure;
for curr_depth = 1:n_depths
    subplot(1,n_depths,curr_depth); hold on;
    
    nonan = ~isnan(use_predicted(:,curr_depth));
    curr_fit = polyfit(use_predicted(nonan,curr_depth),use_measured(nonan,curr_depth),1);
    
    n_bins = 20;
    activity_bounds = linspace(-1,5,n_bins+1);
    
    xlim([-1,5]);ylim([-1,5]);
    plot(use_predicted(:,curr_depth),use_measured(:,curr_depth),'.');
    line([-1,5],[-1,5],'color','k');
    line(xlim,xlim*curr_fit(1)+curr_fit(2),'color','r','linewidth',2);
    axis square;
    xlabel('Predicted');
    ylabel('Measured');
end




% another way - fit for each timepoint independently
use_trials = find( ...
    (trial_contrast_allcat > 0 & trial_side_allcat == 1 & ...
    trial_choice_allcat == -1) | ...
    (trial_contrast_allcat == 0 & trial_choice_allcat == -1) & ...
    move_t > 0 & move_t < 0.5);

% fit in loop - probably dumb way to do it but whatever
pred_fit = nan(length(t_downsample),2,n_depths);
for curr_depth = 1:n_depths
    for curr_t = 1:length(t_downsample)
               
        use_measured = mua_allcat_downsamp_filt_move(use_trials,curr_t,curr_depth);
        use_predicted = mua_allcat_predicted_nlin_move(use_trials,curr_t,curr_depth);
        
        nonan = ~isnan(use_predicted);
        pred_fit(curr_t,:,curr_depth) = ...
            polyfit(use_predicted(nonan),use_measured(nonan),1);                
    end
end

figure; 
subplot(1,2,1); hold on; set(gca,'ColorOrder',copper(n_depths));
plot(t_downsample,squeeze(pred_fit(:,1,:)),'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Multiplicative');
xlabel('Time from movement');
subplot(1,2,2); hold on; set(gca,'ColorOrder',copper(n_depths));
plot(t_downsample,squeeze(pred_fit(:,2,:)),'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Additive');
xlabel('Time from movement');



a = mua_allcat_predicted_nlin_move;
a(use_trials,:,:) = a(use_trials,:,:) .* permute(pred_fit(:,1,:),[2,1,3]);
a(use_trials,:,:) = a(use_trials,:,:) + permute(pred_fit(:,2,:),[2,1,3]);


% Plot striatum
figure; hold on;
for curr_plot = 1:3
    
    switch curr_plot
        case 1
            plot_data = mua_allcat_downsamp_filt_move;
            plot_title = 'Measured';
        case 2
            plot_data = a;
            plot_title = 'Predicted';
        case 3
            plot_data = mua_allcat_downsamp_filt_move - a;
            plot_title = 'Residual';
    end
    
    p_str(curr_plot) = subplot(1,3,curr_plot); hold on;
    for curr_plot_condition = 1:size(plot_conditions,1)
        
        curr_trials = plot_id == curr_plot_condition & use_rxn;
        curr_data = plot_data(curr_trials,:,:);
        
%         % re-align to movement onset
%         t_leeway = 0.5;
%         leeway_samples = round(t_leeway*(sample_rate/downsample_factor));
%         curr_move_idx = move_idx(curr_trials);
%         for i = 1:size(curr_data,1)
%             curr_data(i,:) = circshift(curr_data(i,:),-curr_move_idx(i)+leeway_samples,2);
%         end
        
        curr_data_mean = squeeze(nanmean(curr_data,1));
        
        curr_contrast_side = max(trial_contrast_allcat(curr_trials))*sign(max(trial_side_allcat(curr_trials)));
        contrast_side_idx = find(curr_contrast_side == contrast_side_val);
        curr_col = contrast_side_col(contrast_side_idx,:);
        
        if curr_contrast_side == 0
            switch max(trial_choice_allcat(curr_trials))
                case -1
                    curr_col = 'm';
                case 1
                    curr_col = 'c';
            end
        end
        
        AP_stackplot(curr_data_mean,t_downsample,2,false,curr_col,1:n_depths);
        
    end
    line([0,0],ylim,'color','k');
    xlabel('Time from stim');
    title(plot_title);
end
linkaxes(p_str)


%% Logistic regression on allcat activity

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_Udf_kernel-str_4_depths_long'];

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);
reward_all = cellfun(@(data,use_expts) data(use_expts),reward_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);
reward_all = reward_all(use_animals);

% Get widefield ROIs
n_rois = 200;

% MUA depths
n_depths = size(mua_all{1}{1},3);

% Set up conditions
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_rois,1);
mua_allcat = nan(0,length(t),n_depths,1);
wheel_allcat = nan(0,length(t),1);
reward_allcat = nan(0,length(t),1);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

D_cat = struct('stimulus',cell(1,1),'response',cell(1,1),'repeatNum',cell(1,1));

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    n_align = size(fluor_cat,4);
    
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x,[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
%     wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
    % Concatenate behavioural data
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
    D.response = cell2mat(cellfun(@(x) x.response,D_all{curr_animal},'uni',false));
    D.repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all{curr_animal},'uni',false));
    
    D.day = trial_day;
    
    D_cat.stimulus = vertcat(D_cat.stimulus,D.stimulus);
    D_cat.response = vertcat(D_cat.response,D.response);
    D_cat.repeatNum = vertcat(D_cat.repeatNum,D.repeatNum);
    
    % Get trial ID   
    trial_contrast = max(D.stimulus,[],2);
    [~,side_idx] = max(D.stimulus > 0,[],2);
    trial_side = (side_idx-1.5)*2;
    trial_choice = -(D.response-1.5)*2;
    
    trial_conditions = ...
        [trial_contrast, trial_side, ...
        trial_choice, ones(size(trial_day))];
    [~,trial_id] = ismember(trial_conditions,conditions,'rows');
    
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    reward_allcat = [reward_allcat;vertcat(reward_all{curr_animal}{:})];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];
   
end

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Get maximum wheel velocity in chosen direction
wheel_velocity_allcat = [zeros(size(wheel_allcat,1),1,1),diff(wheel_allcat,[],2)];
[max_speed,max_vel_idx] = max(abs(wheel_velocity_allcat(:,:,1).* ...
    (bsxfun(@times,wheel_velocity_allcat(:,:,1),trial_choice_allcat) > 0)),[],2);
vel_t = t(max_vel_idx);

max_vel = max_speed.*trial_choice_allcat;

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Low-pass filter fluorescence (where does this ~10 Hz crap come from?)
lowpassCutoff = 10; % Hz
[b100s, a100s] = butter(2, lowpassCutoff/((sample_rate)/2), 'low');
fluor_allcat_filt = filter(b100s,a100s,fluor_allcat,[],2);
mua_allcat_filt = filter(b100s,a100s,mua_allcat,[],2);

% re-align to movement onset
fluor_allcat_filt_move = fluor_allcat_filt;
mua_allcat_filt_move = mua_allcat_filt;
t_leeway = 0.5;
leeway_samples = round(t_leeway*(sample_rate));
for i = 1:size(fluor_allcat,1)
    fluor_allcat_filt_move(i,:,:) = circshift(fluor_allcat_filt_move(i,:,:),-move_idx(i)+leeway_samples,2);
    mua_allcat_filt_move(i,:,:) = circshift(mua_allcat_filt_move(i,:,:),-move_idx(i)+leeway_samples,2);    
end


%%% Modeling
use_trials = trial_contrast_allcat == 0.06 & move_t < 0.5;
cvfold = 10;

D_use = structfun(@(x) x(use_trials,:),D_cat,'uni',false);

warning off;
% Fit stim all
use_model = 'AP_test_stim';
g_stim_all = GLM(D_use).setModel(use_model).fit;
behavParameterFit = g_stim_all.parameterFits;

D_use.offset_ZL = g_stim_all.ZL(behavParameterFit, g_stim_all.Zinput(g_stim_all.data));

% Fit stim cross-validated
use_model = 'AP_test_stim';


g_stim = GLM(D_use).setModel(use_model).fitCV(cvfold);
pL = g_stim.p_hat(:,1);    pR = g_stim.p_hat(:,2);
likelihood = pL.*(g_stim.data.response == 1) + pR.*(g_stim.data.response == 2);

loglik_bpt_stim = nanmean(log2(likelihood));


% Fit stim + activity (all time)
use_model = 'AP_test_neur_stimoffset';

loglik_bpt_mua = nan(length(t),n_depths);
for curr_area = 1:n_depths
    for curr_t = 1:length(t)
        
        % Set the activity
        D_use.neur = mua_allcat_filt_move(use_trials,curr_t,curr_area);
        neur_nonan = ~isnan(D_use.neur);
        
        % Pick subset of trials
        D_curr = structfun(@(x) x(neur_nonan,:),D_use,'uni',false);
        
        clear g_act
        g_act = GLM(D_curr).setModel(use_model).fitCV(cvfold);
        pL = g_act.p_hat(:,1);
        pR = g_act.p_hat(:,2);
        likelihood = pL.*(g_act.data.response == 1) + pR.*(g_act.data.response == 2);
        
        frac_correct = sum((pL > 0.5) & g_act.data.response==1)/length(pL);
        
        loglik_bpt_mua(curr_t,curr_area) = nanmean(log2(likelihood));
        
        AP_print_progress_fraction(curr_t,length(t));
    end
    disp(curr_area);
end

loglik_increase_mua = loglik_bpt_mua - loglik_bpt_stim;

figure; hold on;
set(gca,'ColorOrder',copper(n_depths));
plot(t,loglik_increase_mua,'linewidth',2);
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');
xlabel('Time from movement');
ylabel('Loglikelihood increase from stim (bpt)');
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false))







% Logistic regression: fluorescence ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);
n_rois = numel(wf_roi);

roi_mask = cat(3,wf_roi.mask);

U_roi = bsxfun(@rdivide,transpose(reshape(U_master,[],size(U_master,3))'* ...
    reshape(roi_mask,[],size(roi_mask,3))), ...
    sum(reshape(roi_mask,[],size(roi_mask,3)),1)');

n_vs = size(fluor_allcat,3);
fluor_roi = permute(reshape(transpose(U_roi(:,1:n_vs)* ...
    reshape(permute(fluor_allcat(:,:,:,1),[3,2,1,4]),n_vs,[])), ...
    size(fluor_allcat,2),size(fluor_allcat,1),n_rois),[2,1,3]);

% Smooth > derivative > non-negative > normalize
smooth_size = 11;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

fluor_roi_diff = diff(padarray(convn(fluor_roi,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both'),[],2);
fluor_roi_diff(fluor_roi_diff < 0) = 0;
fluor_roi_diff = bsxfun(@rdivide,fluor_roi_diff,nanstd(reshape(fluor_roi_diff,[],1,n_rois),[],1));

% %%%%%%% TEMP: use straight fluorescence (subtract baseline)
% fluor_roi_diff = bsxfun(@minus,fluor_roi(:,1:end-1,:),nanmedian(fluor_roi(:,t < 0,:),2));
% %%%%%%%%%

%%%%%%% TEMP: use L-R derivative
load(wf_roi_fn);
n_rois = numel(wf_roi);

roi_mask = cat(3,wf_roi.mask);

U_roi = bsxfun(@rdivide,transpose(reshape(U_master,[],size(U_master,3))'* ...
    reshape(roi_mask,[],size(roi_mask,3))), ...
    sum(reshape(roi_mask,[],size(roi_mask,3)),1)');

n_vs = size(fluor_allcat,3);
fluor_roi = permute(reshape(transpose(U_roi(:,1:n_vs)* ...
    reshape(permute(fluor_allcat(:,:,:,1),[3,2,1,4]),n_vs,[])), ...
    size(fluor_allcat,2),size(fluor_allcat,1),n_rois),[2,1,3]);

n_rois = size(wf_roi,1);
fluor_roi_hemidiff = fluor_roi(:,:,1:n_rois) - fluor_roi(:,:,n_rois+1:end);
fluor_roi_hemidiff = bsxfun(@minus,fluor_roi_hemidiff(:,1:end-1,:),nanmedian(fluor_roi_hemidiff(:,t < 0,:),2));
% fluor_roi_diff = fluor_roi_hemidiff;

% Smooth > derivative > non-negative > normalize
smooth_size = 11;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

fluor_roi_diff = diff(padarray(convn(fluor_roi,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both'),[],2);
fluor_roi_diff(fluor_roi_diff < 0) = 0;
fluor_roi_diff = fluor_roi_diff(:,:,1:n_rois) - fluor_roi_diff(:,:,n_rois+1:end);
%%%%%%%%%

t_diff =  conv(t,[1,1]/2,'valid');

% Filter/align to movement
fluor_roi_diff_filt = filter(b100s,a100s,fluor_roi_diff,[],2);
fluor_roi_diff_filt_move = fluor_roi_diff_filt;
t_leeway = 0.5;
leeway_samples = round(t_leeway*(sample_rate));
for i = 1:size(fluor_allcat,1)
    fluor_roi_diff_filt_move(i,:,:) = circshift(fluor_roi_diff_filt_move(i,:,:),-move_idx(i)+leeway_samples,2);
end

% Fit stim + activity (all time)
use_model = 'AP_test_neur_stimoffset';

fluor_nonan_trials = ~squeeze(any(any(isnan(fluor_allcat),2),4));
loglik_bpt_fluor_roi = nan(length(t_diff),n_rois);
for curr_area = 1:n_rois
    for curr_t = 1:length(t_diff)
        
       % Set the activity
        D_use.neur = double(fluor_roi_diff_filt_move(use_trials,curr_t,curr_area));
        neur_nonan = ~isnan(D_use.neur);
        
        % Pick subset of trials
        D_curr = structfun(@(x) x(neur_nonan,:),D_use,'uni',false);
        
        clear g_act
        g_act = GLM(D_curr).setModel(use_model).fitCV(cvfold);
        pL = g_act.p_hat(:,1);
        pR = g_act.p_hat(:,2);
        likelihood = pL.*(g_act.data.response==1) + pR.*(g_act.data.response==2);
        
        loglik_bpt_fluor_roi(curr_t,curr_area) = nanmean(log2(likelihood));
        
        AP_print_progress_fraction(curr_t,length(t));
    end
    disp(curr_area);
end

loglik_increase_fluor_roi = loglik_bpt_fluor_roi - loglik_bpt_stim;

figure; hold on;
set(gca,'ColorOrder',jet(n_rois));
plot(t_diff,loglik_increase_fluor_roi,'linewidth',2);
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');
xlabel('Time from movement');
ylabel('Loglikelihood increase from stim (bpt)');
legend({wf_roi.area})


%% Plot activity binned by velocity

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_Udf_kernel-str_4_depths_long'];
% data_fn = 'all_trial_activity_df_msn_earlymove.mat';

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);
reward_all = cellfun(@(data,use_expts) data(use_expts),reward_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);
reward_all = reward_all(use_animals);

% Get widefield components
n_vs = 200;

% MUA depths
n_depths = size(mua_all{1}{1},3);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_vs,1);
mua_allcat = nan(0,length(t),n_depths,1);
wheel_allcat = nan(0,length(t),1);
reward_allcat = nan(0,length(t),1);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    n_align = size(fluor_cat,4);
    
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x,[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
%     wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
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
 
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    reward_allcat = [reward_allcat;vertcat(reward_all{curr_animal}{:})];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];
   
end

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Get maximum wheel velocity in chosen direction
wheel_velocity_allcat = [zeros(size(wheel_allcat,1),1,1),diff(wheel_allcat,[],2)];
[max_speed,max_vel_idx] = max(abs(wheel_velocity_allcat(:,:,1).* ...
    (bsxfun(@times,wheel_velocity_allcat(:,:,1),trial_choice_allcat) > 0)),[],2);
vel_t = t(max_vel_idx);

max_vel = max_speed.*trial_choice_allcat;

%%%%%%%%%%%%%% TEMP: max velocity regardless of final choice
max_speed_t = t > 0 & t < 1;
max_vel = AP_signed_max(wheel_velocity_allcat(:,max_speed_t),2);
% trial_choice_allcat = sign(max_vel);
%%%%%%%%%%%%%%

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Low-pass filter data
lowpassCutoff = 10; % Hz
[b100s, a100s] = butter(2, lowpassCutoff/((sample_rate)/2), 'low');
fluor_allcat_filt = filter(b100s,a100s,fluor_allcat,[],2);
mua_allcat_filt = filter(b100s,a100s,mua_allcat,[],2);

% Get fluorescence in widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi;
n_rois = numel(wf_roi);

roi_mask = cat(3,wf_roi.mask);

U_roi = bsxfun(@rdivide,transpose(reshape(U_master,[],size(U_master,3))'* ...
    reshape(roi_mask,[],size(roi_mask,3))), ...
    sum(reshape(roi_mask,[],size(roi_mask,3)),1)');

fluor_roi = permute(reshape(transpose(U_roi(:,1:n_vs)* ...
    reshape(permute(fluor_allcat_filt(:,:,:,1),[3,2,1,4]),n_vs,[])), ...
    size(fluor_allcat_filt,2),size(fluor_allcat_filt,1),n_rois),[2,1,3]);
fluor_roi_diff = [zeros(size(fluor_roi,1),1,n_rois),diff(fluor_roi,[],2)];

% re-align to movement onset
fluor_roi_diff_move = fluor_roi_diff;
mua_allcat_filt_move = mua_allcat_filt;
wheel_velocity_allcat_move = wheel_velocity_allcat;
t_leeway = 0.5;
leeway_samples = round(t_leeway*(sample_rate));
for i = 1:size(fluor_allcat,1)
    fluor_roi_diff_move(i,:,:) = circshift(fluor_roi_diff_move(i,:,:),-move_idx(i)+leeway_samples,2);
    mua_allcat_filt_move(i,:,:) = circshift(mua_allcat_filt_move(i,:,:),-move_idx(i)+leeway_samples,2);    
    wheel_velocity_allcat_move(i,:,:) = circshift(wheel_velocity_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);    
end

fluor_roi_diff_move = fluor_roi_diff_move(:,:,1:n_rois/2) - ...
    fluor_roi_diff_move(:,:,n_rois/2+1:end);

% Bin maximum velocity
use_trials = move_t > 0 & move_t < 1 & ...
    trial_side_allcat == trial_choice_allcat & trial_contrast_allcat > 0;

% n_prctiles = 9;
% use_vel = prctile(max_vel,100);
% vel_amp_edges = linspace(-use_vel,use_vel,n_prctiles*2);
n_prctiles = 6;
vel_amp_edges = prctile(abs(max_vel),linspace(0,100,n_prctiles));
vel_amp_edges = sort([vel_amp_edges,-vel_amp_edges]);

vel_amp_bins = discretize(max_vel(use_trials),vel_amp_edges);
% don't use the middle bin
vel_amp_bins(vel_amp_bins == n_prctiles) = NaN;

% Plot wheel and activity within bins
[grps_used,max_vel_grp_mean] = grpstats(max_vel(use_trials),vel_amp_bins,{'gname','nanmean'});
grps_used = cellfun(@str2num,grps_used);

use_wheel = wheel_velocity_allcat_move(use_trials,:);
wheel_grp_mean = grpstats(use_wheel,vel_amp_bins);
use_act = mua_allcat_filt_move(use_trials,:,2);
% use_act = fluor_roi_diff_move(use_trials,:,5);
act_grp_mean = grpstats(use_act,vel_amp_bins);

figure; 
col = [0.1*ones(n_prctiles,1),0.1*ones(n_prctiles,1),linspace(1,0.4,n_prctiles)'; ...
    linspace(0.4,1,n_prctiles)',0.1*ones(n_prctiles,1),0.1*ones(n_prctiles,1)];
col_used = col(grps_used,:);

subplot(2,2,1); hold on
set(gca,'ColorOrder',col_used)
plot(t,wheel_grp_mean','linewidth',2)
xlabel('Time from move');
ylabel('Wheel velocity');

subplot(2,2,2); hold on
set(gca,'ColorOrder',col_used)
plot(t,act_grp_mean','linewidth',2)
xlabel('Time from move');
ylabel('Activity');

plot_t = find(t > 0.05,1);
subplot(2,2,3); hold on;
plot(max_vel_grp_mean,max(act_grp_mean,[],2),'k','linewidth',2);
scatter(max_vel_grp_mean,max(act_grp_mean,[],2),80,col_used,'Filled');
xlabel('Wheel velocity');
ylabel('Activity');

plot_t = t < 0.5;
subplot(2,2,4); hold on;
set(gca,'ColorOrder',col_used)
plot(wheel_grp_mean(:,plot_t)',act_grp_mean(:,plot_t)','linewidth',2)
xlim([-max(abs(xlim)),max(abs(xlim))]);
xlabel('Wheel speed');
ylabel('Activity');

% Plot 2 activities by velocity
plot_t = t < 0.1;

use_act_1 = fluor_roi_diff_move(use_trials,:,5);
act_grp_mean_1 = grpstats(use_act_1,vel_amp_bins);

use_act_2 = mua_allcat_filt_move(use_trials,:,2);
act_grp_mean_2 = grpstats(use_act_2,vel_amp_bins);

figure; hold on;
set(gca,'ColorOrder',col_used)
plot(act_grp_mean_1(:,plot_t)',act_grp_mean_2(:,plot_t)','linewidth',2)
xlabel('Activity 1');
ylabel('Activity 2');


%% Predict choice by SVM

n_aligned_depths = 4;

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_Udf_kernel-str_4_depths_long'];

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);
reward_all = cellfun(@(data,use_expts) data(use_expts),reward_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);
reward_all = reward_all(use_animals);

% Get widefield ROIs
n_vs = 200;

% MUA depths
n_depths = size(mua_all{1}{1},3);

% Set up conditions
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_vs,1);
mua_allcat = nan(0,length(t),n_depths,1);
wheel_allcat = nan(0,length(t),1);
reward_allcat = nan(0,length(t),1);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

D_cat = struct('stimulus',cell(1,1),'response',cell(1,1),'repeatNum',cell(1,1));

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    n_align = size(fluor_cat,4);
    
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x,[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    mua_day_max = cell2mat(cellfun(@(x) ...
        repmat(permute(prctile(reshape(permute(x,[1,2,4,3]),[],n_depths),80,1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_max+softnorm);   
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
%     wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
    % Concatenate behavioural data
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
    D.response = cell2mat(cellfun(@(x) x.response,D_all{curr_animal},'uni',false));
    D.repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all{curr_animal},'uni',false));
    
    D.day = trial_day;
    
    D_cat.stimulus = vertcat(D_cat.stimulus,D.stimulus);
    D_cat.response = vertcat(D_cat.response,D.response);
    D_cat.repeatNum = vertcat(D_cat.repeatNum,D.repeatNum);
    
    % Get trial ID   
    trial_contrast = max(D.stimulus,[],2);
    [~,side_idx] = max(D.stimulus > 0,[],2);
    trial_side = (side_idx-1.5)*2;
    trial_choice = -(D.response-1.5)*2;
    
    trial_conditions = ...
        [trial_contrast, trial_side, ...
        trial_choice, ones(size(trial_day))];
    [~,trial_id] = ismember(trial_conditions,conditions,'rows');
    
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    reward_allcat = [reward_allcat;vertcat(reward_all{curr_animal}{:})];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];
   
end

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Get maximum wheel velocity in chosen direction
wheel_velocity_allcat = [zeros(size(wheel_allcat,1),1,1),diff(wheel_allcat,[],2)];
[max_speed,max_vel_idx] = max(abs(wheel_velocity_allcat(:,:,1).* ...
    (bsxfun(@times,wheel_velocity_allcat(:,:,1),trial_choice_allcat) > 0)),[],2);
vel_t = t(max_vel_idx);

max_vel = max_speed.*trial_choice_allcat;

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Low-pass filter fluorescence (where does this ~10 Hz crap come from?)
lowpassCutoff = 10; % Hz
[b100s, a100s] = butter(2, lowpassCutoff/((sample_rate)/2), 'low');
fluor_allcat_filt = filter(b100s,a100s,fluor_allcat,[],2);
mua_allcat_filt = filter(b100s,a100s,mua_allcat,[],2);

% re-align to movement onset
fluor_allcat_filt_move = fluor_allcat_filt;
mua_allcat_filt_move = mua_allcat_filt;
t_leeway = 0.5;
leeway_samples = round(t_leeway*(sample_rate));
for i = 1:size(fluor_allcat,1)
    fluor_allcat_filt_move(i,:,:) = circshift(fluor_allcat_filt_move(i,:,:),-move_idx(i)+leeway_samples,2);
    mua_allcat_filt_move(i,:,:) = circshift(mua_allcat_filt_move(i,:,:),-move_idx(i)+leeway_samples,2);    
end

% Get fluorescence in widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);
n_rois = numel(wf_roi);

roi_mask = cat(3,wf_roi.mask);

U_roi = bsxfun(@rdivide,transpose(reshape(U_master,[],size(U_master,3))'* ...
    reshape(roi_mask,[],size(roi_mask,3))), ...
    sum(reshape(roi_mask,[],size(roi_mask,3)),1)');

smooth_size = 1;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

fluor_roi = permute(reshape(transpose(U_roi(:,1:n_vs)* ...
    reshape(permute(fluor_allcat_filt_move(:,:,:,1),[3,2,1,4]),n_vs,[])), ...
    size(fluor_allcat_filt_move,2),size(fluor_allcat_filt_move,1),n_rois),[2,1,3]);

fluor_roi_diff = [zeros(size(fluor_roi,1),1,n_rois), ...
    diff(padarray(convn(fluor_roi,smWin,'valid'), ...
    [0,floor(size(smWin,2)/2)],'replicate','both'),[],2)];

% % hemidiff
% fluor_roi_diff(fluor_roi_diff < 0) = 0;
% fluor_roi_diff = fluor_roi_diff(:,:,1:size(wf_roi,1)) - ...
%     fluor_roi_diff(:,:,size(wf_roi,1)+1:end);

% (zero negatives, normalize)
fluor_roi_diff(fluor_roi_diff < 0) = 0;
fluor_roi_diff = bsxfun(@rdivide,fluor_roi_diff,nanstd(reshape(fluor_roi_diff,[],1,size(wf_roi,1)),[],1));

% SVM on each time point
activity_cat = cat(3,fluor_roi_diff,mua_allcat_filt_move);
frac_correct_predicted = nan(size(activity_cat,3),length(t));
for curr_area = 10%1:size(activity_cat,3)
    for curr_t = 1:length(t)        
        use_data = activity_cat(:,curr_t,curr_area);
        use_trials = ~isnan(use_data) & move_t > 0 & move_t < 0.5 & trial_contrast_allcat == 0;
        
        SVMModel = fitcsvm(use_data(use_trials),trial_choice_allcat(use_trials));
        CVSVMModel = crossval(SVMModel);
        predicted_choices = kfoldPredict(CVSVMModel);
        frac_correct_predicted(curr_area,curr_t) = ...
            sum(trial_choice_allcat(use_trials) == predicted_choices)/ ...
            sum(use_trials);   
        AP_print_progress_fraction(curr_t,length(t));
    end
end


















