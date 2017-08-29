%% Get behavior times

lick_idx = strcmp({Timeline.hw.inputs.name}, 'beamLickDetector');
lick_trace = Timeline.rawDAQData(:,lick_idx) > 2.5;
lick_times = Timeline.rawDAQTimestamps(find(lick_trace(2:end) & ~lick_trace(1:end-1))+1);
        
lick_bout_cutoff = 1; % seconds between what's considered a lick bout
lick_bout_times = lick_times([1,find(diff(lick_times) >= lick_bout_cutoff)+1]);

min_lick_to_reward_time = 1; % time from lick to reward to define rewarded lick
all_lick_reward_times = bsxfun(@minus,lick_bout_times',reward_t_timeline);
all_lick_reward_times(all_lick_reward_times > 0) = NaN;
lick_to_reward_time = max(all_lick_reward_times,[],2);
rewarded_licks = lick_to_reward_time >= -min_lick_to_reward_time;
if sum(rewarded_licks) ~= length(reward_t_timeline)
   error('Rewarded licks don''t match rewards') 
end

n_complete_stim = min(length(stimOn_t_timeline),length(stimOff_t_timeline));
stimOnTimes = stimOn_t_timeline(1:n_complete_stim);
stimOffTimes = stimOff_t_timeline(1:n_complete_stim);
azimuths = signals_events.trialAzimuthValues(1:n_complete_stim);

epoch_times = reshape([stimOnTimes;stimOffTimes],[],1);
epoch_hit = histcounts(signals_events.hitTimes,epoch_times) > 0;
stim_hit = epoch_hit(1:2:end);

stim_hit_licktime_cell = arrayfun(@(x) lick_times(find(...
    lick_times >= stimOnTimes(x) & ...
    lick_times(x),1)),1:length(stimOnTimes),'uni',false);
stim_hit_licktime = nan(size(stim_hit));
stim_hit_licktime(stim_hit) = [stim_hit_licktime_cell{stim_hit}];

%% Get visual striatum MUA PSTH around events

%use_spikes_idx = ismember(spike_templates,find(templateDepths >= 0 & templateDepths <= 1200));
use_spikes_idx = ismember(spike_templates,find(templateDepths > 3000 & templateDepths < 3200)) & ...
   (ismember(spike_templates,find(msn)));

use_spikes = spike_times_timeline(use_spikes_idx);

%align_times = stimOnTimes(stim_hit & azimuths == 0);
%align_times = stim_hit_licktime(stim_hit & azimuths == 0);
align_times = lick_bout_times(rewarded_licks);

raster_window = [-2,3.5];
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

%% Raster plot by depth aligned to stuff

%align_times = stimOn_t_timeline(~stim_hit & azimuths == 0);
%align_times = stim_hit_licktime(stim_hit & azimuths == 0);
align_times = lick_bout_times(rewarded_licks);

% Group by depth
n_depth_groups = 4;
depth_group_edges = linspace(850,3200,n_depth_groups+1);
depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
depth_group_edges(end) = Inf;
depth_group = discretize(spikeDepths,depth_group_edges);
depth_groups_used = unique(depth_group);

% Create MUA times grouped according to depth
mua_times = cell(n_depth_groups,1);
for curr_depth = 1:n_depth_groups
    use_spikes_idx = (depth_group == curr_depth) & ismember(spike_templates,find(fsi));
    mua_times{curr_depth} = spike_times_timeline(use_spikes_idx);
end

% PSTHs
raster_window = [-2.5,5];
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
figure; 
AP_stackplot(psth_smooth(:,20:end-20)',bins(20:end-20), ...
    10,true,'k',depth_group_centers);
line([0,0],ylim,'linestyle','--','color','k');
ylabel('Depth (\mum)');
xlabel('Time from stim onset (s)')
title('Population raster by depth');

%% PSTH viewer

align_times = stimOn_t_timeline;
%alignIDs = stim_hit(azimuths == 0)*1 + ~stim_hit(azimuths == 0)*2;
%align_times = stim_hit_licktime(stim_hit & azimuths == 0);
%align_times = lick_bout_times(rewarded_licks);

use_spikes_idx = ismember(spike_templates,find(templateDepths > 0 & templateDepths < 3200)) & ...
   (ismember(spike_templates,find(msn)));

raster_window = [-3.5,5];
psthViewer(spike_times_timeline(use_spikes_idx),spike_templates(use_spikes_idx), ...
    align_times,raster_window,alignIDs);

%% Average fluorescence around stuff

% Define the window to get an aligned response to
surround_window = [-2.5,5];

align_times = stimOn_t_timeline(~stim_hit & azimuths == 0);
%align_times = lick_bout_times(rewarded_licks);

% Get the surround time
framerate = 1./nanmedian(diff(frame_t));
surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);

% Don't use times that fall outside of imaging
align_times(align_times + surround_time(1) < frame_t(2) | ...
    align_times + surround_time(2) > frame_t(end)) = [];

% Use closest frames to times
align_surround_times = bsxfun(@plus, align_times', surround_time);
frame_edges = [frame_t,frame_t(end)+1/framerate];
align_frames = discretize(align_surround_times,frame_edges);

% If any aligned V's are NaNs (when does this happen?), don't use
align_frames(any(isnan(align_frames),2),:) = [];

aligned_V = reshape(fV(:,align_frames'), ...
    size(fV,1),size(align_frames,2),size(align_frames,1));

mean_aligned_V = nanmean(aligned_V,3);

% Get and plot the average fluorescence around event
mean_aligned_px = svdFrameReconstruct(U,mean_aligned_V);

AP_image_scroll(mean_aligned_px,surround_time);
warning off; truesize; warning on;

%% Regress stuff to fluorescence

% Skip the first n seconds to do this
skip_seconds = 20;
use_frames = (frame_t > skip_seconds) & (frame_t < (frame_t(end)-skip_seconds));

% Bin choiceworld events by frame
frame_edges = [frame_t,frame_t(end)+1/framerate];
signals_event_trace = [];

% Stim
unique_azimuths = unique(signals_events.trialAzimuthValues);
for trialAzimuth_idx = 1:length(unique_azimuths)        
        curr_azimuth = unique_azimuths(trialAzimuth_idx);       
        use_trials = signals_events.trialAzimuthValues == curr_azimuth;
        align_times = stimOn_t_timeline(use_trials)';
        signals_event_trace = [signals_event_trace;histcounts(align_times,frame_edges)];        
end

% Licks
frame_rewarded_licks_stim1 = histcounts(stim_hit_licktime(stim_hit & azimuths == 0),frame_edges);
frame_rewarded_licks_stim2 = histcounts(stim_hit_licktime(stim_hit & azimuths == 90),frame_edges);
frame_unrewarded_licks = histcounts(lick_bout_times(~rewarded_licks),frame_edges);
signals_event_trace = ...
    [signals_event_trace;frame_rewarded_licks_stim1;frame_rewarded_licks_stim2;frame_unrewarded_licks];

% % Rewards
% water_name = 'rewardEcho';
% water_idx = strcmp({Timeline.hw.inputs.name}, water_name);
% water_times = Timeline.rawDAQTimestamps(find(Timeline.rawDAQData(1:end-1,water_idx) <= 2 & ...
%     Timeline.rawDAQData(2:end,water_idx) > 2) + 1);
% frame_water = histcounts(water_times,frame_edges);
% signals_event_trace = [signals_event_trace;frame_water];

use_svs = 1:50;
kernel_frames = -35*2:35*6;
lambda = 0;
zs = false;
cvfold = 5;

[k,predicted_fluor,explained_var] = ...
    AP_regresskernel(signals_event_trace(:,use_frames), ...
    fV(use_svs,use_frames), ...
    kernel_frames,lambda,zs,cvfold);

% Reshape kernel and convert to pixel space
k_r = permute(reshape(k,size(signals_event_trace,1),length(kernel_frames),length(use_svs)),[3,2,1]);

r_px = zeros(size(U,1),size(U,2),size(k_r,2),size(k_r,3),'single');
for curr_event = 1:size(k_r,3);
    r_px(:,:,:,curr_event) = svdFrameReconstruct(U(:,:,use_svs),k_r(:,:,curr_event));
end

AP_image_scroll(r_px,kernel_frames/framerate);
caxis([prctile(r_px(:),[1,99])]*4);
truesize

%% Regress fluorescence to stuff

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = (frame_t > skip_seconds);

% Bin choiceworld events by frame
frame_edges = [frame_t,frame_t(end)+1/framerate];
signals_event_trace = [];

% Stim
azimuths = unique(signals_events.trialAzimuthValues);
for trialAzimuth_idx = 1:length(azimuths)        
        curr_azimuth = azimuths(trialAzimuth_idx);       
        use_trials = signals_events.trialAzimuthValues == curr_azimuth;
        align_times = stimOn_t_timeline(use_trials)';
        signals_event_trace = [signals_event_trace;histcounts(align_times,frame_edges)];        
end

% Licks
frame_rewarded_licks = histcounts(lick_bout_times(rewarded_licks),frame_edges);
frame_unrewarded_licks = histcounts(lick_bout_times(~rewarded_licks),frame_edges);
signals_event_trace = [signals_event_trace;frame_rewarded_licks;frame_unrewarded_licks];

% Rewards
water_name = 'rewardEcho';
water_idx = strcmp({Timeline.hw.inputs.name}, water_name);
water_times = Timeline.rawDAQTimestamps(find(Timeline.rawDAQData(1:end-1,water_idx) <= 2 & ...
    Timeline.rawDAQData(2:end,water_idx) > 2) + 1);
frame_water = histcounts(water_times,frame_edges);
signals_event_trace = [signals_event_trace;frame_water];

use_svs = 1:50;
kernel_frames = -35*1:7;
lambda = 1e6;
zs = false;
cvfold = 5;

[k,predicted_fluor,explained_var] = ...
    AP_regresskernel(fV(use_svs,use_frames), ...
    signals_event_trace(:,use_frames), ...
    kernel_frames,lambda,zs,cvfold);

% Reshape kernel and convert to pixel space
r = reshape(k,length(use_svs),length(kernel_frames),size(signals_event_trace,1));

r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3);
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(U(:,:,use_svs),r(:,:,curr_spikes));
end

AP_image_scroll(r_px,(kernel_frames)/framerate);
caxis([prctile(r_px(:),[1,99])]*4);
truesize


%% Predict lick bout times from fluorescence

% Skip the first n seconds to do this
skip_seconds = 0;
use_frames = (frame_t > skip_seconds);

% Bin choiceworld events by frame
frame_edges = [frame_t,frame_t(end)+1/framerate];
signals_event_trace = [];
   
% Licks
frame_licks = histcounts(lick_bout_times,frame_edges);
signals_event_trace = [signals_event_trace;frame_licks];

use_svs = 1:50;
kernel_frames = 0:15;
lambda = 1e7;
zs = false;
cvfold = 5;

[k,predicted_licks,explained_var] = ...
    AP_regresskernel(fV(use_svs,use_frames), ...
    signals_event_trace(:,use_frames), ...
    kernel_frames,lambda,zs,cvfold);

% Reshape kernel and convert to pixel space
r = reshape(k,length(use_svs),length(kernel_frames),size(signals_event_trace,1));

r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3);
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(U(:,:,use_svs),r(:,:,curr_spikes));
end

AP_image_scroll(r_px,(kernel_frames)/framerate);
caxis([prctile(r_px(:),[1,99])]*4);
truesize


%% Rasters/PSTHs around licks

lick_id = +rewarded_licks;
rewarded_lick_idx = find(rewarded_licks);
lick_id(rewarded_lick_idx(azimuths(stim_hit) == 0)) = 1;
lick_id(rewarded_lick_idx(azimuths(stim_hit) == 90)) = 2;

%use_spikes_idx = ismember(spike_templates,find(templateDepths >= 0 & templateDepths <= 1500));
use_spikes_idx = ismember(spike_templates,find(templateDepths > 0 & templateDepths < 1200)) & ...
   (ismember(spike_templates,find(msn)));

raster_window = [-1,1];
psthViewer(spike_times_timeline(use_spikes_idx),spike_templates(use_spikes_idx), ...
    lick_bout_times,raster_window,lick_id);




















