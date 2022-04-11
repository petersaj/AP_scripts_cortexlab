%% Plot lick task performance over days

animal = 'AP016';
% just get all days for now (eventually choose, supply date range, etc)
expInfo_path = ['\\zserver.cortexlab.net\Data\expInfo\' animal];
expInfo_dir = dir(expInfo_path);
days = {expInfo_dir(find([expInfo_dir(3:end).isdir])+2).name};

% Get lick task
lick_task_expts = nan(size(days));
for curr_day = 1:length(days)  
    day = days{curr_day};
    % In the event of multiple experiments, check all
    expDay_dir = dir([expInfo_path filesep days{curr_day}]);
    exp_nums = cellfun(@str2num,{expDay_dir(3:end).name});
    use_exp = false(size(exp_nums));
    for curr_exp = 1:length(exp_nums);
        [block_filename, block_exists] = AP_cortexlab_filename(animal,day,exp_nums(curr_exp),'block');
        if ~block_exists
            continue
        end
        % Load the block file
        load(block_filename)
        [~,expDef] = fileparts(block.expDef);
        use_exp(curr_exp) = ~isempty(strfind(expDef,'LickReward'));    
    end
    if any(use_exp)
        lick_task_expts(curr_day) = exp_nums(use_exp);
    end
end

% Initialize the behavior structure
days = days(~isnan(lick_task_expts));
expts = lick_task_expts(~isnan(lick_task_expts));
bhv = struct;

for curr_day = 1:length(days)
    day = days{curr_day};
    experiment = expts(curr_day);
    load_parts = struct;
    AP_load_experiment;

    lick_idx = strcmp({Timeline.hw.inputs.name}, 'beamLickDetector');
    lick_trace = Timeline.rawDAQData(:,lick_idx) > 2.5;
    lick_times = Timeline.rawDAQTimestamps(find(lick_trace(2:end) & ~lick_trace(1:end-1))+1);
    
    epoch_times = reshape([stimOn_times,stimOff_times]',[],1);
    epoch_hit = histcounts(lick_times,epoch_times) > 0;
    
    stim_hit = epoch_hit(1:2:end)';
          
    stim_hit_licktime_cell = arrayfun(@(x) lick_times(find(...
        lick_times >= stimOn_times(x),1)),1:length(stimOn_times),'uni',false);
    stim_hit_licktime = nan(size(stim_hit));
    stim_hit_licktime(stim_hit) = [stim_hit_licktime_cell{stim_hit}];
    stim_lick_delay = stim_hit_licktime - stimOn_times;
    
    % Get lick aligned to stim hit
    surround_interval = [-2,2];
    surround_time = surround_interval(1):Timeline.hw.samplingInterval:surround_interval(2);
    water_surround_lick = bsxfun(@plus,stim_hit_licktime(~isnan(stim_hit_licktime)), ...
        surround_time);   
    water_aligned_lick = interp1(Timeline.rawDAQTimestamps,Timeline.rawDAQData(:,lick_idx), ...
        water_surround_lick);
    
    bhv(curr_day).stim_hit = stim_hit;
    bhv(curr_day).stim_lick_delay = stim_lick_delay;
    bhv(curr_day).water_aligned_lick = water_aligned_lick;
    
end

figure;

subplot(1,3,1); hold on;
plotSpread({bhv.stim_lick_delay});
errorbar(cellfun(@nanmedian,{bhv.stim_lick_delay}), ...
    cellfun(@(x) nanstd(x)./sqrt(sum(~isnan(x))),{bhv.stim_lick_delay}),'r')
xlabel('Day');
ylabel('Time from stim to lick')

subplot(2,3,2); hold on; set(gca,'ColorOrder',copper(length(days)));
frac_stim = arrayfun(@(x) smooth(+bhv(x).stim_hit,20),1:length(days),'uni',false);
for i = 1:length(days)
    plot(frac_stim{i},'linewidth',2);
end
ylim([0 1]);
xlabel('Trial');
ylabel('Fraction hit trials');

subplot(2,3,5);
plot(cellfun(@median,frac_stim),'k','linewidth',2);
xlabel('Day');
ylabel('Median smoothed fraction hit');
ylim([0 1]);

subplot(1,3,3)
day_trials_plotted = arrayfun(@(x) size(bhv(x).water_aligned_lick,1),1:length(days));
imagesc(surround_time,1:sum(day_trials_plotted),1-vertcat(bhv.water_aligned_lick));
colormap(gray);
cumulative_trials = cumsum(day_trials_plotted);
for i = 1:length(days)
    line(xlim,[cumulative_trials(i),cumulative_trials(i)],'color','r');
end
ylabel('Trial');
xlabel('Time from rewarded lick');


%% Get behavior times

lick_idx = strcmp({Timeline.hw.inputs.name}, 'beamLickDetector');
lick_trace = Timeline.rawDAQData(:,lick_idx) > 2.5;
lick_times = Timeline.rawDAQTimestamps(find(lick_trace(2:end) & ~lick_trace(1:end-1))+1);
        
lick_bout_cutoff = 1; % seconds between what's considered a lick bout
lick_bout_starts = lick_times([1,find(diff(lick_times) >= lick_bout_cutoff)+1]);
lick_bout_stops = lick_times([find(diff(lick_times) >= lick_bout_cutoff),length(lick_times)]);

min_lick_to_reward_time = 1; % time from lick to reward to define rewarded lick
all_lick_reward_times = bsxfun(@minus,lick_bout_starts',reward_t_timeline);
all_lick_reward_times(all_lick_reward_times > 0) = NaN;
lick_to_reward_time = max(all_lick_reward_times,[],2);
rewarded_licks = lick_to_reward_time >= -min_lick_to_reward_time;
if sum(rewarded_licks) ~= length(reward_t_timeline)
   error('Rewarded licks don''t match rewards') 
end

azimuths = signals_events.trialAzimuthValues;

epoch_times = reshape([stimOn_times,stimOff_times]',[],1);
epoch_hit = histcounts(lick_times,epoch_times) > 0;
stim_hit = epoch_hit(1:2:end);

stim_hit_licktime_cell = arrayfun(@(x) lick_times(find(...
    lick_times >= stimOn_times(x) & ...
    lick_times(x),1)),1:length(stimOn_times),'uni',false);
stim_hit_licktime = nan(size(stim_hit));
stim_hit_licktime(stim_hit) = [stim_hit_licktime_cell{stim_hit}];


%% Get visual striatum MUA PSTH around events

%use_spikes_idx = ismember(spike_templates,find(template_depths >= 0 & template_depths <= 1200));
use_spikes_idx = ismember(spike_templates,find(template_depths > 2000 & template_depths < 3000)) & ...
   (ismember(spike_templates,find(msn)));

use_spikes = spike_times_timeline(use_spikes_idx);

align_times = stimOn_times(stim_hit & azimuths == 90);
%align_times = stim_hit_licktime(stim_hit & azimuths == 0);
%align_times = lick_bout_starts(rewarded_licks);

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

%align_times = stimOn_times(stim_hit & azimuths == 90);
%align_times = stim_hit_licktime(stim_hit & azimuths == 90);
align_times = lick_bout_starts(~rewarded_licks);

% Group by depth
n_depth_groups = 4;
depth_group_edges = linspace(1000,4000,n_depth_groups+1);
depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
depth_group_edges(end) = Inf;
depth_group = discretize(spike_depths,depth_group_edges);
depth_groups_used = unique(depth_group);

% Create MUA times grouped according to depth
mua_times = cell(n_depth_groups,1);
for curr_depth = 1:n_depth_groups
    use_spikes_idx = (depth_group == curr_depth) & ismember(spike_templates,find(msn));
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
    trace_spacing,false,'k',depth_group_centers);
line([0,0],ylim,'linestyle','--','color','k');
ylabel('Depth (\mum)');
xlabel('Time from stim onset (s)')
title('Population raster by depth');

%% PSTH viewer

align_times = stimOn_times;
alignIDs = azimuths;
%alignIDs = stim_hit(azimuths == 0)*1 + ~stim_hit(azimuths == 0)*2;
%align_times = stim_hit_licktime(stim_hit & azimuths == 0);
%align_times = lick_bout_times(rewarded_licks);

use_spikes_idx = ismember(spike_templates,find(template_depths > 0 & template_depths < 1200)) & ...
   (ismember(spike_templates,find(msn)));

raster_window = [-3.5,5];
psthViewer(spike_times_timeline(use_spikes_idx),spike_templates(use_spikes_idx), ...
    align_times,raster_window,alignIDs);

%% Average fluorescence around stuff

% Define the window to get an aligned response to
surround_window = [-2.5,5];

align_times = stimOn_times(stim_hit & azimuths == 0);
%align_times = lick_bout_times(rewarded_licks);

% Get the surround time
framerate = 1./nanmedian(diff(frame_t));
surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);

% Don't use times that fall outside of imaging
align_times(align_times + surround_time(1) < frame_t(2) | ...
    align_times + surround_time(2) > frame_t(end)) = [];

% Use closest frames to times
align_surround_times = bsxfun(@plus, align_times, surround_time);
frame_edges = [frame_t,frame_t(end)+1/framerate];
align_frames = discretize(align_surround_times,frame_edges);

% If any aligned V's are NaNs (when does this happen?), don't use
align_frames(any(isnan(align_frames),2),:) = [];

aligned_V = reshape(fV(:,align_frames'), ...
    size(fV,1),size(align_frames,2),size(align_frames,1));

mean_aligned_V = nanmean(aligned_V,3);

% Get and plot the average fluorescence around event
mean_aligned_px = svdFrameReconstruct(U,mean_aligned_V);

AP_imscroll(mean_aligned_px,surround_time);
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
        align_times = stimOn_times(use_trials)';
        signals_event_trace = [signals_event_trace;histcounts(align_times,frame_edges)];        
end

% Licks
frame_rewarded_licks_stim1 = histcounts(stim_hit_licktime(stim_hit & azimuths == 0),frame_edges);
frame_rewarded_licks_stim2 = histcounts(stim_hit_licktime(stim_hit & azimuths == 90),frame_edges);
frame_unrewarded_licks = histcounts(lick_bout_starts(~rewarded_licks),frame_edges);
signals_event_trace = ...
    [signals_event_trace;frame_rewarded_licks_stim1;frame_rewarded_licks_stim2;frame_unrewarded_licks];

% Rewards
water_name = 'rewardEcho';
water_idx = strcmp({Timeline.hw.inputs.name}, water_name);
water_times = Timeline.rawDAQTimestamps(find(Timeline.rawDAQData(1:end-1,water_idx) <= 2 & ...
    Timeline.rawDAQData(2:end,water_idx) > 2) + 1);
frame_water = histcounts(water_times,frame_edges);
signals_event_trace = [signals_event_trace;frame_water];

% Quiescence
lick_bout_edges = reshape([lick_bout_starts;lick_bout_stops],[],1);
[~,~,lick_bout_frame_bin] = histcounts(frame_t,lick_bout_edges);
lick_frames = histcounts(lick_times,frame_edges);
quiescence = mod(lick_bout_frame_bin,2) == 0 & ~lick_frames;
quiescence_leeway_time = 1; % give it some leeway
quiescence_leeway_frames = round(quiescence_leeway_time*framerate);
quiescence_leeway = conv(+quiescence,ones(1,quiescence_leeway_frames),'same') ...
    == quiescence_leeway_frames;
signals_event_trace = [signals_event_trace;quiescence_leeway];

use_svs = 1:50;
kernel_frames = -35*2:35*6;
lambda = 0;
zs = [false,false];
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

AP_imscroll(r_px,kernel_frames/framerate);
caxis([prctile(r_px(:),[1,99])]*4);
truesize

% Get map of explained variance
downsample_factor = 10;
spatial_explained_var = AP_spatial_explained_var(U(:,:,use_svs), ...
    fV(use_svs,use_frames),predicted_fluor,downsample_factor);
figure;imagesc(spatial_explained_var);
caxis([-max(abs(spatial_explained_var(:))),max(abs(spatial_explained_var(:)))]);
colormap(colormap_BlueWhiteRed);
colorbar;

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
        align_times = stimOn_times(use_trials)';
        signals_event_trace = [signals_event_trace;histcounts(align_times,frame_edges)];        
end

% Licks
frame_rewarded_licks = histcounts(lick_bout_starts(rewarded_licks),frame_edges);
frame_unrewarded_licks = histcounts(lick_bout_starts(~rewarded_licks),frame_edges);
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
zs = [false,false];
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

AP_imscroll(r_px,(kernel_frames)/framerate);
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
frame_licks = histcounts(lick_bout_starts,frame_edges);
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

AP_imscroll(r_px,(kernel_frames)/framerate);
caxis([prctile(r_px(:),[1,99])]*4);
truesize


%% Rasters/PSTHs around licks

lick_id = +rewarded_licks;
rewarded_lick_idx = find(rewarded_licks);
lick_id(rewarded_lick_idx(azimuths(stim_hit) == 0)) = 1;
lick_id(rewarded_lick_idx(azimuths(stim_hit) == 90)) = 2;

%use_spikes_idx = ismember(spike_templates,find(template_depths >= 0 & template_depths <= 1500));
use_spikes_idx = ismember(spike_templates,find(template_depths > 2500 & template_depths < 3000)) & ...
   (ismember(spike_templates,find(msn)));

raster_window = [-1,1];
psthViewer(spike_times_timeline(use_spikes_idx),spike_templates(use_spikes_idx), ...
    lick_bout_starts,raster_window,lick_id);


%% PSTH around rewarded lick

raster_window = [-5,5];
psth_bin_size = 0.001;

rewarded_stim1_lick_psth = nan(size(templates,1),diff(raster_window)/psth_bin_size);
rewarded_stim2_lick_psth = nan(size(templates,1),diff(raster_window)/psth_bin_size);
nonrewarded_lick_psth = nan(size(templates,1),diff(raster_window)/psth_bin_size);
for curr_template = unique(spike_templates)'
    
    rewarded_stim1_lick_psth(curr_template,:) = psthAndBA( ...
        spike_times_timeline(spike_templates == curr_template), ...
        stim_hit_licktime(stim_hit & azimuths == 0), ...
        raster_window, psth_bin_size);
    
    rewarded_stim2_lick_psth(curr_template,:) = psthAndBA( ...
        spike_times_timeline(spike_templates == curr_template), ...
        stim_hit_licktime(stim_hit & azimuths == 90), ...
        raster_window, psth_bin_size);
    
    nonrewarded_lick_psth(curr_template,:) = psthAndBA( ...
        spike_times_timeline(spike_templates == curr_template), ...
        lick_bout_starts(~rewarded_licks), ...
        raster_window, psth_bin_size);
    
    disp(curr_template);
    
end

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
rewarded_stim1_lick_psth_smooth = conv2(rewarded_stim1_lick_psth, smWin, 'same');
rewarded_stim2_lick_psth_smooth = conv2(rewarded_stim2_lick_psth, smWin, 'same');
nonrewarded_lick_psth_smooth = conv2(nonrewarded_lick_psth, smWin, 'same');

%[Row,Col,Scalar,Model,Residual] = MakeSeparable(psth_smooth');












