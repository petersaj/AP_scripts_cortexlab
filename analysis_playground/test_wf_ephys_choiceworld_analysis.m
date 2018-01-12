%% Batch load responses to passive stim

animal = 'AP026';
protocol = 'stimKalatsky';
experiments = AP_find_experiments(animal,protocol);

% only use experiments with ephys + imaging
experiments = experiments([experiments.imaging] & [experiments.ephys]);

load_parts.cam = false;
load_parts.imaging = true;
load_parts.ephys = false;

batch_vars = struct;
for curr_day = 1:length(experiments);
    
    day = experiments(curr_day).day;
    experiment = experiments(curr_day).experiment;
    
    AP_load_experiment
    
    %%%%%%%%%%%%%%%
    % DO THE STUFF
    %%%%%%%%%%%%%%%
    
    % Set options
    surround_window = [-0.2,3];
    baseline_surround_window = [0,0];
    framerate = 1./median(diff(frame_t));
    surround_samplerate = 1/(framerate*1);
    t_surround = surround_window(1):surround_samplerate:surround_window(2);
    baseline_surround_time = baseline_surround_window(1):surround_samplerate:baseline_surround_window(2);
    
    % Average (time course) responses
    conditions = unique(stimIDs);
    im_stim = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
    for curr_condition_idx = 1:length(conditions)
        curr_condition = conditions(curr_condition_idx);
        
        use_stims = find(stimIDs == curr_condition);
        use_stim_onsets = stimOn_times(use_stims(2:end));
        use_stim_onsets([1,end]) = [];
        
        stim_surround_times = bsxfun(@plus, use_stim_onsets(:), t_surround);
        peri_stim_v = permute(mean(interp1(frame_t,fVdf',stim_surround_times),1),[3,2,1]);
        
        im_stim(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
    end
    
    batch_vars.im_stim{curr_day} = im_stim;
    
    %%%%%%%%%%%%%%%%%%%%
    % THE STUFF IS DONE
    %%%%%%%%%%%%%%%%%%%%
    
    drawnow
    disp(curr_day)
    clearvars -except experiments curr_day animal batch_vars load_parts
    
end

disp('Finished batch.')

% Align images from batch processing
batch_vars_reg = batch_vars;

% % Get ddf
% for curr_day = 1:length(experiments)
%     curr_im = batch_vars_reg.im_stim{curr_day};
%     curr_im(isnan(curr_im)) = 0;
%     curr_im = imgaussfilt(diff(curr_im,[],3),2);
%     curr_im(curr_im < 0) = 0;
%     batch_vars_reg.im_stim{curr_day} = curr_im;
% end

% Align
days = {experiments.day};
[tform_matrix,im_aligned] = AP_align_widefield(animal,days);

for curr_day = 1:length(experiments);
    
    tform = affine2d;
    tform.T = tform_matrix{curr_day};
    
    curr_im = batch_vars_reg.im_stim{curr_day};
    curr_im(isnan(curr_im)) = 0;    

    batch_vars_reg.im_stim{curr_day} = imwarp(curr_im,tform, ...
        'Outputview',imref2d(size(im_aligned(:,:,1))));       

end

a = nanmean(cat(5,batch_vars_reg.im_stim{:}),5);

surround_window = [-0.2,3];
baseline_surround_window = [0,0];
framerate = 35;
surround_samplerate = 1/(framerate*1);
t_surround = surround_window(1):surround_samplerate:surround_window(2);
baseline_surround_time = baseline_surround_window(1):surround_samplerate:baseline_surround_window(2);

f = AP_image_scroll(a,t_surround);
axis image;

t_use = t_surround > 0.1 & t_surround < 0.5;
b = nanmean(a(:,:,t_use,:),3);
figure;imagesc(reshape(b,size(b,1),[],1)); 
colormap(gray); axis image off;
title([animal ': passive stimuli']);


%% Batch get average choiceworld fluorescence

animal = 'AP027';
protocol = 'vanillaChoiceworld';
experiments = AP_find_experiments(animal,protocol);

% only use experiments with ephys + imaging
experiments = experiments([experiments.imaging] & [experiments.ephys]);

load_parts.cam = false;
load_parts.imaging = true;
load_parts.ephys = false;

batch_vars = struct;
for curr_day = 1:length(experiments);
    
    day = experiments(curr_day).day;
    experiment = experiments(curr_day).experiment;
    
    AP_load_experiment
    
    %%%%%%%%%%%%%%%
    % DO THE STUFF
    %%%%%%%%%%%%%%%
    num_stim = min(length(signals_events.trialSideValues),length(stimOn_times));
    stimIDs = signals_events.trialSideValues(1:num_stim).*signals_events.trialContrastValues(1:num_stim);
    stim_onsets = stimOn_times(1:num_stim);
    
    % Discretize the stimIDs by easy/hard/zero
    stimIDs = discretize(stimIDs,[-Inf,-0.125,-0.01,0.01,0.25,Inf],[-2,-1,0,1,2]);
    
    %%%% Get wheel move time
    t_surround = [-0.5,5];
    surround_samples = t_surround/Timeline.hw.samplingInterval;
    
    % Get wheel aligned to stim onset
    rotaryEncoder_idx = strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder');
    t_surround = t_surround(1):Timeline.hw.samplingInterval:t_surround(2);
    pull_times = bsxfun(@plus,stim_onsets,t_surround);
    
    stim_aligned_wheel_raw = interp1(Timeline.rawDAQTimestamps, ...
        Timeline.rawDAQData(:,rotaryEncoder_idx),pull_times);
    stim_aligned_wheel = bsxfun(@minus,stim_aligned_wheel_raw, ...
        nanmedian(stim_aligned_wheel_raw(:,t_surround < 0),2));

    % Define time to first wheel movement
    thresh_displacement = 2;
    [~,wheel_move_sample] = max(abs(stim_aligned_wheel) > thresh_displacement,[],2);
    wheel_move_time = arrayfun(@(x) pull_times(x,wheel_move_sample(x)),1:size(pull_times,1));
    wheel_move_time(wheel_move_sample == 1) = NaN;
    %%%%
    
    
    % Set options
    surround_window = [-0.2,3];
    baseline_surround_window = [0,0];
    framerate = 1./median(diff(frame_t));
    surround_samplerate = 1/(framerate*1);
    t_surround = surround_window(1):surround_samplerate:surround_window(2);
    baseline_surround_time = baseline_surround_window(1):surround_samplerate:baseline_surround_window(2);
    
    % Average (time course) responses
    conditions = unique(stimIDs);
    im_stim_hit = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
    im_stim_miss = nan(size(U,1),size(U,2),length(t_surround),length(conditions));

    for curr_condition_idx = 1:length(conditions)
        curr_condition = conditions(curr_condition_idx);
        
        use_stims = find(stimIDs == curr_condition & signals_events.hitValues(1:num_stim) == 0);
        use_stim_onsets = stim_onsets(use_stims);       
        if length(use_stim_onsets) > 5           
            stim_surround_times = bsxfun(@plus, use_stim_onsets(:), t_surround);
            peri_stim_v = permute(mean(interp1(frame_t,fVdf',stim_surround_times),1),[3,2,1]);           
            im_stim_miss(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
        end
        
        use_stims = find(stimIDs == curr_condition & signals_events.hitValues(1:num_stim) == 1);
        use_stim_onsets = stim_onsets(use_stims);       
        if length(use_stim_onsets) > 5           
            stim_surround_times = bsxfun(@plus, use_stim_onsets(:), t_surround);
            peri_stim_v = permute(mean(interp1(frame_t,fVdf',stim_surround_times),1),[3,2,1]);           
            im_stim_hit(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
        end
        
    end
    
    batch_vars.im_stim_miss{curr_day} = im_stim_miss;
    batch_vars.im_stim_hit{curr_day} = im_stim_hit;

    
    %%%%%%%%%%%%%%%%%%%%
    % THE STUFF IS DONE
    %%%%%%%%%%%%%%%%%%%%
    
    drawnow
    disp(curr_day)
    clearvars -except experiments curr_day animal batch_vars load_parts
    
end

disp('Finished batch.')

% Align

days = {experiments.day};
[tform_matrix,im_aligned] = AP_align_widefield(animal,days);

batch_vars_reg = batch_vars;

% Spatially blur, get the > 0 delta, replace nans
for curr_day = 1:length(experiments)
    curr_im = batch_vars_reg.im_stim_miss{curr_day};
    nan_cond = squeeze(all(all(all(isnan(curr_im),1),2),3));
    curr_im(isnan(curr_im)) = 0;
    curr_im = imgaussfilt(diff(curr_im,[],3),1);
    curr_im(curr_im < 0) = 0;
    curr_im(:,:,:,nan_cond) = NaN;
    batch_vars_reg.im_stim_miss{curr_day} = curr_im;
    
    curr_im = batch_vars_reg.im_stim_hit{curr_day};
    nan_cond = squeeze(all(all(all(isnan(curr_im),1),2),3));
    curr_im(isnan(curr_im)) = 0;
    curr_im = imgaussfilt(diff(curr_im,[],3),1);
    curr_im(curr_im < 0) = 0;
    curr_im(:,:,:,nan_cond) = NaN;
    batch_vars_reg.im_stim_hit{curr_day} = curr_im;
end

% Align across days and replace nans
for curr_day = 1:length(experiments);
    
    tform = affine2d;
    tform.T = tform_matrix{curr_day};
    
    curr_im = batch_vars_reg.im_stim_miss{curr_day};
    nan_cond = squeeze(all(all(all(isnan(curr_im),1),2),3));
    curr_im(isnan(curr_im)) = 0;
    curr_im = imwarp(curr_im,tform, ...
        'Outputview',imref2d(size(im_aligned(:,:,1))));
    curr_im(:,:,:,nan_cond) = NaN;
    batch_vars_reg.im_stim_miss{curr_day} = curr_im;
    
    curr_im = batch_vars_reg.im_stim_hit{curr_day};
    nan_cond = squeeze(all(all(all(isnan(curr_im),1),2),3));
    curr_im(isnan(curr_im)) = 0;
    curr_im = imwarp(curr_im,tform, ...
        'Outputview',imref2d(size(im_aligned(:,:,1)))); 
    curr_im(:,:,:,nan_cond) = NaN;
    batch_vars_reg.im_stim_hit{curr_day} = curr_im;

end

% Mean
avg_hit = nanmean(cat(5,batch_vars_reg.im_stim_hit{:}),5);
avg_miss = nanmean(cat(5,batch_vars_reg.im_stim_miss{:}),5);

surround_window = [-0.2,3];
baseline_surround_window = [0,0];
framerate = 35;
surround_samplerate = 1/(framerate*1);
t_surround = surround_window(1):surround_samplerate:surround_window(2);
baseline_surround_time = baseline_surround_window(1):surround_samplerate:baseline_surround_window(2);

AP_image_scroll(avg_hit,t_surround)
AP_image_scroll(avg_miss,t_surround)


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

use_spikes_idx = ismember(spike_templates,find(templateDepths >= str_depth(2)-600 & templateDepths <= str_depth(2)-300));
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

depth_edges = [500,1500];
sample_rate_factor = 3;

% Define times to align
% SIGNALS - CHOICEWORLD
use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues > 0;% & ~isnan(wheel_move_time);
align_times = reshape(stimOn_times(use_trials(1:length(stimOn_times))),[],1);
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
[roi_trace,roi_mask] = AP_svd_roi(Udf,fVdf,weight_im,retinotopic_map); % weight_im, retinotopic_map, response_im
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
    curr_shuff
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

depth_edges = [500,1500];
sample_rate_factor = 1;

% SIGNALS - CHOICEWORLD
use_trials = signals_events.hitValues == 1;
% align_times = reshape(stimOn_times(use_trials),[],1);
align_times = reshape(signals_events.responseTimes(use_trials),[],1);
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
lambda = 0;
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

spikes_real_aligned = interp1(time_bin_centers,zscore(binned_spikes),t_peri_event);
spikes_pred_aligned = interp1(time_bin_centers,predicted_spikes,t_peri_event);

spikes_real_aligned_mean = grpstats(spikes_real_aligned,stim_conditions);
spikes_pred_aligned_mean = grpstats(spikes_pred_aligned,stim_conditions);

% Plot all responses
figure; hold on
p1 = AP_stackplot(spikes_real_aligned_mean',t_surround,3,false,'k',unique(stim_conditions));
p2 = AP_stackplot(spikes_pred_aligned_mean',t_surround,3,false,'r');
ylabel('Stim');
xlabel('Time from event onset');
legend([p1(1),p2(1)],{'Real','Predicted'});
line([0,0],ylim,'color','k');

% Plot responses by condition and error
response_t = [0,0.3];
response_t_use = t_surround >= response_t(1) & t_surround <= response_t(2);
spikes_real_response = nanmean(spikes_real_aligned(:,response_t),2);
spikes_pred_response = nanmean(spikes_pred_aligned(:,response_t),2);

subplot(1,2,1);


subplot(1,2,2);
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





