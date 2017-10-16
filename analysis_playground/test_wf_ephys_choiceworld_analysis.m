%% Batch template for loading
% this isn't done yet

animal = 'AP024';
protocol = 'vanillaChoiceworld';
experiments = AP_find_experiments(animal,protocol);

% only use experiments with ephys + imaging
experiments = experiments([experiments.imaging] & ~[experiments.ephys]);

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
    
    % Get left/right choice trials (of chosen contrasts)
    frame_edges = [frame_t,frame_t(end)+1/framerate];
    
    left_trials = ismember(signals_events.trialContrastValues,[1,0.5,0.25,0.0125,0.06]) &  ...
        ((ismember(signals_events.trialSideValues,[1]) & ismember(signals_events.hitValues,[1])) | ...
        (ismember(signals_events.trialSideValues,[-1]) & ismember(signals_events.hitValues,[0]))) & ...
        ~signals_events.repeatTrialValues;
    
    right_trials = ismember(signals_events.trialContrastValues,[1,0.5,0.25,0.0125,0.06]) &  ...
        ((ismember(signals_events.trialSideValues,[1]) & ismember(signals_events.hitValues,[0])) | ...
        (ismember(signals_events.trialSideValues,[-1]) & ismember(signals_events.hitValues,[1]))) & ...
        ~signals_events.repeatTrialValues;
        
    % Fix the parameters
    use_window = [-0.2,0.2];
    
    framerate = 1./median(diff(frame_t));
    surround_samplerate = 1/(framerate*1);
    surround_time = use_window(1):surround_samplerate: ...
        use_window(2);
    
    align_surround_times_left = bsxfun(@plus, stimOn_times(left_trials), surround_time);
    align_surround_times_right = bsxfun(@plus, stimOn_times(right_trials), surround_time);
    
    fV_align_left = interp1(frame_t,fVdf',align_surround_times_left);
    fV_align_right = interp1(frame_t,fVdf',align_surround_times_right);
    
    stim_order = [[ones(1,sum(left_trials)),zeros(1,sum(right_trials))]; ...
        [zeros(1,sum(left_trials)),ones(1,sum(right_trials))]];
    
    use_svs = 1:100;
    lambda = 1e6;
    zs = [false,false];
    cvfold = 5;
    
    fV_align_all = [reshape(permute(fV_align_left(:,:,use_svs),[2,3,1]),[],sum(left_trials)), ...
        reshape(permute(fV_align_right(:,:,use_svs),[2,3,1]),[],sum(right_trials))];
    
    [k,predicted_stim,explained_var] = ...
        AP_regresskernel(fV_align_all,stim_order,0,lambda,zs,cvfold);
    
    max_likely_stim = bsxfun(@eq,predicted_stim,max(predicted_stim,[],1));
    correct_decoding = sum(max_likely_stim(:) & stim_order(:))./sum(stim_order(:));
    
    k2 = reshape(k,[],length(use_svs),2);
    k_px_l = svdFrameReconstruct(Udf(:,:,use_svs),k2(:,:,1)');
    k_px_r = svdFrameReconstruct(Udf(:,:,use_svs),k2(:,:,2)');
    
    baseline_time = find(surround_time < 0,1,'last');
    k_px_l_norm = bsxfun(@minus,k_px_l,nanmean(k_px_l(:,:,1:baseline_time),3));
    k_px_r_norm = bsxfun(@minus,k_px_r,nanmean(k_px_r(:,:,1:baseline_time),3));
    
    batch_vars.k_px_l{curr_day} = k_px_l_norm;
    batch_vars.k_px_r{curr_day} = k_px_r_norm;
    
    %%%%%%%%%%%%%%%%%%%%
    % THE STUFF IS DONE
    %%%%%%%%%%%%%%%%%%%%
    
    drawnow
    disp(curr_day)
    clearvars -except experiments curr_day animal batch_vars load_parts
    
end

disp('Finished batch.')

%% Align images from batch processing

batch_vars_reg = batch_vars;

for curr_day = setdiff(1:length(experiments),ref_im_num);
    
    tform = affine2d;
    tform.T = tform_matrix{curr_day};
    
    curr_im = batch_vars_reg.k_px_l{curr_day};
    curr_im(isnan(curr_im)) = 0;    
    batch_vars_reg.k_px_l{curr_day} = imwarp(curr_im,tform, ...
        'Outputview',imref2d(size(avg_im{ref_im_num})));
    
    curr_im = batch_vars_reg.k_px_r{curr_day};
    curr_im(isnan(curr_im)) = 0;    
    batch_vars_reg.k_px_r{curr_day} = imwarp(curr_im,tform, ...
        'Outputview',imref2d(size(avg_im{ref_im_num})));

end

a = nanmean(cat(4,batch_vars_reg.k_px_l{:}),4);
b = nanmean(cat(4,batch_vars_reg.k_px_r{:}),4);
AP_image_scroll([a,b],surround_time)


%% PSTH to choiceworld conditions

stimIDs = signals_events.trialSideValues.*signals_events.trialContrastValues;

use_spikes_idx = ismember(spike_templates,find(templateDepths >= 1000 & templateDepths <= 2000));
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


%% PSTH for left vs. right stim, choose left vs. right stim (by depth)

% Group multiunit by depth
n_depth_groups = 6;
%depth_group_edges = linspace(0,max(channel_positions(:,2)),n_depth_groups+1);
%depth_group_edges = linspace(0,4000,n_depth_groups+1);
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
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
zs = false;
if zs
    trace_spacing = 10;
end
p_rh = AP_stackplot(psth_right_hit',bins,trace_spacing,zs,'k',depth_group_centers);
p_rm = AP_stackplot(psth_right_miss',bins,trace_spacing,zs,'r',depth_group_centers);
p_lh = AP_stackplot(psth_left_hit',bins,trace_spacing,zs,'b',depth_group_centers);
p_lm = AP_stackplot(psth_left_miss',bins,trace_spacing,zs,'m',depth_group_centers);

line([0,0],ylim,'color','k','linestyle','--');
line([0.2,0.2],ylim,'color','k','linestyle','--');

legend([p_rh(1),p_rm(1),p_lh(1),p_lm(1)],{'Stim right hit','Stim right miss','Stim left hit','Stim left miss'})
xlabel('Time from stim');
ylabel('Depth (\mum)');


%% PSTH for left vs. right stim, choose left vs. right stim (errorbars)

stimIDs = signals_events.trialSideValues.*signals_events.trialContrastValues;

use_spikes_idx = ismember(spike_templates,find(templateDepths >= 0 & templateDepths <= 1200));
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
p1 = AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'k');

use_trials = (signals_events.trialSideValues == -1 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == 1 & signals_events.hitValues == 0);
[psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(use_spikes,signals_events.stimOnTimes(use_trials),raster_window,psth_bin_size);
p2 = AP_errorfill(bins,nanmean(binnedArray),nanstd(binnedArray,[],1)./sqrt(sum(~isnan(binnedArray),1)),'r');

legend([p1(1),p2(1)],{'Move left','Move right'});
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


%% PSTH population choose left vs. right stim (move-aligned)

stimIDs = signals_events.trialSideValues.*signals_events.trialContrastValues;

use_spikes_idx = ismember(spike_templates,find(templateDepths >= 0 & templateDepths <= 1200));
%use_spikes_idx = ismember(spike_templates,intersect(find(templateDepths >= 0 & templateDepths <= 1200),find(msn)));

use_spikes = spike_times_timeline(use_spikes_idx);

raster_window = [-1 2];
psth_bin_size = 0.02;

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

% Get wheel movement time for each trial

surround_time = [-0.2,5];
surround_samples = surround_time/Timeline.hw.samplingInterval;

rotaryEncoder_idx = strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder');
surround_time = surround_time(1):Timeline.hw.samplingInterval:surround_time(2);
pull_times = bsxfun(@plus,signals_events.stimOnTimes',surround_time);

stim_aligned_wheel_raw = interp1(Timeline.rawDAQTimestamps, ...
    Timeline.rawDAQData(:,rotaryEncoder_idx),pull_times);
stim_aligned_wheel = bsxfun(@minus,stim_aligned_wheel_raw, ...
    nanmedian(stim_aligned_wheel_raw(:,surround_time < 0),2));

thresh_displacement = 10;
[~,wheel_move_sample] = max(abs(stim_aligned_wheel) > thresh_displacement,[],2);
wheel_move_time = arrayfun(@(x) pull_times(x,wheel_move_sample(x)),1:size(pull_times,1));
wheel_move_time(wheel_move_sample == 1) = NaN;

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

%use_spikes_idx = ismember(spike_templates,find(templateDepths >= 2000 & templateDepths <= 3000));
use_spikes_idx = ismember(spike_templates,intersect(find(templateDepths >= 2000 & templateDepths <= 3000),find(tan)));
use_spikes = spike_times_timeline(use_spikes_idx);

use_templates = unique(spike_templates(use_spikes_idx));

% Get wheel movement time for each trial
surround_time = [-1,2];
surround_samples = surround_time/Timeline.hw.samplingInterval;

rotaryEncoder_idx = strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder');
surround_time = surround_time(1):Timeline.hw.samplingInterval:surround_time(2);
pull_times = bsxfun(@plus,signals_events.stimOnTimes',surround_time);

stim_aligned_wheel_raw = interp1(Timeline.rawDAQTimestamps, ...
    Timeline.rawDAQData(:,rotaryEncoder_idx),pull_times);
stim_aligned_wheel = bsxfun(@minus,stim_aligned_wheel_raw, ...
    nanmedian(stim_aligned_wheel_raw(:,surround_time < 0),2));

thresh_displacement = 10;
[~,wheel_move_sample] = max(abs(stim_aligned_wheel) > thresh_displacement,[],2);
wheel_move_time = arrayfun(@(x) pull_times(x,wheel_move_sample(x)),1:size(pull_times,1));
wheel_move_time(wheel_move_sample == 1) = NaN;

% Plot
go_left_trials = (signals_events.trialSideValues == 1 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == -1 & signals_events.hitValues == 0) & ...
    ~signals_events.repeatTrialValues;

go_right_trials = (signals_events.trialSideValues == -1 & signals_events.hitValues == 1) | ...
    (signals_events.trialSideValues == 1 & signals_events.hitValues == 0) & ...
    ~signals_events.repeatTrialValues;

trial_choice = go_left_trials + 2.*go_right_trials;

raster_window = [-2,2];
psth_bin_size = 0.01;

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

figure; hold on;
plot(nanmean(zscore(template_psth_left,[],2)),'k');
plot(nanmean(zscore(template_psth_right,[],2)),'r');

% Get left/right differences for correct non-repeat trials, plot raster and
% sort by difference
l_r_diff = (sum(template_psth_left,2) - sum(template_psth_right,2))./ ...
    (sum(template_psth_left,2) + sum(template_psth_right,2));

[~,sort_idx] = sort(l_r_diff);
sort_templates = nan(size(templates,1),1);
sort_templates(use_templates(sort_idx)) = 1:length(use_templates);
template_sort = sort_templates(spike_templates(use_spikes_idx));

raster_window = [-1,1];
psthViewer(use_spikes,template_sort, ...
    [wheel_move_time(go_right_trials),wheel_move_time(go_left_trials)], ...
    raster_window,[ones(1,sum(go_right_trials)),2*ones(1,sum(go_left_trials))]);


%% TO DO: average fluorescence around all conditions












