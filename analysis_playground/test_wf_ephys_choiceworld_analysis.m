%% Batch load responses to passive stim

animal = 'AP025';
% protocol = 'vanillaChoiceworld';
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
        use_stim_onsets = stim_onsets(use_stims(2:end));
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

%% Align images from batch processing

batch_vars_reg = batch_vars;

for curr_day = 1:length(experiments)
    curr_im = batch_vars_reg.im_stim{curr_day};
    curr_im(isnan(curr_im)) = 0;
    curr_im = imgaussfilt(diff(curr_im,[],3),2);
    curr_im(curr_im < 0) = 0;
    batch_vars_reg.im_stim{curr_day} = curr_im;
end

for curr_day = setdiff(1:length(experiments),ref_im_num);
    
    tform = affine2d;
    tform.T = tform_matrix{curr_day};
    
    curr_im = batch_vars_reg.im_stim{curr_day};
    curr_im(isnan(curr_im)) = 0;    

    batch_vars_reg.im_stim{curr_day} = imwarp(curr_im,tform, ...
        'Outputview',imref2d(size(avg_im{ref_im_num})));       

end

a = nanmean(cat(5,batch_vars_reg.im_stim{:}),5);

surround_window = [-0.2,3];
baseline_surround_window = [0,0];
framerate = 35;
surround_samplerate = 1/(framerate*1);
t_surround = surround_window(1):surround_samplerate:surround_window(2);
baseline_surround_time = baseline_surround_window(1):surround_samplerate:baseline_surround_window(2);

AP_image_scroll(a,t_surround)


%% Batch get average choiceworld fluorescence

animal = 'AP025';
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
    %     conditions = unique(stimIDs);
    conditions = [-2,-1,0,1,2];
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

avg_im = cell(length(days),1);
for curr_day = 1:length(days)
    data_path = ['\\zserver.cortexlab.net\Data\Subjects\' animal filesep days{curr_day}];
    avg_im{curr_day} = readNPY([data_path filesep 'meanImage_blue.npy']);
end

border_pixels = 20;

%im_align = cellfun(@(x) x(border_pixels:end-border_pixels+1,border_pixels:end-border_pixels+1),avg_im,'uni',false);
% align the left half of the image (without the craniotomy)
%im_align = cellfun(@(x) x(border_pixels:end,1:round(size(x,2)/2)),avg_im,'uni',false);
im_align = cellfun(@(x) imgaussfilt(x(border_pixels:end-border_pixels+1,border_pixels:end-border_pixels+1),3),avg_im,'uni',false);

% Choose reference day
ref_im_num = round(length(im_align)/2);
%ref_im_num = length(avg_im);

disp('Registering average images')
tform_matrix = cell(length(avg_im),1);
tform_matrix{1} = eye(3);

avg_im_reg = nan(size(avg_im{ref_im_num},1),size(avg_im{ref_im_num},2),length(avg_im));
avg_im_reg(:,:,ref_im_num) = avg_im{ref_im_num};

for curr_session = setdiff(1:length(avg_im),ref_im_num)
    
    % This is to do correlation, then affine (if above doesn't work)
    [optimizer, metric] = imregconfig('monomodal');
    optimizer = registration.optimizer.OnePlusOneEvolutionary();
    optimizer.MaximumIterations = 200;
    optimizer.GrowthFactor = 1+1e-6;
    optimizer.InitialRadius = 1e-4;

    %%% for just affine
    tformEstimate_affine = imregtform(im_align{curr_session},im_align{ref_im_num},'affine',optimizer,metric);
    curr_im_reg = imwarp(avg_im{curr_session},tformEstimate_affine,'Outputview',imref2d(size(avg_im{ref_im_num})));
    tform_matrix{curr_session} = tformEstimate_affine.T;
    %%%%
    
    avg_im_reg(:,:,curr_session) = curr_im_reg;
    
    disp(curr_session);
    
end

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
for curr_day = setdiff(1:length(experiments),ref_im_num);
    
    tform = affine2d;
    tform.T = tform_matrix{curr_day};
    
    curr_im = batch_vars_reg.im_stim_miss{curr_day};
    nan_cond = squeeze(all(all(all(isnan(curr_im),1),2),3));
    curr_im(isnan(curr_im)) = 0;
    curr_im = imwarp(curr_im,tform, ...
        'Outputview',imref2d(size(avg_im{ref_im_num})));
    curr_im(:,:,:,nan_cond) = NaN;
    batch_vars_reg.im_stim_miss{curr_day} = curr_im;
    
    curr_im = batch_vars_reg.im_stim_hit{curr_day};
    nan_cond = squeeze(all(all(all(isnan(curr_im),1),2),3));
    curr_im(isnan(curr_im)) = 0;
    curr_im = imwarp(curr_im,tform, ...
        'Outputview',imref2d(size(avg_im{ref_im_num}))); 
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

use_spikes_idx = ismember(spike_templates,find(templateDepths >= 500 & templateDepths <= 1500));
% use_spikes_idx = ismember(spike_templates,intersect(find(templateDepths >= 500 & templateDepths <= 1500),find(msn)));

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
line([0,0],ylim,'linestyle','--','color','k');

% Plot the hit PSTHs on top of each other
figure; hold on
set(gca,'ColorOrder',copper(length(unique(stimIDs))));
plot(bins(20:end-20),stim_psth_hit_smooth(:,20:end-20)','linewidth',2)

%% Plot contrast-dependent response by depth

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
    
%     curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    curr_spike_times = spike_times_timeline((depth_group == curr_depth) & ...
        ismember(spike_templates,find(msn)));
        
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
n_depth_groups = 8;
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
    
%     curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    curr_spike_times = spike_times_timeline((depth_group == curr_depth) & ...
        ismember(spike_templates,find(msn)));
    
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
line([median_move_time,median_move_time],ylim,'color','k','linestyle','--');
line([median_reward_time,median_reward_time],ylim,'color','k','linestyle','--');

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

%use_spikes_idx = ismember(spike_templates,find(templateDepths >= 3000 & templateDepths <= 4000));
use_spikes_idx = ismember(spike_templates,intersect(find(templateDepths >= 1000 & templateDepths <= 2000),find(msn)));

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

% use_spikes_idx = ismember(spike_templates,find(templateDepths >= 0 & templateDepths <= 4000));
use_spikes_idx = ismember(spike_templates,intersect(find(templateDepths >= 2000 & templateDepths <= 3000),find(msn)));
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

% Define times to align
use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues == 0 & signals_events.hitValues == 0;
align_times = reshape(stimOn_times(use_trials(1:length(stimOn_times))),[],1);

t_surround = [-0.5,1.5];
samplerate = framerate*2;
t_surround = t_surround(1):1/samplerate:t_surround(2);
t_peri_event = bsxfun(@plus,align_times,t_surround);

% Draw ROI and align fluorescence
roi_trace = AP_svd_roi(Udf,fVdf,avg_im);
event_aligned_f = interp1(frame_t,roi_trace,t_peri_event);
event_aligned_df = interp1(conv(frame_t,[1,1]/2,'valid'),diff(roi_trace),t_peri_event);
event_aligned_df(event_aligned_df < 0) = 0;

% Pull out MUA at a given depth
use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 0 & templateDepths < 1500)));
% use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 0 & templateDepths < 1500)) &...
%     ismember(spike_templates,find(msn)));
t_peri_event_bins = [t_peri_event - 1/(samplerate*2), ...
    t_peri_event(:,end) + 1/(samplerate*2)];
event_aligned_spikes = reshape([histcounts(use_spikes,reshape(t_peri_event_bins',[],1)),NaN],[],length(align_times))';
event_aligned_spikes(:,length(t_surround)+1) = [];

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


%% SPECIFIC VERSION OF ABOVE (hit vs. miss) 

% Define times to align

% SIGNALS - CHOICEWORLD
% use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues > 0;% & ~isnan(wheel_move_time);
% align_times = reshape(stimOn_times(use_trials(1:length(stimOn_times))),[],1);
% hit_trials = signals_events.hitValues(use_trials) == 1;
% stim_conditions = signals_events.trialContrastValues(use_trials);

% MPEP PASSIVE
align_times = reshape(stimOn_times,[],1);
hit_trials = rand(size(align_times)) > 0.5;
stim_conditions = stimIDs;

interval_surround = [-0.5,1.5];
samplerate = framerate*2;
t_surround = interval_surround(1):1/samplerate:interval_surround(2);
t_peri_event = bsxfun(@plus,align_times,t_surround);

% Draw ROI and align fluorescence
[roi_trace,roi_mask] = AP_svd_roi(Udf,fVdf,weight_im,response_im);
event_aligned_f = interp1(frame_t,roi_trace,t_peri_event);
event_aligned_df = interp1(conv(frame_t,[1,1]/2,'valid'),diff(roi_trace),t_peri_event);
event_aligned_df(event_aligned_df < 0) = 0;

% Pull out MUA at a given depth
use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 500 & templateDepths < 1200)));
% use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 500 & templateDepths < 1200)) &...
%     ismember(spike_templates,find(msn)));
t_peri_event_bins = [t_peri_event - 1/(samplerate*2), ...
    t_peri_event(:,end) + 1/(samplerate*2)];
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
frame_t_df_resample = frame_t_df(1):1/samplerate:frame_t_df(end);
roi_trace_df_resample = interp1(frame_t_df,roi_trace_df,frame_t_df_resample);
roi_trace_df_resample(roi_trace_df < 0) = 0;
spike_bins = [frame_t_df_resample-1/samplerate,frame_t_df_resample(end)+1/samplerate];
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




















