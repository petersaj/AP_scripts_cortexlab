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
    surround_time = surround_window(1):surround_samplerate:surround_window(2);
    baseline_surround_time = baseline_surround_window(1):surround_samplerate:baseline_surround_window(2);
    
    % Average (time course) responses
    conditions = unique(stimIDs);
    im_stim = nan(size(U,1),size(U,2),length(surround_time),length(conditions));
    for curr_condition_idx = 1:length(conditions)
        curr_condition = conditions(curr_condition_idx);
        
        use_stims = find(stimIDs == curr_condition);
        use_stim_onsets = stim_onsets(use_stims(2:end));
        use_stim_onsets([1,end]) = [];
        
        stim_surround_times = bsxfun(@plus, use_stim_onsets(:), surround_time);
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
surround_time = surround_window(1):surround_samplerate:surround_window(2);
baseline_surround_time = baseline_surround_window(1):surround_samplerate:baseline_surround_window(2);

AP_image_scroll(a,surround_time)


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
    surround_time = [-0.5,5];
    surround_samples = surround_time/Timeline.hw.samplingInterval;
    
    % Get wheel aligned to stim onset
    rotaryEncoder_idx = strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder');
    surround_time = surround_time(1):Timeline.hw.samplingInterval:surround_time(2);
    pull_times = bsxfun(@plus,stim_onsets,surround_time);
    
    stim_aligned_wheel_raw = interp1(Timeline.rawDAQTimestamps, ...
        Timeline.rawDAQData(:,rotaryEncoder_idx),pull_times);
    stim_aligned_wheel = bsxfun(@minus,stim_aligned_wheel_raw, ...
        nanmedian(stim_aligned_wheel_raw(:,surround_time < 0),2));

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
    surround_time = surround_window(1):surround_samplerate:surround_window(2);
    baseline_surround_time = baseline_surround_window(1):surround_samplerate:baseline_surround_window(2);
    
    % Average (time course) responses
    %     conditions = unique(stimIDs);
    conditions = [-2,-1,0,1,2];
    im_stim_hit = nan(size(U,1),size(U,2),length(surround_time),length(conditions));
    im_stim_miss = nan(size(U,1),size(U,2),length(surround_time),length(conditions));

    for curr_condition_idx = 1:length(conditions)
        curr_condition = conditions(curr_condition_idx);
        
        use_stims = find(stimIDs == curr_condition & signals_events.hitValues(1:num_stim) == 0);
        use_stim_onsets = stim_onsets(use_stims);       
        if length(use_stim_onsets) > 5           
            stim_surround_times = bsxfun(@plus, use_stim_onsets(:), surround_time);
            peri_stim_v = permute(mean(interp1(frame_t,fVdf',stim_surround_times),1),[3,2,1]);           
            im_stim_miss(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
        end
        
        use_stims = find(stimIDs == curr_condition & signals_events.hitValues(1:num_stim) == 1);
        use_stim_onsets = stim_onsets(use_stims);       
        if length(use_stim_onsets) > 5           
            stim_surround_times = bsxfun(@plus, use_stim_onsets(:), surround_time);
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
surround_time = surround_window(1):surround_samplerate:surround_window(2);
baseline_surround_time = baseline_surround_window(1):surround_samplerate:baseline_surround_window(2);

AP_image_scroll(avg_hit,surround_time)
AP_image_scroll(avg_miss,surround_time)


%% PSTH to choiceworld conditions

stimIDs = signals_events.trialSideValues.*signals_events.trialContrastValues;
stimIDs = discretize(stimIDs,[-Inf,-0.125,-0.01,0.01,0.25,Inf],[-2,-1,0,1,2]);

use_spikes_idx = ismember(spike_templates,find(templateDepths >= 700 & templateDepths <= 1500));
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
use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 1;
align_times = signals_events.stimOnTimes(use_trials(1:length(signals_events.stimOnTimes)))';

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

aligned_V = reshape(fVdf(:,align_frames'), ...
    size(fV,1),size(align_frames,2),size(align_frames,1));

mean_aligned_V = nanmean(aligned_V,3);

% Get and plot the average fluorescence around event
mean_aligned_px = svdFrameReconstruct(Udf,mean_aligned_V);



a = diff(imgaussfilt(mean_aligned_px,2),[],3);
a(a < 0) = 0;
AP_image_scroll(a,surround_time);
warning off; truesize; warning on;











