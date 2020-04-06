% Generate revision figures for ctx-str paper 
% (trials data is prepared in AP_ctx_str_trial_preprocessing)

% NOTE: these are just in order that I wrote them at the moment

%% Load in task data

% Load data

% (task)
% data_fn = 'trial_activity_choiceworld'; % Primary dataset
% data_fn = 'trial_activity_choiceworld_4strdepth'; % Depth-aligned striatum
% exclude_data = true;

% (passive)
data_fn = 'trial_activity_AP_choiceWorldStimPassive_trained';
% data_fn = 'trial_activity_AP_choiceWorldStimPassive_naive';
% data_fn = 'trial_activity_stimKalatsky_naive';
% data_fn = 'trial_activity_stimKalatsky_trained';
exclude_data = false;

% (unused at the moment)
% data_fn = 'trial_activity_choiceworld_wfonly'; % Widefield-only days (no craniotomy, so cleaner)
% exclude_data = true;

AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_recording;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));


%% ~~~~~~~~ Widefield + Striatum + Cortex ephys


%% Correlation between wf/ctx-mua and ctx-mua/str-mua

use_protocol = 'vanillaChoiceworld';
% use_protocol = 'AP_sparseNoise';

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
data_fn = [data_path filesep 'ctx_fluor_mua_corr_' use_protocol];
load(data_fn); 

mua_depth = data(1).cortex_mua_depth{1}; % they're all the same, use 1st
cortex_fluor_corr_cat = cell2mat(horzcat(data.cortex_fluor_corr));
cortex_striatum_corr_cat = cell2mat(permute(horzcat(data.cortex_striatum_corr),[1,3,2]));

figure; 

subplot(1,3,1,'YDir','reverse'); hold on;
plot(cortex_fluor_corr_cat',mua_depth,'color',[0.2,0.8,0.2]);
errorbar(nanmean(cortex_fluor_corr_cat,2), ...
    mua_depth,AP_sem(cortex_fluor_corr_cat,2), ...
    'horizontal','color',[0,0.6,0],'linewidth',2)
xlabel('Correlation');
ylabel('Cortical MUA aligned depth');
title('Cortical fluorescence');

subplot(1,3,2,'YDir','reverse'); hold on;
set(gca,'ColorOrder',copper(4));
plot_str = 1;
plot(permute(cortex_striatum_corr_cat(:,plot_str,:),[1,3,2]),mua_depth,'color',[0.5,0.5,0.5]);
errorbar(nanmean(cortex_striatum_corr_cat(:,plot_str,:),3), ...
    mua_depth,AP_sem(cortex_striatum_corr_cat(:,plot_str,:),3),'k','horizontal','linewidth',2)
xlabel('Correlation');
ylabel('Cortical MUA aligned depth');
title(['Str ' num2str(plot_str) ' multiunit']);

subplot(2,3,3); hold on;
for i = 1:size(cortex_fluor_corr_cat,2)
    plot(cortex_fluor_corr_cat(:,i), ...
        cortex_striatum_corr_cat(:,plot_str,i), ...
        'color',[0.5,0.5,0.5]);
end
plot(nanmean(cortex_fluor_corr_cat,2), ...
    nanmean(cortex_striatum_corr_cat(:,plot_str,:),3),'k','linewidth',2);
xlabel('Fluorescence - cortical MUA correlation');
ylabel('Cortical MUA - striatal MUA correlation')

subplot(2,3,6); hold on;
plot(permute(max(cortex_striatum_corr_cat,[],1),[2,3,1]),'color',[0.5,0.5,0.5]);
errorbar(squeeze(nanmean(max(cortex_striatum_corr_cat,[],1),3)), ...
    squeeze(AP_sem(max(cortex_striatum_corr_cat,[],1),3)),'k','linewidth',2);
xlim([0.5,4.5]);
xlabel('Striatal domain');
ylabel('Cortical MUA max corr');


%% IN PROGRESS: example recording

% AP060 2019-12-06 looks like the best?

animal = 'AP060'; 
day = '2019-12-06'; 
experiment = 3; 
site = 2; % (cortex)
str_align= 'none'; % (cortex)
verbose = true; 
AP_load_experiment;




% Plot CSD
vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
load(vis_ctx_ephys_fn);

curr_animal_idx = strcmp(animal,{vis_ctx_ephys.animal});
curr_day_idx = strcmp(day,vis_ctx_ephys(curr_animal_idx).day);

figure;
imagesc(vis_ctx_ephys(curr_animal_idx).stim_lfp_t{curr_day_idx}, ...
    vis_ctx_ephys(curr_animal_idx).stim_csd_depth{curr_day_idx}, ...
    vis_ctx_ephys(curr_animal_idx).stim_csd{curr_day_idx});
caxis([-max(abs(caxis))/2,max(abs(caxis))/2]);
colormap(brewermap([],'*RdBu'));
ylabel('Depth (\mum)');
xlabel('Time from stim');
colorbar;

% (COPIED FROM ABOVE: PLOT CORTEX MULTIUNIT AND FLUORESCENCE)

%%% DEPTH-ALIGN TEMPLATES, FIND CORTEX BOUNDARY
curr_animal_idx = strcmp(animal,{vis_ctx_ephys.animal});
curr_day_idx = strcmp(day,vis_ctx_ephys(curr_animal_idx).day);
curr_csd_depth = vis_ctx_ephys(curr_animal_idx).stim_csd_depth{curr_day_idx};
curr_csd_depth_aligned = vis_ctx_ephys(curr_animal_idx).stim_csd_depth_aligned{curr_day_idx};
template_depths_aligned = interp1(curr_csd_depth,curr_csd_depth_aligned,template_depths);

% Find cortex end by largest gap between templates
sorted_template_depths = sort([template_depths_aligned]);
[max_gap,max_gap_idx] = max(diff(sorted_template_depths));
ctx_end = sorted_template_depths(max_gap_idx)+1;

ctx_depth = [sorted_template_depths(1),ctx_end];
ctx_units = template_depths_aligned <= ctx_depth(2);

%%% GET FLUORESCENCE AND SPIKES BY DEPTH

% Set binning time
skip_seconds = 60;
spike_binning_t = 1/framerate; % seconds
spike_binning_t_edges = frame_t(1)+skip_seconds:spike_binning_t:frame_t(end)-skip_seconds;
spike_binning_t_centers = spike_binning_t_edges(1:end-1) + diff(spike_binning_t_edges)/2;

% Get fluorescence in pre-drawn ROI
curr_ctx_roi = vis_ctx_ephys(curr_animal_idx).ctx_roi{curr_day_idx};

fVdf_deconv = AP_deconv_wf(fVdf);
fluor_roi = AP_svd_roi(Udf,fVdf_deconv,avg_im,[],curr_ctx_roi);
fluor_roi_interp = interp1(frame_t,fluor_roi,spike_binning_t_centers);

% Set sliding depth window of MUA
depth_corr_range = [-200,1500];
depth_corr_window = 200; % MUA window in microns
depth_corr_window_spacing = 50; % MUA window spacing in microns

depth_corr_bins = [depth_corr_range(1):depth_corr_window_spacing:(depth_corr_range(2)-depth_corr_window); ...
    (depth_corr_range(1):depth_corr_window_spacing:(depth_corr_range(2)-depth_corr_window))+depth_corr_window];
depth_corr_bin_centers = depth_corr_bins(1,:) + diff(depth_corr_bins,[],1)/2;

cortex_mua = zeros(size(depth_corr_bins,2),length(spike_binning_t_centers));
for curr_depth = 1:size(depth_corr_bins,2)
    curr_depth_templates_idx = ...
        find(ctx_units & ...
        template_depths_aligned >= depth_corr_bins(1,curr_depth) & ...
        template_depths_aligned < depth_corr_bins(2,curr_depth));
    
    cortex_mua(curr_depth,:) = histcounts(spike_times_timeline( ...
        ismember(spike_templates,curr_depth_templates_idx)),spike_binning_t_edges);
end

%%% LOAD STRIATUM EPHYS AND GET MUA BY DEPTH
clear load_parts
load_parts.ephys = true;
site = 1; % (striatum is always on probe 1)
str_align = 'kernel';
AP_load_experiment;

striatum_mua = nan(n_aligned_depths,length(spike_binning_t_centers));
for curr_depth = 1:n_aligned_depths
    curr_spike_times = spike_times_timeline(aligned_str_depth_group == curr_depth);
    % Skip if no spikes at this depth
    if isempty(curr_spike_times)
        continue
    end
    striatum_mua(curr_depth,:) = histcounts(curr_spike_times,spike_binning_t_edges);
end

figure;
subplot(5,1,1);
plot(spike_binning_t_centers,fluor_roi_interp,'linewidth',2,'color',[0,0.7,0]);
title('Fluorescence');
subplot(5,1,2:4);
imagesc(spike_binning_t_centers,[],cortex_mua)
caxis([0,10]);
colormap(brewermap([],'Greys'));
title('Cortex MUA');
subplot(5,1,5);
imagesc(spike_binning_t_centers,[],striatum_mua);
caxis([0,10]);
title('Striatum MUA');

linkaxes(get(gcf,'Children'),'x');

plot_t = [128,132];
xlim(plot_t);
















