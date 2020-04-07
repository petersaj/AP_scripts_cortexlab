% Generate revision figures for ctx-str paper
% (trials data is prepared in AP_ctx_str_trial_preprocessing)

% NOTE: these are just in order that I wrote them at the moment

%% [[LOAD DATASETS]]

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


%% ~~~~~~~~ Widefield correlation borders

wf_corr_borders_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_borders\wf_corr_borders.mat';
load(wf_corr_borders_fn);

wf_corr_borders_cat = cell2mat(reshape([wf_corr_borders(:).corr_edges],1,1,[]));
figure;
imagesc(nanmean(wf_corr_borders_cat,3));
axis image off;
caxis([0,max(caxis)])
colormap(brewermap([],'Greys'))
ccf_outline = AP_reference_outline('ccf_aligned',[1,0,0]);
cellfun(@(x) set(x,'linewidth',2),vertcat(ccf_outline{:}));


%% ~~~~~~~~ Probe location variation

%% Plot widefield-estimated/histology-aligned probe location

animal_days = {...
    'AP032','2018-10-26';
    'AP033','2018-10-26';
    'AP034','2018-10-26';
    'AP036','2018-11-14';
    'AP043','2019-12-09';
    'AP060','2019-12-06';
    'AP061','2019-12-09'};

ccf_fig = figure;
plot_x = ceil(sqrt(length(animal_days)));
plot_y = plot_x;

wf_fig = figure;

for curr_animal_day = 1:length(animal_days)
    
    animal = animal_days{curr_animal_day,1};
    day = animal_days{curr_animal_day,2};
    
    probe_angle = 45; % from horizontal
    
    [img_path,img_exists] = AP_cortexlab_filename(animal,day,[],'imaging');
    avg_im = readNPY([img_path filesep 'meanImage_blue.npy']);
    
    avg_im_aligned = AP_align_widefield(avg_im,animal,day);
    
    h = figure('units','normalized','outerposition',[0 0 1 1]);imagesc(avg_im_aligned);
    axis image off;
    colormap(gray);caxis([0,prctile(avg_im_aligned(:),95)]);
    title('Click probe start/end');
    probe_wf = ginput(2);
    close(h)
    
    % Load and invert master CCF tform
    ccf_tform_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\ccf_tform'];
    load(ccf_tform_fn);
    ccf_tform_inverse = invert(ccf_tform);
    
    % Convert aligned widefield probe points to CCF
    % (master alignment to downsampled CCF, need to upsample)
    um2pixel = 20.6;
    probe_ccf_tformed = [probe_wf,[1;1]]*ccf_tform_inverse.T;
    probe_ccf = round(probe_ccf_tformed(:,1:2)*(um2pixel/10));
    
    % Plot probe points on CCF next to average image
    figure(wf_fig);
    
    subplot(length(animal_days),2,((curr_animal_day-1)*2)+1);
    imagesc(avg_im_aligned);
    axis image off;
    colormap(gray);caxis([0,prctile(avg_im_aligned(:),95)]);
    AP_reference_outline('ccf_aligned','r');
    line(probe_wf(:,1),probe_wf(:,2),'color','b','linewidth',1,'linestyle','--');
    title('Widefield');
    
    subplot(length(animal_days),2,((curr_animal_day-1)*2)+2);
    hold on; set(gca,'YDir','reverse');axis image;
    load('C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\cortical_area_boundaries.mat');
    for curr_area_idx =1:length(cortical_area_boundaries)
        p = cellfun(@(outline) plot(outline(:,2),outline(:,1),'color','k'), ...
            cortical_area_boundaries{curr_area_idx},'uni',false);
    end
    line(probe_ccf(:,1),probe_ccf(:,2),'color','r','linewidth',2);
    title('CCF');
    
    % Load in the annotated Allen volume and names
    allen_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
    av = readNPY([allen_path filesep 'annotation_volume_10um_by_index.npy']);
    st = loadStructureTree([allen_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
    
    % Estimate probe entry point by clicked probe start
    depth_start = find(av(probe_ccf(1,2),:,probe_ccf(1,1)) > 1,1);
    probe_entry_ccf = [probe_ccf(1,2),depth_start,probe_ccf(1,1)];
    
    % Estimate height of clicked point from user-set probe angle
    probe_sample_length = pdist2(probe_ccf(1,:),probe_ccf(2,:));
    probe_sample_height = round(probe_sample_length/tand(90-probe_angle));
    probe_air_ccf = [probe_ccf(2,2),depth_start-probe_sample_height,probe_ccf(2,1)];
    
    % Concatenate probe CCF coordinates (in up-down direction);
    probe_ccf = [probe_air_ccf;probe_entry_ccf];
    
    % Get estimated probe vector (widefield)
    r0 = mean(probe_ccf,1);
    xyz = bsxfun(@minus,probe_ccf,r0);
    [~,~,V] = svd(xyz,0);
    probe_direction = V(:,1);
    
    probe_vector_evaluate = [0,sign(probe_direction(2))*1000];
    probe_vector_ccf = round(bsxfun(@plus,bsxfun(@times,probe_vector_evaluate',probe_direction'),r0));
    
    % Get estimated probe vector (histology)
    [probe_filename,probe_filename_exists] = AP_cortexlab_filename(animal,[],[],'probe_ccf');
    load(probe_filename);
    
    use_probe = 1;
    histology_points = probe_ccf(1).points;
    r0 = mean(histology_points,1);
    xyz = bsxfun(@minus,histology_points,r0);
    [~,~,V] = svd(xyz,0);
    histology_probe_direction = V(:,1);
    
    probe_vector_evaluate = [-sign(probe_direction(2))*500,sign(probe_direction(2))*500];
    probe_vector = bsxfun(@plus,bsxfun(@times,probe_vector_evaluate',histology_probe_direction'),r0);
    probe_vector_ccf_histology = round(probe_vector);
    
    % histology done looking A->P so flipped, mirror about bregma
    bregma = allenCCFbregma;
    probe_vector_ccf_histology_mirrored = probe_vector_ccf_histology;
    probe_vector_ccf_histology_mirrored(:,3) = ...
        bregma(3) - (probe_vector_ccf_histology(:,3) - bregma(3));
    
    % Plot estimated probe location on CCF
    % (note the CCF is rotated to allow for dim 1 = x)
    figure(ccf_fig);
    ccf_axes = subplot(plot_y,plot_x,curr_animal_day); hold on
    slice_spacing = 10;
    target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) > 1,[2,1,3]);
    structure_patch = isosurface(target_volume,0);
    structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
    target_structure_color = [0.7,0.7,0.7];
    brain_outline = patch('Vertices',structure_wire.vertices*slice_spacing, ...
        'Faces',structure_wire.faces, ...
        'FaceColor','none','EdgeColor',target_structure_color);
    
    str_id = find(strcmp(st.safe_name,'Caudoputamen'));
    target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) == str_id,[2,1,3]);
    structure_patch = isosurface(target_volume,0);
    structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
    target_structure_color = [0.7,0,0.7];
    striatum_outline = patch('Vertices',structure_wire.vertices*slice_spacing, ...
        'Faces',structure_wire.faces, ...
        'FaceColor','none','EdgeColor',target_structure_color);
    
    axis image vis3d off;
    view([-30,25]);
    cameratoolbar(ccf_fig,'SetCoordSys','y');
    cameratoolbar(ccf_fig,'SetMode','orbit');
    title([animal ' ' day]);
    
    % Plot current probe vector
    probe_wf_line = line(ccf_axes,probe_vector_ccf(:,1),probe_vector_ccf(:,2),probe_vector_ccf(:,3),'color','k','linewidth',3);
    % (make sure histology probe enters on the left)
    if  histology_probe_direction(3) < 0
        probe_histology_line = line(ccf_axes, ...
            probe_vector_ccf_histology(:,1), ...
            probe_vector_ccf_histology(:,2), ...
            probe_vector_ccf_histology(:,3),'color','r','linewidth',3);
    else
        probe_histology_line = line(ccf_axes, ...
            probe_vector_ccf_histology_mirrored(:,1), ...
            probe_vector_ccf_histology_mirrored(:,2), ...
            probe_vector_ccf_histology_mirrored(:,3),'color','r','linewidth',3);
    end
    
    drawnow;
    
    if curr_animal_day == 1
        legend([probe_wf_line,probe_histology_line],{'Widefield-estimated','Histology'},'location','nw');
    end
    
end


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
















