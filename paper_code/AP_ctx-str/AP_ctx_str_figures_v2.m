% Generate figures for ctx-str paper

% Anything that takes a lot of time is done in
% AP_ctx_str_trial_preprocessing and saved for plotting here

% The original scripts here were in test_wf_ephys_choiceworld_analysis


%% Fig 1b left: Example cortex>striatum regression by depth

% Set example experiment to use
animal = 'AP025';
day = '2017-10-04';
experiment = 1;

% Parameters for regression
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 60;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

% Load full data
str_align = 'none';
verbose = true;
AP_load_experiment

% Get time points to query
sample_rate = framerate*regression_params.upsample_factor;
time_bins = frame_t(find(frame_t > ...
    regression_params.skip_seconds,1)):1/sample_rate: ...
    frame_t(find(frame_t-frame_t(end) < ...
    -regression_params.skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Deconvolve and resample V
fVdf_deconv = AP_deconv_wf(fVdf);
fVdf_deconv(isnan(fVdf_deconv)) = 0;
fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';

% Get striatum multiunit in 4 depths for example plot
n_depths = 4;
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
[depth_group_n,depth_group] = histc(spike_depths,depth_group_edges);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

binned_spikes = zeros(n_depths,length(time_bins)-1);
for curr_depth = 1:n_depths
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
end
binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);

% Load lambda from previously estimated and saved
lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
load(lambda_fn);
curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
if any(curr_animal_idx)
    curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
    if any(curr_day_idx)
        lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
    end
end

% Regress fluorescence to spikes
kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
    round(regression_params.kernel_t(2)*sample_rate);

[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
    binned_spikes_std,kernel_frames,lambda, ...
    regression_params.zs,regression_params.cvfold, ...
    false,regression_params.use_constant);

Udf_aligned = single(AP_align_widefield(animal,day,Udf));
k_px = zeros(size(Udf_aligned,1),size(Udf_aligned,2),size(k,2),size(k,3),'single');
for curr_spikes = 1:size(k,3)
    k_px(:,:,:,curr_spikes) = ...
        svdFrameReconstruct(Udf_aligned(:,:,regression_params.use_svs),k(:,:,curr_spikes));
end

% NaN-out depths with no spikes
k_px(:,:,:,~any(binned_spikes,2)) = NaN;

% Keep kernel one frame (t == 0)
k_px_frame = squeeze(k_px(:,:,kernel_frames == 0,:));

% Define ROIs and get fluorescence traces
roi_circle_size = 20;
roi_x = [136,176,142,116]; % roi_x = [131,174,110,51];
roi_y = [299,91,79,87]; % roi_y = [297,96,71,144];
[x,y] = meshgrid(1:size(Udf_aligned,1),1:size(Udf_aligned,2));
roi_mask = cell2mat(arrayfun(@(roi) sqrt((x-roi_x(roi)).^2 + (y-roi_y(roi)).^2) <= ...
    roi_circle_size,permute(1:length(roi_x),[1,3,2]),'uni',false));
roi_trace = AP_svd_roi(Udf_aligned,fVdf_deconv_resample,[],[],roi_mask);

% Plot ROIs
figure;
for i = 1:n_depths
   subplot(n_depths,1,i,'YDir','reverse'); hold on;
   curr_roi_boundaries = cell2mat(bwboundaries(roi_mask(:,:,i)));
   plot(curr_roi_boundaries(:,2),curr_roi_boundaries(:,1),'color',[0,0.8,0],'linewidth',2);
   axis image off;
   AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
end

% Plot overlaid striatal activity and and ROI fluorescence
figure; hold on
AP_stackplot(binned_spikes',time_bin_centers,10,true,'k');
AP_stackplot(roi_trace',time_bin_centers,10,true,[0,0.7,0]);
xlim([185,205]);

% Plot kernel frames
figure;
for i = 1:n_depths
   subplot(n_depths,1,i);
   imagesc(k_px_frame(:,:,i)); hold on;
   axis image off;
   colormap(brewermap([],'*RdBu'));
   caxis([-0.012,0.012]); 
   AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
end

% (UNUSED: was going get spike-triggered average for regression comparison)
% surround_times = [-0.2,0.2];
% 
% surround_frames = round(surround_times(1)*framerate):round(surround_times(2)*framerate);
% sta_t = surround_frames./framerate;
% sta_im = zeros(size(Udf_aligned,1),size(Udf_aligned,2), ...
%     length(surround_frames),size(binned_spikes,1));
% 
% for curr_depth = 1:size(binned_spikes,1)
%     frames_w = repmat(binned_spikes(curr_depth,:)'./ ...
%         sum(binned_spikes(curr_depth,:)),1,length(surround_frames));
%     for curr_sta_frame = 1:length(surround_frames)
%         frames_w(:,curr_sta_frame) = ...
%             circshift(frames_w(:,curr_sta_frame),[surround_frames(curr_sta_frame),0]);
%     end
%     
%     sta_v = fVdf_deconv_resample*frames_w;
%     sta_im(:,:,:,curr_depth) = svdFrameReconstruct(Udf_aligned,sta_v);
% end
% 
% sta_im_max = squeeze(max(sta_im,[],3));


% (UNUSED: was going to plot examples of kernels by 200um depths)
% % Load kernels by depths
% kernel_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
% kernel_fn = ['ephys_kernel_depth'];
% load([kernel_path filesep kernel_fn])
% 
% % Plot templates across striatum depth for example experiment
% curr_animal_idx = strcmp(animal,{ephys_kernel_depth.animal});
% curr_day_idx = strcmp(day,ephys_kernel_depth(curr_animal_idx).days);
% ctx_str_depth_kernels = ephys_kernel_depth(curr_animal_idx).k_px{curr_day_idx};
% 
% figure;
% for i = 1:size(ctx_str_depth_kernels,3)
%     subplot(size(ctx_str_depth_kernels,3),1,i);
%     imagesc(ctx_str_depth_kernels(:,:,i));
%     axis image off;
%     colormap(brewermap([],'*RdBu'));
%     caxis([-0.01,0.01]);
% end


%% Fig 1c-d: Average cortex>striatum domain regression kernels

protocols = {'vanillaChoiceworld','stimSparseNoiseUncorrAsync'};

for protocol = protocols
    
    curr_protocol = cell2mat(protocol);
    
    n_aligned_depths = 4;
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_ephys';
    k_fn = [data_path filesep 'wf_ephys_maps_' curr_protocol '_' num2str(n_aligned_depths) '_depths'];
    load(k_fn);
    
    framerate = 35;
    upsample_factor = 1;
    sample_rate = framerate*upsample_factor;
    kernel_t = [-0.5,0.5];
    kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
    t = kernel_frames/sample_rate;
    
    % Concatenate explained variance
    expl_var_experiment = cell2mat(horzcat(batch_vars.explained_var));
    expl_var_animal = cell2mat(cellfun(@(x) nanmean(cell2mat(x),2),{batch_vars.explained_var},'uni',false));
    figure('Name',curr_protocol);
    errorbar(nanmean(expl_var_experiment,2), ...
        nanstd(expl_var_experiment,[],2)./sqrt(nansum(expl_var_experiment,2)),'k','linewidth',2);
    xlabel('Striatal depth');
    ylabel('Fraction explained variance');
    
    % Concatenate and mean
    % (kernel goes backwards in time - flip to correct)
    k_px_cat = cellfun(@(x) x(:,:,end:-1:1,:),[batch_vars.r_px],'uni',false);
    k_px = nanmean(double(cat(5,k_px_cat{:})),5);
    
    % Get center-of-mass maps
    k_px_positive = k_px;
    k_px_positive(k_px_positive < 0) = 0;
    k_px_com = sum(k_px_positive.*permute(1:n_aligned_depths,[1,3,4,2]),4)./sum(k_px_positive,4);
    k_px_com_colored = nan(size(k_px_com,1),size(k_px_com,2),3,size(k_px_com,3));
    
    jet_basic = jet(255);
    dark_colors = max(jet_basic,[],2) ~= 1;
    jet_alt = interp1(1:255,jet_basic,linspace(find(~dark_colors,1,'first'), ...
        find(~dark_colors,1,'last'),255)) - 0.2;
    use_colormap = jet_alt;
    
    for curr_frame = 1:size(k_px_com,3)
        k_px_com_colored(:,:,:,curr_frame) = ...
            ind2rgb(round(mat2gray(k_px_com(:,:,curr_frame),[1,n_aligned_depths])*255),use_colormap);
    end
    
    % Plot center kernel frames independently at t = 0
    figure('Name',curr_protocol);
    plot_frame = kernel_frames == 0;
    for curr_depth = 1:n_aligned_depths
       subplot(1,n_aligned_depths,curr_depth);
       imagesc(k_px(:,:,plot_frame,curr_depth));
       AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
       axis image off;
       colormap(brewermap([],'*RdBu'));
       caxis([-0.01,0.01]);
    end
    
    % Plot center-of-mass color across time
    plot_t = find(t >= -0.1 & t <= 0.1);
    weight_max = 0.005;
    figure('Name',curr_protocol);
    for t_idx = 1:length(plot_t)
        curr_t = plot_t(t_idx);
        subplot(1,length(plot_t),t_idx);
        p = image(k_px_com_colored(:,:,:,curr_t));
        set(p,'AlphaData', ...
            mat2gray(max(k_px(:,:,curr_t,:),[],4),[0,weight_max]));
        axis image off;
        AP_reference_outline('ccf_aligned','k');
        title(t(curr_t));
    end
    
    drawnow;
    
end


%% Fig 1e: Cortex>striatum projections (Allen connectivity database)

% Define the probe vector manually according to the targeted trajectory
probe_vector_ccf = [520,240,510;520,511,239];

%%% Get the average relative depth of each kernel template

% Load kernels by depths, get depth relative to maximum extent
kernel_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
kernel_fn = ['ephys_kernel_depth'];
load([kernel_path filesep kernel_fn])
total_depths = 1:max(cellfun(@(x) size(x,3),[ephys_kernel_depth.k_px]));
k_px_depths = cellfun(@(x) total_depths(end-size(x,3)+1:end),[ephys_kernel_depth.k_px],'uni',false);

% Load the kernel template matches
n_aligned_depths = 4;
kernel_match_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
kernel_match_fn = ['ephys_kernel_align_' num2str(n_aligned_depths) '_depths.mat'];
load([kernel_match_path filesep kernel_match_fn]);

% Concatenate all relative depths and kernel matches
k_depths = cell2mat(k_px_depths);
k_matches = cell2mat([ephys_kernel_align.kernel_match]')';
k_match_depths_relative = grpstats(k_depths,k_matches)./max(total_depths);

%%% Query allen at each point using targeted trajectory

% Load in the annotated Allen volume and names
allen_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
av = readNPY([allen_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean

% Get probe location per micron
probe_size = pdist2(probe_vector_ccf(1,:),probe_vector_ccf(2,:))*10;
probe_depths = ...
    round([linspace(probe_vector_ccf(1,1)',probe_vector_ccf(2,1)',probe_size); ...
    linspace(probe_vector_ccf(1,2)',probe_vector_ccf(2,2)',probe_size); ...
    linspace(probe_vector_ccf(1,3)',probe_vector_ccf(2,3)',probe_size)]');

% Eliminiate trajectory points that are off the atlas
eliminate_depths = ...
    probe_depths(:,1) < 1 | probe_depths(:,1) > size(av,1) | ...
    probe_depths(:,2) < 1 | probe_depths(:,2) > size(av,2) | ...
    probe_depths(:,3) < 1 | probe_depths(:,3) > size(av,3);
probe_depths(eliminate_depths,:) = [];

% Convert probe depths subscripts to indicies
probe_depths_ind = sub2ind(size(av),probe_depths(:,1),probe_depths(:,2),probe_depths(:,3));

% Get structures that the probe is in
probe_structures = av(probe_depths_ind);

% Get target relative depths through striatum
str_id = find(strcmp(st.safe_name,'Caudoputamen'));
probe_structures_str = probe_structures == str_id;
probe_str = [probe_depths(find(probe_structures_str,1,'first'),:); ...
    probe_depths(find(probe_structures_str,1,'last'),:)];
kernel_depth_ccf = interp1([0,1],probe_str,k_match_depths_relative);

%%%% Just use regular depths?
regular_centers_borders = linspace(0,1,n_aligned_depths*2+1);
kernel_depth_ccf = interp1([0,1],probe_str,regular_centers_borders(2:2:end));

% Plot brain to overlay probes
% (note the CCF is rotated to allow for dim 1 = x)
h = figure; ccf_axes = axes; hold on
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
cameratoolbar(h,'SetCoordSys','y');
cameratoolbar(h,'SetMode','orbit');

scatter3(kernel_depth_ccf(:,1),kernel_depth_ccf(:,2),kernel_depth_ccf(:,3), ...
    100,copper(n_aligned_depths),'filled');

% Mirror all of the locations across the midline and get projections from
% both (because the cortical injections in the database aren't evenly
% distributed, so this is more comprehensive)
kernel_depth_um = round(kernel_depth_ccf*10);
bregma_um = allenCCFbregma*10;
ccf_midline = bregma_um(3);
hemisphere = sign(kernel_depth_um(1,3) - ccf_midline);
str_depths_mirror = [kernel_depth_um(:,1:2),ccf_midline - hemisphere*abs(kernel_depth_um(:,3)-ccf_midline)];

str_depths_query = [kernel_depth_um;str_depths_mirror];
injection_parameters = get_allen_projection(str_depths_query);
injection_coordinates = {injection_parameters.coordinates};

% Standardize injection coordinates by hemisphere (left = contra, right =
% ipsi)
injection_coordinates_standardized = injection_coordinates;
for curr_coord = 1:length(injection_coordinates)
    
    target_hemisphere = sign(ccf_midline - str_depths_query(curr_coord,3));
    injection_coords_ml_offset = abs(injection_coordinates{curr_coord}(:,3) - ccf_midline);
    injection_coordinates_hemisphere = sign(injection_coordinates{curr_coord}(:,3) - ccf_midline);
    
    injection_coords_ipsi = injection_coordinates_hemisphere == target_hemisphere;
    injection_coords_contra = injection_coordinates_hemisphere == -target_hemisphere;
    
    injection_coordinates_standardized{curr_coord}(injection_coords_ipsi,3) = ...
        ccf_midline + injection_coords_ml_offset(injection_coords_ipsi);
    injection_coordinates_standardized{curr_coord}(injection_coords_contra,3) = ...
        ccf_midline - injection_coords_ml_offset(injection_coords_contra);
    
end

% Combine side-standardized injection data across hemispheres
injection_coordinates_bilateral = arrayfun(@(x) ...
    vertcat(injection_coordinates_standardized{x}, ...
    injection_coordinates_standardized{n_aligned_depths+x}), ...
    1:n_aligned_depths,'uni',false);

% Convert points from CCF to widefield
ccf_tform_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\ccf_tform'];
load(ccf_tform_fn);

um2pixel = 20.6;
injection_coordinates_bilateral_wf = cellfun(@(x) ...
    [x(:,[3,1]).*(1/(um2pixel)),ones(size(x,1),1)]*ccf_tform.T, ...
    injection_coordinates_bilateral,'uni',false);

% Get projection density
projection_strength = {injection_parameters.density};

projection_strength_bilateral = arrayfun(@(x) ...
    vertcat(projection_strength{x}', ...
    projection_strength{n_aligned_depths+x}'), ...
    1:n_aligned_depths,'uni',false);

% Load kernel templates for overlay
n_aligned_depths = 4;
kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template_' num2str(n_aligned_depths) '_depths.mat'];
load(kernel_template_fn);

figure; 
colormap(brewermap([],'*RdBu'));
for curr_depth = 1:n_aligned_depths
    subplot(1,n_aligned_depths,curr_depth); hold on; axis image off;
    set(gca,'YDir','reverse');
    
%     imagesc(kernel_template(:,:,curr_depth));
%     caxis([-prctile(abs(kernel_template(:)),99),prctile(abs(kernel_template(:)),99)]);
    
    scatter(injection_coordinates_wf_bilateral{curr_depth}(:,1), ...
        injection_coordinates_wf_bilateral{curr_depth}(:,2), ...
        projection_strength_bilateral{curr_depth}*50 + 10, ...
        'k','filled');
    
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
end

% Color all injection sites and plot overlaid
jet_basic = jet(255);
dark_colors = max(jet_basic,[],2) ~= 1;
jet_alt = max(interp1(1:255,jet_basic,linspace(find(~dark_colors,1,'first'), ...
    find(~dark_colors,1,'last'),255)) - 0.2,0);

cmap = jet_alt;
% depth_color_idx = round(conv(linspace(1,size(cmap,1),n_aligned_depths+1),[0.5,0.5],'valid'));
depth_color_idx = round(linspace(1,size(cmap,1),n_aligned_depths));
plot_colors = cmap(depth_color_idx,:);

figure;
hold on; set(gca,'YDir','reverse'); axis image off;
AP_reference_outline('ccf_aligned','k');
for curr_depth = 1:n_aligned_depths
    % (plot points from both hemispheres)
    scatter(injection_coordinates_wf_bilateral{curr_depth}(:,1), ...
        injection_coordinates_wf_bilateral{curr_depth}(:,2), ...
        projection_strength_bilateral{curr_depth}*100 + 10, ...
        plot_colors(curr_depth,:),'filled');
end


%% Fig 2a: Task performance and reaction time

% Load behavior
bhv_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\bhv_processing\bhv.mat';
load(bhv_fn);

% Exclude bad behavior sessions
exclude_data = true;

bhv_fieldnames = fieldnames(bhv);
experiment_fields = cellfun(@(curr_field) ...
    length([bhv.(curr_field)]) == length([bhv.days]),bhv_fieldnames);

% Load pre-marked experiments to exclude and cut out bad ones
if exclude_data
    exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
    exclude_fn = 'bhv_use_experiments';
    load([exclude_path filesep exclude_fn]);
    
    % Pull out used experiments for the animals loaded
    use_experiments_animals = ismember({bhv.animal},{bhv_use_experiments.animals});
    use_experiments = {bhv_use_experiments(use_experiments_animals).use_experiments}';
    
    % Cut out bad experiments for any experiment data fields
    for curr_field = bhv_fieldnames(experiment_fields)'
        for curr_animal = 1:length(use_experiments)
            bhv(curr_animal).(cell2mat(curr_field)) = ...
                bhv(curr_animal).(cell2mat(curr_field))(use_experiments{curr_animal});
        end
    end
end

conditions = unique(vertcat(bhv.conditions),'rows');
trial_choice_cat = arrayfun(@(x) horzcat(bhv(x).trial_choice{:}),1:length(bhv),'uni',false);
trial_outcome_cat = arrayfun(@(x) horzcat(bhv(x).trial_outcome{:}),1:length(bhv),'uni',false);
trial_side_cat = arrayfun(@(x) horzcat(bhv(x).trial_side{:}),1:length(bhv),'uni',false);
trial_contrast_cat = arrayfun(@(x) horzcat(bhv(x).trial_contrast{:}),1:length(bhv),'uni',false);
trial_condition_cat = cellfun(@(side,contrast) side.*contrast,trial_side_cat,trial_contrast_cat,'uni',false);
trial_wheel_velocity_cat = arrayfun(@(x) vertcat(bhv(x).trial_wheel_velocity{:})',1:length(bhv),'uni',false);
stim_to_move_cat = arrayfun(@(x) horzcat(bhv(x).stim_to_move{:}),1:length(bhv),'uni',false);
stim_to_feedback_cat = arrayfun(@(x) horzcat(bhv(x).stim_to_feedback{:}),1:length(bhv),'uni',false);

% Distinguish early/late movements
go_time = 0.5;
trial_timing = arrayfun(@(animal) cellfun(@(x) 1+(x > go_time), ...
    bhv(animal).stim_to_move,'uni',false),1:length(bhv),'uni',false);
trial_timing_cat = arrayfun(@(animal) ...
    horzcat(trial_timing{animal}{:}),1:length(bhv),'uni',false);

% Plot psychometric 
frac_left = cell2mat(cellfun(@(choice,condition) ...
    grpstats(choice == -1,condition),trial_choice_cat,trial_condition_cat,'uni',false));

figure; hold on; axis square;
plot(conditions,frac_left,'linewidth',2,'color',[0.5,0.5,0.5]);
plot(conditions,nanmean(frac_left,2),'linewidth',5,'color','k');
xlim([-1,1]);
ylim([0,1]);
line([0,0],ylim,'linestyle','--','color','k');
line(xlim,[0.5,0.5],'linestyle','--','color','k');
xlabel('Stimulus side*contrast');
ylabel('Fraction go left');

% Plot psychometric split by early/late movements
frac_left_earlymove = cell2mat(cellfun(@(choice,condition,stim_to_move) ...
    grpstats(choice(stim_to_move < 0.5) == -1,condition(stim_to_move < 0.5)), ...
    trial_choice_cat,trial_condition_cat,stim_to_move_cat,'uni',false));

frac_left_latemove = cell2mat(cellfun(@(choice,condition,stim_to_move) ...
    grpstats(choice(stim_to_move >= 0.5) == -1,condition(stim_to_move >= 0.5)), ...
    trial_choice_cat,trial_condition_cat,stim_to_move_cat,'uni',false));

figure;

subplot(1,3,1); hold on; axis square;
plot(conditions,frac_left,'linewidth',2,'color',[0.5,0.5,0.5]);
plot(conditions,nanmean(frac_left,2),'linewidth',5,'color','k');
xlim([-1,1]);
ylim([0,1]);
line([0,0],ylim,'linestyle','--','color','k');
line(xlim,[0.5,0.5],'linestyle','--','color','k');
xlabel('Condition');
ylabel('Fraction go left');
title('All trials');

subplot(1,3,2); hold on; axis square;
plot(conditions,frac_left_earlymove,'linewidth',2,'color',[0.5,0.5,0.5]);
plot(conditions,nanmean(frac_left_earlymove,2),'linewidth',5,'color','k');
xlim([-1,1]);
ylim([0,1]);
line([0,0],ylim,'linestyle','--','color','k');
line(xlim,[0.5,0.5],'linestyle','--','color','k');
xlabel('Condition');
ylabel('Fraction go left');
title('Early move');

subplot(1,3,3); hold on; axis square;
plot(conditions,frac_left_latemove,'linewidth',2,'color',[0.5,0.5,0.5]);
plot(conditions,nanmean(frac_left_latemove,2),'linewidth',5,'color','k');
xlim([-1,1]);
ylim([0,1]);
line([0,0],ylim,'linestyle','--','color','k');
line(xlim,[0.5,0.5],'linestyle','--','color','k');
xlabel('Condition');
ylabel('Fraction go left');
title('Late move');

% Plot distribution of reaction times
move_time_bins = 0:0.01:1;
day_split = 4;

move_time_centers = move_time_bins(1:end-1) + diff(move_time_bins)/2;
stim_to_move_binned = zeros(day_split,length(move_time_bins)-1,length(bhv));

for curr_animal = 1:length(bhv)
    for curr_day = 1:length(bhv(curr_animal).stim_to_move)
        curr_data = bhv(curr_animal).stim_to_move{curr_day};
        trials_split = round(linspace(1,length(curr_data),day_split+1));
        for curr_split = 1:day_split
            stim_to_move_binned(curr_split,:,curr_animal) = ...
                stim_to_move_binned(curr_split,:,curr_animal) + ...
                histcounts(curr_data(trials_split(curr_split): ...
                trials_split(curr_split+1)),move_time_bins);
        end
    end
end
stim_to_move_binned_norm = bsxfun(@rdivide,stim_to_move_binned, ...
    sum(sum(stim_to_move_binned,1),2));

figure; hold on;
for curr_animal = 1:length(bhv)
    AP_stackplot(stim_to_move_binned_norm(:,:,curr_animal)', ...
        move_time_centers,0.02,false,[0.5,0.5,0.5],1:day_split);
end
AP_stackplot(nanmean(stim_to_move_binned_norm,3)', ...
    move_time_centers,0.02,false,'k',1:day_split);
xlabel('Time to movement onset')
ylabel('Frequency by fraction within day')
axis tight
line([0.5,0.5],ylim,'linestyle','--','color','k');


%% Figure 2b-d: Striatal activity, task regression

% Load data
data_fn = 'trial_activity_choiceworld';
exclude_data = true;
AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(mua_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_animal;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));

% Plot average stimulus-aligned activity in striatum
plot_trials = trial_contrast_allcat > 0 & trial_side_allcat == 1 & trial_choice_allcat == -1;
plot_trials_exp = mat2cell(plot_trials,use_split,1);

mua_allcat_exp = mat2cell(mua_allcat,use_split,length(t),n_depths);
mua_allcat_exp_mean = cell2mat(permute(cellfun(@(x,trials) ...
    permute(nanmean(x(trials,:,:),1),[3,2,1]),mua_allcat_exp,plot_trials_exp,'uni',false),[2,3,1]));

figure;
p = nan(n_depths,1);
for curr_depth = 1:n_depths
    p(curr_depth) = subplot(n_depths,1,curr_depth);
    AP_errorfill(t,nanmean(mua_allcat_exp_mean(curr_depth,:,:),3), ...
        AP_sem(mua_allcat_exp_mean(curr_depth,:,:),3),'k',0.5)
    xlabel('Time from stimulus');
    line([0,0],ylim,'color','k');
end
linkaxes(p);

% Get task>striatum parameters
regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
n_regressors = length(regressor_labels);
t_shifts = {[0,0.5]; ... % stim
    [-0.5,1]; ... % move
    [-0.1,0.5]; ... % go cue
    [-0.5,1]}; % outcome

sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
    round(x(2)*(sample_rate)),t_shifts,'uni',false);
t_shifts = cellfun(@(x) x/sample_rate,sample_shifts,'uni',false);

% Normalize task>striatum kernels (ridiculous loop makes it clearer to read)
mua_taskpred_k_all_norm = {};
for curr_animal = 1:length(mua_taskpred_k_all)
    for curr_exp = 1:length(mua_taskpred_k_all{curr_animal})
        for curr_regressor = 1:n_regressors
            for curr_depth = 1:n_depths
                mua_taskpred_k_all_norm{curr_animal}{curr_exp}{curr_regressor,curr_depth} = ...
                    mua_taskpred_k_all{curr_animal}{curr_exp}{curr_regressor,curr_depth}./ ...
                    ((1/sample_rate)*mua_norm{curr_animal}{curr_exp}(curr_depth));
            end
        end
    end
end

% Average and concatenate task>striatum kernels within animals
task_str_k_animal = cell(n_regressors,n_depths);
for curr_animal = 1:length(mua_taskpred_k_all_norm)
    if isempty(mua_taskpred_k_all_norm{curr_animal})
        continue
    end
    curr_k = cat(3,mua_taskpred_k_all_norm{curr_animal}{:});
    for curr_depth = 1:n_depths
        for curr_regressor = 1:n_regressors
            curr_k_mean = nanmean(cat(3,curr_k{curr_regressor,curr_depth,:}),3);
            task_str_k_animal{curr_regressor,curr_depth} = cat(3, ...
                task_str_k_animal{curr_regressor,curr_depth},curr_k_mean);
        end
    end
end

% Plot task>striatum kernels
figure;
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors
        if curr_regressor == 1
            n_subregressors = 10;
            col = colormap_BlueWhiteRed(n_subregressors/2);
            col(6,:) = [];
        else
            n_subregressors = 2;
            col = lines(n_subregressors);
        end
        
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);  
        
        curr_kernels = task_str_k_animal{curr_regressor,curr_depth};
        for curr_subregressor = 1:n_subregressors
            AP_errorfill(t_shifts{curr_regressor}, ...
                nanmean(curr_kernels(curr_subregressor,:,:),3), ...
                AP_sem(curr_kernels(curr_subregressor,:,:),3), ...
                col(curr_subregressor,:),0.5);
        end
        
        axis off;
        line([0,0],ylim,'color','k');
        
    end
end

linkaxes(p);


% Plot task>striatum regression examples
figure;
for curr_depth = 1:n_depths   
    
    % Set current data (pad trials with NaNs for spacing)
    n_pad = 10;
    curr_data = padarray(mua_allcat(:,:,curr_depth),[0,n_pad],NaN,'post');
    curr_pred_data = padarray(mua_taskpred_allcat(:,:,curr_depth),[0,n_pad],NaN,'post');
    nan_samples = isnan(curr_data) | isnan(curr_pred_data);

    % Smooth
    smooth_filt = ones(1,3)/3;
    curr_data = conv2(curr_data,smooth_filt,'same');
    curr_pred_data = conv2(curr_pred_data,smooth_filt,'same');
    
    % Set common NaNs for R^2    
    curr_data_nonan = curr_data; 
    curr_data_nonan(nan_samples) = NaN;
    
    curr_pred_data_nonan = curr_pred_data; 
    curr_pred_data(nan_samples) = NaN; 
    
    % Get squared error for each trial
    trial_r2 = 1 - (nansum((curr_data_nonan-curr_pred_data_nonan).^2,2)./ ...
        nansum((curr_data_nonan-nanmean(curr_data_nonan,2)).^2,2));
    
    [~,trial_r2_rank] = sort(trial_r2);
    
    plot_prctiles = round(prctile(1:length(trial_r2),linspace(25,75,20)));
    plot_trials = trial_r2_rank(plot_prctiles);
    
    subplot(n_depths,1,curr_depth); hold on;
    plot(reshape(curr_data(plot_trials,:)',[],1),'k');
    plot(reshape(curr_pred_data(plot_trials,:)',[],1),'b');
    
end

% Get R^2 for task regression 
taskpred_r2 = nan(max(split_idx),n_depths);
for curr_exp = 1:max(split_idx)
    curr_data = reshape(permute(mua_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_pred_data = reshape(permute(mua_taskpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_pred_data);
    curr_data(nan_samples) = NaN;
    curr_pred_data(nan_samples) = NaN;

    taskpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_pred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
end
figure;
errorbar(nanmean(taskpred_r2,1),AP_sem(taskpred_r2,1),'k','linewidth',2);
xlabel('Striatum depth');
ylabel('Task explained variance');


%% Fig 2e: Interaction between Str 2 sensory and motor activity

% Load data
data_fn = 'trial_activity_choiceworld';
exclude_data = true;
AP_load_concat_normalize_ctx_str;

% Depth to plot
plot_depth = 2;

% Get split index
trials_allcat = size(mua_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_animal;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));

% Set number of shuffles for significance testing
n_shuff = 10000; 

%%% STIM ACTIVITY

% Set time to average activity
use_stim_t = t > 0.05 & t < 0.15;

% Get stim bins for each trial
stims = unique([0.06,0.125,0.25,0.5,1].*[-1;1]);
stim_bins = trial_contrast_allcat.*trial_side_allcat;
stim_bins_idx = discretize(stim_bins,stims);

% Set contexts for activity
stim_contexts = [trial_choice_allcat == -1, ...
    trial_choice_allcat == 1];

% Get stim-isolated activity by trial
stim_activity = mua_allcat - mua_taskpred_reduced_allcat(:,:,:,1);
stim_trial_activity = nanmean(stim_activity(:,use_stim_t,:),2);

% Split activity by animal, stim, and correct/incorrect
stim_bins = unique(trial_contrast_allcat(trial_contrast_allcat ~= 0).*[-1,1]);
stim_trial_activity_split = cell(max(split_idx),length(stim_bins),size(stim_contexts,2),n_depths);
for curr_exp = 1:max(split_idx)
    for curr_bin = 1:length(stim_bins)
        for curr_context = 1:size(stim_contexts,2)               
            for curr_depth = 1:n_depths
                % Get activity of all group trials (exclude NaNs)
                curr_trials = ...
                    split_idx == curr_exp & ...
                    trial_contrast_allcat.*trial_side_allcat == stim_bins(curr_bin) & ...
                    stim_contexts(:,curr_context);             
                curr_activity = stim_trial_activity(curr_trials,:,curr_depth);
                stim_trial_activity_split{curr_exp,curr_bin,curr_context,curr_depth} = ...
                    curr_activity(~isnan(curr_activity));                
            end
            
        end
    end
end

stim_trial_activity_split_mean = cellfun(@nanmean,stim_trial_activity_split);

% Plot stim activity by movement direction
figure('Name',['Str ' num2str(plot_depth)]);
subplot(1,2,1); hold on;
line_col = [0.6,0,0.6;0,0.6,0];

dot_col = colormap_BlueWhiteRed(5);
dot_col(6,:) = [];

p = nan(size(stim_contexts,2),1);
for curr_context = 1:size(stim_contexts,2)
    p(curr_context) = ...
        errorbar(stims,nanmean(stim_trial_activity_split_mean(:,:,curr_context,plot_depth),1), ...
        AP_sem(stim_trial_activity_split_mean(:,:,curr_context,plot_depth),1),'color',line_col(curr_context,:),'linewidth',3);
    scatter(stims,nanmean(stim_trial_activity_split_mean(:,:,curr_context,plot_depth),1),80,dot_col, ...
        'Filled','MarkerEdgeColor',line_col(curr_context,:),'linewidth',3);
end
xlabel('Stimulus');
ylabel('Activity');
legend(p,{'Move left','Move right'});

% Get condition difference and compare to shuffle
subplot(1,2,2); hold on;

stim_condition_diff = nanmean(stim_trial_activity_split_mean(:,:,1,plot_depth) - ...
    stim_trial_activity_split_mean(:,:,2,plot_depth),1);

stim_condition_shuff_diff = nan(max(split_idx),length(stim_bins),n_shuff);
for curr_exp = 1:max(split_idx)
    for curr_bin = 1:length(stim_bins)
        curr_act = stim_trial_activity_split(curr_exp,curr_bin,:,plot_depth);
        % Shuffle, split, difference
        curr_act_shuff = mat2cell(AP_shake(repmat(vertcat( ...
            curr_act{:}),1,n_shuff),1),cellfun(@length,curr_act),n_shuff);
        stim_condition_shuff_diff(curr_exp,curr_bin,:) = ...
            nanmean(curr_act_shuff{1},1) - nanmean(curr_act_shuff{2},1);            
    end
end
stim_condition_diff_ci = squeeze(prctile(nanmean(stim_condition_shuff_diff,1),[2.5,97.5],3));

plot(stims,stim_condition_diff,'k','linewidth',2);
plot(stims,stim_condition_diff_ci,'k','linewidth',2,'linestyle','--');
xlabel('Stimulus');
ylabel('Condition difference');


%%% MOVEMENT ACTIVITY

% Set time to average activity
use_move_t = t > -0.05 & t < 0.05;

% Get velocity bins for each trial
max_vel = AP_signed_max(wheel_velocity_allcat_move(:,t > 0 & t < 0.2),2);

% Normalize velocity for each day
max_vel_exp = mat2cell(max_vel,trials_recording,1);
max_vel_norm = cell2mat(cellfun(@(x) x./nanmean(abs(x)),max_vel_exp,'uni',false));

n_vel_bins = 4;
vel_edges = prctile(abs(max_vel_norm),linspace(0,100,n_vel_bins+1));
vel_edges = linspace(min(abs(max_vel_norm)),max(abs(max_vel_norm)),n_vel_bins+1);

vel_edges = sort([vel_edges,-vel_edges]);
vel_centers = vel_edges(1:end-1) + diff(vel_edges)/2;

vel_bins = discretize(max_vel_norm,vel_edges);

% Set contexts for activity
move_contexts = [trial_side_allcat == 1 & trial_contrast_allcat > 0, ...
    trial_side_allcat == -1 & trial_contrast_allcat > 0, ...
    trial_contrast_allcat == 0];

% Get move-isolated activity by trial
move_activity = mua_allcat_move - mua_taskpred_reduced_allcat_move(:,:,:,2);
move_trial_activity = nanmean(move_activity(:,use_move_t,:),2);

% Split activity by animal, velocity, and stim side/zero contrast
move_trial_activity_split = cell(max(split_idx),max(vel_bins),size(move_contexts,2),n_depths);
for curr_exp = 1:max(split_idx)
    for curr_bin = 1:max(vel_bins)
        for curr_context = 1:size(move_contexts,2)               
            for curr_depth = 1:n_depths
                % Get activity of all group trials (exclude NaNs)
                curr_trials = ...
                    split_idx == curr_exp & ...
                    vel_bins == curr_bin & ...
                    move_contexts(:,curr_context);              
                curr_activity = move_trial_activity(curr_trials,:,curr_depth);
                move_trial_activity_split{curr_exp,curr_bin,curr_context,curr_depth} = ...
                    curr_activity(~isnan(curr_activity));                
            end
            
        end
    end
end

move_trial_activity_split_mean = cellfun(@nanmean,move_trial_activity_split);

% Plot move activity by stim side
figure('Name',['Str ' num2str(plot_depth)]);
subplot(1,2,1); hold on;
line_col = [0.7,0,0;0,0,0.7;0.5,0.5,0.5];

dot_col = [linspace(0.2,1,n_vel_bins)',0.1*ones(n_vel_bins,1),linspace(0.2,1,n_vel_bins)'; ...
    1,0,0;
    0.1*ones(n_vel_bins,1),linspace(1,0.2,n_vel_bins)',0.1*ones(n_vel_bins,1)];

p = nan(size(move_contexts,2),1);
for curr_context = 1:size(move_contexts,2)
    p(curr_context) = ...
        errorbar(vel_centers,nanmean(move_trial_activity_split_mean(:,:,curr_context,plot_depth),1), ...
        AP_sem(move_trial_activity_split_mean(:,:,curr_context,plot_depth),1),'color',line_col(curr_context,:),'linewidth',3);
    scatter(vel_centers,nanmean(move_trial_activity_split_mean(:,:,curr_context,plot_depth),1),80,dot_col, ...
        'Filled','MarkerEdgeColor',line_col(curr_context,:),'linewidth',3);
end
xlabel('Velocity');
ylabel('Activity');
legend(p,{'Stim right','Stim left','No stim'});
 
% Get condition difference and compare to shuffle
subplot(1,2,2); hold on;

move_condition_diff_1 = nanmean(move_trial_activity_split_mean(:,:,1,plot_depth) - ...
    move_trial_activity_split_mean(:,:,2,plot_depth),1);
move_condition_diff_2 = nanmean(move_trial_activity_split_mean(:,:,2,plot_depth) - ...
    move_trial_activity_split_mean(:,:,3,plot_depth),1);
move_condition_diff_3 = nanmean(move_trial_activity_split_mean(:,:,1,plot_depth) - ...
    move_trial_activity_split_mean(:,:,3,plot_depth),1);

move_condition_shuff_diff_1 = nan(max(split_idx),max(vel_bins),n_shuff);
move_condition_shuff_diff_2 = nan(max(split_idx),max(vel_bins),n_shuff);
move_condition_shuff_diff_3 = nan(max(split_idx),max(vel_bins),n_shuff);
for curr_exp = 1:max(split_idx)
    for curr_bin = 1:max(vel_bins)
        curr_act = move_trial_activity_split(curr_exp,curr_bin,:,plot_depth);
        % Shuffle, split, difference
        curr_act_shuff = mat2cell(AP_shake(repmat(vertcat( ...
            curr_act{:}),1,n_shuff),1),cellfun(@length,curr_act),n_shuff);
        move_condition_shuff_diff_1(curr_exp,curr_bin,:) = ...
            nanmean(curr_act_shuff{1},1) - nanmean(curr_act_shuff{2},1);            
        move_condition_shuff_diff_2(curr_exp,curr_bin,:) = ...
            nanmean(curr_act_shuff{2},1) - nanmean(curr_act_shuff{3},1);            
        move_condition_shuff_diff_3(curr_exp,curr_bin,:) = ...
            nanmean(curr_act_shuff{1},1) - nanmean(curr_act_shuff{3},1);            
    end
end
move_condition_diff_1_ci = squeeze(prctile(nanmean(move_condition_shuff_diff_1,1),[2.5,97.5],3));
move_condition_diff_2_ci = squeeze(prctile(nanmean(move_condition_shuff_diff_2,1),[2.5,97.5],3));
move_condition_diff_3_ci = squeeze(prctile(nanmean(move_condition_shuff_diff_3,1),[2.5,97.5],3));

col = lines(3);
plot(vel_centers,move_condition_diff_1,'color',col(1,:),'linewidth',2);
plot(vel_centers,move_condition_diff_1_ci,'color',col(1,:),'linewidth',2,'linestyle','--');

plot(vel_centers,move_condition_diff_2,'color',col(2,:),'linewidth',2);
plot(vel_centers,move_condition_diff_2_ci,'color',col(2,:),'linewidth',2,'linestyle','--');

plot(vel_centers,move_condition_diff_3,'color',col(3,:),'linewidth',2);
plot(vel_centers,move_condition_diff_3_ci,'color',col(3,:),'linewidth',2,'linestyle','--');
xlabel('Velocity');
ylabel('Condition difference');


%% Fig 3a-d: Cortical activity, task>cortex regression, cortex>striatum regression

% Load data
data_fn = 'trial_activity_choiceworld';
exclude_data = true;
AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(mua_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_animal;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));


%%% Cortical activity

% Get average stim-aligned fluorescence 
plot_trials = trial_contrast_allcat > 0 & trial_side_allcat == 1 & trial_choice_allcat == -1;
plot_trials_exp = mat2cell(plot_trials,use_split,1);

fluor_allcat_deconv_exp = mat2cell(fluor_allcat_deconv,use_split,length(t),n_vs);
plot_px = nanmean(cell2mat(permute(cellfun(@(x,trials) svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    squeeze(nanmean(x(trials,:,:),1))'),fluor_allcat_deconv_exp,plot_trials_exp,'uni',false),[2,3,4,1])),4);

plot_t = [find(t > 0.07,1),find(t > 0.18,1),find(t > 0.3,1),find(t > 0.65,1)];
figure;
for curr_t = 1:length(plot_t)
   subplot(1,length(plot_t),curr_t);
   imagesc(plot_px(:,:,plot_t(curr_t)));
   AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
   axis image off;
   colormap(brewermap([],'Greens'));
   caxis([0,0.01]);
   title([num2str(t(plot_t(curr_t))) 's from stim']);
end

%%% Task>cortex regression results

% Get task>cortex parameters
n_regressors = 4;
t_shifts = {[0,0.5]; ... % stim
    [-0.5,1]; ... % move
    [-0.1,0.5]; ... % go cue
    [-0.5,1]}; % outcome
regressor_labels = {'Stim','Move onset','Go cue','Outcome'};

sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
    round(x(2)*(sample_rate)),t_shifts,'uni',false);
t_shifts = cellfun(@(x) x/sample_rate,sample_shifts,'uni',false);

% Get average task>striatum kernels
regressor_px = cell(n_regressors,1);
for curr_regressor = 1:n_regressors
    curr_k_cell = cellfun(@(x) x(curr_regressor,:),vertcat(fluor_taskpred_k_all{:}),'uni',false);
    curr_k_cell = vertcat(curr_k_cell{:});
    curr_k = permute(cell2mat(permute(arrayfun(@(x) ...
        nanmean(cat(3,curr_k_cell{:,x}),3),1:n_vs,'uni',false),[1,3,2])),[3,2,1]);
    curr_k_px = cell2mat(permute(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:size(curr_k,1)), ...
        curr_k(:,:,x)),1:size(curr_k,3),'uni',false),[1,3,4,2]));
    AP_image_scroll(curr_k_px,t_shifts{curr_regressor});
    axis image;
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'*RdBu'));
    AP_reference_outline('ccf_aligned','k');
    set(gcf,'Name',regressor_labels{curr_regressor});
    
    regressor_px{curr_regressor} = curr_k_px;
end

% Plot max for each kernel across time
regressor_t_max = cellfun(@(x) squeeze(max(x,[],3)),regressor_px,'uni',false);

figure;
max_subregressors = max(cellfun(@(x) size(x,3),regressor_t_max));
max_c = max(abs(cell2mat(cellfun(@(x) x(:),regressor_t_max,'uni',false))));
for curr_regressor = 1:n_regressors
    for curr_subregressor = 1:size(regressor_t_max{curr_regressor},3)
        subplot(n_regressors,max_subregressors, ...
            curr_subregressor+(curr_regressor-1)*max_subregressors);
        imagesc(regressor_t_max{curr_regressor}(:,:,curr_subregressor));
        AP_reference_outline('ccf_aligned','k');
        axis image off; 
        colormap(brewermap([],'Greens'));
        caxis([0,max_c]);
    end
end

%%% Cortex>striatum regression results

% Plot cortex>striatum regression examples
figure;
for curr_depth = 1:n_depths   
    
     % Set current data (pad trials with NaNs for spacing)
    n_pad = 10;
    curr_data = padarray(mua_allcat(:,:,curr_depth),[0,n_pad],NaN,'post');
    curr_pred_data = padarray(mua_ctxpred_allcat(:,:,curr_depth),[0,n_pad],NaN,'post');
    nan_samples = isnan(curr_data) | isnan(curr_pred_data);

    % Smooth
    smooth_filt = ones(1,3)/3;
    curr_data = conv2(curr_data,smooth_filt,'same');
    curr_pred_data = conv2(curr_pred_data,smooth_filt,'same');
    
    % Set common NaNs for R^2    
    curr_data_nonan = curr_data; 
    curr_data_nonan(nan_samples) = NaN;
    
    curr_pred_data_nonan = curr_pred_data; 
    curr_pred_data(nan_samples) = NaN; 
    
    % Get squared error for each trial
    trial_r2 = 1 - (nansum((curr_data_nonan-curr_pred_data_nonan).^2,2)./ ...
        nansum((curr_data_nonan-nanmean(curr_data_nonan,2)).^2,2));
    
    [~,trial_r2_rank] = sort(trial_r2);
    
    plot_prctiles = round(prctile(1:length(trial_r2),linspace(25,75,20)));
    plot_trials = trial_r2_rank(plot_prctiles);
    
    subplot(n_depths,1,curr_depth); hold on;
    plot(reshape(curr_data(plot_trials,:)',[],1),'k');
    plot(reshape(curr_pred_data(plot_trials,:)',[],1),'b');
    
end

% Get R^2 for task regression 
taskpred_r2 = nan(max(split_idx),n_depths);
for curr_exp = 1:max(split_idx)
    curr_data = reshape(permute(mua_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_pred_data = reshape(permute(mua_ctxpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_pred_data);
    curr_data(nan_samples) = NaN;
    curr_pred_data(nan_samples) = NaN;

    taskpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_pred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
end
figure;
errorbar(nanmean(taskpred_r2,1),AP_sem(taskpred_r2,1),'k','linewidth',2);
xlabel('Striatum depth');
ylabel('Cortex explained variance');


%%% Task>cortex>striatum regression results

% Get task>striatum parameters
regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
n_regressors = length(regressor_labels);
t_shifts = {[0,0.5]; ... % stim
    [-0.5,1]; ... % move
    [-0.1,0.5]; ... % go cue
    [-0.5,1]}; % outcome

sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
    round(x(2)*(sample_rate)),t_shifts,'uni',false);
t_shifts = cellfun(@(x) x/sample_rate,sample_shifts,'uni',false);

% Normalize task kernels (ridiculous loop makes it clearer to read)
mua_taskpred_k_all_norm = {};
mua_ctxpred_taskpred_k_all_norm = {};
for curr_animal = 1:length(mua_ctxpred_taskpred_k_all)
    for curr_exp = 1:length(mua_ctxpred_taskpred_k_all{curr_animal})
        for curr_regressor = 1:n_regressors
            for curr_depth = 1:n_depths
                mua_taskpred_k_all_norm{curr_animal}{curr_exp}{curr_regressor,curr_depth} = ...
                    mua_taskpred_k_all{curr_animal}{curr_exp}{curr_regressor,curr_depth}./ ...
                    ((1/sample_rate)*mua_norm{curr_animal}{curr_exp}(curr_depth));
                mua_ctxpred_taskpred_k_all_norm{curr_animal}{curr_exp}{curr_regressor,curr_depth} = ...
                    mua_ctxpred_taskpred_k_all{curr_animal}{curr_exp}{curr_regressor,curr_depth}./ ...
                    ((1/sample_rate)*mua_norm{curr_animal}{curr_exp}(curr_depth));
            end
        end
    end
end

% Average and concatenate task kernels within animals
task_str_k_animal = cell(n_regressors,n_depths);
task_ctx_str_k_animal = cell(n_regressors,n_depths);
for curr_animal = 1:length(mua_ctxpred_taskpred_k_all_norm)
    if isempty(mua_taskpred_k_all_norm{curr_animal})
        continue
    end
    curr_str_k = cat(3,mua_taskpred_k_all_norm{curr_animal}{:});
    curr_ctx_str_k = cat(3,mua_ctxpred_taskpred_k_all_norm{curr_animal}{:});
    for curr_depth = 1:n_depths
        for curr_regressor = 1:n_regressors
            curr_str_k_mean = nanmean(cat(3,curr_str_k{curr_regressor,curr_depth,:}),3);
            task_str_k_animal{curr_regressor,curr_depth} = cat(3, ...
                task_str_k_animal{curr_regressor,curr_depth},curr_str_k_mean);
            
            curr_ctx_str_k_mean = nanmean(cat(3,curr_ctx_str_k{curr_regressor,curr_depth,:}),3);
            task_ctx_str_k_animal{curr_regressor,curr_depth} = cat(3, ...
                task_ctx_str_k_animal{curr_regressor,curr_depth},curr_ctx_str_k_mean);
        end
    end
end

% Plot task>cortex>striatum kernels
figure;
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors
        if curr_regressor == 1
            n_subregressors = 10;
            col = colormap_BlueWhiteRed(n_subregressors/2);
            col(6,:) = [];
        else
            n_subregressors = 2;
            col = lines(n_subregressors);
        end
        
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);  
        
        curr_kernels = task_ctx_str_k_animal{curr_regressor,curr_depth};
        for curr_subregressor = 1:n_subregressors
            AP_errorfill(t_shifts{curr_regressor}, ...
                nanmean(curr_kernels(curr_subregressor,:,:),3), ...
                AP_sem(curr_kernels(curr_subregressor,:,:),3), ...
                col(curr_subregressor,:),0.5);
        end
        
        axis off;
        line([0,0],ylim,'color','k');
        
    end
end
linkaxes(p);


% Plot task>cortex>striatum - task>striatum kernels
figure;
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors
        if curr_regressor == 1
            n_subregressors = 10;
            col = colormap_BlueWhiteRed(n_subregressors/2);
            col(6,:) = [];
        else
            n_subregressors = 2;
            col = lines(n_subregressors);
        end
        
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);  
        
        curr_kernels = task_str_k_animal{curr_regressor,curr_depth} - ...
            task_ctx_str_k_animal{curr_regressor,curr_depth};
        for curr_subregressor = 1:n_subregressors
            AP_errorfill(t_shifts{curr_regressor}, ...
                nanmean(curr_kernels(curr_subregressor,:,:),3), ...
                AP_sem(curr_kernels(curr_subregressor,:,:),3), ...
                col(curr_subregressor,:),0.5);
        end
        
        axis off;
        line([0,0],ylim,'color','k');
        
    end
end
linkaxes(p);



%% Fig 4a: Differences in measured and cortex-predicted striatum preferences

% Load data
data_fn = 'trial_activity_choiceworld';
exclude_data = true;
AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(mua_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_animal;

% Stim/move/outcome rank differences by experiment
mua_exp = mat2cell(mua_allcat,use_split,length(t),n_depths);
mua_ctxpred_exp = mat2cell(mua_ctxpred_allcat,use_split,length(t),n_depths);

mua_taskpred_reduced_exp = mat2cell(mua_taskpred_reduced_allcat,use_split,length(t),n_depths, ...
    size(mua_taskpred_reduced_allcat,4));
mua_ctxpred_taskpred_reduced_exp = mat2cell(mua_ctxpred_taskpred_reduced_allcat,use_split,length(t),n_depths, ...
    size(mua_ctxpred_taskpred_reduced_allcat,4));

move_t_exp = mat2cell(move_t,use_split,1);
move_idx_exp = mat2cell(move_idx,use_split,1);
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
trial_side_allcat_exp = mat2cell(trial_side_allcat,use_split,1);
trial_contrast_allcat_exp = mat2cell(trial_contrast_allcat,use_split,1);
trial_choice_allcat_exp = mat2cell(trial_choice_allcat,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
trial_groups = {'Stim','Move onset','Outcome'};
t_groups = {[0.05,0.15],[-0.05,0.05],[0,0.1]};

act_rank_difference = nan(length(t),n_depths,length(use_split),length(trial_groups));
predicted_act_rank_difference = nan(length(t),n_depths,length(use_split),length(trial_groups));

act_rank_difference_trial = nan(n_depths,length(use_split),length(trial_groups));
predicted_act_rank_difference_trial = nan(n_depths,length(use_split),length(trial_groups));

n_shuff = 10000;
act_rank_difference_trial_shuff = nan(n_depths,length(use_split),length(trial_groups),n_shuff);
predicted_act_rank_difference_trial_shuff = nan(n_depths,length(use_split),length(trial_groups),n_shuff);
act_rank_difference_trial_predshuff = nan(n_depths,length(use_split),length(trial_groups),n_shuff);

for curr_exp = 1:length(mua_exp)
    for curr_group = 1:length(trial_groups)
        
        % Get MUA/predicted using reduced model
        curr_regressor_idx = strcmp(trial_groups{curr_group},regressor_labels);
        curr_mua =  mua_exp{curr_exp} - mua_taskpred_reduced_exp{curr_exp}(:,:,:,curr_regressor_idx);
        curr_mua_ctx = mua_ctxpred_exp{curr_exp} - mua_ctxpred_taskpred_reduced_exp{curr_exp}(:,:,:,curr_regressor_idx);
        
        % Skip if there's no data in this experiment
        if isempty(curr_mua)
            continue
        end
        
        % Set common NaNs
        nan_samples = isnan(curr_mua) | isnan(curr_mua_ctx);
        curr_mua(nan_samples) = NaN;
        curr_mua_ctx(nan_samples) = NaN;
        
        % (movement: align to move onset)
        if any(strfind(lower(trial_groups{curr_group}),'move'))
            t_leeway = -t(1);
            leeway_samples = round(t_leeway*(sample_rate));
            for i = 1:size(curr_mua,1)
                curr_mua(i,:,:) = circshift(curr_mua(i,:,:),-move_idx_exp{curr_exp}(i)+leeway_samples,2);
                curr_mua_ctx(i,:,:) = circshift(curr_mua_ctx(i,:,:),-move_idx_exp{curr_exp}(i)+leeway_samples,2);
            end
        end
        
        % (outcome: align to outcome)
        if any(strfind(lower(trial_groups{curr_group}),'outcome'))
            t_leeway = -t(1);
            leeway_samples = round(t_leeway*(sample_rate));
            for i = 1:size(curr_mua,1)
                curr_mua(i,:,:) = circshift(curr_mua(i,:,:),-outcome_idx_exp{curr_exp}(i)+leeway_samples,2);
                curr_mua_ctx(i,:,:) = circshift(curr_mua_ctx(i,:,:),-outcome_idx_exp{curr_exp}(i)+leeway_samples,2);
            end
        end
        
        % Set trials and grouping to use
        switch trial_groups{curr_group}
            case 'Stim'
                use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} < 0.5 & trial_contrast_allcat_exp{curr_exp} > 0;
                trial_group_1 = trial_side_allcat_exp{curr_exp} == 1;
                trial_group_2 = trial_side_allcat_exp{curr_exp} == -1';
            case 'Move onset'
                use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} < 0.5;
                trial_group_1 = trial_choice_allcat_exp{curr_exp} == -1;
                trial_group_2 = trial_choice_allcat_exp{curr_exp} == 1';
            case 'Outcome'
                use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} < 0.5;
                trial_group_1 = trial_outcome_allcat_exp{curr_exp} == 1;
                trial_group_2 = trial_outcome_allcat_exp{curr_exp} == -1;
        end
        
        act_rank = tiedrank(curr_mua(use_trials,:,:));
        predicted_act_rank = tiedrank(curr_mua_ctx(use_trials,:,:));
        
        act_rank_difference(:,:,curr_exp,curr_group) = squeeze(( ...
            nanmean(act_rank(trial_group_1(use_trials),:,:),1) - ...
            nanmean(act_rank(trial_group_2(use_trials),:,:),1))./max(act_rank,[],1));
        predicted_act_rank_difference(:,:,curr_exp,curr_group) = squeeze(( ...
            nanmean(predicted_act_rank(trial_group_1(use_trials),:,:),1) - ...
            nanmean(predicted_act_rank(trial_group_2(use_trials),:,:),1))./max(predicted_act_rank,[],1));
        
        use_t = t >= t_groups{curr_group}(1) & t <= t_groups{curr_group}(2);
        act_rank_trial = tiedrank(squeeze(nanmean(curr_mua(use_trials,use_t,:),2)));
        predicted_act_rank_trial = tiedrank(squeeze(nanmean(curr_mua_ctx(use_trials,use_t,:),2)));
        
        act_rank_difference_trial(:,curr_exp,curr_group) = ...
            (nanmean(act_rank_trial(trial_group_1(use_trials),:),1) - ...
            nanmean(act_rank_trial(trial_group_2(use_trials),:),1))./max(act_rank_trial,[],1);
        predicted_act_rank_difference_trial(:,curr_exp,curr_group) = ...
            (nanmean(predicted_act_rank_trial(trial_group_1(use_trials),:),1) - ...
            nanmean(predicted_act_rank_trial(trial_group_2(use_trials),:),1))./max(predicted_act_rank_trial,[],1);     
        
        % Shuffle for significance
        shuff_trials = (trial_group_1 | trial_group_2) & use_trials;
        
        for curr_shuff = 1:n_shuff
            % Shuffle for label significance
            trial_group_1_shuff = trial_group_1;
            trial_group_1_shuff(shuff_trials) = AP_shake(trial_group_1_shuff(shuff_trials));
            trial_group_2_shuff = trial_group_2;
            trial_group_2_shuff(shuff_trials) = AP_shake(trial_group_2_shuff(shuff_trials));
            
            act_rank_difference_trial_shuff(:,curr_exp,curr_group,curr_shuff) = ...
                (nanmean(act_rank_trial(trial_group_1_shuff(use_trials),:),1) - ...
                nanmean(act_rank_trial(trial_group_2_shuff(use_trials),:),1))./max(act_rank_trial,[],1);
            predicted_act_rank_difference_trial_shuff(:,curr_exp,curr_group,curr_shuff) = ...
                (nanmean(predicted_act_rank_trial(trial_group_1_shuff(use_trials),:),1) - ...
                nanmean(predicted_act_rank_trial(trial_group_2_shuff(use_trials),:),1))./max(predicted_act_rank_trial,[],1);
        end
        
        
        % Shuffle for measured/predicted difference
        % (build an n_shuff sized matrix of half/half, then shake)
        meas_pred_shuff = AP_shake(cat(3, ...
            repmat(act_rank_trial,1,1,n_shuff), ...
            repmat(predicted_act_rank_trial,1,1,n_shuff)),3);
        
        act_rank_difference_trial_predshuff(:,curr_exp,curr_group,:) = ...
            permute((nanmean(meas_pred_shuff(trial_group_1(use_trials),:,1:n_shuff),1) - ...
            nanmean(meas_pred_shuff(trial_group_2(use_trials),:,1:n_shuff),1))./max(act_rank_trial,[],1) - ...
             (nanmean(meas_pred_shuff(trial_group_1(use_trials),:,n_shuff+1:end),1) - ...
            nanmean(meas_pred_shuff(trial_group_2(use_trials),:,n_shuff+1:end),1))./max(act_rank_trial,[],1),[2,1,3,4]);
        
    end
end

act_rank_difference_mean = squeeze(nanmean(act_rank_difference,3));
predicted_act_rank_difference_mean = squeeze(nanmean(predicted_act_rank_difference,3));

figure;
for curr_group = 1:length(trial_groups)
    p1 = subplot(2,length(trial_groups),curr_group);
    hold on; set(gca,'ColorOrder',copper(n_depths));
    plot(t,act_rank_difference_mean(:,:,curr_group),'linewidth',2);
    axis tight; line([0,0],ylim,'color','k');
    xlabel('Time');
    ylabel('Mean rank difference');
    title([trial_groups{curr_group} ': Measured']);
    
    p2 = subplot(2,length(trial_groups),length(trial_groups)+curr_group);
    hold on; set(gca,'ColorOrder',copper(n_depths));
    plot(t,predicted_act_rank_difference_mean(:,:,curr_group),'linewidth',2);
    axis tight; line([0,0],ylim,'color','k');
    xlabel('Time');
    ylabel('Mean rank difference');
    title([trial_groups{curr_group} ': Predicted']);
    linkaxes([p1,p2],'xy');
end

figure;
p = nan(length(trial_groups),1);
for curr_group = 1:length(trial_groups)
    p(curr_group) = subplot(1,3,curr_group); hold on;
    errorbar(squeeze(nanmean(act_rank_difference_trial(:,:,curr_group),2)), ...
        squeeze(AP_sem(act_rank_difference_trial(:,:,curr_group),2)),'linewidth',2,'color','k');
    errorbar(squeeze(nanmean(predicted_act_rank_difference_trial(:,:,curr_group),2)), ...
        squeeze(AP_sem(predicted_act_rank_difference_trial(:,:,curr_group),2)),'linewidth',2,'color',[0,0.7,0]);
    xlabel('Striatum depth');
    ylabel('Rank difference');
    legend({'Measured','Predicted'});
    title(trial_groups{curr_group});
end
linkaxes(p);

% Get significance from shuffled distribution
trained_predicted_diff_ci = prctile(squeeze(nanmean(act_rank_difference_trial_predshuff,2)),[2.5,97.5],3);
figure; hold on;
set(gca,'ColorOrder',lines(3));
plot(squeeze(nanmean(act_rank_difference_trial-predicted_act_rank_difference_trial,2)),'linewidth',2);
plot(reshape(trained_predicted_diff_ci,n_depths,[]),'linewidth',2,'linestyle','--');
xlabel('Striatum depth');
ylabel('Measured-predicted');



%% Fig S1b: widefield deconvolution
% use AP_deconv_wf_kernelfit to plot examples

load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\gcamp_kernel\gcamp6s_kernel.mat');
gcamp6s_kernel_cat = vertcat(gcamp6s_kernel.regression{:});
gcamp6s_kernel_norm = gcamp6s_kernel_cat./max(abs(gcamp6s_kernel_cat),[],2);
gcamp6s_kernel_mean = nanmean(gcamp6s_kernel_norm,1);

figure; hold on;
plot(fliplr(gcamp6s_kernel.regression_t),gcamp6s_kernel_norm','color',[0.5,0.5,0.5]);
plot(fliplr(gcamp6s_kernel.regression_t),gcamp6s_kernel_mean,'color','k','linewidth',2);
axis tight;
line([0,0],ylim,'color','k','linestyle','--');
line(xlim,[0,0],'color','k','linestyle','--');

xlabel('Lag from fluorescence to spikes');
ylabel('Max-normalized weight');



%% Fig S2b: plot template kernels 

% Load and plot the kernel templates
n_aligned_depths = 4;
kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template_' num2str(n_aligned_depths) '_depths.mat'];
load(kernel_template_fn);
n_kernels = n_aligned_depths;
figure; 
for i = 1:n_kernels
    subplot(n_kernels,1,i);
    imagesc(kernel_template(:,:,i));
    axis image off;
    caxis([-prctile(abs(kernel_template(:)),99),prctile(abs(kernel_template(:)),99)]);
    colormap(brewermap([],'*RdBu'));
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
end



%% Fig S2c: plot cortex>striatum kernels/projections aligned to histology
% use AP_ctx2str_probe to plot examples



%% ~~~~~~~~~~~~~~~~ TESTING BELOW ~~~~~~~~~~~~~~~~~~~~~~~~~~~


%% Overlay task kernels for striatum and ctx-predicted striatum

% Load data
data_fn = 'trial_activity_choiceworld';
exclude_data = true;
AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(mua_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_animal;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));

%%% Task>cortex>striatum regression results

% Get task>striatum parameters
regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
n_regressors = length(regressor_labels);
t_shifts = {[0,0.5]; ... % stim
    [-0.5,1]; ... % move
    [-0.1,0.5]; ... % go cue
    [-0.5,1]}; % outcome

sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
    round(x(2)*(sample_rate)),t_shifts,'uni',false);
t_shifts = cellfun(@(x) x/sample_rate,sample_shifts,'uni',false);

% Normalize task kernels (ridiculous loop makes it clearer to read)
mua_taskpred_k_all_norm = {};
mua_ctxpred_taskpred_k_all_norm = {};
for curr_animal = 1:length(mua_ctxpred_taskpred_k_all)
    for curr_exp = 1:length(mua_ctxpred_taskpred_k_all{curr_animal})
        for curr_regressor = 1:n_regressors
            for curr_depth = 1:n_depths
                mua_taskpred_k_all_norm{curr_animal}{curr_exp}{curr_regressor,curr_depth} = ...
                    mua_taskpred_k_all{curr_animal}{curr_exp}{curr_regressor,curr_depth}./ ...
                    ((1/sample_rate)*mua_norm{curr_animal}{curr_exp}(curr_depth));
                mua_ctxpred_taskpred_k_all_norm{curr_animal}{curr_exp}{curr_regressor,curr_depth} = ...
                    mua_ctxpred_taskpred_k_all{curr_animal}{curr_exp}{curr_regressor,curr_depth}./ ...
                    ((1/sample_rate)*mua_norm{curr_animal}{curr_exp}(curr_depth));
            end
        end
    end
end

% Average and concatenate task kernels within animals
task_str_k_animal = cell(n_regressors,n_depths);
task_ctx_str_k_animal = cell(n_regressors,n_depths);
for curr_animal = 1:length(mua_ctxpred_taskpred_k_all_norm)
    if isempty(mua_taskpred_k_all_norm{curr_animal})
        continue
    end
    curr_str_k = cat(3,mua_taskpred_k_all_norm{curr_animal}{:});
    curr_ctx_str_k = cat(3,mua_ctxpred_taskpred_k_all_norm{curr_animal}{:});
    for curr_depth = 1:n_depths
        for curr_regressor = 1:n_regressors
            curr_str_k_mean = nanmean(cat(3,curr_str_k{curr_regressor,curr_depth,:}),3);
            task_str_k_animal{curr_regressor,curr_depth} = cat(3, ...
                task_str_k_animal{curr_regressor,curr_depth},curr_str_k_mean);
            
            curr_ctx_str_k_mean = nanmean(cat(3,curr_ctx_str_k{curr_regressor,curr_depth,:}),3);
            task_ctx_str_k_animal{curr_regressor,curr_depth} = cat(3, ...
                task_ctx_str_k_animal{curr_regressor,curr_depth},curr_ctx_str_k_mean);
        end
    end
end

% Plot task>striatum kernels with task>cortex>striatum kernel shading
figure;
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors
        if curr_regressor == 1
            n_subregressors = 10;
            col = colormap_BlueWhiteRed(n_subregressors/2);
            col(6,:) = [];
        else
            n_subregressors = 2;
            col = lines(n_subregressors);
        end
        
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);  
        
        curr_task_str_kernels = task_str_k_animal{curr_regressor,curr_depth};
        curr_task_ctx_str_kernels = task_ctx_str_k_animal{curr_regressor,curr_depth};
        for curr_subregressor = 1:n_subregressors
            
            plot(t_shifts{curr_regressor}, ...
                nanmean(curr_task_str_kernels(curr_subregressor,:,:),3), ...
                'color',col(curr_subregressor,:),'linewidth',2);
            
            AP_errorfill(t_shifts{curr_regressor}, ...
                nanmean(curr_task_ctx_str_kernels(curr_subregressor,:,:),3), ...
                AP_sem(curr_task_ctx_str_kernels(curr_subregressor,:,:),3), ...
                col(curr_subregressor,:),0.5,false);
            
        end
        
        axis off;
        line([0,0],ylim,'color','k');
        
    end
end
linkaxes(p);

























