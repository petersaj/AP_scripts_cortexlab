% Generate revision figures for ctx-str paper
% (trials data is prepared in AP_ctx_str_trial_preprocessing)

% NOTE: these are just in order that I wrote them at the moment

%% [[LOAD DATASETS]]

% Load data

% (task)
data_fn = 'trial_activity_choiceworld'; % Primary dataset
% data_fn = 'trial_activity_choiceworld_4strdepth'; % Depth-aligned striatum
exclude_data = false;

% (passive)
% data_fn = 'trial_activity_AP_choiceWorldStimPassive_trained';
% data_fn = 'trial_activity_AP_choiceWorldStimPassive_naive';
% data_fn = 'trial_activity_stimKalatsky_naive';
% data_fn = 'trial_activity_stimKalatsky_trained';
% exclude_data = false;

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


%% CCF domain location / Allen corticostriatal projections

warning('THIS IS NOW SUBBED INTO MAIN FIGURES')

% Load the kernel template matches
n_aligned_depths = 3;
kernel_match_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
kernel_match_fn = ['ephys_kernel_align_' num2str(n_aligned_depths) '_depths.mat'];
load([kernel_match_path filesep kernel_match_fn]);

% Pad and concatenate kernel matches
max_n_kernel_match = max(cell2mat(cellfun(@(x) ...
    cellfun(@length,x),{ephys_kernel_align.kernel_match},'uni',false)));

kernel_match_all = cellfun(@(x) cell2mat(cellfun(@(x) ...
    padarray(x,max_n_kernel_match-length(x),NaN,'pre'),x,'uni',false)), ...
    {ephys_kernel_align.kernel_match},'uni',false);

kernel_match_cat = cell2mat(kernel_match_all);

% (use only depths in at least x% of recordings)
frac_depths = nanmean(~isnan(kernel_match_cat),2);
use_depths = frac_depths > 0.5;
kernel_match_cat_use = kernel_match_cat(use_depths,:);

% Get fraction of each domain for each depth
kernel_match_frac = nan(sum(use_depths),n_aligned_depths);
for curr_depth = 1:n_aligned_depths
    curr_kernel_match = +(kernel_match_cat_use == curr_depth);
    curr_kernel_match(isnan(kernel_match_cat_use)) = NaN;
    kernel_match_frac(:,curr_depth) = nanmean(curr_kernel_match,2);
end

% Get center-of-mass for each domain
kernel_match_com = nan(n_aligned_depths,1);
for curr_depth = 1:n_aligned_depths
    curr_kernel_match = +(kernel_match_cat_use == curr_depth);
    curr_kernel_match(isnan(kernel_match_cat_use)) = NaN;
        kernel_match_com(curr_depth) = ...
            nansum(nanmean(curr_kernel_match,2).*(1:sum(use_depths))') ...
            ./nansum(nanmean(curr_kernel_match,2));
end

% Plot domain locations
str_col = copper(n_aligned_depths);
figure; hold on;
colormap(str_col);
area(kernel_match_frac,'FaceColor','flat');
for curr_depth = 1:n_aligned_depths
    line(repmat(kernel_match_com(curr_depth),2,1),ylim, ...
        'color',min(str_col(curr_depth,:)+0.2,1),'linewidth',3);
end
axis tight;
xlabel('Estimated depth');
ylabel('Fraction domain match');

% Get relative domain COM along trajectory and plot in CCF
kernel_match_com_relative = kernel_match_com./sum(use_depths);

% Define the probe vector manually according to the targeted trajectory
probe_vector_ccf = [520,240,510;520,511,239];

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
kernel_depth_ccf = interp1([0,1],probe_str,kernel_match_com_relative);

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

% Plot brain area outlines at slice
callosum_id = find(strcmp(st.safe_name,'corpus callosum body'));
ventricle_id = find(strcmp(st.safe_name,'lateral ventricle'));

av_slice = permute(av(probe_vector_ccf(1),:,:),[2,3,1]);
slice_brain_outline = bwboundaries(av_slice > 1,'noholes');
slice_str_outline = bwboundaries(av_slice == str_id,'noholes');
slice_callosum_outline = bwboundaries(av_slice == callosum_id,'noholes');
slice_ventricle_outline = bwboundaries(av_slice == ventricle_id,'noholes');

figure; hold on; axis image on; box on; grid on; set(gca,'YDir','reverse');
plot(slice_brain_outline{1}(:,2),slice_brain_outline{1}(:,1),'k','linewidth',2);
cellfun(@(x) plot(x(:,2),x(:,1),'b','linewidth',2),slice_str_outline);
cellfun(@(x) fill(x(:,2),x(:,1),'k','linewidth',2),slice_callosum_outline);
cellfun(@(x) fill(x(:,2),x(:,1),'k','linewidth',2),slice_ventricle_outline);

line(probe_vector_ccf(:,3),probe_vector_ccf(:,2),'linewidth',2,'color','r');
scatter(kernel_depth_ccf(:,3),kernel_depth_ccf(:,2), ...
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
alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment';
ccf_tform_fn = [alignment_path filesep 'ccf_tform.mat'];
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
n_aligned_depths = 3;
kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template_' num2str(n_aligned_depths) '_depths.mat'];
load(kernel_template_fn);

figure; 
colormap(brewermap([],'PRGn'));
for curr_depth = 1:n_aligned_depths
    subplot(n_aligned_depths,1,curr_depth); hold on; axis image off;
    set(gca,'YDir','reverse');
    
    scatter(injection_coordinates_bilateral_wf{curr_depth}(:,1), ...
        injection_coordinates_bilateral_wf{curr_depth}(:,2), ...
        projection_strength_bilateral{curr_depth}*50 + 10, ...
        'k','filled');
    
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
end


    
    

%% Widefield task goodness-of-fit
% NOTE: striatum moved to new fig S8 in v4, what used to be here was bad

error('CHECK THIS! things went back and forth with baseline subtracting and common nans');

% Use raw fluorescence
fluor_exp = vertcat(fluor_all{:});

% Get R^2 for task in cortex ROI
taskpred_r2 = nan(max(split_idx),n_depths);
for curr_exp = 1:max(split_idx)
    
    curr_fluor_deconv = AP_deconv_wf(fluor_exp{curr_exp});
    curr_fluor_deconv_kernelroi = permute(reshape( ...
        AP_svd_roi(U_master(:,:,1:n_vs), ...
        reshape(permute(curr_fluor_deconv,[3,2,1]),n_vs,[]),[],[],kernel_roi.bw), ...
        size(kernel_roi.bw,3),[],size(curr_fluor_deconv,1)),[3,2,1]);
           
    curr_data_baselinesub = reshape(permute(curr_fluor_deconv_kernelroi,[2,1,3]),[],n_depths) - ...
        (nanmean(reshape(curr_fluor_deconv_kernelroi(:,t < 0,:),[],size(curr_fluor_deconv_kernelroi,3)),1));
    curr_taskpred_data = reshape(permute(fluor_kernelroi_taskpred(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    
    % Set common NaNs
    nan_samples = isnan(curr_data_baselinesub) | isnan(curr_taskpred_data);
    curr_data_baselinesub(nan_samples) = NaN;
    curr_taskpred_data(nan_samples) = NaN;

    taskpred_r2(curr_exp,:) = 1 - (nansum((curr_data_baselinesub-curr_taskpred_data).^2,1)./ ...
        nansum((curr_data_baselinesub-nanmean(curr_data_baselinesub,1)).^2,1));

end
figure; hold on;
errorbar(nanmean(taskpred_r2,1),AP_sem(taskpred_r2,1),'b','linewidth',2,'CapSize',0);
xlabel('Cortex ROI');
ylabel('Task explained variance');


% Cortex explained variance

% (spatial explained variance in pixels)
px_taskpred_r2 = nan(size(U_master,1),size(U_master,2),max(split_idx));
for curr_exp = 1:max(split_idx)  
    px_taskpred_r2(:,:,curr_exp) = AP_spatial_explained_var(U_master(:,:,1:n_vs), ...
        reshape(permute(fluor_allcat_deconv(split_idx == curr_exp,:,:),[2,1,3]),[],n_vs)', ...
        reshape(permute(fluor_taskpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_vs)',10);
    AP_print_progress_fraction(curr_exp,max(split_idx));
end

figure;imagesc(nanmedian(px_taskpred_r2,3));
axis image off; 
colormap(brewermap([],'Reds'));
caxis([0,1]); 
c = colorbar;
ylabel(c,'Task R^2');
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);






%% Probe location variation
% (plot widefield-estimated/histology-aligned probe location)

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
    alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment';
    ccf_tform_fn = [alignment_path filesep 'ccf_tform.mat'];
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

%% Striatum cortical kernels pre/post muscimol
error('I don''t think this should be included');

protocols = {'vanillaChoiceworldNoRepeats_pre_muscimol','vanillaChoiceworldNoRepeats_post_muscimol'};

for protocol = protocols 
    protocol = cell2mat(protocol);
    
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
    k_fn = [data_path filesep 'ctx_str_kernels_' protocol];
    load(k_fn);
    
    framerate = 35;
    upsample_factor = 1;
    sample_rate = framerate*upsample_factor;
    kernel_t = [-0.1,0.1];
    kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
    t = kernel_frames/sample_rate;
    
    % Concatenate explained variance
    expl_var_animal = cell2mat(cellfun(@(x) nanmean(horzcat(x{:}),2),ctx_str_expl_var','uni',false));
    figure('Name',protocol);
    errorbar(nanmean(expl_var_animal,2),AP_sem(expl_var_animal,2),'k','linewidth',2);
    xlabel('Striatal depth');
    ylabel('Fraction explained variance');
    
    % Concatenate and mean
    % (kernel is -:+ fluorescence lag, flip to be spike-oriented)
    k_px_timeflipped = cellfun(@(x) cellfun(@(x) x(:,:,end:-1:1,:),x,'uni',false),ctx_str_kernel,'uni',false);
    k_px_animal = cellfun(@(x) nanmean(cat(5,x{:}),5),k_px_timeflipped,'uni',false);
    k_px = nanmean(double(cat(5,k_px_animal{:})),5);
    
    % Get center-of-mass maps
    k_px_positive = k_px;
    k_px_positive(k_px_positive < 0) = 0;
    k_px_com = sum(k_px_positive.*permute(1:n_aligned_depths,[1,3,4,2]),4)./sum(k_px_positive,4);
    k_px_com_colored = nan(size(k_px_com,1),size(k_px_com,2),3,size(k_px_com,3));
    
    use_colormap = min(jet(255)-0.2,1);
    for curr_frame = 1:size(k_px_com,3)
        k_px_com_colored(:,:,:,curr_frame) = ...
            ind2rgb(round(mat2gray(k_px_com(:,:,curr_frame),...
            [1,n_aligned_depths])*size(use_colormap,1)),use_colormap);
    end
    
    % Plot center kernel frames independently at t = 0
    figure('Name',protocol);
    plot_frame = kernel_frames == 0;
    for curr_depth = 1:n_aligned_depths
       subplot(n_aligned_depths,1,curr_depth);
       imagesc(k_px(:,:,plot_frame,curr_depth));
       AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
       axis image off;
       colormap(brewermap([],'PRGn'));
       caxis([-0.01,0.01]);
    end
    
    % Plot center-of-mass color at select time points 
    plot_t = [-0.05:0.025:0.05];
    
    k_px_com_colored_t = ...
        permute(reshape(interp1(t,permute(reshape(k_px_com_colored,[],3,length(t)), ...
        [3,1,2]),plot_t),length(plot_t),size(k_px_com_colored,1), ...
        size(k_px_com_colored,2),3),[2,3,4,1]);
    
    k_px_max = squeeze(max(k_px,[],4));
    k_px_max_t = ...
        permute(reshape(interp1(t,reshape(k_px_max,[],length(t))', ...
        plot_t),length(plot_t),size(k_px_max,1), ...
        size(k_px_max,2)),[2,3,1]);
    
    weight_max = 0.005;
    figure('Name',protocol);
    for t_idx = 1:length(plot_t)
        subplot(1,length(plot_t),t_idx);
        p = image(k_px_com_colored_t(:,:,:,t_idx));
        set(p,'AlphaData', ...
            mat2gray(k_px_max_t(:,:,t_idx),[0,weight_max]));
        axis image off;
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        title([num2str(plot_t(t_idx)),' s']);
    end    
        
    % Plot movie of kernels
    AP_image_scroll(reshape(permute(k_px,[1,4,2,3]),size(k_px,1)*size(k_px,4),size(k_px,2),length(t)),t);
    colormap(brewermap([],'PRGn'));
    caxis([-max(caxis),max(caxis)]);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5],[],[size(k_px,1),size(k_px,2),size(k_px,4),1]);
    axis image off
    
    drawnow;
    
end



%% ~~~~~~~~ BELOW: integrated in to AP_ctx_str_figures_v5


%% Plot kernel ROIs

kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
load(kernel_roi_fn);

figure; hold on
n_depths = size(kernel_roi.bw,3);
str_col = max(hsv(n_depths)-0.2,0);
set(gca,'YDir','reverse');
AP_reference_outline('ccf_aligned','k');
for curr_depths = 1:n_depths
    curr_roi_boundary = cell2mat(bwboundaries(kernel_roi.bw(:,:,curr_depths)));
    patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),str_col(curr_depths,:));
end
axis image off;


%% Task kernel str/ctx correlation

mua_taskpred_catk = cellfun(@(x) cellfun(@(x) ...
    cell2mat(cellfun(@(x) reshape(x,[],size(x,3)),x,'uni',false)), ...
    x,'uni',false),mua_taskpred_k_all,'uni',false);

mua_ctxpred_taskpred_catk = cellfun(@(x) cellfun(@(x) ...
    cell2mat(cellfun(@(x) reshape(x,[],size(x,3)),x,'uni',false)), ...
    x,'uni',false),mua_ctxpred_taskpred_k_all,'uni',false);

ctx_str_taskk_animal = cellfun(@(x,y) [x,y], ...
    mua_taskpred_catk,mua_ctxpred_taskpred_catk,'uni',false);

ctx_str_taskk_corr = nan(4,n_depths,length(ctx_str_taskk_animal));
for curr_animal = 1:length(ctx_str_taskk_animal)
    
    curr_k = cellfun(@(x) reshape(x,[],n_depths), ...
        ctx_str_taskk_animal{curr_animal},'uni',false);
    
    curr_k_str = cat(3,curr_k{:,1});
    curr_k_ctx = cat(3,curr_k{:,2});
     
    % Correlate str/ctx kernels within domain
    ctx_str_taskk_corr(1,:,curr_animal) = ...
        nanmean(cell2mat(cellfun(@(x,y) diag(corr(x,y))', ...
        curr_k(:,1),curr_k(:,2),'uni',false)));
    
    % Correlate str/ctx kernels across domains
    ctx_str_taskk_corr(2,:,curr_animal) = ...
        nanmean(cell2mat(cellfun(@(x,y) ...
        nansum(tril(corr(x),-1)+triu(corr(x),1),1)./(n_depths-1)', ...
        curr_k(:,1),'uni',false)));
       
    % Correlate kernel within task/notask across days within domain
    ctx_str_taskk_corr(3,:,curr_animal) = arrayfun(@(depth) ...
        nanmean(AP_itril(corr(permute(curr_k_str(:,depth,:),[1,3,2])),-1)),1:n_depths);
    ctx_str_taskk_corr(4,:,curr_animal) = arrayfun(@(depth) ...
        nanmean(AP_itril(corr(permute(curr_k_ctx(:,depth,:),[1,3,2])),-1)),1:n_depths);  

end

% Get mean across domains
ctx_str_taskk_corr_strmean = squeeze(nanmean(ctx_str_taskk_corr,2));

% Plot mean and split by domains
figure; 

subplot(2,1,1);hold on; set(gca,'ColorOrder',copper(n_depths));
plot(ctx_str_taskk_corr_strmean,'color',[0.5,0.5,0.5]);
errorbar(nanmean(ctx_str_taskk_corr_strmean,2), ...
    AP_sem(ctx_str_taskk_corr_strmean,2),'k','linewidth',2);
set(gca,'XTick',1:4,'XTickLabelRotation',20,'XTickLabel', ...
    {'Str-ctx within day','Str within day across domains','Str across days','Ctx across days'})
ylabel('Task kernel correlation');
xlim([0.5,4.5]);

subplot(2,1,2);hold on; set(gca,'ColorOrder',copper(n_depths));
errorbar(nanmean(ctx_str_taskk_corr,3), ...
    AP_sem(ctx_str_taskk_corr,3),'linewidth',2)
set(gca,'XTick',1:4,'XTickLabelRotation',20,'XTickLabel', ...
    {'Str-ctx within day','Str within day across domains','Str across days','Ctx across days'})
ylabel('Task kernel correlation');
xlim([0.5,4.5]);
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false))

% (within task-passive v task-task domains statistics)
disp('Str/ctx vs str/str cross-domain:')
curr_p = signrank(squeeze(ctx_str_taskk_corr_strmean(1,:)), ...
    squeeze(ctx_str_taskk_corr_strmean(2,:)));
disp(['All str p = ' num2str(curr_p)]);
for curr_depth = 1:n_depths  
    curr_p = signrank(squeeze(ctx_str_taskk_corr(1,curr_depth,:)), ...
        squeeze(ctx_str_taskk_corr(2,curr_depth,:)));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end

% (within vs across statistics)
disp('Str/ctx-within vs str-across:')
curr_p = signrank(squeeze(ctx_str_taskk_corr_strmean(1,:)), ...
    squeeze(ctx_str_taskk_corr_strmean(3,:)));
disp(['All str ' num2str(curr_depth) ' p = ' num2str(curr_p)]);
for curr_depth = 1:n_depths
    curr_p = signrank(squeeze(ctx_str_taskk_corr(1,curr_depth,:)), ...
        squeeze(ctx_str_taskk_corr(3,curr_depth,:)));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end

% (cross task vs no task statistics)
disp('Str-across vs ctx-across');
curr_p = signrank(squeeze(ctx_str_taskk_corr_strmean(3,:)), ...
    squeeze(ctx_str_taskk_corr_strmean(4,:)));
disp(['All str ' num2str(curr_depth) ' p = ' num2str(curr_p)]);
for curr_depth = 1:n_depths
    curr_p = signrank(squeeze(ctx_str_taskk_corr(3,curr_depth,:)), ...
        squeeze(ctx_str_taskk_corr(4,curr_depth,:)));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end




%% Pre/post learning passive

% data_fns = { ...
%     'trial_activity_AP_choiceWorldStimPassive_naive', ...
%     'trial_activity_AP_choiceWorldStimPassive_trained'};

data_fns = { ...
    'trial_activity_AP_choiceWorldStimPassive_naive', ...
    {'trial_activity_AP_choiceWorldStimPassive_trained', ...
    'trial_activity_AP_lcrGratingPassive_ctxstrephys_str', ...
    'trial_activity_AP_lcrGratingPassive_pre_muscimol'}};

stimIDs = cell(2,1);
mua_training = cell(2,1);
mua_ctxpred_training = cell(2,1);
fluor_training = cell(2,1);
fluor_roi_training = cell(2,1);
fluor_kernelroi_training = cell(2,1);

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;

    % Split by experiment
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    
    % Get trials with movement during stim to exclude
    wheel_thresh = 0.025;
    quiescent_trials = ~any(abs(wheel_allcat(:,t >= -0.5 & t <= 0.5)) > wheel_thresh,2);
    quiescent_trials_exp = mat2cell(quiescent_trials,trials_recording,1);
    
    % Get stim and activity by experiment
    trial_stim_allcat_exp = mat2cell(trial_stim_allcat,trials_recording,1);  
    fluor_deconv_exp = mat2cell(fluor_allcat_deconv,trials_recording,length(t),n_vs);
    fluor_roi_deconv_exp = mat2cell(fluor_roi_deconv,trials_recording,length(t),numel(wf_roi));
    fluor_kernelroi_deconv_exp = mat2cell(fluor_kernelroi_deconv,trials_recording,length(t),n_depths);
    mua_allcat_exp = mat2cell(mua_allcat,trials_recording,length(t),n_depths);
    mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,trials_recording,length(t),n_depths);
    
    % Exclude trials with fluorescence spikes
    % (this is a dirty way to do this but don't have a better alt)
    fluor_spike_thresh = 100;
    fluor_spike_trial = cellfun(@(x) any(any(x > fluor_spike_thresh,2),3), ...
        fluor_kernelroi_deconv_exp,'uni',false);
    
    % Grab data
    stimIDs{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        trial_stim_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    mua_training{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        mua_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    mua_ctxpred_training{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        mua_ctxpred_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_training{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_roi_training{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_roi_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_kernelroi_training{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_kernelroi_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
end

% Get average stim time course (100%R stim)
mua_mean = cellfun(@(act,stim) cell2mat(cellfun(@(act,stim) ...
    nanmean(act(stim == 1,:,:),1),act,stim,'uni',false)), ...
    mua_training,stimIDs,'uni',false);

fluor_kernelroi_mean = cellfun(@(act,stim) cell2mat(cellfun(@(act,stim) ...
    nanmean(act(stim == 1,:,:),1),act,stim,'uni',false)), ...
    fluor_kernelroi_training,stimIDs,'uni',false);

% Get average activity in relevant stim period
stim_avg_t = [0,0.2];
stim_avg_t_idx = t >= stim_avg_t(1) & t <= stim_avg_t(2);

% (only for 100%R stim)
stim_act = cellfun(@(str,ctx,stim) cellfun(@(str,ctx,stim) ...
    [nanmean(str(stim == 1,stim_avg_t_idx,:),2), ...
    nanmean(ctx(stim == 1,stim_avg_t_idx,:),2)], ...
    str,ctx,stim,'uni',false), ...
    mua_training, fluor_kernelroi_training, stimIDs,'uni',false);


% Plot average cortex and striatum stim activity
figure;
p = gobjects(n_depths,3);
for curr_str = 1:n_depths
    
    p(curr_str,1) = subplot(n_depths,3,(curr_str-1)*3+1);
    AP_errorfill(t,nanmean(fluor_kernelroi_mean{1}(:,:,curr_str),1)', ...
        AP_sem(fluor_kernelroi_mean{1}(:,:,curr_str),1)','k');
    AP_errorfill(t,nanmean(fluor_kernelroi_mean{2}(:,:,curr_str),1)', ...
        AP_sem(fluor_kernelroi_mean{2}(:,:,curr_str),1),'r');
    ylabel('Cortex (\DeltaF/F)');
    line(repmat(stim_avg_t(1),2,1),ylim);
    line(repmat(stim_avg_t(2),2,1),ylim);
    
    p(curr_str,2) = subplot(n_depths,3,(curr_str-1)*3+2);
    AP_errorfill(t,nanmean(mua_mean{1}(:,:,curr_str),1)', ...
        AP_sem(mua_mean{1}(:,:,curr_str),1)','k');
    AP_errorfill(t,nanmean(mua_mean{2}(:,:,curr_str),1)', ...
        AP_sem(mua_mean{2}(:,:,curr_str),1)','r');
    xlabel('Time from stim (s)');
    ylabel('Striatum (baseline)');
    title(['Str ' num2str(curr_str)]);
    line(repmat(stim_avg_t(1),2,1),ylim);
    line(repmat(stim_avg_t(2),2,1),ylim);
    
    
    p(curr_str,3) = subplot(n_depths,3,(curr_str-1)*3+3);  hold on;
    
    curr_stim_act_untrained = cell2mat(cellfun(@(x) nanmean(x(:,:,curr_str,1),1),stim_act{1},'uni',false));
    curr_stim_act_trained = cell2mat(cellfun(@(x) nanmean(x(:,:,curr_str,1),1),stim_act{2},'uni',false));
    
    str_col = copper(n_depths);
    
    plot([nanmean(curr_stim_act_untrained(:,2),1),nanmean(curr_stim_act_trained(:,2),1)], ...
       [nanmean(curr_stim_act_untrained(:,1),1),nanmean(curr_stim_act_trained(:,1),1)], ...
       'color',str_col(curr_str,:),'linewidth',2);
    
    errorbar(nanmean(curr_stim_act_untrained(:,2),1), ...
        nanmean(curr_stim_act_untrained(:,1),1), ...
        AP_sem(curr_stim_act_untrained(:,1),1),AP_sem(curr_stim_act_untrained(:,1),1), ...
        AP_sem(curr_stim_act_untrained(:,2),1),AP_sem(curr_stim_act_untrained(:,2),1), ...
        'color',str_col(curr_str,:),'linewidth',2);
    errorbar(nanmean(curr_stim_act_trained(:,2),1), ...
        nanmean(curr_stim_act_trained(:,1),1), ...
        AP_sem(curr_stim_act_trained(:,1),1),AP_sem(curr_stim_act_trained(:,1),1), ...
        AP_sem(curr_stim_act_trained(:,2),1),AP_sem(curr_stim_act_trained(:,2),1), ...
        'color',str_col(curr_str,:),'linewidth',2);

   scatter(nanmean(curr_stim_act_untrained(:,2),1), ...
       nanmean(curr_stim_act_untrained(:,1),1),100, ...
       str_col(curr_str,:),'filled','MarkerEdgeColor',[0,0.7,0],'linewidth',2);
   scatter(nanmean(curr_stim_act_trained(:,2),1), ...
       nanmean(curr_stim_act_trained(:,1),1),100, ...
       str_col(curr_str,:),'filled','MarkerEdgeColor',[0.7,0,0],'linewidth',2);
   
   xlabel('Cortical ROI');
   ylabel('Striatum');
    
end
linkaxes(p(:,1:2),'x');
linkaxes(p(:,3),'xy');


% (Untrained/trained statistics)
stim_act_mean = cellfun(@(x) cell2mat(cellfun(@(x) nanmean(x,1), ...
    x,'uni',false)),stim_act,'uni',false);

disp('Untrained/trained:');
for curr_depth = 1:n_depths
    curr_p = ranksum(stim_act_mean{1}(:,1,curr_depth),stim_act_mean{2}(:,1,curr_depth));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
    
    curr_p = ranksum(stim_act_mean{1}(:,2,curr_depth),stim_act_mean{2}(:,2,curr_depth));
    disp(['Ctx ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end



%% Widefield correlation borders

wf_corr_borders_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_borders\wf_corr_borders.mat';
load(wf_corr_borders_fn);

% Spacing/downsampling (hardcoded - in preprocessing)
px_spacing = 20;
downsample_factor = 10;

% Get average correlation maps
wf_corr_map_recording = [wf_corr_borders(:).corr_map_downsamp];
wf_corr_map_cat = cat(3,wf_corr_map_recording{:});
wf_corr_map_mean = cell(size(wf_corr_map_cat,1),size(wf_corr_map_cat,2));
for curr_y = 1:size(wf_corr_map_cat,1)
    for curr_x = 1:size(wf_corr_map_cat,2)
        wf_corr_map_mean{curr_y,curr_x} = ...
            nanmean(cell2mat(wf_corr_map_cat(curr_y,curr_x,:)),3);
    end
end

% Get average borders
wf_corr_borders_cat = cell2mat(reshape([wf_corr_borders(:).corr_edges],1,1,[]));
wf_corr_borders_mean = nanmean(wf_corr_borders_cat,3);

figure;

[plot_maps_y,plot_maps_x] = ndgrid(4:5:size(wf_corr_map_mean,1)-4,4:5:size(wf_corr_map_mean,2));

% Plot sample correlation map locations
subplot(1,3,1,'YDir','reverse');
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
plot(plot_maps_y(:)*px_spacing,plot_maps_x(:)*px_spacing,'.r','MarkerSize',30);
axis image off;

% Plot sample correlation maps concat
subplot(1,3,2);
imagesc(cell2mat(wf_corr_map_mean(5:5:end,5:5:end)));
axis image off;
caxis([0,1]);
colormap(brewermap([],'Greys'));
c = colorbar;
ylabel(c,'Correlation');

% Plot borders
subplot(1,3,3);
imagesc(wf_corr_borders_mean);
axis image off;
caxis([0,prctile(wf_corr_borders_mean(:),70)])
colormap(brewermap([],'Greys'))
ccf_outline = AP_reference_outline('ccf_aligned',[1,0,0]);
cellfun(@(x) set(x,'linewidth',1),vertcat(ccf_outline{:}));



%% Probe trajectories histology vs widefield-estimated

% Load probe trajectories
histology_probe_ccf_all_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\histology_probe_ccf_all'];
load(histology_probe_ccf_all_fn);
n_animals = length(histology_probe_ccf_all);

% Load in the annotated Allen volume and names
allen_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
av = readNPY([allen_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean

figure;

for curr_animal = 1:n_animals
    
    % Plot average image with drawn probe overlay;
    subplot(n_animals,3,(curr_animal-1)*3+1);
    imagesc(histology_probe_ccf_all(curr_animal).im);
    line(histology_probe_ccf_all(curr_animal).probe_wf(:,1), ...
        histology_probe_ccf_all(curr_animal).probe_wf(:,2), ...
        'color','r','linewidth',1);  
    axis image off;
    caxis([0,prctile(histology_probe_ccf_all(curr_animal).im(:),95)]);
    colormap(gca,'gray');
    
    % Plot brain to overlay probes
    horizontal_axes = subplot(n_animals,3,(curr_animal-1)*3+2,'YDir','reverse'); hold on;
    axis image off;
    coronal_axes = subplot(n_animals,3,(curr_animal-1)*3+3,'YDir','reverse'); hold on;
    axis image off;
    
    % Plot projections
    coronal_outline = bwboundaries(permute((max(av,[],1)) > 1,[2,3,1]));
    horizontal_outline = bwboundaries(permute((max(av,[],2)) > 1,[3,1,2]));
    
    str_id = find(strcmp(st.safe_name,'Caudoputamen'));
    str_coronal_outline = bwboundaries(permute((max(av == str_id,[],1)) > 0,[2,3,1]));
    str_horizontal_outline = bwboundaries(permute((max(av == str_id,[],2)) > 0,[3,1,2]));

    plot(horizontal_axes,horizontal_outline{1}(:,2),horizontal_outline{1}(:,1),'k','linewidth',2);
    plot(horizontal_axes,str_horizontal_outline{1}(:,2),str_horizontal_outline{1}(:,1),'b','linewidth',2);
    plot(horizontal_axes,str_horizontal_outline{2}(:,2),str_horizontal_outline{2}(:,1),'b','linewidth',2);
    axis image off;
    
    plot(coronal_axes,coronal_outline{1}(:,2),coronal_outline{1}(:,1),'k','linewidth',2);
    plot(coronal_axes,str_coronal_outline{1}(:,2),str_coronal_outline{1}(:,1),'b','linewidth',2);
    plot(coronal_axes,str_coronal_outline{2}(:,2),str_coronal_outline{2}(:,1),'b','linewidth',2);
    axis image off;
    
    % Get estimated probe vector (widefield)
    probe_ccf_wf = histology_probe_ccf_all(curr_animal).probe_ccf_wf;
    r0 = mean(probe_ccf_wf,1);
    xyz = bsxfun(@minus,probe_ccf_wf,r0);
    [~,~,V] = svd(xyz,0);
    probe_direction = V(:,1);
    
    probe_vector_evaluate = [0,sign(probe_direction(2))*1000];
    probe_vector_draw_wf = round(bsxfun(@plus,bsxfun(@times,probe_vector_evaluate',probe_direction'),r0));
    
    line(horizontal_axes,probe_vector_draw_wf(:,1),probe_vector_draw_wf(:,3), ...
        'color','r','linewidth',2);
    line(coronal_axes,probe_vector_draw_wf(:,3),probe_vector_draw_wf(:,2), ...
        'color','r','linewidth',2);
    
    % Get estimated probe vector (histology)
    probe_ccf_histology = histology_probe_ccf_all(curr_animal).probe_ccf_histology;
    % (still working on orienting histology: force left hemisphere)
    bregma = allenCCFbregma;
    if sign(probe_ccf_histology(end,3) - probe_ccf_histology(1,3)) == 1
        probe_ccf_histology(:,3) = ...
            probe_ccf_histology(:,3) + ...
            2*(bregma(3) - probe_ccf_histology(:,3));
    end
    
    r0 = mean(probe_ccf_histology,1);
    xyz = bsxfun(@minus,probe_ccf_histology,r0);
    [~,~,V] = svd(xyz,0);
    probe_direction = V(:,1);
    
    probe_vector_evaluate = [-sign(probe_direction(2))*500,sign(probe_direction(2))*500];
    probe_vector_draw_histology = round(bsxfun(@plus,bsxfun(@times,probe_vector_evaluate',probe_direction'),r0));

    line(horizontal_axes,probe_vector_draw_histology(:,1),probe_vector_draw_histology(:,3), ...
        'color',[0,0.7,0],'linewidth',2);
    line(coronal_axes,probe_vector_draw_histology(:,3),probe_vector_draw_histology(:,2), ...
        'color',[0,0.7,0],'linewidth',2);
    
    drawnow;
    
end


%% Probe trajectories estimated from widefield image

% Load estimated probe trajectories
probe_ccf_all_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\probe_ccf_all'];
load(probe_ccf_all_fn);

% Load in the annotated Allen volume and names
allen_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
av = readNPY([allen_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean

% Plot brain to overlay probes
% (note the CCF is rotated to allow for dim 1 = x)
h = figure; 
ccf_axes = subplot(2,2,1); hold on
sagittal_axes = subplot(2,2,2,'YDir','reverse'); hold on;
axis image; grid on;
horizontal_axes = subplot(2,2,3,'YDir','reverse'); hold on;
axis image; grid on;
coronal_axes = subplot(2,2,4,'YDir','reverse'); hold on;
axis image; grid on;

% Plot 1 = 3D
% (Use wire mesh - can add other structures)
slice_spacing = 10;

str_id = find(strcmp(st.safe_name,'Caudoputamen'));
target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) == str_id,[2,1,3]);
structure_patch = isosurface(target_volume,0);
structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
target_structure_color = [0,0,1];
striatum_outline = patch(ccf_axes, ...
    'Vertices',structure_wire.vertices*slice_spacing, ...
    'Faces',structure_wire.faces, ...
    'FaceAlpha',0.5,'FaceColor',target_structure_color,'EdgeColor','none');

% (Use Nick's outline, but rotate to make dim 1 = x)
brainGridData = readNPY([fileparts(which('plotBrainGrid')) filesep 'brainGridData.npy']);
plotBrainGrid(brainGridData(:,[1,3,2]),ccf_axes);

axes(ccf_axes);
axis image vis3d off;
cameratoolbar(h,'SetCoordSys','y');
cameratoolbar(h,'SetMode','orbit');

% Plots 2-4 = projections
coronal_outline = bwboundaries(permute((max(av,[],1)) > 1,[2,3,1]));
horizontal_outline = bwboundaries(permute((max(av,[],2)) > 1,[3,1,2]));
sagittal_outline = bwboundaries(permute((max(av,[],3)) > 1,[2,1,3]));

str_id = find(strcmp(st.safe_name,'Caudoputamen'));
str_coronal_outline = bwboundaries(permute((max(av == str_id,[],1)) > 0,[2,3,1]));
str_horizontal_outline = bwboundaries(permute((max(av == str_id,[],2)) > 0,[3,1,2]));
str_sagittal_outline = bwboundaries(permute((max(av == str_id,[],3)) > 0,[2,1,3]));

plot(sagittal_axes,sagittal_outline{1}(:,2),sagittal_outline{1}(:,1),'k','linewidth',2);
plot(sagittal_axes,str_sagittal_outline{1}(:,2),str_sagittal_outline{1}(:,1),'b','linewidth',2);

plot(horizontal_axes,horizontal_outline{1}(:,2),horizontal_outline{1}(:,1),'k','linewidth',2);
plot(horizontal_axes,str_horizontal_outline{1}(:,2),str_horizontal_outline{1}(:,1),'b','linewidth',2);
plot(horizontal_axes,str_horizontal_outline{2}(:,2),str_horizontal_outline{2}(:,1),'b','linewidth',2);

plot(coronal_axes,coronal_outline{1}(:,2),coronal_outline{1}(:,1),'k','linewidth',2);
plot(coronal_axes,str_coronal_outline{1}(:,2),str_coronal_outline{1}(:,1),'b','linewidth',2);
plot(coronal_axes,str_coronal_outline{2}(:,2),str_coronal_outline{2}(:,1),'b','linewidth',2);

if length(probe_ccf_all) <= 17
    color_set = [brewermap(8,'Dark2');brewermap(8,'Set2')];
    probe_color = color_set(1:length(probe_ccf_all),:);
else
    error('More animals than colors: add another set');
end


for curr_animal = 1:length(probe_ccf_all)
    for curr_day = 1:size(probe_ccf_all{curr_animal},3)

        % Get estimated probe vector
        probe_ccf = probe_ccf_all{curr_animal}(:,:,curr_day);
        r0 = mean(probe_ccf,1);
        xyz = bsxfun(@minus,probe_ccf,r0);
        [~,~,V] = svd(xyz,0);
        probe_direction = V(:,1);
                
        probe_vector_evaluate = [0,sign(probe_direction(2))*1000];
        probe_vector_draw = round(bsxfun(@plus,bsxfun(@times,probe_vector_evaluate',probe_direction'),r0));
        
        % Plot current probe vector
        line(ccf_axes,probe_vector_draw(:,1),probe_vector_draw(:,2),probe_vector_draw(:,3), ...
            'color',probe_color(curr_animal,:),'linewidth',1);
        line(sagittal_axes,probe_vector_draw(:,1),probe_vector_draw(:,2), ...
            'color',probe_color(curr_animal,:),'linewidth',1);
        line(horizontal_axes,probe_vector_draw(:,1),probe_vector_draw(:,3), ...
            'color',probe_color(curr_animal,:),'linewidth',1);
        line(coronal_axes,probe_vector_draw(:,3),probe_vector_draw(:,2), ...
            'color',probe_color(curr_animal,:),'linewidth',1);
        
        drawnow;
        
    end
end

% Put a colormap on the side
cmap_ax = axes('Position',[0.95,0.2,0.2,0.6])
image(permute(probe_color,[1,3,2]));



%% Fluorescence deconvolution and fluor+ctx+str example

% Load and plot kernel
% (flip time to be fluor lag:lead spikes);
kernel_path = fileparts(which('AP_deconv_wf'))
kernel_fn = [kernel_path filesep 'gcamp6s_kernel.mat'];
load(kernel_fn);

gcamp6s_kernel_cat = fliplr(vertcat(gcamp6s_kernel.regression{:}));
gcamp6s_kernel_norm = gcamp6s_kernel_cat./max(abs(gcamp6s_kernel_cat),[],2);
gcamp6s_kernel_mean = nanmean(gcamp6s_kernel_norm)

figure; hold on;
plot(gcamp6s_kernel.regression_t,gcamp6s_kernel_norm','color',[0.5,0.5,0.5]);
plot(gcamp6s_kernel.regression_t,gcamp6s_kernel_mean,'k','linewidth',2);
line(xlim,[0,0],'color','k','linestyle','--');
line([0,0],ylim,'color','k','linestyle','--');


% Load Fluor+ctx-ephys data
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\ctx_deconv_traces');

% Concatenate spikes/deconvolved fluorescence 
% (full recording)
ctx_fluor_full_cat = cellfun(@cell2mat,{ctx_deconv_traces.fluor},'uni',false);
ctx_spikes_full_cat = cellfun(@cell2mat,{ctx_deconv_traces.spikes},'uni',false);
ctx_fluor_deconv_full_cat = cellfun(@cell2mat,{ctx_deconv_traces.fluor_deconv},'uni',false);

% (task)
ctx_spikes_task_cat = cellfun(@(x,protocol) ...
    cell2mat(x(strcmp(protocol,'vanillaChoiceworld'))), ...
    {ctx_deconv_traces.spikes},{ctx_deconv_traces.protocol},'uni',false);
ctx_fluor_deconv_task_cat = cellfun(@(x,protocol) ...
    cell2mat(x(strcmp(protocol,'vanillaChoiceworld'))), ...
    {ctx_deconv_traces.fluor_deconv},{ctx_deconv_traces.protocol},'uni',false);

% (passive)
ctx_spikes_notask_cat = cellfun(@(x,protocol) ...
    cell2mat(x(strcmp(protocol,'AP_sparseNoise'))), ...
    {ctx_deconv_traces.spikes},{ctx_deconv_traces.protocol},'uni',false);
ctx_fluor_deconv_notask_cat = cellfun(@(x,protocol) ...
    cell2mat(x(strcmp(protocol,'AP_sparseNoise'))), ...
    {ctx_deconv_traces.fluor_deconv},{ctx_deconv_traces.protocol},'uni',false);

% Get explained variance
ctx_deconv_full_r2 = cellfun(@(spikes,fluor_deconv) ...
    1-(nansum((spikes-fluor_deconv).^2)./nansum((spikes-nanmean(spikes)).^2)), ...
    ctx_spikes_full_cat,ctx_fluor_deconv_full_cat);

ctx_deconv_task_r2 = cellfun(@(spikes,fluor_deconv) ...
    1-(nansum((spikes-fluor_deconv).^2)./nansum((spikes-nanmean(spikes)).^2)), ...
    ctx_spikes_task_cat,ctx_fluor_deconv_task_cat);

ctx_deconv_notask_r2 = cellfun(@(spikes,fluor_deconv) ...
    1-(nansum((spikes-fluor_deconv).^2)./nansum((spikes-nanmean(spikes)).^2)), ...
    ctx_spikes_notask_cat,ctx_fluor_deconv_notask_cat);


%%% Plot example day

animal = 'AP060';
day = '2019-12-06';
experiment = 1;
plot_t = [100,300];

figure;
disp('Loading example data...');

% Load cortex ephys + imaging
load_parts.ephys = true;
load_parts.imaging = true;
site = 2; % (cortex always probe 2)
str_align = 'none'; % (cortex)
AP_load_experiment;

% Load cortex recording alignment
vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
load(vis_ctx_ephys_fn);

%%% DEPTH-ALIGN TEMPLATES, FIND CORTEX BOUNDARY
curr_animal_idx = strcmp(animal,{vis_ctx_ephys.animal});
curr_day_idx = strcmp(day,vis_ctx_ephys(curr_animal_idx).day);
curr_csd_depth = vis_ctx_ephys(curr_animal_idx).stim_csd_depth{curr_day_idx};
curr_csd_depth_aligned = vis_ctx_ephys(curr_animal_idx).stim_csd_depth_aligned{curr_day_idx};
template_depths_aligned = interp1(curr_csd_depth,curr_csd_depth_aligned,template_depths);
spike_depths_aligned = interp1(curr_csd_depth,curr_csd_depth_aligned,spike_depths);

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

% Plot cortex raster
subplot(6,1,1); hold on;
plot_t_idx = spike_binning_t_centers >= plot_t(1) & ...
    spike_binning_t_centers <= plot_t(2);
plot(spike_binning_t_centers(plot_t_idx), ...
    fluor_roi_interp(plot_t_idx),'linewidth',2,'color',[0,0.7,0]);

subplot(6,1,2,'YDir','reverse'); hold on;
plot_spikes = spike_times_timeline >= plot_t(1) & ...
    spike_times_timeline <= plot_t(2) & ...
    spike_depths_aligned <= ctx_depth(2);
plot(spike_times_timeline(plot_spikes),spike_depths(plot_spikes),'.k');
ylabel('Cortex depth (\mum)');
xlabel('Time (s)');

%%% LOAD STRIATUM EPHYS AND GET MUA BY DEPTH

% Load striatum ephys
load_parts.ephys = true;
load_parts.imaging = false;
site = 1; % (striatum is always on probe 1)
str_align = 'kernel';
AP_load_experiment;

% Plot striatum raster
subplot(6,1,3,'YDir','reverse'); hold on;
plot_spikes = spike_times_timeline >= plot_t(1) & ...
    spike_times_timeline <= plot_t(2) & ...
    spike_depths >= str_depth(1) & spike_depths <= str_depth(2);
plot(spike_times_timeline(plot_spikes),spike_depths(plot_spikes),'.k');
ylabel('Striatum depth (\mum)');
xlabel('Time (s)');

%%%%%% PLOT WHEEL/STIM

% (wheel velocity)
wheel_axes = subplot(6,1,4); hold on;
plot_wheel_idx = Timeline.rawDAQTimestamps >= plot_t(1) & ...
    Timeline.rawDAQTimestamps <= plot_t(2);
plot(wheel_axes,Timeline.rawDAQTimestamps(plot_wheel_idx), ...
    wheel_velocity(plot_wheel_idx),'k','linewidth',2);
ylabel('Wheel velocity');
axis off

% %     (stimuli)
%     % (task)
%     if contains(expDef,'vanilla')
%         stim_col = colormap_BlueWhiteRed(5);
%         [~,trial_contrast_idx] = ...
%             ismember(trial_conditions(:,1).*trial_conditions(:,2),unique(contrasts'.*sides),'rows');
%     elseif strcmp(expDef,'AP_lcrGratingPassive')
%         % (passive)
%         stim_col = [0,0,1;0.5,0.5,0.5;1,0,0];
%         [~,trial_contrast_idx] = ...
%             ismember(trial_conditions(:,1).*trial_conditions(:,2),[-90;0;90],'rows');
%     end
%     stim_lines = arrayfun(@(x) line(wheel_axes,repmat(stimOn_times(x),1,2),ylim(wheel_axes),'color', ...
%             stim_col(trial_contrast_idx(x),:),'linewidth',2), ...
%             find(stimOn_times >= plot_t(1) & stimOn_times <= plot_t(2)));
%
%     % (movement starts)
%     move_col = [0.6,0,0.6;0,0.6,0];
%     [~,trial_choice_idx] = ismember(trial_conditions(:,3),[-1;1],'rows');
%     move_lines = arrayfun(@(x) line(wheel_axes,repmat(wheel_move_time(x),1,2),ylim(wheel_axes),'color', ...
%         move_col(trial_choice_idx(x),:),'linewidth',2), ...
%         find(wheel_move_time >= plot_t(1) & wheel_move_time <= plot_t(2)));
%
%     % (go cues)
%     go_col = [0.8,0.8,0.2];
%     go_cue_times = signals_events.interactiveOnTimes(1:n_trials);
%     go_cue_lines = arrayfun(@(x) line(wheel_axes,repmat(go_cue_times(x),1,2),ylim(wheel_axes),'color', ...
%         go_col,'linewidth',2), ...
%         find(go_cue_times >= plot_t(1) & go_cue_times <= plot_t(2)));
%
%     % (outcomes)
%     outcome_col = [0,0,0.8;0.5,0.5,0.5];
%     reward_lines = arrayfun(@(x) line(wheel_axes,repmat(reward_t_timeline(x),1,2),ylim(wheel_axes),'color', ...
%         outcome_col(1,:),'linewidth',2), ...
%         find(reward_t_timeline >= plot_t(1) & reward_t_timeline <= plot_t(2)));
%     punish_times = signals_events.responseTimes(trial_outcome == -1);
%     punish_lines = arrayfun(@(x) line(wheel_axes,repmat(punish_times(x),1,2),ylim(wheel_axes),'color', ...
%         outcome_col(2,:),'linewidth',2), ...
%         find(punish_times >= plot_t(1) & punish_times <= plot_t(2)));


% Plot the cortex MUA and deconvolved fluorescence match
% Plot example
example_recording = 3;
example_protocol = 1;

% (raw fluorescence)
subplot(6,1,5); hold on;
plot(ctx_deconv_traces(example_recording).t{example_protocol}, ...
    mat2gray(ctx_deconv_traces(example_recording).fluor{example_protocol}), ...
    'color',[0,0.8,0],'linewidth',2);
ylabel('ROI Fluor(\DeltaF/F)');
xlabel('Time(s)')

% (cortex MUA and deconvolved fluorescence)
subplot(6,1,6); hold on;
plot(ctx_deconv_traces(example_recording).t{example_protocol}, ...
    ctx_deconv_traces(example_recording).spikes{example_protocol}, ...
    'color','k','linewidth',2);
plot(ctx_deconv_traces(example_recording).t{example_protocol}, ...
    ctx_deconv_traces(example_recording).fluor_deconv{example_protocol}, ...
    'color',[0,0.5,0],'linewidth',2);
linkaxes(get(gcf,'Children'),'x');
ylabel('Cortex spikes (std)');
xlabel('Time (s)');


curr_axes = flipud(get(gcf,'Children'));
% Link all time axes
linkaxes(curr_axes,'x');
% Link depth axes of raster plots (arbitrary depth, but want same scale)
linkaxes(curr_axes(2:3),'xy');

% (Display average and statistics)
disp(['Deconv explained var: ' num2str(nanmean(ctx_deconv_full_r2)) ...
    ' +/- ' num2str(AP_sem(ctx_deconv_full_r2,2))]);

p = signrank(ctx_deconv_task_r2,ctx_deconv_notask_r2);
disp(['Task vs no-task explained var: ' num2str(p)]);


%% Fluorescence/cortex ephys/striatum ephys correlation

%%% Load correlation data
use_protocol = 'vanillaChoiceworld';
% use_protocol = 'AP_sparseNoise';

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
data_fn = [data_path filesep 'ctx_fluor_mua_corr_' use_protocol];
load(data_fn);

%%% Get average aligned CSD

% Load CSD 
vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
load(vis_ctx_ephys_fn);

% Concatenate for plotting
animal_cat = cellfun(@(x,y) repmat({x},1,length(y)),{vis_ctx_ephys.animal},{vis_ctx_ephys.day},'uni',false);
animal_cat = horzcat(animal_cat{:});
day_cat = [vis_ctx_ephys(:).day];
n_recordings = length(day_cat);

stim_lfp_t = vis_ctx_ephys(1).stim_lfp_t{1};
stim_lfp_depth = vis_ctx_ephys(1).stim_lfp_depth{1};
stim_csd_depth = vis_ctx_ephys(1).stim_csd_depth{1};

stim_lfp_cat = cell2mat(permute([vis_ctx_ephys(:).stim_lfp],[1,3,2]));
stim_csd_cat = cell2mat(permute([vis_ctx_ephys(:).stim_csd],[1,3,2]));
csd_slice_cat = cell2mat([vis_ctx_ephys(:).csd_slice]);

% Get aligned depth for each recording, plot, save
figure;
curr_recording = 1;
stim_csd_aligned_scaled_cat = nan(0);
for curr_animal = 1:length(vis_ctx_ephys)
    for curr_day = 1:length(vis_ctx_ephys(curr_animal).day)
        
        stim_csd = vis_ctx_ephys(curr_animal).stim_csd{curr_day};
        csd_slice = vis_ctx_ephys(curr_animal).csd_slice{curr_day};
        stim_csd_depth = vis_ctx_ephys(curr_animal).stim_csd_depth{curr_day};
        stim_csd_depth_aligned = vis_ctx_ephys(curr_animal).stim_csd_depth_aligned{curr_day};
                
        % Plot aligned CSD
        depth_align_interp = [-500:20:2000];
        stim_csd_aligned = interp1(...
            stim_csd_depth_aligned,stim_csd_cat(:,:,curr_recording),depth_align_interp, ...
            'linear','extrap');
        
        animal = vis_ctx_ephys(curr_animal).animal;
        day = vis_ctx_ephys(curr_animal).day{curr_day};
        stim_lfp_t = vis_ctx_ephys(curr_animal).stim_lfp_t{curr_day};
        stim_lfp_depth = vis_ctx_ephys(curr_animal).stim_lfp_depth{curr_day};
        stim_lfp = vis_ctx_ephys(curr_animal).stim_lfp{curr_day};
        
        subplot(3,n_recordings,curr_recording);
        imagesc(stim_lfp_t,stim_lfp_depth,stim_lfp);
        caxis([-max(abs(caxis))/2,max(abs(caxis))/2]);
        colormap(brewermap([],'PRGn'));
        title({animal,day,'LFP'});
        
        subplot(3,n_recordings,size(stim_csd_cat,3)+curr_recording);
        imagesc(stim_lfp_t,stim_csd_depth,stim_csd);
        caxis([-max(abs(caxis))/2,max(abs(caxis))/2]);
        colormap(brewermap([],'PRGn'));
        title('CSD');
        
        subplot(3,n_recordings,size(stim_csd_cat,3)*2+curr_recording);
        imagesc(stim_lfp_t,stim_csd_depth_aligned,stim_csd_aligned);
        caxis([-max(abs(caxis))/2,max(abs(caxis))/2]);
        colormap(brewermap([],'PRGn'));
        title('Aligned CSD');
                
        % Keep aligned CSD for averaging
        stim_csd_aligned_scaled_cat(:,:,curr_recording) = stim_csd_aligned./min(csd_slice);       
        
        % Iterate recording index (only used for plotting)
        curr_recording = curr_recording + 1;

    end
end

stim_csd_aligned_mean = nanmean(stim_csd_aligned_scaled_cat,3);

%%% Plot correlations

mua_depth = data(1).cortex_mua_depth{1}; % they're all the same, use 1st
cortex_fluor_corr_cat = cell2mat(horzcat(data.cortex_fluor_corr));
cortex_striatum_corr_cat = cell2mat(permute(horzcat(data.cortex_striatum_corr),[1,3,2]));

figure;

subplot(1,4,1,'YDir','reverse');
imagesc(stim_lfp_t,depth_align_interp,stim_csd_aligned_mean);
caxis([-1,1]);
line([0,0],ylim,'color','k');
colormap(brewermap([],'PRGn'));
xlabel('Time from stim');
ylabel('Aligned depth');
title('Average aligned CSD')

subplot(1,4,2,'YDir','reverse'); hold on;
plot(cortex_fluor_corr_cat',mua_depth,'color',[0.2,0.8,0.2]);
errorbar(nanmean(cortex_fluor_corr_cat,2), ...
    mua_depth,AP_sem(cortex_fluor_corr_cat,2), ...
    'horizontal','color',[0,0.6,0],'linewidth',2)
xlabel('Correlation');
ylabel('Cortical MUA aligned depth');
title('Cortical fluorescence');

subplot(1,4,3,'YDir','reverse'); hold on;
set(gca,'ColorOrder',copper(4));
plot_str = 1;
plot(permute(cortex_striatum_corr_cat(:,plot_str,:),[1,3,2]),mua_depth,'color',[0.5,0.5,0.5]);
errorbar(nanmean(cortex_striatum_corr_cat(:,plot_str,:),3), ...
    mua_depth,AP_sem(cortex_striatum_corr_cat(:,plot_str,:),3),'k','horizontal','linewidth',2)
xlabel('Correlation');
ylabel('Cortical MUA aligned depth');
title(['Str ' num2str(plot_str) ' multiunit']);

subplot(2,4,4); hold on;
for i = 1:size(cortex_fluor_corr_cat,2)
    plot(cortex_fluor_corr_cat(:,i), ...
        cortex_striatum_corr_cat(:,plot_str,i), ...
        'color',[0.5,0.5,0.5]);
end
plot(nanmean(cortex_fluor_corr_cat,2), ...
    nanmean(cortex_striatum_corr_cat(:,plot_str,:),3),'k','linewidth',2);
xlabel('Fluorescence - cortical MUA correlation');
ylabel('Cortical MUA - striatal MUA correlation')

subplot(2,4,8); hold on;
plot(permute(max(cortex_striatum_corr_cat,[],1),[2,3,1]),'color',[0.5,0.5,0.5]);
errorbar(squeeze(nanmean(max(cortex_striatum_corr_cat,[],1),3)), ...
    squeeze(AP_sem(max(cortex_striatum_corr_cat,[],1),3)),'k','linewidth',2);
xlim([0.5,3.5]);
xlabel('Striatal domain');
ylabel('Cortical MUA max corr');


% (Fluor-Ctx MUA vs Str-Ctx MUA correlation by depth statistics)
use_str = 1;
use_ctx_str_corr = squeeze(cortex_striatum_corr_cat(:,use_str,:));
fluor_ctx_str_corr = diag(corr(use_ctx_str_corr,cortex_fluor_corr_cat,'rows','complete'));

n_shuff = 10000;
fluor_ctx_str_corr_shuff = nan(size(c,1),n_shuff);
for curr_shuff = 1:n_shuff
    use_ctx_str_corr_circshift = use_ctx_str_corr;
    for i = 1:size(use_ctx_str_corr_circshift,2)
        curr_depth_idx = ~isnan(use_ctx_str_corr_circshift(:,i));
        use_ctx_str_corr_circshift(curr_depth_idx,i) = circshift(use_ctx_str_corr_circshift(curr_depth_idx,i), ...
            randi(sum(curr_depth_idx)));
    end   
    fluor_ctx_str_corr_shuff(:,curr_shuff) = ...
        diag(corr(use_ctx_str_corr_circshift, ...
        cortex_fluor_corr_cat,'rows','complete'));
end
corr_rank = tiedrank([nanmean(fluor_ctx_str_corr,1),nanmean(fluor_ctx_str_corr_shuff)]);
corr_p = 1 - (corr_rank(1)/(n_shuff+1));

disp('Fluor-ctx depth vs Str-ctx depth correlation (circshift stat):')
disp(['p = ' num2str(corr_p) ', r = ' num2str(nanmean(fluor_ctx_str_corr)) ...
    ' +/- SEM ' num2str(AP_sem(fluor_ctx_str_corr,1))])



%% Regress striatum from cortical subregions and striatal domains
% NOTE: this regression is done on trial data rather than the long time
% courses which is what the normal analysis uses. For sanity check, the
% explained using the full data (full) and the trials dataset (trials) are
% compared here but not meant to be included.
%
% Also this takes a long time to run

% Choose depths to run
plot_depth = 1:n_depths;

% Regress kernel ROI activity to striatum domain activity (per recording)
regression_params.use_svs = 1:100;
regression_params.kernel_t = [-0.1,0.1];
regression_params.zs = [false,false];
regression_params.cvfold = 2;
regression_params.use_constant = true;
lambda = 50; % lambda for fluorescence to striatum
domain_lambda = 0; % lambda for other domains to striatum
kernel_frames = floor(regression_params.kernel_t(1)*sample_rate): ...
    ceil(regression_params.kernel_t(2)*sample_rate);

% Set time to use
% (e.g. this type of code can be used to regress only ITI time)
use_t = true(size(t));

% Set regions to zero
% (zero all EXCEPT quadrants)
ctx_zero = true(size(U_master,1),size(U_master,2),6);
ctx_zero(1:260,1:220,1) = false;
ctx_zero(1:260,220:end,2) = false;
ctx_zero(260:end,1:220,3) = false;
ctx_zero(260:end,220:end,4) = false;
ctx_zero(:,220:end,5) = false;
ctx_zero(:,1:220,6) = false;

% (use raw data for trial regression)
mua_exp = vertcat(mua_all{:});
fluor_exp = vertcat(fluor_all{:});
fluor_deconv_exp = cellfun(@AP_deconv_wf,fluor_exp,'uni',false);

mua_ctxtrialpred_exp = cellfun(@(x) nan(size(x)),mua_exp,'uni',false);
mua_ctxtrialpred_k = nan(length(regression_params.use_svs),length(kernel_frames),n_depths,length(use_split));
mua_ctxtrialpred_regionzero_exp = cellfun(@(x) nan([size(x),size(ctx_zero,3)]),mua_exp,'uni',false);
mua_ctxtrialpred_regionzero_k = nan(length(regression_params.use_svs),length(kernel_frames),n_depths,length(use_split));

mua_domaintrialpred_exp = cellfun(@(x) nan(size(x)),mua_exp,'uni',false);
mua_domaintrialpred_k = nan(n_depths-1,length(kernel_frames),n_depths,length(use_split));

for curr_exp = 1:length(trials_recording)
    for curr_depth = plot_depth
        
        curr_mua = reshape(mua_exp{curr_exp}(:,use_t,curr_depth)',1,[]);
        curr_mua_std = curr_mua./nanstd(curr_mua);
        
        curr_fluor = reshape(permute( ...
            fluor_deconv_exp{curr_exp}(:,use_t,:),[2,1,3]),[],n_vs)';
        curr_fluor_full = reshape(permute( ...
            fluor_deconv_exp{curr_exp}(:,:,:),[2,1,3]),[],n_vs)';
        
        curr_fluor_regionzero = nan(size(curr_fluor,1),size(curr_fluor,2),size(ctx_zero,3));
        for curr_zero = 1:size(ctx_zero,3)
            U_master_regionzero = U_master.*~ctx_zero(:,:,curr_zero);
            fVdf_regionzero_altU = ChangeU(U_master(:,:,1:n_vs),curr_fluor,U_master_regionzero(:,:,1:n_vs));
            % (those U's aren't orthonormal, recast back to original Udf_aligned)
            fVdf_regionzero = ChangeU(U_master_regionzero(:,:,1:n_vs),fVdf_regionzero_altU,U_master(:,:,1:n_vs));
            curr_fluor_regionzero(:,:,curr_zero) = fVdf_regionzero;
        end
        
        % Skip if no data
        if all(isnan(curr_mua))
            continue
        end
        
        % Set discontinuities in trial data
        trial_discontinuities = false(size(mua_exp{curr_exp}(:,use_t,curr_depth)));
        trial_discontinuities(:,1) = true;
        trial_discontinuities = reshape(trial_discontinuities',[],1)';
        
        % Full cortex regression
        [ctx_str_k,curr_mua_fluorpred_std,explained_var_trial] = ...
            AP_regresskernel(curr_fluor(regression_params.use_svs,:), ...
            curr_mua_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant,trial_discontinuities);
                
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        ctxpred_spikes = (curr_mua_fluorpred_std - squeeze(ctx_str_k{end})).* ...
            nanstd(curr_mua,[],2) + ...
            nanstd(curr_mua,[],2).*squeeze(ctx_str_k{end});
        
        mua_ctxtrialpred_k(:,:,curr_depth,curr_exp) = ctx_str_k{1};
        mua_ctxtrialpred_exp{curr_exp}(:,use_t,curr_depth) = ...
            reshape(ctxpred_spikes,sum(use_t),[])';
        %         % (apply kernel to full time)
        %         mua_ctxtrialpred_exp{curr_exp}(:,:,curr_depth) = ...
        %             sum(cell2mat(arrayfun(@(x) ...
        %             convn(fluor_allcat_deconv_exp{curr_exp}(:,:,regression_params.use_svs(x)), ...
        %             k_fluor(x,:)','same'),permute(1:length(regression_params.use_svs),[1,3,2]),'uni',false)),3);
        
        % Region-zeroed cortex regression
        for curr_zero = 1:size(ctx_zero,3)
            [ctx_str_k_regionzero,curr_mua_fluorregionzeropred_std,explained_var_trial_regionzero] = ...
                AP_regresskernel(curr_fluor_regionzero(regression_params.use_svs,:,curr_zero), ...
                curr_mua_std,kernel_frames,lambda, ...
                regression_params.zs,regression_params.cvfold, ...
                true,regression_params.use_constant,trial_discontinuities);
            
            % Re-scale the prediction (subtract offset, multiply, add scaled offset)
            ctxpred_spikes_regionzero = (curr_mua_fluorregionzeropred_std - squeeze(ctx_str_k_regionzero{end})).* ...
                nanstd(curr_mua,[],2) + ...
                nanstd(curr_mua,[],2).*squeeze(ctx_str_k_regionzero{end});
                        
            mua_ctxtrialpred_regionzero_k(:,:,curr_depth,curr_exp,curr_zero) = ctx_str_k_regionzero{1};
            mua_ctxtrialpred_regionzero_exp{curr_exp}(:,use_t,curr_depth,curr_zero) = ...
                reshape(ctxpred_spikes_regionzero,sum(use_t),[])';
            %         % (apply kernel to full time)
            %         mua_ctxtrialpred_regionzero_exp{curr_exp}(:,:,curr_depth) = ...
            %             sum(cell2mat(arrayfun(@(x) ...
            %             convn(fluor_allcat_deconv_exp{curr_exp}(:,:,regression_params.use_svs(x)), ...
            %             k_fluorregionzero(x,:)','same'),permute(1:length(regression_params.use_svs),[1,3,2]),'uni',false)),3);
        end
        
        % Regress current domains from other domains
        curr_mua_domains = reshape(permute(mua_exp{curr_exp}(:,use_t,:),[2,1,3]),[],n_depths)';
        curr_mua_domains_std = curr_mua_domains./nanstd(curr_mua_domains,[],2);
        
        [domain_str_k,curr_mua_domainpred_std,explained_var_trial] = ...
            AP_regresskernel(curr_mua_domains_std(setdiff(1:n_depths,curr_depth),:), ...
            curr_mua_std,kernel_frames,domain_lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant,trial_discontinuities);
        
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        curr_mua_domainpred = (curr_mua_domainpred_std - squeeze(domain_str_k{end})).* ...
            nanstd(curr_mua,[],2) + ...
            nanstd(curr_mua,[],2).*squeeze(domain_str_k{end});
        
        mua_domaintrialpred_k(:,:,curr_depth,curr_exp) = domain_str_k{1};
        mua_domaintrialpred_exp{curr_exp}(:,use_t,curr_depth) = ...
            reshape(curr_mua_domainpred,sum(use_t),[])';
      
        
    end
    AP_print_progress_fraction(curr_exp,length(trials_recording));
end

% Plot average cortex->striatum kernel
ctx_str_k_mean = nanmean(cell2mat(permute(vertcat(ctx_str_k_all{:}),[2,3,4,5,1])),5);
ctx_str_k_mean_px = cell2mat(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:100), ...
    ctx_str_k_mean(:,:,x)),permute(1:n_depths,[1,3,4,2]),'uni',false));

mua_ctxtrialpred_k_mean = nanmean(mua_ctxtrialpred_k,4);
mua_ctxtrialpred_k_mean_px = cell2mat(arrayfun(@(x) ...
    svdFrameReconstruct(U_master(:,:,regression_params.use_svs), ...
    mua_ctxtrialpred_k_mean(:,:,x)),permute(1:n_depths,[1,3,4,2]),'uni',false));

mua_ctxtrialpred_regionzero_k_mean = nanmean(mua_ctxtrialpred_regionzero_k,4);
mua_ctxtrialpred_regionzero_k_mean_px = cell2mat(arrayfun(@(x) ...
    svdFrameReconstruct(U_master(:,:,regression_params.use_svs), ...
    reshape(mua_ctxtrialpred_regionzero_k_mean(:,:,x,:),length(regression_params.use_svs),[])), ...
    permute(1:n_depths,[1,3,4,2]),'uni',false));

AP_image_scroll(cat(3,mua_ctxtrialpred_k_mean_px,mua_ctxtrialpred_regionzero_k_mean_px));
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(gca,brewermap([],'PRGn'));
axis image;

% Plot average domain->striautm kernel
figure;
mua_domaintrialpred_k_mean = nanmean(mua_domaintrialpred_k,4);
str_col = max(hsv(n_depths)-0.2,0);
for curr_depth = 1:n_depths
    subplot(n_depths,1,curr_depth); hold on
    set(gca,'ColorOrder',str_col(setdiff(1:n_depths,curr_depth),:));
    plot(kernel_frames,mua_domaintrialpred_k_mean(:,:,curr_depth)','linewidth',2);
    ylabel('Weight');
    title(['Str ' num2str(curr_depth)]);
end


% Get R^2 for task, cortex full, and cortex ROI predictions
mua_ctxpred_exp = vertcat(mua_ctxpred_all{:});

ctxpred_r2 = nan(max(split_idx),n_depths);
ctxtrialpred_r2 = nan(max(split_idx),n_depths);
ctxtrialpred_regionzero_r2 = nan(max(split_idx),n_depths,size(ctx_zero,3));
domaintrialpred_r2 = nan(max(split_idx),n_depths);

for curr_exp = 1:max(split_idx)
    
    curr_data = reshape(permute(mua_exp{curr_exp},[2,1,3]),[],n_depths);
    curr_ctxpred_data = reshape(permute(mua_ctxpred_exp{curr_exp},[2,1,3]),[],n_depths);
    curr_ctxtrialpred_data = reshape(permute(mua_ctxtrialpred_exp{curr_exp},[2,1,3]),[],n_depths);
    curr_ctxtrialpred_regionzero_data = ...
        reshape(permute(mua_ctxtrialpred_regionzero_exp{curr_exp},[2,1,3,4]),[],n_depths,size(ctx_zero,3));
    curr_domaintrialpred_data = reshape(permute(mua_domaintrialpred_exp{curr_exp},[2,1,3]),[],n_depths);
    
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_ctxpred_data) ...
        | isnan(curr_ctxtrialpred_data) | any(isnan(curr_ctxtrialpred_regionzero_data),3) ...
        | isnan(curr_domaintrialpred_data);
    curr_data(nan_samples) = NaN;
    curr_ctxpred_data(nan_samples) = NaN;
    curr_ctxtrialpred_data(nan_samples) = NaN;
    curr_ctxtrialpred_regionzero_data(repmat(nan_samples,1,1,size(ctx_zero,3))) = NaN;
    curr_domaintrialpred_data(nan_samples) = NaN;
    
    ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    ctxtrialpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxtrialpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    ctxtrialpred_regionzero_r2(curr_exp,:,:) = 1 - (nansum((curr_data-curr_ctxtrialpred_regionzero_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    domaintrialpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_domaintrialpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    
end

%%% Cortex subregions explained variance
figure;

% Plot full vs trial (sanity check: they should be ~the same)
subplot(2,3,1); hold on;
errorbar(nanmean(ctxpred_r2,1),AP_sem(ctxpred_r2,1),'.-','linewidth',2,'CapSize',0,'MarkerSize',30);
errorbar(nanmean(ctxtrialpred_r2,1),AP_sem(ctxtrialpred_r2,1),'.-','linewidth',2,'CapSize',0,'MarkerSize',30);
legend({'Ctx full','Ctx trial'});
xlabel('Striatum depth');
ylabel('Explained variance');

% Plot explained variance by subregion
subplot(2,3,2); hold on;
str_col = max(hsv(n_depths)-0.2,0);
set(gca,'ColorOrder',str_col);
errorbar(permute(nanmean(ctxtrialpred_regionzero_r2,1),[3,2,1]), ...
    permute(AP_sem(ctxtrialpred_regionzero_r2,1),[3,2,1]),'.-','linewidth',2,'CapSize',0,'MarkerSize',30);
xlabel('Cortex subregion');
ylabel('Explained variance');
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false));

% Plot explained variance of subregion and domain relative to full
alternate_r2 = cat(3,ctxtrialpred_regionzero_r2,domaintrialpred_r2);
alternate_r2_relative = ...
    (alternate_r2 - ctxtrialpred_r2)./ctxtrialpred_r2;

subplot(2,3,3); hold on;
str_col = max(hsv(n_depths)-0.2,0);
set(gca,'ColorOrder',str_col);
errorbar(permute(nanmean(alternate_r2_relative,1),[3,2,1]), ...
    permute(AP_sem(alternate_r2_relative,1),[3,2,1]),'linewidth',2,'CapSize',0,'Marker','none');
set(gca,'XTick',1:size(alternate_r2,3),'XTickLabel', ...
    [cellfun(@(x) ['Ctx ' num2str(x)],num2cell(1:size(ctx_zero,3)),'uni',false), ...
    'Domains']);
line(xlim,[0,0],'color','k','linestyle','--');
xlabel('Alternate regressor');
ylabel('\DeltaR^2/R^2_{full}')
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false));

% Plot subregions used
for i = 1:size(ctx_zero,3)
   subplot(2,size(ctx_zero,3),size(ctx_zero,3)+i);
   imagesc(~ctx_zero(:,:,i));
   AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
   axis image off;
   colormap(gray);
end


% (Prediction with cortex subsets statistics)
disp('R^2 with cortical subset (1-way anova):');
for curr_depth = 1:n_depths   
    curr_r2 = permute(ctxtrialpred_regionzero_r2_relative(:,curr_depth,:),[3,1,2]);
    [condition_grp,exp_grp] = ndgrid(1:size(curr_r2,1),1:size(curr_r2,2));
    curr_p = anovan(curr_r2(:),condition_grp(:),'display','off');  
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p')]);
end

% (Prediction with domains statistics)
disp('Domain vs. cortex explained variance (signrank):');
for curr_depth = 1:n_depths   
    curr_p = signrank(domaintrialpred_r2_relative(:,curr_depth)); 
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p')]);
end




%%  --- Cortical muscimol


%% ^^^ Cortical spike rate pre/post muscimol
% (use spike rate over all experiments pre/post muscimol)

animal_days = { ...
    'AP052','2019-09-20';
    'AP058','2019-12-06'};

figure;

spike_rate_change_cond = cell(size(animaldays));
for curr_animalday = 1:length(animal_days)
    
    animal = animal_days{curr_animalday,1};
    day = animal_days{curr_animalday,2};
    
    % Load data (first experiment - but spikes throughout used)
    experiment = 1;
    AP_load_experiment
    
    % Estimate boundaries of cortex (the dumb way: first template/gap)
    sorted_template_depths = sort(template_depths);
    ctx_start = sorted_template_depths(1) - 1;
    [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
    ctx_end = sorted_template_depths(max_gap_idx)+1;
    ctx_depth = [ctx_start,ctx_end];
    ctx_templates = template_depths <= ctx_depth(2);
    
    % Set experiments in conditions (1-2 = pre-muscimol, 3-4 = post-muscimol)
    cond_expts = {[1,2],[3,4]};
    
    spike_rate_cond = nan(max(spike_templates),2);
    for curr_cond = 1:2
        
        exp_starts = sync(2).timestamps(sync(2).values == 1);
        exp_stops = sync(2).timestamps(sync(2).values == 0);
        
        curr_exp_start = exp_starts(cond_expts{curr_cond}(1));
        curr_exp_stop = exp_stops(cond_expts{curr_cond}(end));
        
        curr_use_spikes = spike_times >= curr_exp_start & ...
            spike_times <= curr_exp_stop;
        
        spike_rate_cond(:,curr_cond) = ...
            accumarray(spike_templates(curr_use_spikes),1,[max(spike_templates),1])./ ...
            (curr_exp_stop - curr_exp_start);
        
    end
    
    spike_rate_change = (spike_rate_cond(:,2) - spike_rate_cond(:,1))./(spike_rate_cond(:,1)+spike_rate_cond(:,2));
    
    subplot(length(animal_days),2,(curr_animalday-1)*2+1,'YDir','reverse'); hold on;
    plot(reshape([spike_rate_cond,nan(size(template_depths))]',[],1), ...
        reshape(repmat(template_depths,1,3)',[],1),'color',[0.5,0.5,0.5]);
    p1 = plot(spike_rate_cond(:,1),template_depths,'.k','MarkerSize',10);
    p2 = plot(spike_rate_cond(:,2),template_depths,'.r','MarkerSize',10);
    xlabel('Spikes/s')
    ylabel('Depth (\mum)');
    legend([p1,p2],{'Pre-muscimol','Post-muscimol'});
    axis tight;
    xlim([-1,prctile(spike_rate_cond(:,1),95)])
   
    subplot(length(animal_days),2,(curr_animalday-1)*2+2); hold on;
    plot(spike_rate_change,template_depths,'.k','MarkerSize',10);   
    axis tight;
    xlim([-1.1,1.1]);
    line([0,0],ylim);
    set(gca,'YDir','reverse');
    xlabel('(Post-pre)/(pre+post)');
    ylabel('Depth (\mum)');
    title({animal,day,'Change'});
    
    
end



%% ^^^ VFS pre/post musicmol

muscimol_wf_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\muscimol_wf.mat';
load(muscimol_wf_fn);

std_cat = vertcat(muscimol_wf.std);
vfs_cat = vertcat(muscimol_wf.vfs);

std_change = (cat(3,std_cat{:,2})-cat(3,std_cat{:,1}))./(cat(3,std_cat{:,2})+cat(3,std_cat{:,1}));
vfs_softnorm = 0.2;
vfs_change = (abs(cat(3,vfs_cat{:,2}))-abs(cat(3,vfs_cat{:,1})))./(vfs_softnorm+abs(cat(3,vfs_cat{:,2}))+abs(cat(3,vfs_cat{:,1})));

figure; 

% Plot std
subplot(2,3,1);
imagesc(nanmean(cat(3,std_cat{:,1}),3));
axis image off;
caxis([0,0.03]);
colormap(gca,brewermap([],'*Greys'));
AP_reference_outline('ccf_aligned','r');
title('Std (pre-muscimol');
colorbar

subplot(2,3,2);
imagesc(nanmean(cat(3,std_cat{:,2}),3));
axis image off;
caxis([0,0.03]);
colormap(gca,brewermap([],'*Greys'));
AP_reference_outline('ccf_aligned','r');
title('Std (post-muscimol');
colorbar

subplot(2,3,3);
imagesc(nanmean(std_change,3));
axis image off;
caxis([-0.5,0.5]);
colormap(gca,brewermap([],'PRGn'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('(post-pre)/(post+pre)');
colorbar

% Plot VFS
subplot(2,3,4);
imagesc(nanmean(cat(3,vfs_cat{:,1}),3));
axis image off;
caxis([-1,1]);
colormap(gca,brewermap([],'PRGn'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('VFS pre-muscimol');
colorbar

subplot(2,3,5);
imagesc(nanmean(cat(3,vfs_cat{:,2}),3));
axis image off;
caxis([-1,1]);
colormap(gca,brewermap([],'PRGn'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('VFS post-muscimol');
colorbar

subplot(2,3,6);
imagesc(nanmean(vfs_change,3));
axis image off;
caxis([-0.5,0.5]);
colormap(gca,brewermap([],'PRGn'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('(post-pre)/(post+pre)');
colorbar


%% ^^^ Striatal spike rate pre/post muscimol

animals = {'AP045','AP054','AP055','AP053','AP047','AP048'};

spike_rate_cond = cell(size(animals));

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworldNoRepeats';
    experiments = AP_find_experiments(animal,protocol);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
        
    disp(animal);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        
        % Load data (first experiment - but spikes throughout used)
        experiment = experiments(curr_day).experiment(1);
        load_parts.ephys = true;
        AP_load_experiment
        
        % Set experiments in conditions 
        % (1-3 pre-muscimol, 4+ post-muscimol)
        % (assumes all repeated expts/failures were post-muscimol)
        curr_experiments = AP_list_experiments(animal,day);
        cond_expts = {[1:3],[4:length(curr_experiments)]};
        
        for curr_cond = 1:2
                exp_starts = sync(2).timestamps(sync(2).values == 1);
                exp_stops = sync(2).timestamps(sync(2).values == 0);
                
                curr_exp_start = exp_starts(cond_expts{curr_cond}(1));
                curr_exp_stop = exp_stops(cond_expts{curr_cond}(end));
                
                curr_use_spikes = spike_times >= curr_exp_start & ...
                    spike_times <= curr_exp_stop & ~isnan(aligned_str_depth_group);
                
                spike_rate_cond{curr_animal}{curr_day}(:,curr_cond) = ...
                    accumarray(aligned_str_depth_group(curr_use_spikes),1, ...
                    [n_aligned_depths,1])./(curr_exp_stop - curr_exp_start);
        end
         
        clearvars -except animals curr_animal animal experiments curr_day ...
            spike_rate_cond
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
    
end

% Concatenate across recordings and get change
spike_rate_cond_cat = cell2mat(reshape([spike_rate_cond{:}],1,1,[]));
spike_rate_cond_cat_change = ...
    squeeze((spike_rate_cond_cat(:,1,:) - spike_rate_cond_cat(:,2,:))./ ...
    (spike_rate_cond_cat(:,1,:)));
n_depths = size(spike_rate_cond_cat,1);

figure;
plotSpread(spike_rate_cond_cat_change','distributionColors',[0.5,0.5,0.5]);
errorbar((1:n_depths)+0.3,nanmean(spike_rate_cond_cat_change,2), ...
    AP_sem(spike_rate_cond_cat_change,2),'.k', ...
    'MarkerSize',20,'linewidth',2,'linestyle','none');
line(xlim,[0,0],'color','r')
xlabel('(post-pre)/(pre)');
ylabel('Striatal depth');
title('Muscimol change');

% (Condition statistics)
disp('Spike rate pre/post muscimol:')
for curr_depth = 1:n_depths
    curr_p = signrank(squeeze(spike_rate_cond_cat(curr_depth,1,:)), ...
        squeeze(spike_rate_cond_cat(curr_depth,2,:)));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]);
end



%% ^^^ Cortex/striatum passive stim pre/post muscimol

data_fns = { ...
    'trial_activity_AP_lcrGratingPassive_pre_muscimol', ...
    'trial_activity_AP_lcrGratingPassive_post_muscimol'};

stimIDs = cell(2,1);
mua_muscimol = cell(2,1);
mua_ctxpred_muscimol = cell(2,1);
fluor_muscimol = cell(2,1);
fluor_roi_muscimol = cell(2,1);
fluor_kernelroi_muscimol = cell(2,1);

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;

    % Split by experiment
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    
    % Get trials with movement during stim to exclude
    wheel_thresh = 0.025;
    quiescent_trials = ~any(abs(wheel_allcat(:,t >= -0.5 & t <= 0.5)) > wheel_thresh,2);
    quiescent_trials_exp = mat2cell(quiescent_trials,trials_recording,1);
    
    % Get stim and activity by experiment
    trial_stim_allcat_exp = mat2cell(trial_stim_allcat,trials_recording,1);  
    fluor_deconv_exp = mat2cell(fluor_allcat_deconv,trials_recording,length(t),n_vs);
    fluor_roi_deconv_exp = mat2cell(fluor_roi_deconv,trials_recording,length(t),numel(wf_roi));
    fluor_kernelroi_deconv_exp = mat2cell(fluor_kernelroi_deconv,trials_recording,length(t),n_depths);
    mua_allcat_exp = mat2cell(mua_allcat,trials_recording,length(t),n_depths);
    mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,trials_recording,length(t),n_depths);
    
    % Exclude trials with fluorescence spikes
    % (this is a dirty way to do this but don't have a better alt)
    fluor_spike_thresh = 0.5; % deconv df/f threshold (eyeballed)
    fluor_spike_trial = cellfun(@(x) any(any(x > fluor_spike_thresh,2),3), ...
        fluor_kernelroi_deconv_exp,'uni',false);
    
    % Grab data
    stimIDs{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        trial_stim_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    mua_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        mua_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    mua_ctxpred_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        mua_ctxpred_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_roi_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_roi_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_kernelroi_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_kernelroi_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
end

% Plot average fluorescence
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');

use_stim = 1;

fluor_premuscimol_mean = ...
    svdFrameReconstruct(U_master(:,:,1:200), ...
    permute(nanmean(cell2mat(cellfun(@(x,stim) ...
    nanmean(x(stim == use_stim,:,:),1),fluor_muscimol{1},stimIDs{1},'uni',false)),1),[3,2,1]));
fluor_postmuscimol_mean = ...
    svdFrameReconstruct(U_master(:,:,1:200), ...
    permute(nanmean(cell2mat(cellfun(@(x,stim) ...
    nanmean(x(stim == use_stim,:,:),1),fluor_muscimol{2},stimIDs{2},'uni',false)),1),[3,2,1]));

AP_image_scroll([fluor_premuscimol_mean,fluor_postmuscimol_mean]);
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5],[],[size(U_master,1),size(U_master,2),1,2]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
axis image;
colormap(brewermap([],'PRGn'));


figure;
t_stim = t >= 0 & t <= 0.2;

subplot(1,2,1)
imagesc(nanmean(fluor_premuscimol_mean(:,:,t_stim),3));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
c = caxis;
axis image off;
colormap(brewermap([],'PRGn'));
title('Pre-muscimol');

subplot(1,2,2)
imagesc(nanmean(fluor_postmuscimol_mean(:,:,t_stim),3));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
caxis(c);
axis image off;
colormap(brewermap([],'PRGn'));
title('Post-muscimol');

% Get pre/post stim response
use_stim = 1;

mua_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_muscimol{1},stimIDs{1},'uni',false));
mua_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_muscimol{2},stimIDs{2},'uni',false));

mua_ctxpred_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_ctxpred_muscimol{1},stimIDs{1},'uni',false));
mua_ctxpred_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_ctxpred_muscimol{2},stimIDs{2},'uni',false));

fluor_roi_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_roi_muscimol{1},stimIDs{1},'uni',false));
fluor_roi_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_roi_muscimol{2},stimIDs{2},'uni',false));

fluor_kernelroi_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_kernelroi_muscimol{1},stimIDs{1},'uni',false));
fluor_kernelroi_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_kernelroi_muscimol{2},stimIDs{2},'uni',false));

%%% CHANGE OPTION 1
% Fit scaling factor (change) pre-post
n_exps = length(stimIDs{1});
str_muscimol_change = nan(n_exps,n_depths);
ctx_muscimol_change = nan(n_exps,n_depths);
for curr_depth = 1:n_depths
    for curr_exp = 1:n_exps
        str_muscimol_change(curr_exp,curr_depth) = ...
            mua_premuscimol_mean(curr_exp,:,curr_depth)'\ ...
            mua_postmuscimol_mean(curr_exp,:,curr_depth)';
        
        ctx_muscimol_change(curr_exp,curr_depth) = ...
            fluor_kernelroi_premuscimol_mean(curr_exp,:,curr_depth)'\ ...
            fluor_kernelroi_postmuscimol_mean(curr_exp,:,curr_depth)';
    end
end


%%% CHANGE OPTION 2
t_stim = t >= 0 & t <= 0.2;
mua_avg_premuscimol = permute(nanmean(mua_premuscimol_mean(:,t_stim,:),2),[1,3,2]);
mua_avg_postmuscimol = permute(nanmean(mua_postmuscimol_mean(:,t_stim,:),2),[1,3,2]);
str_muscimol_change = (mua_avg_postmuscimol-mua_avg_premuscimol);%./(mua_avg_premuscimol);

fluor_avg_premuscimol = permute(nanmean(fluor_kernelroi_premuscimol_mean(:,t_stim,:),2),[1,3,2]);
fluor_avg_postmuscimol = permute(nanmean(fluor_kernelroi_postmuscimol_mean(:,t_stim,:),2),[1,3,2]);
ctx_muscimol_change = (fluor_avg_postmuscimol-fluor_avg_premuscimol);%./(fluor_avg_premuscimol);


% Plot time courses and change
figure;
p = gobjects(n_depths,3);
for plot_str = 1:n_depths
    
    p(plot_str,1) = subplot(n_depths,3,(plot_str-1)*n_depths+1); hold on;
    AP_errorfill(t,nanmean(mua_premuscimol_mean(:,:,plot_str),1)', ...
        AP_sem(mua_premuscimol_mean(:,:,plot_str),1)','k');
    AP_errorfill(t,nanmean(mua_postmuscimol_mean(:,:,plot_str),1)', ...
        AP_sem(mua_postmuscimol_mean(:,:,plot_str),1)','r');
    xlim([-0.2,1])
    xlabel('Time from stim (s)')
    ylabel(['Str ' num2str(plot_str)]);
    axis square
    
    p(plot_str,2) = subplot(n_depths,3,(plot_str-1)*n_depths+2); hold on;
    AP_errorfill(t,nanmean(fluor_kernelroi_premuscimol_mean(:,:,plot_str),1)', ...
        AP_sem(fluor_kernelroi_premuscimol_mean(:,:,plot_str),1)','k');
    AP_errorfill(t,nanmean(fluor_kernelroi_postmuscimol_mean(:,:,plot_str),1)', ...
        AP_sem(fluor_kernelroi_postmuscimol_mean(:,:,plot_str),1)','r');
    xlim([-0.2,1])
    xlabel('Time from stim (s)')
    ylabel('Cortex ROI');
    axis square
    
    p(plot_str,3) = subplot(n_depths,3,(plot_str-1)*n_depths+3);
    plot(ctx_muscimol_change(:,plot_str),str_muscimol_change(:,plot_str),'.k','MarkerSize',20)
    xlabel(['Cortex ROI (post-pre)']);
    ylabel(['Str ' num2str(plot_str) ' (post-pre)']);
    line(xlim,xlim,'color','k','linestyle','--');

end
linkaxes(p(:,1),'xy');
linkaxes(p(:,2),'xy');

% (str/ctx muscimol statistics)
disp('Striatum/cortex muscimol change correlation:')
for curr_depth = 1:n_depths
    [r,p] = corr(str_muscimol_change(:,curr_depth), ...
        ctx_muscimol_change(:,curr_depth), ...
        'rows','complete','type','Pearson');
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(p) ' r = ' num2str(r)]);
end



%% ^^^ Task performance pre/post muscimol

data_fns = { ...
    'trial_activity_vanillaChoiceworldNoRepeats_pre_muscimol', ...
    'trial_activity_vanillaChoiceworldNoRepeats_post_muscimol'};

frac_orient_right = cell(2,1);
rxn_time = cell(2,1);
move_t_hist = cell(2,1);
for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Split data by recording
    trials_allcat = size(wheel_allcat,1);
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    use_split = trials_recording;
    split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
        [1:length(use_split)]',reshape(use_split,[],1),'uni',false));
    
    % Get psychometric
    stim_conditions = unique(trial_stim_allcat);
    [~,stim_idx] = ismember(trial_stim_allcat,stim_conditions,'rows');
    
    trial_stim_idx_allcat_exp = mat2cell(stim_idx,use_split,1);
    trial_choice_allcat_exp = mat2cell(trial_choice_allcat,use_split,1);
    
    frac_orient_right{curr_data} = cell2mat(cellfun(@(stim,choice) ...
        accumarray(stim,choice == -1,[length(stim_conditions),1],@nanmean,NaN), ...
        trial_stim_idx_allcat_exp,trial_choice_allcat_exp,'uni',false)');
    
    % Reaction time by 
    move_t_exp = mat2cell(move_t,use_split,1);
    rxn_time{curr_data} = cell2mat(cellfun(@(stim,rxn) ...
        accumarray(stim,rxn,[length(stim_conditions),1],@nanmedian,NaN), ...
        trial_stim_idx_allcat_exp,move_t_exp,'uni',false)');
    
    % Get histogram of reaction times
    move_t_bins = -0.2:1/sample_rate:1;
    move_t_bin_centers = move_t_bins(1:end-1) + diff(move_t_bins)./2;
    move_t_bin = mat2cell(discretize(move_t,move_t_bins),trials_recording,1);
    
    move_t_hist{curr_data} = cell2mat(cellfun(@(move_t_bin) ...
        accumarray(move_t_bin(~isnan(move_t_bin)), ...
        1/sum(~isnan(move_t_bin)),[length(move_t_bins)-1,1],@nansum,0), ...
        move_t_bin','uni',false));

end

figure; 

subplot(3,2,1);
AP_errorfill(stim_conditions, ...
    [nanmean(frac_orient_right{1},2),nanmean(frac_orient_right{2},2)], ...
    [AP_sem(frac_orient_right{1},2),AP_sem(frac_orient_right{2},2)],[0,0,0;1,0,0]);
line(xlim,[0.5,0.5],'color','k','linestyle','--');
line([0,0],ylim,'color','k','linestyle','--');
xlabel('Contrast*Side');
ylabel('Reaction time');

subplot(3,2,2);
AP_errorfill(stim_conditions, ...
    nanmean(frac_orient_right{2}-frac_orient_right{1},2), ...
    AP_sem(frac_orient_right{2}-frac_orient_right{1},2),[0,0,0;1,0,0]);
line(xlim,[0,0],'color','k','linestyle','--');
line([0,0],ylim,'color','k','linestyle','--');
xlabel('Contrast*Side');
ylabel('\DeltaReaction time');

subplot(3,2,3);
AP_errorfill(stim_conditions, ...
    [nanmean(rxn_time{1},2),nanmean(rxn_time{2},2)], ...
    [AP_sem(rxn_time{1},2),AP_sem(rxn_time{2},2)],[0,0,0;1,0,0]);
line(xlim,[0.5,0.5],'color','k','linestyle','--');
line([0,0],ylim,'color','k','linestyle','--');
xlabel('Contrast*Side');
ylabel('Fraction orient right');

subplot(3,2,4);
AP_errorfill(stim_conditions, ...
    nanmean(rxn_time{2}-rxn_time{1},2), ...
    AP_sem(rxn_time{2}-rxn_time{1},2),[0,0,0;1,0,0]);
line(xlim,[0,0],'color','k','linestyle','--');
line([0,0],ylim,'color','k','linestyle','--');
xlabel('Contrast*Side');
ylabel('\DeltaFraction orient right');

subplot(3,2,5);
AP_errorfill(move_t_bin_centers, ...
    [nanmean(move_t_hist{1},2),nanmean(move_t_hist{2},2)], ...
    [AP_sem(move_t_hist{1},2),AP_sem(move_t_hist{2},2)],[0,0,0;1,0,0]);
line([0.5,0.5],ylim,'color','k','linestyle','--');
xlabel('Reaction time');
ylabel('Fraction');

subplot(3,2,6);
AP_errorfill(move_t_bin_centers, ...
    nanmean(move_t_hist{2}-move_t_hist{1},2), ...
    AP_sem(move_t_hist{2}-move_t_hist{1},2),[0,0,0;1,0,0]);
line(xlim,[0,0],'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
xlabel('Reaction time');
ylabel('\DeltaFraction');


% (Psychometric stim x condition)
curr_stat_data = permute(cat(3,frac_orient_right{:}),[1,3,2]);
[stim_grp,condition_grp,exp_grp] = ndgrid(1:size(curr_stat_data,1), ...
    1:size(curr_stat_data,2),1:size(curr_stat_data,3));
[curr_p,~,~,terms] = anovan(reshape(curr_stat_data,[],1), ...
        [stim_grp(:),condition_grp(:)],'model','interaction','display','off');
 disp(['Psychomatric condition p = ' num2str(curr_p(2))]);
disp(['Psychomatric stim x condition p = ' num2str(curr_p(3))]);

% (Reaction time stim x condition)
curr_stat_data = permute(cat(3,rxn_time{:}),[1,3,2]);
[stim_grp,condition_grp,exp_grp] = ndgrid(1:size(curr_stat_data,1), ...
    1:size(curr_stat_data,2),1:size(curr_stat_data,3));
[curr_p,~,~,terms] = anovan(reshape(curr_stat_data,[],1), ...
        [stim_grp(:),condition_grp(:)],'model','interaction','display','off');
disp(['Reaction time condition p = ' num2str(curr_p(2))]);
disp(['Reaction time stim x condition p = ' num2str(curr_p(3))]);



%% ^^^ Striatal task trial activity pre/post muscimol

data_fns = { ...
    'trial_activity_vanillaChoiceworldNoRepeats_pre_muscimol', ...
    'trial_activity_vanillaChoiceworldNoRepeats_post_muscimol'};

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Split data by recording
    trials_allcat = size(wheel_allcat,1);
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    use_split = trials_recording;
    split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
        [1:length(use_split)]',reshape(use_split,[],1),'uni',false));
    
    % Plot stim-aligned/sorted measured and predicted striatum activity
    % (correct contra trials)
    switch curr_data
        case 1
            cond_name = 'Pre-muscimol';
        case 2
            cond_name = 'Post-muscimol';
    end
    
    for curr_trial_set = 1:2
        switch curr_trial_set
            case 1
                plot_trials = move_t < Inf & trial_stim_allcat > 0 & trial_choice_allcat == -1;
                figure('Name',[cond_name ': Correct contra trials']);
            case 2
                plot_trials = move_t < Inf & trial_stim_allcat < 0 & trial_choice_allcat == 1;
                figure('Name',[cond_name ': Correct ipsi trials']);
        end
        
        p = gobjects(n_depths,4);
        colormap(brewermap([],'Greys'));
        for curr_depth = 1:n_depths
            
            % Get trials to plot, sort by reaction time
            curr_trials = plot_trials & ~all(isnan(mua_allcat(:,:,curr_depth)),2);
            curr_trials_exp = mat2cell(curr_trials,use_split,1);
            curr_trials_idx = find(curr_trials);
            [~,rxn_sort_idx] = sort(move_t(curr_trials_idx));
            
            sorted_plot_trials = curr_trials_idx(rxn_sort_idx);
            
            curr_plot = mua_allcat(sorted_plot_trials,:,curr_depth);
            curr_taskpred_plot = mua_taskpred_allcat(sorted_plot_trials,:,curr_depth);
            curr_ctxpred_plot = mua_ctxpred_allcat(sorted_plot_trials,:,curr_depth);
            
            % Smooth and plot with stim/move/reward times
            % (as conv(nans-zeroed)./conv(non-nan) to ignore in nans in conv)
            smooth_filt = [50,1]; % (trials x frames)
            
            curr_plot_smooth = conv2(curr_plot,ones(smooth_filt),'same')./ ...
                conv2(~isnan(curr_plot),ones(smooth_filt),'same');
            
            curr_taskpred_plot_smooth = curr_taskpred_plot;
            curr_taskpred_plot_smooth(isnan(curr_taskpred_plot_smooth)) = 0;
            curr_taskpred_plot_smooth = conv2(curr_taskpred_plot_smooth,ones(smooth_filt),'same')./ ...
                conv2(~isnan(curr_taskpred_plot),ones(smooth_filt),'same');
            
            curr_ctxpred_plot_smooth = curr_ctxpred_plot;
            curr_ctxpred_plot_smooth(isnan(curr_ctxpred_plot_smooth)) = 0;
            curr_ctxpred_plot_smooth = conv2(curr_ctxpred_plot_smooth,ones(smooth_filt),'same')./ ...
                conv2(~isnan(curr_ctxpred_plot),ones(smooth_filt),'same');
            
            p(curr_depth,1) = subplot(n_depths,4,1+(curr_depth-1)*4,'YDir','reverse'); hold on
            imagesc(t,[],curr_plot_smooth);
            plot(zeros(size(curr_trials_idx)),1:length(curr_trials_idx),'MarkerSize',1,'color','r');
            plot(move_t(sorted_plot_trials),1:length(curr_trials_idx),'MarkerSize',1,'color',[0.8,0,0.8]);
            %         plot(outcome_t(sorted_plot_trials),1:length(curr_trials_idx),'.','MarkerSize',1,'color','b');
            axis tight;
            xlim([-0.2,1]);
            xlabel('Time from stim');
            ylabel('Trials (rxn sorted)');
            title('Measured');
            
            p(curr_depth,2) = subplot(n_depths,4,2+(curr_depth-1)*4,'YDir','reverse'); hold on
            imagesc(t,[],curr_taskpred_plot_smooth);
            plot(zeros(size(curr_trials_idx)),1:length(curr_trials_idx),'MarkerSize',1,'color','r');
            plot(move_t(sorted_plot_trials),1:length(curr_trials_idx),'MarkerSize',1,'color',[0.8,0,0.8]);
            %         plot(outcome_t(sorted_plot_trials),1:length(curr_trials_idx),'.','MarkerSize',1,'color','b');
            axis tight;
            xlim([-0.2,1]);
            xlabel('Time from stim');
            ylabel('Trials (rxn sorted)');
            title('Task-predicted');
            
            p(curr_depth,3) = subplot(n_depths,4,3+(curr_depth-1)*4,'YDir','reverse'); hold on
            imagesc(t,[],curr_ctxpred_plot_smooth);
            plot(zeros(size(curr_trials_idx)),1:length(curr_trials_idx),'MarkerSize',1,'color','r');
            plot(move_t(sorted_plot_trials),1:length(curr_trials_idx),'MarkerSize',1,'color',[0.8,0,0.8]);
            %         plot(outcome_t(sorted_plot_trials),1:length(curr_trials_idx),'.','MarkerSize',1,'color','b');
            axis tight;
            xlim([-0.2,1]);
            xlabel('Time from stim');
            ylabel('Trials (rxn sorted)');
            title('Cortex-predicted');
            
            % Plot average aligned activity
            % (set alignment shifts)
            t_leeway = -t(1);
            leeway_samples = round(t_leeway*(sample_rate));
            stim_align = zeros(size(trial_stim_allcat));
            move_align = -move_idx + leeway_samples;
            outcome_align = -outcome_idx + leeway_samples;
            use_align = {stim_align,move_align,outcome_align};
            
            align_col = [1,0,0;0.8,0,0.8;0,0,0.8];
            % (split the alignment halfway between median alignment points)
            align_median = cellfun(@(x) -nanmedian(x(plot_trials))/sample_rate,use_align);
            align_break = align_median(1:end-1) + diff(align_median*0.8);
            align_t = {[-0.05,align_break(1)],[align_break(1:2)],[align_break(2),1]};
            
            p(curr_depth,4) = subplot(n_depths,4,4+(curr_depth-1)*4); hold on
            for curr_align = 1:length(use_align)
                curr_mua_align = cell2mat(arrayfun(@(trial) circshift(mua_allcat(trial,:,:), ...
                    use_align{curr_align}(trial),2),transpose(1:size(mua_allcat,1)),'uni',false));
                curr_mua_exp = mat2cell(curr_mua_align(:,:,curr_depth),use_split,length(t));
                curr_mua_exp_mean = cell2mat(cellfun(@(data,trials) nanmean(data(trials,:),1), ...
                    curr_mua_exp,curr_trials_exp,'uni',false));
                
                curr_mua_taskpred_align = cell2mat(arrayfun(@(trial) circshift(mua_taskpred_allcat(trial,:,:), ...
                    use_align{curr_align}(trial),2),transpose(1:size(mua_taskpred_allcat,1)),'uni',false));
                curr_mua_taskpred_exp = mat2cell(curr_mua_taskpred_align(:,:,curr_depth),use_split,length(t));
                curr_mua_taskpred_exp_mean = cell2mat(cellfun(@(data,trials) nanmean(data(trials,:),1), ...
                    curr_mua_taskpred_exp,curr_trials_exp,'uni',false));
                
                curr_mua_ctxpred_align = cell2mat(arrayfun(@(trial) circshift(mua_ctxpred_allcat(trial,:,:), ...
                    use_align{curr_align}(trial),2),transpose(1:size(mua_ctxpred_allcat,1)),'uni',false));
                curr_mua_ctxpred_exp = mat2cell(curr_mua_ctxpred_align(:,:,curr_depth),use_split,length(t));
                curr_mua_ctxpred_exp_mean = cell2mat(cellfun(@(data,trials) nanmean(data(trials,:),1), ...
                    curr_mua_ctxpred_exp,curr_trials_exp,'uni',false));
                
                curr_t_offset = -nanmedian(use_align{curr_align}(plot_trials))/sample_rate;
                curr_t = t + curr_t_offset;
                curr_t_plot = curr_t >= align_t{curr_align}(1) & ...
                    curr_t <= align_t{curr_align}(2);
                
                plot_t = curr_t > align_t{curr_align}(1) & curr_t <= align_t{curr_align}(2);
                
                AP_errorfill(curr_t(plot_t'), ...
                    nanmean(curr_mua_exp_mean(:,plot_t),1)', ...
                    AP_sem(curr_mua_exp_mean(:,plot_t),1)','k');
                
                AP_errorfill(curr_t(plot_t'), ...
                    nanmean(curr_mua_taskpred_exp_mean(:,plot_t),1)', ...
                    AP_sem(curr_mua_taskpred_exp_mean(:,plot_t),1)','b');
                
                AP_errorfill(curr_t(plot_t'), ...
                    nanmean(curr_mua_ctxpred_exp_mean(:,plot_t),1)', ...
                    AP_sem(curr_mua_ctxpred_exp_mean(:,plot_t),1)',[0,0.7,0]);
                
                line(repmat(curr_t_offset,2,1),ylim,'color',align_col(curr_align,:));
            end
            xlabel('~Time from stim');
            ylabel('Striatum depth');
            
        end
        % Link the x-axes, set the c/y-axes same within a row
        linkaxes(p(:),'x');
        
        for curr_row = 1:size(p,1)
            curr_ylim = ylim(p(curr_row,4));
            caxis(p(curr_row,1),[0,curr_ylim(2)]);
            caxis(p(curr_row,2),[0,curr_ylim(2)]);
            caxis(p(curr_row,3),[0,curr_ylim(2)]);
        end
        
        trial_scale = 500;
        t_scale = 0.5;
        y_scale = 1;
        line(p(1,1),min(xlim(p(1,1))) + [0,t_scale],repmat(min(ylim(p(1,1))),2,1),'color','k','linewidth',3);
        line(p(1,4),min(xlim(p(1,4))) + [0,t_scale],repmat(min(ylim(p(1,4))),2,1),'color','k','linewidth',3);
        line(p(1,1),repmat(min(xlim(p(1,1))),2,1),min(ylim(p(1,1))) + [0,trial_scale],'color','k','linewidth',3);
        line(p(1,4),repmat(min(xlim(p(1,4))),2,1),min(ylim(p(1,4))) + [0,y_scale],'color','k','linewidth',3);
        drawnow;
        
    end
    
end


%% ^^^ Striatal task kernels pre/post muscimol

data_fns = { ...
    'trial_activity_vanillaChoiceworldNoRepeats_pre_muscimol', ...
    'trial_activity_vanillaChoiceworldNoRepeats_post_muscimol'};

mua_norm_exp = cell(2,1);
task_str_kernel = cell(2,1);
task_str_ctxpred_kernel = cell(2,1);
task_ctx_roi_kernel = cell(2,1);

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Keep task > str kernels, task > ctx-pred str norm factor
    mua_norm_exp{curr_data} = vertcat(mua_norm{:});
    task_str_kernel{curr_data} = vertcat(mua_taskpred_k_all{:});
    task_str_ctxpred_kernel{curr_data} = vertcat(mua_ctxpred_taskpred_k_all{:});
    
end

% Normalize and concatenate task kernels
mua_task_k_premuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_kernel{1},mua_norm_exp{1},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

mua_ctxpred_task_k_premuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_ctxpred_kernel{1},mua_norm_exp{1},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

mua_task_k_postmuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_kernel{2},mua_norm_exp{2},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

mua_ctxpred_task_k_postmuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_ctxpred_kernel{2},mua_norm_exp{2},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

% Get task>striatum parameters
n_regressors = length(task_regressor_labels);

% Plot task>striatum kernels
stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = [];
move_col = [0.6,0,0.6;0,0.6,0];
go_col = [0.8,0.8,0.2;0.5,0.5,0.5];
outcome_col = [0,0,0.8;0.8,0,0];
task_regressor_cols = {stim_col,move_col,go_col,outcome_col};
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

figure('Name','Pre-muscimol');
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);
        
        curr_kernels = mua_task_k_premuscimol{curr_regressor}(:,:,curr_depth,:);
        n_subregressors = size(mua_task_k_premuscimol{curr_regressor},1);
        col = task_regressor_cols{curr_regressor};
        for curr_subregressor = 1:n_subregressors
            AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
                nanmean(curr_kernels(curr_subregressor,:,:,:),4)', ...
                AP_sem(curr_kernels(curr_subregressor,:,:,:),4)', ...
                col(curr_subregressor,:));
        end
        
        xlabel('Time (s)');
        ylabel('Weight');
        title(task_regressor_labels{curr_regressor});
        line([0,0],ylim,'color','k');
        
    end
end
linkaxes(p);
y_scale = 1;
t_scale = 0.5;
line(min(xlim) + [0,t_scale],repmat(min(ylim),2,1),'color','k','linewidth',3);
line(repmat(min(xlim),2,1),min(ylim) + [0,y_scale],'color','k','linewidth',3);

figure('Name','Post-muscimol');
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);
        
        curr_kernels = mua_task_k_postmuscimol{curr_regressor}(:,:,curr_depth,:);
        n_subregressors = size(mua_task_k_postmuscimol{curr_regressor},1);
        col = task_regressor_cols{curr_regressor};
        for curr_subregressor = 1:n_subregressors
            AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
                nanmean(curr_kernels(curr_subregressor,:,:,:),4)', ...
                AP_sem(curr_kernels(curr_subregressor,:,:,:),4)', ...
                col(curr_subregressor,:),0.5);
        end
        
        xlabel('Time (s)');
        ylabel('Weight');
        title(task_regressor_labels{curr_regressor});
        line([0,0],ylim,'color','k');
        
    end
end
linkaxes(p);
y_scale = 1;
t_scale = 0.5;
line(min(xlim) + [0,t_scale],repmat(min(ylim),2,1),'color','k','linewidth',3);
line(repmat(min(xlim),2,1),min(ylim) + [0,y_scale],'color','k','linewidth',3);

% Plot kernel sums pre/post muscimol
str_k_sum_premuscimol = cellfun(@(x) permute(sum(x,2),[1,3,4,2]),mua_task_k_premuscimol,'uni',false);
str_k_sum_postmuscimol = cellfun(@(x) permute(sum(x,2),[1,3,4,2]),mua_task_k_postmuscimol,'uni',false);

figure;
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors
        if curr_regressor == 1
            x = unique([0.06,0.125,0.25,0.5,1].*[-1;1]);
        else
            x = 1:size(str_k_sum_premuscimol{curr_regressor},1);
        end
        
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);
        hold on
        
        errorbar(x,nanmean(str_k_sum_premuscimol{curr_regressor}(:,curr_depth,:),3), ...
            AP_sem(str_k_sum_premuscimol{curr_regressor}(:,curr_depth,:),3),'k','linewidth',2);
        
        errorbar(x,nanmean(str_k_sum_postmuscimol{curr_regressor}(:,curr_depth,:),3), ...
            AP_sem(str_k_sum_postmuscimol{curr_regressor}(:,curr_depth,:),3),'r','linewidth',2);
        
        axis tight
        xlim(xlim + [-0.2,0.2]);
        
        xlabel('Condition');
        ylabel('Weight sum');
        title(task_regressor_labels{curr_regressor});
        
    end
end

linkaxes(p,'y');

% (Regressor weight sum by condition statistics)
disp('Regressor v condition 2-way anova (only disp condition effect):')
for curr_depth = 1:n_depths
    for curr_regressor = 1:length(task_regressor_labels)
        
        curr_prepost = permute(cat(4, ...
            str_k_sum_premuscimol{curr_regressor}(:,curr_depth,:), ...
            str_k_sum_postmuscimol{curr_regressor}(:,curr_depth,:)),[1,3,4,2]);
        
        [regressor_grp,exp_grp,condition_grp] = ndgrid( ...
            1:size(curr_prepost,1),1:size(curr_prepost,2),1:size(curr_prepost,3));
        
        curr_p = anovan(reshape(curr_prepost,[],1), ...
            [regressor_grp(:),condition_grp(:)],'display','off');
        
        disp(['Str ' num2str(curr_depth) ' ' task_regressor_labels{curr_regressor} ...
            ' p = ' num2str(curr_p(2)')]);
        
    end
end



%% ^^^ Striatal task/cortex explained variance pre/post muscimol

data_fns = { ...
    'trial_activity_vanillaChoiceworldNoRepeats_pre_muscimol', ...
    'trial_activity_vanillaChoiceworldNoRepeats_post_muscimol'};

taskpred_r2_condition = [];
ctxpred_r2_condition = [];
for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
      
    % Use raw data (not normalized or baseline-subtracted) for expl var
    mua_exp = vertcat(mua_all{:});
    mua_taskpred_exp = vertcat(mua_taskpred_all{:});
    mua_ctxpred_exp = vertcat(mua_ctxpred_all{:});
    
    % Get R^2 for task and cortex
    taskpred_r2 = nan(length(mua_exp),n_depths);
    ctxpred_r2 = nan(length(mua_exp),n_depths);
    for curr_exp = 1:length(mua_exp)
        
        curr_str = reshape(permute(mua_exp{curr_exp},[2,1,3]),[],n_depths);
        curr_str_baselinesub = reshape(permute(mua_exp{curr_exp},[2,1,3]),[],n_depths) - ...
            (nanmean(reshape(mua_exp{curr_exp}(:,t < 0,:),[],size(mua_exp{curr_exp},3)),1));
        curr_taskpred = reshape(permute(mua_taskpred_exp{curr_exp},[2,1,3]),[],n_depths);
        curr_ctxpred = reshape(permute(mua_ctxpred_exp{curr_exp},[2,1,3]),[],n_depths);
        
        % Set common NaNs
        nan_samples = isnan(curr_str) | isnan(curr_str_baselinesub) | ...
            isnan(curr_taskpred) | isnan(curr_ctxpred);
        curr_str(nan_samples) = NaN;
        curr_str_baselinesub(nan_samples) = NaN;
        curr_taskpred(nan_samples) = NaN;
        curr_ctxpred(nan_samples) = NaN;
        
        % (task regressed from average baseline-subtracted data)
        taskpred_r2(curr_exp,:) = 1 - (nansum((curr_str_baselinesub-curr_taskpred).^2,1)./ ...
            nansum((curr_str_baselinesub-nanmean(curr_str_baselinesub,1)).^2,1));
        % (cortex regressed from raw data)
        ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_str-curr_ctxpred).^2,1)./ ...
            nansum((curr_str-nanmean(curr_str,1)).^2,1));
        
    end
    
    taskpred_r2_condition(:,:,curr_data) = taskpred_r2;
    ctxpred_r2_condition(:,:,curr_data) = ctxpred_r2;
    
end

% Plot explained variance task vs cortex by experiment
figure; hold on;
line([0,1],[0,1],'color','k','linestyle','--');

str_col = copper(n_depths);
for curr_depth = 1:n_depths
    
    errorbar(squeeze(nanmean(taskpred_r2_condition(:,curr_depth,:),1)), ...
        squeeze(nanmean(ctxpred_r2_condition(:,curr_depth,:),1)), ...
        squeeze(AP_sem(ctxpred_r2_condition(:,curr_depth,:),1)),squeeze(AP_sem(ctxpred_r2_condition(:,curr_depth,:),1)), ...
        squeeze(AP_sem(taskpred_r2_condition(:,curr_depth,:),1)),squeeze(AP_sem(taskpred_r2_condition(:,curr_depth,:),1)), ...
        'color',str_col(curr_depth,:),'linewidth',2);
    
   scatter(nanmean(taskpred_r2_condition(:,curr_depth,1),1), ...
       nanmean(ctxpred_r2_condition(:,curr_depth,1),1),100, ...
       str_col(curr_depth,:),'filled','MarkerEdgeColor',[0,0.7,0],'linewidth',2);
   scatter(nanmean(taskpred_r2_condition(:,curr_depth,2),1), ...
       nanmean(ctxpred_r2_condition(:,curr_depth,2),1),100, ...
       str_col(curr_depth,:),'filled','MarkerEdgeColor',[0.7,0,0],'linewidth',2);
end
xlabel('Task R^2');
ylabel('Cortex R^2');

% (Cortex vs task R2 statistics)
disp('Cortex-task R^2 pre/post-muscimol');
ctx_task_r2_diff = ctxpred_r2_condition - taskpred_r2_condition;
for curr_depth = 1:n_depths
    curr_p = signrank(ctx_task_r2_diff(:,curr_depth,1), ...
        ctx_task_r2_diff(:,curr_depth,2));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end



