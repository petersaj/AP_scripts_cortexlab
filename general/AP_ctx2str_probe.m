% Get Allen projections from the cortex to striatum using probe location

%% Load in probe points
% (skip this if already loaded in from allen_atlas_probe)

% Load saved points
animal = 'AP034';
[probe_filename,probe_filename_exists] = AP_cortexlab_filename(animal,[],[],'probe_histology');
load(probe_filename);

% Get probe vector
histology_points = pointList.pointList{1};
r0 = mean(histology_points,1);
xyz = bsxfun(@minus,histology_points,r0);
[~,~,V] = svd(xyz,0);
histology_probe_direction = V(:,1);

probe_vector_evaluate = [1000,-1000];
probe_vector = bsxfun(@plus,bsxfun(@times,probe_vector_evaluate',histology_probe_direction'),r0);
probe_vector_ccf = round(probe_vector(:,[3,2,1]));

%% Get probe depths per micron within the atlas

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

%% Plot the probe trajectory and Allen atlas in native orientation
% (note the CCF is rotated to allow for dim 1 = x)

figure; hold on
slice_spacing = 10;
target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) > 1,[2,1,3]);
structure_patch = isosurface(target_volume,0);
structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
target_structure_color = [0.7,0.7,0.7];
brain_outline = patch('Vertices',structure_wire.vertices*slice_spacing, ...
    'Faces',structure_wire.faces, ...
    'FaceColor','none','EdgeColor',target_structure_color);

axis image vis3d;
view([-30,25]);
rotate3d;

scatter3(probe_depths(:,1),probe_depths(:,2),probe_depths(:,3),5,'k');
drawnow;

%% Get the Allen projection sites in regular intervals in probe striatum
% NOTE: update this later to use aligned depths instead of even segments

% Convert probe depths subscripts to indicies
probe_depths_ind = sub2ind(size(av),probe_depths(:,1),probe_depths(:,2),probe_depths(:,3));

% Get structures that the probe is in
probe_structures = av(probe_depths_ind);

% Get regular depths from probe that runs through striatum
str_id = find(strcmp(st.safe_name,'Caudoputamen'));
probe_structures_str = probe_structures == str_id;

n_str_depths = 4;
probe_str = [probe_depths(find(probe_structures_str,1,'first'),:); ...
    probe_depths(find(probe_structures_str,1,'last'),:)];
str_depths_um = ...
    round([linspace(probe_str(1,1)',probe_str(2,1)',n_str_depths); ...
    linspace(probe_str(1,2)',probe_str(2,2)',n_str_depths); ...
    linspace(probe_str(1,3)',probe_str(2,3)',n_str_depths)]')*10;

% Plot striatal area on the figure;
scatter3(probe_depths(probe_structures_str,1), ...
    probe_depths(probe_structures_str,2), ...
    probe_depths(probe_structures_str,3),50,jet(sum(probe_structures_str,1)));
drawnow;

% Mirror all of the locations across the midline and get projections from
% both (because the cortical injections in the database aren't evenly
% distributed, so this is more comprehensive)
bregma_um = allenCCFbregma*10;
ccf_midline = bregma_um(3);
hemisphere = sign(str_depths_um(1,3) - ccf_midline);
str_depths_mirror = [str_depths_um(:,1:2),ccf_midline - hemisphere*abs(str_depths_um(:,3)-ccf_midline)];

max_sites = 50;
str_depths_query = [str_depths_um;str_depths_mirror];
injection_parameters = get_allen_projection(str_depths_query,max_sites);
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

% Get relative projection density / injection volumes
% projection_strength = cellfun(@(density,volume) density./volume, ...
%     {injection_parameters.density},{injection_parameters.volume},'uni',false);
% (or currently using: just the density, not sure whether good to norm)
projection_strength = cellfun(@(density,volume) density, ...
    {injection_parameters.density},{injection_parameters.volume},'uni',false);
projection_strength_normalized = cellfun(@(x) mat2gray(x, ...
    [min([projection_strength{:}]),max([projection_strength{:}])]), ...
    projection_strength,'uni',false);

% Plot the injection sites
figure;
hold on; axis image off;
cmap = jet;
plot_colors = repmat(cmap(round(linspace(1,size(cmap,1),size(str_depths_um,1))),:),2,1);
AP_reference_outline('ccf','k');
for curr_coord = 1:length(injection_coordinates)
    scatter(injection_coordinates_standardized{curr_coord}(:,3)-bregma_um(3), ...
        bregma_um(1)-injection_coordinates_standardized{curr_coord}(:,1), ...
        [projection_strength_normalized{curr_coord}*50 + 10], ...
        plot_colors(curr_coord,:),'filled');    
end
text(-4000,3000,'Ipsi','HorizontalAlignment','center')
text(4000,3000,'Contra','HorizontalAlignment','center')
title('Injection sites')

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%% Estimate probe location from average widefield image

probe_angle = 45; % from horizontal

animal = 'AP034';
day = '2018-10-26';
[img_path,img_exists] = AP_cortexlab_filename(animal,day,[],'imaging');
avg_im = readNPY([img_path filesep 'meanImage_blue.npy']);

avg_im_aligned = AP_align_widefield(animal,day,avg_im);

h = figure;imagesc(avg_im_aligned);
axis image off;
colormap(gray);caxis([0,5000]);
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
figure; 

subplot(1,2,1);
imagesc(avg_im_aligned);
axis image off;
colormap(gray);caxis([0,5000]);
AP_reference_outline('ccf_aligned','r');
line(probe_wf(:,1),probe_wf(:,2),'color','b','linewidth',2,'linestyle','--');
title('Widefield');

subplot(1,2,2); hold on; set(gca,'YDir','reverse');axis image;
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

% Get estimated probe vector
r0 = mean(probe_ccf,1);
xyz = bsxfun(@minus,probe_ccf,r0);
[~,~,V] = svd(xyz,0);
probe_direction = V(:,1);

probe_vector_evaluate = [1000,-1000];
probe_vector_ccf = round(bsxfun(@plus,bsxfun(@times,probe_vector_evaluate',probe_direction'),r0));

% Plot estimated probe location on CCF
% (note the CCF is rotated to allow for dim 1 = x)
figure; hold on
slice_spacing = 10;
target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) > 1,[2,1,3]);
structure_patch = isosurface(target_volume,0);
structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
target_structure_color = [0.7,0.7,0.7];
brain_outline = patch('Vertices',structure_wire.vertices*slice_spacing, ...
    'Faces',structure_wire.faces, ...
    'FaceColor','none','EdgeColor',target_structure_color);

axis image vis3d;
view([-30,25]);
rotate3d;

line(probe_ccf(:,1),probe_ccf(:,2),probe_ccf(:,3),'color','r','linewidth',2);
line(probe_vector_ccf(:,1),probe_vector_ccf(:,2),probe_vector_ccf(:,3),'color','b','linewidth',1);

















