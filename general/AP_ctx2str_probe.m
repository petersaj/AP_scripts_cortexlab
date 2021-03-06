% Get Allen projections from the cortex to striatum using probe location

%% Load in probe points
% (skip this if already loaded in from allen_atlas_probe)

% Load saved points
animal = 'AP032';
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

h = figure; hold on
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
cameratoolbar(h,'SetCoordSys','y');
cameratoolbar(h,'SetMode','orbit');

scatter3(probe_depths(:,1),probe_depths(:,2),probe_depths(:,3),15,'k','filled');
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

str_depths_query = [str_depths_um;str_depths_mirror];
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

%% ~~~~~~~~ HISTOLOGY + WIDEFIELD-ALIGNED ~~~~~~~~~~~


%% Compare probe location from average widefield image and histology

probe_angle = 45; % from horizontal

animal = 'AP032';
day = '2018-10-26';
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
figure('Name',animal); 

subplot(1,2,1);
imagesc(avg_im_aligned);
axis image off;
colormap(gray);caxis([0,prctile(avg_im_aligned(:),95)]);
AP_reference_outline('ccf_aligned','r');
line(probe_wf(:,1),probe_wf(:,2),'color','b','linewidth',1,'linestyle','--');
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
h = figure('Name',animal); ccf_axes = axes; hold on
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
title([animal ' ' day]);

% Plot current probe vector
probe_wf_line = line(ccf_axes,probe_vector_ccf(:,1),probe_vector_ccf(:,2),probe_vector_ccf(:,3),'color','k','linewidth',3);
% (make sure histology probe enters on the left)
if  diff(probe_vector_ccf_histology(:,3)) < 0
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


legend([probe_wf_line,probe_histology_line],{'Widefield-estimated','Histology'},'location','nw');


    
%% Directly compare regression and projection maps (from histology)

%%% Define and load
animal = 'AP032';
day = '2018-10-26';
experiment = 1;
verbose = false;

str_align = 'kernel';
AP_load_experiment;

%%% Load lambda from previously estimated and saved
lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
load(lambda_fn);
curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
if any(curr_animal_idx)
    curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
    if any(curr_day_idx)
        lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
    end
end

warning('TEMPORARY LAMBDA')
lambda = 3;

%%% Prepare data for regression
upsample_factor = 2;
sample_rate = (1/median(diff(frame_t)))*upsample_factor;

skip_seconds = 60;
time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Get upsampled dVdf's
use_svs = 1:50;
dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
    diff(fVdf(use_svs,:),[],2)',time_bin_centers)';

% Get striatum depth group by across-experiment alignment
n_depths = n_aligned_depths;
depth_group = aligned_str_depth_group;

binned_spikes = zeros(n_depths,length(time_bin_centers));
for curr_depth = 1:n_depths
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
end

% Get rid of NaNs (if no data?)
binned_spikes(isnan(binned_spikes)) = 0;

kernel_t = [-0.3,0.3];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
zs = [false,true];
cvfold = 10;

%%% Regress MUA from cortex
[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(dfVdf_resample, ...
    binned_spikes,kernel_frames,lambda,zs,cvfold);

% Reshape kernel and convert to pixel space
r = reshape(k,length(use_svs),length(kernel_frames),size(binned_spikes,1));

aUdf = single(AP_align_widefield(Udf,animal,day));
r_px = zeros(size(aUdf,1),size(aUdf,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3)
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(aUdf(:,:,use_svs),r(:,:,curr_spikes));
end

depth_kernels = squeeze(r_px(:,:,kernel_frames == 0,:));

% Load in probe points from histology
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

% Get probe depths per micron within the atlas

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

% Plot the probe trajectory and Allen atlas in native orientation
% (note the CCF is rotated to allow for dim 1 = x)

h = figure; hold on
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
cameratoolbar(h,'SetCoordSys','y');
cameratoolbar(h,'SetMode','orbit');

scatter3(probe_depths(:,1),probe_depths(:,2),probe_depths(:,3),5,'k');
drawnow;

% Get the Allen projection sites in regular intervals in probe striatum
% NOTE: update this later to use aligned depths instead of even segments

% Convert probe depths subscripts to indicies
probe_depths_ind = sub2ind(size(av),probe_depths(:,1),probe_depths(:,2),probe_depths(:,3));

% Get structures that the probe is in
probe_structures = av(probe_depths_ind);

% Get probe in striatum according to CCF
str_id = find(strcmp(st.safe_name,'Caudoputamen'));
probe_structures_str = probe_structures == str_id;
str_depth_ccf = [probe_depths(find(probe_structures_str,1,'first'),:); ...
    probe_depths(find(probe_structures_str,1,'last'),:)];

% OPTION 1: use query points as relative absolute depth in um
%
% [~,idx,~] = unique(spike_templates);
% template_aligned_depth = aligned_str_depth_group(idx)+1;
% template_aligned_depth(isnan(template_aligned_depth)) = 1;
% 
% for curr_depth = 1:n_aligned_depths
%     depth_top = min(spike_depths(aligned_str_depth_group == curr_depth));
%     depth_bottom = max(spike_depths(aligned_str_depth_group == curr_depth));    
%     aligned_depth_relative_ccf = ([depth_top,depth_bottom] - str_depth(2))/10;
% end
% % SANITY CHECK THIS WITH A PLOT
% 
% 
% % % Get striatum length on probe in CCF vs ephys
% % str_length_ccf = pdist2(str_depth_ccf(1,:),str_depth_ccf(2,:))*10;
% % str_length_ephys = diff(str_depth);
% 
% % Get probe vector in striatum using the striatum end
% r0 = str_depth_ccf(2,:);
% xyz = bsxfun(@minus,str_depth_ccf,r0);
% [~,~,V] = svd(xyz,0);
% ccf_probe_direction = V(:,1);
% 
% % Get section of CCF probe corresponding to ephys by distance
% probe_section_evaluate = [aligned_depth_relative_ccf(1),aligned_depth_relative_ccf(2)];
% probe_section_ccf = round(bsxfun(@plus,bsxfun(@times,probe_section_evaluate',ccf_probe_direction'),r0));
% probe_query_um = mean(probe_section_ccf,1)*10;


% OPTION 2: don't use um, instead normalize str distance from what it is to
% what it should be and divide evenly (this won't work if part of the
% striatum isn't recorded from, i.e. the top part is messed up, this seems
% unlikely on the first day though)

[~,idx,~] = unique(spike_templates);
template_aligned_depth = aligned_str_depth_group(idx)+1;
template_aligned_depth(isnan(template_aligned_depth)) = 1;

depth_str_um = nan(n_aligned_depths,3);
for curr_depth = 1:n_aligned_depths
    % (skip if there aren't units in this depth)
    if ~any(aligned_str_depth_group == curr_depth)
        continue
    end
    depth_top = min(spike_depths(aligned_str_depth_group == curr_depth));
    depth_bottom = max(spike_depths(aligned_str_depth_group == curr_depth));  
    depth_center = mean([depth_top,depth_bottom]);
    
    depth_str_relative = (depth_center-str_depth(1))/diff(str_depth);    
    depth_str_um(curr_depth,:) = round((str_depth_ccf(1,:) + diff(str_depth_ccf,[],1)*depth_str_relative)*10);
end
depth_str_um(all(isnan(depth_str_um),2),:) = [];

% Mirror all of the locations across the midline and get projections from
% both (because the cortical injections in the database aren't evenly
% distributed, so this is more comprehensive)
bregma_um = allenCCFbregma*10;
ccf_midline = bregma_um(3);
hemisphere = sign(depth_str_um(1,3) - ccf_midline);
str_depths_mirror = [depth_str_um(:,1:2),ccf_midline - hemisphere*abs(depth_str_um(:,3)-ccf_midline)];

str_depths_query = [depth_str_um;str_depths_mirror];
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
plot_colors = repmat(cmap(round(linspace(1,size(cmap,1),size(depth_str_um,1))),:),2,1);
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

% Convert points from CCF to widefield and plot overlay
ccf_tform_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\ccf_tform'];
load(ccf_tform_fn);

um2pixel = 20.6;
injection_coordinates_wf = cellfun(@(x) ...
    [x(:,[3,1]).*(1/(um2pixel)),ones(size(x,1),1)]*ccf_tform.T, ...
    injection_coordinates_standardized,'uni',false);

figure('Name',[animal ' ' day ': histology-aligned']); 
colormap(brewermap([],'*RdBu'));
used_depths = unique(aligned_str_depth_group(~isnan(aligned_str_depth_group)));
for curr_depth = 1:n_aligned_depths
    subplot(1,n_aligned_depths,curr_depth); hold on; axis image off;
    set(gca,'YDir','reverse');
    
    imagesc(depth_kernels(:,:,curr_depth));
    caxis([-prctile(abs(depth_kernels(:)),99),prctile(abs(depth_kernels(:)),99)]);
    
    curr_inj_idx = find(used_depths == curr_depth);
    curr_mirror_inj_idx = length(used_depths) + curr_inj_idx;
    if ~isempty(curr_inj_idx)
        % (plot points from both hemispheres)
        scatter(injection_coordinates_wf{curr_inj_idx}(:,1), ...
            injection_coordinates_wf{curr_inj_idx}(:,2), ...
            [projection_strength_normalized{curr_inj_idx}*50 + 10], ...
            'k','filled');
        scatter(injection_coordinates_wf{curr_mirror_inj_idx}(:,1), ...
            injection_coordinates_wf{curr_mirror_inj_idx}(:,2), ...
            [projection_strength_normalized{curr_mirror_inj_idx}*50 + 10], ...
            'k','filled');        
    end    
    AP_reference_outline('ccf_aligned','k');
end

%% ~~~~~~~~ WIDEFIELD-ALIGNED ONLY ~~~~~~~~~~~


%% Estimate probe location from average widefield image in batch

probe_angle = 45; % from horizontal

protocol = 'vanillaChoiceworld';
flexible_name = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029', ... % Original 
    'AP045','AP054','AP055','AP053','AP047','AP048', ...        % Muscimol
    'AP043','AP060','AP061'};                                   % Ctx ephys

% Load in the annotated Allen volume and names
allen_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
av = readNPY([allen_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean

% Plot brain to overlay probes
% (note the CCF is rotated to allow for dim 1 = x)
h = figure; 
ccf_axes = subplot(2,2,1); hold on
sagittal_axes = subplot(2,2,2,'YDir','reverse'); hold on;
axis image off;
horizontal_axes = subplot(2,2,3,'YDir','reverse'); hold on;
axis image off;
coronal_axes = subplot(2,2,4,'YDir','reverse'); hold on;
axis image off;

% Plot 1 = 3D
% (Use wire mesh - can add other structures)
slice_spacing = 10;

% target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) > 1,[2,1,3]);
% structure_patch = isosurface(target_volume,0);
% structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
% target_structure_color = [0.7,0.7,0.7];
% brain_outline = patch('Vertices',structure_wire.vertices*slice_spacing, ...
%     'Faces',structure_wire.faces, ...
%     'FaceColor','none','EdgeColor',target_structure_color);

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
view([-30,25]);
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
axis image off;

plot(horizontal_axes,horizontal_outline{1}(:,2),horizontal_outline{1}(:,1),'k','linewidth',2);
plot(horizontal_axes,str_horizontal_outline{1}(:,2),str_horizontal_outline{1}(:,1),'b','linewidth',2);
plot(horizontal_axes,str_horizontal_outline{2}(:,2),str_horizontal_outline{2}(:,1),'b','linewidth',2);
axis image off;

plot(coronal_axes,coronal_outline{1}(:,2),coronal_outline{1}(:,1),'k','linewidth',2);
plot(coronal_axes,str_coronal_outline{1}(:,2),str_coronal_outline{1}(:,1),'b','linewidth',2);
plot(coronal_axes,str_coronal_outline{2}(:,2),str_coronal_outline{2}(:,1),'b','linewidth',2);
axis image off;

if length(animals) <= 17
    probe_color = [brewermap(8,'Set1');brewermap(8,'Set1')];
else
    error('More animals than colors: add another set');
end

probe_vector_ccf = cell(length(animals),1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        
        [img_path,img_exists] = AP_cortexlab_filename(animal,day,[],'imaging');
        avg_im = readNPY([img_path filesep 'meanImage_blue.npy']);
        
        avg_im_aligned = AP_align_widefield(avg_im,animal,day);
        
        h = figure('units','normalized','outerposition',[0 0 1 1]);
        imagesc(avg_im_aligned);
        axis image off;
        colormap(gray);caxis([0,prctile(avg_im_aligned(:),95)]);
        title(['Click probe bottom/top - ' ...
            num2str(curr_animal) '/' num2str(length(animals)) ': ' ...
            num2str(curr_day) '/' num2str(length(experiments))]);
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
        
        %     % Plot probe points on CCF next to average image
        %     figure;
        %
        %     subplot(1,2,1);
        %     imagesc(avg_im_aligned);
        %     axis image off;
        %     colormap(gray);caxis([0,5000]);
        %     AP_reference_outline('ccf_aligned','r');
        %     line(probe_wf(:,1),probe_wf(:,2),'color','b','linewidth',2,'linestyle','--');
        %     title('Widefield');
        %
        %     subplot(1,2,2); hold on; set(gca,'YDir','reverse');axis image;
        %     load('C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\cortical_area_boundaries.mat');
        %     for curr_area_idx =1:length(cortical_area_boundaries)
        %         p = cellfun(@(outline) plot(outline(:,2),outline(:,1),'color','k'), ...
        %             cortical_area_boundaries{curr_area_idx},'uni',false);
        %     end
        %     line(probe_ccf(:,1),probe_ccf(:,2),'color','r','linewidth',2);
        %     title('CCF');
        
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
        
        probe_vector_evaluate = [0,sign(probe_direction(2))*1000];
        probe_vector_draw = round(bsxfun(@plus,bsxfun(@times,probe_vector_evaluate',probe_direction'),r0));
        
        % Plot current probe vector
        line(ccf_axes,probe_vector_draw(:,1),probe_vector_draw(:,2),probe_vector_draw(:,3), ...
            'color',probe_color(curr_animal,:),'linewidth',2);
        line(sagittal_axes,probe_vector_draw(:,1),probe_vector_draw(:,2), ...
            'color',probe_color(curr_animal,:),'linewidth',2);
        line(horizontal_axes,probe_vector_draw(:,1),probe_vector_draw(:,3), ...
            'color',probe_color(curr_animal,:),'linewidth',2);
        line(coronal_axes,probe_vector_draw(:,3),probe_vector_draw(:,2), ...
            'color',probe_color(curr_animal,:),'linewidth',2);
        
        drawnow;
        
        % Store probe vector
        probe_vector_ccf{curr_animal}(:,curr_day) = probe_direction;
        
    end
end



%% Directly compare regression and projection maps (from widefield avg)

probe_angle = 45; % from horizontal

animal = 'AP034';
day = '2018-10-26';
experiment = 1;
verbose = false;
str_align = 'kernel';
[img_path,img_exists] = AP_cortexlab_filename(animal,day,[],'imaging');
avg_im = readNPY([img_path filesep 'meanImage_blue.npy']);

avg_im_aligned = AP_align_widefield(avg_im,animal,day);

h = figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(avg_im_aligned);
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

%%% Load full data and do regression
AP_load_experiment;

%%% Load lambda from previously estimated and saved
lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
load(lambda_fn);
curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
if any(curr_animal_idx)
    curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
    if any(curr_day_idx)
        lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
    end
end

warning('TEMPORARY LAMBDA')
lambda = 3;

%%% Prepare data for regression
upsample_factor = 2;
sample_rate = (1/median(diff(frame_t)))*upsample_factor;

skip_seconds = 60;
time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Get upsampled dVdf's
use_svs = 1:50;
dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
    diff(fVdf(use_svs,:),[],2)',time_bin_centers)';

% Get striatum depth group by across-experiment alignment
n_depths = n_aligned_depths;
depth_group = aligned_str_depth_group;

binned_spikes = zeros(n_depths,length(time_bin_centers));
for curr_depth = 1:n_depths
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
end

% Get rid of NaNs (if no data?)
binned_spikes(isnan(binned_spikes)) = 0;

kernel_t = [-0.3,0.3];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
zs = [false,true];
cvfold = 10;

%%% Regress MUA from cortex
[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(dfVdf_resample, ...
    binned_spikes,kernel_frames,lambda,zs,cvfold);

% Reshape kernel and convert to pixel space
r = reshape(k,length(use_svs),length(kernel_frames),size(binned_spikes,1));

aUdf = single(AP_align_widefield(Udf,animal,day));
r_px = zeros(size(aUdf,1),size(aUdf,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3)
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(aUdf(:,:,use_svs),r(:,:,curr_spikes));
end

depth_kernels = squeeze(r_px(:,:,kernel_frames == 0,:));




% Get probe depths per micron within the atlas

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

% Plot the probe trajectory and Allen atlas in native orientation
% (note the CCF is rotated to allow for dim 1 = x)

h = figure; hold on
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
cameratoolbar(h,'SetCoordSys','y');
cameratoolbar(h,'SetMode','orbit');

scatter3(probe_depths(:,1),probe_depths(:,2),probe_depths(:,3),5,'k');
drawnow;

% Get the Allen projection sites in regular intervals in probe striatum
% NOTE: update this later to use aligned depths instead of even segments

% Convert probe depths subscripts to indicies
probe_depths_ind = sub2ind(size(av),probe_depths(:,1),probe_depths(:,2),probe_depths(:,3));

% Get structures that the probe is in
probe_structures = av(probe_depths_ind);

% Get probe in striatum according to CCF
str_id = find(strcmp(st.safe_name,'Caudoputamen'));
probe_structures_str = probe_structures == str_id;
str_depth_ccf = [probe_depths(find(probe_structures_str,1,'first'),:); ...
    probe_depths(find(probe_structures_str,1,'last'),:)];

% OPTION 1: use query points as relative absolute depth in um
%
% [~,idx,~] = unique(spike_templates);
% template_aligned_depth = aligned_str_depth_group(idx)+1;
% template_aligned_depth(isnan(template_aligned_depth)) = 1;
% 
% for curr_depth = 1:n_aligned_depths
%     depth_top = min(spike_depths(aligned_str_depth_group == curr_depth));
%     depth_bottom = max(spike_depths(aligned_str_depth_group == curr_depth));    
%     aligned_depth_relative_ccf = ([depth_top,depth_bottom] - str_depth(2))/10;
% end
% % SANITY CHECK THIS WITH A PLOT
% 
% 
% % % Get striatum length on probe in CCF vs ephys
% % str_length_ccf = pdist2(str_depth_ccf(1,:),str_depth_ccf(2,:))*10;
% % str_length_ephys = diff(str_depth);
% 
% % Get probe vector in striatum using the striatum end
% r0 = str_depth_ccf(2,:);
% xyz = bsxfun(@minus,str_depth_ccf,r0);
% [~,~,V] = svd(xyz,0);
% ccf_probe_direction = V(:,1);
% 
% % Get section of CCF probe corresponding to ephys by distance
% probe_section_evaluate = [aligned_depth_relative_ccf(1),aligned_depth_relative_ccf(2)];
% probe_section_ccf = round(bsxfun(@plus,bsxfun(@times,probe_section_evaluate',ccf_probe_direction'),r0));
% probe_query_um = mean(probe_section_ccf,1)*10;


% OPTION 2: don't use um, instead normalize str distance from what it is to
% what it should be and divide evenly (this won't work if part of the
% striatum isn't recorded from, i.e. the top part is messed up, this seems
% unlikely on the first day though)

[~,idx,~] = unique(spike_templates);
template_aligned_depth = aligned_str_depth_group(idx)+1;
template_aligned_depth(isnan(template_aligned_depth)) = 1;

depth_str_um = nan(n_aligned_depths,3);
for curr_depth = 1:n_aligned_depths
    % (skip if there aren't units in this depth)
    if ~any(aligned_str_depth_group == curr_depth)
        continue
    end
    depth_top = min(spike_depths(aligned_str_depth_group == curr_depth));
    depth_bottom = max(spike_depths(aligned_str_depth_group == curr_depth));  
    depth_center = mean([depth_top,depth_bottom]);
    
    depth_str_relative = (depth_center-str_depth(1))/diff(str_depth);    
    depth_str_um(curr_depth,:) = round((str_depth_ccf(1,:) + diff(str_depth_ccf,[],1)*depth_str_relative)*10);
end
depth_str_um(all(isnan(depth_str_um),2),:) = [];

% Mirror all of the locations across the midline and get projections from
% both (because the cortical injections in the database aren't evenly
% distributed, so this is more comprehensive)
bregma_um = allenCCFbregma*10;
ccf_midline = bregma_um(3);
hemisphere = sign(depth_str_um(1,3) - ccf_midline);
str_depths_mirror = [depth_str_um(:,1:2),ccf_midline - hemisphere*abs(depth_str_um(:,3)-ccf_midline)];

max_sites = 50;
str_depths_query = [depth_str_um;str_depths_mirror];
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
plot_colors = repmat(cmap(round(linspace(1,size(cmap,1),size(depth_str_um,1))),:),2,1);
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

% Convert points from CCF to widefield and plot overlay
ccf_tform_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\ccf_tform'];
load(ccf_tform_fn);

um2pixel = 20.6;
injection_coordinates_wf = cellfun(@(x) ...
    [x(:,[3,1]).*(1/(um2pixel)),ones(size(x,1),1)]*ccf_tform.T, ...
    injection_coordinates_standardized,'uni',false);

figure('Name',[animal ' ' day ': widefield-aligned']); 
colormap(brewermap([],'*RdBu'));
used_depths = unique(aligned_str_depth_group(~isnan(aligned_str_depth_group)));
for curr_depth = 1:n_aligned_depths
    subplot(1,n_aligned_depths,curr_depth); hold on; axis image off;
    set(gca,'YDir','reverse');
    
    imagesc(depth_kernels(:,:,curr_depth));
    caxis([-prctile(abs(depth_kernels(:)),99),prctile(abs(depth_kernels(:)),99)]);
    
    curr_inj_idx = find(used_depths == curr_depth);
    curr_mirror_inj_idx = length(used_depths) + curr_inj_idx;
    if ~isempty(curr_inj_idx)
        % (plot points from both hemispheres)
        scatter(injection_coordinates_wf{curr_inj_idx}(:,1), ...
            injection_coordinates_wf{curr_inj_idx}(:,2), ...
            [projection_strength_normalized{curr_inj_idx}*50 + 10], ...
            'k','filled');
        scatter(injection_coordinates_wf{curr_mirror_inj_idx}(:,1), ...
            injection_coordinates_wf{curr_mirror_inj_idx}(:,2), ...
            [projection_strength_normalized{curr_mirror_inj_idx}*50 + 10], ...
            'k','filled');        
    end    
    AP_reference_outline('ccf_aligned','k');
end




