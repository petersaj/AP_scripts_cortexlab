%% Load CCF and set paths for slide and slice images

% Load CCF atlas
allen_atlas_path = fileparts(which('template_volume_10um.npy'));
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']);
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

%% Get animal slice path

% Set paths for histology images and directory to save slice/alignment

% % (for writing/debugging)
% im_path = 'C:\Users\Andrew\Desktop\temp_histology\AP077';
% slice_path = [im_path filesep 'slices'];

% New histology
im_path = '\\zinu.cortexlab.net\Subjects\AP118\histology';
slice_path = [im_path filesep 'slices'];

animal = 'AP118';
[probe_ccf_fn,probe_ccf_fn_exists] = AP_cortexlab_filename(animal,[],[],'probe_ccf');
if probe_ccf_fn_exists
    load(probe_ccf_fn);
else
    error('No probe ccf');
end
slice_path = fileparts(probe_ccf_fn);


%% Preprocess slide images to produce slice images

% Set white balance and resize slide images, extract slice images
% (Note: this resizes the images purely for file size reasons - the CCF can
% be aligned to histology no matter what the scaling. If pixel size is
% available in metadata then automatically scales to CCF resolution,
% otherwise user can specify the resize factor as a second argument)

% Set resize factor
% resize_factor = []; % (slides ome.tiff: auto-resize ~CCF size 10um/px)
resize_factor = 1; % (slides tiff: resize factor)

% Set slide or slice images
slice_images = false; % (images are slides - extract individual slices)
% slice_images = true; % (images are already individual slices)

% Preprocess images
AP_process_histology(im_path,resize_factor,slice_images);

% (optional) Rotate, center, pad, flip slice images
AP_rotate_histology(slice_path);

%% Align CCF to slices

% Find CCF slices corresponding to each histology slice
AP_grab_histology_ccf(tv,av,st,slice_path);

% Align CCF slices and histology slices
% (first: automatically, by outline)
AP_auto_align_histology_ccf(slice_path);
% (second: curate manually)
AP_manual_align_histology_ccf(tv,av,st,slice_path);


%% Utilize aligned CCF

% Display aligned CCF over histology slices
AP_view_aligned_histology(st,slice_path);

% Display histology within 3D CCF
AP_view_aligned_histology_volume(tv,av,st,slice_path,1);

% Get probe trajectory from histology, convert to CCF coordinates
AP_get_probe_histology(tv,av,st,slice_path);

% Align histology to electrophysiology
use_probe = 1;
AP_align_probe_histology(st,slice_path, ...
    spike_times,spike_templates,template_depths, ...
    lfp,lfp_channel_positions, ...
    use_probe);

% Extract slices from full-resolution images
% (not worth it at the moment, each slice is 200 MB)
% AP_grab_fullsize_histology_slices(im_path)

% Convert points in histology images to CCF coordinates
ccf_points = AP_histology2ccf(histology_points,slice_path);
% Concatenate points and round to nearest integer coordinate
ccf_points_cat = round(cell2mat(ccf_points));
% Get indicies from subscripts
ccf_points_idx = sub2ind(size(av),ccf_points_cat(:,1),ccf_points_cat(:,2),ccf_points_cat(:,3));
% Find annotated volume (AV) values at points
ccf_points_av = av(ccf_points_idx);
% Get areas from the structure tree (ST) at given AV values
ccf_points_areas = st(ccf_points_areas,:).safe_name;


%% Plot probe

figure; 

% Set up 3D axes
ccf_3d_axes = subplot(1,4,1);
[~, brain_outline] = plotBrainGrid([],ccf_3d_axes);
set(ccf_3d_axes,'ZDir','reverse');
hold(ccf_3d_axes,'on');
axis vis3d equal off manual
view([-30,25]);
caxis([0 300]);
[ap_max,dv_max,ml_max] = size(av);
xlim([-10,ap_max+10])
ylim([-10,ml_max+10])
zlim([-10,dv_max+10])
h = rotate3d(ccf_3d_axes);
h.Enable = 'on';

% Set up 2D axes
ccf_axes = gobjects(3,1);
ccf_axes(1) = subplot(1,4,2,'YDir','reverse');
hold on; axis image off;
ccf_axes(2) = subplot(1,4,3,'YDir','reverse');
hold on; axis image off;
ccf_axes(3) = subplot(1,4,4,'YDir','reverse');
hold on; axis image off;
for curr_view = 1:3
    curr_outline = bwboundaries(squeeze((max(av,[],curr_view)) > 1));
    cellfun(@(x) plot(ccf_axes(curr_view),x(:,2),x(:,1),'k','linewidth',2),curr_outline)
end
linkaxes(ccf_axes);

% Draw probes on all axes
for curr_probe = 1:length(probe_ccf)
    % Plot points and line of best fit
    r0 = mean(probe_ccf(curr_probe).points,1);
    xyz = bsxfun(@minus,probe_ccf(curr_probe).points,r0);
    [~,~,V] = svd(xyz,0);
    histology_probe_direction = V(:,1);
    % (make sure the direction goes down in DV - flip if it's going up)
    if histology_probe_direction(2) < 0
        histology_probe_direction = -histology_probe_direction;
    end   
    
     % Evaluate line of best fit (length of probe to deepest point)
    [~,deepest_probe_idx] = max(probe_ccf.points(:,2));
    probe_deepest_point = probe_ccf.points(deepest_probe_idx,:);
    probe_deepest_point_com_dist = pdist2(r0,probe_deepest_point);
    probe_length_ccf = 3840/10; % mm / ccf voxel size
    
    probe_line_eval = probe_deepest_point_com_dist - [probe_length_ccf,0];
    probe_line = (probe_line_eval'.*histology_probe_direction') + r0;
    
    % Draw probe in 3D view
    line(ccf_3d_axes,probe_line(:,1),probe_line(:,3),probe_line(:,2), ...
        'linewidth',2);
    plot3(ccf_3d_axes,probe_ccf(curr_probe).points(:,1), ...
        probe_ccf(curr_probe).points(:,3), ...
        probe_ccf(curr_probe).points(:,2), ...
        '.','MarkerSize',20);

end

line(ccf_axes(1),probe_line(:,3),probe_line(:,2),'linewidth',2);
line(ccf_axes(2),probe_line(:,3),probe_line(:,1),'linewidth',2);
line(ccf_axes(3),probe_line(:,2),probe_line(:,1),'linewidth',2);
















