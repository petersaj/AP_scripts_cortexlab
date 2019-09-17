function AP_view_aligned_histology_volume(tv,av,st,slice_im_path)
% AP_view_aligned_histology_volume(tv,av,st,slice_im_path)
%
% Plot histology warped onto CCF volume
% Andy Peters (peters.andrew.j@gmail.com)

% Initialize guidata
gui_data = struct;
gui_data.tv = tv;
gui_data.av = av;
gui_data.st = st;

% Load in slice images
gui_data.slice_im_path = slice_im_path;
slice_im_dir = dir([slice_im_path filesep '*.tif']);
slice_im_fn = natsortfiles(cellfun(@(path,fn) [path filesep fn], ...
    {slice_im_dir.folder},{slice_im_dir.name},'uni',false));
gui_data.slice_im = cell(length(slice_im_fn),1);
for curr_slice = 1:length(slice_im_fn)
    gui_data.slice_im{curr_slice} = imread(slice_im_fn{curr_slice});
end

% Load corresponding CCF slices
ccf_slice_fn = [slice_im_path filesep 'histology_ccf.mat'];
load(ccf_slice_fn);
gui_data.histology_ccf = histology_ccf;

% Load histology/CCF alignment
% Manual control point (control points)
ccf_alignment_fn = [slice_im_path filesep 'histology_ccf_alignment.mat'];
if exist(ccf_alignment_fn,'file')
    load(ccf_alignment_fn);
    gui_data.histology_ccf_alignment = histology_ccf_alignment;
end
% Automated outline (affine transform matrix)
auto_ccf_alignment_fn = [slice_im_path filesep 'atlas2histology_tform.mat'];
if exist(auto_ccf_alignment_fn,'file')
    load(auto_ccf_alignment_fn);
    gui_data.histology_ccf_auto_alignment = atlas2histology_tform;
end

% Warp histology to CCF
gui_data.atlas_aligned_histology = cell(length(gui_data.slice_im),1);
for curr_slice = 1:length(gui_data.slice_im)
    curr_av_slice = gui_data.histology_ccf(curr_slice).av_slices;
    curr_av_slice(isnan(curr_av_slice)) = 1;
    curr_slice_im = gui_data.slice_im{curr_slice};
    
    % Manual control points
    if isfield(gui_data,'histology_ccf_alignment')
        tform = fitgeotrans(gui_data.histology_ccf_alignment.histology_control_points{curr_slice}, ...
            gui_data.histology_ccf_alignment.atlas_control_points{curr_slice},'affine');
    end
    
    % Automated outline
    if isfield(gui_data,'histology_ccf_auto_alignment')
        tform = affine2d;
        tform.T = gui_data.histology_ccf_auto_alignment{curr_slice};
        % (transform is CCF -> histology, invert for other direction)
        tform = invert(tform);
    end

    tform_size = imref2d([size(gui_data.histology_ccf(curr_slice).av_slices,1), ...
        size(gui_data.histology_ccf(curr_slice).av_slices,2)]);
    
    gui_data.atlas_aligned_histology{curr_slice} = ...
        imwarp(curr_slice_im,tform,'nearest','OutputView',tform_size);
    
end

% Create figure
gui_fig = figure;

% Set up 3D plot for volume viewing
axes_atlas = axes;
[~, brain_outline] = plotBrainGrid([],axes_atlas);
set(axes_atlas,'YDir','reverse','ZDir','reverse');
hold(axes_atlas,'on');
axis vis3d equal off manual
view([-30,25]);
caxis([0 300]);
[ap_max,dv_max,ml_max] = size(tv);
xlim([-10,ap_max+10])
ylim([-10,ml_max+10])
zlim([-10,dv_max+10])
colormap(hot);

% Turn on rotation by default
h = rotate3d(axes_atlas);
h.Enable = 'on';

% Draw all aligned slices
plot_channel = 1; % Probes are usually in red
histology_surf = gobjects(length(gui_data.slice_im));
for curr_slice = 1:length(gui_data.slice_im)     
        histology_surf(curr_slice) = surface( ...
            gui_data.histology_ccf(curr_slice).plane_ap, ...
            gui_data.histology_ccf(curr_slice).plane_ml, ...
            gui_data.histology_ccf(curr_slice).plane_dv);
        histology_surf(curr_slice).FaceColor = 'texturemap';
        histology_surf(curr_slice).FaceAlpha = 'texturemap';
        histology_surf(curr_slice).EdgeColor = 'none';
        histology_surf(curr_slice).CData = gui_data.atlas_aligned_histology{curr_slice}(:,:,plot_channel);
        histology_surf(curr_slice).AlphaData = gui_data.atlas_aligned_histology{curr_slice}(:,:,plot_channel);
end

























