% Example pipeline for processing histology

% Load CCF atlas
allen_atlas_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']);
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

% Set paths for histology images and directory to save slice/alignment
im_path = 'C:\Users\Andrew\Desktop\temp_histology';
slice_path = [im_path filesep 'slices'];

% Set white balance and resize slide images, extract slice images
AP_process_histology(im_path);

% (optional) Rotate, pad, and center slice images
AP_rotate_histology(slice_path);

% Find CCF slices corresponding to each histology slice
AP_grab_histology_ccf(tv,av,st,slice_path);

% Align CCF slices and histology slices
% (manually, control points)
AP_align_histology_ccf(tv,av,st,slice_path);
% (automatically, by outline)
AP_auto_align_histology_ccf(slice_path)

% Display aligned CCF over histology slices
AP_view_aligned_histology(st,slice_path);

% Display histology within 3D CCF
AP_view_aligned_histology_volume(tv,av,st,slice_path);

% Get probe trajectory from histology, convert to CCF coordinates
AP_get_probe_histology(tv,av,st,slice_path);

% Extract slices from full-resolution images
% (not worth it at the moment, each slice is 200 MB)
% AP_grab_fullsize_histology_slices(im_path)









