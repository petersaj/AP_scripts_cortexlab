% Example pipeline for processing histology

allen_atlas_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']);
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

im_path = 'C:\Users\Andrew\Desktop\test_histology\subset';
slice_path = [im_path filesep 'slices'];

% Set white balance and resize slide images, extract slice images
AP_process_histology(im_path);

% (optional) Rotate, pad, and center slice images
AP_rotate_histology(slice_path);

% Find CCF slices corresponding to each histology slice
AP_grab_histology_ccf(tv,av,st,slice_path);

% Align CCF slices and histology slices
AP_align_histology_ccf(tv,av,st,slice_path);

% Display aligned CCF over histology slices
AP_view_aligned_histology(tv,av,st,slice_path);

% Display histology within 3D CCF
AP_view_aligned_histology_volume(tv,av,st,slice_im_path);

% Get probe trajectory from histology, convert to CCF coordinates
AP_get_probe_histology(slice_path);




