% Example pipeline for processing histology

allen_atlas_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']);
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

im_path = 'C:\Users\Andrew\Desktop\test_histology\subset';
slice_path = [im_path filesep 'slices'];

AP_process_histology(im_path);
AP_rotate_histology(slice_path);
AP_grab_histology_ccf(tv,av,st,slice_path);
AP_align_histology_ccf(tv,av,st,slice_path);
AP_view_aligned_histology(tv,av,st,slice_path);
AP_view_aligned_histology_volume(tv,av,st,slice_im_path);
AP_get_probe_histology(slice_path);

