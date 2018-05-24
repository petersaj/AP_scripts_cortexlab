%% Example widefield movie

% animal = 'AP029'; day = '2017-12-12'; experiment = 1; verbose = true; AP_load_experiment;


AP_wfmovies(U,fV,frame_t)
caxis([-2500,2500])

frames = 435-455;


%% Brain outline movies

allen_atlas_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

slice_spacing = 10;

plot_structure_color = [1,0.6,0.6];
structure_alpha = 1;

use_av = av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) > 1;

structure_3d = isosurface(permute(use_av,[3,1,2]),0);
structure_3d_smoothed = smoothpatch(structure_3d,1,5);

figure('Color','w')
brain_patch = patch('Vertices',structure_3d_smoothed.vertices*slice_spacing, ...
    'Faces',structure_3d_smoothed.faces, ...
    'FaceColor',plot_structure_color,'EdgeColor','none','FaceAlpha',structure_alpha);

set(brain_patch,'DiffuseStrength',0.6,'SpecularStrength',0,'AmbientStrength',0.5);

axis vis3d image off;

set(gca,'Zdir','reverse')

% Set the properties for the mouse photo
view([-159,15]);
h = camlight('headlight');
set(h,'style','infinite');
lighting gouraud












