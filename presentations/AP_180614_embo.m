%% Example widefield movie

% animal = 'AP029'; day = '2017-12-12'; experiment = 1; verbose = true; AP_load_experiment;


AP_wfmovies(U,fV,frame_t)
caxis([-2500,2500])

frames = 435-455;


%% Brain outline movies

embo_folder = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\presentations\180614_embo';

allen_atlas_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

slice_spacing = 10;
structure_alpha = 1;

brain_fig = figure('Color','w','Position', [398,84,1157,898]);

brain_color = [1,0.6,0.6];
brain_av = av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) > 1;
brain_3d = isosurface(permute(brain_av,[3,1,2]),0);
brain_3d_smoothed = smoothpatch(brain_3d,1,10);
brain_patch = patch('Vertices',brain_3d_smoothed.vertices*slice_spacing, ...
    'Faces',brain_3d_smoothed.faces, ...
    'FaceColor',brain_color,'EdgeColor','none','FaceAlpha',structure_alpha, ...
    'DiffuseStrength',0.6,'SpecularStrength',0,'AmbientStrength',0.5);

striatum_id = 574;
striatum_color = [0,0,0.8];
striatum_av = av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) == striatum_id;
striatum_3d = isosurface(permute(striatum_av,[3,1,2]),0);
striatum_3d_smoothed = smoothpatch(striatum_3d,1,10);
striatum_patch = patch('Vertices',striatum_3d_smoothed.vertices*slice_spacing, ...
    'Faces',striatum_3d_smoothed.faces, ...
    'FaceColor',striatum_color,'EdgeColor','none','FaceAlpha',structure_alpha, ...
    'DiffuseStrength',0.6,'SpecularStrength',0,'AmbientStrength',0.5);

axis vis3d image off;
xlim([-100,1480])

set(gca,'Zdir','reverse')

% Set the properties for the mouse photo
view([-180,0]);
h = camlight('headlight');
set(h,'style','infinite');
lighting gouraud

camlight(h,'headlight');

% Rotate brain side to top
cam_steps = [linspace(180,90,30);linspace(0,90,30)];
movie_frames = struct('cdata',[],'colormap',[]);
for i = 1:length(cam_steps)
    view([cam_steps(1,i),cam_steps(2,i)]);
    camlight(h,'headlight');
    movie_frames(i) = getframe(brain_fig);
end
movie_file = [embo_folder filesep 'brain_rotate_side-to-top.avi'];
writerObj = VideoWriter(movie_file);
writerObj.FrameRate = length(cam_steps);
open(writerObj);
writeVideo(writerObj,movie_frames);
close(writerObj);

% Rotate brain top to front
cam_steps = [linspace(90,-90,30);linspace(90,0,30)];
movie_frames = struct('cdata',[],'colormap',[]);
for i = 1:length(cam_steps)
    view([cam_steps(1,i),cam_steps(2,i)]);
    camlight(h,'headlight');
    movie_frames(i) = getframe(brain_fig);
end
movie_file = [embo_folder filesep 'brain_rotate_top-to-front.avi'];
writerObj = VideoWriter(movie_file);
writerObj.FrameRate = length(cam_steps);
open(writerObj);
writeVideo(writerObj,movie_frames);
close(writerObj);

% Fade brain out for striatum from front
alpha_steps = linspace(1,0,20);
view([-90,0]);
camlight(h,'headlight');
movie_frames = struct('cdata',[],'colormap',[]);
for i = 1:length(alpha_steps)
    set(brain_patch,'FaceAlpha',alpha_steps(i))
    movie_frames(i) = getframe(brain_fig);
end
movie_file = [embo_folder filesep 'brain_striatum_fade-out_front.avi'];
writerObj = VideoWriter(movie_file);
writerObj.FrameRate = length(cam_steps);
open(writerObj);
writeVideo(writerObj,movie_frames);
close(writerObj);

% Fade brain out for striatum from top
alpha_steps = linspace(1,0,20);
view([90,90]);
camlight(h,'headlight');
movie_frames = struct('cdata',[],'colormap',[]);
for i = 1:length(alpha_steps)
    set(brain_patch,'FaceAlpha',alpha_steps(i))
    movie_frames(i) = getframe(brain_fig);
end
movie_file = [embo_folder filesep 'brain_striatum_fade-out_top.avi'];
writerObj = VideoWriter(movie_file);
writerObj.FrameRate = length(cam_steps);
open(writerObj);
writeVideo(writerObj,movie_frames);
close(writerObj);






















