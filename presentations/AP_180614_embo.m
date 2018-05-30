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

%% Outline just the cortex

allen_atlas_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

% Get first brain pixel from top-down, get annotation at that point
[~,top_down_depth] = max(av>1, [], 2);
top_down_depth = squeeze(top_down_depth);

[xx,yy] = meshgrid(1:size(top_down_depth,2), 1:size(top_down_depth,1));
top_down_annotation = reshape(av(sub2ind(size(av),yy(:),top_down_depth(:),xx(:))), size(av,1), size(av,3));

% Get all labelled areas
used_areas = unique(top_down_annotation(:));

% Restrict to only cortical areas
structure_id_path = cellfun(@(x) textscan(x(2:end),'%d', 'delimiter',{'/'}),st.structure_id_path);

ctx_path = [997,8,567,688,695,315];
ctx_idx = find(cellfun(@(id) length(id) > length(ctx_path) & ...
    all(id(min(length(id),length(ctx_path))) == ctx_path(min(length(id),length(ctx_path)))),structure_id_path));

ctx_av = ismember(av,ctx_idx);

top_down_ctx = permute(any(ctx_av,2),[1,3,2]);

ctx_boundaries = bwboundaries(top_down_ctx);

figure;
plot(ctx_boundaries{1}(:,2),ctx_boundaries{1}(:,1),'k','linewidth',2);
set(gca,'YDir','reverse');
axis off;


%% Example MUA

animal = 'AP028'; day = '2017-12-16'; experiment = 1; 
load_parts.ephys = true; verbose = true; AP_load_experiment;

% Set upsample value for regression
upsample_factor = 2;
sample_rate = (1/median(diff(frame_t)))*upsample_factor;

% Skip the first/last n seconds to do this
skip_seconds = 60;
time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% (to use aligned striatum depths)
n_depths = n_aligned_depths;
depth_group = aligned_str_depth_group;

binned_spikes = zeros(n_depths,length(time_bins)-1);
for curr_depth = 1:n_depths
    
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
    
end

smooth_size = 8;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
binned_spikes_smooth = conv2(binned_spikes,smWin,'same');

plot_t = time_bin_centers > 690 & time_bin_centers < 705;

figure;
AP_stackplot(binned_spikes_smooth(:,plot_t)',time_bin_centers(plot_t),8,true,'k');


%% Plotting WF areas colored by weight from ctx/str regression

wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

roi_cat = cat(3,wf_roi.mask);

use_t_weights = t > -0.1 & t < 0;

interarea_kernel_mean = nanmean(cat(3,interarea_kernel_all{:}),3);
roi_weights = nanmean(interarea_kernel_mean(:,use_t_weights),2);

figure; hold on
set(gca,'YDir','reverse');
AP_reference_outline('ccf_aligned','k');
for curr_roi = 1:n_rois
    curr_roi_boundary = cell2mat(bwboundaries(roi_cat(:,:,curr_roi)));
    patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),roi_weights(curr_roi));
end
axis image off;
colormap(colormap_BlueWhiteRed);
caxis([-max(abs(caxis)),max(abs(caxis))]);







