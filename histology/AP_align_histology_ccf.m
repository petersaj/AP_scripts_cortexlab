function AP_align_histology_ccf(tv,av,st,slice_im_path)


% Currently: need to have the plane to draw always be normal to the camera
% axis and then the scroll wheel would need to move along that axis,
% actually going to take a lot of work to figure out I think

% at the moment a probe axis is defined straight down - instead, have the
% probe axis run through a plane perpendicular to the camera, then define
% the translation by the axis from the camera?

% can plane just be gotten from campos somehow?


% Meant to be a non-shitty version of AtlasTransformBrowser

% Initialize guidata
gui_data = struct;
gui_data.tv = tv;
gui_data.av = av;
gui_data.st = st;

% Load in slice images
slice_im_dir = dir([slice_im_path filesep '*.tif']);
slice_im_fn = natsortfiles(cellfun(@(path,fn) [path filesep fn], ...
    {slice_im_dir.folder},{slice_im_dir.name},'uni',false));
gui_data.slice_im = cell(length(slice_im_fn),1);
for curr_slice = 1:length(slice_im_fn)
    gui_data.slice_im{curr_slice} = imread(slice_im_fn{curr_slice});
end

% Set up axis for histology image
gui_fig = figure('WindowScrollWheelFcn',@scroll_atlas_slice, ...
    'KeyPressFcn',@keypress);

gui_data.histology_ax = subplot(1,2,1,'YDir','reverse'); 
hold on; axis image off;
gui_data.histology_im_h = image(gui_data.slice_im{1},'Parent',gui_data.histology_ax);
gui_data.curr_histology_slice = 1;

% Set up 3D atlas axis
gui_data.atlas_ax = subplot(1,2,2,'ZDir','reverse','color','k', ...
    'XTick',[],'YTick',[],'ZTick',[]);
hold on
axis vis3d equal manual
view([90,0]);
[ap_max,dv_max,ml_max] = size(tv);
xlim([1,ap_max]);
ylim([1,ml_max]);
zlim([1,dv_max]);
colormap(gui_data.atlas_ax,'gray');
caxis([0,400]);

% Create slice object and first slice point
gui_data.atlas_slice_plot = surface(gui_data.atlas_ax,'EdgeColor','none'); % Slice on 3D atlas
gui_data.atlas_slice_point = camtarget;

% Upload gui data
guidata(gui_fig,gui_data);

% Draw the first slice
update_atlas_slice(gui_fig);

end 



function scroll_atlas_slice(gui_fig,eventdata)

% Get guidata
gui_data = guidata(gui_fig);

% Move slice point along camera -> center axis

% Grab current camera angle
[cam_az,cam_el] = view(gui_data.atlas_ax);

% Camera azimuth is 90 degrees offset from spherical standard
cam_az_sphere = cam_az - 90;
% Camera elevation is reversed (because of CCF orientation)
cam_el_sphere = -cam_el;

[cam_vector_x,cam_vector_y,cam_vector_z] = ...
    sph2cart(deg2rad(cam_az_sphere),deg2rad(cam_el_sphere),eventdata.VerticalScrollCount);
cam_vector = [cam_vector_x,cam_vector_y,cam_vector_z];

% Move slice point
gui_data.atlas_slice_point = gui_data.atlas_slice_point + cam_vector;

% Upload gui data
guidata(gui_fig, gui_data);

% Uodate slice
update_atlas_slice(gui_fig)

end


function keypress(gui_fig,eventdata)

% Get guidata
gui_data = guidata(gui_fig);

disp(eventdata.Key);

switch eventdata.Key
    case 'leftarrow'
        set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [1,0]);
        update_atlas_slice(gui_fig)
    case 'rightarrow'
        set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [-1,0]);
        update_atlas_slice(gui_fig)
    case 'uparrow'
        set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [0,-1]);
        update_atlas_slice(gui_fig)
    case 'downarrow'
        set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [0,1]);
        update_atlas_slice(gui_fig)
    case '1'
        gui_data.curr_histology_slice = ...
            min(gui_data.curr_histology_slice + 1,length(gui_data.slice_im));
        guidata(gui_fig, gui_data);
        update_histology_slice(gui_fig);
    case '2'
        gui_data.curr_histology_slice = max(gui_data.curr_histology_slice - 1,1);
        guidata(gui_fig, gui_data);
        update_histology_slice(gui_fig);
end


end

function update_histology_slice(gui_fig)

% Get guidata
gui_data = guidata(gui_fig);

% Set next histology slice
set(gui_data.histology_im_h,'CData',gui_data.slice_im{gui_data.curr_histology_slice})

% Upload gui data
guidata(gui_fig, gui_data);

end


function update_atlas_slice(gui_fig)

% Get guidata
gui_data = guidata(gui_fig);

% Get plane normal to the camera -> center axis, grab voxels on plane

% Grab current camera angle
[cam_az,cam_el] = view(gui_data.atlas_ax);

% Camera azimuth is 90 degrees offset from spherical standard
cam_az_sphere = cam_az - 90;
% Camera elevation is reversed (because of CCF orientation)
cam_el_sphere = -cam_el;

[cam_vector_x,cam_vector_y,cam_vector_z] = sph2cart(deg2rad(cam_az_sphere),deg2rad(cam_el_sphere),1);
cam_vector = [cam_vector_x,cam_vector_y,cam_vector_z];
plane_offset = -(cam_vector*gui_data.atlas_slice_point');

% Define a plane of points to index
% (the plane grid is defined based on the which cardinal plan is most
% orthogonal to the plotted plane. this is janky but it works)
slice_px_space = 3;

[~,cam_plane] = max(abs(cam_vector./norm(cam_vector)));

switch cam_plane
    
    case 1
        [plane_y,plane_z] = meshgrid(1:slice_px_space:size(gui_data.tv,3),1:slice_px_space:size(gui_data.tv,2));
        plane_x = ...
            (cam_vector(2)*plane_y+cam_vector(3)*plane_z + plane_offset)/ ...
            -cam_vector(1);
        
    case 2
        [plane_x,plane_z] = meshgrid(1:slice_px_space:size(gui_data.tv,1),1:slice_px_space:size(gui_data.tv,2));
        plane_y = ...
            (cam_vector(1)*plane_x+cam_vector(3)*plane_z + plane_offset)/ ...
            -cam_vector(2);
        
    case 3
        [plane_x,plane_y] = meshgrid(1:slice_px_space:size(gui_data.tv,1),1:slice_px_space:size(gui_data.tv,3));
        plane_z = ...
            (cam_vector(1)*plane_x+cam_vector(2)*plane_y + plane_offset)/ ...
            -cam_vector(3);
        
end

% Get the coordiates on the plane
x_idx = round(plane_x);
y_idx = round(plane_y);
z_idx = round(plane_z);

% Find plane coordinates in bounds with the volume
use_xd = x_idx > 0 & x_idx < size(gui_data.tv,1);
use_yd = y_idx > 0 & y_idx < size(gui_data.tv,3);
use_zd = z_idx > 0 & z_idx < size(gui_data.tv,2);
use_idx = use_xd & use_yd & use_zd;

curr_slice_idx = sub2ind(size(gui_data.tv),x_idx(use_idx),z_idx(use_idx),y_idx(use_idx));

% Find plane coordinates that contain brain
curr_slice_isbrain = false(size(use_idx));
curr_slice_isbrain(use_idx) = gui_data.av(curr_slice_idx) > 0;

% Index coordinates in bounds + with brain
grab_pix_idx = sub2ind(size(gui_data.tv),x_idx(curr_slice_isbrain),z_idx(curr_slice_isbrain),y_idx(curr_slice_isbrain));

% Grab pixels from (selected) volume
curr_slice = nan(size(use_idx));
curr_slice(curr_slice_isbrain) = gui_data.tv(grab_pix_idx);

% Update the slice display
set(gui_data.atlas_slice_plot,'XData',plane_x,'YData',plane_y,'ZData',plane_z,'CData',curr_slice);

% Upload gui_data
guidata(gui_fig, gui_data);

end



%% Storing here for now - transform
function asdf

example_slice_fn = 'C:\Users\Andrew\Desktop\test_histology\processed\ps_slices\12.tif';
example_slice = imread(example_slice_fn);
figure;imshow(example_slice);


allen_atlas_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
ccf_fig = allenAtlasBrowser(tv,av,st);


% tv_slice/av_slice = get ccf fig
[moving_points,fixed_points] = cpselect(example_slice,rescale(tv_slice),'Wait',true);



% this gui isn't good and need to have something that updates on every
% click, because otherwise no idea how it'll turn out

tform = fitgeotrans(moving_points,fixed_points,'affine');
R = imref2d(size(tv_slice));
example_slice_warp = imwarp(example_slice, tform, 'OutputView',R);


tform = fitgeotrans(fixed_points,moving_points,'affine');
R = imref2d(size(example_slice));
av_slice_warp = imwarp(av_slice, tform,'nearest','OutputView',R);


figure;

% slice -> ccf
slice_boundaries = bwboundaries(round(conv2(av_slice,ones(3)./9,'same')) ~= av_slice,4);

subplot(1,2,1);
imshow(example_slice_warp);

for i = 1:length(slice_boundaries)
   line(slice_boundaries{i}(:,2),slice_boundaries{i}(:,1), ...
       'color',[0.5,0.5,0.5],'linewidth',1); 
end
title('Slice -> CCF');

% ccf -> slice
slice_boundaries = bwboundaries(round(conv2(av_slice_warp,ones(3)./9,'same')) ~= av_slice_warp,4);

subplot(1,2,2);
imshow(example_slice);

for i = 1:length(slice_boundaries)
   line(slice_boundaries{i}(:,2),slice_boundaries{i}(:,1), ...
       'color',[0.5,0.5,0.5],'linewidth',1); 
end
title('CCF -> Slice');

% caxis([1,1305]);
% colormap(allen_ccf_colormap)
% axis image off





end




















