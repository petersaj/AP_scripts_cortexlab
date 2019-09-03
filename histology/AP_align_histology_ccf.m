function AP_align_histology_ccf(tv,av,st,slice_im_path)

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

% Make figure with axes for 1) image, 2) atlas, 3) overlay
gui_fig = figure('WindowScrollWheelFcn',@scroll_slice, ...
    'KeyPressFcn',@keypress);

gui_data.slice_ax = subplot(1,3,1,'YDir','reverse'); 
hold on; axis image off; title('Histology');
gui_data.slice_im_h = image(gui_data.slice_im{1},'Parent',gui_data.slice_ax); 

% gui_data.atlas_ax = subplot(1,3,2,'YDir','reverse'); 
% hold on; axis image off; colormap(gui_data.atlas_ax,gray); title('Atlas');
% caxis([0,max(tv(:))]);
% gui_data.curr_atlas_slice = 500;
% gui_data.atlas_im_h = imagesc(permute(gui_data.tv(gui_data.curr_atlas_slice,:,:),[2,3,1]),'Parent',gui_data.atlas_ax); 

gui_data.slice_atlas_overlay_ax = subplot(1,3,3,'YDir','reverse'); 
hold on; axis image off; title('Aligned atlas over histology');
gui_data.slice_atlas_overlay_im_h = image(gui_data.slice_im{1},'Parent',gui_data.slice_atlas_overlay_ax); 








% Testing 3d atlas?
bregma = allenCCFbregma;
probe_ref_top = [bregma(1),bregma(3),0];

gui_data.atlas_vector = ...
    [bregma(1),bregma(3),0; ...
    bregma(1),bregma(3),size(gui_data.tv,2)];

gui_data.atlas_ax = subplot(1,3,2,'ZDir','reverse');
hold on
axis vis3d equal off manual
view([90,0]);
[ap_max,dv_max,ml_max] = size(tv);
xlim([-10,ap_max+10])
ylim([-10,ml_max+10])
zlim([-10,dv_max+10])

gui_data.slice_plot = surface(gui_data.atlas_ax,'EdgeColor','none'); % Slice on 3D atlas



% Upload gui data
guidata(gui_fig, gui_data);


update_slice_3d(gui_fig);



end 



function scroll_slice(gui_fig,eventdata)

% Get guidata
gui_data = guidata(gui_fig);

gui_data.atlas_vector = gui_data.atlas_vector + ...
    [eventdata.VerticalScrollCount,0,0; ...
    eventdata.VerticalScrollCount,0,0];

% Upload gui data
guidata(gui_fig, gui_data);

update_slice_3d(gui_fig)

end


function keypress(gui_fig,eventdata)

% Get guidata
gui_data = guidata(gui_fig);

disp(eventdata.Key);

switch eventdata.Key
    case 'leftarrow'
        set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [1,0]);
        update_slice_3d(gui_fig)
    case 'rightarrow'
        set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [-1,0]);
        update_slice_3d(gui_fig)
    case 'rightarrow'
        set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [-1,0]);
        update_slice_3d(gui_fig)
    case 'rightarrow'
        set(gui_data.atlas_ax,'View',get(gui_data.atlas_ax,'View') + [-1,0]);
        update_slice_3d(gui_fig)
end


end


function update_slice_3d(gui_fig,varargin)

% attempting modifying allen_ccf_npx

% Get guidata
gui_data = guidata(gui_fig);

% Get current position of camera
curr_campos = campos;

% Get constant atlas vector
atlas_vector = diff(gui_data.atlas_vector,[],1);

% Get probe-camera vector
probe_camera_vector = gui_data.atlas_vector(1,:) - curr_campos;

% Get the vector to plot the plane in (along with probe vector)
plot_vector = cross(probe_camera_vector,atlas_vector);

% Get the normal vector of the plane
normal_vector = cross(plot_vector,atlas_vector);

% Get the plane offset through the probe
plane_offset = -(normal_vector*gui_data.atlas_vector(1,:)');

% Define a plane of points to index
% (the plane grid is defined based on the which cardinal plan is most
% orthogonal to the plotted plane. this is janky but it works)
slice_px_space = 3;
%[~,cam_plane] = max(abs((campos - camtarget)./norm(campos - camtarget)));

[~,cam_plane] = max(abs(normal_vector./norm(normal_vector)));

switch cam_plane
    
    case 1
        [plane_y,plane_z] = meshgrid(1:slice_px_space:size(gui_data.tv,3),1:slice_px_space:size(gui_data.tv,2));
        plane_x = ...
            (normal_vector(2)*plane_y+normal_vector(3)*plane_z + plane_offset)/ ...
            -normal_vector(1);
        
    case 2
        [plane_x,plane_z] = meshgrid(1:slice_px_space:size(gui_data.tv,1),1:slice_px_space:size(gui_data.tv,2));
        plane_y = ...
            (normal_vector(1)*plane_x+normal_vector(3)*plane_z + plane_offset)/ ...
            -normal_vector(2);
        
    case 3
        [plane_x,plane_y] = meshgrid(1:slice_px_space:size(gui_data.tv,1),1:slice_px_space:size(gui_data.tv,3));
        plane_z = ...
            (normal_vector(1)*plane_x+normal_vector(2)*plane_y + plane_offset)/ ...
            -normal_vector(3);
        
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

% Index coordinates in bounds + with brain
grab_pix_idx = sub2ind(size(gui_data.tv),x_idx,z_idx,y_idx);

% Grab pixels from (selected) volume
curr_slice = gui_data.tv(grab_pix_idx);
colormap(gui_data.atlas_ax,'gray');
caxis([0,255]);



% Update the slice display
set(gui_data.slice_plot,'XData',plane_x,'YData',plane_y,'ZData',plane_z,'CData',curr_slice);

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




















