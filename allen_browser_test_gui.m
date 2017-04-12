function allen_browser_test_gui(tv,av,st,bregma)
% allen_browser_test_gui(tv,av,st)

% If not already loaded in, load in atlas
if nargin < 4
    cd('C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF')
    tv = readNPY('template_volume_10um.npy'); % grey-scale "background signal intensity"
    av = readNPY('annotation_volume_10um_by_index.npy'); % the number at each pixel labels the area, see note below
    st = loadStructureTree('structure_tree_safe.csv'); % a table of what all the labels mean
    bregma = allenCCFbregma();
end

% Set up the gui and axes
probe_atlas_gui = figure('Toolbar','none','color','w');
colormap(gray);

axes_3d = axes;
hold(axes_3d,'on');
axis vis3d
set(gca, 'ZDir', 'reverse')
axis equal
axis off
view([-30,25]);

% (coordinates in: plot [ap,ml,dv], volume [ap,dv,ml])
[ap_max,dv_max,ml_max] = size(tv);
axis manual
xlim([-10,ap_max+10])
ylim([-10,ml_max+10])
zlim([-10,dv_max+10])

% Plot outline of the brain
% (with downsampled polygons)
% slice_spacing = 10;
% target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) > 1,[3,1,2]);
% structure_patch = isosurface(target_volume,0);
% structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
% target_structure_color = [0.7,0.7,0.7];
% cortex_wire_h = patch('Vertices',structure_wire.vertices*slice_spacing, ...
%     'Faces',structure_wire.faces, ...
%     'FaceColor','none','EdgeColor',target_structure_color);
% (with Nick's method)
gridIn3D(double(av > 1),0.5,80,bregma);


% Outline structures
structure_alpha = 0.2;
structure_patch_h = gobjects;

target_structure = 'CP';
st_idx = find(strcmp(st.acronym,target_structure));
target_structure_color = hex2dec(reshape(st.color_hex_triplet{st_idx},3,[]))./255;
slice_spacing = 10;
structure_patch = isosurface(permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) == st_idx,[3,1,2]),0);
structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);

structure_patch_h(1) = patch('Vertices',structure_patch.vertices*slice_spacing, ...
    'Faces',structure_patch.faces, ...
    'FaceColor',target_structure_color,'EdgeColor','none','FaceAlpha',structure_alpha);

target_structure = 'SNr';
st_idx = find(strcmp(st.acronym,target_structure));
target_structure_color = hex2dec(reshape(st.color_hex_triplet{st_idx},3,[]))./255;
slice_spacing = 10;
structure_patch = isosurface(permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) == st_idx,[3,1,2]),0);
structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);

structure_patch_h(2) = patch('Vertices',structure_patch.vertices*slice_spacing, ...
    'Faces',structure_patch.faces, ...
    'FaceColor',target_structure_color,'EdgeColor','none','FaceAlpha',structure_alpha);

% Set up the probe
probe_ref_top = [ap_max,bregma(3),0];
probe_ref_bottom = [0,bregma(3),0];
probe_vector = [probe_ref_top',probe_ref_bottom'];
probe_ref_line = line(probe_vector(1,:),probe_vector(2,:),probe_vector(3,:),'linewidth',4,'color','k','linestyle','--');

probe_angle_ml = 10;
probe_angle_dv = 45;

rotate(probe_ref_line,[0,0,1],probe_angle_ml,probe_ref_top)
rotate(probe_ref_line,[0,1,0],probe_angle_dv,probe_ref_top)

% set up probe-angle plane
slice_plot = update_slice(surface('EdgeColor','none'),tv,av,probe_ref_line);

% Create and store guidata
gui_data = struct;

gui_data.tv = tv;
gui_data.av = av;
gui_data.st = st;
gui_data.bregma = bregma;

%gui_data.handles.cortex_wire_h = cortex_wire_h; % commented until Nick's
%function has handles too
gui_data.handles.structure_patch_h = structure_patch_h;

gui_data.handles.axes_3d = axes_3d;
gui_data.handles.slice_plot = slice_plot;
gui_data.handles.probe_ref_line = probe_ref_line;



% Set functions for mouse wheel / button press and pass necessary handles
set(probe_atlas_gui,'KeyPressFcn',@key_press); 



% Upload gui_data
guidata(probe_atlas_gui, gui_data);

end


function key_press(probe_atlas_gui,eventdata)

% Get guidata
gui_data = guidata(probe_atlas_gui);

switch eventdata.Key
    
    case 'uparrow'
        
        ap_offset = -20;
        
        set(gui_data.handles.probe_ref_line,'XData',get(gui_data.handles.probe_ref_line,'XData') + ap_offset);
        gui_data.handles.slice_plot = update_slice(gui_data.handles.slice_plot, ...
            gui_data.tv,gui_data.av,gui_data.handles.probe_ref_line);
        
    case 'downarrow'
        
        ap_offset = 20;
        
        set(gui_data.handles.probe_ref_line,'XData',get(gui_data.handles.probe_ref_line,'XData') + ap_offset);
        gui_data.handles.slice_plot = update_slice(gui_data.handles.slice_plot, ...
            gui_data.tv,gui_data.av,gui_data.handles.probe_ref_line);
        
    case 'rightarrow'
        
        ml_offset = 20;
        
        set(gui_data.handles.probe_ref_line,'YData',get(gui_data.handles.probe_ref_line,'YData') + ml_offset);
        gui_data.handles.slice_plot = update_slice(gui_data.handles.slice_plot, ...
            gui_data.tv,gui_data.av,gui_data.handles.probe_ref_line);
        
    case 'leftarrow'
        
        ml_offset = -20;
        
        set(gui_data.handles.probe_ref_line,'YData',get(gui_data.handles.probe_ref_line,'YData') + ml_offset);
        gui_data.handles.slice_plot = update_slice(gui_data.handles.slice_plot, ...
            gui_data.tv,gui_data.av,gui_data.handles.probe_ref_line);
        
    case 'm'
        % Toggle face visibility of target structure
        curr_structure_FaceAlpha = [gui_data.handles.structure_patch_h.FaceAlpha];
        if any(curr_structure_FaceAlpha)
            new_structure_FaceAlpha = 0.2;
        else
            new_structure_FaceAlpha = 0;
        end
        [gui_data.handles.structure_patch_h.FaceAlpha] = ...
            deal(new_structure_FaceAlpha);
        
    case 's'
        % Toggle slice visibility
        switch gui_data.handles.slice_plot.Visible
            case 'on'
                gui_data.handles.slice_plot.Visible = 'off';
            case 'off'
                gui_data.handles.slice_plot.Visible = 'on';
        end
        
    case 'u'
        % Update slice
        gui_data.handles.slice_plot = update_slice( ...
            gui_data.handles.slice_plot,gui_data.tv,gui_data.av,gui_data.handles.probe_ref_line);
        
    case 'r'
        % Toggle 3D rotation
        rotate3d(probe_atlas_gui)
        %(this is annoyingly complicated: need to restore key-press
        %functionality after turning on 3D rotation)       
        hManager = uigetmodemanager(probe_atlas_gui);
        [hManager.WindowListenerHandles.Enabled] = deal(false);
        set(probe_atlas_gui,'KeyPressFcn',@key_press); 
        
end

% Upload gui_data
guidata(probe_atlas_gui, gui_data);

end

function slice_plot = update_slice(slice_plot,tv,av,probe_ref_line)

% Get current position of camera
curr_campos = campos;

% Get probe vector
probe_ref_top = [probe_ref_line.XData(1),probe_ref_line.YData(1),probe_ref_line.ZData(1)];
probe_ref_bottom = [probe_ref_line.XData(2),probe_ref_line.YData(2),probe_ref_line.ZData(2)];
probe_vector = probe_ref_top - probe_ref_bottom;

% Get probe-camera vector
probe_camera_vector = probe_ref_top - curr_campos;

% Get the vector to plot the plane in (along with probe vector)
plot_vector = cross(probe_camera_vector,probe_vector);

% Get the normal vector of the plane
normal_vector = cross(plot_vector,probe_vector);

% Get the plane offset through the probe
plane_offset = -(normal_vector*probe_ref_top');

% Define a plane within the boundaries orthogonal to viewing axis on probe
[plane_x,plane_y] = meshgrid(1:5:size(tv,1),1:5:size(tv,3));

plane_z = ...
    (normal_vector(1)*plane_x+normal_vector(2)*plane_y + plane_offset)/ ...
    -normal_vector(3);

% Get the coordiates on the plane
x_idx = plane_x;
y_idx = plane_y;
z_idx = round(plane_z);

% Find plane coordinates in bounds with the volume
use_zd = z_idx > 0 & z_idx < size(tv,2);
curr_slice_idx = sub2ind(size(tv),x_idx(use_zd),z_idx(use_zd),y_idx(use_zd));

%(coordinates in: plot [ap,ml,dv], volume [ap,dv,ml])

% Find plane coordinates that contain brain
curr_slice_isbrain = false(size(use_zd));
curr_slice_isbrain(use_zd) = av(curr_slice_idx) > 1;

% Index coordinates in bounds + with brain
grab_pix_idx = sub2ind(size(tv),x_idx(curr_slice_isbrain),z_idx(curr_slice_isbrain),y_idx(curr_slice_isbrain));

% Grab pixels from volume
curr_slice = nan(size(use_zd));
curr_slice(curr_slice_isbrain) = tv(grab_pix_idx);

% Update the slice display
set(slice_plot,'XData',plane_x,'YData',plane_y,'ZData',plane_z,'CData',curr_slice);

end






