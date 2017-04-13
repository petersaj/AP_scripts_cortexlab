function allen_browser_test_gui(tv,av,st,bregma)
% allen_browser_test_gui(tv,av,st,bregma)
%
% This gui is for looking at trajectories in the brain with the Allen CCF
% 
% Coordinates in: plot [ap,ml,dv], volume [ap,dv,ml]
%
% TO ADD: 
% - mouse over structure names
% - live manupulator update
% - bregma-lambda scaling and angle adjustment
% - choose structures to show in mesh

% Initialize gui_data structure
gui_data = struct;

% If not already loaded in, load in atlas
if nargin < 4
    cd('C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF')
    tv = readNPY('template_volume_10um.npy'); % grey-scale "background signal intensity"
    av = readNPY('annotation_volume_10um_by_index.npy'); % the number at each pixel labels the area, see note below
    st = loadStructureTree('structure_tree_safe.csv'); % a table of what all the labels mean
    bregma = allenCCFbregma();
end

% Set up the gui and axes
probe_atlas_gui = figure('Toolbar','none','Menubar','none','color','w', ...
    'Name','Atlas-probe viewer');
colormap(gray);

axes_3d = axes;
hold(axes_3d,'on');
axis vis3d
set(gca, 'ZDir', 'reverse')
axis equal
axis off
view([-30,25]);
caxis([0 300]);

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
brain_outline = gridIn3D(double(av > 1),0.5,80,bregma);

% Set up the probe
probe_ref_top = [bregma(1),bregma(3),0];
probe_ref_bottom = [bregma(1),bregma(3),size(tv,2)];
probe_vector = [probe_ref_top',probe_ref_bottom'];
probe_ref_line = line(probe_vector(1,:),probe_vector(2,:),probe_vector(3,:),'linewidth',4,'color','k','linestyle','--');

% Store data
gui_data.tv = tv; % Intensity atlas
gui_data.av = av; % Annotated atlas
gui_data.st = st; % Labels table
gui_data.bregma = bregma; % Bregma for external referencing
gui_data.structure_plot_idx = []; % Plotted structures

%Store handles
gui_data.handles.cortex_outline = brain_outline; 
gui_data.handles.structure_patch = []; % Plotted structures
gui_data.handles.axes_3d = axes_3d; % Axes with 3D atlas
gui_data.handles.slice_plot = surface('EdgeColor','none'); % Slice on 3D atlas
gui_data.handles.probe_ref_line = probe_ref_line; % Probe reference line on 3D atlas

% Set functions for key presses
set(probe_atlas_gui,'KeyPressFcn',@key_press); 

% Upload gui_data
guidata(probe_atlas_gui, gui_data);

% Display the first slice
update_slice(probe_atlas_gui);

end

function key_press(probe_atlas_gui,eventdata)

% Get guidata
gui_data = guidata(probe_atlas_gui);

switch eventdata.Key
    
    case 'uparrow'
        
        ap_offset = -20;
        
        set(gui_data.handles.probe_ref_line,'XData',get(gui_data.handles.probe_ref_line,'XData') + ap_offset);
        update_slice(probe_atlas_gui);
        
    case 'downarrow'
        
        ap_offset = 20;
        
        set(gui_data.handles.probe_ref_line,'XData',get(gui_data.handles.probe_ref_line,'XData') + ap_offset);
        update_slice(probe_atlas_gui);
        
    case 'rightarrow'
        
        ml_offset = 20;
        
        set(gui_data.handles.probe_ref_line,'YData',get(gui_data.handles.probe_ref_line,'YData') + ml_offset);
        update_slice(probe_atlas_gui);
        
    case 'leftarrow'
        
        ml_offset = -20;
        
        set(gui_data.handles.probe_ref_line,'YData',get(gui_data.handles.probe_ref_line,'YData') + ml_offset);
        update_slice(probe_atlas_gui);
        
    case 'b' 
        % Toggle brain outline visibility
        current_visibility = gui_data.handles.cortex_outline{1}{1}{1}.Visible;
        switch current_visibility; case 'on'; new_visibility = 'off'; case 'off'; new_visibility = 'on'; end;        
        % Pull out all handles (stored as nested cells)
        brain_h_direction = cellfun(@(x) [x{:}],[gui_data.handles.cortex_outline{:}],'uni',false);
        set(horzcat(brain_h_direction{:}),'Visible',new_visibility);
        
    case 'm'
        % Toggle plotted structure visibility
        if ~isempty(gui_data.structure_plot_idx)
            current_visibility = get(gui_data.handles.structure_patch(1),'Visible');
            switch current_visibility; case 'on'; new_visibility = 'off'; case 'off'; new_visibility = 'on'; end;
            set(gui_data.handles.structure_patch,'Visible',new_visibility);
        end
        
    case 's'
        % Toggle slice visibility
        current_visibility = gui_data.handles.slice_plot(1).Visible;
        switch current_visibility; case 'on'; new_visibility = 'off'; case 'off'; new_visibility = 'on'; end;
        set(gui_data.handles.slice_plot,'Visible',new_visibility);
                
    case 'r'
        % Toggle 3D rotation
        h = rotate3d(gui_data.handles.axes_3d);
        switch h.Enable
            case 'off'
                h.Enable = 'on';
                % Update the slice whenever a rotation is completed  
                h.ActionPostCallback = @update_slice;
                %(need to restore key-press functionality with rotation)
                hManager = uigetmodemanager(probe_atlas_gui);
                [hManager.WindowListenerHandles.Enabled] = deal(false);
                set(probe_atlas_gui,'KeyPressFcn',@key_press);                                           
            case 'on'
                h.Enable = 'off';
        end
        
    case 'a' 
        % Set probe angle
        set_probe_angle(probe_atlas_gui);
        
    case 'add'
        % Add structure(s) to display
        slice_spacing = 10;

        % Prompt for which structures to show
        plot_structures = listdlg('PromptString','Select a structure to add:','ListString',gui_data.st.safe_name);
        
        for curr_plot_structure = plot_structures
            % If this label isn't used, don't plot
            if ~any(reshape(gui_data.av( ...
                    1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end),[],1) == curr_plot_structure)
                disp(['"' gui_data.st.safe_name{curr_plot_structure} '" is not parsed in the atlas'])
                continue
            end
                       
            gui_data.structure_plot_idx(end+1) = curr_plot_structure;
            
            plot_structure_color = hex2dec(reshape(gui_data.st.color_hex_triplet{curr_plot_structure},3,[]))./255;
            structure_3d = isosurface(permute(gui_data.av(1:slice_spacing:end, ...
                1:slice_spacing:end,1:slice_spacing:end) == curr_plot_structure,[3,1,2]),0);
            
            structure_alpha = 0.2;
            gui_data.handles.structure_patch(end+1) = patch('Vertices',structure_3d.vertices*slice_spacing, ...
                'Faces',structure_3d.faces, ...
                'FaceColor',plot_structure_color,'EdgeColor','none','FaceAlpha',structure_alpha);
        end
        
    case 'subtract'
        % Remove structure(s) already plotted
        remove_structures = listdlg('PromptString','Select a structure to remove:', ...
            'ListString',gui_data.st.safe_name(gui_data.structure_plot_idx));        
        delete(gui_data.handles.structure_patch(remove_structures))
        gui_data.structure_plot_idx(remove_structures) = [];
        gui_data.handles.structure_patch(remove_structures) = [];
        
end

% Upload gui_data
guidata(probe_atlas_gui, gui_data);

end

function update_slice(probe_atlas_gui,varargin)

% Get guidata
gui_data = guidata(probe_atlas_gui);

% Get current position of camera
curr_campos = campos;

% Get probe vector
probe_ref_top = [gui_data.handles.probe_ref_line.XData(1), ...
    gui_data.handles.probe_ref_line.YData(1),gui_data.handles.probe_ref_line.ZData(1)];
probe_ref_bottom = [gui_data.handles.probe_ref_line.XData(2), ...
    gui_data.handles.probe_ref_line.YData(2),gui_data.handles.probe_ref_line.ZData(2)];
probe_vector = probe_ref_top - probe_ref_bottom;

% Get probe-camera vector
probe_camera_vector = probe_ref_top - curr_campos;

% Get the vector to plot the plane in (along with probe vector)
plot_vector = cross(probe_camera_vector,probe_vector);

% Get the normal vector of the plane
normal_vector = cross(plot_vector,probe_vector);

% Get the plane offset through the probe
plane_offset = -(normal_vector*probe_ref_top');

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

% Find plane coordinates that contain brain
curr_slice_isbrain = false(size(use_idx));
curr_slice_isbrain(use_idx) = gui_data.av(curr_slice_idx) > 1;

% Index coordinates in bounds + with brain
grab_pix_idx = sub2ind(size(gui_data.tv),x_idx(curr_slice_isbrain),z_idx(curr_slice_isbrain),y_idx(curr_slice_isbrain));

% Grab pixels from volume
curr_slice = nan(size(use_idx));
curr_slice(curr_slice_isbrain) = gui_data.tv(grab_pix_idx);

% Update the slice display
set(gui_data.handles.slice_plot,'XData',plane_x,'YData',plane_y,'ZData',plane_z,'CData',curr_slice);

% Upload gui_data
guidata(probe_atlas_gui, gui_data);

end

function set_probe_angle(probe_atlas_gui,varargin)

% Get guidata
gui_data = guidata(probe_atlas_gui);

% Prompt for angles
prompt_text = {'Probe ML angle (relative to lambda -> bregma)', ....
    'Probe DV angle (relative to horizontal)'};
probe_angles = ...
    (cellfun(@str2num,inputdlg(prompt_text,'Enter probe angles',1,{'0','90'}))/360)*2*pi;

% Update the probe trajectory reference (reset to bregma)
[ap_max,dv_max,ml_max] = size(gui_data.tv);

max_ref_length = sqrt(sum(([ap_max,dv_max,ml_max].^2)));

[x,y,z] = sph2cart(pi-probe_angles(1),probe_angles(2),max_ref_length);

probe_ref_top = [gui_data.bregma(1),gui_data.bregma(3),0];
probe_ref_bottom = probe_ref_top + [x,y,z];
probe_ref_vector = [probe_ref_top;probe_ref_bottom];

set(gui_data.handles.probe_ref_line,'XData',probe_ref_vector(:,1), ...
    'YData',probe_ref_vector(:,2), ...
    'ZData',probe_ref_vector(:,3));

% Update the slice
update_slice(probe_atlas_gui);

% Upload gui_data
guidata(probe_atlas_gui, gui_data);

end



