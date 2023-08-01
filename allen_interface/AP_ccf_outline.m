%% Draw CCF structures
% Draw 3-view outline of the brain and selected structures

%% Load atlas

allen_atlas_path = fileparts(which('template_volume_10um.npy'));
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean

%% Set up axes, plot brain outline

figure('Color','w');
h = tiledlayout('flow','Tilespacing','none');

% Set up 2D axes
ccf_axes(1) = nexttile;
ccf_axes(2) = nexttile;
ccf_axes(3) = nexttile;

set(ccf_axes,'YDir','reverse')
axis(ccf_axes,'equal','off')
hold(ccf_axes,'on');

% Set up 3D axes
ccf_3d_axes = nexttile;
set(ccf_3d_axes,'ZDir','reverse');
hold(ccf_3d_axes,'on');
axis(ccf_3d_axes,'vis3d','equal','off','manual');
view([-30,25]);
axis tight;
h = rotate3d(ccf_3d_axes);
h.Enable = 'on';

% Plot 3D/2D brain outlines
% (2D)
for curr_view = 1:3
    curr_outline = bwboundaries(squeeze((max(av,[],curr_view)) > 1));
    % (only plot largest outline)
    [~,curr_outline_idx] = max(cellfun(@length,curr_outline));
    curr_outline_reduced = reducepoly(curr_outline{curr_outline_idx});
    plot(ccf_axes(curr_view), ...
        curr_outline_reduced(:,2), ...
        curr_outline_reduced(:,1),'k','linewidth',2);
    % (draw 1mm scalebar)
    %     line(ccf_axes(curr_view),[0,0],[0,100],'color','k','linewidth',2);
end
linkaxes(ccf_axes);

% (3D)
% wireframe
% [~, brain_outline] = plotBrainGrid([],ccf_3d_axes);

% mesh
slice_spacing = 5;
brain_volume = ...
    bwmorph3(bwmorph3(av(1:slice_spacing:end, ...
    1:slice_spacing:end,1:slice_spacing:end)>1,'majority'),'majority');
brain_outline_patchdata = isosurface(permute(brain_volume,[3,1,2]),0.5);
brain_outline = patch( ...
    ccf_3d_axes, ...
    'Vertices',brain_outline_patchdata.vertices*slice_spacing, ...
    'Faces',brain_outline_patchdata.faces, ...
    'FaceColor',[0.7,0.7,0.7],'EdgeColor','none','FaceAlpha',0.1);

%% Get area in search, draw

% Prompt for which structures to show (only structures which are
% labelled in the slice-spacing downsampled annotated volume)
structure_search = lower(inputdlg('Search structures'));
structure_match = find(contains(lower(st.safe_name),structure_search));

selected_structure = listdlg('PromptString','Select a structure to plot:', ...
    'ListString',st.safe_name(structure_match),'ListSize',[520,500], ...
    'SelectionMode','single');

plot_structure = structure_match(selected_structure);

% Get all areas within and below the selected hierarchy level
plot_structure_id = st.structure_id_path{plot_structure};
plot_ccf_idx = find(cellfun(@(x) contains(x,plot_structure_id), ...
    st.structure_id_path));

% Get structure color and volume
slice_spacing = 5;
plot_structure_color = hex2dec(reshape(st.color_hex_triplet{plot_structure},2,[])')./255;
plot_ccf_volume = ismember(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end),plot_ccf_idx);

% Plot 2D structure
for curr_view = 1:3
    curr_outline = bwboundaries(squeeze((max(plot_ccf_volume,[],curr_view))));
    cellfun(@(x) plot(ccf_axes(curr_view),x(:,2)*slice_spacing, ...
        x(:,1)*slice_spacing,'color',plot_structure_color,'linewidth',2),curr_outline)
end

% Plot 3D structure
structure_3d = isosurface(permute(plot_ccf_volume,[3,1,2]),0);
structure_alpha = 0.2;
patch(ccf_3d_axes, ...
    'Vertices',structure_3d.vertices*slice_spacing, ...
    'Faces',structure_3d.faces, ...
    'FaceColor',plot_structure_color,'EdgeColor','none','FaceAlpha',structure_alpha);

%% Draw probes (from AP_histology)

probe_color = 'r';

[probe_ccf_file,probe_ccf_path] = uigetfile('*.mat','Pick probe histology file');
load(fullfile(probe_ccf_path,probe_ccf_file));

% Loop through probes and draw
for curr_probe = 1:length(probe_ccf)

    % Get line of best fit through mean of marked points
    probe_coords_mean = mean(probe_ccf(curr_probe).points,1);
    xyz = bsxfun(@minus,probe_ccf(curr_probe).points,probe_coords_mean);
    [~,~,V] = svd(xyz,0);
    histology_probe_direction = V(:,1);

    % (make sure the direction goes down in DV - flip if it's going up)
    if histology_probe_direction(2) < 0
        histology_probe_direction = -histology_probe_direction;
    end

    % Evaluate line of best fit (length of probe to deepest point)
    [~,deepest_probe_idx] = max(probe_ccf(curr_probe).points(:,2));
    probe_deepest_point = probe_ccf(curr_probe).points(deepest_probe_idx,:);
    probe_deepest_point_com_dist = pdist2(probe_coords_mean,probe_deepest_point);
    probe_length_ccf = 3840/10; % mm / ccf voxel size

    probe_line_eval = probe_deepest_point_com_dist - [probe_length_ccf,0];
    probe_line = (probe_line_eval'.*histology_probe_direction') + probe_coords_mean;

    % Draw probe in 3D view
    line(ccf_3d_axes,probe_line(:,1),probe_line(:,3),probe_line(:,2), ...
        'linewidth',2,'color',probe_color)

    % Draw probes on coronal + saggital
    line(ccf_axes(1),probe_line(:,3),probe_line(:,2),'linewidth',2,'color',probe_color);
    line(ccf_axes(2),probe_line(:,3),probe_line(:,1),'linewidth',2,'color',probe_color);
    line(ccf_axes(3),probe_line(:,2),probe_line(:,1),'linewidth',2,'color',probe_color);

    % Draw probe start/end on horizontal
    plot(ccf_axes(2), probe_line(1,3),probe_line(1,1), ...
        'o','MarkerSize',5,'color',probe_color);
    plot(ccf_axes(2), probe_line(end,3),probe_line(end,1), ...
        '.','MarkerSize',20,'color',probe_color);

end

%% Draw probes (from Neuropixels Trajectory Explorer)

probe_color = 'b';

[probe_ccf_file,probe_ccf_path] = uigetfile('*.mat','Pick probe NTE file');
load(fullfile(probe_ccf_path,probe_ccf_file));

% Loop through probes and draw
for curr_probe = 1:length(probe_positions_ccf)

    % Probe line is directly stored
    probe_line = probe_positions_ccf{curr_probe}';

    % Draw probe in 3D view
    line(ccf_3d_axes,probe_line(:,1),probe_line(:,3),probe_line(:,2), ...
        'linewidth',2,'color',probe_color)

    % Draw probes on coronal + saggital
    line(ccf_axes(1),probe_line(:,3),probe_line(:,2),'linewidth',2,'color',probe_color);
    line(ccf_axes(2),probe_line(:,3),probe_line(:,1),'linewidth',2,'color',probe_color);
    line(ccf_axes(3),probe_line(:,2),probe_line(:,1),'linewidth',2,'color',probe_color);

    % Draw probe start/end on horizontal
    plot(ccf_axes(2), probe_line(1,3),probe_line(1,1), ...
        'o','MarkerSize',5,'color',probe_color);
    plot(ccf_axes(2), probe_line(end,3),probe_line(end,1), ...
        '.','MarkerSize',20,'color',probe_color);

end














