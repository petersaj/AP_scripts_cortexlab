%% Draw CCF structures
% Draw 3-view outline of the brain and selected structures
%
% TO DO: CCF DV axis is scaled??? correct for that

%% Load atlas

allen_atlas_path = fileparts(which('template_volume_10um.npy'));
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean

%% Set up axes, plot brain outline

% Set up axes
figure;
ccf_axes = gobjects(3,1);
ccf_axes(1) = subplot(1,3,1,'YDir','reverse');
hold on; axis image off;
ccf_axes(2) = subplot(1,3,2,'YDir','reverse');
hold on; axis image off;
ccf_axes(3) = subplot(1,3,3,'YDir','reverse');
hold on; axis image off;

% Plot brain projections
for curr_view = 1:3
    curr_outline = bwboundaries(squeeze((max(av,[],curr_view)) > 1));
    cellfun(@(x) plot(ccf_axes(curr_view),x(:,2),x(:,1),'k','linewidth',2),curr_outline)
    % (draw 1mm scalebar)
%     line(ccf_axes(curr_view),[0,0],[0,100],'color','k','linewidth',2);
end
linkaxes(ccf_axes);


%% Get area in hierarchy, draw

% Bring up hierarchical selector
plot_structure = hierarchicalSelect(st);

% Get all areas within and below the selected hierarchy level
plot_structure_id = st.structure_id_path{plot_structure};
plot_ccf_idx = find(cellfun(@(x) contains(x,plot_structure_id), ...
    st.structure_id_path));

% Get structure CCF color
plot_structure_color = hex2dec(reshape(st.color_hex_triplet{plot_structure},2,[])')./255;

% Get structure volume
plot_ccf_volume = ismember(av,plot_ccf_idx);

for curr_view = 1:3
    curr_outline = bwboundaries(squeeze((max(plot_ccf_volume,[],curr_view))));
    cellfun(@(x) plot(ccf_axes(curr_view),x(:,2),x(:,1),'color',plot_structure_color,'linewidth',2),curr_outline)
end


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

% plot the structure
slice_spacing = 5;
plot_structure_color = hex2dec(reshape(st.color_hex_triplet{plot_structure},2,[])')./255;

% Get structure volume
plot_ccf_volume = ismember(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end),plot_ccf_idx);

for curr_view = 1:3
    curr_outline = bwboundaries(squeeze((max(plot_ccf_volume,[],curr_view))));
    cellfun(@(x) plot(ccf_axes(curr_view),x(:,2)*slice_spacing, ...
        x(:,1)*slice_spacing,'color',plot_structure_color,'linewidth',2),curr_outline)
end


