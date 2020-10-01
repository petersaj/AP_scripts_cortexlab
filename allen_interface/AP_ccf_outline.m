%% TO GENERALIZE: Draw outline of CCF structures

%% Load atlas

allen_atlas_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean


%% Outline of VisAM in CCF

% grab slice from gui into curr_ccf_slice
curr_ccf_slice = get(gco,'CData');


am_id = find(contains(st.safe_name,'Anteromedial visual area'));
whitematter_id = find(contains(st.safe_name,{'callosum','commissure', ...
    'radiation','capsule','bundle','ventricle'}));

figure; hold on;
set(gca,'YDir','reverse');
imagesc(curr_ccf_slice);

am_outline = bwboundaries(ismember(curr_ccf_slice,am_id));
whitematter = bwboundaries(ismember(curr_ccf_slice,whitematter_id));
brain_outline = bwboundaries(curr_ccf_slice > 0);

cellfun(@(x) plot(x(:,2),x(:,1),'k','linewidth',2),am_outline);
cellfun(@(x) plot(x(:,2),x(:,1),'k','linewidth',2),whitematter);
cellfun(@(x) plot(x(:,2),x(:,1),'k','linewidth',2),brain_outline);

axis tight image off
