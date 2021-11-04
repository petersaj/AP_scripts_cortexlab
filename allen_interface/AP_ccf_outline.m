%% Draw CCF structures (IN PROGRESS)
%
% This is disorganized code to draw CCF regions
% Turn this into a function 


%% Load atlas

allen_atlas_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean

%% Out line of coronal/horizontal

figure;

% Plot brain to overlay probes
horizontal_axes = subplot(1,3,1,'YDir','reverse'); hold on;
axis image off;
coronal_axes = subplot(1,3,2,'YDir','reverse'); hold on;
axis image off;
sagittal_axes = subplot(1,3,3,'YDir','reverse'); hold on;
axis image off;

% Plot projections
coronal_outline = bwboundaries(permute((max(av,[],1)) > 1,[2,3,1]));
horizontal_outline = bwboundaries(permute((max(av,[],2)) > 1,[3,1,2]));
sagittal_outline = bwboundaries(permute((max(av,[],3)) > 1,[2,1,3]));

str_id = find(strcmp(st.safe_name,'Caudoputamen'));
str_coronal_outline = bwboundaries(permute((max(av == str_id,[],1)) > 0,[2,3,1]));
str_horizontal_outline = bwboundaries(permute((max(av == str_id,[],2)) > 0,[3,1,2]));
str_sagittal_outline = bwboundaries(permute((max(av == str_id,[],3)) > 0,[2,1,3]));

vstr_id = find(strcmp(st.safe_name,'Nucleus accumbens'));
vstr_coronal_outline = bwboundaries(permute((max(av == vstr_id,[],1)) > 0,[2,3,1]));

cellfun(@(x) plot(horizontal_axes,x(:,2),x(:,1),'k','linewidth',2),horizontal_outline)
cellfun(@(x) plot(horizontal_axes,x(:,2),x(:,1),'b','linewidth',2),str_horizontal_outline)
axis image off;

cellfun(@(x) plot(coronal_axes,x(:,2),x(:,1),'k','linewidth',2),coronal_outline)
cellfun(@(x) plot(coronal_axes,x(:,2),x(:,1),'b','linewidth',2),str_coronal_outline)
cellfun(@(x) plot(coronal_axes,x(:,2),x(:,1),'color',[0.8,0,0.8],'linewidth',2),vstr_coronal_outline)
axis image off;

cellfun(@(x) plot(sagittal_axes,x(:,2),x(:,1),'k','linewidth',2),sagittal_outline)
cellfun(@(x) plot(sagittal_axes,x(:,2),x(:,1),'b','linewidth',2),str_sagittal_outline)
axis image off;

%% Get area in hierarchy, draw

% Bring up hierarchical selector
plot_structure = hierarchicalSelect(st);


%% Outline of VisAM in CCF

% TO USE THIS: open allen_ccf_npx (might also work on allenAtlasBrowser,
% but that seems broken at the moment?), move to desired annotated slice,
% click once, then the code below will pull that slice and draw outlines

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
