% function AP_vfs_ccf_align
% 
% Not really sure where to put this code at the moment or what to do with
% it, but this aligns a widefield visual field sign map to the allen CCF
% boundaries, requires imaged visual field sign 'imaged_vfs'

%% Make top-down boundaries from scratch (I save these and load)

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

plot_areas = intersect(used_areas,ctx_idx);

bregma = allenCCFbregma;

% Get outlines of all areas
top_down_cortical_area_boundaries = cell(size(plot_areas));
for curr_area_idx = 1:length(plot_areas)
    top_down_cortical_area_boundaries{curr_area_idx} = bwboundaries(top_down_annotation == plot_areas(curr_area_idx));
end

%% Label the visual areas according to visual field sign

al_idx = find(cellfun(@(name) strcmp(name,'Anterolateral visual area layer 1'),st.safe_name(used_areas)));
am_idx = find(cellfun(@(name) strcmp(name,'Anteromedial visual area layer 1'),st.safe_name(used_areas)));
lm_idx = find(cellfun(@(name) strcmp(name,'Lateral visual area layer 1'),st.safe_name(used_areas)));
v1_idx = find(cellfun(@(name) strcmp(name,'Primary visual area layer 1'),st.safe_name(used_areas)));
p_idx = find(cellfun(@(name) strcmp(name,'Posterolateral visual area layer 1'),st.safe_name(used_areas)));
pm_idx = find(cellfun(@(name) strcmp(name,'posteromedial visual area layer 1'),st.safe_name(used_areas)));
li_idx = find(cellfun(@(name) strcmp(name,'Laterointermediate area layer 1'),st.safe_name(used_areas)));
rl_idx = find(cellfun(@(name) strcmp(name,'Rostrolateral area layer 1'),st.safe_name(used_areas)));

ccf_vfs = zeros(size(top_down_annotation));
ccf_vfs(ismember(top_down_annotation,used_areas([v1_idx,am_idx,al_idx,li_idx]))) = 1;
ccf_vfs(ismember(top_down_annotation,used_areas([p_idx,pm_idx,rl_idx,lm_idx]))) = -1;

%% Align the imaged VFS to the (downsampled) Allen VFS

% Hard-code microns per pixel
um2pixel = 20.6;

ccf_vfs_d = imresize(ccf_vfs,10/um2pixel,'nearest');

[optimizer, metric] = imregconfig('monomodal');
optimizer = registration.optimizer.OnePlusOneEvolutionary();
optimizer.MaximumIterations = 200;
optimizer.GrowthFactor = 1+1e-6;
optimizer.InitialRadius = 1e-4;

tformEstimate_affine = imregtform(imaged_vfs,ccf_vfs_d,'affine',optimizer,metric);
curr_im_reg = imwarp(imaged_vfs,tformEstimate_affine,'Outputview',imref2d(size(ccf_vfs_d)));
tform_matrix = tformEstimate_affine.T;

%% Plot the CCF boundaries over the aligned VFS

figure; 
colormap(colormap_BlueWhiteRed);

subplot(1,2,1);
imagesc(ccf_vfs_d); hold on; axis image off;
for curr_area_idx =1:length(top_down_cortical_area_boundaries)
    cellfun(@(outline) plot((outline(:,2)*10)/um2pixel, ...
        (outline(:,1)*10)/um2pixel,'k'),top_down_cortical_area_boundaries{curr_area_idx},'uni',false);
end

subplot(1,2,2);
imagesc(curr_im_reg); hold on; axis image off;
for curr_area_idx =1:length(top_down_cortical_area_boundaries)
    cellfun(@(outline) plot((outline(:,2)*10)/um2pixel, ...
        (outline(:,1)*10)/um2pixel,'k'),top_down_cortical_area_boundaries{curr_area_idx},'uni',false);
end

%% ALTERNATELY: Align the (downsampled) Allen VFS to the imaged VFS

ccf_vfs_d = imresize(ccf_vfs,10/um2pixel,'nearest');

[optimizer, metric] = imregconfig('monomodal');
optimizer = registration.optimizer.OnePlusOneEvolutionary();
optimizer.MaximumIterations = 200;
optimizer.GrowthFactor = 1+1e-6;
optimizer.InitialRadius = 1e-4;

tformEstimate_affine = imregtform(ccf_vfs_d,imaged_vfs,'affine',optimizer,metric);
tform_matrix = tformEstimate_affine.T;

boundaries_tform_long = cellfun(@(areas) cellfun(@(coords) ...
    [fliplr(coords*(10/um2pixel)),ones(size(coords,1),1)]*tform_matrix,areas,'uni',false), ...
    top_down_cortical_area_boundaries,'uni',false);

cortical_boundaries_aligned = cellfun(@(areas) cellfun(@(coords) ...
    coords(:,[2,1]),areas,'uni',false),boundaries_tform_long,'uni',false);

figure; 
imagesc(imaged_vfs); hold on; axis image off;
colormap(colormap_BlueWhiteRed);
for curr_area_idx = 1:length(cortical_boundaries_aligned)
    cellfun(@(outline) plot(outline(:,2),outline(:,1),'k'), ...
        cortical_boundaries_aligned{curr_area_idx},'uni',false);
end







