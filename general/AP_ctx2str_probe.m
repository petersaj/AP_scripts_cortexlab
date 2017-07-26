%% Get projections from the cortex to striatum 
% Uses probe location from allen_atlas_probe

% Get probe location per micron
probe_size = pdist2(probe_vector_ccf(1,:),probe_vector_ccf(2,:));
probe_depths = ...
    round([linspace(probe_vector_ccf(1,1)',probe_vector_ccf(2,1)',probe_size); ...
    linspace(probe_vector_ccf(1,2)',probe_vector_ccf(2,2)',probe_size); ...
    linspace(probe_vector_ccf(1,3)',probe_vector_ccf(2,3)',probe_size)]');

% Load in the annotated Allen volume
av_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\annotation_volume_10um_by_index.npy';
av = readNPY(av_fn);

% Convert probe depths subscripts to indicies
probe_depths_ind = sub2ind(size(av),round(probe_depths(:,1)/10),round(probe_depths(:,2)/10),round(probe_depths(:,3)/10));

% Get structures that the probe is in
probe_structures = av(probe_depths_ind);

% Get regular depths from probe that runs through striatum
str_id = 573; % just hardcode for now
n_str_depths = 6;
probe_str = [probe_depths(find(probe_structures == str_id,1,'first'),:); ...
    probe_depths(find(probe_structures == str_id,1,'last'),:)];
str_depths = ...
    round([linspace(probe_str(1,1)',probe_str(2,1)',n_str_depths); ...
    linspace(probe_str(1,2)',probe_str(2,2)',n_str_depths); ...
    linspace(probe_str(1,3)',probe_str(2,3)',n_str_depths)]');

% Mirror all of the locations across the midline and get projections from
% both (because the cortical injections in the database aren't evenly
% distributed, so this is more comprehensive)
ccf_midline = 5700;
hemisphere = sign(str_depths(1,3) - ccf_midline);
str_depths_mirror = [str_depths(:,1:2),ccf_midline - hemisphere*abs(str_depths(:,3)-ccf_midline)];

max_sites = 50;
str_depths_query = [str_depths;str_depths_mirror];
injection_coordinates = get_allen_projection(str_depths_query,max_sites);

% Standardize injection coordinates by hemisphere (left = contra, right =
% ipsi)
injection_coordinates_standardized = injection_coordinates;
for curr_coord = 1:length(injection_coordinates)
    
    target_hemisphere = sign(str_depths_query(curr_coord,3) - ccf_midline);
    injection_coords_ml_offset = abs(injection_coordinates{curr_coord}(:,3) - ccf_midline);
    injection_coordinates_hemisphere = sign(injection_coordinates{curr_coord}(:,3) - ccf_midline);
    
    injection_coords_ipsi = injection_coordinates_hemisphere == target_hemisphere;
    injection_coords_contra = injection_coordinates_hemisphere == -target_hemisphere;
    
    injection_coordinates_standardized{curr_coord}(injection_coords_ipsi,3) = ...
        ccf_midline + injection_coords_ml_offset(injection_coords_ipsi);
    injection_coordinates_standardized{curr_coord}(injection_coords_contra,3) = ...
        ccf_midline - injection_coords_ml_offset(injection_coords_contra);
    
end

% Plot the injection sites
figure; 

% By site
plot_colors = repmat(jet(size(str_depths,1)),2,1);
dorsal_brain = permute((max(av,[],2)) > 1,[3,1,2]);
[pt1,pt2] = ind2sub(size(dorsal_brain),find(dorsal_brain,1));
dorsal_brain_outline = bwtraceboundary(dorsal_brain,[pt1,pt2],'N');

subplot(1,2,1,'YDir','reverse'); hold on; axis equal; axis off;
ylim([1,size(dorsal_brain,2)]);
xlim([1,size(dorsal_brain,1)]);
plot(dorsal_brain_outline(:,1),dorsal_brain_outline(:,2),'k');
line(repmat(ccf_midline/10,1,2),ylim,'color','k');
for curr_coord = 1:length(injection_coordinates)
    plot(injection_coordinates_standardized{curr_coord}(:,3)/10, ...
        injection_coordinates_standardized{curr_coord}(:,1)/10, ...
        '.','color',plot_colors(curr_coord,:),'MarkerSize',10);    
end
text(ccf_midline/10-500,100,'Contralateral','HorizontalAlignment','center')
text(ccf_midline/10+500,100,'Ipsilateral','HorizontalAlignment','center')
title('Injection sites')

% As binned COM
bin_size = 1;
bin_edges = {0:bin_size:size(av,1),0:bin_size:size(av,3)};
n_inj = nan(length(bin_edges{1}),length(bin_edges{2}),length(injection_coordinates));
for curr_coord = 1:length(injection_coordinates)    
    n_inj(:,:,curr_coord) = hist3(injection_coordinates_standardized{curr_coord}(:,[1,3])/10,'Edges',bin_edges);
end

n_inj_combined_hemisphere = n_inj(:,:,1:size(str_depths,1)) + n_inj(:,:,size(str_depths,1)+1:end);

gauss_sigma = 30;
n_inj_blur = imgaussfilt(n_inj_combined_hemisphere,repmat(gauss_sigma,1,2));

r_px_com = sum(bsxfun(@times,n_inj_blur,permute(1:size(str_depths,1),[1,3,2])),3)./sum(n_inj_blur,3);

r_px_com_col = ind2rgb(round(mat2gray(r_px_com)*255),jet(255));

subplot(1,2,2,'YDir','reverse'); hold on; axis equal; axis off;
plot(dorsal_brain_outline(:,1),dorsal_brain_outline(:,2),'k');
line(repmat(ccf_midline/10,1,2),ylim,'color','k');
ylim([1,size(dorsal_brain,2)]);
xlim([1,size(dorsal_brain,1)]);
p = imagesc(r_px_com_col);
set(p,'AlphaData',mat2gray(max(n_inj_blur,[],3), ...
     [0,double(prctile(reshape(max(n_inj_blur,[],3),[],1),99))]));
text(ccf_midline/10-500,100,'Contralateral','HorizontalAlignment','center')
text(ccf_midline/10+500,100,'Ipsilateral','HorizontalAlignment','center')
title('Injection site center-of-mass')






