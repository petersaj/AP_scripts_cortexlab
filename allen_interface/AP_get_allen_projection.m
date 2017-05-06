%% Load in from Allen projection data

% Download spatial search
allen_temp_dir = 'C:\data_temp\AllenAPI';
seed_point = [4950,3300,6800];
seed_example_url = sprintf('http://api.brain-map.org/api/v2/data/query.xml?criteria=service::mouse_connectivity_target_spatial[seed_point$eq%d,%d,%d][start_row$eq25][num_rows$eq50][primary_structure_only$eqtrue][injection_structures$eqIsocortex]',seed_point);
seed_table_fn = [allen_temp_dir filesep 'seed_experiments'];
urlwrite(seed_example_url,seed_table_fn);
seed_table_raw = xml2struct(seed_table_fn);
projection_experiments = seed_table_raw.Response.objects.object;
% Delete downloaded data
delete(seed_table_fn);

% Get coordinates of all injections
injection_coordinates = nan(length(projection_experiments),3);
for curr_experiment = 1:length(projection_experiments)
    % Injection coordinates are /10 to fit with Allen CCF
    injection_coordinates(curr_experiment,:) = ...
        cellfun(@(x) str2num(x.Text), ...
        projection_experiments{curr_experiment}.injection_dash_coordinates.injection_dash_coordinate);
end

% Plot the injection points on a dorsal outline of the brain
dorsal_brain = permute((max(av,[],2)) > 1,[3,1,2]);
[pt1,pt2] = ind2sub(size(dorsal_brain),find(dorsal_brain,1));

dorsal_brain_outline = bwtraceboundary(dorsal_brain,[pt1,pt2],'N');

figure; axes('YDir','reverse'); hold on; axis equal; axis off;
plot(dorsal_brain_outline(:,1),dorsal_brain_outline(:,2),'k');
plot(injection_coordinates(:,3)/10,injection_coordinates(:,1)/10,'.r')

% Bin the injection points and make a heatmap
bin_size = 1;
n_inj = hist3(injection_coordinates(:,[1,3])/10,'Edges',{0:bin_size:size(av,1),0:bin_size:size(av,3)});
gauss_sigma = 30;
n_inj_blur = imgaussfilt(n_inj,repmat(gauss_sigma,1,2));
figure; axes('YDir','reverse'); hold on; colormap(flipud(gray));
imagesc(n_inj_blur);
plot(dorsal_brain_outline(:,1),dorsal_brain_outline(:,2),'k');












