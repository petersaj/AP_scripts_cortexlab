function injection_coordinates = AP_get_allen_projection(projection_coordinates,max_coords)
% injection_coordinates = AP_get_allen_projection(projection_coordinates)
%
% Load in from Allen projection data
% INPUTS
% projection_coordinates - n x 3 array of coordinates to query
% max_coords - the maximum number of injection coordinates to return
% OUTPUTS
% injection_coordinates - cell array of cortical injection sites

allen_temp_dir = 'C:\data_temp\AllenAPI';

% Loop through all input coordinates
disp('Getting injection coordinates...')
injection_coordinates = cell(size(projection_coordinates,1),1);
for curr_coord = 1:size(projection_coordinates,1);
    
    disp(['Coordinate ' num2str(curr_coord) '/' num2str(size(projection_coordinates,1))]);

    % Download spatial search
    seed_point = projection_coordinates(curr_coord,:);
    seed_url = sprintf('http://api.brain-map.org/api/v2/data/query.xml?criteria=service::mouse_connectivity_target_spatial[seed_point$eq%d,%d,%d][num_rows$eq%d][primary_structure_only$eqtrue][injection_structures$eqIsocortex]',seed_point,max_coords);
    seed_table_fn = [allen_temp_dir filesep 'seed_experiments'];
    urlwrite(seed_url,seed_table_fn);
    seed_table_raw = xml2struct(seed_table_fn);
    projection_experiments = seed_table_raw.Response.objects.object;
    
    % Delete downloaded data
    delete(seed_table_fn);
    
    % Get coordinates of all injections
    injection_coordinates{curr_coord} = nan(length(projection_experiments),3);
    for curr_experiment = 1:length(projection_experiments)
        % Injection coordinates are /10 to fit with Allen CCF
        injection_coordinates{curr_coord}(curr_experiment,:) = ...
            cellfun(@(x) str2num(x.Text), ...
            projection_experiments{curr_experiment}.injection_dash_coordinates.injection_dash_coordinate);
    end    
    
end
disp('Done');







