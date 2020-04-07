function injection_parameters = get_allen_projection(projection_coordinates)
% injection_coordinates = AP_get_allen_projection(projection_coordinates)
%
% Load in from Allen projection data
% INPUTS
% projection_coordinates - n x 3 array of coordinates to query
% OUTPUTS
% injection_coordinates - cell array of cortical injection sites

allen_temp_dir = 'C:\data_temp\AllenAPI';
if ~exist(allen_temp_dir,'dir')
    mkdir(allen_temp_dir)
end

% Loop through all input coordinates
disp('Getting injection coordinates...')
injection_parameters = struct('coordinates',cell(size(projection_coordinates,1),1), ...
    'volume',cell(size(projection_coordinates,1),1), ...
    'density',cell(size(projection_coordinates,1),1));
for curr_coord = 1:size(projection_coordinates,1)
    
    disp(['Coordinate ' num2str(curr_coord) '/' num2str(size(projection_coordinates,1))]);

    % Download spatial search
    seed_point = projection_coordinates(curr_coord,:);
    seed_url = sprintf('http://api.brain-map.org/api/v2/data/query.xml?criteria=service::mouse_connectivity_target_spatial[seed_point$eq%d,%d,%d][primary_structure_only$eqtrue][injection_structures$eqIsocortex]',seed_point);
    seed_table_fn = [allen_temp_dir filesep 'seed_experiments'];
    urlwrite(seed_url,seed_table_fn);
    seed_table_raw = xml2struct(seed_table_fn);
    projection_experiments = seed_table_raw.Response.objects.object;
    
    % Delete downloaded data
    delete(seed_table_fn);
    
    % Get coordinates, injection volume, and projection density
    for curr_experiment = 1:length(projection_experiments)
        % Injection coordinates are /10 to fit with Allen CCF
        injection_parameters(curr_coord).coordinates(curr_experiment,:) = ...
            cellfun(@(x) str2num(x.Text), ...
            projection_experiments{curr_experiment}.injection_dash_coordinates.injection_dash_coordinate);
        % Injection volume
        injection_parameters(curr_coord).volume(curr_experiment) = ...
            str2num(projection_experiments{curr_experiment}.injection_dash_volume.Text);
        % Projection density
        injection_parameters(curr_coord).density(curr_experiment) = ...
            str2num(projection_experiments{curr_experiment}.density.Text);
    end    
        
end
disp('Done');







