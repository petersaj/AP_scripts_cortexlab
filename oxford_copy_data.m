%% Copy all data to local drive for UCL to Oxford move (not raw ephys)
%
% Just data from operant task and unused corticostriatal imaging
% (from shutdown_copy_data.m)


%% Set animals

corticostriatal_animals = {'AP075','AP077','AP079','AP085','AP087', ...
    'AP089','AP090','AP091','AP092','AP093','AP094','AP095','AP096','AP097'};

operant_animals = {'AP100','AP101','AP103','AP104','AP105','AP106', ...
    'AP107','AP108','AP109','AP110','AP111','AP112','AP113','AP114', ...
    'AP115'};

animals = [ ...
    corticostriatal_animals, ...
    operant_animals, ...
    ];

%% Set local and server locations

% Drive to copy
local_drive = 'D:\CortexlabData';

server_location = cell(0);
% (zserver: different files used to be split across folders)
server_location{end+1} = '\\zserver.cortexlab.net\Data\Subjects';
server_location{end+1} = '\\zserver.cortexlab.net\Data\expInfo';
server_location{end+1} = '\\zserver.cortexlab.net\Data\trodes';
server_location{end+1} = '\\zserver.cortexlab.net\Data\EyeCamera';
% (servers after zserver filled)
server_location{end+1} = '\\zubjects.cortexlab.net\Subjects';
server_location{end+1} = '\\znas.cortexlab.net\Subjects';
server_location{end+1} = '\\zinu.cortexlab.net\Subjects';
server_location{end+1} = '\\zaru.cortexlab.net\Subjects';

% Check that servers are accessible (login needed on restart)
warning on
for curr_location = 1:length(server_location)
   if ~exist(server_location{curr_location},'dir')
       error([server_location{curr_location} ' not available']);
   end
end

%% Copy from server to local (but NOT raw ephys)

% Loop through animals
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    fprintf('Copying %s... \n', animal);
    tic;
    
    % Loop through server locations
    for curr_server_location = 1:length(server_location)
        
        curr_path = [server_location{curr_server_location} filesep animal];
        curr_dir = dir(curr_path);
        day_paths = cellfun(@(x) ~isempty(regexp(x,'\d\d\d\d-\d\d-\d\d')),{curr_dir.name}) &...
            [curr_dir.isdir];
        curr_days = {curr_dir(day_paths).name};
        
        for curr_day_idx = 1:length(curr_days)
            curr_day = curr_days{curr_day_idx};
            curr_day_dir = dir([curr_path filesep curr_day]);
            
            % Copy everything not ephys
            copy_idx = setdiff(find(~contains({curr_day_dir.name},{'ephys'})),[1,2]);          
            for curr_copy = copy_idx
                curr_server_pull = [curr_path filesep curr_day filesep curr_day_dir(curr_copy).name];
                curr_local_drop = [local_drive filesep animal filesep curr_day filesep curr_day_dir(curr_copy).name];
                fprintf('%s \n ---> %s',curr_server_pull,curr_local_drop)
                copyfile(curr_server_pull,curr_local_drop);
            end
            
            % If ephys, only copy kilosort output
            ephys_copy_idx = strmatch('ephys',{curr_day_dir.name});
            if any(ephys_copy_idx)
                ephys_copy_dir = dir([curr_path filesep curr_day filesep curr_day_dir(ephys_copy_idx).name]);
                ephys_copy_dir_idx = setdiff(find(~contains({ephys_copy_dir.name},'.dat')),[1,2]);
                for curr_copy = ephys_copy_dir_idx
                    curr_server_pull = [ephys_copy_dir(curr_copy).folder filesep ephys_copy_dir(curr_copy).name];
                    curr_local_drop = [local_drive filesep animal filesep curr_day filesep 'ephys' filesep ephys_copy_dir(curr_copy).name];
                    
                    if ~exist(fileparts(curr_local_drop),'dir')
                        mkdir(fileparts(curr_local_drop));
                    end
                    
                    fprintf('%s \n ---> %s',curr_server_pull,curr_local_drop)
                    copyfile(curr_server_pull,curr_local_drop);            
                end
            end
        end
    end

    fprintf('Done %s (%d/%d animals) \n', animal,curr_animal,length(animals));
    toc;
    
end








