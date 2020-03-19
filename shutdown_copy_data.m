%% Copy all data to local drive for shutdown prep (not raw ephys)


%% Set animals

trained_animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
naive_animals = {'AP032','AP033','AP034','AP035','AP036'};
muscimol_animals = {'AP045','AP054','AP055','AP053','AP047','AP048'};
cortical_animals = {'AP043','AP060','AP061'};
corticostriatal_animals = {'AP063','AP068'};

animals = [trained_animals, ...
    naive_animals, ...
    muscimol_animals, ...
    cortical_animals, ...
    corticostriatal_animals];

%% Set local and server locations

% Local path
local_drive = 'E:\full_data_backup';

% List servers
server1 = '\\zserver.cortexlab.net';
server2 = '\\zubjects.cortexlab.net';
server3 = '\\znas.cortexlab.net';

% Check that servers are accessible (login needed on restart)
if ~exist([server1 filesep 'Data'])
    error('Zserver not available');
end
if ~exist([server2 filesep 'Subjects'])
    error('Zubjects not available');
end
if ~exist([server3 filesep 'Subjects'])
    error('Znas not available');
end

% List all folders to check
server_location = cell(0);
server_location{end+1} = [server3 filesep 'Subjects'];
server_location{end+1} = [server2 filesep 'Subjects'];
server_location{end+1} = [server1 filesep 'Data' filesep 'Subjects'];
server_location{end+1} = [server1 filesep 'Data' filesep 'expInfo'];
server_location{end+1} = [server1 filesep 'Data' filesep 'trodes'];
server_location{end+1} = [server1 filesep 'Data' filesep 'EyeCamera'];


%% Copy from server to local (but NOT raw ephys)

% Loop through animals
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
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
                sprintf('%s \n ---> %s',curr_server_pull,curr_local_drop)
                copyfile(curr_server_pull,curr_local_drop);
            end
            
            % If ephys, copy everything not .dat
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
                    
                    sprintf('%s \n ---> %s',curr_server_pull,curr_local_drop)
                    copyfile(curr_server_pull,curr_local_drop);            
                end
            end
        end
    end
end











