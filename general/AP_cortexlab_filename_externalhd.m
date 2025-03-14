function [filename,file_exists] = AP_cortexlab_filename_externalhd(animal,day,experiment,file,site,recording)
% [filename,file_exists] = AP_cortexlab_filename_externalhd(animal,day,experiment,file,site,recording)
%
% This is an absolute mess because of lab-wide inconsistency
%
% file types:
% expInfo
% timeline
% block
% parameters
% protocol
% eyecam
% eyecam_processed
% eyecam_dlc
% facecam
% facecam_processed
% facecam_movement
% hardware
% imaging
% ephys
% ephys_ks1
% ephys_dir
% ephys_ap
% probe_ccf
%
% EXTERNALHD VERSION: looking for data on external hd, after leaving lab so
% no server access

% Make inputs strings if they're numbers
if isnumeric(experiment)
    experiment = num2str(experiment);
end

if isnumeric(day)
    experiment = num2str(day);
end

% Site = multiple probes
if exist('site','var') && ~isempty(site)
    if isnumeric(site)
        site = num2str(site);
    end
    site_dir = [filesep 'site' site];
else
    site = [];
    site_dir = [];
end

% Recording = separated recordings from single probe
if exist('recording','var') && ~isempty(recording)
    if isnumeric(recording)
        recording = num2str(recording);
    end
    recording_dir = [filesep 'experiment' recording];
else
    recording = [];
    recording_dir = [];
end

% List all folders to check
server_location = cell(0);
% (zserver: different files used to be split across folders)
server_location{end+1} = '\\qnap-ap001.dpag.ox.ac.uk\APlab\Users\Andy_Peters\cortexlab_corticostriatal_data';

% Check that servers are accessible (login needed on restart)
warning on
for curr_location = 1:length(server_location)
   if ~exist(server_location{curr_location},'dir')
       error([server_location{curr_location} ' not available']);
   end
end

switch file
    
    case 'expInfo'
        filepattern = [animal filesep day filesep];
        [filename,file_exists] = check_locations(filepattern,server_location);
        if file_exists
            filename = fileparts(filename{1});
        end
        
    case 'timeline'
        filepattern = [animal filesep day filesep experiment ...
            filesep day '_' experiment '_' animal '_Timeline.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
    case 'block'
        filepattern = [animal filesep day filesep experiment ...
            filesep day '_' experiment '_' animal '_Block.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
    case 'parameters'
        filepattern = [animal filesep day filesep experiment ...
            filesep day '_' experiment '_' animal '_parameters.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
    case 'protocol'
        filepattern = [animal filesep day filesep experiment filesep 'Protocol.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
    case 'eyecam'
        filepattern = [animal filesep day filesep experiment filesep 'eye.mj2'];
        [filename,file_exists] = check_locations(filepattern,server_location);

    case 'eyecam_processed'
        filepattern = [animal filesep day filesep experiment filesep 'eye_proc.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);

    case 'eyecam_dlc'
        filepattern = [animal filesep day filesep experiment filesep 'eyeDLC*.csv'];
        [filename,file_exists] = check_locations(filepattern,server_location);

    case 'facecam'
        filepattern = [animal filesep day filesep experiment filesep 'face.mj2'];
        [filename,file_exists] = check_locations(filepattern,server_location);

    case 'facecam_processed'
        filepattern = [animal filesep day filesep experiment filesep 'face_proc.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
    case 'facecam_movement'
        % (output from AP_mouse_movie_movement)
        filepattern = [animal filesep day filesep experiment filesep 'facecam_movement.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
    case 'hardware'
        filepattern = [animal filesep day filesep experiment ...
            filesep day '_' experiment '_' animal '_hardwareInfo.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
    case 'imaging'
        filepattern = [animal filesep day filesep filesep 'svd*'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        if file_exists && iscell(filename)
            filename = fileparts(filename{1});
        elseif file_exists && isstr(filename)
            filename = fileparts(filename);
        end
        
    case 'ephys_dir'
        % (the path where the ephys data is kept)
        filepattern = [animal filesep day filesep 'ephys' site_dir];
        [filename,file_exists] = check_locations(filepattern,server_location);
        if file_exists
            filename = fileparts(filename{1});
        end
        
    case 'ephys_ap'
        % (the raw action potential band data file)
        
        % Old open ephys
        filepattern = [animal filesep day filesep 'ephys' site_dir filesep 'experiment1_10*-0_0.dat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
        % New open ephys
        if ~file_exists
            filepattern = [animal filesep day filesep ...
                'ephys' site_dir filesep 'experiment1' filesep 'recording1' ...
                filesep 'continuous'  filesep 'Neuropix-3a-100.0' filesep 'continuous.dat'];
            [filename,file_exists] = check_locations(filepattern,server_location);
        end
        
    case 'ephys'
        % (folder with kilosort/phy outputs)
        
        kilosort_version = 2; % (kilosort 2 by default)
        
        % Drop the kilosort version in the base workspace
        assignin('base','kilosort_version',kilosort_version);
        
        filepattern = [animal filesep day filesep 'ephys'  filesep 'kilosort2' filesep site_dir filesep recording_dir];
        [filename,file_exists] = check_locations(filepattern,server_location);
        if file_exists
            filename = fileparts(filename{1});
        end
        
    case 'ephys_ks1'
        % folder with kilosort/phy outputs
        
        kilosort_version = 1; % (kilosort 2 by default)
        
        % Drop the kilosort version in the base workspace
        assignin('base','kilosort_version',kilosort_version);
        
        filepattern = [animal filesep day filesep 'ephys' filesep 'kilosort' filesep site_dir filesep recording_dir];
        [filename,file_exists] = check_locations(filepattern,server_location);
        if file_exists
            filename = fileparts(filename{1});
        end
        
    case 'probe_ccf'
        % Histology probe location from AP histology
        % (sometimes upper/lowecase "Histology" folder)
        filepattern = [animal filesep '*istology' filesep 'slices' filesep 'probe_ccf.mat'];
        [filename,file_exists] = check_locations(filepattern,server_location);
        
end
end


function [filename,file_exists] = check_locations(filepattern,server_location)

% Loop through all server locations and look for file
for curr_location = 1:length(server_location)
    curr_filename = [server_location{curr_location} filesep filepattern];
    curr_filename_dir = dir(curr_filename);
    file_exists = ~isempty(curr_filename_dir);
    
    if file_exists
        % If found, break and use this filename
        if length(curr_filename_dir) == 1
            filename = [curr_filename_dir.folder filesep curr_filename_dir.name];
        else
            filename = cellfun(@(path,fn) [path filesep fn], ...
                {curr_filename_dir.folder},{curr_filename_dir.name},'uni',false);
        end
        break
    else
        % If not found, clear filename
        filename = [];
    end
end

end





