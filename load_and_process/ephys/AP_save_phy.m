function AP_save_phy(confirm)
% AP_save_phy(confirm)
%
% Move clustering data from local phy directory to server
%
% confirm - ask for confirmation to overwrite (default true)

if ~exist('confirm','var') || isempty(confirm)
    confirm = true;
end

% Local phy directory
local_phy_dir = 'C:\data_temp\phy';

% Get animal and day from header
header_filename = [local_phy_dir filesep 'dat_params.txt'];

header_fid = fopen(header_filename);
header_info = textscan(header_fid,'%s %s', 'delimiter',{' = '});
fclose(header_fid);

raw_path_idx = strcmp(header_info{1},'raw_path');
raw_path = header_info{2}{raw_path_idx};

animal_day = regexp(raw_path,'\\Subjects\\(\w*)\\(\d\d\d\d-\d\d-\d\d)','tokens');
animal = animal_day{1}{1};
day = animal_day{1}{2};

% Get server kilosort directory
[ephys_path,ephys_path_exists] = AP_cortexlab_filename(animal,day,[],'ephys');

% TEMPORARY: if kilosort 1 directory, switch to kilosort 2
% (once everything is sorted, AP_cortexlab_filename should default to using
% kilosort2 if it's available and kilosort if it's not)
if ~contains(ephys_path,'kilosort2')
    warning('Temporary: replacing kilosort with kilosort 2')
    ephys_path = strrep(ephys_path,'kilosort','kilosort2');
    if ~exist(ephys_path,'dir')
        error('No kilosort2 directory');
    end
end

% User confirm animal and day overwrite
if confirm
    user_confirm = strcmp(input(['Overwrite phy sorting ' animal ' ' day ' (y/n)? '],'s'),'y');
    if ~user_confirm
        disp('Not saving');
        return
    end
end

% Set standard files to copy (cluster labels and spike clusters)
cluster_files = dir([local_phy_dir filesep 'cluster_*']);
phy_files = [{'spike_clusters.npy'},{cluster_files.name}];

% Move files to server
for curr_file = 1:length(phy_files)
    local_filename = [local_phy_dir filesep phy_files{curr_file}];
    server_filename = [ephys_path filesep phy_files{curr_file}];
    copy_success = copyfile(local_filename,server_filename);
    if copy_success
        disp(['Overwrote ' server_filename]);
    else
        disp(['COULD NOT COPY ' server_filename]);
    end
end

disp('Saved phy clustering to server');




