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
local_phy_dir = 'G:\data_temp\phy';

% Get animal and day from header
header_filename = [local_phy_dir filesep 'dat_params.txt'];

header_fid = fopen(header_filename);
header_info = textscan(header_fid,'%s %s', 'delimiter',{' = '});
fclose(header_fid);

raw_path_idx = strcmp(header_info{1},'raw_path');
raw_path = header_info{2}{raw_path_idx};

animal_day_site = regexp(raw_path,'\\Subjects\\(\w*)\\(\d\d\d\d-\d\d-\d\d)\\ephys([^'']*)','tokens');
animal = animal_day_site{1}{1};
day = animal_day_site{1}{2};
site_dir = animal_day_site{1}{3};

% Get server kilosort directory
[ephys_path,ephys_path_exists] = AP_cortexlab_filename(animal,day,[],'ephys');
kilosort_path = [ephys_path filesep site_dir];

% Check that the path exists (it must because it's where kilosort lives -
% if it doesn't it means the path is wrong)
if ~exist(kilosort_path,'dir')
    error(['Kilosort path nonexistant: ' kilosort_path])
end

% User confirm animal and day overwrite
if confirm
    confirm_text = ['Save phy sorting for ' animal ' ' day ' to ' kilosort_path ' (y/n)? '];
    confirm_text_slash = strrep(confirm_text,'\','\\');
    user_confirm = strcmp(input(confirm_text_slash,'s'),'y');
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
    server_filename = [kilosort_path filesep phy_files{curr_file}];
    copy_success = copyfile(local_filename,server_filename);
    if copy_success
        disp(['Overwrote ' server_filename]);
    else
        disp(['COULD NOT COPY ' server_filename]);
    end
end

disp('Saved phy clustering to server');

% Delete local phy files
delete_confirm = strcmp(input('Delete local phy (y/n)? ','s'),'y');
if delete_confirm
    rmdir(local_phy_dir,'s');
    mkdir(local_phy_dir);
    disp('Cleared local phy data');
end


