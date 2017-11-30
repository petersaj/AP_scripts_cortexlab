function AP_prepare_phy(animal,day,car)
% AP_prepare_phy(data_filename)
%
% Clears whatever is currently in local phy folder
% Sets up params file
% Loads raw and kilosorted data to local file
% TO DO:
% CAR - true does common average referencing

%% Get parts of raw filename

[ephys_path,ephys_exists] = AP_cortexlab_filename(animal,day,[],'ephys');

[ephys_raw_path,ephys_raw_exists] = AP_cortexlab_filename(animal,day,[],'ephysraw');
ap_filename = 'experiment1_100-0_0.dat';

if ~ephys_exists || ~ephys_raw_exists
    error('No ephys data')
end

%% Clear out local phy folder

local_phy_dir = 'C:\data_temp\phy';
if ~exist(local_phy_dir)
    mkdir(local_phy_dir)
elseif length(dir(local_phy_dir)) > 2
    rmdir(local_phy_dir,'s');
    mkdir(local_phy_dir);
end

%% Copy kilosorted data to local

disp('Copying kilosorted data to local...')
copyfile(ephys_path,local_phy_dir);
disp('Done.');

%% Copy raw to local
disp('Copying raw data to local (NOT DOING CAR AT THE MOMENT)...')

% get the expected filename from params
header_path = [local_phy_dir filesep 'params.py'];
header_fid = fopen(header_path);
header_info = textscan(header_fid,'%s %s', 'delimiter',{' = '});
fclose(header_fid);
header = struct;
for i = 1:length(header_info{1})
    header.(header_info{1}{i}) = header_info{2}{i};
end
    
copyfile([ephys_raw_path ap_filename],[local_phy_dir filesep header.dat_path(2:end-1)]);
disp('Done.');



































