function AP_prepare_phy(animal,day,site,car)
% AP_prepare_phy(animal,day,site,car)
%
% Clears whatever is currently in local phy folder
% Sets up params file
% Loads raw and kilosorted data to local file
% car - if true, does common-average referencing

%% Set defaults

if ~exist('site','var')
    site = [];
end

if ~exist('car','var') || isempty(car)
    car = false;
end

%% Get parts of raw filename

[ephys_path,ephys_exists] = AP_cortexlab_filename(animal,day,[],'ephys',site);
[ephys_raw_ap_filename,ephys_raw_ap_exists] = AP_cortexlab_filename(animal,day,[],'ephys_ap',site);

if ~ephys_exists || ~ephys_raw_ap_exists
    error('No ephys data')
end

%% TEMPORARY: replace directory kilosort with kilosort2
% (once everything is sorted, AP_cortexlab_filename should default to using
% kilosort2 if it's available and kilosort if it's not)

warning('Temporary: replacing kilosort with kilosort 2')
ephys_path = strrep(ephys_path,'kilosort','kilosort2');
if ~exist(ephys_path,'dir')
    error('No kilosort2 directory');
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
tic
copyfile(ephys_path,local_phy_dir);
toc
disp('Done.');

%% Copy raw to local
disp('Copying raw data to local...')

% get the expected filename from params
header_path = [local_phy_dir filesep 'params.py'];
header_fid = fopen(header_path);
header_info = textscan(header_fid,'%s %s', 'delimiter',{' = '});
fclose(header_fid);
header = struct;
for i = 1:length(header_info{1})
    header.(header_info{1}{i}) = header_info{2}{i};
end
    
tic
local_ap_filename = [local_phy_dir filesep header.dat_path(2:end-1)];
copyfile(ephys_raw_ap_filename,local_ap_filename);
toc

if car
    disp('Doing common-average referencing...')
    % If 'CAR' appears in the filename, get rid of it
    car_text = strfind(local_ap_filename,'_CAR');
    noncar_ap_filename = local_ap_filename;
    noncar_ap_filename(car_text:car_text+3) = [];
    movefile(local_ap_filename,noncar_ap_filename);

    % Subtract common median across AP-band channels (hardcode channels?)
    n_chan = 384;
    medianTrace = applyCARtoDat(noncar_ap_filename,n_chan);
end

disp('Done.');


