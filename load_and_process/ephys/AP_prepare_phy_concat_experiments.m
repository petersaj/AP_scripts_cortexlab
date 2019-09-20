function AP_prepare_phy_concat_experiments(animal,day,site,staging)
% AP_prepare_phy_concat_experiments(animal,day,site,staging)
%
% Prepares data for Phy
% - clears local phy folder
% - sets up params file
% - copies kilosorted data from server to local
% - copies AND CONCATENATES AP band recordings 
% -- NOTE: if this is present on HDD, then copies from there, otherwise
% copies, CAR's, and concatenates from server
%
% staging - puts in staging folder (false by default)
%
%

%% Set defaults

if ~exist('site','var')
    site = [];
end

if ~exist('staging','var') || isempty(staging)
    staging = false;
end

%% Get parts of raw filename

% Get ephys and kilosorted directory
[ephys_path,ephys_exists] = AP_cortexlab_filename(animal,day,[],'ephys_dir',site);
[ks_path,ks_exists] = AP_cortexlab_filename(animal,day,[],'ephys',site);

if ~ephys_exists || ~ks_exists 
    error('No kilosorted ephys data')
end

% Get experiments (if turned off between)
ephys_exp_paths = dir([ephys_path filesep 'experiment*']);

% Get OE filenames (check for multiple experiments, do separately)
exp_rec_dirs = cellfun(@(x) [x filesep 'recording1'],{ephys_exp_paths.name},'uni',false');
ap_data_filename = cellfun(@(exp_rec_dir) ...
    [ephys_path filesep exp_rec_dir filesep 'continuous' filesep 'Neuropix-3a-100.0' filesep 'continuous.dat'], ...
    exp_rec_dirs,'uni',false);

%% Clear out local phy folder

if ~staging
    local_phy_dir = 'C:\data_temp\phy';
elseif staging
    local_phy_dir = 'C:\data_temp\phy_staging';
end

if ~exist(local_phy_dir)
    mkdir(local_phy_dir)
elseif length(dir(local_phy_dir)) > 2
    rmdir(local_phy_dir,'s');
    mkdir(local_phy_dir);
end

disp(['Using local phy directory: ' local_phy_dir]);

%% Copy kilosorted data to local

disp('Copying kilosorted data to local...')
tic
copyfile(ks_path,local_phy_dir);
toc
disp('Done.');

%% Copy raw to local

% Get the expected filename from params
header_path = [local_phy_dir filesep 'params.py'];
header_fid = fopen(header_path);
header_info = textscan(header_fid,'%s %s', 'delimiter',{' = '});
fclose(header_fid);
header = struct;
for i = 1:length(header_info{1})
    header.(header_info{1}{i}) = header_info{2}{i};
end

% Set local AP data name
local_ap_filename = [local_phy_dir filesep header.dat_path(2:end-1)];

% Copy AP data from all experiments to local HDD
% (set HDD path and clear out)
hdd_phy_path = 'E:\data_temp\phy';
rmdir(hdd_phy_path,'s');
mkdir(hdd_phy_path);

% Check whether CAR/concatenated data is on HDD
% (this is part of the pipeline from
% AP_preprocess_phase3_newOE_concat_experiments)
hdd_storage_path = ['E:\local_data\' animal filesep day];
hdd_storage_fn = [hdd_storage_path filesep header.dat_path(2:end-1)];

if exist(hdd_storage_fn,'file')
    % If the file is on the HDD, copy from there to SDD
    disp('Copying preprocessed AP file from HDD to SDD...')
    copyfile(hdd_storage_fn,local_ap_filename);
    
else
    % Otherwise, copy/CAR/concatenate from server
    % (set temp HDD filename)
    ap_temp_filename = cellfun(@(recording_name) ...
        [hdd_phy_path filesep animal '_' day  ...
        '_' 'ephys_apband_' recording_name '.dat'], ...
        {ephys_exp_paths.name},'uni',false);
    
    % (copy AP files to HDD)
    disp('Copying AP files from server to HDD...')
    tic
    disp(datestr(now,'yyyy-mm-dd HH:MM'));
    for curr_recording = 1:length(ephys_exp_paths)
        copyfile(ap_data_filename{curr_recording},ap_temp_filename{curr_recording});
        AP_print_progress_fraction(curr_recording,length(ephys_exp_paths));
    end
    disp(datestr(now,'yyyy-mm-dd HH:MM'));
    toc
    
    % (CAR/copy/concatenate AP files into SDD phy folder)
    % (CAR not necessary, but read anyway and self-contained so might as well)
    disp('Car/copy/concatenating AP files HDD to SDD...')
    tic
    disp(datestr(now,'yyyy-mm-dd HH:MM'));
    n_chan = 384;
    AP_applyCARtoDat(ap_temp_filename,n_chan,local_ap_filename);
    disp(datestr(now,'yyyy-mm-dd HH:MM'));
    toc

end

disp('Done.');

%% Make a sound to indicate finished

fs = 8192;
x = 1:fs*0.05;
y = [sin(5000*x/fs),sin(7500*x/fs),sin(10000*x/fs)]*0.2;
sound(y);






