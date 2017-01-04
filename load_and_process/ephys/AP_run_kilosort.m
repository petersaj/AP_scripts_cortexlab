function AP_run_kilosort(data_filename,input_board,sample_rate)
% AP_run_kilosort(data_filename,input_board,sample_rate)
%
% data_filename = .dat flat binary file of all channels together
% 
% New version (old in _old)
% Runs kilosort (modified from master_file_example_MOVEME)

%% Assume the data is already copied locally for now
[data_path,data_file,data_ext] = fileparts(data_filename);

%% Run Kilosort

% Set path
data_path = ['C:\temp'];

% Load channel map
ops.chanMap = 'C:\Users\Andrew\Documents\CarandiniHarrisLab\data\kilosort_temp\forPRBimecP3opt3.mat';
load(ops.chanMap);

% Run config script to get options
AP_kilosort_config_IMEC_P3O2

tic; % start timer
%
if ops.GPU     
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end

if strcmp(ops.datatype , 'openEphys')
   ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
end
%
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)

% AutoMerge. rez2Phy will use for clusters the new 5th column of st3 if you run this)
%     rez = merge_posthoc2(rez);

% save matlab results file
save(fullfile(ops.root,  'rez.mat'), 'rez', '-v7.3');

% remove temporary file
delete(ops.fproc);

%% Convert kilosort output to npy

fprintf('Time %3.0fs. Done. saving results... \n', toc);

phy_dir = [local_path filesep 'phy'];
mkdir(phy_dir);
rezToPhy(rez, phy_dir);

fprintf('Time %3.0fs. Done. \n', toc);

%% Copy output files to original path (should be basket)

disp('Copying sorted data to original path...')

% phy folder
local_phy_path = [local_path filesep 'phy'];
copyfile(local_phy_path,data_path);

% ks_results file
ks_results_filename = [data_file '_ks_results.mat'];
local_ks_results_filename = [local_path filesep ks_results_filename];
copyfile(local_ks_results_filename,[data_path filesep ks_results_filename]);

%% Delete temporary data on local SSD

disp('Deleting temporary files...')

delete(local_data_filename);
rmdir(local_phy_path,'s');
delete(local_ks_results_filename);

disp('Done.')




