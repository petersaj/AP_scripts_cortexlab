function AP_run_kilosort(data_filename,sample_rate)
% AP_run_kilosort(data_filename,sample_rate)
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
ops.chanMap = 'C:\Users\Andrew\Documents\CarandiniHarrisLab\data\kilosort_temp\channel_maps\forPRBimecP3opt3.mat';
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

% save matlab and phy results
ks_results_dir = ['C:\Users\Andrew\Documents\CarandiniHarrisLab\data\kilosort_temp\kilosort_results'];
mkdir(ks_results_dir);

save([ks_results_dir filesep 'rez.mat'], 'rez', '-v7.3');
rezToPhy(rez, ks_results_dir);


%% Copy output files to original path (should be basket)

disp('Copying sorted data to original path...')

disp('MAKE THIS')


%% Delete temporary data on local SSD

disp('Deleting temporary files...')

disp('MAKE THIS')
%delete(ops.fproc);
disp('Done.')




