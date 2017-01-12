function AP_run_kilosort(data_filename,sample_rate)
% AP_run_kilosort(data_filename,sample_rate)
%
% data_filename = .dat flat binary file of all channels together
% 
% New version (old version in _old)
% Runs kilosort (modified from master_file_example_MOVEME)

%% Assume the data is already copied locally for now
[data_path,data_file,data_ext] = fileparts(data_filename);

%% Run Kilosort

% Load channel map
ops.chanMap = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\kilosort_channelmaps\forPRBimecP3opt3.mat';
load(ops.chanMap);

% Run config script to get options
AP_kilosort_config_IMEC_P3O2

tic;
if ops.GPU     
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end

if strcmp(ops.datatype , 'openEphys')
   ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
end

% Run kilosort
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)

% AutoMerge. rez2Phy will use for clusters the new 5th column of st3 if you run this)
%     rez = merge_posthoc2(rez);

% save matlab and phy results
ks_results_dir = [data_path filesep 'results'];
mkdir(ks_results_dir);

save([ks_results_dir filesep 'rez.mat'], 'rez', '-v7.3');
rezToPhy(rez, ks_results_dir);


%% Copy output files to original path (should be basket)

disp('Copying sorted data to original path...')

disp('MAKE THIS')


%% Delete temporary data on local SSD

disp('Deleting temporary files...')

delete(ops.fproc);
disp('ALSO DELETE THE REST OF THE LOCAL RESULTS AFTER THEY''RE MOVED')

disp('Done.')




