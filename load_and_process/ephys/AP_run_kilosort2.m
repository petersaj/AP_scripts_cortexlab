function AP_run_kilosort2(data_filename,sample_rate,save_path,t_range)
% AP_run_kilosort2(data_filename,sample_rate,save_path,t_range)
%
% data_filename = .dat flat binary file of all channels together
% save_path = where the output files are saved
% t_range = data to use/truncate if recording issue (default [0,inf])
% 
% New version (old version in _old)
% Runs kilosort (modified from master_file_example_MOVEME)

%% Assume the data is already copied locally for now
[data_path,data_file,data_ext] = fileparts(data_filename);

%% Run Kilosort

% Run config script to get options
AP_kilosort2_config_IMEC_P3O2

tic;
if ops.GPU     
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end

% (not used in kilosort2?)
% if strcmp(ops.datatype , 'openEphys')
%    ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
% end

% Run kilosort 2
rez = preprocessDataSub(ops);
rez = clusterSingleBatches(rez);
rez = learnAndSolve8b(rez);
rez = splitAllClusters(rez);

% Convert results to phy, save
ks_results_dir = [save_path filesep 'results'];
mkdir(ks_results_dir);

save([ks_results_dir filesep 'rez.mat'], 'rez', '-v7.3');
rezToPhy(rez, ks_results_dir);


%% Delete temporary data on local SSD

disp('Deleting temporary files...')

delete(ops.fproc);

disp('Done.')




