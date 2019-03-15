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

%%% Run kilosort 2

% Set directory to save results
ks_results_dir = [save_path filesep 'results'];
mkdir(ks_results_dir);

% preprocess data to create temp_wh.dat
rez = preprocessDataSub(ops);

% time-reordering as a function of drift
rez = clusterSingleBatches(rez);
save(fullfile(ks_results_dir, 'rez.mat'), 'rez', '-v7.3');

% main tracking and template matching algorithm
rez = learnAndSolve8b(rez);

% final merges
rez = find_merges(rez, 1);

% final splits by SVD
rez = splitAllClusters(rez, 1);

% final splits by amplitudes
rez = splitAllClusters(rez, 0);

% decide on cutoff
rez = set_cutoff(rez);

fprintf('found %d good units \n', sum(rez.good>0))

% write to Phy
fprintf('Saving results to Phy  \n')
rezToPhy(rez, ks_results_dir);

%%% Save results 

% discard features in final rez file (too slow to save)
rez.cProj = [];
rez.cProjPC = [];

% save final results as rez2
fprintf('Saving final results in rez2  \n')
save([ks_results_dir filesep 'rez2.mat'], 'rez', '-v7.3');


%% Delete temporary data on local SSD

disp('Deleting temporary files...')

delete(ops.fproc);

disp('Done.')




