function AP_preprocess_neuropixels(animal,day)
% AP_preprocess_neuropixels(animal,day)
%
% Spike sorts with Pykilosort (into 'pykilosort' folder)
%
% Currently assumes Neuropixels Phase 3A recorded with Open Ephys

%% Get paths and filenames

[ephys_path,ephys_exists] = AP_cortexlab_filename(animal,day,[],'ephys_dir');

if ~ephys_exists
    error([animal ' ' day ': No ephys data found']);
end

save_paths = {[ephys_path filesep 'pykilosort']};
data_paths = {ephys_path};

% Check for multiple sites (assume sites are marked as site#)
data_path_dir = dir([data_paths{1} filesep 'site*']);
if ~isempty(data_path_dir)
    data_paths = cellfun(@(x) [data_paths{1} filesep x],{data_path_dir.name},'uni',false);
    save_paths = cellfun(@(x) [save_paths{1} filesep x],{data_path_dir.name},'uni',false);
end

for curr_site = 1:length(data_paths)  
    
    % Get experiments (if turned off between)
    curr_data_path = data_paths{curr_site};   
    ephys_exp_paths = dir([curr_data_path filesep 'experiment*']);
       
    for curr_exp = 1:length(ephys_exp_paths)
        
        % Update save path with experiment (only if more than one, legacy)
        if length(ephys_exp_paths) == 1
            curr_save_path = save_paths{curr_site};
        elseif length(ephys_exp_paths) > 1
            curr_save_path = [save_paths{curr_site} filesep ephys_exp_paths(curr_exp).name];
        end
        
        % Make save path
        if ~exist(curr_save_path,'dir')
            mkdir(curr_save_path)
        end
              
        % Get OE filenames (check for multiple experiments, do separately)
        exp_rec_dir = [ephys_exp_paths(curr_exp).name filesep 'recording1'];
        ap_data_filename = [curr_data_path filesep exp_rec_dir filesep 'continuous' filesep 'Neuropix-3a-100.0' filesep 'continuous.dat'];
        lfp_data_filename = [curr_data_path filesep exp_rec_dir filesep 'continuous' filesep 'Neuropix-3a-100.1' filesep 'continuous.dat'];
        sync_filename = [curr_data_path filesep filesep exp_rec_dir filesep 'events' filesep 'Neuropix-3a-100.0' filesep 'TTL_1' filesep 'channel_states.npy' ];
        sync_timestamps_filename = [curr_data_path filesep exp_rec_dir filesep 'events' filesep 'Neuropix-3a-100.0' filesep 'TTL_1' filesep 'timestamps.npy' ];
        messages_filename = [curr_data_path filesep exp_rec_dir filesep 'sync_messages.txt'];
        settings_filename = [curr_data_path filesep exp_rec_dir filesep 'structure.oebin'];
        
        
        %% Get and save recording parameters
        
        % The gains and filter cuts aren't recorded anymore?!
        ap_gain = {500};
        lfp_gain = {125};
        filter_cut = {300};
        
        % (0.195x for int16 to uV? how's this change with gain, just another x?)
        
        % Hard-coded parameters
        n_channels = 384;
        ap_sample_rate = 30000;
        lfp_sample_rate = 2500;
        
        params = {'raw_path',['''' curr_data_path '''']; ...
            'n_channels',num2str(n_channels); ...
            'ap_sample_rate',num2str(ap_sample_rate); ... % this should be 30000 AP, 2500 LFP
            'lfp_sample_rate',num2str(lfp_sample_rate);
            'ap_gain',num2str(ap_gain{1}); ...
            'lfp_gain',num2str(lfp_gain{1})
            'filter_cutoff',num2str(filter_cut{1})};
        
        param_filename = [curr_save_path filesep 'dat_params.txt'];
        
        formatSpec = '%s = %s\r\n';
        fid = fopen(param_filename,'w');
        for curr_param = 1:size(params,1)
            fprintf(fid,formatSpec,params{curr_param,:});
        end
        fclose(fid);
        
        
        %% Get/save digital input events
        
        % Get processor start information 
        % (assume 1 = software, 2 = AP, 3 = LFP)
        % This way of getting start times is garbage too, assume AP
        messages_id = fopen(messages_filename);
        messages_text = textscan(messages_id,'%*s %d@%dHz','delimiter',{'time: '});
        fclose(messages_id);
        
        % Processor start time: to be subtracted from sync so start = 0
        % (this is usually 0 but not if hit record directly from preview,
        % so either does nothing or corrects offset)
        start_time_sec = double(messages_text{1}(2)/messages_text{2}(2));
        
        % Get/save digital input event times,
        sync_data = readNPY(sync_filename);
        sync_timestamps_int = readNPY(sync_timestamps_filename);
        sync_timestamps = double(sync_timestamps_int)/ap_sample_rate;
        
        sync_channels = unique(abs(sync_data));
        sync = struct('timestamps',cell(size(sync_channels)),'values',cell(size(sync_channels)));
        for curr_sync = 1:length(sync_channels)
            sync_events = abs(sync_data) == (sync_channels(curr_sync));
            sync(curr_sync).timestamps = sync_timestamps(sync_events) - start_time_sec;
            sync(curr_sync).values = sign(sync_data(sync_events)) == 1;
        end
        
        sync_save_filename = [curr_save_path filesep 'sync.mat'];
        save(sync_save_filename,'sync');
        
        %% Run kilosort
        
        % Set up local directory and clear out
        ssd_kilosort_path = 'G:\data_temp\kilosort';
        
        % Clear out local kilosort directories
        if exist(ssd_kilosort_path,'dir')
            rmdir(ssd_kilosort_path,'s')
        end
        mkdir(ssd_kilosort_path);
        
        % Copy AP-band data locally
        disp('Copying AP data to local drive...')
        apband_local_filename = [ssd_kilosort_path filesep animal '_' day  '_' 'ephys_apband.dat'];
        copyfile(ap_data_filename,apband_local_filename);
        disp('Done');
        
        % Set up python
        % (set pykilosort python environment)
        py_version = pyenv('Version','C:\Users\Andrew\Anaconda3\envs\pyks2\pythonw.exe');
        % (add pykilosort environment paths to windows system path)
        pre_pykilosort_syspath = getenv('PATH');
        py_env_paths = {
            fullfile(char(py_version.Home),'Library','bin'); ...
            fullfile(char(py_version.Home),'Scripts')};
        run_pykilosort_syspath = strjoin( ...
            unique(cat(1,py_env_paths, ...
            split(pre_pykilosort_syspath,pathsep)),'stable'),pathsep);
        setenv('PATH',run_pykilosort_syspath);

        % Run pykilosort
        % (directly on raw data - pykilosort does local common average
        pykilosort_output_path = fullfile(ssd_kilosort_path,'pykilosort');
        pyrunfile('AP_run_pykilosort.py', ...
            data_filename = apband_local_filename, ...
            pykilosort_output_path = pykilosort_output_path);

        % Revert system paths to pre-pykilosort
        % (just in case alternate python environments used elsewhere)
        setenv('PATH',pre_pykilosort_syspath);

        
        %% Copy kilosort results to server
                
        disp('Copying sorted data to server...');

        kilosort_results_path = fullfile(pykilosort_output_path,'output');
        copyfile(pykilosort_results_path,curr_save_path);
        
        %% Copy kilosort results and raw data to phy folder for clustering
        
        local_phy_path = ['G:\data_temp\phy_staging' filesep animal '_' day];
        mkdir(local_phy_path)
        
        % Move the cleaned data into the phy directory
        [~,ap_file,ap_ext] = fileparts(apband_local_filename);
        movefile(apband_local_filename,fullfile(local_phy_path,[ap_file ap_ext]))
        
        % Move the results into the phy directory
        movefile(fullfile(kilosort_results_path,'*'),local_phy_path)
        
        % Copy the dat_params file into the phy directory
        copyfile(param_filename,fullfile(local_phy_path,'dat_params.txt'));
        
        %% Delete all temporary local data
        rmdir(ssd_kilosort_path,'s');
        mkdir(ssd_kilosort_path);
        
    end
    
end

fprintf('Done preprocessing Neuropixels: %s %s',animal,day);


