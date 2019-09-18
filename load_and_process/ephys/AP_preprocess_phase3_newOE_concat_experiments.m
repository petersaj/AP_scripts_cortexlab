function AP_preprocess_phase3_newOE_concat_experiments(animal,day,t_range)
% AP_preprocess_phase3(animal,day,t_range)
%
% t_range = specify data time range to use
% (now using kilosort2, putting into a 'kilosort2' folder)
%
% Copied from AP_preprocess_phase3_newOE:
% doesn't accomodate old OE, done when started multiple experiments
% (stopped and started ephys recordings, which creates multiple experiment
% folders)

%% Get paths and filenames

[ephys_path,ephys_exists] = AP_cortexlab_filename(animal,day,[],'ephys_dir');

if ~ephys_exists
    error([animal ' ' day ': No ephys data found']);
end

save_paths = {[ephys_path filesep 'kilosort2']};
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
    
    % Set and make save path
    curr_save_path = save_paths{curr_site};
            if ~exist(curr_save_path,'dir')
            mkdir(curr_save_path)
        end
              
        % Get OE filenames (check for multiple experiments, do separately)
        exp_rec_dirs = cellfun(@(x) [x filesep 'recording1'],{ephys_exp_paths.name},'uni',false');
        ap_data_filename = cellfun(@(exp_rec_dir) ...
            [curr_data_path filesep exp_rec_dir filesep 'continuous' filesep 'Neuropix-3a-100.0' filesep 'continuous.dat'], ...
            exp_rec_dirs,'uni',false);
        sync_filename = cellfun(@(exp_rec_dir) ...
            [curr_data_path filesep filesep exp_rec_dir filesep 'events' filesep 'Neuropix-3a-100.0' filesep 'TTL_1' filesep 'channel_states.npy' ], ...
            exp_rec_dirs,'uni',false);
        sync_timestamps_filename = cellfun(@(exp_rec_dir) ...
            [curr_data_path filesep exp_rec_dir filesep 'events' filesep 'Neuropix-3a-100.0' filesep 'TTL_1' filesep 'timestamps.npy' ], ...
            exp_rec_dirs,'uni',false);
        messages_filename = cellfun(@(exp_rec_dir) ...
            [curr_data_path filesep exp_rec_dir filesep 'sync_messages.txt'], ...
            exp_rec_dirs,'uni',false);
        
        
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
                
        % DOING: need to get size of AP recordings and add to sync
        n_ap_samples = nan(length(ap_data_filename),1);
        for curr_recording = 1:length(ephys_exp_paths)
            d = dir(ap_data_filename{curr_recording});
            n_ap_samples(curr_recording) = d.bytes/n_channels/2;
        end
        
        % Get/save digital input event times
        for curr_recording = 1:length(ephys_exp_paths)
            sync_data = readNPY(sync_filename{curr_recording});
            sync_timestamps_int = readNPY(sync_timestamps_filename{curr_recording});
     
            % If first recording, initialize sync
            if curr_recording == 1
                sync_channels = unique(abs(sync_data));
                sync = struct('timestamps',cell(size(sync_channels)),'values',cell(size(sync_channels)));
            end
            
            % If not first recording, add offset 
            % (assumes 0-idx - no idea if that's true)
            if curr_recording > 1
                sync_timestamps_int = sync_timestamps_int + ...
                    n_ap_samples(curr_recording  - 1) - 1;
            end
  
            sync_timestamps = double(sync_timestamps_int)/ap_sample_rate;          
            
            for curr_sync = 1:length(sync_channels)
                sync_events = abs(sync_data) == (sync_channels(curr_sync));
                sync(curr_sync).timestamps = [sync(curr_sync).timestamps;sync_timestamps(sync_events)];
                sync(curr_sync).values = [sync(curr_sync).values;sign(sync_data(sync_events)) == 1];
            end            
        end
        
        sync_save_filename = [curr_save_path filesep 'sync.mat'];
        save(sync_save_filename,'sync');
        
        %% Run kilosort
        
        % Set up local directory and clear out
        ssd_kilosort_path = 'C:\data_temp\kilosort';
        hdd_kilosort_path = 'E:\data_temp\kilosort';
        
        % Clear out local kilosort directories
        rmdir(ssd_kilosort_path,'s');
        mkdir(ssd_kilosort_path);
        
        rmdir(hdd_kilosort_path,'s');
        mkdir(hdd_kilosort_path);
        
        % Clear out whatever's currently in phy (usually not enough room)
        local_phy_path = 'C:\data_temp\phy';
        rmdir(local_phy_path,'s');
        mkdir(local_phy_path);
        
        % Copy AP data locally
        disp('Copying AP data to local drive...')
        ap_temp_filename = cellfun(@(recording_name) ...
            [ssd_kilosort_path filesep animal '_' day  ...
                '_' 'ephys_apband_' recording_name '.dat'], ...
                {ephys_exp_paths.name},'uni',false);
        for curr_recording = 1:length(ephys_exp_paths)
            copyfile(ap_data_filename{curr_recording},ap_temp_filename{curr_recording});
            AP_print_progress_fraction(curr_recording,length(ephys_exp_paths));
        end
        disp('Done');
                     
        % Common average reference data from server to HDD
        % (input multiple filenames, concatenate into one file)
        ops.NchanTOT = n_channels;
        ap_car_filename = [hdd_kilosort_path filesep animal '_' day '_' 'ephys_apband_CAR.dat'];
        AP_applyCARtoDat(ap_temp_filename,ops.NchanTOT,ap_car_filename);
                
        % Delete the un-CAR'd data
        for curr_recording = 1:length(ephys_exp_paths)
            delete(ap_temp_filename{curr_recording});
        end
        
        % Kilosort 1 (old)
        %     AP_run_kilosort(ap_temp_car_filename,ap_sample_rate);
        
        % Kilosort 2
        % Set default time range to be [0,inf]
        % NOTE: should probably change this to be 5-10s cut on each end? Need
        % to be able to get total time in seconds from file
        if ~exist('t_range','var')
            t_range = [0,inf];
        end
        AP_run_kilosort2(ap_car_filename,ap_sample_rate,ssd_kilosort_path,t_range);
        
        %% Copy kilosort results to server
        
        disp('Copying sorted data to server...');
        ks_results_path = [ssd_kilosort_path filesep 'results'];
        copyfile(ks_results_path,curr_save_path);
        
        %% Copy kilosort results and raw data to phy folder for clustering
        % (only if clustering immediately afterwards)
        
        %     local_phy_path = 'C:\data_temp\phy';
        %
        %     % Clear out whatever's currently in phy
        %     rmdir(local_phy_path,'s');
        %     mkdir(local_phy_path);
        %
        %     % Copy the CAR'd data
        %     [~,ap_file,ap_ext] = fileparts(ap_temp_car_filename_hdd);
        %     movefile(ap_temp_car_filename_hdd,[local_phy_path filesep ap_file ap_ext])
        %
        %     % Copy the results
        %     movefile([ks_results_path filesep '*'],local_phy_path)
        
        %% Delete all temporarly local data
        rmdir(ssd_kilosort_path,'s');
        mkdir(ssd_kilosort_path);
        
    
end

disp('Done processing phase 3 data.');


