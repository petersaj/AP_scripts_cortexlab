function AP_preprocess_phase3(animal,day,t_range)
% AP_preprocess_phase3(animal,day,t_range)
%
% t_range = specify data time range to use
% (now using kilosort2, putting into a 'kilosort2' folder)

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
    
    curr_data_path = data_paths{curr_site};
    curr_save_path = save_paths{curr_site};
    
    if ~exist(curr_save_path,'dir')
        mkdir(curr_save_path)
    end
    
    % Switch file formats between old/new versions
    % (do this by looking for a folder called 'continuous' created
    % automatically by open ephys)
    if exist([curr_data_path filesep 'experiment1'],'dir')
        oe_file_structure = 'new';
    else
        oe_file_structure = 'old';
    end
    
    switch oe_file_structure
        
        case 'old'
            % Filenames are semi-hardcoded in open ephys convention
            ap_data_dir = dir([curr_data_path filesep 'experiment*_10*-0_0.dat']);
            lfp_data_dir = dir([curr_data_path filesep 'experiment*_10*-1_0.dat']);
            sync_dir = dir([curr_data_path filesep 'experiment*_all_channels_0.events']);
            messages_dir = dir([curr_data_path filesep 'experiment*_messages_0.events']);
            settings_dir = dir([curr_data_path filesep 'settings*.xml']);
            
            ap_data_filename = [curr_data_path filesep ap_data_dir.name];
            lfp_data_filename = [curr_data_path filesep lfp_data_dir.name];
            sync_filename = [curr_data_path filesep sync_dir.name];
            messages_filename = [curr_data_path filesep messages_dir.name];
            settings_filename = [curr_data_path filesep settings_dir.name];
            
        case 'new'
            % and then they totally changed everything (2018-07-27)
            exp_rec_dir = ['experiment1' filesep 'recording1'];
            ap_data_filename = [curr_data_path filesep exp_rec_dir filesep 'continuous' filesep 'Neuropix-3a-100.0' filesep 'continuous.dat'];
            lfp_data_filename = [curr_data_path filesep exp_rec_dir filesep 'continuous' filesep 'Neuropix-3a-100.1' filesep 'continuous.dat'];
            sync_filename = [curr_data_path filesep filesep exp_rec_dir filesep 'events' filesep 'Neuropix-3a-100.0' filesep 'TTL_1' filesep 'channel_states.npy' ];
            sync_timestamps_filename = [curr_data_path filesep exp_rec_dir filesep 'events' filesep 'Neuropix-3a-100.0' filesep 'TTL_1' filesep 'timestamps.npy' ];
            messages_filename = [curr_data_path filesep exp_rec_dir filesep 'sync_messages.txt'];
            settings_filename = [curr_data_path filesep exp_rec_dir filesep 'structure.oebin'];
            
    end
    
    %% Get and save recording parameters
    
    switch oe_file_structure
        
        case 'old'           
            % Get index of electrophysiology channels in recordings
            ephys_settings = xml2struct(settings_filename);
            
            % Get sample rate, gain, cutoff frequency (separate numbers and suffixes)
            ap_gain = textscan(ephys_settings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.EDITOR.NEUROPIXELS.Attributes.apGainValue,'%d%s');
            lfp_gain = textscan(ephys_settings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.EDITOR.NEUROPIXELS.Attributes.lfpGainValue,'%d%s');
            filter_cut = textscan(ephys_settings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.EDITOR.NEUROPIXELS.Attributes.filterCut,'%d%s');
            
        case 'new'
            % The gains and filter cuts aren't recorded anymore?!
            ap_gain = {500};
            lfp_gain = {125}; 
            filter_cut = {300};
    end
    
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
    
    switch oe_file_structure
        
        case 'old'         
            % Get experiment start time (these messages are saved in a super dumb way
            % that are impossible to parse generally, so this is messy)
            messages_id = fopen(messages_filename);
            messages_text = textscan(messages_id,'%*d %s %s', 'delimiter',{': '});
            fclose(messages_id);
            
            start_time_idx = strcmp(messages_text{1},'start time');
            start_time = str2num(messages_text{2}{start_time_idx}(1:strfind(messages_text{2}{start_time_idx},'@')-1));
            start_time_freq = str2num(messages_text{2}{start_time_idx}(strfind(messages_text{2}{start_time_idx},'@')+1: ...
                strfind(messages_text{2}{start_time_idx},'Hz')-1));
            start_time_sec = start_time/start_time_freq;
            
            % Get/save digital input event times,
            [sync_data, sync_timestamps, sync_info] = load_open_ephys_data_faster(sync_filename);
            sync_channels = unique(sync_data);
            sync = struct('timestamps',cell(size(sync_channels)),'values',cell(size(sync_channels)));
            for curr_sync = 1:length(sync_channels)
                sync_events = sync_data == (sync_channels(curr_sync));
                sync(curr_sync).timestamps = sync_timestamps(sync_events);
                sync(curr_sync).values = logical(sync_info.eventId(sync_events));
                
                % correct for experiment start time (not always necessary??)
                % as far as I can tell it's random whether this is needed or not: if
                % it's not then you get negative numbers at first, so maybe check for
                % those and then it can be automated? it's not a good sign that it's
                % variable though... I should probably just switch to spikeglx
                if sync(curr_sync).timestamps(1) - start_time_sec > 0
                    sync(curr_sync).timestamps = sync(curr_sync).timestamps - start_time_sec;
                end
                
            end
            
            sync_save_filename = [curr_save_path filesep 'sync.mat'];
            save(sync_save_filename,'sync');
            
        case 'new'
            
            % This way of getting start times is garbage too, assume AP
            % band is listed third (1 = software, 2 = AP, 3 = LFP)
            messages_id = fopen(messages_filename);
            messages_text = textscan(messages_id,'%*s %d@%dHz','delimiter',{'time: '});
            fclose(messages_id);
            
            start_time_sec = double(messages_text{1}(2)/messages_text{2}(2));
            
            % Get/save digital input event times,
            sync_data = readNPY(sync_filename);
            sync_timestamps_int = readNPY(sync_timestamps_filename);
            sync_timestamps = double(sync_timestamps_int)/ap_sample_rate;
            
            sync_channels = unique(abs(sync_data));
            sync = struct('timestamps',cell(size(sync_channels)),'values',cell(size(sync_channels)));
            for curr_sync = 1:length(sync_channels)
                sync_events = abs(sync_data) == (sync_channels(curr_sync));
                
                % Subtract system start time from sync
                % (Kilosort starts spikes at 0, open ephys sometimes doesn't
                % start sync time at 0 if record direct from preview)
                sync(curr_sync).timestamps = sync_timestamps(sync_events) - start_time_sec;
                sync(curr_sync).values = sign(sync_data(sync_events)) == 1;               
            end
            
            sync_save_filename = [curr_save_path filesep 'sync.mat'];
            save(sync_save_filename,'sync');
            
    end
    
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
    ap_temp_filename = [ssd_kilosort_path filesep animal '_' day  '_' 'ephys_apband.dat'];
    copyfile(ap_data_filename,ap_temp_filename);
    disp('Done');
    
    % Common average reference data from server to HDD
    ops.NchanTOT = n_channels;
    medianTrace = applyCARtoDat(ap_temp_filename,ops.NchanTOT,hdd_kilosort_path);         
    [~,ap_car_name,ap_car_ext] = fileparts(ap_temp_filename);
    ap_temp_car_filename_hdd = [hdd_kilosort_path filesep ap_car_name '_CAR' ap_car_ext];
    
    % Delete the un-CAR'd data
    delete(ap_temp_filename);
    
    % Kilosort 1 (old)
%     AP_run_kilosort(ap_temp_car_filename,ap_sample_rate);

    % Kilosort 2
    % Set default time range to be [0,inf]
    % NOTE: should probably change this to be 5-10s cut on each end? Need
    % to be able to get total time in seconds from file
    if ~exist('t_range','var')
        t_range = [0,inf];
    end
    AP_run_kilosort2(ap_temp_car_filename_hdd,ap_sample_rate,ssd_kilosort_path,t_range);
    
    %% Copy kilosort results to server
    
    disp('Copying sorted data to server...');
    ks_results_path = [ssd_kilosort_path filesep 'results'];
    copyfile(ks_results_path,curr_save_path);
    
    %% Copy kilosort results and raw data to phy folder for clustering
    % (only if clustering immediately afterwards)
    
    local_phy_path = 'C:\data_temp\phy';
    
    % Clear out whatever's currently in phy
    rmdir(local_phy_path,'s');
    mkdir(local_phy_path);
    
    % Copy the CAR'd data
    [~,ap_file,ap_ext] = fileparts(ap_temp_car_filename_hdd);
    movefile(ap_temp_car_filename_hdd,[local_phy_path filesep ap_file ap_ext])
    
    % Copy the results
    movefile([ks_results_path filesep '*'],local_phy_path)
    
    %% Delete all temporarly local data
    rmdir(ssd_kilosort_path,'s');
    mkdir(ssd_kilosort_path);
    
end

disp('Done processing phase 3 data.');


