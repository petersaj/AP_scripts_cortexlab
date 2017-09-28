function AP_preprocess_phase3(animal,day)
% AP_preprocess_phase3(animal,day)

%% Get paths and filenames

data_paths =  ...
    {['\\zserver.cortexlab.net\Data\Subjects\' animal filesep day '\ephys']};
save_paths =  ...
    {['\\basket.cortexlab.net\data\ajpeters\' animal filesep day '\ephys']};

% Check for multiple sites (based on the data path containing only folders)
data_path_dir = dir(data_paths{1});
if all([data_path_dir(3:end).isdir])
    data_paths = cellfun(@(x) [data_paths{1} filesep x],{data_path_dir(3:end).name},'uni',false);
    save_paths = cellfun(@(x) [save_paths{1} filesep x],{data_path_dir(3:end).name},'uni',false);
end

for curr_site = 1:length(data_paths)
    
    curr_data_path = data_paths{curr_site};
    curr_save_path = save_paths{curr_site};
    
    if ~exist(curr_save_path,'dir')
        mkdir(curr_save_path)
    end
    
    % Filenames are semi-hardcoded in open ephys convention
    ap_data_dir = dir([curr_data_path filesep 'experiment*_100-0_0.dat']);
    lfp_data_dir = dir([curr_data_path filesep 'experiment*_100-1_0.dat']);
    sync_dir = dir([curr_data_path filesep 'experiment*_all_channels_0.events']);
    messages_dir = dir([curr_data_path filesep 'experiment*_messages_0.events']);
    settings_dir = dir([curr_data_path filesep 'settings*.xml']);
    
    ap_data_filename = [curr_data_path filesep ap_data_dir.name];
    lfp_data_filename = [curr_data_path filesep lfp_data_dir.name];
    sync_filename = [curr_data_path filesep sync_dir.name];
    messages_filename = [curr_data_path filesep messages_dir.name];
    settings_filename = [curr_data_path filesep settings_dir.name];
    
    %% Get and save recording parameters
    
    % Get index of electrophysiology channels in recordings
    ephys_settings = xml2struct(settings_filename);
    
    % Get sample rate, gain, cutoff frequency (separate numbers and suffixes)
    ap_gain = textscan(ephys_settings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.EDITOR.NEUROPIXELS.Attributes.apGainValue,'%d%s');
    lfp_gain = textscan(ephys_settings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.EDITOR.NEUROPIXELS.Attributes.lfpGainValue,'%d%s');
    filter_cut = textscan(ephys_settings.SETTINGS.SIGNALCHAIN.PROCESSOR{1}.EDITOR.NEUROPIXELS.Attributes.filterCut,'%d%s');
    
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
        sync(curr_sync).timestamps = sync(curr_sync).timestamps - start_time_sec;
    end
    
    sync_save_filename = [curr_save_path filesep 'sync.mat'];
    save(sync_save_filename,'sync');
    
    %% Run kilosort
    
    % Set up local directory and clear out
    local_kilosort_path = 'C:\data_temp\kilosort';
    rmdir(local_kilosort_path,'s');
    mkdir(local_kilosort_path);
    
    % Clear out whatever's currently in phy (usually not enough room)
    local_phy_path = 'C:\data_temp\phy';
    rmdir(local_phy_path,'s');
    mkdir(local_phy_path);
    
    % Copy data locally
    disp('Copying data to local drive...')
    ap_temp_filename = [local_kilosort_path filesep animal '_' day  '_' 'ephys_apband.dat'];
    if ~exist(local_kilosort_path,'dir')
        mkdir(local_kilosort_path)
    end
    copyfile(ap_data_filename,ap_temp_filename);
    disp('Done');
    
    % Subtract common median across AP-band channels (hardcode channels?)
    ops.NchanTOT = 384;
    medianTrace = applyCARtoDat(ap_temp_filename, ops.NchanTOT);
    ap_temp_car_filename = [ap_temp_filename(1:end-4) '_CAR.dat'];
    
    % Get rid of the original non-CAR (usually not enough disk space)
    delete(ap_temp_filename);
    
    % Run kilosort on CAR data
    AP_run_kilosort(ap_temp_car_filename,ap_sample_rate);
    
    
    %% Copy kilosort results to basket
    
    disp('Copying sorted data to basket...');
    ks_results_path = [local_kilosort_path filesep 'results'];
    copyfile(ks_results_path,curr_save_path);
    
    %% Copy kilosort results and raw data to phy folder for clustering
    
    local_phy_path = 'C:\data_temp\phy';
    
    % Clear out whatever's currently in phy
    rmdir(local_phy_path,'s');
    mkdir(local_phy_path);
    
    % Copy the CAR'd data
    [~,ap_file,ap_ext] = fileparts(ap_temp_car_filename);
    movefile(ap_temp_car_filename,[local_phy_path filesep ap_file ap_ext])
    
    % Copy the results
    movefile([ks_results_path filesep '*'],local_phy_path)
    
    %% Delete all temporarly local data
    rmdir(local_kilosort_path,'s');
    mkdir(local_kilosort_path);
    
end

disp('Done processing phase 3 data.');


