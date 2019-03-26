function ephys_sync_samples = AP_sync_from_channel(raw_data_filename,processed_data_path)
% Get sync pulse times from electrophysiology recording
% (sync trace is the last channel recorded)
%
% INPUTS: 
% raw_data_filename = filename of the .dat recording file
% processed_data_path = general processed path, appends animal/day
%
% OUTPUTS: 
% saves sync channel in analysis folder if not already saved
% ephys_sync_samples = pulse start samples in trace

% Load meta file to get number of channels
[raw_data_path,raw_data_file,raw_data_ext] = fileparts(raw_data_filename);
ephys_header_filename = [raw_data_path filesep raw_data_file '.meta'];
ephys_header_id = fopen(ephys_header_filename);
ephys_header = textscan(ephys_header_id,'%s %s', 'delimiter',{' ='});
fclose(ephys_header_id);

nChans_idx = strcmp(ephys_header{1},'nChans');
nChans = str2num(ephys_header{2}{nChans_idx});

% Get last two paths of raw data (assume animal/day)
data_paths = strread(raw_data_path,'%s','delimiter','\\');
processed_data_path_expt = [processed_data_path filesep  ...
    data_paths{end-1} filesep data_paths{end}];

if ~exist(processed_data_path_expt,'dir');
    mkdir(processed_data_path_expt)
end

% Pull out and save sync channel if not already done
ephys_sync_filename = [processed_data_path_expt filesep raw_data_file ...
    '_ch' num2str(nChans) ,'.mat'];
if ~exist(ephys_sync_filename,'file');
    
    [raw_data_path,raw_data_file,raw_data_ext] = fileparts(raw_data_filename);
    
    disp('Sync trace not yet pulled. Copying raw data to local SSD...')
    local_ssd_path = 'C:\Users\Andy\Documents\CarandiniHarrisLab\data';
    temp_data_filename = [local_ssd_path filesep raw_data_file raw_data_ext];
    copyfile(raw_data_filename,temp_data_filename);
    
    disp('Pulling and saving sync channel...')
    AP_save_sync_chan(temp_data_filename,ephys_sync_filename,nChans,nChans);
    
    % Delete temporary local data
    delete(temp_data_filename);
    
end

% Get sync pulse samples from sync trace:
% 1) load pulled/saved trace
ephys_sync = load(ephys_sync_filename);
% 2) abs and rescale trace (arbitrary offset and scale)
ephys_sync = mat2gray(abs(ephys_sync.dat));
ephys_sync_samples = find(ephys_sync(1:end-1) <= 0.5 & ephys_sync(2:end) > 0.5);




















