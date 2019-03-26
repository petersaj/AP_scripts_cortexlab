function AP_save_phy
% AP_save_phy
%
% Move clustering data from local phy directory to server

% Local phy directory
local_phy_dir = 'C:\data_temp\phy';

% Get animal and day from header
header_filename = [local_phy_dir filesep 'dat_params.txt'];

header_fid = fopen(header_filename);
header_info = textscan(header_fid,'%s %s', 'delimiter',{' = '});
fclose(header_fid);

raw_path_idx = strcmp(header_info{1},'raw_path');
raw_path = header_info{2}{raw_path_idx};

animal_day = regexp(raw_path,'\\Subjects\\(\w*)\\(\d\d\d\d-\d\d-\d\d)','tokens');
animal = animal_day{1}{1};
day = animal_day{1}{2};

% Get server kilosort directory
[ephys_path,ephys_path_exists] = AP_cortexlab_filename(animal,day,[],'ephys');

% TEMPORARY: replace directory kilosort with kilosort2
% (once everything is sorted, AP_cortexlab_filename should default to using
% kilosort2 if it's available and kilosort if it's not)
warning('Temporary: replacing kilosort with kilosort 2')
ephys_path = strrep(ephys_path,'kilosort','kilosort2');
if ~exist(ephys_path,'dir')
    error('No kilosort2 directory');
end

% User confirm animal and day overwrite
user_confirm = strcmp(input(['Overwrite phy sorting ' animal ' ' day ' (y/n)? '],'s'),'y');
if ~user_confirm
    disp('Not saving');
    return
end

% Set standard files to copy
phy_files = {'cluster_group.tsv','spike_clusters.npy'};

% User select save mean cluster waveforms
user_mean_waveforms = strcmp(input(['Save mean cluster waveforms (y/n)? '],'s'),'y');

% Make cluster mean waveforms if selected
% cluster x time x channel like templates, rows go in order of
% unique(spike_clusters)
if user_mean_waveforms
    
    n_waveforms = 200; % number of waveforms to average
    
    spike_times = readNPY([local_phy_dir filesep 'spike_times.npy']);
    spike_clusters = readNPY([local_phy_dir filesep 'spike_clusters.npy']);
    
    dat_dir = dir([local_phy_dir filesep '*.dat']);
    
    gwfparams.dataDir = local_phy_dir;
    gwfparams.fileName = [local_phy_dir filesep dat_dir.name];
    gwfparams.dataType = 'int16';
    gwfparams.nCh = 384;
    gwfparams.wfWin = [-40 80];
    gwfparams.nWf = n_waveforms;
    gwfparams.spikeTimes = spike_times;
    gwfparams.spikeClusters = spike_clusters;
    
    wf = getWaveForms(gwfparams);
    
    cluster_mean_waveforms = permute(wf.waveFormsMean,[1,3,2]);
    
    % (baseline subtract and normalize: do this on load)
    % cluster_mean_waveforms = (templates-templates(:,1,:))./abs(max(max(templates,[],2),[],3));
    
    save([local_phy_dir filesep 'cluster_mean_waveforms.mat'],'cluster_mean_waveforms');
    
    phy_files{end+1} = 'cluster_mean_waveforms.mat';
    
end

% Move files to server
for curr_file = 1:length(phy_files)
    local_filename = [local_phy_dir filesep phy_files{curr_file}];
    server_filename = [ephys_path filesep phy_files{curr_file}];
    copyfile(local_filename,server_filename);
    disp(['Overwrote ' server_filename]);
end

disp('Saved phy clustering to server');




