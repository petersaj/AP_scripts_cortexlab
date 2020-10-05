function AP_clean_dat(ap_filename,n_channels,ttl_path,output_filename)
% AP_clean_dat(ap_filename,n_channels,ttl_path,output_path)
%
% ap_filename - filename of AP band DAT file
% n_channels - number of channels (neuropixels: 384)
% ttl_path - path of TTL sync data
% output_filename - filename to write cleaned AP DAT file
%
% Remove artifacts from raw ephys AP band DAT file and save new file

% Memory map AP data
n_bytes = 2; % default is int16 = 2 bytes
ap_fileinfo = dir(ap_filename);
n_ap_samples = ap_fileinfo.bytes/n_bytes/n_channels;
ap_memmap = memmapfile(ap_filename,'Format',{'int16',[n_channels,n_ap_samples],'ap'});

% Get light onsets
sync_filename = [ttl_path filesep 'channel_states.npy'];
sync_timestamps_filename = [ttl_path filesep 'timestamps.npy'];

sync_data = double(readNPY(sync_filename));
sync_timestamps_int = double(readNPY(sync_timestamps_filename));

led_sync_idx = 3; % hardcoded: input channel for LED TTL
light_on_idx = sync_data == led_sync_idx;
light_off_idx = sync_data == -led_sync_idx;

% Set chunks to process
% (lght order resets between experiments:
% split chunks with experiment gap to not mix light colors)
chunk_size = 1000000;
ap_samples_split = unique([1:chunk_size:n_ap_samples,n_ap_samples]);

% (to exclude multiple experiments in one chunk to avoid light order
% switching: not needed anymore since not using average light artifact)
% ap_sample_rate = 30000; % always the same, just hardcode
% expt_gap_time = 2; % minimum time to define experiment gap
% expt_gap_samples = round(expt_gap_time)*ap_sample_rate;
% light_on_samples = sync_timestamps_int(light_on_idx);
% expt_gaps = light_on_samples(find(diff(light_on_samples) > expt_gap_samples)+1);
% expt_gaps_leeway = expt_gaps - round(ap_sample_rate*0.5);
% ap_samples_split = unique([ap_samples_split,expt_gaps_leeway']);

% Open file for writing cleaned DAT
fidOut = fopen(output_filename, 'w');

% Loop through chunks and remove artifacts
for curr_chunk = 1:length(ap_samples_split)-1
   
    % Pull data chunk
    use_samples = ap_samples_split(curr_chunk): ...
        ap_samples_split(curr_chunk+1);

    curr_ap_raw = ap_memmap.Data.ap(:,use_samples);

    % Get current light syncs
    curr_sync_idx = ...
        sync_timestamps_int >= use_samples(1) & ...
        sync_timestamps_int <= use_samples(end);
    
    light_on = sync_timestamps_int(light_on_idx & curr_sync_idx)- use_samples(1) + 1;
    light_off = sync_timestamps_int(light_off_idx & curr_sync_idx) - use_samples(1) + 1;
    
    % Subtract median within each channel
    curr_ap_medsub = curr_ap_raw - nanmedian(curr_ap_raw,2);
    
    % Subtract moving median across channels (local median reference)
    n_channels_movmed = 20;
    curr_ap_medsub_localref = curr_ap_medsub - ...
        movmedian(curr_ap_medsub,n_channels_movmed,1);
    
    % Zero-out samples around light on/off (often small artifacts)
    n_light_samples_remove = 5;
    zero_light_samples = reshape(sort([light_on;light_off] + ...
        [-n_light_samples_remove:n_light_samples_remove]),[],1);
    zero_light_samples(zero_light_samples < 1 | ...
        zero_light_samples > size(curr_ap_raw,2)) = [];
    
    curr_ap_medsub_localref_lightzero = curr_ap_medsub_localref;
    curr_ap_medsub_localref_lightzero(:,zero_light_samples) = 0;

    % Write cleaned DAT to new file
    fwrite(fidOut, curr_ap_medsub_localref_lightzero, 'int16');
    
%     % Unused tested method: remove median light artifact
%     % (barely does better than moving median subtraction, and with light
%     % on/off interpolation doesn't make any difference)
%     % Get current light syncs (light alternates) 
%     blue_on = light_on(1:2:end);
%     blue_off = light_off(1:2:end);
%     violet_on = light_on(2:2:end);
%     violet_off = light_off(2:2:end);
%     
%     % Get light on/off time average
%     light_on_mean = mean(light_off - light_on);
%     light_off_mean = mean(light_on(2:end) - light_off(1:end-1));
%     
%     light_surround_samples = [-round(light_off_mean/2):round(light_on_mean+(light_off_mean/2))];
%    
%     blue_on_pull_t = blue_on + light_surround_samples;
%     blue_on_pull_t(any(blue_on_pull_t < 0 | blue_on_pull_t > size(curr_ap_raw,2),2),:) = [];
%     violet_on_pull_t = violet_on + light_surround_samples;
%     violet_on_pull_t(any(violet_on_pull_t < 0 | violet_on_pull_t > size(curr_ap_raw,2),2),:) = [];
%     
%     ap_blue = reshape(curr_ap_raw(:,blue_on_pull_t'),n_channels,length(light_surround_samples),[]);
%     ap_violet = reshape(curr_ap_raw(:,violet_on_pull_t'),n_channels,length(light_surround_samples),[]);
%     
%     % subtract baseline
%     baseline_sample = -1;
%     ap_blue_lightsubtract = nanmedian(ap_blue - ap_blue(:,light_surround_samples == baseline_sample,:),3);
%     ap_violet_lightsubtract = nanmedian(ap_violet - ap_violet(:,light_surround_samples == baseline_sample,:),3);
%         
%     curr_ap_lightsub = curr_ap_raw;
%     curr_ap_lightsub(:,blue_on_pull_t') = curr_ap_lightsub(:,blue_on_pull_t') - repmat(ap_blue_lightsubtract,1,size(blue_on_pull_t,1));
%     curr_ap_lightsub(:,violet_on_pull_t') = curr_ap_lightsub(:,violet_on_pull_t') - repmat(ap_violet_lightsubtract,1,size(violet_on_pull_t,1));
%     
%     curr_ap_lightsub_medsub = curr_ap_lightsub - nanmedian(curr_ap_lightsub,2);
%     
%     n_channels_median = 10;
%     curr_ap_lightsub_medsub_car = curr_ap_lightsub_medsub - ...
%         movmedian(curr_ap_lightsub_medsub,n_channels_median,1);

    AP_print_progress_fraction(curr_chunk,length(ap_samples_split)-1);
end

% Close all open files
fclose('all');







