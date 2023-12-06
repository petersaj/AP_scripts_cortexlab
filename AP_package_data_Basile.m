
%% Trained data

animals = {'AP100','AP101','AP104','AP105','AP106'};

for curr_animal = animals

    animal = cell2mat(curr_animal);

    external_drive_path = 'E:\Cortexlab_data';

    %% Find data on external hard drive

    % Initialize pathname, add to it with each server location
    days_combined = {};
    days_pathnames_combined = {};

    expInfo_path = fullfile(external_drive_path,animal);
    expInfo_dir = dir(expInfo_path);
    day_paths = cellfun(@(x) ~isempty(regexp(x,'\d\d\d\d-\d\d-\d\d')),{expInfo_dir.name}) &...
        [expInfo_dir.isdir];
    curr_days = {expInfo_dir(day_paths).name};
    curr_days_pathname = cellfun(@(x) [expInfo_path filesep x],curr_days,'uni',false);

    days_combined = [days_combined,curr_days];
    days_pathnames_combined = [days_pathnames_combined,curr_days_pathname];

    [days,unique_day_idx] = unique(days_combined);
    days_pathnames = days_pathnames_combined(unique_day_idx);

    ephys_recorded = false(size(days));
    for curr_day = 1:length(days)
        if exist(fullfile(external_drive_path,animal,days{curr_day},'ephys'),'dir')
            ephys_recorded(curr_day) = true;
        end
    end

    ephys_days = days(ephys_recorded);

    %% Get spike times and templates

    for curr_day_idx = 1:length(ephys_days)

        day = ephys_days{curr_day_idx};

        curr_day_path = fullfile(external_drive_path,animal,day);
        ephys_path = fullfile(curr_day_path,'ephys','kilosort2');

        %%%%%%%%% Load data (copied from AP_load_experiment)

        % Read header information
        header_path = [ephys_path filesep 'dat_params.txt'];
        header_fid = fopen(header_path);
        header_info = textscan(header_fid,'%s %s', 'delimiter',{' = '});
        fclose(header_fid);

        header = struct;
        for i = 1:length(header_info{1})
            header.(header_info{1}{i}) = header_info{2}{i};
        end

        % Load spike data
        if isfield(header,'sample_rate')
            ephys_sample_rate = str2num(header.sample_rate);
        elseif isfield(header,'ap_sample_rate')
            ephys_sample_rate = str2num(header.ap_sample_rate);
        end
        spike_times = double(readNPY([ephys_path filesep 'spike_times.npy']))./ephys_sample_rate;
        spike_templates_0idx = readNPY([ephys_path filesep 'spike_templates.npy']);
        templates_whitened = readNPY([ephys_path filesep 'templates.npy']);
        channel_positions = readNPY([ephys_path filesep 'channel_positions.npy']);
        channel_map = readNPY([ephys_path filesep 'channel_map.npy']);
        winv = readNPY([ephys_path filesep 'whitening_mat_inv.npy']);
        template_amplitudes = readNPY([ephys_path filesep 'amplitudes.npy']);

        % Default channel map/positions are from end: make from surface
        % (hardcode this: kilosort2 drops channels)
        max_depth = 3840;
        channel_positions(:,2) = max_depth - channel_positions(:,2);

        % Unwhiten templates
        templates = zeros(size(templates_whitened));
        for t = 1:size(templates_whitened,1)
            templates(t,:,:) = squeeze(templates_whitened(t,:,:))*winv;
        end

        % Get the waveform of all templates (channel with largest amplitude)
        [~,max_site] = max(max(abs(templates),[],2),[],3);
        templates_max = nan(size(templates,1),size(templates,2));
        for curr_template = 1:size(templates,1)
            templates_max(curr_template,:) = ...
                templates(curr_template,:,max_site(curr_template));
        end
        waveforms = templates_max;

        % Get depth of each template
        % (get min-max range for each channel)
        template_chan_amp = squeeze(range(templates,2));
        % (zero-out low amplitude channels)
        template_chan_amp_thresh = max(template_chan_amp,[],2)*0.5;
        template_chan_amp_overthresh = template_chan_amp.*(template_chan_amp >= template_chan_amp_thresh);
        % (get center-of-mass on thresholded channel amplitudes)
        template_depths = sum(template_chan_amp_overthresh.*channel_positions(:,2)',2)./sum(template_chan_amp_overthresh,2);

        % Get the depth of each spike (templates are zero-indexed)
        spike_depths = template_depths(spike_templates_0idx+1);

        % Get trough-to-peak time for each template
        templates_max_signfix = bsxfun(@times,templates_max, ...
            sign(abs(min(templates_max,[],2)) - abs(max(templates_max,[],2))));

        [~,waveform_trough] = min(templates_max,[],2);
        [~,waveform_peak_rel] = arrayfun(@(x) ...
            max(templates_max(x,waveform_trough(x):end),[],2), ...
            transpose(1:size(templates_max,1)));
        waveform_peak = waveform_peak_rel + waveform_trough;

        templateDuration = waveform_peak - waveform_trough;
        templateDuration_us = (templateDuration/ephys_sample_rate)*1e6;

        % Keep all units
        good_templates_idx = unique(spike_templates_0idx);
        good_templates = ismember(0:size(templates,1)-1,good_templates_idx);

        % Throw out all non-good template data
        templates = templates(good_templates,:,:);
        template_depths = template_depths(good_templates);
        waveforms = waveforms(good_templates,:);
        templateDuration = templateDuration(good_templates);
        templateDuration_us = templateDuration_us(good_templates);

        % Throw out all non-good spike data
        good_spike_idx = ismember(spike_templates_0idx,good_templates_idx);
        spike_times = spike_times(good_spike_idx);
        spike_templates_0idx = spike_templates_0idx(good_spike_idx);
        template_amplitudes = template_amplitudes(good_spike_idx);
        spike_depths = spike_depths(good_spike_idx);

        % Rename the spike templates according to the remaining templates
        % (and make 1-indexed from 0-indexed)
        new_spike_idx = nan(max(spike_templates_0idx)+1,1);
        new_spike_idx(good_templates_idx+1) = 1:length(good_templates_idx);
        spike_templates = new_spike_idx(spike_templates_0idx+1);

        %%%%%%%%% Get good single units from bombcell

        % JF code for loading good single units from bombcell
        % unitType: 0 = noise, 1 = good, 2 = multiunit
        curr_qmetrics_path = fullfile(curr_day_path,'ephys','qMetrics');
        load(fullfile(curr_qmetrics_path, 'qMetric.mat'))
        load(fullfile(curr_qmetrics_path, 'param.mat'))
        clearvars unitType;

        % DEFAULT CHANGE: eliminate amplitude cutoff
        % (for one recording it got rid of almost all cells, and
        % cells under amplitude cutoff still look good)
        param.minAmplitude = 0;

        % (classify good cells)
        unitType = nan(length(qMetric.percSpikesMissing), 1);
        unitType( ...
            qMetric.nPeaks > param.maxNPeaks | ...
            qMetric.nTroughs > param.maxNTroughs | ...
            qMetric.somatic ~= param.somatic | ...
            qMetric.spatialDecaySlope <=  param.minSpatialDecaySlope | ...
            qMetric.waveformDuration < param.minWvDuration |...
            qMetric.waveformDuration > param.maxWvDuration  | ...
            qMetric.waveformBaseline >= param.maxWvBaselineFraction) = 0;
        unitType( ...
            any(qMetric.percSpikesMissing <= param.maxPercSpikesMissing, 2)' & ...
            qMetric.nSpikes > param.minNumSpikes & ...
            any(qMetric.Fp <= param.maxRPVviolations, 2)' & ...
            qMetric.rawAmplitude > param.minAmplitude & isnan(unitType)') = 1;
        unitType(isnan(unitType)') = 2;

        % (some upwards waveforms not caught in .somatic? remove)
        upward_waveforms = max(waveforms,[],2) > abs(min(waveforms,[],2));

        % Templates already 1/re-indexed, grab good ones
        good_templates = unitType == 1 & ~upward_waveforms;
        good_templates_idx = find(good_templates);

        % Throw out all non-good template data
        templates = templates(good_templates,:,:);
        template_depths = template_depths(good_templates);
        waveforms = waveforms(good_templates,:);
        templateDuration = templateDuration(good_templates);
        templateDuration_us = templateDuration_us(good_templates);

        % Throw out all non-good spike data
        good_spike_idx = ismember(spike_templates,good_templates_idx);
        spike_times = spike_times(good_spike_idx);
        spike_templates_0idx = spike_templates_0idx(good_spike_idx);
        template_amplitudes = template_amplitudes(good_spike_idx);
        spike_depths = spike_depths(good_spike_idx);

        % Rename the spike templates according to the remaining templates
        % (and make 1-indexed from 0-indexed)
        new_spike_idx = nan(max(spike_templates_0idx)+1,1);
        new_spike_idx(unique(spike_templates_0idx)+1) = 1:length(unique(spike_templates_0idx));
        spike_templates = new_spike_idx(spike_templates_0idx+1);

        %%%%%%%%% Load timeline and associated inputs

        % Get protocols within day
        experiments_dir = dir(curr_day_path);
        experiments_num_idx = cellfun(@(x) ~isempty(x), regexp({experiments_dir.name},'^\d*$'));
        experiment_num = sort(cellfun(@str2num,{experiments_dir(experiments_num_idx).name}));

        exp_protocol = cell(size(experiment_num));
        for curr_exp = 1:length(experiment_num)
            block_dir = dir(fullfile(curr_day_path,num2str(experiment_num(curr_exp)),'*Block*'));
            block_filename = fullfile(block_dir.folder,block_dir.name);
            load(block_filename)
            [~,expDef] = fileparts(block.expDef);
            exp_protocol{curr_exp} = expDef;
        end

        % Load timeline and flipper from task
        load_exp = experiment_num(find(strcmp(exp_protocol,'AP_stimWheelRight'),1,'last'));
        timeline_dir = dir(fullfile(curr_day_path,num2str(experiment_num(load_exp)),'*Timeline.mat'));
        load(fullfile(timeline_dir.folder,timeline_dir.name));

        flipper_name = 'flipper';
        flipper_idx = strcmp({Timeline.hw.inputs.name}, flipper_name);
        flipper_thresh = 2; % TTL threshold
        flipper_trace = Timeline.rawDAQData(:,flipper_idx) > flipper_thresh;
        flipper_flip = find((~flipper_trace(1:end-1) & flipper_trace(2:end)) | ...
            (flipper_trace(1:end-1) & ~flipper_trace(2:end)))+1;
        flipper_flip_times_timeline = Timeline.rawDAQTimestamps(flipper_flip)';

        % Get photodiode flips (compensate for screen flicker)
        photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
        % (define stim screen on from photodiode - sometimes sample-length
        % offset maybe because of backlight onset delay)
        stimScreen_on = Timeline.rawDAQData(:,photodiode_idx) > 0.2;
        stimScreen_on_t = Timeline.rawDAQTimestamps(stimScreen_on);
        photodiode_thresh = 2; % old: max(Timeline.rawDAQData(:,photodiode_idx))/2
        photodiode_trace = Timeline.rawDAQData(stimScreen_on,photodiode_idx) > photodiode_thresh;
        % (medfilt because photodiode can be intermediate value when backlight
        % coming on)
        % (OLD: this worked fine until photodiode bug: sometimes gray)
        %     photodiode_trace_medfilt = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
        %         photodiode_idx),3) > photodiode_thresh;
        %     photodiode_flip = find((~photodiode_trace_medfilt(1:end-1) & photodiode_trace_medfilt(2:end)) | ...
        %         (photodiode_trace_medfilt(1:end-1) & ~photodiode_trace_medfilt(2:end)))+1;
        %     photodiode_flip_times = stimScreen_on_t(photodiode_flip)';
        % (NEW: accomodating photodiode bug flipping sometimes to gray)
        photodiode_trace_medfilt = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
            photodiode_idx),3);
        photodiode_diff_thresh = range(Timeline.rawDAQData(:,photodiode_idx))*0.2;
        photodiode_diff_t = 50; % time (in ms) to get delayed differential
        photodiode_diff_samples = round(Timeline.hw.daqSampleRate/1000*photodiode_diff_t);
        photodiode_diff_filt = [1,zeros(1,photodiode_diff_samples),-1];
        photodiode_trace_diff = abs(conv(photodiode_trace_medfilt,photodiode_diff_filt,'valid')) > ...
            photodiode_diff_thresh;
        photodiode_flip = find(~photodiode_trace_diff(1:end-1) & ...
            photodiode_trace_diff(2:end))+ photodiode_diff_samples + 1;
        photodiode_flip_times = stimScreen_on_t(photodiode_flip)';

        % Get stim on times by closest photodiode flip
        signals_events = block.events;
        n_trials = length(signals_events.endTrialTimes);
        [~,closest_stimOn_photodiode] = ...
            arrayfun(@(x) min(abs(signals_events.stimOnTimes(x) - ...
            photodiode_flip_times)), ...
            1:n_trials);
        stimOn_times = photodiode_flip_times(closest_stimOn_photodiode);

        timeline_reward_idx = strcmp({Timeline.hw.inputs.name}, 'rewardEcho');
        reward_thresh = max(Timeline.rawDAQData(:,timeline_reward_idx))/2;
        reward_trace = Timeline.rawDAQData(:,timeline_reward_idx) > reward_thresh;
        reward_t_timeline = Timeline.rawDAQTimestamps(find(reward_trace(2:end) & ~reward_trace(1:end-1))+1)';

        %%%%%%%%% Convert spike times to timeline

        % These are the digital channels going into the FPGA
        photodiode_sync_idx = 1;
        acqLive_sync_idx = 2;
        led_sync_idx = 3;
        flipper_sync_idx = 4;

        % Load sync/photodiode
        load(([ephys_path filesep 'sync.mat']));

        % Get sync points for alignment
        % Get experiment index by finding numbered folders
        experiment_idx = load_exp == experiment_num;

        % Get flipper experiment differences by long delays
        % (note: this is absolute difference, if recording stopped and
        % started then the clock starts over again, although I thought it
        % wasn't supposed to when I grab the concatenated sync, so
        % something might be wrong)
        flip_diff_thresh = 1; % time between flips to define experiment gap (s)
        flipper_expt_idx = [1;find(abs(diff(sync(flipper_sync_idx).timestamps)) > ...
            flip_diff_thresh)+1;length(sync(flipper_sync_idx).timestamps)+1];

        flipper_flip_times_ephys = sync(flipper_sync_idx).timestamps( ...
            flipper_expt_idx(find(experiment_idx)):flipper_expt_idx(find(experiment_idx)+1)-1);

        % If same number of flips in ephys/timeline, use all
        if length(flipper_flip_times_ephys) ~= length(flipper_flip_times_timeline)
            error('mismatch flipper ephys/timeline');
        end
        sync_timeline = flipper_flip_times_timeline;
        sync_ephys = flipper_flip_times_ephys;

        % Get spike times in timeline time
        spike_times_timeline = interp1(sync_ephys,sync_timeline,spike_times,'linear','extrap');

        %%%%%%%% Save spike times/templates and events
        curr_save_dir = fullfile('D:\Basile_data',strjoin({animal,day},'_'));
        mkdir(curr_save_dir);

        writeNPY(spike_times_timeline,[curr_save_dir filesep 'spike_times.npy']);
        writeNPY(spike_templates,[curr_save_dir filesep 'spike_templates.npy']);
        writeNPY(stimOn_times,[curr_save_dir filesep 'stimOn_times.npy']);
        writeNPY(reward_t_timeline,[curr_save_dir filesep 'reward_times.npy']);

    end

    % Prep next loop
    clearvars -except animals curr_animal

end

%% Naive data

animals = {'AP116','AP117','AP118','AP119'};

for curr_animal = animals

    animal = cell2mat(curr_animal);

    external_drive_path = 'D:\Subjects';

    %% Find data on external hard drive

    % Initialize pathname, add to it with each server location
    days_combined = {};
    days_pathnames_combined = {};

    expInfo_path = fullfile(external_drive_path,animal);
    expInfo_dir = dir(expInfo_path);
    day_paths = cellfun(@(x) ~isempty(regexp(x,'\d\d\d\d-\d\d-\d\d')),{expInfo_dir.name}) &...
        [expInfo_dir.isdir];
    curr_days = {expInfo_dir(day_paths).name};
    curr_days_pathname = cellfun(@(x) [expInfo_path filesep x],curr_days,'uni',false);

    days_combined = [days_combined,curr_days];
    days_pathnames_combined = [days_pathnames_combined,curr_days_pathname];

    [days,unique_day_idx] = unique(days_combined);
    days_pathnames = days_pathnames_combined(unique_day_idx);

    ephys_recorded = false(size(days));
    for curr_day = 1:length(days)
        if exist(fullfile(external_drive_path,animal,days{curr_day},'ephys'),'dir')
            ephys_recorded(curr_day) = true;
        end
    end

    ephys_days = days(ephys_recorded);

    %% Get spike times and templates

    for curr_day_idx = 1:length(ephys_days)

        day = ephys_days{curr_day_idx};

        curr_day_path = fullfile(external_drive_path,animal,day);
        ephys_path = fullfile(curr_day_path,'ephys','kilosort2');

        %%%%%%%%% Load data (copied from AP_load_experiment)

        % Read header information
        header_path = [ephys_path filesep 'dat_params.txt'];
        header_fid = fopen(header_path);
        header_info = textscan(header_fid,'%s %s', 'delimiter',{' = '});
        fclose(header_fid);

        header = struct;
        for i = 1:length(header_info{1})
            header.(header_info{1}{i}) = header_info{2}{i};
        end

        % Load spike data
        if isfield(header,'sample_rate')
            ephys_sample_rate = str2num(header.sample_rate);
        elseif isfield(header,'ap_sample_rate')
            ephys_sample_rate = str2num(header.ap_sample_rate);
        end
        spike_times = double(readNPY([ephys_path filesep 'spike_times.npy']))./ephys_sample_rate;
        spike_templates_0idx = readNPY([ephys_path filesep 'spike_templates.npy']);
        templates_whitened = readNPY([ephys_path filesep 'templates.npy']);
        channel_positions = readNPY([ephys_path filesep 'channel_positions.npy']);
        channel_map = readNPY([ephys_path filesep 'channel_map.npy']);
        winv = readNPY([ephys_path filesep 'whitening_mat_inv.npy']);
        template_amplitudes = readNPY([ephys_path filesep 'amplitudes.npy']);

        % Default channel map/positions are from end: make from surface
        % (hardcode this: kilosort2 drops channels)
        max_depth = 3840;
        channel_positions(:,2) = max_depth - channel_positions(:,2);

        % Unwhiten templates
        templates = zeros(size(templates_whitened));
        for t = 1:size(templates_whitened,1)
            templates(t,:,:) = squeeze(templates_whitened(t,:,:))*winv;
        end

        % Get the waveform of all templates (channel with largest amplitude)
        [~,max_site] = max(max(abs(templates),[],2),[],3);
        templates_max = nan(size(templates,1),size(templates,2));
        for curr_template = 1:size(templates,1)
            templates_max(curr_template,:) = ...
                templates(curr_template,:,max_site(curr_template));
        end
        waveforms = templates_max;

        % Get depth of each template
        % (get min-max range for each channel)
        template_chan_amp = squeeze(abs(diff(prctile(templates,[0,100],2),[],2)));
        % (zero-out low amplitude channels)
        template_chan_amp_thresh = max(template_chan_amp,[],2)*0.5;
        template_chan_amp_overthresh = template_chan_amp.*(template_chan_amp >= template_chan_amp_thresh);
        % (get center-of-mass on thresholded channel amplitudes)
        template_depths = sum(template_chan_amp_overthresh.*channel_positions(:,2)',2)./sum(template_chan_amp_overthresh,2);

        % Get the depth of each spike (templates are zero-indexed)
        spike_depths = template_depths(spike_templates_0idx+1);

        % Get trough-to-peak time for each template
        templates_max_signfix = bsxfun(@times,templates_max, ...
            sign(abs(min(templates_max,[],2)) - abs(max(templates_max,[],2))));

        [~,waveform_trough] = min(templates_max,[],2);
        [~,waveform_peak_rel] = arrayfun(@(x) ...
            max(templates_max(x,waveform_trough(x):end),[],2), ...
            transpose(1:size(templates_max,1)));
        waveform_peak = waveform_peak_rel + waveform_trough;

        templateDuration = waveform_peak - waveform_trough;
        templateDuration_us = (templateDuration/ephys_sample_rate)*1e6;

        % Keep all units
        good_templates_idx = unique(spike_templates_0idx);
        good_templates = ismember(0:size(templates,1)-1,good_templates_idx);

        % Throw out all non-good template data
        templates = templates(good_templates,:,:);
        template_depths = template_depths(good_templates);
        waveforms = waveforms(good_templates,:);
        templateDuration = templateDuration(good_templates);
        templateDuration_us = templateDuration_us(good_templates);

        % Throw out all non-good spike data
        good_spike_idx = ismember(spike_templates_0idx,good_templates_idx);
        spike_times = spike_times(good_spike_idx);
        spike_templates_0idx = spike_templates_0idx(good_spike_idx);
        template_amplitudes = template_amplitudes(good_spike_idx);
        spike_depths = spike_depths(good_spike_idx);

        % Rename the spike templates according to the remaining templates
        % (and make 1-indexed from 0-indexed)
        new_spike_idx = nan(max(spike_templates_0idx)+1,1);
        new_spike_idx(good_templates_idx+1) = 1:length(good_templates_idx);
        spike_templates = new_spike_idx(spike_templates_0idx+1);

        %%%%%%%%% OLD: good units from manual (bombcell not run yet)

         cluster_filepattern = [ephys_path filesep 'cluster_group*'];
         cluster_filedir = dir(cluster_filepattern);
         if ~isempty(cluster_filedir)
             cluster_filename = [ephys_path filesep cluster_filedir.name];
             fid = fopen(cluster_filename);
             cluster_groups = textscan(fid,'%d%s','HeaderLines',1);
             fclose(fid);
         end


        % Check that all used spike templates have a label
        spike_templates_0idx_unique = unique(spike_templates_0idx);
        if ~all(ismember(spike_templates_0idx_unique,uint32(cluster_groups{1}))) || ...
                ~all(ismember(cluster_groups{2},{'good','mua','noise'}))
            warning([animal ' ' day ': not all templates labeled']);
        end
        
        % Define good units from labels
        good_templates_idx = uint32(cluster_groups{1}( ...
            strcmp(cluster_groups{2},'good') | strcmp(cluster_groups{2},'mua')));
        good_templates = ismember(0:size(templates,1)-1,good_templates_idx);        

%         %%%%%%%%% Get good single units from bombcell
% 
%         % JF code for loading good single units from bombcell
%         % unitType: 0 = noise, 1 = good, 2 = multiunit
%         curr_qmetrics_path = fullfile(curr_day_path,'ephys','qMetrics');
%         load(fullfile(curr_qmetrics_path, 'qMetric.mat'))
%         load(fullfile(curr_qmetrics_path, 'param.mat'))
%         clearvars unitType;
% 
%         % DEFAULT CHANGE: eliminate amplitude cutoff
%         % (for one recording it got rid of almost all cells, and
%         % cells under amplitude cutoff still look good)
%         param.minAmplitude = 0;
% 
%         % (classify good cells)
%         unitType = nan(length(qMetric.percSpikesMissing), 1);
%         unitType( ...
%             qMetric.nPeaks > param.maxNPeaks | ...
%             qMetric.nTroughs > param.maxNTroughs | ...
%             qMetric.somatic ~= param.somatic | ...
%             qMetric.spatialDecaySlope <=  param.minSpatialDecaySlope | ...
%             qMetric.waveformDuration < param.minWvDuration |...
%             qMetric.waveformDuration > param.maxWvDuration  | ...
%             qMetric.waveformBaseline >= param.maxWvBaselineFraction) = 0;
%         unitType( ...
%             any(qMetric.percSpikesMissing <= param.maxPercSpikesMissing, 2)' & ...
%             qMetric.nSpikes > param.minNumSpikes & ...
%             any(qMetric.Fp <= param.maxRPVviolations, 2)' & ...
%             qMetric.rawAmplitude > param.minAmplitude & isnan(unitType)') = 1;
%         unitType(isnan(unitType)') = 2;
% 
%         % (some upwards waveforms not caught in .somatic? remove)
%         upward_waveforms = max(waveforms,[],2) > abs(min(waveforms,[],2));
% 
%         % Templates already 1/re-indexed, grab good ones
%         good_templates = unitType == 1 & ~upward_waveforms;
%         good_templates_idx = find(good_templates);

        %%%%%%%%

        % Throw out all non-good template data
        templates = templates(good_templates,:,:);
        template_depths = template_depths(good_templates);
        waveforms = waveforms(good_templates,:);
        templateDuration = templateDuration(good_templates);
        templateDuration_us = templateDuration_us(good_templates);

        % Throw out all non-good spike data
        good_spike_idx = ismember(spike_templates,good_templates_idx);
        spike_times = spike_times(good_spike_idx);
        spike_templates_0idx = spike_templates_0idx(good_spike_idx);
        template_amplitudes = template_amplitudes(good_spike_idx);
        spike_depths = spike_depths(good_spike_idx);

        % Rename the spike templates according to the remaining templates
        % (and make 1-indexed from 0-indexed)
        new_spike_idx = nan(max(spike_templates_0idx)+1,1);
        new_spike_idx(unique(spike_templates_0idx)+1) = 1:length(unique(spike_templates_0idx));
        spike_templates = new_spike_idx(spike_templates_0idx+1);

        %%%%%%%%% Load timeline and associated inputs

        % Get protocols within day
        experiments_dir = dir(curr_day_path);
        experiments_num_idx = cellfun(@(x) ~isempty(x), regexp({experiments_dir.name},'^\d*$'));
        experiment_num = sort(cellfun(@str2num,{experiments_dir(experiments_num_idx).name}));

        block_filename = cell(size(experiment_num));
        exp_protocol = cell(size(experiment_num));
        for curr_exp = 1:length(experiment_num)
            block_dir = dir(fullfile(curr_day_path,num2str(experiment_num(curr_exp)),'*Block*'));
            block_filename{curr_exp} = fullfile(block_dir.folder,block_dir.name);
            load(block_filename{curr_exp})
            [~,expDef] = fileparts(block.expDef);
            exp_protocol{curr_exp} = expDef;
        end

        % Load timeline and flipper from task
        load_exp = experiment_num(find(strcmp(exp_protocol,'AP_lcrGratingPassive'),1,'last'));
        timeline_dir = dir(fullfile(curr_day_path,num2str(experiment_num(load_exp)),'*Timeline.mat'));
        load(fullfile(timeline_dir.folder,timeline_dir.name));

        flipper_name = 'flipper';
        flipper_idx = strcmp({Timeline.hw.inputs.name}, flipper_name);
        flipper_thresh = 2; % TTL threshold
        flipper_trace = Timeline.rawDAQData(:,flipper_idx) > flipper_thresh;
        flipper_flip = find((~flipper_trace(1:end-1) & flipper_trace(2:end)) | ...
            (flipper_trace(1:end-1) & ~flipper_trace(2:end)))+1;
        flipper_flip_times_timeline = Timeline.rawDAQTimestamps(flipper_flip)';

        % Get photodiode flips (compensate for screen flicker)
        photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
        % (define stim screen on from photodiode - sometimes sample-length
        % offset maybe because of backlight onset delay)
        stimScreen_on = Timeline.rawDAQData(:,photodiode_idx) > 0.2;
        stimScreen_on_t = Timeline.rawDAQTimestamps(stimScreen_on);
        photodiode_thresh = 2; % old: max(Timeline.rawDAQData(:,photodiode_idx))/2
        photodiode_trace = Timeline.rawDAQData(stimScreen_on,photodiode_idx) > photodiode_thresh;
        % (medfilt because photodiode can be intermediate value when backlight
        % coming on)
        % (OLD: this worked fine until photodiode bug: sometimes gray)
        %     photodiode_trace_medfilt = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
        %         photodiode_idx),3) > photodiode_thresh;
        %     photodiode_flip = find((~photodiode_trace_medfilt(1:end-1) & photodiode_trace_medfilt(2:end)) | ...
        %         (photodiode_trace_medfilt(1:end-1) & ~photodiode_trace_medfilt(2:end)))+1;
        %     photodiode_flip_times = stimScreen_on_t(photodiode_flip)';
        % (NEW: accomodating photodiode bug flipping sometimes to gray)
        photodiode_trace_medfilt = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
            photodiode_idx),3);
        photodiode_diff_thresh = abs(diff(prctile(Timeline.rawDAQData(:,photodiode_idx),[0,100])))*0.2;
        photodiode_diff_t = 50; % time (in ms) to get delayed differential
        photodiode_diff_samples = round(Timeline.hw.daqSampleRate/1000*photodiode_diff_t);
        photodiode_diff_filt = [1,zeros(1,photodiode_diff_samples),-1];
        photodiode_trace_diff = abs(conv(photodiode_trace_medfilt,photodiode_diff_filt,'valid')) > ...
            photodiode_diff_thresh;
        photodiode_flip = find(~photodiode_trace_diff(1:end-1) & ...
            photodiode_trace_diff(2:end))+ photodiode_diff_samples + 1;
        photodiode_flip_times = stimScreen_on_t(photodiode_flip)';

        % Get stim on times
        load(block_filename{load_exp})
        signals_events = block.events;
        stim_azimuth = [signals_events.stimAzimuthValues];
        % (quick and dirty sanity check: stim off/on + 1 at start)
        if length(photodiode_flip_times) == length(stim_azimuth)*2+1
            stimOn_times = photodiode_flip_times(2:2:end);
        end
      
        %%%%%%%%% Convert spike times to timeline

        % These are the digital channels going into the FPGA
        photodiode_sync_idx = 1;
        acqLive_sync_idx = 2;
        led_sync_idx = 3;
        flipper_sync_idx = 4;

        % Load sync/photodiode
        load(([ephys_path filesep 'sync.mat']));

        % Get sync points for alignment
        % Get experiment index by finding numbered folders
        experiment_idx = load_exp == experiment_num;

        % Get flipper experiment differences by long delays
        % (note: this is absolute difference, if recording stopped and
        % started then the clock starts over again, although I thought it
        % wasn't supposed to when I grab the concatenated sync, so
        % something might be wrong)
        flip_diff_thresh = 1; % time between flips to define experiment gap (s)
        flipper_expt_idx = [1;find(abs(diff(sync(flipper_sync_idx).timestamps)) > ...
            flip_diff_thresh)+1;length(sync(flipper_sync_idx).timestamps)+1];

        flipper_flip_times_ephys = sync(flipper_sync_idx).timestamps( ...
            flipper_expt_idx(find(experiment_idx)):flipper_expt_idx(find(experiment_idx)+1)-1);

        % If same number of flips in ephys/timeline, use all
        if length(flipper_flip_times_ephys) ~= length(flipper_flip_times_timeline)
            error('mismatch flipper ephys/timeline');
        end
        sync_timeline = flipper_flip_times_timeline;
        sync_ephys = flipper_flip_times_ephys;

        % Get spike times in timeline time
        spike_times_timeline = interp1(sync_ephys,sync_timeline,spike_times,'linear','extrap');

        %%%%%%%% Save spike times/templates and events
        curr_save_dir = fullfile('D:\Basile_data',strjoin({animal,day},'_'));
        mkdir(curr_save_dir);

        writeNPY(spike_times_timeline,[curr_save_dir filesep 'spike_times.npy']);
        writeNPY(spike_templates,[curr_save_dir filesep 'spike_templates.npy']);
        writeNPY(stimOn_times,[curr_save_dir filesep 'stimOn_times.npy']);
        writeNPY(stim_azimuth,[curr_save_dir filesep 'stim_azimuth.npy']);

    end

    % Prep next loop
    clearvars -except animals curr_animal

end








