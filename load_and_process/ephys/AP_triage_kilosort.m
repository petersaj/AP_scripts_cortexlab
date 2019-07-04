function AP_triage_kilosort(animal,day,site)
% AP_triage_kilosort(animal,day,site)
%
% Automatic triage of kilosort units (by template and waveform)
% Creates cluster file in kilosort location (cluster_AP_triage.tsv)
%
% If animal = 'local', works with local phy directory

%% Set paths

if ~exist('site','var')
    site = [];
end

if strcmpi(animal,'local')
    kilosort_path = 'C:\data_temp\phy';
    local_dat = dir([kilosort_path filesep '*.dat']);
    % (if multiple use first, sometimes a CAR file might exist)
    ephys_ap_filename = [local_dat(1).folder filesep local_dat(1).name];
else
    [kilosort_path,kilsort_exists] = AP_cortexlab_filename(animal,day,[],'ephys',site);
    [ephys_ap_filename,ephys_ap_exists] = AP_cortexlab_filename(animal,day,[],'ephys_ap',site);
    
    if ~kilsort_exists || ~ephys_ap_exists
        error('Ephys files not found on server');
    end
end
    

%% Triage by template

% Read header information
header_path = [kilosort_path filesep 'dat_params.txt'];
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
templates_whitened = readNPY([kilosort_path filesep 'templates.npy']);
channel_positions = readNPY([kilosort_path filesep 'channel_positions.npy']);
winv = readNPY([kilosort_path filesep 'whitening_mat_inv.npy']);

% Default channel map/positions are from end: make from surface
channel_positions(:,2) = max(channel_positions(:,2)) - channel_positions(:,2);

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

% (there's a NaN sometimes?)
waveforms(isnan(waveforms)) = 0;

% Triage units the same was as in kilosort 2

% Get SVD of waveforms (first component used)
[u,~,v] = svd(waveforms','econ');

% Standardize sign of first component to be negative
flip_u = -sign(AP_signed_max(u(:,1),1));
u(:,1) = u(:,1)*flip_u;
v(:,1) = v(:,1)*flip_u;

% Get trough/post-trough peak
[waveform_trough,waveform_trough_t] = min(waveforms,[],2);
[waveform_post_peak,waveform_post_peak_t] = arrayfun(@(x) ...
    max(waveforms(x,waveform_trough_t(x):end),[],2), ...
    transpose(1:size(waveforms,1)));
trough_peak_t = (waveform_post_peak_t/ephys_sample_rate)*1e6;

fwhm_trough = (sum(waveforms < (waveform_trough/2),2)/ephys_sample_rate)*1e6;
fwhm_peak = (sum(waveforms > (waveform_post_peak/2),2)/ephys_sample_rate)*1e6;

trough_peak_t_fwhm = trough_peak_t./(fwhm_trough + fwhm_peak);

% Get number of channels with 50% of the max range
template_channel_amp = squeeze(range(templates,2));
amp_thresh = max(template_channel_amp,[],2)*0.5;
large_amp_channel_n = sum(template_channel_amp > amp_thresh,2);

% Get the first timepoint > 10% of the max
waveform_deviate_check = 24; % the first sample can be weird
[~,waveform_deviate_t] = max(abs(waveforms(:,waveform_deviate_check:end)./ ...
    max(abs(waveforms(:,waveform_deviate_check:end)),[],2)) > 0.1,[],2);
waveform_deviate_t = waveform_deviate_t + waveform_deviate_check;

% Set automatic cutoffs and corresponding bad templates

v_cutoff = 0; % upward going spikes might be axons
v_bad = v(:,1) < v_cutoff;

large_amp_channel_n_cutoff = 14; % too many large channels
large_amp_channel_bad = large_amp_channel_n > large_amp_channel_n_cutoff;

fwhm_trough_cutoff = 600; % too wide a trough thickness
fwhm_trough_bad = fwhm_trough > fwhm_trough_cutoff;

waveform_deviate_t_cutoff = 28; % deviates from baseline too early
waveform_deviate_t_bad = waveform_deviate_t < waveform_deviate_t_cutoff;

trough_peak_t_fwhm_cutoff = 3; % non-proportional peak-trough time
trough_peak_t_fwhm_bad = trough_peak_t_fwhm > trough_peak_t_fwhm_cutoff;

bad_template_cutoffs = [v_cutoff,large_amp_channel_n_cutoff, ...
    fwhm_trough_cutoff, waveform_deviate_t_cutoff, ...
    trough_peak_t_fwhm_cutoff];
bad_template_values = [v(:,1),large_amp_channel_n,fwhm_trough, ...
    waveform_deviate_t,trough_peak_t_fwhm];
bad_templates = [v_bad,large_amp_channel_bad,fwhm_trough_bad, ...
    waveform_deviate_t_bad,trough_peak_t_fwhm_bad];
bad_template_labels = {'SVD','Large amp channel','Trough FWHM', ...
    'Waveform deviation t','Trough-peak/FWHM'};

% Plot template triage
figure('Name','Kilosort 2 template triage');

for curr_bad = 1:size(bad_templates,2)
    subplot(size(bad_templates,2),2,1+(curr_bad-1)*2); hold on;
    plot(waveforms(bad_templates(:,curr_bad),:)'./ ...
        max(abs(waveforms(bad_templates(:,curr_bad),:)),[],2)','k')
    title(bad_template_labels{curr_bad});
    
    subplot(size(bad_templates,2),2,2+(curr_bad-1)*2); hold on;
    plot(bad_template_values(:,curr_bad),'.k');
    line(xlim,repmat(bad_template_cutoffs(curr_bad),[1,2]));
    xlabel('Template');
    ylabel(bad_template_labels{curr_bad});
end


%% Triage by waveforms

% Load samples where spikes happen (misnamed "spike times")
spike_samples = readNPY([kilosort_path filesep 'spike_times.npy']);
% Load template belonging to each spike, 1-idx
spike_templates = readNPY([kilosort_path filesep 'spike_templates.npy']) + 1;
% Load channels used by kilosort (it drops channels with no activity), 1-idx
used_channels_idx = readNPY([kilosort_path filesep 'channel_map.npy']) + 1;

% Get AP band file info
ap_dat_dir = dir(ephys_ap_filename);

% Set data format
pull_spike_samples = -40:41;
n_channels = 384;
ephys_datatype = 'int16';

% Load AP-band data (by memory map)
% (get bytes (uint8) per sample (open ephys int16))
dataTypeNBytes = numel(typecast(cast(0, ephys_datatype), 'uint8')); % determine number of bytes per sample
% (get samples per channel)
n_samples = ap_dat_dir.bytes/(n_channels*dataTypeNBytes);  % Number of samples per channel
% (memory map file)
ap_data = memmapfile(ephys_ap_filename,'Format',{ephys_datatype,[n_channels,n_samples],'data'});

% Loop through templates, get average waveforms
max_pull_spikes = 200;

waveform_baseline_std = nan(size(templates,1),size(templates,3));
waveforms_mean = nan(size(templates));
disp('Getting mean waveforms...')
for curr_template = 1:size(templates,1)
    
    % (pull spikes evenly distributed from first to last)
    curr_spikes_idx = find(spike_templates == curr_template);
    
    if isempty(curr_spikes_idx)
        continue
    end
    
    curr_pull_spikes = unique(round(linspace(1,length(curr_spikes_idx),max_pull_spikes)));
    curr_spike_samples = spike_samples(curr_spikes_idx(curr_pull_spikes));
    curr_spike_samples_pull = double(curr_spike_samples) + pull_spike_samples;
    
    out_of_bounds_spikes = any(curr_spike_samples_pull < 1,2) | ...
        any(curr_spike_samples_pull > size(ap_data.data.data,2),2);
    curr_spike_samples_pull(out_of_bounds_spikes,:) = [];
    
    % (grab from file)
    curr_spike_waveforms = reshape(ap_data.data.data(:,reshape(curr_spike_samples_pull',[],1)),n_channels,length(pull_spike_samples),[]);
    
    % (subtract median across channels)
    curr_spike_waveforms_car = curr_spike_waveforms - nanmedian(curr_spike_waveforms,1);
    curr_spike_waveforms_car_sub = curr_spike_waveforms_car - curr_spike_waveforms_car(:,1,:);
    
    % (get "baseline" std)
    baseline_samples = 1:20;
    waveform_baseline_std(curr_template,:) = ...
        std(double(reshape(curr_spike_waveforms_car_sub(used_channels_idx, ...
        baseline_samples,:),length(used_channels_idx),[])),[],2)';   
    
    % (get mean waveform)
    waveforms_mean(curr_template,:,:) = ...
        permute(nanmean(curr_spike_waveforms_car_sub(used_channels_idx,:,:),3),[3,2,1]);
    
    AP_print_progress_fraction(curr_template,size(templates,1)); 
end
disp('done.')

% Get amplitude of peak relative to baseline std
[max_amp,max_channel] = max(max(abs(waveforms_mean),[],2),[],3);
max_channel_std = arrayfun(@(x) waveform_baseline_std(x,max_channel(x)), ...
    1:size(waveforms_mean,1))';
amp_std_ratio = max_amp./max_channel_std;

% Get correlation between mean and template waveform on biggest channels
amp_thresh = 0.5;

template_mean_corr = nan(size(templates,1),1);
for curr_template = 1:size(templates,1)
    % Get mean waveform channels with > threshold amplitude
    mean_channel_amp = squeeze(range(waveforms_mean(curr_template,:,:),2));
    mean_thresh = max(mean_channel_amp)*amp_thresh;
    mean_superthresh_channels = mean_channel_amp > mean_thresh;
    
    corr_grid = corrcoef([permute(templates(curr_template,:,mean_superthresh_channels),[2,3,1]), ...
        permute(waveforms_mean(curr_template,:,mean_superthresh_channels),[2,3,1])]);
    template_mean_corr(curr_template) = nanmean(diag(corr_grid( ...
        1:sum(mean_superthresh_channels),sum(mean_superthresh_channels)+1:end)));
end

% Set automatic cutoffs and corresponding bad waveforms

amp_std_ratio_cutoff = 2; % really small units are questionable
bad_amp_std_ratio = amp_std_ratio < amp_std_ratio_cutoff;

template_mean_corr_cutoff = 0.8; % uncorrelated templates and waveforms suggests something's weird
bad_template_mean_corr = template_mean_corr < template_mean_corr_cutoff;

bad_waveform_cutoffs = [amp_std_ratio_cutoff,template_mean_corr_cutoff];
bad_waveform_values = [amp_std_ratio,template_mean_corr];
bad_waveforms = [bad_amp_std_ratio,bad_template_mean_corr];
bad_waveform_labels = {'Amplitude/baseline std','Template/waveform corr'};

% Plot waveform triage
figure('Name','Kilosort 2 waveform triage');

for curr_bad = 1:size(bad_waveforms,2)
    subplot(size(bad_waveforms,2),2,1+(curr_bad-1)*2); hold on;
    plot(waveforms(bad_waveforms(:,curr_bad),:)'./ ...
        max(abs(waveforms(bad_waveforms(:,curr_bad),:)),[],2)','k')
    title(bad_waveform_labels{curr_bad});
    
    subplot(size(bad_waveforms,2),2,2+(curr_bad-1)*2); hold on;
    plot(bad_waveform_values(:,curr_bad),'.k');
    line(xlim,repmat(bad_waveform_cutoffs(curr_bad),[1,2]));
    xlabel('Template');
    ylabel(bad_waveform_labels{curr_bad});
end


% UNUSED: this is to plot waveforms and templates side-by-side
% % Plot mean waveforms and templates for low correlations
% figure; 
% 
% waveform_axes = subplot(1,2,1); hold on;
% set(waveform_axes,'YDir','reverse'); axis off
% waveform_lines = arrayfun(@(x) plot(waveform_axes,0,0,'k','linewidth',2),1:size(templates,3));
% waveform_xscale = 3;
% waveform_yscale = 0.05;
% 
% template_axes = subplot(1,2,2); hold on;
% set(template_axes,'YDir','reverse'); axis off;
% template_lines = arrayfun(@(x) plot(template_axes,0,0,'k','linewidth',2),1:size(templates,3));
% template_xscale = 3;
% template_yscale = 0.5;
% 
% linkaxes([waveform_axes,template_axes],'y');
% 
% for plot_template = find(template_mean_corr < 0.8)'
%     
%     % Get mean waveform channels with > threshold amplitude
%     mean_channel_amp = squeeze(range(waveforms_mean(plot_template,:,:),2));
%     mean_thresh = max(mean_channel_amp)*amp_thresh;
%     mean_superthresh_channels = mean_channel_amp > mean_thresh;
%     [~,max_channel] = max(max(abs(waveforms_mean(plot_template,:,:)),[],2),[],3);
%     
%     % Plot mean waveforms
%     waveform_y = permute(waveforms_mean(plot_template,:,:),[3,2,1]);
%     waveform_y = -waveform_y*waveform_yscale + channel_positions(:,2);
%     waveform_x = (1:size(waveforms_mean,2)) + channel_positions(:,1)*waveform_xscale;
%     
%     arrayfun(@(ch) set(waveform_lines(ch),'XData',waveform_x(ch,:),'YData',waveform_y(ch,:)),1:length(used_channels_idx));
%     arrayfun(@(ch) set(waveform_lines(ch),'Color','r'),find(mean_superthresh_channels));
%     arrayfun(@(ch) set(waveform_lines(ch),'Color','k'),find(~mean_superthresh_channels));
%     
%     title(waveform_axes,[num2str(plot_template) ', corr: ' num2str(template_mean_corr(plot_template))]);
%     
%     % Plot template waveforms
%     template_y = permute(templates(plot_template,:,:),[3,2,1]);
%     template_y = -template_y*template_yscale + channel_positions(:,2);
%     template_x = (1:size(templates,2)) + channel_positions(:,1)*template_xscale;
%     
%     arrayfun(@(ch) set(template_lines(ch),'XData',template_x(ch,:),'YData',template_y(ch,:)),1:size(templates,3));
%     arrayfun(@(ch) set(template_lines(ch),'Color','r'),find(mean_superthresh_channels));
%     arrayfun(@(ch) set(template_lines(ch),'Color','k'),find(~mean_superthresh_channels));
%         
%     % Set ylim
%     yrange = range(channel_positions(:,2))*0.03.*[-1,1];
%     ylim(template_axes,channel_positions(max_channel,2) + yrange);
%     
%     waitforbuttonpress
%     
% end



%% Save

% Set final good units
triage_good_units = ~any([bad_templates,bad_waveforms],2);

% Save classifications in CSV
label_text = cell(size(templates,1),2);
label_text(:,1) = num2cell((1:size(templates,1))-1);
label_text(triage_good_units,2) = {'good'};
label_text(any(bad_templates,2) & ~any(bad_waveforms,2),2) = {'bad_template'};
label_text(~any(bad_templates,2) & any(bad_waveforms,2),2) = {'bad_waveform'};
label_text(any(bad_templates,2) & any(bad_waveforms,2),2) = {'bad_template_waveform'};

% Write
triage_label = 'cluster_AP_triage';
triage_label_filename = [kilosort_path filesep triage_label '.tsv'];

fid = fopen(triage_label_filename,'w');
fprintf(fid,'%s\t%s\n','cluster_id','AP_triage');
for curr_template = 1:size(templates,1)
    fprintf(fid,'%d\t%s\n',label_text{curr_template,:});
end
fclose(fid);

disp('Saved triaged templates')





