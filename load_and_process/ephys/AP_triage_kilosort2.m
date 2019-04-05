function AP_triage_kilosort2
% Automatic triage of kilosort 2 units

local_phy_dir = 'C:\data_temp\phy';

% Read header information
header_path = [local_phy_dir filesep 'dat_params.txt'];
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
spike_times = double(readNPY([local_phy_dir filesep 'spike_times.npy']))./ephys_sample_rate;
spike_templates_0idx = readNPY([local_phy_dir filesep 'spike_templates.npy']);
templates_whitened = readNPY([local_phy_dir filesep 'templates.npy']);
channel_positions = readNPY([local_phy_dir filesep 'channel_positions.npy']);
channel_map = readNPY([local_phy_dir filesep 'channel_map.npy']);
winv = readNPY([local_phy_dir filesep 'whitening_mat_inv.npy']);
template_amplitudes = readNPY([local_phy_dir filesep 'amplitudes.npy']);

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
[u,s,v] = svd(waveforms','econ');

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

bad_cutoffs = [v_cutoff,large_amp_channel_n_cutoff, ...
    fwhm_trough_cutoff, waveform_deviate_t_cutoff, ...
    trough_peak_t_fwhm_cutoff];
bad_values = [v(:,1),large_amp_channel_n,fwhm_trough, ...
    waveform_deviate_t,trough_peak_t_fwhm];
bad_templates = [v_bad,large_amp_channel_bad,fwhm_trough_bad, ...
    waveform_deviate_t_bad,trough_peak_t_fwhm_bad];
bad_labels = {'SVD','Large amp channel','Trough FWHM', ...
    'Waveform deviation t','Trough-peak/FWHM'};

% Set good templates
triage_good_templates = ~any(bad_templates,2);

% Plot template triage
figure('Name','Kilosort 2 template triage');

for curr_bad = 1:size(bad_templates,2)
    subplot(size(bad_templates,2),2,1+(curr_bad-1)*2); hold on;
    plot(waveforms(bad_templates(:,curr_bad),:)'./ ...
        max(abs(waveforms(bad_templates(:,curr_bad),:)),[],2)','k')
    title(bad_labels{curr_bad});
    
    subplot(size(bad_templates,2),2,2+(curr_bad-1)*2); hold on;
    plot(bad_values(:,curr_bad),'.k');
    line(xlim,repmat(bad_cutoffs(curr_bad),[1,2]));
    xlabel('Template');
    ylabel(bad_labels{curr_bad});
end

% Save classifications in CSV
label_text = cell(size(templates,1),2);
label_text(:,1) = num2cell((1:size(templates,1))-1);
label_text(triage_good_templates,2) = {'good'};
label_text(~triage_good_templates,2) = {'triaged'};

% Write
triage_label = 'AP_triage';
triage_label_filename = [local_phy_dir filesep triage_label '.tsv'];

fid = fopen(triage_label_filename,'w');
fprintf(fid,'%s\t%s\n','cluster_id','Triage');
for curr_template = 1:size(templates,1)
    fprintf(fid,'%d\t%s\n',label_text{curr_template,:});
end
fclose(fid);

disp('Saved triaged templates')





