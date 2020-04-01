%% Test analysis for wf + cortex ephys + striatum ephys

%% Testing: get PSTHs

% Load experiment
animal = 'AP061';
day = '2019-12-09';
experiment = 1;
site = 2;
lfp_channel = 'all';
verbose = true;
AP_load_experiment;


% Find cortex end by largest gap between templates
sorted_template_depths = sort([template_depths]);
[max_gap,max_gap_idx] = max(diff(sorted_template_depths));
ctx_end = sorted_template_depths(max_gap_idx)+1;

% Get LFP median-subtracted correlation
lfp_medsub_corr = ...
    corrcoef((movmedian(zscore(double(lfp),[],2),10,1) - ...
    nanmedian(zscore(double(lfp),[],2),1))');

% Get average correlation along the diagonal
on_diag = 0; % (can also use -x:x)
off_diag = 1:20;

lfp_medsub_corr_diag = nanmean(cell2mat(arrayfun(@(x) ...
    padarray(diag(lfp_medsub_corr,off_diag(x)), ...
    abs(off_diag(x)),nan,'pre'), ...
    1:length(off_diag),'uni',false)),2);

% Get average LFP correlation within cortex units
lfp_ctx_units_channels = lfp_channel_positions > min(template_depths) & ...
    lfp_channel_positions < ctx_end;
lfp_medsub_corr_diag_ctx = nanmedian(lfp_medsub_corr_diag(lfp_ctx_units_channels));

% Define the cortex as reaching 70% ctx lfp corr average after cortex
lfp_corr_thresh = lfp_medsub_corr_diag_ctx*0.9;
lfp_corr_rise_idx = ...
    find(lfp_medsub_corr_diag(find(lfp_ctx_units_channels,1):-1:1) > lfp_corr_thresh,1);
ctx_start = lfp_channel_positions(find(lfp_ctx_units_channels,1) - lfp_corr_rise_idx);

% Get final cortex borders
ctx_depth = [ctx_start,ctx_end];

% Plot cortex start/end
figure('Name',[animal ' ' day]);

subplot(1,3,1);
norm_spike_n = mat2gray(log10(accumarray(spike_templates,1)+1));
gscatter(norm_spike_n,template_depths,[],'k',[],10);
xlim([0,1])
set(gca,'YDir','reverse');
xlabel('Norm log_{10} spike rate');
ylabel('Depth (\mum)');
ylim([0,max(channel_positions(:,2))])
line(xlim,repmat(ctx_start,2,1),'color','b','linewidth',2)
line(xlim,repmat(ctx_end,2,1),'color','b','linewidth',2)

subplot(1,3,2:3);
imagesc(lfp_channel_positions,lfp_channel_positions,lfp_medsub_corr)
colormap(brewermap([],'*RdBu'));
caxis([-1,1]);
line(xlim,repmat(ctx_start,2,1),'color','k','linewidth',2)
line(xlim,repmat(ctx_end,2,1),'color','k','linewidth',2)
line(repmat(ctx_start,2,1),ylim,'color','k','linewidth',2)
line(repmat(ctx_end,2,1),ylim,'color','k','linewidth',2)
title('LFP med sub corr');

% Set binning time
skip_seconds = 60;
spike_binning_t = 1/framerate; % seconds
spike_binning_t_edges = frame_t(1)+skip_seconds:spike_binning_t:frame_t(end)-skip_seconds;
spike_binning_t_centers = spike_binning_t_edges(1:end-1) + diff(spike_binning_t_edges)/2;

% Get cortical MUA in sliding window
depth_corr_window = 200; % MUA window in microns
depth_corr_window_spacing = 50; % MUA window spacing in microns

depth_corr_bins = [ctx_depth(1):depth_corr_window_spacing:(ctx_depth(2)-depth_corr_window); ...
    (ctx_depth(1):depth_corr_window_spacing:(ctx_depth(2)-depth_corr_window))+depth_corr_window];
depth_corr_bin_centers = depth_corr_bins(1,:) + diff(depth_corr_bins,[],1)/2;

cortex_mua = zeros(size(depth_corr_bins,2),length(spike_binning_t_centers));
for curr_depth = 1:size(depth_corr_bins,2)
    curr_depth_templates_idx = ...
        find(template_depths >= depth_corr_bins(1,curr_depth) & ...
        template_depths < depth_corr_bins(2,curr_depth));
    
    cortex_mua(curr_depth,:) = histcounts(spike_times_timeline( ...
        ismember(spike_templates,curr_depth_templates_idx)),spike_binning_t_edges);
end

% Get cortical SUA
cortex_units = find(template_depths >= ctx_depth(1) & template_depths <= ctx_depth(2));
cortex_sua = zeros(length(cortex_units),length(spike_binning_t_centers));
for curr_unit = 1:length(cortex_units)
    cortex_sua(curr_unit,:) = histcounts(spike_times_timeline( ...
        ismember(spike_templates,cortex_units(curr_unit))),spike_binning_t_edges);
end

% Load striatal spikes and bin by depth

clear load_parts
load_parts.ephys = true;
site = 1;
AP_load_experiment;

% Set binning time
skip_seconds = 60;
spike_binning_t = 1/framerate; % seconds
spike_binning_t_edges = frame_t(1)+skip_seconds:spike_binning_t:frame_t(end)-skip_seconds;
spike_binning_t_centers = spike_binning_t_edges(1:end-1) + diff(spike_binning_t_edges)/2;

depth_corr_window = 200; % MUA window in microns
depth_corr_window_spacing = 50; % MUA window spacing in microns

depth_corr_bins = [str_depth(1):depth_corr_window_spacing:(str_depth(2)-depth_corr_window); ...
    (str_depth(1):depth_corr_window_spacing:(str_depth(2)-depth_corr_window))+depth_corr_window];
depth_corr_bin_centers = depth_corr_bins(1,:) + diff(depth_corr_bins,[],1)/2;

striatum_mua = zeros(size(depth_corr_bins,2),length(spike_binning_t_edges)-1);
for curr_depth = 1:size(depth_corr_bins,2)
    curr_depth_templates_idx = ...
        find(template_depths >= depth_corr_bins(1,curr_depth) & ...
        template_depths < depth_corr_bins(2,curr_depth));
    
    striatum_mua(curr_depth,:) = histcounts(spike_times_timeline( ...
        ismember(spike_templates,curr_depth_templates_idx)),spike_binning_t_edges);
end

ctx_str_corr = 1-pdist2(cortex_mua,striatum_mua,'correlation');
figure;
subplot(4,4,[1,2,3,5,6,7,9,10,11])
imagesc(ctx_str_corr);
colormap(brewermap([],'*RdBu'));
caxis([-1,1]);
title('MUA correlation');

subplot(4,4,[13,14,15]);
plot(nanmean(ctx_str_corr,1),'k','linewidth',2);
xlabel('Striatal depth');
subplot(4,4,[4,8,12]);
plot(nanmean(ctx_str_corr,2),1:size(ctx_str_corr,1),'k','linewidth',2);
set(gca,'YDir','reverse');
ylabel('Cortical depth');
set(gca,'YAxisLocation','right')

drawnow;

% Get fluorescence in ROI
fVdf_deconv = AP_deconv_wf(fVdf);
fluor_roi = AP_svd_roi(Udf,fVdf_deconv,avg_im);
fluor_roi_interp = interp1(frame_t,fluor_roi,spike_binning_t_centers);

% Plot stim-aligned responses

% Set options
surround_window = [-0.5,3];

surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);

stimIDs = trial_conditions(:,1).*trial_conditions(:,2);
use_stims = find(stimIDs > 0);
use_stimOn_times = stimOn_times(use_stims);
psth_times = bsxfun(@plus, use_stimOn_times(:), surround_time);

fluor_psth = permute(interp1(spike_binning_t_centers,fluor_roi_interp',psth_times),[3,2,1]);
cortex_psth = permute(interp1(spike_binning_t_centers,cortex_mua',psth_times),[3,2,1]);
striatum_psth = permute(interp1(spike_binning_t_centers,striatum_mua',psth_times),[3,2,1]);

baseline_samples = surround_time < 0;

fluor_psth_norm_mean = ...
    nanmean(fluor_psth - nanmean(fluor_psth(:,baseline_samples,:),2),3)./ ...
    nanstd(reshape(fluor_psth(:,baseline_samples,:),size(fluor_psth,1),[]),[],2);
cortex_psth_norm_mean = ...
    nanmean(cortex_psth - nanmean(cortex_psth(:,baseline_samples,:),2),3)./ ...
    nanstd(reshape(cortex_psth(:,baseline_samples,:),size(cortex_psth,1),[]),[],2);
striatum_psth_norm_mean = ...
    nanmean(striatum_psth - nanmean(striatum_psth(:,baseline_samples,:),2),3)./ ...
    nanstd(reshape(striatum_psth(:,baseline_samples,:),size(striatum_psth,1),[]),[],2);

% Regress cortex to multiunit
use_svs = 1:100;
kernel_t = [-0.5,0.5];
kernel_frames = round(kernel_t(1)*framerate):round(kernel_t(2)*framerate);
lambda = 2;
zs = [false,false];
cvfold = 5;
return_constant = false;
use_constant = true;

fVdf_deconv = AP_deconv_wf(fVdf);
fVdf_deconv(isnan(fVdf_deconv)) = 0;
fVdf_deconv_resample = interp1(frame_t,fVdf_deconv(use_svs,:)',spike_binning_t_centers)';
        
cortex_mua_std = cortex_mua./nanstd(cortex_mua,[],2);
cortex_mua_std(isnan(cortex_mua_std)) = 0;
[k,cortex_mua_fluorpred,cortex_mua_explained_var] = ...
    AP_regresskernel(fVdf_deconv_resample, ...
    cortex_mua_std,kernel_frames,lambda,zs,cvfold,return_constant,use_constant);

striatum_mua_std = striatum_mua./nanstd(striatum_mua,[],2);
striatum_mua_std(isnan(striatum_mua_std)) = 0;
[k,striatum_mua_fluorpred,striatum_mua_explained_var] = ...
    AP_regresskernel(fVdf_deconv_resample, ...
    striatum_mua_std,kernel_frames,lambda,zs,cvfold,return_constant,use_constant);

cortex_fluorpred_psth = permute(interp1(spike_binning_t_centers,cortex_mua_fluorpred',psth_times),[3,2,1]);
striatum_fluorpred_psth = permute(interp1(spike_binning_t_centers,striatum_mua_fluorpred',psth_times),[3,2,1]);

cortex_fluorpred_psth_norm_mean = ...
    nanmean(cortex_fluorpred_psth - nanmean(cortex_fluorpred_psth(:,baseline_samples,:),2),3)./ ...
    nanstd(reshape(cortex_fluorpred_psth(:,baseline_samples,:),size(cortex_fluorpred_psth,1),[]),[],2);
striatum_fluorpred_psth_norm_mean = ...
    nanmean(striatum_fluorpred_psth - nanmean(striatum_fluorpred_psth(:,baseline_samples,:),2),3)./ ...
    nanstd(reshape(striatum_fluorpred_psth(:,baseline_samples,:),size(striatum_fluorpred_psth,1),[]),[],2);

figure; 
subplot(2,4,1);
imagesc(cortex_psth_norm_mean);
caxis([-5,5]);
colormap(brewermap([],'*RdBu'));
title('Spikes')
subplot(2,4,2);
imagesc(cortex_fluorpred_psth_norm_mean);
caxis([-5,5]);
colormap(brewermap([],'*RdBu'));
title('Fluor-pred spikes')
subplot(2,4,3);
imagesc(cortex_psth_norm_mean-cortex_fluorpred_psth_norm_mean);
caxis([-5,5]);
colormap(brewermap([],'*RdBu'));
title('Difference')
subplot(2,4,4);
plot(cortex_mua_explained_var.total,1:size(cortex_mua,1),'k','linewidth',2);
set(gca,'YDir','reverse');
xlabel('Fluorescence explained variance')

subplot(2,4,5);
imagesc(striatum_psth_norm_mean);
caxis([-5,5]);
colormap(brewermap([],'*RdBu'));
subplot(2,4,6);
imagesc(striatum_fluorpred_psth_norm_mean);
caxis([-5,5]);
colormap(brewermap([],'*RdBu'));
subplot(2,4,7);
imagesc(striatum_psth_norm_mean-striatum_fluorpred_psth_norm_mean);
caxis([-5,5]);
colormap(brewermap([],'*RdBu'));
subplot(2,4,8);
plot(striatum_mua_explained_var.total,1:size(striatum_mua,1),'k','linewidth',2);
set(gca,'YDir','reverse');
xlabel('Fluorescence explained variance')

% Note: the above is normalized within measured and prediction instead of
% being normalized the same across both, looks more accurate when
% normalizing within

% Regress cortex MUA to striatum MUA
[k,striatum_mua_ctxmuapred,striatum_mua_explained_var] = ...
    AP_regresskernel(cortex_mua_std, ...
    striatum_mua_std,kernel_frames,lambda,zs,cvfold,return_constant,use_constant);

% % Regress cortex SUA to striatum MUA (lots of units makes this take long
% % time - probably not enough memory for long enough recording)
% cortex_sua_std = cortex_sua./nanstd(cortex_sua,[],2);
% cortex_sua_std(isnan(cortex_sua_std)) = 0;
% [k,striatum_mua_ctxsuapred,striatum_sua_explained_var] = ...
%     AP_regresskernel(cortex_sua_std, ...
%     striatum_mua_std,kernel_frames,lambda,zs,cvfold,return_constant,use_constant);

striatum_ctxmuapred_psth = permute(interp1(spike_binning_t_centers,striatum_mua_ctxmuapred',psth_times),[3,2,1]);
str_std_psth = permute(interp1(spike_binning_t_centers,striatum_mua_std',psth_times),[3,2,1]);

figure;
subplot(2,3,1);
imagesc(nanmean(str_std_psth,3));
caxis([-5,5]);
colormap(brewermap([],'*RdBu'));
title('Striatum')

subplot(2,3,2);
imagesc(nanmean(striatum_fluorpred_psth,3));
caxis([-5,5]);
colormap(brewermap([],'*RdBu'));
title('WF-pred striatum')

subplot(2,3,5)
imagesc(nanmean(str_std_psth,3) - nanmean(striatum_fluorpred_psth,3));
caxis([-5,5]);
colormap(brewermap([],'*RdBu'));
title('Diff')

subplot(2,3,3);
imagesc(nanmean(striatum_ctxmuapred_psth,3));
caxis([-5,5]);
colormap(brewermap([],'*RdBu'));
title('AM ephys-pred striatum');

subplot(2,3,6)
imagesc(nanmean(str_std_psth,3) - nanmean(striatum_ctxmuapred_psth,3));
caxis([-5,5]);
colormap(brewermap([],'*RdBu'));
title('Diff')

%% Testing: get CSD to identify layer

% Load experiment (with one LFP channel)
animal = 'AP061';
day = '2019-12-12';
experiment = 4;
site = 2;
lfp_channel = 1;
str_align = 'none';
load_parts.ephys = true;
verbose = true;
AP_load_experiment;


% Load in each channel, light-fix, get average stim response
% (loading script sets up a bunch of this already): load in each channel
use_stims = trial_conditions(:,2) == 90;
use_stimOn_times = stimOn_times(use_stims);

stim_surround_window = [-0.1,0.6];
stim_surround_t = stim_surround_window(1):1/lfp_sample_rate:stim_surround_window(2);

pull_times = use_stimOn_times + stim_surround_t;

lfp_stim_align_avg = nan(n_channels,length(stim_surround_t));

for curr_lfp_channel = 1:n_channels
    
    lfp = double(lfp_memmap{lfp_load_file}.Data.lfp(curr_lfp_channel,lfp_load_start_rel:lfp_load_stop_rel));
        
    % Pull out LFP around light on
    use_blue_on = blue_on >= lfp_t_timeline(1) & blue_on <= lfp_t_timeline(end);
    blue_on_pull_t = blue_on(use_blue_on) + light_surround_t;
    blue_on_lfp = interp1(lfp_t_timeline,lfp',blue_on_pull_t);
    
    use_violet_on = violet_on >= lfp_t_timeline(1) & violet_on <= lfp_t_timeline(end);
    violet_on_pull_t = violet_on(use_violet_on) + light_surround_t;
    violet_on_lfp = interp1(lfp_t_timeline,lfp',violet_on_pull_t);
    
    % Subtract baseline
    baseline_t = find(light_surround_t < 0,1,'last');
    blue_on_lfp_baselinesub = blue_on_lfp - blue_on_lfp(:,baseline_t,:);
    violet_on_lfp_baselinesub = violet_on_lfp - violet_on_lfp(:,baseline_t,:);
    
    % Get rolling median (allow light artifact to change slightly)
    n_light = 500;
    blue_on_lfp_baselinesub_med = movmedian(blue_on_lfp_baselinesub,n_light,1);
    violet_on_lfp_baselinesub_med = movmedian(violet_on_lfp_baselinesub,n_light,1);
    
    % Interpolate out the artifact to remove
    n_lfp_channels = size(lfp,1);
    blue_light_remove = interp1( ...
        reshape(permute(blue_on_pull_t,[2,1]),[],1), ...
        reshape(permute(blue_on_lfp_baselinesub_med,[2,1,3]),[],n_lfp_channels), ...
        reshape(lfp_t_timeline,[],1))';
    violet_light_remove = interp1( ...
        reshape(permute(violet_on_pull_t,[2,1]),[],1), ...
        reshape(permute(violet_on_lfp_baselinesub_med,[2,1,3]),[],n_lfp_channels), ...
        reshape(lfp_t_timeline,[],1))';
    
    % Zero-out any NaNs (e.g. remove nothing)
    blue_light_remove(isnan(blue_light_remove)) = 0;
    violet_light_remove(isnan(violet_light_remove)) = 0;
    
    % Remove the artifact
    lfp_lightfix = lfp - (blue_light_remove + violet_light_remove);
    
    % Get average stim response
    lfp_stim_align = interp1(lfp_t_timeline,lfp_lightfix,pull_times);
    
    baseline_time = stim_surround_t < 0;
    lfp_stim_align_avg(curr_lfp_channel,:) = nanmean((lfp_stim_align - nanmean(lfp_stim_align(:,baseline_time),2)),1);
      
    AP_print_progress_fraction(curr_lfp_channel,n_channels);
    
end

% Average channels at same depth and smooth
n_channel_smooth = 10;
[lfp_stim_depth,lfp_stim_align_avg_depth] = grpstats(lfp_stim_align_avg,lfp_channel_positions,{'gname','mean'});
lfp_stim_depth = cellfun(@str2num,lfp_stim_depth);
lfp_stim_align_avg_depth_smooth = movmean(lfp_stim_align_avg_depth,n_channel_smooth,1);

% Interpolate LFP for line plots
plot_lfp_depths = 0:200:max(lfp_stim_depth);
plot_lfp = interp1(lfp_stim_depth,lfp_stim_align_avg_depth_smooth,plot_lfp_depths);
yscale = -0.2; % flip because YDir is reversed

figure;

% Plot template depths
p1 = subplot(1,3,1);
norm_spike_n = mat2gray(log10(accumarray(spike_templates,1)+1));
scatter(norm_spike_n,template_depths,20,'k','Filled');
ylabel('Depth (\mum)');
xlabel('Norm spike rate');
ylim([min(lfp_channel_positions),max(lfp_channel_positions)]);
set(gca,'YDir','reverse');

% Plot LFP stim response (all channels)
p2 = subplot(1,3,2);
imagesc(stim_surround_t,lfp_stim_depth,lfp_stim_align_avg_depth_smooth);
line([0,0],ylim,'color','k');
colormap(brewermap([],'*RdBu'));
caxis([-max(abs(caxis)),max(abs(caxis))]);

AP_stackplot(plot_lfp',stim_surround_t,plot_lfp_depths,yscale,'k');

% Estimate LFP CSD
lfp_stim_align_avg_depth_smooth_d1 = movmean(diff(lfp_stim_align_avg_depth_smooth,1,1),n_channel_smooth,1);
lfp_stim_align_avg_depth_smooth_d2 = movmean(diff(lfp_stim_align_avg_depth_smooth_d1,1,1),n_channel_smooth,1);
lfp_stim_align_avg_depth_smooth_csd = -lfp_stim_align_avg_depth_smooth_d2; 

lfp_stim_depth_d2 = conv(conv(lfp_stim_depth,[1,1]/2,'valid'),[1,1]/2,'valid');

p3 = subplot(1,3,3);
imagesc(stim_surround_t,lfp_stim_depth_d2,lfp_stim_align_avg_depth_smooth_d2);
line([0,0],ylim,'color','k');
colormap(brewermap([],'*RdBu'));
caxis([-max(abs(caxis)),max(abs(caxis))]);

AP_stackplot(plot_lfp',stim_surround_t,plot_lfp_depths,yscale,'k');

% Link axes
linkaxes([p1,p2,p3],'y')
linkaxes([p2,p3],'x')


%% Get cortex widefield/ephys relationship by depth (batch)

% TO DO HERE: 
% - cortex-fluorescence correlation looks good
% - make same plot for cortex-striatum MUA (needs running, then new plot)

animals = {'AP043','AP060','AP061'};
data = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    protocol = 'vanillaChoiceworld';
%     protocol = 'AP_lcrGratingPassive';   
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    for curr_day = 1:length(experiments)
        
        % Load experiment
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(1); % (only one doubled experiment, and the second one was worse)
        site = 2; % cortex probe is always site 2
        lfp_channel = 'all';
        str_align = 'none'; % (because this is on cortex)
        verbose = false;
        AP_load_experiment
        
        %%% FIND CORTEX BOUNDARIES      
        
        % Find cortex end by largest gap between templates
        sorted_template_depths = sort([template_depths]);
        [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
        ctx_end = sorted_template_depths(max_gap_idx)+1;
        
        % Get LFP median-subtracted correlation
        lfp_medsub_corr = ...
            corrcoef((movmedian(double(lfp - nanmedian(lfp,2)),15,1) - ...
            nanmedian(double(lfp - nanmedian(lfp,2)),1))');
        
        % Get average correlation along the diagonal
        on_diag = 0; % (can also use -x:x)
        off_diag = 1:20;
        
        lfp_medsub_corr_diag = nanmean(cell2mat(arrayfun(@(x) ...
            padarray(diag(lfp_medsub_corr,off_diag(x)), ...
            abs(off_diag(x)),nan,'pre'), ...
            1:length(off_diag),'uni',false)),2);
        
        % Get average LFP correlation within cortex units
        lfp_ctx_units_channels = lfp_channel_positions > min(template_depths) & ...
            lfp_channel_positions < ctx_end;
        lfp_medsub_corr_diag_ctx = nanmedian(lfp_medsub_corr_diag(lfp_ctx_units_channels));
        
%         % Define the cortex start as reaching 70% ctx lfp corr average after cortex
%         lfp_corr_thresh = lfp_medsub_corr_diag_ctx*0.9;
%         lfp_corr_rise_idx = ...
%             find(lfp_medsub_corr_diag(find(lfp_ctx_units_channels,1):-1:1) > lfp_corr_thresh,1);
%         ctx_start = lfp_channel_positions(find(lfp_ctx_units_channels,1) - lfp_corr_rise_idx);        

%         % Define the cortex start as a zero-crossing in correlation
%         % (assume at least first n channels are outside brain)
%         outside_channels = 50:100;
%         lfp_outside_corr = nanmean(lfp_medsub_corr(:,outside_channels),2);
%         lfp_outside_corr_flip = lfp_outside_corr(1:end-1) > 0 & ...
%             lfp_outside_corr(2:end) < 0;
%         ctx_start = lfp_channel_positions(find(lfp_outside_corr_flip > 0,1,'last'));
        
        % Define the cortex start as the first template
        % (the dumbest way and not good if superficial damage)
        ctx_start = sorted_template_depths(1) - 1;
        
        % Get LFP power in brain-ish band
        lfp_beta_power = medfilt1(nanmean(lfp_power(lfp_power_freq > 10 & lfp_power_freq < 60,:),1),10);
        lfp_beta_power_thresh = 0.2*prctile(lfp_beta_power,95);
        lfp_beta_power_flip = lfp_beta_power(1:end-1) < lfp_beta_power_thresh & ...
            lfp_beta_power(2:end) > lfp_beta_power_thresh;
        
        % Get final cortex borders
        ctx_depth = [ctx_start,ctx_end];
        
        % Plot cortex start/end
        figure('Name',[animal ' ' day]);
        
        p1 = subplot(1,4,1);
        norm_spike_n = mat2gray(log10(accumarray(spike_templates,1)+1));
        gscatter(norm_spike_n,template_depths,[],'k',[],10);
        xlim([0,1])
        set(gca,'YDir','reverse');
        xlabel('Norm log_{10} spike rate');
        ylabel('Depth (\mum)');
        line(xlim,repmat(ctx_start,2,1),'color','b','linewidth',2)
        line(xlim,repmat(ctx_end,2,1),'color','b','linewidth',2)
        
        p2 = subplot(1,4,2,'YDir','reverse'); hold on;
        plot(lfp_beta_power,lfp_channel_positions,'k','linewidth',2);
        title('LFP 30-40 Hz');
        xlabel('log_{10} power')
        line(xlim,repmat(ctx_start,2,1),'color','b','linewidth',2)
        line(xlim,repmat(ctx_end,2,1),'color','b','linewidth',2)       
        
        p3 = subplot(1,4,3:4);
        imagesc(lfp_channel_positions,lfp_channel_positions,lfp_medsub_corr)
        colormap(brewermap([],'*RdBu'));
        caxis([-1,1]);
        line(xlim,repmat(ctx_start,2,1),'color','k','linewidth',2)
        line(xlim,repmat(ctx_end,2,1),'color','k','linewidth',2)
        line(repmat(ctx_start,2,1),ylim,'color','k','linewidth',2)
        line(repmat(ctx_end,2,1),ylim,'color','k','linewidth',2)
        title('LFP med sub corr');
        
        ylim([0,max(channel_positions(:,2))])
        linkaxes([p1,p2,p3],'y');
        
        drawnow;
                        
        %%% GET FLUORESCENCE AND SPIKES BY DEPTH        
        
        % Set binning time
        skip_seconds = 60;       
        spike_binning_t = 1/framerate; % seconds
        spike_binning_t_edges = frame_t(1)+skip_seconds:spike_binning_t:frame_t(end)-skip_seconds;
        spike_binning_t_centers = spike_binning_t_edges(1:end-1) + diff(spike_binning_t_edges)/2;
        
        % Get fluorescence in ROI
        fVdf_deconv = AP_deconv_wf(fVdf);
        fluor_roi = AP_svd_roi(Udf,fVdf_deconv,avg_im);
        fluor_roi_interp = interp1(frame_t,fluor_roi,spike_binning_t_centers);
        
        % Get single unit spiking
        ctx_units = find(template_depths >= ctx_depth(1) & template_depths <= ctx_depth(2));
        sua_spikes = zeros(length(ctx_units),length(spike_binning_t_centers));
        for curr_unit_idx = 1:length(ctx_units)
            sua_spikes(curr_unit_idx,:) = histcounts(spike_times_timeline( ...
                spike_templates == ctx_units(curr_unit_idx)),spike_binning_t_edges);
        end
        sua_spikes_std = sua_spikes./std(sua_spikes,[],2);
        sua_spikes_std(isnan(sua_spikes_std)) = 0;
        
        % Set sliding depth window of MUA
        depth_corr_window = 200; % MUA window in microns
        depth_corr_window_spacing = 50; % MUA window spacing in microns
        
        depth_corr_bins = [ctx_depth(1):depth_corr_window_spacing:(ctx_depth(2)-depth_corr_window); ...
            (ctx_depth(1):depth_corr_window_spacing:(ctx_depth(2)-depth_corr_window))+depth_corr_window];
        depth_corr_bin_centers = depth_corr_bins(1,:) + diff(depth_corr_bins,[],1)/2;
        
        cortex_mua = zeros(size(depth_corr_bins,2),length(spike_binning_t_centers));
        for curr_depth = 1:size(depth_corr_bins,2)
            curr_depth_templates_idx = ...
                find(template_depths >= depth_corr_bins(1,curr_depth) & ...
                template_depths < depth_corr_bins(2,curr_depth));
            
            cortex_mua(curr_depth,:) = histcounts(spike_times_timeline( ...
                ismember(spike_templates,curr_depth_templates_idx)),spike_binning_t_edges);
        end
        
        %%% LOAD STRIATUM EPHYS AND GET MUA BY DEPTH        
        clear load_parts
        load_parts.ephys = true;
        site = 1; % (striatum is always on probe 1)
        str_align = 'kernel';
        AP_load_experiment;
        
        striatum_mua = nan(n_aligned_depths,length(spike_binning_t_centers));
        for curr_depth = 1:n_aligned_depths
            curr_spike_times = spike_times_timeline(aligned_str_depth_group == curr_depth);
            % Skip if no spikes at this depth
            if isempty(curr_spike_times)
                continue
            end
            striatum_mua(curr_depth,:) = histcounts(curr_spike_times,spike_binning_t_edges);
        end
        
        %%% CROSS-CORRELATION OF FLUORESCENCE AND SPIKES
        maxlag = round(framerate*0.5);
        
        cortex_fluor_xcorr = nan(size(cortex_mua,1),maxlag*2+1);
        for curr_ctx_depth = 1:size(cortex_mua,1)
            [cortex_fluor_xcorr(curr_ctx_depth,:),lags] = ...
                xcorr(fluor_roi_interp,cortex_mua(curr_ctx_depth,:),maxlag,'unbiased');
        end
        
        cortex_striatum_xcorr = nan(size(cortex_mua,1),maxlag*2+1,n_aligned_depths);
        for curr_ctx_depth = 1:size(cortex_mua,1)
            for curr_str_depth = 1:n_aligned_depths
            [cortex_striatum_xcorr(curr_ctx_depth,:,curr_str_depth),lags] = ...
                xcorr(striatum_mua(curr_str_depth,:),cortex_mua(curr_ctx_depth,:),maxlag,'unbiased');
            end
        end
        
        %%% REGULAR CORRELATION
        cortex_fluor_corr = 1-pdist2(cortex_mua,fluor_roi_interp,'correlation');
        cortex_striatum_corr = 1-pdist2(cortex_mua,striatum_mua,'correlation');      
                
%         %%%%%%%%%%% TESTING
%         
%         figure('Name',[animal ' ' day]);
%         subplot(1,3,1);
%         imagesc(avg_im);
%         colormap(gray);
%         axis image off;
%         drawnow; 
%         
%         % Get and plot fluorescence kernel for cortex spikes        
%         use_svs = 1:100;
%         kernel_t = [-0.5,0.5];
%         kernel_frames = round(kernel_t(1)*framerate):round(kernel_t(2)*framerate);
%         lambda = 10;
%         zs = [false,false];
%         cvfold = 5;
%         return_constant = false;
%         use_constant = true;
%         
%         for curr_plot = 1:2
%             
%             switch curr_plot
%                 case 1
%                     mua_std = cortex_mua./nanstd(cortex_mua,[],2);
%                     mua_std(isnan(mua_std)) = 0;
%                 case 2
%                     mua_std = striatum_mua./nanstd(striatum_mua,[],2);
%                     mua_std(isnan(mua_std)) = 0;
%             end
%             n_depths = size(mua_std,1);
%             
%             fVdf_deconv_resample = interp1(frame_t,fVdf_deconv(use_svs,:)',spike_binning_t_centers)';
%             
%             [k,predicted_spikes,explained_var] = ...
%                 AP_regresskernel(fVdf_deconv_resample, ...
%                 mua_std,kernel_frames,lambda,zs,cvfold,return_constant,use_constant);
%             
%             r_px = zeros(size(Udf,1),size(Udf,2),size(k,2),size(k,3),'single');
%             for curr_spikes = 1:size(k,3)
%                 r_px(:,:,:,curr_spikes) = svdFrameReconstruct(Udf(:,:,use_svs),k(:,:,curr_spikes));
%             end
%             
%             % Get center of mass for each pixel
%             % (get max r for each pixel, filter out big ones)
%             r_px_max = squeeze(max(r_px,[],3));
%             r_px_max(isnan(r_px_max)) = 0;
%             % for i = 1:n_depths
%             %     r_px_max(:,:,i) = medfilt2(r_px_max(:,:,i),[10,10]);
%             % end
%             r_px_max_norm = bsxfun(@rdivide,r_px_max, ...
%                 permute(max(reshape(r_px_max,[],n_depths),[],1),[1,3,2]));
%             r_px_max_norm(isnan(r_px_max_norm)) = 0;
%             r_px_com = sum(bsxfun(@times,r_px_max_norm,permute(1:n_depths,[1,3,2])),3)./sum(r_px_max_norm,3);
%             
%             % Plot map of cortical pixel by preferred depth of probe
%             r_px_com_col = ind2rgb(round(mat2gray(r_px_com,[1,n_depths])*255),jet(255));
%             
%             a1 = subplot(1,3,curr_plot+1,'YDir','reverse'); hold on;
%             imagesc(avg_im); colormap(gray); caxis([0,prctile(avg_im(:),99.7)]);
%             axis off; axis image;
%             p = imagesc(r_px_com_col);
%             axis off; axis image;
%             set(p,'AlphaData',mat2gray(max(r_px_max_norm,[],3), ...
%                 [0,double(prctile(reshape(max(r_px_max_norm,[],3),[],1),95))]));
%             
%             drawnow;
%             
%         end
%         
%         %%%%%%%%%%%%%
        
        
        % Package data to save
        data(curr_animal,curr_day).cortex_mua_depth = depth_corr_bin_centers-ctx_depth(1);        
%         data(curr_animal,curr_day).cortex_fluor_xcorr = max(cortex_fluor_xcorr,[],2);
%         data(curr_animal,curr_day).cortex_str_xcorr = permute(max(cortex_striatum_xcorr,[],2),[1,3,2]);
        data(curr_animal,curr_day).cortex_fluor_corr = cortex_fluor_corr;
        data(curr_animal,curr_day).cortex_striatum_corr = cortex_striatum_corr;
        
        % Clear variables for next experiment
        clearvars -except animals curr_animal animal ...
            protocol flexible_name experiments curr_exp ...
            data
               
    end
    
end

use_data = cellfun(@(x) ~isempty(x),{data.cortex_mua_depth});

% Concatenate data
% % (raw)
% mua_depth_cat = round(horzcat(data(:).mua_depth))';
% fluor_mua_xcorr_cat = vertcat(data(:).fluor_mua_xcorr);

% (0-1 normalized)
mua_depth_norm = cellfun(@(x) round(mat2gray(x)'*10)/10, ...
    {data(use_data).cortex_mua_depth},'uni',false)';
cortex_fluor_corr_norm = cellfun(@mat2gray, ...
    {data(use_data).cortex_fluor_corr},'uni',false)';

% Plot MUA 
figure; 

subplot(1,2,1); hold on;
for curr_exp = 1:sum(use_data)
  plot(cortex_fluor_corr_norm{curr_exp},mua_depth_norm{curr_exp},'color',[0.5,0.5,0.5]);
end
[mua_corr_mean,mua_corr_sem,mua_group_depth] = grpstats( ...
    vertcat(cortex_fluor_corr_norm{:}),vertcat(mua_depth_norm{:}),{'mean','sem','gname'});
mua_group_depth = cellfun(@str2num,mua_group_depth);
errorbar(mua_corr_mean,mua_group_depth,mua_corr_sem,'horizontal','linewidth',2,'color','k');
ylabel('Cortical depth (max normalized)');
xlabel('MUA-Fluorescence correlation (max normalized)');
set(gca,'YDir','reverse');

% Plot cortex/striatum MUA correlation
cortex_striatum_corr_norm = cellfun(@(x) x, ...
    {data(use_data).cortex_striatum_corr},'uni',false)';

subplot(1,2,2); hold on;
set(gca,'ColorOrder',copper(4));
% for curr_exp = 1:sum(use_data)
%   plot(cortex_striatum_corr_norm{curr_exp},mua_depth_norm{curr_exp});
% end

cortex_striatum_corr_norm_cat = vertcat(cortex_striatum_corr_norm{:});
for curr_depth = 1:size(cortex_striatum_corr_norm_cat,2)
    [mua_corr_mean,mua_corr_sem,mua_group_depth] = grpstats( ...
        cortex_striatum_corr_norm_cat(:,curr_depth),vertcat(mua_depth_norm{:}),{'mean','sem','gname'});
    mua_group_depth = cellfun(@str2num,mua_group_depth);
    errorbar(mua_corr_mean,mua_group_depth,mua_corr_sem,'horizontal','linewidth',2);
end
set(gca,'YDir','reverse');










