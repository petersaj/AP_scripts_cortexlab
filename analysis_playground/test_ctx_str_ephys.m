%% Test analysis for wf + cortex ephys + striatum ephys

%% (TESTING)

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





%% Get cortex widefield/ephys relationship by depth (batch)
% (IN PROGRESS: still need to fix depth estimation)
% (hopefully just get rid of regression and use correlation)

animals = {'AP043','AP060','AP061'};
data = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    protocol = 'vanillaChoiceworld';
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
        depth_corr_window = 100; % MUA window in microns
        depth_corr_window_spacing = 50; % MUA window spacing in microns
        
        depth_corr_bins = [ctx_depth(1):depth_corr_window_spacing:(ctx_depth(2)-depth_corr_window); ...
            (ctx_depth(1):depth_corr_window_spacing:(ctx_depth(2)-depth_corr_window))+depth_corr_window];
        depth_corr_bin_centers = depth_corr_bins(1,:) + diff(depth_corr_bins,[],1)/2;
        
        mua_spikes = zeros(size(depth_corr_bins,2),length(spike_binning_t_edges)-1);
        for curr_depth = 1:size(depth_corr_bins,2)
            curr_depth_templates_idx = ...
                find(template_depths >= depth_corr_bins(1,curr_depth) & ...
                template_depths < depth_corr_bins(2,curr_depth));
            
            mua_spikes(curr_depth,:) = histcounts(spike_times_timeline( ...
                ismember(spike_templates,curr_depth_templates_idx)),spike_binning_t_edges);
        end
        
        %%% CROSS-CORRELATION OF FLUORESCENCE AND SPIKES
        maxlag = round(framerate*0.5);
        
        mua_xcorr = nan(size(mua_spikes,1),maxlag*2+1);
        for curr_depth = 1:size(mua_spikes,1)
            [mua_xcorr(curr_depth,:),lags] = xcorr(fluor_roi_interp,mua_spikes(curr_depth,:),maxlag,'unbiased');
        end
        
        sua_xcorr = nan(size(sua_spikes,1),maxlag*2+1);
        for curr_unit = 1:size(sua_spikes,1)
            [sua_xcorr(curr_unit,:),lags] = xcorr(fluor_roi_interp,sua_spikes(curr_unit,:),maxlag,'unbiased');
        end        
        
        % Package data to save
        data(curr_animal,curr_day).mua_depth = depth_corr_bin_centers-ctx_depth(1);
        data(curr_animal,curr_day).sua_depth = template_depths(ctx_units)-ctx_depth(1);
        
        data(curr_animal,curr_day).fluor_mua_xcorr = max(mua_xcorr,[],2);
        data(curr_animal,curr_day).fluor_sua_xcorr = max(sua_xcorr,[],2);
        
        % Clear variables for next experiment
        clearvars -except animals curr_animal animal ...
            protocol flexible_name experiments curr_exp ...
            data
               
    end
    
end

% Plot fluorescence/MUA correlation for each experiment
figure; hold on;
for curr_exp = 1:numel(data)
    plot(data(curr_exp).mua_depth,data(curr_exp).fluor_mua_xcorr);
end

% Concatenate data
mua_depth_cat = horzcat(data(:).mua_depth)';
fluor_mua_xcorr_cat = vertcat(data(:).fluor_mua_xcorr);

% Plot MUA 
figure; 
subplot(2,2,1);
[mua_xcorr_mean,mua_xcorr_sem,mua_group_depth] = grpstats(fluor_mua_xcorr_cat,mua_depth_cat,{'mean','sem','gname'});
mua_group_depth = cellfun(@str2num,mua_group_depth);
errorbar(mua_xcorr_mean,mua_group_depth,mua_xcorr_sem,'horizontal','linewidth',2);
ylabel('Cortical depth');
xlabel('MUA xcorr');
set(gca,'YDir','reverse');

% Plot SUA
sua_depth_groups = 0:100:1400;
sua_depth_group_centers = sua_depth_groups(1:end-1)+diff(sua_depth_groups)/2;
use_data = ~cellfun(@isempty,{data.sua_depth});
sua_depth_cat = {data(use_data).sua_depth};

fluor_sua_xcorr_cat = {data(use_data).fluor_sua_xcorr};
fluor_sua_xcorr_binned = ...
    cell2mat(cellfun(@(sua_depth,fluor_sua_xcorr) ...
    accumarray(discretize(sua_depth,sua_depth_groups), ...
    fluor_sua_xcorr,[length(sua_depth_groups)-1,1],@nanmean,NaN), ...
    sua_depth_cat,fluor_sua_xcorr_cat,'uni',false));

subplot(2,2,3);hold on;
plot(vertcat(data(:).fluor_sua_xcorr),vertcat(data(:).sua_depth),'.k');
errorbar(nanmean(fluor_sua_xcorr_binned,2), ...
    sua_depth_group_centers,AP_sem(fluor_sua_xcorr_binned,2), ...
    'horizontal','linewidth',2,'color','r');
set(gca,'YDir','reverse');
ylabel('Cortical depth');
xlabel('SUA xcorr');











