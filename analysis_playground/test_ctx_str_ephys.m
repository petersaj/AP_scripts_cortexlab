%% Test analysis for recording cortex + striatum (+wf) simultaneous ephys


%% (TESTING) 

% Load experiment
animal = 'AP061';
day = '2019-12-10';
experiment = 1;
site = 2;
lfp_channel = 'all';
verbose = true;
AP_load_experiment;

% Get LFP median-subtracted correlation
lfp_medsub_corr = ...
    corrcoef((movmedian(zscore(double(lfp),[],2),10,1) - ...
    nanmedian(zscore(double(lfp),[],2),1))');
lfp_medsub_corr(triu(true(length(lfp_medsub_corr)),0)) = nan;

% Find cortex start by dim2 correlation falling below 70% of range
% (have a buffer from the start: channels across air/saline border have a
% big shift in correlation too)
saline_start_buffer = 500; % microns
saline_start_buffer_channel = find(lfp_channel_positions > saline_start_buffer,1);

lfp_medsub_corr_avg_2 = smooth(nanmean(lfp_medsub_corr,2),10);

corr_drop = find(lfp_medsub_corr_avg_2(saline_start_buffer_channel:end) > ...
    min(lfp_medsub_corr_avg_2(saline_start_buffer_channel:end)) + ...
    range(lfp_medsub_corr_avg_2(saline_start_buffer_channel:end))*0.7,1,'last');
ctx_start = lfp_channel_positions(saline_start_buffer_channel + corr_drop);

% Find cortex end by largest gap between templates
sorted_template_depths = sort([template_depths]);
[max_gap,max_gap_idx] = max(diff(sorted_template_depths));
ctx_end = sorted_template_depths(max_gap_idx)+1;

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


%%%%%%%%% TESTING: REGRESSION/PREDICTED RESPONSES

% Set upsample value for regression
upsample_factor = 1;
sample_rate = (1/median(diff(frame_t)))*upsample_factor;

% Skip the first/last n seconds to do this
skip_seconds = 60;
time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% (to group multiunit by depth within striatum)
n_depths = round(diff(ctx_depth)/200);
% n_depths = 1;
depth_group_edges = round(linspace(ctx_depth(1),ctx_depth(2),n_depths+1));
[depth_group_n,depth_group] = histc(spike_depths,depth_group_edges);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

binned_spikes = zeros(n_depths,length(time_bins)-1);
for curr_depth = 1:n_depths  
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);    
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);  
end

binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
binned_spikes_std(isnan(binned_spikes_std)) = 0;

use_svs = 1:200;
kernel_t = [-0.5,0.5];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
lambda = 2;
zs = [false,false];
cvfold = 5;
return_constant = false;
use_constant = true;

fVdf_resample = interp1(frame_t,fVdf(use_svs,:)',time_bin_centers)';
dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf(use_svs,:),[],2)',time_bin_centers)';
        
fVdf_deconv = AP_deconv_wf(fVdf);
fVdf_deconv(isnan(fVdf_deconv)) = 0;
fVdf_deconv_resample = interp1(frame_t,fVdf_deconv(use_svs,:)',time_bin_centers)';
        
% % TO USE fV
% [k,predicted_spikes,explained_var] = ...
%     AP_regresskernel(fVdf_resample, ...
%     binned_spikes_std,kernel_frames,lambda,zs,cvfold);
% % TO USE dfV
% [k,predicted_spikes,explained_var] = ...
%     AP_regresskernel(dfVdf_resample, ...
%     binned_spikes_std,kernel_frames,lambda,zs,cvfold,return_constant,use_constant);
% TO USE DECONV
[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(fVdf_deconv_resample, ...
    binned_spikes_std,kernel_frames,lambda,zs,cvfold,return_constant,use_constant);

% Convert kernel to pixel space
r_px = zeros(size(Udf,1),size(Udf,2),size(k,2),size(k,3),'single');
for curr_spikes = 1:size(k,3)
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(Udf(:,:,use_svs),k(:,:,curr_spikes));
end

AP_image_scroll(r_px,kernel_frames/framerate);
caxis([-prctile(r_px(:),99.9),prctile(r_px(:),99.9)])
colormap(colormap_BlueWhiteRed);
axis image;

% Get center of mass for each pixel
% (get max r for each pixel, filter out big ones)
r_px_max = squeeze(max(r_px,[],3));
r_px_max(isnan(r_px_max)) = 0;
% for i = 1:n_depths
%     r_px_max(:,:,i) = medfilt2(r_px_max(:,:,i),[10,10]);
% end
r_px_max_norm = bsxfun(@rdivide,r_px_max, ...
    permute(max(reshape(r_px_max,[],n_depths),[],1),[1,3,2]));
r_px_max_norm(isnan(r_px_max_norm)) = 0;
r_px_com = sum(bsxfun(@times,r_px_max_norm,permute(1:n_depths,[1,3,2])),3)./sum(r_px_max_norm,3);

% Plot map of cortical pixel by preferred depth of probe
r_px_com_col = ind2rgb(round(mat2gray(r_px_com,[1,n_depths])*255),jet(255));
figure;
a1 = axes('YDir','reverse');
imagesc(avg_im); colormap(gray); caxis([0,prctile(avg_im(:),99.7)]);
axis off; axis image;
a2 = axes('Visible','off'); 
p = imagesc(r_px_com_col);
axis off; axis image;
set(p,'AlphaData',mat2gray(max(r_px_max_norm,[],3), ...
     [0,double(prctile(reshape(max(r_px_max_norm,[],3),[],1),95))]));
set(gcf,'color','w');

c1 = colorbar('peer',a1,'Visible','off');
c2 = colorbar('peer',a2);
ylabel(c2,'Depth (\mum)');
colormap(c2,jet);
set(c2,'YDir','reverse');
set(c2,'YTick',linspace(0,1,n_depths));
set(c2,'YTickLabel',round(linspace(depth_group_edges(1),depth_group_edges(end),n_depths)));
 

% Plot measured/widefield-predicted spikes

% Set options
surround_window = [-0.5,3];
baseline_window = [-0.1,0];

surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);
baseline_surround_time = baseline_window(1):surround_samplerate:baseline_window(2);

% (passive)
use_stims = find(stimIDs == 3);
% (choiceworld)
% stimIDs = trial_conditions(:,1).*trial_conditions(:,2);
% use_stims = find(stimIDs > 0);

use_stimOn_times = stimOn_times(use_stims);

stim_surround_times = bsxfun(@plus, use_stimOn_times(:), surround_time);
stim_baseline_surround_times = bsxfun(@plus, use_stimOn_times(:), baseline_surround_time);

psth_measured_stim = permute(interp1(time_bin_centers,binned_spikes_std',stim_surround_times),[3,2,1]);
psth_predicted_stim = permute(interp1(time_bin_centers,predicted_spikes',stim_surround_times),[3,2,1]);

psth_measured_baseline = permute(nanmean(interp1(time_bin_centers,binned_spikes_std',stim_baseline_surround_times),2),[3,2,1]);
psth_predicted_baseline = permute(nanmean(interp1(time_bin_centers,predicted_spikes',stim_baseline_surround_times),2),[3,2,1]);

psth_measured = nanmean(psth_measured_stim - psth_measured_baseline,3);
psth_predicted = nanmean(psth_predicted_stim - psth_predicted_baseline,3);

figure; hold on;
AP_stackplot(psth_measured',surround_time,3,false,'k');
AP_stackplot(psth_predicted',surround_time,3,false,[0,0.7,0]);
% plot(surround_time,psth_measured,'k','linewidth',2);
% plot(surround_time,psth_predicted,'r','linewidth',2);



%%%%%%%%% TESTING: FLUORESCENCE BY DEPTH

%%% Get sliding MUA
depth_corr_window = 100; % MUA window in microns
depth_corr_window_spacing = 50; % MUA window spacing in microns

depth_corr_bins = [ctx_depth(1):depth_corr_window_spacing:(ctx_depth(2)-depth_corr_window); ...
    (ctx_depth(1):depth_corr_window_spacing:(ctx_depth(2)-depth_corr_window))+depth_corr_window];
depth_corr_bin_centers = depth_corr_bins(1,:) + diff(depth_corr_bins,[],1)/2;

skip_seconds = 60;

spike_binning_t = 1/framerate; % seconds
spike_binning_t_edges = frame_t(1)+skip_seconds:spike_binning_t:frame_t(end)-skip_seconds;
spike_binning_t_centers = spike_binning_t_edges(1:end-1) + diff(spike_binning_t_edges)/2;

binned_spikes_depth = zeros(size(depth_corr_bins,2),length(spike_binning_t_edges)-1);
for curr_depth = 1:size(depth_corr_bins,2)
    curr_depth_templates_idx = ...
        find(template_depths >= depth_corr_bins(1,curr_depth) & ...
        template_depths < depth_corr_bins(2,curr_depth));
    
    binned_spikes_depth(curr_depth,:) = histcounts(spike_times_timeline( ...
        ismember(spike_templates,curr_depth_templates_idx)),spike_binning_t_edges);
end

fluor_roi = AP_svd_roi(Udf,fVdf_deconv,avg_im);
fluor_roi_interp = interp1(frame_t,fluor_roi,spike_binning_t_centers);

fluor_mua_corr = corrcoef([fluor_roi_interp',binned_spikes_depth']);

figure;
subplot(1,2,1);
imagesc(fluor_mua_corr);
caxis([-1,1])
colormap(colormap_BlueWhiteRed);
axis image;

subplot(1,2,2);
plot(depth_corr_bin_centers,fluor_mua_corr(1,2:end),'k','linewidth',2);
ylabel('Fluor-MUA correlation');
xlabel('Cortex depth (\mum)');

% Regress from multiunit to fluorescence
binned_spikes_depth_std = binned_spikes_depth./nanstd(binned_spikes_depth,[],2);
binned_spikes_depth_std(isnan(binned_spikes_depth_std)) = 0;
fluor_roi_interp_std = fluor_roi_interp./nanstd(fluor_roi_interp);

n_lambdas = 10;
lambda_range = [0,300];
lambdas = linspace(lambda_range(1),lambda_range(2),n_lambdas)';
explained_var_lambdas = nan(n_lambdas,1);

figure;
h = plot(lambdas,explained_var_lambdas,'k','linewidth',2);drawnow;
for curr_lambda = 1:n_lambdas
    [~,~,curr_explained_var] = ...
        AP_regresskernel(binned_spikes_depth_std, ...
        fluor_roi_interp_std,kernel_frames,lambdas(curr_lambda),zs,cvfold,return_constant,use_constant);
    explained_var_lambdas(curr_lambda) = curr_explained_var.total;
    set(h,'YData',explained_var_lambdas);
    drawnow;
end

all_lambdas = lambda_range(1):lambda_range(end);
explained_var_lambdas_interp = interp1(lambdas,explained_var_lambdas,all_lambdas,'pchip');
plot(all_lambdas,explained_var_lambdas_interp,'r');

[~,lambda_idx] = max(explained_var_lambdas_interp);
lambda = all_lambdas(lambda_idx);

fluor_roi_interp_std = fluor_roi_interp./nanstd(fluor_roi_interp);
[mua_k,mua_predicted_fluor,explained_var] = ...
    AP_regresskernel(binned_spikes_depth_std, ...
    fluor_roi_interp_std,kernel_frames,lambda,zs,cvfold,return_constant,use_constant);

% Regress from single units to fluorescence
ctx_units = find(template_depths >= ctx_depth(1) & template_depths <= ctx_depth(2));
unit_spikes = zeros(length(ctx_units),length(spike_binning_t_centers));
for curr_unit_idx = 1:length(ctx_units)
    unit_spikes(curr_unit_idx,:) = histcounts(spike_times_timeline( ...
        spike_templates == ctx_units(curr_unit_idx)),spike_binning_t_edges);
end
unit_spikes_std = unit_spikes./std(unit_spikes,[],2);
unit_spikes_std(isnan(unit_spikes_std)) = 0;

n_lambdas = 10;
lambda_range = [0,300];
lambdas = linspace(lambda_range(1),lambda_range(2),n_lambdas)';
explained_var_lambdas = nan(n_lambdas,1);

figure; hold on;
h = plot(lambdas,explained_var_lambdas,'k','linewidth',2);drawnow;
for curr_lambda = 1:n_lambdas
    [~,~,curr_explained_var] = ...
        AP_regresskernel(unit_spikes_std, ...
        fluor_roi_interp_std,kernel_frames,lambdas(curr_lambda),zs,cvfold,return_constant,use_constant);
    explained_var_lambdas(curr_lambda) = curr_explained_var.total;
    set(h,'YData',explained_var_lambdas);
    drawnow;
end

all_lambdas = lambda_range(1):lambda_range(end);
explained_var_lambdas_interp = interp1(lambdas,explained_var_lambdas,all_lambdas,'pchip');
plot(all_lambdas,explained_var_lambdas_interp,'r');

[~,lambda_idx] = max(explained_var_lambdas_interp);
lambda = all_lambdas(lambda_idx);

[sua_k,sua_predicted_fluor,explained_var] = ...
    AP_regresskernel(unit_spikes_std, ...
    fluor_roi_interp_std,kernel_frames,lambda,zs,cvfold,return_constant,use_constant);

% Plot MUA/SUA kernels
figure;
subplot(2,2,1);
imagesc([],depth_corr_bin_centers-ctx_depth(1),mua_k);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
ylabel('MUA depth from ctx surface');
title('MUA > fluorescence kernel');

subplot(2,2,2);
[~,sort_idx] = sort(template_depths(ctx_units));
imagesc(sua_k(sort_idx,:));
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
ylabel('MUA depth (sorted)');
title('MUA > fluorescence kernel');

subplot(2,2,3);
plot(max(mua_k,[],2),depth_corr_bin_centers-ctx_depth(1),'k','linewidth',2);
ylabel('MUA depth from ctx surface');
set(gca,'YDir','reverse');

subplot(2,2,4);
plot(max(sua_k,[],2),template_depths(ctx_units)-ctx_depth(1),'.k');
ylabel('MUA depth from ctx surface');
set(gca,'YDir','reverse');





%% Get cortex widefield/ephys relationship by depth (batch)

animals = {'AP043','AP060','AP061'};
data = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    protocol = 'vanillaChoiceworld';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    for curr_exp = 1:length(experiments)
        
        % Load experiment
        day = experiments(curr_exp).day;
        % (only one doubled experiment, and the second one was worse)
        experiment = experiments(curr_exp).experiment(1);
        site = 2; % cortex probe is always site 2
        lfp_channel = 'all';
        verbose = false;
        AP_load_experiment
          
        %%% FIND CORTEX BOUNDARIES
        
        % Get LFP median-subtracted correlation
        lfp_medsub_corr = ...
            corrcoef((movmedian(zscore(double(lfp),[],2),10,1) - ...
            nanmedian(zscore(double(lfp),[],2),1))');
        lfp_medsub_corr(triu(true(length(lfp_medsub_corr)),0)) = nan;
        
        % Find cortex start by dim2 correlation falling below 70% of range
        % (have a buffer from the start: channels across air/saline border have a
        % big shift in correlation too)
        saline_start_buffer = 500; % microns
        saline_start_buffer_channel = find(lfp_channel_positions > saline_start_buffer,1);
        
        lfp_medsub_corr_avg_2 = smooth(nanmean(lfp_medsub_corr,2),10);
        
        corr_drop = find(lfp_medsub_corr_avg_2(saline_start_buffer_channel:end) > ...
            min(lfp_medsub_corr_avg_2(saline_start_buffer_channel:end)) + ...
            range(lfp_medsub_corr_avg_2(saline_start_buffer_channel:end))*0.7,1,'last');
        ctx_start = lfp_channel_positions(saline_start_buffer_channel + corr_drop);       
        
        % Find cortex end by largest gap between templates
        sorted_template_depths = sort([template_depths]);
        [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
        ctx_end = sorted_template_depths(max_gap_idx)+1;
        
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
        
        ctx_depth = [ctx_start,ctx_end];
        
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
        maxlag = framerate*0.5;
        
        mua_xcorr = nan(size(mua_spikes,1),maxlag*2+1);
        for curr_depth = 1:size(mua_spikes,1)
            [mua_xcorr(curr_depth,:),lags] = xcorr(fluor_roi_interp,mua_spikes(curr_depth,:),maxlag,'unbiased');
        end
        
        sua_xcorr = nan(size(sua_spikes,1),maxlag*2+1);
        for curr_unit = 1:size(sua_spikes,1)
            [sua_xcorr(curr_unit,:),lags] = xcorr(fluor_roi_interp,sua_spikes(curr_unit,:),maxlag,'unbiased');
        end
    
        % Regression setup
        sua_spikes_std = sua_spikes./nanstd(sua_spikes,[],2);
        sua_spikes_std(isnan(sua_spikes_std)) = 0;
        
        mua_spikes_std = mua_spikes./nanstd(mua_spikes,[],2);
        mua_spikes_std(isnan(mua_spikes_std)) = 0;
        
        fluor_roi_interp_std = fluor_roi_interp./nanstd(fluor_roi_interp);
        
        kernel_frames = round(-0.5*framerate):round(0.5*framerate);
        zs = [false,false];
        cvfold = 5;
        return_constant = false;
        use_constant = true;
        
        n_lambdas = 10;
        lambda_range = [0,300];
        lambdas = linspace(lambda_range(1),lambda_range(2),n_lambdas)';
        explained_var_lambdas = nan(n_lambdas,1);
        
        % Regress multiunit to fluorescence
        figure;
        h = plot(lambdas,explained_var_lambdas,'k','linewidth',2);drawnow;
        for curr_lambda = 1:n_lambdas
            [~,~,curr_explained_var] = ...
                AP_regresskernel(mua_spikes_std, ...
                fluor_roi_interp_std,kernel_frames,lambdas(curr_lambda),zs,cvfold,return_constant,use_constant);
            explained_var_lambdas(curr_lambda) = curr_explained_var.total;
            set(h,'YData',explained_var_lambdas);
            drawnow;
        end
        
        all_lambdas = lambda_range(1):lambda_range(end);
        explained_var_lambdas_interp = interp1(lambdas,explained_var_lambdas,all_lambdas,'pchip');
        plot(all_lambdas,explained_var_lambdas_interp,'r');
        
        [~,lambda_idx] = max(explained_var_lambdas_interp);
        lambda = all_lambdas(lambda_idx);
        
        [mua_k,mua_predicted_fluor,~] = ...
            AP_regresskernel(binned_spikes_depth_std, ...
            fluor_roi_interp_std,kernel_frames,lambda,zs,cvfold,return_constant,use_constant);
        
        % Regress single units to fluorescence
        figure; hold on;
        h = plot(lambdas,explained_var_lambdas,'k','linewidth',2);drawnow;
        for curr_lambda = 1:n_lambdas
            [~,~,curr_explained_var] = ...
                AP_regresskernel(sua_spikes_std, ...
                fluor_roi_interp_std,kernel_frames,lambdas(curr_lambda),zs,cvfold,return_constant,use_constant);
            explained_var_lambdas(curr_lambda) = curr_explained_var.total;
            set(h,'YData',explained_var_lambdas);
            drawnow;
        end
        
        all_lambdas = lambda_range(1):lambda_range(end);
        explained_var_lambdas_interp = interp1(lambdas,explained_var_lambdas,all_lambdas,'pchip');
        plot(all_lambdas,explained_var_lambdas_interp,'r');
        
        [~,lambda_idx] = max(explained_var_lambdas_interp);
        lambda = all_lambdas(lambda_idx);
        
        [sua_k,sua_predicted_fluor,~] = ...
            AP_regresskernel(sua_spikes_std, ...
            fluor_roi_interp_std,kernel_frames,lambda,zs,cvfold,return_constant,use_constant);
        
        % Gather and save metrics
        data(curr_animal,curr_day).mua_depth = depth_corr_bin_centers-ctx_depth(1);
        data(curr_animal,curr_day).sua_depth = template_depths(ctx_units)-ctx_depth(1);
        
        data(curr_animal,curr_day).fluor_mua_xcorr = max(mua_xcorr,[],2);
        data(curr_animal,curr_day).fluor_sua_xcorr = max(sua_xcorr,[],2);
        
        data(curr_animal,curr_day).fluor_mua_weight = max(mua_k,[],2);
        data(curr_animal,curr_day).fluor_sua_weight = max(sua_k,[],2);
        
        % Plot MUA/SUA kernels
        figure;
        subplot(2,2,1);
        imagesc([],depth_corr_bin_centers-ctx_depth(1),mua_k);
        caxis([-max(abs(caxis)),max(abs(caxis))]);
        colormap(brewermap([],'*RdBu'));
        ylabel('MUA depth from ctx surface');
        title('MUA > fluorescence kernel');
        
        subplot(2,2,2);
        [~,sort_idx] = sort(template_depths(ctx_units));
        imagesc(sua_k(sort_idx,:));
        caxis([-max(abs(caxis)),max(abs(caxis))]);
        colormap(brewermap([],'*RdBu'));
        ylabel('MUA depth (sorted)');
        title('MUA > fluorescence kernel');
        
        subplot(2,2,3); hold on;
        plot(max(sua_k,[],2),template_depths(ctx_units)-ctx_depth(1),'.k');
        plot(max(mua_k,[],2),depth_corr_bin_centers-ctx_depth(1),'r','linewidth',2);
        ylabel('Depth from ctx surface');
        xlabel('Max weight');
        set(gca,'YDir','reverse');
        legend({'SUA','MUA'});

        subplot(2,2,4); hold on;
        plot(max(sua_xcorr,[],2),template_depths(ctx_units)-ctx_depth(1),'.k');
        plot(max(mua_xcorr,[],2),depth_corr_bin_centers-ctx_depth(1),'r','linewidth',2);
        ylabel('Depth from ctx surface');
        xlabel('Max xcorr');
        set(gca,'YDir','reverse');
        legend({'SUA','MUA'});
        
        drawnow;        
        
        % Clear variables for next experiment
        clearvars -except animals curr_animal animal ...
            protocol flexible_name experiments curr_exp ...
            data
               
    end
    
end










