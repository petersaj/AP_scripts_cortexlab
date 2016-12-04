%% Define experiment

animal = 'AP011';
day = '2016-11-07';
experiment = '6';

% ugggghhhhh mpep
mpep_animal = ['M111111_' animal];


%% Load experiment info

% Load timeline
timeline_filename = get_cortexlab_filename(mpep_animal,day,experiment,'timeline');
load(timeline_filename);
timeline_sample_rate = Timeline.hw.daqSampleRate;

% Get camera times
timeline_cam2_idx = strcmp({Timeline.hw.inputs.name}, 'cam2');
cam2_samples = find(Timeline.rawDAQData(1:end-1,timeline_cam2_idx) <= 2 & ...
    Timeline.rawDAQData(2:end,timeline_cam2_idx) > 2);
cam2_time = cam2_samples./timeline_sample_rate;

% Load in protocol
protocol_filename = get_cortexlab_filename(mpep_animal,day,experiment,'protocol');
load(protocol_filename);

% Get stimulus onsets and parameters

photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
thresh = max(Timeline.rawDAQData(:,photodiode_idx))/2;
photodiode_trace = Timeline.rawDAQData(:,photodiode_idx) > thresh;
photodiode_flip = find((~photodiode_trace(1:end-1) & photodiode_trace(2:end)) | ...
    (photodiode_trace(1:end-1) & ~photodiode_trace(2:end)))+1;

photodiode = struct('timestamps',[],'values',[]);
photodiode.timestamps = Timeline.rawDAQTimestamps(photodiode_flip)';
photodiode.values = photodiode_trace(photodiode_flip);

photodiode_onsets = photodiode.timestamps(photodiode.values == 1);

refresh_rate_cutoff = 1/5;
stim_onsets = photodiode_onsets( ...
    [1;find(diff(photodiode_onsets) > refresh_rate_cutoff) + 1]);

stimIDs = zeros(size(stim_onsets));
for q = 1:size(Protocol.seqnums,1)
    stimIDs(Protocol.seqnums(q,:)) = q;
end

%% Load imaging data

data_path = ['\\zserver.cortexlab.net\Data\Subjects\' mpep_animal filesep day];
experiment_path = [data_path filesep num2str(experiment)];

if strmatch(day,'2016-05-31')
    % this day was messed up because of light shield falling off
    frame_t = cam2_time(1:size(V,2))';
else
    frame_t = readNPY([experiment_path filesep 'svdTemporalComponents_cam2.timestamps.npy']);
    
end

Fs = 1./nanmedian(diff(frame_t));

U = readUfromNPY([data_path filesep 'svdSpatialComponents_cam2.npy']);
V = readVfromNPY([experiment_path filesep 'svdTemporalComponents_cam2.npy']);
fV = detrendAndFilt(V, Fs);

avg_im = readNPY([data_path filesep 'meanImage_cam2.npy']);

%% Load ephys data

data_path = ['\\basket.cortexlab.net\data\ajpeters\' animal filesep day filesep 'ephys' filesep num2str(experiment)];

% Load sync/photodiode
load(([data_path filesep 'sync.mat']));

% Read header information
header_path = [data_path filesep 'dat_params.txt'];
header_fid = fopen(header_path);
header_info = textscan(header_fid,'%s %s', 'delimiter',{' = '});
fclose(header_fid);

header = struct;
for i = 1:length(header_info{1})
    header.(header_info{1}{i}) = header_info{2}{i};
end

% Load spike data 
ephys_sample_rate = str2num(header.sample_rate);
spike_times = double(readNPY([data_path filesep 'spike_times.npy']))./ephys_sample_rate;
spike_templates = readNPY([data_path filesep 'spike_templates.npy']);
templates = readNPY([data_path filesep 'templates.npy']);
channel_positions = readNPY([data_path filesep 'channel_positions.npy']);
channel_map = readNPY([data_path filesep 'channel_map.npy']);
winv = readNPY([data_path filesep 'whitening_mat_inv.npy']);
template_amplitudes = readNPY([data_path filesep 'amplitudes.npy']);

% Get stim onset times
% Check that sync matches photodiode number
if length(sync.timestamps) == length(photodiode.timestamps)
    refresh_rate_cutoff = 1/10;
    stim_onset_idx_ephys = [1,find(diff(sync.timestamps) > refresh_rate_cutoff) + 1];
    stim_onset_idx_timeline = [1;find(diff(photodiode.timestamps) > refresh_rate_cutoff) + 1];
else
    warning(['Ephys vs. Timeline photodiode = ' num2str(length(sync.timestamps) - length(photodiode.timestamps))]);
    
    refresh_rate_cutoff = 1/10;
    stim_onset_idx_ephys = [1,find(diff(sync.timestamps) > refresh_rate_cutoff) + 1];
    stim_onset_idx_timeline = [1;find(diff(photodiode.timestamps) > refresh_rate_cutoff) + 1];
    if length(stim_onset_idx_ephys) == length(stim_onset_idx_timeline)
        stim_onset_idx = stim_onset_idx_ephys;
        warning('But same number of stimuli');
    else
        error('Ephys vs. Timeline photodiode: different number of stimuli')
    end
end

% Get the spike times in timeline time
%spike_times_timeline = AP_clock_fix(spike_times,sync.timestamps,photodiode.timestamps);
spike_times_timeline = AP_clock_fix(spike_times,sync.timestamps(stim_onset_idx_ephys),photodiode.timestamps(stim_onset_idx_timeline));

% PSTHs
raster_window = [-1,3];
psthViewer(spike_times_timeline,spike_templates, ...
    photodiode.timestamps(stim_onset_idx_timeline),raster_window, ...
    ones(size(stim_onset_idx_timeline)));

psth_bin_size = 0.001;
[psth,bins,rasterX,rasterY,spikeCounts] = psthAndBA( ...
    spike_times_timeline, ...
    photodiode.timestamps(stim_onset_idx_timeline), ...
    raster_window, psth_bin_size);

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
psth_smooth = conv2(psth, smWin, 'same');
figure; hold on;
plot(bins(20:end-20),psth_smooth(20:end-20)','k','linewidth',2);
line([0,0],ylim,'linestyle','--','color','k');
ylabel('Population spikes');
xlabel('Time from stim onset')


%% Spike-triggered averaging of widefield (single templates)

curr_template = 134;
surround_times = [-2,2];

framerate = 1./median(diff(cam2_time));
surround_frames = round(surround_times(1)*framerate):round(surround_times(2)*framerate);
sta_t = surround_frames./framerate;

curr_spike_times = spike_times_timeline(spike_templates == curr_template);
curr_spike_times(curr_spike_times + surround_times(1) < frame_t(2) | ...
    curr_spike_times + surround_times(2) > frame_t(end)) = [];

% Discretize spikes into frames and count spikes per frame
frame_edges = [frame_t,frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(curr_spike_times,frame_edges);

% Prepare weighting matrix for each frame
frames_w = repmat(frame_spikes'./sum(frame_spikes),1,length(surround_frames));
for curr_sta_frame = 1:length(surround_frames);
    frames_w(:,curr_sta_frame) = ...
        circshift(frames_w(:,curr_sta_frame),[surround_frames(curr_sta_frame),0]);
end

sta_v = fV*frames_w;
sta_im = svdFrameReconstruct(U,sta_v);

% Normalize the image
sta_im_norm = bsxfun(@rdivide,bsxfun(@minus,sta_im,mean(sta_im,3)),std(sta_im,[],3));

% Draw the movie
AP_image_scroll(sta_im_norm,sta_t);
caxis([-std(sta_im_norm(:))*3,std(sta_im_norm(:))*3])
set(gcf,'Name',num2str(curr_template))



[col,row] = MakeSeparable(reshape(sta_im,[],size(sta_im,3)));

[col,row] = MakeSeparable(sta_v);


%% Spike-triggered averaging of widefield (all templates)

surround_times = [-2,2];

framerate = 1./median(diff(cam2_time));
surround_frames = round(surround_times(1)*framerate):round(surround_times(2)*framerate);
sta_t = surround_frames./framerate;

use_templates = unique(spike_templates);
sta_v_all = zeros(size(U,3),length(surround_frames),length(use_templates));

for curr_template_idx = 1:length(use_templates)
    
    curr_template = use_templates(curr_template_idx);
    
    curr_spike_times = spike_times_timeline(spike_templates == curr_template);
    curr_spike_times(curr_spike_times + surround_times(1) < frame_t(2) | ...
        curr_spike_times + surround_times(2) > frame_t(end)) = [];
    
    if isempty(curr_spike_times)
        continue
    end
    
    % Discretize spikes into frames and count spikes per frame
    frame_edges = [frame_t,frame_t(end)+1/framerate];
    [frame_spikes,~,spike_frames] = histcounts(curr_spike_times,frame_edges);
    
    % Prepare weighting matrix for each frame
    frames_w = repmat(frame_spikes'./sum(frame_spikes),1,length(surround_frames));
    for curr_sta_frame = 1:length(surround_frames);
        frames_w(:,curr_sta_frame) = ...
            circshift(frames_w(:,curr_sta_frame),[surround_frames(curr_sta_frame),0]);
    end
    
    sta_v_all(:,:,curr_template_idx) = fV*frames_w;
    disp(curr_template);
end

% This doesn't really work
% Get PCA of all STA V's
[coeff,score,latent] = princomp(reshape(sta_v_all,[],size(sta_v_all,3)));
sta_v_pca = reshape(score,size(sta_v_all,1),size(sta_v_all,2),size(sta_v_all,3));


% sta_im = svdFrameReconstruct(U,sta_v);
% 
% % Normalize the image
% sta_im_norm = bsxfun(@rdivide,bsxfun(@minus,sta_im,mean(sta_im,3)),std(sta_im,[],3));
% 
% % Draw the movie
% AP_image_scroll(sta_im_norm,sta_t);
% caxis([-std(sta_im_norm(:))*3,std(sta_im_norm(:))*3])
% set(gcf,'Name',num2str(curr_template));

%% STA for multiunit

use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 900 & templateDepths < Inf)-1));
%use_spikes = spike_times_timeline(ismember(spike_templates,205));

%use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 0 & templateDepths < 400)-1) & ...
%    ismember(spike_templates,use_templates(use_template_narrow))-1);

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);

surround_times = [-5,5];
framerate = 1./median(diff(frame_t));
surround_frames = round(surround_times(1)*framerate):round(surround_times(2)*framerate);
sta_t = surround_frames./framerate;

sta_im = zeros(size(U,1),size(U,2),length(surround_frames));

% Prepare weighting matrix for each frame
frames_w = repmat(frame_spikes'./sum(frame_spikes),1,length(surround_frames));
for curr_sta_frame = 1:length(surround_frames);
    frames_w(:,curr_sta_frame) = ...
        circshift(frames_w(:,curr_sta_frame),[surround_frames(curr_sta_frame),0]);
end

sta_v = fV*frames_w;
sta_im= svdFrameReconstruct(U,sta_v);

% Normalize the image
% sta_im_norm = bsxfun(@rdivide,bsxfun(@minus,sta_im,mean(sta_im,3)),std(sta_im,[],3) + ...
%     prctile(reshape(std(sta_im,[],3),[],1),50));

% Draw the movie
AP_image_scroll(sta_im,sta_t);

%% STA for templates
template_sta = zeros(size(U,1),size(U,2),length(good_templates));

for curr_template_idx = 1:length(good_templates)
    
    curr_template = good_templates(curr_template_idx);
    
    use_spikes = spike_times_timeline(spike_templates == curr_template);
        
    frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
    [frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
    
    surround_times = [0,0.5];
    framerate = 1./median(diff(frame_t));
    surround_frames = round(surround_times(1)*framerate):round(surround_times(2)*framerate);
    sta_t = surround_frames./framerate;
        
    % Prepare weighting matrix for each frame
    frames_w = repmat(frame_spikes'./sum(frame_spikes),1,length(surround_frames));
    for curr_sta_frame = 1:length(surround_frames);
        frames_w(:,curr_sta_frame) = ...
            circshift(frames_w(:,curr_sta_frame),[surround_frames(curr_sta_frame),0]);
    end
    
    sta_v = fV*frames_w;
    sta_im = svdFrameReconstruct(U,sta_v);
    
    template_sta(:,:,curr_template_idx) = max(sta_im,[],3);
     
    disp(curr_template_idx)
end

% Rearrange STAs by depth
[~,sort_idx] = sort(templateDepths(good_templates+1));
template_sta = template_sta(:,:,sort_idx);

% Plot
AP_image_scroll(template_sta,good_templates(sort_idx));


%% STA for templates by clustered area

use_templates = good_templates;

template_sta = zeros(size(kgrp,1),size(kgrp,2),length(use_templates));
cluster_sta = zeros(size(cluster_trace,1),length(use_templates));

for curr_template_idx = 1:length(use_templates)
    
    curr_template = good_templates(curr_template_idx);
    
    use_spikes = spike_times_timeline(spike_templates == curr_template);
        
    frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
    [frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
    
    surround_times = [0,0.5];
    framerate = 1./median(diff(frame_t));
    surround_frames = round(surround_times(1)*framerate):round(surround_times(2)*framerate);
    sta_t = surround_frames./framerate;
        
    % Prepare weighting matrix for each frame
    frames_w = repmat(frame_spikes'./sum(frame_spikes),1,length(surround_frames));
    for curr_sta_frame = 1:length(surround_frames);
        frames_w(:,curr_sta_frame) = ...
            circshift(frames_w(:,curr_sta_frame),[surround_frames(curr_sta_frame),0]);
    end
    
    sta_cluster = cluster_trace*frames_w;
    sta_cluster_max = max(sta_cluster,[],2);
    cluster_sta(:,curr_template_idx) = sta_cluster_max;
    
    sta_im = zeros(size(kgrp));
    for curr_cluster = 1:size(cluster_trace,1)
        sta_im(kgrp == curr_cluster) = sta_cluster_max(curr_cluster);
    end
    template_sta(:,:,curr_template_idx) = sta_im;
    
    disp(curr_template_idx)
end

% Rearrange STAs by depth
[~,sort_idx] = sort(templateDepths(use_templates+1));
template_sta = template_sta(:,:,sort_idx);

% Plot
AP_image_scroll(template_sta,use_templates(sort_idx));

% PCA results
cluster_sta_nonan = cluster_sta;
cluster_sta_nonan(isnan(cluster_sta_nonan)) = 0;
cluster_sta_nonan_zscore = zscore(cluster_sta_nonan,[],1);
[coeff,score,latent] = pca(cluster_sta_nonan_zscore');
%[coeff,score,latent] = pca(cluster_sta_nonan');

% Plot PCs
figure; colormap(gray);
n_subplots = ceil(sqrt(size(cluster_trace,1)));
for curr_pc = 1:size(coeff,2)
    curr_pc_map = zeros(size(kgrp,1),size(kgrp,2));
    for curr_cluster = 1:size(cluster_trace,1)
        curr_pc_map(kgrp == curr_cluster) = coeff(curr_cluster,curr_pc);
    end
    subplot(n_subplots,n_subplots,curr_pc);
    imagesc(curr_pc_map);
    axis off
    title(['PC ' num2str(curr_pc)]);
    colorbar
end

% Plot PC scores
figure;
scatter3(score(:,1),score(:,2),score(:,3),50,templateDepths(use_templates+1),'filled');
xlabel('PC1');ylabel('PC2');zlabel('PC3');
c = colorbar;
ylabel(c,'Depth (\mum)');
axis square

% Plot first three cluster STAs
figure;
scatter3(cluster_sta_nonan(1,:),cluster_sta_nonan(2,:), ...
    cluster_sta_nonan(3,:),50,templateDepths(use_templates+1),'filled');
xlabel('PC1');ylabel('PC2');zlabel('PC3');
c = colorbar;
ylabel(c,'Depth (\mum)');
axis square

% Plot template by depth and colored by PC score
x_vals = rand(size(use_templates));
figure;scatter(x_vals,templateDepths(use_templates+1),50,'.k')

%% Correlation of templates with clustered area traces

use_templates = good_templates;

template_corr = zeros(size(kgrp,1),size(kgrp,2),length(use_templates));
cluster_corr = zeros(size(cluster_trace,1),length(use_templates));

for curr_template_idx = 1:length(use_templates)
    
    curr_template = use_templates(curr_template_idx);
    
    use_spikes = spike_times_timeline(spike_templates == curr_template);
    
    frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
    [frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
    
%     % Smooth spikes
%     smooth_time = round(35*1);
%     smooth_filt = ones(1,smooth_time)./smooth_time;
%     frame_spikes = conv2(frame_spikes,smooth_filt,'same');
    
    corr_lags = 5;%35*1;
    cluster_xcorr = nan(size(cluster_trace,1),corr_lags*2+1);
        
    % Max cross correlation
    for curr_cluster = 1:size(cluster_trace,1)
        cluster_xcorr(curr_cluster,:) = xcorr(cluster_trace(curr_cluster,:),frame_spikes,corr_lags,'coeff');
    end
    
    cluster_max_xcorr = zeros(size(kgrp));
    for curr_cluster = 1:size(cluster_trace,1)
        cluster_max_xcorr(kgrp == curr_cluster) = max(cluster_xcorr(curr_cluster,:));
    end
    
    cluster_corr(:,curr_template_idx) = max(cluster_xcorr,[],2);
    template_corr(:,:,curr_template_idx) = cluster_max_xcorr;
    
    disp(curr_template_idx);
end

% Rearrange STAs by depth
[~,sort_idx] = sort(templateDepths(use_templates+1));
template_corr = template_corr(:,:,sort_idx);

% Plot correlation by area
template_corr(isnan(template_corr)) = 0;
AP_image_scroll(template_corr,use_templates(sort_idx));
caxis([-0.2 0.2]);
colormap(redblue)

figure;

% Plot vs depth
subplot(1,2,1);
scatter3(cluster_corr(1,:),cluster_corr(2,:), ...
    cluster_corr(3,:),50,templateDepths(use_templates+1),'filled');
xlabel('ROI 1');ylabel('ROI 2');zlabel('ROI 3');
c = colorbar;
ylabel(c,'Depth (\mum)');
axis square

% Plot vs waveform duration
subplot(1,2,2);
scatter3(cluster_corr(1,:),cluster_corr(2,:), ...
    cluster_corr(3,:),50,templateDuration_us(use_templates+1),'filled');
xlabel('ROI 1');ylabel('ROI 2');zlabel('ROI 3');
c = colorbar;
ylabel(c,'Waveform duration (\mus)');
axis square

% Plot all depth vs. ROI separately
figure;
subplot(1,3,1)
plot(cluster_corr(1,:),templateDepths(use_templates+1),'.k','MarkerSize',10)
set(gca,'YDir','reverse');
ylabel('Depth (\mum)');
xlabel('Correlation with fluorescence')
title('ROI 1')
subplot(1,3,2)
plot(cluster_corr(2,:),templateDepths(use_templates+1),'.k','MarkerSize',10)
set(gca,'YDir','reverse');
ylabel('Depth (\mum)');
xlabel('Correlation with fluorescence')
title('ROI 2')
subplot(1,3,3)
plot(cluster_corr(3,:),templateDepths(use_templates+1),'.k','MarkerSize',10)
set(gca,'YDir','reverse');
ylabel('Depth (\mum)');
xlabel('Correlation with fluorescence')
title('ROI 3')


%% Group multiunit by depth, get STAs

% Group by depth
n_depth_groups = 6;
depth_group_edges = linspace(0,max(templateDepths),n_depth_groups+1);
depth_group_edges(end) = Inf;

[depth_group_n,depth_group] = histc(spikeDepths,depth_group_edges);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges(1:end-1)+diff(depth_group_edges);

surround_times = [-2,2];
framerate = 1./median(diff(frame_t));
surround_frames = round(surround_times(1)*framerate):round(surround_times(2)*framerate);
sta_t = surround_frames./framerate;

sta_im = zeros(size(U,1),size(U,2),length(surround_frames),length(unique(depth_group)));

for curr_depth = unique(depth_group)'
    
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    curr_spike_times(curr_spike_times + surround_times(1) < frame_t(2) | ...
        curr_spike_times + surround_times(2) > frame_t(end)) = [];
    
    % Discretize spikes into frames and count spikes per frame
    frame_edges = [frame_t,frame_t(end)+1/framerate];
    [frame_spikes,~,spike_frames] = histcounts(curr_spike_times,frame_edges);
    
    % Prepare weighting matrix for each frame
    frames_w = repmat(frame_spikes'./sum(frame_spikes),1,length(surround_frames));
    for curr_sta_frame = 1:length(surround_frames);
        frames_w(:,curr_sta_frame) = ...
            circshift(frames_w(:,curr_sta_frame),[surround_frames(curr_sta_frame),0]);
    end
    
    sta_v = fV*frames_w;
    sta_im(:,:,:,curr_depth) = svdFrameReconstruct(U,sta_v);

    disp(curr_depth)
    
end

% Normalize the image
% sta_im_norm = bsxfun(@rdivide,bsxfun(@minus,sta_im,mean(sta_im,3)),std(sta_im,[],3) + ...
%     prctile(reshape(std(sta_im,[],3),[],1),50));

sta_im_norm = bsxfun(@rdivide,bsxfun(@minus,sta_im,px_mean),px_std);

% Draw the movie
AP_image_scroll(sta_im,sta_t);

% Plot ROI traces across depth (after drawing an ROI)
trace_norm = bsxfun(@minus,roi.trace,min(roi.trace,[],2));
trace_norm = bsxfun(@rdivide,trace_norm,max(trace_norm,[],2));

figure; hold on
set(gca,'ColorOrder',copper(n_depth_groups));
plot(sta_t,roi.trace','linewidth',2);
%plot(sta_t,trace_norm','linewidth',2);
legend(arrayfun(@(x) ['Depth: ' num2str(round(depth_group_centers(x))) ...
    ' \mum'],1:n_depth_groups,'uni',false),'location','nw');
xlabel('Time from spike (s)');
ylabel('Fluorescence');


%% Group multiunit by waveform duration, get STAs

% Group waveforms by waveforms categorized elsewhere
waveform_groups = ismember(spike_templates,use_templates(use_template_narrow))+1;

surround_times = [-2,2];
framerate = 1./median(diff(frame_t));
surround_frames = round(surround_times(1)*framerate):round(surround_times(2)*framerate);
sta_t = surround_frames./framerate;

sta_im = zeros(size(U,1),size(U,2),length(surround_frames),2);

for curr_group = 1:2
    
    curr_spike_times = spike_times_timeline(waveform_groups == curr_group);
    curr_spike_times(curr_spike_times + surround_times(1) < frame_t(2) | ...
        curr_spike_times + surround_times(2) > frame_t(end)) = [];
    
    % Discretize spikes into frames and count spikes per frame
    frame_edges = [frame_t,frame_t(end)+1/framerate];
    [frame_spikes,~,spike_frames] = histcounts(curr_spike_times,frame_edges);
    
    % Prepare weighting matrix for each frame
    frames_w = repmat(frame_spikes'./sum(frame_spikes),1,length(surround_frames));
    for curr_sta_frame = 1:length(surround_frames);
        frames_w(:,curr_sta_frame) = ...
            circshift(frames_w(:,curr_sta_frame),[surround_frames(curr_sta_frame),0]);
    end
    
    sta_v = fV*frames_w;
    sta_im(:,:,:,curr_group) = svdFrameReconstruct(U,sta_v);

    disp(curr_group)
    
end

% Normalize the image
sta_im_norm = zscore(sta_im,[],3);

% Draw the movie
AP_image_scroll(sta_im_norm,sta_t);
caxis([-std(sta_im_norm(:))*3,std(sta_im_norm(:))*3])

% Plot ROI traces across depth (after drawing an ROI)
trace_norm = bsxfun(@minus,roi.trace,min(roi.trace,[],2));
trace_norm = bsxfun(@rdivide,trace_norm,max(trace_norm,[],2));

figure; hold on
set(gca,'ColorOrder',copper(2));
plot(sta_t,roi.trace','linewidth',2);
%plot(sta_t,trace_norm','linewidth',2);
% is this right?
probe_end = 1300;
legend({'Short waveforms','Long waveforms'});
xlabel('Time from spike (s)');
ylabel('Fluorescence');


%% Get fluorescence traces/multiunit crosscorr

% Choose ROI
h = figure;
imagesc(avg_im);
colormap(gray);
caxis([0 prctile(avg_im(:),90)]);
roiMask = roipoly;
close(h);

% Get fluorescence across session in ROI
U_roi = reshape(U(repmat(roiMask,1,1,size(U,3))),[],size(U,3));
roi_trace = nanmean(U_roi*fV);

% Get population spikes per frame
use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 1000 & templateDepths < Inf)-1));

framerate = 1./nanmedian(diff(frame_t));
frame_edges = [frame_t,frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);

figure;hold on;
plot(zscore(roi_trace),'k','linewidth',1);
plot(zscore(frame_spikes),'r','linewidth',1);
xlabel('Frames');
legend({'Fluorescence','Spikes'});

corr_lags = 100;
[~,lags] = xcorr(ones(size(frame_spikes)),corr_lags);
lags_t = lags./framerate;
figure;
subplot(2,1,1); hold on;
plot(lags_t,xcov(frame_spikes,corr_lags,'coeff'),'k','linewidth',2)
plot(lags_t,xcov(roi_trace,corr_lags,'coeff'),'r','linewidth',2)
line([0,0],ylim,'linestyle','--','color','k');
legend({'Spikes','Fluorescence'});
xlabel('Lag (s)')
ylabel('Autocorrelation')
subplot(2,1,2);
plot(lags_t,xcov(roi_trace,frame_spikes,corr_lags),'k','linewidth',2)
line([0,0],ylim,'linestyle','--','color','k');
xlabel('Lag (s)');
ylabel('Cross-correlation')

%% Get fluorescence traces with single unit

n_rois = 1;

% Choose ROI
h = figure;
imagesc(avg_im);
set(gca,'YDir','normal');
colormap(gray);
caxis([0 prctile(avg_im(:),90)]);
roiMask = false(size(avg_im,1),size(avg_im,2),n_rois);
for i = 1:n_rois
    roiMask(:,:,i) = roipoly;
end
close(h);

% Get fluorescence across session in ROI
roi_trace = nan(n_rois,size(fV,2));
for i = 1:n_rois
    U_roi = reshape(U(repmat(roiMask(:,:,i),1,1,size(U,3))),[],size(U,3));
    roi_trace(i,:) = nanmean(U_roi*fVn);
end

% Get population spikes per frame
framerate = 1./nanmedian(diff(frame_t));
frame_edges = [frame_t,frame_t(end)+1/framerate];

frame_spikes = zeros(max(spike_templates),length(frame_t));
for curr_template_idx = 1:max(spike_templates)+1;
    frame_spikes(curr_template_idx,:) = histcounts(spike_times_timeline(spike_templates == curr_template_idx-1),frame_edges);
end

spike_roi_corr_all = corrcoef([roi_trace(1,:)',frame_spikes']);
spike_roi_corr = spike_roi_corr_all(1,2:end);
spike_roi_corr(isnan(spike_roi_corr)) = 0;

[~,sort_idx] = sort(spike_roi_corr,'descend');
use_templates = [sort_idx(1:20),sort_idx(end-19:end)]-1;
% 
% figure; hold on
% for i = 1:length(use_templates)
%     curr_spikes = spike_times_timeline(spike_templates == use_templates(i));
%     plot(curr_spikes,i,'.k');
% end

figure; 
p1 = subplot(4,1,1); hold on;

roi_trace_plot = bsxfun(@plus,zscore(roi_trace,[],2),5*transpose(1:size(roi_trace,1)));
plot(roi_trace_plot','k','linewidth',2);

p2 = subplot(4,1,2:4);
imagesc(1-frame_spikes(use_templates+1,:));
colormap(gray);caxis([0 1]);
linkaxes([p1,p2],'x');

%% Figure for Kenneth

%use_time = [393,396];
use_time = [380,400];
%use_time = [380,434];

n_rois = 8;

% Choose ROI
h = figure;
imagesc(avg_im);
set(gca,'YDir','normal');
colormap(gray);
caxis([0 10000]);
roiMask = false(size(avg_im,1),size(avg_im,2),n_rois);
for i = 1:n_rois
    roiMask(:,:,i) = roipoly;
end
close(h);

% Get fluorescence across session in ROI
use_frames = frame_t > use_time(1) & frame_t < use_time(2);

roi_trace = nan(n_rois,sum(use_frames));
for i = 1:n_rois
    U_roi = reshape(U(repmat(roiMask(:,:,i),1,1,size(U,3))),[],size(U,3));
    roi_trace(i,:) = nanmean(U_roi*fV(:,use_frames));
end

roi_trace_ipsi = roi_trace(1:2:end,:);
roi_trace_contra = roi_trace(2:2:end,:);


% Raster by depth (color coded by nucleus)

unique_depths = sort(unique(templateDepths));

nucleus_colors = lines(5);
nucleus_borders = [0,220,360,500,1150,max(templateDepths)];
depth_nucleus = discretize(unique_depths,nucleus_borders);

bin_edges = use_time(1):0.001:use_time(end);
bin_centers = bin_edges(1:end-1) + (diff(bin_edges)/2);
raster_x_time = cell(length(unique_depths),1);
raster_y = cell(length(unique_depths),1);
binned_spikes_depth = zeros(length(unique_depths),length(corr_edges)-1);
for curr_depth = 1:length(unique_depths);
    curr_binned_spikes = histcounts(spike_times_timeline( ...
        ismember(spike_templates,find(templateDepths == unique_depths(curr_depth))-1)), ...
        bin_edges);
    
    [raster_x,raster_y{curr_depth}] = ...
        rasterize(find(curr_binned_spikes));
    
    raster_x(~isnan(raster_x)) = bin_centers(raster_x(~isnan(raster_x)));
    raster_x_time{curr_depth} = raster_x - use_time(1);
end

raster_height = 10;
raster_y_depth = cellfun(@(raster_y,depth) (raster_y*raster_height) + ...
    depth,raster_y,num2cell(unique_depths),'uni',false);

figure; hold on
for curr_depth = 1:length(unique_depths);
    curr_color = nucleus_colors(depth_nucleus(curr_depth),:);
    plot(raster_x_time{curr_depth},raster_y_depth{curr_depth},'color',curr_color);    
end

set(gca,'YDir','Reverse');


% Smoothed population rate by nucleus

nucleus_borders = [0,220,360,500,1150,max(unique_depths)];
template_nucleus = discretize(templateDepths,nucleus_borders);

bin_edges = use_time(1):0.001:use_time(end);
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
binned_spikes_nucleus = zeros(length(nucleus_borders)-1,length(bin_edges)-1);
for curr_nucleus = 1:length(nucleus_borders)-1;
    binned_spikes_nucleus(curr_nucleus,:) = histcounts(spike_times_timeline( ...
        ismember(spike_templates,find(template_nucleus == curr_nucleus)-1)), ...
        bin_edges);
end

smooth_size = 100;
gw = gausswin(round(smooth_size*6),3)';
smWin = gw./sum(gw);
smoothed_spikes_nucleus = conv2(binned_spikes_nucleus, smWin, 'same');

% Plot fluorescence and population rates in nuclei together

figure; 

p1 = subplot(2,1,1); hold on;

trace_spacing = 5*transpose(1:size(roi_trace_ipsi,1));

roi_trace_ipsi_plot = bsxfun(@plus,zscore(roi_trace_ipsi,[],2),trace_spacing);
h1 = plot(frame_t(use_frames),roi_trace_ipsi_plot','k','linewidth',1);

roi_trace_contra_plot = bsxfun(@plus,zscore(roi_trace_contra,[],2),trace_spacing);
h2 = plot(frame_t(use_frames),roi_trace_contra_plot','r','linewidth',1);

set(gca,'YTick',trace_spacing);
set(gca,'YTickLabel',{'Retrosplenial','Visual','Barrel','Motor'});
legend([h1(1),h2(1)],'Ipsi','Contra')

p2 = subplot(2,1,2);

trace_spacing = fliplr(10*(1:size(binned_spikes_nucleus,1)));
plot(bin_centers-bin_centers(1),bsxfun(@plus,zscore(smoothed_spikes_nucleus,[],2)',trace_spacing),'linewidth',2);
set(gca,'YTick',sort(trace_spacing));
set(gca,'YTickLabel',fliplr({'MD','CM','AV','VL','VPL'}));

linkaxes([p1,p2],'x');

%%%%%%%%%%%%%%%%%%%%


%%%%%% HOLDING FOR LATER?
% Define bad templates, manual for now based on PCA
[coeff,score,latent] = pca(waveforms');
[~,sort_idx] = sort(coeff(:,1),'descend');
good_templates = sort_idx(1:400);


% Temporary: autocorrelation of VL
use_depths = [600,900];

bin_edges = 0:0.001:frame_t(end);

for curr_depth = 1:length(unique_depths);
    vl_spikes = histcounts(spike_times_timeline( ...
        ismember(spike_templates,find(templateDepths > use_depths(1) & templateDepths < use_depths(2))-1)), ...
        bin_edges);
end
figure;
[vl_autocorr,lags] = xcorr(vl_spikes,100);
plot(lags,vl_autocorr);
xlabel('Lag (ms)');
ylabel('Autocorrelation')


%%%%%% Correlation plot for delineating areas


% Correlation plot 

corr_edges = 0:0.01:frame_t(end);

binned_spikes = zeros(max(spike_templates),length(corr_edges)-1);
for curr_template_idx = 1:max(spike_templates)+1;
    binned_spikes(curr_template_idx,:) = histcounts(spike_times_timeline(spike_templates == curr_template_idx-1),corr_edges);
end

[~,sort_idx] = sort(templateDepths,'ascend');
binned_spikes_depthsort = binned_spikes(sort_idx,:);

figure;imagesc(corrcoef(binned_spikes_depthsort'));
colormap(hot);

% Nick suggestion: correlation plot by depth MUA

unique_depths = sort(unique(templateDepths));

n_depth_groups = 30;
depth_group_edges = linspace(min(templateDepths),max(templateDepths),n_depth_groups+1);
depth_group = discretize(templateDepths,depth_group_edges);
depth_group_centers = grpstats(templateDepths,depth_group);
unique_depths = unique(depth_group);

corr_edges = 0:0.01:frame_t(end);

binned_spikes_depth = zeros(length(unique_depths),length(corr_edges)-1);
for curr_depth = 1:length(unique_depths);
    binned_spikes_depth(curr_depth,:) = histcounts(spike_times_timeline( ...
        ismember(spike_templates,find(depth_group == unique_depths(curr_depth))-1)), ...
        corr_edges);
end

figure;imagesc(depth_group_centers,depth_group_centers,corrcoef(binned_spikes_depth'));
colormap(hot);


%% TEMPORARY KENNETH FIGURE: STA from thalamic nuclei

unique_depths = sort(unique(templateDepths));

nucleus_borders = [0,220,360,500,1150,max(unique_depths)];
template_nucleus = discretize(templateDepths,nucleus_borders);
n_nuclei = length(nucleus_borders) - 1;

bin_edges = [frame_t,frame_t(end)+(1/framerate)];%frame_t(1):0.01:frame_t(end);
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
binned_spikes_nucleus = zeros(length(nucleus_borders)-1,length(bin_edges)-1);
for curr_nucleus = 1:length(nucleus_borders)-1;
    binned_spikes_nucleus(curr_nucleus,:) = histcounts(spike_times_timeline( ...
        ismember(spike_templates,find(template_nucleus == curr_nucleus)-1)), ...
        bin_edges);
end

corr_lags = 100;
v_xcorr = nan(size(V,1),corr_lags*2+1,n_nuclei);
for curr_nucleus = 1:n_nuclei
    for curr_u = 1:size(U,3)
        v_xcorr(curr_u,:,curr_nucleus) = xcorr(fV(curr_u,:),binned_spikes_nucleus(curr_nucleus,:),corr_lags);
    end
end

lags_t = (-corr_lags:corr_lags)/framerate;

svd_xcorr = nan(size(U,1),size(U,2),corr_lags*2+1,n_nuclei);
for curr_nucleus = 1:n_nuclei
    svd_xcorr(:,:,:,curr_nucleus) = svdFrameReconstruct(U,v_xcorr(:,:,curr_nucleus));
end

% Normalize the image
svd_xcorr_norm = bsxfun(@rdivide,bsxfun(@minus,svd_xcorr,mean(svd_xcorr,3)),std(svd_xcorr,[],3) + ...
    prctile(reshape(std(svd_xcorr,[],3),[],1),50));

% Draw the movie
AP_image_scroll(svd_xcorr,lags_t);

% Plot ROI traces across depth (after drawing an ROI)
trace_spacing = fliplr(5*(1:size(roi.trace,1)));
roi_trace_plot = bsxfun(@plus,zscore(roi.trace,[],2)',trace_spacing);

figure; hold on
set(gca,'ColorOrder',lines(n_nuclei));
plot(lags_t,roi_trace_plot,'linewidth',2);
xlabel('Time from spike (s)');
ylabel('Fluorescence');
set(gca,'YTick',sort(trace_spacing));
set(gca,'YTickLabel',fliplr({'Habenula?','CM','AV','VL','VPL?'}));

%%%%%%%%%%%%%% proper STA

surround_times = [-2,2];
framerate = 1./median(diff(cam2_time));
surround_frames = round(surround_times(1)*framerate):round(surround_times(2)*framerate);
sta_t = surround_frames./framerate;

sta_im = nan(size(U,1),size(U,2),length(surround_frames),n_nuclei);

for curr_nucleus = 1:n_nuclei
    
    % Prepare weighting matrix for each frame
    frame_spikes = binned_spikes_nucleus(curr_nucleus,:);
    frames_w = repmat(frame_spikes'./sum(frame_spikes),1,length(surround_frames));
    for curr_sta_frame = 1:length(surround_frames);
        frames_w(:,curr_sta_frame) = ...
            circshift(frames_w(:,curr_sta_frame),[surround_frames(curr_sta_frame),0]);
    end
    
    sta_v = fV(:,1:end-1)*frames_w;
    sta_im(:,:,:,curr_nucleus) = svdFrameReconstruct(U,sta_v);
    
end

% Normalize the image
sta_im_norm = bsxfun(@rdivide,bsxfun(@minus,sta_im,mean(sta_im,3)),std(sta_im,[],3) + ...
    prctile(reshape(std(sta_im,[],3),[],1),50));

% Draw the movie
AP_image_scroll(sta_im_norm,sta_t);
caxis([-std(sta_im_norm(:))*3,std(sta_im_norm(:))*3])

%%%%%%%%%%%%% average image triggered on all
surround_times = [-2,2];
framerate = 1./median(diff(cam2_time));
surround_frames = round(surround_times(1)*framerate):round(surround_times(2)*framerate);
sta_t = surround_frames./framerate;

% Prepare weighting matrix for each frame
frame_spikes = sum(binned_spikes_nucleus(curr_nucleus,:),1);
frames_w = repmat(frame_spikes'./sum(frame_spikes),1,length(surround_frames));
for curr_sta_frame = 1:length(surround_frames);
    frames_w(:,curr_sta_frame) = ...
        circshift(frames_w(:,curr_sta_frame),[surround_frames(curr_sta_frame),0]);
end

sta_v = fV(:,1:end-1)*frames_w;
sta_im = svdFrameReconstruct(U,sta_v);

% Normalize the image
sta_im_norm = bsxfun(@rdivide,bsxfun(@minus,sta_im,mean(sta_im,3)),std(sta_im,[],3) + ...
    prctile(reshape(std(sta_im,[],3),[],1),50));

% Draw ROIs over image
figure; hold on;
imagesc(nanmean(sta_im_norm,3));

for curr_roi = 1:size(roiMask,3);
    first_nonzero = find(roiMask(:,:,curr_roi) > 0,1);
    [y_nonzero x_nonzero] = ind2sub([size(U,1),size(U,2)],first_nonzero);
    roi_perim = bwtraceboundary(roiMask(:,:,curr_roi),[y_nonzero x_nonzero],'N');
    plot(roi_perim(:,2),roi_perim(:,1),'b','linewidth',2);
end
    

%% Get fluorescence traces/multiunit crosscorr in parts within recording

% Choose ROI
h = figure;
imagesc(avg_im);
set(gca,'YDir','normal');
colormap(gray);
caxis([0 10000]);
roiMask = roipoly;
close(h);

% Get fluorescence across session in ROI
U_roi = reshape(U(repmat(roiMask,1,1,size(U,3))),[],size(U,3));
roi_trace = nanmean(U_roi*fV);

% Get population spikes per frame
framerate = 1./nanmedian(diff(frame_t));
frame_edges = [frame_t,frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(spike_times_timeline,frame_edges);

n_split = 4;

split_frame_edges = linspace(1,size(fV,2),n_split+1);
split_frames = discretize(1:size(fV,2),split_frame_edges);

corr_lags = 100;
split_xcorr = nan(n_split,corr_lags*2+1);
for curr_split = 1:n_split;
    roi_trace_split = roi_trace(split_frames == curr_split);
    frame_spikes_split = frame_spikes(:,split_frames == curr_split);    
    [split_xcorr(curr_split,:),lags] = xcov(roi_trace_split,frame_spikes_split,corr_lags);   
end

figure; hold on
plot(lags./framerate,split_xcorr','linewidth',2)
xlabel('Lag (s)')
ylabel('Cross-covariance');
legend(cellfun(@(x) ['Part ' num2str(x)],num2cell(1:n_split),'uni',false));


%% Get kernel between spikes and fluorescence

% Skip the first n seconds to do this
skip_seconds = 30;
use_frames = frame_t > skip_seconds;

% Choose ROI
h = figure;
imagesc(avg_im);
set(gca,'YDir','reverse');
colormap(gray);
caxis([0 prctile(avg_im(:),90)]);
title('Pick ROI to define kernel')
drawnow
roiMask = roipoly;
close(h);

% Get fluorescence across session in ROI
U_roi = reshape(U(repmat(roiMask,1,1,size(U,3))),[],size(U,3));
roi_trace_full = nanmean(U_roi*fV);
roi_trace = roi_trace_full(use_frames);

% Get population spikes per frame
framerate = 1./nanmedian(diff(frame_t));
frame_edges = [frame_t,frame_t(end)+1/framerate];

use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 0 & templateDepths < Inf)-1));

[frame_spikes_full,~,spike_frames] = histcounts(use_spikes,frame_edges);
frame_spikes = frame_spikes_full(use_frames);

corr_lags = 100;
[~,lags] = xcorr(ones(size(frame_spikes)),corr_lags);
lags_t = lags./framerate;
figure;
subplot(2,1,1); hold on;
plot(lags_t,xcov(frame_spikes,corr_lags,'coeff'),'k','linewidth',2)
plot(lags_t,xcov(roi_trace,corr_lags,'coeff'),'r','linewidth',2)
line([0,0],ylim,'linestyle','--','color','k');
legend({'Spikes','Fluorescence'});
xlabel('Lag (s)')
ylabel('Autocorrelation')
subplot(2,1,2);
plot(lags_t,xcov(roi_trace,frame_spikes,corr_lags),'k','linewidth',2)
line([0,0],ylim,'linestyle','--','color','k');
xlabel('Lag (s)');
ylabel('Cross-correlation')

% Non-normalized
x_nonorm = ifft((fft(roi_trace).*conj(fft(frame_spikes)))); % unnormalized

% I think this is what Krumin suggested? (not sure if he meant mean)
x = ifft((fft(roi_trace).*conj(fft(frame_spikes)))./(mean(fft(frame_spikes)).*mean(conj(fft(frame_spikes)))));

% This looks like from Nauhaus 2012?
x_autonorm = ifft((fft(roi_trace).*conj(fft(frame_spikes)))./(fft(frame_spikes).*conj(fft(frame_spikes))));

plot_frames = 35*10;

figure;

t_shift = [frame_t(end-plot_frames+1:end)-frame_t(end)-1/framerate,frame_t(1:plot_frames)-frame_t(1)];

p1 = subplot(2,1,1);
plot(t_shift,[x_nonorm(end-plot_frames+1:end),x_nonorm(1:plot_frames)],'k','linewidth',2);
xlabel('Time (s)');
ylabel('Impulse response')
title('Non-normalized: ifft(F*S'')');

p2 = subplot(2,1,2);
plot(t_shift,[x_autonorm(end-plot_frames+1:end),x_autonorm(1:plot_frames)],'k','linewidth',2);
xlabel('Time (s)');
ylabel('Impluse response');
title('Normalized: ifft(F*S''/S*S'')');

linkaxes([p1,p2],'x')

% Use the corrected impulse response for convolving kernel
gcamp_kernel = x_autonorm(1:plot_frames);
frame_spikes_conv_full = conv(frame_spikes_full,gcamp_kernel);
frame_spikes_conv = frame_spikes_conv_full(1:length(frame_spikes_full));

figure;hold on;
plot(frame_t,zscore(roi_trace_full),'k','linewidth',1);
plot(frame_t,zscore(frame_spikes_conv),'b','linewidth',1);
xlabel('Time (s)');
legend({'Fluorescence','Spikes conv'});

%% Make kernel that looks like GCaMP6f

use_t = 1:100;
t = frame_t(use_t);

event_trace = double(x_autonorm(use_t)./max(x_autonorm(use_t)));
starting = [0,0.2,1,0.1,0.5,0.2];
estimates = fminsearch(@(x) AP_fit_gcamp_kernel(x,frame_t(use_t),event_trace),starting);

t0 = estimates(1);
tau_on = estimates(2);
A1 = estimates(3);
tau1 = estimates(4);
A2 = estimates(5);
tau2 = estimates(6);

fitted_curve = (1-exp(-(t-t0)./tau_on)).*(A1*exp(-(t-t0)./tau1) + A2*exp(-(t-t0)./tau2));
fitted_curve(fitted_curve < 0) = 0;

figure; hold on;
plot(frame_t(use_t),event_trace,'k');
plot(frame_t(use_t),fitted_curve,'r');

% Here's a reasonable-ish one
t0 = 2.85;
tau_on = 1.46;
A1 = -866;
tau1 = 0.05;
A2 = 59.87;
tau2 = 0.12;

fitted_curve = (1-exp(-(t-t0)./tau_on)).*(A1*exp(-(t-t0)./tau1) + A2*exp(-(t-t0)./tau2));
fitted_curve(fitted_curve < 0) = 0;

gcamp_kernel = fitted_curve;


%% Get STA-equivalent via xcorr

skip_frames = 35*10;

%use_spikes = spike_times_timeline;
use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 0 & templateDepths < 500)-1));
%use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 0 & templateDepths < 400)-1) & ...
%    (ismember(spike_templates,use_templates(use_template_narrow))));

framerate = 1./nanmedian(diff(frame_t));

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);

corr_lags = 35*1;
v_xcorr = nan(size(fV,1),corr_lags*2+1);

for curr_u = 1:size(U,3)  
    v_xcorr(curr_u,:) = xcov(fV(curr_u,skip_frames:end)-mean(fV(curr_u,skip_frames:end)), ...
        frame_spikes(skip_frames:end) - mean(frame_spikes(skip_frames:end)),corr_lags,'biased');
end

lags_t = (-corr_lags:corr_lags)/framerate;

svd_xcorr = svdFrameReconstruct(U,v_xcorr);

% Normalize the image
svd_xcorr_norm = bsxfun(@rdivide,svd_xcorr,px_std*std(frame_spikes(skip_frames:end)));

% Draw the movie
AP_image_scroll(svd_xcorr_norm,lags_t);



%% Get zero-lag correlation in SVD space

skip_frames = 35*10;

%use_spikes = spike_times_timeline;
use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 0 & templateDepths < 500)-1));
%use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 0 & templateDepths < 400)-1) & ...
%    (ismember(spike_templates,use_templates(use_template_narrow))));

framerate = 1./nanmedian(diff(frame_t));

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);

v_cov = mean(bsxfun(@times,bsxfun(@minus,fV(:,skip_frames:end),mean(fV(:,skip_frames:end),2)), ...
    (frame_spikes(skip_frames:end)-mean(frame_spikes(skip_frames:end)))),2);

svd_corr = bsxfun(@rdivide,svdFrameReconstruct(U,v_cov),px_std.*std(frame_spikes(skip_frames:end)));

figure;imagesc(svd_corr);colormap(gray);




%% For last cell: make baseline std based shifted times
% ??

framerate = 1./nanmedian(diff(frame_t));

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(spike_times_timeline,frame_edges);

frame_spikes = circshift(frame_spikes,[0,randi(length(frame_spikes))]);

corr_lags = 100;
v_xcorr = nan(size(fV,1),corr_lags*2+1);

for curr_u = 1:size(U,3)  
    v_xcorr(curr_u,:) = xcorr(fV(curr_u,:),frame_spikes,corr_lags);
end

lags_t = (-corr_lags:corr_lags)/framerate;

svd_xcorr = svdFrameReconstruct(U,v_xcorr);

svd_std = nanstd(svd_xcorr,[],3);
svd_mean = nanmean(svd_xcorr,3);



%% PCA of activity in population and correlation with fluorescence
% THIS SEEMS MESSED UP DO THIS AGAIN

framerate = 1./median(diff(frame_t));
frame_edges = [frame_t,frame_t(end)+1/framerate];

curr_templates = intersect(find(templateDepths > 0 & templateDepths < Inf)-1,use_templates(msn));

frame_spikes = zeros(length(curr_templates),length(frame_t));
for curr_template_idx = 1:length(curr_templates);
    curr_template = curr_templates(curr_template_idx);
    frame_spikes(curr_template_idx,:) = histcounts( ...
        spike_times_timeline(spike_templates == curr_template),frame_edges);
end

% sort spike population by type of activity
frame_spikes_norm = zscore(frame_spikes,[],2);
[coeff,score,latent] = pca(frame_spikes_norm');
use_metric = coeff(:,1);
[~,sort_idx] = sort(use_metric,'descend');

% Plot sort index by depth
figure;plot(use_metric,templateDepths(curr_templates+1),'.k','MarkerSize',15);
xlabel('PC1 coefficient')
ylabel('Depth (\mum)')
set(gca,'YDir','reverse');

% Choose ROI
h = figure;
imagesc(avg_im);
colormap(gray);
caxis([0 prctile(avg_im(:),90)]);
roiMask = roipoly;
close(h);

% Get fluorescence across session in ROI
U_roi = reshape(U(repmat(roiMask,1,1,size(U,3))),[],size(U,3));
roi_trace = nanmean(U_roi*fV);

frame_spikes_1 = nanmean(frame_spikes(sort_idx(1:10),:),1);
frame_spikes_2 = nanmean(frame_spikes(sort_idx(end-10:end),:),1);
%use_spikes = nanmean(frame_spikes,1);

corr_lags = 35*10;
[~,lags] = xcorr(ones(size(frame_spikes_1)),corr_lags);
lags_t = lags./framerate;
figure;
subplot(2,1,1); hold on;
plot(lags_t,xcov(frame_spikes_1,corr_lags,'coeff'),'k','linewidth',2)
plot(lags_t,xcov(frame_spikes_2,corr_lags,'coeff'),'b','linewidth',2)
plot(lags_t,xcov(roi_trace,corr_lags,'coeff'),'r','linewidth',2)
line([0,0],ylim,'linestyle','--','color','k');
xlabel('Lag (s)')
ylabel('Autocorrelation')
legend({'High PC1 spikes','Low PC1 spikes','Fluorescence'});
subplot(2,1,2); hold on;
plot(lags_t,xcov(roi_trace,frame_spikes_1,corr_lags),'k','linewidth',2)
plot(lags_t,xcov(roi_trace,frame_spikes_2,corr_lags),'b','linewidth',2)
line([0,0],ylim,'linestyle','--','color','k');
xlabel('Lag (s)');
ylabel('Cross-correlation')
legend({'High PC1 spikes','Low PC1 spikes'});


% STA-equivalent via xcorr
v_xcorr_1 = nan(size(fV,1),corr_lags*2+1);
v_xcorr_2 = nan(size(fV,1),corr_lags*2+1);

for curr_u = 1:size(U,3)  
    v_xcorr_1(curr_u,:) = xcorr(fV(curr_u,:),frame_spikes_1,corr_lags);
    v_xcorr_2(curr_u,:) = xcorr(fV(curr_u,:),frame_spikes_2,corr_lags);
end

lags_t = (-corr_lags:corr_lags)/framerate;

svd_xcorr_1 = svdFrameReconstruct(U,v_xcorr_1);
svd_xcorr_2 = svdFrameReconstruct(U,v_xcorr_2);

% Normalize the image
svd_xcorr_norm_1 = bsxfun(@rdivide,bsxfun(@minus,svd_xcorr_1,mean(svd_xcorr_1,3)),std(svd_xcorr_1,[],3) + ...
    prctile(reshape(std(svd_xcorr_1,[],3),[],1),50));

svd_xcorr_norm_2 = bsxfun(@rdivide,bsxfun(@minus,svd_xcorr_2,mean(svd_xcorr_2,3)),std(svd_xcorr_2,[],3) + ...
    prctile(reshape(std(svd_xcorr_2,[],3),[],1),50));

% Draw the movie
AP_image_scroll([svd_xcorr_norm_1,svd_xcorr_norm_2],lags_t);


% Get correlations and PCA by different spike bin widths
bin_time = 0.02; % seconds

corr_edges = min(spike_times_timeline):bin_time:max(spike_times_timeline);
binned_spikes = zeros(length(curr_templates),length(corr_edges)-1);
for curr_template_idx = 1:length(curr_templates);
    curr_spikes = spike_times_timeline(spike_templates == curr_templates(curr_template_idx));
    binned_spikes(curr_template_idx,:) = histcounts(curr_spikes,corr_edges);
end

binned_corr = corrcoef(binned_spikes(sort_idx,:)');
figure;
subplot(5,1,1:4);
imagesc(binned_corr);
axis square
colormap(hot);
caxis([0, max(AP_itril(binned_corr,-2))]);
h = colorbar;
xlabel('Tempate (sorted by PCA)')
ylabel('Tempate (sorted by PCA)')
title(['Bin time: ' num2str(bin_time)]);

% Sort by PCA on differently binned spikes
binned_spikes_norm = zscore(binned_spikes,[],2);
[coeff,score,latent] = pca(binned_spikes_norm');
use_metric = coeff(:,1);
[~,sort_idx] = sort(use_metric,'descend');

% Plot sort index by depth
figure;plot(use_metric,templateDepths(curr_templates+1),'.k','MarkerSize',15);
xlabel('PC1 coefficient')
ylabel('Depth (\mum)')
set(gca,'YDir','reverse');

use_spikes_1 = nanmean(binned_spikes(sort_idx(1:20),:),1);
use_spikes_2 = nanmean(binned_spikes(sort_idx(end-19:end,:),:),1);
[x,lags] = xcov(use_spikes_1,use_spikes_2,3/bin_time);
x(lags == 0) = NaN; % why is this happening that it's huge here?
subplot(5,1,5);
plot(lags*bin_time,x,'k');
ylabel('Cross correlation (top/bottom 20)')
xlabel('Lag group 2 (sec)');

frame_spikes_1 = nanmean(frame_spikes(sort_idx(1:20),:),1);
frame_spikes_2 = nanmean(frame_spikes(sort_idx(end-19:end),:),1);

corr_lags = 100;
[~,lags] = xcorr(ones(size(frame_spikes_1)),corr_lags);
lags_t = lags./framerate;
figure;
subplot(2,1,1); hold on;
plot(lags_t,xcov(frame_spikes_1,corr_lags,'coeff'),'k','linewidth',2)
plot(lags_t,xcov(frame_spikes_2,corr_lags,'coeff'),'b','linewidth',2)
plot(lags_t,xcov(roi_trace,corr_lags,'coeff'),'r','linewidth',2)
line([0,0],ylim,'linestyle','--','color','k');
xlabel('Lag (s)')
ylabel('Autocorrelation')
legend({'High PC1 spikes','Low PC1 spikes','Fluorescence'});
subplot(2,1,2); hold on;
plot(lags_t,xcov(roi_trace,frame_spikes_1,corr_lags),'k','linewidth',2)
plot(lags_t,xcov(roi_trace,frame_spikes_2,corr_lags),'b','linewidth',2)
line([0,0],ylim,'linestyle','--','color','k');
xlabel('Lag (s)');
ylabel('Cross-correlation')
legend({'High PC1 spikes','Low PC1 spikes'});


% NOTE: this didn't really look like anything in particular
% Possibly interesting: if one follows the other, is it maybe cyclic?
% would jPCA give interesting components? actually, what to mark
% churchland's cross correlations look like?
% smooth_length = 30;
% smooth_filt = ones(1,smooth_length)/smooth_length;
% binned_spikes_smooth = conv2(binned_spikes,smooth_filt,'same');
% jPC_scores = AP_jPCA(zscore(binned_spikes_smooth,[],2)',size(binned_spikes,2),'rotation',[]);
% 
% 
% figure; hold on;
% xlim([min(jPC_scores(:,1)),max(jPC_scores(:,1))]);
% ylim([min(jPC_scores(:,2)),max(jPC_scores(:,2))]);
% trail_points = 30;
% plot_buffer = nan(trail_points,2);
% h = plot(plot_buffer(:,1),plot_buffer(:,2));
% for i = 1:size(binned_spikes,2)
%     plot_buffer = circshift(plot_buffer,[-1,0]);
%     plot_buffer(end,:) = jPC_scores(i,:);
%     set(h,'XData',plot_buffer(:,1),'YData',plot_buffer(:,2));
%     drawnow;
% end



%% Correlation between spikes and clustered cortical areas
% (do after clustering)

framerate = 1./nanmedian(diff(frame_t));

%use_spikes = spike_times_timeline;
%use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 800 & templateDepths < Inf)-1));
%use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 0 & templateDepths < 400)-1) & ...
%    (ismember(spike_templates,use_templates(use_template_narrow)));

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);

% % Smooth spikes
% smooth_time = round(35*1);
% smooth_filt = ones(1,smooth_time)./smooth_time;
% frame_spikes = conv2(frame_spikes,smooth_filt,'same');

corr_lags = 35;
cluster_xcorr = nan(size(cluster_trace,1),corr_lags*2+1);

%cluster_trace_zscore = zscore(cluster_trace,[],2);

% Max cross correlation
for curr_cluster = 1:size(cluster_trace,1) 
    cluster_xcorr(curr_cluster,:) = xcorr(cluster_trace(curr_cluster,:),frame_spikes,corr_lags,'coeff');
end

cluster_max_xcorr = zeros(size(kgrp));
for curr_cluster = 1:size(cluster_trace,1) 
   cluster_max_xcorr(kgrp == curr_cluster) = max(cluster_xcorr(curr_cluster,:)); 
end
figure;imagesc(cluster_max_xcorr);colormap(gray);


%% Correlation between spikes (or any trace) and pixel fluorescence

% Use the corrected impulse response for convolving kernel
use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 0 & templateDepths < 500)-1));
%use_spikes = spike_times_timeline(ismember(spike_templates,427));
%use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 0 & templateDepths < Inf)-1) & ...
%    ismember(spike_templates,use_templates(use_template_narrow))-1);

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);

if exist('gcamp_kernel','var')
    frame_spikes_conv_full = conv(frame_spikes,gcamp_kernel);
    frame_spikes_conv = frame_spikes_conv_full(1:length(frame_spikes));
end

% Set the trace to use
%use_trace = interp1(facecam_t(~isnan(facecam_t)),facecam.proc.data.groom.motion(~isnan(facecam_t)),frame_t);
%use_trace = frame_spikes_conv;
use_trace = frame_spikes;
%use_trace = smooth(frame_spikes,35);
%use_trace = wheel_speed;
%use_trace = interp1(t,beta_power,frame_t);

% Get cross correlation with all pixels
corr_lags = round(framerate*3);

framerate = 1./nanmedian(diff(frame_t));
v_xcorr = nan(size(fV,1),corr_lags*2+1);
for curr_u = 1:size(U,3)  
    v_xcorr(curr_u,:) = xcorr(fV(curr_u,:),use_trace,corr_lags);
end
lags_t = (-corr_lags:corr_lags)/framerate;
svd_xcorr = svdFrameReconstruct(U,v_xcorr);

% Draw the movie
AP_image_scroll(svd_xcorr,lags_t);

% Just correlation with V, and correlation in pixel space
v_corr = nan(size(fV,1),1);
for curr_u = 1:size(U,3)  
    curr_corr = corrcoef(fV(curr_u,:),use_trace);
    v_corr(curr_u,:) = curr_corr(2);
end
svd_corr = svdFrameReconstruct(U,v_corr);
figure;imagesc(svd_corr);


% Correlation in pixel space
U_downsample_factor = 10;
maxlags = round(35/2);

% Make mask
h = figure;
imagesc(avg_im);
set(gca,'YDir','reverse');
colormap(gray);
caxis([0 prctile(avg_im(:),90)]);
title('Draw mask to use pixels')
roiMask = roipoly;
delete(h);
roiMaskd = imresize(roiMask,1/U_downsample_factor,'bilinear') > 0;

Ud = imresize(U,1/U_downsample_factor,'bilinear');
Ud_flat = reshape(Ud,[],size(U,3));
svd_corr = nan(size(Ud,1),size(Ud,2));
for curr_px = find(roiMaskd)';
    curr_trace = Ud_flat(curr_px,:)*fV;
    %curr_corr = corrcoef(curr_trace,use_trace);
    %svd_corr(curr_px) = curr_corr(2);
    curr_corr = xcorr(curr_trace,use_trace,maxlags,'coeff');
    svd_corr(curr_px) = max(curr_corr);
    disp(curr_px/find(roiMaskd,1,'last'));
end
figure;imagesc(svd_corr);colormap(gray)


%% Get spike trace by depth, get contribution to fluorescence
% dirty: needs roi_trace, and set the thing to correlate below

unique_depths = sort(unique(templateDepths));

n_depth_groups = 8;
depth_group_edges = linspace(min(templateDepths),max(templateDepths),n_depth_groups+1);
depth_group = discretize(templateDepths,depth_group_edges);
depth_group_centers = grpstats(templateDepths,depth_group);
unique_depths = unique(depth_group);

frame_edges = [frame_t,frame_t(end)+1/framerate];

depth_mua = zeros(length(unique_depths),length(frame_edges)-1);
for curr_depth = 1:length(unique_depths);
    use_spikes = spike_times_timeline(ismember(spike_templates, ...
        find(depth_group == unique_depths(curr_depth))-1));
    depth_mua(curr_depth,:) = histcounts(use_spikes,frame_edges);
end

% Use the corrected impulse response for convolving kernel
depth_mua_conv_full = conv2(depth_mua,gcamp_kernel);
depth_mua_conv = depth_mua_conv_full(:,1:length(frame_t));

% Or, just smooth traces
smooth_kernel = ones(1,35)/35;
depth_mua_smooth = conv2(depth_mua,smooth_kernel,'same');

fluor_mua_corr = 1-pdist2(roi_trace,depth_mua,'correlation');

figure;
plot(depth_group_centers,fluor_mua_corr,'k','linewidth',2)
xlabel('Depth (\mum)');
ylabel('Correlation of conv spikes with fluorescence'); 


%% Spatiotemporal correlation-fixed spatial kernel for spikes

use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 0 & templateDepths < Inf)-1));
%use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 0 & templateDepths < 400)-1) & ...
%    ismember(spike_templates,use_templates(use_template_narrow))-1);

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
% frame_spikes_conv_full = conv(frame_spikes,gcamp_kernel);
% frame_spikes_conv = frame_spikes_conv_full(1:length(frame_spikes));

spike_trace = frame_spikes;

U_downsample_factor = 1;
Ud = imresize(U,1/U_downsample_factor,'bilinear');
Ud_flat = reshape(Ud,[],size(U,3));

% Get time-shifted spike trace
surround_frames = 35*3;
surround_t = (-surround_frames:surround_frames)/framerate;
use_spikes_shift = zeros(surround_frames*2+1,length(spike_trace));
for curr_frame_idx = 1:surround_frames*2+1
   curr_shift = curr_frame_idx - surround_frames;
   use_spikes_shift(curr_frame_idx,:) = ...
       circshift(spike_trace,[0,curr_shift]);
end

use_svs = 1:2000;
k = use_spikes_shift*fV(use_svs,:)'*(diag(1./dataSummary_n.dataSummary.Sv(use_svs)))*Ud_flat(:,use_svs)';
k2 = reshape(k',size(Ud,1),size(Ud,2),surround_frames*2+1);
k2_blur = imgaussfilt(k2,1);
AP_image_scroll(k2_blur,surround_t);

figure;imagesc(max(k2_blur,[],3));
colormap(gray);

%% Spatiotemporal correlation-fixed spatial kernel for spikes (templates)

use_templates = good_templates;

U_downsample_factor = 1;
Ud = imresize(U,1/U_downsample_factor,'bilinear');
Ud_flat = reshape(Ud,[],size(U,3));

k2_all = zeros(size(Ud,1),size(Ud,2),length(use_templates));

for curr_template_idx = 1:length(use_templates)
    
    curr_template = use_templates(curr_template_idx);
    
    use_spikes = spike_times_timeline(spike_templates == curr_template);
    
    %use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 0 & templateDepths < 400)-1) & ...
    %    ismember(spike_templates,use_templates(use_template_narrow))-1);
    
    frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
    [frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
    % frame_spikes_conv_full = conv(frame_spikes,gcamp_kernel);
    % frame_spikes_conv = frame_spikes_conv_full(1:length(frame_spikes));
    
    spike_trace = frame_spikes;
    
    % Get time-shifted spike trace
    surround_frames = 7;
    surround_t = (-surround_frames:surround_frames)/framerate;
    use_spikes_shift = zeros(surround_frames*2+1,length(spike_trace));
    for curr_frame_idx = 1:surround_frames*2+1
        curr_shift = curr_frame_idx - surround_frames;
        use_spikes_shift(curr_frame_idx,:) = ...
            circshift(spike_trace,[0,curr_shift]);
    end
    
    use_svs = 1:2000;
    k = use_spikes_shift*fV(use_svs,:)'*(diag(1./dataSummary_n.dataSummary.Sv(use_svs)))*Ud_flat(:,use_svs)';
    k2 = reshape(k',size(Ud,1),size(Ud,2),surround_frames*2+1);
    k2_all(:,:,curr_template_idx) = max(k2,[],3);
    
    disp(curr_template_idx);
  
end

% Rearrange STAs by depth
[~,sort_idx] = sort(templateDepths(good_templates+1));
k2_all = k2_all(:,:,sort_idx);

AP_image_scroll(k2_all);


%% Reduced rank regression (Kenneth) by depth

n_depths = 8;
use_depths = linspace(0,1300,n_depths+1);

canonU = zeros(size(U,1),size(U,2),n_depths);
for i = 1:n_depths;
    
    use_spikes = spike_times_timeline(ismember(spike_templates, ...
        find(templateDepths > use_depths(i) & templateDepths < use_depths(i+1))-1));
    
    frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
    [frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
    
    %frame_spikes_conv_full = conv(frame_spikes,gcamp_kernel);
    %frame_spikes_conv = frame_spikes_conv_full(1:length(frame_spikes));
    
    %Y = frame_spikes';
    Y = smooth(frame_spikes,35);
    %Y = frame_spikes_conv';
    [a,b,R2,~] = CanonCor2(Y,fV');
    
    % [a b R2 V] = CanonCor2(Y, X)
    %
    % the approximation of Y based on the first n projections is:
    % Y = X * b(:,1:n) *a'(:,1:n);
    %
    % the nth variable for the ith case gives the approximation
    % Y(i,:)' = a(:,n) * b(:,n)' * X(i,:)'
    
    canonU(:,:,i) = sum(bsxfun(@times,U,permute(b(:,1),[2,3,1])),3)/size(U,3).*a(1);
    
    disp(i);
    
end

canonU_blur = imgaussfilt(canonU,3);
AP_image_scroll(canonU_blur);colormap(colormap_blueblackred);


%% Get spikes to widefield transformation and/or deconvolve V

% I don't think any of this makes sense 

use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 1000 & templateDepths < 1200)-1));

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);

shifts = -35:35;
frame_spikes_shift = zeros(length(shifts),length(frame_spikes));
for i = 1:length(shifts)
    frame_spikes_shift(i,:) = circshift(frame_spikes,[0,shifts(i)]);
end

frame_spikes_shift_conv_full = conv2(frame_spikes_shift,gcamp_kernel);
frame_spikes_shift_conv = frame_spikes_shift_conv_full(:,1:length(frame_spikes));

sv_weights = frame_spikes_shift_conv'\fV';

m = svdFrameReconstruct(U,sv_weights');


%% Get linear combination of pixel > pixel corr to spike > pixel corr
% this is probably dumb, didnt' really work

% Use the corrected impulse response for convolving kernel
use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 0 & templateDepths < 300)-1));
%use_spikes = spike_times_timeline(ismember(spike_templates,find(templateDepths > 0 & templateDepths < 400)-1) & ...
%    ismember(spike_templates,use_templates(use_template_narrow))-1);

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);

if exist('gcamp_kernel','var')
    frame_spikes_conv_full = conv(frame_spikes,gcamp_kernel);
    frame_spikes_conv = frame_spikes_conv_full(1:length(frame_spikes));
end

% Set the trace to use
use_trace = frame_spikes_conv;

% Correlation in pixel space
U_downsample_factor = 10;
maxlags = round(35/2);

% Make mask
h = figure;
imagesc(avg_im);
set(gca,'YDir','reverse');
colormap(gray);
caxis([0 prctile(avg_im(:),90)]);
title('Draw mask to use pixels')
roiMask = roipoly;
delete(h);
roiMaskd = imresize(roiMask,1/U_downsample_factor,'bilinear') > 0;

Ud = imresize(U,1/U_downsample_factor,'bilinear');
Ud_flat = reshape(Ud,[],size(U,3));
svd_corr = nan(size(Ud,1),size(Ud,2));
for curr_px = find(roiMaskd)';
    curr_trace = Ud_flat(curr_px,:)*fV;
    %curr_corr = corrcoef(curr_trace,use_trace);
    %svd_corr(curr_px) = curr_corr(2);
    curr_corr = xcorr(curr_trace,use_trace,maxlags,'coeff');
    svd_corr(curr_px) = max(curr_corr);
    disp(curr_px/find(roiMaskd,1,'last'));
end
figure;imagesc(svd_corr);colormap(gray)

px_trace = Ud_flat(roiMaskd(:),:)*fV;
px_corr = corrcoef(px_trace');

x = px_corr\svd_corr(roiMaskd);

r = zeros(size(Ud,1),size(Ud,2));
r(roiMaskd) = x;
figure;imagesc(r);colormap(gray);













