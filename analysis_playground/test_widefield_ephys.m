%% Define experiment

animal = 'AP003';
day = '2016-06-02';
experiment = '3';

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
protocol_filename = get_cortexlab_filename(mpep_animal,day,experiment,'protocol','8digit');
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

flipped_banks = true;

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

% Load clusters, if they exist
cluster_filename = [data_path filesep 'cluster_groups.csv'];
if exist(cluster_filename,'file')
    fid = fopen(cluster_filename);
    cluster_groups = textscan(fid,'%d%s','HeaderLines',1);
    fclose(fid);
end

% Flip channel map and positions if banks are reversed
if flipped_banks
    channel_map = [channel_map(61:end);channel_map(1:60)];
    channel_positions = [channel_positions(61:end,:);channel_positions(1:60,:)];
end

% Default channel map/positions are from end: make from surface
channel_positions(:,2) = max(channel_positions(:,2)) - channel_positions(:,2);

template_abs = permute(max(abs(templates),[],2),[3,1,2]);
[~,max_channel_idx] =  max(template_abs,[],1);
template_depths = channel_positions(max_channel_idx,2);
spike_depths = template_depths(spike_templates+1);


% Load LFP
n_channels = str2num(header.n_channels);
lfp_filename = [data_path filesep 'lfp.dat'];
fid = fopen(lfp_filename);
lfp_all = fread(fid,[n_channels,inf],'int16');
fclose(fid);
% eliminate non-connected channels and sort by position (surface to deep)
lfp = lfp_all(flipud(channel_map)+1,:);
% get time of LFP sample points (NOTE: this is messy, based off of sample
% rate and knowing what kwik2dat does, not sure how accurate)
sample_rate = str2num(header.sample_rate);
lfp_cutoff = str2num(header.lfp_cutoff);
lfp_downsamp = (sample_rate/lfp_cutoff)/2;
lfp_t = ([1:size(lfp,2)]*lfp_downsamp)/sample_rate;


% Get the spike times in timeline time
%spike_times_timeline = AP_clock_fix(spike_times,sync.timestamps,photodiode.timestamps);
spike_times_timeline = AP_clock_fix(spike_times,sync.timestamps(stim_onset_idx_ephys),photodiode.timestamps(stim_onset_idx_timeline));
lfp_t_timeline = AP_clock_fix(lfp_t,sync.timestamps(stim_onset_idx_ephys),photodiode.timestamps(stim_onset_idx_timeline));

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

curr_template = 43;
surround_times = [-2,2];

framerate = 1./median(diff(frame_t));
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
AP_imscroll(sta_im_norm,sta_t);
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
% AP_imscroll(sta_im_norm,sta_t);
% caxis([-std(sta_im_norm(:))*3,std(sta_im_norm(:))*3])
% set(gcf,'Name',num2str(curr_template));

%% STA for multiunit

% use_spikes = spike_times_timeline(spike_templates == 279);

use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 500 & template_depths < 1500)));

% use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 1000 & template_depths < 2000)) &...
%     ismember(spike_templates,find(msn)));

%use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 400)) & ...
%    ismember(spike_templates,use_templates(use_template_narrow)));

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);

surround_times = [-0.5,0.5];
framerate = 1./median(diff(frame_t));
surround_frames = round(surround_times(1)*framerate):round(surround_times(2)*framerate);
sta_t = surround_frames./framerate;

sta_im = zeros(size(U,1),size(U,2),length(surround_frames));

% Prepare weighting matrix for each frame
frames_w = repmat(frame_spikes'./sum(frame_spikes),1,length(surround_frames));
for curr_sta_frame = 1:length(surround_frames)
    frames_w(:,curr_sta_frame) = ...
        circshift(frames_w(:,curr_sta_frame),[surround_frames(curr_sta_frame),0]);
end

sta_v = fV*frames_w;
sta_im = svdFrameReconstruct(U,sta_v);

% Normalize the image
%sta_im_norm = bsxfun(@rdivide,bsxfun(@minus,sta_im,px_mean),px_std);

% Draw the movie
AP_imscroll(sta_im,sta_t);
axis image;

%% STA for multiunit by depth

% Skip the first/last n seconds to do this
skip_seconds = 60;
use_frames = (frame_t > skip_seconds & frame_t < frame_t(end)-skip_seconds);
%use_frames = (frame_t > skip_seconds) & (frame_t < max(frame_t)/2);
%use_frames = (frame_t > max(frame_t)/2);

% Group multiunit by depth
n_depth_groups = 6;
%depth_group_edges = linspace(0,max(channel_positions(:,2)),n_depth_groups+1);
%depth_group_edges = linspace(700,3500,n_depth_groups+1);
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
depth_group_edges_use = depth_group_edges;
%depth_group_edges_use = [3500 Inf];

[depth_group_n,depth_group] = histc(spike_depths,depth_group_edges_use);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges_use(1:end-1)+(diff(depth_group_edges_use)/2);

framerate = 1./median(diff(frame_t));
frame_spikes = zeros(length(depth_group_edges_use)-1,length(frame_t));
for curr_depth = 1:length(depth_group_edges_use)-1
    
    %     curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    curr_spike_times = spike_times_timeline((depth_group == curr_depth) & ...
        ismember(spike_templates,find(msn)));
    
    % Discretize spikes into frames and count spikes per frame
    frame_edges = [frame_t,frame_t(end)+1/framerate];
    frame_spikes(curr_depth,:) = histcounts(curr_spike_times,frame_edges);
    
end

% Nope: array too big
% use_spikes = ismember(spike_templates,find(msn)) & (spike_times_timeline > frame_t(1)+skip_seconds) & ...
%     (spike_times_timeline < frame_t(end)-skip_seconds);
% [avgPeriEventV, periEventV] = eventLockedAvgSVD(U,fV,frame_t,spike_times_timeline(use_spikes),depth_group(use_spikes),[-3,3])

surround_times = [-0.5,1];
framerate = 1./median(diff(frame_t));
surround_frames = round(surround_times(1)*framerate):round(surround_times(2)*framerate);
sta_t = surround_frames./framerate;

sta_im = zeros(size(U,1),size(U,2),length(surround_frames));

% Prepare weighting matrix for each frame
frames_w = repmat(bsxfun(@rdivide,frame_spikes,sum(frame_spikes,2))',1,1,length(surround_frames));
for curr_sta_frame = 1:length(surround_frames);
    frames_w(:,:,curr_sta_frame) = ...
        circshift(frames_w(:,:,curr_sta_frame),[surround_frames(curr_sta_frame),0]);
end

sta_im = nan(size(U,1),size(U,2),length(surround_frames),size(frame_spikes,1));
for curr_depth = 1:size(frame_spikes,1)
   curr_surround_v = fVdf_deconv*permute(frames_w(:,curr_depth,:),[1,3,2]); 
   sta_im(:,:,:,curr_depth) = svdFrameReconstruct(Udf,curr_surround_v);
end

% Normalize the image
%sta_im_norm = bsxfun(@rdivide,bsxfun(@minus,sta_im,px_mean),px_std);

% Draw the movie
AP_imscroll(sta_im,sta_t);

% Plot map of cortical pixel by preferred depth of probe
r_px_max = squeeze(max(sta_im,[],3));
r_px_max_zeronan = r_px_max;
r_px_max_zeronan(isnan(r_px_max_zeronan)) = 0;
r_px_max_norm = reshape(zscore(reshape(r_px_max_zeronan,[],size(r_px_max,3)),[],1),size(r_px_max));
r_px_com = nansum(bsxfun(@times,r_px_max_norm,permute(1:n_depth_groups,[1,3,2])),3)./nansum(r_px_max_norm,3);

r_px_com_col = ind2rgb(round(mat2gray(r_px_com,[1,n_depth_groups])*255),jet(255));
figure;
a1 = axes('YDir','reverse');
imagesc(avg_im); colormap(gray); caxis([0,prctile(avg_im(:),99)]);
axis off; 
a2 = axes('Visible','off'); 
p = imagesc(r_px_com_col);
axis off; 
set(p,'AlphaData',mat2gray(max(r_px_max,[],3), ...
     [0,double(prctile(reshape(max(r_px_max,[],3),[],1),99.5))]));
set(gcf,'color','w');

c1 = colorbar('peer',a1,'Visible','off');
c2 = colorbar('peer',a2);
ylabel(c2,'Depth (\mum)');
colormap(c2,jet);
set(c2,'YDir','reverse');
set(c2,'YTick',linspace(0,1,6));
set(c2,'YTickLabel',linspace(depth_group_edges(1),depth_group_edges(end),6));


%% STA for templates
template_sta = zeros(size(U,1),size(U,2),length(good_templates_idx));

for curr_template_idx = 1:length(good_templates_idx)
    
    curr_template = good_templates_idx(curr_template_idx) + 1;
    
    use_spikes = spike_times_timeline(spike_templates == curr_template);
        
    frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
    [frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
    
    surround_times = [0,0.5];
    framerate = 1./median(diff(frame_t));
    surround_frames = round(surround_times(1)*framerate):round(surround_times(2)*framerate);
    sta_t = surround_frames./framerate;
        
    % Prepare weighting matrix for each frame
    frames_w = repmat(frame_spikes'./sum(frame_spikes),1,length(surround_frames));
    for curr_sta_frame = 1:length(surround_frames)
        frames_w(:,curr_sta_frame) = ...
            circshift(frames_w(:,curr_sta_frame),[surround_frames(curr_sta_frame),0]);
    end
    
    sta_v = fV*frames_w;
    sta_im = svdFrameReconstruct(U,sta_v);
    
    template_sta(:,:,curr_template_idx) = max(sta_im,[],3);
     
    disp(curr_template_idx)
end

% Rearrange STAs by depth
[~,sort_idx] = sort(template_depths(good_templates_idx));
template_sta = template_sta(:,:,sort_idx);

% Plot
AP_imscroll(template_sta,good_templates_idx(sort_idx));


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
[~,sort_idx] = sort(template_depths(use_templates+1));
template_sta = template_sta(:,:,sort_idx);

% Plot
AP_imscroll(template_sta,use_templates(sort_idx));

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
scatter3(score(:,1),score(:,2),score(:,3),50,template_depths(use_templates+1),'filled');
xlabel('PC1');ylabel('PC2');zlabel('PC3');
c = colorbar;
ylabel(c,'Depth (\mum)');
axis square

% Plot first three cluster STAs
figure;
scatter3(cluster_sta_nonan(1,:),cluster_sta_nonan(2,:), ...
    cluster_sta_nonan(3,:),50,template_depths(use_templates+1),'filled');
xlabel('PC1');ylabel('PC2');zlabel('PC3');
c = colorbar;
ylabel(c,'Depth (\mum)');
axis square

% Plot template by depth and colored by PC score
x_vals = rand(size(use_templates));
figure;scatter(x_vals,template_depths(use_templates+1),50,'.k')

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
[~,sort_idx] = sort(template_depths(use_templates+1));
template_corr = template_corr(:,:,sort_idx);

% Plot correlation by area
template_corr(isnan(template_corr)) = 0;
AP_imscroll(template_corr,use_templates(sort_idx));
caxis([-0.2 0.2]);
colormap(redblue)

figure;

% Plot vs depth
subplot(1,2,1);
scatter3(cluster_corr(1,:),cluster_corr(2,:), ...
    cluster_corr(3,:),50,template_depths(use_templates+1),'filled');
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
plot(cluster_corr(1,:),template_depths(use_templates+1),'.k','MarkerSize',10)
set(gca,'YDir','reverse');
ylabel('Depth (\mum)');
xlabel('Correlation with fluorescence')
title('ROI 1')
subplot(1,3,2)
plot(cluster_corr(2,:),template_depths(use_templates+1),'.k','MarkerSize',10)
set(gca,'YDir','reverse');
ylabel('Depth (\mum)');
xlabel('Correlation with fluorescence')
title('ROI 2')
subplot(1,3,3)
plot(cluster_corr(3,:),template_depths(use_templates+1),'.k','MarkerSize',10)
set(gca,'YDir','reverse');
ylabel('Depth (\mum)');
xlabel('Correlation with fluorescence')
title('ROI 3')


%% Group multiunit by depth, get STAs

% Group by depth
n_depth_groups = 5;
depth_group_edges = linspace(0,max(template_depths),n_depth_groups+1);
depth_group_edges(end) = Inf;

[depth_group_n,depth_group] = histc(spike_depths,depth_group_edges);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

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
AP_imscroll(sta_im,sta_t);

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
AP_imscroll(sta_im_norm,sta_t);
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
use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 1000 & template_depths < Inf)));

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
use_time = [0,620];

n_rois = 2;

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

unique_depths = sort(unique(template_depths));

nucleus_colors = lines(5);
nucleus_borders = [0,220,360,500,1150,max(template_depths)];
depth_nucleus = discretize(unique_depths,nucleus_borders);

bin_edges = use_time(1):0.001:use_time(end);
bin_centers = bin_edges(1:end-1) + (diff(bin_edges)/2);
raster_x_time = cell(length(unique_depths),1);
raster_y = cell(length(unique_depths),1);
binned_spikes_depth = zeros(length(unique_depths),length(corr_edges)-1);
for curr_depth = 1:length(unique_depths);
    curr_binned_spikes = histcounts(spike_times_timeline( ...
        ismember(spike_templates,find(template_depths == unique_depths(curr_depth)))), ...
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
template_nucleus = discretize(template_depths,nucleus_borders);

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
        ismember(spike_templates,find(template_depths > use_depths(1) & template_depths < use_depths(2)))), ...
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

[~,sort_idx] = sort(template_depths,'ascend');
binned_spikes_depthsort = binned_spikes(sort_idx,:);

figure;imagesc(corrcoef(binned_spikes_depthsort'));
colormap(hot);

% Nick suggestion: correlation plot by depth MUA

unique_depths = sort(unique(template_depths));

n_depth_groups = 30;
depth_group_edges = linspace(min(template_depths),max(template_depths),n_depth_groups+1);
depth_group = discretize(template_depths,depth_group_edges);
depth_group_centers = grpstats(template_depths,depth_group);
unique_depths = unique(depth_group);

corr_edges = 0:0.01:frame_t(end);

binned_spikes_depth = zeros(length(unique_depths),length(corr_edges)-1);
for curr_depth = 1:length(unique_depths);
    binned_spikes_depth(curr_depth,:) = histcounts(spike_times_timeline( ...
        ismember(spike_templates,find(depth_group == unique_depths(curr_depth)))), ...
        corr_edges);
end

figure;imagesc(depth_group_centers,depth_group_centers,corrcoef(binned_spikes_depth'));
colormap(hot);


%% TEMPORARY KENNETH FIGURE: STA from thalamic nuclei

unique_depths = sort(unique(template_depths));

nucleus_borders = [0,220,360,500,1150,max(unique_depths)];
template_nucleus = discretize(template_depths,nucleus_borders);
n_nuclei = length(nucleus_borders) - 1;

bin_edges = [frame_t,frame_t(end)+(1/framerate)];%frame_t(1):0.01:frame_t(end);
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
binned_spikes_nucleus = zeros(length(nucleus_borders)-1,length(bin_edges)-1);
for curr_nucleus = 1:length(nucleus_borders)-1;
    binned_spikes_nucleus(curr_nucleus,:) = histcounts(spike_times_timeline( ...
        ismember(spike_templates,find(template_nucleus == curr_nucleus))), ...
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
AP_imscroll(svd_xcorr,lags_t);

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
AP_imscroll(sta_im_norm,sta_t);
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
skip_seconds = 10;
use_frames = frame_t > skip_seconds;

% Get fluorescence in ROI
roi_trace_full = AP_svd_roi(Udf,fVdf,avg_im);
roi_trace = roi_trace_full(use_frames);

% Get population spikes per frame
framerate = 1./nanmedian(diff(frame_t));
frame_edges = [frame_t,frame_t(end)+1/framerate];

use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 1500)));

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
soft_reg_factor = 1e6;
x_autonorm = ifft((fft(roi_trace).*conj(fft(frame_spikes)))./(soft_reg_factor+fft(frame_spikes).*conj(fft(frame_spikes))));

plot_frames = 35*5;

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


%% Temporal autocorrelation-fixed STA
 
% Skip the first n seconds to do this
skip_seconds = 30;
use_frames = frame_t > skip_seconds;
use_t = frame_t(use_frames);

%use_spikes = spike_times_timeline;
use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 1500)));
%use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 400)) & ...
%    (ismember(spike_templates,use_templates(use_template_narrow))));

framerate = 1./nanmedian(diff(frame_t));

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
frame_spikes = frame_spikes(use_frames);

fluor_spikes_corr = bsxfun(@times,fft(fV(:,use_frames),[],2),conj(fft(frame_spikes)));
spikes_autocorr = fft(frame_spikes).*conj(fft(frame_spikes));

v_autonorm = ifft(bsxfun(@rdivide,fluor_spikes_corr,spikes_autocorr),[],2);

plot_frames = 35*2;
t_shift = [use_t(end-plot_frames+1:end)-use_t(end)-1/framerate,use_t(1:plot_frames)-use_t(1)];
v_autonorm_shift = [v_autonorm(:,end-plot_frames+1:end),v_autonorm(:,1:plot_frames)];

tfun = svdFrameReconstruct(U,v_autonorm_shift);
% tfun_norm = bsxfun(@rdivide,bsxfun(@minus,tfun,px_mean),px_std);

AP_imscroll(tfun);



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


%% XCORR between pixels and spikes

skip_frames = 35*10;

%use_spikes = spike_times_timeline;
use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < Inf)));
%use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 400)) & ...
%    (ismember(spike_templates,use_templates(use_template_narrow))));

framerate = 1./nanmedian(diff(frame_t));

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);

corr_lags = round(35*1);
v_xcorr = nan(size(fV,1),corr_lags*2+1);

for curr_u = 1:size(U,3)  
    v_xcorr(curr_u,:) = xcov(fV(curr_u,skip_frames:end)-mean(fV(curr_u,skip_frames:end)), ...
        frame_spikes(skip_frames:end) - mean(frame_spikes(skip_frames:end)),corr_lags,'biased');
    disp(curr_u);
end

lags_t = (-corr_lags:corr_lags)/framerate;

svd_xcorr = svdFrameReconstruct(U,v_xcorr);

% Normalize the image
svd_xcorr_norm = bsxfun(@rdivide,svd_xcorr,px_std*std(frame_spikes(skip_frames:end)));

% Set NaNs to zeros for better visualization
svd_xcorr_norm(isnan(svd_xcorr_norm)) = 0;

% Draw the movie
AP_imscroll(svd_xcorr_norm,lags_t);


%% Get one-shift-lag correlation in SVD space

lag = 0; % frames
skip_frames = 35*10;

%use_spikes = spike_times_timeline;
use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 600 & template_depths < 800)));
%use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 400)) & ...
%    (ismember(spike_templates,use_templates(use_template_narrow))));

% Discretize spikes into frames and count spikes per frame
frame_edges = [frame_t,frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);

frame_spikes = circshift(frame_spikes,[0,lag]);

v_cov = mean(bsxfun(@times,bsxfun(@minus,fV(:,skip_frames:end),mean(fV(:,skip_frames:end),2)), ...
    (frame_spikes(skip_frames:end)-mean(frame_spikes(skip_frames:end)))),2);

svd_corr = bsxfun(@rdivide,svdFrameReconstruct(U,v_cov),px_std.*std(frame_spikes(skip_frames:end)));

% Set NaNs to zeros for better visualization
svd_corr(isnan(svd_corr)) = 0;

figure;
imagesc(svd_corr);
colormap(colormap_blueblackred);
axis off;


%% Get one-shift-lag correlation in SVD space at different depths

lag = 3; % frames
skip_frames = 35*10;

% Group by depth
n_depth_groups = 8;
depth_group_edges = linspace(0,max(template_depths),n_depth_groups+1);
depth_group_edges(end) = Inf;

[depth_group_n,depth_group] = histc(spike_depths,depth_group_edges);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

framerate = 1./median(diff(frame_t));

svd_corr = zeros(size(U,1),size(U,2),length(unique(depth_group)));

for curr_depth = unique(depth_group)'
    
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    curr_spike_times(curr_spike_times < frame_t(2) | ...
        curr_spike_times > frame_t(end)) = [];
    
    % Discretize spikes into frames and count spikes per frame
    frame_edges = [frame_t,frame_t(end)+1/framerate];
    [frame_spikes,~,spike_frames] = histcounts(curr_spike_times,frame_edges);
        
    frame_spikes = circshift(frame_spikes,[0,lag]);
    
    v_cov = mean(bsxfun(@times,bsxfun(@minus,fV(:,skip_frames:end),mean(fV(:,skip_frames:end),2)), ...
        (frame_spikes(skip_frames:end)-mean(frame_spikes(skip_frames:end)))),2);
    
    svd_corr(:,:,curr_depth) = bsxfun(@rdivide,svdFrameReconstruct(U,v_cov),px_std.*std(frame_spikes(skip_frames:end)));

    disp(curr_depth);
    
end

% Set NaNs to zeros for better visualization
svd_corr(isnan(svd_corr)) = 0;

AP_imscroll(svd_corr);


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

curr_templates = intersect(find(template_depths > 0 & template_depths < Inf),use_templates(msn));

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
figure;plot(use_metric,template_depths(curr_templates+1),'.k','MarkerSize',15);
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
AP_imscroll([svd_xcorr_norm_1,svd_xcorr_norm_2],lags_t);


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
figure;plot(use_metric,template_depths(curr_templates+1),'.k','MarkerSize',15);
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
%use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 800 & template_depths < Inf)));
%use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 400)) & ...
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
use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 2500 & template_depths < 2700)));
%use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < Inf)) & ...
%    ismember(spike_templates,use_templates(use_template_narrow)));

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
AP_imscroll(svd_xcorr,lags_t);

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

n_depth_groups = 6;


% Choose ROI
h = figure;
imagesc(avg_im);
set(gca,'YDir','reverse');
colormap(gray);
caxis([0 prctile(avg_im(:),90)]);
roiMask = roipoly;
close(h);

% Get fluorescence across session in ROI
U_roi = reshape(U(repmat(roiMask,1,1,size(U,3))),[],size(U,3));
roi_trace = nanmean(U_roi*fV);

unique_depths = sort(unique(template_depths));

depth_group_edges = linspace(min(template_depths),max(template_depths),n_depth_groups+1);
depth_group = discretize(template_depths,depth_group_edges);
depth_group_centers = grpstats(template_depths,depth_group);
unique_depths = unique(depth_group);

frame_edges = [frame_t,frame_t(end)+1/framerate];

depth_mua = zeros(length(unique_depths),length(frame_edges)-1);
for curr_depth = 1:length(unique_depths);
    use_spikes = spike_times_timeline(ismember(spike_templates, ...
        find(depth_group == unique_depths(curr_depth))-1));
    depth_mua(curr_depth,:) = histcounts(use_spikes,frame_edges);
end

% Use the corrected impulse response for convolving kernel
if exist('gcamp_kernel','var')
    depth_mua_conv_full = conv2(depth_mua,gcamp_kernel);
    depth_mua_conv = depth_mua_conv_full(:,1:length(frame_t));
end;

% Or, just smooth traces
smooth_kernel = ones(1,35)/35;
depth_mua_smooth = conv2(depth_mua,smooth_kernel,'same');

fluor_mua_corr = 1-pdist2(roi_trace,depth_mua_smooth,'correlation');

figure;
plot(depth_group_centers,fluor_mua_corr,'k','linewidth',2)
xlabel('Depth (\mum)');
ylabel('Correlation of conv spikes with fluorescence'); 


%% Spatiotemporal correlation-fixed spatial kernel for spikes

use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 2500 & template_depths < 2700)));
%use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 400)) & ...
%    ismember(spike_templates,use_templates(use_template_narrow)));

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
lambda = 1e-1;%mean(1./dataSummary_n.dataSummary.Sv(use_svs))/2;
% S is already in fV so should probably take out first for this to work
k = use_spikes_shift*fV(use_svs,:)'*(diag(lambda+1./dataSummary_n.dataSummary.Sv(use_svs)))*Ud_flat(:,use_svs)';
k2 = reshape(k',size(Ud,1),size(Ud,2),surround_frames*2+1);
AP_imscroll(k2,surround_t);


%% Spatiotemporal correlation-fixed spatial kernel for spikes (templates)

use_templates = good_templates;

U_downsample_factor = 1;
Ud = imresize(U,1/U_downsample_factor,'bilinear');
Ud_flat = reshape(Ud,[],size(U,3));

k2_all = zeros(size(Ud,1),size(Ud,2),length(use_templates));

for curr_template_idx = 1:length(use_templates)
    
    curr_template = use_templates(curr_template_idx);
    
    use_spikes = spike_times_timeline(spike_templates == curr_template);
    
    %use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 400)) & ...
    %    ismember(spike_templates,use_templates(use_template_narrow)));
    
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
[~,sort_idx] = sort(template_depths(good_templates));
k2_all = k2_all(:,:,sort_idx);

AP_imscroll(k2_all);


%% Reduced rank regression (Kenneth) by depth

n_depths = 1;
use_depths = linspace(0,1300,n_depths+1);

canonU = zeros(size(U,1),size(U,2),n_depths);
for i = 1:n_depths;
    
    use_spikes = spike_times_timeline(ismember(spike_templates, ...
        find(template_depths > use_depths(i) & template_depths < use_depths(i+1))));
    
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
AP_imscroll(canonU_blur);colormap(colormap_blueblackred);


%% Get spikes to widefield transformation and/or deconvolve V

% I don't think any of this makes sense 

use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 1000 & template_depths < 1200)));

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
use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 300)));
%use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 400)) & ...
%    ismember(spike_templates,use_templates(use_template_narrow)));

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




%% Spatiotemporal autocorr correct fluorescence-to-spike kernel

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = frame_t > skip_seconds;
use_t = frame_t(use_frames);

%use_spikes = spike_times_timeline;
use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 2500 & template_depths < 2700)));
%use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 400)) & ...
%    (ismember(spike_templates,use_templates(use_template_narrow))));

framerate = 1./nanmedian(diff(frame_t));

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
frame_spikes = frame_spikes(use_frames);

spikes_fluor_corr = bsxfun(@times,fft(frame_spikes),conj(fft(fV(:,use_frames),[],2)));
fluor_autocorr = fft(fV(:,use_frames),[],2).*conj(fft(fV(:,use_frames),[],2));
fluor_autocorr_buffer = bsxfun(@plus,fluor_autocorr,std(fluor_autocorr(:)));

v_autonorm = ifft(bsxfun(@rdivide,spikes_fluor_corr,fluor_autocorr_buffer),[],2);

plot_frames = 35*2;
t_shift = [use_t(end-plot_frames+1:end)-use_t(end)-1/framerate,use_t(1:plot_frames)-use_t(1)];
v_autonorm_shift = [v_autonorm(:,end-plot_frames+1:end),v_autonorm(:,1:plot_frames)];

tfun = svdFrameReconstruct(U,v_autonorm_shift);
%tfun_norm = bsxfun(@rdivide,bsxfun(@minus,tfun,px_mean),px_std);

AP_imscroll(tfun);



% TEMP - DO THIS FOR ROI




% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = frame_t > skip_seconds;
use_t = frame_t(use_frames);

%use_spikes = spike_times_timeline;
use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < Inf)));
%use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 400)) & ...
%    (ismember(spike_templates,use_templates(use_template_narrow))));

framerate = 1./nanmedian(diff(frame_t));

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
frame_spikes = frame_spikes(use_frames);

spikes_fluor_corr = bsxfun(@times,fft(frame_spikes),conj(fft(roi_trace(use_frames),[],2)));
fluor_autocorr = fft(roi_trace(use_frames),[],2).*conj(fft(roi_trace(use_frames),[],2));

v_autonorm = ifft(bsxfun(@rdivide,spikes_fluor_corr,fluor_autocorr+std(fluor_autocorr)),[],2);

plot_frames = 35*2;
t_shift = [use_t(end-plot_frames+1:end)-use_t(end)-1/framerate,use_t(1:plot_frames)-use_t(1)];
v_autonorm_shift = [v_autonorm(:,end-plot_frames+1:end),v_autonorm(:,1:plot_frames)];



%% Make fake spike train based on given fluor kernel

% Choose ROI
h = figure;
imagesc(avg_im);
set(gca,'YDir','reverse');
colormap(gray);
caxis([0 prctile(avg_im(:),90)]);
roiMask = roipoly;
close(h);

% Get fluorescence across session in ROI
U_roi = reshape(U(repmat(roiMask,1,1,size(U,3))),[],size(U,3));
roi_trace = nanmean(U_roi*fV);
roi_trace_slope = [0;diff(smooth(roi_trace,2))];
roi_trace_slope(roi_trace_slope < 0) = 0;
fake_frame_spikes = poissrnd(mat2gray(roi_trace_slope)*10)';

jitter_time = median(diff(frame_t));
fake_spike_times = cell2mat(arrayfun(@(x) ...
    repmat(frame_t(x),1,fake_frame_spikes(x)), ...
    1:length(fake_frame_spikes),'uni',false));
fake_spike_times = fake_spike_times + (rand(1,length(fake_spike_times))-0.5)*jitter_time;



%% Fake spike regression: test lambdas

frame_spikes = fake_frame_spikes;

skip_seconds = 4;
use_frames = (frame_t > skip_seconds);

% Downsample and get pixels to compare across
U_downsample_factor = 10;
Ud = imresize(U,1/U_downsample_factor,'bilinear');
px = reshape(svdFrameReconstruct(Ud,fV),[],size(fV,2));

% Coherence
[px_coherence,f] = mscohere(px(:,use_frames)',frame_spikes(use_frames)',hanning(round(framerate*10)),round(framerate*5),[],framerate);
coherence_map = reshape(sum(px_coherence,1),size(Ud,1),size(Ud,2));

% Correlation
maxlags = round(35/2);
correlation_map = nan(size(Ud,1),size(Ud,2));
for curr_px = 1:size(px,1)
    curr_corr = xcorr(px(curr_px,:),frame_spikes,maxlags,'coeff');
    correlation_map(curr_px) = max(curr_corr);
end

% Regression
use_lambdas = [0,1e4,1e5,1e6,1e7,1e8];
regression_maps = nan(size(Ud,1),size(Ud,2),length(use_lambdas));
explained_var_all = nan(length(use_lambdas),1);
for curr_lambda = 1:length(use_lambdas)
    
    use_svs = 1:50;
    kernel_frames = -10:10;
    downsample_factor = 1;
    lambda = use_lambdas(curr_lambda);
    zs = false;
    cvfold = 3;
    
    kernel_frames_downsample = round(downsample(kernel_frames,downsample_factor)/downsample_factor);
    
    [k,predicted_spikes,explained_var] = ...
        AP_regresskernel(downsample(fV(use_svs,use_frames)',downsample_factor)', ...
        downsample(frame_spikes(:,use_frames)',downsample_factor)',kernel_frames_downsample,lambda,zs,cvfold);
    
    r = reshape(k,length(use_svs),length(kernel_frames_downsample),size(frame_spikes,1));
    regression_maps(:,:,curr_lambda) = max(svdFrameReconstruct(Ud(:,:,use_svs),r),[],3);
    explained_var_all(curr_lambda) = explained_var.total;
end

% Plot
sq_plot = 3;%ceil(sqrt(length(use_lambdas) + 3));

figure; colormap(gray);

subplot(sq_plot+1,sq_plot,1);
imagesc(avg_im); hold on;
mask_im = imagesc(padarray(roiMask,[0,0,2],'post'));
set(mask_im,'AlphaData',0.2);
axis off;
title('Real');

subplot(sq_plot+1,sq_plot,2);
imagesc(coherence_map);
axis off;
title('Coherence')
 
subplot(sq_plot+1,sq_plot,3);
imagesc(correlation_map);
axis off;
title('Correlation')

for curr_regression_map = 1:size(regression_maps,3)
    subplot(sq_plot+1,sq_plot,3+curr_regression_map);
    imagesc(regression_maps(:,:,curr_regression_map));
    axis off;
    title(['Regression ' num2str(use_lambdas(curr_regression_map))]);
end

subplot(sq_plot+1,1,sq_plot+1);
plot(explained_var_all,'k','linewidth',2);
ylabel('Explained variance');
set(gca,'XTick',1:length(use_lambdas));
set(gca,'XTickLabel',cellfun(@num2str,num2cell(use_lambdas),'uni',false));
xlabel('\lambda');

%% Fake spike regression: test number of SVs

frame_spikes = fake_frame_spikes;

skip_seconds = 4;
use_frames = (frame_t > skip_seconds);

% Downsample and get pixels to compare across
U_downsample_factor = 10;
Ud = imresize(U,1/U_downsample_factor,'bilinear');
px = reshape(svdFrameReconstruct(Ud,fV),[],size(fV,2));

% Coherence
[px_coherence,f] = mscohere(px(:,use_frames)',frame_spikes(use_frames)',hanning(round(framerate*10)),round(framerate*5),[],framerate);
coherence_map = reshape(sum(px_coherence,1),size(Ud,1),size(Ud,2));

% Correlation
maxlags = round(35/2);
correlation_map = nan(size(Ud,1),size(Ud,2));
for curr_px = 1:size(px,1)
    curr_corr = xcorr(px(curr_px,:),frame_spikes,maxlags,'coeff');
    correlation_map(curr_px) = max(curr_corr);
end

% Regression
use_svs_all = [25,50,75,100,200,300];
regression_maps = nan(size(Ud,1),size(Ud,2),length(use_svs_all));
explained_var_all = nan(length(use_svs_all),1);
for curr_sv = 1:length(use_svs_all)
    
    use_svs = 1:use_svs_all(curr_sv);
    kernel_frames = -10:10;
    downsample_factor = 1;
    lambda = 0;
    zs = false;
    cvfold = 3;
    
    kernel_frames_downsample = round(downsample(kernel_frames,downsample_factor)/downsample_factor);
    
    [k,predicted_spikes,explained_var] = ...
        AP_regresskernel(downsample(fV(use_svs,use_frames)',downsample_factor)', ...
        downsample(frame_spikes(:,use_frames)',downsample_factor)',kernel_frames_downsample,lambda,zs,cvfold);
    
    r = reshape(k,length(use_svs),length(kernel_frames_downsample),size(frame_spikes,1));
    regression_maps(:,:,curr_sv) = max(svdFrameReconstruct(Ud(:,:,use_svs),r),[],3);
    explained_var_all(curr_sv) = explained_var.total;
end

% Plot
sq_plot = 3;%ceil(sqrt(length(use_lambdas) + 3));

figure; colormap(gray);

subplot(sq_plot+1,sq_plot,1);
imagesc(avg_im); hold on;
mask_im = imagesc(padarray(roiMask,[0,0,2],'post'));
set(mask_im,'AlphaData',0.2);
axis off;
title('Real');

subplot(sq_plot+1,sq_plot,2);
imagesc(coherence_map);
axis off;
title('Coherence')
 
subplot(sq_plot+1,sq_plot,3);
imagesc(correlation_map);
axis off;
title('Correlation')

for curr_regression_map = 1:size(regression_maps,3)
    subplot(sq_plot+1,sq_plot,3+curr_regression_map);
    imagesc(regression_maps(:,:,curr_regression_map));
    axis off;
    title(['Vs ' num2str(use_svs_all(curr_regression_map))]);
end

subplot(sq_plot+1,1,sq_plot+1);
plot(explained_var_all,'k','linewidth',2);
ylabel('Explained variance');
set(gca,'XTick',1:length(use_svs_all));
set(gca,'XTickLabel',cellfun(@num2str,num2cell(use_svs_all),'uni',false));
xlabel('Vs');


%% Matrix division method of fluorescence -> spike kernel (SVD space)
% NOTE: the lambda scaling matrix is ridge regression
% lambda is a scalar to tell what level of V values shouldn't contribute
% (i.e. anything lower than sqrt(lambda))
% this corresponds, but not exactly, to truncating the Vs

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = frame_t > skip_seconds;


% Group multiunit by depth
n_depth_groups = 4;
depth_group_edges = linspace(800,max(template_depths),n_depth_groups+1);
depth_group_edges(end) = Inf;

[depth_group_n,depth_group] = histc(spike_depths,depth_group_edges);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

framerate = 1./median(diff(frame_t));
frame_spikes = zeros(n_depth_groups,length(frame_t));
for curr_depth = 1:n_depth_groups'
    
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    curr_spike_times(curr_spike_times + surround_times(1) < frame_t(2) | ...
        curr_spike_times + surround_times(2) > frame_t(end)) = [];
    
    % Discretize spikes into frames and count spikes per frame
    frame_edges = [frame_t,frame_t(end)+1/framerate];
    frame_spikes(curr_depth,:) = histcounts(curr_spike_times,frame_edges);
    
end

use_svs = 1:500;
kernel_frames = -10:5;

n_svs = length(use_svs);
n_kernel_frames = length(kernel_frames);

fluor_design = repmat(fV(use_svs,use_frames)',[1,1,n_kernel_frames]);

% Temporally shift each page
for curr_kernel_frame = 1:n_kernel_frames;
    fluor_design(:,:,curr_kernel_frame) = ...
        circshift(fluor_design(:,:,curr_kernel_frame),[kernel_frames(curr_kernel_frame),0,0]);
end

fluor_design = reshape(fluor_design,[],size(fluor_design,2)*size(fluor_design,3));

% Ridge regression for reducing noise: add offsets to design matrix to penalize k
lambda = std(fV(500,:))^2;
ridge_matrix = lambda*eye(size(fluor_design,2));

fluor_gpu = gpuArray([bsxfun(@minus,fluor_design,mean(fluor_design,1));ridge_matrix]);
spikes_gpu = gpuArray([bsxfun(@minus,frame_spikes(:,use_frames), ...
    mean(frame_spikes(:,use_frames),2))';zeros(size(fluor_design,2),size(frame_spikes,1))]);

% Use the pseudoinverse (pinv doesn't work on gpu) - looks the same though
k = gather(inv(fluor_gpu'*fluor_gpu)*fluor_gpu'*spikes_gpu);

% Reshape kernel and convert to pixel space
r = reshape(k,n_svs,n_kernel_frames,size(frame_spikes,1));

r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3);
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(U(:,:,use_svs),r(:,:,curr_spikes));
end

AP_imscroll(r_px,kernel_frames/framerate);

clear fluor_design

%% Regression from fluor to spikes (AP_regresskernel) MUA

% Skip the first n seconds to do this
skip_seconds = 60;
use_frames = (frame_t > skip_seconds & frame_t < (frame_t(end) - skip_seconds));
%use_frames = (frame_t > skip_seconds) & (frame_t < max(frame_t)/2);
%use_frames = (frame_t > max(frame_t)/2);

use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 600)));
% use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 1500)) &...
%     ismember(spike_templates,find(msn)));
% use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 1300 & template_depths < 2500)) &...
%     ismember(spike_templates,find(msn)) & ismember(spike_templates,find(l_r_diff < 0.5)));

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
frame_spikes = single(frame_spikes);
 
use_svs = 1:50;
kernel_frames = -17:17;
downsample_factor = 1;
lambda = 20;
zs = [false,true];
cvfold = 5;

kernel_frames_downsample = round(downsample(kernel_frames,downsample_factor)/downsample_factor);

% TO USE fV
% [k,predicted_spikes,explained_var] = ...
%     AP_regresskernel(downsample(fVdf(use_svs,use_frames)',downsample_factor)', ...
%     downsample(frame_spikes(:,use_frames)',downsample_factor)',kernel_frames_downsample,lambda,zs,cvfold);
% TO USE dfV
use_frames_idx = find(use_frames);
[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(downsample(diff(fVdf(use_svs,use_frames_idx),[],2)',downsample_factor)', ...
    downsample(frame_spikes(:,use_frames_idx(2:end))',downsample_factor)',kernel_frames_downsample,lambda,zs,cvfold);

% Reshape kernel and convert to pixel space
r = reshape(k,length(use_svs),length(kernel_frames_downsample),size(frame_spikes,1));

r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3);
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(Udf(:,:,use_svs),r(:,:,curr_spikes));
end

AP_imscroll(r_px,(kernel_frames_downsample*downsample_factor)/framerate);
caxis([prctile(r_px(:),[1,99])]*4);
truesize

%% Regression from fluor (PIXELS, DF) to spikes (AP_regresskernel) MUA

% Skip the first n seconds to do this
skip_seconds = 60;
use_frames = (frame_t > skip_seconds & frame_t < (frame_t(end) - skip_seconds));
%use_frames = (frame_t > skip_seconds) & (frame_t < max(frame_t)/2);
%use_frames = (frame_t > max(frame_t)/2);

use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 550 & template_depths < 1000)));
% use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 2000 & template_depths < 2500)) &...
%     ismember(spike_templates,find(msn)));
% use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 1300 & template_depths < 2500)) &...
%     ismember(spike_templates,find(msn)) & ismember(spike_templates,find(l_r_diff < 0.5)));

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
frame_spikes = single(frame_spikes);
 
use_svs = 1:50;
kernel_frames = -10:5;
downsample_factor = 1;
lambda = 10;
zs = [false,false];
cvfold = 5;

downsample_factor = 20;
Udfd = imresize(Udf,1/downsample_factor,'bilinear');
f = svdFrameReconstruct(Udfd,fVdf);
f = diff(f,1,3);
f(f < 0) = 0;

[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(reshape(f(:,:,use_frames),[],sum(use_frames)), ...
    frame_spikes(:,use_frames),kernel_frames,lambda,zs,cvfold);

% Reshape kernel and convert to pixel space
r = reshape(k,size(Udfd,1),size(Udfd,2),[]);

r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3);
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(Udf(:,:,use_svs),r(:,:,curr_spikes));
end

AP_imscroll(r_px,kernel_frames/framerate);
caxis([prctile(r_px(:),[1,99])]*4);
truesize

%% Evaluating prediction with Q (or L?)
% Q = -1/2(integral)f(t)^2 dt + (sum over s) f(ts)
% L = (sum over s) log f(ts) - (integral) f df

%f_real = zscore(frame_spikes(:,use_frames),[],2);
f_real = frame_spikes(:,use_frames);
f_pred = predicted_spikes;
T = frame_t(find(use_frames,1,'last')) - frame_t(find(use_frames,1,'first'));

Q = (-1/2)*(sum(f_real.^2)/T) + sum(f_pred);
Q_norm = Q/(T*sum(f_real/T).^2);




%% Regression fluor -> MUA: get lambda via cross-validation (preset lambdas)

% Skip the first n seconds to do this
skip_seconds = 60;
use_frames = (frame_t > skip_seconds & frame_t < frame_t(end)-skip_seconds);

% Group multiunit by depth
n_depths = 6;
% depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
depth_group_edges = [0,2000];
depth_group_edges_use = depth_group_edges;

[depth_group_n,depth_group] = histc(spike_depths,depth_group_edges_use);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges_use(1:end-1)+(diff(depth_group_edges_use)/2);

% Use spikes at specific depth
use_spikes = spike_times_timeline(depth_group ~= 0);
frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
frame_spikes = single(frame_spikes);

use_svs = 1:50;
kernel_t = [-0.3,0.3];
kernel_frames = round(kernel_t(1)*framerate):round(kernel_t(2)*framerate);
zs = [false,true]; % MUA has to be z-scored if large lambda
cvfold = 5;
use_frames_idx = find(use_frames);

lambdas = [0,logspace(0,2,10)];
explained_var_lambdas = nan(size(lambdas));

for curr_lambda_idx = 1:length(lambdas)
    
    curr_lambda = lambdas(curr_lambda_idx);
    
    [~,predicted_spikes,explained_var] = ...
        AP_regresskernel(conv2(diff(fVdf(use_svs,use_frames_idx),[],2),[1,1]/2,'valid'), ...
        frame_spikes(:,use_frames_idx(2:end-1)),kernel_frames,curr_lambda,zs,cvfold);
    
    explained_var_lambdas(curr_lambda_idx) = explained_var.total;    
        
    AP_print_progress_fraction(curr_lambda_idx,length(lambdas));
    
end

figure;
semilogx(lambdas,explained_var_lambdas,'k','linewidth',2);
xlabel('\lambda');
ylabel('Explained variance');

%% Regression fluor -> MUA: get lambda via cross-validation (auto stop)

% Skip the first n seconds to do this
skip_seconds = 60;
use_frames = (frame_t > skip_seconds & frame_t < frame_t(end)-skip_seconds);

% Group multiunit by depth
n_depths = 6;
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
depth_group_edges_use = depth_group_edges;

[depth_group_n,depth_group] = histc(spike_depths,depth_group_edges_use);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges_use(1:end-1)+(diff(depth_group_edges_use)/2);

% Use spikes at specific depth
use_spikes = spike_times_timeline(depth_group ~= 0);
frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
frame_spikes = single(frame_spikes);

lambda_start = 1e3;
n_reduce_lambda = 3;

use_svs = 1:50;
kernel_t = [-0.3,0.3]; 
kernel_frames = round(kernel_t(1)*framerate):round(kernel_t(2)*framerate);
zs = [false,true]; % MUA has to be z-scored if large lambda
cvfold = 5;
use_frames_idx = find(use_frames);

update_lambda = n_reduce_lambda;
curr_lambda = lambda_start;
lambdas = 0;
explained_var_lambdas = -Inf;

while update_lambda;
    
    % TO USE dfV
    [~,~,explained_var] = ...
        AP_regresskernel(conv2(diff(fVdf(use_svs,use_frames_idx),[],2),[1,1]/2,'valid'), ...
        frame_spikes(:,use_frames_idx(2:end-1)),kernel_frames,curr_lambda,zs,cvfold);
    
    lambdas(end+1) = curr_lambda;
    explained_var_lambdas(end+1) = explained_var.total;
    
    if explained_var_lambdas(end) > explained_var_lambdas(end-1)
        curr_lambda = curr_lambda*10;
    else
        lambdas(end) = [];
        explained_var_lambdas(end) = [];
        curr_lambda = curr_lambda/2;
        update_lambda = update_lambda-1;
    end
        AP_print_progress_fraction(update_lambda,update_lambda);
end

lambda = lambdas(end);
disp(['Best lambda = ' num2str(lambda) ', Frac var = ' num2str(explained_var_lambdas(end))]);

%% Regression fluor -> MUA: get lambda via cross-validation (auto stop, NEW)
% (used to be AP_estimate_lambda, then subbed that for something else)
% (autostop is bad - too sensitive to noise)

% upsample_factor = 0.2;
% sample_rate = (1/median(diff(frame_t)))*upsample_factor;
% 
% % Skip the first n seconds to do this
% skip_seconds = 60;
% time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
% time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
% 
% % Use all spikes in striatum
% use_spikes = spike_times_timeline(ismember(spike_templates, ...
%     find(template_depths > str_depth(1) & template_depths <= str_depth(2))));
% binned_spikes = single(histcounts(use_spikes,time_bins));
% 
% use_svs = 1:50;
% kernel_t = [-0.3,0.3];
% kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
% zs = [false,true]; % MUA has to be z-scored if large lambda
% cvfold = 50;
% 
% % Resample and get derivative of V
% dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
%     diff(fVdf(use_svs,:),[],2)',time_bin_centers)';
% 
% n_update_lambda = 1;
% lambda_range = [3,8]; % ^10
% n_lambdas = 100;
% 
% figure; hold on;
% set(gca,'XScale','log');
% xlabel('\lambda');
% ylabel('Explained variance');
% drawnow;
% for curr_update_lambda = 1:n_update_lambda
%   
%     lambdas = logspace(lambda_range(1),lambda_range(2),n_lambdas);
%     explained_var_lambdas = nan(n_lambdas,1);
% 
%     curr_plot = plot(lambdas,explained_var_lambdas,'linewidth',2);
%     
%     for curr_lambda_idx = 1:length(lambdas)
%         
%         curr_lambda = lambdas(curr_lambda_idx);
%         
%         [~,predicted_spikes,explained_var] = ...
%             AP_regresskernel(dfVdf_resample, ...
%             binned_spikes,kernel_frames,curr_lambda,zs,cvfold);
%         
%         explained_var_lambdas(curr_lambda_idx) = explained_var.total;
%         
%         set(curr_plot,'YData',explained_var_lambdas);
%         drawnow;
%         
%     end        
% 
%     lambda_bin_size = diff(lambda_range)/n_lambdas;
%     explained_var_lambdas_smoothed = smooth(explained_var_lambdas,3);
%     [best_lambda_explained_var,best_lambda_idx] = max(explained_var_lambdas_smoothed);
%     best_lambda = lambdas(best_lambda_idx);
%     lambda_range = log10([best_lambda,best_lambda]) + ...
%         [-lambda_bin_size,lambda_bin_size];
%         
% end
% 
% plot(lambdas,explained_var_lambdas_smoothed,'r');
% line(xlim,repmat(best_lambda_explained_var,1,2),'color','k');
% line(repmat(best_lambda,1,2),ylim,'color','k');
% 
% lambda = best_lambda;
% 
% disp(['Best lambda = ' num2str(lambda) ', Frac var = ' num2str(explained_var_lambdas(best_lambda_idx))]);


%% Regression fluor -> templates: get lambda via cross-validation (auto stop)

% Skip the first n seconds to do this
skip_seconds = 60;
use_frames = (frame_t > skip_seconds & frame_t < frame_t(end)-skip_seconds);

% (to use aligned striatum depths)
depth_group = aligned_str_depth_group;

use_templates = unique(spike_templates(~isnan(depth_group)));
binned_spikes = zeros(length(use_templates),length(time_bins)-1);
for curr_template_idx = 1:length(use_templates)    
    curr_template = use_templates(curr_template_idx);   
    curr_spike_times = spike_times_timeline(spike_templates == curr_template);   
    binned_spikes(curr_template_idx,:) = histcounts(curr_spike_times,time_bins);    
end

lambda_start = 1e3;
n_reduce_lambda = 3;

use_svs = 1:50;
kernel_frames = 0;
zs = [false,true]; % MUA has to be z-scored if large lambda
cvfold = 5;
use_frames_idx = find(use_frames);

update_lambda = n_reduce_lambda;
curr_lambda = lambda_start;
lambdas = 0;
explained_var_lambdas = -Inf;

while update_lambda;
    
    % TO USE dfV
    [~,~,explained_var] = ...
        AP_regresskernel(conv2(diff(fVdf(use_svs,use_frames_idx),[],2),[1,1]/2,'valid'), ...
        binned_spikes(:,use_frames_idx(2:end-1)),kernel_frames,curr_lambda,zs,cvfold);
    
    lambdas(end+1) = curr_lambda;
    explained_var_lambdas(end+1) = nanmedian(explained_var.total);
    
    if explained_var_lambdas(end) > explained_var_lambdas(end-1)
        curr_lambda = curr_lambda*10;
    else
        lambdas(end) = [];
        explained_var_lambdas(end) = [];
        curr_lambda = curr_lambda/2;
        update_lambda = update_lambda-1;
    end
        
end

lambda = lambdas(end);
disp(['Best lambda = ' num2str(lambda) ', Frac var = ' num2str(explained_var_lambdas(end))]);

%% Regression from fluor to spikes (AP_regresskernel) MUA depth
% NOTE that the timing is more carefully controlled in the next cell

% Skip the first/last n seconds to do this
skip_seconds = 60;
use_frames = (frame_t > skip_seconds & frame_t < frame_t(end)-skip_seconds);
%use_frames = (frame_t > skip_seconds) & (frame_t < max(frame_t)/2);
%use_frames = (frame_t > max(frame_t)/2);

% Group multiunit by depth
n_depth_groups = 6;
% depth_group_edges = linspace(0,max(channel_positions(:,2)),n_depth_groups+1);
%depth_group_edges = linspace(700,3500,n_depth_groups+1);
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
depth_group_edges_use = depth_group_edges;
%depth_group_edges_use = [3500 Inf];

[depth_group_n,depth_group] = histc(spike_depths,depth_group_edges_use);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges_use(1:end-1)+(diff(depth_group_edges_use)/2);

framerate = 1./median(diff(frame_t));
frame_spikes = zeros(length(depth_group_edges_use)-1,length(frame_t));
for curr_depth = 1:length(depth_group_edges_use)-1
    
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
%     curr_spike_times = spike_times_timeline((depth_group == curr_depth) & ...
%         ismember(spike_templates,find(msn)));

    % Discretize spikes into frames and count spikes per frame
    frame_edges = [frame_t,frame_t(end)+1/framerate];
    frame_spikes(curr_depth,:) = histcounts(curr_spike_times,frame_edges);
    
end

use_svs = 1:50;
kernel_frames = -17:17;
downsample_factor = 1;
lambda = 20;
zs = [false,true];
cvfold = 5;

use_frames_idx = find(use_frames);

% TO USE fV
% [k,predicted_spikes,explained_var] = ...
%     AP_regresskernel(fVdf(use_svs,use_frames), ...
%     frame_spikes(:,use_frames),kernel_frames,lambda,zs,cvfold);
% TO USE dfV
% [k,predicted_spikes,explained_var] = ...
%     AP_regresskernel(diff(fVdf(use_svs,use_frames_idx),[],2), ...
%    frame_spikes(:,use_frames_idx(2:end)),kernel_frames,lambda,zs,cvfold);
% TO USE deconv
[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(diff(fVdf_deconv(use_svs,use_frames_idx),[],2), ...
   frame_spikes(:,use_frames_idx(2:end)),kernel_frames,lambda,zs,cvfold,false,true);

% Reshape kernel and convert to pixel space
r = reshape(k,length(use_svs),length(kernel_frames),size(frame_spikes,1));

r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3);
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(Udf(:,:,use_svs),r(:,:,curr_spikes));
end

AP_imscroll(r_px,kernel_frames*downsample_factor/framerate);
caxis([-prctile(r_px(:),99.9),prctile(r_px(:),99.9)])
colormap(colormap_BlueWhiteRed);
axis image;

% Get center of mass for each pixel 
% r_px_max = squeeze(sqrt(sum(r_px.^2,3)));
r_px_max = squeeze(max(r_px,[],3));
r_px_max_zeronan = r_px_max;
r_px_max_zeronan(isnan(r_px_max_zeronan)) = 0;
r_px_max_norm = bsxfun(@rdivide,bsxfun(@minus,r_px_max_zeronan,min(r_px_max_zeronan,[],3)), ...
    max(bsxfun(@minus,r_px_max_zeronan,min(r_px_max_zeronan,[],3)),[],3));
r_px_com = sum(bsxfun(@times,r_px_max_norm,permute(1:n_depth_groups,[1,3,2])),3)./sum(r_px_max_norm,3);

% Plot map of cortical pixel by preferred depth of probe
r_px_com_col = ind2rgb(round(mat2gray(r_px_com,[1,n_depth_groups])*255),jet(255));
figure;
a1 = axes('YDir','reverse');
imagesc(avg_im); colormap(gray); caxis([0,prctile(avg_im(:),99.7)]);
axis off; axis image;
a2 = axes('Visible','off'); 
p = imagesc(r_px_com_col);
axis off; axis image;
set(p,'AlphaData',mat2gray(max(r_px_max,[],3), ...
     [0,double(prctile(reshape(max(r_px_max,[],3),[],1),99))]));
set(gcf,'color','w');

c1 = colorbar('peer',a1,'Visible','off');
c2 = colorbar('peer',a2);
ylabel(c2,'Depth (\mum)');
colormap(c2,jet);
set(c2,'YDir','reverse');
set(c2,'YTick',linspace(0,1,6));
set(c2,'YTickLabel',linspace(depth_group_edges(1),depth_group_edges(end),6));

%% Regression from fluor to spikes (AP_regresskernel) MUA depth - RESAMPLE

% % (to estimate lambda value)
% plot_lambda_estimate = true;
% AP_estimate_lambda;

% Set upsample value for regression
upsample_factor = 1;
sample_rate = (1/median(diff(frame_t)))*upsample_factor;

% Skip the first/last n seconds to do this
skip_seconds = 60;
time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% % (to group multiunit by depth from top)
% n_depths = 10;
% depth_group_edges = round(linspace(min(template_depths),max(template_depths),n_depths+1));
% [depth_group_n,depth_group] = histc(spike_depths,depth_group_edges);
% depth_groups_used = unique(depth_group);
% depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

% % (to group multiunit by depth within striatum)
% % n_depths = round(diff(str_depth)/200);
% n_depths = 4;
% depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
% [depth_group_n,depth_group] = histc(spike_depths,depth_group_edges);
% depth_groups_used = unique(depth_group);
% depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

% % (to use aligned striatum depths)
% n_depths = n_aligned_depths;
% depth_group = aligned_str_depth_group;

% (for manual depth)
depth_group_edges = [str_depth(2),4000];
n_depths = length(depth_group_edges) - 1;
[depth_group_n,depth_group] = histc(spike_depths,depth_group_edges);

binned_spikes = zeros(n_depths,length(time_bins)-1);
for curr_depth = 1:n_depths
    
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
%     curr_spike_times = spike_times_timeline((depth_group == curr_depth) & ...
%         ismember(spike_templates,find(msn)));
    
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
    
end

binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
binned_spikes_std(isnan(binned_spikes_std)) = 0;

use_svs = 1:200;
kernel_t = [-0.5,0.5];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
lambda = 5;
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
r_px = AP_svdFrameReconstruct(Udf(:,:,use_svs),k);

AP_imscroll(r_px,kernel_frames/framerate);
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
 
%% Regression from fluor to spikes (AP_regresskernel) MUA depth - RESAMPLE, ZS TEMPLATES

upsample_factor = 1;
sample_rate = (1/median(diff(frame_t)))*upsample_factor;

% Skip the first/last n seconds to do this
skip_seconds = 60;

time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% % (to group multiunit by depth from top)
% n_depths = 6;
% depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
% [depth_group_n,depth_group] = histc(spike_depths,depth_group_edges);
% depth_groups_used = unique(depth_group);
% depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

% (to use aligned striatum depths)
n_depths = n_aligned_depths;
depth_group = aligned_str_depth_group;

% % (for manual depth)
% depth_group_edges = [0,1500];
% n_depths = length(depth_group_edges) - 1;
% [depth_group_n,depth_group] = histc(spike_depths,depth_group_edges);

[~,depth_group_idx] = unique(spike_templates);
template_depth_group = depth_group(depth_group_idx);

binned_spikes_template = zeros(size(templates,1),length(time_bins)-1);
for curr_template = 1:size(templates,1)    
    curr_spike_times = spike_times_timeline(spike_templates == curr_template);   
    binned_spikes_template(curr_template,:) = histcounts(curr_spike_times,time_bins);    
end

binned_spikes_template_zs = zscore(binned_spikes_template,[],2);

binned_spikes = zeros(n_depths,length(time_bins)-1);
for curr_depth = 1:n_depths   
    binned_spikes(curr_depth,:) = nanmean(binned_spikes_template_zs(template_depth_group == curr_depth,:),1);   
end

% Get rid of NaNs (if no data?)
binned_spikes(isnan(binned_spikes)) = 0;

use_svs = 1:50;
kernel_t = [-0.3,0.3];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
lambda = 4; % (COMMENT OUT TO USE LAMBDA FROM ESTIMATION ABOVE)
zs = [false,true];
cvfold = 5;

fVdf_resample = interp1(frame_t,fVdf(use_svs,:)',time_bin_centers)';

% TO USE fV
% [k,predicted_spikes,explained_var] = ...
%     AP_regresskernel(fVdf_resample, ...
%     binned_spikes,kernel_frames,lambda,zs,cvfold);
% TO USE dfV
[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(conv2(diff(fVdf_resample,[],2),[1,1]/2,'valid'), ...
    binned_spikes(:,2:end-1),kernel_frames,lambda,zs,cvfold);

% Reshape kernel and convert to pixel space
r = reshape(k,length(use_svs),length(kernel_frames),size(binned_spikes,1));

r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3);
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(Udf(:,:,use_svs),r(:,:,curr_spikes));
end

AP_imscroll(r_px,kernel_frames/framerate);
caxis([-prctile(r_px(:),99.9),prctile(r_px(:),99.9)])
colormap(colormap_BlueWhiteRed);
axis image;

% Get center of mass for each pixel
% (get max r for each pixel, filter out big ones)
r_px_max = squeeze(max(r_px,[],3)).^3;
r_px_max(isnan(r_px_max)) = 0;
for i = 1:n_depths
    r_px_max(:,:,i) = medfilt2(r_px_max(:,:,i),[10,10]);
end
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
set(c2,'YTickLabel',linspace(depth_group_edges(1),depth_group_edges(end),n_depths));


%% Regression from fluor to spikes (AP_regresskernel) templates - RESAMPLE

upsample_factor = 1;
sample_rate = (1/median(diff(frame_t)))*upsample_factor;

% Skip the first/last n seconds to do this
skip_seconds = 60;

time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

use_spikes = spike_times_timeline > time_bins(1) & spike_times_timeline < time_bins(end);
template_spike_n = accumarray(spike_templates(use_spikes),1,[size(templates,1),1]);
use_templates = find(template_spike_n > 0);
use_templates = 525;

binned_spikes = zeros(length(use_templates),length(time_bins)-1);
for curr_template_idx = 1:length(use_templates)    
    curr_template = use_templates(curr_template_idx);   
    curr_spike_times = spike_times_timeline(spike_templates == curr_template);   
    binned_spikes(curr_template_idx,:) = histcounts(curr_spike_times,time_bins);    
end
binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);

dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
    diff(fVdf,[],2)',time_bin_centers)';

fVdf_deconv = AP_deconv_wf(fVdf);
fVdf_deconv(isnan(fVdf_deconv)) = 0;
fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';

use_svs = 1:200;
kernel_t = [-0.2,0.2];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
lambda = 2; % (comment out to use lambda from above)
zs = [false,false];
cvfold = 5;

% TO USE fV
% [k,predicted_spikes,explained_var] = ...
%     AP_regresskernel(fVdf_resample(use_svs,:), ...
%     binned_spikes,kernel_frames,lambda,zs,cvfold);
% TO USE dfV
% [k,predicted_spikes,explained_var] = ...
%     AP_regresskernel(dfVdf_resample(use_svs,:), ...
%     binned_spikes_std,kernel_frames,lambda*10,zs,cvfold);
% TO USE DECONV
[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(fVdf_deconv_resample(use_svs,:), ...
    binned_spikes_std, ...
    kernel_frames,lambda,zs,cvfold);

% Reshape kernel and convert to pixel space
r = reshape(k,length(use_svs),length(kernel_frames),size(binned_spikes_std,1));

r_px = arrayfun(@(x) svdFrameReconstruct(Udf(:,:,use_svs),r(:,:,x)),1:size(r,3),'uni',false);
r_px = cat(4,r_px{:});

% (kernel at t = 0)
r_px_t0 = squeeze(r_px(:,:,kernel_frames == 0,:));

% (kernel t = 0:1 mean)
r_px_t01 = squeeze(nanmean(r_px(:,:,ismember(kernel_frames,[0,1]),:),3));

% Plot map of cortical pixel by preferred depth of probe
r_px_com = sum(bsxfun(@times,r_px_t0,permute(template_depths(use_templates),[2,3,1])),3)./sum(r_px_t0,3);

depth_range_min = min(template_depths(use_templates));
depth_range_max = max(template_depths(use_templates));

r_px_binary_frac_aligned = AP_align_widefield(animal,day,r_px_binary_frac);
r_px_com_aligned = AP_align_widefield(animal,day,r_px_com);
avg_im_aligned = AP_align_widefield(animal,day,avg_im);

r_px_com_col = ind2rgb(round(mat2gray(r_px_com_aligned,[depth_range_min,depth_range_max])*255),jet(255));
figure;
a1 = axes('YDir','reverse');
imagesc(avg_im_aligned); colormap(gray); caxis([0,prctile(avg_im(:),99.7)]);
axis off; axis image;
a2 = axes('Visible','off'); 
p = imagesc(r_px_com_col);
axis off; axis image;
set(p,'AlphaData',mat2gray(r_px_binary_frac_aligned, ...
     [0,double(prctile(reshape(r_px_binary_frac_aligned,[],1),90))]));
set(gcf,'color','w');

c1 = colorbar('peer',a1,'Visible','off');
c2 = colorbar('peer',a2);
ylabel(c2,'Depth (\mum)');
colormap(c2,jet);
set(c2,'YDir','reverse');
set(c2,'YTick',linspace(0,1,6));
set(c2,'YTickLabel',linspace(depth_range_min,depth_range_max,6));
 
AP_reference_outline('ccf_aligned','b');AP_reference_outline('retinotopy','m');

% Trying out a thing: correlation of kernel across depths
[~,sort_idx] = sort(template_depths(use_templates));
r_px_binary_reshape = reshape(r_px_binary(:,:,sort_idx),[],size(r_px_binary,3));
r_px_binary_reshape_medfilt = medfilt1(+r_px_binary_reshape',20)';
figure;imagesc(corrcoef(r_px_binary_reshape_medfilt));
axis square; colormap(gray);
xlabel('Template');
ylabel('Template');


%% Matrix division method of fluorescence -> spike kernel (pixel space)
% this doesn't make sense that it doesn't work at all...

U_downsample_factor = 15;

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

% Super downsample the images, do it in pixel space
Ud = imresize(U,1/U_downsample_factor,'bilinear');
Udpx = reshape(Ud(repmat(roiMaskd,1,1,size(Ud,3))),[],size(Ud,3));

% Reconstruct all data
px = Udpx*fV;

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = frame_t > skip_seconds;
use_t = frame_t(use_frames);

use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 800 & template_depths < Inf)));

framerate = 1./nanmedian(diff(frame_t));

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);


kernel_frames = -6:3;
lambda = 1e4;

k = AP_regresskernel(px(:,use_frames),frame_spikes(:,use_frames),kernel_frames,lambda);

% Reshape kernel
k = reshape(k,size(px,1),length(kernel_frames));
r = zeros(size(Ud,1),size(Ud,2),length(kernel_frames));
r(repmat(roiMaskd,1,1,length(kernel_frames))) = k;

AP_imscroll(r,kernel_frames/framerate);


%% Hemo corr check: spike-fluorescence coherence/kernel

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = frame_t > skip_seconds;

use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 500)));
frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);

% Hemo correct at different bands
hemo_freq = ...
    [0.1,2; ...
    0.1,3; ...
    0.2,3; ...    
    6,12; ...
    7,9; ...
    7,13];

fV_allhemos = cell(size(hemo_freq,1),1);
for i = 1:size(hemo_freq,1)
    curr_hemo = HemoCorrectLocal(Un,Vn_th,Vh_Un,framerate,hemo_freq(i,:),3);   
    fV_allhemos{i} = curr_hemo;
    disp(i);
end

% Get fluorescence from ROI, get correlation between un/corrected
h = figure;
imagesc(avg_im);
set(gca,'YDir','reverse');
colormap(gray);
caxis([0 prctile(avg_im(:),90)]);
roiMask = roipoly;
close(h);

U_roi_c = reshape(U(repmat(roiMask,1,1,size(U,3))),[],size(U,3));
roi_trace_c = cell2mat(cellfun(@(x) nanmean(U_roi_c*x)',fV_allhemos','uni',false));

U_roi_n = reshape(Un(repmat(roiMask,1,1,size(U,3))),[],size(U,3));
roi_trace_n = nanmean(U_roi_n*Vn);

[c,f] = mscohere(roi_trace_c(use_frames,:),frame_spikes(use_frames)', ...
    [],[],[],framerate);

[n,f] = mscohere(roi_trace_n(use_frames)',frame_spikes(use_frames)', ...
   [],[],[],framerate);

figure;
subplot(3,1,1:2);
hold on
plot(f,conv2(bsxfun(@minus,c,n),ones(50,1)/50,'same'));
ylabel('\DeltaCoherence (corrected - uncorrected)');
xlabel('Frequency');
legend(cellfun(@num2str,mat2cell(hemo_freq,ones(size(hemo_freq,1),1),2),'uni',false));
line(xlim,[0,0],'color','k')
title([animal ' ' day ' ' experiment])

subplot(3,1,3);
plot(sum(bsxfun(@minus,c,n),1),'k','linewidth',2);
set(gca,'XTick',1:size(hemo_freq,1));
set(gca,'XTickLabel',(cellfun(@num2str,mat2cell(hemo_freq,ones(size(hemo_freq,1),1),2),'uni',false)));
xlabel('Correction frequency range')
ylabel('Summed \DeltaCoherence');

% Impulse response
c_kernel = ifft(bsxfun(@rdivide,bsxfun(@times,fft(roi_trace_c(use_frames,:)), ...
    conj(fft(frame_spikes(use_frames)'))),fft(frame_spikes(use_frames)').*conj(fft(frame_spikes(use_frames)'))));
plot_frames = [round(framerate*-0.5),round(framerate*6)];
t_shift = [frame_t(end+plot_frames(1)+1:end)-frame_t(end)-1/framerate,frame_t(1:plot_frames(2))-frame_t(1)];
c_kernel_shift = [c_kernel(end+plot_frames(1)+1:end,:);c_kernel(1:plot_frames(2),:)];
n_kernel = ifft(bsxfun(@rdivide,bsxfun(@times,fft(roi_trace_n(use_frames)'), ...
    conj(fft(frame_spikes(use_frames)'))),fft(frame_spikes(use_frames)').*conj(fft(frame_spikes(use_frames)'))));
n_kernel_shift = [n_kernel(end+plot_frames(1)+1:end,:);n_kernel(1:plot_frames(2),:)];

figure;
subplot(2,1,1); hold on;
plot(t_shift,c_kernel_shift);
plot(t_shift,n_kernel_shift,'k');
line([t_shift(1) t_shift(end)],[0,0],'color','k')
xlabel('Time from spike (s)')
ylabel('Fluorescence');
xlim([t_shift(1) t_shift(end)])
title([animal ' ' day ' ' experiment])
legend(cellfun(@num2str,mat2cell(hemo_freq,ones(size(hemo_freq,1),1),2),'uni',false));

subplot(2,1,2); hold on;
plot(t_shift,bsxfun(@minus,c_kernel_shift,n_kernel_shift));
line([t_shift(1) t_shift(end)],[0,0],'color','k')
ylabel('\DeltaFluorescence (corrected - uncorrected)')
xlabel('Time from spike (s)');
xlim([t_shift(1) t_shift(end)])

% Frequency power plot
L = sum(use_frames);
NFFT = 2^nextpow2(L);
[Pn,F] = pwelch(single(roi_trace_n(use_frames))',[],[],NFFT,framerate);
[Pc,F] = pwelch(single(roi_trace_c(use_frames,:)),[],[],NFFT,framerate);

figure; hold on;
plot(F,log10(conv2(Pc,ones(50,1)/50,'same')));
plot(F,log10(smooth(Pn,50)),'k')
xlabel('Frequency');
ylabel('Log Power');
legend(cellfun(@num2str,mat2cell(hemo_freq,ones(size(hemo_freq,1),1),2),'uni',false));
title([animal ' ' day ' ' experiment])

% % Get coherence between spikes and all pixels
% 
% Super downsample the images, do it in pixel space
U_downsample_factor = 10;
Ud = imresize(U,1/U_downsample_factor,'bilinear');

% Reconstruct all data
px = reshape(svdFrameReconstruct(Ud,fV),[],size(fV,2));

[px_coherence,f] = mscohere(px(:,use_frames)',frame_spikes(use_frames)',hanning(round(framerate*10)),round(framerate*5),[],framerate);
px_coherence_sum = reshape(sum(px_coherence,1),size(Ud,1),size(Ud,2));

figure;imagesc(px_coherence_sum);colormap(gray);caxis([0,max(caxis)]);axis off;

px_coherence_allf = reshape(permute(px_coherence,[2,3,1]),size(Ud,1),size(Ud,2),[]);
AP_imscroll(px_coherence_allf,f);


%% Fluor-spike regression (add spike history)

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = (frame_t > skip_seconds);

use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 700)));

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
frame_spikes = single(frame_spikes);

use_svs = 1:500;
fluor_kernel_frames = -6:3;
spike_kernel_frames = -5:-1;
lambda = 10;

[k,predicted_spikes,explained_var] = ...
    AP_regresskernel({fV(use_svs,use_frames),frame_spikes(:,use_frames)}, ...
    frame_spikes(:,use_frames),{fluor_kernel_frames,spike_kernel_frames},lambda);

% Reshape kernel and convert to pixel space
px_regressors_idx = 1:length(use_svs)*length(fluor_kernel_frames);
spike_regressors_idx = px_regressors_idx(end)+1:px_regressors_idx(end)+length(spike_kernel_frames);
r = reshape(k(px_regressors_idx),length(use_svs),length(kernel_frames),size(frame_spikes,1));

r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3);
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(U(:,:,use_svs),r(:,:,curr_spikes));
end

AP_imscroll(r_px,kernel_frames/framerate);
caxis([prctile(r_px(:),[1,99])]*2)

figure;plot(spike_kernel_frames/framerate,k(spike_regressors_idx),'k')
ylabel('Weight');
xlabel('Spike lag (s)');

%% Fluor/past spike/stim (sparse noise) -> spike regression (MUA)

[Uy,Ux,nSV] = size(U);

myScreenInfo.windowPtr = NaN; % so we can call the stimulus generation and it won't try to display anything
stimNum = 1;
ss = eval([Protocol.xfile(1:end-2) '(myScreenInfo, Protocol.pars(:,stimNum));']);
stim_screen = cat(3,ss.ImageTextures{:});
ny = size(stim_screen,1);
nx = size(stim_screen,2);

warning('for now just grab stim times from retinotopy code')
% Interpolate stim screen for all times
stim_screen_interp = abs(single(interp1(stim_times,reshape(stim_screen,[],length(stim_times))',frame_t,'nearest')'));
stim_screen_interp(isnan(stim_screen_interp)) = 0;

% Use only the onsets of stim
stim_screen_interp = single([zeros(size(stim_screen_interp,1),1),diff(stim_screen_interp,[],2)] == 1);

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = (frame_t > skip_seconds);

use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 800)));
%use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 800)) & ...
%     ismember(spike_templates,find(narrow)));

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
frame_spikes = single(frame_spikes);

use_svs = 1:200;
fluor_kernel_frames = -40:10;
spike_kernel_frames = -5:-1;
stim_kernel_frames = -5:10;
lambda = 1000;
zs = true;

[k,predicted_spikes,explained_var] = ...
    AP_regresskernel({fV(use_svs,use_frames),frame_spikes(:,use_frames),stim_screen_interp(:,use_frames)}, ...
    frame_spikes(:,use_frames),{fluor_kernel_frames,spike_kernel_frames,stim_kernel_frames},lambda,zs);

% Reshape kernel and convert to pixel space
px_regressors_idx = 1:length(use_svs)*length(fluor_kernel_frames);
spike_regressors_idx = px_regressors_idx(end)+1:px_regressors_idx(end)+length(spike_kernel_frames);
stim_regressors_idx = spike_regressors_idx(end)+1:spike_regressors_idx(end)+ ...
    size(stim_screen_interp,1)*length(stim_kernel_frames);

r = reshape(k(px_regressors_idx),length(use_svs),length(fluor_kernel_frames),size(frame_spikes,1));
r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3);
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(U(:,:,use_svs),r(:,:,curr_spikes));
end

AP_imscroll(r_px,fluor_kernel_frames/framerate);
caxis([prctile(r_px(:),[1,99.5])]*2)

figure;plot(spike_kernel_frames/framerate,k(spike_regressors_idx),'k','linewidth',2)
ylabel('Weight');
xlabel('Spike history (s)');
title('Spike regressors');

stim_r = reshape(k(stim_regressors_idx),ny,nx,length(stim_kernel_frames));
AP_imscroll(stim_r,stim_kernel_frames/framerate);

%% Fluor/stim (sparse noise abs and signed) -> spike regression (templates)

[Uy,Ux,nSV] = size(U);

myScreenInfo.windowPtr = NaN; % so we can call the stimulus generation and it won't try to display anything
stimNum = 1;
ss = eval([Protocol.xfile(1:end-2) '(myScreenInfo, Protocol.pars(:,stimNum));']);
stim_screen = cat(3,ss.ImageTextures{:});
ny = size(stim_screen,1);
nx = size(stim_screen,2);

warning('for now just grab stim times from retinotopy code')
% Interpolate stim screen for all times
stim_screen_interp = single(interp1(stim_times,reshape(stim_screen,[],length(stim_times))',frame_t,'nearest')');
stim_screen_interp(isnan(stim_screen_interp)) = 0;

% Get onsets of stims
stim_screen_interp([false(size(stim_screen_interp,1),1), ...
    ~(stim_screen_interp(:,1:end-1) == 0 & stim_screen_interp(:,2:end) ~= 0)]) = 0;

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = (frame_t > skip_seconds);


frame_spikes = zeros(length(good_templates),length(frame_t),'single');
for curr_template_idx = 1:length(good_templates)
    
    curr_template = good_templates(curr_template_idx);
    
    use_spikes = spike_times_timeline(spike_templates == curr_template);
    
    %use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 0 & template_depths < 400)) & ...
    %    ismember(spike_templates,use_templates(use_template_narrow)));
    
    frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
    [curr_frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
    % frame_spikes_conv_full = conv(frame_spikes,gcamp_kernel);
    % frame_spikes_conv = frame_spikes_conv_full(1:length(frame_spikes));
    
    frame_spikes(curr_template_idx,:) = curr_frame_spikes;
    
end

use_svs = 1:200;
fluor_kernel_frames = -40:10;
stim_kernel_frames = -5:10;
lambda = 1000;

[k,predicted_spikes,explained_var] = ...
    AP_regresskernel({fV(use_svs,use_frames),stim_screen_interp(:,use_frames),abs(stim_screen_interp(:,use_frames))}, ...
    frame_spikes(:,use_frames),{fluor_kernel_frames,stim_kernel_frames,stim_kernel_frames},lambda);

% Reshape kernel and convert to pixel space
px_regressors_idx = 1:length(use_svs)*length(fluor_kernel_frames);
stim_regressors_idx = px_regressors_idx(end)+1:px_regressors_idx(end)+ ...
    size(stim_screen_interp,1)*length(stim_kernel_frames);
abs_stim_regressors_idx = stim_regressors_idx(end)+1:stim_regressors_idx(end)+ ...
    size(stim_screen_interp,1)*length(stim_kernel_frames);

r = reshape(k(px_regressors_idx,:),length(use_svs),length(fluor_kernel_frames),size(frame_spikes,1));
r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3);
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(U(:,:,use_svs),r(:,:,curr_spikes));
end

AP_imscroll(r_px,fluor_kernel_frames/framerate);
caxis([prctile(r_px(:),[1,99.5])]*2)

stim_r = reshape(k(stim_regressors_idx,:),ny,nx,length(stim_kernel_frames),size(frame_spikes,1));
AP_imscroll(stim_r,stim_kernel_frames/framerate);

abs_stim_r = reshape(k(abs_stim_regressors_idx,:),ny,nx,length(stim_kernel_frames),size(frame_spikes,1));
AP_imscroll(abs_stim_r,stim_kernel_frames/framerate);

%% Fluor/stim (stimID) -> spike regression (MUA by depth)

stim_regressors = zeros(max(unique(stimIDs)),length(frame_t),'single');
unique_stimIDs = unique(stimIDs);
for curr_stimID_idx = 1:length(unique_stimIDs)
    frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
    stim_regressors(curr_stimID_idx,:) = histcounts(stim_onsets(stimIDs == ...
        unique_stimIDs(curr_stimID_idx)),frame_edges); 
end

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = (frame_t > skip_seconds);

% Group multiunit by depth
n_depth_groups = 10;
depth_group_edges = linspace(0,double(max(channel_positions(:,2))),n_depth_groups+1);
depth_group_edges_use = depth_group_edges;
%depth_group_edges = [500,2000,Inf];
depth_group_edges_use(end) = Inf;

[depth_group_n,depth_group] = histc(spike_depths,depth_group_edges_use);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges_use(1:end-1)+(diff(depth_group_edges_use)/2);

framerate = 1./median(diff(frame_t));
frame_spikes = zeros(length(depth_group_edges_use),length(frame_t));
for curr_depth = 1:length(depth_group_edges_use)-1   
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);

    % Discretize spikes into frames and count spikes per frame
    frame_edges = [frame_t,frame_t(end)+1/framerate];
    frame_spikes(curr_depth,:) = histcounts(curr_spike_times,frame_edges);
end

use_svs = 1:200;
fluor_kernel_frames = -30:10;
stim_kernel_frames = -5:10;
lambda = 1000;

[k,predicted_spikes,explained_var] = ...
    AP_regresskernel({fV(use_svs,use_frames),stim_regressors(:,use_frames)}, ...
    frame_spikes(:,use_frames),{fluor_kernel_frames,stim_kernel_frames},lambda,true);

% Reshape kernel and convert to pixel space
px_regressors_idx = 1:length(use_svs)*length(fluor_kernel_frames);
stim_regressors_idx = px_regressors_idx(end)+1:px_regressors_idx(end)+ ...
    size(stim_regressors,1)*length(stim_kernel_frames);

r = reshape(k(px_regressors_idx,:),length(use_svs),length(fluor_kernel_frames),size(frame_spikes,1));
r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3);
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(U(:,:,use_svs),r(:,:,curr_spikes));
end

AP_imscroll(r_px,fluor_kernel_frames/framerate);
caxis([prctile(r_px(:),[1,99.5])]*2)

stim_r = reshape(k(stim_regressors_idx,:),size(stim_regressors,1),length(stim_kernel_frames),size(frame_spikes,1));
AP_imscroll(stim_r);
xlabel('Time from spike');
ylabel('Weight');
title('Stimuli');
legend(cellfun(@num2str,num2cell(unique(stimIDs))))


%% Fluor/stim (stimID) -> spike regression (MUA)

stim_regressors = zeros(max(unique(stimIDs)),length(frame_t),'single');
for curr_stimID = unique(stimIDs)'
    frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
    stim_regressors(curr_stimID,:) = histcounts(stim_onsets(stimIDs == curr_stimID),frame_edges); 
end

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = (frame_t > skip_seconds);

use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths >= 2400 & template_depths <= 2600)));

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
frame_spikes = single(frame_spikes);

use_svs = 1:50;
fluor_kernel_frames = -10:5;
stim_kernel_frames = -5:20;
lambdas = [0,0];
cv_fold = 5;

[k,predicted_spikes,explained_var] = ...
    AP_regresskernel({fV(use_svs,use_frames),stim_regressors(:,use_frames)}, ...
    frame_spikes(:,use_frames),{fluor_kernel_frames,stim_kernel_frames},lambdas,false,cv_fold);

% (stim only)
% [k,predicted_spikes,explained_var] = ...
%     AP_regresskernel(stim_regressors(:,use_frames), ...
%     frame_spikes(:,use_frames),stim_kernel_frames,0,false,cv_fold);

% Reshape kernel and convert to pixel space
px_regressors_idx = 1:length(use_svs)*length(fluor_kernel_frames);
stim_regressors_idx = px_regressors_idx(end)+1:px_regressors_idx(end)+ ...
    size(stim_regressors,1)*length(stim_kernel_frames);

r = reshape(k(px_regressors_idx),length(use_svs),length(fluor_kernel_frames),size(frame_spikes,1));
r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3);
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(U(:,:,use_svs),r(:,:,curr_spikes));
end

AP_imscroll(r_px,fluor_kernel_frames/framerate);
caxis([prctile(r_px(:),[1,99.5])]*2)

stim_r = reshape(k(stim_regressors_idx),size(stim_regressors,1),length(stim_kernel_frames));
figure;plot(stim_kernel_frames/framerate,stim_r','linewidth',2);
xlabel('Time from spike');
ylabel('Weight');
title('Stimuli');
legend(cellfun(@num2str,num2cell(unique(stimIDs))))

%% Fluor/stim (sparse noise) -> spike regression (templates separately)

[Uy,Ux,nSV] = size(U);

myScreenInfo.windowPtr = NaN; % so we can call the stimulus generation and it won't try to display anything
stimNum = 1;
ss = eval([Protocol.xfile(1:end-2) '(myScreenInfo, Protocol.pars(:,stimNum));']);
stim_screen = cat(3,ss.ImageTextures{:});
ny = size(stim_screen,1);
nx = size(stim_screen,2);

warning('for now just grab stim times from retinotopy code')
% Interpolate stim screen for all times
stim_screen_interp = abs(single(interp1(stim_times,reshape(stim_screen,[],length(stim_times))',frame_t,'nearest')'));
stim_screen_interp(isnan(stim_screen_interp)) = 0;

% Use only the onsets of stim
stim_screen_interp = single([zeros(size(stim_screen_interp,1),1),diff(stim_screen_interp,[],2)] == 1);

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = (frame_t > skip_seconds);

frame_spikes = zeros(size(templates,1),length(frame_t),'single');
for curr_template = 1:size(templates,1)
       
    use_spikes = spike_times_timeline(spike_templates == curr_template);
    
    frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
    [curr_frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
    
    frame_spikes(curr_template,:) = curr_frame_spikes;
    
end
    
use_svs = 1:200;
fluor_kernel_frames = -30:10;
stim_kernel_frames = 0:5;
lambdas = [2782559,0];
zs = false;
cvfold = 5;

[k,predicted_spikes,explained_var] = ...
    AP_regresskernel({fV(use_svs,use_frames),stim_screen_interp(:,use_frames)}, ...
    frame_spikes(:,use_frames),{fluor_kernel_frames,stim_kernel_frames},lambdas,zs,cvfold);

% Reshape kernel and convert to pixel space
px_regressors_idx = 1:length(use_svs)*length(fluor_kernel_frames);
stim_regressors_idx = px_regressors_idx(end)+1:px_regressors_idx(end)+ ...
    size(stim_screen_interp,1)*length(stim_kernel_frames)*size(frame_spikes,1);

r = reshape(k(px_regressors_idx,:),length(use_svs),length(fluor_kernel_frames),size(frame_spikes,1));
r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3);
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(U(:,:,use_svs),r(:,:,curr_spikes));
end

AP_imscroll(r_px,fluor_kernel_frames/framerate);
caxis([prctile(r_px(:),[1,99])]*2)

stim_r = reshape(k(stim_regressors_idx),ny,nx,length(stim_kernel_frames),size(frame_spikes,1));
AP_imscroll(stim_r,stim_kernel_frames);


%% Fluor/stim (stimID) -> spike regression (templates separately)

use_templates = find(template_depths >= 2400 & template_depths <= 2600);

stim_regressors = zeros(max(unique(stimIDs)),length(frame_t),'single');
for curr_stimID = unique(stimIDs)'
    frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
    stim_regressors(curr_stimID,:) = histcounts(stim_onsets(stimIDs == curr_stimID),frame_edges); 
end

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = (frame_t > skip_seconds);

frame_spikes = zeros(length(use_templates),length(frame_t),'single');
for curr_template_idx = 1:length(use_templates)
    
    curr_template = use_templates(curr_template_idx);
    
    use_spikes = spike_times_timeline(spike_templates == curr_template);
       
    frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
    [curr_frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);

    frame_spikes(curr_template_idx,:) = curr_frame_spikes;
    
end
    
use_svs = 1:200;
fluor_kernel_frames = -30:10;
stim_kernel_frames = 0:10;
lambdas = [2782559,0];
zs = false;
cvfold = 5;

[k,predicted_spikes,explained_var] = ...
    AP_regresskernel({fV(use_svs,use_frames),stim_regressors(:,use_frames)}, ...
    frame_spikes(:,use_frames),{fluor_kernel_frames,stim_kernel_frames},lambdas,zs,cvfold);

% Reshape kernel and convert to pixel space
px_regressors_idx = 1:length(use_svs)*length(fluor_kernel_frames);
stim_regressors_idx = px_regressors_idx(end)+1:px_regressors_idx(end)+ ...
    size(stim_regressors,1)*length(stim_kernel_frames);

r = reshape(k(px_regressors_idx,:),length(use_svs),length(fluor_kernel_frames),size(frame_spikes,1));
r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3);
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(U(:,:,use_svs),r(:,:,curr_spikes));
end

AP_imscroll(r_px,fluor_kernel_frames/framerate);
caxis([prctile(r_px(:),[1,99])]*2);

stim_r = reshape(k(stim_regressors_idx,:),size(stim_regressors,1),length(stim_kernel_frames),size(frame_spikes,1));
AP_imscroll(stim_r);
axis on
ylabel('Stim');
xlabel('Frame');

%% Regress spikes to fluorescence

% Skip the first n seconds to do this
skip_seconds = 30;
use_frames = (frame_t > skip_seconds & frame_t < (frame_t(end) - skip_seconds));
%use_frames = (frame_t > skip_seconds) & (frame_t < max(frame_t)/2);
%use_frames = (frame_t > max(frame_t)/2);

use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 500 & template_depths < 1500)));
% use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 2000 & template_depths < 2500)) &...
%     ismember(spike_templates,find(msn)));

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
frame_spikes = single(frame_spikes);

use_svs = 1:50;
kernel_frames = -200:200;
downsample_factor = 1;
lambda = 0;%3e5;
zs = [false,false];

kernel_frames_downsample = round(downsample(kernel_frames,downsample_factor)/downsample_factor);

[k,predicted_fluor,explained_var] = ...
    AP_regresskernel( ...
    downsample(frame_spikes(:,use_frames)',downsample_factor)', ...
    downsample(fV(use_svs,use_frames)',downsample_factor)', ...
    kernel_frames_downsample,lambda,zs);

figure;plot(kernel_frames_downsample*downsample_factor/framerate,k)

fluor_kernel = svdFrameReconstruct(U(:,:,use_svs),k');
AP_imscroll(fluor_kernel,kernel_frames_downsample*downsample_factor/framerate);
truesize;

% Get map of explained variance
downsample_factor = 10;
spatial_explained_var = AP_spatial_explained_var(U(:,:,use_svs), ...
    fV(use_svs,use_frames),predicted_fluor,downsample_factor);
figure;imagesc(spatial_explained_var);
caxis([-max(abs(spatial_explained_var(:))),max(abs(spatial_explained_var(:)))]);
colormap(colormap_BlueWhiteRed);
colorbar;


%% Regression from MUA by depth to fluorescence

% Skip the first n seconds to do this
skip_seconds = 60;
use_frames = (frame_t > skip_seconds & frame_t < (frame_t(end) - skip_seconds));

% Group multiunit by depth
n_depth_groups = 8;
depth_group_edges = linspace(min(template_depths),max(template_depths),n_depth_groups+1);
% depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
depth_group_edges_use = depth_group_edges;
%depth_group_edges_use = [400,1500,2000,2300,3000,4000];
depth_group_edges_use(end) = Inf;

[depth_group_n,depth_group] = histc(spike_depths,depth_group_edges_use);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges_use(1:end-1)+(diff(depth_group_edges_use)/2);

framerate = 1./median(diff(frame_t));
frame_spikes = zeros(length(depth_group_edges_use)-1,length(frame_t));
for curr_depth = 1:length(depth_group_edges_use)-1
    
%     curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    curr_spike_times = spike_times_timeline(depth_group == curr_depth & ...
        ismember(spike_templates,find(msn)));
    
    % Discretize spikes into frames and count spikes per frame
    frame_edges = [frame_t,frame_t(end)+1/framerate];
    frame_spikes(curr_depth,:) = histcounts(curr_spike_times,frame_edges);
    
end

use_svs = 1:50;
kernel_frames = -35*3:35*3;
downsample_factor = 1;
lambda = 0;
zs = [true,false];
cvfold = 5;

kernel_frames_downsample = round(downsample(kernel_frames,downsample_factor)/downsample_factor);

[k,predicted_fluor,explained_var] = ...
    AP_regresskernel( ...
    downsample(frame_spikes(:,use_frames)',downsample_factor)', ...
    downsample(fVdf(use_svs,use_frames)',downsample_factor)', ...
    kernel_frames_downsample,lambda,zs,cvfold);

k = permute(reshape(k,size(frame_spikes,1),length(kernel_frames_downsample),[]),[3,2,1]);

fluor_kernel = arrayfun(@(x) svdFrameReconstruct(Udf(:,:,use_svs),k(:,:,x)),1:size(k,3),'uni',false);
AP_imscroll(cat(4,fluor_kernel{:}),kernel_frames_downsample*downsample_factor/framerate);
truesize;

% Get map of explained variance
downsample_factor = 10;
spatial_explained_var = AP_spatial_explained_var(Udf(:,:,use_svs), ...
    fVdf(use_svs,use_frames),predicted_fluor,downsample_factor);
figure;imagesc(spatial_explained_var);
caxis([-max(abs(spatial_explained_var(:))),max(abs(spatial_explained_var(:)))]);
colormap(colormap_BlueWhiteRed);
colorbar;

%% Regression from MUA by depth to fluorescence - expl var map

% Skip the first n seconds to do this
skip_seconds = 60;
use_frames = (frame_t > skip_seconds & frame_t < (frame_t(end) - skip_seconds));

% Group multiunit by depth
n_depth_groups = 8;
depth_group_edges = linspace(min(template_depths),max(template_depths),n_depth_groups+1);
% depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
depth_group_edges_use = depth_group_edges;
depth_group_edges_use(end) = Inf;

[depth_group_n,depth_group] = histc(spike_depths,depth_group_edges_use);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges_use(1:end-1)+(diff(depth_group_edges_use)/2);

framerate = 1./median(diff(frame_t));
frame_spikes = zeros(length(depth_group_edges_use)-1,length(frame_t));
for curr_depth = 1:length(depth_group_edges_use)-1
    
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
%     curr_spike_times = spike_times_timeline(depth_group == curr_depth & ...
%         ismember(spike_templates,find(msn)));

    % Discretize spikes into frames and count spikes per frame
    frame_edges = [frame_t,frame_t(end)+1/framerate];
    frame_spikes(curr_depth,:) = histcounts(curr_spike_times,frame_edges);
    
end

use_svs = 1:50;
kernel_frames = -35*3:35*3;
lambda = 0;
zs = [false,false];
cvfold = 5;

spatial_explained_var = nan(size(U,1),size(U,2),n_depth_groups);
for curr_depth = 1:n_depth_groups
    
    [k,predicted_fluor,explained_var] = ...
        AP_regresskernel( ...
        frame_spikes(curr_depth,use_frames), ...
        fVdf(use_svs,use_frames), ...
        kernel_frames,lambda,zs,cvfold);
    
    % Get map of explained variance
    downsample_factor = 10;
    spatial_explained_var(:,:,curr_depth) = ...
        AP_spatial_explained_var(Udf(:,:,use_svs), ...
        fVdf(use_svs,use_frames),predicted_fluor,downsample_factor);
    
    disp(curr_depth);
    
end

% Get center of mass for each pixel 
r_px_max = spatial_explained_var;
%r_px_max = reshape(zscore(reshape(r_px_max,[],size(r_px_max,3)),[],1),size(r_px_max));
r_px_com = sum(bsxfun(@times,r_px_max,permute(1:n_depth_groups,[1,3,2])),3)./sum(r_px_max,3);

% Plot map of cortical pixel by preferred depth of probe
r_px_com_col = ind2rgb(round(mat2gray(r_px_com,[1,n_depth_groups])*255),jet(255));
figure;
a1 = axes('YDir','reverse');
imagesc(avg_im); colormap(gray); caxis([0,prctile(avg_im(:),99)]);
axis off;
a2 = axes('Visible','off');
p = imagesc(r_px_com_col);
axis off;
set(p,'AlphaData',mat2gray(max(r_px_max,[],3), ...
    [0,double(prctile(reshape(max(r_px_max,[],3),[],1),100))]));
set(gcf,'color','w');

c1 = colorbar('peer',a1,'Visible','off');
c2 = colorbar('peer',a2);
ylabel(c2,'Depth (\mum)');
colormap(c2,jet);
set(c2,'YDir','reverse');
set(c2,'YTick',linspace(0,1,6));
set(c2,'YTickLabel',linspace(depth_group_edges(1),depth_group_edges(end),6));

%% TEST: face to spikes (MUA) 

face_interp = interp1(facecam_t(~isnan(facecam_t)),facecam.proc.data.face.motionSVD(~isnan(facecam_t),:),frame_t);

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = (frame_t > skip_seconds);

use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths > 2400 & template_depths < Inf)));

frame_edges = [frame_t,frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
frame_spikes = single(frame_spikes);

use_face_svs = 1:50;
kernel_frames = -30:30;
downsample_factor = 1;
lambda = 1e3;
zs = false;
cvfold = 3;

kernel_frames_downsample = round(downsample(kernel_frames,downsample_factor)/downsample_factor);

[k,predicted_spikes,explained_var] = ...
    AP_regresskernel( ...
    downsample(face_interp(use_frames,use_face_svs),downsample_factor)', ...
    downsample(frame_spikes(:,use_frames)',downsample_factor)', ...   
    kernel_frames_downsample,lambda,zs,cvfold);

k = reshape(k,length(use_face_svs),length(kernel_frames_downsample),size(frame_spikes,1));

face_k = reshape(reshape(facecam.proc.data.face.motionMask(:,:,use_face_svs),[],length(use_face_svs))*k, ...
    facecam.proc.data.face.nY,facecam.proc.data.face.nX,[]);

AP_imscroll(face_k,kernel_frames_downsample*downsample_factor/framerate)


%% TEST: face to spikes (MUA by depth) 

face_interp = interp1(facecam_t(~isnan(facecam_t)),facecam.proc.data.face.motionSVD(~isnan(facecam_t),:),frame_t);

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = (frame_t > skip_seconds);
%use_frames = (frame_t > skip_seconds) & (frame_t < max(frame_t)/2);
%use_frames = (frame_t > max(frame_t)/2);

% Group multiunit by depth
n_depth_groups = 6;
%depth_group_edges = linspace(2400,3820,n_depth_groups+1);
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
depth_group_edges_use = depth_group_edges;
%depth_group_edges_use = [400,1500,2000,2300,3000,4000];
depth_group_edges_use(end) = Inf;

[depth_group_n,depth_group] = histc(spike_depths,depth_group_edges_use);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges_use(1:end-1)+(diff(depth_group_edges_use)/2);

framerate = 1./median(diff(frame_t));
frame_spikes = zeros(length(depth_group_edges_use)-1,length(frame_t));
for curr_depth = 1:length(depth_group_edges_use)-1
    
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
%     curr_spike_times = spike_times_timeline((depth_group == curr_depth) & ...
%         ismember(spike_templates,find(msn)));
    
    % Discretize spikes into frames and count spikes per frame
    frame_edges = [frame_t,frame_t(end)+1/framerate];
    frame_spikes(curr_depth,:) = histcounts(curr_spike_times,frame_edges);
    
end

use_face_svs = 1:50;
kernel_frames = -20:20;
downsample_factor = 1;
lambda = 1e6;
zs = [false,false];
cvfold = 5;

kernel_frames_downsample = round(downsample(kernel_frames,downsample_factor)/downsample_factor);

[k,predicted_spikes,explained_var] = ...
    AP_regresskernel( ...
    downsample(face_interp(use_frames,use_face_svs),downsample_factor)', ...
    downsample(frame_spikes(:,use_frames)',downsample_factor)', ...   
    kernel_frames_downsample,lambda,zs);

k = reshape(k,length(use_face_svs),length(kernel_frames_downsample),size(frame_spikes,1));

face_k = arrayfun(@(x) reshape(reshape(facecam.proc.data.face.motionMask(:,:,use_face_svs),[],length(use_face_svs))*k(:,:,x), ...
    facecam.proc.data.face.nY,facecam.proc.data.face.nX,[]),1:size(frame_spikes,1),'uni',false);

face_k = cat(4,face_k{:});

AP_imscroll(face_k,kernel_frames_downsample*downsample_factor/framerate)


%% Regression stim ID to spikes, then fluor to residual spikes

unique_stimIDs = unique(stimIDs);
stim_regressors = ones(length(unique_stimIDs),length(frame_t),'single');
for curr_stimID_idx = 1:length(unique_stimIDs);
    curr_stimID = unique_stimIDs(curr_stimID_idx);
    frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
    stim_regressors(curr_stimID_idx,:) = histcounts(stim_onsets(stimIDs == curr_stimID),frame_edges); 
end

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = (frame_t > skip_seconds);

use_spikes = spike_times_timeline(ismember(spike_templates,find(template_depths >= 0 & template_depths <= 1500)));

frame_edges = [frame_t(1),mean([frame_t(2:end);frame_t(1:end-1)],1),frame_t(end)+1/framerate];
[frame_spikes,~,spike_frames] = histcounts(use_spikes,frame_edges);
frame_spikes = single(frame_spikes);

% Regress stim ID alone to spikes 
stim_kernel_frames = -35*5:35*5;
lambda = 0;
cv_fold = 5;

[stim_k,predicted_spikes,explained_var] = ...
    AP_regresskernel(stim_regressors(:,use_frames), ...
    frame_spikes(:,use_frames),stim_kernel_frames,lambda,false,cv_fold);

stim_r = reshape(stim_k,size(stim_regressors,1),length(stim_kernel_frames));
figure; hold on;
set(gca,'ColorOrder',copper(size(stim_regressors,1)));
plot(stim_kernel_frames/framerate,stim_r','linewidth',2);
xlabel('Time from spike');
ylabel('Weight');
title('Stimuli');
legend(cellfun(@num2str,num2cell(unique(stimIDs)),'uni',false));

% Remove the predicted component from stimuli
residual_frame_spikes = frame_spikes(:,use_frames) - predicted_spikes;

use_svs = 1:50;
fluor_kernel_frames = -20:20;
lambdas = 1e6;
cv_fold = 5;

[fluor_k,predicted_spikes,explained_var] = ...
    AP_regresskernel(fV(use_svs,use_frames), ...
    residual_frame_spikes,fluor_kernel_frames,lambdas,false,cv_fold);

% Reshape kernel and convert to pixel space
px_regressors_idx = 1:length(use_svs)*length(fluor_kernel_frames);

r = reshape(fluor_k,length(use_svs),length(fluor_kernel_frames),size(frame_spikes,1));
r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3);
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(U(:,:,use_svs),r(:,:,curr_spikes));
end

AP_imscroll(r_px,fluor_kernel_frames/framerate);
caxis([prctile(r_px(:),[1,99.5])]*2)


%% Predict wheel movement from MUA (velocity)

% Set upsample value for regression
upsample_factor = 2;
sample_rate = (1/median(diff(frame_t)))*upsample_factor;

% Skip the first/last n seconds to do this
skip_seconds = 60;
time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Resample wheel
wheel_t = conv(Timeline.rawDAQTimestamps,[1,1]/2,'valid');
wheel_velocity_resample = interp1(wheel_t,wheel_velocity,time_bin_centers);

% (to use aligned striatum depths)
n_depths = n_aligned_depths;
depth_group = aligned_str_depth_group;

binned_spikes = zeros(n_depths,length(time_bins)-1);
for curr_depth = 1:n_depths   
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);    
end
% Get rid of NaNs (if no data?)
binned_spikes_norm = bsxfun(@rdivide,binned_spikes,nanstd(binned_spikes,[],2));   
binned_spikes_norm(isnan(binned_spikes_norm)) = 0;

kernel_samples = [-70:1:70];
lambda = 1e2;
cv = 5;

[curr_kernel,curr_predicted_signal,curr_expl_var] = AP_regresskernel( ...
        binned_spikes_norm,wheel_velocity_resample, ...
        kernel_samples,lambda,[false,false],cv);

figure; 
subplot(2,1,1); hold on;
plot(time_bin_centers,wheel_velocity_resample/max(abs(wheel_velocity_resample)),'k');
plot(time_bin_centers,curr_predicted_signal/max(abs(curr_predicted_signal)),'b');    
ylabel('Wheel velocity (max-normalized)');
xlabel('Time (s)');
legend({'Measured','Predicted'});

subplot(2,1,2); hold on; set(gca,'ColorOrder',copper(4));
curr_kernel_reshape = reshape(curr_kernel,size(binned_spikes,1),[]);
plot(kernel_samples/sample_rate,curr_kernel_reshape','linewidth',2);
axis tight
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:4),'uni',false)')
ylabel('Weight');
xlabel('Time offset (s)');


%% Predict wheel movement from MUA (speed)

% Set upsample value for regression
upsample_factor = 2;
sample_rate = (1/median(diff(frame_t)))*upsample_factor;

% Skip the first/last n seconds to do this
skip_seconds = 60;
time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Resample wheel
wheel_t = conv(Timeline.rawDAQTimestamps,[1,1]/2,'valid');
wheel_velocity_resample = abs(interp1(wheel_t,wheel_velocity,time_bin_centers));

% (to use aligned striatum depths)
n_depths = n_aligned_depths;
depth_group = aligned_str_depth_group;

binned_spikes = zeros(n_depths,length(time_bins)-1);
for curr_depth = 1:n_depths   
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);    
end
% Get rid of NaNs (if no data?)
binned_spikes(isnan(binned_spikes)) = 0;
binned_spikes_norm = bsxfun(@rdivide,binned_spikes,nanstd(binned_spikes,[],2));   

kernel_samples = [-70:1:70];
lambda = 1e2;
cv = 5;

[curr_kernel,curr_predicted_signal,curr_expl_var] = AP_regresskernel( ...
        binned_spikes_norm,wheel_velocity_resample, ...
        kernel_samples,lambda,[false,false],cv);

figure; 
subplot(2,1,1); hold on;
plot(time_bin_centers,wheel_velocity_resample/max(abs(wheel_velocity_resample)),'k');
plot(time_bin_centers,curr_predicted_signal/max(abs(curr_predicted_signal)),'b');    
ylabel('Wheel speed (max-normalized)');
xlabel('Time (s)');
legend({'Measured','Predicted'});

subplot(2,1,2); hold on; set(gca,'ColorOrder',copper(4));
curr_kernel_reshape = reshape(curr_kernel,size(binned_spikes,1),[]);
plot(kernel_samples/sample_rate,curr_kernel_reshape','linewidth',2);
axis tight
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:4),'uni',false)')
ylabel('Weight');
xlabel('Time offset (s)');

%% Predict wheel movement from MUA (+/- separately)

% Set upsample value for regression
upsample_factor = 2;
sample_rate = (1/median(diff(frame_t)))*upsample_factor;

% Skip the first/last n seconds to do this
skip_seconds = 60;
time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Resample wheel
wheel_t = conv(Timeline.rawDAQTimestamps,[1,1]/2,'valid');
wheel_velocity_resample = interp1(wheel_t,wheel_velocity,time_bin_centers);

% (to use aligned striatum depths)
n_depths = n_aligned_depths;
depth_group = aligned_str_depth_group;

binned_spikes = zeros(n_depths,length(time_bins)-1);
for curr_depth = 1:n_depths   
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);    
end
% Get rid of NaNs (if no data?)
binned_spikes(isnan(binned_spikes)) = 0;
binned_spikes_norm = bsxfun(@rdivide,binned_spikes,nanstd(binned_spikes,[],2));   

kernel_samples = [-70:1:70];
lambda = 1e2;
cv = 5;

wheel_velocity_separate = zeros(2,size(wheel_velocity_resample,2));
wheel_velocity_separate(1,wheel_velocity_resample > 0) = wheel_velocity_resample(wheel_velocity_resample > 0);
wheel_velocity_separate(2,wheel_velocity_resample < 0) = wheel_velocity_resample(wheel_velocity_resample < 0);

[curr_kernel,curr_predicted_signal,curr_expl_var] = AP_regresskernel( ...
        binned_spikes_norm,wheel_velocity_separate, ...
        kernel_samples,lambda,[false,false],cv);
    
curr_kernel_reshape = reshape(curr_kernel,size(binned_spikes,1),[],2);
    
figure;
subplot(3,1,1); hold on;
plot(time_bin_centers,wheel_velocity_resample/max(abs(wheel_velocity_resample)),'k');
plot(time_bin_centers,curr_predicted_signal'/max(abs(curr_predicted_signal(:))));    
ylabel('Wheel velocity (max-normalized)');
xlabel('Time (s)');
legend({'Measured','Predicted rightward','Predicted leftward'});

subplot(3,1,2); hold on; set(gca,'ColorOrder',copper(4));
plot(kernel_samples/sample_rate,curr_kernel_reshape(:,:,1)','linewidth',2);
axis tight
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:4),'uni',false)')
ylabel('Weight');
xlabel('Time offset (s)');

subplot(3,1,3); hold on; set(gca,'ColorOrder',copper(4));
plot(kernel_samples/sample_rate,curr_kernel_reshape(:,:,2)','linewidth',2);
axis tight
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:4),'uni',false)')
ylabel('Weight');
xlabel('Time offset (s)');


%% Predict wheel movement from fluor

% Set upsample value for regression
upsample_factor = 2;
sample_rate = (1/median(diff(frame_t)))*upsample_factor;

% Skip the first/last n seconds to do this
skip_seconds = 60;
time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Resample wheel
wheel_t = conv(Timeline.rawDAQTimestamps,[1,1]/2,'valid');
wheel_velocity_resample = interp1(wheel_t,wheel_velocity,time_bin_centers);

% resample fV 
use_svs = 1:50;
fVdf_resample = interp1(frame_t,fVdf(use_svs,:)',time_bin_centers)';
dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf(use_svs,:),[],2)',time_bin_centers)';   

kernel_samples = [-40:1:20];
lambda = 1e6;
cv = 5;

[curr_kernel,curr_predicted_signal,curr_expl_var] = AP_regresskernel( ...
        dfVdf_resample,wheel_velocity_resample, ...
        kernel_samples,lambda,[false,false],cv);

figure; 
subplot(2,1,1); hold on;
plot(time_bin_centers,wheel_velocity_resample/max(abs(wheel_velocity_resample)),'k');
plot(time_bin_centers,curr_predicted_signal/max(abs(curr_predicted_signal)),'b');    
ylabel('Wheel velocity (max-normalized)');
xlabel('Time (s)');
legend({'Measured','Predicted'});

curr_kernel_reshape = reshape(curr_kernel,size(dfVdf_resample,1),[]);

if ~exist('aUdf','var')
    aUdf = single(AP_align_widefield(animal,day,Udf));
end

curr_kernel_px = svdFrameReconstruct(aUdf(:,:,use_svs),curr_kernel_reshape);
AP_imscroll(curr_kernel_px,kernel_samples/sample_rate);
axis image
caxis([-prctile(abs(curr_kernel_px(:)),99),prctile(abs(curr_kernel_px(:)),99)]);
colormap(colormap_BlueWhiteRed);


%% Spikes raster and fluorescence in timeframe

% Set time to plot
plot_t = [1,500];

% Align U's, deconvolve widefield
use_components = 1:200;
if ~exist('aUdf','var')
    aUdf = AP_align_widefield(Udf,animal,day);
end
if ~exist('fVdf_deconv','var')
    fVdf_deconv = AP_deconv_wf(fVdf);
end

raster_fig = figure;

% (wheel velocity)
wheel_axes = subplot(6,1,6);
plot_wheel_idx = Timeline.rawDAQTimestamps >= plot_t(1) & ...
    Timeline.rawDAQTimestamps <= plot_t(2);
plot(wheel_axes,Timeline.rawDAQTimestamps(plot_wheel_idx), ...
    wheel_velocity(plot_wheel_idx),'k','linewidth',2);
ylabel('Wheel velocity');

% (striatum raster)
raster_axes = subplot(6,1,3:5,'YDir','reverse'); hold on;
plot_spikes = spike_times_timeline >= plot_t(1) & ...
    spike_times_timeline <= plot_t(2) & ...
    spike_depths >= str_depth(1) & spike_depths <= str_depth(2);
plot(raster_axes,spike_times_timeline(plot_spikes),spike_depths(plot_spikes),'.k');
ylabel('Depth (\mum)');
xlabel('Time (s)');

% (fluorescence from select ROIs)
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
roi_trace = AP_svd_roi(aUdf(:,:,use_components),fVdf_deconv(use_components,:),[],[],cat(3,wf_roi.mask));

plot_rois = [3,7,9];
fluor_spacing = []; % (use default)
fluor_axes = subplot(6,1,1:2); hold on;
plot_fluor_idx = frame_t >= plot_t(1) & frame_t <= plot_t(2);
AP_stackplot(roi_trace(plot_rois,plot_fluor_idx)', ...
    frame_t(plot_fluor_idx),fluor_spacing,[],[0,0.7,0],{wf_roi(plot_rois).area});

linkaxes([wheel_axes,raster_axes,fluor_axes],'x');







