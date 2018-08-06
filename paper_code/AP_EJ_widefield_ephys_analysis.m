%% AP analysis for EJ's 3-5Hz oscillation paper
% Simultaneous widefield and cortical electrophysiology
% Q: are 3-5Hz oscillations that EJ sees represented in spiking
% 
% All animals are tetO-GC6S mice
% Experiments are concatenated and include: spares noise, full screen
% flickering stimuli, gray screens spontaneous, and screens off spontaneous

%% Load and concantenate data from an animal/day

% % Experiment 1
% animal = 'AP026';
% day = '2017-12-09';

% % Experiment 2
% animal = 'AP027';
% day = '2017-12-09';

% Experiment 3
animal = 'AP024';
day = '2018-07-02';

% Find all experiments for that day
% (get folders with only a number - those're the experiment folders)
experiments_dir = dir(AP_cortexlab_filename(animal,day,[],'expInfo'));
experiments_num_idx = cellfun(@(x) ~isempty(x), regexp({experiments_dir.name},'^\d*$'));
experiments = cellfun(@str2num,{experiments_dir(experiments_num_idx).name});

% Loop through experiments, collate data
frame_t_all = cell(size(experiments));
fVdf_all = cell(size(experiments));
spike_times_timeline_all = cell(size(experiments));

for curr_exp = 1:length(experiments)
    experiment = experiments(curr_exp);
    AP_load_experiment;
    
    frame_t_all{curr_exp} = frame_t;
    fVdf_all{curr_exp} = fVdf;
    spike_times_timeline_all{curr_exp} = spike_times_timeline;
    
    AP_print_progress_fraction(curr_exp,length(experiments));
end

frame_t_stitch = cellfun(@(frame_t,frame_t_add) frame_t + frame_t_add, ...
    frame_t_all, num2cell(cumsum([0,cellfun(@(x) x(end)+1/framerate,frame_t_all(1:end-1))])),'uni',false);

% Concatenate all data
frame_t = cat(2,frame_t_all{:});
fVdf = cat(2,fVdf_all{:});
spike_times_timeline = cat(2,spike_times_timeline_all{:});


%% Choose depths to pool spikes for MUA

% Get MUA correlation by depth
n_corr_groups = 40;
depth_group_edges = linspace(0,max(channel_positions(:,2)),n_corr_groups+1);
depth_group = discretize(templateDepths,depth_group_edges);
depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);
unique_depths = 1:length(depth_group_edges)-1;

spike_binning = 0.01; % seconds
corr_edges = spike_times_timeline(1):spike_binning:spike_times_timeline(end);
corr_centers = corr_edges(1:end-1) + diff(corr_edges);

binned_spikes_depth = zeros(length(unique_depths),length(corr_edges)-1);
for curr_depth = 1:length(unique_depths)
    binned_spikes_depth(curr_depth,:) = histcounts(spike_times_timeline( ...
        ismember(spike_templates,find(depth_group == unique_depths(curr_depth)))), ...
        corr_edges);
end

mua_corr = corrcoef(binned_spikes_depth');

% Plot templates by depth and MUA correlation by depth, choose borders

figure;

template_depth_ax = subplot(1,2,1); hold on;
plotSpread(template_depth_ax,templateDepths,'distributionColors','k');
ylim([0,4000]);
set(gca,'YDir','reverse'); 
ylabel('Depth');
title('Templates by depth');
axis square;

mua_corr_ax = subplot(1,2,2); hold on;
imagesc(depth_group_centers,depth_group_centers,mua_corr);
caxis([0,prctile(mua_corr(:),95)]);
colormap(hot)
axis image;
ylim([0,4000]);
xlim([0,4000]);
set(gca,'YDir','reverse'); 
xlabel('Depth');
ylabel('Depth');
title('MUA correlation by depth')

[~,mua_top] = ginput(1);
line(template_depth_ax,xlim(template_depth_ax),repmat(mua_top,1,2),'color','b','linewidth',2);
line(mua_corr_ax,xlim(mua_corr_ax),repmat(mua_top,1,2),'color','b','linewidth',2);
line(mua_corr_ax,repmat(mua_top,1,2),ylim(mua_corr_ax),'color','b','linewidth',2);
drawnow;

[~,mua_bottom] = ginput(1);
line(template_depth_ax,xlim(template_depth_ax),repmat(mua_bottom,1,2),'color','b','linewidth',2);
line(mua_corr_ax,xlim(mua_corr_ax),repmat(mua_bottom,1,2),'color','b','linewidth',2);
line(mua_corr_ax,repmat(mua_bottom,1,2),ylim(mua_corr_ax),'color','b','linewidth',2);
drawnow;

mua_borders = [mua_top,mua_bottom];


%% Bin spikes to get MUA

% Discretize spikes into frames and count spikes per frame
use_spikes = spikeDepths >= mua_borders(1) & spikeDepths <= mua_borders(2);

frame_edges = [frame_t,frame_t(end)+1/framerate];
frame_spikes = histcounts(spike_times_timeline(use_spikes),frame_edges);

%% Get fluorescence -> MUA regression map, draw ROI for fluorescence

use_svs = 1:50;
kernel_frames = -30:30;
downsample_factor = 1;
lambda = 3.5;
zs = [false,true];
cvfold = 5;

% Regress MUA from fluorescence
[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(fVdf(use_svs,:), ...
    frame_spikes,kernel_frames,lambda,zs,cvfold);

% Reshape kernel and convert to pixel space
r = reshape(k,length(use_svs),length(kernel_frames),size(frame_spikes,1));

r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3)
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(Udf(:,:,use_svs),r(:,:,curr_spikes));
end

% AP_image_scroll(r_px,kernel_frames*downsample_factor/framerate);
% caxis([-prctile(r_px(:),99.9),prctile(r_px(:),99.9)])
% colormap(colormap_BlueWhiteRed);
% axis image;

max_weights = max(abs(r_px),[],3);

% Draw ROI and get fluorescence using average image and weights
[fluor_trace,fluor_mask] = AP_svd_roi(Udf,fVdf,avg_im,max_weights);

% Show the average image and regression weights with ROI drawn
first_nonzero = find(fluor_mask > 0,1);
[y_nonzero x_nonzero] = ind2sub([size(U,1),size(U,2)],first_nonzero);
roi_perim = bwtraceboundary(fluor_mask,[y_nonzero x_nonzero],'N');

figure;colormap(gray);

subplot(1,2,1); hold on;
imagesc(avg_im); 
axis image off;
caxis([0,prctile(avg_im(:),99)]);
plot(roi_perim(:,2),roi_perim(:,1),'b','linewidth',2);
set(gca,'YDir','reverse');
title('Average image')

subplot(1,2,2); hold on;
imagesc(max_weights); 
axis image off;
caxis([0,prctile(max_weights(:),99)]);
plot(roi_perim(:,2),roi_perim(:,1),'b','linewidth',2);
set(gca,'YDir','reverse');
title('MUA regression weights')


%% Plot MUA and fluorescence together

figure; hold on
plot(frame_t,zscore(frame_spikes),'k','linewidth',2);
plot(frame_t,zscore(fluor_trace),'color',[0,0.5,0],'linewidth',2);
plot(conv(frame_t,[1,1]/2,'valid'),zscore(diff(fluor_trace)),'color',[0,0,0.8],'linewidth',2);

xlabel('Time');
ylabel('Z-scored activity');
legend({'MUA','Fluorescence','Fluorescence derivative'});

%% MUA/fluorescence cross-correlation and coherence

% Cross-correlation
corr_lag = 5; % in seconds
corr_lag_samples = round(corr_lag*framerate);
[mua_fluor_xcorr,lags] = xcorr(fluor_trace,frame_spikes,corr_lag_samples,'coeff');
lags_t = lags./framerate;

% Coherence
[mua_fluor_coherence,coherence_f] = mscohere(frame_spikes',fluor_trace', ...
    hamming(round(framerate*5)),round(framerate*2.5),[],framerate);

% Plot
figure;

subplot(2,1,1);
plot(lags_t,mua_fluor_xcorr,'k','linewidth',2)
line([0,0],ylim,'linestyle','--','color','r');
xlabel('MUA Lag (s)');
ylabel('Normalized cross-correlation')

subplot(2,1,2);
plot(coherence_f,mua_fluor_coherence,'k','linewidth',2)
xlabel('Freqency');
ylabel('Coherence');


%% MUA/fluorescence spectral correlation

% Spectrogram settings
spect_overlap = 80;
window_length = 3; % in seconds
window_length_samples = round(window_length/(1/framerate));
N = window_length_samples; % window length
df = framerate/N; % frequency increment

% Time of traces
use_t = frame_t;
framerate = 1./median(diff(use_t));

% Power spectrum of MUA
use_trace = frame_spikes(1:end-1);

[fluor_spect,spect_f,spect_t] = spectrogram(use_trace,window_length_samples, ...
    round(spect_overlap/100*window_length_samples),[],framerate);
fluor_spect_norm = (fluor_spect/framerate).*conj(fluor_spect/framerate);  % framerate is used to normalize the FFT amplitudes
fluor_spect_power = fluor_spect_norm*2*df;

% Power spectrum of fluorescence
use_trace = fluor_trace;

[mua_spect,spect_f,spect_t] = spectrogram(use_trace,window_length_samples, ...
    round(spect_overlap/100*window_length_samples),[],framerate);
mua_spect_norm = (mua_spect/framerate).*conj(mua_spect/framerate);  % framerate is used to normalize the FFT amplitudes
mua_spect_power = mua_spect_norm*2*df;

% Correlate power spectra of MUA and fluorescence
spectra_corr = mat2cell(corrcoef([mua_spect_power',fluor_spect_power']),...
    repmat(length(spect_f),1,2),repmat(length(spect_f),1,2));

figure; colormap(hot);
c = [min(reshape(cell2mat(spectra_corr),[],1)), ...
    max(reshape(cell2mat(spectra_corr),[],1))];

subplot(1,3,1); 
imagesc(spect_f,spect_f,spectra_corr{1,1});
xlabel('MUA frequency');
ylabel('MUA frequency');
axis square; caxis(c);

subplot(1,3,2); 
imagesc(spect_f,spect_f,spectra_corr{2,2});
xlabel('Fluor frequency');
ylabel('Fluor frequency');
axis square; caxis(c);

subplot(1,3,3); 
imagesc(spect_f,spect_f,spectra_corr{1,2});
xlabel('Fluor frequency');
ylabel('MUA frequency');
axis square; caxis(c);












