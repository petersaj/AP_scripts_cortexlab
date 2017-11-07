%% Fluorescence ROIs and striatal MUAs: example traces and regression

%animal = 'AP025'; day = '2017-10-01'; experiment = 3; AP_load_experiment;

% (Get the striatum map in 3 depths before this)

r_px_max_rescale = mat2gray(r_px_max,[0,2e-11]);

[roi_trace1,roi_mask1] = AP_svd_roi(Udf,fVdf,avg_im,r_px_max_rescale(:,:,1));
[roi_trace2,roi_mask2] = AP_svd_roi(Udf,fVdf,avg_im,r_px_max_rescale(:,:,2));
[roi_trace3,roi_mask3] = AP_svd_roi(Udf,fVdf,avg_im,r_px_max_rescale(:,:,3));

roi_trace = [roi_trace1;roi_trace2;roi_trace3];

frametime = median(diff(frame_t));
roi_trace_df = interp1(frame_t(1:end-1)+frametime/2,diff(roi_trace,[],2)',frame_t)';
roi_trace_df(:,1) = 0;
roi_trace_df(:,end) = 0;

% Plot the ROIs
figure;
imagesc(avg_im); caxis([0,40000]);
hold on; axis image off;
colormap(gray);
roi_boundaries = bwboundaries(roi_mask1+roi_mask2+roi_mask3);
for i = 1:length(roi_boundaries)
    plot(roi_boundaries{i}(:,2),roi_boundaries{i}(:,1),'m','linewidth',3);
end

% Group multiunit by depth
n_depth_groups = 3;
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
depth_group_edges_use = depth_group_edges;

[depth_group_n,depth_group] = histc(spikeDepths,depth_group_edges_use);
depth_groups_used = unique(depth_group);
depth_group_centers = depth_group_edges_use(1:end-1)+(diff(depth_group_edges_use)/2);

time_bins = [frame_t,frame_t(end)+frametime];

binned_spikes = zeros(length(depth_group_edges_use)-1,length(time_bins)-1);
for curr_depth = 1:length(depth_group_edges_use)-1
    
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
%     curr_spike_times = spike_times_timeline((depth_group == curr_depth) & ...
%         ismember(spike_templates,find(msn)));
    
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
    
end

smooth_size = 35;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
binned_spikes_smooth = conv2(binned_spikes,smWin,'same');

t_plot = [143,173];
t_use = frame_t > t_plot(1) & frame_t < t_plot(2);

% Plot example traces
figure;
p1 = subplot(1,3,1);
AP_stackplot(roi_trace(:,t_use)',frame_t(t_use),5,true,[0,0.7,0]);
xlabel('Time (s)');
ylabel('F')

p2 = subplot(1,3,2);
roi_trace_df_smooth = conv2(roi_trace_df,smWin,'same');
AP_stackplot(roi_trace_df_smooth(:,t_use)',frame_t(t_use),5,true,[0,0.7,0]);
xlabel('Time (s)');
ylabel('dF')

p3 = subplot(1,3,3);
AP_stackplot(binned_spikes_smooth(:,t_use)',frame_t(t_use),5,true,'k');
xlabel('Time (s)');
ylabel('Spikes')

linkaxes([p1,p2,p3],'x');
axis tight;

% Use ROI fluorescence to predict spikes
kernel_frames = -17:17;
lambda = 10;
zs = [true,true];
cvfold = 5;

[k_f,predicted_spikes_f,explained_var_f] = ...
    AP_regresskernel(roi_trace, ...
    binned_spikes(2,:),kernel_frames,lambda,zs,cvfold);
predicted_spikes_f_smooth = conv2(predicted_spikes_f,smWin,'same');
k_f_r = reshape(k_f,size(roi_trace,1),[]);

[k_df,predicted_spikes_df,explained_var_df] = ...
    AP_regresskernel(roi_trace_df, ...
    binned_spikes(2,:),kernel_frames,lambda,zs,cvfold);
k_df_r = reshape(k_df,size(roi_trace,1),[]);
predicted_spikes_df_smooth = conv2(predicted_spikes_df,smWin,'same');

% Plot the regression stuff
figure;
subplot(2,2,1);
plot(-kernel_frames/framerate,k_f_r','linewidth',2);
line([0,0],ylim,'color','k');
xlabel('Time from spike (s)');
ylabel('Weight');
title('F regression')
axis tight;
legend({'ROI 1','ROI 2','ROI 3'});

subplot(2,2,2);
plot(-kernel_frames/framerate,k_df_r','linewidth',2);
line([0,0],ylim,'color','k');
xlabel('Time from spike (s)');
ylabel('Weight');
title('dF regression')
axis tight;
legend({'ROI 1','ROI 2','ROI 3'});

subplot(2,2,3); hold on;
plot(frame_t(t_use),zscore(binned_spikes_smooth(2,t_use)),'k','linewidth',2)
plot(frame_t(t_use),zscore(predicted_spikes_f_smooth(t_use)),'r','linewidth',2)
xlabel('Time (s)');
ylabel('Spikes');
legend({'Real','Predicted'});
axis tight;

subplot(2,2,4); hold on;
plot(frame_t(t_use),zscore(binned_spikes_smooth(2,t_use)),'k','linewidth',2)
plot(frame_t(t_use),zscore(predicted_spikes_df_smooth(t_use)),'r','linewidth',2)
xlabel('Time (s)');
ylabel('Spikes');
legend({'Real','Predicted'});
axis tight;

% Scatter plot (/heatmap?) of predicted vs. real
predicted_spikes_rescale = (predicted_spikes_df*std(binned_spikes(2,:)) + mean(binned_spikes(2,:)));
figure; hold on;

[pred_spikes_grp,pred_spikes_quartiles] = grpstats(predicted_spikes_rescale,binned_spikes(2,:),{'gname',@(x) prctile(x,[25;50;75])});

pred_spikes_grp = cellfun(@str2num,pred_spikes_grp);
plot(pred_spikes_grp,pred_spikes_quartiles(:,2),'r');
plot([pred_spikes_grp,fliplr(pred_spikes_grp)],[pred_spikes_quartiles(:,1),fliplr(pred_spikes_quartiles(:,3))],'b');

min_spikes = min([binned_spikes(2,:),predicted_spikes_rescale]);
max_spikes = max([binned_spikes(2,:),predicted_spikes_rescale]);
xlim([min_spikes,max_spikes]);
ylim([min_spikes,max_spikes]);
line([0,max_spikes],[0,max_spikes],'color','k');
xlabel('Real spikes');
ylabel('Predicted spikes');
axis square;


%% Example regression weights in time
% (Get the striatum map in 3 depths before this)

plot_weights = r_px(:,:,end:-1:1,2);
kernel_frame_times = -fliplr(kernel_frames)*downsample_factor/framerate;

plot_times = [-0.1,0.1];
plot_kernel_frames = kernel_frame_times >= plot_times(1) & kernel_frame_times <= plot_times(2);
plot_weights_reshape = reshape(plot_weights(:,:,plot_kernel_frames),size(plot_weights,1),[]);
figure;imagesc(plot_weights_reshape);
caxis([-prctile(plot_weights(:),99.9),prctile(plot_weights(:),99.9)])
colormap(colormap_BlueWhiteRed); colorbar;
axis image off; 
title(sprintf('%0.5g,',round(1000*kernel_frame_times(plot_kernel_frames))/1000));

figure
imagesc(r_px_max(:,:,2));
caxis([-prctile(plot_weights(:),99.9),prctile(plot_weights(:),99.9)])
colormap(colormap_BlueWhiteRed);
axis image off; 


%% Example map
% (Get the striatum map in 6 depths before this)

% Plot the maximum weights by depth
r_px_max_rescale = mat2gray(r_px_max,[0,3e-11]);
r_px_max_rescale_reshape = reshape(permute(r_px_max_rescale,[1,3,2]),[],size(r_px_max,2));
figure;
imagesc(r_px_max_rescale_reshape);
caxis([-1 1]);
colormap(colormap_BlueWhiteRed);
axis image off;

% Plot the explained variance by depth 
figure;plot(explained_var.total,depth_group_centers,'k','linewidth',2)
set(gca,'YDir','reverse')
xlabel('Explained variance')
ylabel('Depth')
axis tight;

% Plot the mean firing rate and variance by depth
spike_rate = nanmean(binned_spikes,2);
spike_var = var(binned_spikes,[],2);
figure; hold on
p = plotyy(depth_group_centers,spike_rate,depth_group_centers,spike_var);
axis(p(1),'tight')
axis(p(2),'tight')
xlabel('Depth');
ylabel(p(1),'Rate');
ylabel(p(2),'Variance');

%% Maps of example cells
% 6 (tan), 8/95 (msn), 53/119 (fsi)
% these are saved, easier enough to just produce on the fly



