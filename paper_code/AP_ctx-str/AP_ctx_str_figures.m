%% Fig 1b: Example traces

warning('Probably not best example, check others');

% Load and align
str_align = 'kernel';
animal = 'AP028'; 
day = '2017-12-16'; 
experiment = 1; 
verbose = false; 
AP_load_experiment;

avg_im_aligned = AP_align_widefield(animal,day,avg_im);
Udf_aligned = single(AP_align_widefield(animal,day,Udf));

% Define ROIs and get fluorescence traces
roi_circle_size = 20;
roi_x = [131,174,110,51];
roi_y = [297,96,71,144];
[x,y] = meshgrid(1:size(avg_im_aligned,1),1:size(avg_im_aligned,2));
roi_mask = cell2mat(arrayfun(@(roi) sqrt((x-roi_x(roi)).^2 + (y-roi_y(roi)).^2) <= ...
    roi_circle_size,permute(1:length(roi_x),[1,3,2]),'uni',false));
roi_trace = AP_svd_roi(Udf_aligned,fVdf,[],[],roi_mask);

roi_trace_deriv = diff(roi_trace,[],2);
roi_trace_deriv(roi_trace_deriv < 0) = 0;
frame_t_deriv = conv(frame_t,[1,1]/2,'valid');

% Bin spikes by aligned depth
n_depths = n_aligned_depths;
depth_group = aligned_str_depth_group;

time_bins = [frame_t_deriv,frame_t_deriv(end)+1/framerate];
binned_spikes = zeros(n_depths,length(time_bins)-1);
for curr_depth = 1:n_depths
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
end

% Plot ROIs and traces
figure;
subplot(1,6,1);

roi_boundaries = bwboundaries(sum(roi_mask,3));
imagesc(avg_im_aligned);colormap(gray);
caxis([0,prctile(avg_im_aligned(:),99)]);
axis image;
AP_reference_outline('ccf_aligned','r');
p = cellfun(@(x) plot(x(:,2),x(:,1),'b','linewidth',2),roi_boundaries);

subplot(1,6,2:6); hold on;
p1 = AP_stackplot(bsxfun(@rdivide,binned_spikes,std(binned_spikes,[],2))', ...
    frame_t_deriv,10,false,'k');
p2 = AP_stackplot(bsxfun(@rdivide,roi_trace_deriv,std(roi_trace_deriv,[],2))', ...
    frame_t_deriv,10,false,[0,0.7,0]);
xlabel('Time (seconds)');
ylabel('Activity (std)');
legend([p1(1),p2(1)],{'MUA','\DeltaFluorescence'});
xlim([177,200]);

%% Fig 1b: Average regression maps




































