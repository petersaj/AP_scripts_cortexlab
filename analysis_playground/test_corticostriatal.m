%% Testing for corticostriatal-specific imaging

%% Get explained variance in sliding window

% Load data
animal = 'AP063';
day = '2020-03-15';
experiments = AP_list_experiments(animal,day);
experiment = experiments(find(contains({experiments.protocol},'AP_lcrGratingPassive'),1,'last')).experiment;
verbose = true;
AP_load_experiment;


% Get striatum MUA in sliding window
skip_seconds = 60;
spike_binning_t = 1/framerate; % seconds
spike_binning_t_edges = frame_t(1)+skip_seconds:spike_binning_t:frame_t(end)-skip_seconds;
spike_binning_t_centers = spike_binning_t_edges(1:end-1) + diff(spike_binning_t_edges)/2;

mua_window = 500; % MUA window in microns
mua_window_spacing = 250; % MUA window spacing in microns

mua_bins = [str_depth(1):mua_window_spacing:(str_depth(2)-mua_window); ...
    (str_depth(1):mua_window_spacing:(str_depth(2)-mua_window))+mua_window];
mua_bin_centers = mua_bins(1,:) + diff(mua_bins,[],1)/2;
n_depths = length(mua_bin_centers);

striatum_mua = zeros(size(mua_bins,2),length(spike_binning_t_edges)-1);
for curr_depth = 1:size(mua_bins,2)
    curr_depth_templates_idx = ...
        find(template_depths >= mua_bins(1,curr_depth) & ...
        template_depths < mua_bins(2,curr_depth));
    
    striatum_mua(curr_depth,:) = histcounts(spike_times_timeline( ...
        ismember(spike_templates,curr_depth_templates_idx)),spike_binning_t_edges);
end

striatum_mua_std = striatum_mua./nanstd(striatum_mua,[],2);
striatum_mua_std(isnan(striatum_mua_std)) = 0;



use_svs = 1:200;
kernel_t = [-0.5,0.5];
kernel_frames = round(kernel_t(1)*framerate):round(kernel_t(2)*framerate);
lambda = 20;
zs = [false,false];
cvfold = 5;
return_constant = false;
use_constant = true;
        
fVdf_deconv = AP_deconv_wf(fVdf);
fVdf_deconv(isnan(fVdf_deconv)) = 0;
fVdf_deconv_resample = interp1(frame_t,fVdf_deconv(use_svs,:)',spike_binning_t_centers)';
        
[k,striatum_mua_std_predicted,explained_var] = ...
    AP_regresskernel(fVdf_deconv_resample, ...
    striatum_mua_std,kernel_frames,lambda,zs,cvfold,return_constant,use_constant);

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
r_px_max_norm = bsxfun(@rdivide,r_px_max, ...
    permute(max(reshape(r_px_max,[],n_depths),[],1),[1,3,2]));
r_px_max_norm(isnan(r_px_max_norm)) = 0;
r_px_com = sum(bsxfun(@times,r_px_max_norm,permute(1:n_depths,[1,3,2])),3)./sum(r_px_max_norm,3);

% Plot map of cortical pixel by preferred depth of probe
r_px_com_col = ind2rgb(round(mat2gray(r_px_com,[1,n_depths])*255),jet(255));
figure;
a1 = axes('YDir','reverse');
imagesc(avg_im); colormap(gray); caxis([0,prctile(avg_im(:),99)]);
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
set(c2,'YTickLabel',round(linspace(mua_bins(1),mua_bins(end),n_depths)));

% Plot MUA correlation and explained variance by cortex
figure;
subplot(4,1,1); 
plot(mua_bin_centers,explained_var.total,'k','linewidth',2);
line(xlim,[0,0],'color','r');

subplot(4,1,2:4);
imagesc(mua_bin_centers,mua_bin_centers,corrcoef(striatum_mua_std'));
colormap(hot);
caxis([0,1]);


% Plot measured and predicted stim-aligned activity

% Set options
surround_window = [-0.5,3];
baseline_window = [-0.1,0];

surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);
baseline_surround_time = baseline_window(1):surround_samplerate:baseline_window(2);

% (passive)
% use_stims = find(stimIDs == 3);
% (choiceworld)
stimIDs = trial_conditions(:,1).*trial_conditions(:,2);
use_stims = find(stimIDs > 0);

use_stimOn_times = stimOn_times(use_stims);

stim_surround_times = bsxfun(@plus, use_stimOn_times(:), surround_time);
stim_baseline_surround_times = bsxfun(@plus, use_stimOn_times(:), baseline_surround_time);

psth_measured_stim = permute(interp1(spike_binning_t_centers, ...
    striatum_mua_std',stim_surround_times),[3,2,1]);
psth_predicted_stim = permute(interp1(spike_binning_t_centers, ...
    striatum_mua_std_predicted',stim_surround_times),[3,2,1]);

psth_measured_baseline = permute(nanmean(interp1(spike_binning_t_centers, ...
    striatum_mua_std',stim_baseline_surround_times),2),[3,2,1]);
psth_predicted_baseline = permute(nanmean(interp1(spike_binning_t_centers, ...
    striatum_mua_std_predicted',stim_baseline_surround_times),2),[3,2,1]);

psth_measured = nanmean(psth_measured_stim - psth_measured_baseline,3);
psth_predicted = nanmean(psth_predicted_stim - psth_predicted_baseline,3);

figure; hold on;
AP_stackplot(psth_measured',surround_time,2,false,'k');
AP_stackplot(psth_predicted',surround_time,2,false,'r');
line([0,0],ylim,'color','k','linestyle','--');



%% Load and average stim-aligned image

animals = {'AP063','AP064','AP066','AP068'};
day = '2020-02-03';
experiment = 1;
verbose = true; 

im_stim_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    AP_load_experiment;
    
    % Set options
    surround_window = [-0.5,3];
    baseline_window = [-0.1,0];
    
    surround_samplerate = 1/(framerate*1);
    surround_time = surround_window(1):surround_samplerate:surround_window(2);
    baseline_surround_time = baseline_window(1):surround_samplerate:baseline_window(2);
    
    % Average (time course) responses
    use_vs = 1:size(U,3);
    
    conditions = unique(stimIDs);
    im_stim = nan(size(U,1),size(U,2),length(surround_time),length(conditions));
    for curr_condition_idx = 1:length(conditions)
        curr_condition = conditions(curr_condition_idx);
        
        use_stims = find(stimIDs == curr_condition);
        use_stimOn_times = stimOn_times(use_stims(2:end));
        use_stimOn_times([1,end]) = [];
        
        stim_surround_times = bsxfun(@plus, use_stimOn_times(:), surround_time);
        stim_baseline_surround_times = bsxfun(@plus, use_stimOn_times(:), baseline_surround_time);
        
        peri_stim_v = permute(interp1(frame_t,fVdf',stim_surround_times),[3,2,1]);
        baseline_v = permute(nanmean(interp1(frame_t,fVdf',stim_baseline_surround_times),2),[3,2,1]);
        
        stim_v_mean = nanmean(peri_stim_v - baseline_v,3);
        
        im_stim(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf(:,:,use_vs), ...
            stim_v_mean(use_vs,:));
    end
    
    im_stim_all{curr_animal} = im_stim;
    
end

im_stim_cat = cat(5,im_stim_all{:});
AP_image_scroll(nanmean(im_stim_cat(:,:,:,2,:),5));
axis image










