% Generate figures for ctx-str paper

% Anything that takes a lot of time is done in
% AP_ctx_str_trial_preprocessing and saved for plotting here

% The original scripts here were in test_wf_ephys_choiceworld_analysis

%% Fig 1 b/c: Cortical activity during task

data_fn = ['trial_activity_choiceworld'];
exclude_data = true;

[fluor_allcat,mua_allcat,wheel_allcat,reward_allcat,D_allcat] = ...
    AP_load_concat_normalize_ctx_str(data_fn,exclude_data);

n_vs = size(fluor_allcat,3);

% Get trial information
trial_contrast_allcat = max(D_allcat.stimulus,[],2);
[~,side_idx] = max(D_allcat.stimulus > 0,[],2);
trial_side_allcat = (side_idx-1.5)*2;
trial_choice_allcat = -(D_allcat.response-1.5)*2;

% Get time (make this be saved in trial data)
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

%%% Plot map of cortical activity at time points on correct visual trials
rxn_time_use = [0.1,0.3];

vis_correct_L_trials = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == -1 & ...
    move_t >= rxn_time_use(1) & ...
    move_t <= rxn_time_use(2);

vis_correct_L_px = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    diff(squeeze(nanmean(fluor_allcat(vis_correct_L_trials,:,:),1))',[],2));

t_diff = conv(t,[1,1]/2,'valid');
AP_image_scroll(vis_correct_L_px,t_diff);
axis image;
caxis([0,max(abs(caxis))])
colormap(brewermap([],'BuGn'));
AP_reference_outline('ccf_aligned','k');

plot_t = [find(t_diff > 0.06,1),find(t_diff > 0.120,1), ...
    find(t_diff > 0.3,1),find(t_diff > 0.8,1)];

figure;
colormap(brewermap([],'BuGn'));
for curr_t = 1:length(plot_t)
    subplot(1,length(plot_t),curr_t);
    imagesc(vis_correct_L_px(:,:,plot_t(curr_t)));
    caxis([0,max(abs(caxis))]);
    axis image off;
    AP_reference_outline('ccf_aligned','k');
    title([num2str(round(t_diff(plot_t(curr_t))*1000)) ' ms after stim']);
end

%%% Plot ROI PSTH's
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];

contrast_side_col = colormap_BlueWhiteRed(5);
contrast_side_col(6,:) = 0;
contrast_side_val = unique(sort([-contrasts,contrasts]))';

plot_conditions = ...
    [contrasts,contrasts; ...
    -ones(1,6),-1,ones(1,5); ...
    ones(1,6),-ones(1,6)]';

[~,plot_id] = ismember( ...
    [trial_contrast_allcat,trial_side_allcat,trial_choice_allcat], ...
    plot_conditions,'rows');

% Get fluorescence in widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);
n_rois = numel(wf_roi);

roi_mask = cat(3,wf_roi.mask);

fluor_roi = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat,[3,2,1]),n_vs,[]),[],[],roi_mask), ...
    size(roi_mask,3),[],size(fluor_allcat,1)),[3,2,1]);

% (diff, > 0, normalize)
fluor_roi_diff = diff(fluor_roi,[],2);
fluor_roi_diff(fluor_roi_diff < 0) = 0;
fluor_roi_diff = bsxfun(@rdivide,fluor_roi_diff,nanstd(reshape(fluor_roi_diff,[],1,size(wf_roi,1)),[],1));

t_diff =  conv(t,[1,1]/2,'valid');

figure;
for curr_align = 1:2
    p(curr_align) = subplot(1,2,curr_align); hold on;
    for curr_plot_condition = 1:size(plot_conditions,1)
        
        curr_trials = plot_id == curr_plot_condition & ...
            move_t >= rxn_time_use(1) & ...
            move_t <= rxn_time_use(2);

        curr_data = fluor_roi_diff(curr_trials,:,:);
        
        if curr_align == 2
            % Re-align to movement onset
            t_leeway = -t_diff(1);
            leeway_samples = round(t_leeway*sample_rate);
            curr_move_idx = move_idx(curr_trials);
            for i = 1:size(curr_data,1)
                curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
            end
        end
        
        curr_data_mean = squeeze(nanmean(curr_data,1));
        
        curr_contrast_side = max(trial_contrast_allcat(curr_trials))*sign(max(trial_side_allcat(curr_trials)));
        contrast_side_idx = find(curr_contrast_side == contrast_side_val);
        curr_col = contrast_side_col(contrast_side_idx,:);
                
        if curr_contrast_side == 0
            switch max(trial_choice_allcat(curr_trials))
                case -1
                    curr_col = 'm';
                case 1
                    curr_col = 'c';
            end
        end
        
        AP_stackplot(curr_data_mean,t_diff,3.5,false,curr_col,{wf_roi.area},true);
        
    end
    line([0,0],ylim,'color','k');
    switch curr_align
        case 1
            xlabel('Time from stim');
        case 2
            xlabel('Time from move');
    end   
end
linkaxes(p,'y');

% Plot widefield ROIs
figure; hold on
set(gca,'YDir','reverse');
AP_reference_outline('ccf_aligned','k');
for curr_roi = 1:n_rois
    curr_roi_boundary = cell2mat(bwboundaries(roi_mask(:,:,curr_roi)));
    patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),[0,0.7,0]);
    
    text(nanmean(curr_roi_boundary(:,2)),nanmean(curr_roi_boundary(:,1)), ...
        wf_roi(curr_roi).area,'FontSize',10,'HorizontalAlignment','center','color','w')
end
axis image off;

%% Fig 1 d/e: Task -> cortex regression

% Load regression results
regression_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\trial_activity_choiceworld_task_regression';
load(regression_fn);

% Get number of V's/depths
n_vs = size(fluor_allcat_predicted,3);
n_depths = size(mua_allcat_predicted,3);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Get move onset index
[~,move_idx] = max(abs(wheel(:,:,1)) > 2,[],2);
move_t = t_downsample_diff(move_idx)';

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Set which components to use
use_svs = 200;

% Get time shifts in samples
sample_shifts = cellfun(@(x) round(x(1)*(sample_rate/downsample_factor)): ...
    round(x(2)*(sample_rate/downsample_factor)),t_shifts,'uni',false);

% Plot fluorescence kernels
for curr_regressor = 1:length(regressors)
    curr_k_v = permute(cell2mat(permute(fluor_kernel(curr_regressor,1:use_svs),[1,3,2])),[3,2,1]);
    curr_k_px = reshape(svdFrameReconstruct(U_master(:,:,1:use_svs), ...
        reshape(curr_k_v,use_svs,[])),size(U_master,1),size(U_master,2),[],size(curr_k_v,3));
    AP_image_scroll(curr_k_px,sample_shifts{curr_regressor}/(sample_rate/downsample_factor));
    axis image
    caxis([0,prctile(abs(curr_k_px(:)),100)]);
    colormap(brewermap([],'BuGn'))
    AP_reference_outline('ccf_aligned','k');
    set(gcf,'Name',regressor_labels{curr_regressor});
end

% Plot fluorescence constant
curr_k_v = permute(cell2mat(permute(fluor_kernel(length(regressors)+1,1:use_svs),[1,3,2])),[3,2,1]);
curr_k_px = reshape(svdFrameReconstruct(U_master(:,:,1:use_svs), ...
    reshape(curr_k_v,use_svs,[])),size(U_master,1),size(U_master,2),[],size(curr_k_v,3));
figure;imagesc(curr_k_px);
axis image off
caxis([-prctile(abs(curr_k_px(:)),100),prctile(abs(curr_k_px(:)),100)]);
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k');
title('Constant');

% Plot fluorescence ROI kernels
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);
n_rois = numel(wf_roi);

roi_mask = cat(3,wf_roi.mask);

U_roi = bsxfun(@rdivide,transpose(reshape(U_master,[],size(U_master,3))'* ...
    reshape(roi_mask,[],size(roi_mask,3))), ...
    sum(reshape(roi_mask,[],size(roi_mask,3)),1)');

for curr_regressor = 1:length(regressors)
    curr_k_v = permute(cell2mat(permute(fluor_kernel(curr_regressor,1:use_svs),[1,3,2])),[3,2,1]);
    
    curr_k_roi = reshape(AP_svd_roi(U_master(:,:,1:n_vs), ...
        reshape(curr_k_v,n_vs,[]),[],[],roi_mask)', ...
        size(curr_k_v,2),size(curr_k_v,3),n_rois);
    
    if size(curr_k_roi,2) > 1
        curr_col = colormap_BlueWhiteRed(size(curr_k_roi,2)/2);
        curr_col(size(curr_k_roi,2)/2+1,:) = [];
    else
        curr_col = 'k';
    end
    
    figure; hold on;
    for curr_subk = 1:size(curr_k_roi,2)
        AP_stackplot(squeeze(curr_k_roi(:,curr_subk,:)), ...
            sample_shifts{curr_regressor}/(sample_rate/downsample_factor), ...
            range(curr_k_roi(:))*1.2,false,curr_col(curr_subk,:),{wf_roi.area},true);
    end
    line([0,0],ylim,'color','k');
    title(regressor_labels{curr_regressor});
end

% Plot ROI actual and predicted
rxn_time_use = [0.1,0.3];
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];

contrast_side_col = colormap_BlueWhiteRed(5);
contrast_side_col(6,:) = 0;
contrast_side_val = unique(sort([-contrasts,contrasts]))';

plot_conditions = ...
    [contrasts,contrasts; ...
    -ones(1,6),-1,ones(1,5); ...
    ones(1,6),-ones(1,6)]';

[~,plot_id] = ismember( ...
    [trial_contrast_allcat,trial_side_allcat,trial_choice_allcat], ...
    plot_conditions,'rows');

% Get real/predicted fluorescence in widefield ROIs
fluor_downsample_diff_roi = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_downsamp_diff,[3,2,1]),n_vs,[]),[],[],roi_mask), ...
    size(roi_mask,3),[],size(fluor_allcat_downsamp_diff,1)),[3,2,1]);

fluor_predicted_roi = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_predicted,[3,2,1]),n_vs,[]),[],[],roi_mask), ...
    size(roi_mask,3),[],size(fluor_allcat_predicted,1)),[3,2,1]);

figure;
for curr_plot = 1:3
    
    switch curr_plot
        case 1
            plot_data = fluor_downsample_diff_roi;
            plot_title = 'Measured';
        case 2
            plot_data = fluor_predicted_roi;
            plot_title = 'Predicted';
        case 3
            plot_data = fluor_downsample_diff_roi - fluor_predicted_roi;
            plot_title = 'Residual';
    end
    
    p_ctx(curr_plot) = subplot(1,3,curr_plot); hold on;
    for curr_plot_condition = 1:size(plot_conditions,1)
        
        curr_trials = plot_id == curr_plot_condition & ...
            move_t >= rxn_time_use(1) & ...
            move_t <= rxn_time_use(2);
        curr_data = plot_data(curr_trials,:,:);
        
        % Re-align to movement onset
        t_leeway = -t_downsample_diff(1);
        leeway_samples = round(t_leeway*(sample_rate/downsample_factor));
        curr_move_idx = move_idx(curr_trials);
        for i = 1:size(curr_data,1)
            curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
        end
        
        curr_data_mean = squeeze(nanmean(curr_data,1));
        
        curr_contrast_side = max(trial_contrast_allcat(curr_trials))*sign(max(trial_side_allcat(curr_trials)));
        contrast_side_idx = find(curr_contrast_side == contrast_side_val);
        curr_col = contrast_side_col(contrast_side_idx,:);
        
        if curr_contrast_side == 0
            switch max(trial_choice_allcat(curr_trials))
                case -1
                    curr_col = 'm';
                case 1
                    curr_col = 'c';
            end
        end
        
        AP_stackplot(curr_data_mean,t_downsample_diff,6e-3,false,curr_col,{wf_roi.area});
        
    end
    line([0,0],ylim,'color','k');
    xlabel('Time from stim');
    title(plot_title);   
end
linkaxes(p_ctx)

% Get max stim kernel and contrast slope maps
stim_regressor = strcmp(regressor_labels,'Stim');
k_v_stim = permute(cell2mat(permute(fluor_kernel(stim_regressor,1:use_svs),[1,3,2])),[3,2,1]);
k_px_stim = reshape(svdFrameReconstruct(U_master(:,:,1:use_svs), ...
    reshape(k_v_stim,use_svs,[])),size(U_master,1),size(U_master,2),[],size(k_v_stim,3));
k_px_stim_maxt = squeeze(max(k_px_stim,[],3));

contrasts = [0.06,0.125,0.25,0.5,1];

contrast_k_L = AP_regresskernel(contrasts, ...
    reshape(k_px_stim_maxt(:,:,5:-1:1),[],5),0,[],[false,false],1,true,true);
contrast_slope_L = reshape(contrast_k_L{1},size(U_master,1),size(U_master,2));

contrast_k_R = AP_regresskernel(contrasts, ...
    reshape(k_px_stim_maxt(:,:,6:10),[],5),0,[],[false,false],1,true,true);
contrast_slope_R = reshape(contrast_k_R{1},size(U_master,1),size(U_master,2));

contrast_slope_diff = contrast_slope_L - contrast_slope_R;

% Get max move kernel and move ipsi-contra differences
move_onset_regressor = strcmp(regressor_labels,'Move onset');
k_v_move_onset = permute(cell2mat(permute(fluor_kernel(move_onset_regressor,1:use_svs),[1,3,2])),[3,2,1]);
k_px_move_onset = reshape(svdFrameReconstruct(U_master(:,:,1:use_svs), ...
    reshape(k_v_move_onset,use_svs,[])),size(U_master,1),size(U_master,2),[],size(k_v_move_onset,3));
k_px_move_onset_maxt = squeeze(max(k_px_move_onset,[],3));

k_px_move_onset_maxt_diff = k_px_move_onset_maxt(:,:,1)-k_px_move_onset_maxt(:,:,2);

figure;
subplot(2,2,1);
imagesc(max(k_px_stim_maxt,[],3))
axis image off;
AP_reference_outline('ccf_aligned','k');
caxis([0,max(k_px_stim_maxt(:))])
colormap(gca,brewermap([],'Purples'));
title('Max stim kernel weight');

subplot(2,2,2);
imagesc(contrast_slope_diff); 
axis image off;
AP_reference_outline('ccf_aligned','k');
caxis([-max(k_px_stim_maxt(:)),max(k_px_stim_maxt(:))])
colormap(gca,brewermap([],'RdBu'));
title('Contrast slope R-L');

subplot(2,2,3);
imagesc(max(k_px_move_onset_maxt,[],3));
axis image off;
AP_reference_outline('ccf_aligned','k');
caxis([0,max(k_px_move_onset_maxt(:))])
colormap(gca,brewermap([],'Purples'));
title('Max move onset kernel weight');

subplot(2,2,4);
imagesc(k_px_move_onset_maxt_diff);
axis image off;
AP_reference_outline('ccf_aligned','k');
caxis([-max(k_px_move_onset_maxt(:)),max(k_px_move_onset_maxt(:))])
colormap(gca,brewermap([],'*RdBu'));
title('Move L-R');


% Plot contrast response function in ROIs 
contrast_sides = sort(reshape([0.06,0.125,0.25,0.5,1].*[-1,1]',[],1));

stim_max_k_roi = cell2mat(arrayfun(@(x) ....
    nanmean(reshape(reshape(repmat(roi_mask(:,:,x),1,1,10).* ...
    k_px_stim_maxt,[],10),[],10),1),1:size(roi_mask,3),'uni',false)');

figure;
AP_stackplot(stim_max_k_roi',contrast_sides, ...
    range(stim_max_k_roi(:))*1.2,false,[0,0.7,0],{wf_roi.area},true);
xlabel('Contrast*Side');
title('Stim kernel maximum');




%% %%%%%%%%%%%%%%% UP-TO-DATE ABOVE %%%%%%%%%%%%%%%%%%



%% Fig 1a: Example average widefield

animal = 'AP028';
day = '2017-12-16'; 
[img_path,img_exists] = AP_cortexlab_filename(animal,day,[],'imaging');
avg_im = readNPY([img_path filesep 'meanImage_blue.npy']);


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
% (half-done: variable names etc will change when new code finished)

n_aligned_depths = 4;

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_ephys';

k_fn = [data_path filesep 'wf_ephys_maps_concat_' num2str(n_aligned_depths) '_depths_kernel'];
load(k_fn);

t = batch_vars(1).t{1};

% Concatenate and mean
% (kernel goes backwards in time - flip to correct)
k_px_trained_cat = cellfun(@(x) x(:,:,end:-1:1,:),[batch_vars(1:6).k_px],'uni',false);
k_px_naive_cat = cellfun(@(x) x(:,:,end:-1:1,:),[batch_vars(7:11).k_px],'uni',false);

k_px_trained = nanmean(double(cat(5,k_px_trained_cat{:})),5);
k_px_naive = nanmean(double(cat(5,k_px_naive_cat{:})),5);

AP_image_scroll(k_px_trained,t);
caxis([-max(abs(caxis)),max(abs(caxis))]);
axis image;
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k');
set(gcf,'Name','Trained');

AP_image_scroll(k_px_naive,t);
caxis([-max(abs(caxis)),max(abs(caxis))]);
axis image;
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k');
set(gcf,'Name','Naive');

% Get center-of-mass maps
k_px_trained_positive = k_px_trained;
k_px_trained_positive(k_px_trained_positive < 0) = 0;
k_px_trained_com = sum(k_px_trained_positive.*permute(1:n_aligned_depths,[1,3,4,2]),4)./sum(k_px_trained_positive,4);
k_px_trained_com_colored = nan(size(k_px_trained_com,1),size(k_px_trained_com,2),3,size(k_px_trained_com,3));
for curr_frame = 1:size(k_px_trained_com,3)
    k_px_trained_com_colored(:,:,:,curr_frame) = ...
        ind2rgb(round(mat2gray(k_px_trained_com(:,:,curr_frame),[1,n_aligned_depths])*255),jet(255));
end

k_px_naive_positive = k_px_naive;
k_px_naive_positive(k_px_naive_positive < 0) = 0;
k_px_naive_com = sum(k_px_naive_positive.*permute(1:n_aligned_depths,[1,3,4,2]),4)./sum(k_px_naive_positive,4);
k_px_naive_com_colored = nan(size(k_px_naive_com,1),size(k_px_naive_com,2),3,size(k_px_naive_com,3));
for curr_frame = 1:size(k_px_naive_com,3)
    k_px_naive_com_colored(:,:,:,curr_frame) = ...
        ind2rgb(round(mat2gray(k_px_naive_com(:,:,curr_frame),[1,n_aligned_depths])*255),jet(255));
end

% Whiten relative to weight, plot movie
k_px_trained_com_colored_weighted = k_px_trained_com_colored + ...
    ((1-(permute(max(abs(k_px_trained_positive),[],4)./ ...
    prctile(abs(k_px_trained_positive(:)),100),[1,2,4,3]))) .* ...
    (ones(1,1,3,length(t)) - k_px_trained_com_colored));

k_px_naive_com_colored_weighted = k_px_naive_com_colored + ...
    ((1-(permute(max(abs(k_px_naive_positive),[],4)./ ...
    prctile(abs(k_px_naive_positive(:)),100),[1,2,4,3]))) .* ...
    (ones(1,1,3,length(t)) - k_px_naive_com_colored));

AP_image_scroll(k_px_trained_com_colored_weighted,t,true);
axis image;
AP_reference_outline('ccf_aligned','k');
set(gcf,'Name','Trained');

AP_image_scroll(k_px_naive_com_colored_weighted,t,true);
axis image;
AP_reference_outline('ccf_aligned','k');
set(gcf,'Name','Naive');

% Plot at t == 0 transparency relative to weight
k_px_trained_t0 = squeeze(k_px_trained(:,:,t == 0,:));
k_px_trained_weight_norm = mat2gray(max(k_px_trained_t0,[],3), ...
    [0,prctile(reshape(max(k_px_trained_t0,[],3),[],1),95)]);

k_px_naive_t0 = squeeze(k_px_naive(:,:,t == 0,:));
k_px_naive_weight_norm = mat2gray(max(k_px_naive_t0,[],3), ...
    [0,prctile(reshape(max(k_px_naive_t0,[],3),[],1),95)]);

figure;
subplot(1,2,1);
p = imagesc(k_px_trained_com_colored(:,:,:,t == 0));
axis image off;
set(p,'AlphaData',k_px_trained_weight_norm);
AP_reference_outline('ccf_aligned','k');
title('Trained  (t = 0)');

subplot(1,2,2);
p = imagesc(k_px_naive_com_colored(:,:,:,t == 0));
axis image off;
set(p,'AlphaData',k_px_naive_weight_norm);
AP_reference_outline('ccf_aligned','k');
title('Naive (t = 0)');


%% Fig 1b: Allen projection maps vs regression maps?
% (maybe get average centroid of str 1/2/3/4 then get one map?)
% (copy code from AP_ctx2str_probe)

warning('This should probably be the average vector of wf-estimated');
probe_vector_ccf = [520,240,510;520,511,239];

%%% Get the average relative depth of each kernel template

% Load kernels by depths, get depth relative to maximum extent
kernel_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
kernel_fn = ['ephys_kernel_depth'];
load([kernel_path filesep kernel_fn])
total_depths = 1:max(cellfun(@(x) size(x,3),[ephys_kernel_depth.k_px]));
k_px_depths = cellfun(@(x) total_depths(end-size(x,3)+1:end),[ephys_kernel_depth.k_px],'uni',false);

% Load the kernel template matches
n_aligned_depths = 4;
kernel_match_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
kernel_match_fn = ['ephys_kernel_align_' num2str(n_aligned_depths) '_depths.mat'];
load([kernel_match_path filesep kernel_match_fn]);

% Concatenate all relative depths and kernel matches
k_depths = cell2mat(k_px_depths);
k_matches = cell2mat([ephys_kernel_align.kernel_match]')';
k_match_depths_relative = grpstats(k_depths,k_matches)./max(total_depths);

%%% Query allen at each point using targeted trajectory

% Load in the annotated Allen volume and names
allen_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
av = readNPY([allen_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean

% Get probe location per micron
probe_size = pdist2(probe_vector_ccf(1,:),probe_vector_ccf(2,:))*10;
probe_depths = ...
    round([linspace(probe_vector_ccf(1,1)',probe_vector_ccf(2,1)',probe_size); ...
    linspace(probe_vector_ccf(1,2)',probe_vector_ccf(2,2)',probe_size); ...
    linspace(probe_vector_ccf(1,3)',probe_vector_ccf(2,3)',probe_size)]');

% Eliminiate trajectory points that are off the atlas
eliminate_depths = ...
    probe_depths(:,1) < 1 | probe_depths(:,1) > size(av,1) | ...
    probe_depths(:,2) < 1 | probe_depths(:,2) > size(av,2) | ...
    probe_depths(:,3) < 1 | probe_depths(:,3) > size(av,3);
probe_depths(eliminate_depths,:) = [];

% Convert probe depths subscripts to indicies
probe_depths_ind = sub2ind(size(av),probe_depths(:,1),probe_depths(:,2),probe_depths(:,3));

% Get structures that the probe is in
probe_structures = av(probe_depths_ind);

% Get target relative depths through striatum
str_id = find(strcmp(st.safe_name,'Caudoputamen'));
probe_structures_str = probe_structures == str_id;
probe_str = [probe_depths(find(probe_structures_str,1,'first'),:); ...
    probe_depths(find(probe_structures_str,1,'last'),:)];
kernel_depth_ccf = interp1([0,1],probe_str,k_match_depths_relative);

%%%% Just use regular depths?
regular_centers_borders = linspace(0,1,n_aligned_depths*2+1);
kernel_depth_ccf = interp1([0,1],probe_str,regular_centers_borders(2:2:end));

% Plot brain to overlay probes
% (note the CCF is rotated to allow for dim 1 = x)
h = figure; ccf_axes = axes; hold on
slice_spacing = 10;
target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) > 1,[2,1,3]);
structure_patch = isosurface(target_volume,0);
structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
target_structure_color = [0.7,0.7,0.7];
brain_outline = patch('Vertices',structure_wire.vertices*slice_spacing, ...
    'Faces',structure_wire.faces, ...
    'FaceColor','none','EdgeColor',target_structure_color);

str_id = find(strcmp(st.safe_name,'Caudoputamen'));
target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) == str_id,[2,1,3]);
structure_patch = isosurface(target_volume,0);
structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
target_structure_color = [0.7,0,0.7];
striatum_outline = patch('Vertices',structure_wire.vertices*slice_spacing, ...
    'Faces',structure_wire.faces, ...
    'FaceColor','none','EdgeColor',target_structure_color);

axis image vis3d off;
view([-30,25]);
cameratoolbar(h,'SetCoordSys','y');
cameratoolbar(h,'SetMode','orbit');

scatter3(kernel_depth_ccf(:,1),kernel_depth_ccf(:,2),kernel_depth_ccf(:,3), ...
    100,copper(n_aligned_depths),'filled');

% Mirror all of the locations across the midline and get projections from
% both (because the cortical injections in the database aren't evenly
% distributed, so this is more comprehensive)
kernel_depth_um = round(kernel_depth_ccf*10);
bregma_um = allenCCFbregma*10;
ccf_midline = bregma_um(3);
hemisphere = sign(kernel_depth_um(1,3) - ccf_midline);
str_depths_mirror = [kernel_depth_um(:,1:2),ccf_midline - hemisphere*abs(kernel_depth_um(:,3)-ccf_midline)];

max_sites = 50;
str_depths_query = [kernel_depth_um;str_depths_mirror];
injection_parameters = get_allen_projection(str_depths_query,max_sites);
injection_coordinates = {injection_parameters.coordinates};

% Standardize injection coordinates by hemisphere (left = contra, right =
% ipsi)
injection_coordinates_standardized = injection_coordinates;
for curr_coord = 1:length(injection_coordinates)
    
    target_hemisphere = sign(ccf_midline - str_depths_query(curr_coord,3));
    injection_coords_ml_offset = abs(injection_coordinates{curr_coord}(:,3) - ccf_midline);
    injection_coordinates_hemisphere = sign(injection_coordinates{curr_coord}(:,3) - ccf_midline);
    
    injection_coords_ipsi = injection_coordinates_hemisphere == target_hemisphere;
    injection_coords_contra = injection_coordinates_hemisphere == -target_hemisphere;
    
    injection_coordinates_standardized{curr_coord}(injection_coords_ipsi,3) = ...
        ccf_midline + injection_coords_ml_offset(injection_coords_ipsi);
    injection_coordinates_standardized{curr_coord}(injection_coords_contra,3) = ...
        ccf_midline - injection_coords_ml_offset(injection_coords_contra);
    
end

% Get relative projection density / injection volumes
% projection_strength = cellfun(@(density,volume) density./volume, ...
%     {injection_parameters.density},{injection_parameters.volume},'uni',false);
% (or currently using: just the density, not sure whether good to norm)
projection_strength = cellfun(@(density,volume) density, ...
    {injection_parameters.density},{injection_parameters.volume},'uni',false);
projection_strength_normalized = cellfun(@(x) mat2gray(x, ...
    [min([projection_strength{:}]),max([projection_strength{:}])]), ...
    projection_strength,'uni',false);

% Convert points from CCF to widefield
ccf_tform_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\ccf_tform'];
load(ccf_tform_fn);

um2pixel = 20.6;
injection_coordinates_wf = cellfun(@(x) ...
    [x(:,[3,1]).*(1/(um2pixel)),ones(size(x,1),1)]*ccf_tform.T, ...
    injection_coordinates_standardized,'uni',false);

% Load kernel templates for overlay
n_aligned_depths = 4;
kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template_' num2str(n_aligned_depths) '_depths.mat'];
load(kernel_template_fn);

figure; 
colormap(brewermap([],'*RdBu'));
for curr_depth = 1:n_aligned_depths
    subplot(1,n_aligned_depths,curr_depth); hold on; axis image off;
    set(gca,'YDir','reverse');
    
    imagesc(kernel_template(:,:,curr_depth));
    caxis([-prctile(abs(kernel_template(:)),99),prctile(abs(kernel_template(:)),99)]);
    
    % (plot points from both hemispheres)
    scatter(injection_coordinates_wf{curr_depth}(:,1), ...
        injection_coordinates_wf{curr_depth}(:,2), ...
        projection_strength_normalized{curr_depth}*50 + 10, ...
        'k','filled');
    scatter(injection_coordinates_wf{n_aligned_depths + curr_depth}(:,1), ...
        injection_coordinates_wf{n_aligned_depths + curr_depth}(:,2), ...
        projection_strength_normalized{n_aligned_depths + curr_depth}*50 + 10, ...
        'k','filled');
    
    AP_reference_outline('ccf_aligned','k');
end


%% Fig 2a: Behavior psychometric 
% (from AP_vanillaChoiceworld_behavior - currently no eliminations)

% Load behavior
bhv_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\bhv_processing\bhv.mat';
load(bhv_fn);

animals = {bhv.animal};

conditions = unique(vertcat(bhv.conditions),'rows');
trial_choice_cat = arrayfun(@(x) horzcat(bhv(x).trial_choice{:}),1:length(bhv),'uni',false);
trial_outcome_cat = arrayfun(@(x) horzcat(bhv(x).trial_outcome{:}),1:length(bhv),'uni',false);
trial_side_cat = arrayfun(@(x) horzcat(bhv(x).trial_side{:}),1:length(bhv),'uni',false);
trial_contrast_cat = arrayfun(@(x) horzcat(bhv(x).trial_contrast{:}),1:length(bhv),'uni',false);
trial_condition_cat = cellfun(@(side,contrast) side.*contrast,trial_side_cat,trial_contrast_cat,'uni',false);
trial_wheel_velocity_cat = arrayfun(@(x) vertcat(bhv(x).trial_wheel_velocity{:})',1:length(bhv),'uni',false);
stim_to_move_cat = arrayfun(@(x) horzcat(bhv(x).stim_to_move{:}),1:length(bhv),'uni',false);
stim_to_feedback_cat = arrayfun(@(x) horzcat(bhv(x).stim_to_feedback{:}),1:length(bhv),'uni',false);

% Distinguish early/late movements
go_time = 0.5;
trial_timing = arrayfun(@(animal) cellfun(@(x) 1+(x > go_time), ...
    bhv(animal).stim_to_move,'uni',false),1:length(bhv),'uni',false);
trial_timing_cat = arrayfun(@(animal) ...
    horzcat(trial_timing{animal}{:}),1:length(bhv),'uni',false);

% Plot psychometric 
frac_left = cell2mat(cellfun(@(choice,condition) ...
    grpstats(choice == -1,condition),trial_choice_cat,trial_condition_cat,'uni',false));

frac_left_earlymove = cell2mat(cellfun(@(choice,condition,stim_to_move) ...
    grpstats(choice(stim_to_move < 0.5) == -1,condition(stim_to_move < 0.5)), ...
    trial_choice_cat,trial_condition_cat,stim_to_move_cat,'uni',false));

frac_left_latemove = cell2mat(cellfun(@(choice,condition,stim_to_move) ...
    grpstats(choice(stim_to_move >= 0.5) == -1,condition(stim_to_move >= 0.5)), ...
    trial_choice_cat,trial_condition_cat,stim_to_move_cat,'uni',false));

figure;

subplot(1,3,1); hold on; axis square;
plot(conditions,frac_left,'linewidth',2,'color',[0.5,0.5,0.5]);
plot(conditions,nanmean(frac_left,2),'linewidth',5,'color','k');
xlim([-1,1]);
ylim([0,1]);
line([0,0],ylim,'linestyle','--','color','k');
line(xlim,[0.5,0.5],'linestyle','--','color','k');
xlabel('Condition');
ylabel('Fraction go left');
title('All trials');

subplot(1,3,2); hold on; axis square;
plot(conditions,frac_left_earlymove,'linewidth',2,'color',[0.5,0.5,0.5]);
plot(conditions,nanmean(frac_left_earlymove,2),'linewidth',5,'color','k');
xlim([-1,1]);
ylim([0,1]);
line([0,0],ylim,'linestyle','--','color','k');
line(xlim,[0.5,0.5],'linestyle','--','color','k');
xlabel('Condition');
ylabel('Fraction go left');
title('Early move');

subplot(1,3,3); hold on; axis square;
plot(conditions,frac_left_latemove,'linewidth',2,'color',[0.5,0.5,0.5]);
plot(conditions,nanmean(frac_left_latemove,2),'linewidth',5,'color','k');
xlim([-1,1]);
ylim([0,1]);
line([0,0],ylim,'linestyle','--','color','k');
line(xlim,[0.5,0.5],'linestyle','--','color','k');
xlabel('Condition');
ylabel('Fraction go left');
title('Late move');


%% Fig 2b/c: Cortical/striatal activity during task

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
data_fn = ['trial_activity_choiceworld'];

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);
reward_all = cellfun(@(data,use_expts) data(use_expts),reward_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);
reward_all = reward_all(use_animals);

% Get number of widefield components and MUA depths
n_vs = size(fluor_all{1}{1},3);
n_depths = size(mua_all{1}{1},3);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_vs,1);
mua_allcat = nan(0,length(t),n_depths,1);
wheel_allcat = nan(0,length(t),1);
reward_allcat = nan(0,length(t),1);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

D_cat = struct('stimulus',cell(1,1),'response',cell(1,1),'repeatNum',cell(1,1));

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
        
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths),'uni',false));
    
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x,[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
%     wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
    % Concatenate behavioural data
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
    D.response = cell2mat(cellfun(@(x) x.response,D_all{curr_animal},'uni',false));
    D.repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all{curr_animal},'uni',false));
    
    D.day = trial_day;
    
    D_cat.stimulus = vertcat(D_cat.stimulus,D.stimulus);
    D_cat.response = vertcat(D_cat.response,D.response);
    D_cat.repeatNum = vertcat(D_cat.repeatNum,D.repeatNum);
    
    % Get trial ID   
    trial_contrast = max(D.stimulus,[],2);
    [~,side_idx] = max(D.stimulus > 0,[],2);
    trial_side = (side_idx-1.5)*2;
    trial_choice = -(D.response-1.5)*2;
    
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    reward_allcat = [reward_allcat;vertcat(reward_all{curr_animal}{:})];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];
   
end

% Get max velocity in chosen direction
wheel_velocity_allcat = [zeros(size(wheel_allcat,1),1),diff(wheel_allcat,[],2)];
[max_speed,max_vel_idx] = max(abs(wheel_velocity_allcat(:,:,1).* ...
    (bsxfun(@times,wheel_velocity_allcat(:,:,1),trial_choice_allcat) > 0)),[],2);
vel_t = t(max_vel_idx);
max_vel = max_speed.*trial_choice_allcat;

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Reaction times to plot
rxn_time_bins = {[0.2,0.5]};

move_align = true;
normalize_px = true;

% Get major trial types
vis_L_trials_hit = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == -1;

vis_R_trials_hit = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == -1 & ...
    trial_choice_allcat == 1;

vis_L_trials_miss = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == -1 & ...
    trial_choice_allcat == -1;

vis_R_trials_miss = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == 1;

zero_L_trials = ...
    trial_contrast_allcat == 0 & ...
    trial_choice_allcat == -1;

zero_R_trials = ...
    trial_contrast_allcat == 0 & ...
    trial_choice_allcat == 1;

trial_types = ...
    [vis_L_trials_hit, vis_R_trials_hit, ...
    vis_L_trials_miss, vis_R_trials_miss, ...
    zero_L_trials, zero_R_trials];

%%% Get and plot fluorescence
px_trial_types = nan(size(U_master,1),size(U_master,2), ...
    size(fluor_allcat,2)-1,size(trial_types,2),length(rxn_time_bins));
for curr_rxn = 1:length(rxn_time_bins)
    for curr_trial_type = 1:size(trial_types,2)
        
        curr_trials = trial_types(:,curr_trial_type) & ...
            move_t > rxn_time_bins{curr_rxn}(1) & ...
            move_t < rxn_time_bins{curr_rxn}(2);
        
        % Fluorescence derivative
        curr_data = fluor_allcat(curr_trials,:,:);
        
%         % Smoothed Fluorescence derivative
%         smooth_factor = 3;
%         curr_data = convn(fluor_allcat(curr_trials,:,:), ...
%             ones(1,smooth_factor)/smooth_factor,'same');
        
        if move_align
            % Re-align to movement onset
            t_leeway = 0.5;
            leeway_samples = round(t_leeway*sample_rate);
            curr_move_idx = move_idx(curr_trials);
            for i = 1:size(curr_data,1)
                curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
            end
        end
        
        curr_data_mean = squeeze(nanmean(curr_data,1))';       
        curr_px = diff(svdFrameReconstruct(U_master(:,:,1:n_vs),curr_data_mean),[],3);
%         curr_px(curr_px < 0) = 0;
        
        px_trial_types(:,:,:,curr_trial_type,curr_rxn) = curr_px;
    end
end

t_diff = conv(t,[1,1]/2,'valid');

% Flip and combine trial types
% 1) visual hit (contra), 2) visual miss (ipsi), 3) zero (L)
px_combined = cat(4, ...
    (px_trial_types(:,:,:,1,:) + AP_reflect_widefield(px_trial_types(:,:,:,2,:)))./2, ...
    (px_trial_types(:,:,:,3,:) + AP_reflect_widefield(px_trial_types(:,:,:,4,:)))./2, ...
    (px_trial_types(:,:,:,5,:) + AP_reflect_widefield(px_trial_types(:,:,:,6,:)))./2);

px_combined_hemidiff = px_combined - AP_reflect_widefield(px_combined);

if normalize_px
    % Normalize by dividing by max of each frame
    px_dims = size(px_combined);
    px_combined = bsxfun(@rdivide,px_combined,max(abs(reshape( ...
         px_combined,[px_dims(1)*px_dims(2),1,px_dims(3:end)])),[],1));
     
    px_combined_hemidiff = bsxfun(@rdivide,px_combined_hemidiff,max(abs(reshape( ...
         px_combined_hemidiff,[px_dims(1)*px_dims(2),1,px_dims(3:end)])),[],1));
end

% Plot
plot_rxn = 1;
AP_image_scroll(cat(4,px_combined(:,:,:,:,plot_rxn),px_combined_hemidiff(:,:,:,:,plot_rxn)),t_diff);
axis image; caxis([-1,1]); 
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k');

% Plot all concatenated
px_dims = size(px_combined);
AP_image_scroll([reshape(permute(px_combined,[1,4,2,5,3]), ...
    [],size(px_combined,2)*size(px_combined,5),size(px_combined,3)), ...
    reshape(permute(px_combined_hemidiff,[1,4,2,5,3]), ...
    [],size(px_combined_hemidiff,2)*size(px_combined_hemidiff,5), ...
    size(px_combined_hemidiff,3))],t_diff);
axis image;
AP_reference_outline('ccf_aligned','k',[], ...
    [size(px_combined,1),size(px_combined,2),3,2]);
caxis([-1,1]);
colormap(brewermap([],'*RdBu'))

%%% Get and plot MUA 
mua_trial_types = nan(n_depths, ...
    size(mua_allcat,2),size(trial_types,2),length(rxn_time_bins));
for curr_rxn = 1:length(rxn_time_bins)
    for curr_trial_type = 1:size(trial_types,2)
        
        curr_trials = trial_types(:,curr_trial_type) & ...
            move_t > rxn_time_bins{curr_rxn}(1) & ...
            move_t < rxn_time_bins{curr_rxn}(2);
        
        curr_data = mua_allcat(curr_trials,:,:);
        
        if move_align
            % Re-align to movement onset
            t_leeway = 0.5;
            leeway_samples = round(t_leeway*sample_rate);
            curr_move_idx = move_idx(curr_trials);
            for i = 1:size(curr_data,1)
                curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
            end
        end
        
        mua_trial_types(:,:,curr_trial_type,curr_rxn) = squeeze(nanmean(curr_data,1))';              
    end
end

% mua_trial_types_pseudodiff = cat(3,mua_trial_types(:,:,[1,3,5],:), ...
%     mua_trial_types(:,:,[1,3,5],:) -  mua_trial_types(:,:,[2,4,6],:));
mua_trial_types_pseudodiff = cat(3,mua_trial_types(:,:,[1,3,5],:), ...
    mua_trial_types(:,:,[2,4,6],:));

figure;
for curr_rxn = 1:length(rxn_time_bins)
    for curr_trial_type = 1:size(trial_types,2)/2
        
        subplot(size(trial_types,2)/2,length(rxn_time_bins)*2, ...
            curr_rxn + length(rxn_time_bins)*2*(curr_trial_type-1)); 
        hold on; set(gca,'ColorOrder',copper(n_depths));
        plot(t, ...
            mua_trial_types_pseudodiff(:,:,curr_trial_type,curr_rxn)', ...
            'linewidth',2);
        ylim([min(mua_trial_types(:)),max(mua_trial_types(:))]);
        line([0,0],ylim,'color','k')
        
        subplot(size(trial_types,2)/2,length(rxn_time_bins)*2, ...
            curr_rxn + length(rxn_time_bins) + length(rxn_time_bins)*2*(curr_trial_type-1)); 
        hold on; set(gca,'ColorOrder',copper(n_depths));
        plot(t, ...
            mua_trial_types_pseudodiff(:,:,curr_trial_type+3,curr_rxn)', ...
            'linewidth',2);
       ylim([min(mua_trial_types(:)),max(mua_trial_types(:))]);
        line([0,0],ylim,'color','k')
        
    end
end

%%% Get and plot wheel 
wheel_trial_types = nan(...
    size(wheel_velocity_allcat,2),size(trial_types,2),length(rxn_time_bins));
for curr_rxn = 1:length(rxn_time_bins)
    for curr_trial_type = 1:size(trial_types,2)
        
        curr_trials = trial_types(:,curr_trial_type) & ...
            move_t > rxn_time_bins{curr_rxn}(1) & ...
            move_t < rxn_time_bins{curr_rxn}(2);
        
        curr_data = wheel_velocity_allcat(curr_trials,:);
        if move_align
            % Re-align to movement onset
            t_leeway = 0.5;
            leeway_samples = round(t_leeway*sample_rate);
            curr_move_idx = move_idx(curr_trials);
            for i = 1:size(curr_data,1)
                curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
            end
        end
        wheel_trial_types(:,curr_trial_type,curr_rxn) = squeeze(nanmean(curr_data,1))';
               
    end
end

figure;
for curr_rxn = 1:length(rxn_time_bins)
    subplot(1,length(rxn_time_bins),curr_rxn); hold on;
    plot(t,wheel_trial_types(:,:,curr_rxn),'linewidth',2);
    ylim([-max(abs(wheel_trial_types(:))),max(abs(wheel_trial_types(:)))]);
    line([0,0],ylim,'color','k')
end


%% Fig 3b/c: regression from task events to cortex/striatum
% (clean this up, save everything necessary after regression)

%%% Load and prepare regression

% Load regression results
regression_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\all_trial_activity_regressed.mat';
load(regression_fn);

% Get number of V's/depths
n_vs = size(fluor_allcat_predicted,3);
n_depths = size(mua_allcat_predicted,3);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Get move onset index
[~,move_idx] = max(abs(wheel(:,:,1)) > 2,[],2);
move_t = t_downsample_diff(move_idx)';

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Set which components to use
use_svs = 200;

% Get time shifts in samples
sample_shifts = cellfun(@(x) round(x(1)*(sample_rate/downsample_factor)): ...
    round(x(2)*(sample_rate/downsample_factor)),t_shifts,'uni',false);

%%% Plot kernels

% Plot fluorescence kernels
for curr_regressor = 1:length(regressors)
    curr_k_v = permute(cell2mat(permute(fluor_kernel(curr_regressor,1:use_svs),[1,3,2])),[3,2,1]);
    curr_k_px = reshape(svdFrameReconstruct(U_master(:,:,1:use_svs), ...
        reshape(curr_k_v,use_svs,[])),size(U_master,1),size(U_master,2),[],size(curr_k_v,3));
    AP_image_scroll(curr_k_px,sample_shifts{curr_regressor}/(sample_rate/downsample_factor));
    axis image
    caxis([-prctile(abs(curr_k_px(:)),100),prctile(abs(curr_k_px(:)),100)]);
    colormap(brewermap([],'*RdBu'))
    AP_reference_outline('ccf_aligned','k');
    set(gcf,'Name',regressor_labels{curr_regressor});
end
% Plot fluorescence constant
curr_k_v = permute(cell2mat(permute(fluor_kernel(length(regressors)+1,1:use_svs),[1,3,2])),[3,2,1]);
curr_k_px = reshape(svdFrameReconstruct(U_master(:,:,1:use_svs), ...
    reshape(curr_k_v,use_svs,[])),size(U_master,1),size(U_master,2),[],size(curr_k_v,3));
figure;imagesc(curr_k_px);
axis image off
caxis([-prctile(abs(curr_k_px(:)),100),prctile(abs(curr_k_px(:)),100)]);
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k');
title('Constant');

% Plot fluorescence ROI kernels
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);
n_rois = numel(wf_roi);

roi_mask = cat(3,wf_roi.mask);

U_roi = bsxfun(@rdivide,transpose(reshape(U_master,[],size(U_master,3))'* ...
    reshape(roi_mask,[],size(roi_mask,3))), ...
    sum(reshape(roi_mask,[],size(roi_mask,3)),1)');

for curr_regressor = 1:length(regressors)
    curr_k_v = permute(cell2mat(permute(fluor_kernel(curr_regressor,1:use_svs),[1,3,2])),[3,2,1]);
    
    curr_k_roi = reshape(AP_svd_roi(U_master(:,:,1:n_vs), ...
        reshape(curr_k_v,n_vs,[]),[],[],roi_mask)', ...
        size(curr_k_v,2),size(curr_k_v,3),n_rois);
    
    if size(curr_k_roi,2) > 1
        curr_col = colormap_BlueWhiteRed(size(curr_k_roi,2)/2);
        curr_col(size(curr_k_roi,2)/2+1,:) = [];
    else
        curr_col = 'k'
    end
    
    figure; hold on;
    for curr_subk = 1:size(curr_k_roi,2)
        AP_stackplot(squeeze(curr_k_roi(:,curr_subk,:)), ...
            sample_shifts{curr_regressor}/(sample_rate/downsample_factor), ...
            1.5e-3,false,curr_col(curr_subk,:));
    end
    title(regressor_labels{curr_regressor});
end

% Plot MUA kernels
figure;
for curr_regressor = 1:length(regressors)
    
    curr_k_cat = permute(cat(3,mua_kernel{curr_regressor,:}),[2,1,3]);
    
    if size(curr_k_cat,2) > 1
        curr_col = colormap_BlueWhiteRed(size(curr_k_cat,2)/2);
        curr_col(size(curr_k_cat,2)/2+1,:) = [];
    else
        curr_col = 'k';
    end
    
    figure; hold on;
    for curr_subk = 1:size(curr_k_cat,2)
        AP_stackplot(squeeze(curr_k_cat(:,curr_subk,:)), ...
            sample_shifts{curr_regressor}/(sample_rate/downsample_factor), ...
            1,false,curr_col(curr_subk,:));
    end
    title(regressor_labels{curr_regressor});
end


%% Fig 4a: Passive stim responses in cortex/striatum
% (TO DO: normalize this the same as during task so they're directly
% comparable)
% (Do this for both types of passive stim)

% saving new trial data here in progress: 
% C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
data_fn = ['trial_activity_passive_fullscreen_naive'];

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% % Load use experiments and cut out bad ones
% exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
% exclude_fn{1} = 'bhv_use_experiments';
% % exclude_fn{2} = 'expl_var_use_experiments';
% use_experiments_all = {};
% for curr_exclude = 1:length(exclude_fn)
%     curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
%     use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
% end
% use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
%     1:size(use_experiments_all,2),'uni',false);
% 
% D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
% fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
% mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
% wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);

% Get widefield ROIs
n_vs = 200;

% MUA depths
n_depths = size(mua_all{1}{1},3);

% Set up conditions
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_vs,1);
mua_allcat = nan(0,length(t),n_depths,1);
wheel_allcat = nan(0,length(t),1);
reward_allcat = nan(0,length(t),1);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

D_cat = struct('stimulus',cell(1,1),'response',cell(1,1),'repeatNum',cell(1,1));

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    n_align = size(fluor_cat,4);
    
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    t_std = t > -0.5 & t < 0.5;
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x(:,t_std,:,1),[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
    
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
    
    % Concatenate stim
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
    D_cat.stimulus = vertcat(D_cat.stimulus,D.stimulus);
   
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    
end

%%% REMOVE MOVE TRIALS FOR PASSIVE
move_trials = any(abs(wheel_allcat(:,t >= -0.5 & t <= 2) > 2),2);
fluor_allcat(move_trials,:,:) = [];
mua_allcat(move_trials,:,:) = [];
wheel_allcat(move_trials,:,:) = [];
D_cat.stimulus(move_trials) = [];
%%%

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

%%% Get and plot fluorescence
px_stim = nan(size(U_master,1),size(U_master,2), ...
    size(fluor_allcat,2)-1,length(unique(D_cat.stimulus)));
for curr_stim = unique(D_cat.stimulus)'
    
    curr_trials = D_cat.stimulus == curr_stim;
    
    %             % Straight fluorescence
    %             curr_data = fluor_allcat(curr_trials,:,:,plot_align);
    %             curr_baseline = nanmean(fluor_allcat(curr_trials,t > -0.2 & t < 0,:,1),2);
    %             curr_data = squeeze(nanmean(bsxfun(@minus,curr_data,curr_baseline),1))';
    %             curr_px = svdFrameReconstruct(U_master(:,:,1:n_rois),curr_data(:,1:end-1));
    
    %         % Fluorescence derivative
    %         curr_data = squeeze(nanmean(fluor_allcat(curr_trials,:,:),1))';
    %         curr_px = diff(svdFrameReconstruct(U_master(:,:,1:n_rois),curr_data),[],3);
    %         curr_px(curr_px < 0) = 0;
    
    % Smoothed Fluorescence derivative
    smooth_factor = 3;
    curr_data = convn(fluor_allcat(curr_trials,:,:), ...
        ones(1,smooth_factor)/smooth_factor,'same');
    
    curr_data_mean = squeeze(nanmean(curr_data,1))';   
    curr_px = diff(svdFrameReconstruct(U_master(:,:,1:n_vs),curr_data_mean),[],3);
%     curr_px(curr_px < 0) = 0;
    
    px_stim(:,:,:,curr_stim) = curr_px;
end

t_diff = conv(t,[1,1]/2,'valid');

% Choose stim to plot
plot_stim = unique(D_cat.stimulus);
% plot_stim = [2,5,8,10];

% Plot all together
AP_image_scroll(reshape(permute(px_stim(:,:,:,plot_stim),[1,2,4,3]), ...
    size(px_stim,1),[],size(px_stim,3)),t_diff);
axis image; caxis([-1.5e-3,1.5e-3]); 
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k',[], ...
    [size(px_stim,1),size(px_stim,2),1,length(plot_stim)]);

% Plot striatum
plot_cols = lines(length(plot_stim));
figure; hold on;
p = gobjects(n_depths,length(plot_stim));
for curr_plot_condition = 1:length(plot_stim)  
    curr_trials = D_cat.stimulus == plot_stim(curr_plot_condition);
    curr_data_mean = squeeze(nanmean(mua_allcat(curr_trials,:,:),1));   
    curr_col = plot_cols(curr_plot_condition,:);   
    p(:,curr_plot_condition) = AP_stackplot(curr_data_mean,t,1,false,curr_col,1:n_depths);    
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');
legend(p(1,:),cellfun(@num2str,num2cell(unique(D_cat.stimulus)),'uni',false));

% Get fluorescence in widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);
n_rois = numel(wf_roi);

roi_mask = cat(3,wf_roi.mask);

U_roi = bsxfun(@rdivide,transpose(reshape(U_master,[],size(U_master,3))'* ...
    reshape(roi_mask,[],size(roi_mask,3))), ...
    sum(reshape(roi_mask,[],size(roi_mask,3)),1)');

smooth_size = 1;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

fluor_roi = permute(reshape(transpose(U_roi(:,1:n_vs)* ...
    reshape(permute(fluor_allcat(:,:,:,1),[3,2,1,4]),n_vs,[])), ...
    size(fluor_allcat,2),size(fluor_allcat,1),n_rois),[2,1,3]);

fluor_roi_diff = diff(padarray(convn(fluor_roi,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both'),[],2);

% (zero negatives, normalize)
fluor_roi_diff(fluor_roi_diff < 0) = 0;
fluor_roi_diff = bsxfun(@rdivide,fluor_roi_diff,nanstd(reshape(fluor_roi_diff,[],1,size(wf_roi,1)),[],1));

t_diff =  conv(t,[1,1]/2,'valid');

figure; hold on;
plot_cols = lines(length(plot_stim));
for curr_plot_condition = 1:length(plot_stim)
    
    curr_trials = D_cat.stimulus == plot_stim(curr_plot_condition);
    curr_data = fluor_roi_diff(curr_trials,:,:);
   
    curr_data_mean = squeeze(nanmean(curr_data,1));

    curr_col = plot_cols(curr_plot_condition,:);   
    AP_stackplot(curr_data_mean,t_diff,3.5,false,curr_col,{wf_roi.area});
    
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');
title('Widefield ROI');



% Plot tuning curves?
condition = reshape([0.06,0.125,0.25,0.5,1]'.*[-1,1],[],1);
[~,condition_sort_idx] = sort(condition);

t_ctx = [0.05,0.15];
ctx_mean = squeeze(nanmean(fluor_roi_diff(:,t >= t_ctx(1) & t <= t_ctx(2),:),2));
ctx_condition_mean = grpstats(ctx_mean,D_cat.stimulus,'nanmean');

t_str = [0.05,0.15];
str_mean = squeeze(nanmean(mua_allcat(:,t >= t_str(1) & t <= t_str(2),:),2));
str_condition_mean = grpstats(str_mean,D_cat.stimulus,'nanmean');

figure; 
subplot(1,2,1); hold on;
set(gca,'ColorOrder',jet(n_rois));
plot(condition(condition_sort_idx),ctx_condition_mean(condition_sort_idx,:),'linewidth',2)
ylabel('Fluorescence (std)');
xlabel('Contrast*side');

subplot(1,2,2); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot(condition(condition_sort_idx),str_condition_mean(condition_sort_idx,:),'linewidth',2)
ylabel('Spikes (std)');
xlabel('Contrast*side');





