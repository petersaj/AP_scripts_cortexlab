% Generate figures for ctx-str paper

% Anything that takes a lot of time is done in
% AP_ctx_str_trial_preprocessing and saved for plotting here

% The original scripts here were in test_wf_ephys_choiceworld_analysis

%% Fig 1a: Behavior psychometric 
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

%% Fig 1 b/c: Cortical activity during task

data_fn = ['trial_activity_choiceworld'];
exclude_data = true;

[fluor_allcat,fluor_roi_diff,mua_allcat,wheel_allcat,reward_allcat,D_allcat] = ...
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
t_diff =  conv(t,[1,1]/2,'valid');

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

%% Fig 1 d/e: Task -> cortical activity regression

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

% Plot cortex kernels
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

% Plot cortex constant
curr_k_v = permute(cell2mat(permute(fluor_kernel(length(regressors)+1,1:use_svs),[1,3,2])),[3,2,1]);
curr_k_px = reshape(svdFrameReconstruct(U_master(:,:,1:use_svs), ...
    reshape(curr_k_v,use_svs,[])),size(U_master,1),size(U_master,2),[],size(curr_k_v,3));
figure;imagesc(curr_k_px);
axis image off
caxis([-prctile(abs(curr_k_px(:)),100),prctile(abs(curr_k_px(:)),100)]);
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k');
title('Constant');

% Plot cortex ROI kernels
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

% Plot cortex ROI actual and predicted
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

fluor_downsample_diff_roi = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_downsamp_diff,[3,2,1]),n_vs,[]),[],[],roi_mask), ...
    size(roi_mask,3),[],size(fluor_allcat_downsamp_diff,1)),[3,2,1]);

fluor_predicted_roi = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_predicted,[3,2,1]),n_vs,[]),[],[],roi_mask), ...
    size(roi_mask,3),[],size(fluor_allcat_predicted,1)),[3,2,1]);

figure;
p = [];
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
    
    p(curr_plot) = subplot(1,3,curr_plot); hold on;
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
linkaxes(p);

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

%% Fig 2b: Example traces/kernels

warning('Probably not best example, check others');

% Load and align
str_align = 'kernel';

% animal = 'AP028'; 
% day = '2017-12-20'; % 16/20 (blood in 16, a little in 20)

animal = 'AP027';
day = '2017-11-25';

experiment = 1; 
verbose = false; 
AP_load_experiment;

avg_im_aligned = AP_align_widefield(animal,day,avg_im);
Udf_aligned = single(AP_align_widefield(animal,day,Udf));

% Smoothing window
smooth_size = 9; % MUST BE ODD
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

% Define ROIs and get fluorescence traces
roi_circle_size = 10;
roi_x = [131,160,128,51]; % roi_x = [131,174,110,51];
roi_y = [297,131,59,144]; % roi_y = [297,96,71,144];
[x,y] = meshgrid(1:size(avg_im_aligned,1),1:size(avg_im_aligned,2));
roi_mask = cell2mat(arrayfun(@(roi) sqrt((x-roi_x(roi)).^2 + (y-roi_y(roi)).^2) <= ...
    roi_circle_size,permute(1:length(roi_x),[1,3,2]),'uni',false));
roi_trace = AP_svd_roi(Udf_aligned,fVdf,[],[],roi_mask);

roi_trace_deriv = diff(convn(roi_trace,smWin,'same'),[],2);
roi_trace_deriv(roi_trace_deriv < 0) = 0;
frame_t_deriv = conv(frame_t,[1,1]/2,'valid');

% % Bin spikes by aligned depth
% n_depths = n_aligned_depths;
% depth_group = aligned_str_depth_group;
% 
% time_bins = [frame_t_deriv,frame_t_deriv(end)+1/framerate];
% binned_spikes = zeros(n_depths,length(time_bins)-1);
% for curr_depth = 1:n_depths
%     curr_spike_times = spike_times_timeline(depth_group == curr_depth);
%     binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
% end

% Bin spikes evenly across striatum

n_depths = 4;
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
depth_group = discretize(spikeDepths,depth_group_edges);
use_depths = 1:n_depths;

% depth_group_edges = [str_depth(1),str_depth(1)+500, ...
%     str_depth(2)-500,str_depth(2)];
% depth_group = discretize(spikeDepths,depth_group_edges);
% use_depths = [1,max(depth_group)];

% n_depths = round(diff(str_depth)/200);
% depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
% depth_group = discretize(spikeDepths,depth_group_edges);
% use_depths = 1:n_depths;

time_bins = [frame_t_deriv,frame_t_deriv(end)+1/framerate];
binned_spikes = zeros(length(use_depths),length(time_bins)-1);
for curr_depth = 1:length(use_depths)
    curr_spike_times = spike_times_timeline(depth_group == use_depths(curr_depth));
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
end
binned_spikes_smoothed = convn(binned_spikes,smWin,'same');

% Get spike-triggered average
skip_seconds = 60;
surround_times = [-0.2,0.2];

framerate = 1./median(diff(frame_t));
skip_frames = round(skip_seconds*framerate);
surround_frames = round(surround_times(1)*framerate):round(surround_times(2)*framerate);
sta_t = surround_frames./framerate;

sta_im = zeros(size(Udf_aligned,1),size(Udf_aligned,2), ...
    length(surround_frames),size(binned_spikes,1));

for curr_depth = 1:size(binned_spikes,1)
    frames_w = repmat(binned_spikes(curr_depth,skip_frames:end-skip_frames)'./ ...
        sum(binned_spikes(curr_depth,skip_frames:end-skip_frames)),1,length(surround_frames));
    for curr_sta_frame = 1:length(surround_frames)
        frames_w(:,curr_sta_frame) = ...
            circshift(frames_w(:,curr_sta_frame),[surround_frames(curr_sta_frame),0]);
    end
    
    sta_v = diff(fVdf(:,skip_frames:end-skip_frames),[],2)*frames_w;
    sta_im(:,:,:,curr_depth) = svdFrameReconstruct(Udf_aligned,sta_v);
end

sta_im_max = squeeze(max(sta_im,[],3));

% Regress fluorescence to spikes
use_svs = 1:50;
skip_seconds = 60;
upsample_factor = 1;
kernel_t = [-0.3,0.3];

lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
load(lambda_fn);
curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
if any(curr_animal_idx)
    curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
    if any(curr_day_idx)
        lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
    end
end

large_lambda = lambda*10;

kernel_frames = round(kernel_t(1)*framerate): ...
    round(kernel_t(2)*framerate);
zs = [false,true];
cvfold = 10;

[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(diff(fVdf(use_svs,:),[],2), ...
    binned_spikes,kernel_frames,large_lambda,zs,cvfold);

Udf_aligned = single(AP_align_widefield(animal,day,Udf));
k_px = zeros(size(Udf_aligned,1),size(Udf_aligned,2),size(k,2),size(k,3),'single');
for curr_spikes = 1:size(k,3)
    k_px(:,:,:,curr_spikes) = ...
        svdFrameReconstruct(Udf_aligned(:,:,use_svs),k(:,:,curr_spikes));
end

k_px_max = squeeze(max(k_px,[],3));

% Plot STA and regression
figure;
for curr_depth = 1:size(binned_spikes,1)
    h = subplot(2,size(binned_spikes,1),curr_depth);
    imagesc(sta_im_max(:,:,curr_depth));
    caxis([0,max(abs(caxis))]);
    colormap(h,brewermap([],'BuGn'));
    AP_reference_outline('ccf_aligned','k');
    axis image off;
    
    h = subplot(2,size(binned_spikes,1),size(binned_spikes,1) + curr_depth);
    imagesc(k_px_max(:,:,curr_depth));
    caxis([0,max(abs(caxis))]);
    colormap(h,brewermap([],'BuGn'));
    AP_reference_outline('ccf_aligned','k');
    axis image off;
end

% Plot ROIs and traces
figure;

subplot(1,3,1);
roi_boundaries = bwboundaries(sum(roi_mask,3));
imagesc(avg_im_aligned);colormap(gray);
caxis([0,prctile(avg_im_aligned(:),99)]);
axis image off;
AP_reference_outline('ccf_aligned','r');
p = cellfun(@(x) plot(x(:,2),x(:,1),'b','linewidth',2),roi_boundaries);

p1 = subplot(2,3,2:3);
AP_stackplot(bsxfun(@rdivide,roi_trace_deriv,std(roi_trace_deriv,[],2))', ...
    frame_t_deriv,5,false,[0,0.7,0]);
title('\DeltaFluorescence')
xlabel('Time (seconds)');
ylabel('Activity (std)');

p2 = subplot(2,3,5:6);
AP_stackplot(bsxfun(@rdivide,binned_spikes_smoothed,std(binned_spikes,[],2))', ...
    frame_t_deriv,5,false,'k');
title('MUA')
xlabel('Time (seconds)');
ylabel('Activity (std)');

linkaxes([p1,p2],'xy');
xlim([80,107]);


%% Fig 2c,e: Average regression maps (concatenated protocols)

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

% Plot kernels at t = 0;
figure;
colormap(brewermap([],'*RdBu'));
for curr_depth = 1:n_aligned_depths
    subplot(2,n_aligned_depths,curr_depth)
    imagesc(k_px_trained(:,:,t == 0,curr_depth));
    caxis([-max(caxis),max(caxis)]*0.8)
    axis image off;
    AP_reference_outline('ccf_aligned','k');
    
    subplot(2,n_aligned_depths,n_aligned_depths + curr_depth)
    imagesc(k_px_naive(:,:,t == 0,curr_depth));
    caxis([-max(caxis),max(caxis)]*0.8)
    axis image off;
    AP_reference_outline('ccf_aligned','k');
end

% Get center-of-mass maps
jet_basic = jet(255);
dark_colors = max(jet_basic,[],2) ~= 1;
jet_alt = interp1(1:255,jet_basic,linspace(find(~dark_colors,1,'first'), ...
    find(~dark_colors,1,'last'),255)) - 0.2;
use_colormap = jet_alt;

k_px_trained_positive = k_px_trained;
k_px_trained_positive(k_px_trained_positive < 0) = 0;
k_px_trained_com = sum(k_px_trained_positive.*permute(1:n_aligned_depths,[1,3,4,2]),4)./sum(k_px_trained_positive,4);
k_px_trained_com_colored = nan(size(k_px_trained_com,1),size(k_px_trained_com,2),3,size(k_px_trained_com,3));
for curr_frame = 1:size(k_px_trained_com,3)
    k_px_trained_com_colored(:,:,:,curr_frame) = ...
        ind2rgb(round(mat2gray(k_px_trained_com(:,:,curr_frame),[1,n_aligned_depths])*255),use_colormap);
end

k_px_naive_positive = k_px_naive;
k_px_naive_positive(k_px_naive_positive < 0) = 0;
k_px_naive_com = sum(k_px_naive_positive.*permute(1:n_aligned_depths,[1,3,4,2]),4)./sum(k_px_naive_positive,4);
k_px_naive_com_colored = nan(size(k_px_naive_com,1),size(k_px_naive_com,2),3,size(k_px_naive_com,3));
for curr_frame = 1:size(k_px_naive_com,3)
    k_px_naive_com_colored(:,:,:,curr_frame) = ...
        ind2rgb(round(mat2gray(k_px_naive_com(:,:,curr_frame),[1,n_aligned_depths])*255),use_colormap);
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

% Plot at all multiple time points, transparency relative to weight
plot_t = find(t >= -0.05 & t <= 0.05);
weight_max = max(k_px_trained(:))*0.8;
figure;
for t_idx = 1:length(plot_t)
    curr_t = plot_t(t_idx);  
    subplot(1,length(plot_t),t_idx);
    p = image(k_px_trained_com_colored(:,:,:,curr_t));
    set(p,'AlphaData', ...
        mat2gray(max(k_px_trained(:,:,curr_t,:),[],4),[0,weight_max]));
    axis image off;
    AP_reference_outline('ccf_aligned','k');
    title(t(curr_t));
end

%% Fig 2f: Average regression maps (protocols separately)

protocols = {'vanillaChoiceworld', ...
    'stimSparseNoiseUncorrAsync', ...
    'stimKalatsky'};

for protocol = protocols
    
    curr_protocol = cell2mat(protocol);
    
    n_aligned_depths = 4;
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_ephys';
    k_fn = [data_path filesep 'wf_ephys_maps_' curr_protocol '_' num2str(n_aligned_depths) '_depths_kernel'];
    load(k_fn);
    
    framerate = 35;
    upsample_factor = 2;
    sample_rate = framerate*upsample_factor;
    kernel_t = [-0.3,0.3];
    kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
    t = kernel_frames/sample_rate;
    
    % Concatenate and mean
    % (kernel goes backwards in time - flip to correct)
    k_px_cat = cellfun(@(x) x(:,:,end:-1:1,:),[batch_vars.r_px],'uni',false);
    k_px = nanmean(double(cat(5,k_px_cat{:})),5);
    
    % Get center-of-mass maps
    k_px_positive = k_px;
    k_px_positive(k_px_positive < 0) = 0;
    k_px_com = sum(k_px_positive.*permute(1:n_aligned_depths,[1,3,4,2]),4)./sum(k_px_positive,4);
    k_px_com_colored = nan(size(k_px_com,1),size(k_px_com,2),3,size(k_px_com,3));
    
    jet_basic = jet(255);
    dark_colors = max(jet_basic,[],2) ~= 1;
    jet_alt = interp1(1:255,jet_basic,linspace(find(~dark_colors,1,'first'), ...
        find(~dark_colors,1,'last'),255)) - 0.2;
    use_colormap = jet_alt;
    
    for curr_frame = 1:size(k_px_com,3)
        k_px_com_colored(:,:,:,curr_frame) = ...
            ind2rgb(round(mat2gray(k_px_com(:,:,curr_frame),[1,n_aligned_depths])*255),use_colormap);
    end
    
    % Plot at t = 0, transparency relative to weight
    figure;
    p = image(k_px_com_colored(:,:,:,t == 0));
    % weight_max = max(k_px(:))*0.8;
    weight_max = 0.01;
    set(p,'AlphaData', ...
        mat2gray(max(k_px(:,:,t == 0,:),[],4),[0,weight_max]));
    axis image off;
    AP_reference_outline('ccf_aligned','k');
    title([curr_protocol ': t = 0']);
    drawnow;
    
end


%% Fig 2c: Allen projection maps
% (code from AP_ctx2str_probe)

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
    
%     imagesc(kernel_template(:,:,curr_depth));
%     caxis([-prctile(abs(kernel_template(:)),99),prctile(abs(kernel_template(:)),99)]);
    
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

%% Fig 3?: Striatal activity during task

data_fn = ['trial_activity_choiceworld'];
exclude_data = true;

[fluor_allcat,fluor_roi_diff,mua_allcat,wheel_allcat,reward_allcat,D_allcat] = ...
    AP_load_concat_normalize_ctx_str(data_fn,exclude_data);

n_depths = size(fluor_allcat,3);

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

% Get trial types to plot
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

% Plot MUA
figure;
for curr_align = 1:2
    p(curr_align) = subplot(1,2,curr_align); hold on;
    for curr_plot_condition = 1:size(plot_conditions,1)
        
        curr_trials = plot_id == curr_plot_condition & ...
            move_t >= rxn_time_use(1) & ...
            move_t <= rxn_time_use(2);

        curr_data = mua_allcat(curr_trials,:,:);
        
        if curr_align == 2
            % Re-align to movement onset
            t_leeway = -t(1);
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
        
        AP_stackplot(curr_data_mean,t,3.5,false,curr_col);
        
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

%% Fig 3a: Task -> striatal activity regression

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

% Get time shifts in samples
sample_shifts = cellfun(@(x) round(x(1)*(sample_rate/downsample_factor)): ...
    round(x(2)*(sample_rate/downsample_factor)),t_shifts,'uni',false);

% Plot striatum kernels
figure;
p = [];
for curr_regressor = 1:length(regressors)
    p(curr_regressor) = subplot(1,length(regressors),curr_regressor); hold on;
    
    curr_k_v = permute(cell2mat(permute(mua_kernel(curr_regressor,:),[1,3,2])),[2,1,3]);
    
    if size(curr_k_v,2) > 1
        curr_col = colormap_BlueWhiteRed(size(curr_k_v,2)/2);
        curr_col(size(curr_k_v,2)/2+1,:) = [];
    else
        curr_col = 'k';
    end
    
    for curr_subk = 1:size(curr_k_v,2)
        AP_stackplot(squeeze(curr_k_v(:,curr_subk,:)), ...
            sample_shifts{curr_regressor}/(sample_rate/downsample_factor), ...
            1,false,curr_col(curr_subk,:),1:4,true);
    end
    line([0,0],ylim,'color','k');
    title(regressor_labels{curr_regressor});
end
linkaxes(p,'xy')

% Plot actual and predicted
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

figure;
p = [];
for curr_plot = 1:3
    
    switch curr_plot
        case 1
            plot_data = mua_allcat_downsamp;
            plot_title = 'Measured';
        case 2
            plot_data = mua_allcat_predicted;
            plot_title = 'Predicted';
        case 3
            plot_data = mua_allcat_downsamp - mua_allcat_predicted;
            plot_title = 'Residual';
    end
    
    p(curr_plot) = subplot(1,3,curr_plot); hold on;
    for curr_plot_condition = 1:size(plot_conditions,1)
        
        curr_trials = plot_id == curr_plot_condition & ...
            move_t >= rxn_time_use(1) & ...
            move_t <= rxn_time_use(2);
        curr_data = plot_data(curr_trials,:,:);
        
%         % Re-align to movement onset
%         t_leeway = -t_downsample_diff(1);
%         leeway_samples = round(t_leeway*(sample_rate/downsample_factor));
%         curr_move_idx = move_idx(curr_trials);
%         for i = 1:size(curr_data,1)
%             curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
%         end
        
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
        
        AP_stackplot(curr_data_mean,t_downsample_diff,2,false,curr_col);
        
    end
    axis tight
    line([0,0],ylim,'color','k');
    xlabel('Time from stim');
    title(plot_title);   
end
linkaxes(p)

% Get max stim kernel and contrast slope maps
stim_regressor = strcmp(regressor_labels,'Stim');
k_v_stim = permute(cell2mat(permute(mua_kernel(stim_regressor,:),[1,3,2])),[3,2,1]);
k_v_stim_maxt = squeeze(max(k_v_stim,[],2));

contrasts = [0.06,0.125,0.25,0.5,1];
contrast_sides = [fliplr(contrasts)*-1,contrasts];

contrast_k_L = AP_regresskernel(contrasts, ...
    reshape(k_v_stim_maxt(:,5:-1:1),[],5),0,[],[false,false],1,true,true);
contrast_slope_L = reshape(contrast_k_L{1},[],1);

contrast_k_R = AP_regresskernel(contrasts, ...
    reshape(k_v_stim_maxt(:,6:10),[],5),0,[],[false,false],1,true,true);
contrast_slope_R = reshape(contrast_k_R{1},[],1);

contrast_slope_diff = contrast_slope_R - contrast_slope_L;

% Get max move kernel and move ipsi-contra differences
move_onset_regressor = strcmp(regressor_labels,'Move onset');
k_v_move_onset = permute(cell2mat(permute(mua_kernel(move_onset_regressor,:),[1,3,2])),[3,2,1]);
k_v_move_onset_maxt = squeeze(max(k_v_move_onset,[],2));

k_px_move_onset_maxt_diff = k_v_move_onset_maxt(:,1)-k_v_move_onset_maxt(:,2);

% Get max reward kernel
reward_regressor = strcmp(regressor_labels,'Reward');
k_v_reward = permute(cell2mat(permute(mua_kernel(reward_regressor,:),[1,3,2])),[3,2,1]);
k_v_reward_maxt = squeeze(max(k_v_reward,[],2));

% Plot max kernels and differences
figure;
subplot(1,3,1); hold on;
plot(max(k_v_stim_maxt,[],2),'k','linewidth',2);
plot(max(k_v_move_onset_maxt,[],2),'color',[0.6,0,0.6],'linewidth',2);
plot(max(k_v_reward_maxt,[],2),'color',[0,0,0.6],'linewidth',2);
xlabel('Striatum depth');
ylabel('Maximum kernel weight');
legend({'Stim','Move onset','Reward'});
axis square;

subplot(1,3,2); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot(contrast_sides,k_v_stim_maxt','linewidth',2);
xlabel('Contrast*Side');
ylabel('Max stim kernel');
line(xlim,[0,0],'color','k','linestyle','--');
axis square;

subplot(1,3,3); hold on;
set(gca,'ColorOrder',copper(n_depths));
plot(k_px_move_onset_maxt_diff,'linewidth',2);
xlabel('Striatum depth');
ylabel('Move L - R');
line(xlim,[0,0],'color','k','linestyle','--');
axis square;

%% Fig 4?: stim/no-stim move?



%% Fig 4?: stim responses: task > 0.5s reaction

data_fn = ['trial_activity_choiceworld'];
exclude_data = true;

[fluor_allcat,fluor_roi_diff,mua_allcat,wheel_allcat,reward_allcat,D_allcat] = ...
    AP_load_concat_normalize_ctx_str(data_fn,exclude_data);

n_vs = size(fluor_allcat,3);
n_rois = size(fluor_roi_diff,3);
n_depths = size(mua_allcat,3);

% Get ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Get trial information
trial_contrast_allcat = max(D_allcat.stimulus,[],2);
[~,side_idx] = max(D_allcat.stimulus > 0,[],2);
trial_side_allcat = (side_idx-1.5)*2;
trial_choice_allcat = -(D_allcat.response-1.5)*2;

trial_contrast_side = trial_contrast_allcat.*trial_side_allcat;

% Get time (make this be saved in trial data)
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);
t_diff = conv(t,[1,1]/2,'valid');

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Stim and time to plot
unique_stim = unique(trial_contrast_allcat.*trial_side_allcat);
use_t = t > 0 & t < 0.5;
use_t_diff = t_diff > 0 & t_diff < 0.5;

% Get average long reaction time stim responses
fluor_stim_mean = nan(length(unique_stim),size(fluor_allcat,2),n_vs);
fluor_roi_stim_mean = nan(length(unique_stim),size(fluor_roi_diff,2),n_rois);
mua_stim_mean = nan(length(unique_stim),size(mua_allcat,2),n_depths);
for curr_stim_idx = 1:length(unique_stim)
    curr_trials = move_t > 0.5 & ...
        (trial_contrast_allcat.*trial_side_allcat) == unique_stim(curr_stim_idx);
    fluor_stim_mean(curr_stim_idx,:,:) = nanmean(fluor_allcat(curr_trials,:,:),1);
    fluor_roi_stim_mean(curr_stim_idx,:,:) = nanmean(fluor_roi_diff(curr_trials,:,:),1);
    mua_stim_mean(curr_stim_idx,:,:) = nanmean(mua_allcat(curr_trials,:,:),1);
end

figure('Name','Trained task > 0.5 reaction time');

% Plot 100/R fluorescence
subplot(2,3,1);
plot_px = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    permute(fluor_stim_mean(end,:,:),[3,2,1]));
plot_px_max = max(diff(plot_px(:,:,use_t),[],3),[],3);
imagesc(plot_px_max);
caxis([0,max(abs(caxis))]);
colormap(brewermap([],'BuGn'));
AP_reference_outline('ccf_aligned','k');
axis image off;

% Plot 100/R and contrast curve for visual cortex and striatal depths
subplot(2,3,2); hold on;
plot_rois = [1:3,7];
plot(t_diff,squeeze(fluor_roi_stim_mean(end,:,plot_rois)),'linewidth',2);
axis tight;
line([0,0],ylim,'color','k');
xlim([-0.2,0.5]);
legend({wf_roi(plot_rois).area});
ylabel('Fluorescence');

subplot(2,3,3); hold on; set(gca,'ColorOrder',copper(n_depths));
plot(t,squeeze(mua_stim_mean(end,:,:)),'linewidth',2);
axis tight;
line([0,0],ylim,'color','k');
xlim([-0.2,0.5]);
legend(cellfun(@num2str,num2cell(1:n_depths),'uni',false));
ylabel('MUA');

subplot(2,3,5); hold on;
fluor_roi_stim_max = squeeze(max(fluor_roi_stim_mean(:,use_t_diff,:),[],2));
plot(unique_stim,fluor_roi_stim_max(:,plot_rois),'linewidth',2);
xlabel('Contrast*Side');
ylabel('Fluorescence');

subplot(2,3,6); hold on; set(gca,'ColorOrder',copper(n_depths));
mua_stim_max = squeeze(max(mua_stim_mean(:,use_t,:),[],2));
plot(unique_stim,mua_stim_max,'linewidth',2);
xlabel('Contrast*Side');
ylabel('MUA');

%% Fig 4?: stim responses: passive choiceworld trained

data_fn = ['trial_activity_passive_choiceworld'];
exclude_data = false;

[fluor_allcat,fluor_roi_diff,mua_allcat,wheel_allcat,reward_allcat,D_allcat] = ...
    AP_load_concat_normalize_ctx_str(data_fn,exclude_data);

n_vs = size(fluor_allcat,3);
n_rois = size(fluor_roi_diff,3);
n_depths = size(mua_allcat,3);

% Get ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Get time (make this be saved in trial data)
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);
t_diff = conv(t,[1,1]/2,'valid');

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Stim and time to plot
unique_stim = unique(D_allcat.stimulus);
use_t = t > 0 & t < 0.5;
use_t_diff = t_diff > 0 & t_diff < 0.5;

% Get average long reaction time stim responses
fluor_stim_mean = nan(length(unique_stim),size(fluor_allcat,2),n_vs);
fluor_roi_stim_mean = nan(length(unique_stim),size(fluor_roi_diff,2),n_rois);
mua_stim_mean = nan(length(unique_stim),size(mua_allcat,2),n_depths);
for curr_stim_idx = 1:length(unique_stim)
    curr_trials = move_t > 0.5 & ...
        D_allcat.stimulus == unique_stim(curr_stim_idx);
    fluor_stim_mean(curr_stim_idx,:,:) = nanmean(fluor_allcat(curr_trials,:,:),1);
    fluor_roi_stim_mean(curr_stim_idx,:,:) = nanmean(fluor_roi_diff(curr_trials,:,:),1);
    mua_stim_mean(curr_stim_idx,:,:) = nanmean(mua_allcat(curr_trials,:,:),1);
end

figure('Name','Trained passive choiceworld');

% Plot 100/R fluorescence
subplot(2,3,1);
plot_px = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    permute(fluor_stim_mean(end,:,:),[3,2,1]));
plot_px_max = max(diff(plot_px(:,:,use_t),[],3),[],3);
imagesc(plot_px_max);
caxis([0,max(abs(caxis))]);
colormap(brewermap([],'BuGn'));
AP_reference_outline('ccf_aligned','k');
axis image off;

% Plot 100/R and contrast curve for visual cortex and striatal depths
subplot(2,3,2); hold on;
plot_rois = [1:3,7];
plot(t_diff,squeeze(fluor_roi_stim_mean(end,:,plot_rois)),'linewidth',2);
axis tight;
line([0,0],ylim,'color','k');
xlim([-0.2,0.5]);
legend({wf_roi(plot_rois).area});
ylabel('Fluorescence');

subplot(2,3,3); hold on; set(gca,'ColorOrder',copper(n_depths));
plot(t,squeeze(mua_stim_mean(end,:,:)),'linewidth',2);
axis tight;
line([0,0],ylim,'color','k');
xlim([-0.2,0.5]);
legend(cellfun(@num2str,num2cell(1:n_depths),'uni',false));
ylabel('MUA');

subplot(2,3,5); hold on;
fluor_roi_stim_max = squeeze(max(fluor_roi_stim_mean(:,use_t_diff,:),[],2));
plot(unique_stim,fluor_roi_stim_max(:,plot_rois),'linewidth',2);
xlabel('Contrast*Side');
ylabel('Fluorescence');

subplot(2,3,6); hold on; set(gca,'ColorOrder',copper(n_depths));
mua_stim_max = squeeze(max(mua_stim_mean(:,use_t,:),[],2));
plot(unique_stim,mua_stim_max,'linewidth',2);
xlabel('Contrast*Side');
ylabel('MUA');


%% Fig 4?: stim responses: passive choiceworld naive

data_fn = ['trial_activity_passive_choiceworld_naive'];
exclude_data = false;

[fluor_allcat,fluor_roi_diff,mua_allcat,wheel_allcat,reward_allcat,D_allcat] = ...
    AP_load_concat_normalize_ctx_str(data_fn,exclude_data);

n_vs = size(fluor_allcat,3);
n_rois = size(fluor_roi_diff,3);
n_depths = size(mua_allcat,3);

% Get ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Get time (make this be saved in trial data)
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);
t_diff = conv(t,[1,1]/2,'valid');

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Stim and time to plot
unique_stim = unique(D_allcat.stimulus);
use_t = t > 0 & t < 0.5;
use_t_diff = t_diff > 0 & t_diff < 0.5;

% Get average long reaction time stim responses
fluor_stim_mean = nan(length(unique_stim),size(fluor_allcat,2),n_vs);
fluor_roi_stim_mean = nan(length(unique_stim),size(fluor_roi_diff,2),n_rois);
mua_stim_mean = nan(length(unique_stim),size(mua_allcat,2),n_depths);
for curr_stim_idx = 1:length(unique_stim)
    curr_trials = move_t > 0.5 & ...
        D_allcat.stimulus == unique_stim(curr_stim_idx);
    fluor_stim_mean(curr_stim_idx,:,:) = nanmean(fluor_allcat(curr_trials,:,:),1);
    fluor_roi_stim_mean(curr_stim_idx,:,:) = nanmean(fluor_roi_diff(curr_trials,:,:),1);
    mua_stim_mean(curr_stim_idx,:,:) = nanmean(mua_allcat(curr_trials,:,:),1);
end

figure('Name','Naive passive choiceworld');

% Plot 100/R fluorescence
subplot(2,3,1);
plot_px = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    permute(fluor_stim_mean(end,:,:),[3,2,1]));
plot_px_max = max(diff(plot_px(:,:,use_t),[],3),[],3);
imagesc(plot_px_max);
caxis([0,max(abs(caxis))]);
colormap(brewermap([],'BuGn'));
AP_reference_outline('ccf_aligned','k');
axis image off;

% Plot 100/R and contrast curve for visual cortex and striatal depths
subplot(2,3,2); hold on;
plot_rois = [1:3,7];
plot(t_diff,squeeze(fluor_roi_stim_mean(end,:,plot_rois)),'linewidth',2);
axis tight;
line([0,0],ylim,'color','k');
xlim([-0.2,0.5]);
legend({wf_roi(plot_rois).area});
ylabel('Fluorescence');

subplot(2,3,3); hold on; set(gca,'ColorOrder',copper(n_depths));
plot(t,squeeze(mua_stim_mean(end,:,:)),'linewidth',2);
axis tight;
line([0,0],ylim,'color','k');
xlim([-0.2,0.5]);
legend(cellfun(@num2str,num2cell(1:n_depths),'uni',false));
ylabel('MUA');

subplot(2,3,5); hold on;
fluor_roi_stim_max = squeeze(max(fluor_roi_stim_mean(:,use_t_diff,:),[],2));
plot(unique_stim,fluor_roi_stim_max(:,plot_rois),'linewidth',2);
xlabel('Contrast*Side');
ylabel('Fluorescence');

subplot(2,3,6); hold on; set(gca,'ColorOrder',copper(n_depths));
mua_stim_max = squeeze(max(mua_stim_mean(:,use_t,:),[],2));
plot(unique_stim,mua_stim_max,'linewidth',2);
xlabel('Contrast*Side');
ylabel('MUA');

%% Fig 4?: stim responses: passive fullscreen trained

data_fn = ['trial_activity_passive_fullscreen'];
exclude_data = false;

[fluor_allcat,fluor_roi_diff,mua_allcat,wheel_allcat,reward_allcat,D_allcat] = ...
    AP_load_concat_normalize_ctx_str(data_fn,exclude_data);

n_vs = size(fluor_allcat,3);
n_rois = size(fluor_roi_diff,3);
n_depths = size(mua_allcat,3);

% Get ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Get time (make this be saved in trial data)
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);
t_diff = conv(t,[1,1]/2,'valid');

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Stim and time to plot
unique_stim = unique(D_allcat.stimulus);
use_t = t > 0 & t < 0.5;
use_t_diff = t_diff > 0 & t_diff < 0.5;

% Get average long reaction time stim responses
fluor_stim_mean = nan(length(unique_stim),size(fluor_allcat,2),n_vs);
fluor_roi_stim_mean = nan(length(unique_stim),size(fluor_roi_diff,2),n_rois);
mua_stim_mean = nan(length(unique_stim),size(mua_allcat,2),n_depths);
for curr_stim_idx = 1:length(unique_stim)
    curr_trials = move_t > 0.5 & ...
        D_allcat.stimulus == unique_stim(curr_stim_idx);
    fluor_stim_mean(curr_stim_idx,:,:) = nanmean(fluor_allcat(curr_trials,:,:),1);
    fluor_roi_stim_mean(curr_stim_idx,:,:) = nanmean(fluor_roi_diff(curr_trials,:,:),1);
    mua_stim_mean(curr_stim_idx,:,:) = nanmean(mua_allcat(curr_trials,:,:),1);
end

figure('Name','Trained passive fullscreen');

% Plot 100/R fluorescence
subplot(2,3,1);
plot_px = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    permute(fluor_stim_mean(end,:,:),[3,2,1]));
plot_px_max = max(diff(plot_px(:,:,use_t),[],3),[],3);
imagesc(plot_px_max);
caxis([0,max(abs(caxis))]);
colormap(brewermap([],'BuGn'));
AP_reference_outline('ccf_aligned','k');
axis image off;

% Plot 100/R and contrast curve for visual cortex and striatal depths
subplot(2,3,2); hold on;
plot_rois = [1:3,7];
plot(t_diff,squeeze(fluor_roi_stim_mean(end,:,plot_rois)),'linewidth',2);
axis tight;
line([0,0],ylim,'color','k');
xlim([-0.2,0.5]);
legend({wf_roi(plot_rois).area});
ylabel('Fluorescence');

subplot(2,3,3); hold on; set(gca,'ColorOrder',copper(n_depths));
plot(t,squeeze(mua_stim_mean(end,:,:)),'linewidth',2);
axis tight;
line([0,0],ylim,'color','k');
xlim([-0.2,0.5]);
legend(cellfun(@num2str,num2cell(1:n_depths),'uni',false));
ylabel('MUA');

subplot(2,3,5); hold on;
fluor_roi_stim_max = squeeze(max(fluor_roi_stim_mean(:,use_t_diff,:),[],2));
plot(unique_stim,fluor_roi_stim_max(:,plot_rois),'linewidth',2);
set(gca,'XTick',1:3,'XTickLabel',{'Left','Center','Right'});
xlabel('Stim');
ylabel('Fluorescence');

subplot(2,3,6); hold on; set(gca,'ColorOrder',copper(n_depths));
mua_stim_max = squeeze(max(mua_stim_mean(:,use_t,:),[],2));
plot(unique_stim,mua_stim_max,'linewidth',2);
set(gca,'XTick',1:3,'XTickLabel',{'Left','Center','Right'});
xlabel('Stim');
ylabel('MUA');

%% Fig 4?: stim responses: passive fullscreen naive

data_fn = ['trial_activity_passive_fullscreen_naive'];
exclude_data = false;

[fluor_allcat,fluor_roi_diff,mua_allcat,wheel_allcat,reward_allcat,D_allcat] = ...
    AP_load_concat_normalize_ctx_str(data_fn,exclude_data);

n_vs = size(fluor_allcat,3);
n_rois = size(fluor_roi_diff,3);
n_depths = size(mua_allcat,3);

% Get ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Get time (make this be saved in trial data)
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);
t_diff = conv(t,[1,1]/2,'valid');

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Stim and time to plot
unique_stim = unique(D_allcat.stimulus);
use_t = t > 0 & t < 0.5;
use_t_diff = t_diff > 0 & t_diff < 0.5;

% Get average long reaction time stim responses
fluor_stim_mean = nan(length(unique_stim),size(fluor_allcat,2),n_vs);
fluor_roi_stim_mean = nan(length(unique_stim),size(fluor_roi_diff,2),n_rois);
mua_stim_mean = nan(length(unique_stim),size(mua_allcat,2),n_depths);
for curr_stim_idx = 1:length(unique_stim)
    curr_trials = move_t > 0.5 & ...
        D_allcat.stimulus == unique_stim(curr_stim_idx);
    fluor_stim_mean(curr_stim_idx,:,:) = nanmean(fluor_allcat(curr_trials,:,:),1);
    fluor_roi_stim_mean(curr_stim_idx,:,:) = nanmean(fluor_roi_diff(curr_trials,:,:),1);
    mua_stim_mean(curr_stim_idx,:,:) = nanmean(mua_allcat(curr_trials,:,:),1);
end

figure('Name','Naive passive fullscreen');

% Plot 100/R fluorescence
subplot(2,3,1);
plot_px = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    permute(fluor_stim_mean(end,:,:),[3,2,1]));
plot_px_max = max(diff(plot_px(:,:,use_t),[],3),[],3);
imagesc(plot_px_max);
caxis([0,max(abs(caxis))]);
colormap(brewermap([],'BuGn'));
AP_reference_outline('ccf_aligned','k');
axis image off;

% Plot 100/R and contrast curve for visual cortex and striatal depths
subplot(2,3,2); hold on;
plot_rois = [1:3,7];
plot(t_diff,squeeze(fluor_roi_stim_mean(end,:,plot_rois)),'linewidth',2);
axis tight;
line([0,0],ylim,'color','k');
xlim([-0.2,0.5]);
legend({wf_roi(plot_rois).area});
ylabel('Fluorescence');

subplot(2,3,3); hold on; set(gca,'ColorOrder',copper(n_depths));
plot(t,squeeze(mua_stim_mean(end,:,:)),'linewidth',2);
axis tight;
line([0,0],ylim,'color','k');
xlim([-0.2,0.5]);
legend(cellfun(@num2str,num2cell(1:n_depths),'uni',false));
ylabel('MUA');

subplot(2,3,5); hold on;
fluor_roi_stim_max = squeeze(max(fluor_roi_stim_mean(:,use_t_diff,:),[],2));
plot(unique_stim,fluor_roi_stim_max(:,plot_rois),'linewidth',2);
set(gca,'XTick',1:3,'XTickLabel',{'Left','Center','Right'});
xlabel('Stim');
ylabel('Fluorescence');

subplot(2,3,6); hold on; set(gca,'ColorOrder',copper(n_depths));
mua_stim_max = squeeze(max(mua_stim_mean(:,use_t,:),[],2));
plot(unique_stim,mua_stim_max,'linewidth',2);
set(gca,'XTick',1:3,'XTickLabel',{'Left','Center','Right'});
xlabel('Stim');
ylabel('MUA');





