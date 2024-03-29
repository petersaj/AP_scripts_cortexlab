% Generate figures for ctx-str paper

% Anything that takes a lot of time is done in
% AP_ctx_str_trial_preprocessing and saved for plotting here

% The original scripts here were in test_wf_ephys_choiceworld_analysis

%% Fig 1a: Behavior psychometric 
% (from AP_vanillaChoiceworld_behavior - currently no eliminations)

% Load behavior
bhv_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\bhv_processing\bhv.mat';
load(bhv_fn);

% Exclude bad behavior sessions
exclude_data = true;

bhv_fieldnames = fieldnames(bhv);
experiment_fields = cellfun(@(curr_field) ...
    length([bhv.(curr_field)]) == length([bhv.days]),bhv_fieldnames);

% Load pre-marked experiments to exclude and cut out bad ones
if exclude_data
    exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
    exclude_fn{1} = 'bhv_use_experiments';
    % exclude_fn{2} = 'expl_var_use_experiments';
    use_experiments_all = {};
    for curr_exclude = 1:length(exclude_fn)
        curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
        use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
    end
    use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
        1:size(use_experiments_all,2),'uni',false)';
    
    % Cut out bad experiments for any experiment data fields
    for curr_field = bhv_fieldnames(experiment_fields)'
        for curr_animal = 1:length(use_experiments)
            bhv(curr_animal).(cell2mat(curr_field)) = ...
                bhv(curr_animal).(cell2mat(curr_field))(use_experiments{curr_animal});
        end
    end
end


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

%% Fig 1 b/c: Cortical activity and regression during task

% Load data
data_fn = ['trial_activity_choiceworld_DECONVTEST'];
exclude_data = true;
AP_load_concat_normalize_ctx_str;

n_vs = size(fluor_allcat_deconv,3);
n_depths = size(mua_allcat,3);

% Load regression results
regression_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\trial_activity_choiceworld_task_regression';
load(regression_fn);

% Get trial information
trial_contrast_allcat = max(D_allcat.stimulus,[],2);
[~,side_idx] = max(D_allcat.stimulus > 0,[],2);
trial_side_allcat = (side_idx-1.5)*2;
trial_contrastside_allcat = trial_contrast_allcat.*trial_side_allcat;
trial_choice_allcat = -(D_allcat.response-1.5)*2;

% % Get time (make this be saved in trial data)
% framerate = 35;
% raster_window = [-0.5,3];
% upsample_factor = 3;
% sample_rate = (framerate*upsample_factor);
% t = raster_window(1):1/sample_rate:raster_window(2);
sample_rate = 1/mean(diff(t));

% Get reaction time
[move_trial,move_idx] = max(abs(wheel_allcat) > 0.02,[],2);
move_t = t(move_idx)';
move_t(~move_trial) = Inf;

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Load widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);
n_rois = numel(wf_roi);

% Get predicted fluorescence in ROIs
fluor_predicted_roi = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_predicted,[3,2,1]),n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
    n_rois,[],size(fluor_allcat_predicted,1)),[3,2,1]);

fluor_predicted_reduced_roi = arrayfun(@(x) ...
    permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_predicted_reduced(:,:,:,x),[3,2,1]), ...
    n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
    n_rois,[],size(fluor_allcat_predicted,1)),[3,2,1]), ...
    1:size(fluor_allcat_predicted_reduced,4),'uni',false);

% Plot map of cortical activity at time points on correct visual trials
plot_rxn_time = [0.1,0.3];

vis_correct_L_trials = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == -1 & ...
    move_t >= plot_rxn_time(1) & ...
    move_t <= plot_rxn_time(2);

vis_correct_L_px = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    squeeze(nanmean(fluor_allcat_deconv(vis_correct_L_trials,:,:),1))');

AP_imscroll(vis_correct_L_px,t);
axis image;
caxis([0,max(abs(caxis))])
colormap(brewermap([],'BuGn'));
AP_reference_outline('ccf_aligned','k');

plot_t = [find(t > 0.06,1),find(t > 0.120,1), ...
    find(t > 0.3,1),find(t > 0.8,1)];

figure;
colormap(brewermap([],'BuGn'));
for curr_t = 1:length(plot_t)
    subplot(1,length(plot_t),curr_t);
    imagesc(vis_correct_L_px(:,:,plot_t(curr_t)));
    caxis([0,max(abs(caxis))]);
    axis image off;
    AP_reference_outline('ccf_aligned','k');
    title([num2str(round(t(plot_t(curr_t))*1000)) ' ms after stim']);
end

% Plot measured and predicted relating to each kernel
% % (stim = V1, move onset/ongoing = SMl, beep = PPC, reward = SMf)
% plot_areas = {'V1p_L','RSPa_L','SMl_L','PPC_L','SMf_L'};
% plot_reduction = [1,2,3,4,5]; % (which reduced variable to plot)
% plot_conditions = {unique(trial_contrastside_allcat),[-1,1],1:3,[1,3],0:1};
% plot_conditions_compare = { ...
%     trial_contrastside_allcat, ...  
%     trial_choice_allcat, ...
%     discretize(move_t,linspace(0,0.5,4)), ...
%     discretize(move_t,[0.2,0.3,0.6,0.7]), ...
%     trial_choice_allcat == -trial_side_allcat};
% 
% stim_col = colormap_BlueWhiteRed(5);
% stim_col(6,:) = 0;
% plot_cols = { ...
%     stim_col, ...
%     [1,0,0;0,0,1], ...
%     copper(3), ...
%     [0.5,0.5,0;0.6,0,0.6], ...
%     [0,0,0;0,0,0.7]};

% (without move ongoing)
% (stim = V1, move onset = SMl, beep = PPC, reward = SMf)
plot_areas = {'V1p_L','RSPa_L','SMl_L','PPC_L','SMf_L'};
plot_reduction = [1,2,2,3,4]; % (which reduced variable to plot)
plot_conditions = {unique(trial_contrastside_allcat),[-1,1],1:3,[1,3],0:1};
plot_conditions_compare = { ...
    trial_contrastside_allcat, ...  
    trial_choice_allcat, ...
    discretize(move_t,linspace(0,0.5,4)), ...
    discretize(move_t,[0.2,0.3,0.6,0.7]), ...
    trial_choice_allcat == -trial_side_allcat};

stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = 0;
plot_cols = { ...
    stim_col, ...
    [1,0,0;0,0,1], ...
    copper(3), ...
    [0.5,0.5,0;0.6,0,0.6], ...
    [0,0,0;0,0,0.7]};

figure;
for curr_area = 1:length(plot_areas)
    area_idx = strcmp(plot_areas{curr_area},{wf_roi.area});

    p1 = subplot(5,length(plot_areas),curr_area);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(fluor_roi_deconv(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel([wf_roi(area_idx).area ' Measured']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p2 = subplot(5,length(plot_areas),curr_area+length(plot_areas));
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(fluor_predicted_roi(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel([wf_roi(area_idx).area ' Predicted']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p3 = subplot(5,length(plot_areas),curr_area+length(plot_areas)*2);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(fluor_roi_deconv(curr_trials,:,area_idx) - ...
            fluor_predicted_roi(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel([wf_roi(area_idx).area ' Residual']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p4 = subplot(5,length(plot_areas),curr_area+length(plot_areas)*3);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean( ...
            fluor_predicted_reduced_roi{plot_reduction(curr_area)}( ...
            curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel([wf_roi(area_idx).area ' Reduced (' ...
        regressor_labels{plot_reduction(curr_area)} ')']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p5 = subplot(5,length(plot_areas),curr_area+length(plot_areas)*4);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean( fluor_roi_deconv(curr_trials,:,area_idx) - ...
            fluor_predicted_reduced_roi{plot_reduction(curr_area)}( ...
            curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel([wf_roi(area_idx).area ' Reduced residual (' ...
        regressor_labels{plot_reduction(curr_area)} ')']);
    axis tight
    line([0,0],ylim,'color','k');
       
    linkaxes([p1,p2,p3,p4,p5],'xy');
    
end

% Plot max weights aross regressors
max_k_px = nan(size(U_master,1),size(U_master,2),length(regressors));
for curr_regressor = 1:length(regressors)
    curr_k_v = permute(cell2mat(permute(fluor_kernel(curr_regressor,1:n_vs),[1,3,2])),[3,2,1]);
    curr_k_px = reshape(svdFrameReconstruct(U_master(:,:,1:n_vs), ...
        reshape(curr_k_v,n_vs,[])),size(U_master,1),size(U_master,2),[],size(curr_k_v,3));
    max_k_px(:,:,curr_regressor) = max(max(curr_k_px,[],3),[],4);
end

figure;
max_weight = prctile(max_k_px(:),99);
for curr_regressor = 1:length(regressors)
    subplot(1,length(regressors),curr_regressor);
    colormap(brewermap([],'Purples'));
    imagesc(max_k_px(:,:,curr_regressor));
    AP_reference_outline('ccf_aligned','k');
    title(regressor_labels{curr_regressor});
    axis image off;
    caxis([0,max_weight]);
end

% Plot long regression example?
plot_areas = {'V1p_L','SMl_L'};
figure;
for curr_area = 1:length(plot_areas)
    plot_area_idx = strcmp(plot_areas{curr_area},{wf_roi.area});
    long_trace = reshape(permute(fluor_roi_deconv(:,:,plot_area_idx),[2,1]),[],1);
    long_trace_predicted = reshape(permute(fluor_predicted_roi(:,:,plot_area_idx),[2,1]),[],1);
    
    % (get rid of NaNs)
    nan_samples = isnan(long_trace) | isnan(long_trace_predicted);
    long_trace(nan_samples) = [];
    long_trace_predicted(nan_samples) = [];
    
    % (correlate measured and predicted in chunks)
    n_chunks = 1000;
    chunk_boundaries = round(linspace(1,length(long_trace),n_chunks+1));
    chunk_corr = nan(n_chunks,1);
    for curr_chunk = 1:n_chunks
        curr_corr = corrcoef(long_trace(chunk_boundaries(curr_chunk):chunk_boundaries(curr_chunk+1)), ...
            long_trace_predicted(chunk_boundaries(curr_chunk):chunk_boundaries(curr_chunk+1)));
        chunk_corr(curr_chunk) = curr_corr(2);
    end
    
    % (plot chunk percentiles)
    plot_prctiles = [0,25,50,75,100];
    [chunk_corr_sorted,sort_idx] = sort(chunk_corr);
    chunk_corr_prctile = prctile(chunk_corr,plot_prctiles);
    for curr_prctile = 1:length(chunk_corr_prctile)
        subplot(length(chunk_corr_prctile),length(plot_areas), ...
            (curr_prctile-1)*length(plot_areas)+curr_area); hold on;
        curr_plot_chunk_idx = sort_idx(find(chunk_corr_sorted >= chunk_corr_prctile(curr_prctile),1));
        curr_chunk_plot = chunk_boundaries(curr_plot_chunk_idx):chunk_boundaries(curr_plot_chunk_idx+1);
        plot([1:length(curr_chunk_plot)]/sample_rate,long_trace(curr_chunk_plot),'color',[0,0.6,0]);
        plot([1:length(curr_chunk_plot)]/sample_rate,long_trace_predicted(curr_chunk_plot),'color',[0.6,0,0.6]);
        ylabel([num2str(plot_prctiles(curr_prctile)) 'th percentile']);
        if curr_prctile == 1
            title(plot_areas{curr_area})
        end
    end
end

% Plot widefield ROIs
figure; hold on
set(gca,'YDir','reverse');
AP_reference_outline('ccf_aligned','k');
for curr_roi = 1:n_rois
    curr_roi_boundary = cell2mat(bwboundaries(wf_roi(curr_roi).mask));
    patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),[0,0.7,0]);
    
    text(nanmean(curr_roi_boundary(:,2)),nanmean(curr_roi_boundary(:,1)), ...
        wf_roi(curr_roi).area,'FontSize',10,'HorizontalAlignment','center','color','w')
end
axis image off;


% Storing this here for now? Side differences, but move this later

% % Get max stim kernel and contrast slope maps
% stim_regressor = strcmp(regressor_labels,'Stim');
% k_v_stim = permute(cell2mat(permute(fluor_kernel(stim_regressor,1:n_vs),[1,3,2])),[3,2,1]);
% k_px_stim = reshape(svdFrameReconstruct(U_master(:,:,1:n_vs), ...
%     reshape(k_v_stim,n_vs,[])),size(U_master,1),size(U_master,2),[],size(k_v_stim,3));
% k_px_stim_maxt = squeeze(max(k_px_stim,[],3));
% 
% contrasts = [0.06,0.125,0.25,0.5,1];
% 
% contrast_k_L = AP_regresskernel(contrasts, ...
%     reshape(k_px_stim_maxt(:,:,5:-1:1),[],5),0,[],[false,false],1,true,true);
% contrast_slope_L = reshape(contrast_k_L{1},size(U_master,1),size(U_master,2));
% 
% contrast_k_R = AP_regresskernel(contrasts, ...
%     reshape(k_px_stim_maxt(:,:,6:10),[],5),0,[],[false,false],1,true,true);
% contrast_slope_R = reshape(contrast_k_R{1},size(U_master,1),size(U_master,2));
% 
% contrast_slope_diff = contrast_slope_L - contrast_slope_R;
% 
% % Get max move kernel and move ipsi-contra differences
% move_onset_regressor = strcmp(regressor_labels,'Move onset');
% k_v_move_onset = permute(cell2mat(permute(fluor_kernel(move_onset_regressor,1:n_vs),[1,3,2])),[3,2,1]);
% k_px_move_onset = reshape(svdFrameReconstruct(U_master(:,:,1:n_vs), ...
%     reshape(k_v_move_onset,n_vs,[])),size(U_master,1),size(U_master,2),[],size(k_v_move_onset,3));
% k_px_move_onset_maxt = squeeze(max(k_px_move_onset,[],3));
% 
% k_px_move_onset_maxt_diff = k_px_move_onset_maxt(:,:,1)-k_px_move_onset_maxt(:,:,2);
% 
% figure;
% subplot(2,2,1);
% imagesc(max(k_px_stim_maxt,[],3))
% axis image off;
% AP_reference_outline('ccf_aligned','k');
% caxis([0,max(k_px_stim_maxt(:))])
% colormap(gca,brewermap([],'Purples'));
% title('Max stim kernel weight');
% 
% subplot(2,2,2);
% imagesc(contrast_slope_diff); 
% axis image off;
% AP_reference_outline('ccf_aligned','k');
% caxis([-max(k_px_stim_maxt(:)),max(k_px_stim_maxt(:))])
% colormap(gca,brewermap([],'RdBu'));
% title('Contrast slope R-L');
% 
% subplot(2,2,3);
% imagesc(max(k_px_move_onset_maxt,[],3));
% axis image off;
% AP_reference_outline('ccf_aligned','k');
% caxis([0,max(k_px_move_onset_maxt(:))])
% colormap(gca,brewermap([],'Purples'));
% title('Max move onset kernel weight');
% 
% subplot(2,2,4);
% imagesc(k_px_move_onset_maxt_diff);
% axis image off;
% AP_reference_outline('ccf_aligned','k');
% caxis([-max(k_px_move_onset_maxt(:)),max(k_px_move_onset_maxt(:))])
% colormap(gca,brewermap([],'*RdBu'));
% title('Move L-R');
% 
% % Plot contrast response function in ROIs 
% contrast_sides = sort(reshape([0.06,0.125,0.25,0.5,1].*[-1,1]',[],1));
% 
% stim_max_k_roi = cell2mat(arrayfun(@(x) ....
%     nanmean(reshape(reshape(repmat(roi_mask(:,:,x),1,1,10).* ...
%     k_px_stim_maxt,[],10),[],10),1),1:size(roi_mask,3),'uni',false)');
% 
% figure;
% AP_stackplot(stim_max_k_roi',contrast_sides, ...
%     range(stim_max_k_roi(:))*1.2,false,[0,0.7,0],{wf_roi.area},true);
% xlabel('Contrast*Side');
% title('Stim kernel maximum');



%% Fig 2b: Example traces/kernels

% Load and align
str_align = 'kernel';

% animal = 'AP028'; 
% day = '2017-12-20'; % 16/20 (blood in 16, a little in 20)

% animal = 'AP027';
% day = '2017-11-25'; % some of those pco interger type artifacts

animal = 'AP025';
day = '2017-10-04';

experiment = 1; 
verbose = true; 
AP_load_experiment;

avg_im_aligned = AP_align_widefield(animal,day,avg_im);
Udf_aligned = single(AP_align_widefield(animal,day,Udf));

% Parameters for regression
n_aligned_depths = 4;
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 60;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.2,0.2];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

% Get lambda
lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
load(lambda_fn);
curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
if any(curr_animal_idx)
    curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
    if any(curr_day_idx)
        lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
    end
end

% Get time points to bin
sample_rate = framerate*regression_params.upsample_factor;
time_bins = frame_t(find(frame_t > ...
    regression_params.skip_seconds,1)):1/sample_rate: ...
    frame_t(find(frame_t-frame_t(end) < ...
    -regression_params.skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Get upsampled dVdf's
deriv_smooth_factor = 3;
dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'),diff(convn(fVdf, ...
    ones(1,deriv_smooth_factor)/deriv_smooth_factor,'same'),[],2)',time_bin_centers)';

% % %%%%%%%%%%% TESTING DECONV
% 
% % Deconvolve from the cortical recordings kernel
% load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\gcamp_kernel\gcamp6s_kernel.mat');
% 
% gcamp6s_kernel_cat = vertcat(gcamp6s_kernel.regression{:});
% gcamp6s_kernel = nanmean(gcamp6s_kernel_cat./max(gcamp6s_kernel_cat,[],2),1);
% 
% fVdf_deconv = convn(fVdf,gcamp6s_kernel,'same');
%
% dfVdf_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
% 
% % %%%%%%%%%%%%

% Get striatum depth group by across-experiment alignment
n_depths = 4;
depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
depth_group = discretize(spike_depths,depth_group_edges);
use_depths = 1:n_depths;

binned_spikes = zeros(n_depths,length(time_bin_centers));
for curr_depth = 1:n_depths
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
end

% Define ROIs and get fluorescence traces
roi_circle_size = 10;
roi_x = [136,176,142,116]; % roi_x = [131,174,110,51];
roi_y = [299,91,79,87]; % roi_y = [297,96,71,144];
[x,y] = meshgrid(1:size(avg_im_aligned,1),1:size(avg_im_aligned,2));
roi_mask = cell2mat(arrayfun(@(roi) sqrt((x-roi_x(roi)).^2 + (y-roi_y(roi)).^2) <= ...
    roi_circle_size,permute(1:length(roi_x),[1,3,2]),'uni',false));
roi_trace_deriv = AP_svd_roi(Udf_aligned,dfVdf_resample,[],[],roi_mask);

% Get spike-triggered average
surround_times = [-0.2,0.2];

surround_frames = round(surround_times(1)*framerate):round(surround_times(2)*framerate);
sta_t = surround_frames./framerate;
sta_im = zeros(size(Udf_aligned,1),size(Udf_aligned,2), ...
    length(surround_frames),size(binned_spikes,1));

for curr_depth = 1:size(binned_spikes,1)
    frames_w = repmat(binned_spikes(curr_depth,:)'./ ...
        sum(binned_spikes(curr_depth,:)),1,length(surround_frames));
    for curr_sta_frame = 1:length(surround_frames)
        frames_w(:,curr_sta_frame) = ...
            circshift(frames_w(:,curr_sta_frame),[surround_frames(curr_sta_frame),0]);
    end
    
    sta_v = dfVdf_resample*frames_w;
    sta_im(:,:,:,curr_depth) = svdFrameReconstruct(Udf_aligned,sta_v);
end

sta_im_max = squeeze(max(sta_im,[],3));

% Regress fluorescence to spikes
kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
    round(regression_params.kernel_t(2)*sample_rate);

binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
binned_spikes_std(isnan(binned_spikes_std)) = 0;

[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(dfVdf_resample(regression_params.use_svs,:), ...
    binned_spikes_std,kernel_frames,lambda, ...
    regression_params.zs,regression_params.cvfold, ...
    false,regression_params.use_constant);

Udf_aligned = single(AP_align_widefield(animal,day,Udf));
k_px = zeros(size(Udf_aligned,1),size(Udf_aligned,2),size(k,2),size(k,3),'single');
for curr_spikes = 1:size(k,3)
    k_px(:,:,:,curr_spikes) = ...
        svdFrameReconstruct(Udf_aligned(:,:,regression_params.use_svs),k(:,:,curr_spikes));
end

k_px_max = squeeze(max(k_px,[],3));
k_px_t0 = squeeze(k_px(:,:,kernel_frames == 0,:));

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
    imagesc(k_px_t0(:,:,curr_depth));
    caxis([0,prctile(k_px_t0(:),99)]);
    colormap(h,brewermap([],'Purples'));
    AP_reference_outline('ccf_aligned','k');
    axis image off;
end

% Plot ROIs and traces

% Smoothing window
smooth_size = 9; % MUST BE ODD
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

binned_spikes_smoothed = convn(binned_spikes,smWin,'same');
predicted_spikes_smoothed = convn(predicted_spikes,smWin,'same');

figure;

subplot(1,3,1);
roi_boundaries = bwboundaries(sum(roi_mask,3));
imagesc(avg_im_aligned);colormap(gray);
caxis([0,prctile(avg_im_aligned(:),99)]);
axis image off;
AP_reference_outline('ccf_aligned','r');
p = cellfun(@(x) plot(x(:,2),x(:,1),'b','linewidth',2),roi_boundaries);

p1 = subplot(2,3,2:3); hold on;
AP_stackplot(roi_trace_deriv',time_bin_centers,5,true,[0,0.7,0]);
title('\DeltaFluorescence')
xlabel('Time (seconds)');
ylabel('Activity (std)');

p2 = subplot(2,3,5:6); hold on;
AP_stackplot(binned_spikes_smoothed',time_bin_centers,5,true,'k');
AP_stackplot(predicted_spikes_smoothed',time_bin_centers,5,true,[0.6,0,0.6]);
title('MUA')
xlabel('Time (seconds)');
ylabel('Activity (std)');

linkaxes([p1,p2],'xy');
axis tight;
xlim([395,433]);


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

AP_imscroll(k_px_trained,t);
caxis([-max(abs(caxis)),max(abs(caxis))]);
axis image;
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k');
set(gcf,'Name','Trained');

AP_imscroll(k_px_naive,t);
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

AP_imscroll(k_px_trained_com_colored_weighted,t,true);
axis image;
AP_reference_outline('ccf_aligned','k');
set(gcf,'Name','Trained');

AP_imscroll(k_px_naive_com_colored_weighted,t,true);
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
plot_t = find(t >= -0.06 & t <= 0.06);
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
    k_fn = [data_path filesep 'wf_ephys_maps_' curr_protocol '_' num2str(n_aligned_depths) '_depths_kernel_DECONVTEST'];
    load(k_fn);
    
    framerate = 35;
    upsample_factor = 1;
    sample_rate = framerate*upsample_factor;
    kernel_t = [-0.5,0.5];
    kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
    t = kernel_frames/sample_rate;
    
    % Concatenate explained variance
    expl_var_experiment = cell2mat(horzcat(batch_vars.explained_var));
    expl_var_animal = cell2mat(cellfun(@(x) nanmean(cell2mat(x),2),{batch_vars.explained_var},'uni',false));
    figure;errorbar(nanmean(expl_var_experiment,2), ...
        nanstd(expl_var_experiment,[],2)./sqrt(nansum(expl_var_experiment,2)),'k','linewidth',2);
    
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
    
    % Plot at kernel_frames = 1, transparency relative to weight
    figure;
    plot_frame = kernel_frames == 0;
    p = image(k_px_com_colored(:,:,:,plot_frame));
    % weight_max = max(k_px(:))*0.8;
    weight_max = 0.005;
    set(p,'AlphaData', ...
        mat2gray(max(k_px(:,:,plot_frame,:),[],4),[0,weight_max]));
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

%% Fig 3?: (OLD - remove?) Striatal activity during task

data_fn = ['trial_activity_choiceworld'];
exclude_data = true;

[t,fluor_allcat_deconv,fluor_roi_deconv,mua_allcat,wheel_allcat,reward_allcat,D_allcat] = ...
    AP_load_concat_normalize_ctx_str(data_fn,exclude_data);

n_depths = size(fluor_allcat_deconv,3);

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
[move_trial,move_idx] = max(abs(wheel_allcat) > 0.02,[],2);
move_t = t(move_idx)';
move_t(~move_trial) = Inf;

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


%% Fig 3?: (OLD - remove?) Striatal activity and regression during task

% Load data
data_fn = ['trial_activity_choiceworld_framerate'];
exclude_data = true;
[t,fluor_allcat_deconv,fluor_roi_deconv,mua_allcat,wheel_allcat,reward_allcat,D_allcat] = ...
    AP_load_concat_normalize_ctx_str(data_fn,exclude_data);

% Load regression results
regression_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\trial_activity_choiceworld_task_regression';
load(regression_fn);

% Get trial information
trial_contrast_allcat = max(D_allcat.stimulus,[],2);
[~,side_idx] = max(D_allcat.stimulus > 0,[],2);
trial_side_allcat = (side_idx-1.5)*2;
trial_choice_allcat = -(D_allcat.response-1.5)*2;

% Get number of V's/depths
n_vs = size(fluor_allcat_predicted,3);
n_depths = size(mua_allcat_predicted,3);

% Get time
sample_rate = 1/mean(diff(t));

% Get move onset index
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Get time shifts in samples
sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
    round(x(2)*(sample_rate)),t_shifts,'uni',false);

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
            sample_shifts{curr_regressor}/(sample_rate), ...
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
            plot_data = mua_allcat;
            plot_title = 'Measured';
        case 2
            plot_data = mua_allcat_predicted;
            plot_title = 'Predicted';
        case 3
            plot_data = mua_allcat - mua_allcat_predicted;
            plot_title = 'Residual';
    end
    
    p(curr_plot) = subplot(1,3,curr_plot); hold on;
    for curr_plot_condition = 1:size(plot_conditions,1)
        
        curr_trials = plot_id == curr_plot_condition & ...
            move_t >= rxn_time_use(1) & ...
            move_t <= rxn_time_use(2);
        curr_data = plot_data(curr_trials,:,:);
        
%         % Re-align to movement onset
%         t_leeway = -t_downsample_deriv(1);
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
        
        AP_stackplot(curr_data_mean,t,2,false,curr_col);
        
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

%% Fig 3?: (NEW - replace above?) Striatal activity and regression during task

% Load data
data_fn = ['trial_activity_choiceworld_DECONVTEST'];
exclude_data = true;

AP_load_concat_normalize_ctx_str;
n_vs = size(fluor_allcat_deconv,3);
n_depths = size(mua_allcat,3);

% Load regression results
regression_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\trial_activity_choiceworld_task_regression';
load(regression_fn);

% Get trial information
trial_contrast_allcat = max(D_allcat.stimulus,[],2);
[~,side_idx] = max(D_allcat.stimulus > 0,[],2);
trial_side_allcat = (side_idx-1.5)*2;
trial_contrastside_allcat = trial_contrast_allcat.*trial_side_allcat;
trial_choice_allcat = -(D_allcat.response-1.5)*2;

% % Get time
sample_rate = 1/mean(diff(t));

% Get reaction time
[move_trial,move_idx] = max(abs(wheel_allcat) > 0.02,[],2);
move_t = t(move_idx)';
move_t(~move_trial) = Inf;

% Plot striatal activity on correct visual trials
plot_rxn_time = [0.1,0.3];

vis_correct_L_trials = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == -1 & ...
    move_t >= plot_rxn_time(1) & ...
    move_t <= plot_rxn_time(2);

vis_correct_L_mua = squeeze(nanmean(mua_allcat(vis_correct_L_trials,:,:),1));

figure; hold on
set(gca,'ColorOrder',copper(n_depths));
plot(t,vis_correct_L_mua,'linewidth',2)
axis tight; 
line([0,0],ylim);
ylabel({'MUA (std)'});
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false));
xlabel('Time from stim onset');

% Plot measured and predicted relating to each kernel
% plot_areas = [1,2,3,1,4];
% plot_reduction = [1,2,3,4,5]; % (which reduced variable to plot)
% plot_conditions = {unique(trial_contrastside_allcat),[-1,1],1:3,[1,3],0:1};
% plot_conditions_compare = { ...
%     trial_contrastside_allcat, ...  
%     trial_choice_allcat, ...
%     discretize(move_t,linspace(0,0.5,4)), ...
%     discretize(move_t,[0.2,0.3,0.6,0.7]), ...
%     trial_choice_allcat == -trial_side_allcat};
% 
% stim_col = colormap_BlueWhiteRed(5);
% stim_col(6,:) = 0;
% plot_cols = { ...
%     stim_col, ...
%     [1,0,0;0,0,1], ...
%     copper(3), ...
%     [0.5,0.5,0;0.6,0,0.6], ...
%     [0,0,0;0,0,0.7]};

% (without move ongoing kernel)
plot_areas = [1,2,3,1,4];
plot_reduction = [1,2,2,3,4]; % (which reduced variable to plot)
plot_conditions = {unique(trial_contrastside_allcat),[-1,1],1:3,[1,3],0:1};
plot_conditions_compare = { ...
    trial_contrastside_allcat, ...  
    trial_choice_allcat, ...
    discretize(move_t,linspace(0,0.5,4)), ...
    discretize(move_t,[0.2,0.3,0.6,0.7]), ...
    trial_choice_allcat == -trial_side_allcat};

stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = 0;
plot_cols = { ...
    stim_col, ...
    [1,0,0;0,0,1], ...
    copper(3), ...
    [0.5,0.5,0;0.6,0,0.6], ...
    [0,0,0;0,0,0.7]};


figure;
for curr_area = 1:length(plot_areas)
    area_idx = plot_areas(curr_area);

    p1 = subplot(6,length(plot_areas),curr_area);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(mua_allcat(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Str ' num2str(area_idx) ' Measured']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p2 = subplot(6,length(plot_areas),curr_area+length(plot_areas));
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(mua_allcat_predicted(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Str ' num2str(area_idx) ' Predicted']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p3 = subplot(6,length(plot_areas),curr_area+length(plot_areas)*2);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(mua_allcat(curr_trials,:,area_idx) - ...
            mua_allcat_predicted(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Str ' num2str(area_idx) ' Residual']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p4 = subplot(6,length(plot_areas),curr_area+length(plot_areas)*3);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean( ...
            mua_allcat_predicted_reduced( ...
            curr_trials,:,area_idx,plot_reduction(curr_area)),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Str ' num2str(area_idx) ' Reduced (' ...
        regressor_labels{plot_reduction(curr_area)} ')']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p5 = subplot(6,length(plot_areas),curr_area+length(plot_areas)*4);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(mua_allcat(curr_trials,:,area_idx) - ...
            mua_allcat_predicted_reduced( ...
            curr_trials,:,area_idx,plot_reduction(curr_area)),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Str ' num2str(area_idx) ' Reduced residual (' ...
        regressor_labels{plot_reduction(curr_area)} ')']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p6 = subplot(6,length(plot_areas),curr_area+length(plot_areas)*5);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(predicted_mua_std_allcat(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Str ' num2str(area_idx) ' Cortex-predicted']);
    axis tight
    line([0,0],ylim,'color','k');
       
    linkaxes([p1,p2,p3,p4,p5,p6],'xy');
    
end

% Plot max weights aross regressors
max_k = nan(n_depths,length(regressors));
for curr_regressor = 1:length(regressors)
    curr_k = permute(cell2mat(permute(mua_kernel(curr_regressor,:),[1,3,2])),[3,2,1]);
    max_k(:,curr_regressor) = max(max(curr_k,[],2),[],3);
end

figure; hold on;
set(gca,'YDir','reverse', ...
    'ColorOrder',[0,0,0;0,0.8,0;0,0.5,0;0.7,0,0.7;0,0,0.7]);
plot(max_k,1:n_depths,'linewidth',2)
legend(regressor_labels)
ylabel('Striatum depth');
xlabel('Maximum kernel weight');
set(gca,'YTick',1:n_depths)

% Plot max weights over time across regressors
regressor_t = cellfun(@(x) (round(x(1)*(sample_rate)): ...
    round(x(2)*(sample_rate)))/sample_rate,t_shifts,'uni',false);

figure;
for curr_regressor = 1:length(regressors)
    curr_k = permute(cell2mat(permute(mua_kernel(curr_regressor,:),[1,3,2])),[3,2,1]);
    curr_k_maxt = max(curr_k,[],3);
    
    subplot(1,length(regressors),curr_regressor); hold on;
    set(gca,'ColorOrder',copper(n_depths));
    plot(regressor_t{curr_regressor},curr_k_maxt','linewidth',2);   
    title(regressor_labels{curr_regressor});
end

% Plot long regression example?
plot_areas = 1:4;
figure;
for curr_area = 1:length(plot_areas)
    plot_area_idx = plot_areas(curr_area);
    long_trace = reshape(permute(mua_allcat(:,:,plot_area_idx),[2,1]),[],1);
    long_trace_predicted = reshape(permute(mua_allcat_predicted(:,:,plot_area_idx),[2,1]),[],1);
    
%     % Filter long trace
%     lowpassCutoff = 6; % Hz
%     [b100s, a100s] = butter(2, lowpassCutoff/((sample_rate)/2), 'low');
%     long_trace = filter(b100s,a100s,long_trace,[],1);
        
    % (correlate measured and predicted in chunks)
    n_chunks = 1000;
    chunk_boundaries = round(linspace(1,length(long_trace),n_chunks+1));
    chunk_corr = nan(n_chunks,1);
    for curr_chunk = 1:n_chunks
%         curr_corr = corrcoef(long_trace(chunk_boundaries(curr_chunk):chunk_boundaries(curr_chunk+1)), ...
%             long_trace_predicted(chunk_boundaries(curr_chunk):chunk_boundaries(curr_chunk+1)));
%         chunk_corr(curr_chunk) = curr_corr(2);

        curr_corr = sum((long_trace(chunk_boundaries(curr_chunk):chunk_boundaries(curr_chunk+1))- ...
            long_trace_predicted(chunk_boundaries(curr_chunk):chunk_boundaries(curr_chunk+1))).^2);
        chunk_corr(curr_chunk) = curr_corr;       
    end
    
    % (plot chunk percentiles)
    plot_prctiles = [0,25,50,75,100];
    [chunk_corr_sorted,sort_idx] = sort(chunk_corr);
    chunk_corr_prctile = prctile(chunk_corr,plot_prctiles);
    for curr_prctile = 1:length(chunk_corr_prctile)
        subplot(length(chunk_corr_prctile),length(plot_areas), ...
            (curr_prctile-1)*length(plot_areas)+curr_area); hold on;
        curr_plot_chunk_idx = sort_idx(find(chunk_corr_sorted >= chunk_corr_prctile(curr_prctile),1));
        curr_chunk_plot = chunk_boundaries(curr_plot_chunk_idx):chunk_boundaries(curr_plot_chunk_idx+1);
        plot([1:length(curr_chunk_plot)]/sample_rate,long_trace(curr_chunk_plot),'color',[0,0.6,0]);
        plot([1:length(curr_chunk_plot)]/sample_rate,long_trace_predicted(curr_chunk_plot),'color',[0.6,0,0.6]);
        ylabel([num2str(plot_prctiles(curr_prctile)) 'th percentile']);
        if curr_prctile == 1
            title(plot_areas(curr_area))
        end
    end
end


%% ? ctx-predicted str stuff

% Load data
data_fn = ['trial_activity_choiceworld_DECONVTEST'];
exclude_data = true;
AP_load_concat_normalize_ctx_str;

n_vs = size(fluor_allcat_deconv,3);
n_depths = size(mua_allcat,3);

% Load regression results
regression_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\trial_activity_choiceworld_task_regression';
load(regression_fn)

% Get trial information
trial_contrast_allcat = max(D_allcat.stimulus,[],2);
[~,side_idx] = max(D_allcat.stimulus > 0,[],2);
trial_side_allcat = (side_idx-1.5)*2;
trial_contrastside_allcat = trial_contrast_allcat.*trial_side_allcat;
trial_choice_allcat = -(D_allcat.response-1.5)*2;

% % Get time (make this be saved in trial data)
sample_rate = 1/mean(diff(t));

% Get reaction time
[move_trial,move_idx] = max(abs(wheel_allcat) > 0.02,[],2);
move_t = t(move_idx)';
move_t(~move_trial) = Inf;

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');


%%%%%% Get cortical activity in ROIs similar to kernel

% Load widefield kernel ROIs
kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
load(kernel_roi_fn);

% Get predicted fluorescence in ROIs
fluor_kernel_roi_bw = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],kernel_roi.bw), ...
    size(kernel_roi.bw,3),[],size(fluor_allcat_deconv,1)),[3,2,1]);

fluor_kernel_roi_weighted = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],kernel_roi.max_weighted), ...
    size(kernel_roi.max_weighted,3),[],size(fluor_allcat_deconv,1)),[3,2,1]);


%%%%%% Plot measured and predicted relating to each kernel
plot_areas = [1,2,3,1,4];
plot_reduction = [1,2,3,4,5]; % (which reduced variable to plot)
plot_conditions = {unique(trial_contrastside_allcat),[-1,1],1:3,[1,3],0:1};
plot_conditions_compare = { ...
    trial_contrastside_allcat, ...  
    trial_choice_allcat, ...
    discretize(move_t,linspace(0,0.5,4)), ...
    discretize(move_t,[0.2,0.3,0.6,0.7]), ...
    trial_choice_allcat == -trial_side_allcat};

stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = 0;
plot_cols = { ...
    stim_col, ...
    [1,0,0;0,0,1], ...
    copper(3), ...
    [0.5,0.5,0;0.6,0,0.6], ...
    [0,0,0;0,0,0.7]};

figure;
for curr_area = 1:length(plot_areas)
    area_idx = plot_areas(curr_area);
    
    p1 = subplot(4,length(plot_areas),curr_area);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(mua_allcat(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Str ' num2str(area_idx) ' Measured']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p2 = subplot(4,length(plot_areas),curr_area+length(plot_areas)*1);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(fluor_kernel_roi_bw(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['ROI']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p3 = subplot(4,length(plot_areas),curr_area+length(plot_areas)*2);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(fluor_kernel_roi_weighted(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Weighted ROI']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p4 = subplot(4,length(plot_areas),curr_area+length(plot_areas)*3);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(predicted_mua_std_allcat(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Experiment kernel predicted']);
    axis tight
    line([0,0],ylim,'color','k');
    
    linkaxes([p1,p2,p3,p4],'x');
    
end


% Plot residuals
figure;
for curr_area = 1:length(plot_areas)
    area_idx = plot_areas(curr_area);
    
    p1 = subplot(3,length(plot_areas),curr_area);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(mua_allcat(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Str ' num2str(area_idx) ' Measured']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p2 = subplot(3,length(plot_areas),curr_area+length(plot_areas)*1);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(mua_allcat(curr_trials,:,area_idx) - ...
            mua_allcat_predicted(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Str ' num2str(area_idx) ' Task residual']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p3 = subplot(3,length(plot_areas),curr_area+length(plot_areas)*2);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(mua_allcat(curr_trials,:,area_idx) - ...
            predicted_mua_std_allcat(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Str ' num2str(area_idx) ' Ctx residual']);
    axis tight
    line([0,0],ylim,'color','k');
    
    linkaxes([p1,p2,p3],'xy');
    
end


%%%%%% Plot measured and predicted relating to each kernel (move-aligned)

% Re-align activity to movement onset
move_onset_regressor_idx = strcmp('Move onset',regressor_labels);
predicted_mua_std_allcat_move = predicted_mua_std_allcat - predicted_mua_allcat_predicted_reduced(:,:,:,move_onset_regressor_idx);
mua_allcat_move = mua_allcat - mua_allcat_predicted_reduced(:,:,:,move_onset_regressor_idx);
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
for i = 1:size(fluor_allcat_deconv,1)
    predicted_mua_std_allcat_move(i,:,:) = circshift(predicted_mua_std_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
    mua_allcat_move(i,:,:) = circshift(mua_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
end

% Plot measured and predicted relating to each kernel
plot_areas = [1,2,3,4];

plot_conditions = ...
    [1,-1,1; ...
    -1,1,1; ...
    0,-1,1; ...
    0,1,1];
plot_condition_col = ... 
    [1,0,0; ...
    0,0,1; ...
    1,0,1; ...
    0,1,1];
trial_conditions = ...
    [sign(trial_contrastside_allcat), trial_choice_allcat, move_t < 0.5];
[~,trial_id] = ismember(trial_conditions,plot_conditions,'rows');

plot_conditions = {[1:4],[1:4],[1:4],[1:4]};
plot_conditions_compare = { ...
    trial_id, ...  
    trial_id, ...
    trial_id, ...
    trial_id};

stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = 0;
plot_cols = { ...
    plot_condition_col, ...
    plot_condition_col, ...
    plot_condition_col, ...
    plot_condition_col};

% (reward)
% plot_condition_col = ... 
%     [0,0,1; ...
%     0,0,0];
% 
% plot_conditions = {[1,0],[1,0],[1,0],[1,0]};
% plot_conditions_compare = { ...
%     trial_side_allcat == -trial_choice_allcat & move_t < 0.5, ...  
%     trial_side_allcat == -trial_choice_allcat & move_t < 0.5, ...  
%     trial_side_allcat == -trial_choice_allcat & move_t < 0.5, ...  
%     trial_side_allcat == -trial_choice_allcat & move_t < 0.5};
% 
% plot_cols = { ...
%     plot_condition_col, ...
%     plot_condition_col, ...
%     plot_condition_col, ...
%     plot_condition_col};

figure;
for curr_area = 1:length(plot_areas)
    area_idx = plot_areas(curr_area);
    
    p1 = subplot(2,length(plot_areas),curr_area);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(mua_allcat_move(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Str ' num2str(area_idx) ' Measured']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p2 = subplot(2,length(plot_areas),curr_area+length(plot_areas)*1);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(predicted_mua_std_allcat_move(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Experiment kernel predicted']);
    axis tight
    line([0,0],ylim,'color','k');
    
    linkaxes([p1,p2],'xy');
    
end



% %%%% SVM
% 
% % L/R stim SVM on each time point
% % (subtract all but stim)
% disp('Predicting stim side...');
% stim_regressor_idx = strcmp('Stim',regressor_labels);
% predicted_mua_std_allcat_stim = predicted_mua_std_allcat - predicted_mua_allcat_predicted_reduced(:,:,:,stim_regressor_idx);
% mua_allcat_stim = mua_allcat - mua_allcat_predicted_reduced(:,:,:,stim_regressor_idx);
% 
% activity_cat = cat(3,mua_allcat_stim,predicted_mua_std_allcat_stim);
% stim_frac_correct_predicted = nan(size(activity_cat,3),length(t));
% for curr_area = 1:size(activity_cat,3)
%     for curr_t = 1:length(t)        
%         use_data = activity_cat(:,curr_t,curr_area);
%         use_predict = trial_side_allcat;
%         
%         use_trials = ~isnan(use_data) & move_t < 0.5 & trial_contrast_allcat > 0;
%         
%         SVMModel = fitcsvm(use_data(use_trials),use_predict(use_trials));
%         CVSVMModel = crossval(SVMModel);
%         predicted_val = kfoldPredict(CVSVMModel);
%         stim_frac_correct_predicted(curr_area,curr_t) = ...
%             sum(use_predict(use_trials) == predicted_val)/ ...
%             sum(use_trials);           
%         AP_print_progress_fraction(curr_t,length(t));
%     end
%     AP_print_progress_fraction(curr_area,size(activity_cat,3));
% end
% 
% % L/R move SVM on each time point
% % (realign to movement, subtract all but move onset)
% disp('Predicting move direction...');
% move_onset_regressor_idx = strcmp('Move onset',regressor_labels);
% predicted_mua_std_allcat_move = predicted_mua_std_allcat - predicted_mua_allcat_predicted_reduced(:,:,:,move_onset_regressor_idx);
% mua_allcat_move = mua_allcat - mua_allcat_predicted_reduced(:,:,:,move_onset_regressor_idx);
% t_leeway = -t(1);
% leeway_samples = round(t_leeway*(sample_rate));
% for i = 1:size(fluor_allcat_deriv,1)
%     predicted_mua_std_allcat_move(i,:,:) = circshift(predicted_mua_std_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
%     mua_allcat_move(i,:,:) = circshift(mua_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
% end
% 
% activity_cat = cat(3,mua_allcat_move,predicted_mua_std_allcat_move);
% move_frac_correct_predicted = nan(size(activity_cat,3),length(t));
% for curr_area = 1:size(activity_cat,3)
%     for curr_t = 1:length(t)        
%         use_data = activity_cat(:,curr_t,curr_area);
%         use_predict = trial_choice_allcat;
%         
%         use_trials = ~isnan(use_data) & move_t < 0.5;
%         
%         SVMModel = fitcsvm(use_data(use_trials),use_predict(use_trials));
%         CVSVMModel = crossval(SVMModel);
%         predicted_val = kfoldPredict(CVSVMModel);
%         move_frac_correct_predicted(curr_area,curr_t) = ...
%             sum(use_predict(use_trials) == predicted_val)/ ...
%             sum(use_trials);           
%         AP_print_progress_fraction(curr_t,length(t));
%     end
%     AP_print_progress_fraction(curr_area,size(activity_cat,3));
% end
% 
% % Reward/no reward SVM on each time point
% % (subtract all but reward)
% disp('Predicting reward...');
% reward_regressor_idx = strcmp('Reward',regressor_labels);
% predicted_mua_std_allcat_reward = predicted_mua_std_allcat - predicted_mua_allcat_predicted_reduced(:,:,:,reward_regressor_idx);
% mua_allcat_reward = mua_allcat - mua_allcat_predicted_reduced(:,:,:,reward_regressor_idx);
% 
% activity_cat = cat(3,mua_allcat_reward,predicted_mua_std_allcat_reward);
% reward_frac_correct_predicted = nan(size(activity_cat,3),length(t));
% for curr_area = 1:size(activity_cat,3)
%     for curr_t = 1:length(t)        
%         use_data = activity_cat(:,curr_t,curr_area);
%         use_predict = trial_side_allcat == -trial_choice_allcat;
%         
%         use_trials = ~isnan(use_data) & move_t < 0.5;
%         
%         SVMModel = fitcsvm(use_data(use_trials),use_predict(use_trials));
%         CVSVMModel = crossval(SVMModel);
%         predicted_val = kfoldPredict(CVSVMModel);
%         reward_frac_correct_predicted(curr_area,curr_t) = ...
%             sum(use_predict(use_trials) == predicted_val)/ ...
%             sum(use_trials);           
%         AP_print_progress_fraction(curr_t,length(t));
%     end
%     AP_print_progress_fraction(curr_area,size(activity_cat,3));
% end
% 
% % Plot SVM predictions 
% figure; 
% subplot(1,2,1);
% hold on; set(gca,'ColorOrder',[copper(n_depths);cool(n_depths)]);
% plot(t,stim_frac_correct_predicted','linewidth',2);
% line([0,0],ylim,'color','k');
% xlabel('Time from stim')
% ylabel('Fraction correct predicted');
% title('Stim side');
% 
% subplot(1,2,2);
% hold on; set(gca,'ColorOrder',[copper(n_depths);cool(n_depths)]);
% plot(t,move_frac_correct_predicted','linewidth',2);
% line([0,0],ylim,'color','k');
% xlabel('Time from move')
% ylabel('Fraction correct predicted');
% title('Move direction');


%%%% Get activity related to different predictions

%%%%% TESTING
% check the ROIs drawn from the kernels
fluor_kernel_roi_bw = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],kernel_roi.bw), ...
    size(kernel_roi.bw,3),[],size(fluor_allcat_deconv,1)),[3,2,1]);

fluor_kernel_roi_bw_reduced = arrayfun(@(x) ...
    permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_predicted_reduced(:,:,:,x),[3,2,1]), ...
    n_vs,[]),[],[],kernel_roi.bw), ...
    size(kernel_roi.bw,3),[],size(fluor_allcat_predicted,1)),[3,2,1]), ...
    1:size(fluor_allcat_predicted_reduced,4),'uni',false);
fluor_kernel_roi_bw_reduced = cat(4,fluor_kernel_roi_bw_reduced{:});
%%%%%%


trial_groups = {'Stim','Move onset','Reward'};
figure;
for curr_group = 1:length(trial_groups)
    
    % Get MUA/predicted using reduced model
    curr_regressor_idx = strcmp(trial_groups{curr_group},regressor_labels);
    curr_mua =  mua_allcat - mua_allcat_predicted_reduced(:,:,:,curr_regressor_idx);
%     curr_predicted_mua = predicted_mua_std_allcat - predicted_mua_allcat_predicted_reduced(:,:,:,curr_regressor_idx);
    curr_predicted_mua = fluor_kernel_roi_bw - fluor_kernel_roi_bw_reduced(:,:,:,curr_regressor_idx);

    % (if movement, align to move onset)
    if any(strfind(lower(trial_groups{curr_group}),'move'))
        t_leeway = -t(1);
        leeway_samples = round(t_leeway*(sample_rate));
        for i = 1:size(fluor_allcat_deconv,1)
            curr_mua(i,:,:) = circshift(curr_mua(i,:,:),-move_idx(i)+leeway_samples,2);
            curr_predicted_mua(i,:,:) = circshift(curr_predicted_mua(i,:,:),-move_idx(i)+leeway_samples,2);
        end
    end
    
    % Set trials and grouping to use
    switch trial_groups{curr_group}
        case 'Stim'
            use_trials = move_t > 0 & move_t < 1 & trial_contrast_allcat > 0;
            trial_group = trial_side_allcat;
            group1 = 1;
            group2 = -1;
        case 'Move onset'
            use_trials = move_t > 0 & move_t < 1;
            trial_group = trial_choice_allcat;
            group1 = -1;
            group2 = 1;
        case 'Reward'
            use_trials = move_t > 0 & move_t < 1;
            trial_group = trial_choice_allcat == -trial_side_allcat;
            group1 = 1;
            group2 = 0;
    end    
    
    act_rank = tiedrank(curr_mua(use_trials,:,:));
    predicted_act_rank = tiedrank(curr_predicted_mua(use_trials,:,:));
    
    act_rank_difference = squeeze(( ...
        nanmean(act_rank(trial_group(use_trials) == group1,:,:),1) - ...
        nanmean(act_rank(trial_group(use_trials) == group2,:,:),1))./sum(use_trials)); 
    predicted_act_rank_difference = squeeze(( ...
        nanmean(predicted_act_rank(trial_group(use_trials) == group1,:,:),1) - ...
        nanmean(predicted_act_rank(trial_group(use_trials) == group2,:,:),1))./sum(use_trials));    
    
    p1 = subplot(2,length(trial_groups),curr_group); 
    hold on; set(gca,'ColorOrder',copper(n_depths));
    plot(t,act_rank_difference,'linewidth',2);
    axis tight; line([0,0],ylim,'color','k');
    xlabel('Time');
    ylabel('Mean rank difference');
    title([trial_groups{curr_group} ': Measured']);
    
    p2 = subplot(2,length(trial_groups),length(trial_groups)+curr_group); 
    hold on; set(gca,'ColorOrder',copper(n_depths));
    plot(t,predicted_act_rank_difference,'linewidth',2);
    axis tight; line([0,0],ylim,'color','k');
    xlabel('Time');
    ylabel('Mean rank difference');
    title([trial_groups{curr_group} ': Predicted']);
    linkaxes([p1,p2],'xy');        
    
%     %%% (testing shuffle)    
%     n_shuff = 100;
%     rank_diff_shuff = nan(size(act_rank,2),size(act_rank,3),n_shuff);
%     rank_diff_predicted_shuff = nan(size(act_rank,2),size(act_rank,3),n_shuff);
%     for curr_shuff = 1:n_shuff
%         trial_group_shuff = AP_shake(trial_group(use_trials));
%         rank_diff_shuff(:,:,curr_shuff) = squeeze(( ...
%             nanmean(act_rank(trial_group_shuff == group1,:,:),1) - ...
%             nanmean(act_rank(trial_group_shuff == group2,:,:),1))./sum(use_trials));
%         rank_diff_predicted_shuff(:,:,curr_shuff) = squeeze(( ...
%             nanmean(predicted_act_rank(trial_group_shuff == group1,:,:),1) - ...
%             nanmean(predicted_act_rank(trial_group_shuff == group2,:,:),1))./sum(use_trials));
%         AP_print_progress_fraction(curr_shuff,n_shuff);
%     end
%     
%     rank_diff_rank = permute(tiedrank(permute(cat(3,act_rank_difference,rank_diff_shuff),[3,1,2])),[2,3,1]);
%     rank_diff_predicted_rank = permute(tiedrank(permute(cat(3,predicted_act_rank_difference,rank_diff_predicted_shuff),[3,1,2])),[2,3,1]);
%     
%     rank_diff_p = rank_diff_rank(:,:,1)./size(rank_diff_rank,3);
%     rank_diff_predicted_p = rank_diff_predicted_rank(:,:,1)./size(rank_diff_predicted_rank,3);
%     
%     max_plot = max(act_rank_difference(:));
%     spacing = max_plot/10;
%     
%     rank_diff_p_plot = +(rank_diff_p > 0.95);
%     rank_diff_p_plot(~rank_diff_p_plot) = NaN;
%     rank_diff_p_plot = +rank_diff_p_plot.*(1:n_depths)*spacing+max_plot;  
%     plot(p1,t,rank_diff_p_plot,'.','MarkerSize',5);
%     axis tight;
%     
%     rank_diff_predicted_p_plot = +(rank_diff_predicted_p > 0.95);
%     rank_diff_predicted_p_plot(~rank_diff_predicted_p_plot) = NaN;
%     rank_diff_predicted_p_plot = +rank_diff_predicted_p_plot.*(1:n_depths)*spacing+max_plot;   
%     plot(p2,t,rank_diff_predicted_p_plot,'.','MarkerSize',5);
%     axis tight;
    
end

warning('this is fucked up: negative MUA in reduced')
% (predicting stim v no stim)
trial_groups = {'Move L','Move R'};
figure;
for curr_group = 1:length(trial_groups)
    
    % Get MUA/predicted using reduced model
    curr_regressor_idx = strcmp('Move onset',regressor_labels);
    curr_mua =  mua_allcat - mua_allcat_predicted_reduced(:,:,:,curr_regressor_idx);
    curr_predicted_mua = predicted_mua_std_allcat - predicted_mua_allcat_predicted_reduced(:,:,:,curr_regressor_idx);
%     curr_predicted_mua = fluor_kernel_roi_bw - fluor_kernel_roi_bw_reduced(:,:,:,curr_regressor_idx);

    % (if movement, align to move onset)
    if any(strfind(lower(trial_groups{curr_group}),'move'))
        t_leeway = -t(1);
        leeway_samples = round(t_leeway*(sample_rate));
        for i = 1:size(fluor_allcat_deconv,1)
            curr_mua(i,:,:) = circshift(curr_mua(i,:,:),-move_idx(i)+leeway_samples,2);
            curr_predicted_mua(i,:,:) = circshift(curr_predicted_mua(i,:,:),-move_idx(i)+leeway_samples,2);
        end
    end
    
    % Set trials and grouping to use
    switch trial_groups{curr_group}
        case 'Move L'
            use_trials = move_t > 0.1 & move_t < 0.5 & ...
                trial_choice_allcat == -1 & ...
                (trial_contrast_allcat == 0 | trial_side_allcat == 1);
            trial_group = trial_contrast_allcat > 0;
            group1 = 1;
            group2 = 0;
        case 'Move R'
            use_trials = move_t > 0.1 & move_t < 0.5 & ...
                trial_choice_allcat == 1 & ...
                (trial_contrast_allcat == 0 | trial_side_allcat == -1);
            trial_group = trial_contrast_allcat > 0;
            group1 = 1;
            group2 = 0;
    end    
    
    act_rank = tiedrank(curr_mua(use_trials,:,:));    
    predicted_act_rank = tiedrank(curr_predicted_mua(use_trials,:,:));
%     act_rank = curr_mua(use_trials,:,:);    
%     predicted_act_rank = curr_predicted_mua(use_trials,:,:);
    
    act_rank_difference = squeeze(( ...
        nanmean(act_rank(trial_group(use_trials) == group1,:,:),1) - ...
        nanmean(act_rank(trial_group(use_trials) == group2,:,:),1))./sum(use_trials));    
    predicted_act_rank_difference = squeeze(( ...
        nanmean(predicted_act_rank(trial_group(use_trials) == group1,:,:),1) - ...
        nanmean(predicted_act_rank(trial_group(use_trials) == group2,:,:),1))./sum(use_trials));    
    
    p1 = subplot(2,length(trial_groups),curr_group); 
    hold on; set(gca,'ColorOrder',copper(n_depths));
    plot(t,act_rank_difference,'linewidth',2);
    axis tight; line([0,0],ylim,'color','k');
    xlabel('Time');
    ylabel('Mean rank difference');
    title([trial_groups{curr_group} ': Measured']);
    
    p2 = subplot(2,length(trial_groups),length(trial_groups)+curr_group); 
    hold on; set(gca,'ColorOrder',copper(n_depths));
    plot(t,predicted_act_rank_difference,'linewidth',2);
    axis tight; line([0,0],ylim,'color','k');
    xlabel('Time');
    ylabel('Mean rank difference');
    title([trial_groups{curr_group} ': Predicted']);
    linkaxes([p1,p2],'xy');        
    
end


% (same, for fluor)
trial_groups = {'Stim','Move onset','Reward'};
figure;
for curr_group = 1:length(trial_groups)
    
    curr_regressor_idx = strcmp(trial_groups{curr_group},regressor_labels);
    
    fluor_predicted_roi_reduced = permute(reshape( ...
        AP_svd_roi(U_master(:,:,1:n_vs), ...
        reshape(permute(fluor_allcat_predicted_reduced(:,:,:,curr_regressor_idx),[3,2,1]), ...
        n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
        n_rois,[],size(fluor_allcat_predicted,1)),[3,2,1]);
    
    curr_fluor = fluor_roi_deconv - fluor_predicted_roi_reduced;
%     curr_fluor = curr_fluor(:,:,1:size(wf_roi,1)) - curr_fluor(:,:,size(wf_roi,1)+1:end);
    
    % (if movement, align to move onset)
    if strcmp(trial_groups{curr_group},'Move onset')
        t_leeway = -t(1);
        leeway_samples = round(t_leeway*(sample_rate));
        for i = 1:size(fluor_allcat_deconv,1)
            curr_fluor(i,:,:) = circshift(curr_fluor(i,:,:),-move_idx(i)+leeway_samples,2);
        end
    end
    
    % Set trials and grouping (-1,1) to use
    switch trial_groups{curr_group}
        case 'Stim'
            use_trials = move_t < 1 & trial_contrast_allcat > 0;
            trial_group = trial_side_allcat;
            group1 = 1;
            group2 = -1;
        case 'Move onset'
            use_trials = move_t < 1;
            trial_group = trial_choice_allcat;
            group1 = -1;
            group2 = 1;
        case 'Reward'
            use_trials = move_t < 1;
            trial_group = trial_choice_allcat == -trial_side_allcat;
            group1 = 1;
            group2 = 0;
    end
    
    act_rank = tiedrank(curr_fluor(use_trials,:,:));
    act_rank_difference = squeeze(( ...
        nanmean(act_rank(trial_group(use_trials) == group1,:,:),1) - ...
        nanmean(act_rank(trial_group(use_trials) == group2,:,:),1))./sum(use_trials));   
    
    p1 = subplot(1,length(trial_groups),curr_group); 
    hold on; set(gca,'ColorOrder',jet(n_rois));
    plot(t,act_rank_difference,'linewidth',2);
    axis tight; line([0,0],ylim,'color','k');
    xlabel('Time');
    ylabel('Mean rank difference');
    title([trial_groups{curr_group}]);
   
    
    
%     %%% (testing shuffle)    
%     n_shuff = 1000;
%     rank_diff_shuff = nan(size(act_rank,2),size(act_rank,3),n_shuff);
%     for curr_shuff = 1:n_shuff
%         trial_group_shuff = AP_shake(trial_group(use_trials));
%         rank_diff_shuff(:,:,curr_shuff) = squeeze(( ...
%             nanmean(act_rank(trial_group_shuff == group1,:,:),1) - ...
%             nanmean(act_rank(trial_group_shuff == group2,:,:),1))./sum(use_trials));  
%         AP_print_progress_fraction(curr_shuff,n_shuff);
%     end
%     
%     rank_diff_rank = permute(tiedrank(permute(cat(3,act_rank_difference,rank_diff_shuff),[3,1,2])),[2,3,1]);
%     rank_diff_p = rank_diff_rank(:,:,1)./size(rank_diff_rank,3);
%     
%     max_plot = max(act_rank_difference(:));
%     spacing = max_plot/10;
%     
%     rank_diff_p_plot = +(rank_diff_p > 0.95);
%     rank_diff_p_plot(~rank_diff_p_plot) = NaN;
%     rank_diff_p_plot = +rank_diff_p_plot.*(1:n_rois)*spacing+max_plot;  
%     plot(p1,t,rank_diff_p_plot,'.','MarkerSize',5);
%     axis tight;
    
   
    
    
end






% Plot long regression example?
plot_areas = 1:4;
figure;
for curr_area = 1:length(plot_areas)
    plot_area_idx = plot_areas(curr_area);
    long_trace = reshape(permute(mua_allcat(:,:,plot_area_idx),[2,1]),[],1);
    long_trace_ctx_predicted = reshape(permute(predicted_mua_std_allcat(:,:,plot_area_idx),[2,1]),[],1);
    long_trace_task_predicted = reshape(permute(mua_allcat_predicted(:,:,plot_area_idx),[2,1]),[],1);
    
    % Filter long trace
    lowpassCutoff = 6; % Hz
    [b100s, a100s] = butter(2, lowpassCutoff/((sample_rate)/2), 'low');
    long_trace = filter(b100s,a100s,long_trace,[],1);
        
    % (correlate measured and predicted in chunks)
    n_chunks = 1000;
    chunk_boundaries = round(linspace(1,length(long_trace),n_chunks+1));
    chunk_corr = nan(n_chunks,1);
    for curr_chunk = 1:n_chunks
%         curr_corr = corrcoef(long_trace(chunk_boundaries(curr_chunk):chunk_boundaries(curr_chunk+1)), ...
%             long_trace_predicted(chunk_boundaries(curr_chunk):chunk_boundaries(curr_chunk+1)));
%         chunk_corr(curr_chunk) = curr_corr(2);

        curr_corr = sum((long_trace(chunk_boundaries(curr_chunk):chunk_boundaries(curr_chunk+1))- ...
            long_trace_ctx_predicted(chunk_boundaries(curr_chunk):chunk_boundaries(curr_chunk+1))).^2);
        chunk_corr(curr_chunk) = curr_corr;       
    end
    
    % (plot chunk percentiles)
    plot_prctiles = [0,25,50,75,100];
    [chunk_corr_sorted,sort_idx] = sort(chunk_corr);
    chunk_corr_prctile = prctile(chunk_corr,plot_prctiles);
    for curr_prctile = 1:length(chunk_corr_prctile)
        subplot(length(chunk_corr_prctile),length(plot_areas), ...
            (curr_prctile-1)*length(plot_areas)+curr_area); hold on;
        curr_plot_chunk_idx = sort_idx(find(chunk_corr_sorted >= chunk_corr_prctile(curr_prctile),1));
        curr_chunk_plot = chunk_boundaries(curr_plot_chunk_idx):chunk_boundaries(curr_plot_chunk_idx+1);
        plot([1:length(curr_chunk_plot)]/sample_rate,long_trace(curr_chunk_plot),'color',[0,0.6,0]);
        plot([1:length(curr_chunk_plot)]/sample_rate,long_trace_ctx_predicted(curr_chunk_plot),'color',[0.6,0,0.6]);
        plot([1:length(curr_chunk_plot)]/sample_rate,long_trace_task_predicted(curr_chunk_plot),'color',[0.3,0.3,0.8]);
        ylabel([num2str(plot_prctiles(curr_prctile)) 'th percentile']);
        if curr_prctile == 1
            title(plot_areas(curr_area))
        end
    end
end

% Plot all on top of each other 
figure; hold on;
curr_depth = 4;
plot(reshape(mua_allcat(:,:,curr_depth)',[],1));
plot(reshape(predicted_mua_std_allcat(:,:,curr_depth)',[],1));
plot(reshape(mua_allcat_predicted(:,:,curr_depth)',[],1));



%% Fig 4?: stim responses: task > 0.5s reaction

data_fn = ['trial_activity_choiceworld_FLUORTASKTEST'];
exclude_data = true;

AP_load_concat_normalize_ctx_str;

n_vs = size(fluor_allcat_deconv,3);
n_rois = size(fluor_roi_deconv,3);
n_depths = size(mua_allcat,3);

% Get trial information
trial_contrast_allcat = max(D_allcat.stimulus,[],2);
[~,side_idx] = max(D_allcat.stimulus > 0,[],2);
trial_side_allcat = (side_idx-1.5)*2;
trial_choice_allcat = -(D_allcat.response-1.5)*2;

trial_contrast_side = trial_contrast_allcat.*trial_side_allcat;

% Get reaction time
[move_trial,move_idx] = max(abs(wheel_allcat) > 0.02,[],2);
move_t = t(move_idx)';
move_t(~move_trial) = Inf;

% Stim and time to plot
unique_stim = unique(trial_contrast_side);

% Get average long reaction time stim responses
fluor_stim_mean = nan(length(unique_stim),size(fluor_allcat_deconv,2),n_vs);
fluor_roi_stim_mean = nan(length(unique_stim),size(fluor_roi_deconv,2),n_rois);
mua_stim_mean = nan(length(unique_stim),size(mua_allcat,2),n_depths);
mua_ctxpred_stim_mean = nan(length(unique_stim),size(mua_ctxpred_allcat,2),n_depths);
for curr_stim_idx = 1:length(unique_stim)
    curr_trials = move_t > 0.5 & ...
        trial_contrast_side == unique_stim(curr_stim_idx);
    fluor_stim_mean(curr_stim_idx,:,:) = nanmean(fluor_allcat_deconv(curr_trials,:,:),1);
    fluor_roi_stim_mean(curr_stim_idx,:,:) = nanmean(fluor_roi_deconv(curr_trials,:,:),1);
    mua_stim_mean(curr_stim_idx,:,:) = nanmean(mua_allcat(curr_trials,:,:),1);
    mua_ctxpred_stim_mean(curr_stim_idx,:,:) = nanmean(mua_ctxpred_allcat(curr_trials,:,:),1);
end

% Plot 100/R fluorescence, curves, and tuning curves
figure('Name','Trained task > 0.5 reaction time');

use_t = t > 0.05 & t < 0.15;

subplot(2,3,1);
plot_px = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    permute(fluor_stim_mean(end,:,:),[3,2,1]));
plot_px_max = mean(plot_px(:,:,use_t),3);
imagesc(plot_px_max);
caxis([0,max(abs(caxis))]);
colormap(brewermap([],'BuGn'));
AP_reference_outline('ccf_aligned','k');
axis image off;

fp1 = subplot(2,3,2); hold on;
plot_rois = [1:3,7];
plot(t,squeeze(fluor_roi_stim_mean(end,:,plot_rois)),'linewidth',2);
axis tight;
line([0,0],ylim,'color','k');
xlim([-0.2,0.5]);
legend({wf_roi(plot_rois).area});
ylabel('Fluorescence');

mp1 = subplot(2,3,3); hold on; set(gca,'ColorOrder',copper(n_depths));
plot(t,squeeze(mua_stim_mean(end,:,:)),'linewidth',2);
axis tight;
line([0,0],ylim,'color','k');
xlim([-0.2,0.5]);
legend(cellfun(@num2str,num2cell(1:n_depths),'uni',false));
ylabel('MUA');

fp2 = subplot(2,3,4); hold on;
fluor_roi_stim_max = squeeze(max(fluor_roi_stim_mean(:,use_t,:),[],2));
plot(unique_stim,fluor_roi_stim_max(:,plot_rois),'linewidth',2);
xlabel('Contrast*Side');
ylabel('Fluorescence');

mp2 = subplot(2,3,5); hold on; set(gca,'ColorOrder',copper(n_depths));
mua_stim_max = squeeze(max(mua_stim_mean(:,use_t,:),[],2));
plot(unique_stim,mua_stim_max,'linewidth',2);
xlabel('Contrast*Side');
ylabel('MUA');

mp3 = subplot(2,3,6); hold on; set(gca,'ColorOrder',copper(n_depths));
predicted_mua_stim_max = squeeze(max(mua_ctxpred_stim_mean(:,use_t,:),[],2));
plot(unique_stim,predicted_mua_stim_max,'linewidth',2);
xlabel('Contrast*Side');
ylabel('Ctx-predicted MUA');

linkaxes([fp1,fp2],'y');
linkaxes([mp1,mp2,mp3],'y');


% Plot MUA and fluor-predicted MUA
figure;
spacing = max(mua_stim_mean(:));

p1 = subplot(1,3,1); hold on; col = lines(length(unique_stim));
for curr_stim = 1:length(unique_stim)
    AP_stackplot(permute(mua_stim_mean(curr_stim,:,:),[2,3,1]),t, ...
        spacing,false,col(curr_stim,:),1:n_depths,true);
end
title('Measured');

p2 = subplot(1,3,2); hold on; col = lines(length(unique_stim));
for curr_stim = 1:length(unique_stim)
    AP_stackplot(permute(mua_ctxpred_stim_mean(curr_stim,:,:),[2,3,1]),t, ...
        spacing,false,col(curr_stim,:),1:n_depths,true);
end
title('Predicted');

p3 = subplot(1,3,3); hold on; col = lines(length(unique_stim));
for curr_stim = 1:length(unique_stim)
    AP_stackplot(permute(mua_stim_mean(curr_stim,:,:) - ...
        mua_ctxpred_stim_mean(curr_stim,:,:),[2,3,1]),t, ...
        spacing,false,col(curr_stim,:),1:n_depths,true);
end
title('Measured - predicted');

linkaxes([p1,p2,p3]);



% Rank differences by experiment
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_animal;

mua_exp = mat2cell(mua_allcat,use_split,length(t),n_depths);
mua_ctxpred_exp = mat2cell(mua_ctxpred_allcat,use_split,length(t),n_depths);

move_t_exp = mat2cell(move_t,use_split,1);
trial_side_allcat_exp = mat2cell(trial_side_allcat,use_split,1);
trial_contrast_allcat_exp = mat2cell(trial_contrast_allcat,use_split,1);

stim_exp = mat2cell(trial_side_allcat.*trial_contrast_allcat,use_split,1);

trial_groups = {'Stim'};
t_groups = {[0.05,0.15]};

act_rank_difference = nan(length(t),n_depths,length(use_split),length(trial_groups));
predicted_act_rank_difference = nan(length(t),n_depths,length(use_split),length(trial_groups));

act_rank_difference_trial = nan(n_depths,length(use_split),length(trial_groups));
predicted_act_rank_difference_trial = nan(n_depths,length(use_split),length(trial_groups));

for curr_exp = 1:length(mua_exp)
    for curr_group = 1:length(trial_groups)
        
        % Get MUA/predicted using reduced model
        curr_mua =  mua_exp{curr_exp};
        curr_mua_ctx = mua_ctxpred_exp{curr_exp};
        
        % Skip if there's no data in this experiment
        if isempty(curr_mua)
            continue
        end
        
        % Set common NaNs
        nan_samples = isnan(curr_mua) | isnan(curr_mua_ctx);
        curr_mua(nan_samples) = NaN;
        curr_mua_ctx(nan_samples) = NaN;       
        
        % Set trials and grouping to use
        switch trial_groups{curr_group}
            case 'Stim'
                use_trials = move_t_exp{curr_exp} > 0.5 & move_t_exp{curr_exp} < 1 & trial_contrast_allcat_exp{curr_exp} > 0;
                trial_group_1 = trial_side_allcat_exp{curr_exp} == 1;
                trial_group_2 = trial_side_allcat_exp{curr_exp} == -1;
        end
        
        act_rank = tiedrank(curr_mua(use_trials,:,:));
        predicted_act_rank = tiedrank(curr_mua_ctx(use_trials,:,:));
        
        act_rank_difference(:,:,curr_exp,curr_group) = squeeze(( ...
            nanmean(act_rank(trial_group_1(use_trials),:,:),1) - ...
            nanmean(act_rank(trial_group_2(use_trials),:,:),1))./max(act_rank,[],1));
        predicted_act_rank_difference(:,:,curr_exp,curr_group) = squeeze(( ...
            nanmean(predicted_act_rank(trial_group_1(use_trials),:,:),1) - ...
            nanmean(predicted_act_rank(trial_group_2(use_trials),:,:),1))./max(predicted_act_rank,[],1));
        
        use_t = t >= t_groups{curr_group}(1) & t <= t_groups{curr_group}(2);
        act_rank_trial = tiedrank(squeeze(nanmean(curr_mua(use_trials,use_t,:),2)));
        predicted_act_rank_trial = tiedrank(squeeze(nanmean(curr_mua_ctx(use_trials,use_t,:),2)));
        
        act_rank_difference_trial(:,curr_exp,curr_group) = ...
            (nanmean(act_rank_trial(trial_group_1(use_trials),:),1) - ...
            nanmean(act_rank_trial(trial_group_2(use_trials),:),1))./max(act_rank_trial,[],1);
        predicted_act_rank_difference_trial(:,curr_exp,curr_group) = ...
            (nanmean(predicted_act_rank_trial(trial_group_1(use_trials),:),1) - ...
            nanmean(predicted_act_rank_trial(trial_group_2(use_trials),:),1))./max(predicted_act_rank_trial,[],1);
      
    end
end

act_rank_difference_mean = squeeze(nanmean(act_rank_difference,3));
predicted_act_rank_difference_mean = squeeze(nanmean(predicted_act_rank_difference,3));

figure;
for curr_group = 1:length(trial_groups)
    p1 = subplot(2,length(trial_groups),curr_group);
    hold on; set(gca,'ColorOrder',copper(n_depths));
    plot(t,act_rank_difference_mean(:,:,curr_group),'linewidth',2);
    axis tight; line([0,0],ylim,'color','k');
    xlabel('Time');
    ylabel('Mean rank difference');
    title([trial_groups{curr_group} ': Measured']);
    
    p2 = subplot(2,length(trial_groups),length(trial_groups)+curr_group);
    hold on; set(gca,'ColorOrder',copper(n_depths));
    plot(t,predicted_act_rank_difference_mean(:,:,curr_group),'linewidth',2);
    axis tight; line([0,0],ylim,'color','k');
    xlabel('Time');
    ylabel('Mean rank difference');
    title([trial_groups{curr_group} ': Predicted']);
    linkaxes([p1,p2],'xy');
end

figure;

subplot(1,2,1); hold on;
errorbar(squeeze(nanmean(act_rank_difference_trial,2)), ...
    AP_sem(act_rank_difference_trial,2),'linewidth',2);
errorbar(squeeze(nanmean(predicted_act_rank_difference_trial,2)), ...
    AP_sem(predicted_act_rank_difference_trial,2),'linewidth',2);
set(gca,'XTick',1:n_depths);
xlim([0.5,n_depths+0.5]);
line(xlim,[0,0],'color','k');
xlabel('Striatum depth');
ylabel('Stim rank difference');
legend({'Measured','Cortex-predicted'});

% Plot the prediction difference by stimulus
used_exp = cellfun(@(x) ~isempty(x),stim_exp);
mua_ctxpred_diff_exp = ...
    cell2mat(permute(cellfun(@(mua,ctxpred,stim) cell2mat(arrayfun(@(depth) ...
    grpstats(mua(:,:,depth)-ctxpred(:,:,depth),stim,@nanmean), ...
    permute(1:n_depths,[1,3,2]),'uni',false)), ...
    mua_exp(used_exp),mua_ctxpred_exp(used_exp),stim_exp(used_exp),'uni',false),[2,3,4,1]));

t_use = t > 0.05 & t < 0.15;
mua_ctxpred_diff_exp_mean_t = squeeze(nanmean(mua_ctxpred_diff_exp(:,t_use,:,:),2));

subplot(1,2,2); hold on
set(gca,'ColorOrder',copper(n_depths));
errorbar(nanmean(mua_ctxpred_diff_exp_mean_t,3), ...
    AP_sem(mua_ctxpred_diff_exp_mean_t,3),'linewidth',2);
line(xlim,[0,0],'color','k');
xlabel('Stim');
ylabel('Cortex-predicted striatum error');
legend(cellfun(@num2str,num2cell(1:n_depths),'uni',false))





%% Fig 4?: passive stim responses (OLD)

compare_stim_1 = [3:4];
compare_stim_2 = [1:2];

data_fn = 'trial_activity_passive_choiceworld_trained_DECONVTEST';
% data_fn = 'trial_activity_passive_choiceworld_naive_DECONVTEST';
% data_fn = 'trial_activity_passive_fullscreen_trained_DECONVTEST';
% data_fn = 'trial_activity_passive_fullscreen_naive_DECONVTEST';

exclude_data = true;

AP_load_concat_normalize_ctx_str;

n_vs = size(fluor_allcat_deconv,3);
n_rois = size(fluor_roi_deconv,3);
n_depths = size(mua_allcat,3);

% Get stim values
unique_stim = unique(D_allcat.stimulus);

% Get average long reaction time stim responses
fluor_stim_mean = nan(length(unique_stim),size(fluor_allcat_deconv,2),n_vs);
fluor_roi_stim_mean = nan(length(unique_stim),size(fluor_roi_deconv,2),n_rois);
mua_stim_mean = nan(length(unique_stim),size(mua_allcat,2),n_depths);
mua_ctxpred_stim_mean = nan(length(unique_stim),size(mua_ctxpred_allcat,2),n_depths);
for curr_stim_idx = 1:length(unique_stim)
    curr_trials =  ...
        D_allcat.stimulus == unique_stim(curr_stim_idx);
    fluor_stim_mean(curr_stim_idx,:,:) = nanmean(fluor_allcat_deconv(curr_trials,:,:),1);
    fluor_roi_stim_mean(curr_stim_idx,:,:) = nanmean(fluor_roi_deconv(curr_trials,:,:),1);
    mua_stim_mean(curr_stim_idx,:,:) = nanmean(mua_allcat(curr_trials,:,:),1);
    mua_ctxpred_stim_mean(curr_stim_idx,:,:) = nanmean(mua_ctxpred_allcat(curr_trials,:,:),1);
end

% Plot example fluorescence, curves, and tuning curves
use_t = t > 0.05 & t < 0.15;

figure;
subplot(2,3,1);
plot_px = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    permute(fluor_stim_mean(end,:,:),[3,2,1]));
plot_px_max = mean(plot_px(:,:,use_t),3);
imagesc(plot_px_max);
caxis([0,max(abs(caxis))]);
colormap(brewermap([],'BuGn'));
AP_reference_outline('ccf_aligned','k');
axis image off;

fp1 = subplot(2,3,2); hold on;
plot_rois = [1:3,7];
plot(t,squeeze(fluor_roi_stim_mean(end,:,plot_rois)),'linewidth',2);
axis tight;
line([0,0],ylim,'color','k');
xlim([-0.2,0.5]);
legend({wf_roi(plot_rois).area});
ylabel('Fluorescence');

mp1 = subplot(2,3,3); hold on; set(gca,'ColorOrder',copper(n_depths));
plot(t,squeeze(mua_stim_mean(end,:,:)),'linewidth',2);
axis tight;
line([0,0],ylim,'color','k');
xlim([-0.2,0.5]);
legend(cellfun(@num2str,num2cell(1:n_depths),'uni',false));
ylabel('MUA');

fp2 = subplot(2,3,4); hold on;
fluor_roi_stim_max = squeeze(max(fluor_roi_stim_mean(:,use_t,:),[],2));
plot(unique_stim,fluor_roi_stim_max(:,plot_rois),'linewidth',2);
xlabel('Stim');
ylabel('Fluorescence');

mp2 = subplot(2,3,5); hold on; set(gca,'ColorOrder',copper(n_depths));
mua_stim_max = squeeze(max(mua_stim_mean(:,use_t,:),[],2));
plot(unique_stim,mua_stim_max,'linewidth',2);
xlabel('Stim');
ylabel('MUA');

mp3 = subplot(2,3,6); hold on; set(gca,'ColorOrder',copper(n_depths));
predicted_mua_stim_max = squeeze(max(mua_ctxpred_stim_mean(:,use_t,:),[],2));
plot(unique_stim,predicted_mua_stim_max,'linewidth',2);
xlabel('Stim');
ylabel('Ctx-predicted MUA');

linkaxes([fp1,fp2],'y');
linkaxes([mp1,mp2,mp3],'y');


% Plot MUA and fluor-predicted MUA
figure;
spacing = max(mua_stim_mean(:));

p1 = subplot(1,3,1); hold on; col = lines(length(unique_stim));
for curr_stim = 1:length(unique_stim)
    AP_stackplot(permute(mua_stim_mean(curr_stim,:,:),[2,3,1]),t, ...
        spacing,false,col(curr_stim,:),1:n_depths,true);
end
title('Measured');

p2 = subplot(1,3,2); hold on; col = lines(length(unique_stim));
for curr_stim = 1:length(unique_stim)
    AP_stackplot(permute(mua_ctxpred_stim_mean(curr_stim,:,:),[2,3,1]),t, ...
        spacing,false,col(curr_stim,:),1:n_depths,true);
end
title('Predicted');

p3 = subplot(1,3,3); hold on; col = lines(length(unique_stim));
for curr_stim = 1:length(unique_stim)
    AP_stackplot(permute(mua_stim_mean(curr_stim,:,:) - ...
        mua_ctxpred_stim_mean(curr_stim,:,:),[2,3,1]),t, ...
        spacing,false,col(curr_stim,:),1:n_depths,true);
end
title('Measured - predicted');

linkaxes([p1,p2,p3]);


% Rank differences by experiment
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_animal;

mua_exp = mat2cell(mua_allcat,use_split,length(t),n_depths);
mua_ctxpred_exp = mat2cell(mua_ctxpred_allcat,use_split,length(t),n_depths);

stim_exp = mat2cell(D_allcat.stimulus,use_split,1);

trial_groups = {'Stim'};
t_groups = {[0.05,0.15]};

act_rank_difference = nan(length(t),n_depths,length(use_split),length(trial_groups));
predicted_act_rank_difference = nan(length(t),n_depths,length(use_split),length(trial_groups));

act_rank_difference_trial = nan(n_depths,length(use_split),length(trial_groups));
predicted_act_rank_difference_trial = nan(n_depths,length(use_split),length(trial_groups));

for curr_exp = 1:length(mua_exp)
    for curr_group = 1:length(trial_groups)
        
        % Get MUA/predicted using reduced model
        curr_mua =  mua_exp{curr_exp};
        curr_mua_ctx = mua_ctxpred_exp{curr_exp};
        
        % Set common NaNs
        nan_samples = isnan(curr_mua) | isnan(curr_mua_ctx);
        curr_mua(nan_samples) = NaN;
        curr_mua_ctx(nan_samples) = NaN;       
        
        % Set trials and grouping to use
        switch trial_groups{curr_group}
            case 'Stim'
                use_trials = true(size(stim_exp{curr_exp}));                
                trial_group_1 = ismember(stim_exp{curr_exp},compare_stim_1);
                trial_group_2 = ismember(stim_exp{curr_exp},compare_stim_2);  
        end
        
        act_rank = tiedrank(curr_mua(use_trials,:,:));
        predicted_act_rank = tiedrank(curr_mua_ctx(use_trials,:,:));
        
        act_rank_difference(:,:,curr_exp,curr_group) = squeeze(( ...
            nanmean(act_rank(trial_group_1(use_trials),:,:),1) - ...
            nanmean(act_rank(trial_group_2(use_trials),:,:),1))./max(act_rank,[],1));
        predicted_act_rank_difference(:,:,curr_exp,curr_group) = squeeze(( ...
            nanmean(predicted_act_rank(trial_group_1(use_trials),:,:),1) - ...
            nanmean(predicted_act_rank(trial_group_2(use_trials),:,:),1))./max(predicted_act_rank,[],1));
        
        use_t = t >= t_groups{curr_group}(1) & t <= t_groups{curr_group}(2);
        act_rank_trial = tiedrank(squeeze(nanmean(curr_mua(use_trials,use_t,:),2)));
        predicted_act_rank_trial = tiedrank(squeeze(nanmean(curr_mua_ctx(use_trials,use_t,:),2)));
        
        act_rank_difference_trial(:,curr_exp,curr_group) = ...
            (nanmean(act_rank_trial(trial_group_1(use_trials),:),1) - ...
            nanmean(act_rank_trial(trial_group_2(use_trials),:),1))./max(act_rank_trial,[],1);
        predicted_act_rank_difference_trial(:,curr_exp,curr_group) = ...
            (nanmean(predicted_act_rank_trial(trial_group_1(use_trials),:),1) - ...
            nanmean(predicted_act_rank_trial(trial_group_2(use_trials),:),1))./max(predicted_act_rank_trial,[],1);
      
    end
end

act_rank_difference_mean = squeeze(nanmean(act_rank_difference,3));
predicted_act_rank_difference_mean = squeeze(nanmean(predicted_act_rank_difference,3));

figure;
for curr_group = 1:length(trial_groups)
    p1 = subplot(2,length(trial_groups),curr_group);
    hold on; set(gca,'ColorOrder',copper(n_depths));
    plot(t,act_rank_difference_mean(:,:,curr_group),'linewidth',2);
    axis tight; line([0,0],ylim,'color','k');
    xlabel('Time');
    ylabel('Mean rank difference');
    title([trial_groups{curr_group} ': Measured']);
    
    p2 = subplot(2,length(trial_groups),length(trial_groups)+curr_group);
    hold on; set(gca,'ColorOrder',copper(n_depths));
    plot(t,predicted_act_rank_difference_mean(:,:,curr_group),'linewidth',2);
    axis tight; line([0,0],ylim,'color','k');
    xlabel('Time');
    ylabel('Mean rank difference');
    title([trial_groups{curr_group} ': Predicted']);
    linkaxes([p1,p2],'xy');
end

figure; 

subplot(1,2,1); hold on;
errorbar(squeeze(nanmean(act_rank_difference_trial,2)), ...
    AP_sem(act_rank_difference_trial,2),'linewidth',2);
errorbar(squeeze(nanmean(predicted_act_rank_difference_trial,2)), ...
    AP_sem(predicted_act_rank_difference_trial,2),'linewidth',2);
set(gca,'XTick',1:n_depths);
xlim([0.5,n_depths+0.5]);
line(xlim,[0,0],'color','k');
xlabel('Striatum depth');
ylabel('Stim rank difference');
legend({'Measured','Cortex-predicted'});

% Plot the prediction difference by stimulus
mua_ctxpred_diff_exp = ...
    cell2mat(permute(cellfun(@(mua,ctxpred,stim) cell2mat(arrayfun(@(depth) ...
    grpstats(mua(:,:,depth)-ctxpred(:,:,depth),stim,@nanmean), ...
    permute(1:n_depths,[1,3,2]),'uni',false)), ...
    mua_exp,mua_ctxpred_exp,stim_exp,'uni',false),[2,3,4,1]));

t_use = t > 0 & t < 0.2;
mua_ctxpred_diff_exp_mean_t = squeeze(nanmean(mua_ctxpred_diff_exp(:,t_use,:,:),2));

subplot(1,2,2); hold on
set(gca,'ColorOrder',copper(n_depths));
errorbar(nanmean(mua_ctxpred_diff_exp_mean_t,3), ...
    AP_sem(mua_ctxpred_diff_exp_mean_t,3),'linewidth',2);
line(xlim,[0,0],'color','k');
xlabel('Stim');
ylabel('Cortex-predicted striatum error');
legend(cellfun(@num2str,num2cell(1:n_depths),'uni',false))


%% Passive choiceworld: combined

data_fns = { ...
    'trial_activity_choiceworld_FLUORTASKTEST', ...
    'trial_activity_passive_choiceworld_trained_DECONVTEST', ...
    'trial_activity_passive_choiceworld_naive_DECONVTEST'};

n_t = 88;
n_depths = 4;
n_animals = 6;

act_rank_difference = nan(n_t,n_depths,n_animals,length(data_fns));
predicted_act_rank_difference = nan(n_t,n_depths,n_animals,length(data_fns));
act_rank_difference_trial = nan(n_depths,n_animals,length(data_fns));
predicted_act_rank_difference_trial = nan(n_depths,n_animals,length(data_fns));

n_shuff = 10000;
act_rank_difference_trial_shuff = nan(n_depths,n_animals,length(data_fns),n_shuff);
predicted_act_rank_difference_trial_shuff = nan(n_depths,n_animals,length(data_fns),n_shuff);
act_rank_difference_trial_predshuff = nan(n_depths,n_animals,length(data_fns),n_shuff);

for curr_group = 1:length(data_fns)    
    
    clearvars -except data_fns curr_group ...
        act_rank_difference ...
        predicted_act_rank_difference ...
        act_rank_difference_trial ...
        predicted_act_rank_difference_trial ...
        n_shuff ...
        act_rank_difference_trial_shuff ...
        predicted_act_rank_difference_trial_shuff ...
        act_rank_difference_trial_predshuff
    
    data_fn = data_fns{curr_group};
    exclude_data = true;
    AP_load_concat_normalize_ctx_str;
    
    switch curr_group
        case 1 % trained behavior
            % Get trial information
            stim = D_allcat.stimulus(:,2) - D_allcat.stimulus(:,1);
            unique_stim = unique(stim);
            
            compare_stim_1 = unique_stim(unique_stim > 0);
            compare_stim_2 = unique_stim(unique_stim < 0);
            
%             compare_stim_1 = [1];
%             compare_stim_2 = [-1];
            
%             compare_stim_1 = [0.125];
%             compare_stim_2 = [-0.125];

        case 2 % trained passive
            stim = D_allcat.stimulus;
            
            compare_stim_1 = [3:4];
            compare_stim_2 = [1:2];

%             compare_stim_1 = [4];
%             compare_stim_2 = [1];
            
        case 3 % naive passive
            stim = D_allcat.stimulus;
            
            compare_stim_1 = [6:10];
            compare_stim_2 = [1:5];
            
%             compare_stim_1 = [7,10];
%             compare_stim_2 = [1,4];
            
%             compare_stim_1 = [10];
%             compare_stim_2 = [1];
            
    end
           
    % Get trials with movment during stimulus to eliminate
    move_trial = any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > 0.02,2);
    
    % Split data
    trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
    trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
    use_split = trials_animal;
    
    mua_allcat_exp = mat2cell(mua_allcat,use_split,length(t),n_depths);
    mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,use_split,length(t),n_depths);
    
    stim_exp = mat2cell(stim,use_split,1);
    move_trial_exp = mat2cell(move_trial,use_split,1);
    
    % Get rank differences
    t_stim = [0.05,0.15];

    for curr_exp = 1:length(mua_allcat_exp)
        
        % Get MUA/predicted using reduced model
        curr_mua =  mua_allcat_exp{curr_exp};
        curr_mua_ctx = mua_ctxpred_allcat_exp{curr_exp};
        
        % Skip if no data
        if isempty(curr_mua(:))
            continue
        end
        
        % Set common NaNs
        nan_samples = isnan(curr_mua) | isnan(curr_mua_ctx);
        curr_mua(nan_samples) = NaN;
        curr_mua_ctx(nan_samples) = NaN;
        
        % Set trials and grouping to use
        use_trials = ~move_trial_exp{curr_exp};
        trial_group_1 = ismember(stim_exp{curr_exp},compare_stim_1);
        trial_group_2 = ismember(stim_exp{curr_exp},compare_stim_2);
        
        act_rank = tiedrank(curr_mua(use_trials,:,:));
        predicted_act_rank = tiedrank(curr_mua_ctx(use_trials,:,:));
        
        act_rank_difference(:,:,curr_exp,curr_group) = squeeze(( ...
            nanmean(act_rank(trial_group_1(use_trials),:,:),1) - ...
            nanmean(act_rank(trial_group_2(use_trials),:,:),1))./max(act_rank,[],1));
        predicted_act_rank_difference(:,:,curr_exp,curr_group) = squeeze(( ...
            nanmean(predicted_act_rank(trial_group_1(use_trials),:,:),1) - ...
            nanmean(predicted_act_rank(trial_group_2(use_trials),:,:),1))./max(predicted_act_rank,[],1));
        
        use_t = t >= t_stim(1) & t <= t_stim(2);
        act_rank_trial = tiedrank(squeeze(nanmean(curr_mua(use_trials,use_t,:),2)));
        predicted_act_rank_trial = tiedrank(squeeze(nanmean(curr_mua_ctx(use_trials,use_t,:),2)));
        
        act_rank_difference_trial(:,curr_exp,curr_group) = ...
            (nanmean(act_rank_trial(trial_group_1(use_trials),:),1) - ...
            nanmean(act_rank_trial(trial_group_2(use_trials),:),1))./max(act_rank_trial,[],1);
        predicted_act_rank_difference_trial(:,curr_exp,curr_group) = ...
            (nanmean(predicted_act_rank_trial(trial_group_1(use_trials),:),1) - ...
            nanmean(predicted_act_rank_trial(trial_group_2(use_trials),:),1))./max(predicted_act_rank_trial,[],1);

        % Shuffle for significance
        shuff_trials = (trial_group_1 | trial_group_2) & use_trials;
        
        for curr_shuff = 1:n_shuff
            % Shuffle for label significance
            trial_group_1_shuff = trial_group_1;
            trial_group_1_shuff(shuff_trials) = AP_shake(trial_group_1_shuff(shuff_trials));
            trial_group_2_shuff = trial_group_2;
            trial_group_2_shuff(shuff_trials) = AP_shake(trial_group_2_shuff(shuff_trials));
            
            act_rank_difference_trial_shuff(:,curr_exp,curr_group,curr_shuff) = ...
                (nanmean(act_rank_trial(trial_group_1_shuff(use_trials),:),1) - ...
                nanmean(act_rank_trial(trial_group_2_shuff(use_trials),:),1))./max(act_rank_trial,[],1);
            predicted_act_rank_difference_trial_shuff(:,curr_exp,curr_group,curr_shuff) = ...
                (nanmean(predicted_act_rank_trial(trial_group_1_shuff(use_trials),:),1) - ...
                nanmean(predicted_act_rank_trial(trial_group_2_shuff(use_trials),:),1))./max(predicted_act_rank_trial,[],1);
        end
        
        % Shuffle for measured/predicted difference
        % (build an n_shuff sized matrix of half/half, then shake)
        meas_pred_shuff = AP_shake(cat(3, ...
            repmat(act_rank_trial,1,1,n_shuff), ...
            repmat(predicted_act_rank_trial,1,1,n_shuff)),3);
        
        act_rank_difference_trial_predshuff(:,curr_exp,curr_group,:) = ...
            permute((nanmean(meas_pred_shuff(trial_group_1(use_trials),:,1:n_shuff),1) - ...
            nanmean(meas_pred_shuff(trial_group_2(use_trials),:,1:n_shuff),1))./max(act_rank_trial,[],1) - ...
             (nanmean(meas_pred_shuff(trial_group_1(use_trials),:,n_shuff+1:end),1) - ...
            nanmean(meas_pred_shuff(trial_group_2(use_trials),:,n_shuff+1:end),1))./max(act_rank_trial,[],1),[2,1,3,4]);
        
    end
  
    AP_print_progress_fraction(curr_group,length(data_fns));
    
end


figure; 
p = nan(3,length(data_fn));
for curr_group = 1:length(data_fns)
    
    p(1,curr_group) = subplot(3,length(data_fns),curr_group+length(data_fns)*0); hold on;
%     set(gca,'ColorOrder',copper(n_depths));
%     plot(t,nanmean(act_rank_difference(:,:,:,curr_data_group),3),'linewidth',2);
%     
    col = copper(n_depths);
    for curr_depth = 1:n_depths
        AP_errorfill(t,nanmean(act_rank_difference(:,curr_depth,:,curr_group),3), ...
            AP_sem(act_rank_difference(:,curr_depth,:,curr_group),3),col(curr_depth,:));
    end
    xlabel('Time from stim');
    ylabel('Rank difference');
    
    p(2,curr_group) = subplot(3,length(data_fns),curr_group+length(data_fns)*1); hold on;
%     set(gca,'ColorOrder',copper(n_depths));
%     plot(t,nanmean(predicted_act_rank_difference(:,:,:,curr_data_group),3),'linewidth',2);

    col = copper(n_depths);
    for curr_depth = 1:n_depths
        AP_errorfill(t,nanmean(predicted_act_rank_difference(:,curr_depth,:,curr_group),3), ...
            AP_sem(predicted_act_rank_difference(:,curr_depth,:,curr_group),3),col(curr_depth,:));
    end
    xlabel('Time from stim');
    ylabel('Rank difference');
    
    p(3,curr_group) = subplot(3,length(data_fns),curr_group+length(data_fns)*2); hold on;
    
    errorbar(squeeze(nanmean(act_rank_difference_trial(:,:,curr_group),2)), ...
        AP_sem(act_rank_difference_trial(:,:,curr_group),2),'linewidth',2,'color','k');
    errorbar(squeeze(nanmean(predicted_act_rank_difference_trial(:,:,curr_group),2)), ...
        AP_sem(predicted_act_rank_difference_trial(:,:,curr_group),2),'linewidth',2,'color',[0,0.7,0]);
    set(gca,'XTick',1:n_depths);
    xlim([0.5,n_depths+0.5]);
    line(xlim,[0,0],'color','k');
    xlabel('Striatum depth');
    ylabel('Stim rank difference');
    legend({'Measured','Cortex-predicted'});
    title(data_fn);
    drawnow;
    
end

linkaxes(p(1:2,:),'xy');
linkaxes(p(3,:),'xy');

% Plot significance
figure; 

% Cross-group significance
bhv_passive_diff_ci = prctile(squeeze(nanmean(AP_shake(cat(3, ...
    repmat(act_rank_difference_trial(:,:,1),1,1,n_shuff/2), ...
    repmat(act_rank_difference_trial(:,:,3),1,1,n_shuff/2)),3) - ...
    AP_shake(cat(3, ...
    repmat(act_rank_difference_trial(:,:,1),1,1,n_shuff/2), ...
    repmat(act_rank_difference_trial(:,:,3),1,1,n_shuff/2)),3),2)),[2.5,97.5],2);

subplot(1,length(data_fns)+2,1); hold on;
plot(nanmean(act_rank_difference_trial(:,:,1) - act_rank_difference_trial(:,:,3),2),'r','linewidth',2);
plot(bhv_passive_diff_ci,'k','linewidth',2,'linestyle','--');
title('Trained-naive');

predicted_bhv_passive_diff_ci = prctile(squeeze(nanmean(AP_shake(cat(3, ...
    repmat(predicted_act_rank_difference_trial(:,:,1),1,1,n_shuff/2), ...
    repmat(predicted_act_rank_difference_trial(:,:,3),1,1,n_shuff/2)),3) - ...
    AP_shake(cat(3, ...
    repmat(predicted_act_rank_difference_trial(:,:,1),1,1,n_shuff/2), ...
    repmat(predicted_act_rank_difference_trial(:,:,3),1,1,n_shuff/2)),3),2)),[2.5,97.5],2);

subplot(1,length(data_fns)+2,2); hold on;
plot(nanmean(predicted_act_rank_difference_trial(:,:,1) - predicted_act_rank_difference_trial(:,:,3),2),'r','linewidth',2);
plot(bhv_passive_diff_ci,'k','linewidth',2,'linestyle','--');
title('Trained-naive (predicted)');

% Significance of measured vs predicted
meas_pred_diff_ci = prctile(squeeze(nanmean(act_rank_difference_trial_predshuff,2)),[2.5,97.5],3);
for curr_group = 1:length(data_fns)
    subplot(1,length(data_fns)+2,2+curr_group); hold on;
    plot(squeeze(nanmean(act_rank_difference_trial(:,:,curr_group),2)) - ...
        squeeze(nanmean(predicted_act_rank_difference_trial(:,:,1),2)),'r','linewidth',2);
    plot(reshape(permute(meas_pred_diff_ci(:,curr_group,:),[1,3,2]),n_depths,[]),'k','linewidth',2,'linestyle','--');
    xlabel('Striatum depth');
    ylabel('Measured - Predicted');
    title(data_fns{curr_group});
end


%% Passive full screen: combined

data_fns = { ...
    'trial_activity_passive_fullscreen_trained_DECONVTEST', ...
    'trial_activity_passive_fullscreen_naive_DECONVTEST'};

n_t = 88;
n_depths = 4;
n_animals = 6;

act_rank_difference = nan(n_t,n_depths,n_animals,length(data_fns));
predicted_act_rank_difference = nan(n_t,n_depths,n_animals,length(data_fns));
act_rank_difference_trial = nan(n_depths,n_animals,length(data_fns));
predicted_act_rank_difference_trial = nan(n_depths,n_animals,length(data_fns));

n_shuff = 1000;
act_rank_difference_trial_stimshuff = nan(n_depths,n_animals,length(data_fns),n_shuff);
predicted_act_rank_difference_trial_stimshuff = nan(n_depths,n_animals,length(data_fns),n_shuff);
act_rank_difference_trial_predshuff = nan(n_depths,n_animals,length(data_fns),n_shuff);

for curr_group = 1:length(data_fns)    
    
    clearvars -except data_fns curr_group ...
        act_rank_difference ...
        predicted_act_rank_difference ...
        act_rank_difference_trial ...
        predicted_act_rank_difference_trial ...
        n_shuff ...
        act_rank_difference_trial_stimshuff ...
        predicted_act_rank_difference_trial_stimshuff ...
        act_rank_difference_trial_predshuff
    
    data_fn = data_fns{curr_group};
    exclude_data = true;
    AP_load_concat_normalize_ctx_str;
       
    stim = D_allcat.stimulus;
    compare_stim_1 = 3;
    compare_stim_2 = 2;
        
    % Get trials with movment during stimulus to eliminate
    move_trial = any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > 0.02,2);
    
    % Split data
    trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
    trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
    use_split = trials_animal;
    
    mua_allcat_exp = mat2cell(mua_allcat,use_split,length(t),n_depths);
    mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,use_split,length(t),n_depths);
    
    stim_exp = mat2cell(stim,use_split,1);
    move_trial_exp = mat2cell(move_trial,use_split,1);
    
    % Get rank differences
    t_stim = [0.05,0.15];

    for curr_exp = 1:length(mua_allcat_exp)
        
        % Get MUA/predicted using reduced model
        curr_mua =  mua_allcat_exp{curr_exp};
        curr_mua_ctx = mua_ctxpred_allcat_exp{curr_exp};
        
        % Skip if no data
        if isempty(curr_mua(:))
            continue
        end
        
        % Set common NaNs
        nan_samples = isnan(curr_mua) | isnan(curr_mua_ctx);
        curr_mua(nan_samples) = NaN;
        curr_mua_ctx(nan_samples) = NaN;
        
        % Set trials and grouping to use
        use_trials = ~move_trial_exp{curr_exp};
        trial_group_1 = ismember(stim_exp{curr_exp},compare_stim_1);
        trial_group_2 = ismember(stim_exp{curr_exp},compare_stim_2);
        
        act_rank = tiedrank(curr_mua(use_trials,:,:));
        predicted_act_rank = tiedrank(curr_mua_ctx(use_trials,:,:));
        
        act_rank_difference(:,:,curr_exp,curr_group) = squeeze(( ...
            nanmean(act_rank(trial_group_1(use_trials),:,:),1) - ...
            nanmean(act_rank(trial_group_2(use_trials),:,:),1))./max(act_rank,[],1));
        predicted_act_rank_difference(:,:,curr_exp,curr_group) = squeeze(( ...
            nanmean(predicted_act_rank(trial_group_1(use_trials),:,:),1) - ...
            nanmean(predicted_act_rank(trial_group_2(use_trials),:,:),1))./max(predicted_act_rank,[],1));
        
        use_t = t >= t_stim(1) & t <= t_stim(2);
        act_rank_trial = tiedrank(squeeze(nanmean(curr_mua(use_trials,use_t,:),2)));
        predicted_act_rank_trial = tiedrank(squeeze(nanmean(curr_mua_ctx(use_trials,use_t,:),2)));
        
        act_rank_difference_trial(:,curr_exp,curr_group) = ...
            (nanmean(act_rank_trial(trial_group_1(use_trials),:),1) - ...
            nanmean(act_rank_trial(trial_group_2(use_trials),:),1))./max(act_rank_trial,[],1);
        predicted_act_rank_difference_trial(:,curr_exp,curr_group) = ...
            (nanmean(predicted_act_rank_trial(trial_group_1(use_trials),:),1) - ...
            nanmean(predicted_act_rank_trial(trial_group_2(use_trials),:),1))./max(predicted_act_rank_trial,[],1);
        
        % Shuffle for label significance
        shuff_trials = (trial_group_1 | trial_group_2) & use_trials;
       
        for curr_shuff = 1:n_shuff
            trial_group_1_shuff = trial_group_1;
            trial_group_1_shuff(shuff_trials) = AP_shake(trial_group_1(shuff_trials));
            trial_group_2_shuff = trial_group_2;
            trial_group_2_shuff(shuff_trials) = AP_shake(trial_group_2(shuff_trials));
            
            act_rank_difference_trial_stimshuff(:,curr_exp,curr_group,curr_shuff) = ...
                (nanmean(act_rank_trial(trial_group_1_shuff(use_trials),:),1) - ...
                nanmean(act_rank_trial(trial_group_2_shuff(use_trials),:),1))./max(act_rank_trial,[],1);
            predicted_act_rank_difference_trial_stimshuff(:,curr_exp,curr_group,curr_shuff) = ...
                (nanmean(predicted_act_rank_trial(trial_group_1_shuff(use_trials),:),1) - ...
                nanmean(predicted_act_rank_trial(trial_group_2_shuff(use_trials),:),1))./max(predicted_act_rank_trial,[],1);           
        end
        
        % Shuffle for measured/predicted difference
        % (build an n_shuff sized matrix of half/half, then shake)
        meas_pred_shuff = AP_shake(cat(3, ...
            repmat(act_rank_trial,1,1,n_shuff), ...
            repmat(predicted_act_rank_trial,1,1,n_shuff)),3);
        
        act_rank_difference_trial_predshuff(:,curr_exp,curr_group,:) = ...
            permute((nanmean(meas_pred_shuff(trial_group_1(use_trials),:,1:n_shuff),1) - ...
            nanmean(meas_pred_shuff(trial_group_2(use_trials),:,1:n_shuff),1))./max(act_rank_trial,[],1) - ...
            (nanmean(meas_pred_shuff(trial_group_1(use_trials),:,n_shuff+1:end),1) - ...
            nanmean(meas_pred_shuff(trial_group_2(use_trials),:,n_shuff+1:end),1))./max(act_rank_trial,[],1),[2,1,3,4]);
        
    end
  
    AP_print_progress_fraction(curr_group,length(data_fns));
    
end


figure; 
p = nan(3,length(data_fn));
for curr_group = 1:length(data_fns)
    
    p(1,curr_group) = subplot(3,length(data_fns),curr_group+length(data_fns)*0); hold on;
%     set(gca,'ColorOrder',copper(n_depths));
%     plot(t,nanmean(act_rank_difference(:,:,:,curr_data_group),3),'linewidth',2);
%     
    col = copper(n_depths);
    for curr_depth = 1:n_depths
        AP_errorfill(t,nanmean(act_rank_difference(:,curr_depth,:,curr_group),3), ...
            AP_sem(act_rank_difference(:,curr_depth,:,curr_group),3),col(curr_depth,:));
    end
    xlabel('Time from stim');
    ylabel('Rank difference');
    
    p(2,curr_group) = subplot(3,length(data_fns),curr_group+length(data_fns)*1); hold on;
%     set(gca,'ColorOrder',copper(n_depths));
%     plot(t,nanmean(predicted_act_rank_difference(:,:,:,curr_data_group),3),'linewidth',2);

    col = copper(n_depths);
    for curr_depth = 1:n_depths
        AP_errorfill(t,nanmean(predicted_act_rank_difference(:,curr_depth,:,curr_group),3), ...
            AP_sem(predicted_act_rank_difference(:,curr_depth,:,curr_group),3),col(curr_depth,:));
    end
    xlabel('Time from stim');
    ylabel('Rank difference');
    
    p(3,curr_group) = subplot(3,length(data_fns),curr_group+length(data_fns)*2); hold on;
    
    errorbar(squeeze(nanmean(act_rank_difference_trial(:,:,curr_group),2)), ...
        AP_sem(act_rank_difference_trial(:,:,curr_group),2),'linewidth',2);
    errorbar(squeeze(nanmean(predicted_act_rank_difference_trial(:,:,curr_group),2)), ...
        AP_sem(predicted_act_rank_difference_trial(:,:,curr_group),2),'linewidth',2);
    set(gca,'XTick',1:n_depths);
    xlim([0.5,n_depths+0.5]);
    line(xlim,[0,0],'color','k');
    xlabel('Striatum depth');
    ylabel('Stim rank difference');
    legend({'Measured','Cortex-predicted'});
    title(data_fn);
    drawnow;
    
end

linkaxes(p(1:2,:),'xy');
linkaxes(p(3,:),'xy');

% Plot significance
figure; 

% Cross-group significance
bhv_passive_diff_ci = prctile(squeeze(nanmean(AP_shake(cat(3, ...
    repmat(act_rank_difference_trial(:,:,1),1,1,n_shuff/2), ...
    repmat(act_rank_difference_trial(:,:,2),1,1,n_shuff/2)),3) - ...
    AP_shake(cat(3, ...
    repmat(act_rank_difference_trial(:,:,1),1,1,n_shuff/2), ...
    repmat(act_rank_difference_trial(:,:,2),1,1,n_shuff/2)),3),2)),[2.5,97.5],2);

subplot(1,length(data_fns)+1,1); hold on;
plot(nanmean(act_rank_difference_trial(:,:,1) - act_rank_difference_trial(:,:,2),2),'r','linewidth',2);
plot(bhv_passive_diff_ci,'k','linewidth',2,'linestyle','--');
title('Trained-naive');

% Significance of measured vs predicted
meas_pred_diff_ci = prctile(squeeze(nanmean(act_rank_difference_trial_predshuff,2)),[2.5,97.5],3);
for curr_group = 1:length(data_fns)
    subplot(1,length(data_fns)+1,1+curr_group); hold on;
    plot(squeeze(nanmean(act_rank_difference_trial(:,:,curr_group),2)) - ...
        squeeze(nanmean(predicted_act_rank_difference_trial(:,:,1),2)),'r','linewidth',2);
    plot(reshape(permute(meas_pred_diff_ci(:,curr_group,:),[1,3,2]),n_depths,[]),'k','linewidth',2,'linestyle','--');
    xlabel('Striatum depth');
    ylabel('Measured - Predicted');
    title(data_fns{curr_group});
end



%% Fig 4?: stim/no-stim move?

% Load data
data_fn = ['trial_activity_choiceworld_framerate'];
exclude_data = true;
AP_load_concat_normalize_ctx_str

% Load regression results
regression_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\trial_activity_choiceworld_task_regression';
load(regression_fn);

n_vs = size(fluor_allcat_deconv,3);
n_rois = size(fluor_roi_deconv,3);
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

sample_rate = 1/mean(diff(t));

% Get reaction time
[move_trial,move_idx] = max(abs(wheel_allcat) > 0.02,[],2);
move_t = t(move_idx)';
move_t(~move_trial) = Inf;

% Get maximum wheel velocity in chosen direction
wheel_velocity_allcat = [zeros(size(wheel_allcat,1),1,1),diff(wheel_allcat,[],2)];
[max_speed,max_vel_idx] = max(abs(wheel_velocity_allcat(:,:,1).* ...
    (bsxfun(@times,wheel_velocity_allcat,trial_choice_allcat) > 0)),[],2);
max_vel = max_speed.*trial_choice_allcat;




% Plot movement onset regressors across the cortex and striatum?
move_onset_regressor_idx = strcmp('Move onset',regressor_labels);
regressor_t = cellfun(@(x) (round(x(1)*(sample_rate)): ...
    round(x(2)*(sample_rate)))/sample_rate,t_shifts,'uni',false);

ctx_v_move_onset_kernel = permute(cell2mat(permute( ...
    fluor_kernel(move_onset_regressor_idx,1:n_vs),[1,3,2])),[3,2,1]);
ctx_px_move_onset_kernel = reshape(svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    reshape(ctx_v_move_onset_kernel,n_vs,[])),size(U_master,1),size(U_master,2),[],2);
ctx_roi_move_onset_kernel = reshape(AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(ctx_v_move_onset_kernel,n_vs,[]),[],[],cat(3,wf_roi.mask)),[], ...
    size(ctx_v_move_onset_kernel,2),size(ctx_v_move_onset_kernel,3));

str_move_onset_kernel = ...
    permute(cell2mat(permute(mua_kernel(move_onset_regressor_idx,:),[1,3,2])),[3,2,1]);

figure;
p1 = subplot(2,2,1); hold on; set(gca,'ColorOrder',jet(n_rois));
plot(regressor_t{2},ctx_roi_move_onset_kernel(:,:,2)','linewidth',2)
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');
legend({wf_roi.area});
title('Move onset L kernel');

p2 = subplot(2,2,3); hold on; set(gca,'ColorOrder',jet(n_rois));
plot(regressor_t{2},-diff(ctx_roi_move_onset_kernel,[],3)','linewidth',2)
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');
title('Move onset L-R kernel');

p3 = subplot(2,2,2); hold on; set(gca,'ColorOrder',copper(n_depths));
plot(regressor_t{2},str_move_onset_kernel(:,:,2)','linewidth',2)
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');
title('Move onset L kernel');

p4 = subplot(2,2,4); hold on; set(gca,'ColorOrder',copper(n_depths));
plot(regressor_t{2},-diff(str_move_onset_kernel,[],3)','linewidth',2)
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');
title('Move onset L-R kernel');

linkaxes([p1,p2],'xy');
linkaxes([p3,p4],'xy');

% Get predicted fluorescence in ROIs from reduced model w/o move onset
move_onset_regressor_idx = strcmp('Move onset',regressor_labels);
fluor_predicted_roi_reduced = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_predicted_reduced(:,:,:,move_onset_regressor_idx),[3,2,1]), ...
    n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
    n_rois,[],size(fluor_allcat_predicted,1)),[3,2,1]);

% Re-align activity to movement onset
fluor_predicted_roi_reduced_move = fluor_roi_deconv - fluor_predicted_roi_reduced;
mua_allcat_move = mua_allcat-mua_taskpred_reduced_allcat(:,:,:,2);
wheel_velocity_allcat_move = wheel_velocity_allcat;
t_leeway = 0.5;
leeway_samples = round(t_leeway*(sample_rate));
for i = 1:size(fluor_allcat_deconv,1)
    fluor_predicted_roi_reduced_move(i,:,:) = circshift(fluor_predicted_roi_reduced_move(i,:,:),-move_idx(i)+leeway_samples,2);
    mua_allcat_move(i,:,:) = circshift(mua_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
    wheel_velocity_allcat_move(i,:,:) = circshift(wheel_velocity_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);    
end

% Bin maximum velocity
use_rxn = move_t > 0 & move_t < 1;

vis_trials = use_rxn & ...
    trial_side_allcat == -trial_choice_allcat & trial_contrast_allcat > 0;
zero_trials = use_rxn & trial_contrast_allcat == 0;

n_prctiles = 6;
vel_amp_edges = prctile(abs(max_vel),linspace(0,100,n_prctiles));
vel_amp_edges = sort([vel_amp_edges,-vel_amp_edges]);

vel_amp_bins = discretize(max_vel,vel_amp_edges);
% don't use the middle bin
vel_amp_bins(vel_amp_bins == n_prctiles) = NaN;

% Plot wheel and activity within bins
[vis_grps_used,vis_max_vel_grp_mean] = ...
    grpstats(max_vel(vis_trials),vel_amp_bins(vis_trials),{'gname','nanmean'});
vis_grps_used = cellfun(@str2num,vis_grps_used);

[zero_grps_used,zero_max_vel_grp_mean] = ...
    grpstats(max_vel(zero_trials),vel_amp_bins(zero_trials),{'gname','nanmean'});
zero_grps_used = cellfun(@str2num,zero_grps_used);

use_activity = mua_allcat_move(:,:,2);
% use_activity = fluor_predicted_roi_reduced_move(:,:,7);

vis_wheel_grp_mean = grpstats(wheel_velocity_allcat_move(vis_trials,:),vel_amp_bins(vis_trials));
vis_activity_grp_mean = grpstats(use_activity(vis_trials,:),vel_amp_bins(vis_trials));

zero_wheel_grp_mean = grpstats(wheel_velocity_allcat_move(zero_trials,:),vel_amp_bins(zero_trials));
zero_activity_grp_mean = grpstats(use_activity(zero_trials,:),vel_amp_bins(zero_trials));

figure; 
col = [0.1*ones(n_prctiles,1),0.1*ones(n_prctiles,1),linspace(1,0.4,n_prctiles)'; ...
    linspace(0.4,1,n_prctiles)',0.1*ones(n_prctiles,1),0.1*ones(n_prctiles,1)];
col_used = col(unique([vis_grps_used;zero_grps_used]),:);

subplot(2,4,1); hold on
set(gca,'ColorOrder',col_used)
plot(t,vis_wheel_grp_mean','linewidth',2)
xlabel('Time from move');
ylabel('Wheel velocity');
axis tight
line([0,0],ylim,'color','k');

subplot(2,4,2); hold on
set(gca,'ColorOrder',col_used)
plot(t,vis_activity_grp_mean','linewidth',2)
xlabel('Time from move');
ylabel('Activity');
axis tight
line([0,0],ylim,'color','k');

subplot(2,4,3); hold on;
set(gca,'ColorOrder',col_used)
plot_t = t < 0.5;
plot(vis_wheel_grp_mean(:,plot_t)',vis_activity_grp_mean(:,plot_t)','linewidth',2)
xlim([-max(abs(xlim)),max(abs(xlim))]);
xlabel('Wheel speed');
ylabel('Activity');
axis tight
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');

subplot(2,4,5); hold on
set(gca,'ColorOrder',col_used)
plot(t,zero_wheel_grp_mean','linewidth',2)
xlabel('Time from move');
ylabel('Wheel velocity');
axis tight
line([0,0],ylim,'color','k');

subplot(2,4,6); hold on
set(gca,'ColorOrder',col_used)
plot(t,zero_activity_grp_mean','linewidth',2)
xlabel('Time from move');
ylabel('Activity');
axis tight
line([0,0],ylim,'color','k');

subplot(2,4,7); hold on;
set(gca,'ColorOrder',col_used)
plot_t = t < 0.5;
plot(zero_wheel_grp_mean(:,plot_t)',zero_activity_grp_mean(:,plot_t)','linewidth',2)
xlim([-max(abs(xlim)),max(abs(xlim))]);
xlabel('Wheel speed');
ylabel('Activity');
axis tight
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');

subplot(1,4,4); hold on;
p1 = plot(vis_max_vel_grp_mean,max(vis_activity_grp_mean(:,plot_t),[],2),'k','linewidth',2);
scatter(vis_max_vel_grp_mean,max(vis_activity_grp_mean(:,plot_t),[],2),80,col_used, ...
    'Filled','MarkerEdgeColor','k','linewidth',2);

p2 = plot(zero_max_vel_grp_mean,max(zero_activity_grp_mean(:,plot_t),[],2),'color',[0.7,0.7,0.7],'linewidth',2);
scatter(zero_max_vel_grp_mean,max(zero_activity_grp_mean(:,plot_t),[],2),80,col_used, ...
    'Filled','MarkerEdgeColor',[0.7,0.7,0.7],'linewidth',2);
xlabel('Max wheel velocity');
ylabel('Max activity');

legend([p1,p2],{'Visual','Zero-contrast'});
axis square;





% Summarize this: do this in V-space and make a pixel map? maybe compare
% visual and zero and then speed modulation?
use_rxn = move_t > 0 & move_t < 1;

move_onset_regressor_idx = strcmp('Move onset',regressor_labels);
fluor_move_activity = fluor_allcat_deconv - fluor_allcat_predicted_reduced(:,:,:,move_onset_regressor_idx);
mua_allcat_move = mua_allcat-mua_allcat_predicted_reduced(:,:,:,move_onset_regressor_idx);

t_leeway = 0.5;
leeway_samples = round(t_leeway*(sample_rate));
for i = 1:size(fluor_allcat_deconv,1)
    fluor_move_activity(i,:,:) = circshift(fluor_move_activity(i,:,:),-move_idx(i)+leeway_samples,2);  
    mua_allcat_move(i,:,:) = circshift(mua_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);  
end

vis_L_trials = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == -1 & ...
    use_rxn;

vis_R_trials = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == -1 & ...
    trial_choice_allcat == 1 & ...
    use_rxn;

zero_L_trials = ...
    trial_contrast_allcat == 0 & ...
    trial_choice_allcat == -1 & ...
    use_rxn;

zero_R_trials = ...
    trial_contrast_allcat == 0 & ...
    trial_choice_allcat == 1 & ...
    use_rxn;

trial_types = [vis_L_trials,vis_R_trials,zero_L_trials,zero_R_trials];
trial_types_labels = {'Vis L','Vis R','Zero L','Zero R'};
compare_types = {[1,2],[3,4],[1,3],[2,4]};

trial_types_mua = nan(n_depths,size(mua_allcat_move,2),size(trial_types,2));
trial_types_px = nan(size(U_master,1),size(U_master,2),size(fluor_move_activity,2),size(trial_types,2));
for curr_trial_type = 1:size(trial_types,2)   
    trial_types_mua(:,:,curr_trial_type) = ...
        permute(nanmean(mua_allcat_move(trial_types(:,curr_trial_type),:,:),1),[3,2,1]);   
    
    trial_types_px(:,:,:,curr_trial_type) = ....
        svdFrameReconstruct(U_master(:,:,1:n_vs), ...
        permute(nanmean(fluor_move_activity(trial_types(:,curr_trial_type),:,:),1),[3,2,1]));    
end

figure;
p1 = []; p2 = [];
for curr_plot = 1:4
    p1(curr_plot) = subplot(2,4,curr_plot);
    hold on; set(gca,'ColorOrder',copper(n_depths))
    plot(t,trial_types_mua(:,:,curr_plot)','linewidth',2);
    axis tight;
    line([0,0],ylim,'color','k');
    title(trial_types_labels{curr_plot});
end
for curr_plot = 1:length(compare_types)
    p2(curr_plot) = subplot(2,4,4+curr_plot);
    hold on; set(gca,'ColorOrder',copper(n_depths))
    plot(t,trial_types_mua(:,:,compare_types{curr_plot}(1))' - ...
        trial_types_mua(:,:,compare_types{curr_plot}(2))','linewidth',2);
    axis tight;
    line([0,0],ylim,'color','k');
    line(xlim,[0,0],'color','k');
    title([trial_types_labels{compare_types{curr_plot}(1)} '-' ...
        trial_types_labels{compare_types{curr_plot}(2)}]);
end
linkaxes([p1,p2],'xy');

trial_types_px_max = squeeze(max(trial_types_px(:,:,t < 0,:),[],3));
cmax = max(abs(trial_types_px_max(:)));
figure;
p1 = []; p2 = [];
for curr_plot = 1:4
    p1(curr_plot) = subplot(2,4,curr_plot);
    imagesc(trial_types_px_max(:,:,curr_plot));
    colormap(gca,brewermap([],'BuGn'));
    caxis([0,cmax]);
    axis image off;
    AP_reference_outline('ccf_aligned','k');
    title(trial_types_labels{curr_plot});
end
for curr_plot = 1:length(compare_types)
    p2(curr_plot) = subplot(2,4,4+curr_plot);
    imagesc(trial_types_px_max(:,:,compare_types{curr_plot}(1)) - ...
        trial_types_px_max(:,:,compare_types{curr_plot}(2)));
    colormap(gca,brewermap([],'*RdBu'));
    caxis([-cmax,cmax]);
    axis image off;
    AP_reference_outline('ccf_aligned','k');
    title(trial_types_labels{curr_plot});
    title([trial_types_labels{compare_types{curr_plot}(1)} '-' ...
        trial_types_labels{compare_types{curr_plot}(2)}]);
end























