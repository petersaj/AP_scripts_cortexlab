% (this is to work on AP_ctx_str_figures and related code)

%% Load choiceworld trial data (** NEEDED FOR BELOW **)

% Load data
% data_fn = 'trial_activity_choiceworld_FLUORTASKTEST';
% data_fn = 'trial_activity_choiceworld_MOVELTEST';
% data_fn = 'trial_activity_choiceworld_FWDCTXTEST';
% data_fn = 'trial_activity_choiceworld_STIMRTEST';
% data_fn = 'trial_activity_choiceworld_MOVEONGOINGTEST';
data_fn = 'trial_activity_choiceworld_EARLYMOVETEST';

exclude_data = true;
AP_load_concat_normalize_ctx_str;

% Get trial information
trial_contrast_allcat = max(D_allcat.stimulus,[],2);
[~,side_idx] = max(D_allcat.stimulus > 0,[],2);
trial_side_allcat = (side_idx-1.5)*2;
trial_contrastside_allcat = trial_contrast_allcat.*trial_side_allcat;
trial_choice_allcat = -(D_allcat.response-1.5)*2;

trial_outcome_allcat = D_allcat.outcome;
 
% % Get time (make this be saved in trial data)
sample_rate = 1/mean(diff(t));

% Get reaction time
[move_trial,move_idx] = max(abs(wheel_allcat) > 0.02,[],2);
move_t = t(move_idx)';
move_t(~move_trial) = Inf;

% Get outcome time
[~,outcome_idx] = max(any(outcome_allcat,3),[],2);
outcome_t = t(outcome_idx)';

% Get wheel velocity
wheel_velocity_allcat = wheel_allcat;
[max_speed,max_vel_idx] = max(abs(wheel_velocity_allcat(:,:,1).* ...
    (bsxfun(@times,wheel_velocity_allcat,trial_choice_allcat) > 0)),[],2);
max_vel = max_speed.*trial_choice_allcat;

% Make move-aligned data
wheel_velocity_allcat_move = wheel_velocity_allcat;
mua_allcat_move = mua_allcat;
mua_ctxpred_allcat_move = mua_ctxpred_allcat;
mua_taskpred_allcat_move = mua_taskpred_allcat;
mua_taskpred_reduced_allcat_move = mua_taskpred_reduced_allcat;
mua_ctxpred_taskpred_allcat_move = mua_taskpred_allcat;
mua_ctxpred_taskpred_reduced_allcat_move = mua_taskpred_reduced_allcat;

t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
for i = 1:size(mua_allcat,1)
    wheel_velocity_allcat_move(i,:,:) = circshift(wheel_velocity_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
    mua_allcat_move(i,:,:) = circshift(mua_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
    mua_ctxpred_allcat_move(i,:,:) = circshift(mua_ctxpred_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
    mua_taskpred_allcat_move(i,:,:) = circshift(mua_taskpred_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
    mua_taskpred_reduced_allcat_move(i,:,:,:) = circshift(mua_taskpred_reduced_allcat_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
    mua_ctxpred_taskpred_allcat_move(i,:,:) = circshift(mua_ctxpred_taskpred_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
    mua_ctxpred_taskpred_reduced_allcat_move(i,:,:,:) = circshift(mua_ctxpred_taskpred_reduced_allcat_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
end

disp('Finished loading trials')

%% Load concat task regression results (only cortex used)

% Load regression results
regression_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\trial_activity_choiceworld_task_regression';
load(regression_fn);
 
% Get task regression in ROIs
fluor_roi_deconv_taskpred = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_predicted,[3,2,1]),n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
    size(cat(3,wf_roi.mask),3),[],size(fluor_allcat_predicted,1)),[3,2,1]);

fluor_roi_deconv_taskpred_reduced = arrayfun(@(x) ...
    permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_predicted_reduced(:,:,:,x),[3,2,1]), ...
    n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
    size(cat(3,wf_roi.mask),3),[],size(fluor_allcat_predicted,1)),[3,2,1]), ...
    1:size(fluor_allcat_predicted_reduced,4),'uni',false);

fluor_roi_deconv_taskpred_reduced = cat(4,fluor_roi_deconv_taskpred_reduced{:});

disp('Finished loading concatenated regression')


%% Trial-based ctx/task->str prediction accuracy

trial_contrastside_allcat_early = trial_contrastside_allcat;
trial_contrastside_allcat_early(move_t > 0.5) = NaN;

%%%%%% Plot measured and predicted relating to each kernel
plot_areas = [1,2,3,1,4];
plot_reduction = [1,2,2,3,4]; % (which reduced variable to plot)
plot_conditions = {unique(trial_contrastside_allcat),[-1,1],1:3,[1,3],0:1};
plot_conditions_compare = { ...
    trial_contrastside_allcat_early, ...  
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
    
    p1 = subplot(3,length(plot_areas),curr_area);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(mua_allcat_move(curr_trials,:,area_idx) - ...
            mua_taskpred_reduced_allcat_move(curr_trials,:,area_idx,plot_reduction(area_idx)),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Str ' num2str(area_idx) ' Measured']);
    axis tight
    line([0,0],ylim,'color','k');       
    
    p2 = subplot(3,length(plot_areas),curr_area+length(plot_areas)*1);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(mua_ctxpred_allcat_move(curr_trials,:,area_idx) - ...
            mua_ctxpred_taskpred_reduced_allcat_move(curr_trials,:,area_idx,plot_reduction(area_idx)),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Cortex predicted']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p3 = subplot(3,length(plot_areas),curr_area+length(plot_areas)*2);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(mua_taskpred_allcat_move(curr_trials,:,area_idx) - ...
            mua_taskpred_reduced_allcat_move(curr_trials,:,area_idx,plot_reduction(area_idx)),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Task predicted']);
    axis tight
    line([0,0],ylim,'color','k');
    
    linkaxes([p1,p2,p3],'xy');
    
end

% Plot residuals
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
        curr_data_mean = nanmean(mua_allcat(curr_trials,:,area_idx) - ...
            mua_ctxpred_allcat(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Str ' num2str(area_idx) ' Ctx residual']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p3 = subplot(4,length(plot_areas),curr_area+length(plot_areas)*2);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(mua_allcat(curr_trials,:,area_idx) - ...
            mua_taskpred_allcat(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Str ' num2str(area_idx) ' Task residual']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p4 = subplot(4,length(plot_areas),curr_area+length(plot_areas)*3);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(mua_ctxpred_allcat(curr_trials,:,area_idx) - ...
            mua_taskpred_allcat(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Str ' num2str(area_idx) ' Ctx - Task']);
    axis tight
    line([0,0],ylim,'color','k');
    
    linkaxes([p1,p2,p3,p4],'xy');
    
end



% Plot R^2 (normalized data together)
r2_trials = move_t > 0 & move_t < 0.5;

predicted_points = +~any(isnan(cat(4,mua_allcat(r2_trials,:,:),mua_ctxpred_allcat(r2_trials,:,:),mua_taskpred_allcat(r2_trials,:,:))),4);
predicted_points(~predicted_points) = NaN;

sse_total_t = nansum((mua_allcat(r2_trials,:,:).*predicted_points-nanmean(mua_allcat(r2_trials,:,:).*predicted_points,1)).^2,1);
sse_residual_ctx_t = nansum((mua_allcat(r2_trials,:,:).*predicted_points-mua_ctxpred_allcat(r2_trials,:,:).*predicted_points).^2,1);
sse_residual_task_t = nansum((mua_allcat(r2_trials,:,:).*predicted_points-mua_taskpred_allcat(r2_trials,:,:).*predicted_points).^2,1);
explained_var_ctx_t = 1-(sse_residual_ctx_t./sse_total_t);
explained_var_task_t = 1-(sse_residual_task_t./sse_total_t);

sse_total = nansum((reshape(mua_allcat(r2_trials,:,:).*predicted_points,[],n_depths)- ...
    nanmean(reshape(mua_allcat(r2_trials,:,:).*predicted_points,[],n_depths),1)).^2,1);
sse_residual_ctx = nansum((reshape(mua_allcat(r2_trials,:,:).*predicted_points,[],n_depths)- ...
    reshape(mua_ctxpred_allcat(r2_trials,:,:).*predicted_points,[],n_depths)).^2,1);
sse_residual_task = nansum((reshape(mua_allcat(r2_trials,:,:).*predicted_points,[],n_depths)- ...
    reshape(mua_taskpred_allcat(r2_trials,:,:).*predicted_points,[],n_depths)).^2,1);
explained_var_ctx = 1-(sse_residual_ctx./sse_total);
explained_var_task = 1-(sse_residual_task./sse_total);

figure;
for curr_depth = 1:n_depths
    subplot(1,n_depths,curr_depth); hold on;
    p1 = plot(t,explained_var_ctx_t(:,:,curr_depth),'b','linewidth',2);
    line(xlim,repmat(explained_var_ctx(curr_depth),1,2),'color','b','linewidth',2);
    
    p2 = plot(t,explained_var_task_t(:,:,curr_depth),'r','linewidth',2);
    line(xlim,repmat(explained_var_task(curr_depth),1,2),'color','r','linewidth',2);
    
    ylabel('R^2');
    xlabel('Time from stim');
    legend([p1,p2],{'Cortex','Task'})
    ylim([0,1]);
end



% Plot R^2 (recordings separate)
r2_trials = move_t > 0 & move_t < 1;

% split everthing up by animal or recording
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));

use_split = trials_animal;

r2_trials_exp = mat2cell(r2_trials,use_split,1);

move_idx_exp = cellfun(@(x,trials) x(trials,:,:), ...
    mat2cell(move_idx,use_split,1),r2_trials_exp,'uni',false);
mua_exp = cellfun(@(x,trials) x(trials,:,:), ...
    mat2cell(mua_allcat,use_split,length(t),n_depths),r2_trials_exp,'uni',false);
mua_ctxpred_exp = cellfun(@(x,trials) x(trials,:,:), ...
    mat2cell(mua_ctxpred_allcat,use_split,length(t),n_depths),r2_trials_exp,'uni',false);
mua_taskpred_exp = cellfun(@(x,trials) x(trials,:,:), ...
    mat2cell(mua_taskpred_allcat,use_split,length(t),n_depths),r2_trials_exp,'uni',false);

mua_taskpred_reduced_exp = cellfun(@(x,trials) x(trials,:,:,:), ...
    mat2cell(mua_taskpred_reduced_allcat,use_split,length(t),n_depths, ...
    size(mua_taskpred_reduced_allcat,4)),r2_trials_exp,'uni',false);

% remove_regressor = 2;
% mua_exp = cellfun(@(x,y) x-y(:,:,:,remove_regressor),mua_exp,mua_taskpred_reduced_exp,'uni',false);
% mua_ctxpred_exp = cellfun(@(x,y) x-y(:,:,:,remove_regressor),mua_ctxpred_exp,mua_taskpred_reduced_exp,'uni',false);
% mua_taskpred_exp = cellfun(@(x,y) x-y(:,:,:,remove_regressor),mua_taskpred_exp,mua_taskpred_reduced_exp,'uni',false);

for curr_exp = 1:length(mua_exp)

    % Set common nans
    nan_points = any(isnan(cat(4,mua_exp{curr_exp}, ...
        mua_ctxpred_exp{curr_exp},mua_taskpred_exp{curr_exp})),4);
    mua_exp{curr_exp}(nan_points) = NaN;
    mua_ctxpred_exp{curr_exp}(nan_points) = NaN;
    mua_taskpred_exp{curr_exp}(nan_points) = NaN;
     
%     % Low-pass filter data
%     lowpassCutoff = 6; % Hz
%     [b100s, a100s] = butter(2, lowpassCutoff/((sample_rate)/2), 'low');
%     mua_exp{curr_exp}(nan_points) = 0;
%     mua_exp{curr_exp} = permute(filtfilt(b100s,a100s,permute(mua_exp{curr_exp},[2,1,3])),[2,1,3]);
%     mua_exp{curr_exp}(nan_points) = NaN;
    
%     % (to align to movement)    
%     t_leeway = t(1);
%     leeway_samples = round(-t_leeway*sample_rate);
%     for i = 1:size(move_idx_exp{curr_exp},1)
%         mua_exp{curr_exp}(i,:,:) = circshift(mua_exp{curr_exp}(i,:,:),-move_idx_exp{curr_exp}(i)+leeway_samples,2);
%         mua_ctxpred_exp{curr_exp}(i,:,:) = circshift(mua_ctxpred_exp{curr_exp}(i,:,:),-move_idx_exp{curr_exp}(i)+leeway_samples,2);
%         mua_taskpred_exp{curr_exp}(i,:,:) = circshift(mua_taskpred_exp{curr_exp}(i,:,:),-move_idx_exp{curr_exp}(i)+leeway_samples,2);
%     end
end

exp_sse_total = cell2mat(cellfun(@(measured) ...
    nansum((reshape(measured,[],n_depths)-nanmean(reshape(measured,[],n_depths),1)).^2,1),mua_exp,'uni',false));
exp_sse_residual_ctx = cell2mat(cellfun(@(measured,ctxpred) ...
    nansum((reshape(measured,[],n_depths)-reshape(ctxpred,[],n_depths)).^2,1),mua_exp,mua_ctxpred_exp,'uni',false));
exp_sse_residual_task = cell2mat(cellfun(@(measured,taskpred) ...
    nansum((reshape(measured,[],n_depths)-reshape(taskpred,[],n_depths)).^2,1),mua_exp,mua_taskpred_exp,'uni',false));

exp_explained_var_ctx = 1-(exp_sse_residual_ctx./exp_sse_total);
exp_explained_var_task = 1-(exp_sse_residual_task./exp_sse_total);


exp_sse_total_t = cell2mat(cellfun(@(measured) ...
    nansum((measured-nanmean(measured,1)).^2,1),mua_exp,'uni',false));
exp_sse_residual_ctx_t = cell2mat(cellfun(@(measured,ctxpred) ...
    nansum((measured-ctxpred).^2,1),mua_exp,mua_ctxpred_exp,'uni',false));
exp_sse_residual_task_t = cell2mat(cellfun(@(measured,taskpred) ...
    nansum((measured-taskpred).^2,1),mua_exp,mua_taskpred_exp,'uni',false));

exp_explained_var_ctx_t = 1-(exp_sse_residual_ctx_t./exp_sse_total_t);
exp_explained_var_task_t = 1-(exp_sse_residual_task_t./exp_sse_total_t);

exp_explained_var_ctx_t(exp_explained_var_ctx_t < 0) = 0;
exp_explained_var_task_t(exp_explained_var_task_t < 0) = 0;

figure;
for curr_depth = 1:n_depths
    p1 = subplot(3,n_depths,curr_depth); 
    imagesc(t,[],exp_explained_var_ctx_t(:,:,curr_depth)); caxis([0,1]); colormap(hot);
    title('Expl var cortex');
    ylabel('Experiment');
    
    p2 = subplot(3,n_depths,curr_depth+n_depths); 
    imagesc(t,[],exp_explained_var_task_t(:,:,curr_depth)); caxis([0,1]); colormap(hot);
    title('Expl var task');
    ylabel('Experiment');
    
    p3 = subplot(3,n_depths,curr_depth+n_depths*2);     
    hold on;
    plot(t,nanmean(exp_explained_var_ctx_t(:,:,curr_depth),1),'b','linewidth',2);
    plot(t,nanmean(exp_explained_var_task_t(:,:,curr_depth),1),'r','linewidth',2);
    plot(t,nanmean(exp_explained_var_ctx_t(:,:,curr_depth),1) - ...
        nanmean(exp_explained_var_task_t(:,:,curr_depth),1),'k','linewidth',1);
    line(xlim,repmat(nanmean(exp_explained_var_ctx(:,curr_depth)),2,1),'color','b','linewidth',2);
    line(xlim,repmat(nanmean(exp_explained_var_task(:,curr_depth)),2,1),'color','r','linewidth',2);
    ylabel('R^2');
    xlabel('Time');
    legend({'Cortex','Task','Cortex-Task'})
    
    linkaxes([p1,p2,p3],'x');
end



% Get rank differences for stim/move/reward
regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
trial_groups = {'Stim','Move onset','Outcome'};
figure;
for curr_group = 1:length(trial_groups)
    
    % Get MUA/predicted using reduced model
    curr_regressor_idx = strcmp(trial_groups{curr_group},regressor_labels);
    curr_mua =  mua_allcat - mua_taskpred_reduced_allcat(:,:,:,curr_regressor_idx);
    curr_mua_ctx = mua_ctxpred_allcat - mua_ctxpred_taskpred_reduced_allcat(:,:,:,curr_regressor_idx);
    
    % Set common NaNs
    nan_samples = isnan(curr_mua) | isnan(curr_mua_ctx);
    curr_mua(nan_samples) = NaN;
    curr_mua_ctx(nan_samples) = NaN;

    % (movement: align to move onset)
    if any(strfind(lower(trial_groups{curr_group}),'move'))
        t_leeway = -t(1);
        leeway_samples = round(t_leeway*(sample_rate));
        for i = 1:size(fluor_allcat_deconv,1)
            curr_mua(i,:,:) = circshift(curr_mua(i,:,:),-move_idx(i)+leeway_samples,2);
            curr_mua_ctx(i,:,:) = circshift(curr_mua_ctx(i,:,:),-move_idx(i)+leeway_samples,2);
        end
    end
    
    % (outcome: align to outcome)
    if any(strfind(lower(trial_groups{curr_group}),'outcome'))       
        t_leeway = -t(1);
        leeway_samples = round(t_leeway*(sample_rate));
        for i = 1:size(fluor_allcat_deconv,1)
            curr_mua(i,:,:) = circshift(curr_mua(i,:,:),-outcome_idx(i)+leeway_samples,2);
            curr_mua_ctx(i,:,:) = circshift(curr_mua_ctx(i,:,:),-outcome_idx(i)+leeway_samples,2);
        end
    end
    
    % Set trials and grouping to use
    switch trial_groups{curr_group}
        case 'Stim'
            use_trials = move_t > 0 & move_t < 0.5 & trial_contrast_allcat > 0;
            trial_group = trial_side_allcat;
            group1 = 1;
            group2 = -1;
        case 'Move onset'
            use_trials = move_t > 0 & move_t < 0.5;
            trial_group = trial_choice_allcat;
            group1 = -1;
            group2 = 1;
        case 'Outcome'
            use_trials = move_t > 0 & move_t < 0.5;
            trial_group = trial_choice_allcat == -trial_side_allcat;
            group1 = 1;
            group2 = 0;
    end    
    
    act_rank = tiedrank(curr_mua(use_trials,:,:));
    predicted_act_rank = tiedrank(curr_mua_ctx(use_trials,:,:));
    
    act_rank_difference = squeeze(( ...
        nanmean(act_rank(trial_group(use_trials) == group1,:,:),1) - ...
        nanmean(act_rank(trial_group(use_trials) == group2,:,:),1))./max(act_rank,[],1)); 
    predicted_act_rank_difference = squeeze(( ...
        nanmean(predicted_act_rank(trial_group(use_trials) == group1,:,:),1) - ...
        nanmean(predicted_act_rank(trial_group(use_trials) == group2,:,:),1))./max(predicted_act_rank,[],1));     
    
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


% Rank differences by experiment
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_animal;

mua_exp = mat2cell(mua_allcat,use_split,length(t),n_depths);
mua_ctxpred_exp = mat2cell(mua_ctxpred_allcat,use_split,length(t),n_depths);

mua_taskpred_reduced_exp = mat2cell(mua_taskpred_reduced_allcat,use_split,length(t),n_depths, ...
    size(mua_taskpred_reduced_allcat,4));
mua_ctxpred_taskpred_reduced_exp = mat2cell(mua_ctxpred_taskpred_reduced_allcat,use_split,length(t),n_depths, ...
    size(mua_ctxpred_taskpred_reduced_allcat,4));

move_t_exp = mat2cell(move_t,use_split,1);
move_idx_exp = mat2cell(move_idx,use_split,1);
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
trial_side_allcat_exp = mat2cell(trial_side_allcat,use_split,1);
trial_contrast_allcat_exp = mat2cell(trial_contrast_allcat,use_split,1);
trial_choice_allcat_exp = mat2cell(trial_choice_allcat,use_split,1);

regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
trial_groups = {'Stim','Move onset','Outcome'};
t_groups = {[0,0.3],[-0.2,0.2],[0,0.2]};

act_rank_difference = nan(length(t),n_depths,length(use_split),length(trial_groups));
predicted_act_rank_difference = nan(length(t),n_depths,length(use_split),length(trial_groups));

act_rank_difference_trial = nan(n_depths,length(use_split),length(trial_groups));
predicted_act_rank_difference_trial = nan(n_depths,length(use_split),length(trial_groups));

for curr_exp = 1:length(mua_exp)
    for curr_group = 1:length(trial_groups)
        
        % Get MUA/predicted using reduced model
        curr_regressor_idx = strcmp(trial_groups{curr_group},regressor_labels);
        curr_mua =  mua_exp{curr_exp} - mua_taskpred_reduced_exp{curr_exp}(:,:,:,curr_regressor_idx);
        curr_mua_ctx = mua_ctxpred_exp{curr_exp} - mua_ctxpred_taskpred_reduced_exp{curr_exp}(:,:,:,curr_regressor_idx);
        
        % Skip if there's no data in this experiment
        if isempty(curr_mua)
            continue
        end
        
        % Set common NaNs
        nan_samples = isnan(curr_mua) | isnan(curr_mua_ctx);
        curr_mua(nan_samples) = NaN;
        curr_mua_ctx(nan_samples) = NaN;
        
        % (movement: align to move onset)
        if any(strfind(lower(trial_groups{curr_group}),'move'))
            t_leeway = -t(1);
            leeway_samples = round(t_leeway*(sample_rate));
            for i = 1:size(curr_mua,1)
                curr_mua(i,:,:) = circshift(curr_mua(i,:,:),-move_idx_exp{curr_exp}(i)+leeway_samples,2);
                curr_mua_ctx(i,:,:) = circshift(curr_mua_ctx(i,:,:),-move_idx_exp{curr_exp}(i)+leeway_samples,2);
            end
        end
        
        % (outcome: align to outcome)
        if any(strfind(lower(trial_groups{curr_group}),'outcome'))
            t_leeway = -t(1);
            leeway_samples = round(t_leeway*(sample_rate));
            for i = 1:size(curr_mua,1)
                curr_mua(i,:,:) = circshift(curr_mua(i,:,:),-outcome_idx_exp{curr_exp}(i)+leeway_samples,2);
                curr_mua_ctx(i,:,:) = circshift(curr_mua_ctx(i,:,:),-outcome_idx_exp{curr_exp}(i)+leeway_samples,2);
            end
        end
        
        % Set trials and grouping to use
        switch trial_groups{curr_group}
            case 'Stim'
                use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} < 0.5 & trial_contrast_allcat_exp{curr_exp} > 0;
                trial_group = trial_side_allcat_exp{curr_exp};
                group1 = 1;
                group2 = -1;
            case 'Move onset'
                use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} < 0.5;
                trial_group = trial_choice_allcat_exp{curr_exp};
                group1 = -1;
                group2 = 1;                
            case 'Outcome'
                use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} < 0.5;
                trial_group = trial_choice_allcat_exp{curr_exp} == -trial_side_allcat_exp{curr_exp};
                group1 = 1;
                group2 = 0;
        end
        
        act_rank = tiedrank(curr_mua(use_trials,:,:));
        predicted_act_rank = tiedrank(curr_mua_ctx(use_trials,:,:));
        
        act_rank_difference(:,:,curr_exp,curr_group) = squeeze(( ...
            nanmean(act_rank(trial_group(use_trials) == group1,:,:),1) - ...
            nanmean(act_rank(trial_group(use_trials) == group2,:,:),1))./max(act_rank,[],1));
        predicted_act_rank_difference(:,:,curr_exp,curr_group) = squeeze(( ...
            nanmean(predicted_act_rank(trial_group(use_trials) == group1,:,:),1) - ...
            nanmean(predicted_act_rank(trial_group(use_trials) == group2,:,:),1))./max(predicted_act_rank,[],1));
        
        use_t = t >= t_groups{curr_group}(1) & t <= t_groups{curr_group}(2);
        act_rank_trial = tiedrank(squeeze(nanmean(curr_mua(use_trials,use_t,:),2)));
        predicted_act_rank_trial = tiedrank(squeeze(nanmean(curr_mua_ctx(use_trials,use_t,:),2)));
        
        act_rank_difference_trial(:,curr_exp,curr_group) = ...
            (nanmean(act_rank_trial(trial_group(use_trials) == group1,:),1) - ...
            nanmean(act_rank_trial(trial_group(use_trials) == group2,:),1))./max(act_rank_trial,[],1);
        predicted_act_rank_difference_trial(:,curr_exp,curr_group) = ...
            (nanmean(predicted_act_rank_trial(trial_group(use_trials) == group1,:),1) - ...
            nanmean(predicted_act_rank_trial(trial_group(use_trials) == group2,:),1))./max(predicted_act_rank_trial,[],1);
        
%         n_shuff = 100;
%         shuff_rank_difference = nan(length(t),n_depths,n_shuff);
%         shuff_predicted_rank_difference = nan(length(t),n_depths,n_shuff);
%         for curr_shuff = 1:n_shuff
%             curr_shake_group = AP_shake(trial_group(use_trials));
%             shuff_rank_difference(:,:,curr_shuff) = squeeze(( ...
%                 nanmean(act_rank(curr_shake_group == group1,:,:),1) - ...
%                 nanmean(act_rank(curr_shake_group == group2,:,:),1))./max(act_rank,[],1));
%             shuff_predicted_rank_difference(:,:,curr_shuff) = squeeze(( ...
%                 nanmean(predicted_act_rank(curr_shake_group == group1,:,:),1) - ...
%                 nanmean(predicted_act_rank(curr_shake_group == group2,:,:),1))./max(predicted_act_rank,[],1));
%             AP_print_progress_fraction(curr_shuff,n_shuff);
%         end
        
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
p1 = subplot(1,2,1);
errorbar(squeeze(nanmean(act_rank_difference_trial,2)), ...
    squeeze(nanstd(act_rank_difference_trial,[],2)./ ...
    sqrt(sum(~isnan(act_rank_difference_trial),2))),'linewidth',2);
xlabel('Striatum depth');
ylabel('Rank difference');
legend(trial_groups);
title('Striatum');

p2 = subplot(1,2,2);
errorbar(squeeze(nanmean(predicted_act_rank_difference_trial,2)), ...
    squeeze(nanstd(predicted_act_rank_difference_trial,[],2)./ ...
    sqrt(sum(~isnan(predicted_act_rank_difference_trial),2))),'linewidth',2);
xlabel('Striatum depth');
ylabel('Rank difference');
legend(trial_groups);
title('Cortex-predicted striatum');

linkaxes([p1,p2]);



% Plot long regression example
plot_areas = 1:4;
figure;
for curr_area = 1:length(plot_areas)
    plot_area_idx = plot_areas(curr_area);
    long_trace = reshape(permute(mua_allcat(:,:,plot_area_idx),[2,1]),[],1);
    long_trace_ctx_predicted = reshape(permute(mua_ctxpred_allcat(:,:,plot_area_idx),[2,1]),[],1);
    long_trace_task_predicted = reshape(permute(mua_taskpred_allcat(:,:,plot_area_idx),[2,1]),[],1);
    
    % Low-pass filter data
    lowpassCutoff = 6; % Hz
    [b100s, a100s] = butter(2, lowpassCutoff/((sample_rate)/2), 'low');
    long_trace(~isnan(long_trace)) = filtfilt(b100s,a100s,long_trace(~isnan(long_trace)));
    
    % (correlate measured and predicted in chunks)
    chunk_t = 60; %s
    chunk_samples = round(chunk_t*sample_rate);
    chunk_boundaries = 1:chunk_samples:length(long_trace);
    n_chunks = length(chunk_boundaries)-1;
    chunk_corr = nan(n_chunks,1);
    for curr_chunk = 1:length(chunk_boundaries)-1
        chunk_samples = chunk_boundaries(curr_chunk):chunk_boundaries(curr_chunk+1);
        nonan_chunk_samples = chunk_samples(~isnan(long_trace(chunk_samples)) & ...
            ~isnan(long_trace_ctx_predicted(chunk_samples)));
        
        curr_corr = corrcoef(long_trace(nonan_chunk_samples), ...
            long_trace_ctx_predicted(nonan_chunk_samples));
        chunk_corr(curr_chunk) = curr_corr(2);

%         curr_corr = sum((long_trace(nonan_chunk_samples)- ...
%             long_trace_ctx_predicted(nonan_chunk_samples)).^2);       
%         chunk_corr(curr_chunk) = curr_corr;       
    end
    
    % (plot chunk percentiles)
    plot_prctiles = [25,50,75];
    [chunk_corr_sorted,sort_idx] = sort(chunk_corr);
    chunk_corr_prctile = prctile(chunk_corr,plot_prctiles);
    for curr_prctile = 1:length(chunk_corr_prctile)
        subplot(length(chunk_corr_prctile),length(plot_areas), ...
            (curr_prctile-1)*length(plot_areas)+curr_area); hold on;
        curr_plot_chunk_idx = sort_idx(find(chunk_corr_sorted >= chunk_corr_prctile(curr_prctile),1));
        curr_chunk_plot = chunk_boundaries(curr_plot_chunk_idx):chunk_boundaries(curr_plot_chunk_idx+1);
        plot([1:length(curr_chunk_plot)]/sample_rate,long_trace(curr_chunk_plot),'color','k');
        plot([1:length(curr_chunk_plot)]/sample_rate,long_trace_ctx_predicted(curr_chunk_plot),'color',[0,0.6,0]);
        plot([1:length(curr_chunk_plot)]/sample_rate,long_trace_task_predicted(curr_chunk_plot),'color',[0.3,0.3,0.8]);
        ylabel([num2str(plot_prctiles(curr_prctile)) 'th percentile']);
        if curr_prctile == 1
            title(plot_areas(curr_area))
        end
        xlabel('Time (s)')
    end
end


% Plot all on top of each other 
figure; hold on;
curr_depth = 2;
plot(reshape(mua_allcat(:,:,curr_depth)',[],1));
plot(reshape(mua_ctxpred_allcat(:,:,curr_depth)',[],1));
plot(reshape(mua_taskpred_allcat(:,:,curr_depth)',[],1));





% Coherence (concatenated)
coherence_trials = move_t > 0 & move_t < 0.5;

mua_ctxpred_coherence = cell(n_depths,1);
mua_taskpred_coherence = cell(n_depths,1);
mua_ctxtaskpred_coherence = cell(n_depths,1);
for curr_depth = 1:n_depths    
    measured_data = double(reshape(mua_allcat(coherence_trials,:,curr_depth)',[],1));
    ctx_predicted_data = double(reshape(mua_ctxpred_allcat(coherence_trials,:,curr_depth)',[],1));
    task_predicted_data = double(reshape(mua_taskpred_allcat(coherence_trials,:,curr_depth)',[],1));
    
    nan_points = isnan(measured_data) | isnan(ctx_predicted_data) | isnan(task_predicted_data);
    
    [mua_ctxpred_coherence{curr_depth},coherence_f] = mscohere( ...
        measured_data(~nan_points), ctx_predicted_data(~nan_points), ...
        hamming(round(sample_rate*5)),round(sample_rate*2.5),[],sample_rate);   
    
    [mua_taskpred_coherence{curr_depth},coherence_f] = mscohere( ...
        measured_data(~nan_points), task_predicted_data(~nan_points), ...
        hamming(round(sample_rate*5)),round(sample_rate*2.5),[],sample_rate);
    
    [mua_ctxtaskpred_coherence{curr_depth},coherence_f] = mscohere( ...
        ctx_predicted_data(~nan_points), task_predicted_data(~nan_points), ...
        hamming(round(sample_rate*5)),round(sample_rate*2.5),[],sample_rate);
    
    AP_print_progress_fraction(curr_depth,n_depths);
end
mua_ctxpred_coherence = horzcat(mua_ctxpred_coherence{:});
mua_taskpred_coherence = horzcat(mua_taskpred_coherence{:});
mua_ctxtaskpred_coherence = horzcat(mua_ctxtaskpred_coherence{:});

figure; 
subplot(2,2,1); hold on; set(gca,'ColorOrder',copper(n_depths)); set(gca,'XScale','log')
plot(coherence_f,mua_ctxpred_coherence,'linewidth',2);
xlabel('Frequency'); ylabel('Coherence with measured');
title('Cortex predicted');
subplot(2,2,2); hold on; set(gca,'ColorOrder',copper(n_depths)); set(gca,'XScale','log')
plot(coherence_f,mua_taskpred_coherence,'linewidth',2);
xlabel('Frequency'); ylabel('Coherence with measured');
title('Task predicted');
subplot(2,2,3); hold on; set(gca,'ColorOrder',copper(n_depths)); set(gca,'XScale','log')
plot(coherence_f,mua_ctxpred_coherence-mua_taskpred_coherence,'linewidth',2);
xlabel('Frequency'); ylabel('Coherence with measured');
title('Cortex - task predicted');
subplot(2,2,4); hold on; set(gca,'ColorOrder',copper(n_depths)); set(gca,'XScale','log')
plot(coherence_f,mua_ctxtaskpred_coherence,'linewidth',2);
xlabel('Frequency'); ylabel('Coherence ctx/task');
title('Cortex predicted / task predicted');



% Coherence (by experiment)
coherence_trials = move_t > 0 & move_t < 0.5;

% split everthing up by animal or recording
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));

use_split = trials_recording;

coherence_trials_exp = mat2cell(coherence_trials,trials_recording,1);
move_idx_exp = mat2cell(move_idx,trials_recording,1);
mua_exp = mat2cell(mua_allcat,trials_recording,length(t),n_depths);
mua_ctxpred_exp = mat2cell(mua_ctxpred_allcat,trials_recording,length(t),n_depths);
mua_taskpred_exp = mat2cell(mua_taskpred_allcat,trials_recording,length(t),n_depths);

exp_mua_ctxpred_coherence = nan(length(use_split),length(coherence_f),n_depths);
exp_mua_taskpred_coherence = nan(length(use_split),length(coherence_f),n_depths);
for curr_exp = 1:length(mua_exp)
    
    % Set common nans
    nan_points = any(isnan(cat(4,mua_exp{curr_exp}, ...
        mua_ctxpred_exp{curr_exp},mua_taskpred_exp{curr_exp})),4);
    mua_exp{curr_exp}(nan_points) = NaN;
    mua_ctxpred_exp{curr_exp}(nan_points) = NaN;
    mua_taskpred_exp{curr_exp}(nan_points) = NaN;
        
    for curr_depth = 1:n_depths
        
        measured_data = double(reshape(mua_exp{curr_exp}(coherence_trials_exp{curr_exp},:,curr_depth)',[],1));
        ctx_predicted_data = double(reshape(mua_ctxpred_exp{curr_exp}(coherence_trials_exp{curr_exp},:,curr_depth)',[],1));
        task_predicted_data = double(reshape(mua_taskpred_exp{curr_exp}(coherence_trials_exp{curr_exp},:,curr_depth)',[],1));
        
        nan_points = isnan(measured_data) | isnan(ctx_predicted_data) | isnan(task_predicted_data);
        
        if all(nan_points)
            continue
        end
        
        [exp_mua_ctxpred_coherence(curr_exp,:,curr_depth),coherence_f] = mscohere( ...
            measured_data(~nan_points), ctx_predicted_data(~nan_points), ...
            hamming(round(sample_rate*5)),round(sample_rate*2.5),[],sample_rate);
        
        [exp_mua_taskpred_coherence(curr_exp,:,curr_depth),coherence_f] = mscohere( ...
            measured_data(~nan_points), task_predicted_data(~nan_points), ...
            hamming(round(sample_rate*5)),round(sample_rate*2.5),[],sample_rate);
           
    end
    AP_print_progress_fraction(curr_exp,length(mua_exp));
end

figure; 
subplot(1,3,1); hold on; set(gca,'ColorOrder',copper(n_depths)); set(gca,'XScale','log')
plot(coherence_f,squeeze(nanmean(exp_mua_ctxpred_coherence,1)),'linewidth',2);
xlabel('Frequency'); ylabel('Coherence with measured');
title('Cortex predicted');
subplot(1,3,2); hold on; set(gca,'ColorOrder',copper(n_depths)); set(gca,'XScale','log')
plot(coherence_f,squeeze(nanmean(exp_mua_taskpred_coherence,1)),'linewidth',2);
xlabel('Frequency'); ylabel('Coherence with measured');
title('Task predicted');
subplot(1,3,3); hold on; set(gca,'ColorOrder',copper(n_depths)); set(gca,'XScale','log')
plot(coherence_f,squeeze(nanmean(exp_mua_ctxpred_coherence,1))-squeeze(nanmean(exp_mua_taskpred_coherence,1)),'linewidth',2);
xlabel('Frequency'); ylabel('Coherence with measured');
title('Cortex - task predicted');





% % Apply empirical static nonlinearity
% figure;
% mua_allcat_predicted_nlin = mua_allcat;
% for curr_depth = 1:n_depths       
%     measured_data = double(reshape(mua_allcat(:,:,curr_depth)',[],1));
%     predicted_data = double(reshape(mua_ctxpred_allcat(:,:,curr_depth)',[],1));
%     predicted_data(predicted_data < 0) = 0;
%     
%     n_bins = 100;
%     activity_bounds = linspace(0,5,n_bins+1);
%     activity_bin_centers = conv2(activity_bounds,[1,1]/2,'valid');
%   
%     measured_bins = discretize(measured_data,activity_bounds);
%     predicted_bins = discretize(predicted_data,activity_bounds);
%     
%     measured_data_binmean = accumarray(predicted_bins(~isnan(predicted_bins)), ...
%         measured_data(~isnan(predicted_bins)),[n_bins,1],@nanmedian,nan);
%     predicted_data_binmean = accumarray(predicted_bins(~isnan(predicted_bins)),...
%         predicted_data(~isnan(predicted_bins)),[n_bins,1],@nanmedian,nan);
%     
%     % smooth out the measured data binmean
%     measured_data_binmean_smooth = medfilt1(measured_data_binmean,10);
%     
%     predicted_data_nlin = nan(size(predicted_data));
%     predicted_data_nlin(~isnan(predicted_bins)) = measured_data_binmean_smooth(predicted_bins(~isnan(predicted_bins)));
%     
%     predicted_data_nlin_bins = discretize(predicted_data_nlin,activity_bounds);
%     predicted_data_nlin_binmean = accumarray( ...
%         predicted_bins(~isnan(predicted_bins)), ...
%         predicted_data_nlin(~isnan(predicted_bins)),[n_bins,1],@mean,nan);
%     
%     mua_allcat_predicted_nlin(:,:,curr_depth) = ...
%         reshape(predicted_data_nlin, ...
%         size(mua_allcat,2),size(mua_allcat,1))';
%     
%     subplot(1,n_depths,curr_depth); hold on;
%     plot(predicted_data,measured_data,'.')
%     plot(predicted_data_binmean,measured_data_binmean,'linewidth',2);
%     plot(predicted_data_binmean,measured_data_binmean_smooth,'linewidth',2);
%     plot(predicted_data_nlin_binmean,measured_data_binmean,'linewidth',2);
%     xlim([-2,6]);ylim([-2,6]);
%     line(xlim,ylim,'color','k');
%     xlabel('Predicted')
%     ylabel('Measured')
%     axis square;
% end

%% Trial-based str vs str ROI

% %%%%% TESTING individual exp basis? this didn't work at all

% kernel_t = (round(regression_params.kernel_t(1)*sample_rate): ...
%             round(regression_params.kernel_t(2)*sample_rate))/sample_rate;
%         
% ctx_str_k_avg = fliplr(nanmean(cell2mat(permute(vertcat(ctx_str_k_all{:}),[2,3,4,1])),4));
% ctx_str_k_px_avg = reshape(svdFrameReconstruct( ...
%     U_master(:,:,1:size(ctx_str_k_avg,1)), ...
%     reshape(ctx_str_k_avg,size(ctx_str_k_avg,1),[])), ...
%     size(U_master,1),size(U_master,2),[],size(ctx_str_k_avg,3));
% 
% AP_image_scroll(ctx_str_k_px_avg,kernel_t);
% axis image; caxis([-max(abs(caxis)),max(abs(caxis))]);
% colormap(brewermap([],'*RdBu'));
% AP_reference_outline('ccf_aligned','k');
% 
% kernel_roi_t = nan(size(mua_allcat));
% for curr_depth = 1:n_depths
%     temp_fluor = reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]);
%     temp_fluor(isnan(temp_fluor)) = 0;
%     kernel_roi_t(:,:,curr_depth) = ...
%         reshape(sum(convn(temp_fluor,fliplr(ctx_str_k_avg(:,:,curr_depth)),'same'),1),length(t),size(mua_allcat,1))';
%     AP_print_progress_fraction(curr_depth,n_depths);
% end

% 
% mua_ctxpred_allcat_test = [];
% 
% for curr_animal = 1:length(mua_all)
%     for curr_day = 1:length(mua_all{curr_animal})
%         deconv_fluorescence = [];
%         for curr_depth = 1:n_depths
%             
%             test_k = ctx_str_k_all{curr_animal}{curr_day}(:,:,curr_depth);
%             test_fluor = AP_deconv_wf(fluor_all{curr_animal}{curr_day});
%             test_str = mua_all{curr_animal}{curr_day}(:,:,curr_depth);
%             test_ctx_str = mua_ctxpred_all{curr_animal}{curr_day}(:,:,curr_depth);
%             
%             %%%% COPIED FROM AP_deconv_wf
%             % Get valid time from kernel (assume time-symmetric)
%             t_rel = 1:size(test_fluor,2);
%             t_rel_conv_valid = round(conv(t_rel,ones(1,size(test_k,2))./size(test_k,2),'valid'));
%             
%             % Remove NaNs to make convolution faster (put back in later)
%             fluorescence_nonan = test_fluor;
%             fluorescence_nonan(isnan(fluorescence_nonan)) = 0;
%             
%             % Deconvolve widefield, set invalid points to NaN
%             dims = 1:ndims(test_fluor);
%             dims_permute = circshift(dims,-1);
%             [~,dims_unpermute] = sort(dims_permute);
%             
%             deconv_fluorescence(:,:,curr_depth) = permute(interp1(t_rel_conv_valid, ...
%                 permute(convn(fluorescence_nonan(:,:,1:50),permute(test_k,[3,2,1]),'valid'), ...
%                 dims_permute),t_rel),dims_unpermute);      
%         end
%         mua_ctxpred_allcat_test = [mua_ctxpred_allcat_test;deconv_fluorescence];        
%     end
%     AP_print_progress_fraction(curr_animal,length(mua_all));
% end
% 
% %%%%%%

% Load widefield kernel ROIs
kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
load(kernel_roi_fn);

kernel_roi_bw_hemidiff = kernel_roi.bw - AP_reflect_widefield(kernel_roi.bw);

% Get predicted fluorescence in ROIs
fluor_kernel_roi_bw = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],kernel_roi.bw), ...
    size(kernel_roi.bw,3),[],size(fluor_allcat_deconv,1)),[3,2,1]);

fluor_kernel_roi_bw_hemidiff = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],kernel_roi_bw_hemidiff), ...
    size(kernel_roi_bw_hemidiff,3),[],size(fluor_allcat_deconv,1)),[3,2,1]);

fluor_kernel_roi_weighted = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],kernel_roi.max_weighted), ...
    size(kernel_roi.max_weighted,3),[],size(fluor_allcat_deconv,1)),[3,2,1]);


fluor_kernel_roi_bw_move = fluor_kernel_roi_bw;
fluor_kernel_roi_bw_hemidiff_move = fluor_kernel_roi_bw_hemidiff;
fluor_kernel_roi_weighted_move = fluor_kernel_roi_weighted;

t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
for i = 1:size(mua_allcat,1)
    fluor_kernel_roi_bw_move(i,:,:) = circshift(fluor_kernel_roi_bw_move(i,:,:),-move_idx(i)+leeway_samples,2);
    fluor_kernel_roi_bw_hemidiff_move(i,:,:) = circshift(fluor_kernel_roi_bw_hemidiff_move(i,:,:),-move_idx(i)+leeway_samples,2);
    fluor_kernel_roi_weighted_move(i,:,:) = circshift(fluor_kernel_roi_weighted_move(i,:,:),-move_idx(i)+leeway_samples,2);
end


trial_contrastside_allcat_early = trial_contrastside_allcat;
trial_contrastside_allcat_early(move_t < 0 | move_t > 0.5) = NaN;

%%%%%% Plot measured and predicted relating to each kernel
% plot_areas = [1,2,3,1,4];
plot_areas = [2,2,2,2,2];
plot_conditions = {unique(trial_contrastside_allcat),[-1,1],1:3,[1,3],0:1};
plot_conditions_compare = { ...
    trial_contrastside_allcat_early, ...  
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
        curr_data_mean = nanmean(mua_allcat_move(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Str ' num2str(area_idx) ' Measured']);
    axis tight
    line([0,0],ylim,'color','k');       
    
    p2 = subplot(4,length(plot_areas),curr_area+length(plot_areas)*1);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(fluor_kernel_roi_bw_move(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['BW kernel ROI']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p3 = subplot(4,length(plot_areas),curr_area+length(plot_areas)*2);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(fluor_kernel_roi_bw_hemidiff_move(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['BW hemidiff kernel ROI']);
    axis tight
    line([0,0],ylim,'color','k');
    
    p4 = subplot(4,length(plot_areas),curr_area+length(plot_areas)*3);
    hold on; set(gca,'ColorOrder',plot_cols{curr_area});
    for curr_condition = reshape(plot_conditions{curr_area},1,[])
        curr_trials = plot_conditions_compare{curr_area} == curr_condition;
        curr_data_mean = nanmean(fluor_kernel_roi_weighted_move(curr_trials,:,area_idx),1);
        plot(t,curr_data_mean,'linewidth',2);
    end
    ylabel(['Weighted kernel ROI']);
    axis tight
    line([0,0],ylim,'color','k');
    
    linkaxes([p1,p2,p3,p4],'x');
    
end


%% Striatum: get activity by (correct/incorrect) stim and (visual/zero) velocity (experiment-separated)

plot_depth = 2;
min_trials = 5; % minimum trials per group to keep

% Split data
trials_allcat = size(mua_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_animal;

mua_allcat_exp = mat2cell(mua_allcat,use_split,length(t),n_depths);
mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,use_split,length(t),n_depths);
mua_taskpred_allcat_exp = mat2cell(mua_taskpred_allcat,use_split,length(t),n_depths);

mua_taskpred_reduced_allcat_exp = mat2cell(mua_taskpred_reduced_allcat,use_split,length(t),n_depths, ...
    size(mua_taskpred_reduced_allcat,4));
mua_ctxpred_taskpred_reduced_allcat_exp = mat2cell(mua_ctxpred_taskpred_reduced_allcat,use_split,length(t),n_depths, ...
    size(mua_ctxpred_taskpred_reduced_allcat,4));

wheel_velocity_allcat_exp = mat2cell(wheel_velocity_allcat,use_split,length(t),1);

move_t_exp = mat2cell(move_t,use_split,1);
move_idx_exp = mat2cell(move_idx,use_split,1);
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);

trial_side_allcat_exp = mat2cell(trial_side_allcat,use_split,1);
trial_contrast_allcat_exp = mat2cell(trial_contrast_allcat,use_split,1);
trial_choice_allcat_exp = mat2cell(trial_choice_allcat,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

% Get velocity and bins
% (to use max velocity regardless of final choice)
max_vel = AP_signed_max(wheel_velocity_allcat_move(:,t > 0 & t < 0.5),2);
% (to use summed velocity regardless of final choice)
% max_vel = sum(wheel_velocity_allcat_move(:,t > 0 & t < 0.5),2);
% (to use maximum cumulative velocity regardless of final choice)
% max_vel = AP_signed_max(cumsum(wheel_velocity_allcat_move(:,t > 0 & t < 0.5),2),2);

n_vel_bins = 5;
vel_amp_edges = prctile(abs(max_vel),linspace(0,100,n_vel_bins+1));
% vel_amp_edges = linspace(prctile(abs(max_vel),10),prctile(abs(max_vel),90),n_vel_bins);
vel_amp_edges = sort([vel_amp_edges,-vel_amp_edges]);

vel_amp_bins = discretize(max_vel,vel_amp_edges);
vel_amp_bins(vel_amp_bins == n_vel_bins+1) = NaN;

max_vel_exp = mat2cell(max_vel,use_split,1);
vel_amp_bins_exp = mat2cell(vel_amp_bins,use_split,1);

% Loop through experiments, get context activity
stim_activity_context = nan(10,length(t),2,length(use_split));
move_activity_context = nan(n_vel_bins*2+1,length(t),3,length(use_split));
wheel_context = nan(n_vel_bins*2+1,length(t),3,length(use_split));

for curr_exp = 1:length(use_split)
        
    use_rxn = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} < 1;
       
    % Stim activity
    use_activity = (mua_allcat_exp{curr_exp}(:,:,plot_depth) - mua_taskpred_reduced_allcat_exp{curr_exp}(:,:,plot_depth,1));
%     use_activity = (mua_allcat_exp{curr_exp}(:,:,plot_depth) - mua_taskpred_reduced_allcat_exp{curr_exp}(:,:,plot_depth,1)) - ...
%         (mua_ctxpred_allcat_exp{curr_exp}(:,:,plot_depth) - mua_ctxpred_taskpred_reduced_allcat_exp{curr_exp}(:,:,plot_depth,1));
    
    % Stim contexts
    correct_trials = use_rxn & ...
        trial_side_allcat_exp{curr_exp} == -trial_choice_allcat_exp{curr_exp} & trial_contrast_allcat_exp{curr_exp} > 0;
    incorrect_trials = use_rxn & ...
        trial_side_allcat_exp{curr_exp} == trial_choice_allcat_exp{curr_exp} & trial_contrast_allcat_exp{curr_exp} > 0;
    
    % (if not enough incorrect/zero trials, skip them)
    stims = unique([0.06,0.125,0.25,0.5,1].*[-1;1]);
    stim_bins = trial_contrast_allcat_exp{curr_exp}.*trial_side_allcat_exp{curr_exp};
    stim_bins_idx = discretize(stim_bins,stims);
    correct_trials_n = accumarray(stim_bins_idx(correct_trials),1,[length(stims),1],@sum,0);
    incorrect_trials_n = accumarray(stim_bins_idx(incorrect_trials),1,[length(stims),1],@sum,0);
       
    correct_trials(ismember(stim_bins_idx,find(correct_trials_n < min_trials))) = 0;
    incorrect_trials(ismember(stim_bins_idx,find(incorrect_trials_n < min_trials))) = 0;
    
    % Bin and average by stim/context   
    [correct_stim_used,correct_activity_grp_mean] = ...
        grpstats(use_activity(correct_trials,:),stim_bins(correct_trials),{'gname','nanmean'});
    [incorrect_stim_used,incorrect_activity_grp_mean] = ...
        grpstats(use_activity(incorrect_trials,:),stim_bins(incorrect_trials),{'gname','nanmean'});
    
    correct_stim_used_idx = ismember(stims,cellfun(@str2num,correct_stim_used));
    incorrect_stim_used_idx = ismember(stims,cellfun(@str2num,incorrect_stim_used));
    
    stim_activity_context(correct_stim_used_idx,:,1,curr_exp) = correct_activity_grp_mean;
    stim_activity_context(incorrect_stim_used_idx,:,2,curr_exp) = incorrect_activity_grp_mean;
       
    
    % Move activity
    use_activity = (mua_allcat_exp{curr_exp}(:,:,plot_depth) - mua_taskpred_reduced_allcat_exp{curr_exp}(:,:,plot_depth,2));
%     use_activity = (mua_allcat_exp{curr_exp}(:,:,plot_depth) - mua_taskpred_reduced_allcat_exp{curr_exp}(:,:,plot_depth,2)) - ...
%         (mua_ctxpred_allcat_exp{curr_exp}(:,:,plot_depth) - mua_ctxpred_taskpred_reduced_allcat_exp{curr_exp}(:,:,plot_depth,2));
    
    t_leeway = -t(1);
    leeway_samples = round(t_leeway*(sample_rate));
    wheel_velocity_allcat_exp_move = wheel_velocity_allcat_exp{curr_exp};
    for i = 1:use_split(curr_exp)
        use_activity(i,:,:) = circshift(use_activity(i,:,:),-move_idx_exp{curr_exp}(i)+leeway_samples,2);
        wheel_velocity_allcat_exp_move(i,:,:) = circshift(wheel_velocity_allcat_exp_move(i,:,:),-move_idx_exp{curr_exp}(i)+leeway_samples,2);
    end
    
    % Move contexts
    binned_trials = ~isnan(vel_amp_bins_exp{curr_exp});
    
    vis_correct_trials = binned_trials & use_rxn & ...
        trial_outcome_allcat_exp{curr_exp} == 1 & trial_contrast_allcat_exp{curr_exp} > 0;
    vis_incorrect_trials = binned_trials & use_rxn & ...
        trial_outcome_allcat_exp{curr_exp} == -1 & trial_contrast_allcat_exp{curr_exp} > 0;
    zero_trials = binned_trials & use_rxn & trial_contrast_allcat_exp{curr_exp} == 0;    
    
    % (if not enough incorrect/zero trials, skip them)
    vis_correct_trials_n =  accumarray(vel_amp_bins_exp{curr_exp}(vis_correct_trials),1,[n_vel_bins*2+1,1],@sum,0);
    vis_incorrect_trials_n =  accumarray(vel_amp_bins_exp{curr_exp}(vis_incorrect_trials),1,[n_vel_bins*2+1,1],@sum,0);
    zero_trials_n =  accumarray(vel_amp_bins_exp{curr_exp}(zero_trials),1,[n_vel_bins*2+1,1],@sum,0);
    
    vis_correct_trials(ismember(vel_amp_bins_exp{curr_exp},find(vis_correct_trials_n < min_trials))) = 0;
    vis_incorrect_trials(ismember(vel_amp_bins_exp{curr_exp},find(vis_incorrect_trials_n < min_trials))) = 0;
    zero_trials(ismember(vel_amp_bins_exp{curr_exp},find(zero_trials_n < min_trials))) = 0;

    % Bin and average by velocity/context
    [vis_correct_grps_used,vis_correct_max_vel_grp_mean] = ...
        grpstats(max_vel_exp{curr_exp}(vis_correct_trials),vel_amp_bins_exp{curr_exp}(vis_correct_trials),{'gname','nanmean'});
    vis_correct_grps_used = cellfun(@str2num,vis_correct_grps_used);
    
    [vis_incorrect_grps_used,vis_incorrect_max_vel_grp_mean] = ...
        grpstats(max_vel_exp{curr_exp}(vis_incorrect_trials),vel_amp_bins_exp{curr_exp}(vis_incorrect_trials),{'gname','nanmean'});
    vis_incorrect_grps_used = cellfun(@str2num,vis_incorrect_grps_used);
    
    [zero_grps_used,zero_max_vel_grp_mean] = ...
        grpstats(max_vel_exp{curr_exp}(zero_trials),vel_amp_bins_exp{curr_exp}(zero_trials),{'gname','nanmean'});
    zero_grps_used = cellfun(@str2num,zero_grps_used);
    
    vis_correct_wheel_grp_mean = grpstats(wheel_velocity_allcat_exp_move(vis_correct_trials,:),vel_amp_bins_exp{curr_exp}(vis_correct_trials));
    vis_correct_activity_grp_mean = grpstats(use_activity(vis_correct_trials,:),vel_amp_bins_exp{curr_exp}(vis_correct_trials));
    
    vis_incorrect_wheel_grp_mean = grpstats(wheel_velocity_allcat_exp_move(vis_incorrect_trials,:),vel_amp_bins_exp{curr_exp}(vis_incorrect_trials));
    vis_incorrect_activity_grp_mean = grpstats(use_activity(vis_incorrect_trials,:),vel_amp_bins_exp{curr_exp}(vis_incorrect_trials));
    
    zero_wheel_grp_mean = grpstats(wheel_velocity_allcat_exp_move(zero_trials,:),vel_amp_bins_exp{curr_exp}(zero_trials));
    zero_activity_grp_mean = grpstats(use_activity(zero_trials,:),vel_amp_bins_exp{curr_exp}(zero_trials));
    
    % Store
    move_activity_context(vis_correct_grps_used,:,1,curr_exp) = vis_correct_activity_grp_mean;
    move_activity_context(vis_incorrect_grps_used,:,2,curr_exp) = vis_incorrect_activity_grp_mean;
    move_activity_context(zero_grps_used,:,3,curr_exp) = zero_activity_grp_mean;
    
    wheel_context(vis_correct_grps_used,:,1,curr_exp) = vis_correct_wheel_grp_mean;
    wheel_context(vis_incorrect_grps_used,:,2,curr_exp) = vis_incorrect_wheel_grp_mean;
    wheel_context(zero_grps_used,:,3,curr_exp) = zero_wheel_grp_mean;
     
end

stim_activity_context_mean = nanmean(stim_activity_context,4);
move_activity_context_mean = nanmean(move_activity_context,4);
wheel_context_mean = nanmean(wheel_context,4);

use_stim_t = t > 0 & t < 0.2;
use_move_t = t > -0.1 & t < 0.2;
stim_activity_context_max_t = squeeze(nanmean(stim_activity_context(:,use_stim_t,:,:),2));
move_activity_context_max_t = squeeze(nanmean(move_activity_context(:,use_move_t,:,:),2));

% Plot stim activity
figure;
col = colormap_BlueWhiteRed(5);
col(6,:) = [];

p1 = subplot(1,3,1); hold on
set(gca,'ColorOrder',col)
plot(t,stim_activity_context_mean(:,:,1)','linewidth',2)
xlabel('Time from stim');
ylabel('Activity');
axis tight
line([0,0],ylim,'color','k');

p2 = subplot(1,3,2); hold on
set(gca,'ColorOrder',col)
plot(t,stim_activity_context_mean(:,:,2)','linewidth',2)
xlabel('Time from stim');
ylabel('Activity');
axis tight
line([0,0],ylim,'color','k');

subplot(1,3,3); hold on;
l1 = errorbar(stims,nanmean(stim_activity_context_max_t(:,1,:),3), ...
    AP_sem(stim_activity_context_max_t(:,1,:),3),'k','linewidth',2);
scatter(stims,nanmean(stim_activity_context_max_t(:,1,:),3),80,col, ...
    'Filled','MarkerEdgeColor','k','linewidth',2);

l2 = errorbar(stims,nanmean(stim_activity_context_max_t(:,2,:),3), ...
    AP_sem(stim_activity_context_max_t(:,2,:),3),'color',[0.7,0.7,0.7],'linewidth',2);
scatter(stims,nanmean(stim_activity_context_max_t(:,2,:),3),80,col, ...
    'Filled','MarkerEdgeColor',[0.7,0.7,0.7],'linewidth',2);
xlabel('Stim');
ylabel('Max activity');

legend([l1,l2],{'Correct','Incorrect'});
axis square;

linkaxes([p1,p2]);


% Plot move activity 
figure;
col = [linspace(0.2,1,n_vel_bins)',0.1*ones(n_vel_bins,1),linspace(0.2,1,n_vel_bins)'; ...
    1,0,0;
    0.1*ones(n_vel_bins,1),linspace(1,0.2,n_vel_bins)',0.1*ones(n_vel_bins,1)];

subplot(3,3,1); hold on
set(gca,'ColorOrder',col)
plot(t,wheel_context_mean(:,:,1)','linewidth',2)
xlabel('Time from move');
ylabel('Wheel velocity');
axis tight
line([0,0],ylim,'color','k');

p1 = subplot(3,3,2); hold on
set(gca,'ColorOrder',col)
plot(t,move_activity_context_mean(:,:,1)','linewidth',2)
xlabel('Time from move');
ylabel('Activity');
axis tight
line([0,0],ylim,'color','k');

subplot(3,3,4); hold on
set(gca,'ColorOrder',col)
plot(t,wheel_context_mean(:,:,2)'','linewidth',2)
xlabel('Time from move');
ylabel('Wheel velocity');
axis tight
line([0,0],ylim,'color','k');

p2 = subplot(3,3,5); hold on
set(gca,'ColorOrder',col)
plot(t,move_activity_context_mean(:,:,2)','linewidth',2)
xlabel('Time from move');
ylabel('Activity');
axis tight
line([0,0],ylim,'color','k');

subplot(3,3,7); hold on
set(gca,'ColorOrder',col)
plot(t,wheel_context_mean(:,:,3)'','linewidth',2)
xlabel('Time from move');
ylabel('Wheel velocity');
axis tight
line([0,0],ylim,'color','k');

p3 = subplot(3,3,8); hold on
set(gca,'ColorOrder',col)
plot(t,move_activity_context_mean(:,:,3)','linewidth',2)
xlabel('Time from move');
ylabel('Activity');
axis tight
line([0,0],ylim,'color','k');


subplot(1,3,3); hold on;
l1 = errorbar(1:(n_vel_bins)*2+1,nanmean(move_activity_context_max_t(:,1,:),3), ...
    AP_sem(move_activity_context_max_t(:,1,:),3),'k','linewidth',2);
scatter(1:(n_vel_bins)*2+1,nanmean(move_activity_context_max_t(:,1,:),3),80,col, ...
    'Filled','MarkerEdgeColor','k','linewidth',2);

l2 = errorbar(1:(n_vel_bins)*2+1,nanmean(move_activity_context_max_t(:,2,:),3), ...
    AP_sem(move_activity_context_max_t(:,2,:),3),'color',[0.7,0,0],'linewidth',2);
scatter(1:(n_vel_bins)*2+1,nanmean(move_activity_context_max_t(:,2,:),3),80,col, ...
    'Filled','MarkerEdgeColor',[0.7,0,0],'linewidth',2);

l3 = errorbar(1:(n_vel_bins)*2+1,nanmean(move_activity_context_max_t(:,3,:),3), ...
    AP_sem(move_activity_context_max_t(:,3,:),3),'color',[0.7,0.7,0.7],'linewidth',2);
scatter(1:(n_vel_bins)*2+1,nanmean(move_activity_context_max_t(:,3,:),3),80,col, ...
    'Filled','MarkerEdgeColor',[0.7,0.7,0.7],'linewidth',2);

xlabel('Max velocity');
ylabel('Max activity');

legend([l1,l2,l3],{'Visual correct','Visual incorrect','Zero-contrast'});
axis square;

linkaxes([p1,p2,p3]);


%% Cortex: get activity by (correct/incorrect) stim and (visual/zero) velocity (experiment-separated)

plot_roi = 7;
min_trials = 5; % minimum trials per group to keep

% Split data
trials_all = size(mua_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_animal;

fluor_roi_deconv_exp = mat2cell(fluor_roi_deconv,use_split,length(t),n_rois);
fluor_roi_deconv_taskpred_reduced_exp = mat2cell(fluor_roi_taskpred_reduced,use_split,length(t),n_rois, ...
    size(fluor_roi_taskpred_reduced,4));

wheel_velocity_allcat_exp = mat2cell(wheel_velocity_allcat,use_split,length(t),1);

move_t_exp = mat2cell(move_t,use_split,1);
move_idx_exp = mat2cell(move_idx,use_split,1);
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);

trial_side_allcat_exp = mat2cell(trial_side_allcat,use_split,1);
trial_contrast_allcat_exp = mat2cell(trial_contrast_allcat,use_split,1);
trial_choice_allcat_exp = mat2cell(trial_choice_allcat,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

% Get velocity and bins
% (to use max velocity regardless of final choice)
max_vel = AP_signed_max(wheel_velocity_allcat_move(:,t > 0 & t < 0.5),2);
% (to use summed velocity regardless of final choice)
% max_vel = sum(wheel_velocity_allcat_move(:,t > 0 & t < 0.5),2);
% (to use maximum cumulative velocity regardless of final choice)
% max_vel = AP_signed_max(cumsum(wheel_velocity_allcat_move(:,t > 0 & t < 0.5),2),2);

n_vel_bins = 5;
vel_amp_edges = prctile(abs(max_vel),linspace(0,100,n_vel_bins+1));
% vel_amp_edges = linspace(prctile(abs(max_vel),10),prctile(abs(max_vel),90),n_vel_bins);
vel_amp_edges = sort([vel_amp_edges,-vel_amp_edges]);

vel_amp_bins = discretize(max_vel,vel_amp_edges);
vel_amp_bins(vel_amp_bins == n_vel_bins+1) = NaN;

max_vel_exp = mat2cell(max_vel,use_split,1);
vel_amp_bins_exp = mat2cell(vel_amp_bins,use_split,1);

% Loop through experiments, get context activity
stim_activity_context = nan(10,length(t),2,length(use_split));
move_activity_context = nan(n_vel_bins*2+1,length(t),3,length(use_split));
wheel_context = nan(n_vel_bins*2+1,length(t),3,length(use_split));

for curr_exp = 1:length(use_split)
        
    use_rxn = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} < 1;
       
    % Stim activity
    use_activity = fluor_roi_deconv_exp{curr_exp}(:,:,plot_roi) - fluor_roi_deconv_taskpred_reduced_exp{curr_exp}(:,:,plot_roi,1);

    % Stim contexts
    correct_trials = use_rxn & ...
        trial_side_allcat_exp{curr_exp} == -trial_choice_allcat_exp{curr_exp} & trial_contrast_allcat_exp{curr_exp} > 0;
    incorrect_trials = use_rxn & ...
        trial_side_allcat_exp{curr_exp} == trial_choice_allcat_exp{curr_exp} & trial_contrast_allcat_exp{curr_exp} > 0;
    
    % (if not enough incorrect/zero trials, skip them)
    stims = unique([0.06,0.125,0.25,0.5,1].*[-1;1]);
    stim_bins = trial_contrast_allcat_exp{curr_exp}.*trial_side_allcat_exp{curr_exp};
    stim_bins_idx = discretize(stim_bins,stims);
    correct_trials_n = accumarray(stim_bins_idx(correct_trials),1,[length(stims),1],@sum,0);
    incorrect_trials_n = accumarray(stim_bins_idx(incorrect_trials),1,[length(stims),1],@sum,0);
       
    correct_trials(ismember(stim_bins_idx,find(correct_trials_n < min_trials))) = 0;
    incorrect_trials(ismember(stim_bins_idx,find(incorrect_trials_n < min_trials))) = 0;
    
    % Bin and average by stim/context   
    [correct_stim_used,correct_activity_grp_mean] = ...
        grpstats(use_activity(correct_trials,:),stim_bins(correct_trials),{'gname','nanmean'});
    [incorrect_stim_used,incorrect_activity_grp_mean] = ...
        grpstats(use_activity(incorrect_trials,:),stim_bins(incorrect_trials),{'gname','nanmean'});
    
    correct_stim_used_idx = ismember(stims,cellfun(@str2num,correct_stim_used));
    incorrect_stim_used_idx = ismember(stims,cellfun(@str2num,incorrect_stim_used));
    
    stim_activity_context(correct_stim_used_idx,:,1,curr_exp) = correct_activity_grp_mean;
    stim_activity_context(incorrect_stim_used_idx,:,2,curr_exp) = incorrect_activity_grp_mean;
       
    
    % Move activity
    use_activity = fluor_roi_deconv_exp{curr_exp}(:,:,plot_roi) - fluor_roi_deconv_taskpred_reduced_exp{curr_exp}(:,:,plot_roi,4);

    t_leeway = -t(1);
    leeway_samples = round(t_leeway*(sample_rate));
    wheel_velocity_allcat_exp_move = wheel_velocity_allcat_exp{curr_exp};
    for i = 1:use_split(curr_exp)
        use_activity(i,:,:) = circshift(use_activity(i,:,:),-move_idx_exp{curr_exp}(i)+leeway_samples,2);
        wheel_velocity_allcat_exp_move(i,:,:) = circshift(wheel_velocity_allcat_exp_move(i,:,:),-move_idx_exp{curr_exp}(i)+leeway_samples,2);
    end
    
    % Move contexts
    binned_trials = ~isnan(vel_amp_bins_exp{curr_exp});
    
    vis_correct_trials = binned_trials & use_rxn & ...
        trial_outcome_allcat_exp{curr_exp} == 1 & trial_contrast_allcat_exp{curr_exp} > 0;
    vis_incorrect_trials = binned_trials & use_rxn & ...
        trial_outcome_allcat_exp{curr_exp} == -1 & trial_contrast_allcat_exp{curr_exp} > 0;
    zero_trials = binned_trials & use_rxn & trial_contrast_allcat_exp{curr_exp} == 0;    
    
    % (if not enough incorrect/zero trials, skip them)
    vis_correct_trials_n =  accumarray(vel_amp_bins_exp{curr_exp}(vis_correct_trials),1,[n_vel_bins*2+1,1],@sum,0);
    vis_incorrect_trials_n =  accumarray(vel_amp_bins_exp{curr_exp}(vis_incorrect_trials),1,[n_vel_bins*2+1,1],@sum,0);
    zero_trials_n =  accumarray(vel_amp_bins_exp{curr_exp}(zero_trials),1,[n_vel_bins*2+1,1],@sum,0);
    
    vis_correct_trials(ismember(vel_amp_bins_exp{curr_exp},find(vis_correct_trials_n < min_trials))) = 0;
    vis_incorrect_trials(ismember(vel_amp_bins_exp{curr_exp},find(vis_incorrect_trials_n < min_trials))) = 0;
    zero_trials(ismember(vel_amp_bins_exp{curr_exp},find(zero_trials_n < min_trials))) = 0;

    % Bin and average by velocity/context
    [vis_correct_grps_used,vis_correct_max_vel_grp_mean] = ...
        grpstats(max_vel_exp{curr_exp}(vis_correct_trials),vel_amp_bins_exp{curr_exp}(vis_correct_trials),{'gname','nanmean'});
    vis_correct_grps_used = cellfun(@str2num,vis_correct_grps_used);
    
    [vis_incorrect_grps_used,vis_incorrect_max_vel_grp_mean] = ...
        grpstats(max_vel_exp{curr_exp}(vis_incorrect_trials),vel_amp_bins_exp{curr_exp}(vis_incorrect_trials),{'gname','nanmean'});
    vis_incorrect_grps_used = cellfun(@str2num,vis_incorrect_grps_used);
    
    [zero_grps_used,zero_max_vel_grp_mean] = ...
        grpstats(max_vel_exp{curr_exp}(zero_trials),vel_amp_bins_exp{curr_exp}(zero_trials),{'gname','nanmean'});
    zero_grps_used = cellfun(@str2num,zero_grps_used);
    
    vis_correct_wheel_grp_mean = grpstats(wheel_velocity_allcat_exp_move(vis_correct_trials,:),vel_amp_bins_exp{curr_exp}(vis_correct_trials));
    vis_correct_activity_grp_mean = grpstats(use_activity(vis_correct_trials,:),vel_amp_bins_exp{curr_exp}(vis_correct_trials));
    
    vis_incorrect_wheel_grp_mean = grpstats(wheel_velocity_allcat_exp_move(vis_incorrect_trials,:),vel_amp_bins_exp{curr_exp}(vis_incorrect_trials));
    vis_incorrect_activity_grp_mean = grpstats(use_activity(vis_incorrect_trials,:),vel_amp_bins_exp{curr_exp}(vis_incorrect_trials));
    
    zero_wheel_grp_mean = grpstats(wheel_velocity_allcat_exp_move(zero_trials,:),vel_amp_bins_exp{curr_exp}(zero_trials));
    zero_activity_grp_mean = grpstats(use_activity(zero_trials,:),vel_amp_bins_exp{curr_exp}(zero_trials));
    
    % Store
    move_activity_context(vis_correct_grps_used,:,1,curr_exp) = vis_correct_activity_grp_mean;
    move_activity_context(vis_incorrect_grps_used,:,2,curr_exp) = vis_incorrect_activity_grp_mean;
    move_activity_context(zero_grps_used,:,3,curr_exp) = zero_activity_grp_mean;
    
    wheel_context(vis_correct_grps_used,:,1,curr_exp) = vis_correct_wheel_grp_mean;
    wheel_context(vis_incorrect_grps_used,:,2,curr_exp) = vis_incorrect_wheel_grp_mean;
    wheel_context(zero_grps_used,:,3,curr_exp) = zero_wheel_grp_mean;
     
end

stim_activity_context_mean = nanmean(stim_activity_context,4);
move_activity_context_mean = nanmean(move_activity_context,4);
wheel_context_mean = nanmean(wheel_context,4);

use_stim_t = t > 0 & t < 0.2;
use_move_t = t > -0.1 & t < 0.2;
stim_activity_context_max_t = squeeze(nanmean(stim_activity_context(:,use_stim_t,:,:),2));
move_activity_context_max_t = squeeze(nanmean(move_activity_context(:,use_move_t,:,:),2));

% Plot stim activity
figure('Name',wf_roi(plot_roi).area);
col = colormap_BlueWhiteRed(5);
col(6,:) = [];

p1 = subplot(1,3,1); hold on
set(gca,'ColorOrder',col)
plot(t,stim_activity_context_mean(:,:,1)','linewidth',2)
xlabel('Time from stim');
ylabel('Activity');
axis tight
line([0,0],ylim,'color','k');

p2 = subplot(1,3,2); hold on
set(gca,'ColorOrder',col)
plot(t,stim_activity_context_mean(:,:,2)','linewidth',2)
xlabel('Time from stim');
ylabel('Activity');
axis tight
line([0,0],ylim,'color','k');

subplot(1,3,3); hold on;
l1 = errorbar(stims,nanmean(stim_activity_context_max_t(:,1,:),3), ...
    AP_sem(stim_activity_context_max_t(:,1,:),3),'k','linewidth',2);
scatter(stims,nanmean(stim_activity_context_max_t(:,1,:),3),80,col, ...
    'Filled','MarkerEdgeColor','k','linewidth',2);

l2 = errorbar(stims,nanmean(stim_activity_context_max_t(:,2,:),3), ...
    AP_sem(stim_activity_context_max_t(:,2,:),3),'color',[0.7,0.7,0.7],'linewidth',2);
scatter(stims,nanmean(stim_activity_context_max_t(:,2,:),3),80,col, ...
    'Filled','MarkerEdgeColor',[0.7,0.7,0.7],'linewidth',2);
xlabel('Stim');
ylabel('Max activity');

legend([l1,l2],{'Correct','Incorrect'});
axis square;

linkaxes([p1,p2]);


% Plot move activity 
figure('Name',wf_roi(plot_roi).area);
col = [linspace(0.2,1,n_vel_bins)',0.1*ones(n_vel_bins,1),linspace(0.2,1,n_vel_bins)'; ...
    1,0,0;
    0.1*ones(n_vel_bins,1),linspace(1,0.2,n_vel_bins)',0.1*ones(n_vel_bins,1)];

subplot(3,3,1); hold on
set(gca,'ColorOrder',col)
plot(t,wheel_context_mean(:,:,1)','linewidth',2)
xlabel('Time from move');
ylabel('Wheel velocity');
axis tight
line([0,0],ylim,'color','k');

p1 = subplot(3,3,2); hold on
set(gca,'ColorOrder',col)
plot(t,move_activity_context_mean(:,:,1)','linewidth',2)
xlabel('Time from move');
ylabel('Activity');
axis tight
line([0,0],ylim,'color','k');

subplot(3,3,4); hold on
set(gca,'ColorOrder',col)
plot(t,wheel_context_mean(:,:,2)'','linewidth',2)
xlabel('Time from move');
ylabel('Wheel velocity');
axis tight
line([0,0],ylim,'color','k');

p2 = subplot(3,3,5); hold on
set(gca,'ColorOrder',col)
plot(t,move_activity_context_mean(:,:,2)','linewidth',2)
xlabel('Time from move');
ylabel('Activity');
axis tight
line([0,0],ylim,'color','k');

subplot(3,3,7); hold on
set(gca,'ColorOrder',col)
plot(t,wheel_context_mean(:,:,3)'','linewidth',2)
xlabel('Time from move');
ylabel('Wheel velocity');
axis tight
line([0,0],ylim,'color','k');

p3 = subplot(3,3,8); hold on
set(gca,'ColorOrder',col)
plot(t,move_activity_context_mean(:,:,3)','linewidth',2)
xlabel('Time from move');
ylabel('Activity');
axis tight
line([0,0],ylim,'color','k');


subplot(1,3,3); hold on;
l1 = errorbar(1:(n_vel_bins)*2+1,nanmean(move_activity_context_max_t(:,1,:),3), ...
    AP_sem(move_activity_context_max_t(:,1,:),3),'k','linewidth',2);
scatter(1:(n_vel_bins)*2+1,nanmean(move_activity_context_max_t(:,1,:),3),80,col, ...
    'Filled','MarkerEdgeColor','k','linewidth',2);

l2 = errorbar(1:(n_vel_bins)*2+1,nanmean(move_activity_context_max_t(:,2,:),3), ...
    AP_sem(move_activity_context_max_t(:,2,:),3),'color',[0.7,0,0],'linewidth',2);
scatter(1:(n_vel_bins)*2+1,nanmean(move_activity_context_max_t(:,2,:),3),80,col, ...
    'Filled','MarkerEdgeColor',[0.7,0,0],'linewidth',2);

l3 = errorbar(1:(n_vel_bins)*2+1,nanmean(move_activity_context_max_t(:,3,:),3), ...
    AP_sem(move_activity_context_max_t(:,3,:),3),'color',[0.7,0.7,0.7],'linewidth',2);
scatter(1:(n_vel_bins)*2+1,nanmean(move_activity_context_max_t(:,3,:),3),80,col, ...
    'Filled','MarkerEdgeColor',[0.7,0.7,0.7],'linewidth',2);

xlabel('Max velocity');
ylabel('Max activity');

legend([l1,l2,l3],{'Visual correct','Visual incorrect','Zero-contrast'});
axis square;

linkaxes([p1,p2,p3]);


%% Average task/ctx -> striatum kernels

regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
n_regressors = length(regressor_labels);
t_shifts = {[0,0.5]; ... % stim
    [-0.5,1]; ... % move
    [-0.1,0.5]; ... % go cue
    [-0.5,1]}; % outcome

% regressor_labels = {'Stim','Move onset','Outcome'};
% n_regressors = length(regressor_labels);
% t_shifts = {[0,0.5]; ... % stim
%     [-0.5,1]; ... % move
%     [-0.5,1]}; % outcome

% regressor_labels = {'Stim','Move onset','Move ongoing','Move offset','Go cue','Outcome'};
% n_regressors = length(regressor_labels);
% t_shifts = {[0,0.5]; ... % stim
%     [-0.5,0.1]; ... % move onset
%     [0.1,0.1]; ... % move ongoing
%     [0.1,0.5]; ... % move offset
%     [-0.1,0.5]; ... % go cue
%     [-0.5,1]}; % outcome

sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
    round(x(2)*(sample_rate)),t_shifts,'uni',false);
t_shifts = cellfun(@(x) x/sample_rate,sample_shifts,'uni',false);


figure('Name','Striatum');
p = nan(n_depths,n_regressors);
task_regressors_avg = cell(n_depths,n_regressors);
for curr_regressor = 1:n_regressors
    
    if curr_regressor == 1
        col = colormap_BlueWhiteRed(5);
        col(6,:) = [];
    else
        col = lines(2);
    end
    
    for curr_depth = 1:n_depths
        % Concatenate, normalize by day's std, average
        curr_regressors = cellfun(@(x,mua_std) ...
            x{curr_regressor,curr_depth}./mua_std(curr_depth), ...
            vertcat(mua_taskpred_k_all{:}),vertcat(mua_day_std{:}),'uni',false);

        task_regressors_avg{curr_depth,curr_regressor} = ...
            nanmean(cat(3,curr_regressors{:}),3);
        
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors); 
        hold on; set(gca,'ColorOrder',col);
        plot(t_shifts{curr_regressor},task_regressors_avg{curr_depth,curr_regressor}','linewidth',2);
        title(['Regressor ' num2str(curr_regressor) ' depth ' num2str(curr_depth)])
    end
end
linkaxes(p,'xy');


figure('Name','Cortex-predicted striatum');
p = nan(n_depths,n_regressors);
ctxpred_task_regressors_avg = cell(n_depths,n_regressors);
for curr_regressor = 1:n_regressors
    
    if curr_regressor == 1
        col = colormap_BlueWhiteRed(5);
        col(6,:) = [];
    else
        col = lines(2);
    end
    
    for curr_depth = 1:n_depths
        % Concatenate, normalize by day's std, average
        curr_regressors = cellfun(@(x,mua_std) ...
            x{curr_regressor,curr_depth}./mua_std(curr_depth), ...
            vertcat(mua_ctxpred_taskpred_k_all{:}),vertcat(mua_day_std{:}),'uni',false);

        ctxpred_task_regressors_avg{curr_depth,curr_regressor} = ...
            nanmean(cat(3,curr_regressors{:}),3);
        
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors); 
        hold on; set(gca,'ColorOrder',col);
        plot(t_shifts{curr_regressor},ctxpred_task_regressors_avg{curr_depth,curr_regressor}','linewidth',2);
        title(['Regressor ' num2str(curr_regressor) ' depth ' num2str(curr_depth)])
    end
end
linkaxes(p,'xy');


figure('Name','Striatum - Cortex-predicted striatum');
p = nan(n_depths,n_regressors);
for curr_regressor = 1:n_regressors
    
    if curr_regressor == 1
        col = colormap_BlueWhiteRed(5);
        col(6,:) = [];
    else
        col = lines(2);
    end
    
    for curr_depth = 1:n_depths
        
        curr_task_regressor = task_regressors_avg{curr_depth,curr_regressor};
        curr_ctxpred_task_regressor = ctxpred_task_regressors_avg{curr_depth,curr_regressor};
        
        % get a scaling factor?
%         scale_factor = curr_ctxpred_task_regressor(:)\curr_task_regressor(:);
        
        curr_plot_regressors = curr_task_regressor - curr_ctxpred_task_regressor;
        
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors); 
        hold on; set(gca,'ColorOrder',col);
        plot(t_shifts{curr_regressor},curr_plot_regressors','linewidth',2);
        title(['Regressor ' num2str(curr_regressor) ' depth ' num2str(curr_depth)])
    end
end
linkaxes(p,'xy');


kernel_t = fliplr(-((round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate))/sample_rate));

% Get average ctx->str kernels (flip in time to be relative to event)
ctx_str_k_avg = fliplr(nanmean(cell2mat(permute(vertcat(ctx_str_k_all{:}),[2,3,4,1])),4));
ctx_str_k_px_avg = reshape(svdFrameReconstruct( ...
    U_master(:,:,1:size(ctx_str_k_avg,1)), ...
    reshape(ctx_str_k_avg,size(ctx_str_k_avg,1),[])), ...
    size(U_master,1),size(U_master,2),[],size(ctx_str_k_avg,3));

AP_image_scroll(ctx_str_k_px_avg,kernel_t);
axis image; caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k');
set(gcf,'Name','Ctx->Str');

% Get average ctx->wheel kernels (flip in time to be relative to event)
ctx_wheel_k_avg = fliplr(nanmean(cell2mat(permute(vertcat(ctx_wheel_k_all{:}),[2,3,4,1])),4));
ctx_wheel_k_px_avg = reshape(svdFrameReconstruct( ...
    U_master(:,:,1:size(ctx_wheel_k_avg,1)), ...
    reshape(ctx_wheel_k_avg,size(ctx_wheel_k_avg,1),[])), ...
    size(U_master,1),size(U_master,2),[],size(ctx_wheel_k_avg,3));

AP_image_scroll(ctx_wheel_k_px_avg,kernel_t);
axis image; caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k');
set(gcf,'Name','Ctx->Wheel');

AP_image_scroll(ctx_wheel_k_px_avg-AP_reflect_widefield(ctx_wheel_k_px_avg),kernel_t);
axis image; caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k');
set(gcf,'Name','Ctx->Wheel symmetric');


%% Average task -> cortex kernels

n_regressors = 4;

t_shifts = {[0,0.5]; ... % stim
    [-0.5,1]; ... % move
    [-0.1,0.5]; ... % go cue
    [-0.5,1]}; % outcome

sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
    round(x(2)*(sample_rate)),t_shifts,'uni',false);
t_shifts = cellfun(@(x) x/sample_rate,sample_shifts,'uni',false);
regressor_labels = {'Stim','Move onset','Go cue','Outcome'};

for curr_regressor = 1:n_regressors
    curr_k_cell = cellfun(@(x) x(curr_regressor,:),vertcat(fluor_taskpred_k_all{:}),'uni',false);
    curr_k_cell = vertcat(curr_k_cell{:});
    curr_k = permute(cell2mat(permute(arrayfun(@(x) ...
        nanmean(cat(3,curr_k_cell{:,x}),3),1:n_vs,'uni',false),[1,3,2])),[3,2,1]);
    curr_k_px = cell2mat(permute(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:size(curr_k,1)), ...
        curr_k(:,:,x)),1:size(curr_k,3),'uni',false),[1,3,4,2]));
    AP_image_scroll(curr_k_px,t_shifts{curr_regressor});
    axis image;
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'*RdBu'));
    AP_reference_outline('ccf_aligned','k');
    set(gcf,'Name',regressor_labels{curr_regressor});
end


%% Regression to wheel velocity

kernel_t = (round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate))/sample_rate;

% Get average ctx->str kernels (flip in time to be relative to event)
ctx_str_k_avg = fliplr(nanmean(cell2mat(permute(vertcat(ctx_str_k_all{:}),[2,3,4,1])),4));
ctx_str_k_px_avg = reshape(svdFrameReconstruct( ...
    U_master(:,:,1:size(ctx_str_k_avg,1)), ...
    reshape(ctx_str_k_avg,size(ctx_str_k_avg,1),[])), ...
    size(U_master,1),size(U_master,2),[],size(ctx_str_k_avg,3));

AP_image_scroll(ctx_str_k_px_avg,kernel_t);
axis image; caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k');

% Get average ctx->wheel kernels (flip in time to be relative to event)
ctx_wheel_k_avg = fliplr(nanmean(cell2mat(permute(vertcat(ctx_wheel_k_all{:}),[2,3,4,1])),4));
ctx_wheel_k_px_avg = reshape(svdFrameReconstruct( ...
    U_master(:,:,1:size(ctx_wheel_k_avg,1)), ...
    reshape(ctx_wheel_k_avg,size(ctx_wheel_k_avg,1),[])), ...
    size(U_master,1),size(U_master,2),[],size(ctx_wheel_k_avg,3));

AP_image_scroll(ctx_wheel_k_px_avg,kernel_t);
axis image; caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k');

AP_image_scroll(ctx_wheel_k_px_avg-AP_reflect_widefield(ctx_wheel_k_px_avg),kernel_t);
axis image; caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k');

% Plot binned measured/cortex-predicted wheel velocity
wheel_ctxpred_allcat = cell2mat(vertcat(wheel_ctxpred_all{:}));
figure; hold on;
plot(reshape(wheel_velocity_allcat',[],1));
plot(reshape(wheel_ctxpred_allcat',[],1));

n_vel_bins = 10;
vel_amp_edges = linspace(prctile(abs(reshape(wheel_velocity_allcat',[],1)),0), ...
    prctile(abs(reshape(wheel_velocity_allcat',[],1)),100),n_vel_bins);
vel_amp_edges = sort([vel_amp_edges,-vel_amp_edges]);

vel_amp_bins = discretize(wheel_velocity_allcat,vel_amp_edges);
vel_amp_bins(vel_amp_bins == n_vel_bins) = NaN; 

[vel_grp_mean,vel_grp_sem] = ...
    grpstats(reshape(wheel_velocity_allcat',[],1), ...
    reshape(vel_amp_bins',[],1),{'nanmean','nanstd'});

[predicted_vel_grp_mean,predicted_vel_grp_sem] = ...
    grpstats(reshape(wheel_ctxpred_allcat',[],1), ...
    reshape(vel_amp_bins',[],1),{'nanmean','nanstd'});

figure; hold on;
errorbar(vel_grp_mean,predicted_vel_grp_mean,predicted_vel_grp_sem, ...
    predicted_vel_grp_sem,vel_grp_sem,vel_grp_sem,'k','linewidth',2)
xlabel('Measured wheel velocity');
ylabel('Cortex-predicted wheel velocity');
line(ylim,ylim);
axis tight
line(xlim,[0,0])
line([0,0],ylim)

%  Apply a scalar, which makes it look much better?
wheel_scale = predicted_vel_grp_mean\vel_grp_mean;
wheel_ctxpred_allcat_scaled = wheel_ctxpred_allcat.*wheel_scale;

[predicted_vel_scaled_grp_mean,predicted_vel_scaled_grp_sem] = ...
    grpstats(reshape(wheel_ctxpred_allcat_scaled',[],1), ...
    reshape(vel_amp_bins',[],1),{'nanmean','nanstd'});
errorbar(vel_grp_mean,predicted_vel_scaled_grp_mean,predicted_vel_scaled_grp_sem, ...
    predicted_vel_scaled_grp_sem,vel_grp_sem,vel_grp_sem,'r','linewidth',2)

% Get trial choice by velocity / cortex-predicted velocity
[~,outcome_idx] = max(any(outcome_allcat,3),[],2);
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
wheel_velocity_allcat_outcome = wheel_velocity_allcat;
wheel_ctxpred_allcat_outcome = wheel_ctxpred_allcat;
for i = 1:size(mua_allcat,1)
    wheel_velocity_allcat_outcome(i,:,:) = circshift(wheel_velocity_allcat_outcome(i,:,:),-outcome_idx(i)+leeway_samples,2);
    wheel_ctxpred_allcat_outcome(i,:,:) = circshift(wheel_ctxpred_allcat_outcome(i,:,:),-outcome_idx(i)+leeway_samples,2);
end
    
% Find best t to distinguish (I guess delay between response and outcome)
[~,response_idx] = min(max(wheel_velocity_allcat_outcome(trial_choice_allcat == -1,:),[],1) - ...
    min(wheel_velocity_allcat_outcome(trial_choice_allcat == 1,:),[],1));

response_idx_leeway = response_idx-1:response_idx+1;
velocity_choice = sign(AP_signed_max(wheel_velocity_allcat_outcome(:,response_idx_leeway),2));
velocity_ctxpred_choice = sign(AP_signed_max(wheel_ctxpred_allcat_outcome(:,response_idx_leeway),2));

disp(['Velocity choice prediction: ' num2str(nanmean(velocity_choice == trial_choice_allcat)) ...
    ', Velocity ctxpred choice prediction: ' num2str(nanmean(velocity_ctxpred_choice == trial_choice_allcat))]);


stim_correctchoice = nanmean(trial_choice_allcat(trial_contrast_allcat > 0) == ...
    -trial_side_allcat(trial_contrast_allcat > 0));
velocity_correctchoice = nanmean(velocity_choice == trial_choice_allcat);
velocity_ctxpred_correctchoice = nanmean(velocity_ctxpred_choice == trial_choice_allcat);
velocity_correctchoice_zerocontrast = nanmean(velocity_choice(trial_contrast_allcat == 0) == ...
    trial_choice_allcat(trial_contrast_allcat == 0));
velocity_ctxpred_correctchoice_zerocontrast = nanmean(velocity_ctxpred_choice(trial_contrast_allcat == 0) == ...
    trial_choice_allcat(trial_contrast_allcat == 0));


figure;plot([stim_correctchoice,velocity_correctchoice,velocity_ctxpred_correctchoice ...
    velocity_correctchoice_zerocontrast,velocity_ctxpred_correctchoice_zerocontrast],'linewidth',2);
ylabel('Fraction correct choice prediction');
set(gca,'XTick',1:5,'XTickLabel',{'Stim','Vel','Vel-ctxpred','Vel zero','Vel-ctxpred zero'});


% (was trying to get choice prediction from mua

r = mua_allcat_move - mua_taskpred_reduced_allcat_move(:,:,:,2);

a = linspace(0,10,100);
b = nan(size(a));
for i = 1:length(a)
    mua_choice = -sign(nanmean(r(:,18:20,2),2)-a(i));
    b(i) = nanmean(mua_choice == trial_choice_allcat);
end

figure;plot(a,b);
xlabel('Threshold');
ylabel('Fraction correct');






