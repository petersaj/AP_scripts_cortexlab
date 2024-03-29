% (this is to work on AP_ctx_str_figures and related code)

%% Load in task data

% Load data
% data_fn = 'trial_activity_choiceworld';
% data_fn = 'trial_activity_choiceworld_stimxmove';
% data_fn = 'trial_activity_choiceworld_wfonly';
% data_fn = 'trial_activity_choiceworld_oldregressors';
% exclude_data = true;

% data_fn = 'trial_activity_AP_choiceWorldStimPassive_trained';
% data_fn = 'trial_activity_AP_choiceWorldStimPassive_naive';
% data_fn = 'trial_activity_stimKalatsky_naive';
% data_fn = 'trial_activity_stimKalatsky_trained';
% data_fn = 'trial_activity_AP_lcrGratingPassive_pre_muscimol';
% data_fn = 'trial_activity_AP_lcrGratingPassive_post_muscimol';
% data_fn = 'trial_activity_vanillaChoiceworldNoRepeats_pre_muscimol';
data_fn = 'trial_activity_vanillaChoiceworldNoRepeats_post_muscimol';

exclude_data = false;

AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_recording;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));


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


% Stim/move/outcome rank differences by experiment
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
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
trial_groups = {'Stim','Move onset','Outcome'};
t_groups = {[0.05,0.15],[-0.05,0.05],[0,0.1]};

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
                trial_group = trial_outcome_allcat_exp{curr_exp} == 1;
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

%% Str/Ctx-str rank biases (subset of above)

% Stim/move/outcome rank differences by experiment
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
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
trial_groups = {'Stim','Move onset','Outcome'};
t_groups = {[0.05,0.15],[-0.05,0.05],[0,0.1]};

act_rank_difference = nan(length(t),n_depths,length(use_split),length(trial_groups));
predicted_act_rank_difference = nan(length(t),n_depths,length(use_split),length(trial_groups));

act_rank_difference_trial = nan(n_depths,length(use_split),length(trial_groups));
predicted_act_rank_difference_trial = nan(n_depths,length(use_split),length(trial_groups));

n_shuff = 1000;
act_rank_difference_trial_shuff = nan(n_depths,length(use_split),length(trial_groups),n_shuff);
predicted_act_rank_difference_trial_shuff = nan(n_depths,length(use_split),length(trial_groups),n_shuff);
act_rank_difference_trial_predshuff = nan(n_depths,length(use_split),length(trial_groups),n_shuff);

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
                trial_group_1 = trial_side_allcat_exp{curr_exp} == 1;
                trial_group_2 = trial_side_allcat_exp{curr_exp} == -1';
            case 'Move onset'
                use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} < 0.5;
                trial_group_1 = trial_choice_allcat_exp{curr_exp} == -1;
                trial_group_2 = trial_choice_allcat_exp{curr_exp} == 1';
            case 'Outcome'
                use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} < 0.5;
                trial_group_1 = trial_outcome_allcat_exp{curr_exp} == 1;
                trial_group_2 = trial_outcome_allcat_exp{curr_exp} == -1;
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
p = nan(length(trial_groups),1);
for curr_group = 1:length(trial_groups)
    p(curr_group) = subplot(1,3,curr_group); hold on;
    errorbar(squeeze(nanmean(act_rank_difference_trial(:,:,curr_group),2)), ...
        squeeze(AP_sem(act_rank_difference_trial(:,:,curr_group),2)),'linewidth',2,'color','k');
    errorbar(squeeze(nanmean(predicted_act_rank_difference_trial(:,:,curr_group),2)), ...
        squeeze(AP_sem(predicted_act_rank_difference_trial(:,:,curr_group),2)),'linewidth',2,'color',[0,0.7,0]);
    xlabel('Striatum depth');
    ylabel('Rank difference');
    legend({'Measured','Predicted'});
    title(trial_groups{curr_group});
end
linkaxes(p);

% Get significance from shuffled distribution
trained_predicted_diff_ci = prctile(squeeze(nanmean(act_rank_difference_trial_predshuff,2)),[2.5,97.5],3);
figure; hold on;
set(gca,'ColorOrder',lines(3));
plot(squeeze(nanmean(act_rank_difference_trial-predicted_act_rank_difference_trial,2)),'linewidth',2);
plot(reshape(trained_predicted_diff_ci,n_depths,[]),'linewidth',2,'linestyle','--');
xlabel('Striatum depth');
ylabel('Measured-predicted');




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
% AP_imscroll(ctx_str_k_px_avg,kernel_t);
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
plot_areas = [1,2,3,1,4];
% plot_areas = [2,2,2,2,2];
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
min_trials = 2; % minimum trials per group to keep

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

% (to use mousecam movement signed by choice)
% max_vel = sum(movement_allcat_move(:,t > 0 & t < 0.5),2).*trial_choice_allcat;

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
wheel_bin = nan(n_vel_bins*2+1,3,length(use_split));
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
    
    wheel_bin(vis_correct_grps_used,1,curr_exp) = vis_correct_max_vel_grp_mean;
    wheel_bin(vis_incorrect_grps_used,2,curr_exp) = vis_incorrect_max_vel_grp_mean;
    wheel_bin(zero_grps_used,3,curr_exp) = zero_max_vel_grp_mean;
    
end

stim_activity_context_mean = nanmean(stim_activity_context,4);
move_activity_context_mean = nanmean(move_activity_context,4);
wheel_context_mean = nanmean(wheel_context,4);
wheel_bin_mean = nanmean(wheel_bin,3);

use_stim_t = t > 0.05 & t < 0.15;
use_move_t = t > -0.05 & t < 0.05;
stim_activity_context_max_t = squeeze(nanmean(stim_activity_context(:,use_stim_t,:,:),2));
move_activity_context_max_t = squeeze(nanmean(move_activity_context(:,use_move_t,:,:),2));

% Plot stim activity
figure('Name',['Str ' num2str(plot_depth)]);
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
figure('Name',['Str ' num2str(plot_depth)]);
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
plot(t,wheel_context_mean(:,:,2)','linewidth',2)
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
plot(t,wheel_context_mean(:,:,3)','linewidth',2)
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
l1 = errorbar(wheel_bin_mean(:,1),nanmean(move_activity_context_max_t(:,1,:),3), ...
    AP_sem(move_activity_context_max_t(:,1,:),3),'k','linewidth',2);
scatter(wheel_bin_mean(:,1),nanmean(move_activity_context_max_t(:,1,:),3),80,col, ...
    'Filled','MarkerEdgeColor','k','linewidth',2);

l2 = errorbar(wheel_bin_mean(:,2),nanmean(move_activity_context_max_t(:,2,:),3), ...
    AP_sem(move_activity_context_max_t(:,2,:),3),'color',[0.7,0,0],'linewidth',2);
scatter(wheel_bin_mean(:,2),nanmean(move_activity_context_max_t(:,2,:),3),80,col, ...
    'Filled','MarkerEdgeColor',[0.7,0,0],'linewidth',2);

l3 = errorbar(wheel_bin_mean(:,3),nanmean(move_activity_context_max_t(:,3,:),3), ...
    AP_sem(move_activity_context_max_t(:,3,:),3),'color',[0.7,0.7,0.7],'linewidth',2);
scatter(wheel_bin_mean(:,3),nanmean(move_activity_context_max_t(:,3,:),3),80,col, ...
    'Filled','MarkerEdgeColor',[0.7,0.7,0.7],'linewidth',2);

xlabel('Max velocity');
ylabel('Max activity');

legend([l1,l2,l3],{'Visual correct','Visual incorrect','Zero-contrast'});
axis square;

linkaxes([p1,p2,p3]);


%% Cortex: get activity by (correct/incorrect) stim and (visual/zero) velocity (experiment-separated)

plot_roi = 7;
min_trials = 2; % minimum trials per group to keep

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

% (to use mousecam movement signed by choice)
% max_vel = sum(movement_allcat_move(:,t > 0 & t < 0.5),2).*trial_choice_allcat;

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
wheel_bin = nan(n_vel_bins*2+1,3,length(use_split));
for curr_exp = 1:length(use_split)
    
    use_rxn = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} < 0.5;
    
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
    use_activity = fluor_roi_deconv_exp{curr_exp}(:,:,plot_roi) - fluor_roi_deconv_taskpred_reduced_exp{curr_exp}(:,:,plot_roi,2);
    
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
    
    wheel_bin(vis_correct_grps_used,1,curr_exp) = vis_correct_max_vel_grp_mean;
    wheel_bin(vis_incorrect_grps_used,2,curr_exp) = vis_incorrect_max_vel_grp_mean;
    wheel_bin(zero_grps_used,3,curr_exp) = zero_max_vel_grp_mean;
    
end

stim_activity_context_mean = nanmean(stim_activity_context,4);
move_activity_context_mean = nanmean(move_activity_context,4);
wheel_context_mean = nanmean(wheel_context,4);
wheel_bin_mean = nanmean(wheel_bin,3);

use_stim_t = t > 0.05 & t < 0.1;
use_move_t = t > -0.05 & t < 0.05;
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
l1 = errorbar(wheel_bin_mean(:,1),nanmean(move_activity_context_max_t(:,1,:),3), ...
    AP_sem(move_activity_context_max_t(:,1,:),3),'k','linewidth',2);
scatter(wheel_bin_mean(:,1),nanmean(move_activity_context_max_t(:,1,:),3),80,col, ...
    'Filled','MarkerEdgeColor','k','linewidth',2);

l2 = errorbar(wheel_bin_mean(:,2),nanmean(move_activity_context_max_t(:,2,:),3), ...
    AP_sem(move_activity_context_max_t(:,2,:),3),'color',[0.7,0,0],'linewidth',2);
scatter(wheel_bin_mean(:,2),nanmean(move_activity_context_max_t(:,2,:),3),80,col, ...
    'Filled','MarkerEdgeColor',[0.7,0,0],'linewidth',2);

l3 = errorbar(wheel_bin_mean(:,3),nanmean(move_activity_context_max_t(:,3,:),3), ...
    AP_sem(move_activity_context_max_t(:,3,:),3),'color',[0.7,0.7,0.7],'linewidth',2);
scatter(wheel_bin_mean(:,3),nanmean(move_activity_context_max_t(:,3,:),3),80,col, ...
    'Filled','MarkerEdgeColor',[0.7,0.7,0.7],'linewidth',2);

xlabel('Max velocity');
ylabel('Max activity');

legend([l1,l2,l3],{'Visual correct','Visual incorrect','Zero-contrast'});
axis square;

linkaxes([p1,p2,p3]);


%% Ctx/Str: plot activity by stim (move L/R) and move (stim L/R)

% Activity and areas to plot

use_act = mua_allcat;
use_act_taskreduced = mua_taskpred_reduced_allcat;
use_act_move = mua_allcat_move;
use_act_taskreduced_move = mua_taskpred_reduced_allcat_move;
plot_area = 1;
n_areas = size(use_act,3);

% use_act = fluor_roi_deconv;
% use_act_taskreduced = fluor_roi_taskpred_reduced;
% use_act_move = fluor_roi_deconv_move;
% use_act_taskreduced_move = fluor_roi_taskpred_reduced;
% plot_area = 1;
% n_areas = size(use_act,3);

% Get split index
trials_allcat = size(mua_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_animal;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));

% Set number of shuffles for significance testing
n_shuff = 10000; 

%%% STIM ACTIVITY

% Set time to average activity
use_stim_t = t > 0.05 & t < 0.15;

% Get stim bins for each trial
stims = unique([0.06,0.125,0.25,0.5,1].*[-1;1]);
stim_bins = trial_contrast_allcat.*trial_side_allcat;
stim_bins_idx = discretize(stim_bins,stims);

% Set contexts for activity
stim_contexts = [trial_choice_allcat == -1, ...
    trial_choice_allcat == 1];

% Get stim-isolated activity by trial
stim_activity = use_act - use_act_taskreduced(:,:,:,1);
stim_trial_activity = nanmean(stim_activity(:,use_stim_t,:),2);

% Split activity by animal, stim, and correct/incorrect
stim_bins = unique(trial_contrast_allcat(trial_contrast_allcat ~= 0).*[-1,1]);
stim_trial_activity_split = cell(max(split_idx),length(stim_bins),size(stim_contexts,2),n_areas);
for curr_exp = 1:max(split_idx)
    for curr_bin = 1:length(stim_bins)
        for curr_context = 1:size(stim_contexts,2)               
            for curr_depth = 1:n_areas
                % Get activity of all group trials (exclude NaNs)
                curr_trials = ...
                    split_idx == curr_exp & ...
                    trial_contrast_allcat.*trial_side_allcat == stim_bins(curr_bin) & ...
                    stim_contexts(:,curr_context);             
                curr_activity = stim_trial_activity(curr_trials,:,curr_depth);
                stim_trial_activity_split{curr_exp,curr_bin,curr_context,curr_depth} = ...
                    curr_activity(~isnan(curr_activity));                
            end
            
        end
    end
end

stim_trial_activity_split_mean = cellfun(@nanmean,stim_trial_activity_split);

% Plot stim activity by movement direction
figure('Name',['Area ' num2str(plot_area)]);
subplot(1,2,1); hold on;
line_col = [0.6,0,0.6;0,0.6,0];

dot_col = colormap_BlueWhiteRed(5);
dot_col(6,:) = [];

p = nan(size(stim_contexts,2),1);
for curr_context = 1:size(stim_contexts,2)
    p(curr_context) = ...
        errorbar(stims,nanmean(stim_trial_activity_split_mean(:,:,curr_context,plot_area),1), ...
        AP_sem(stim_trial_activity_split_mean(:,:,curr_context,plot_area),1),'color',line_col(curr_context,:),'linewidth',3);
    scatter(stims,nanmean(stim_trial_activity_split_mean(:,:,curr_context,plot_area),1),80,dot_col, ...
        'Filled','MarkerEdgeColor',line_col(curr_context,:),'linewidth',3);
end
xlabel('Stimulus');
ylabel('Activity');
legend(p,{'Move left','Move right'});

% Get condition difference and compare to shuffle
subplot(1,2,2); hold on;

stim_condition_diff = nanmean(stim_trial_activity_split_mean(:,:,1,plot_area) - ...
    stim_trial_activity_split_mean(:,:,2,plot_area),1);

stim_condition_shuff_diff = nan(max(split_idx),length(stim_bins),n_shuff);
for curr_exp = 1:max(split_idx)
    for curr_bin = 1:length(stim_bins)
        curr_act = stim_trial_activity_split(curr_exp,curr_bin,:,plot_area);
        % Shuffle, split, difference
        curr_act_shuff = mat2cell(AP_shake(repmat(vertcat( ...
            curr_act{:}),1,n_shuff),1),cellfun(@length,curr_act),n_shuff);
        stim_condition_shuff_diff(curr_exp,curr_bin,:) = ...
            nanmean(curr_act_shuff{1},1) - nanmean(curr_act_shuff{2},1);            
    end
end
stim_condition_diff_ci = squeeze(prctile(nanmean(stim_condition_shuff_diff,1),[2.5,97.5],3));

plot(stims,stim_condition_diff,'k','linewidth',2);
plot(stims,stim_condition_diff_ci,'k','linewidth',2,'linestyle','--');
xlabel('Stimulus');
ylabel('Condition difference');


%%% MOVEMENT ACTIVITY

% Set time to average activity
use_move_t = t > -0.05 & t < 0.05;

% Get velocity bins for each trial
max_vel = AP_signed_max(wheel_velocity_allcat_move(:,t > 0 & t < 0.2),2);

% Normalize velocity for each day
max_vel_exp = mat2cell(max_vel,trials_recording,1);
max_vel_norm = cell2mat(cellfun(@(x) x./nanmean(abs(x)),max_vel_exp,'uni',false));

n_vel_bins = 4;
vel_edges = prctile(abs(max_vel_norm),linspace(0,100,n_vel_bins+1));
vel_edges = linspace(min(abs(max_vel_norm)),max(abs(max_vel_norm)),n_vel_bins+1);

vel_edges = sort([vel_edges,-vel_edges]);
vel_centers = vel_edges(1:end-1) + diff(vel_edges)/2;

vel_bins = discretize(max_vel_norm,vel_edges);

% Set contexts for activity
move_contexts = [trial_side_allcat == 1 & trial_contrast_allcat > 0, ...
    trial_side_allcat == -1 & trial_contrast_allcat > 0, ...
    trial_contrast_allcat == 0];

% Get move-isolated activity by trial
move_activity = use_act_move - use_act_taskreduced_move(:,:,:,2);
move_trial_activity = nanmean(move_activity(:,use_move_t,:),2);

% Split activity by animal, velocity, and stim side/zero contrast
move_trial_activity_split = cell(max(split_idx),max(vel_bins),size(move_contexts,2),n_areas);
for curr_exp = 1:max(split_idx)
    for curr_bin = 1:max(vel_bins)
        for curr_context = 1:size(move_contexts,2)               
            for curr_area = 1:n_areas
                % Get activity of all group trials (exclude NaNs)
                curr_trials = ...
                    split_idx == curr_exp & ...
                    vel_bins == curr_bin & ...
                    move_contexts(:,curr_context);              
                curr_activity = move_trial_activity(curr_trials,:,curr_area);
                move_trial_activity_split{curr_exp,curr_bin,curr_context,curr_area} = ...
                    curr_activity(~isnan(curr_activity));                
            end
            
        end
    end
end

move_trial_activity_split_mean = cellfun(@nanmean,move_trial_activity_split);

% Plot move activity by stim side
figure('Name',['Area ' num2str(plot_area)]);
subplot(1,2,1); hold on;
line_col = [0.7,0,0;0,0,0.7;0.5,0.5,0.5];

dot_col = [linspace(0.2,1,n_vel_bins)',0.1*ones(n_vel_bins,1),linspace(0.2,1,n_vel_bins)'; ...
    1,0,0;
    0.1*ones(n_vel_bins,1),linspace(1,0.2,n_vel_bins)',0.1*ones(n_vel_bins,1)];

p = nan(size(move_contexts,2),1);
for curr_context = 1:size(move_contexts,2)
    p(curr_context) = ...
        errorbar(vel_centers,nanmean(move_trial_activity_split_mean(:,:,curr_context,plot_area),1), ...
        AP_sem(move_trial_activity_split_mean(:,:,curr_context,plot_area),1),'color',line_col(curr_context,:),'linewidth',3);
    scatter(vel_centers,nanmean(move_trial_activity_split_mean(:,:,curr_context,plot_area),1),80,dot_col, ...
        'Filled','MarkerEdgeColor',line_col(curr_context,:),'linewidth',3);
end
xlabel('Velocity');
ylabel('Activity');
legend(p,{'Stim right','Stim left','No stim'});
 
% Get condition difference and compare to shuffle
subplot(1,2,2); hold on;

move_condition_diff_1 = nanmean(move_trial_activity_split_mean(:,:,1,plot_area) - ...
    move_trial_activity_split_mean(:,:,2,plot_area),1);
move_condition_diff_2 = nanmean(move_trial_activity_split_mean(:,:,2,plot_area) - ...
    move_trial_activity_split_mean(:,:,3,plot_area),1);
move_condition_diff_3 = nanmean(move_trial_activity_split_mean(:,:,1,plot_area) - ...
    move_trial_activity_split_mean(:,:,3,plot_area),1);

move_condition_shuff_diff_1 = nan(max(split_idx),max(vel_bins),n_shuff);
move_condition_shuff_diff_2 = nan(max(split_idx),max(vel_bins),n_shuff);
move_condition_shuff_diff_3 = nan(max(split_idx),max(vel_bins),n_shuff);
for curr_exp = 1:max(split_idx)
    for curr_bin = 1:max(vel_bins)
        curr_act = move_trial_activity_split(curr_exp,curr_bin,:,plot_area);
        % Shuffle, split, difference
        curr_act_shuff = mat2cell(AP_shake(repmat(vertcat( ...
            curr_act{:}),1,n_shuff),1),cellfun(@length,curr_act),n_shuff);
        move_condition_shuff_diff_1(curr_exp,curr_bin,:) = ...
            nanmean(curr_act_shuff{1},1) - nanmean(curr_act_shuff{2},1);            
        move_condition_shuff_diff_2(curr_exp,curr_bin,:) = ...
            nanmean(curr_act_shuff{2},1) - nanmean(curr_act_shuff{3},1);            
        move_condition_shuff_diff_3(curr_exp,curr_bin,:) = ...
            nanmean(curr_act_shuff{1},1) - nanmean(curr_act_shuff{3},1);            
    end
end
move_condition_diff_1_ci = squeeze(prctile(nanmean(move_condition_shuff_diff_1,1),[2.5,97.5],3));
move_condition_diff_2_ci = squeeze(prctile(nanmean(move_condition_shuff_diff_2,1),[2.5,97.5],3));
move_condition_diff_3_ci = squeeze(prctile(nanmean(move_condition_shuff_diff_3,1),[2.5,97.5],3));

col = lines(3);
plot(vel_centers,move_condition_diff_1,'color',col(1,:),'linewidth',2);
plot(vel_centers,move_condition_diff_1_ci,'color',col(1,:),'linewidth',2,'linestyle','--');

plot(vel_centers,move_condition_diff_2,'color',col(2,:),'linewidth',2);
plot(vel_centers,move_condition_diff_2_ci,'color',col(2,:),'linewidth',2,'linestyle','--');

plot(vel_centers,move_condition_diff_3,'color',col(3,:),'linewidth',2);
plot(vel_centers,move_condition_diff_3_ci,'color',col(3,:),'linewidth',2,'linestyle','--');
xlabel('Velocity');
ylabel('Condition difference');


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

% regressor_labels = {'Stim','Move onset','Move x Stim','Go cue','Outcome'};
% n_regressors = length(regressor_labels);
% t_shifts = {[0,0.5]; ... % stim
%     [-0.5,1]; ... % move
%     [-0.5,1]; ... % move
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
            x{curr_regressor,curr_depth}./((1/sample_rate)*mua_std(curr_depth)), ...
            vertcat(mua_taskpred_k_all{:}),vertcat(mua_norm{:}),'uni',false);
        
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
            x{curr_regressor,curr_depth}./((1/sample_rate)*mua_std(curr_depth)), ...
            vertcat(mua_ctxpred_taskpred_k_all{:}),vertcat(mua_norm{:}),'uni',false);
        
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

AP_imscroll(ctx_str_k_px_avg,kernel_t);
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

AP_imscroll(ctx_wheel_k_px_avg,kernel_t);
axis image; caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k');
set(gcf,'Name','Ctx->Wheel');

AP_imscroll(ctx_wheel_k_px_avg-AP_reflect_widefield(ctx_wheel_k_px_avg),kernel_t);
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
regressor_labels = {'Stim','Move onset','Go cue','Outcome'};

% n_regressors = 5;
% t_shifts = {[0,0.5]; ... % stim
%     [-0.5,1]; ... % move
%     [-0.5,1]; ... % move
%     [-0.1,0.5]; ... % go cue
%     [-0.5,1]}; % outcome
% regressor_labels = {'Stim','Move onset','Move x Stim','Go cue','Outcome'};

sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
    round(x(2)*(sample_rate)),t_shifts,'uni',false);
t_shifts = cellfun(@(x) x/sample_rate,sample_shifts,'uni',false);

regressor_px = cell(n_regressors,1);
for curr_regressor = 1:n_regressors
    curr_k_cell = cellfun(@(x) x(curr_regressor,:),vertcat(fluor_taskpred_k_all{:}),'uni',false);
    curr_k_cell = vertcat(curr_k_cell{:});
    curr_k = permute(cell2mat(permute(arrayfun(@(x) ...
        nanmean(cat(3,curr_k_cell{:,x}),3),1:n_vs,'uni',false),[1,3,2])),[3,2,1]);
    curr_k_px = cell2mat(permute(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:size(curr_k,1)), ...
        curr_k(:,:,x)),1:size(curr_k,3),'uni',false),[1,3,4,2]));
    AP_imscroll(curr_k_px,t_shifts{curr_regressor});
    axis image;
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'*RdBu'));
    AP_reference_outline('ccf_aligned','k');
    set(gcf,'Name',regressor_labels{curr_regressor});
    
    regressor_px{curr_regressor} = curr_k_px;
end

regressor_t_max = cellfun(@(x) squeeze(max(x,[],3)),regressor_px,'uni',false);

figure;
max_subregressors = max(cellfun(@(x) size(x,3),regressor_t_max));
max_c = max(abs(cell2mat(cellfun(@(x) x(:),regressor_t_max,'uni',false))));
for curr_regressor = 1:n_regressors
    for curr_subregressor = 1:size(regressor_t_max{curr_regressor},3)
        subplot(n_regressors,max_subregressors, ...
            curr_subregressor+(curr_regressor-1)*max_subregressors);
        imagesc(regressor_t_max{curr_regressor}(:,:,curr_subregressor));
        AP_reference_outline('ccf_aligned','k');
        axis image off;
        colormap(brewermap([],'Greens'));
        caxis([0,max_c]);
    end
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

AP_imscroll(ctx_str_k_px_avg,kernel_t);
axis image; caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k');

% Get average ctx->wheel kernels (flip in time to be relative to event)
ctx_wheel_k_avg = fliplr(nanmean(cell2mat(permute(vertcat(ctx_wheel_k_all{:}),[2,3,4,1])),4));
ctx_wheel_k_px_avg = reshape(svdFrameReconstruct( ...
    U_master(:,:,1:size(ctx_wheel_k_avg,1)), ...
    reshape(ctx_wheel_k_avg,size(ctx_wheel_k_avg,1),[])), ...
    size(U_master,1),size(U_master,2),[],size(ctx_wheel_k_avg,3));

AP_imscroll(ctx_wheel_k_px_avg,kernel_t);
axis image; caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k');

AP_imscroll(ctx_wheel_k_px_avg-AP_reflect_widefield(ctx_wheel_k_px_avg),kernel_t);
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


%% Std of px (necessary for cell below) (actually crap)
error('This doesn''t make sense for time point analyses')

% In chunks of n frames
chunk_size = 10000;

use_v = reshape(permute(fluor_allcat_deconv,[2,1,3]),[],n_vs)';
px_mean = svdFrameReconstruct(U_master(:,:,1:n_vs),nanmean(use_v,2));

n_frames = size(use_v,2);
frame_chunks = discretize(1:n_frames,linspace(1,n_frames,round(n_frames/chunk_size)));

px_std_sq = zeros(size(U_master,1),size(U_master,2));

for curr_chunk = 1:max(frame_chunks)
    
    curr_im = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
        use_v(:,frame_chunks == curr_chunk));
    px_std_sq = px_std_sq + sum((curr_im - px_mean).^2,3);
    clear curr_im
    
    AP_print_progress_fraction(curr_chunk,max(frame_chunks));
end

px_std = sqrt(px_std_sq./(n_frames-1));


%% Timepoint Ctx/Str correlation (experiment-separated)

plot_depth = 1;


% % (to use fluorescence hemidiff)
% mirror_matrix = reshape(U_master(:,:,1:n_vs),[],n_vs)'* ...
%     reshape(AP_reflect_widefield(U_master(:,:,1:n_vs)),[],n_vs);
% fluor_allcat_deconv_mirror = reshape(transpose( ...
%     mirror_matrix*reshape(fluor_allcat_deconv,[],n_vs)'),size(fluor_allcat_deconv));
% fluor_allcat_deconv_hemidiff = fluor_allcat_deconv - fluor_allcat_deconv_mirror;


% Split data
trials_allcat = size(mua_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_animal;

mua_allcat_exp = mat2cell(mua_allcat,use_split,length(t),n_depths);
mua_taskpred_allcat_exp = mat2cell(mua_taskpred_allcat,use_split,length(t),n_depths);

mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,use_split,length(t),n_depths);
mua_ctxpred_taskpred_allcat_exp = mat2cell(mua_ctxpred_taskpred_allcat,use_split,length(t),n_depths);

mua_taskpred_reduced_allcat_exp = mat2cell(mua_taskpred_reduced_allcat,use_split,length(t),n_depths, ...
    size(mua_taskpred_reduced_allcat,4));
mua_ctxpred_taskpred_reduced_allcat_exp = mat2cell(mua_ctxpred_taskpred_reduced_allcat,use_split,length(t),n_depths, ...
    size(mua_ctxpred_taskpred_reduced_allcat,4));

fluor_allcat_deconv_exp = mat2cell(fluor_allcat_deconv,use_split,length(t),n_vs);
fluor_taskpred_allcat_exp = mat2cell(fluor_taskpred_allcat,use_split,length(t),n_vs);
fluor_taskpred_reduced_allcat_exp = mat2cell(fluor_taskpred_reduced_allcat,use_split,length(t),n_vs, ...
    size(fluor_taskpred_reduced_allcat,4));

fluor_roi_deconv_exp = mat2cell(fluor_roi_deconv,use_split,length(t),n_rois);
fluor_roi_taskpred_exp = mat2cell(fluor_roi_taskpred,use_split,length(t),n_rois);
fluor_roi_taskpred_reduced_exp = mat2cell(fluor_roi_taskpred_reduced,use_split,length(t),n_rois, ...
    size(fluor_roi_taskpred_reduced,4));

trial_side_allcat_exp = mat2cell(trial_side_allcat,use_split,1);
trial_contrast_allcat_exp = mat2cell(trial_contrast_allcat,use_split,1);
trial_choice_allcat_exp = mat2cell(trial_choice_allcat,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

move_t_exp = mat2cell(move_t,use_split,1);
move_idx_exp = mat2cell(move_idx,use_split,1);

act_t_corr = nan(length(t),length(t),n_rois,length(use_split));
act_t_corr_ctxpred = nan(length(t),length(t),n_rois,length(use_split));
for curr_exp = 1:length(use_split)
    
    %     use_trials = true(size(mua_allcat_exp{curr_exp},1),1);
    
    use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} <= 0.5;
    
    %     use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} <= 0.5 & ...
    %         trial_contrast_allcat_exp{curr_exp} > 0 & ...
    %         trial_side_allcat_exp{curr_exp} == 1 & ...
    %         trial_choice_allcat_exp{curr_exp} == -1;
    
    %     use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} <= 0.5 & ...
    %         trial_contrast_allcat_exp{curr_exp} == 0;
    
    %     % Set activity to use
    
    use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    use_act_mua_ctxpred = mua_ctxpred_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    use_act_fluor = fluor_roi_deconv_exp{curr_exp}(use_trials,:,:);
    
    %     use_reduced = 1;
    %     use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_taskpred_reduced_allcat_exp{curr_exp}(use_trials,:,plot_depth,2);
    %     use_act_mua_ctxpred = mua_ctxpred_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_ctxpred_taskpred_reduced_allcat_exp{curr_exp}(use_trials,:,plot_depth,use_reduced);
    %     use_act_fluor = fluor_roi_deconv_exp{curr_exp}(use_trials,:,:) - fluor_roi_taskpred_reduced_exp{curr_exp}(use_trials,:,:,use_reduced);
    
    %     use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_taskpred_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    %     use_act_mua_ctxpred = mua_ctxpred_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_ctxpred_taskpred_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    %     use_act_fluor = fluor_roi_deconv_exp{curr_exp}(use_trials,:,:) - fluor_roi_taskpred_exp{curr_exp}(use_trials,:,:);
    
    %     use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_ctxpred_taskpred_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    %     use_act_mua_ctxpred = zeros(size(use_act_mua));
    %     use_act_fluor = fluor_roi_deconv_exp{curr_exp}(use_trials,:,:) - fluor_roi_taskpred_exp{curr_exp}(use_trials,:,:);
    
    %     use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_ctxpred_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    %     use_act_mua_ctxpred = zeros(size(use_act_mua));
    %     use_act_fluor = fluor_roi_deconv_exp{curr_exp}(use_trials,:,:);
    
    %     use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_ctxpred_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    %     use_act_mua_ctxpred = zeros(size(use_act_mua));
    %     use_act_fluor = fluor_roi_deconv_exp{curr_exp}(use_trials,:,:) - fluor_roi_taskpred_reduced_exp{curr_exp}(use_trials,:,:,1);
    
    %     use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_taskpred_reduced_allcat_exp{curr_exp}(use_trials,:,plot_depth,2);
    %     use_act_mua_ctxpred = mua_ctxpred_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_ctxpred_taskpred_reduced_allcat_exp{curr_exp}(use_trials,:,plot_depth,2);
    %     use_act_fluor = fluor_roi_deconv_exp{curr_exp}(use_trials,:,:) - fluor_roi_taskpred_reduced_exp{curr_exp}(use_trials,:,:,2);
    
%     %      (to align to movement)
%     t_leeway = -t(1);
%     leeway_samples = round(t_leeway*(sample_rate));
%     move_idx_use = move_idx_exp{curr_exp}(use_trials);
%     for i = 1:sum(use_trials)
%         use_act_mua(i,:,:) = circshift(use_act_mua(i,:,:),-move_idx_use(i)+leeway_samples,2);
%         use_act_mua_ctxpred(i,:,:) = circshift(use_act_mua_ctxpred(i,:,:),-move_idx_use(i)+leeway_samples,2);
%         use_act_fluor(i,:,:) = circshift(use_act_fluor(i,:,:),-move_idx_use(i)+leeway_samples,2);
%     end
    
    % Set NaNs to mean to ignore in correlation
    nan_samples = isnan(use_act_mua) | isnan(use_act_mua_ctxpred) | any(isnan(use_act_fluor),3);
    use_act_mua(nan_samples) = AP_index_ans(nan_samples.*nanmean(use_act_mua,1),nan_samples);
    use_act_mua_ctxpred(nan_samples) = AP_index_ans(nan_samples.*nanmean(use_act_mua_ctxpred,1),nan_samples);
    use_act_fluor(repmat(nan_samples,1,1,size(use_act_fluor,3))) = ...
        AP_index_ans(nan_samples.*nanmean(use_act_fluor,1), ...
        repmat(nan_samples,1,1,size(use_act_fluor,3)));
    
    % Get correlation of all time points
    act_t_corr(:,:,:,curr_exp) = cell2mat(arrayfun(@(x) 1-pdist2(use_act_mua', ...
        use_act_fluor(:,:,x)','correlation'),permute(1:n_rois,[1,3,2]),'uni',false));
    act_t_corr_ctxpred(:,:,:,curr_exp) = cell2mat(arrayfun(@(x) 1-pdist2(use_act_mua_ctxpred', ...
        use_act_fluor(:,:,x)','correlation'),permute(1:n_rois,[1,3,2]),'uni',false));
    
end
act_t_corr_mean = nanmean(act_t_corr,4);
act_t_corr_ctxpred_mean = nanmean(act_t_corr_ctxpred,4);

figure;
p = nan(4,n_rois);
for curr_roi = 1:n_rois
    
    p(1,curr_roi) = subplot(4,n_rois,curr_roi);
    imagesc(t,t,act_t_corr_mean(:,:,curr_roi)); axis square;
    caxis([-0.2,0.2])
    colormap(brewermap([],'*RdBu'));
    line(xlim,[0,0],'color','k');
    line([0,0],ylim,'color','k');
    line(xlim,ylim,'color','k')
    title(wf_roi(curr_roi).area);
    
    p(2,curr_roi) = subplot(4,n_rois,n_rois + curr_roi);
    imagesc(t,t,act_t_corr_mean(:,:,curr_roi) - act_t_corr_mean(:,:,curr_roi)');
    axis square;
    caxis([-0.2,0.2])
    colormap(brewermap([],'*RdBu'));
    line(xlim,[0,0],'color','k');
    line([0,0],ylim,'color','k');
    line(xlim,ylim,'color','k');
    
    p(3,curr_roi) = subplot(4,n_rois,n_rois*2 + curr_roi);
    imagesc(t,t,act_t_corr_ctxpred_mean(:,:,curr_roi)); axis square;
    caxis([-0.2,0.2])
    colormap(brewermap([],'*RdBu'));
    line(xlim,[0,0],'color','k');
    line([0,0],ylim,'color','k');
    line(xlim,ylim,'color','k');
    title(wf_roi(curr_roi).area);
    
    p(4,curr_roi) = subplot(4,n_rois,n_rois*3 + curr_roi);
    imagesc(t,t,act_t_corr_ctxpred_mean(:,:,curr_roi) - act_t_corr_ctxpred_mean(:,:,curr_roi)');
    axis square;
    caxis([-0.2,0.2])
    colormap(brewermap([],'*RdBu'));
    line(xlim,[0,0],'color','k');
    line([0,0],ylim,'color','k');
    line(xlim,ylim,'color','k');
    
end
linkaxes(p,'xy');


% Correlation by V's
use_t = t > -0.5 & t < 0.5;

muafluor_corr_px = nan(size(U_master,1),size(U_master,2),sum(use_t),length(use_split));
mua2fluor_corr_px = nan(size(U_master,1),size(U_master,2),sum(use_t),length(use_split));
fluor2mua_corr_px = nan(size(U_master,1),size(U_master,2),sum(use_t),length(use_split));
for curr_exp = 1:length(use_split)
    
    %     use_trials = true(size(mua_allcat_exp{curr_exp},1),1);
    
    use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} <= 0.5;
    
    %     use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} <= 0.5 & ...
    %         trial_contrast_allcat_exp{curr_exp} == 0;
    
    %     use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} <= 0.5 & ...
    %         trial_contrast_allcat_exp{curr_exp} == 0 & ...
    %         trial_choice_allcat_exp{curr_exp} == -1;
    
    %     use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} <= 0.5 & ...
    %         trial_contrast_allcat_exp{curr_exp} > 0 & ...
    %         trial_side_allcat_exp{curr_exp} == 1 & ...
    %         trial_choice_allcat_exp{curr_exp} == -1;
    
    %     use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} <= 0.5 & ...
    %         trial_contrast_allcat_exp{curr_exp} > 0;
    
    % Set activity to use
    
    use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    use_act_fluor = fluor_allcat_deconv_exp{curr_exp}(use_trials,:,:);
    
    %         use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_taskpred_reduced_allcat_exp{curr_exp}(use_trials,:,plot_depth,1);
    %         use_act_fluor = fluor_allcat_deconv_exp{curr_exp}(use_trials,:,:) - fluor_taskpred_reduced_allcat_exp{curr_exp}(use_trials,:,:,1);
    
    %         use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_taskpred_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    %         use_act_fluor = fluor_allcat_deconv_exp{curr_exp}(use_trials,:,:) - fluor_taskpred_allcat_exp{curr_exp}(use_trials,:,:);
    
    %         use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_ctxpred_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    %         use_act_fluor = fluor_allcat_deconv_exp{curr_exp}(use_trials,:,:);
    
    %     use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_ctxpred_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    %     use_act_fluor = fluor_allcat_deconv_exp{curr_exp}(use_trials,:,:) - fluor_taskpred_allcat_exp{curr_exp}(use_trials,:,:);
    
    %     use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_ctxpred_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    %     use_act_fluor = fluor_allcat_deconv_exp{curr_exp}(use_trials,:,:) - fluor_taskpred_reduced_allcat_exp{curr_exp}(use_trials,:,:,1);
    
    %     use_act_mua = (mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - ...
    %         mua_taskpred_reduced_allcat_exp{curr_exp}(use_trials,:,plot_depth,2)) - ...
    %         (mua_ctxpred_allcat_exp{curr_exp}(use_trials,:,plot_depth) - ...
    %         mua_ctxpred_taskpred_reduced_allcat_exp{curr_exp}(use_trials,:,plot_depth,2));
    %     use_act_fluor = fluor_allcat_deconv_exp{curr_exp}(use_trials,:,:);
    
    % (to align to movement)
    t_leeway = -t(1);
    leeway_samples = round(t_leeway*(sample_rate));
    move_idx_use = move_idx_exp{curr_exp}(use_trials);
    for i = 1:sum(use_trials)
        use_act_mua(i,:,:) = circshift(use_act_mua(i,:,:),-move_idx_use(i)+leeway_samples,2);
        use_act_fluor(i,:,:) = circshift(use_act_fluor(i,:,:),-move_idx_use(i)+leeway_samples,2);
    end
    
    % Skip of no MUA in this recording
    if ~any(use_act_mua(:))
        continue
    end
    
    % Set NaNs to mean to ignore in correlation
    nan_samples = isnan(use_act_mua) | any(isnan(use_act_fluor),3);
    use_act_mua(nan_samples) = AP_index_ans(nan_samples.*nanmean(use_act_mua,1),nan_samples);
    use_act_fluor(repmat(nan_samples,1,1,size(use_act_fluor,3))) = ...
        AP_index_ans(nan_samples.*nanmean(use_act_fluor,1), ...
        repmat(nan_samples,1,1,size(use_act_fluor,3)));
    
    % Get std for correlation calculation (t-by-t, otherwise huge)
    px_mean = svdFrameReconstruct(U_master(:,:,1:n_vs),permute(nanmean(use_act_fluor,1),[3,2,1]));
    px_std_sq = zeros(size(U_master,1),size(U_master,2),length(t));
    for curr_t = find(use_t)
        curr_im = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
            permute(use_act_fluor(:,curr_t,:),[3,1,2]));
        px_std_sq(:,:,curr_t) = px_std_sq(:,:,curr_t) + ...
            sum((curr_im - px_mean(:,:,curr_t)).^2,3);
        clear curr_im
        %         AP_print_progress_fraction(curr_t,length(t));
    end
    px_std = sqrt(px_std_sq./(sum(use_trials)-1));
    
    mua_std = nanstd(use_act_mua,[],1);
    
    combined_std = px_std(:,:,use_t).*permute(mua_std(:,use_t),[1,3,2]);
    
    % Get MUA-fluor covariance
    mua_fluor_cov = cell2mat(arrayfun(@(x) ...
        ((use_act_mua(:,use_t)-nanmean(use_act_mua(:,use_t),1))'* ...
        (use_act_fluor(:,use_t,x)-nanmean(use_act_fluor(:,use_t,x),1)))./(sum(use_trials)-1), ...
        permute(1:n_vs,[1,3,2]),'uni',false));
    
    % Get forward and backward covariance diagonals
    on_diag = 0; % (can also use -x:x)
    off_diag = 1:2;
    
    muafluor_cov = cell2mat(arrayfun(@(v) nanmean(cell2mat(arrayfun(@(x) ...
        padarray(diag(mua_fluor_cov(:,:,v),on_diag(x)), ...
        abs(on_diag(x)),nan,'pre'), ...
        1:length(on_diag),'uni',false)),2),1:n_vs,'uni',false));
    
    mua2fluor_cov = cell2mat(arrayfun(@(v) nanmean(cell2mat(arrayfun(@(x) ...
        padarray(diag(mua_fluor_cov(:,:,v),off_diag(x)), ...
        abs(off_diag(x)),nan,'pre'), ...
        1:length(off_diag),'uni',false)),2),1:n_vs,'uni',false));
    
    fluor2mua_cov = cell2mat(arrayfun(@(v) nanmean(cell2mat(arrayfun(@(x) ...
        padarray(diag(mua_fluor_cov(:,:,v),-off_diag(x)), ...
        abs(off_diag(x)),nan,'pre'), ...
        1:length(off_diag),'uni',false)),2),1:n_vs,'uni',false));
    
    % Convert into pixels, get correlation
    muafluor_corr_px(:,:,:,curr_exp) = svdFrameReconstruct(U_master(:,:,1:n_vs),muafluor_cov')./combined_std;
    mua2fluor_corr_px(:,:,:,curr_exp) = svdFrameReconstruct(U_master(:,:,1:n_vs),mua2fluor_cov')./combined_std;
    fluor2mua_corr_px(:,:,:,curr_exp) = svdFrameReconstruct(U_master(:,:,1:n_vs),fluor2mua_cov')./combined_std;
    
    %     % Convert into pixels, get correlation (testing soft)
    %     softdenom = nanmean(combined_std(:));
    %     muafluor_corr_px(:,:,:,curr_exp) = svdFrameReconstruct(U_master(:,:,1:n_vs),muafluor_cov')./(combined_std+softdenom);
    %     mua2fluor_corr_px(:,:,:,curr_exp) = svdFrameReconstruct(U_master(:,:,1:n_vs),mua2fluor_cov')./(combined_std+softdenom);
    %     fluor2mua_corr_px(:,:,:,curr_exp) = svdFrameReconstruct(U_master(:,:,1:n_vs),fluor2mua_cov')./(combined_std+softdenom);
    
    AP_print_progress_fraction(curr_exp,length(use_split));
end

muafluor_corr_px_mean = nanmean(muafluor_corr_px,4);
mua2fluor_corr_px_mean = nanmean(mua2fluor_corr_px,4);
fluor2mua_corr_px_mean = nanmean(fluor2mua_corr_px,4);

AP_imscroll([muafluor_corr_px_mean,fluor2mua_corr_px_mean,mua2fluor_corr_px_mean],t(use_t));
axis image
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k',[],[size(U_master,1),size(U_master,2),1,3]);



%%%%% I THINK I'M DOING CORR WRONG? taken from px corr viewer:

use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth);
use_act_fluor = fluor_allcat_deconv_exp{curr_exp}(use_trials,:,:);

test_t = 20;
test_c = [squeeze(use_act_fluor(:,test_t,:))';use_act_mua(:,test_t)'];

Ur = reshape(U_master(:,:,1:n_vs), size(U_master,1)*size(U_master,2),n_vs);
covV = cov(test_c');
varP = dot((Ur*covV(1:n_vs,:))', Ur');
varS = covV(end,end);


covP = Ur(pixelInd,:)*covV(1:n_vs,1:n_vs)*Ur';
stdPxPy = varS.^0.5 * varP.^0.5;
corrMat = covP./stdPxPy;





%% Timepoint Ctx/Str correlation (passive)

plot_depth = 2;

% Split data
trials_allcat = size(mua_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_animal;

mua_allcat_exp = mat2cell(mua_allcat,use_split,length(t),n_depths);

mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,use_split,length(t),n_depths);



fluor_allcat_deconv_exp = mat2cell(fluor_allcat_deconv,use_split,length(t),n_vs);

fluor_roi_deconv_exp = mat2cell(fluor_roi_deconv,use_split,length(t),n_rois);

stim_exp = mat2cell(D_allcat.stimulus,use_split,1);

wheel_thresh = 0.025;
move_trial = any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > wheel_thresh,2);
move_trial_exp =  mat2cell(move_trial,use_split,1);


act_t_corr = nan(length(t),length(t),n_rois,length(use_split));
act_t_corr_ctxpred = nan(length(t),length(t),n_rois,length(use_split));
for curr_exp = 1:length(use_split)
    
    use_trials = ~move_trial_exp{curr_exp} & ismember(stim_exp{curr_exp},[6:10]);
    
    %     % Set activity to use
    use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    use_act_mua_ctxpred = mua_ctxpred_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    use_act_fluor = fluor_roi_deconv_exp{curr_exp}(use_trials,:,:);
    
    % Set NaNs to mean to ignore in correlation
    nan_samples = isnan(use_act_mua) | isnan(use_act_mua_ctxpred) | any(isnan(use_act_fluor),3);
    use_act_mua(nan_samples) = AP_index_ans(nan_samples.*nanmean(use_act_mua,1),nan_samples);
    use_act_mua_ctxpred(nan_samples) = AP_index_ans(nan_samples.*nanmean(use_act_mua_ctxpred,1),nan_samples);
    use_act_fluor(repmat(nan_samples,1,1,size(use_act_fluor,3))) = ...
        AP_index_ans(nan_samples.*nanmean(use_act_fluor,1), ...
        repmat(nan_samples,1,1,size(use_act_fluor,3)));
    
    % Get correlation of all time points
    act_t_corr(:,:,:,curr_exp) = cell2mat(arrayfun(@(x) 1-pdist2(use_act_mua', ...
        use_act_fluor(:,:,x)','correlation'),permute(1:n_rois,[1,3,2]),'uni',false));
    act_t_corr_ctxpred(:,:,:,curr_exp) = cell2mat(arrayfun(@(x) 1-pdist2(use_act_mua_ctxpred', ...
        use_act_fluor(:,:,x)','correlation'),permute(1:n_rois,[1,3,2]),'uni',false));
    
end
act_t_corr_mean = nanmean(act_t_corr,4);
act_t_corr_ctxpred_mean = nanmean(act_t_corr_ctxpred,4);

figure;
p = nan(4,n_rois);
for curr_roi = 1:n_rois
    
    p(1,curr_roi) = subplot(4,n_rois,curr_roi);
    imagesc(t,t,act_t_corr_mean(:,:,curr_roi)); axis square;
    caxis([-0.2,0.2])
    colormap(brewermap([],'*RdBu'));
    line(xlim,[0,0],'color','k');
    line([0,0],ylim,'color','k');
    line(xlim,ylim,'color','k')
    title(wf_roi(curr_roi).area);
    
    p(2,curr_roi) = subplot(4,n_rois,n_rois + curr_roi);
    imagesc(t,t,act_t_corr_mean(:,:,curr_roi) - act_t_corr_mean(:,:,curr_roi)');
    axis square;
    caxis([-0.2,0.2])
    colormap(brewermap([],'*RdBu'));
    line(xlim,[0,0],'color','k');
    line([0,0],ylim,'color','k');
    line(xlim,ylim,'color','k');
    
    p(3,curr_roi) = subplot(4,n_rois,n_rois*2 + curr_roi);
    imagesc(t,t,act_t_corr_ctxpred_mean(:,:,curr_roi)); axis square;
    caxis([-0.2,0.2])
    colormap(brewermap([],'*RdBu'));
    line(xlim,[0,0],'color','k');
    line([0,0],ylim,'color','k');
    line(xlim,ylim,'color','k');
    title(wf_roi(curr_roi).area);
    
    p(4,curr_roi) = subplot(4,n_rois,n_rois*3 + curr_roi);
    imagesc(t,t,act_t_corr_ctxpred_mean(:,:,curr_roi) - act_t_corr_ctxpred_mean(:,:,curr_roi)');
    axis square;
    caxis([-0.2,0.2])
    colormap(brewermap([],'*RdBu'));
    line(xlim,[0,0],'color','k');
    line([0,0],ylim,'color','k');
    line(xlim,ylim,'color','k');
    
end
linkaxes(p,'xy');


%% Px-Str correlation (as above, but bootstrapped to reduce noise?)

plot_depth = 2;
n_boot = 100;

% Split data
trials_allcat = size(mua_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_animal;

mua_allcat_exp = mat2cell(mua_allcat,use_split,length(t),n_depths);
mua_taskpred_allcat_exp = mat2cell(mua_taskpred_allcat,use_split,length(t),n_depths);

mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,use_split,length(t),n_depths);
mua_ctxpred_taskpred_allcat_exp = mat2cell(mua_ctxpred_taskpred_allcat,use_split,length(t),n_depths);

mua_taskpred_reduced_allcat_exp = mat2cell(mua_taskpred_reduced_allcat,use_split,length(t),n_depths, ...
    size(mua_taskpred_reduced_allcat,4));
mua_ctxpred_taskpred_reduced_allcat_exp = mat2cell(mua_ctxpred_taskpred_reduced_allcat,use_split,length(t),n_depths, ...
    size(mua_ctxpred_taskpred_reduced_allcat,4));

fluor_allcat_deconv_exp = mat2cell(fluor_allcat_deconv,use_split,length(t),n_vs);
fluor_taskpred_allcat_exp = mat2cell(fluor_taskpred_allcat,use_split,length(t),n_vs);
fluor_taskpred_reduced_allcat_exp = mat2cell(fluor_taskpred_reduced_allcat,use_split,length(t),n_vs, ...
    size(fluor_taskpred_reduced_allcat,4));

fluor_roi_deconv_exp = mat2cell(fluor_roi_deconv,use_split,length(t),n_rois);
fluor_roi_taskpred_exp = mat2cell(fluor_roi_taskpred,use_split,length(t),n_rois);
fluor_roi_taskpred_reduced_exp = mat2cell(fluor_roi_taskpred_reduced,use_split,length(t),n_rois, ...
    size(fluor_roi_taskpred_reduced,4));

trial_side_allcat_exp = mat2cell(trial_side_allcat,use_split,1);
trial_contrast_allcat_exp = mat2cell(trial_contrast_allcat,use_split,1);
trial_choice_allcat_exp = mat2cell(trial_choice_allcat,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

move_t_exp = mat2cell(move_t,use_split,1);
move_idx_exp = mat2cell(move_idx,use_split,1);


% Correlation by ROIs
act_t_corr = nan(length(t),length(t),n_rois,length(use_split));
act_t_corr_ctxpred = nan(length(t),length(t),n_rois,length(use_split));
for curr_exp = 1:length(use_split)
    
    use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} <= 0.5;
    
    %     use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} <= 0.5 & ...
    %         trial_contrast_allcat_exp{curr_exp} > 0 & ...
    %         trial_side_allcat_exp{curr_exp} == 1 & ...
    %         trial_choice_allcat_exp{curr_exp} == -1;
    
    %         use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} <= 0.5 & ...
    %             trial_contrast_allcat_exp{curr_exp} == 0;
    
    % Set activity to use
    
    %         use_reduced = 2;
    %         use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_taskpred_reduced_allcat_exp{curr_exp}(use_trials,:,plot_depth,2);
    %         use_act_mua_ctxpred = mua_ctxpred_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_ctxpred_taskpred_reduced_allcat_exp{curr_exp}(use_trials,:,plot_depth,use_reduced);
    %         use_act_fluor = fluor_roi_deconv_exp{curr_exp}(use_trials,:,:) - fluor_roi_taskpred_reduced_exp{curr_exp}(use_trials,:,:,use_reduced);
    
    use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_taskpred_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    use_act_mua_ctxpred = mua_ctxpred_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_ctxpred_taskpred_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    use_act_fluor = fluor_roi_deconv_exp{curr_exp}(use_trials,:,:) - fluor_roi_taskpred_exp{curr_exp}(use_trials,:,:);
    
    %     use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_ctxpred_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    %     use_act_mua_ctxpred = zeros(size(use_act_mua));
    %     use_act_fluor = fluor_roi_deconv_exp{curr_exp}(use_trials,:,:);
    
    %     use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_ctxpred_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    %     use_act_mua_ctxpred = zeros(size(use_act_mua));
    %     use_act_fluor = fluor_roi_deconv_exp{curr_exp}(use_trials,:,:) - fluor_roi_taskpred_reduced_exp{curr_exp}(use_trials,:,:,1);
    
    %     use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_taskpred_reduced_allcat_exp{curr_exp}(use_trials,:,plot_depth,2);
    %     use_act_mua_ctxpred = mua_ctxpred_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_ctxpred_taskpred_reduced_allcat_exp{curr_exp}(use_trials,:,plot_depth,2);
    %     use_act_fluor = fluor_roi_deconv_exp{curr_exp}(use_trials,:,:) - fluor_roi_taskpred_reduced_exp{curr_exp}(use_trials,:,:,2);
    
    if isempty(use_act_mua(:))
        continue
    end
    
    % (to align to movement)
    t_leeway = -t(1);
    leeway_samples = round(t_leeway*(sample_rate));
    move_idx_use = move_idx_exp{curr_exp}(use_trials);
    for i = 1:sum(use_trials)
        use_act_mua(i,:,:) = circshift(use_act_mua(i,:,:),-move_idx_use(i)+leeway_samples,2);
        use_act_mua_ctxpred(i,:,:) = circshift(use_act_mua_ctxpred(i,:,:),-move_idx_use(i)+leeway_samples,2);
        use_act_fluor(i,:,:) = circshift(use_act_fluor(i,:,:),-move_idx_use(i)+leeway_samples,2);
    end
    
    % Set NaNs to mean to ignore in correlation
    nan_samples = isnan(use_act_mua) | isnan(use_act_mua_ctxpred) | any(isnan(use_act_fluor),3);
    use_act_mua(nan_samples) = AP_index_ans(nan_samples.*nanmean(use_act_mua,1),nan_samples);
    use_act_mua_ctxpred(nan_samples) = AP_index_ans(nan_samples.*nanmean(use_act_mua_ctxpred,1),nan_samples);
    use_act_fluor(repmat(nan_samples,1,1,size(use_act_fluor,3))) = ...
        AP_index_ans(nan_samples.*nanmean(use_act_fluor,1), ...
        repmat(nan_samples,1,1,size(use_act_fluor,3)));
    
    act_t_corr_boot = nan(length(t),length(t),n_rois,n_boot);
    act_t_corr_ctxpred_boot = nan(length(t),length(t),n_rois,n_boot);
    for curr_boot = 1:n_boot
        
        curr_boot_trials = randi(sum(use_trials),sum(use_trials),1);
        use_act_mua_boot = use_act_mua(curr_boot_trials,:);
        use_act_mua_cxpred_boot = use_act_mua_ctxpred(curr_boot_trials,:);
        use_act_fluor_boot = use_act_fluor(curr_boot_trials,:,:);
        
        % Get correlation of all time points
        act_t_corr_boot(:,:,:,curr_boot) = cell2mat(arrayfun(@(x) 1-pdist2(use_act_mua_boot', ...
            use_act_fluor_boot(:,:,x)','correlation'),permute(1:n_rois,[1,3,2]),'uni',false));
        act_t_corr_ctxpred_boot(:,:,:,curr_boot) = cell2mat(arrayfun(@(x) 1-pdist2(use_act_mua_cxpred_boot', ...
            use_act_fluor_boot(:,:,x)','correlation'),permute(1:n_rois,[1,3,2]),'uni',false));
        
        AP_print_progress_fraction(curr_boot,n_boot);
        
    end
    act_t_corr(:,:,:,curr_exp) = nanmean(act_t_corr_boot,4);
    act_t_corr_ctxpred(:,:,:,curr_exp) = nanmean(act_t_corr_ctxpred_boot,4);
    
end
act_t_corr_mean = nanmean(act_t_corr,4);
act_t_corr_ctxpred_mean = nanmean(act_t_corr_ctxpred,4);

figure;
p = nan(4,n_rois);
for curr_roi = 1:n_rois
    
    p(1,curr_roi) = subplot(4,n_rois,curr_roi);
    imagesc(t,t,act_t_corr_mean(:,:,curr_roi)); axis square;
    caxis([-0.2,0.2])
    colormap(brewermap([],'*RdBu'));
    line(xlim,[0,0],'color','k');
    line([0,0],ylim,'color','k');
    line(xlim,ylim,'color','k')
    title(wf_roi(curr_roi).area);
    
    p(2,curr_roi) = subplot(4,n_rois,n_rois + curr_roi);
    imagesc(t,t,act_t_corr_mean(:,:,curr_roi) - act_t_corr_mean(:,:,curr_roi)');
    axis square;
    caxis([-0.2,0.2])
    colormap(brewermap([],'*RdBu'));
    line(xlim,[0,0],'color','k');
    line([0,0],ylim,'color','k');
    line(xlim,ylim,'color','k')
    
    p(3,curr_roi) = subplot(4,n_rois,n_rois*2 + curr_roi);
    imagesc(t,t,act_t_corr_ctxpred_mean(:,:,curr_roi)); axis square;
    caxis([-0.2,0.2])
    colormap(brewermap([],'*RdBu'));
    line(xlim,[0,0],'color','k');
    line([0,0],ylim,'color','k');
    line(xlim,ylim,'color','k')
    title(wf_roi(curr_roi).area);
    
    p(4,curr_roi) = subplot(4,n_rois,n_rois*3 + curr_roi);
    imagesc(t,t,act_t_corr_ctxpred_mean(:,:,curr_roi) - act_t_corr_ctxpred_mean(:,:,curr_roi)');
    axis square;
    caxis([-0.2,0.2])
    colormap(brewermap([],'*RdBu'));
    line(xlim,[0,0],'color','k');
    line([0,0],ylim,'color','k');
    line(xlim,ylim,'color','k')
    
end
linkaxes(p,'xy');
drawnow;



% Correlation by V's (boostrapping this does essentially nothing?)
use_t = t > -0.1 & t < 0.1;

muafluor_corr_px = nan(size(U_master,1),size(U_master,2),sum(use_t),length(use_split));
mua2fluor_corr_px = nan(size(U_master,1),size(U_master,2),sum(use_t),length(use_split));
fluor2mua_corr_px = nan(size(U_master,1),size(U_master,2),sum(use_t),length(use_split));
for curr_exp = 1:length(use_split)
    
    %     use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} <= 0.5 & ...
    %         trial_contrast_allcat_exp{curr_exp} == 0 & ...
    %         trial_choice_allcat_exp{curr_exp} == -1;
    
    %     use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} <= 0.5 & ...
    %         trial_contrast_allcat_exp{curr_exp} > 0 & ...
    %         trial_side_allcat_exp{curr_exp} == 1 & ...
    %         trial_choice_allcat_exp{curr_exp} == -1;
    
    use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} <= 0.5 & ...
        trial_contrast_allcat_exp{curr_exp} > 0;
    
    % Set activity to use
    %     use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_taskpred_reduced_allcat_exp{curr_exp}(use_trials,:,plot_depth,1);
    %     use_act_fluor = fluor_allcat_deconv_exp{curr_exp}(use_trials,:,:) - fluor_taskpred_reduced_allcat_exp{curr_exp}(use_trials,:,:,1);
    
    %     use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_taskpred_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    %     use_act_fluor = fluor_allcat_deconv_exp{curr_exp}(use_trials,:,:) - fluor_taskpred_allcat_exp{curr_exp}(use_trials,:,:);
    
    use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_ctxpred_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    use_act_fluor = fluor_allcat_deconv_exp{curr_exp}(use_trials,:,:);
    
    %     use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_ctxpred_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    %     use_act_fluor = fluor_allcat_deconv_exp{curr_exp}(use_trials,:,:) - fluor_taskpred_allcat_exp{curr_exp}(use_trials,:,:);
    
    %     use_act_mua = mua_allcat_exp{curr_exp}(use_trials,:,plot_depth) - mua_ctxpred_allcat_exp{curr_exp}(use_trials,:,plot_depth);
    %     use_act_fluor = fluor_allcat_deconv_exp{curr_exp}(use_trials,:,:) - fluor_taskpred_reduced_allcat_exp{curr_exp}(use_trials,:,:,1);
    
    % (to align to movement)
    t_leeway = -t(1);
    leeway_samples = round(t_leeway*(sample_rate));
    move_idx_use = move_idx_exp{curr_exp}(use_trials);
    for i = 1:sum(use_trials)
        use_act_mua(i,:,:) = circshift(use_act_mua(i,:,:),-move_idx_use(i)+leeway_samples,2);
        use_act_fluor(i,:,:) = circshift(use_act_fluor(i,:,:),-move_idx_use(i)+leeway_samples,2);
    end
    
    % Skip of no MUA in this recording
    if ~any(use_act_mua(:))
        continue
    end
    
    % Set NaNs to mean to ignore in correlation
    nan_samples = isnan(use_act_mua) | any(isnan(use_act_fluor),3);
    use_act_mua(nan_samples) = AP_index_ans(nan_samples.*nanmean(use_act_mua,1),nan_samples);
    use_act_fluor(repmat(nan_samples,1,1,size(use_act_fluor,3))) = ...
        AP_index_ans(nan_samples.*nanmean(use_act_fluor,1), ...
        repmat(nan_samples,1,1,size(use_act_fluor,3)));
    
    %%% TESTING: bootstrap to denoise?
    muafluor_corr_px_boot = nan(size(U_master,1),size(U_master,2),sum(use_t),n_boot);
    mua2fluor_corr_px_boot = nan(size(U_master,1),size(U_master,2),sum(use_t),n_boot);
    fluor2mua_corr_px_boot = nan(size(U_master,1),size(U_master,2),sum(use_t),n_boot);
    for curr_boot = 1:n_boot
        
        curr_boot_trials = randi(sum(use_trials),sum(use_trials),1);
        use_act_mua_boot = use_act_mua(curr_boot_trials,:);
        use_act_fluor_boot = use_act_fluor(curr_boot_trials,:,:);
        
        % Get std for correlation calculation (t-by-t, otherwise huge)
        px_mean = svdFrameReconstruct(U_master(:,:,1:n_vs),permute(nanmean(use_act_fluor_boot,1),[3,2,1]));
        px_std_sq = zeros(size(U_master,1),size(U_master,2),length(t));
        for curr_t = find(use_t)
            curr_im = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
                permute(use_act_fluor_boot(:,curr_t,:),[3,1,2]));
            px_std_sq(:,:,curr_t) = px_std_sq(:,:,curr_t) + ...
                sum((curr_im - px_mean(:,:,curr_t)).^2,3);
            clear curr_im
            %         AP_print_progress_fraction(curr_t,length(t));
        end
        px_std = sqrt(px_std_sq./(sum(use_trials)-1));
        
        mua_std = nanstd(use_act_mua_boot,[],1);
        
        combined_std = px_std(:,:,use_t).*permute(mua_std(:,use_t),[1,3,2]);
        
        % Get MUA-fluor covariance
        mua_fluor_cov = cell2mat(arrayfun(@(x) ...
            ((use_act_mua_boot(:,use_t)-nanmean(use_act_mua_boot(:,use_t),1))'* ...
            (use_act_fluor_boot(:,use_t,x)-nanmean(use_act_fluor_boot(:,use_t,x),1)))./(sum(use_trials)-1), ...
            permute(1:n_vs,[1,3,2]),'uni',false));
        
        % Get forward and backward covariance diagonals
        on_diag = 0; % (can also use -x:x)
        off_diag = 1:2;
        
        muafluor_cov = cell2mat(arrayfun(@(v) nanmean(cell2mat(arrayfun(@(x) ...
            padarray(diag(mua_fluor_cov(:,:,v),on_diag(x)), ...
            abs(on_diag(x)),nan,'pre'), ...
            1:length(on_diag),'uni',false)),2),1:n_vs,'uni',false));
        
        mua2fluor_cov = cell2mat(arrayfun(@(v) nanmean(cell2mat(arrayfun(@(x) ...
            padarray(diag(mua_fluor_cov(:,:,v),off_diag(x)), ...
            abs(off_diag(x)),nan,'pre'), ...
            1:length(off_diag),'uni',false)),2),1:n_vs,'uni',false));
        
        fluor2mua_cov = cell2mat(arrayfun(@(v) nanmean(cell2mat(arrayfun(@(x) ...
            padarray(diag(mua_fluor_cov(:,:,v),-off_diag(x)), ...
            abs(off_diag(x)),nan,'pre'), ...
            1:length(off_diag),'uni',false)),2),1:n_vs,'uni',false));
        
        % Convert into pixels, get correlation
        muafluor_corr_px_boot(:,:,:,curr_boot) = svdFrameReconstruct(U_master(:,:,1:n_vs),muafluor_cov')./combined_std;
        mua2fluor_corr_px_boot(:,:,:,curr_boot) = svdFrameReconstruct(U_master(:,:,1:n_vs),mua2fluor_cov')./combined_std;
        fluor2mua_corr_px_boot(:,:,:,curr_boot) = svdFrameReconstruct(U_master(:,:,1:n_vs),fluor2mua_cov')./combined_std;
        
        AP_print_progress_fraction(curr_boot,n_boot)
    end
    
    muafluor_corr_px(:,:,:,curr_exp) = nanmean(muafluor_corr_px_boot,4);
    mua2fluor_corr_px(:,:,:,curr_exp) = nanmean(mua2fluor_corr_px_boot,4);
    fluor2mua_corr_px(:,:,:,curr_exp) = nanmean(fluor2mua_corr_px_boot,4);
    
    %     AP_print_progress_fraction(curr_exp,length(use_split));
end

muafluor_corr_px_mean = nanmean(muafluor_corr_px,4);
mua2fluor_corr_px_mean = nanmean(mua2fluor_corr_px,4);
fluor2mua_corr_px_mean = nanmean(fluor2mua_corr_px,4);

AP_imscroll([muafluor_corr_px_mean,fluor2mua_corr_px_mean,mua2fluor_corr_px_mean],t(use_t));
axis image
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k',[],[size(U_master,1),size(U_master,2),1,3]);


%% ~~~~~~~~ EXPLORATORY ANALYSIS AGAIN (FOR UNITS AND OTHER) ~~~~~~~~~~~~


%% Get PSTH for each unit

% AP_cellraster({stimOn_times,wheel_move_time,reward_t_timeline});

outcome_time = signals_events.responseTimes(1:n_trials)';

% (all trials)
AP_cellraster({stimOn_times,wheel_move_time,outcome_time}, ...
    {trial_conditions(1:n_trials,1).*trial_conditions(1:n_trials,2), ...
    trial_choice(1:n_trials),trial_outcome(1:n_trials)});

% % (to plot subset of trials)
% use_trials = stim_to_move > 0.5;
% AP_cellraster({stimOn_times(use_trials),wheel_move_time(use_trials),outcome_time(use_trials)}, ...
%     {trial_conditions(use_trials,1).*trial_conditions(use_trials,2), ...
%     trial_choice(use_trials),trial_outcome(use_trials)});



%% Regress task to ephys units

%%%% Copied from AP_ctx_str_trial_preprocessing (removed some wf stuff)

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

% Get event-aligned activity
raster_window = [-0.5,2];
upsample_factor = 1;
raster_sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):raster_sample_rate:raster_window(2);

% Get align times
use_align = stimOn_times;
use_align(isnan(use_align)) = 0;

t_peri_event = bsxfun(@plus,use_align,t);
t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];

%%% Trial-align wheel velocity
event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
    wheel_velocity,t_peri_event);

%%% Trial-align facecam movement
event_aligned_movement = interp1(facecam_t(~isnan(facecam_t)), ...
    frame_movement(~isnan(facecam_t)),t_peri_event);

%%% Trial-align outcome (reward page 1, punish page 2)
% (note incorrect outcome imprecise from signals, but looks good)
event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);

event_aligned_outcome(trial_outcome == 1,:,1) = ...
    (cell2mat(arrayfun(@(x) ...
    histcounts(reward_t_timeline,t_bins(x,:)), ...
    find(trial_outcome == 1)','uni',false))) > 0;

event_aligned_outcome(trial_outcome == -1,:,2) = ...
    (cell2mat(arrayfun(@(x) ...
    histcounts(signals_events.responseTimes,t_bins(x,:)), ...
    find(trial_outcome == -1)','uni',false))) > 0;

% Pick trials to keep
use_trials = ...
    trial_outcome ~= 0 & ...
    ~signals_events.repeatTrialValues(1:n_trials) & ...
    stim_to_feedback < 1.5;

% Get behavioural data
D = struct;
D.stimulus = zeros(sum(use_trials),2);

L_trials = signals_events.trialSideValues(1:n_trials) == -1;
R_trials = signals_events.trialSideValues(1:n_trials) == 1;

D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);

D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
D.repeatNum = ones(sum(use_trials),1);

D.outcome = reshape(trial_outcome(use_trials),[],1);

%%% Regress task to cortex/striatum/cortex-predicted striatum

% Get time points to bin
sample_rate = framerate*regression_params.upsample_factor;
time_bins = frame_t(find(frame_t > ...
    regression_params.skip_seconds,1)):1/sample_rate: ...
    frame_t(find(frame_t-frame_t(end) < ...
    -regression_params.skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Get reaction time for building regressors
wheel_thresh = 0.025;
[move_trial,move_idx] = max(abs(event_aligned_wheel) > wheel_thresh,[],2);
move_idx(~move_trial) = NaN;
move_t = nan(size(move_idx));
move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';

% Stim regressors
unique_stim = unique(contrasts(contrasts > 0).*sides');
stim_contrastsides = ...
    signals_events.trialSideValues(1:length(stimOn_times)).* ...
    signals_events.trialContrastValues(1:length(stimOn_times));

stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
for curr_stim = 1:length(unique_stim)
    curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
    stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
end

% Stim move regressors (one for each stim when it starts to move)
stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
for curr_stim = 1:length(unique_stim)
    
    % (find the first photodiode flip after the stim azimuth has
    % moved past a threshold)
    
    curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
        stim_contrastsides == unique_stim(curr_stim));
    
    azimuth_move_threshold = 5; % degrees to consider stim moved
    stim_move_times_signals = ...
        signals_events.stimAzimuthTimes( ...
        abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
    curr_stim_move_times_signals = arrayfun(@(x) ...
        stim_move_times_signals(find(stim_move_times_signals > ...
        curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
    
    curr_stim_move_times_photodiode = arrayfun(@(x) ...
        photodiode_flip_times(find(photodiode_flip_times > ...
        curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));
    
    stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
    
end

% Stim center regressors (one for each stim when it's stopped during reward)
unique_contrasts = unique(contrasts(contrasts > 0));
stim_contrasts = ...
    signals_events.trialContrastValues(1:length(stimOn_times));

stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
for curr_contrast = 1:length(unique_contrasts)
    
    % (find the last photodiode flip before the reward)
    curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ...
        stim_contrasts == unique_contrasts(curr_contrast));
    
    curr_reward_times = arrayfun(@(x) ...
        reward_t_timeline(find(reward_t_timeline > ...
        curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
    
    curr_prereward_photodiode_times = arrayfun(@(x) ...
        photodiode_flip_times(find(photodiode_flip_times < ...
        curr_reward_times(x),1,'last')),1:length(curr_reward_times));
    
    stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
    
end

% Move onset regressors (L/R)
move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
    find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == -1));
move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
    find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == 1));

move_onset_regressors = zeros(2,length(time_bin_centers));
move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);

% Move onset x stim regressors (one for each contrast/side)
move_onset_stim_time_absolute = arrayfun(@(curr_stim) ...
    arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
    find(~isnan(move_idx) & stim_contrastsides' == unique_stim(curr_stim))), ...
    1:length(unique_stim),'uni',false);

move_onset_stim_regressors = zeros(10,length(time_bin_centers));
for curr_stim = 1:length(unique_stim)
    move_onset_stim_regressors(curr_stim,:) = ...
        histcounts(move_onset_stim_time_absolute{curr_stim},time_bins);
end

% Move ongoing regressors (L/R choice for duration of movement)
wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);

move_stopped_t = 0.5;
move_stopped_samples = round(sample_rate*move_stopped_t);
wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.02), ...
    ones(1,move_stopped_samples),'full') > 0;
wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);

move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
    find(time_bin_centers > x): ...
    find(time_bin_centers > x & ~wheel_moving,1), ...
    move_time_L_absolute','uni',false));
move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
    find(time_bin_centers > x): ...
    find(time_bin_centers > x & ~wheel_moving,1), ...
    move_time_R_absolute','uni',false));

move_ongoing_regressors = zeros(2,length(time_bin_centers));
move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
move_ongoing_regressors(2,move_ongoing_R_samples) = 1;

% Go cue regressors - separate for early/late move
% (using signals timing - not precise but looks good)
if length(signals_events.interactiveOnTimes) ~= length(move_t)
    error('Different number of interactive ons and move times')
end

go_cue_regressors = zeros(2,length(time_bin_centers));
go_cue_regressors(1,:) = histcounts( ...
    signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
go_cue_regressors(2,:) = histcounts( ...
    signals_events.interactiveOnTimes(move_t > 0.5),time_bins);

% Outcome regressors
% (using signals timing - not precise but looks good)
outcome_regressors = zeros(2,length(time_bin_centers));

outcome_regressors(1,:) = histcounts( ...
    reward_t_timeline,time_bins);
outcome_regressors(2,:) = histcounts( ...
    signals_events.responseTimes(trial_outcome == -1),time_bins);

% Concatenate selected regressors, set parameters
regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
regressor_labels = {'Stim','Move onset','Go cue','Outcome'};

t_shifts = {[0,0.5]; ... % stim
    [-0.5,1]; ... % move
    [-0.1,0.5]; ... % go cue
    [-0.5,1]}; % outcome

sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
    round(x(2)*(sample_rate)),t_shifts,'uni',false);
lambda = 0;
zs = [false,false];
cvfold = 5;
use_constant = false;
return_constant = false;


%%%%%%% New: regression task -> striatum units

%%% Trial-align striatum units
event_aligned_unit = nan(length(stimOn_times),length(t),max(spike_templates));
t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
for curr_unit = 1:max(spike_templates)
    curr_spikes = spike_times_timeline(spike_templates == curr_unit);
    event_aligned_unit(:,:,curr_unit) = cell2mat(arrayfun(@(x) ...
        histcounts(curr_spikes,t_bins(x,:)), ...
        [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
    AP_print_progress_fraction(curr_unit,max(spike_templates));
end

%%% Binned unit activity across the experiment
binned_spikes = zeros(max(spike_templates),length(time_bin_centers));
for curr_unit = 1:max(spike_templates)
    curr_spike_times = spike_times_timeline(spike_templates == curr_unit);
    binned_spikes(curr_unit,:) = histcounts(curr_spike_times,time_bins);
    AP_print_progress_fraction(curr_unit,max(spike_templates));
end
binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
binned_spikes_std(isnan(binned_spikes_std)) = 0;

%%% Task > striatal unit regression
baseline = nanmean(reshape(event_aligned_unit(:,t < 0,:),[], ...
    size(event_aligned_unit,3))*raster_sample_rate)';
binned_spikes_baselinesubtracted = binned_spikes - baseline;

[unit_taskpred_k,~,unit_expl_var,~] = ...
    AP_regresskernel(regressors,binned_spikes_baselinesubtracted, ...
    sample_shifts,lambda,zs,cvfold,return_constant,use_constant);

% Plot unit explained variance by depth
used_spikes = spike_times_timeline > time_bins(1) & ...
    spike_times_timeline < time_bins(end);
norm_spike_n = mat2gray(log(accumarray(spike_templates(used_spikes),1,[size(templates,1),1])+1));

use_partial_var = 2;
[~,max_regressor_idx] = max(unit_expl_var.partial(:,:,use_partial_var),[],2);
regressor_cols = [01,0,0;1,0,1;0.8,0.7,0.2;0,0,1];

figure;

subplot(1,length(regressors)+2,1,'YDir','reverse'); hold on;
scatter3(norm_spike_n,template_depths, ...
    1:max(spike_templates),20, ...
    regressor_cols(max_regressor_idx,:),'filled')
xlabel('Normalized n spikes');
ylabel('Depth (\mum)');
title({animal,day,'Best regressor'})

subplot(1,length(regressors)+2,2,'YDir','reverse'); hold on;
curr_expl_var = unit_expl_var.total;
curr_expl_var_dotsize = 100*mat2gray(curr_expl_var,[0,1]) + 1;
curr_expl_var_dotsize(isnan(curr_expl_var) | curr_expl_var < -1 | curr_expl_var >= 1) = NaN;
scatter3(norm_spike_n,template_depths, ...
    1:max(spike_templates),curr_expl_var_dotsize, ...
    regressor_cols(max_regressor_idx,:),'filled')
xlabel('Normalized n spikes');
ylabel('Depth (\mum)');
title('Total expl var');

% (plot units size-scaled by explained variance)
for curr_regressor = 1:length(regressors)
    subplot(1,length(regressors)+2,curr_regressor+2,'YDir','reverse');
    hold on;
    curr_expl_var = unit_expl_var.partial(:,curr_regressor,use_partial_var);
    curr_nan = isnan(curr_expl_var) | curr_expl_var < 0 | curr_expl_var >= 1;
    curr_expl_var(curr_nan) = NaN;
    curr_expl_var = mat2gray(curr_expl_var);
    curr_expl_var(curr_nan) = NaN;
    scatter3(norm_spike_n,template_depths, ...
        1:max(spike_templates),curr_expl_var*50+1, ...
        regressor_cols(curr_regressor,:),'filled')
    title(regressor_labels{curr_regressor});
    xlabel('Normalized n spikes');
    ylabel('Depth (\mum)');
end







% %%%%% try just doing a comparison of aligned peaks?
% % event_aligned_unit
%
% % Stim-aligned
% stim_aligned_unit = event_aligned_unit;
%
% % Move-aligned
% [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
%
% realign_idx = move_idx;
% t_leeway = -t(1);
% leeway_samples = round(t_leeway*(sample_rate));
% curr_realigned = event_aligned_unit;
% for curr_unit = 1:max(spike_templates)
%     for curr_trial = 1:size(curr_realigned,1)
%         curr_realigned(curr_trial,:,curr_unit) = ...
%             circshift(curr_realigned(curr_trial,:,curr_unit),-realign_idx(curr_trial)+leeway_samples,2);
%     end
% end
%
% move_aligned_unit = curr_realigned;
%
% % Outcome-aligned
% [~,outcome_idx] = max(any(event_aligned_outcome,3),[],2);
%
% realign_idx = outcome_idx;
% t_leeway = -t(1);
% leeway_samples = round(t_leeway*(sample_rate));
% curr_realigned = event_aligned_unit;
% for curr_unit = 1:max(spike_templates)
%     for curr_trial = 1:size(curr_realigned,1)
%         curr_realigned(curr_trial,:,curr_unit) = ...
%             circshift(curr_realigned(curr_trial,:,curr_unit),-realign_idx(curr_trial)+leeway_samples,2);
%     end
% end
%
% outcome_aligned_unit = curr_realigned;
%
% % Get trial-trial correlation
% smooth_size = 10;
% gw = gausswin(smooth_size,4)';
% smWin = gw./sum(gw);
%
% stim_aligned_unit_smoothed = convn(padarray(stim_aligned_unit, ...
%     [0,floor(length(smWin)/2)],'replicate','both'), ...
%     smWin,'valid');
% move_aligned_unit_smoothed = convn(padarray(move_aligned_unit, ...
%     [0,floor(length(smWin)/2)],'replicate','both'), ...
%     smWin,'valid');
% outcome_aligned_unit_smoothed = convn(padarray(outcome_aligned_unit, ...
%     [0,floor(length(smWin)/2)],'replicate','both'), ...
%     smWin,'valid');
%
% trial_corr = nan(3,max(spike_templates));
% for curr_unit = 1:max(spike_templates)
%     trial_corr(1,curr_unit) = nanmean(AP_itril(squareform(1-pdist(stim_aligned_unit_smoothed(:,:,curr_unit),'correlation')),-1));
%     trial_corr(2,curr_unit) = nanmean(AP_itril(squareform(1-pdist(move_aligned_unit_smoothed(:,:,curr_unit),'correlation')),-1));
%     trial_corr(3,curr_unit) = nanmean(AP_itril(squareform(1-pdist(outcome_aligned_unit_smoothed(:,:,curr_unit),'correlation')),-1));
% end
%
% unit_std = std(reshape(stim_aligned_unit_smoothed,[],max(spike_templates)),[],1);
% trial_corr_norm = trial_corr.^2./sum(trial_corr.^2,1).*unit_std;
%
% % Plot units by depth size-scaled by trial-trial correlation
% norm_spike_n = mat2gray(log(accumarray(spike_templates,1)+1));
%
% figure;
% for curr_align = 1:size(trial_corr,1)
%     subplot(1,size(trial_corr,1),curr_align,'YDir','reverse');
%     hold on;
%     curr_corr = trial_corr_norm(curr_align,:);
%     curr_nan = isnan(curr_corr) | curr_corr < 0 | curr_corr > 1;
%     curr_corr(curr_nan) = NaN;
%     curr_corr = mat2gray(curr_corr);
%     curr_corr(curr_nan) = NaN;
%     scatter3(norm_spike_n,template_depths, ...
%         1:max(spike_templates),curr_corr*50+1,'k','filled')
% end

%% Regress task to ephys units (BATCH)

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

unit_kernel_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'none';
        AP_load_experiment;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Get align times
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        
        %%% Trial-align wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Trial-align facecam movement
        event_aligned_movement = interp1(facecam_t(~isnan(facecam_t)), ...
            frame_movement(~isnan(facecam_t)),t_peri_event);
        
        %%% Trial-align outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_bins(x,:)), ...
            find(trial_outcome == 1),'uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_bins(x,:)), ...
            find(trial_outcome == -1),'uni',false))) > 0;
        
        % Pick trials to keep
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials)' & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials');
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials');
        
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%% Regress task to cortex/striatum/cortex-predicted striatum
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times))'.* ...
            signals_events.trialContrastValues(1:length(stimOn_times))';
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end
        
        % Stim move regressors (one for each stim when it starts to move)
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
            
        end
        
        % Stim center regressors (one for each stim when it's stopped during reward)
        unique_contrasts = unique(contrasts(contrasts > 0));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ...
                abs(stim_contrastsides) == unique_contrasts(curr_contrast));
            
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
            
        end
        
        % Move onset regressors (L/R)
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move onset x stim regressors (one for each contrast/side)
        move_onset_stim_time_absolute = arrayfun(@(curr_stim) ...
            arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & stim_contrastsides == unique_stim(curr_stim))), ...
            1:length(unique_stim),'uni',false);
        
        move_onset_stim_regressors = zeros(10,length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            move_onset_stim_regressors(curr_stim,:) = ...
                histcounts(move_onset_stim_time_absolute{curr_stim},time_bins);
        end
        
        % Move ongoing regressors (L/R choice for duration of movement)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        
        move_stopped_t = 0.5;
        move_stopped_samples = round(sample_rate*move_stopped_t);
        wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.02), ...
            ones(1,move_stopped_samples),'full') > 0;
        wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);
        
        move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_L_absolute','uni',false));
        move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_R_absolute','uni',false));
        
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
        move_ongoing_regressors(2,move_ongoing_R_samples) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        if length(signals_events.interactiveOnTimes) ~= length(move_t)
            error('Different number of interactive ons and move times')
        end
        
        go_cue_regressors = zeros(2,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        outcome_regressors = zeros(2,length(time_bin_centers));
        
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate selected regressors, set parameters
        regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
        
        t_shifts = {[0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [-0.1,0.5]; ... % go cue
            [-0.5,1]}; % outcome
        
        sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;
        
        %%%%%%% Regress task -> units
        
        %%% Trial-align units
        event_aligned_unit = nan(length(stimOn_times),length(t),max(spike_templates));
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_unit = 1:max(spike_templates)
            curr_spikes = spike_times_timeline(spike_templates == curr_unit);
            event_aligned_unit(:,:,curr_unit) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Binned unit activity across the experiment
        binned_spikes = nan(max(spike_templates),length(time_bin_centers));
        for curr_unit = 1:max(spike_templates)
            curr_spike_times = spike_times_timeline(spike_templates == curr_unit);
            if isempty(curr_spike_times)
                continue
            end
            binned_spikes(curr_unit,:) = histcounts(curr_spike_times,time_bins);
        end
        
        %%% Task > striatal units
        baseline = nanmean(reshape(event_aligned_unit(:,t < 0,:),[], ...
            size(event_aligned_unit,3))*raster_sample_rate)';
        binned_spikes_baselinesubtracted = binned_spikes - baseline;
        
        [unit_taskpred_k,binned_spikes_taskpred,unit_expl_var,~] = ...
            AP_regresskernel(regressors,binned_spikes_baselinesubtracted, ...
            sample_shifts,lambda,zs,cvfold,return_constant,use_constant);
        
        % Get normalized total number of spikes
        used_spikes = spike_times_timeline > time_bins(1) & ...
            spike_times_timeline < time_bins(end);
        spike_rate = accumarray(spike_templates(used_spikes),1,[size(templates,1),1])./(time_bins(end)-time_bins(1));
        
        % Package
        unit_kernel_all(curr_animal,curr_day).template_depths = template_depths;
        unit_kernel_all(curr_animal,curr_day).spike_rate = spike_rate;
        unit_kernel_all(curr_animal,curr_day).unit_expl_var_total = unit_expl_var.total;
        unit_kernel_all(curr_animal,curr_day).unit_expl_var_partial = unit_expl_var.partial;
        unit_kernel_all(curr_animal,curr_day).unit_taskpred_k = unit_taskpred_k;
        
        % Clear
        clearvars -except regression_params animals curr_animal animal days curr_day experiments unit_kernel_all
        AP_print_progress_fraction(curr_day,length(experiments));
    end
end

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
save_fn = 'unit_kernel_all.mat';
save([save_path filesep save_fn],'unit_kernel_all');

%% Plot single unit kernel explained variance

% Load the unit kernel results
unit_kernel_dir = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
unit_kernel_fn = 'unit_kernel_all_triaged.mat';
load([unit_kernel_dir filesep unit_kernel_fn]);
regressor_labels = {'Stim','Move','Go cue','Outcome'};
regressor_cols = [01,0,0;1,0,1;0.8,0.7,0.2;0,0,1];

% Get estimation of end of striatum for each recording
ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
ephys_align_fn = ['ephys_depth_align.mat'];
load([ephys_align_path filesep ephys_align_fn]);

% Plot all days of a given animal
curr_animal = 1;
for curr_day = 1:size(unit_kernel_all,2)
    
    if isempty(unit_kernel_all(curr_animal,curr_day).unit_expl_var_total)
        continue
    end
    
    curr_str_depths = ephys_depth_align(curr_animal).str_depth(curr_day,:);
    
    use_partial = 2;
    [~,max_regressor_idx] = max(unit_kernel_all(curr_animal,curr_day).unit_expl_var_partial(:,:,use_partial),[],2);
    
    figure;
    
    subplot(1,length(regressor_labels)+2,1,'YDir','reverse'); hold on;
    curr_expl_var = unit_kernel_all(curr_animal,curr_day).unit_expl_var_total;
    scatter3(log10(unit_kernel_all(curr_animal,curr_day).spike_rate), ...
        unit_kernel_all(curr_animal,curr_day).template_depths, ...
        1:length(unit_kernel_all(curr_animal,curr_day).template_depths),20, ...
        regressor_cols(max_regressor_idx,:),'filled')
    line(xlim,repmat(curr_str_depths(1),1,2),'color','k')
    line(xlim,repmat(curr_str_depths(2),1,2),'color','k')
    xlabel('Normalized n spikes');
    ylabel('Depth (\mum)');
    title({curr_animal,curr_day,'Best regressor'})
    
    subplot(1,length(regressor_labels)+2,2,'YDir','reverse'); hold on;
    curr_expl_var = unit_kernel_all(curr_animal,curr_day).unit_expl_var_total;
    curr_expl_var_dotsize = 100*rescale(curr_expl_var.*(curr_expl_var > 0),0,1) + 1;
    scatter3(log10(unit_kernel_all(curr_animal,curr_day).spike_rate), ...
        unit_kernel_all(curr_animal,curr_day).template_depths, ...
        1:length(unit_kernel_all(curr_animal,curr_day).template_depths),curr_expl_var_dotsize, ...
        regressor_cols(max_regressor_idx,:),'filled')
    line(xlim,repmat(curr_str_depths(1),1,2),'color','k')
    line(xlim,repmat(curr_str_depths(2),1,2),'color','k')
    xlabel('Normalized n spikes');
    ylabel('Depth (\mum)');
    title('Total expl var');
    
    % (plot units size-scaled by explained variance)
    for curr_regressor = 1:length(regressor_labels)
        subplot(1,length(regressor_labels)+2,curr_regressor+2,'YDir','reverse');
        hold on;
        curr_expl_var = unit_kernel_all(curr_animal,curr_day).unit_expl_var_partial(:,curr_regressor,use_partial);
        curr_expl_var_dotsize = 100*rescale(curr_expl_var.*(curr_expl_var > 0),0,1) + 1;
        
        scatter3(log10(unit_kernel_all(curr_animal,curr_day).spike_rate), ...
            unit_kernel_all(curr_animal,curr_day).template_depths, ...
            1:length(unit_kernel_all(curr_animal,curr_day).template_depths), ...
            curr_expl_var_dotsize, ...
            regressor_cols(curr_regressor,:),'filled')
        line(xlim,repmat(curr_str_depths(1),1,2),'color','k')
        line(xlim,repmat(curr_str_depths(2),1,2),'color','k')
        title(regressor_labels{curr_regressor});
        xlabel('Normalized n spikes');
        ylabel('Depth (\mum)');
    end
    
end




%% Plot all unit regression explained variance

unit_kernel_dir = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
% unit_kernel_fn = 'unit_kernel_all.mat';
unit_kernel_fn = 'unit_kernel_all_triaged.mat';
load([unit_kernel_dir filesep unit_kernel_fn]);

regressor_labels = {'Stim','Move','Go cue','Outcome'};
regressor_cols = [01,0,0;1,0,1;0.8,0.7,0.2;0,0,1];

% Concatenate all units (raw distance)
animaldays = reshape(arrayfun(@(x) ~isempty(unit_kernel_all(x).template_depths), ...
    1:numel(unit_kernel_all)),size(unit_kernel_all));

use_partial = 2;

template_depths_cat = AP_index_ans(reshape({unit_kernel_all.template_depths}, ...
    size(unit_kernel_all))',animaldays');
unit_expl_var_total_cat = AP_index_ans(reshape({unit_kernel_all.unit_expl_var_total}, ...
    size(unit_kernel_all))',animaldays');
unit_expl_var_partial_cat = AP_index_ans(reshape({unit_kernel_all.unit_expl_var_partial}, ...
    size(unit_kernel_all))',animaldays');
spike_rate_cat = AP_index_ans(reshape({unit_kernel_all.spike_rate}, ...
    size(unit_kernel_all))',animaldays');

[~,max_regressor_idx_cat] = cellfun(@(x) max(x(:,:,use_partial),[],2),unit_expl_var_partial_cat,'uni',false);

h = figure;

subplot(1,length(regressor_labels)+3,1,'YDir','reverse'); hold on;
scatter(log10(cell2mat(spike_rate_cat)), ...
    cell2mat(template_depths_cat),20, ...
    regressor_cols(cell2mat(max_regressor_idx_cat),:),'filled')
xlabel('log10(spike rate)');
ylabel('Depth (\mum)');
title({'All units','Best regressor'})

subplot(1,length(regressor_labels)+3,2,'YDir','reverse'); hold on;
curr_expl_var = cell2mat(unit_expl_var_total_cat);
curr_expl_var_dotsize = 100*mat2gray(curr_expl_var,[0,1]) + 1;
curr_expl_var_dotsize(isnan(curr_expl_var) | curr_expl_var < -1 | curr_expl_var >= 1) = NaN;
scatter(log10(cell2mat(spike_rate_cat)), ...
    cell2mat(template_depths_cat),curr_expl_var_dotsize, ...
    regressor_cols(cell2mat(max_regressor_idx_cat),:),'filled')
xlabel('log10(spike rate)');
ylabel('Depth (\mum)');
title({'All regressors'})

% (plot units size-scaled by explained variance)
for curr_regressor = 1:length(regressor_labels)
    for curr_norm = 1:2
        
        subplot(2,length(regressor_labels)+3,curr_regressor+2+((length(regressor_labels)+3)*(curr_norm-1)),'YDir','reverse');
        hold on;
        curr_expl_var = cellfun(@(x) x(:,curr_regressor,use_partial),unit_expl_var_partial_cat,'uni',false);
        
        switch curr_norm
            case 1
                % Normalize within experiment
                curr_expl_var_norm = cell2mat(cellfun(@(x) x./nanmax(x(x < 1)),curr_expl_var,'uni',false));
                curr_expl_var_dotsize = 100*mat2gray(curr_expl_var_norm,[0,1]) + 1;
                curr_expl_var_dotsize(isnan(curr_expl_var_norm) | curr_expl_var_norm < -1 | curr_expl_var_norm >= 1) = NaN;
                
                curr_expl_var_norm(curr_expl_var_norm < -1 | curr_expl_var_norm >= 1) = NaN;
                curr_expl_var_norm(curr_expl_var_norm < 0) = 0;
                norm_type = 'experiment normalized';
            case 2
                % Normalize across experiments
                curr_expl_var_norm = cell2mat(curr_expl_var);
                curr_expl_var_norm = curr_expl_var_norm./ ...
                    max(curr_expl_var_norm(curr_expl_var_norm < 1));
                curr_expl_var_dotsize = 100*mat2gray(curr_expl_var_norm,[0,1]) + 1;
                curr_expl_var_dotsize(isnan(curr_expl_var_norm) | curr_expl_var_norm < -1 | curr_expl_var_norm >= 1) = NaN;
                
                curr_expl_var_norm(curr_expl_var_norm < -1 | curr_expl_var_norm >= 1) = NaN;
                curr_expl_var_norm(curr_expl_var_norm < 0) = 0;
                norm_type = 'concat normalized';
        end
        
        scatter(log10(cell2mat(spike_rate_cat)), ...
            cell2mat(template_depths_cat), ...
            curr_expl_var_dotsize, ...
            regressor_cols(curr_regressor,:),'filled')
        title({regressor_labels{curr_regressor},norm_type});
        xlabel('log10(spike rate)');
        ylabel('Depth (\mum)');
        
        subplot(2,length(regressor_labels)+3,curr_norm*(length(regressor_labels)+3),'YDir','reverse');
        hold on;
        depth_bins = linspace(min(cell2mat(template_depths_cat)), ...
            max(cell2mat(template_depths_cat)),round(range(cell2mat(template_depths_cat))/400));
        depth_bin_centers = depth_bins(1:end-1) + diff(depth_bins)./2;
        curr_depth_groups = discretize(cell2mat(template_depths_cat),depth_bins);
        curr_expl_var_norm_depth = accumarray(curr_depth_groups, ...
            curr_expl_var_norm,[length(depth_bin_centers),1],@nanmean);
        plot(rescale(curr_expl_var_norm_depth,0,1),depth_bin_centers,'color',regressor_cols(curr_regressor,:),'linewidth',2);
        title({'Binned explained var',norm_type});
        xlabel('Normalized explained var');
        ylabel('Depth (\mum)')
        
    end
end

linkaxes(get(h,'Children'),'y');


% Concatenate all units (normalized distance)

% Get estimation of end of striatum for each recording
ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
ephys_align_fn = ['ephys_depth_align.mat'];
load([ephys_align_path filesep ephys_align_fn]);
str_depth_cat = vertcat(ephys_depth_align(1:6).str_depth);


animaldays = reshape(arrayfun(@(x) ~isempty(unit_kernel_all(x).template_depths), ...
    1:numel(unit_kernel_all)),size(unit_kernel_all));

str_units = cellfun(@(unit_depths,str_depths) ...
    unit_depths >= str_depths(1) & unit_depths <= str_depths(2), ...
    AP_index_ans(reshape({unit_kernel_all.template_depths}, ...
    size(unit_kernel_all))',animaldays'), ...
    mat2cell(str_depth_cat,ones(size(str_depth_cat,1),1),2),'uni',false);

use_partial = 2;

template_depths_cat = cellfun(@(unit_depths,str_depths,str_units) ...
    unit_depths(str_units) - str_depths(2), ...
    AP_index_ans(reshape({unit_kernel_all.template_depths}, ...
    size(unit_kernel_all))',animaldays'), ...
    mat2cell(str_depth_cat,ones(size(str_depth_cat,1),1),2), ...
    str_units,'uni',false);
unit_expl_var_total_cat = cellfun(@(x,str_units) x(str_units), ...
    AP_index_ans(reshape({unit_kernel_all.unit_expl_var_total}, ...
    size(unit_kernel_all))',animaldays'),str_units,'uni',false);
unit_expl_var_partial_cat = cellfun(@(x,str_units) x(str_units,:,:), ...
    AP_index_ans(reshape({unit_kernel_all.unit_expl_var_partial}, ...
    size(unit_kernel_all))',animaldays'),str_units,'uni',false);
spike_rate_cat = cellfun(@(x,str_units) x(str_units), ...
    AP_index_ans(reshape({unit_kernel_all.spike_rate}, ...
    size(unit_kernel_all))',animaldays'),str_units,'uni',false);

[~,max_regressor_idx_cat] = cellfun(@(x) max(x(:,:,use_partial),[],2),unit_expl_var_partial_cat,'uni',false);

h = figure;

subplot(1,length(regressor_labels)+3,1,'YDir','reverse'); hold on;
scatter(log10(cell2mat(spike_rate_cat)), ...
    cell2mat(template_depths_cat),20, ...
    regressor_cols(cell2mat(max_regressor_idx_cat),:),'filled')
xlabel('log10(spike rate)');
ylabel('Depth from str end (\mum)');
title({'All units','Best regressor'})
line(xlim,[0,0],'color','k','linewidth',2);

subplot(1,length(regressor_labels)+3,2,'YDir','reverse'); hold on;
curr_expl_var = cell2mat(unit_expl_var_total_cat);
curr_expl_var_dotsize = 100*rescale(curr_expl_var.*(curr_expl_var > 0),0,1) + 1;
curr_expl_var_dotsize(isnan(curr_expl_var) | curr_expl_var < -1 | curr_expl_var > 1) = NaN;
scatter(log10(cell2mat(spike_rate_cat)), ...
    cell2mat(template_depths_cat),curr_expl_var_dotsize, ...
    regressor_cols(cell2mat(max_regressor_idx_cat),:),'filled')
xlabel('log10(spike rate)');
ylabel('Depth from str end (\mum)');
title({'All regressors'})
line(xlim,[0,0],'color','k','linewidth',2);

% (plot units size-scaled by explained variance)
for curr_regressor = 1:length(regressor_labels)
    for curr_norm = 1:2
        
        p1 = subplot(2,length(regressor_labels)+3,curr_regressor+2+((length(regressor_labels)+3)*(curr_norm-1)),'YDir','reverse');
        hold on;
        curr_expl_var = cellfun(@(x) x(:,curr_regressor,use_partial),unit_expl_var_partial_cat,'uni',false);
        
        switch curr_norm
            case 1
                % Normalize within experiment
                curr_expl_var_norm = cell2mat(cellfun(@(x) x./nanmax(x(x < 1)),curr_expl_var,'uni',false));
                curr_expl_var_dotsize = 100*rescale(curr_expl_var_norm.*(curr_expl_var_norm > 0),0,1) + 1;
                norm_type = 'experiment normalized';
            case 2
                % Normalize across experiments
                curr_expl_var_cat = cell2mat(curr_expl_var);
                curr_expl_var_norm = curr_expl_var_cat./ ...
                    max(curr_expl_var_cat(curr_expl_var_cat < 1));
                curr_expl_var_dotsize = 100*rescale(curr_expl_var_norm.*(curr_expl_var_norm > 0),0,1) + 1;
                norm_type = 'concat normalized';
        end
        
        scatter(log10(cell2mat(spike_rate_cat)), ...
            cell2mat(template_depths_cat), ...
            curr_expl_var_dotsize, ...
            regressor_cols(curr_regressor,:),'filled');
        title({regressor_labels{curr_regressor},norm_type});
        xlabel('log10(spike rate)');
        ylabel('Depth from str end (\mum)');
        line(xlim,[0,0],'color','k','linewidth',2);
        
        p2 = subplot(2,length(regressor_labels)+3,curr_norm*(length(regressor_labels)+3),'YDir','reverse');
        hold on;
        depth_bins = linspace(min(cell2mat(template_depths_cat)), ...
            max(cell2mat(template_depths_cat)),round(range(cell2mat(template_depths_cat))/400));
        depth_bin_centers = depth_bins(1:end-1) + diff(depth_bins)./2;
        curr_depth_groups = discretize(cell2mat(template_depths_cat),depth_bins);
        
        rate_cutoff = -0.5;
        use_units = log10(cell2mat(spike_rate_cat)) > rate_cutoff;
        line(p1,repmat(rate_cutoff,1,2),ylim(p1),'color','k');
        
        curr_expl_var_norm_depth = accumarray(curr_depth_groups(use_units), ...
            curr_expl_var_norm(use_units),[length(depth_bin_centers),1],@nanmean);
        plot(rescale(curr_expl_var_norm_depth,0,1),depth_bin_centers,'color',regressor_cols(curr_regressor,:),'linewidth',2);
        title({'Binned explained var',norm_type});
        xlabel('Normalized explained var');
        ylabel('Depth from str end (\mum)')
        line(xlim,[0,0],'color','k','linewidth',2);
        
    end
end

% (draw striatum starts)
for curr_plot = 1:2
    subplot(2,length(regressor_labels)+3,curr_plot*(length(regressor_labels)+3))
    str_start_line = nan(size(str_depth_cat,1),1);
    for i = 1:size(str_depth_cat,1)
        str_start_line(i) = line(xlim,repmat(-diff(str_depth_cat(i,:)),1,2),'color',[0.8,0.8,0.8]);
    end
    set(gca,'Children',circshift(get(gca,'Children'),-size(str_depth_cat,1)));
end

linkaxes(get(h,'Children'),'y');


%% Plot all unit kernels (to show reward response is different in top and bottom?)

unit_kernel_dir = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
% unit_kernel_fn = 'unit_kernel_all.mat';
unit_kernel_fn = 'unit_kernel_all_triaged.mat';
load([unit_kernel_dir filesep unit_kernel_fn]);

animaldays = reshape(arrayfun(@(x) ~isempty(unit_kernel_all(x).template_depths), ...
    1:numel(unit_kernel_all)),size(unit_kernel_all));

regressor_labels = {'Stim','Move','Go cue','Outcome'};
regressor_cols = [01,0,0;1,0,1;0.8,0.7,0.2;0,0,1];

% Get estimation of end of striatum for each recording
ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
ephys_align_fn = ['ephys_depth_align.mat'];
load([ephys_align_path filesep ephys_align_fn]);
str_depth_cat = vertcat(ephys_depth_align(1:6).str_depth);

% Get units in striatum
str_units = cellfun(@(unit_depths,str_depths) ...
    unit_depths >= str_depths(1) & unit_depths <= str_depths(2), ...
    AP_index_ans(reshape({unit_kernel_all.template_depths}, ...
    size(unit_kernel_all))',animaldays'), ...
    mat2cell(str_depth_cat,ones(size(str_depth_cat,1),1),2),'uni',false);

% Concatenate unit kernel data
use_partial = 2;

template_depths_cat = cellfun(@(unit_depths,str_depths,str_units) ...
    unit_depths(str_units) - str_depths(2), ...
    AP_index_ans(reshape({unit_kernel_all.template_depths}, ...
    size(unit_kernel_all))',animaldays'), ...
    mat2cell(str_depth_cat,ones(size(str_depth_cat,1),1),2), ...
    str_units,'uni',false);

unit_expl_var_total_cat = cellfun(@(x,str_units) x(str_units), ...
    AP_index_ans(reshape({unit_kernel_all.unit_expl_var_total}, ...
    size(unit_kernel_all))',animaldays'),str_units,'uni',false);

unit_expl_var_partial_cat = cellfun(@(x,str_units) x(str_units,:,:), ...
    AP_index_ans(reshape({unit_kernel_all.unit_expl_var_partial}, ...
    size(unit_kernel_all))',animaldays'),str_units,'uni',false);

spike_rate_cat = cellfun(@(x,str_units) x(str_units), ...
    AP_index_ans(reshape({unit_kernel_all.spike_rate}, ...
    size(unit_kernel_all))',animaldays'),str_units,'uni',false);

unit_taskpred_kernel_allcat = cellfun(@(x,str_units) ...
    cellfun(@(x) x(:,:,str_units),x,'uni',false), ...
    AP_index_ans(reshape({unit_kernel_all.unit_taskpred_k}, ...
    size(unit_kernel_all))',animaldays'),str_units,'uni',false);

unit_taskpred_kernel_cat = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, unit_taskpred_kernel_allcat,'uni',false),[2,3,1])), ...
    1:length(regressor_labels),'uni',false);

[~,max_regressor_idx_cat] = cellfun(@(x) max(x(:,:,use_partial),[],2),unit_expl_var_partial_cat,'uni',false);


% Just blunt for now:
all_depths = cell2mat(template_depths_cat);
all_log_spike_rate = log10(cell2mat(spike_rate_cat));

spike_cutoff_units = all_log_spike_rate > -0.5;

curr_regressor = 4;
a = unit_taskpred_kernel_cat{curr_regressor}./std(unit_taskpred_kernel_cat{curr_regressor},[],2);

kernel_mean = nanmean(a,3);
figure;hold on
plot(kernel_mean')

% kernel_superficial = nanmean(a(:,:,spike_cutoff_units & all_depths < -1500),3);
% kernel_deep = nanmean(a(:,:,spike_cutoff_units & all_depths >= -1500),3);
%
% figure;hold on
% plot(kernel_superficial')
% plot(kernel_deep')




%% Cortex -> striatum unit regression

upsample_factor = 1;
sample_rate = (1/median(diff(frame_t)))*upsample_factor;

% Skip the first/last n seconds to do this
skip_seconds = 60;

time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

use_spikes = spike_times_timeline > time_bins(1) & spike_times_timeline < time_bins(end);
template_spike_n = accumarray(spike_templates(use_spikes),1,[size(templates,1),1]);
use_templates = find(template_spike_n > 0);

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

use_svs = 1:50;
kernel_t = [-0.2,0.2];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
lambda = 10; % (comment out to use lambda from above)
zs = [false,false];
cvfold = 5;

% Regress cortex to units (using deconvolved fluorescence)
[k,predicted_spikes,explained_var] = ...
    AP_regresskernel(fVdf_deconv_resample(use_svs,:), ...
    binned_spikes_std, ...
    kernel_frames,lambda,zs,cvfold);

% Reshape kernel and convert to pixel space, get single-frame map
k_px = arrayfun(@(x) svdFrameReconstruct(Udf(:,:,use_svs),k(:,:,x)),1:size(k,3),'uni',false);
k_px = cat(4,k_px{:});

% (kernel t = 0:1 mean)
k_px_t01 = squeeze(nanmean(k_px(:,:,ismember(kernel_frames,[0,1]),:),3));

ctx_str_unit_map = nan(size(Udf,1),size(Udf,2),size(templates,3));
ctx_str_unit_map(:,:,use_templates) = k_px_t01;

AP_imscroll(ctx_str_unit_map)
axis image;
caxis([-0.005,0.005])
colormap(brewermap([],'*RdBu'));


%% (after unit kernel and getting ctx>str maps)

% Load unit kernel
unit_kernel_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\unit_kernel_all.mat';
load(unit_kernel_fn);

regressor_labels = {'Stim','Move','Go cue','Outcome'};
regressor_cols = [01,0,0;1,0,1;0.8,0.7,0.2;0,0,1];

% Get striatum depths in recording
ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
ephys_align_fn = ['ephys_depth_align.mat'];
load([ephys_align_path filesep ephys_align_fn]);
curr_animal_idx = strcmp(animal,{ephys_depth_align.animal});
curr_day_idx = strcmp(day,ephys_depth_align(curr_animal_idx).day);

curr_str_depths = ephys_depth_align(curr_animal_idx).str_depth(curr_day_idx,:);
curr_str_units = template_depths > curr_str_depths(1) & template_depths < curr_str_depths(2);

% Correlation between map and depth vs. explained variance similarity
depth_dist = abs(template_depths - template_depths');

use_partial = 2;
unit_expl_var = unit_kernel_all(curr_animal_idx,curr_day_idx).unit_expl_var_partial(:,:,use_partial);
unit_expl_var_norm = rescale(unit_expl_var,0,1,'InputMin',0,'InputMax',1);
unit_expl_var_norm(isnan(unit_expl_var) | unit_expl_var < -1 | unit_expl_var >= 1) = NaN;

regressor_dist = squareform(pdist(unit_expl_var_norm,'correlation'));
map_corr = corrcoef(reshape(ctx_str_unit_map,[],size(templates,1)));

depth_dist_str = AP_itril(depth_dist(curr_str_units,curr_str_units),-1);
regressor_dist_str = AP_itril(regressor_dist(curr_str_units,curr_str_units),-1);
map_corr_str = AP_itril(map_corr(curr_str_units,curr_str_units),-1);

% Plot by min-max bins
dist_bins = linspace(0,1,20);
dist_bin_centers = dist_bins(1:end-1) + diff(dist_bins)/2;

depth_dist_bin = discretize(mat2gray(depth_dist_str,prctile(depth_dist_str,[0,100])),dist_bins);
map_corr_depth = accumarray(depth_dist_bin,map_corr_str,[length(dist_bin_centers),1],@nanmean);

regressor_dist_bin = discretize(mat2gray(regressor_dist_str,prctile(regressor_dist_str,[0,100])),dist_bins);
map_corr_regressor = accumarray(regressor_dist_bin,map_corr_str,[length(dist_bin_centers),1],@nanmean);

figure; hold on
plot(dist_bin_centers,map_corr_depth,'k','linewidth',2);
plot(dist_bin_centers,map_corr_regressor,'r','linewidth',2);
legend({'Depth','Regressor explained variance'});
ylabel('Cortex map correlation');
xlabel('Unit distance');


%% Additive fit between MUA and cortex (quick and dirty)

% Choose split for data
trials_allcat = size(mua_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_animal;

% (make outcome-aligned data, here for now)
mua_allcat_outcome = mua_allcat;
fluor_roi_deconv_outcome = fluor_roi_deconv;

t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
for i = 1:size(mua_allcat,1)
    mua_allcat_outcome(i,:,:) = circshift(mua_allcat_outcome(i,:,:),-outcome_idx(i)+leeway_samples,2);
    fluor_roi_deconv_outcome(i,:,:) = circshift(fluor_roi_deconv_outcome(i,:,:),-outcome_idx(i)+leeway_samples,2);
end

% Split MUA and fluorescence
mua_allcat_exp = mat2cell(mua_allcat,use_split,length(t),n_depths);
fluor_roi_deconv_exp = mat2cell(fluor_roi_deconv,use_split,length(t),n_rois);

mua_allcat_move_exp = mat2cell(mua_allcat_move,use_split,length(t),n_depths);
fluor_roi_deconv_move_exp = mat2cell(fluor_roi_deconv_move,use_split,length(t),n_rois);

mua_allcat_outcome_exp = mat2cell(mua_allcat_outcome,use_split,length(t),n_depths);
fluor_roi_deconv_outcome_exp = mat2cell(fluor_roi_deconv_outcome,use_split,length(t),n_rois);

% Get average activity within window
trial_groups = {'Stim','Move onset'};
t_baseline = [-0.1,-0.05];
t_event = {[0.05,0.15],[-0.05,0.05],[0,0.1]};


% NOT LOOPING FOR NOW - JUST HARDCODING PLOTS

t_col = colormap_BlueWhiteRed(length(t));
t_col(length(t)+1,:) = [];
t_col = permute(reshape(t_col',3,[],2),[2,1,3]);
t_col(:,:,1) = flipud(t_col(:,:,1));

% STIM
curr_event = 1;
figure; subplot(1,2,1); hold on;
trial_type_fit = nan(length(t),2);
for curr_t = 1:length(t)
    
    curr_mua_baseline = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        mua_allcat_exp,'uni',false);
    
    curr_mua_event = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        mua_allcat_exp,'uni',false);
    
    curr_fluor_baseline = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        fluor_roi_deconv_exp,'uni',false);
    
    curr_fluor_event = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        fluor_roi_deconv_exp,'uni',false);
    
    
    curr_depth = 1;
    curr_roi = 3;
    
    % (concatenate all data)
    curr_mua_baseline_cat = cell2mat(cellfun(@(x) x(:,:,curr_depth),curr_mua_baseline,'uni',false));
    curr_mua_event_cat = cell2mat(cellfun(@(x) x(:,:,curr_depth),curr_mua_event,'uni',false));
    
    curr_mua_baseline_cat = curr_mua_baseline_cat(trial_stim_allcat < 0);
    curr_mua_event_cat = curr_mua_event_cat(trial_stim_allcat > 0);
    
    curr_fluor_baseline_cat = cell2mat(cellfun(@(x) x(:,:,curr_roi),curr_fluor_baseline,'uni',false));
    curr_fluor_event_cat = cell2mat(cellfun(@(x) x(:,:,curr_roi),curr_fluor_event,'uni',false));
    
    curr_fluor_baseline_cat = curr_fluor_baseline_cat(trial_stim_allcat < 0);
    curr_fluor_event_cat = curr_fluor_event_cat(trial_stim_allcat > 0);
    
    % (bin data by fluorescence)
    fluor_bin_range = prctile([curr_fluor_baseline_cat;curr_fluor_event_cat],[5,95]);
    fluor_bin_edges = linspace(fluor_bin_range(1),fluor_bin_range(2),10);
    fluor_bin_centers = fluor_bin_edges(1:end-1) + diff(fluor_bin_edges)./2;
    
    fluor_bins_baseline = discretize(curr_fluor_baseline_cat,fluor_bin_edges);
    fluor_bins_event = discretize(curr_fluor_event_cat,fluor_bin_edges);
    
    curr_mua_baseline_binned = ...
        accumarray(fluor_bins_baseline(~isnan(fluor_bins_baseline)), ...
        curr_mua_baseline_cat(~isnan(fluor_bins_baseline)), ...
        [length(fluor_bin_centers),1],@nanmean,NaN);
    
    curr_mua_event_binned = ...
        accumarray(fluor_bins_event(~isnan(fluor_bins_event)), ...
        curr_mua_event_cat(~isnan(fluor_bins_event)), ...
        [length(fluor_bin_centers),1],@nanmean,NaN);
    
    %     % (to fit multiplicative and additive)
    %     trial_type_fit(curr_t,:) = ...
    %         [curr_mua_baseline_binned,ones(size(curr_mua_baseline_binned))]\curr_mua_event_binned;
    
    % (to fit only additive);
    trial_type_fit(curr_t,:) = ...
        ones(size(curr_mua_baseline_binned))\(curr_mua_event_binned - curr_mua_baseline_binned);
    
    % Plot
    plot(fluor_bin_centers,curr_mua_baseline_binned,'color',t_col(curr_t,:,1),'linewidth',2);
    plot(fluor_bin_centers,curr_mua_event_binned,'color',t_col(curr_t,:,2),'linewidth',2);
    xlabel([wf_roi(curr_roi).area ' fluorescence']);
    ylabel(['Str ' num2str(curr_depth) ' MUA']);
end

subplot(1,2,2); hold on
plot(t,trial_type_fit,'linewidth',2);
xlabel('Time from event');
ylabel('Additive offset');


% MOVE
curr_event = 2;
figure; subplot(1,2,1); hold on;
trial_type_fit = nan(length(t),2);
for curr_t = 1:length(t)
    
    curr_mua_baseline = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        mua_allcat_move_exp,'uni',false);
    
    curr_mua_event = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        mua_allcat_move_exp,'uni',false);
    
    curr_fluor_baseline = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        fluor_roi_deconv_move_exp,'uni',false);
    
    curr_fluor_event = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        fluor_roi_deconv_move_exp,'uni',false);
    
    
    curr_depth = 2;
    curr_roi = 7;
    
    % (concatenate all data)
    curr_mua_baseline_cat = cell2mat(cellfun(@(x) x(:,:,curr_depth),curr_mua_baseline,'uni',false));
    curr_mua_event_cat = cell2mat(cellfun(@(x) x(:,:,curr_depth),curr_mua_event,'uni',false));
    
    curr_mua_baseline_cat = curr_mua_baseline_cat(trial_choice_allcat == 1);
    curr_mua_event_cat = curr_mua_event_cat(trial_choice_allcat == -1);
    
    curr_fluor_baseline_cat = cell2mat(cellfun(@(x) x(:,:,curr_roi),curr_fluor_baseline,'uni',false));
    curr_fluor_event_cat = cell2mat(cellfun(@(x) x(:,:,curr_roi),curr_fluor_event,'uni',false));
    
    curr_fluor_baseline_cat = curr_fluor_baseline_cat(trial_choice_allcat == 1);
    curr_fluor_event_cat = curr_fluor_event_cat(trial_choice_allcat == -1);
    
    % (bin data by fluorescence)
    fluor_bin_range = prctile([curr_fluor_baseline_cat;curr_fluor_event_cat],[5,95]);
    fluor_bin_edges = linspace(fluor_bin_range(1),fluor_bin_range(2),10);
    fluor_bin_centers = fluor_bin_edges(1:end-1) + diff(fluor_bin_edges)./2;
    
    fluor_bins_baseline = discretize(curr_fluor_baseline_cat,fluor_bin_edges);
    fluor_bins_event = discretize(curr_fluor_event_cat,fluor_bin_edges);
    
    curr_mua_baseline_binned = ...
        accumarray(fluor_bins_baseline(~isnan(fluor_bins_baseline)), ...
        curr_mua_baseline_cat(~isnan(fluor_bins_baseline)), ...
        [length(fluor_bin_centers),1],@nanmean,NaN);
    
    curr_mua_event_binned = ...
        accumarray(fluor_bins_event(~isnan(fluor_bins_event)), ...
        curr_mua_event_cat(~isnan(fluor_bins_event)), ...
        [length(fluor_bin_centers),1],@nanmean,NaN);
    
    %     % (to fit multiplicative and additive)
    %     trial_type_fit(curr_t,:) = ...
    %         [curr_mua_baseline_binned,ones(size(curr_mua_baseline_binned))]\curr_mua_event_binned;
    
    % (to fit only additive);
    trial_type_fit(curr_t,:) = ...
        ones(size(curr_mua_baseline_binned))\(curr_mua_event_binned - curr_mua_baseline_binned);
    
    % Plot
    plot(fluor_bin_centers,curr_mua_baseline_binned,'color',t_col(curr_t,:,1),'linewidth',2);
    plot(fluor_bin_centers,curr_mua_event_binned,'color',t_col(curr_t,:,2),'linewidth',2);
    xlabel([wf_roi(curr_roi).area ' fluorescence']);
    ylabel(['Str ' num2str(curr_depth) ' MUA']);
end

subplot(1,2,2); hold on
plot(t,trial_type_fit,'linewidth',2);
xlabel('Time from event');
ylabel('Additive offset');


% OUTCOME
curr_event = 3;
figure; subplot(1,2,1); hold on;
trial_type_fit = nan(length(t),2);
for curr_t = 1:length(t)
    
    curr_mua_baseline = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        mua_allcat_outcome_exp,'uni',false);
    
    curr_mua_event = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        mua_allcat_outcome_exp,'uni',false);
    
    curr_fluor_baseline = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        fluor_roi_deconv_outcome_exp,'uni',false);
    
    curr_fluor_event = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        fluor_roi_deconv_outcome_exp,'uni',false);
    
    
    curr_depth = 4;
    curr_roi = 10;
    
    % (concatenate all data)
    curr_mua_baseline_cat = cell2mat(cellfun(@(x) x(:,:,curr_depth),curr_mua_baseline,'uni',false));
    curr_mua_event_cat = cell2mat(cellfun(@(x) x(:,:,curr_depth),curr_mua_event,'uni',false));
    
    curr_mua_baseline_cat = curr_mua_baseline_cat(trial_outcome_allcat == -1);
    curr_mua_event_cat = curr_mua_event_cat(trial_outcome_allcat == 1);
    
    curr_fluor_baseline_cat = cell2mat(cellfun(@(x) x(:,:,curr_roi),curr_fluor_baseline,'uni',false));
    curr_fluor_event_cat = cell2mat(cellfun(@(x) x(:,:,curr_roi),curr_fluor_event,'uni',false));
    
    curr_fluor_baseline_cat = curr_fluor_baseline_cat(trial_outcome_allcat == -1);
    curr_fluor_event_cat = curr_fluor_event_cat(trial_outcome_allcat == 1);
    
    % (bin data by fluorescence)
    fluor_bin_range = prctile([curr_fluor_baseline_cat;curr_fluor_event_cat],[5,95]);
    fluor_bin_edges = linspace(fluor_bin_range(1),fluor_bin_range(2),10);
    fluor_bin_centers = fluor_bin_edges(1:end-1) + diff(fluor_bin_edges)./2;
    
    fluor_bins_baseline = discretize(curr_fluor_baseline_cat,fluor_bin_edges);
    fluor_bins_event = discretize(curr_fluor_event_cat,fluor_bin_edges);
    
    curr_mua_baseline_binned = ...
        accumarray(fluor_bins_baseline(~isnan(fluor_bins_baseline)), ...
        curr_mua_baseline_cat(~isnan(fluor_bins_baseline)), ...
        [length(fluor_bin_centers),1],@nanmean,NaN);
    
    curr_mua_event_binned = ...
        accumarray(fluor_bins_event(~isnan(fluor_bins_event)), ...
        curr_mua_event_cat(~isnan(fluor_bins_event)), ...
        [length(fluor_bin_centers),1],@nanmean,NaN);
    
    %     % (to fit multiplicative and additive)
    %     trial_type_fit(curr_t,:) = ...
    %         [curr_mua_baseline_binned,ones(size(curr_mua_baseline_binned))]\curr_mua_event_binned;
    
    % (to fit only additive);
    trial_type_fit(curr_t,:) = ...
        ones(size(curr_mua_baseline_binned))\(curr_mua_event_binned - curr_mua_baseline_binned);
    
    % Plot
    plot(fluor_bin_centers,curr_mua_baseline_binned,'color',t_col(curr_t,:,1),'linewidth',2);
    plot(fluor_bin_centers,curr_mua_event_binned,'color',t_col(curr_t,:,2),'linewidth',2);
    xlabel([wf_roi(curr_roi).area ' fluorescence']);
    ylabel(['Str ' num2str(curr_depth) ' MUA']);
end

subplot(1,2,2); hold on
plot(t,trial_type_fit,'linewidth',2);
xlabel('Time from event');
ylabel('Additive offset');



%% (Additive model as above, but using kernel-based ROIs)

% Get fluorescence in kernel-based ROIs
kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
load(kernel_roi_fn);
fluor_kernel_roi_bw = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],kernel_roi.bw), ...
    size(kernel_roi.bw,3),[],size(fluor_allcat_deconv,1)),[3,2,1]);

% (align to move)
fluor_kernel_roi_bw_move = fluor_kernel_roi_bw;
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
for i = 1:size(mua_allcat,1)
    fluor_kernel_roi_bw_move(i,:,:) = circshift(fluor_kernel_roi_bw_move(i,:,:),-move_idx(i)+leeway_samples,2);
end



% Choose split for data
trials_allcat = size(mua_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_animal;

% (make outcome-aligned data, here for now)
mua_allcat_outcome = mua_allcat;
fluor_kernel_roi_bw_outcome = fluor_kernel_roi_bw;

t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
for i = 1:size(mua_allcat,1)
    mua_allcat_outcome(i,:,:) = circshift(mua_allcat_outcome(i,:,:),-outcome_idx(i)+leeway_samples,2);
    fluor_kernel_roi_bw_outcome(i,:,:) = circshift(fluor_kernel_roi_bw_outcome(i,:,:),-outcome_idx(i)+leeway_samples,2);
end

% Split MUA and fluorescence
mua_allcat_exp = mat2cell(mua_allcat,use_split,length(t),n_depths);
fluor_kernel_roi_bw_exp = mat2cell(fluor_kernel_roi_bw,use_split,length(t),n_depths);

mua_allcat_move_exp = mat2cell(mua_allcat_move,use_split,length(t),n_depths);
fluor_kernel_roi_bw_move_exp = mat2cell(fluor_kernel_roi_bw_move,use_split,length(t),n_depths);

mua_allcat_outcome_exp = mat2cell(mua_allcat_outcome,use_split,length(t),n_depths);
fluor_kernel_roi_bw_outcome_exp = mat2cell(fluor_kernel_roi_bw_outcome,use_split,length(t),n_depths);

% Get average activity within window
trial_groups = {'Stim','Move onset'};
t_baseline = [-0.1,-0.05];
t_event = {[0.05,0.15],[-0.05,0.05],[0,0.1]};


% NOT LOOPING FOR NOW - JUST HARDCODING PLOTS

t_col = colormap_BlueWhiteRed(length(t));
t_col(length(t)+1,:) = [];
t_col = permute(reshape(t_col',3,[],2),[2,1,3]);
t_col(:,:,1) = flipud(t_col(:,:,1));

% STIM
curr_event = 1;
figure; subplot(1,2,1); hold on;
trial_type_fit = nan(length(t),2);
for curr_t = 1:length(t)
    
    curr_mua_baseline = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        mua_allcat_exp,'uni',false);
    
    curr_mua_event = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        mua_allcat_exp,'uni',false);
    
    curr_fluor_baseline = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        fluor_kernel_roi_bw_exp,'uni',false);
    
    curr_fluor_event = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        fluor_kernel_roi_bw_exp,'uni',false);
    
    
    curr_depth = 1;
    curr_roi = 1;
    
    % (concatenate all data)
    curr_mua_baseline_cat = cell2mat(cellfun(@(x) x(:,:,curr_depth),curr_mua_baseline,'uni',false));
    curr_mua_event_cat = cell2mat(cellfun(@(x) x(:,:,curr_depth),curr_mua_event,'uni',false));
    
    curr_mua_baseline_cat = curr_mua_baseline_cat(trial_contrast_allcat > 0 & trial_side_allcat == -1);
    curr_mua_event_cat = curr_mua_event_cat(trial_contrast_allcat > 0 & trial_side_allcat == 1);
    
    curr_fluor_baseline_cat = cell2mat(cellfun(@(x) x(:,:,curr_roi),curr_fluor_baseline,'uni',false));
    curr_fluor_event_cat = cell2mat(cellfun(@(x) x(:,:,curr_roi),curr_fluor_event,'uni',false));
    
    curr_fluor_baseline_cat = curr_fluor_baseline_cat(trial_contrast_allcat > 0 & trial_side_allcat == -1);
    curr_fluor_event_cat = curr_fluor_event_cat(trial_contrast_allcat > 0 & trial_side_allcat == 1);
    
    % (bin data by fluorescence)
    fluor_bin_range = prctile([curr_fluor_baseline_cat;curr_fluor_event_cat],[5,95]);
    fluor_bin_edges = linspace(fluor_bin_range(1),fluor_bin_range(2),10);
    fluor_bin_centers = fluor_bin_edges(1:end-1) + diff(fluor_bin_edges)./2;
    
    fluor_bins_baseline = discretize(curr_fluor_baseline_cat,fluor_bin_edges);
    fluor_bins_event = discretize(curr_fluor_event_cat,fluor_bin_edges);
    
    curr_mua_baseline_binned = ...
        accumarray(fluor_bins_baseline(~isnan(fluor_bins_baseline)), ...
        curr_mua_baseline_cat(~isnan(fluor_bins_baseline)), ...
        [length(fluor_bin_centers),1],@nanmean,NaN);
    
    curr_mua_event_binned = ...
        accumarray(fluor_bins_event(~isnan(fluor_bins_event)), ...
        curr_mua_event_cat(~isnan(fluor_bins_event)), ...
        [length(fluor_bin_centers),1],@nanmean,NaN);
    
    %     % (to fit multiplicative and additive)
    %     trial_type_fit(curr_t,:) = ...
    %         [curr_mua_baseline_binned,ones(size(curr_mua_baseline_binned))]\curr_mua_event_binned;
    
    % (to fit only additive);
    trial_type_fit(curr_t,:) = ...
        ones(size(curr_mua_baseline_binned))\(curr_mua_event_binned - curr_mua_baseline_binned);
    
    % Plot
    plot(fluor_bin_centers,curr_mua_baseline_binned,'color',t_col(curr_t,:,1),'linewidth',2);
    plot(fluor_bin_centers,curr_mua_event_binned,'color',t_col(curr_t,:,2),'linewidth',2);
    ylabel(['Str ' num2str(curr_depth) ' MUA']);
end

subplot(1,2,2); hold on
plot(t,trial_type_fit,'linewidth',2);
xlabel('Time from event');
ylabel('Additive offset');


% MOVE
curr_event = 2;
figure; subplot(1,2,1); hold on;
trial_type_fit = nan(length(t),2);
for curr_t = 1:length(t)
    
    curr_mua_baseline = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        mua_allcat_move_exp,'uni',false);
    
    curr_mua_event = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        mua_allcat_move_exp,'uni',false);
    
    curr_fluor_baseline = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        fluor_kernel_roi_bw_move_exp,'uni',false);
    
    curr_fluor_event = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        fluor_kernel_roi_bw_move_exp,'uni',false);
    
    
    curr_depth = 2;
    curr_roi = 2;
    
    % (concatenate all data)
    curr_mua_baseline_cat = cell2mat(cellfun(@(x) x(:,:,curr_depth),curr_mua_baseline,'uni',false));
    curr_mua_event_cat = cell2mat(cellfun(@(x) x(:,:,curr_depth),curr_mua_event,'uni',false));
    
    curr_mua_baseline_cat = curr_mua_baseline_cat(trial_choice_allcat == 1);
    curr_mua_event_cat = curr_mua_event_cat(trial_choice_allcat == -1);
    
    curr_fluor_baseline_cat = cell2mat(cellfun(@(x) x(:,:,curr_roi),curr_fluor_baseline,'uni',false));
    curr_fluor_event_cat = cell2mat(cellfun(@(x) x(:,:,curr_roi),curr_fluor_event,'uni',false));
    
    curr_fluor_baseline_cat = curr_fluor_baseline_cat(trial_choice_allcat == 1);
    curr_fluor_event_cat = curr_fluor_event_cat(trial_choice_allcat == -1);
    
    % (bin data by fluorescence)
    fluor_bin_range = prctile([curr_fluor_baseline_cat;curr_fluor_event_cat],[5,95]);
    fluor_bin_edges = linspace(fluor_bin_range(1),fluor_bin_range(2),10);
    fluor_bin_centers = fluor_bin_edges(1:end-1) + diff(fluor_bin_edges)./2;
    
    fluor_bins_baseline = discretize(curr_fluor_baseline_cat,fluor_bin_edges);
    fluor_bins_event = discretize(curr_fluor_event_cat,fluor_bin_edges);
    
    curr_mua_baseline_binned = ...
        accumarray(fluor_bins_baseline(~isnan(fluor_bins_baseline)), ...
        curr_mua_baseline_cat(~isnan(fluor_bins_baseline)), ...
        [length(fluor_bin_centers),1],@nanmean,NaN);
    
    curr_mua_event_binned = ...
        accumarray(fluor_bins_event(~isnan(fluor_bins_event)), ...
        curr_mua_event_cat(~isnan(fluor_bins_event)), ...
        [length(fluor_bin_centers),1],@nanmean,NaN);
    
    %     % (to fit multiplicative and additive)
    %     trial_type_fit(curr_t,:) = ...
    %         [curr_mua_baseline_binned,ones(size(curr_mua_baseline_binned))]\curr_mua_event_binned;
    
    % (to fit only additive);
    trial_type_fit(curr_t,:) = ...
        ones(size(curr_mua_baseline_binned))\(curr_mua_event_binned - curr_mua_baseline_binned);
    
    % Plot
    plot(fluor_bin_centers,curr_mua_baseline_binned,'color',t_col(curr_t,:,1),'linewidth',2);
    plot(fluor_bin_centers,curr_mua_event_binned,'color',t_col(curr_t,:,2),'linewidth',2);
    ylabel(['Str ' num2str(curr_depth) ' MUA']);
end

subplot(1,2,2); hold on
plot(t,trial_type_fit,'linewidth',2);
xlabel('Time from event');
ylabel('Additive offset');


% OUTCOME
curr_event = 3;
figure; subplot(1,2,1); hold on;
trial_type_fit = nan(length(t),2);
for curr_t = 1:length(t)
    
    curr_mua_baseline = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        mua_allcat_outcome_exp,'uni',false);
    
    curr_mua_event = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        mua_allcat_outcome_exp,'uni',false);
    
    curr_fluor_baseline = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        fluor_kernel_roi_bw_outcome_exp,'uni',false);
    
    curr_fluor_event = cellfun(@(x) nanmean(x(:,curr_t,:),2), ...
        fluor_kernel_roi_bw_outcome_exp,'uni',false);
    
    
    curr_depth = 4;
    curr_roi = 4;
    
    % (concatenate all data)
    curr_mua_baseline_cat = cell2mat(cellfun(@(x) x(:,:,curr_depth),curr_mua_baseline,'uni',false));
    curr_mua_event_cat = cell2mat(cellfun(@(x) x(:,:,curr_depth),curr_mua_event,'uni',false));
    
    curr_mua_baseline_cat = curr_mua_baseline_cat(trial_outcome_allcat == -1);
    curr_mua_event_cat = curr_mua_event_cat(trial_outcome_allcat == 1);
    
    curr_fluor_baseline_cat = cell2mat(cellfun(@(x) x(:,:,curr_roi),curr_fluor_baseline,'uni',false));
    curr_fluor_event_cat = cell2mat(cellfun(@(x) x(:,:,curr_roi),curr_fluor_event,'uni',false));
    
    curr_fluor_baseline_cat = curr_fluor_baseline_cat(trial_outcome_allcat == -1);
    curr_fluor_event_cat = curr_fluor_event_cat(trial_outcome_allcat == 1);
    
    % (bin data by fluorescence)
    fluor_bin_range = prctile([curr_fluor_baseline_cat;curr_fluor_event_cat],[5,95]);
    fluor_bin_edges = linspace(fluor_bin_range(1),fluor_bin_range(2),10);
    fluor_bin_centers = fluor_bin_edges(1:end-1) + diff(fluor_bin_edges)./2;
    
    fluor_bins_baseline = discretize(curr_fluor_baseline_cat,fluor_bin_edges);
    fluor_bins_event = discretize(curr_fluor_event_cat,fluor_bin_edges);
    
    curr_mua_baseline_binned = ...
        accumarray(fluor_bins_baseline(~isnan(fluor_bins_baseline)), ...
        curr_mua_baseline_cat(~isnan(fluor_bins_baseline)), ...
        [length(fluor_bin_centers),1],@nanmean,NaN);
    
    curr_mua_event_binned = ...
        accumarray(fluor_bins_event(~isnan(fluor_bins_event)), ...
        curr_mua_event_cat(~isnan(fluor_bins_event)), ...
        [length(fluor_bin_centers),1],@nanmean,NaN);
    
    %     % (to fit multiplicative and additive)
    %     trial_type_fit(curr_t,:) = ...
    %         [curr_mua_baseline_binned,ones(size(curr_mua_baseline_binned))]\curr_mua_event_binned;
    
    % (to fit only additive);
    trial_type_fit(curr_t,:) = ...
        ones(size(curr_mua_baseline_binned))\(curr_mua_event_binned - curr_mua_baseline_binned);
    
    % Plot
    plot(fluor_bin_centers,curr_mua_baseline_binned,'color',t_col(curr_t,:,1),'linewidth',2);
    plot(fluor_bin_centers,curr_mua_event_binned,'color',t_col(curr_t,:,2),'linewidth',2);
    ylabel(['Str ' num2str(curr_depth) ' MUA']);
end

subplot(1,2,2); hold on
plot(t,trial_type_fit,'linewidth',2);
xlabel('Time from event');
ylabel('Additive offset');

%% MUA (measured, task-pred, ctx-pred) v fluorescence by condition

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_contrast_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set windows to average activity
timeavg_labels = {'Pre-stim','Stim','Move onset','Outcome'};
timeavg_t = {[-0.2,-0.1],[0.05,0.15],[-0.05,0.05],[0.05,0.15]};
timeavg_align = {stim_align,stim_align,move_align,outcome_align};

% Set activity percentiles and bins
act_prctile = [10,90];
n_act_bins = 5;

% Set areas and conditions
plot_str_ctx = {[1,3],[2,7],[3,8],[4,10]};
trial_condition_groups = ...
    {[sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1], ...
    [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
    [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
    [trial_outcome_allcat == 1, trial_outcome_allcat == -1]};

% Loop across area pairs, plot binned predicted v measured activity
for curr_str_ctx = 1:length(plot_str_ctx)
    
    plot_str = plot_str_ctx{curr_str_ctx}(1);
    plot_ctx = plot_str_ctx{curr_str_ctx}(2);
    trial_conditions = trial_condition_groups{curr_str_ctx};
    
    %     use_reduced = 4;
    
    str_v_ctx_fig = figure('color','w');
    for curr_mua = 1:3
        
        % (set striatum activity to use)
        switch curr_mua
            case 1
                curr_str_act_allcat = single(mua_allcat);
                
                %                 curr_str_act_reduced = mua_taskpred_reduced_allcat(:,:,:,use_reduced);
                %                 curr_str_act_reduced(isnan(curr_str_act_reduced)) = 0;
                %                 curr_str_act_allcat = single(mua_allcat - curr_str_act_reduced);
                mua_label = 'Measured';
            case 2
                curr_str_act_allcat = single(mua_taskpred_allcat);
                
                %                 curr_str_act_reduced = mua_taskpred_reduced_allcat(:,:,:,use_reduced);
                %                 curr_str_act_reduced(isnan(curr_str_act_reduced)) = 0;
                %                 curr_str_act_allcat = single(mua_taskpred_allcat - curr_str_act_reduced);
                mua_label = 'Task-pred';
            case 3
                curr_str_act_allcat = single(mua_ctxpred_allcat);
                
                %                 curr_str_act_reduced = mua_ctxpred_taskpred_reduced_allcat(:,:,:,use_reduced);
                %                 curr_str_act_reduced(isnan(curr_str_act_reduced)) = 0;
                %                 curr_str_act_allcat = single(mua_ctxpred_allcat - curr_str_act_reduced);
                mua_label = 'Ctx-pred';
        end
        
        for curr_timeavg = 1:length(timeavg_labels)
            
            % (set cortex activity to use)
            curr_ctx_act_allcat = fluor_roi_deconv;
            %             curr_ctx_act_reduced = fluor_roi_taskpred_reduced(:,:,:,use_reduced);
            %             curr_ctx_act_reduced(isnan(curr_ctx_act_reduced)) = 0;
            %             curr_ctx_act_allcat = fluor_roi_deconv - curr_ctx_act_reduced;
            
            % (re-align and split activity and conditions)
            curr_ctx_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift(curr_ctx_act_allcat(trial,:,:), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_ctx_act_allcat,1)),'uni',false)), ...
                use_split,length(t),size(curr_ctx_act_allcat,3));
            
            curr_str_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift(curr_str_act_allcat(trial,:,:), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_str_act_allcat,1)),'uni',false)), ...
                use_split,length(t),size(curr_str_act_allcat,3));
            
            trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
            
            % (get average activity within window)
            curr_event_t = t > timeavg_t{curr_timeavg}(1) & t < timeavg_t{curr_timeavg}(2);
            curr_ctx_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_ctx_act,'uni',false);
            curr_str_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_str_act,'uni',false);
            
            % (bin measured data across percentile range)
            bin_range = prctile(cell2mat(cellfun(@(x) x(:,plot_ctx),curr_ctx_act_avg,'uni',false)),act_prctile);
            bin_edges = linspace(bin_range(1),bin_range(2),n_act_bins+1);
            bin_centers = bin_edges(1:end-1) + diff(bin_edges)./2;
            
            trial_bins = cellfun(@(x) discretize(x,bin_edges),curr_ctx_act_avg,'uni',false);
            total_bins = cellfun(@(x) discretize(x,[bin_edges(1),bin_edges(end)]),curr_ctx_act_avg,'uni',false);
            
            % (get trials to use: no NaNs in either time series or bins)
            nonan_trials = cellfun(@(ctx_act,str_act,trial_bins) ...
                squeeze(~any(isnan(ctx_act(:,curr_event_t,plot_ctx)),2)) & ...
                squeeze(~any(isnan(str_act(:,curr_event_t,plot_str)),2)) & ...
                ~isnan(trial_bins), ...
                curr_ctx_act,curr_str_act,trial_bins,'uni',false);
            
            % (get average binned activity for measured/predicted by condition)
            ctx_act_binmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_ctx_act_avg,trial_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            str_act_binmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_str_act_avg,trial_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            % (get the average total activity)
            ctx_act_totalmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [1,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_ctx_act_avg,total_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            str_act_totalmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [1,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_str_act_avg,total_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            % Plot binned predicted v measured
            figure(str_v_ctx_fig);
            subplot(3,length(timeavg_labels), ...
                sub2ind(fliplr([3,length(timeavg_labels)]),curr_timeavg,curr_mua));
            hold on;
            
            errorbar( ...
                squeeze(nanmean(ctx_act_binmean(:,plot_ctx,:,:),3)), ...
                squeeze(nanmean(str_act_binmean(:,plot_str,:,:),3)), ...
                squeeze(AP_sem(str_act_binmean(:,plot_str,:,:),3)), ...
                'linewidth',2,'CapSize',0);
            errorbar( ...
                squeeze(nanmean(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
                squeeze(nanmean(str_act_totalmean(:,plot_str,:,:),3)), ...
                squeeze(AP_sem(str_act_totalmean(:,plot_str,:,:),3)), ...
                squeeze(AP_sem(str_act_totalmean(:,plot_str,:,:),3)), ...
                squeeze(AP_sem(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
                squeeze(AP_sem(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
                '.','linewidth',3,'CapSize',0);
            
            xlabel(['Ctx (' wf_roi(plot_ctx).area ')']);
            ylabel([mua_label ' Str (' num2str(plot_str) ')'])
            title(timeavg_labels{curr_timeavg});
            
        end
        
    end
    % Link axes of all plots
    linkaxes(get(str_v_ctx_fig,'Children'));
end


%% MUA (meas/ctx-pred overlay) v fluorescence (across-epoch condition)

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_contrast_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set windows to average activity
timeavg_labels = {'Pre-stim','Stim','Move onset','Outcome'};
timeavg_t = {[-0.2,-0.1],[0.05,0.15],[-0.05,0.05],[0.05,0.15]};
timeavg_align = {stim_align,stim_align,move_align,outcome_align};

% Set activity percentiles and bins
act_prctile = [10,90];
n_act_bins = 5;

% Set areas and conditions
% plot_str_ctx = {[1,3],[2,7],[3,8],[4,10]};
plot_str_ctx = {[1,1],[2,2],[3,3],[4,4]};

trial_condition_groups = ...
    {[sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1], ...
    [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
    [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
    [trial_outcome_allcat == 1, trial_outcome_allcat == -1]};

% trial_condition_groups = ...
%     {[sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1], ...
%     [sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1], ...
%     [sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1], ...
%     [sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1]};

% trial_condition_groups = ...
%     {[trial_choice_allcat == -1,trial_choice_allcat == 1], ...
%     [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
%     [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
%     [trial_choice_allcat == -1,trial_choice_allcat == 1]};

% trial_condition_groups = ...
%     {[sign(trial_contrastside_allcat) == 1 & trial_choice_allcat == -1, ...
%     sign(trial_contrastside_allcat) == 1 & trial_choice_allcat == 1], ...
%     ...
%     [trial_choice_allcat == -1 & sign(trial_contrastside_allcat) == 1, ... 
%     trial_choice_allcat == -1 & sign(trial_contrastside_allcat) == 0, ...
%     trial_choice_allcat == -1 & sign(trial_contrastside_allcat) == -1], ...
%     ...
%     [trial_choice_allcat == -1 & sign(trial_contrastside_allcat) == 1, ... 
%     trial_choice_allcat == -1 & sign(trial_contrastside_allcat) == 0, ...
%     trial_choice_allcat == -1 & sign(trial_contrastside_allcat) == -1], ...
%     ...
%     [trial_choice_allcat == -1 & trial_outcome_allcat == 1, ... 
%     trial_choice_allcat == 1 & trial_outcome_allcat == 1]};

% Loop across area pairs, plot binned predicted v measured activity
str_v_ctx_fig = figure('color','w');
for curr_str_ctx = 1:length(plot_str_ctx)
    
    plot_str = plot_str_ctx{curr_str_ctx}(1);
    plot_ctx = plot_str_ctx{curr_str_ctx}(2);
    trial_conditions = trial_condition_groups{curr_str_ctx};

    for curr_mua = 1:2
        
        % (set striatum activity to use)
        switch curr_mua
            case 1
                curr_str_act_allcat = single(mua_allcat);
                mua_label = 'Measured';
            case 2
                curr_str_act_allcat = single(mua_ctxpred_allcat);
                mua_label = 'Ctx-pred';
        end
        
        for curr_timeavg = 1:length(timeavg_labels)
            
            % (set cortex activity to use)
%             curr_ctx_act_allcat = fluor_roi_deconv;
            curr_ctx_act_allcat = fluor_kernelroi_deconv;   
            
            % (re-align and split activity and conditions)
            curr_ctx_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift(curr_ctx_act_allcat(trial,:,:), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_ctx_act_allcat,1)),'uni',false)), ...
                use_split,length(t),size(curr_ctx_act_allcat,3));
            
            curr_str_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift(curr_str_act_allcat(trial,:,:), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_str_act_allcat,1)),'uni',false)), ...
                use_split,length(t),size(curr_str_act_allcat,3));
            
            trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
            
            % (get average activity within window)
            curr_event_t = t > timeavg_t{curr_timeavg}(1) & t < timeavg_t{curr_timeavg}(2);
            curr_ctx_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_ctx_act,'uni',false);
            curr_str_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_str_act,'uni',false);
            
            % (bin measured data across percentile range)
            bin_range = prctile(cell2mat(cellfun(@(x) x(:,plot_ctx),curr_ctx_act_avg,'uni',false)),act_prctile);
            bin_edges = linspace(bin_range(1),bin_range(2),n_act_bins+1);
            bin_centers = bin_edges(1:end-1) + diff(bin_edges)./2;
            
            trial_bins = cellfun(@(x) discretize(x,bin_edges),curr_ctx_act_avg,'uni',false);
            total_bins = cellfun(@(x) discretize(x,[bin_edges(1),bin_edges(end)]),curr_ctx_act_avg,'uni',false);
            
            % (get trials to use: no NaNs in either time series or bins)
            nonan_trials = cellfun(@(ctx_act,str_act,trial_bins) ...
                squeeze(~any(isnan(ctx_act(:,curr_event_t,plot_ctx)),2)) & ...
                squeeze(~any(isnan(str_act(:,curr_event_t,plot_str)),2)) & ...
                ~isnan(trial_bins), ...
                curr_ctx_act,curr_str_act,trial_bins,'uni',false);
            
            % (get average binned activity for measured/predicted by condition)
            ctx_act_binmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_ctx_act_avg,trial_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            str_act_binmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_str_act_avg,trial_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            % (get the average total activity)
            ctx_act_totalmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [1,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_ctx_act_avg,total_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            str_act_totalmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [1,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_str_act_avg,total_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            % Plot binned predicted v measured           
            subplot(length(plot_str_ctx),length(timeavg_labels), ...
                sub2ind(fliplr([length(plot_str_ctx),length(timeavg_labels)]),curr_timeavg,curr_str_ctx));
            hold on;
            col = lines(size(trial_conditions,2));
            switch curr_mua
                case 1
                    errorbar( ...
                        squeeze(nanmean(ctx_act_binmean(:,plot_ctx,:,:),3)), ...
                        squeeze(nanmean(str_act_binmean(:,plot_str,:,:),3)), ...
                        squeeze(AP_sem(str_act_binmean(:,plot_str,:,:),3)), ...
                        'linewidth',2,'CapSize',0);
%                     errorbar( ...
%                         squeeze(nanmean(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
%                         squeeze(nanmean(str_act_totalmean(:,plot_str,:,:),3)), ...
%                         squeeze(AP_sem(str_act_totalmean(:,plot_str,:,:),3)), ...
%                         squeeze(AP_sem(str_act_totalmean(:,plot_str,:,:),3)), ...
%                         squeeze(AP_sem(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
%                         squeeze(AP_sem(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
%                         '.','linewidth',3,'CapSize',0);
                case 2
                    for curr_cond = 1:size(trial_conditions,2)
                        AP_errorfill( ...
                            squeeze(nanmean(ctx_act_binmean(:,plot_ctx,:,curr_cond),3)), ...
                            squeeze(nanmean(str_act_binmean(:,plot_str,:,curr_cond),3)), ...
                            squeeze(AP_sem(str_act_binmean(:,plot_str,:,curr_cond),3)), ...
                            col(curr_cond,:),0.5,false);
                    end
            end
            xlabel(['Ctx (' num2str(plot_ctx) ')']);
            ylabel([mua_label ' Str (' num2str(plot_str) ')'])
            title(timeavg_labels{curr_timeavg});
            
        end
        
    end
end

% Link axes of all plots
linkaxes(get(str_v_ctx_fig,'Children'));


%% MUA (meas/ctx-pred overlay) v fluorescence (task-reduction, epoch-specific condition)

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_contrast_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set windows to average activity
timeavg_labels = {'Stim','Move onset','Outcome'};
timeavg_t = {[0.05,0.15],[-0.05,0.05],[0.05,0.15]};
timeavg_align = {stim_align,move_align,outcome_align};
timeavg_trial_conditions = ...
    {[sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1], ...
    [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
    [trial_outcome_allcat == 1, trial_outcome_allcat == -1]};
timeavg_task_reduction = [1,2,4];

% Set activity percentiles and bins
act_prctile = [10,90];
n_act_bins = 5;

% Set areas and conditions
plot_str_ctx = {[1,1],[2,2],[3,3],[4,4]};

% Loop across area pairs, plot binned predicted v measured activity
str_v_ctx_fig = figure('color','w');
for curr_str_ctx = 1:length(plot_str_ctx)
    
    plot_str = plot_str_ctx{curr_str_ctx}(1);
    plot_ctx = plot_str_ctx{curr_str_ctx}(2);

    for curr_mua = 1:2
        
        % (set striatum activity to use)
        switch curr_mua
            case 1
                curr_str_act_allcat = single(mua_allcat);
                curr_str_act_taskpred_reduced_allcat = single(mua_taskpred_reduced_allcat);
                mua_label = 'Measured';
            case 2
                curr_str_act_allcat = single(mua_ctxpred_allcat);
                curr_str_act_taskpred_reduced_allcat = single(mua_ctxpred_taskpred_reduced_allcat);
                mua_label = 'Ctx-pred';
        end
        
        for curr_timeavg = 1:length(timeavg_labels)
            
            trial_conditions = timeavg_trial_conditions{curr_timeavg};
            curr_task_reduction = timeavg_task_reduction(curr_timeavg);
            
            % (set cortex activity to use)
%             curr_ctx_act_allcat = fluor_roi_deconv;
            curr_ctx_act_allcat = fluor_kernelroi_deconv;   
            curr_ctx_act_taskpred_reduced_allcat = fluor_kernelroi_taskpred_reduced;
            
            % (re-align and split activity and conditions)
            curr_ctx_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_ctx_act_allcat(trial,:,:) - curr_ctx_act_taskpred_reduced_allcat(trial,:,:,curr_task_reduction), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_ctx_act_allcat,1)),'uni',false)), ...
                use_split,length(t),size(curr_ctx_act_allcat,3));
            
            curr_str_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_str_act_allcat(trial,:,:) - curr_str_act_taskpred_reduced_allcat(trial,:,:,curr_task_reduction), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_str_act_allcat,1)),'uni',false)), ...
                use_split,length(t),size(curr_str_act_allcat,3));
            
            trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
            
            % (get average activity within window)
            curr_event_t = t > timeavg_t{curr_timeavg}(1) & t < timeavg_t{curr_timeavg}(2);
            curr_ctx_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_ctx_act,'uni',false);
            curr_str_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_str_act,'uni',false);
            
            % (bin measured data across percentile range)
            bin_range = prctile(cell2mat(cellfun(@(x) x(:,plot_ctx),curr_ctx_act_avg,'uni',false)),act_prctile);
            bin_edges = linspace(bin_range(1),bin_range(2),n_act_bins+1);
            bin_centers = bin_edges(1:end-1) + diff(bin_edges)./2;
            
            trial_bins = cellfun(@(x) discretize(x,bin_edges),curr_ctx_act_avg,'uni',false);
            total_bins = cellfun(@(x) discretize(x,[bin_edges(1),bin_edges(end)]),curr_ctx_act_avg,'uni',false);
            
            % (get trials to use: no NaNs in either time series or bins)
            nonan_trials = cellfun(@(ctx_act,str_act,trial_bins) ...
                squeeze(~any(isnan(ctx_act(:,curr_event_t,plot_ctx)),2)) & ...
                squeeze(~any(isnan(str_act(:,curr_event_t,plot_str)),2)) & ...
                ~isnan(trial_bins), ...
                curr_ctx_act,curr_str_act,trial_bins,'uni',false);
            
            % (get average binned activity for measured/predicted by condition)
            ctx_act_binmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_ctx_act_avg,trial_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            str_act_binmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_str_act_avg,trial_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            % (get the average total activity)
            ctx_act_totalmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [1,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_ctx_act_avg,total_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            str_act_totalmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [1,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_str_act_avg,total_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            % Plot binned predicted v measured           
            subplot(length(plot_str_ctx),length(timeavg_labels), ...
                sub2ind(fliplr([length(plot_str_ctx),length(timeavg_labels)]),curr_timeavg,curr_str_ctx));
            hold on;
            col = lines(size(trial_conditions,2));
            switch curr_mua
                case 1
                    errorbar( ...
                        squeeze(nanmean(ctx_act_binmean(:,plot_ctx,:,:),3)), ...
                        squeeze(nanmean(str_act_binmean(:,plot_str,:,:),3)), ...
                        squeeze(AP_sem(str_act_binmean(:,plot_str,:,:),3)), ...
                        'linewidth',2,'CapSize',0);
%                     errorbar( ...
%                         squeeze(nanmean(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
%                         squeeze(nanmean(str_act_totalmean(:,plot_str,:,:),3)), ...
%                         squeeze(AP_sem(str_act_totalmean(:,plot_str,:,:),3)), ...
%                         squeeze(AP_sem(str_act_totalmean(:,plot_str,:,:),3)), ...
%                         squeeze(AP_sem(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
%                         squeeze(AP_sem(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
%                         '.','linewidth',3,'CapSize',0);
                case 2
                    for curr_cond = 1:size(trial_conditions,2)
                        AP_errorfill( ...
                            squeeze(nanmean(ctx_act_binmean(:,plot_ctx,:,curr_cond),3)), ...
                            squeeze(nanmean(str_act_binmean(:,plot_str,:,curr_cond),3)), ...
                            squeeze(AP_sem(str_act_binmean(:,plot_str,:,curr_cond),3)), ...
                            col(curr_cond,:),0.5,false);
                    end
            end
            xlabel(['Ctx (' num2str(plot_ctx) ')']);
            ylabel([' Str (' num2str(plot_str) ')'])
            title([timeavg_labels{curr_timeavg} '(' task_regressor_labels{curr_task_reduction} '-reduced)']);
            
        end
        
    end
end

% Link axes of all plots
linkaxes(get(str_v_ctx_fig,'Children'));



%% MUA (passive) v fluorescence by stim (overlay measured/ctx-pred)

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(fluor_allcat,1),1);

% Set windows to average activity
timeavg_labels = {'Pre-stim','Stim'};
timeavg_t = {[-0.2,-0.1],[0.05,0.15]};
timeavg_align = {stim_align,stim_align};

% Set activity percentiles and bins
act_prctile = [10,90];
n_act_bins = 5;

% Set areas and conditions
% plot_str_ctx = {[1,3],[2,7]};
plot_str_ctx = {[1,1],[2,2]};
% trial_condition_groups = ...
%     {[sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1], ...
%     [sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1]};
% trial_condition_groups = ...
%     {[sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1], ...
%     [sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1]};
trial_condition_groups = ...
    {[trial_stim_allcat == 1, trial_stim_allcat == 2, trial_stim_allcat == 3], ...
    [trial_stim_allcat == 1, trial_stim_allcat == 2, trial_stim_allcat == 3]};

% Loop across area pairs, plot binned predicted v measured activity
str_v_ctx_fig = figure('color','w');
for curr_str_ctx = 1:length(plot_str_ctx)
    
    plot_str = plot_str_ctx{curr_str_ctx}(1);
    plot_ctx = plot_str_ctx{curr_str_ctx}(2);
    trial_conditions = trial_condition_groups{curr_str_ctx};

    for curr_mua = 1:2
        
        % (set striatum activity to use)
        switch curr_mua
            case 1
                curr_str_act_allcat = single(mua_allcat);
                mua_label = 'Measured';
            case 2
                curr_str_act_allcat = single(mua_ctxpred_allcat);
                mua_label = 'Ctx-pred';
        end
        
        for curr_timeavg = 1:length(timeavg_labels)
            
            % (set cortex activity to use)
%             curr_ctx_act_allcat = fluor_roi_deconv;    
            curr_ctx_act_allcat = fluor_kernelroi_deconv;    
            
            % (re-align and split activity and conditions)
            curr_ctx_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift(curr_ctx_act_allcat(trial,:,:), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_ctx_act_allcat,1)),'uni',false)), ...
                use_split,length(t),size(curr_ctx_act_allcat,3));
            
            curr_str_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift(curr_str_act_allcat(trial,:,:), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_str_act_allcat,1)),'uni',false)), ...
                use_split,length(t),size(curr_str_act_allcat,3));
            
            trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
            
            % (get average activity within window)
            curr_event_t = t > timeavg_t{curr_timeavg}(1) & t < timeavg_t{curr_timeavg}(2);
            curr_ctx_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_ctx_act,'uni',false);
            curr_str_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_str_act,'uni',false);
            
            % (bin measured data across percentile range)
            bin_range = prctile(cell2mat(cellfun(@(x) x(:,plot_ctx),curr_ctx_act_avg,'uni',false)),act_prctile);
            bin_edges = linspace(bin_range(1),bin_range(2),n_act_bins+1);
            bin_centers = bin_edges(1:end-1) + diff(bin_edges)./2;
            
            trial_bins = cellfun(@(x) discretize(x,bin_edges),curr_ctx_act_avg,'uni',false);
            total_bins = cellfun(@(x) discretize(x,[bin_edges(1),bin_edges(end)]),curr_ctx_act_avg,'uni',false);
            
            % (get trials to use: no NaNs in either time series or bins)
            nonan_trials = cellfun(@(ctx_act,str_act,trial_bins) ...
                squeeze(~any(isnan(ctx_act(:,curr_event_t,plot_ctx)),2)) & ...
                squeeze(~any(isnan(str_act(:,curr_event_t,plot_str)),2)) & ...
                ~isnan(trial_bins), ...
                curr_ctx_act,curr_str_act,trial_bins,'uni',false);
            
            % (get average binned activity for measured/predicted by condition)
            ctx_act_binmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_ctx_act_avg,trial_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            str_act_binmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_str_act_avg,trial_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            % (get the average total activity)
            ctx_act_totalmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [1,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_ctx_act_avg,total_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            str_act_totalmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [1,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_str_act_avg,total_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            % Plot binned predicted v measured           
            subplot(length(plot_str_ctx),length(timeavg_labels), ...
                sub2ind(fliplr([length(plot_str_ctx),length(timeavg_labels)]),curr_timeavg,curr_str_ctx));
            hold on;
            col = lines(size(trial_conditions,2));
            switch curr_mua
                case 1
                    errorbar( ...
                        squeeze(nanmean(ctx_act_binmean(:,plot_ctx,:,:),3)), ...
                        squeeze(nanmean(str_act_binmean(:,plot_str,:,:),3)), ...
                        squeeze(AP_sem(str_act_binmean(:,plot_str,:,:),3)), ...
                        'linewidth',2,'CapSize',0);
%                     errorbar( ...
%                         squeeze(nanmean(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
%                         squeeze(nanmean(str_act_totalmean(:,plot_str,:,:),3)), ...
%                         squeeze(AP_sem(str_act_totalmean(:,plot_str,:,:),3)), ...
%                         squeeze(AP_sem(str_act_totalmean(:,plot_str,:,:),3)), ...
%                         squeeze(AP_sem(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
%                         squeeze(AP_sem(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
%                         '.','linewidth',3,'CapSize',0);
                case 2
                    for curr_cond = 1:size(trial_conditions,2)
                        AP_errorfill( ...
                            squeeze(nanmean(ctx_act_binmean(:,plot_ctx,:,curr_cond),3)), ...
                            squeeze(nanmean(str_act_binmean(:,plot_str,:,curr_cond),3)), ...
                            squeeze(AP_sem(str_act_binmean(:,plot_str,:,curr_cond),3)), ...
                            col(curr_cond,:),0.5,false);
                    end
            end
            xlabel(['Ctx (' wf_roi(plot_ctx).area ')']);
            ylabel([mua_label ' Str (' num2str(plot_str) ')'])
            title(timeavg_labels{curr_timeavg});
            
        end
        
    end
    % Link axes of all plots
    linkaxes(get(str_v_ctx_fig,'Children'));
end


%% Measured vs predicted (across-epochs conditions)

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_contrast_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set windows to average activity
timeavg_labels = {'Pre-stim','Stim','Move onset','Outcome'};
timeavg_t = {[-0.2,-0.1],[0.05,0.15],[-0.05,0.05],[0.05,0.15]};
timeavg_align = {stim_align,stim_align,move_align,outcome_align};

% timeavg_labels = {'Pre-stim','Stim'};
% timeavg_t = {[-0.2,-0.1],[0.05,0.15]};
% timeavg_align = {stim_align,stim_align};

% Set activity percentiles and bins
act_prctile = [10,90];
n_act_bins = 5;

% Set areas and conditions

% plot_areas = [1,3,7,10];
plot_areas = [1,2,3,4];

trial_condition_groups = ...
    {[sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1], ...
    [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
    [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
    [trial_outcome_allcat == 1, trial_outcome_allcat == -1]};

% trial_condition_groups = ...
%     {[sign(trial_contrastside_allcat) == 1 & trial_choice_allcat == -1, ...
%     sign(trial_contrastside_allcat) == 1 & trial_choice_allcat == 1], ...
%     ...
%     [trial_choice_allcat == -1 & sign(trial_contrastside_allcat) == 1, ... 
%     trial_choice_allcat == -1 & sign(trial_contrastside_allcat) == 0, ...
%     trial_choice_allcat == -1 & sign(trial_contrastside_allcat) == -1], ...
%     ...
%     [trial_choice_allcat == -1 & sign(trial_contrastside_allcat) == 1, ... 
%     trial_choice_allcat == -1 & sign(trial_contrastside_allcat) == 0, ...
%     trial_choice_allcat == -1 & sign(trial_contrastside_allcat) == -1], ...
%     ...
%     [trial_choice_allcat == -1 & trial_outcome_allcat == 1, ... 
%     trial_choice_allcat == 1 & trial_outcome_allcat == 1]};

% trial_condition_groups = ...
%     {[sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1], ...
%     [sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1], ...
%     [sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1], ...
%     [sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1]};

% trial_condition_groups = ...
%     {[trial_stim_allcat == 1, trial_stim_allcat == 2, trial_stim_allcat == 3], ...
%     [trial_stim_allcat == 1, trial_stim_allcat == 2, trial_stim_allcat == 3], ...
%     [trial_stim_allcat == 1, trial_stim_allcat == 2, trial_stim_allcat == 3], ...
%     [trial_stim_allcat == 1, trial_stim_allcat == 2, trial_stim_allcat == 3]};

% Loop across area pairs, plot binned predicted v measured activity
% (cortex)
curr_act_allcat = fluor_roi_deconv;
curr_act_pred_allcat = fluor_roi_taskpred;
% (striatum)
% curr_act_allcat = mua_allcat - mua_taskpred_reduced_allcat(:,:,:,1);
% % curr_act_pred_allcat = mua_taskpred_allcat;
% curr_act_pred_allcat = mua_ctxpred_allcat - mua_ctxpred_taskpred_reduced_allcat(:,:,:,1);

measured_v_pred_fig = figure('color','w');
for curr_area_idx = 1:length(plot_areas)
    
    plot_area = plot_areas(curr_area_idx);
    trial_conditions = trial_condition_groups{curr_area_idx};
    
    % Set up the summary values
    curr_act_pred_diff = nan(size(trial_conditions,2),length(use_split),length(timeavg_labels));
    curr_act_pred_rank_diff = nan(size(trial_conditions,2),length(use_split),length(timeavg_labels));
    n_shuff = 1000;
    curr_act_pred_diff_shuff = nan(size(trial_conditions,2),length(use_split),n_shuff,length(timeavg_labels));
    curr_expl_var = nan(size(trial_conditions,2),length(use_split),length(timeavg_labels));
    
    for curr_timeavg = 1:length(timeavg_labels)
        
        % (re-align and split activity)
        curr_act = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift(curr_act_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        curr_act_pred = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift(curr_act_pred_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
        
        % (get average activity within window)
        curr_event_t = t > timeavg_t{curr_timeavg}(1) & t < timeavg_t{curr_timeavg}(2);
        curr_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act,'uni',false);
        curr_act_pred_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_pred,'uni',false);

        % (bin predicted data across percentile range)
        bin_range = prctile(cell2mat(curr_act_pred_avg),act_prctile);
        bin_edges = linspace(bin_range(1),bin_range(2),n_act_bins+1);
        bin_centers = bin_edges(1:end-1) + diff(bin_edges)./2;
        
        trial_bins = cellfun(@(x) discretize(x,bin_edges),curr_act_pred_avg,'uni',false);
        total_bins = cellfun(@(x) discretize(x,[bin_edges(1),bin_edges(end)]),curr_act_pred_avg,'uni',false);
        
        % (get trials to use: no NaNs in either time series or bins)
        use_trials = cellfun(@(act,act_taskpred,trial_bins) ...
            squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
            squeeze(~any(isnan(act_taskpred(:,curr_event_t)),2)) & ...
            ~isnan(trial_bins), ...
            curr_act,curr_act_pred,trial_bins,'uni',false);
        
        % (get average binned activity)
        act_binmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,trial_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_binmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,trial_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));      

        % (get the average total activity)       
        act_totalmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [1,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,total_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_totalmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [1,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,total_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        % Get average act-pred difference (amplitude and rank) (ALL TRIALS)
        curr_act_pred_diff(:,:,curr_timeavg) = ...
            cell2mat(arrayfun(@(cond) cellfun(@(act,pred,trial_cond,use_trials) ...
            nanmean(act(trial_cond(:,cond)) - ...
            pred(trial_cond(:,cond))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp,use_trials), ...
            1:size(trial_conditions,2),'uni',false))';
        
        curr_act_avg_rank = cellfun(@(x) tiedrank(x)./max(tiedrank(x)),curr_act_avg,'uni',false);
        curr_act_pred_avg_rank = cellfun(@(x) tiedrank(x)./max(tiedrank(x)),curr_act_pred_avg,'uni',false);
        curr_act_pred_rank_diff(:,:,curr_timeavg) = ...
            cell2mat(arrayfun(@(cond) cellfun(@(act,pred,trial_cond,use_trials) ...
            nanmean(act(trial_cond(:,cond)) - ...
            pred(trial_cond(:,cond))), ...
            curr_act_avg_rank,curr_act_pred_avg_rank,trial_conditions_exp,use_trials), ...
            1:size(trial_conditions,2),'uni',false))';
        
        % Get condition-shuffled act-pred difference
        % (unused at the moment)
        trial_conditions_exp_shuff = mat2cell(AP_shake(repmat(trial_conditions,1,1,n_shuff),2), ...
            use_split,size(trial_conditions,2),n_shuff);        
        curr_act_pred_diff_shuff(:,:,:,curr_timeavg) = ...
            cell2mat(arrayfun(@(shuff) ...
            cell2mat(arrayfun(@(cond) cellfun(@(act,pred,trial_cond,use_trials) ...
            nanmean(act(trial_cond(:,cond,shuff)) - ...
            pred(trial_cond(:,cond,shuff))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp_shuff,use_trials), ...
            1:size(trial_conditions,2),'uni',false))',permute(1:n_shuff,[1,3,2]),'uni',false));
        
        % Get the explained variance by condition
        curr_expl_var(:,:,curr_timeavg) = ...
            cell2mat(cellfun(@(act,pred,trial_cond,use_trials) cell2mat(arrayfun(@(condition) ...
            1 - (nansum((act(use_trials & trial_cond(:,condition)) - ...
            pred(use_trials & trial_cond(:,condition))).^2)./ ...
            (nansum((act(use_trials & trial_cond(:,condition)) - ...
            nanmean(act(use_trials & trial_cond(:,condition)),1)).^2,1))), ...
            [1:size(trial_conditions,2)]','uni',false)), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp,use_trials,'uni',false)');
        
        % Plot binned predicted v measured
        figure(measured_v_pred_fig);
        subplot(length(plot_areas),length(timeavg_labels)+3, ...
            sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+3]),curr_timeavg,curr_area_idx));
        hold on;
        
        errorbar( ...
            squeeze(nanmean(act_pred_binmean,2)), ...
            squeeze(nanmean(act_binmean,2)), ...
            squeeze(AP_sem(act_binmean,2)), ...
            'linewidth',2,'CapSize',0);
        errorbar( ...
            squeeze(nanmean(act_pred_totalmean,2)), ...
            squeeze(nanmean(act_totalmean,2)), ...
            squeeze(AP_sem(act_totalmean,2)), ...
            squeeze(AP_sem(act_totalmean,2)), ...
            squeeze(AP_sem(act_pred_totalmean,2)), ...
            squeeze(AP_sem(act_pred_totalmean,2)), ...
            '.','linewidth',3,'CapSize',0);
        xlabel(['Predicted (' num2str(plot_area) ')']);
        ylabel(['Measured (' num2str(plot_area) ')'])
        ylim(xlim); axis square;
        title(timeavg_labels{curr_timeavg});
        
    end
    
    % Plot measured - predicted (amplitude and rank)
    subplot(length(plot_areas),length(timeavg_labels)+3, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+3]),length(timeavg_labels)+1,curr_area_idx));
     hold on;
    errorbar( ...
        permute(nanmean(curr_act_pred_diff,2),[3,1,2]), ...
        permute(AP_sem(curr_act_pred_diff,2),[3,1,2]),'linewidth',2,'CapSize',0);
    ylabel('Act-pred');
    set(gca,'XTick',1:4,'XTickLabels',timeavg_labels,'XTickLabelRotation',45)
    xlim([0.5,length(timeavg_labels)+0.5]);
    
    subplot(length(plot_areas),length(timeavg_labels)+3, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+3]),length(timeavg_labels)+2,curr_area_idx));
     hold on;
    errorbar( ...
        permute(nanmean(curr_act_pred_rank_diff,2),[3,1,2]), ...
        permute(AP_sem(curr_act_pred_rank_diff,2),[3,1,2]),'linewidth',2,'CapSize',0);
    ylabel('Act-pred (rank)');
    set(gca,'XTick',1:4,'XTickLabels',timeavg_labels,'XTickLabelRotation',45)
    xlim([0.5,length(timeavg_labels)+0.5]);
    
    % Get/plot measured - predicted significance
    subplot(length(plot_areas),length(timeavg_labels)+3, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+3]),length(timeavg_labels)+1,curr_area_idx));
    real_diff = squeeze(nanmean(curr_act_pred_diff(1,:,:) - curr_act_pred_diff(2,:,:),2));
    shuff_diff_ci = permute(prctile(nanmean(curr_act_pred_diff_shuff(1,:,:,:) - ...
        curr_act_pred_diff_shuff(2,:,:,:),2),[2.5,97.5],3),[4,3,2,1]);
    sig_diff = real_diff < shuff_diff_ci(:,1) | real_diff > shuff_diff_ci(:,2);
    plot(find(sig_diff),repmat(max(ylim),sum(sig_diff),1),'*','color','k','MarkerSize',5)
  
    % Plot explained variance
    subplot(length(plot_areas),length(timeavg_labels)+3, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+3]),length(timeavg_labels)+3,curr_area_idx));
    hold on;
    errorbar( ...
        permute(nanmean(curr_expl_var,2),[3,1,2]), ...
        permute(AP_sem(curr_expl_var,2),[3,1,2]),'linewidth',2,'CapSize',0);
    ylabel('R^2');
    set(gca,'XTick',1:4,'XTickLabels',timeavg_labels,'XTickLabelRotation',45)
    xlim([0.5,length(timeavg_labels)+0.5]);
    
end

% Link axes of all predicted v measured and all explained var
all_axes = get(measured_v_pred_fig,'Children');
meas_pred_axes = all_axes(setdiff(1:length(all_axes), ...
    [1:length(timeavg_labels)+3:length(all_axes), ...
    2:length(timeavg_labels)+3:length(all_axes), ...
    3:length(timeavg_labels)+3:length(all_axes)]));
linkaxes(meas_pred_axes)
linkaxes(all_axes(intersect(1:length(all_axes),1:length(timeavg_labels)+3:length(all_axes))));
linkaxes(all_axes(intersect(1:length(all_axes),2:length(timeavg_labels)+3:length(all_axes))));
linkaxes(all_axes(intersect(1:length(all_axes),3:length(timeavg_labels)+3:length(all_axes))));

for curr_axes = 1:length(meas_pred_axes)
   line(meas_pred_axes(curr_axes),xlim(meas_pred_axes(curr_axes)), ...
       xlim(meas_pred_axes(curr_axes)),'color','k'); 
end


%% Measured vs predicted (task-reduction, epoch-specific condition)

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_contrast_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set windows to average activity
timeavg_labels = {'Stim','Move onset','Outcome'};
timeavg_t = {[0.05,0.15],[-0.05,0.05],[0.05,0.15]};
timeavg_align = {stim_align,move_align,outcome_align};
timeavg_trial_conditions = ...
    {[sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1], ...
    [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
    [trial_outcome_allcat == 1, trial_outcome_allcat == -1]};
timeavg_task_reduction = [1,2,4];

% timeavg_labels = {'Pre-stim','Stim'};
% timeavg_t = {[-0.2,-0.1],[0.05,0.15]};
% timeavg_align = {stim_align,stim_align};

% Set activity percentiles and bins
act_prctile = [10,90];
n_act_bins = 5;

% Set areas and conditions
plot_areas = [1,2,3,4];

% Loop across area pairs, plot binned predicted v measured activity
curr_act_allcat = mua_allcat;
curr_act_taskpred_reduced_allcat = mua_taskpred_reduced_allcat;
curr_act_pred_allcat = mua_ctxpred_allcat;
curr_act_pred_taskpred_reduced_allcat = mua_ctxpred_taskpred_reduced_allcat;

measured_v_pred_fig = figure('color','w');
for curr_area_idx = 1:length(plot_areas)
    
    plot_area = plot_areas(curr_area_idx);
    
    % Set up the summary values
    curr_act_pred_diff = nan(2,length(use_split),length(timeavg_labels));
    curr_act_pred_rank_condition_diff = nan(length(use_split),length(timeavg_labels));
    n_shuff = 1000;
    curr_act_pred_rank_condition_diff_shuff = nan(length(use_split),n_shuff,length(timeavg_labels));
    
    for curr_timeavg = 1:length(timeavg_labels)
        
        trial_conditions = timeavg_trial_conditions{curr_timeavg};
        curr_task_reduction = timeavg_task_reduction(curr_timeavg);
        
        % (re-align and split activity)
        curr_act = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_allcat(trial,:,plot_area) - curr_act_taskpred_reduced_allcat(trial,:,plot_area,curr_task_reduction), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        curr_act_pred = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_pred_allcat(trial,:,plot_area) - curr_act_pred_taskpred_reduced_allcat(trial,:,plot_area,curr_task_reduction), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
        
        % (get average activity within window)
        curr_event_t = t > timeavg_t{curr_timeavg}(1) & t < timeavg_t{curr_timeavg}(2);
        curr_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act,'uni',false);
        curr_act_pred_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_pred,'uni',false);

        % (bin predicted data across percentile range)
        bin_range = prctile(cell2mat(curr_act_pred_avg),act_prctile);
        bin_edges = linspace(bin_range(1),bin_range(2),n_act_bins+1);
        bin_centers = bin_edges(1:end-1) + diff(bin_edges)./2;
        
        trial_bins = cellfun(@(x) discretize(x,bin_edges),curr_act_pred_avg,'uni',false);
        total_bins = cellfun(@(x) discretize(x,[bin_edges(1),bin_edges(end)]),curr_act_pred_avg,'uni',false);
        
        % (get trials to use: no NaNs in either time series or bins)
        use_trials = cellfun(@(act,act_taskpred,trial_bins) ...
            squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
            squeeze(~any(isnan(act_taskpred(:,curr_event_t)),2)) & ...
            ~isnan(trial_bins), ...
            curr_act,curr_act_pred,trial_bins,'uni',false);
        
        % (get average binned activity)
        act_binmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,trial_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_binmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,trial_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));      

        % (get the average total activity)       
        act_totalmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [1,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,total_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_totalmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [1,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,total_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        % Get average act-pred difference (ALL TRIALS)
        curr_act_pred_diff(:,:,curr_timeavg) = ...
            cell2mat(arrayfun(@(cond) cellfun(@(act,pred,trial_cond,use_trials) ...
            nanmean(act(trial_cond(:,cond)) - pred(trial_cond(:,cond))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp,use_trials), ...
            1:size(trial_conditions,2),'uni',false))';
        
        % Get condition  difference (assume 2 conditions)
        curr_act_avg_rank = cellfun(@(x) tiedrank(x)./max(tiedrank(x)),curr_act_avg,'uni',false);
        curr_act_pred_avg_rank = cellfun(@(x) tiedrank(x)./max(tiedrank(x)),curr_act_pred_avg,'uni',false);        
        curr_act_pred_condition_diff(:,curr_timeavg) = ...
            cellfun(@(act,pred,trial_cond,use_trials) ...       
            nanmean(act(trial_cond(:,1)) - pred(trial_cond(:,1))) - ...
            nanmean(act(trial_cond(:,2)) - pred(trial_cond(:,2))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp,use_trials);
        
        % Get condition-shuffled act-pred difference
        trial_conditions_exp_shuff = cellfun(@(x) AP_shake(repmat(x,1,1,n_shuff),1), ...
            trial_conditions_exp,'uni',false);   
                
        curr_act_pred_condition_diff_shuff(:,:,curr_timeavg) = ...
            cell2mat(arrayfun(@(shuff) ...
            cellfun(@(act,pred,trial_cond,use_trials) ...
            nanmean(act(trial_cond(:,1,shuff)) - pred(trial_cond(:,1,shuff))) - ...
            nanmean(act(trial_cond(:,2,shuff)) - pred(trial_cond(:,2,shuff))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp_shuff,use_trials), ...
            1:n_shuff,'uni',false));
        
        % Plot binned predicted v measured
        figure(measured_v_pred_fig);
        subplot(length(plot_areas),length(timeavg_labels)+2, ...
            sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),curr_timeavg,curr_area_idx));
        hold on;
        
        errorbar( ...
            squeeze(nanmean(act_pred_binmean,2)), ...
            squeeze(nanmean(act_binmean,2)), ...
            squeeze(AP_sem(act_binmean,2)), ...
            'linewidth',2,'CapSize',0);
        errorbar( ...
            squeeze(nanmean(act_pred_totalmean,2)), ...
            squeeze(nanmean(act_totalmean,2)), ...
            squeeze(AP_sem(act_totalmean,2)), ...
            squeeze(AP_sem(act_totalmean,2)), ...
            squeeze(AP_sem(act_pred_totalmean,2)), ...
            squeeze(AP_sem(act_pred_totalmean,2)), ...
            '.','linewidth',3,'CapSize',0);
        xlabel(['Predicted (' num2str(plot_area) ')']);
        ylabel(['Measured (' num2str(plot_area) ')'])
        ylim(xlim); axis square;
        title([timeavg_labels{curr_timeavg} ' (' task_regressor_labels{curr_task_reduction} '-reduced)']);
        
    end
    
    % Plot measured - predicted
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+1,curr_area_idx));
     hold on;
    errorbar( ...
        permute(nanmean(curr_act_pred_diff,2),[3,1,2]), ...
        permute(AP_sem(curr_act_pred_diff,2),[3,1,2]),'linewidth',2,'CapSize',0);
    ylabel('Meas - Pred');
    set(gca,'XTick',1:4,'XTickLabels',timeavg_labels,'XTickLabelRotation',45)
    xlim([0.5,length(timeavg_labels)+0.5]);
    
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+2,curr_area_idx));
     hold on;
    errorbar( ...
        nanmean(curr_act_pred_condition_diff,1), ...
        AP_sem(curr_act_pred_condition_diff,1),'linewidth',2,'CapSize',0);
    ylabel('Condition difference');
    set(gca,'XTick',1:4,'XTickLabels',timeavg_labels,'XTickLabelRotation',45)
    xlim([0.5,length(timeavg_labels)+0.5]);
    
    % Get and plot significance
    real_diff = nanmean(curr_act_pred_condition_diff,1); 
    alpha = [5/2/3,100-(5/2/3)];
    shuff_diff_ci = permute(nanmean(prctile(curr_act_pred_condition_diff_shuff,alpha,2),1),[2,3,1]);
    sig_diff = real_diff < shuff_diff_ci(1,:) | real_diff > shuff_diff_ci(2,:);
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+1,curr_area_idx));
    plot(find(sig_diff),repmat(max(ylim),sum(sig_diff),1),'*','color','k','MarkerSize',5);
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+2,curr_area_idx));
    plot(shuff_diff_ci','r');
    plot(find(sig_diff),repmat(max(ylim),sum(sig_diff),1),'*','color','k','MarkerSize',5);
    
end

% Link axes of all predicted v measured and all explained var
all_axes = get(measured_v_pred_fig,'Children');
meas_pred_axes = all_axes(setdiff(1:length(all_axes), ...
    [1:length(timeavg_labels)+2:length(all_axes), ...
    2:length(timeavg_labels)+2:length(all_axes)]));
linkaxes(meas_pred_axes)
linkaxes(all_axes(1:length(timeavg_labels)+2:length(all_axes)));
linkaxes(all_axes(2:length(timeavg_labels)+2:length(all_axes)));

for curr_axes = 1:length(meas_pred_axes)
   line(meas_pred_axes(curr_axes),xlim(meas_pred_axes(curr_axes)), ...
       xlim(meas_pred_axes(curr_axes)),'color','k'); 
end
for curr_axes = 1:length(timeavg_labels)+2:length(all_axes)
   line(all_axes(curr_axes),xlim(all_axes(curr_axes)),[0,0],'color','k');
end
for curr_axes = 2:length(timeavg_labels)+2:length(all_axes)
   line(all_axes(curr_axes),xlim(all_axes(curr_axes)),[0,0],'color','k');
end

%% Measured vs predicted (passive)

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(mua_allcat,1),1);

% Set windows to average activity
timeavg_labels = {'Stim'};
timeavg_t = {[0.05,0.15]};
timeavg_align = {stim_align};
timeavg_trial_conditions = ...
    {[sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1]};
% timeavg_trial_conditions = ...
%     {[trial_stim_allcat == 1,trial_stim_allcat == 2,trial_stim_allcat == 3]};

% Set activity percentiles and bins
act_prctile = [10,90];
n_act_bins = 5;

% Set areas and conditions
plot_areas = [1,2,3,4];

% Loop across area pairs, plot binned predicted v measured activity
curr_act_allcat = mua_allcat;
curr_act_pred_allcat = mua_ctxpred_allcat;

measured_v_pred_fig = figure('color','w');
for curr_area_idx = 1:length(plot_areas)
    
    plot_area = plot_areas(curr_area_idx);
    
    % Set up the summary values
    curr_act_pred_diff = nan(size(timeavg_trial_conditions{1},2),length(use_split),length(timeavg_labels));
    curr_act_pred_rank_condition_diff = nan(length(use_split),length(timeavg_labels));
    n_shuff = 1000;
    curr_act_pred_rank_condition_diff_shuff = nan(length(use_split),n_shuff,length(timeavg_labels));
    
    for curr_timeavg = 1:length(timeavg_labels)
        
        trial_conditions = timeavg_trial_conditions{curr_timeavg};
        
        % (re-align and split activity)
        curr_act = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        curr_act_pred = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_pred_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
        
        % (get average activity within window)
        curr_event_t = t > timeavg_t{curr_timeavg}(1) & t < timeavg_t{curr_timeavg}(2);
        curr_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act,'uni',false);
        curr_act_pred_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_pred,'uni',false);

        % (bin predicted data across percentile range)
        bin_range = prctile(cell2mat(curr_act_pred_avg),act_prctile);
        bin_edges = linspace(bin_range(1),bin_range(2),n_act_bins+1);
        bin_centers = bin_edges(1:end-1) + diff(bin_edges)./2;
        
        trial_bins = cellfun(@(x) discretize(x,bin_edges),curr_act_pred_avg,'uni',false);
        total_bins = cellfun(@(x) discretize(x,[bin_edges(1),bin_edges(end)]),curr_act_pred_avg,'uni',false);
        
        % (get trials to use: no NaNs in either time series or bins)
        use_trials = cellfun(@(act,act_taskpred,trial_bins) ...
            squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
            squeeze(~any(isnan(act_taskpred(:,curr_event_t)),2)) & ...
            ~isnan(trial_bins), ...
            curr_act,curr_act_pred,trial_bins,'uni',false);
        
        % (get average binned activity)
        act_binmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,trial_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_binmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,trial_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));      

        % (get the average total activity)       
        act_totalmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [1,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,total_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_totalmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [1,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,total_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        % Get average act-pred difference (ALL TRIALS)
        curr_act_pred_diff(:,:,curr_timeavg) = ...
            cell2mat(arrayfun(@(cond) cellfun(@(act,pred,trial_cond,use_trials) ...
            nanmean(act(trial_cond(:,cond)) - pred(trial_cond(:,cond))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp,use_trials), ...
            1:size(trial_conditions,2),'uni',false))';
        
        % Get condition  difference (assume 2 conditions)
        curr_act_avg_rank = cellfun(@(x) tiedrank(x)./max(tiedrank(x)),curr_act_avg,'uni',false);
        curr_act_pred_avg_rank = cellfun(@(x) tiedrank(x)./max(tiedrank(x)),curr_act_pred_avg,'uni',false);        
        curr_act_pred_condition_diff(:,curr_timeavg) = ...
            cellfun(@(act,pred,trial_cond,use_trials) ...       
            nanmean(act(trial_cond(:,1)) - pred(trial_cond(:,1))) - ...
            nanmean(act(trial_cond(:,2)) - pred(trial_cond(:,2))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp,use_trials);
        
        % Get condition-shuffled act-pred difference
        trial_conditions_exp_shuff = cellfun(@(x) AP_shake(repmat(x,1,1,n_shuff),1), ...
            trial_conditions_exp,'uni',false);   
                
        curr_act_pred_condition_diff_shuff(:,:,curr_timeavg) = ...
            cell2mat(arrayfun(@(shuff) ...
            cellfun(@(act,pred,trial_cond,use_trials) ...
            nanmean(act(trial_cond(:,1,shuff)) - pred(trial_cond(:,1,shuff))) - ...
            nanmean(act(trial_cond(:,2,shuff)) - pred(trial_cond(:,2,shuff))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp_shuff,use_trials), ...
            1:n_shuff,'uni',false));
        
        % Plot binned predicted v measured
        figure(measured_v_pred_fig);
        subplot(length(plot_areas),length(timeavg_labels)+2, ...
            sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),curr_timeavg,curr_area_idx));
        hold on;
        
        errorbar( ...
            squeeze(nanmean(act_pred_binmean,2)), ...
            squeeze(nanmean(act_binmean,2)), ...
            squeeze(AP_sem(act_binmean,2)), ...
            'linewidth',2,'CapSize',0);
        errorbar( ...
            squeeze(nanmean(act_pred_totalmean,2)), ...
            squeeze(nanmean(act_totalmean,2)), ...
            squeeze(AP_sem(act_totalmean,2)), ...
            squeeze(AP_sem(act_totalmean,2)), ...
            squeeze(AP_sem(act_pred_totalmean,2)), ...
            squeeze(AP_sem(act_pred_totalmean,2)), ...
            '.','linewidth',3,'CapSize',0);
        xlabel(['Predicted (' num2str(plot_area) ')']);
        ylabel(['Measured (' num2str(plot_area) ')'])
        ylim(xlim); axis square;
        title([timeavg_labels{curr_timeavg}]);
        
    end
    
    % Plot measured - predicted
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+1,curr_area_idx));
     hold on;
    errorbar( ...
        permute(nanmean(curr_act_pred_diff,2),[3,1,2]), ...
        permute(AP_sem(curr_act_pred_diff,2),[3,1,2]),'linewidth',2,'CapSize',0);
    ylabel('Meas - Pred');
    xlabel('Stim');
    set(gca,'XTick',1:size(timeavg_trial_conditions{1},2));
    xlim([0.5,size(timeavg_trial_conditions{1},2)+0.5]);
    
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+2,curr_area_idx));
     hold on;
    errorbar( ...
        nanmean(curr_act_pred_condition_diff,1), ...
        AP_sem(curr_act_pred_condition_diff,1),'linewidth',2,'CapSize',0);
    ylabel('Condition difference');
    xlabel('Stim');
    set(gca,'XTick',1:size(timeavg_trial_conditions{1},2));
    xlim([0.5,size(timeavg_trial_conditions{1},2)+0.5]);
    
    % Get and plot significance
    real_diff = nanmean(curr_act_pred_condition_diff,1); 
    alpha = [5/2/3,100-(5/2/3)];
    shuff_diff_ci = permute(nanmean(prctile(curr_act_pred_condition_diff_shuff,alpha,2),1),[2,3,1]);
    sig_diff = real_diff < shuff_diff_ci(1,:) | real_diff > shuff_diff_ci(2,:);
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+1,curr_area_idx));
    plot(find(sig_diff),repmat(max(ylim),sum(sig_diff),1),'*','color','k','MarkerSize',5);
    
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+2,curr_area_idx));
    plot(shuff_diff_ci','r');
    plot(find(sig_diff),repmat(max(ylim),sum(sig_diff),1),'*','color','k','MarkerSize',5);
    
end

% Link axes of all predicted v measured and all explained var
all_axes = get(measured_v_pred_fig,'Children');
meas_pred_axes = all_axes(setdiff(1:length(all_axes), ...
    [1:length(timeavg_labels)+2:length(all_axes), ...
    2:length(timeavg_labels)+2:length(all_axes)]));
linkaxes(meas_pred_axes)
linkaxes(all_axes(1:length(timeavg_labels)+2:length(all_axes)));
linkaxes(all_axes(2:length(timeavg_labels)+2:length(all_axes)));

for curr_axes = 1:length(meas_pred_axes)
   line(meas_pred_axes(curr_axes),xlim(meas_pred_axes(curr_axes)), ...
       xlim(meas_pred_axes(curr_axes)),'color','k'); 
end
for curr_axes = 1:length(timeavg_labels)+2:length(all_axes)
   line(all_axes(curr_axes),xlim(all_axes(curr_axes)),[0,0],'color','k');
end
for curr_axes = 2:length(timeavg_labels)+2:length(all_axes)
   line(all_axes(curr_axes),xlim(all_axes(curr_axes)),[0,0],'color','k');
end

%% Str/Ctx->Str nonlinearity (concatenated)

% Apply empirical static nonlinearity
figure;
mua_allcat_predicted_nlin = mua_allcat;
for curr_depth = 1:n_depths
    measured_data = double(reshape(mua_allcat(:,:,curr_depth)',[],1));
    predicted_data = double(reshape(mua_ctxpred_allcat(:,:,curr_depth)',[],1));
    predicted_data(predicted_data < 0) = 0;

    n_bins = 100;
    activity_bounds = linspace(0,5,n_bins+1);
    activity_bin_centers = conv2(activity_bounds,[1,1]/2,'valid');

    measured_bins = discretize(measured_data,activity_bounds);
    predicted_bins = discretize(predicted_data,activity_bounds);

    measured_data_binmean = accumarray(predicted_bins(~isnan(predicted_bins)), ...
        measured_data(~isnan(predicted_bins)),[n_bins,1],@nanmedian,nan);
    predicted_data_binmean = accumarray(predicted_bins(~isnan(predicted_bins)),...
        predicted_data(~isnan(predicted_bins)),[n_bins,1],@nanmedian,nan);

    % smooth out the measured data binmean
    measured_data_binmean_smooth = medfilt1(measured_data_binmean,10);

    predicted_data_nlin = nan(size(predicted_data));
    predicted_data_nlin(~isnan(predicted_bins)) = measured_data_binmean_smooth(predicted_bins(~isnan(predicted_bins)));

    predicted_data_nlin_bins = discretize(predicted_data_nlin,activity_bounds);
    predicted_data_nlin_binmean = accumarray( ...
        predicted_bins(~isnan(predicted_bins)), ...
        predicted_data_nlin(~isnan(predicted_bins)),[n_bins,1],@mean,nan);

    mua_allcat_predicted_nlin(:,:,curr_depth) = ...
        reshape(predicted_data_nlin, ...
        size(mua_allcat,2),size(mua_allcat,1))';

    subplot(1,n_depths,curr_depth); hold on;
    plot(predicted_data,measured_data,'.')
    plot(predicted_data_binmean,measured_data_binmean,'linewidth',2);
    plot(predicted_data_binmean,measured_data_binmean_smooth,'linewidth',2);
    plot(predicted_data_nlin_binmean,measured_data_binmean,'linewidth',2);
    xlim([-2,6]);ylim([-2,6]);
    line(xlim,ylim,'color','k');
    xlabel('Predicted')
    ylabel('Measured')
    axis square;
end

%% Plot ctx L/ctx R/str L
% (raw/reduced)
% (by stim or grouped)

% Get flipped kernel ROIs and fluorescence
kernel_roiR = AP_reflect_widefield(kernel_roi);

fluor_kernelroiR_deconv = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],kernel_roiR), ...
    size(kernel_roiR,3),[],size(fluor_allcat_deconv,1)),[3,2,1]);

fluor_kernelroiR_deconv_std = ...
    cell2mat(cellfun(@(x) repmat(nanstd(reshape(x,[],1,size(kernel_roi,3)),[],1),size(x,1),1), ...
    mat2cell(fluor_kernelroiR_deconv,n_trials_day,length(t),size(kernel_roi,3)),'uni',false));

fluor_kernelroiR_deconv = fluor_kernelroiR_deconv./fluor_kernelroiR_deconv_std;

fluor_kernelroiR_taskpred_reduced = cell2mat(permute(arrayfun(@(x) ...
        permute(reshape( ...
        AP_svd_roi(U_master(:,:,1:n_vs), ...
        reshape(permute(fluor_taskpred_reduced_allcat(:,:,:,x),[3,2,1]), ...
        n_vs,[]),[],[],kernel_roiR), ...
        size(kernel_roi,3),[],size(fluor_taskpred_reduced_allcat,1)),[3,2,1]), ...
        1:size(fluor_taskpred_reduced_allcat,4),'uni',false),[1,3,4,2]))./fluor_kernelroiR_deconv_std;

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_contrast_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Pick alignments, time averages, conditions
use_align = stim_align;
use_t_avg = [0.05,0.15]; % stim = 50:150, move = -50:50, outcome = 50:150
% use_conditions =  ...
%     [sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1];
% use_conditions =  ...
%     [trial_contrastside_allcat <= -0.5,trial_contrastside_allcat >= 0.5];
use_condition_labels = {'Stim R','Stim L'};

use_conditions =  ...
    [trial_contrastside_allcat > 0,trial_contrastside_allcat == 0];
use_condition_labels = {'Stim R','No stim'};
use_reduction = 1;

use_t_avg_idx = t > use_t_avg(1) & t < use_t_avg(2);

% Re-align activity
curr_strL_act = cell2mat(arrayfun(@(trial) circshift( ...
    mua_allcat(trial,:,:) - mua_taskpred_reduced_allcat(trial,:,:,use_reduction), ...
    use_align(trial),2),transpose(1:length(use_align)),'uni',false));

curr_ctxL_act = cell2mat(arrayfun(@(trial) circshift( ...
    fluor_kernelroi_deconv(trial,:,:) - fluor_kernelroi_taskpred_reduced(trial,:,:,use_reduction), ...
    use_align(trial),2),transpose(1:length(use_align)),'uni',false));

curr_ctxR_act = cell2mat(arrayfun(@(trial) circshift( ...
    fluor_kernelroiR_deconv(trial,:,:) - fluor_kernelroiR_taskpred_reduced(trial,:,:,use_reduction), ...
    use_align(trial),2),transpose(1:length(use_align)),'uni',false));

% Time-average activity
curr_strL_act_avg = permute(nanmean(curr_strL_act(:,use_t_avg_idx,:),2),[1,3,2]);
curr_ctxL_act_avg = permute(nanmean(curr_ctxL_act(:,use_t_avg_idx,:),2),[1,3,2]);
curr_ctxR_act_avg = permute(nanmean(curr_ctxR_act(:,use_t_avg_idx,:),2),[1,3,2]);

% Plot
plot_trials = any(use_conditions,2);
plot_conditions = use_condition_labels(sum(use_conditions(plot_trials,:).*[1,2],2))';

diff(use_conditions(plot_trials,:),[],2);

curr_depth = 1;

figure('Name',['Str ' num2str(curr_depth)]); 

% (for all separately)
% curr_data = [curr_strL_act_avg(plot_trials,curr_depth), ...
%     curr_ctxL_act_avg(plot_trials,curr_depth),curr_ctxR_act_avg(plot_trials,curr_depth)];
% gplotmatrix(curr_data,[],plot_conditions,{'r','b'},'.',2,'on',[],{'Str_L','Ctx_L','Ctx_R'});
% (for L-R and L+R)
curr_data = [curr_strL_act_avg(plot_trials,curr_depth), ...
    curr_ctxL_act_avg(plot_trials,curr_depth) - curr_ctxR_act_avg(plot_trials,curr_depth), ...
    curr_ctxL_act_avg(plot_trials,curr_depth) + curr_ctxR_act_avg(plot_trials,curr_depth)];
gplotmatrix(curr_data,[],plot_conditions,{'r','b'},'.',2,'on',[],{'Str_L','Ctx_L-Ctx_R','Ctx_L+Ctx_R'});
legend(use_condition_labels);


%% Str L vs Ctx R

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_contrast_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set windows to average activity
timeavg_labels = {'Stim','Move onset','Outcome'};
timeavg_t = {[0.05,0.15],[-0.05,0.05],[0.05,0.15]};
timeavg_align = {stim_align,move_align,outcome_align};
timeavg_trial_conditions = ...
    {[sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1], ...
    [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
    [trial_outcome_allcat == 1, trial_outcome_allcat == -1]};
timeavg_task_reduction = [1,2,4];

% Set activity percentiles and bins
act_prctile = [10,90];
n_act_bins = 5;

% Set areas and conditions
plot_str_ctx = {[1,1],[2,2],[3,3],[4,4]};

% Loop across area pairs, plot binned predicted v measured activity
str_v_ctx_fig = figure('color','w');
for curr_str_ctx = 1:length(plot_str_ctx)
    
    plot_str = plot_str_ctx{curr_str_ctx}(1);
    plot_ctx = plot_str_ctx{curr_str_ctx}(2);

    for curr_mua = 1:2
        
        % (set striatum activity to use)
        switch curr_mua
            case 1
                curr_str_act_allcat = single(mua_allcat);
                curr_str_act_taskpred_reduced_allcat = single(mua_taskpred_reduced_allcat);
                mua_label = 'Measured';
            case 2
                curr_str_act_allcat = single(mua_ctxpred_allcat);
                curr_str_act_taskpred_reduced_allcat = single(mua_ctxpred_taskpred_reduced_allcat);
                mua_label = 'Ctx-pred';
        end
        
        for curr_timeavg = 1:length(timeavg_labels)
            
            trial_conditions = timeavg_trial_conditions{curr_timeavg};
            curr_task_reduction = timeavg_task_reduction(curr_timeavg);
            
            % (set cortex activity to use)
%             curr_ctx_act_allcat = fluor_roi_deconv;
            curr_ctx_act_allcat = fluor_kernelroiR_deconv;   
            curr_ctx_act_taskpred_reduced_allcat = fluor_kernelroiR_taskpred_reduced;
            
            % (re-align and split activity and conditions)
            curr_ctx_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_ctx_act_allcat(trial,:,:) - curr_ctx_act_taskpred_reduced_allcat(trial,:,:,curr_task_reduction), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_ctx_act_allcat,1)),'uni',false)), ...
                use_split,length(t),size(curr_ctx_act_allcat,3));
            
            curr_str_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_str_act_allcat(trial,:,:) - curr_str_act_taskpred_reduced_allcat(trial,:,:,curr_task_reduction), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_str_act_allcat,1)),'uni',false)), ...
                use_split,length(t),size(curr_str_act_allcat,3));
            
            trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
            
            % (get average activity within window)
            curr_event_t = t > timeavg_t{curr_timeavg}(1) & t < timeavg_t{curr_timeavg}(2);
            curr_ctx_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_ctx_act,'uni',false);
            curr_str_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_str_act,'uni',false);
            
            % (bin measured data across percentile range)
            bin_range = prctile(cell2mat(cellfun(@(x) x(:,plot_ctx),curr_ctx_act_avg,'uni',false)),act_prctile);
            bin_edges = linspace(bin_range(1),bin_range(2),n_act_bins+1);
            bin_centers = bin_edges(1:end-1) + diff(bin_edges)./2;
            
            trial_bins = cellfun(@(x) discretize(x,bin_edges),curr_ctx_act_avg,'uni',false);
            total_bins = cellfun(@(x) discretize(x,[bin_edges(1),bin_edges(end)]),curr_ctx_act_avg,'uni',false);
            
            % (get trials to use: no NaNs in either time series or bins)
            nonan_trials = cellfun(@(ctx_act,str_act,trial_bins) ...
                squeeze(~any(isnan(ctx_act(:,curr_event_t,plot_ctx)),2)) & ...
                squeeze(~any(isnan(str_act(:,curr_event_t,plot_str)),2)) & ...
                ~isnan(trial_bins), ...
                curr_ctx_act,curr_str_act,trial_bins,'uni',false);
            
            % (get average binned activity for measured/predicted by condition)
            ctx_act_binmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_ctx_act_avg,trial_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            str_act_binmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_str_act_avg,trial_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            % (get the average total activity)
            ctx_act_totalmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [1,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_ctx_act_avg,total_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            str_act_totalmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [1,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_str_act_avg,total_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            % Plot binned predicted v measured           
            subplot(length(plot_str_ctx),length(timeavg_labels), ...
                sub2ind(fliplr([length(plot_str_ctx),length(timeavg_labels)]),curr_timeavg,curr_str_ctx));
            hold on;
            col = lines(size(trial_conditions,2));
            switch curr_mua
                case 1
                    errorbar( ...
                        squeeze(nanmean(ctx_act_binmean(:,plot_ctx,:,:),3)), ...
                        squeeze(nanmean(str_act_binmean(:,plot_str,:,:),3)), ...
                        squeeze(AP_sem(str_act_binmean(:,plot_str,:,:),3)), ...
                        'linewidth',2,'CapSize',0);
%                     errorbar( ...
%                         squeeze(nanmean(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
%                         squeeze(nanmean(str_act_totalmean(:,plot_str,:,:),3)), ...
%                         squeeze(AP_sem(str_act_totalmean(:,plot_str,:,:),3)), ...
%                         squeeze(AP_sem(str_act_totalmean(:,plot_str,:,:),3)), ...
%                         squeeze(AP_sem(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
%                         squeeze(AP_sem(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
%                         '.','linewidth',3,'CapSize',0);
                case 2
                    for curr_cond = 1:size(trial_conditions,2)
                        AP_errorfill( ...
                            squeeze(nanmean(ctx_act_binmean(:,plot_ctx,:,curr_cond),3)), ...
                            squeeze(nanmean(str_act_binmean(:,plot_str,:,curr_cond),3)), ...
                            squeeze(AP_sem(str_act_binmean(:,plot_str,:,curr_cond),3)), ...
                            col(curr_cond,:),0.5,false);
                    end
            end
            xlabel(['Ctx (' num2str(plot_ctx) ')']);
            ylabel([' Str (' num2str(plot_str) ')'])
            title([timeavg_labels{curr_timeavg} '(' task_regressor_labels{curr_task_reduction} '-reduced)']);
            
        end
        
    end
end

% Link axes of all plots
linkaxes(get(str_v_ctx_fig,'Children'));

%% Str L vs Ctx R (binned by Ctx L)

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_contrast_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set windows to average activity
% timeavg_labels = {'Stim','Move onset','Outcome'};
% timeavg_t = {[0.05,0.15],[-0.05,0.05],[0.05,0.15]};
% timeavg_align = {stim_align,move_align,outcome_align};
% timeavg_trial_conditions = ...
%     {[sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1], ...
%     [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
%     [trial_outcome_allcat == 1, trial_outcome_allcat == -1]};
% timeavg_task_reduction = [1,2,4];

timeavg_labels = {'Stim'};
timeavg_t = {[0.05,0.15]};
timeavg_align = {stim_align};
timeavg_trial_conditions = ...
    {[trial_contrastside_allcat > 0,trial_contrastside_allcat == 0]};
timeavg_task_reduction = [1];

% Set activity percentiles and bins
act_prctile = [10,90];
n_act_bins = 5;

% Set areas and conditions
plot_str_ctx = {[1,1]};

% Loop across area pairs, plot binned predicted v measured activity
str_v_ctx_fig = figure('Name', ...
    [timeavg_labels{curr_timeavg} '(' task_regressor_labels{curr_task_reduction} '-reduced)'], ...
    'color','w');
for curr_str_ctx = 1:length(plot_str_ctx)
    
    plot_str = plot_str_ctx{curr_str_ctx}(1);
    plot_ctx = plot_str_ctx{curr_str_ctx}(2);

    for curr_mua = 1:2
        
        % (set striatum activity to use)
        switch curr_mua
            case 1
                curr_str_act_allcat = single(mua_allcat);
                curr_str_act_taskpred_reduced_allcat = single(mua_taskpred_reduced_allcat);
                mua_label = 'Measured';
            case 2
                curr_str_act_allcat = single(mua_ctxpred_allcat);
                curr_str_act_taskpred_reduced_allcat = single(mua_ctxpred_taskpred_reduced_allcat);
                mua_label = 'Ctx-pred';
        end
        
        for curr_timeavg = 1:length(timeavg_labels)
            
            trial_conditions = timeavg_trial_conditions{curr_timeavg};
            curr_task_reduction = timeavg_task_reduction(curr_timeavg);
            
            % (set cortex activity to use)
%             curr_ctx_act_allcat = fluor_roi_deconv;
            curr_ctx_act_allcat = fluor_kernelroi_deconv;   
            curr_ctx_act_taskpred_reduced_allcat = fluor_kernelroi_taskpred_reduced;
            
            curr_ctx_contra_act_allcat = fluor_kernelroiR_deconv;
            curr_ctx_contra_act_taskpred_reduced_allcat = fluor_kernelroiR_taskpred_reduced;
            
            % (re-align and split activity and conditions)
            curr_ctx_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_ctx_act_allcat(trial,:,:) - curr_ctx_act_taskpred_reduced_allcat(trial,:,:,curr_task_reduction), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_ctx_act_allcat,1)),'uni',false)), ...
                use_split,length(t),size(curr_ctx_act_allcat,3));
            
            curr_ctxR_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_ctx_contra_act_allcat(trial,:,:) - curr_ctx_contra_act_taskpred_reduced_allcat(trial,:,:,curr_task_reduction), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_ctx_contra_act_allcat,1)),'uni',false)), ...
                use_split,length(t),size(curr_ctx_contra_act_allcat,3));
            
            curr_str_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_str_act_allcat(trial,:,:) - curr_str_act_taskpred_reduced_allcat(trial,:,:,curr_task_reduction), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_str_act_allcat,1)),'uni',false)), ...
                use_split,length(t),size(curr_str_act_allcat,3));
            
            trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
            
            % (get average activity within window)
            curr_event_t = t > timeavg_t{curr_timeavg}(1) & t < timeavg_t{curr_timeavg}(2);
            curr_ctx_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_ctx_act,'uni',false);
            curr_ctxR_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_ctxR_act,'uni',false);
            curr_str_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_str_act,'uni',false);
            
            % (bin measured data across percentile range)
            bin_range = prctile(cell2mat(cellfun(@(x) x(:,plot_ctx),curr_ctx_act_avg,'uni',false)),act_prctile);
            bin_edges = linspace(bin_range(1),bin_range(2),n_act_bins+1);
            bin_centers = bin_edges(1:end-1) + diff(bin_edges)./2;
            
            trial_bins = cellfun(@(x) discretize(x,bin_edges),curr_ctx_act_avg,'uni',false);
            total_bins = cellfun(@(x) discretize(x,[bin_edges(1),bin_edges(end)]),curr_ctx_act_avg,'uni',false);
            
            % (bin Ctx_R)
            trial_bins_R = cellfun(@(x) discretize(x,bin_edges),curr_ctxR_act_avg,'uni',false);
            
            % (get trials to use: no NaNs in either time series or bins)
            nonan_trials = cellfun(@(ctx_act,str_act,trial_bins) ...
                squeeze(~any(isnan(ctx_act(:,curr_event_t,plot_ctx)),2)) & ...
                squeeze(~any(isnan(str_act(:,curr_event_t,plot_str)),2)) & ...
                ~isnan(trial_bins), ...
                curr_ctx_act,curr_str_act,trial_bins,'uni',false);      
            
            % Plot binned predicted v measured
            % (for each bin of Ctx L)
            for curr_bin = 1:n_act_bins                          
                
                % Get Str and Ctx_R activity within this Ctx_L bin
                % (get average binned activity for measured/predicted by condition)
                curr_use_trials = cellfun(@(nonan_trials,binL_trials,binR_trials) ...
                    nonan_trials & binL_trials == curr_bin & ~isnan(binR_trials), ...
                    nonan_trials,trial_bins,trial_bins_R,'uni',false);
                
                ctxR_act_binmean = cell2mat(arrayfun(@(condition) ...
                    cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                    accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                    act(use_trials(:,area) & trial_cond(:,condition),area), ...
                    [n_act_bins,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                    curr_ctxR_act_avg,trial_bins_R,trial_conditions_exp,curr_use_trials,'uni',false),[2,3,1])), ...
                    permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
                
                str_act_binmean = cell2mat(arrayfun(@(condition) ...
                    cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                    accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                    act(use_trials(:,area) & trial_cond(:,condition),area), ...
                    [n_act_bins,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                    curr_str_act_avg,trial_bins_R,trial_conditions_exp,curr_use_trials,'uni',false),[2,3,1])), ...
                    permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));           
                
                subplot(length(plot_str_ctx),n_act_bins, ...
                    sub2ind([n_act_bins,length(plot_str_ctx)],curr_bin,curr_str_ctx));
                hold on;
                col = lines(size(trial_conditions,2));
                switch curr_mua
                    case 1
                        % (as in figure)
%                         errorbar( ...
%                             squeeze(nanmean(ctxR_act_binmean(:,plot_ctx,:,:),3)), ...
%                             squeeze(nanmean(str_act_binmean(:,plot_str,:,:),3)), ...
%                             squeeze(AP_sem(str_act_binmean(:,plot_str,:,:),3)), ...
%                             'linewidth',2,'CapSize',0);
                        
                        % (to just make scatter plot)
                        for curr_condition = 1:2
                            curr_ctxR_scatter = cell2mat(cellfun(@(act,trials,cond) ...
                                act(trials(:,plot_ctx) & cond(:,curr_condition),plot_ctx), ...
                                curr_ctxR_act_avg,curr_use_trials,trial_conditions_exp,'uni',false));
                            curr_str_scatter = cell2mat(cellfun(@(act,trials,cond) ...
                                act(trials(:,plot_str) & cond(:,curr_condition),plot_str), ...
                                curr_str_act_avg,curr_use_trials,trial_conditions_exp,'uni',false));
                            plot(curr_ctxR_scatter,curr_str_scatter,'.','MarkerSize',5);
                        end
                    
       
                    case 2
                        % (don't plot ctx>str now)
%                         for curr_cond = 1:size(trial_conditions,2)
%                             AP_errorfill( ...
%                                 squeeze(nanmean(ctxR_act_binmean(:,plot_ctx,:,curr_cond),3)), ...
%                                 squeeze(nanmean(str_act_binmean(:,plot_str,:,curr_cond),3)), ...
%                                 squeeze(AP_sem(str_act_binmean(:,plot_str,:,curr_cond),3)), ...
%                                 col(curr_cond,:),0.5,false);
%                         end
                end
                xlabel(['Ctx_R (' num2str(plot_ctx) ')']);
                ylabel([' Str (' num2str(plot_str) ')'])
                title(sprintf('Ctx_L bin %0.2f:%0.2f',bin_edges(curr_bin),bin_edges(curr_bin+1)));
                
            end
            
        end
        
    end
end

% Link axes of all plots
linkaxes(get(str_v_ctx_fig,'Children'));

%% Ctx -> Str regression concat and timepoint-specific 
% (this is somewhere in old scripts but easier to redo)

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [0,0];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
    round(regression_params.kernel_t(2)*sample_rate);

lambda = 10;

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_contrast_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set trials/depth/alignment to use
use_trials = move_t < 0.5; % & trial_contrastside_allcat > 0;
curr_depth = 1;
use_align = stim_align;

% Re-align activity
curr_str_act = cell2mat(arrayfun(@(trial) circshift( ...
    mua_allcat(trial,:,:), ...
    use_align(trial),2),transpose(1:length(use_align)),'uni',false));

curr_ctx_act = cell2mat(arrayfun(@(trial) circshift( ...
    fluor_allcat_deconv(trial,:,:), ...
    use_align(trial),2),transpose(1:length(use_align)),'uni',false));

% Regress cortex to striatum (all time together)
str_data = reshape(permute(curr_str_act(use_trials,:,:),[2,1,3]),[],n_depths)';
ctx_data = reshape(permute(curr_ctx_act(use_trials,:,:),[2,1,3]),[],n_vs)';

[ctx_str_k,ctxpred_spikes_std,explained_var] = ...
    AP_regresskernel(ctx_data(regression_params.use_svs,:), ...
    str_data(curr_depth,:),kernel_frames,lambda, ...
    regression_params.zs,regression_params.cvfold, ...
    true,regression_params.use_constant);

k = svdFrameReconstruct(U_master(:,:,regression_params.use_svs),ctx_str_k{1});

% Regress cortex to striatum (timepoint-by-timepoint)
k_t = cell(2,1);
k_t{1} = nan(size(U_master,1),size(U_master,2),length(t));
k_t{2} = nan(length(t),1);
expl_var_t = nan(length(t),1);
for curr_t = 1:length(t)
    
    str_data = reshape(permute(curr_str_act(use_trials,curr_t,:),[2,1,3]),[],n_depths)';
    ctx_data = reshape(permute(curr_ctx_act(use_trials,curr_t,:),[2,1,3]),[],n_vs)';
        
    [ctx_str_k,ctxpred_spikes_std,explained_var_t] = ...
        AP_regresskernel(ctx_data(regression_params.use_svs,:), ...
        str_data(curr_depth,:),kernel_frames,lambda, ...
        regression_params.zs,regression_params.cvfold, ...
        true,regression_params.use_constant);
    
    k_t{1}(:,:,curr_t) = svdFrameReconstruct(U_master(:,:,regression_params.use_svs),ctx_str_k{1});
    k_t{2}(curr_t) = ctx_str_k{2};
    expl_var_t(curr_t) = explained_var_t.total;
    
    AP_print_progress_fraction(curr_t,length(t))
    
end

% Plot kernels
AP_imscroll(k,{'Concat time'})
caxis([-max(caxis),max(caxis)]);
colormap(brewermap([],'*RdBu'));
axis image
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

AP_imscroll(k_t{1},t)
caxis([-max(caxis),max(caxis)]);
colormap(brewermap([],'*RdBu'));
axis image
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

figure;plot(t,k_t{2},'k','linewidth',2)
ylabel('Offset');
xlabel('Time');
line(xlim,[0,0],'color','k');
line([0,0],ylim,'color','k');

% Plot explained variance; 
figure; hold on
plot(t,expl_var_t,'k','linewidth',2);
line(xlim,repmat(explained_var.total,2,1),'color','r','linewidth',2);
ylabel('Explained variance');
xlabel('Time');
legend({'Concat time','Timepoint'})



%% Ctx L only to Str (hacky)

% Load bregma and master CCF tform
bregma = allenCCFbregma;
bregma(3) = bregma(3) + 0.5;
ccf_tform_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\ccf_tform'];
load(ccf_tform_fn);

um2pixel = 20.6;
bregma_resize = bregma*(10/um2pixel);
bregma_align = [bregma_resize([3,1]),1]*ccf_tform.T;

U_master_L = U_master;
U_master_L(:,round(bregma_align(1)):end,:) = 0;

U_tform_matrix = reshape(U_master(:,:,1:n_vs),[],n_vs)'* ...
    reshape(U_master_L(:,:,1:n_vs),[],n_vs);
fluor_allcat_deconv_L = permute(reshape(transpose( ...
    U_tform_matrix*reshape(permute(fluor_allcat_deconv,[2,1,3]),[],n_vs)'), ...
    size(permute(fluor_allcat_deconv,[2,1,3]))),[2,1,3]);

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [0,0];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
    round(regression_params.kernel_t(2)*sample_rate);

lambda = 10;

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_contrast_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set trials/depth/alignment to use
use_trials = move_t < 0.5; % & trial_contrastside_allcat > 0;
curr_depth = 1;
use_align = stim_align;

% Re-align activity
curr_str_act = cell2mat(arrayfun(@(trial) circshift( ...
    mua_allcat(trial,:,:), ...
    use_align(trial),2),transpose(1:length(use_align)),'uni',false));

curr_ctx_act = cell2mat(arrayfun(@(trial) circshift( ...
    fluor_allcat_deconv_L(trial,:,:), ...
    use_align(trial),2),transpose(1:length(use_align)),'uni',false));

% Regress cortex to striatum (all time together)
str_data = reshape(permute(curr_str_act(use_trials,:,:),[2,1,3]),[],n_depths)';
ctx_data = reshape(permute(curr_ctx_act(use_trials,:,:),[2,1,3]),[],n_vs)';

[ctx_str_k,ctxpred_spikes_std,explained_var] = ...
    AP_regresskernel(ctx_data(regression_params.use_svs,:), ...
    str_data(curr_depth,:),kernel_frames,lambda, ...
    regression_params.zs,regression_params.cvfold, ...
    true,regression_params.use_constant);

k = svdFrameReconstruct(U_master_L(:,:,regression_params.use_svs),ctx_str_k{1});

% Regress cortex to striatum (timepoint-by-timepoint)
k_t = cell(2,1);
k_t{1} = nan(size(U_master,1),size(U_master,2),length(t));
k_t{2} = nan(length(t),1);
expl_var_t = nan(length(t),1);
for curr_t = 1:length(t)
    
    str_data = reshape(permute(curr_str_act(use_trials,curr_t,:),[2,1,3]),[],n_depths)';
    ctx_data = reshape(permute(curr_ctx_act(use_trials,curr_t,:),[2,1,3]),[],n_vs)';
        
    [ctx_str_k,ctxpred_spikes_std,explained_var_t] = ...
        AP_regresskernel(ctx_data(regression_params.use_svs,:), ...
        str_data(curr_depth,:),kernel_frames,lambda, ...
        regression_params.zs,regression_params.cvfold, ...
        true,regression_params.use_constant);
    
    k_t{1}(:,:,curr_t) = svdFrameReconstruct(U_master_L(:,:,regression_params.use_svs),ctx_str_k{1});
    k_t{2}(curr_t) = ctx_str_k{2};
    expl_var_t(curr_t) = explained_var_t.total;
    
    AP_print_progress_fraction(curr_t,length(t))
    
end

% Plot kernels
AP_imscroll(k,{'Concat time'})
caxis([-max(caxis),max(caxis)]);
colormap(brewermap([],'*RdBu'));
axis image
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

AP_imscroll(k_t{1},t)
caxis([-max(caxis),max(caxis)]);
colormap(brewermap([],'*RdBu'));
axis image
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

figure;plot(t,k_t{2},'k','linewidth',2)
ylabel('Offset');
xlabel('Time');
line(xlim,[0,0],'color','k');
line([0,0],ylim,'color','k');

% Plot explained variance; 
figure; hold on
plot(t,expl_var_t,'k','linewidth',2);
line(xlim,repmat(explained_var.total,2,1),'color','r','linewidth',2);
ylabel('Explained variance');
xlabel('Time');
legend({'Concat time','Timepoint'})


%% Ctx -> Str regression: constant vs time-varying ctx weights (1 L/R ROIs)

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_contrast_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

use_align = stim_align;
use_trials = move_t < 0.5;
plot_str = 1; % 1,2
plot_ctx_L = 3; % 3,7
plot_ctx_R = plot_ctx_L + size(wf_roi,1);

% Re-align activity
curr_str = cell2mat(arrayfun(@(trial) circshift( ...
    mua_allcat(trial,:,plot_str), ...
    use_align(trial),2),find(use_trials),'uni',false));

curr_ctx_L = cell2mat(arrayfun(@(trial) circshift( ...
    fluor_roi_deconv(trial,:,plot_ctx_L), ...
    use_align(trial),2),find(use_trials),'uni',false));

curr_ctx_R = cell2mat(arrayfun(@(trial) circshift( ...
    fluor_roi_deconv(trial,:,plot_ctx_R), ...
    use_align(trial),2),find(use_trials),'uni',false));


% % Re-align activity (reduced)
% use_reduction = 1;
% 
% curr_str = cell2mat(arrayfun(@(trial) circshift( ...
%     mua_allcat(trial,:,plot_str) - mua_taskpred_reduced_allcat(trial,:,plot_str,use_reduction), ...
%     use_align(trial),2),find(use_trials),'uni',false));
% 
% curr_ctx_L = cell2mat(arrayfun(@(trial) circshift( ...
%     fluor_roi_deconv(trial,:,plot_ctx_L) - fluor_roi_taskpred_reduced(trial,:,plot_ctx_L,use_reduction), ...
%     use_align(trial),2),find(use_trials),'uni',false));
% 
% curr_ctx_R = cell2mat(arrayfun(@(trial) circshift( ...
%     fluor_roi_deconv(trial,:,plot_ctx_R) - fluor_roi_taskpred_reduced(trial,:,plot_ctx_R,use_reduction), ...
%     use_align(trial),2),find(use_trials),'uni',false));

% Regress Str from Ctx L,R and get cross-validated expl var
k = nan(length(t),3);
ev = nan(length(t),2);
predicted_str_t = nan(sum(use_trials),length(t));
for curr_t = 1:length(t)
        
    [kLR,predicted_signals,evLR] = ...
        AP_regresskernel([curr_ctx_L(:,curr_t),curr_ctx_R(:,curr_t)]', ...
        curr_str(:,curr_t)',0,0,[false,false],10,1,1,[]);   
    
    predicted_str_t(:,curr_t) = predicted_signals;
    k(curr_t,:) = cell2mat(kLR)';
    ev(curr_t,1) = evLR.total;
    
end

figure; 
plot(t,k,'linewidth',2);
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');
legend({'Ctx_L','Ctx_R','Offset'});
xlabel('Time');
ylabel('Weight');
title('Str regression');


% Kenneth asked: fit a*Ctx_L + b*Ctx_R + Xc in stages
% (get static a and b first, estimate Xc from residual, refit a/b)
% (compare this to a changing a and b)
% (wants this towards: does changing weights buy anything?)

% Regress Str from Ctx L,R, estimate Xc from residual
[k,predicted_str_noxc,expl_var] = ...
    AP_regresskernel([reshape(curr_ctx_L',[],1),reshape(curr_ctx_R',[],1)]', ...
    reshape(curr_str',[],1)',0,0,[false,false],10,1,1,[]);

predicted_str_allt_noxc = reshape(predicted_str_noxc,length(t),[])';

xc = grpstats(curr_str - reshape(predicted_str_noxc,[],sum(use_trials))', ...
    trial_contrastside_allcat(use_trials),@nanmean);

% Regress (Str - Xc) from Ctx L,R
[~,contrastside_idx] = ismember(trial_contrastside_allcat(use_trials), ...
    unique(trial_contrastside_allcat),'rows');
xc_trial = xc(contrastside_idx,:);

curr_str_xc = curr_str - xc_trial;

[k_xc,predicted_str_minusxc,expl_var_xc] = ...
    AP_regresskernel([reshape(curr_ctx_L',[],1),reshape(curr_ctx_R',[],1)]', ...
    reshape(curr_str_xc',[],1)',0,0,[false,false],10,1,1,[]);

predicted_str_xc = predicted_str_minusxc + reshape(xc_trial',[],1)';

predicted_str_allt = reshape(predicted_str_xc,length(t),[])';

% Get R^2 for both predictions
sse_total = nansum((reshape(curr_str',[],1) - nanmean(reshape(curr_str',[],1))).^2);
sse_residual_t = nansum((reshape(curr_str',[],1) - reshape(predicted_str_t',[],1)).^2);
sse_residual_allt_noxc = nansum((reshape(curr_str',[],1) - reshape(predicted_str_allt_noxc',[],1)).^2);
sse_residual_allt = nansum((reshape(curr_str',[],1) - reshape(predicted_str_allt',[],1)).^2);

expl_var_t = 1 - (sse_residual_t/sse_total);
expl_var_allt_noxc = 1 - (sse_residual_allt_noxc/sse_total);
expl_var_allt = 1 - (sse_residual_allt/sse_total);

figure; 
plot([expl_var_t,expl_var_allt_noxc,expl_var_allt],'k','linewidth',2);
set(gca,'XTick',1:3,'XTickLabel',{'Timepoint','All time (no X_c)','All time (with X_c)'});
ylabel('Explained variance');
xlim([0.5,3.5])

% Plot Xc
col = colormap_BlueWhiteRed(length(unique(abs(trial_contrastside_allcat(trial_contrastside_allcat ~= 0)))));
if ~any(trial_contrastside_allcat == 0)
    col(ceil(size(col,1)/2),:) = [];
end

figure; hold on
set(gca,'ColorOrder',col);
plot(t,xc','linewidth',2);
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');



%% Ctx -> Str regression: constant vs time-varying ctx weights (Vs)

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_contrast_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

use_align = stim_align;
use_trials = move_t < 0.5;
plot_str = 3; % 1,2

use_vs = 1:50;
lambda = 10;

% Re-align activity
curr_str = cell2mat(arrayfun(@(trial) circshift( ...
    mua_allcat(trial,:,plot_str), ...
    use_align(trial),2),find(use_trials),'uni',false));

curr_ctx = cell2mat(arrayfun(@(trial) circshift( ...
    fluor_allcat_deconv(trial,:,use_vs), ...
    use_align(trial),2),find(use_trials),'uni',false));

% Regress Str from Ctx and get cross-validated expl var
k = nan(length(use_vs),length(t));
ev = nan(length(t),2);
predicted_str_t = nan(sum(use_trials),length(t));
for curr_t = 1:length(t)
        
    [kLR,predicted_signals,evLR] = ...
        AP_regresskernel(permute(curr_ctx(:,curr_t,:),[3,1,2]), ...
        curr_str(:,curr_t)',0,lambda,[false,false],10,1,1,[]);   
    
    predicted_str_t(:,curr_t) = predicted_signals;
    k(:,curr_t) = kLR{1};
    
end

figure; 
AP_imscroll(svdFrameReconstruct(U_master(:,:,use_vs),k));
caxis([-max(caxis),max(caxis)]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
axis image off


% Kenneth asked: fit a*Ctx_L + b*Ctx_R + Xc in stages
% (get static a and b first, estimate Xc from residual, refit a/b)
% (compare this to a changing a and b)
% (wants this towards: does changing weights buy anything?)

% Regress Str from Ctx L,R, estimate Xc from residual
[k,predicted_str_noxc,expl_var] = ...
    AP_regresskernel(reshape(permute(curr_ctx,[2,1,3]),[],length(use_vs))', ...
    reshape(curr_str',[],1)',0,lambda,[false,false],10,1,1,[]);

predicted_str_allt_noxc = reshape(predicted_str_noxc,length(t),[])';

xc = grpstats(curr_str - reshape(predicted_str_noxc,[],sum(use_trials))', ...
    trial_contrastside_allcat(use_trials),@nanmean);

% Regress (Str - Xc) from Ctx L,R
[~,contrastside_idx] = ismember(trial_contrastside_allcat(use_trials), ...
    unique(trial_contrastside_allcat),'rows');
xc_trial = xc(contrastside_idx,:);

curr_str_xc = curr_str - xc_trial;

[k_xc,predicted_str_minusxc,expl_var_xc] = ...
    AP_regresskernel(reshape(permute(curr_ctx,[2,1,3]),[],length(use_vs))', ...
    reshape(curr_str_xc',[],1)',0,0,[false,false],10,1,1,[]);

predicted_str_xc = predicted_str_minusxc + reshape(xc_trial',[],1)';

predicted_str_allt = reshape(predicted_str_xc,length(t),[])';

figure; 
imagesc(svdFrameReconstruct(U_master(:,:,use_vs),k_xc{1}));
caxis([-max(caxis),max(caxis)]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
axis image off


% Get R^2 for both predictions
sse_total = nansum((reshape(curr_str',[],1) - nanmean(reshape(curr_str',[],1))).^2);
sse_residual_t = nansum((reshape(curr_str',[],1) - reshape(predicted_str_t',[],1)).^2);
sse_residual_allt_noxc = nansum((reshape(curr_str',[],1) - reshape(predicted_str_allt_noxc',[],1)).^2);
sse_residual_allt = nansum((reshape(curr_str',[],1) - reshape(predicted_str_allt',[],1)).^2);

expl_var_t = 1 - (sse_residual_t/sse_total);
expl_var_allt_noxc = 1 - (sse_residual_allt_noxc/sse_total);
expl_var_allt = 1 - (sse_residual_allt/sse_total);

figure; 
plot([expl_var_t,expl_var_allt_noxc,expl_var_allt],'k','linewidth',2);
set(gca,'XTick',1:3,'XTickLabel',{'Timepoint','All time (no X_c)','All time (with X_c)'});
ylabel('Explained variance');
xlim([0.5,3.5])

% Plot Xc
col = colormap_BlueWhiteRed(length(unique(abs(trial_contrastside_allcat(trial_contrastside_allcat ~= 0)))));
if ~any(trial_contrastside_allcat == 0)
    col(ceil(size(col,1)/2),:) = [];
end

figure; hold on
set(gca,'ColorOrder',col);
plot(t,xc','linewidth',2);
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');
xlabel('Time');
ylabel('X_c');

%% Ctx -> Str regression: iterating Xc 

use_trials = move_t < 0.5;

curr_str = mua_allcat(use_trials,:,1);
curr_ctx_L = fluor_roi_deconv(use_trials,:,3);
curr_ctx_R = fluor_roi_deconv(use_trials,:,13);

% curr_ctx_L = fluor_kernelroi_deconv(use_trials,:,1);
% curr_ctx_R = fluor_kernelroiR_deconv(use_trials,:,1);

% Initialize Xc
xc = zeros(length(unique(trial_contrastside_allcat)),length(t));

% Regress (Str - Xc) from Ctx L,R
n_iter = 10;
expl_var_allt = nan(n_iter,1);
for curr_iter = 1:n_iter
    figure;imagesc(xc);
    
    [~,contrastside_idx] = ismember(trial_contrastside_allcat(use_trials), ...
        unique(trial_contrastside_allcat),'rows');
    xc_trial = xc(contrastside_idx,:);
    
    curr_str_xc = curr_str - xc_trial;
    
    [k_xc,predicted_str_minusxc,expl_var_xc] = ...
        AP_regresskernel([reshape(curr_ctx_L',[],1),reshape(curr_ctx_R',[],1)]', ...
        reshape(curr_str_xc',[],1)',0,0,[false,false],10,1,1,[]);
    
    predicted_str_xc = predicted_str_minusxc + reshape(xc_trial',[],1)';
    
    predicted_str_allt = reshape(predicted_str_xc,length(t),[])';
    
    xc = xc + grpstats(curr_str - reshape(predicted_str_xc,[],sum(use_trials))', ...
        trial_contrastside_allcat(use_trials),@nanmean);
    
    % Get R^2
    sse_total = nansum((reshape(curr_str',[],1) - nanmean(reshape(curr_str',[],1))).^2);
    sse_residual_allt = nansum((reshape(curr_str',[],1) - reshape(predicted_str_allt',[],1)).^2);
    
    expl_var_allt(curr_iter) = 1 - (sse_residual_allt/sse_total);
    
end

figure; 
plot(expl_var_allt,'k','linewidth',2);
xlabel('Xc iteration');
ylabel('Explained variance');


%% Ctx -> Str: w*ctx + Xc, w from pre-stim or stim

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_contrast_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

%%%%%% TESTING NEW ALIGNMENT: 
% max stim response for each contrast (better than stim onset)
str1_stimresponses = grpstats(mua_allcat(:,:,1),trial_contrastside_allcat,@nanmean);
[~,stimresponse_max_idx] = max(str1_stimresponses(:,t < 0.5),[],2);
stimresponse_max_idx(1:6) = stimresponse_max_idx(7);
[~,stim_idx] = ismember(trial_contrastside_allcat,unique(trial_contrastside_allcat),'rows');
stim_max_align = -stimresponse_max_idx(stim_idx) + leeway_samples;
%%%%%%

use_align = stim_align;
plot_str = 1; % 1,2
plot_ctx_L = 3; % 3,7
plot_ctx_R = plot_ctx_L + size(wf_roi,1);

use_trials = move_t < 0.5 & ...
    ~any(isnan(mua_allcat(:,:,plot_str)),2) & ...
    ~any(isnan(fluor_roi_deconv(:,:,plot_ctx_L)),2);

% Re-align activity
curr_str = cell2mat(arrayfun(@(trial) circshift( ...
    mua_allcat(trial,:,plot_str), ...
    use_align(trial),2),find(use_trials),'uni',false));

curr_ctx_L = cell2mat(arrayfun(@(trial) circshift( ...
    fluor_roi_deconv(trial,:,plot_ctx_L), ...
    use_align(trial),2),find(use_trials),'uni',false));

% % Re-align activity (reduced)
% use_reduction = 1;
% 
% curr_str = cell2mat(arrayfun(@(trial) circshift( ...
%     mua_allcat(trial,:,plot_str) - mua_taskpred_reduced_allcat(trial,:,plot_str,use_reduction), ...
%     use_align(trial),2),find(use_trials),'uni',false));
% 
% curr_ctx_L = cell2mat(arrayfun(@(trial) circshift( ...
%     fluor_roi_deconv(trial,:,plot_ctx_L) - fluor_roi_taskpred_reduced(trial,:,plot_ctx_L,use_reduction), ...
%     use_align(trial),2),find(use_trials),'uni',false));


% % Regress Str = W*Ctx + Xc
% % NOTE THIS CROSS-VAL IS JANKY AT THE MOMENT:
% % it's using the average W across cv (should use proper train/test
% % combinations for each cv's W)
% 
% t_shifts = 0;
% lambdas = 0;
% zs = [false,false];
% cvfold = 10;
% return_constant = true;
% use_constant = true; 
% 
% expl_var = nan(length(t));
% expl_var_xc = nan(length(t));
% for curr_t = 1:length(t)
%     
%     % Get K from one time point
%     k_ctx = ...
%         AP_regresskernel(curr_ctx_L(:,curr_t)',curr_str(:,curr_t)', ...
%         t_shifts,lambdas,zs,cvfold,return_constant,use_constant);
%     
%     % Apply K to all time points
%     ctx_str = reshape(cell2mat(k_ctx)'* ...
%         [reshape(curr_ctx_L',[],1),ones(numel(curr_ctx_L),1)]',length(t),[])';
%     
%     % Estimate offset by contrast
%     xc = grpstats(curr_str - ctx_str, ...
%         trial_contrastside_allcat(use_trials),@nanmean);
%     
%     [~,contrastside_idx] = ismember(trial_contrastside_allcat(use_trials), ...
%         unique(trial_contrastside_allcat),'rows');
%     xc_trial = xc(contrastside_idx,:);
% 
%     % Get R^2 for each time point
%     sse_total = nansum((curr_str - nanmean(curr_str,1)).^2);
%     sse_residual = nansum((curr_str - ctx_str).^2);
%     sse_residual_xc = nansum((curr_str - (ctx_str + xc_trial)).^2);
%     
%     expl_var(curr_t,:) = 1 - (sse_residual./sse_total);
%     expl_var_xc(curr_t,:) = 1 - (sse_residual_xc./sse_total);
%     
%     AP_print_progress_fraction(curr_t,length(t));
%     
% end
% 
% figure;
% 
% subplot(6,1,1:4);
% imagesc(t,t,expl_var_xc - expl_var);
% caxis([-max(caxis),max(caxis)]);
% line([0,0],ylim,'color','k');
% line(xlim,[0,0],'color','k');
% ylabel('Trained weights');
% xlabel('Tested fit');
% title('\DeltaExpl var: with X_c - without X_c');
% c = colorbar; ylabel(c,'Fraction expl var');
% colormap(brewermap([],'*RdBu'));
% axis square
% 
% subplot(6,1,5);
% plot(t,nanmean(expl_var_xc - expl_var,1),'k','linewidth',2);
% xlabel('Time');
% ylabel('\DeltaExpl var')
% 
% stim_t = 21;
% line(repmat(t(stim_t),2,1),ylim,'color','r');
% subplot(6,1,6);
% plot(t,expl_var_xc(:,stim_t) - expl_var(:,stim_t),'k','linewidth',2);
% xlabel('Trained weights');
% ylabel('\DeltaExpl var');
% ylim([0,max(ylim)]);


%%%%%%%% WAS WORKING HERE: 
% doing non-cv should have a hot diagonal (weights from t are best for t)
% but it didn't, looked like something was wrong
% (and anyway all this is pointless if prediction error can't be fixed)

% Alternate Method: 
% str = W*fluor + Xc + mu for a time point 
% (regressor matrix is fluor and trial x contrast X matrix)
% apply W and mu to all other time points, estimate Xc from residual 

cvfold = 10;

% Set up Xr regressors (trial x contrast)
stim_used = unique(curr_contrastside);
Xr = zeros(sum(use_trials),length(stim_used));
for curr_stim_idx = 1:length(stim_used)
    Xr(curr_contrastside == ...
        stim_used(curr_stim_idx),curr_stim_idx) = 1;
end

curr_contrastside = trial_contrastside_allcat(use_trials);
[~,contrastside_idx] = ismember(curr_contrastside,stim_used,'rows');

expl_var_cv = nan(1,length(t));
expl_var_xc = nan(length(t));
for curr_t = 1:length(t)
    
    % (doing regression manually here to correctly cross-weight cv)   
    cv_set = AP_shake(round(linspace(1,cvfold,sum(use_trials)))');
    
    ctx_str_cv = nan(size(curr_str,1),1);
    ctx_str_xc = nan(size(curr_str));
    for curr_cv = 1:cvfold
        
        train_trials = cv_set ~= curr_cv;
        test_trials = cv_set == curr_cv;
        
        % (cancel CV)
        train_trials = true(size(cv_set));
        test_trials = true(size(cv_set));
        
        % Do standard cross-validated regression 
        k_cv = [curr_ctx_L(train_trials,curr_t),ones(sum(train_trials),1),Xr(train_trials,:)]\ ...
            curr_str(train_trials,curr_t);
        
        ctx_str_cv(test_trials) = [curr_ctx_L(test_trials,curr_t),ones(sum(test_trials),1),Xr(test_trials,:)]*k_cv;
                
        % Do cross-time weights
        % (get W/mu from timepoint train trials)
        k = [curr_ctx_L(train_trials,curr_t),ones(sum(train_trials),1)]\ ...
            curr_str(train_trials,curr_t);
        
        % (apply W/mu to all time test sets)
        ctx_str_test = reshape([reshape(curr_ctx_L(test_trials,:)',[],1), ...
            ones(numel(curr_ctx_L(test_trials,:)),1)]*k(1:2),length(t),[])';
        
        % (estimate Xc on all test sets from residual)
        Xc = Xr(test_trials,:)\(curr_str(test_trials,:) - ctx_str_test);  
        Xc_trial = Xc(contrastside_idx(test_trials),:);
        
        ctx_str_xc(test_trials,:) = ctx_str_test + Xc_trial;
        
    end
        
    % Get R^2 for each time point
    sse_total = nansum((curr_str - nanmean(curr_str,1)).^2);
    sse_residual_cv = nansum((curr_str(:,curr_t) - ctx_str_cv).^2);
    sse_residual_xc = nansum((curr_str - ctx_str_xc).^2);
        
    expl_var_cv(curr_t) = 1 - (sse_residual_cv./sse_total(curr_t));
    expl_var_xc(curr_t,:) = 1 - (sse_residual_xc./sse_total);
    
    AP_print_progress_fraction(curr_t,length(t));
    
end

figure;

subplot(6,1,1:4);
imagesc(t,t,expl_var_xc);
caxis([-max(abs(caxis)),max(abs(caxis))]);
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');
line(xlim,xlim,'color','k');
ylabel('Trained time');
xlabel('Tested time');
c = colorbar; ylabel(c,'Expl var');
colormap(brewermap([],'*RdBu'));
axis square

subplot(6,1,5); hold on;
plot(t,nanmean(expl_var_cv,1),'k','linewidth',2);
plot(t,nanmean(expl_var_xc,1),'r','linewidth',2);
xlabel('Time');
ylabel('Expl var')
legend({'Cross-val','Cross weight'});

stim_t = find(t > 0,1);
line(repmat(t(stim_t),2,1),ylim,'color','r');
subplot(6,1,6);
plot(t,expl_var_xc(:,stim_t) - expl_var_cv(stim_t),'k','linewidth',2);
xlabel('Trained weights');
ylabel('\DeltaExpl var');
line(repmat(t(stim_t),2,1),ylim,'color','r');
line(xlim,[0,0],'color','k');


%% Ctx -> Str: compare models with different ROIs and timings

plot_str = 2; 

% % Set alignment shifts
% t_leeway = -t(1);
% leeway_samples = round(t_leeway*(sample_rate));
% stim_align = zeros(size(trial_contrast_allcat));
% move_align = -move_idx + leeway_samples;
% outcome_align = -outcome_idx + leeway_samples;
% 
% use_align = stim_align;
% 
% use_trials = move_t < 0.5 & ...
%     ~any(isnan(mua_allcat(:,:,plot_str)),2) & ...
%     ~any(isnan(fluor_roi_deconv(:,:,plot_ctx_L)),2);
% 
% % Re-align activity
% curr_str = cell2mat(arrayfun(@(trial) circshift( ...
%     mua_allcat(trial,:,plot_str), ...
%     use_align(trial),2),find(use_trials),'uni',false));
% 
% curr_ctx = cell2mat(arrayfun(@(trial) circshift( ...
%     fluor_roi_deconv(trial,:,:), ...
%     use_align(trial),2),find(use_trials),'uni',false));


% Set discontinuities in trial data
trial_discontinuities = [true(size(mua_allcat,1),1), ...
    false(size(mua_allcat,1),size(mua_allcat,2)-1)];

% All ROIs together
[k,ctxpred_spikes,explained_var] = ...
    AP_regresskernel( ...
    reshape(permute(fluor_roi_deconv,[2,1,3]),[],n_rois)', ...
    reshape(mua_allcat(:,:,plot_str)',[],1)', ...
    0,0,[false,false],10,1,1,reshape(trial_discontinuities',[],1)');

expl_var_allrois_nolag = explained_var.total;

% All L ROIs together
[k,ctxpred_spikes,explained_var] = ...
    AP_regresskernel( ...
    reshape(permute(fluor_roi_deconv(:,:,1:n_rois/2),[2,1,3]),[],n_rois/2)', ...
    reshape(mua_allcat(:,:,plot_str)',[],1)', ...
    0,0,[false,false],10,1,1,reshape(trial_discontinuities',[],1)');

expl_var_Lrois_nolag = explained_var.total;

% Single ROIs, no time shift
expl_var_singleroi_nolag = nan(n_rois,1);
for curr_roi = 1:n_rois
    
    [k,ctxpred_spikes,explained_var] = ...
        AP_regresskernel( ...
        reshape(permute(fluor_roi_deconv(:,:,curr_roi),[2,1,3]),[],1)', ...
        reshape(mua_allcat(:,:,plot_str)',[],1)', ...
        0,0,[false,false],10,1,1,reshape(trial_discontinuities',[],1)');
    
    expl_var_singleroi_nolag(curr_roi) = explained_var.total;
    AP_print_progress_fraction(curr_roi,n_rois);
    
end

% Single ROIs, time shift
t_shifts = -5:5;
expl_var_singleroi_lags = nan(n_rois,1);
for curr_roi = 1:n_rois
    
    [k,ctxpred_spikes,explained_var] = ...
        AP_regresskernel( ...
        reshape(permute(fluor_roi_deconv(:,:,curr_roi),[2,1,3]),[],1)', ...
        reshape(mua_allcat(:,:,plot_str)',[],1)', ...
        t_shifts,0,[false,false],10,1,1,reshape(trial_discontinuities',[],1)');
    
    expl_var_singleroi_lags(curr_roi) = explained_var.total;
    AP_print_progress_fraction(curr_roi,n_rois);
    
end

% Left ROIs (adding from most to least explained variance)
[~,L_roi_rank] = sort(expl_var_singleroi_lags(1:n_rois/2),'descend');

t_shifts = -5:5;
expl_var_multroiL_lags = nan(n_rois/2,1);
for curr_roi = 1:n_rois/2
    
    [k,ctxpred_spikes,explained_var] = ...
        AP_regresskernel( ...
        reshape(permute(fluor_roi_deconv(:,:,L_roi_rank(1:curr_roi)),[2,1,3]),[],curr_roi)', ...
        reshape(mua_allcat(:,:,plot_str)',[],1)', ...
        t_shifts,0,[false,false],10,1,1,reshape(trial_discontinuities',[],1)');
    
    expl_var_multroiL_lags(curr_roi) = explained_var.total;
    AP_print_progress_fraction(curr_roi,n_rois/2);
    
end

% All ROIs (adding from most to least explained variance)
[~,roi_rank] = sort(expl_var_singleroi_lags,'descend');

t_shifts = -5:5;
expl_var_multroi_lags = nan(n_rois,1);
for curr_roi = 1:n_rois
    
    [k,ctxpred_spikes,explained_var] = ...
        AP_regresskernel( ...
        reshape(permute(fluor_roi_deconv(:,:,roi_rank(1:curr_roi)),[2,1,3]),[],curr_roi)', ...
        reshape(mua_allcat(:,:,plot_str)',[],1)', ...
        t_shifts,0,[false,false],10,1,1,reshape(trial_discontinuities',[],1)');
    
    expl_var_multroi_lags(curr_roi) = explained_var.total;
    AP_print_progress_fraction(curr_roi,n_rois);
    
end


% Plot explained variance
figure; hold on;

plot(expl_var_singleroi_nolag(roi_rank));
plot(expl_var_singleroi_lags(roi_rank));
plot(find(ismember(roi_rank,L_roi_rank)),expl_var_multroiL_lags)
plot(expl_var_multroi_lags);
line(xlim,repmat(expl_var_Lrois_nolag,2,1),'color','k');
line(xlim,repmat(expl_var_allrois_nolag,2,1),'color','r');

ylabel('Explained variance');
xlabel('ROI (sorted by individual explained variance');
set(gca,'XTick',1:n_rois,'XTickLabel',{wf_roi(roi_rank).area})
legend({'Single ROI, no lag','Single ROI, lag', ...
    'Cumulative L ROIs, lag','Cumulative all ROIs, lag', ...
    'Left rois, no lag', 'All rois, no lag'})
title(['Str ' num2str(plot_str)]);






%%%%%%% try adding in Xc?

[k,ctxpred_spikes,explained_var] = ...
    AP_regresskernel( ...
    reshape(permute(fluor_roi_deconv,[2,1,3]),[],n_rois)', ...
    reshape(mua_allcat(:,:,plot_str)',[],1)', ...
    t_shifts,0,[false,false],10,1,1,reshape(trial_discontinuities',[],1)');

curr_contrastside = trial_contrastside_allcat;
[~,contrastside_idx] = ismember(curr_contrastside,stim_used,'rows');

Xc = grpstats(reshape(reshape(mua_allcat(:,:,plot_str)',[],1)' - ...
    ctxpred_spikes,length(t),[])',trial_contrastside_allcat);
Xc_trial = Xc(contrastside_idx,:);

curr_mua = reshape(mua_allcat(:,:,plot_str)',[],1)';
nonan = ~isnan(curr_mua) & ~isnan(ctxpred_spikes);

sse_total = sum((curr_mua(nonan) - nanmean(curr_mua(nonan))).^2);
sse_residual = sum((curr_mua(nonan) - ctxpred_spikes(nonan)).^2);
ev = 1-(sse_residual./sse_total);

ctxpred_spikes_xc = ctxpred_spikes + reshape(Xc_trial',[],1)';
sse_residual_xc = sum((curr_mua(nonan) - ctxpred_spikes_xc(nonan)).^2);
ev_xc = 1-(sse_residual_xc./sse_total);

%% Str vs ctx-predicted str with flipped binnings and ctx/ctx+task models
% (like fig 4g)

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(mua_allcat,1),1);
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

%%%%%% TESTING NEW ALIGNMENT: 
% max stim response for each contrast (better than stim onset)
str1_stimresponses = grpstats(mua_allcat(:,:,1),trial_contrastside_allcat,@nanmean);
[~,stimresponse_max_idx] = max(str1_stimresponses(:,t < 0.5),[],2);
stimresponse_max_idx(1:6) = stimresponse_max_idx(7);
[~,stim_idx] = ismember(trial_contrastside_allcat,unique(trial_contrastside_allcat),'rows');
stim_max_align = -stimresponse_max_idx(stim_idx) + leeway_samples;
%%%%%%

% Set windows to average activity

% timeavg_labels = {'Stim','Move onset','Outcome'};
% timeavg_t = {[0.05,0.15],[-0.05,0.05],[0.05,0.15]};
% timeavg_align = {stim_align,move_align,outcome_align};
% timeavg_trial_conditions = ...
%     {[sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1], ...
%     [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
%     [trial_outcome_allcat == 1, trial_outcome_allcat == -1]};
% timeavg_task_reduction = [1,2,4];

timeavg_labels = {'Stim'};
timeavg_t = {[0.05,0.15]};
timeavg_align = {stim_align};
timeavg_trial_conditions = ...
    {[trial_contrastside_allcat > 0,trial_contrastside_allcat < 0]};
timeavg_task_reduction = [1];

% timeavg_labels = {'Pre-stim','Stim','Post-stim','Post-stim+'};
% timeavg_t = {[-0.15,-0.05],[0.05,0.15],[0.2,0.3],[0.8,0.9]};
% timeavg_align = {stim_align,stim_align,stim_align,stim_align};
% timeavg_trial_conditions = ...
%     {[trial_contrastside_allcat > 0,trial_contrastside_allcat < 0], ...
%     [trial_contrastside_allcat > 0,trial_contrastside_allcat < 0], ...
%     [trial_contrastside_allcat > 0,trial_contrastside_allcat < 0], ...
%     [trial_contrastside_allcat > 0,trial_contrastside_allcat < 0]};
% timeavg_task_reduction = [1,1,1,1];

% timeavg_labels = {'Move'};
% timeavg_t = {[-0.05,0.05]};
% timeavg_align = {move_align};
% timeavg_trial_conditions = ...
%     {[trial_choice_allcat == -1,trial_choice_allcat == 1]};
% timeavg_task_reduction = [2];

% timeavg_labels = {'Go cue'};
% timeavg_t = {[0.55,0.65]};
% timeavg_align = {stim_align};
% timeavg_trial_conditions = ...
%     {[move_t < 0.5,move_t >= 0.5]};
% timeavg_task_reduction = [1];

% timeavg_labels = {'Outcome'};
% timeavg_t = {[-0.05,0.05]};
% timeavg_align = {outcome_align};
% timeavg_trial_conditions = ...
%     {[trial_outcome_allcat == 1,trial_outcome_allcat == -1]};
% timeavg_task_reduction = [4];

% Set activity percentiles and bins
act_prctile = [10,90];
n_act_bins = 5;

% Set areas and conditions
% plot_areas = [1,2,3,4];
plot_areas = [1];

% Loop across area pairs, plot binned predicted v measured activity
curr_act_allcat = mua_allcat;
curr_act_taskpred_reduced_allcat = mua_taskpred_reduced_allcat;

% (ctx-predicted)
curr_act_pred_allcat = mua_ctxpred_allcat;
curr_act_pred_taskpred_reduced_allcat = mua_ctxpred_taskpred_reduced_allcat;

% % (task-predicted)
% curr_act_pred_allcat = mua_taskpred_allcat;
% curr_act_pred_taskpred_reduced_allcat = mua_taskpred_reduced_allcat;

% % (ctx + task predicted)
% task_fix = mua_taskpred_allcat - mua_ctxpred_taskpred_allcat;
% task_fix_reduced = mua_taskpred_reduced_allcat - mua_ctxpred_taskpred_reduced_allcat;
% curr_act_pred_allcat = mua_ctxpred_allcat + task_fix;
% curr_act_pred_taskpred_reduced_allcat = mua_ctxpred_taskpred_reduced_allcat + task_fix_reduced;

% % (ctx + task predicted (by event))
% fix_reduced = 1;
% task_fix = (mua_taskpred_allcat - mua_taskpred_reduced_allcat(:,:,:,fix_reduced)) - ...
%     (mua_ctxpred_taskpred_allcat - mua_ctxpred_taskpred_reduced_allcat(:,:,:,fix_reduced));
% curr_act_pred_allcat = mua_ctxpred_allcat + task_fix;

% (ctx & task predicted) (in trial_activity_choiceworld_ctxtaskpred)
% curr_act_pred_allcat = mua_ctxtaskpred_allcat;

for curr_area_idx = 1:length(plot_areas)   
    
    plot_area = plot_areas(curr_area_idx);
    measured_v_pred_fig = figure('color','w','Name',['Str ' num2str(plot_area)]);
    
    for curr_timeavg = 1:length(timeavg_labels)
        
        trial_conditions = timeavg_trial_conditions{curr_timeavg};
        curr_task_reduction = timeavg_task_reduction(curr_timeavg);
        
        % (re-align and split activity)
        act_title = timeavg_labels{curr_timeavg};
        
        curr_act = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        curr_act_pred = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_pred_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
%         % (re-align, REDUCE, and split activity)
%         act_title = [timeavg_labels{curr_timeavg} ' (' task_regressor_labels{curr_task_reduction} '-reduced)'];
%         
%         curr_act = mat2cell(...
%             cell2mat(arrayfun(@(trial) circshift( ...
%             curr_act_allcat(trial,:,plot_area) - curr_act_taskpred_reduced_allcat(trial,:,plot_area,curr_task_reduction), ...
%             timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_allcat,1)),'uni',false)), ...
%             use_split,length(t));
%         
%         curr_act_pred = mat2cell(...
%             cell2mat(arrayfun(@(trial) circshift( ...
%             curr_act_pred_allcat(trial,:,plot_area) - curr_act_pred_taskpred_reduced_allcat(trial,:,plot_area,curr_task_reduction), ...
%             timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
%             use_split,length(t));
        
        trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
        
        % (get average activity within window)
        curr_event_t = t >= timeavg_t{curr_timeavg}(1) & t <= timeavg_t{curr_timeavg}(2);
        curr_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act,'uni',false);
        curr_act_pred_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_pred,'uni',false);

        % (bin predicted data across percentile range)
        act_bin_range = prctile(cell2mat(curr_act_avg),act_prctile);
        act_bin_edges = linspace(act_bin_range(1),act_bin_range(2),n_act_bins+1);
        act_bin_centers = act_bin_edges(1:end-1) + diff(act_bin_edges)./2;       
        act_trial_bins = cellfun(@(x) discretize(x,act_bin_edges),curr_act_avg,'uni',false);
        
        % (bin predicted data across percentile range)
        pred_bin_range = prctile(cell2mat(curr_act_pred_avg),act_prctile);
        pred_bin_edges = linspace(pred_bin_range(1),pred_bin_range(2),n_act_bins+1);
        pred_bin_centers = pred_bin_edges(1:end-1) + diff(pred_bin_edges)./2;        
        pred_trial_bins = cellfun(@(x) discretize(x,pred_bin_edges),curr_act_pred_avg,'uni',false); 
        
        % (get activity binned by measured)        
        act_use_trials = cellfun(@(act,act_taskpred,trial_bins) ...
            squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
            squeeze(~any(isnan(act_taskpred(:,curr_event_t)),2)) & ...
            ~isnan(trial_bins), ...
            curr_act,curr_act_pred,act_trial_bins,'uni',false);
        
        act_actbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,act_trial_bins,trial_conditions_exp,act_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_actbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,act_trial_bins,trial_conditions_exp,act_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        % (get activity binned by predicted)        
        pred_use_trials = cellfun(@(act,act_taskpred,trial_bins) ...
            squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
            squeeze(~any(isnan(act_taskpred(:,curr_event_t)),2)) & ...
            ~isnan(trial_bins), ...
            curr_act,curr_act_pred,pred_trial_bins,'uni',false);
        
        act_predbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_predbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        % Plot binned measured, predicted, and error (by predicted bins)
        measured_pred_fig = figure('color','w','Name', ...
            ['Str ' num2str(plot_area)' ', ' timeavg_labels{curr_timeavg}]);
        n_col_bins = n_act_bins + 2;
        
        [binned_act_pred_t,binned_act_pred_grp] = grpstats(cell2mat(curr_act_pred), ...
            [cell2mat(pred_trial_bins),trial_conditions],{'nanmean','gname'});
        binned_act_pred_grp = cellfun(@str2num,binned_act_pred_grp);        
        
        [binned_act_t,binned_act_grp] = grpstats(cell2mat(curr_act), ...
            [cell2mat(pred_trial_bins),trial_conditions],{'nanmean','gname'});
        binned_act_grp = cellfun(@str2num,binned_act_grp);
        
        binned_act_t_error = binned_act_t - binned_act_pred_t;
        
        % (plot predicted data)
        for curr_cond = 1:size(trial_conditions,2)           
            subplot(3,size(trial_conditions,2), ...
                sub2ind(fliplr([3,size(trial_conditions,2)]),curr_cond,1)); hold on;
            set(gca,'ColorOrder',[brewermap(n_col_bins,'*Greens')]);
            plot(t,binned_act_pred_t(binned_act_pred_grp(:,curr_cond + 1) == 1,:)','linewidth',2);
            xlabel('Time'); ylabel('Predicted data'); 
            title(['Condition ' num2str(curr_cond)]);
            line(repmat(timeavg_t{curr_timeavg}(1),2,1),ylim,'color','k');
            line(repmat(timeavg_t{curr_timeavg}(2),2,1),ylim,'color','k');            
        end
 
        % (plot measured data)        
        for curr_cond = 1:size(trial_conditions,2)           
            subplot(3,size(trial_conditions,2), ...
                sub2ind(fliplr([3,size(trial_conditions,2)]),curr_cond,2)); hold on;
            set(gca,'ColorOrder',[brewermap(n_col_bins,'*Greys')]);
            plot(t,binned_act_t(binned_act_grp(:,curr_cond + 1) == 1,:)','linewidth',2);
            xlabel('Time'); ylabel('Measured data'); 
            title(['Condition ' num2str(curr_cond)]);
            line(repmat(timeavg_t{curr_timeavg}(1),2,1),ylim,'color','k');
            line(repmat(timeavg_t{curr_timeavg}(2),2,1),ylim,'color','k');            
        end        
        
        % (plot error)      
        for curr_cond = 1:size(trial_conditions,2)           
            subplot(3,size(trial_conditions,2), ...
                sub2ind(fliplr([3,size(trial_conditions,2)]),curr_cond,3)); hold on;
            set(gca,'ColorOrder',[brewermap(n_col_bins,'*OrRd')]);
            plot(t,binned_act_t_error(binned_act_grp(:,curr_cond + 1) == 1,:)','linewidth',2);
            xlabel('Time'); ylabel('Prediction error'); 
            title(['Condition ' num2str(curr_cond)]);
            line(repmat(timeavg_t{curr_timeavg}(1),2,1),ylim,'color','k');
            line(repmat(timeavg_t{curr_timeavg}(2),2,1),ylim,'color','k');
            line(xlim,[0,0],'color','k')
        end          
        
        linkaxes(get(measured_pred_fig,'Children'),'xy');       
        
        % Plot measured v predicted in bins
        figure(measured_v_pred_fig)
        
        % (measured vs binned predicted)
        subplot(2,length(timeavg_labels), ...
            sub2ind(fliplr([2,length(timeavg_labels)]),curr_timeavg,1));
        hold on;
        
        errorbar( ...
            squeeze(nanmean(act_pred_predbinmean,2)), ...
            squeeze(nanmean(act_predbinmean,2)), ...
            squeeze(AP_sem(act_predbinmean,2)), ...
            'linewidth',2,'CapSize',0);
        xlabel(['Predicted (' num2str(plot_area) ')']);
        ylabel(['Measured (' num2str(plot_area) ')'])
        ylim(xlim); axis square;
        title(act_title);
        line(xlim,xlim,'color','k');
              
        % (predicted vs binned measured)
        subplot(2,length(timeavg_labels), ...
            sub2ind(fliplr([2,length(timeavg_labels)]),curr_timeavg,2));
        hold on;
        
        errorbar( ...
            squeeze(nanmean(act_actbinmean,2)), ...
            squeeze(nanmean(act_pred_actbinmean,2)), ...
            squeeze(AP_sem(act_pred_actbinmean,2)), ...
            'linewidth',2,'CapSize',0);
        xlabel(['Measured (' num2str(plot_area) ')']);
        ylabel(['Predicted (' num2str(plot_area) ')'])
        ylim(xlim); axis square;
        title([act_title]);
        line(xlim,xlim,'color','k');
        
    end
end

linkaxes(get(measured_v_pred_fig,'Children'),'xy');


%% Str vs ctx-predicted str (dots = measured, lines = model)
% (like fig 4g)

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(mua_allcat,1),1);
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% %%%%%% TESTING NEW ALIGNMENT: 
% % max stim response for each contrast (better than stim onset)
% str1_stimresponses = grpstats(mua_allcat(:,:,1),trial_contrastside_allcat,@nanmean);
% [~,stimresponse_max_idx] = max(str1_stimresponses(:,t < 0.5),[],2);
% stimresponse_max_idx(1:6) = stimresponse_max_idx(7);
% [~,stim_idx] = ismember(trial_contrastside_allcat,unique(trial_contrastside_allcat),'rows');
% stim_max_align = -stimresponse_max_idx(stim_idx) + leeway_samples;
% %%%%%%

% Set windows to average activity

timeavg_labels = {'Stim','Move onset','Outcome'};
timeavg_t = {[0.05,0.15],[-0.05,0.05],[0.05,0.15]};
timeavg_align = {stim_align,move_align,outcome_align};
timeavg_trial_conditions = ...
    {[sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1], ...
    [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
    [trial_outcome_allcat == 1, trial_outcome_allcat == -1]};

% timeavg_labels = {'Stim'};
% timeavg_t = {[0.05,0.15]};
% timeavg_align = {stim_align};
% timeavg_trial_conditions = ...
%     {[trial_contrastside_allcat > 0,trial_contrastside_allcat < 0]};

% timeavg_labels = {'Stim'};
% timeavg_t = {[0.05,0.15]};
% timeavg_align = {stim_align};
% timeavg_trial_conditions = ...
%     {[trial_stim_allcat == 1,trial_stim_allcat == 2,trial_stim_allcat == 3]};

% timeavg_labels = {'Pre-stim','Stim','Post-stim','Post-stim+'};
% timeavg_t = {[-0.15,-0.05],[0.05,0.15],[0.2,0.3],[0.8,0.9]};
% timeavg_align = {stim_align,stim_align,stim_align,stim_align};
% timeavg_trial_conditions = ...
%     {[trial_contrastside_allcat > 0,trial_contrastside_allcat < 0], ...
%     [trial_contrastside_allcat > 0,trial_contrastside_allcat < 0], ...
%     [trial_contrastside_allcat > 0,trial_contrastside_allcat < 0], ...
%     [trial_contrastside_allcat > 0,trial_contrastside_allcat < 0]};

% timeavg_labels = {'Move'};
% timeavg_t = {[-0.05,0.05]};
% timeavg_align = {move_align};
% timeavg_trial_conditions = ...
%     {[trial_choice_allcat == -1,trial_choice_allcat == 1]};

% timeavg_labels = {'Go cue'};
% timeavg_t = {[0.55,0.65]};
% timeavg_align = {stim_align};
% timeavg_trial_conditions = ...
%     {[move_t < 0.5,move_t >= 0.5]};

% timeavg_labels = {'Outcome'};
% timeavg_t = {[-0.05,0.05]};
% timeavg_align = {outcome_align};
% timeavg_trial_conditions = ...
%     {[trial_outcome_allcat == 1,trial_outcome_allcat == -1]};

% Set activity percentiles and bins
act_prctile = [10,90];
n_act_bins = 5;

% Set areas and conditions
plot_areas = [1,2,3,4];
% plot_areas = [1];

% Loop across area pairs, plot binned predicted v measured activity
curr_act_allcat = mua_allcat;

% (ctx-predicted)
curr_act_pred_allcat = mua_ctxpred_allcat;

% (fix by all events)
task_fix = mua_taskpred_allcat - mua_ctxpred_taskpred_allcat;
curr_act_pred_fix_allcat = mua_ctxpred_allcat + task_fix;

% (fix by specific events)
% fix_reduced = 2;
% task_fix = (mua_taskpred_allcat - mua_taskpred_reduced_allcat(:,:,:,fix_reduced)) - ...
%     (mua_ctxpred_taskpred_allcat - mua_ctxpred_taskpred_reduced_allcat(:,:,:,fix_reduced));
% curr_act_pred_fix_allcat = mua_ctxpred_allcat + task_fix;

% % (for passive)
% % DIRTY! need to get offset for each experiment
% stim_act = grpstats(curr_act_allcat(:,:,1),trial_contrastside_allcat);
% stim_act_pred = grpstats(curr_act_pred_allcat(:,:,1),trial_contrastside_allcat);
% stim_act_offset = stim_act - stim_act_pred;
% [~,contrastside_idx] = ismember(trial_contrastside_allcat,unique(trial_contrastside_allcat),'rows');
% task_fix = stim_act_offset(contrastside_idx,:);
% curr_act_pred_fix_allcat = mua_ctxpred_allcat + task_fix;

for curr_area_idx = 1:length(plot_areas)   
    
    plot_area = plot_areas(curr_area_idx);
    measured_v_pred_fig = figure('color','w','Name',['Str ' num2str(plot_area)]);
    
    for curr_timeavg = 1:length(timeavg_labels)
        
        trial_conditions = timeavg_trial_conditions{curr_timeavg};
        
        % (re-align and split activity)
        act_title = timeavg_labels{curr_timeavg};
        
        curr_act = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        curr_act_pred = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_pred_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        curr_act_predfix = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_pred_fix_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_fix_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
        
        % (get average activity within window)
        curr_event_t = t >= timeavg_t{curr_timeavg}(1) & t <= timeavg_t{curr_timeavg}(2);
        curr_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act,'uni',false);
        curr_act_pred_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_pred,'uni',false);
        curr_act_predfix_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_predfix,'uni',false);
        
        % (bin predicted data across percentile range)
        pred_bin_range = prctile(cell2mat(curr_act_pred_avg),act_prctile);
        pred_bin_edges = linspace(pred_bin_range(1),pred_bin_range(2),n_act_bins+1);
        pred_bin_centers = pred_bin_edges(1:end-1) + diff(pred_bin_edges)./2;        
        pred_trial_bins = cellfun(@(x) discretize(x,pred_bin_edges),curr_act_pred_avg,'uni',false);         
        
        % (get activity binned by predicted)        
        pred_use_trials = cellfun(@(act,act_taskpred,trial_bins) ...
            squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
            squeeze(~any(isnan(act_taskpred(:,curr_event_t)),2)) & ...
            ~isnan(trial_bins), ...
            curr_act,curr_act_pred,pred_trial_bins,'uni',false);
        
        act_predbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_predbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        % (get "fixed" predicted activity binned by predicted)
        act_predfix_predbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_predfix_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        % Plot binned measured, predicted, and error (by predicted bins)
        measured_pred_fig = figure('color','w','Name', ...
            ['Str ' num2str(plot_area)' ', ' timeavg_labels{curr_timeavg}]);
        n_col_bins = n_act_bins + 2;
        
        [binned_act_pred_t,binned_act_pred_grp] = grpstats(cell2mat(curr_act_pred), ...
            [cell2mat(pred_trial_bins),trial_conditions],{'nanmean','gname'});
        binned_act_pred_grp = cellfun(@str2num,binned_act_pred_grp);        
        
        [binned_act_t,binned_act_grp] = grpstats(cell2mat(curr_act), ...
            [cell2mat(pred_trial_bins),trial_conditions],{'nanmean','gname'});
        binned_act_grp = cellfun(@str2num,binned_act_grp);
        
        binned_act_t_error = binned_act_t - binned_act_pred_t;
        
        % (plot predicted data)
        for curr_cond = 1:size(trial_conditions,2)           
            subplot(3,size(trial_conditions,2), ...
                sub2ind(fliplr([3,size(trial_conditions,2)]),curr_cond,1)); hold on;
            set(gca,'ColorOrder',[brewermap(n_col_bins,'*Greens')]);
            plot(t,binned_act_pred_t(binned_act_pred_grp(:,curr_cond + 1) == 1,:)','linewidth',2);
            xlabel('Time'); ylabel('Predicted data'); 
            title(['Condition ' num2str(curr_cond)]);
            line(repmat(timeavg_t{curr_timeavg}(1),2,1),ylim,'color','k');
            line(repmat(timeavg_t{curr_timeavg}(2),2,1),ylim,'color','k');            
        end
 
        % (plot measured data)        
        for curr_cond = 1:size(trial_conditions,2)           
            subplot(3,size(trial_conditions,2), ...
                sub2ind(fliplr([3,size(trial_conditions,2)]),curr_cond,2)); hold on;
            set(gca,'ColorOrder',[brewermap(n_col_bins,'*Greys')]);
            plot(t,binned_act_t(binned_act_grp(:,curr_cond + 1) == 1,:)','linewidth',2);
            xlabel('Time'); ylabel('Measured data'); 
            title(['Condition ' num2str(curr_cond)]);
            line(repmat(timeavg_t{curr_timeavg}(1),2,1),ylim,'color','k');
            line(repmat(timeavg_t{curr_timeavg}(2),2,1),ylim,'color','k');            
        end        
        
        % (plot error)      
        for curr_cond = 1:size(trial_conditions,2)           
            subplot(3,size(trial_conditions,2), ...
                sub2ind(fliplr([3,size(trial_conditions,2)]),curr_cond,3)); hold on;
            set(gca,'ColorOrder',[brewermap(n_col_bins,'*OrRd')]);
            plot(t,binned_act_t_error(binned_act_grp(:,curr_cond + 1) == 1,:)','linewidth',2);
            xlabel('Time'); ylabel('Prediction error'); 
            title(['Condition ' num2str(curr_cond)]);
            line(repmat(timeavg_t{curr_timeavg}(1),2,1),ylim,'color','k');
            line(repmat(timeavg_t{curr_timeavg}(2),2,1),ylim,'color','k');
            line(xlim,[0,0],'color','k')
        end          
        
        linkaxes(get(measured_pred_fig,'Children'),'xy');       
        
        % Plot measured v predicted in bins
        figure(measured_v_pred_fig)
        
        % (measured vs binned predicted)
        subplot(1,length(timeavg_labels), ...
            sub2ind(fliplr([1,length(timeavg_labels)]),curr_timeavg,1));
        hold on;
        set(gca,'ColorOrder',[1,0,0;0,0,1;0,1,0]);
        
        errorbar( ...
            squeeze(nanmean(act_pred_predbinmean,2)), ...
            squeeze(nanmean(act_predbinmean,2)), ...
            squeeze(AP_sem(act_predbinmean,2)), ...
            '.','MarkerSize',30,'linewidth',2,'CapSize',20);
%         plot( ...
%             squeeze(nanmean(act_pred_predbinmean,2)), ...
%             squeeze(nanmean(act_predfix_predbinmean,2)),'linewidth',2);        

cond_col = {'r','b','g'};
for curr_cond = 1:size(trial_conditions,2)
   AP_errorfill( ...
        squeeze(nanmean(act_pred_predbinmean(:,:,curr_cond),2)), ...
        squeeze(nanmean(act_predfix_predbinmean(:,:,curr_cond),2)), ...
        squeeze(AP_sem(act_predfix_predbinmean(:,:,curr_cond),2)),cond_col{curr_cond},0.5,true);
end
        
        xlabel(['Predicted (' num2str(plot_area) ')']);
        ylabel(['Measured (' num2str(plot_area) ')'])
        ylim(xlim); axis square;
        title(act_title);
        line(xlim,xlim,'color','k');
        legend({'Measured cond 1','Measured cond 2','','"Fixed" cond 1','','"Fixed" cond 2'});
        
        
    end
end

linkaxes(get(measured_v_pred_fig,'Children'),'xy');





















