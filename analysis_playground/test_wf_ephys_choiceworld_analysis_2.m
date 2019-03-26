% (this is to work on AP_ctx_str_figures and related code)

%% Load choiceworld trial data (** NEEDED FOR BELOW **)

% Load data
% data_fn = 'trial_activity_choiceworld';
% data_fn = 'trial_activity_choiceworld_MOVELTEST';
% data_fn = 'trial_activity_choiceworld_FWDCTXTEST';
% data_fn = 'trial_activity_choiceworld_STIMRTEST';
% data_fn = 'trial_activity_choiceworld_MOVEONGOINGTEST';
% data_fn = 'trial_activity_choiceworld_EARLYMOVETEST';
% data_fn = 'trial_activity_choiceworld_REGRESSRESIDUAL';
% data_fn = 'trial_activity_choiceworld_CTXSTRMOVEONLY';
% data_fn = 'trial_activity_choiceworld_MOVExSTIM';
% data_fn = 'trial_activity_choiceworld_MSN';
% data_fn = 'trial_activity_choiceworld_FSI';
data_fn = 'trial_activity_choiceworld_200umdepth';

exclude_data = true;
AP_load_concat_normalize_ctx_str;


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
% max_vel = AP_signed_max(wheel_velocity_allcat_move(:,t > 0 & t < 0.5),2);

% (to use summed velocity regardless of final choice)
% max_vel = sum(wheel_velocity_allcat_move(:,t > 0 & t < 0.5),2);

% (to use maximum cumulative velocity regardless of final choice)
% max_vel = AP_signed_max(cumsum(wheel_velocity_allcat_move(:,t > 0 & t < 0.5),2),2);

% (to use mousecam movement signed by choice)
max_vel = sum(movement_allcat_move(:,t > 0 & t < 0.5),2).*trial_choice_allcat;

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
    AP_image_scroll(curr_k_px,t_shifts{curr_regressor});
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

plot_depth = 2;


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
    
    %     use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} <= 0.5;
    
    use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} <= 0.5 & ...
        trial_contrast_allcat_exp{curr_exp} > 0 & ...
        trial_side_allcat_exp{curr_exp} == 1 & ...
        trial_choice_allcat_exp{curr_exp} == -1;
    
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
    
    %      (to align to movement)
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

AP_image_scroll([muafluor_corr_px_mean,fluor2mua_corr_px_mean,mua2fluor_corr_px_mean],t(use_t));
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


move_trial = any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > 0.02,2);
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

AP_image_scroll([muafluor_corr_px_mean,fluor2mua_corr_px_mean,mua2fluor_corr_px_mean],t(use_t));
axis image
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k',[],[size(U_master,1),size(U_master,2),1,3]);


%% ~~~~~~~~ EXPLORATORY ANALYSIS AGAIN FOR UNITS ~~~~~~~~~~~~~~


%% Get PSTH for each unit

n_units = size(templates,1);

% PSTH timing
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t = conv2(t_bins,[1,1]/2,'valid');

% Event to align
use_align = wheel_move_time';
% use_align = stimOn_times;
use_align(isnan(use_align)) = [];
t_peri_event = use_align + t_bins;

% Get PSTHs
unit_psth = nan(n_units,length(t_bins)-1);
for curr_unit = 1:n_units
    curr_spikes = spike_times_timeline(spike_templates == curr_unit);
    curr_spikes_binned = cell2mat(arrayfun(@(x) ...
        histcounts(curr_spikes,t_peri_event(x,:)), ...
        [1:size(t_peri_event,1)]','uni',false))./psth_bin_size;
    unit_psth(curr_unit,:) = sum(curr_spikes_binned,1);
end

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
unit_psth_smoothed = convn(unit_psth,smWin,'same');

% PCA of PSTHs
[coeff,score,latent] = pca(zscore(unit_psth_smoothed,[],2)');

depth_bins = linspace(0,4000,20);
depth_centers = depth_bins(1:end-1) + diff(depth_bins)./2;
depth_group = discretize(template_depths,depth_bins);

plot_coeffs = 1:4;

figure;
subplot(1,2,1); hold on;
plot(t,score(:,plot_coeffs),'linewidth',2);
line([0,0],ylim,'color','k');
xlabel('Time');
ylabel('Activity');

subplot(1,2,2,'YDir','reverse'); hold on;
for curr_coeff = plot_coeffs
    curr_depth_coeff = accumarray(depth_group,coeff(:,curr_coeff), ...
        [length(depth_bins)-1,1],@nanmean,nan);
    plot(smooth(curr_depth_coeff,3),depth_centers,'linewidth',2);
end
xlabel('Coeff loading');
ylabel('Depth');


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
[move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
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
unit_expl_var = nan(length(regressors)+1,max(spike_templates));
unit_taskpred_k = cell(length(regressors)+return_constant,max(spike_templates));
unit_taskpred = nan(size(event_aligned_unit));
unit_taskpred_reduced = ...
    repmat(nan(size(event_aligned_unit)),1,1,1,length(regressors));
for curr_unit = 1:max(spike_templates)
    
    baseline = nanmean(reshape(event_aligned_unit(:,t < 0,curr_unit),[],1)*raster_sample_rate);
    activity = single(binned_spikes(curr_unit,:)) - baseline;
    
    % Skip if nothing in this depth
    if ~any(activity(:))
        continue
    end
    
    [task_kernel,activity_predicted,expl_var,activity_predicted_reduced] = ...
        AP_regresskernel(regressors,activity,sample_shifts, ...
        lambda,zs,cvfold,return_constant,use_constant);
    
    unit_expl_var(1,curr_unit) = expl_var.total;
    unit_expl_var(2:end,curr_unit) = expl_var.partial;
    
    unit_taskpred_k(:,curr_unit) = task_kernel;
    
    unit_taskpred(:,:,curr_unit) = ...
        interp1(time_bin_centers,activity_predicted',t_peri_event)./raster_sample_rate;
    
    unit_taskpred_reduced(:,:,curr_unit,:) = cell2mat(arrayfun(@(x) ...
        interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
        t_peri_event)./raster_sample_rate,permute(1:length(regressors),[1,3,4,2]),'uni',false));
    
    AP_print_progress_fraction(curr_unit,max(spike_templates));
    
end


%%%%% try just doing a comparison of aligned peaks?





