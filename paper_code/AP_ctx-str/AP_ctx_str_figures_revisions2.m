%% Code for organizing second round/final revisions



%% Celltype stim response (TANS are stim, not move) [UNUSED]

mua_exp = vertcat(mua_all{:});

% (move aligned)
move_idx_exp = mat2cell(move_idx,use_split,1);
mua_exp_movealign = vertcat(mua_all{:});
for curr_exp = 1:length(mua_exp_movealign)
   for curr_trial = 1:size(mua_exp_movealign{curr_exp},1)
       mua_exp_movealign{curr_exp}(curr_trial,:,:) = ...
           circshift(mua_exp_movealign{curr_exp}(curr_trial,:,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
   end
end

trial_stim_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
move_t_exp = mat2cell(move_t,use_split,1);
trial_choice_allcat_exp = mat2cell(trial_choice_allcat,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

% Get SUA average response by stim
unique_stim = unique(trial_stim_allcat);

sua_stim_avg = nan(length(good_units_allcat),length(t),length(unique_stim));
for curr_stim_idx = 1:length(unique_stim)
    curr_stim = unique_stim(curr_stim_idx);
    sua_stim_avg(:,:,curr_stim_idx) = cell2mat(cellfun(@(act,stim,rxn,choice,outcome) ...
        squeeze(nanmean(act(stim == curr_stim & rxn < 0.5 & outcome == 1,:,:),1)), ...
        mua_exp,trial_stim_allcat_exp,move_t_exp, ...
        trial_choice_allcat_exp,trial_outcome_allcat_exp,'uni',false)')';
end

% Get MUA time-window response by stim and domain/celltype
stim_avg_t = [0,0.2];
stim_avg_t_idx = t >= stim_avg_t(1) & t <= stim_avg_t(2);

sua_stim_avg_t = squeeze(nanmean(sua_stim_avg(:,stim_avg_t_idx,:),2));

grp_matrix = cat(3,permute(repmat([celltype_allcat(good_units_allcat), ...
    domain_aligned_allcat(good_units_allcat), ...
    recordings_allcat(good_units_allcat)],1,1,length(unique_stim)),[1,3,2]), ...
    repmat(1:length(unique_stim),sum(good_units_allcat),1));

sua_stim_avg_t_grp = accumarray(reshape(grp_matrix,[],size(grp_matrix,3)), ...
    reshape(sua_stim_avg_t(good_units_allcat,:),[],1), ...
    [max(celltype_allcat),n_aligned_depths,max(recordings_allcat),length(unique_stim)], ...
    @nanmean,NaN);

% Plot stim response by domain/celltype
plot_celltypes = 1:3;
stim_col = brewermap(11,'*RdBu');
stim_col(6,:) = [0.5,0.5,0.5];

figure;
for curr_depth = 1:n_aligned_depths
    for curr_celltype = plot_celltypes
        
        subplot(n_aligned_depths,length(plot_celltypes), ...
            (curr_depth-1)*length(plot_celltypes)+curr_celltype); hold on;
        set(gca,'ColorOrder',stim_col);
        
        curr_cells = domain_aligned_allcat == curr_depth & ...
            celltype_allcat == curr_celltype & good_units_allcat;
        
        curr_stim_avg = squeeze(nanmean(sua_stim_avg(curr_cells,:,:),1));
        plot(t,curr_stim_avg);
    end
end

% Plot stim spikes within time by stim
celltype_col = ...
    [0.9,0.4,0.6; ...
    0.4,0.5,0.6; ...
    0.5,0.3,0.1; ...
    1,0.5,0];

figure;
for curr_depth = 1:n_aligned_depths
    subplot(n_aligned_depths,1,curr_depth);
    hold on;
    for curr_celltype = plot_celltypes
        errorbar(unique_stim, ...
            squeeze(nanmean(sua_stim_avg_t_grp(curr_celltype,curr_depth,:,:),3)), ...
            squeeze(AP_sem(sua_stim_avg_t_grp(curr_celltype,curr_depth,:,:),3)), ...
            'linewidth',2,'color',celltype_col(curr_celltype,:));
    end
end



%% Cortex/task explained striatal variance - spike thinned [UNUSED]
% DOESN'T SEEM LIKE THIS WORKS YET, ALSO RECTIFYING MEANS THE RATE ISN'T
% EXACTLY THE SAME

% Use raw data (not normalized or baseline-subtracted) for expl var
mua_exp = vertcat(mua_all{:});
mua_taskpred_exp = vertcat(mua_taskpred_all{:});
mua_ctxpred_exp = vertcat(mua_ctxpred_all{:});

% Thin spikes to equalize rate across domains
mua_exp_thinned = cellfun(@(x) nan(size(x)),mua_exp,'uni',false);
mua_taskpred_exp_thinned = cellfun(@(x) nan(size(x)),mua_taskpred_exp,'uni',false);
mua_ctxpred_exp_thinned = cellfun(@(x) nan(size(x)),mua_ctxpred_exp,'uni',false);
for curr_exp = 1:max(split_idx)
    str_rates = nanmean(reshape(mua_exp{curr_exp},[],n_depths));
    str_var = var(reshape(mua_exp{curr_exp},[],n_depths));
    if isnan(str_rates(1))
        continue
    end
    
    for curr_depth = 1:n_depths
        if all(isnan(reshape(mua_exp{curr_exp}(:,:,curr_depth),[],1)))
            continue
        end
        
        % Pull off spikes until the variance is equal
        curr_mua_thinned = mua_exp{curr_exp}(:,:,curr_depth);
        curr_thin = zeros(size(curr_mua_thinned));
        while var(curr_mua_thinned(:)) >= min(str_var)
            curr_thin = curr_thin + poissrnd(100, ...
                size(mua_exp{curr_exp}(:,:,curr_depth)));
            curr_mua_thinned = max(0,mua_exp{curr_exp}(:,:,curr_depth)-curr_thin);            
        end
        
        mua_exp_thinned{curr_exp}(:,:,curr_depth) = ...
            max(0,mua_exp{curr_exp}(:,:,curr_depth) - curr_thin);
        mua_taskpred_exp_thinned{curr_exp}(:,:,curr_depth) = ...
            max(0,mua_taskpred_exp{curr_exp}(:,:,curr_depth) - curr_thin);
        mua_ctxpred_exp_thinned{curr_exp}(:,:,curr_depth) = ...
            max(0,mua_ctxpred_exp{curr_exp}(:,:,curr_depth) - curr_thin);       
    end    
end


% Get R^2 for task and cortex
taskpred_r2 = nan(max(split_idx),n_depths);
ctxpred_r2 = nan(max(split_idx),n_depths);
for curr_exp = 1:max(split_idx)
       
    curr_data = reshape(permute(mua_exp_thinned{curr_exp},[2,1,3]),[],n_depths);
    curr_data_baselinesub = reshape(permute(mua_exp_thinned{curr_exp},[2,1,3]),[],n_depths) - ...
        (nanmean(reshape(mua_exp_thinned{curr_exp}(:,t < 0,:),[],size(mua_exp_thinned{curr_exp},3)),1));
    curr_taskpred_data = reshape(permute(mua_taskpred_exp_thinned{curr_exp},[2,1,3]),[],n_depths);
    curr_ctxpred_data = reshape(permute(mua_ctxpred_exp_thinned{curr_exp},[2,1,3]),[],n_depths);
       
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_data_baselinesub) | ...
        isnan(curr_taskpred_data) | isnan(curr_ctxpred_data);
    curr_data(nan_samples) = NaN;
    curr_data_baselinesub(nan_samples) = NaN;
    curr_taskpred_data(nan_samples) = NaN;
    curr_ctxpred_data(nan_samples) = NaN;

    % (task regressed from average baseline-subtracted data)
    taskpred_r2(curr_exp,:) = 1 - (nansum((curr_data_baselinesub-curr_taskpred_data).^2,1)./ ...
        nansum((curr_data_baselinesub-nanmean(curr_data_baselinesub,1)).^2,1));
    % (cortex regressed from raw data)
    ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    
end

figure; hold on;
errorbar(nanmean(taskpred_r2,1),AP_sem(taskpred_r2,1),'b','linewidth',2,'CapSize',0);
errorbar(nanmean(ctxpred_r2,1),AP_sem(ctxpred_r2,1),'color',[0,0.6,0],'linewidth',2,'CapSize',0);
xlabel('Striatum depth');
ylabel('Explained variance');
legend({'Task','Cortex'});


% Plot explained variance task vs cortex by experiment
figure; hold on;
str_col = max(hsv(n_depths)-0.2,0);
for curr_str = 1:n_depths
    errorbar(squeeze(nanmean(taskpred_r2(:,curr_str),1)), ...
        squeeze(nanmean(ctxpred_r2(:,curr_str),1)), ...
        squeeze(AP_sem(ctxpred_r2(:,curr_str),1)),squeeze(AP_sem(ctxpred_r2(:,curr_str),1)), ...
        squeeze(AP_sem(taskpred_r2(:,curr_str),1)),squeeze(AP_sem(taskpred_r2(:,curr_str),1)), ...
        'color','k','linewidth',2);
    
    scatter(taskpred_r2(:,curr_str),ctxpred_r2(:,curr_str),10, ...
        str_col(curr_str,:),'filled');
    scatter(nanmean(taskpred_r2(:,curr_str),1), ...
        nanmean(ctxpred_r2(:,curr_str),1),80, ...
        str_col(curr_str,:),'filled','MarkerEdgeColor','k','linewidth',2);
end
axis tight equal;
line(xlim,xlim,'color','k','linestyle','--');
xlabel('Task R^2');
ylabel('Cortex R^2');
legend({'DMS','DCS','DLS'})

% (Task R2 statistics)
disp('Task R^2 values:');
for curr_depth = 1:n_depths
    disp(['Str ' num2str(curr_depth) ...
        ' R2 = ' num2str(nanmean(taskpred_r2(:,curr_depth))) '+/- ' ...
        num2str(AP_sem(taskpred_r2(:,curr_depth),1))]); 
end

% (Cortex R2 statistics)
disp('Cortex R^2 values:');
for curr_depth = 1:n_depths
    disp(['Str ' num2str(curr_depth) ...
        ' R2 = ' num2str(nanmean(ctxpred_r2(:,curr_depth))) '+/- ' ...
        num2str(AP_sem(ctxpred_r2(:,curr_depth),1))]); 
end

% (Cortex vs task R2 statistics)
disp('Cortex vs Task R^2 signrank:');
for curr_depth = 1:n_depths
    curr_p = signrank(ctxpred_r2(:,curr_depth), ...
        taskpred_r2(:,curr_depth));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end





% Plot explained variance task vs cortex R2 by experiment
mua_var_exp = log10(cell2mat(cellfun(@(x) var(reshape(x,[],n_depths)),mua_exp_thinned,'uni',false)));

figure; hold on;
str_col = max(hsv(n_depths)-0.2,0);
for curr_str = 1:n_depths
    errorbar(squeeze(nanmean(mua_var_exp(:,curr_str),1)), ...
        squeeze(nanmean(ctxpred_r2(:,curr_str),1)), ...
        squeeze(AP_sem(ctxpred_r2(:,curr_str),1)),squeeze(AP_sem(ctxpred_r2(:,curr_str),1)), ...
        squeeze(AP_sem(mua_var_exp(:,curr_str),1)),squeeze(AP_sem(mua_var_exp(:,curr_str),1)), ...
        'color','k','linewidth',2);
    
    scatter(mua_var_exp(:,curr_str),ctxpred_r2(:,curr_str),10, ...
        str_col(curr_str,:),'filled');
    scatter(nanmean(mua_var_exp(:,curr_str),1), ...
        nanmean(ctxpred_r2(:,curr_str),1),80, ...
        str_col(curr_str,:),'filled','MarkerEdgeColor','k','linewidth',2);
end
xlabel('log_{10}(striatum variance)');
ylabel('Cortex R^2');
legend({'DMS','DCS','DLS'})



%% TRYING AGAIN: just plot variance by explained variance? (fixed below)

% Use raw data (not normalized or baseline-subtracted) for expl var
mua_exp = vertcat(mua_all{:});
mua_taskpred_exp = vertcat(mua_taskpred_all{:});
mua_ctxpred_exp = vertcat(mua_ctxpred_all{:});

% Get R^2 for task and cortex
taskpred_r2 = nan(max(split_idx),n_depths);
ctxpred_r2 = nan(max(split_idx),n_depths);
for curr_exp = 1:max(split_idx)
       
    curr_data = reshape(permute(mua_exp{curr_exp},[2,1,3]),[],n_depths);
    curr_data_baselinesub = reshape(permute(mua_exp{curr_exp},[2,1,3]),[],n_depths) - ...
        (nanmean(reshape(mua_exp{curr_exp}(:,t < 0,:),[],size(mua_exp{curr_exp},3)),1));
    curr_taskpred_data = reshape(permute(mua_taskpred_exp{curr_exp},[2,1,3]),[],n_depths);
    curr_ctxpred_data = reshape(permute(mua_ctxpred_exp{curr_exp},[2,1,3]),[],n_depths);
       
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_data_baselinesub) | ...
        isnan(curr_taskpred_data) | isnan(curr_ctxpred_data);
    curr_data(nan_samples) = NaN;
    curr_data_baselinesub(nan_samples) = NaN;
    curr_taskpred_data(nan_samples) = NaN;
    curr_ctxpred_data(nan_samples) = NaN;

    % (task regressed from average baseline-subtracted data)
    taskpred_r2(curr_exp,:) = 1 - (nansum((curr_data_baselinesub-curr_taskpred_data).^2,1)./ ...
        nansum((curr_data_baselinesub-nanmean(curr_data_baselinesub,1)).^2,1));
    % (cortex regressed from raw data)
    ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    
end

% Plot explained variance task vs cortex R2 by experiment
mua_var_exp_log = log10(cell2mat(cellfun(@(x) var(reshape(x,[],n_depths)),mua_exp,'uni',false)));

figure; hold on;
str_col = max(hsv(n_depths)-0.2,0);
for curr_str = 1:n_depths
    errorbar(squeeze(nanmean(mua_var_exp_log(:,curr_str),1)), ...
        squeeze(nanmean(ctxpred_r2(:,curr_str),1)), ...
        squeeze(AP_sem(ctxpred_r2(:,curr_str),1)),squeeze(AP_sem(ctxpred_r2(:,curr_str),1)), ...
        squeeze(AP_sem(mua_var_exp_log(:,curr_str),1)),squeeze(AP_sem(mua_var_exp_log(:,curr_str),1)), ...
        'color','k','linewidth',2);
       
    scatter(nanmean(mua_var_exp_log(:,curr_str),1), ...
        nanmean(ctxpred_r2(:,curr_str),1),80, ...
        str_col(curr_str,:),'filled','MarkerEdgeColor','k','linewidth',2);
    scatter(mua_var_exp_log(:,curr_str),ctxpred_r2(:,curr_str),10, ...
        str_col(curr_str,:),'filled');
end
xlabel('log_{10}(striatum variance)');
ylabel('Cortex R^2');
legend({'DMS','DCS','DLS'})

% (variance vs. explained variance stats)
[exp_grp,domain_grp] = ndgrid(1:size(mua_var_exp_log,1),1:size(mua_var_exp_log,2));
use_points = ~isnan(mua_var_exp_log);
[h,atab,ctab,stats] = aoctool(mua_var_exp_log(use_points), ...
    ctxpred_r2(use_points),domain_grp(use_points));



%% ~~~~~~~~ BELOW: integrated in to AP_ctx_str_figures_v6
 
%% Plot Cortex R2 vs variance (low DMS expl var is due to low var)

% Load data
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\str_ctxpred.mat')

% Concatenate data across cohorts and animals
mua_animal = [str_ctxpred.str]';
mua_exp = vertcat(mua_animal{:});

mua_ctxpred_animal = [str_ctxpred.str_ctxpred]';
mua_ctxpred_exp = vertcat(mua_ctxpred_animal{:});

ctxpred_r2 = cellfun(@(mua,mua_ctxpred) ...
    1 - (nansum((mua-mua_ctxpred).^2,2)./(nansum((mua-nanmean(mua,2)).^2,2))), ...
    mua_exp,mua_ctxpred_exp,'uni',false);

ctxpred_r2_task = horzcat(ctxpred_r2{:,1})';
ctxpred_r2_passive = horzcat(ctxpred_r2{:,2})';

figure;
n_depths = size(ctxpred_r2_task,2);
str_col = max(hsv(n_depths)-0.2,0);


% Get and plot striatal variance by task vs. passive
mua_var_exp = cellfun(@(x) var(x,[],2),mua_exp,'uni',false);
mua_var_task = log10(horzcat(mua_var_exp{:,1})');
mua_var_passive = log10(horzcat(mua_var_exp{:,2})');


% Plot cortex R2 vs striatal variance
figure; hold on
for curr_str = 1:n_depths
    errorbar(squeeze(nanmean(mua_var_task(:,curr_str),1)), ...
        squeeze(nanmean(ctxpred_r2_task(:,curr_str),1)), ...
        squeeze(AP_sem(ctxpred_r2_task(:,curr_str),1)), ...
        squeeze(AP_sem(ctxpred_r2_task(:,curr_str),1)), ...
        squeeze(AP_sem(mua_var_task(:,curr_str),1)), ...
        squeeze(AP_sem(mua_var_task(:,curr_str),1)), ...
        'color','k','linewidth',2);
    
    scatter(mua_var_task(:,curr_str),ctxpred_r2_task(:,curr_str),10, ...
        str_col(curr_str,:),'filled');
    scatter(nanmean(mua_var_task(:,curr_str),1), ...
        nanmean(ctxpred_r2_task(:,curr_str),1),80, ...
        str_col(curr_str,:),'filled','MarkerEdgeColor','k','linewidth',2);
end
xlabel('log(Task variance)');
ylabel('Cortex task R^2');
legend({'DMS','DCS','DLS'})

% (variance vs. explained variance stats)
[exp_grp,domain_grp] = ndgrid(1:size(mua_var_task,1),1:size(mua_var_task,2));
use_points = ~isnan(mua_var_task);
[h,atab,ctab,stats] = aoctool(mua_var_task(use_points), ...
    ctxpred_r2_task(use_points),domain_grp(use_points));


%% Naive cortical kernels: compare to trained

% Load Master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Load trained data
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\str_ctxpred.mat')
% Concatenate kernels and convert to pixels
% (and flip in time so it's fluorescence lead:lag spikes)
ctx_str_k_trained_animal = [str_ctxpred.ctx_str_k]';
ctx_str_k_trained_animal = cellfun(@(x) cellfun(@(x) ...
    flip(AP_svdFrameReconstruct(U_master(:,:,1:100),x),3),x,'uni',false), ...
    ctx_str_k_trained_animal,'uni',false);

ctx_str_k_px_trained_cat = vertcat(ctx_str_k_trained_animal{:}); 
ctx_str_k_px_notask_trained_mean = nanmean(cat(5,ctx_str_k_px_trained_cat{:,2}),5);

% Load naive data
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\str_ctxpred_naive.mat')
% Concatenate kernels and convert to pixels
% (and flip in time so it's fluorescence lead:lag spikes)
ctx_str_k_naive_animal = [str_ctxpred.ctx_str_k]';
ctx_str_k_naive_animal = cellfun(@(x) cellfun(@(x) ...
    flip(AP_svdFrameReconstruct(U_master(:,:,1:100),x),3),x,'uni',false), ...
    ctx_str_k_naive_animal,'uni',false);

ctx_str_k_px_naive_cat = vertcat(ctx_str_k_naive_animal{:}); 
ctx_str_k_px_notask_naive_mean = nanmean(cat(5,ctx_str_k_px_naive_cat{:,1}),5);


% Get time
framerate = 35;
upsample_factor = 1;
sample_rate = framerate*upsample_factor;
kernel_t = [-0.1,0.1];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
t = kernel_frames/sample_rate;

% Get mean kernels and plot
n_depths = size(ctx_str_k_px_notask_trained_mean,4);

AP_imscroll([ctx_str_k_px_notask_trained_mean,ctx_str_k_px_notask_naive_mean]);
axis image;
colormap(brewermap([],'PRGn'));
caxis([-max(abs(caxis)),max(abs(caxis))]);
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5],[],[size(U_master,1),size(U_master,2),1,2]);

figure;
for curr_depth = 1:n_depths
    subplot(n_depths,2,(curr_depth-1)*2+1);
    imagesc(ctx_str_k_px_notask_trained_mean(:,:,t == 0,curr_depth));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'PRGn'));
    axis image off;
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title('Task');
    
    subplot(n_depths,2,(curr_depth-1)*2+2);
    imagesc(ctx_str_k_px_notask_naive_mean(:,:,t == 0,curr_depth));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'PRGn'));
    axis image off;
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title('Passive');
end

figure;
for curr_depth = 1:n_depths
    subplot(n_depths,1,curr_depth);
    imagesc(reshape(ctx_str_k_px_notask_trained_mean(:,:,:,curr_depth), ...
        size(ctx_str_k_px_notask_trained_mean,1),[]));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'PRGn'));
    axis image off;
    title('Task');
end

figure;
for curr_depth = 1:n_depths
    subplot(n_depths,1,curr_depth);
    imagesc(reshape(ctx_str_k_px_notask_naive_mean(:,:,:,curr_depth), ...
        size(ctx_str_k_px_notask_naive_mean,1),[]));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'PRGn'));
    axis image off;
    title('Passive');
end


% Correlate trained/passive kernels within/across domain
ctx_str_k_trained_animalmean = cellfun(@(x) ...
    reshape(nanmean(cat(5,x{:,2}),5),[],n_depths),ctx_str_k_trained_animal,'uni',false);

ctx_str_k_naive_animalmean = cellfun(@(x) ...
    reshape(nanmean(cat(5,x{:,1}),5),[],n_depths),ctx_str_k_naive_animal,'uni',false);

ctx_str_k_corr = corr([ ...
    horzcat(ctx_str_k_trained_animalmean{:}), ...
    horzcat(ctx_str_k_naive_animalmean{:})],'type','Pearson');

ctx_str_k_corr_tril = tril(ctx_str_k_corr,-1);
ctx_str_k_corr_tril(triu(true(size(ctx_str_k_corr)))) = NaN;

corr_grp = ...
    [repmat(1:n_depths,1, ...
    length(ctx_str_k_trained_animalmean)+ ...
    length(ctx_str_k_naive_animalmean)); ...
    1*ones(1,length(ctx_str_k_trained_animalmean)*n_depths), ...
    2*ones(1,length(ctx_str_k_naive_animalmean)*n_depths)];

k_corr_trained_withindomain = ctx_str_k_corr_tril(...
    corr_grp(2,:) == 1 & corr_grp(2,:)' == 1 & ...
    corr_grp(1,:) == corr_grp(1,:)');

k_corr_trained_acrossdomain = ctx_str_k_corr_tril(...
    corr_grp(2,:) == 1 & corr_grp(2,:)' == 1 & ...
    corr_grp(1,:) ~= corr_grp(1,:)');

k_corr_naive_withindomain = ctx_str_k_corr_tril(...
    corr_grp(2,:) == 2 & corr_grp(2,:)' == 2 & ...
    corr_grp(1,:) == corr_grp(1,:)');

k_corr_naive_acrossdomain = ctx_str_k_corr_tril(...
    corr_grp(2,:) == 2 & corr_grp(2,:)' == 2 & ...
    corr_grp(1,:) ~= corr_grp(1,:)');

k_corr_trainednaive_withindomain = ctx_str_k_corr_tril(...
    corr_grp(2,:) == 1 & corr_grp(2,:)' == 2 & ...
    corr_grp(1,:) == corr_grp(1,:)');

k_corr_trainednaive_acrossdomain = ctx_str_k_corr_tril(...
    corr_grp(2,:) == 1 & corr_grp(2,:)' == 2 & ...
    corr_grp(1,:) ~= corr_grp(1,:)');

k_corr_grp = {k_corr_trained_withindomain, ...
    k_corr_trainednaive_withindomain, ...
    k_corr_trainednaive_acrossdomain};

figure; hold on;
distributionPlot(k_corr_grp,'showMM',0,'color',[0.5,0.5,0.5])
plot(cellfun(@nanmean,k_corr_grp),'k','linewidth',2);
ylabel('Kernel spatiotemporal correlation');
set(gca,'XTickLabel',{'Trained within domain', ...
    'Trained-naive within domain', ...
    'Trained-naive across domain'}, ...
    'XTickLabelRotation',45);

% (cross task vs no task statistics)
curr_p = ranksum(k_corr_grp{1},k_corr_grp{2});
disp(['Trained within domain vs. trained-naive within domain p = ' ...
    num2str(curr_p)])
curr_p = ranksum(k_corr_grp{1},k_corr_grp{3});
disp(['Trained within domain vs. trained-naive across domain p = ' ...
    num2str(curr_p)])


%% Depth fluor/MUA correlation using depth-specific deconv kernel

%%% Load and plot depth-specific deconv kernel
deconv_k_depth_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\gcamp6s_kernel_ctxdepth.mat';
load(deconv_k_depth_fn)

gcamp6s_kernel_cat = fliplr(cat(3,gcamp6s_kernel_ctxdepth.regression{:}));
gcamp6s_kernel_norm = gcamp6s_kernel_cat./max(abs(gcamp6s_kernel_cat),[],2);
gcamp6s_kernel_mean = nanmean(gcamp6s_kernel_norm,3);

depth_col = lines(size(gcamp6s_kernel_mean,1));

figure; 

subplot(1,2,1); hold on;
p = AP_errorfill(gcamp6s_kernel_ctxdepth.regression_t, ...
    permute(nanmean(gcamp6s_kernel_norm,3),[2,1,3]), ...
    permute(nanstd(gcamp6s_kernel_norm,[],3),[2,1,3]), ...
    depth_col);
line(xlim,[0,0],'color','k','linestyle','--');
line([0,0],ylim,'color','k','linestyle','--');
set(p(2),'linestyle','--');

% (depth kernel statistics)
disp('Kernel 2-way ANOVA:')
[depth_grp,t_grp,exp_grp] = ndgrid(1:size(gcamp6s_kernel_norm,1), ...
    1:size(gcamp6s_kernel_norm,2),1:size(gcamp6s_kernel_norm,3));
curr_p = anovan(gcamp6s_kernel_norm(:),[depth_grp(:),t_grp(:)],'display','off');
disp(['depth p = ' num2str(curr_p(1))]);
disp(['time p = ' num2str(curr_p(2))]);

%%% Load and plot correlation data
use_protocol = 'vanillaChoiceworld';

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
data_fn = [data_path filesep 'ctx_fluor_mua_corr_' use_protocol '_ctxdepthkernel'];
load(data_fn);

mua_depth = data(1).cortex_mua_depth{1}; % they're all the same, use 1st
cortex_fluor_corr_cat = cell2mat(permute(horzcat( ...
    data.cortex_fluor_corr),[1,3,2]));

subplot(1,2,2); hold on;
p = AP_errorfill(mua_depth, ...
    nanmean(cortex_fluor_corr_cat,3), ...
    AP_sem(cortex_fluor_corr_cat,3), ...
    depth_col);
xlim([0,1400]);
xlabel('Cortical MUA aligned depth');
ylabel('Correlation');
legend(p,{'Superficial kernel','Deep kernel'});
set(p(2),'linestyle','--');

% (depth kernel correlation statistics)
disp('Depth correlation 2-way ANOVA:')
[depth_grp,kernel_grp,exp_grp] = ndgrid(1:size(cortex_fluor_corr_cat,1), ...
    1:size(cortex_fluor_corr_cat,2),1:size(cortex_fluor_corr_cat,3));
curr_p = anovan(cortex_fluor_corr_cat(:),[depth_grp(:),kernel_grp(:)],'display','off');
disp(['depth p = ' num2str(curr_p(1))]);
disp(['kernel p = ' num2str(curr_p(2))]);







