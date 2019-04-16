% Generate figures for ctx-str paper

% Anything that takes a lot of time is done in
% AP_ctx_str_trial_preprocessing and saved for plotting here

% The original scripts here were in test_wf_ephys_choiceworld_analysis

% (this is addendum to AP_ctx_str_figures_v2)

%% Figure 1: Example raster plot

animal = 'AP025'; 
day = '2017-10-01'; 
experiment = 1; 
str_align = 'none'; 
verbose = true; 
AP_load_experiment;

% Load task->unit regression
unit_kernel_dir = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
unit_kernel_fn = 'unit_kernel_all_triaged.mat';
load([unit_kernel_dir filesep unit_kernel_fn]);
regressor_labels = {'Stim','Move','Go cue','Outcome'};
regressor_cols = [01,0,0;1,0,1;0.8,0.7,0.2;0,0,1];

curr_animal = 2;
curr_day = 3;
use_partial = 2;
[~,max_regressor_idx] = max(unit_kernel_all(curr_animal,curr_day).unit_expl_var_partial(:,:,use_partial),[],2);



figure; 

% Plot wheel velocity
wheel_axes = subplot(6,1,1);
plot(wheel_axes,Timeline.rawDAQTimestamps,wheel_velocity,'k','linewidth',2);
ylabel('Wheel velocity');

% Plot raster
raster_axes = subplot(6,1,2:6,'YDir','reverse'); hold on;
% (time to plot spikes)
raster_plot_t = [0,60*10]; % time from start of experiment
plot_spikes = spike_times_timeline >= raster_plot_t(1) & ...
    spike_times_timeline <= raster_plot_t(2) & ...
    spike_depths >= str_depth(1) & spike_depths <= str_depth(2);
% plot(raster_axes,spike_times_timeline(plot_spikes),spike_depths(plot_spikes),'.k');
scatter(raster_axes,spike_times_timeline(plot_spikes),spike_depths(plot_spikes), ...
    5,regressor_cols(max_regressor_idx(spike_templates(plot_spikes)),:),'filled');
ylabel('Depth (\mum)');
xlabel('Time (s)')

linkaxes([wheel_axes,raster_axes],'x');

% Plot all stimuli
stim_col = colormap_BlueWhiteRed(5);
[~,trial_contrast_idx] = ...
    ismember(trial_conditions(:,1).*trial_conditions(:,2),unique(contrasts'.*sides),'rows');
stim_lines = arrayfun(@(x) line(raster_axes,repmat(stimOn_times(x),1,2),ylim,'color', ...
    stim_col(trial_contrast_idx(x),:),'linewidth',2),1:length(stimOn_times));

% Plot movement starts
move_col = [0.6,0,0.6;0,0.6,0];
[~,trial_choice_idx] = ismember(trial_conditions(:,3),[-1;1],'rows');
move_lines = arrayfun(@(x) line(raster_axes,repmat(wheel_move_time(x),1,2),ylim,'color', ...
    move_col(trial_choice_idx(x),:),'linewidth',2),find(~isnan(wheel_move_time))');

% Plot go cues
go_col = [0.8,0.8,0.2];
go_cue_times = signals_events.interactiveOnTimes(1:n_trials);
go_cue_lines = arrayfun(@(x) line(raster_axes,repmat(go_cue_times(x),1,2),ylim,'color', ...
    go_col,'linewidth',2,'linestyle','--'),1:length(go_cue_times));

% Plot outcomes
outcome_col = [0,0,0.8;0.5,0.5,0.5];
reward_lines = arrayfun(@(x) line(raster_axes,repmat(reward_t_timeline(x),1,2),ylim,'color', ...
    outcome_col(1,:),'linewidth',2,'linestyle','--'),1:length(reward_t_timeline));
punish_times = signals_events.responseTimes(trial_outcome == -1);
punish_lines = arrayfun(@(x) line(raster_axes,repmat(punish_times(x),1,2),ylim,'color', ...
    outcome_col(2,:),'linewidth',2,'linestyle','--'),1:length(punish_times));

% Set x limits
xlim(raster_axes,[134,154]);

% Write legend
[~,unique_contrasts_h] = unique(trial_contrast_idx);
[~,unique_move_h] = unique(trial_choice_idx(trial_choice_idx > 0));
legend([stim_lines(unique_contrasts_h),move_lines(unique_move_h), ...
    go_cue_lines(1),reward_lines(1),punish_lines(1)], ...
    [cellfun(@(x) ['Stim ' num2str(x)],num2cell(unique(contrasts'.*sides)),'uni',false); ...
    {'Move L';'Move R';'Go cue';'Reward';'Punish'}]);




%% Figure 1: Task -> striatum unit regression

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

% Plot example recording
curr_animal = 2;
curr_day = 3;

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
    
    p1 = subplot(1,length(regressor_labels)+3,curr_regressor+2,'YDir','reverse');
    hold on;
    curr_expl_var = cellfun(@(x) x(:,curr_regressor,use_partial),unit_expl_var_partial_cat,'uni',false);
  
    % Normalize across experiments
    curr_expl_var_cat = cell2mat(curr_expl_var);
    curr_expl_var_norm = curr_expl_var_cat./ ...
        max(curr_expl_var_cat(curr_expl_var_cat < 1));
    curr_expl_var_dotsize = 100*rescale(curr_expl_var_norm.*(curr_expl_var_norm > 0),0,1) + 1;
    norm_type = 'concat normalized';    
    
    scatter(log10(cell2mat(spike_rate_cat)), ...
        cell2mat(template_depths_cat), ...
        curr_expl_var_dotsize, ...
        regressor_cols(curr_regressor,:),'filled');
    title({regressor_labels{curr_regressor},norm_type});
    xlabel('log10(spike rate)');
    ylabel('Depth from str end (\mum)');
    line(xlim,[0,0],'color','k','linewidth',2);
    
    p2 = subplot(1,length(regressor_labels)+3,length(regressor_labels)+3,'YDir','reverse');
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

% (draw striatum starts)
subplot(1,length(regressor_labels)+3,length(regressor_labels)+3)
str_start_line = nan(size(str_depth_cat,1),1);
for i = 1:size(str_depth_cat,1)
    str_start_line(i) = line(xlim,repmat(-diff(str_depth_cat(i,:)),1,2),'color',[0.8,0.8,0.8]);
end
set(gca,'Children',circshift(get(gca,'Children'),-size(str_depth_cat,1)));

linkaxes(get(h,'Children'),'y');


%% Figure 1: Striatal unit regression examples

% (animal 2 day 3, AP025 2017-10-01, run unit regression first)

smooth_size = 11;
gw = gausswin(smooth_size,5)';
smWin = gw./sum(gw);
% smWin = 1;

binned_spikes_smoothed = conv2(binned_spikes,smWin,'same');
binned_spikes_taskpred_smoothed = conv2(binned_spikes_taskpred,smWin,'same');

figure;
example_units = [237,150,170,66];
p = nan(length(example_units),1);
for curr_unit_idx = 1:length(example_units)
    curr_unit = example_units(curr_unit_idx);
    
    p(curr_unit_idx) = subplot(length(example_units),1,curr_unit_idx);
    hold on;
    plot(binned_spikes_smoothed(curr_unit,:),'k');
    plot(binned_spikes_taskpred_smoothed(curr_unit,:),'r');
    title(curr_unit);
    xlabel('Time');
    ylabel('Spikes');
end

linkaxes(p,'x');

x_bounds = [10877,12010];
xlim(p(1),x_bounds);


%% Figure 3: Compare kernels with domain maps

% Load data
data_fn = 'trial_activity_choiceworld';
exclude_data = true;
AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(mua_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_animal;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));



%%% Task>cortex regression results

% Get task>cortex parameters
n_regressors = 4;
t_shifts = {[0,0.5]; ... % stim
    [-0.5,1]; ... % move
    [-0.1,0.5]; ... % go cue
    [-0.5,1]}; % outcome
regressor_labels = {'Stim','Move onset','Go cue','Outcome'};

sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
    round(x(2)*(sample_rate)),t_shifts,'uni',false);
t_shifts = cellfun(@(x) x/sample_rate,sample_shifts,'uni',false);

% Get average task>striatum kernels
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

% Plot max for each kernel across time
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





% Get average cortex->striatum kernel
ctx_str_k_mean = nanmean(cell2mat(permute(vertcat(ctx_str_k_all{:}),[2,3,4,5,1])),5);
ctx_str_k_mean_px = cell2mat(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:50), ...
    ctx_str_k_mean(:,:,x)),permute(1:n_depths,[1,3,4,2]),'uni',false));

ctx_str_kernel_frames_t = [-0.5,0.5];
ctx_str_kernel_frames = round(ctx_str_kernel_frames_t(1)*sample_rate): ...
    round(ctx_str_kernel_frames_t(2)*sample_rate);
ctx_str_kernel_t = ctx_str_kernel_frames./sample_rate;

% (use t = 0)
ctx_str_kernel_t0 = squeeze(ctx_str_k_mean_px(:,:,ctx_str_kernel_t == 0,:));

% Get correlation between task->cortex and cortex->striatum regressor
ctx_task_str_corr = cellfun(@(k) cell2mat(arrayfun(@(x) 1 - pdist2(reshape(k(:,:,x),[],1)', ...
    reshape(ctx_str_kernel_t0,[],n_depths)','correlation'), ...
    transpose(1:size(k,3)),'uni',false)),regressor_t_max,'uni',false); 

figure;
for curr_regressor = 1:n_regressors
    if curr_regressor == 1
        col = colormap_BlueWhiteRed(size(ctx_task_str_corr{curr_regressor},1)/2);
        col(median(1:size(col,1)),:) = [];
    else
        col = lines(size(ctx_task_str_corr{curr_regressor}));
    end
    
    subplot(n_regressors,1,curr_regressor); hold on;
    set(gca,'ColorOrder',col);
    bar(ctx_task_str_corr{curr_regressor}');    
    set(gca,'XTick',1:4);
    xlabel('Striatum domain')
    ylabel('Kernel correlation')
    title(regressor_labels{curr_regressor});
end









