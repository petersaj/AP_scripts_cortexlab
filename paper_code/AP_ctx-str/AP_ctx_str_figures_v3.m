% Generate figures for ctx-str paper 
% CURRENTLY BUILDING FOR COMPREHENSIVE v10
% (trials data is prepared in AP_ctx_str_trial_preprocessing)

%% Load in task data

% Load data
% data_fn = 'trial_activity_choiceworld'; % Primary dataset
data_fn = 'trial_activity_choiceworld_strdepth'; % Depth-aligned striatum
% data_fn = 'trial_activity_choiceworld_stimxmove'; % Including stim x move regressors (no real difference)
% data_fn = 'trial_activity_choiceworld_wfonly'; % Widefield-only days (no craniotomy, so cleaner)
% data_fn = 'trial_activity_choiceworld_oldregressors'; % Longer timeshifts and go early/late and outcome correct/incorrect
exclude_data = true;

% data_fn = 'trial_activity_AP_choiceWorldStimPassive_trained';
% data_fn = 'trial_activity_AP_choiceWorldStimPassive_naive';
% data_fn = 'trial_activity_stimKalatsky_naive';
% data_fn = 'trial_activity_stimKalatsky_trained';
% exclude_data = false;

AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_animal;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));


%% Fig 1a: Behavior psychometric

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
    exclude_fn = 'bhv_use_experiments';
    load([exclude_path filesep exclude_fn]);
    
    % Pull out used experiments for the animals loaded
    use_experiments_animals = ismember({bhv.animal},{bhv_use_experiments.animals});
    use_experiments = {bhv_use_experiments(use_experiments_animals).use_experiments}';
    
    % Cut out bad experiments for any experiment data fields
    for curr_field = bhv_fieldnames(experiment_fields)'
        for curr_animal = 1:length(use_experiments)
            bhv(curr_animal).(cell2mat(curr_field)) = ...
                bhv(curr_animal).(cell2mat(curr_field))(use_experiments{curr_animal});
        end
    end
end

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

figure; hold on; axis square;
plot(conditions,frac_left,'linewidth',2,'color',[0.5,0.5,0.5]);
plot(conditions,nanmean(frac_left,2),'linewidth',5,'color','k');
xlim([-1,1]);
ylim([0,1]);
line([0,0],ylim,'linestyle','--','color','k');
line(xlim,[0.5,0.5],'linestyle','--','color','k');
xlabel('Stimulus side*contrast');
ylabel('Fraction go left');


%% Fig 1b: Example recording

animal = 'AP025'; 
day = '2017-10-01'; 
experiment = 1; 
str_align = 'none'; 
verbose = true; 
AP_load_experiment;

%%% Plot example data

% Align U's, deconvolve widefield
use_components = 1:200;
aUdf = AP_align_widefield(animal,day,Udf);
fVdf_deconv = AP_deconv_wf(fVdf);

% Set time to plot
plot_t = [134,152];

raster_fig = figure;

% (wheel velocity)
wheel_axes = subplot(6,1,6);
plot_wheel_idx = Timeline.rawDAQTimestamps >= plot_t(1) & ...
    Timeline.rawDAQTimestamps <= plot_t(2);
plot(wheel_axes,Timeline.rawDAQTimestamps(plot_wheel_idx), ...
    wheel_velocity(plot_wheel_idx),'k','linewidth',2);
ylabel('Wheel velocity');

% (stimuli)
stim_col = colormap_BlueWhiteRed(5);
[~,trial_contrast_idx] = ...
    ismember(trial_conditions(:,1).*trial_conditions(:,2),unique(contrasts'.*sides),'rows');
stim_lines = arrayfun(@(x) line(wheel_axes,repmat(stimOn_times(x),1,2),ylim(wheel_axes),'color', ...
    stim_col(trial_contrast_idx(x),:),'linewidth',2), ...
    find(stimOn_times >= plot_t(1) & stimOn_times <= plot_t(2)));

% (movement starts)
move_col = [0.6,0,0.6;0,0.6,0];
[~,trial_choice_idx] = ismember(trial_conditions(:,3),[-1;1],'rows');
move_lines = arrayfun(@(x) line(wheel_axes,repmat(wheel_move_time(x),1,2),ylim(wheel_axes),'color', ...
    move_col(trial_choice_idx(x),:),'linewidth',2), ...
    find(wheel_move_time >= plot_t(1) & wheel_move_time <= plot_t(2)));

% (go cues)
go_col = [0.8,0.8,0.2];
go_cue_times = signals_events.interactiveOnTimes(1:n_trials);
go_cue_lines = arrayfun(@(x) line(wheel_axes,repmat(go_cue_times(x),1,2),ylim(wheel_axes),'color', ...
    go_col,'linewidth',2,'linestyle','--'), ...
    find(go_cue_times >= plot_t(1) & go_cue_times <= plot_t(2)));

% (outcomes)
outcome_col = [0,0,0.8;0.5,0.5,0.5];
reward_lines = arrayfun(@(x) line(wheel_axes,repmat(reward_t_timeline(x),1,2),ylim(wheel_axes),'color', ...
    outcome_col(1,:),'linewidth',2,'linestyle','--'), ...
    find(reward_t_timeline >= plot_t(1) & reward_t_timeline <= plot_t(2)));
punish_times = signals_events.responseTimes(trial_outcome == -1);
punish_lines = arrayfun(@(x) line(wheel_axes,repmat(punish_times(x),1,2),ylim(wheel_axes),'color', ...
    outcome_col(2,:),'linewidth',2,'linestyle','--'), ...
    find(punish_times >= plot_t(1) & punish_times <= plot_t(2)));

% (striatum raster)
raster_axes = subplot(6,1,3:5,'YDir','reverse'); hold on;
plot_spikes = spike_times_timeline >= plot_t(1) & ...
    spike_times_timeline <= plot_t(2) & ...
    spike_depths >= str_depth(1) & spike_depths <= str_depth(2);
plot(raster_axes,spike_times_timeline(plot_spikes),spike_depths(plot_spikes),'.k');
ylabel('Depth (\mum)');
xlabel('Time (s)')

% (fluorescence from select ROIs)
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
roi_trace = AP_svd_roi(aUdf(:,:,use_components),fVdf_deconv(use_components,:),[],[],cat(3,wf_roi.mask));

plot_rois = [1,7,9];
fluor_spacing = 70;
fluor_axes = subplot(6,1,1:2); hold on;
plot_fluor_idx = frame_t >= plot_t(1) & frame_t <= plot_t(2);
AP_stackplot(roi_trace(plot_rois,plot_fluor_idx)', ...
    frame_t(plot_fluor_idx),fluor_spacing,false,[0,0.7,0],{wf_roi(plot_rois).area});

linkaxes([wheel_axes,raster_axes,fluor_axes],'x');

% % Write legend
% [~,unique_contrasts_h] = unique(trial_contrast_idx);
% [~,unique_move_h] = unique(trial_choice_idx(trial_choice_idx > 0));
% legend([stim_lines(unique_contrasts_h),move_lines(unique_move_h), ...
%     go_cue_lines(1),reward_lines(1),punish_lines(1)], ...
%     [cellfun(@(x) ['Stim ' num2str(x)],num2cell(unique(contrasts'.*sides)),'uni',false); ...
%     {'Move L';'Move R';'Go cue';'Reward';'Punish'}]);

% Plot ROIs
figure; hold on
set(gca,'YDir','reverse');
AP_reference_outline('ccf_aligned','k');
for curr_roi = plot_rois
    curr_roi_boundary = cell2mat(bwboundaries(wf_roi(curr_roi).mask()));
    patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),[0,0.8,0]);   
    text(nanmean(curr_roi_boundary(:,2)),nanmean(curr_roi_boundary(:,1)), ...
        wf_roi(curr_roi).area,'FontSize',12,'HorizontalAlignment','center')
end
axis image off;

%% Fig 1x: Average stim-aligned cortex

% Get average stim-aligned fluorescence 
plot_trials = trial_contrast_allcat > 0 & trial_side_allcat == 1 & trial_choice_allcat == -1;
plot_trials_exp = mat2cell(plot_trials,use_split,1);

fluor_allcat_deconv_exp = mat2cell(fluor_allcat_deconv,use_split,length(t),n_vs);
plot_px = nanmean(cell2mat(permute(cellfun(@(x,trials) svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    squeeze(nanmean(x(trials,:,:),1))'),fluor_allcat_deconv_exp,plot_trials_exp,'uni',false),[2,3,4,1])),4);

plot_t = [find(t > 0.07,1),find(t > 0.18,1),find(t > 0.3,1),find(t > 0.65,1)];
figure;
for curr_t = 1:length(plot_t)
   subplot(1,length(plot_t),curr_t);
   imagesc(plot_px(:,:,plot_t(curr_t)));
   AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
   axis image off;
   colormap(brewermap([],'Greens'));
   caxis([0,0.01]);
   title([num2str(t(plot_t(curr_t))) 's from stim']);
end

%% Fig 1c: Task > cortex kernels

% Get task>cortex parameters
n_regressors = length(task_regressor_labels);
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

% Get average task > cortex kernels (pixels and ROIs)
regressor_px = cell(n_regressors,1);
regressor_roi = cell(n_regressors,1);
for curr_regressor = 1:n_regressors
    curr_k_cell = cellfun(@(x) x(curr_regressor,:),vertcat(fluor_taskpred_k_all{:}),'uni',false);
    curr_k_cell = vertcat(curr_k_cell{:});
    curr_k = nanmean(cat(4,curr_k_cell{:}),4); 
    curr_k_px = cell2mat(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:n_vs), ...
        permute(curr_k(x,:,:),[3,2,1])),permute(1:size(curr_k,1),[1,3,4,2]),'uni',false));   
    AP_image_scroll(curr_k_px,task_regressor_t_shifts{curr_regressor});
    axis image;
    caxis([-max(caxis),max(caxis)])
    colormap(crameri('cork'));
    AP_reference_outline('ccf_aligned','k');
    set(gcf,'Name',task_regressor_labels{curr_regressor});
    
    regressor_px{curr_regressor} = curr_k_px;
    regressor_roi{curr_regressor} = permute(cell2mat(arrayfun(@(x) AP_svd_roi(U_master(:,:,1:n_vs), ...
        permute(curr_k(x,:,:),[3,2,1]),[],[],cat(3,wf_roi.mask)), ...
        permute(1:size(curr_k,1),[1,3,2]),'uni',false)),[3,2,1]);
end

% Plot max for each kernel group
regressor_t_max = cellfun(@(x) squeeze(max(x,[],3)),regressor_px,'uni',false);

figure; hold on;
max_c = max(abs(cell2mat(cellfun(@(x) x(:),regressor_t_max,'uni',false))));
for curr_regressor = 1:n_regressors
    subplot(1,n_regressors,curr_regressor)
    imagesc(max(regressor_t_max{curr_regressor},[],3));
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    axis image off;
    colormap(brewermap([],'Greens'));
    caxis([0,max_c]);
    title(task_regressor_labels{curr_regressor});
end

% Plot max for each kernel separately
figure;
max_subregressors = max(cellfun(@(x) size(x,3),regressor_t_max));
max_c = max(abs(cell2mat(cellfun(@(x) x(:),regressor_t_max,'uni',false))));
for curr_regressor = 1:n_regressors
    for curr_subregressor = 1:size(regressor_t_max{curr_regressor},3)
        subplot(n_regressors,max_subregressors, ...
            curr_subregressor+(curr_regressor-1)*max_subregressors);
        imagesc(regressor_t_max{curr_regressor}(:,:,curr_subregressor));
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        axis image off; 
        colormap(brewermap([],'Greens'));
        caxis([0,max_c]);
    end
end

% Get ROI traces for each subregressor
plot_rois = [1,7,9,10];

stim_col = colormap_BlueWhiteRed(5);
move_col = [0.6,0,0.6;0,0.6,0];
go_col = [0.8,0.8,0.2];
outcome_col = [0,0,0.8];
task_regressor_cols = {stim_col,move_col,go_col,outcome_col};
figure;
p = nan(length(plot_rois),n_regressors);
for curr_roi_idx = 1:length(plot_rois)
    curr_roi = plot_rois(curr_roi_idx);    
    for curr_regressor = 1:n_regressors
        p(curr_roi_idx,curr_regressor) = ...
            subplot(length(plot_rois),n_regressors, ...
            sub2ind([n_regressors,length(plot_rois)],curr_regressor,curr_roi_idx)); hold on;
        set(gca,'ColorOrder',task_regressor_cols{curr_regressor});
        plot(task_regressor_t_shifts{curr_regressor}, ...
            regressor_roi{curr_regressor}(:,:,curr_roi)','linewidth',2);
        ylabel([wf_roi(curr_roi).area ' \DeltaF/F_0']);
        xlabel('Time (s)');
        line([0,0],ylim,'color','k');
    end
    linkaxes(p(curr_roi_idx,:),'y');
end


%% Fig 1e: Task > striatum explained variance histogram

% Load the unit kernel results
unit_kernel_dir = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
unit_kernel_fn = 'unit_kernel_all.mat';
load([unit_kernel_dir filesep unit_kernel_fn]);
task_regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
task_regressor_cols = [1,0,0;1,0,1;0.8,0.7,0.2;0,0,1];

% Get estimation of end of striatum for each recording
ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
ephys_align_fn = ['ephys_depth_align.mat'];
load([ephys_align_path filesep ephys_align_fn]);

% Concatenate all units and re-depth by end of striatum

% (get estimation of end of striatum for each recording)
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

% Plot all units explained by kernels and expl var histogram by depth
use_partial = 2;
rate_cutoff = 0.1;
use_units_rate = cell2mat(spike_rate_cat) > rate_cutoff;

template_depths_allcat = cell2mat(template_depths_cat);

depth_bin_um = 400;
depth_bins = linspace(min(template_depths_allcat), ...
    max(template_depths_allcat),round(range(template_depths_allcat)/depth_bin_um));
depth_bin_centers = depth_bins(1:end-1) + diff(depth_bins)./2;
unit_depth_bins = discretize(template_depths_allcat,depth_bins);

unit_expl_var_partial_allcat = cell2mat(cellfun(@(x) x(:,:,use_partial), ...
    unit_expl_var_partial_cat,'uni',false));

figure; 
for curr_regressor = 1:length(task_regressor_labels)
    subplot(1,length(task_regressor_labels)+1,curr_regressor);
    hold on; set(gca,'YDir','reverse');
    plot_units = use_units_rate & unit_expl_var_partial_allcat(:,curr_regressor) > 0;
    scatter(unit_expl_var_partial_allcat(plot_units,curr_regressor), ...
        template_depths_allcat(plot_units),10,task_regressor_cols(curr_regressor,:),'filled');
    xlabel('Expl var');
    ylabel('Depth (\mum)');
    
    subplot(1,length(task_regressor_labels)+1,length(task_regressor_labels)+1);   
    hold on; set(gca,'YDir','reverse');
    curr_expl_var_norm = unit_expl_var_partial_allcat(:,curr_regressor);
    curr_expl_var_norm(curr_expl_var_norm < 0) = 0;
    curr_expl_var_norm_depth = accumarray(unit_depth_bins(plot_units), ...
        curr_expl_var_norm(plot_units,:), ...
        [length(depth_bin_centers),1],@nanmean,NaN);   
    plot(curr_expl_var_norm_depth,depth_bin_centers, ...
        'linewidth',2,'color',task_regressor_cols(curr_regressor,:));
    xlabel('Mean expl var (> 0)');
    ylabel('Depth (\mum)');
end

% (draw striatum starts)
for curr_plot = 1:length(task_regressor_labels)+1
    subplot(1,length(task_regressor_labels)+1,curr_plot)
    str_start_line = nan(size(str_depth_cat,1),1);
    for i = 1:size(str_depth_cat,1)
        str_start_line(i) = line(xlim,repmat(-diff(str_depth_cat(i,:)),1,2),'color',[0.8,0.8,0.8]);
    end
    set(gca,'Children',circshift(get(gca,'Children'),-size(str_depth_cat,1)));
end

%% Fig 2: BRING CODE OVER FROM V2


%% Fig 3a-b: Striatal domain activity and task regression

% Plot average stimulus-aligned activity in striatum
plot_trials = trial_contrast_allcat > 0 & trial_side_allcat == 1 & trial_choice_allcat == -1;
plot_trials_exp = mat2cell(plot_trials,use_split,1);

mua_allcat_exp = mat2cell(mua_allcat,use_split,length(t),n_depths);
mua_allcat_exp_mean = cell2mat(permute(cellfun(@(x,trials) ...
    permute(nanmean(x(trials,:,:),1),[3,2,1]),mua_allcat_exp,plot_trials_exp,'uni',false),[2,3,1]));

figure;
p = nan(n_depths,1);
for curr_depth = 1:n_depths
    p(curr_depth) = subplot(n_depths,1,curr_depth);
    AP_errorfill(t,nanmean(mua_allcat_exp_mean(curr_depth,:,:),3), ...
        AP_sem(mua_allcat_exp_mean(curr_depth,:,:),3),'k',0.5)
    xlabel('Time from stimulus');
    ylabel('Spikes (std)');
    line([0,0],ylim,'color','k');
    title(['Striatum depth ' num2str(curr_depth)]);
end
linkaxes(p);

% Get task>striatum parameters
n_regressors = length(task_regressor_labels);
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

% Normalize task > striatum kernels across experiments with mua_norm
mua_taskpred_k_all_norm = cellfun(@(kernel_animal,mua_norm_animal) ...
    cellfun(@(kernel_exp,mua_norm_exp) ...
    cellfun(@(kernel_regressor) ...
    kernel_regressor./(mua_norm_exp/sample_rate), ...
    kernel_exp,'uni',false),kernel_animal,mua_norm_animal,'uni',false), ...
    mua_taskpred_k_all,mua_norm,'uni',false);

% Average and concatenate task>striatum kernels within animals
task_str_k_animal = cell(n_regressors,1);
for curr_animal = 1:length(mua_taskpred_k_all_norm)
    if isempty(mua_taskpred_k_all_norm{curr_animal})
        continue
    end
    curr_k = cat(2,mua_taskpred_k_all_norm{curr_animal}{:});
    for curr_regressor = 1:n_regressors
        curr_k_mean = nanmean(cat(4,curr_k{curr_regressor,:}),4);        
        task_str_k_animal{curr_regressor} = cat(4, ...
            task_str_k_animal{curr_regressor},curr_k_mean);
    end
end

% Plot task>striatum kernels
stim_col = colormap_BlueWhiteRed(5);
move_col = [0.6,0,0.6;0,0.6,0];
go_col = [0.8,0.8,0.2];
outcome_col = [0,0,0.8];
task_regressor_cols = {stim_col,move_col,go_col,outcome_col};

figure;
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors  
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);  
        
        curr_kernels = permute(task_str_k_animal{curr_regressor}(:,:,curr_depth,:),[1,2,4,3]);
        n_subregressors = size(task_str_k_animal{curr_regressor},1);
        col = task_regressor_cols{curr_regressor};
        for curr_subregressor = 1:n_subregressors
            AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
                nanmean(curr_kernels(curr_subregressor,:,:),3), ...
                AP_sem(curr_kernels(curr_subregressor,:,:),3), ...
                col(curr_subregressor,:),0.5);
        end
        
        xlabel('Time (s)');
        ylabel('Weight');
        title(task_regressor_labels{curr_regressor});
        line([0,0],ylim,'color','k');
        
    end
end
linkaxes(p);


% Plot task>striatum regression examples
figure;
for curr_depth = 1:n_depths   
    
    % Set current data (pad trials with NaNs for spacing)
    n_pad = 10;
    curr_data = padarray(mua_allcat(:,:,curr_depth),[0,n_pad],NaN,'post');
    curr_pred_data = padarray(mua_taskpred_allcat(:,:,curr_depth),[0,n_pad],NaN,'post');
    nan_samples = isnan(curr_data) | isnan(curr_pred_data);

    % Smooth
    smooth_filt = ones(1,3)/3;
    curr_data = conv2(curr_data,smooth_filt,'same');
    curr_pred_data = conv2(curr_pred_data,smooth_filt,'same');
    
    % Set common NaNs for R^2    
    curr_data_nonan = curr_data; 
    curr_data_nonan(nan_samples) = NaN;
    
    curr_pred_data_nonan = curr_pred_data; 
    curr_pred_data(nan_samples) = NaN; 
    
    % Get squared error for each trial
    trial_r2 = 1 - (nansum((curr_data_nonan-curr_pred_data_nonan).^2,2)./ ...
        nansum((curr_data_nonan-nanmean(curr_data_nonan,2)).^2,2));
    
    [~,trial_r2_rank] = sort(trial_r2);
    
    plot_prctiles = round(prctile(1:length(trial_r2),linspace(25,75,20)));
    plot_trials = trial_r2_rank(plot_prctiles);
    
    subplot(n_depths,1,curr_depth); hold on;
    plot(reshape(curr_data(plot_trials,:)',[],1),'k');
    plot(reshape(curr_pred_data(plot_trials,:)',[],1),'b');
    
end

% Get R^2 for task regression 
taskpred_r2 = nan(max(split_idx),n_depths);
for curr_exp = 1:max(split_idx)
    curr_data = reshape(permute(mua_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_pred_data = reshape(permute(mua_taskpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_pred_data);
    curr_data(nan_samples) = NaN;
    curr_pred_data(nan_samples) = NaN;

    taskpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_pred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
end
figure;
errorbar(nanmean(taskpred_r2,1),AP_sem(taskpred_r2,1),'k','linewidth',2);
xlabel('Striatum depth');
ylabel('Task explained variance');


%% Fig 3x: ctx->str ?????????????????????




%% Fig 4x: Plot kernel ROIs

% Plot kernel ROIs
figure;
for curr_roi = 1:n_depths
    subplot(n_depths,1,curr_roi,'YDir','reverse');
    
    curr_roi_boundary = bwboundaries(kernel_roi(:,:,curr_roi));
    for curr_boundary = 1:length(curr_roi_boundary)
        patch(curr_roi_boundary{curr_boundary}(:,2), ...
            curr_roi_boundary{curr_boundary}(:,1),[0,0.8,0]);
    end
    
    AP_reference_outline('ccf_aligned','k');
    axis image off;
    title(['ROI for Str ' num2str(curr_roi)])
end


%% Fig 4x: Striatum v Cortex by condition

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


%% Fig 4x: Cortex-predicted striatum error

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


%% Fig 4x: Cortex-predicted striatum error (passive choiceworld stim: trained/naive)
% run separately for trained/naive then combine posthoc with copyobj

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

% Set activity percentiles and bins
act_prctile = [10,90];
n_act_bins = 5;

% Set areas and conditions
plot_areas = [1,2,3,4];

% Loop across area pairs, plot binned predicted v measured activity
curr_act_allcat = mua_allcat;
curr_act_pred_allcat = mua_ctxpred_allcat;

measured_v_pred_fig = figure('Name',data_fn,'color','w');
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
    set(gca,'XTick',1:2,'XTickLabel',{'Contra','Ipsi'});
    xlim([0.5,size(timeavg_trial_conditions{1},2)+0.5]);
    
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+2,curr_area_idx));
     hold on;
    errorbar( ...
        nanmean(curr_act_pred_condition_diff,1), ...
        AP_sem(curr_act_pred_condition_diff,1),'linewidth',2,'CapSize',0);
    ylabel('Condition difference');
    xlabel('Stim');
    set(gca,'XTick',1:2,'XTickLabel',{'Contra','Ipsi'});
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

%% Fig 4x: Cortex-predicted striatum error (passive fullscreen stim: trained/naive)
% run separately for trained/naive then combine posthoc with copyobj

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(mua_allcat,1),1);

% Set windows to average activity
timeavg_labels = {'Stim'};
timeavg_t = {[0.05,0.15]};
timeavg_align = {stim_align};
timeavg_trial_conditions = ...
    {[trial_stim_allcat == 1,trial_stim_allcat == 2,trial_stim_allcat == 3]};

% Set activity percentiles and bins
act_prctile = [10,90];
n_act_bins = 5;

% Set areas and conditions
plot_areas = [1,2,3,4];

% Loop across area pairs, plot binned predicted v measured activity
curr_act_allcat = mua_allcat;
curr_act_pred_allcat = mua_ctxpred_allcat;

measured_v_pred_fig = figure('Name',data_fn,'color','w');
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
            nanmean(act(trial_cond(:,2)) - pred(trial_cond(:,2))) - ...
            nanmean(act(trial_cond(:,3)) - pred(trial_cond(:,3))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp,use_trials);
        
        % Get condition-shuffled act-pred difference
        trial_conditions_exp_shuff = cellfun(@(x) AP_shake(repmat(x,1,1,n_shuff),1), ...
            trial_conditions_exp,'uni',false);   
                
        curr_act_pred_condition_diff_shuff(:,:,curr_timeavg) = ...
            cell2mat(arrayfun(@(shuff) ...
            cellfun(@(act,pred,trial_cond,use_trials) ...
            nanmean(act(trial_cond(:,2,shuff)) - pred(trial_cond(:,2,shuff))) - ...
            nanmean(act(trial_cond(:,3,shuff)) - pred(trial_cond(:,3,shuff))), ...
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
    set(gca,'XTick',1:3,'XTickLabel',{'Ipsi','Central','Contra'});
    xlim([0.5,size(timeavg_trial_conditions{1},2)+0.5]);
    
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+2,curr_area_idx));
     hold on;
    errorbar( ...
        nanmean(curr_act_pred_condition_diff,1), ...
        AP_sem(curr_act_pred_condition_diff,1),'linewidth',2,'CapSize',0);
    ylabel('Condition difference');
    xlabel('Stim');
    set(gca,'XTick',1:3,'XTickLabel',{'Ipsi','Central','Contra'});
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


%% ~~~~~~~~~~~~~ UNUSED ~~~~~~~~~~~

%% Spatial map of explained variance in cortex 

use_t = true(size(t));
spatial_downsample = 10;

correct_trials = trial_side_allcat == -trial_choice_allcat;
incorrect_trials = trial_side_allcat == trial_choice_allcat;

visual_trials = trial_contrast_allcat > 0;
zero_trials = trial_contrast_allcat == 0;

ctx_expl_var_correct = AP_spatial_explained_var(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv(correct_trials,use_t,:),[2,1,3]),[],n_vs)', ...
    reshape(permute(fluor_taskpred_allcat(correct_trials,use_t,:),[2,1,3]),[],n_vs)',spatial_downsample);
ctx_expl_var_incorrect = AP_spatial_explained_var(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv(incorrect_trials,use_t,:),[2,1,3]),[],n_vs)', ...
    reshape(permute(fluor_taskpred_allcat(incorrect_trials,use_t,:),[2,1,3]),[],n_vs)',spatial_downsample);

ctx_expl_var_visual = AP_spatial_explained_var(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv(visual_trials,use_t,:),[2,1,3]),[],n_vs)', ...
    reshape(permute(fluor_taskpred_allcat(visual_trials,use_t,:),[2,1,3]),[],n_vs)',spatial_downsample);
ctx_expl_var_zero = AP_spatial_explained_var(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv(zero_trials,use_t,:),[2,1,3]),[],n_vs)', ...
    reshape(permute(fluor_taskpred_allcat(zero_trials,use_t,:),[2,1,3]),[],n_vs)',spatial_downsample);

figure;
subplot(2,3,1);
imagesc(ctx_expl_var_correct);
axis image off; 
colormap(brewermap([],'*RdBu'))
caxis([-1,1]); c = colorbar; ylabel(c,'Variance explained');
AP_reference_outline('ccf_aligned','k');
title('Correct trials');

subplot(2,3,2);
imagesc(ctx_expl_var_incorrect);
axis image off; 
colormap(brewermap([],'*RdBu'))
caxis([-1,1]); c = colorbar; ylabel(c,'Variance explained');
AP_reference_outline('ccf_aligned','k');
title('Incorrect trials');

subplot(2,3,3);
imagesc(ctx_expl_var_correct - ctx_expl_var_incorrect);
axis image off; 
colormap(brewermap([],'*RdBu'))
caxis([-1,1]); c = colorbar; ylabel(c,'Variance explained');
AP_reference_outline('ccf_aligned','k');
title('Correct - incorrect');

subplot(2,3,4);
imagesc(ctx_expl_var_visual);
axis image off; 
colormap(brewermap([],'*RdBu'))
caxis([-1,1]); c = colorbar; ylabel(c,'Variance explained');
AP_reference_outline('ccf_aligned','k');
title('Visual trials');

subplot(2,3,5);
imagesc(ctx_expl_var_zero);
axis image off; 
colormap(brewermap([],'*RdBu'))
caxis([-1,1]); c = colorbar; ylabel(c,'Variance explained');
AP_reference_outline('ccf_aligned','k');
title('Zero trials');

subplot(2,3,6);
imagesc(ctx_expl_var_visual - ctx_expl_var_zero);
axis image off; 
colormap(brewermap([],'*RdBu'))
caxis([-1,1]); c = colorbar; ylabel(c,'Variance explained');
AP_reference_outline('ccf_aligned','k');
title('Visual - zero');


%% Compare kernels with domain maps

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

%% Example data with task regression overlay

animal = 'AP025'; 
day = '2017-10-01'; 
experiment = 1; 
str_align = 'none'; 
verbose = true; 
AP_load_experiment;

% Align U deconvolve V, set components to use
use_components = 1:200;
aUdf = AP_align_widefield(animal,day,Udf);
fVdf_deconv = AP_deconv_wf(fVdf);

% Regress task events to fluorescence
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

% Regress task -> fluor
event_aligned_V_deconv = ...
    interp1(frame_t,fVdf_deconv(use_components,:)',t_peri_event);

baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,:),[],size(event_aligned_V_deconv,3)))';
activity = interp1(frame_t,fVdf_deconv(use_components,:)',time_bin_centers)' - baseline;

[~,fluor_taskpred_short,~,~] = ...
    AP_regresskernel(regressors,activity,sample_shifts, ...
    lambda,zs,cvfold,return_constant,use_constant);

% (interpolate the task-predicted fluorescence and add back baseline)
fluor_taskpred = interp1(time_bin_centers,fluor_taskpred_short',frame_t)' + baseline;

%%% Plot example data

% Set time to plot
plot_t = [134,152];

raster_fig = figure;

% (wheel velocity)
wheel_axes = subplot(6,1,1);
plot_wheel_idx = Timeline.rawDAQTimestamps >= plot_t(1) & ...
    Timeline.rawDAQTimestamps <= plot_t(2);
plot(wheel_axes,Timeline.rawDAQTimestamps(plot_wheel_idx), ...
    wheel_velocity(plot_wheel_idx),'k','linewidth',2);
ylabel('Wheel velocity');

% (stimuli)
stim_col = colormap_BlueWhiteRed(5);
[~,trial_contrast_idx] = ...
    ismember(trial_conditions(:,1).*trial_conditions(:,2),unique(contrasts'.*sides),'rows');
stim_lines = arrayfun(@(x) line(wheel_axes,repmat(stimOn_times(x),1,2),ylim(wheel_axes),'color', ...
    stim_col(trial_contrast_idx(x),:),'linewidth',2), ...
    find(stimOn_times >= plot_t(1) & stimOn_times <= plot_t(2)));

% (movement starts)
move_col = [0.6,0,0.6;0,0.6,0];
[~,trial_choice_idx] = ismember(trial_conditions(:,3),[-1;1],'rows');
move_lines = arrayfun(@(x) line(wheel_axes,repmat(wheel_move_time(x),1,2),ylim(wheel_axes),'color', ...
    move_col(trial_choice_idx(x),:),'linewidth',2), ...
    find(wheel_move_time >= plot_t(1) & wheel_move_time <= plot_t(2)));

% (go cues)
go_col = [0.8,0.8,0.2];
go_cue_times = signals_events.interactiveOnTimes(1:n_trials);
go_cue_lines = arrayfun(@(x) line(wheel_axes,repmat(go_cue_times(x),1,2),ylim(wheel_axes),'color', ...
    go_col,'linewidth',2,'linestyle','--'), ...
    find(go_cue_times >= plot_t(1) & go_cue_times <= plot_t(2)));

% (outcomes)
outcome_col = [0,0,0.8;0.5,0.5,0.5];
reward_lines = arrayfun(@(x) line(wheel_axes,repmat(reward_t_timeline(x),1,2),ylim(wheel_axes),'color', ...
    outcome_col(1,:),'linewidth',2,'linestyle','--'), ...
    find(reward_t_timeline >= plot_t(1) & reward_t_timeline <= plot_t(2)));
punish_times = signals_events.responseTimes(trial_outcome == -1);
punish_lines = arrayfun(@(x) line(wheel_axes,repmat(punish_times(x),1,2),ylim(wheel_axes),'color', ...
    outcome_col(2,:),'linewidth',2,'linestyle','--'), ...
    find(punish_times >= plot_t(1) & punish_times <= plot_t(2)));

% (striatum raster)
raster_axes = subplot(6,1,2:4,'YDir','reverse'); hold on;
plot_spikes = spike_times_timeline >= plot_t(1) & ...
    spike_times_timeline <= plot_t(2) & ...
    spike_depths >= str_depth(1) & spike_depths <= str_depth(2);
plot(raster_axes,spike_times_timeline(plot_spikes),spike_depths(plot_spikes),'.k');
ylabel('Depth (\mum)');
xlabel('Time (s)')

% (fluorescence from select ROIs)
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
roi_trace = AP_svd_roi(aUdf(:,:,use_components),fVdf_deconv(use_components,:),[],[],cat(3,wf_roi.mask));
roi_trace_taskpred = AP_svd_roi(aUdf(:,:,use_components),fluor_taskpred,[],[],cat(3,wf_roi.mask));

plot_rois = [1,7,9];
fluor_spacing = 70;
fluor_axes = subplot(6,1,5:6); hold on;
plot_fluor_idx = frame_t >= plot_t(1) & frame_t <= plot_t(2);
AP_stackplot(roi_trace(plot_rois,plot_fluor_idx)', ...
    frame_t(plot_fluor_idx),fluor_spacing,false,[0,0.7,0],{wf_roi(plot_rois).area});
AP_stackplot(roi_trace_taskpred(plot_rois,plot_fluor_idx)', ...
    frame_t(plot_fluor_idx),fluor_spacing,false,'b',{wf_roi(plot_rois).area});


linkaxes([wheel_axes,raster_axes,fluor_axes],'x');


% % Write legend
% [~,unique_contrasts_h] = unique(trial_contrast_idx);
% [~,unique_move_h] = unique(trial_choice_idx(trial_choice_idx > 0));
% legend([stim_lines(unique_contrasts_h),move_lines(unique_move_h), ...
%     go_cue_lines(1),reward_lines(1),punish_lines(1)], ...
%     [cellfun(@(x) ['Stim ' num2str(x)],num2cell(unique(contrasts'.*sides)),'uni',false); ...
%     {'Move L';'Move R';'Go cue';'Reward';'Punish'}]);


% Plot fluorescence at regular intervals within time range

% plot_frames_idx = [4721,4833,4924,5008]; % (across 4 trials)
plot_frames_idx = [4827,4831,4834,4838]; % (within one trial)

plot_frames = svdFrameReconstruct(aUdf(:,:,use_components),fVdf_deconv(use_components,plot_frames_idx));
plot_frames_taskpred = svdFrameReconstruct(aUdf(:,:,use_components),fluor_taskpred(:,plot_frames_idx));

wf_fig = figure;
for curr_frame = 1:length(plot_frames_idx)
    subplot(2,length(plot_frames_idx),curr_frame);
    imagesc(plot_frames(:,:,curr_frame));
    colormap(crameri('cork'));
    caxis([-0.03,0.03]);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title(sprintf('Deconv - %0.2fs',frame_t(plot_frames_idx(curr_frame))));
    axis image off;
    
    subplot(2,length(plot_frames_idx),length(plot_frames_idx) + curr_frame);
    imagesc(plot_frames_taskpred(:,:,curr_frame));
    colormap(crameri('cork'));
    caxis([-0.03,0.03]);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title(sprintf('Task predicted - %0.2fs',frame_t(plot_frames_idx(curr_frame))));
    axis image off;
end

% (draw lines on the ROI plot where the frames were taken from)
frame_lines = arrayfun(@(x) line(fluor_axes,repmat(frame_t(x),1,2), ...
    ylim(fluor_axes),'color','k','linewidth',2),plot_frames_idx);


% Plot ROIs
figure; hold on
set(gca,'YDir','reverse');
AP_reference_outline('ccf_aligned','k');
for curr_roi = plot_rois
    curr_roi_boundary = cell2mat(bwboundaries(wf_roi(curr_roi).mask()));
    patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),[0,0.8,0]);   
    text(nanmean(curr_roi_boundary(:,2)),nanmean(curr_roi_boundary(:,1)), ...
        wf_roi(curr_roi).area,'FontSize',12,'HorizontalAlignment','center')
end
axis image off;


%% Striatal unit regression examples

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






