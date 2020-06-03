%% ~~~~~~~~~~~~ Load

%% Load SUA data

disp('Loading SUA data');

data_fn = 'G:\JF_single_cell_data\trial_activity_choiceworld.mat';
load(data_fn);

% Find fields with experiment data (cell arrays with length animals)
data_struct_fieldnames = fieldnames(trial_data_all);
experiment_fields = cellfun(@(curr_field) ...
    length(trial_data_all.(curr_field)) == length(trial_data_all.animals) && ...
    iscell(trial_data_all.(curr_field)) && ...
    any(cellfun(@(x) iscell(x),trial_data_all.(curr_field))),data_struct_fieldnames);

% Load pre-marked experiments to exclude and cut out bad ones
if exist('exclude_data','var') && exclude_data
    exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
    exclude_fn = 'bhv_use_experiments';
    load([exclude_path filesep exclude_fn]);
    
    % Pull out used experiments for the animals loaded
    use_experiments_animals = ismember(trial_data_all.animals,{bhv_use_experiments.animals});
    use_experiments = {bhv_use_experiments(use_experiments_animals).use_experiments}';
    
    % (old, when multiple kinds of exclusions
%     exclude_fn{1} = 'bhv_use_experiments';
%     % exclude_fn{2} = 'expl_var_use_experiments';
%     use_experiments_all = {};
%     for curr_exclude = 1:length(exclude_fn)
%         curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
%         use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
%     end
%     use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
%         1:size(use_experiments_all,2),'uni',false)';
    
    % Cut out bad experiments for any experiment data fields
    if ~isempty(use_experiments)
        for curr_field = data_struct_fieldnames(experiment_fields)'
            trial_data_all.(cell2mat(curr_field)) = cellfun(@(x,expts) ...
                x(expts),trial_data_all.(cell2mat(curr_field)), ...
                use_experiments,'uni',false);
        end
    end
end

% If any animals don't have any data - throw away
nodata_animals = cellfun(@(x) isempty(x),trial_data_all.trial_info_all);
trial_data_all.animals(nodata_animals) = [];
for curr_field = data_struct_fieldnames(experiment_fields)'
            trial_data_all.(cell2mat(curr_field))(nodata_animals) = [];
end

% Unpack data structure into workspace then throw away
arrayfun(@(x) assignin('base',cell2mat(x),trial_data_all.(cell2mat(x))),data_struct_fieldnames);
clear trial_data_all

% Get sample rate and set "baseline" time
sample_rate = 1/mean(diff(t));
t_baseline = t < 0;

% Get if this is a task dataset
task_dataset = exist('outcome_all','var');

% Concatenate trial info data
trial_info_fields = fieldnames(trial_info_all{end}{end});
trial_info_allcat = cell2struct(arrayfun(@(curr_field) ...
    cell2mat(cellfun(@(x) x(curr_field), ...
    cellfun(@struct2cell,vertcat(trial_info_all{:}),'uni',false))), ...
    1:length(trial_info_fields),'uni',false),trial_info_fields,2);

% Concatenate wheel
wheel_allcat = cell2mat(vertcat(wheel_all{:}));

% % (movement from mousecam if exists: normalize and concatenate)
% if exist('movement_all','var')
%     movement_all_norm = cellfun(@(x) cellfun(@(x) x./nanstd(x(:)), ...
%         x,'uni',false),movement_all,'uni',false);
%     movement_allcat = cell2mat(vertcat(movement_all_norm{:}));
% end

%%% Get task/stim-relevant

if task_dataset
    
    % Get trial information
    trial_stim_allcat = trial_info_allcat.stimulus; % old: diff(trial_info_allcat.stimulus,[],2);
    trial_choice_allcat = -(trial_info_allcat.response-1.5)*2;    
    trial_outcome_allcat = trial_info_allcat.outcome;
    
    % Get reaction time and t index for movement onset
    move_t = trial_info_allcat.stim_to_move;
    [~,move_idx] = min(abs(move_t - t),[],2);
    
    % Get outcome time
    outcome_allcat = cell2mat(vertcat(outcome_all{:}));
    [~,outcome_idx] = max(any(outcome_allcat,3),[],2);
    outcome_t = t(outcome_idx)';
    
    % Get wheel velocity
    wheel_velocity_allcat = wheel_allcat;
    [max_speed,max_vel_idx] = max(abs(wheel_velocity_allcat(:,:,1).* ...
        (bsxfun(@times,wheel_velocity_allcat,trial_choice_allcat) > 0)),[],2);
    max_vel = max_speed.*trial_choice_allcat;
    wheel_velocity_allcat_move = wheel_velocity_allcat;
    
    t_leeway = -t(1);
    leeway_samples = round(t_leeway*(sample_rate));
    for i = 1:size(wheel_velocity_allcat,1)
        wheel_velocity_allcat_move(i,:,:) = circshift(wheel_velocity_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
    end 
   
elseif isfield(trial_info_allcat,'stimulus')  
    
    if length(unique(trial_info_allcat.stimulus)) == 3
        % Use stim IDs
        trial_stim_allcat = trial_info_allcat.stimulus;
        % For both mpep and signals, ID is 1 = left, 2 = center, 3 = right
        trial_stim_allcat(trial_stim_allcat == 1) = -1;
        trial_stim_allcat(trial_stim_allcat == 2) = 0;
        trial_stim_allcat(trial_stim_allcat == 3) = 1;
    else
        % Passive choiceworld uses stimID = side*contrast
        trial_stim_allcat = trial_info_allcat.stimulus;
    end

end


% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_recording;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));

% Concatenate depths
depth_allcat = cell2mat(cellfun(@(x) x(~isnan(x)), ...
    vertcat(allDepths{:}),'uni',false));

% Concatenate groups and make logical vectors
groups_allcat = cell2mat(cellfun(@(x) x(~isnan(x)), ...
    vertcat(allGroups{:}),'uni',false));

% (groups are cell type * depth: msn, fsi, tan, th, ??)
n_domains = 3; % hardcoded: I think not stored
n_celltypes = max(groups_allcat)./n_domains;
celltype_allcat = ceil(groups_allcat./n_domains);
celltype_labels = {'MSN','FSI','TAN','TH','Other'};
if length(celltype_labels) ~= n_celltypes
    error('Mismatching number of celltypes and labels')
end
domain_allcat = mod(groups_allcat,n_domains) + ...
    n_domains*(mod(groups_allcat,n_domains) == 0);

% Get maps for all cells
% (load master U)
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');
% (use t = 0 kernels for all cells, hardcoded at the moment)
use_k_frame = 5;
ctx_str_k_frame = cell2mat(cellfun(@(x) reshape(x(:,use_k_frame,:),100,[]), ...
    vertcat(ctx_str_k_all{:})','uni',false));
ctx_str_k_px = svdFrameReconstruct(U_master(:,:,1:100),ctx_str_k_frame);


% Deconvolve fluorescence
fluor_deconv_allcat_exp = cellfun(@(x) AP_deconv_wf(x),vertcat(fluor_all{:}),'uni',false);


disp('Finished.')


%% Align by depth and domain

%%% Align by relative depth 
ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
ephys_align_fn = ['ephys_depth_align.mat'];
load([ephys_align_path filesep ephys_align_fn]);

% Pull out relative depths for each cell
depth_aligned = cell(size(mua_all));
for curr_animal = 1:length(mua_all)
    curr_animal_idx = strcmp(animals{curr_animal},{ephys_depth_align.animal});
    for curr_day = 1:length(mua_all{curr_animal_idx})
        % (relative to end)
        depth_aligned{curr_animal}{curr_day} = ...
            allDepths{curr_animal}{curr_day}(~isnan(allDepths{curr_animal}{curr_day})) - ...
            ephys_depth_align(curr_animal_idx).str_depth(curr_day,2);
%         % (relative to start:end)
%         depth_aligned{curr_animal}{curr_day} = ...
%             (allDepths{curr_animal}{curr_day}(~isnan(allDepths{curr_animal}{curr_day})) - ...
%             ephys_depth_align(curr_animal_idx).str_depth(curr_day,1))./ ...
%             diff(ephys_depth_align(curr_animal_idx).str_depth(curr_day,:));
    end
end

depth_aligned_allcat = cell2mat(horzcat(depth_aligned{:})');

%%% Align by domain

% Load domain alignment
n_domains = 3;
ephys_kernel_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
ephys_kernel_align_fn = ['ephys_kernel_align_' num2str(n_domains) '_depths.mat'];
load([ephys_kernel_align_path filesep ephys_kernel_align_fn]);

% Loop through cells, pull out domains
domain_aligned = cell(size(mua_all));
for curr_animal = 1:length(mua_all)
    curr_animal_idx = strcmp(animals{curr_animal},{ephys_kernel_align.animal});
    for curr_day = 1:length(mua_all{curr_animal_idx})
        
        str_depth = ephys_depth_align(curr_animal_idx).str_depth(curr_day,:);
        kernel_match = ephys_kernel_align(curr_animal_idx).kernel_match{curr_day};

        % Get depth groups
        n_depths = round(diff(str_depth)/200);
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        
        % Get domain depth boundaries
        kernel_match_boundary_idx = unique([1;find(diff(kernel_match) ~= 0)+1;n_depths+1]);      
        kernel_match_depth_edges = depth_group_edges(kernel_match_boundary_idx);
        % (extend first and last forever - sometimes cells a little past
        % the border but classified as str probably because changed border
        % calculation slightly)
        kernel_match_depth_edges(1) = -Inf;
        kernel_match_depth_edges(end) = Inf;
        
        kernel_match_idx = kernel_match(kernel_match_boundary_idx(1:end-1));
        
        % Assign neurons to domains
        domain_aligned{curr_animal}{curr_day} = ...
            discretize(allDepths{curr_animal}{curr_day}(~isnan(allDepths{curr_animal}{curr_day})), ...
            kernel_match_depth_edges,kernel_match_idx);
 
    end
end

domain_aligned_allcat = cell2mat(horzcat(domain_aligned{:})');

% Plot for sanity check
figure; plot3(depth_allcat,depth_aligned_allcat,domain_aligned_allcat,'.k')
xlabel('Raw depth');
ylabel('Aligned depth');
zlabel('Aligned domain');
axis vis3d


%% Plot type/domain counts

domain_type_count= accumarray([domain_aligned_allcat,celltype_allcat],1);

figure;
bar(domain_type_count,'stacked');
xlabel('Domain');
ylabel('Count');
legend(celltype_labels);

%% ~~~~~~~~~~~~ Activity

%% Rasters (all cells)

trial_stim_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
move_t_exp = mat2cell(move_t,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

% Sort cells by max activity time
mua_allcat_stimalign_exp = vertcat(mua_all{:});

act_mean = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1 & ...
    trial_rand(1:size(act,1)) < 0.5,:,:),1)),mua_allcat_stimalign_exp, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false)')';

[~,max_idx] = max(AP_deconv_wf(act_mean,true),[],2);
[~,sort_idx] = sort(max_idx);

figure;imagesc(t,[],zscore(AP_deconv_wf(act_mean(sort_idx,:),true),[],2));
line([0,0],ylim,'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
caxis([-3,3]);
colormap(brewermap([],'*RdBu'));

% Pad and plot all rasters
plot_trials_1 = cellfun(@(stim,rxn,outcome) ...
    stim > 0, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false);
plot_trials_2 = cellfun(@(stim,rxn,outcome) ...
    stim < 0, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false);

max_trials = max(cellfun(@(x,y) sum(x)+sum(y),plot_trials_1,plot_trials_2));
mua_allcat_stimalign_exp_pad = cell2mat(permute(cellfun(@(act,trials_1,trials_2) ...
    padarray(act([find(trials_1);find(trials_2)],:,:), ...
    [max_trials-(sum(trials_1)+sum(trials_2)),0,0],NaN,'post'), ...
    mua_allcat_stimalign_exp,plot_trials_1,plot_trials_2,'uni',false),[2,3,1]));

cell_string = arrayfun(@(x) ...
    [' Sort idx: ' num2str(x) ...
    ', Cell: ' num2str(sort_idx(x)) ...
    ', Domain: ' num2str(domain_aligned_allcat(x)), ...
    ', ' celltype_labels{celltype_allcat(x)}],1:length(sort_idx),'uni',false);
AP_image_scroll(mua_allcat_stimalign_exp_pad(:,:,sort_idx),cell_string);
line(repmat(find(t >= 0,1),2,1),ylim,'linestyle','--','color','r');
colormap(brewermap([],'Greys'));
caxis([0,50]);



%% Rasters (select cells) 

% Set cells and alignment to plot
plot_cells = move_cells & ismember(celltype_allcat,[1,2]) & domain_aligned_allcat == 1;
plot_align = 2; % (stim,move,outcome)

% Get events by experiment
trial_stim_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
move_t_exp = mat2cell(move_t,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

% Get aligned activity
% (stim)
mua_allcat_stimalign_exp = vertcat(mua_all{:});

% (move aligned)
move_idx_exp = mat2cell(move_idx,use_split,1);
mua_allcat_movealign_exp = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_movealign_exp)
   for curr_trial = 1:size(mua_allcat_movealign_exp{curr_exp},1)
       mua_allcat_movealign_exp{curr_exp}(curr_trial,:,:) = ...
           circshift(mua_allcat_movealign_exp{curr_exp}(curr_trial,:,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
   end
end

% (outcome aligned)
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
mua_allcat_outcomealign_exp = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_outcomealign_exp)
    for curr_trial = 1:size(mua_allcat_outcomealign_exp{curr_exp},1)
        mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:,:) = ...
            circshift(mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:,:), ...
            -outcome_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
    end
end


% Pad and plot rasters
plot_trials_1 = cellfun(@(stim,rxn,outcome) ...
    stim > 0, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false);
plot_trials_2 = cellfun(@(stim,rxn,outcome) ...
    stim < 0, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false);

switch plot_align
    case 1
        plot_act = mua_allcat_stimalign_exp;
    case 2
        plot_act = mua_allcat_movealign_exp;
    case 3
        plot_act = mua_allcat_outcomealign_exp;
end

max_trials = max(cellfun(@(x,y) sum(x)+sum(y),plot_trials_1,plot_trials_2));
act_exp_pad = cell2mat(permute(cellfun(@(act,trials_1,trials_2) ...
    padarray(act([find(trials_1);find(trials_2)],:,:), ...
    [max_trials-(sum(trials_1)+sum(trials_2)),0,0],NaN,'post'), ...
    plot_act,plot_trials_1,plot_trials_2,'uni',false),[2,3,1]));

% act_exp_pad = cell2mat(permute(cellfun(@(act_stim,act_move,act_outcome,trials_1,trials_2) ...
%     padarray( ...
%     [act_stim([find(trials_1);find(trials_2)],:,:), ...
%     act_move([find(trials_1);find(trials_2)],:,:), ...
%     act_outcome([find(trials_1);find(trials_2)],:,:)], ...
%     [max_trials-(sum(trials_1)+sum(trials_2)),0,0],NaN,'post'), ...
%     mua_allcat_stimalign_exp,mua_allcat_movealign_exp,mua_allcat_outcomealign_exp, ...
%     plot_trials_1,plot_trials_2,'uni',false),[2,3,1]));

plot_cells_idx = find(plot_cells);
cell_string = arrayfun(@(x) ...
    [' Sort idx: ' num2str(x) ...
    ', Cell: ' num2str(plot_cells_idx(x)) ...
    ', Domain: ' num2str(domain_aligned_allcat(plot_cells_idx(x))), ...
    ', ' celltype_labels{celltype_allcat(plot_cells_idx(x))}],1:length(plot_cells_idx),'uni',false);
AP_image_scroll(act_exp_pad(:,:,plot_cells),cell_string);
line(repmat(find(t >= 0,1),2,1),ylim,'linestyle','--','color','r');
line(repmat(find(t >= 0,1),2,1)+length(t),ylim,'linestyle','--','color','r');
line(repmat(find(t >= 0,1),2,1)+length(t)*2,ylim,'linestyle','--','color','r');

colormap(brewermap([],'Greys'));
caxis([0,50]);


%% Stim/move/reward sort and plot

plot_cells = true(size(celltype_allcat));
% plot_cells = celltype_allcat == 1;

trial_stim_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
move_t_exp = mat2cell(move_t,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

% Align activity

% (stim aligned)
mua_allcat_stimalign_exp = vertcat(mua_all{:});

% (move aligned)
move_idx_exp = mat2cell(move_idx,use_split,1);
mua_allcat_movealign_exp = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_movealign_exp)
   for curr_trial = 1:size(mua_allcat_movealign_exp{curr_exp},1)
       mua_allcat_movealign_exp{curr_exp}(curr_trial,:,:) = ...
           circshift(mua_allcat_movealign_exp{curr_exp}(curr_trial,:,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
   end
end

% (outcome aligned)
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
mua_allcat_outcomealign_exp = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_outcomealign_exp)
    for curr_trial = 1:size(mua_allcat_outcomealign_exp{curr_exp},1)
        mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:,:) = ...
            circshift(mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:,:), ...
            -outcome_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
    end
end

% Split trials in 2 using random vector
trial_rand = rand(1000,1);

% Get sorting order from stim-aligned
sort_act = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1 & ...
    trial_rand(1:size(act,1)) < 0.5,:,:),1)),mua_allcat_stimalign_exp, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false)')';
sort_act_smooth = AP_deconv_wf(sort_act(plot_cells,:),true);

[~,max_idx] = max(sort_act_smooth,[],2);
[~,sort_idx] = sort(max_idx);

figure;imagesc(t,[],zscore(sort_act_smooth(sort_idx,:),[],2));
colormap(brewermap([],'*RdBu'));
caxis([-3,3]);
ylabel('Cells (sorted)');
xlabel('Time from stim');

for curr_event = 1:3
    switch curr_event
        case 1
            curr_trials_1 = cellfun(@(stim,rxn,outcome) ...
                stim > 0 & rxn < 0.5 & outcome == 1, ...
                trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false);
            curr_trials_2 = cellfun(@(stim,rxn,outcome) ...
                stim < 0 & rxn < 0.5 & outcome == 1, ...
                trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false);
            act_labels = {'Ipsi stim','Contra stim','Difference'};
        case 2
            curr_trials_1 = cellfun(@(stim,rxn,outcome) ...
                rxn < 0.5, ...
                trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false);
            curr_trials_2 = cellfun(@(stim,rxn,outcome) ...
                rxn > 0.5, ...
                trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false);
            act_labels = {'Move early','Move late','Difference'};
        case 3
            curr_trials_1 = cellfun(@(stim,rxn,outcome) ...
                rxn < 0.5 & outcome == 1, ...
                trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false);
            curr_trials_2 = cellfun(@(stim,rxn,outcome) ...
                rxn < 0.5 & outcome == -1, ...
                trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false);
            act_labels = {'Reward','No reward','Difference'};
    end
    
    
  
    act_1 = cell2mat(cellfun(@(act,curr_trials) ...
        squeeze(nanmean(act(curr_trials & ...
        trial_rand(1:size(act,1)) >= 0.5,:,:),1)),mua_allcat_stimalign_exp, ...
        curr_trials_1,'uni',false)')';
    
    act_2 = cell2mat(cellfun(@(act,curr_trials) ...
        squeeze(nanmean(act(curr_trials & ...
        trial_rand(1:size(act,1)) >= 0.5,:,:),1)),mua_allcat_stimalign_exp, ...
        curr_trials_2,'uni',false)')';
    
    act_1_smooth = AP_deconv_wf(act_1(plot_cells,:),true);
    act_2_smooth = AP_deconv_wf(act_2(plot_cells,:),true);
    
    act_1_smooth(isnan(act_1_smooth)) = 0;
    act_2_smooth(isnan(act_2_smooth)) = 0;
    
    figure;
    
    subplot(1,3,1);
    imagesc(t,[],zscore(act_1_smooth(sort_idx,:),[],2));
    caxis([-3,3])
    line([0,0],ylim,'color','k','linestyle','--');
    line([0.5,0.5],ylim,'color','k','linestyle','--');
    xlabel('Time');
    ylabel('Cell (sorted)');
    title(act_labels{1});
    
    subplot(1,3,2);
    imagesc(t,[],zscore(act_2_smooth(sort_idx,:),[],2));
    caxis([-3,3])
    line([0,0],ylim,'color','k','linestyle','--');
    line([0.5,0.5],ylim,'color','k','linestyle','--');
    xlabel('Time');
    ylabel('Cell (sorted)');
    title(act_labels{2});
    
    subplot(1,3,3);
    imagesc(t,[],imgaussfilt(zscore(act_1_smooth(sort_idx,:)-act_2_smooth(sort_idx,:),[],2),3));
    caxis([-3,3])
    line([0,0],ylim,'color','k','linestyle','--');
    line([0.5,0.5],ylim,'color','k','linestyle','--');
    xlabel('Time');
    ylabel('Cell (sorted)');
    title(act_labels{3});
    
    colormap(brewermap([],'*RdBu'));
    
end

linkaxes(get(gcf,'Children'),'xy');



%% Pull out stim/move/reward cells

trial_stim_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
move_t_exp = mat2cell(move_t,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

% (stim aligned)
mua_allcat_stimalign_exp = vertcat(mua_all{:});

% (move aligned)
move_idx_exp = mat2cell(move_idx,use_split,1);
mua_allcat_movealign_exp = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_movealign_exp)
   for curr_trial = 1:size(mua_allcat_movealign_exp{curr_exp},1)
       mua_allcat_movealign_exp{curr_exp}(curr_trial,:,:) = ...
           circshift(mua_allcat_movealign_exp{curr_exp}(curr_trial,:,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
   end
end

% (outcome aligned)
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
mua_allcat_outcomealign_exp = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_outcomealign_exp)
    for curr_trial = 1:size(mua_allcat_outcomealign_exp{curr_exp},1)
        mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:,:) = ...
            circshift(mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:,:), ...
            -outcome_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
    end
end

% Get activity by conditions
stim_contra_act = AP_deconv_wf(cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim > 0 & rxn < 0.5,:,:),1)), ...
    mua_allcat_stimalign_exp,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false)')',true);

stim_ipsi_act = AP_deconv_wf(cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim <= 0 & rxn < 0.5,:,:),1)), ...
    mua_allcat_stimalign_exp,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false)')',true);

move_early_act = AP_deconv_wf(cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim <= 0 & rxn < 0.5,:,:),1)), ...
    mua_allcat_stimalign_exp,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false)')',true);

move_late_act = AP_deconv_wf(cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim <= 0 & rxn >= 0.5,:,:),1)), ...
    mua_allcat_stimalign_exp,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false)')',true);

reward_act = AP_deconv_wf(cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(outcome == 1,:,:),1)), ...
    mua_allcat_outcomealign_exp,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false)')',true);

no_reward_act = AP_deconv_wf(cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(outcome == -1,:,:),1)), ...
    mua_allcat_outcomealign_exp,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false)')',true);


% Get max aligned activity
% (unused at the moment)
use_trials = cellfun(@(stim,rxn,outcome) ...
    stim > 0 & rxn < 0.5 & outcome == 1, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false);

avg_act_aligned = AP_deconv_wf( ...
    cell2mat(cellfun(@(use_trials,stim_aligned,move_aligned,outcome_aligned) ...
    permute(cat(1,nanmean(stim_aligned(use_trials,:,:),1), ...
    nanmean(move_aligned(use_trials,:,:),1), ...
    nanmean(outcome_aligned(use_trials,:,:),1)),[3,2,1]), ...
    use_trials,mua_allcat_stimalign_exp, ...
    mua_allcat_movealign_exp,mua_allcat_outcomealign_exp,'uni',false)),true);

[~,max_aligned_idx] = max(max(avg_act_aligned,[],2),[],3);

% % Pull out cells
% % (first criteria)
% use_t = t > 0 & t < 0.5;
% [~,stim_max_idx] = max(stim_contra_act,[],2);
% stim_cells = ismember(stim_max_idx,find(use_t)) & ...
%     (max(stim_contra_act(:,use_t),[],2) > ...
%     2*max(stim_ipsi_act(:,use_t),[],2)) & max_aligned_idx == 1;
% 
% use_t = t > 0 & t < 0.5;
% [~,move_max_idx] = max(move_early_act,[],2);
% move_cells = ismember(move_max_idx,find(use_t)) & ...
%     (max(move_early_act(:,use_t),[],2) > ...
%     2*max(move_late_act(:,use_t),[],2)) & max_aligned_idx == 2;
% 
% use_t = t > 0 & t < 0.5;
% [~,reward_max_idx] = max(reward_act,[],2);
% reward_cells = ismember(reward_max_idx,find(use_t)) & ...
%     (max(reward_act(:,use_t),[],2) > ...
%     2*max(no_reward_act(:,use_t),[],2)) & max_aligned_idx == 3;

% Pull out cells
% NOTE: these criteria have changed a lot
use_t = t > 0 & t < 0.5;
[~,stim_max_idx] = max(stim_contra_act,[],2);
stim_cells = ismember(stim_max_idx,find(use_t)) & ...
    (max(stim_contra_act(:,use_t),[],2) > ...
    max(stim_ipsi_act(:,use_t),[],2)) & trial_corr_max_align == 1;% & ...
%     taskpred_partial_r2_allcat(:,1) > 0;

use_t = t > 0 & t < 0.5;
[~,move_max_idx] = max(move_early_act,[],2);
move_cells = ismember(move_max_idx,find(use_t)) & ...
    (max(move_early_act(:,use_t),[],2) > ...
    max(move_late_act(:,use_t),[],2)) & trial_corr_max_align == 2;% & ...
%     taskpred_partial_r2_allcat(:,2) > 0;

use_t = t > 0 & t < 0.5;
[~,reward_max_idx] = max(reward_act,[],2);
reward_cells = ismember(reward_max_idx,find(use_t)) & ...
    (max(reward_act(:,use_t),[],2) > ...
    max(no_reward_act(:,use_t),[],2)) & trial_corr_max_align == 3;% & ...
%     taskpred_partial_r2_allcat(:,4) > 0;


% Plot cells by condition
figure; 
subplot(1,3,1);
imagesc(t,[],zscore([stim_contra_act(stim_cells,:);stim_ipsi_act(stim_cells,:)],[],2));
line([0,0],ylim,'color','k','linestyle','--');
line(xlim,repmat(sum(stim_cells),2,1),'color','k');
xlabel('Time from stim');
title('Stim cells');
caxis([-4,4])
colormap(brewermap([],'*RdBu'));

subplot(1,3,2);
imagesc(t,[],zscore([move_early_act(move_cells,:);move_late_act(move_cells,:)],[],2));
line([0,0],ylim,'color','k','linestyle','--');
line(xlim,repmat(sum(move_cells),2,1),'color','k');
xlabel('Time from stim');
title('Move cells');
caxis([-4,4])
colormap(brewermap([],'*RdBu'));

subplot(1,3,3);
imagesc(t,[],zscore([reward_act(reward_cells,:);no_reward_act(reward_cells,:)],[],2));
line([0,0],ylim,'color','k','linestyle','--');
line(xlim,repmat(sum(reward_cells),2,1),'color','k');
xlabel('Time from outcome');
title('Reward cells');
caxis([-4,4])
colormap(brewermap([],'*RdBu'));


%% Plot cell response count/frac by type/domain (after above) 

domain_type_count = ...
    accumarray([domain_aligned_allcat,celltype_allcat, ...
    stim_cells+1,move_cells+1,reward_cells+1],1);

cols = zeros(2,2,2,3);
cols(2,:,:,1) = 1;
cols(:,2,:,2) = 1;
cols(:,:,2,3) = 1;

figure;
for curr_domain = 1:n_domains
    for curr_type = 1:n_celltypes
        subplot(n_domains,n_celltypes, ...
            (curr_domain-1)*n_celltypes+curr_type);
        
        curr_counts = squeeze(domain_type_count(curr_domain,curr_type,:,:,:));
        % hack the colors
        curr_counts(curr_counts == 0) = 0.01;
        
        pie(curr_counts)

        title([num2str(curr_domain) ' ' celltype_labels{curr_type}]);
    end
end
legend({'Nothing','Stim','Move','Stim+Move','Reward','Move+Reward','Stim+Move+Reward'})
colormap(lines)

figure; 
subplot(1,2,1); hold on
plot(rand(sum(stim_cells),1),depth_aligned_allcat(stim_cells),'.r','MarkerSize',10);
plot(rand(sum(move_cells),1),depth_aligned_allcat(move_cells),'.m','MarkerSize',10);
plot(rand(sum(reward_cells),1),depth_aligned_allcat(reward_cells),'.b','MarkerSize',10);
ylabel('Aligned depth');
set(gca,'YDir','reverse');

subplot(1,2,2); hold on
depth_bin_edges = linspace(min(depth_aligned_allcat),max(depth_aligned_allcat),6);
depth_bin_centers = conv(depth_bin_edges,[1,1]./2,'valid');
plot(histcounts(depth_aligned_allcat(stim_cells),depth_bin_edges)./sum(stim_cells), ...
    depth_bin_centers,'r','linewidth',2);
plot(histcounts(depth_aligned_allcat(move_cells),depth_bin_edges)./sum(move_cells), ...
    depth_bin_centers,'m','linewidth',2);
plot(histcounts(depth_aligned_allcat(reward_cells),depth_bin_edges)./sum(reward_cells), ...
    depth_bin_centers,'b','linewidth',2);
ylabel('Aligned depth');
xlabel('Fraction')
set(gca,'YDir','reverse');
legend({'Stim','Move','Reward'});

linkaxes(get(gcf,'Children'),'y');



%% Split average trial correlation for all cells

trial_stim_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
move_t_exp = mat2cell(move_t,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

% Align activity

% (stim aligned)
mua_allcat_stimalign_exp = vertcat(mua_all{:});

% (move aligned)
move_idx_exp = mat2cell(move_idx,use_split,1);
mua_allcat_movealign_exp = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_movealign_exp)
   for curr_trial = 1:size(mua_allcat_movealign_exp{curr_exp},1)
       mua_allcat_movealign_exp{curr_exp}(curr_trial,:,:) = ...
           circshift(mua_allcat_movealign_exp{curr_exp}(curr_trial,:,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
   end
end

% (outcome aligned)
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
mua_allcat_outcomealign_exp = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_outcomealign_exp)
    for curr_trial = 1:size(mua_allcat_outcomealign_exp{curr_exp},1)
        mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:,:) = ...
            circshift(mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:,:), ...
            -outcome_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
    end
end

% Split trials in 2 and get correlation for each alignment
n_split = 5;
n_cells = sum(cellfun(@(x) size(x,3),vertcat(mua_all{:})));
trial_corr = nan(n_cells,3);
for curr_align = 1:3
    switch curr_align
        case 1
            curr_act = mua_allcat_stimalign_exp;
        case 2
            curr_act = mua_allcat_movealign_exp;
        case 3
            curr_act = mua_allcat_outcomealign_exp;
    end
    
    trial_rand = rand(1000,1);
    rand_split = linspace(0,1,n_split+1);
    
    act_split = cell2mat(arrayfun(@(x) cell2mat(cellfun(@(act,stim,rxn,outcome) ...
        squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1 & ...
        trial_rand(1:size(act,1)) >= rand_split(x) & ...
        trial_rand(1:size(act,1)) <= rand_split(x+1),:,:),1)),curr_act, ...
        trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false)')', ...
        permute(1:n_split,[1,3,2]),'uni',false));
    
    act_split_smooth = AP_deconv_wf(act_split,true);
    
    act_split_corr = arrayfun(@(x) ...
        nanmedian(AP_itril(corr(permute(act_split_smooth(x,:,:),[2,3,1])),-1)), ...
        1:size(act_split,1));
        
    trial_corr(:,curr_align) = act_split_corr;
    
end

[trial_corr_max,trial_corr_max_align] = max(trial_corr,[],2);


%% Trial-trial correlation for all cells

trial_stim_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
move_t_exp = mat2cell(move_t,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

% Align activity

% (stim aligned)
mua_allcat_stimalign_exp = vertcat(mua_all{:});

% (move aligned)
move_idx_exp = mat2cell(move_idx,use_split,1);
mua_allcat_movealign_exp = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_movealign_exp)
   for curr_trial = 1:size(mua_allcat_movealign_exp{curr_exp},1)
       mua_allcat_movealign_exp{curr_exp}(curr_trial,:,:) = ...
           circshift(mua_allcat_movealign_exp{curr_exp}(curr_trial,:,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
   end
end

% (outcome aligned)
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
mua_allcat_outcomealign_exp = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_outcomealign_exp)
    for curr_trial = 1:size(mua_allcat_outcomealign_exp{curr_exp},1)
        mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:,:) = ...
            circshift(mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:,:), ...
            -outcome_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
    end
end

% Smooth data for trial-trial correlations
mua_allcat_stimalign_exp_smooth = cellfun(@(x) AP_deconv_wf(x,true),mua_allcat_stimalign_exp,'uni',false);
mua_allcat_movealign_exp_smooth = cellfun(@(x) AP_deconv_wf(x,true),mua_allcat_movealign_exp,'uni',false);
mua_allcat_outcomealign_exp_smooth = cellfun(@(x) AP_deconv_wf(x,true),mua_allcat_outcomealign_exp,'uni',false);

trial_corr = cell2mat(cellfun(@(stim_act,move_act,outcome_act,stim,rxn,outcome) cell2mat(arrayfun(@(neuron) ...
        cat(2, ...
        nanmean(AP_itril(corr(stim_act(stim > 0 & outcome == 1,:,neuron)'),-1)), ...
        nanmean(AP_itril(corr(move_act(stim > 0 & outcome == 1,:,neuron)'),-1)), ...
        nanmean(AP_itril(corr(outcome_act(stim > 0 & outcome == 1,:,neuron)'),-1))), ...
        transpose(1:size(stim_act,3)),'uni',false)), ...
        mua_allcat_stimalign_exp_smooth,mua_allcat_movealign_exp_smooth,mua_allcat_outcomealign_exp_smooth, ...
        trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false));

[trial_corr_max,trial_corr_max_align] = max(trial_corr,[],2);



%% ~~~~~~~~~~~~ Maps

%% Plot average maps by domain/celltype

% Plot average within each domain AND cell type
figure;
for curr_depth = 1:n_domains
    for curr_celltype = 1:n_celltypes
        subplot(n_domains,n_celltypes, ...
            (curr_depth-1)*n_celltypes+curr_celltype)
        imagesc(nanmean(ctx_str_k_px(:,:,domain_aligned_allcat == curr_depth & ...
            celltype_allcat == curr_celltype),3));
        axis image
        colormap(brewermap([],'*RdBu'));
        caxis([-max(abs(caxis)),max(abs(caxis))]);
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        axis off
        title(['Str ' num2str(curr_depth) ' ' celltype_labels{curr_celltype}]);
    end
end

%% Plot maps by condition 1/2

curr_domain = 2;

% max-normalize maps?
ctx_str_k_px_norm = ctx_str_k_px./max(max(ctx_str_k_px,[],1),[],2);

figure;

for curr_celltype = 1:n_celltypes
    
    curr_k_1 = nanmean(ctx_str_k_px_norm(:,:, ...
        stim_cells & ...
        domain_aligned_allcat == curr_domain & ...
        celltype_allcat == curr_celltype),3);
    
     curr_k_2 = nanmean(ctx_str_k_px_norm(:,:, ...
        move_cells & ...
        domain_aligned_allcat == curr_domain & ...
        celltype_allcat == curr_celltype),3);
    
    subplot(3,n_celltypes,0*n_celltypes+curr_celltype)
    imagesc(curr_k_1);
    axis image
    colormap(brewermap([],'*RdBu'));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
%     caxis([-2e-4,2e-4]);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    axis off
    title(['Cond 1: Str ' num2str(curr_domain) ' ' celltype_labels{curr_celltype}]);
    
    subplot(3,n_celltypes,1*n_celltypes+curr_celltype)
    imagesc(curr_k_2);
    axis image
    colormap(brewermap([],'*RdBu'));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
%     caxis([-2e-4,2e-4]);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    axis off
    title(['Cond 2: Str ' num2str(curr_domain) ' ' celltype_labels{curr_celltype}]);
    
    subplot(3,n_celltypes,2*n_celltypes+curr_celltype)
    imagesc(curr_k_1 - curr_k_2);
    axis image
    colormap(brewermap([],'*RdBu'));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
%     caxis([-2e-4,2e-4]);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    axis off
    title(['Diff: Str ' num2str(curr_domain) ' ' celltype_labels{curr_celltype}]);
    
end

%% >>>>>> HERE: regress maps from different time points

% Set regression parameters
regression_params.use_svs = 1:100;
regression_params.kernel_t = [-0.1,0.1];
regression_params.zs = [false,false];
regression_params.cvfold = 2;
regression_params.use_constant = true;
lambda = 200;
kernel_frames = floor(regression_params.kernel_t(1)*sample_rate): ...
    ceil(regression_params.kernel_t(2)*sample_rate);

% Pick time to use for regression
use_t = t < 0;
% use_t = true(size(t));

% Regress cortex to striatum
mua_ctxtrialpred_exp = cellfun(@(x) nan(size(x)),mua_allcat_exp,'uni',false);
mua_ctxtrialpred_k = cell(size(mua_ctxtrialpred_exp));
for curr_exp = 1:length(trials_recording)
        
        n_cells = size(mua_allcat_exp{curr_exp},3);
        n_vs = size(fluor_deconv_allcat_exp{curr_exp},3);
    
        curr_mua = reshape(permute(mua_allcat_exp{curr_exp}(:,use_t,:), ...
            [2,1,3]),[],n_cells)';
        curr_mua_std = curr_mua./nanstd(curr_mua,[],2);
        
        curr_fluor = reshape(permute( ...
            fluor_deconv_allcat_exp{curr_exp}(:,use_t,:),[2,1,3]),[],n_vs)';
        curr_fluor_full = reshape(permute( ...
            fluor_deconv_allcat_exp{curr_exp}(:,:,:),[2,1,3]),[],n_vs)';

        % Set discontinuities in trial data
        trial_discontinuities = false(size(mua_allcat_exp{curr_exp}(:,use_t,1)));
        trial_discontinuities(:,1) = true;
        trial_discontinuities = reshape(trial_discontinuities',[],1)';
        
        % Full cortex regression
        [k_fluor,curr_mua_fluorpred_std,explained_var_trial] = ...
            AP_regresskernel(curr_fluor(regression_params.use_svs,:), ...
            curr_mua_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant,trial_discontinuities);
        %  NOTE: THIS ISN'T USABLE, NEED CONSTANT ALSO
        curr_mua_fluorpred = curr_mua_fluorpred_std.*nanstd(curr_mua,[],2);
        
        mua_ctxtrialpred_k{curr_exp} = k_fluor;
     
    AP_print_progress_fraction(curr_exp,length(trials_recording));
end

% Get trial-based maps for all cells
% (use t = 0 kernels for all cells, hardcoded at the moment)
use_k_frame = 5;
ctx_str_trialk_frame = cell2mat(cellfun(@(x) reshape(x(:,use_k_frame,:),100,[]), ...
    mua_ctxtrialpred_k','uni',false));
ctx_str_trialk_px = svdFrameReconstruct(U_master(:,:,1:100),ctx_str_trialk_frame);

% Max-normalize this for now
ctx_str_trialk_px = ctx_str_trialk_px./max(max(ctx_str_trialk_px,[],1),[],2);


%% Plot trial-regressed maps by condition 1/2

curr_domain = 2;

figure;

for curr_celltype = 1:n_celltypes
    
    curr_k_1 = nanmean(ctx_str_trialk_px(:,:, ...
        stim_cells & ...
        domain_aligned_allcat == curr_domain & ...
        celltype_allcat == curr_celltype),3);
    
     curr_k_2 = nanmean(ctx_str_trialk_px(:,:, ...
        move_cells & ...
        domain_aligned_allcat == curr_domain & ...
        celltype_allcat == curr_celltype),3);
    
    subplot(3,n_celltypes,0*n_celltypes+curr_celltype)
    imagesc(curr_k_1);
    axis image
    colormap(brewermap([],'*RdBu'));
    caxis([-0.5,0.5])
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    axis off
    title(['Cond 1: Str ' num2str(curr_domain) ' ' celltype_labels{curr_celltype}]);
    
    subplot(3,n_celltypes,1*n_celltypes+curr_celltype)
    imagesc(curr_k_2);
    axis image
    colormap(brewermap([],'*RdBu'));
    caxis([-0.5,0.5]);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    axis off
    title(['Cond 2: Str ' num2str(curr_domain) ' ' celltype_labels{curr_celltype}]);
    
    subplot(3,n_celltypes,2*n_celltypes+curr_celltype)
    imagesc(curr_k_1 - curr_k_2);
    axis image
    colormap(brewermap([],'*RdBu'));
    caxis([-0.5,0.5]);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    axis off
    title(['Diff: Str ' num2str(curr_domain) ' ' celltype_labels{curr_celltype}]);
    
end



%% ~~~~~~~~~~~~ Regressed activity

%% Get explained total variance for all cells

% Get data by experiment
mua_allcat_exp = vertcat(mua_all{:});
mua_taskpred_allcat_exp = vertcat(mua_taskpred_all{:});
mua_taskpred_reduced_allcat_exp = vertcat(mua_taskpred_reduced_all{:});
mua_ctxpred_allcat_exp = vertcat(mua_ctxpred_all{:});

% Get R^2 for task, cortex full, and cortex ROI predictions
taskpred_r2 = cell(size(mua_all));
taskpred_partial_r2 = cell(size(mua_all));
ctxpred_r2 = cell(size(mua_all));

for curr_exp = 1:length(mua_allcat_exp)
       
    curr_data = ...
        reshape(permute(mua_allcat_exp{curr_exp},[2,1,3]),[],size(mua_allcat_exp{curr_exp},3));
    curr_data_baselinesub = ...
        reshape(permute(mua_allcat_exp{curr_exp},[2,1,3]),[],size(mua_allcat_exp{curr_exp},3)) - ...
        nanmean(reshape(mua_allcat_exp{curr_exp}(:,t < 0,:),[],size(mua_allcat_exp{curr_exp},3)),1);
    
    curr_taskpred_data = ...
        reshape(permute(mua_taskpred_allcat_exp{curr_exp},[2,1,3]),[],size(mua_allcat_exp{curr_exp},3));
    curr_taskpred_reduced_data = ...
        reshape(permute(mua_taskpred_reduced_allcat_exp{curr_exp},[2,1,3,4]),[], ...
        size(mua_allcat_exp{curr_exp},3),size(mua_taskpred_reduced_allcat_exp{curr_exp},4));
    curr_ctxpred_data = ...
        reshape(permute(mua_ctxpred_allcat_exp{curr_exp},[2,1,3]),[],size(mua_allcat_exp{curr_exp},3));
    
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_data_baselinesub) | ...
        isnan(curr_taskpred_data) | isnan(curr_ctxpred_data);
    curr_data(nan_samples) = NaN;
    curr_data_baselinesub(nan_samples) = NaN;
    curr_taskpred_data(nan_samples) = NaN;
    curr_taskpred_reduced_data(repmat(nan_samples,1,1,size(curr_taskpred_reduced_data,3))) = NaN;
    curr_ctxpred_data(nan_samples) = NaN;

    % (note: task regression done on baseline subtracted)
    taskpred_r2{curr_exp} = (1 - (nansum((curr_data_baselinesub-curr_taskpred_data).^2,1)./ ...
        nansum((curr_data_baselinesub-nanmean(curr_data_baselinesub,1)).^2,1)))';
    taskpred_partial_r2{curr_exp} = permute((1 - (nansum((curr_data_baselinesub-curr_taskpred_data).^2,1)./ ...
        nansum((curr_data_baselinesub-curr_taskpred_reduced_data).^2,1))),[2,3,1]);
        
    % (cortex regression done on raw data)
    ctxpred_r2{curr_exp} = (1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1)))';

end

taskpred_r2_allcat = cell2mat(taskpred_r2);
taskpred_partial_r2_allcat = cell2mat(taskpred_partial_r2);
ctxpred_r2_allcat = cell2mat(ctxpred_r2);

taskpred_r2_allcat(isinf(taskpred_r2_allcat)) = NaN;
taskpred_partial_r2_allcat(isinf(taskpred_partial_r2_allcat)) = NaN;
ctxpred_r2_allcat(isinf(ctxpred_r2_allcat)) = NaN;

figure;plot(taskpred_r2_allcat,ctxpred_r2_allcat,'.k')
xlim([-0.1,1])
ylim([-0.1,1])
line(xlim,xlim);
xlabel('Task R^2');
ylabel('Cortex R^2');


%% Get explained trial-averaged variance for all cells

warning('Weird use of time: > x% trials not nan');

% Get data by experiment
mua_exp = vertcat(mua_all{:});
mua_taskpred_exp = vertcat(mua_taskpred_all{:});
mua_ctxpred_exp = vertcat(mua_ctxpred_all{:});

% Get R^2 for task, cortex full, and cortex ROI predictions
taskpred_trialavg_r2 = cell(size(mua_all));
ctxpred_trialavg_r2 = cell(size(mua_all));

for curr_exp = 1:length(mua_allcat_exp)
       
    % Current data
    curr_data = mua_exp{curr_exp};
    curr_data_baselinesub = mua_exp{curr_exp} - ...
        nanmean(reshape(mua_exp{curr_exp}(:,t < 0,:),[],1,size(mua_exp{curr_exp},3)),1);
    curr_taskpred_data = mua_taskpred_exp{curr_exp};
    curr_ctxpred_data = mua_ctxpred_exp{curr_exp};
    
    % Set common nans
    % NOTE: weird at the moment, using time points with > 50% of trials
    use_trials = all(~all(isnan(curr_taskpred_data),2) & ~all(isnan(curr_ctxpred_data),2),3);
    nan_frac_thresh = 0.5; % if frac trials nan, don't use timepoint
    nan_t = any(nanmean(isnan(curr_taskpred_data(use_trials,:,:)),1) > nan_frac_thresh | ...
        any(isnan(curr_ctxpred_data(use_trials,:,:)),1) > nan_frac_thresh,3);
    
    curr_data(:,nan_t,:) = NaN;
    curr_data_baselinesub(:,nan_t,:) = NaN;
    curr_taskpred_data(:,nan_t,:) = NaN;
    curr_ctxpred_data(:,nan_t,:) = NaN;
    
    % Trial-average data
    curr_data_avg = permute(nanmean(curr_data(use_trials,:,:),1),[3,2,1]);
    curr_data_baselinesub_avg = permute(nanmean(curr_data_baselinesub(use_trials,:,:),1),[3,2,1]);
    curr_taskpred_data_avg = permute(nanmean(curr_taskpred_data(use_trials,:,:),1),[3,2,1]);
    curr_ctxpred_data_avg = permute(nanmean(curr_ctxpred_data(use_trials,:,:),1),[3,2,1]);
    
    % Trial-average R^2
    % (note: task regression done on baseline subtracted)
    taskpred_trialavg_r2{curr_exp} = (1 - (nansum((curr_data_baselinesub_avg-curr_taskpred_data_avg).^2,2)./ ...
        nansum((curr_data_baselinesub_avg-nanmean(curr_data_baselinesub_avg,2)).^2,2)));
    % (cortex regression done on raw data)
    ctxpred_trialavg_r2{curr_exp} = (1 - (nansum((curr_data_avg-curr_ctxpred_data_avg).^2,2)./ ...
        nansum((curr_data_avg-nanmean(curr_data_avg,2)).^2,2)));

end

taskpred_trialavg_r2_allcat = cell2mat(taskpred_trialavg_r2);
ctxpred_trialavg_r2_allcat = cell2mat(ctxpred_trialavg_r2);

figure;plot(taskpred_trialavg_r2_allcat,ctxpred_trialavg_r2_allcat,'.k')
xlim([-0.1,1])
ylim([-0.1,1])
line(xlim,xlim);
xlabel('Task R^2');
ylabel('Cortex R^2');


%% Get partial explained trial-averaged variance for all cells

warning('Weird use of time: > x% trials not nan');

% Get data by experiment
mua_exp = vertcat(mua_all{:});
mua_taskpred_exp = vertcat(mua_taskpred_all{:});
mua_taskpred_reduced_exp = vertcat(mua_taskpred_reduced_all{:});

% Get R^2 for task, cortex full, and cortex ROI predictions
taskpred_trialavg_partial_r2 = cell(size(mua_all));

for curr_exp = 1:length(mua_allcat_exp)
       
    % Current data
    curr_data_baselinesub = mua_exp{curr_exp} - ...
        nanmean(reshape(mua_exp{curr_exp}(:,t < 0,:),[],1,size(mua_exp{curr_exp},3)),1);
    curr_taskpred_data = mua_taskpred_exp{curr_exp};
    curr_taskpred_reduced_data = mua_taskpred_reduced_exp{curr_exp};
    
    % Set common nans
    % NOTE: weird at the moment, using time points with > 50% of trials
    use_trials = all(~all(isnan(curr_taskpred_data),2),3);
    nan_frac_thresh = 0.5; % if frac trials nan, don't use timepoint
    nan_t = any(nanmean(isnan(curr_taskpred_data(use_trials,:,:)),1) > nan_frac_thresh,3);
    
    curr_data_baselinesub(:,nan_t,:) = NaN;
    curr_taskpred_data(:,nan_t,:) = NaN;
    curr_taskpred_reduced_data(:,nan_t,:,:) = NaN;
    
    % Trial-average data
    curr_data_baselinesub_avg = permute(nanmean(curr_data_baselinesub(use_trials,:,:),1),[3,2,1]);
    curr_taskpred_data_avg = permute(nanmean(curr_taskpred_data(use_trials,:,:),1),[3,2,1]);
    curr_taskpred_reduced_data_avg = permute(nanmean(curr_taskpred_reduced_data(use_trials,:,:,:),1),[3,2,4,1]);
    
    % Trial-average R^2
    % partial variance = 1-(sse residual full model / sse residual reduced model)
    curr_sse_full = nansum((curr_data_baselinesub_avg - curr_taskpred_data_avg).^2,2);
    curr_sse_reduced = permute(nansum((curr_data_baselinesub_avg - ...
        curr_taskpred_reduced_data_avg).^2,2),[1,3,2]);
    
    taskpred_trialavg_partial_r2{curr_exp} = 1 - (curr_sse_full./curr_sse_reduced);

end

taskpred_trialavg_partial_r2_allcat = cell2mat(taskpred_trialavg_partial_r2);

figure;plot3(taskpred_trialavg_partial_r2_allcat(:,1), ...
    taskpred_trialavg_partial_r2_allcat(:,2), ...
    taskpred_trialavg_partial_r2_allcat(:,4),'.k')
xlabel('Stim R^2');
ylabel('Move R^2');
zlabel('Reward R^2');
xlim([-0.1,1]);
ylim(xlim);zlim(xlim);
axis vis3d;



%% Explained variance by type/task-responsiveness

figure;
curr_domain = 1;
for curr_celltype = 1:n_celltypes
    
    curr_ctxpred_r2_stim = ctxpred_r2_allcat(stim_cells & ...
        domain_aligned_allcat == curr_domain & celltype_allcat == curr_celltype);
    
    curr_ctxpred_r2_notstim = ctxpred_r2_allcat(~stim_cells & ...
        domain_aligned_allcat == curr_domain & celltype_allcat == curr_celltype);
    
    curr_taskpred_r2_stim = taskpred_r2_allcat(stim_cells & ...
        domain_aligned_allcat == curr_domain & celltype_allcat == curr_celltype);
    
    curr_taskpred_r2_notstim = taskpred_r2_allcat(~stim_cells & ...
        domain_aligned_allcat == curr_domain & celltype_allcat == curr_celltype);
    
    curr_plot_data = ...
        {curr_ctxpred_r2_stim;curr_ctxpred_r2_notstim; ...
        curr_taskpred_r2_stim;curr_taskpred_r2_notstim};
    curr_plot_grp = cellfun(@(x,grp) repmat(grp,length(x),1), ...
        curr_plot_data,num2cell(1:length(curr_plot_data))','uni',false);
    
    subplot(1,n_celltypes, ...
        (curr_domain-1)*n_celltypes+curr_celltype)
    
    boxplot(vertcat(curr_plot_data{:}),vertcat(curr_plot_grp{:}));
    set(gca,'XTickLabel',{'Stim: ctx','Not stim: ctx','Stim: task','Not stim: task'});
    set(gca,'XTickLabelRotation',45);
    ylabel('Trial-average R^2');
    title([num2str(curr_domain) ' ' celltype_labels{curr_celltype}]);
    
end

%% Explained trial-average variance by type/task-responsiveness

figure;
curr_domain = 1;
for curr_celltype = 1:n_celltypes
    
    curr_ctxpred_r2_stim = ctxpred_trialavg_r2_allcat(stim_cells & ...
        domain_aligned_allcat == curr_domain & celltype_allcat == curr_celltype);
    
    curr_ctxpred_r2_notstim = ctxpred_trialavg_r2_allcat(~stim_cells & ...
        domain_aligned_allcat == curr_domain & celltype_allcat == curr_celltype);
    
    curr_taskpred_r2_stim = taskpred_trialavg_r2_allcat(stim_cells & ...
        domain_aligned_allcat == curr_domain & celltype_allcat == curr_celltype);
    
    curr_taskpred_r2_notstim = taskpred_trialavg_r2_allcat(~stim_cells & ...
        domain_aligned_allcat == curr_domain & celltype_allcat == curr_celltype);
    
    curr_plot_data = ...
        {curr_ctxpred_r2_stim;curr_ctxpred_r2_notstim; ...
        curr_taskpred_r2_stim;curr_taskpred_r2_notstim};
    curr_plot_grp = cellfun(@(x,grp) repmat(grp,length(x),1), ...
        curr_plot_data,num2cell(1:length(curr_plot_data))','uni',false);
    
    subplot(1,n_celltypes, ...
        (curr_domain-1)*n_celltypes+curr_celltype)
    
    boxplot(vertcat(curr_plot_data{:}),vertcat(curr_plot_grp{:}));
    set(gca,'XTickLabel',{'Stim: ctx','Not stim: ctx','Stim: task','Not stim: task'});
    set(gca,'XTickLabelRotation',45);
    ylabel('R^2');
    title([num2str(curr_domain) ' ' celltype_labels{curr_celltype}]);
    
end

%% Task kernels of domain/celltype

% Concatenate task kernels
n_regressors = length(task_regressor_labels);
mua_taskpred_k_allcat = arrayfun(@(x) cell2mat(permute(...
    cellfun(@(k) k(x),vertcat(mua_taskpred_k_all{:})),[2,3,1])),1:4,'uni',false);
mua_ctxpred_taskpred_k_allcat = arrayfun(@(x) cell2mat(permute(...
    cellfun(@(k) k(x),vertcat(mua_ctxpred_taskpred_k_all{:})),[2,3,1])),1:4,'uni',false);

% Set kernel colors and times
stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = [];
move_col = [0.6,0,0.6;0,0.6,0];
go_col = [0.8,0.8,0.2;0.5,0.5,0.5];
outcome_col = [0,0,0.8;0.8,0,0];
task_regressor_cols = {stim_col,move_col,go_col,outcome_col};
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

% Get baseline of each cell to normalize
t_baseline = t < 0.5;
mua_baseline = ...
    cell2mat(cellfun(@(x) nanmean(reshape(x(:,t_baseline,:),[],size(x,3)),1)', ...
    mua_allcat_stimalign_exp,'uni',false));

% Plot each regressor by domain/celltype
for curr_regressor = 1:n_regressors
    
    figure('Name',task_regressor_labels{curr_regressor});
    for curr_domain = 1:n_domains
        for curr_celltype = 1:n_celltypes
            
            curr_cells = domain_aligned_allcat == curr_domain & ...
                celltype_allcat == curr_celltype;
            
            subplot(n_domains,n_celltypes, ...
                (curr_domain-1)*n_celltypes+curr_celltype);
            hold on;
            set(gca,'ColorOrder',task_regressor_cols{curr_regressor});
            
            softnorm = 0.2;
            curr_k = mua_taskpred_k_allcat{curr_regressor}(:,:,curr_cells)./ ...
                (softnorm+permute(mua_baseline(curr_cells),[2,3,1]));
            
            plot(task_regressor_t_shifts{curr_regressor}, ...
                nanmean(curr_k,3)', ...
                'linewidth',2);
            
            title([num2str(curr_domain) ' ' celltype_labels{curr_celltype}]);
            
        end
    end
    linkaxes(get(gcf,'Children'),'xy');
    
end





%% Task kernels of select cells

% Plot kernels from cell groups
plot_regressor = 1;
plot_cells_1 = stim_cells & domain_aligned_allcat == 1  & celltype_allcat == 1;
plot_cells_2 = ~stim_cells & domain_aligned_allcat == 1  & celltype_allcat == 1;
% plot_cells_1 = stim_cells & celltype_allcat == 1;
% plot_cells_2 = ~stim_cells & celltype_allcat == 1;

% Concatenate task kernels
n_regressors = length(task_regressor_labels);
mua_taskpred_k_allcat = arrayfun(@(x) cell2mat(permute(...
    cellfun(@(k) k(x),vertcat(mua_taskpred_k_all{:})),[2,3,1])),1:4,'uni',false);
mua_ctxpred_taskpred_k_allcat = arrayfun(@(x) cell2mat(permute(...
    cellfun(@(k) k(x),vertcat(mua_ctxpred_taskpred_k_all{:})),[2,3,1])),1:4,'uni',false);

% Set kernel colors and times
stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = [];
move_col = [0.6,0,0.6;0,0.6,0];
go_col = [0.8,0.8,0.2;0.5,0.5,0.5];
outcome_col = [0,0,0.8;0.8,0,0];
task_regressor_cols = {stim_col,move_col,go_col,outcome_col};
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

% Plot average within groups
figure; 
subplot(2,2,1); hold on;
set(gca,'ColorOrder',task_regressor_cols{plot_regressor});
plot(task_regressor_t_shifts{plot_regressor}, ...
    nanmean(mua_taskpred_k_allcat{plot_regressor}(:,:,plot_cells_1),3)');
title('Group 1');

subplot(2,2,2); hold on;
set(gca,'ColorOrder',task_regressor_cols{plot_regressor});
plot(task_regressor_t_shifts{plot_regressor}, ...
    nanmean(mua_taskpred_k_allcat{plot_regressor}(:,:,plot_cells_2),3)');
title('Group 2');

subplot(2,2,3); hold on;
set(gca,'ColorOrder',task_regressor_cols{plot_regressor});
plot(task_regressor_t_shifts{plot_regressor}, ...
    nanmean(mua_ctxpred_taskpred_k_allcat{plot_regressor}(:,:,plot_cells_1),3)');
title('Group 1 ctxpred');

subplot(2,2,4); hold on;
set(gca,'ColorOrder',task_regressor_cols{plot_regressor});
plot(task_regressor_t_shifts{plot_regressor}, ...
    nanmean(mua_ctxpred_taskpred_k_allcat{plot_regressor}(:,:,plot_cells_2),3)');
title('Group 2 ctxpred');

linkaxes(get(gcf,'Children'),'xy');

























