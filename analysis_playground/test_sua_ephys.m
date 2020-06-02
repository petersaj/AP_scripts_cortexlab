%% ~~~~~~~~~~~~ Load

%% Load SUA data

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
n_celltypes = max(groups_allcat)./n_aligned_depths;
celltype_allcat = ceil(groups_allcat./n_aligned_depths);
celltype_labels = {'MSN','FSI','TAN','TH','Other'};
if length(celltype_labels) ~= n_celltypes
    error('Mismatching number of celltypes and labels')
end
domain_allcat = mod(groups_allcat,n_aligned_depths) + ...
    n_aligned_depths*(mod(groups_allcat,n_aligned_depths) == 0);

% Get maps for all cells
% (load master U)
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');
% (use t = 0 kernels for all cells, hardcoded at the moment)
use_k_frame = 5;
ctx_str_k_frame = cell2mat(cellfun(@(x) reshape(x(:,use_k_frame,:),100,[]), ...
    vertcat(ctx_str_k_all{:})','uni',false));
ctx_str_k_px = svdFrameReconstruct(U_master(:,:,1:100),ctx_str_k_frame);


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
n_aligned_depths = 3;
ephys_kernel_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
ephys_kernel_align_fn = ['ephys_kernel_align_' num2str(n_aligned_depths) '_depths.mat'];
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
    ', Cell: ' num2str(x) ...
    ', Domain: ' num2str(domain_aligned_allcat(x)), ...
    ', ' celltype_labels{celltype_allcat(x)}],1:length(sort_idx),'uni',false);
AP_image_scroll(mua_allcat_stimalign_exp_pad(:,:,sort_idx),cell_string);
line(repmat(find(t >= 0,1),2,1),ylim,'linestyle','--','color','r');
colormap(brewermap([],'Greys'));
caxis([0,50]);



%% Rasters (select cells) 

% Set cells to plot
plot_cells = find(stim_cells & celltype_allcat == 1 & domain_allcat == 1);
% plot_cells = find(stim_cells);

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
       mua_allcat_movealign_exp{curr_exp}(curr_trial,:) = ...
           circshift(mua_allcat_movealign_exp{curr_exp}(curr_trial,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples);
   end
end

% (outcome aligned)
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
mua_allcat_outcomealign_exp = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_outcomealign_exp)
    for curr_trial = 1:size(mua_allcat_outcomealign_exp{curr_exp},1)
        mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:) = ...
            circshift(mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:), ...
            -outcome_idx_exp{curr_exp}(curr_trial)+leeway_samples);
    end
end


% Pad and plot rasters
plot_trials_1 = cellfun(@(stim,rxn,outcome) ...
    stim > 0, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false);
plot_trials_2 = cellfun(@(stim,rxn,outcome) ...
    stim < 0, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false);

max_trials = max(cellfun(@(x,y) sum(x)+sum(y),plot_trials_1,plot_trials_2));
act_exp_pad = cell2mat(permute(cellfun(@(act,trials_1,trials_2) ...
    padarray(act([find(trials_1);find(trials_2)],:,:), ...
    [max_trials-(sum(trials_1)+sum(trials_2)),0,0],NaN,'post'), ...
    mua_allcat_stimalign_exp,plot_trials_1,plot_trials_2,'uni',false),[2,3,1]));

% act_exp_pad = cell2mat(permute(cellfun(@(act_stim,act_move,act_outcome,trials_1,trials_2) ...
%     padarray( ...
%     [act_stim([find(trials_1);find(trials_2)],:,:), ...
%     act_move([find(trials_1);find(trials_2)],:,:), ...
%     act_outcome([find(trials_1);find(trials_2)],:,:)], ...
%     [max_trials-(sum(trials_1)+sum(trials_2)),0,0],NaN,'post'), ...
%     mua_allcat_stimalign_exp,mua_allcat_movealign_exp,mua_allcat_outcomealign_exp, ...
%     plot_trials_1,plot_trials_2,'uni',false),[2,3,1]));

cell_string = arrayfun(@(x) ...
    [' Sort idx: ' num2str(x) ...
    ', Cell: ' num2str(plot_cells(x)) ...
    ', Domain: ' num2str(domain_aligned_allcat(plot_cells(x))), ...
    ', ' celltype_labels{celltype_allcat(plot_cells(x))}],1:length(plot_cells),'uni',false);
AP_image_scroll(act_exp_pad(:,:,plot_cells),cell_string);
line(repmat(find(t >= 0,1),2,1),ylim,'linestyle','--','color','r');
line(repmat(find(t >= 0,1),2,1)+length(t),ylim,'linestyle','--','color','r');
line(repmat(find(t >= 0,1),2,1)+length(t)*2,ylim,'linestyle','--','color','r');

colormap(brewermap([],'Greys'));
caxis([0,50]);



%% Stim: sort and plot 

trial_stim_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
move_t_exp = mat2cell(move_t,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

% sort half, plot other half (early/late move)

% (stim aligned)
mua_allcat_stimalign_exp = vertcat(mua_all{:});

% (move aligned)
move_idx_exp = mat2cell(move_idx,use_split,1);
mua_allcat_movealign_exp = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_movealign_exp)
   for curr_trial = 1:size(mua_allcat_movealign_exp{curr_exp},1)
       mua_allcat_movealign_exp{curr_exp}(curr_trial,:) = ...
           circshift(mua_allcat_movealign_exp{curr_exp}(curr_trial,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples);
   end
end

trial_rand = rand(1000,1);

b1 = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1 & ...
    trial_rand(1:size(act,1)) < 0.5,:,:),1)),mua_allcat_stimalign_exp, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false)')';

b2 = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim > 0 & rxn < 0.5 & ...
    trial_rand(1:size(act,1)) >= 0.5,:,:),1)),mua_allcat_stimalign_exp, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false)')';

b3 = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim < 0 & rxn < 0.5 & ...
    trial_rand(1:size(act,1)) >= 0.5,:,:),1)),mua_allcat_stimalign_exp, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false)')';

c1 = AP_deconv_wf(b1,true);
c2 = AP_deconv_wf(b2,true);
c3 = AP_deconv_wf(b3,true);

[~,max_idx] = max(c1,[],2);
[~,sort_idx] = sort(max_idx);

figure;
subplot(1,4,1);
imagesc(t,[],zscore(c1(sort_idx,:),[],2));
caxis([-3,3])
colormap(brewermap([],'*RdBu'));
title('Sorted half');
subplot(1,4,2);
imagesc(t,[],zscore(c2(sort_idx,:),[],2));
caxis([-3,3])
colormap(brewermap([],'*RdBu'));
title('Unsorted half (early move)');
subplot(1,4,3);
imagesc(t,[],zscore(c3(sort_idx,:),[],2));
caxis([-3,3])
colormap(brewermap([],'*RdBu'));
title('Unsorted half (late move)');
subplot(1,4,4);
imagesc(t,[],imgaussfilt(zscore(c2(sort_idx,:),[],2) - zscore(c3(sort_idx,:),[],2),2));
caxis([-1,1])
colormap(brewermap([],'*RdBu'));
title('Unsorted contra-ipsi (smoothed)');

linkaxes(get(gcf,'Children'),'xy');

% Get stim activity for each cells?
use_t = t > 0 & t < 0.5;
act_diff = zscore(c2-c3,[],2) - zscore(c3,[],2);
stim_act = squeeze(nanmean(act_diff(:,use_t,:),2));


% Pull out cells (by max condition difference > thresh)
use_t = t > 0 & t < 0.2;

stim_cells = ismember(max_idx,find(use_t)) & ...
    max(c2(:,use_t),[],2) > 2*max(c3(:,use_t),[],2);

figure; 
imagesc(zscore([c2(stim_cells,:),c3(stim_cells,:)],[],2));
caxis([-3,3])
colormap(brewermap([],'*RdBu'));
title('Contra/Ipsi stim cells');


%% Stim: sort and plot celltype/domain

plot_cells = celltype_allcat == 1 & domain_allcat == 1;

trial_stim_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
move_t_exp = mat2cell(move_t,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

% sort half, plot other half (early/late move)

% (stim aligned)
mua_allcat_stimalign_exp = vertcat(mua_all{:});

% (move aligned)
move_idx_exp = mat2cell(move_idx,use_split,1);
mua_allcat_movealign_exp = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_movealign_exp)
   for curr_trial = 1:size(mua_allcat_movealign_exp{curr_exp},1)
       mua_allcat_movealign_exp{curr_exp}(curr_trial,:) = ...
           circshift(mua_allcat_movealign_exp{curr_exp}(curr_trial,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples);
   end
end

trial_rand = rand(1000,1);

b1 = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1 & ...
    trial_rand(1:size(act,1)) < 0.5,:,:),1)),mua_allcat_stimalign_exp, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false)')';

b2 = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim > 0 & rxn < 0.5 & ...
    trial_rand(1:size(act,1)) >= 0.5,:,:),1)),mua_allcat_stimalign_exp, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false)')';

b3 = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim < 0 & rxn < 0.5 & ...
    trial_rand(1:size(act,1)) >= 0.5,:,:),1)),mua_allcat_stimalign_exp, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false)')';

c1 = AP_deconv_wf(b1(plot_cells,:),true);
c2 = AP_deconv_wf(b2(plot_cells,:),true);
c3 = AP_deconv_wf(b3(plot_cells,:),true);

[~,max_idx] = max(c1,[],2);
[~,sort_idx] = sort(max_idx);

figure;
subplot(1,4,1);
imagesc(t,[],zscore(c1(sort_idx,:),[],2));
caxis([-3,3])
colormap(brewermap([],'*RdBu'));
title('Sorted half');
subplot(1,4,2);
imagesc(t,[],zscore(c2(sort_idx,:),[],2));
caxis([-3,3])
colormap(brewermap([],'*RdBu'));
title('Unsorted half (early move)');
subplot(1,4,3);
imagesc(t,[],zscore(c3(sort_idx,:),[],2));
caxis([-3,3])
colormap(brewermap([],'*RdBu'));
title('Unsorted half (late move)');
subplot(1,4,4);
imagesc(t,[],imgaussfilt(zscore(c2(sort_idx,:),[],2) - zscore(c3(sort_idx,:),[],2),2));
caxis([-1,1])
colormap(brewermap([],'*RdBu'));
title('Unsorted contra-ipsi (smoothed)');

linkaxes(get(gcf,'Children'),'xy');



%% Move: sort and plot 

trial_stim_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
move_t_exp = mat2cell(move_t,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

% sort half, plot other half (early/late move)

% (stim aligned)
mua_allcat_stimalign_exp = vertcat(mua_all{:});

% (move aligned)
move_idx_exp = mat2cell(move_idx,use_split,1);
mua_allcat_movealign_exp = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_movealign_exp)
   for curr_trial = 1:size(mua_allcat_movealign_exp{curr_exp},1)
       mua_allcat_movealign_exp{curr_exp}(curr_trial,:) = ...
           circshift(mua_allcat_movealign_exp{curr_exp}(curr_trial,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples);
   end
end

trial_rand = rand(1000,1);

b1 = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1 & ...
    trial_rand(1:size(act,1)) < 0.5,:,:),1)),mua_allcat_stimalign_exp, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false)')';

b2 = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim < 0 & rxn < 0.5 & outcome == 1 & ...
    trial_rand(1:size(act,1)) >= 0.5,:,:),1)),mua_allcat_stimalign_exp, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false)')';

b3 = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim < 0 & rxn >= 0.5 & outcome == 1 & ...
    trial_rand(1:size(act,1)) >= 0.5,:,:),1)),mua_allcat_stimalign_exp, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false)')';

c1 = AP_deconv_wf(b1,true);
c2 = AP_deconv_wf(b2,true);
c3 = AP_deconv_wf(b3,true);

[~,max_idx] = max(c1,[],2);
[~,sort_idx] = sort(max_idx);

figure;
subplot(1,4,1);
imagesc(t,[],zscore(c1(sort_idx,:),[],2));
caxis([-3,3])
colormap(brewermap([],'*RdBu'));
title('Sorted half');
subplot(1,4,2);
imagesc(t,[],zscore(c2(sort_idx,:),[],2));
caxis([-3,3])
colormap(brewermap([],'*RdBu'));
title('Unsorted half (early move)');
subplot(1,4,3);
imagesc(t,[],zscore(c3(sort_idx,:),[],2));
caxis([-3,3])
colormap(brewermap([],'*RdBu'));
title('Unsorted half (late move)');
subplot(1,4,4);
imagesc(t,[],imgaussfilt(zscore(c2(sort_idx,:),[],2) - zscore(c3(sort_idx,:),[],2),2));
caxis([-1,1])
colormap(brewermap([],'*RdBu'));
title('Unsorted contra-ipsi (smoothed)');

linkaxes(get(gcf,'Children'),'xy');


use_t = t > 0 & t < 0.5;
act_diff = zscore(c2,[],2) - zscore(c3,[],2);
move_act = squeeze(nanmean(act_diff(:,use_t,:),2));


% Pull out cells (by max condition difference > thresh)
use_t = t > 0 & t < 0.5;

move_cells = ismember(max_idx,find(use_t)) & ...
    max(c2(:,use_t),[],2) > 2*max(c3(:,use_t),[],2);

figure; 
imagesc(zscore([c2(move_cells,:),c3(move_cells,:)],[],2));
caxis([-3,3])
colormap(brewermap([],'*RdBu'));
title('Early/late move cells');


%% Reward: sort and plot 

trial_stim_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
move_t_exp = mat2cell(move_t,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

% sort half, plot other half (early/late move)

% (outcome aligned)
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
mua_allcat_outcomealign_exp = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_outcomealign_exp)
    for curr_trial = 1:size(mua_allcat_outcomealign_exp{curr_exp},1)
        mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:) = ...
            circshift(mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:), ...
            -outcome_idx_exp{curr_exp}(curr_trial)+leeway_samples);
    end
end

trial_rand = rand(1000,1);

b1 = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1 & ...
    trial_rand(1:size(act,1)) < 0.5,:,:),1)),mua_allcat_outcomealign_exp, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false)')';

b2 = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(outcome == 1 & ...
    trial_rand(1:size(act,1)) >= 0.5,:,:),1)),mua_allcat_outcomealign_exp, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false)')';

b3 = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(outcome == -1 & ...
    trial_rand(1:size(act,1)) >= 0.5,:,:),1)),mua_allcat_outcomealign_exp, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false)')';

c1 = AP_deconv_wf(b1,true);
c2 = AP_deconv_wf(b2,true);
c3 = AP_deconv_wf(b3,true);

[~,max_idx] = max(c1,[],2);
[~,sort_idx] = sort(max_idx);

figure;
subplot(1,4,1);
imagesc(t,[],zscore(c1(sort_idx,:),[],2));
caxis([-3,3])
colormap(brewermap([],'*RdBu'));
title('Sorted half');
subplot(1,4,2);
imagesc(t,[],zscore(c2(sort_idx,:),[],2));
caxis([-3,3])
colormap(brewermap([],'*RdBu'));
title('Unsorted half (rewarded)');
subplot(1,4,3);
imagesc(t,[],zscore(c3(sort_idx,:),[],2));
caxis([-3,3])
colormap(brewermap([],'*RdBu'));
title('Unsorted half (unrewarded)');
subplot(1,4,4);
imagesc(t,[],imgaussfilt(zscore(c2(sort_idx,:),[],2) - zscore(c3(sort_idx,:),[],2),2));
caxis([-1,1])
colormap(brewermap([],'*RdBu'));
title('Unsorted contra-ipsi (smoothed)');

linkaxes(get(gcf,'Children'),'xy');


use_t = t > 0.5 & t < 1;
act_diff = zscore(c2,[],2) - zscore(c3,[],2);
reward_act = squeeze(nanmean(act_diff(:,use_t,:),2));


% Pull out cells (by max condition difference > thresh)
use_t = t > 0 & t < 0.5;

reward_cells = ismember(max_idx,find(use_t)) & ...
    max(c2(:,use_t),[],2) > 2*max(c3(:,use_t),[],2);

figure; 
imagesc(zscore([c2(reward_cells,:),c3(reward_cells,:)],[],2));
caxis([-3,3])
colormap(brewermap([],'*RdBu'));
title('Early/late move cells');

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
       mua_allcat_movealign_exp{curr_exp}(curr_trial,:) = ...
           circshift(mua_allcat_movealign_exp{curr_exp}(curr_trial,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples);
   end
end

% (outcome aligned)
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
mua_allcat_outcomealign_exp = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_outcomealign_exp)
    for curr_trial = 1:size(mua_allcat_outcomealign_exp{curr_exp},1)
        mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:) = ...
            circshift(mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:), ...
            -outcome_idx_exp{curr_exp}(curr_trial)+leeway_samples);
    end
end


stim_contra_act = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim > 0 & rxn < 0.5,:,:),1)), ...
    mua_allcat_stimalign_exp,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false)')';

stim_ipsi_act = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim < 0 & rxn < 0.5,:,:),1)), ...
    mua_allcat_stimalign_exp,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false)')';

move_early_act = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim < 0 & rxn < 0.5,:,:),1)), ...
    mua_allcat_stimalign_exp,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false)')';

move_late_act = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim < 0 & rxn > 0.5,:,:),1)), ...
    mua_allcat_stimalign_exp,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false)')';

reward_act = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(outcome == 1,:,:),1)), ...
    mua_allcat_outcomealign_exp,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false)')';

no_reward_act = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(outcome == -1,:,:),1)), ...
    mua_allcat_outcomealign_exp,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false)')';


% Pull out cells (by max condition difference > thresh)
use_t = t > 0 & t < 0.5;

stim_cells = (nanmean(stim_contra_act(:,use_t),2) > 2*nanmean(stim_ipsi_act(:,use_t),2)) & trial_corr_max_align == 1;





%% Count cells by type/domain (after above) 

domain_type_count = ...
    accumarray([domain_aligned_allcat,celltype_allcat, ...
    stim_cells+1,move_cells+1,reward_cells+1],1);

cols = zeros(2,2,2,3);
cols(2,:,:,1) = 1;
cols(:,2,:,2) = 1;
cols(:,:,2,3) = 1;

figure;
for curr_domain = 1:n_aligned_depths
    for curr_type = 1:n_celltypes
        subplot(n_aligned_depths,n_celltypes, ...
            (curr_domain-1)*n_celltypes+curr_type);
        
        curr_counts = squeeze(domain_type_count(curr_domain,curr_type,:,:,:));
        % hack the colors
        curr_counts(curr_counts == 0) = 0.01;
        
        pie(curr_counts)

        title([num2str(curr_domain) ' ' celltype_labels{curr_type}]);
    end
end
legend({'Nothing','Stim','Move','Stim+Move','Reward','Move+Reward','Stim+Move+Reward'})
colormap(hot)

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



%% Stim/move/reward sort and plot

plot_cells = celltype_allcat == 1;

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
       mua_allcat_movealign_exp{curr_exp}(curr_trial,:) = ...
           circshift(mua_allcat_movealign_exp{curr_exp}(curr_trial,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples);
   end
end

% (outcome aligned)
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
mua_allcat_outcomealign_exp = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_outcomealign_exp)
    for curr_trial = 1:size(mua_allcat_outcomealign_exp{curr_exp},1)
        mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:) = ...
            circshift(mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:), ...
            -outcome_idx_exp{curr_exp}(curr_trial)+leeway_samples);
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
title('Sorted trials');

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

%% Get quick trial-trial corr

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
       mua_allcat_movealign_exp{curr_exp}(curr_trial,:) = ...
           circshift(mua_allcat_movealign_exp{curr_exp}(curr_trial,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples);
   end
end

% (outcome aligned)
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
mua_allcat_outcomealign_exp = vertcat(mua_all{:});
for curr_exp = 1:length(mua_allcat_outcomealign_exp)
   for curr_trial = 1:size(mua_allcat_outcomealign_exp{curr_exp},1)
       mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:) = ...
           circshift(mua_allcat_outcomealign_exp{curr_exp}(curr_trial,:), ...
           -outcome_idx_exp{curr_exp}(curr_trial)+leeway_samples);
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


%% ~~~~~~~~~~~~ Maps

%% Plot average maps by domain/celltype

% Plot average within each domain AND cell type
figure;
for curr_depth = 1:n_aligned_depths
    for curr_celltype = 1:n_celltypes
        subplot(n_aligned_depths,n_celltypes, ...
            (curr_depth-1)*n_celltypes+curr_celltype)
        imagesc(nanmean(ctx_str_k_px(:,:,domain_allcat == curr_depth & ...
            celltype_allcat == curr_celltype),3));
        axis image
        colormap(brewermap([],'*RdBu'));
        caxis([-max(abs(caxis)),max(abs(caxis))]);
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        axis off
        title(['Str ' num2str(curr_depth) ' ' celltype_labels{curr_celltype}]);
    end
end

%% Plot maps by stim/no stim

curr_domain = 1;

figure;

for curr_celltype = 1:n_celltypes
    
    curr_k_stim = nanmean(ctx_str_k_px(:,:, ...
        stim_cells & ...
        domain_allcat == curr_domain & ...
        celltype_allcat == curr_celltype),3);
    
     curr_k_notstim = nanmean(ctx_str_k_px(:,:, ...
        ~stim_cells & ...
        domain_allcat == curr_domain & ...
        celltype_allcat == curr_celltype),3);
    
    subplot(3,n_celltypes,0*n_celltypes+curr_celltype)
    imagesc(curr_k_stim);
    axis image
    colormap(brewermap([],'*RdBu'));
    caxis([-2e-4,2e-4]);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    axis off
    title(['Stim: Str ' num2str(curr_domain) ' ' celltype_labels{curr_celltype}]);
    
    subplot(3,n_celltypes,1*n_celltypes+curr_celltype)
    imagesc(curr_k_notstim);
    axis image
    colormap(brewermap([],'*RdBu'));
    caxis([-2e-4,2e-4]);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    axis off
    title(['Not stim: Str ' num2str(curr_domain) ' ' celltype_labels{curr_celltype}]);
    
    subplot(3,n_celltypes,2*n_celltypes+curr_celltype)
    imagesc(curr_k_stim - curr_k_notstim);
    axis image
    colormap(brewermap([],'*RdBu'));
    caxis([-2e-4,2e-4]);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    axis off
    title(['Diff: Str ' num2str(curr_domain) ' ' celltype_labels{curr_celltype}]);
    
end


%% ~~~~~~~~~~~~ Regressed activity

%% Get explained variance for all cells

% Get data by experiment
mua_allcat_exp = vertcat(mua_all{:});
mua_taskpred_allcat_exp = vertcat(mua_taskpred_all{:});
mua_ctxpred_allcat_exp = vertcat(mua_ctxpred_all{:});

% Get R^2 for task, cortex full, and cortex ROI predictions
taskpred_r2 = cell(size(mua_all));
ctxpred_r2 = cell(size(mua_all));

for curr_exp = 1:length(mua_allcat_exp)
       
    curr_data = ...
        reshape(permute(mua_allcat_exp{curr_exp},[2,1,3]),[],size(mua_allcat_exp{curr_exp},3));
    curr_taskpred_data = ...
        reshape(permute(mua_taskpred_allcat_exp{curr_exp},[2,1,3]),[],size(mua_allcat_exp{curr_exp},3));
    curr_ctxpred_data = ...
        reshape(permute(mua_ctxpred_allcat_exp{curr_exp},[2,1,3]),[],size(mua_allcat_exp{curr_exp},3));
    
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_taskpred_data) | isnan(curr_ctxpred_data);
    curr_data(nan_samples) = NaN;
    curr_taskpred_data(nan_samples) = NaN;
    curr_ctxpred_data(nan_samples) = NaN;

    taskpred_r2{curr_exp} = (1 - (nansum((curr_data-curr_taskpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1)))';
    ctxpred_r2{curr_exp} = (1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1)))';

end

taskpred_r2_allcat = cell2mat(taskpred_r2);
ctxpred_r2_allcat = cell2mat(ctxpred_r2);

taskpred_r2_allcat(isinf(taskpred_r2_allcat)) = NaN;
ctxpred_r2_allcat(isinf(ctxpred_r2_allcat)) = NaN;

figure;plot(taskpred_r2_allcat,ctxpred_r2_allcat,'.k')
xlim([-0.1,1])
ylim([-0.1,1])
line(xlim,xlim);
xlabel('Task R^2');
ylabel('Cortex R^2');


%% Explained variance by type/task-responsiveness

figure;
curr_domain = 1;
for curr_celltype = 1:n_celltypes
    
    curr_ctxpred_r2_stim = ctxpred_r2_allcat(stim_cells & ...
        domain_allcat == curr_domain & celltype_allcat == curr_celltype);
    
    curr_ctxpred_r2_notstim = ctxpred_r2_allcat(~stim_cells & ...
        domain_allcat == curr_domain & celltype_allcat == curr_celltype);
    
    curr_taskpred_r2_stim = taskpred_r2_allcat(stim_cells & ...
        domain_allcat == curr_domain & celltype_allcat == curr_celltype);
    
    curr_taskpred_r2_notstim = taskpred_r2_allcat(~stim_cells & ...
        domain_allcat == curr_domain & celltype_allcat == curr_celltype);
    
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


%% Task kernels

% IN PROGRESS

% Concatenate task kernels
n_regressors = length(task_regressor_labels);
mua_taskpred_k_allcat = arrayfun(@(x) cell2mat(permute(...
    cellfun(@(k) k(x),vertcat(mua_taskpred_k_all{:})),[2,3,1])),1:4,'uni',false);
mua_ctxpred_taskpred_k_allcat = arrayfun(@(x) cell2mat(permute(...
    cellfun(@(k) k(x),vertcat(mua_ctxpred_taskpred_k_all{:})),[2,3,1])),1:4,'uni',false);

% Plot task>striatum kernels
stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = [];
move_col = [0.6,0,0.6;0,0.6,0];
go_col = [0.8,0.8,0.2;0.5,0.5,0.5];
outcome_col = [0,0,0.8;0.8,0,0];
task_regressor_cols = {stim_col,move_col,go_col,outcome_col};
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

% Plot kernels from cell groups
plot_regressor = 1;
plot_cells_1 = stim_cells & domain_allcat == 1  & celltype_allcat == 1;
plot_cells_2 = ~stim_cells & domain_allcat == 1  & celltype_allcat == 1;
% plot_cells_1 = stim_cells & celltype_allcat == 1;
% plot_cells_2 = ~stim_cells & celltype_allcat == 1;

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

























