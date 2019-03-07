% AP_load_concat_normalize_ctx_str(data_fn,exclude_data)
%
% (not a function any more, used to be this)
% [t,fluor_allcat_deriv,fluor_roi_deriv,mua_allcat,wheel_allcat,reward_allcat,D_allcat] = ...
%     AP_load_concat_normalize_ctx_str(data_fn,exclude_data)
%
% Loads, concatenates, and normalizes trial data for the ctx-str project
% data_fn - filename (e.g. trial_activity_choiceworld)
% exclude_data - true/false at the moment (for bad behavior days)
% (can in future extend this to wonky recording days?)

% Load data into structure
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
data_struct = load([data_path filesep data_fn]);

% Find fields with experiment data (cell arrays with length animals)
data_struct_fieldnames = fieldnames(data_struct);
experiment_fields = cellfun(@(curr_field) ...
    length(data_struct.(curr_field)) == length(data_struct.animals) && ...
    all(cellfun(@(x) iscell(x),data_struct.(curr_field))),data_struct_fieldnames);

% Load pre-marked experiments to exclude and cut out bad ones
if exist('exclude_data','var') && exclude_data
    exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
    exclude_fn = 'bhv_use_experiments';
    load([exclude_path filesep exclude_fn]);
    
    % Pull out used experiments for the animals loaded
    use_experiments_animals = ismember(data_struct.animals,{bhv_use_experiments.animals});
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
            data_struct.(cell2mat(curr_field)) = cellfun(@(x,expts) ...
                x(expts),data_struct.(cell2mat(curr_field)), ...
                use_experiments,'uni',false);
        end
    end
end

% Unpack data structure into workspace
arrayfun(@(x) assignin('base',cell2mat(x),data_struct.(cell2mat(x))),data_struct_fieldnames);
clear data_struct

% Get if this is a task dataset
task_dataset = exist('outcome_all','var');

% Get number of widefield components and MUA depths
n_vs = size(fluor_all{end}{end},3);
n_depths = size(mua_all{end}{end},3);

% Concatenate behavior data
D_fields = fieldnames(D_all{end}{end});
D_allcat = cell2struct(arrayfun(@(curr_field) ...
    cell2mat(cellfun(@(x) x(curr_field), ...
    cellfun(@struct2cell,vertcat(D_all{:}),'uni',false))), ...
    1:length(D_fields),'uni',false),D_fields,2);

wheel_allcat = cell2mat(vertcat(wheel_all{:}));
if task_dataset
    outcome_allcat = cell2mat(vertcat(outcome_all{:}));
end

% (movement from mousecam if exists: normalize and concatenate)
if exist('movement_all','var')
    movement_all_norm = cellfun(@(x) cellfun(@(x) x./nanstd(x(:)), ...
        x,'uni',false),movement_all,'uni',false);
    movement_allcat = cell2mat(vertcat(movement_all_norm{:}));
end

% Concatenate cortex data
fluor_allcat = cell2mat(vertcat(fluor_all{:}));
if task_dataset
    fluor_taskpred_allcat = cell2mat(vertcat(fluor_taskpred_all{:}));
    fluor_taskpred_reduced_allcat = cell2mat(vertcat(fluor_taskpred_reduced_all{:}));
end

% Normalize MUA by experiment and concatenate
t_baseline = t < 0;
mua_day_baseline = cellfun(@(mua_animal) cellfun(@(mua_day) ...
    nanmean(reshape(mua_day(:,t_baseline,:),[],1,n_depths)), ...
    mua_animal,'uni',false),mua_all,'uni',false);
mua_day_std = cellfun(@(mua_animal) cellfun(@(mua_day) ...
    nanstd(reshape(mua_day(:,t_baseline,:),[],1,n_depths)), ...
    mua_animal,'uni',false),mua_all,'uni',false);

% Choose what to normalize MUA by (std or baseline)
mua_norm = mua_day_baseline;

% (NaN-out any trials that are all zeros)
mua_nan_trials = +any(cell2mat(vertcat(mua_all{:})),2);
mua_nan_trials(~mua_nan_trials) = NaN;

mua_allcat = cell2mat(cellfun(@(mua,mua_baseline,mua_std) ...
    (mua-mua_baseline)./mua_std,vertcat(mua_all{:}), ...
    vertcat(mua_day_baseline{:}),vertcat(mua_norm{:}),'uni',false)) ...
    .*mua_nan_trials;

mua_ctxpred_allcat = cell2mat(cellfun(@(mua,mua_baseline,mua_std) ...
    (mua-mua_baseline)./mua_std,vertcat(mua_ctxpred_all{:}), ...
    vertcat(mua_day_baseline{:}),vertcat(mua_norm{:}),'uni',false)) ...
    .*mua_nan_trials;

if task_dataset
    % Don't subtract baseline for task-predicted: done before regression
    mua_taskpred_allcat = cell2mat(cellfun(@(mua,mua_baseline,mua_std) ...
        (mua)./mua_std,vertcat(mua_taskpred_all{:}), ...
        vertcat(mua_day_baseline{:}),vertcat(mua_norm{:}),'uni',false)) ...
        .*mua_nan_trials;
    mua_taskpred_reduced_allcat = cell2mat(cellfun(@(mua,mua_baseline,mua_std) ...
        (mua)./mua_std,vertcat(mua_taskpred_reduced_all{:}), ...
        vertcat(mua_day_baseline{:}),vertcat(mua_norm{:}),'uni',false)) ...
        .*mua_nan_trials;
    
    mua_ctxpred_taskpred_allcat = cell2mat(cellfun(@(mua,mua_baseline,mua_std) ...
        (mua)./mua_std,vertcat(mua_ctxpred_taskpred_all{:}), ...
        vertcat(mua_day_baseline{:}),vertcat(mua_norm{:}),'uni',false)) ...
        .*mua_nan_trials;
    mua_ctxpred_taskpred_reduced_allcat = cell2mat(cellfun(@(mua,mua_baseline,mua_std) ...
        (mua)./mua_std,vertcat(mua_ctxpred_taskpred_reduced_all{:}), ...
        vertcat(mua_day_baseline{:}),vertcat(mua_norm{:}),'uni',false)) ...
        .*mua_nan_trials;
end

% Deconvolve fluorescence, subtract baseline
fluor_allcat_deconv = AP_deconv_wf(fluor_allcat);
fluor_allcat_deconv_baseline = nanmean(reshape(fluor_allcat_deconv(:,t_baseline,:),[],1,n_vs));
fluor_allcat_deconv = fluor_allcat_deconv - fluor_allcat_deconv_baseline;

% Get fluorescence ROIs
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1); % (left ROIs only)
n_rois = numel(wf_roi);

fluor_roi_deconv = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
    n_rois,[],size(fluor_allcat_deconv,1)),[3,2,1]);

if task_dataset
    fluor_roi_taskpred = permute(reshape( ...
        AP_svd_roi(U_master(:,:,1:n_vs), ...
        reshape(permute(fluor_taskpred_allcat,[3,2,1]),n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
        n_rois,[],size(fluor_taskpred_allcat,1)),[3,2,1]);
    
    fluor_roi_taskpred_reduced = cell2mat(permute(arrayfun(@(x) ...
        permute(reshape( ...
        AP_svd_roi(U_master(:,:,1:n_vs), ...
        reshape(permute(fluor_taskpred_reduced_allcat(:,:,:,x),[3,2,1]), ...
        n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
        size(cat(3,wf_roi.mask),3),[],size(fluor_taskpred_reduced_allcat,1)),[3,2,1]), ...
        1:size(fluor_taskpred_reduced_allcat,4),'uni',false),[1,3,4,2]));
end

%%%% OLD: smoothed fluorescence derivative
%
% % Get fluorescence derivative and interpolate to match original t
% deriv_smooth = 3;
%
% t_smooth_diff = conv(conv(t,ones(1,deriv_smooth)/deriv_smooth,'valid'),[1,1]/2,'valid');
% fluor_allcat_deriv = permute(interp1(t_smooth_diff, permute(diff( ...
%     convn(fluor_allcat,ones(1,deriv_smooth)/deriv_smooth,'valid'),[],2), ...
%     [2,1,3]),t,'linear','extrap'),[2,1,3]);
%
% fluor_roi_deriv = permute(reshape( ...
%     AP_svd_roi(U_master(:,:,1:n_vs), ...
%     reshape(permute(fluor_allcat_deriv,[3,2,1]),n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
%     n_rois,[],size(fluor_allcat_deriv,1)),[3,2,1]);

% If loading task data, get task-relevant variables and movement info
if task_dataset
    
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
    
    % Make move-aligned MUA
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
    
    % Make move-aligned fluorescence
    fluor_allcat_deconv_move = fluor_allcat_deconv;
    fluor_taskpred_reduced_allcat_move = fluor_taskpred_reduced_allcat;
    fluor_roi_deconv_move = fluor_roi_deconv;
    fluor_roi_taskpred_move = fluor_roi_taskpred;
    fluor_roi_taskpred_reduced_move = fluor_roi_taskpred_reduced;
    
    t_leeway = -t(1);
    leeway_samples = round(t_leeway*(sample_rate));
    for i = 1:size(mua_allcat,1)
        fluor_allcat_deconv_move(i,:,:,:) = circshift(fluor_allcat_deconv_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
        fluor_taskpred_reduced_allcat_move(i,:,:,:) = circshift(fluor_taskpred_reduced_allcat_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
        fluor_roi_deconv_move(i,:,:,:) = circshift(fluor_roi_deconv_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
        fluor_roi_taskpred_move(i,:,:,:) = circshift(fluor_roi_taskpred_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
        fluor_roi_taskpred_reduced_move(i,:,:,:) = circshift(fluor_roi_taskpred_reduced_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
    end
    
    % (move-align mousecam movement if it exists)
    if exist('movement_allcat','var')
        movement_allcat_move = movement_allcat;
        
        t_leeway = -t(1);
        leeway_samples = round(t_leeway*(sample_rate));
        for i = 1:size(mua_allcat,1)
            movement_allcat_move(i,:,:,:) = circshift(movement_allcat_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
        end
    end
        
end

disp('Finished loading trials')







