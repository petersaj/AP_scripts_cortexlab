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
    for curr_field = data_struct_fieldnames(experiment_fields)'
        data_struct.(cell2mat(curr_field)) = cellfun(@(x,expts) ...
            x(expts),data_struct.(cell2mat(curr_field)), ...
            use_experiments,'uni',false);
    end   
end

% Unpack data structure into workspace
arrayfun(@(x) assignin('base',cell2mat(x),data_struct.(cell2mat(x))),data_struct_fieldnames);
clear data_struct

% Get number of widefield components and MUA depths
n_vs = size(fluor_all{1}{1},3);
n_depths = size(mua_all{1}{1},3);

% Concatenate data
D_fields = fieldnames(D_all{1}{1});
D_allcat = cell2struct(arrayfun(@(curr_field) ...
    cell2mat(cellfun(@(x) x(curr_field), ...
    cellfun(@struct2cell,vertcat(D_all{:}),'uni',false))), ...
    1:length(D_fields),'uni',false),D_fields,2);

fluor_allcat = cell2mat(vertcat(fluor_all{:}));
wheel_allcat = cell2mat(vertcat(wheel_all{:}));
outcome_allcat = cell2mat(vertcat(outcome_all{:}));

% Normalize MUA by experiment and concatenate
t_baseline = t < 0;
mua_day_baseline = cellfun(@(mua_animal) cellfun(@(mua_day) ...
    nanmean(reshape(mua_day(:,t_baseline,:),[],1,n_depths)),mua_animal,'uni',false),mua_all,'uni',false);
mua_day_std = cellfun(@(mua_animal) cellfun(@(mua_day) ...
    nanstd(reshape(mua_day(:,t_baseline,:),[],1,n_depths)),mua_animal,'uni',false),mua_all,'uni',false);

% (NaN-out any trials that are all zeros)
mua_nan_trials = +any(cell2mat(vertcat(mua_all{:})),2);
mua_nan_trials(~mua_nan_trials) = NaN;

mua_allcat = cell2mat(cellfun(@(mua,mua_baseline,mua_std) ...
    (mua-mua_baseline)./mua_std,vertcat(mua_all{:}), ...
    vertcat(mua_day_baseline{:}),vertcat(mua_day_std{:}),'uni',false)) ...
    .*mua_nan_trials;

mua_ctxpred_allcat = cell2mat(cellfun(@(mua,mua_baseline,mua_std) ...
    (mua-mua_baseline)./mua_std,vertcat(mua_ctxpred_all{:}), ...
    vertcat(mua_day_baseline{:}),vertcat(mua_day_std{:}),'uni',false)) ...
    .*mua_nan_trials;

% Don't subtract baseline for task-predicted: done before regression
mua_taskpred_allcat = cell2mat(cellfun(@(mua,mua_baseline,mua_std) ...
    (mua)./mua_std,vertcat(mua_taskpred_all{:}), ...
    vertcat(mua_day_baseline{:}),vertcat(mua_day_std{:}),'uni',false)) ...
    .*mua_nan_trials;
mua_taskpred_reduced_allcat = cell2mat(cellfun(@(mua,mua_baseline,mua_std) ...
    (mua)./mua_std,vertcat(mua_taskpred_reduced_all{:}), ...
    vertcat(mua_day_baseline{:}),vertcat(mua_day_std{:}),'uni',false)) ...
    .*mua_nan_trials;

mua_ctxpred_taskpred_allcat = cell2mat(cellfun(@(mua,mua_baseline,mua_std) ...
    (mua)./mua_std,vertcat(mua_ctxpred_taskpred_all{:}), ...
    vertcat(mua_day_baseline{:}),vertcat(mua_day_std{:}),'uni',false)) ...
    .*mua_nan_trials;
mua_ctxpred_taskpred_reduced_allcat = cell2mat(cellfun(@(mua,mua_baseline,mua_std) ...
    (mua)./mua_std,vertcat(mua_ctxpred_taskpred_reduced_all{:}), ...
    vertcat(mua_day_baseline{:}),vertcat(mua_day_std{:}),'uni',false)) ...
    .*mua_nan_trials;

% Deconvolve fluorescence, subtract baseline
fluor_allcat_deconv = AP_deconv_wf(fluor_allcat);
fluor_allcat_deconv_baseline = nanmean(reshape(fluor_allcat_deconv(:,t_baseline,:),[],1,n_vs));
fluor_allcat_deconv = fluor_allcat_deconv - fluor_allcat_deconv_baseline;

% Get fluorescence ROIs
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
% wf_roi = wf_roi(:,1); % (left ROIs)
n_rois = numel(wf_roi); % (left and right ROIs)

fluor_roi_deconv = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
    n_rois,[],size(fluor_allcat_deconv,1)),[3,2,1]);



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
% 







