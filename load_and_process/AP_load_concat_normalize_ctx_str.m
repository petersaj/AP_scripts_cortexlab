function [fluor_allcat,fluor_roi_diff,mua_allcat,wheel_allcat,reward_allcat,D_allcat] = ...
    AP_load_concat_normalize_ctx_str(data_fn,exclude_data)
% [fluor_allcat,fluor_roi_diff,mua_allcat,wheel_allcat,reward_allcat,D_allcat] = ...
%     AP_load_concat_normalize_ctx_str(data_fn,exclude_data)
%
% Loads, concatenates, and normalizes trial data for the ctx-str project
% data_fn - filename (e.g. trial_activity_choiceworld)
% exclude_data - true/false at the moment (for bad behavior days)
% (can in future extend this to wonky recording days?)

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time (should save this in the data)
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Load pre-marked experiments to exclude and cut out bad ones
if exist('exclude_data','var') && isempty(exclude_data)
    if exclude_data
        exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
        exclude_fn{1} = 'bhv_use_experiments';
        % exclude_fn{2} = 'expl_var_use_experiments';
        use_experiments_all = {};
        for curr_exclude = 1:length(exclude_fn)
            curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
            use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
        end
        use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
            1:size(use_experiments_all,2),'uni',false);
        
        D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
        fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
        mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
        wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);
        if exist('reward_all','var')
            reward_all = cellfun(@(data,use_expts) data(use_expts),reward_all,use_experiments','uni',false);
        end
    end
end

% If any animals have empty data, remove them
use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);
if exist('reward_all','var')
    reward_all = reward_all(use_animals);
end

% Get number of widefield components and MUA depths
n_vs = size(fluor_all{1}{1},3);
n_depths = size(mua_all{1}{1},3);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_vs,1);
mua_allcat = nan(0,length(t),n_depths,1);
wheel_allcat = nan(0,length(t),1);
reward_allcat = nan(0,length(t),1);
D_allcat = struct('stimulus',cell(1,1),'response',cell(1,1),'repeatNum',cell(1,1));

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
    t_baseline = t < 0;
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths),'uni',false));
    
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x,[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths),'uni',false));
    
    softnorm = 20;
    mua_cat_norm = (mua_cat_raw_smoothed-mua_day_baseline)./(mua_day_std+softnorm);
    
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
    
    % Concatenate behavioural data
    stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));    
    if exist('reward_all','var')
        response = cell2mat(cellfun(@(x) x.response,D_all{curr_animal},'uni',false));
        repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all{curr_animal},'uni',false));  
    end

    % Concatenate to across-day data
    fluor_allcat = [fluor_allcat;fluor_cat];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    D_allcat.stimulus = vertcat(D_allcat.stimulus,stimulus);
    if exist('reward_all','var')
        reward_allcat = [reward_allcat;vertcat(reward_all{curr_animal}{:})];
        D_allcat.response = vertcat(D_allcat.response,response);
        D_allcat.repeatNum = vertcat(D_allcat.repeatNum,repeatNum);
    end
    
end

% Get fluorescence in ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);
n_rois = numel(wf_roi);

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

roi_mask = cat(3,wf_roi.mask);

fluor_roi = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat,[3,2,1]),n_vs,[]),[],[],roi_mask), ...
    size(roi_mask,3),[],size(fluor_allcat,1)),[3,2,1]);

% (diff, normalize)
fluor_roi_diff_raw = diff(fluor_roi,[],2);
fluor_roi_std = nanstd(reshape(fluor_roi_diff_raw,[],1,size(fluor_roi_diff_raw,3)),[],1);
fluor_roi_diff = fluor_roi_diff_raw./fluor_roi_std;












