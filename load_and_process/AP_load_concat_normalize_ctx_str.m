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

%%%%%%%%%% TO ADD: 
% mua_ctxpred_all mua_taskpred_all mua_taskpred_reduced_all


% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
load([data_path filesep data_fn]);
n_animals = length(D_all);

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
        predicted_mua_std_all = cellfun(@(data,use_expts) data(use_expts),predicted_mua_std_all,use_experiments','uni',false);
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
if exist('mua_ctxpred_all','var')
    mua_ctxpred_all = mua_ctxpred_all(use_animals);
end
if exist('mua_taskpred_all','var')
    mua_taskpred_all = mua_taskpred_all(use_animals);
    mua_taskpred_reduced_all = mua_taskpred_reduced_all(use_animals);
end

% Get number of widefield components and MUA depths
n_vs = size(fluor_all{1}{1},3);
n_depths = size(mua_all{1}{1},3);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_vs,1);
mua_allcat = nan(0,length(t),n_depths,1);
mua_ctxpred_allcat = nan(0,length(t),n_depths,1);
mua_taskpred_allcat = nan(0,length(t),n_depths,1);
mua_taskpred_reduced_allcat = nan(0,length(t),n_depths,4);
wheel_allcat = nan(0,length(t),1);
reward_allcat = nan(0,length(t),1);
D_allcat = struct('stimulus',cell(1,1),'response',cell(1,1),'repeatNum',cell(1,1));

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
   
    % Normalize and concatenate MUA
    
    % (filter to align with widefield and predictions)   
    sample_rate = 1/mean(diff(t));
    lowpassCutoff = 6; % Hz
    [b100s, a100s] = butter(2, lowpassCutoff/((sample_rate)/2), 'low');
    filter_pad = 10;
    curr_mua_pad = cellfun(@(x) filter(b100s,a100s, ...
        padarray(x,[0,filter_pad,0],'both','replicate'),[],2), ...
        mua_all{curr_animal},'uni',false);
    curr_mua = cellfun(@(x) x(:,filter_pad+1:end-filter_pad,:),curr_mua_pad,'uni',false);
    
    % (NaN-out trials with no MUA)
    for curr_day = 1:length(curr_mua)
        curr_unfilled_trials = +any(mua_all{curr_animal}{curr_day},2);
        curr_unfilled_trials(curr_unfilled_trials == 0) = NaN;
        curr_mua{curr_day} = curr_mua{curr_day}.*curr_unfilled_trials;
    end
    
    % (subtract baseline, divide by std across recordings)
    t_baseline = t < 0;
    curr_mua_day_baseline = cellfun(@(x) nanmean(reshape(x(:,t_baseline,:),[],1,n_depths),1),curr_mua,'uni',false); 
    curr_mua_day_std = cellfun(@(x) nanstd(reshape(x,[],1,n_depths),[],1),curr_mua,'uni',false);
    curr_mua_norm = cellfun(@(raw,baseline,std) (raw-baseline)./std, ...
        curr_mua,curr_mua_day_baseline,curr_mua_day_std,'uni',false);
    mua_cat_norm = cat(1,curr_mua_norm{:});
    
    % Concatenated cortex-predicted MUA (normalize by measured MUA)
    if exist('mua_ctxpred_all','var')
        curr_mua_ctxpred = mua_ctxpred_all{curr_animal};
        curr_mua_ctxpred_norm = cellfun(@(raw,baseline,std) (raw-baseline)./std, ...
            curr_mua_ctxpred,curr_mua_day_baseline,curr_mua_day_std,'uni',false);
        mua_ctxpred_cat = cat(1,curr_mua_ctxpred_norm{:});
    end
    
    % Concatenate task-predicted MUA (normalize by measured MUA)
    if exist('mua_taskpred_all','var')
        curr_mua_taskpred = mua_taskpred_all{curr_animal};
        curr_mua_taskpred_norm = cellfun(@(raw,baseline,std) (raw-baseline)./std, ...
            curr_mua_taskpred,curr_mua_day_baseline,curr_mua_day_std,'uni',false);
       mua_taskpred_cat = cat(1,curr_mua_taskpred_norm{:});
        
        curr_mua_taskpred_reduced = mua_taskpred_reduced_all{curr_animal};
        curr_mua_taskpred_reduced_norm = cellfun(@(raw,baseline,std) (raw-baseline)./std, ...
            curr_mua_taskpred_reduced,curr_mua_day_baseline,curr_mua_day_std,'uni',false);
        mua_taskpred_reduced_cat = cat(1,curr_mua_taskpred_reduced_norm{:});
    end
    
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
    if exist('mua_ctxpred_all','var')
        mua_ctxpred_allcat = [mua_ctxpred_allcat;mua_ctxpred_cat];
    end
    if exist('mua_taskpred_all','var')
        mua_taskpred_allcat = [mua_taskpred_allcat;mua_taskpred_cat];
        mua_taskpred_reduced_allcat = [mua_taskpred_reduced_allcat;mua_taskpred_reduced_cat];
    end
    
end

% Deconvolve fluorescence
fluor_allcat_deriv = AP_deconv_wf(fluor_allcat);

% Get fluorescence ROIs
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);
n_rois = numel(wf_roi);

fluor_roi_deriv = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deriv,[3,2,1]),n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
    n_rois,[],size(fluor_allcat_deriv,1)),[3,2,1]);

% I have no idea what's right here, can't do < 0 because don't know offset?
% t_baseline = t < -0.1;
% fluor_roi_deriv_baseline = nanmean(nanmean(fluor_roi_deriv(:,t_baseline,:),2),1);
% for curr_roi = 1:n_rois
%     fluor_roi_deriv(fluor_roi_deriv(:,:,curr_roi) < ...
%         fluor_roi_deriv_baseline(curr_roi)) = fluor_roi_deriv_baseline(curr_roi);
% end
% %%%%%%%%%%

% (storing filter here - implement later in MUA?
% % Filter and std-normalize spikes
% lowpassCutoff = 6; % Hz
% [b100s, a100s] = butter(2, lowpassCutoff/((sample_rate)/2), 'low');
% binned_spikes_filt = filter(b100s,a100s,binned_spikes,[],2);
% binned_spikes_filt_std = binned_spikes_filt./nanstd(binned_spikes_filt,[],2);
% binned_spikes_filt_std(isnan(binned_spikes_filt_std)) = 0;



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


