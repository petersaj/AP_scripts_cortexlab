% AP_load_concat_normalize_ctx_str(data_fn,exclude_data)
%
% Loads, concatenates, and normalizes trial data for the ctx-str project
% data_fn - base filename(s) (cell array if > 1, will concatenate data)
% exclude_data - true/false at the moment (for bad behavior days)

%% Load in and unpack data

% (turn on warnings)
warning on;

% Load data (saved as structure trial_data_all)
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';

if ischar(data_fn)
    % Single dataset
    disp(['Loading ' data_fn '...']);
    load([trial_data_path filesep data_fn]);
    
elseif iscell(data_fn)
    % Multiple datasets to merge (this is really dirty)
    
    preload_vars = who;
    
    % (load in split)
    clear trial_data_all_split
    for curr_load_data = 1:length(data_fn)
        disp(['Loading ' data_fn{curr_load_data} '...']);
        temp_data{curr_load_data} = load([trial_data_path filesep data_fn{curr_load_data}],'trial_data_all');
%         trial_data_all_split(curr_data) = temp_data.trial_data_all;
    end
    
    % (find shared fields - need to iterate intersect if > 2 datasets)
    split_fieldnames = cellfun(@(x) fieldnames(x.trial_data_all),temp_data,'uni',false);
    intersect_fieldnames = intersect(split_fieldnames{1},split_fieldnames{2});
    if length(data_fn) > 2
        for curr_intersect = 3:length(data_fn)
            intersect_fieldnames = intersect(intersect_fieldnames,split_fieldnames{curr_intersect});
        end
    end
    
    % (initialize combined structure)
    trial_data_all = cell2struct(cell(size(intersect_fieldnames)),intersect_fieldnames);
    data_animals = arrayfun(@(x) temp_data{x}.trial_data_all.animals,1:length(data_fn),'uni',false);
    trial_data_all.animals = horzcat(data_animals{:});
    
    % (concatenate experiment fields)
    experiment_fields = cellfun(@(curr_field) ...
        length(temp_data{1}.trial_data_all.(curr_field)) == ...
        length(temp_data{1}.trial_data_all.animals) && ...
        iscell(temp_data{1}.trial_data_all.(curr_field)) && ...
        any(cellfun(@(x) iscell(x),temp_data{1}.trial_data_all.(curr_field))),intersect_fieldnames);

    for curr_field = intersect_fieldnames(experiment_fields)'
       for curr_load_data = 1:length(data_fn)
          trial_data_all.(curr_field{:}) = ...
              cat(1,trial_data_all.(curr_field{:}), ...
              reshape(temp_data{curr_load_data}.trial_data_all.(curr_field{:}),[],1));
       end
    end    
    
    % (grab non-experiment fields from the first dataset)
    % (NOTE: this assumes they're the same)
    for curr_field = intersect_fieldnames(~experiment_fields & ...
            ~strcmp(intersect_fieldnames,'animals'))'
        trial_data_all.(curr_field{:}) = ...
            temp_data{1}.trial_data_all.(curr_field{:});
    end
        
    % (if mixing protocols: ensure stimIDs convered into contrast*side)
    for curr_animal = 1:length(trial_data_all.trial_info_all)
        for curr_day = 1:length(trial_data_all.trial_info_all{curr_animal})
            
            curr_stim = trial_data_all.trial_info_all{curr_animal}{curr_day}.stimulus;
            if length(unique(curr_stim)) == 3 && all(unique(curr_stim) == [1;2;3])
                % Stim [1,2,3] are coded stim IDs
                curr_stim_recode = curr_stim;
                % For both mpep and signals, ID is 1 = left, 2 = center, 3 = right
                curr_stim_recode(curr_stim == 1) = -1;
                curr_stim_recode(curr_stim == 2) = 0;
                curr_stim_recode(curr_stim == 3) = 1;
                
                trial_data_all.trial_info_all{curr_animal}{curr_day}.stimulus = ...
                    curr_stim_recode;
            end
            
        end
    end
    
    % clear temp variables
    clearvars('-except',preload_vars{:},'data_fn','trial_data_all')
    
else
    error('Unrecognized data_fn type')
end

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

%% Get task/stim-relevant

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
    
    if length(unique(trial_info_allcat.stimulus)) == 3 && ...
            all(unique(trial_info_allcat.stimulus) == [1;2;3])
        % Stim [1,2,3] are coded stim IDs
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

%% Cortical fluorescence

% Get number of widefield ROIs 
n_vs = size(fluor_all{end}{end},3);

% Concatenate cortex data
fluor_allcat = cell2mat(vertcat(fluor_all{:}));
if task_dataset
    fluor_taskpred_allcat = cell2mat(vertcat(fluor_taskpred_all{:}));
    fluor_taskpred_reduced_allcat = cell2mat(vertcat(fluor_taskpred_reduced_all{:}));
end

% Deconvolve fluorescence, subtract average baseline
fluor_allcat_deconv = AP_deconv_wf(fluor_allcat);
fluor_allcat_deconv_baseline = nanmean(reshape(fluor_allcat_deconv(:,t_baseline,:),[],1,n_vs));
fluor_allcat_deconv = fluor_allcat_deconv - fluor_allcat_deconv_baseline;

% Get fluorescence ROIs
% (by cortical area)
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);
fluor_roi_deconv = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
    n_rois,[],size(fluor_allcat_deconv,1)),[3,2,1]);

% (from striatum kernels)
kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
load(kernel_roi_fn);
n_kernel_rois = size(kernel_roi.bw,3);
fluor_kernelroi_deconv = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],kernel_roi.bw), ...
    size(kernel_roi.bw,3),[],size(fluor_allcat_deconv,1)),[3,2,1]);

% Normalize ROI fluorescence by baseline and std
t_baseline = t < 0;
n_trials_day = cellfun(@(x) size(x,1),vertcat(fluor_all{:}));

% (baseline of every trial)
fluor_roi_trial_baseline = cellfun(@(fluor) ...
    nanmean(fluor(:,t_baseline,:),2), ...
    mat2cell(fluor_roi_deconv,n_trials_day,length(t),n_rois),'uni',false);
fluor_kernelroi_trial_baseline = cellfun(@(fluor) ...
    nanmean(fluor(:,t_baseline,:),2), ...
    mat2cell(fluor_kernelroi_deconv,n_trials_day,length(t),n_kernel_rois),'uni',false);

% (std of mean-subtracted baseline)
fluor_roi_day_baseline_std = ...
    cellfun(@(x,baseline) repmat(nanstd(reshape(x(:,t_baseline,:) - baseline,[],1,n_rois),[],1),size(x,1),1), ...
    mat2cell(fluor_roi_deconv,n_trials_day,length(t),n_rois),fluor_roi_trial_baseline,'uni',false);
fluor_kernelroi_day_baseline_std = ...
    cellfun(@(x,baseline) repmat(nanstd(reshape(x(:,t_baseline,:) - baseline,[],1,n_kernel_rois),[],1),size(x,1),1), ...
    mat2cell(fluor_kernelroi_deconv,n_trials_day,length(t),n_kernel_rois),fluor_kernelroi_trial_baseline,'uni',false);

% (std of whole day)
fluor_roi_day_std = ...
    cellfun(@(x,baseline) repmat(nanstd(reshape(x(:,:,:) - baseline,[],1,n_rois),[],1),size(x,1),1), ...
    mat2cell(fluor_roi_deconv,n_trials_day,length(t),n_rois),fluor_roi_trial_baseline,'uni',false);
fluor_kernelroi_day_std = ...
    cellfun(@(x,baseline) repmat(nanstd(reshape(x(:,:,:) - baseline,[],1,n_kernel_rois),[],1),size(x,1),1), ...
    mat2cell(fluor_kernelroi_deconv,n_trials_day,length(t),n_kernel_rois),fluor_kernelroi_trial_baseline,'uni',false);

% Normalize fluor

% % (don't normalize)
% fluor_roi_norm = 1;
% fluor_kernelroi_norm = 1;

% % % (use day std)
% fluor_roi_norm = fluor_roi_day_std;
% fluor_kernelroi_norm = fluor_kernelroi_day_std;

% (use day baseline std)
fluor_roi_norm = fluor_roi_day_baseline_std;
fluor_kernelroi_norm = fluor_kernelroi_day_baseline_std;

fluor_roi_deconv = (fluor_roi_deconv - cell2mat(fluor_roi_trial_baseline))./cell2mat(fluor_roi_norm);
fluor_kernelroi_deconv = (fluor_kernelroi_deconv - cell2mat(fluor_kernelroi_trial_baseline))./cell2mat(fluor_kernelroi_norm);

if task_dataset
    % Get task-predicted activity
    fluor_roi_taskpred = permute(reshape( ...
        AP_svd_roi(U_master(:,:,1:n_vs), ...
        reshape(permute(fluor_taskpred_allcat,[3,2,1]),n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
        n_rois,[],size(fluor_taskpred_allcat,1)),[3,2,1])./cell2mat(fluor_roi_norm);
    fluor_roi_taskpred_reduced = cell2mat(permute(arrayfun(@(x) ...
        permute(reshape( ...
        AP_svd_roi(U_master(:,:,1:n_vs), ...
        reshape(permute(fluor_taskpred_reduced_allcat(:,:,:,x),[3,2,1]), ...
        n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
        n_rois,[],size(fluor_taskpred_reduced_allcat,1)),[3,2,1]), ...
        1:size(fluor_taskpred_reduced_allcat,4),'uni',false),[1,3,4,2]))./cell2mat(fluor_roi_norm);
    
    fluor_kernelroi_taskpred = permute(reshape( ...
        AP_svd_roi(U_master(:,:,1:n_vs), ...
        reshape(permute(fluor_taskpred_allcat,[3,2,1]),n_vs,[]),[],[],kernel_roi.bw), ...
        size(kernel_roi.bw,3),[],size(fluor_taskpred_allcat,1)),[3,2,1])./cell2mat(fluor_kernelroi_norm);   
    fluor_kernelroi_taskpred_reduced = cell2mat(permute(arrayfun(@(x) ...
        permute(reshape( ...
        AP_svd_roi(U_master(:,:,1:n_vs), ...
        reshape(permute(fluor_taskpred_reduced_allcat(:,:,:,x),[3,2,1]), ...
        n_vs,[]),[],[],kernel_roi.bw), ...
        size(kernel_roi.bw,3),[],size(fluor_taskpred_reduced_allcat,1)),[3,2,1]), ...
        1:size(fluor_taskpred_reduced_allcat,4),'uni',false),[1,3,4,2]))./cell2mat(fluor_kernelroi_norm);
    
    % Make move-aligned fluorescence
    fluor_allcat_deconv_move = fluor_allcat_deconv;
    fluor_taskpred_reduced_allcat_move = fluor_taskpred_reduced_allcat;
    fluor_roi_deconv_move = fluor_roi_deconv;
    fluor_roi_taskpred_move = fluor_roi_taskpred;
    fluor_roi_taskpred_reduced_move = fluor_roi_taskpred_reduced;
    fluor_kernelroi_deconv_move = fluor_kernelroi_deconv;
    
    t_leeway = -t(1);
    leeway_samples = round(t_leeway*(sample_rate));
    for i = 1:size(fluor_allcat_deconv,1)
        fluor_allcat_deconv_move(i,:,:,:) = circshift(fluor_allcat_deconv_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
        fluor_taskpred_reduced_allcat_move(i,:,:,:) = circshift(fluor_taskpred_reduced_allcat_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
        fluor_roi_deconv_move(i,:,:,:) = circshift(fluor_roi_deconv_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
        fluor_roi_taskpred_move(i,:,:,:) = circshift(fluor_roi_taskpred_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
        fluor_roi_taskpred_reduced_move(i,:,:,:) = circshift(fluor_roi_taskpred_reduced_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
        fluor_kernelroi_deconv_move(i,:,:,:) = circshift(fluor_kernelroi_deconv_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
    end
end

%% Striatal multiunit

if exist('mua_all','var')
    
    % Get number of striatum depths
    n_depths = size(mua_all{end}{end},3);
    
    % Normalize MUA
    
    % (day baseline)
    mua_day_baseline = cellfun(@(mua_animal) cellfun(@(mua_day) ...
        nanmean(reshape(mua_day(:,t_baseline,:),1,[],size(mua_day,3)),2), ...
        mua_animal,'uni',false),mua_all,'uni',false);   
    
    % (baseline of every trial)
    mua_trial_baseline = cellfun(@(mua_animal) cellfun(@(mua_day) ...
        nanmean(mua_day(:,t_baseline,:),2), ...
        mua_animal,'uni',false),mua_all,'uni',false);    
    
    % (std of trial-baseline subtracted activity within day)
    mua_day_std = cellfun(@(mua_animal,baseline) cellfun(@(mua_day,baseline) ...
        nanstd(reshape(mua_day - baseline,[],1,n_depths)), ...
        mua_animal,baseline,'uni',false),mua_all,mua_trial_baseline,'uni',false);
    
    % (std of trial-baseline subtracted baseline within day)
    mua_day_baseline_std = cellfun(@(mua_animal,baseline) cellfun(@(mua_day,baseline) ...
        nanstd(reshape(mua_day(:,t_baseline,:) - baseline,[],1,n_depths)), ...
        mua_animal,baseline,'uni',false),mua_all,mua_trial_baseline,'uni',false);
    
    % (set baseline and norm)
    mua_use_baseline = mua_trial_baseline;
    mua_norm = mua_day_baseline_std; % mua_day_baseline_std or mua_day_std
    
    % (NaN-out days with no spikes)
    mua_nan_trials = cell2mat(cellfun(@(x) ...
        +repmat(any(reshape(x,[],1,n_depths),1),size(x,1),1),vertcat(mua_all{:}),'uni',false));
    mua_nan_trials(~mua_nan_trials) = NaN;
    
    mua_allcat = cell2mat(cellfun(@(mua,mua_baseline,mua_norm) ...
        (mua-mua_baseline)./mua_norm,vertcat(mua_all{:}), ...
        vertcat(mua_use_baseline{:}),vertcat(mua_norm{:}),'uni',false)) ...
        .*mua_nan_trials;
    
    mua_ctxpred_allcat = cell2mat(cellfun(@(mua,mua_baseline,mua_norm) ...
        (mua-mua_baseline)./mua_norm,vertcat(mua_ctxpred_all{:}), ...
        vertcat(mua_use_baseline{:}),vertcat(mua_norm{:}),'uni',false)) ...
        .*mua_nan_trials;
    
    if task_dataset
        % Get task-predicted activity
        % (don't subtract baseline for task-predicted: done before regression)
        mua_taskpred_allcat = cell2mat(cellfun(@(mua,mua_baseline,mua_norm) ...
            (mua)./mua_norm,vertcat(mua_taskpred_all{:}), ...
            vertcat(mua_use_baseline{:}),vertcat(mua_norm{:}),'uni',false)) ...
            .*mua_nan_trials;
        mua_taskpred_reduced_allcat = cell2mat(cellfun(@(mua,mua_baseline,mua_norm) ...
            (mua)./mua_norm,vertcat(mua_taskpred_reduced_all{:}), ...
            vertcat(mua_use_baseline{:}),vertcat(mua_norm{:}),'uni',false)) ...
            .*mua_nan_trials;
        
        mua_ctxpred_taskpred_allcat = cell2mat(cellfun(@(mua,mua_baseline,mua_norm) ...
            (mua)./mua_norm,vertcat(mua_ctxpred_taskpred_all{:}), ...
            vertcat(mua_use_baseline{:}),vertcat(mua_norm{:}),'uni',false)) ...
            .*mua_nan_trials;
        mua_ctxpred_taskpred_reduced_allcat = cell2mat(cellfun(@(mua,mua_baseline,mua_norm) ...
            (mua)./mua_norm,vertcat(mua_ctxpred_taskpred_reduced_all{:}), ...
            vertcat(mua_use_baseline{:}),vertcat(mua_norm{:}),'uni',false)) ...
            .*mua_nan_trials;
        
        % Make move-aligned MUA
        mua_allcat_move = mua_allcat;
        mua_ctxpred_allcat_move = mua_ctxpred_allcat;
        mua_taskpred_allcat_move = mua_taskpred_allcat;
        mua_taskpred_reduced_allcat_move = mua_taskpred_reduced_allcat;
        mua_ctxpred_taskpred_allcat_move = mua_taskpred_allcat;
        mua_ctxpred_taskpred_reduced_allcat_move = mua_taskpred_reduced_allcat;
        
        t_leeway = -t(1);
        leeway_samples = round(t_leeway*(sample_rate));
        for i = 1:size(mua_allcat,1)
            mua_allcat_move(i,:,:) = circshift(mua_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
            mua_ctxpred_allcat_move(i,:,:) = circshift(mua_ctxpred_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
            mua_taskpred_allcat_move(i,:,:) = circshift(mua_taskpred_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
            mua_taskpred_reduced_allcat_move(i,:,:,:) = circshift(mua_taskpred_reduced_allcat_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
            mua_ctxpred_taskpred_allcat_move(i,:,:) = circshift(mua_ctxpred_taskpred_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
            mua_ctxpred_taskpred_reduced_allcat_move(i,:,:,:) = circshift(mua_ctxpred_taskpred_reduced_allcat_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
        end      
    end
    
end


%% Finish

disp('Finished loading trials')







