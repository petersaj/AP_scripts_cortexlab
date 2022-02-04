% AP_load_trials_wf(data_fn,exclude_data)
%
% Loads widefield trial data
% from 'operant_grab_trial_data'
% ASSUMES DECONVOLUTION ALREADY DONE
% data_fn - base filename(s) (cell array if > 1, will concatenate data)
% exclude_data - true/false at the moment (for bad behavior days)

%% Load in and unpack data

% (turn on warnings)
warning on;

% Load data (saved as structure trial_data_all)
if ~exist('trial_data_path','var')
    % If no path specified, assume paper data
    trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
end

if ischar(data_fn)
    % Single dataset
    disp(['Loading ' data_fn '...']);
    load([trial_data_path filesep data_fn]);
    
elseif iscell(data_fn)
    % Multiple datasets to merge (this is really dirty)
    
    preload_vars = who;
    
    % (load in all datasets)
    clear temp_data
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

% Concatenate cortex data, subtract baseline
fluor_allcat = cell2mat(vertcat(fluor_all{:}));
if task_dataset
    % (concatenate)
    fluor_taskpred_allcat = cell2mat(vertcat(fluor_taskpred_all{:}));
    fluor_taskpred_reduced_allcat = cell2mat(vertcat(fluor_taskpred_reduced_all{:}));
    
%     % (subtract baseline on each trial)
%     fluor_taskpred_allcat = fluor_taskpred_allcat - ...
%         nanmean(fluor_taskpred_allcat(:,t_baseline,:),2);
%     fluor_taskpred_reduced_allcat = fluor_taskpred_reduced_allcat - ...
%         nanmean(fluor_taskpred_reduced_allcat(:,t_baseline,:,:),2);
end

% % Deconvolve fluorescence, subtract baseline for each trial
%%%% CHANGED: DECONV DONE IN LOADING SCRIPT (to allow resampling)
% fluor_allcat_deconv = AP_deconv_wf(fluor_allcat) - ...
%     nanmean(AP_deconv_wf(fluor_allcat(:,t_baseline,:)),2);
fluor_allcat_deconv = fluor_allcat;
clear fluor_allcat;

% Get fluorescence ROIs
% (by cortical area)
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);
fluor_roi_deconv = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
    n_rois,[],size(fluor_allcat_deconv,1)),[3,2,1]);

if task_dataset
    % Get task-predicted activity
    fluor_roi_taskpred = permute(reshape( ...
        AP_svd_roi(U_master(:,:,1:n_vs), ...
        reshape(permute(fluor_taskpred_allcat,[3,2,1]),n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
        n_rois,[],size(fluor_taskpred_allcat,1)),[3,2,1]);
    fluor_roi_taskpred_reduced = cell2mat(permute(arrayfun(@(x) ...
        permute(reshape( ...
        AP_svd_roi(U_master(:,:,1:n_vs), ...
        reshape(permute(fluor_taskpred_reduced_allcat(:,:,:,x),[3,2,1]), ...
        n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
        n_rois,[],size(fluor_taskpred_reduced_allcat,1)),[3,2,1]), ...
        1:size(fluor_taskpred_reduced_allcat,4),'uni',false),[1,3,4,2]));
    
%     % Make move-aligned fluorescence
%     fluor_allcat_deconv_move = fluor_allcat_deconv;
%     fluor_taskpred_reduced_allcat_move = fluor_taskpred_reduced_allcat;
%     fluor_roi_deconv_move = fluor_roi_deconv;
%     fluor_roi_taskpred_move = fluor_roi_taskpred;
%     fluor_roi_taskpred_reduced_move = fluor_roi_taskpred_reduced;
%     
%     t_leeway = -t(1);
%     leeway_samples = round(t_leeway*(sample_rate));
%     for i = 1:size(fluor_allcat_deconv,1)
%         fluor_allcat_deconv_move(i,:,:,:) = circshift(fluor_allcat_deconv_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
%         fluor_taskpred_reduced_allcat_move(i,:,:,:) = circshift(fluor_taskpred_reduced_allcat_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
%         fluor_roi_deconv_move(i,:,:,:) = circshift(fluor_roi_deconv_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
%         fluor_roi_taskpred_move(i,:,:,:) = circshift(fluor_roi_taskpred_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
%         fluor_roi_taskpred_reduced_move(i,:,:,:) = circshift(fluor_roi_taskpred_reduced_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
%     end

end

%% Finish

disp('Finished loading trials')







