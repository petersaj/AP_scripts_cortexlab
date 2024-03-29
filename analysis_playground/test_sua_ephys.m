%% ~~~~~~~~~~~~ Load development

%% --> INTEGRATE THIS: bugfixes from Julie
% TANS now just in one group
% High prop ISI UINs now not included (as others) and are instead short
% waveform long prop isi cells that behave like TANs but with short
% waveforms

% SUAbugsBgone
% with some good old for loops
%clear all;
%load('E:\analysis\wf_ephys_choiceworld\paper\data\trial_activity_naive.mat')
removeMissingSpikesUnits = true;
for iAnimal = 1:size(trial_data_all.goodUnits, 1)
    for iRecording = 1:size(trial_data_all.goodUnits{iAnimal}, 1)

        % 1. if flag remove spikes missing, add these to the goodUnits
        if removeMissingSpikesUnits
            trial_data_all.goodUnits{iAnimal}{iRecording} = trial_data_all.goodUnits{iAnimal}{iRecording} ...
                & nanmin(trial_data_all.percent_missing_ndtr{iAnimal}{iRecording}(2:end, :)) < 30;
        end

        % 2. re-classify to: add noise guys to uins, remove large post spike suppression FSI and UINs, and short waveform guys from TANs. 
        if ~isempty(trial_data_all.acg{iAnimal}{iRecording})
         for iCell = 1:size(trial_data_all.acg{iAnimal}{iRecording}, 1)
            pss_allcat2temp = find(trial_data_all.acg{iAnimal}{iRecording}(iCell, 500:1000) >= nanmean(trial_data_all.acg{iAnimal}{iRecording}(iCell, 600:900)));
            pss_allcat2(iCell) = pss_allcat2temp(1);
         end
        allDepths = zeros(size(trial_data_all.allGroups{iAnimal}{iRecording},2),1);
        allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [1,4,7,10,13,16]))=1; 
        allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [2,5,8,11,14,17]))=2;
        allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [3,6,9,12,15,18]))=3;
        largePssShorties = find(trial_data_all.templateDuration{iAnimal}{iRecording} < 400 & pss_allcat2' > 40); 
        
        fsi = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 < 40 & trial_data_all.prop_long_isi{iAnimal}{iRecording} <= 0.1;
        uin = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 < 40 & trial_data_all.prop_long_isi{iAnimal}{iRecording} > 0.1;
        tan = trial_data_all.templateDuration{iAnimal}{iRecording} > 400 & pss_allcat2 >= 40;
        msn = trial_data_all.templateDuration{iAnimal}{iRecording} > 400 & pss_allcat2 < 40;
        shortDurLongPss = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 >= 40;
        
        trial_data_all.allGroups{iAnimal}{iRecording}(msn) = allDepths(msn);
        trial_data_all.allGroups{iAnimal}{iRecording}(fsi) = allDepths(fsi) + 3;
        trial_data_all.allGroups{iAnimal}{iRecording}(tan) = allDepths(tan) + 6;
        trial_data_all.allGroups{iAnimal}{iRecording}(uin) = allDepths(uin) + 9;
        trial_data_all.allGroups{iAnimal}{iRecording}(shortDurLongPss) = allDepths(shortDurLongPss) + 12;
        
        clearvars pss_allcat2
        end
    end
end






%% Load SUA data (single dataset)

disp('Loading SUA data');

% data_fn = 'G:\JF_single_cell_data\trial_activity_choiceworld.mat';
% data_fn = 'G:\JF_single_cell_data\trial_activity_ctx_task.mat';
% data_fn = 'G:\JF_single_cell_data\trial_activity_muscimol_task.mat';

data_fn = 'G:\JF_single_cell_data\trial_activity_naive';

load(data_fn);
exclude_data = false;

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

% BUGFIX? If any days don't have data - throw away
% (this isn't present in MUA dataset, don't know where this is from)
for curr_animal = 1:length(trial_data_all.animals)
    nodata_days = cellfun(@(x) isempty(x),trial_data_all.trial_info_all{curr_animal});
    if any(nodata_days)
        for curr_field = data_struct_fieldnames(experiment_fields)'
            trial_data_all.(cell2mat(curr_field)){curr_animal}(nodata_days) = [];
        end
    end
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
    
    if length(unique(trial_info_allcat.stimulus)) == 3 && ...
            all(unique(trial_info_allcat.stimulus) == [1;2;3])
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

% Concatenate good units
good_units_exp = cellfun(@transpose,vertcat(goodUnits{:}),'uni',false);
good_units_allcat = cell2mat(vertcat(goodUnits{:})')';

% (groups are cell type * depth: msn, fsi, tan, th, ??)
n_aligned_depths = 3; % hardcoded: I think not stored
n_celltypes = max(groups_allcat)./n_aligned_depths;
if n_celltypes ~= 6
    error('Incorrect number of celltypes')
end
celltype_allcat = ceil(groups_allcat./n_aligned_depths);
% NOTE: Group 3+6 are both TAN (from relaxing metrics)
% COMBINING THESE MANUALLY
celltype_labels = {'MSN','FSI','TAN','UIN','Other'};
celltype_allcat(celltype_allcat == 6) = 3;
n_celltypes = 5;

domain_allcat = mod(groups_allcat,n_aligned_depths) + ...
    n_aligned_depths*(mod(groups_allcat,n_aligned_depths) == 0);

% Get maps for all cells
% (load master U)
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');
% (use t = 0 kernels for all cells, hardcoded at the moment)
use_k_frame = 5;
ctx_str_k_frame = cell2mat(cellfun(@(x) reshape(x(:,use_k_frame,:),100,[]), ...
    vertcat(ctx_str_k_all{:})','uni',false));
ctx_str_k_px = svdFrameReconstruct(U_master(:,:,1:100),ctx_str_k_frame(:,good_units_allcat));


% Deconvolve fluorescence
fluor_deconv_allcat_exp = cellfun(@(x) AP_deconv_wf(x),vertcat(fluor_all{:}),'uni',false);


% Get animals/days/cell number
% (note this includes cortical cells, which are NaNs in allGroups)
animals_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp)), ...
    vertcat(allAnimals{:}),vertcat(allGroups{:}),'uni',false));
days_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp)), ...
    vertcat(allDays{:}),vertcat(allGroups{:}),'uni',false));
neurons_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp))', ...
    vertcat(allNeurons{:}),vertcat(allGroups{:}),'uni',false));
recordings_allcat = cell2mat(cellfun(@(x,grp) x*ones(sum(~isnan(grp)),1), ...    
    num2cell((1:length(vertcat(allDays{:})))'),vertcat(allGroups{:}),'uni',false));


disp('Finished.')


%% Load SUA data (combine datasets)


data_fn = { ...
    'trial_activity_choiceworld.mat', ... % original 
    'trial_activity_ctx_task.mat', ...    % + cortex ephys
    'trial_activity_muscimol_task.mat'};  % muscimol group

% (turn on warnings)
warning on;

% Load data (saved as structure trial_data_all)
trial_data_path = 'G:\JF_single_cell_data\';

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
    
    if length(unique(trial_info_allcat.stimulus)) == 3 && ...
            all(unique(trial_info_allcat.stimulus) == [1;2;3])
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

% Concatenate good units
good_units_exp = cellfun(@transpose,vertcat(goodUnits{:}),'uni',false);
good_units_allcat = cell2mat(vertcat(goodUnits{:})')';

% (groups are cell type * depth: msn, fsi, tan, th, ??)
n_aligned_depths = 3; % hardcoded: I think not stored
n_celltypes = max(groups_allcat)./n_aligned_depths;
if n_celltypes ~= 6
    error('Incorrect number of celltypes')
end
celltype_allcat = ceil(groups_allcat./n_aligned_depths);
% NOTE: Group 3+6 are both TAN (from relaxing metrics)
% COMBINING THESE MANUALLY
celltype_labels = {'MSN','FSI','TAN','UIN','Other'};
celltype_allcat(celltype_allcat == 6) = 3;
n_celltypes = 5;

domain_allcat = mod(groups_allcat,n_aligned_depths) + ...
    n_aligned_depths*(mod(groups_allcat,n_aligned_depths) == 0);

% Get maps for all cells
% (load master U)
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');
% (use t = 0 kernels for all cells, hardcoded at the moment)
use_k_frame = 5;
ctx_str_k_frame = cell2mat(cellfun(@(x) reshape(x(:,use_k_frame,:),100,[]), ...
    vertcat(ctx_str_k_all{:})','uni',false));
ctx_str_k_px = svdFrameReconstruct(U_master(:,:,1:100),ctx_str_k_frame(:,good_units_allcat));


% Deconvolve fluorescence
fluor_deconv_allcat_exp = cellfun(@(x) AP_deconv_wf(x),vertcat(fluor_all{:}),'uni',false);


% Get animals/days/cell number
% (note this includes cortical cells, which are NaNs in allGroups)
animals_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp)), ...
    vertcat(allAnimals{:}),vertcat(allGroups{:}),'uni',false));
days_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp)), ...
    vertcat(allDays{:}),vertcat(allGroups{:}),'uni',false));
neurons_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp))', ...
    vertcat(allNeurons{:}),vertcat(allGroups{:}),'uni',false));
recordings_allcat = cell2mat(cellfun(@(x,grp) x*ones(sum(~isnan(grp)),1), ...    
    num2cell((1:length(vertcat(allDays{:})))'),vertcat(allGroups{:}),'uni',false));


disp('Finished.')



%% Load SUA data (combine datasets, w/ classification bugfix)

data_fn = { ...
    'trial_activity_choiceworld.mat', ... % original 
    'trial_activity_ctx_task.mat', ...    % + cortex ephys
    'trial_activity_muscimol_task.mat'};  % muscimol group

% (turn on warnings)
warning on;

% Load data (saved as structure trial_data_all)
trial_data_path = 'G:\JF_single_cell_data\';

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


% Julie's bugfixes
% (combines 2 TAN categories, integrates old 'other' with UINs, puts short
% waveforms that look like TANs into new 'other')
removeMissingSpikesUnits = true;
for iAnimal = 1:size(trial_data_all.goodUnits, 1)
    for iRecording = 1:size(trial_data_all.goodUnits{iAnimal}, 1)

        % 1. if flag remove spikes missing, add these to the goodUnits
        if removeMissingSpikesUnits
            trial_data_all.goodUnits{iAnimal}{iRecording} = trial_data_all.goodUnits{iAnimal}{iRecording} ...
                & nanmin(trial_data_all.percent_missing_ndtr{iAnimal}{iRecording}(2:end, :)) < 30;
        end

        % 2. re-classify to: add noise guys to uins, remove large post spike suppression FSI and UINs, and short waveform guys from TANs. 
        if ~isempty(trial_data_all.acg{iAnimal}{iRecording})
         for iCell = 1:size(trial_data_all.acg{iAnimal}{iRecording}, 1)
            pss_allcat2temp = find(trial_data_all.acg{iAnimal}{iRecording}(iCell, 500:1000) >= nanmean(trial_data_all.acg{iAnimal}{iRecording}(iCell, 600:900)));
            pss_allcat2(iCell) = pss_allcat2temp(1);
         end
        allDepths = zeros(size(trial_data_all.allGroups{iAnimal}{iRecording},2),1);
        allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [1,4,7,10,13,16]))=1; 
        allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [2,5,8,11,14,17]))=2;
        allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [3,6,9,12,15,18]))=3;
        largePssShorties = find(trial_data_all.templateDuration{iAnimal}{iRecording} < 400 & pss_allcat2' > 40); 
        
        fsi = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 < 40 & trial_data_all.prop_long_isi{iAnimal}{iRecording} <= 0.1;
        uin = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 < 40 & trial_data_all.prop_long_isi{iAnimal}{iRecording} > 0.1;
        tan = trial_data_all.templateDuration{iAnimal}{iRecording} > 400 & pss_allcat2 >= 40;
        msn = trial_data_all.templateDuration{iAnimal}{iRecording} > 400 & pss_allcat2 < 40;
        shortDurLongPss = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 >= 40;
        
        trial_data_all.allGroups{iAnimal}{iRecording}(msn) = allDepths(msn);
        trial_data_all.allGroups{iAnimal}{iRecording}(fsi) = allDepths(fsi) + 3;
        trial_data_all.allGroups{iAnimal}{iRecording}(tan) = allDepths(tan) + 6;
        trial_data_all.allGroups{iAnimal}{iRecording}(uin) = allDepths(uin) + 9;
        trial_data_all.allGroups{iAnimal}{iRecording}(shortDurLongPss) = allDepths(shortDurLongPss) + 12;
        
        clearvars pss_allcat2
        end
    end
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
    
    if length(unique(trial_info_allcat.stimulus)) == 3 && ...
            all(unique(trial_info_allcat.stimulus) == [1;2;3])
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

% Concatenate good units
good_units_exp = cellfun(@transpose,vertcat(goodUnits{:}),'uni',false);
good_units_allcat = cell2mat(vertcat(goodUnits{:})')';

% (groups are cell type * depth: msn, fsi, tan, th, ??)
n_aligned_depths = 3; % hardcoded: I think not stored
n_celltypes = max(groups_allcat)./n_aligned_depths;
if n_celltypes ~= 5
    error('Incorrect number of celltypes')
end
celltype_allcat = ceil(groups_allcat./n_aligned_depths);
celltype_labels = {'MSN','FSI','TAN','UIN','Short bursty TAN-like'};

domain_allcat = mod(groups_allcat,n_aligned_depths) + ...
    n_aligned_depths*(mod(groups_allcat,n_aligned_depths) == 0);

% Get maps for all cells
% (load master U)
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');
% (use t = 0 kernels for all cells, hardcoded at the moment)
use_k_frame = 5;
ctx_str_k_frame = cell2mat(cellfun(@(x) reshape(x(:,use_k_frame,:),100,[]), ...
    vertcat(ctx_str_k_all{:})','uni',false));
ctx_str_k_px = svdFrameReconstruct(U_master(:,:,1:100),ctx_str_k_frame(:,good_units_allcat));


% Deconvolve fluorescence, baseline subtract, get kernel ROIs
fluor_deconv_allcat_exp = cellfun(@(x) AP_deconv_wf(x),vertcat(fluor_all{:}),'uni',false);

n_vs = size(fluor_all{end}{end},3);
kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
load(kernel_roi_fn);
fluor_kernelroi_deconv_exp = cellfun(@(x) ...
    permute(AP_deconv_wf(AP_svd_roi(U_master(:,:,1:n_vs), ...
    permute(x-nanmean(x(:,t < 0,:),2),[3,2,1]),[],[],kernel_roi.bw)),[3,2,1]), ...
    vertcat(fluor_all{:}),'uni',false);


% Get animals/days/cell number
% (note this includes cortical cells, which are NaNs in allGroups)
animals_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp)), ...
    vertcat(allAnimals{:}),vertcat(allGroups{:}),'uni',false));
days_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp)), ...
    vertcat(allDays{:}),vertcat(allGroups{:}),'uni',false));
neurons_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp))', ...
    vertcat(allNeurons{:}),vertcat(allGroups{:}),'uni',false));
recordings_allcat = cell2mat(cellfun(@(x,grp) x*ones(sum(~isnan(grp)),1), ...    
    num2cell((1:length(vertcat(allDays{:})))'),vertcat(allGroups{:}),'uni',false));


disp('Finished.')


%% Get fluorescence in kernel ROIs (put in above)

% Get number of widefield components 
n_vs = size(fluor_all{end}{end},3);

% Load kernel ROIs
kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
load(kernel_roi_fn);
n_kernel_rois = size(kernel_roi.bw,3);

% Get deconvolved/baseline-subtracted kernel ROI fluorescence
fluor_kernelroi_deconv_exp = cellfun(@(x) ...
    permute(AP_deconv_wf(AP_svd_roi(U_master(:,:,1:n_vs), ...
    permute(x-nanmean(x(:,t < 0,:),2),[3,2,1]),[],[],kernel_roi.bw)),[3,2,1]), ...
    vertcat(fluor_all{:}),'uni',false);



%% Align by depth and domain

%%% Align by relative depth 
ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
ephys_align_fn = ['ephys_depth_align.mat'];
load([ephys_align_path filesep ephys_align_fn]);

% Pull out relative depths for each cell
depth_aligned = cell(size(mua_all));
for curr_animal = 1:length(mua_all)
    curr_animal_idx = strcmp(animals{curr_animal},{ephys_depth_align.animal});
    for curr_day = 1:length(mua_all{curr_animal})
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
    depth_curr_animal_idx = strcmp(animals{curr_animal},{ephys_depth_align.animal});
    kernel_curr_animal_idx = strcmp(animals{curr_animal},{ephys_kernel_align.animal});
    for curr_day = 1:length(mua_all{curr_animal})
        
        str_depth = ephys_depth_align(depth_curr_animal_idx).str_depth(curr_day,:);
        kernel_match = ephys_kernel_align(kernel_curr_animal_idx).kernel_match{curr_day};

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

include_celltypes = [1:4]; % 5 is probably incorrectly split units
include_units = good_units_allcat & ...
    ismember(celltype_allcat,include_celltypes);

domain_type_count = accumarray([ ...
    domain_aligned_allcat(include_units), ...
    celltype_allcat(include_units)],1);

celltype_col = ...
    [0.9,0.4,0.6; ...
    0.4,0.5,0.6; ...
    0.5,0.3,0.1; ...
    1,0.5,0];

figure; hold on;
set(gca,'ColorOrder',celltype_col);
bar(domain_type_count,'stacked');
xlabel('Domain');
ylabel('Count');
legend(celltype_labels);

%% >>>> Final: Load task dataset


data_fn = { ...
    'trial_activity_choiceworld.mat', ... % original 
    'trial_activity_ctx_task.mat', ...    % + cortex ephys
    'trial_activity_muscimol_task.mat'};  % muscimol group

% (turn on warnings)
warning on;

% Load data (saved as structure trial_data_all)
trial_data_path = 'G:\JF_single_cell_data\';

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


% Julie's bugfixes
% (combines 2 TAN categories, integrates old 'other' with UINs, puts short
% waveforms that look like TANs into new 'other')
removeMissingSpikesUnits = true;
for iAnimal = 1:size(trial_data_all.goodUnits, 1)
    for iRecording = 1:size(trial_data_all.goodUnits{iAnimal}, 1)

        % 1. if flag remove spikes missing, add these to the goodUnits
        if removeMissingSpikesUnits
            trial_data_all.goodUnits{iAnimal}{iRecording} = trial_data_all.goodUnits{iAnimal}{iRecording} ...
                & nanmin(trial_data_all.percent_missing_ndtr{iAnimal}{iRecording}(2:end, :)) < 30;
        end

        % 2. re-classify to: add noise guys to uins, remove large post spike suppression FSI and UINs, and short waveform guys from TANs. 
        if ~isempty(trial_data_all.acg{iAnimal}{iRecording})
         for iCell = 1:size(trial_data_all.acg{iAnimal}{iRecording}, 1)
            pss_allcat2temp = find(trial_data_all.acg{iAnimal}{iRecording}(iCell, 500:1000) >= nanmean(trial_data_all.acg{iAnimal}{iRecording}(iCell, 600:900)));
            pss_allcat2(iCell) = pss_allcat2temp(1);
         end
        allDepths = zeros(size(trial_data_all.allGroups{iAnimal}{iRecording},2),1);
        allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [1,4,7,10,13,16]))=1; 
        allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [2,5,8,11,14,17]))=2;
        allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [3,6,9,12,15,18]))=3;
        largePssShorties = find(trial_data_all.templateDuration{iAnimal}{iRecording} < 400 & pss_allcat2' > 40); 
        
        fsi = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 < 40 & trial_data_all.prop_long_isi{iAnimal}{iRecording} <= 0.1;
        uin = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 < 40 & trial_data_all.prop_long_isi{iAnimal}{iRecording} > 0.1;
        tan = trial_data_all.templateDuration{iAnimal}{iRecording} > 400 & pss_allcat2 >= 40;
        msn = trial_data_all.templateDuration{iAnimal}{iRecording} > 400 & pss_allcat2 < 40;
        shortDurLongPss = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 >= 40;
        
        trial_data_all.allGroups{iAnimal}{iRecording}(msn) = allDepths(msn);
        trial_data_all.allGroups{iAnimal}{iRecording}(fsi) = allDepths(fsi) + 3;
        trial_data_all.allGroups{iAnimal}{iRecording}(tan) = allDepths(tan) + 6;
        trial_data_all.allGroups{iAnimal}{iRecording}(uin) = allDepths(uin) + 9;
        trial_data_all.allGroups{iAnimal}{iRecording}(shortDurLongPss) = allDepths(shortDurLongPss) + 12;
        
        clearvars pss_allcat2
        end
    end
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
    
    if length(unique(trial_info_allcat.stimulus)) == 3 && ...
            all(unique(trial_info_allcat.stimulus) == [1;2;3])
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

% Concatenate good units
good_units_exp = cellfun(@transpose,vertcat(goodUnits{:}),'uni',false);
good_units_allcat = cell2mat(vertcat(goodUnits{:})')';

% (groups are cell type * depth: msn, fsi, tan, th, ??)
n_aligned_depths = 3; % hardcoded: I think not stored
n_celltypes = max(groups_allcat)./n_aligned_depths;
if n_celltypes ~= 5
    error('Incorrect number of celltypes')
end
celltype_allcat = ceil(groups_allcat./n_aligned_depths);
celltype_labels = {'MSN','FSI','TAN','UIN','Short bursty TAN-like'};

domain_allcat = mod(groups_allcat,n_aligned_depths) + ...
    n_aligned_depths*(mod(groups_allcat,n_aligned_depths) == 0);

% Get maps for all cells
% (load master U)
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');
% (use t = 0 kernels for all cells, hardcoded at the moment)
use_k_frame = 5;
ctx_str_k_frame = cell2mat(cellfun(@(x) reshape(x(:,use_k_frame,:),100,[]), ...
    vertcat(ctx_str_k_all{:})','uni',false));
ctx_str_k_px = svdFrameReconstruct(U_master(:,:,1:100),ctx_str_k_frame(:,good_units_allcat));


% Deconvolve fluorescence and get kernel ROIs
fluor_deconv_allcat_exp = cellfun(@(x) AP_deconv_wf(x),vertcat(fluor_all{:}),'uni',false);

n_vs = size(fluor_all{end}{end},3);
kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
load(kernel_roi_fn);
fluor_kernelroi_deconv_exp = cellfun(@(x) ...
    permute(AP_deconv_wf(AP_svd_roi(U_master(:,:,1:n_vs), ...
    permute(x-nanmean(x(:,t < 0,:),2),[3,2,1]),[],[],kernel_roi.bw)),[3,2,1]), ...
    vertcat(fluor_all{:}),'uni',false);


% Get animals/days/cell number
% (note this includes cortical cells, which are NaNs in allGroups)
animals_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp)), ...
    vertcat(allAnimals{:}),vertcat(allGroups{:}),'uni',false));
days_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp)), ...
    vertcat(allDays{:}),vertcat(allGroups{:}),'uni',false));
neurons_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp))', ...
    vertcat(allNeurons{:}),vertcat(allGroups{:}),'uni',false));
recordings_allcat = cell2mat(cellfun(@(x,grp) x*ones(sum(~isnan(grp)),1), ...    
    num2cell((1:length(vertcat(allDays{:})))'),vertcat(allGroups{:}),'uni',false));


%%% Align by relative depth 
ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
ephys_align_fn = ['ephys_depth_align.mat'];
load([ephys_align_path filesep ephys_align_fn]);

% Pull out relative depths for each cell
depth_aligned = cell(size(mua_all));
for curr_animal = 1:length(mua_all)
    curr_animal_idx = strcmp(animals{curr_animal},{ephys_depth_align.animal});
    for curr_day = 1:length(mua_all{curr_animal})
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
    depth_curr_animal_idx = strcmp(animals{curr_animal},{ephys_depth_align.animal});
    kernel_curr_animal_idx = strcmp(animals{curr_animal},{ephys_kernel_align.animal});
    for curr_day = 1:length(mua_all{curr_animal})
        
        str_depth = ephys_depth_align(depth_curr_animal_idx).str_depth(curr_day,:);
        kernel_match = ephys_kernel_align(kernel_curr_animal_idx).kernel_match{curr_day};

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



disp('Finished.')


%% >>>> Final: Load passive naive dataset


data_fn = 'trial_activity_naive.mat';

% (turn on warnings)
warning on;

% Load data (saved as structure trial_data_all)
trial_data_path = 'G:\JF_single_cell_data\';

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


% Julie's bugfixes
% (combines 2 TAN categories, integrates old 'other' with UINs, puts short
% waveforms that look like TANs into new 'other')
removeMissingSpikesUnits = true;
for iAnimal = 1:size(trial_data_all.goodUnits, 1)
    for iRecording = 1:size(trial_data_all.goodUnits{iAnimal}, 1)

        % 1. if flag remove spikes missing, add these to the goodUnits
        if removeMissingSpikesUnits
            trial_data_all.goodUnits{iAnimal}{iRecording} = trial_data_all.goodUnits{iAnimal}{iRecording} ...
                & nanmin(trial_data_all.percent_missing_ndtr{iAnimal}{iRecording}(2:end, :)) < 30;
        end

        % 2. re-classify to: add noise guys to uins, remove large post spike suppression FSI and UINs, and short waveform guys from TANs. 
        if ~isempty(trial_data_all.acg{iAnimal}{iRecording})
         for iCell = 1:size(trial_data_all.acg{iAnimal}{iRecording}, 1)
            pss_allcat2temp = find(trial_data_all.acg{iAnimal}{iRecording}(iCell, 500:1000) >= nanmean(trial_data_all.acg{iAnimal}{iRecording}(iCell, 600:900)));
            pss_allcat2(iCell) = pss_allcat2temp(1);
         end
        allDepths = zeros(size(trial_data_all.allGroups{iAnimal}{iRecording},2),1);
        allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [1,4,7,10,13,16]))=1; 
        allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [2,5,8,11,14,17]))=2;
        allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [3,6,9,12,15,18]))=3;
        largePssShorties = find(trial_data_all.templateDuration{iAnimal}{iRecording} < 400 & pss_allcat2' > 40); 
        
        fsi = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 < 40 & trial_data_all.prop_long_isi{iAnimal}{iRecording} <= 0.1;
        uin = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 < 40 & trial_data_all.prop_long_isi{iAnimal}{iRecording} > 0.1;
        tan = trial_data_all.templateDuration{iAnimal}{iRecording} > 400 & pss_allcat2 >= 40;
        msn = trial_data_all.templateDuration{iAnimal}{iRecording} > 400 & pss_allcat2 < 40;
        shortDurLongPss = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 >= 40;
        
        trial_data_all.allGroups{iAnimal}{iRecording}(msn) = allDepths(msn);
        trial_data_all.allGroups{iAnimal}{iRecording}(fsi) = allDepths(fsi) + 3;
        trial_data_all.allGroups{iAnimal}{iRecording}(tan) = allDepths(tan) + 6;
        trial_data_all.allGroups{iAnimal}{iRecording}(uin) = allDepths(uin) + 9;
        trial_data_all.allGroups{iAnimal}{iRecording}(shortDurLongPss) = allDepths(shortDurLongPss) + 12;
        
        clearvars pss_allcat2
        end
    end
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

% BUGFIX? If any days don't have data - throw away
% (this isn't present in MUA dataset, don't know where this is from)
for curr_animal = 1:length(trial_data_all.animals)
    nodata_days = cellfun(@(x) isempty(x),trial_data_all.trial_info_all{curr_animal});
    if any(nodata_days)
        for curr_field = data_struct_fieldnames(experiment_fields)'
            trial_data_all.(cell2mat(curr_field)){curr_animal}(nodata_days) = [];
        end
    end
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
    
    if length(unique(trial_info_allcat.stimulus)) == 3 && ...
            all(unique(trial_info_allcat.stimulus) == [1;2;3])
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

% Concatenate good units
good_units_exp = cellfun(@transpose,vertcat(goodUnits{:}),'uni',false);
good_units_allcat = cell2mat(vertcat(goodUnits{:})')';

% (groups are cell type * depth: msn, fsi, tan, th, ??)
n_aligned_depths = 3; % hardcoded: I think not stored
n_celltypes = max(groups_allcat)./n_aligned_depths;
if n_celltypes ~= 5
    error('Incorrect number of celltypes')
end
celltype_allcat = ceil(groups_allcat./n_aligned_depths);
celltype_labels = {'MSN','FSI','TAN','UIN','Short bursty TAN-like'};

domain_allcat = mod(groups_allcat,n_aligned_depths) + ...
    n_aligned_depths*(mod(groups_allcat,n_aligned_depths) == 0);

% Get maps for all cells
% (load master U)
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');
% (use t = 0 kernels for all cells, hardcoded at the moment)
use_k_frame = 5;
ctx_str_k_frame = cell2mat(cellfun(@(x) reshape(x(:,use_k_frame,:),100,[]), ...
    vertcat(ctx_str_k_all{:})','uni',false));
ctx_str_k_px = svdFrameReconstruct(U_master(:,:,1:100),ctx_str_k_frame(:,good_units_allcat));


% Deconvolve fluorescence and get kernel ROIs
fluor_deconv_allcat_exp = cellfun(@(x) AP_deconv_wf(x),vertcat(fluor_all{:}),'uni',false);

n_vs = size(fluor_all{end}{end},3);
kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
load(kernel_roi_fn);
fluor_kernelroi_deconv_exp = cellfun(@(x) ...
    permute(AP_deconv_wf(AP_svd_roi(U_master(:,:,1:n_vs), ...
    permute(x-nanmean(x(:,t < 0,:),2),[3,2,1]),[],[],kernel_roi.bw)),[3,2,1]), ...
    vertcat(fluor_all{:}),'uni',false);


% Get animals/days/cell number
% (note this includes cortical cells, which are NaNs in allGroups)
animals_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp)), ...
    vertcat(allAnimals{:}),vertcat(allGroups{:}),'uni',false));
days_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp)), ...
    vertcat(allDays{:}),vertcat(allGroups{:}),'uni',false));
neurons_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp))', ...
    vertcat(allNeurons{:}),vertcat(allGroups{:}),'uni',false));
recordings_allcat = cell2mat(cellfun(@(x,grp) x*ones(sum(~isnan(grp)),1), ...    
    num2cell((1:length(vertcat(allDays{:})))'),vertcat(allGroups{:}),'uni',false));


%%% Align by relative depth 
ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
ephys_align_fn = ['ephys_depth_align.mat'];
load([ephys_align_path filesep ephys_align_fn]);

% Pull out relative depths for each cell
depth_aligned = cell(size(mua_all));
for curr_animal = 1:length(mua_all)
    curr_animal_idx = strcmp(animals{curr_animal},{ephys_depth_align.animal});
    for curr_day = 1:length(mua_all{curr_animal})
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
    depth_curr_animal_idx = strcmp(animals{curr_animal},{ephys_depth_align.animal});
    kernel_curr_animal_idx = strcmp(animals{curr_animal},{ephys_kernel_align.animal});
    for curr_day = 1:length(mua_all{curr_animal})
        
        str_depth = ephys_depth_align(depth_curr_animal_idx).str_depth(curr_day,:);
        kernel_match = ephys_kernel_align(kernel_curr_animal_idx).kernel_match{curr_day};

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



disp('Finished.')



%% >>>> Final: Load passive trained dataset


data_fn = { ...
    'trial_activity_trainedPassive.mat', ... % original 
    'trial_activity_ctx_passive.mat', ...    % + cortex ephys
    'trial_activity_muscimol_passive.mat'};  % muscimol group

% (turn on warnings)
warning on;

% Load data (saved as structure trial_data_all)
trial_data_path = 'G:\JF_single_cell_data\';

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


% Julie's bugfixes
% (combines 2 TAN categories, integrates old 'other' with UINs, puts short
% waveforms that look like TANs into new 'other')
removeMissingSpikesUnits = true;
for iAnimal = 1:size(trial_data_all.goodUnits, 1)
    for iRecording = 1:size(trial_data_all.goodUnits{iAnimal}, 1)

        % 1. if flag remove spikes missing, add these to the goodUnits
        if removeMissingSpikesUnits
            trial_data_all.goodUnits{iAnimal}{iRecording} = trial_data_all.goodUnits{iAnimal}{iRecording} ...
                & nanmin(trial_data_all.percent_missing_ndtr{iAnimal}{iRecording}(2:end, :)) < 30;
        end

        % 2. re-classify to: add noise guys to uins, remove large post spike suppression FSI and UINs, and short waveform guys from TANs. 
        if ~isempty(trial_data_all.acg{iAnimal}{iRecording})
         for iCell = 1:size(trial_data_all.acg{iAnimal}{iRecording}, 1)
            pss_allcat2temp = find(trial_data_all.acg{iAnimal}{iRecording}(iCell, 500:1000) >= nanmean(trial_data_all.acg{iAnimal}{iRecording}(iCell, 600:900)));
            pss_allcat2(iCell) = pss_allcat2temp(1);
         end
        allDepths = zeros(size(trial_data_all.allGroups{iAnimal}{iRecording},2),1);
        allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [1,4,7,10,13,16]))=1; 
        allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [2,5,8,11,14,17]))=2;
        allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [3,6,9,12,15,18]))=3;
        largePssShorties = find(trial_data_all.templateDuration{iAnimal}{iRecording} < 400 & pss_allcat2' > 40); 
        
        fsi = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 < 40 & trial_data_all.prop_long_isi{iAnimal}{iRecording} <= 0.1;
        uin = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 < 40 & trial_data_all.prop_long_isi{iAnimal}{iRecording} > 0.1;
        tan = trial_data_all.templateDuration{iAnimal}{iRecording} > 400 & pss_allcat2 >= 40;
        msn = trial_data_all.templateDuration{iAnimal}{iRecording} > 400 & pss_allcat2 < 40;
        shortDurLongPss = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 >= 40;
        
        trial_data_all.allGroups{iAnimal}{iRecording}(msn) = allDepths(msn);
        trial_data_all.allGroups{iAnimal}{iRecording}(fsi) = allDepths(fsi) + 3;
        trial_data_all.allGroups{iAnimal}{iRecording}(tan) = allDepths(tan) + 6;
        trial_data_all.allGroups{iAnimal}{iRecording}(uin) = allDepths(uin) + 9;
        trial_data_all.allGroups{iAnimal}{iRecording}(shortDurLongPss) = allDepths(shortDurLongPss) + 12;
        
        clearvars pss_allcat2
        end
    end
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

% BUGFIX? If any days don't have data - throw away
% (this isn't present in MUA dataset, don't know where this is from)
for curr_animal = 1:length(trial_data_all.animals)
    nodata_days = cellfun(@(x) isempty(x),trial_data_all.trial_info_all{curr_animal});
    if any(nodata_days)
        for curr_field = data_struct_fieldnames(experiment_fields)'
            trial_data_all.(cell2mat(curr_field)){curr_animal}(nodata_days) = [];
        end
    end
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
    
    if length(unique(trial_info_allcat.stimulus)) == 3 && ...
            all(unique(trial_info_allcat.stimulus) == [1;2;3])
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

% Concatenate good units
good_units_exp = cellfun(@transpose,vertcat(goodUnits{:}),'uni',false);
good_units_allcat = cell2mat(vertcat(goodUnits{:})')';

% (groups are cell type * depth: msn, fsi, tan, th, ??)
n_aligned_depths = 3; % hardcoded: I think not stored
n_celltypes = max(groups_allcat)./n_aligned_depths;
if n_celltypes ~= 5
    error('Incorrect number of celltypes')
end
celltype_allcat = ceil(groups_allcat./n_aligned_depths);
celltype_labels = {'MSN','FSI','TAN','UIN','Short bursty TAN-like'};

domain_allcat = mod(groups_allcat,n_aligned_depths) + ...
    n_aligned_depths*(mod(groups_allcat,n_aligned_depths) == 0);

% Get maps for all cells
% (load master U)
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');
% (use t = 0 kernels for all cells, hardcoded at the moment)
use_k_frame = 5;
ctx_str_k_frame = cell2mat(cellfun(@(x) reshape(x(:,use_k_frame,:),100,[]), ...
    vertcat(ctx_str_k_all{:})','uni',false));
ctx_str_k_px = svdFrameReconstruct(U_master(:,:,1:100),ctx_str_k_frame(:,good_units_allcat));


% Deconvolve fluorescence and get kernel ROIs
fluor_deconv_allcat_exp = cellfun(@(x) AP_deconv_wf(x),vertcat(fluor_all{:}),'uni',false);

n_vs = size(fluor_all{end}{end},3);
kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
load(kernel_roi_fn);
fluor_kernelroi_deconv_exp = cellfun(@(x) ...
    permute(AP_deconv_wf(AP_svd_roi(U_master(:,:,1:n_vs), ...
    permute(x-nanmean(x(:,t < 0,:),2),[3,2,1]),[],[],kernel_roi.bw)),[3,2,1]), ...
    vertcat(fluor_all{:}),'uni',false);


% Get animals/days/cell number
% (note this includes cortical cells, which are NaNs in allGroups)
animals_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp)), ...
    vertcat(allAnimals{:}),vertcat(allGroups{:}),'uni',false));
days_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp)), ...
    vertcat(allDays{:}),vertcat(allGroups{:}),'uni',false));
neurons_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp))', ...
    vertcat(allNeurons{:}),vertcat(allGroups{:}),'uni',false));
recordings_allcat = cell2mat(cellfun(@(x,grp) x*ones(sum(~isnan(grp)),1), ...    
    num2cell((1:length(vertcat(allDays{:})))'),vertcat(allGroups{:}),'uni',false));


%%% Align by relative depth 
ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
ephys_align_fn = ['ephys_depth_align.mat'];
load([ephys_align_path filesep ephys_align_fn]);

% Pull out relative depths for each cell
depth_aligned = cell(size(mua_all));
for curr_animal = 1:length(mua_all)
    curr_animal_idx = strcmp(animals{curr_animal},{ephys_depth_align.animal});
    for curr_day = 1:length(mua_all{curr_animal})
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
    depth_curr_animal_idx = strcmp(animals{curr_animal},{ephys_depth_align.animal});
    kernel_curr_animal_idx = strcmp(animals{curr_animal},{ephys_kernel_align.animal});
    for curr_day = 1:length(mua_all{curr_animal})
        
        str_depth = ephys_depth_align(depth_curr_animal_idx).str_depth(curr_day,:);
        kernel_match = ephys_kernel_align(kernel_curr_animal_idx).kernel_match{curr_day};

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



disp('Finished.')

%% >>>> Final: Load/plot passive naive/trained datasets

% data_fns = { ...
%     'trial_activity_naive_sua', ...
%     {'trial_activity_trainedPassive_sua.mat', ... % original
%     'trial_activity_ctx_passive_sua.mat', ...     % + cortex ephys
%     'trial_activity_muscimol_passive_sua.mat'}};  % muscimol group

% (leave out original group - 1s instead of 0.5s)
data_fns = { ...
    'trial_activity_naive_sua', ...
    'trial_activity_muscimol_passive_sua.mat'};  % muscimol group

stim_act_celltype_training = cell(size(data_fns));
stim_act_celltype_group = cell(size(data_fns));

for curr_data = 1:length(data_fns)
    
    preload_vars = who;
    
    data_fn = data_fns{curr_data};
    
    % (turn on warnings)
    warning on;
    
    % Load data (saved as structure trial_data_all)
    trial_data_path = 'G:\JF_single_cell_data\';
    
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
    
    
    % Julie's bugfixes
    % (combines 2 TAN categories, integrates old 'other' with UINs, puts short
    % waveforms that look like TANs into new 'other')
    removeMissingSpikesUnits = true;
    for iAnimal = 1:size(trial_data_all.goodUnits, 1)
        for iRecording = 1:size(trial_data_all.goodUnits{iAnimal}, 1)
            
            % 1. if flag remove spikes missing, add these to the goodUnits
            if removeMissingSpikesUnits
                trial_data_all.goodUnits{iAnimal}{iRecording} = trial_data_all.goodUnits{iAnimal}{iRecording} ...
                    & nanmin(trial_data_all.percent_missing_ndtr{iAnimal}{iRecording}(2:end, :)) < 30;
            end
            
            % 2. re-classify to: add noise guys to uins, remove large post spike suppression FSI and UINs, and short waveform guys from TANs.
            if ~isempty(trial_data_all.acg{iAnimal}{iRecording})
                for iCell = 1:size(trial_data_all.acg{iAnimal}{iRecording}, 1)
                    pss_allcat2temp = find(trial_data_all.acg{iAnimal}{iRecording}(iCell, 500:1000) >= nanmean(trial_data_all.acg{iAnimal}{iRecording}(iCell, 600:900)));
                    pss_allcat2(iCell) = pss_allcat2temp(1);
                end
                allDepths = zeros(size(trial_data_all.allGroups{iAnimal}{iRecording},2),1);
                allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [1,4,7,10,13,16]))=1;
                allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [2,5,8,11,14,17]))=2;
                allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [3,6,9,12,15,18]))=3;
                largePssShorties = find(trial_data_all.templateDuration{iAnimal}{iRecording} < 400 & pss_allcat2' > 40);
                
                fsi = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 < 40 & trial_data_all.prop_long_isi{iAnimal}{iRecording} <= 0.1;
                uin = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 < 40 & trial_data_all.prop_long_isi{iAnimal}{iRecording} > 0.1;
                tan = trial_data_all.templateDuration{iAnimal}{iRecording} > 400 & pss_allcat2 >= 40;
                msn = trial_data_all.templateDuration{iAnimal}{iRecording} > 400 & pss_allcat2 < 40;
                shortDurLongPss = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 >= 40;
                
                trial_data_all.allGroups{iAnimal}{iRecording}(msn) = allDepths(msn);
                trial_data_all.allGroups{iAnimal}{iRecording}(fsi) = allDepths(fsi) + 3;
                trial_data_all.allGroups{iAnimal}{iRecording}(tan) = allDepths(tan) + 6;
                trial_data_all.allGroups{iAnimal}{iRecording}(uin) = allDepths(uin) + 9;
                trial_data_all.allGroups{iAnimal}{iRecording}(shortDurLongPss) = allDepths(shortDurLongPss) + 12;
                
                clearvars pss_allcat2
            end
        end
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
    
    % BUGFIX? If any days don't have data - throw away
    % (this isn't present in MUA dataset, don't know where this is from)
    for curr_animal = 1:length(trial_data_all.animals)
        nodata_days = cellfun(@(x) isempty(x),trial_data_all.trial_info_all{curr_animal});
        if any(nodata_days)
            for curr_field = data_struct_fieldnames(experiment_fields)'
                trial_data_all.(cell2mat(curr_field)){curr_animal}(nodata_days) = [];
            end
        end
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
        
        if length(unique(trial_info_allcat.stimulus)) == 3 && ...
            all(unique(trial_info_allcat.stimulus) == [1;2;3])
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
    
    % Concatenate good units
    good_units_exp = cellfun(@transpose,vertcat(goodUnits{:}),'uni',false);
    good_units_allcat = cell2mat(vertcat(goodUnits{:})')';
    
    % (groups are cell type * depth: msn, fsi, tan, th, ??)
    n_aligned_depths = 3; % hardcoded: I think not stored
    n_celltypes = max(groups_allcat)./n_aligned_depths;
    if n_celltypes ~= 5
        error('Incorrect number of celltypes')
    end
    celltype_allcat = ceil(groups_allcat./n_aligned_depths);
    celltype_labels = {'MSN','FSI','TAN','UIN','Short bursty TAN-like'};
    
    domain_allcat = mod(groups_allcat,n_aligned_depths) + ...
        n_aligned_depths*(mod(groups_allcat,n_aligned_depths) == 0);
    
    % Get maps for all cells
    % (load master U)
    load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');
    % (use t = 0 kernels for all cells, hardcoded at the moment)
    use_k_frame = 5;
    ctx_str_k_frame = cell2mat(cellfun(@(x) reshape(x(:,use_k_frame,:),100,[]), ...
        vertcat(ctx_str_k_all{:})','uni',false));
    ctx_str_k_px = svdFrameReconstruct(U_master(:,:,1:100),ctx_str_k_frame(:,good_units_allcat));
    
    
    % Deconvolve fluorescence and get kernel ROIs
    fluor_deconv_allcat_exp = cellfun(@(x) AP_deconv_wf(x),vertcat(fluor_all{:}),'uni',false);
    
    n_vs = size(fluor_all{end}{end},3);
    kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
    load(kernel_roi_fn);
    fluor_kernelroi_deconv_exp = cellfun(@(x) ...
        permute(AP_deconv_wf(AP_svd_roi(U_master(:,:,1:n_vs), ...
        permute(x-nanmean(x(:,t < 0,:),2),[3,2,1]),[],[],kernel_roi.bw)),[3,2,1]), ...
        vertcat(fluor_all{:}),'uni',false);
    
    
    % Get animals/days/cell number
    % (note this includes cortical cells, which are NaNs in allGroups)
    animals_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp)), ...
        vertcat(allAnimals{:}),vertcat(allGroups{:}),'uni',false));
    days_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp)), ...
        vertcat(allDays{:}),vertcat(allGroups{:}),'uni',false));
    neurons_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp))', ...
        vertcat(allNeurons{:}),vertcat(allGroups{:}),'uni',false));
    recordings_allcat = cell2mat(cellfun(@(x,grp) x*ones(sum(~isnan(grp)),1), ...
        num2cell((1:length(vertcat(allDays{:})))'),vertcat(allGroups{:}),'uni',false));
    
    
    %%% Align by relative depth
    ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
    ephys_align_fn = ['ephys_depth_align.mat'];
    load([ephys_align_path filesep ephys_align_fn]);
    
    % Pull out relative depths for each cell
    depth_aligned = cell(size(mua_all));
    for curr_animal = 1:length(mua_all)
        curr_animal_idx = strcmp(animals{curr_animal},{ephys_depth_align.animal});
        for curr_day = 1:length(mua_all{curr_animal})
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
        depth_curr_animal_idx = strcmp(animals{curr_animal},{ephys_depth_align.animal});
        kernel_curr_animal_idx = strcmp(animals{curr_animal},{ephys_kernel_align.animal});
        for curr_day = 1:length(mua_all{curr_animal})
            
            str_depth = ephys_depth_align(depth_curr_animal_idx).str_depth(curr_day,:);
            kernel_match = ephys_kernel_align(kernel_curr_animal_idx).kernel_match{curr_day};
            
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
    
       
    %%%%% Average stim activity
    
    % Split data by experiment
    mua_exp = vertcat(mua_all{:});
    wheel_exp = vertcat(wheel_all{:});
    stim_exp = mat2cell(trial_stim_allcat,use_split,1);
    
    % Get quiescent trials
    wheel_thresh = 0.025;
    quiescent_trials = cellfun(@(wheel) ...
        ~any(abs(wheel(:,t >= -0.5 & t <= 0.5)) > wheel_thresh,2), ...
        wheel_exp,'uni',false);
    
    % Baseline subtract spikes
    mua_exp_baselinesub = cellfun(@(act,quiescent) act - ...
        nanmean(reshape(act(quiescent,t < 0,:),[],1,size(act,3)),1), ...
        mua_exp,quiescent_trials,'uni',false);   
    %%%%%%%%%%%% TESTING;
%     mua_exp_baselinesub = mua_exp; warning('TESTING no baselinesub');
    
    % Get average activity with no wheel movement and 100%R stim
    stim_act = cellfun(@(stim,quiescent,act) ...
        permute(nanmean(act(stim == 1 & quiescent,:,:),1), ...
        [3,2,1]),stim_exp,quiescent_trials,mua_exp_baselinesub,'uni',false);
    
    % Get average activity for each celltype/domain/recording
    [stim_act_celltype,group_char] = grpstats(cell2mat(stim_act), ...
        [celltype_allcat,domain_aligned_allcat,recordings_allcat],...
        {'mean','gname'});
    
    group = cellfun(@str2num,group_char);
    
    stim_act_celltype_training{curr_data} = stim_act_celltype;
    stim_act_celltype_group{curr_data} = group;
    
    clearvars('-except',preload_vars{:},'t','n_aligned_depths','n_celltypes');
    
end

% Plot average stim timecourse before/after training
figure;
p = reshape(tight_subplot(n_aligned_depths,n_celltypes),[n_celltypes,n_aligned_depths])';
training_col = [0,0,0;1,0,0];
for curr_training = 1:2
    for curr_depth = 1:n_aligned_depths
        for curr_celltype = 1:n_celltypes
            
            curr_data = stim_act_celltype_training{curr_training}( ...
                all(stim_act_celltype_group{curr_training}(:,1:2) == ...
                [curr_celltype,curr_depth],2),:);
            
            axes(p(curr_depth,curr_celltype));
            AP_errorfill(t',nanmean(curr_data,1)', ...
                AP_sem(curr_data,1)',training_col(curr_training,:));
            
            xlim([-0.2,1])
            
        end
    end
end

for curr_celltype = 1:n_celltypes
    linkaxes(p(:,curr_celltype),'y');
end

% Get average stim activity before/after training
stim_avg_t = [0,0.2];
stim_avg_t_idx = t >= stim_avg_t(1) & t <= stim_avg_t(2);
stim_act_celltype_training_avg = cellfun(@(x) ...
    nanmean(x(:,stim_avg_t_idx),2),stim_act_celltype_training,'uni',false);

stim_act_avg_grp = cellfun(@(grp,act) ...
    accumarray(grp,act,max(grp,[],1),@nanmean,NaN), ...
    stim_act_celltype_group,stim_act_celltype_training_avg,'uni',false);


% (statistics)
disp('Pre/post ranksum (MSN,FSI,TAN,UIN):');
for curr_depth = 1:3
    curr_p = nan(1,4);
    for curr_celltype = 1:4
        curr_p(curr_celltype) = ...
            ranksum(reshape(stim_act_avg_grp{1}(curr_celltype,curr_depth,:),[],1), ...
            reshape(stim_act_avg_grp{2}(curr_celltype,curr_depth,:),[],1));
    end
    disp(['Str ' num2str(curr_depth) ': ' num2str(curr_p)]);
end


% Plot change by celltype
celltype_col = ...
    [0.9,0.4,0.6; ...
    0.4,0.5,0.6; ...
    0.5,0.3,0.1; ...
    1,0.5,0];

figure;
for curr_depth = 1:n_aligned_depths
    for curr_celltype = 1:4
        curr_untrained = squeeze(stim_act_avg_grp{1}(curr_celltype,curr_depth,:));
        curr_trained = squeeze(stim_act_avg_grp{2}(curr_celltype,curr_depth,:));
        
        curr_data = cellfun(@(x) squeeze(x(curr_celltype,curr_depth,:)), ...
            stim_act_avg_grp,'uni',false);
        
        subplot(n_aligned_depths,4,(curr_depth-1)*4+curr_celltype);
        plotSpread(curr_data,'distributionColors',min(celltype_col(curr_celltype,:)+0.2,1))
        errorbar(cellfun(@nanmean,curr_data),cellfun(@AP_sem,curr_data), ...
            'color',celltype_col(curr_celltype,:),'linewidth',2);
        set(gca,'XTickLabel',{'Untrained','Trained'}, ...
            'XTickLabelRotation',45)
        ylabel('Spikes/s');
    end
end



figure;
for curr_depth = 1:n_aligned_depths
    subplot(n_aligned_depths,1,curr_depth); hold on;
    
    for curr_celltype = 1:4
        curr_untrained = squeeze(stim_act_avg_grp{1}(curr_celltype,curr_depth,:));
        curr_trained = squeeze(stim_act_avg_grp{2}(curr_celltype,curr_depth,:));
        
        curr_data = cellfun(@(x) squeeze(x(curr_celltype,curr_depth,:)), ...
            stim_act_avg_grp,'uni',false);
                
        errorbar(cellfun(@nanmean,curr_data),cellfun(@AP_sem,curr_data), ...
            'color',celltype_col(curr_celltype,:),'linewidth',2);
    end
end


figure;
for curr_depth = 1:n_aligned_depths
    subplot(n_aligned_depths,1,curr_depth); hold on;
    
    for curr_celltype = 1:4
        curr_untrained = squeeze(stim_act_avg_grp{1}(curr_celltype,curr_depth,:));
        curr_trained = squeeze(stim_act_avg_grp{2}(curr_celltype,curr_depth,:));
        
        curr_data = cellfun(@(x) squeeze(x(curr_celltype,curr_depth,:)), ...
            stim_act_avg_grp,'uni',false);
        
        curr_data_norm = cellfun(@(x) x./nanmean(curr_data{1}),curr_data,'uni',false);
        
        errorbar(cellfun(@nanmean,curr_data_norm),cellfun(@AP_sem,curr_data_norm), ...
            'color',celltype_col(curr_celltype,:),'linewidth',2);
    end
end


%% >>>> TEST: Load/plot passive pre/post muscimol

data_fns = { ...
    'trial_activity_muscimol_passive_sua', ...
    'trial_activity_postmuscimol_passive_sua'};

% % (leave out original group - 1s instead of 0.5s)
% data_fns = { ...
%     'trial_activity_naive', ...
%     {'trial_activity_ctx_passive.mat', ...     % + cortex ephys
%     'trial_activity_muscimol_passive.mat'}};  % muscimol group

stim_act_celltype_training = cell(size(data_fns));
stim_act_celltype_group = cell(size(data_fns));

for curr_data = 1:length(data_fns)
    
    preload_vars = who;
    
    data_fn = data_fns{curr_data};
    
    % (turn on warnings)
    warning on;
    
    % Load data (saved as structure trial_data_all)
    trial_data_path = 'G:\JF_single_cell_data\';
    
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
    
    
    % Julie's bugfixes
    % (combines 2 TAN categories, integrates old 'other' with UINs, puts short
    % waveforms that look like TANs into new 'other')
    removeMissingSpikesUnits = true;
    for iAnimal = 1:size(trial_data_all.goodUnits, 1)
        for iRecording = 1:size(trial_data_all.goodUnits{iAnimal}, 1)
            
            % 1. if flag remove spikes missing, add these to the goodUnits
            if removeMissingSpikesUnits
                trial_data_all.goodUnits{iAnimal}{iRecording} = trial_data_all.goodUnits{iAnimal}{iRecording} ...
                    & nanmin(trial_data_all.percent_missing_ndtr{iAnimal}{iRecording}(2:end, :)) < 30;
            end
            
            % 2. re-classify to: add noise guys to uins, remove large post spike suppression FSI and UINs, and short waveform guys from TANs.
            if ~isempty(trial_data_all.acg{iAnimal}{iRecording})
                for iCell = 1:size(trial_data_all.acg{iAnimal}{iRecording}, 1)
                    pss_allcat2temp = find(trial_data_all.acg{iAnimal}{iRecording}(iCell, 500:1000) >= nanmean(trial_data_all.acg{iAnimal}{iRecording}(iCell, 600:900)));
                    pss_allcat2(iCell) = pss_allcat2temp(1);
                end
                allDepths = zeros(size(trial_data_all.allGroups{iAnimal}{iRecording},2),1);
                allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [1,4,7,10,13,16]))=1;
                allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [2,5,8,11,14,17]))=2;
                allDepths(ismember(trial_data_all.allGroups{iAnimal}{iRecording}, [3,6,9,12,15,18]))=3;
                largePssShorties = find(trial_data_all.templateDuration{iAnimal}{iRecording} < 400 & pss_allcat2' > 40);
                
                fsi = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 < 40 & trial_data_all.prop_long_isi{iAnimal}{iRecording} <= 0.1;
                uin = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 < 40 & trial_data_all.prop_long_isi{iAnimal}{iRecording} > 0.1;
                tan = trial_data_all.templateDuration{iAnimal}{iRecording} > 400 & pss_allcat2 >= 40;
                msn = trial_data_all.templateDuration{iAnimal}{iRecording} > 400 & pss_allcat2 < 40;
                shortDurLongPss = trial_data_all.templateDuration{iAnimal}{iRecording} <= 400 & pss_allcat2 >= 40;
                
                trial_data_all.allGroups{iAnimal}{iRecording}(msn) = allDepths(msn);
                trial_data_all.allGroups{iAnimal}{iRecording}(fsi) = allDepths(fsi) + 3;
                trial_data_all.allGroups{iAnimal}{iRecording}(tan) = allDepths(tan) + 6;
                trial_data_all.allGroups{iAnimal}{iRecording}(uin) = allDepths(uin) + 9;
                trial_data_all.allGroups{iAnimal}{iRecording}(shortDurLongPss) = allDepths(shortDurLongPss) + 12;
                
                clearvars pss_allcat2
            end
        end
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
    
    % BUGFIX? If any days don't have data - throw away
    % (this isn't present in MUA dataset, don't know where this is from)
    for curr_animal = 1:length(trial_data_all.animals)
        nodata_days = cellfun(@(x) isempty(x),trial_data_all.trial_info_all{curr_animal});
        if any(nodata_days)
            for curr_field = data_struct_fieldnames(experiment_fields)'
                trial_data_all.(cell2mat(curr_field)){curr_animal}(nodata_days) = [];
            end
        end
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
        
        if length(unique(trial_info_allcat.stimulus)) == 3 && ...
            all(unique(trial_info_allcat.stimulus) == [1;2;3])
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
    
    % Concatenate good units
    good_units_exp = cellfun(@transpose,vertcat(goodUnits{:}),'uni',false);
    good_units_allcat = cell2mat(vertcat(goodUnits{:})')';
    
    % (groups are cell type * depth: msn, fsi, tan, th, ??)
    n_aligned_depths = 3; % hardcoded: I think not stored
    n_celltypes = max(groups_allcat)./n_aligned_depths;
    if n_celltypes ~= 5
        error('Incorrect number of celltypes')
    end
    celltype_allcat = ceil(groups_allcat./n_aligned_depths);
    celltype_labels = {'MSN','FSI','TAN','UIN','Short bursty TAN-like'};
    
    domain_allcat = mod(groups_allcat,n_aligned_depths) + ...
        n_aligned_depths*(mod(groups_allcat,n_aligned_depths) == 0);
    
    % Get maps for all cells
    % (load master U)
    load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');
    % (use t = 0 kernels for all cells, hardcoded at the moment)
    use_k_frame = 5;
    ctx_str_k_frame = cell2mat(cellfun(@(x) reshape(x(:,use_k_frame,:),100,[]), ...
        vertcat(ctx_str_k_all{:})','uni',false));
    ctx_str_k_px = svdFrameReconstruct(U_master(:,:,1:100),ctx_str_k_frame(:,good_units_allcat));
    
    
    % Deconvolve fluorescence and get kernel ROIs
    fluor_deconv_allcat_exp = cellfun(@(x) AP_deconv_wf(x),vertcat(fluor_all{:}),'uni',false);
    
    n_vs = size(fluor_all{end}{end},3);
    kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
    load(kernel_roi_fn);
    fluor_kernelroi_deconv_exp = cellfun(@(x) ...
        permute(AP_deconv_wf(AP_svd_roi(U_master(:,:,1:n_vs), ...
        permute(x-nanmean(x(:,t < 0,:),2),[3,2,1]),[],[],kernel_roi.bw)),[3,2,1]), ...
        vertcat(fluor_all{:}),'uni',false);
    
    
    % Get animals/days/cell number
    % (note this includes cortical cells, which are NaNs in allGroups)
    animals_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp)), ...
        vertcat(allAnimals{:}),vertcat(allGroups{:}),'uni',false));
    days_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp)), ...
        vertcat(allDays{:}),vertcat(allGroups{:}),'uni',false));
    neurons_allcat = cell2mat(cellfun(@(x,grp) x(~isnan(grp))', ...
        vertcat(allNeurons{:}),vertcat(allGroups{:}),'uni',false));
    recordings_allcat = cell2mat(cellfun(@(x,grp) x*ones(sum(~isnan(grp)),1), ...
        num2cell((1:length(vertcat(allDays{:})))'),vertcat(allGroups{:}),'uni',false));
    
    
    %%% Align by relative depth
    ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
    ephys_align_fn = ['ephys_depth_align.mat'];
    load([ephys_align_path filesep ephys_align_fn]);
    
    % Pull out relative depths for each cell
    depth_aligned = cell(size(mua_all));
    for curr_animal = 1:length(mua_all)
        curr_animal_idx = strcmp(animals{curr_animal},{ephys_depth_align.animal});
        for curr_day = 1:length(mua_all{curr_animal})
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
        depth_curr_animal_idx = strcmp(animals{curr_animal},{ephys_depth_align.animal});
        kernel_curr_animal_idx = strcmp(animals{curr_animal},{ephys_kernel_align.animal});
        for curr_day = 1:length(mua_all{curr_animal})
            
            str_depth = ephys_depth_align(depth_curr_animal_idx).str_depth(curr_day,:);
            kernel_match = ephys_kernel_align(kernel_curr_animal_idx).kernel_match{curr_day};
            
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
    
       
    %%%%% Average stim activity
    
    % Split data by experiment
    mua_exp = vertcat(mua_all{:});
    wheel_exp = vertcat(wheel_all{:});
    stim_exp = mat2cell(trial_stim_allcat,use_split,1);
    
    % Get quiescent trials
    wheel_thresh = 0.025;
    quiescent_trials = cellfun(@(wheel) ...
        ~any(abs(wheel(:,t >= -0.5 & t <= 0.5)) > wheel_thresh,2), ...
        wheel_exp,'uni',false);
    
    % Baseline subtract spikes
    mua_exp_baselinesub = cellfun(@(act,quiescent) act - ...
        nanmean(reshape(act(quiescent,t < 0,:),[],1,size(act,3)),1), ...
        mua_exp,quiescent_trials,'uni',false);   
    %%%%%%%%%%%% TESTING;
%     mua_exp_baselinesub = mua_exp; warning('TESTING no baselinesub');
    
    % Get average activity with no wheel movement and 100%R stim
    stim_act = cellfun(@(stim,quiescent,act) ...
        permute(nanmean(act(stim == 1 & quiescent,:,:),1), ...
        [3,2,1]),stim_exp,quiescent_trials,mua_exp_baselinesub,'uni',false);
    
    % Get average activity for each celltype/domain/recording
    [stim_act_celltype,group_char] = grpstats(cell2mat(stim_act), ...
        [celltype_allcat,domain_aligned_allcat,recordings_allcat],...
        {'mean','gname'});
    
    group = cellfun(@str2num,group_char);
    
    stim_act_celltype_training{curr_data} = stim_act_celltype;
    stim_act_celltype_group{curr_data} = group;
    
    clearvars('-except',preload_vars{:},'t','n_aligned_depths','n_celltypes');
    
end

% Plot average stim timecourse before/after training
figure;
p = reshape(tight_subplot(n_aligned_depths,n_celltypes),[n_celltypes,n_aligned_depths])';
training_col = [0,0,0;1,0,0];
for curr_training = 1:2
    for curr_depth = 1:n_aligned_depths
        for curr_celltype = 1:n_celltypes
            
            curr_data = stim_act_celltype_training{curr_training}( ...
                all(stim_act_celltype_group{curr_training}(:,1:2) == ...
                [curr_celltype,curr_depth],2),:);
            
            axes(p(curr_depth,curr_celltype));
            AP_errorfill(t',nanmean(curr_data,1)', ...
                AP_sem(curr_data,1)',training_col(curr_training,:));
            
            xlim([-0.2,1])
            
        end
    end
end

for curr_celltype = 1:n_celltypes
    linkaxes(p(:,curr_celltype),'y');
end

% Get average stim activity before/after training
stim_avg_t = [0,0.2];
stim_avg_t_idx = t >= stim_avg_t(1) & t <= stim_avg_t(2);
stim_act_celltype_training_avg = cellfun(@(x) ...
    nanmean(x(:,stim_avg_t_idx),2),stim_act_celltype_training,'uni',false);

stim_act_avg_grp = cellfun(@(grp,act) ...
    accumarray(grp,act,max(grp,[],1),@nanmean,NaN), ...
    stim_act_celltype_group,stim_act_celltype_training_avg,'uni',false);


% (statistics)
disp('Pre/post ranksum (MSN,FSI,TAN,UIN):');
for curr_depth = 1:3
    curr_p = nan(1,4);
    for curr_celltype = 1:4
        curr_p(curr_celltype) = ...
            ranksum(reshape(stim_act_avg_grp{1}(curr_celltype,curr_depth,:),[],1), ...
            reshape(stim_act_avg_grp{2}(curr_celltype,curr_depth,:),[],1));
    end
    disp(['Str ' num2str(curr_depth) ': ' num2str(curr_p)]);
end


% Plot change by celltype
celltype_col = ...
    [0.9,0.4,0.6; ...
    0.4,0.5,0.6; ...
    0.5,0.3,0.1; ...
    1,0.5,0];

figure;
for curr_depth = 1:n_aligned_depths
    for curr_celltype = 1:4
        curr_untrained = squeeze(stim_act_avg_grp{1}(curr_celltype,curr_depth,:));
        curr_trained = squeeze(stim_act_avg_grp{2}(curr_celltype,curr_depth,:));
        
        curr_data = cellfun(@(x) squeeze(x(curr_celltype,curr_depth,:)), ...
            stim_act_avg_grp,'uni',false);
        
        subplot(n_aligned_depths,4,(curr_depth-1)*4+curr_celltype);
        plotSpread(curr_data,'distributionColors',min(celltype_col(curr_celltype,:)+0.2,1))
        errorbar(cellfun(@nanmean,curr_data),cellfun(@AP_sem,curr_data), ...
            'color',celltype_col(curr_celltype,:),'linewidth',2);
        set(gca,'XTickLabel',{'Untrained','Trained'}, ...
            'XTickLabelRotation',45)
        ylabel('Spikes/s');
    end
end



figure;
for curr_depth = 1:n_aligned_depths
    subplot(n_aligned_depths,1,curr_depth); hold on;
    
    for curr_celltype = 1:4
        curr_untrained = squeeze(stim_act_avg_grp{1}(curr_celltype,curr_depth,:));
        curr_trained = squeeze(stim_act_avg_grp{2}(curr_celltype,curr_depth,:));
        
        curr_data = cellfun(@(x) squeeze(x(curr_celltype,curr_depth,:)), ...
            stim_act_avg_grp,'uni',false);
                
        errorbar(cellfun(@nanmean,curr_data),cellfun(@AP_sem,curr_data), ...
            'color',celltype_col(curr_celltype,:),'linewidth',2);
    end
end


figure;
for curr_depth = 1:n_aligned_depths
    subplot(n_aligned_depths,1,curr_depth); hold on;
    
    for curr_celltype = 1:4
        curr_untrained = squeeze(stim_act_avg_grp{1}(curr_celltype,curr_depth,:));
        curr_trained = squeeze(stim_act_avg_grp{2}(curr_celltype,curr_depth,:));
        
        curr_data = cellfun(@(x) squeeze(x(curr_celltype,curr_depth,:)), ...
            stim_act_avg_grp,'uni',false);
        
        curr_data_norm = cellfun(@(x) x./nanmean(curr_data{1}),curr_data,'uni',false);
        
        errorbar(cellfun(@nanmean,curr_data_norm),cellfun(@AP_sem,curr_data_norm), ...
            'color',celltype_col(curr_celltype,:),'linewidth',2);
    end
end



%% ~~~~~~~~~~~~ Activity

%% Average activity of each celltype aligned by depth

mua_allcat_exp = vertcat(mua_all{:});

trial_stim_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
move_t_exp = mat2cell(move_t,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

act_mean = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:,:),1)), ...
    mua_allcat_exp,trial_stim_allcat_exp,move_t_exp, ...
    trial_outcome_allcat_exp,'uni',false)')';

figure;

subplot(1,3,1);
curr_cells_idx = find(good_units_allcat & ismember(celltype_allcat,[1]));
[~,sort_idx] = sort(depth_aligned_allcat(curr_cells_idx));
imagesc(zscore(act_mean(curr_cells_idx(sort_idx),:),[],2));
caxis([-3,3]);
colormap(brewermap([],'*RdBu'));
title('MSN');

subplot(1,3,2);
curr_cells_idx = find(good_units_allcat & ismember(celltype_allcat,[2]));
[~,sort_idx] = sort(depth_aligned_allcat(curr_cells_idx));
imagesc(zscore(act_mean(curr_cells_idx(sort_idx),:),[],2));
caxis([-3,3]);
colormap(brewermap([],'*RdBu'));
title('FSI');

subplot(1,3,3);
curr_cells_idx = find(good_units_allcat & ismember(celltype_allcat,[3,6]));
[~,sort_idx] = sort(depth_aligned_allcat(curr_cells_idx));
imagesc(zscore(act_mean(curr_cells_idx(sort_idx),:),[],2));
caxis([-3,3]);
colormap(brewermap([],'*RdBu'));
title('TAN');


%% Trial-average stim-aligned activity across cells (by domain/celltype)

mua_exp = vertcat(mua_all{:});

trial_stim_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
move_t_exp = mat2cell(move_t,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

 % Split trials in 2 using random vector to sort/plot
trial_rand = rand(max(cellfun(@(x) size(x,1),mua_exp)),1);

act_mean_sort = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1 & ...
    trial_rand(1:size(act,1)) < 0.5,:,:),1)), ...
    mua_exp,trial_stim_allcat_exp,move_t_exp, ...
    trial_outcome_allcat_exp,'uni',false)')';

act_mean_plot = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1 & ...
    trial_rand(1:size(act,1)) >= 0.5,:,:),1)), ...
    mua_exp,trial_stim_allcat_exp,move_t_exp, ...
    trial_outcome_allcat_exp,'uni',false)')';

% Plot average within each domain cell type
figure;
for curr_depth = 1:n_aligned_depths
    for curr_celltype = 1:n_celltypes
        
        curr_cells = domain_aligned_allcat == curr_depth & ...
            celltype_allcat == curr_celltype & good_units_allcat & ...
            any(act_mean_plot,2);
        curr_cells_idx = find(curr_cells);
        [~,max_idx] = max(act_mean_sort(curr_cells,:),[],2);
        [~,sort_idx] = sort(max_idx);
        
        curr_act_sorted = act_mean_plot(curr_cells_idx(sort_idx),:);
        curr_act_sorted_norm = curr_act_sorted./max(curr_act_sorted,[],2);
        
        % smooth MSNs across cells: too many to plot accurately
        if curr_celltype == 1
            n_smooth = 5;
            curr_act_sorted_norm = convn(curr_act_sorted_norm, ...
                ones(n_smooth,1)./n_smooth,'same');
        end
        

        subplot(n_aligned_depths,n_celltypes, ...
            (curr_depth-1)*n_celltypes+curr_celltype);
        imagesc(t,[],curr_act_sorted_norm);
        colormap(brewermap([],'Greys'));
        caxis([0,max(abs(caxis))]);
        axis off
        title(['Str ' num2str(curr_depth) ' ' celltype_labels{curr_celltype}]);        
        
    end
end

%% Concatenate trial-average and multi-aligned activity

% Concatenate multi-aligned trial-averaged activity
mua_exp = vertcat(mua_all{:});

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));

% Align activity by move/outcome (natively stim-aligned)

% (get indicies for alignments - used later)
stim_align = zeros(size(trial_stim_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;
use_align = {stim_align,move_align,outcome_align};

% (move aligned)
move_idx_exp = mat2cell(move_idx,use_split,1);
mua_exp_movealign = vertcat(mua_all{:});
for curr_exp = 1:length(mua_exp_movealign)
   for curr_trial = 1:size(mua_exp_movealign{curr_exp},1)
       mua_exp_movealign{curr_exp}(curr_trial,:,:) = ...
           circshift(mua_exp_movealign{curr_exp}(curr_trial,:,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
   end
end

% (outcome aligned)
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
mua_exp_outcomealign = vertcat(mua_all{:});
for curr_exp = 1:length(mua_exp_outcomealign)
    for curr_trial = 1:size(mua_exp_outcomealign{curr_exp},1)
        mua_exp_outcomealign{curr_exp}(curr_trial,:,:) = ...
            circshift(mua_exp_outcomealign{curr_exp}(curr_trial,:,:), ...
            -outcome_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
    end
end

% Get average activity across alignments
use_align_labels = {'Stim','Move onset','Outcome'};

act_mean = nan(length(good_units_allcat),length(t),length(use_align_labels));
for curr_align = 1:length(use_align_labels)
    
    switch curr_align
        case 1
            curr_mua = mua_exp;
        case 2
            curr_mua = mua_exp_movealign;
        case 3
            curr_mua = mua_exp_outcomealign;
    end
    
    act_mean(:,:,curr_align) = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
        squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:,:),1)), ...
        curr_mua,trial_stim_allcat_exp,move_t_exp, ...
        trial_outcome_allcat_exp,'uni',false)')';
    
end

% Concatenate multi-aligned averages
% (split the alignment halfway between median alignment points)
align_median = cellfun(@(x) -nanmedian(x)/sample_rate,use_align);
align_break = align_median(1:end-1) + diff(align_median*0.8);
align_t = {[-0.05,align_break(1)],[align_break(1:2)],[align_break(2),1]};

align_samples = cell(size(align_t));
for curr_align = 1:length(align_samples)    
    curr_t_offset = -nanmedian(use_align{curr_align})/sample_rate;
    curr_t = t + curr_t_offset;
    curr_t_use = curr_t >= align_t{curr_align}(1) & ...
        curr_t <= align_t{curr_align}(2);
    align_samples{curr_align} = curr_t_use;
end

act_mean_multialign = cell2mat(arrayfun(@(x) ...
    act_mean(:,align_samples{x},x),1:size(act_mean,3),'uni',false));


% Get pairwise across-recording activity correlations by celltype
% (this is a super dirty way to do it but I'm in a hurry)
n_recordings = max(recordings_allcat);
celltype_act_corr = cell(3,n_recordings);
celltype_act_minrate = cell(3,n_recordings);
for curr_depth = 1:n_aligned_depths
    for curr_celltype = 1:3
        
        % Get all pairwise activity correlations
        curr_cells = ...
            domain_aligned_allcat == curr_depth & ...
            celltype_allcat == curr_celltype;        
        curr_act_corr = corr(act_mean_multialign(curr_cells,:)');
        
        curr_rate = nanmean(act_mean_multialign(curr_cells,:),2);
        curr_act_minrate = min(curr_rate,curr_rate');
        
        for curr_recording = 1:n_recordings
            
            % Store
            if any(recordings_allcat(curr_cells) == curr_recording)
                curr_act_corr_use = ...
                    curr_act_corr(recordings_allcat(curr_cells) == curr_recording, ...
                    recordings_allcat(curr_cells) ~= curr_recording);
                
                curr_act_minrate_use = ...
                    curr_act_minrate(recordings_allcat(curr_cells) == curr_recording, ...
                    recordings_allcat(curr_cells) ~= curr_recording);
            else
                curr_act_corr_use = NaN;
                curr_act_minrate_use = NaN;
            end
         
            for curr_depth_compare = 1:2
                celltype_act_corr{curr_celltype,curr_recording} = ...
                   cat(1,celltype_act_corr{curr_celltype,curr_recording}, ...
                   curr_act_corr_use(:));
               
               celltype_act_minrate{curr_celltype,curr_recording} = ...
                   cat(1,celltype_act_minrate{curr_celltype,curr_recording}, ...
                   curr_act_minrate_use(:));
            end           
        end
    end
end

figure; hold on;
fr_bins = linspace(0,10,6);%logspace(-2,1,10);
fr_bin_centers = fr_bins(1:end-1)+diff(fr_bins)./2;
set(gca,'ColorOrder',cool(length(fr_bin_centers)));
for curr_fr_bin = 1:length(fr_bin_centers)
    celltype_act_corr_mean_fr_bin = cellfun(@(r,fr) ...
        nanmean(r(fr > fr_bins(curr_fr_bin) & ...
        fr <= fr_bins(curr_fr_bin+1))), ...
        celltype_act_corr,celltype_act_minrate);
    
    errorbar(nanmean(celltype_act_corr_mean_fr_bin,2)', ...
        AP_sem(celltype_act_corr_mean_fr_bin,2),'linewidth',1);
    drawnow;
end
errorbar(nanmean(celltype_act_corr_mean,2)', ...
    AP_sem(celltype_act_corr_mean,2),'k','linewidth',2);
xlim([0.5,3.5]);
set(gca,'XTick',1:3,'XTickLabel',{'MSN','FSI','TAN'});
ylabel('Pairwise activity correlation (across recordings)');
legend([cellfun(@(x) [num2str(x) 'sp/s'], ...
    num2cell(fr_bin_centers),'uni',false),'All rates']);



%% (quick test: fraction of response carried by fraction of cells)

curr_cells_idx = find(good_units_allcat & ismember(celltype_allcat,[1]) & domain_aligned_allcat == 1);

act_mean = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:,:),1)), ...
    mua_allcat_stimalign_exp,trial_stim_allcat_exp,move_t_exp, ...
    trial_outcome_allcat_exp,'uni',false)')';

use_t = t > 0 & t < 0.15;
sort_act = nanmean(act_mean(curr_cells_idx,use_t),2);
[~,sort_idx] = sort(sort_act,'ascend');

baseline_act = nanmean(act_mean(curr_cells_idx,t < 0),2);

a = act_mean(curr_cells_idx(sort_idx),:);
b = cumsum(a,1);

figure;
subplot(1,3,1);
imagesc(b);
subplot(1,3,2);
imagesc(b./max(b,[],1));
subplot(1,3,3); hold on;
plot((1:length(curr_cells_idx))./length(curr_cells_idx), ...
    nanmean(b(:,t < 0),2)./max(nanmean(b(:,t < 0),2)));
plot((1:length(curr_cells_idx))./length(curr_cells_idx), ...
    nanmean(b(:,use_t),2)./max(nanmean(b(:,use_t),2)));
line([0,1],[0,1],'color','k','linestyle','--');
xlabel('Fraction of cells');
ylabel('Fracion of total spikes');
legend({'Baseline','Stim'});






%% Rasters (all cells)

trial_stim_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
move_t_exp = mat2cell(move_t,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

% Sort cells by max activity time
mua_allcat_stimalign_exp = vertcat(mua_all{:});

act_mean = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:,:),1)), ...
    mua_allcat_stimalign_exp,trial_stim_allcat_exp,move_t_exp, ...
    trial_outcome_allcat_exp,'uni',false)')';

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
mua_allcat_exp_pad = cell2mat(permute(cellfun(@(act,trials_1,trials_2) ...
    padarray(act([find(trials_1);find(trials_2)],:,:), ...
    [max_trials-(sum(trials_1)+sum(trials_2)),0,0],NaN,'post'), ...
    mua_allcat_stimalign_exp,plot_trials_1,plot_trials_2,'uni',false),[2,3,1]));

cell_string = arrayfun(@(x) ...
    {[' Sort idx: ' num2str(x) ...
    ', Cell: ' num2str(sort_idx(x)) ...
    ', Domain: ' num2str(domain_aligned_allcat(sort_idx(x))) ...
    ', ' celltype_labels{celltype_allcat(sort_idx(x))}],...
    ['Animal: ' num2str(animals_allcat(sort_idx(x))) ...
    ', Day: ' num2str(days_allcat(sort_idx(x))) ...
    ', Neuron: ' num2str(neurons_allcat(sort_idx(x)))]}, ...
    1:length(sort_idx),'uni',false);

AP_imscroll(mua_allcat_exp_pad(:,:,sort_idx),cell_string);
line(repmat(find(t >= 0,1),2,1),ylim,'linestyle','--','color','r');
colormap(brewermap([],'Greys'));
caxis([0,50]);



%% Rasters (select cells) 

% Set cells and alignment to plot
% plot_cells = celltype_allcat == 1 & domain_aligned_allcat == 1;
plot_cells = good_units_allcat & ismember(celltype_allcat,[2]) & domain_aligned_allcat == 2;
% plot_cells = (stim_cells) & good_units_allcat & ismember(celltype_allcat,[1]) & domain_aligned_allcat == 1;
% plot_cells = ~stim_cells & ~move_cells & ~reward_cells & ismember(celltype_allcat,[1]) & domain_aligned_allcat == 1;
% plot_cells = stim_cells;
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

% (get experiment/neuron index)


plot_cells_idx = find(plot_cells);
cell_string = arrayfun(@(x) ...
    {[' Sort idx: ' num2str(x) ...
    ', Cell: ' num2str(plot_cells_idx(x)) ...
    ', Domain: ' num2str(domain_aligned_allcat(plot_cells_idx(x))) ...
    ', ' celltype_labels{celltype_allcat(plot_cells_idx(x))}],...
    ['Animal: ' num2str(animals_allcat(plot_cells_idx(x))) ...
    ', Day: ' num2str(days_allcat(plot_cells_idx(x))) ...
    ', Neuron: ' num2str(neurons_allcat(plot_cells_idx(x)))]}, ...
    1:length(plot_cells_idx),'uni',false);
AP_imscroll(act_exp_pad(:,:,plot_cells),cell_string);
line(repmat(find(t >= 0,1),2,1),ylim,'linestyle','--','color','r');
line(repmat(find(t >= 0,1),2,1)+length(t),ylim,'linestyle','--','color','r');
line(repmat(find(t >= 0,1),2,1)+length(t)*2,ylim,'linestyle','--','color','r');

colormap(brewermap([],'Greys'));
caxis([0,50]);


%% Sort and plot by different alignments

% Pick cells to plot
plot_cells = cellfun(@(good,celltype,domain) good & ismember(celltype,[4]), ...
    good_units_exp,celltype_exp,domain_exp,'uni',false);

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

% Get average cell activity across alignments
curr_str_act_stim_avg = cell2mat(cellfun(@(act,curr_trials,cells) ...
    permute(nanmean(act(curr_trials,:,cells),1),[3,2,1]), ...
    mua_allcat_stimalign_exp,plot_trials,plot_cells,'uni',false));
     
curr_str_act_move_avg = cell2mat(cellfun(@(act,curr_trials,cells) ...
    permute(nanmean(act(curr_trials,:,cells),1),[3,2,1]), ...
    mua_allcat_movealign_exp,plot_trials,plot_cells,'uni',false));

curr_str_act_outcome_avg = cell2mat(cellfun(@(act,curr_trials,cells) ...
    permute(nanmean(act(curr_trials,:,cells),1),[3,2,1]), ...
    mua_allcat_outcomealign_exp,plot_trials,plot_cells,'uni',false));

% Get sort index
% [~,sort_idx] = sort(depth_aligned_allcat(cell2mat(plot_cells)));

align_act = cell2mat(cellfun(@(act,curr_trials,cells) ...
    permute(nanmean(act(curr_trials,:,cells),1),[3,2,1]), ...
    mua_allcat_stimalign_exp,plot_trials,plot_cells,'uni',false));
[~,max_idx] = max([curr_str_act_stim_avg],[],2);
[~,sort_idx] = sort(max_idx);

% Plot
figure; hold on; set(gca,'YDir','reverse');

align_col = [1,0,0;0.8,0,0.8;0,0,0.8];
align_t = {[-0.05,0.15],[0.15,0.6],[0.6,1]};
for curr_align = 1:3
    
    switch curr_align
        case 1
            curr_act = curr_str_act_stim_avg;
            curr_shift = zeros(size(trial_stim_allcat));           
        case 2            
            curr_act = curr_str_act_move_avg;
            t_leeway = -t(1);
            leeway_samples = round(t_leeway*(sample_rate));
            curr_shift = -move_idx + leeway_samples;
        case 3            
            curr_act = curr_str_act_outcome_avg;
            t_leeway = -t(1);
            leeway_samples = round(t_leeway*(sample_rate));
            curr_shift = -outcome_idx + leeway_samples;
    end
    
    curr_t_offset = -nanmean(curr_shift(cell2mat(plot_trials)))/sample_rate;   
    curr_t = t + curr_t_offset;
    curr_t_plot = curr_t >= align_t{curr_align}(1) & ...
        curr_t <= align_t{curr_align}(2);
    
    curr_act_norm = curr_act./ ...
        max([curr_str_act_stim_avg,curr_str_act_move_avg,curr_str_act_outcome_avg],[],2);
    plot_t = curr_t > align_t{curr_align}(1) & curr_t < align_t{curr_align}(2);
    
    imagesc(curr_t(plot_t),[],curr_act_norm(sort_idx,plot_t));
    line(repmat(curr_t_offset,2,1),[0,length(sort_idx)],'color',align_col(curr_align,:));
    colormap(gca,brewermap([],'Greys'));
    
end

axis tight



%% Stim/move/reward sort and plot

% plot_cells = true(size(celltype_allcat));
plot_cells = good_units_allcat & ismember(celltype_allcat,[2]);

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

figure;imagesc(t,[],sort_act_smooth(sort_idx,:)./max(sort_act_smooth(sort_idx,:),[],2));
colormap(brewermap([],'*RdBu'));
caxis([-1,1]);
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

%% Trying to get task condition activity

% Get task activity for each cell
% (stim)
use_t = t > 0.05 & t < 0.15;

stim_act_cond = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    [squeeze(nanmean(nanmean(act(stim > 0 & rxn > 0.5,use_t,:),2),1)), ...
    squeeze(nanmean(nanmean(act(stim <= 0 & rxn > 0.5,use_t,:),2),1))], ...
    mua_allcat_stimalign_exp,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false));

vis_act = (stim_act_cond(:,1) - stim_act_cond(:,2));

% (move)
use_t = t > 0 & t < 0.5;

move_act_cond = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    [squeeze(nanmean(nanmean(act(stim <= 0 & rxn < 0.5,use_t,:),2),1)), ...
    squeeze(nanmean(nanmean(act(stim <= 0 & rxn >= 0.5,use_t,:),2),1))], ...
    mua_allcat_stimalign_exp_smooth,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false));

move_act = (move_act_cond(:,1) - move_act_cond(:,2));

% (reward)
use_t = t > 0 & t < 0.5;

reward_act_cond = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    [squeeze(nanmean(nanmean(act(stim > 0 & rxn < 0.5,use_t,:),2),1)), ...
    squeeze(nanmean(nanmean(act(stim <= 0 & rxn < 0.5,use_t,:),2),1))], ...
    mua_allcat_outcomealign_exp_smooth,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false));

reward_act = (reward_act_cond(:,1) - reward_act_cond(:,2))./ ...
    (sum(abs(reward_act_cond),2));





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

% (OR condition act)
mua_allcat_stimalign_exp_smooth = cellfun(@(x) AP_deconv_wf(x,true),mua_allcat_stimalign_exp,'uni',false);
mua_allcat_movealign_exp_smooth = cellfun(@(x) AP_deconv_wf(x,true),mua_allcat_movealign_exp,'uni',false);
mua_allcat_outcomealign_exp_smooth = cellfun(@(x) AP_deconv_wf(x,true),mua_allcat_outcomealign_exp,'uni',false);

use_t = t > 0 & t < 0.5;

stim_act_cond = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    [squeeze(max(nanmean(act(stim > 0 & rxn < 0.5,use_t,:),1),[],2)), ...
    squeeze(max(nanmean(act(stim <= 0 & rxn < 0.5,use_t,:),1),[],2))], ...
    mua_allcat_stimalign_exp_smooth,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false));

move_act_cond = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    [squeeze(max(nanmean(act(stim <= 0 & rxn < 0.5,use_t,:),1),[],2)), ...
    squeeze(max(nanmean(act(stim <= 0 & rxn >= 0.5,use_t,:),1),[],2))], ...
    mua_allcat_stimalign_exp_smooth,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false));

outcome_act_cond = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    [squeeze(max(nanmean(act(stim > 0 & rxn < 0.5,use_t,:),1),[],2)), ...
    squeeze(max(nanmean(act(stim <= 0 & rxn < 0.5,use_t,:),1),[],2))], ...
    mua_allcat_outcomealign_exp_smooth,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false));

max_act_diff = [diff(fliplr(stim_act_cond),[],2), ...
    diff(fliplr(move_act_cond),[],2), ...
    diff(fliplr(outcome_act_cond),[],2)];
[~,max_act_diff_idx] = max(max_act_diff,[],2);

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

% ATTEMPTS AT CRITERIA FOR PULLING OUT CELLS

% % Pull out cells
% use_t = t > 0 & t < 0.5;
% [~,stim_max_idx] = max(stim_contra_act,[],2);
% stim_cells = ismember(stim_max_idx,find(use_t)) & ...
%     (max(stim_contra_act(:,use_t),[],2) > ...
%     1.5*max(stim_ipsi_act(:,use_t),[],2)) & max_aligned_idx == 1;
% 
% use_t = t > 0 & t < 0.5;
% [~,move_max_idx] = max(move_early_act,[],2);
% move_cells = ismember(move_max_idx,find(use_t)) & ...
%     (max(move_early_act(:,use_t),[],2) > ...
%     1.5*max(move_late_act(:,use_t),[],2)) & max_aligned_idx == 2;
% 
% use_t = t > 0 & t < 0.5;
% [~,reward_max_idx] = max(reward_act,[],2);
% reward_cells = ismember(reward_max_idx,find(use_t)) & ...
%     (max(reward_act(:,use_t),[],2) > ...
%     1.5*max(no_reward_act(:,use_t),[],2)) & max_aligned_idx == 3;


% Pull out cells
% (BEST SO FAR)
use_t = t > 0 & t < 0.5;
[~,stim_max_idx] = max(stim_contra_act,[],2);
stim_cells = ismember(stim_max_idx,find(use_t)) & ...
    (max(stim_contra_act(:,use_t),[],2) > ...
    1.5*max(stim_ipsi_act(:,use_t),[],2)) & trial_corr_max_align == 1;% & ...
%     taskpred_partial_r2_allcat(:,1) > 0;

use_t = t > 0 & t < 0.5;
[~,move_max_idx] = max(move_early_act,[],2);
move_cells = ismember(move_max_idx,find(use_t)) & ...
    (max(move_early_act(:,use_t),[],2) > ...
    1.5*max(move_late_act(:,use_t),[],2)) & trial_corr_max_align == 2;% & ...
%     taskpred_partial_r2_allcat(:,2) > 0;

use_t = t > 0 & t < 0.5;
[~,reward_max_idx] = max(reward_act,[],2);
reward_cells = ismember(reward_max_idx,find(use_t)) & ...
    (max(reward_act(:,use_t),[],2) > ...
    1.5*max(no_reward_act(:,use_t),[],2)) & trial_corr_max_align == 3;% & ...
%     taskpred_partial_r2_allcat(:,4) > 0;


% % Pull out cells
% % (trying corr align OR max align?)
% use_t = t > 0 & t < 0.5;
% [~,stim_max_idx] = max(stim_contra_act,[],2);
% stim_cells = ismember(stim_max_idx,find(use_t)) & ...
%     (max(stim_contra_act(:,use_t),[],2) > ...
%     1.5*max(stim_ipsi_act(:,use_t),[],2)) & ...
%     (trial_corr_max_align == 1 | max_act_diff_idx == 1);
% 
% use_t = t > 0 & t < 0.5;
% [~,move_max_idx] = max(move_early_act,[],2);
% move_cells = ismember(move_max_idx,find(use_t)) & ...
%     (max(move_early_act(:,use_t),[],2) > ...
%     1.5*max(move_late_act(:,use_t),[],2)) & ...
%     (trial_corr_max_align == 2 | max_act_diff_idx == 2);
% 
% use_t = t > 0 & t < 0.5;
% [~,reward_max_idx] = max(reward_act,[],2);
% reward_cells = ismember(reward_max_idx,find(use_t)) & ...
%     (max(reward_act(:,use_t),[],2) > ...
%     1.5*max(no_reward_act(:,use_t),[],2)) & ...
%     (trial_corr_max_align == 3 | max_act_diff_idx == 3);

% % Pull out cells
% stim_cells = stim_act_cond(:,1) > 1.5.*stim_act_cond(:,2) & ...
%     max_act_diff_idx == 1;
% 
% move_cells = move_act_cond(:,1) > 1.5.*move_act_cond(:,2) & ...
%     max_act_diff_idx == 2;
% 
% reward_cells = outcome_act_cond(:,1) > 1.5.*outcome_act_cond(:,2) & ...
%     max_act_diff_idx == 3;


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

% (select trials)
trial_corr = cell2mat(cellfun(@(stim_act,move_act,outcome_act,stim,rxn,outcome) cell2mat(arrayfun(@(neuron) ...
    cat(2, ...
    nanmean(AP_itril(corr(stim_act(stim > 0 & outcome == 1,:,neuron)'),-1)), ...
    nanmean(AP_itril(corr(move_act(stim > 0 & outcome == 1,:,neuron)'),-1)), ...
    nanmean(AP_itril(corr(outcome_act(stim > 0 & outcome == 1,:,neuron)'),-1))), ...
    transpose(1:size(stim_act,3)),'uni',false)), ...
    mua_allcat_stimalign_exp_smooth,mua_allcat_movealign_exp_smooth,mua_allcat_outcomealign_exp_smooth, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false));

% % (all trials)
% trial_corr = cell2mat(cellfun(@(stim_act,move_act,outcome_act,stim,rxn,outcome) cell2mat(arrayfun(@(neuron) ...
%     cat(2, ...
%     nanmean(AP_itril(corr(stim_act(:,:,neuron)'),-1)), ...
%     nanmean(AP_itril(corr(move_act(:,:,neuron)'),-1)), ...
%     nanmean(AP_itril(corr(outcome_act(:,:,neuron)'),-1))), ...
%     transpose(1:size(stim_act,3)),'uni',false)), ...
%     mua_allcat_stimalign_exp_smooth,mua_allcat_movealign_exp_smooth,mua_allcat_outcomealign_exp_smooth, ...
%     trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false));

[trial_corr_max,trial_corr_max_align] = max(trial_corr,[],2);


% % (testing)
% plot_cells = cellfun(@(good,celltype,domain) good & ismember(celltype,[3,6]) & domain == 1, ...
%     good_units_exp,celltype_exp,domain_exp,'uni',false);





%% Correlation of cortical + unit activity?


% correlation with individual kernel
mua_exp = vertcat(mua_all{:});
mua_ctxpred_exp = vertcat(mua_ctxpred_all{:});

use_t = t < 0;
str_ctx_baseline_corr = cell2mat(cellfun(@(str,ctx) ...
    diag(corr(reshape(permute(str(:,use_t,:),[2,1,3]),[],size(str,3)), ...
    reshape(permute(ctx(:,use_t,:),[2,1,3]),[],size(ctx,3)),'type','spearman','rows','complete')),...
    mua_exp,mua_ctxpred_exp,'uni',false));

use_t = t > 0 & t < 0.5;
str_ctx_stim_corr = cell2mat(cellfun(@(str,ctx) ...
    diag(corr(reshape(permute(str(:,use_t,:),[2,1,3]),[],size(str,3)), ...
    reshape(permute(ctx(:,use_t,:),[2,1,3]),[],size(ctx,3)),'type','spearman','rows','complete')),...
    mua_exp,mua_ctxpred_exp,'uni',false));

% Plot correlation distribution
figure; 
p1 = distributionPlot( ...
    max(str_ctx_baseline_corr(good_units_allcat & celltype_allcat == 1,:),[],2), ...
    'groups',domain_aligned_allcat(good_units_allcat & celltype_allcat == 1), ...
    'histOri','left','xValues',[1:3]-0.25,'distWidth',0.5,'color','k');
p2 = distributionPlot( ...
    max(str_ctx_stim_corr(good_units_allcat & celltype_allcat == 1,:),[],2), ...
    'groups',domain_aligned_allcat(good_units_allcat & celltype_allcat == 1), ...
    'histOri','right','xValues',[1:3]+0.25,'distWidth',0.5,'color','r');
xlabel('Domain');
ylabel('Unit-cortex correlation');
legend([p1{1}(1),p2{1}(1)],{'Baseline','Stim period'});

linkaxes(get(gcf,'Children'),'xy');

% Plot cumulative correlation distribution
figure; 
curr_celltype = 1;
subplot(1,2,1);hold on
set(gca,'ColorOrder',copper(n_aligned_depths));
for curr_domain = 1:3
    plot_cells = ...
        good_units_allcat & domain_aligned_allcat == curr_domain & celltype_allcat == curr_celltype & ...
        ~isnan(str_ctx_baseline_corr) & ~isnan(str_ctx_stim_corr);
    plot(sort(str_ctx_baseline_corr(plot_cells)),[1:sum(plot_cells)]./sum(plot_cells),'linewidth',2);
end
title('Baseline');
subplot(1,2,2);hold on
set(gca,'ColorOrder',copper(n_aligned_depths));
for curr_domain = 1:3
    plot_cells = ...
        good_units_allcat & domain_aligned_allcat == curr_domain & celltype_allcat == curr_celltype & ...
        ~isnan(str_ctx_baseline_corr) & ~isnan(str_ctx_stim_corr);
    plot(sort(str_ctx_stim_corr(plot_cells)),[1:sum(plot_cells)]./sum(plot_cells),'linewidth',2);
end
title('Stim period');


% Get sum of activity in striatum
str_act_baseline_std = cell2mat(cellfun(@(x) std(reshape(x(:,t < 0,:),[],size(x,3)),[],1)', ...
    mua_exp,'uni',false));
str_act_baseline_rate = cell2mat(cellfun(@(x) mean(reshape(x(:,t < 0,:),[],size(x,3)),1)', ...
    mua_exp,'uni',false));
str_act_trial_rate = cell2mat(cellfun(@(x) mean(reshape(x(:,t > 0 & t < 0.5,:),[],size(x,3)),1)', ...
    mua_exp,'uni',false));


ctx_act_baseline_std = cell2mat(cellfun(@(x) std(reshape(x(:,t > 0 & t < 0.5,:),[],size(x,3)),[],1), ...
    fluor_kernelroi_deconv_exp,'uni',false));
ctx_act_trial_rate = cell2mat(cellfun(@(x) mean(reshape(x(:,t > 0 & t < 0.5,:),[],size(x,3)),1), ...
    fluor_kernelroi_deconv_exp,'uni',false));
ctx_std_change = ctx_act_trial_rate./ctx_act_baseline_std;

softnorm = 0.1;
str_rate_change = (str_act_trial_rate-str_act_baseline_rate)./(str_act_baseline_rate + softnorm);
str_std_change = (str_act_trial_rate)./(str_act_baseline_std + softnorm);


figure; hold on;
hist_edges = [linspace(0,50,100)];
histogram(str_rate_change(good_units_allcat & celltype_allcat == 1 & ...
    domain_aligned_allcat == 1),hist_edges,'normalization','probability')
histogram(str_rate_change(good_units_allcat & celltype_allcat == 1 & ...
    domain_aligned_allcat == 2),hist_edges,'normalization','probability')
histogram(str_rate_change(good_units_allcat & celltype_allcat == 1 & ...
    domain_aligned_allcat == 3),hist_edges,'normalization','probability')



% Plot std change of ctx vs str?
figure; hold on;
hist_edges = [linspace(0,2,20)];
distributionPlot(ctx_std_change,'histOpt',0,'divFactor',hist_edges, ...
    'histOri','left','xValues',[1:3]-0.25,'distWidth',0.5,'color',[0,0.7,0])
distributionPlot(str_std_change(good_units_allcat & celltype_allcat == 1), ...
    'groups',domain_aligned_allcat(good_units_allcat & celltype_allcat == 1), ...
    'histOpt',0,'divFactor',hist_edges,'histOri','right','xValues',[1:3]+0.25,'distWidth',0.5,'color','k')
ylabel('Std change');
xlabel('Domain');



%% ~~~~~~~~~~~~ Maps

%% QUICK TEST: map tsne

downsamp_factor = 1/10;
ctx_str_k_px_downsamp_reshape = reshape(imresize(ctx_str_k_px, ...
    downsamp_factor,'bilinear'),[],size(ctx_str_k_px,3));

use_cells = domain_aligned_allcat == 1 & celltype_allcat == 1;
y = tsne(ctx_str_k_px_downsamp_reshape(:,use_cells)','Distance','Correlation');
figure;gscatter(y(:,1),y(:,2));




%% Get relative map weight in ROI

ctx_str_k_px_maxnorm = ctx_str_k_px./nanmax(nanmax(ctx_str_k_px,[],1),[],2);

AP_imscroll(ctx_str_k_px_maxnorm);
axis image;
caxis([-1,1])
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

% (get ROI and sort)
[~,sort_idx] = sort(roi.trace,'descend');
AP_imscroll(ctx_str_k_px_maxnorm(:,:,sort_idx));
axis image;
caxis([-1,1])
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);


%% Plot average maps by domain/celltype

% Plot average within each domain AND cell type
figure;
for curr_depth = 1:n_aligned_depths
    for curr_celltype = 1:n_celltypes
        subplot(n_aligned_depths,n_celltypes, ...
            (curr_depth-1)*n_celltypes+curr_celltype)
        imagesc(nanmean(ctx_str_k_px(:,:, ...
            domain_aligned_allcat(good_units_allcat) == curr_depth & ...
            celltype_allcat(good_units_allcat) == curr_celltype),3));
        axis image
        colormap(brewermap([],'*RdBu'));
        caxis([-max(abs(caxis)),max(abs(caxis))]);
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        axis off
        title(['Str ' num2str(curr_depth) ' ' celltype_labels{curr_celltype}]);
    end
end



%%%% TESTING





% Correlate maps
curr_depth = 1;

msn_maps = ctx_str_k_px(:,:, ...
    domain_aligned_allcat(good_units_allcat) == curr_depth & ...
    celltype_allcat(good_units_allcat) == 1);

fsi_maps = ctx_str_k_px(:,:, ...
    domain_aligned_allcat(good_units_allcat) == curr_depth & ...
    celltype_allcat(good_units_allcat) == 2);

tan_maps = ctx_str_k_px(:,:, ...
    domain_aligned_allcat(good_units_allcat) == curr_depth & ...
    ismember(celltype_allcat(good_units_allcat),[3,6]));

% (cell-cell map corr?)
map_corr = mat2cell(corr([reshape(msn_maps,[],size(msn_maps,3)), ...
    reshape(fsi_maps,[],size(fsi_maps,3)), ...
    reshape(tan_maps,[],size(tan_maps,3))]), ...
    [size(msn_maps,3),size(fsi_maps,3),size(tan_maps,3)], ...
    [size(msn_maps,3),size(fsi_maps,3),size(tan_maps,3)]);

nanmean(AP_itril(map_corr{1,1},-1))

across_corr = nanmean(map_corr{1,2}(:));

map_corr_cat = cell2mat(map_corr);
n_shuff = 100;
across_corr_shuff = nan(n_shuff,1);
for curr_shuff = 1:n_shuff
    curr_perm = randperm(length(map_corr_cat));
    curr_map_corr_cat_shuff = map_corr_cat(curr_perm,curr_perm);
    
    curr_map_corr_shuff = mat2cell(curr_map_corr_cat_shuff, ...
        [size(msn_maps,3),size(fsi_maps,3),size(tan_maps,3)], ...
        [size(msn_maps,3),size(fsi_maps,3),size(tan_maps,3)]);
    
    across_corr_shuff(curr_shuff) = nanmean(curr_map_corr_shuff{1,2}(:));
    AP_print_progress_fraction(curr_shuff,n_shuff);    
end

corr_rank = tiedrank([across_corr;across_corr_shuff]);
corr_p = (corr_rank(1)/(n_shuff+1));



% (cell-template map corr?)


% Load kernel templates
n_aligned_depths = 3;
kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template_' num2str(n_aligned_depths) '_depths.mat'];
load(kernel_template_fn);

figure; hold on;
str_col = max(hsv(n_aligned_depths)-0.2,0);
set(gca,'ColorOrder',str_col);
for curr_depth = 1:n_aligned_depths
    
    use_cells = good_units_allcat & domain_aligned_allcat == curr_depth & ...
        ismember(celltype_allcat,[1,2,3]);
    use_maps = ctx_str_k_px(:,:, ...
        domain_aligned_allcat(good_units_allcat) == curr_depth & ...
        ismember(celltype_allcat(good_units_allcat),[1,2,3]));
    
    map_corr = corr(reshape(use_maps,[],size(use_maps,3)), ...
        reshape(kernel_template(:,:,curr_depth),[],1),'type','Pearson');
    
    curr_p = anovan(map_corr,[recordings_allcat(use_cells),celltype_allcat(use_cells)],'display','off');
    disp(curr_p(2))
    a = accumarray([recordings_allcat(use_cells),celltype_allcat(use_cells)],map_corr,[],@nanmean);
    
    errorbar(nanmean(a,1),AP_sem(a,1),'linewidth',2)
    
end
set(gca,'XTick',1:n_aligned_depths,'XTickLabel',{'MSN','FSI','TAN'});
ylabel('Cell-domain kernel correlation')

% Plot kernel template
figure
for curr_depth = 1:n_aligned_depths
    subplot(n_aligned_depths,1,curr_depth);
    imagesc(kernel_template(:,:,curr_depth))
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'*RdBu'));
    axis image off;
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title('Domain kernel')
end




% Plot average within each celltype/domain
figure;
for curr_depth = 1:n_aligned_depths
    for curr_celltype = 1:n_celltypes
        subplot(n_aligned_depths,n_celltypes, ...
            (curr_depth-1)*n_celltypes+curr_celltype)
        imagesc(nanmean(ctx_str_k_px(:,:, ...
            domain_aligned_allcat(good_units_allcat) == curr_depth & ...
            celltype_allcat(good_units_allcat) == curr_celltype),3));
        axis image
        colormap(brewermap([],'*RdBu'));
        caxis([-max(abs(caxis)),max(abs(caxis))]);
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        axis off
        title(['Str ' num2str(curr_depth) ' ' celltype_labels{curr_celltype}]);
    end
end

% Correlate maps for each celltype within and across domains

% Load kernel templates
n_aligned_depths = 3;
kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template_' num2str(n_aligned_depths) '_depths.mat'];
load(kernel_template_fn);

for curr_depth = 1:n_aligned_depths
    for curr_celltype = 1:3
        
        use_cells = good_units_allcat & domain_aligned_allcat == curr_depth & ...
            ismember(celltype_allcat,[1,2,3]);
        use_maps = ctx_str_k_px(:,:, ...
            domain_aligned_allcat(good_units_allcat) == curr_depth & ...
            ismember(celltype_allcat(good_units_allcat),[1,2,3]));
        
    end
end

% Get map correlations for each celltype within/across domain
% (this is a super dirty way to do it but I'm in a hurry)
n_recordings = max(recordings_allcat);
celltype_map_corr = cell(3,2,n_recordings);
for curr_depth = 1:n_aligned_depths
    for curr_celltype = 1:3
        for curr_recording = 1:n_recordings
            
            use_maps = ctx_str_k_px(:,:, ...
                domain_aligned_allcat(good_units_allcat) == curr_depth & ...
                celltype_allcat(good_units_allcat) == curr_celltype & ...
                recordings_allcat(good_units_allcat) == curr_recording);
            
            if ~isempty(use_maps)
                map_corr = corr(reshape(use_maps,[],size(use_maps,3)), ...
                reshape(kernel_template,[],n_aligned_depths),'type','Pearson');
            else 
               map_corr = nan('single');
            end
    
            [cell_id,depth_group_id] = ndgrid(1:size(map_corr,1),1:size(map_corr,2));
            depth_group = (depth_group_id ~= curr_depth)+1;
            
            for curr_depth_compare = 1:2
                celltype_map_corr{curr_celltype,curr_depth_compare,curr_recording} = ...
                   cat(1,celltype_map_corr{curr_celltype,curr_depth_compare,curr_recording}, ...
                   reshape(map_corr(depth_group == curr_depth_compare),[],1));       
            end
        end
    end
end

celltype_map_corr_mean = cellfun(@nanmean,celltype_map_corr);

figure; hold on;
errorbar(nanmean(celltype_map_corr_mean,3)', ...
    AP_sem(celltype_map_corr_mean,3)','linewidth',2);
xlim([0.5,2.5]);
set(gca,'XTick',1:2,'XTickLabel',{'Within domain','Across domain'});
ylabel('Cell-domain kernel correlation')
legend({'MSN','FSI','TAN'});

% (Map correlation vs. celltype statistics)
disp('Cell-domain map correlation by celltype (2-way anova):');
[celltype_grp,domain_grp,exp_grp] = ndgrid(1:size(celltype_map_corr_mean,1), ...
    1:size(celltype_map_corr_mean,2),1:size(celltype_map_corr_mean,3));
p = anovan(celltype_map_corr_mean(:),...
    [celltype_grp(:),domain_grp(:)],'model','interaction','display','off');
disp(['p(cell type) = ' num2str(p(1))]);
disp(['p(cell type) = ' num2str(p(2))]);
disp(['p(cell type) = ' num2str(p(3))]);


%% Plot maps by condition 1/2

curr_domain = 1;

% max-normalize maps?
ctx_str_k_px_norm = ctx_str_k_px./max(max(ctx_str_k_px,[],1),[],2);

figure;

for curr_celltype = 1:n_celltypes
    
    curr_k_1 = nanmean(ctx_str_k_px_norm(:,:, ...
        stim_cells & ...
        domain_aligned_allcat == curr_domain & ...
        celltype_allcat == curr_celltype),3);
    
     curr_k_2 = nanmean(ctx_str_k_px_norm(:,:, ...
        ~stim_cells & ~move_cells & ~reward_cells & ...
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



%% STA maps from different time points

% Pick time to use for regression
% use_t = t < 0;
% use_t = t > 0 & t < 0.5;
% use_t = t < 0 | t > 0.5; 
use_t = true(size(t));

% Regress cortex to striatum
mua_ctx_sta = cell(size(mua_allcat_stimalign_exp));
for curr_exp = 1:length(trials_recording)
        
        n_cells = size(mua_allcat_stimalign_exp{curr_exp},3);
        n_vs = size(fluor_deconv_allcat_exp{curr_exp},3);
    
        curr_mua = reshape(permute(mua_allcat_stimalign_exp{curr_exp}(:,use_t,:), ...
            [2,1,3]),[],n_cells)';
        
        curr_fluor = fluor_deconv_allcat_exp{curr_exp};
        curr_fluor_baselinesub = curr_fluor - nanmean(curr_fluor(:,t < 0,:),2);
        curr_fluor_flat = reshape(permute(curr_fluor_baselinesub(:,use_t,:),[2,1,3]),[],n_vs)';

        %%% Get spike-triggered average      
        mua_ctx_sta{curr_exp} = (curr_fluor_flat*curr_mua')./sum(curr_mua',1);
     
        AP_print_progress_fraction(curr_exp,length(trials_recording));
end

% Get pixels
mua_ctx_sta_cat = horzcat(mua_ctx_sta{:});
mua_ctx_sta_px = svdFrameReconstruct(U_master(:,:,1:n_vs),mua_ctx_sta_cat(:,good_units_allcat));

% Plot STA top, kernel bottom
norm_sta_kernel = [(mua_ctx_sta_px./max(max(mua_ctx_sta_px,[],1),[],2)); ...
    (ctx_str_k_px./max(max(ctx_str_k_px,[],1),[],2))];
AP_imscroll(norm_sta_kernel);
axis image off
colormap(brewermap([],'*RdBu'));
caxis([-1,1])
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);



%% Regress maps from different time points

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
% use_t = t > 0 & t < 0.5;
% use_t = t > 1;
% use_t = t < 0 | t > 0.5;
% use_t = true(size(t));

mua_exp = vertcat(mua_all{:});

% Regress cortex to striatum
mua_ctxtrialpred_exp = cellfun(@(x) nan(size(x)),mua_exp,'uni',false);
mua_ctxtrialpred_k = cell(size(mua_ctxtrialpred_exp));
for curr_exp = 1:length(trials_recording)
        
        n_cells = size(mua_exp{curr_exp},3);
        n_vs = size(fluor_deconv_allcat_exp{curr_exp},3);
    
        curr_mua = reshape(permute(mua_exp{curr_exp}(:,use_t,:), ...
            [2,1,3]),[],n_cells)';
        curr_mua_std = curr_mua./nanstd(curr_mua,[],2);
        
        curr_fluor = reshape(permute( ...
            fluor_deconv_allcat_exp{curr_exp}(:,use_t,:),[2,1,3]),[],n_vs)';
        curr_fluor_full = reshape(permute( ...
            fluor_deconv_allcat_exp{curr_exp}(:,:,:),[2,1,3]),[],n_vs)';

        % Set discontinuities in trial data
        trial_discontinuities = false(size(mua_exp{curr_exp}(:,use_t,1)));
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

curr_domain = 1;

figure;

for curr_celltype = 1:n_celltypes
    
    curr_k_1 = nanmean(ctx_str_trialk_px(:,:, ...
        stim_cells & ...
        domain_aligned_allcat == curr_domain & ...
        celltype_allcat == curr_celltype),3);
    
     curr_k_2 = nanmean(ctx_str_trialk_px(:,:, ...
        ~stim_cells& ~move_cells & ~reward_cells & ...
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

%% ~~~~~~~~~~~~ Activity vs. maps

%% Correlate activity/map with average activity/map (all together)

avg_act = AP_deconv_wf(cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:,:),1)),mua_allcat_stimalign_exp, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false)')',true);

curr_celltype = celltype_allcat == 1;
curr_cells = good_units_allcat & domain_aligned_allcat == 2 & curr_celltype;

avg_act_domain = nanmean(zscore(avg_act(curr_cells,:),[],2),1);
avg_map_domain = nanmean(ctx_str_k_px(:,:,curr_cells(good_units_allcat)),3);

act_corr = 1-pdist2(avg_act_domain,avg_act,'correlation');
map_corr = 1-pdist2(avg_map_domain(:)',reshape(ctx_str_k_px,[],size(ctx_str_k_px,3))','correlation');

figure; 
subplot(4,3,1);
plot(t,avg_act_domain,'k','linewidth',2);
title('Average activity');

subplot(4,3,2);
imagesc(avg_map_domain);
axis image off
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('Average map');

subplot(4,3,[4,7,10]);
plot(act_corr(curr_celltype),depth_aligned_allcat(curr_celltype),'.k')
set(gca,'YDir','reverse');
xlabel('Activity correlation');
ylabel('Depth');

subplot(4,3,[5,8,11]);
plot(map_corr(curr_celltype(good_units_allcat)),depth_aligned_allcat(good_units_allcat & curr_celltype),'.k')
set(gca,'YDir','reverse');
xlabel('Map correlation');
ylabel('Depth');

subplot(4,3,3);
plot(act_corr(good_units_allcat & curr_celltype),map_corr(curr_celltype(good_units_allcat)),'.k')
xlabel('Activity correlation');
ylabel('Map correlation');
line([-1,1],[-1,1])

%% Correlate activity/map with average activity/map (within domain)

figure;

for curr_domain = 1:n_aligned_depths
    
    avg_act = AP_deconv_wf(cell2mat(cellfun(@(act,stim,rxn,outcome) ...
        squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:,:),1)),mua_allcat_stimalign_exp, ...
        trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false)')',true);
    
    curr_celltype = celltype_allcat == 1;
    curr_cells = good_units_allcat & domain_aligned_allcat == 1 & curr_celltype;
    
    avg_act_domain = nanmean(zscore(avg_act(curr_cells,:),[],2),1);
    avg_map_domain = nanmean(ctx_str_k_px(:,:,curr_cells(good_units_allcat)),3);
    
    act_corr = 1-pdist2(avg_act_domain,avg_act,'correlation');
    map_corr = 1-pdist2(avg_map_domain(:)',reshape(ctx_str_k_px,[],size(ctx_str_k_px,3))','correlation');
    
    hist_edges = linspace(-1,1,20);
    
    act_corr_hist = histcounts(act_corr(curr_cells),hist_edges);
    map_corr_hist = histcounts(map_corr(curr_cells(good_units_allcat)),hist_edges);
    
    
    
    
    
    subplot(n_aligned_depths,4,1);
    plot(t,avg_act_domain,'k','linewidth',2);
    title('Average activity');
    
    subplot(n_aligned_depths,4,2);
    imagesc(avg_map_domain);
    axis image off
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'*RdBu'));
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title('Average map');
    
    subplot(n_aligned_depths,4,3);
    plot(act_corr(curr_cells),map_corr(curr_cells(good_units_allcat)),'.k')
    xlabel('Activity correlation');
    ylabel('Map correlation');
    line([-1,1],[-1,1]);
    line(xlim,[-1,1],'linestyle','--');
    line([-1,1],ylim,'linestyle','--');
    
    subplot(n_aligned_depths,4,4); hold on;
    histogram(act_corr(curr_cells),hist_edges,'normalization','probability');
    histogram(map_corr(curr_cells(good_units_allcat)),hist_edges,'normalization','probability');
    xlabel('Correlation');
    ylabel('Fraction');
    legend({'Activity','Map'});
    
end


%% TESTING from above

plot_cells = domain_aligned_allcat == 1 & celltype_allcat == 1;
plot_cells_idx = find(plot_cells);

[~,sort_idx] = sort(map_corr(plot_cells));
AP_imscroll(ctx_str_k_px(:,:,plot_cells_idx(sort_idx)));
axis image off
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);



%% Correlate activity/map with average activity/map (by experiment)

% Get cortical activity
fluor_allcat = cell2mat(vertcat(fluor_all{:}));
fluor_allcat_deconv = AP_deconv_wf(fluor_allcat);
% fluor_allcat_deconv_baseline = nanmean(reshape(fluor_allcat_deconv(:,t_baseline,:),[],1,n_vs));
% fluor_allcat_deconv = fluor_allcat_deconv - fluor_allcat_deconv_baseline;

kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
load(kernel_roi_fn);
n_kernel_rois = size(kernel_roi.bw,3);
fluor_kernelroi_deconv = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],kernel_roi.bw), ...
    size(kernel_roi.bw,3),[],size(fluor_allcat_deconv,1)),[3,2,1]);

% (stim aligned)
fluor_kernelroi_deconv_stimalign_exp = mat2cell(fluor_kernelroi_deconv,use_split,length(t),n_kernel_rois);

% (move aligned)
move_idx_exp = mat2cell(move_idx,use_split,1);
fluor_kernelroi_deconv_movealign_exp = fluor_kernelroi_deconv_stimalign_exp;
for curr_exp = 1:length(fluor_kernelroi_deconv_movealign_exp)
   for curr_trial = 1:size(fluor_kernelroi_deconv_movealign_exp{curr_exp},1)
       fluor_kernelroi_deconv_movealign_exp{curr_exp}(curr_trial,:,:) = ...
           circshift(fluor_kernelroi_deconv_movealign_exp{curr_exp}(curr_trial,:,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
   end
end

% (outcome aligned)
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
fluor_kernelroi_deconv_outcomealign_exp = fluor_kernelroi_deconv_stimalign_exp;
for curr_exp = 1:length(fluor_kernelroi_deconv_outcomealign_exp)
    for curr_trial = 1:size(fluor_kernelroi_deconv_outcomealign_exp{curr_exp},1)
        fluor_kernelroi_deconv_outcomealign_exp{curr_exp}(curr_trial,:,:) = ...
            circshift(fluor_kernelroi_deconv_outcomealign_exp{curr_exp}(curr_trial,:,:), ...
            -outcome_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
    end
end

% Split striatal maps / indicies by experiment
n_cells_exp = cellfun(@(x) length(x),good_units_exp);
n_good_cells_exp = cellfun(@(x) sum(x),vertcat(goodUnits{:}));
ctx_str_k_px_exp = reshape(mat2cell(ctx_str_k_px, ...
    size(U_master,1),size(U_master,2),n_good_cells_exp),[],1);

domain_exp = mat2cell(domain_aligned_allcat,n_cells_exp,1);
celltype_exp = mat2cell(celltype_allcat,n_cells_exp,1);

% Loop through domains and get activity/map correlations
curr_celltype = 1;

avg_map_exp = cell(length(n_cells_exp),n_aligned_depths);
avg_act_exp = cell(length(n_cells_exp),n_aligned_depths);
avg_ctxact_exp = cell(length(n_cells_exp),n_aligned_depths);

avg_map_corr_exp = cell(length(n_cells_exp),n_aligned_depths);
avg_act_corr_exp = cell(length(n_cells_exp),n_aligned_depths);
avg_ctxact_corr_exp = cell(length(n_cells_exp),n_aligned_depths);

for curr_domain = 1:n_aligned_depths
    
    switch curr_domain
        case 1
            str_act_align = mua_allcat_stimalign_exp;
            ctx_act_align = fluor_kernelroi_deconv_stimalign_exp;
        case 2
            str_act_align = mua_allcat_movealign_exp;
            ctx_act_align = fluor_kernelroi_deconv_movealign_exp;
        case 3
            str_act_align = mua_allcat_outcomealign_exp;
            ctx_act_align = fluor_kernelroi_deconv_outcomealign_exp;
    end
    
    avg_map_exp(:,curr_domain) = cellfun(@(map,domain,celltype,good_units) ...
        nanmean(map(:,:,domain(good_units) == curr_domain & ...
        celltype(good_units) == curr_celltype),3), ...
        ctx_str_k_px_exp,domain_exp,celltype_exp,good_units_exp,'uni',false);
    
    avg_act_exp(:,curr_domain) = cellfun(@(act,stim,rxn,outcome,domain,celltype,good_units) ...
        permute(nanmean(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:, ...
        good_units & domain == curr_domain & celltype == curr_celltype),1),3),[3,2,1]), ...
        str_act_align, ...
        trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp, ...
        domain_exp,celltype_exp,good_units_exp,'uni',false);
    
    avg_ctxact_exp(:,curr_domain) = cellfun(@(act,stim,rxn,outcome,domain,celltype,good_units) ...
        permute(nanmean(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:,curr_domain),1),3),[3,2,1]), ...
        ctx_act_align, ...
        trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp, ...
        domain_exp,celltype_exp,good_units_exp,'uni',false);
    
    % Map/avg map correlation
    avg_map_corr_exp(:,curr_domain) = ...
        cellfun(@(avg_map,map,domain,celltype,good_units) 1-pdist2( ...
        reshape(map(:,:,domain(good_units) == curr_domain & ...
        celltype(good_units) == curr_celltype),size(map,1)*size(map,2),[])', ...
        avg_map(:)','correlation'),avg_map_exp(:,curr_domain),ctx_str_k_px_exp, ...
        domain_exp,celltype_exp,good_units_exp,'uni',false);
    
    % Act/avg act correlation
    avg_act_corr_exp(:,curr_domain) = ...
        cellfun(@(avg_act,act,stim,rxn,outcome,domain,celltype,good_units) ...
        1-pdist2(permute(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:, ...
        good_units & domain == curr_domain & celltype == curr_celltype),1),[3,2,1]), ...
        avg_act,'correlation'),avg_act_exp(:,curr_domain),str_act_align, ...
        trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp, ...
        domain_exp,celltype_exp,good_units_exp,'uni',false);
    
    % Act/avg ctx act correlation
    avg_ctxact_corr_exp(:,curr_domain) = ...
        cellfun(@(avg_act,act,stim,rxn,outcome,domain,celltype,good_units) ...
        1-pdist2(permute(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:, ...
        good_units & domain == curr_domain & celltype == curr_celltype),1),[3,2,1]), ...
        avg_act,'correlation'),avg_ctxact_exp(:,curr_domain),str_act_align, ...
        trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp, ...
        domain_exp,celltype_exp,good_units_exp,'uni',false);
    
end

figure;
for curr_domain = 1:n_aligned_depths
    
    % Plot average activity
    % (normalize activity: subtract baseline, divide std)
    % (!! assume first domain is stim-aligned !!)

    curr_act = cell2mat(cellfun(@(act,act_stimalign,baseline) ...
        (act-nanmean(act_stimalign(:,t < 0,:),2))./nanstd(act,[],2), ...
        avg_act_exp(:,curr_domain),avg_act_exp(:,1),'uni',false));
    curr_ctxact = cell2mat(cellfun(@(act,act_stimalign,baseline) ...
        (act-nanmean(act_stimalign(:,t < 0,:),2))./nanstd(act,[],2), ...
        avg_ctxact_exp(:,curr_domain),avg_ctxact_exp(:,1),'uni',false));
 
    subplot(n_aligned_depths,4,(curr_domain-1)*4+1); hold on;
    AP_errorfill(t,nanmean(curr_act,1)',AP_sem(curr_act,1)','k');
    AP_errorfill(t,nanmean(curr_ctxact,1)',AP_sem(curr_ctxact,1)',[0,0.7,0]);
    line([0,0],ylim,'color','r');
    ylabel('Activity');
    title(['Domain ' num2str(curr_domain)]);
    
    % Plot average map
    subplot(n_aligned_depths,4,(curr_domain-1)*4+2);
    imagesc(nanmean(cat(3,avg_map_exp{:,curr_domain}),3));
    axis image off
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'*RdBu'));
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    
    % Plot map/act correlation points
    subplot(n_aligned_depths,4,(curr_domain-1)*4+3); hold on;
    plot(cell2mat(avg_map_corr_exp(:,curr_domain)), ...
        cell2mat(avg_act_corr_exp(:,curr_domain)),'.k');
    plot(cell2mat(avg_map_corr_exp(:,curr_domain)), ...
        cell2mat(avg_ctxact_corr_exp(:,curr_domain)),'.','color',[0,0.7,0]);
    line([0,0],[-1,1],'linestyle','--');
    line([-1,1],[0,0],'linestyle','--');
    xlabel('Map correlation');
    ylabel('Activity correlation');
        
    % Plot correlation histogram
    hist_edges = linspace(-1,1,20);    
    subplot(n_aligned_depths,4,(curr_domain-1)*4+4); hold on;
    histogram(cell2mat(avg_act_corr_exp(:,curr_domain)),hist_edges,'normalization','probability');
    histogram(cell2mat(avg_ctxact_corr_exp(:,curr_domain)),hist_edges,'normalization','probability');
    histogram(cell2mat(avg_map_corr_exp(:,curr_domain)),hist_edges,'normalization','probability');
    xlabel('Correlation');
    ylabel('Fraction');
    legend({'Avg activity','Avg ctx activity','Maps'});
    
end

% TESTING 
figure; 
for curr_domain = 1:n_aligned_depths

    % Plot map/act correlation points
    subplot(n_aligned_depths,2,(curr_domain-1)*2+1); hold on;
    plot(cell2mat(avg_act_corr_exp(:,curr_domain)), ...
        cell2mat(avg_ctxact_corr_exp(:,curr_domain)),'.k');
    line([0,0],[-1,1],'linestyle','--');
    line([-1,1],[0,0],'linestyle','--');
    xlabel('Str activity correlation');
    ylabel('Ctx activity correlation');
    
    subplot(n_aligned_depths,2,(curr_domain-1)*2+2); hold on;
    plot([cell2mat(avg_act_corr_exp(:,curr_domain)), ...
        cell2mat(avg_ctxact_corr_exp(:,curr_domain))]','color',[0.5,0.5,0.5]);
    errorbar(nanmean([cell2mat(avg_act_corr_exp(:,curr_domain)), ...
        cell2mat(avg_ctxact_corr_exp(:,curr_domain))],1), ...
        AP_sem([cell2mat(avg_act_corr_exp(:,curr_domain)), ...
        cell2mat(avg_ctxact_corr_exp(:,curr_domain))],1),'color','k','linewidth',2);
    set(gca,'XTick',1:2,'XTickLabel',{'Str act','Ctx act'});
    xlim([0.5,2.5])
    ylabel('Correlation');
    
end


%% Correlate activity/map with average activity/map (across experiments)

% Get cortical activity
fluor_allcat = cell2mat(vertcat(fluor_all{:}));
fluor_allcat_deconv = AP_deconv_wf(fluor_allcat);
% fluor_allcat_deconv_baseline = nanmean(reshape(fluor_allcat_deconv(:,t_baseline,:),[],1,n_vs));
% fluor_allcat_deconv = fluor_allcat_deconv - fluor_allcat_deconv_baseline;

kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
load(kernel_roi_fn);
n_kernel_rois = size(kernel_roi.bw,3);
fluor_kernelroi_deconv = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],kernel_roi.bw), ...
    size(kernel_roi.bw,3),[],size(fluor_allcat_deconv,1)),[3,2,1]);

% (stim aligned)
fluor_kernelroi_deconv_stimalign_exp = mat2cell(fluor_kernelroi_deconv,use_split,length(t),n_kernel_rois);

% (move aligned)
move_idx_exp = mat2cell(move_idx,use_split,1);
fluor_kernelroi_deconv_movealign_exp = fluor_kernelroi_deconv_stimalign_exp;
for curr_exp = 1:length(fluor_kernelroi_deconv_movealign_exp)
   for curr_trial = 1:size(fluor_kernelroi_deconv_movealign_exp{curr_exp},1)
       fluor_kernelroi_deconv_movealign_exp{curr_exp}(curr_trial,:,:) = ...
           circshift(fluor_kernelroi_deconv_movealign_exp{curr_exp}(curr_trial,:,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
   end
end

% (outcome aligned)
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
fluor_kernelroi_deconv_outcomealign_exp = fluor_kernelroi_deconv_stimalign_exp;
for curr_exp = 1:length(fluor_kernelroi_deconv_outcomealign_exp)
    for curr_trial = 1:size(fluor_kernelroi_deconv_outcomealign_exp{curr_exp},1)
        fluor_kernelroi_deconv_outcomealign_exp{curr_exp}(curr_trial,:,:) = ...
            circshift(fluor_kernelroi_deconv_outcomealign_exp{curr_exp}(curr_trial,:,:), ...
            -outcome_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
    end
end

% Split striatal maps / indicies by experiment
n_cells_exp = cellfun(@(x) length(x),good_units_exp);
n_good_cells_exp = cellfun(@(x) sum(x),vertcat(goodUnits{:}));
ctx_str_k_px_exp = reshape(mat2cell(ctx_str_k_px, ...
    size(U_master,1),size(U_master,2),n_good_cells_exp),[],1);

domain_exp = mat2cell(domain_aligned_allcat,n_cells_exp,1);
celltype_exp = mat2cell(celltype_allcat,n_cells_exp,1);

% Loop through domains and get activity/map correlations
curr_celltype = 3;

avg_map = cell(1,n_aligned_depths);
avg_act = cell(1,n_aligned_depths);
avg_ctxact = cell(1,n_aligned_depths);

avg_map_corr_exp = cell(length(n_cells_exp),n_aligned_depths);
avg_act_corr_exp = cell(length(n_cells_exp),n_aligned_depths);
avg_ctxact_corr_exp = cell(length(n_cells_exp),n_aligned_depths);

for curr_domain = 1:n_aligned_depths
    
    switch curr_domain
        case 1
            str_act_align = mua_allcat_stimalign_exp;
            ctx_act_align = fluor_kernelroi_deconv_stimalign_exp;
        case 2
            str_act_align = mua_allcat_movealign_exp;
            ctx_act_align = fluor_kernelroi_deconv_movealign_exp;
        case 3
            str_act_align = mua_allcat_outcomealign_exp;
            ctx_act_align = fluor_kernelroi_deconv_outcomealign_exp;
    end
    
    curr_units = good_units_allcat & domain_aligned_allcat == curr_domain & ...
        celltype_allcat == curr_celltype;
    curr_good_units = domain_aligned_allcat(good_units_allcat) == curr_domain & ...
        celltype_allcat(good_units_allcat) == curr_celltype;
    
    avg_map{curr_domain} = nanmean(ctx_str_k_px(:,:,curr_good_units),3);
    avg_act{curr_domain} = nanmean(cell2mat(cellfun(@(act,stim,rxn,outcome,domain,celltype,good_units) ...
        permute(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:, ...
        good_units & domain == curr_domain & celltype == curr_celltype),1),[3,2,1]), ...
        str_act_align,trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp, ...
        domain_exp,celltype_exp,good_units_exp,'uni',false)),1);
    avg_ctxact{curr_domain} = nanmean(cell2mat(cellfun(@(act,stim,rxn,outcome,domain,celltype,good_units) ...
        permute(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:,curr_domain),1),[3,2,1]), ...
        ctx_act_align,trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp, ...
        domain_exp,celltype_exp,good_units_exp,'uni',false)),1);
    
    % Map/avg map correlation
    avg_map_corr_exp(:,curr_domain) = ...
        cellfun(@(map,domain,celltype,good_units) 1-pdist2( ...
        reshape(map(:,:,domain(good_units) == curr_domain & ...
        celltype(good_units) == curr_celltype),size(map,1)*size(map,2),[])', ...
        avg_map{curr_domain}(:)','correlation'),ctx_str_k_px_exp, ...
        domain_exp,celltype_exp,good_units_exp,'uni',false);
    
    % Act/avg act correlation
    avg_act_corr_exp(:,curr_domain) = ...
        cellfun(@(act,stim,rxn,outcome,domain,celltype,good_units) ...
        1-pdist2(permute(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:, ...
        good_units & domain == curr_domain & celltype == curr_celltype),1),[3,2,1]), ...
        avg_act{curr_domain},'correlation'),str_act_align, ...
        trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp, ...
        domain_exp,celltype_exp,good_units_exp,'uni',false);
    
    % Act/avg ctx act correlation
    avg_ctxact_corr_exp(:,curr_domain) = ...
        cellfun(@(act,stim,rxn,outcome,domain,celltype,good_units) ...
        1-pdist2(permute(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:, ...
        good_units & domain == curr_domain & celltype == curr_celltype),1),[3,2,1]), ...
        avg_ctxact{curr_domain},'correlation'),str_act_align, ...
        trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp, ...
        domain_exp,celltype_exp,good_units_exp,'uni',false);
    
end

figure;
for curr_domain = 1:n_aligned_depths
    
    % Plot average activity
    % (normalize activity: subtract baseline, divide std)
    % (!! assume first domain is stim-aligned !!)

    curr_act = cell2mat(cellfun(@(act,act_stimalign,baseline) ...
        (act-nanmean(act_stimalign(:,t < 0,:),2))./nanstd(act,[],2), ...
        avg_act(:,curr_domain),avg_act(:,1),'uni',false));
    curr_ctxact = cell2mat(cellfun(@(act,act_stimalign,baseline) ...
        (act-nanmean(act_stimalign(:,t < 0,:),2))./nanstd(act,[],2), ...
        avg_ctxact(:,curr_domain),avg_ctxact(:,1),'uni',false));
 
    subplot(n_aligned_depths,4,(curr_domain-1)*4+1); hold on;
    AP_errorfill(t,nanmean(curr_act,1)',AP_sem(curr_act,1)','k');
    AP_errorfill(t,nanmean(curr_ctxact,1)',AP_sem(curr_ctxact,1)',[0,0.7,0]);
    line([0,0],ylim,'color','r');
    ylabel('Activity');
    title(['Domain ' num2str(curr_domain)]);
    
    % Plot average map
    subplot(n_aligned_depths,4,(curr_domain-1)*4+2);
    imagesc(nanmean(cat(3,avg_map_exp{:,curr_domain}),3));
    axis image off
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'*RdBu'));
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    
    % Plot map/act correlation points
    subplot(n_aligned_depths,4,(curr_domain-1)*4+3); hold on;
    plot(cell2mat(avg_map_corr_exp(:,curr_domain)), ...
        cell2mat(avg_act_corr_exp(:,curr_domain)),'.k');
    plot(cell2mat(avg_map_corr_exp(:,curr_domain)), ...
        cell2mat(avg_ctxact_corr_exp(:,curr_domain)),'.','color',[0,0.7,0]);
    line([0,0],[-1,1],'linestyle','--');
    line([-1,1],[0,0],'linestyle','--');
    xlabel('Map correlation');
    ylabel('Activity correlation');
        
    % Plot correlation histogram
    hist_edges = linspace(-1,1,20);    
    subplot(n_aligned_depths,4,(curr_domain-1)*4+4); hold on;
    histogram(cell2mat(avg_act_corr_exp(:,curr_domain)),hist_edges,'normalization','probability');
    histogram(cell2mat(avg_ctxact_corr_exp(:,curr_domain)),hist_edges,'normalization','probability');
    histogram(cell2mat(avg_map_corr_exp(:,curr_domain)),hist_edges,'normalization','probability');
    xlabel('Correlation');
    ylabel('Fraction');
    legend({'Avg activity','Avg ctx activity','Maps'});
    
end

% TESTING 
figure; 
for curr_domain = 1:n_aligned_depths

    % Plot map/act correlation points
    subplot(n_aligned_depths,2,(curr_domain-1)*2+1); hold on;
    plot(cell2mat(avg_act_corr_exp(:,curr_domain)), ...
        cell2mat(avg_ctxact_corr_exp(:,curr_domain)),'.k');
    line([0,0],[-1,1],'linestyle','--');
    line([-1,1],[0,0],'linestyle','--');
    xlabel('Str activity correlation');
    ylabel('Ctx activity correlation');
    
    subplot(n_aligned_depths,2,(curr_domain-1)*2+2); hold on;
    plot([cell2mat(avg_act_corr_exp(:,curr_domain)), ...
        cell2mat(avg_ctxact_corr_exp(:,curr_domain))]','color',[0.5,0.5,0.5]);
    errorbar(nanmean([cell2mat(avg_act_corr_exp(:,curr_domain)), ...
        cell2mat(avg_ctxact_corr_exp(:,curr_domain))],1), ...
        AP_sem([cell2mat(avg_act_corr_exp(:,curr_domain)), ...
        cell2mat(avg_ctxact_corr_exp(:,curr_domain))],1),'color','k','linewidth',2);
    set(gca,'XTick',1:2,'XTickLabel',{'Str act','Ctx act'});
    xlim([0.5,2.5])
    ylabel('Correlation');
    
end




%% Sort by activity or map correlation

curr_domain = 1;
curr_celltype = 1;

% Map/avg map correlation
avg_map_corr_exp = cellfun(@(map,domain,celltype) ...
    1-pdist2( ...
    reshape(map,size(map,1)*size(map,2),[])', ...
    nanmean(reshape(map(:,:,domain == curr_domain & celltype == curr_celltype),size(map,1)*size(map,2),[]),2)', ...
    'correlation'),ctx_str_k_px_exp,domain_exp,celltype_exp,'uni',false);

% Act/avg act correlation
avg_act_corr_exp = cellfun(@(act,stim,rxn,outcome,domain,celltype) ...
    1-pdist2( ...
    permute(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:,:),1),[3,2,1]), ...
    permute(nanmean(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:,domain == curr_domain & celltype == curr_celltype),1),3),[3,2,1]), ...
    'correlation'),mua_allcat_stimalign_exp, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp, ...
    domain_exp,celltype_exp,'uni',false);

% Act/avg ctxact correlation
avg_ctxact_corr_exp = cellfun(@(act,ctx_act,stim,rxn,outcome,domain,celltype) ...
    1-pdist2( ...
    permute(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:,:),1),[3,2,1]), ...
    permute(nanmean(ctx_act(stim > 0 & rxn < 0.5 & outcome == 1,:,curr_domain),1),[3,2,1]), ...
    'correlation'),mua_allcat_stimalign_exp,fluor_kernelroi_deconv_stimalign_exp, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp, ...
    domain_exp,celltype_exp,'uni',false);


use_measure = cell2mat(avg_ctxact_corr_exp);
plot_cells_idx = find(good_units_allcat & celltype_allcat == 1 & domain_aligned_allcat == 2);
[~,sort_idx_cells] = sort(use_measure(plot_cells_idx));
sort_idx = plot_cells_idx(sort_idx_cells);





%% Task responses vs map

% Get task activity for each cell
% (stim)
stim_act_cond = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    [squeeze(nanmean(nanmean(act(stim > 0,t > 0 & t < 0.15,:),2),1)), ...
    squeeze(nanmean(nanmean(act(stim > 0,t < 0,:),2),1))], ...
    mua_allcat_stimalign_exp,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false));
vis_act = (stim_act_cond(:,1) - stim_act_cond(:,2))./ ...
    (sum(abs(stim_act_cond),2));

% (move)
move_act_cond = cell2mat(cellfun(@(act_move,act_stim,stim,rxn,outcome) ...
    [squeeze(nanmean(nanmean(act_move(stim <= 0,t > 0 & t < 0.2,:),2),1)), ...
    squeeze(nanmean(nanmean(act_stim(stim <= 0,t < 0,:),2),1))], ...
    mua_allcat_movealign_exp_smooth,mua_allcat_stimalign_exp_smooth, ...
    trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false));
move_act = (move_act_cond(:,1) - move_act_cond(:,2))./ ...
    (sum(abs(move_act_cond),2));

% (reward)
use_t = t > 0 & t < 0.5;

reward_act_cond = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    [squeeze(nanmean(nanmean(act(outcome == 1,use_t,:),2),1)), ...
    squeeze(nanmean(nanmean(act(outcome == -1,use_t,:),2),1))], ...
    mua_allcat_outcomealign_exp_smooth,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false));

reward_act = (reward_act_cond(:,1) - reward_act_cond(:,2))./ ...
    (sum(abs(reward_act_cond),2));



% Get correlation between map and average map instead
k_corr = 1-pdist2(reshape(ctx_str_k_px,[],size(ctx_str_k_px,3))', ...
    reshape(kernel_roi.max_weighted,[],n_aligned_depths)','correlation');


% Plot weights vs visual activity
figure;

curr_domain = 1;
curr_act = vis_act;

plot_cells = good_units_allcat & domain_aligned_allcat == curr_domain & ...
    celltype_allcat == 1 & abs(curr_act) ~= 1 & ~isnan(curr_act);

subplot(1,2,1);
plot(curr_act(plot_cells),k_corr(plot_cells(good_units_allcat),curr_domain),'.k');
xlabel('Visual activity');
ylabel('ROI weight');
line(repmat(nanmean(vis_act_ctx),2,1),ylim,'color',[0,0.7,0],'linestyle','--','linewidth',2)
legend({'Str neuron','Cortex mean'});

subplot(1,2,2); hold on;
plot(sort(curr_act(plot_cells)),[1:sum(plot_cells)]./sum(plot_cells),'k','linewidth',2);
plot(sort(k_corr(plot_cells(good_units_allcat),curr_domain)), ...
    [1:sum(plot_cells)]./sum(plot_cells),'r','linewidth',2);
line(repmat(nanmean(vis_act_ctx),2,1),ylim,'color',[0,0.7,0],'linestyle','--','linewidth',2);
xlabel('Visual activity');
ylabel('Cumulative sum');
legend({'Visual activity (str neuron)','Cortex weight','Cortex mean'})



% Plot type of activity for each domain
figure; 
p1 = subplot(1,2,1); hold on; set(gca,'ColorOrder',copper(3));
xlabel('Condition activity');
ylabel('Cumulative sum');
p2 = subplot(1,2,2); hold on; set(gca,'ColorOrder',copper(3));
xlabel('Map correlation');
ylabel('Cumulative sum');
for curr_domain = 1:n_aligned_depths
    switch curr_domain
        case 1
            curr_act = vis_act;
        case 2
            curr_act = move_act;
        case 3
            curr_act = reward_act;
    end
    plot_cells = good_units_allcat & domain_aligned_allcat == curr_domain & ...
        celltype_allcat == 1 & abs(curr_act) ~= 1 & ~isnan(curr_act);
    plot(p1,sort(curr_act(plot_cells)),[1:sum(plot_cells)]./sum(plot_cells),'linewidth',2); 
    plot(p2,sort(k_corr(plot_cells(good_units_allcat),curr_domain)), ...
    [1:sum(plot_cells)]./sum(plot_cells),'linewidth',2);
end
legend(p1,{'DMS/Vis','DCS/Move','DLS/Reward'});
legend(p2,{'DMS','DCS','DLS'});



% Histogram of response types?
figure; hold on;
for act_cond = 1:3
    switch act_cond
        case 1
            curr_act = vis_act;
        case 2
            curr_act = move_act;
        case 3
            curr_act = reward_act;
    end
    subplot(1,3,act_cond); hold on; set(gca,'ColorOrder',copper(3));
    for curr_domain = 1:n_aligned_depths
        plot_cells = good_units_allcat & domain_aligned_allcat == curr_domain & ...
            celltype_allcat == 1 & abs(curr_act) ~= 1 & ~isnan(curr_act);
        hist_edges = [linspace(0,1,20)];
        histogram(curr_act(plot_cells),hist_edges,'normalization','probability')
    end
end


% Cumulative distribution of response types?
plot_celltype = 1;
figure; hold on;
for act_cond = 1:3
    switch act_cond
        case 1
            curr_act = vis_act;
        case 2
            curr_act = move_act;
        case 3
            curr_act = reward_act;
    end
    subplot(1,3,act_cond); hold on; set(gca,'ColorOrder',copper(3));
    for curr_domain = 1:n_aligned_depths
        plot_cells = good_units_allcat & domain_aligned_allcat == curr_domain & ...
            celltype_allcat == plot_celltype & abs(curr_act) ~= 1 & ~isnan(curr_act);
        plot(sort(curr_act(plot_cells)),[1:sum(plot_cells)]./sum(plot_cells),'linewidth',2);
    end
end



AP_imscroll(ctx_str_k_px_maxnorm(:,:,plot_cells(good_units_allcat)));
axis image;
caxis([-1,1])
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);


%% [compare above] task respones in cortex?
% (this is garbage)

curr_domain = 1;
stim_act_cond_ctx = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    [squeeze(nanmean(nanmean(act(stim > 0,t > 0 & t < 0.15,curr_domain),2),1)), ...
    squeeze(nanmean(nanmean(act(stim > 0,t < 0,curr_domain),2),1))], ...
    fluor_kernelroi_deconv_exp,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false));

vis_act_ctx = stim_act_cond_ctx(:,1)./ctx_act_baseline_std(:,curr_domain);


vis_act = stim_act_cond(:,1)./str_act_baseline_std;


curr_domain = 2;
move_act_cond_ctx = cell2mat(cellfun(@(act_move,act_stim,stim,rxn,outcome) ...
    [squeeze(nanmean(nanmean(act_move(stim <= 0,t > 0 & t < 0.2,curr_domain),2),1)), ...
    squeeze(nanmean(nanmean(act_stim(stim <= 0,t < 0,curr_domain),2),1))], ...
    fluor_kernelroi_deconv_movealign_exp,fluor_kernelroi_deconv_exp, ...
    trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false));
move_act_ctx = move_act_cond_ctx(:,1)./ctx_act_baseline_std(:,curr_domain);

move_act = move_act_cond(:,1)./str_act_baseline_std;


figure;
vis_act_rel = vis_act(good_units_allcat & celltype_allcat == 1 & ...
    domain_aligned_allcat == 1 & ~isinf(vis_act))./nanmean(vis_act_ctx);
move_act_rel = move_act(good_units_allcat & celltype_allcat == 1 & ...
    domain_aligned_allcat == 2 & ~isinf(move_act))./nanmean(move_act_ctx);
subplot(1,2,1); hold on;
histogram(vis_act_rel,hist_edges,'normalization','probability');
histogram(move_act_rel,hist_edges,'normalization','probability');

subplot(1,2,2);hold on;
plot(sort(vis_act_rel),[1:length(vis_act_rel)]./length(vis_act_rel),'linewidth',2); 
plot(sort(move_act_rel),[1:length(move_act_rel)]./length(move_act_rel),'linewidth',2); 





%% TESTING: visual response vs. map

% Get cortical activity in kernel ROI

% Concatenate cortex data
fluor_allcat = cell2mat(vertcat(fluor_all{:}));
% fluor_allcat_deconv = AP_deconv_wf(fluor_allcat);
% fluor_allcat_deconv_baseline = nanmean(reshape(fluor_allcat_deconv(:,t_baseline,:),[],1,n_vs));
fluor_allcat_deconv = fluor_allcat_deconv - fluor_allcat_deconv_baseline;

kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
load(kernel_roi_fn);
n_kernel_rois = size(kernel_roi.bw,3);
fluor_kernelroi_deconv = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],kernel_roi.bw), ...
    size(kernel_roi.bw,3),[],size(fluor_allcat_deconv,1)),[3,2,1]);
fluor_kernelroi_deconv_exp = mat2cell(fluor_kernelroi_deconv,use_split,length(t),n_kernel_rois);

curr_roi = kernel_roi.bw(:,:,1);

% Get correlation between fluorescence and every cell
curr_domain = 1;
str_ctx_corr = cell2mat(cellfun(@(spikes,fluor) 1-pdist2( ...
    reshape(permute(spikes,[2,1,3]),[],size(spikes,3))', ...
    reshape(fluor(:,:,curr_domain)',[],1)','correlation'), ...
    mua_allcat_stimalign_exp,fluor_kernelroi_deconv_exp,'uni',false));

% Get weight in kernel ROI
ctx_str_k_px_maxnorm = ctx_str_k_px./nanmax(nanmax(abs(ctx_str_k_px),[],1),[],2);
ctx_str_k_px_maxnorm_roi = ...
    curr_roi(:)'* ...
    reshape(ctx_str_k_px_maxnorm,[],size(ctx_str_k_px_maxnorm,3))./ ...
    sum(curr_roi(:));

ctx_str_trialk_px_roi =  ...
    curr_roi(:)'* ...
    reshape(ctx_str_trialk_px,[],size(ctx_str_trialk_px,3))./ ...
    sum(curr_roi(:));


% Visual response (different from above)
use_t = t > 0.05 & t < 0.15;

stim_act_cond = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    [squeeze(nanmean(nanmean(act(stim > 0 & rxn < 0.5,use_t,:),2),1)), ...
    squeeze(nanmean(nanmean(act(stim <= 0 & rxn < 0.5,use_t,:),2),1))], ...
    mua_allcat_stimalign_exp,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false));

vis_act = (stim_act_cond(:,1) - stim_act_cond(:,2))./ ...
    (sum(stim_act_cond,2));

% Visual response of the cortex
use_t = t > 0 & t < 0.1;

curr_domain = 1;

stim_act_cond_ctx = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
    [squeeze(nanmean(nanmean(act(stim > 0 & rxn < 0.5,use_t,curr_domain),2),1)), ...
    squeeze(nanmean(nanmean(act(stim <= 0 & rxn < 0.5,use_t,curr_domain),2),1))], ...
    fluor_kernelroi_deconv_exp,trial_stim_allcat_exp, ...
    move_t_exp,trial_outcome_allcat_exp,'uni',false));

vis_act_ctx = (stim_act_cond_ctx(:,1) - stim_act_cond_ctx(:,2))./ ...
    (sum(abs(stim_act_cond_ctx),2));


% Plot in domain 1
figure;
plot_cells = good_units_allcat & domain_aligned_allcat == 1 & celltype_allcat == 1 & abs(vis_act) ~= 1 & ~isnan(vis_act);
subplot(2,2,1);
plot(vis_act(plot_cells),ctx_str_k_px_maxnorm_roi(plot_cells(good_units_allcat)),'.k');
xlabel('Visual activity');
ylabel('ROI weight');
line(repmat(nanmean(vis_act_ctx),2,1),ylim,'color',[0,0.7,0],'linestyle','--','linewidth',2)
legend({'Str neuron','Cortex mean'});

subplot(2,2,2); hold on;
plot(sort(vis_act(plot_cells)),[1:sum(plot_cells)]./sum(plot_cells),'k','linewidth',2);
plot(sort(ctx_str_k_px_maxnorm_roi(plot_cells(good_units_allcat))),[1:sum(plot_cells)]./sum(plot_cells),'r','linewidth',2);
plot(sort(vis_act_ctx),[1:length(vis_act_ctx)]./length(vis_act_ctx),'color',[0,0.7,0],'linewidth',2);
line(repmat(nanmean(vis_act_ctx),2,1),ylim,'color',[0,0.7,0],'linestyle','--','linewidth',2)
xlabel('Visual activity');
ylabel('Cumulative sum');
legend({'Visual activity (str neuron)','Cortex weight','Visual activity (cortex ROI)','Cortex mean'})

subplot(2,2,3);
plot(vis_act(plot_cells),str_ctx_corr(plot_cells),'.k');
xlabel('Visual activity');
ylabel('Ctx-str corr');
line(repmat(nanmean(vis_act_ctx),2,1),ylim,'color',[0,0.7,0],'linestyle','--','linewidth',2)
legend({'Str neuron','Cortex mean'});

subplot(2,2,4); hold on;
plot(sort(vis_act(plot_cells)),[1:sum(plot_cells)]./sum(plot_cells),'k','linewidth',2);
plot(sort(str_ctx_corr(plot_cells)),[1:sum(plot_cells)]./sum(plot_cells),'r','linewidth',2);
plot(sort(vis_act_ctx),[1:length(vis_act_ctx)]./length(vis_act_ctx),'color',[0,0.7,0],'linewidth',2);
line(repmat(nanmean(vis_act_ctx),2,1),ylim,'color',[0,0.7,0],'linestyle','--','linewidth',2)
xlabel('Visual activity');
ylabel('Cumulative sum');
legend({'Visual activity (str neuron)','Str-ctx corr','Visual activity (cortex ROI)','Cortex mean'})


% Plot things sorted

[~,sort_idx_full] = sort(vis_act);
a = find(plot_cells);
b = ismember(sort_idx_full,a);
sort_idx = sort_idx_full(b);

figure;plot(vis_act(sort_idx),'.k')
xlabel('Sorted cell');
ylabel('Visual activity');

cell_string = arrayfun(@(x) ...
    {[' Sort idx: ' num2str(x) ...
    ', Cell: ' num2str(sort_idx(x)) ...
    ', Domain: ' num2str(domain_aligned_allcat(sort_idx(x))) ...
    ', ' celltype_labels{celltype_allcat(sort_idx(x))}],...
    ['Animal: ' num2str(animals_allcat(sort_idx(x))) ...
    ', Day: ' num2str(days_allcat(sort_idx(x))) ...
    ', Neuron: ' num2str(neurons_allcat(sort_idx(x)))]}, ...
    1:length(sort_idx),'uni',false);

AP_imscroll(mua_allcat_stimalign_exp_pad(:,:,sort_idx),cell_string);
line(repmat(find(t >= 0,1),2,1),ylim,'linestyle','--','color','r');
colormap(brewermap([],'Greys'));
caxis([0,50]);


AP_imscroll(ctx_str_k_px_maxnorm(:,:,sort_idx));
axis image;
caxis([-1,1])
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);


AP_imscroll(mua_ctx_sta_px_maxnorm(:,:,sort_idx));
axis image;
caxis([-1,1])
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);


AP_imscroll(ctx_str_trialk_px(:,:,sort_idx));
axis image;
caxis([-1,1])
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);



% Kenneth suggested: sort activity by fluor activity to show that it's
% correlated with spikes?


% Pad and select trials of fluor like spikes
plot_trials_1 = cellfun(@(stim,rxn,outcome) ...
    stim > 0, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false);
plot_trials_2 = cellfun(@(stim,rxn,outcome) ...
    stim < 0, ...
    trial_stim_allcat_exp,move_t_exp,trial_outcome_allcat_exp,'uni',false);

max_trials = max(cellfun(@(x,y) sum(x)+sum(y),plot_trials_1,plot_trials_2));

curr_roi = 1;
fluor_kernelroi_deconv_exp_pad = cell2mat(permute(cellfun(@(act,trials_1,trials_2) ...
    padarray(act([find(trials_1);find(trials_2)],:,curr_roi), ...
    [max_trials-(sum(trials_1)+sum(trials_2)),0,0],NaN,'post'), ...
    fluor_kernelroi_deconv_exp,plot_trials_1,plot_trials_2,'uni',false),[2,3,1]));



% get exp index for cell
curr_cell = 10244;

a = cellfun(@(x) size(x,3),mua_allcat_stimalign_exp_smooth);
b = cumsum(a);
curr_exp = find(b >= curr_cell,1,'first');

[~,max_idx] = max(fluor_kernelroi_deconv_exp_pad(:,t > -0.5,curr_exp),[],2);
[~,trial_sort_idx] = sort(max_idx);

figure; hold on;
imagesc(fluor_kernelroi_deconv_exp_pad(trial_sort_idx,:,curr_exp));
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));

[y,x] = find(mua_allcat_stimalign_exp_pad(trial_sort_idx,:,curr_cell) > 0);
plot(x,y,'.k');
title(num2cell(curr_cell))





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
    ylabel('R^2');
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
    ylabel('Trial-average R^2');
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
    mua_allcat_exp,'uni',false));

% Plot each regressor by domain/celltype
for curr_regressor = 1:n_regressors
    
    figure('Name',task_regressor_labels{curr_regressor});
    for curr_domain = 1:n_aligned_depths
        for curr_celltype = 1:n_celltypes
            
            curr_cells = good_units_allcat  &...
                domain_aligned_allcat == curr_domain & ...
                celltype_allcat == curr_celltype;
            
            subplot(n_aligned_depths,n_celltypes, ...
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
plot_cells_1 = good_units_allcat & stim_cells & domain_aligned_allcat == 1  & celltype_allcat == 1;
plot_cells_2 = good_units_allcat & ~stim_cells & ~move_cells & ~reward_cells & domain_aligned_allcat == 1  & celltype_allcat == 1;
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


%% ~~~~~~ Paper figures/claims: single units

%% Fig 4c: Average activity within and across cells by domain/celltype

mua_exp = vertcat(mua_all{:});

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));

% Align activity by move/outcome (natively stim-aligned)

% (get indicies for alignments - used later)
stim_align = zeros(size(trial_stim_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;
use_align = {stim_align,move_align,outcome_align};

% (move aligned)
move_idx_exp = mat2cell(move_idx,use_split,1);
mua_exp_movealign = vertcat(mua_all{:});
for curr_exp = 1:length(mua_exp_movealign)
   for curr_trial = 1:size(mua_exp_movealign{curr_exp},1)
       mua_exp_movealign{curr_exp}(curr_trial,:,:) = ...
           circshift(mua_exp_movealign{curr_exp}(curr_trial,:,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
   end
end

% (outcome aligned)
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
mua_exp_outcomealign = vertcat(mua_all{:});
for curr_exp = 1:length(mua_exp_outcomealign)
    for curr_trial = 1:size(mua_exp_outcomealign{curr_exp},1)
        mua_exp_outcomealign{curr_exp}(curr_trial,:,:) = ...
            circshift(mua_exp_outcomealign{curr_exp}(curr_trial,:,:), ...
            -outcome_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
    end
end

% Get average activity across alignments
use_align_labels = {'Stim','Move onset','Outcome'};

% (split trials for sort/plot)
trial_rand = rand(max(cellfun(@(x) size(x,1),mua_exp)),1);

act_mean_sort = nan(length(good_units_allcat),length(t),length(use_align_labels));
act_mean_plot = nan(length(good_units_allcat),length(t),length(use_align_labels));
for curr_align = 1:length(use_align_labels)
    
    switch curr_align
        case 1
            curr_mua = mua_exp;
        case 2
            curr_mua = mua_exp_movealign;
        case 3
            curr_mua = mua_exp_outcomealign;
    end
    
    act_mean_sort(:,:,curr_align) = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
        squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1 & ...
        trial_rand(1:size(act,1)) < 0.5,:,:),1)), ...
        curr_mua,trial_stim_allcat_exp,move_t_exp, ...
        trial_outcome_allcat_exp,'uni',false)')';
    
    act_mean_plot(:,:,curr_align) = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
        squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1 & ...
        trial_rand(1:size(act,1)) >= 0.5,:,:),1)), ...
        curr_mua,trial_stim_allcat_exp,move_t_exp, ...
        trial_outcome_allcat_exp,'uni',false)')';  
    
end


% Set align breakpoints
align_col = [1,0,0;0.8,0,0.8;0,0,0.8];
% (split the alignment halfway between median alignment points)
align_median = cellfun(@(x) -nanmedian(x)/sample_rate,use_align);
align_break = align_median(1:end-1) + diff(align_median*0.8);
align_t = {[-0.05,align_break(1)],[align_break(1:2)],[align_break(2),1]};

% Plot aligned_activity (all cells)
figure;
for curr_depth = 1:n_aligned_depths
    for curr_celltype = 1:n_celltypes
        
        curr_cells = domain_aligned_allcat == curr_depth & ...
            celltype_allcat == curr_celltype & good_units_allcat & ...
            any(any(act_mean_plot,3),2);
        curr_cells_idx = find(curr_cells);
        [~,max_idx] = max(act_mean_sort(curr_cells,:,1),[],2);
        [~,sort_idx] = sort(max_idx,'descend');
   
        for curr_align = 1:length(use_align)
            
            subplot(n_aligned_depths,n_celltypes, ...
            (curr_depth-1)*n_celltypes+curr_celltype); hold on;
            
            curr_t_offset = -nanmedian(use_align{curr_align})/sample_rate;
            curr_t = t + curr_t_offset;
            curr_t_plot = curr_t >= align_t{curr_align}(1) & ...
                curr_t <= align_t{curr_align}(2);
            
            curr_act_sorted = act_mean_plot(curr_cells_idx(sort_idx),:,curr_align);
            curr_act_sorted_norm = curr_act_sorted./max(curr_act_sorted,[],2);
        
            % smooth if too many cells to plot accurately
            if sum(curr_cells) > 500
                n_smooth = 5;
                curr_act_sorted_norm = convn(curr_act_sorted_norm, ...
                    ones(n_smooth,1)./n_smooth,'same');
            end

            plot_t = curr_t > align_t{curr_align}(1) & curr_t <= align_t{curr_align}(2);
            
            imagesc(curr_t(plot_t),[],curr_act_sorted_norm(:,plot_t));
            caxis([0,1]);
            axis tight;
            line(repmat(curr_t_offset,2,1),ylim,'color',align_col(curr_align,:));
            colormap(gca,brewermap([],'Greys'));
            ylabel('Neuron (sorted)');
            xlabel('~Time from stim');
            title(['Str ' num2str(curr_depth) ' ' celltype_labels{curr_celltype}]);
            
        end       
    end
end

% Plot aligned_activity (average across cells)
figure;
for curr_depth = 1:n_aligned_depths
    for curr_celltype = 1:n_celltypes
        
        curr_cells = domain_aligned_allcat == curr_depth & ...
            celltype_allcat == curr_celltype & good_units_allcat & ...
            any(any(act_mean_plot,3),2);   
   
        for curr_align = 1:length(use_align)
            
            subplot(n_aligned_depths,n_celltypes, ...
            (curr_depth-1)*n_celltypes+curr_celltype); hold on;
            
            curr_t_offset = -nanmedian(use_align{curr_align})/sample_rate;
            curr_t = t + curr_t_offset;
            curr_t_plot = curr_t >= align_t{curr_align}(1) & ...
                curr_t <= align_t{curr_align}(2);
            
            curr_act = act_mean_plot(curr_cells,:,curr_align);
            
            plot_t = curr_t > align_t{curr_align}(1) & curr_t <= align_t{curr_align}(2);
            
            AP_errorfill(curr_t(plot_t)', ...
                nanmean(curr_act(:,plot_t),1)', ...
                AP_sem(curr_act(:,plot_t),1)','k');
            
            axis tight;
            line(repmat(curr_t_offset,2,1),ylim,'color',align_col(curr_align,:));
            colormap(gca,brewermap([],'Greys'));
            ylabel('Spikes/s');
            xlabel('~Time from stim');
            title(['Str ' num2str(curr_depth) ' ' celltype_labels{curr_celltype}]);
            
        end       
    end
end


%% Fig 4d: SUA-celltype/cortex activity correlation

% Get trial params by experiment
trial_stim_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
move_t_exp = mat2cell(move_t,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));

% Align activity by move/outcome (natively stim-aligned)

% (get indicies for alignments - used later)
stim_align = zeros(size(trial_stim_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;
use_align = {stim_align,move_align,outcome_align};

% (stim aligned)
mua_exp = vertcat(mua_all{:});

% (move aligned)
move_idx_exp = mat2cell(move_idx,use_split,1);
mua_exp_movealign = vertcat(mua_all{:});
for curr_exp = 1:length(mua_exp_movealign)
   for curr_trial = 1:size(mua_exp_movealign{curr_exp},1)
       mua_exp_movealign{curr_exp}(curr_trial,:,:) = ...
           circshift(mua_exp_movealign{curr_exp}(curr_trial,:,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
   end
end

move_idx_exp = mat2cell(move_idx,use_split,1);
fluor_kernelroi_deconv_exp_movealign = fluor_kernelroi_deconv_exp;
for curr_exp = 1:length(fluor_kernelroi_deconv_exp_movealign)
   for curr_trial = 1:size(fluor_kernelroi_deconv_exp_movealign{curr_exp},1)
       fluor_kernelroi_deconv_exp_movealign{curr_exp}(curr_trial,:,:) = ...
           circshift(fluor_kernelroi_deconv_exp_movealign{curr_exp}(curr_trial,:,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
   end
end

% (outcome aligned)
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
mua_exp_outcomealign = vertcat(mua_all{:});
for curr_exp = 1:length(mua_exp_outcomealign)
    for curr_trial = 1:size(mua_exp_outcomealign{curr_exp},1)
        mua_exp_outcomealign{curr_exp}(curr_trial,:,:) = ...
            circshift(mua_exp_outcomealign{curr_exp}(curr_trial,:,:), ...
            -outcome_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
    end
end

outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
fluor_kernelroi_deconv_exp_outcomealign = fluor_kernelroi_deconv_exp;
for curr_exp = 1:length(fluor_kernelroi_deconv_exp_outcomealign)
    for curr_trial = 1:size(fluor_kernelroi_deconv_exp_outcomealign{curr_exp},1)
        fluor_kernelroi_deconv_exp_outcomealign{curr_exp}(curr_trial,:,:) = ...
            circshift(fluor_kernelroi_deconv_exp_outcomealign{curr_exp}(curr_trial,:,:), ...
            -outcome_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
    end
end

% Get average activity across alignments
use_align_labels = {'Stim','Move onset','Outcome'};

str_unit_act_mean = nan(length(good_units_allcat),length(t),length(use_align_labels));
for curr_align = 1:length(use_align_labels)
    
    switch curr_align
        case 1
            curr_mua = mua_exp;
        case 2
            curr_mua = mua_exp_movealign;
        case 3
            curr_mua = mua_exp_outcomealign;
    end
    
    str_unit_act_mean(:,:,curr_align) = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
        squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:,:),1)), ...
        curr_mua,trial_stim_allcat_exp,move_t_exp, ...
        trial_outcome_allcat_exp,'uni',false)')';
    
end

n_recordings = length(fluor_kernelroi_deconv_exp);
ctx_kernelroi_act_mean = nan(n_aligned_depths,length(t),n_recordings,length(use_align_labels));
for curr_align = 1:length(use_align_labels)
    
    switch curr_align
        case 1
            curr_act = fluor_kernelroi_deconv_exp;
        case 2
            curr_act = fluor_kernelroi_deconv_exp_movealign;
        case 3
            curr_act = fluor_kernelroi_deconv_exp_outcomealign;
    end
    
    ctx_kernelroi_act_mean(:,:,:,curr_align) = cell2mat(permute(cellfun(@(act,stim,rxn,outcome) ...
        squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:,:),1))', ...
        curr_act,trial_stim_allcat_exp,move_t_exp, ...
        trial_outcome_allcat_exp,'uni',false),[2,3,1]));
    
end

% Concatenate multi-aligned averages
% (split the alignment halfway between median alignment points)
align_median = cellfun(@(x) -nanmedian(x)/sample_rate,use_align);
align_break = align_median(1:end-1) + diff(align_median*0.8);
align_t = {[-0.05,align_break(1)],[align_break(1:2)],[align_break(2),1]};

align_samples = cell(size(align_t));
for curr_align = 1:length(align_samples)    
    curr_t_offset = -nanmedian(use_align{curr_align})/sample_rate;
    curr_t = t + curr_t_offset;
    curr_t_use = curr_t >= align_t{curr_align}(1) & ...
        curr_t <= align_t{curr_align}(2);
    align_samples{curr_align} = curr_t_use;
end

str_unit_actmean_multialign = cell2mat(arrayfun(@(x) ...
    str_unit_act_mean(:,align_samples{x},x),1:size(str_unit_act_mean,3),'uni',false));

ctx_kernelroi_actmean_multialign = cell2mat(arrayfun(@(x) ...
    ctx_kernelroi_act_mean(:,align_samples{x},:,x),1:size(ctx_kernelroi_act_mean,4),'uni',false));


% Get cell-average correlations (across recording)
n_recordings = max(recordings_allcat);
celltype_act_corr = nan(size(str_unit_actmean_multialign,1),3);
for curr_depth = 1:n_aligned_depths
    for curr_celltype = 1:3
        for curr_celltype_compare = 1:3
            for curr_recording = 1:n_recordings
                
                % Get cells in domain in/out of recording
                curr_cells_currday = ...
                    good_units_allcat & ...
                    domain_aligned_allcat == curr_depth & ...
                    celltype_allcat == curr_celltype & ...
                    recordings_allcat == curr_recording;
                
                curr_cells_otherday = ...
                    good_units_allcat & ...
                    domain_aligned_allcat == curr_depth & ...
                    celltype_allcat == curr_celltype_compare & ...
                    recordings_allcat ~= curr_recording;
                act_mean_multialign_otherday_avg = ...
                    nanmean(str_unit_actmean_multialign(curr_cells_otherday,:),1);
                
                if any(curr_cells_currday)
                    % Correlate day's cells with other-day average
                    celltype_act_corr(curr_cells_currday,curr_celltype_compare) = ...
                        corr(str_unit_actmean_multialign(curr_cells_currday,:)', ...
                        act_mean_multialign_otherday_avg');
                end
            end
        end
    end
end

% Get cell-cortex ROI correlation (across recording)
n_recordings = max(recordings_allcat);
celltype_ctxact_corr = nan(size(str_unit_actmean_multialign,1),1);
for curr_depth = 1:n_aligned_depths
    for curr_celltype = 1:3
        for curr_recording = 1:n_recordings
            
            % Get cells in domain and cortex out of recording
            curr_cells_currday = ...
                good_units_allcat & ...
                domain_aligned_allcat == curr_depth & ...
                celltype_allcat == curr_celltype & ...
                recordings_allcat == curr_recording;
            
            ctx_act_otherday = ...
                nanmean(ctx_kernelroi_actmean_multialign(curr_depth,:, ...
                setdiff(1:n_recordings,curr_recording)),3);
            
            if any(curr_cells_currday)
                % Correlate day's cells with other-day average
                celltype_ctxact_corr(curr_cells_currday) = ...
                    corr(str_unit_actmean_multialign(curr_cells_currday,:)', ...
                    ctx_act_otherday');
            end
        end
    end
end



% Concatenate celltype and fluorescence correlation
celltype_allcorr = [celltype_act_corr,celltype_ctxact_corr];

% Get correlation mean and in FR percentile bins
fr = nanmean(str_unit_actmean_multialign,2);
n_prctiles = 4;
fr_prctilebin_edges = linspace(0,100,n_prctiles+1);
fr_prctilebin_centers = fr_prctilebin_edges(1:end-1)+diff(fr_prctilebin_edges)./2;

fr_prctilebins = nan(size(fr));
for curr_celltype = 1:3
    curr_units = good_units_allcat & celltype_allcat == curr_celltype;
    curr_fr = fr(curr_units);
    curr_fr_prctilebin_edges = prctile(curr_fr,fr_prctilebin_edges);
    fr_prctilebins(curr_units) = discretize(curr_fr,curr_fr_prctilebin_edges);
end

use_units = ~any(isnan(celltype_allcorr),2) & ~isnan(fr_prctilebins);

% Make grouping variables and average by group
grp_matrix = [celltype_allcat(use_units), ...
    fr_prctilebins(use_units), ...
    recordings_allcat(use_units)];
grp_matrix_repmat = reshape(permute(repmat( ...
    grp_matrix,1,1,size(celltype_allcorr,2)),[1,3,2]),[],size(grp_matrix,2));

[~,celltype_compare_repmat] = ndgrid(1:size(celltype_allcorr,1),1:size(celltype_allcorr,2));
celltype_compare_grp = reshape(celltype_compare_repmat(use_units,:),[],1);

celltype_actcorr_frbins = accumarray( ...
    [grp_matrix_repmat,celltype_compare_grp], ...
    reshape(celltype_allcorr(use_units,:),[],1), ...
    [3,length(fr_prctilebin_centers),max(recordings_allcat),size(celltype_allcorr,2)],@nanmean,NaN);
celltype_fr_frbins = accumarray(grp_matrix,fr(use_units), ...
    [3,length(fr_prctilebin_centers),max(recordings_allcat)],@nanmean,NaN);


figure;
celltype_col = ...
    [0.9,0.4,0.6; ...
    0.4,0.5,0.6; ...
    0.5,0.3,0.1; ...
    1,0.5,0];
ctx_col = [0,0.7,0];
p = gobjects(2,3);
for curr_celltype = 1:3
    
    curr_fr = permute(squeeze(celltype_fr_frbins(curr_celltype,:,:,:)),[1,3,2]);
    curr_actcorr = permute(squeeze(celltype_actcorr_frbins(curr_celltype,:,:,:)),[1,3,2]);
    
    p(1,curr_celltype) = subplot(2,3,curr_celltype); hold on;
    set(gca,'ColorOrder',[celltype_col(1:3,:);ctx_col],'YScale','log');
    errorbar(nanmean(curr_actcorr,3), ...
        repmat(nanmean(curr_fr,3),1,size(celltype_allcorr,2)), ...
        AP_sem(curr_actcorr,3),'horizontal','linewidth',2);
    title(celltype_labels{curr_celltype})
    legend([celltype_labels(1:3),'Cortex ROI']);
    xlabel('Avg act corr');
    ylabel('Firing rate');

    p(2,curr_celltype) = subplot(2,3,3+curr_celltype); hold on;
    set(gca,'ColorOrder',[celltype_col(1:3,:);ctx_col],'XScale','log');
    errorbar(repmat(nanmean(curr_fr,3),1,size(celltype_allcorr,2)), ...
        nanmean(curr_actcorr,3),AP_sem(curr_actcorr,3),'linewidth',2);
    title(celltype_labels{curr_celltype})
    legend([celltype_labels(1:3),'Cortex ROI']);
    xlabel('Firing rate');
    ylabel('Avg act corr');
    axis tight
    curr_xlim = xlim;
    xlim([curr_xlim(1)*0.7,curr_xlim(2)*1.3]);
    
end
linkaxes(p(1,:),'xy');
linkaxes(p(2,:),'y');


% Plot histogram of firing rates
fr = nanmean(str_unit_actmean_multialign,2);
fr_bin_edges = [logspace(0,2,50),Inf];

figure; hold on
celltype_col = ...
    [0.9,0.4,0.6; ...
    0.4,0.5,0.6; ...
    0.5,0.3,0.1; ...
    1,0.5,0];
set(gca,'ColorOrder',celltype_col,'XScale','log');
histogram(fr(good_units_allcat & celltype_allcat == 1),fr_bin_edges,'Normalization','probability');
histogram(fr(good_units_allcat & celltype_allcat == 2),fr_bin_edges,'Normalization','probability');
histogram(fr(good_units_allcat & celltype_allcat == 3),fr_bin_edges,'Normalization','probability');
legend(celltype_labels(1:3));
xlabel('Firing rate');
ylabel('Fraction');



% (statistics comparing celltypes/cortex)

disp('MSN/FSI/CTX ANOVA:')
for curr_celltype = 1:3   
    curr_actcorr = permute(celltype_actcorr_frbins(curr_celltype,:,:,:),[2,4,3,1]);
    [fr_grp,compare_grp,recording_grp] = ndgrid( ...
        1:size(curr_actcorr,1), ...
        1:size(curr_actcorr,2), ...
        1:size(curr_actcorr,3));
    use_comparison = [1,2,4]; % (msn v fsi v cortex)
    p = anovan(reshape(curr_actcorr(:,use_comparison,:),[],1), ...
        [reshape(fr_grp(:,use_comparison,:),[],1), ...
        reshape(compare_grp(:,use_comparison,:),[],1)], ...
        'model','interaction','display','off');
    disp([celltype_labels{curr_celltype} ' v (msn/fsi/ctx): p(fr,type,interaction) = ' num2str(p')])    
end

disp('TAN ANOVA:')
for curr_celltype = 1:2
    curr_actcorr = permute(celltype_actcorr_frbins(curr_celltype,:,:,:),[2,4,3,1]);
    [fr_grp,compare_grp,recording_grp] = ndgrid( ...
        1:size(curr_actcorr,1), ...
        1:size(curr_actcorr,2), ...
        1:size(curr_actcorr,3));
    use_comparison = [curr_celltype,3]; % (msn v fsi v cortex)
    p = anovan(reshape(curr_actcorr(:,use_comparison,:),[],1), ...
        [reshape(fr_grp(:,use_comparison,:),[],1), ...
        reshape(compare_grp(:,use_comparison,:),[],1)], ...
        'model','interaction','display','off');
    disp([celltype_labels{curr_celltype} ' v (tan): p(fr,type,interaction) = ' num2str(p')])
end


%% UNUSED Fig 4x: MSNs < FSIs < TANs activity correlation due to sparseness

% Concatenate multi-aligned trial-averaged activity
mua_exp = vertcat(mua_all{:});

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));

% Align activity by move/outcome (natively stim-aligned)

% (get indicies for alignments - used later)
stim_align = zeros(size(trial_stim_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;
use_align = {stim_align,move_align,outcome_align};

% (move aligned)
move_idx_exp = mat2cell(move_idx,use_split,1);
mua_exp_movealign = vertcat(mua_all{:});
for curr_exp = 1:length(mua_exp_movealign)
   for curr_trial = 1:size(mua_exp_movealign{curr_exp},1)
       mua_exp_movealign{curr_exp}(curr_trial,:,:) = ...
           circshift(mua_exp_movealign{curr_exp}(curr_trial,:,:), ...
           -move_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
   end
end

% (outcome aligned)
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
mua_exp_outcomealign = vertcat(mua_all{:});
for curr_exp = 1:length(mua_exp_outcomealign)
    for curr_trial = 1:size(mua_exp_outcomealign{curr_exp},1)
        mua_exp_outcomealign{curr_exp}(curr_trial,:,:) = ...
            circshift(mua_exp_outcomealign{curr_exp}(curr_trial,:,:), ...
            -outcome_idx_exp{curr_exp}(curr_trial)+leeway_samples,2);
    end
end

% Get average activity across alignments
use_align_labels = {'Stim','Move onset','Outcome'};

act_mean = nan(length(good_units_allcat),length(t),length(use_align_labels));
for curr_align = 1:length(use_align_labels)
    
    switch curr_align
        case 1
            curr_mua = mua_exp;
        case 2
            curr_mua = mua_exp_movealign;
        case 3
            curr_mua = mua_exp_outcomealign;
    end
    
    act_mean(:,:,curr_align) = cell2mat(cellfun(@(act,stim,rxn,outcome) ...
        squeeze(nanmean(act(stim > 0 & rxn < 0.5 & outcome == 1,:,:),1)), ...
        curr_mua,trial_stim_allcat_exp,move_t_exp, ...
        trial_outcome_allcat_exp,'uni',false)')';
    
end

% Concatenate multi-aligned averages
% (split the alignment halfway between median alignment points)
align_median = cellfun(@(x) -nanmedian(x)/sample_rate,use_align);
align_break = align_median(1:end-1) + diff(align_median*0.8);
align_t = {[-0.05,align_break(1)],[align_break(1:2)],[align_break(2),1]};

align_samples = cell(size(align_t));
for curr_align = 1:length(align_samples)    
    curr_t_offset = -nanmedian(use_align{curr_align})/sample_rate;
    curr_t = t + curr_t_offset;
    curr_t_use = curr_t >= align_t{curr_align}(1) & ...
        curr_t <= align_t{curr_align}(2);
    align_samples{curr_align} = curr_t_use;
end

act_mean_multialign = cell2mat(arrayfun(@(x) ...
    act_mean(:,align_samples{x},x),1:size(act_mean,3),'uni',false));


% Get pairwise across-recording activity correlations by celltype
% (this is a super dirty way to do it but I'm in a hurry)
n_recordings = max(recordings_allcat);
celltype_act_corr = cell(3,n_recordings);
celltype_act_minrate = cell(3,n_recordings);
for curr_depth = 1:n_aligned_depths
    for curr_celltype = 1:3
        
        % Get all pairwise activity correlations
        curr_cells = ...
            domain_aligned_allcat == curr_depth & ...
            celltype_allcat == curr_celltype;        
        curr_act_corr = corr(act_mean_multialign(curr_cells,:)');
        
        curr_rate = nanmean(act_mean_multialign(curr_cells,:),2);
        curr_act_minrate = min(curr_rate,curr_rate');
        
        for curr_recording = 1:n_recordings
            
            % Store
            if any(recordings_allcat(curr_cells) == curr_recording)
                curr_act_corr_use = ...
                    curr_act_corr(recordings_allcat(curr_cells) == curr_recording, ...
                    recordings_allcat(curr_cells) ~= curr_recording);
                
                curr_act_minrate_use = ...
                    curr_act_minrate(recordings_allcat(curr_cells) == curr_recording, ...
                    recordings_allcat(curr_cells) ~= curr_recording);
            else
                curr_act_corr_use = NaN;
                curr_act_minrate_use = NaN;
            end
         
            for curr_depth_compare = 1:2
                celltype_act_corr{curr_celltype,curr_recording} = ...
                   cat(1,celltype_act_corr{curr_celltype,curr_recording}, ...
                   curr_act_corr_use(:));
               
               celltype_act_minrate{curr_celltype,curr_recording} = ...
                   cat(1,celltype_act_minrate{curr_celltype,curr_recording}, ...
                   curr_act_minrate_use(:));
            end           
        end
    end
end

figure; hold on;
fr_bins = linspace(0,10,6);%logspace(-2,1,10);
fr_bin_centers = fr_bins(1:end-1)+diff(fr_bins)./2;
set(gca,'ColorOrder',cool(length(fr_bin_centers)));
for curr_fr_bin = 1:length(fr_bin_centers)
    celltype_act_corr_mean_fr_bin = cellfun(@(r,fr) ...
        nanmean(r(fr > fr_bins(curr_fr_bin) & ...
        fr <= fr_bins(curr_fr_bin+1))), ...
        celltype_act_corr,celltype_act_minrate);
    
    errorbar(nanmean(celltype_act_corr_mean_fr_bin,2)', ...
        AP_sem(celltype_act_corr_mean_fr_bin,2),'linewidth',1);
    drawnow;
end
errorbar(nanmean(celltype_act_corr_mean,2)', ...
    AP_sem(celltype_act_corr_mean,2),'k','linewidth',2);
xlim([0.5,3.5]);
set(gca,'XTick',1:3,'XTickLabel',{'MSN','FSI','TAN'});
ylabel('Pairwise activity correlation (across recordings)');
legend([cellfun(@(x) [num2str(x) 'sp/s'], ...
    num2cell(fr_bin_centers),'uni',false),'All rates']);



%% Fig 4e: Cortical maps are similar within domain across celltypes

% Plot average within each celltype/domain
figure;
for curr_depth = 1:n_aligned_depths
    for curr_celltype = 1:n_celltypes
        subplot(n_aligned_depths,n_celltypes, ...
            (curr_depth-1)*n_celltypes+curr_celltype)
        imagesc(nanmean(ctx_str_k_px(:,:, ...
            domain_aligned_allcat(good_units_allcat) == curr_depth & ...
            celltype_allcat(good_units_allcat) == curr_celltype),3));
        axis image
        colormap(brewermap([],'*RdBu'));
        caxis([-max(abs(caxis)),max(abs(caxis))]);
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        axis off
        title(['Str ' num2str(curr_depth) ' ' celltype_labels{curr_celltype}]);
    end
end

% Correlate maps for each celltype within and across domains

% Load kernel templates
n_aligned_depths = 3;
kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template_' num2str(n_aligned_depths) '_depths.mat'];
load(kernel_template_fn);

for curr_depth = 1:n_aligned_depths
    for curr_celltype = 1:3
        
        use_cells = good_units_allcat & domain_aligned_allcat == curr_depth & ...
            ismember(celltype_allcat,[1,2,3]);
        use_maps = ctx_str_k_px(:,:, ...
            domain_aligned_allcat(good_units_allcat) == curr_depth & ...
            ismember(celltype_allcat(good_units_allcat),[1,2,3]));
        
    end
end

% Get map correlations for each celltype within/across domain
% (this is a super dirty way to do it but I'm in a hurry)
n_recordings = max(recordings_allcat);
celltype_map_corr = cell(3,2,n_recordings);
for curr_depth = 1:n_aligned_depths
    for curr_celltype = 1:3
        for curr_recording = 1:n_recordings
            
            use_maps = ctx_str_k_px(:,:, ...
                domain_aligned_allcat(good_units_allcat) == curr_depth & ...
                celltype_allcat(good_units_allcat) == curr_celltype & ...
                recordings_allcat(good_units_allcat) == curr_recording);
            
            if ~isempty(use_maps)
                map_corr = corr(reshape(use_maps,[],size(use_maps,3)), ...
                reshape(kernel_template,[],n_aligned_depths),'type','Pearson');
            else 
               map_corr = nan('single');
            end
    
            [cell_id,depth_group_id] = ndgrid(1:size(map_corr,1),1:size(map_corr,2));
            depth_group = (depth_group_id ~= curr_depth)+1;
            
            for curr_depth_compare = 1:2
                celltype_map_corr{curr_celltype,curr_depth_compare,curr_recording} = ...
                   cat(1,celltype_map_corr{curr_celltype,curr_depth_compare,curr_recording}, ...
                   reshape(map_corr(depth_group == curr_depth_compare),[],1));       
            end
        end
    end
end

celltype_map_corr_mean = cellfun(@nanmean,celltype_map_corr);

figure; hold on;
errorbar(nanmean(celltype_map_corr_mean,3)', ...
    AP_sem(celltype_map_corr_mean,3)','linewidth',2);
xlim([0.5,2.5]);
set(gca,'XTick',1:2,'XTickLabel',{'Within domain','Across domain'});
ylabel('Cell-domain kernel correlation')
legend({'MSN','FSI','TAN'});

% (Map correlation vs. celltype statistics)
disp('Cell-domain map correlation by celltype (2-way anova):');
[celltype_grp,domain_grp,exp_grp] = ndgrid(1:size(celltype_map_corr_mean,1), ...
    1:size(celltype_map_corr_mean,2),1:size(celltype_map_corr_mean,3));
p = anovan(celltype_map_corr_mean(:),...
    [celltype_grp(:),domain_grp(:)],'model','interaction','display','off');
disp(['p(cell type) = ' num2str(p(1))]);
disp(['p(cell type) = ' num2str(p(2))]);
disp(['p(cell type) = ' num2str(p(3))]);


%% Cell type changes with training

mua_exp = vertcat(mua_all{:});
wheel_exp = vertcat(wheel_all{:});
stim_exp = mat2cell(trial_stim_allcat,use_split,1);

% Get average activity with no wheel movement and 100%R stim
wheel_thresh = 0.025;
stim_act = cellfun(@(stim,wheel,act) ...
    permute(nanmean(act(stim == 1 & ...
    ~any(abs(wheel(:,t >= -0.5 & t <= 0.5)) > wheel_thresh,2),:,:),1), ...
    [3,2,1]),stim_exp,wheel_exp,mua_exp,'uni',false);

% Get average activity for each celltype/domain/recording
[stim_act_celltype,group_char] = grpstats(cell2mat(stim_act), ...
    [celltype_allcat,domain_aligned_allcat,recordings_allcat],...
    {'mean','gname'});

group = cellfun(@str2num,group_char);

figure
p = reshape(tight_subplot(n_aligned_depths,n_celltypes),[n_celltypes,n_aligned_depths])';
for curr_depth = 1:n_aligned_depths
    for curr_celltype = 1:n_celltypes
        plot(p(curr_depth,curr_celltype),t, ...
            nanmean(stim_act_celltype( ...
            all(group(:,1:2) == [curr_celltype,curr_depth],2),:),1));
    end
end
linkaxes(get(gcf,'Children'),'xy')






