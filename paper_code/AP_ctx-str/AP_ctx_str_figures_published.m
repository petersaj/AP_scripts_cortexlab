% Generate figures for ctx-str paper 
% (data is prepared in AP_ctx_str_trial_preprocessing)
%
% Code blocks with preceeding symbols: run code block at top with
% corresponding symbol to load data
%
% If no symbol: code is loaded within code block


%% ~~~~~~~~~~~~~ Load data associated with each symbol


%% ++ Depth-aligned striatum

data_fn = 'trial_activity_choiceworld_16strdepth'; % Depth-aligned striatum

AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_recording;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));

%% ** Task (combined: original,ctxephys,pre-muscimol)

data_fn = { ...
    'trial_activity_choiceworld'... 
    'trial_activity_vanillaChoiceworld_ctxstrephys_str'...
    'trial_activity_vanillaChoiceworldNoRepeats_pre_muscimol'};

AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_recording;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));


%% @@ Task single unit/celltypes (combined: original,ctxephys,pre-muscimol)

data_fn = { ...
    'trial_activity_choiceworld_sua.mat', ... % original 
    'trial_activity_ctx_task_sua.mat', ...    % + cortex ephys
    'trial_activity_muscimol_task_sua.mat'};  % muscimol group

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

% (groups are cell type * depth: msn, fsi, tan, uin, narrow tan-like)
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


%% ~~~~~~~~~~~~~ Main figures (some EDF included)


%% Fig 1a-d: Example recording

animal = 'AP025'; 
day = '2017-10-01'; 
experiment = 1; 
str_align = 'none'; 
verbose = true; 
AP_load_experiment;

%%% Plot example data

% Align U's, deconvolve widefield
use_components = 1:200;
aUdf = AP_align_widefield(Udf,animal,day);
fVdf_deconv = AP_deconv_wf(fVdf);
% Set time to plot
plot_t = [15,45];

raster_fig = figure;

% (wheel velocity)
wheel_axes = subplot(6,1,6);
plot_wheel_idx = Timeline.rawDAQTimestamps >= plot_t(1) & ...
    Timeline.rawDAQTimestamps <= plot_t(2);
plot(wheel_axes,Timeline.rawDAQTimestamps(plot_wheel_idx), ...
    wheel_velocity(plot_wheel_idx),'k','linewidth',2);
ylabel('Wheel velocity');
axis off

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
    go_col,'linewidth',2), ...
    find(go_cue_times >= plot_t(1) & go_cue_times <= plot_t(2)));

% (outcomes)
outcome_col = [0,0,0.8;0.5,0.5,0.5];
reward_lines = arrayfun(@(x) line(wheel_axes,repmat(reward_t_timeline(x),1,2),ylim(wheel_axes),'color', ...
    outcome_col(1,:),'linewidth',2), ...
    find(reward_t_timeline >= plot_t(1) & reward_t_timeline <= plot_t(2)));
punish_times = signals_events.responseTimes(trial_outcome == -1);
punish_lines = arrayfun(@(x) line(wheel_axes,repmat(punish_times(x),1,2),ylim(wheel_axes),'color', ...
    outcome_col(2,:),'linewidth',2), ...
    find(punish_times >= plot_t(1) & punish_times <= plot_t(2)));

% (striatum raster)
raster_axes = subplot(6,1,3:5,'YDir','reverse'); hold on;
plot_spikes = spike_times_timeline >= plot_t(1) & ...
    spike_times_timeline <= plot_t(2) & ...
    spike_depths >= str_depth(1) & spike_depths <= str_depth(2);
plot(raster_axes,spike_times_timeline(plot_spikes),spike_depths(plot_spikes),'.k');
ylabel('Depth (\mum)');
xlabel('Time (s)');
depth_scale = 1000;
line(repmat(min(xlim),2,1),[min(ylim),min(ylim) + depth_scale],'color','k','linewidth',3);
axis off

% (fluorescence from select ROIs)
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
roi_trace = AP_svd_roi(aUdf(:,:,use_components),fVdf_deconv(use_components,:),[],[],cat(3,wf_roi.mask));

plot_rois = [1,11,7,17,9,19,10,20];
fluor_spacing = []; % (use default)
fluor_axes = subplot(6,1,1:2); hold on;
plot_fluor_idx = frame_t >= plot_t(1) & frame_t <= plot_t(2);
AP_stackplot(roi_trace(plot_rois,plot_fluor_idx)', ...
    frame_t(plot_fluor_idx),fluor_spacing,[],[0,0.7,0],{wf_roi(plot_rois).area});

y_scale = 0.02;
t_scale = 2;
line([min(xlim),min(xlim) + t_scale],repmat(min(ylim),2,1),'color','k','linewidth',3);
line(repmat(min(xlim),2,1),[min(ylim),min(ylim) + y_scale],'color','k','linewidth',3);
axis off

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

%% ** Fig 1e: Average cortical fluorescence during trial

% Get average stim-aligned fluorescence 
plot_trials = move_t < 0.5 & trial_stim_allcat > 0 & trial_choice_allcat == -1;
plot_trials_exp = mat2cell(plot_trials,use_split,1);

fluor_allcat_deconv_exp = mat2cell(fluor_allcat_deconv,use_split,length(t),n_vs);

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_stim_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set windows to average activity
use_align_labels = {'Stim','Move onset','Outcome'};
use_align = {stim_align,move_align,outcome_align};
plot_t = {[0,0.1],[0,0.1],[0,0.1]};

figure;
p = tight_subplot(1,length(cell2mat(plot_t)));
curr_plot_idx = 0;
for curr_align = 1:length(use_align)
        
    % (re-align activity)
    curr_ctx_act = cellfun(@(act,trials,shift) cell2mat(arrayfun(@(trial) ...
        circshift(act(trial,:,:),shift(trial),2), ...
        find(trials),'uni',false)), ...
        fluor_allcat_deconv_exp,plot_trials_exp, ...
        mat2cell(use_align{curr_align},use_split,1),'uni',false);
    
    curr_ctx_act_mean = ...
        permute(nanmean(cell2mat(cellfun(@(x) nanmean(x,1), ...
        curr_ctx_act,'uni',false)),1),[3,2,1]);
    
    curr_ctx_act_mean_t = interp1(t,curr_ctx_act_mean',plot_t{curr_align})';
    curr_ctx_act_mean_t_px = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
        curr_ctx_act_mean_t);
    
    for curr_plot_t = 1:length(plot_t{curr_align})
        curr_plot_idx = curr_plot_idx+1;
        axes(p(curr_plot_idx));
        imagesc(curr_ctx_act_mean_t_px(:,:,curr_plot_t));
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        axis image off;
        colormap(brewermap([],'Greens'));
        caxis([0,0.01]);
        title([use_align_labels{curr_align} ': ' num2str(plot_t{curr_align}(curr_plot_t)) ' sec']);
        colorbar
    end
    
end

% Matteo-requested: darken colormap to see visual activity better
brighten(-0.75);



%% ++ Fig 1f: Striatum multiunit by depth

% Plot depths present in > 50% of recordings
frac_depths = nanmean(cell2mat(cellfun(@(x) ...
    nanmean(~squeeze(all(all(isnan(x),1),2)),2), ...
    vertcat(mua_all{:})','uni',false)),2);
use_depths = frac_depths > 0.5;

% Plot average stimulus-aligned activity in striatum
plot_trials = move_t < 0.5 & trial_stim_allcat > 0 & trial_choice_allcat == -1;
plot_trials_exp = mat2cell(plot_trials,use_split,1);

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_stim_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Get average activity across alignments
use_align_labels = {'Stim','Move onset','Outcome'};
use_align = {stim_align,move_align,outcome_align};
act_align_avg = nan(n_depths,length(t),length(use_align));

for curr_align = 1:length(use_align)
    
    curr_act_align = cell2mat(arrayfun(@(trial) circshift(mua_allcat(trial,:,:), ...
        use_align{curr_align}(trial),2),transpose(1:size(mua_allcat,1)),'uni',false));
    
    act_align_avg(:,:,curr_align) = ....
        permute(nanmean(cell2mat(cellfun(@(act,trials) ...
        nanmean(act(trials,:,:),1), ...
        mat2cell(curr_act_align,use_split,length(t),n_depths), ...
        plot_trials_exp,'uni',false)),1),[3,2,1]);
    
end

% Plot aligned_activity
align_col = [1,0,0;0.8,0,0.8;0,0,0.8];
% (split the alignment halfway between median alignment points)
align_median = cellfun(@(x) -nanmedian(x(plot_trials))/sample_rate,use_align);
align_break = align_median(1:end-1) + diff(align_median*0.8);
align_t = {[-0.05,align_break(1)],[align_break(1:2)],[align_break(2),1]};

figure; hold on; set(gca,'YDir','reverse');

for curr_align = 1:length(use_align)
    
    curr_t_offset = -nanmedian(use_align{curr_align}(plot_trials))/sample_rate;   
    curr_t = t + curr_t_offset;
    curr_t_plot = curr_t >= align_t{curr_align}(1) & ...
        curr_t <= align_t{curr_align}(2);
    
    curr_act_norm = act_align_avg(:,:,curr_align)./ ...
        max(max(act_align_avg,[],2),[],3);
    plot_t = curr_t > align_t{curr_align}(1) & curr_t <= align_t{curr_align}(2);
    
    imagesc(curr_t(plot_t),[],curr_act_norm(use_depths,plot_t));
    line(repmat(curr_t_offset,2,1),[0.5,sum(use_depths)+0.5],'color',align_col(curr_align,:));
    colormap(gca,brewermap([],'Greys')); 
    
end
xlabel('~Time from stim');
ylabel('Striatum depth');

axis tight;


%% Fig 2a,b,c: Cortex > striatum kernels by depth

% Load kernels by depths
kernel_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
kernel_fn = ['ephys_kernel_depth'];
load([kernel_path filesep kernel_fn])

% (use only trained animals)
naive_animals = {'AP032','AP033','AP034','AP035','AP036'};
use_animals = ~ismember({ephys_kernel_depth.animal},naive_animals);
k_px_cat = horzcat(ephys_kernel_depth(use_animals).k_px)';

% Load, concatenate, mean STA
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
k_fn = [data_path filesep 'ctx_str_kernels_vanillaChoiceworld_16strdepths'];
load(k_fn);
ctx_str_sta_mean = nanmean(cell2mat(permute(horzcat(ctx_str_sta{:}),[1,3,4,2])),4);
ctx_str_sta_mean_norm = ctx_str_sta_mean./max(max(ctx_str_sta_mean,[],1),[],2);

% Pad and concatenate kernels
max_depths = max(cellfun(@(x) size(x,3),k_px_cat));
k_px_cat_pad = cell2mat(permute(cellfun(@(x) padarray(x, ...
    [0,0,max_depths-size(x,3)],NaN,'pre'),k_px_cat,'uni',false),[2,3,4,1]));
k_px_cat_pad_mean = nanmean(k_px_cat_pad,4);
k_px_cat_pad_mean_norm = k_px_cat_pad_mean./max(max(k_px_cat_pad_mean,[],1),[],2);

frac_depths = squeeze(nanmean(all(all(~isnan(k_px_cat_pad),1),2),4));
use_depths = frac_depths > 0.5;

% Plot STA and kernels
if size(k_px_cat_pad_mean_norm,3) ~= size(ctx_str_sta_mean,3)
    error('Different n depths');
end
figure;
subplot(2,1,1);
imagesc(reshape(ctx_str_sta_mean_norm(:,:,use_depths),size(ctx_str_sta_mean_norm,1),[]));
caxis([-1,1]);
colormap(brewermap([],'PRGn'));
axis image off;
subplot(2,1,2);
imagesc(reshape(k_px_cat_pad_mean_norm(:,:,use_depths),size(k_px_cat_pad,1),[]));
caxis([-1,1]);
colormap(brewermap([],'PRGn'));
axis image off;

% Plot STA and kernels (rainbow colormap)
if size(k_px_cat_pad_mean_norm,3) ~= size(ctx_str_sta_mean,3)
    error('Different n depths');
end

use_colormap = flipud(min(jet(sum(use_depths))-0.2,1));

sta_colored = 1-reshape(repmat(ctx_str_sta_mean_norm(:,:,use_depths),1,1,1,3).* ...
    (1-permute(use_colormap,[3,4,1,2])),size(ctx_str_sta_mean_norm,1),[],3);
k_px_colored = 1-reshape(repmat(k_px_cat_pad_mean_norm(:,:,use_depths),1,1,1,3).* ...
    (1-permute(use_colormap,[3,4,1,2])),size(k_px_cat_pad_mean_norm,1),[],3);

figure;
subplot(2,1,1);
image(sta_colored);
caxis([-1,1]);
colormap(brewermap([],'PRGn'));
axis image off;
subplot(2,1,2);
imagesc(k_px_colored);
caxis([-1,1]);
colormap(brewermap([],'PRGn'));
axis image off;


% Plot center-of-mass color
n_aligned_depths = sum(use_depths);
k_px_norm = k_px_cat_pad_mean_norm(:,:,use_depths);

k_px_com = sum(k_px_norm.*permute(1:n_aligned_depths,[1,3,2]),3)./sum(k_px_norm,3);

use_colormap = flipud(min(jet(255)-0.2,1));
k_px_com_colored = ...
    ind2rgb(round(mat2gray(k_px_com,...
    [1,n_aligned_depths])*size(use_colormap,1)),use_colormap);

figure;
p = image(k_px_com_colored);
set(p,'AlphaData',mat2gray(max(k_px_norm,[],3),[0,1]));
axis image off;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);



%% Fig 2d: Ternary plot for map-template correlations

% Load kernel templates for overlay
n_aligned_depths = 3;
kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template_' num2str(n_aligned_depths) '_depths.mat'];
load(kernel_template_fn);

% Load the kernel template matches
n_aligned_depths = 3;
kernel_match_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
kernel_match_fn = ['ephys_kernel_align_' num2str(n_aligned_depths) '_depths.mat'];
load([kernel_match_path filesep kernel_match_fn]);

% (naive animals not included in this plot)
naive_animals = {'AP032','AP033','AP034','AP035','AP036'};
use_animals = ~ismember({ephys_kernel_align.animal},naive_animals);
kernel_corr_cat = cell2mat([ephys_kernel_align(use_animals).kernel_corr]');
kernel_match_cat = cell2mat(horzcat(ephys_kernel_align(use_animals).kernel_match)');

% Get ternary coordinates of map correlations
kernel_corr_cat_offset = kernel_corr_cat+1;
kernel_corr_cat_norm = kernel_corr_cat_offset./sum(kernel_corr_cat_offset,2) - ...
    (kernel_corr_cat_offset./sum(kernel_corr_cat_offset,2))./3;

ternary_tform = [0,1;sqrt(3/4),-0.5;-sqrt(3/4),-0.5];
ternary_limits = (ternary_tform\(eye(3)*max(max(kernel_corr_cat_norm,[],1))))';

kernel_corr_ternary = (ternary_tform\kernel_corr_cat_norm')';

% Plot ternary correlation coorinates and templates
figure;

subplot(n_aligned_depths+1,1,1); hold on;
use_depths = ~isnan(kernel_match_cat);
str_col = max(hsv(n_aligned_depths)-0.2,0);
patch(ternary_limits(:,1),ternary_limits(:,2),'w','linewidth',1);
scatter(kernel_corr_ternary(use_depths,1), ...
    kernel_corr_ternary(use_depths,2), ...
    5,str_col(kernel_match_cat(use_depths),:),'filled');
axis image off;

for curr_template = 1:n_aligned_depths
    subplot(n_aligned_depths+1,1,curr_template+1);
    imagesc(kernel_template(:,:,curr_template));
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    axis image off;
    colormap(brewermap([],'PRGn'));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
end



%% Fig 2e,h: Striatum domain locations and corticostriatal projections

% Set striatum domain colors
str_col = max(hsv(n_aligned_depths)-0.2,0);

% Load the kernel template matches
n_aligned_depths = 3;
kernel_match_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
kernel_match_fn = ['ephys_kernel_align_' num2str(n_aligned_depths) '_depths.mat'];
load([kernel_match_path filesep kernel_match_fn]);

% Pad and concatenate kernel matches
max_n_kernel_match = max(cell2mat(cellfun(@(x) ...
    cellfun(@length,x),{ephys_kernel_align.kernel_match},'uni',false)));

kernel_match_all = cellfun(@(x) cell2mat(cellfun(@(x) ...
    padarray(x,max_n_kernel_match-length(x),NaN,'pre'),x,'uni',false)), ...
    {ephys_kernel_align.kernel_match},'uni',false);

kernel_match_cat = cell2mat(kernel_match_all);

% (use only depths in at least x% of recordings)
frac_depths = nanmean(~isnan(kernel_match_cat),2);
use_depths = frac_depths > 0.5;
kernel_match_cat_use = kernel_match_cat(use_depths,:);

% Get fraction of each domain for each depth
kernel_match_frac = nan(sum(use_depths),n_aligned_depths);
for curr_depth = 1:n_aligned_depths
    curr_kernel_match = +(kernel_match_cat_use == curr_depth);
    curr_kernel_match(isnan(kernel_match_cat_use)) = NaN;
    kernel_match_frac(:,curr_depth) = nanmean(curr_kernel_match,2);
end

% Get center-of-mass for each domain
kernel_match_com = nan(n_aligned_depths,1);
for curr_depth = 1:n_aligned_depths
    curr_kernel_match = +(kernel_match_cat_use == curr_depth);
    curr_kernel_match(isnan(kernel_match_cat_use)) = NaN;
        kernel_match_com(curr_depth) = ...
            nansum(nanmean(curr_kernel_match,2).*(1:sum(use_depths))') ...
            ./nansum(nanmean(curr_kernel_match,2));
end

% Get relative domain COM along trajectory and plot in CCF
kernel_match_com_relative = kernel_match_com./sum(use_depths);

% Define the probe vector manually according to the targeted trajectory
probe_vector_ccf = [520,240,510;520,511,239];

% Load in the annotated Allen volume and names
allen_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
av = readNPY([allen_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean

% Get probe location per micron
probe_size = pdist2(probe_vector_ccf(1,:),probe_vector_ccf(2,:))*10;
probe_depths = ...
    round([linspace(probe_vector_ccf(1,1)',probe_vector_ccf(2,1)',probe_size); ...
    linspace(probe_vector_ccf(1,2)',probe_vector_ccf(2,2)',probe_size); ...
    linspace(probe_vector_ccf(1,3)',probe_vector_ccf(2,3)',probe_size)]');

% Eliminiate trajectory points that are off the atlas
eliminate_depths = ...
    probe_depths(:,1) < 1 | probe_depths(:,1) > size(av,1) | ...
    probe_depths(:,2) < 1 | probe_depths(:,2) > size(av,2) | ...
    probe_depths(:,3) < 1 | probe_depths(:,3) > size(av,3);
probe_depths(eliminate_depths,:) = [];

% Convert probe depths subscripts to indicies
probe_depths_ind = sub2ind(size(av),probe_depths(:,1),probe_depths(:,2),probe_depths(:,3));

% Get structures that the probe is in
probe_structures = av(probe_depths_ind);

% Get target relative depths through striatum
str_id = find(strcmp(st.safe_name,'Caudoputamen'));
probe_structures_str = probe_structures == str_id;
probe_str = [probe_depths(find(probe_structures_str,1,'first'),:); ...
    probe_depths(find(probe_structures_str,1,'last'),:)];
kernel_depth_ccf = interp1([0,1],probe_str,kernel_match_com_relative);

% Plot brain to overlay probes
% (note the CCF is rotated to allow for dim 1 = x)
h = figure; ccf_axes = axes; hold on
slice_spacing = 10;
target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) > 1,[2,1,3]);
structure_patch = isosurface(target_volume,0);
structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
target_structure_color = [0.7,0.7,0.7];
brain_outline = patch('Vertices',structure_wire.vertices*slice_spacing, ...
    'Faces',structure_wire.faces, ...
    'FaceColor','none','EdgeColor',target_structure_color);

str_id = find(strcmp(st.safe_name,'Caudoputamen'));
target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) == str_id,[2,1,3]);
structure_patch = isosurface(target_volume,0);
structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
target_structure_color = [0.7,0,0.7];
striatum_outline = patch('Vertices',structure_wire.vertices*slice_spacing, ...
    'Faces',structure_wire.faces, ...
    'FaceColor','none','EdgeColor',target_structure_color);

axis image vis3d off;
view([-30,25]);
cameratoolbar(h,'SetCoordSys','y');
cameratoolbar(h,'SetMode','orbit');

scatter3(kernel_depth_ccf(:,1),kernel_depth_ccf(:,2),kernel_depth_ccf(:,3), ...
    100,str_col,'filled');

% Plot brain area outlines at slice
callosum_id = find(strcmp(st.safe_name,'corpus callosum body'));
ventricle_id = find(strcmp(st.safe_name,'lateral ventricle'));

av_slice = permute(av(probe_vector_ccf(1),:,:),[2,3,1]);
slice_brain_outline = bwboundaries(av_slice > 1,'noholes');
slice_str_outline = bwboundaries(av_slice == str_id,'noholes');
slice_callosum_outline = bwboundaries(av_slice == callosum_id,'noholes');
slice_ventricle_outline = bwboundaries(av_slice == ventricle_id,'noholes');

% Plot mean domain CCF locations
figure; 
subplot(1,2,1);
hold on; axis image on; box on; grid on; set(gca,'YDir','reverse');
plot(slice_brain_outline{1}(:,2),slice_brain_outline{1}(:,1),'k','linewidth',2);
cellfun(@(x) plot(x(:,2),x(:,1),'b','linewidth',2),slice_str_outline);
cellfun(@(x) fill(x(:,2),x(:,1),'k','linewidth',2),slice_callosum_outline);
cellfun(@(x) fill(x(:,2),x(:,1),'k','linewidth',2),slice_ventricle_outline);

line(probe_vector_ccf(:,3),probe_vector_ccf(:,2),'linewidth',2,'color','r');
scatter(kernel_depth_ccf(:,3),kernel_depth_ccf(:,2), ...
    100,str_col,'filled');

% Plot fraction domains by depth
subplot(1,2,2); hold on;
colormap(str_col);
area(kernel_match_frac,'FaceColor','flat');
for curr_depth = 1:n_aligned_depths
    line(repmat(kernel_match_com(curr_depth),2,1),ylim, ...
        'color',min(str_col(curr_depth,:)+0.2,1),'linewidth',3);
end
axis tight;
xlabel('Estimated depth');
ylabel('Fraction domain match');

% Mirror all of the locations across the midline and get projections from
% both (because the cortical injections in the database aren't evenly
% distributed, so this is more comprehensive)
kernel_depth_um = round(kernel_depth_ccf*10);
bregma_um = allenCCFbregma*10;
ccf_midline = bregma_um(3);
hemisphere = sign(kernel_depth_um(1,3) - ccf_midline);
str_depths_mirror = [kernel_depth_um(:,1:2),ccf_midline - hemisphere*abs(kernel_depth_um(:,3)-ccf_midline)];

str_depths_query = [kernel_depth_um;str_depths_mirror];
injection_parameters = get_allen_projection(str_depths_query);
injection_coordinates = {injection_parameters.coordinates};

% Standardize injection coordinates by hemisphere (left = contra, right =
% ipsi)
injection_coordinates_standardized = injection_coordinates;
for curr_coord = 1:length(injection_coordinates)
    
    target_hemisphere = sign(ccf_midline - str_depths_query(curr_coord,3));
    injection_coords_ml_offset = abs(injection_coordinates{curr_coord}(:,3) - ccf_midline);
    injection_coordinates_hemisphere = sign(injection_coordinates{curr_coord}(:,3) - ccf_midline);
    
    injection_coords_ipsi = injection_coordinates_hemisphere == target_hemisphere;
    injection_coords_contra = injection_coordinates_hemisphere == -target_hemisphere;
    
    injection_coordinates_standardized{curr_coord}(injection_coords_ipsi,3) = ...
        ccf_midline + injection_coords_ml_offset(injection_coords_ipsi);
    injection_coordinates_standardized{curr_coord}(injection_coords_contra,3) = ...
        ccf_midline - injection_coords_ml_offset(injection_coords_contra);
    
end

% Combine side-standardized injection data across hemispheres
injection_coordinates_bilateral = arrayfun(@(x) ...
    vertcat(injection_coordinates_standardized{x}, ...
    injection_coordinates_standardized{n_aligned_depths+x}), ...
    1:n_aligned_depths,'uni',false);

% Convert points from CCF to widefield
alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment';
ccf_tform_fn = [alignment_path filesep 'ccf_tform.mat'];
load(ccf_tform_fn);

um2pixel = 20.6;
injection_coordinates_bilateral_wf = cellfun(@(x) ...
    [x(:,[3,1]).*(1/(um2pixel)),ones(size(x,1),1)]*ccf_tform.T, ...
    injection_coordinates_bilateral,'uni',false);

% Get projection density
projection_strength = {injection_parameters.density};

projection_strength_bilateral = arrayfun(@(x) ...
    vertcat(projection_strength{x}', ...
    projection_strength{n_aligned_depths+x}'), ...
    1:n_aligned_depths,'uni',false);

% Load kernel templates for overlay
n_aligned_depths = 3;
kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template_' num2str(n_aligned_depths) '_depths.mat'];
load(kernel_template_fn);

figure; 
colormap(brewermap([],'PRGn'));
for curr_depth = 1:n_aligned_depths
    subplot(n_aligned_depths,1,curr_depth); hold on; axis image off;
    set(gca,'YDir','reverse');
    
    scatter(injection_coordinates_bilateral_wf{curr_depth}(:,1), ...
        injection_coordinates_bilateral_wf{curr_depth}(:,2), ...
        projection_strength_bilateral{curr_depth}*50 + 10, ...
        'k','filled');
    
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
end

% Plot heatmap of injection sites
figure;
colormap(brewermap([],'Greys'));
for curr_depth = 1:n_aligned_depths
    subplot(n_aligned_depths,1,curr_depth); hold on; axis image off;
    set(gca,'YDir','reverse');
    
    bin_width = 10;
    x_bins = 0:bin_width:500;
    y_bins = 0:bin_width:500;
    [bin_n,~,~,bin_x,bin_y] = histcounts2( ...
        injection_coordinates_bilateral_wf{curr_depth}(:,1), ...
        injection_coordinates_bilateral_wf{curr_depth}(:,2),x_bins,y_bins);
    
    x_bin_centers = x_bins(1:end-1) + diff(x_bins)./2;
    y_bin_centers = y_bins(1:end-1) + diff(x_bins)./2;
    
    % (tried this - not used)
%     projection_strength_bin = accumarray( ...
%         [bin_x,bin_y], ...
%         projection_strength_bilateral{curr_depth},[length(x_bin_centers), ...
%         length(y_bin_centers)],[],0);
    
    imagesc(y_bin_centers,x_bin_centers,imgaussfilt(bin_n',2));
        
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
end



%% Fig 2f,g,i,j: Average cortex > striatum domain kernels

% Load data
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\str_ctxpred.mat')

% Load Master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Get time
framerate = 35;
upsample_factor = 1;
sample_rate = framerate*upsample_factor;
kernel_t = [-0.1,0.1];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
t = kernel_frames/sample_rate;

% Concatenate kernels and convert to pixels
% (and flip in time so it's fluorescence lead:lag spikes)
ctx_str_k_animal = [str_ctxpred.ctx_str_k]';
ctx_str_k_px_animal = cellfun(@(x) cellfun(@(x) ...
    flip(AP_svdFrameReconstruct(U_master(:,:,1:100),x),3),x,'uni',false), ...
    ctx_str_k_animal,'uni',false);


% Get mean kernels and plot
ctx_str_k_px_cat = vertcat(ctx_str_k_px_animal{:}); 

ctx_str_k_px_task_mean = nanmean(cat(5,ctx_str_k_px_cat{:,1}),5);
ctx_str_k_px_notask_mean = nanmean(cat(5,ctx_str_k_px_cat{:,2}),5);

n_depths = size(ctx_str_k_px_task_mean,4);

AP_image_scroll([ctx_str_k_px_task_mean,ctx_str_k_px_notask_mean]);
axis image;
colormap(brewermap([],'PRGn'));
caxis([-max(abs(caxis)),max(abs(caxis))]);
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5],[],[size(U_master,1),size(U_master,2),1,2]);

figure;
for curr_depth = 1:n_depths
    subplot(n_depths,2,(curr_depth-1)*2+1);
    imagesc(ctx_str_k_px_task_mean(:,:,t == 0,curr_depth));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'PRGn'));
    axis image off;
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title('Task');
    
    subplot(n_depths,2,(curr_depth-1)*2+2);
    imagesc(ctx_str_k_px_notask_mean(:,:,t == 0,curr_depth));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'PRGn'));
    axis image off;
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title('Passive');
end

figure;
for curr_depth = 1:n_depths
    subplot(n_depths,1,curr_depth);
    imagesc(reshape(ctx_str_k_px_task_mean(:,:,:,curr_depth), ...
        size(ctx_str_k_px_task_mean,1),[]));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'PRGn'));
    axis image off;
    title('Task');
end

figure;
for curr_depth = 1:n_depths
    subplot(n_depths,1,curr_depth);
    imagesc(reshape(ctx_str_k_px_notask_mean(:,:,:,curr_depth), ...
        size(ctx_str_k_px_notask_mean,1),[]));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'PRGn'));
    axis image off;
    title('Passive');
end


% Get within- across-condition kernel correlation
corr_type = 'Pearson';
task_notask_k_corr = nan(4,n_depths,length(ctx_str_k_px_animal));
for curr_animal = 1:length(ctx_str_k_px_animal)
    
    curr_px = cellfun(@(x) reshape(x,[],n_depths), ...
        ctx_str_k_px_animal{curr_animal},'uni',false);
    
    curr_px_task = cat(3,curr_px{:,1});
    curr_px_notask = cat(3,curr_px{:,2});
     
    % Correlate kernel task/notask within domain
    task_notask_k_corr(1,:,curr_animal) = ...
        nanmean(cell2mat(cellfun(@(x,y) diag(corr(x,y,'type',corr_type))', ...
        curr_px(:,1),curr_px(:,2),'uni',false)));
    
    % Correlate kernel task across domains
    task_notask_k_corr(2,:,curr_animal) = ...
        nanmean(cell2mat(cellfun(@(x,y) ...
        nansum(tril(corr(x,'type',corr_type),-1)+triu(corr(x,'type',corr_type),1),1)./(n_depths-1)', ...
        curr_px(:,1),'uni',false)));
       
    % Correlate kernel within task/notask across days within domain
    task_notask_k_corr(3,:,curr_animal) = arrayfun(@(depth) ...
        nanmean(AP_itril(corr(permute(curr_px_task(:,depth,:),[1,3,2]),'type',corr_type),-1)),1:n_depths);
    task_notask_k_corr(4,:,curr_animal) = arrayfun(@(depth) ...
        nanmean(AP_itril(corr(permute(curr_px_notask(:,depth,:),[1,3,2]),'type',corr_type),-1)),1:n_depths);  

end

% Get mean across domains
task_notask_k_corr_strmean = squeeze(nanmean(task_notask_k_corr,2));

% Plot mean and split by domains
str_col = max(hsv(n_depths)-0.2,0);

figure; 

subplot(2,1,1);hold on; set(gca,'ColorOrder',str_col);
plot(task_notask_k_corr_strmean,'color',[0.5,0.5,0.5]);
errorbar(nanmean(task_notask_k_corr_strmean,2), ...
    AP_sem(task_notask_k_corr_strmean,2),'k','linewidth',2);
set(gca,'XTick',1:4,'XTickLabelRotation',20,'XTickLabel', ...
    {'Task-no task within day','Task within day across domains','Task across days','No task across days'})
ylabel('Spatiotemporal correlation');
xlim([0.5,4.5]);

subplot(2,1,2);hold on; set(gca,'ColorOrder',str_col);
errorbar(nanmean(task_notask_k_corr,3), ...
    AP_sem(task_notask_k_corr,3),'linewidth',2)
set(gca,'XTick',1:4,'XTickLabelRotation',20,'XTickLabel', ...
    {'Task-no task within day','Task within day across domains','Task across days','No task across days'})
ylabel('Spatiotemporal correlation');
xlim([0.5,4.5]);
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false))

% (within task-passive v task-task domains statistics)
disp('Task/passive vs task/task cross-domain:')
curr_p = signrank(squeeze(task_notask_k_corr_strmean(1,:)), ...
    squeeze(task_notask_k_corr_strmean(2,:)));
disp(['All str p = ' num2str(curr_p)]);
for curr_depth = 1:n_depths  
    curr_p = signrank(squeeze(task_notask_k_corr(1,curr_depth,:)), ...
        squeeze(task_notask_k_corr(2,curr_depth,:)));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end

% (within vs across statistics)
disp('Task/passive-within vs task-across:')
curr_p = signrank(squeeze(task_notask_k_corr_strmean(1,:)), ...
    squeeze(task_notask_k_corr_strmean(3,:)));
disp(['All str p = ' num2str(curr_p)]);
for curr_depth = 1:n_depths
    curr_p = signrank(squeeze(task_notask_k_corr(1,curr_depth,:)), ...
        squeeze(task_notask_k_corr(3,curr_depth,:)));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end

% (cross task vs no task statistics)
disp('Task-across vs passive-across');
curr_p = signrank(squeeze(task_notask_k_corr_strmean(3,:)), ...
    squeeze(task_notask_k_corr_strmean(4,:)));
disp(['All str p = ' num2str(curr_p)]);
for curr_depth = 1:n_depths
    curr_p = signrank(squeeze(task_notask_k_corr(3,curr_depth,:)), ...
        squeeze(task_notask_k_corr(4,curr_depth,:)));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end


% Total weight sum in time
ctx_str_k_px_sum = cellfun(@(x) ...
    permute(sum(abs(reshape(x,[],size(x,3),size(x,4))),1),[2,3,1]), ...
    ctx_str_k_px_cat,'uni',false);
% (min-max norm)
ctx_str_k_px_sum_norm = cellfun(@(x) ...
    (x-min(x,[],1))./max(x-min(x,[],1)), ...
    ctx_str_k_px_sum,'uni',false);

% (average across domain)
ctx_str_k_px_sum_norm_mean = permute(cell2mat(permute(cellfun(@(x) ...
    nanmean(x,2),ctx_str_k_px_sum_norm,'uni',false),[3,1,2])),[1,3,2]);
figure;
AP_errorfill(t,nanmean(ctx_str_k_px_sum_norm_mean,3), ...
    AP_sem(ctx_str_k_px_sum_norm_mean,3),[0,0,0;1,0,0]);
xlabel('Ctx lead:lag str (s)');
ylabel('Sum(abs(W))');
line([0,0],ylim,'color','k');

% (for each domain)
figure;
task_sum_norm_tdiff = nan(n_depths,size(ctx_str_k_px_sum,1));
for curr_depth = 1:n_depths
    subplot(n_depths,1,curr_depth); hold on;
    
    % Min/max-normalized sum of weights    
    curr_task_sum_norm = cell2mat(cellfun(@(x) ...
        x(:,curr_depth), ...
        ctx_str_k_px_sum_norm(:,1)','uni',false));
    curr_notask_sum_norm = cell2mat(cellfun(@(x) ...
        x(:,curr_depth), ...
        ctx_str_k_px_sum_norm(:,2)','uni',false));
    
    % Get weight difference t < 0 and t > 0
    task_sum_norm_tdiff(curr_depth,:) = ...
        nanmean(curr_task_sum_norm(t < 0,:) - ...
        flipud(curr_task_sum_norm(t > 0,:)),1);
    
    AP_errorfill(t,nanmean(curr_task_sum_norm,2),AP_sem(curr_task_sum_norm,2),'k');
    AP_errorfill(t,nanmean(curr_notask_sum_norm,2),AP_sem(curr_notask_sum_norm,2),'r');
    xlabel('Ctx lead:lag str (s)');
    ylabel('Sum(abs(W))');
    line([0,0],ylim,'color','k');
    
end

% (cross task vs no task statistics)
disp('Weight t<0 vs t>0:')
curr_p = signrank(nanmean(task_sum_norm_tdiff,1));
disp(['All str  p = ' num2str(curr_p)]);

for curr_depth = 1:n_depths
    curr_p = signrank(task_sum_norm_tdiff(curr_depth,:));
    disp(['Str ' num2str(curr_depth) ' weight t <0 vs > 0 = ' num2str(curr_p)]); 
end



%%%% MAYBE DON'T DO THE BELOW
k_px = ctx_str_k_px_task_mean;
k_px = k_px./max(max(k_px,[],1),[],2);


% Plot center-of-mass color at select time points
k_px_com = sum(k_px.*permute(1:n_depths,[1,3,4,2]),4)./sum(k_px,4);
k_px_com_colored = nan(size(k_px_com,1),size(k_px_com,2),3,size(k_px_com,3));

use_colormap = min(jet(255)-0.2,1);
for curr_frame = 1:size(k_px_com,3)
    k_px_com_colored(:,:,:,curr_frame) = ...
        ind2rgb(round(mat2gray(k_px_com(:,:,curr_frame),...
        [1,n_depths])*size(use_colormap,1)),use_colormap);
end

plot_t = [-0.05:0.025:0.05];
k_px_com_colored_t = ...
    permute(reshape(interp1(t,permute(reshape(k_px_com_colored,[],3,length(t)), ...
    [3,1,2]),plot_t),length(plot_t),size(k_px_com_colored,1), ...
    size(k_px_com_colored,2),3),[2,3,4,1]);

k_px_max = squeeze(max(k_px,[],4));
k_px_max_t = ...
    permute(reshape(interp1(t,reshape(k_px_max,[],length(t))', ...
    plot_t),length(plot_t),size(k_px_max,1), ...
    size(k_px_max,2)),[2,3,1]);

weight_max = 1;
figure;
for t_idx = 1:length(plot_t)
    subplot(1,length(plot_t),t_idx);
    p = image(k_px_com_colored_t(:,:,:,t_idx));
    set(p,'AlphaData', ...
        mat2gray(k_px_max_t(:,:,t_idx),[0,weight_max]));
    axis image off;
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title([num2str(plot_t(t_idx)),' s']);
end


%% ** Fig 3a-e, EDFig5,7a-c: Striatal domain activity and task regression

% Plot stim-aligned/sorted measured and predicted striatum activity
% (correct contra trials)
for curr_trial_set = 1:2
    switch curr_trial_set
        case 1
            plot_trials = move_t < Inf & trial_stim_allcat > 0 & trial_choice_allcat == -1;
            figure('Name','Correct contra trials');
        case 2
            plot_trials = move_t < Inf & trial_stim_allcat < 0 & trial_choice_allcat == 1;
            figure('Name','Correct ipsi trials');
    end
    
    p = gobjects(n_depths,4);
    colormap(brewermap([],'Greys'));
    for curr_depth = 1:n_depths
        
        % Get trials to plot, sort by reaction time
        curr_trials = plot_trials & ~all(isnan(mua_allcat(:,:,curr_depth)),2);
        curr_trials_exp = mat2cell(curr_trials,use_split,1);
        curr_trials_idx = find(curr_trials);
        % (sort by reaction time)
        [~,rxn_sort_idx] = sort(move_t(curr_trials_idx));
        % (sort by reaction time and contrast)
%         [~,rxn_sort_idx] = sortrows([trial_stim_allcat(curr_trials),move_t(curr_trials)]);
        
        sorted_plot_trials = curr_trials_idx(rxn_sort_idx);
        
        curr_plot = mua_allcat(sorted_plot_trials,:,curr_depth);
        curr_taskpred_plot = mua_taskpred_allcat(sorted_plot_trials,:,curr_depth);
        curr_ctxpred_plot = mua_ctxpred_allcat(sorted_plot_trials,:,curr_depth);
        
        % Smooth and plot with stim/move/reward times
        % (as conv(nans-zeroed)./conv(non-nan) to ignore in nans in conv)
        smooth_filt = [100,1]; % (trials x frames)
        
        curr_plot_smooth = conv2(curr_plot,ones(smooth_filt),'same')./ ...
            conv2(~isnan(curr_plot),ones(smooth_filt),'same');
        
        curr_taskpred_plot_smooth = curr_taskpred_plot;
        curr_taskpred_plot_smooth(isnan(curr_taskpred_plot_smooth)) = 0;
        curr_taskpred_plot_smooth = conv2(curr_taskpred_plot_smooth,ones(smooth_filt),'same')./ ...
            conv2(~isnan(curr_taskpred_plot),ones(smooth_filt),'same');
        
        curr_ctxpred_plot_smooth = curr_ctxpred_plot;
        curr_ctxpred_plot_smooth(isnan(curr_ctxpred_plot_smooth)) = 0;
        curr_ctxpred_plot_smooth = conv2(curr_ctxpred_plot_smooth,ones(smooth_filt),'same')./ ...
            conv2(~isnan(curr_ctxpred_plot),ones(smooth_filt),'same');
        
        p(curr_depth,1) = subplot(n_depths,4,1+(curr_depth-1)*4,'YDir','reverse'); hold on
        imagesc(t,[],curr_plot_smooth);
        plot(zeros(size(curr_trials_idx)),1:length(curr_trials_idx),'MarkerSize',1,'color','r');
        plot(move_t(sorted_plot_trials),1:length(curr_trials_idx),'MarkerSize',1,'color',[0.8,0,0.8]);
        %         plot(outcome_t(sorted_plot_trials),1:length(curr_trials_idx),'.','MarkerSize',1,'color','b');
        axis tight;
        xlim([-0.2,1]);
        xlabel('Time from stim');
        ylabel('Trials (rxn sorted)');
        title('Measured');
        
        p(curr_depth,2) = subplot(n_depths,4,2+(curr_depth-1)*4,'YDir','reverse'); hold on
        imagesc(t,[],curr_taskpred_plot_smooth);
        plot(zeros(size(curr_trials_idx)),1:length(curr_trials_idx),'MarkerSize',1,'color','r');
        plot(move_t(sorted_plot_trials),1:length(curr_trials_idx),'MarkerSize',1,'color',[0.8,0,0.8]);
        %         plot(outcome_t(sorted_plot_trials),1:length(curr_trials_idx),'.','MarkerSize',1,'color','b');
        axis tight;
        xlim([-0.2,1]);
        xlabel('Time from stim');
        ylabel('Trials (rxn sorted)');
        title('Task-predicted');
        
        p(curr_depth,3) = subplot(n_depths,4,3+(curr_depth-1)*4,'YDir','reverse'); hold on
        imagesc(t,[],curr_ctxpred_plot_smooth);
        plot(zeros(size(curr_trials_idx)),1:length(curr_trials_idx),'MarkerSize',1,'color','r');
        plot(move_t(sorted_plot_trials),1:length(curr_trials_idx),'MarkerSize',1,'color',[0.8,0,0.8]);
        %         plot(outcome_t(sorted_plot_trials),1:length(curr_trials_idx),'.','MarkerSize',1,'color','b');
        axis tight;
        xlim([-0.2,1]);
        xlabel('Time from stim');
        ylabel('Trials (rxn sorted)');
        title('Cortex-predicted');
        
        % Plot average aligned activity
        % (set alignment shifts)
        t_leeway = -t(1);
        leeway_samples = round(t_leeway*(sample_rate));
        stim_align = zeros(size(trial_stim_allcat));
        move_align = -move_idx + leeway_samples;
        outcome_align = -outcome_idx + leeway_samples;
        use_align = {stim_align,move_align,outcome_align};
        
        align_col = [1,0,0;0.8,0,0.8;0,0,0.8];
        % (split the alignment halfway between median alignment points)
        align_median = cellfun(@(x) -nanmedian(x(plot_trials))/sample_rate,use_align);
        align_break = align_median(1:end-1) + diff(align_median*0.8);
        align_t = {[-0.05,align_break(1)],[align_break(1:2)],[align_break(2),1]};
        
        p(curr_depth,4) = subplot(n_depths,4,4+(curr_depth-1)*4); hold on
        for curr_align = 1:length(use_align)
            curr_mua_align = cell2mat(arrayfun(@(trial) circshift(mua_allcat(trial,:,:), ...
                use_align{curr_align}(trial),2),transpose(1:size(mua_allcat,1)),'uni',false));
            curr_mua_exp = mat2cell(curr_mua_align(:,:,curr_depth),use_split,length(t));
            curr_mua_exp_mean = cell2mat(cellfun(@(data,trials) nanmean(data(trials,:),1), ...
                curr_mua_exp,curr_trials_exp,'uni',false));
            
            curr_mua_taskpred_align = cell2mat(arrayfun(@(trial) circshift(mua_taskpred_allcat(trial,:,:), ...
                use_align{curr_align}(trial),2),transpose(1:size(mua_taskpred_allcat,1)),'uni',false));
            curr_mua_taskpred_exp = mat2cell(curr_mua_taskpred_align(:,:,curr_depth),use_split,length(t));
            curr_mua_taskpred_exp_mean = cell2mat(cellfun(@(data,trials) nanmean(data(trials,:),1), ...
                curr_mua_taskpred_exp,curr_trials_exp,'uni',false));
            
            curr_mua_ctxpred_align = cell2mat(arrayfun(@(trial) circshift(mua_ctxpred_allcat(trial,:,:), ...
                use_align{curr_align}(trial),2),transpose(1:size(mua_ctxpred_allcat,1)),'uni',false));
            curr_mua_ctxpred_exp = mat2cell(curr_mua_ctxpred_align(:,:,curr_depth),use_split,length(t));
            curr_mua_ctxpred_exp_mean = cell2mat(cellfun(@(data,trials) nanmean(data(trials,:),1), ...
                curr_mua_ctxpred_exp,curr_trials_exp,'uni',false));
            
            curr_t_offset = -nanmedian(use_align{curr_align}(plot_trials))/sample_rate;
            curr_t = t + curr_t_offset;
            curr_t_plot = curr_t >= align_t{curr_align}(1) & ...
                curr_t <= align_t{curr_align}(2);
            
            plot_t = curr_t > align_t{curr_align}(1) & curr_t <= align_t{curr_align}(2);
            
            AP_errorfill(curr_t(plot_t'), ...
                nanmean(curr_mua_exp_mean(:,plot_t),1)', ...
                AP_sem(curr_mua_exp_mean(:,plot_t),1)','k');
            
            AP_errorfill(curr_t(plot_t'), ...
                nanmean(curr_mua_taskpred_exp_mean(:,plot_t),1)', ...
                AP_sem(curr_mua_taskpred_exp_mean(:,plot_t),1)','b');
            
            AP_errorfill(curr_t(plot_t'), ...
                nanmean(curr_mua_ctxpred_exp_mean(:,plot_t),1)', ...
                AP_sem(curr_mua_ctxpred_exp_mean(:,plot_t),1)',[0,0.7,0]);
            
            line(repmat(curr_t_offset,2,1),ylim,'color',align_col(curr_align,:));
        end
        xlabel('~Time from stim');
        ylabel('Striatum depth');
        
    end
    % Link the x-axes, set the c/y-axes same within a row
    linkaxes(p(:),'x');
    
    for curr_row = 1:size(p,1)
        curr_ylim = ylim(p(curr_row,4));
        caxis(p(curr_row,1),[0,curr_ylim(2)]);
        caxis(p(curr_row,2),[0,curr_ylim(2)]);
        caxis(p(curr_row,3),[0,curr_ylim(2)]);
    end
    
    trial_scale = 500;
    t_scale = 0.5;
    y_scale = 1;
    line(p(1,1),min(xlim(p(1,1))) + [0,t_scale],repmat(min(ylim(p(1,1))),2,1),'color','k','linewidth',3);
    line(p(1,4),min(xlim(p(1,4))) + [0,t_scale],repmat(min(ylim(p(1,4))),2,1),'color','k','linewidth',3);
    line(p(1,1),repmat(min(xlim(p(1,1))),2,1),min(ylim(p(1,1))) + [0,trial_scale],'color','k','linewidth',3);
    line(p(1,4),repmat(min(xlim(p(1,4))),2,1),min(ylim(p(1,4))) + [0,y_scale],'color','k','linewidth',3);
    
end

% Get task>striatum parameters
n_regressors = length(task_regressor_labels);

% Normalize task > striatum kernels across experiments with mua_norm
mua_taskpred_k_allcat_norm = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    vertcat(mua_taskpred_k_all{:}),vertcat(mua_norm{:}),'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

mua_ctxpred_taskpred_k_allcat_norm = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    vertcat(mua_ctxpred_taskpred_k_all{:}),vertcat(mua_norm{:}),'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

% Plot task>striatum kernels
stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = [];
move_col = [0.6,0,0.6;1,0.6,0];
go_col = [0,0,0;0.5,0.5,0.5];
outcome_col = [0.2,0.8,1;0,0,0];
task_regressor_cols = {stim_col,move_col,go_col,outcome_col};
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

figure('Name','Task > Striatum');
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors  
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);  
        
        curr_kernels = mua_taskpred_k_allcat_norm{curr_regressor}(:,:,curr_depth,:);
        n_subregressors = size(mua_taskpred_k_allcat_norm{curr_regressor},1);
        col = task_regressor_cols{curr_regressor};
        for curr_subregressor = 1:n_subregressors
            AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
                nanmean(curr_kernels(curr_subregressor,:,:,:),4)', ...
                AP_sem(curr_kernels(curr_subregressor,:,:,:),4)', ...
                col(curr_subregressor,:),0.5);
        end
        
        xlabel('Time (s)');
        ylabel('Weight');
        title(task_regressor_labels{curr_regressor});
        line([0,0],ylim,'color','k');
        
    end
end
linkaxes(p);
y_scale = 1;
t_scale = 0.5;
line(min(xlim) + [0,t_scale],repmat(min(ylim),2,1),'color','k','linewidth',3);
line(repmat(min(xlim),2,1),min(ylim) + [0,y_scale],'color','k','linewidth',3);

figure('Name','Task > Cortex (str)');
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors  
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);  
        
        curr_kernels = mua_ctxpred_taskpred_k_allcat_norm{curr_regressor}(:,:,curr_depth,:);
        n_subregressors = size(mua_ctxpred_taskpred_k_allcat_norm{curr_regressor},1);
        col = task_regressor_cols{curr_regressor};
        for curr_subregressor = 1:n_subregressors
            AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
                nanmean(curr_kernels(curr_subregressor,:,:,:),4)', ...
                AP_sem(curr_kernels(curr_subregressor,:,:,:),4)', ...
                col(curr_subregressor,:),0.5);
        end
        
        xlabel('Time (s)');
        ylabel('Weight');
        title(task_regressor_labels{curr_regressor});
        line([0,0],ylim,'color','k');
        
    end
end
linkaxes(p);
y_scale = 1;
t_scale = 0.5;
line(min(xlim) + [0,t_scale],repmat(min(ylim),2,1),'color','k','linewidth',3);
line(repmat(min(xlim),2,1),min(ylim) + [0,y_scale],'color','k','linewidth',3);


% Task kernel str/ctx correlation
mua_taskpred_catk = cellfun(@(x) cellfun(@(x) ...
    cell2mat(cellfun(@(x) reshape(x,[],size(x,3)),x,'uni',false)), ...
    x,'uni',false),mua_taskpred_k_all,'uni',false);

mua_ctxpred_taskpred_catk = cellfun(@(x) cellfun(@(x) ...
    cell2mat(cellfun(@(x) reshape(x,[],size(x,3)),x,'uni',false)), ...
    x,'uni',false),mua_ctxpred_taskpred_k_all,'uni',false);

ctx_str_taskk_animal = cellfun(@(x,y) [x,y], ...
    mua_taskpred_catk,mua_ctxpred_taskpred_catk,'uni',false);

ctx_str_taskk_corr = nan(4,n_depths,length(ctx_str_taskk_animal));
for curr_animal = 1:length(ctx_str_taskk_animal)
    
    curr_k = cellfun(@(x) reshape(x,[],n_depths), ...
        ctx_str_taskk_animal{curr_animal},'uni',false);
    
    curr_k_str = cat(3,curr_k{:,1});
    curr_k_ctx = cat(3,curr_k{:,2});
     
    % Correlate str/ctx kernels within domain
    ctx_str_taskk_corr(1,:,curr_animal) = ...
        nanmean(cell2mat(cellfun(@(x,y) diag(corr(x,y))', ...
        curr_k(:,1),curr_k(:,2),'uni',false)));
    
    % Correlate str/ctx kernels across domains
    ctx_str_taskk_corr(2,:,curr_animal) = ...
        nanmean(cell2mat(cellfun(@(x,y) ...
        nansum(tril(corr(x),-1)+triu(corr(x),1),1)./(n_depths-1)', ...
        curr_k(:,1),'uni',false)));
       
    % Correlate kernel within task/notask across days within domain
    ctx_str_taskk_corr(3,:,curr_animal) = arrayfun(@(depth) ...
        nanmean(AP_itril(corr(permute(curr_k_str(:,depth,:),[1,3,2])),-1)),1:n_depths);
    ctx_str_taskk_corr(4,:,curr_animal) = arrayfun(@(depth) ...
        nanmean(AP_itril(corr(permute(curr_k_ctx(:,depth,:),[1,3,2])),-1)),1:n_depths);  

end

% Get mean across domains
ctx_str_taskk_corr_strmean = squeeze(nanmean(ctx_str_taskk_corr,2));

% Plot mean and split by domains
figure; 

subplot(2,1,1);hold on; set(gca,'ColorOrder',copper(n_depths));
plot(ctx_str_taskk_corr_strmean,'color',[0.5,0.5,0.5]);
errorbar(nanmean(ctx_str_taskk_corr_strmean,2), ...
    AP_sem(ctx_str_taskk_corr_strmean,2),'k','linewidth',2);
set(gca,'XTick',1:4,'XTickLabelRotation',20,'XTickLabel', ...
    {'Str-ctx within day','Str within day across domains','Str across days','Ctx across days'})
ylabel('Task kernel correlation');
xlim([0.5,4.5]);

subplot(2,1,2);hold on; set(gca,'ColorOrder',copper(n_depths));
errorbar(nanmean(ctx_str_taskk_corr,3), ...
    AP_sem(ctx_str_taskk_corr,3),'linewidth',2)
set(gca,'XTick',1:4,'XTickLabelRotation',20,'XTickLabel', ...
    {'Str-ctx within day','Str within day across domains','Str across days','Ctx across days'})
ylabel('Task kernel correlation');
xlim([0.5,4.5]);
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false))

% (within task-passive v task-task domains statistics)
disp('Str/ctx vs str/str cross-domain:')
curr_p = signrank(squeeze(ctx_str_taskk_corr_strmean(1,:)), ...
    squeeze(ctx_str_taskk_corr_strmean(2,:)));
disp(['All str p = ' num2str(curr_p)]);
for curr_depth = 1:n_depths  
    curr_p = signrank(squeeze(ctx_str_taskk_corr(1,curr_depth,:)), ...
        squeeze(ctx_str_taskk_corr(2,curr_depth,:)));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end

% (within vs across statistics)
disp('Str/ctx-within vs str-across:')
curr_p = signrank(squeeze(ctx_str_taskk_corr_strmean(1,:)), ...
    squeeze(ctx_str_taskk_corr_strmean(3,:)));
disp(['All str ' num2str(curr_depth) ' p = ' num2str(curr_p)]);
for curr_depth = 1:n_depths
    curr_p = signrank(squeeze(ctx_str_taskk_corr(1,curr_depth,:)), ...
        squeeze(ctx_str_taskk_corr(3,curr_depth,:)));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end

% (cross task vs no task statistics)
disp('Str-across vs ctx-across');
curr_p = signrank(squeeze(ctx_str_taskk_corr_strmean(3,:)), ...
    squeeze(ctx_str_taskk_corr_strmean(4,:)));
disp(['All str ' num2str(curr_depth) ' p = ' num2str(curr_p)]);
for curr_depth = 1:n_depths
    curr_p = signrank(squeeze(ctx_str_taskk_corr(3,curr_depth,:)), ...
        squeeze(ctx_str_taskk_corr(4,curr_depth,:)));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end


%% ** Fig 3f: cortex vs task explained variance

% Use raw data (not normalized or baseline-subtracted) for expl var
mua_exp = vertcat(mua_all{:});
mua_taskpred_exp = vertcat(mua_taskpred_all{:});
mua_ctxpred_exp = vertcat(mua_ctxpred_all{:});

% Get R^2 for task and cortex
taskpred_r2 = nan(max(split_idx),n_depths);
ctxpred_r2 = nan(max(split_idx),n_depths);
for curr_exp = 1:max(split_idx)
       
    curr_data = reshape(permute(mua_exp{curr_exp},[2,1,3]),[],n_depths);
    curr_data_baselinesub = reshape(permute(mua_exp{curr_exp},[2,1,3]),[],n_depths) - ...
        (nanmean(reshape(mua_exp{curr_exp}(:,t < 0,:),[],size(mua_exp{curr_exp},3)),1));
    curr_taskpred_data = reshape(permute(mua_taskpred_exp{curr_exp},[2,1,3]),[],n_depths);
    curr_ctxpred_data = reshape(permute(mua_ctxpred_exp{curr_exp},[2,1,3]),[],n_depths);
       
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_data_baselinesub) | ...
        isnan(curr_taskpred_data) | isnan(curr_ctxpred_data);
    curr_data(nan_samples) = NaN;
    curr_data_baselinesub(nan_samples) = NaN;
    curr_taskpred_data(nan_samples) = NaN;
    curr_ctxpred_data(nan_samples) = NaN;

    % (task regressed from average baseline-subtracted data)
    taskpred_r2(curr_exp,:) = 1 - (nansum((curr_data_baselinesub-curr_taskpred_data).^2,1)./ ...
        nansum((curr_data_baselinesub-nanmean(curr_data_baselinesub,1)).^2,1));
    % (cortex regressed from raw data)
    ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    
end

figure; hold on;
errorbar(nanmean(taskpred_r2,1),AP_sem(taskpred_r2,1),'b','linewidth',2,'CapSize',0);
errorbar(nanmean(ctxpred_r2,1),AP_sem(ctxpred_r2,1),'color',[0,0.6,0],'linewidth',2,'CapSize',0);
xlabel('Striatum depth');
ylabel('Task explained variance');
legend({'Task','Cortex'});


% Plot explained variance task vs cortex by experiment
figure; hold on;
str_col = max(hsv(n_depths)-0.2,0);
for curr_str = 1:n_depths
    errorbar(squeeze(nanmean(taskpred_r2(:,curr_str),1)), ...
        squeeze(nanmean(ctxpred_r2(:,curr_str),1)), ...
        squeeze(AP_sem(ctxpred_r2(:,curr_str),1)),squeeze(AP_sem(ctxpred_r2(:,curr_str),1)), ...
        squeeze(AP_sem(taskpred_r2(:,curr_str),1)),squeeze(AP_sem(taskpred_r2(:,curr_str),1)), ...
        'color','k','linewidth',2);
    
    scatter(taskpred_r2(:,curr_str),ctxpred_r2(:,curr_str),10, ...
        str_col(curr_str,:),'filled');
    scatter(nanmean(taskpred_r2(:,curr_str),1), ...
        nanmean(ctxpred_r2(:,curr_str),1),80, ...
        str_col(curr_str,:),'filled','MarkerEdgeColor','k','linewidth',2);
end
axis tight equal;
line(xlim,xlim,'color','k','linestyle','--');
xlabel('Task R^2');
ylabel('Cortex R^2');
legend({'DMS','DCS','DLS'})

% (Task R2 statistics)
disp('Task R^2 values:');
for curr_depth = 1:n_depths
    disp(['Str ' num2str(curr_depth) ...
        ' R2 = ' num2str(nanmean(taskpred_r2(:,curr_depth))) '+/- ' ...
        num2str(AP_sem(taskpred_r2(:,curr_depth),1))]); 
end

% (Cortex R2 statistics)
disp('Cortex R^2 values:');
for curr_depth = 1:n_depths
    disp(['Str ' num2str(curr_depth) ...
        ' R2 = ' num2str(nanmean(ctxpred_r2(:,curr_depth))) '+/- ' ...
        num2str(AP_sem(ctxpred_r2(:,curr_depth),1))]); 
end

% (Cortex vs task R2 statistics)
disp('Cortex vs Task R^2 signrank:');
for curr_depth = 1:n_depths
    curr_p = signrank(ctxpred_r2(:,curr_depth), ...
        taskpred_r2(:,curr_depth));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end


%% @@ Fig 4b,c: Average activity within and across cells by domain/celltype

% Get trial params and data by experiment
trial_stim_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
move_t_exp = mat2cell(move_t,use_split,1);
trial_outcome_allcat_exp = mat2cell(trial_outcome_allcat,use_split,1);

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

% (striatum)
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

% (cortex)
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


% Set align breakpoints
align_col = [1,0,0;0.8,0,0.8;0,0,0.8];
% (split the alignment halfway between median alignment points)
align_median = cellfun(@(x) -nanmedian(x)/sample_rate,use_align);
align_break = align_median(1:end-1) + diff(align_median*0.8);
align_t = {[-0.05,align_break(1)],[align_break(1:2)],[align_break(2),1]};

% Plot aligned_activity (all cells)
plot_celltypes = 1:4;

figure;
for curr_depth = 1:n_aligned_depths
    for curr_celltype = plot_celltypes
        
        curr_cells = domain_aligned_allcat == curr_depth & ...
            celltype_allcat == curr_celltype & good_units_allcat & ...
            any(any(act_mean_plot,3),2);
        curr_cells_idx = find(curr_cells);
        [~,max_idx] = max(act_mean_sort(curr_cells,:,1),[],2);
        [~,sort_idx] = sort(max_idx,'descend');
   
        for curr_align = 1:length(use_align)
            
            subplot(n_aligned_depths,length(plot_celltypes), ...
                (curr_depth-1)*length(plot_celltypes)+curr_celltype); hold on;
            
            curr_t_offset = -nanmedian(use_align{curr_align})/sample_rate;
            curr_t = t + curr_t_offset;
            curr_t_plot = curr_t >= align_t{curr_align}(1) & ...
                curr_t <= align_t{curr_align}(2);
            
            curr_act_sorted = act_mean_plot(curr_cells_idx(sort_idx),:,curr_align);
            curr_act_sorted_norm = curr_act_sorted./max(curr_act_sorted,[],2);
        
            % smooth if too many cells to plot accurately
            % (taken out: safer bet to rasterize in illustrator)
%             if sum(curr_cells) > 500
%                 n_smooth = 5;
%                 curr_act_sorted_norm = convn(curr_act_sorted_norm, ...
%                     ones(n_smooth,1)./n_smooth,'same');
%             end

            plot_t = curr_t > align_t{curr_align}(1) & curr_t <= align_t{curr_align}(2);
            
            % (scalebar)
            if sum(curr_cells) > 500
                line([0,0],[0,200],'color','g','linewidth',2);
            else
                line([0,0],[0,20],'color','g','linewidth',2);
            end
            
            imagesc(curr_t(plot_t),[],curr_act_sorted_norm(:,plot_t));
            caxis([0,1]);
            axis tight off;
            line(repmat(curr_t_offset,2,1),ylim,'color',align_col(curr_align,:));
            colormap(gca,brewermap([],'Greys'));
            ylabel('Neuron (sorted)');
            xlabel('~Time from stim');
%             title(['Str ' num2str(curr_depth) ' ' celltype_labels{curr_celltype}]);
            
        end       
    end
end

% Plot aligned_activity (average across cells)
figure;
for curr_depth = 1:n_aligned_depths
    for curr_celltype = plot_celltypes
        
        curr_cells = domain_aligned_allcat == curr_depth & ...
            celltype_allcat == curr_celltype & good_units_allcat & ...
            any(any(act_mean_plot,3),2);   
   
        for curr_align = 1:length(use_align)
            
            subplot(n_aligned_depths,length(plot_celltypes), ...
                (curr_depth-1)*length(plot_celltypes)+curr_celltype); hold on;
            
            curr_t_offset = -nanmedian(use_align{curr_align})/sample_rate;
            curr_t = t + curr_t_offset;
            curr_t_plot = curr_t >= align_t{curr_align}(1) & ...
                curr_t <= align_t{curr_align}(2);
            
            curr_act = act_mean_plot(curr_cells,:,curr_align);
            
            plot_t = curr_t > align_t{curr_align}(1) & curr_t <= align_t{curr_align}(2);
            
            AP_errorfill(curr_t(plot_t)', ...
                nanmean(curr_act(:,plot_t),1)', ...
                AP_sem(curr_act(:,plot_t),1)','k',[],[],1);
            % (Matteo's suggestion plot decile: looks crazy)
%             plot(curr_t(plot_t),prctile(curr_act(:,plot_t),linspace(10,100,90),1)','color',[0.5,0.5,0.5]);
%             plot(curr_t(plot_t),nanmean(curr_act(:,plot_t),1),'color','k');
            
            % (scalebar on first plot)
            if curr_depth == 1 && curr_celltype == 1 && curr_align == 1
               line([0,0.5],repmat(min(ylim),2,1),'color','g','linewidth',2); 
            end
            
            axis tight;
            line(repmat(curr_t_offset,2,1),ylim,'color',align_col(curr_align,:));
            colormap(gca,brewermap([],'Greys'));
            ylabel('Spikes/s');
%             xlabel('~Time from stim');
%             title(['Str ' num2str(curr_depth) ' ' celltype_labels{curr_celltype}]);
            set(gca,'FontSize',9,'FontName','Myriad Pro');
            set(gca,'XTick',[]);
            
        end       
    end
end


% Plot average cortical kernel ROI activity
figure;
for curr_depth = 1:n_aligned_depths
    for curr_align = 1:length(use_align)
        
        subplot(n_aligned_depths,1,curr_depth); hold on;
        
        curr_t_offset = -nanmedian(use_align{curr_align})/sample_rate;
        curr_t = t + curr_t_offset;
        curr_t_plot = curr_t >= align_t{curr_align}(1) & ...
            curr_t <= align_t{curr_align}(2);
        
        curr_act = permute(ctx_kernelroi_act_mean(curr_depth,:,:,curr_align),[3,2,1,4]);
        
        plot_t = curr_t > align_t{curr_align}(1) & curr_t <= align_t{curr_align}(2);
        
        AP_errorfill(curr_t(plot_t)', ...
            nanmean(curr_act(:,plot_t),1)', ...
            AP_sem(curr_act(:,plot_t),1)',[0,0.7,0],[],[],1);
                
    end
end
linkaxes(get(gcf,'Children'),'xy');




%% @@ Fig 4d: SUA-celltype/cortex activity correlation

include_celltypes = 1:4;

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
    for curr_celltype = include_celltypes
        for curr_celltype_compare = include_celltypes
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
    for curr_celltype = include_celltypes
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
for curr_celltype = include_celltypes
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
    [length(include_celltypes),length(fr_prctilebin_centers),max(recordings_allcat),size(celltype_allcorr,2)],@nanmean,NaN);
celltype_fr_frbins = accumarray(grp_matrix,fr(use_units), ...
    [length(include_celltypes),length(fr_prctilebin_centers),max(recordings_allcat)],@nanmean,NaN);


figure;
celltype_col = ...
    [0.9,0.4,0.6; ...
    0.4,0.5,0.6; ...
    0.5,0.3,0.1; ...
    1,0.5,0];
ctx_col = [0,0.7,0];
p = gobjects(2,3);
for curr_celltype = include_celltypes
    
    curr_fr = permute(squeeze(celltype_fr_frbins(curr_celltype,:,:,:)),[1,3,2]);
    curr_actcorr = permute(squeeze(celltype_actcorr_frbins(curr_celltype,:,:,:)),[1,3,2]);
    
    p(1,curr_celltype) = subplot(2,length(include_celltypes),curr_celltype); hold on;
    set(gca,'ColorOrder',[celltype_col(include_celltypes,:);ctx_col],'YScale','log');
    errorbar(nanmean(curr_actcorr,3), ...
        repmat(nanmean(curr_fr,3),1,size(celltype_allcorr,2)), ...
        AP_sem(curr_actcorr,3),'horizontal','linewidth',2);
    title(celltype_labels{curr_celltype})
    legend([celltype_labels(include_celltypes),'Cortex ROI']);
    xlabel('Avg act corr');
    ylabel('Firing rate');

    p(2,curr_celltype) = subplot(2,length(include_celltypes),length(include_celltypes)+curr_celltype); hold on;
    set(gca,'ColorOrder',[celltype_col(include_celltypes,:);ctx_col],'XScale','log');
    errorbar(repmat(nanmean(curr_fr,3),1,size(celltype_allcorr,2)), ...
        nanmean(curr_actcorr,3),AP_sem(curr_actcorr,3),'linewidth',2);
    title(celltype_labels{curr_celltype})
    legend([celltype_labels(include_celltypes),'Cortex ROI']);
    xlabel('Firing rate');
    ylabel('Avg act corr');
    axis tight
    curr_xlim = xlim;
    xlim([curr_xlim(1)*0.7,curr_xlim(2)*1.3]);
    
end
linkaxes(p(1,:),'xy');
linkaxes(p(2,:),'y');

% (statistics comparing celltypes/cortex)
disp('MSN/FSI - FSI/MSN v CTX ANOVA:')
for curr_celltype = 1:3   
    curr_actcorr = permute(celltype_actcorr_frbins(curr_celltype,:,:,:),[2,4,3,1]);
    [fr_grp,compare_grp,recording_grp] = ndgrid( ...
        1:size(curr_actcorr,1), ...
        1:size(curr_actcorr,2), ...
        1:size(curr_actcorr,3));
    use_comparison = setdiff([1,2,4],curr_celltype); % (msn v fsi v cortex)
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


% (if UINs included:)
disp('MSN/FSI - FSI/MSN v CTX ANOVA:')
for curr_celltype = 4   
    curr_actcorr = permute(celltype_actcorr_frbins(curr_celltype,:,:,:),[2,4,3,1]);
    [fr_grp,compare_grp,recording_grp] = ndgrid( ...
        1:size(curr_actcorr,1), ...
        1:size(curr_actcorr,2), ...
        1:size(curr_actcorr,3));
    use_comparison = setdiff([1,2,4],curr_celltype); % (msn v fsi v cortex)
    p = anovan(reshape(curr_actcorr(:,use_comparison,:),[],1), ...
        [reshape(fr_grp(:,use_comparison,:),[],1), ...
        reshape(compare_grp(:,use_comparison,:),[],1)], ...
        'model','interaction','display','off');
    disp([celltype_labels{curr_celltype} ' v (msn/fsi/ctx): p(fr,type,interaction) = ' num2str(p')])    
end

disp('TAN ANOVA:')
for curr_celltype = 4
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


%% @@ ^^ Fig 4d shuffle statistics

% Stats for difference in activity correlation between cell types
% (compared to shuffling 2 cell type labels within day)

shuffle_celltypes = [2,3];

n_shuff = 1000;
shuffle_diff = nan(3,4,n_shuff);
for curr_shuff = 1:n_shuff
    
    if curr_shuff == 1
        % Get the real value on the first shuffle
        celltype_allcat_shuff = celltype_allcat;
    else
        % Shuffle labels across chosen pair within day
        celltype_exp = mat2cell(celltype_allcat,cellfun(@(x) size(x,3),mua_exp),1);
        celltype_exp_shuff = celltype_exp;
        for curr_exp = 1:length(celltype_exp_shuff)
            curr_shuff_idx = ismember(celltype_exp_shuff{curr_exp},shuffle_celltypes);
            if sum(curr_shuff_idx) > 1
                celltype_exp_shuff{curr_exp}(curr_shuff_idx) = ...
                    AP_shake(celltype_exp_shuff{curr_exp}(curr_shuff_idx));
            end
        end
        celltype_allcat_shuff = cell2mat(celltype_exp_shuff);
    end
    
    % Get cell-average correlations (across recording)
    n_recordings = max(recordings_allcat);
    celltype_act_corr = nan(size(str_unit_actmean_multialign,1),3);
    for curr_depth = 1:n_aligned_depths
        for curr_celltype = include_celltypes
            for curr_celltype_compare = include_celltypes
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
                        celltype_allcat_shuff == curr_celltype_compare & ...
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
    
    % Concatenate celltype and fluorescence correlation
    celltype_allcorr = [celltype_act_corr];
    
    % Get correlation mean and in FR percentile bins
    fr = nanmean(str_unit_actmean_multialign,2);
    n_prctiles = 4;
    fr_prctilebin_edges = linspace(0,100,n_prctiles+1);
    fr_prctilebin_centers = fr_prctilebin_edges(1:end-1)+diff(fr_prctilebin_edges)./2;
    
    fr_prctilebins = nan(size(fr));
    for curr_celltype = include_celltypes
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
        [length(include_celltypes),length(fr_prctilebin_centers),max(recordings_allcat),size(celltype_allcorr,2)],@nanmean,NaN);
    celltype_fr_frbins = accumarray(grp_matrix,fr(use_units), ...
        [length(include_celltypes),length(fr_prctilebin_centers),max(recordings_allcat)],@nanmean,NaN);

    curr_shuffle_diff = ...
        nanmean(celltype_actcorr_frbins(:,:,:,shuffle_celltypes(1)) - ...
        celltype_actcorr_frbins(:,:,:,shuffle_celltypes(2)),3);    
    shuffle_diff(:,:,curr_shuff) = curr_shuffle_diff;
    
    AP_print_progress_fraction(curr_shuff,n_shuff);
    
end

% Get pvalue (first value is the real difference, others are shuffled)
corr_rank = tiedrank(permute(nanmean(shuffle_diff,2),[3,1,2]));
corr_p = 1 - (corr_rank(1,:)/(n_shuff+1));

disp(['Shuffle: ' (celltype_labels{shuffle_celltypes})]);
for curr_celltype = 1:3
    disp([celltype_labels{curr_celltype} ' , p = ' num2str(corr_p(curr_celltype))]);
end



%% Fig 5a,b: Untrained/trained passive cortex/striatum mua

data_fns = { ...
    'trial_activity_AP_choiceWorldStimPassive_naive', ...
    {'trial_activity_AP_choiceWorldStimPassive_trained', ...
    'trial_activity_AP_lcrGratingPassive_ctxstrephys_str', ...
    'trial_activity_AP_lcrGratingPassive_pre_muscimol'}};

% % (leave out original dataset: 1s instead of 0.5s stim)
% data_fns = { ...
%     'trial_activity_AP_choiceWorldStimPassive_naive', ...
%     {'trial_activity_AP_lcrGratingPassive_ctxstrephys_str', ...
%     'trial_activity_AP_lcrGratingPassive_pre_muscimol'}};

stimIDs = cell(2,1);
mua_training = cell(2,1);
mua_ctxpred_training = cell(2,1);
fluor_training = cell(2,1);
fluor_roi_training = cell(2,1);
fluor_kernelroi_training = cell(2,1);

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;

    % Split by experiment
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    
    % Get trials with movement during stim to exclude
    wheel_thresh = 0.025;
    quiescent_trials = ~any(abs(wheel_allcat(:,t >= -0.5 & t <= 0.5)) > wheel_thresh,2);
    quiescent_trials_exp = mat2cell(quiescent_trials,trials_recording,1);
    
    % Get stim and activity by experiment
    trial_stim_allcat_exp = mat2cell(trial_stim_allcat,trials_recording,1);  
    fluor_deconv_exp = mat2cell(fluor_allcat_deconv,trials_recording,length(t),n_vs);
    fluor_roi_deconv_exp = mat2cell(fluor_roi_deconv,trials_recording,length(t),numel(wf_roi));
    fluor_kernelroi_deconv_exp = mat2cell(fluor_kernelroi_deconv,trials_recording,length(t),n_depths);
    mua_allcat_exp = mat2cell(mua_allcat,trials_recording,length(t),n_depths);
    mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,trials_recording,length(t),n_depths);
    
    % Exclude trials with fluorescence spikes
    % (this is a dirty way to do this but don't have a better alt)
    fluor_spike_thresh = 100;
    fluor_spike_trial = cellfun(@(x) any(any(x > fluor_spike_thresh,2),3), ...
        fluor_kernelroi_deconv_exp,'uni',false);
    
    % Grab data
    stimIDs{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        trial_stim_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    mua_training{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        mua_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    mua_ctxpred_training{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        mua_ctxpred_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_training{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_roi_training{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_roi_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_kernelroi_training{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_kernelroi_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
end

% Get average stim time course (100%R stim)
mua_mean = cellfun(@(act,stim) cell2mat(cellfun(@(act,stim) ...
    nanmean(act(stim == 1,:,:),1),act,stim,'uni',false)), ...
    mua_training,stimIDs,'uni',false);

fluor_kernelroi_mean = cellfun(@(act,stim) cell2mat(cellfun(@(act,stim) ...
    nanmean(act(stim == 1,:,:),1),act,stim,'uni',false)), ...
    fluor_kernelroi_training,stimIDs,'uni',false);

% Get average activity in relevant stim period
stim_avg_t = [0,0.2];
stim_avg_t_idx = t >= stim_avg_t(1) & t <= stim_avg_t(2);

% (only for 100%R stim)
stim_act = cellfun(@(str,ctx,stim) cellfun(@(str,ctx,stim) ...
    [nanmean(str(stim == 1,stim_avg_t_idx,:),2), ...
    nanmean(ctx(stim == 1,stim_avg_t_idx,:),2)], ...
    str,ctx,stim,'uni',false), ...
    mua_training, fluor_kernelroi_training, stimIDs,'uni',false);


% Plot average cortex and striatum stim activity
figure;
p = gobjects(n_depths,3);
for curr_str = 1:n_depths
    
    p(curr_str,1) = subplot(n_depths,3,(curr_str-1)*3+1);
    AP_errorfill(t,nanmean(fluor_kernelroi_mean{1}(:,:,curr_str),1)', ...
        AP_sem(fluor_kernelroi_mean{1}(:,:,curr_str),1)','k');
    AP_errorfill(t,nanmean(fluor_kernelroi_mean{2}(:,:,curr_str),1)', ...
        AP_sem(fluor_kernelroi_mean{2}(:,:,curr_str),1),'r');
    ylabel('Cortex (\DeltaF/F)');
    line(repmat(stim_avg_t(1),2,1),ylim);
    line(repmat(stim_avg_t(2),2,1),ylim);
    
    p(curr_str,2) = subplot(n_depths,3,(curr_str-1)*3+2);
    AP_errorfill(t,nanmean(mua_mean{1}(:,:,curr_str),1)', ...
        AP_sem(mua_mean{1}(:,:,curr_str),1)','k');
    AP_errorfill(t,nanmean(mua_mean{2}(:,:,curr_str),1)', ...
        AP_sem(mua_mean{2}(:,:,curr_str),1)','r');
    xlabel('Time from stim (s)');
    ylabel('Striatum (baseline)');
    title(['Str ' num2str(curr_str)]);
    line(repmat(stim_avg_t(1),2,1),ylim);
    line(repmat(stim_avg_t(2),2,1),ylim);
    
    
    p(curr_str,3) = subplot(n_depths,3,(curr_str-1)*3+3);  hold on;
    
    curr_stim_act_untrained = cell2mat(cellfun(@(x) nanmean(x(:,:,curr_str,1),1),stim_act{1},'uni',false));
    curr_stim_act_trained = cell2mat(cellfun(@(x) nanmean(x(:,:,curr_str,1),1),stim_act{2},'uni',false));
    
    str_col = copper(n_depths);
    
    plot([nanmean(curr_stim_act_untrained(:,2),1),nanmean(curr_stim_act_trained(:,2),1)], ...
       [nanmean(curr_stim_act_untrained(:,1),1),nanmean(curr_stim_act_trained(:,1),1)], ...
       'color',str_col(curr_str,:),'linewidth',2);
      
    scatter(curr_stim_act_untrained(:,2), ...
        curr_stim_act_untrained(:,1),10, ...
        [0.2,0.2,0.2],'filled');
    scatter(curr_stim_act_trained(:,2), ...
        curr_stim_act_trained(:,1),10, ...
        [1,0.2,0.2],'filled');

   errorbar(nanmean(curr_stim_act_untrained(:,2),1), ...
        nanmean(curr_stim_act_untrained(:,1),1), ...
        AP_sem(curr_stim_act_untrained(:,1),1),AP_sem(curr_stim_act_untrained(:,1),1), ...
        AP_sem(curr_stim_act_untrained(:,2),1),AP_sem(curr_stim_act_untrained(:,2),1), ...
        '.','color','k','linewidth',2,'MarkerSize',20);
    errorbar(nanmean(curr_stim_act_trained(:,2),1), ...
        nanmean(curr_stim_act_trained(:,1),1), ...
        AP_sem(curr_stim_act_trained(:,1),1),AP_sem(curr_stim_act_trained(:,1),1), ...
        AP_sem(curr_stim_act_trained(:,2),1),AP_sem(curr_stim_act_trained(:,2),1), ...
        '.','color','r','linewidth',2,'MarkerSize',20);
   
   xlabel('Cortical ROI');
   ylabel('Striatum');
    
end
linkaxes(p(:,1),'xy');
linkaxes(p(:,2),'xy');
linkaxes(p(:,3),'xy');


% Plot kernel ROIs
kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
load(kernel_roi_fn);

figure; hold on
n_depths = size(kernel_roi.bw,3);
str_col = max(hsv(n_depths)-0.2,0);
set(gca,'YDir','reverse');
AP_reference_outline('ccf_aligned','k');
for curr_depths = 1:n_depths
    curr_roi_boundary = cell2mat(bwboundaries(kernel_roi.bw(:,:,curr_depths)));
    patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),str_col(curr_depths,:));
end
axis image off;


% (Untrained/trained statistics)
stim_act_mean = cellfun(@(x) cell2mat(cellfun(@(x) nanmean(x,1), ...
    x,'uni',false)),stim_act,'uni',false);

disp('Untrained/trained ranksum:');
for curr_depth = 1:n_depths
    curr_p = ranksum(stim_act_mean{1}(:,1,curr_depth),stim_act_mean{2}(:,1,curr_depth));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
    
    curr_p = ranksum(stim_act_mean{1}(:,2,curr_depth),stim_act_mean{2}(:,2,curr_depth));
    disp(['Ctx ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end

%% Fig 5c: Untrained/trained passive striatum cell types

data_fns = { ...
    'trial_activity_naive_sua', ...
    {'trial_activity_trainedPassive_sua.mat', ... % original
    'trial_activity_ctx_passive_sua.mat', ...     % + cortex ephys
    'trial_activity_muscimol_passive_sua.mat'}};  % muscimol group

stim_act_celltype_training = cell(size(data_fns));
stim_act_celltype_group = cell(size(data_fns));

for curr_data = 1:length(data_fns)
    
    preload_vars = who;
    
    data_fn = data_fns{curr_data};
    
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
%     %%%%%%%%%%%% TO NOT BASELINE-SUBTRACT
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

%% Fig 5d: Naive cortical kernels: compare to trained

% Load Master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Load trained data
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\str_ctxpred.mat')
% Concatenate kernels and convert to pixels
% (and flip in time so it's fluorescence lead:lag spikes)
ctx_str_k_trained_animal = [str_ctxpred.ctx_str_k]';
ctx_str_k_trained_animal = cellfun(@(x) cellfun(@(x) ...
    flip(AP_svdFrameReconstruct(U_master(:,:,1:100),x),3),x,'uni',false), ...
    ctx_str_k_trained_animal,'uni',false);

ctx_str_k_px_trained_cat = vertcat(ctx_str_k_trained_animal{:}); 
ctx_str_k_px_notask_trained_mean = nanmean(cat(5,ctx_str_k_px_trained_cat{:,2}),5);

% Load naive data
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\str_ctxpred_naive.mat')
% Concatenate kernels and convert to pixels
% (and flip in time so it's fluorescence lead:lag spikes)
ctx_str_k_naive_animal = [str_ctxpred.ctx_str_k]';
ctx_str_k_naive_animal = cellfun(@(x) cellfun(@(x) ...
    flip(AP_svdFrameReconstruct(U_master(:,:,1:100),x),3),x,'uni',false), ...
    ctx_str_k_naive_animal,'uni',false);

ctx_str_k_px_naive_cat = vertcat(ctx_str_k_naive_animal{:}); 
ctx_str_k_px_notask_naive_mean = nanmean(cat(5,ctx_str_k_px_naive_cat{:,1}),5);


% Get time
framerate = 35;
upsample_factor = 1;
sample_rate = framerate*upsample_factor;
kernel_t = [-0.1,0.1];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
t = kernel_frames/sample_rate;

% Get mean kernels and plot
n_depths = size(ctx_str_k_px_notask_trained_mean,4);

AP_image_scroll([ctx_str_k_px_notask_trained_mean,ctx_str_k_px_notask_naive_mean]);
axis image;
colormap(brewermap([],'PRGn'));
caxis([-max(abs(caxis)),max(abs(caxis))]);
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5],[],[size(U_master,1),size(U_master,2),1,2]);

figure;
for curr_depth = 1:n_depths
    subplot(n_depths,2,(curr_depth-1)*2+1);
    imagesc(ctx_str_k_px_notask_trained_mean(:,:,t == 0,curr_depth));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'PRGn'));
    axis image off;
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title('Task');
    
    subplot(n_depths,2,(curr_depth-1)*2+2);
    imagesc(ctx_str_k_px_notask_naive_mean(:,:,t == 0,curr_depth));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'PRGn'));
    axis image off;
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title('Passive');
end

figure;
for curr_depth = 1:n_depths
    subplot(n_depths,1,curr_depth);
    imagesc(reshape(ctx_str_k_px_notask_trained_mean(:,:,:,curr_depth), ...
        size(ctx_str_k_px_notask_trained_mean,1),[]));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'PRGn'));
    axis image off;
    title('Task');
end

figure;
for curr_depth = 1:n_depths
    subplot(n_depths,1,curr_depth);
    imagesc(reshape(ctx_str_k_px_notask_naive_mean(:,:,:,curr_depth), ...
        size(ctx_str_k_px_notask_naive_mean,1),[]));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'PRGn'));
    axis image off;
    title('Passive');
end


% Correlate trained/passive kernels within/across domain
ctx_str_k_trained_animalmean = cellfun(@(x) ...
    reshape(nanmean(cat(5,x{:,2}),5),[],n_depths),ctx_str_k_trained_animal,'uni',false);

ctx_str_k_naive_animalmean = cellfun(@(x) ...
    reshape(nanmean(cat(5,x{:,1}),5),[],n_depths),ctx_str_k_naive_animal,'uni',false);

ctx_str_k_corr = corr([ ...
    horzcat(ctx_str_k_trained_animalmean{:}), ...
    horzcat(ctx_str_k_naive_animalmean{:})],'type','Pearson');

ctx_str_k_corr_tril = tril(ctx_str_k_corr,-1);
ctx_str_k_corr_tril(triu(true(size(ctx_str_k_corr)))) = NaN;

corr_grp = ...
    [repmat(1:n_depths,1, ...
    length(ctx_str_k_trained_animalmean)+ ...
    length(ctx_str_k_naive_animalmean)); ...
    1*ones(1,length(ctx_str_k_trained_animalmean)*n_depths), ...
    2*ones(1,length(ctx_str_k_naive_animalmean)*n_depths)];

k_corr_trained_withindomain = ctx_str_k_corr_tril(...
    corr_grp(2,:) == 1 & corr_grp(2,:)' == 1 & ...
    corr_grp(1,:) == corr_grp(1,:)');

k_corr_trained_acrossdomain = ctx_str_k_corr_tril(...
    corr_grp(2,:) == 1 & corr_grp(2,:)' == 1 & ...
    corr_grp(1,:) ~= corr_grp(1,:)');

k_corr_naive_withindomain = ctx_str_k_corr_tril(...
    corr_grp(2,:) == 2 & corr_grp(2,:)' == 2 & ...
    corr_grp(1,:) == corr_grp(1,:)');

k_corr_naive_acrossdomain = ctx_str_k_corr_tril(...
    corr_grp(2,:) == 2 & corr_grp(2,:)' == 2 & ...
    corr_grp(1,:) ~= corr_grp(1,:)');

k_corr_trainednaive_withindomain = ctx_str_k_corr_tril(...
    corr_grp(2,:) == 1 & corr_grp(2,:)' == 2 & ...
    corr_grp(1,:) == corr_grp(1,:)');

k_corr_trainednaive_acrossdomain = ctx_str_k_corr_tril(...
    corr_grp(2,:) == 1 & corr_grp(2,:)' == 2 & ...
    corr_grp(1,:) ~= corr_grp(1,:)');

k_corr_grp = {k_corr_trained_withindomain, ...
    k_corr_trainednaive_withindomain, ...
    k_corr_trainednaive_acrossdomain};

figure; hold on;
distributionPlot(k_corr_grp,'showMM',0,'color',[0.5,0.5,0.5])
plot(cellfun(@nanmean,k_corr_grp),'k','linewidth',2);
ylabel('Kernel spatiotemporal correlation');
set(gca,'XTickLabel',{'Trained within domain', ...
    'Trained-naive within domain', ...
    'Trained-naive across domain'}, ...
    'XTickLabelRotation',45);

% (cross task vs no task statistics)
curr_p = ranksum(k_corr_grp{1},k_corr_grp{2});
disp(['Trained within domain vs. trained-naive within domain p = ' ...
    num2str(curr_p)])
curr_p = ranksum(k_corr_grp{1},k_corr_grp{3});
disp(['Trained within domain vs. trained-naive across domain p = ' ...
    num2str(curr_p)])


%% ~~~~~~~~~~~~~ Extended data figures


%% ** EDFig1b,c,d: Psychometric and reaction time

figure;

% Plot psychometric
stim_conditions = unique(trial_stim_allcat);
[~,stim_idx] = ismember(trial_stim_allcat,stim_conditions,'rows');

trial_stim_idx_allcat_exp = mat2cell(stim_idx,trials_animal,1);
trial_choice_allcat_exp = mat2cell(trial_choice_allcat,trials_animal,1);

frac_orient_right = cell2mat(cellfun(@(stim,choice) ...
    accumarray(stim,choice == -1,[length(stim_conditions),1],@nanmean,NaN), ...
    trial_stim_idx_allcat_exp,trial_choice_allcat_exp,'uni',false)');

subplot(1,3,1); hold on;
axis square;
AP_errorfill(stim_conditions,nanmean(frac_orient_right,2), ...
    AP_sem(frac_orient_right,2),'k');
line([0,0],ylim,'color','k','linestyle','--');
line(xlim,[0.5,0.5],'color','k','linestyle','--');
xlabel('Stimulus side*contrast');
ylabel('Fraction orient right');

% Plot reaction times by stim
subplot(1,3,2); hold on;
axis square;
rxn_time_stim = cell2mat(cellfun(@(stim,move_t) ...
    accumarray(stim,move_t,[length(stim_conditions),1],@nanmedian,NaN), ...
    trial_stim_idx_allcat_exp,mat2cell(move_t,trials_animal,1),'uni',false)');
AP_errorfill(stim_conditions,nanmean(rxn_time_stim,2), ...
    AP_sem(rxn_time_stim,2),'k');
ylim([0,0.6])
line([0,0],ylim,'color','k','linestyle','--');
line(xlim,[0.5,0.5],'color','k','linestyle','--');
xlabel('Stimulus side*contrast')
ylabel('Median reaction time');

% Plot reaction time by trial percentile within session
n_trial_prctile = 4;
trial_bin = arrayfun(@(x) min(floor(linspace(1,n_trial_prctile+1,x)), ...
    n_trial_prctile)',trials_recording,'uni',false);

move_t_bins = -0.2:1/sample_rate:1;
move_t_bin_centers = move_t_bins(1:end-1) + diff(move_t_bins)./2;
move_t_bin = mat2cell(discretize(move_t,move_t_bins),trials_recording,1);

move_t_hist = cell2mat(permute(cellfun(@(trial_bin,move_t_bin) ...
    accumarray([trial_bin(~isnan(move_t_bin)),move_t_bin(~isnan(move_t_bin))], ...
    1,[n_trial_prctile,length(move_t_bins)-1],@nansum,0), ...
    trial_bin,move_t_bin,'uni',false),[2,3,1]));

move_t_frac = move_t_hist./sum(move_t_hist,2);

n_recording_animal = cellfun(@length,wheel_all);
move_t_frac_animals = mat2cell(move_t_frac,n_trial_prctile,length(move_t_bin_centers),n_recording_animal);
move_t_frac_animalmean = cell2mat(cellfun(@(x) nanmean(x,3), ...
    move_t_frac_animals,'uni',false));

subplot(1,3,3);
AP_errorfill(move_t_bin_centers,nanmean(move_t_frac_animalmean,3)', ...
    AP_sem(move_t_frac_animalmean,3)',copper(n_trial_prctile));
xlabel('Time from stim onset');
ylabel('Fraction of reaction times');
line([0,0],ylim,'color','k');
title('Trial quartiles');


%% EDFig2a,b,c,d: Widefield alignment

% Show average images across days for one animal
animal = 'AP025';

protocol = 'vanillaChoiceworld';
experiments = AP_find_experiments(animal,protocol);
experiments = experiments([experiments.imaging] & [experiments.ephys]);

% (load retinotopy for that animal)
retinotopy_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment\retinotopy';
retinotopy_fn = [retinotopy_path filesep animal '_retinotopy'];
load(retinotopy_fn);

figure;
for curr_day = 1:length(experiments)
    
    day = experiments(curr_day).day;
    
    [img_path,img_exists] = AP_cortexlab_filename(animal,day,[],'imaging');
    avg_im = readNPY([img_path filesep 'meanImage_purple.npy']);
        
    subplot(2,length(experiments),curr_day);
    imagesc(avg_im);
    axis image off;
    colormap(gca,gray);
    caxis([0,40000]);
    title([animal ' Day ' num2str(curr_day)]);
    
    curr_day_idx = strcmp(day,{retinotopy.day});
    subplot(2,length(experiments),length(experiments)+curr_day);
    imagesc(retinotopy(curr_day_idx).vfs);
    axis image off;
    colormap(gca,brewermap([],'*RdBu'));
    caxis([-1,1]);
    
end

% Plot retinotopy for all animals
retinotopy_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\retinotopy\retinotopy';
load(retinotopy_fn);

figure;
sq_subplots = ceil(sqrt(length(retinotopy)));
for curr_animal = 1:length(retinotopy)
   subplot(sq_subplots,sq_subplots,curr_animal);
   imagesc(retinotopy(curr_animal).vfs);
   axis image off
   colormap(brewermap([],'*RdBu'));
   caxis([-1,1]);
   title(retinotopy(curr_animal).animal);
end


% Show colored CCF and master retinotopy with aligned CCF
% (these steps are normally done in AP_vfs_ccf_align)

% Load CCF av and st
ccf_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
av = readNPY([ccf_path filesep 'annotation_volume_10um_by_index.npy']); 
st = loadStructureTree([ccf_path filesep 'structure_tree_safe_2017.csv']);

% Get first brain pixel from top-down, get annotation at that point
[~,top_down_depth] = max(av>1, [], 2);
top_down_depth = squeeze(top_down_depth);

[xx,yy] = meshgrid(1:size(top_down_depth,2), 1:size(top_down_depth,1));
top_down_annotation = reshape(av(sub2ind(size(av),yy(:),top_down_depth(:),xx(:))), size(av,1), size(av,3));

% Get all labelled areas
used_areas = unique(top_down_annotation(:));

% Restrict to only cortical areas
structure_id_path = cellfun(@(x) textscan(x(2:end),'%d', 'delimiter',{'/'}),st.structure_id_path);

ctx_path = [997,8,567,688,695,315];
ctx_idx = find(cellfun(@(id) length(id) > length(ctx_path) & ...
    all(id(min(length(id),length(ctx_path))) == ctx_path(min(length(id),length(ctx_path)))),structure_id_path));

plot_areas = intersect(used_areas,ctx_idx);

bregma = allenCCFbregma;

% Get outlines of all areas
top_down_cortical_area_boundaries = cell(size(plot_areas));
for curr_area_idx = 1:length(plot_areas)
    top_down_cortical_area_boundaries{curr_area_idx} = bwboundaries(top_down_annotation == plot_areas(curr_area_idx));
end

% Color CCF by VFS
a_idx = find(cellfun(@(name) strcmp(name,'Anterior area layer 1'),st.safe_name(used_areas)));
al_idx = find(cellfun(@(name) strcmp(name,'Anterolateral visual area layer 1'),st.safe_name(used_areas)));
am_idx = find(cellfun(@(name) strcmp(name,'Anteromedial visual area layer 1'),st.safe_name(used_areas)));
lm_idx = find(cellfun(@(name) strcmp(name,'Lateral visual area layer 1'),st.safe_name(used_areas)));
v1_idx = find(cellfun(@(name) strcmp(name,'Primary visual area layer 1'),st.safe_name(used_areas)));
p_idx = find(cellfun(@(name) strcmp(name,'Posterolateral visual area layer 1'),st.safe_name(used_areas)));
pm_idx = find(cellfun(@(name) strcmp(name,'posteromedial visual area layer 1'),st.safe_name(used_areas)));
li_idx = find(cellfun(@(name) strcmp(name,'Laterointermediate area layer 1'),st.safe_name(used_areas)));
rl_idx = find(cellfun(@(name) strcmp(name,'Rostrolateral area layer 1'),st.safe_name(used_areas)));

ccf_vfs = zeros(size(top_down_annotation));
ccf_vfs(ismember(top_down_annotation,used_areas([v1_idx,am_idx,al_idx,li_idx]))) = 1;
ccf_vfs(ismember(top_down_annotation,used_areas([a_idx,p_idx,pm_idx,rl_idx,lm_idx]))) = -1;

% Load master retinotopy
alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment';
master_vfs_fn = [alignment_path filesep 'master_vfs.mat'];
load(master_vfs_fn);

% Plot CCF VFS, and master VFS
figure;
subplot(1,2,1,'YDir','reverse'); hold on;
imagesc(ccf_vfs);
cellfun(@(area) plot(area(:,2),area(:,1),'color',[0.5,0.5,0.5]), ...
    vertcat(top_down_cortical_area_boundaries{:}),'uni',false);
axis image off;
colormap(brewermap([],'*RdBu'));
caxis([-1,1]);    
title('CCF');

subplot(1,2,2);
imagesc(master_vfs);
axis image off;
colormap(brewermap([],'*RdBu'));
caxis([-1,1]);
title('Combined');
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);




%% EDFig2e,f: Widefield correlation borders

wf_corr_borders_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_borders\wf_corr_borders.mat';
load(wf_corr_borders_fn);

% Spacing/downsampling (hardcoded - in preprocessing)
px_spacing = 20;
downsample_factor = 10;

% Get average correlation maps
wf_corr_map_recording = [wf_corr_borders(:).corr_map_downsamp];
wf_corr_map_cat = cat(3,wf_corr_map_recording{:});
wf_corr_map_mean = cell(size(wf_corr_map_cat,1),size(wf_corr_map_cat,2));
for curr_y = 1:size(wf_corr_map_cat,1)
    for curr_x = 1:size(wf_corr_map_cat,2)
        wf_corr_map_mean{curr_y,curr_x} = ...
            nanmean(cell2mat(wf_corr_map_cat(curr_y,curr_x,:)),3);
    end
end

% Get average borders
wf_corr_borders_cat = cell2mat(reshape([wf_corr_borders(:).corr_edges],1,1,[]));
wf_corr_borders_mean = nanmean(wf_corr_borders_cat,3);

figure;

[plot_maps_y,plot_maps_x] = ndgrid(4:5:size(wf_corr_map_mean,1)-4,4:5:size(wf_corr_map_mean,2));

%%%%%%%%%%%%%%%%%%%% This is unused: k-means of correlated pixels
% (CCF lines lie at reflection points and correlated points are tangent to
% that, that's probably too confusing for anyone not used to this and I
% don't know how to pull out reflection points, so leaving out)

% (get widefield inside brain)
h = figure('visible','off');
brain_outlines = AP_reference_outline('ccf_aligned');
brain_outlines = vertcat(brain_outlines{:});

% (get brain pixels in widefield using aligned CCF borders)
[m,n] = size(wf_corr_borders_mean);
brain_pixels = false(m,n);
for curr_region = 1:length(brain_outlines)
    curr_bw = poly2mask(brain_outlines{curr_region}.XData, ...
        brain_outlines{curr_region}.YData,m,n);
    brain_pixels(curr_bw) = true;
end
brain_pixels = imclose(brain_pixels,strel('disk',5'));
close(h);

corr_map_cat = cat(3,wf_corr_map_mean{:});
corr_map_cat_flat = reshape(corr_map_cat,[],size(corr_map_cat,3));
corr_map_cat_flat(isnan(corr_map_cat_flat)) = 0;

n_grps = 6;
brain_pixels_downsamp = imresize(brain_pixels,1/downsample_factor,'nearest');
kidx = kmeans(corr_map_cat_flat(brain_pixels_downsamp(:),:),n_grps);

kidx_brain = nan(size(brain_pixels_downsamp));
kidx_brain(brain_pixels_downsamp(:)) = kidx;

kidx_upsample = imresize(kidx_brain,downsample_factor,'nearest');

figure;imagesc(kidx_upsample)
AP_reference_outline('ccf_aligned','k');
axis image off;
colormap([1,1,1;brewermap([],'Set3')]);
caxis([0,n_grps]);
title('Pixel correlation k-means');
%%%%%%%%%%%%%%%%%%%% 

% Plot sample correlation map locations
subplot(1,3,1,'YDir','reverse');
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
plot(plot_maps_y(:)*px_spacing,plot_maps_x(:)*px_spacing,'.r','MarkerSize',30);
axis image off;

% Plot sample correlation maps concat
subplot(1,3,2);
imagesc(cell2mat(wf_corr_map_mean(5:5:end,5:5:end)));
axis image off;
caxis([0,1]);
colormap(brewermap([],'Greys'));
c = colorbar;
ylabel(c,'Correlation');

% Plot borders
subplot(1,3,3);
imagesc(wf_corr_borders_mean);
axis image off;
caxis([0,prctile(wf_corr_borders_mean(:),70)])
colormap(brewermap([],'Greys'))
ccf_outline = AP_reference_outline('ccf_aligned',[1,0,0]);
cellfun(@(x) set(x,'linewidth',1),vertcat(ccf_outline{:}));


%% EDFig3a: Probe trajectories histology vs widefield-estimated

% Load probe trajectories
histology_probe_ccf_all_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\histology_probe_ccf_all'];
load(histology_probe_ccf_all_fn);
n_animals = length(histology_probe_ccf_all);

% Load in the annotated Allen volume and names
allen_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
av = readNPY([allen_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean

figure;

for curr_animal = 1:n_animals
    
    % Plot average image with drawn probe overlay;
    subplot(n_animals,3,(curr_animal-1)*3+1);
    imagesc(histology_probe_ccf_all(curr_animal).im);
    line(histology_probe_ccf_all(curr_animal).probe_wf(:,1), ...
        histology_probe_ccf_all(curr_animal).probe_wf(:,2), ...
        'color','r','linewidth',1);  
    axis image off;
    caxis([0,prctile(histology_probe_ccf_all(curr_animal).im(:),95)]);
    colormap(gca,'gray');
    
    % Plot brain to overlay probes
    horizontal_axes = subplot(n_animals,3,(curr_animal-1)*3+2,'YDir','reverse'); hold on;
    axis image off;
    coronal_axes = subplot(n_animals,3,(curr_animal-1)*3+3,'YDir','reverse'); hold on;
    axis image off;
    
    % Plot projections
    coronal_outline = bwboundaries(permute((max(av,[],1)) > 1,[2,3,1]));
    horizontal_outline = bwboundaries(permute((max(av,[],2)) > 1,[3,1,2]));
    
    str_id = find(strcmp(st.safe_name,'Caudoputamen'));
    str_coronal_outline = bwboundaries(permute((max(av == str_id,[],1)) > 0,[2,3,1]));
    str_horizontal_outline = bwboundaries(permute((max(av == str_id,[],2)) > 0,[3,1,2]));
    
    vstr_id = find(strcmp(st.safe_name,'Nucleus accumbens'));
    vstr_coronal_outline = bwboundaries(permute((max(av == vstr_id,[],1)) > 0,[2,3,1]));
 
    cellfun(@(x) plot(horizontal_axes,x(:,2),x(:,1),'k','linewidth',2),horizontal_outline)
    cellfun(@(x) plot(horizontal_axes,x(:,2),x(:,1),'b','linewidth',2),str_horizontal_outline)
    axis image off;
    
    cellfun(@(x) plot(coronal_axes,x(:,2),x(:,1),'k','linewidth',2),coronal_outline)
    cellfun(@(x) plot(coronal_axes,x(:,2),x(:,1),'b','linewidth',2),str_coronal_outline)
    cellfun(@(x) plot(coronal_axes,x(:,2),x(:,1),'color',[0.8,0,0.8],'linewidth',2),vstr_coronal_outline)
    axis image off;
    
    % Get estimated probe vector (widefield)
    probe_ccf_wf = histology_probe_ccf_all(curr_animal).probe_ccf_wf;
    r0 = mean(probe_ccf_wf,1);
    xyz = bsxfun(@minus,probe_ccf_wf,r0);
    [~,~,V] = svd(xyz,0);
    probe_direction = V(:,1);
    
    probe_vector_evaluate = [0,sign(probe_direction(2))*1000];
    probe_vector_draw_wf = round(bsxfun(@plus,bsxfun(@times,probe_vector_evaluate',probe_direction'),r0));
    
    line(horizontal_axes,probe_vector_draw_wf(:,1),probe_vector_draw_wf(:,3), ...
        'color','r','linewidth',2);
    line(coronal_axes,probe_vector_draw_wf(:,3),probe_vector_draw_wf(:,2), ...
        'color','r','linewidth',2);
    
    % Get estimated probe vector (histology)
    probe_ccf_histology = histology_probe_ccf_all(curr_animal).probe_ccf_histology;
    % (still working on orienting histology: force left hemisphere)
    bregma = allenCCFbregma;
    if sign(probe_ccf_histology(end,3) - probe_ccf_histology(1,3)) == 1
        probe_ccf_histology(:,3) = ...
            probe_ccf_histology(:,3) + ...
            2*(bregma(3) - probe_ccf_histology(:,3));
    end
    
    r0 = mean(probe_ccf_histology,1);
    xyz = bsxfun(@minus,probe_ccf_histology,r0);
    [~,~,V] = svd(xyz,0);
    probe_direction = V(:,1);
    
    probe_vector_evaluate = [-sign(probe_direction(2))*500,sign(probe_direction(2))*500];
    probe_vector_draw_histology = round(bsxfun(@plus,bsxfun(@times,probe_vector_evaluate',probe_direction'),r0));

    line(horizontal_axes,probe_vector_draw_histology(:,1),probe_vector_draw_histology(:,3), ...
        'color',[0,0.7,0],'linewidth',2);
    line(coronal_axes,probe_vector_draw_histology(:,3),probe_vector_draw_histology(:,2), ...
        'color',[0,0.7,0],'linewidth',2);
    
    drawnow;
    
end


%% EDFig3b: Probe trajectories estimated from widefield image

% Load estimated probe trajectories
probe_ccf_all_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\probe_ccf_all'];
load(probe_ccf_all_fn);

% Load in the annotated Allen volume and names
allen_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
av = readNPY([allen_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean

% Plot brain to overlay probes
% (note the CCF is rotated to allow for dim 1 = x)
h = figure; 
ccf_axes = subplot(2,2,1); hold on
sagittal_axes = subplot(2,2,2,'YDir','reverse'); hold on;
axis image; grid on;
horizontal_axes = subplot(2,2,3,'YDir','reverse'); hold on;
axis image; grid on;
coronal_axes = subplot(2,2,4,'YDir','reverse'); hold on;
axis image; grid on;

% Plot 1 = 3D
% (Use wire mesh - can add other structures)
slice_spacing = 10;

str_id = find(strcmp(st.safe_name,'Caudoputamen'));
target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) == str_id,[2,1,3]);
structure_patch = isosurface(target_volume,0);
structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
target_structure_color = [0,0,1];
striatum_outline = patch(ccf_axes, ...
    'Vertices',structure_wire.vertices*slice_spacing, ...
    'Faces',structure_wire.faces, ...
    'FaceAlpha',0.5,'FaceColor',target_structure_color,'EdgeColor','none');

% (Use Nick's outline, but rotate to make dim 1 = x)
brainGridData = readNPY([fileparts(which('plotBrainGrid')) filesep 'brainGridData.npy']);
plotBrainGrid(brainGridData(:,[1,3,2]),ccf_axes);

axes(ccf_axes);
axis image vis3d off;
cameratoolbar(h,'SetCoordSys','y');
cameratoolbar(h,'SetMode','orbit');

% Plots 2-4 = projections
coronal_outline = bwboundaries(permute((max(av,[],1)) > 1,[2,3,1]));
horizontal_outline = bwboundaries(permute((max(av,[],2)) > 1,[3,1,2]));
sagittal_outline = bwboundaries(permute((max(av,[],3)) > 1,[2,1,3]));

str_id = find(strcmp(st.safe_name,'Caudoputamen'));
str_coronal_outline = bwboundaries(permute((max(av == str_id,[],1)) > 0,[2,3,1]));
str_horizontal_outline = bwboundaries(permute((max(av == str_id,[],2)) > 0,[3,1,2]));
str_sagittal_outline = bwboundaries(permute((max(av == str_id,[],3)) > 0,[2,1,3]));

vstr_id = find(strcmp(st.safe_name,'Nucleus accumbens'));
vstr_coronal_outline = bwboundaries(permute((max(av == vstr_id,[],1)) > 0,[2,3,1]));

cellfun(@(x) plot(sagittal_axes,x(:,2),x(:,1),'k','linewidth',2),sagittal_outline)
cellfun(@(x) plot(sagittal_axes,x(:,2),x(:,1),'b','linewidth',2),str_sagittal_outline)

cellfun(@(x) plot(horizontal_axes,x(:,2),x(:,1),'k','linewidth',2),horizontal_outline)
cellfun(@(x) plot(horizontal_axes,x(:,2),x(:,1),'b','linewidth',2),str_horizontal_outline)

cellfun(@(x) plot(coronal_axes,x(:,2),x(:,1),'k','linewidth',2),coronal_outline)
cellfun(@(x) plot(coronal_axes,x(:,2),x(:,1),'b','linewidth',2),str_coronal_outline)
cellfun(@(x) plot(coronal_axes,x(:,2),x(:,1),'color',[0.8,0,0.8],'linewidth',2),vstr_coronal_outline)

if length(probe_ccf_all) <= 17
    color_set = [brewermap(8,'Dark2');brewermap(8,'Set2')];
    probe_color = color_set(1:length(probe_ccf_all),:);
else
    error('More animals than colors: add another set');
end

for curr_animal = 1:length(probe_ccf_all)
    for curr_day = 1:size(probe_ccf_all{curr_animal},3)

        % Get estimated probe vector
        probe_ccf = probe_ccf_all{curr_animal}(:,:,curr_day);
        r0 = mean(probe_ccf,1);
        xyz = bsxfun(@minus,probe_ccf,r0);
        [~,~,V] = svd(xyz,0);
        probe_direction = V(:,1);
                
        probe_vector_evaluate = [0,sign(probe_direction(2))*1000];
        probe_vector_draw = round(bsxfun(@plus,bsxfun(@times,probe_vector_evaluate',probe_direction'),r0));
        
        % Plot current probe vector
        line(ccf_axes,probe_vector_draw(:,1),probe_vector_draw(:,2),probe_vector_draw(:,3), ...
            'color',probe_color(curr_animal,:),'linewidth',1);
        line(sagittal_axes,probe_vector_draw(:,1),probe_vector_draw(:,2), ...
            'color',probe_color(curr_animal,:),'linewidth',1);
        line(horizontal_axes,probe_vector_draw(:,1),probe_vector_draw(:,3), ...
            'color',probe_color(curr_animal,:),'linewidth',1);
        line(coronal_axes,probe_vector_draw(:,3),probe_vector_draw(:,2), ...
            'color',probe_color(curr_animal,:),'linewidth',1);
        
        drawnow;
        
    end
end

% Put a colormap on the side
cmap_ax = axes('Position',[0.95,0.2,0.2,0.6]);
image(permute(probe_color,[1,3,2]));




%% EDFig3d: Striatum border estimation (histology and electrophysiology)

animal = 'AP032';

% Get probe areas estimated from histology
load(['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\histology\' animal '\processed\probe_areas']);

% Get and plot striatum boundaries
ephys_align_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
ephys_align_fn = ['ephys_depth_align'];
load([ephys_align_path filesep ephys_align_fn])

% Parameters from batch
n_corr_groups = 40;
depth_group_edges = linspace(0,3820,n_corr_groups+1);
depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);

curr_animal = find(strcmp({ephys_depth_align.animal},animal));

% Plot first day and histology-estimated boundaries
str_id = 574;
str_start_histology = probe_areas.y(find(probe_areas.av == str_id,1,'first'));
str_end_histology = probe_areas.y(find(probe_areas.av == str_id,1,'last'));

curr_day = 1;
mua_corr = ephys_depth_align(curr_animal).mua_corr{curr_day};
template_depths = ephys_depth_align(curr_animal).template_depths{curr_day};
str_depth = ephys_depth_align(curr_animal).str_depth(curr_day,:);

figure;
imagesc(depth_group_centers,depth_group_centers,mua_corr);
caxis([0,1]); colormap(hot); axis square;
line(xlim,[str_start_histology,str_start_histology],'color','w')
line([str_start_histology,str_start_histology],ylim,'color','w')
line(xlim,[str_end_histology,str_end_histology],'color','w')
line([str_end_histology,str_end_histology],ylim,'color','w')
xlabel('Depth (\mum)');
ylabel('Depth (\mum)');
title(['Day ' num2str(curr_day) ' (histology-approximated)']);
    
% Plot all days and ephys-estimated boundaries
n_days = length(ephys_depth_align(curr_animal).mua_corr);
figure;
for curr_day = 1:n_days
    mua_corr = ephys_depth_align(curr_animal).mua_corr{curr_day};
    template_depths = ephys_depth_align(curr_animal).template_depths{curr_day};
    str_depth = ephys_depth_align(curr_animal).str_depth(curr_day,:);
    
    subplot(1,n_days,curr_day);
    imagesc(depth_group_centers,depth_group_centers,mua_corr);
    caxis([0,1]); colormap(hot); axis square;
    line(xlim,[str_depth(1),str_depth(1)],'color','w');
    line([str_depth(1),str_depth(1)],ylim,'color','w');
    line(xlim,[str_depth(2),str_depth(2)],'color','w');
    line([str_depth(2),str_depth(2)],ylim,'color','w');
    xlabel('Depth (\mum)');
    ylabel('Depth (\mum)');
    title(['Day ' num2str(curr_day)]);
end
set(gcf,'Name',ephys_depth_align(curr_animal).animal);


%% EDFig4a,b: Fluorescence deconvolution and fluor+ctx+str example

% Load and plot kernel
% (flip time to be fluor lag:lead spikes);
kernel_path = fileparts(which('AP_deconv_wf'));
kernel_fn = [kernel_path filesep 'gcamp6s_kernel.mat'];
load(kernel_fn);

gcamp6s_kernel_cat = fliplr(vertcat(gcamp6s_kernel.regression{:}));
gcamp6s_kernel_norm = gcamp6s_kernel_cat./max(abs(gcamp6s_kernel_cat),[],2);
gcamp6s_kernel_mean = nanmean(gcamp6s_kernel_norm);

figure; hold on;
plot(gcamp6s_kernel.regression_t,gcamp6s_kernel_norm','color',[0.5,0.5,0.5]);
plot(gcamp6s_kernel.regression_t,gcamp6s_kernel_mean,'k','linewidth',2);
line(xlim,[0,0],'color','k','linestyle','--');
line([0,0],ylim,'color','k','linestyle','--');


% Load Fluor+ctx-ephys data
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\ctx_deconv_traces');

% Concatenate spikes/deconvolved fluorescence 
% (full recording)
ctx_fluor_full_cat = cellfun(@cell2mat,{ctx_deconv_traces.fluor},'uni',false);
ctx_spikes_full_cat = cellfun(@cell2mat,{ctx_deconv_traces.spikes},'uni',false);
ctx_fluor_deconv_full_cat = cellfun(@cell2mat,{ctx_deconv_traces.fluor_deconv},'uni',false);

% (task)
ctx_spikes_task_cat = cellfun(@(x,protocol) ...
    cell2mat(x(strcmp(protocol,'vanillaChoiceworld'))), ...
    {ctx_deconv_traces.spikes},{ctx_deconv_traces.protocol},'uni',false);
ctx_fluor_deconv_task_cat = cellfun(@(x,protocol) ...
    cell2mat(x(strcmp(protocol,'vanillaChoiceworld'))), ...
    {ctx_deconv_traces.fluor_deconv},{ctx_deconv_traces.protocol},'uni',false);

% (passive)
ctx_spikes_notask_cat = cellfun(@(x,protocol) ...
    cell2mat(x(strcmp(protocol,'AP_sparseNoise'))), ...
    {ctx_deconv_traces.spikes},{ctx_deconv_traces.protocol},'uni',false);
ctx_fluor_deconv_notask_cat = cellfun(@(x,protocol) ...
    cell2mat(x(strcmp(protocol,'AP_sparseNoise'))), ...
    {ctx_deconv_traces.fluor_deconv},{ctx_deconv_traces.protocol},'uni',false);

% Get explained variance
ctx_deconv_full_r2 = cellfun(@(spikes,fluor_deconv) ...
    1-(nansum((spikes-fluor_deconv).^2)./nansum((spikes-nanmean(spikes)).^2)), ...
    ctx_spikes_full_cat,ctx_fluor_deconv_full_cat);

ctx_deconv_task_r2 = cellfun(@(spikes,fluor_deconv) ...
    1-(nansum((spikes-fluor_deconv).^2)./nansum((spikes-nanmean(spikes)).^2)), ...
    ctx_spikes_task_cat,ctx_fluor_deconv_task_cat);

ctx_deconv_notask_r2 = cellfun(@(spikes,fluor_deconv) ...
    1-(nansum((spikes-fluor_deconv).^2)./nansum((spikes-nanmean(spikes)).^2)), ...
    ctx_spikes_notask_cat,ctx_fluor_deconv_notask_cat);


%%% Plot example day

animal = 'AP060';
day = '2019-12-06';
experiment = 1;
plot_t = [100,300];

figure;
disp('Loading example data...');

% Load cortex ephys + imaging
load_parts.ephys = true;
load_parts.imaging = true;
site = 2; % (cortex always probe 2)
str_align = 'none'; % (cortex)
AP_load_experiment;

% Load cortex recording alignment
vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
load(vis_ctx_ephys_fn);

%%% DEPTH-ALIGN TEMPLATES, FIND CORTEX BOUNDARY
curr_animal_idx = strcmp(animal,{vis_ctx_ephys.animal});
curr_day_idx = strcmp(day,vis_ctx_ephys(curr_animal_idx).day);
curr_csd_depth = vis_ctx_ephys(curr_animal_idx).stim_csd_depth{curr_day_idx};
curr_csd_depth_aligned = vis_ctx_ephys(curr_animal_idx).stim_csd_depth_aligned{curr_day_idx};
template_depths_aligned = interp1(curr_csd_depth,curr_csd_depth_aligned,template_depths);
spike_depths_aligned = interp1(curr_csd_depth,curr_csd_depth_aligned,spike_depths);

% Find cortex end by largest gap between templates
sorted_template_depths = sort([template_depths_aligned]);
[max_gap,max_gap_idx] = max(diff(sorted_template_depths));
ctx_end = sorted_template_depths(max_gap_idx)+1;

ctx_depth = [sorted_template_depths(1),ctx_end];
ctx_units = template_depths_aligned <= ctx_depth(2);

%%% GET FLUORESCENCE AND SPIKES BY DEPTH

% Set binning time
skip_seconds = 60;
spike_binning_t = 1/framerate; % seconds
spike_binning_t_edges = frame_t(1)+skip_seconds:spike_binning_t:frame_t(end)-skip_seconds;
spike_binning_t_centers = spike_binning_t_edges(1:end-1) + diff(spike_binning_t_edges)/2;

% Get fluorescence in pre-drawn ROI
curr_ctx_roi = vis_ctx_ephys(curr_animal_idx).ctx_roi{curr_day_idx};

fVdf_deconv = AP_deconv_wf(fVdf);
fluor_roi = AP_svd_roi(Udf,fVdf_deconv,avg_im,[],curr_ctx_roi);
fluor_roi_interp = interp1(frame_t,fluor_roi,spike_binning_t_centers);

% Plot cortex raster
subplot(6,1,1); hold on;
plot_t_idx = spike_binning_t_centers >= plot_t(1) & ...
    spike_binning_t_centers <= plot_t(2);
plot(spike_binning_t_centers(plot_t_idx), ...
    fluor_roi_interp(plot_t_idx),'linewidth',2,'color',[0,0.7,0]);

subplot(6,1,2,'YDir','reverse'); hold on;
plot_spikes = spike_times_timeline >= plot_t(1) & ...
    spike_times_timeline <= plot_t(2) & ...
    spike_depths_aligned <= ctx_depth(2);
plot(spike_times_timeline(plot_spikes),spike_depths(plot_spikes),'.k');
ylabel('Cortex depth (\mum)');
xlabel('Time (s)');

%%% LOAD STRIATUM EPHYS AND GET MUA BY DEPTH

% Load striatum ephys
load_parts.ephys = true;
load_parts.imaging = false;
site = 1; % (striatum is always on probe 1)
str_align = 'kernel';
AP_load_experiment;

% Plot striatum raster
subplot(6,1,3,'YDir','reverse'); hold on;
plot_spikes = spike_times_timeline >= plot_t(1) & ...
    spike_times_timeline <= plot_t(2) & ...
    spike_depths >= str_depth(1) & spike_depths <= str_depth(2);
plot(spike_times_timeline(plot_spikes),spike_depths(plot_spikes),'.k');
ylabel('Striatum depth (\mum)');
xlabel('Time (s)');

%%%%%% PLOT WHEEL/STIM

% (wheel velocity)
wheel_axes = subplot(6,1,4); hold on;
plot_wheel_idx = Timeline.rawDAQTimestamps >= plot_t(1) & ...
    Timeline.rawDAQTimestamps <= plot_t(2);
plot(wheel_axes,Timeline.rawDAQTimestamps(plot_wheel_idx), ...
    wheel_velocity(plot_wheel_idx),'k','linewidth',2);
ylabel('Wheel velocity');
axis off

% %     (stimuli)
%     % (task)
%     if contains(expDef,'vanilla')
%         stim_col = colormap_BlueWhiteRed(5);
%         [~,trial_contrast_idx] = ...
%             ismember(trial_conditions(:,1).*trial_conditions(:,2),unique(contrasts'.*sides),'rows');
%     elseif strcmp(expDef,'AP_lcrGratingPassive')
%         % (passive)
%         stim_col = [0,0,1;0.5,0.5,0.5;1,0,0];
%         [~,trial_contrast_idx] = ...
%             ismember(trial_conditions(:,1).*trial_conditions(:,2),[-90;0;90],'rows');
%     end
%     stim_lines = arrayfun(@(x) line(wheel_axes,repmat(stimOn_times(x),1,2),ylim(wheel_axes),'color', ...
%             stim_col(trial_contrast_idx(x),:),'linewidth',2), ...
%             find(stimOn_times >= plot_t(1) & stimOn_times <= plot_t(2)));
%
%     % (movement starts)
%     move_col = [0.6,0,0.6;0,0.6,0];
%     [~,trial_choice_idx] = ismember(trial_conditions(:,3),[-1;1],'rows');
%     move_lines = arrayfun(@(x) line(wheel_axes,repmat(wheel_move_time(x),1,2),ylim(wheel_axes),'color', ...
%         move_col(trial_choice_idx(x),:),'linewidth',2), ...
%         find(wheel_move_time >= plot_t(1) & wheel_move_time <= plot_t(2)));
%
%     % (go cues)
%     go_col = [0.8,0.8,0.2];
%     go_cue_times = signals_events.interactiveOnTimes(1:n_trials);
%     go_cue_lines = arrayfun(@(x) line(wheel_axes,repmat(go_cue_times(x),1,2),ylim(wheel_axes),'color', ...
%         go_col,'linewidth',2), ...
%         find(go_cue_times >= plot_t(1) & go_cue_times <= plot_t(2)));
%
%     % (outcomes)
%     outcome_col = [0,0,0.8;0.5,0.5,0.5];
%     reward_lines = arrayfun(@(x) line(wheel_axes,repmat(reward_t_timeline(x),1,2),ylim(wheel_axes),'color', ...
%         outcome_col(1,:),'linewidth',2), ...
%         find(reward_t_timeline >= plot_t(1) & reward_t_timeline <= plot_t(2)));
%     punish_times = signals_events.responseTimes(trial_outcome == -1);
%     punish_lines = arrayfun(@(x) line(wheel_axes,repmat(punish_times(x),1,2),ylim(wheel_axes),'color', ...
%         outcome_col(2,:),'linewidth',2), ...
%         find(punish_times >= plot_t(1) & punish_times <= plot_t(2)));


% Plot the cortex MUA and deconvolved fluorescence match
% Plot example
example_recording = 3;
example_protocol = 1;

% (raw fluorescence)
subplot(6,1,5); hold on;
plot(ctx_deconv_traces(example_recording).t{example_protocol}, ...
    mat2gray(ctx_deconv_traces(example_recording).fluor{example_protocol}), ...
    'color',[0,0.8,0],'linewidth',2);
ylabel('ROI Fluor(\DeltaF/F)');
xlabel('Time(s)')

% (cortex MUA and deconvolved fluorescence)
subplot(6,1,6); hold on;
plot(ctx_deconv_traces(example_recording).t{example_protocol}, ...
    ctx_deconv_traces(example_recording).spikes{example_protocol}, ...
    'color','k','linewidth',2);
plot(ctx_deconv_traces(example_recording).t{example_protocol}, ...
    ctx_deconv_traces(example_recording).fluor_deconv{example_protocol}, ...
    'color',[0,0.5,0],'linewidth',2);
linkaxes(get(gcf,'Children'),'x');
ylabel('Cortex spikes (std)');
xlabel('Time (s)');


curr_axes = flipud(get(gcf,'Children'));
% Link all time axes
linkaxes(curr_axes,'x');
% Link depth axes of raster plots (arbitrary depth, but want same scale)
linkaxes(curr_axes(2:3),'xy');

% (Display average and statistics)
disp(['Deconv explained var: ' num2str(nanmean(ctx_deconv_full_r2)) ...
    ' +/- ' num2str(AP_sem(ctx_deconv_full_r2,2))]);

p = signrank(ctx_deconv_task_r2,ctx_deconv_notask_r2);
disp(['Task vs no-task explained var: ' num2str(p)]);


%% EDFig4c,d: Fluorescence/cortex ephys/striatum ephys correlation

%%% Load correlation data
use_protocol = 'vanillaChoiceworld';
% use_protocol = 'AP_sparseNoise';

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
data_fn = [data_path filesep 'ctx_fluor_mua_corr_' use_protocol];
load(data_fn);

%%% Get average aligned CSD

% Load CSD 
vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
load(vis_ctx_ephys_fn);

% Concatenate for plotting
animal_cat = cellfun(@(x,y) repmat({x},1,length(y)),{vis_ctx_ephys.animal},{vis_ctx_ephys.day},'uni',false);
animal_cat = horzcat(animal_cat{:});
day_cat = [vis_ctx_ephys(:).day];
n_recordings = length(day_cat);

stim_lfp_t = vis_ctx_ephys(1).stim_lfp_t{1};
stim_lfp_depth = vis_ctx_ephys(1).stim_lfp_depth{1};
stim_csd_depth = vis_ctx_ephys(1).stim_csd_depth{1};

stim_lfp_cat = cell2mat(permute([vis_ctx_ephys(:).stim_lfp],[1,3,2]));
stim_csd_cat = cell2mat(permute([vis_ctx_ephys(:).stim_csd],[1,3,2]));
csd_slice_cat = cell2mat([vis_ctx_ephys(:).csd_slice]);

% Get aligned depth for each recording, plot, save
figure;
curr_recording = 1;
stim_csd_aligned_scaled_cat = nan(0);
for curr_animal = 1:length(vis_ctx_ephys)
    for curr_day = 1:length(vis_ctx_ephys(curr_animal).day)
        
        stim_csd = vis_ctx_ephys(curr_animal).stim_csd{curr_day};
        csd_slice = vis_ctx_ephys(curr_animal).csd_slice{curr_day};
        stim_csd_depth = vis_ctx_ephys(curr_animal).stim_csd_depth{curr_day};
        stim_csd_depth_aligned = vis_ctx_ephys(curr_animal).stim_csd_depth_aligned{curr_day};
                
        % Plot aligned CSD
        depth_align_interp = [-500:20:2000];
        stim_csd_aligned = interp1(...
            stim_csd_depth_aligned,stim_csd_cat(:,:,curr_recording),depth_align_interp, ...
            'linear','extrap');
        
        animal = vis_ctx_ephys(curr_animal).animal;
        day = vis_ctx_ephys(curr_animal).day{curr_day};
        stim_lfp_t = vis_ctx_ephys(curr_animal).stim_lfp_t{curr_day};
        stim_lfp_depth = vis_ctx_ephys(curr_animal).stim_lfp_depth{curr_day};
        stim_lfp = vis_ctx_ephys(curr_animal).stim_lfp{curr_day};
        
        subplot(3,n_recordings,curr_recording);
        imagesc(stim_lfp_t,stim_lfp_depth,stim_lfp);
        caxis([-max(abs(caxis))/2,max(abs(caxis))/2]);
        colormap(brewermap([],'*RdBu'));
        title({animal,day,'LFP'});
        
        subplot(3,n_recordings,size(stim_csd_cat,3)+curr_recording);
        imagesc(stim_lfp_t,stim_csd_depth,stim_csd);
        caxis([-max(abs(caxis))/2,max(abs(caxis))/2]);
        colormap(brewermap([],'*RdBu'));
        title('CSD');
        
        subplot(3,n_recordings,size(stim_csd_cat,3)*2+curr_recording);
        imagesc(stim_lfp_t,stim_csd_depth_aligned,stim_csd_aligned);
        caxis([-max(abs(caxis))/2,max(abs(caxis))/2]);
        colormap(brewermap([],'*RdBu'));
        title('Aligned CSD');
                
        % Keep aligned CSD for averaging
        stim_csd_aligned_scaled_cat(:,:,curr_recording) = stim_csd_aligned./min(csd_slice);       
        
        % Iterate recording index (only used for plotting)
        curr_recording = curr_recording + 1;

    end
end

stim_csd_aligned_mean = nanmean(stim_csd_aligned_scaled_cat,3);

%%% Plot correlations

mua_depth = data(1).cortex_mua_depth{1}; % they're all the same, use 1st
cortex_fluor_corr_cat = cell2mat(horzcat(data.cortex_fluor_corr));
cortex_striatum_corr_cat = cell2mat(permute(horzcat(data.cortex_striatum_corr),[1,3,2]));

figure;

subplot(1,2,1,'YDir','reverse');
imagesc(stim_lfp_t,depth_align_interp,stim_csd_aligned_mean);
caxis([-1,1]);
line([0,0],ylim,'color','k');
colormap(brewermap([],'*RdBu'));
xlabel('Time from stim');
ylabel('Aligned depth');
title('Average aligned CSD')
ylim([0,1400]);

subplot(1,2,2); hold on;
p1 = AP_errorfill(mua_depth, ...
    nanmean(cortex_fluor_corr_cat,2), ...
    AP_sem(cortex_fluor_corr_cat,2),[0,0.6,0]);

plot_str = 1;
p2 = AP_errorfill(mua_depth, ...
    nanmean(cortex_striatum_corr_cat(:,plot_str,:),3), ...
    AP_sem(cortex_striatum_corr_cat(:,plot_str,:),3),'k');
xlabel('Cortical MUA aligned depth');
ylabel('Correlation');
legend([p1,p2],{'Cortical fluorescence', ...
    ['Str ' num2str(plot_str) ' multiunit']});
xlim([0,1400]);

% (Fluor-Ctx MUA vs Str-Ctx MUA correlation by depth statistics)
use_str = 1;
use_ctx_str_corr = squeeze(cortex_striatum_corr_cat(:,use_str,:));
fluor_ctx_str_corr = diag(corr(use_ctx_str_corr,cortex_fluor_corr_cat,'rows','complete'));

n_shuff = 10000;
fluor_ctx_str_corr_shuff = nan(size(fluor_ctx_str_corr,1),n_shuff);
for curr_shuff = 1:n_shuff
    use_ctx_str_corr_circshift = use_ctx_str_corr;
    for i = 1:size(use_ctx_str_corr_circshift,2)
        curr_depth_idx = ~isnan(use_ctx_str_corr_circshift(:,i));
        use_ctx_str_corr_circshift(curr_depth_idx,i) = circshift(use_ctx_str_corr_circshift(curr_depth_idx,i), ...
            randi(sum(curr_depth_idx)));
    end   
    fluor_ctx_str_corr_shuff(:,curr_shuff) = ...
        diag(corr(use_ctx_str_corr_circshift, ...
        cortex_fluor_corr_cat,'rows','complete'));
end
corr_rank = tiedrank([nanmean(fluor_ctx_str_corr,1),nanmean(fluor_ctx_str_corr_shuff)]);
corr_p = 1 - (corr_rank(1)/(n_shuff+1));

disp('Fluor-ctx depth vs Str-ctx depth correlation (circshift stat):')
disp(['p = ' num2str(corr_p) ', r = ' num2str(nanmean(fluor_ctx_str_corr)) ...
    ' +/- SEM ' num2str(AP_sem(fluor_ctx_str_corr,1))])

%% EDFig 4e:  Depth fluor/MUA correlation using depth-specific deconv kernel

%%% Load and plot depth-specific deconv kernel
deconv_k_depth_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\gcamp6s_kernel_ctxdepth.mat';
load(deconv_k_depth_fn)

gcamp6s_kernel_cat = fliplr(cat(3,gcamp6s_kernel_ctxdepth.regression{:}));
gcamp6s_kernel_norm = gcamp6s_kernel_cat./max(abs(gcamp6s_kernel_cat),[],2);
gcamp6s_kernel_mean = nanmean(gcamp6s_kernel_norm,3);

depth_col = lines(size(gcamp6s_kernel_mean,1));

figure; 

subplot(1,2,1); hold on;
p = AP_errorfill(gcamp6s_kernel_ctxdepth.regression_t, ...
    permute(nanmean(gcamp6s_kernel_norm,3),[2,1,3]), ...
    permute(nanstd(gcamp6s_kernel_norm,[],3),[2,1,3]), ...
    depth_col);
line(xlim,[0,0],'color','k','linestyle','--');
line([0,0],ylim,'color','k','linestyle','--');
set(p(2),'linestyle','--');

% (depth kernel statistics)
disp('Kernel 2-way ANOVA:')
[depth_grp,t_grp,exp_grp] = ndgrid(1:size(gcamp6s_kernel_norm,1), ...
    1:size(gcamp6s_kernel_norm,2),1:size(gcamp6s_kernel_norm,3));
curr_p = anovan(gcamp6s_kernel_norm(:),[depth_grp(:),t_grp(:)],'display','off');
disp(['depth p = ' num2str(curr_p(1))]);
disp(['time p = ' num2str(curr_p(2))]);

%%% Load and plot correlation data
use_protocol = 'vanillaChoiceworld';

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
data_fn = [data_path filesep 'ctx_fluor_mua_corr_' use_protocol '_ctxdepthkernel'];
load(data_fn);

mua_depth = data(1).cortex_mua_depth{1}; % they're all the same, use 1st
cortex_fluor_corr_cat = cell2mat(permute(horzcat( ...
    data.cortex_fluor_corr),[1,3,2]));

subplot(1,2,2); hold on;
p = AP_errorfill(mua_depth, ...
    nanmean(cortex_fluor_corr_cat,3), ...
    AP_sem(cortex_fluor_corr_cat,3), ...
    depth_col);
xlim([0,1400]);
xlabel('Cortical MUA aligned depth');
ylabel('Correlation');
legend(p,{'Superficial kernel','Deep kernel'});
set(p(2),'linestyle','--');

% (depth kernel correlation statistics)
disp('Depth correlation 2-way ANOVA:')
[depth_grp,kernel_grp,exp_grp] = ndgrid(1:size(cortex_fluor_corr_cat,1), ...
    1:size(cortex_fluor_corr_cat,2),1:size(cortex_fluor_corr_cat,3));
curr_p = anovan(cortex_fluor_corr_cat(:),[depth_grp(:),kernel_grp(:)],'display','off');
disp(['depth p = ' num2str(curr_p(1))]);
disp(['kernel p = ' num2str(curr_p(2))]);



%% EDFig4f: Latency between VisAM and DMS spikes
% (takes a while to run - not really worth preprocessing/saving though)

animals = {'AP043','AP060','AP061'};

% Load ephys alignment
vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
load(vis_ctx_ephys_fn);

lag_t = 1000; % ms
ctx_str_xcorr = struct;
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    protocol = 'vanillaChoiceworld';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        % Load experiment
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(1); % (only one doubled experiment, and the second one was worse)
        site = 2; % cortex probe is always site 2
        str_align = 'none'; % (because this is on cortex)
        load_parts.ephys = true;
        AP_load_experiment
        
        %%%%% CORTEX
        
        % Get depth group boundaries
        % (these are just manual based on the CSD)
        n_aligned_depths = 2;
        ctx_depth_edges = linspace(0,1200,n_aligned_depths+1);
        
        % Depth-align templates
        curr_animal_idx = strcmp(animal,{vis_ctx_ephys.animal});
        curr_day_idx = strcmp(day,vis_ctx_ephys(curr_animal_idx).day);
        curr_csd_depth = vis_ctx_ephys(curr_animal_idx).stim_csd_depth{curr_day_idx};
        curr_csd_depth_aligned = vis_ctx_ephys(curr_animal_idx).stim_csd_depth_aligned{curr_day_idx};
        
        template_depths_aligned = interp1(curr_csd_depth,curr_csd_depth_aligned,template_depths);
        spike_depths_aligned = template_depths_aligned(spike_templates);
        
        % Find cortex end by largest gap between templates
        sorted_template_depths = sort([template_depths_aligned]);
        [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
        ctx_end = sorted_template_depths(max_gap_idx)+1;
        ctx_depth = [sorted_template_depths(1),ctx_end];
        ctx_units = template_depths_aligned <= ctx_depth(2);
        
        % Assign templates to depth groups
        ctx_spike_depths = spike_depths_aligned;
        ctx_spike_depths(spike_depths_aligned < ctx_depth(1) | spike_depths_aligned > ctx_depth(2)) = NaN;
        ctx_depth_group = discretize(ctx_spike_depths,ctx_depth_edges);
        
        % Set binning time
        skip_seconds = 60;
        spike_binning_t = 0.001; % seconds
        spike_binning_t_edges = spike_times_timeline(1)+skip_seconds: ...
            spike_binning_t:spike_times_timeline(end)-skip_seconds;
        spike_binning_t_centers = spike_binning_t_edges(1:end-1) + diff(spike_binning_t_edges)/2;
            
        cortex_mua = nan(size(ctx_depth_group,2),length(spike_binning_t_centers));
        for curr_depth = 1:n_aligned_depths
            curr_spike_times = spike_times_timeline(ctx_depth_group == curr_depth);
            % (skip if no spikes at this depth)
            if isempty(curr_spike_times)
                continue
            end
            cortex_mua(curr_depth,:) = histcounts(curr_spike_times,spike_binning_t_edges);
        end               
        
        
        %%%%% STRIATUM
       
        site = 1; % (striatum is always on probe 1)
        n_aligned_depths = 3;
        str_align = 'kernel';
        AP_load_experiment;
        
        striatum_mua = nan(n_aligned_depths,length(spike_binning_t_centers));
        for curr_depth = 1:n_aligned_depths
            curr_spike_times = spike_times_timeline(aligned_str_depth_group == curr_depth);
            % Skip if no spikes at this depth
            if isempty(curr_spike_times)
                continue
            end
            striatum_mua(curr_depth,:) = histcounts(curr_spike_times,spike_binning_t_edges);
        end
        
        %%% Cortex-striatum correlation
              
        % Get cross-correlations
        curr_ctx_str_xcorr = nan(size(cortex_mua,1),lag_t*2+1);
        for curr_ctx = 1:size(cortex_mua,1)
            curr_ctx_str_xcorr(curr_ctx,:) = ...
                xcorr(striatum_mua(1,:),cortex_mua(curr_ctx,:),lag_t,'coeff');
        end       
        ctx_str_xcorr(curr_animal).ctx_str{curr_day} = curr_ctx_str_xcorr;
        ctx_str_xcorr(curr_animal).ctx_layer{curr_day} = ...
            xcorr(cortex_mua(1,:),cortex_mua(2,:),lag_t,'coeff');
        
        % Get autocorrelation-normalized impulse response (Nauhaus 2012)         
        curr_ctx_str_xcorr_autonorm = nan(size(cortex_mua,1),lag_t*2+1);
        for curr_ctx = 1:size(cortex_mua,1)
            soft_reg_factor = 1e4;
            curr_autonorm = ifft((fft(striatum_mua(1,:)).* ...
                conj(fft(cortex_mua(curr_ctx,:))))./ ...
                (soft_reg_factor+fft(cortex_mua(curr_ctx,:)).* ...
                conj(fft(cortex_mua(curr_ctx,:)))));
            curr_ctx_str_xcorr_autonorm(curr_ctx,:) = ...
                [curr_autonorm(end-lag_t:end),curr_autonorm(1:lag_t)];
        end
        ctx_str_xcorr(curr_animal).ctx_str_autonorm{curr_day} = curr_ctx_str_xcorr_autonorm;
        
        curr_autonorm = ifft((fft(cortex_mua(1,:)).* ...
            conj(fft(cortex_mua(2,:))))./ ...
            (soft_reg_factor+fft(cortex_mua(2,:)).* ...
            conj(fft(cortex_mua(2,:)))));
        ctx_str_xcorr(curr_animal).ctx_layer_autonorm{curr_day} = ...
                [curr_autonorm(end-lag_t:end),curr_autonorm(1:lag_t)];  
        
        % Clear variables for next experiment
        clearvars('-except',preload_vars{:});
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
    
end

% Concatenate (no normalizing?)
ctx_str_cat = cell2mat(permute([ctx_str_xcorr(:).ctx_str],[1,3,2]));
ctx_str_autonorm_cat = cell2mat(permute([ctx_str_xcorr(:).ctx_str_autonorm],[1,3,2]));

ctx_layer_cat = cell2mat([ctx_str_xcorr(:).ctx_layer]');
ctx_layer_autonorm_cat = cell2mat([ctx_str_xcorr(:).ctx_layer_autonorm]');

% Plot
lags = -lag_t:lag_t;

figure;
subplot(2,2,1);
AP_errorfill(lags,nanmean(ctx_str_cat,3)',AP_sem(ctx_str_cat,3)');
line([0,0],ylim,'color','k');
xlabel('Lag (ms)')
ylabel('Correlation (norm)')
title('Deep cortex-DMS xcorr (autocorr norm)');
subplot(2,2,2);
AP_errorfill(lags,nanmean(ctx_str_autonorm_cat,3)',AP_sem(ctx_str_autonorm_cat,3)');
line([0,0],ylim,'color','k');
xlabel('Lag (ms)')
ylabel('Impulse response')
title('Deep cortex-DMS SC''/CC''');
subplot(2,2,3);
AP_errorfill(lags,nanmean(ctx_layer_cat,1)',AP_sem(ctx_layer_cat,1)','k');
line([0,0],ylim,'color','k');
xlabel('Lag (ms)')
ylabel('Correlation (norm)')
title('Cortex deep-superficial xcorr (autocorr norm)');
subplot(2,2,4);
AP_errorfill(lags,nanmean(ctx_layer_autonorm_cat,1)',AP_sem(ctx_layer_autonorm_cat,1)','k');
line([0,0],ylim,'color','k');
xlabel('Lag (ms)')
ylabel('Impulse response')
title('Cortex deep-superficial SC''/CC''');

linkaxes(get(gcf,'Children'),'x');


%% ** EDFig6a: Stimulus activity vs choice

% Get stim activity by stim/choice/depth/experiment
stim_avg_t = [0,0.2];
stim_avg_t_idx = t >= stim_avg_t(1) & t <= stim_avg_t(2);
stim_act = permute(nanmean(mua_allcat(:,stim_avg_t_idx,:),2),[1,3,2]);
stims = unique(trial_stim_allcat);

stim_act_grp = nan(length(stims),2,n_depths,length(use_split));
for curr_depth = 1:n_depths
    for curr_exp = 1:length(use_split)
        for curr_stim = 1:length(stims)
            
            stim_act_grp(curr_stim,1,curr_depth,curr_exp) = ...
                nanmean(stim_act( ...
                move_t < 0.5 & ...
                split_idx == curr_exp & ...
                trial_stim_allcat == stims(curr_stim) & ...
                trial_choice_allcat == -1,curr_depth));
            
            stim_act_grp(curr_stim,2,curr_depth,curr_exp) = ...
                nanmean(stim_act( ...
                move_t < 0.5 & ...
                split_idx == curr_exp & ...
                trial_stim_allcat == stims(curr_stim) & ...
                trial_choice_allcat == 1,curr_depth));
            
        end
    end
end

% Plot stim responses by condition
move_col = [0.6,0,0.6;1,0.6,0];
figure;
for curr_depth = 1:n_depths
   subplot(n_depths,1,curr_depth); hold on;
   set(gca,'ColorOrder',move_col);
    errorbar(repmat(stims,1,2), ...
        nanmean(stim_act_grp(:,:,curr_depth,:),4), ...
        AP_sem(stim_act_grp(:,:,curr_depth,:),4),'linewidth',2);
    xlabel('Contrast*side');
    ylabel(['Str ' num2str(curr_depth)']);
end

% (stim response by choice statistics)
disp('Stim v choice 2-way anova:')
for curr_depth = 1:n_depths
    [stim_grp,choice_grp,exp_grp] = ndgrid(1:size(stim_act_grp,1), ...
        1:size(stim_act_grp,2),1:size(stim_act_grp,4));
    [curr_p,~,~,terms] = anovan(reshape(stim_act_grp(:,:,curr_depth,:),[],1), ...
        [stim_grp(:),choice_grp(:)],'model','interaction','display','off');
    disp(['Str ' num2str(curr_depth) ' choice: p = ' num2str(curr_p(2)')]);
    disp(['Str ' num2str(curr_depth) ' stim x choice: p = ' num2str(curr_p(2)')]); 
end


%% ** EDFig6b,c: Task > cortex kernels (go cue)

% Get task>cortex parameters
n_regressors = length(task_regressor_labels);
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

% Concatenate and average task>cortex kernels
fluor_taskpred_k_cat = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(k) k{regressor}, ...
    vertcat(fluor_taskpred_k_all{:}),'uni',false),[2,3,4,1])), ...
    1:n_regressors,'uni',false);

fluor_taskpred_k_avg = cellfun(@(x) nanmean(x,4), fluor_taskpred_k_cat,'uni',false);
fluor_taskpred_k_avg_px = cellfun(@(x) ...
    AP_svdFrameReconstruct(U_master(:,:,1:n_vs),permute(x,[3,2,1])), ...
    fluor_taskpred_k_avg,'uni',false);

% Plot kernel max across time
max_subregressors = max(cellfun(@(x) size(x,1),fluor_taskpred_k_avg));

figure;
for curr_regressor = 1:length(task_regressor_labels)

    c_max = max(fluor_taskpred_k_avg_px{curr_regressor}(:));  
    for curr_subregressor = 1:size(fluor_taskpred_k_avg_px{curr_regressor},4)
        subplot(length(task_regressor_labels),max_subregressors, ...
            (curr_regressor-1)*max_subregressors+curr_subregressor);
        imagesc(max(fluor_taskpred_k_avg_px{curr_regressor}(:,:,:,curr_subregressor),[],3));
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        axis image off;
        colormap(brewermap([],'Greys'));
        caxis([0,c_max]);
        title([task_regressor_labels{curr_regressor}]);       
    end
    
end

%% ** EDFig7d: Task explained variance cortex vs striatum

% Use raw data (not normalized or baseline-subtracted) for expl var
mua_exp = vertcat(mua_all{:});
mua_taskpred_exp = vertcat(mua_taskpred_all{:});
mua_ctxpred_exp = vertcat(mua_ctxpred_all{:});
mua_ctxpred_taskpred_exp = vertcat(mua_ctxpred_taskpred_all{:});

% Get R^2 for task and cortex
task_str_r2 = nan(max(split_idx),n_depths);
task_ctx_r2 = nan(max(split_idx),n_depths);
for curr_exp = 1:max(split_idx)
       
    curr_str_data_baselinesub = reshape(permute(mua_exp{curr_exp},[2,1,3]),[],n_depths) - ...
        (nanmean(reshape(mua_exp{curr_exp}(:,t < 0,:),[],size(mua_exp{curr_exp},3)),1));
    curr_str_taskpred_data = reshape(permute(mua_taskpred_exp{curr_exp},[2,1,3]),[],n_depths);
    
    curr_ctx_data_baselinesub = reshape(permute(mua_ctxpred_exp{curr_exp},[2,1,3]),[],n_depths) - ...
        (nanmean(reshape(mua_ctxpred_exp{curr_exp}(:,t < 0,:),[],size(mua_ctxpred_exp{curr_exp},3)),1));
    curr_ctx_taskpred_data = reshape(permute(mua_ctxpred_taskpred_exp{curr_exp},[2,1,3]),[],n_depths);
       
    % Set common NaNs
    nan_samples = isnan(curr_str_data_baselinesub) | isnan(curr_str_taskpred_data) | ...
        isnan(curr_ctx_data_baselinesub) | isnan(curr_ctx_taskpred_data);
    curr_str_data_baselinesub(nan_samples) = NaN;
    curr_str_taskpred_data(nan_samples) = NaN;
    curr_ctx_data_baselinesub(nan_samples) = NaN;
    curr_ctx_taskpred_data(nan_samples) = NaN;

    task_str_r2(curr_exp,:) = 1 - (nansum((curr_str_data_baselinesub-curr_str_taskpred_data).^2,1)./ ...
        nansum((curr_str_data_baselinesub-nanmean(curr_str_data_baselinesub,1)).^2,1));
    task_ctx_r2(curr_exp,:) = 1 - (nansum((curr_ctx_data_baselinesub-curr_ctx_taskpred_data).^2,1)./ ...
        nansum((curr_ctx_data_baselinesub-nanmean(curr_ctx_data_baselinesub,1)).^2,1));
    
end

figure; hold on;
errorbar(nanmean(task_str_r2,1),AP_sem(task_str_r2,1),'b','linewidth',2,'CapSize',0);
errorbar(nanmean(task_ctx_r2,1),AP_sem(task_ctx_r2,1),'color',[0,0.6,0],'linewidth',2,'CapSize',0);
xlabel('Striatum depth');
ylabel('Task explained variance');
legend({'Task','Cortex'});

% Plot explained variance task vs cortex by experiment
figure; hold on;
str_col = max(hsv(n_depths)-0.2,0);
for curr_str = 1:n_depths
    errorbar(squeeze(nanmean(task_str_r2(:,curr_str),1)), ...
        squeeze(nanmean(task_ctx_r2(:,curr_str),1)), ...
        squeeze(AP_sem(task_ctx_r2(:,curr_str),1)),squeeze(AP_sem(task_ctx_r2(:,curr_str),1)), ...
        squeeze(AP_sem(task_str_r2(:,curr_str),1)),squeeze(AP_sem(task_str_r2(:,curr_str),1)), ...
        'color','k','linewidth',2);
    
    scatter(task_str_r2(:,curr_str),task_ctx_r2(:,curr_str),10, ...
        str_col(curr_str,:),'filled');
    scatter(nanmean(task_str_r2(:,curr_str),1), ...
        nanmean(task_ctx_r2(:,curr_str),1),80, ...
        str_col(curr_str,:),'filled','MarkerEdgeColor','k','linewidth',2);
end
axis tight equal;
line(xlim,xlim,'color','k','linestyle','--');
xlabel('Task striatum R^2');
ylabel('Task cortex R^2');
legend({'DMS','DCS','DLS'})

% (Task-striatum R2 statistics)
disp('Task R^2 values:');
for curr_depth = 1:n_depths
    disp(['Str ' num2str(curr_depth) ...
        ' R2 = ' num2str(nanmean(task_str_r2(:,curr_depth))) '+/- ' ...
        num2str(AP_sem(task_str_r2(:,curr_depth),1))]); 
end

% (Task-cortex R2 statistics)
disp('Cortex R^2 values:');
for curr_depth = 1:n_depths
    disp(['Str ' num2str(curr_depth) ...
        ' R2 = ' num2str(nanmean(task_ctx_r2(:,curr_depth))) '+/- ' ...
        num2str(AP_sem(task_ctx_r2(:,curr_depth),1))]); 
end

% (Cortex vs striatum task R2 statistics)
disp('Cortex vs Task R^2 signrank:');
for curr_depth = 1:n_depths
    curr_p = signrank(task_ctx_r2(:,curr_depth), ...
        task_str_r2(:,curr_depth));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end

disp('Cortex vs Task R^2 correlation:');
[r,p] = corr(task_str_r2(:,curr_depth), ...
    task_ctx_r2(:,curr_depth), ...
    'rows','complete','type','pearson');
disp(['r = ' num2str(r) ' p = ' num2str(p)]);



%% ** EDFig8a: Striatum prediction: cortical subregions, striatal domains
% NOTE: this regression is done on trial data rather than the long time
% courses which is what the normal analysis uses. For sanity check, the
% explained using the full data (full) and the trials dataset (trials) are
% compared here but not meant to be included.
%
% Also this takes a long time to run

% Choose depths to run
plot_depth = 1:n_depths;

% Regress kernel ROI activity to striatum domain activity (per recording)
regression_params.use_svs = 1:100;
regression_params.kernel_t = [-0.1,0.1];
regression_params.zs = [false,false];
regression_params.cvfold = 2;
regression_params.use_constant = true;
lambda = 50; % lambda for fluorescence to striatum
domain_lambda = 0; % lambda for other domains to striatum
kernel_frames = floor(regression_params.kernel_t(1)*sample_rate): ...
    ceil(regression_params.kernel_t(2)*sample_rate);

% Set time to use
% (e.g. this type of code can be used to regress only ITI time)
use_t = true(size(t));

% Set regions to zero
% (zero all EXCEPT quadrants)
ctx_zero = true(size(U_master,1),size(U_master,2),6);
ctx_zero(1:260,1:220,1) = false;
ctx_zero(1:260,220:end,2) = false;
ctx_zero(260:end,1:220,3) = false;
ctx_zero(260:end,220:end,4) = false;
ctx_zero(:,220:end,5) = false;
ctx_zero(:,1:220,6) = false;

% (use raw data for trial regression)
mua_exp = vertcat(mua_all{:});
fluor_exp = vertcat(fluor_all{:});
fluor_deconv_exp = cellfun(@AP_deconv_wf,fluor_exp,'uni',false);

mua_ctxtrialpred_exp = cellfun(@(x) nan(size(x)),mua_exp,'uni',false);
mua_ctxtrialpred_k = nan(length(regression_params.use_svs),length(kernel_frames),n_depths,length(use_split));
mua_ctxtrialpred_regionzero_exp = cellfun(@(x) nan([size(x),size(ctx_zero,3)]),mua_exp,'uni',false);
mua_ctxtrialpred_regionzero_k = nan(length(regression_params.use_svs),length(kernel_frames),n_depths,length(use_split));

mua_domaintrialpred_exp = cellfun(@(x) nan(size(x)),mua_exp,'uni',false);
mua_domaintrialpred_k = nan(n_depths-1,length(kernel_frames),n_depths,length(use_split));

for curr_exp = 1:length(trials_recording)
    for curr_depth = plot_depth
        
        curr_mua = reshape(mua_exp{curr_exp}(:,use_t,curr_depth)',1,[]);
        curr_mua_std = curr_mua./nanstd(curr_mua);
        
        curr_fluor = reshape(permute( ...
            fluor_deconv_exp{curr_exp}(:,use_t,:),[2,1,3]),[],n_vs)';
        curr_fluor_full = reshape(permute( ...
            fluor_deconv_exp{curr_exp}(:,:,:),[2,1,3]),[],n_vs)';
        
        curr_fluor_regionzero = nan(size(curr_fluor,1),size(curr_fluor,2),size(ctx_zero,3));
        for curr_zero = 1:size(ctx_zero,3)
            U_master_regionzero = U_master.*~ctx_zero(:,:,curr_zero);
            fVdf_regionzero_altU = ChangeU(U_master(:,:,1:n_vs),curr_fluor,U_master_regionzero(:,:,1:n_vs));
            % (those U's aren't orthonormal, recast back to original Udf_aligned)
            fVdf_regionzero = ChangeU(U_master_regionzero(:,:,1:n_vs),fVdf_regionzero_altU,U_master(:,:,1:n_vs));
            curr_fluor_regionzero(:,:,curr_zero) = fVdf_regionzero;
        end
        
        % Skip if no data
        if all(isnan(curr_mua))
            continue
        end
        
        % Set discontinuities in trial data
        trial_discontinuities = false(size(mua_exp{curr_exp}(:,use_t,curr_depth)));
        trial_discontinuities(:,1) = true;
        trial_discontinuities = reshape(trial_discontinuities',[],1)';
        
        % Full cortex regression
        [ctx_str_k,curr_mua_fluorpred_std,explained_var_trial] = ...
            AP_regresskernel(curr_fluor(regression_params.use_svs,:), ...
            curr_mua_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant,trial_discontinuities);
                
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        ctxpred_spikes = (curr_mua_fluorpred_std - squeeze(ctx_str_k{end})).* ...
            nanstd(curr_mua,[],2) + ...
            nanstd(curr_mua,[],2).*squeeze(ctx_str_k{end});
        
        mua_ctxtrialpred_k(:,:,curr_depth,curr_exp) = ctx_str_k{1};
        mua_ctxtrialpred_exp{curr_exp}(:,use_t,curr_depth) = ...
            reshape(ctxpred_spikes,sum(use_t),[])';
        %         % (apply kernel to full time)
        %         mua_ctxtrialpred_exp{curr_exp}(:,:,curr_depth) = ...
        %             sum(cell2mat(arrayfun(@(x) ...
        %             convn(fluor_allcat_deconv_exp{curr_exp}(:,:,regression_params.use_svs(x)), ...
        %             k_fluor(x,:)','same'),permute(1:length(regression_params.use_svs),[1,3,2]),'uni',false)),3);
        
        % Region-zeroed cortex regression
        for curr_zero = 1:size(ctx_zero,3)
            [ctx_str_k_regionzero,curr_mua_fluorregionzeropred_std,explained_var_trial_regionzero] = ...
                AP_regresskernel(curr_fluor_regionzero(regression_params.use_svs,:,curr_zero), ...
                curr_mua_std,kernel_frames,lambda, ...
                regression_params.zs,regression_params.cvfold, ...
                true,regression_params.use_constant,trial_discontinuities);
            
            % Re-scale the prediction (subtract offset, multiply, add scaled offset)
            ctxpred_spikes_regionzero = (curr_mua_fluorregionzeropred_std - squeeze(ctx_str_k_regionzero{end})).* ...
                nanstd(curr_mua,[],2) + ...
                nanstd(curr_mua,[],2).*squeeze(ctx_str_k_regionzero{end});
                        
            mua_ctxtrialpred_regionzero_k(:,:,curr_depth,curr_exp,curr_zero) = ctx_str_k_regionzero{1};
            mua_ctxtrialpred_regionzero_exp{curr_exp}(:,use_t,curr_depth,curr_zero) = ...
                reshape(ctxpred_spikes_regionzero,sum(use_t),[])';
            %         % (apply kernel to full time)
            %         mua_ctxtrialpred_regionzero_exp{curr_exp}(:,:,curr_depth) = ...
            %             sum(cell2mat(arrayfun(@(x) ...
            %             convn(fluor_allcat_deconv_exp{curr_exp}(:,:,regression_params.use_svs(x)), ...
            %             k_fluorregionzero(x,:)','same'),permute(1:length(regression_params.use_svs),[1,3,2]),'uni',false)),3);
        end
        
        % Regress current domains from other domains
        curr_mua_domains = reshape(permute(mua_exp{curr_exp}(:,use_t,:),[2,1,3]),[],n_depths)';
        curr_mua_domains_std = curr_mua_domains./nanstd(curr_mua_domains,[],2);
        
        [domain_str_k,curr_mua_domainpred_std,explained_var_trial] = ...
            AP_regresskernel(curr_mua_domains_std(setdiff(1:n_depths,curr_depth),:), ...
            curr_mua_std,kernel_frames,domain_lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant,trial_discontinuities);
        
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        curr_mua_domainpred = (curr_mua_domainpred_std - squeeze(domain_str_k{end})).* ...
            nanstd(curr_mua,[],2) + ...
            nanstd(curr_mua,[],2).*squeeze(domain_str_k{end});
        
        mua_domaintrialpred_k(:,:,curr_depth,curr_exp) = domain_str_k{1};
        mua_domaintrialpred_exp{curr_exp}(:,use_t,curr_depth) = ...
            reshape(curr_mua_domainpred,sum(use_t),[])';
      
        
    end
    AP_print_progress_fraction(curr_exp,length(trials_recording));
end

% Plot average cortex->striatum kernel
ctx_str_k_mean = nanmean(cell2mat(permute(vertcat(ctx_str_k_all{:}),[2,3,4,5,1])),5);
ctx_str_k_mean_px = cell2mat(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:100), ...
    ctx_str_k_mean(:,:,x)),permute(1:n_depths,[1,3,4,2]),'uni',false));

mua_ctxtrialpred_k_mean = nanmean(mua_ctxtrialpred_k,4);
mua_ctxtrialpred_k_mean_px = cell2mat(arrayfun(@(x) ...
    svdFrameReconstruct(U_master(:,:,regression_params.use_svs), ...
    mua_ctxtrialpred_k_mean(:,:,x)),permute(1:n_depths,[1,3,4,2]),'uni',false));

mua_ctxtrialpred_regionzero_k_mean = nanmean(mua_ctxtrialpred_regionzero_k,4);
mua_ctxtrialpred_regionzero_k_mean_px = cell2mat(arrayfun(@(x) ...
    svdFrameReconstruct(U_master(:,:,regression_params.use_svs), ...
    reshape(mua_ctxtrialpred_regionzero_k_mean(:,:,x,:),length(regression_params.use_svs),[])), ...
    permute(1:n_depths,[1,3,4,2]),'uni',false));

AP_image_scroll(cat(3,mua_ctxtrialpred_k_mean_px,mua_ctxtrialpred_regionzero_k_mean_px));
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(gca,brewermap([],'PRGn'));
axis image;

% Plot average domain->striautm kernel
figure;
mua_domaintrialpred_k_mean = nanmean(mua_domaintrialpred_k,4);
str_col = max(hsv(n_depths)-0.2,0);
for curr_depth = 1:n_depths
    subplot(n_depths,1,curr_depth); hold on
    set(gca,'ColorOrder',str_col(setdiff(1:n_depths,curr_depth),:));
    plot(kernel_frames,mua_domaintrialpred_k_mean(:,:,curr_depth)','linewidth',2);
    ylabel('Weight');
    title(['Str ' num2str(curr_depth)]);
end


% Get R^2 for task, cortex full, and cortex ROI predictions
mua_ctxpred_exp = vertcat(mua_ctxpred_all{:});

ctxpred_r2 = nan(max(split_idx),n_depths);
ctxtrialpred_r2 = nan(max(split_idx),n_depths);
ctxtrialpred_regionzero_r2 = nan(max(split_idx),n_depths,size(ctx_zero,3));
domaintrialpred_r2 = nan(max(split_idx),n_depths);

for curr_exp = 1:max(split_idx)
    
    curr_data = reshape(permute(mua_exp{curr_exp},[2,1,3]),[],n_depths);
    curr_ctxpred_data = reshape(permute(mua_ctxpred_exp{curr_exp},[2,1,3]),[],n_depths);
    curr_ctxtrialpred_data = reshape(permute(mua_ctxtrialpred_exp{curr_exp},[2,1,3]),[],n_depths);
    curr_ctxtrialpred_regionzero_data = ...
        reshape(permute(mua_ctxtrialpred_regionzero_exp{curr_exp},[2,1,3,4]),[],n_depths,size(ctx_zero,3));
    curr_domaintrialpred_data = reshape(permute(mua_domaintrialpred_exp{curr_exp},[2,1,3]),[],n_depths);
    
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_ctxpred_data) ...
        | isnan(curr_ctxtrialpred_data) | any(isnan(curr_ctxtrialpred_regionzero_data),3) ...
        | isnan(curr_domaintrialpred_data);
    curr_data(nan_samples) = NaN;
    curr_ctxpred_data(nan_samples) = NaN;
    curr_ctxtrialpred_data(nan_samples) = NaN;
    curr_ctxtrialpred_regionzero_data(repmat(nan_samples,1,1,size(ctx_zero,3))) = NaN;
    curr_domaintrialpred_data(nan_samples) = NaN;
    
    ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    ctxtrialpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxtrialpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    ctxtrialpred_regionzero_r2(curr_exp,:,:) = 1 - (nansum((curr_data-curr_ctxtrialpred_regionzero_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    domaintrialpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_domaintrialpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    
end

%%% Cortex subregions explained variance
figure;

% Plot full vs trial (sanity check: they should be ~the same)
subplot(2,3,1); hold on;
errorbar(nanmean(ctxpred_r2,1),AP_sem(ctxpred_r2,1),'.-','linewidth',2,'CapSize',0,'MarkerSize',30);
errorbar(nanmean(ctxtrialpred_r2,1),AP_sem(ctxtrialpred_r2,1),'.-','linewidth',2,'CapSize',0,'MarkerSize',30);
legend({'Ctx full','Ctx trial'});
xlabel('Striatum depth');
ylabel('Explained variance');

% Plot explained variance by subregion
subplot(2,3,2); hold on;
str_col = max(hsv(n_depths)-0.2,0);
set(gca,'ColorOrder',str_col);
errorbar(permute(nanmean(ctxtrialpred_regionzero_r2,1),[3,2,1]), ...
    permute(AP_sem(ctxtrialpred_regionzero_r2,1),[3,2,1]),'.-','linewidth',2,'CapSize',0,'MarkerSize',30);
xlabel('Cortex subregion');
ylabel('Explained variance');
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false));

% Plot explained variance of subregion and domain relative to full
alternate_r2 = cat(3,ctxtrialpred_regionzero_r2,domaintrialpred_r2);
alternate_r2_relative = ...
    (alternate_r2 - ctxtrialpred_r2)./ctxtrialpred_r2;

subplot(2,3,3); hold on;
str_col = max(hsv(n_depths)-0.2,0);
set(gca,'ColorOrder',str_col);
errorbar(permute(nanmean(alternate_r2_relative,1),[3,2,1]), ...
    permute(AP_sem(alternate_r2_relative,1),[3,2,1]),'linewidth',2,'CapSize',0,'Marker','none');
set(gca,'XTick',1:size(alternate_r2,3),'XTickLabel', ...
    [cellfun(@(x) ['Ctx ' num2str(x)],num2cell(1:size(ctx_zero,3)),'uni',false), ...
    'Domains']);
line(xlim,[0,0],'color','k','linestyle','--');
xlabel('Alternate regressor');
ylabel('\DeltaR^2/R^2_{full}')
legend(cellfun(@(x) ['Str ' num2str(x)],num2cell(1:n_depths),'uni',false));

% Plot subregions used
for i = 1:size(ctx_zero,3)
   subplot(2,size(ctx_zero,3),size(ctx_zero,3)+i);
   imagesc(~ctx_zero(:,:,i));
   AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
   axis image off;
   colormap(gray);
end


% (Prediction with cortex subsets statistics)
disp('R^2 with cortical subset (2-way anova):');
for curr_depth = 1:n_depths   
    curr_r2 = permute(ctxtrialpred_regionzero_r2_relative(:,curr_depth,:),[3,1,2]);
    [condition_grp,exp_grp] = ndgrid(1:size(curr_r2,1),1:size(curr_r2,2));
    curr_p = anovan(curr_r2(:),condition_grp(:),'display','off');  
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p')]);
end

% (Prediction with domains statistics)
disp('Domain vs. cortex explained variance (signrank):');
for curr_depth = 1:n_depths   
    curr_p = signrank(domaintrialpred_r2_relative(:,curr_depth)); 
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p')]);
end

%% EDFig8b: Passive fluor+ctx+str example

%%% Plot example day

animal = 'AP060';
day = '2019-12-06';
experiment = 2;
plot_t = [100,300];

figure;
disp('Loading example data...');

% Load cortex ephys + imaging
load_parts.ephys = true;
load_parts.imaging = true;
site = 2; % (cortex always probe 2)
str_align = 'none'; % (cortex)
AP_load_experiment;

% Load cortex recording alignment
vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
load(vis_ctx_ephys_fn);

%%% DEPTH-ALIGN TEMPLATES, FIND CORTEX BOUNDARY
curr_animal_idx = strcmp(animal,{vis_ctx_ephys.animal});
curr_day_idx = strcmp(day,vis_ctx_ephys(curr_animal_idx).day);
curr_csd_depth = vis_ctx_ephys(curr_animal_idx).stim_csd_depth{curr_day_idx};
curr_csd_depth_aligned = vis_ctx_ephys(curr_animal_idx).stim_csd_depth_aligned{curr_day_idx};
template_depths_aligned = interp1(curr_csd_depth,curr_csd_depth_aligned,template_depths);
spike_depths_aligned = interp1(curr_csd_depth,curr_csd_depth_aligned,spike_depths);

% Find cortex end by largest gap between templates
sorted_template_depths = sort([template_depths_aligned]);
[max_gap,max_gap_idx] = max(diff(sorted_template_depths));
ctx_end = sorted_template_depths(max_gap_idx)+1;

ctx_depth = [sorted_template_depths(1),ctx_end];
ctx_units = template_depths_aligned <= ctx_depth(2);

%%% GET FLUORESCENCE AND SPIKES BY DEPTH

% Set binning time
skip_seconds = 60;
spike_binning_t = 1/framerate; % seconds
spike_binning_t_edges = frame_t(1)+skip_seconds:spike_binning_t:frame_t(end)-skip_seconds;
spike_binning_t_centers = spike_binning_t_edges(1:end-1) + diff(spike_binning_t_edges)/2;

% Get fluorescence in pre-drawn ROI
curr_ctx_roi = vis_ctx_ephys(curr_animal_idx).ctx_roi{curr_day_idx};

fVdf_deconv = AP_deconv_wf(fVdf);
fluor_roi = AP_svd_roi(Udf,fVdf_deconv,avg_im,[],curr_ctx_roi);
fluor_roi_interp = interp1(frame_t,fluor_roi,spike_binning_t_centers);

% Plot cortex raster
subplot(4,1,1); hold on;
plot_t_idx = spike_binning_t_centers >= plot_t(1) & ...
    spike_binning_t_centers <= plot_t(2);
plot(spike_binning_t_centers(plot_t_idx), ...
    fluor_roi_interp(plot_t_idx),'linewidth',2,'color',[0,0.7,0]);

subplot(4,1,2,'YDir','reverse'); hold on;
plot_spikes = spike_times_timeline >= plot_t(1) & ...
    spike_times_timeline <= plot_t(2) & ...
    spike_depths_aligned <= ctx_depth(2);
plot(spike_times_timeline(plot_spikes),spike_depths(plot_spikes),'.k');
ylabel('Cortex depth (\mum)');
xlabel('Time (s)');

%%% LOAD STRIATUM EPHYS AND GET MUA BY DEPTH

% Load striatum ephys
load_parts.ephys = true;
load_parts.imaging = false;
site = 1; % (striatum is always on probe 1)
str_align = 'kernel';
AP_load_experiment;

% Plot striatum raster
subplot(4,1,3,'YDir','reverse'); hold on;
plot_spikes = spike_times_timeline >= plot_t(1) & ...
    spike_times_timeline <= plot_t(2) & ...
    spike_depths >= str_depth(1) & spike_depths <= str_depth(2);
plot(spike_times_timeline(plot_spikes),spike_depths(plot_spikes),'.k');
ylabel('Striatum depth (\mum)');
xlabel('Time (s)');

%%%%%% PLOT WHEEL/STIM

% (wheel velocity)
wheel_axes = subplot(4,1,4); hold on;
plot_wheel_idx = Timeline.rawDAQTimestamps >= plot_t(1) & ...
    Timeline.rawDAQTimestamps <= plot_t(2);
plot(wheel_axes,Timeline.rawDAQTimestamps(plot_wheel_idx), ...
    wheel_velocity(plot_wheel_idx),'k','linewidth',2);
ylabel('Wheel velocity');
axis off

curr_axes = flipud(get(gcf,'Children'));
% Link all time axes
linkaxes(curr_axes,'x');
% Link depth axes of raster plots (arbitrary depth, but want same scale)
linkaxes(curr_axes(2:3),'xy');


%% EDFig8c-e: Cortex-explained variance task vs. passive

% Load data
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\str_ctxpred.mat')

% Concatenate data across cohorts and animals
mua_animal = [str_ctxpred.str]';
mua_exp = vertcat(mua_animal{:});

mua_ctxpred_animal = [str_ctxpred.str_ctxpred]';
mua_ctxpred_exp = vertcat(mua_ctxpred_animal{:});

ctxpred_r2 = cellfun(@(mua,mua_ctxpred) ...
    1 - (nansum((mua-mua_ctxpred).^2,2)./(nansum((mua-nanmean(mua,2)).^2,2))), ...
    mua_exp,mua_ctxpred_exp,'uni',false);

ctxpred_r2_task = horzcat(ctxpred_r2{:,1})';
ctxpred_r2_passive = horzcat(ctxpred_r2{:,2})';

figure;
n_depths = size(ctxpred_r2_task,2);
str_col = max(hsv(n_depths)-0.2,0);

% Plot cortex-explained variance task vs. passive by experiment
subplot(1,3,1); hold on;
for curr_str = 1:n_depths
    errorbar(squeeze(nanmean(ctxpred_r2_passive(:,curr_str),1)), ...
        squeeze(nanmean(ctxpred_r2_task(:,curr_str),1)), ...
        squeeze(AP_sem(ctxpred_r2_task(:,curr_str),1)), ...
        squeeze(AP_sem(ctxpred_r2_task(:,curr_str),1)), ...
        squeeze(AP_sem(ctxpred_r2_passive(:,curr_str),1)), ...
        squeeze(AP_sem(ctxpred_r2_passive(:,curr_str),1)), ...
        'color','k','linewidth',2);
    
    scatter(ctxpred_r2_passive(:,curr_str),ctxpred_r2_task(:,curr_str),10, ...
        str_col(curr_str,:),'filled');
    scatter(nanmean(ctxpred_r2_passive(:,curr_str),1), ...
        nanmean(ctxpred_r2_task(:,curr_str),1),80, ...
        str_col(curr_str,:),'filled','MarkerEdgeColor','k','linewidth',2);
end
axis tight equal
line(xlim,xlim,'color','k','linestyle','--');
xlabel('Cortex passive R^2');
ylabel('Cortex task R^2');
legend({'DMS','DCS','DLS'})

% Get and plot striatal variance by task vs. passive
mua_var_exp = cellfun(@(x) var(x,[],2),mua_exp,'uni',false);
mua_var_task = log10(horzcat(mua_var_exp{:,1})');
mua_var_passive = log10(horzcat(mua_var_exp{:,2})');

subplot(1,3,2); hold on;
for curr_str = 1:n_depths
    errorbar(squeeze(nanmean(mua_var_passive(:,curr_str),1)), ...
        squeeze(nanmean(mua_var_task(:,curr_str),1)), ...
        squeeze(AP_sem(mua_var_task(:,curr_str),1)), ...
        squeeze(AP_sem(mua_var_task(:,curr_str),1)), ...
        squeeze(AP_sem(mua_var_passive(:,curr_str),1)), ...
        squeeze(AP_sem(mua_var_passive(:,curr_str),1)), ...
        'color','k','linewidth',2);
    
    scatter(mua_var_passive(:,curr_str),mua_var_task(:,curr_str),10, ...
        str_col(curr_str,:),'filled');
    scatter(nanmean(mua_var_passive(:,curr_str),1), ...
        nanmean(mua_var_task(:,curr_str),1),80, ...
        str_col(curr_str,:),'filled','MarkerEdgeColor','k','linewidth',2);
end
axis tight equal
line(xlim,xlim,'color','k','linestyle','--');
xlabel('log(Passive variance)');
ylabel('log(Task variance)');
legend({'DMS','DCS','DLS'})

% Plot cortex R2 vs striatal variance
subplot(1,3,3); hold on;
for curr_str = 1:n_depths
    errorbar(squeeze(nanmean(mua_var_task(:,curr_str),1)), ...
        squeeze(nanmean(ctxpred_r2_task(:,curr_str),1)), ...
        squeeze(AP_sem(ctxpred_r2_task(:,curr_str),1)), ...
        squeeze(AP_sem(ctxpred_r2_task(:,curr_str),1)), ...
        squeeze(AP_sem(mua_var_task(:,curr_str),1)), ...
        squeeze(AP_sem(mua_var_task(:,curr_str),1)), ...
        'color','k','linewidth',2);
    
    scatter(mua_var_task(:,curr_str),ctxpred_r2_task(:,curr_str),10, ...
        str_col(curr_str,:),'filled');
    scatter(nanmean(mua_var_task(:,curr_str),1), ...
        nanmean(ctxpred_r2_task(:,curr_str),1),80, ...
        str_col(curr_str,:),'filled','MarkerEdgeColor','k','linewidth',2);
end
axis square;
xlabel('log(Task variance)');
ylabel('Cortex task R^2');
legend({'DMS','DCS','DLS'})


% (Cortex-passive R2 statistics)
disp('Cortex (passive) R^2 values:');
for curr_depth = 1:n_depths
    disp(['Str ' num2str(curr_depth) ...
        ' R2 = ' num2str(nanmean(ctxpred_r2_passive(:,curr_depth))) '+/- ' ...
        num2str(AP_sem(ctxpred_r2_passive(:,curr_depth),1))]); 
end

% (Cortex-task R2 statistics)
disp('Cortex (task) R^2 values:');
for curr_depth = 1:n_depths
    disp(['Str ' num2str(curr_depth) ...
        ' R2 = ' num2str(nanmean(ctxpred_r2_task(:,curr_depth))) '+/- ' ...
        num2str(AP_sem(ctxpred_r2_task(:,curr_depth),1))]); 
end

% (Cortex vs task R2 statistics)
disp('Cortex task vs passive R^2 signrank:');
for curr_depth = 1:n_depths
    curr_p = signrank(ctxpred_r2_task(:,curr_depth), ...
        ctxpred_r2_passive(:,curr_depth));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end

% (Cortex vs task variance statistics)
disp('Cortex task vs variance R^2 signrank:');
for curr_depth = 1:n_depths
    curr_p = signrank(mua_var_task(:,curr_depth), ...
        mua_var_passive(:,curr_depth));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]); 
end

% (R2 vs variance statistics)
[exp_grp,domain_grp] = ndgrid(1:size(mua_var_task,1),1:size(mua_var_task,2));
use_points = ~isnan(mua_var_task);
[h,atab,ctab,stats] = aoctool(mua_var_task(use_points), ...
    ctxpred_r2_task(use_points),domain_grp(use_points));



%% EDFig9a Cortical spike rate pre/post muscimol
% (use spike rate over all experiments pre/post muscimol)

animal_days = { ...
    'AP052','2019-09-20';
    'AP058','2019-12-06'};

figure;

spike_rate_change_cond = cell(size(animal_days));
for curr_animalday = 1:length(animal_days)
    
    animal = animal_days{curr_animalday,1};
    day = animal_days{curr_animalday,2};
    
    % Load data (first experiment - but spikes throughout used)
    experiment = 1;
    AP_load_experiment
    
    % Estimate boundaries of cortex (the dumb way: first template/gap)
    sorted_template_depths = sort(template_depths);
    ctx_start = sorted_template_depths(1) - 1;
    [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
    ctx_end = sorted_template_depths(max_gap_idx)+1;
    ctx_depth = [ctx_start,ctx_end];
    ctx_templates = template_depths <= ctx_depth(2);
    
    % Set experiments in conditions (1-2 = pre-muscimol, 3-4 = post-muscimol)
    cond_expts = {[1,2],[3,4]};
    
    spike_rate_cond = nan(max(spike_templates),2);
    for curr_cond = 1:2
        
        exp_starts = sync(2).timestamps(sync(2).values == 1);
        exp_stops = sync(2).timestamps(sync(2).values == 0);
        
        curr_exp_start = exp_starts(cond_expts{curr_cond}(1));
        curr_exp_stop = exp_stops(cond_expts{curr_cond}(end));
        
        curr_use_spikes = spike_times >= curr_exp_start & ...
            spike_times <= curr_exp_stop;
        
        spike_rate_cond(:,curr_cond) = ...
            accumarray(spike_templates(curr_use_spikes),1,[max(spike_templates),1])./ ...
            (curr_exp_stop - curr_exp_start);
        
    end
    
    spike_rate_change = (spike_rate_cond(:,2) - spike_rate_cond(:,1))./(spike_rate_cond(:,1)+spike_rate_cond(:,2));
    
    subplot(length(animal_days),2,(curr_animalday-1)*2+1,'YDir','reverse'); hold on;
    plot(reshape([spike_rate_cond,nan(size(template_depths))]',[],1), ...
        reshape(repmat(template_depths,1,3)',[],1),'color',[0.5,0.5,0.5]);
    p1 = plot(spike_rate_cond(:,1),template_depths,'.k','MarkerSize',10);
    p2 = plot(spike_rate_cond(:,2),template_depths,'.r','MarkerSize',10);
    xlabel('Spikes/s')
    ylabel('Depth (\mum)');
    legend([p1,p2],{'Pre-muscimol','Post-muscimol'});
    axis tight;
    xlim([-1,prctile(spike_rate_cond(:,1),95)])
   
    subplot(length(animal_days),2,(curr_animalday-1)*2+2); hold on;
    plot(spike_rate_change,template_depths,'.k','MarkerSize',10);   
    axis tight;
    xlim([-1.1,1.1]);
    line([0,0],ylim);
    set(gca,'YDir','reverse');
    xlabel('(Post-pre)/(pre+post)');
    ylabel('Depth (\mum)');
    title({animal,day,'Change'});
    
    
end



%% EDFig9b VFS pre/post musicmol

muscimol_wf_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\muscimol_wf.mat';
load(muscimol_wf_fn);

std_cat = vertcat(muscimol_wf.std);
vfs_cat = vertcat(muscimol_wf.vfs);

std_change = (cat(3,std_cat{:,2})-cat(3,std_cat{:,1}))./(cat(3,std_cat{:,2})+cat(3,std_cat{:,1}));
vfs_softnorm = 0.2;
vfs_change = (abs(cat(3,vfs_cat{:,2}))-abs(cat(3,vfs_cat{:,1})))./(vfs_softnorm+abs(cat(3,vfs_cat{:,2}))+abs(cat(3,vfs_cat{:,1})));

figure; 

% Plot std
subplot(2,3,1);
imagesc(nanmean(cat(3,std_cat{:,1}),3));
axis image off;
caxis([0,0.03]);
colormap(gca,brewermap([],'*Greys'));
AP_reference_outline('ccf_aligned','r');
title('Std (pre-muscimol');
colorbar

subplot(2,3,2);
imagesc(nanmean(cat(3,std_cat{:,2}),3));
axis image off;
caxis([0,0.03]);
colormap(gca,brewermap([],'*Greys'));
AP_reference_outline('ccf_aligned','r');
title('Std (post-muscimol');
colorbar

subplot(2,3,3);
imagesc(nanmean(std_change,3));
axis image off;
caxis([-0.5,0.5]);
colormap(gca,brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('(post-pre)/(post+pre)');
colorbar

% Plot VFS
subplot(2,3,4);
imagesc(nanmean(cat(3,vfs_cat{:,1}),3));
axis image off;
caxis([-1,1]);
colormap(gca,brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('VFS pre-muscimol');
colorbar

subplot(2,3,5);
imagesc(nanmean(cat(3,vfs_cat{:,2}),3));
axis image off;
caxis([-1,1]);
colormap(gca,brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('VFS post-muscimol');
colorbar

subplot(2,3,6);
imagesc(nanmean(vfs_change,3));
axis image off;
caxis([-0.5,0.5]);
colormap(gca,brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('(post-pre)/(post+pre)');
colorbar


%% EDFig9c Striatal spike rate pre/post muscimol

animals = {'AP045','AP054','AP055','AP053','AP047','AP048'};

spike_rate_cond = cell(size(animals));

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworldNoRepeats';
    experiments = AP_find_experiments(animal,protocol);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
        
    disp(animal);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        
        % Load data (first experiment - but spikes throughout used)
        experiment = experiments(curr_day).experiment(1);
        load_parts.ephys = true;
        AP_load_experiment
        
        % Set experiments in conditions 
        % (1-3 pre-muscimol, 4+ post-muscimol)
        % (assumes all repeated expts/failures were post-muscimol)
        curr_experiments = AP_list_experiments(animal,day);
        cond_expts = {[1:3],[4:length(curr_experiments)]};
        
        for curr_cond = 1:2
                exp_starts = sync(2).timestamps(sync(2).values == 1);
                exp_stops = sync(2).timestamps(sync(2).values == 0);
                
                curr_exp_start = exp_starts(cond_expts{curr_cond}(1));
                curr_exp_stop = exp_stops(cond_expts{curr_cond}(end));
                
                curr_use_spikes = spike_times >= curr_exp_start & ...
                    spike_times <= curr_exp_stop & ~isnan(aligned_str_depth_group);
                
                spike_rate_cond{curr_animal}{curr_day}(:,curr_cond) = ...
                    accumarray(aligned_str_depth_group(curr_use_spikes),1, ...
                    [n_aligned_depths,1])./(curr_exp_stop - curr_exp_start);
        end
         
        clearvars -except animals curr_animal animal experiments curr_day ...
            spike_rate_cond
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
    
end

% Concatenate across recordings and get change
spike_rate_cond_cat = cell2mat(reshape([spike_rate_cond{:}],1,1,[]));
spike_rate_cond_cat_change = ...
    squeeze((spike_rate_cond_cat(:,1,:) - spike_rate_cond_cat(:,2,:))./ ...
    (spike_rate_cond_cat(:,1,:)));
n_depths = size(spike_rate_cond_cat,1);

figure;
plotSpread(spike_rate_cond_cat_change','distributionColors',[0.5,0.5,0.5]);
errorbar((1:n_depths)+0.3,nanmean(spike_rate_cond_cat_change,2), ...
    AP_sem(spike_rate_cond_cat_change,2),'.k', ...
    'MarkerSize',20,'linewidth',2,'linestyle','none');
line(xlim,[0,0],'color','r')
xlabel('(post-pre)/(pre)');
ylabel('Striatal depth');
title('Muscimol change');

% (Condition statistics)
disp('Spike rate pre/post muscimol signrank:')
for curr_depth = 1:n_depths
    curr_p = signrank(squeeze(spike_rate_cond_cat(curr_depth,1,:)), ...
        squeeze(spike_rate_cond_cat(curr_depth,2,:)));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(curr_p)]);
end


%% EDFig9d Cortex/striatum passive stim pre/post muscimol

data_fns = { ...
    'trial_activity_AP_lcrGratingPassive_pre_muscimol', ...
    'trial_activity_AP_lcrGratingPassive_post_muscimol'};

stimIDs = cell(2,1);
mua_muscimol = cell(2,1);
mua_ctxpred_muscimol = cell(2,1);
fluor_muscimol = cell(2,1);
fluor_roi_muscimol = cell(2,1);
fluor_kernelroi_muscimol = cell(2,1);

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;

    % Split by experiment
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    
    % Get trials with movement during stim to exclude
    wheel_thresh = 0.025;
    quiescent_trials = ~any(abs(wheel_allcat(:,t >= -0.5 & t <= 0.5)) > wheel_thresh,2);
    quiescent_trials_exp = mat2cell(quiescent_trials,trials_recording,1);
    
    % Get stim and activity by experiment
    trial_stim_allcat_exp = mat2cell(trial_stim_allcat,trials_recording,1);  
    fluor_deconv_exp = mat2cell(fluor_allcat_deconv,trials_recording,length(t),n_vs);
    fluor_roi_deconv_exp = mat2cell(fluor_roi_deconv,trials_recording,length(t),numel(wf_roi));
    fluor_kernelroi_deconv_exp = mat2cell(fluor_kernelroi_deconv,trials_recording,length(t),n_depths);
    mua_allcat_exp = mat2cell(mua_allcat,trials_recording,length(t),n_depths);
    mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,trials_recording,length(t),n_depths);
    
    % Exclude trials with fluorescence spikes
    % (this is a dirty way to do this but don't have a better alt)
    fluor_spike_thresh = 0.015; % deconv df/f threshold (eyeballed)
    fluor_spike_trial = cellfun(@(x) any(any(x > fluor_spike_thresh,2),3), ...
        fluor_kernelroi_deconv_exp,'uni',false);
    
    % Grab data
    stimIDs{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        trial_stim_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    mua_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        mua_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    mua_ctxpred_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        mua_ctxpred_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_roi_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_roi_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_kernelroi_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_kernelroi_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
end

% Plot average fluorescence
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');

use_stim = 1;

fluor_premuscimol_mean = ...
    svdFrameReconstruct(U_master(:,:,1:200), ...
    permute(nanmean(cell2mat(cellfun(@(x,stim) ...
    nanmean(x(stim == use_stim,:,:),1),fluor_muscimol{1},stimIDs{1},'uni',false)),1),[3,2,1]));
fluor_postmuscimol_mean = ...
    svdFrameReconstruct(U_master(:,:,1:200), ...
    permute(nanmean(cell2mat(cellfun(@(x,stim) ...
    nanmean(x(stim == use_stim,:,:),1),fluor_muscimol{2},stimIDs{2},'uni',false)),1),[3,2,1]));

AP_image_scroll([fluor_premuscimol_mean,fluor_postmuscimol_mean]);
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5],[],[size(U_master,1),size(U_master,2),1,2]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
axis image;
colormap(brewermap([],'PRGn'));


figure;
t_stim = t >= 0 & t <= 0.2;

subplot(1,2,1)
imagesc(nanmean(fluor_premuscimol_mean(:,:,t_stim),3));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
c = caxis;
axis image off;
colormap(brewermap([],'PRGn'));
title('Pre-muscimol');

subplot(1,2,2)
imagesc(nanmean(fluor_postmuscimol_mean(:,:,t_stim),3));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
caxis(c);
axis image off;
colormap(brewermap([],'PRGn'));
title('Post-muscimol');

% Get pre/post stim response
use_stim = 1;

mua_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_muscimol{1},stimIDs{1},'uni',false));
mua_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_muscimol{2},stimIDs{2},'uni',false));

mua_ctxpred_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_ctxpred_muscimol{1},stimIDs{1},'uni',false));
mua_ctxpred_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_ctxpred_muscimol{2},stimIDs{2},'uni',false));

fluor_roi_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_roi_muscimol{1},stimIDs{1},'uni',false));
fluor_roi_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_roi_muscimol{2},stimIDs{2},'uni',false));

fluor_kernelroi_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_kernelroi_muscimol{1},stimIDs{1},'uni',false));
fluor_kernelroi_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_kernelroi_muscimol{2},stimIDs{2},'uni',false));

%%% MUSCIMOL CHANGE OPTION 1
% Fit scaling factor (change) pre-post
n_exps = length(stimIDs{1});
str_muscimol_change = nan(n_exps,n_depths);
ctx_muscimol_change = nan(n_exps,n_depths);
for curr_depth = 1:n_depths
    for curr_exp = 1:n_exps
        str_muscimol_change(curr_exp,curr_depth) = ...
            mua_premuscimol_mean(curr_exp,:,curr_depth)'\ ...
            mua_postmuscimol_mean(curr_exp,:,curr_depth)';
        
        ctx_muscimol_change(curr_exp,curr_depth) = ...
            fluor_kernelroi_premuscimol_mean(curr_exp,:,curr_depth)'\ ...
            fluor_kernelroi_postmuscimol_mean(curr_exp,:,curr_depth)';
    end
end

% %%% MUSCIMOL CHANGE OPTION 2
% % Get difference in average stimulus response
% t_stim = t >= 0 & t <= 0.2;
% mua_avg_premuscimol = permute(nanmean(mua_premuscimol_mean(:,t_stim,:),2),[1,3,2]);
% mua_avg_postmuscimol = permute(nanmean(mua_postmuscimol_mean(:,t_stim,:),2),[1,3,2]);
% str_muscimol_change = (mua_avg_postmuscimol-mua_avg_premuscimol);%./(mua_avg_premuscimol);
% 
% fluor_avg_premuscimol = permute(nanmean(fluor_kernelroi_premuscimol_mean(:,t_stim,:),2),[1,3,2]);
% fluor_avg_postmuscimol = permute(nanmean(fluor_kernelroi_postmuscimol_mean(:,t_stim,:),2),[1,3,2]);
% ctx_muscimol_change = (fluor_avg_postmuscimol-fluor_avg_premuscimol);%./(fluor_avg_premuscimol);

% Plot time courses and change
figure;
p = gobjects(n_depths,3);
for plot_str = 1:n_depths
    
    p(plot_str,1) = subplot(n_depths,3,(plot_str-1)*n_depths+1); hold on;
    AP_errorfill(t,nanmean(mua_premuscimol_mean(:,:,plot_str),1)', ...
        AP_sem(mua_premuscimol_mean(:,:,plot_str),1)','k');
    AP_errorfill(t,nanmean(mua_postmuscimol_mean(:,:,plot_str),1)', ...
        AP_sem(mua_postmuscimol_mean(:,:,plot_str),1)','r');
    xlim([-0.2,1])
    xlabel('Time from stim (s)')
    ylabel(['Str ' num2str(plot_str)]);
    axis square
    
    p(plot_str,2) = subplot(n_depths,3,(plot_str-1)*n_depths+2); hold on;
    AP_errorfill(t,nanmean(fluor_kernelroi_premuscimol_mean(:,:,plot_str),1)', ...
        AP_sem(fluor_kernelroi_premuscimol_mean(:,:,plot_str),1)','k');
    AP_errorfill(t,nanmean(fluor_kernelroi_postmuscimol_mean(:,:,plot_str),1)', ...
        AP_sem(fluor_kernelroi_postmuscimol_mean(:,:,plot_str),1)','r');
    xlim([-0.2,1])
    xlabel('Time from stim (s)')
    ylabel('Cortex ROI');
    axis square
    
    p(plot_str,3) = subplot(n_depths,3,(plot_str-1)*n_depths+3);
    plot(ctx_muscimol_change(:,plot_str),str_muscimol_change(:,plot_str),'.k','MarkerSize',20)
    xlabel(['Cortex ROI (post-pre)']);
    ylabel(['Str ' num2str(plot_str) ' (post-pre)']);
    line(xlim,xlim,'color','k','linestyle','--');

end
linkaxes(p(:,1),'xy');
linkaxes(p(:,2),'xy');

% (str/ctx muscimol statistics)
disp('Striatum/cortex muscimol change correlation:')
for curr_depth = 1:n_depths
    [r,p] = corr(str_muscimol_change(:,curr_depth), ...
        ctx_muscimol_change(:,curr_depth), ...
        'rows','complete','type','pearson');
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(p) ' r = ' num2str(r)]);
end


disp('Striatum/cortex muscimol signrank:')
for curr_depth = 1:n_depths
    p = signrank(mua_avg_premuscimol(:,curr_depth),mua_avg_postmuscimol(:,curr_depth));
    disp(['Str ' num2str(curr_depth) ' p = ' num2str(p)]);
    p = signrank(fluor_avg_premuscimol(:,curr_depth),fluor_avg_postmuscimol(:,curr_depth));
    disp(['Ctx ' num2str(curr_depth) ' p = ' num2str(p)]);
end


%% EDFig9e Task performance pre/post muscimol

data_fns = { ...
    'trial_activity_vanillaChoiceworldNoRepeats_pre_muscimol', ...
    'trial_activity_vanillaChoiceworldNoRepeats_post_muscimol'};


performance_mean = cell(2,1);
rxn_time_median = cell(2,1);

frac_orient_right = cell(2,1);
rxn_time = cell(2,1);
move_t_hist = cell(2,1);

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Split data by recording
    trials_allcat = size(wheel_allcat,1);
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    use_split = trials_recording;
    split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
        [1:length(use_split)]',reshape(use_split,[],1),'uni',false));
    
    % Get all fraction correct and reaction times
    performance_mean{curr_data} = cellfun(@(x) nansum(x==1)./nansum(x~=0),mat2cell(trial_outcome_allcat,use_split,1));
    rxn_time_median{curr_data} = cellfun(@nanmedian,mat2cell(move_t,use_split,1));
    
    % Get psychometric
    stim_conditions = unique(trial_stim_allcat);
    [~,stim_idx] = ismember(trial_stim_allcat,stim_conditions,'rows');
    
    trial_stim_idx_allcat_exp = mat2cell(stim_idx,use_split,1);
    trial_choice_allcat_exp = mat2cell(trial_choice_allcat,use_split,1);
    
    frac_orient_right{curr_data} = cell2mat(cellfun(@(stim,choice) ...
        accumarray(stim,choice == -1,[length(stim_conditions),1],@nanmean,NaN), ...
        trial_stim_idx_allcat_exp,trial_choice_allcat_exp,'uni',false)');
    
    % Reaction time by stim
    move_t_exp = mat2cell(move_t,use_split,1);
    rxn_time{curr_data} = cell2mat(cellfun(@(stim,rxn) ...
        accumarray(stim,rxn,[length(stim_conditions),1],@nanmedian,NaN), ...
        trial_stim_idx_allcat_exp,move_t_exp,'uni',false)');
    
    % Get histogram of reaction times
    move_t_bins = -0.2:1/sample_rate:1;
    move_t_bin_centers = move_t_bins(1:end-1) + diff(move_t_bins)./2;
    move_t_bin = mat2cell(discretize(move_t,move_t_bins),trials_recording,1);
    
    move_t_hist{curr_data} = cell2mat(cellfun(@(move_t_bin) ...
        accumarray(move_t_bin(~isnan(move_t_bin)), ...
        1/sum(~isnan(move_t_bin)),[length(move_t_bins)-1,1],@nansum,0), ...
        move_t_bin','uni',false));

end

figure; 

subplot(3,2,1);
AP_errorfill(stim_conditions, ...
    [nanmean(frac_orient_right{1},2),nanmean(frac_orient_right{2},2)], ...
    [AP_sem(frac_orient_right{1},2),AP_sem(frac_orient_right{2},2)],[0,0,0;1,0,0]);
line(xlim,[0.5,0.5],'color','k','linestyle','--');
line([0,0],ylim,'color','k','linestyle','--');
xlabel('Contrast*Side');
ylabel('Fraction orient right');

subplot(3,2,2);
AP_errorfill(stim_conditions, ...
    nanmean(frac_orient_right{2}-frac_orient_right{1},2), ...
    AP_sem(frac_orient_right{2}-frac_orient_right{1},2),[0,0,0;1,0,0]);
line(xlim,[0,0],'color','k','linestyle','--');
line([0,0],ylim,'color','k','linestyle','--');
xlabel('Contrast*Side');
ylabel('\DeltaFraction orient right');

subplot(3,2,3);
AP_errorfill(stim_conditions, ...
    [nanmean(rxn_time{1},2),nanmean(rxn_time{2},2)], ...
    [AP_sem(rxn_time{1},2),AP_sem(rxn_time{2},2)],[0,0,0;1,0,0]);
line(xlim,[0.5,0.5],'color','k','linestyle','--');
line([0,0],ylim,'color','k','linestyle','--');
xlabel('Contrast*Side');
ylabel('Reaction time');

subplot(3,2,4);
AP_errorfill(stim_conditions, ...
    nanmean(rxn_time{2}-rxn_time{1},2), ...
    AP_sem(rxn_time{2}-rxn_time{1},2),[0,0,0;1,0,0]);
line(xlim,[0,0],'color','k','linestyle','--');
line([0,0],ylim,'color','k','linestyle','--');
xlabel('Contrast*Side');
ylabel('\DeltaReaction time');

subplot(3,2,5);
AP_errorfill(move_t_bin_centers, ...
    [nanmean(move_t_hist{1},2),nanmean(move_t_hist{2},2)], ...
    [AP_sem(move_t_hist{1},2),AP_sem(move_t_hist{2},2)],[0,0,0;1,0,0]);
line([0.5,0.5],ylim,'color','k','linestyle','--');
xlabel('Reaction time');
ylabel('Fraction');

subplot(3,2,6);
AP_errorfill(move_t_bin_centers, ...
    nanmean(move_t_hist{2}-move_t_hist{1},2), ...
    AP_sem(move_t_hist{2}-move_t_hist{1},2),[0,0,0;1,0,0]);
line(xlim,[0,0],'color','k','linestyle','--');
line([0.5,0.5],ylim,'color','k','linestyle','--');
xlabel('Reaction time');
ylabel('\DeltaFraction');

% Values and statistics
disp('Performance mean:')
disp(num2str(cellfun(@nanmean,performance_mean)'));
disp('Performance sem:')
disp(num2str(cellfun(@(x) AP_sem(x,1),performance_mean)'));
disp('Reaction time mean:')
disp(num2str(cellfun(@nanmean,rxn_time_median)'));
disp('Reaction time sem:')
disp(num2str(cellfun(@(x) AP_sem(x,1),rxn_time_median)'));

disp('Stim/condition 2-way anova:');

% (Psychometric stim x condition)
curr_stat_data = permute(cat(3,frac_orient_right{:}),[1,3,2]);
[stim_grp,condition_grp,exp_grp] = ndgrid(1:size(curr_stat_data,1), ...
    1:size(curr_stat_data,2),1:size(curr_stat_data,3));
[curr_p,~,~,terms] = anovan(reshape(curr_stat_data,[],1), ...
        [stim_grp(:),condition_grp(:)],'model','interaction','display','off');
 disp(['Psychomatric condition p = ' num2str(curr_p(2))]);
disp(['Psychomatric stim x condition p = ' num2str(curr_p(3))]);

% (Reaction time stim x condition)
curr_stat_data = permute(cat(3,rxn_time{:}),[1,3,2]);
[stim_grp,condition_grp,exp_grp] = ndgrid(1:size(curr_stat_data,1), ...
    1:size(curr_stat_data,2),1:size(curr_stat_data,3));
[curr_p,~,~,terms] = anovan(reshape(curr_stat_data,[],1), ...
        [stim_grp(:),condition_grp(:)],'model','interaction','display','off');
disp(['Reaction time condition p = ' num2str(curr_p(2))]);
disp(['Reaction time stim x condition p = ' num2str(curr_p(3))]);


%% EDFig9f Striatal task kernels pre/post muscimol

data_fns = { ...
    'trial_activity_vanillaChoiceworldNoRepeats_pre_muscimol', ...
    'trial_activity_vanillaChoiceworldNoRepeats_post_muscimol'};

mua_norm_exp = cell(2,1);
task_str_kernel = cell(2,1);
task_str_ctxpred_kernel = cell(2,1);
task_ctx_roi_kernel = cell(2,1);

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Keep task > str kernels, task > ctx-pred str norm factor
    mua_norm_exp{curr_data} = vertcat(mua_norm{:});
    task_str_kernel{curr_data} = vertcat(mua_taskpred_k_all{:});
    task_str_ctxpred_kernel{curr_data} = vertcat(mua_ctxpred_taskpred_k_all{:});
    
end

% Normalize and concatenate task kernels
mua_task_k_premuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_kernel{1},mua_norm_exp{1},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

mua_ctxpred_task_k_premuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_ctxpred_kernel{1},mua_norm_exp{1},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

mua_task_k_postmuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_kernel{2},mua_norm_exp{2},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

mua_ctxpred_task_k_postmuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_ctxpred_kernel{2},mua_norm_exp{2},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

% Get task>striatum parameters
n_regressors = length(task_regressor_labels);

% Plot task>striatum kernels
stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = [];
move_col = [0.6,0,0.6;0,0.6,0];
go_col = [0.8,0.8,0.2;0.5,0.5,0.5];
outcome_col = [0,0,0.8;0.8,0,0];
task_regressor_cols = {stim_col,move_col,go_col,outcome_col};
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

figure('Name','Pre-muscimol');
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);
        
        curr_kernels = mua_task_k_premuscimol{curr_regressor}(:,:,curr_depth,:);
        n_subregressors = size(mua_task_k_premuscimol{curr_regressor},1);
        col = task_regressor_cols{curr_regressor};
        for curr_subregressor = 1:n_subregressors
            AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
                nanmean(curr_kernels(curr_subregressor,:,:,:),4)', ...
                AP_sem(curr_kernels(curr_subregressor,:,:,:),4)', ...
                col(curr_subregressor,:));
        end
        
        xlabel('Time (s)');
        ylabel('Weight');
        title(task_regressor_labels{curr_regressor});
        line([0,0],ylim,'color','k');
        
    end
end
linkaxes(p);
y_scale = 1;
t_scale = 0.5;
line(min(xlim) + [0,t_scale],repmat(min(ylim),2,1),'color','k','linewidth',3);
line(repmat(min(xlim),2,1),min(ylim) + [0,y_scale],'color','k','linewidth',3);

figure('Name','Post-muscimol');
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);
        
        curr_kernels = mua_task_k_postmuscimol{curr_regressor}(:,:,curr_depth,:);
        n_subregressors = size(mua_task_k_postmuscimol{curr_regressor},1);
        col = task_regressor_cols{curr_regressor};
        for curr_subregressor = 1:n_subregressors
            AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
                nanmean(curr_kernels(curr_subregressor,:,:,:),4)', ...
                AP_sem(curr_kernels(curr_subregressor,:,:,:),4)', ...
                col(curr_subregressor,:),0.5);
        end
        
        xlabel('Time (s)');
        ylabel('Weight');
        title(task_regressor_labels{curr_regressor});
        line([0,0],ylim,'color','k');
        
    end
end
linkaxes(p);
y_scale = 1;
t_scale = 0.5;
line(min(xlim) + [0,t_scale],repmat(min(ylim),2,1),'color','k','linewidth',3);
line(repmat(min(xlim),2,1),min(ylim) + [0,y_scale],'color','k','linewidth',3);

% Plot kernel sums pre/post muscimol
str_k_sum_premuscimol = cellfun(@(x) permute(sum(x,2),[1,3,4,2]),mua_task_k_premuscimol,'uni',false);
str_k_sum_postmuscimol = cellfun(@(x) permute(sum(x,2),[1,3,4,2]),mua_task_k_postmuscimol,'uni',false);

figure;
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors
        if curr_regressor == 1
            x = unique([0.06,0.125,0.25,0.5,1].*[-1;1]);
        else
            x = 1:size(str_k_sum_premuscimol{curr_regressor},1);
        end
        
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);
        hold on
        
        errorbar(x,nanmean(str_k_sum_premuscimol{curr_regressor}(:,curr_depth,:),3), ...
            AP_sem(str_k_sum_premuscimol{curr_regressor}(:,curr_depth,:),3),'k','linewidth',2);
        
        errorbar(x,nanmean(str_k_sum_postmuscimol{curr_regressor}(:,curr_depth,:),3), ...
            AP_sem(str_k_sum_postmuscimol{curr_regressor}(:,curr_depth,:),3),'r','linewidth',2);
        
        axis tight
        xlim(xlim + [-0.2,0.2]);
        
        xlabel('Condition');
        ylabel('Weight sum');
        title(task_regressor_labels{curr_regressor});
        
    end
end

linkaxes(p,'y');

% (Regressor weight sum by condition statistics)
disp('Regressor v condition 2-way anova (only disp condition effect):')
for curr_depth = 1:n_depths
    for curr_regressor = 1:length(task_regressor_labels)
        
        curr_prepost = permute(cat(4, ...
            str_k_sum_premuscimol{curr_regressor}(:,curr_depth,:), ...
            str_k_sum_postmuscimol{curr_regressor}(:,curr_depth,:)),[1,3,4,2]);
        
        [regressor_grp,exp_grp,condition_grp] = ndgrid( ...
            1:size(curr_prepost,1),1:size(curr_prepost,2),1:size(curr_prepost,3));
        
        curr_p = anovan(reshape(curr_prepost,[],1), ...
            [regressor_grp(:),condition_grp(:)],'display','off');
        
        disp(['Str ' num2str(curr_depth) ' ' task_regressor_labels{curr_regressor} ...
            ' p = ' num2str(curr_p(2)')]);
        
    end
end


%% @@ EDFig10b: Firing rate histogram by cell type

include_celltypes = 1:4;

% Get average firing rates of all cells
unit_fr = cell2mat(cellfun(@(x) ...
    nanmean(reshape(x,[],size(x,3)),1)', ...
    vertcat(mua_all{:}),'uni',false));

% Plot histogram of firing rates
fr_bin_edges = [logspace(-2,2,30)];

figure; hold on
celltype_col = ...
    [0.9,0.4,0.6; ...
    0.4,0.5,0.6; ...
    0.5,0.3,0.1; ...
    1,0.5,0];
set(gca,'ColorOrder',celltype_col,'XScale','log');
for curr_celltype = include_celltypes
    histogram(unit_fr(good_units_allcat & celltype_allcat == ...
        curr_celltype),fr_bin_edges, ...
        'Normalization','probability','EdgeColor','none');
end
legend(celltype_labels(1:4));
xlabel('Firing rate');
ylabel('Fraction');


%% @@ EDFig10c: Counts of cells by type and domain

include_celltypes = [1:4]; % 5 is the narrow tan-like cells
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
set(gca,'ColorOrder',celltype_col,'XTick',[]);
bar(domain_type_count,'stacked');
xlabel('Domain');
ylabel('Count');
legend(celltype_labels);


%% EDFig 10d: Example rasters

example_msns = ...
    [187,1,8;
    525,2,4;
    132,1,1];
example_fsis = ...
    [412,6,1;
    256,2,3;
    264,4,2];
example_tans = ...
    [364,5,1; ...
    290,5,1; ...
    167,4,4];
example_uins = ...
    [325,4,5;
    482,3,2;
    44,1,2];

example_units = cat(3,example_msns,example_fsis,example_tans,example_uins);

% Get animals/days using the depth alignment file
alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
ephys_depth_align_fn = [alignment_path filesep 'ephys_depth_align'];
load(ephys_depth_align_fn);

% Set up subplots
figure;
n_depths = size(example_units,1);
n_celltypes = size(example_units,3);
p = reshape(tight_subplot(n_depths*2,n_celltypes),[n_celltypes,n_depths*2])';
drawnow;

% Set raster time bins
raster_window = [-0.2,0.9];
psth_bin_size = 0.001;
raster_t_bins = raster_window(1):psth_bin_size:raster_window(2);
raster_t = raster_t_bins(1:end-1) + diff(raster_t_bins)./2;

for curr_depth = 1:n_depths
    for curr_celltype = 1:n_celltypes
        
        preload_vars = who;
        
        curr_animal = example_units(curr_depth,2,curr_celltype);
        curr_day = example_units(curr_depth,3,curr_celltype);
        curr_unit = example_units(curr_depth,1,curr_celltype);
        
        animal = ephys_depth_align(curr_animal).animal;
        day = ephys_depth_align(curr_animal).day{curr_day};
        experiment = 1;
        load_parts.ephys = true;
        AP_load_experiment
        
        % Get alignment times to use
        outcome_time = signals_events.responseTimes';
        plot_trials = ...
            (trial_conditions(1:n_trials,1).*trial_conditions(1:n_trials,2)) > 0 & ...
            trial_choice(1:n_trials) == -1 & ...
            trial_outcome(1:n_trials) == 1;
        switch curr_depth
            case 1
                use_align = stimOn_times(plot_trials);
            case 2
                if curr_celltype ~= 3
                    use_align = wheel_move_time(plot_trials);
                else
                    % (if TAN in str2, stim-align)
                    use_align = stimOn_times(plot_trials);
                    
                end
            case 3
                use_align = outcome_time(plot_trials);
        end        
        
        t_peri_event = use_align + raster_t_bins;
        % (handle NaNs by setting rows with NaN times to 0)
        t_peri_event(any(isnan(t_peri_event),2),:) = 0;
        
        % Bin spikes (use only spikes within time range, big speed-up)
        curr_spikes_idx = ismember(spike_templates,curr_unit);
        curr_raster_spike_times = spike_times_timeline(curr_spikes_idx);
        curr_raster_spike_times(curr_raster_spike_times < min(t_peri_event(:)) | ...
            curr_raster_spike_times > max(t_peri_event(:))) = [];
        
        if ~any(diff(reshape(t_peri_event',[],1)) < 0)
            % (if no backward time jumps, can do long bin and cut out in-between, faster)
            curr_raster_continuous = reshape([histcounts(curr_raster_spike_times, ...
                reshape(t_peri_event',[],1)),NaN],size(t_peri_event'))';
            curr_raster = curr_raster_continuous(:,1:end-1);
        else
            % (otherwise, bin trial-by-trial)
            curr_raster = cell2mat(arrayfun(@(x) ...
                histcounts(curr_raster_spike_times,t_peri_event(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false));
        end
        
        % Get smoothed PSTH
        smooth_size = 51;
        gw = gausswin(smooth_size,3)';
        smWin = gw./sum(gw);
        bin_t = mean(diff(raster_t));
        
        curr_psth = nanmean(curr_raster,1);
        curr_smoothed_psth = conv2(padarray(curr_psth, ...
            [0,floor(length(smWin)/2)],'replicate','both'), ...
            smWin,'valid')./bin_t;       
               
        % Plot PSTH and raster
        axes(p(curr_depth*2-1,curr_celltype));
        plot(raster_t,curr_smoothed_psth,'k','linewidth',1);
        axis tight off
        line([0,0],ylim);
        % (scalebar on first plot)
        if curr_depth == 1 & curr_celltype == 1
            line([-0.1,-0.1],[10,20],'color','r','linewidth',1);
            line([-0.1,0.1],[10,10],'color','r','linewidth',1);
        end
        
        axes(p(curr_depth*2,curr_celltype));
        [raster_y,raster_x] = find(curr_raster);
        plot(raster_t(raster_x),raster_y,'.k');
        axis tight off    
        % (scalebar on first plot)
        if curr_depth == 1 & curr_celltype == 1
            line([-0.1,-0.1],[10,30],'color','r','linewidth',1);
        end
        
        drawnow;

        clearvars('-except',preload_vars{:})
        
    end
end

linkaxes(p(:),'xy')
set(p,'XLim',raster_window);



























