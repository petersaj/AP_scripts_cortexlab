%% Test analysis for imaging with operant task (over learning)
% (note: this includes things brought from test_corticostriatal and
% test_learning)

%% Get average rewarded and other movements

[wheel_vel,wheel_move,wheel_vel_split] = AP_parse_wheel(wheel_position,Timeline.hw.daqSampleRate);

wheel_split_n = cellfun(@length,wheel_vel_split);
wheel_pos_split = mat2cell(wheel_position,wheel_split_n);

wheel_split_n = cellfun(@length,wheel_vel_split)';
split_move = cellfun(@(x) x(1),mat2cell(wheel_move,wheel_split_n));

max_move = max(wheel_split_n(split_move == 1));
vel_pad = cell2mat(cellfun(@(x) padarray(x,max_move-length(x),NaN,'post'), ...
    wheel_vel_split(split_move == 1),'uni',false)');
pos_pad_raw = cell2mat(cellfun(@(x) padarray(x,max_move-length(x),NaN,'post'), ...
    wheel_pos_split(split_move == 1),'uni',false)'); 
pos_pad = pos_pad_raw - pos_pad_raw(1,:);


rewarded_epochs = cellfun(@(t) any(ismember(t,reward_t_timeline)), ...
    mat2cell(Timeline.rawDAQTimestamps',wheel_split_n));

rewarded_movements = rewarded_epochs(split_move == 1);

% Plot average rewarded and non-rewarded movements
figure; hold on
plot(nanmean(vel_pad(:,~rewarded_movements),2),'k');
plot(nanmean(vel_pad(:,rewarded_movements),2),'b');
legend({'Unrewarded','Rewarded'});
xlabel('Time from movement onset');
ylabel('Wheel velocity');


%% ~~~~~~~~~ GRAB & SAVE BATCH  ~~~~~~~~~

%% Passive - corticostriatal

clear all
disp('Passive trial activity (corticostriatal)')

% (AP089/90/91: not variable quiescence)
animals = {'AP092','AP093','AP094','AP095','AP096','AP097'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    experiments_full = AP_find_experiments(animal,protocol);
    imaged_experiments = [experiments_full.imaging];

    % Use only days with operant training
    operant_protocol = 'AP_stimWheelRight';
    operant_experiments = AP_find_experiments(animal,operant_protocol,false);
    experiment_operant = ismember({experiments_full.day},{operant_experiments.day});
    
    % Don't use days with retinotopy (those were often muscimol days)
    retinotopy_protocol = 'AP_kalatsky';
    retinotopy_experiments = AP_find_experiments(animal,retinotopy_protocol,false);
    experiment_retintopy = ismember({experiments_full.day},{retinotopy_experiments.day});
    
    use_experiments = [experiments_full.imaging] & experiment_operant & ...
        ~experiment_retintopy;
    experiments = experiments_full(use_experiments);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        AP_load_experiment;
        
        % Pull out trial data
        AP_ctx_str_grab_trial_data;
        
        % Store trial data into master structure
        trial_data_fieldnames = fieldnames(trial_data);
        for curr_trial_data_field = trial_data_fieldnames'
            trial_data_all.(cell2mat(curr_trial_data_field)){curr_animal,1}{curr_day,1} = ...
                trial_data.(cell2mat(curr_trial_data_field));
        end
        
        % Store general info
        trial_data_all.animals = animals;
        trial_data_all.t = t;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
    
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = ['trial_activity_passive_corticostriatal'];
save([save_path filesep save_fn],'-v7.3');
disp(['Saved: ' save_path filesep save_fn])

%% Task - corticostriatal

clear all
disp('Task trial activity (corticostriatal)')

% (AP089/90/91: not variable quiescence)
animals = {'AP092','AP093','AP094','AP095','AP096','AP097'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_stimWheelRight';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & ~[experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        AP_load_experiment;
        
        % Pull out trial data
        test_cstr_grab_trial_data;
        
        % Store trial data into master structure
        trial_data_fieldnames = fieldnames(trial_data);
        for curr_trial_data_field = trial_data_fieldnames'
            trial_data_all.(cell2mat(curr_trial_data_field)){curr_animal,1}{curr_day,1} = ...
                trial_data.(cell2mat(curr_trial_data_field));
        end
        
        % Store general info
        trial_data_all.animals = animals;
        trial_data_all.t = t;
        trial_data_all.task_regressor_labels = task_regressor_labels;
        trial_data_all.task_regressor_sample_shifts = task_regressor_sample_shifts;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = ['trial_activity_task_corticostriatal'];
save([save_path filesep save_fn],'-v7.3');


%% Passive - tetO-GC6s

clear all
disp('Passive trial activity (tetO)')

animals = {'AP100','AP101','AP103','AP104','AP105','AP106'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    experiments_full = AP_find_experiments(animal,protocol);
    imaged_experiments = [experiments_full.imaging];

    % Use only days with operant training
    operant_protocol = 'AP_stimWheelRight';
    operant_experiments = AP_find_experiments(animal,operant_protocol,false);
    experiment_operant = ismember({experiments_full.day},{operant_experiments.day});
    
    % Use days before muscimol experiment starts
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    muscimol_start_day = muscimol(muscimol_animal_idx).day{1};
    muscimol_experiments = datenum({experiments_full.day})' >= datenum(muscimol_start_day);
    
    % (set experiments to use)
    use_experiments = [experiments_full.imaging] & experiment_operant & ...
        ~muscimol_experiments;
    experiments = experiments_full(use_experiments);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        AP_load_experiment;
        
        % Pull out trial data
        AP_ctx_str_grab_trial_data;
        
        % Store trial data into master structure
        trial_data_fieldnames = fieldnames(trial_data);
        for curr_trial_data_field = trial_data_fieldnames'
            trial_data_all.(cell2mat(curr_trial_data_field)){curr_animal,1}{curr_day,1} = ...
                trial_data.(cell2mat(curr_trial_data_field));
        end
        
        % Store general info
        trial_data_all.animals = animals;
        trial_data_all.t = t;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
    
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = ['trial_activity_passive_teto'];
save([save_path filesep save_fn],'-v7.3');
disp(['Saved: ' save_path filesep save_fn])


%% Task - tetO-GC6s

clear all
disp('Task trial activity (tetO)')

animals = {'AP100','AP101','AP103','AP104','AP105','AP106'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_stimWheelRight';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & ~[experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        AP_load_experiment;
        
        % Pull out trial data
        test_cstr_grab_trial_data;
        
        % Store trial data into master structure
        trial_data_fieldnames = fieldnames(trial_data);
        for curr_trial_data_field = trial_data_fieldnames'
            trial_data_all.(cell2mat(curr_trial_data_field)){curr_animal,1}{curr_day,1} = ...
                trial_data.(cell2mat(curr_trial_data_field));
        end
        
        % Store general info
        trial_data_all.animals = animals;
        trial_data_all.t = t;
        trial_data_all.task_regressor_labels = task_regressor_labels;
        trial_data_all.task_regressor_sample_shifts = task_regressor_sample_shifts;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = ['trial_activity_task_teto'];
save([save_path filesep save_fn],'-v7.3');



%% ~~~~~~~~~ BATCH ANALYSIS ~~~~~~~~~

%% >> Passive 

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_passive_corticostriatal';
AP_load_concat_normalize_ctx_str;

% Get animal and day index for each trial
trial_animal = cell2mat(arrayfun(@(x) ...
    x*ones(size(vertcat(wheel_all{x}{:}),1),1), ...
    [1:length(wheel_all)]','uni',false));
trial_day = cell2mat(cellfun(@(x) cell2mat(cellfun(@(curr_day,x) ...
    curr_day*ones(size(x,1),1),num2cell(1:length(x))',x,'uni',false)), ...
    wheel_all,'uni',false));

% Get trials with movement during stim to exclude
wheel_thresh = 0.025;
quiescent_trials = ~any(abs(wheel_allcat(:,t >= -0.5 & t <= 0.5)) > wheel_thresh,2);

% Get average fluorescence by animal, day, stim
stim_unique = unique(trial_stim_allcat);
stim_v_avg = cell(size(animals));
for curr_animal = 1:length(animals)        
    for curr_day = 1:max(trial_day(trial_animal == curr_animal))
        for curr_stim_idx = 1:length(stim_unique)
            use_trials = quiescent_trials & ...
                trial_animal == curr_animal & ...
                trial_day == curr_day & ...
                trial_stim_allcat == stim_unique(curr_stim_idx);
            stim_v_avg{curr_animal}(:,:,curr_day,curr_stim_idx) = ...
                permute(nanmean(fluor_allcat_deconv(use_trials,:,:),1),[3,2,1]);
        end       
    end
end


% (average stim response for each day)
use_stim = 3;
min_days = min(cellfun(@(x) size(x,3),stim_v_avg));
stim_v_avg_dayavg = nanmean(cell2mat(permute(cellfun(@(x) x(:,:,1:min_days,use_stim), ...
    stim_v_avg,'uni',false),[1,3,4,2])),4);
stim_px_avg_dayavg = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    stim_v_avg_dayavg);
AP_image_scroll(stim_px_avg_dayavg,t);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
% (as above but within each mouse)
use_t = t > 0.05 & t < 0.2;
stim_px_avg_day = AP_svdFrameReconstruct( ...
    U_master(:,:,1:n_vs),cell2mat(permute(cellfun(@(x) ...
    squeeze(nanmean(x(:,use_t,1:min_days,use_stim),2)), ...
    stim_v_avg,'uni',false),[1,3,2])));
AP_image_scroll(stim_px_avg_day);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

% (average stim response "pre/post learning")
naive_days = 1:3;
% expert_days = 15:17;
expert_days = 6:7;

stim_avg_naive = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    nanmean(cell2mat(permute(cellfun(@(x) ...
    squeeze(nanmean(x(:,:,naive_days,:),3)), ...
    stim_v_avg,'uni',false),[1,3,4,2])),4));
stim_avg_expert = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    nanmean(cell2mat(permute(cellfun(@(x) ...
    squeeze(nanmean(x(:,:,expert_days,:),3)), ...
    stim_v_avg,'uni',false),[1,3,4,2])),4));

use_t = t > 0.05 & t < 0.15;
use_stim = 3;
stim_avg_learning = cat(3,nanmean(stim_avg_naive(:,:,use_t,use_stim),3), ...
    nanmean(stim_avg_expert(:,:,use_t,use_stim),3));

figure;
subplot(1,3,1);
imagesc(stim_avg_learning(:,:,1));
axis image off;
AP_reference_outline('ccf_aligned','k');
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
title('Naive');

subplot(1,3,2);
imagesc(stim_avg_learning(:,:,2));
axis image off;
AP_reference_outline('ccf_aligned','k');
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
title('Expert');

subplot(1,3,3);
imagesc(diff(stim_avg_learning,[],3));
axis image off;
AP_reference_outline('ccf_aligned','k');
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
title('Difference');

%% >> Task

% Task
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_task_teto';

AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_recording;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));



%% Task kernels

% Get task>cortex parameters
n_regressors = length(task_regressor_labels);
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

% Get average task > cortex kernels (V's and ROIs)
regressor_v = cell(n_regressors,1);
for curr_regressor = 1:n_regressors
    regressor_v{curr_regressor} = nanmean(cell2mat(cellfun(@(x) x{curr_regressor}, ...
        permute(vertcat(fluor_taskpred_k_all{:}),[2,3,4,1]),'uni',false)),4);   
    AP_print_progress_fraction(curr_regressor,n_regressors);
end

% Get regressor pixels
regressor_px = cellfun(@(v) cell2mat(arrayfun(@(subregressor) ...
    svdFrameReconstruct(U_master(:,:,1:n_vs),permute(v(subregressor,:,:),[3,2,1])), ...
    permute(1:size(v,1),[1,3,4,2]),'uni',false)),regressor_v,'uni',false);

plot_regressor = 4;
AP_image_scroll(regressor_px{plot_regressor})
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
axis image;












