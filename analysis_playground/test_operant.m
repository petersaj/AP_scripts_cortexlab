%% Test analysis for imaging with operant task (over learning)
% (note: this includes things brought from test_corticostriatal and
% test_learning)

%% ~~~~~~~~~ SINGLE-SESSION ~~~~~~~~~

%% Passive PSTH

% Get quiescent trials
wheel_window = [0,0.5];
wheel_window_t = wheel_window(1):1/Timeline.hw.daqSampleRate:wheel_window(2);
wheel_window_t_peri_event = bsxfun(@plus,stimOn_times,wheel_window_t);
event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
    wheel_velocity,wheel_window_t_peri_event);
quiescent_trials = ~any(abs(event_aligned_wheel) > 0,2);

% PSTH viewer
AP_cellraster(stimOn_times(quiescent_trials),stimIDs(quiescent_trials))


%% Operant PSTH

% Standard events
AP_cellraster({stimOn_times,wheel_move_time,reward_t_timeline})


% Align rewarded and non-rewarded movements
[wheel_vel,wheel_move,wheel_vel_split] = AP_parse_wheel(wheel_position,Timeline.hw.daqSampleRate);
rewarded_epochs = cellfun(@(t) any(ismember(t,reward_t_timeline)), ...
    mat2cell(Timeline.rawDAQTimestamps',cellfun(@length,wheel_vel_split)));
rewarded_movements = wheel_vel_split(rewarded_epochs == 1);




%% Get average rewarded and other movements

[wheel_vel,wheel_move,wheel_vel_split] = AP_parse_wheel(wheel_position,Timeline.hw.daqSampleRate);

wheel_split_n = cellfun(@length,wheel_vel_split);
wheel_pos_split = mat2cell(wheel_position,wheel_split_n);

wheel_split_n = cellfun(@length,wheel_vel_split)';
split_move = cellfun(@(x) x(1),mat2cell(wheel_move,wheel_split_n));

wheel_gain = 8; % (deg/mm - in signals protocol)

max_move = max(wheel_split_n(split_move == 1));
vel_pad = cell2mat(cellfun(@(x) padarray(x,max_move-length(x),NaN,'post'), ...
    wheel_vel_split(split_move == 1),'uni',false)');
pos_pad_raw = cell2mat(cellfun(@(x) padarray(x,max_move-length(x),NaN,'post'), ...
    wheel_pos_split(split_move == 1),'uni',false)'); 
pos_pad = (pos_pad_raw - pos_pad_raw(1,:))*wheel_gain;

rewarded_epochs = cellfun(@(t) any(ismember(t,reward_t_timeline)), ...
    mat2cell(Timeline.rawDAQTimestamps',wheel_split_n));

rewarded_movements = rewarded_epochs(split_move == 1);

% Plot average rewarded and non-rewarded movements
figure; 

subplot(2,1,1);hold on
plot(nanmean(pos_pad(:,~rewarded_movements),2),'k');
plot(nanmean(pos_pad(:,rewarded_movements),2),'b');
legend({'Unrewarded','Rewarded'});
xlabel('Time from movement onset');
ylabel('Wheel position');
xlim([0,1000]);

subplot(2,1,2); hold on
plot(nanmean(vel_pad(:,~rewarded_movements),2),'k');
plot(nanmean(vel_pad(:,rewarded_movements),2),'b');
legend({'Unrewarded','Rewarded'});
xlabel('Time from movement onset');
ylabel('Wheel velocity');
xlim([0,1000]);


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
    experiments = AP_find_experiments(animal,protocol);
    
    % Get days with muscimol
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    muscimol_start_day = muscimol(muscimol_animal_idx).day{1};
    muscimol_experiments = datenum({experiments.day})' >= datenum(muscimol_start_day);
    
    % Set experiments to use (imaging, not muscimol)
    experiments = experiments([experiments.imaging] & ~muscimol_experiments);
    
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
    
    % Get days with muscimol
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    muscimol_start_day = muscimol(muscimol_animal_idx).day{1};
    muscimol_experiments = datenum({experiments.day})' >= datenum(muscimol_start_day);
    
    % Set experiments to use (imaging, not muscimol)
    experiments = experiments([experiments.imaging] & ~muscimol_experiments);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        AP_load_experiment;
        
        % Pull out trial data
        operant_grab_trial_data;
        
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
    experiments = AP_find_experiments(animal,protocol);
    
    % Get days with muscimol
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    muscimol_start_day = muscimol(muscimol_animal_idx).day{1};
    muscimol_experiments = datenum({experiments.day})' >= datenum(muscimol_start_day);
    
    % Set experiments to use (imaging, not muscimol)
    experiments = experiments([experiments.imaging] & ~muscimol_experiments);
    
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
       
    % Get days with muscimol
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    muscimol_start_day = muscimol(muscimol_animal_idx).day{1};
    muscimol_experiments = datenum({experiments.day})' >= datenum(muscimol_start_day);
    
    % Set experiments to use (imaging, not muscimol)
    experiments = experiments([experiments.imaging] & ~muscimol_experiments);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        AP_load_experiment;
        
        % Pull out trial data
        operant_grab_trial_data;
        
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

%% Muscimol: Passive - tetO-GC6s

clear all
disp('Passive trial activity (tetO)')

animals = {'AP100','AP101','AP103','AP104','AP105','AP106'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);
    
    % Get days with muscimol
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    muscimol_experiments = ismember({experiments.day},muscimol(muscimol_animal_idx).day);
    
    % Set experiments to use (imaging, not muscimol)
    experiments = experiments([experiments.imaging] & muscimol_experiments);
    
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
        
        % Store muscimol info
        muscimol_day_idx = ismember(muscimol(muscimol_animal_idx).day,day);
        trial_data_all.muscimol_area{curr_animal}{curr_day} = ...
            muscimol(muscimol_animal_idx).area{muscimol_day_idx};
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
    
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = ['trial_activity_passive_teto_muscimol'];
save([save_path filesep save_fn],'-v7.3');
disp(['Saved: ' save_path filesep save_fn])


%% ~~~~~~~~~ BATCH ANALYSIS ~~~~~~~~~

%% >> Passive 

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_passive_tetO';
% data_fn = 'trial_activity_passive_corticostriatal';

AP_load_concat_normalize_ctx_str;

% Get animal and day index for each trial
trial_animal = cell2mat(arrayfun(@(x) ...
    x*ones(size(vertcat(wheel_all{x}{:}),1),1), ...
    [1:length(wheel_all)]','uni',false));
trial_day = cell2mat(cellfun(@(x) cell2mat(cellfun(@(curr_day,x) ...
    curr_day*ones(size(x,1),1),num2cell(1:length(x))',x,'uni',false)), ...
    wheel_all,'uni',false));

% Get trials with movement during stim to exclude
quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > 0,2);

% Get average fluorescence by animal, day, stim
stim_unique = unique(trial_stim_allcat);
stim_v_avg = cell(size(animals));
stim_roi_avg = cell(size(animals));
stim_roi = cell(size(animals));
for curr_animal = 1:length(animals)        
    for curr_day = 1:max(trial_day(trial_animal == curr_animal))
        for curr_stim_idx = 1:length(stim_unique)
            use_trials = quiescent_trials & ...
                trial_animal == curr_animal & ...
                trial_day == curr_day & ...
                trial_stim_allcat == stim_unique(curr_stim_idx);
            stim_v_avg{curr_animal}(:,:,curr_day,curr_stim_idx) = ...
                permute(nanmean(fluor_allcat_deconv(use_trials,:,:),1),[3,2,1]);
            stim_roi_avg{curr_animal}(:,:,curr_day,curr_stim_idx) = ...
                permute(nanmean(fluor_roi_deconv(use_trials,:,:),1),[3,2,1]);
            
            stim_roi{curr_animal}{curr_day}{curr_stim_idx} = ...
                fluor_roi_deconv(use_trials,:,:);
        end       
    end
end

% (average stim response for each day)
use_stim = 3;

n_days = cellfun(@length,trial_info_all);
max_days = max(n_days);
stim_v_avg_dayavg = nanmean(cell2mat(permute(cellfun(@(x) ...
    padarray(x,[0,0,max_days-size(x,3),0],NaN,'post'), ...
    stim_v_avg,'uni',false),[1,3,4,5,2])),5);
stim_px_avg_dayavg = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    stim_v_avg_dayavg(:,:,:,use_stim));
AP_image_scroll(stim_px_avg_dayavg,t);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

% (plot time average by day)
use_t = t > 0.05 & t < 0.2;
stim_px_avg_dayavg_tavg = squeeze(nanmean(stim_px_avg_dayavg(:,:,use_t,:),3));

min_days = min(n_days);
c = repmat(prctile(stim_px_avg_dayavg_tavg(:),95),1,2).*[-1,1];
figure;
for curr_day = 1:min_days
    subplot(1,min_days,curr_day);
    imagesc(stim_px_avg_dayavg_tavg(:,:,curr_day));
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    caxis(c)
    axis image off;
    colormap(brewermap([],'PrGn'));
    if curr_day <= 3
        title('Passive')
    else
       title(['Training ' num2str(curr_day-3)]) 
    end
end


% (average stim response "pre/post learning")
naive_days = 1:3;
unlearned_days = 4:8;
expert_days = 9:13;

stim_avg_naive = nanmean(stim_px_avg_dayavg(:,:,:,naive_days),4);
stim_avg_unlearned = nanmean(stim_px_avg_dayavg(:,:,:,unlearned_days),4);
stim_avg_learned = nanmean(stim_px_avg_dayavg(:,:,:,expert_days),4);

AP_image_scroll([stim_avg_naive,stim_avg_unlearned,stim_avg_learned]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;

use_t = t > 0.05 & t < 0.2;
stim_learning_tavg = cat(3,nanmean(stim_avg_naive(:,:,use_t),3), ...
    nanmean(stim_avg_unlearned(:,:,use_t),3),nanmean(stim_avg_learned(:,:,use_t),3));
c = repmat(prctile(stim_learning_tavg(:),95),1,2).*[-1,1];

figure; 
for i = 1:3
    subplot(1,3,i)
    imagesc(stim_learning_tavg(:,:,i));
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    caxis(c)
    axis image off;
    colormap(brewermap([],'PrGn'));
    switch i
        case 1
            title('Naive');
        case 2
            title('Early training');
        case 3
            title('Late training');
    end
end


% Plot stim response in single animal over days
% (exclude days 1-3: pre-behavior);
use_animal = 2;
curr_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),stim_v_avg{use_animal}(:,:,4:end,3));
AP_image_scroll(curr_px,t);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;

%% ROI-ROI correlation

use_stim = 3;
stim_roi_avg_t = cellfun(@(x) cellfun(@(x) ...
    squeeze(nanmean(x{use_stim}(:,use_t,:),2)),x,'uni',false),stim_roi,'uni',false);

stim_roi_avg_t_corr = cellfun(@(x) cellfun(@corrcoef,x,'uni',false),stim_roi_avg_t,'uni',false);
v1_avg_corr = cellfun(@(x) cell2mat(cellfun(@(x) x(1,:),x,'uni',false)'),stim_roi_avg_t_corr,'uni',false);

figure; 

subplot(2,1,1); hold on;
for i = 1:length(v1_avg_corr);plot(v1_avg_corr{i}(:,3));end
ylabel('V1-AM corr');
ylim([-0.2,1]);

subplot(2,1,2); hold on;
for i = 1:length(v1_avg_corr);plot(v1_avg_corr{i}(:,7));end
ylabel('V1-FRm corr');
ylim([-0.2,1]);


curr_animal = 1;
k = nan(6,n_days(curr_animal));
explained_var = nan(n_days(curr_animal),1);
for curr_day = 1:n_days(curr_animal)
    [k(:,curr_day),~,curr_expl_var] = ...
        AP_regresskernel(stim_roi_avg_t{curr_animal}{curr_day}(:,[1,2,3,4,5,6])', ...
        stim_roi_avg_t{curr_animal}{curr_day}(:,7)',0,0,[false,false],5,false,true);
    explained_var(curr_day) = curr_expl_var.total;
end




%% Pixel-behavior correlation 
% (run behavior first)

% Get average stim (skip first 3 days: no behavior)
stim_t = [0.05,0.2];
stim_px = cellfun(@(x) ...
    squeeze(AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    nanmean(x(:,t >= stim_t(1) & t <= stim_t(2),4:end,3),2))), ...
    stim_v_avg,'uni',false);

stim_roi_stimavg = cellfun(@(x) ...
    squeeze(nanmean(x(:,t >= stim_t(1) & t <= stim_t(2),4:end,3),2)), ...
    stim_roi_avg,'uni',false);

% Correlate behavior measure with pixels

stim_move_t_median = cellfun(@(x) cellfun(@nanmedian,x),{bhv.stim_move_t},'uni',false);
stim_reward_t_median = cellfun(@(x) cellfun(@nanmedian,x),{bhv.stim_reward_t},'uni',false);
stim_move_reward_t_median = cellfun(@(x,y) cellfun(@(x,y) nanmedian(y-x),x,y),{bhv.stim_move_t},{bhv.stim_reward_t},'uni',false);



stim_surround_t = bhv(curr_animal).stim_surround_t;
poststim_move_frac_max = cellfun(@(x) cellfun(@(x) ...
    max(nanmean(x(:,stim_surround_t > 0),1)),x), ...
    {bhv.stim_surround_wheel},'uni',false);

use_bhv = poststim_move_frac_max;

[r_long,p_long] = cellfun(@(px,bhv) corr(reshape(px,[],size(px,3))', ...
    bhv(1:size(px,3))'),stim_px,use_bhv,'uni',false);

r = cell2mat(permute(cellfun(@(x) reshape(x,size(U_master(:,:,1))),r_long,'uni',false),[1,3,2]));
p = cell2mat(permute(cellfun(@(x) reshape(x,size(U_master(:,:,1))),p_long,'uni',false),[1,3,2]));

r_mean = nanmean(r,3);
p_mean = nanmean(p,3);

figure;
imagesc(r_mean,'AlphaData',1-p_mean);
axis image off
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
title('Kernel:Dprime corr')
AP_reference_outline('ccf_aligned','k');


% For more data: correlate all trials instead of avg
use_bhv = stim_move_t_median';

trial_bhv = cell(size(animals));
for curr_animal = 1:length(animals)
    curr_bhv = [nan(1,3),use_bhv{curr_animal}(1:n_days(curr_animal)-3)];
    trial_bhv{curr_animal} = cellfun(@(x,y) repmat(x,y,1), num2cell(curr_bhv)', ...
        cellfun(@(x) size(x,1),wheel_all{curr_animal},'uni',false),'uni',false);   
end
trial_bhv_allcat = cell2mat(vertcat(trial_bhv{:}));

use_trials = quiescent_trials & trial_stim_allcat == 1;

U_downsample_factor = 10;
Ud = imresize(U_master,1/U_downsample_factor,'bilinear');
trial_px = AP_svdFrameReconstruct(Ud(:,:,1:n_vs), ...
    permute(fluor_allcat_deconv(use_trials,:,:),[3,2,1]));

corr_px = nan(size(Ud,1),size(Ud,2),length(t));
corr_px_p = nan(size(Ud,1),size(Ud,2),length(t));
for curr_t = 1:length(t)
    [r,p] = corr(reshape(trial_px(:,:,curr_t,:),[],sum(use_trials))', ...
        trial_bhv_allcat(use_trials),'rows','complete');
    corr_px(:,:,curr_t) = reshape(r,size(Ud,1),size(Ud,2));
end
AP_image_scroll(corr_px,t);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
axis image;


% (testing: xcorr roi and move_t)
use_roi = 1;
a = cell2mat(cellfun(@(x,y) xcorr(x(use_roi,:),y(1:size(x,2)),7,'coeff'), ...
    stim_roi_stimavg,poststim_move_frac_max,'uni',false)')';
figure; hold on
plot(a);
plot(nanmean(a,2),'k','linewidth',2);
xlabel('Lag (move_t relative to roi activity)');
ylabel('xcorr');


%% >> Task

% Task
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_task_teto';
% data_fn = 'trial_activity_task_corticostriatal';

AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_animal;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));


%% Average trial activity

% (package back into animal/day)
n_days = cellfun(@length,wheel_all);
fluor_deconv_animal = mat2cell(mat2cell(fluor_allcat_deconv,trials_recording,length(t),n_vs),n_days);
move_t_animal = mat2cell(mat2cell(move_t,trials_recording),n_days);

% Average activity within animal/day
max_days = max(n_days);
fluor_deconv_cat = nan(n_vs,length(t),max_days,length(animals));
for curr_animal = 1:length(animals)
    fluor_deconv_cat(:,:,1:n_days(curr_animal),curr_animal) = ...
        permute(cell2mat(cellfun(@(x) nanmean(x,1), ...
        fluor_deconv_animal{curr_animal},'uni',false)),[3,2,1]);
end

% Plot average for each day
curr_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),nanmean(fluor_deconv_cat,4));
AP_image_scroll(curr_px,t);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'))

% Plot animal across days
use_animal = 6;
curr_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),fluor_deconv_cat(:,:,:,use_animal));
AP_image_scroll(curr_px,t);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'))




%% Task kernels

% Get task>cortex parameters
n_regressors = length(task_regressor_labels);
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

% Get average task > cortex kernels separately for each day
n_days = cellfun(@length,trial_info_all);
max_days = max(n_days);

% Make grid of kernel/day/mouse
fluor_taskpred_k_cat = cell(n_regressors,max_days,length(animals));
for curr_animal = 1:length(animals)
    fluor_taskpred_k_cat(:,1:n_days(curr_animal),curr_animal) = ...
        horzcat(fluor_taskpred_k_all{curr_animal}{:});
end

ctx_wheel_k_cat = cell(max_days,length(animals));
for curr_animal = 1:length(animals)
    ctx_wheel_k_cat(1:n_days(curr_animal),curr_animal) = ...
        horzcat(ctx_wheel_k_all{curr_animal});
end

% Average regressors by day
k_day_mean = cell(n_regressors,max_days);
for curr_regressor = 1:n_regressors
    for curr_day = 1:max_days
        k_day_mean{curr_regressor,curr_day} = permute(nanmean(cat(4, ...
            fluor_taskpred_k_cat{curr_regressor,curr_day,:}),4),[3,2,1]);
    end
end

ctx_wheel_k_day_mean = cell(1,max_days);
for curr_day = 1:max_days
    ctx_wheel_k_day_mean{curr_day} = nanmean(cat(4, ...
        ctx_wheel_k_cat{curr_day,:}),4);
end


% Get regressor pixels
k_day_mean_px = cellfun(@(x) AP_svdFrameReconstruct(U_master(:,:,1:n_vs),x), ...
    k_day_mean,'uni',false);

ctx_wheel_k_day_mean_px = cellfun(@(x) AP_svdFrameReconstruct(U_master(:,:,1:size(x,1)),x), ...
    ctx_wheel_k_day_mean,'uni',false);

% Plot average kernels
% (by day)
curr_regressor = 2;
curr_subregressor = 1;
curr_regressor_plot = cell2mat(permute(cellfun(@(x) x(:,:,:,curr_subregressor), ...
    k_day_mean_px(curr_regressor,:),'uni',false),[1,3,4,2]));

AP_image_scroll(curr_regressor_plot);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
axis image;

figure;imagesc(reshape(squeeze(sum(curr_regressor_plot,3)),size(U_master,1),[]));
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image off;


% (by unlearned vs learned)
unlearned_days = 1:5;
learned_days = 6:10;
curr_regressor_unlearned = nanmean(curr_regressor_plot(:,:,:,unlearned_days),4);
curr_regressor_learned = nanmean(curr_regressor_plot(:,:,:,learned_days),4);
AP_image_scroll([curr_regressor_unlearned,curr_regressor_learned,curr_regressor_learned-curr_regressor_unlearned]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;


% Plot kernel over days in single animal
use_animal = 4;
curr_k = cell2mat(cellfun(@(x) x{1},fluor_taskpred_k_all{use_animal},'uni',false));
curr_k_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),permute(curr_k,[3,2,1]));
AP_image_scroll(curr_k_px);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;



%% Within-animal: correlate reaction time to trial activity?

% (package back into animal/day)
n_days = cellfun(@length,wheel_all);
fluor_deconv_animal = mat2cell(mat2cell(fluor_allcat_deconv,trials_recording,length(t),n_vs),n_days);
move_t_animal = mat2cell(mat2cell(move_t,trials_recording),n_days);

% Plot time to movement for each animal/day
figure('Name','Time to move'); p = tight_subplot(length(animals),1);
for curr_animal = 1:length(animals)
    subplot(p(curr_animal));
    
    plot(vertcat(move_t_animal{curr_animal}{:}),'.k');
    
    curr_n_trials = cumsum(cellfun(@length,move_t_animal{curr_animal}));
    for i = 1:length(curr_n_trials)
        line([curr_n_trials,curr_n_trials],ylim,'color','r','linewidth',2);
    end
    ylabel(animals(curr_animal));
    set(gca,'XTick',[]);
end


%% >> Muscimol: passive

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_passive_teto_muscimol';

AP_load_concat_normalize_ctx_str;

% (package deconv back into animal/day)
fluor_deconv_exp = mat2cell(fluor_allcat_deconv,trials_recording,length(t),n_vs);
fluor_roi_exp = mat2cell(fluor_roi_deconv,trials_recording,length(t),n_rois);
muscimol_area_cat = horzcat(muscimol_area{:})';

quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > 0,2);
quiescent_trials_exp = mat2cell(quiescent_trials,trials_recording,1);
trial_stim_exp = mat2cell(trial_stim_allcat,trials_recording,1);

fluor_deconv_exp_stimavg = cell2mat(cellfun(@(act,quiescent,stim) ...
    nanmean(act(quiescent & stim == 1,:,:),1), ...
    fluor_deconv_exp,quiescent_trials_exp,trial_stim_exp,'uni',false));
fluor_roi_exp_stimavg = cell2mat(cellfun(@(act,quiescent,stim) ...
    nanmean(act(quiescent & stim == 1,:,:),1), ...
    fluor_roi_exp,quiescent_trials_exp,trial_stim_exp,'uni',false));

unique_muscimol_area = unique(muscimol_area_cat);
fluor_muscimol_avg = nan(n_vs,length(t),length(unique_muscimol_area));
for curr_area = 1:length(unique_muscimol_area)
    curr_exps = strcmp(muscimol_area_cat,unique_muscimol_area{curr_area});
    fluor_muscimol_avg(:,:,curr_area) = ...
        permute(nanmean(fluor_deconv_exp_stimavg(curr_exps,:,:),1),[3,2,1]);
end

fluor_muscimol_avg_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    fluor_muscimol_avg);

AP_image_scroll(fluor_muscimol_avg_px)
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;

% Plot time averages and ROIs
use_t = t > 0.05 & t < 0.2;
fluor_roi_stimavg = squeeze(nanmean(fluor_roi_exp_stimavg(:,use_t,:),2));
fluor_muscimol_px_stimavg = squeeze(nanmean(fluor_muscimol_avg_px(:,:,use_t,:),3));

figure;
for curr_area = 1:length(unique_muscimol_area)
    subplot(1,length(unique_muscimol_area)+1,curr_area);
    imagesc(fluor_muscimol_px_stimavg(:,:,curr_area));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'PrGn'));
    axis image off;
    AP_reference_outline('ccf_aligned','k');
    title(unique_muscimol_area{curr_area});
end

subplot(1,length(unique_muscimol_area)+1,length(unique_muscimol_area)+1);
plot_rois = [1,7];
area_col = lines(length(unique_muscimol_area));
gscatter(fluor_roi_stimavg(:,plot_rois(1)),fluor_roi_stimavg(:,plot_rois(2)), ...
    muscimol_area_cat);
axis square;
xlabel(wf_roi(plot_rois(1)).area);
ylabel(wf_roi(plot_rois(2)).area);



% (testing: by trial)
trial_muscimol = cellfun(@(x,y) repmat({x},y,1),muscimol_area_cat, ...
    num2cell(trials_recording),'uni',false);
trial_muscimol = vertcat(trial_muscimol{:});

a = squeeze(nanmean(fluor_roi_deconv(:,use_t,:),2));
figure;gscatter(a(:,1),a(:,7),trial_muscimol)








