%% Test analysis for imaging with operant task (over learning)
% (note: this includes things brought from test_corticostriatal and
% test_learning)

%% ~~~~~~~~~ SINGLE-SESSION ~~~~~~~~~

%% Load example session

% animal = 'AP100';
% day = '2021-05-03';
% day = '2021-05-07';

animal = 'AP104';
day = '2021-06-06';
% day = '2021-06-08';
% day = '2021-06-11';
% day = '2021-06-13';

experiment = 1;
verbose = true;
load_parts.imaging = false;
AP_load_experiment;

t = Timeline.rawDAQTimestamps;

%% Get rewarded/unrewarded movement times/epochs

% Plot velocity/reward/stim
wheel_velocity_move = wheel_velocity;
wheel_velocity_move(~wheel_move) = NaN;

figure; hold on;
plot(t,wheel_velocity,'k');
plot(t,wheel_velocity_move,'r');
plot(reward_t_timeline,0,'.b','MarkerSize',10);
plot(t,stimOn_epochs*0.1,'g');

% Get real and null distribution of stim-to-move times
% (get gaps between iti movements and prior movements)
use_wheel_move_iti = wheel_move_iti_idx(wheel_move_iti_idx > 1); % no gap for first movement
iti_move_gaps = wheel_starts(use_wheel_move_iti) - wheel_stops(use_wheel_move_iti - 1);
% (get gaps between stim movements and prior movements)
use_wheel_move_stim = wheel_move_stim_idx(wheel_move_stim_idx > 1); % no gap for first movement
wheel_move_preresponse_stop = wheel_stops(use_wheel_move_stim-1);
% (make null distribution: positive prior movement + gap relative to stim)
% (from 2nd stim, in case no move before first)
n_shuff = 10000;
stim_to_move_shuff = ...
    stimOn_times(end-length(use_wheel_move_stim)+1:end) ...
    - (wheel_move_preresponse_stop + ...
    reshape(randsample(iti_move_gaps,length(wheel_move_preresponse_stop)*n_shuff,true),[],n_shuff));
stim_to_move_shuff(stim_to_move_shuff < 0) = NaN; % positives only, negatives would've been resets

% (TESTING: get move start time after and end time before reward)
iti_move_gaps = wheel_starts(wheel_move_iti_idx(2:length(wheel_move_iti_idx))) - ...
    wheel_starts(wheel_move_iti_idx(2:length(wheel_move_iti_idx))-1);

wheel_move_preresponse_stop = wheel_stops(wheel_move_stim_idx(2:length(wheel_move_stim_idx))-1);

r = wheel_starts(wheel_move_stim_idx(2:end)) - wheel_move_preresponse_stop;


a = wheel_starts(wheel_move_stim_idx(2:length(wheel_move_stim_idx))-1);
b = wheel_stops(wheel_move_stim_idx(2:length(wheel_move_stim_idx))-2);
c = a - b;
figure; hold on;
histogram(stim_to_move_shuff(:),0:0.02:1,'EdgeColor','none','normalization','probability')
histogram(stim_to_move,0:0.02:1,'EdgeColor','none','normalization','probability')


% (align movement to prestim wheel stop)
figure;subplot(1,2,1);
% use_align = wheel_move_preresponse_stop;
use_align = wheel_stops(wheel_move_stim_idx-2);
surround_t_centers = -10:0.1:10;
surround_times = use_align + surround_t_centers;
surround_move = interp1(t,wheel_move,surround_times,'previous');
imagesc(surround_move);
colormap(brewermap([],'Greys'));

% (align movement to second-last stop vs stim)
subplot(1,2,2);
use_align = stimOn_times;
surround_t_centers = -10:0.1:10;
surround_times = use_align + surround_t_centers;
surround_move = interp1(t,wheel_move,surround_times,'previous');
imagesc(surround_move);
colormap(brewermap([],'Greys'));


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

% Get ITI movements (when stim is off)
[wheel_velocity,wheel_move,wheel_velocity_split] = AP_parse_wheel(wheel_position,Timeline.hw.daqSampleRate);
t_movesplit = mat2cell(Timeline.rawDAQTimestamps',cellfun(@length,wheel_velocity_split));

stimOff_times = signals_events.stimOffTimes';
stim_on_epochs = interp1([stimOn_times;stimOff_times], ...
    [ones(size(stimOn_times));zeros(size(stimOff_times))], ...
    Timeline.rawDAQTimestamps','previous','extrap');

iti_move_epochs = cellfun(@(move,stim) all(move) & ~any(stim), ...
    mat2cell(wheel_move,cellfun(@length,wheel_velocity_split)), ...
    mat2cell(stim_on_epochs,cellfun(@length,wheel_velocity_split)));
iti_move_starts = cellfun(@(x) x(1), t_movesplit(iti_move_epochs));
    
% PSTH with rewarded & ITI move starts
[~,rxn_sort_idx] = sort(stim_to_move);
[~,iti_move_sort_idx] = sort(cellfun(@(x) sum(abs(x)),t_movesplit(iti_move_epochs)));

AP_cellraster({stimOn_times,wheel_move_time,iti_move_starts,reward_t_timeline}, ...
    {rxn_sort_idx,rxn_sort_idx,iti_move_sort_idx,1:length(reward_t_timeline)})


%% Get average rewarded and other movements

wheel_split_n = cellfun(@length,wheel_velocity_split);
wheel_pos_split = mat2cell(wheel_position,wheel_split_n);

wheel_split_n = cellfun(@length,wheel_velocity_split)';
split_move = cellfun(@(x) x(1),mat2cell(wheel_move,wheel_split_n));

wheel_gain = 8; % (deg/mm - in signals protocol)

max_move = max(wheel_split_n(split_move == 1));
vel_pad = cell2mat(cellfun(@(x) padarray(x,max_move-length(x),NaN,'post'), ...
    wheel_velocity_split(split_move == 1),'uni',false)');
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

%% Stim-surround movement

t = Timeline.rawDAQTimestamps;

% (probability of any movement around stim)
stim_surround_t_centers = -10:0.1:10;
stim_surround_times = stimOn_times + stim_surround_t_centers;
stim_surround_move = interp1(t,wheel_move,stim_surround_times,'previous');

figure;
subplot(1,3,1);
imagesc(stim_surround_t_centers,[],stim_surround_move);
ylabel('Trial order');

subplot(1,3,2);
trial_quiescence = signals_events.trialQuiescenceValues(1:n_trials);
[~,sort_idx] = sort(trial_quiescence);
imagesc(stim_surround_t_centers,[],stim_surround_move(sort_idx,:));
hold on;
plot(-trial_quiescence(sort_idx),1:length(sort_idx),'.r')
plot(stim_to_feedback(sort_idx),1:length(sort_idx),'.c')
ylabel('Sorted by quiescence')

subplot(1,3,3);
[~,sort_idx] = sort(stim_to_move);
imagesc(stim_surround_t_centers,[],stim_surround_move(sort_idx,:));
hold on;
plot(-trial_quiescence(sort_idx),1:length(sort_idx),'.r')
plot(stim_to_feedback(sort_idx),1:length(sort_idx),'.c')
ylabel('Sorted by rxn')

colormap(brewermap([],'Greys'));

%% Move-surround movement (and stim)

t = Timeline.rawDAQTimestamps;

% wheel start surround
surround_t_centers = -10:0.1:10;
surround_times = wheel_starts + surround_t_centers;
surround_move = interp1(t,wheel_move,surround_times,'previous');

figure;
subplot(1,2,1);
imagesc(surround_t_centers,[],surround_move(wheel_move_iti_idx,:));
title('ITI movements')

subplot(1,2,2);
imagesc(surround_t_centers,[],surround_move(wheel_move_stim_idx,:));
hold on;
plot(-stim_to_move,1:length(stim_to_move),'.r');
title('Stim movements')

colormap(brewermap([],'Greys'));

% wheel stop surround
surround_t_centers = 0:0.1:5;
surround_times = wheel_stops + surround_t_centers;
surround_move = interp1(t,wheel_move,surround_times,'previous');

figure;
subplot(1,2,1);
imagesc(surround_t_centers,[],surround_move(wheel_move_iti_idx(2:end)-1,:));
title('ITI movements')

subplot(1,2,2);
imagesc(surround_t_centers,[],surround_move(wheel_move_stim_idx-1,:));
hold on;
title('Stim movements')

colormap(brewermap([],'Greys'));



%% Move time regression (Cox proportional hazards)

t = Timeline.rawDAQTimestamps;

% Plot velocity/reward/stim
wheel_velocity_move = wheel_velocity;
wheel_velocity_move(~wheel_move) = NaN;

figure; hold on;
plot(t,wheel_velocity,'k');
plot(t,wheel_velocity_move,'r');
plot(reward_t_timeline,0,'.b','MarkerSize',10);
plot(t,stimOn_epochs*0.1,'--g');


% (carried over from AP_operant_behavior)
t_since_reward = t - ...
    interp1(reward_t_timeline,reward_t_timeline,t,'previous','extrap');

t_since_move = t - ...
    interp1(wheel_stops,wheel_stops,t,'previous','extrap');

t_since_stim = t - ...
    interp1(stimOn_times,stimOn_times,t,'previous','extrap');

% r = [t_since_reward;t_since_move;t_since_stim];
% s = diff([0;wheel_move(1:1000:end)]');s = s == 1;
% t_lags = [-500:500];
% [k,ps,ev] = AP_regresskernel(r(:,1:1000:end),s,t_lags,0);
% figure;plot(t_lags,k');

a = wheel_starts(2:end)-wheel_stops(1:end-1);
b = a(wheel_move_stim_idx-1);

r_idx = (wheel_move_stim_idx-1)+[-10:-1];
use_r = all(r_idx>0,2);

% [k,ps,ev] = AP_regresskernel(a(r_idx(use_r,:)),b(use_r),0,0,[],5,0,1);

x = [a(r_idx(use_r,:)),stim_to_move(use_r),ones(sum(use_r),1)]\b(use_r);


a = diff([0;wheel_move]) == 1;
ar = t_since_reward(a);
am = t_since_move(a);
as = t_since_stim(a);
figure; hold on;
plot3(ar,am,as,'.k')
plot3(ar(wheel_move_stim_idx),am(wheel_move_stim_idx),as(wheel_move_stim_idx),'or')

% Regress time from stim to move relative to previous gaps
use_last_gaps = 5;
wheel_gaps = wheel_starts(2:end)-wheel_stops(1:end-1);
wheel_gaps_regressor_idx = (wheel_move_stim_idx-1)+[-use_last_gaps:-1];
use_move_regressors = all(wheel_gaps_regressor_idx>0,2);
wheel_gaps_regressors = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors,:));

[b,logl,H,stats] = coxphfit(wheel_gaps_regressors,stim_to_move(use_move_regressors));
figure;errorbar(stats.beta,stats.se);
xlabel('Previous move gap');
ylabel('Cox \beta');


% Regress time from move to post-stim move relative to pre-stim gap/stim
use_last_gaps = 5;
wheel_gaps = wheel_starts(2:end)-wheel_stops(1:end-1);
wheel_gaps_regressor_idx = (wheel_move_stim_idx-1)+[-use_last_gaps:0];
use_move_regressors = all(wheel_gaps_regressor_idx>0,2);

wheel_gap_predict = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors,end));
wheel_gaps_regressors = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors,1:end-1));
wheel_stim_regressor = stimOn_times(use_move_regressors) - wheel_stops(wheel_move_stim_idx(use_move_regressors)-1);

[b,logl_stim,H,stats] = coxphfit([wheel_gaps_regressors,wheel_stim_regressor],wheel_gap_predict);
figure;errorbar(stats.beta,stats.se);
xlabel('Previous move gap');
ylabel('Cox \beta');

[b_nostim,logl_nostim,H_nostim,stats_nostim] = coxphfit([wheel_gaps_regressors],wheel_gap_predict);

disp(logl_stim-logl_nostim);


% Regress time from move to ANY move (with/without stim info)
use_last_gaps = 5;
wheel_gaps = wheel_starts(2:end)-wheel_stops(1:end-1);
wheel_gaps_regressor_idx = (1:length(wheel_gaps))'+[-use_last_gaps:0];
use_move_regressors = all(wheel_gaps_regressor_idx>0,2);

wheel_gap_predict = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors,end));
wheel_gaps_regressors = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors,1:end-1));

[b,logl,H,stats] = coxphfit([wheel_gaps_regressors],wheel_gap_predict);
figure;errorbar(stats.beta,stats.se);
xlabel('Previous move gap');
ylabel('Cox \beta');

a = stimOn_times - wheel_stops(wheel_move_stim_idx-1);
idx = ismember(2:length(wheel_starts),wheel_move_stim_idx)';
m = zeros(size(idx));
m(idx) = a;

[b_stim,logl_stim,H,stats_stim] = coxphfit([wheel_gaps_regressors,m(use_move_regressors)],wheel_gap_predict);
figure;errorbar(stats_stim.beta,stats_stim.se);
xlabel('Previous move gap');
ylabel('Cox \beta');

disp(logl_stim-logl);


% Regress ITI gap time move only (sanity check?)
use_last_gaps = 5;
wheel_gaps = wheel_starts(2:end)-wheel_stops(1:end-1);
wheel_gaps_regressor_idx = (1:length(wheel_gaps))'+[-use_last_gaps:0];
use_move_regressors = all(wheel_gaps_regressor_idx>0,2);

iti_idx = ismember(2:length(wheel_starts),wheel_move_iti_idx)';
stim_idx = ismember(2:length(wheel_starts),wheel_move_stim_idx)';

wheel_gap_predict = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors & iti_idx,end));
wheel_gaps_regressors = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors & iti_idx,1:end-1));

[b,logl,H,stats] = coxphfit([wheel_gaps_regressors],wheel_gap_predict);
figure;errorbar(stats.beta,stats.se);
xlabel('Previous move gap');
ylabel('Cox \beta');




% (predict ITI movements, test on stim movements)
use_last_gaps = 5;
wheel_gaps = wheel_starts(2:end)-wheel_stops(1:end-1);
wheel_gaps_regressor_idx = (1:length(wheel_gaps))'+[-use_last_gaps:0];
use_move_regressors = all(wheel_gaps_regressor_idx>0,2);

iti_idx = ismember(2:length(wheel_starts),wheel_move_iti_idx)';
stim_idx = ismember(2:length(wheel_starts),wheel_move_stim_idx)';

% -> fit on ITI
r = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors & iti_idx,1:end-1));
s = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors & iti_idx,end));
[k,ps,ev] = AP_regresskernel(r',s',0,0,[],1,1,1);
disp(ev)

% -> test on stim
r_test = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors & stim_idx,1:end-1));
s_test = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors & stim_idx,end));

s_test_pred = [r_test,ones(size(r_test,1),1)]*cell2mat(k);
r2 = 1 - (nansum((s_test-s_test_pred).^2,1)./ ...
        nansum((s_test-nanmean(s_test,1)).^2,1));
disp(r2)




% (predict ITI movements, test on stim movements)
use_last_gaps = 1;
wheel_gaps = wheel_starts(2:end)-wheel_stops(1:end-1);
wheel_gaps_regressor_idx = (1:length(wheel_gaps))'+[-use_last_gaps:0];
use_move_regressors = all(wheel_gaps_regressor_idx>0,2);

iti_idx = ismember(2:length(wheel_starts),wheel_move_iti_idx)';
stim_idx = ismember(2:length(wheel_starts),wheel_move_stim_idx)';

% -> fit on ITI
r = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors & stim_idx,1:end-1));
s = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors & stim_idx,end));
[k,ps,ev] = AP_regresskernel(r',s',0,0,[]);
disp(ev)

% -> test on stim
r_test = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors & stim_idx,1:end-1));
s_test = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors & stim_idx,end));

s_test_pred = [r_test,ones(size(r_test,1),1)]*cell2mat(k);
r2 = 1 - (nansum((s_test-s_test_pred).^2,1)./ ...
        nansum((s_test-nanmean(s_test,1)).^2,1));
disp(r2)






%% ~~~~~~~~~ PIP'S NAIVE EXPERIMENTS  ~~~~~~~~~

% PC052/3/4/5

% ephys path is something different
ephys_path = '\\128.40.224.65\Subjects\PC052\2021-09-13\ephys\kilosort\imec0';

% params has different name
header_path = [ephys_path filesep 'params.py'];

% sync file is a little different
sync_long = load('\\128.40.224.65\Subjects\PC052\2021-09-13\ephys\kilosort\imec0\sync.mat');
sync(4).timestamps = (find(diff(sync_long.sync) ~= 0)+1)'./ephys_sample_rate;

% (also looks like good/bad sort not done yet)

animal = 'PC052';day = '2021-09-13';experiment = 4;verbose = true;AP_load_experiment;




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


%% Passive - tetO

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


%% Task - tetO

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

%% Muscimol: Passive - tetO

clear all
disp('Muscimol: Passive trial activity (tetO)')

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

%% Muscimol: Task - tetO

clear all
disp('Muscimol: Task trial activity (tetO)')

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
save_fn = ['trial_activity_task_teto_muscimol'];
save([save_path filesep save_fn],'-v7.3');

%% Muscimol: Passive - corticostriatal

clear all
disp('Muscimol: Passive trial activity (cstr)')

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
save_fn = ['trial_activity_passive_cstr_muscimol'];
save([save_path filesep save_fn],'-v7.3');
disp(['Saved: ' save_path filesep save_fn])


%% ~~~~~~~~~ BATCH ANALYSIS ~~~~~~~~~

%% >> Passive 

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_passive_tetO';
% data_fn = 'trial_activity_passive_corticostriatal';

AP_load_trials_wf;

% Get animal and day index for each trial
trial_animal = cell2mat(arrayfun(@(x) ...
    x*ones(size(vertcat(wheel_all{x}{:}),1),1), ...
    [1:length(wheel_all)]','uni',false));
trial_day = cell2mat(cellfun(@(x) cell2mat(cellfun(@(curr_day,x) ...
    curr_day*ones(size(x,1),1),num2cell(1:length(x))',x,'uni',false)), ...
    wheel_all,'uni',false));

trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));

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
use_t = t >= 0.05 & t <= 0.1;
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
unlearned_days = 4:6;
expert_days = 6:10;

stim_avg_naive = nanmean(stim_px_avg_dayavg(:,:,:,naive_days),4);
stim_avg_unlearned = nanmean(stim_px_avg_dayavg(:,:,:,unlearned_days),4);
stim_avg_learned = nanmean(stim_px_avg_dayavg(:,:,:,expert_days),4);

AP_image_scroll([stim_avg_naive,stim_avg_unlearned,stim_avg_learned]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;

use_t = t >= 0.05 & t <= 0.2;
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
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
axis image;
set(gcf,'name',animals{use_animal});


%% Time-average ROI activity by day

curr_act = fluor_roi_deconv(:,:,7);
% curr_act = fluor_roi_deconv(:,:,7) - fluor_roi_deconv(:,:,17);

use_trials = trial_stim_allcat == 1 & quiescent_trials;

use_t = t >= 0.07 & t <= 0.17;
curr_act_timeavg = nanmean(double(curr_act(:,use_t)),2);

roi_animal_day_avg = accumarray([trial_animal(use_trials),trial_day(use_trials)], ...
    curr_act_timeavg(use_trials),[max(trial_animal),max(trial_day)], ...
    @nanmean,NaN);
figure; hold on;
plot([1:max(trial_day)]-3,roi_animal_day_avg');
plot([1:max(trial_day)]-3,nanmean(roi_animal_day_avg,1),'k','linewidth',2);

%%%% split days
n_daysplit = 4;
day_split_idx = cell2mat(arrayfun(@(x) ...
    min(floor(linspace(1,n_daysplit+1,x)),n_daysplit)', ...
    trials_recording,'uni',false));
roi_animal_daysplit_avg = accumarray( ...
    [trial_animal(use_trials),trial_day(use_trials),day_split_idx(use_trials)], ...
    curr_act_timeavg(use_trials),[max(trial_animal),max(trial_day),n_daysplit], ...
    @nanmean,NaN);
figure; hold on;
roi_animal_day_avg_long = reshape(permute( ...
    padarray(roi_animal_daysplit_avg,[0,0,1],NaN,'post'),[3,2,1]),[],length(animals));
% plot([1:size(roi_animal_day_avg_long,1)]/(n_daysplit+1),roi_animal_day_avg_long);
errorbar([1:size(roi_animal_day_avg_long,1)]/(n_daysplit+1) - 3, ...
    nanmean(roi_animal_day_avg_long,2),AP_sem(roi_animal_day_avg_long,2),'k','linewidth',2);



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

% Average response given performance
curr_stim_v_avg_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    nanmean(cell2mat(permute(cellfun(@(x,y) x(:,:,find(y > 0.2)+3,3), ...
    stim_v_avg,move_prepost_max_ratio_animal','uni',false),[1,3,2])),3));
AP_image_scroll(curr_stim_v_avg_px,t)
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image off;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);


% Get average stim (skip first 3 days: no behavior)
stim_t = [0.05,0.2];
stim_px = cellfun(@(x) ...
    squeeze(AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    nanmean(x(:,t >= stim_t(1) & t <= stim_t(2),4:end,3),2))), ...
    stim_v_avg,'uni',false);

stim_roi_stimavg = cellfun(@(x) ...
    squeeze(nanmean(x(:,t >= stim_t(1) & t <= stim_t(2),4:end,3),2)), ...
    stim_roi_avg,'uni',false);

% Plot pixels for each mouse
figure;
for curr_animal = 1:length(animals)
    subplot(1,length(animals),curr_animal);
    imagesc(nanmean(stim_px{curr_animal},3));
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'PrGn'));
    axis image off;
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title(animals{curr_animal});
end

% Correlate behavior measure with pixels

stim_move_t_median = cellfun(@(x) cellfun(@nanmedian,x),{bhv.stim_move_t},'uni',false);
stim_reward_t_median = cellfun(@(x) cellfun(@nanmedian,x),{bhv.stim_reward_t},'uni',false);
stim_move_reward_t_median = cellfun(@(x,y) cellfun(@(x,y) nanmedian(y-x),x,y),{bhv.stim_move_t},{bhv.stim_reward_t},'uni',false);


stim_surround_t = bhv(curr_animal).stim_surround_t;
poststim_move_frac_max = cellfun(@(x) cellfun(@(x) ...
    max(nanmean(x(:,stim_surround_t > 0),1)),x), ...
    {bhv.stim_surround_wheel},'uni',false);

use_bhv = stim_move_t_median;

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

%% Average passive by behavior

% Load and prep behavior
data_fn_parse = regexp(data_fn,'_','split');
bhv_fn = [trial_data_path filesep 'bhv_' data_fn_parse{end}];
load(bhv_fn);

% Load muscimol injection info
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
muscimol_fn = [data_path filesep 'muscimol.mat'];
load(muscimol_fn);

% Use days before muscimol
use_days = cell(size(bhv));
for curr_animal = 1:length(bhv)
    muscimol_animal_idx = ismember({muscimol.animal},bhv(curr_animal).animal);
    if isempty(muscimol_animal_idx)
        continue
    end
    
    pre_muscimol_day_idx = datenum(bhv(curr_animal).day) < ...
        datenum(muscimol(muscimol_animal_idx).day(1));
    muscimol_day_idx = ismember(datenum(bhv(curr_animal).day), ...
        datenum(muscimol(muscimol_animal_idx).day));
    
    use_days{curr_animal} = pre_muscimol_day_idx;    
end

% Stim to move time (NOT SAME NUM TRIALS - REWARDED ONLY)
stim_move_t = cell2mat(cellfun(@(x,y) cell2mat(x(y)), ...
    {bhv.stim_move_t},use_days,'uni',false))';

% Pre/post move probability
stim_surround_wheel_avg = cell2mat(cellfun(@(x,y) ...
    cell2mat(cellfun(@(x) nanmean(x,1),x(y)','uni',false)), ...
    {bhv.stim_surround_wheel},use_days,'uni',false)');

stim_surround_t = bhv(1).stim_surround_t;
move_prestim_max = max(stim_surround_wheel_avg(:,stim_surround_t<0,:),[],2);
move_poststim_max = max(stim_surround_wheel_avg(:,stim_surround_t>=0,:),[],2);
move_prepost_max_ratio = ...
    (move_poststim_max-move_prestim_max)./(move_poststim_max+move_prestim_max);

move_prepost_max_ratio_padcat = ...
    cell2mat(cellfun(@(x) padarray(x,[3,0],NaN,'pre'), ...
    mat2cell(move_prepost_max_ratio,cellfun(@sum,use_days)),'uni',false));

% Group stim responses by behavioral performance
stim_v_avg_cat = cell2mat(permute(stim_v_avg,[1,3,2,4])); 

move_prepost_max_ratio_padcat_grp = ...
    discretize(move_prepost_max_ratio_padcat,[-Inf,0.2,Inf]);
move_prepost_max_ratio_padcat_grp(isnan(move_prepost_max_ratio_padcat_grp)) = 0;
stim_v_avg_cat_grp = reshape(grpstats( ...
    reshape(permute(stim_v_avg_cat,[1,2,4,3]),[],size(stim_v_avg_cat,3))', ...
    move_prepost_max_ratio_padcat_grp)',n_vs,length(t),length(stim_unique),[]);

curr_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),stim_v_avg_cat_grp);
AP_image_scroll(reshape(permute(curr_px,[1,2,4,3,5]), ...
    size(curr_px,1),[],size(curr_px,3),size(curr_px,5)),t);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5],[], ...
    [size(U_master,1),size(U_master,2),1,size(curr_px,5)]);

% Plot time-average
use_t = t >= 0.07 & t <= 0.17;
curr_px_avg = squeeze(nanmean(curr_px(:,:,use_t,:,:),3));
c = prctile(curr_px_avg(:),99.5).*[-1,1];
figure;
for curr_stim = 1:length(stim_unique)
    for curr_learn = 1:size(curr_px_avg,3)
        subplot(length(stim_unique),size(curr_px_avg,3), ...
            size(curr_px_avg,3)*(curr_stim-1) + curr_learn);
        imagesc(curr_px_avg(:,:,curr_stim,curr_learn));
        caxis(c);
        colormap(brewermap([],'PrGn'));
        axis image off
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    end
end

% Get fluorescence in ROIs
stim_v_avg_cat_grp_roi = ...
    AP_svd_roi(U_master(:,:,1:n_vs),stim_v_avg_cat_grp,[],[], ...
    cat(3,wf_roi.mask));
    





%% >> (Choiceworld task: post-learn only)

trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\paper_unused';
data_fn = 'trial_activity_choiceworld_wfonly';

AP_load_trials_wf;


%% >> Task

% Task
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_task_teto';
% data_fn = 'trial_activity_task_corticostriatal';

AP_load_trials_wf;

% Get animal and day index for each trial
trial_animal = cell2mat(arrayfun(@(x) ...
    x*ones(size(vertcat(wheel_all{x}{:}),1),1), ...
    [1:length(wheel_all)]','uni',false));
trial_day = cell2mat(cellfun(@(x) cell2mat(cellfun(@(curr_day,x) ...
    curr_day*ones(size(x,1),1),num2cell(1:length(x))',x,'uni',false)), ...
    wheel_all,'uni',false));

% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_animal;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));


% Load and prep behavior
data_fn_parse = regexp(data_fn,'_','split');
bhv_fn = [trial_data_path filesep 'bhv_' data_fn_parse{end}];
load(bhv_fn);

% Load muscimol injection info
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
muscimol_fn = [data_path filesep 'muscimol.mat'];
load(muscimol_fn);

% Use days before muscimol
use_days = cell(size(bhv));
for curr_animal = 1:length(bhv)
    muscimol_animal_idx = ismember({muscimol.animal},bhv(curr_animal).animal);
    if isempty(muscimol_animal_idx)
        continue
    end
    
    pre_muscimol_day_idx = datenum(bhv(curr_animal).day) < ...
        datenum(muscimol(muscimol_animal_idx).day(1));
    muscimol_day_idx = ismember(datenum(bhv(curr_animal).day), ...
        datenum(muscimol(muscimol_animal_idx).day));
    
    use_days{curr_animal} = pre_muscimol_day_idx;    
end

% Stim to move time (NOT SAME NUM TRIALS - REWARDED ONLY)
stim_move_t = cell2mat(cellfun(@(x,y) cell2mat(x(y)), ...
    {bhv.stim_move_t},use_days,'uni',false))';

% Pre/post move probability
stim_surround_wheel_avg = cell2mat(cellfun(@(x,y) ...
    cell2mat(cellfun(@(x) nanmean(x,1),x(y)','uni',false)), ...
    {bhv.stim_surround_wheel},use_days,'uni',false)');

stim_surround_t = bhv(1).stim_surround_t;
move_prestim_max = max(stim_surround_wheel_avg(:,stim_surround_t<0,:),[],2);
move_poststim_max = max(stim_surround_wheel_avg(:,stim_surround_t>=0,:),[],2);
move_prepost_max_ratio = ...
    (move_poststim_max-move_prestim_max)./(move_poststim_max+move_prestim_max);

move_prepost_max_ratio_animal = mat2cell(move_prepost_max_ratio,cellfun(@sum,use_days));

move_prepost_max_ratio_allcat = cell2mat(cellfun(@(x,y) repmat(x,y,1), ...
    num2cell(move_prepost_max_ratio),num2cell(trials_recording),'uni',false));



%% Average trial activity by day

% % (move-align fluor?)
% fluor_allcat_deconv_move = fluor_allcat_deconv;
% fluor_roi_deconv_move = fluor_roi_deconv;
% t_leeway = -t(1);
% leeway_samples = round(t_leeway*(sample_rate));
% for i = 1:size(fluor_allcat_deconv,1)
%     fluor_allcat_deconv_move(i,:,:,:) = circshift(fluor_allcat_deconv_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
%     fluor_roi_deconv_move(i,:,:,:) = circshift(fluor_roi_deconv_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
% end

% (package back into animal/day)
n_days = cellfun(@length,wheel_all);
fluor_deconv_animal = mat2cell(mat2cell(fluor_allcat_deconv,trials_recording,length(t),n_vs),n_days);
fluor_taskpred_reduced_animal = mat2cell(mat2cell(fluor_taskpred_reduced_allcat, ...
    trials_recording,length(t),n_vs,length(task_regressor_labels)),n_days);

move_t_animal = mat2cell(mat2cell(move_t,trials_recording),n_days);

% Average activity within animal/day
max_days = max(n_days);
fluor_deconv_cat = nan(n_vs,length(t),max_days,length(animals));
for curr_animal = 1:length(animals)
    % (all trials)
    use_trials = cellfun(@(x) true(size(x)),move_t_animal{curr_animal},'uni',false);
    % (reaction-time limited)
%         use_trials = cellfun(@(x) x > 0.1 & x < 0.5,move_t_animal{curr_animal},'uni',false);
    % (trial subset from start or end)
    %     use_trials = cellfun(@(x) find(x > 0.05,10,'first'),move_t_animal{curr_animal},'uni',false);
    
    % (full activity)
%     fluor_deconv_cat(:,:,1:n_days(curr_animal),curr_animal) = ...
%         permute(cell2mat(cellfun(@(x,trials) nanmean(x(trials,:,:),1), ...
%         fluor_deconv_animal{curr_animal},use_trials,'uni',false)),[3,2,1]);
    
    % (reduced activity)
    fluor_deconv_cat(:,:,1:n_days(curr_animal),curr_animal) = ...
        permute(cell2mat(cellfun(@(act,act_reduced,trials) ...
        nanmean(act(trials,:,:) - act_reduced(trials,:,:,1),1), ...
        fluor_deconv_animal{curr_animal},fluor_taskpred_reduced_animal{curr_animal}, ...
        use_trials,'uni',false)),[3,2,1]);
        
end

% Plot average for each day
curr_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),nanmean(fluor_deconv_cat,4));
AP_image_scroll(curr_px,t);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'))
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

% Plot animal across days
use_animal = 5;
curr_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),fluor_deconv_cat(:,:,:,use_animal));
AP_image_scroll(curr_px,t);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'))





%% Average trial activity by rxn/bhv

% Average all activity by reaction time depending on learned
% rxn_bins = [-0.05:0.02:1]';
% rxn_bins = [-Inf,Inf];
% rxn_bins = [0.1,0.5];
rxn_bins = [0.1:0.1:0.5];
rxn_bin_centers = rxn_bins(1:end-1) + diff(rxn_bins)./2;
move_t_discretize = discretize(move_t,rxn_bins);
fluor_rxn = nan(n_vs,length(t),length(rxn_bins)-1);
for curr_rxn = 1:length(rxn_bins)-1
        use_trials = move_t_discretize == curr_rxn & ...
            trial_outcome_allcat == 1 & ...
            move_prepost_max_ratio_allcat >= 0.2;
        
%     use_trials = move_t_discretize == curr_rxn & ...
%         trial_outcome_allcat == 1 & ...
%         trial_stim_allcat == 1;
     
%     % (raw activity)
%     fluor_rxn(:,:,curr_rxn) = ...
%         permute(nanmean(fluor_allcat_deconv(use_trials,:,:),1),[3,2,1]);
    
    % (activity - taskpred_reduced)
        fluor_rxn(:,:,curr_rxn) = ...
            permute(nanmean(fluor_allcat_deconv(use_trials,:,:) - ...
            fluor_taskpred_reduced_allcat(use_trials,:,:,1),1),[3,2,1]);
        
end

curr_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),fluor_rxn);
AP_image_scroll(curr_px,t);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'))
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);


%% Reduced activity by behavior (animal day average)

% Get average fluorescence by animal, day, stim
stim_unique = unique(trial_stim_allcat);
stim_v_avg = cell(size(animals));
for curr_animal = 1:length(animals)        
    for curr_day = 1:max(trial_day(trial_animal == curr_animal))
        for curr_stim_idx = 1:length(stim_unique)
            use_trials = ...
                move_t >= 0.1 & ...
                trial_animal == curr_animal & ...
                trial_day == curr_day & ...
                trial_stim_allcat == stim_unique(curr_stim_idx);
            stim_v_avg{curr_animal}(:,:,curr_day,curr_stim_idx) = ...
                permute(nanmean( ...
                fluor_allcat_deconv(use_trials,:,:) - ...
                fluor_taskpred_reduced_allcat(use_trials,:,:,1),1),[3,2,1]);
        end       
    end
end

% Group stim responses by behavioral performance
stim_v_avg_cat = cell2mat(permute(stim_v_avg,[1,3,2,4])); 

move_prepost_max_ratio_grp = ...
    discretize(move_prepost_max_ratio,[-Inf,0.2,Inf]);
move_prepost_max_ratio_grp(isnan(move_prepost_max_ratio)) = 0;
stim_v_avg_cat_grp = reshape(grpstats( ...
    reshape(stim_v_avg_cat,[],size(stim_v_avg_cat,3))', ...
    move_prepost_max_ratio_grp)',n_vs,length(t),[]);

curr_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),stim_v_avg_cat_grp);
AP_image_scroll(curr_px,t);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

% Plot time-average
use_t = t >= 0.07 & t <= 0.17;
curr_px_avg = squeeze(nanmean(curr_px(:,:,use_t,:,:),3));
c = prctile(curr_px_avg(:),99.5).*[-1,1];
figure;
for curr_learn = 1:size(curr_px_avg,3)
    subplot(1,size(curr_px_avg,3),curr_learn);
    imagesc(curr_px_avg(:,:,curr_learn));
    caxis(c);
    colormap(brewermap([],'PrGn'));
    axis image off
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
end




%% Average trial activity by wheel velocity

use_trials = trial_outcome_allcat == 1 & move_prepost_max_ratio_allcat > 0.2 & move_t > 0.1 & move_t < 0.2;

wheel_vel_summed = sum(abs(wheel_allcat(:,t > 0 & t < 0.5)),2);

wheel_vel_bins = prctile(wheel_vel_summed(use_trials),[0:10:100]);
wheel_vel_discretize = discretize(wheel_vel_summed,wheel_vel_bins);

fluor_binned = nan(n_vs,length(t),length(wheel_vel_bins)-1);
for curr_bin = 1:length(wheel_vel_bins)-1
    curr_trials = use_trials & wheel_vel_discretize == curr_bin;
    % (raw activity)
    fluor_binned(:,:,curr_bin) = ...
        permute(nanmean(fluor_allcat_deconv(curr_trials,:,:),1),[3,2,1]);
    % (task reduced activity)
%     fluor_binned(:,:,curr_bin) = ...
%         permute(nanmean(fluor_allcat_deconv(curr_trials,:,:) - ...
%         fluor_taskpred_reduced_allcat(curr_trials,:,:,1),1),[3,2,1]);
end

curr_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),fluor_binned);
AP_image_scroll(curr_px,t);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'))
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);


%% Average ROI activity by reaction time

use_trials = trial_outcome_allcat == 1 & move_prepost_max_ratio_allcat > 0.2;

% curr_act = fluor_roi_deconv(use_trials,:,7);
% curr_act = fluor_roi_taskpred(use_trials,:,1);
% curr_act = fluor_roi_deconv(use_trials,:,7) - fluor_roi_deconv(use_trials,:,17);
% curr_act = fluor_roi_taskpred(use_trials,:,4) - fluor_roi_taskpred(use_trials,:,14);
curr_act = fluor_roi_deconv(use_trials,:,7) - fluor_roi_taskpred_reduced(use_trials,:,7,1);

rxn_bins = [-0.05:0.01:0.5]';
rxn_bin_centers = rxn_bins(1:end-1) + diff(rxn_bins)./2;
move_t_discretize = discretize(move_t(use_trials),rxn_bins);

% By each time point
curr_act_rxn = accumarray( ...
    {reshape(repmat(1:length(t),sum(~isnan(move_t_discretize)),1),[],1), ...
    reshape(repmat(move_t_discretize(~isnan(move_t_discretize)),1,length(t)),[],1)}, ...
    reshape(curr_act(~isnan(move_t_discretize),:),[],1), ...
    [length(t),length(rxn_bin_centers)],@nanmean,single(NaN));

figure;imagesc(t,rxn_bin_centers,curr_act_rxn');
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
set(gca,'YDir','normal')
xlabel('Activity time from stim');
ylabel('Reaction time');

%% Time-average ROI activity by day/behavior

% curr_act = fluor_roi_deconv(:,:,7);
% curr_act = fluor_roi_deconv(:,:,7) - fluor_roi_deconv(:,:,17);
curr_act = fluor_roi_deconv(:,:,7) - fluor_roi_taskpred_reduced(:,:,7,1);

use_t = t >= 0.07 & t <= 0.17;
curr_act_timeavg = nanmean(double(curr_act(:,use_t)),2);

use_trials = true(size(move_t));
% use_trials = move_t >= 0.1;

roi_animal_day_avg = accumarray([trial_animal(use_trials),trial_day(use_trials)], ...
    curr_act_timeavg(use_trials),[max(trial_animal),max(trial_day)], ...
    @nanmean,NaN);
figure; hold on;
plot(roi_animal_day_avg');
plot(nanmean(roi_animal_day_avg,1),'k','linewidth',2);

%%%% split days
n_daysplit = 6;
day_split_idx = cell2mat(arrayfun(@(x) ...
    min(floor(linspace(1,n_daysplit+1,x)),n_daysplit)', ...
    trials_recording,'uni',false));
roi_animal_daysplit_avg = accumarray( ...
    [trial_animal(use_trials),trial_day(use_trials),day_split_idx(use_trials)], ...
    curr_act_timeavg(use_trials),[max(trial_animal),max(trial_day),n_daysplit], ...
    @nanmean,NaN);
figure; hold on;
roi_animal_daysplit_avg_long = reshape(permute( ...
    padarray(roi_animal_daysplit_avg,[0,0,1],NaN,'post'),[3,2,1]),[],length(animals));
% plot([1:size(roi_animal_day_avg_long,1)]/(n_daysplit+1),roi_animal_day_avg_long);
errorbar([1:size(roi_animal_daysplit_avg_long,1)]/(n_daysplit+1), ...
    nanmean(roi_animal_daysplit_avg_long,2),AP_sem(roi_animal_daysplit_avg_long,2),'k','linewidth',2);

% Split by pre/post move ratio
moveratio_bins = [-1:0.1:1];
moveratio_bin_centers = moveratio_bins(1:end-1)+diff(moveratio_bins)/2;
moveratio_discretize = discretize(move_prepost_max_ratio_allcat,moveratio_bins);

rxn_bins = [0:0.05:0.7];
rxn_bin_centers = rxn_bins(1:end-1) + diff(rxn_bins)./2;
rxn_discretize = discretize(move_t,rxn_bins);

use_trials = ~isnan(rxn_discretize);

curr_act_binned = accumarray( ...
    [moveratio_discretize(use_trials),rxn_discretize(use_trials)], ...
    curr_act_timeavg(use_trials), ...
    [length(moveratio_bin_centers),length(rxn_bin_centers)],@nanmean,NaN);


figure;imagesc(rxn_bin_centers,moveratio_bin_centers,curr_act_binned);
xlabel('Reaction time');
ylabel('Move ratio');
colormap(hot);



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

% Get average kernel given behavior
move_prepost_max_ratio_padcat = ...
    cell2mat(cellfun(@(x) padarray(x,max_days-length(x),nan,'post'), ...
    move_prepost_max_ratio_animal,'uni',false)');
use_prepost_ratio = move_prepost_max_ratio_padcat >= 0.2;
% use_prepost_ratio = move_prepost_max_ratio_padcat < 0.2 & circshift(move_prepost_max_ratio_padcat,[-1,0]) > 0.2;
% use_prepost_ratio = move_prepost_max_ratio_padcat > 0.2 & circshift(move_prepost_max_ratio_padcat,[1,0]) < 0.2;
curr_k = squeeze(fluor_taskpred_k_cat(1,:,:));
curr_k_avg_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    permute(nanmean(vertcat(curr_k{use_prepost_ratio}),1),[3,2,1]));
AP_image_scroll(curr_k_avg_px);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
axis image;

% Get kernel weights in ROIs
fluor_taskpred_k_cat_roi = cellfun(@(x) ...
    AP_svd_roi(U_master(:,:,1:n_vs),permute(x,[3,2,1]),[],[],cat(3,wf_roi.mask)), ...
    fluor_taskpred_k_cat,'uni',false);

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
curr_regressor = 1;
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
curr_k = cell2mat(cellfun(@(x) x{1}(1,:,:),fluor_taskpred_k_all{use_animal},'uni',false));
curr_k_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),permute(curr_k,[3,2,1]));
AP_image_scroll(curr_k_px);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;
set(gcf,'name',animals{use_animal})


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

%% >> Muscimol: task

trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_task_teto_muscimol';

AP_load_trials_wf;

% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_animal;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));

% Get task>cortex parameters
n_regressors = length(task_regressor_labels);
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

% Average kernels by muscimol area
muscimol_area_cat = horzcat(muscimol_area{:})';
unique_muscimol_area = unique(muscimol_area_cat);

muscimol_k = cell(n_regressors,1);
for curr_area = 1:length(unique_muscimol_area)
    curr_exps = cellfun(@(x) strcmp(x,unique_muscimol_area{curr_area}), ...
        muscimol_area,'uni',false)';
    curr_k = cellfun(@(x,y) x(y),fluor_taskpred_k_all,curr_exps,'uni',false);
    curr_k_cat = vertcat(curr_k{:});
    for curr_regressor = 1:n_regressors
        muscimol_k{curr_regressor}(:,:,:,curr_area) = ...
            nanmean(cell2mat(permute(cellfun(@(x) x{curr_regressor}, ...
            curr_k_cat,'uni',false),[2,3,4,1])),4);
    end
end

curr_regressor = 1;
curr_subregressor = 1;
curr_k_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ... 
    permute(muscimol_k{curr_regressor}(curr_subregressor,:,:,:),[3,2,4,1]));
AP_image_scroll(curr_k_px);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

% Plot kernel over days in single animal
use_animal = 3;
curr_regressor = 1;
curr_subregressor = 1;

curr_k = cell2mat(cellfun(@(x) x{curr_regressor}(curr_subregressor,:,:), ...
    fluor_taskpred_k_all{use_animal},'uni',false));
curr_k_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),permute(curr_k,[3,2,1]));
AP_image_scroll(curr_k_px);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;


%% >> Muscimol: passive

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
% data_fn = 'trial_activity_passive_teto_muscimol';
data_fn = 'trial_activity_passive_cstr_muscimol';

AP_load_trials_wf;

% (package deconv back into animal/day)
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));

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

unique_muscimol_area = flipud(unique(muscimol_area_cat));
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
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

% Plot time averages and ROIs
use_t = t >= 0.07 & t <= 0.17;
fluor_roi_stimavg = squeeze(nanmean(fluor_roi_exp_stimavg(:,use_t,:),2));
fluor_muscimol_px_stimavg = squeeze(nanmean(fluor_muscimol_avg_px(:,:,use_t,:),3));

figure;
c = prctile(fluor_muscimol_px_stimavg(:),100).*[-1,1];
for curr_area = 1:length(unique_muscimol_area)
    subplot(1,length(unique_muscimol_area)+1,curr_area);
    imagesc(fluor_muscimol_px_stimavg(:,:,curr_area));
    caxis(c);
    colormap(brewermap([],'PrGn'));
    axis image off;
    AP_reference_outline('ccf_aligned','k');
    title(unique_muscimol_area{curr_area});
end

subplot(1,length(unique_muscimol_area)+1,length(unique_muscimol_area)+1);
plot_rois = [1,7];
area_col = lines(length(unique_muscimol_area));
gscatter(fluor_roi_stimavg(:,plot_rois(1)),fluor_roi_stimavg(:,plot_rois(2)), ...
    muscimol_area_cat,area_col);
axis square;
xlabel(wf_roi(plot_rois(1)).area);
ylabel(wf_roi(plot_rois(2)).area);



% (testing: by trial)
trial_muscimol = cellfun(@(x,y) repmat({x},y,1),muscimol_area_cat, ...
    num2cell(trials_recording),'uni',false);
trial_muscimol = vertcat(trial_muscimol{:});

a = squeeze(nanmean(fluor_roi_deconv(:,use_t,:),2));
figure;gscatter(a(:,1),a(:,7),trial_muscimol)








