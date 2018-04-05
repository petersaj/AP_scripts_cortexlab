%% Single behavior (OLD)

animal = 'AP017';
protocol = 'vanillaChoiceworld';
experiments = AP_find_experiments(animal,protocol);
bhv = struct;

for curr_day = 1:length(experiments)
    
    day = experiments(curr_day).day;
    experiment_num = experiments(curr_day).experiment;
    
    if length(experiment_num) > 1
        warning('NOT USING ALL EXPERIMENTS AT THE MOMENT')
    end
    
    for curr_experiment = length(experiment_num);
        
        experiment = experiment_num(curr_experiment);
        
        [block_filename, block_exists] = AP_cortexlab_filename(animal,day,experiment,'block');
        
        % Load the block file
        load(block_filename)
        
        % Time of session (in minutes)
        session_duration = block.duration/60;
        
        % Trial counts
        n_trials = length(block.paramsValues);
        total_water = sum(block.outputs.rewardValues);
        
        % Wheel movements/biases
        wheel_movement = diff(block.inputs.wheelValues);
        left_wheel_movement = abs(wheel_movement.*(wheel_movement < 0));
        right_wheel_movement = abs(wheel_movement.*(wheel_movement > 0));
        wheel_bias = (sum(right_wheel_movement)-sum(left_wheel_movement))/ ...
            (sum(right_wheel_movement)+sum(left_wheel_movement));
        
        % Performance (note that this excludes repeat on incorrect trials)
        performance = block.events.sessionPerformanceValues(:,end-10:end);
        
        % Get whether all contrasts were used
        use_all_contrasts = all(block.events.useContrastsValues(end-10:end));
        
        % Store in behavior structure
        bhv.session_duration(curr_day) = session_duration;
        bhv.n_trials(curr_day) = n_trials;
        bhv.total_water(curr_day) = total_water;
        bhv.wheel_movement(curr_day) = sum(abs(wheel_movement));
        bhv.wheel_bias(curr_day) = wheel_bias;
        bhv.conditions = performance(1,:);
        bhv.n_trials_condition(curr_day,:) = performance(2,:);
        bhv.go_left_trials(curr_day,:) = performance(end,:);
        bhv.use_all_contrasts(curr_day) = use_all_contrasts;
        
        disp(['Behavior: day ' num2str(curr_day) '/' num2str(length(experiments))]);
    end
    
end

% Plot summary
day_num = cellfun(@(x) datenum(x),{experiments.day});
day_labels = cellfun(@(x) x(6:end),{experiments.day},'uni',false);
figure('Name',animal)

% Trials and water
subplot(2,2,1);
yyaxis left
plot(day_num,bhv.n_trials./bhv.session_duration,'linewidth',2);
ylabel('Trials/min');
yyaxis right
plot(day_num,bhv.total_water,'linewidth',2);
ylabel('Total water');
xlabel('Session');
set(gca,'XTick',day_num);
set(gca,'XTickLabel',day_labels);
set(gca,'XTickLabelRotation',90);

imaging_days = day_num([experiments.imaging]);
for i = 1:length(imaging_days)
    line(repmat(imaging_days(i),1,2),ylim,'color','b');
end

ephys_days = day_num([experiments.ephys]);
for i = 1:length(ephys_days)
    line(repmat(ephys_days(i),1,2),ylim,'color','r','linestyle','--');
end

% Wheel movement and bias
subplot(2,2,2);
yyaxis left
plot(day_num,bhv.wheel_movement./bhv.session_duration,'linewidth',2);
ylabel('Wheel movement / min');
yyaxis right
plot(day_num,bhv.wheel_bias,'linewidth',2);
ylim([-1,1]);
ylabel('Wheel bias');
xlabel('Session');
set(gca,'XTick',day_num);
set(gca,'XTickLabel',day_labels);
set(gca,'XTickLabelRotation',90);

imaging_days = day_num([experiments.imaging]);
for i = 1:length(imaging_days)
    line(repmat(imaging_days(i),1,2),ylim,'color','b');
end

ephys_days = day_num([experiments.ephys]);
for i = 1:length(ephys_days)
    line(repmat(ephys_days(i),1,2),ylim,'color','r','linestyle','--');
end

% Psychometric over all days
subplot(2,2,3);
imagesc(bhv.conditions,1:size(bhv.go_left_trials),bhv.go_left_trials./bhv.n_trials_condition);colormap(redblue)
color_step = diff(caxis)/size(colormap,1);
colormap([repmat(0.5,3);colormap]);
caxis([0-color_step*(255/size(colormap,1)),1]);
c = colorbar;
ylabel(c,'Go left (frac)');
xlabel('Condition');
ylabel('Session');
set(gca,'YTick',1:length(experiments));
set(gca,'YTickLabel',day_labels);
axis square;
hold on;
if any([experiments.imaging])
    plot(0,find([experiments.imaging]),'.k');
end
if any([experiments.ephys])
    plot(0,find([experiments.ephys]),'ok');
end

% Psychometric of combined days that use all contrasts
subplot(2,2,4); hold on;
combine_days = bhv.use_all_contrasts;
combined_psychometric = sum(bhv.go_left_trials(combine_days,:))./sum(bhv.n_trials_condition(combine_days,:),1);
plot(bhv.conditions,combined_psychometric,'color','k','linewidth',2);
combine_days_performance = bhv.go_left_trials(combine_days,:)./bhv.n_trials_condition(combine_days,:);
errorbar(bhv.conditions,nanmean(combine_days_performance,1), ...
    nanstd(combine_days_performance,[],1)./sqrt(sum(~isnan(combine_days_performance))),'r');
xlim([-1,1]);
ylim([0,1]);
axis square;
xlabel('Condition');
ylabel('Fraction go left (combined)');
legend({'Sum','Average'});

%% Batch behavior

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

bhv = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = false;
    load_parts.ephys = false;   
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment              
        
        % Time of session (in minutes)
        session_duration = block.duration/60;
        
        % Water
        total_water = sum(block.outputs.rewardValues);
 
        % Performance (note that this excludes repeat on incorrect trials)
        performance = block.events.sessionPerformanceValues(:,end-10:end);
        
        % Trial conditions and behavior times (pad to the same number)
        n_trials = length(block.paramsValues);
        trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
        stim_to_move = wheel_move_time - stimOn_times';
        stim_to_feedback = padarray(signals_events.responseTimes, ...
            [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
            padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');    
        
        % Define trials to use 
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % L/R choice
        go_right = (signals_events.trialSideValues(1:n_trials) == 1 & trial_outcome(1:n_trials) == -1) | ...
            (signals_events.trialSideValues(1:n_trials) == -1 & trial_outcome(1:n_trials) == 1);
        go_left = (signals_events.trialSideValues(1:n_trials) == -1 & trial_outcome(1:n_trials) == -1) | ...
            (signals_events.trialSideValues(1:n_trials) == 1 & trial_outcome(1:n_trials) == 1);
        
        trial_choice = go_right - go_left;
        
        if any(~ismember(trial_choice(use_trials),[-1,1]))
            error('Non-valid choice')
        end
               
        % Max wheel speed around choice
        % CHANGE THIS: make it be stim -> feedback for each trial
        use_move_times = reshape(wheel_move_time(use_trials),[],1);       
        interval_surround = [0,0.6];
        sample_t = 0.01;
        t_surround = interval_surround(1):sample_t:interval_surround(2);
        t_peri_move = bsxfun(@plus,use_move_times,t_surround);       
        move_aligned_wheel = interp1(conv2(Timeline.rawDAQTimestamps,[0.5,0.5],'valid'), ...
            wheel_velocity,t_peri_move);
        min_wheel_velocity = min(move_aligned_wheel,[],2);
        max_wheel_velocity = max(move_aligned_wheel,[],2);
        trial_wheel_velocity = nan(sum(use_trials),1);
        trial_wheel_velocity(trial_choice(use_trials) == -1) = min_wheel_velocity(trial_choice(use_trials) == -1);
        trial_wheel_velocity(trial_choice(use_trials) == 1) = max_wheel_velocity(trial_choice(use_trials) == 1);
        
        % Store in behavior structure
        bhv(curr_animal).session_duration(curr_day) = session_duration;
        bhv(curr_animal).n_trials(curr_day) = n_trials;
        bhv(curr_animal).total_water(curr_day) = total_water;
        bhv(curr_animal).conditions = performance(1,:);
                        
        bhv(curr_animal).trial_side{curr_day} = signals_events.trialSideValues(use_trials);
        bhv(curr_animal).trial_contrast{curr_day} = signals_events.trialContrastValues(use_trials);
        bhv(curr_animal).trial_choice{curr_day} = trial_choice(use_trials);
        bhv(curr_animal).trial_outcome{curr_day} = trial_outcome(use_trials);
        bhv(curr_animal).stim_to_move{curr_day} = stim_to_move(use_trials);
        bhv(curr_animal).stim_to_feedback{curr_day} = stim_to_feedback(use_trials);
        bhv(curr_animal).trial_wheel_velocity{curr_day} = trial_wheel_velocity;
        
        clearvars -except animals protocol experiments load_parts curr_animal curr_day animal bhv 
    end
    
    AP_print_progress_fraction(curr_animal,length(animals));
    
end

conditions = unique(vertcat(bhv.conditions),'rows');
trial_choice_cat = arrayfun(@(x) horzcat(bhv(x).trial_choice{:}),1:length(bhv),'uni',false);
trial_outcome_cat = arrayfun(@(x) horzcat(bhv(x).trial_outcome{:}),1:length(bhv),'uni',false);
trial_side_cat = arrayfun(@(x) horzcat(bhv(x).trial_side{:}),1:length(bhv),'uni',false);
trial_contrast_cat = arrayfun(@(x) horzcat(bhv(x).trial_contrast{:}),1:length(bhv),'uni',false);
trial_condition_cat = cellfun(@(side,contrast) side.*contrast,trial_side_cat,trial_contrast_cat,'uni',false);
trial_wheel_velocity_cat = arrayfun(@(x) vertcat(bhv(x).trial_wheel_velocity{:})',1:length(bhv),'uni',false);
stim_to_move_cat = arrayfun(@(x) horzcat(bhv(x).stim_to_move{:}),1:length(bhv),'uni',false);
stim_to_feedback_cat = arrayfun(@(x) horzcat(bhv(x).stim_to_feedback{:}),1:length(bhv),'uni',false);

% Distinguish early/late movements
go_time = 0.5;
trial_timing = arrayfun(@(animal) cellfun(@(x) 1+(x > go_time), ...
    bhv(animal).stim_to_move,'uni',false),1:length(bhv),'uni',false);
trial_timing_cat = arrayfun(@(animal) ...
    horzcat(trial_timing{animal}{:}),1:length(bhv),'uni',false);

% Get easy performance / anti-lapse rate
performance_contrasts = [0.5,1];

l_correct = cellfun(@(contrast,side,choice) cellfun(@(contrast,side,choice) ...
    nanmean(choice(ismember(contrast,performance_contrasts) & side == -1) == 1), ...
    contrast,side,choice), ...
    {bhv.trial_contrast},{bhv.trial_side},{bhv.trial_choice},'uni',false);

r_correct = cellfun(@(contrast,side,choice) cellfun(@(contrast,side,choice) ...
    nanmean(choice(ismember(contrast,performance_contrasts) & side == 1) == -1), ...
    contrast,side,choice), ...
    {bhv.trial_contrast},{bhv.trial_side},{bhv.trial_choice},'uni',false);

thresh_correct = 0.85; % threshold for lapse rate 
use_experiment = cellfun(@(l,r) min(l,r) > thresh_correct,l_correct,r_correct,'uni',false);

% Plot psychometric 
frac_left = cell2mat(cellfun(@(choice,condition) ...
    grpstats(choice == -1,condition),trial_choice_cat,trial_condition_cat,'uni',false));

frac_left_earlymove = cell2mat(cellfun(@(choice,condition,stim_to_move) ...
    grpstats(choice(stim_to_move < 0.5) == -1,condition(stim_to_move < 0.5)), ...
    trial_choice_cat,trial_condition_cat,stim_to_move_cat,'uni',false));

frac_left_latemove = cell2mat(cellfun(@(choice,condition,stim_to_move) ...
    grpstats(choice(stim_to_move >= 0.5) == -1,condition(stim_to_move >= 0.5)), ...
    trial_choice_cat,trial_condition_cat,stim_to_move_cat,'uni',false));

figure;

subplot(1,3,1); hold on;
plot(conditions,frac_left,'linewidth',2,'color',[0.5,0.5,0.5]);
plot(conditions,nanmean(frac_left,2),'linewidth',5,'color','k');
xlim([-1,1]);
ylim([0,1]);
line([0,0],ylim,'linestyle','--','color','k');
line(xlim,[0.5,0.5],'linestyle','--','color','k');
xlabel('Condition');
ylabel('Fraction go left');
title('All trials');

subplot(1,3,2); hold on;
plot(conditions,frac_left_earlymove,'linewidth',2,'color',[0.5,0.5,0.5]);
plot(conditions,nanmean(frac_left_earlymove,2),'linewidth',5,'color','k');
xlim([-1,1]);
ylim([0,1]);
line([0,0],ylim,'linestyle','--','color','k');
line(xlim,[0.5,0.5],'linestyle','--','color','k');
xlabel('Condition');
ylabel('Fraction go left');
title('Early move');

subplot(1,3,3); hold on;
plot(conditions,frac_left_latemove,'linewidth',2,'color',[0.5,0.5,0.5]);
plot(conditions,nanmean(frac_left_latemove,2),'linewidth',5,'color','k');
xlim([-1,1]);
ylim([0,1]);
line([0,0],ylim,'linestyle','--','color','k');
line(xlim,[0.5,0.5],'linestyle','--','color','k');
xlabel('Condition');
ylabel('Fraction go left');
title('Late move');

% Plot wheel velocity by condition
[velocity_condition_hit_earlymove_grp,group] = cellfun(@(wheel,outcome,condition,stim_to_move) ...
    grpstats(wheel(outcome == 1 & stim_to_move < 0.5),condition(outcome == 1 & stim_to_move < 0.5), ...
    {'mean','gname'}),trial_wheel_velocity_cat,trial_outcome_cat,trial_condition_cat,stim_to_move_cat,'uni',false);
group = cellfun(@(x) cellfun(@(x) str2num(x),x),group,'uni',false);
all_group = unique(vertcat(group{:}),'rows');
velocity_condition_hit_earlymove = nan(size(all_group,1),length(bhv));
for curr_animal = 1:length(bhv)
    [~,group_idx] = intersect(all_group,group{curr_animal},'rows');
    velocity_condition_hit_earlymove(group_idx,curr_animal) = velocity_condition_hit_earlymove_grp{curr_animal};
end

[velocity_condition_miss_earlymove_grp,group] = cellfun(@(wheel,outcome,condition,stim_to_move) ...
    grpstats(wheel(outcome == -1 & stim_to_move < 0.5),condition(outcome == -1 & stim_to_move < 0.5), ...
    {'mean','gname'}),trial_wheel_velocity_cat,trial_outcome_cat,trial_condition_cat,stim_to_move_cat,'uni',false);
group = cellfun(@(x) cellfun(@(x) str2num(x),x),group,'uni',false);
all_group = unique(vertcat(group{:}),'rows');
velocity_condition_miss_earlymove = nan(size(all_group,1),length(bhv));
for curr_animal = 1:length(bhv)
    [~,group_idx] = intersect(all_group,group{curr_animal},'rows');
    velocity_condition_miss_earlymove(group_idx,curr_animal) = velocity_condition_miss_earlymove_grp{curr_animal};
end

figure; hold on;
plot(conditions,abs(velocity_condition_miss_earlymove),'color',[0.8,0.6,0.6],'linewidth',2);
plot(conditions,abs(velocity_condition_hit_earlymove),'color',[0.6,0.8,0.6],'linewidth',2);
plot(conditions,nanmean(abs(velocity_condition_miss_earlymove),2),'color',[0.8,0,0],'linewidth',2);
plot(conditions,nanmean(abs(velocity_condition_hit_earlymove),2),'color',[0,0.8,0],'linewidth',2);
ylabel('Max wheel velocity towards choice');
xlabel('Condition');

% Plot distribution of stim to move and feedback times
move_time_bins = 0:0.01:1;
day_split = 4;

move_time_centers = move_time_bins(1:end-1) + diff(move_time_bins)/2;
stim_to_move_binned = zeros(day_split,length(move_time_bins)-1,length(bhv));

for curr_animal = 1:length(bhv)
    for curr_day = 1:length(bhv(curr_animal).stim_to_move)
        curr_data = bhv(curr_animal).stim_to_move{curr_day};
        trials_split = round(linspace(1,length(curr_data),day_split+1));
        for curr_split = 1:day_split
            stim_to_move_binned(curr_split,:,curr_animal) = ...
                stim_to_move_binned(curr_split,:,curr_animal) + ...
                histcounts(curr_data(trials_split(curr_split): ...
                trials_split(curr_split+1)),move_time_bins);
        end
    end
end
stim_to_move_binned_norm = bsxfun(@rdivide,stim_to_move_binned, ...
    sum(sum(stim_to_move_binned,1),2));

figure; hold on;
for curr_animal = 1:length(bhv)
    AP_stackplot(stim_to_move_binned_norm(:,:,curr_animal)', ...
        move_time_centers,0.02,false,[0.5,0.5,0.5],1:day_split);
end
AP_stackplot(nanmean(stim_to_move_binned_norm,3)', ...
    move_time_centers,0.02,false,'k',1:day_split);
xlabel('Time to movement onset')
ylabel('Frequency by fraction within day')
axis tight
line([0.5,0.5],ylim,'linestyle','--','color','k');
   
% Plot stim to move / feedback by condition and success pooling days
% (separate early/late movements)

[group,stim_to_move_prebeep] = arrayfun(@(x) ...
    grpstats(stim_to_move_cat{x}(trial_timing_cat{x} == 1)', ...
    [trial_condition_cat{x}(trial_timing_cat{x} == 1)', ...
    trial_outcome_cat{x}(trial_timing_cat{x} == 1)'], ...
    {'gname','nanmedian'}),1:length(bhv),'uni',false);
group = cellfun(@(x) cellfun(@(x) str2num(x),x),group,'uni',false);
all_group = unique(vertcat(group{:}),'rows');
stim_to_move_prebeep_cat = nan(size(all_group,1),length(bhv));
for curr_animal = 1:length(bhv)
    [~,group_idx] = intersect(all_group,group{curr_animal},'rows');
    stim_to_move_prebeep_cat(group_idx,curr_animal) = stim_to_move_prebeep{curr_animal};
end

[group,stim_to_move_postbeep] = arrayfun(@(x) ...
    grpstats(stim_to_move_cat{x}(trial_timing_cat{x} == 2)', ...
    [trial_condition_cat{x}(trial_timing_cat{x} == 2)', ...
    trial_outcome_cat{x}(trial_timing_cat{x} == 2)'], ...
    {'gname','nanmedian'}),1:length(bhv),'uni',false);
group = cellfun(@(x) cellfun(@(x) str2num(x),x),group,'uni',false);
all_group = unique(vertcat(group{:}),'rows');
stim_to_move_postbeep_cat = nan(size(all_group,1),length(bhv));
for curr_animal = 1:length(bhv)
    [~,group_idx] = intersect(all_group,group{curr_animal},'rows');
    stim_to_move_postbeep_cat(group_idx,curr_animal) = stim_to_move_postbeep{curr_animal};
end

[group,stim_to_feedback_prebeep] = arrayfun(@(x) ...
    grpstats(stim_to_feedback_cat{x}(trial_timing_cat{x} == 1)', ...
    [trial_condition_cat{x}(trial_timing_cat{x} == 1)', ...
    trial_outcome_cat{x}(trial_timing_cat{x} == 1)'], ...
    {'gname','nanmedian'}),1:length(bhv),'uni',false);
group = cellfun(@(x) cellfun(@(x) str2num(x),x),group,'uni',false);
all_group = unique(vertcat(group{:}),'rows');
stim_to_feedback_prebeep_cat = nan(size(all_group,1),length(bhv));
for curr_animal = 1:length(bhv)
    [~,group_idx] = intersect(all_group,group{curr_animal},'rows');
    stim_to_feedback_prebeep_cat(group_idx,curr_animal) = stim_to_feedback_prebeep{curr_animal};
end

[group,stim_to_feedback_postbeep] = arrayfun(@(x) ...
    grpstats(stim_to_feedback_cat{x}(trial_timing_cat{x} == 2)', ...
    [trial_condition_cat{x}(trial_timing_cat{x} == 2)', ...
    trial_outcome_cat{x}(trial_timing_cat{x} == 2)'], ...
    {'gname','nanmedian'}),1:length(bhv),'uni',false);
group = cellfun(@(x) cellfun(@(x) str2num(x),x),group,'uni',false);
all_group = unique(vertcat(group{:}),'rows');
stim_to_feedback_postbeep_cat = nan(size(all_group,1),length(bhv));
for curr_animal = 1:length(bhv)
    [~,group_idx] = intersect(all_group,group{curr_animal},'rows');
    stim_to_feedback_postbeep_cat(group_idx,curr_animal) = stim_to_feedback_postbeep{curr_animal};
end

figure;
subplot(2,2,1); hold on;
plot(all_group(all_group(:,2) == -1,1),stim_to_move_prebeep_cat(all_group(:,2) == -1,:),'color',[0.8,0.6,0.6],'linewidth',2);
plot(all_group(all_group(:,2) == 1,1),stim_to_move_prebeep_cat(all_group(:,2) == 1,:),'color',[0.6,0.8,0.6],'linewidth',2);
p1 = plot(all_group(all_group(:,2) == -1,1),nanmean(stim_to_move_prebeep_cat(all_group(:,2) == -1,:),2),'color',[0.8,0,0],'linewidth',5);
p2 = plot(all_group(all_group(:,2) == 1,1),nanmean(stim_to_move_prebeep_cat(all_group(:,2) == 1,:),2),'color',[0,0.8,0],'linewidth',5);
ylabel('Stim to move time (s)');
xlabel('Condition');
title('Pre-beep movement')
legend([p1,p2],{'Miss','Hit'})
axis tight
ylim([0,1])
line(xlim,[0.5,0.5],'linestyle','--','color','k')

subplot(2,2,2); hold on;
plot(all_group(all_group(:,2) == -1,1),stim_to_feedback_prebeep_cat(all_group(:,2) == -1,:),'color',[0.8,0.6,0.6],'linewidth',2);
plot(all_group(all_group(:,2) == 1,1),stim_to_feedback_prebeep_cat(all_group(:,2) == 1,:),'color',[0.6,0.8,0.6],'linewidth',2);
p1 = plot(all_group(all_group(:,2) == -1,1),nanmean(stim_to_feedback_prebeep_cat(all_group(:,2) == -1,:),2),'color',[0.8,0,0],'linewidth',5);
p2 = plot(all_group(all_group(:,2) == 1,1),nanmean(stim_to_feedback_prebeep_cat(all_group(:,2) == 1,:),2),'color',[0,0.8,0],'linewidth',5);
ylabel('Stim to feedback time (s)');
xlabel('Condition');
title('Pre-beep movement')
legend([p1,p2],{'Miss','Hit'})
axis tight
ylim([0,1])
line(xlim,[0.5,0.5],'linestyle','--','color','k')

subplot(2,2,3); hold on;
plot(all_group(all_group(:,2) == -1,1),stim_to_move_postbeep_cat(all_group(:,2) == -1,:),'color',[0.8,0.6,0.6],'linewidth',2);
plot(all_group(all_group(:,2) == 1,1),stim_to_move_postbeep_cat(all_group(:,2) == 1,:),'color',[0.6,0.8,0.6],'linewidth',2);
p1 = plot(all_group(all_group(:,2) == -1,1),nanmean(stim_to_move_postbeep_cat(all_group(:,2) == -1,:),2),'color',[0.8,0,0],'linewidth',5);
p2 = plot(all_group(all_group(:,2) == 1,1),nanmean(stim_to_move_postbeep_cat(all_group(:,2) == 1,:),2),'color',[0,0.8,0],'linewidth',5);
ylabel('Stim to move time (s)');
xlabel('Condition');
title('Post-beep movement')
legend([p1,p2],{'Miss','Hit'})
axis tight
ylim([0,1])
line(xlim,[0.5,0.5],'linestyle','--','color','k')

subplot(2,2,4); hold on;
plot(all_group(all_group(:,2) == -1,1),stim_to_feedback_postbeep_cat(all_group(:,2) == -1,:),'color',[0.8,0.6,0.6],'linewidth',2);
plot(all_group(all_group(:,2) == 1,1),stim_to_feedback_postbeep_cat(all_group(:,2) == 1,:),'color',[0.6,0.8,0.6],'linewidth',2);
p1 = plot(all_group(all_group(:,2) == -1,1),nanmean(stim_to_feedback_postbeep_cat(all_group(:,2) == -1,:),2),'color',[0.8,0,0],'linewidth',5);
p2 = plot(all_group(all_group(:,2) == 1,1),nanmean(stim_to_feedback_postbeep_cat(all_group(:,2) == 1,:),2),'color',[0,0.8,0],'linewidth',5);
ylabel('Stim to feedback time (s)');
xlabel('Condition');
title('Post-beep movement')
legend([p1,p2],{'Miss','Hit'})
axis tight
ylim([0,1])
line(xlim,[0.5,0.5],'linestyle','--','color','k')

% Get numbers of trials each day/animal (contrast/side/choice)
% [animal,contrast,side,choice,timing]
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];

% Re-number all conditions
conditions = combvec(contrasts,sides,choices,timings);
conditions_str = cellfun(@(x) sprintf('%+g/',x),mat2cell(conditions,size(conditions,1),ones(1,size(conditions,2))),'uni',false);
n_conditions = size(conditions,2);

condition_count = arrayfun(@(x) nan(n_conditions,length(bhv(x).session_duration)),1:length(bhv),'uni',false);
for curr_animal = 1:length(bhv);
    for curr_day = 1:length(bhv(curr_animal).session_duration)
        trial_conditions = ...
            [bhv(curr_animal).trial_contrast{curr_day}; bhv(curr_animal).trial_side{curr_day}; ...
            bhv(curr_animal).trial_choice{curr_day}; trial_timing{curr_animal}{curr_day}];
        [~,trial_id] = ismember(trial_conditions',conditions','rows');
        curr_counts = histcounts(trial_id,'BinLimits',[1,n_conditions],'BinMethod','integers');
        condition_count{curr_animal}(:,curr_day) = curr_counts;
    end
end

figure; colormap(gray);
max_count = max(reshape(horzcat(condition_count{:}),[],1));
for curr_animal = 1:length(bhv)
    subplot(1,length(bhv)+1,curr_animal);
    imagesc(condition_count{curr_animal});
    caxis([0,max_count]);
    title(animals{curr_animal});
    if curr_animal == 1
        set(gca,'YTick',1:n_conditions,'YTickLabel',conditions_str)
        ylabel(['contrast/side/choice/timing, max ' num2str(max_count)])
        xlabel('day');
    else
        axis off
    end
end
subplot(1,length(bhv)+1,length(bhv)+1);
condition_count_cat = cell2mat(cellfun(@(x) sum(x,2),condition_count,'uni',false));
plot(condition_count_cat,1:n_conditions);
axis tight;
set(gca,'YDir','reverse','YTick',1:n_conditions,'YTickLabel',conditions_str,'YAxisLocation','Right')
xlabel('# Trials')
legend(animals);

    
    
    
    
    
    

