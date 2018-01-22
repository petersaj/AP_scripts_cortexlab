%% Single behavior

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
        
        % Trial counts (subtract 1 because last trial always incomplete)
        n_trials = length(block.paramsValues) - 1;
        total_water = sum(block.outputs.rewardValues);
 
        % Performance (note that this excludes repeat on incorrect trials)
        performance = block.events.sessionPerformanceValues(:,end-10:end);
        
        % Trial conditions and behavior times
        % (only use trials that were completed and not repeat trials)
        use_trials = any([signals_events.hitValues;signals_events.missValues],1) & ~signals_events.repeatTrialValues;
        
        trial_conditions = signals_events.trialSideValues(use_trials).*signals_events.trialContrastValues(use_trials);
        trial_hit = signals_events.hitValues(use_trials);
        stim_to_move = wheel_move_time(use_trials) - stimOn_times(use_trials)';
        stim_to_feedback = signals_events.responseTimes(use_trials) - stimOn_times(use_trials)';        
        
        % Store in behavior structure
        bhv(curr_animal).session_duration(curr_day) = session_duration;
        bhv(curr_animal).n_trials(curr_day) = n_trials;
        bhv(curr_animal).total_water(curr_day) = total_water;
        bhv(curr_animal).conditions = performance(1,:);
        bhv(curr_animal).n_trials_condition(curr_day,:) = performance(2,:);
        bhv(curr_animal).go_left_trials(curr_day,:) = performance(end,:);
        
        bhv(curr_animal).trial_conditions{curr_day} = trial_conditions;
        bhv(curr_animal).trial_hit{curr_day} = trial_hit;
        bhv(curr_animal).stim_to_move{curr_day} = stim_to_move;
        bhv(curr_animal).stim_to_feedback{curr_day} = stim_to_feedback;
        
        clearvars -except animals protocol experiments load_parts curr_animal curr_day animal bhv 
    end
    
    AP_print_progress_fraction(curr_animal,length(animals));
    
end

% Plot psychometric pooling days
conditions = unique(vertcat(bhv.conditions),'rows');
pooled_psychometric = cell2mat(cellfun(@(go_left,n_trials) sum(go_left,1)./sum(n_trials,1), ...
    {bhv.go_left_trials},{bhv.n_trials_condition},'uni',false)');

figure; hold on;
plot(conditions,pooled_psychometric','linewidth',2,'color',[0.5,0.5,0.5]);
plot(conditions,nanmean(pooled_psychometric,1),'linewidth',5','color','k');

xlim([-1,1]);
ylim([0,1]);
line([0,0],ylim,'linestyle','--','color','k');
line(xlim,[0.5,0.5],'linestyle','--','color','k');
axis square;
xlabel('Condition');
ylabel('Fraction go left');
title('Pooling across days');

% Plot stim to move / feedback by condition and success pooling days
[group,stim_to_move] = arrayfun(@(x) ...
    grpstats(horzcat(bhv(x).stim_to_move{:})', ...
    [horzcat(bhv(x).trial_conditions{:})',horzcat(bhv(x).trial_hit{:})'], ...
    {'gname','nanmedian'}),1:length(bhv),'uni',false);
group = cellfun(@(x) cellfun(@(x) str2num(x),x),group,'uni',false);
if ~any(reshape(bsxfun(@eq,cat(3,group{:}),group{1}),[],1))
    error('Different conditions across animals')
end
group = group{1};
stim_to_move_cat = horzcat(stim_to_move{:});

[group,stim_to_feedback] = arrayfun(@(x) ...
    grpstats(horzcat(bhv(x).stim_to_feedback{:})', ...
    [horzcat(bhv(x).trial_conditions{:})',horzcat(bhv(x).trial_hit{:})'], ...
    {'gname','nanmedian'}),1:length(bhv),'uni',false);
group = cellfun(@(x) cellfun(@(x) str2num(x),x),group,'uni',false);
if ~any(reshape(bsxfun(@eq,cat(3,group{:}),group{1}),[],1))
    error('Different conditions across animals')
end
group = group{1};
stim_to_feedback_cat = horzcat(stim_to_feedback{:});

figure;
subplot(1,2,1); hold on;
plot(group(group(:,2) == 1,1),stim_to_move_cat(group(:,2) == 0,:),'color',[0.8,0.6,0.6],'linewidth',2);
plot(group(group(:,2) == 1,1),stim_to_move_cat(group(:,2) == 1,:),'color',[0.6,0.8,0.6],'linewidth',2);
p1 = plot(group(group(:,2) == 1,1),nanmean(stim_to_move_cat(group(:,2) == 0,:),2),'color',[0.8,0,0],'linewidth',5);
p2 = plot(group(group(:,2) == 1,1),nanmean(stim_to_move_cat(group(:,2) == 1,:),2),'color',[0,0.8,0],'linewidth',5);
ylabel('Stim to move time (s)');
xlabel('Condition');
legend([p1,p2],{'Miss','Hit'})
axis tight square
line(xlim,[0.5,0.5],'linestyle','--','color','k')

subplot(1,2,2); hold on;
plot(group(group(:,2) == 1,1),stim_to_feedback_cat(group(:,2) == 0,:),'color',[0.8,0.6,0.6],'linewidth',2);
plot(group(group(:,2) == 1,1),stim_to_feedback_cat(group(:,2) == 1,:),'color',[0.6,0.8,0.6],'linewidth',2);
p1 = plot(group(group(:,2) == 1,1),nanmean(stim_to_feedback_cat(group(:,2) == 0,:),2),'color',[0.8,0,0],'linewidth',5);
p2 = plot(group(group(:,2) == 1,1),nanmean(stim_to_feedback_cat(group(:,2) == 1,:),2),'color',[0,0.8,0],'linewidth',5);
ylabel('Stim to feedback time (s)');
xlabel('Condition');
legend([p1,p2],{'Miss','Hit'})
axis tight square
line(xlim,[0.5,0.5],'linestyle','--','color','k')






    
    
    
    
    
    
    
    
    
    
    

