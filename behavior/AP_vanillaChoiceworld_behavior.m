% Load in and plot behavior

animal = 'AP029';
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












