%% Get and plot behavior within mice (vanillaChoiceworld)
% animals = {'AP063','AP064','AP066','AP068','AP071','AP085','AP086','AP087'};
% animals = {'AP047','AP048','AP077','AP079'};
animals = {'AP089','AP090','AP091'};
protocol = 'vanillaChoiceworldBiasNoCue';
flexible_name = false;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    bhv = struct;
            
    if isempty(experiments)
        disp(['No behavior data: ' animal]);
        continue
    end
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment_num = experiments(curr_day).experiment;
         
        % If multiple experiments, only use the last one 
        % (usually multiple happens if mess ups and final one is good)
        for curr_experiment = length(experiment_num)
            
            experiment = experiment_num(curr_experiment);
            
            [block_filename, block_exists] = AP_cortexlab_filename(animal,day,experiment,'block');
            
            % Load the block file
            load(block_filename)
            
            % Get protocol name
            [~,curr_protocol] = fileparts(block.expDef);
            
            % Time of session (in minutes)
            session_duration = block.duration/60;
            
            % Trial counts
            n_trials = length(block.paramsValues);
            total_water = sum(block.outputs.rewardValues);
            
            % Estimate reaction time
            % (evenly resample velocity - not even natively);
            wheel_resample_rate = 1000;
            wheel_t_resample = block.inputs.wheelTimes(1):1/wheel_resample_rate:block.inputs.wheelTimes(end);
            wheel_values_resample = interp1(block.inputs.wheelTimes,block.inputs.wheelValues,wheel_t_resample);
            
            wheel_smooth_t = 0.05; % seconds
            wheel_smooth_samples = round(wheel_smooth_t*wheel_resample_rate);
            wheel_velocity = interp1(conv(wheel_t_resample,[1,1]/2,'valid'), ...
                diff(smooth(wheel_values_resample,wheel_smooth_samples)),wheel_t_resample)';
            
            wheel_thresh = 0.025;
            wheel_starts = wheel_t_resample(abs(wheel_velocity(1:end-1)) < wheel_thresh & ...
                abs(wheel_velocity(2:end)) > wheel_thresh);
            
            response_trials = 1:length(block.events.responseValues);
            trial_wheel_starts = arrayfun(@(x) ...
                wheel_starts(find(wheel_starts > block.events.stimOnTimes(x),1)), ...
                response_trials);
            trial_move_t = trial_wheel_starts - block.events.stimOnTimes(response_trials);
            
            % Wheel movements/biases
            left_wheel_velocity = abs(wheel_velocity.*(wheel_velocity < 0));
            right_wheel_velocity = abs(wheel_velocity.*(wheel_velocity > 0));
            wheel_bias = (nansum(right_wheel_velocity)-nansum(left_wheel_velocity))/ ...
                (nansum(right_wheel_velocity)+nansum(left_wheel_velocity));
            
            % Get reaction times for each stimulus
            trial_stim = block.events.trialContrastValues(response_trials).*block.events.trialSideValues(response_trials);
            stim_list = unique(reshape(unique(block.events.contrastsValues).*[-1;1],[],1));
            [~,trial_stim_idx] = ismember(trial_stim,stim_list);
            stim_rxn_time = accumarray(trial_stim_idx(response_trials)',trial_move_t',[11,1],@nanmedian,nan);
            
            % Performance (note that this excludes repeat on incorrect trials)
            performance = block.events.sessionPerformanceValues(:,end-10:end);
            
            % Hit/miss recorded for last trial, circshift to align
            hitValues = circshift(block.events.hitValues,[0,-1]);
            missValues = circshift(block.events.missValues,[0,-1]);
            frac_correct = nanmean(hitValues(trial_stim ~= 0));   
            
            % Dprime (loglinear approach to normalizing)
%             dprime = norminv(nanmean(hitValues(trial_stim > 0))) - ...
%                 norminv(nanmean(missValues(trial_stim < 0)));
            dprime = ...
                norminv((sum(hitValues(trial_stim > 0))+0.5)./(sum(trial_stim > 0) + 1)) - ...
                norminv((sum(missValues(trial_stim < 0))+0.5)./(sum(trial_stim < 0) + 1));

            % Get whether all contrasts were used
            use_all_contrasts = all(block.events.useContrastsValues(end-10:end));
            
            % Store in behavior structure
            bhv.protocol{curr_day} = curr_protocol;
            bhv.session_duration(curr_day) = session_duration;
            bhv.n_trials(curr_day) = n_trials;
            bhv.total_water(curr_day) = total_water;
            bhv.wheel_velocity(curr_day) = nansum(abs(wheel_velocity));
            bhv.wheel_bias(curr_day) = wheel_bias;
            bhv.conditions = performance(1,:);
            bhv.n_trials_condition(curr_day,:) = performance(2,:);
            bhv.go_left_trials(curr_day,:) = performance(end,:);
            bhv.use_all_contrasts(curr_day) = use_all_contrasts;
            bhv.stim_rxn_time(curr_day,:) = stim_rxn_time;
            
            bhv.frac_correct(curr_day) = frac_correct;
            bhv.dprime(curr_day) = dprime;           
            
            AP_print_progress_fraction(curr_day,length(experiments));
        end
        
    end
    
    % Plot summary
    day_num = cellfun(@(x) datenum(x),{experiments.day});
    day_labels = cellfun(@(day,protocol) [day(6:end)], ...
        {experiments.day},bhv.protocol,'uni',false);
    
    [unique_protocols,~,protocol_idx] = unique(bhv.protocol);
    protocol_col = hsv(length(unique_protocols));
    
    figure('Name',animal)
    
    % Trials and water
    subplot(3,3,1); hold on;
    yyaxis left
    % plot(day_num,bhv.n_trials./bhv.session_duration,'linewidth',2);
    % ylabel('Trials/min');
    plot(day_num,bhv.n_trials,'linewidth',2);
    ylabel('Trials');
    yyaxis right
    plot(day_num,bhv.total_water,'linewidth',2);
    ylabel('Total water');
    xlabel('Session');
    set(gca,'XTick',day_num);
    set(gca,'XTickLabel',day_labels);
    set(gca,'XTickLabelRotation',90);
    
    protocol_plot = gscatter(day_num,zeros(size(day_num)),[bhv.protocol]');
    
    imaging_days = day_num([experiments.imaging]);
    for i = 1:length(imaging_days)
        line(repmat(imaging_days(i),1,2),ylim,'color','k');
    end
    
    ephys_days = day_num([experiments.ephys]);
    for i = 1:length(ephys_days)
        line(repmat(ephys_days(i),1,2),ylim,'color','r','linestyle','--');
    end
    
    % Frac correct and d'
    subplot(3,3,2); hold on;
    yyaxis left
    % plot(day_num,bhv.n_trials./bhv.session_duration,'linewidth',2);
    % ylabel('Trials/min');
    plot(day_num,bhv.frac_correct,'linewidth',2);
    ylabel('Frac correct');
    yyaxis right
    plot(day_num,bhv.dprime,'linewidth',2);
    ylabel('d''');
    xlabel('Session');
    set(gca,'XTick',day_num);
    set(gca,'XTickLabel',day_labels);
    set(gca,'XTickLabelRotation',90);
        
    imaging_days = day_num([experiments.imaging]);
    for i = 1:length(imaging_days)
        line(repmat(imaging_days(i),1,2),ylim,'color','k');
    end
    
    ephys_days = day_num([experiments.ephys]);
    for i = 1:length(ephys_days)
        line(repmat(ephys_days(i),1,2),ylim,'color','r','linestyle','--');
    end
      
    % Wheel movement and bias
    subplot(3,3,3);
    yyaxis left
    % plot(day_num,bhv.wheel_velocity./bhv.session_duration,'linewidth',2);
    % ylabel('Wheel movement / min');
    plot(day_num,bhv.wheel_velocity,'linewidth',2);
    ylabel('Wheel movement');
    yyaxis right
    plot(day_num,bhv.wheel_bias,'linewidth',2);
    ylim([-1,1]);
    line(xlim,[0,0]);
    ylabel('Wheel bias');
    xlabel('Session');
    set(gca,'XTick',day_num);
    set(gca,'XTickLabel',day_labels);
    set(gca,'XTickLabelRotation',90);
    
    imaging_days = day_num([experiments.imaging]);
    for i = 1:length(imaging_days)
        line(repmat(imaging_days(i),1,2),ylim,'color','k');
    end
    
    ephys_days = day_num([experiments.ephys]);
    for i = 1:length(ephys_days)
        line(repmat(ephys_days(i),1,2),ylim,'color','r','linestyle','--');
    end
    
    % Psychometric over all days
    subplot(3,2,3);
    im = imagesc(bhv.conditions,1:size(bhv.go_left_trials),bhv.go_left_trials./bhv.n_trials_condition);
    set(im,'AlphaData',~isnan(get(im,'CData')));
    set(gca,'color',[0.5,0.5,0.5]);
    colormap(brewermap([],'*RdBu'));
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
    
    % Reaction time over days
    subplot(3,2,4);
    im = imagesc(bhv.conditions,1:size(bhv.go_left_trials),bhv.stim_rxn_time);
    set(im,'AlphaData',~isnan(get(im,'CData')));
    set(gca,'color',[0.5,0.5,0.5]);
    colormap(brewermap([],'*RdBu'));
    caxis([0.2,0.8])
    c = colorbar;
    ylabel(c,'Reaction time');
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
    subplot(3,2,5); hold on;
    combine_days = bhv.use_all_contrasts;
    combine_days_performance = bhv.go_left_trials(combine_days,:)./bhv.n_trials_condition(combine_days,:);
    errorbar(bhv.conditions,nanmean(combine_days_performance,1), ...
        nanstd(combine_days_performance,[],1)./sqrt(sum(~isnan(combine_days_performance))),'k','linewidth',2);
    xlim([-1,1]);
    ylim([0,1]);
    line(xlim,[0.5,0.5],'color','k','linestyle','--');
    line([0,0],ylim,'color','k','linestyle','--');
    axis square;
    xlabel('Condition');
    ylabel('Fraction go left');
    
    % Reaction time of combined days that use all contrasts
    subplot(3,2,6); hold on;
    combine_days = bhv.use_all_contrasts;
    combine_days_rxn_time = bhv.stim_rxn_time(combine_days,:);
    errorbar(bhv.conditions,nanmean(combine_days_rxn_time,1), ...
        nanstd(combine_days_rxn_time,[],1)./sqrt(sum(~isnan(combine_days_rxn_time))),'k','linewidth',2);
    xlim([-1,1]);
    line(xlim,[0.5,0.5],'color','k','linestyle','--');
    axis square;
    xlabel('Condition');
    ylabel('Reaction time');
    
    drawnow;
    
end

%% Get and plot single mouse behavior (milesChoiceworld/vanillaChoiceworld)
animals = {'AP060','AP061'};
protocol = 'choiceworld';
flexible_name = true;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    bhv = struct;
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment_num = experiments(curr_day).experiment;
        
        % If multiple experiments, only use the last one 
        % (usually multiple happens if mess ups and final one is good)
        for curr_experiment = length(experiment_num)
            
            experiment = experiment_num(curr_experiment);
            
            [block_filename, block_exists] = AP_cortexlab_filename(animal,day,experiment,'block');
            
            % Load the block file
            load(block_filename)
            
            % Get protocol name
            [~,curr_protocol] = fileparts(block.expDef);
            
            % Get relevant parameters (differently named)
            if contains(curr_protocol,'vanilla')
                stimOn_times = block.events.stimOnTimes;
                trial_contrast = block.events.trialContrastValues;
                trial_side = block.events.trialSideValues;
            else
                stimOn_times = block.events.stimulusOnTimes;
                trial_contrast = max([block.events.contrastRightValues; ...
                    block.events.contrastLeftValues],[],1);
                trial_side = ...
                    (block.events.contrastRightValues > 0) - ...
                    (block.events.contrastLeftValues > 0);
            end
            
            % Time of session (in minutes)
            session_duration = block.duration/60;
            
            % Trial counts
            n_trials = length(block.paramsValues);
            total_water = sum(block.outputs.rewardValues);
            
            % Estimate reaction time
            % (evenly resample velocity - not even natively);
            wheel_resample_rate = 1000;
            wheel_t_resample = block.inputs.wheelTimes(1):1/wheel_resample_rate:block.inputs.wheelTimes(end);
            wheel_values_resample = interp1(block.inputs.wheelTimes,block.inputs.wheelValues,wheel_t_resample);
            
            wheel_smooth_t = 0.05; % seconds
            wheel_smooth_samples = round(wheel_smooth_t*wheel_resample_rate);
            wheel_velocity = interp1(conv(wheel_t_resample,[1,1]/2,'valid'), ...
                diff(smooth(wheel_values_resample,wheel_smooth_samples)),wheel_t_resample)';
            
            wheel_thresh = 0.025;
            wheel_starts = wheel_t_resample(abs(wheel_velocity(1:end-1)) < wheel_thresh & ...
                abs(wheel_velocity(2:end)) > wheel_thresh);
            
            response_trials = 1:length(block.events.responseValues);
            trial_wheel_starts = arrayfun(@(x) ...
                wheel_starts(find(wheel_starts > stimOn_times(x),1)), ...
                response_trials);
            trial_move_t = trial_wheel_starts - stimOn_times(response_trials);
            
            % Wheel movements/biases
            left_wheel_velocity = abs(wheel_velocity.*(wheel_velocity < 0));
            right_wheel_velocity = abs(wheel_velocity.*(wheel_velocity > 0));
            wheel_bias = (nansum(right_wheel_velocity)-nansum(left_wheel_velocity))/ ...
                (nansum(right_wheel_velocity)+nansum(left_wheel_velocity));
            
            % Get reaction times for each stimulus
            full_stim_set = sort(unique([1,0.5,0.25,0.125,0.06,0]'*[1,-1]));
            
            trial_stim = trial_contrast(response_trials).*trial_side(response_trials);
            [~,trial_stim_idx] = ismember(trial_stim,full_stim_set);
            stim_rxn_time = accumarray(trial_stim_idx(response_trials)',trial_move_t',[11,1],@nanmedian,nan);
            
            % Performance and reaction time
            if contains(curr_protocol,'vanilla')
              
                performance = block.events.sessionPerformanceValues(:,end-10:end);
            else
                % Not sure whether this is right yet
                non_repeat_trials = block.events.repeatNumValues(response_trials) == 1;
                                
                n_trials_subset = accumarray(trial_stim_idx(non_repeat_trials)', ...
                    1,[11,1],@sum,nan);
                go_left_subset = accumarray(trial_stim_idx(non_repeat_trials)', ...
                    block.events.responseValues(non_repeat_trials)' == 1,[11,1],@sum,nan);
                                              
                performance = nan(3,11);
                performance(1,:) = full_stim_set;
                performance(2,:) = n_trials_subset;
                performance(3,:) = go_left_subset;
            end
            
            % Get whether all contrasts were used
            use_all_contrasts = all(block.events.useContrastsValues(end-10:end));
            
            % Store in behavior structure
            bhv.protocol{curr_day} = curr_protocol;
            bhv.session_duration(curr_day) = session_duration;
            bhv.n_trials(curr_day) = n_trials;
            bhv.total_water(curr_day) = total_water;
            bhv.wheel_velocity(curr_day) = nansum(abs(wheel_velocity));
            bhv.wheel_bias(curr_day) = wheel_bias;
            bhv.conditions = performance(1,:);
            bhv.n_trials_condition(curr_day,:) = performance(2,:);
            bhv.go_left_trials(curr_day,:) = performance(end,:);
            bhv.use_all_contrasts(curr_day) = use_all_contrasts;
            bhv.stim_rxn_time(curr_day,:) = stim_rxn_time;
            
            AP_print_progress_fraction(curr_day,length(experiments));
        end
        
    end
    
    % Plot summary
    day_num = cellfun(@(x) datenum(x),{experiments.day});
    day_labels = cellfun(@(day,protocol) [protocol(19:end) ' ' day(6:end)], ...
        {experiments.day},bhv.protocol,'uni',false);
    figure('Name',animal)
    
    % Trials and water
    subplot(3,2,1);
    yyaxis left
    % plot(day_num,bhv.n_trials./bhv.session_duration,'linewidth',2);
    % ylabel('Trials/min');
    plot(day_num,bhv.n_trials,'linewidth',2);
    ylabel('Trials');
    yyaxis right
    plot(day_num,bhv.total_water,'linewidth',2);
    ylabel('Total water');
    xlabel('Session');
    set(gca,'XTick',day_num);
    set(gca,'XTickLabel',day_labels);
    set(gca,'XTickLabelRotation',90);
    
    imaging_days = day_num([experiments.imaging]);
    for i = 1:length(imaging_days)
        line(repmat(imaging_days(i),1,2),ylim,'color','k');
    end
    
    ephys_days = day_num([experiments.ephys]);
    for i = 1:length(ephys_days)
        line(repmat(ephys_days(i),1,2),ylim,'color','r','linestyle','--');
    end
    
    % Wheel movement and bias
    subplot(3,2,2);
    yyaxis left
    % plot(day_num,bhv.wheel_velocity./bhv.session_duration,'linewidth',2);
    % ylabel('Wheel movement / min');
    plot(day_num,bhv.wheel_velocity,'linewidth',2);
    ylabel('Wheel movement');
    yyaxis right
    plot(day_num,bhv.wheel_bias,'linewidth',2);
    ylim([-1,1]);
    line(xlim,[0,0]);
    ylabel('Wheel bias');
    xlabel('Session');
    set(gca,'XTick',day_num);
    set(gca,'XTickLabel',day_labels);
    set(gca,'XTickLabelRotation',90);
    
    imaging_days = day_num([experiments.imaging]);
    for i = 1:length(imaging_days)
        line(repmat(imaging_days(i),1,2),ylim,'color','k');
    end
    
    ephys_days = day_num([experiments.ephys]);
    for i = 1:length(ephys_days)
        line(repmat(ephys_days(i),1,2),ylim,'color','r','linestyle','--');
    end
    
    % Psychometric over all days
    subplot(3,2,3);
    im = imagesc(bhv.conditions,1:size(bhv.go_left_trials),bhv.go_left_trials./bhv.n_trials_condition);
    set(im,'AlphaData',~isnan(get(im,'CData')));
    set(gca,'color',[0.5,0.5,0.5]);
    colormap(brewermap([],'*RdBu'));
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
    
    % Reaction time over days
    subplot(3,2,4);
    im = imagesc(bhv.conditions,1:size(bhv.go_left_trials),bhv.stim_rxn_time);
    set(im,'AlphaData',~isnan(get(im,'CData')));
    set(gca,'color',[0.5,0.5,0.5]);
    colormap(brewermap([],'*RdBu'));
    caxis([0.2,0.8])
    c = colorbar;
    ylabel(c,'Reaction time');
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
    subplot(3,2,5); hold on;
    combine_days = bhv.use_all_contrasts;
    combine_days_performance = bhv.go_left_trials(combine_days,:)./bhv.n_trials_condition(combine_days,:);
    errorbar(bhv.conditions,nanmean(combine_days_performance,1), ...
        nanstd(combine_days_performance,[],1)./sqrt(sum(~isnan(combine_days_performance))),'k','linewidth',2);
    xlim([-1,1]);
    ylim([0,1]);
    line(xlim,[0.5,0.5],'color','k','linestyle','--');
    line([0,0],ylim,'color','k','linestyle','--');
    axis square;
    xlabel('Condition');
    ylabel('Fraction go left');
    
    % Reaction time of combined days that use all contrasts
    subplot(3,2,6); hold on;
    combine_days = bhv.use_all_contrasts;
    combine_days_rxn_time = bhv.stim_rxn_time(combine_days,:);
    errorbar(bhv.conditions,nanmean(combine_days_rxn_time,1), ...
        nanstd(combine_days_rxn_time,[],1)./sqrt(sum(~isnan(combine_days_rxn_time))),'k','linewidth',2);
    xlim([-1,1]);
    line(xlim,[0.5,0.5],'color','k','linestyle','--');
    axis square;
    xlabel('Condition');
    ylabel('Reaction time');
    
    drawnow;
    
end

%% Get and plot single mouse behavior (ChoiceWorld (old) + vanillaChoiceWorld)
clear all;
animal = 'AP085';

protocol = 'ChoiceWorld';
choiceworld_experiments = AP_find_experiments(animal,protocol);
[choiceworld_experiments.protocol] = deal('ChoiceWorld');

protocol = 'vanillaChoiceworld';
flexible_name = true;
vanillachoiceworld_experiments = AP_find_experiments(animal,protocol,flexible_name);
[vanillachoiceworld_experiments.protocol] = deal('vanillaChoiceworld');

experiments = [choiceworld_experiments;vanillachoiceworld_experiments];

% Sort days (can be out of order if search across servers)
[~,sort_idx] = sort({experiments.day});
experiments = experiments(sort_idx);

bhv = struct;

for curr_day = 1:length(experiments)
    
    day = experiments(curr_day).day;
    experiment_num = experiments(curr_day).experiment;
    
    if length(experiment_num) > 1
        warning('(excluding > 1 exp/day)');disp(newline);
    end
    
    for curr_experiment = length(experiment_num)
        
        experiment = experiment_num(curr_experiment);
        
        [block_filename, block_exists] = AP_cortexlab_filename(animal,day,experiment,'block');
        
        % Load the block file
        load(block_filename)
        
        % Time of session (in minutes)
        session_duration = block.duration/60;
        
        switch experiments(curr_day).protocol       
            case 'ChoiceWorld'
                % Trial counts
                n_trials = length([block.trial.responseMadeTime]);
                n_rewards = sum([block.trial.feedbackType] == 1);
                
                % Time to response
                response_times = [block.trial(1:n_trials).responseMadeTime] - ...
                    [block.trial(1:n_trials).stimulusCueStartedTime];
                
                % Resample velocity over even/standard time intervals
                wheel_resample_rate = 1000;
                wheel_t_resample = block.inputSensorPositionTimes(1):1/wheel_resample_rate:block.inputSensorPositionTimes(end);
                wheel_values_resample = interp1(block.inputSensorPositionTimes,block.inputSensorPositions,wheel_t_resample);
                
            case 'vanillaChoiceworld'
                % Trial counts
                n_trials = length([block.events.responseTimes]);
                n_rewards = sum(block.events.hitValues);
                
                % Time to response
                response_times = [block.events.responseTimes(1:n_trials)] - ...
                    [block.events.stimOnTimes(1:n_trials)];
                
                % Resample velocity over even/standard time intervals
                wheel_resample_rate = 1000;
                wheel_t_resample = block.inputs.wheelTimes(1):1/wheel_resample_rate:block.inputs.wheelTimes(end);
                wheel_values_resample = interp1(block.inputs.wheelTimes,block.inputs.wheelValues,wheel_t_resample);
        end
             
        % Get wheel velocity
        wheel_smooth_t = 0.05; % seconds
        wheel_smooth_samples = round(wheel_smooth_t*wheel_resample_rate);
        wheel_velocity = interp1(conv(wheel_t_resample,[1,1]/2,'valid'), ...
            diff(smooth(wheel_values_resample,wheel_smooth_samples)),wheel_t_resample)';
        
        wheel_thresh = 0.025;        
        wheel_starts = wheel_t_resample(abs(wheel_velocity(1:end-1)) < wheel_thresh & ...
            abs(wheel_velocity(2:end)) > wheel_thresh);       
    
        % Wheel movements/biases
        left_wheel_velocity = abs(wheel_velocity.*(wheel_velocity < 0));
        right_wheel_velocity = abs(wheel_velocity.*(wheel_velocity > 0));
        wheel_bias = (nansum(right_wheel_velocity)-nansum(left_wheel_velocity))/ ...
            (nansum(right_wheel_velocity)+nansum(left_wheel_velocity));
                
        % Store in behavior structure
        bhv.rigname{curr_day} = block.rigName;
        bhv.session_duration(curr_day) = session_duration;
        bhv.n_trials(curr_day) = n_trials;
        bhv.n_rewards(curr_day) = n_rewards;
        bhv.response_time_med(curr_day) = median(response_times);
        bhv.wheel_velocity(curr_day) = nansum(abs(wheel_velocity));
        bhv.wheel_bias(curr_day) = wheel_bias;
     
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
end

% Plot summary
% (to have spaces on days without training)
day_num = cellfun(@(x) datenum(x),{experiments.day});
day_labels = cellfun(@(curr_day,curr_rig) [num2str(curr_day) ' ' curr_rig], ...
    num2cell(1:length(experiments)),bhv.rigname,'uni',false);

% % (to only plot days with training)
% day_num = 1:length(experiments);
% day_labels = 1:length(experiments);

figure('Name',animal)

% Trials and water
subplot(1,2,1);
yyaxis left
plot(day_num,bhv.n_trials,'linewidth',2);
ylabel('Trials');
yyaxis right
plot(day_num,bhv.n_rewards./bhv.n_trials,'linewidth',2);
line(xlim,[0.5,0.5])
ylim([0,1]);
ylabel('Frac correct');
xlabel('Session');
set(gca,'XTick',day_num);
set(gca,'XTickLabel',day_labels);
set(gca,'XTickLabelRotation',90);

yyaxis left; hold on
choiceworld_days = day_num(strcmp({experiments.protocol},'ChoiceWorld'));
scatter(choiceworld_days,repmat(0,1,length(choiceworld_days)),20,'r','filled');

% Wheel movement and bias
subplot(1,2,2,'YScale','log'); hold on
plot(day_num,bhv.response_time_med,'linewidth',2);
line(xlim,[0.5,0.5])
line(xlim,[1,1])
ylabel('Median response time');
xlabel('Session');
set(gca,'XTick',day_num);
set(gca,'XTickLabel',day_labels);
set(gca,'XTickLabelRotation',90);
choiceworld_days = day_num(strcmp({experiments.protocol},'ChoiceWorld'));
scatter(choiceworld_days,repmat(0.5,1,length(choiceworld_days)),20,'r','filled');


%% Get and plot learning rate across mice?

% % To use historical data on saved list
% load('training_list.mat');
% 
% valid_task = find(cellfun(@ischar,training_list(:,2)));
% use_tasks = {'Burgess AFC','Burgess NFC','vanillaChoiceworld','Multisensory'};
% use_task_mice = valid_task(ismember(training_list(valid_task,2),use_tasks));
% 
% animals = training_list(use_task_mice,1);

% To enter animals manually
animals = {'AP063','AP064','AP066','AP068','AP071','AP085','AP086','AP087'};

bhv = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    protocol = 'choiceworld';
    experiments = AP_find_experiments(animal,protocol,true);
    
%     % If there's more than 20 days, just load the first 20
%     if length(experiments) > 20
%         experiments = experiments(1:20);
%     end
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        % (if more than one of the same experiment - use the last)
        experiment = max(experiments(curr_day).experiment);
                
        [block_filename, block_exists] = AP_cortexlab_filename(animal,day,experiment,'block');
        
        % Load the block file
        load(block_filename)
        
        % Skip session if no 'duration' (because errored out?)
        if ~isfield(block,'duration')
            continue
        end
        
        % Get information from block
        
        % Time of session (in minutes)
        session_duration = block.duration/60;  
        
        if isfield(block,'expType')
            % Old ChoiceWorld/pre-Signals
            
            % Protocol
            [~,curr_protocol] = fileparts(block.expType);
            
            % Trial counts
            n_trials = length(block.trial);
            n_rewards = length(block.rewardDeliveredSizes);
            
            % Resample position over even/standard time intervals
            wheel_resample_rate = 1000;
            wheel_t_resample = block.inputSensorPositionTimes(1):1/wheel_resample_rate:block.inputSensorPositionTimes(end);
            wheel_values_resample = interp1(block.inputSensorPositionTimes,block.inputSensorPositions,wheel_t_resample);
            
        elseif isfield(block,'expDef')
            % Signals
            
            % Protocol
            [~,curr_protocol] = fileparts(block.expDef);
            
            % Trial counts
            n_trials = length(block.paramsValues);
            n_rewards = length(block.outputs.rewardValues);
            
            % Resample position over even/standard time intervals
            wheel_resample_rate = 1000;
            wheel_t_resample = block.inputs.wheelTimes(1):1/wheel_resample_rate:block.inputs.wheelTimes(end);
            wheel_values_resample = interp1(block.inputs.wheelTimes,block.inputs.wheelValues,wheel_t_resample);
        end        
       
        % (if there's only one trial - skip this day)
        if n_trials == 1
            continue
        end
   
        % Get wheel velocity
        wheel_smooth_t = 0.05; % seconds
        wheel_smooth_samples = round(wheel_smooth_t*wheel_resample_rate);
        wheel_velocity = interp1(conv(wheel_t_resample,[1,1]/2,'valid'), ...
            diff(smooth(wheel_values_resample,wheel_smooth_samples)),wheel_t_resample)';
        
        wheel_thresh = 0.025;
        wheel_starts = wheel_t_resample(abs(wheel_velocity(1:end-1)) < wheel_thresh & ...
            abs(wheel_velocity(2:end)) > wheel_thresh);
        
        % Wheel movements/biases
        left_wheel_velocity = abs(wheel_velocity.*(wheel_velocity < 0));
        right_wheel_velocity = abs(wheel_velocity.*(wheel_velocity > 0));
        wheel_bias = (nansum(right_wheel_velocity)-nansum(left_wheel_velocity))/ ...
            (nansum(right_wheel_velocity)+nansum(left_wheel_velocity));
        
%          % > Trial contrasts (one contrast per side)
%         trial_contrasts = cell2mat(arrayfun(@(x) ...
%             block.trial(x).condition.visCueContrast,1:length(block.trial),'uni',false));
%         trial_outcome = [block.trial(:).feedbackType];
%         
%         % Calculate trial info
%         
%         % Trial performance (100/50 on one side only)
%         left_easy_trials = trial_contrasts(1,:) >= 0.5 & trial_contrasts(2,:) == 0;
%         right_easy_trials = trial_contrasts(2,:) >= 0.5 & trial_contrasts(1,:) == 0;
%         
%         frac_left_correct = nanmean(trial_outcome(left_easy_trials(1:length(trial_outcome))) == 1);
%         frac_right_correct = nanmean(trial_outcome(right_easy_trials(1:length(trial_outcome))) == 1);
%         
%         frac_correct = [frac_left_correct;frac_right_correct];
        
        % Store in behavior structure
        bhv(curr_animal).name = animal;
        bhv(curr_animal).date(curr_day) = block.startDateTime;
        bhv(curr_animal).expDef{curr_day} = curr_protocol;
        bhv(curr_animal).session_duration(curr_day) = session_duration;
        bhv(curr_animal).n_trials(curr_day) = n_trials;
        bhv(curr_animal).n_rewards(curr_day) = n_rewards;
        bhv(curr_animal).wheel_velocity(curr_day) = nansum(abs(wheel_velocity));
        bhv(curr_animal).wheel_bias(curr_day) = wheel_bias;
%         bhv(curr_animal).frac_correct(:,curr_day) = frac_correct;
               
    end
    
    AP_print_progress_fraction(curr_animal,length(animals));
    
end

% Get rid of any empty data
% bhv(cellfun(@isempty,{bhv.name})) = [];

% % Plot summary (individual)
% figure;
% for curr_animal = 1:length(animals)
%     
%     % % (to have spaces on days without training)
%     % day_num = cellfun(@(x) datenum(x),{experiments.day});
%     % day_labels = cellfun(@(x) x(6:end),{experiments.day},'uni',false);
%     % (to only plot days with training)
%     day_num = 1:length(bhv(curr_animal).session_duration);
%     day_labels = 1:length(bhv(curr_animal).session_duration);
%     
%     % Trials and water
%     subplot(length(animals),3,(curr_animal-1)*3+1);
%     yyaxis left
%     plot(day_num,bhv(curr_animal).n_trials./bhv(curr_animal).session_duration,'linewidth',2);
%     ylabel('Trials/min');
%     yyaxis right
%     plot(day_num,bhv(curr_animal).n_rewards,'linewidth',2);
%     ylabel('Number of rewards');
%     xlabel('Session');
%     set(gca,'XTick',day_num);
%     set(gca,'XTickLabel',day_labels);
%     set(gca,'XTickLabelRotation',0);
%     
%     % Performance (on easy trials)
%     subplot(length(animals),3,(curr_animal-1)*3+2);
%     plot(day_num,bhv(curr_animal).frac_correct,'linewidth',2);
%     ylabel('Fraction correct');
%     legend({'Left','Right'})
%     xlabel('Session');
%     set(gca,'XTick',day_num);
%     set(gca,'XTickLabel',day_labels);
%     set(gca,'XTickLabelRotation',0);
%     
%     % Wheel movement and bias
%     subplot(length(animals),3,(curr_animal-1)*3+3);
%     yyaxis left
%     plot(day_num,bhv(curr_animal).wheel_velocity./bhv(curr_animal).session_duration,'linewidth',2);
%     ylabel('Wheel movement / min');
%     yyaxis right
%     plot(day_num,bhv(curr_animal).wheel_bias,'linewidth',2);
%     ylim([-1,1]);
%     ylabel('Wheel bias');
%     xlabel('Session');
%     set(gca,'XTick',day_num);
%     set(gca,'XTickLabel',day_labels);
%     set(gca,'XTickLabelRotation',0);   
% end

% % Plot summary (overlaid)
% figure;
% for curr_animal = 1:length(bhv)
%     
%     if isempty(bhv(curr_animal).date)
%        continue 
%     end
%     
%     % % (to have spaces on days without training)
%     % day_num = cellfun(@(x) datenum(x),{experiments.day});
%     % day_labels = cellfun(@(x) x(6:end),{experiments.day},'uni',false);
%     % (to only plot days with training)
%     day_num = 1:length(bhv(curr_animal).session_duration);
%     
%     % Trials and water
%     subplot(1,3,1); hold on;
%     plot(day_num,bhv(curr_animal).n_trials./bhv(curr_animal).session_duration,'k');
%     ylabel('Trials/min');
%     xlabel('Session');
%     set(gca,'XTick',day_num);
%     
%     % Performance (on easy trials)
%     subplot(1,3,2); hold on;
%     plot(day_num,bhv(curr_animal).frac_correct(1,:),'color','b');
%     plot(day_num,bhv(curr_animal).frac_correct(2,:),'color','r');
%     ylabel('Fraction correct');
%     xlabel('Session');
%     set(gca,'XTick',day_num);
%     
%     % Wheel movement and bias
%     subplot(1,3,3); hold on;
%     plot(day_num,bhv(curr_animal).wheel_velocity./bhv(curr_animal).session_duration,'k');
%     ylabel('Wheel movement / min');
%     xlabel('Session');
%     set(gca,'XTick',day_num);
%     
% end

% Pad and concatenate data
max_sessions = max(cellfun(@length,{bhv.date}));

trial_rate = cellfun(@(trials,duration) trials./duration,{bhv.n_trials},{bhv.session_duration},'uni',false);
trial_rate_pad = cell2mat(cellfun(@(x) padarray(x,[0,max_sessions - length(x)],NaN,'post'),trial_rate','uni',false));

n_rewards = cellfun(@(rewards,duration) rewards./duration,{bhv.n_rewards},{bhv.session_duration},'uni',false);
n_rewards_pad = cell2mat(cellfun(@(x) padarray(x,[0,max_sessions - length(x)],NaN,'post'),n_rewards','uni',false));

wheel_velocity = cellfun(@(wheel_velocity,duration) wheel_velocity./duration,{bhv.wheel_velocity},{bhv.session_duration},'uni',false);
wheel_velocity_pad = cell2mat(cellfun(@(x) padarray(x,[0,max_sessions - length(x)],NaN,'post'),wheel_velocity','uni',false));

training_dates = cellfun(@(x) x(1),{bhv.date});

figure;
subplot(2,3,1); hold on;
plot(trial_rate_pad','color',[0.5,0.5,0.5]);
plot(nanmedian(trial_rate_pad,1),'r','linewidth',3);
xlabel('Session');
ylabel('Trials/min');

subplot(2,3,2); hold on;
plot(n_rewards_pad','color',[0.5,0.5,0.5]);
plot(nanmedian(n_rewards_pad,1),'r','linewidth',3);
xlabel('Session');
ylabel('Rewards/min');

subplot(2,3,3); hold on;
plot(wheel_velocity_pad','color',[0.5,0.5,0.5]);
plot(nanmedian(wheel_velocity_pad,1),'r','linewidth',3);
xlabel('Session');
ylabel('Velocity/min');

avg_days = 5:7;
reward_avg = nanmean(n_rewards_pad(:,avg_days),2);
wheel_avg = nanmean(wheel_velocity_pad(:,avg_days),2);
subplot(2,1,2); hold on;
[~,date_sort_idx] = sort(training_dates);
set(gca,'XTick',1:length(bhv))
set(gca,'XTickLabel',{bhv.name})
set(gca,'XTickLabelRotation',90)

yyaxis left
plot(1:length(bhv),reward_avg(date_sort_idx));
ylabel(['Rewards/min days ' num2str(avg_days)])

yyaxis right
plot(1:length(bhv),wheel_avg(date_sort_idx));
ylabel(['Velocity/min days ' num2str(avg_days)])


%% Batch load behavior

% animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
% animals = {'AP063','AP064','AP066','AP068','AP071','AP085','AP086','AP087'};
animals = {'AP040','AP041','AP045','AP047','AP048'};

bhv = struct;

for curr_animal = 1:length(animals)
     
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    experiments = experiments([experiments.imaging] & ~[experiments.ephys]);
    
    % Restrict to experiments also with passive data
    passive_protocol = 'AP_lcrGratingPassive';
    passive_experiments = AP_find_experiments(animal,passive_protocol,true);
    passive_experiments = passive_experiments([passive_experiments.imaging]);
    use_experiments = ismember({experiments.day},{passive_experiments.day});
    experiments = experiments(use_experiments);
    
    load_parts.cam = false;
    load_parts.imaging = false;
    load_parts.ephys = false;   
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
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
            ~signals_events.repeatTrialValues(1:n_trials);
        
        % L/R choice
        go_right = (signals_events.trialSideValues(1:n_trials) == 1 & trial_outcome(1:n_trials) == -1) | ...
            (signals_events.trialSideValues(1:n_trials) == -1 & trial_outcome(1:n_trials) == 1);
        go_left = (signals_events.trialSideValues(1:n_trials) == -1 & trial_outcome(1:n_trials) == -1) | ...
            (signals_events.trialSideValues(1:n_trials) == 1 & trial_outcome(1:n_trials) == 1);
        
        trial_choice = go_right - go_left;
        
        if any(~ismember(trial_choice(use_trials),[-1,1]))
            error('Non-valid choice')
        end
               
        % Store in behavior structure
        bhv(curr_animal).animal = animal;
        bhv(curr_animal).days{curr_day} = day;
       
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
        
        clearvars('-except',preload_vars{:});
    end
    
    AP_print_progress_fraction(curr_animal,length(animals));
    
end

%% Get Dprime at >= 50% contrast (from above)

hit_rate_right_norm = cellfun(@(contrast,side,outcome) cellfun(@(contrast,side,outcome) ...
    (sum(outcome(contrast >= 0.5 & side == 1) == 1)+0.5)/ ...
    (sum(contrast > 0 & side == 1)+1),contrast,side,outcome), ...
    {bhv.trial_contrast},{bhv.trial_side},{bhv.trial_outcome},'uni',false);
miss_rate_left_norm = cellfun(@(contrast,side,outcome) cellfun(@(contrast,side,outcome) ...
    (sum(outcome(contrast >= 0.5 & side == -1) == -1)+0.5)/ ...
    (sum(contrast > 0 & side == -1)+1),contrast,side,outcome), ...
    {bhv.trial_contrast},{bhv.trial_side},{bhv.trial_outcome},'uni',false);

dprime = cellfun(@(hit,miss) norminv(hit)-norminv(miss), ...
    hit_rate_right_norm,miss_rate_left_norm,'uni',false);


%% Get side-slope (for ignoring side)
% (I don't think this worked - log(prob/1-prob) doesn't make straight line)
figure; hold on

stims_all = unique([-1,1]'*[0,6,25,12.5,25,50,100]/100)';
trials_n_all = cell(length(bhv),1);
go_left_n_all = cell(length(bhv),1);
zero_left_frac = nan(length(bhv),1);
stim_side_slopediff = nan(length(bhv),1);
for curr_animal = 1:length(bhv)
    for curr_day = 1:length(bhv(curr_animal).days)
    
    trial_side = bhv(curr_animal).trial_side{curr_day};
    trial_contrast = bhv(curr_animal).trial_contrast{curr_day};
    trial_choice = bhv(curr_animal).trial_choice{curr_day};
    
    [stims,stims_n,go_left_n] = ...
        grpstats(trial_choice == -1,trial_side.*trial_contrast, ...
        {'gname','numel','sum'});
    stims = cellfun(@str2num,stims);
    
    trials_n_all{curr_animal}(curr_day,ismember(stims_all,stims)) = stims_n;
    go_left_n_all{curr_animal}(curr_day,ismember(stims_all,stims)) = go_left_n;
    
    end
        
    % Fit lines to left/right side
    % (combine all days with all stim)
    use_days = all(trials_n_all{curr_animal},2);
    trials_n_usedays = sum(trials_n_all{curr_animal}(use_days,:),1);
    go_left_n_usedays = sum(go_left_n_all{curr_animal}(use_days,:),1);
    zero_left_frac(curr_animal) = ...
        sum(go_left_n_all{curr_animal}(use_days,stims_all == 0),1)./ ...
        sum(trials_n_all{curr_animal}(use_days,stims_all == 0),1);
    
    % (loglinear-normalize)
    soft_left_hitrate = (go_left_n_usedays+0.5)./(trials_n_usedays+1);
    log_frac_norm = log(soft_left_hitrate./(1-soft_left_hitrate));
    
    right_stim_fit = polyfit(stims_all(stims_all >= 0), ...
        log_frac_norm(stims_all >= 0),1);
    left_stim_fit = polyfit(stims_all(stims_all <= 0), ...
        log_frac_norm(stims_all <= 0),1);
    stim_side_slopediff(curr_animal) = ...
        right_stim_fit(1)-left_stim_fit(1);
    
    plot(stims_all,go_left_n_usedays./trials_n_usedays,'linewidth',2);
    
    AP_print_progress_fraction(curr_animal,length(bhv));
end

legend(num2str(stim_side_slopediff))




%% Plot behavior (from above)

animals = {bhv.animal};

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
min_correct = cellfun(@(l,r) min(l,r),l_correct,r_correct,'uni',false);
use_experiments = cellfun(@(x) x > thresh_correct,min_correct,'uni',false);

trial_choice_cat_use = arrayfun(@(x) horzcat(bhv(x).trial_choice{use_experiments{x}}),1:length(bhv),'uni',false);
trial_side_cat_use = arrayfun(@(x) horzcat(bhv(x).trial_side{use_experiments{x}}),1:length(bhv),'uni',false);
trial_contrast_cat_use = arrayfun(@(x) horzcat(bhv(x).trial_contrast{use_experiments{x}}),1:length(bhv),'uni',false);

trial_condition_cat_use = cellfun(@(side,contrast) side.*contrast,trial_side_cat_use,trial_contrast_cat_use,'uni',false);

frac_left_use = cell2mat(cellfun(@(choice,condition) ...
    grpstats(choice == -1,condition),trial_choice_cat_use,trial_condition_cat_use,'uni',false));

animal_experiments = cell2mat(cellfun(@(mouse,data) repmat(mouse,1,length(data)), ...
    num2cell(1:length(bhv)),l_correct,'uni',false));

figure;
subplot(1,2,1);
p1 = scatter(1:length(animal_experiments),horzcat(min_correct{:}), ...
    40,animal_experiments,'Filled','MarkerEdgeColor','k');
axis tight;
p2 = line(xlim,[thresh_correct,thresh_correct],'linewidth',3);
xlabel('Experiment');
ylabel(['Min performance easy contrasts']);
legend([p1(1),p2],{'Mouse','Performance threshold'});

subplot(1,2,2); hold on;
plot(conditions,frac_left,'linewidth',2,'color',[0.5,0,0]);
p1 = plot(conditions,nanmean(frac_left,2),'linewidth',5,'color','r');
plot(conditions,frac_left_use,'linewidth',2,'color',[0.5,0.5,0.5]);
p2 = plot(conditions,nanmean(frac_left_use,2),'linewidth',5,'color','k');
xlim([-1,1]);
ylim([0,1]);
line([0,0],ylim,'linestyle','--','color','k');
line(xlim,[0.5,0.5],'linestyle','--','color','k');
xlabel('Condition');
ylabel('Fraction go left');
legend([p1,p2],{'All','> Threshold'},'location','nw');

% Plot wheel velocity by correct/incorrect
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

% Plot wheel velocity by left/right
[velocity_condition_moveL_earlymove_grp,group] = cellfun(@(wheel,choice,condition,stim_to_move) ...
    grpstats(wheel(choice == -1 & stim_to_move < 0.5),condition(choice == -1 & stim_to_move < 0.5), ...
    {'mean','gname'}),trial_wheel_velocity_cat,trial_choice_cat,trial_condition_cat,stim_to_move_cat,'uni',false);
group = cellfun(@(x) cellfun(@(x) str2num(x),x),group,'uni',false);
all_group = unique(vertcat(group{:}),'rows');
velocity_condition_hit_earlymove = nan(size(all_group,1),length(bhv));
for curr_animal = 1:length(bhv)
    [~,group_idx] = intersect(all_group,group{curr_animal},'rows');
    velocity_condition_moveL_earlymove(group_idx,curr_animal) = velocity_condition_moveL_earlymove_grp{curr_animal};
end

[velocity_condition_moveR_earlymove_grp,group] = cellfun(@(wheel,choice,condition,stim_to_move) ...
    grpstats(wheel(choice == 1 & stim_to_move < 0.5),condition(choice == 1 & stim_to_move < 0.5), ...
    {'mean','gname'}),trial_wheel_velocity_cat,trial_choice_cat,trial_condition_cat,stim_to_move_cat,'uni',false);
group = cellfun(@(x) cellfun(@(x) str2num(x),x),group,'uni',false);
all_group = unique(vertcat(group{:}),'rows');
velocity_condition_hit_earlymove = nan(size(all_group,1),length(bhv));
for curr_animal = 1:length(bhv)
    [~,group_idx] = intersect(all_group,group{curr_animal},'rows');
    velocity_condition_moveR_earlymove(group_idx,curr_animal) = velocity_condition_moveR_earlymove_grp{curr_animal};
end

figure; hold on;
plot(conditions,nanmean(abs(velocity_condition_moveL_earlymove),2),'color',[0.6,0,0.6],'linewidth',2);
plot(conditions,nanmean(abs(velocity_condition_moveR_earlymove),2),'color',[0,0.6,0],'linewidth',2);
ylabel('Max wheel velocity towards choice');
xlabel('Condition');

legend({'Move L','Move R'});

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
   
% Plot distribution of stim to move and feedback times (correct/incorrect)
move_time_bins = 0:0.01:1;
day_split = 4;

move_time_centers = move_time_bins(1:end-1) + diff(move_time_bins)/2;
stim_to_move_binned_correct = zeros(day_split,length(move_time_bins)-1,length(bhv));
stim_to_move_binned_incorrect = zeros(day_split,length(move_time_bins)-1,length(bhv));
stim_to_move_binned_zero = zeros(day_split,length(move_time_bins)-1,length(bhv));

for curr_animal = 1:length(bhv)
    for curr_day = 1:length(bhv(curr_animal).stim_to_move)
        curr_data = bhv(curr_animal).stim_to_move{curr_day}( ...
            bhv(curr_animal).trial_side{curr_day} == -bhv(curr_animal).trial_choice{curr_day} & ...
            bhv(curr_animal).trial_contrast{curr_day} ~= 0);
        trials_split = round(linspace(1,length(curr_data),day_split+1));
        for curr_split = 1:day_split
            stim_to_move_binned_correct(curr_split,:,curr_animal) = ...
                stim_to_move_binned_correct(curr_split,:,curr_animal) + ...
                histcounts(curr_data(trials_split(curr_split): ...
                trials_split(curr_split+1)),move_time_bins);
        end
        
        curr_data = bhv(curr_animal).stim_to_move{curr_day}( ...
            bhv(curr_animal).trial_side{curr_day} == bhv(curr_animal).trial_choice{curr_day} & ...
            bhv(curr_animal).trial_contrast{curr_day} ~= 0);
        trials_split = round(linspace(1,length(curr_data),day_split+1));
        for curr_split = 1:day_split
            stim_to_move_binned_incorrect(curr_split,:,curr_animal) = ...
                stim_to_move_binned_incorrect(curr_split,:,curr_animal) + ...
                histcounts(curr_data(trials_split(curr_split): ...
                trials_split(curr_split+1)),move_time_bins);
        end
        
        curr_data = bhv(curr_animal).stim_to_move{curr_day}( ...
            bhv(curr_animal).trial_contrast{curr_day} == 0);
        trials_split = round(linspace(1,length(curr_data),day_split+1));
        for curr_split = 1:day_split
            stim_to_move_binned_zero(curr_split,:,curr_animal) = ...
                stim_to_move_binned_zero(curr_split,:,curr_animal) + ...
                histcounts(curr_data(trials_split(curr_split): ...
                trials_split(curr_split+1)),move_time_bins);
        end
    end
end
stim_to_move_binned_correct_norm = bsxfun(@rdivide,stim_to_move_binned_correct, ...
    sum(sum(stim_to_move_binned_correct,1),2));
stim_to_move_binned_incorrect_norm = bsxfun(@rdivide,stim_to_move_binned_incorrect, ...
    sum(sum(stim_to_move_binned_incorrect,1),2)); 
stim_to_move_binned_zero_norm = bsxfun(@rdivide,stim_to_move_binned_zero, ...
    sum(sum(stim_to_move_binned_zero,1),2)); 

figure; hold on;
for curr_animal = 1:length(bhv)
    AP_stackplot(stim_to_move_binned_correct_norm(:,:,curr_animal)', ...
        move_time_centers,0.02,false,[0.5,1,0.5],1:day_split);
end
for curr_animal = 1:length(bhv)
    AP_stackplot(stim_to_move_binned_incorrect_norm(:,:,curr_animal)', ...
        move_time_centers,0.02,false,[1,0.5,0.5],1:day_split);
end
for curr_animal = 1:length(bhv)
    AP_stackplot(stim_to_move_binned_zero_norm(:,:,curr_animal)', ...
        move_time_centers,0.02,false,[0.5,0.5,0.5],1:day_split);
end

p1 = AP_stackplot(nanmean(stim_to_move_binned_correct_norm,3)', ...
    move_time_centers,0.02,false,[0,0.8,0],1:day_split);
p2 = AP_stackplot(nanmean(stim_to_move_binned_incorrect_norm,3)', ...
    move_time_centers,0.02,false,[0.8,0,0],1:day_split);
p3 = AP_stackplot(nanmean(stim_to_move_binned_zero_norm,3)', ...
    move_time_centers,0.02,false,[0,0,0],1:day_split);

xlabel('Time to movement onset')
ylabel('Frequency by fraction within day')
legend([p1(1),p2(1),p3(1)],{'Correct','Incorrect','Zero'});
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
condition_ids = combvec(contrasts,sides,choices,timings);
conditions_str = cellfun(@(x) sprintf('%+g/',x),mat2cell(condition_ids,size(condition_ids,1),ones(1,size(condition_ids,2))),'uni',false);
n_conditions = size(condition_ids,2);

condition_count = arrayfun(@(x) nan(n_conditions,length(bhv(x).session_duration)),1:length(bhv),'uni',false);
for curr_animal = 1:length(bhv);
    for curr_day = 1:length(bhv(curr_animal).session_duration)
        trial_conditions = ...
            [bhv(curr_animal).trial_contrast{curr_day}; bhv(curr_animal).trial_side{curr_day}; ...
            bhv(curr_animal).trial_choice{curr_day}; trial_timing{curr_animal}{curr_day}];
        [~,trial_id] = ismember(trial_conditions',condition_ids','rows');
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

    
%% Plot psychometric and reaction time of only included data

% Load behavior
bhv_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\bhv_processing\bhv.mat';
load(bhv_fn);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);
use_animals = cellfun(@any,use_experiments);

for curr_animal = 1:length(bhv)
    % messy but whatever: chage animal/conditions field then put back
    temp_name = bhv(curr_animal).animal;
    temp_conditions = bhv(curr_animal).conditions;
    bhv(curr_animal).animal = false(size(use_experiments{curr_animal}));
    
    bhv(curr_animal) = structfun(@(x) x(use_experiments{curr_animal}), ...
        bhv(curr_animal),'uni',false);
    
    bhv(curr_animal).animal = temp_name;
    bhv(curr_animal).conditions = temp_conditions;
end
bhv(~use_animals) = [];

animals = {bhv.animal};

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



%% Batch load and save TOTAL behavior (to look for biases on zero-contrast)

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
        bhv(curr_animal).animal = animal;
        bhv(curr_animal).days{curr_day} = day;
       
        bhv(curr_animal).session_duration(curr_day) = session_duration;
        bhv(curr_animal).n_trials(curr_day) = n_trials;
        bhv(curr_animal).total_water(curr_day) = total_water;
        bhv(curr_animal).conditions = performance(1,:);
                        
        bhv(curr_animal).trial_side{curr_day} = signals_events.trialSideValues;
        bhv(curr_animal).trial_contrast{curr_day} = signals_events.trialContrastValues;
        bhv(curr_animal).trial_choice{curr_day} = trial_choice;
        bhv(curr_animal).trial_outcome{curr_day} = trial_outcome;
        bhv(curr_animal).stim_to_move{curr_day} = stim_to_move;
        bhv(curr_animal).stim_to_feedback{curr_day} = stim_to_feedback;
        bhv(curr_animal).trial_wheel_velocity{curr_day} = trial_wheel_velocity;
        
        clearvars -except animals protocol experiments load_parts curr_animal curr_day animal bhv 
    end
    
    AP_print_progress_fraction(curr_animal,length(animals));
    
end

% Save behavior
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\bhv_processing';
save_fn = 'bhv_alltrials';
save([save_path filesep save_fn],'bhv');

    
%% Check for history bias in zero-contrast trials    
    
% Load behavior
bhv_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\bhv_processing\bhv_alltrials.mat';
load(bhv_fn);

animals = {bhv.animal};

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


% Predict responses from linear combination of past choices?
% (do this smarter with log-likelihood later?)
n_history_trials = 10;
history_kernel = nan(n_history_trials,length(animals));
history_choice_mean = nan(n_history_trials,length(animals),2);
history_correct_prediction = nan(length(animals),1);
for curr_animal = 1:length(animals)
    zero_trials = find(trial_condition_cat{curr_animal} == 0)';
    zero_trials(zero_trials < n_history_trials+1) = [];
    
    zero_choices = trial_choice_cat{curr_animal}(zero_trials)';
    history_trials = bsxfun(@plus,zero_trials,-n_history_trials:-1);
    history_choices = trial_choice_cat{curr_animal}(history_trials);
    history_outcome = trial_outcome_cat{curr_animal}(history_trials);
    
    % (use only correct trials)
    history_choices = history_choices.*(history_outcome == 1);
    
    history_kernel(:,curr_animal) = history_choices\zero_choices;
    history_choice_mean(:,curr_animal,1) = ...
        nanmean(history_choices(zero_choices == -1,:) == -1,1);
    history_choice_mean(:,curr_animal,2) = ...
        nanmean(history_choices(zero_choices == 1,:) == -1,1);
    
    history_prediction = history_choices*history_kernel(:,curr_animal) < 0.5;
    history_correct_prediction(curr_animal) = ...
        sum(history_prediction & zero_choices == -1)./sum(zero_choices ~= 0);
end

figure; 
subplot(1,3,1);hold on;
plot(-n_history_trials:-1,history_kernel);
plot(-n_history_trials:-1,nanmean(history_kernel,2),'k','linewidth',2);
xlabel('Past trial');
ylabel('Weight');
line(xlim,[0,0],'color','k')

subplot(1,3,2); 
plot(history_correct_prediction,'k','linewidth',2);
xlabel('Animal');
ylabel('Fraction correctly predicted');
line(xlim,[0.5,0.5],'color','k')

subplot(1,3,3); hold on;
plot(-n_history_trials:-1,history_choice_mean(:,:,1),'k');
plot(-n_history_trials:-1,history_choice_mean(:,:,2),'r');
plot(-n_history_trials:-1,nanmean(history_choice_mean(:,:,1),2),'k','linewidth',2);
plot(-n_history_trials:-1,nanmean(history_choice_mean(:,:,2),2),'r','linewidth',2);
xlabel('Past trial');
ylabel('Fraction move L');
line(xlim,[0.5,0.5],'color','k')


% As above, but all animals concatenated
n_history_trials = 40;
history_kernel = nan(n_history_trials,1);
history_choice_mean = nan(n_history_trials,2);
history_correct_prediction = nan;

a = horzcat(trial_condition_cat{:});
b = horzcat(trial_choice_cat{:});
% (use only correct trials)
c = horzcat(trial_outcome_cat{:});
b = b.*(c == 1);

zero_trials = find(a == 0)';
zero_trials(zero_trials < n_history_trials+1) = [];

zero_choices = b(zero_trials)';
history_trials = bsxfun(@plus,zero_trials,-n_history_trials:-1);
history_choices = b(history_trials);

history_kernel = history_choices\zero_choices;
history_choice_mean(:,1) = ...
    nanmean(history_choices(zero_choices == -1,:) == -1,1);
history_choice_mean(:,2) = ...
    nanmean(history_choices(zero_choices == 1,:) == -1,1);

history_prediction = history_choices*history_kernel < 0.5;
history_correct_prediction = ...
    sum(history_prediction & zero_choices == -1)./sum(zero_choices ~= 0);

figure; 
subplot(1,3,1);hold on;
plot(-n_history_trials:-1,history_kernel,'k','linewidth',2);
xlabel('Past trial');
ylabel('Weight');
line(xlim,[0,0],'color','k')

subplot(1,3,2); 
plot(history_correct_prediction,'.k','MarkerSize',20);
xlabel('Animal');
ylabel('Fraction correctly predicted');
line(xlim,[0.5,0.5],'color','k')

subplot(1,3,3); hold on;
plot(-n_history_trials:-1,history_choice_mean(:,1),'k','linewidth',2);
plot(-n_history_trials:-1,history_choice_mean(:,2),'r','linewidth',2);
xlabel('Past trial');
ylabel('Fraction move L');
line(xlim,[0.5,0.5],'color','k')


%% (testing) Get side-ignoring

clear all;
animal = 'AP085';

protocol = 'ChoiceWorld';
choiceworld_experiments = AP_find_experiments(animal,protocol);
[choiceworld_experiments.protocol] = deal('ChoiceWorld');

protocol = 'vanillaChoiceworld';
flexible_name = true;
vanillachoiceworld_experiments = AP_find_experiments(animal,protocol,flexible_name);
[vanillachoiceworld_experiments.protocol] = deal('vanillaChoiceworld');

experiments = [choiceworld_experiments;vanillachoiceworld_experiments];

% Sort days (can be out of order if search across servers)
[~,sort_idx] = sort({experiments.day});
experiments = experiments(sort_idx);

bhv = struct;

for curr_day = 1:length(experiments)
    
    day = experiments(curr_day).day;
    experiment_num = experiments(curr_day).experiment;
    
    if length(experiment_num) > 1
        warning('(excluding > 1 exp/day)');disp(newline);
    end
    
    for curr_experiment = length(experiment_num)
        
        experiment = experiment_num(curr_experiment);
        
        [block_filename, block_exists] = AP_cortexlab_filename(animal,day,experiment,'block');
        
        % Load the block file
        load(block_filename)
        
        % Time of session (in minutes)
        session_duration = block.duration/60;
        
        switch experiments(curr_day).protocol       
            case 'ChoiceWorld'
                % Trial counts
                n_trials = length([block.trial.responseMadeTime]);
                n_rewards = sum([block.trial.feedbackType] == 1);
                
                % Time to response
                response_times = [block.trial(1:n_trials).responseMadeTime] - ...
                    [block.trial(1:n_trials).stimulusCueStartedTime];
                
                % Resample velocity over even/standard time intervals
                wheel_resample_rate = 1000;
                wheel_t_resample = block.inputSensorPositionTimes(1):1/wheel_resample_rate:block.inputSensorPositionTimes(end);
                wheel_values_resample = interp1(block.inputSensorPositionTimes,block.inputSensorPositions,wheel_t_resample);
                
                % (No stim slope at the moment)
                stim_side_slopediff = NaN;
                               
            case 'vanillaChoiceworld'
                % Trial counts
                n_trials = length([block.events.responseTimes]);
                n_rewards = sum(block.events.hitValues);
                
                % Time to response
                response_times = [block.events.responseTimes(1:n_trials)] - ...
                    [block.events.stimOnTimes(1:n_trials)];
                
                % Resample velocity over even/standard time intervals
                wheel_resample_rate = 1000;
                wheel_t_resample = block.inputs.wheelTimes(1):1/wheel_resample_rate:block.inputs.wheelTimes(end);
                wheel_values_resample = interp1(block.inputs.wheelTimes,block.inputs.wheelValues,wheel_t_resample);
                
                
                % Fit lines to left/right performance, get slope difference
                performance = block.events.sessionPerformanceValues(:,end-10:end);                
                stim_conditions = performance(1,:);
                frac_orient_right = performance(3,:)./performance(2,:);
                if ~any(isnan(frac_orient_right))
                    right_stim_fit = polyfit(stim_conditions(stim_conditions >= 0), ...
                        frac_orient_right(stim_conditions >= 0),1);
                    left_stim_fit = polyfit(stim_conditions(stim_conditions <= 0), ...
                        frac_orient_right(stim_conditions <= 0),1);
                    stim_side_slopediff = right_stim_fit(1)-left_stim_fit(1);
                else
                    stim_side_slopediff = NaN;
                end       
                
        end   
             
        % Get wheel velocity
        wheel_smooth_t = 0.05; % seconds
        wheel_smooth_samples = round(wheel_smooth_t*wheel_resample_rate);
        wheel_velocity = interp1(conv(wheel_t_resample,[1,1]/2,'valid'), ...
            diff(smooth(wheel_values_resample,wheel_smooth_samples)),wheel_t_resample)';
        
        wheel_thresh = 0.025;        
        wheel_starts = wheel_t_resample(abs(wheel_velocity(1:end-1)) < wheel_thresh & ...
            abs(wheel_velocity(2:end)) > wheel_thresh);       
    
        % Wheel movements/biases
        left_wheel_velocity = abs(wheel_velocity.*(wheel_velocity < 0));
        right_wheel_velocity = abs(wheel_velocity.*(wheel_velocity > 0));
        wheel_bias = (nansum(right_wheel_velocity)-nansum(left_wheel_velocity))/ ...
            (nansum(right_wheel_velocity)+nansum(left_wheel_velocity));
                
        % Store in behavior structure
        bhv.rigname{curr_day} = block.rigName;
        bhv.session_duration(curr_day) = session_duration;
        bhv.n_trials(curr_day) = n_trials;
        bhv.n_rewards(curr_day) = n_rewards;
        bhv.response_time_med(curr_day) = median(response_times);
        bhv.wheel_velocity(curr_day) = nansum(abs(wheel_velocity));
        bhv.wheel_bias(curr_day) = wheel_bias;
        bhv.stim_side_slopediff(curr_day) = stim_side_slopediff;
     
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
end

    

%% --------- AP_stimWheelRight (aka one-sided choiceworld)

%% Get and plot mice behavior

% animals = {'AP089','AP090','AP091','AP092','AP093','AP094'};
% animals = {'AP092','AP093','AP094'};
% animals = {'AP095','AP096','AP097'};
animals = {'AP096'};
protocol = 'AP_stimWheelRight';
flexible_name = false;
bhv = struct;

for curr_animal = 1:length(animals)
    
    preload_vars = who;
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol,flexible_name);
            
    if isempty(experiments)
        disp(['No behavior data: ' animal]);
        continue
    end
    
    disp([animal ', day:'])
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment_num = experiments(curr_day).experiment;
         
        % If multiple experiments, only use the last one 
        % (usually multiple happens if mess ups and final one is good)
        for curr_experiment = length(experiment_num)
            
            experiment = experiment_num(curr_experiment);
            
            [block_filename, block_exists] = AP_cortexlab_filename(animal,day,experiment,'block');
            
            % Load the block file
            load(block_filename)
            
            % Get protocol name
            [~,curr_protocol] = fileparts(block.expDef);
            
            % Time of session (in minutes)
            session_duration = block.duration/60;
            
            % Trial counts
            n_trials = length(block.paramsValues);
            total_water = sum(block.outputs.rewardValues);

            % Get hit trial index (first is dud because of setup)
            rewarded_trials = block.events.hitValues(2:end) == 1;
            reward_times = block.events.hitTimes(find(rewarded_trials)+1);
            
            % Estimate reaction time
            % (evenly resample velocity - not even natively);
            wheel_resample_rate = 1000;
            wheel_t_resample = block.inputs.wheelTimes(1):1/wheel_resample_rate:block.inputs.wheelTimes(end);
            wheel_values_resample = interp1(block.inputs.wheelTimes,block.inputs.wheelValues,wheel_t_resample);
            
            [wheel_velocity,wheel_move] = AP_parse_wheel(wheel_values_resample,wheel_resample_rate);
            
            wheel_starts = wheel_t_resample(diff(wheel_move) == 1);
            wheel_stops = wheel_t_resample(diff(wheel_move) == -1);
                 
            response_trials = 1:length(block.events.responseValues);
            trial_wheel_starts = arrayfun(@(x) ...
                wheel_starts(find(wheel_starts > block.events.stimOnTimes(x),1)), ...
                response_trials);
            trial_move_t = trial_wheel_starts - block.events.stimOnTimes(response_trials);
            stim_move_t = trial_move_t(rewarded_trials);

            % (TO DO: get probability of movement within window around
            % stim)
            stim_surround_t_edges = -5:0.1:5;
            stim_surround_t_centers = stim_surround_t_edges(1:end-1) + diff(stim_surround_t_edges)/2;
            stim_surround_times = block.events.stimOnTimes' + stim_surround_t_edges;
            stim_surround_move = cell2mat(arrayfun(@(x) ...
                histcounts(wheel_starts,stim_surround_times(x,:)), ...
                1:size(stim_surround_times,1),'uni',false)');
            % (need to do some kind of hazard rate here, maybe like what
            % Ruediger did, I think there was another one like this too)
            
            % (temp: any movement, not just starts)
            stim_surround_t_centers = -10:0.1:10;
            stim_surround_times = block.events.stimOnTimes' + stim_surround_t_centers;
            stim_surround_move = interp1(wheel_t_resample,wheel_move, ...
                stim_surround_times,'previous');

            % (testing/todo: find  movement breaks during the iti and see
            % what time to restart movement normally is?
            % (from stimOff to stimOn)?
            %             block.events.stimOnTimes(2:end)
%             wheel_stops = wheel_t_resample(abs(wheel_velocity(1:end-1)) > wheel_thresh & ...
%                 abs(wheel_velocity(2:end)) < wheel_thresh);
% 
%             wheel_move = abs(wheel_velocity) > wheel_thresh;
%             
%             poststim_move_t = arrayfun(@(x) wheel_t_resample(find(wheel_move & ...
%                 wheel_t_resample' > block.events.stimOnTimes(x),1,'first')), ...
%                 find(response_trials)) - block.events.stimOnTimes(response_trials);
%             
%             prestim_move_t = arrayfun(@(x) wheel_t_resample(find(wheel_move & ...
%                 wheel_t_resample' < block.events.stimOnTimes(x),1,'last')), ...
%                 find(response_trials)) - block.events.stimOnTimes(response_trials);

            % (or another way: align by breaks within ITI, eg if mouse took
            % a break from seconds 5-7, what's the profile of moving again?
            % it's probably shallow, while the stim-related one is strong)
            test_bin = [5,7];
            trial_iti_time = wheel_t_resample - ...
                interp1(block.events.stimOnTimes,block.events.stimOnTimes, ...
                wheel_t_resample,'previous');
            trial_iti_time_use = trial_iti_time >= test_bin(1) & ...
                trial_iti_time < test_bin(2);
            %(not done - I don't think necessary)



            % Time to reward
            stim_reward_t = reward_times - ...
                block.events.stimOnTimes(rewarded_trials);

            % Resetting quiescence period for each trial
            if isfield('trialQuiescencevalues',block.events)
                quiescence_t = block.events.trialQuiescenceValues(response_trials);
            else
                quiescence_t = 0.5*ones(size(block.events.stimOnTimes));
            end
            
            % ITI time (including quiescence resets)
            iti_t = block.events.stimOnTimes(response_trials(2:end)) - ...
                block.events.stimOffTimes(response_trials(1:end-1));
            
            % Wheel movements/biases
            left_wheel_velocity = abs(wheel_velocity.*(wheel_velocity < 0));
            right_wheel_velocity = abs(wheel_velocity.*(wheel_velocity > 0));
            wheel_bias = (nansum(right_wheel_velocity)-nansum(left_wheel_velocity))/ ...
                (nansum(right_wheel_velocity)+nansum(left_wheel_velocity));

            % Store in behavior structure
            bhv(curr_animal).animal = animal;
            bhv(curr_animal).day{curr_day} = day;
            bhv(curr_animal).protocol{curr_day} = curr_protocol;
            bhv(curr_animal).session_duration(curr_day) = session_duration;
            bhv(curr_animal).n_trials(curr_day) = n_trials;
            bhv(curr_animal).total_water(curr_day) = total_water;
            bhv(curr_animal).wheel_velocity(curr_day) = nansum(abs(wheel_velocity));
            bhv(curr_animal).stim_move_t{curr_day} = stim_move_t;
            bhv(curr_animal).stim_reward_t{curr_day} = stim_reward_t;
            bhv(curr_animal).quiescence_t{curr_day} = quiescence_t;
            bhv(curr_animal).iti_t{curr_day} = iti_t;
            bhv(curr_animal).wheel_bias(curr_day) = wheel_bias;

            bhv(curr_animal).stim_surround_t = stim_surround_t_centers;
            bhv(curr_animal).stim_surround_wheel{curr_day} = stim_surround_move;
            
            AP_print_progress_fraction(curr_day,length(experiments));
        end
        
    end
    
    % Plot summary
    day_num = cellfun(@(x) datenum(x),{experiments.day});
    day_labels = cellfun(@(day,protocol) [day(6:end)], ...
        {experiments.day},bhv(curr_animal).protocol,'uni',false);
    
    [unique_protocols,~,protocol_idx] = unique(bhv(curr_animal).protocol);
    protocol_col = hsv(length(unique_protocols));
    
    figure('Name',animal)
    
    % Trials and water
    subplot(2,3,1); hold on;
    yyaxis left
    % plot(day_num,bhv(curr_animal).n_trials./bhv(curr_animal).session_duration,'linewidth',2);
    % ylabel('Trials/min');
    plot(day_num,bhv(curr_animal).n_trials,'linewidth',2);
    ylabel('Trials');
    yyaxis right
    plot(day_num,bhv(curr_animal).total_water,'linewidth',2);
    ylabel('Total water');
    xlabel('Session');
    set(gca,'XTick',day_num);
    set(gca,'XTickLabel',day_labels);
    set(gca,'XTickLabelRotation',90);
    
    protocol_plot = gscatter(day_num,zeros(size(day_num)),[bhv(curr_animal).protocol]');
    legend off;
    
    imaging_days = day_num([experiments.imaging]);
    for i = 1:length(imaging_days)
        line(repmat(imaging_days(i),1,2),ylim,'color','k');
    end
    
    ephys_days = day_num([experiments.ephys]);
    for i = 1:length(ephys_days)
        line(repmat(ephys_days(i),1,2),ylim,'color','r','linestyle','--');
    end
    
    % Wheel movement and bias
    subplot(2,3,2);
    yyaxis left
    % plot(day_num,bhv(curr_animal).wheel_velocity./bhv(curr_animal).session_duration,'linewidth',2);
    % ylabel('Wheel movement / min');
    plot(day_num,bhv(curr_animal).wheel_velocity,'linewidth',2);
    ylabel('Wheel movement');
    yyaxis right
    plot(day_num,bhv(curr_animal).wheel_bias,'linewidth',2);
    ylim([-1,1]);
    line(xlim,[0,0]);
    ylabel('Wheel bias');
    xlabel('Session');
    set(gca,'XTick',day_num);
    set(gca,'XTickLabel',day_labels);
    set(gca,'XTickLabelRotation',90);
    
    imaging_days = day_num([experiments.imaging]);
    for i = 1:length(imaging_days)
        line(repmat(imaging_days(i),1,2),ylim,'color','k');
    end
    
    ephys_days = day_num([experiments.ephys]);
    for i = 1:length(ephys_days)
        line(repmat(ephys_days(i),1,2),ylim,'color','r','linestyle','--');
    end

    % Stim-to-reward time
    subplot(2,3,3);
    plot(day_num,cellfun(@nanmedian,bhv(curr_animal).stim_reward_t),'k','linewidth',2);
    ylabel('Stim to reward time (s)');
    xlabel('Session')
    set(gca,'XTick',day_num);
    set(gca,'XTickLabel',day_labels);
    set(gca,'XTickLabelRotation',90);

    % Move fraction heatmap
    stim_surround_wheel_cat = ...
        cell2mat(cellfun(@(x) nanmean(x,1),bhv(curr_animal).stim_surround_wheel,'uni',false)');
    subplot(2,3,4);
    imagesc(bhv(curr_animal).stim_surround_t,day_num,stim_surround_wheel_cat);
    colormap(gca,brewermap([],'Greys'));
    caxis([0,1]);
    xlabel('Time from stim (s)');
    set(gca,'YTick',day_num);
    set(gca,'YTickLabel',day_labels);
    
    % Move fraction lineplot
    subplot(2,3,5); hold on;
    set(gca,'ColorOrder',copper(size(stim_surround_wheel_cat,1)));
    plot(bhv(curr_animal).stim_surround_t,stim_surround_wheel_cat');
    ylabel([0,1]);
    ylabel('Fraction move');
    xlabel('Time from stim (s)');
    
    % Move fraction pre/post stim
    t_pre = bhv(curr_animal).stim_surround_t > -4 & bhv(curr_animal).stim_surround_t < -2;
    t_post = bhv(curr_animal).stim_surround_t > 0 & bhv(curr_animal).stim_surround_t < 2;
    move_pre_stim = ...
        cell2mat(cellfun(@(x) nanmean(reshape(x(:,t_pre),[],1)), ...
        bhv(curr_animal).stim_surround_wheel,'uni',false)');
    move_post_stim = ...
        cell2mat(cellfun(@(x) nanmean(reshape(x(:,t_post),[],1)), ...
        bhv(curr_animal).stim_surround_wheel,'uni',false)');
    
    subplot(2,3,6); hold on;
    plot(day_num,move_pre_stim,'-k','linewidth',2);
    plot(day_num,move_post_stim,'-r','linewidth',2);
    ylabel({'Fraction moving','(2s pre,post)'})
    ylim([0,1]);
    set(gca,'XTick',day_num);
    set(gca,'XTickLabel',day_labels);
    set(gca,'XTickLabelRotation',90);
    
    drawnow;
    clearvars('-except',preload_vars{:});
    
end

% End if only one mouse
if length(animals) == 1
    return
end

% If more than one mouse, plot group
% (only for min days across mice)
min_days = min(cellfun(@length,{bhv.n_trials}));

% Movement around stim
stim_surround_wheel_avg = ...
    cell2mat(permute(cellfun(@(x) cell2mat(cellfun(@(x) ...
    nanmean(x,1),x(1:min_days),'uni',false)'), ...
    {bhv.stim_surround_wheel},'uni',false),[1,3,2]));

stim_surround_t = bhv(1).stim_surround_t;
figure; hold on
set(gca,'ColorOrder',copper(min_days));
plot(stim_surround_t,nanmean(stim_surround_wheel_avg,3)','linewidth',2);
xlabel('Time from stim (s)');
ylabel('Fraction moving');


    




    




    





    



    



    



    

