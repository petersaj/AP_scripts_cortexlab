%% AP_stimWheelRight (aka one-sided choiceworld)


%% ~~~~~~~ Grab behavior

%% Plot behavior (for checking ongoing)

% animals = {'AP107'};
animals = {'AP113','AP114','AP115'};

% protocol = 'AP_stimWheel';
protocol = 'AP_stimWheelLeftReverse';
flexible_name = true;
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
            
            % Load (behavior/timeline only)
            load_parts.imaging = false;
            AP_load_experiment;
            
            % Get protocol name
            [~,curr_protocol] = fileparts(block.expDef);
            
            % Time of session (in minutes)
            session_duration = block.duration/60;
            
            % Trial counts (use only full trials with a response)
            total_water = sum(block.outputs.rewardValues);
            
            % Hit fraction
            frac_hit = nanmean(trial_outcome == 1);
            
            % (probability of any movement around stim)
            stim_surround_t_centers = -10:0.1:10;
            stim_surround_times = stimOn_times + stim_surround_t_centers;
            stim_surround_move = interp1(Timeline.rawDAQTimestamps, ...
                +wheel_move,stim_surround_times,'previous');
            
            % Resetting quiescence period for each trial
            if isfield(block.events,'trialQuiescenceValues')
                quiescence_t = block.events.trialQuiescenceValues(1:n_trials);
            else
                quiescence_t = 0.5*ones(size(block.events.stimOnTimes));
            end
            
            % ITI time (including quiescence resets)
            iti_t = block.events.stimOnTimes(2:n_trials) - ...
                block.events.stimOffTimes(1:n_trials-1);
            
            % Wheel movements/biases
            left_wheel_velocity = abs(wheel_velocity.*(wheel_velocity < 0));
            right_wheel_velocity = abs(wheel_velocity.*(wheel_velocity > 0));
            wheel_bias = (nansum(right_wheel_velocity)-nansum(left_wheel_velocity))/ ...
                (nansum(right_wheel_velocity)+nansum(left_wheel_velocity));
            
            % Store in behavior structure
            % (general timings)
            bhv(curr_animal).animal = animal;
            bhv(curr_animal).day{curr_day,1} = day;
            bhv(curr_animal).protocol{curr_day,1} = curr_protocol;
            bhv(curr_animal).session_duration(curr_day,1) = session_duration;
            bhv(curr_animal).n_trials(curr_day,1) = n_trials;
            bhv(curr_animal).frac_hit(curr_day,1) = frac_hit;
            bhv(curr_animal).total_water(curr_day,1) = total_water;
            bhv(curr_animal).wheel_velocity(curr_day,1) = nansum(abs(wheel_velocity));
            bhv(curr_animal).stim_move_t{curr_day,1} = stim_to_move;
            bhv(curr_animal).stim_feedback_t{curr_day,1} = stim_to_feedback;
            bhv(curr_animal).quiescence_t{curr_day,1} = quiescence_t';
            bhv(curr_animal).iti_t{curr_day,1} = iti_t';
            bhv(curr_animal).wheel_bias(curr_day,1) = wheel_bias;
            
            % (stim-aligned movement)
            bhv(curr_animal).stim_surround_t = stim_surround_t_centers;
            bhv(curr_animal).stim_surround_wheel{curr_day,1} = stim_surround_move;
            
            AP_print_progress_fraction(curr_day,length(experiments));
        end
        
    end
    
    % Plot summary
    day_num = cellfun(@(x) datenum(x),{experiments.day}');
    day_labels = cellfun(@(day,protocol) [day(6:end)], ...
        {experiments.day}',bhv(curr_animal).protocol,'uni',false);
    
    [unique_protocols,~,protocol_idx] = unique(bhv(curr_animal).protocol);
    protocol_col = hsv(length(unique_protocols));
    
    figure('Name',animal)
    
    % Trials and water
    subplot(2,3,1); hold on;
    yyaxis left
    % plot(day_num,bhv(curr_animal).n_trials./bhv(curr_animal).session_duration,'linewidth',2);
    % ylabel('Trials/min');
    plot(day_num,bhv(curr_animal).n_trials./bhv(curr_animal).session_duration,'linewidth',2);
    ylabel('Trials/min');
    yyaxis right
    plot(day_num,bhv(curr_animal).total_water,'linewidth',2);
    protocol_plot = gscatter(day_num,zeros(size(day_num)),[bhv(curr_animal).protocol]);
    legend off;
    ylabel('Total water');
    xlabel('Day');
    set(gca,'XTick',day_num);
    set(gca,'XTickLabel',day_labels);
    set(gca,'XTickLabelRotation',90);
    
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
    xlabel('Day');
    set(gca,'XTick',day_num);
    set(gca,'XTickLabel',day_labels);
    set(gca,'XTickLabelRotation',90);
    
    % Stim to move/reward time
    subplot(2,3,3); hold on;
    plot(day_num,cellfun(@nanmedian,bhv(curr_animal).stim_move_t),'k','linewidth',2);
    plot(day_num,cellfun(@nanmedian,bhv(curr_animal).stim_feedback_t),'b','linewidth',2);
    ylabel('Stim to (move/reward) time (s)');
    xlabel('Session')
    set(gca,'YScale','log')
    set(gca,'XTick',day_num);
    set(gca,'XTickLabel',day_labels);
    set(gca,'XTickLabelRotation',90);
    
    % Stim-aligned move fraction
    stim_surround_wheel_cat = ...
        cell2mat(cellfun(@(x) nanmean(x,1),bhv(curr_animal).stim_surround_wheel,'uni',false));
    subplot(2,4,5);
    imagesc(bhv(curr_animal).stim_surround_t,[],stim_surround_wheel_cat);
    colormap(gca,brewermap([],'Greys'));
    caxis([0,1]);
    xlabel('Time from stim (s)');
    set(gca,'YTick',1:length(day_num));
    set(gca,'YTickLabel',day_labels);
    
    subplot(2,4,6); hold on;
    set(gca,'ColorOrder',copper(size(stim_surround_wheel_cat,1)));
    plot(bhv(curr_animal).stim_surround_t,stim_surround_wheel_cat');
    ylim([0,1]);
    ylabel('Fraction move');
    xlabel('Time from stim (s)');
    
    % Reaction time histogram
    rxn_bins = [0:0.05:1];
    rxn_bin_centers = rxn_bins(1:end-1) + diff(rxn_bins)./2;
    rxn_hist = cell2mat(cellfun(@(x) ...
        histcounts(x,rxn_bins,'normalization','probability'), ...
        bhv(curr_animal).stim_move_t,'uni',false));
    
    subplot(2,4,7);
    imagesc(rxn_bin_centers,[],rxn_hist);
    colormap(gca,brewermap([],'Greys'));
    caxis([0,0.5]);
    set(gca,'YTick',1:length(day_num));
    set(gca,'YTickLabel',day_labels);
    xlabel('Reaction time');
    
    subplot(2,4,8); hold on;
    set(gca,'ColorOrder',copper(size(stim_surround_wheel_cat,1)));
    plot(rxn_bin_centers,rxn_hist');
    xlabel('Reaction time');
    ylabel('Fraction');
    
    drawnow;
    
    clearvars('-except',preload_vars{:});
    
end


%% Grab and save behavior

% early corticostriatal mice (fixed task timings)
% animal_group = 'cstr_fixtime';
% animals =  {'AP089','AP090','AP091'};

% % corticostriatal mice
% animal_group = 'cstr';
% animals = {'AP092','AP093','AP094','AP095','AP096','AP097'};

% % tetO mice
% animal_group = 'teto';
% animals = {'AP100','AP101','AP103','AP104','AP105', ...
%     'AP106','AP107','AP108','AP109','AP111','AP112'};
% protocol = 'AP_stimWheelRight';
% flexible_name = false;
% bhv = struct;

% reversal mice
animal_group = 'reversal';
animals = {'AP107','AP113','AP114','AP115'};
protocol = 'AP_stimWheel';
flexible_name = true;
bhv = struct;

% Flags
plot_flag = false;
save_bhv = true;

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
            
            % Load (behavior/timeline only)
            load_parts.imaging = false;
            AP_load_experiment;
            
            % Get protocol name
            [~,curr_protocol] = fileparts(block.expDef);
            
            % Time of session (in minutes)
            session_duration = block.duration/60;
            
            % Trial counts (use only full trials with a response)
            total_water = sum(block.outputs.rewardValues);
            
            % Hit fraction
            frac_hit = nanmean(trial_outcome == 1);
            
            % (probability of any movement around stim)
            stim_surround_t_centers = -10:0.1:10;
            stim_surround_times = stimOn_times + stim_surround_t_centers;
            stim_surround_move = interp1(Timeline.rawDAQTimestamps, ...
                +wheel_move,stim_surround_times,'previous');
            
            % Resetting quiescence period for each trial
            if isfield(block.events,'trialQuiescenceValues')
                quiescence_t = block.events.trialQuiescenceValues(1:n_trials);
            else
                quiescence_t = 0.5*ones(size(block.events.stimOnTimes));
            end
            
            % ITI time for each trial
            iti_t = block.events.trialITIValues;
            
            % Wheel movements/biases
            % mm in clicks from +hw.DaqRotaryEncoder, lilrig encoder = 100
            wheel_click2mm = 0.4869;
            wheel_mm = sum(abs(diff(wheel_position)))*wheel_click2mm;

            left_wheel_velocity = abs(wheel_velocity.*(wheel_velocity < 0));
            right_wheel_velocity = abs(wheel_velocity.*(wheel_velocity > 0));
            wheel_bias = (nansum(right_wheel_velocity)-nansum(left_wheel_velocity))/ ...
                (nansum(right_wheel_velocity)+nansum(left_wheel_velocity));
            
            %%%%%%%%% STIM RESPONSE VS. NULL
            
            t = Timeline.rawDAQTimestamps';
            
            % Get quiescence reset threshold
            % (only stored in script! whoops, this was dumb)
            expDef_fn = [fileparts(block_filename) filesep day '_' ...
                num2str(experiment) '_' animal '_expDef.m'];
            if ~exist(expDef_fn,'file')
                error('%s %s: no expDef.m',animal,day)
            end
            expDef_text = fileread(expDef_fn);
            [~,quiescThreshold_txt] = regexp(expDef_text, ...
                'quiescThreshold = (\d*)','match','tokens');
            quiescThreshold = str2num(quiescThreshold_txt{1}{1});
            
            % Get quiescence reset times
            %             % (real: from new trial onset)
            %             quiescence_reset_t = AP_find_stimWheel_quiescence;
            % (extrapolated: from response to response)
            quiescence_reset_t_extrap = AP_extrap_stimWheel_quiescence;
            
            alt_stimOn_times = cell(n_trials,1);
            alt_stimOn_trialparams = cell(n_trials,1);
            % (skip trial 1: no ITI and bad quiescence watch)
            for curr_trial = 2:n_trials
                
                % Pull out current trial times (last to next response)
                curr_trial_t_idx = t >= signals_events.responseTimes(curr_trial-1) & ...
                    t <= signals_events.responseTimes(curr_trial);
                curr_trial_t = t(curr_trial_t_idx);
                
                % Get quiescence reset times for different ITIs
                param_timestep = 0.1; % (hardcoded in expDef)
                possible_iti = max([block.paramsValues.itiMin]):param_timestep:max([block.paramsValues.itiMax]);
                possible_quiescence = max([block.paramsValues.quiescenceMin]):param_timestep:max([block.paramsValues.quiescenceMax]);
                
                t_from_quiescence_reset_trialitis = nan(length(curr_trial_t),length(possible_iti));
                for curr_possible_iti = 1:length(possible_iti)
                    curr_possible_itiend = signals_events.responseTimes(curr_trial-1) + ...
                        possible_iti(curr_possible_iti);
                    curr_quiescence_resets = sort([quiescence_reset_t_extrap;curr_possible_itiend]);
                    
                    t_from_quiescence_reset_trialitis(:,curr_possible_iti) = ...
                        curr_trial_t - interp1(curr_quiescence_resets, ...
                        curr_quiescence_resets,curr_trial_t,'previous','extrap');
                end
                
                % Find alternate stim times which would have given same first move
                
                % (getting possible iti + quiescence crosses)
                alt_iti_reached = ((t(curr_trial_t_idx) - curr_trial_t(1)) > possible_iti);
                alt_quiescence_reached = ...
                    t_from_quiescence_reset_trialitis > permute(possible_quiescence,[1,3,2]);
                
                % (get possible stim times as iti x quiescence grid)
                [alt_stim_value,alt_stim_idx] = max( ...
                    permute(alt_iti_reached & alt_quiescence_reached,[2,3,1]),[],3);
                alt_stim_t = curr_trial_t(alt_stim_idx);
                alt_stimOn_times_all = alt_stim_t(alt_stim_value);
                
                % (get alt stim times that would have resulted in the same
                % first movement since that's the measured value)
                stim_leeway = 0.1;
                curr_wheel_move_alt_stim_idx = ...
                    arrayfun(@(stim) find(wheel_starts > stim-stim_leeway,1,'first'), ...
                    alt_stimOn_times_all);
                use_alt_stimOn_times = ...
                    curr_wheel_move_alt_stim_idx == wheel_move_stim_idx(curr_trial);
                
                % (make sure that real parameters = real stim time: missed
                % wheel clicks sometimes give non-reproducible traces, in
                % which case the trial shouldn't be used for stats)
                curr_quiescence_idx = find(possible_quiescence == quiescence_t(curr_trial));
                curr_iti_idx = find(possible_iti == iti_t(curr_trial-1));
                curr_block_stimOn = signals_events.stimOnTimes(curr_trial);
                curr_alt_stim_offset = curr_block_stimOn - ...
                    alt_stim_t(curr_iti_idx,curr_quiescence_idx);
                if curr_alt_stim_offset > 0.01
                    continue
                end
                
                % (apply the block vs actual stim on time delay to all
                % times - note this is regularly backwards in time??)
                curr_stim_pd_offset = stimOn_times(curr_trial) - curr_block_stimOn;
                alt_stimOn_times_all_pd = alt_stimOn_times_all + curr_stim_pd_offset;
                
                % (store alternate stim times)
                alt_stimOn_times{curr_trial} = alt_stimOn_times_all_pd(use_alt_stimOn_times);
                
%                 % (trial plot)
%                 figure; hold on;
%                 t_plot_scale = 0.1;
%                 plot(t(curr_trial_t_idx),wheel_velocity(curr_trial_t_idx),'k')
%                 plot(t(curr_trial_t_idx),[0;diff(wheel_position(curr_trial_t_idx))]*0.1,'g')
%                 plot(alt_stimOn_times_all,0,'ob');
%                 plot(alt_stimOn_times{curr_trial},0,'.r');
%                 line(repmat(curr_trial_t(1)+signals_events.trialITIValues(curr_trial-1),2,1),ylim);
%                 line(xlim,repmat(signals_events.trialQuiescenceValues(curr_trial),2,1)*t_plot_scale,'color','m');
%                 line(repmat(curr_block_stimOn,1,2),ylim,'color','r','linestyle','--');
%                 line(repmat(stimOn_times(curr_trial),1,2),ylim,'color','k','linestyle','--');
%                 drawnow;
                
            end
            
            % Get would-be reaction time after alt stim times
            stim_leeway = 0.1;
            wheel_move_alt_stim_idx = ...
                arrayfun(@(stim) find(wheel_starts > stim-stim_leeway,1,'first'), ...
                cell2mat(alt_stimOn_times));
            
            alt_stim_to_move = ...
                mat2cell(wheel_starts(wheel_move_alt_stim_idx) - cell2mat(alt_stimOn_times), ...
                cellfun(@length,alt_stimOn_times));
            
            % Compare same stim/alt-stim trials (sub-sample alt-stim)
            use_alt_trials = cellfun(@(x) ~isempty(x),alt_stimOn_times);
            surround_t_centers = 0:0.01:2;
            
            % (day total)
            curr_trials = use_alt_trials;
            curr_stim_surround_times = stimOn_times(curr_trials) + surround_t_centers;
            curr_stim_surround_move = interp1(t,+wheel_move,curr_stim_surround_times,'previous');
            
            n_alt = 1000; % random subset alt trials for same n
            curr_alt_subset = cell2mat(cellfun(@(x) randsample(x,n_alt,true)', ...
                alt_stimOn_times(curr_trials),'uni',false));
            
            curr_alt_stim_surround_times = permute(curr_alt_subset,[1,3,2]) + surround_t_centers;
            curr_alt_stim_surround_move = interp1(t,+wheel_move,curr_alt_stim_surround_times,'previous');
            
            stim_move_max = max(nanmean(curr_stim_surround_move,1));
            alt_stim_move_max = nanmean(max(nanmean(curr_alt_stim_surround_move,1),[],2));
            
            stim_response_idx = (stim_move_max-alt_stim_move_max)./ ...
                (stim_move_max+alt_stim_move_max);
            
            % (daysplit)
            n_daysplit = 4;
            day_split_idx = min(floor(linspace(1,n_daysplit+1,n_trials)),n_daysplit)';
            
            stim_move_max_daysplit = nan(1,n_daysplit);
            alt_stim_move_max_daysplit = nan(1,n_daysplit);
            for curr_daysplit = 1:n_daysplit
                
                curr_trials = day_split_idx == curr_daysplit & use_alt_trials;
                
                curr_stim_surround_times = stimOn_times(curr_trials) + surround_t_centers;
                curr_stim_surround_move = interp1(t,+wheel_move,curr_stim_surround_times,'previous');
                
                n_alt = 1000; % random subset alt trials for same n
                curr_alt_subset = cell2mat(cellfun(@(x) randsample(x,n_alt,true)', ...
                    alt_stimOn_times(curr_trials),'uni',false));
                
                curr_alt_stim_surround_times = permute(curr_alt_subset,[1,3,2]) + surround_t_centers;
                curr_alt_stim_surround_move = interp1(t,+wheel_move,curr_alt_stim_surround_times,'previous');
                
                stim_move_max_daysplit(curr_daysplit) = max(nanmean(curr_stim_surround_move,1));
                alt_stim_move_max_daysplit(curr_daysplit) = nanmean(max(nanmean(curr_alt_stim_surround_move,1),[],2));
                
            end
            
            stim_response_idx_daysplit = (stim_move_max_daysplit-alt_stim_move_max_daysplit)./ ...
                (stim_move_max_daysplit+alt_stim_move_max_daysplit);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Store in behavior structure
            % (general values)
            bhv(curr_animal).animal = animal;
            bhv(curr_animal).day{curr_day,1} = day;
            bhv(curr_animal).protocol{curr_day,1} = curr_protocol;
            bhv(curr_animal).session_duration(curr_day,1) = session_duration;
            bhv(curr_animal).n_trials(curr_day,1) = n_trials;
            bhv(curr_animal).frac_hit(curr_day,1) = frac_hit;
            bhv(curr_animal).total_water(curr_day,1) = total_water;
            bhv(curr_animal).wheel_mm(curr_day,1) = wheel_mm;
            bhv(curr_animal).wheel_velocity(curr_day,1) = nansum(abs(wheel_velocity));
            bhv(curr_animal).stim_move_t{curr_day,1} = stim_to_move;
            bhv(curr_animal).stim_feedback_t{curr_day,1} = stim_to_feedback;
            bhv(curr_animal).quiescence_t{curr_day,1} = quiescence_t';
            bhv(curr_animal).iti_t{curr_day,1} = iti_t';
            bhv(curr_animal).wheel_bias(curr_day,1) = wheel_bias;
            
            % (stim-aligned movement)
            bhv(curr_animal).stim_surround_t = stim_surround_t_centers;
            bhv(curr_animal).stim_surround_wheel{curr_day,1} = stim_surround_move;
            
            % (stim vs null response)
            bhv(curr_animal).alt_stim_move_t{curr_day,1} = alt_stim_to_move;
            bhv(curr_animal).alt_stim_trialparams{curr_day,1} = alt_stimOn_trialparams;
            
            bhv(curr_animal).stim_move_max{curr_day,1} = stim_move_max;
            bhv(curr_animal).alt_stim_move_max{curr_day,1} = alt_stim_move_max;
            bhv(curr_animal).stim_response_idx(curr_day,1) = stim_response_idx;
            
            bhv(curr_animal).stim_move_max_daysplit{curr_day,1} = stim_move_max_daysplit;
            bhv(curr_animal).alt_stim_move_max_daysplit{curr_day,1} = alt_stim_move_max_daysplit;
            bhv(curr_animal).stim_response_idx_daysplit{curr_day,1} = stim_response_idx_daysplit;
            
            AP_print_progress_fraction(curr_day,length(experiments));
        end
        
    end
    
    if plot_flag
        
        % Plot summary
        day_num = cellfun(@(x) datenum(x),{experiments.day}');
        day_labels = cellfun(@(day,protocol) [day(6:end)], ...
            {experiments.day}',bhv(curr_animal).protocol,'uni',false);
        
        [unique_protocols,~,protocol_idx] = unique(bhv(curr_animal).protocol);
        protocol_col = hsv(length(unique_protocols));
        
        figure('Name',animal)
        
        % Trials and water
        subplot(2,3,1); hold on;
        yyaxis left
        % plot(day_num,bhv(curr_animal).n_trials./bhv(curr_animal).session_duration,'linewidth',2);
        % ylabel('Trials/min');
        plot(day_num,bhv(curr_animal).n_trials./bhv(curr_animal).session_duration,'linewidth',2);
        ylabel('Trials/min');
        yyaxis right
        plot(day_num,bhv(curr_animal).total_water,'linewidth',2);
        ylabel('Total water');
        xlabel('Session');
        set(gca,'XTick',day_num);
        set(gca,'XTickLabel',day_labels);
        set(gca,'XTickLabelRotation',90);
        
        protocol_plot = gscatter(day_num,zeros(size(day_num)),[bhv(curr_animal).protocol]);
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
        plot(day_num,cellfun(@nanmedian,bhv(curr_animal).stim_feedback_t),'k','linewidth',2);
        ylabel('Stim to reward time (s)');
        xlabel('Session')
        set(gca,'XTick',day_num);
        set(gca,'XTickLabel',day_labels);
        set(gca,'XTickLabelRotation',90);
        
        % Move fraction heatmap
        stim_surround_wheel_cat = ...
            cell2mat(cellfun(@(x) nanmean(x,1),bhv(curr_animal).stim_surround_wheel,'uni',false));
        subplot(2,3,4);
        imagesc(bhv(curr_animal).stim_surround_t,[],stim_surround_wheel_cat);
        colormap(gca,brewermap([],'Greys'));
        caxis([0,1]);
        xlabel('Time from stim (s)');
        set(gca,'YTick',1:length(day_num));
        set(gca,'YTickLabel',day_labels);
        
        % Move fraction lineplot
        subplot(2,3,5); hold on;
        set(gca,'ColorOrder',copper(size(stim_surround_wheel_cat,1)));
        plot(bhv(curr_animal).stim_surround_t,stim_surround_wheel_cat');
        ylim([0,1]);
        ylabel('Fraction move');
        xlabel('Time from stim (s)');
        
        % Move fraction pre/post stim
        stim_surround_t = bhv(curr_animal).stim_surround_t;
        move_prestim_max = ...
            cell2mat(cellfun(@(x) max(nanmean(x(:,stim_surround_t < 0),1)), ...
            bhv(curr_animal).stim_surround_wheel,'uni',false));
        move_poststim_max = ...
            cell2mat(cellfun(@(x) max(nanmean(x(:,stim_surround_t > 0),1)), ...
            bhv(curr_animal).stim_surround_wheel,'uni',false));
        move_prepost_max_ratio = ...
            (move_poststim_max-move_prestim_max)./(move_poststim_max+move_prestim_max);
        
        subplot(2,3,6); hold on;
        plot(day_num,move_prepost_max_ratio,'k','linewidth',2);
        ylabel({'Post:pre max move ratio'});
        ylim([-1,1]);
        line(xlim,[0,0],'color','k','linestyle','--');
        set(gca,'XTick',day_num);
        set(gca,'XTickLabel',day_labels);
        set(gca,'XTickLabelRotation',90);
        
        drawnow;
    end
    
    clearvars('-except',preload_vars{:});
    
end

if save_bhv
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    save_fn = [data_path filesep 'bhv_' animal_group];
    save(save_fn,'bhv');
    disp(['Saved ' save_fn]);
end


%% (load data saved after above)

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';

% % (load separate)
bhv_fn = [data_path filesep 'bhv_teto'];
% bhv_fn = [data_path filesep 'bhv_cstr'];
% bhv_fn = [data_path filesep 'bhv_reversal'];

load(bhv_fn);

% % (load and combine)
% bhv_fn1 = [data_path filesep 'bhv_teto'];
% bhv1 = load(bhv_fn1);
% bhv_fn2 = [data_path filesep 'bhv_cstr'];
% bhv2 = load(bhv_fn2);
% bhv = [bhv1.bhv,bhv2.bhv];
% clearvars bhv_fn1 bhv1 bhv_fn2 bhv2

animals = {bhv.animal};



%% ~~~~~~~ Plot behavior

%% Trial behavior (pre-muscimol or muscimol)

% Load muscimol injection info
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
muscimol_fn = [data_path filesep 'muscimol.mat'];
load(muscimol_fn);

% Use days before muscimol
use_days = cell(size(bhv));
for curr_animal = 1:length(bhv)
    muscimol_animal_idx = ismember({muscimol.animal},bhv(curr_animal).animal);
    if ~any(muscimol_animal_idx)
        use_days{curr_animal} = true(length(bhv(curr_animal).day),1);
        continue
    end
    
    pre_muscimol_day_idx = datenum(bhv(curr_animal).day) < ...
        datenum(muscimol(muscimol_animal_idx).day(1));
    muscimol_day_idx = ismember(datenum(bhv(curr_animal).day), ...
        datenum(muscimol(muscimol_animal_idx).day));
    
    use_days{curr_animal} = pre_muscimol_day_idx;
end

stim_move_t_cat = cellfun(@(x,y) cell2mat(x(y)),{bhv.stim_move_t},use_days,'uni',false);
stim_feedback_t_cat = cellfun(@(x,y) cell2mat(x(y)),{bhv.stim_feedback_t},use_days,'uni',false);

trials_cumsum = cellfun(@(x,y) cumsum(cellfun(@length,x(y))),{bhv(:).stim_move_t},use_days,'uni',false);

figure;
for curr_animal = 1:length(animals)
    subplot(length(animals),1,curr_animal); hold on;
    plot(stim_move_t_cat{curr_animal},'.k');
    plot(stim_feedback_t_cat{curr_animal},'.b');
    set(gca,'YScale','log')
    ylim([0.05,10])
    xlim([0,trials_cumsum{curr_animal}(end)]);
    for i = 1:length(trials_cumsum{curr_animal})
        line(repmat(trials_cumsum{curr_animal}(i),2,1), ...
            ylim,'color','r','linewidth',2)
    end
    ylabel(animals{curr_animal});
    
    muscimol_animal_idx = ismember({muscimol.animal},bhv(curr_animal).animal);
    if any(muscimol_animal_idx)
        muscimol_day_idx = ismember(datenum(bhv(curr_animal).day(use_days{curr_animal})), ...
            datenum(muscimol(muscimol_animal_idx).day));
        bhv_muscimol_day_idx = ismember(datenum(muscimol(muscimol_animal_idx).day), ...
            datenum(bhv(curr_animal).day(use_days{curr_animal})));
        if any(muscimol_day_idx)
            trial_n_center = diff([0,trials_cumsum{curr_animal}])./2 + ...
                [0,trials_cumsum{curr_animal}(1:end-1)];
            text(trial_n_center(muscimol_day_idx),repmat(12,1,sum(muscimol_day_idx)), ...
                muscimol(muscimol_animal_idx).area(bhv_muscimol_day_idx), ...
                'color','r','HorizontalAlign','center');
        end
    end
    
end


%% Histograms of reaction times

rxn_bins = [-Inf,0:0.02:0.5,Inf];
rxn_bin_centers = rxn_bins(1:end-1) + diff(rxn_bins)./2;

% Load muscimol injection info
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
muscimol_fn = [data_path filesep 'muscimol.mat'];
load(muscimol_fn);

% Use days before muscimol
use_days = cell(size(bhv));
for curr_animal = 1:length(bhv)
    muscimol_animal_idx = ismember({muscimol.animal},bhv(curr_animal).animal);
    if ~any(muscimol_animal_idx)
        use_days{curr_animal} = true(length(bhv(curr_animal).day),1);
        continue
    end
    muscimol_day_idx = datenum(bhv(curr_animal).day) >= ...
        datenum(muscimol(muscimol_animal_idx).day(1));
    use_days{curr_animal} = ~muscimol_day_idx;
end

max_days = max(cellfun(@sum,use_days));

% Animal/day
% (match null n just single sample)
alt_rxn_matched = ...
    cellfun(@(x) cellfun(@(x) cell2mat(cellfun(@(x) ...
    datasample(x,min(length(x),1)),x,'uni',false)),x,'uni',false), ...
    {bhv.alt_stim_move_t},'uni',false);
figure;
h = reshape(tight_subplot(length(bhv),max_days),max_days,length(bhv))';
for curr_animal = 1:length(bhv)
    for curr_day = find(use_days{curr_animal})'
        subplot(h(curr_animal,curr_day)); hold on
        histogram(bhv(curr_animal).stim_move_t{curr_day},rxn_bins,'EdgeColor','none','normalization','pdf')
        histogram(alt_rxn_matched{curr_animal}{curr_day}, ...
            rxn_bins,'EdgeColor','none','normalization','pdf')
        if curr_day == 1
           ylabel(bhv(curr_animal).animal); 
        end
        if curr_animal == 1
           title(sprintf('Day %d',curr_day)); 
        end
    end
end

% Day (animals overlaid)
figure;
for curr_animal = 1:length(bhv)
    for curr_day = 1:min(length(bhv(curr_animal).stim_move_t),10)
        subplot(max_days,1,curr_day); hold on;
        histogram(bhv(curr_animal).stim_move_t{curr_day},rxn_bins,'EdgeColor','none','normalization','pdf')
    end
end

% Total (single day)
plot_day = 6;

animal_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(vertcat(x{plot_day}),rxn_bins,'normalization','pdf'), ...
    use_days,{bhv.stim_move_t},'uni',false)');

animal_alt_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(cell2mat(vertcat(x{plot_day})),rxn_bins,'normalization','pdf'), ...
    use_days,{bhv.alt_stim_move_t},'uni',false)');

figure; hold on;
histogram('BinEdges',rxn_bins,'BinCounts',nanmean(animal_rxn,1),'EdgeColor','none')
histogram('BinEdges',rxn_bins,'BinCounts',nanmean(animal_alt_rxn,1),'EdgeColor','none')
xlabel('Reaction time');
ylabel('PDF');
legend({'Measured','Null'});
title(sprintf('Day %d',plot_day));

% Total (single day - 1:1 rxn)
plot_day = 1;

animal_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(vertcat(x{plot_day}),rxn_bins,'normalization','pdf'), ...
    use_days,{bhv.stim_move_t},'uni',false)');

animal_alt_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(cell2mat(vertcat(x{plot_day})),rxn_bins,'normalization','pdf'), ...
    use_days,{bhv.alt_stim_move_t},'uni',false)');

figure; hold on;
histogram('BinEdges',rxn_bins,'BinCounts',nanmean(animal_rxn,1),'EdgeColor','none')
histogram('BinEdges',rxn_bins,'BinCounts',nanmean(animal_alt_rxn,1),'EdgeColor','none')
xlabel('Reaction time');
ylabel('PDF');
legend({'Measured','Null'});
title(sprintf('Day %d',plot_day));

% Total (single animal)
plot_animal = 10;

animal_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(vertcat(x{:}),rxn_bins,'normalization','pdf'), ...
    use_days(plot_animal),{bhv(plot_animal).stim_move_t},'uni',false)');

animal_alt_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(cell2mat(vertcat(x{:})),rxn_bins,'normalization','pdf'), ...
    use_days(plot_animal),{bhv(plot_animal).alt_stim_move_t},'uni',false)');

figure; hold on;
histogram('BinEdges',rxn_bins,'BinCounts',nanmean(animal_rxn,1),'EdgeColor','none')
histogram('BinEdges',rxn_bins,'BinCounts',nanmean(animal_alt_rxn,1),'EdgeColor','none')
xlabel('Reaction time');
ylabel('PDF');
legend({'Measured','Null'});
title(sprintf('Animal %s',animals{plot_animal}));


% Total
animal_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(vertcat(x{use}),rxn_bins,'normalization','pdf'), ...
    use_days,{bhv.stim_move_t},'uni',false)');

animal_alt_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(cell2mat(vertcat(x{use})),rxn_bins,'normalization','pdf'), ...
    use_days,{bhv.alt_stim_move_t},'uni',false)');

figure; hold on;
histogram('BinEdges',rxn_bins,'BinCounts',nanmean(animal_rxn,1),'EdgeColor','none')
histogram('BinEdges',rxn_bins,'BinCounts',nanmean(animal_alt_rxn,1),'EdgeColor','none')
xlabel('Reaction time');
ylabel('PDF');
legend({'Measured','Null (all)'});

% Total (trial/alt number-matched)
alt_rxn_matched = ...
    cellfun(@(x) cellfun(@(x) cell2mat(cellfun(@(x) ...
    datasample(x,min(length(x),1)),x,'uni',false)),x,'uni',false), ...
    {bhv.alt_stim_move_t},'uni',false);

animal_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(vertcat(x{use}),rxn_bins,'normalization','pdf'), ...
    cellfun(@(x) x&([1:length(x)]' > 6),use_days,'uni',false),{bhv.stim_move_t},'uni',false)');

animal_alt_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(vertcat(x{use}),rxn_bins,'normalization','pdf'), ...
    cellfun(@(x) x&([1:length(x)]' > 6),use_days,'uni',false),alt_rxn_matched,'uni',false)');

figure;
subplot(2,1,1); hold on;
histogram('BinEdges',rxn_bins,'BinCounts',nanmean(animal_rxn,1),'EdgeColor','none')
histogram('BinEdges',rxn_bins,'BinCounts',nanmean(animal_alt_rxn,1),'EdgeColor','none')
xlabel('Reaction time');
ylabel('PDF');
legend({'Measured','Null (n trial matched)'});
subplot(2,1,2);
plot(rxn_bin_centers,nanmean(animal_rxn - animal_alt_rxn,1));
line(xlim,[0,0],'color','r');
ylabel('Distribution difference');
linkaxes(get(gcf,'children'),'xy');


% %%%%% TESTING STATS rxn within window
% THIS IS GOOD, USE THIS
max_days = max(cellfun(@sum,use_days));
r_frac_p_cat = nan(max_days,length(animals));
for curr_animal = 1:length(animals)
    for curr_day = find(use_days{curr_animal})'
        
        n_sample = 10000;
        use_trials = ~cellfun(@isempty,bhv(curr_animal).alt_stim_move_t{curr_day});
%         use_trials = ~cellfun(@isempty,bhv(curr_animal).alt_stim_move_t{curr_day}) & ...
%             linspace(0,1,length(bhv(curr_animal).stim_move_t{curr_day}))' < 0.9;
%           use_trials = ~cellfun(@isempty,bhv(curr_animal).alt_stim_move_t{curr_day}) & ...
%             (length(bhv(curr_animal).stim_move_t{curr_day}):-1:1)' > 10;
        
        r = bhv(curr_animal).stim_move_t{curr_day}(use_trials);
        ar = cell2mat(cellfun(@(x) datasample(x,n_sample)', ...
            bhv(curr_animal).alt_stim_move_t{curr_day}(use_trials),'uni',false));
        
        rxn_window = [0.1,0.2];
        
        r_frac = nanmean(r >= rxn_window(1) & r <= rxn_window(2));
        ar_frac = nanmean(ar >= rxn_window(1) & ar <= rxn_window(2),1);
        
        r_frac_rank = tiedrank([r_frac,ar_frac]);
        r_frac_p = r_frac_rank(1)./(n_sample+1);
        
        r_frac_p_cat(curr_day,curr_animal) = r_frac_p;
        
    end
end
figure; 
imagesc(r_frac_p_cat > 0.99,'AlphaData',~isnan(r_frac_p_cat));
set(gca,'Color',[0.5,0.5,0.5])
title('rxn window')
xlabel('Animal');
ylabel('Day');
[~,learned_day] = max(r_frac_p_cat,[],1);


% %%% TESTING STATS: interquartile range
% % (less good, and needed to exclude last 10% trials to get reasonable)
% max_days = max(cellfun(@sum,use_days));
% p = nan(max_days,length(animals));
% for curr_animal = 1:length(animals)
%     for curr_day = find(use_days{curr_animal})'
%         
%         n_sample = 10000;
% %         use_trials = ~cellfun(@isempty,bhv(curr_animal).alt_stim_move_t{curr_day});
%         use_trials = ~cellfun(@isempty,bhv(curr_animal).alt_stim_move_t{curr_day}) & ...
%             linspace(0,1,length(bhv(curr_animal).stim_move_t{curr_day}))' < 0.9;
% %         use_trials = ~cellfun(@isempty,bhv(curr_animal).alt_stim_move_t{curr_day}) & ...
% %             (length(bhv(curr_animal).stim_move_t{curr_day}):-1:1)' > 10;
%         
%         r = bhv(curr_animal).stim_move_t{curr_day}(use_trials);
%         ar = cell2mat(cellfun(@(x) datasample(x,n_sample)', ...
%             bhv(curr_animal).alt_stim_move_t{curr_day}(use_trials),'uni',false));
%         
% %         figure;
% %         subplot(1,2,1);hold on;
% %         plot(r,'.k');plot(ar(:,1),'.r');
% %         subplot(1,2,2);hold on;
% %         plot(sort(ar(:,1:100:end)),linspace(0,1,length(r)),'r');
% %         plot(sort(r),linspace(0,1,length(r)),'k','linewidth',2);
% %         line(xlim,[0.25,0.25]);
% %         line(xlim,[0.75,0.75]);
%         
%         r_frac_rank = tiedrank(iqr([r,ar],1));
% %         r_frac_rank = tiedrank([diff(prctile([r,ar],[20,80]))]);
% %         r_frac_rank = tiedrank([nanmean(diff(prctile([r,ar],[25,50,75])))]);
%         
%         r_frac_p = r_frac_rank(1)./(n_sample+1);
% 
%         p(curr_day,curr_animal) = r_frac_p;
%         
%     end
% end
% figure;
% imagesc(p < 0.01,'AlphaData',~isnan(p));
% set(gca,'Color',[0.5,0.5,0.5])
% xlabel('Animal');
% ylabel('Day');
% title('iqr');


% %%% TESTING KS STAT (multisample)
% max_days = max(cellfun(@sum,use_days));
% p = nan(length(animals),max_days);
% for curr_animal = 1:length(animals)
%     for curr_day = find(use_days{curr_animal})'
%         
%         use_rxn = ~cellfun(@isempty, ...
%             bhv(curr_animal).alt_stim_move_t{curr_day});
%         
%         n_sample = 1000;
%         curr_rxn = bhv(curr_animal).stim_move_t{curr_day}(use_rxn);
%         curr_alt_rxn = cell2mat(cellfun(@(x) ...
%             datasample(x,n_sample*2)', ...
%             bhv(curr_animal).alt_stim_move_t{curr_day}(use_rxn),'uni',false));
%         
%         curr_ksstat_null = nan(n_sample,1);
%         curr_ksstat_test = nan(n_sample,1);
%         for curr_sample = 1:n_sample
%             [~,~,curr_ksstat_null(curr_sample)] = ...
%                 kstest2(curr_alt_rxn(:,curr_sample), ...
%                 curr_alt_rxn(:,n_sample+curr_sample));
%             [~,~,curr_ksstat_test(curr_sample)] = ...
%                 kstest2(curr_alt_rxn(:,curr_sample), ...
%                 curr_rxn);
%         end
%         
%         curr_ksstat_rank = tiedrank([nanmean(curr_ksstat_test);curr_ksstat_null]);
%         curr_ksstat_p = curr_ksstat_rank(1)./(n_sample+1);
%         
%         p(curr_animal,curr_day) = curr_ksstat_p;
%         
%     end
% end
% figure;imagesc(p > 0.95);
% title('ks stat')

% Median reaction time (padded)
rxn_median_pad = cell2mat(cellfun(@(x) padarray(cellfun(@(x) ...
    nanmedian(x),x), ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_move_t},use_days,'uni',false),'uni',false));

% Fraction rxn within boundary
max_days = max(cellfun(@sum,use_days));

rxn_frac_window = [0.1,0.26];
rxn_frac_pad = cell2mat(cellfun(@(x) padarray(cellfun(@(x) ...
    nanmean(x > rxn_frac_window(1) & x < rxn_frac_window(2)),x), ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_move_t},use_days,'uni',false),'uni',false));

alt_rxn_frac_pad = cell2mat(cellfun(@(x) padarray(cellfun(@(x) ...
    nanmean(cell2mat(x) > rxn_frac_window(1) & cell2mat(x) < rxn_frac_window(2)),x), ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.alt_stim_move_t},use_days,'uni',false),'uni',false));

figure; hold on;
plot(rxn_frac_pad);
errorbar(nanmean(rxn_frac_pad,2),AP_sem(rxn_frac_pad,2),'k','linewidth',2);
errorbar(nanmean(alt_rxn_frac_pad,2),AP_sem(alt_rxn_frac_pad,2),'r','linewidth',2);

% Day-split reaction time fraction

n_daysplit = 4;

% Index: within-day split, day, animal
trial_split_idx = cellfun(@(x,use_days,animal_num) ...
    cellfun(@(x,day) ...
    [min(floor(linspace(1,n_daysplit+1,length(x))),n_daysplit)', ...
    repmat(day,length(x),1),repmat(animal_num,length(x),1)], ...
    x(use_days),num2cell(1:sum(use_days))','uni',false), ...
    {bhv.stim_move_t},use_days,num2cell(1:length(bhv)),'uni',false);

stim_move_t_cat = cell2mat(cellfun(@(x,use_days) cell2mat(x(use_days)), ...
    {bhv.stim_move_t},use_days,'uni',false)');
% alt_stim_move_t_cat = cell2mat(cellfun(@(x,use_days) cell2mat(x(use_days)), ...
%     alt_rxn_matched,use_days,'uni',false)');
alt_stim_move_t_cat = cell2mat(cellfun(@(x,use_days) ...
    cell2mat(cellfun(@(x) cellfun(@nanmedian,x),x(use_days),'uni',false)), ...
    {bhv.alt_stim_move_t},use_days,'uni',false)');

stim_move_t_stimtime = accumarray(cell2mat(cat(1,trial_split_idx{:})), ...
    stim_move_t_cat > rxn_frac_window(1) & ...
    stim_move_t_cat < rxn_frac_window(2)', ...
    [n_daysplit,max_days,length(bhv)],@nanmean,NaN);
alt_stim_move_t_stimtime = accumarray(cell2mat(cat(1,trial_split_idx{:})), ...
    alt_stim_move_t_cat > rxn_frac_window(1) & ...
    alt_stim_move_t_cat < rxn_frac_window(2)', ...
    [n_daysplit,max_days,length(bhv)],@nanmean,NaN);

figure; hold on;
stim_move_t_stimtime_long = reshape(padarray(stim_move_t_stimtime,[1,0,0],NaN,'post'),[],length(animals));
alt_stim_move_t_stimtime_long = reshape(padarray(alt_stim_move_t_stimtime,[1,0,0],NaN,'post'),[],length(animals));

plot([1:size(stim_move_t_stimtime_long,1)]/(n_daysplit+1),stim_move_t_stimtime_long,'color',[0.5,0.5,0.5]);
p1 = errorbar([1:size(stim_move_t_stimtime_long,1)]/(n_daysplit+1), ...
    nanmean(stim_move_t_stimtime_long,2),AP_sem(stim_move_t_stimtime_long,2),'k','linewidth',2);

plot([1:size(alt_stim_move_t_stimtime_long,1)]/(n_daysplit+1),alt_stim_move_t_stimtime_long,'color',[1,0.5,0.5]);
p2 = errorbar([1:size(alt_stim_move_t_stimtime_long,1)]/(n_daysplit+1), ...
    nanmean(alt_stim_move_t_stimtime_long,2),AP_sem(alt_stim_move_t_stimtime_long,2),'r','linewidth',2);

xlabel('Day');
ylabel('Frac rxn 100-250 ms');
legend([p1,p2],{'Measured','Null'},'location','nw');



%% Non-muscimol behavior

% Load muscimol injection info
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
muscimol_fn = [data_path filesep 'muscimol.mat'];
load(muscimol_fn);

% Use days before muscimol
use_days = cell(size(bhv));
for curr_animal = 1:length(bhv)
    muscimol_animal_idx = ismember({muscimol.animal},bhv(curr_animal).animal);
    if ~any(muscimol_animal_idx)
        use_days{curr_animal} = true(length(bhv(curr_animal).day),1);
        continue
    end
    muscimol_day_idx = datenum(bhv(curr_animal).day) >= ...
        datenum(muscimol(muscimol_animal_idx).day(1));
    use_days{curr_animal} = ~muscimol_day_idx;
end

% (get max days for padding)
max_days = max(cellfun(@sum,use_days));

% Time to movement
stim_move_t_allpad = cell2mat(cellfun(@(x) padarray(cellfun(@nanmedian,x), ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_move_t},use_days,'uni',false),'uni',false));
stim_move_t_std_allpad = cell2mat(cellfun(@(x) padarray(cellfun(@nanstd,x), ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_move_t},use_days,'uni',false),'uni',false));

% Time to reward
stim_feedback_t_allpad = cell2mat(cellfun(@(x) padarray(cellfun(@nanmedian,x), ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_feedback_t},use_days,'uni',false),'uni',false));
stim_feedback_t_std_allpad = cell2mat(cellfun(@(x) padarray(cellfun(@nanstd,x), ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_feedback_t},use_days,'uni',false),'uni',false));

% Fraction hits
frac_hit_allpad = cell2mat(cellfun(@(x) padarray(x, ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.frac_hit},use_days,'uni',false),'uni',false));

% Stim vs null response index
stim_response_idx_allpad = cell2mat(cellfun(@(x) padarray(x, ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_response_idx},use_days,'uni',false),'uni',false));

% Movement around stim
stim_surround_wheel_avg = cell2mat(cellfun(@(x) ...
    padarray(cell2mat(cellfun(@(x) nanmean(x,1),x,'uni',false)), ...
    [max_days-length(x),0],NaN,'post'), ...
    permute( ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_surround_wheel},use_days,'uni',false), ...
    [1,3,2]),'uni',false));

% Movement probability pre/post stim
stim_surround_t = bhv(1).stim_surround_t;
move_prestim_max = permute(max(stim_surround_wheel_avg(:,stim_surround_t<0,:),[],2),[3,1,2]);
move_poststim_max = permute(max(stim_surround_wheel_avg(:,stim_surround_t>=0,:),[],2),[3,1,2]);
move_prepost_max_ratio = ...
    (move_poststim_max-move_prestim_max)./(move_poststim_max+move_prestim_max);

%%% (unused)
% Trials per minute
trials_per_min = cell2mat(cellfun(@(x) padarray(x,[max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,y,use_days) x(use_days)./y(use_days), ....
    {bhv.n_trials},{bhv.session_duration},use_days,'uni',false),'uni',false));

% Total ITI (initial + resets)
iti_allpad = cell2mat(cellfun(@(x) padarray(cellfun(@nanmedian,x), ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.iti_t},use_days,'uni',false),'uni',false));

% Plot
figure;

subplot(2,2,1);
hold on;
errorbar(nanmean(stim_move_t_allpad,2),AP_sem(stim_move_t_allpad,2),'k','linewidth',2);
errorbar(nanmean(stim_feedback_t_allpad,2),AP_sem(stim_feedback_t_allpad,2),'b','linewidth',2);
ylabel('Time (s)');
xlabel('Session');
legend({'Stim to move','Stim to reward'});
set(gca,'YScale','log');
set(gca,'YTick',[0.25,0.5,1,2,4]);

subplot(2,2,2);
hold on;
stim_surround_t = bhv(1).stim_surround_t;
set(gca,'ColorOrder',copper(size(stim_surround_wheel_avg,1)));
plot(stim_surround_t,nanmean(stim_surround_wheel_avg,3)','linewidth',2);
xlabel('Time from stim (s)');
ylabel('Fraction moving');

subplot(2,2,3);
hold on;
plot(move_prepost_max_ratio');
errorbar(nanmean(move_prepost_max_ratio,1),AP_sem(move_prepost_max_ratio,1),'k','linewidth',2);
line(xlim,[0,0],'linestyle','--','color','k');
ylabel('Pre/post stim move ratio');
xlabel('Session');

subplot(2,2,4);
hold on;
plot(stim_response_idx_allpad);
errorbar(nanmean(stim_response_idx_allpad,2),AP_sem(stim_response_idx_allpad,2),'k','linewidth',2);
line(xlim,[0,0],'linestyle','--','color','k');
ylabel('Stim response index');
xlabel('Session');


% Stim response index daysplit
stim_move_max_plotcat = reshape(permute(padarray(cell2mat(permute(cellfun(@(x,use_days) ...
    padarray(vertcat(x{use_days}),[max_days-sum(use_days),0],NaN,'post'), ...
    {bhv.stim_move_max_daysplit},use_days,'uni',false),[1,3,2])),[0,1],NaN,'post'),[2,1,3]),[],length(bhv));

alt_stim_move_max_plotcat = reshape(permute(padarray(cell2mat(permute(cellfun(@(x,use_days) ...
    padarray(vertcat(x{use_days}),[max_days-sum(use_days),0],NaN,'post'), ...
    {bhv.alt_stim_move_max_daysplit},use_days,'uni',false),[1,3,2])),[0,1],NaN,'post'),[2,1,3]),[],length(bhv));

stim_response_idx_plotcat = reshape(permute(padarray(cell2mat(permute(cellfun(@(x,use_days) ...
    padarray(vertcat(x{use_days}),[max_days-sum(use_days),0],NaN,'post'), ...
    {bhv.stim_response_idx_daysplit},use_days,'uni',false),[1,3,2])),[0,1],NaN,'post'),[2,1,3]),[],length(bhv));

figure; hold on;
n_daysplit = length(bhv(1).stim_response_idx_daysplit{1});
daysplit_x = (1:(n_daysplit+1)*max_days)/(n_daysplit+1);
errorbar(daysplit_x,nanmean(stim_move_max_plotcat,2), ...
    AP_sem(stim_move_max_plotcat,2),'r','linewidth',2);
errorbar(daysplit_x,nanmean(alt_stim_move_max_plotcat,2), ...
    AP_sem(alt_stim_move_max_plotcat,2),'b','linewidth',2);
errorbar(daysplit_x,nanmean(stim_response_idx_plotcat,2), ...
    AP_sem(stim_response_idx_plotcat,2),'k','linewidth',2);
ylabel('Stim response index');
xlabel('Day');
legend({'Stim move max','Null move max','Stim response idx'},'location','se');


%% Muscimol behavior
% (grouped by injection location)

% Load muscimol injection info
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
muscimol_fn = [data_path filesep 'muscimol.mat'];
load(muscimol_fn);

muscimol_animals = find(ismember({muscimol.animal},animals));
muscimol_areas = flipud(unique(vertcat(muscimol(muscimol_animals).area)));

stim_surround_t = bhv(1).stim_surround_t;
muscimol_stim_surround_wheel = repmat( ...
    {nan(length(animals),length(stim_surround_t))}, ...
    length(muscimol_areas),1);

stim_move_t_muscimol_all = cell(length(muscimol_areas),length(animals));
alt_stim_move_t_muscimol_all = cell(length(muscimol_areas),length(animals));

stim_move_t_muscimol = nan(length(muscimol_areas),length(animals));
stim_feedback_t_muscimol = nan(length(muscimol_areas),length(animals));
move_feedback_t_muscimol = nan(length(muscimol_areas),length(animals));

rxn_frac_muscimol = nan(length(muscimol_areas),length(animals));
stim_response_idx_muscimol = nan(length(muscimol_areas),length(animals));

wheel_velocity_muscimol = nan(length(muscimol_areas),length(animals));
wheel_bias_muscimol = nan(length(muscimol_areas),length(animals));

figure;
for curr_area_idx = 1:length(muscimol_areas)
    curr_area = muscimol_areas{curr_area_idx};
    for curr_muscimol_animal_idx = 1:length(muscimol_animals)
        
        % Get index in bhv to match animal in muscimol
        curr_bhv_animal_idx = find(strcmp( ...
            muscimol(muscimol_animals(curr_muscimol_animal_idx)).animal, ...
            {bhv.animal}));
        
        % Get days of current animal with muscimol in current area
        muscimol_day_idx = find(strcmp(curr_area, ...
            muscimol(muscimol_animals(curr_muscimol_animal_idx)).area));
        
        % (if this mouse didn't have this experiment)
        if isempty(muscimol_day_idx)
            continue
        end
        
        % If muscimol, use last available day (doubled once when didn't work)
        if ~strcmp(curr_area,'washout')
            use_day = muscimol(muscimol_animals(curr_muscimol_animal_idx)).day{muscimol_day_idx(end)};
        else
            use_day = muscimol(muscimol_animals(curr_muscimol_animal_idx)).day(muscimol_day_idx);
        end
        bhv_day_idx = ismember(bhv(curr_bhv_animal_idx).day,use_day);
        
        curr_stim_surround_wheel = nanmean(cell2mat(cellfun(@(x) nanmean(x,1), ...
            bhv(curr_bhv_animal_idx).stim_surround_wheel(bhv_day_idx),'uni',false)),1);
        
        muscimol_stim_surround_wheel{curr_area_idx}(curr_muscimol_animal_idx,:) = ...
            curr_stim_surround_wheel;
        
        stim_move_t_muscimol_all{curr_area_idx,curr_muscimol_animal_idx} = ...
            cat(1,bhv(curr_bhv_animal_idx).stim_move_t{bhv_day_idx});
        alt_stim_move_t_muscimol_all{curr_area_idx,curr_muscimol_animal_idx} = ...
            cell2mat(cat(1,bhv(curr_bhv_animal_idx).alt_stim_move_t{bhv_day_idx}));
        
        stim_move_t_muscimol(curr_area_idx,curr_muscimol_animal_idx) = ...
            nanmedian(cat(1,bhv(curr_bhv_animal_idx).stim_move_t{bhv_day_idx}));
        stim_feedback_t_muscimol(curr_area_idx,curr_muscimol_animal_idx) = ...
            nanmedian(cat(1,bhv(curr_bhv_animal_idx).stim_feedback_t{bhv_day_idx}));
        move_feedback_t_muscimol(curr_area_idx,curr_muscimol_animal_idx) = ...
            nanmedian(cat(1,bhv(curr_bhv_animal_idx).stim_feedback_t{bhv_day_idx}) - ...
            cat(1,bhv(curr_bhv_animal_idx).stim_move_t{bhv_day_idx}));
        
        rxn_frac_muscimol(curr_area_idx,curr_muscimol_animal_idx) = ...
            nanmedian(cellfun(@(x) nanmean(x > 0.1 & x < 0.25), ...
            bhv(curr_bhv_animal_idx).stim_move_t(bhv_day_idx)));
        
        stim_response_idx_muscimol(curr_area_idx,curr_muscimol_animal_idx) = ...
            nanmedian(cat(1,bhv(curr_bhv_animal_idx).stim_response_idx(bhv_day_idx)));
        
        wheel_velocity_muscimol(curr_area_idx,curr_muscimol_animal_idx) = ...
            nanmean(bhv(curr_bhv_animal_idx).wheel_velocity(bhv_day_idx)./ ...
            bhv(curr_bhv_animal_idx).session_duration(bhv_day_idx));
        wheel_bias_muscimol(curr_area_idx,curr_muscimol_animal_idx) = ...
            nanmean(bhv(curr_bhv_animal_idx).wheel_bias(bhv_day_idx));
        
    end
    
    % Plot for each animal and muscimol area:
    % stim-surround movement
    subplot(length(muscimol_areas)+1,2,(curr_area_idx-1)*2+1); hold on;
    plot(stim_surround_t,muscimol_stim_surround_wheel{curr_area_idx}');
    plot(stim_surround_t,nanmean(muscimol_stim_surround_wheel{curr_area_idx},1),'k','linewidth',2);
    ylim([0,1]);
    title(curr_area);
    % reaction-time distribution
    subplot(length(muscimol_areas)+1,2,(curr_area_idx-1)*2+2); hold on;
    rxn_bins = [-Inf,0:0.01:1,Inf];
    cellfun(@(x) histogram(x,rxn_bins,'EdgeColor','none','normalization','pdf'), ...
        stim_move_t_muscimol_all(curr_area_idx,:))
    title(curr_area);
    
end

% Plot means overlaid
subplot(length(muscimol_areas)+1,2,length(muscimol_areas)*2+1); hold on;
plot(stim_surround_t, ...
    cell2mat(cellfun(@(x) nanmean(x,1),muscimol_stim_surround_wheel,'uni',false))', ...
    'linewidth',2);
legend(muscimol_areas);

subplot(length(muscimol_areas)+1,2,length(muscimol_areas)*2+2); hold on;
rxn_bins = [0:0.01:1];
rxn_bin_centers = rxn_bins(1:end-1)+(diff(rxn_bins)/2);
plot(rxn_bin_centers, ...
    nanmean(cell2mat(permute(cellfun(@(x) histcounts(x,rxn_bins,'normalization','probability'), ...
    stim_move_t_muscimol_all,'uni',false),[1,3,2])),3)','linewidth',2);
legend(muscimol_areas);

% Plot histograms for stim/null
rxn_bins = [0:0.01:1];
rxn_bin_centers = rxn_bins(1:end-1)+(diff(rxn_bins)/2);
stim_move_t_muscimol_hist = ...
    cell2mat(permute(cellfun(@(x) histcounts(x,rxn_bins,'normalization','probability'), ...
    stim_move_t_muscimol_all,'uni',false),[1,3,2]));
alt_stim_move_t_muscimol_hist = ...
    cell2mat(permute(cellfun(@(x) histcounts(x,rxn_bins,'normalization','probability'), ...
    alt_stim_move_t_muscimol_all,'uni',false),[1,3,2]));

figure;
for curr_area = 1:length(muscimol_areas)
    subplot(length(muscimol_areas),1,curr_area); hold on;
    histogram('BinEdges',rxn_bins,'BinCounts', ...
        squeeze(nanmean(stim_move_t_muscimol_hist(curr_area,:,:),3)),'EdgeColor','none')
    histogram('BinEdges',rxn_bins,'BinCounts', ...
        squeeze(nanmean(alt_stim_move_t_muscimol_hist(curr_area,:,:),3)),'EdgeColor','none')
    title(muscimol_areas{curr_area});
end
linkaxes(get(gcf,'Children'),'xy')

% Movement probability pre/post stim
stim_surround_t = bhv(1).stim_surround_t;
move_prestim_max = cellfun(@(x) max(x(:,stim_surround_t<0,:),[],2),muscimol_stim_surround_wheel,'uni',false);
move_poststim_max = cellfun(@(x) max(x(:,stim_surround_t>0,:),[],2),muscimol_stim_surround_wheel,'uni',false);
move_prepost_max_ratio = cell2mat(cellfun(@(x,y) (y-x)./(x+y),move_prestim_max,move_poststim_max,'uni',false)')';

% Plot average timings
figure;
subplot(1,7,1); hold on;
plot(stim_move_t_muscimol,'linewidth',2)
set(gca,'YScale','log')
set(gca,'XTick',1:length(muscimol_areas),'XTickLabel',muscimol_areas);
ylabel('Stim move t');

subplot(1,7,2); hold on;
plot(stim_feedback_t_muscimol,'linewidth',2)
set(gca,'YScale','log')
set(gca,'XTick',1:length(muscimol_areas),'XTickLabel',muscimol_areas);
ylabel('Stim reward t');

subplot(1,7,3); hold on;
plot(move_feedback_t_muscimol,'linewidth',2)
set(gca,'YScale','log')
set(gca,'XTick',1:length(muscimol_areas),'XTickLabel',muscimol_areas);
ylabel('Move feedback t');

subplot(1,7,4); hold on;
plot(wheel_velocity_muscimol,'linewidth',2)
set(gca,'YScale','log')
set(gca,'XTick',1:length(muscimol_areas),'XTickLabel',muscimol_areas);
ylabel('Wheel velocity/s');

subplot(1,7,5); hold on;
plot(move_prepost_max_ratio,'linewidth',2)
set(gca,'XTick',1:length(muscimol_areas),'XTickLabel',muscimol_areas);
ylabel('Post/post ratio');

subplot(1,7,6); hold on;
plot(rxn_frac_muscimol,'linewidth',2)
set(gca,'XTick',1:length(muscimol_areas),'XTickLabel',muscimol_areas);
ylabel('Rxn stim frac');

subplot(1,7,7); hold on;
plot(stim_response_idx_muscimol,'linewidth',2)
set(gca,'XTick',1:length(muscimol_areas),'XTickLabel',muscimol_areas);
ylabel('Stim response idx');


% Plot pre/post changes
figure;
subplot(1,4,1); hold on;
gscatter(move_prepost_max_ratio(:),wheel_velocity_muscimol(:), ...
    reshape(repmat(muscimol_areas,1,length(animals)),[],1),lines(length(muscimol_areas)),[],20);
xlabel('Pre/post move');
ylabel('Velocity/s');
axis square;

subplot(1,4,2); hold on;
washout_idx = strcmp('washout',muscimol_areas);
wheel_velocity_muscimol_norm = ...
    (wheel_velocity_muscimol-wheel_velocity_muscimol(washout_idx,:))./ ...
    (wheel_velocity_muscimol+wheel_velocity_muscimol(washout_idx,:));
move_prepost_max_ratio_diff = move_prepost_max_ratio - move_prepost_max_ratio(washout_idx,:);
gscatter(move_prepost_max_ratio_diff(:),wheel_velocity_muscimol_norm(:), ...
    reshape(repmat(muscimol_areas,1,length(animals)),[],1),lines(length(muscimol_areas)),[],20);
xlabel('Diff. Pre/post move');
ylabel('Norm. Velocity/s');
axis square;

subplot(1,4,3); hold on;
hold on;
washout_idx = strcmp('washout',muscimol_areas);
wheel_velocity_muscimol_norm = ...
    (wheel_velocity_muscimol-wheel_velocity_muscimol(washout_idx,:))./ ...
    (wheel_velocity_muscimol+wheel_velocity_muscimol(washout_idx,:));
stim_response_idx_muscimol_diff = stim_response_idx_muscimol - stim_response_idx_muscimol(washout_idx,:);
gscatter(stim_response_idx_muscimol_diff(:),wheel_velocity_muscimol_norm(:), ...
    reshape(repmat(muscimol_areas,1,length(animals)),[],1),lines(length(muscimol_areas)),[],20);
xlabel('Diff. stim response idx');
ylabel('Norm. Velocity/s');
axis square;

subplot(1,4,4); hold on;
hold on;
washout_idx = strcmp('washout',muscimol_areas);
rxn_frac_diff = rxn_frac_muscimol - rxn_frac_muscimol(washout_idx,:);
stim_response_idx_muscimol_diff = stim_response_idx_muscimol - stim_response_idx_muscimol(washout_idx,:);
gscatter(stim_response_idx_muscimol_diff(:),rxn_frac_diff(:), ...
    reshape(repmat(muscimol_areas,1,length(animals)),[],1),lines(length(muscimol_areas)),[],20);
xlabel('Diff. stim response idx');
ylabel('Diff. rxn frac');
axis square;


%% Learning within- or across-day: split trials within day

n_daysplit = 4;

% Load muscimol injection info
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
muscimol_fn = [data_path filesep 'muscimol.mat'];
load(muscimol_fn);

% Use days before muscimol
use_days = cell(size(bhv));
for curr_animal = 1:length(bhv)
    muscimol_animal_idx = ismember({muscimol.animal},bhv(curr_animal).animal);
    if ~any(muscimol_animal_idx)
        use_days{curr_animal} = true(length(bhv(curr_animal).day),1);
        continue
    end
    muscimol_day_idx = datenum(bhv(curr_animal).day) >= ...
        datenum(muscimol(muscimol_animal_idx).day(1));
    use_days{curr_animal} = ~muscimol_day_idx;
end

% (get max days for padding)
max_days = max(cellfun(@sum,use_days));

% Time to movement: bin, concat, plot
trial_split_idx = cellfun(@(x,use_days,animal_num) ...
    cellfun(@(x,day) ...
    [min(floor(linspace(1,n_daysplit+1,length(x))),n_daysplit)', ...
    repmat(day,length(x),1),repmat(animal_num,length(x),1)], ...
    x(use_days),num2cell(1:sum(use_days))','uni',false), ...
    {bhv.stim_move_t},use_days,num2cell(1:length(bhv)),'uni',false);

stim_move_t_cat = cell2mat(cellfun(@(x,use_days) cell2mat(x(use_days)), ...
    {bhv.stim_move_t},use_days,'uni',false)');
stim_move_t_daysplit = accumarray(cell2mat(cat(1,trial_split_idx{:})), ...
    stim_move_t_cat,[n_daysplit,max_days,length(bhv)],@nanmedian,NaN);

figure; hold on;
stim_move_t_daysplit_long = reshape(padarray(stim_move_t_daysplit,[1,0,0],NaN,'post'),[],length(animals));
plot([1:size(stim_move_t_daysplit_long,1)]/(n_daysplit+1),stim_move_t_daysplit_long);
set(gca,'YScale','log')
errorbar([1:size(stim_move_t_daysplit_long,1)]/(n_daysplit+1), ...
    nanmean(stim_move_t_daysplit_long,2),AP_sem(stim_move_t_daysplit_long,2),'k','linewidth',2);
set(gca,'YScale','log')
linkaxes(get(gcf,'Children'),'xy')
xlabel('Day');
ylabel('Stim to move (s)');

% Movement around stim
stim_surround_wheel_avg = cell2mat(cellfun(@(x) ...
    padarray(cell2mat(cellfun(@(x) nanmean(x,1),x,'uni',false)), ...
    [max_days-length(x),0],NaN,'post'), ...
    permute( ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_surround_wheel},use_days,'uni',false), ...
    [1,3,2]),'uni',false));

stim_surround_wheel_cat = cell2mat(cellfun(@(x,use_days) cell2mat(x(use_days)), ...
    {bhv.stim_surround_wheel},use_days,'uni',false)');

stim_surround_wheel_grp = accumarray(reshape(cat(3, ...
    repmat(1:size(stim_surround_wheel_cat,2),size(stim_surround_wheel_cat,1),1), ...
    repmat(permute(cell2mat(cat(1,trial_split_idx{:})),[1,3,2]),1,size(stim_surround_wheel_cat,2))),[],4), ...
    stim_surround_wheel_cat(:),[],@nanmean,NaN);

stim_surround_t = bhv(1).stim_surround_t;
stim_move_prestim_daysplit = squeeze(max(stim_surround_wheel_grp(stim_surround_t < 0,:,:,:),[],1));
stim_move_poststim_daysplit = squeeze(max(stim_surround_wheel_grp(stim_surround_t >= 0,:,:,:),[],1));
stim_move_ratio_daysplit = (stim_move_poststim_daysplit - stim_move_prestim_daysplit)./ ...
    (stim_move_poststim_daysplit + stim_move_prestim_daysplit);

figure; hold on
stim_move_ratio_daysplit_long = ...
    reshape(padarray(stim_move_ratio_daysplit,[1,0,0],NaN,'post'),[],length(animals));
plot([1:size(stim_move_ratio_daysplit_long,1)]/(n_daysplit+1),stim_move_ratio_daysplit_long);
errorbar([1:size(stim_move_ratio_daysplit_long,1)]/(n_daysplit+1), ...
    nanmean(stim_move_ratio_daysplit_long,2),AP_sem(stim_move_ratio_daysplit_long,2),'k','linewidth',2);
linkaxes(get(gcf,'Children'),'xy')
xlabel('Day');
ylabel('Stim response ratio');








