%% Notes
%
%
% Batch scripts to save preprocessed data here, saved to:
% C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning


%% ~~~~~~~~~~~~~~ Behavior

%% Grab and save behavior

animal_group = 'teto';
animals = { ...
    'AP100','AP101','AP103','AP104','AP105', ...
    'AP106','AP107','AP108','AP109','AP111', ...
    'AP113','AP114','AP115'};

protocol = 'AP_stimWheelRight';
flexible_name = false;
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
            
            % Time of session (in minutes)curr_day
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

            % Time between prior move stop and post-stim move start
            prestim_move_pause = ...
                wheel_starts(wheel_move_stim_idx(2:end)) - ...
                wheel_stops(wheel_move_stim_idx(2:end)-1);

            % Time between all movements
            all_move_pause = wheel_starts(2:end)-wheel_stops(1:end-1);

            % Non-stim movements per second
            wheel_iti_move_rate = arrayfun(@(x) ...
                sum(wheel_starts > signals_events.responseTimes(x) & ...
                wheel_starts < stimOn_times(x+1))./ ...
                (stimOn_times(x+1) - signals_events.responseTimes(x)), ...
                1:n_trials-1);


            wheel_iti_move_rate = arrayfun(@(x) ...
                sum(wheel_starts > signals_events.responseTimes(x) & ...
                wheel_starts < stimOn_times(x+1))./ ...
                (stimOn_times(x+1) - signals_events.responseTimes(x)), ...
                1:n_trials-1);

            wheel_iti_move_rate = arrayfun(@(x) ...
                sum(wheel_starts > signals_events.responseTimes(x) & ...
                wheel_starts < stimOn_times(x+1))./ ...
                (stimOn_times(x+1) - signals_events.responseTimes(x)), ...
                1:n_trials-1);

            nonresponse_move_rate = ...
                (length(wheel_starts)-length(wheel_move_response_idx))./ ...
                block.duration;

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
            % (general timings)
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

            bhv(curr_animal).prestim_move_pause{curr_day,1} = prestim_move_pause;
            bhv(curr_animal).all_move_pause{curr_day,1} = all_move_pause;
            bhv(curr_animal).wheel_iti_move_rate{curr_day,1} = wheel_iti_move_rate;
            bhv(curr_animal).nonresponse_move_rate(curr_day,1) = nonresponse_move_rate;

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


%% Find "learned" days
% Days when response to stimulus was significantly different from null

% Load behavior from above
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
bhv_fn = [data_path filesep 'bhv_teto'];
load(bhv_fn);

% Define "learned" days by median reaction time
% (exclude rxn < 0.1: too fast for stim, doesn't change with learning)
% (keeping window defined, but not used for stats)
rxn_window = [0.1,0.2];
n_sample = 10000;
learned_days = cell(size(bhv));
for curr_animal = 1:length(bhv)

    curr_n_days = length(bhv(curr_animal).stim_move_t);
    learned_days{curr_animal} = nan(curr_n_days,1);

    for curr_day = 1:curr_n_days

        % (use trials only with alt - no alt if couldn't recapture stim
        % time with params, happens with some skipped wheel clicks)
        use_trials = ~cellfun(@isempty,bhv(curr_animal).alt_stim_move_t{curr_day});
        
        curr_rxn = bhv(curr_animal).stim_move_t{curr_day}(use_trials);
        curr_alt_rxn = cell2mat(cellfun(@(x) datasample(x,n_sample)', ...
            bhv(curr_animal).alt_stim_move_t{curr_day}(use_trials),'uni',false));     
        
        % Get measured and null medians (exclude rxn < 0.1)
        rxn_stat = nanmedian(curr_rxn.*AP_nanout(curr_rxn < 0.1),1);
        alt_rxn_stat = nanmedian(curr_alt_rxn.*AP_nanout(curr_alt_rxn < 0.1),1);
        
        rxn_stat_rank = tiedrank([rxn_stat,alt_rxn_stat]);
        rxn_stat_p = rxn_stat_rank(1)./(n_sample+1);
        
        % (null rejected at 5%)
        rxn_window_frac_h = rxn_stat_p < 0.05;

        learned_days{curr_animal}(curr_day) = rxn_window_frac_h;
    end
end

% Plot to check
learned_days_padcat = AP_padcatcell(learned_days);

figure; 
imagesc(learned_days_padcat,'AlphaData',~isnan(learned_days_padcat));
set(gca,'Color',[0.5,0.5,0.5])
title('Sig. reaction time stat days');
set(gca,'XTick',1:length({bhv.animal}),'XTickLabel',{bhv.animal})
xlabel('Animal');
ylabel('Day');

% Put learned days into behavior structure and save
[bhv.learned_days] = learned_days{:};
[bhv.learned_days_rxn_window] = deal(rxn_window);
save(bhv_fn,'bhv');
disp(['Saved learned days ' bhv_fn]);

%% Find weekend/other gaps in training days
% (not used for anything - just to check)

% Load behavior
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
bhv_fn = [data_path filesep 'bhv_teto'];
load(bhv_fn);
animals = {bhv.animal};

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

% Get number of days from start
training_day_absolute = ...
    cellfun(@(day,use) datenum(day(use)) - datenum(day(1)), ...
    {bhv.day},use_days,'uni',false);


%% ~~~~~~~~~~~~~~ Widefield ROIs

%% Draw ROIs

% Set ROIs to draw
roi_areas = {'V1','V1binoc','AM','RSP','M2post','mPFC','SM'};

% Load reference image
% (this was manually made: RGB image of right stim pre/post and central)
wf_roi_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\wf_processing\wf_rois';
wf_ref_fn = 'wf_roi_ref.fig';
openfig([wf_roi_path filesep wf_ref_fn]);

wf_roi = struct('area',cell(length(roi_areas),2),'mask',cell(length(roi_areas),2));
for curr_area = 1:length(roi_areas)
    
    % Get ROI from left hemisphere
    title(['Draw ' roi_areas{curr_area} '_L']);
    curr_mask_L = roipoly;
    wf_roi(curr_area,1).area = [roi_areas{curr_area} '_L'];
    wf_roi(curr_area,1).mask = curr_mask_L;
    
    % Reflect ROI to right hemisphere
    curr_mask_R = AP_reflect_widefield(curr_mask_L) > 0;
    wf_roi(curr_area,2).area = [roi_areas{curr_area} '_R'];
    wf_roi(curr_area,2).mask = curr_mask_R;
    
    % Draw ROIs
    curr_roi_L = cell2mat(bwboundaries(curr_mask_L));
    curr_roi_R = cell2mat(bwboundaries(curr_mask_R));
    plot(curr_roi_L(:,2),curr_roi_L(:,1),'m','linewidth',2);
    plot(curr_roi_R(:,2),curr_roi_R(:,1),'m','linewidth',2);
    drawnow;
    
end

close(gcf);
% Save (commented so don't accidentally re-write)
% wf_roi_fn = [wf_roi_path filesep 'wf_roi'];
% save(wf_roi_fn,'wf_roi');
% disp('Saved new widefield ROIs');

%% Plot widefield ROIs

wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

roi_cat = cat(3,wf_roi.mask);
roi_col = [autumn(size(wf_roi,1));winter(size(wf_roi,1))];
% roi_col = repmat([0.5,0.7,0.5],n_rois,1);

figure; hold on
set(gca,'YDir','reverse');
AP_reference_outline('ccf_aligned',[0.7,0.7,0.7]);
for curr_roi = 1:n_rois
    curr_roi_boundary = cell2mat(bwboundaries(roi_cat(:,:,curr_roi)));
    patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),roi_col(curr_roi,:), ...
        'EdgeColor','none');
%     
%     text(nanmean(curr_roi_boundary(:,2)),nanmean(curr_roi_boundary(:,1)), ...
%         wf_roi(curr_roi).area,'FontSize',9,'HorizontalAlignment','center')
end
axis image off;

%% ~~~~~~~~~~~~~~ Facecam ROIs

%% Collect sample facecam frame from all experiments

% Overwrite, or add to existing data
overwrite_flag = false;

animals = { ...
    'AP100','AP101','AP103','AP104','AP105', ...
    'AP106','AP107','AP108','AP109','AP111', ...
    'AP113','AP114','AP115'};

% Init structure to overwrite, or load 
if overwrite_flag
    facecam_align = struct;
else
    facecam_processing_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\facecam_processing';
    facecam_align_fn = fullfile(facecam_processing_path,'facecam_align.mat');
    load(facecam_align_fn);
end

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);

    % If adding rather than overwriting, get days without data
    if ~overwrite_flag && length(facecam_align) >= curr_animal
        load_days = find(~ismember({experiments.day}, ...
            [facecam_align(curr_animal).day(cellfun(@ischar,[facecam_align(curr_animal).day]))]));
        if ~any(load_days)
            continue
        end
    else
        load_days = find(true(1,length(experiments)));
    end
    
    disp(['Loading ' animal]);
    
    for curr_day = load_days
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment (only cam)
        load_parts.cam = true;
        AP_load_experiment       

        if ~facecam_exists
            continue
        end

        % Grab sample frame (last frame)
        vr = VideoReader(facecam_fn);
        grab_frame = vr.NumFrames;
        facecam_sample_frame = read(vr,grab_frame);
        
        facecam_align(curr_animal).animal = animal;
        facecam_align(curr_animal).day{curr_day} = day;
        facecam_align(curr_animal).im{curr_day} = facecam_sample_frame;

        % Prep for next loop
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars('-except',preload_vars{:});

    end
end
disp('Done loading all');

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\facecam_processing';
save_fn = 'facecam_align';
save([save_path filesep save_fn],'facecam_align','-v7.3');
disp(['Saved: ' save_path filesep save_fn])


%% Align facecam frames (nose/eye control point)

% Overwrite, or add to existing data
overwrite_flag = false;

% Load facecam sample frames
facecam_processing_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\facecam_processing';
facecam_align_fn = fullfile(facecam_processing_path,'facecam_align.mat');
load(facecam_align_fn);

% Get days without tform to align (or all if overwriting)
if ~overwrite_flag
    im_cat = {facecam_align.im};
    tform_cat = {facecam_align.tform};
    tform_cat_pad = cellfun(@(im,tform) ...
        [tform;cell(length(im)-length(tform),1)],im_cat,tform_cat,'uni',false);
    align_days = find(cellfun(@(im,tform) ~isempty(im) & isempty(tform), ...
        horzcat(im_cat{:}),vertcat(tform_cat_pad{:})'));
else
    align_days = 1:length(length([facecam_align.im]));
end

% Plot images and select control points
im_unaligned = cellfun(@double,[facecam_align.im],'uni',false);

target_im = im_unaligned{1};
target_size = size(target_im);

figure;
target_ax = subplot(1,2,1);
imagesc(target_im);
axis image off; hold on;
source_ax = subplot(1,2,2);
source_h = imagesc([]);
axis image off; hold on;

title(target_ax,'Click: nose end, eye center');
target_ctrl_points = ginput(2);
plot(target_ax,target_ctrl_points(:,1),target_ctrl_points(:,2),'.r','MarkerSize',20);
title(target_ax,'');

im_aligned = nan(target_size(1),target_size(2),length(im_unaligned));
source_ctrl_points = cell(length(im_unaligned),1);
if overwrite_flag
    cam_tform = cell(length(im_unaligned),1);
else
    cam_tform = vertcat(tform_cat_pad{:});
end
for curr_im = align_days
    source_im = im_unaligned{curr_im};
    if isempty(source_im)
        continue
    end

    % Click control points
    title(source_ax,{['Click: nose eye'], ...
        [sprintf('%d/%d',curr_im,length(im_unaligned))]});
    set(source_h,'CData',source_im);
    source_ctrl_points{curr_im} = ginput(2);

    % Store tform
    cam_tform{curr_im} = fitgeotrans(source_ctrl_points{curr_im}, ...
        target_ctrl_points,'nonreflectivesimilarity');
    tform_size = imref2d(target_size);
    im_aligned(:,:,curr_im) = ...
        imwarp(source_im,cam_tform{curr_im},'OutputView',tform_size);

end
close(gcf);

% Plot aligned
AP_imscroll(im_aligned); axis image

% Save transform (into original struct)
% (package back into animals)
n_days_animal = cellfun(@length,{facecam_align.day});
cam_tform_animal = mat2cell(cam_tform,n_days_animal);
[facecam_align.tform] = cam_tform_animal{:};

save(facecam_align_fn,'facecam_align');
disp(['Saved: ' facecam_align_fn]);


%% Load and align (manual check)

% Load facecam align
facecam_processing_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\facecam_processing';
facecam_align_fn = fullfile(facecam_processing_path,'facecam_align.mat');
load(facecam_align_fn);

% Align facecams
im_unaligned_cat = cellfun(@double,[facecam_align.im],'uni',false);
im_tform_cat = cat(1,facecam_align.tform);
im_ref = im_unaligned_cat{1};

im_aligned = nan(size(im_ref,1),size(im_ref,2),length(im_unaligned_cat));
for curr_im =1:length(im_unaligned_cat)
    if isempty(im_unaligned_cat{curr_im})
        continue
    end
    tform_size = imref2d(size(im_ref));
    im_aligned(:,:,curr_im) = ...
        imwarp(im_unaligned_cat{curr_im},im_tform_cat{curr_im}, ...
        'OutputView',tform_size);
end

day_label_split = cellfun(@(animal,day) cellfun(@(day) ...
    sprintf('%s %s',animal,day),day,'uni',false), ...
    {facecam_align.animal},{facecam_align.day},'uni',false);
day_label = horzcat(day_label_split{:});

im_unaligned_padcat = AP_padcatcell(im_unaligned_cat);
AP_imscroll([ ...
    im_unaligned_padcat(1:size(im_aligned,1),1:size(im_aligned,2),:), ...
    im_aligned],day_label); axis image;

% Plot single day all mice to check identity
figure;
h = tiledlayout('flow');
for x = 1:length(facecam_align)
    nexttile
    imagesc(facecam_align(x).im{2}); axis image off;
    title(sprintf('Animal %d: %s',x,facecam_align(x).animal));
end


%% Draw facecam whisker ROI and align to each recording

% Load facecam align
facecam_processing_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\facecam_processing';
facecam_align_fn = fullfile(facecam_processing_path,'facecam_align.mat');
load(facecam_align_fn);

% Align facecams
im_unaligned_cat = cellfun(@double,[facecam_align.im],'uni',false);
im_tform_cat = cat(1,facecam_align.tform);
im_ref = im_unaligned_cat{1};

im_aligned = nan(size(im_ref,1),size(im_ref,2),length(im_unaligned_cat));
for curr_im =1:length(im_unaligned_cat)
    if isempty(im_unaligned_cat{curr_im})
        continue
    end
    tform_size = imref2d(size(im_ref));
    im_aligned(:,:,curr_im) = ...
        imwarp(im_unaligned_cat{curr_im},im_tform_cat{curr_im}, ...
        'OutputView',tform_size);
end

% Plot average image, draw whisker ROI
figure;imagesc(nanmean(im_aligned,3));axis image off;
title('Draw whisker line ROI');
whisker_line = drawline;
master_whisker_mask = createMask(whisker_line);
close(gcf);

figure;
image(imoverlay(mat2gray(nanmean(im_aligned,3)),master_whisker_mask,'r'));
axis image off
title('Master whisker mask');

% Align whisker mask to individual day
whisker_mask = cell(size(im_unaligned_cat));
for curr_im =1:length(im_unaligned_cat)
    if isempty(im_unaligned_cat{curr_im})
        continue
    end

    tform_size = imref2d(size(im_unaligned_cat{curr_im}));
    whisker_mask{curr_im} = ...
        imwarp(master_whisker_mask,invert(im_tform_cat{curr_im}), ...
        'OutputView',tform_size); 
end

% Package into original structure
n_days_animal = cellfun(@length,{facecam_align.day});
whisker_mask_animal = mat2cell(whisker_mask,1,n_days_animal);
[facecam_align.whisker_mask] = whisker_mask_animal{:};

% Plot through all days to check the aligned whisker ROIs
whisker_mask_cat = [facecam_align.whisker_mask];
figure; h = image(im_unaligned_cat{1}); axis image off;
for curr_im =1:length(im_unaligned_cat)
    if isempty(im_unaligned_cat{curr_im})
        continue
    end
    set(h,'CData', ...
        imoverlay(mat2gray(im_unaligned_cat{curr_im}), ...
        whisker_mask_cat{curr_im},'r'));
    pause(0.1);
end

warning('Saftey: saving turned off');
% % Save whisker mask (into original struct)
% save(facecam_align_fn,'facecam_align');
% disp(['Saved: ' facecam_align_fn]);


%% ~~~~~~~~~~~~~~ Widefield data

%% Passive widefield

disp('Passive trial activity')

% Initialize
clear all
trial_data_all = struct;

animals = { ...
    'AP100','AP101','AP103','AP104','AP105', ...
    'AP106','AP107','AP108','AP109','AP111', ...
    'AP113','AP114','AP115'};

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    passive_experiments = AP_find_experiments(animal,protocol);
    
    % Get days after muscimol starts (to exclude)
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    if any(muscimol_animal_idx)
        muscimol_start_day = muscimol(muscimol_animal_idx).day{1};
        muscimol_experiments = datenum({passive_experiments.day})' >= datenum(muscimol_start_day);
    else
        muscimol_experiments = false(size({passive_experiments.day}));
    end

    % Get days with standard task (to exclude other tasks)
    task_protocol = 'AP_stimWheelRight';
    task_experiments = AP_find_experiments(animal,task_protocol);
    standard_task_experiments =  datenum({passive_experiments.day})' <= ...
        datenum(task_experiments(end).day);
    
    % Set experiments to use (imaging, not muscimol, standard task)
    experiments = passive_experiments([passive_experiments.imaging] & ...
        ~muscimol_experiments & standard_task_experiments);
    
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
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
    
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = fullfile(save_path,'trial_activity_passive_teto');
save(save_fn,'-v7.3');
disp(['Saved: ' save_fn])


%% Task widefield

disp('Task trial activity')

% Initialize
clear all
trial_data_all = struct;

animals = { ...
    'AP100','AP101','AP103','AP104','AP105', ...
    'AP106','AP107','AP108','AP109','AP111', ...
    'AP113','AP114','AP115'};

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_stimWheelRight';
    experiments = AP_find_experiments(animal,protocol);
    
    % Get days with muscimol
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    if any(muscimol_animal_idx)
        muscimol_start_day = muscimol(muscimol_animal_idx).day{1};
        muscimol_experiments = datenum({experiments.day})' >= datenum(muscimol_start_day);
    else
        muscimol_experiments = false(size({experiments.day}));
    end
    
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

        % (turned regression off)
%         trial_data_all.task_regressor_labels = task_regressor_labels;
%         trial_data_all.task_regressor_sample_shifts = task_regressor_sample_shifts;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end

    % Save (after each animal: matlab often crashes during this)
    save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    save_fn = fullfile(save_path,'trial_activity_task_teto');
    save(save_fn,'trial_data_all','-v7.3');
    disp(['Saved ' save_fn]);

end

disp('Finished loading all')


%% Passive long-term post-learning

disp('Long-term passive trial activity')

% Initialize
clear all
trial_data_all = struct;

animals = {'AP113','AP114','AP115'};

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    passive_experiments = AP_find_experiments(animal,protocol);    

    % Get days with task
    % (these mice were trained in right, then left task)
    task_right_protocol = 'AP_stimWheelRight';
    task_right_experiments = AP_find_experiments(animal,task_right_protocol);

    task_left_protocol = 'AP_stimWheelLeftReverse';
    task_left_experiments = AP_find_experiments(animal,task_left_protocol);

    % Set experiments to use:
    % - Last passive with left task
    % - Passive after task is finished
    last_task_left_experiment = datenum({passive_experiments.day}) == ....
        datenum(task_left_experiments(end).day);

    task_last_day = max(vertcat(datenum({task_right_experiments.day}),datenum({task_left_experiments.day})));
    post_training_experiments = datenum({passive_experiments.day}) > task_last_day;

    experiments = passive_experiments( ...
        last_task_left_experiment | ...
        post_training_experiments);

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
        
        % Store day relative to last task day (right and left-hand tasks)
        trial_data_all.last_task_right_day{curr_animal} = datenum(task_right_experiments(end).day);
        trial_data_all.last_task_left_day{curr_animal} = datenum(task_left_experiments(end).day);
        trial_data_all.recording_day{curr_animal}(curr_day) = datenum(day);

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
save_fn = fullfile(save_path,'trial_activity_passive_teto_postlearn');
save(save_fn,'-v7.3');
disp(['Saved: ' save_fn])




%% ~~~~~~~~~~~~~~ Muscimol + widefield data

%% Muscimol - passive

clear all
disp('Muscimol: Passive trial activity (tetO)')

animals = {'AP100','AP105','AP106','AP107','AP108'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);
    
    % Get days with V1 muscimol and washout
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    
    v1_muscimol_idx = find(strcmpi(muscimol(muscimol_animal_idx).area,'v1'));
    washout_idx = find(strcmpi(muscimol(muscimol_animal_idx).area,'washout'));
    v1_washout_idx = washout_idx(find(washout_idx > v1_muscimol_idx,1,'first'));
    
    % Set experiments to use (V1 muscimol and subsequent washout)
    use_muscimol_experiments = ismember({experiments.day}, ...
        muscimol(muscimol_animal_idx).day([v1_muscimol_idx,v1_washout_idx]));
    experiments = experiments(use_muscimol_experiments);
    
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


%% Muscimol - task

clear all
disp('Muscimol: Task trial activity (tetO)')

animals = {'AP100','AP105','AP106','AP107','AP108'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_stimWheelRight';
    experiments = AP_find_experiments(animal,protocol);
    
    % Get days with V1 muscimol and washout
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    
    v1_muscimol_idx = find(strcmpi(muscimol(muscimol_animal_idx).area,'v1'));
    washout_idx = find(strcmpi(muscimol(muscimol_animal_idx).area,'washout'));
    v1_washout_idx = washout_idx(find(washout_idx > v1_muscimol_idx,1,'first'));
    
    % Set experiments to use (V1 muscimol and subsequent washout)
    use_muscimol_experiments = ismember({experiments.day}, ...
        muscimol(muscimol_animal_idx).day([v1_muscimol_idx,v1_washout_idx]));
    experiments = experiments(use_muscimol_experiments);
    
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


%% ~~~~~~~~~~~~~~ Electrophysiology data

%% Passive ephys

clear all
disp('Passive trial activity (ephys)')

animals = {'AP100','AP101','AP104','AP105','AP106'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);
    
    % Set experiments to use (ephys)
    experiments = experiments([experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        % (load LFP to find cortex start)
        lfp_channel = 'all';
        ephys_align = 'cortex';
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
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
    
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = ['trial_activity_passive_ephys'];
save([save_path filesep save_fn],'-v7.3');
disp(['Saved: ' save_path filesep save_fn])


%% Task ephys

clear all
disp('Task trial activity (ephys)')

animals = {'AP100','AP101','AP104','AP105','AP106'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_stimWheelRight';
    experiments = AP_find_experiments(animal,protocol);
    
    % Set experiments to use (ephys)
    experiments = experiments([experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        % (load LFP to find cortex start)
        lfp_channel = 'all';
        ephys_align = 'cortex';
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

        % (not saving by depth anymore)
%         trial_data_all.mua_depth_edges = depth_group_edges;

        % (turned regression off)
%         trial_data_all.task_regressor_labels = task_regressor_labels;
%         trial_data_all.task_regressor_sample_shifts = task_regressor_sample_shifts;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = ['trial_activity_task_ephys'];
save([save_path filesep save_fn],'-v7.3');



%% Single-unit (task + passive)
% (save baseline + response for move (task) and stim (passive))

animals = {'AP100','AP101','AP104','AP105','AP106'};

% Initialize save variable
single_unit_data_all = struct;

for curr_animal = 1:length(animals)

    animal = animals{curr_animal};

    % Get passive ephys experiments
    passive_protocol = 'AP_lcrGratingPassive';
    passive_experiments = AP_find_experiments(animal,passive_protocol);
    passive_experiments = passive_experiments([passive_experiments.ephys]);

    % Get task ephys experiments
    task_protocol = 'AP_stimWheelRight';
    task_experiments = AP_find_experiments(animal,task_protocol);
    task_experiments = task_experiments([task_experiments.ephys]);

    % (make sure they're the same length and days)
    if length(passive_experiments) ~= length(task_experiments) || ...
            ~all([passive_experiments.day] == [task_experiments.day])
        error('Mismatched passive/task experiments')
    end

    for curr_day = 1:length(task_experiments)
        day = task_experiments(curr_day).day;
        for curr_exptype_cell = {'task','passive'}
            curr_exptype = cell2mat(curr_exptype_cell);

            switch curr_exptype
                case 'task'
                    experiment = task_experiments(curr_day).experiment(end);
                case 'passive'
                    experiment = passive_experiments(curr_day).experiment(end);
            end

            %%%% LOAD DATA, APPLY JF GOOD SINGLE UNITS

            % (load LFP to find cortex start)
            preload_vars = who;
            lfp_channel = 'all';
            ephys_align = 'cortex';
            ephys_quality_control = false;
            AP_load_experiment;

            % JF code for loading good single units from bombcell
            % unitType: 0 = noise, 1 = good, 2 = multiunit
            ephysDirPath = AP_cortexlab_filename(animal, day, experiment, 'ephys_dir');
            qMetric_fn = fullfile(ephysDirPath, 'qMetrics');
            load(fullfile(qMetric_fn, 'qMetric.mat'))
            load(fullfile(qMetric_fn, 'param.mat'))
            clearvars unitType;

            % DEFAULT CHANGE: eliminate amplitude cutoff
            % (for one recording it got rid of almost all cells, and
            % cells under amplitude cutoff still look good)
            param.minAmplitude = 0;

            % (classify good cells)
            unitType = nan(length(qMetric.percSpikesMissing), 1);
            unitType( ...
                qMetric.nPeaks > param.maxNPeaks | ...
                qMetric.nTroughs > param.maxNTroughs | ...
                qMetric.somatic ~= param.somatic | ...
                qMetric.spatialDecaySlope <=  param.minSpatialDecaySlope | ...
                qMetric.waveformDuration < param.minWvDuration |...
                qMetric.waveformDuration > param.maxWvDuration  | ...
                qMetric.waveformBaseline >= param.maxWvBaselineFraction) = 0;
            unitType( ...
                any(qMetric.percSpikesMissing <= param.maxPercSpikesMissing, 2)' & ...
                qMetric.nSpikes > param.minNumSpikes & ...
                any(qMetric.Fp <= param.maxRPVviolations, 2)' & ...
                qMetric.rawAmplitude > param.minAmplitude & isnan(unitType)') = 1;
            unitType(isnan(unitType)') = 2;

            % (some upwards waveforms not caught in .somatic? remove)
            upward_waveforms = max(waveforms,[],2) > abs(min(waveforms,[],2));

            % Templates already 1/re-indexed, grab good ones
            good_templates = unitType == 1 & ~upward_waveforms;
            good_templates_idx = find(good_templates);

            % Throw out all non-good template data
            templates = templates(good_templates,:,:);
            template_depths = template_depths(good_templates);
            waveforms = waveforms(good_templates,:);
            templateDuration = templateDuration(good_templates);
            templateDuration_us = templateDuration_us(good_templates);

            % Throw out all non-good spike data
            good_spike_idx = ismember(spike_templates,good_templates_idx);
            spike_times = spike_times(good_spike_idx);
            spike_templates_0idx = spike_templates_0idx(good_spike_idx);
            template_amplitudes = template_amplitudes(good_spike_idx);
            spike_depths = spike_depths(good_spike_idx);
            spike_times_timeline = spike_times_timeline(good_spike_idx);

            % Rename the spike templates according to the remaining templates
            % (and make 1-indexed from 0-indexed)
            new_spike_idx = nan(max(spike_templates_0idx)+1,1);
            new_spike_idx(unique(spike_templates_0idx)+1) = 1:length(unique(spike_templates_0idx));
            spike_templates = new_spike_idx(spike_templates_0idx+1);

            % Get raw waveforms from JF qMetrics
            % (convert to uV: raw * 0.195 (OE legacy) * 2.34 (uV/bit)
            waveforms_raw_all_full = cat(3,qMetric.rawWaveforms.spkMapMean);
            [~,waveforms_raw_all_maxchan] = max(max(abs(waveforms_raw_all_full),[],2),[],1);
            waveforms_raw = cell2mat(arrayfun(@(x) ...
                waveforms_raw_all_full(waveforms_raw_all_maxchan(x),:,x), ...
                good_templates_idx,'uni',false)).* 0.195 * 2.34;

            %%%% GET RESPONSES ALIGNED TO EVENT
            switch curr_exptype
                case 'task'
                    % Task: rewardable delay-period movement onsets
                    % (selection copied from operant_grab_trial_data.m)
                    wheel_moves_deg = arrayfun(@(x) wheel_position_deg( ...
                        Timeline.rawDAQTimestamps >= wheel_starts(x) & ...
                        Timeline.rawDAQTimestamps <= wheel_stops(x)) - ...
                        wheel_position_deg(find(Timeline.rawDAQTimestamps >= wheel_starts(x),1)), ...
                        1:length(wheel_starts),'uni',false);

                    deg_reward = -90;
                    deg_punish = 90;

                    wheel_moves_deg_rewardlimit = find(cellfun(@(x) ...
                        any(x <= deg_reward) && ...
                        ~any(x >= deg_punish),wheel_moves_deg));

                    move_rewardable = ...
                        intersect(wheel_move_nostim_idx',wheel_moves_deg_rewardlimit);

                    use_align = wheel_starts(move_rewardable);

                case 'passive'
                    % Passive: stim onsets (right-side & quiescent)
                    wheel_window = [0,0.5];
                    wheel_window_t = wheel_window(1):1/Timeline.hw.daqSampleRate:wheel_window(2);
                    wheel_window_t_peri_event = bsxfun(@plus,stimOn_times,wheel_window_t);
                    event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
                        +wheel_move,wheel_window_t_peri_event,'previous');
                    quiescent_trials = ~any(abs(event_aligned_wheel) > 0,2);

                    use_align = stimOn_times(quiescent_trials & stimIDs == 3);

            end

            % Get spike rate within baseline/response window
            baseline_window = [-0.5,-0.3];
            switch curr_exptype
                case 'task'
                    response_window = [-0.1,0.1];
                case 'passive'
                    response_window = [0,0.2];
            end

            event_fr = nan(length(use_align),size(templates,1),2);
            for curr_event = 1:length(use_align)

                curr_use_spikes_baseline = ~isnan(discretize(spike_times_timeline, ...
                    use_align(curr_event) + baseline_window));
                curr_use_spikes_response = ~isnan(discretize(spike_times_timeline, ...
                    use_align(curr_event) + response_window));

                curr_baseline_fr = accumarray(spike_templates(curr_use_spikes_baseline),1, ...
                    [size(templates,1),1],@sum)./diff(baseline_window);
                curr_response_fr = accumarray(spike_templates(curr_use_spikes_response),1, ...
                    [size(templates,1),1],@sum)./diff(response_window);

                event_fr(curr_event,:,1) = curr_baseline_fr;
                event_fr(curr_event,:,2) = curr_response_fr;

            end

            % Get PSTH for each unit
            raster_window = [-0.5,1];
            raster_sample_rate = 50;

            raster_sample_time = 1/raster_sample_rate;
            t = raster_window(1):raster_sample_time:raster_window(2);
            t_peri_stim = bsxfun(@plus,use_align,t);
            t_peri_stim_bins = [t_peri_stim-raster_sample_time/2,t_peri_stim(:,end)+raster_sample_time/2];

            unit_psth = nan(size(templates,1),length(t));
            for curr_unit = 1:size(templates,1)
                curr_spikes = spike_times_timeline(spike_templates == curr_unit);
                curr_psth = cell2mat(arrayfun(@(x) ...
                    histcounts(curr_spikes,t_peri_stim_bins(x,:)), ...
                    [1:size(t_peri_stim_bins,1)]','uni',false))*raster_sample_rate;
                unit_psth(curr_unit,:) = nanmean(curr_psth,1);
            end

            % Store rate, psth, and event index info
            switch curr_exptype
                case 'task'
                    single_unit_data_all(curr_animal).task_move_fr{curr_day} = ...
                        event_fr;
                    single_unit_data_all(curr_animal).task_move_psth{curr_day} = ...
                        unit_psth;

                case 'passive'
                    single_unit_data_all(curr_animal).passive_stim_fr{curr_day} = ...
                        event_fr;
                    single_unit_data_all(curr_animal).passive_stim_psth{curr_day} = ...
                        unit_psth;
            end

            % Get and store area of each unit
            probe_area_boundary_starts = cellfun(@(x) x(1),probe_area_boundaries);
            [~,area_sort] = sort(probe_area_boundary_starts);

            unit_area = probe_areas(area_sort(discretize(template_depths, ...
                [-Inf;probe_area_boundary_starts(area_sort(2:end));Inf])));

            single_unit_data_all(curr_animal).unit_area{curr_day} = unit_area;

            % Store waveform and spike width for each unit (both mine and Julie's)
            single_unit_data_all(curr_animal).waveforms{curr_day} = waveforms;
            single_unit_data_all(curr_animal).waveforms_raw{curr_day} = waveforms_raw;

            single_unit_data_all(curr_animal).waveform_duration_JF{curr_day} = ...
                qMetric.waveformDuration(good_templates)';
            single_unit_data_all(curr_animal).waveform_duration_AP{curr_day} = ...
                templateDuration_us;

            % Prep next loop
            clearvars('-except',preload_vars{:});
        end

         AP_print_progress_fraction(curr_day,length(task_experiments));

    end
end

disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = 'single_unit_data_all';
save([save_path filesep save_fn],'single_unit_data_all','-v7.3');
disp(['Saved: ' save_path filesep save_fn])


% (JUST STORING HERE - raster-making for debugging/sanity-check)
% raster_window = [-0.5,0.5];
% raster_sample_rate = 50;
% 
% raster_sample_time = 1/raster_sample_rate;
% t = raster_window(1):raster_sample_time:raster_window(2);
% 
% t_peri_event = bsxfun(@plus,use_align,t);
% t_peri_event_bins = [t_peri_event-raster_sample_time/2,t_peri_event(:,end)+raster_sample_time/2];
% 
% event_raster = nan(size(t_peri_event_bins,1),length(t),size(templates,1));
% for curr_unit = 1:size(templates,1)
%     curr_spikes = spike_times_timeline(spike_templates == curr_unit);
%     event_raster(:,:,curr_unit) = cell2mat(arrayfun(@(x) ...
%         histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
%         [1:size(t_peri_event_bins,1)]','uni',false))*raster_sample_rate;
% end
% AP_imscroll(event_raster(:,:,event_fr_diff_sig));
% AP_imscroll(event_raster(:,:,~event_fr_diff_sig));
% 
% [~,sig_sort] = sort(event_fr_diff(event_fr_diff_sig),'descend');
% [~,nonsig_sort] = sort(event_fr_diff(~event_fr_diff_sig),'descend');
% 
% event_raster_avg_sig = permute(nanmean(event_raster(:,:,event_fr_diff_sig),1),[3,2,1]);
% event_raster_avg_nonsig = permute(nanmean(event_raster(:,:,~event_fr_diff_sig),1),[3,2,1]);
% 
% figure;
% imagesc(t,[],[event_raster_avg_sig(sig_sort,:);event_raster_avg_nonsig(nonsig_sort,:)]);
% yline(length(sig_sort),'linewidth',3);
% xline(0);
% colormap(AP_colormap('BWR'))
% caxis(max(abs(caxis))*[-1,1]);
% xlabel('Time from event');
% ylabel('Unit');


%% Passive ephys (naive)

clear all
disp('Passive trial activity (ephys)')

animals = {'AP116','AP117','AP118','AP119'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);
    
    % Set experiments to use (ephys)
    experiments = experiments([experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        % (load LFP to find cortex start)
        lfp_channel = 'all';
        ephys_align = 'cortex';
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
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
    
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = ['trial_activity_passive_ephys_naive'];
save([save_path filesep save_fn],'-v7.3');
disp(['Saved: ' save_path filesep save_fn])
















