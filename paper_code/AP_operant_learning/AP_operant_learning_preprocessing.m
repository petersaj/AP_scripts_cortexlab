%% Notes
%
%
% Batch scripts to save preprocessed data here, saved to:
% C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning


%% ~~~~~~~~~~~~~~ Behavior

%% Grab and save behavior

animal_group = 'teto';
animals = {'AP100','AP101','AP103','AP104','AP105', ...
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

% Define "learned" days by fraction of reaction times within window
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
        
        rxn_window_frac = nanmean(curr_rxn >= rxn_window(1) & curr_rxn <= rxn_window(2));
        alt_rxn_window_frac = nanmean(curr_alt_rxn >= rxn_window(1) & curr_alt_rxn <= rxn_window(2),1);
        
        rxn_window_frac_rank = tiedrank([rxn_window_frac,alt_rxn_window_frac]);
        rxn_window_frac_p = rxn_window_frac_rank(1)./(n_sample+1);
        
        % (null rejected at 1%)
        rxn_window_frac_h = rxn_window_frac_p > 0.99;

        learned_days{curr_animal}(curr_day) = rxn_window_frac_h;
    end
end

% Plot to check
learned_days_padcat = AP_padcatcell(learned_days);

figure; 
imagesc(learned_days_padcat,'AlphaData',~isnan(learned_days_padcat));
set(gca,'Color','r')
title(sprintf('Frac reaction times %.2f-%.2f',rxn_window(1),rxn_window(2)));
xlabel('Animal');
ylabel('Day');

% Put learned days into behavior structure and save
[bhv.learned_days] = learned_days{:};
save(bhv_fn,'bhv');
disp(['Saved learned days ' bhv_fn]);


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

figure; hold on
set(gca,'YDir','reverse');
AP_reference_outline('ccf_aligned','k');%AP_reference_outline('retinotopy','m');
for curr_roi = 1:n_rois
    curr_roi_boundary = cell2mat(bwboundaries(roi_cat(:,:,curr_roi)));
    patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),roi_col(curr_roi,:));
    
    text(nanmean(curr_roi_boundary(:,2)),nanmean(curr_roi_boundary(:,1)), ...
        wf_roi(curr_roi).area,'FontSize',12,'HorizontalAlignment','center')
end
axis image off;















