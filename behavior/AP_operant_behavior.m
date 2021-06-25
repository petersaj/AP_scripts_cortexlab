%% --------- AP_stimWheelRight (aka one-sided choiceworld)

%% Get and plot mice behavior

% early corticostriatal mice
% animals =  {'AP089','AP090','AP091'};

% corticostriatal mice
animals = {'AP092','AP093','AP094','AP095','AP096','AP097'};

% tetO mice
% animals = {'AP100','AP101','AP103','AP104'};

% (current mice)
% animals = {'AP105','AP106'};

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

%             % (to get probability of move starts around stim)
%             stim_surround_t_edges = -5:0.1:5;
%             stim_surround_t_centers = stim_surround_t_edges(1:end-1) + diff(stim_surround_t_edges)/2;
%             stim_surround_times = block.events.stimOnTimes' + stim_surround_t_edges;
%             stim_surround_move = cell2mat(arrayfun(@(x) ...
%                 histcounts(wheel_starts,stim_surround_times(x,:)), ...
%                 1:size(stim_surround_times,1),'uni',false)');
%             % (need to do some kind of hazard rate here, maybe like what
%             % Ruediger did, I think there was another one like this too)
            
            % (probability of any movement around stim)
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

%             % (or another way: align by breaks within ITI, eg if mouse took
%             % a break from seconds 5-7, what's the profile of moving again?
%             % it's probably shallow, while the stim-related one is strong)
%             test_bin = [5,7];
%             trial_iti_time = wheel_t_resample - ...
%                 interp1(block.events.stimOnTimes,block.events.stimOnTimes, ...
%                 wheel_t_resample,'previous');
%             trial_iti_time_use = trial_iti_time >= test_bin(1) & ...
%                 trial_iti_time < test_bin(2);
            %(not done - I don't think necessary)



            % Time to reward
            stim_reward_t = reward_times - ...
                block.events.stimOnTimes(rewarded_trials);

            % Resetting quiescence period for each trial
            if isfield(block.events,'trialQuiescenceValues')
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
    plot(day_num,bhv(curr_animal).n_trials./bhv(curr_animal).session_duration,'linewidth',2);
    ylabel('Trials/min');
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
    t = bhv(curr_animal).stim_surround_t;
    move_prestim_max = ...
        cell2mat(cellfun(@(x) max(nanmean(x(:,t < 0),1)), ...
        bhv(curr_animal).stim_surround_wheel,'uni',false)');
    move_poststim_max = ...
        cell2mat(cellfun(@(x) max(nanmean(x(:,t > 0),1)), ...
        bhv(curr_animal).stim_surround_wheel,'uni',false)');
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
    clearvars('-except',preload_vars{:});
    
end

% End if only one mouse
if length(animals) == 1
    return
end

% If more than one mouse, plot group
% (only for min days across mice)
min_days = min(cellfun(@length,{bhv.n_trials}));
max_days = max(cellfun(@length,{bhv.n_trials}));

% Time to movement
stim_move_t_allpad = cell2mat(cellfun(@(x) padarray(cellfun(@nanmedian,x), ...
    [0,max_days-length(x)],NaN,'post'),{bhv.stim_move_t}','uni',false));
stim_move_t_std_allpad = cell2mat(cellfun(@(x) padarray(cellfun(@nanstd,x), ...
    [0,max_days-length(x)],NaN,'post'),{bhv.stim_move_t}','uni',false));

% Time to reward
stim_reward_t_allpad = cell2mat(cellfun(@(x) padarray(cellfun(@nanmedian,x), ...
    [0,max_days-length(x)],NaN,'post'),{bhv.stim_reward_t}','uni',false));
stim_reward_t_std_allpad = cell2mat(cellfun(@(x) padarray(cellfun(@nanstd,x), ...
    [0,max_days-length(x)],NaN,'post'),{bhv.stim_reward_t}','uni',false));

% Movement around stim
stim_surround_wheel_avg = cell2mat(cellfun(@(x) ...
    padarray(cell2mat(cellfun(@(x) nanmean(x,1),x,'uni',false)'), ...
    [max_days-length(x),0],NaN,'post'), ...
    permute({bhv.stim_surround_wheel},[1,3,2]),'uni',false));

stim_surround_t = bhv(1).stim_surround_t;
figure; hold on
set(gca,'ColorOrder',copper(size(stim_surround_wheel_avg,1)));
plot(stim_surround_t,nanmean(stim_surround_wheel_avg,3)','linewidth',2);
xlabel('Time from stim (s)');
ylabel('Fraction moving');


%% Muscimol behavior (after loading above)

% Load muscimol injection info
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
muscimol_fn = [data_path filesep 'muscimol.mat'];
load(muscimol_fn);

muscimol_animals = find(ismember({muscimol.animal},animals));
muscimol_areas = unique(vertcat(muscimol(muscimol_animals).area));

stim_surround_t = bhv(1).stim_surround_t;
muscimol_stim_surround_wheel = repmat( ...
    {nan(length(muscimol_animals),length(stim_surround_t))}, ...
    length(muscimol_areas),1);

figure;
for curr_area_idx = 1:length(muscimol_areas)
    curr_area = muscimol_areas{curr_area_idx};
    for curr_animal_idx = 1:length(muscimol_animals)
        muscimol_day_idx = find(strcmp(curr_area, ...
            muscimol(muscimol_animals(curr_animal_idx)).area));
        % Use last available day (doubled once when didn't work)
        use_day = muscimol(muscimol_animals(curr_animal_idx)).day{muscimol_day_idx(end)};
        bhv_day_idx = ismember(bhv(curr_animal_idx).day,use_day);
        
        curr_stim_surround_wheel = nanmean( ...
            bhv(curr_animal_idx).stim_surround_wheel{bhv_day_idx},1);
        
        muscimol_stim_surround_wheel{curr_area_idx}(curr_animal_idx,:) = ...
            curr_stim_surround_wheel;
    end   
    
    % Plot each animal and each condition
    subplot(length(muscimol_areas)+1,1,curr_area_idx); hold on;
    plot(stim_surround_t,muscimol_stim_surround_wheel{curr_area_idx}');
    plot(stim_surround_t,nanmean(muscimol_stim_surround_wheel{curr_area_idx},1),'k','linewidth',2);
    ylim([0,1]);
    title(curr_area);
end

% Plot means overlaid
subplot(length(muscimol_areas)+1,1,length(muscimol_areas)+1); hold on;
plot(stim_surround_t, ...
    cell2mat(cellfun(@(x) nanmean(x,1),muscimol_stim_surround_wheel,'uni',false))', ...
    'linewidth',2);
legend(muscimol_areas);







    




    





    



    



    



    

