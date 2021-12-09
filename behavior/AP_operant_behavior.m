%% AP_stimWheelRight (aka one-sided choiceworld)

%%%%% CHANGE HOW THIS IS DONE: LOAD MUSCIMOL INFO DIRECTLY INTO GRAB?


%% ~~~~~~~ Grab behavior

%% Get and plot mice behavior

% early corticostriatal mice (fixed task timings)
% animal_group = 'cstr_fixtime';
% animals =  {'AP089','AP090','AP091'};

% % corticostriatal mice
% animal_group = 'cstr';
% animals = {'AP092','AP093','AP094','AP095','AP096','AP097'};

% % tetO mice
% animal_group = 'teto';
% animals = {'AP100','AP101','AP103','AP104','AP105','AP106'};

% Andrada's tetO mice
animal_group = 'teto_AM';
animals = {'AP107','AP108','AP109'};

protocol = 'AP_stimWheelRight';
flexible_name = false;
bhv = struct;

% Save flag
save_bhv = false;

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
            n_trials = length(block.events.responseValues);
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
            
            %%%%%%%%%% TESTING
            
            %             %%% (Real vs null reaction times
            %             % Get real and null distribution of stim-to-move times
            %             % (get gaps between iti movements and prior movements)
            %             use_wheel_move_iti = wheel_move_iti_idx(wheel_move_iti_idx > 1); % no gap for first movement
            %             iti_move_gaps = wheel_starts(use_wheel_move_iti) - wheel_stops(use_wheel_move_iti - 1);
            %             % (get gaps between rewarded movements and prior movements)
            %             use_wheel_move_response = wheel_move_response_idx(wheel_move_response_idx > 1); % no gap for first movement
            %             wheel_move_preresponse_stop = wheel_stops(use_wheel_move_response-1);
            %             % (make null distribution: positive prior movement + gap relative to stim)
            %             % (from 2nd stim, in case no move before first)
            %             n_shuff = 10000;
            %             stim_to_move_shuff = ...
            %                 stimOn_times(end-length(use_wheel_move_response)+1:end) ...
            %                 - (wheel_move_preresponse_stop + ...
            %                 reshape(randsample(iti_move_gaps,length(wheel_move_preresponse_stop)*n_shuff,true),[],n_shuff));
            %             stim_to_move_shuff(stim_to_move_shuff < 0) = NaN; % positives only, negatives would've been resets
            %             % (p-value of reaction time)
            %             stim_to_move_shuff_med = nanmedian(stim_to_move_shuff,1);
            %             stim_to_move_med = nanmedian(stim_to_move);
            %             stim_to_move_rank = tiedrank([stim_to_move_med,stim_to_move_shuff_med]);
            %             stim_to_move_p = stim_to_move_rank(1)./max(stim_to_move_rank);
            %             % (reaction time relative to random)
            %             stim_to_move_shuff_med_total = nanmedian(stim_to_move_shuff(:));
            %             stim_to_move_shuffcomp = (stim_to_move_med - stim_to_move_shuff_med_total)./ ...
            %                 (stim_to_move_med + stim_to_move_shuff_med_total);
            
            
%             % Cox proportional regression
%             % Regress time from move to post-stim move relative to pre-stim gap/stim
%             use_last_gaps = 5;
%             wheel_gaps = wheel_starts(2:end)-wheel_stops(1:end-1);
%             wheel_gaps_regressor_idx = (wheel_move_stim_idx-1)+[-use_last_gaps:0];
%             use_move_regressors = all(wheel_gaps_regressor_idx>0,2);
%             
%             wheel_gap_predict = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors,end));
%             wheel_gaps_regressors = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors,1:end-1));
%             wheel_stim_regressor = stimOn_times(use_move_regressors) - wheel_stops(wheel_move_stim_idx(use_move_regressors)-1);
%             
%             [b,logl_stim,H,stats] = coxphfit([wheel_gaps_regressors,wheel_stim_regressor],wheel_gap_predict);
%             [b_nostim,logl_nostim,H_nostim,stats_nostim] = coxphfit([wheel_gaps_regressors],wheel_gap_predict);
%             
%             logdiff = logl_stim-logl_nostim;           
%                         
%             bhv(curr_animal).b{curr_day} = b;
%             bhv(curr_animal).logdiff(curr_day) = logdiff;
            
            
            % Store in behavior structure
            bhv(curr_animal).animal = animal;
            bhv(curr_animal).day{curr_day} = day;
            bhv(curr_animal).protocol{curr_day} = curr_protocol;
            bhv(curr_animal).session_duration(curr_day) = session_duration;
            bhv(curr_animal).n_trials(curr_day) = n_trials;
            bhv(curr_animal).frac_hit(curr_day) = frac_hit;
            bhv(curr_animal).total_water(curr_day) = total_water;
            bhv(curr_animal).wheel_velocity(curr_day) = nansum(abs(wheel_velocity));
            bhv(curr_animal).stim_move_t{curr_day} = stim_to_move;
            bhv(curr_animal).stim_reward_t{curr_day} = stim_to_feedback(trial_outcome == 1);
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
    stim_surround_t = bhv(curr_animal).stim_surround_t;
    move_prestim_max = ...
        cell2mat(cellfun(@(x) max(nanmean(x(:,stim_surround_t < 0),1)), ...
        bhv(curr_animal).stim_surround_wheel,'uni',false)');
    move_poststim_max = ...
        cell2mat(cellfun(@(x) max(nanmean(x(:,stim_surround_t > 0),1)), ...
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

if save_bhv
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    save_fn = [data_path filesep 'bhv_' animal_group];
    save(save_fn,'bhv');
    disp(['Saved ' save_fn]);
end


%% (load data saved after above)

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';

% (load separate)
bhv_fn = [data_path filesep 'bhv_teto'];
% bhv_fn = [data_path filesep 'bhv_corticostriatal'];
load(bhv_fn);

% % (load and combine)
% bhv_fn1 = [data_path filesep 'bhv_teto'];
% bhv1 = load(bhv_fn1);
% bhv_fn2 = [data_path filesep 'bhv_corticostriatal'];
% bhv2 = load(bhv_fn2);
% bhv = [bhv1.bhv,bhv2.bhv];
% clearvars bhv1 bhv2

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
    if isempty(muscimol_animal_idx)
        continue
    end
    
    pre_muscimol_day_idx = datenum(bhv(curr_animal).day) < ...
        datenum(muscimol(muscimol_animal_idx).day(1));
    muscimol_day_idx = ismember(datenum(bhv(curr_animal).day), ...
        datenum(muscimol(muscimol_animal_idx).day));
    
    use_days{curr_animal} = pre_muscimol_day_idx;
end

stim_move_t_cat = cellfun(@(x,y) cell2mat(x(y)),{bhv(:).stim_move_t},use_days,'uni',false);
stim_reward_t_cat = cellfun(@(x,y) cell2mat(x(y)),{bhv(:).stim_reward_t},use_days,'uni',false);
move_reward_t_cat = cellfun(@(x,y) y-x,stim_move_t_cat,stim_reward_t_cat,'uni',false);

trials_cumsum = cellfun(@(x,y) cumsum(cellfun(@length,x(y))),{bhv(:).stim_move_t},use_days,'uni',false);

figure;
for curr_animal = 1:length(animals)
    subplot(length(animals),1,curr_animal); hold on;
    plot(stim_move_t_cat{curr_animal},'.k');
    plot(stim_reward_t_cat{curr_animal},'.b');
    plot(move_reward_t_cat{curr_animal},'.m');
    set(gca,'YScale','log')
    ylim([0.05,10])
    xlim([0,trials_cumsum{curr_animal}(end)]);
    for i = 1:length(trials_cumsum{curr_animal})
        line(repmat(trials_cumsum{curr_animal}(i),2,1), ...
            ylim,'color','r','linewidth',2)
    end
    ylabel(animals{curr_animal});
    
    muscimol_animal_idx = ismember({muscimol.animal},bhv(curr_animal).animal);
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

% TESTING: fit slope each day
p = cellfun(@(x,y) cellfun(@(x,y) polyfit(1:length(x),log10(x),1),x(y),'uni',false), ...
    {bhv(:).stim_move_t},use_days,'uni',false);
p_offset = cellfun(@(x) cellfun(@(x) x(2),x),p,'uni',false);
p_slope = cellfun(@(x) cellfun(@(x) x(1),x),p,'uni',false);


%% TESTING: histograms of reaction times

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
    muscimol_day_idx = datenum(bhv(curr_animal).day) >= ...
        datenum(muscimol(muscimol_animal_idx).day(1));
    use_days{curr_animal} = ~muscimol_day_idx;
end

max_days = max(cellfun(@sum,use_days));

% Animal/day
figure;
for curr_animal = 1:length(bhv)
    for curr_day = find(use_days{curr_animal})'
        subplot(max_days,length(bhv),(curr_day-1)*length(bhv)+curr_animal); hold on;
        histogram(bhv(curr_animal).stim_move_t{curr_day},0:0.02:1,'EdgeColor','none','normalization','probability')
    end
end

% Day (animals overlaid)
figure;
for curr_animal = 1:length(bhv)
    for curr_day = 1:10
        subplot(10,1,curr_day); hold on;
        histogram(bhv(curr_animal).stim_move_t{curr_day},0:0.02:1,'EdgeColor','none','normalization','probability')
    end
end

% Animals (days combined)
figure; hold on;

cellfun(@(x) histogram(cell2mat(x),0:0.005:1,'EdgeColor','none', ...
    'normalization','probability'),{bhv.stim_move_t},'uni',false);
cellfun(@(x) line(repmat(mode(round(cell2mat(x)*100)/100),2,1),ylim,'linewidth',2),{bhv.stim_move_t});

% Total
figure;
histogram(cell2mat([bhv.stim_move_t]),0:0.005:1,'EdgeColor','none','normalization','probability');

% Fraction rxn within boundary
b = cellfun(@(x) cellfun(@(x) nanmean(x > 0.1 & x < 0.25),x),{bhv.stim_move_t},'uni',false);
figure; hold on;
cellfun(@plot,b);



% Day-split reaction time fraction

n_daysplit = 4;

% Index: within-day split, day, animal
trial_split_idx = cellfun(@(x,use_days,animal_num) ...
    cellfun(@(x,day) ...
    [min(floor(linspace(1,n_daysplit+1,length(x))),n_daysplit); ...
    repmat(day,1,length(x));repmat(animal_num,1,length(x))], ...
    x(use_days),num2cell(1:sum(use_days)),'uni',false), ...
    {bhv.stim_move_t},use_days,num2cell(1:length(bhv)),'uni',false);

stim_move_t_cat = cell2mat(cellfun(@(x,use_days) cell2mat(x(use_days)), ...
    {bhv.stim_move_t},use_days,'uni',false));

stim_move_t_stimtime = accumarray(cell2mat([trial_split_idx{:}])', ...
    stim_move_t_cat > 0.1 & stim_move_t_cat < 0.25',[n_daysplit,max_days,length(bhv)],@nanmean,NaN);

figure; hold on;
stim_move_t_stimtime_long = reshape(padarray(stim_move_t_stimtime,[1,0,0],NaN,'post'),[],length(animals));
plot([1:size(stim_move_t_stimtime_long,1)]/(n_daysplit+1),stim_move_t_stimtime_long);
errorbar([1:size(stim_move_t_stimtime_long,1)]/(n_daysplit+1), ...
    nanmean(stim_move_t_stimtime_long,2),AP_sem(stim_move_t_stimtime_long,2),'k','linewidth',2);
xlabel('Day');
ylabel('Rxn in stim window');




%% Non-muscimol behavior

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
    muscimol_day_idx = datenum(bhv(curr_animal).day) >= ...
        datenum(muscimol(muscimol_animal_idx).day(1));
    use_days{curr_animal} = ~muscimol_day_idx;
end

% (get max days for padding)
max_days = max(cellfun(@sum,use_days));

% Time to movement
stim_move_t_allpad = cell2mat(cellfun(@(x) padarray(cellfun(@nanmedian,x), ...
    [0,max_days-length(x)],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_move_t},use_days,'uni',false)','uni',false));
stim_move_t_std_allpad = cell2mat(cellfun(@(x) padarray(cellfun(@nanstd,x), ...
    [0,max_days-length(x)],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_move_t},use_days,'uni',false)','uni',false));

% Time to reward
stim_reward_t_allpad = cell2mat(cellfun(@(x) padarray(cellfun(@nanmedian,x), ...
    [0,max_days-length(x)],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_reward_t},use_days,'uni',false)','uni',false));
stim_reward_t_std_allpad = cell2mat(cellfun(@(x) padarray(cellfun(@nanstd,x), ...
    [0,max_days-length(x)],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_reward_t},use_days,'uni',false)','uni',false));

% Movement around stim
stim_surround_wheel_avg = cell2mat(cellfun(@(x) ...
    padarray(cell2mat(cellfun(@(x) nanmean(x,1),x,'uni',false)'), ...
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
trials_per_min = cell2mat(cellfun(@(x) padarray(x,[0,max_days-length(x)],NaN,'post'), ...
    cellfun(@(x,y,use_days) x(use_days)./y(use_days), ....
    {bhv.n_trials},{bhv.session_duration},use_days,'uni',false),'uni',false)');

% Total ITI (initial + resets)
iti_allpad = cell2mat(cellfun(@(x) padarray(cellfun(@nanmedian,x), ...
    [0,max_days-length(x)],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.iti_t},use_days,'uni',false)','uni',false));

% Plot
figure;

subplot(1,3,1);
hold on;
errorbar(nanmean(stim_move_t_allpad,1),AP_sem(stim_move_t_allpad,1),'k','linewidth',2);
errorbar(nanmean(stim_reward_t_allpad,1),AP_sem(stim_reward_t_allpad,1),'b','linewidth',2);
ylabel('Time (s)');
xlabel('Session');
legend({'Stim to move','Stim to reward'});
set(gca,'YScale','log');
set(gca,'YTick',[0.25,0.5,1,2,4]);

subplot(1,3,2);
hold on;
stim_surround_t = bhv(1).stim_surround_t;
set(gca,'ColorOrder',copper(size(stim_surround_wheel_avg,1)));
plot(stim_surround_t,nanmean(stim_surround_wheel_avg,3)','linewidth',2);
xlabel('Time from stim (s)');
ylabel('Fraction moving');

subplot(1,3,3);
hold on;
plot(move_prepost_max_ratio');
errorbar(nanmean(move_prepost_max_ratio,1),AP_sem(move_prepost_max_ratio,1),'k','linewidth',2);
line(xlim,[0,0],'linestyle','--','color','k');
ylabel('Pre/post stim move ratio');
xlabel('Session');


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
    {nan(length(muscimol_animals),length(stim_surround_t))}, ...
    length(muscimol_areas),1);

stim_move_t_muscimol = nan(length(muscimol_areas),length(animals));
stim_reward_t_muscimol = nan(length(muscimol_areas),length(animals));
move_reward_t_muscimol = nan(length(muscimol_areas),length(animals));

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
            bhv(curr_bhv_animal_idx).stim_surround_wheel(bhv_day_idx),'uni',false)'),1);
        
        muscimol_stim_surround_wheel{curr_area_idx}(curr_muscimol_animal_idx,:) = ...
            curr_stim_surround_wheel;
        
        stim_move_t_muscimol(curr_area_idx,curr_muscimol_animal_idx) = ...
            nanmedian([bhv(curr_bhv_animal_idx).stim_move_t{bhv_day_idx}]);
        stim_reward_t_muscimol(curr_area_idx,curr_muscimol_animal_idx) = ...
            nanmedian([bhv(curr_bhv_animal_idx).stim_reward_t{bhv_day_idx}]);
        move_reward_t_muscimol(curr_area_idx,curr_muscimol_animal_idx) = ...
            nanmedian([bhv(curr_bhv_animal_idx).stim_reward_t{bhv_day_idx}] - ...
            [bhv(curr_bhv_animal_idx).stim_move_t{bhv_day_idx}]);
        
        
        wheel_velocity_muscimol(curr_area_idx,curr_muscimol_animal_idx) = ...
            nanmean(bhv(curr_bhv_animal_idx).wheel_velocity(bhv_day_idx)./ ...
            bhv(curr_bhv_animal_idx).session_duration(bhv_day_idx));
        wheel_bias_muscimol(curr_area_idx,curr_muscimol_animal_idx) = ...
            nanmean(bhv(curr_bhv_animal_idx).wheel_bias(bhv_day_idx));
        
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

% Movement probability pre/post stim
stim_surround_t = bhv(1).stim_surround_t;
move_prestim_max = cellfun(@(x) max(x(:,stim_surround_t<0,:),[],2),muscimol_stim_surround_wheel,'uni',false);
move_poststim_max = cellfun(@(x) max(x(:,stim_surround_t>0,:),[],2),muscimol_stim_surround_wheel,'uni',false);
move_prepost_max_ratio = cell2mat(cellfun(@(x,y) (y-x)./(x+y),move_prestim_max,move_poststim_max,'uni',false)')';

figure;
subplot(1,5,1); hold on;
plot(stim_move_t_muscimol,'linewidth',2)
set(gca,'YScale','log')
set(gca,'XTick',1:length(muscimol_areas),'XTickLabel',muscimol_areas);
ylabel('Stim move t');

subplot(1,5,2); hold on;
plot(stim_reward_t_muscimol,'linewidth',2)
set(gca,'YScale','log')
set(gca,'XTick',1:length(muscimol_areas),'XTickLabel',muscimol_areas);
ylabel('Stim reward t');

subplot(1,5,3); hold on;
plot(move_reward_t_muscimol,'linewidth',2)
set(gca,'YScale','log')
set(gca,'XTick',1:length(muscimol_areas),'XTickLabel',muscimol_areas);
ylabel('Move reward t');

subplot(1,5,4); hold on;
plot(wheel_velocity_muscimol,'linewidth',2)
set(gca,'YScale','log')
set(gca,'XTick',1:length(muscimol_areas),'XTickLabel',muscimol_areas);
ylabel('Wheel velocity/s');

subplot(1,5,5); hold on;
plot(move_prepost_max_ratio,'linewidth',2)
set(gca,'XTick',1:length(muscimol_areas),'XTickLabel',muscimol_areas);
ylabel('Post/post ratio');



figure;
subplot(1,2,1); hold on;
gscatter(move_prepost_max_ratio(:),wheel_velocity_muscimol(:), ...
    reshape(repmat(muscimol_areas,1,length(animals)),[],1),lines(length(muscimol_areas)),[],20);
xlabel('Pre/post move');
ylabel('Velocity/s');
axis square;

subplot(1,2,2); hold on;
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
    if isempty(muscimol_animal_idx)
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
    [min(floor(linspace(1,n_daysplit+1,length(x))),n_daysplit); ...
    repmat(day,1,length(x));repmat(animal_num,1,length(x))], ...
    x(use_days),num2cell(1:sum(use_days)),'uni',false), ...
    {bhv.stim_move_t},use_days,num2cell(1:length(bhv)),'uni',false);

stim_move_t_cat = cell2mat(cellfun(@(x,use_days) cell2mat(x(use_days)), ...
    {bhv.stim_move_t},use_days,'uni',false));
stim_move_t_daysplit = accumarray(cell2mat([trial_split_idx{:}])', ...
    stim_move_t_cat',[n_daysplit,max_days,length(bhv)],@nanmedian,NaN);

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
    padarray(cell2mat(cellfun(@(x) nanmean(x,1),x,'uni',false)'), ...
    [max_days-length(x),0],NaN,'post'), ...
    permute( ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_surround_wheel},use_days,'uni',false), ...
    [1,3,2]),'uni',false));

stim_surround_wheel_cat = cell2mat(cellfun(@(x,use_days) cell2mat(x(use_days)'), ...
    {bhv.stim_surround_wheel},use_days,'uni',false)');

stim_surround_wheel_grp = accumarray(reshape(cat(3, ...
    repmat(1:size(stim_surround_wheel_cat,2),size(stim_surround_wheel_cat,1),1), ...
    repmat(permute(cell2mat([trial_split_idx{:}])',[1,3,2]),1,size(stim_surround_wheel_cat,2))),[],4), ...
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








