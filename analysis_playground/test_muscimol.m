%% Test analysis for muscimol during recording
% (AP040 and AP041 at the moment)

%% Plot psychometrics and reaction times (only one choiceworld)

animal = 'AP041';
protocol = 'vanillaChoiceworld';
experiments = AP_find_experiments(animal,protocol);

% experiments = experiments([experiments.imaging] & [experiments.ephys]);
experiments = experiments([experiments.imaging] & ~[experiments.ephys]);

init_array = cell(size(experiments));
bhv = struct('move_t',init_array,'stim_contrastsides',init_array,'trial_choice',init_array);

for curr_day = 1:length(experiments)
    
    load_parts.cam = false;
    load_parts.imaging = false;
    load_parts.ephys = false;
    
    day = experiments(curr_day).day;
    experiment = experiments(curr_day).experiment(end);
    
    AP_load_experiment
    framerate = 35;
    
    % Get event-aligned activity
    raster_window = [-0.5,2];
    upsample_factor = 1;
    raster_sample_rate = 1/(framerate*upsample_factor);
    t = raster_window(1):raster_sample_rate:raster_window(2);
    
    % Get align times
    use_align = stimOn_times;
    use_align(isnan(use_align)) = 0;
    
    t_peri_event = bsxfun(@plus,use_align,t);
    t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
    
    %%% Trial-align wheel velocity
    event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
        wheel_velocity,t_peri_event);
    
    % Get reaction time for building regressors
    [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.03,[],2);
    move_idx(~move_trial) = NaN;
    move_t = nan(size(move_idx));
    move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
    
    % Get stim
    unique_stim = unique(contrasts(contrasts > 0).*sides');
    stim_contrastsides = ...
        signals_events.trialSideValues(1:length(stimOn_times))'.* ...
        signals_events.trialContrastValues(1:length(stimOn_times))';
    
    % Package
    bhv(curr_day).day = day;
    bhv(curr_day).move_t = move_t;
    bhv(curr_day).stim_contrastsides = stim_contrastsides;
    bhv(curr_day).trial_choice = trial_choice;
    
    AP_print_progress_fraction(curr_day,length(experiments));
    
end

% Plot all psychometrics and reaction times overlaid
figure; hold on
psychometric_ax = subplot(1,2,1); hold on;
line(psychometric_ax,[-1,1],[0.5,0.5],'linestyle','--','color','k');
line(psychometric_ax,[0,0],[0,1],'linestyle','--','color','k');

rxn_ax = subplot(1,2,2); hold on;
line(rxn_ax,[-1,1],[0.5,0.5],'linestyle','--','color','k');
line(rxn_ax,[0,0],[0,1],'linestyle','--','color','k');

% Plot psychometric and number of trials over days
unique_stim_contrastsides = unique(vertcat(bhv.stim_contrastsides));
frac_left_all = nan(length(bhv),11);
rxn_all = nan(length(bhv),11);
for i = 1:length(bhv)
    n_trials = length(bhv(i).move_t);
    frac_left = grpstats(bhv(i).trial_choice(1:n_trials) == -1,bhv(i).stim_contrastsides(1:n_trials));
    plot(psychometric_ax,unique(bhv(i).stim_contrastsides),frac_left,'k','linewidth',2);
    
    rxn = grpstats(bhv(i).move_t(1:n_trials),bhv(i).stim_contrastsides(1:n_trials),{'nanmedian'});
    plot(rxn_ax,unique(bhv(i).stim_contrastsides),rxn,'k','linewidth',2);
    
    curr_stim_contrastsides = unique(bhv(i).stim_contrastsides(1:n_trials));
    [~,curr_stim_contrastsides_idx] = ismember(curr_stim_contrastsides,unique_stim_contrastsides);
    frac_left_all(i,curr_stim_contrastsides_idx) = frac_left;
    rxn_all(i,curr_stim_contrastsides_idx) = rxn;
end

n_trials_all = cellfun(@length,{bhv.move_t});

figure; 

subplot(2,1,1);
plot(n_trials_all,'k','linewidth',2);
xlabel('Session');
ylabel('Number of trials');

subplot(2,1,2);
imagesc(frac_left_all');
colormap(brewermap([],'*RdBu'));

%% Plot psychometrics and reaction times (pre/post muscimol choiceworld)

animal = 'AP054';
protocol = 'vanillaChoiceworld';
flexible_name = true;
experiments = AP_find_experiments(animal,protocol,flexible_name);

experiments = experiments([experiments.imaging] & [experiments.ephys]);

init_array = cell(size(experiments));

curr_day = 4;
day = experiments(curr_day).day;

if length(experiments(curr_day).experiment) ~= 2
    error('~= 2 experiments')
end

bhv = struct('frac_go_left',init_array,'rxn_time',init_array);
for curr_exp = 1:2
    
    experiment = experiments(curr_day).experiment(curr_exp);
    
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
    
    % Get whether all contrasts were used
    use_all_contrasts = all(block.events.useContrastsValues(end-10:end));
    
    
    % Package final 
    bhv(curr_exp).frac_go_left = performance(3,:)./performance(2,:);
    bhv(curr_exp).rxn_time = stim_rxn_time;
    
end

figure; 
subplot(1,2,1); hold on;
plot(performance(1,:),bhv(1).frac_go_left,'k','linewidth',2);
plot(performance(1,:),bhv(2).frac_go_left,'r','linewidth',2);
line(xlim,[0.5,0.5],'color','k','linestyle','--');
line([0,0],[0,1],'color','k','linestyle','--');

subplot(1,2,2); hold on;
plot(performance(1,:),bhv(1).rxn_time,'k','linewidth',2);
plot(performance(1,:),bhv(2).rxn_time,'r','linewidth',2);
line([-1,1],[0.5,0.5],'color','k','linestyle','--');




