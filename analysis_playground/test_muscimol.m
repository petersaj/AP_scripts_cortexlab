%% Test analysis for muscimol during recording
% (AP040 and AP041 at the moment)

%% Plot psychometrics and reaction times

animal = 'AP041';
protocol = 'vanillaChoiceworld';
experiments = AP_find_experiments(animal,protocol);

experiments = experiments([experiments.imaging] & [experiments.ephys]);

init_array = cell(size(experiments));
bhv = struct('move_t',init_array,'stim_contrastsides',init_array,'trial_choice',init_array);

for curr_day = 1:length(experiments)
    
    load_parts.cam = false;
    load_parts.imaging = false;
    load_parts.ephys = false;
    
    day = experiments(curr_day).day;
    experiment = experiments(curr_day).experiment;
    
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
    bhv(curr_day).move_t = move_t;
    bhv(curr_day).stim_contrastsides = stim_contrastsides;
    bhv(curr_day).trial_choice = trial_choice;
    
    AP_print_progress_fraction(curr_day,length(experiments));
    
end


figure; hold on
psychometric_ax = subplot(1,2,1); hold on;
line(psychometric_ax,[-1,1],[0.5,0.5],'linestyle','--','color','k');
line(psychometric_ax,[0,0],[0,1],'linestyle','--','color','k');

rxn_ax = subplot(1,2,2); hold on;
line(rxn_ax,[-1,1],[0.5,0.5],'linestyle','--','color','k');
line(rxn_ax,[0,0],[0,1],'linestyle','--','color','k');

for i = 1:length(bhv)
    n_trials = length(bhv(i).move_t);
    frac_left = grpstats(bhv(i).trial_choice(1:n_trials) == -1,bhv(i).stim_contrastsides(1:n_trials));
    plot(psychometric_ax,unique(bhv(i).stim_contrastsides),frac_left,'k','linewidth',2);
    
    rxn = grpstats(bhv(i).move_t(1:n_trials),bhv(i).stim_contrastsides(1:n_trials),{'nanmedian'});
    plot(rxn_ax,unique(bhv(i).stim_contrastsides),rxn,'k','linewidth',2);
end














