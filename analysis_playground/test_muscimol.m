%% ~~~~~~ Test analysis for muscimol during recording ~~~~~~
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

animal = 'AP053';
protocol = 'vanillaChoiceworld';
flexible_name = true;
experiments = AP_find_experiments(animal,protocol,flexible_name);

experiments = experiments([experiments.imaging] & [experiments.ephys]);

init_array = cell(size(experiments));

curr_day = 1;%length(experiments);
day = experiments(curr_day).day;

if length(experiments(curr_day).experiment) ~= 2
    error('~= 2 experiments')
end

bhv = struct('frac_orient_right',init_array,'rxn_time',init_array);
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
    bhv(curr_exp).frac_orient_right = performance(3,:)./performance(2,:);
    bhv(curr_exp).rxn_time = stim_rxn_time;
    
end

figure;
subplot(1,2,1); hold on;
plot(performance(1,:),bhv(1).frac_orient_right,'k','linewidth',2);
plot(performance(1,:),bhv(2).frac_orient_right,'r','linewidth',2);
line(xlim,[0.5,0.5],'color','k','linestyle','--');
line([0,0],[0,1],'color','k','linestyle','--');

subplot(1,2,2); hold on;
plot(performance(1,:),bhv(1).rxn_time,'k','linewidth',2);
plot(performance(1,:),bhv(2).rxn_time,'r','linewidth',2);
line([-1,1],[0.5,0.5],'color','k','linestyle','--');

%% Plot psychometrics and reaction times (pre/post muscimol choiceworld, BATCH)

animals = {'AP045','AP054','AP055','AP053','AP047','AP048'};

bhv = struct('frac_orient_right',cell(length(animals),1), ...
    'rxn_time',cell(length(animals),1));

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        
        if length(experiments(curr_day).experiment) ~= 2
            error('~= 2 experiments')
        end
        
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
            stim_rxn_time_gocue = accumarray(trial_stim_idx(response_trials)',trial_move_t' > 0.5,[11,1],@nanmean,nan);
            
            % Performance (note that this excludes repeat on incorrect trials)
            performance = block.events.sessionPerformanceValues(:,end-10:end);
            
            % Get whether all contrasts were used
            use_all_contrasts = all(block.events.useContrastsValues(end-10:end));
            
            
            % Package final
            bhv(curr_animal).frac_orient_right{curr_day}(curr_exp,:) = performance(3,:)./performance(2,:);
            bhv(curr_animal).rxn_time{curr_day}(curr_exp,:) = stim_rxn_time;
            bhv(curr_animal).rxn_time_gocue{curr_day}(curr_exp,:) = stim_rxn_time_gocue;
            
        end
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
end

% Save for quick loading later
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\tests\muscimol_test';
save_fn = [save_path filesep 'muscimol_bhv.mat'];
save(save_fn,'bhv');
disp(['Saved ' save_fn]);

% Plot
frac_orient_right_cat = cell2mat(permute(cat(2,bhv.frac_orient_right),[1,3,2]));
rxn_time_cat = cell2mat(permute(cat(2,bhv.rxn_time),[1,3,2]));
rxn_time_gocue_cat = cell2mat(permute(cat(2,bhv.rxn_time_gocue),[1,3,2]));

contrastsides = unique([1,0.5,0.25,0.125,0.06,0].*[-1;1]);

figure;
subplot(2,2,1); hold on;
plot(contrastsides,squeeze(frac_orient_right_cat(1,:,:)),'color',[0.5,0.5,0.5],'linewidth',1);
plot(contrastsides,squeeze(frac_orient_right_cat(2,:,:)),'color',[1,0.5,0.5],'linewidth',1);
plot(contrastsides,nanmean(frac_orient_right_cat(1,:,:),3),'color','k','linewidth',3);
plot(contrastsides,nanmean(frac_orient_right_cat(2,:,:),3),'color','r','linewidth',3);
line(xlim,[0.5,0.5],'color','k','linestyle','--');
line([0,0],[0,1],'color','k','linestyle','--');
xlabel('Side*Contrast');
ylabel('Fraction orient right')

subplot(2,2,2); hold on;
plot(contrastsides,squeeze(rxn_time_cat(1,:,:)),'color',[0.5,0.5,0.5],'linewidth',1);
plot(contrastsides,squeeze(rxn_time_cat(2,:,:)),'color',[1,0.5,0.5],'linewidth',1);
plot(contrastsides,nanmean(rxn_time_cat(1,:,:),3),'color','k','linewidth',3);
plot(contrastsides,nanmean(rxn_time_cat(2,:,:),3),'color','r','linewidth',3);
line([-1,1],[0.5,0.5],'color','k','linestyle','--');
xlabel('Side*Contrast');
ylabel('Reaction time');

subplot(2,2,3); hold on;
plot(contrastsides,squeeze(diff(frac_orient_right_cat,[],1)),'color',[0.5,0.5,0.5],'linewidth',1);
plot(contrastsides,nanmean(diff(frac_orient_right_cat,[],1),3),'color','k','linewidth',3);
line(xlim,[0,0],'color','k','linestyle','--');
line([0,0],[-1,1],'color','k','linestyle','--');
xlabel('Side*Contrast');
ylabel('\Delta Fraction orient right')

subplot(2,2,4); hold on;
plot(contrastsides,squeeze(diff(rxn_time_cat,[],1)),'color',[0.5,0.5,0.5],'linewidth',1);
plot(contrastsides,nanmean(diff(rxn_time_cat,[],1),3),'color','k','linewidth',3);
line([-1,1],[0,0],'color','k','linestyle','--');
xlabel('Side*Contrast');
ylabel('\Delta Reaction time');


%% Get VFS pre/post musicmol

animals = {'AP045','AP054','AP055','AP053','AP047','AP048'};

vfs = cell(length(animals),1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_sparseNoise';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        
        if length(experiments(curr_day).experiment) ~= 2
            warning([animal ' ' day '~= 2 experiments, using first/last']);
            experiments(curr_day).experiment = experiments(curr_day).experiment([1,end]);
        end
        
        for curr_exp = 1:2
            
            experiment = experiments(curr_day).experiment(curr_exp);
            
            AP_load_experiment;
            lilrig_retinotopy;
            vfs_aligned = AP_align_widefield(vfs_median,animal,day);
            
            % Package final
            vfs{curr_animal}{curr_day,curr_exp} = vfs_aligned;
            
        end
        AP_print_progress_fraction(curr_day,length(experiments));
    end
end

% Plot pre/post for each animal and day
for curr_animal = 1:length(vfs)
    figure('Name',animals{curr_animal});
    for curr_day = 1:size(vfs{curr_animal},1)
        subplot(size(vfs{curr_animal},1),2,(curr_day-1)*2+1);
        imagesc(vfs{curr_animal}{curr_day,1});
        axis image off;
        colormap(brewermap([],'*RdBu'));
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        
        subplot(size(vfs{curr_animal},1),2,(curr_day-1)*2+2);
        imagesc(vfs{curr_animal}{curr_day,2});
        axis image off;
        colormap(brewermap([],'*RdBu'));
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    end
end

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\tests\muscimol_test';
save_fn = 'muscimol_vfs.mat';
save([save_path filesep save_fn],'vfs');
disp(['Saved ' [save_path filesep save_fn]]);



%% ~~~~~~~~ Non-trial data


%% Cortical: spike rate and CSD
% NOTE: this doesn't make sense for AP058 in barrel cortex (visual striatum
% below, but I think the muscimol blew out visual cortex too because those
% responses are gone)

animal_days = { ...
    'AP052','2019-09-20';
    'AP058','2019-12-06'};

figure;

for curr_animalday = 1:length(animal_days)
    
    animal = animal_days{curr_animalday,1};
    day = animal_days{curr_animalday,2};
    curr_protocols = AP_list_experiments(animal,day);
    
    protocol = 'AP_lcrGratingPassive';
    curr_use_exp = [curr_protocols(strmatch(protocol, ...
        {curr_protocols.protocol})).experiment];
    
    % (2 expts: 1 pre-muscimol, 2 post-muscimol)
    for curr_exp = 1:length(curr_use_exp)
        
        experiment = curr_use_exp(curr_exp);
        lfp_channel = 1;
        AP_load_experiment;
        
        % Estimate boundaries of cortex (the dumb way: first template/gap)
        sorted_template_depths = sort(template_depths);
        ctx_start = sorted_template_depths(1) - 1;
        [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
        ctx_end = sorted_template_depths(max_gap_idx)+1;
        ctx_depth = [ctx_start,ctx_end];
        
        % Get stim-aligned average LFP
        % (craniotomies on the right, use left stimuli)
        use_stims = trial_conditions(:,2) == -90;
        use_stimOn_times = stimOn_times(use_stims);
        
        stim_surround_window = [-0.05,0.1];
        stim_lfp_t = stim_surround_window(1):1/lfp_sample_rate:stim_surround_window(2);
        
        pull_times = use_stimOn_times + stim_lfp_t;
        
        lfp_stim_align_avg = nan(n_channels,length(stim_lfp_t));
        
        for curr_lfp_channel = 1:n_channels
            
            lfp = double(lfp_memmap{lfp_load_file}.Data.lfp(curr_lfp_channel,lfp_load_start_rel:lfp_load_stop_rel));
            
            % Get average stim response
            lfp_stim_align = interp1(lfp_t_timeline,lfp,pull_times);
            
            baseline_time = stim_lfp_t < 0;
            lfp_stim_align_avg(curr_lfp_channel,:) = nanmean((lfp_stim_align - nanmean(lfp_stim_align(:,baseline_time),2)),1);
            
            AP_print_progress_fraction(curr_lfp_channel,n_channels);
            
        end
        
        % Average channels at same depth
        lfp_connected_channels = channel_map_full.connected;
        [stim_lfp_depth,stim_lfp] = grpstats( ...
            lfp_stim_align_avg(lfp_connected_channels,:), ...
            lfp_channel_positions(lfp_connected_channels),{'gname','mean'});
        stim_lfp_depth = cellfun(@str2num,stim_lfp_depth);
        
        % Calculate CSD (2nd derivative, smooth at each step)
        n_channel_smooth = 10;
        stim_lfp_smooth = movmean(stim_lfp,n_channel_smooth,1);
        stim_lfp_smooth_d1 = movmean(diff(stim_lfp_smooth,1,1),n_channel_smooth,1);
        stim_lfp_smooth_d2 = movmean(diff(stim_lfp_smooth_d1,1,1),n_channel_smooth,1);
        stim_csd = -stim_lfp_smooth_d2;
        stim_csd_depth = conv(conv(stim_lfp_depth,[1,1]/2,'valid'),[1,1]/2,'valid');
        
        % Get CSD profile at time slice
        csd_slice_t = [0.04,0.06]; % time to use after stim (s)
        csd_slice_samples = stim_lfp_t >= csd_slice_t(1) & stim_lfp_t <= csd_slice_t(2);
        csd_slice = nanmean(stim_csd(:,csd_slice_samples),2);
        
        % Get spike rate in current experiment
        curr_exp_start = sync(2).timestamps(find(experiment_idx)*2 - 1);
        curr_exp_stop = sync(2).timestamps(find(experiment_idx)*2);
        curr_use_spikes = spike_times >= curr_exp_start & ...
            spike_times <= curr_exp_stop;
        curr_spike_rate = accumarray(spike_templates(curr_use_spikes),1, ...
            [max(spike_templates),1])./(curr_exp_stop - curr_exp_start);
        
        % Plot CSD and slice
        switch curr_exp
            case 1
                cond_label = 'pre-muscimol';
            case 2
                cond_label = 'post-muscimol';
        end
        
        subplot(length(animal_days),4,(curr_animalday-1)*4+(curr_exp-1)*2+1,'YDir','reverse'); hold on;
        plot(log10(curr_spike_rate),template_depths,'.k','MarkerSize',10);
        line(xlim,repmat(ctx_start,2,1));
        line(xlim,repmat(ctx_end,2,1));
        xlabel('Log spike rate');
        ylabel('Depth (\mum)');
        title({animal,day,['CSD ' cond_label]});
        
        subplot(length(animal_days),4,(curr_animalday-1)*4+(curr_exp-1)*2+2);
        imagesc(stim_lfp_t,stim_csd_depth,stim_csd);
        colormap(brewermap([],'*RdBu'));
        caxis([-max(abs(caxis)),max(abs(caxis))]);
        line([0,0],ylim,'color','k');
        xlabel('Time from stim (s)');
        ylabel('Depth (\mum)');
        title({animal,day,['Spike rate ' cond_label]});
        
        linkaxes(get(gcf,'Children'),'y');
        drawnow;
        
    end
    
end


%% Cortical spike rate pre/post muscimol
% (use spike rate over all experiments pre/post muscimol)

animal_days = { ...
    'AP052','2019-09-20';
    'AP058','2019-12-06'};

figure;

for curr_animalday = 1:length(animal_days)
    
    animal = animal_days{curr_animalday,1};
    day = animal_days{curr_animalday,2};
    
    % Load data (first experiment - but spikes throughout used)
    experiment = 1;
    AP_load_experiment
    
    % Estimate boundaries of cortex (the dumb way: first template/gap)
    sorted_template_depths = sort(template_depths);
    ctx_start = sorted_template_depths(1) - 1;
    [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
    ctx_end = sorted_template_depths(max_gap_idx)+1;
    ctx_depth = [ctx_start,ctx_end];
    ctx_templates = template_depths <= ctx_depth(2);
    
    % Set experiments in conditions (1-2 = pre-muscimol, 3-4 = post-muscimol)
    cond_expts = {[1,2],[3,4]};
    
    spike_rate_cond = nan(max(spike_templates),2);
    for curr_cond = 1:2
        
        exp_starts = sync(2).timestamps(sync(2).values == 1);
        exp_stops = sync(2).timestamps(sync(2).values == 0);
        
        curr_exp_start = exp_starts(cond_expts{curr_cond}(1));
        curr_exp_stop = exp_stops(cond_expts{curr_cond}(end));
        
        curr_use_spikes = spike_times >= curr_exp_start & ...
            spike_times <= curr_exp_stop;
        
        spike_rate_cond(:,curr_cond) = ...
            accumarray(spike_templates(curr_use_spikes),1,[max(spike_templates),1])./ ...
            (curr_exp_stop - curr_exp_start);
        
    end
    
    spike_rate_change = (spike_rate_cond(:,2) - spike_rate_cond(:,1))./(spike_rate_cond(:,1)+spike_rate_cond(:,2));
    
    subplot(length(animal_days),3,(curr_animalday-1)*3+1,'YDir','reverse'); hold on;
    plot(spike_rate_cond(:,1),template_depths,'.k','MarkerSize',10);
    xlabel('Spikes/s')
    ylabel('Depth (\mum)');
    title({animal,day,'Pre-muscimol'});
    
    subplot(length(animal_days),3,(curr_animalday-1)*3+2,'YDir','reverse'); hold on;
    plot(spike_rate_cond(:,2),template_depths,'.k','MarkerSize',10);
    xlabel('Spikes/s')
    ylabel('Depth (\mum)');
    title({animal,day,'Post-muscimol'});
    
    subplot(length(animal_days),3,(curr_animalday-1)*3+3); hold on;
    plot(spike_rate_change,template_depths,'.k','MarkerSize',10);
    line([0,0],ylim);
    xlim([-1.1,1.1])
    set(gca,'YDir','reverse');
    xlabel('(Post-pre)/(pre+post)');
    ylabel('Depth (\mum)');
    title({animal,day,'Change'});
    
end

%% Striatal spike rate pre/post muscimol

animals = {'AP045','AP054','AP055','AP053','AP047','AP048'};

spike_rate_cond = cell(size(animals));

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworldNoRepeats';
    experiments = AP_find_experiments(animal,protocol);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
        
    disp(animal);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        
        % Load data (first experiment - but spikes throughout used)
        experiment = experiments(curr_day).experiment(1);
        load_parts.ephys = true;
        AP_load_experiment
        
        % Set experiments in conditions 
        % (1-3 pre-muscimol, 4+ post-muscimol)
        % (assumes all repeated expts/failures were post-muscimol)
        curr_experiments = AP_list_experiments(animal,day);
        cond_expts = {[1:3],[4:length(curr_experiments)]};
        
        for curr_cond = 1:2
                exp_starts = sync(2).timestamps(sync(2).values == 1);
                exp_stops = sync(2).timestamps(sync(2).values == 0);
                
                curr_exp_start = exp_starts(cond_expts{curr_cond}(1));
                curr_exp_stop = exp_stops(cond_expts{curr_cond}(end));
                
                curr_use_spikes = spike_times >= curr_exp_start & ...
                    spike_times <= curr_exp_stop & ~isnan(aligned_str_depth_group);
                
                spike_rate_cond{curr_animal}{curr_day}(:,curr_cond) = ...
                    accumarray(aligned_str_depth_group(curr_use_spikes),1, ...
                    [n_aligned_depths,1])./(curr_exp_stop - curr_exp_start);
        end
         
        clearvars -except animals curr_animal animal experiments curr_day ...
            spike_rate_cond
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
    
end

% Concatenate across recordings and get change
spike_rate_cond_cat = cell2mat(reshape([spike_rate_cond{:}],1,1,[]));
spike_rate_cond_cat_change = ...
    squeeze((spike_rate_cond_cat(:,1,:) - spike_rate_cond_cat(:,2,:))./ ...
    (spike_rate_cond_cat(:,1,:) + spike_rate_cond_cat(:,2,:)));
n_depths = size(spike_rate_cond_cat,1);

figure;
set(gca,'YDir','reverse'); hold on;
plot(spike_rate_cond_cat_change,1:n_depths,'color',[0.5,0.5,0.5]);
errorbar(nanmean(spike_rate_cond_cat_change,2),1:n_depths, ...
    AP_sem(spike_rate_cond_cat_change,2),'k','horizontal','linewidth',2);
line([0,0],ylim,'color','r')
xlim([-1.1,1.1]);
ylim([1-0.2,n_depths+0.2])
xlabel('(post-pre)/(post+pre)');
ylabel('Striatal depth');
title('Muscimol change');


%% Cortex > striatum domain kernels pre/post muscimol

protocols = {'vanillaChoiceworldNoRepeats_pre_muscimol','vanillaChoiceworldNoRepeats_post_muscimol'};

for protocol = protocols 
    protocol = cell2mat(protocol);
    
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
    k_fn = [data_path filesep 'ctx_str_kernels_' protocol];
    load(k_fn);
    
    framerate = 35;
    upsample_factor = 1;
    sample_rate = framerate*upsample_factor;
    kernel_t = [-0.5,0.5];
    kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
    t = kernel_frames/sample_rate;
    
    % Concatenate explained variance
    expl_var_animal = cell2mat(cellfun(@(x) nanmean(horzcat(x{:}),2),ctx_str_expl_var','uni',false));
    figure('Name',protocol);
    errorbar(nanmean(expl_var_animal,2),AP_sem(expl_var_animal,2),'k','linewidth',2);
    xlabel('Striatal depth');
    ylabel('Fraction explained variance');
    
    % Concatenate and mean
    % (kernel is -:+ fluorescence lag, flip to be spike-oriented)
    k_px_timeflipped = cellfun(@(x) cellfun(@(x) x(:,:,end:-1:1,:),x,'uni',false),ctx_str_kernel,'uni',false);
    k_px_animal = cellfun(@(x) nanmean(cat(5,x{:}),5),k_px_timeflipped,'uni',false);
    k_px = nanmean(double(cat(5,k_px_animal{:})),5);
    
    % Get center-of-mass maps
    k_px_positive = k_px;
    k_px_positive(k_px_positive < 0) = 0;
    k_px_com = sum(k_px_positive.*permute(1:n_aligned_depths,[1,3,4,2]),4)./sum(k_px_positive,4);
    k_px_com_colored = nan(size(k_px_com,1),size(k_px_com,2),3,size(k_px_com,3));
    
    use_colormap = min(jet(255)-0.2,1);
    for curr_frame = 1:size(k_px_com,3)
        k_px_com_colored(:,:,:,curr_frame) = ...
            ind2rgb(round(mat2gray(k_px_com(:,:,curr_frame),...
            [1,n_aligned_depths])*size(use_colormap,1)),use_colormap);
    end
    
    % Plot center kernel frames independently at t = 0
    figure('Name',protocol);
    plot_frame = kernel_frames == 0;
    for curr_depth = 1:n_aligned_depths
       subplot(n_aligned_depths,1,curr_depth);
       imagesc(k_px(:,:,plot_frame,curr_depth));
       AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
       axis image off;
       colormap(brewermap([],'*RdBu'));
       caxis([-0.01,0.01]);
    end
    
    % Plot center-of-mass color at select time points 
    plot_t = [-0.05:0.025:0.05];
    
    k_px_com_colored_t = ...
        permute(reshape(interp1(t,permute(reshape(k_px_com_colored,[],3,length(t)), ...
        [3,1,2]),plot_t),length(plot_t),size(k_px_com_colored,1), ...
        size(k_px_com_colored,2),3),[2,3,4,1]);
    
    k_px_max = squeeze(max(k_px,[],4));
    k_px_max_t = ...
        permute(reshape(interp1(t,reshape(k_px_max,[],length(t))', ...
        plot_t),length(plot_t),size(k_px_max,1), ...
        size(k_px_max,2)),[2,3,1]);
    
    weight_max = 0.005;
    figure('Name',protocol);
    for t_idx = 1:length(plot_t)
        subplot(1,length(plot_t),t_idx);
        p = image(k_px_com_colored_t(:,:,:,t_idx));
        set(p,'AlphaData', ...
            mat2gray(k_px_max_t(:,:,t_idx),[0,weight_max]));
        axis image off;
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        title([num2str(plot_t(t_idx)),' s']);
    end    
        
    % Plot movie of kernels
    AP_imscroll(reshape(permute(k_px,[1,4,2,3]),size(k_px,1)*size(k_px,4),size(k_px,2),length(t)),t);
    colormap(brewermap([],'*RdBu'));
    caxis([-max(caxis),max(caxis)]);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5],[],[size(k_px,1),size(k_px,2),size(k_px,4),1]);
    axis image off
    
    drawnow;
    
end

%% VFS pre/post musicmol

muscimol_wf_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\muscimol_wf.mat';
load(muscimol_wf_fn);

std_cat = vertcat(muscimol_wf.std);
vfs_cat = vertcat(muscimol_wf.vfs);

figure; 

subplot(1,2,1);
imagesc(nanmean(cat(3,vfs_cat{:,1}),3));
axis image off;
caxis([-1,1]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('VFS pre-muscimol');

subplot(1,2,2);
imagesc(nanmean(cat(3,vfs_cat{:,2}),3));
axis image off;
caxis([-1,1]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('VFS post-muscimol');



%% ~~~~~~~~ After preprocessing/saving trial activity ~~~~~~~~

%% Passive gratings pre/post muscimol

data_fns = { ...
    'trial_activity_AP_lcrGratingPassive_pre_muscimol', ...
    'trial_activity_AP_lcrGratingPassive_post_muscimol'};

stimIDs = cell(2,1);
mua_muscimol = cell(2,1);
mua_ctxpred_muscimol = cell(2,1);
fluor_muscimol = cell(2,1);
fluor_roi_muscimol = cell(2,1);

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Split by experiment
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    
    % Get trials with movement during stim to exclude
    wheel_thresh = 0.025;
    quiescent_trials = ~any(abs(wheel_allcat(:,t >= -0.5 & t <= 0.5)) > wheel_thresh,2);
    quiescent_trials_exp = mat2cell(quiescent_trials,trials_recording,1);
    
    % Get stim by experiment
    trial_stim_allcat_exp = mat2cell(trial_stim_allcat,trials_recording,1);
    
    % Get fluorescence in ROIs deconv but not baseline-subtracted
    fluor_deconv_exp = mat2cell(fluor_allcat_deconv,trials_recording,length(t),n_vs);
    fluor_roi_deconv_exp = mat2cell(fluor_roi_deconv,trials_recording,length(t),numel(wf_roi));
    mua_allcat_exp = mat2cell(mua_allcat,trials_recording,length(t),n_depths);
    mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,trials_recording,length(t),n_depths);
    
    % Grab data
    stimIDs{curr_data} = cellfun(@(data,trials) data(trials,:,:), ...
        trial_stim_allcat_exp,quiescent_trials_exp,'uni',false);
    
    mua_muscimol{curr_data} = cellfun(@(data,trials) data(trials,:,:), ...
        mua_allcat_exp,quiescent_trials_exp,'uni',false);
    
    mua_ctxpred_muscimol{curr_data} = cellfun(@(data,trials) data(trials,:,:), ...
        mua_ctxpred_allcat_exp,quiescent_trials_exp,'uni',false);
    
    fluor_muscimol{curr_data} = cellfun(@(data,trials) data(trials,:,:), ...
        fluor_deconv_exp,quiescent_trials_exp,'uni',false);
    
    fluor_roi_muscimol{curr_data} = cellfun(@(data,trials) data(trials,:,:), ...
        fluor_roi_deconv_exp,quiescent_trials_exp,'uni',false);
    
end

% Plot average fluorescence
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');

use_stim = 1;

fluor_premuscimol_mean = ...
    svdFrameReconstruct(U_master(:,:,1:200), ...
    permute(nanmean(cell2mat(cellfun(@(x,stim) ...
    nanmean(x(stim == use_stim,:,:),1),fluor_muscimol{1},stimIDs{1},'uni',false)),1),[3,2,1]));
fluor_postmuscimol_mean = ...
    svdFrameReconstruct(U_master(:,:,1:200), ...
    permute(nanmean(cell2mat(cellfun(@(x,stim) ...
    nanmean(x(stim == use_stim,:,:),1),fluor_muscimol{2},stimIDs{2},'uni',false)),1),[3,2,1]));

AP_imscroll([fluor_premuscimol_mean,fluor_postmuscimol_mean]);
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5],[],[size(U_master,1),size(U_master,2),1,2]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
axis image;
colormap(brewermap([],'*RdBu'));


figure;
t_stim = t >= 0.05 & t <= 0.15;

subplot(1,2,1)
imagesc(nanmean(fluor_premuscimol_mean(:,:,t_stim),3));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
c = caxis;
axis image off;
colormap(brewermap([],'*RdBu'));
title('Pre-muscimol');

subplot(1,2,2)
imagesc(nanmean(fluor_postmuscimol_mean(:,:,t_stim),3));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
caxis(c);
axis image off;
colormap(brewermap([],'*RdBu'));
title('Post-muscimol');

% Get pre/post stim response
use_stim = 1;

mua_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_muscimol{1},stimIDs{1},'uni',false));
mua_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_muscimol{2},stimIDs{2},'uni',false));

mua_ctxpred_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_ctxpred_muscimol{1},stimIDs{1},'uni',false));
mua_ctxpred_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_ctxpred_muscimol{2},stimIDs{2},'uni',false));

fluor_roi_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_roi_muscimol{1},stimIDs{1},'uni',false));
fluor_roi_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_roi_muscimol{2},stimIDs{2},'uni',false));

% Plot all str responses
figure;
for curr_str = 1:n_depths
    subplot(n_depths,1,curr_str);
    AP_errorfill(t,nanmean(mua_premuscimol_mean(:,:,curr_str),1), ...
        AP_sem(mua_premuscimol_mean(:,:,curr_str),1),'k',1,false);
    AP_errorfill(t,nanmean(mua_postmuscimol_mean(:,:,curr_str),1), ...
        AP_sem(mua_postmuscimol_mean(:,:,curr_str),1),'r',1,false);
    xlabel('Time from stim (s)');
    ylabel('Spikes (std)');
    title(['Str ' num2str(curr_str)]);
end
linkaxes(get(gcf,'Children'))

% Plot ctx responses
figure;
plot_ctx = [1,3,7];
for curr_ctx_idx = 1:length(plot_ctx)
    curr_ctx = plot_ctx(curr_ctx_idx);
    subplot(length(plot_ctx),1,curr_ctx_idx);
    AP_errorfill(t,nanmean(fluor_roi_premuscimol_mean(:,:,curr_ctx),1), ...
        AP_sem(fluor_roi_premuscimol_mean(:,:,curr_ctx),1),'k',1,false);
    AP_errorfill(t,nanmean(fluor_roi_postmuscimol_mean(:,:,curr_ctx),1), ...
        AP_sem(fluor_roi_postmuscimol_mean(:,:,curr_ctx),1),'r',1,false);
    ylabel(wf_roi(curr_ctx).area);
end
linkaxes(get(gcf,'Children'))

% Plot pre/post muscimol and repsonse change for pair of str/ctx
plot_str = 1;
plot_ctx = 3;

t_stim = t >= 0.05 & t <= 0.15;
mua_avg_premuscimol = permute(nanmean(mua_premuscimol_mean(:,t_stim,:),2),[1,3,2]);
mua_avg_postmuscimol = permute(nanmean(mua_postmuscimol_mean(:,t_stim,:),2),[1,3,2]);
mua_avg_postpre_change = (mua_avg_postmuscimol-mua_avg_premuscimol)./(mua_avg_premuscimol);

mua_ctxpred_avg_premuscimol = permute(nanmean(mua_ctxpred_premuscimol_mean(:,t_stim,:),2),[1,3,2]);
mua_ctxpred_avg_postmuscimol = permute(nanmean(mua_ctxpred_postmuscimol_mean(:,t_stim,:),2),[1,3,2]);
mua_ctxpred_avg_postpre_change = (mua_ctxpred_avg_postmuscimol-mua_ctxpred_avg_premuscimol)./(mua_ctxpred_avg_premuscimol);

fluor_avg_premuscimol = permute(nanmean(fluor_roi_premuscimol_mean(:,t_stim,:),2),[1,3,2]);
fluor_avg_postmuscimol = permute(nanmean(fluor_roi_postmuscimol_mean(:,t_stim,:),2),[1,3,2]);
fluor_avg_postpre_change = (fluor_avg_postmuscimol-fluor_avg_premuscimol)./(fluor_avg_premuscimol);

figure;
subplot(1,3,1);hold on;
AP_errorfill(t,nanmean(mua_premuscimol_mean(:,:,plot_str),1), ...
    AP_sem(mua_premuscimol_mean(:,:,plot_str),1),'k');
AP_errorfill(t,nanmean(mua_postmuscimol_mean(:,:,plot_str),1), ...
    AP_sem(mua_postmuscimol_mean(:,:,plot_str),1),'r');
xlim([-0.2,1])
xlabel('Time from stim (s)')
ylabel(['Str ' num2str(plot_str)]);
axis square

subplot(1,3,2);hold on;
AP_errorfill(t,nanmean(fluor_roi_premuscimol_mean(:,:,plot_ctx),1), ...
    AP_sem(fluor_roi_premuscimol_mean(:,:,plot_ctx),1),'k');
AP_errorfill(t,nanmean(fluor_roi_postmuscimol_mean(:,:,plot_ctx),1), ...
    AP_sem(fluor_roi_postmuscimol_mean(:,:,plot_ctx),1),'r');
xlim([-0.2,1])
xlabel('Time from stim (s)')
ylabel(wf_roi(plot_ctx).area);
axis square

subplot(1,3,3);
plot(fluor_avg_postpre_change(:,plot_ctx),mua_avg_postpre_change(:,plot_str),'.k','MarkerSize',20)
xlabel([wf_roi(plot_ctx).area  ' (post-pre)/pre']);
ylabel(['Str ' num2str(plot_str) ' (post-pre)/pre']);
line([-1,1],[-1,1],'color','k');
line([-1,1],[0,0],'color','k','linestyle','--');
line([0,0],[-1,1],'color','k','linestyle','--');
axis square;

nonan_points = ~isnan(mua_avg_postpre_change(:,plot_str)) & ...
    ~isnan(fluor_avg_postpre_change(:,plot_ctx));
[r,p] = corrcoef(fluor_avg_postpre_change(nonan_points,plot_ctx), ...
    mua_avg_postpre_change(nonan_points,plot_str));


%%%%%%% TESTING

% fluor_premuscimol_mean = cell2mat(permute(cellfun(@(x) ...
%     svdFrameReconstruct(U_master(:,:,1:200),permute(x,[3,2,1])), ...
%     cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1), ...
%     fluor_muscimol{1},stimIDs{1},'uni',false),'uni',false),[2,3,4,1]));
% 
% fluor_postmuscimol_mean = cell2mat(permute(cellfun(@(x) ...
%     svdFrameReconstruct(U_master(:,:,1:200),permute(x,[3,2,1])), ...
%     cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1), ...
%     fluor_muscimol{2},stimIDs{2},'uni',false),'uni',false),[2,3,4,1]));
% 
% fluor_avg_premuscimol = squeeze(nanmean(fluor_premuscimol_mean(:,:,t_stim,:),3));
% fluor_avg_postmuscimol = squeeze(nanmean(fluor_postmuscimol_mean(:,:,t_stim,:),3));

%% Striatal average task activity pre/post muscimol

data_fns = { ...
    'trial_activity_vanillaChoiceworldNoRepeats_pre_muscimol', ...
    'trial_activity_vanillaChoiceworldNoRepeats_post_muscimol'};

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Split data by recording
    trials_allcat = size(wheel_allcat,1);
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    use_split = trials_recording;
    split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
        [1:length(use_split)]',reshape(use_split,[],1),'uni',false));
    
    % Plot stim-aligned/sorted measured and predicted striatum activity
    % (correct contra trials)      
    switch curr_data
        case 1
            cond_name = 'Pre-muscimol';
        case 2
            cond_name = 'Post-muscimol';
    end
    
    for curr_trial_set = 1:2
        switch curr_trial_set
            case 1
                plot_trials = move_t < Inf & trial_stim_allcat > 0 & trial_choice_allcat == -1;
                figure('Name',[cond_name ': Correct contra trials']);
            case 2
                plot_trials = move_t < Inf & trial_stim_allcat < 0 & trial_choice_allcat == 1;
                figure('Name',[cond_name ': Correct ipsi trials']);
        end
        
        p = gobjects(n_depths,4);
        colormap(brewermap([],'Greys'));
        for curr_depth = 1:n_depths
            
            % Get trials to plot, sort by reaction time
            curr_trials = plot_trials & ~all(isnan(mua_allcat(:,:,curr_depth)),2);
            curr_trials_idx = find(curr_trials);
            [~,rxn_sort_idx] = sort(move_t(curr_trials_idx));
            
            sorted_plot_trials = curr_trials_idx(rxn_sort_idx);
            
            curr_plot = mua_allcat(sorted_plot_trials,:,curr_depth);
            curr_taskpred_plot = mua_taskpred_allcat(sorted_plot_trials,:,curr_depth);
            curr_ctxpred_plot = mua_ctxpred_allcat(sorted_plot_trials,:,curr_depth);
            
            % Smooth and plot with stim/move/reward times
            % (as conv(nans-zeroed)./conv(non-nan) to ignore in nans in conv)
            smooth_filt = [50,1]; % (trials x frames)
            
            curr_plot_smooth = conv2(curr_plot,ones(smooth_filt),'same')./ ...
                conv2(~isnan(curr_plot),ones(smooth_filt),'same');
            
            curr_taskpred_plot_smooth = curr_taskpred_plot;
            curr_taskpred_plot_smooth(isnan(curr_taskpred_plot_smooth)) = 0;
            curr_taskpred_plot_smooth = conv2(curr_taskpred_plot_smooth,ones(smooth_filt),'same')./ ...
                conv2(~isnan(curr_taskpred_plot),ones(smooth_filt),'same');
            
            curr_ctxpred_plot_smooth = curr_ctxpred_plot;
            curr_ctxpred_plot_smooth(isnan(curr_ctxpred_plot_smooth)) = 0;
            curr_ctxpred_plot_smooth = conv2(curr_ctxpred_plot_smooth,ones(smooth_filt),'same')./ ...
                conv2(~isnan(curr_ctxpred_plot),ones(smooth_filt),'same');
            
            p(curr_depth,1) = subplot(n_depths,4,1+(curr_depth-1)*4,'YDir','reverse'); hold on
            imagesc(t,[],curr_plot_smooth);
            plot(zeros(size(curr_trials_idx)),1:length(curr_trials_idx),'MarkerSize',1,'color','r');
            plot(move_t(sorted_plot_trials),1:length(curr_trials_idx),'MarkerSize',1,'color',[0.8,0,0.8]);
            %         plot(outcome_t(sorted_plot_trials),1:length(curr_trials_idx),'.','MarkerSize',1,'color','b');
            axis tight;
            xlim([-0.2,1]);
            xlabel('Time from stim');
            ylabel('Trials (rxn sorted)');
            title('Measured');
            
            p(curr_depth,2) = subplot(n_depths,4,2+(curr_depth-1)*4,'YDir','reverse'); hold on
            imagesc(t,[],curr_taskpred_plot_smooth);
            plot(zeros(size(curr_trials_idx)),1:length(curr_trials_idx),'MarkerSize',1,'color','r');
            plot(move_t(sorted_plot_trials),1:length(curr_trials_idx),'MarkerSize',1,'color',[0.8,0,0.8]);
            %         plot(outcome_t(sorted_plot_trials),1:length(curr_trials_idx),'.','MarkerSize',1,'color','b');
            axis tight;
            xlim([-0.2,1]);
            xlabel('Time from stim');
            ylabel('Trials (rxn sorted)');
            title('Task-predicted');
            
            p(curr_depth,3) = subplot(n_depths,4,3+(curr_depth-1)*4,'YDir','reverse'); hold on
            imagesc(t,[],curr_ctxpred_plot_smooth);
            plot(zeros(size(curr_trials_idx)),1:length(curr_trials_idx),'MarkerSize',1,'color','r');
            plot(move_t(sorted_plot_trials),1:length(curr_trials_idx),'MarkerSize',1,'color',[0.8,0,0.8]);
            %         plot(outcome_t(sorted_plot_trials),1:length(curr_trials_idx),'.','MarkerSize',1,'color','b');
            axis tight;
            xlim([-0.2,1]);
            xlabel('Time from stim');
            ylabel('Trials (rxn sorted)');
            title('Cortex-predicted');
            
            % Split and average trials by animal
            curr_trials_exp = mat2cell(curr_trials,use_split,1);
            
            curr_mua_exp = mat2cell(mua_allcat(:,:,curr_depth),use_split,length(t));
            curr_mua_exp_mean = cell2mat(cellfun(@(data,trials) nanmean(data(trials,:),1), ...
                curr_mua_exp,curr_trials_exp,'uni',false));
            
            curr_taskpred_mua_exp = mat2cell(mua_taskpred_allcat(:,:,curr_depth),use_split,length(t));
            curr_taskpred_mua_exp_mean = cell2mat(cellfun(@(data,trials) nanmean(data(trials,:),1), ...
                curr_taskpred_mua_exp,curr_trials_exp,'uni',false));
            
            curr_ctxpred_mua_exp = mat2cell(mua_ctxpred_allcat(:,:,curr_depth),use_split,length(t));
            curr_ctxpred_mua_exp_mean = cell2mat(cellfun(@(data,trials) nanmean(data(trials,:),1), ...
                curr_ctxpred_mua_exp,curr_trials_exp,'uni',false));
            
            % Plot PSTH (measured, task-predicted, cortex-predicted);
            p(curr_depth,4) = subplot(n_depths,4,4+(curr_depth-1)*4); hold on
            p1 = AP_errorfill(t,nanmean(curr_mua_exp_mean,1), ...
                AP_sem(curr_mua_exp_mean,1),'k',0.5);
            p2 = AP_errorfill(t,nanmean(curr_taskpred_mua_exp_mean,1), ...
                AP_sem(curr_taskpred_mua_exp_mean,1),[0,0,0.7],0.5);
            p3 = AP_errorfill(t,nanmean(curr_ctxpred_mua_exp_mean,1), ...
                AP_sem(curr_ctxpred_mua_exp_mean,1),[0,0.7,0],0.5);
            xlim([-0.2,1])
            line([0,0],ylim,'color','r');
            line(repmat(median(move_t(sorted_plot_trials)),1,2),ylim,'color',[0.8,0,0.8],'linestyle','--');
            line(repmat(median(outcome_t(sorted_plot_trials)),1,2),ylim,'color','b','linestyle','--');
            xlabel('Time from stim');
            ylabel('Spikes (std)');
            legend([p1,p2,p3],{'Measured','Task-predicted','Cortex-predicted'});
            
        end
        % Link the x-axes, set the c/y-axes same within a row
        linkaxes(p(:),'x');
        
        for curr_row = 1:size(p,1)
            curr_ylim = ylim(p(curr_row,4));
            caxis(p(curr_row,1),[0,curr_ylim(2)]);
            caxis(p(curr_row,2),[0,curr_ylim(2)]);
            caxis(p(curr_row,3),[0,curr_ylim(2)]);
        end
        
        trial_scale = 500;
        t_scale = 0.5;
        y_scale = 1;
        line(p(1,1),min(xlim(p(1,1))) + [0,t_scale],repmat(min(ylim(p(1,1))),2,1),'color','k','linewidth',3);
        line(p(1,4),min(xlim(p(1,4))) + [0,t_scale],repmat(min(ylim(p(1,4))),2,1),'color','k','linewidth',3);
        line(p(1,1),repmat(min(xlim(p(1,1))),2,1),min(ylim(p(1,1))) + [0,trial_scale],'color','k','linewidth',3);
        line(p(1,4),repmat(min(xlim(p(1,4))),2,1),min(ylim(p(1,4))) + [0,y_scale],'color','k','linewidth',3);
        drawnow;
        
    end
    
end


%% Get task stim kernels pre/post muscimol

data_fns = { ...
    'trial_activity_vanillaChoiceworldNoRepeats_pre_muscimol', ...
    'trial_activity_vanillaChoiceworldNoRepeats_post_muscimol'};

mua_norm_exp = cell(2,1);
task_str_kernel = cell(2,1);
task_str_ctxpred_kernel = cell(2,1);
task_ctx_roi_kernel = cell(2,1);

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Get average task > cortex kernels (ROIs)
    n_regressors = length(task_regressor_labels);
    regressor_roi = cell(n_regressors,1);
    for curr_regressor = 1:n_regressors
        
        curr_k = cell2mat(cellfun(@(x) x{curr_regressor}, ...
            permute(vertcat(fluor_taskpred_k_all{:}),[2,3,4,1]),'uni',false));
        
        curr_k_roi = nan(size(curr_k,1),size(curr_k,2),n_rois,size(curr_k,4));
        for curr_subregressor = 1:size(curr_k,1)
            for curr_exp = 1:size(curr_k,4)
                curr_k_roi(curr_subregressor,:,:,curr_exp) = ...
                    permute(AP_svd_roi(U_master(:,:,1:n_vs), ...
                    permute(curr_k(curr_subregressor,:,:,curr_exp),[3,2,1]), ...
                    [],[],cat(3,wf_roi.mask)),[3,2,1]);
            end
        end
        
        regressor_roi{curr_regressor} = curr_k_roi;
        
        AP_print_progress_fraction(curr_regressor,n_regressors);
    end
    
    % Keep task > str kernels, task > ctx-pred str norm factor
    mua_norm_exp{curr_data} = vertcat(mua_norm{:});
    task_str_kernel{curr_data} = vertcat(mua_taskpred_k_all{:});
    task_str_ctxpred_kernel{curr_data} = vertcat(mua_ctxpred_taskpred_k_all{:});
    task_ctx_roi_kernel{curr_data} = regressor_roi;
    
end

% Normalize and concatenate task kernels
mua_task_k_premuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_kernel{1},mua_norm_exp{1},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

mua_ctxpred_task_k_premuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_ctxpred_kernel{1},mua_norm_exp{1},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

mua_task_k_postmuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_kernel{2},mua_norm_exp{2},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

mua_ctxpred_task_k_postmuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_ctxpred_kernel{2},mua_norm_exp{2},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

% Get task>striatum parameters
n_regressors = length(task_regressor_labels);

% Plot task>striatum kernels
stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = [];
move_col = [0.6,0,0.6;0,0.6,0];
go_col = [0.8,0.8,0.2;0.5,0.5,0.5];
outcome_col = [0,0,0.8;0.8,0,0];
task_regressor_cols = {stim_col,move_col,go_col,outcome_col};
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

figure('Name','Pre-muscimol');
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);
        
        curr_kernels = mua_task_k_premuscimol{curr_regressor}(:,:,curr_depth,:);
        n_subregressors = size(mua_task_k_premuscimol{curr_regressor},1);
        col = task_regressor_cols{curr_regressor};
        for curr_subregressor = 1:n_subregressors
            AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
                nanmean(curr_kernels(curr_subregressor,:,:,:),4), ...
                AP_sem(curr_kernels(curr_subregressor,:,:,:),4), ...
                col(curr_subregressor,:));
        end
        
        xlabel('Time (s)');
        ylabel('Weight');
        title(task_regressor_labels{curr_regressor});
        line([0,0],ylim,'color','k');
        
    end
end
linkaxes(p);
y_scale = 1;
t_scale = 0.5;
line(min(xlim) + [0,t_scale],repmat(min(ylim),2,1),'color','k','linewidth',3);
line(repmat(min(xlim),2,1),min(ylim) + [0,y_scale],'color','k','linewidth',3);

figure('Name','Post-muscimol');
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);
        
        curr_kernels = mua_task_k_postmuscimol{curr_regressor}(:,:,curr_depth,:);
        n_subregressors = size(mua_task_k_postmuscimol{curr_regressor},1);
        col = task_regressor_cols{curr_regressor};
        for curr_subregressor = 1:n_subregressors
            AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
                nanmean(curr_kernels(curr_subregressor,:,:,:),4), ...
                AP_sem(curr_kernels(curr_subregressor,:,:,:),4), ...
                col(curr_subregressor,:),0.5);
        end
        
        xlabel('Time (s)');
        ylabel('Weight');
        title(task_regressor_labels{curr_regressor});
        line([0,0],ylim,'color','k');
        
    end
end
linkaxes(p);
y_scale = 1;
t_scale = 0.5;
line(min(xlim) + [0,t_scale],repmat(min(ylim),2,1),'color','k','linewidth',3);
line(repmat(min(xlim),2,1),min(ylim) + [0,y_scale],'color','k','linewidth',3);

% Plot kernel sums pre/post muscimol
str_k_sum_premuscimol = cellfun(@(x) permute(sum(x,2),[1,3,4,2]),mua_task_k_premuscimol,'uni',false);
str_k_sum_postmuscimol = cellfun(@(x) permute(sum(x,2),[1,3,4,2]),mua_task_k_postmuscimol,'uni',false);

figure;
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors
        if curr_regressor == 1
            x = unique([0.06,0.125,0.25,0.5,1].*[-1;1]);
        else
            x = 1:size(str_k_sum_premuscimol{curr_regressor},1);
        end
        
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);
        hold on
        
        errorbar(x,nanmean(str_k_sum_premuscimol{curr_regressor}(:,curr_depth,:),3), ...
            AP_sem(str_k_sum_premuscimol{curr_regressor}(:,curr_depth,:),3),'k','linewidth',2);
        
        errorbar(x,nanmean(str_k_sum_postmuscimol{curr_regressor}(:,curr_depth,:),3), ...
            AP_sem(str_k_sum_postmuscimol{curr_regressor}(:,curr_depth,:),3),'r','linewidth',2);
        
        axis tight
        xlim(xlim + [-0.2,0.2]);
        
        xlabel('Condition');
        ylabel('Weight sum');
        title(task_regressor_labels{curr_regressor});
        
    end
end

linkaxes(p,'y');




%% Cortex-explained striatal variance pre/post muscimol

data_fns = { ...
    'trial_activity_vanillaChoiceworldNoRepeats_pre_muscimol', ...
    'trial_activity_vanillaChoiceworldNoRepeats_post_muscimol'};

figure;
for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Split data by recording
    trials_allcat = size(wheel_allcat,1);
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    use_split = trials_recording;
    split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
        [1:length(use_split)]',reshape(use_split,[],1),'uni',false));
    
    % Get R^2 for task, cortex full, and cortex ROI predictions
    taskpred_r2 = nan(max(split_idx),n_depths);
    ctxpred_r2 = nan(max(split_idx),n_depths);
    ctxroipred_r2 = nan(max(split_idx),n_depths);
    for curr_exp = 1:max(split_idx)
        
        curr_data = reshape(permute(mua_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
        curr_taskpred_data = reshape(permute(mua_taskpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
        curr_ctxpred_data = reshape(permute(mua_ctxpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
        
        % Set common NaNs
        nan_samples = isnan(curr_data) | isnan(curr_taskpred_data) | isnan(curr_ctxpred_data);
        curr_data(nan_samples) = NaN;
        curr_taskpred_data(nan_samples) = NaN;
        curr_ctxpred_data(nan_samples) = NaN;
        
        taskpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_taskpred_data).^2,1)./ ...
            nansum((curr_data-nanmean(curr_data,1)).^2,1));
        ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
            nansum((curr_data-nanmean(curr_data,1)).^2,1));
    end
    
    subplot(2,1,curr_data);
    errorbar(nanmean(taskpred_r2,1),AP_sem(taskpred_r2,1),'b','linewidth',2,'CapSize',0);
    errorbar(nanmean(ctxpred_r2,1),AP_sem(ctxpred_r2,1),'color',[0,0.6,0],'linewidth',2,'CapSize',0);
    xlabel('Striatum depth');
    ylabel('Task explained variance');
    legend({'Task','Cortex'});
    
    
    
end


%% Offset between measured and predicted pre/post muscimol

data_fns = { ...
    'trial_activity_vanillaChoiceworldNoRepeats_pre_muscimol', ...
    'trial_activity_vanillaChoiceworldNoRepeats_post_muscimol'};

ctx_v_error = cell(2,2);
measured_pred_error_fig = figure;
measured_v_pred_fig = figure;
for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Split data by recording
    trials_allcat = size(wheel_allcat,1);
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    use_split = trials_recording;
    split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
        [1:length(use_split)]',reshape(use_split,[],1),'uni',false));
    
    % Set alignment shifts
    t_leeway = -t(1);
    leeway_samples = round(t_leeway*(sample_rate));
    stim_align = zeros(size(mua_allcat,1),1);
    move_align = -move_idx + leeway_samples;
    outcome_align = -outcome_idx + leeway_samples;
    
    % Set windows to average activity
    timeavg_labels = {'Stim','Move early','Move late','Outcome'};
    timeavg_t = {[0.05,0.15],[-0.05,0.05],[0.05,0.15],[0.05,0.15]};
    timeavg_align = {stim_align,move_align,move_align,outcome_align};
    timeavg_trial_conditions = ...
        {[trial_stim_allcat > 0,trial_stim_allcat == 0,trial_stim_allcat < 0], ...
        [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
        [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
        [trial_outcome_allcat == 1, trial_outcome_allcat == -1]};
    timeavg_condition_colors = ...
        {[1,0,0;0,0,0;0,0,1], ...
        [0.6,0,0.6;0,0.6,0], ...
        [0.6,0,0.6;0,0.6,0], ...
        [0,0,0.7;0,0,0]};
    timeavg_event_offset = {'Stim','Move onset','Move onset','Outcome'};
    
    % Set activity percentiles and bins
    act_prctile = [5,95];
    n_act_bins = 5;
    
    % Set areas and conditions
    plot_areas = [1,2,3,4];
    
    % Loop across area pairs, plot binned predicted v measured activity
    curr_act_allcat = mua_allcat;
    curr_act_pred_allcat = mua_ctxpred_allcat;

    % (old: to use all events)
    % task_fix = mua_taskpred_allcat - mua_ctxpred_taskpred_allcat;
    % curr_act_pred_fix_allcat = curr_act_pred_allcat + task_fix;
    
    % (to use individual events)
    task_fix = (mua_taskpred_allcat - mua_taskpred_reduced_allcat) - ...
        (mua_ctxpred_taskpred_allcat - mua_ctxpred_taskpred_reduced_allcat);
    curr_act_pred_fix_allcat = curr_act_pred_allcat + task_fix;
    
    for curr_area_idx = 1:length(plot_areas)
        
        plot_area = plot_areas(curr_area_idx);
        
        for curr_timeavg = 1:length(timeavg_labels)
            
            trial_conditions = timeavg_trial_conditions{curr_timeavg};
            curr_event_offset_idx = find(strcmp(timeavg_event_offset{curr_timeavg},task_regressor_labels));
            curr_act_pred_fix_allcat_event = ...
                curr_act_pred_fix_allcat(:,:,:,curr_event_offset_idx);
            
            % (re-align and split activity)
            act_title = timeavg_labels{curr_timeavg};
            
            curr_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_act_allcat(trial,:,plot_area), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_allcat,1)),'uni',false)), ...
                use_split,length(t));
            
            curr_act_pred = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_act_pred_allcat(trial,:,plot_area), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
                use_split,length(t));
            
            curr_act_predfix = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_act_pred_fix_allcat_event(trial,:,plot_area), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_fix_allcat_event,1)),'uni',false)), ...
                use_split,length(t));
            
            trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
            
            % (get average activity within window)
            curr_event_t = t >= timeavg_t{curr_timeavg}(1) & t <= timeavg_t{curr_timeavg}(2);
            curr_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act,'uni',false);
            curr_act_pred_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_pred,'uni',false);
            curr_act_predfix_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_predfix,'uni',false);

            %%% Plot total error in prediction           
            curr_meas_pred_error = cellfun(@(x,y) x-y,curr_act_avg,curr_act_pred_avg,'uni',false);
            curr_meas_pred_error_cond = cell2mat(cellfun(@(x,cond) arrayfun(@(curr_cond) ...
                nanmean(x(cond(:,curr_cond))),1:size(cond,2)),curr_meas_pred_error,trial_conditions_exp,'uni',false));
                   
            figure(measured_pred_error_fig);
            
            subplot(length(plot_areas),length(timeavg_labels), ...
                sub2ind(fliplr([length(plot_areas),length(timeavg_labels)]),curr_timeavg,curr_area_idx));
            hold on;
            
            switch curr_data
                case 1
                    col = 'k';
                case 2
                    col = 'r';
            end
            
            errorbar(nanmean(curr_meas_pred_error_cond,1),AP_sem(curr_meas_pred_error_cond,1),col,'linewidth',2);
            ylabel('Measured-Predicted');
            xlabel('Condition');
            
            drawnow;
            
            %%% Get and plot offset for Str 1 stim only
            
            if curr_area_idx == 1 &&  curr_timeavg == 1
                                
                % (bin predicted data across percentile range)
                pred_bin_edges = prctile(cell2mat(curr_act_pred_avg),linspace(act_prctile(1),act_prctile(2),n_act_bins+1));
                pred_bin_centers = pred_bin_edges(1:end-1) + diff(pred_bin_edges)./2;
                pred_trial_bins = cellfun(@(x) discretize(x,pred_bin_edges),curr_act_pred_avg,'uni',false);
                
                % (get activity binned by predicted)
                pred_use_trials = cellfun(@(act,act_pred,trial_bins) ...
                    squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
                    squeeze(~any(isnan(act_pred(:,curr_event_t)),2)) & ...
                    ~isnan(trial_bins), ...
                    curr_act,curr_act_pred,pred_trial_bins,'uni',false);
                
                act_predbinmean = cell2mat(arrayfun(@(condition) ...
                    cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                    accumarray(bins(use_trials & trial_cond(:,condition)), ...
                    act(use_trials & trial_cond(:,condition)), ...
                    [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                    curr_act_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
                    permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
                
                act_pred_predbinmean = cell2mat(arrayfun(@(condition) ...
                    cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                    accumarray(bins(use_trials & trial_cond(:,condition)), ...
                    act(use_trials & trial_cond(:,condition)), ...
                    [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                    curr_act_pred_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
                    permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
                
                % (get "fixed" predicted activity binned by predicted)
                act_predfix_predbinmean = cell2mat(arrayfun(@(condition) ...
                    cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                    accumarray(bins(use_trials & trial_cond(:,condition)), ...
                    act(use_trials & trial_cond(:,condition)), ...
                    [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                    curr_act_predfix_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
                    permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
                
                % Get condition difference significance from shuffle
                n_shuff = 1000;
                trial_conditions_shuff = cellfun(@(x) AP_shake(x,1), ...
                    repmat(trial_conditions_exp,1,n_shuff),'uni',false);
                
                act_predbinmean_condshuff = cell2mat(arrayfun(@(curr_shuff) ...
                    cell2mat(arrayfun(@(condition) ...
                    cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                    accumarray(bins(use_trials & trial_cond(:,condition)), ...
                    act(use_trials & trial_cond(:,condition)), ...
                    [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                    curr_act_avg,pred_trial_bins,trial_conditions_shuff(:,curr_shuff), ...
                    pred_use_trials,'uni',false)'), ...
                    permute(1:size(trial_conditions,2),[1,3,2]),'uni',false)), ...
                    permute(1:n_shuff,[1,3,4,2]),'uni',false));
                
                act_predfix_predbinmean_condshuff = cell2mat(arrayfun(@(curr_shuff) ...
                    cell2mat(arrayfun(@(condition) ...
                    cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                    accumarray(bins(use_trials & trial_cond(:,condition)), ...
                    act(use_trials & trial_cond(:,condition)), ...
                    [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                    curr_act_predfix_avg,pred_trial_bins,trial_conditions_shuff(:,curr_shuff), ...
                    pred_use_trials,'uni',false)'), ...
                    permute(1:size(trial_conditions,2),[1,3,2]),'uni',false)), ...
                    permute(1:n_shuff,[1,3,4,2]),'uni',false));
                
                % (measured data: null = no difference between conditions)
                cond_combos = nchoosek(1:size(trial_conditions,2),2);
                cond_sig_diff = false(n_act_bins,size(cond_combos,1));
                predfix_cond_sig_diff = false(n_act_bins,size(cond_combos,1));
                for curr_cond_combo = 1:size(cond_combos,1)
                    curr_combo_diff = nanmean(act_predbinmean(:,:,cond_combos(curr_cond_combo,1)) - ...
                        act_predbinmean(:,:,cond_combos(curr_cond_combo,2)),2);
                    shuff_prctile = squeeze(prctile(nanmean(act_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,1),:) - ...
                        act_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,2),:),2),[2.5,97.5],4));
                    cond_sig_diff(:,curr_cond_combo) = ...
                        curr_combo_diff < shuff_prctile(:,1) | ...
                        curr_combo_diff > shuff_prctile(:,2);
                    
                    predfix_curr_combo_diff = nanmean(act_predfix_predbinmean(:,:,cond_combos(curr_cond_combo,1)) - ...
                        act_predfix_predbinmean(:,:,cond_combos(curr_cond_combo,2)),2);
                    predfix_shuff_prctile = squeeze(prctile(nanmean(act_predfix_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,1),:) - ...
                        act_predfix_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,2),:),2),[2.5,97.5],4));
                    predfix_cond_sig_diff(:,curr_cond_combo) = ...
                        predfix_curr_combo_diff < predfix_shuff_prctile(:,1) | ...
                        predfix_curr_combo_diff > predfix_shuff_prctile(:,2);
                end
                
                % Plot binned measured, predicted, and error (by predicted bins)
                measured_pred_fig = figure('color','w','Name', ...
                    ['Str ' num2str(plot_area)' ', ' timeavg_labels{curr_timeavg}]);
                n_col_bins = n_act_bins + 2;
                
                [binned_act_pred_t,binned_act_pred_grp] = grpstats(cell2mat(curr_act_pred), ...
                    [cell2mat(pred_trial_bins),trial_conditions],{'nanmean','gname'});
                binned_act_pred_grp = cellfun(@str2num,binned_act_pred_grp);
                
                [binned_act_t,binned_act_grp] = grpstats(cell2mat(curr_act), ...
                    [cell2mat(pred_trial_bins),trial_conditions],{'nanmean','gname'});
                binned_act_grp = cellfun(@str2num,binned_act_grp);
                
                % (plot predicted data timecourse)
                for curr_cond = 1:size(trial_conditions,2)
                    subplot(2,size(trial_conditions,2), ...
                        sub2ind(fliplr([2,size(trial_conditions,2)]),curr_cond,1)); hold on;
                    
                    curr_mean = nanmean(cell2mat(cellfun(@(act,use_trials,cond) ...
                        act(use_trials & cond(:,curr_cond),:), ...
                        curr_act_pred,pred_use_trials,trial_conditions_exp,'uni',false)),1);
                    curr_std = nanstd(cell2mat(cellfun(@(act,use_trials,cond) ...
                        act(use_trials & cond(:,curr_cond),:), ...
                        curr_act_pred,pred_use_trials,trial_conditions_exp,'uni',false)),[],1);
                    
                    AP_errorfill(t,curr_mean,curr_std,'k',0.5,true);
                    
                    set(gca,'ColorOrder',brewermap(n_col_bins,'*Greens'));
                    plot(t,binned_act_pred_t(binned_act_pred_grp(:,curr_cond + 1) == 1,:)','linewidth',2);
                    xlabel('Time'); ylabel('Predicted data');
                    title(['Condition ' num2str(curr_cond)]);
                    line(repmat(timeavg_t{curr_timeavg}(1),2,1),ylim,'color','k');
                    line(repmat(timeavg_t{curr_timeavg}(2),2,1),ylim,'color','k');
                end
                
                % (plot measured data timecourse)
                for curr_cond = 1:size(trial_conditions,2)
                    subplot(2,size(trial_conditions,2), ...
                        sub2ind(fliplr([2,size(trial_conditions,2)]),curr_cond,2)); hold on;
                    set(gca,'ColorOrder',[brewermap(n_col_bins,'*Greys')]);
                    plot(t,binned_act_t(binned_act_grp(:,curr_cond + 1) == 1,:)','linewidth',2);
                    xlabel('Time'); ylabel('Measured data');
                    title(['Condition ' num2str(curr_cond)]);
                    line(repmat(timeavg_t{curr_timeavg}(1),2,1),ylim,'color','k');
                    line(repmat(timeavg_t{curr_timeavg}(2),2,1),ylim,'color','k');
                end
                
                linkaxes(get(measured_pred_fig,'Children'),'xy');
                
                % Plot measured v binned predicted
                figure(measured_v_pred_fig);
                subplot(1,2,curr_data); hold on;
                set(gca,'ColorOrder',timeavg_condition_colors{curr_timeavg});
                
                fill_cols = min(timeavg_condition_colors{curr_timeavg} + 0.5,1);
                for curr_cond = 1:size(trial_conditions,2)
                    AP_errorfill( ...
                        squeeze(nanmean(act_pred_predbinmean(:,:,curr_cond),2)), ...
                        squeeze(nanmean(act_predfix_predbinmean(:,:,curr_cond),2)), ...
                        squeeze(AP_sem(act_predfix_predbinmean(:,:,curr_cond),2)), ...
                        fill_cols(curr_cond,:),1,false);
                end
                
                errorbar( ...
                    squeeze(nanmean(act_pred_predbinmean,2)), ...
                    squeeze(nanmean(act_predbinmean,2)), ...
                    squeeze(AP_sem(act_predbinmean,2)), ...
                    'linestyle','none','linewidth',2,'CapSize',10);
                
                ctx_col = brewermap(n_col_bins,'*Greens');
                scatter( ...
                    reshape(squeeze(nanmean(act_pred_predbinmean,2)),[],1), ...
                    reshape(squeeze(nanmean(act_predbinmean,2)),[],1), ...
                    60,repmat(ctx_col(1:n_act_bins,:),size(trial_conditions,2),1),'filled');
                
                xlabel(['Predicted (' num2str(plot_area) ')']);
                ylabel(['Measured (' num2str(plot_area) ')'])
                axis square tight;
                xlim(xlim + [-0.1,0.1]);
                title(act_title);
                
                % (plot significant measured condition differences)
                % (* and o = significant in both measured and "fixed" predicted)
                curr_ylim = max(ylim);
                for curr_cond_combo = 1:size(cond_combos,1)
                    % (plot * for measured condition differences)
                    sig_y = curr_ylim + 0.1*curr_cond_combo;
                    sig_x = pred_bin_centers;
                    if any(cond_sig_diff(:,curr_cond_combo))
                        plot(sig_x(cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                            '*','MarkerSize',10,'color', ...
                            timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,1),:));
                        plot(sig_x(cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                            '*','MarkerSize',5,'color', ...
                            timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,2),:));
                    end
                    % (plot o for [predicted condition differences)
                    if any(predfix_cond_sig_diff(:,curr_cond_combo))
                        plot(sig_x(predfix_cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                            'o','MarkerSize',15,'color', ...
                            timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,1),:));
                        plot(sig_x(predfix_cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                            'o','MarkerSize',10,'color', ...
                            timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,2),:));
                    end
                end   
                drawnow;
                
                %%%%%%%%%% TESTING %%%%%%%%%%%%%
                % Store the cortical activity and error
%                 use_cond = 1;                
%                 ctx_v_error{curr_data,1} = ...
%                     cellfun(@(x,cond) nanmean(x(cond(:,use_cond)),1),curr_act_ctx_avg,trial_conditions_exp);
%                 ctx_v_error{curr_data,2} = ...
%                     cellfun(@(x,cond) nanmean(x(cond(:,use_cond)),1),curr_meas_pred_error,trial_conditions_exp);
                % Store error but in bins
                use_cond = 1;     
                curr_offset = nanmean(act_predbinmean(:,:,1) - act_predbinmean(:,:,3),1);
                ctx_v_error{curr_data,1} = ...
                    cellfun(@(x,cond,pred_use_trials) ...
                    nanmean(x(pred_use_trials & cond(:,use_cond)),1), ...
                    curr_act_pred_avg,trial_conditions_exp,pred_use_trials);
                ctx_v_error{curr_data,2} = ...
                    curr_offset';
                
            end
        end
    end
end


% Plot striatum/cortex offset relative to cortical activity
ctx_change = (ctx_v_error{2,1} - ctx_v_error{1,1})./(abs(ctx_v_error{2,1}) + abs(ctx_v_error{1,1}));
error_change = (ctx_v_error{2,2} - ctx_v_error{1,2})./(abs(ctx_v_error{2,2}) + abs(ctx_v_error{1,2}));
figure; 
plot(ctx_change,error_change,'.k','MarkerSize',20);
line([-1,1],[0,0],'color','k','linestyle','--');
line([0,0],[-1,1],'color','k','linestyle','--');
line([-1,1],[-1,1],'color','k');
xlabel('Cortex (pred striatum) change');
ylabel('Error change');

nonan_points = ~isnan(ctx_change) & ~isnan(error_change);
[r,p] = corrcoef(ctx_change(nonan_points),error_change(nonan_points));
disp(r(2));
disp(p(2));

% HAVEN'T GOTTEN ABOVE TO LOOK RIGHT - NOT SURE WHAT'S RIGHT TO DO HERE
% above is biased based on which conditions mouse was shown, maybe just use
% std or something?


%% Cortex/striatum offset pre/post muscimol (passive)

data_fns = { ...
    'trial_activity_AP_lcrGratingPassive_pre_muscimol', ...
    'trial_activity_AP_lcrGratingPassive_post_muscimol'};

trial_cond_all = cell(2,1);
trial_avg_all = cell(2,1);
trial_avg_pred_all = cell(2,1);

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Split data by recording
    trials_allcat = size(wheel_allcat,1);
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    use_split = trials_recording;
    split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
        [1:length(use_split)]',reshape(use_split,[],1),'uni',false));
       
    % Set alignment shifts
    t_leeway = -t(1);
    leeway_samples = round(t_leeway*(sample_rate));
    stim_align = zeros(size(mua_allcat,1),1);
    
    % Get trials with movement during stim to exclude
    wheel_thresh = 0.025;
    quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > wheel_thresh,2);
    
    % Set windows to average activity
    timeavg_labels = {'Stim'};
    timeavg_t = {[0.05,0.15]};
    timeavg_align = {stim_align};
    timeavg_trial_conditions = ...
        {[trial_stim_allcat > 0 & quiescent_trials, ...
        trial_stim_allcat < 0 & quiescent_trials]};
    timeavg_condition_colors = ...
        {[1,0,0;0,0,1]};
    
    % Set activity percentiles and bins
    act_prctile = [5,95];
    n_act_bins = 5;
    
    % Set areas and conditions
    plot_areas = [1];
    
    % Loop across area pairs, plot binned predicted v measured activity
    curr_act_allcat = mua_allcat;
    
    % (ctx-predicted)
    curr_act_pred_allcat = mua_ctxpred_allcat;
        
    % (fix by average stim response within experiment)
    mua_allcat_exp = mat2cell(mua_allcat,use_split,length(t),n_depths);
    mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,use_split,length(t),n_depths);
    trial_contrastside_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
    trial_conditions_exp = mat2cell(timeavg_trial_conditions{1},use_split, ...
        size(timeavg_trial_conditions{1},2));
    
    str_stim_exp = nan(length(unique(trial_stim_allcat)),length(t),n_depths,length(use_split));
    ctx_stim_exp = nan(length(unique(trial_stim_allcat)),length(t),n_depths,length(use_split));
    curr_act_pred_fix_allcat_exp = mua_ctxpred_allcat_exp;
    for curr_exp = 1:length(use_split)
        
        curr_stim = unique(trial_contrastside_allcat_exp{curr_exp});
        for curr_stim_idx = 1:length(curr_stim)
            curr_trials = any(trial_conditions_exp{curr_exp},2) & ...
                trial_contrastside_allcat_exp{curr_exp} == curr_stim(curr_stim_idx);
            curr_act_stim_avg = nanmean(mua_allcat_exp{curr_exp}(curr_trials,:,:),1);
            curr_act_pred_stim_avg = nanmean(mua_ctxpred_allcat_exp{curr_exp}(curr_trials,:,:),1);
            
            str_stim_exp(curr_stim_idx,:,:,curr_exp) = curr_act_stim_avg;
            ctx_stim_exp(curr_stim_idx,:,:,curr_exp) = curr_act_pred_stim_avg;
            
            curr_stim_fix = curr_act_stim_avg - curr_act_pred_stim_avg;
            curr_act_pred_fix_allcat_exp{curr_exp}(curr_trials,:,:) = ...
                curr_act_pred_fix_allcat_exp{curr_exp}(curr_trials,:,:) + curr_stim_fix;
        end
    end
    curr_act_pred_fix_allcat = cell2mat(curr_act_pred_fix_allcat_exp);
    
    % Plot stim response
    stim_col = colormap_BlueWhiteRed(5);
    stim_col(6,:) = 0.5;
    stim_col_contrastside = unique([0,0.06,0.125,0.25,0.5,1].*[-1;1]);
    [~,used_stim_idx] = ismember(unique(trial_stim_allcat), ...
        stim_col_contrastside,'rows');
    
    str_stim_avg = nanmean(str_stim_exp,4);
    figure('color','w');
    for curr_area_idx = 1:length(plot_areas)
        curr_area = plot_areas(curr_area_idx);
        subplot(length(plot_areas),2,length(plot_areas)*(curr_area_idx-1)+1);
        hold on;
        for curr_stim_idx = 1:size(str_stim_avg,1)
            AP_errorfill(t, ...
                nanmean(str_stim_exp(curr_stim_idx,:,curr_area,:),4), ...
                AP_sem(str_stim_exp(curr_stim_idx,:,curr_area,:),4), ...
                stim_col(used_stim_idx(curr_stim_idx),:),0.5,true);
        end
        line([0,0],ylim,'color','k');
        ylabel(['Measured (' num2str(curr_area) ')']);
        xlabel('Time from stim');
        
        subplot(length(plot_areas),2,length(plot_areas)*(curr_area_idx-1)+2);
        hold on;
        for curr_stim_idx = 1:size(str_stim_avg,1)
            AP_errorfill(t, ...
                nanmean(ctx_stim_exp(curr_stim_idx,:,curr_area,:),4), ...
                AP_sem(ctx_stim_exp(curr_stim_idx,:,curr_area,:),4), ...
                stim_col(used_stim_idx(curr_stim_idx),:),0.5,true);
        end
        line([0,0],ylim,'color','k');
        ylabel(['Predicted (' num2str(curr_area) ')']);
        xlabel('Time from stim');
    end
    linkaxes(get(gcf,'Children'),'xy');
    
    measured_v_pred_fig = figure('color','w');
    for curr_area_idx = 1:length(plot_areas)
        
        plot_area = plot_areas(curr_area_idx);
        
        for curr_timeavg = 1:length(timeavg_labels)
            
            trial_conditions = timeavg_trial_conditions{curr_timeavg};
            
            % (re-align and split activity)
            act_title = timeavg_labels{curr_timeavg};
            
            curr_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_act_allcat(trial,:,plot_area), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_allcat,1)),'uni',false)), ...
                use_split,length(t));
            
            curr_act_pred = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_act_pred_allcat(trial,:,plot_area), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
                use_split,length(t));
            
            curr_act_predfix = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_act_pred_fix_allcat(trial,:,plot_area), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_fix_allcat,1)),'uni',false)), ...
                use_split,length(t));
            
            trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
            
            % (get average activity within window)
            curr_event_t = t >= timeavg_t{curr_timeavg}(1) & t <= timeavg_t{curr_timeavg}(2);
            curr_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act,'uni',false);
            curr_act_pred_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_pred,'uni',false);
            curr_act_predfix_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_predfix,'uni',false);
            
            % (store average activity and conditions)
            trial_cond_all{curr_data} = trial_conditions_exp;
            trial_avg_all{curr_data} = curr_act_avg;
            trial_avg_pred_all{curr_data} = curr_act_pred_avg;
            
            % (bin predicted data across percentile range)
            pred_bin_edges = prctile(cell2mat(curr_act_pred_avg),linspace(act_prctile(1),act_prctile(2),n_act_bins+1));
            pred_bin_centers = pred_bin_edges(1:end-1) + diff(pred_bin_edges)./2;
            pred_trial_bins = cellfun(@(x) discretize(x,pred_bin_edges),curr_act_pred_avg,'uni',false);
            
            % (get activity binned by predicted)
            pred_use_trials = cellfun(@(act,act_pred,trial_bins) ...
                squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
                squeeze(~any(isnan(act_pred(:,curr_event_t)),2)) & ...
                ~isnan(trial_bins), ...
                curr_act,curr_act_pred,pred_trial_bins,'uni',false);
            
            act_predbinmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                accumarray(bins(use_trials & trial_cond(:,condition)), ...
                act(use_trials & trial_cond(:,condition)), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                curr_act_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
                permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
            
            act_pred_predbinmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                accumarray(bins(use_trials & trial_cond(:,condition)), ...
                act(use_trials & trial_cond(:,condition)), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                curr_act_pred_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
                permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
            
            % (get "fixed" predicted activity binned by predicted)
            act_predfix_predbinmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                accumarray(bins(use_trials & trial_cond(:,condition)), ...
                act(use_trials & trial_cond(:,condition)), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                curr_act_predfix_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
                permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
            
            % Get condition difference significance from shuffle
            n_shuff = 1000;
            trial_conditions_shuff = cellfun(@(x) AP_shake(x,1), ...
                repmat(trial_conditions_exp,1,n_shuff),'uni',false);
            
            act_predbinmean_condshuff = cell2mat(arrayfun(@(curr_shuff) ...
                cell2mat(arrayfun(@(condition) ...
                cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                accumarray(bins(use_trials & trial_cond(:,condition)), ...
                act(use_trials & trial_cond(:,condition)), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                curr_act_avg,pred_trial_bins,trial_conditions_shuff(:,curr_shuff), ...
                pred_use_trials,'uni',false)'), ...
                permute(1:size(trial_conditions,2),[1,3,2]),'uni',false)), ...
                permute(1:n_shuff,[1,3,4,2]),'uni',false));
            
            act_predfix_predbinmean_condshuff = cell2mat(arrayfun(@(curr_shuff) ...
                cell2mat(arrayfun(@(condition) ...
                cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                accumarray(bins(use_trials & trial_cond(:,condition)), ...
                act(use_trials & trial_cond(:,condition)), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                curr_act_predfix_avg,pred_trial_bins,trial_conditions_shuff(:,curr_shuff), ...
                pred_use_trials,'uni',false)'), ...
                permute(1:size(trial_conditions,2),[1,3,2]),'uni',false)), ...
                permute(1:n_shuff,[1,3,4,2]),'uni',false));
            
            % (measured data: null = no difference between conditions)
            cond_combos = nchoosek(1:size(trial_conditions,2),2);
            cond_sig_diff = false(n_act_bins,size(cond_combos,1));
            predfix_cond_sig_diff = false(n_act_bins,size(cond_combos,1));
            for curr_cond_combo = 1:size(cond_combos,1)
                curr_combo_diff = nanmean(act_predbinmean(:,:,cond_combos(curr_cond_combo,1)) - ...
                    act_predbinmean(:,:,cond_combos(curr_cond_combo,2)),2);
                shuff_prctile = squeeze(prctile(nanmean(act_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,1),:) - ...
                    act_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,2),:),2),[2.5,97.5],4));
                cond_sig_diff(:,curr_cond_combo) = ...
                    curr_combo_diff < shuff_prctile(:,1) | ...
                    curr_combo_diff > shuff_prctile(:,2);
                
                predfix_curr_combo_diff = nanmean(act_predfix_predbinmean(:,:,cond_combos(curr_cond_combo,1)) - ...
                    act_predfix_predbinmean(:,:,cond_combos(curr_cond_combo,2)),2);
                predfix_shuff_prctile = squeeze(prctile(nanmean(act_predfix_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,1),:) - ...
                    act_predfix_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,2),:),2),[2.5,97.5],4));
                predfix_cond_sig_diff(:,curr_cond_combo) = ...
                    predfix_curr_combo_diff < predfix_shuff_prctile(:,1) | ...
                    predfix_curr_combo_diff > predfix_shuff_prctile(:,2);
            end
            
            % Plot measured v predicted in bins
            figure(measured_v_pred_fig)
            
            % (measured vs binned predicted)
            subplot(length(plot_areas),length(timeavg_labels), ...
                sub2ind(fliplr([length(plot_areas),length(timeavg_labels)]),curr_timeavg,curr_area_idx));
            hold on;
            set(gca,'ColorOrder',timeavg_condition_colors{curr_timeavg});
            
            % Plot "fixed" curves based on average stim difference
%             fill_cols = min(timeavg_condition_colors{curr_timeavg} + 0.5,1);
%             for curr_cond = 1:size(trial_conditions,2)
%                 AP_errorfill( ...
%                     squeeze(nanmean(act_pred_predbinmean(:,:,curr_cond),2)), ...
%                     squeeze(nanmean(act_predfix_predbinmean(:,:,curr_cond),2)), ...
%                     squeeze(AP_sem(act_predfix_predbinmean(:,:,curr_cond),2)), ...
%                     fill_cols(curr_cond,:),1,false);
%             end
            
            % Plot real curves
            errorbar( ...
                squeeze(nanmean(act_pred_predbinmean,2)), ...
                squeeze(nanmean(act_predbinmean,2)), ...
                squeeze(AP_sem(act_predbinmean,2)), ...
                'linestyle','-','linewidth',2,'CapSize',10);
            
            n_col_bins = n_act_bins + 2;
            ctx_col = brewermap(n_col_bins,'*Greens');
            scatter( ...
                reshape(squeeze(nanmean(act_pred_predbinmean,2)),[],1), ...
                reshape(squeeze(nanmean(act_predbinmean,2)),[],1), ...
                60,repmat(ctx_col(1:n_act_bins,:),size(trial_conditions,2),1),'filled');
            
            xlabel(['Predicted (' num2str(plot_area) ')']);
            ylabel(['Measured (' num2str(plot_area) ')'])
            axis square;
            title(act_title);
            
            % (plot significant measured condition differences)
            % (* and o = significant in both measured and "fixed" predicted)
            curr_ylim = max(ylim);
            for curr_cond_combo = 1:size(cond_combos,1)
                % (plot * for measured condition differences)
                sig_y = curr_ylim + 0.1*curr_cond_combo;
                sig_x = pred_bin_centers;
                if any(cond_sig_diff(:,curr_cond_combo))
                    plot(sig_x(cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                        '*','MarkerSize',10,'color', ...
                        timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,1),:));
                    plot(sig_x(cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                        '*','MarkerSize',5,'color', ...
                        timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,2),:));
                end
                % (plot o for [predicted condition differences)
                if any(predfix_cond_sig_diff(:,curr_cond_combo))
                    plot(sig_x(predfix_cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                        'o','MarkerSize',15,'color', ...
                        timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,1),:));
                    plot(sig_x(predfix_cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                        'o','MarkerSize',10,'color', ...
                        timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,2),:));
                end
            end
            
            drawnow;
            
        end
    end
    
    linkaxes(get(measured_v_pred_fig,'Children'),'xy');
    
end

% (testing out here)

trial_cond_all;
trial_avg_all;
trial_avg_pred_all;

% average activity? [exp,stim,muscimol]
a = cellfun(@(act,cond) ...
    cellfun(@(act,cond) [nanmean(act(cond(:,1))),nanmean(act(cond(:,2)))], ...
    act,cond,'uni',false), ...
    trial_avg_all,trial_cond_all,'uni',false);
a2 = cell2mat(permute(horzcat(a{:}),[1,3,2]));

b = cellfun(@(act,cond) ...
    cellfun(@(act,cond) [nanmean(act(cond(:,1))),nanmean(act(cond(:,2)))], ...
    act,cond,'uni',false), ...
    trial_avg_pred_all,trial_cond_all,'uni',false);
b2 = cell2mat(permute(horzcat(b{:}),[1,3,2]));

% testing cat across trials?

a = cell2mat(vertcat(trial_avg_all{1}));
b = cell2mat(vertcat(trial_avg_pred_all{1}));
c = cell2mat(vertcat(trial_cond_all{1}));



















