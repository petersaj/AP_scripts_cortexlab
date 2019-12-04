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

animal = 'AP047';
protocol = 'vanillaChoiceworld';
flexible_name = true;
experiments = AP_find_experiments(animal,protocol,flexible_name);

% experiments = experiments([experiments.imaging] & [experiments.ephys]);

init_array = cell(size(experiments));

curr_day = length(experiments);
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

%% Plot psychometrics and reaction times (pre/post muscimol choiceworld, BATCH)

animals = {'AP045','AP054','AP055'};

bhv = struct('frac_go_left',cell(length(animals),1), ...
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
            
            % Performance (note that this excludes repeat on incorrect trials)
            performance = block.events.sessionPerformanceValues(:,end-10:end);
            
            % Get whether all contrasts were used
            use_all_contrasts = all(block.events.useContrastsValues(end-10:end));
            
            
            % Package final
            bhv(curr_animal).frac_go_left{curr_day}(curr_exp,:) = performance(3,:)./performance(2,:);
            bhv(curr_animal).rxn_time{curr_day}(curr_exp,:) = stim_rxn_time;
            
        end  
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
end

frac_go_left_cat = cell2mat(permute(cat(2,bhv.frac_go_left),[1,3,2]));
rxn_time_cat = cell2mat(permute(cat(2,bhv.rxn_time),[1,3,2]));

figure; 
subplot(1,2,1); hold on;
plot(performance(1,:),squeeze(frac_go_left_cat(1,:,:)),'k','linewidth',2);
plot(performance(1,:),squeeze(frac_go_left_cat(2,:,:)),'r','linewidth',2);
line(xlim,[0.5,0.5],'color','k','linestyle','--');
line([0,0],[0,1],'color','k','linestyle','--');

subplot(1,2,2); hold on;
plot(performance(1,:),squeeze(rxn_time_cat(1,:,:)),'k','linewidth',2);
plot(performance(1,:),squeeze(rxn_time_cat(2,:,:)),'r','linewidth',2);
line([-1,1],[0.5,0.5],'color','k','linestyle','--');

%% Get VFS pre/post musicmol

animals = {'AP045','AP054','AP055','AP053'};

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
            error('~= 2 experiments')
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

vfs_cat = vertcat(vfs{:});

% Plot pre/post mean
vfs_pre_mean = nanmean(cat(3,vfs_cat{:,1}),3);
vfs_post_mean = nanmean(cat(3,vfs_cat{:,2}),3);

figure; 
subplot(1,2,1);
imagesc(vfs_pre_mean);
axis image off
colormap(brewermap([],'*RdBu'));
caxis([-1,1])
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('Pre-muscimol');

subplot(1,2,2);
imagesc(vfs_post_mean);
axis image off
colormap(brewermap([],'*RdBu'));
caxis([-1,1])
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('Post-muscimol');


% Save this somewhere
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\tests\muscimol_test';
save_fn = 'muscimol_vfs.mat';
save([save_path filesep save_fn],'vfs');

%% ~~~~~~~~ After preprocessing/saving trial activity ~~~~~~~~

%% Get cortex/striatum activity pre/post muscimol

data_fns = { ...
    'trial_activity_AP_lcrGratingPassive_pre_muscimol', ...
    'trial_activity_AP_lcrGratingPassive_post_muscimol'};

mua_stim_act_exp = [];

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Split data by experiment
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    
    trial_stim_allcat_exp = mat2cell(trial_stim_allcat,trials_recording,1);
    mua_allcat_exp = mat2cell(mua_allcat, ...
        trials_recording,size(mua_allcat,2),size(mua_allcat,3));
    mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat, ...
        trials_recording,size(mua_ctxpred_allcat,2),size(mua_ctxpred_allcat,3));
    fluor_roi_deconv_exp = mat2cell(fluor_roi_deconv, ...
        trials_recording,size(fluor_roi_deconv,2),size(fluor_roi_deconv,3));
    
    % Get trials with movement during stim to exclude
    wheel_thresh = 0.025;
    quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > wheel_thresh,2);
    quiescent_trials_exp = mat2cell(quiescent_trials,trials_recording,1);
    
    % Get stimulus activity for each recording
    use_stim_time = t >= 0.05 & t <= 0.15;
    
    mua_stim_act_exp(:,:,curr_data) = cell2mat(cellfun(@(act,stim,quiescent) ...
        squeeze(nanmean(nanmean(act(stim == 3 & quiescent,use_stim_time,:,:),2),1)), ...
        mua_allcat_exp,trial_stim_allcat_exp,quiescent_trials_exp,'uni',false)');
    
    mua_ctxpred_stim_act_exp(:,:,curr_data) = cell2mat(cellfun(@(act,stim,quiescent) ...
        squeeze(nanmean(nanmean(act(stim == 3 & quiescent,use_stim_time,:,:),2),1)), ...
        mua_ctxpred_allcat_exp,trial_stim_allcat_exp,quiescent_trials_exp,'uni',false)');
    
    fluor_roi_stim_act_exp(:,:,curr_data) = cell2mat(cellfun(@(act,stim,quiescent) ...
        squeeze(nanmean(nanmean(act(stim == 3 & quiescent,use_stim_time,:,:),2),1)), ...
        fluor_roi_deconv_exp,trial_stim_allcat_exp,quiescent_trials_exp,'uni',false)');
    
    clearvars -except data_fns curr_data mua_stim_act_exp mua_ctxpred_stim_act_exp fluor_roi_stim_act_exp
    
    AP_print_progress_fraction(curr_data,length(data_fns));
    
end

% Nothing plotted at the moment












