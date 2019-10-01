%% Notes
%
% AP_ctx_str_ephys_alignment - beforehand: gets striatal domains to group MUA
%
% Batch scripts to save preprocessed data here, saved to:
% C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data


%% Cortex -> kernel-aligned striatum map (protocols separately)

disp('Cortex -> striatum regression maps across protocols');

n_aligned_depths = 4;

% Parameters for regression
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

protocols = {'vanillaChoiceworld','stimSparseNoiseUncorrAsync'};

for protocol = protocols
    protocol = cell2mat(protocol);
    
    animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
    
    init_array = cellfun(@(x) cell(0,0),animals','uni',false);
    ctx_str_kernel = init_array;
    ctx_str_expl_var = init_array;
    
    for curr_animal = 1:length(animals)
        
        animal = animals{curr_animal};
        disp(['Loading ' animal ' ' protocol]);
        
        % Only use days with choiceworld (sometimes recorded in cortex, no bhv)
        behavior_protocol = 'vanillaChoiceworld';
        behavior_experiments = AP_find_experiments(animal,behavior_protocol);
        
        curr_experiments = AP_find_experiments(animal,protocol);
        
        behavior_day = ismember({curr_experiments.day},{behavior_experiments.day});
        
        experiments = curr_experiments([curr_experiments.imaging] & [curr_experiments.ephys] & behavior_day);
        
        % Skip if this animal doesn't have this experiment
        if isempty(experiments)
            continue
        end
        
        disp(animal);
        
        load_parts.cam = false;
        load_parts.imaging = true;
        load_parts.ephys = true;
        
        for curr_day = 1:length(experiments)
            
            day = experiments(curr_day).day;
            experiment = experiments(curr_day).experiment(end);
            
            % Load data and align striatum by depth
            str_align = 'kernel';
            AP_load_experiment;
            
            %%% Load lambda from previously estimated and saved
            lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
            load(lambda_fn);
            curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
            if any(curr_animal_idx)
                curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
                if any(curr_day_idx)
                    lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
                end
            end
            
            %%% Prepare data for regression
            
            % Get time points to bin
            sample_rate = framerate*regression_params.upsample_factor;
            time_bins = frame_t(find(frame_t > ...
                regression_params.skip_seconds,1)):1/sample_rate: ...
                frame_t(find(frame_t-frame_t(end) < ...
                -regression_params.skip_seconds,1,'last'));
            time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
            
            % Deconvolve fluorescence
            fVdf_deconv = AP_deconv_wf(fVdf);
            fVdf_deconv(isnan(fVdf_deconv)) = 0;
            fVdf_deconv_resample = interp1(frame_t,fVdf_deconv(regression_params.use_svs,:)',time_bin_centers)';
            
            % Get striatum depth group by across-experiment alignment
            n_depths = n_aligned_depths;
            depth_group = aligned_str_depth_group;
            
            binned_spikes = zeros(n_depths,length(time_bin_centers));
            for curr_depth = 1:n_depths
                curr_spike_times = spike_times_timeline(depth_group == curr_depth);
                binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
            end
            
            binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
            
            %%% Regress MUA from cortex
            kernel_frames = floor(regression_params.kernel_t(1)*sample_rate): ...
                ceil(regression_params.kernel_t(2)*sample_rate);
            
            [k,ctxpred_spikes_std,explained_var] = ...
                AP_regresskernel(fVdf_deconv_resample, ...
                binned_spikes_std,kernel_frames,lambda, ...
                regression_params.zs,regression_params.cvfold, ...
                false,regression_params.use_constant);
            
            % Reshape kernel and convert to pixel space            
            aUdf = AP_align_widefield(Udf,animal,day);
            k_px = zeros(size(aUdf,1),size(aUdf,2),size(k,2),size(k,3),'single');
            for curr_spikes = 1:size(k,3)
                k_px(:,:,:,curr_spikes) = svdFrameReconstruct(aUdf(:,:,regression_params.use_svs),k(:,:,curr_spikes));
            end    
            
            % Store
            ctx_str_kernel{curr_animal}{curr_day} = k_px;
            ctx_str_expl_var{curr_animal}{curr_day} = explained_var.total;
            
            AP_print_progress_fraction(curr_day,length(experiments));
            clearvars -except regression_params n_aligned_depths ...
                animals animal curr_animal protocol ...
                experiments curr_day animal load_parts ...
                ctx_str_kernel ctx_str_expl_var
            
        end
        
        disp(['Finished ' animal]);
        
    end
    
    % Save
    save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data'];
    save([save_path filesep 'ctx_str_kernels_' protocol],'-v7.3');
    warning('saving -v7.3');
    disp(['Finished ' protocol]);
    
end


%% Choiceworld task -> striatum unit regression 

disp('Task -> Striatal unit regression')

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

unit_kernel_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'none';
        AP_load_experiment;
        
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
        
        %%% Trial-align facecam movement
        event_aligned_movement = interp1(facecam_t(~isnan(facecam_t)), ...
            frame_movement(~isnan(facecam_t)),t_peri_event);
        
        %%% Trial-align outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_bins(x,:)), ...
            find(trial_outcome == 1),'uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_bins(x,:)), ...
            find(trial_outcome == -1),'uni',false))) > 0;
        
        % Pick trials to keep
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials)' & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials');
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials');
        
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%% Regress task to cortex/striatum/cortex-predicted striatum
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.025,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Build regressors (only a subset of these are used)
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times))'.* ...
            signals_events.trialContrastValues(1:length(stimOn_times))';
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end
        
        % Stim move regressors (one for each stim when it starts to move)
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
            
        end
        
        % Stim center regressors (one for each stim when it's stopped during reward)
        unique_contrasts = unique(contrasts(contrasts > 0));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ...
                abs(stim_contrastsides) == unique_contrasts(curr_contrast));
            
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
            
        end
        
        % Move onset regressors (L/R)
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move onset x stim regressors (one for each contrast/side)
        move_onset_stim_time_absolute = arrayfun(@(curr_stim) ...
            arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & stim_contrastsides == unique_stim(curr_stim))), ...
            1:length(unique_stim),'uni',false);
        
        move_onset_stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            move_onset_stim_regressors(curr_stim,:) = ...
                histcounts(move_onset_stim_time_absolute{curr_stim},time_bins);
        end
        
        % Move ongoing regressors (L/R choice for duration of movement)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        
        move_stopped_t = 0.5;
        move_stopped_samples = round(sample_rate*move_stopped_t);
        wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.025), ...
            ones(1,move_stopped_samples),'full') > 0;
        wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);
        
        move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_L_absolute','uni',false));
        move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_R_absolute','uni',false));
        
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
        move_ongoing_regressors(2,move_ongoing_R_samples) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        % (for go cue only on late move trials)
        go_cue_regressors = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        % (for go cue with early/late move trials)
%         go_cue_regressors = zeros(1,length(time_bin_centers));
%         go_cue_regressors(1,:) = histcounts( ...
%             signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
%         go_cue_regressors(2,:) = histcounts( ...
%             signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        % (regressors for hit only)
        outcome_regressors = histcounts(reward_t_timeline,time_bins);
        % (regressors for both hit and miss)
%         outcome_regressors = zeros(2,length(time_bin_centers));
%         outcome_regressors(1,:) = histcounts( ...
%             reward_t_timeline,time_bins);
%         outcome_regressors(2,:) = histcounts( ...
%             signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate selected regressors, set parameters        
        task_regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        task_regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
                       
        task_t_shifts = { ...
            [0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [0,0.5]; ... % go cue
            [0,0.5]}; % outcome      
        
        task_regressor_sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),task_t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;        
        
        %%%%%%% Regress task -> units
        
        %%% Trial-align units
        event_aligned_unit = nan(length(stimOn_times),length(t),max(spike_templates));
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_unit = 1:max(spike_templates)
            curr_spikes = spike_times_timeline(spike_templates == curr_unit);
            event_aligned_unit(:,:,curr_unit) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Binned unit activity across the experiment
        binned_spikes = nan(max(spike_templates),length(time_bin_centers));
        for curr_unit = 1:max(spike_templates)
            curr_spike_times = spike_times_timeline(spike_templates == curr_unit);
            if isempty(curr_spike_times)
                continue
            end
            binned_spikes(curr_unit,:) = histcounts(curr_spike_times,time_bins);
        end
        
        %%% Task > striatal units
        baseline = nanmean(reshape(event_aligned_unit(:,t < 0,:),[], ...
            size(event_aligned_unit,3))*raster_sample_rate)';
        binned_spikes_baselinesubtracted = binned_spikes - baseline;
        
        [unit_taskpred_k,binned_spikes_taskpred,unit_expl_var,~] = ...
            AP_regresskernel(task_regressors,binned_spikes_baselinesubtracted, ...
            task_regressor_sample_shifts,lambda,zs,cvfold,return_constant,use_constant);
        
        % Get normalized total number of spikes
        used_spikes = spike_times_timeline > time_bins(1) & ...
            spike_times_timeline < time_bins(end);
        spike_rate = accumarray(spike_templates(used_spikes),1,[size(templates,1),1])./(time_bins(end)-time_bins(1));
        
        % Package
        unit_kernel_all(curr_animal,curr_day).template_depths = template_depths;
        unit_kernel_all(curr_animal,curr_day).spike_rate = spike_rate;
        unit_kernel_all(curr_animal,curr_day).unit_expl_var_total = unit_expl_var.total;
        unit_kernel_all(curr_animal,curr_day).unit_expl_var_partial = unit_expl_var.partial;
        unit_kernel_all(curr_animal,curr_day).unit_taskpred_k = unit_taskpred_k;
        
        % Clear
        clearvars -except ...
            task_regressor_labels task_regressor_sample_shifts ...
            regression_params animals curr_animal animal days curr_day experiments unit_kernel_all
        AP_print_progress_fraction(curr_day,length(experiments));
    end
end

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = 'unit_kernel_all.mat';
save([save_path filesep save_fn]);

disp('Finished');


%% Choiceworld trial activity (striatum depth)

clear all
disp('Choiceworld trial activity (striatum depth)')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);

mua_taskpred_k_all = cell(1,1);
mua_taskpred_all = cell(1,1);
mua_taskpred_reduced_all = cell(1,1);
mua_taskpred_expl_var_total_all = cell(1,1);
mua_taskpred_expl_var_partial_all = cell(1,1);

mua_ctxpred_taskpred_k_all = cell(1,1);
mua_ctxpred_taskpred_all = cell(1,1);
mua_ctxpred_taskpred_reduced_all = cell(1,1);
mua_ctxpred_taskpred_expl_var_total_all = cell(1,1);
mua_ctxpred_taskpred_expl_var_partial_all = cell(1,1);

fluor_taskpred_k_all = cell(1,1);
fluor_taskpred_all = cell(1,1);
fluor_taskpred_reduced_all = cell(1,1);
fluor_taskpred_expl_var_total_all = cell(1,1);
fluor_taskpred_expl_var_partial_all = cell(1,1);

wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
movement_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'depth';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = AP_align_widefield(Udf,animal,day);
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Get align times
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        
        %%% Trial-align cortex
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% Trial-align striatum
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            if isempty(curr_spikes)
                continue
            end
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regress cortex to striatum
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Deconvolve fluoresence
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        
        % Bin spikes across the experiment
        binned_spikes = nan(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            % Skip if no spikes at this depth
            if isempty(curr_spike_times)
                continue
            end
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        % Regress cortex to striatum
        [ctx_str_k,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Recast the k's into the master U
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
        
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        ctxpred_spikes = (ctxpred_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,ctxpred_spikes',t_peri_event)./raster_sample_rate;
        
        % Regress cortex to wheel velocity/speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Recast the k's into the master U
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%% Trial-align wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Trial-align facecam movement
        event_aligned_movement = interp1(facecam_t(~isnan(facecam_t)), ...
            frame_movement(~isnan(facecam_t)),t_peri_event);
        
        %%% Trial-align outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_peri_event_bins(x,:)), ...
            find(trial_outcome == 1),'uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_peri_event_bins(x,:)), ...
            find(trial_outcome == -1),'uni',false))) > 0;
        
        % Pick trials to keep
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials)' & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials)' == -1;
        R_trials = signals_events.trialSideValues(1:n_trials)' == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = 3-(abs((trial_choice(use_trials)+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%% Regress task to cortex/striatum/cortex-predicted striatum
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.025,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Build regressors (only a subset of these are used)
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times))'.* ...
            signals_events.trialContrastValues(1:length(stimOn_times))';
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end
        
        % Stim move regressors (one for each stim when it starts to move)
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
            
        end
        
        % Stim center regressors (one for each stim when it's stopped during reward)
        unique_contrasts = unique(contrasts(contrasts > 0));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ...
                abs(stim_contrastsides) == unique_contrasts(curr_contrast));
            
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
            
        end
        
        % Move onset regressors (L/R)
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move onset x stim regressors (one for each contrast/side)
        move_onset_stim_time_absolute = arrayfun(@(curr_stim) ...
            arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & stim_contrastsides == unique_stim(curr_stim))), ...
            1:length(unique_stim),'uni',false);
        
        move_onset_stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            move_onset_stim_regressors(curr_stim,:) = ...
                histcounts(move_onset_stim_time_absolute{curr_stim},time_bins);
        end
        
        % Move ongoing regressors (L/R choice for duration of movement)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        
        move_stopped_t = 0.5;
        move_stopped_samples = round(sample_rate*move_stopped_t);
        wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.025), ...
            ones(1,move_stopped_samples),'full') > 0;
        wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);
        
        move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_L_absolute','uni',false));
        move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_R_absolute','uni',false));
        
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
        move_ongoing_regressors(2,move_ongoing_R_samples) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        % (for go cue only on late move trials)
%         go_cue_regressors = histcounts( ...
%             signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        % (for go cue with early/late move trials)
        go_cue_regressors = zeros(1,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        % (regressors for hit only)
%         outcome_regressors = histcounts(reward_t_timeline,time_bins);
        % (regressors for both hit and miss)
        outcome_regressors = zeros(2,length(time_bin_centers));
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate selected regressors, set parameters
        
        task_regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        task_regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
                       
        task_t_shifts = { ...
            [0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [0,0.5]; ... % go cue
            [0,0.5]}; % outcome       
        
        % (old extended timings)
        %         task_t_shifts = {[0,0.5]; ... % stim
        %             [-0.5,1]; ... % move
        %             [-0.1,0.5]; ... % go cue
        %             [-0.5,1]}; % outcome

        % (to include stim x move)
        %         task_regressors = {stim_regressors;move_onset_regressors;move_onset_stim_regressors;go_cue_regressors;outcome_regressors};
        %         task_regressor_labels = {'Stim','Move','Stim x move','Go cue','Outcome'};
        %
        %         task_t_shifts = { ...
        %             [0,0.5]; ... % stim
        %             [-0.5,1]; ... % move
        %             [-0.5,1]; ... % stim x move
        %             [0,0.5]; ... % go cue
        %             [0,0.5]}; % outcome
        
        task_regressor_sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),task_t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;
        
        % Regression task -> MUA
        baseline = nanmean(reshape(event_aligned_mua(:,t < 0,:),[], ...
            size(event_aligned_mua,3))*raster_sample_rate,1)';
        activity = single(binned_spikes) - baseline;
        
        [mua_taskpred_k,mua_taskpred_long,mua_taskpred_expl_var,mua_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        mua_taskpred = ...
            interp1(time_bin_centers,mua_taskpred_long',t_peri_event)./raster_sample_rate;
        
        mua_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,mua_taskpred_reduced_long(:,:,x)', ...
            t_peri_event)./raster_sample_rate,permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        % Regression task -> MUA-ctxpred
        baseline = nanmean(reshape(event_aligned_mua_ctxpred(:,t < 0,:),[], ...
            size(event_aligned_mua_ctxpred,3))*raster_sample_rate,1)';
        activity = single(ctxpred_spikes) - baseline;
        
        [mua_ctxpred_taskpred_k,mua_ctxpred_taskpred_long,mua_ctxpred_taskpred_expl_var,mua_ctxpred_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        mua_ctxpred_taskpred = ...
            interp1(time_bin_centers,mua_ctxpred_taskpred_long',t_peri_event)./raster_sample_rate;
        
        mua_ctxpred_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,mua_ctxpred_taskpred_reduced_long(:,:,x)', ...
            t_peri_event)./raster_sample_rate,permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        % Regression task -> (master U, deconvolved) fluor
        event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
        fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);
        
        baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,:),[],size(event_aligned_V_deconv,3)))';
        activity = single(fVdf_deconv_resample_recast(use_components,:))-baseline;
        
        [fluor_taskpred_k,fluor_taskpred_long,fluor_taskpred_expl_var,fluor_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        fluor_taskpred = ...
            interp1(time_bin_centers,fluor_taskpred_long',t_peri_event);
        
        fluor_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,fluor_taskpred_reduced_long(:,:,x)', ...
            t_peri_event),permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        mua_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_taskpred_k;
        mua_taskpred_all{curr_animal,1}{curr_day,1} = mua_taskpred(use_trials,:,:,:);
        mua_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_taskpred_reduced(use_trials,:,:,:);
        mua_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = mua_taskpred_expl_var.total;
        mua_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = mua_taskpred_expl_var.partial;
        
        mua_ctxpred_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_k;
        mua_ctxpred_taskpred_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred(use_trials,:,:,:);
        mua_ctxpred_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_reduced(use_trials,:,:,:);
        mua_ctxpred_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_expl_var.total;
        mua_ctxpred_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_expl_var.partial;
        
        fluor_taskpred_k_all{curr_animal,1}{curr_day,1} = fluor_taskpred_k;
        fluor_taskpred_all{curr_animal,1}{curr_day,1} = fluor_taskpred(use_trials,:,:,:);
        fluor_taskpred_reduced_all{curr_animal,1}{curr_day,1} = fluor_taskpred_reduced(use_trials,:,:,:);
        fluor_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.total;
        fluor_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.partial;
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        movement_all{curr_animal,1}{curr_day,1} = event_aligned_movement(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars -except curr_animal animal protocol experiments load_parts curr_day ...
            n_aligned_depths regression_params animals t ...
            ...
            fluor_all ...
            mua_all ...
            mua_ctxpred_all ...
            ...
            mua_taskpred_k_all ...
            mua_taskpred_all ...
            mua_taskpred_reduced_all ...
            mua_taskpred_expl_var_total_all ...
            mua_taskpred_expl_var_partial_all ...
            ...
            mua_ctxpred_taskpred_k_all ...
            mua_ctxpred_taskpred_all ...
            mua_ctxpred_taskpred_reduced_all ...
            mua_ctxpred_taskpred_expl_var_total_all ...
            mua_ctxpred_taskpred_expl_var_partial_all ...
            ...
            fluor_taskpred_k_all ...
            fluor_taskpred_all ...
            fluor_taskpred_reduced_all ...
            fluor_taskpred_expl_var_total_all ...
            fluor_taskpred_expl_var_partial_all ...
            ...
            wheel_ctxpred_all ...
            ctx_str_k_all ...
            ctx_wheel_k_all ...
            wheel_all ...
            movement_all ...
            outcome_all ...
            D_all ...
            ...
            task_regressor_labels ...
            task_regressor_sample_shifts
        
    end
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_taskpred_expl_var_total_all ...
    mua_taskpred_expl_var_partial_all ...
    ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    mua_ctxpred_taskpred_expl_var_total_all ...
    mua_ctxpred_taskpred_expl_var_partial_all ...
    ...
    fluor_taskpred_k_all ...
    fluor_taskpred_all ...
    fluor_taskpred_reduced_all ...
    fluor_taskpred_expl_var_total_all ...
    fluor_taskpred_expl_var_partial_all ...
    ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    movement_all ...
    outcome_all ...
    D_all ...
    ...
    task_regressor_labels ...
    task_regressor_sample_shifts

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_' num2str(n_aligned_depths) 'strdepth'];
save([save_path filesep save_fn],'-v7.3');


%% Choiceworld trial activity (striatum domain)

clear all
disp('Choiceworld trial activity (striatum domain)')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);

mua_taskpred_k_all = cell(1,1);
mua_taskpred_all = cell(1,1);
mua_taskpred_reduced_all = cell(1,1);
mua_taskpred_expl_var_total_all = cell(1,1);
mua_taskpred_expl_var_partial_all = cell(1,1);

mua_ctxpred_taskpred_k_all = cell(1,1);
mua_ctxpred_taskpred_all = cell(1,1);
mua_ctxpred_taskpred_reduced_all = cell(1,1);
mua_ctxpred_taskpred_expl_var_total_all = cell(1,1);
mua_ctxpred_taskpred_expl_var_partial_all = cell(1,1);

fluor_taskpred_k_all = cell(1,1);
fluor_taskpred_all = cell(1,1);
fluor_taskpred_reduced_all = cell(1,1);
fluor_taskpred_expl_var_total_all = cell(1,1);
fluor_taskpred_expl_var_partial_all = cell(1,1);

wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
movement_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = AP_align_widefield(Udf,animal,day);
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Get align times
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        
        %%% Trial-align cortex
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% Trial-align striatum
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            if isempty(curr_spikes)
                continue
            end
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regress cortex to striatum
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Deconvolve fluoresence
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        
        % Bin spikes across the experiment
        binned_spikes = nan(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            % Skip if no spikes at this depth
            if isempty(curr_spike_times)
                continue
            end
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        % Regress cortex to striatum
        [ctx_str_k,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Recast the k's into the master U
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
        
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        ctxpred_spikes = (ctxpred_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,ctxpred_spikes',t_peri_event)./raster_sample_rate;
        
        % Regress cortex to wheel velocity/speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Recast the k's into the master U
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%% Trial-align wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Trial-align facecam movement
        event_aligned_movement = interp1(facecam_t(~isnan(facecam_t)), ...
            frame_movement(~isnan(facecam_t)),t_peri_event);
        
        %%% Trial-align outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_peri_event_bins(x,:)), ...
            find(trial_outcome == 1),'uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_peri_event_bins(x,:)), ...
            find(trial_outcome == -1),'uni',false))) > 0;
        
        % Pick trials to keep
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials)' & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials)' == -1;
        R_trials = signals_events.trialSideValues(1:n_trials)' == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = 3-(abs((trial_choice(use_trials)+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%% Regress task to cortex/striatum/cortex-predicted striatum
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.025,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Build regressors (only a subset of these are used)
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times))'.* ...
            signals_events.trialContrastValues(1:length(stimOn_times))';
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end
        
        % Stim move regressors (one for each stim when it starts to move)
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
            
        end
        
        % Stim center regressors (one for each stim when it's stopped during reward)
        unique_contrasts = unique(contrasts(contrasts > 0));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ...
                abs(stim_contrastsides) == unique_contrasts(curr_contrast));
            
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
            
        end
        
        % Move onset regressors (L/R)
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move onset x stim regressors (one for each contrast/side)
        move_onset_stim_time_absolute = arrayfun(@(curr_stim) ...
            arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & stim_contrastsides == unique_stim(curr_stim))), ...
            1:length(unique_stim),'uni',false);
        
        move_onset_stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            move_onset_stim_regressors(curr_stim,:) = ...
                histcounts(move_onset_stim_time_absolute{curr_stim},time_bins);
        end
        
        % Move ongoing regressors (L/R choice for duration of movement)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        
        move_stopped_t = 0.5;
        move_stopped_samples = round(sample_rate*move_stopped_t);
        wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.025), ...
            ones(1,move_stopped_samples),'full') > 0;
        wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);
        
        move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_L_absolute','uni',false));
        move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_R_absolute','uni',false));
        
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
        move_ongoing_regressors(2,move_ongoing_R_samples) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        % (for go cue only on late move trials)
%         go_cue_regressors = histcounts( ...
%             signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        % (for go cue with early/late move trials)
        go_cue_regressors = zeros(1,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        % (regressors for hit only)
%         outcome_regressors = histcounts(reward_t_timeline,time_bins);
        % (regressors for both hit and miss)
        outcome_regressors = zeros(2,length(time_bin_centers));
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate selected regressors, set parameters
        
        task_regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        task_regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
                       
        task_t_shifts = { ...
            [0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [0,0.5]; ... % go cue
            [0,0.5]}; % outcome       
        
        % (old extended timings)
        %         task_t_shifts = {[0,0.5]; ... % stim
        %             [-0.5,1]; ... % move
        %             [-0.1,0.5]; ... % go cue
        %             [-0.5,1]}; % outcome

        % (to include stim x move)
        %         task_regressors = {stim_regressors;move_onset_regressors;move_onset_stim_regressors;go_cue_regressors;outcome_regressors};
        %         task_regressor_labels = {'Stim','Move','Stim x move','Go cue','Outcome'};
        %
        %         task_t_shifts = { ...
        %             [0,0.5]; ... % stim
        %             [-0.5,1]; ... % move
        %             [-0.5,1]; ... % stim x move
        %             [0,0.5]; ... % go cue
        %             [0,0.5]}; % outcome
        
        task_regressor_sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),task_t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;
        
        % Regression task -> MUA
        baseline = nanmean(reshape(event_aligned_mua(:,t < 0,:),[], ...
            size(event_aligned_mua,3))*raster_sample_rate,1)';
        activity = single(binned_spikes) - baseline;
        
        [mua_taskpred_k,mua_taskpred_long,mua_taskpred_expl_var,mua_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        mua_taskpred = ...
            interp1(time_bin_centers,mua_taskpred_long',t_peri_event)./raster_sample_rate;
        
        mua_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,mua_taskpred_reduced_long(:,:,x)', ...
            t_peri_event)./raster_sample_rate,permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        % Regression task -> MUA-ctxpred
        baseline = nanmean(reshape(event_aligned_mua_ctxpred(:,t < 0,:),[], ...
            size(event_aligned_mua_ctxpred,3))*raster_sample_rate,1)';
        activity = single(ctxpred_spikes) - baseline;
        
        [mua_ctxpred_taskpred_k,mua_ctxpred_taskpred_long,mua_ctxpred_taskpred_expl_var,mua_ctxpred_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        mua_ctxpred_taskpred = ...
            interp1(time_bin_centers,mua_ctxpred_taskpred_long',t_peri_event)./raster_sample_rate;
        
        mua_ctxpred_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,mua_ctxpred_taskpred_reduced_long(:,:,x)', ...
            t_peri_event)./raster_sample_rate,permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        % Regression task -> (master U, deconvolved) fluor
        event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
        fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);
        
        baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,:),[],size(event_aligned_V_deconv,3)))';
        activity = single(fVdf_deconv_resample_recast(use_components,:))-baseline;
        
        [fluor_taskpred_k,fluor_taskpred_long,fluor_taskpred_expl_var,fluor_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        fluor_taskpred = ...
            interp1(time_bin_centers,fluor_taskpred_long',t_peri_event);
        
        fluor_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,fluor_taskpred_reduced_long(:,:,x)', ...
            t_peri_event),permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        mua_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_taskpred_k;
        mua_taskpred_all{curr_animal,1}{curr_day,1} = mua_taskpred(use_trials,:,:,:);
        mua_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_taskpred_reduced(use_trials,:,:,:);
        mua_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = mua_taskpred_expl_var.total;
        mua_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = mua_taskpred_expl_var.partial;
        
        mua_ctxpred_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_k;
        mua_ctxpred_taskpred_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred(use_trials,:,:,:);
        mua_ctxpred_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_reduced(use_trials,:,:,:);
        mua_ctxpred_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_expl_var.total;
        mua_ctxpred_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_expl_var.partial;
        
        fluor_taskpred_k_all{curr_animal,1}{curr_day,1} = fluor_taskpred_k;
        fluor_taskpred_all{curr_animal,1}{curr_day,1} = fluor_taskpred(use_trials,:,:,:);
        fluor_taskpred_reduced_all{curr_animal,1}{curr_day,1} = fluor_taskpred_reduced(use_trials,:,:,:);
        fluor_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.total;
        fluor_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.partial;
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        movement_all{curr_animal,1}{curr_day,1} = event_aligned_movement(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars -except curr_animal animal protocol experiments load_parts curr_day ...
            n_aligned_depths regression_params animals t ...
            ...
            fluor_all ...
            mua_all ...
            mua_ctxpred_all ...
            ...
            mua_taskpred_k_all ...
            mua_taskpred_all ...
            mua_taskpred_reduced_all ...
            mua_taskpred_expl_var_total_all ...
            mua_taskpred_expl_var_partial_all ...
            ...
            mua_ctxpred_taskpred_k_all ...
            mua_ctxpred_taskpred_all ...
            mua_ctxpred_taskpred_reduced_all ...
            mua_ctxpred_taskpred_expl_var_total_all ...
            mua_ctxpred_taskpred_expl_var_partial_all ...
            ...
            fluor_taskpred_k_all ...
            fluor_taskpred_all ...
            fluor_taskpred_reduced_all ...
            fluor_taskpred_expl_var_total_all ...
            fluor_taskpred_expl_var_partial_all ...
            ...
            wheel_ctxpred_all ...
            ctx_str_k_all ...
            ctx_wheel_k_all ...
            wheel_all ...
            movement_all ...
            outcome_all ...
            D_all ...
            ...
            task_regressor_labels ...
            task_regressor_sample_shifts
        
    end
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_taskpred_expl_var_total_all ...
    mua_taskpred_expl_var_partial_all ...
    ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    mua_ctxpred_taskpred_expl_var_total_all ...
    mua_ctxpred_taskpred_expl_var_partial_all ...
    ...
    fluor_taskpred_k_all ...
    fluor_taskpred_all ...
    fluor_taskpred_reduced_all ...
    fluor_taskpred_expl_var_total_all ...
    fluor_taskpred_expl_var_partial_all ...
    ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    movement_all ...
    outcome_all ...
    D_all ...
    ...
    task_regressor_labels ...
    task_regressor_sample_shifts

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld'];
save([save_path filesep save_fn],'-v7.3');


%% Choiceworld trial activity (wf-only days)

clear all
disp('Choiceworld trial activity (wf-only days)')

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

% Initialize all saved variables for indexing
fluor_all = cell(1,1);

fluor_taskpred_k_all = cell(1,1);
fluor_taskpred_all = cell(1,1);
fluor_taskpred_reduced_all = cell(1,1);
fluor_taskpred_expl_var_total_all = cell(1,1);
fluor_taskpred_expl_var_partial_all = cell(1,1);

wheel_all = cell(1,1);
movement_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & ~[experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = AP_align_widefield(Udf,animal,day);
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Get align times
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        
        %%% Trial-align cortex
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Deconvolve fluoresence
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        
        %%% Trial-align wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Trial-align outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_peri_event_bins(x,:)), ...
            find(trial_outcome == 1),'uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_peri_event_bins(x,:)), ...
            find(trial_outcome == -1),'uni',false))) > 0;
        
        % Pick trials to keep
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials)' & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials)' == -1;
        R_trials = signals_events.trialSideValues(1:n_trials)' == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = 3-(abs((trial_choice(use_trials)+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%% Regress task to cortex/striatum/cortex-predicted striatum
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.025,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Build regressors (only a subset of these are used)
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times))'.* ...
            signals_events.trialContrastValues(1:length(stimOn_times))';
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end
        
        % Stim move regressors (one for each stim when it starts to move)
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
            
        end
        
        % Stim center regressors (one for each stim when it's stopped during reward)
        unique_contrasts = unique(contrasts(contrasts > 0));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ...
                abs(stim_contrastsides) == unique_contrasts(curr_contrast));
            
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
            
        end
        
        % Move onset regressors (L/R)
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move onset x stim regressors (one for each contrast/side)
        move_onset_stim_time_absolute = arrayfun(@(curr_stim) ...
            arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & stim_contrastsides == unique_stim(curr_stim))), ...
            1:length(unique_stim),'uni',false);
        
        move_onset_stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            move_onset_stim_regressors(curr_stim,:) = ...
                histcounts(move_onset_stim_time_absolute{curr_stim},time_bins);
        end
        
        % Move ongoing regressors (L/R choice for duration of movement)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        
        move_stopped_t = 0.5;
        move_stopped_samples = round(sample_rate*move_stopped_t);
        wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.025), ...
            ones(1,move_stopped_samples),'full') > 0;
        wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);
        
        move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_L_absolute','uni',false));
        move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_R_absolute','uni',false));
        
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
        move_ongoing_regressors(2,move_ongoing_R_samples) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        % (for go cue only on late move trials)
%         go_cue_regressors = histcounts( ...
%             signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        % (for go cue with early/late move trials)
        go_cue_regressors = zeros(1,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        % (regressors for hit only)
%         outcome_regressors = histcounts(reward_t_timeline,time_bins);
        % (regressors for both hit and miss)
        outcome_regressors = zeros(2,length(time_bin_centers));
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate selected regressors, set parameters
        task_regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        task_regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
        
        task_t_shifts = { ...
            [0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [0,0.5]; ... % go cue
            [0,0.5]}; % outcome
        
        task_regressor_sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),task_t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;
        
        % Regression task -> (master U, deconvolved) fluor
        event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
        fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);
        
        baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,:),[],size(event_aligned_V_deconv,3)))';
        activity = single(fVdf_deconv_resample_recast(use_components,:))-baseline;
        
        [fluor_taskpred_k,fluor_taskpred_long,fluor_taskpred_expl_var,fluor_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        fluor_taskpred = ...
            interp1(time_bin_centers,fluor_taskpred_long',t_peri_event);
        
        fluor_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,fluor_taskpred_reduced_long(:,:,x)', ...
            t_peri_event),permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        
        fluor_taskpred_k_all{curr_animal,1}{curr_day,1} = fluor_taskpred_k;
        fluor_taskpred_all{curr_animal,1}{curr_day,1} = fluor_taskpred(use_trials,:,:,:);
        fluor_taskpred_reduced_all{curr_animal,1}{curr_day,1} = fluor_taskpred_reduced(use_trials,:,:,:);
        fluor_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.total;
        fluor_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.partial;
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars -except curr_animal animal protocol experiments load_parts curr_day ...
            n_aligned_depths regression_params animals t ...
            ...
            fluor_all ...
            ...
            mua_taskpred_k_all ...
            mua_taskpred_all ...
            mua_taskpred_reduced_all ...
            mua_taskpred_expl_var_total_all ...
            mua_taskpred_expl_var_partial_all ...
            ...
            mua_ctxpred_taskpred_k_all ...
            mua_ctxpred_taskpred_all ...
            mua_ctxpred_taskpred_reduced_all ...
            mua_ctxpred_taskpred_expl_var_total_all ...
            mua_ctxpred_taskpred_expl_var_partial_all ...
            ...
            fluor_taskpred_k_all ...
            fluor_taskpred_all ...
            fluor_taskpred_reduced_all ...
            fluor_taskpred_expl_var_total_all ...
            fluor_taskpred_expl_var_partial_all ...
            ...
            wheel_ctxpred_all ...
            ctx_str_k_all ...
            ctx_wheel_k_all ...
            wheel_all ...
            movement_all ...
            outcome_all ...
            D_all ...
            ...
            task_regressor_labels ...
            task_regressor_sample_shifts
        
    end
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_taskpred_expl_var_total_all ...
    mua_taskpred_expl_var_partial_all ...
    ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    mua_ctxpred_taskpred_expl_var_total_all ...
    mua_ctxpred_taskpred_expl_var_partial_all ...
    ...
    fluor_taskpred_k_all ...
    fluor_taskpred_all ...
    fluor_taskpred_reduced_all ...
    fluor_taskpred_expl_var_total_all ...
    fluor_taskpred_expl_var_partial_all ...
    ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    movement_all ...
    outcome_all ...
    D_all ...
    ...
    task_regressor_labels ...
    task_regressor_sample_shifts

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_wfonly'];
save([save_path filesep save_fn],'-v7.3');


%% Passive trial activity (all protocols and animal groups)

clear all

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

trained_animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
naive_animals = {'AP032','AP033','AP034','AP035','AP036'};
protocols = {'AP_choiceWorldStimPassive','stimKalatsky'};

for animal_group = {'trained','naive'}
    animal_group = cell2mat(animal_group);
    
    switch animal_group
        case 'trained'
            animals = trained_animals;
        case 'naive'
            animals = naive_animals;
    end
    
    for protocol = protocols
        protocol = cell2mat(protocol);        
        
        disp(['Passive trial activity: ' animal_group ' ' protocol])
        
        % Initialize all saved variables for indexing and common structure
        init_array = cellfun(@(x) cell(0,0),animals','uni',false);
        
        fluor_all = init_array;
        mua_all = init_array;
        mua_ctxpred_all = init_array;
        ctx_str_k_all = init_array;
        wheel_all = init_array;
        D_all = init_array;
        
        for curr_animal = 1:length(animals)
            
            animal = animals{curr_animal};
            
            if strcmp(animal_group,'trained')
                % Only use days with choiceworld (sometimes recorded in cortex, no bhv)
                bhv_protocol = 'vanillaChoiceworld';
                behavior_experiments = AP_find_experiments(animal,bhv_protocol);
                passive_experiments = AP_find_experiments(animal,protocol);
                behavior_day = ismember({passive_experiments.day},{behavior_experiments.day});
                experiments = passive_experiments([passive_experiments.imaging] & [passive_experiments.ephys] & behavior_day);
            elseif strcmp(animal_group,'naive')
                experiments = AP_find_experiments(animal,protocol);
                experiments = experiments([experiments.imaging] & [experiments.ephys]);
            end
            
            load_parts.cam = false;
            load_parts.imaging = true;
            load_parts.ephys = true;
            
            disp(['Loading ' animal]);
            
            for curr_day = 1:length(experiments)
                
                day = experiments(curr_day).day;
                experiment = experiments(curr_day).experiment(end);
                
                % Load experiment
                str_align = 'kernel';
                AP_load_experiment;
                
                % Prepare fluorescence
                % Convert U to master U
                load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
                Udf_aligned = AP_align_widefield(Udf,animal,day);
                fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
                
                % Set components to keep
                use_components = 1:200;
                
                % (aligned striatum depths)
                n_depths = n_aligned_depths;
                depth_group = aligned_str_depth_group;
                
                % Get event-aligned activity
                raster_window = [-0.5,2];
                upsample_factor = 1;
                raster_sample_rate = 1/(framerate*upsample_factor);
                t = raster_window(1):raster_sample_rate:raster_window(2);
                
                % Align (only to stim)
                use_align = stimOn_times;
                use_align(isnan(use_align)) = 0;
                
                t_peri_event = bsxfun(@plus,use_align,t);
                
                %%% Fluorescence
                event_aligned_V = ...
                    interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
                
                %%% MUA
                event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
                t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
                for curr_depth = 1:n_depths
                    curr_spikes = spike_times_timeline(depth_group == curr_depth);
                    % (for only msns in depth group)
                    %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
                    %                     ismember(spike_templates,find(msn)));
                    
                    event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                        histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                        [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
                end
                
                %%% Regressed fluorescence -> MUA
                
                % Get time points to bin
                sample_rate = framerate*regression_params.upsample_factor;
                time_bins = frame_t(find(frame_t > ...
                    regression_params.skip_seconds,1)):1/sample_rate: ...
                    frame_t(find(frame_t-frame_t(end) < ...
                    -regression_params.skip_seconds,1,'last'));
                time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
                
                % Get upsampled dVdf's
                dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
                    diff(fVdf,[],2)',time_bin_centers)';
                
                % Deconvolve fluorescence
                fVdf_deconv = AP_deconv_wf(fVdf);
                fVdf_deconv(isnan(fVdf_deconv)) = 0;
                fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
                
                binned_spikes = zeros(n_depths,length(time_bin_centers));
                for curr_depth = 1:n_depths
                    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
                    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
                end
                binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
                binned_spikes_std(isnan(binned_spikes_std)) = 0;
                
                % Load lambda from previously estimated and saved
                lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
                load(lambda_fn);
                curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
                if any(curr_animal_idx)
                    curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
                    if any(curr_day_idx)
                        lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
                    end
                end
                
                kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
                    round(regression_params.kernel_t(2)*sample_rate);
                
                [ctx_str_k,ctxpred_spikes_std,explained_var] = ...
                    AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
                    binned_spikes_std,kernel_frames,lambda, ...
                    regression_params.zs,regression_params.cvfold, ...
                    true,regression_params.use_constant);
                
                % Cast the k's into the master to save small V's
                ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
                    reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
                    size(ctx_str_k{1}));
                
                % Re-scale the prediction (subtract offset, multiply, add scaled offset)
                ctxpred_spikes = (ctxpred_spikes_std - squeeze(ctx_str_k{end})).* ...
                    nanstd(binned_spikes,[],2) + ...
                    nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
                
                event_aligned_mua_ctxpred = ...
                    interp1(time_bin_centers,ctxpred_spikes',t_peri_event)./raster_sample_rate;
                
                %%% Wheel velocity
                event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
                    wheel_velocity,t_peri_event);
                
                % Pick trials to use
                use_trials = true(size(stimIDs));
                
                % Get stim info
                D = struct;
                D.stimulus = stimIDs;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Store everything
                fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
                mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
                
                ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
                mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
                
                wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
                
                D_all{curr_animal,1}{curr_day,1} = D;
                
                AP_print_progress_fraction(curr_day,length(experiments));
            end
            
            clearvars -except ...
                trained_animals naive_animals protocols ...
                protocol animal_group curr_animal ...
                n_aligned_depths regression_params animals t ...
                fluor_all ...
                mua_all ...
                mua_ctxpred_all ...
                ctx_str_k_all ...
                wheel_all ...
                D_all
            
        end
        
        clearvars -except ...
            trained_animals naive_animals protocols ...
            protocol animal_group...
            n_aligned_depths regression_params animals t ...
            fluor_all ...
            mua_all ...
            mua_ctxpred_all ...
            ctx_str_k_all ...
            wheel_all ...
            D_all
        
        disp('Finished loading all')
        
        save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
        save_fn = ['trial_activity_' protocol '_' animal_group];
        save([save_path filesep save_fn],'-v7.3');
        
    end
end



%% ~~~~~~~~~~~~~ Widefield ROIs ~~~~~~~~~~~~~

%% (lock on)
if false

%% Make reference image for drawing widefield ROIs

data_fn = 'trial_activity_choiceworld';
exclude_data = true;
AP_load_concat_normalize_ctx_str;

%%% Generate reference image for widefield ROIs
rxn_time_use = [0.1,0.3];

use_trials = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == -1 & ...
    move_t >= rxn_time_use(1) & ...
    move_t <= rxn_time_use(2);

plot_px = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    diff(squeeze(nanmean(fluor_allcat(use_trials,:,:),1))',[],2));

use_t = t >= 0 & t <= 0.8;
plot_px_use_t = double(plot_px(:,:,use_t));
plot_px_use_t(plot_px_use_t < 0) = 0;
plot_px_use_t = bsxfun(@rdivide,plot_px_use_t,max(max(plot_px_use_t,[],1),[],2));

plot_px_com = sum((plot_px_use_t.*permute(1:sum(use_t),[1,3,2])),3)./sum(plot_px_use_t,3);

t_leeway = 20;
plot_px_colored = ...
    ind2rgb(round(mat2gray(plot_px_com,[t_leeway,sum(use_t)-t_leeway])*255),jet(255));

plot_px_alpha = mat2gray(max(plot_px_use_t,[],3), ...
    [0,prctile(reshape(max(plot_px_use_t,[],3),[],1),99)]);

f = figure;

h = image(plot_px_colored);
set(h,'AlphaData',plot_px_alpha);

axis image off
AP_reference_outline('ccf_aligned','k');
title('Widefield ROI reference image')

wf_roi_ref = plot_px_colored;
wf_roi_ref_alpha = plot_px_alpha;

wf_roi_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois';
wf_ref_fn = 'wf_roi_ref';
savefig(f,[wf_roi_path filesep wf_ref_fn]);


%% Draw widefield ROIs

% Set ROIs to draw
roi_areas = {'V1p','V1c','AM','RSPa','RSPp','PPC','FRm','FRa','SMl','SMf'};

% Load reference image
wf_roi_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois';
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

wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
save(wf_roi_fn,'wf_roi');
disp('Saved new widefield ROIs');

%% Create ROIs from kernel templates

% Load kernel templates
n_aligned_depths = 4;
kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template_' num2str(n_aligned_depths) '_depths.mat'];
load(kernel_template_fn);

kernel_bw = false(size(kernel_template));

% Set cutoff by std, get rid of small islands, keep only ipsi side
kernel_std_prctile = 1;
smallest_area = 2000;

for curr_kernel = 1:size(kernel_bw,3)
    curr_kernel_cutoff = std(abs(reshape(kernel_template(:,:,curr_kernel),[],1)));
    curr_kernel_bw = kernel_template(:,:,curr_kernel) > kernel_std_prctile*curr_kernel_cutoff;
    kernel_bw(:,:,curr_kernel) = bwareaopen(curr_kernel_bw,smallest_area);
end

bregma = allenCCFbregma;
ccf_tform_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\ccf_tform'];
load(ccf_tform_fn);

um2pixel = 20.6;
bregma_resize = bregma*(10/um2pixel);
bregma_align = [bregma_resize([3,1]),1]*ccf_tform.T;
kernel_bw(:,round(bregma_align(1)):end,:) = false;

% Make a max-weighted version of each kernel
kernel_max_weighted = kernel_template./ ...
    max(abs(reshape(kernel_template,1,[],size(kernel_template,3))),[],2);

% Plot the template kernels and ROIs
figure;
for i = 1:n_aligned_depths
    subplot(3,n_aligned_depths,i);
    imagesc(kernel_template(:,:,i));
    AP_reference_outline('ccf_aligned','k');
    axis image off;
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(gca,brewermap([],'*RdBu'));
    title('Kernel template')
    
    subplot(3,n_aligned_depths,i+n_aligned_depths);
    imagesc(kernel_bw(:,:,i));
    AP_reference_outline('ccf_aligned','r');
    axis image off;
    colormap(gca,gray);
    title('ROI BW')
    
    subplot(3,n_aligned_depths,i+n_aligned_depths*2);
    imagesc(kernel_max_weighted(:,:,i));
    AP_reference_outline('ccf_aligned','k');
    axis image off;
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(gca,brewermap([],'*RdBu'));
    title('ROI weighted')
end

% Get the average full kernel
n_aligned_depths = 4;
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_ephys';
k_fn = [data_path filesep 'wf_ephys_maps_vanillaChoiceworld_' num2str(n_aligned_depths) '_depths'];
load(k_fn);
k_px_trained_cat = cellfun(@(x) x(:,:,end:-1:1,:),[batch_vars(1:6).r_px],'uni',false);
k_px_trained = nanmean(double(cat(5,k_px_trained_cat{:})),5);

% Save kernel ROIs
kernel_roi = struct;
kernel_roi.bw = kernel_bw;
kernel_roi.max_weighted = kernel_max_weighted;
kernel_roi.kernel = k_px_trained;
kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
save(kernel_roi_fn,'kernel_roi');
disp('Saved kernel ROIs');

%% (lock off)
end

%% ~~~~~~~~~~~ UNUSED ~~~~~~~~~~~~~~~

%% (lock on)
if false

%% Passive choiceworld trial activity (trained)

clear all
disp('Passive choiceworld trial activity (trained)')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);
wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    % Only use days with choiceworld (sometimes recorded in cortex, no bhv)
    protocol = 'vanillaChoiceworld';
    behavior_experiments = AP_find_experiments(animal,protocol);
    
    protocol = 'AP_choiceWorldStimPassive';
    passive_experiments = AP_find_experiments(animal,protocol);
    
    behavior_day = ismember({passive_experiments.day},{behavior_experiments.day});
    
    experiments = passive_experiments([passive_experiments.imaging] & [passive_experiments.ephys] & behavior_day);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Align (only to stim)
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf,[],2)',time_bin_centers)';
        
        %%%% TESTING: DECONV
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        %%%%
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        %         [~,predicted_spikes_std,explained_var] = ...
        %             AP_regresskernel(dfVdf_resample, ...
        %             binned_spikes_std,kernel_frames,lambda, ...
        %             regression_params.zs,regression_params.cvfold, ...
        %             false,regression_params.use_constant);
        
        %%%% TESTING DECONV
        [ctx_str_k,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Shove the k's into the master to save small V's
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
        
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        ctxpred_spikes = (ctxpred_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,ctxpred_spikes',t_peri_event)./raster_sample_rate;
        
        %%%%
        
        %%%% TESTING CTX->WHEEL VELOCITY
        % Resample wheel, predict both velocity and speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Shove the k's into the master to save small V's
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%%%
        
        %%% Wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        % Pick trials to use
        use_trials = true(size(stimIDs));
        
        % Get stim info
        D = struct;
        D.stimulus = stimIDs;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);
        
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except curr_animal ...
        n_aligned_depths regression_params animals t ...
        fluor_all ...
        mua_all ...
        mua_ctxpred_all ...
        wheel_ctxpred_all ...
        ctx_str_k_all ...
        ctx_wheel_k_all ...
        wheel_all ...
        D_all
    
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    D_all

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_passive_choiceworld_trained_DECONVTEST'];
save([save_path filesep save_fn],'-v7.3');

%% Passive choiceworld trial activity (naive)

clear all
disp('Passive choiceworld trial activity (naive)')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP032','AP033','AP034','AP035','AP036'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);
wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    protocol = 'AP_choiceWorldStimPassive';
    experiments = AP_find_experiments(animal,protocol);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Align (only to stim)
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf,[],2)',time_bin_centers)';
        
        %%%% TESTING: DECONV
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        %%%%
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        %         [~,predicted_spikes_std,explained_var] = ...
        %             AP_regresskernel(dfVdf_resample, ...
        %             binned_spikes_std,kernel_frames,lambda, ...
        %             regression_params.zs,regression_params.cvfold, ...
        %             false,regression_params.use_constant);
        
        %%%% TESTING DECONV
        [ctx_str_k,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Shove the k's into the master to save small V's
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
        
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        ctxpred_spikes = (ctxpred_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,ctxpred_spikes',t_peri_event)./raster_sample_rate;
        
        %%%%
        
        %%%% TESTING CTX->WHEEL VELOCITY
        % Resample wheel, predict both velocity and speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Shove the k's into the master to save small V's
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%%%
        
        %%% Wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        % Pick trials to use
        use_trials = true(size(stimIDs));
        
        % Get stim info
        D = struct;
        D.stimulus = stimIDs;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);
        
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except curr_animal ...
        n_aligned_depths regression_params animals t ...
        fluor_all ...
        mua_all ...
        mua_ctxpred_all ...
        wheel_ctxpred_all ...
        ctx_str_k_all ...
        ctx_wheel_k_all ...
        wheel_all ...
        D_all
    
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    D_all

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_passive_choiceworld_naive_DECONVTEST'];
save([save_path filesep save_fn],'-v7.3');


%% Passive fullscreen trial activity (trained)

clear all
disp('Passive fullscreen trial activity (trained)')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);
wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    % Only use days with choiceworld (sometimes recorded in cortex, no bhv)
    protocol = 'vanillaChoiceworld';
    behavior_experiments = AP_find_experiments(animal,protocol);
    
    protocol = 'stimKalatsky';
    passive_experiments = AP_find_experiments(animal,protocol);
    
    behavior_day = ismember({passive_experiments.day},{behavior_experiments.day});
    
    experiments = passive_experiments([passive_experiments.imaging] & [passive_experiments.ephys] & behavior_day);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Align (only to stim)
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf,[],2)',time_bin_centers)';
        
        %%%% TESTING: DECONV
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        %%%%
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        %         [~,predicted_spikes_std,explained_var] = ...
        %             AP_regresskernel(dfVdf_resample, ...
        %             binned_spikes_std,kernel_frames,lambda, ...
        %             regression_params.zs,regression_params.cvfold, ...
        %             false,regression_params.use_constant);
        
        %%%% TESTING DECONV
        [ctx_str_k,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Shove the k's into the master to save small V's
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
        
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        ctxpred_spikes = (ctxpred_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,ctxpred_spikes',t_peri_event)./raster_sample_rate;
        
        %%%%
        
        %%%% TESTING CTX->WHEEL VELOCITY
        % Resample wheel, predict both velocity and speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Shove the k's into the master to save small V's
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%%%
        
        %%% Wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        % Pick trials to use
        use_trials = true(size(stimIDs));
        
        % Get stim info
        D = struct;
        D.stimulus = stimIDs;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);
        
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except curr_animal ...
        n_aligned_depths regression_params animals t ...
        fluor_all ...
        mua_all ...
        mua_ctxpred_all ...
        wheel_ctxpred_all ...
        ctx_str_k_all ...
        ctx_wheel_k_all ...
        wheel_all ...
        D_all
    
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    D_all

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_passive_fullscreen_trained_DECONVTEST'];
save([save_path filesep save_fn],'-v7.3');



%% Passive fullscreen trial activity (naive)

clear all
disp('Passive fullscreen trial activity (naive)')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP032','AP033','AP034','AP035','AP036'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);
wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    protocol = 'stimKalatsky';
    experiments = AP_find_experiments(animal,protocol);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Align (only to stim)
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf,[],2)',time_bin_centers)';
        
        %%%% TESTING: DECONV
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        %%%%
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        %         [~,predicted_spikes_std,explained_var] = ...
        %             AP_regresskernel(dfVdf_resample, ...
        %             binned_spikes_std,kernel_frames,lambda, ...
        %             regression_params.zs,regression_params.cvfold, ...
        %             false,regression_params.use_constant);
        
        %%%% TESTING DECONV
        [ctx_str_k,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Shove the k's into the master to save small V's
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
        
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        ctxpred_spikes = (ctxpred_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,ctxpred_spikes',t_peri_event)./raster_sample_rate;
        
        %%%%
        
        %%%% TESTING CTX->WHEEL VELOCITY
        % Resample wheel, predict both velocity and speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Shove the k's into the master to save small V's
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%%%
        
        %%% Wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        % Pick trials to use
        use_trials = true(size(stimIDs));
        
        % Get stim info
        D = struct;
        D.stimulus = stimIDs;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);
        
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except curr_animal ...
        n_aligned_depths regression_params animals t ...
        fluor_all ...
        mua_all ...
        mua_ctxpred_all ...
        wheel_ctxpred_all ...
        ctx_str_k_all ...
        ctx_wheel_k_all ...
        wheel_all ...
        D_all
    
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    D_all

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_passive_fullscreen_naive_DECONVTEST'];
save([save_path filesep save_fn],'-v7.3');



%% Regress task events to cortex/striatum activity (all concatenated)

data_fn = ['trial_activity_choiceworld_DECONVTEST'];
exclude_data = true;

AP_load_concat_normalize_ctx_str;

n_vs = size(fluor_allcat_deconv,3);
n_depths = size(mua_allcat,3);

% Get trial information
trial_contrast_allcat = max(D_allcat.stimulus,[],2);
[~,side_idx] = max(D_allcat.stimulus > 0,[],2);
trial_side_allcat = (side_idx-1.5)*2;
trial_choice_allcat = -(D_allcat.response-1.5)*2;

% Get timing
sample_rate = 1./mean(diff(t));

% Get reaction time
[move_trial,move_idx] = max(abs(wheel_allcat) > 0.025,[],2);
move_t = t(move_idx)';
move_t(~move_trial) = Inf;

% Get maximum wheel velocity in chosen direction
wheel_velocity_allcat = [zeros(size(wheel_allcat,1),1,1),diff(wheel_allcat,[],2)];
[max_speed,max_vel_idx] = max(abs(wheel_velocity_allcat(:,:,1).* ...
    (bsxfun(@times,wheel_velocity_allcat(:,:,1),trial_choice_allcat) > 0)),[],2);
vel_t = t(max_vel_idx);

max_vel = max_speed.*trial_choice_allcat;

%%% Regression (separate regressors to get partial explained variance)

% Stim regressors
contrasts = unique(trial_contrast_allcat(trial_contrast_allcat > 0));
contrastsides = sort([-contrasts;contrasts]);

stim_regressors = zeros(size(wheel_allcat,1),size(wheel_allcat,2),length(contrastsides));
for curr_condition_idx = 1:length(contrastsides)
    stim_regressors(trial_contrast_allcat.*trial_side_allcat == ...
        contrastsides(curr_condition_idx),find(t > 0,1),curr_condition_idx) = 1;
end

% Move onset regressors (L/R)
move_onset_regressors = zeros(size(wheel_allcat,1),size(wheel_allcat,2),2);
for curr_trial = find(move_trial)'
    
    % To use binary
    if trial_choice_allcat(curr_trial) == -1
        move_onset_regressors(curr_trial,move_idx(curr_trial),1) = 1;
    elseif trial_choice_allcat(curr_trial) == 1
        move_onset_regressors(curr_trial,move_idx(curr_trial),2) = 1;
    end
    
    %     % To fold in maximum velocity in chosen direction
    %     if trial_choice_allcat(curr_trial) == -1
    %         move_onset_regressors(curr_trial,move_idx(curr_trial),1) = abs(max_vel(curr_trial));
    %     elseif trial_choice_allcat(curr_trial) == 1
    %         move_onset_regressors(curr_trial,move_idx(curr_trial),2) = abs(max_vel(curr_trial));
    %     end
    
end

% Move ongoing regressor - (0.5s buffer from move onset, L/R)
wheel_vel = [zeros(size(wheel_allcat,1),1),diff(wheel_allcat,[],2)];
move_ongoing_regressors = zeros(size(wheel_allcat,1),size(wheel_allcat,2),2);
move_ongoing_regressors(:,:,1) = abs(wheel_allcat) > 2 & wheel_vel < -2;
move_ongoing_regressors(:,:,2) = abs(wheel_allcat) > 2 & wheel_vel > 2;

move_onset_buffer_t = 0.5;
move_onset_buffer_samples = ones(1,round(move_onset_buffer_t*sample_rate));
move_onset_buffer = padarray( ...
    convn(any(move_onset_regressors,3),move_onset_buffer_samples,'valid'), ...
    [0,length(move_onset_buffer_samples)-1],0,'pre');

move_ongoing_regressors = +(move_ongoing_regressors & ~move_onset_buffer);

% % Wheel regressors - separate leftward/rightward velocity
% wheel_vel_norm = [zeros(size(wheel_allcat,1),1),diff(wheel_allcat(:,:,1),[],2)];
% wheel_vel_norm = wheel_vel_norm/prctile(abs(wheel_vel_norm(:)),95);
% wheel_vel_norm(abs(wheel_vel_norm) > 1) = sign(wheel_vel_norm(abs(wheel_vel_norm) > 1));
% wheel_regressors = abs(repmat(wheel_vel_norm,1,1,2).*cat(3,wheel_vel_norm < 0,wheel_vel_norm > 0));

% Go cue regressors - separate for early/late move
go_cue_regressors = zeros(size(wheel_allcat,1),size(wheel_allcat,2));
go_cue_regressors(move_t <= 0.5,find(t > 0.5,1),1) = 1;
go_cue_regressors(move_t > 0.5,find(t > 0.5,1),2) = 1;

% Outcome regressors
% (using signals timing - not precise but looks good)
outcome_regressors = outcome_allcat;

% Day regressors
% (allow a separate offset for each day)
day_trials = cellfun(@(x) size(x,1),vertcat(fluor_all{:}));
if sum(day_trials) ~= size(fluor_allcat,1)
    error('Incorrect day trial count');
end
day_regressors = zeros(sum(day_trials),length(t),length(day_trials));
day_trials_cumulative = [0;cumsum(day_trials)];
for curr_day = 1:length(day_trials)
    day_regressors((day_trials_cumulative(curr_day)+1): ...
        day_trials_cumulative(curr_day+1),:,curr_day) = 1;
end

% Time in trial
trial_time_regressors = zeros(size(wheel_allcat));
trial_time_regressors(:,find(t > 0,1)) = 1;

% Concatenate regressors, set parameters
task_regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
task_regressor_labels = {'Stim','Move onset','Go cue','Outcome'};

t_shifts = {[0,0.5]; ... % stim
    [-0.5,1]; ... % move
    [-0.1,0.5]; ... % go cue
    [-0.5,1]; ... % outcome
    };

sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
    round(x(2)*(sample_rate)),t_shifts,'uni',false);
lambda = 0;
zs = [false,false];
cvfold = 5;
use_constant = false;
return_constant = false;


%%% Do regression on fluor
disp('Regressing task to cortex');
fluor_allcat_predicted = nan(size(fluor_allcat_deconv));
fluor_allcat_predicted_reduced = ...
    repmat(nan(size(fluor_allcat_deconv)),1,1,1,length(task_regressors));
fluor_kernel = cell(length(task_regressors)+return_constant,n_vs);
for curr_v = 1:n_vs
    
    baseline = nanmean(reshape(fluor_allcat_deconv(:,t < 0,curr_v),[],1));
    activity = fluor_allcat_deconv(:,:,curr_v) - baseline;
    
    activity_reshape = reshape(activity',[],1)';
    regressors_reshape = cellfun(@(x) ...
        reshape(permute(x(:,:,:),[2,1,3]),[],size(x,3))',task_regressors,'uni',false);
    
    warning off
    [kernel,activity_predicted_reshape,expl_var,activity_predicted_reduced_reshape] = AP_regresskernel( ...
        cellfun(@(x) x(:,~isnan(activity_reshape)),regressors_reshape,'uni',false), ...
        activity_reshape(~isnan(activity_reshape)),sample_shifts,lambda,zs,cvfold,return_constant,use_constant);
    warning on
    
    % Reshape full predictions
    activity_predicted_transpose = permute(nan(size(activity)),[2,1]);
    activity_predicted_transpose(~isnan(activity')) = activity_predicted_reshape;
    
    % Reshape reduced predictions
    activity_predicted_reduced_transpose = permute(repmat(nan(size(activity)),1,1,length(task_regressors)),[2,1,3]);
    activity_predicted_reduced_transpose(repmat(~isnan(activity'),1,1,length(task_regressors))) = ...
        activity_predicted_reduced_reshape;
    
    % Store predicted/reduced predicted/kernels
    fluor_allcat_predicted(:,:,curr_v) = ...
        permute(activity_predicted_transpose,[2,1,3]);
    fluor_allcat_predicted_reduced(:,:,curr_v,:) = ...
        permute(activity_predicted_reduced_transpose,[2,1,4,3]);
    fluor_kernel(:,curr_v) = kernel;
    
    AP_print_progress_fraction(curr_v,n_vs);
    
end

%%% Do regression on MUA
disp('Regressing task to striatum');
mua_allcat_predicted = nan(size(mua_allcat));
mua_allcat_predicted_reduced = ...
    repmat(nan(size(mua_allcat)),1,1,1,length(task_regressors));
mua_kernel = cell(length(task_regressors)+return_constant,n_depths);
for curr_depth = 1:n_depths
    
    % Convert MUA to single for GPU memory (fluorescence is already single)
    baseline = nanmean(reshape(mua_allcat(:,t < 0,curr_depth),[],1));
    activity = single(mua_allcat(:,:,curr_depth)) - baseline;
    
    activity_reshape = reshape(activity',[],1)';
    regressors_reshape = cellfun(@(x) ...
        reshape(permute(x(:,:,:),[2,1,3]),[],size(x,3))',task_regressors,'uni',false);
    
    [kernel,activity_predicted_reshape,expl_var,activity_predicted_reduced_reshape] = AP_regresskernel( ...
        cellfun(@(x) x(:,~isnan(activity_reshape)),regressors_reshape,'uni',false), ...
        activity_reshape(~isnan(activity_reshape)),sample_shifts,lambda,zs,cvfold,return_constant,use_constant);
    
    % Reshape full predictions
    activity_predicted_transpose = permute(nan(size(activity)),[2,1]);
    activity_predicted_transpose(~isnan(activity')) = activity_predicted_reshape;
    
    % Reshape reduced predictions
    activity_predicted_reduced_transpose = permute(repmat(nan(size(activity)),1,1,length(task_regressors)),[2,1,3]);
    activity_predicted_reduced_transpose(repmat(~isnan(activity'),1,1,length(task_regressors))) = ...
        activity_predicted_reduced_reshape;
    
    % Store predicted/reduced predicted/kernels
    mua_allcat_predicted(:,:,curr_depth) = ...
        permute(activity_predicted_transpose,[2,1,3]);
    mua_allcat_predicted_reduced(:,:,curr_depth,:) = ...
        permute(activity_predicted_reduced_transpose,[2,1,4,3]);
    mua_kernel(:,curr_depth) = kernel;
    
    AP_print_progress_fraction(curr_depth,n_depths);
    
end

%%% Do regression on ctx-predicted striatum
disp('Regressing task to cortex-predicted striatum');
predicted_mua_allcat_predicted = nan(size(mua_allcat));
predicted_mua_allcat_predicted_reduced = ...
    repmat(nan(size(mua_allcat)),1,1,1,length(task_regressors));
predicted_mua_kernel = cell(length(task_regressors)+return_constant,n_depths);
for curr_depth = 1:n_depths
    
    % Convert MUA to single for GPU memory (fluorescence is already single)
    baseline = nanmean(reshape(mua_ctxpred_allcat(:,t < 0,curr_depth),[],1));
    activity = single(mua_ctxpred_allcat(:,:,curr_depth)) - baseline;
    
    activity_reshape = reshape(activity',[],1)';
    regressors_reshape = cellfun(@(x) ...
        reshape(permute(x(:,:,:),[2,1,3]),[],size(x,3))',task_regressors,'uni',false);
    
    [kernel,activity_predicted_reshape,expl_var,activity_predicted_reduced_reshape] = AP_regresskernel( ...
        cellfun(@(x) x(:,~isnan(activity_reshape)),regressors_reshape,'uni',false), ...
        activity_reshape(~isnan(activity_reshape)),sample_shifts,lambda,zs,cvfold,return_constant);
    
    % Reshape full predictions
    activity_predicted_transpose = permute(nan(size(activity)),[2,1]);
    activity_predicted_transpose(~isnan(activity')) = activity_predicted_reshape;
    
    % Reshape reduced predictions
    activity_predicted_reduced_transpose = permute(repmat(nan(size(activity)),1,1,length(task_regressors)),[2,1,3]);
    activity_predicted_reduced_transpose(repmat(~isnan(activity'),1,1,length(task_regressors))) = ...
        activity_predicted_reduced_reshape;
    
    % Store predicted/reduced predicted/kernels
    predicted_mua_allcat_predicted(:,:,curr_depth) = ...
        permute(activity_predicted_transpose,[2,1,3]);
    predicted_mua_allcat_predicted_reduced(:,:,curr_depth,:) = ...
        permute(activity_predicted_reduced_transpose,[2,1,4,3]);
    predicted_mua_kernel(:,curr_depth) = kernel;
    
    AP_print_progress_fraction(curr_depth,n_depths);
    
end



% Save regression results
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_task_regression'];
save([save_path filesep save_fn], ...
    'regressors','regressor_labels','t_shifts', ...
    'fluor_allcat_predicted','fluor_allcat_predicted_reduced','fluor_kernel', ...
    'mua_allcat_predicted','mua_allcat_predicted_reduced','mua_kernel', ...
    'predicted_mua_allcat_predicted','predicted_mua_allcat_predicted_reduced','predicted_mua_kernel','-v7.3');
disp('Saved regression');


%% Cortex -> kernel-aligned striatum maps (deriv, concat experiments)
clear all
disp('Cortex -> kernel-aligned striatum regression');

% Parameters for regression
n_aligned_depths = 4;
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 60;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.2,0.2];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029', ...
    'AP032','AP033','AP034','AP035','AP036'};

batch_vars = struct;
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    % Find experiments
    % (use only behavior days because cortical recordings afterwards)
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    days = {experiments([experiments.imaging] & [experiments.ephys]).day};
    if isempty(experiments)
        % (if no behavior days then it was a naive mouse - use passive expt)
        protocol = 'AP_choiceWorldStimPassive';
        experiments = AP_find_experiments(animal,protocol);
        days = {experiments([experiments.imaging] & [experiments.ephys]).day};
    end
    
    disp(animal);
    
    for curr_day = 1:length(days)
        
        day = days{curr_day};
        
        % Find all protocols for that day, use non-multiples
        protocols = AP_list_experiments(animal,day);
        curr_experiments = [protocols(~[protocols.multiple]).experiment];
        
        % Loop through experiments, collate data
        time_bin_centers_all = cell(size(curr_experiments));
        dfVdf_resample_all = cell(size(curr_experiments));
        binned_spikes_all = cell(size(curr_experiments));
        
        for curr_exp = 1:length(curr_experiments)
            experiment = curr_experiments(curr_exp);
            AP_load_experiment;
            
            % Get time points to query
            sample_rate = framerate*regression_params.upsample_factor;
            time_bins = frame_t(find(frame_t > ...
                regression_params.skip_seconds,1)):1/sample_rate: ...
                frame_t(find(frame_t-frame_t(end) < ...
                -regression_params.skip_seconds,1,'last'));
            time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
            time_bin_centers_all{curr_exp} = time_bin_centers;
            
            % Get upsampled dVdf's
            dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
                diff(fVdf(regression_params.use_svs,:),[],2)',time_bin_centers)';
            dfVdf_resample_all{curr_exp} = dfVdf_resample;
            
            % Get striatum depth group by across-experiment alignment
            n_depths = n_aligned_depths;
            depth_group = aligned_str_depth_group;
            
            binned_spikes = zeros(n_depths,length(time_bin_centers));
            for curr_depth = 1:n_depths
                curr_spike_times = spike_times_timeline(depth_group == curr_depth);
                binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
            end
            
            binned_spikes_all{curr_exp} = binned_spikes;
        end
        
        % Concatenate all data
        time_bin_centers = cat(2,time_bin_centers_all{:});
        dfVdf_resample = cat(2,dfVdf_resample_all{:});
        binned_spikes = cat(2,binned_spikes_all{:});
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        % Regress fluorescence to spikes
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;
        
        [k,ctxpred_spikes,explained_var] = ...
            AP_regresskernel(dfVdf_resample, ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant);
        
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        k_px = zeros(size(Udf_aligned,1),size(Udf_aligned,2),size(k,2),size(k,3),'single');
        for curr_spikes = 1:size(k,3)
            k_px(:,:,:,curr_spikes) = ...
                svdFrameReconstruct(Udf_aligned(:,:,regression_params.use_svs),k(:,:,curr_spikes));
        end
        
        % NaN-out depths with no spikes
        k_px(:,:,:,~any(binned_spikes,2)) = NaN;
        
        % Store kernel pixels
        batch_vars(curr_animal).animal = animal;
        batch_vars(curr_animal).regression_params = regression_params;
        batch_vars(curr_animal).day{curr_day} = day;
        batch_vars(curr_animal).t{curr_day} = kernel_frames/sample_rate;
        batch_vars(curr_animal).k_px{curr_day} = k_px;
        batch_vars(curr_animal).explained_var{curr_day} = explained_var.total;
        
        AP_print_progress_fraction(curr_day,length(days));
        clearvars -except n_aligned_depths regression_params animals animal curr_animal days curr_day batch_vars
    end
    disp(['Finished ' animal]);
end

% Save
save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_ephys'];
save([save_path filesep 'wf_ephys_maps_concat_' num2str(n_aligned_depths) '_depths_kernel'],'batch_vars','-v7.3');
warning('saving -v7.3');
disp('Finished batch');



%% Choiceworld trial activity (old, derivative? and other stuff)
clear all
disp('Choiceworld trial activity')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 60;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.2,0.2];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

fluor_all = cell(length(animals),1);
mua_all = cell(length(animals),1);
predicted_mua_std_all = cell(length(animals),1);
wheel_all = cell(length(animals),1);
reward_all = cell(length(animals),1);
D_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    fluor_all{curr_animal} = cell(length(experiments),1);
    mua_all{curr_animal} = cell(length(experiments),1);
    predicted_mua_std_all{curr_animal} = cell(length(experiments),1);
    wheel_all{curr_animal} = cell(length(experiments),1);
    D_all{curr_animal} = cell(length(experiments),1);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,3];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Align (only to stim)
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf(regression_params.use_svs,:),[],2)',time_bin_centers)';
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;
        
        [~,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(dfVdf_resample, ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant);
        
        event_aligned_predicted_mua_std = ...
            interp1(time_bin_centers,ctxpred_spikes_std',t_peri_event);
        
        %%% Wheel
        event_aligned_wheel_raw = interp1(Timeline.rawDAQTimestamps, ...
            wheel_position,t_peri_event);
        event_aligned_wheel = bsxfun(@minus,event_aligned_wheel_raw, ...
            nanmedian(event_aligned_wheel_raw(:,t < 0),2));
        
        %%% Reward
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        event_aligned_reward = (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_peri_event_bins(x,:)), ...
            [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate) > 0;
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        % D.response = (trial_choice(use_trials)'+1)/2+1;
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);
        
        % Store activity and stim
        fluor_all{curr_animal}{curr_day} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal}{curr_day} = event_aligned_mua(use_trials,:,:,:);
        predicted_mua_std_all{curr_animal}{curr_day} = event_aligned_predicted_mua_std(use_trials,:,:,:);
        wheel_all{curr_animal}{curr_day} = event_aligned_wheel(use_trials,:,:);
        reward_all{curr_animal}{curr_day} = event_aligned_reward(use_trials,:,:);
        D_all{curr_animal}{curr_day} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except n_aligned_depths regression_params animals curr_animal ...
        t fluor_all mua_all predicted_mua_std_all wheel_all reward_all D_all
    
end
clearvars -except n_aligned_depths t fluor_all mua_all predicted_mua_std_all wheel_all reward_all D_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld'];
save([save_path filesep save_fn],'-v7.3');


%% Pull out and save activity in each trial (choiceworld - frame sampling rate)
clear all
disp('Getting trial activity choiceworld (at framerate)')

n_aligned_depths = 4;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

fluor_all = cell(length(animals),1);
mua_all = cell(length(animals),1);
wheel_all = cell(length(animals),1);
reward_all = cell(length(animals),1);
D_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    fluor_all{curr_animal} = cell(length(experiments),1);
    mua_all{curr_animal} = cell(length(experiments),1);
    wheel_all{curr_animal} = cell(length(experiments),1);
    D_all{curr_animal} = cell(length(experiments),1);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,3];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Align (only to stim)
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        % Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        % MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        % Wheel
        event_aligned_wheel_raw = interp1(Timeline.rawDAQTimestamps, ...
            wheel_position,t_peri_event);
        event_aligned_wheel = bsxfun(@minus,event_aligned_wheel_raw, ...
            nanmedian(event_aligned_wheel_raw(:,t < 0),2));
        
        % Reward
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        event_aligned_reward = (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_peri_event_bins(x,:)), ...
            [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate) > 0;
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        % D.response = (trial_choice(use_trials)'+1)/2+1;
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);
        
        % Store activity and stim
        fluor_all{curr_animal}{curr_day} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal}{curr_day} = event_aligned_mua(use_trials,:,:,:);
        wheel_all{curr_animal}{curr_day} = event_aligned_wheel(use_trials,:,:);
        reward_all{curr_animal}{curr_day} = event_aligned_reward(use_trials,:,:);
        D_all{curr_animal}{curr_day} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except n_aligned_depths animals curr_animal t fluor_all mua_all wheel_all reward_all D_all
    
end
clearvars -except n_aligned_depths t fluor_all mua_all wheel_all reward_all D_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_framerate'];
save([save_path filesep save_fn],'-v7.3');


%% Pull out and save activity in each trial (choiceworld)
clear all
disp('Getting trial activity choiceworld')

n_aligned_depths = 4;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

fluor_all = cell(length(animals),1);
mua_all = cell(length(animals),1);
wheel_all = cell(length(animals),1);
reward_all = cell(length(animals),1);
D_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    fluor_all{curr_animal} = cell(length(experiments),1);
    mua_all{curr_animal} = cell(length(experiments),1);
    wheel_all{curr_animal} = cell(length(experiments),1);
    D_all{curr_animal} = cell(length(experiments),1);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,3];
        upsample_factor = 3;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Align (only to stim)
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        % Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        % MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        % Wheel
        event_aligned_wheel_raw = interp1(Timeline.rawDAQTimestamps, ...
            wheel_position,t_peri_event);
        event_aligned_wheel = bsxfun(@minus,event_aligned_wheel_raw, ...
            nanmedian(event_aligned_wheel_raw(:,t < 0),2));
        
        % Reward
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        event_aligned_reward = (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_peri_event_bins(x,:)), ...
            [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate) > 0;
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        % D.response = (trial_choice(use_trials)'+1)/2+1;
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);
        
        % Store activity and stim
        fluor_all{curr_animal}{curr_day} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal}{curr_day} = event_aligned_mua(use_trials,:,:,:);
        wheel_all{curr_animal}{curr_day} = event_aligned_wheel(use_trials,:,:);
        reward_all{curr_animal}{curr_day} = event_aligned_reward(use_trials,:,:);
        D_all{curr_animal}{curr_day} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except n_aligned_depths animals curr_animal t fluor_all mua_all wheel_all reward_all D_all
    
end
clearvars -except n_aligned_depths t fluor_all mua_all wheel_all reward_all D_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld'];
save([save_path filesep save_fn],'-v7.3');


%% Choiceworld trial activity (old, not including task regression)

clear all
disp('Choiceworld trial activity')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);
mua_taskpred_k_all = cell(1,1);
mua_taskpred_all = cell(1,1);
mua_taskpred_reduced_all = cell(1,1);
mua_ctxpred_taskpred_k_all = cell(1,1);
mua_ctxpred_taskpred_all = cell(1,1);
mua_ctxpred_taskpred_reduced_all = cell(1,1);
wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Align (only to stim)
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf,[],2)',time_bin_centers)';
        
        %%%% TESTING: DECONV
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        %%%%
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        %         [~,predicted_spikes_std,explained_var] = ...
        %             AP_regresskernel(dfVdf_resample, ...
        %             binned_spikes_std,kernel_frames,lambda, ...
        %             regression_params.zs,regression_params.cvfold, ...
        %             false,regression_params.use_constant);
        
        %%%% TESTING DECONV
        [ctx_str_k,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Shove the k's into the master to save small V's
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
        
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        ctxpred_spikes = (ctxpred_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,ctxpred_spikes',t_peri_event)./raster_sample_rate;
        
        %%%%
        
        %%%% TESTING CTX->WHEEL VELOCITY
        % Resample wheel, predict both velocity and speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Shove the k's into the master to save small V's
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%%%
        
        %%% Wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_peri_event_bins(x,:)), ...
            find(trial_outcome == 1)','uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_peri_event_bins(x,:)), ...
            find(trial_outcome == -1)','uni',false))) > 0;
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);
        
        D.outcome = trial_outcome(use_trials);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING: TASK REGRESSION
        
        % Get reaction time
        [move_trial,move_idx] = max(abs(event_aligned_wheel(use_trials,:)) > 0.025,[],2);
        move_t = t(move_idx)';
        move_t(~move_trial) = Inf;
        
        %%% Regression (separate regressors to get partial explained variance)
        
        % Stim regressors
        contrastsides = diff(D.stimulus,[],2);
        unique_contrastsides = sort(unique(contrastsides(contrastsides ~= 0)));
        
        stim_regressors = zeros(sum(use_trials),length(t),length(unique_contrastsides));
        for curr_condition_idx = 1:length(unique_contrastsides)
            stim_regressors(contrastsides == ...
                contrastsides(curr_condition_idx),find(t > 0,1),curr_condition_idx) = 1;
        end
        
        % Move onset regressors (L/R)
        move_onset_regressors = zeros(sum(use_trials),length(t),2);
        for curr_trial = 1:sum(use_trials)
            % To use binary
            if D.response(curr_trial) == 1
                move_onset_regressors(curr_trial,move_idx(curr_trial),1) = 1;
            elseif D.response(curr_trial) == 2
                move_onset_regressors(curr_trial,move_idx(curr_trial),2) = 1;
            end
        end
        
        % Go cue regressors - separate for early/late move
        go_cue_regressors = zeros(sum(use_trials),length(t),2);
        go_cue_regressors(move_t <= 0.5,find(t > 0.5,1),1) = 1;
        go_cue_regressors(move_t > 0.5,find(t > 0.5,1),2) = 1;
        
        % Outcome regressors
        outcome_regressors = event_aligned_outcome(use_trials,:,:);
        
        task_regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        task_regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
        
        t_shifts = {[0,1]; ... % stim
            [-0.5,1]; ... % move onset
            [0,0.5]; ... % go cue
            [0,1]}; % outcome
        
        sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = true;
        return_constant = true;
        
        % Regression task->MUA
        use_mua_task_regress = event_aligned_mua(use_trials,:,:);
        mua_taskpred_k = cell(length(task_regressors)+return_constant,n_depths);
        mua_taskpred = nan(size(use_mua_task_regress));
        mua_taskpred_reduced = ...
            repmat(nan(size(use_mua_task_regress)),1,1,1,length(task_regressors));
        for curr_depth = 1:n_depths
            
            activity = single(use_mua_task_regress(:,:,curr_depth));
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
            
            activity_reshape = reshape(activity',[],1)';
            regressors_reshape = cellfun(@(x) ...
                reshape(permute(x,[2,1,3]),[],size(x,3))',task_regressors,'uni',false);
            
            [task_kernel,activity_predicted_reshape,expl_var,activity_predicted_reduced_reshape] = AP_regresskernel( ...
                cellfun(@(x) x(:,~isnan(activity_reshape)),regressors_reshape,'uni',false), ...
                activity_reshape(~isnan(activity_reshape)),sample_shifts,lambda,zs,cvfold,return_constant,use_constant);
            
            % Reshape full predictions
            activity_predicted_transpose = permute(nan(size(activity)),[2,1]);
            activity_predicted_transpose(~isnan(activity')) = activity_predicted_reshape;
            
            % Reshape reduced predictions
            activity_predicted_reduced_transpose = permute(repmat(nan(size(activity)),1,1,length(task_regressors)),[2,1,3]);
            activity_predicted_reduced_transpose(repmat(~isnan(activity'),1,1,length(task_regressors))) = ...
                activity_predicted_reduced_reshape;
            
            % Store predicted/reduced predicted/kernels
            mua_taskpred_k(:,curr_depth) = task_kernel;
            mua_taskpred(:,:,curr_depth) = ...
                permute(activity_predicted_transpose,[2,1,3]);
            mua_taskpred_reduced(:,:,curr_depth,:) = ...
                permute(activity_predicted_reduced_transpose,[2,1,4,3]);
        end
        
        % Regression task->MUA-ctxpred
        use_mua_task_regress = event_aligned_mua_ctxpred(use_trials,:,:);
        mua_ctxpred_taskpred_k = cell(length(task_regressors)+return_constant,n_depths);
        mua_ctxpred_taskpred = nan(size(use_mua_task_regress));
        mua_ctxpred_taskpred_reduced = ...
            repmat(nan(size(use_mua_task_regress)),1,1,1,length(task_regressors));
        for curr_depth = 1:n_depths
            
            activity = single(use_mua_task_regress(:,:,curr_depth));
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
            
            activity_reshape = reshape(activity',[],1)';
            regressors_reshape = cellfun(@(x) ...
                reshape(permute(x,[2,1,3]),[],size(x,3))',task_regressors,'uni',false);
            
            [task_kernel,activity_predicted_reshape,expl_var,activity_predicted_reduced_reshape] = AP_regresskernel( ...
                cellfun(@(x) x(:,~isnan(activity_reshape)),regressors_reshape,'uni',false), ...
                activity_reshape(~isnan(activity_reshape)),sample_shifts,lambda,zs,cvfold,return_constant,use_constant);
            
            % Reshape full predictions
            activity_predicted_transpose = permute(nan(size(activity)),[2,1]);
            activity_predicted_transpose(~isnan(activity')) = activity_predicted_reshape;
            
            % Reshape reduced predictions
            activity_predicted_reduced_transpose = permute(repmat(nan(size(activity)),1,1,length(task_regressors)),[2,1,3]);
            activity_predicted_reduced_transpose(repmat(~isnan(activity'),1,1,length(task_regressors))) = ...
                activity_predicted_reduced_reshape;
            
            % Store predicted/reduced predicted/kernels
            mua_ctxpred_taskpred_k(:,curr_depth) = task_kernel;
            mua_ctxpred_taskpred(:,:,curr_depth) = ...
                permute(activity_predicted_transpose,[2,1,3]);
            mua_ctxpred_taskpred_reduced(:,:,curr_depth,:) = ...
                permute(activity_predicted_reduced_transpose,[2,1,4,3]);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        mua_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_taskpred_k;
        mua_taskpred_all{curr_animal,1}{curr_day,1} = mua_taskpred;
        mua_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_taskpred_reduced;
        
        mua_ctxpred_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_k;
        mua_ctxpred_taskpred_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred;
        mua_ctxpred_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_reduced;
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except curr_animal ...
        n_aligned_depths regression_params animals t ...
        fluor_all ...
        mua_all ...
        mua_ctxpred_all ...
        mua_taskpred_k_all ...
        mua_taskpred_all ...
        mua_taskpred_reduced_all ...
        mua_ctxpred_taskpred_k_all ...
        mua_ctxpred_taskpred_all ...
        mua_ctxpred_taskpred_reduced_all ...
        wheel_ctxpred_all ...
        ctx_str_k_all ...
        ctx_wheel_k_all ...
        wheel_all ...
        outcome_all ...
        D_all
    
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    outcome_all ...
    D_all

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_DECONVTEST'];
save([save_path filesep save_fn],'-v7.3');


%% Passive fullscreen trial activity (trained)

clear all
disp('Passive fullscreen trial activity (trained)')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 60;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.2,0.2];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

fluor_all = cell(length(animals),1);
mua_all = cell(length(animals),1);
predicted_mua_std_all = cell(length(animals),1);
wheel_all = cell(length(animals),1);
D_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    % Only use days with choiceworld (sometimes recorded in cortex, no bhv)
    protocol = 'vanillaChoiceworld';
    behavior_experiments = AP_find_experiments(animal,protocol);
    
    protocol = 'stimKalatsky';
    passive_experiments = AP_find_experiments(animal,protocol);
    
    behavior_day = ismember({passive_experiments.day},{behavior_experiments.day});
    
    experiments = passive_experiments([passive_experiments.imaging] & [passive_experiments.ephys] & behavior_day);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    fluor_all{curr_animal} = cell(length(experiments),1);
    mua_all{curr_animal} = cell(length(experiments),1);
    predicted_mua_std_all{curr_animal} = cell(length(experiments),1);
    wheel_all{curr_animal} = cell(length(experiments),1);
    D_all{curr_animal} = cell(length(experiments),1);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        % Set components to use
        use_components = 1:200;
        
        % Group multiunit by depth
        % (evenly across recorded striatum)
        %         n_depths = 6;
        %         depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        %         depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
        %         depth_group = discretize(spike_depths,depth_group_edges);
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,3];
        raster_sample_rate = 1/(framerate*regression_params.upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        use_align = stimOn_times;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            % (for all spikes in depth group)
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf(regression_params.use_svs,:),[],2)',time_bin_centers)';
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;
        
        [~,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(dfVdf_resample, ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant);
        
        event_aligned_predicted_mua_std = ...
            interp1(time_bin_centers,ctxpred_spikes_std',t_peri_event);
        
        %%% Wheel
        event_aligned_wheel_raw = interp1(Timeline.rawDAQTimestamps, ...
            wheel_position,t_peri_event);
        event_aligned_wheel = bsxfun(@minus,event_aligned_wheel_raw, ...
            nanmedian(event_aligned_wheel_raw(:,t < 0),2));
        
        % Get stim info
        D = struct;
        D.stimulus = stimIDs;
        
        % Store activity and stim
        fluor_all{curr_animal}{curr_day} = event_aligned_V;
        mua_all{curr_animal}{curr_day} = event_aligned_mua;
        predicted_mua_std_all{curr_animal}{curr_day} = event_aligned_predicted_mua_std;
        wheel_all{curr_animal}{curr_day} = event_aligned_wheel;
        D_all{curr_animal}{curr_day} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
    
    clearvars -except n_aligned_depths regression_params animals curr_animal ...
        t fluor_all mua_all predicted_mua_std_all wheel_all D_all
    
end
clearvars -except n_aligned_depths t fluor_all mua_all predicted_mua_std_all wheel_all D_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_passive_fullscreen'];
save([save_path filesep save_fn],'-v7.3');

%% Passive choiceworld trial activity (trained)

clear all
disp('Passive choiceworld trial activity (trained)')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 60;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.2,0.2];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

fluor_all = cell(length(animals),1);
mua_all = cell(length(animals),1);
predicted_mua_std_all = cell(length(animals),1);
wheel_all = cell(length(animals),1);
D_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    % Only use days with choiceworld (sometimes recorded in cortex, no bhv)
    protocol = 'vanillaChoiceworld';
    behavior_experiments = AP_find_experiments(animal,protocol);
    
    protocol = 'AP_choiceWorldStimPassive';
    passive_experiments = AP_find_experiments(animal,protocol);
    
    behavior_day = ismember({passive_experiments.day},{behavior_experiments.day});
    
    experiments = passive_experiments([passive_experiments.imaging] & [passive_experiments.ephys] & behavior_day);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    fluor_all{curr_animal} = cell(length(experiments),1);
    mua_all{curr_animal} = cell(length(experiments),1);
    predicted_mua_std_all{curr_animal} = cell(length(experiments),1);
    wheel_all{curr_animal} = cell(length(experiments),1);
    D_all{curr_animal} = cell(length(experiments),1);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        % Set components to use
        use_components = 1:200;
        
        % Group multiunit by depth
        % (evenly across recorded striatum)
        %         n_depths = 6;
        %         depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        %         depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
        %         depth_group = discretize(spike_depths,depth_group_edges);
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,3];
        raster_sample_rate = 1/(framerate*regression_params.upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        use_align = stimOn_times;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            % (for all spikes in depth group)
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf(regression_params.use_svs,:),[],2)',time_bin_centers)';
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;
        
        [~,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(dfVdf_resample, ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant);
        
        event_aligned_predicted_mua_std = ...
            interp1(time_bin_centers,ctxpred_spikes_std',t_peri_event);
        
        %%% Wheel
        event_aligned_wheel_raw = interp1(Timeline.rawDAQTimestamps, ...
            wheel_position,t_peri_event);
        event_aligned_wheel = bsxfun(@minus,event_aligned_wheel_raw, ...
            nanmedian(event_aligned_wheel_raw(:,t < 0),2));
        
        % Get stim info
        D = struct;
        D.stimulus = sign(conditions(stimIDs,1)).*conditions(stimIDs,2);
        
        % Store activity and stim
        fluor_all{curr_animal}{curr_day} = event_aligned_V;
        mua_all{curr_animal}{curr_day} = event_aligned_mua;
        predicted_mua_std_all{curr_animal}{curr_day} = event_aligned_predicted_mua_std;
        wheel_all{curr_animal}{curr_day} = event_aligned_wheel;
        D_all{curr_animal}{curr_day} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
    
    clearvars -except n_aligned_depths regression_params animals curr_animal ...
        t fluor_all mua_all predicted_mua_std_all wheel_all D_all
    
end
clearvars -except n_aligned_depths t fluor_all mua_all predicted_mua_std_all wheel_all D_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_passive_choiceworld'];
save([save_path filesep save_fn],'-v7.3');


%% Passive fullscreen trial activity (naive)

clear all
disp('Passive fullscreen trial activity (naive)')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 60;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.2,0.2];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP032','AP033','AP034','AP035','AP036'};

fluor_all = cell(length(animals),1);
mua_all = cell(length(animals),1);
predicted_mua_std_all = cell(length(animals),1);
wheel_all = cell(length(animals),1);
D_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    protocol = 'stimKalatsky';
    experiments = AP_find_experiments(animal,protocol);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    fluor_all{curr_animal} = cell(length(experiments),1);
    mua_all{curr_animal} = cell(length(experiments),1);
    predicted_mua_std_all{curr_animal} = cell(length(experiments),1);
    wheel_all{curr_animal} = cell(length(experiments),1);
    D_all{curr_animal} = cell(length(experiments),1);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        % Set components to use
        use_components = 1:200;
        
        % Group multiunit by depth
        % (evenly across recorded striatum)
        %         n_depths = 6;
        %         depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        %         depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
        %         depth_group = discretize(spike_depths,depth_group_edges);
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,3];
        raster_sample_rate = 1/(framerate*regression_params.upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        use_align = stimOn_times;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            % (for all spikes in depth group)
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf(regression_params.use_svs,:),[],2)',time_bin_centers)';
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;
        
        [~,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(dfVdf_resample, ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant);
        
        event_aligned_predicted_mua_std = ...
            interp1(time_bin_centers,ctxpred_spikes_std',t_peri_event);
        
        %%% Wheel
        event_aligned_wheel_raw = interp1(Timeline.rawDAQTimestamps, ...
            wheel_position,t_peri_event);
        event_aligned_wheel = bsxfun(@minus,event_aligned_wheel_raw, ...
            nanmedian(event_aligned_wheel_raw(:,t < 0),2));
        
        % Get stim info
        D = struct;
        D.stimulus = stimIDs;
        
        % Store activity and stim
        fluor_all{curr_animal}{curr_day} = event_aligned_V;
        mua_all{curr_animal}{curr_day} = event_aligned_mua;
        predicted_mua_std_all{curr_animal}{curr_day} = event_aligned_predicted_mua_std;
        wheel_all{curr_animal}{curr_day} = event_aligned_wheel;
        D_all{curr_animal}{curr_day} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
    
    clearvars -except n_aligned_depths regression_params animals curr_animal ...
        t fluor_all mua_all predicted_mua_std_all wheel_all D_all
    
end
clearvars -except n_aligned_depths t fluor_all mua_all predicted_mua_std_all wheel_all D_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_passive_fullscreen_naive'];
save([save_path filesep save_fn],'-v7.3');


%% Passive choiceworld trial activity (naive)
clear all
disp('Passive choiceworld trial activity (naive)')

n_aligned_depths = 4;

animals = {'AP032','AP033','AP034','AP035','AP036'};

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 60;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.2,0.2];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

fluor_all = cell(length(animals),1);
mua_all = cell(length(animals),1);
predicted_mua_std_all = cell(length(animals),1);
wheel_all = cell(length(animals),1);
D_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    protocol = 'AP_choiceWorldStimPassive';
    experiments = AP_find_experiments(animal,protocol);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    fluor_all{curr_animal} = cell(length(experiments),1);
    mua_all{curr_animal} = cell(length(experiments),1);
    predicted_mua_std_all{curr_animal} = cell(length(experiments),1);
    wheel_all{curr_animal} = cell(length(experiments),1);
    D_all{curr_animal} = cell(length(experiments),1);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        % Set components to use
        use_components = 1:200;
        
        % Group multiunit by depth
        % (evenly across recorded striatum)
        %         n_depths = 6;
        %         depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        %         depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
        %         depth_group = discretize(spike_depths,depth_group_edges);
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,3];
        raster_sample_rate = 1/(framerate*regression_params.upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        use_align = stimOn_times;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            % (for all spikes in depth group)
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf(regression_params.use_svs,:),[],2)',time_bin_centers)';
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;
        
        [~,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(dfVdf_resample, ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant);
        
        event_aligned_predicted_mua_std = ...
            interp1(time_bin_centers,ctxpred_spikes_std',t_peri_event);
        
        %%% Wheel
        event_aligned_wheel_raw = interp1(Timeline.rawDAQTimestamps, ...
            wheel_position,t_peri_event);
        event_aligned_wheel = bsxfun(@minus,event_aligned_wheel_raw, ...
            nanmedian(event_aligned_wheel_raw(:,t < 0),2));
        
        % Get stim info
        D = struct;
        D.stimulus = sign(conditions(stimIDs,1)).*conditions(stimIDs,2);
        
        % Store activity and stim
        fluor_all{curr_animal}{curr_day} = event_aligned_V;
        mua_all{curr_animal}{curr_day} = event_aligned_mua;
        predicted_mua_std_all{curr_animal}{curr_day} = event_aligned_predicted_mua_std;
        wheel_all{curr_animal}{curr_day} = event_aligned_wheel;
        D_all{curr_animal}{curr_day} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
    
    clearvars -except n_aligned_depths regression_params animals curr_animal ...
        t fluor_all mua_all predicted_mua_std_all wheel_all D_all
    
end
clearvars -except n_aligned_depths t fluor_all mua_all predicted_mua_std_all wheel_all D_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_passive_choiceworld_naive'];
save([save_path filesep save_fn],'-v7.3');

%% (lock off)
end

%% TEST

%% Choiceworld trial activity (striatum domain) - No WF deconv

clear all
disp('Choiceworld trial activity (striatum domain)')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);

mua_taskpred_k_all = cell(1,1);
mua_taskpred_all = cell(1,1);
mua_taskpred_reduced_all = cell(1,1);
mua_taskpred_expl_var_total_all = cell(1,1);
mua_taskpred_expl_var_partial_all = cell(1,1);

mua_ctxpred_taskpred_k_all = cell(1,1);
mua_ctxpred_taskpred_all = cell(1,1);
mua_ctxpred_taskpred_reduced_all = cell(1,1);
mua_ctxpred_taskpred_expl_var_total_all = cell(1,1);
mua_ctxpred_taskpred_expl_var_partial_all = cell(1,1);

fluor_taskpred_k_all = cell(1,1);
fluor_taskpred_all = cell(1,1);
fluor_taskpred_reduced_all = cell(1,1);
fluor_taskpred_expl_var_total_all = cell(1,1);
fluor_taskpred_expl_var_partial_all = cell(1,1);

wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
movement_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = AP_align_widefield(Udf,animal,day);
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Get align times
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        
        %%% Trial-align cortex
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% Trial-align striatum
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            if isempty(curr_spikes)
                continue
            end
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regress cortex to striatum
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Resample fluorescence (don't deconvolve)
        fVdf_resample = interp1(frame_t,fVdf',time_bin_centers)';
        
        % Deconvolve fluoresence
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        
        % Bin spikes across the experiment
        binned_spikes = nan(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            % Skip if no spikes at this depth
            if isempty(curr_spike_times)
                continue
            end
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        % Regress cortex to striatum
        [ctx_str_k,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Recast the k's into the master U
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
        
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        ctxpred_spikes = (ctxpred_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,ctxpred_spikes',t_peri_event)./raster_sample_rate;
        
        % Regress cortex to wheel velocity/speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Recast the k's into the master U
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%% Trial-align wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Trial-align facecam movement
        event_aligned_movement = interp1(facecam_t(~isnan(facecam_t)), ...
            frame_movement(~isnan(facecam_t)),t_peri_event);
        
        %%% Trial-align outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_peri_event_bins(x,:)), ...
            find(trial_outcome == 1),'uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_peri_event_bins(x,:)), ...
            find(trial_outcome == -1),'uni',false))) > 0;
        
        % Pick trials to keep
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials)' & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials)' == -1;
        R_trials = signals_events.trialSideValues(1:n_trials)' == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = 3-(abs((trial_choice(use_trials)+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%% Regress task to cortex/striatum/cortex-predicted striatum
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.025,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Build regressors (only a subset of these are used)
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times))'.* ...
            signals_events.trialContrastValues(1:length(stimOn_times))';
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end
        
        % Stim move regressors (one for each stim when it starts to move)
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
            
        end
        
        % Stim center regressors (one for each stim when it's stopped during reward)
        unique_contrasts = unique(contrasts(contrasts > 0));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ...
                abs(stim_contrastsides) == unique_contrasts(curr_contrast));
            
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
            
        end
        
        % Move onset regressors (L/R)
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move onset x stim regressors (one for each contrast/side)
        move_onset_stim_time_absolute = arrayfun(@(curr_stim) ...
            arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & stim_contrastsides == unique_stim(curr_stim))), ...
            1:length(unique_stim),'uni',false);
        
        move_onset_stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            move_onset_stim_regressors(curr_stim,:) = ...
                histcounts(move_onset_stim_time_absolute{curr_stim},time_bins);
        end
        
        % Move ongoing regressors (L/R choice for duration of movement)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        
        move_stopped_t = 0.5;
        move_stopped_samples = round(sample_rate*move_stopped_t);
        wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.025), ...
            ones(1,move_stopped_samples),'full') > 0;
        wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);
        
        move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_L_absolute','uni',false));
        move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_R_absolute','uni',false));
        
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
        move_ongoing_regressors(2,move_ongoing_R_samples) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        % (for go cue only on late move trials)
%         go_cue_regressors = histcounts( ...
%             signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        % (for go cue with early/late move trials)
        go_cue_regressors = zeros(1,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        % (regressors for hit only)
%         outcome_regressors = histcounts(reward_t_timeline,time_bins);
        % (regressors for both hit and miss)
        outcome_regressors = zeros(2,length(time_bin_centers));
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate selected regressors, set parameters
        
        task_regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        task_regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
                       
        task_t_shifts = { ...
            [0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [0,0.5]; ... % go cue
            [0,0.5]}; % outcome       
        
        % (old extended timings)
        %         task_t_shifts = {[0,0.5]; ... % stim
        %             [-0.5,1]; ... % move
        %             [-0.1,0.5]; ... % go cue
        %             [-0.5,1]}; % outcome

        % (to include stim x move)
        %         task_regressors = {stim_regressors;move_onset_regressors;move_onset_stim_regressors;go_cue_regressors;outcome_regressors};
        %         task_regressor_labels = {'Stim','Move','Stim x move','Go cue','Outcome'};
        %
        %         task_t_shifts = { ...
        %             [0,0.5]; ... % stim
        %             [-0.5,1]; ... % move
        %             [-0.5,1]; ... % stim x move
        %             [0,0.5]; ... % go cue
        %             [0,0.5]}; % outcome
        
        task_regressor_sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),task_t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;
        
        % Regression task -> MUA
        baseline = nanmean(reshape(event_aligned_mua(:,t < 0,:),[], ...
            size(event_aligned_mua,3))*raster_sample_rate,1)';
        activity = single(binned_spikes) - baseline;
        
        [mua_taskpred_k,mua_taskpred_long,mua_taskpred_expl_var,mua_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        mua_taskpred = ...
            interp1(time_bin_centers,mua_taskpred_long',t_peri_event)./raster_sample_rate;
        
        mua_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,mua_taskpred_reduced_long(:,:,x)', ...
            t_peri_event)./raster_sample_rate,permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        % Regression task -> MUA-ctxpred
        baseline = nanmean(reshape(event_aligned_mua_ctxpred(:,t < 0,:),[], ...
            size(event_aligned_mua_ctxpred,3))*raster_sample_rate,1)';
        activity = single(ctxpred_spikes) - baseline;
        
        [mua_ctxpred_taskpred_k,mua_ctxpred_taskpred_long,mua_ctxpred_taskpred_expl_var,mua_ctxpred_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        mua_ctxpred_taskpred = ...
            interp1(time_bin_centers,mua_ctxpred_taskpred_long',t_peri_event)./raster_sample_rate;
        
        mua_ctxpred_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,mua_ctxpred_taskpred_reduced_long(:,:,x)', ...
            t_peri_event)./raster_sample_rate,permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        % Regression task -> (master U, deconvolved) fluor
        event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
        fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);
        
        baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,:),[],size(event_aligned_V_deconv,3)))';
        activity = single(fVdf_deconv_resample_recast(use_components,:))-baseline;
        
        [fluor_taskpred_k,fluor_taskpred_long,fluor_taskpred_expl_var,fluor_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        fluor_taskpred = ...
            interp1(time_bin_centers,fluor_taskpred_long',t_peri_event);
        
        fluor_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,fluor_taskpred_reduced_long(:,:,x)', ...
            t_peri_event),permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        mua_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_taskpred_k;
        mua_taskpred_all{curr_animal,1}{curr_day,1} = mua_taskpred(use_trials,:,:,:);
        mua_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_taskpred_reduced(use_trials,:,:,:);
        mua_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = mua_taskpred_expl_var.total;
        mua_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = mua_taskpred_expl_var.partial;
        
        mua_ctxpred_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_k;
        mua_ctxpred_taskpred_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred(use_trials,:,:,:);
        mua_ctxpred_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_reduced(use_trials,:,:,:);
        mua_ctxpred_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_expl_var.total;
        mua_ctxpred_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_expl_var.partial;
        
        fluor_taskpred_k_all{curr_animal,1}{curr_day,1} = fluor_taskpred_k;
        fluor_taskpred_all{curr_animal,1}{curr_day,1} = fluor_taskpred(use_trials,:,:,:);
        fluor_taskpred_reduced_all{curr_animal,1}{curr_day,1} = fluor_taskpred_reduced(use_trials,:,:,:);
        fluor_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.total;
        fluor_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.partial;
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        movement_all{curr_animal,1}{curr_day,1} = event_aligned_movement(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars -except curr_animal animal protocol experiments load_parts curr_day ...
            n_aligned_depths regression_params animals t ...
            ...
            fluor_all ...
            mua_all ...
            mua_ctxpred_all ...
            ...
            mua_taskpred_k_all ...
            mua_taskpred_all ...
            mua_taskpred_reduced_all ...
            mua_taskpred_expl_var_total_all ...
            mua_taskpred_expl_var_partial_all ...
            ...
            mua_ctxpred_taskpred_k_all ...
            mua_ctxpred_taskpred_all ...
            mua_ctxpred_taskpred_reduced_all ...
            mua_ctxpred_taskpred_expl_var_total_all ...
            mua_ctxpred_taskpred_expl_var_partial_all ...
            ...
            fluor_taskpred_k_all ...
            fluor_taskpred_all ...
            fluor_taskpred_reduced_all ...
            fluor_taskpred_expl_var_total_all ...
            fluor_taskpred_expl_var_partial_all ...
            ...
            wheel_ctxpred_all ...
            ctx_str_k_all ...
            ctx_wheel_k_all ...
            wheel_all ...
            movement_all ...
            outcome_all ...
            D_all ...
            ...
            task_regressor_labels ...
            task_regressor_sample_shifts
        
    end
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_taskpred_expl_var_total_all ...
    mua_taskpred_expl_var_partial_all ...
    ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    mua_ctxpred_taskpred_expl_var_total_all ...
    mua_ctxpred_taskpred_expl_var_partial_all ...
    ...
    fluor_taskpred_k_all ...
    fluor_taskpred_all ...
    fluor_taskpred_reduced_all ...
    fluor_taskpred_expl_var_total_all ...
    fluor_taskpred_expl_var_partial_all ...
    ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    movement_all ...
    outcome_all ...
    D_all ...
    ...
    task_regressor_labels ...
    task_regressor_sample_shifts

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_nodeconv'];
save([save_path filesep save_fn],'-v7.3');

