% Batch scripts to save preprocessed data here, saved to: 
% C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data
%
% Things lower down can be dependent on things higher up

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
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        % Wheel
        event_aligned_wheel_raw = interp1(Timeline.rawDAQTimestamps, ...
            wheel_position,t_peri_event);
        event_aligned_wheel = bsxfun(@minus,event_aligned_wheel_raw, ...
            nanmedian(event_aligned_wheel_raw(:,t < 0),2));
        
        % Reward
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];      
        event_aligned_reward = (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_bins(x,:)), ...
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
    
    clearvars -except n_aligned_depths animals curr_animal fluor_all mua_all wheel_all reward_all D_all
    
end
clearvars -except n_aligned_depths fluor_all mua_all wheel_all reward_all D_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld'];
save([save_path filesep save_fn],'-v7.3');

%% Pull out and save activity in each trial (passive fullscreen - trained)
clear all
disp('Getting trial activity passive fullscreen trained')

n_aligned_depths = 4;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

fluor_all = cell(length(animals),1);
mua_all = cell(length(animals),1);
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
        %         depth_group = discretize(spikeDepths,depth_group_edges);
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,3];
        upsample_factor = 3;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);      
        
        use_align = stimOn_times;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        % Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        % MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            % (for all spikes in depth group)
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        % Wheel
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
        wheel_all{curr_animal}{curr_day} = event_aligned_wheel;
        D_all{curr_animal}{curr_day} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
    
    clearvars -except n_aligned_depths animals curr_animal fluor_all mua_all wheel_all D_all
    
end
clearvars -except n_aligned_depths fluor_all mua_all wheel_all D_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_passive_fullscreen'];
save([save_path filesep save_fn],'-v7.3');

%% Pull out and save activity in each trial (passive choiceworld - trained)
clear all
disp('Getting trial activity passive choiceworld trained')

n_aligned_depths = 4;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

fluor_all = cell(length(animals),1);
mua_all = cell(length(animals),1);
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
        %         depth_group = discretize(spikeDepths,depth_group_edges);
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,3];
        upsample_factor = 3;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);      
        
        use_align = stimOn_times;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        % Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        % MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            % (for all spikes in depth group)
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        % Wheel
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
        wheel_all{curr_animal}{curr_day} = event_aligned_wheel;
        D_all{curr_animal}{curr_day} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
    
    clearvars -except n_aligned_depths animals curr_animal fluor_all mua_all wheel_all D_all
    
end
clearvars -except n_aligned_depths fluor_all mua_all wheel_all D_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_passive_choiceworld'];
save([save_path filesep save_fn],'-v7.3');

%% Pull out and save activity in each trial (passive fullscreen - naive)
clear all
disp('Getting trial activity passive fullscreen naive')

n_aligned_depths = 4;

animals = {'AP032','AP033','AP034','AP035','AP036'};

fluor_all = cell(length(animals),1);
mua_all = cell(length(animals),1);
wheel_all = cell(length(animals),1);
reward_all = cell(length(animals),1);
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
        %         depth_group = discretize(spikeDepths,depth_group_edges);
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,3];
        upsample_factor = 3;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);      
        
        use_align = stimOn_times;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        % Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        % MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            % (for all spikes in depth group)
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        % Wheel
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
        wheel_all{curr_animal}{curr_day} = event_aligned_wheel;
        D_all{curr_animal}{curr_day} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
    
    clearvars -except n_aligned_depths animals curr_animal fluor_all mua_all wheel_all D_all
    
end
clearvars -except n_aligned_depths fluor_all mua_all wheel_all D_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_passive_fullscreen_naive'];
save([save_path filesep save_fn],'-v7.3');

%% Pull out and save activity in each trial (passive choiceworld - naive)
clear all
disp('Getting trial activity passive choiceworld naive')

n_aligned_depths = 4;

animals = {'AP032','AP033','AP034','AP035','AP036'};

fluor_all = cell(length(animals),1);
mua_all = cell(length(animals),1);
wheel_all = cell(length(animals),1); 
reward_all = cell(length(animals),1);  
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
        %         depth_group = discretize(spikeDepths,depth_group_edges);
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,3];
        upsample_factor = 3;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);      
        
        use_align = stimOn_times;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        % Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        % MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            % (for all spikes in depth group)
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        % Wheel
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
        wheel_all{curr_animal}{curr_day} = event_aligned_wheel;
        D_all{curr_animal}{curr_day} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
    
    clearvars -except n_aligned_depths animals curr_animal fluor_all mua_all wheel_all D_all
    
end
clearvars -except n_aligned_depths fluor_all mua_all wheel_all D_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_passive_choiceworld_naive'];
save([save_path filesep save_fn],'-v7.3');

%% Regress task events to cortex/striatum activity

%%%% TO DO HERE: low-pass filter?

% Load trial choiceworld data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
data_fn = ['trial_activity_choiceworld'];

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);
reward_all = cellfun(@(data,use_expts) data(use_expts),reward_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);
reward_all = reward_all(use_animals);

% Get widefield ROIs
n_rois = 200;

% MUA depths
n_depths = size(mua_all{1}{1},3);

% Set up conditions
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_rois,1);
mua_allcat = nan(0,length(t),n_depths,1);
wheel_allcat = nan(0,length(t),1);
reward_allcat = nan(0,length(t),1);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

D_cat = struct('stimulus',cell(1,1),'response',cell(1,1),'repeatNum',cell(1,1));

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    n_align = size(fluor_cat,4);
    
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x,[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
%     wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
    % Concatenate behavioural data
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
    D.response = cell2mat(cellfun(@(x) x.response,D_all{curr_animal},'uni',false));
    D.repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all{curr_animal},'uni',false));
    
    D.day = trial_day;
    
    D_cat.stimulus = vertcat(D_cat.stimulus,D.stimulus);
    D_cat.response = vertcat(D_cat.response,D.response);
    D_cat.repeatNum = vertcat(D_cat.repeatNum,D.repeatNum);
    
    % Get trial ID   
    trial_contrast = max(D.stimulus,[],2);
    [~,side_idx] = max(D.stimulus > 0,[],2);
    trial_side = (side_idx-1.5)*2;
    trial_choice = -(D.response-1.5)*2;
    
    trial_conditions = ...
        [trial_contrast, trial_side, ...
        trial_choice, ones(size(trial_day))];
    [~,trial_id] = ismember(trial_conditions,conditions,'rows');
    
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    reward_allcat = [reward_allcat;vertcat(reward_all{curr_animal}{:})];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];
   
end

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Get maximum wheel velocity in chosen direction
wheel_velocity_allcat = [zeros(size(wheel_allcat,1),1,1),diff(wheel_allcat,[],2)];
[max_speed,max_vel_idx] = max(abs(wheel_velocity_allcat(:,:,1).* ...
    (bsxfun(@times,wheel_velocity_allcat(:,:,1),trial_choice_allcat) > 0)),[],2);
vel_t = t(max_vel_idx);

max_vel = max_speed.*trial_choice_allcat;

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Downsample (otherwise it's too much for regression) and d(smooth(fluor))
downsample_factor = 4;
t_downsample = linspace(t(1),t(end),round(length(t)/downsample_factor));
t_diff =  conv(t,[1,1]/2,'valid');
t_downsample_diff = conv(t_downsample,[1,1]/2,'valid');

smooth_factor = downsample_factor;
fluor_allcat_downsamp_diff = permute(interp1(t_diff,permute(diff(convn( ...
    fluor_allcat,ones(1,smooth_factor)/smooth_factor,'same'),[],2),[2,1,3,4]),t_downsample_diff),[2,1,3,4]);

mua_allcat_downsamp = permute(interp1(t,permute(mua_allcat,[2,1,3,4]),t_downsample_diff),[2,1,3,4]);

wheel = interp1(t,wheel_allcat(:,:,1)',t_downsample_diff)';

%%% Regression (separate regressors to get partial explained variance)

% Stim regressors
contrasts = unique(trial_contrast_allcat(trial_contrast_allcat > 0));
contrastsides = sort([-contrasts;contrasts]);

stim_regressors = zeros(size(wheel,1),size(wheel,2),length(contrastsides));
for curr_condition_idx = 1:length(contrastsides)
    stim_regressors(trial_contrast_allcat.*trial_side_allcat == ...
        contrastsides(curr_condition_idx),find(t_downsample_diff > 0,1),curr_condition_idx) = 1;
end

% Move onset regressors (L/R)
[~,move_idx] = max(abs(wheel(:,:,1)) > 2,[],2);
move_onset_regressors = zeros(size(wheel,1),size(wheel,2),2);
for curr_trial = 1:size(move_onset_regressors,1)
    
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

% Move ongoing regressor - 1 if ongoing movement in that direction
wheel_vel = [zeros(size(wheel,1),1),diff(wheel(:,:,1),[],2)];
move_ongoing_regressors = zeros(size(wheel,1),size(wheel,2),2);
move_ongoing_regressors(:,:,1) = wheel_vel < -2;
move_ongoing_regressors(:,:,2) = wheel_vel > 2;

% Wheel regressors - separate leftward/rightward velocity
wheel_vel_norm = [zeros(size(wheel,1),1),diff(wheel(:,:,1),[],2)];
wheel_vel_norm = wheel_vel_norm/prctile(abs(wheel_vel_norm(:)),95);
wheel_vel_norm(abs(wheel_vel_norm) > 1) = sign(wheel_vel_norm(abs(wheel_vel_norm) > 1));
wheel_regressors = abs(repmat(wheel_vel_norm,1,1,2).*cat(3,wheel_vel_norm < 0,wheel_vel_norm > 0)); 

% Go cue regressors - separate for early/late move
go_cue_regressors = zeros(size(wheel,1),size(wheel,2));
go_cue_regressors(move_t <= 0.5,find(t_downsample_diff > 0.5,1),1) = 1;
go_cue_regressors(move_t > 0.5,find(t_downsample_diff > 0.5,1),2) = 1;

% Reward regressors
reward_allcat_downsamp = permute(interp1(t,permute(reward_allcat,[2,1,3,4]),t_downsample_diff,'nearest'),[2,1,3,4]);

reward_allcat_regressor = zeros(size(reward_allcat,1),length(t_downsample_diff));
for curr_trial = 1:size(reward_allcat,1)
   curr_reward = find(reward_allcat(curr_trial,:));
   for i = curr_reward
       curr_reward_t = t(i);
       reward_allcat_regressor(curr_trial,find(t_downsample_diff > curr_reward_t,1)) = 1;
   end
end

regressors = {stim_regressors;move_onset_regressors;move_ongoing_regressors;go_cue_regressors;reward_allcat_regressor};
regressor_labels = {'Stim','Move onset','Move ongoing','Go cue','Reward'};

% Set regression parameters
regress_trials = true(size(move_t));
t_shifts = {[-0.1,0.6];[-0.2,0.5];[0,0];[-0.05,0.3];[-0.1,0.5]};

sample_shifts = cellfun(@(x) round(x(1)*(sample_rate/downsample_factor)): ...
    round(x(2)*(sample_rate/downsample_factor)),t_shifts,'uni',false);
lambda = 0;
zs = [false,false];
cvfold = 5;
return_constant = true;

%%% Do regression on fluor
disp('Regressing task to cortex');
fluor_allcat_predicted = nan(size(fluor_allcat_downsamp_diff));

fluor_kernel = cell(length(regressors)+1,n_rois);
fluor_expl_var = cell(length(regressors),n_rois);
for curr_roi = 1:n_rois
    
    activity = fluor_allcat_downsamp_diff(regress_trials,:,curr_roi,1);
    
    activity_reshape = reshape(activity',[],1)';
    regressors_reshape = cellfun(@(x) ...
        reshape(permute(x(regress_trials,:,:),[2,1,3]),[],size(x,3))',regressors,'uni',false);
    
    [kernel,activity_predicted_reshape,expl_var,activity_predicted_reduced_reshape] = AP_regresskernel( ...
        cellfun(@(x) x(:,~isnan(activity_reshape)),regressors_reshape,'uni',false), ...
        activity_reshape(~isnan(activity_reshape)),sample_shifts,lambda,zs,cvfold,return_constant);
    
    activity_predicted = nan(size(activity))';
    % (to use full model)
    activity_predicted(~isnan(activity')) = activity_predicted_reshape;
    % (to use reduced model)
%     activity_predicted(~isnan(activity')) = activity_predicted_reduced_reshape(:,:,1);
    activity_predicted = activity_predicted';    
    fluor_allcat_predicted(regress_trials,:,curr_roi,1) = activity_predicted;    
    
    fluor_kernel(:,curr_roi) = kernel;
    
    AP_print_progress_fraction(curr_roi,n_rois);
    
end

%%% Do regression on MUA
disp('Regressing task to striatum');
mua_allcat_predicted = nan(size(mua_allcat_downsamp));

mua_kernel = cell(length(regressors)+1,n_depths);
mua_expl_var = cell(length(regressors),n_depths);
for curr_depth = 1:n_depths
    
    activity = mua_allcat_downsamp(regress_trials,:,curr_depth,1);
    
    activity_reshape = reshape(activity',[],1)';
    regressors_reshape = cellfun(@(x) ...
        reshape(permute(x(regress_trials,:,:),[2,1,3]),[],size(x,3))',regressors,'uni',false);
    
    [kernel,activity_predicted_reshape,expl_var,activity_predicted_reduced_reshape] = AP_regresskernel( ...
        cellfun(@(x) x(:,~isnan(activity_reshape)),regressors_reshape,'uni',false), ...
        activity_reshape(~isnan(activity_reshape)),sample_shifts,lambda,zs,cvfold,return_constant);
    
    activity_predicted = nan(size(activity))';
    % (to use full model)
    activity_predicted(~isnan(activity')) = activity_predicted_reshape;
    % (to use reduced model)
%     activity_predicted(~isnan(activity')) = activity_predicted_reduced_reshape(:,:,1);
    activity_predicted = activity_predicted';    
    mua_allcat_predicted(regress_trials,:,curr_depth,1) = activity_predicted;
    
    mua_kernel(:,curr_depth) = kernel;
    
    AP_print_progress_fraction(curr_depth,n_depths);
    
end

% Save regression results
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_task_regression'];
save([save_path filesep save_fn], ...
    'trial_contrast_allcat','trial_side_allcat','trial_choice_allcat', ...
    'downsample_factor','t_downsample_diff','wheel', ...
    'fluor_allcat_downsamp_diff','mua_allcat_downsamp', ...
    'regressors','regressor_labels','t_shifts','regress_trials', ...
    'fluor_allcat_predicted','fluor_kernel', ...
    'mua_allcat_predicted','mua_kernel','-v7.3');
disp('Saved regression');

%% Make reference image for drawing widefield ROIs


% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
data_fn = ['trial_activity_choiceworld'];

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);
reward_all = cellfun(@(data,use_expts) data(use_expts),reward_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);
reward_all = reward_all(use_animals);

% Get number of widefield components and MUA depths
n_vs = size(fluor_all{1}{1},3);
n_depths = size(mua_all{1}{1},3);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_vs,1);
mua_allcat = nan(0,length(t),n_depths,1);
wheel_allcat = nan(0,length(t),1);
reward_allcat = nan(0,length(t),1);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

D_cat = struct('stimulus',cell(1,1),'response',cell(1,1),'repeatNum',cell(1,1));

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
        
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths),'uni',false));
    
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x,[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
%     wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
    % Concatenate behavioural data
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
    D.response = cell2mat(cellfun(@(x) x.response,D_all{curr_animal},'uni',false));
    D.repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all{curr_animal},'uni',false));
    
    D.day = trial_day;
    
    D_cat.stimulus = vertcat(D_cat.stimulus,D.stimulus);
    D_cat.response = vertcat(D_cat.response,D.response);
    D_cat.repeatNum = vertcat(D_cat.repeatNum,D.repeatNum);
    
    % Get trial ID   
    trial_contrast = max(D.stimulus,[],2);
    [~,side_idx] = max(D.stimulus > 0,[],2);
    trial_side = (side_idx-1.5)*2;
    trial_choice = -(D.response-1.5)*2;
    
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    reward_allcat = [reward_allcat;vertcat(reward_all{curr_animal}{:})];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];
   
end

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

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

t_diff = conv(t,[1,1]/2,'valid');

use_t = t_diff >= 0 & t_diff <= 0.8;
plot_px_use_t = double(plot_px(:,:,use_t));
plot_px_use_t(plot_px_use_t < 0) = 0;
plot_px_use_t = bsxfun(@rdivide,plot_px_use_t,max(max(plot_px_use_t,[],1),[],2));

plot_px_com = sum((plot_px_use_t.*permute(1:sum(use_t),[1,3,2])),3)./sum(plot_px_use_t,3);

t_leeway = 20;
plot_px_colored = ...
    ind2rgb(round(mat2gray(plot_px_com,[t_leeway,sum(use_t)-t_leeway])*255),jet(255));

plot_px_alpha = mat2gray(max(plot_px_use_t,[],3), ...
    [0,prctile(reshape(max(plot_px_use_t,[],3),[],1),99)]);
figure;

h = image(plot_px_colored);
set(h,'AlphaData',plot_px_alpha);

axis image off
AP_reference_outline('ccf_aligned','k');

wf_roi_ref = plot_px_colored;

wf_roi_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois';
wf_ref_fn = 'wf_roi_ref';
save([wf_roi_path filesep wf_ref_fn],'wf_roi_ref');

%% Draw widefield ROIs


alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment';
load([alignment_path filesep 'animal_wf_tform']);
im_size = animal_wf_tform(1).im_size;

roi_areas = {'V1','AM','RSPp','PPC','FRm','FRl','SMl','SMf'};

wf_roi = struct('area',cell(length(roi_areas),2),'mask',cell(length(roi_areas),2));

% Draw all ROIs on left hemisphere
for curr_area = 1:length(roi_areas)
    curr_roi_name = [roi_areas{curr_area} '_L'];
    disp(curr_roi_name);
    wf_roi(curr_area).area = curr_roi_name;
    [~,wf_roi(curr_area).mask] = AP_svd_roi(nan(im_size),[],'master');
end

% Reflect all ROIs to right hemisphere
for curr_area = 1:length(roi_areas)
    curr_roi_name = [roi_areas{curr_area} '_R'];
    wf_roi(curr_area,2).area = curr_roi_name;
    
    L_roi = wf_roi(curr_area,1).mask;
    R_roi = AP_reflect_widefield(L_roi) > 0;
    
    wf_roi(curr_area,2).mask = R_roi;
end

wf_roi_plot = sum(cat(3,wf_roi(:).mask),3);
figure;imagesc(wf_roi_plot);
colormap(flipud(gray));
AP_reference_outline('ccf_aligned','r');AP_reference_outline('retinotopy','b');
axis image off;
title('New ROIs')

wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
save(wf_roi_fn,'wf_roi');
disp('Saved new widefield ROIs');










