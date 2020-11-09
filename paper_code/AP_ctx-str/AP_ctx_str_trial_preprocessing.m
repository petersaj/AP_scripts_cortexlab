%% Notes
%
% AP_ctx_str_ephys_alignment - beforehand: gets striatal domains to group MUA
%
% Batch scripts to save preprocessed data here, saved to:
% C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data


%% ~~~~~~~~~~~~~ PRIMARY RECORDINGS ~~~~~~~~~~~~~

%% Get widefield correlation boundaries
clear all
animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';

wf_corr_borders = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    disp(animal);
    
    experiments = AP_find_experiments(animal,protocol);
    % (wf-only days: no craniotomy)
    experiments = experiments([experiments.imaging] & ~[experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = false;
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        fVdf_deconv = AP_deconv_wf(fVdf);
        aUdf = AP_align_widefield(Udf,animal,day);       

        % Get V covariance
        Ur = reshape(aUdf, size(aUdf,1)*size(aUdf,2),[]); % P x S
        covV = cov(fVdf_deconv'); % S x S % this is the only one that takes some time really
        varP = dot((Ur*covV)', Ur'); % 1 x P
        
        ySize = size(aUdf,1); xSize = size(aUdf,2);

        px_spacing = 20;
        use_y = 1:px_spacing:size(aUdf,1);
        use_x = 1:px_spacing:size(aUdf,2);
        corr_map = cell(length(use_y),length(use_x));
        for curr_x_idx = 1:length(use_x)
            curr_x = use_x(curr_x_idx);
            for curr_y_idx = 1:length(use_y)
                curr_y = use_y(curr_y_idx);
                
                pixel = [curr_y,curr_x];
                pixelInd = sub2ind([ySize, xSize], pixel(1), pixel(2));
                
                covP = Ur(pixelInd,:)*covV*Ur'; % 1 x P
                stdPxPy = varP(pixelInd).^0.5 * varP.^0.5; % 1 x P
                corrMat = reshape(covP./stdPxPy,ySize,xSize); % 1 x P
                
                corr_map{curr_y_idx,curr_x_idx} = corrMat;
            end
        end       
        
        % Get downsampled correlation maps
        downsample_factor = 10;
        corr_map_downsamp = cellfun(@(x) ...
            imresize(x,1/downsample_factor,'bilinear'),corr_map,'uni',false);
                
%         % OLD: Correlation map edge detection
%         corr_map_cat = cat(3,corr_map{:});
%         corr_map_edge = corr_map_cat-imgaussfilt(corr_map_cat,20);
%         corr_map_edge_norm = reshape(zscore(reshape(corr_map_edge,[], ...
%             size(corr_map_edge,3)),[],1),size(corr_map_edge));      
%         corr_edges = nanmean(corr_map_edge,3);
        
        % NEW: Correlation map edge detection
        corr_map_cat = cat(3,corr_map{:});
        corr_map_edge = corr_map_cat-imfilter(corr_map_cat,fspecial('disk',20));
        corr_map_edge_norm = reshape(zscore(reshape(corr_map_edge,[], ...
            size(corr_map_edge,3)),[],1),size(corr_map_edge));      
        corr_edges = nanmean(corr_map_edge,3);
        % (get rid of high frequency of edges - vasculature/bone)
        corr_edges_highpass = corr_edges - imfilter(corr_edges,fspecial('disk',10));
        corr_edges_lowpass = corr_edges - corr_edges_highpass;
        
        wf_corr_borders(curr_animal).corr_map_downsamp{curr_day} = corr_map_downsamp;
        wf_corr_borders(curr_animal).corr_edges{curr_day} = corr_edges_lowpass;
        
        clearvars -except animals animal curr_animal protocol experiments curr_day animal wf_corr_borders load_parts
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
end

% Save widefield correlation borders
save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_borders'];
save_fn = [save_path filesep 'wf_corr_borders'];
save(save_fn,'wf_corr_borders');
disp(['Saved ' save_fn]);


%% Striatum cortical kernels (depth-aligned)

clear all

disp('Cortex -> striatum regression maps across protocols');

n_aligned_depths = 16;

% Parameters for regression
regression_params.use_svs = 1:100;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.1,0.1];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

protocols = {'vanillaChoiceworld'};

for protocol = protocols
    protocol = cell2mat(protocol);
    
    animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
    
    init_array = cellfun(@(x) cell(0,0),animals','uni',false);
    ctx_str_kernel = init_array;
    ctx_str_expl_var = init_array;
    ctx_str_sta = init_array;
    
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
            str_align = 'depth';
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
            fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
            
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
                AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
                binned_spikes_std,kernel_frames,lambda, ...
                regression_params.zs,regression_params.cvfold, ...
                false,regression_params.use_constant);
            
            % Reshape kernel and convert to pixel space
            aUdf = AP_align_widefield(Udf,animal,day);
            k_px = zeros(size(aUdf,1),size(aUdf,2),size(k,2),size(k,3),'single');
            for curr_spikes = 1:size(k,3)
                k_px(:,:,:,curr_spikes) = svdFrameReconstruct(aUdf(:,:,regression_params.use_svs),k(:,:,curr_spikes));
            end
            
            %%% Get spike-triggered average
            sta_v = fVdf_deconv_resample*binned_spikes_std';
            sta_px = svdFrameReconstruct(aUdf,sta_v);
            
            % Store
            ctx_str_kernel{curr_animal}{curr_day} = k_px;
            ctx_str_expl_var{curr_animal}{curr_day} = explained_var.total;
            ctx_str_sta{curr_animal}{curr_day} = sta_px;
            
            AP_print_progress_fraction(curr_day,length(experiments));
            clearvars -except regression_params n_aligned_depths ...
                animals animal curr_animal protocol ...
                experiments curr_day animal load_parts ...
                ctx_str_kernel ctx_str_expl_var ctx_str_sta
            
        end
        
        disp(['Finished ' animal]);
        
    end
    
    % Save
    save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data'];
    save_fn = [save_path filesep 'ctx_str_kernels_' protocol '_' num2str(n_aligned_depths) 'strdepths'];
    save(save_fn,'-v7.3');
    disp(['Saved ' save_fn]);
    
end


%% Striatum cortical kernels (kernel-aligned)

clear all

disp('Cortex -> striatum regression maps across protocols');

% Parameters for regression
regression_params.use_svs = 1:100;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.1,0.1];
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


%% Choiceworld trial activity (widefield-only days)

clear all
disp('Choiceworld trial activity (widefield-only days)')

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize save variable
trial_data_all = struct;

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
        
        % Pull out trial data
        AP_ctx_str_grab_trial_data;
        
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
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars -except animals curr_animal animal protocol experiments curr_day ...
            trial_data_all
        
    end
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_filename = [save_path filesep 'trial_activity_choiceworld_wfonly'];
save(save_filename,'-v7.3');
disp(['Saved ' save_filename]);


%% Choiceworld trial activity (striatum depth)

clear all
disp('Choiceworld trial activity (striatum depth)')

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
n_aligned_depths = 16;

% Initialize save variable
trial_data_all = struct;

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
        
        % Pull out trial data
        AP_ctx_str_grab_trial_data;
        
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
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars -except animals n_aligned_depths ...
            curr_animal animal protocol experiments curr_day ...
            trial_data_all
        
    end
end

clearvars -except trial_data_all n_aligned_depths
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_' num2str(n_aligned_depths) 'strdepth'];
save([save_path filesep save_fn],'-v7.3');


%% Choiceworld trial activity (striatum domain)

clear all
disp('Choiceworld trial activity (striatum domain)')

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize save variable
trial_data_all = struct;

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
        
        % Pull out trial data
        AP_ctx_str_grab_trial_data;
        
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
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars -except animals curr_animal animal protocol experiments curr_day ...
            trial_data_all
        
    end
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_filename = [save_path filesep 'trial_activity_choiceworld'];
save(save_filename,'-v7.3');
disp(['Saved ' save_filename]);


%% Passive trial activity (all protocols and animal groups)

clear all

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
        
        % Initialize save variable
        trial_data_all = struct;
        
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
            
            disp(['Loading ' animal]);
            
            for curr_day = 1:length(experiments)
                
                day = experiments(curr_day).day;
                experiment = experiments(curr_day).experiment(end);
                
                % Load experiment
                str_align = 'kernel';
                AP_load_experiment;
                
                % Pull out trial data
                AP_ctx_str_grab_trial_data;
                
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
            end
            
            % Clear for next loop
            clearvars -except ...
                trained_animals naive_animals protocols ...
                protocol animal_group animals curr_animal ...
                animals curr_animal animal protocol experiments curr_day ...
                trial_data_all
            
        end
        
        clearvars -except ...
            trained_animals naive_animals protocols ...
            protocol animal_group animals...
            trial_data_all
        disp('Finished loading all')
        
        save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
        save_fn = ['trial_activity_' protocol '_' animal_group];
        save([save_path filesep save_fn],'-v7.3');
        
    end
end



%% ~~~~~~~~~~~~~ Widefield ROIs ~~~~~~~~~~~~~

%% [[lock on]]
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
    n_aligned_depths = 3;
    kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template_' num2str(n_aligned_depths) '_depths.mat'];
    load(kernel_template_fn);
    
    kernel_bw = false(size(kernel_template));
    
    % Set cutoff by fraction of max weight, keep large blobs ipsi
    bregma = allenCCFbregma;
    ccf_tform_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\ccf_tform'];
    load(ccf_tform_fn);
    
    um2pixel = 20.6;
    bregma_resize = bregma*(10/um2pixel);
    bregma_align = [bregma_resize([3,1]),1]*ccf_tform.T;
    ipsi_side = true(size(kernel_template(:,:,1)));
    ipsi_side(:,round(bregma_align(1)):end,:) = false;
    
    frac_max_weight = 0.75; % zero pixels > max weight * this
    min_px = 100; % get rid of small islands
    for curr_depth = 1:n_aligned_depths
        kernel_bw(:,:,curr_depth) = ...
            bwareaopen(kernel_template(:,:,curr_depth) > ...
            max(reshape(kernel_template(:,:,curr_depth), ...
            [],1),[],1)*frac_max_weight,min_px).*ipsi_side;
    end
    
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
    
    %     % Get the average full kernel
    %     n_aligned_depths = 4;
    %     data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_ephys';
    %     k_fn = [data_path filesep 'wf_ephys_maps_vanillaChoiceworld_' num2str(n_aligned_depths) '_depths'];
    %     load(k_fn);
    %     k_px_trained_cat = cellfun(@(x) x(:,:,end:-1:1,:),[batch_vars(1:6).r_px],'uni',false);
    %     k_px_trained = nanmean(double(cat(5,k_px_trained_cat{:})),5);
    
    % Save kernel ROIs
    kernel_roi = struct;
    kernel_roi.bw = kernel_bw;
    kernel_roi.max_weighted = kernel_max_weighted;
    %     kernel_roi.kernel = k_px_trained;
    kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
    save(kernel_roi_fn,'kernel_roi');
    disp('Saved kernel ROIs');
    
    %% [[lock off]]
end

%% ~~~~~~~~~~~ MUSCIMOL ~~~~~~~~~~~~~

%% Widefield VFS and std dev (pre/post muscimol)

clear all

animals = {'AP045','AP054','AP055','AP053','AP047','AP048'};

% Load master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');

muscimol_wf = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    disp(animal);
    
    % Get standard deviation from task and retinotopy from sparse noise
    task_experiments = AP_find_experiments(animal,'vanillaChoiceworldNoRepeats');
    task_experiments = task_experiments([task_experiments.imaging] & [task_experiments.ephys]);
    
    retinotopy_experiments = AP_find_experiments(animal,'AP_sparseNoise');
    retinotopy_experiments = retinotopy_experiments([retinotopy_experiments.imaging] & [retinotopy_experiments.ephys]);
    
    for curr_day = 1:length(task_experiments)
        
        day = task_experiments(curr_day).day;
        
        if length(task_experiments(curr_day).experiment) ~= 2
            warning([animal ' ' day '~= 2 experiments, using first/last']);
            task_experiments(curr_day).experiment = task_experiments(curr_day).experiment([1,end]);
        end
        
        if length(retinotopy_experiments(curr_day).experiment) ~= 2
            warning([animal ' ' day '~= 2 experiments, using first/last']);
            retinotopy_experiments(curr_day).experiment = retinotopy_experiments(curr_day).experiment([1,end]);
        end
        
        for curr_exp = 1:2
            
            % Get retinotopy (from sparse noise)
            experiment = retinotopy_experiments(curr_day).experiment(curr_exp);
            load_parts.imaging = true;
            AP_load_experiment;
            lilrig_retinotopy;
            close(gcf);
            vfs_aligned = AP_align_widefield(vfs_median,animal,day);
            
            % Get standard deviation (from task)
            % (do this in master-U space)
            experiment = task_experiments(curr_day).experiment(curr_exp);
            load_parts.imaging = true;
            AP_load_experiment;
            
            px_std_sq = zeros(size(Udf,1),size(Udf,2));
            px_mean = svdFrameReconstruct(Udf,nanmean(fVdf,2));
            
            skip_frames = 35*10;
            n_frames = size(fVdf,2) - skip_frames*2 + 1;
            chunk_size = 5000;
            frame_chunks = unique([skip_frames:chunk_size:size(fVdf,2)-skip_frames, ...
                size(fVdf,2)-skip_frames]);
            
            for curr_chunk = 1:length(frame_chunks)-1
                curr_im = svdFrameReconstruct(Udf,fVdf(:,frame_chunks(curr_chunk): ...
                    frame_chunks(curr_chunk+1)));
                px_std_sq = px_std_sq + sum((bsxfun(@minus,curr_im,px_mean).^2./n_frames),3);
                clear curr_im
            end
            px_std = sqrt(px_std_sq);
            px_std_aligned = AP_align_widefield(px_std,animal,day);
            avg_im_aligned = AP_align_widefield(avg_im,animal,day);
            
            % Package final
            muscimol_wf(curr_animal).vfs{curr_day,curr_exp} = vfs_aligned;
            muscimol_wf(curr_animal).std{curr_day,curr_exp} = px_std_aligned;
            muscimol_wf(curr_animal).avg{curr_day,curr_exp} = avg_im_aligned;
            
        end
        
        AP_print_progress_fraction(curr_day,length(task_experiments));
        
    end
end

% Save
save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data'];
save_fn = 'muscimol_wf.mat';
save([save_path filesep save_fn],'muscimol_wf');
disp(['Saved ' [save_path filesep save_fn]]);


%% Striatum cortical kernels (pre/post muscimol)

clear all

disp('Cortex -> striatum regression maps across protocols');

n_aligned_depths = 3;

% Parameters for regression
regression_params.use_svs = 1:100;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.1,0.1];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

exp_conditions = {'pre_muscimol','post_muscimol'};

for curr_exp_condition = 1:2
    
    protocols = {'vanillaChoiceworldNoRepeats','AP_sparseNoise'};
    
    for protocol = protocols
        protocol = cell2mat(protocol);
        
        animals = {'AP045','AP054','AP055','AP053','AP047','AP048'};
        
        init_array = cellfun(@(x) cell(0,0),animals','uni',false);
        ctx_str_kernel = init_array;
        ctx_str_expl_var = init_array;
        
        for curr_animal = 1:length(animals)
            
            animal = animals{curr_animal};
            disp(['Loading ' animal ' ' protocol]);
            
            experiments = AP_find_experiments(animal,protocol);
            experiments = experiments([experiments.imaging] & [experiments.ephys]);
            
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
                
                % Assumption: if more than 2 experiments (pre/post), use first
                % last (because the failures always hapened post)
                n_experiments = length(experiments(curr_day).experiment);
                if n_experiments ~= 2
                    warning([animal ' ' day ': ' num2str(n_experiments) ' experiments, using first/last'])
                end
                condition_experiments = experiments(curr_day).experiment([1,end]);
                experiment = condition_experiments(curr_exp_condition);
                
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
                    ctx_str_kernel ctx_str_expl_var ...
                    exp_conditions curr_exp_condition ...
                    
            end
            
            disp(['Finished ' animal]);
            
        end
        
        % Save
        save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data'];
        save_fn = [save_path filesep 'ctx_str_kernels_' protocol '_' exp_conditions{curr_exp_condition}];
        save(save_fn,'-v7.3');
        warning('saving -v7.3');
        disp(['Saved ' save_fn]);
        
    end
    
end

%% Choiceworld trial activity (pre/post muscimol)

clear all

animals = {'AP045','AP054','AP055','AP053','AP047','AP048'};
protocol = 'vanillaChoiceworldNoRepeats';
exp_conditions = {'pre_muscimol','post_muscimol'};

for curr_exp_condition = 1:2
    
    % Initialize save variable
    trial_data_all = struct;
    
    for curr_animal = 1:length(animals)
        
        animal = animals{curr_animal};
        
        experiments = AP_find_experiments(animal,protocol);
        experiments = experiments([experiments.imaging] & [experiments.ephys]);
        
        disp(['Loading ' animal]);
        
        for curr_day = 1:length(experiments)
            
            day = experiments(curr_day).day;
            % Assumption: if more than 2 experiments (pre/post), use first
            % last (because the failures always hapened post)
            n_experiments = length(experiments(curr_day).experiment);
            if n_experiments ~= 2
                warning([animal ' ' day ': ' num2str(n_experiments) ' experiments, using first/last'])
            end
            condition_experiments = experiments(curr_day).experiment([1,end]);
            experiment = condition_experiments(curr_exp_condition);
            
            % Load experiment
            str_align = 'kernel';
            AP_load_experiment;
            
            % Pull out trial data
            AP_ctx_str_grab_trial_data;
            
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
            
            AP_print_progress_fraction(curr_day,length(experiments));
            
            % Clear for next loop
            clearvars -except ...
                exp_conditions curr_exp_condition ...
                animals curr_animal animal protocol experiments curr_day ...
                trial_data_all
            
        end
    end
    
    clearvars -except ...
        exp_conditions curr_exp_condition animals protocol...
        trial_data_all
    disp('Finished loading all')
    
    save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
    save_fn = ['trial_activity_' protocol '_' exp_conditions{curr_exp_condition}];
    save([save_path filesep save_fn],'-v7.3');
    
end



%% Passive trial activity (pre/post muscimol)

clear all

animals = {'AP045','AP054','AP055','AP053','AP047','AP048'};
protocol = 'AP_lcrGratingPassive';
exp_conditions = {'pre_muscimol','post_muscimol'};

for curr_exp_condition = 1:2
    
    % Initialize save variable
    trial_data_all = struct;
    
    for curr_animal = 1:length(animals)
        
        animal = animals{curr_animal};
        
        experiments = AP_find_experiments(animal,protocol);
        experiments = experiments([experiments.imaging] & [experiments.ephys]);
        
        disp(['Loading ' animal]);
        
        for curr_day = 1:length(experiments)
            
            day = experiments(curr_day).day;
            % Assumption: if more than 2 experiments (pre/post), use first
            % last (because the failures always hapened post)
            n_experiments = length(experiments(curr_day).experiment);
            if n_experiments ~= 2
                warning([animal ' ' day ': ' num2str(n_experiments) ' experiments, using first/last'])
            end
            condition_experiments = experiments(curr_day).experiment([1,end]);
            experiment = condition_experiments(curr_exp_condition);
            
            % Load experiment
            str_align = 'kernel';
            AP_load_experiment;
            
            % Pull out trial data
            AP_ctx_str_grab_trial_data;
            
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
            clearvars -except ...
                exp_conditions curr_exp_condition ...
                animals curr_animal animal protocol experiments curr_day ...
                trial_data_all
            
            AP_print_progress_fraction(curr_day,length(experiments));
        end
    end
    
    clearvars -except ...
        exp_conditions curr_exp_condition animals protocol ...
        trial_data_all
    disp('Finished loading all')
    
    save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
    save_fn = ['trial_activity_' protocol '_' exp_conditions{curr_exp_condition}];
    save([save_path filesep save_fn],'-v7.3');
    disp(['Saved ' save_fn]);
    
end


%% ~~~~~~~~~~~ CORTICAL EPHYS ~~~~~~~~~~~~~

%% Align LFP and get correlation between fluorescence and spikes:

%% [[lock on]]

if false
    
    %% 1) Get stim-aligned LFP
    
    animals = {'AP043','AP060','AP061'};
    
    vis_ctx_ephys = struct;
    
    for curr_animal = 1:length(animals)
        
        animal = animals{curr_animal};
        
        protocol = 'AP_lcrGratingPassive';
        flexible_name = true;
        experiments = AP_find_experiments(animal,protocol);
        experiments = experiments([experiments.imaging] & [experiments.ephys]);
        
        for curr_day = 1:length(experiments)
            
            % Load experiment
            day = experiments(curr_day).day;
            experiment = experiments(curr_day).experiment(end);
            site = 2; % cortex probe is always site 2
            lfp_channel = 1; % initialize the single-channel variables
            str_align = 'none'; % (because this is on cortex)
            load_parts.ephys = true;
            verbose = false;
            AP_load_experiment
            
            % Estimate boundaries of cortex (the dumb way: first template/gap)
            sorted_template_depths = sort(template_depths);
            ctx_start = sorted_template_depths(1) - 1;
            [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
            ctx_end = sorted_template_depths(max_gap_idx)+1;
            ctx_depth = [ctx_start,ctx_end];
            
            % Load in each channel, light-fix, get average stim response
            % (loading script sets up a bunch of this already): load in each channel
            use_stims = trial_conditions(:,2) == 90;
            use_stimOn_times = stimOn_times(use_stims);
            
            stim_surround_window = [-0.05,0.1];
            stim_surround_t = stim_surround_window(1):1/lfp_sample_rate:stim_surround_window(2);
            
            pull_times = use_stimOn_times + stim_surround_t;
            
            lfp_stim_align_avg = nan(n_channels,length(stim_surround_t));
            
            for curr_lfp_channel = 1:n_channels
                
                lfp = double(lfp_memmap{lfp_load_file}.Data.lfp(curr_lfp_channel,lfp_load_start_rel:lfp_load_stop_rel));
                
                % Pull out LFP around light on
                use_blue_on = blue_on >= lfp_t_timeline(1) & blue_on <= lfp_t_timeline(end);
                blue_on_pull_t = blue_on(use_blue_on) + light_surround_t;
                blue_on_lfp = interp1(lfp_t_timeline,lfp',blue_on_pull_t);
                
                use_violet_on = violet_on >= lfp_t_timeline(1) & violet_on <= lfp_t_timeline(end);
                violet_on_pull_t = violet_on(use_violet_on) + light_surround_t;
                violet_on_lfp = interp1(lfp_t_timeline,lfp',violet_on_pull_t);
                
                % Subtract baseline
                baseline_t = find(light_surround_t < 0,1,'last');
                blue_on_lfp_baselinesub = blue_on_lfp - blue_on_lfp(:,baseline_t,:);
                violet_on_lfp_baselinesub = violet_on_lfp - violet_on_lfp(:,baseline_t,:);
                
                % Get rolling median (allow light artifact to change slightly)
                n_light = 500;
                blue_on_lfp_baselinesub_med = movmedian(blue_on_lfp_baselinesub,n_light,1);
                violet_on_lfp_baselinesub_med = movmedian(violet_on_lfp_baselinesub,n_light,1);
                
                % Interpolate out the artifact to remove
                n_lfp_channels = size(lfp,1);
                blue_light_remove = interp1( ...
                    reshape(permute(blue_on_pull_t,[2,1]),[],1), ...
                    reshape(permute(blue_on_lfp_baselinesub_med,[2,1,3]),[],n_lfp_channels), ...
                    reshape(lfp_t_timeline,[],1))';
                violet_light_remove = interp1( ...
                    reshape(permute(violet_on_pull_t,[2,1]),[],1), ...
                    reshape(permute(violet_on_lfp_baselinesub_med,[2,1,3]),[],n_lfp_channels), ...
                    reshape(lfp_t_timeline,[],1))';
                
                % Zero-out any NaNs (e.g. remove nothing)
                blue_light_remove(isnan(blue_light_remove)) = 0;
                violet_light_remove(isnan(violet_light_remove)) = 0;
                
                % Remove the artifact
                lfp_lightfix = lfp - (blue_light_remove + violet_light_remove);
                
                % Additional 300Hz filter: artifact blips sometimes still there
                freqCutoff = 300; % Hz
                [b100s, a100s] = butter(2,freqCutoff/(lfp_sample_rate/2),'low');
                lfp_lightfix_filt = single(filtfilt(b100s,a100s,double(lfp_lightfix)')');
                
                % Get average stim response
                lfp_stim_align = interp1(lfp_t_timeline,lfp_lightfix_filt,pull_times);
                
                baseline_time = stim_surround_t < 0;
                lfp_stim_align_avg(curr_lfp_channel,:) = nanmean((lfp_stim_align - nanmean(lfp_stim_align(:,baseline_time),2)),1);
                
                AP_print_progress_fraction(curr_lfp_channel,n_channels);
                
            end
            
            % Average channels at same depth
            lfp_connected_channels = channel_map_full.connected;
            [lfp_stim_depth,lfp_stim_align_avg_depth] = grpstats( ...
                lfp_stim_align_avg(lfp_connected_channels,:), ...
                lfp_channel_positions(lfp_connected_channels),{'gname','mean'});
            lfp_stim_depth = cellfun(@str2num,lfp_stim_depth);
            
            % Package
            vis_ctx_ephys(curr_animal).animal = animal;
            vis_ctx_ephys(curr_animal).day{curr_day} = day;
            
            vis_ctx_ephys(curr_animal).stim_lfp_t{curr_day} = stim_surround_t;
            
            vis_ctx_ephys(curr_animal).stim_lfp{curr_day} = lfp_stim_align_avg_depth;
            vis_ctx_ephys(curr_animal).stim_lfp_depth{curr_day} = lfp_stim_depth;
            
        end
    end
    
    save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
    save_fn = [save_path filesep 'vis_ctx_ephys.mat'];
    save(save_fn,'vis_ctx_ephys');
    disp(['Saved ' save_fn]);
    
    
    %% 2) Get CSD
    
    % Load stim LFP
    vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
    load(vis_ctx_ephys_fn);
    
    for curr_animal = 1:length(vis_ctx_ephys)
        for curr_day = 1:length(vis_ctx_ephys(curr_animal).stim_lfp)
            
            stim_lfp_t = vis_ctx_ephys(curr_animal).stim_lfp_t{curr_day};
            stim_lfp = vis_ctx_ephys(curr_animal).stim_lfp{curr_day};
            stim_lfp_depth = vis_ctx_ephys(curr_animal).stim_lfp_depth{curr_day};
            
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
            
            % Plot CSD and slice
            figure;
            subplot(1,3,1:2);
            imagesc(stim_lfp_t,stim_csd_depth,stim_csd)
            colormap(brewermap([],'*RdBu'));
            caxis([-max(abs(caxis)),max(abs(caxis))]);
            line(repmat(csd_slice_t(1),2,1),ylim,'color','k');
            line(repmat(csd_slice_t(2),2,1),ylim,'color','k');
            subplot(1,3,3,'YDir','reverse'); hold on;
            plot(csd_slice,stim_csd_depth,'k','linewidth',2)
            linkaxes(get(gcf,'Children'),'y')
            
            % Store
            vis_ctx_ephys(curr_animal).stim_csd{curr_day} = stim_csd;
            vis_ctx_ephys(curr_animal).csd_slice{curr_day} = csd_slice;
            vis_ctx_ephys(curr_animal).stim_csd_depth{curr_day} = stim_csd_depth;
            
        end
    end
    
    % Save csd into struct
    save(vis_ctx_ephys_fn,'vis_ctx_ephys');
    disp(['Saved ' vis_ctx_ephys_fn]);
    
    
    %% 3) Align CSD across recordings
    
    % Load CSD
    vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
    load(vis_ctx_ephys_fn);
    
    % Concatenate for plotting
    animal_cat = cellfun(@(x,y) repmat({x},1,length(y)),{vis_ctx_ephys.animal},{vis_ctx_ephys.day},'uni',false);
    animal_cat = horzcat(animal_cat{:});
    day_cat = [vis_ctx_ephys(:).day];
    n_recordings = length(day_cat);
    
    stim_lfp_t = vis_ctx_ephys(1).stim_lfp_t{1};
    stim_lfp_depth = vis_ctx_ephys(1).stim_lfp_depth{1};
    stim_csd_depth = vis_ctx_ephys(1).stim_csd_depth{1};
    
    stim_lfp_cat = cell2mat(permute([vis_ctx_ephys(:).stim_lfp],[1,3,2]));
    stim_csd_cat = cell2mat(permute([vis_ctx_ephys(:).stim_csd],[1,3,2]));
    csd_slice_cat = cell2mat([vis_ctx_ephys(:).csd_slice]);
    
    % Anchor points: biggest sink,superficial source
    
    % Get median distance between anchor points across recordings
    csd_slice_cat = cell2mat([vis_ctx_ephys(:).csd_slice]);
    
    [~,csd_slice_cat_sink_idx] = min(csd_slice_cat,[],1);
    [~,csd_slice_cat_source_idx_rel] = arrayfun(@(x) ...
        max(csd_slice_cat(csd_slice_cat_sink_idx(x):-1:1,x)), ...
        1:length(csd_slice_cat_sink_idx));
    csd_slice_cat_source_idx = csd_slice_cat_sink_idx - csd_slice_cat_source_idx_rel - 1;
    
    csd_slice_cat_sink_depth = vis_ctx_ephys(1).stim_csd_depth{1}(csd_slice_cat_sink_idx);
    csd_slice_cat_source_depth = vis_ctx_ephys(1).stim_csd_depth{1}(csd_slice_cat_source_idx);
    sink_source_dist = median(csd_slice_cat_sink_depth - csd_slice_cat_source_depth);
    
    % Get aligned depth for each recording, plot, save
    figure;
    curr_recording = 1;
    stim_csd_aligned_scaled_cat = nan(0);
    for curr_animal = 1:length(vis_ctx_ephys)
        for curr_day = 1:length(vis_ctx_ephys(curr_animal).day)
            
            stim_csd = vis_ctx_ephys(curr_animal).stim_csd{curr_day};
            csd_slice = vis_ctx_ephys(curr_animal).csd_slice{curr_day};
            stim_csd_depth = vis_ctx_ephys(curr_animal).stim_csd_depth{curr_day};
            
            % Get anchor points
            [~,csd_slice_sink_idx] = min(csd_slice);
            [~,csd_slice_source_idx_rel] = max(csd_slice(csd_slice_sink_idx:-1:1));
            csd_slice_source_idx = csd_slice_sink_idx - csd_slice_source_idx_rel - 1;
            
            csd_slice_sink_depth = stim_csd_depth(csd_slice_sink_idx);
            csd_slice_source_depth = stim_csd_depth(csd_slice_source_idx);
            
            % Scale relative depth: 0 at first source, constant distance to sink
            stim_csd_depth_aligned = ...
                (stim_csd_depth - csd_slice_source_depth)/ ...
                (csd_slice_sink_depth - csd_slice_source_depth)* ...
                sink_source_dist;
            
            % Plot aligned CSD
            depth_align_interp = [-500:20:2000];
            stim_csd_aligned = interp1(...
                stim_csd_depth_aligned,stim_csd_cat(:,:,curr_recording),depth_align_interp, ...
                'linear','extrap');
            
            animal = vis_ctx_ephys(curr_animal).animal;
            day = vis_ctx_ephys(curr_animal).day{curr_day};
            stim_lfp_t = vis_ctx_ephys(curr_animal).stim_lfp_t{curr_day};
            stim_lfp_depth = vis_ctx_ephys(curr_animal).stim_lfp_depth{curr_day};
            stim_lfp = vis_ctx_ephys(curr_animal).stim_lfp{curr_day};
            
            subplot(3,n_recordings,curr_recording);
            imagesc(stim_lfp_t,stim_lfp_depth,stim_lfp);
            caxis([-max(abs(caxis))/2,max(abs(caxis))/2]);
            colormap(brewermap([],'*RdBu'));
            title({animal,day,'LFP'});
            
            subplot(3,n_recordings,size(stim_csd_cat,3)+curr_recording);
            imagesc(stim_lfp_t,stim_csd_depth,stim_csd);
            caxis([-max(abs(caxis))/2,max(abs(caxis))/2]);
            colormap(brewermap([],'*RdBu'));
            title('CSD');
            
            subplot(3,n_recordings,size(stim_csd_cat,3)*2+curr_recording);
            imagesc(stim_lfp_t,stim_csd_depth_aligned,stim_csd_aligned);
            caxis([-max(abs(caxis))/2,max(abs(caxis))/2]);
            colormap(brewermap([],'*RdBu'));
            title('Aligned CSD');
            
            % Keep aligned CSD for averaging
            stim_csd_aligned_scaled_cat(:,:,curr_recording) = stim_csd_aligned./min(csd_slice);
            
            % Store aligned depth
            vis_ctx_ephys(curr_animal).stim_csd_depth_aligned{curr_day} = stim_csd_depth_aligned;
            
            % Iterate recording index (only used for plotting)
            curr_recording = curr_recording + 1;
            
        end
    end
    
    % Plot average aligned CSD
    figure;
    imagesc(stim_lfp_t,depth_align_interp,nanmean(stim_csd_aligned_scaled_cat,3));
    caxis([-1,1]);
    line([0,0],ylim,'color','k');
    colormap(brewermap([],'*RdBu'));
    xlabel('Time from stim');
    ylabel('Aligned depth');
    title('Average aligned CSD')
    
    
    % Save alignment into struct
    save(vis_ctx_ephys_fn,'vis_ctx_ephys');
    disp(['Saved ' vis_ctx_ephys_fn]);
    
    
    %% 4) Draw cortical widefield ROIs for each cortical recording
    
    % Load stim LFP
    vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
    load(vis_ctx_ephys_fn);
    
    % Draw cortical probe ROI for each recording
    for curr_animal = 1:length(vis_ctx_ephys)
        for curr_day = 1:length(vis_ctx_ephys(curr_animal).day)
            animal = vis_ctx_ephys(curr_animal).animal;
            day = vis_ctx_ephys(curr_animal).day{curr_day};
            [img_path,img_exists] = AP_cortexlab_filename(animal,day,[],'imaging');
            avg_im_blue = readNPY([img_path filesep 'meanImage_blue.npy']);
            
            h = figure;imagesc(avg_im_blue);
            axis image off;
            caxis([0,prctile(avg_im_blue(:),99)]);
            colormap(gray)
            title(['Draw ROI: ' animal ' ' day]);
            curr_mask = roipoly;
            vis_ctx_ephys(curr_animal).ctx_roi{curr_day} = curr_mask;
            close(h);
        end
    end
    
    % Save ROIs into struct
    save(vis_ctx_ephys_fn,'vis_ctx_ephys');
    disp(['Saved ' vis_ctx_ephys_fn]);
    
    %% [[lock off]]
    
end

%% Get activity correlation cortex wf/mua and ctx ephys/str mua (task)

animals = {'AP043','AP060','AP061'};

% Load ephys alignment
vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
load(vis_ctx_ephys_fn);

data = struct;
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    protocol = 'vanillaChoiceworld';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    for curr_day = 1:length(experiments)
        
        % Load experiment
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(1); % (only one doubled experiment, and the second one was worse)
        site = 2; % cortex probe is always site 2
        lfp_channel = 'all';
        str_align = 'none'; % (because this is on cortex)
        verbose = false;
        AP_load_experiment
        
        %%% DEPTH-ALIGN TEMPLATES, FIND CORTEX BOUNDARY
        curr_animal_idx = strcmp(animal,{vis_ctx_ephys.animal});
        curr_day_idx = strcmp(day,vis_ctx_ephys(curr_animal_idx).day);
        curr_csd_depth = vis_ctx_ephys(curr_animal_idx).stim_csd_depth{curr_day_idx};
        curr_csd_depth_aligned = vis_ctx_ephys(curr_animal_idx).stim_csd_depth_aligned{curr_day_idx};
        template_depths_aligned = interp1(curr_csd_depth,curr_csd_depth_aligned,template_depths);
        
        % Find cortex end by largest gap between templates
        sorted_template_depths = sort([template_depths_aligned]);
        [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
        ctx_end = sorted_template_depths(max_gap_idx)+1;
        
        ctx_depth = [sorted_template_depths(1),ctx_end];
        ctx_units = template_depths_aligned <= ctx_depth(2);
        
        %%% GET FLUORESCENCE AND SPIKES BY DEPTH
        
        % Set binning time
        skip_seconds = 60;
        spike_binning_t = 1/framerate; % seconds
        spike_binning_t_edges = frame_t(1)+skip_seconds:spike_binning_t:frame_t(end)-skip_seconds;
        spike_binning_t_centers = spike_binning_t_edges(1:end-1) + diff(spike_binning_t_edges)/2;
        
        % Get fluorescence in pre-drawn ROI
        curr_ctx_roi = vis_ctx_ephys(curr_animal_idx).ctx_roi{curr_day_idx};
        
        fVdf_deconv = AP_deconv_wf(fVdf);
        fluor_roi = AP_svd_roi(Udf,fVdf_deconv,avg_im,[],curr_ctx_roi);
        fluor_roi_interp = interp1(frame_t,fluor_roi,spike_binning_t_centers);
        
        % Set sliding depth window of MUA
        depth_corr_range = [-200,1500];
        depth_corr_window = 200; % MUA window in microns
        depth_corr_window_spacing = 50; % MUA window spacing in microns
        
        depth_corr_bins = [depth_corr_range(1):depth_corr_window_spacing:(depth_corr_range(2)-depth_corr_window); ...
            (depth_corr_range(1):depth_corr_window_spacing:(depth_corr_range(2)-depth_corr_window))+depth_corr_window];
        depth_corr_bin_centers = depth_corr_bins(1,:) + diff(depth_corr_bins,[],1)/2;
        
        cortex_mua = zeros(size(depth_corr_bins,2),length(spike_binning_t_centers));
        for curr_depth = 1:size(depth_corr_bins,2)
            curr_depth_templates_idx = ...
                find(ctx_units & ...
                template_depths_aligned >= depth_corr_bins(1,curr_depth) & ...
                template_depths_aligned < depth_corr_bins(2,curr_depth));
            
            cortex_mua(curr_depth,:) = histcounts(spike_times_timeline( ...
                ismember(spike_templates,curr_depth_templates_idx)),spike_binning_t_edges);
        end
        
        % Plot cortex units and depth bins
        norm_spike_n = mat2gray(log10(accumarray(spike_templates,1)+1));
        figure; set(gca,'YDir','reverse'); hold on;
        plot(norm_spike_n(ctx_units),template_depths_aligned(ctx_units),'.b','MarkerSize',20);
        plot(norm_spike_n(~ctx_units),template_depths_aligned(~ctx_units),'.k','MarkerSize',20);
        line(xlim,repmat(depth_corr_range(1),2,1),'color','r','linewidth',2);
        line(xlim,repmat(depth_corr_range(2),2,1),'color','r','linewidth',2);
        for i = 1:length(depth_corr_bin_centers)
            line(xlim,repmat(depth_corr_bin_centers(i),2,1),'color','r');
        end
        xlabel('Norm firing rate');
        ylabel('Aligned depth');
        title([animal ' ' day]);
        drawnow;
        
        %%% LOAD STRIATUM EPHYS AND GET MUA BY DEPTH
        clear load_parts
        load_parts.ephys = true;
        site = 1; % (striatum is always on probe 1)
        str_align = 'kernel';
        AP_load_experiment;
        
        striatum_mua = nan(n_aligned_depths,length(spike_binning_t_centers));
        for curr_depth = 1:n_aligned_depths
            curr_spike_times = spike_times_timeline(aligned_str_depth_group == curr_depth);
            % Skip if no spikes at this depth
            if isempty(curr_spike_times)
                continue
            end
            striatum_mua(curr_depth,:) = histcounts(curr_spike_times,spike_binning_t_edges);
        end
        
        %%% REGULAR CORRELATION
        cortex_fluor_corr = 1-pdist2(cortex_mua,fluor_roi_interp,'correlation');
        cortex_striatum_corr = 1-pdist2(cortex_mua,striatum_mua,'correlation');
        
        %%% PACKAGE AND SAVE
        data(curr_animal).cortex_mua_depth{curr_day} = depth_corr_bin_centers;
        data(curr_animal).cortex_fluor_corr{curr_day} = cortex_fluor_corr;
        data(curr_animal).cortex_striatum_corr{curr_day} = cortex_striatum_corr;
        
        % Clear variables for next experiment
        clearvars -except animals vis_ctx_ephys ...
            curr_animal animal ...
            protocol flexible_name experiments curr_exp ...
            data
        
    end
    
end

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
data_fn = [data_path filesep 'ctx_fluor_mua_corr_' protocol];
save(data_fn,'data');
disp(['Saved ' data_fn]);

%% Get activity correlation cortex wf/mua and ctx ephys/str mua (passive)

animals = {'AP043','AP060','AP061'};

% Load ephys alignment
vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
load(vis_ctx_ephys_fn);

data = struct;
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    protocol = 'AP_sparseNoise';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    for curr_day = 1:length(experiments)
        
        % Load experiment
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(1); % (only one doubled experiment, and the second one was worse)
        site = 2; % cortex probe is always site 2
        lfp_channel = 'all';
        str_align = 'none'; % (because this is on cortex)
        verbose = false;
        AP_load_experiment
        
        %%% DEPTH-ALIGN TEMPLATES, FIND CORTEX BOUNDARY
        curr_animal_idx = strcmp(animal,{vis_ctx_ephys.animal});
        curr_day_idx = strcmp(day,vis_ctx_ephys(curr_animal_idx).day);
        curr_csd_depth = vis_ctx_ephys(curr_animal_idx).stim_csd_depth{curr_day_idx};
        curr_csd_depth_aligned = vis_ctx_ephys(curr_animal_idx).stim_csd_depth_aligned{curr_day_idx};
        template_depths_aligned = interp1(curr_csd_depth,curr_csd_depth_aligned,template_depths);
        
        % Find cortex end by largest gap between templates
        sorted_template_depths = sort([template_depths_aligned]);
        [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
        ctx_end = sorted_template_depths(max_gap_idx)+1;
        
        ctx_depth = [sorted_template_depths(1),ctx_end];
        ctx_units = template_depths_aligned <= ctx_depth(2);
        
        %%% GET FLUORESCENCE AND SPIKES BY DEPTH
        
        % Set binning time
        skip_seconds = 60;
        spike_binning_t = 1/framerate; % seconds
        spike_binning_t_edges = frame_t(1)+skip_seconds:spike_binning_t:frame_t(end)-skip_seconds;
        spike_binning_t_centers = spike_binning_t_edges(1:end-1) + diff(spike_binning_t_edges)/2;
        
        % Get fluorescence in pre-drawn ROI
        curr_ctx_roi = vis_ctx_ephys(curr_animal_idx).ctx_roi{curr_day_idx};
        
        fVdf_deconv = AP_deconv_wf(fVdf);
        fluor_roi = AP_svd_roi(Udf,fVdf_deconv,avg_im,[],curr_ctx_roi);
        fluor_roi_interp = interp1(frame_t,fluor_roi,spike_binning_t_centers);
        
        % Set sliding depth window of MUA
        depth_corr_range = [-200,1500];
        depth_corr_window = 200; % MUA window in microns
        depth_corr_window_spacing = 50; % MUA window spacing in microns
        
        depth_corr_bins = [depth_corr_range(1):depth_corr_window_spacing:(depth_corr_range(2)-depth_corr_window); ...
            (depth_corr_range(1):depth_corr_window_spacing:(depth_corr_range(2)-depth_corr_window))+depth_corr_window];
        depth_corr_bin_centers = depth_corr_bins(1,:) + diff(depth_corr_bins,[],1)/2;
        
        cortex_mua = zeros(size(depth_corr_bins,2),length(spike_binning_t_centers));
        for curr_depth = 1:size(depth_corr_bins,2)
            curr_depth_templates_idx = ...
                find(ctx_units & ...
                template_depths_aligned >= depth_corr_bins(1,curr_depth) & ...
                template_depths_aligned < depth_corr_bins(2,curr_depth));
            
            cortex_mua(curr_depth,:) = histcounts(spike_times_timeline( ...
                ismember(spike_templates,curr_depth_templates_idx)),spike_binning_t_edges);
        end
        
        % Plot cortex units and depth bins
        norm_spike_n = mat2gray(log10(accumarray(spike_templates,1)+1));
        figure; set(gca,'YDir','reverse'); hold on;
        plot(norm_spike_n(ctx_units),template_depths_aligned(ctx_units),'.b','MarkerSize',20);
        plot(norm_spike_n(~ctx_units),template_depths_aligned(~ctx_units),'.k','MarkerSize',20);
        line(xlim,repmat(depth_corr_range(1),2,1),'color','r','linewidth',2);
        line(xlim,repmat(depth_corr_range(2),2,1),'color','r','linewidth',2);
        for i = 1:length(depth_corr_bin_centers)
            line(xlim,repmat(depth_corr_bin_centers(i),2,1),'color','r');
        end
        xlabel('Norm firing rate');
        ylabel('Aligned depth');
        title([animal ' ' day]);
        drawnow;
        
        %%% LOAD STRIATUM EPHYS AND GET MUA BY DEPTH
        clear load_parts
        load_parts.ephys = true;
        site = 1; % (striatum is always on probe 1)
        str_align = 'kernel';
        AP_load_experiment;
        
        striatum_mua = nan(n_aligned_depths,length(spike_binning_t_centers));
        for curr_depth = 1:n_aligned_depths
            curr_spike_times = spike_times_timeline(aligned_str_depth_group == curr_depth);
            % Skip if no spikes at this depth
            if isempty(curr_spike_times)
                continue
            end
            striatum_mua(curr_depth,:) = histcounts(curr_spike_times,spike_binning_t_edges);
        end
        
        %%% REGULAR CORRELATION
        cortex_fluor_corr = 1-pdist2(cortex_mua,fluor_roi_interp,'correlation');
        cortex_striatum_corr = 1-pdist2(cortex_mua,striatum_mua,'correlation');
        
        %%% PACKAGE AND SAVE
        data(curr_animal).cortex_mua_depth{curr_day} = depth_corr_bin_centers;
        data(curr_animal).cortex_fluor_corr{curr_day} = cortex_fluor_corr;
        data(curr_animal).cortex_striatum_corr{curr_day} = cortex_striatum_corr;
        
        % Clear variables for next experiment
        clearvars -except animals vis_ctx_ephys ...
            curr_animal animal ...
            protocol flexible_name experiments curr_exp ...
            data
        
    end
    
end

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
data_fn = [data_path filesep 'ctx_fluor_mua_corr_' protocol];
save(data_fn,'data');

%% Get activity correlation cortex wf/mua and ctx ephys/str mua (task depth deconv)

animals = {'AP043','AP060','AP061'};

% Load ephys alignment
vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
load(vis_ctx_ephys_fn);

data = struct;
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    protocol = 'vanillaChoiceworld';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    for curr_day = 1:length(experiments)
        
        % Set variables to keep
        preload_vars = who;
        
        % Load experiment
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(1); % (only one doubled experiment, and the second one was worse)
        site = 2; % cortex probe is always site 2
        lfp_channel = 'all';
        str_align = 'none'; % (because this is on cortex)
        verbose = false;
        AP_load_experiment
        
        %%% DEPTH-ALIGN TEMPLATES, FIND CORTEX BOUNDARY
        curr_animal_idx = strcmp(animal,{vis_ctx_ephys.animal});
        curr_day_idx = strcmp(day,vis_ctx_ephys(curr_animal_idx).day);
        curr_csd_depth = vis_ctx_ephys(curr_animal_idx).stim_csd_depth{curr_day_idx};
        curr_csd_depth_aligned = vis_ctx_ephys(curr_animal_idx).stim_csd_depth_aligned{curr_day_idx};
        template_depths_aligned = interp1(curr_csd_depth,curr_csd_depth_aligned,template_depths);
        
        % Find cortex end by largest gap between templates
        sorted_template_depths = sort([template_depths_aligned]);
        [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
        ctx_end = sorted_template_depths(max_gap_idx)+1;
        
        ctx_depth = [sorted_template_depths(1),ctx_end];
        ctx_units = template_depths_aligned <= ctx_depth(2);
        
        %%% GET FLUORESCENCE AND SPIKES BY DEPTH
        
        % Set binning time
        skip_seconds = 60;
        spike_binning_t = 1/framerate; % seconds
        spike_binning_t_edges = frame_t(1)+skip_seconds:spike_binning_t:frame_t(end)-skip_seconds;
        spike_binning_t_centers = spike_binning_t_edges(1:end-1) + diff(spike_binning_t_edges)/2;
        
        % Get fluorescence in pre-drawn ROI
        curr_ctx_roi = vis_ctx_ephys(curr_animal_idx).ctx_roi{curr_day_idx};
        
        % Deconvolve fluorescence with cortex kernels by depth
        % (copied from AP_deconv_wf but manually using other kernel)
        
        % (load depth-specific deconv kernel, average)
        deconv_k_depth_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\gcamp6s_kernel_ctxdepth.mat';
        load(deconv_k_depth_fn)
        kernel_cat = cat(3,gcamp6s_kernel_ctxdepth.regression{:});
        kernel_mean = nanmean(kernel_cat./max(abs(kernel_cat),[],2),3);
        kernel = kernel_mean./sum(kernel_mean.^2,2);              
        
        % (loop through depth kernels and convolve)
        fluor_roi_interp = nan(size(kernel,1),length(spike_binning_t_centers));
        for curr_kernel_depth = 1:size(kernel,1)
            
            % Remove NaNs for convolution (put them back after)
            activity = fVdf;
            activity_nonan = activity;
            activity_nonan(isnan(activity_nonan)) = 0;           
            
            % Deconvolve to same size with replication padding
            convolved_activity = convn(padarray(activity_nonan, ...
                [0,floor(length(kernel(curr_kernel_depth,:))/2)],'replicate','both'), ...
                kernel(curr_kernel_depth,:),'valid');
            
            % Put original NaNs back in
            convolved_activity(isnan(activity)) = NaN;
            
            % Get ROI fluorescence and interpolate time
            fluor_roi = AP_svd_roi(Udf,convolved_activity,avg_im,[],curr_ctx_roi);
            fluor_roi_interp(curr_kernel_depth,:) = ...
                interp1(frame_t,fluor_roi,spike_binning_t_centers);
                    
        end             
        
        % Set sliding depth window of MUA
        depth_corr_range = [-200,1500];
        depth_corr_window = 200; % MUA window in microns
        depth_corr_window_spacing = 50; % MUA window spacing in microns
        
        depth_corr_bins = [depth_corr_range(1):depth_corr_window_spacing:(depth_corr_range(2)-depth_corr_window); ...
            (depth_corr_range(1):depth_corr_window_spacing:(depth_corr_range(2)-depth_corr_window))+depth_corr_window];
        depth_corr_bin_centers = depth_corr_bins(1,:) + diff(depth_corr_bins,[],1)/2;
        
        cortex_mua = zeros(size(depth_corr_bins,2),length(spike_binning_t_centers));
        for curr_depth = 1:size(depth_corr_bins,2)
            curr_depth_templates_idx = ...
                find(ctx_units & ...
                template_depths_aligned >= depth_corr_bins(1,curr_depth) & ...
                template_depths_aligned < depth_corr_bins(2,curr_depth));
            
            cortex_mua(curr_depth,:) = histcounts(spike_times_timeline( ...
                ismember(spike_templates,curr_depth_templates_idx)),spike_binning_t_edges);
        end
        
        % Plot cortex units and depth bins
        norm_spike_n = mat2gray(log10(accumarray(spike_templates,1)+1));
        figure; set(gca,'YDir','reverse'); hold on;
        plot(norm_spike_n(ctx_units),template_depths_aligned(ctx_units),'.b','MarkerSize',20);
        plot(norm_spike_n(~ctx_units),template_depths_aligned(~ctx_units),'.k','MarkerSize',20);
        line(xlim,repmat(depth_corr_range(1),2,1),'color','r','linewidth',2);
        line(xlim,repmat(depth_corr_range(2),2,1),'color','r','linewidth',2);
        for i = 1:length(depth_corr_bin_centers)
            line(xlim,repmat(depth_corr_bin_centers(i),2,1),'color','r');
        end
        xlabel('Norm firing rate');
        ylabel('Aligned depth');
        title([animal ' ' day]);
        drawnow;
        
        %%% LOAD STRIATUM EPHYS AND GET MUA BY DEPTH
        clear load_parts
        load_parts.ephys = true;
        site = 1; % (striatum is always on probe 1)
        str_align = 'kernel';
        AP_load_experiment;
        
        striatum_mua = nan(n_aligned_depths,length(spike_binning_t_centers));
        for curr_depth = 1:n_aligned_depths
            curr_spike_times = spike_times_timeline(aligned_str_depth_group == curr_depth);
            % Skip if no spikes at this depth
            if isempty(curr_spike_times)
                continue
            end
            striatum_mua(curr_depth,:) = histcounts(curr_spike_times,spike_binning_t_edges);
        end
        
        %%% REGULAR CORRELATION
        cortex_fluor_corr = 1-pdist2(cortex_mua,fluor_roi_interp,'correlation');
        cortex_striatum_corr = 1-pdist2(cortex_mua,striatum_mua,'correlation');
        
        %%% PACKAGE AND SAVE
        data(curr_animal).cortex_mua_depth{curr_day} = depth_corr_bin_centers;
        data(curr_animal).cortex_fluor_corr{curr_day} = cortex_fluor_corr;
        data(curr_animal).cortex_striatum_corr{curr_day} = cortex_striatum_corr;
        
        % Clear variables for next experiment
        clearvars('-except',preload_vars{:})
        
    end
    
end

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
data_fn = [data_path filesep 'ctx_fluor_mua_corr_' protocol '_ctxdepthkernel'];
save(data_fn,'data');
disp(['Saved ' data_fn]);


%% Task trial activity (ctx/str MUA separately)
clear all

animals = {'AP043','AP060','AP061'};
protocol = 'vanillaChoiceworld';
recording_site = {'ctxstrephys_str','ctxstrephys_ctx'};

% Load ephys alignment
vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
load(vis_ctx_ephys_fn);

for curr_site = 1:2
    
    % Initialize save variable
    trial_data_all = struct;
    
    for curr_animal = 1:length(animals)
        
        animal = animals{curr_animal};
        
        experiments = AP_find_experiments(animal,protocol);
        experiments = experiments([experiments.imaging] & [experiments.ephys]);
        
        disp(['Loading ' animal]);
        
        for curr_day = 1:length(experiments)
            
            day = experiments(curr_day).day;
            experiment = experiments(curr_day).experiment(1);
            
            % Load experiment
            site = curr_site;
            switch site
                case 1
                    % Site 1 = striautm
                    str_align = 'kernel';
                case 2
                    % Site 2 = cortex
                    str_align = 'none';
            end
            AP_load_experiment;
            
            % If cortical experiment - depth-align cortical MUA
            if site == 2
                % Get depth group boundaries
                % (these are just manual based on the CSD)
                n_aligned_depths = 2;
                ctx_depth_edges = linspace(0,1200,n_aligned_depths+1);
                
                % Depth-align templates
                curr_animal_idx = strcmp(animal,{vis_ctx_ephys.animal});
                curr_day_idx = strcmp(day,vis_ctx_ephys(curr_animal_idx).day);
                curr_csd_depth = vis_ctx_ephys(curr_animal_idx).stim_csd_depth{curr_day_idx};
                curr_csd_depth_aligned = vis_ctx_ephys(curr_animal_idx).stim_csd_depth_aligned{curr_day_idx};
                
                template_depths_aligned = interp1(curr_csd_depth,curr_csd_depth_aligned,template_depths);
                spike_depths_aligned = template_depths_aligned(spike_templates);
                
                % Find cortex end by largest gap between templates
                sorted_template_depths = sort([template_depths_aligned]);
                [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
                ctx_end = sorted_template_depths(max_gap_idx)+1;
                ctx_depth = [sorted_template_depths(1),ctx_end];
                ctx_units = template_depths_aligned <= ctx_depth(2);
                
                % Assign templates to depth groups
                % (NOTE: variable still called 'str' to make it easier to
                % fit into other code)
                ctx_spike_depths = spike_depths_aligned;
                ctx_spike_depths(spike_depths_aligned < ctx_depth(1) | spike_depths_aligned > ctx_depth(2)) = NaN;
                aligned_str_depth_group = discretize(ctx_spike_depths,ctx_depth_edges);
            end
            
            % Pull out trial data
            AP_ctx_str_grab_trial_data;
            
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
            
            % Clear for next loop
            clearvars -except ...
                vis_ctx_ephys ...
                recording_site curr_site ...
                animals curr_animal animal protocol experiments curr_day ...
                trial_data_all
            
            AP_print_progress_fraction(curr_day,length(experiments));
        end
    end
    
    clearvars -except ...
        vis_ctx_ephys ...
        recording_site curr_site animals protocol ...
        trial_data_all
    
    disp('Finished loading all')
    
    save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
    save_fn = ['trial_activity_' protocol '_' recording_site{curr_site}];
    save([save_path filesep save_fn],'-v7.3');
    disp(['Saved ' save_fn]);
    
end



%% Passive trial activity (ctx/str MUA separately)
clear all

animals = {'AP043','AP060','AP061'};
protocol = 'AP_lcrGratingPassive';
recording_site = {'ctxstrephys_str','ctxstrephys_ctx'};

% Load ephys alignment
vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
load(vis_ctx_ephys_fn);

for curr_site = 1:2
    
    % Initialize save variable
    trial_data_all = struct;
    
    for curr_animal = 1:length(animals)
        
        animal = animals{curr_animal};
        
        experiments = AP_find_experiments(animal,protocol);
        experiments = experiments([experiments.imaging] & [experiments.ephys]);
        
        disp(['Loading ' animal]);
        
        for curr_day = 1:length(experiments)
            
            day = experiments(curr_day).day;
            experiment = experiments(curr_day).experiment(1);
            
            % Load experiment
            site = curr_site;
            switch site
                case 1
                    % Site 1 = striautm
                    str_align = 'kernel';
                case 2
                    % Site 2 = cortex
                    str_align = 'none';
            end
            AP_load_experiment;
            
            % If cortical experiment - depth-align cortical MUA
            if site == 2
                % Get depth group boundaries
                % (these are just manual based on the CSD)
                n_aligned_depths = 2;
                ctx_depth_edges = linspace(0,1200,n_aligned_depths+1);
                
                % Depth-align templates
                curr_animal_idx = strcmp(animal,{vis_ctx_ephys.animal});
                curr_day_idx = strcmp(day,vis_ctx_ephys(curr_animal_idx).day);
                curr_csd_depth = vis_ctx_ephys(curr_animal_idx).stim_csd_depth{curr_day_idx};
                curr_csd_depth_aligned = vis_ctx_ephys(curr_animal_idx).stim_csd_depth_aligned{curr_day_idx};
                
                template_depths_aligned = interp1(curr_csd_depth,curr_csd_depth_aligned,template_depths);
                spike_depths_aligned = template_depths_aligned(spike_templates);
                
                % Find cortex end by largest gap between templates
                sorted_template_depths = sort([template_depths_aligned]);
                [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
                ctx_end = sorted_template_depths(max_gap_idx)+1;
                ctx_depth = [sorted_template_depths(1),ctx_end];
                ctx_units = template_depths_aligned <= ctx_depth(2);
                
                % Assign templates to depth groups
                % (NOTE: variable still called 'str' to make it easier to
                % fit into other code)
                ctx_spike_depths = spike_depths_aligned;
                ctx_spike_depths(spike_depths_aligned < ctx_depth(1) | spike_depths_aligned > ctx_depth(2)) = NaN;
                aligned_str_depth_group = discretize(ctx_spike_depths,ctx_depth_edges);
            end
            
            % Pull out trial data
            AP_ctx_str_grab_trial_data;
            
            % Store trial data into master structure
            trial_data_fieldnames = fieldnames(trial_data);
            for curr_trial_data_field = trial_data_fieldnames'
                trial_data_all.(cell2mat(curr_trial_data_field)){curr_animal,1}{curr_day,1} = ...
                    trial_data.(cell2mat(curr_trial_data_field));
            end
            
            % Store general info
            trial_data_all.animals = animals;
            trial_data_all.t = t;
            
            % Clear for next loop
            clearvars -except ...
                vis_ctx_ephys ...
                recording_site curr_site ...
                animals curr_animal animal protocol experiments curr_day ...
                trial_data_all
            
            AP_print_progress_fraction(curr_day,length(experiments));
        end
    end
    
    clearvars -except ...
        vis_ctx_ephys ...
        recording_site curr_site animals protocol ...
        trial_data_all
    
    disp('Finished loading all')
    
    save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
    save_fn = ['trial_activity_' protocol '_' recording_site{curr_site}];
    save([save_path filesep save_fn],'-v7.3');
    disp(['Saved ' save_fn]);
    
end


%% [[unused]] Cortex fluorescence > multiunit kernels 
% (did this first and then realized the one below makes more sense

clear all

% Parameters for regression
regression_params.use_svs = 1:100;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.1,0.1];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

% Load Master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');

% Load ephys alignment
vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
load(vis_ctx_ephys_fn);

animals = {'AP043','AP060','AP061'};
task_protocol = 'vanillaChoiceworld';
retinotopy_protocol = 'AP_sparseNoise';
site = 2;

ctx_fluorpred = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    disp(animal);
    
    % Get days with both task and retinotopy protocol
    % and that have imaging + ephys
    task_experiments = AP_find_experiments(animal,task_protocol);
    task_experiments = task_experiments([task_experiments.imaging] & [task_experiments.ephys]);
    
    retinotopy_experiments = AP_find_experiments(animal,retinotopy_protocol);
    retinotopy_experiments = retinotopy_experiments([retinotopy_experiments.imaging] & [retinotopy_experiments.ephys]);
    
    use_days = intersect({task_experiments.day},{retinotopy_experiments.day});
    % (use the first experiment of each protocol)
    use_experiments = ...
        [cellfun(@(x) x(1), {task_experiments(ismember({task_experiments.day},use_days)).experiment}'), ...
        cellfun(@(x) x(1), {retinotopy_experiments(ismember({retinotopy_experiments.day},use_days)).experiment}')];
    
    for curr_day = 1:length(use_days)
        for curr_exp = 1:size(use_experiments,2)
            
            % Set variables to keep
            preload_vars = who;
            
            % Load data and align striatum by depth
            day = use_days{curr_day};
            experiment = use_experiments(curr_day,curr_exp);
            str_align = 'none'; % (because it's cortex)
            load_parts.cam = false;
            load_parts.imaging = true;
            load_parts.ephys = true;
            AP_load_experiment;
            
            % Convert U to master U
            Udf_aligned = AP_align_widefield(Udf,animal,day);
            
            % Deconvolve fluoresence (for use in regression)
            fVdf_deconv = AP_deconv_wf(fVdf);
            fVdf_deconv(isnan(fVdf_deconv)) = 0;
            
            % Get time points to bin
            sample_rate = framerate*regression_params.upsample_factor;
            time_bins = frame_t(find(frame_t > ...
                regression_params.skip_seconds,1)):1/sample_rate: ...
                frame_t(find(frame_t-frame_t(end) < ...
                -regression_params.skip_seconds,1,'last'));
            time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
            
            % Resample deconvolved fluorescence
            fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
                      
            %%% Align cortical recording 
            
            % Get depth group boundaries
            % (these are just manual based on the CSD)
            n_aligned_depths = 2;
            ctx_depth_edges = linspace(0,1200,n_aligned_depths+1);
            
            % Depth-align templates
            curr_animal_idx = strcmp(animal,{vis_ctx_ephys.animal});
            curr_day_idx = strcmp(day,vis_ctx_ephys(curr_animal_idx).day);
            curr_csd_depth = vis_ctx_ephys(curr_animal_idx).stim_csd_depth{curr_day_idx};
            curr_csd_depth_aligned = vis_ctx_ephys(curr_animal_idx).stim_csd_depth_aligned{curr_day_idx};
            
            template_depths_aligned = interp1(curr_csd_depth,curr_csd_depth_aligned,template_depths);
            spike_depths_aligned = template_depths_aligned(spike_templates);
            
            % Find cortex end by largest gap between templates
            sorted_template_depths = sort([template_depths_aligned]);
            [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
            ctx_end = sorted_template_depths(max_gap_idx)+1;
            ctx_depth = [sorted_template_depths(1),ctx_end];
            ctx_units = template_depths_aligned <= ctx_depth(2);
            
            % Assign templates to depth groups
            ctx_spike_depths = spike_depths_aligned;
            ctx_spike_depths(spike_depths_aligned < ctx_depth(1) | spike_depths_aligned > ctx_depth(2)) = NaN;
            aligned_ctx_depth_group = discretize(ctx_spike_depths,ctx_depth_edges);
                        
            %%%
            
            % Bin spikes to match widefield frames
            binned_spikes = nan(n_aligned_depths,length(time_bin_centers));
            for curr_depth = 1:n_aligned_depths
                curr_spike_times = spike_times_timeline(aligned_ctx_depth_group == curr_depth);
                
                % (skip if no spikes at this depth)
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
            [ctx_fluor_k,ctxpred_spikes_std,explained_var] = ...
                AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
                binned_spikes_std,kernel_frames,lambda, ...
                regression_params.zs,regression_params.cvfold, ...
                true,regression_params.use_constant);
            
            % Recast the k's into the master U
            ctx_fluor_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
                reshape(ctx_fluor_k{1},size(ctx_fluor_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
                size(ctx_fluor_k{1}));
            
            % Re-scale the prediction (subtract offset, multiply, add scaled offset)
            ctxpred_spikes = (ctxpred_spikes_std - squeeze(ctx_fluor_k{end})).* ...
                nanstd(binned_spikes,[],2) + ...
                nanstd(binned_spikes,[],2).*squeeze(ctx_fluor_k{end});
            
            % Store
            ctx_fluorpred.animal{curr_animal} = animal;
            ctx_fluorpred.ctx{curr_animal}{curr_day,curr_exp} = binned_spikes;
            ctx_fluorpred.ctx_fluorpred{curr_animal}{curr_day,curr_exp} = ctxpred_spikes;
            ctx_fluorpred.ctx_fluor_k{curr_animal}{curr_day,curr_exp} = ctx_fluor_k_recast;
            
            % Reset for next experiment
            clearvars('-except',preload_vars{:})
            
        end
        AP_print_progress_fraction(curr_day,length(use_days));
    end
    
    % Save after each animal
    save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data'];
    save_fn = [save_path filesep 'ctx_fluorpred'];
    disp(['Saving ' save_fn]);
    save(save_fn,'ctx_fluorpred','-v7.3');
    disp('Saved.');
    
end

disp('Done.');

%% Multiunit and deconvolved fluorescence traces across protocols
% Copied from AP_deconv_wf_kernelfit: 
% that just saves the kernels in the github repo, 
% this saves the traces in the data folder to get explained variance

% Find recordings
animals = {'AP043','AP060','AP061'};

recordings = struct('animal',{},'day',{});

recording_idx = 1;
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    protocol = 'vanillaChoiceworld';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    for curr_day = 1:length(experiments)        
        recordings(recording_idx).animal = animal;
        recordings(recording_idx).day = experiments(curr_day).day; 
        recordings(recording_idx).site = 2; % Cortex probe is always 2
        recording_idx = recording_idx+1;
    end
    
end

% Loop through all recordings
ctx_deconv_traces = struct;
for curr_recording = 1:length(recordings)
       
    % Set varibles to keep
    preload_vars = who;
    
    animal = recordings(curr_recording).animal;
    day = recordings(curr_recording).day;
    site = recordings(curr_recording).site;
    
    % Load and concantenate data from an animal/day

    % Find all experiments for that day
    % (get folders with only a number - those're the experiment folders)
    experiments_dir = dir(AP_cortexlab_filename(animal,day,[],'expInfo'));
    experiments_num_idx = cellfun(@(x) ~isempty(x), regexp({experiments_dir.name},'^\d*$'));
    experiments = cellfun(@str2num,{experiments_dir(experiments_num_idx).name});
    
    % Loop through experiments, collate data
    upsample_factor = 1;
    skip_seconds = 60;
    
    protocols_all = cell(size(experiments));
    time_bin_centers_all = cell(size(experiments));
    fVdf_resample_all = cell(size(experiments));
    binned_spikes_all = cell(size(experiments));
    
    for curr_exp = 1:length(experiments)
        
        experiment = experiments(curr_exp);
        str_align = 'none';
        AP_load_experiment;
        
        protocols_all{curr_exp} = expDef;
        
        % Set cortical units
        
        % Find cortex end by largest gap between templates
        sorted_template_depths = sort([template_depths]);
        [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
        ctx_end = sorted_template_depths(max_gap_idx);
        
        ctx_depth = [sorted_template_depths(1),ctx_end];
        ctx_units = template_depths <= ctx_depth(2);       
        
        % Resample and concatenate data
        % Get time points to query
        sample_rate = framerate*upsample_factor;
        time_bins = frame_t(find(frame_t >= ...
            skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) <= ...
            -skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        time_bin_centers_all{curr_exp} = time_bin_centers;
        
        % Get upsampled dVdf's
        fVdf_resample_all{curr_exp} = interp1(frame_t,fVdf',time_bin_centers)';
        
        % Get binned cortex MUA  
        curr_exp_spikes = spike_depths >= ctx_depth(1) & spike_depths <= ctx_depth(2);        
        binned_spikes_all{curr_exp} = histcounts(spike_times_timeline(curr_exp_spikes),time_bins);   
        
    end
    
    % Concatenate all data
    time_bin_centers = cat(2,time_bin_centers_all{:});
    fVdf_resample = cat(2,fVdf_resample_all{:});
    binned_spikes = cat(2,binned_spikes_all{:});
    
    % Std-normalize spikes
    binned_spikes_std = binned_spikes/nanstd(binned_spikes);
    
 
    % Get fluorescence from pre-drawn ROI
    
    % Load ephys alignment with ROI (created in AP_ctx_str_trial_preprocessing)
    vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
    load(vis_ctx_ephys_fn);
    
    curr_animal_idx = strcmp(animal,{vis_ctx_ephys.animal});
    curr_day_idx = strcmp(day,vis_ctx_ephys(curr_animal_idx).day);
    curr_roi = vis_ctx_ephys(curr_animal_idx).ctx_roi{curr_day_idx};
         
    % Get fluorescence using average image and weights
    [fluor_trace,fluor_mask] = AP_svd_roi(Udf,fVdf_resample,[],[],curr_roi);

    % Get kernel from wf > mua and mua > "deconv" wf
    kernel_t = [-0.5,0.5];
    kernel_frames = floor(kernel_t(1)*sample_rate):ceil(kernel_t(2)*sample_rate);  
    zs = [false,false];
    cvfold = 5;
    return_constant = true;
    use_constant = true;
    lambda = 0;

    % Regression from ROI fluorescence to multiunit
    [curr_gcamp_kernel,predicted_spikes,explained_var_spikes] = ...
        AP_regresskernel(fluor_trace, ...
        binned_spikes_std,kernel_frames,lambda,zs, ...
        cvfold,return_constant,use_constant);

    % Store traces   
    ctx_deconv_traces(curr_recording).protocol = protocols_all';
    ctx_deconv_traces(curr_recording).t = time_bin_centers_all';
    % (split traces back into protocols)
    protocol_split_n = cellfun(@length,time_bin_centers_all);
    ctx_deconv_traces(curr_recording).fluor = ...
        mat2cell(fluor_trace',protocol_split_n,1);
    ctx_deconv_traces(curr_recording).spikes = ...
        mat2cell(binned_spikes_std',protocol_split_n,1);
    ctx_deconv_traces(curr_recording).fluor_deconv = ...
        mat2cell(predicted_spikes',protocol_split_n,1);
    
    AP_print_progress_fraction(curr_recording,length(recordings));
    
    % Reset for next recording
    clearvars('-except',preload_vars{:})
    
end

% Save
save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data'];
save_fn = [save_path filesep 'ctx_deconv_traces'];
disp(['Saving ' save_fn]);
save(save_fn,'ctx_deconv_traces','-v7.3');
disp('Saved.');


%% ~~~~~~~~~~~ CORTICOSTRIATAL  ~~~~~~~~~~~~~

%% Choiceworld trial activity (striatum depth)
clear all

lambda = 20; % Set lambda - not done empirically yet

animals = {'AP063','AP068'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(1);
        
        % Load experiment
        n_aligned_depths = 4;
        str_align = 'depth';
        AP_load_experiment;
        
        % Pull out trial data
        AP_ctx_str_grab_trial_data;
        
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
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars -except animals curr_animal animal protocol experiments curr_day ...
            trial_data_all lambda
        
    end
end

clearvars -except trial_data_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_corticostriatal'];
save([save_path filesep save_fn],'-v7.3');

%% Passive trial activity trained (striatum depth)
clear all

lambda = 20; % Set lambda - not done empirically yet

animals = {'AP063','AP068'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        n_aligned_depths = 4;
        str_align = 'depth';
        AP_load_experiment;
        
        % Pull out trial data
        AP_ctx_str_grab_trial_data;
        
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
        clearvars -except animals curr_animal animal protocol experiments curr_day ...
            trial_data_all lambda
        
    end
end

clearvars -except trial_data_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_AP_lcrGratingPassive_corticostriatal'];
save([save_path filesep save_fn],'-v7.3');

%% Passive trial activity naive (widefield only)
clear all

lambda = 20; % (only for wheel - unused)

animals = {'AP063','AP068'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & ~[experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        AP_load_experiment;
        
        % Pull out trial data
        AP_ctx_str_grab_trial_data;
        
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
        clearvars -except animals curr_animal animal protocol experiments curr_day ...
            trial_data_all lambda
        
    end
end

clearvars -except trial_data_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_AP_lcrGratingPassive_corticostriatal_naive'];
save([save_path filesep save_fn],'-v7.3');



%% ~~~~~~~~~~~ COMBINED/OTHER  ~~~~~~~~~~~~~


%% Histology vs widefield-estimated probe trajectory
% Copied/modified from AP_ctx2str_probe
% (this generates figures, but only to check: plotted in fig script)

% Set recordings with histology
recordings = struct('animal',{},'day',{});

recordings(1).animal = 'AP032';
recordings(1).day = '2018-10-26';

recordings(2).animal = 'AP033';
recordings(2).day = '2018-10-26';

recordings(3).animal = 'AP034';
recordings(3).day = '2018-10-26';

recordings(4).animal = 'AP036';
recordings(4).day = '2018-11-14';

recordings(5).animal = 'AP043';
recordings(5).day = '2019-12-09';

recordings(6).animal = 'AP060';
recordings(6).day = '2019-12-06';

recordings(7).animal = 'AP061';
recordings(7).day = '2019-12-09';

histology_probe_ccf_all = struct;
for curr_recording = 1:length(recordings)
    
    probe_angle = 45; % from horizontal
    
    animal = recordings(curr_recording).animal;
    day = recordings(curr_recording).day;
    [img_path,img_exists] = AP_cortexlab_filename(animal,day,[],'imaging');
    avg_im = readNPY([img_path filesep 'meanImage_blue.npy']);
    
    avg_im_aligned = AP_align_widefield(avg_im,animal,day);
    
    h = figure('units','normalized','outerposition',[0 0 1 1]);imagesc(avg_im_aligned);
    axis image off;
    colormap(gray);caxis([0,prctile(avg_im_aligned(:),95)]);
    title('Click probe start/end');
    probe_wf = ginput(2);
    close(h)
    
    % Load and invert master CCF tform
    ccf_tform_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\ccf_tform'];
    load(ccf_tform_fn);
    ccf_tform_inverse = invert(ccf_tform);
    
    % Convert aligned widefield probe points to CCF
    % (master alignment to downsampled CCF, need to upsample)
    um2pixel = 20.6;
    probe_ccf_tformed = [probe_wf,[1;1]]*ccf_tform_inverse.T;
    probe_ccf_wf = round(probe_ccf_tformed(:,1:2)*(um2pixel/10));
    
    % Plot probe points on CCF next to average image
    figure('Name',animal);
    
    subplot(1,2,1);
    imagesc(avg_im_aligned);
    axis image off;
    colormap(gray);caxis([0,prctile(avg_im_aligned(:),95)]);
    AP_reference_outline('ccf_aligned','r');
    line(probe_wf(:,1),probe_wf(:,2),'color','b','linewidth',1,'linestyle','--');
    title('Widefield');
    
    subplot(1,2,2); hold on; set(gca,'YDir','reverse');axis image;
    load('C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\cortical_area_boundaries.mat');
    for curr_area_idx =1:length(cortical_area_boundaries)
        p = cellfun(@(outline) plot(outline(:,2),outline(:,1),'color','k'), ...
            cortical_area_boundaries{curr_area_idx},'uni',false);
    end
    line(probe_ccf_wf(:,1),probe_ccf_wf(:,2),'color','r','linewidth',2);
    title('CCF');
    
    % Load in the annotated Allen volume and names
    allen_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
    av = readNPY([allen_path filesep 'annotation_volume_10um_by_index.npy']);
    st = loadStructureTree([allen_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
    
    % Estimate probe entry point by clicked probe start
    depth_start = find(av(probe_ccf_wf(1,2),:,probe_ccf_wf(1,1)) > 1,1);
    probe_entry_ccf = [probe_ccf_wf(1,2),depth_start,probe_ccf_wf(1,1)];
    
    % Estimate height of clicked point from user-set probe angle
    probe_sample_length = pdist2(probe_ccf_wf(1,:),probe_ccf_wf(2,:));
    probe_sample_height = round(probe_sample_length/tand(90-probe_angle));
    probe_air_ccf = [probe_ccf_wf(2,2),depth_start-probe_sample_height,probe_ccf_wf(2,1)];
    
    % Concatenate probe CCF coordinates (in up-down direction);
    probe_ccf_wf = [probe_air_ccf;probe_entry_ccf];
    
    % Get estimated probe vector (widefield)
    r0 = mean(probe_ccf_wf,1);
    xyz = bsxfun(@minus,probe_ccf_wf,r0);
    [~,~,V] = svd(xyz,0);
    probe_direction = V(:,1);
    
    probe_vector_evaluate = [0,sign(probe_direction(2))*1000];
    probe_vector_ccf = round(bsxfun(@plus,bsxfun(@times,probe_vector_evaluate',probe_direction'),r0));
    
    % Get estimated probe vector (histology)
    [probe_filename,probe_filename_exists] = AP_cortexlab_filename(animal,[],[],'probe_ccf');
    load(probe_filename);
    
    use_probe = 1;
    histology_points = probe_ccf(use_probe).points;
    r0 = mean(histology_points,1);
    xyz = bsxfun(@minus,histology_points,r0);
    [~,~,V] = svd(xyz,0);
    histology_probe_direction = V(:,1);
    
    probe_vector_evaluate = [-sign(probe_direction(2))*500,sign(probe_direction(2))*500];
    probe_vector = bsxfun(@plus,bsxfun(@times,probe_vector_evaluate',histology_probe_direction'),r0);
    probe_vector_ccf_histology = round(probe_vector);
    
    % histology done looking A->P so flipped, mirror about bregma
    bregma = allenCCFbregma;
    probe_vector_ccf_histology_mirrored = probe_vector_ccf_histology;
    probe_vector_ccf_histology_mirrored(:,3) = ...
        bregma(3) - (probe_vector_ccf_histology(:,3) - bregma(3));
    
    % Plot estimated probe location on CCF
    % (note the CCF is rotated to allow for dim 1 = x)
    h = figure('Name',animal); ccf_axes = axes; hold on
    slice_spacing = 10;
    target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) > 1,[2,1,3]);
    structure_patch = isosurface(target_volume,0);
    structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
    target_structure_color = [0.7,0.7,0.7];
    brain_outline = patch('Vertices',structure_wire.vertices*slice_spacing, ...
        'Faces',structure_wire.faces, ...
        'FaceColor','none','EdgeColor',target_structure_color);
    
    str_id = find(strcmp(st.safe_name,'Caudoputamen'));
    target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) == str_id,[2,1,3]);
    structure_patch = isosurface(target_volume,0);
    structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
    target_structure_color = [0.7,0,0.7];
    striatum_outline = patch('Vertices',structure_wire.vertices*slice_spacing, ...
        'Faces',structure_wire.faces, ...
        'FaceColor','none','EdgeColor',target_structure_color);
    
    axis image vis3d off;
    view([-30,25]);
    cameratoolbar(h,'SetCoordSys','y');
    cameratoolbar(h,'SetMode','orbit');
    title([animal ' ' day]);
    
    % Plot current probe vector
    probe_wf_line = line(ccf_axes,probe_vector_ccf(:,1),probe_vector_ccf(:,2),probe_vector_ccf(:,3),'color','k','linewidth',3);
    % (make sure histology probe enters on the left)
    if  diff(probe_vector_ccf_histology(:,3)) < 0
        probe_histology_line = line(ccf_axes, ...
            probe_vector_ccf_histology(:,1), ...
            probe_vector_ccf_histology(:,2), ...
            probe_vector_ccf_histology(:,3),'color','r','linewidth',3);
    else
        probe_histology_line = line(ccf_axes, ...
            probe_vector_ccf_histology_mirrored(:,1), ...
            probe_vector_ccf_histology_mirrored(:,2), ...
            probe_vector_ccf_histology_mirrored(:,3),'color','r','linewidth',3);
    end
    drawnow;   
    
    legend([probe_wf_line,probe_histology_line],{'Widefield-estimated','Histology'},'location','nw');
    
    % Store 
    histology_probe_ccf_all(curr_recording).im = avg_im_aligned;
    histology_probe_ccf_all(curr_recording).probe_wf = probe_wf;
    histology_probe_ccf_all(curr_recording).probe_ccf_wf = probe_ccf_wf;
    histology_probe_ccf_all(curr_recording).probe_ccf_histology = probe_ccf(use_probe).points;
    
    clearvars -except recordings curr_recording histology_probe_ccf_all
    
end

% Save
save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data'];
save_fn = [save_path filesep 'histology_probe_ccf_all'];
disp(['Saving ' save_fn]);
save(save_fn,'histology_probe_ccf_all');



%% Estimate probe trajectory from widefield image
% Copied/modified from AP_ctx2str_probe
% (this generates figures, but only to check: plotted in fig script)

probe_angle = 45; % from horizontal

protocol = 'vanillaChoiceworld';
flexible_name = true;

animals = { ...
    'AP024','AP025','AP026','AP027','AP028','AP029', ... % Original 
    'AP045','AP054','AP055','AP053','AP047','AP048', ... % Muscimol
    'AP043','AP060','AP061'};                            % Ctx ephys

% Load in the annotated Allen volume and names
allen_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
av = readNPY([allen_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean

% Plot brain to overlay probes
% (note the CCF is rotated to allow for dim 1 = x)
h = figure; 
ccf_axes = subplot(2,2,1); hold on
sagittal_axes = subplot(2,2,2,'YDir','reverse'); hold on;
axis image off;
horizontal_axes = subplot(2,2,3,'YDir','reverse'); hold on;
axis image off;
coronal_axes = subplot(2,2,4,'YDir','reverse'); hold on;
axis image off;

% Plot 1 = 3D
% (Use wire mesh - can add other structures)
slice_spacing = 10;

str_id = find(strcmp(st.safe_name,'Caudoputamen'));
target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) == str_id,[2,1,3]);
structure_patch = isosurface(target_volume,0);
structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
target_structure_color = [0,0,1];
striatum_outline = patch(ccf_axes, ...
    'Vertices',structure_wire.vertices*slice_spacing, ...
    'Faces',structure_wire.faces, ...
    'FaceAlpha',0.5,'FaceColor',target_structure_color,'EdgeColor','none');

% (Use Nick's outline, but rotate to make dim 1 = x)
brainGridData = readNPY([fileparts(which('plotBrainGrid')) filesep 'brainGridData.npy']);
plotBrainGrid(brainGridData(:,[1,3,2]),ccf_axes);

axes(ccf_axes);
axis image vis3d off;
view([-30,25]);
cameratoolbar(h,'SetCoordSys','y');
cameratoolbar(h,'SetMode','orbit');

% Plots 2-4 = projections
coronal_outline = bwboundaries(permute((max(av,[],1)) > 1,[2,3,1]));
horizontal_outline = bwboundaries(permute((max(av,[],2)) > 1,[3,1,2]));
sagittal_outline = bwboundaries(permute((max(av,[],3)) > 1,[2,1,3]));

str_id = find(strcmp(st.safe_name,'Caudoputamen'));
str_coronal_outline = bwboundaries(permute((max(av == str_id,[],1)) > 0,[2,3,1]));
str_horizontal_outline = bwboundaries(permute((max(av == str_id,[],2)) > 0,[3,1,2]));
str_sagittal_outline = bwboundaries(permute((max(av == str_id,[],3)) > 0,[2,1,3]));

plot(sagittal_axes,sagittal_outline{1}(:,2),sagittal_outline{1}(:,1),'k','linewidth',2);
plot(sagittal_axes,str_sagittal_outline{1}(:,2),str_sagittal_outline{1}(:,1),'b','linewidth',2);
axis image off;

plot(horizontal_axes,horizontal_outline{1}(:,2),horizontal_outline{1}(:,1),'k','linewidth',2);
plot(horizontal_axes,str_horizontal_outline{1}(:,2),str_horizontal_outline{1}(:,1),'b','linewidth',2);
plot(horizontal_axes,str_horizontal_outline{2}(:,2),str_horizontal_outline{2}(:,1),'b','linewidth',2);
axis image off;

plot(coronal_axes,coronal_outline{1}(:,2),coronal_outline{1}(:,1),'k','linewidth',2);
plot(coronal_axes,str_coronal_outline{1}(:,2),str_coronal_outline{1}(:,1),'b','linewidth',2);
plot(coronal_axes,str_coronal_outline{2}(:,2),str_coronal_outline{2}(:,1),'b','linewidth',2);
axis image off;

if length(animals) <= 17
    probe_color = [brewermap(8,'Set1');brewermap(8,'Set1')];
else
    error('More animals than colors: add another set');
end

probe_ccf_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        
        [img_path,img_exists] = AP_cortexlab_filename(animal,day,[],'imaging');
        avg_im = readNPY([img_path filesep 'meanImage_blue.npy']);
        
        avg_im_aligned = AP_align_widefield(avg_im,animal,day);
        
        h = figure('units','normalized','outerposition',[0 0 1 1]);
        imagesc(avg_im_aligned);
        axis image off;
        colormap(gray);caxis([0,prctile(avg_im_aligned(:),95)]);
        title(['Click probe bottom/top - ' ...
            num2str(curr_animal) '/' num2str(length(animals)) ': ' ...
            num2str(curr_day) '/' num2str(length(experiments))]);
        probe_wf = ginput(2);
        close(h)
        
        % Load and invert master CCF tform
        ccf_tform_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\ccf_tform'];
        load(ccf_tform_fn);
        ccf_tform_inverse = invert(ccf_tform);
        
        % Convert aligned widefield probe points to CCF
        % (master alignment to downsampled CCF, need to upsample)
        um2pixel = 20.6;
        probe_ccf_tformed = [probe_wf,[1;1]]*ccf_tform_inverse.T;
        probe_ccf = round(probe_ccf_tformed(:,1:2)*(um2pixel/10));
        
        % Estimate probe entry point by clicked probe start
        depth_start = find(av(probe_ccf(1,2),:,probe_ccf(1,1)) > 1,1);
        probe_entry_ccf = [probe_ccf(1,2),depth_start,probe_ccf(1,1)];
        
        % Estimate height of clicked point from user-set probe angle
        probe_sample_length = pdist2(probe_ccf(1,:),probe_ccf(2,:));
        probe_sample_height = round(probe_sample_length/tand(90-probe_angle));
        probe_air_ccf = [probe_ccf(2,2),depth_start-probe_sample_height,probe_ccf(2,1)];
        
        % Concatenate probe CCF coordinates (in up-down direction);
        probe_ccf = [probe_air_ccf;probe_entry_ccf];
        
        % Get estimated probe vector
        r0 = mean(probe_ccf,1);
        xyz = bsxfun(@minus,probe_ccf,r0);
        [~,~,V] = svd(xyz,0);
        probe_direction = V(:,1);
        
        probe_vector_evaluate = [0,sign(probe_direction(2))*1000];
        probe_vector_draw = round(bsxfun(@plus,bsxfun(@times,probe_vector_evaluate',probe_direction'),r0));
        
        % Plot current probe vector
        line(ccf_axes,probe_vector_draw(:,1),probe_vector_draw(:,2),probe_vector_draw(:,3), ...
            'color',probe_color(curr_animal,:),'linewidth',2);
        line(sagittal_axes,probe_vector_draw(:,1),probe_vector_draw(:,2), ...
            'color',probe_color(curr_animal,:),'linewidth',2);
        line(horizontal_axes,probe_vector_draw(:,1),probe_vector_draw(:,3), ...
            'color',probe_color(curr_animal,:),'linewidth',2);
        line(coronal_axes,probe_vector_draw(:,3),probe_vector_draw(:,2), ...
            'color',probe_color(curr_animal,:),'linewidth',2);
        
        drawnow;
        
        % Store probe vector
        probe_ccf_all{curr_animal}(:,:,curr_day) = probe_ccf;
        
    end
end

% Save
save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data'];
save_fn = [save_path filesep 'probe_ccf_all'];
disp(['Saving ' save_fn]);
save(save_fn,'probe_ccf_all');



%% Striatum cortical kernels (master U, paired task/retinotopy)

clear all

disp('Cortex/striatum, activity/kernels, task/notask');

% Parameters for regression
regression_params.use_svs = 1:100;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.1,0.1];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

% Load Master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');

n_cohorts = 3;
str_ctxpred = struct;

for curr_cohort = 1:n_cohorts
    
    switch curr_cohort
        case 1
            animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
            task_protocol = 'vanillaChoiceworld';
            retinotopy_protocol = 'stimSparseNoiseUncorrAsync';
            site = [];
        case 2
            animals = {'AP045','AP054','AP055','AP053','AP047','AP048'};
            task_protocol = 'vanillaChoiceworldNoRepeats';
            retinotopy_protocol = 'AP_sparseNoise';
            site = [];
        case 3
            animals = {'AP043','AP060','AP061'};
            task_protocol = 'vanillaChoiceworld';
            retinotopy_protocol = 'AP_sparseNoise';
            site = 1;
    end
        
    for curr_animal = 1:length(animals)
        
        animal = animals{curr_animal};
        disp(animal);
        
        % Get days with both task and retinotopy protocol
        % and that have imaging + ephys
        task_experiments = AP_find_experiments(animal,task_protocol);
        task_experiments = task_experiments([task_experiments.imaging] & [task_experiments.ephys]);
        
        retinotopy_experiments = AP_find_experiments(animal,retinotopy_protocol);
        retinotopy_experiments = retinotopy_experiments([retinotopy_experiments.imaging] & [retinotopy_experiments.ephys]);
        
        use_days = intersect({task_experiments.day},{retinotopy_experiments.day});
        % (use the first experiment of each protocol)
        use_experiments = ...
            [cellfun(@(x) x(1), {task_experiments(ismember({task_experiments.day},use_days)).experiment}'), ...
            cellfun(@(x) x(1), {retinotopy_experiments(ismember({retinotopy_experiments.day},use_days)).experiment}')];
        
        for curr_day = 1:length(use_days)
            for curr_exp = 1:size(use_experiments,2)
                
                % Set variables to keep
                preload_vars = who;
                
                % Load data and align striatum by depth
                day = use_days{curr_day};
                experiment = use_experiments(curr_day,curr_exp);
                str_align = 'kernel';
                load_parts.cam = false;
                load_parts.imaging = true;
                load_parts.ephys = true;
                AP_load_experiment;
                
                % Convert U to master U
                Udf_aligned = AP_align_widefield(Udf,animal,day);
                
                % Deconvolve fluoresence (for use in regression)
                fVdf_deconv = AP_deconv_wf(fVdf);
                fVdf_deconv(isnan(fVdf_deconv)) = 0;
                
                % Get time points to bin
                sample_rate = framerate*regression_params.upsample_factor;
                time_bins = frame_t(find(frame_t > ...
                    regression_params.skip_seconds,1)):1/sample_rate: ...
                    frame_t(find(frame_t-frame_t(end) < ...
                    -regression_params.skip_seconds,1,'last'));
                time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
                
                % Resample deconvolved fluorescence
                fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
                
                % Bin spikes to match widefield frames
                binned_spikes = nan(n_aligned_depths,length(time_bin_centers));
                for curr_depth = 1:n_aligned_depths
                    curr_spike_times = spike_times_timeline(aligned_str_depth_group == curr_depth);
                    
                    % (skip if no spikes at this depth)
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
                
                % Store
                str_ctxpred(curr_cohort).animal{curr_animal} = animal;
                str_ctxpred(curr_cohort).str{curr_animal}{curr_day,curr_exp} = binned_spikes;
                str_ctxpred(curr_cohort).str_ctxpred{curr_animal}{curr_day,curr_exp} = ctxpred_spikes;
                str_ctxpred(curr_cohort).ctx_str_k{curr_animal}{curr_day,curr_exp} = ctx_str_k_recast;
                
                % Reset for next experiment
                clearvars('-except',preload_vars{:})
                
            end
            AP_print_progress_fraction(curr_day,length(use_days));
        end
        
        % Save after each animal
        save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data'];
        save_fn = [save_path filesep 'str_ctxpred'];
        disp(['Saving ' save_fn]);
        save(save_fn,'str_ctxpred','-v7.3');
        disp('Saved.');
        
    end
end
disp('Done.');

%% Striatum cortical kernels (naive: master U, paired task/retinotopy)

clear all

disp('Cortex/striatum, activity/kernels (naive)');

% Parameters for regression
regression_params.use_svs = 1:100;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.1,0.1];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

% Load Master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');

n_cohorts = 1;
str_ctxpred = struct;

for curr_cohort = 1:n_cohorts
    
    switch curr_cohort
        case 1
            animals = {'AP032','AP033','AP034','AP035','AP036'};
            retinotopy_protocol = 'stimSparseNoiseUncorrAsync';
    end
        
    for curr_animal = 1:length(animals)
        
        animal = animals{curr_animal};
        disp(animal);
        
        % Get days with imaging + ephys
        retinotopy_experiments = AP_find_experiments(animal,retinotopy_protocol);
        retinotopy_experiments = retinotopy_experiments([retinotopy_experiments.imaging] & [retinotopy_experiments.ephys]);
        
        use_days = {retinotopy_experiments.day};
        % (use the last experiment of each protocol)
        use_experiments = ...
            cellfun(@(x) x(end), {retinotopy_experiments(ismember({retinotopy_experiments.day},use_days)).experiment}');
        
        for curr_day = 1:length(use_days)
            for curr_exp = 1:size(use_experiments,2)
                
                % Set variables to keep
                preload_vars = who;
                
                % Load data and align striatum by depth
                day = use_days{curr_day};
                experiment = use_experiments(curr_day,curr_exp);
                str_align = 'kernel';
                load_parts.cam = false;
                load_parts.imaging = true;
                load_parts.ephys = true;
                AP_load_experiment;
                
                % Convert U to master U
                Udf_aligned = AP_align_widefield(Udf,animal,day);
                
                % Deconvolve fluoresence (for use in regression)
                fVdf_deconv = AP_deconv_wf(fVdf);
                fVdf_deconv(isnan(fVdf_deconv)) = 0;
                
                % Get time points to bin
                sample_rate = framerate*regression_params.upsample_factor;
                time_bins = frame_t(find(frame_t > ...
                    regression_params.skip_seconds,1)):1/sample_rate: ...
                    frame_t(find(frame_t-frame_t(end) < ...
                    -regression_params.skip_seconds,1,'last'));
                time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
                
                % Resample deconvolved fluorescence
                fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
                
                % Bin spikes to match widefield frames
                binned_spikes = nan(n_aligned_depths,length(time_bin_centers));
                for curr_depth = 1:n_aligned_depths
                    curr_spike_times = spike_times_timeline(aligned_str_depth_group == curr_depth);
                    
                    % (skip if no spikes at this depth)
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
                
                % Store
                str_ctxpred(curr_cohort).animal{curr_animal} = animal;
                str_ctxpred(curr_cohort).str{curr_animal}{curr_day,curr_exp} = binned_spikes;
                str_ctxpred(curr_cohort).str_ctxpred{curr_animal}{curr_day,curr_exp} = ctxpred_spikes;
                str_ctxpred(curr_cohort).ctx_str_k{curr_animal}{curr_day,curr_exp} = ctx_str_k_recast;
                
                % Reset for next experiment
                clearvars('-except',preload_vars{:})
                
            end
            AP_print_progress_fraction(curr_day,length(use_days));
        end
        
        % Save after each animal
        save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data'];
        save_fn = [save_path filesep 'str_ctxpred_naive'];
        disp(['Saving ' save_fn]);
        save(save_fn,'str_ctxpred','-v7.3');
        disp('Saved.');
        
    end
end
disp('Done.');






