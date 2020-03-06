%% Notes
%
% AP_ctx_str_ephys_alignment - beforehand: gets striatal domains to group MUA
%
% Batch scripts to save preprocessed data here, saved to:
% C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data

%% ~~~~~~~~~~~~~ Main experiment ~~~~~~~~~~~~~


%% Cortex -> kernel-aligned striatum map (protocols separately)

disp('Cortex -> striatum regression maps across protocols');

n_aligned_depths = 4;

% Parameters for regression
regression_params.use_svs = 1:100;
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


%% Choiceworld trial activity (striatum depth)

clear all
disp('Choiceworld trial activity (striatum depth)')

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
            trial_data_all         
        
    end
end

clearvars -except trial_data_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_4strdepth'];
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
        n_aligned_depths = 4;
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

%% ~~~~~~~~~~~ MUSCIMOL ~~~~~~~~~~~~~

%% Cortex -> kernel-aligned striatum map (protocols separately) (MUSCIMOL)

disp('Cortex -> striatum regression maps across protocols');

n_aligned_depths = 4;

% Parameters for regression
regression_params.use_svs = 1:100;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
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


%% Choiceworld trial activity (MUSCIMOL)
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
        exp_conditions curr_exp_condition ...
        trial_data_all           
    disp('Finished loading all')
    
    save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
    save_fn = ['trial_activity_' protocol '_' exp_conditions{curr_exp_condition}];
    save([save_path filesep save_fn],'-v7.3');
    
end



%% Passive trial activity (MUSCIMOL)
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
                exp_conditions curr_exp_condition ...
                trial_data_all
    disp('Finished loading all')
    
    save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
    save_fn = ['trial_activity_' protocol '_' exp_conditions{curr_exp_condition}];
    save([save_path filesep save_fn],'-v7.3');
    
end


%% ~~~~~~~~~~~ CORTICAL EPHYS ~~~~~~~~~~~~~

animals = {'AP043','AP060','AP061'};









