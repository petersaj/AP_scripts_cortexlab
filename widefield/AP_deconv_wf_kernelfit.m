%% Create kernel to deconvolve tetO-GCaMP6s to spikes
%
% Updated from AP_deconv_wf_kernelfit_old: 
% now new mice doing task, all recordings in AM
% NOTE: cortical ROIs and ephys alignment done in AP_ctx_str_trial_preprocessing



%% Find recordings

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


%% Initialize variables

% GCaMP kernel (to go from widefield fluorescence to spikes)
gcamp_regression_kernel = cell(size(recordings));

% Spikes kernel (to go from spikes to deconvolved widefield fluorescence)
spikes_regression_kernel = cell(size(recordings));


%% Loop through all recordings
disp('Loading and getting kernels:')

for curr_recording = 1:length(recordings)
       
    %% Clear workspace, set current recording
    clearvars -except curr_recording recordings ...
        gcamp_regression_kernel spikes_regression_kernel
    
    animal = recordings(curr_recording).animal;
    day = recordings(curr_recording).day;
    site = recordings(curr_recording).site;
    
    %% Load and concantenate data from an animal/day

    % Find all experiments for that day
    % (get folders with only a number - those're the experiment folders)
    experiments_dir = dir(AP_cortexlab_filename(animal,day,[],'expInfo'));
    experiments_num_idx = cellfun(@(x) ~isempty(x), regexp({experiments_dir.name},'^\d*$'));
    experiments = cellfun(@str2num,{experiments_dir(experiments_num_idx).name});
    
    % Loop through experiments, collate data
    skip_seconds = 60;
    
    time_bin_centers_all = cell(size(experiments));
    fVdf_resample_all = cell(size(experiments));
    binned_spikes_all = cell(size(experiments));
    
    for curr_exp = 1:length(experiments)
        
        experiment = experiments(curr_exp);
        str_align = 'none';
        AP_load_experiment;
        
        %% Set cortical units
        
        % Find cortex end by largest gap between templates
        sorted_template_depths = sort([template_depths]);
        [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
        ctx_end = sorted_template_depths(max_gap_idx);
        
        ctx_depth = [sorted_template_depths(1),ctx_end];
        ctx_units = template_depths <= ctx_depth(2);       
        
        %% Resample and concatenate data
        % Get time points to query
        sample_rate = framerate;
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
    
 
    %% Get fluorescence from pre-drawn ROI
    
    % Load ephys alignment with ROI (created in AP_ctx_str_trial_preprocessing)
    vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
    load(vis_ctx_ephys_fn);
    
    curr_animal_idx = strcmp(animal,{vis_ctx_ephys.animal});
    curr_day_idx = strcmp(day,vis_ctx_ephys(curr_animal_idx).day);
    curr_roi = vis_ctx_ephys(curr_animal_idx).ctx_roi{curr_day_idx};
         
    % Get fluorescence using average image and weights
    [fluor_trace,fluor_mask] = AP_svd_roi(Udf,fVdf_resample,[],[],curr_roi);

    %% Get kernel from wf > mua and mua > "deconv" wf

    kernel_t = [-0.5,0.5];
    kernel_frames = floor(kernel_t(1)*sample_rate):ceil(kernel_t(2)*sample_rate);  
    zs = [false,false];
    cvfold = 5;
    return_constant = true;
    use_constant = true;
    lambda = 0;

    % Regression from widefield to multiunit
    [curr_gcamp_kernel,predicted_spikes,explained_var_spikes] = ...
        AP_regresskernel(fluor_trace, ...
        binned_spikes_std,kernel_frames,lambda,zs, ...
        cvfold,return_constant,use_constant);
    
    % Regression from multiunit predicted widefield
    [curr_spikes_kernel,predicted_deconv_fluor,explained_var_fluor] = ...
        AP_regresskernel(binned_spikes_std, ...
        predicted_spikes,kernel_frames,lambda,zs, ...
        cvfold,return_constant,use_constant);
    
    % Get and store kernel in ROI
    gcamp_regression_kernel_t = kernel_frames/framerate;
    gcamp_regression_kernel{curr_recording} = curr_gcamp_kernel{1};
    
    spikes_regression_kernel_t = kernel_frames/framerate;
    spikes_regression_kernel{curr_recording} = curr_spikes_kernel{1};
   
    AP_print_progress_fraction(curr_recording,length(recordings));
    
end

%% Plot all kernels

% GCaMP regression kernel
figure;
subplot(1,2,1);
plot(gcamp_regression_kernel_t, vertcat(gcamp_regression_kernel{:})');
subplot(1,2,2);
plot(spikes_regression_kernel_t, vertcat(spikes_regression_kernel{:})');


%% Save GCaMP kernel for deconvolution

% Save the GCaMP kernels for deconvolution
gcamp6s_kernel.regression_t = gcamp_regression_kernel_t;
gcamp6s_kernel.regression = gcamp_regression_kernel;

gcamp6s_kernel.spikes_regression_t = spikes_regression_kernel_t;
gcamp6s_kernel.spikes_regression = spikes_regression_kernel;

save_dir = 'C:\Github\AP_scripts_cortexlab\widefield';
save_fn = [save_dir filesep 'gcamp6s_kernel.mat'];
save(save_fn,'gcamp6s_kernel');
disp(['Saved GCaMP6s deconvolution kernel: ' save_fn]);


