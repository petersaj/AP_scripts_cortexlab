%% AP analysis for EJ's 3-5Hz oscillation paper
% Simultaneous widefield and cortical electrophysiology
% Q: are 3-5Hz oscillations that EJ sees represented in spiking
%
% All animals are tetO-GC6S mice
% Experiments are concatenated and include: spares noise, full screen
% flickering stimuli, gray screens spontaneous, and screens off spontaneous


%% ~~~~ WIDEFIELD/SPIKE/LFP COHERENCE ~~~~

%% >>>> PREPROCESS AND SAVE DATA <<<<

%% Set experiments

% Initialize data structure and save filename
save_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\EJ_ctx_widefield_ephys\revision2\ctx_wf_traces.mat';
ctx_wf_traces = struct;
warning on;

% Set recordings to use
recordings = struct('animal',{},'day',{});

recordings(end+1).animal = 'AP060';
recordings(end).day = '2019-12-07';
recordings(end).genotype = 'tetO-GC6s';
recordings(end).site = 2;
recordings(end).kilosort_version = 2;

recordings(end+1).animal = 'AP043';
recordings(end).day = '2019-12-07';
recordings(end).genotype = 'tetO-GC6s';
recordings(end).site = 2;
recordings(end).kilosort_version = 2;

recordings(end+1).animal = 'AP026';
recordings(end).day = '2017-12-09';
recordings(end).genotype = 'tetO-GC6s';
recordings(end).kilosort_version = 1;

recordings(end+1).animal = 'AP027';
recordings(end).day = '2017-12-09';
recordings(end).genotype = 'tetO-GC6s';
recordings(end).kilosort_version = 1;

recordings(end+1).animal = 'AP024';
recordings(end).day = '2018-07-02';
recordings(end).genotype = 'tetO-GC6s';
recordings(end).kilosort_version = 1;

recordings(end+1).animal = 'AP029';
recordings(end).day = '2018-08-15';
recordings(end).genotype = 'tetO-GC6s';
recordings(end).kilosort_version = 1;

% (these are excluded)

% % Widefield signal weird, also most of brain isn't visible 
% recordings(end+1).animal = 'AP004';
% recordings(end).day = '2016-07-28';
% recordings(end).genotype = 'EMX-GC6s';
% 
% % Epileptiform, The regression location of this is Rsp but probe in Vis?
% recordings(end+1).animal = 'AP006';
% recordings(end).day = '2016-10-01';
% recordings(end).genotype = 'EMX-GC6f';
% 
% % Ok, but probe is in Cg and kernel is in Rsp?
% recordings(end+1).animal = 'AP007';
% recordings(end).day = '2016-12-18';
% recordings(end).genotype = 'SNAP25-GC6s';


%% Loop through all recordings

for curr_recording = 1:length(recordings)
       
    %% Clear workspace, set current recording
    clearvars -except ...
        save_fn ctx_wf_traces...
        curr_recording recordings
    
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
    upsample_factor = 1;
    skip_seconds = 60;
    
    time_bin_centers_all = cell(size(experiments));
    fVdf_resample_all = cell(size(experiments));
    binned_spikes_all = cell(size(experiments));
    lfp_all = cell(size(experiments));
    
    lfp_200Hz_all = cell(size(experiments));
    lfp_200Hz_t_all = cell(size(experiments));
    
    for curr_exp = 1:length(experiments)
        
        experiment = experiments(curr_exp);
        % Try loading with current conventions, otherwise load old
        try
            kilosort_version = recordings(curr_recording).kilosort_version;
            AP_load_experiment;
        catch me
            AP_load_old;
        end
        
        %% On first experiment: choose templates for cortex MUA
        if curr_exp == 1
            % Get MUA correlation by depth
            n_corr_groups = 40;
            depth_group_edges = linspace(0,max(channel_positions(:,2)),n_corr_groups+1);
            depth_group = discretize(template_depths,depth_group_edges);
            depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);
            unique_depths = 1:length(depth_group_edges)-1;
            
            spike_binning = 0.01; % seconds
            corr_edges = spike_times_timeline(1):spike_binning:spike_times_timeline(end);
            corr_centers = corr_edges(1:end-1) + diff(corr_edges);
            
            binned_spikes_depth = zeros(length(unique_depths),length(corr_edges)-1);
            for curr_depth = 1:length(unique_depths)
                binned_spikes_depth(curr_depth,:) = histcounts(spike_times_timeline( ...
                    ismember(spike_templates,find(depth_group == unique_depths(curr_depth)))), ...
                    corr_edges);
            end
            
            mua_corr = corrcoef(binned_spikes_depth');
            
            % Plot templates by depth and MUA correlation by depth, choose borders
            
            figure('Name',animal);
            
            template_depth_ax = subplot(1,2,1); hold on;
            plotSpread(template_depth_ax,template_depths,'distributionColors','k');
            ylim([0,max(channel_positions(:,2))]);
            set(gca,'YDir','reverse');
            ylabel('Depth');
            title('Templates by depth');
            axis square;
            
            mua_corr_ax = subplot(1,2,2); hold on;
            imagesc(depth_group_centers,depth_group_centers,mua_corr);
            caxis([0,prctile(mua_corr(:),95)]);
            colormap(hot)
            axis image;
            ylim([0,max(channel_positions(:,2))]);
            xlim([0,max(channel_positions(:,2))]);
            set(gca,'YDir','reverse');
            xlabel('Depth');
            ylabel('Depth');
            title('MUA correlation by depth')
            
            [~,mua_top] = ginput(1);
            line(template_depth_ax,xlim(template_depth_ax),repmat(mua_top,1,2),'color','b','linewidth',2);
            line(mua_corr_ax,xlim(mua_corr_ax),repmat(mua_top,1,2),'color','b','linewidth',2);
            line(mua_corr_ax,repmat(mua_top,1,2),ylim(mua_corr_ax),'color','b','linewidth',2);
            drawnow;
            
            [~,mua_bottom] = ginput(1);
            line(template_depth_ax,xlim(template_depth_ax),repmat(mua_bottom,1,2),'color','b','linewidth',2);
            line(mua_corr_ax,xlim(mua_corr_ax),repmat(mua_bottom,1,2),'color','b','linewidth',2);
            line(mua_corr_ax,repmat(mua_bottom,1,2),ylim(mua_corr_ax),'color','b','linewidth',2);
            drawnow;
            
            mua_borders = [mua_top,mua_bottom];
        end
        
        %% Load LFP channel (from middle of MUA-defined region)
        % (copied and modified from AP_load_experiment)
        
        % Find the channel closest to the middle of the defined region
        cortex_midpoint = mua_borders(1) + diff(mua_borders)/2;
        [~,lfp_channel] = min(abs(channel_positions(:,2) - cortex_midpoint));
        
        % Load LFP
        n_channels = str2num(header.n_channels);
        %lfp_filename = [ephys_path filesep 'lfp.dat']; (this is old)
        [data_path,data_path_exists] = AP_cortexlab_filename(animal,day,experiment,'ephys_dir',site);
        lfp_dir = dir([data_path filesep 'experiment*-1_0.dat']);
        if ~isempty(lfp_dir)
            lfp_filename = [data_path filesep lfp_dir.name];
        else
            % For one dataset new open ephys: just do this manually
            lfp_filename = [data_path filesep 'experiment1\recording1\continuous\Neuropix-3a-100.1\continuous.dat'];
        end
        if ~exist(lfp_filename)
            error('LFP filename incorrect')
        end
            
        lfp_sample_rate = str2num(header.lfp_sample_rate);
        lfp_cutoff = str2num(header.filter_cutoff);     
        
        % Get acqLive times for current experiment
        experiment_ephys_starts = sync(acqLive_sync_idx).timestamps(sync(acqLive_sync_idx).values == 1);
        experiment_ephys_stops = sync(acqLive_sync_idx).timestamps(sync(acqLive_sync_idx).values == 0);
        acqlive_ephys_currexpt = [experiment_ephys_starts(experiment_idx), ...
            experiment_ephys_stops(experiment_idx)];
        
        % Load single LFP channel within experiment bounds
        n_bytes = 2; % LFP = int16 = 2 bytes
        lfp_fileinfo = dir(lfp_filename); 
        n_lfp_samples = lfp_fileinfo.bytes/n_bytes/n_channels;
        lfp_memmap = memmapfile(lfp_filename, 'Format', {'int16', [n_channels n_lfp_samples], 'lfp'});
        lfp_load_start = round((lfp_sample_rate*acqlive_ephys_currexpt(1)));
        lfp_load_stop = round((lfp_sample_rate*acqlive_ephys_currexpt(2)));
        lfp = double(lfp_memmap.Data.lfp(lfp_channel,lfp_load_start:lfp_load_stop));
  
        % Get LFP times and convert to timeline time
        lfp_load_start_t = lfp_load_start/lfp_sample_rate;
        lfp_t = [0:size(lfp,2)-1]/lfp_sample_rate + lfp_load_start_t;
        lfp_t_timeline = interp1(acqlive_ephys_currexpt,acqLive_timeline,lfp_t,'linear','extrap');
        
        %%% Remove light artifact
        
        % Get light times (assume blue/violet alternate)
        light_t_timeline = interp1(sync_ephys,sync_timeline,sync(3).timestamps,'linear','extrap');
        light_on = light_t_timeline(sync(3).values == 1);
        light_off = light_t_timeline(sync(3).values == 0);
        
        light_on_median = median(light_off - light_on);
        light_off_median = median(light_on(2:end) - light_off(1:end-1));
        
        blue_on = light_on(1:2:end);
        blue_off = light_off(1:2:end);
        violet_on = light_on(2:2:end);
        violet_off = light_off(2:2:end);
        
        light_surround_t = [-(light_off_median/2):1/lfp_sample_rate:(light_on_median+(light_off_median/2))];
        
        lfp_lightfix = lfp;
        
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
        blue_on_lfp_baselinesub_med = movmedian(blue_on_lfp_baselinesub,n_light,1,'omitnan');
        violet_on_lfp_baselinesub_med = movmedian(violet_on_lfp_baselinesub,n_light,1,'omitnan');
                
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
          
                
        %% Resample and concatenate data
        % Get time points to query
        sample_rate = framerate*upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        time_bin_centers_all{curr_exp} = time_bin_centers;
        
        % Get upsampled dVdf's
        fVdf_resample_all{curr_exp} = interp1(frame_t,fVdf',time_bin_centers)';
        
        % Get binned cortex MUA  
        curr_exp_spikes = spike_depths >= mua_borders(1) & spike_depths <= mua_borders(2);        
        binned_spikes_all{curr_exp} = histcounts(spike_times_timeline(curr_exp_spikes),time_bins);   
        
        % Filter and downsample LFP (framerate)
        lowpassCutoff = framerate/2; % Hz
        [b100s, a100s] = butter(2, lowpassCutoff/(lfp_sample_rate/2), 'low');
        lfp_lightfix_lowpass = filter(b100s,a100s,lfp_lightfix);  
        lfp_lightfix_lowpass_resample = interp1(lfp_t_timeline,lfp_lightfix_lowpass,time_bin_centers);
        lfp_all{curr_exp} = lfp_lightfix_lowpass_resample;
        
        % Filter and downsample LFP (200 Hz)
        lowpassCutoff = 200/2; % Hz
        [b100s, a100s] = butter(2, lowpassCutoff/(lfp_sample_rate/2), 'low');
        lfp_lightfix_200Hz_lowpass = filter(b100s,a100s,lfp_lightfix);  
        lfp_200Hz_t = time_bin_centers(1):(1/200):time_bin_centers(end);
        lfp_lightfix_200Hz_lowpass_resample = interp1(lfp_t_timeline,lfp_lightfix_200Hz_lowpass,lfp_200Hz_t);
        lfp_200Hz_all{curr_exp} = lfp_lightfix_200Hz_lowpass_resample;
        lfp_200Hz_t_all{curr_exp} = lfp_200Hz_t;

    end
 
    
     %% Concatenate all data for analyzing
    time_bin_centers = cat(2,time_bin_centers_all{:});
    fVdf_resample = cat(2,fVdf_resample_all{:});
    binned_spikes = cat(2,binned_spikes_all{:});
    lfp = cat(2,lfp_all{:});

    
    %% Get fluorescence -> MUA regression map, draw ROI for fluorescence
    
    use_svs = 1:50;
    kernel_frames = -10:10; % specified as frames, t gave slight difference
    downsample_factor = 1;
    zs = [false,false];
    cvfold = 5;
    return_constant = true;
    use_constant = true;
    
    %%% Just set lambda, this step is now only used to pick an ROI
    %%% (uncomment if a good map is wanted)
    best_lambda = 10;
    
%     % Find optimal regression lambda
%     disp('Finding optimal regression lambda...')
%     lambda_range = [-1,1]; % ^10
%     n_lambdas = 30;
%     
%     lambdas = logspace(lambda_range(1),lambda_range(2),n_lambdas)';
%     explained_var_lambdas = nan(n_lambdas,1);
%     
%     figure('Name',animal); hold on;
%     plot_expl_var = plot(lambdas,explained_var_lambdas,'linewidth',2);
%     xlabel('\lambda');
%     ylabel('Explained variance');
%     drawnow;
%     for curr_lambda_idx = 1:length(lambdas)
%         
%         curr_lambda = lambdas(curr_lambda_idx);
%         
%         [~,~,explained_var] = ...
%         AP_regresskernel(fVdf_resample(use_svs,:), ...
%         binned_spikes./nanstd(binned_spikes),kernel_frames,curr_lambda, ...
%         zs,cvfold,return_constant,use_constant);   
%         
%         explained_var_lambdas(curr_lambda_idx) = explained_var.total;
%         
%         set(plot_expl_var,'YData',explained_var_lambdas);
%         drawnow
%         
%         AP_print_progress_fraction(curr_lambda_idx,length(lambdas));
%     end
%     
%     lambda_bin_size = diff(lambda_range)/n_lambdas;
%     explained_var_lambdas_smoothed = smooth(explained_var_lambdas,5);
%     [best_lambda_explained_var,best_lambda_idx] = max(explained_var_lambdas_smoothed);
%     best_lambda = lambdas(best_lambda_idx);
%     lambda_range = log10([best_lambda,best_lambda]) + ...
%         [-lambda_bin_size,lambda_bin_size];
%     
%     plot(lambdas,explained_var_lambdas_smoothed,'linewidth',2);
%     line(repmat(best_lambda,1,2),ylim,'color','k');
    
    % Regress using optimal lambda
    [curr_gcamp_kernel,predicted_spikes,explained_var] = ...
        AP_regresskernel(fVdf_resample(use_svs,:), ...
        binned_spikes./nanstd(binned_spikes),kernel_frames,best_lambda,zs, ...
        cvfold,return_constant,use_constant);   
    
    % Convert kernel to pixel space
    kernel_px = zeros(size(U,1),size(U,2),size(curr_gcamp_kernel{1},2),size(curr_gcamp_kernel{1},3),'single');
    for curr_spikes = 1:size(curr_gcamp_kernel{1},3)
        kernel_px(:,:,:,curr_spikes) = svdFrameReconstruct(Udf(:,:,use_svs),curr_gcamp_kernel{1}(:,:,curr_spikes));
    end
    
    % AP_imscroll(r_px,kernel_frames*downsample_factor/framerate);
    % caxis([-prctile(r_px(:),99.9),prctile(r_px(:),99.9)])
    % colormap(colormap_BlueWhiteRed);
    % axis image;
   
    max_weights = max(abs(kernel_px),[],3);
    
    % Draw ROI and get fluorescence using average image and weights
    [fluor_trace,fluor_mask] = AP_svd_roi(Udf,fVdf_resample,avg_im,max_weights);
    
    % Show the average image and regression weights with ROI drawn
    first_nonzero = find(fluor_mask > 0,1);
    [y_nonzero x_nonzero] = ind2sub([size(U,1),size(U,2)],first_nonzero);
    roi_perim = bwtraceboundary(fluor_mask,[y_nonzero x_nonzero],'N');
    
    figure('Name',animal);
    colormap(gray);
    
    subplot(1,2,1); hold on;
    imagesc(avg_im);
    axis image off;
    caxis([0,prctile(avg_im(:),99)]);
    plot(roi_perim(:,2),roi_perim(:,1),'b','linewidth',2);
    set(gca,'YDir','reverse');
    title('Average image')
    
    subplot(1,2,2); hold on;
    imagesc(max_weights);
    axis image off;
    caxis([0,prctile(max_weights(:),99)]);
    plot(roi_perim(:,2),roi_perim(:,1),'b','linewidth',2);
    set(gca,'YDir','reverse');
    title('MUA regression weights')
    
    drawnow;  
    
    %% Save data
    ctx_wf_traces(curr_recording).time_bin_centers = time_bin_centers_all;
    ctx_wf_traces(curr_recording).fluor_trace_all = ...
        mat2cell(fluor_trace,1,cellfun(@length,time_bin_centers_all));
    ctx_wf_traces(curr_recording).binned_spikes_all = binned_spikes_all;
    ctx_wf_traces(curr_recording).lfp_all = lfp_all;
    ctx_wf_traces(curr_recording).lfp_200Hz_all = lfp_200Hz_all;
    ctx_wf_traces(curr_recording).lfp_200Hz_t_all = lfp_200Hz_t_all;
    
    save(save_fn,'ctx_wf_traces')
    disp(['Saved ' save_fn])
    
end

%% >>>> LOAD, ANALYZE, PLOT <<<<

%% Load data from above
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\EJ_ctx_widefield_ephys\revision2';
data_fn = [data_path filesep 'ctx_wf_traces.mat'];
load(data_fn);

%% Initialize variables

n_recordings = length(ctx_wf_traces);

% Cross-correlation
mua_fluor_xcorr = cell(1,n_recordings);
lfp_fluor_xcorr = cell(1,n_recordings);
mua_lfp_xcorr = cell(1,n_recordings);

% Coherence
mua_fluor_coherence = cell(1,n_recordings);
lfp_fluor_coherence = cell(1,n_recordings);
mua_lfp_coherence = cell(1,n_recordings);

% Spectral correlation
mua_fluor_spectra_corr = cell(1,n_recordings);
lfp_fluor_spectra_corr = cell(1,n_recordings);
mua_lfp_spectra_corr = cell(1,n_recordings);

% Power band correlation (LFP gamma/delta, MUA 3-6 Hz, Fluor 3-6 Hz)
lfp_mua_fluor_power_corr = cell(1,n_recordings);

%% Set windowing options
framerate = 35;
window_length = 2; % in seconds
window_overlap = 1; % in seconds
window_length_samples = round(window_length/(1/framerate));
window_overlap_samples = round(window_overlap/(1/framerate));

%% Loop through all recordings and analyze

for curr_recording = 1:n_recordings
    
    % Concatenate traces from experiments
    time_bin_centers = horzcat(ctx_wf_traces(curr_recording).time_bin_centers{:});
    fluor_trace = AP_deconv_wf(horzcat(ctx_wf_traces(curr_recording).fluor_trace_all{:}));
    binned_spikes = horzcat(ctx_wf_traces(curr_recording).binned_spikes_all{:});
    lfp = horzcat(ctx_wf_traces(curr_recording).lfp_all{:});
    lfp_200Hz = horzcat(ctx_wf_traces(curr_recording).lfp_200Hz_all{:});
    
    % (median-center LFP)
    lfp = lfp - nanmean(lfp);
    lfp_200Hz = lfp_200Hz - nanmean(lfp_200Hz);
    
    %% Cross-correlation and coherence
    
    % Cross-correlation
    corr_lag = 5; % in seconds
    corr_lag_samples = round(corr_lag*framerate);
    [mua_fluor_xcorr{curr_recording},lags] = xcorr(fluor_trace,binned_spikes,corr_lag_samples,'coeff');
    [lfp_fluor_xcorr{curr_recording},lags] = xcorr(fluor_trace,lfp,corr_lag_samples,'coeff');
    [mua_lfp_xcorr{curr_recording},lags] = xcorr(lfp,binned_spikes,corr_lag_samples,'coeff');
    lags_t = lags./framerate;
    
    % Coherence
    [mua_fluor_coherence{curr_recording},coherence_f] = mscohere(binned_spikes',fluor_trace', ...
        hamming(window_length_samples),window_overlap_samples,[],framerate);
    [lfp_fluor_coherence{curr_recording},coherence_f] = mscohere(lfp',fluor_trace', ...
        hamming(window_length_samples),window_overlap_samples,[],framerate);
    [mua_lfp_coherence{curr_recording},coherence_f] = mscohere(binned_spikes',lfp', ...
        hamming(window_length_samples),window_overlap_samples,[],framerate);
    
    %% Spectral correlation
        
    % Power spectra
    [mua_spect,spect_f,spect_t,mua_spect_power] = spectrogram(binned_spikes,window_length_samples, ...
        window_overlap_samples,[],framerate);    
    [lfp_spect,spect_f,spect_t,lfp_spect_power] = spectrogram(lfp,window_length_samples, ...
        window_overlap_samples,[],framerate);  
    [fluor_spect,spect_f,spect_t,fluor_spect_power] = spectrogram(fluor_trace,window_length_samples, ...
        window_overlap_samples,[],framerate);

%     % (to plot the power spectra)
%     figure; 
%     p1 = subplot(1,3,1);
%     imagesc(spect_t,spect_f,mua_spect_power);
%     colormap(hot); set(gca,'YDir','normal');
%     p2 = subplot(1,3,2);
%     imagesc(spect_t,spect_f,lfp_spect_power);
%     colormap(hot); set(gca,'YDir','normal');
%     p3 = subplot(1,3,3);
%     imagesc(spect_t,spect_f,fluor_spect_power);
%     colormap(hot); set(gca,'YDir','normal');
%     linkaxes([p1,p2,p3]);
    
    % Set time bins to use (if splitting high/low firing rate)
    spect_t_binsize = mean(diff(spect_t));
    spect_t_bins = [spect_t - spect_t_binsize/2,spect_t(end) + spect_t_binsize/2];
    spect_t_mua_bins = discretize([1:length(binned_spikes)]./framerate,spect_t_bins);
    spect_t_binned_spikes = accumarray(spect_t_mua_bins(~isnan(spect_t_mua_bins))', ...
        binned_spikes(~isnan(spect_t_mua_bins))')';
    
    highpassCutoff = 0.5; % Hz
    [b100s, a100s] = butter(2, highpassCutoff/(framerate/2), 'high');
    spect_t_binned_spikes_highpass = filter(b100s,a100s,spect_t_binned_spikes')';
    
    median_spikes = median(spect_t_binned_spikes_highpass);
    % (NOTE: SET THIS MANUALLY)
%     warning('Using <= median spikes time points');
%     use_t = spect_t_binned_spikes_highpass <= median_spikes;
%     warning('Using > median spikes time points')
%     use_t = spect_t_binned_spikes_highpass > median_spikes;
    warning('Using all time points')
    use_t = true(size(spect_t));
    
    % Spectral correlation
    mua_fluor_spectra_corr_full = mat2cell(corrcoef([mua_spect_power(:,use_t)',fluor_spect_power(:,use_t)']),...
        repmat(length(spect_f),1,2),repmat(length(spect_f),1,2));
    mua_fluor_spectra_corr{curr_recording} = mua_fluor_spectra_corr_full([1,4,2]);
    
    lfp_fluor_spectra_corr_full = mat2cell(corrcoef([lfp_spect_power(:,use_t)',fluor_spect_power(:,use_t)']),...
        repmat(length(spect_f),1,2),repmat(length(spect_f),1,2));
    lfp_fluor_spectra_corr{curr_recording} = lfp_fluor_spectra_corr_full([1,4,2]);
    
    mua_lfp_spectra_corr_full = mat2cell(corrcoef([mua_spect_power(:,use_t)',lfp_spect_power(:,use_t)']),...
        repmat(length(spect_f),1,2),repmat(length(spect_f),1,2));
    mua_lfp_spectra_corr{curr_recording} = mua_lfp_spectra_corr_full([1,4,2]);
    
    
    %%%%%%%%% NEW THING:
    % correlate gamma/delta LFP with 3-6 Hz in widefield
    
    % Power spectrum of 200 Hz LFP
    window_length_samples_200Hz = round(window_length/(1/200));
    window_overlap_samples_200Hz = round(window_overlap/(1/200));
    
    [lfp200Hz_spect,spect_f_200Hz,spect_t_200Hz,lfp200Hz_spect_power] = ...
        spectrogram(lfp_200Hz,window_length_samples_200Hz, ...
        window_overlap_samples_200Hz,[],200);
    
    delta_f = spect_f >= 1 & spect_f <= 4;
    rest_f = spect_f >= 3 & spect_f <= 6;
    gamma_f = spect_f >= 20 & spect_f <= 80;
    
    delta_f_200Hz = spect_f_200Hz >= 1 & spect_f_200Hz <= 4;
    rest_f_200Hz = spect_f_200Hz >= 3 & spect_f_200Hz <= 6;
    gamma_f_200Hz = spect_f_200Hz >= 20 & spect_f_200Hz <= 80;
    
    lfp200Hz_gamma_delta_ratio = ...
        nanmean(lfp200Hz_spect_power(gamma_f_200Hz,:),1)./ ...
        nanmean(lfp200Hz_spect_power(delta_f_200Hz,:),1);
    
    mua_3_6_power = nanmean(mua_spect_power(rest_f,:),1);
    fluor_3_6_power = nanmean(fluor_spect_power(rest_f,:),1);
    
    lfp_mua_fluor_power_corr{curr_recording} = ...
        corrcoef([lfp200Hz_gamma_delta_ratio',mua_3_6_power',fluor_3_6_power']);
    
end

%% Plot figures

% Cross-correlation/coherence
figure;
for curr_recording = 1:n_recordings
    
    p1 = subplot(3,2,1); hold on;
    plot(lags_t,mua_fluor_xcorr{curr_recording},'linewidth',1)
    xlabel('MUA Lag (s)');
    ylabel('MUA-fluor norm. cross-correlation')
    
    p2 = subplot(3,2,2); hold on;
    plot(coherence_f,mua_fluor_coherence{curr_recording},'linewidth',1)
    xlabel('Frequency');
    ylabel('MUA-fluor Coherence');
    
    p3 = subplot(3,2,3); hold on;
    plot(lags_t,lfp_fluor_xcorr{curr_recording},'linewidth',1)
    xlabel('LFP Lag (s)');
    ylabel('LFP-fluor norm. cross-correlation')
    
    p4 = subplot(3,2,4); hold on;
    plot(coherence_f,lfp_fluor_coherence{curr_recording},'linewidth',1)
    xlabel('Frequency');
    ylabel('LFP-Fluor Coherence');
    
    p5 = subplot(3,2,5); hold on;
    plot(lags_t,mua_lfp_xcorr{curr_recording},'linewidth',1)
    xlabel('MUA Lag (s)');
    ylabel('MUA-LFP norm. cross-correlation')
    
    p6 = subplot(3,2,6); hold on;
    plot(coherence_f,mua_lfp_coherence{curr_recording},'linewidth',1)
    xlabel('Frequency');
    ylabel('MUA-LFP Coherence');
    
end

subplot(3,2,1); 
mua_fluor_xcorr_mean = nanmean(vertcat(mua_fluor_xcorr{:}),1);
plot(lags_t,mua_fluor_xcorr_mean,'k','linewidth',2);

subplot(3,2,2); 
mua_fluor_coherence_mean = nanmean(horzcat(mua_fluor_coherence{:}),2);
plot(coherence_f,mua_fluor_coherence_mean,'k','linewidth',2);

subplot(3,2,3); 
lfp_fluor_xcorr_mean = nanmean(vertcat(lfp_fluor_xcorr{:}),1);
plot(lags_t,lfp_fluor_xcorr_mean,'k','linewidth',2);

subplot(3,2,4); 
lfp_fluor_coherence_mean = nanmean(horzcat(lfp_fluor_coherence{:}),2);
plot(coherence_f,lfp_fluor_coherence_mean,'k','linewidth',2);

subplot(3,2,5); 
mua_lfp_xcorr_mean = nanmean(vertcat(mua_lfp_xcorr{:}),1);
plot(lags_t,mua_lfp_xcorr_mean,'k','linewidth',2);

subplot(3,2,6); 
mua_lfp_coherence_mean = nanmean(horzcat(mua_lfp_coherence{:}),2);
plot(coherence_f,mua_lfp_coherence_mean,'k','linewidth',2);

line(p1,[0,0],ylim(p1),'linestyle','--','color','r');
line(p3,[0,0],ylim(p3),'linestyle','--','color','r');
line(p5,[0,0],ylim(p5),'linestyle','--','color','r');

legend(p1,[repmat({'Mouse'},1,n_recordings),'Average']);
axis(p1,'tight');
axis(p3,'tight');
axis(p5,'tight');

% Spectral correlation 
figure; colormap(hot);

mua_fluor_spectra_corr_mean = arrayfun(@(curr_spect) ...
    nanmean(cell2mat(permute(cellfun(@(s) s{curr_spect}, ...
    mua_fluor_spectra_corr,'uni',false),[1,3,2])),3),1:3,'uni',false);

lfp_fluor_spectra_corr_mean = arrayfun(@(curr_spect) ...
    nanmean(cell2mat(permute(cellfun(@(s) s{curr_spect}, ...
    lfp_fluor_spectra_corr,'uni',false),[1,3,2])),3),1:3,'uni',false);

mua_lfp_spectra_corr_mean = arrayfun(@(curr_spect) ...
    nanmean(cell2mat(permute(cellfun(@(s) s{curr_spect}, ...
    mua_lfp_spectra_corr,'uni',false),[1,3,2])),3),1:3,'uni',false);

subplot(3,3,1);
imagesc(spect_f,spect_f,mua_fluor_spectra_corr_mean{1});
xlabel('MUA frequency');
ylabel('MUA frequency');
axis square; caxis([0,1]);
title('Mean');
c = colorbar;
ylabel(c,'Correlation');

subplot(3,3,2);
imagesc(spect_f,spect_f,mua_fluor_spectra_corr_mean{2});
xlabel('Fluor frequency');
ylabel('Fluor frequency');
axis square; caxis([0,1]);
title('Mean');
c = colorbar;
ylabel(c,'Correlation');

subplot(3,3,3);
imagesc(spect_f,spect_f,mua_fluor_spectra_corr_mean{3});
xlabel('MUA frequency');
ylabel('Fluor frequency');
axis square; caxis([0,1]);
title('Mean');
c = colorbar;
ylabel(c,'Correlation');

subplot(3,3,4);
imagesc(spect_f,spect_f,lfp_fluor_spectra_corr_mean{1});
xlabel('LFP frequency');
ylabel('LFP frequency');
axis square; caxis([0,1]);
title('Mean');
c = colorbar;
ylabel(c,'Correlation');

subplot(3,3,5);
imagesc(spect_f,spect_f,lfp_fluor_spectra_corr_mean{2});
xlabel('Fluor frequency');
ylabel('Fluor frequency');
axis square; caxis([0,1]);
title('Mean');
c = colorbar;
ylabel(c,'Correlation');

subplot(3,3,6);
imagesc(spect_f,spect_f,lfp_fluor_spectra_corr_mean{3});
xlabel('LFP frequency');
ylabel('Fluor frequency');
axis square; caxis([0,1]);
title('Mean');
c = colorbar;
ylabel(c,'Correlation');

subplot(3,3,7);
imagesc(spect_f,spect_f,mua_lfp_spectra_corr_mean{1});
xlabel('MUA frequency');
ylabel('MUA frequency');
axis square; caxis([0,1]);
title('Mean');
c = colorbar;
ylabel(c,'Correlation');

subplot(3,3,8);
imagesc(spect_f,spect_f,mua_lfp_spectra_corr_mean{2});
xlabel('LFP frequency');
ylabel('LFP frequency');
axis square; caxis([0,1]);
title('Mean');
c = colorbar;
ylabel(c,'Correlation');

subplot(3,3,9);
imagesc(spect_f,spect_f,mua_lfp_spectra_corr_mean{3});
xlabel('MUA frequency');
ylabel('LFP frequency');
axis square; caxis([0,1]);
title('Mean');
c = colorbar;
ylabel(c,'Correlation');

% Power band correlation 
% (order of variables: LFP gamma/delta, MUA 3-6Hz, Fluor 3-6Hz)
lfp_mua_fluor_power_corr_cat = cat(3,lfp_mua_fluor_power_corr{:});
figure;
subplot(1,2,1);
imagesc(nanmean(lfp_mua_fluor_power_corr_cat,3));
caxis([-1,1]); axis image;
set(gca,'YTick',1:3,'YTickLabel',{'LFP \gamma/\delta','MUA 3-6Hz','LFP 3-6Hz'});
set(gca,'XTick',1:3,'XTickLabel',{'LFP \gamma/\delta','MUA 3-6Hz','LFP 3-6Hz'});
c = colorbar; 
ylabel(c,'Correlation');
colormap(brewermap([],'*RdBu'));
subplot(1,2,2); hold on;
plot_corr = [squeeze(lfp_mua_fluor_power_corr_cat(1,2,:)), ...
    squeeze(lfp_mua_fluor_power_corr_cat(1,3,:)), ...
    squeeze(lfp_mua_fluor_power_corr_cat(2,3,:))];
plot(plot_corr','color',[0.5,0.5,0.5]);
plot(nanmean(plot_corr,1),'linewidth',2,'color','k');
set(gca,'XTick',1:3,'XTickLabel',{'LFP \gamma/\delta:MUA 3-6Hz', ...
    'LFP \gamma/\delta:Fluor 3-6Hz','MUA 3-6Hz:Fluor 3-6Hz'});
ylabel('Correlation');
xlim([0.5,3.5]);
line(xlim,[0,0],'linestyle','--','color','r');

%% ~~~~ MUSCIMOL ~~~~

%% >>>> PREPROCESS AND SAVE DATA <<<<

% Initialize data structure and save filename
save_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\EJ_ctx_widefield_ephys\revision2\muscimol_traces.mat';
muscimol_traces = struct;

% Set recordings to use
recordings = struct('animal',{},'day',{});

recordings(end+1).animal = 'AP040';
recordings(end).day = '2019-09-12';

recordings(end+1).animal = 'AP041';
recordings(end).day = '2019-09-12';

recordings(end+1).animal = 'AP054';
recordings(end).day = '2019-10-27';

recordings(end+1).animal = 'AP055';
recordings(end).day = '2019-10-27';

recordings(end+1).animal = 'AP045';
recordings(end).day = '2019-11-03';

recordings(end+1).animal = 'AP053';
recordings(end).day = '2019-11-22';

recordings(end+1).animal = 'AP047';
recordings(end).day = '2019-11-29';

recordings(end+1).animal = 'AP048';
recordings(end).day = '2019-12-02';

for curr_recording = 1:length(recordings)
    
    animal = recordings(curr_recording).animal;
    day = recordings(curr_recording).day;
    
    experiments = AP_list_experiments(animal,day);
    use_experiments = [experiments(strcmp({experiments.protocol},'AP_sparseNoise')).experiment];
    if length(use_experiments) ~= 2
        error([animal ' ' day ' ~= 2 experiments'])
    end
    
    for curr_experiment_idx = 1:length(use_experiments)
        
        % Load experiment
        experiment = use_experiments(curr_experiment_idx);
        load_parts.imaging = true;
        AP_load_experiment;
        
        avg_im_aligned = AP_align_widefield(avg_im,animal,day);
        Udf_aligned = AP_align_widefield(Udf,animal,day);
        
        % Draw ROI over craniotomy on first experiment
        if curr_experiment_idx == 1
            [muscimol_trace,muscimol_mask] = AP_svd_roi(Udf_aligned,fVdf,avg_im_aligned);
            intact_mask = fliplr(muscimol_mask);
            intact_trace = AP_svd_roi(Udf_aligned,fVdf,[],[],intact_mask);
        else
            muscimol_trace = AP_svd_roi(Udf_aligned,fVdf,[],[],muscimol_mask);
            intact_trace = AP_svd_roi(Udf_aligned,fVdf,[],[],intact_mask);
        end
               
        % Save
        if curr_experiment_idx == 1
            muscimol_traces(curr_recording).muscimolhemi_premuscimol = muscimol_trace;
            muscimol_traces(curr_recording).intacthemi_premuscimol = intact_trace;
        elseif curr_experiment_idx == 2
            muscimol_traces(curr_recording).muscimolhemi_postmuscimol = muscimol_trace;
            muscimol_traces(curr_recording).intacthemi_postmuscimol = intact_trace;
        end
        
    end
 
end

save(save_fn,'muscimol_traces');
disp(['Saved ' save_fn]);


%% >>>> LOAD, ANALYZE, PLOT <<<<

% Load data
muscimol_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\EJ_ctx_widefield_ephys\revision2\muscimol_traces.mat';
load(muscimol_fn);
framerate = 35;

muscimol_power_pre = cell(length(muscimol_traces),1);
intact_power_pre = cell(length(muscimol_traces),1);
muscimol_power_post = cell(length(muscimol_traces),1);
intact_power_post = cell(length(muscimol_traces),1);

figure; 
for curr_recording = 1:length(muscimol_traces)
    
    muscimol_trace_pre = muscimol_traces(curr_recording).muscimolhemi_premuscimol;
    intact_trace_pre = muscimol_traces(curr_recording).intacthemi_premuscimol;
    
    muscimol_trace_post = muscimol_traces(curr_recording).muscimolhemi_postmuscimol;
    intact_trace_post = muscimol_traces(curr_recording).intacthemi_postmuscimol;
    
    Fs = framerate;
    L = length(muscimol_trace_pre);
    freqs = 0.1:0.2:15;
    [P_muscimol_pre,F] = pwelch(double(muscimol_trace_pre)',[],[],freqs,Fs);
    [P_intact_pre,F] = pwelch(double(intact_trace_pre)',[],[],freqs,Fs);
    [P_muscimol_post,F] = pwelch(double(muscimol_trace_post)',[],[],freqs,Fs);
    [P_intact_post,F] = pwelch(double(intact_trace_post)',[],[],freqs,Fs);
    
    muscimol_power_pre{curr_recording} = P_muscimol_pre;
    intact_power_pre{curr_recording} = P_intact_pre;
    muscimol_power_post{curr_recording} = P_muscimol_post;
    intact_power_post{curr_recording} = P_intact_post;
    
    subplot(length(muscimol_traces),1,curr_recording); hold on
    plot_t = (0:length(muscimol_trace_post)-1)/framerate;
    plot(plot_t,muscimol_trace_post);
    plot(plot_t,intact_trace_post);
    xlabel('Time (s)');
    ylabel('\DeltaF/F');
    if curr_recording == 1
       legend({'Muscimol side','Intact side'}); 
    end

end

figure; 
p1 = subplot(2,2,1); hold on;
plot(F,log10(vertcat(intact_power_pre{:}))','color',[0.5,0.5,0.5]);
plot(F,log10(vertcat(muscimol_power_pre{:}))','color',[1,0.5,0.5]);
plot(F,log10(nanmean(vertcat(intact_power_pre{:})',2)),'k','linewidth',2);
plot(F,log10(nanmean(vertcat(muscimol_power_pre{:})',2)),'r','linewidth',2);
xlabel('Frequency');
ylabel('Log_10 power');
title('Pre-muscimol');

p2 = subplot(2,2,2); hold on;
plot(F,log10(vertcat(intact_power_post{:}))','color',[0.5,0.5,0.5]);
plot(F,log10(vertcat(muscimol_power_post{:}))','color',[1,0.5,0.5]);
plot(F,log10(nanmean(vertcat(intact_power_post{:})',2)),'k','linewidth',2);
plot(F,log10(nanmean(vertcat(muscimol_power_post{:})',2)),'r','linewidth',2);
xlabel('Frequency');
ylabel('Log_10 power');
title('Post-muscimol');

p3 = subplot(2,2,3); hold on;
plot(F,[vertcat(muscimol_power_pre{:})./vertcat(intact_power_pre{:})],'color',[0.5,0.5,0.5]);
plot(F,nanmean([vertcat(muscimol_power_pre{:})./vertcat(intact_power_pre{:})]',2),'k','linewidth',2);
line(xlim,[1,1],'color','r','linestyle','--');
xlabel('Frequency');
ylabel('Muscimol/intact ratio')
title('Pre-muscimol');

p4 = subplot(2,2,4); hold on;
plot(F,[vertcat(muscimol_power_post{:})./vertcat(intact_power_post{:})]','color',[0.5,0.5,0.5]);
plot(F,nanmean([vertcat(muscimol_power_post{:})./vertcat(intact_power_post{:})]',2),'k','linewidth',2);
line(xlim,[1,1],'color','r','linestyle','--');
xlabel('Frequency');
ylabel('Muscimol/intact ratio')
title('Post-muscimol');

linkaxes([p1,p2]);
linkaxes([p3,p4]);











