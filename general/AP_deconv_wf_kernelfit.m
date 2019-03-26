%% Create kernel to deconvolve tetO-GCaMP6s to spikes
% (adapted from AP_EJ_widefield_ephys_paper.m)
%
% Simultaneous widefield and cortical electrophysiology
% All animals are tetO-GC6S mice
% Experiments are concatenated and include: spares noise, full screen
% flickering stimuli, gray screens spontaneous, and screens off spontaneous

%% Set experiments

recordings = struct('animal',{},'day',{});

% Recording 1
recordings(1).animal = 'AP026';
recordings(1).day = '2017-12-09';

% Recording 2
recordings(2).animal = 'AP027';
recordings(2).day = '2017-12-09';

% Recording 3
recordings(3).animal = 'AP024';
recordings(3).day = '2018-07-02';

% Recording 4
recordings(4).animal = 'AP029';
recordings(4).day = '2018-08-15';

%% Initialize variables to keep

% Cross-correlation
mua_fluor_xcorr = cell(size(recordings));

% Coherence
mua_fluor_coherence = cell(size(recordings));

% Spectra correlations
spectra_corr = cell(size(recordings));

% GCaMP kernel (this is for me)
gcamp_regression_kernel = cell(size(recordings));
gcamp_fft_kernel = cell(size(recordings));

%% Loop through all recordings

for curr_recording = 1:length(recordings)
       
    %% Clear workspace, set current recording
    clearvars -except curr_recording recordings ...
        mua_fluor_xcorr mua_fluor_coherence spectra_corr ...
        gcamp_regression_kernel gcamp_fft_kernel
    
    animal = recordings(curr_recording).animal;
    day = recordings(curr_recording).day;
    
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
    
    for curr_exp = 1:length(experiments)
        
        experiment = experiments(curr_exp);
        AP_load_experiment;
        
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
            ylim([0,4000]);
            set(gca,'YDir','reverse');
            ylabel('Depth');
            title('Templates by depth');
            axis square;
            
            mua_corr_ax = subplot(1,2,2); hold on;
            imagesc(depth_group_centers,depth_group_centers,mua_corr);
            caxis([0,prctile(mua_corr(:),95)]);
            colormap(hot)
            axis image;
            ylim([0,4000]);
            xlim([0,4000]);
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
        
        %% Resample and concatenate data
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
        curr_exp_spikes = spike_depths >= mua_borders(1) & spike_depths <= mua_borders(2);        
        binned_spikes_all{curr_exp} = histcounts(spike_times_timeline(curr_exp_spikes),time_bins);   
        
    end
    
    % Concatenate all data
    time_bin_centers = cat(2,time_bin_centers_all{:});
    fVdf_resample = cat(2,fVdf_resample_all{:});
    binned_spikes = cat(2,binned_spikes_all{:});
    
    % Std-normalize spikes
    binned_spikes_std = binned_spikes/nanstd(binned_spikes);
    
 
    %% Get fluorescence -> MUA regression map, draw ROI for fluorescence
    
    use_svs = 1:50;
    kernel_t = [-0.5,0.5];
    kernel_frames = floor(kernel_t(1)*sample_rate):ceil(kernel_t(2)*sample_rate);
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
        binned_spikes_std,kernel_frames,best_lambda,zs, ...
        cvfold,return_constant,use_constant);   
    
    % Convert kernel to pixel space
    kernel_px = zeros(size(U,1),size(U,2),size(curr_gcamp_kernel{1},2),size(curr_gcamp_kernel{1},3),'single');
    for curr_spikes = 1:size(curr_gcamp_kernel{1},3)
        kernel_px(:,:,:,curr_spikes) = svdFrameReconstruct(Udf(:,:,use_svs),curr_gcamp_kernel{1}(:,:,curr_spikes));
    end
    
    % AP_image_scroll(r_px,kernel_frames*downsample_factor/framerate);
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

    %% Regress from fluorescence trace to MUA
    
    zs = [false,false];
    cvfold = 5;
    return_constant = true;
    use_constant = true;

    % Find optimal regression lambda
    lambda_range = [0,20];
    n_lambdas = 60;
    
    lambdas = linspace(lambda_range(1),lambda_range(2),n_lambdas)';
    explained_var_lambdas = nan(n_lambdas,1);
    
    figure('Name',animal); p = axes; hold on;
    plot_expl_var = plot(p,lambdas,explained_var_lambdas,'k');
    xlabel('\lambda');
    ylabel('Explained variance');
    drawnow;
    for curr_lambda_idx = 1:length(lambdas)
        
        curr_lambda = lambdas(curr_lambda_idx);
        
        [~,predicted_spikes,explained_var] = ...
            AP_regresskernel(fluor_trace, ...
            binned_spikes_std,kernel_frames,curr_lambda, ...
            zs,cvfold,return_constant,use_constant);
          
        explained_var_lambdas(curr_lambda_idx) = explained_var.total;        
        
        set(plot_expl_var,'YData',explained_var_lambdas);
        drawnow
        
    end
    
    [best_lambda_explained_var,best_lambda_idx] = max(explained_var_lambdas);
    best_lambda = lambdas(best_lambda_idx);
    line(repmat(best_lambda,1,2),ylim,'color','r','linewidth',2);
    
    % Do regression with the best lambda
    [curr_gcamp_kernel,predicted_spikes,explained_var] = ...
        AP_regresskernel(fluor_trace, ...
        binned_spikes_std,kernel_frames,best_lambda,zs, ...
        cvfold,return_constant,use_constant);
    
    % Get and store kernel in ROI
    gcamp_regression_kernel_t = kernel_frames/framerate;
    gcamp_regression_kernel{curr_recording} = curr_gcamp_kernel{1};
    
    %% Plot all traces together (std normalized)
    
    plot_t = [1:length(time_bin_centers)]/sample_rate;
 
    fluor_trace_diff = interp1(conv(plot_t,[1,1]/2,'valid'),diff(fluor_trace),plot_t);
    fVdf_deconv = convn(fVdf_resample,gcamp_regression_kernel{curr_recording},'same');
    dfVdf_roi = AP_svd_roi(Udf,fVdf_deconv,[],[],fluor_mask);
    
    figure('Name',animal);
    hold on
    plot(plot_t,binned_spikes_std,'linewidth',2);
    plot(plot_t,predicted_spikes,'linewidth',2);
    plot(plot_t,fluor_trace/nanstd(fluor_trace),'linewidth',2);
    plot(plot_t,fluor_trace_diff./nanstd(fluor_trace_diff),'linewidth',2);
    plot(plot_t,dfVdf_roi./nanstd(dfVdf_roi),'linewidth',2);
    
    xlabel('Time');
    ylabel('Activity (std)');
    legend({'MUA','Predicted MUA','Fluorescence','Fluorescence derivative','Fluorescence deconv'});
    
    drawnow;    
    
    
    %% MUA/fluorescence cross-correlation and coherence
    
    % Cross-correlation
    corr_lag = 5; % in seconds
    corr_lag_samples = round(corr_lag*framerate);
    [mua_fluor_xcorr{curr_recording},lags] = xcorr(fluor_trace,binned_spikes,corr_lag_samples,'coeff');
    lags_t = lags./framerate;
    
    % Coherence
    [mua_fluor_coherence{curr_recording},coherence_f] = mscohere(binned_spikes',fluor_trace', ...
        hamming(round(framerate*5)),round(framerate*2.5),[],framerate);
    
    % Plot
    figure('Name',animal);
    
    subplot(2,1,1);
    plot(lags_t,mua_fluor_xcorr{curr_recording},'k','linewidth',2)
    line([0,0],ylim,'linestyle','--','color','r');
    xlabel('MUA Lag (s)');
    ylabel('Normalized cross-correlation')
    
    subplot(2,1,2);
    plot(coherence_f,mua_fluor_coherence{curr_recording},'k','linewidth',2)
    xlabel('Freqency');
    ylabel('Coherence');
    
    drawnow;
    
    %% MUA/fluorescence spectral correlation
    
    % Spectrogram settings
    spect_overlap = 80;
    window_length = 3; % in seconds
    window_length_samples = round(window_length/(1/framerate));
    N = window_length_samples; % window length
    df = framerate/N; % frequency increment
    
    % Power spectrum of MUA
    use_trace = binned_spikes(1:end-1);
    
    [fluor_spect,spect_f,spect_t] = spectrogram(use_trace,window_length_samples, ...
        round(spect_overlap/100*window_length_samples),[],framerate);
    fluor_spect_norm = (fluor_spect/framerate).*conj(fluor_spect/framerate);  % framerate is used to normalize the FFT amplitudes
    fluor_spect_power = fluor_spect_norm*2*df;
    
    % Power spectrum of fluorescence
    use_trace = fluor_trace;
    
    [mua_spect,spect_f,spect_t] = spectrogram(use_trace,window_length_samples, ...
        round(spect_overlap/100*window_length_samples),[],framerate);
    mua_spect_norm = (mua_spect/framerate).*conj(mua_spect/framerate);  % framerate is used to normalize the FFT amplitudes
    mua_spect_power = mua_spect_norm*2*df;
    
    % Correlate power spectra of MUA and fluorescence
    spectra_corr_full = mat2cell(corrcoef([mua_spect_power',fluor_spect_power']),...
        repmat(length(spect_f),1,2),repmat(length(spect_f),1,2));
    spectra_corr{curr_recording} = spectra_corr_full([1,4,2]);
    
    figure('Name',animal);
    colormap(hot);
    c = [min(reshape(cell2mat(spectra_corr{curr_recording}),[],1)), ...
        max(reshape(cell2mat(spectra_corr{curr_recording}),[],1))];
    
    subplot(1,3,1);
    imagesc(spect_f,spect_f,spectra_corr{curr_recording}{1});
    xlabel('MUA frequency');
    ylabel('MUA frequency');
    axis square; caxis(c);
    
    subplot(1,3,2);
    imagesc(spect_f,spect_f,spectra_corr{curr_recording}{2});
    xlabel('Fluor frequency');
    ylabel('Fluor frequency');
    axis square; caxis(c);
    
    subplot(1,3,3);
    imagesc(spect_f,spect_f,spectra_corr{curr_recording}{3});
    xlabel('Fluor frequency');
    ylabel('MUA frequency');
    axis square; caxis(c);
    
    drawnow;
    
    %% Get the kernel from spikes to fluorescence with normalized ifft       
    
    % Non-normalized
    x_nonorm = ifft((fft(fluor_trace).*conj(fft(binned_spikes)))); % unnormalized
    
    % Normalized like in Nauhaus 2012?
    soft_reg_factor = 0;
    x_autonorm = ifft((fft(fluor_trace).*conj(fft(binned_spikes)))./(soft_reg_factor+fft(binned_spikes).*conj(fft(binned_spikes))));
    
    plot_frames = 35*5;
    
    figure;
    t_shift = [frame_t(end-plot_frames+1:end)-frame_t(end)-1/framerate,frame_t(1:plot_frames)-frame_t(1)];
    
    p1 = subplot(2,1,1); hold on;
    plot(t_shift,[x_nonorm(end-plot_frames+1:end),x_nonorm(1:plot_frames)],'k','linewidth',2);
    xlabel('Time (s)');
    ylabel('Impulse response')
    title('Non-normalized: ifft(F*S'')');
    
    p2 = subplot(2,1,2); hold on;
    plot(t_shift,[x_autonorm(end-plot_frames+1:end),x_autonorm(1:plot_frames)],'k','linewidth',2);
    xlabel('Time (s)');
    ylabel('Impluse response');
    title('Normalized: ifft(F*S''/S*S'')');
    
    linkaxes([p1,p2],'x')
    
    % Use the corrected impulse response for convolving kernel
    gcamp_fft_kernel_frames = 1:round(framerate*20);
    gcamp_fft_kernel_t = (gcamp_fft_kernel_frames-1)/framerate;
    gcamp_fft_kernel{curr_recording} = x_autonorm(gcamp_fft_kernel_frames);        
    plot(p2,gcamp_fft_kernel_t,gcamp_fft_kernel{curr_recording},'r');
    drawnow;
    
end

%% Plot summary figures

% Cross-correlation/coherence
figure;
for curr_recording = 1:length(recordings)
    
    p1 = subplot(2,1,1); hold on;
    plot(lags_t,mua_fluor_xcorr{curr_recording},'linewidth',2)
    xlabel('MUA Lag (s)');
    ylabel('Normalized cross-correlation')
    
    p2 = subplot(2,1,2); hold on;
    plot(coherence_f,mua_fluor_coherence{curr_recording},'linewidth',2)
    xlabel('Freqency');
    ylabel('Coherence');
    
end
line(p1,[0,0],ylim(p1),'linestyle','--','color','r');
legend(p1,{recordings.animal});
axis(p1,'tight');

% Spectral correlation (individual)
figure; colormap(hot);
for curr_recording = 1:length(recordings)
        
    subplot(length(recordings),3,1+3*(curr_recording-1));
    imagesc(spect_f,spect_f,spectra_corr{curr_recording}{1});
    xlabel('MUA frequency');
    ylabel('MUA frequency');
    axis square; caxis([0,1]);
    title(recordings(curr_recording).animal);
    c = colorbar;
    ylabel(c,'Correlation');
    
    subplot(length(recordings),3,2+3*(curr_recording-1));
    imagesc(spect_f,spect_f,spectra_corr{curr_recording}{2});
    xlabel('Fluor frequency');
    ylabel('Fluor frequency');
    axis square; caxis([0,1]);
    c = colorbar;
    ylabel(c,'Correlation');
    
    title(recordings(curr_recording).animal);
    subplot(length(recordings),3,3+3*(curr_recording-1));
    imagesc(spect_f,spect_f,spectra_corr{curr_recording}{3});
    xlabel('Fluor frequency');
    ylabel('MUA frequency');
    axis square; caxis([0,1]);
    title(recordings(curr_recording).animal);
    c = colorbar;
    ylabel(c,'Correlation');
end

% Spectral correlation (mean)
figure; colormap(hot);

spectra_corr_mean = arrayfun(@(curr_spect) ...
    nanmean(cell2mat(permute(cellfun(@(s) s{curr_spect}, ...
    spectra_corr,'uni',false),[1,3,2])),3),1:3,'uni',false);

subplot(1,3,1);
imagesc(spect_f,spect_f,spectra_corr_mean{1});
xlabel('MUA frequency');
ylabel('MUA frequency');
axis square; caxis([0,1]);
title('Mean');
c = colorbar;
ylabel(c,'Correlation');

subplot(1,3,2);
imagesc(spect_f,spect_f,spectra_corr_mean{2});
xlabel('Fluor frequency');
ylabel('Fluor frequency');
axis square; caxis([0,1]);
title('Mean');
c = colorbar;
ylabel(c,'Correlation');

subplot(1,3,3);
imagesc(spect_f,spect_f,spectra_corr_mean{3});
xlabel('Fluor frequency');
ylabel('MUA frequency');
axis square; caxis([0,1]);
title(recordings(curr_recording).animal);
title('Mean');
c = colorbar;
ylabel(c,'Correlation');

% GCaMP fft kernel
gcamp_kernel_norm = cellfun(@(x) (x-x(1))/max(x-x(1)),gcamp_fft_kernel,'uni',false);
figure; hold on
plot(gcamp_fft_kernel_t,vertcat(gcamp_kernel_norm{:})','linewidth',2);
legend({recordings.animal});
title('GCaMP kernel');

% GCaMP regression kernel
figure;
plot(gcamp_regression_kernel_t, vertcat(gcamp_regression_kernel{:})');

%% Save GCaMP kernel for deconvolution

% Save the GCaMP kernels for deconvolution
gcamp6s_kernel.fft_t = gcamp_fft_kernel_t;
gcamp6s_kernel.fft = gcamp_fft_kernel;

gcamp6s_kernel.regression_t = gcamp_regression_kernel_t;
gcamp6s_kernel.regression = gcamp_regression_kernel;

save('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\gcamp_kernel\gcamp6s_kernel.mat', ...
    'gcamp6s_kernel');
disp('Saved GCaMP6s deconvolution kernel');


