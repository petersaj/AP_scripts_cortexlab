%% AP analysis for EJ's 3-5Hz oscillation paper
% Simultaneous widefield and cortical electrophysiology
% Q: are 3-5Hz oscillations that EJ sees represented in spiking
%
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
    frame_t_all = cell(size(experiments));
    fVdf_all = cell(size(experiments));
    
    spike_times_timeline_all = cell(size(experiments));
    spikeDepths_all = cell(size(experiments));
    spike_templates_all = cell(size(experiments));
    
    disp('Loading experiments...');
    for curr_exp = 1:length(experiments)
        experiment = experiments(curr_exp);
        AP_load_experiment;
        
        frame_t_all{curr_exp} = frame_t;
        fVdf_all{curr_exp} = fVdf;
        
        curr_exp_spikes = spike_times_timeline > frame_t(1) & ...
            spike_times_timeline < frame_t(end);
        spike_times_timeline_all{curr_exp} = spike_times_timeline(curr_exp_spikes);
        spikeDepths_all{curr_exp} = spikeDepths(curr_exp_spikes);
        spike_templates_all{curr_exp} = spike_templates(curr_exp_spikes);
        
        AP_print_progress_fraction(curr_exp,length(experiments));
    end
    
    stitch_t = cumsum([0,cellfun(@(x) ...
        (x(end)-x(1))+1/framerate,frame_t_all(1:end-1))]);
    
    frame_t_stitch = cellfun(@(frame_t,stitch_t) ...
        frame_t-frame_t(1) + stitch_t, ...
        frame_t_all, num2cell(stitch_t),'uni',false);
    
    spike_times_timeline_all_stitch = cellfun(@(frame_t,spike_times_timeline,stitch_t) ...
        spike_times_timeline-frame_t(1) + stitch_t, ...
        frame_t_all,spike_times_timeline_all,num2cell(stitch_t),'uni',false);
    
    % Concatenate all data
    frame_t = cat(2,frame_t_stitch{:});
    fVdf = cat(2,fVdf_all{:});
    spike_times_timeline = cat(1,spike_times_timeline_all_stitch{:});
    spikeDepths = cat(1,spikeDepths_all{:});
    spike_templates = cat(1,spike_templates_all{:});
    
    
    %% Choose depths to pool spikes for MUA
    
    % Get MUA correlation by depth
    n_corr_groups = 40;
    depth_group_edges = linspace(0,max(channel_positions(:,2)),n_corr_groups+1);
    depth_group = discretize(templateDepths,depth_group_edges);
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
    plotSpread(template_depth_ax,templateDepths,'distributionColors','k');
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
    
    
    %% Bin spikes to get MUA
    
    % Discretize spikes into frames and count spikes per frame
    use_spikes = spikeDepths >= mua_borders(1) & spikeDepths <= mua_borders(2);
    
    frame_edges = [frame_t,frame_t(end)+1/framerate];
    frame_spikes = histcounts(spike_times_timeline(use_spikes),frame_edges);
    
    %% Get fluorescence -> MUA regression map, draw ROI for fluorescence
    
    use_svs = 1:50;
%     kernel_t = [-0.3,0.3];
%     kernel_frames = round(kernel_t(1)*framerate):round(kernel_t(2)*framerate);
    kernel_frames = -10:10;
    downsample_factor = 1;
    zs = [false,false];
    cvfold = 5;
    return_constant = false;
    use_constant = true;
    
    % Find optimal regression lambda
    disp('Finding optimal regression lambda...')
    lambda_range = [-0.5,1.5]; % ^10
    n_lambdas = 30;
    
    lambdas = logspace(lambda_range(1),lambda_range(2),n_lambdas)';
    explained_var_lambdas = nan(n_lambdas,1);
    
    figure('Name',animal); hold on;
    plot_expl_var = plot(lambdas,explained_var_lambdas,'linewidth',2);
    xlabel('\lambda');
    ylabel('Explained variance');
    drawnow;
    for curr_lambda_idx = 1:length(lambdas)
        
        curr_lambda = lambdas(curr_lambda_idx);
        
        [~,~,explained_var] = ...
        AP_regresskernel(fVdf(use_svs,:), ...
        frame_spikes./nanstd(frame_spikes),kernel_frames,curr_lambda, ...
        zs,cvfold,return_constant,use_constant);   
        
        explained_var_lambdas(curr_lambda_idx) = explained_var.total;
        
        set(plot_expl_var,'YData',explained_var_lambdas);
        drawnow
        
        AP_print_progress_fraction(curr_lambda_idx,length(lambdas));
    end
    
    lambda_bin_size = diff(lambda_range)/n_lambdas;
    explained_var_lambdas_smoothed = smooth(explained_var_lambdas,5);
    [best_lambda_explained_var,best_lambda_idx] = max(explained_var_lambdas_smoothed);
    best_lambda = lambdas(best_lambda_idx);
    lambda_range = log10([best_lambda,best_lambda]) + ...
        [-lambda_bin_size,lambda_bin_size];
    
    plot(lambdas,explained_var_lambdas_smoothed,'linewidth',2);
    line(repmat(best_lambda,1,2),ylim,'color','k');
    
    % Regress using optimal lambda
    [curr_gcamp_kernel,predicted_spikes,explained_var] = ...
        AP_regresskernel(fVdf(use_svs,:), ...
        frame_spikes./nanstd(frame_spikes),kernel_frames,best_lambda,zs, ...
        cvfold,return_constant,use_constant);   
    
    % Convert kernel to pixel space
    kernel_px = zeros(size(U,1),size(U,2),size(curr_gcamp_kernel,2),size(curr_gcamp_kernel,3),'single');
    for curr_spikes = 1:size(curr_gcamp_kernel,3)
        kernel_px(:,:,:,curr_spikes) = svdFrameReconstruct(Udf(:,:,use_svs),curr_gcamp_kernel(:,:,curr_spikes));
    end
    
    % AP_image_scroll(r_px,kernel_frames*downsample_factor/framerate);
    % caxis([-prctile(r_px(:),99.9),prctile(r_px(:),99.9)])
    % colormap(colormap_BlueWhiteRed);
    % axis image;
   
    max_weights = max(abs(kernel_px),[],3);
    
    % Draw ROI and get fluorescence using average image and weights
    [fluor_trace,fluor_mask] = AP_svd_roi(Udf,fVdf,avg_im,max_weights);
    
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
    
    % Get and store kernel in ROI
    gcamp_regression_kernel_t = kernel_frames/framerate;
    gcamp_regression_kernel{curr_recording} = ...
        AP_svd_roi(Udf(:,:,use_svs),curr_gcamp_kernel,[],[],fluor_mask);
    
    
    %% Plot MUA and fluorescence together
    
    fluor_trace_diff = diff(fluor_trace);
    
    figure('Name',animal);
    hold on
    plot(frame_t,frame_spikes/nanstd(abs(frame_spikes)),'k','linewidth',2);
    plot(frame_t,fluor_trace/nanstd(abs(fluor_trace)),'color',[0,0.5,0],'linewidth',2);
    plot(conv(frame_t,[1,1]/2,'valid'),fluor_trace_diff./nanstd(abs(fluor_trace_diff)),'color',[0,0,0.8],'linewidth',2);
    
    xlabel('Time');
    ylabel('Activity (std)');
    legend({'MUA','Fluorescence','Fluorescence derivative'});
    
    drawnow;    
    
    
    %% MUA/fluorescence cross-correlation and coherence
    
    % Cross-correlation
    corr_lag = 5; % in seconds
    corr_lag_samples = round(corr_lag*framerate);
    [mua_fluor_xcorr{curr_recording},lags] = xcorr(fluor_trace,frame_spikes,corr_lag_samples,'coeff');
    lags_t = lags./framerate;
    
    % Coherence
    [mua_fluor_coherence{curr_recording},coherence_f] = mscohere(frame_spikes',fluor_trace', ...
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
    use_trace = frame_spikes(1:end-1);
    
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
    x_nonorm = ifft((fft(fluor_trace).*conj(fft(frame_spikes)))); % unnormalized
    
    % Normalized like in Nauhaus 2012?
    soft_reg_factor = 0;
    x_autonorm = ifft((fft(fluor_trace).*conj(fft(frame_spikes)))./(soft_reg_factor+fft(frame_spikes).*conj(fft(frame_spikes))));
    
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


