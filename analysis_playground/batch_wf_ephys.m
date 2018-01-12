%% Get cortex > spike prediction kernels across recordings

animal = 'AP028';
protocol = 'vanillaChoiceworld';
experiments = AP_find_experiments(animal,protocol);

% only use experiments with ephys + imaging
experiments = experiments([experiments.imaging] & [experiments.ephys]);

load_parts.cam = false;
load_parts.imaging = true;
load_parts.ephys = true;

batch_vars = struct;
for curr_day = 1:length(experiments);
    
    day = experiments(curr_day).day;
    experiment = experiments(curr_day).experiment;
    
    AP_load_experiment
    
    %%%%%%%%%%%%%%%
    % DO THE STUFF
    %%%%%%%%%%%%%%%    
    
    sample_rate = (1/median(diff(frame_t)))*1;
    
    % Skip the first/last n seconds to do this
    skip_seconds = 60;
    
    time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
    time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
    
    % Group multiunit by depth
    n_depth_groups = 18;
    depth_group_edges = linspace(0,max(channel_positions(:,2)),n_depth_groups+1);
    % depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
    depth_group_edges_use = depth_group_edges;
    
    [depth_group_n,depth_group] = histc(spikeDepths,depth_group_edges_use);
    depth_groups_used = unique(depth_group);
    depth_group_centers = depth_group_edges_use(1:end-1)+(diff(depth_group_edges_use)/2);
    
    binned_spikes = zeros(length(depth_group_edges_use)-1,length(time_bins)-1);
    for curr_depth = 1:length(depth_group_edges_use)-1
        
        curr_spike_times = spike_times_timeline(depth_group == curr_depth);
        %     curr_spike_times = spike_times_timeline((depth_group == curr_depth) & ...
        %         ismember(spike_templates,find(tan)));
        %     curr_spike_times = spike_times_timeline(spike_templates == 279);
        
        binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        
    end
    
    use_svs = 1:50;
    kernel_frames = -35:17;
    downsample_factor = 1;
    lambda = 2e5;
    zs = [false,true];
    cvfold = 5;
    
    fVdf_resample = interp1(frame_t,fVdf(use_svs,:)',time_bin_centers)';
    
    % TO USE fV
    % [k,predicted_spikes,explained_var] = ...
    %     AP_regresskernel(fVdf_resample, ...
    %     binned_spikes,kernel_frames,lambda,zs,cvfold);
    % TO USE dfV
    [k,predicted_spikes,explained_var] = ...
        AP_regresskernel(conv2(diff(fVdf_resample,[],2),[1,1]/2,'valid'), ...
        binned_spikes(:,2:end-1),kernel_frames,lambda,zs,cvfold);
    
    % Reshape kernel and convert to pixel space
    r = reshape(k,length(use_svs),length(kernel_frames),size(binned_spikes,1));
    
    r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
    for curr_spikes = 1:size(r,3);
        r_px(:,:,:,curr_spikes) = svdFrameReconstruct(Udf(:,:,use_svs),r(:,:,curr_spikes));
    end
    
    % Get center of mass for each pixel
    % r_px_max = squeeze(sqrt(sum(r_px.^2,3)));
    r_px_max = squeeze(max(r_px,[],3));
    r_px_max_zeronan = r_px_max;
    r_px_max_zeronan(isnan(r_px_max_zeronan)) = 0;
    r_px_max_norm = bsxfun(@rdivide,bsxfun(@minus,r_px_max_zeronan,min(r_px_max_zeronan,[],3)), ...
        max(bsxfun(@minus,r_px_max_zeronan,min(r_px_max_zeronan,[],3)),[],3));
    r_px_com = sum(bsxfun(@times,r_px_max_norm,permute(1:n_depth_groups,[1,3,2])),3)./sum(r_px_max_norm,3);
    
    r_px_weight = max(r_px_max,[],3);

    batch_vars.r_px_com{curr_day} = r_px_com;
    batch_vars.r_px_weight{curr_day} = r_px_weight;

    %%%%%%%%%%%%%%%%%%%%
    % THE STUFF IS DONE
    %%%%%%%%%%%%%%%%%%%%
    
    drawnow
    AP_print_progress_fraction(curr_day,length(experiments));
    clearvars -except experiments curr_day animal batch_vars load_parts
    
end

disp('Finished batch.')

% Align images from batch processing
batch_vars_reg = batch_vars;

% Align
days = {experiments.day};
[tform_matrix,im_aligned] = AP_align_widefield(animal,days);

for curr_day = 1:length(experiments);
    
    tform = affine2d;
    tform.T = tform_matrix{curr_day};
    
    curr_im = batch_vars_reg.r_px_com{curr_day};
    curr_im(isnan(curr_im)) = 0;    
    batch_vars_reg.r_px_com{curr_day} = imwarp(curr_im,tform, ...
        'Outputview',imref2d(size(im_aligned(:,:,1))));       
    
    curr_im = batch_vars_reg.r_px_weight{curr_day};
    curr_im(isnan(curr_im)) = 0;    
    batch_vars_reg.r_px_weight{curr_day} = imwarp(curr_im,tform, ...
        'Outputview',imref2d(size(im_aligned(:,:,1)))); 

end

% Plot map of cortical pixel by preferred depth of probe
n_depth_groups = 18;
r_px_com_col = cellfun(@(x) ind2rgb(round(mat2gray(x,[1,n_depth_groups])*255), ...
    jet(255)),batch_vars_reg.r_px_com,'uni',false);

figure;
for curr_day = 1:length(experiments)    
    curr_plot = subplot(2,ceil(length(experiments)/2),curr_day,'Visible','off');
    p = imagesc(curr_plot,r_px_com_col{curr_day});
    axis off; axis image;
    set(p,'AlphaData',mat2gray(max(batch_vars_reg.r_px_weight{curr_day},[],3), ...
        [0,double(prctile(reshape(max(batch_vars_reg.r_px_weight{curr_day},[],3),[],1),99))]));
    set(gcf,'color','w');    
end


%% Get ephys properties and visual modulation index across spikes

animal = 'AP028';
protocol = 'stimKalatsky';
experiments = AP_find_experiments(animal,protocol);

% only use experiments with ephys + imaging
experiments = experiments([experiments.imaging] & [experiments.ephys]);

load_parts.cam = false;
load_parts.imaging = false;
load_parts.ephys = true;

batch_vars = struct;
for curr_day = 1:length(experiments);
    
    day = experiments(curr_day).day;
    experiment = experiments(curr_day).experiment;
    
    AP_load_experiment
    
    %%%%%%%%%%%%%%%
    % DO THE STUFF
    %%%%%%%%%%%%%%%
    
    n_depth_groups = 18;
    
    %%% Store template information
    batch_vars.templateDepths{curr_day} = templateDepths;
   
    %%% Get rate over time
    depth_group_edges = linspace(0,max(channel_positions(:,2)),n_depth_groups+1);
    depth_group = discretize(templateDepths,depth_group_edges);
    depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);
    unique_depths = 1:length(depth_group_edges)-1;
    
    spike_binning = 60*5; % seconds
    corr_edges = spike_times_timeline(1):spike_binning:spike_times_timeline(end);
    corr_centers = corr_edges(1:end-1) + diff(corr_edges);
    
    binned_spikes_depth = zeros(length(unique_depths),length(corr_edges)-1);
    for curr_depth = 1:length(unique_depths);
        binned_spikes_depth(curr_depth,:) = histcounts(spike_times_timeline( ...
            ismember(spike_templates,find(depth_group == unique_depths(curr_depth)))), ...
            corr_edges);
    end
    
    batch_vars.spike_rate{curr_day} = binned_spikes_depth/spike_binning;
    
    %%% Get correlation of MUA and LFP   
    n_corr_groups = 40;
    depth_group_edges = linspace(0,max(channel_positions(:,2)),n_corr_groups+1);
    depth_group = discretize(templateDepths,depth_group_edges);
    depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);
    unique_depths = 1:length(depth_group_edges)-1;
    
    spike_binning = 0.01; % seconds
    corr_edges = spike_times_timeline(1):spike_binning:spike_times_timeline(end);
    corr_centers = corr_edges(1:end-1) + diff(corr_edges);
    
    binned_spikes_depth = zeros(length(unique_depths),length(corr_edges)-1);
    for curr_depth = 1:length(unique_depths);
        binned_spikes_depth(curr_depth,:) = histcounts(spike_times_timeline( ...
            ismember(spike_templates,find(depth_group == unique_depths(curr_depth)))), ...
            corr_edges);
    end
    
    mua_corr = corrcoef(binned_spikes_depth');
    
    channel_depth_grp = discretize(channel_positions(:,2),depth_group_edges);
    lfp_depth_mean = grpstats(lfp,channel_depth_grp);
    lfp_depth_mean_mediansub = bsxfun(@minus,lfp_depth_mean,nanmedian(lfp_depth_mean,1));
    lfp_corr = corrcoef(lfp_depth_mean_mediansub');
    
    batch_vars.corr_group_centers = depth_group_centers;
    batch_vars.mua_corr{curr_day} = mua_corr;
    batch_vars.lfp_corr{curr_day} = lfp_corr;
   
    %%% Get visual modulation
    align_times = stimOn_times(ismember(stimIDs,[1,2,3]));
    
    % Group by depth
    n_str_groups = 10;
    depth_group_edges = linspace(str_depth(1),str_depth(2),n_str_groups+1);
    depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
    depth_group_edges(end) = Inf;
    depth_group = discretize(spikeDepths,depth_group_edges);
    depth_groups_used = unique(depth_group);
    
    % Create MUA times grouped according to depth
    mua_times = cell(n_str_groups,1);
    for curr_depth = 1:n_str_groups
        mua_times{curr_depth} = spike_times_timeline(depth_group == curr_depth);
        %     mua_times{curr_depth} = spike_times_timeline(depth_group == curr_depth & ...
        %         ismember(spike_templates,find(msn)));
    end
    
    % PSTHs
    raster_window = [-0.2,0.2];
    psth_bin_size = 0.2;
    
    depth_psth = nan(n_str_groups,diff(raster_window)/psth_bin_size,length(unique(stimIDs)));
    unique_stim = unique(stimIDs);
    for curr_stim_idx = 1:length(unique_stim)
        curr_stim = unique_stim(curr_stim_idx);
        for curr_depth = 1:n_str_groups
            [depth_psth(curr_depth,:,curr_stim_idx),bins,rasterX,rasterY,spikeCounts] = psthAndBA( ...
                mua_times{curr_depth}, ...
                stimOn_times(ismember(stimIDs,curr_stim)), ...
                raster_window, psth_bin_size);
        end
    end
        
    baseline_rate = depth_psth(:,1,:);
    vis_rate = depth_psth(:,2,:);
    vis_modulation = permute((vis_rate-baseline_rate)./(vis_rate+baseline_rate),[1,3,2]);
    
    batch_vars.vis_modulation{curr_day} = vis_modulation;
    
    %%%%%%%%%%%%%%%%%%%%
    % THE STUFF IS DONE
    %%%%%%%%%%%%%%%%%%%%
    
    drawnow
    AP_print_progress_fraction(curr_day,length(experiments));
    clearvars -except experiments curr_day animal batch_vars load_parts
    
end

disp('Finished batch.')

% Plot basic stuff
figure;
% Templates by depth
plotSpread(batch_vars.templateDepths);set(gca,'YDir','reverse');
for curr_day = 1:length(experiments)
    line(curr_day+[-0.5,0.5],repmat(batch_vars.str_depth(curr_day,1),1,2),'color','g')
    line(curr_day+[-0.5,0.5],repmat(batch_vars.str_depth(curr_day,2),1,2),'color','g')
end
axis off; ylim([0 3840])

figure; colormap(hot);
% Spike rate
spike_rate_range = prctile(reshape(cat(2,batch_vars.spike_rate{:}),[],1),[5,95]);
for curr_day = 1:length(experiments)
    subplot(3,length(experiments),length(experiments)*0+curr_day);
    imagesc(batch_vars.spike_rate{curr_day});
    axis off;
    caxis(spike_rate_range);
end
% Ephys correlation
mua_corr_range = prctile(reshape(cat(2,batch_vars.mua_corr{:}),[],1),[5,95]);
lfp_corr_range = prctile(reshape(cat(2,batch_vars.lfp_corr{:}),[],1),[5,95]);
for curr_day = 1:length(experiments)
   subplot(3,length(experiments),length(experiments)*1+curr_day); 
   imagesc(batch_vars.corr_group_centers,batch_vars.corr_group_centers,batch_vars.mua_corr{curr_day});
   axis square off; caxis(mua_corr_range);
   line(xlim,repmat(batch_vars.str_depth(curr_day,1),1,2),'color','g')
   line(repmat(batch_vars.str_depth(curr_day,1),1,2),ylim,'color','g')
   line(xlim,repmat(batch_vars.str_depth(curr_day,2),1,2),'color','g')
   line(repmat(batch_vars.str_depth(curr_day,2),1,2),ylim,'color','g')
   
   subplot(3,length(experiments),length(experiments)*2+curr_day);
   imagesc(batch_vars.corr_group_centers,batch_vars.corr_group_centers,batch_vars.lfp_corr{curr_day});
   axis square off; caxis(lfp_corr_range);
end

% Plot average visual modulation
vis_modulation = cat(3,batch_vars.vis_modulation{:});
figure;plot(nanmean(vis_modulation,3));
legend({'Left','Center','Right'})

%% Get widefield > spike map in striatum

animal = 'AP029';
protocol = 'vanillaChoiceworld';
experiments = AP_find_experiments(animal,protocol);

% only use experiments with ephys + imaging
experiments = experiments([experiments.imaging] & [experiments.ephys]);

load_parts.cam = false;
load_parts.imaging = true;
load_parts.ephys = true;

batch_vars = struct;
for curr_day = 1:length(experiments);
    
    day = experiments(curr_day).day;
    experiment = experiments(curr_day).experiment;
    
    AP_load_experiment
    
    %%%%%%%%%%%%%%%
    % DO THE STUFF
    %%%%%%%%%%%%%%%    
       
    %%% Get the boundaries of the striatum
     n_corr_groups = 40;
    depth_group_edges = linspace(0,max(channel_positions(:,2)),n_corr_groups+1);
    depth_group = discretize(templateDepths,depth_group_edges);
    depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);
    unique_depths = 1:length(depth_group_edges)-1;
    
    spike_binning = 0.01; % seconds
    corr_edges = spike_times_timeline(1):spike_binning:spike_times_timeline(end);
    corr_centers = corr_edges(1:end-1) + diff(corr_edges);
    
    binned_spikes_depth = zeros(length(unique_depths),length(corr_edges)-1);
    for curr_depth = 1:length(unique_depths);
        binned_spikes_depth(curr_depth,:) = histcounts(spike_times_timeline( ...
            ismember(spike_templates,find(depth_group == unique_depths(curr_depth)))), ...
            corr_edges);
    end
    
    mua_corr = corrcoef(binned_spikes_depth');
    
    % start of striatum: look for ventricle (the largest gap)
    sorted_template_depths = sort([0;templateDepths]);
    [~,max_gap_idx] = max(diff(sorted_template_depths));
    str_start = sorted_template_depths(max_gap_idx+1)-1;
    % end of striatum: biggest drop in MUA correlation near end
    groups_back = 10;
    mua_corr_end = mua_corr(end-groups_back+1:end,end-groups_back+1:end);
    mua_corr_end(triu(true(length(mua_corr_end)),0)) = nan;
    median_corr = nanmedian(mua_corr_end,2);
    [x,max_corr_drop] = min(diff(median_corr));
    str_end = depth_group_centers(end-groups_back+max_corr_drop-1);    
    
    str_depth = [str_start,str_end];    
    
    %%% Get the widefield > spikes map
    sample_rate = (1/median(diff(frame_t)))*1;
    
    % Skip the first/last n seconds to do this
    skip_seconds = 60;
    
    time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
    time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
    
    % Group multiunit by depth
    n_depth_groups = 10;
    depth_group_edges = linspace(str_depth(1),str_depth(2),n_depth_groups+1);
    % depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
    depth_group_edges_use = depth_group_edges;
    
    [depth_group_n,depth_group] = histc(spikeDepths,depth_group_edges_use);
    depth_groups_used = unique(depth_group);
    depth_group_centers = depth_group_edges_use(1:end-1)+(diff(depth_group_edges_use)/2);
    
    binned_spikes = zeros(length(depth_group_edges_use)-1,length(time_bins)-1);
    for curr_depth = 1:length(depth_group_edges_use)-1
        
        curr_spike_times = spike_times_timeline(depth_group == curr_depth);
        %     curr_spike_times = spike_times_timeline((depth_group == curr_depth) & ...
        %         ismember(spike_templates,find(tan)));
        %     curr_spike_times = spike_times_timeline(spike_templates == 279);
        
        binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        
    end
    
    use_svs = 1:50;
    kernel_frames = -35:17;
    downsample_factor = 1;
    lambda = 2e5;
    zs = [false,true];
    cvfold = 5;
    
    fVdf_resample = interp1(frame_t,fVdf(use_svs,:)',time_bin_centers)';
    
    % TO USE fV
    % [k,predicted_spikes,explained_var] = ...
    %     AP_regresskernel(fVdf_resample, ...
    %     binned_spikes,kernel_frames,lambda,zs,cvfold);
    % TO USE dfV
    [k,predicted_spikes,explained_var] = ...
        AP_regresskernel(conv2(diff(fVdf_resample,[],2),[1,1]/2,'valid'), ...
        binned_spikes(:,2:end-1),kernel_frames,lambda,zs,cvfold);
    
    % Reshape kernel and convert to pixel space
    r = reshape(k,length(use_svs),length(kernel_frames),size(binned_spikes,1));
    
    r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
    for curr_spikes = 1:size(r,3);
        r_px(:,:,:,curr_spikes) = svdFrameReconstruct(Udf(:,:,use_svs),r(:,:,curr_spikes));
    end
    
    % Get center of mass for each pixel
    % r_px_max = squeeze(sqrt(sum(r_px.^2,3)));
    r_px_max = squeeze(max(r_px,[],3));
    r_px_max_zeronan = r_px_max;
    r_px_max_zeronan(isnan(r_px_max_zeronan)) = 0;
    r_px_max_norm = bsxfun(@rdivide,bsxfun(@minus,r_px_max_zeronan,min(r_px_max_zeronan,[],3)), ...
        max(bsxfun(@minus,r_px_max_zeronan,min(r_px_max_zeronan,[],3)),[],3));
    r_px_com = sum(bsxfun(@times,r_px_max_norm,permute(1:n_depth_groups,[1,3,2])),3)./sum(r_px_max_norm,3);
    
    r_px_weight = max(r_px_max,[],3);

    batch_vars.r_px_com{curr_day} = r_px_com;
    batch_vars.r_px_weight{curr_day} = r_px_weight;

    %%%%%%%%%%%%%%%%%%%%
    % THE STUFF IS DONE
    %%%%%%%%%%%%%%%%%%%%
    
    drawnow
    AP_print_progress_fraction(curr_day,length(experiments));
    clearvars -except experiments curr_day animal batch_vars load_parts
    
end

disp('Finished batch.')

% Align images from batch processing
batch_vars_reg = batch_vars;

% Align
days = {experiments.day};
[tform_matrix,im_aligned] = AP_align_widefield(animal,days);

for curr_day = 1:length(experiments);
    
    tform = affine2d;
    tform.T = tform_matrix{curr_day};
    
    curr_im = batch_vars_reg.r_px_com{curr_day};
    curr_im(isnan(curr_im)) = 0;    
    batch_vars_reg.r_px_com{curr_day} = imwarp(curr_im,tform, ...
        'Outputview',imref2d(size(im_aligned(:,:,1))));       
    
    curr_im = batch_vars_reg.r_px_weight{curr_day};
    curr_im(isnan(curr_im)) = 0;    
    batch_vars_reg.r_px_weight{curr_day} = imwarp(curr_im,tform, ...
        'Outputview',imref2d(size(im_aligned(:,:,1)))); 

end

% Plot map of cortical pixel by preferred depth of probe
n_depth_groups = 10;
r_px_com_col = cellfun(@(x) ind2rgb(round(mat2gray(x,[1,n_depth_groups])*255), ...
    jet(255)),batch_vars_reg.r_px_com,'uni',false);

figure;
for curr_day = 1:length(experiments)    
    curr_plot = subplot(2,ceil(length(experiments)/2),curr_day,'Visible','off');
    p = imagesc(curr_plot,r_px_com_col{curr_day});
    axis off; axis image;
    set(p,'AlphaData',mat2gray(max(batch_vars_reg.r_px_weight{curr_day},[],3), ...
        [0,double(prctile(reshape(max(batch_vars_reg.r_px_weight{curr_day},[],3),[],1),99))]));
    set(gcf,'color','w');    
end

%% Get striatum responses during choiceworld

animal = 'AP025';
protocol = 'vanillaChoiceworld';
experiments = AP_find_experiments(animal,protocol);

% only use experiments with ephys + imaging
experiments = experiments([experiments.imaging] & [experiments.ephys]);

load_parts.cam = false;
load_parts.imaging = false;
load_parts.ephys = true;

batch_vars = struct;
for curr_day = 1:length(experiments);
    
    day = experiments(curr_day).day;
    experiment = experiments(curr_day).experiment;
    
    AP_load_experiment
    
    %%%%%%%%%%%%%%%
    % DO THE STUFF
    %%%%%%%%%%%%%%%
    
    % Group multiunit by depth
    n_depth_groups = 6;
    depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
    depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
    
    depth_group = discretize(spikeDepths,depth_group_edges);
    
    raster_window = [-0.5,2.5];
    psth_bin_size = 0.001;
    smooth_size = 50;
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    
    psth_right_hit = nan(n_depth_groups,diff(raster_window)/psth_bin_size);
    psth_right_miss = nan(n_depth_groups,diff(raster_window)/psth_bin_size);
    psth_left_hit = nan(n_depth_groups,diff(raster_window)/psth_bin_size);
    psth_left_miss = nan(n_depth_groups,diff(raster_window)/psth_bin_size);
    
    for curr_depth = 1:n_depth_groups
        
        curr_spike_times = spike_times_timeline(depth_group == curr_depth);
        %     curr_spike_times = spike_times_timeline((depth_group == curr_depth) & ...
        %         ismember(spike_templates,find(msn)));
        
        use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 1;
        if sum(use_trials) > 0
            [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,stimOn_times(use_trials),raster_window,psth_bin_size);
            psth_smooth = conv2(psth,smWin,'same');
            psth_right_hit(curr_depth,:) = psth_smooth;
        end
        
        use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues > 0 & signals_events.missValues == 1;
        if sum(use_trials) > 0
            [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,stimOn_times(use_trials),raster_window,psth_bin_size);
            psth_smooth = conv2(psth,smWin,'same');
            psth_right_miss(curr_depth,:) = psth_smooth;
        end
        
        use_trials = signals_events.trialSideValues == -1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 1;
        if sum(use_trials) > 0
            [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,stimOn_times(use_trials),raster_window,psth_bin_size);
            psth_smooth = conv2(psth,smWin,'same');
            psth_left_hit(curr_depth,:) = psth_smooth;
        end
        
        use_trials = signals_events.trialSideValues == -1 & signals_events.trialContrastValues > 0 & signals_events.missValues == 1;
        if sum(use_trials) > 0
            [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,stimOn_times(use_trials),raster_window,psth_bin_size);
            psth_smooth = conv2(psth,smWin,'same');
            psth_left_miss(curr_depth,:) = psth_smooth;
        end
        
    end
    
    batch_vars.psth_right_hit{curr_day} = psth_right_hit;
    batch_vars.psth_right_miss{curr_day} = psth_right_miss;
    batch_vars.psth_left_hit{curr_day} = psth_left_hit;
    batch_vars.psth_left_miss{curr_day} = psth_left_miss;
    
    %%%%%%%%%%%%%%%%%%%%
    % THE STUFF IS DONE
    %%%%%%%%%%%%%%%%%%%%
    
    drawnow
    AP_print_progress_fraction(curr_day,length(experiments));
    clearvars -except experiments curr_day animal batch_vars load_parts
    
end

disp('Finished batch.')

psth_right_hit_all = nanmean(cat(3,batch_vars.psth_right_hit{:}),3);
psth_left_hit_all = nanmean(cat(3,batch_vars.psth_left_hit{:}),3);

figure; hold on;
p_rh = AP_stackplot(psth_right_hit_all',[],5,true,'k');
p_lh = AP_stackplot(psth_left_hit_all',[],5,true,'r');


%% Get striatum responses during passive

animal = 'AP025';
% protocol = 'AP_choiceWorldStimPassive';
protocol = 'stimKalatsky';
experiments = AP_find_experiments(animal,protocol);

% only use experiments with ephys + imaging
experiments = experiments([experiments.imaging] & [experiments.ephys]);

load_parts.cam = false;
load_parts.imaging = false;
load_parts.ephys = true;

batch_vars = struct;
for curr_day = 1:length(experiments);
    
    day = experiments(curr_day).day;
    experiment = experiments(curr_day).experiment;
    
    AP_load_experiment
    
    %%%%%%%%%%%%%%%
    % DO THE STUFF
    %%%%%%%%%%%%%%%
    
    % Group multiunit by depth
    n_depth_groups = 6;
    depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
    depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
    
    depth_group = discretize(spikeDepths,depth_group_edges);
    
    raster_window = [-0.5,2.5];
    psth_bin_size = 0.001;
    smooth_size = 50;
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    
    depth_group = discretize(spikeDepths,depth_group_edges);
    
    stim_psth = nan(n_depth_groups,diff(raster_window)/psth_bin_size,length(unique(stimIDs)));
    unique_stim = unique(stimIDs);
    for curr_stim_idx = 1:length(unique_stim)
        curr_stim = unique_stim(curr_stim_idx);
        for curr_depth = 1:n_depth_groups
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            
            [psth,bins,rasterX,rasterY,spikeCounts] = psthAndBA( ...
                curr_spike_times, ...
                stimOn_times(ismember(stimIDs,curr_stim)), ...
                raster_window, psth_bin_size);
            
            psth_smooth = conv2(psth,smWin,'same');
            stim_psth(curr_depth,:,curr_stim) = psth_smooth;
        end
    end
    
    batch_vars.stim_psth{curr_day} = stim_psth;
    
    %%%%%%%%%%%%%%%%%%%%
    % THE STUFF IS DONE
    %%%%%%%%%%%%%%%%%%%%
    
    drawnow
    AP_print_progress_fraction(curr_day,length(experiments));
    clearvars -except experiments curr_day animal batch_vars load_parts
    
end

stim_psth_avg = nanmean(cat(5,batch_vars.stim_psth{:}),5);

figure; hold on;
for curr_stim = 1:size(stim_psth_avg,3)
    p = AP_stackplot(stim_psth_avg(:,:,curr_stim)',[],5,true);
end

%% Meta-batch: striatum responses during choiceworld

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';
batch_vars = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    % only use experiments with ephys + imaging
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = false;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        
        %%%%%%%%%%%%%%%
        % DO THE STUFF
        %%%%%%%%%%%%%%%
        
        % Group multiunit by depth
        n_depth_groups = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
        depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
        
        depth_group = discretize(spikeDepths,depth_group_edges);
        
        raster_window = [-0.5,2.5];
        psth_bin_size = 0.001;
        smooth_size = 50;
        gw = gausswin(smooth_size,3)';
        smWin = gw./sum(gw);
        
        psth_right_hit = nan(n_depth_groups,diff(raster_window)/psth_bin_size);
        psth_right_miss = nan(n_depth_groups,diff(raster_window)/psth_bin_size);
        psth_left_hit = nan(n_depth_groups,diff(raster_window)/psth_bin_size);
        psth_left_miss = nan(n_depth_groups,diff(raster_window)/psth_bin_size);
        
        for curr_depth = 1:n_depth_groups
            
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            %     curr_spike_times = spike_times_timeline((depth_group == curr_depth) & ...
            %         ismember(spike_templates,find(msn)));
            
            use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 1;
            if sum(use_trials) > 0
                [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,stimOn_times(use_trials),raster_window,psth_bin_size);
                psth_smooth = conv2(psth,smWin,'same');
                psth_right_hit(curr_depth,:) = psth_smooth;
            end
            
            use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues > 0 & signals_events.missValues == 1;
            if sum(use_trials) > 0
                [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,stimOn_times(use_trials),raster_window,psth_bin_size);
                psth_smooth = conv2(psth,smWin,'same');
                psth_right_miss(curr_depth,:) = psth_smooth;
            end
            
            use_trials = signals_events.trialSideValues == -1 & signals_events.trialContrastValues > 0 & signals_events.hitValues == 1;
            if sum(use_trials) > 0
                [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,stimOn_times(use_trials),raster_window,psth_bin_size);
                psth_smooth = conv2(psth,smWin,'same');
                psth_left_hit(curr_depth,:) = psth_smooth;
            end
            
            use_trials = signals_events.trialSideValues == -1 & signals_events.trialContrastValues > 0 & signals_events.missValues == 1;
            if sum(use_trials) > 0
                [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,stimOn_times(use_trials),raster_window,psth_bin_size);
                psth_smooth = conv2(psth,smWin,'same');
                psth_left_miss(curr_depth,:) = psth_smooth;
            end
            
        end
        
        batch_vars(curr_animal).psth_right_hit{curr_day} = psth_right_hit;
        batch_vars(curr_animal).psth_right_miss{curr_day} = psth_right_miss;
        batch_vars(curr_animal).psth_left_hit{curr_day} = psth_left_hit;
        batch_vars(curr_animal).psth_left_miss{curr_day} = psth_left_miss;
        
        %%%%%%%%%%%%%%%%%%%%
        % THE STUFF IS DONE
        %%%%%%%%%%%%%%%%%%%%
        
        drawnow
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal])
    
end

psth_right_hit_all = nanmean(cat(3,batch_vars.psth_right_hit{:}),3);
psth_left_hit_all = nanmean(cat(3,batch_vars.psth_left_hit{:}),3);

figure; hold on;
p_rh = AP_stackplot(psth_right_hit_all',[],5,true,'k');
p_lh = AP_stackplot(psth_left_hit_all',[],5,true,'r');











