%% Get widefield area boundaries in batch

clear all
animal = 'AP027';
protocol = 'vanillaChoiceworld';
experiments = AP_find_experiments(animal,protocol);
experiments = experiments([experiments.imaging]);

load_parts.cam = false;
load_parts.imaging = true;
load_parts.ephys = false;

batch_vars = struct;
for curr_day = 1:length(experiments);
    
    day = experiments(curr_day).day;
    experiment = experiments(curr_day).experiment;
    
    AP_load_experiment
    
    %%%%%%%%%%%%%%%
    % DO THE STUFF
    %%%%%%%%%%%%%%%
    
    % Get V covariance
    Ur = reshape(U, size(U,1)*size(U,2),[]); % P x S
    covV = cov(fV'); % S x S % this is the only one that takes some time really
    varP = dot((Ur*covV)', Ur'); % 1 x P
    
    ySize = size(U,1); xSize = size(U,2);
    
    px_spacing = 20;
    use_y = 1:px_spacing:size(U,1);
    use_x = 1:px_spacing:size(U,2);
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
       
    % Correlation map edge detection
    corr_map_norm = mat2gray(cat(3,corr_map{:}),[0.5,1]);
    corr_map_edge = imgaussfilt(corr_map_norm,5)-imgaussfilt(corr_map_norm,20);
    corr_edges = nanmean(corr_map_edge,3);
        
    batch_vars.corr_edges{curr_day} = corr_edges;
    
    %%%%%%%%%%%%%%%%%%%%
    % THE STUFF IS DONE
    %%%%%%%%%%%%%%%%%%%%
    
    drawnow
    AP_print_progress_fraction(curr_day,length(experiments))
    clearvars -except experiments curr_day animal batch_vars load_parts
    
end

disp('Finished batch.')

% Align images from batch processing
days = {experiments.day};
corr_edges_aligned = AP_align_widefield(animal,days,batch_vars.corr_edges);

wf_borders_fig = figure('Name',animal);
imagesc(nanmean(corr_edges_aligned,3));
axis image off; colormap(gray); caxis([0,0.05])

fn = ['\\basket.cortexlab.net\data\ajpeters\wf_borders' filesep animal '_wf_borders'];
saveas(wf_borders_fig,fn);

%% Batch load responses to passive stim

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
% protocol = 'stimKalatsky';
protocol = 'AP_choiceWorldStimPassive';

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    % Move on if this animal doesn't have this experiment
    if isempty(experiments)
        continue
    end
    
    disp(animal);
    
    experiments = experiments([experiments.imaging]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = false;
    
    batch_vars = struct;
    for curr_day = 1:length(experiments);
        
        % Load the experiment
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;        
        AP_load_experiment
             
        % Set options
        surround_window = [-0.5,5];
        framerate = 1./median(diff(frame_t));
        surround_samplerate = 1/(framerate*1);
        t_surround = surround_window(1):surround_samplerate:surround_window(2);
        
        % Average (time course) responses
        conditions = unique(stimIDs);
        im_stim = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        for curr_condition_idx = 1:length(conditions)
            curr_condition = conditions(curr_condition_idx);
            
            use_stims = find(stimIDs == curr_condition);
            use_stim_onsets = stimOn_times(use_stims(2:end));
            use_stim_onsets([1,end]) = [];
            
            stim_surround_times = bsxfun(@plus, use_stim_onsets(:), t_surround);
            peri_stim_v = permute(mean(interp1(frame_t,fVdf',stim_surround_times),1),[3,2,1]);
            
            im_stim(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
        end
        
        batch_vars.im_stim{curr_day} = im_stim;        
        
        % Prepare for next loop
        AP_print_progress_fraction(curr_day,length(experiments))
        clearvars -except animals protocol curr_animal experiments curr_day animal batch_vars load_parts
        
    end
    
    %     % Get ddf
    %     ddf_im = batch_vars.im_stim;
    %     for curr_day = 1:length(experiments)
    %         curr_im = ddf_im{curr_day};
    %         curr_im(isnan(curr_im)) = 0;
    %         curr_im = imgaussfilt(diff(curr_im,[],3),2);
    %         curr_im(curr_im < 0) = 0;
    %         ddf_im{curr_day} = curr_im;
    %     end
    
    % Align
    days = {experiments.day};
    im_aligned = AP_align_widefield(animal,days,batch_vars.im_stim);
    im_aligned_average = nanmean(im_aligned,5);
    
    %     % Plot
    %     surround_window = [-0.5,5];
    %     t_surround = linspace(surround_window(1),surround_window(2),size(im_aligned_average,3));
    %     AP_image_scroll(im_aligned_average,t_surround);
    %     colormap(gray); axis image off;
    %     title([animal ': passive stimuli']);
    
    % Save
    save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\passive'];
    save([save_path filesep animal '_' protocol],'im_aligned_average');    
    
    disp(['Finished ' animal]);
    
end

disp('Finished batch');

%% Batch get average choiceworld fluorescence

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);    
    
    disp(animal);
    
    experiments = experiments([experiments.imaging]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = false;
    
    batch_vars = struct;
    % (initialize these because iteratively added below)
    batch_vars.im_stim_miss = [];
    batch_vars.im_stim_hit = [];    
    batch_vars.n_im_stim_miss = [];
    batch_vars.n_im_stim_hit = [];
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        
        use_stim = true(1,min(length(signals_events.trialSideValues),length(stimOn_times)));
        stimIDs = signals_events.trialSideValues(use_stim).*signals_events.trialContrastValues(use_stim);
        stim_onsets = stimOn_times(use_stim);
        
        %         % (to only use certain stimIDs)
        %         use_stim = ismember(stimIDs,[-1,-0.125,0.125,1]);
        %         stimIDs = stimIDs(use_stim);
        %         stim_onsets = stim_onsets(use_stim);
        %         % (to discretize the stimIDs by easy/hard/zero)
        %         stimIDs = discretize(stimIDs,[-Inf,-0.125,-0.01,0.01,0.25,Inf],[-2,-1,0,1,2]);
        
        %%%% Get wheel move time
        t_surround = [-0.5,5];
        surround_samples = t_surround/Timeline.hw.samplingInterval;
        
        % Get wheel aligned to stim onset
        rotaryEncoder_idx = strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder');
        t_surround = t_surround(1):Timeline.hw.samplingInterval:t_surround(2);
        pull_times = bsxfun(@plus,stim_onsets,t_surround);
        
        stim_aligned_wheel_raw = interp1(Timeline.rawDAQTimestamps, ...
            Timeline.rawDAQData(:,rotaryEncoder_idx),pull_times);
        stim_aligned_wheel = bsxfun(@minus,stim_aligned_wheel_raw, ...
            nanmedian(stim_aligned_wheel_raw(:,t_surround < 0),2));
        
        % Define time to first wheel movement
        thresh_displacement = 2;
        [~,wheel_move_sample] = max(abs(stim_aligned_wheel) > thresh_displacement,[],2);
        wheel_move_time = arrayfun(@(x) pull_times(x,wheel_move_sample(x)),1:size(pull_times,1));
        wheel_move_time(wheel_move_sample == 1) = NaN;
        %%%%
        
        % Set options
        surround_window = [-0.5,5];
        baseline_surround_window = [0,0];
        framerate = 1./median(diff(frame_t));
        surround_samplerate = 1/(framerate*1);
        t_surround = surround_window(1):surround_samplerate:surround_window(2);
        baseline_surround_time = baseline_surround_window(1):surround_samplerate:baseline_surround_window(2);
        
        % Average (time course) responses
        conditions = unique(stimIDs);
        im_stim_hit = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        im_stim_miss = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        
        for curr_condition_idx = 1:length(conditions)
            curr_condition = conditions(curr_condition_idx);
            
            use_stims = find(stimIDs == curr_condition & signals_events.hitValues(use_stim) == 0);
            use_stim_onsets = stim_onsets(use_stims);
            if length(use_stim_onsets) > 5
                stim_surround_times = bsxfun(@plus, use_stim_onsets(:), t_surround);
                peri_stim_v = permute(mean(interp1(frame_t,fVdf',stim_surround_times),1),[3,2,1]);
                im_stim_miss(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
            end
            
            use_stims = find(stimIDs == curr_condition & signals_events.hitValues(use_stim) == 1);
            use_stim_onsets = stim_onsets(use_stims);
            if length(use_stim_onsets) > 5
                stim_surround_times = bsxfun(@plus, use_stim_onsets(:), t_surround);
                peri_stim_v = permute(mean(interp1(frame_t,fVdf',stim_surround_times),1),[3,2,1]);
                im_stim_hit(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
            end
            
        end
        
        % Align and average as it goes, otherwise the variable's too big
        im_stim_miss_align = AP_align_widefield(animal,day,im_stim_miss);
        im_stim_hit_align = AP_align_widefield(animal,day,im_stim_hit);
        
        batch_vars.im_stim_miss = nansum(cat(5,batch_vars.im_stim_miss,im_stim_miss_align),5);
        batch_vars.im_stim_hit = nansum(cat(5,batch_vars.im_stim_hit,im_stim_hit_align),5);
        
        % Count conditions to divide at end
        batch_vars.n_im_stim_miss = sum(cat(5,batch_vars.n_im_stim_miss,any(any(any(im_stim_miss_align,1),2),3)),5);
        batch_vars.n_im_stim_hit = sum(cat(5,batch_vars.n_im_stim_hit,any(any(any(im_stim_hit_align,1),2),3)),5);
        
        % Prepare for next loop
        AP_print_progress_fraction(curr_day,length(experiments))
        clearvars -except animals protocol curr_animal experiments curr_day animal batch_vars load_parts
        
    end
    
    % Divide sum to get average
    im_stim_miss_avg = bsxfun(@rdivide,batch_vars.im_stim_miss,batch_vars.n_im_stim_miss);
    im_stim_hit_avg = bsxfun(@rdivide,batch_vars.im_stim_hit,batch_vars.n_im_stim_hit);
    
    % Save
    save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
    save([save_path filesep animal '_im_stim_miss_avg'],'im_stim_miss_avg');
    save([save_path filesep animal 'im_stim_hit_avg'],'im_stim_hit_avg');
    
    disp(['Finished ' animal]);
    
end

disp('Finished batch.')


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

animal = 'AP028';
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

animal = 'AP028';
protocol = 'AP_choiceWorldStimPassive';
% protocol = 'stimKalatsky';
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

psth_right_hit_all = cellfun(@(x) nanmean(cat(3,x{:}),3),{batch_vars(:).psth_right_hit},'uni',false);
psth_right_miss_all = cellfun(@(x) nanmean(cat(3,x{:}),3),{batch_vars(:).psth_right_miss},'uni',false);

psth_left_hit_all = cellfun(@(x) nanmean(cat(3,x{:}),3),{batch_vars(:).psth_left_hit},'uni',false);
psth_left_miss_all = cellfun(@(x) nanmean(cat(3,x{:}),3),{batch_vars(:).psth_left_miss},'uni',false);

right_hit = cat(3,psth_right_hit_all{:});
right_hit_norm = nanmean(bsxfun(@rdivide,right_hit,nanmedian(right_hit(:,1:400,:),2))-1,3);

right_miss = cat(3,psth_right_miss_all{:});
right_miss_norm = nanmean(bsxfun(@rdivide,right_miss,nanmedian(right_miss(:,1:400,:),2))-1,3);

left_hit = cat(3,psth_left_hit_all{:});
left_hit_norm = nanmean(bsxfun(@rdivide,left_hit,nanmedian(left_hit(:,1:400,:),2))-1,3);

left_miss = cat(3,psth_left_miss_all{:});
left_miss_norm = nanmean(bsxfun(@rdivide,left_miss,nanmedian(left_miss(:,1:400,:),2))-1,3);

figure; hold on;
raster_window = [-0.5,2.5];
psth_bin_size = 0.001;
t = raster_window(1)+psth_bin_size/2:psth_bin_size:raster_window(2)-psth_bin_size/2;
p_rh = AP_stackplot(right_hit_norm',t,3,false,'k');
p_rm = AP_stackplot(right_miss_norm',t,3,false,'r');
p_lh = AP_stackplot(left_hit_norm',t,3,false,'b');
p_lm = AP_stackplot(left_miss_norm',t,3,false,'m');

line([0,0],ylim,'linestyle','--','color','k');
title('Average across animals')
ylabel('Baseline-normalized response by depth')
xlabel('Time from stim onset')
legend([p_rh(1),p_rm(1),p_lh(1),p_lm(1)],{'Right hit','Right miss','Left hit','Left miss'})




















