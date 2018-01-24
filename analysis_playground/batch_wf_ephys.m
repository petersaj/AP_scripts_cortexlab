%% Batch widefield area boundaries

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
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
                     
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
        
    % Align images from batch processing
    days = {experiments.day};
    corr_edges_aligned = AP_align_widefield(animal,days,batch_vars.corr_edges);
    
    wf_borders_fig = figure('Name',animal);
    imagesc(nanmean(corr_edges_aligned,3));
    axis image off; colormap(gray); caxis([0,0.05])
    
    fn = ['\\basket.cortexlab.net\data\ajpeters\wf_borders' filesep animal '_wf_borders'];
    saveas(wf_borders_fig,fn);
   
    AP_print_progress_fraction(curr_animal,length(animals))
    
end

disp('Finished batch.')


%% Batch widefield responses to passive stim

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% protocol = 'stimKalatsky';
protocol = 'AP_choiceWorldStimPassive';

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    % Skip if this animal doesn't have this experiment
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

%% Batch widefield choiceworld 

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
    batch_vars.im_stim_earlymove_hit = [];
    batch_vars.im_stim_latemove_hit = [];
    
    batch_vars.im_stim_earlymove_miss = [];
    batch_vars.im_stim_latemove_miss = [];
    
    batch_vars.n_im_stim_earlymove_hit = [];
    batch_vars.n_im_stim_latemove_hit = [];
    
    batch_vars.n_im_stim_earlymove_miss = [];
    batch_vars.n_im_stim_latemove_miss = [];

    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        conditions = unique(block.events.sessionPerformanceValues(1,:));

        % Define trials to use
        n_trials = length(block.paramsValues);
        trial_conditions = signals_events.trialSideValues(1:n_trials).*signals_events.trialContrastValues(1:n_trials);
        trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
        stim_to_move = padarray(wheel_move_time - stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
        stim_to_feedback = padarray(signals_events.responseTimes, ...
            [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
            padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');    
      
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Set options
        surround_window = [-0.5,5];
        framerate = 1./median(diff(frame_t));
        surround_samplerate = 1/(framerate*1);
        t_surround = surround_window(1):surround_samplerate:surround_window(2);
        
        % Average (time course) responses
        im_stim_earlymove_hit = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        im_stim_latemove_hit = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        
        im_stim_earlymove_miss = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        im_stim_latemove_miss = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        
        for curr_condition_idx = 1:length(conditions)
            curr_condition = conditions(curr_condition_idx);
            
            use_stims = use_trials & trial_conditions == curr_condition & trial_outcome == 1 & stim_to_move < 0.5;
            use_stim_onsets = stimOn_times(use_stims);           
            if length(use_stim_onsets) > 5
                stim_surround_times = bsxfun(@plus, use_stim_onsets(:), t_surround);
                peri_stim_v = permute(mean(interp1(frame_t,fVdf',stim_surround_times),1),[3,2,1]);
                im_stim_earlymove_hit(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
            end      
            
            use_stims = use_trials & trial_conditions == curr_condition & trial_outcome == 1 & stim_to_move >= 0.5;
            use_stim_onsets = stimOn_times(use_stims);           
            if length(use_stim_onsets) > 5
                stim_surround_times = bsxfun(@plus, use_stim_onsets(:), t_surround);
                peri_stim_v = permute(mean(interp1(frame_t,fVdf',stim_surround_times),1),[3,2,1]);
                im_stim_latemove_hit(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
            end      
                        
            use_stims = use_trials & trial_conditions == curr_condition & trial_outcome == -1 & stim_to_move < 0.5;
            use_stim_onsets = stimOn_times(use_stims);  
            if length(use_stim_onsets) > 5
                stim_surround_times = bsxfun(@plus, use_stim_onsets(:), t_surround);
                peri_stim_v = permute(mean(interp1(frame_t,fVdf',stim_surround_times),1),[3,2,1]);
                im_stim_earlymove_miss(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
            end
            
            use_stims = use_trials & trial_conditions == curr_condition & trial_outcome == -1 & stim_to_move >= 0.5;
            use_stim_onsets = stimOn_times(use_stims);  
            if length(use_stim_onsets) > 5
                stim_surround_times = bsxfun(@plus, use_stim_onsets(:), t_surround);
                peri_stim_v = permute(mean(interp1(frame_t,fVdf',stim_surround_times),1),[3,2,1]);
                im_stim_latemove_miss(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
            end
            
        end
        
        % Align and average as it goes, otherwise the variable's too big
        im_stim_earlymove_hit_align = AP_align_widefield(animal,day,im_stim_earlymove_hit);
        im_stim_latemove_hit_align = AP_align_widefield(animal,day,im_stim_latemove_hit);
        
        im_stim_earlymove_miss_align = AP_align_widefield(animal,day,im_stim_earlymove_miss);
        im_stim_latemove_miss_align = AP_align_widefield(animal,day,im_stim_latemove_miss);
        
        batch_vars.im_stim_earlymove_hit = nansum(cat(5,batch_vars.im_stim_earlymove_hit,im_stim_earlymove_hit_align),5);
        batch_vars.im_stim_latemove_hit = nansum(cat(5,batch_vars.im_stim_latemove_hit,im_stim_latemove_hit_align),5);
        
        batch_vars.im_stim_earlymove_miss = nansum(cat(5,batch_vars.im_stim_earlymove_miss,im_stim_earlymove_miss_align),5);
        batch_vars.im_stim_latemove_miss = nansum(cat(5,batch_vars.im_stim_latemove_miss,im_stim_latemove_miss_align),5);
        
        % Count conditions to divide at end
        batch_vars.n_im_stim_earlymove_hit = ...
            sum(cat(5,batch_vars.n_im_stim_earlymove_hit,any(any(any(im_stim_earlymove_hit_align,1),2),3)),5);
        batch_vars.n_im_stim_latemove_hit = ...
            sum(cat(5,batch_vars.n_im_stim_latemove_hit,any(any(any(im_stim_latemove_hit_align,1),2),3)),5);
        
        batch_vars.n_im_stim_earlymove_miss = ...
            sum(cat(5,batch_vars.n_im_stim_earlymove_miss,any(any(any(im_stim_earlymove_miss_align,1),2),3)),5);
        batch_vars.n_im_stim_latemove_miss = ...
            sum(cat(5,batch_vars.n_im_stim_latemove_miss,any(any(any(im_stim_latemove_miss_align,1),2),3)),5);
        
        % Prepare for next loop
        AP_print_progress_fraction(curr_day,length(experiments))
        clearvars -except animals protocol curr_animal experiments curr_day animal batch_vars load_parts
        
    end
    
    % Divide sum to get average
    im_stim_earlymove_hit_avg = bsxfun(@rdivide,batch_vars.im_stim_earlymove_hit,batch_vars.n_im_stim_earlymove_hit);
    im_stim_latemove_hit_avg = bsxfun(@rdivide,batch_vars.im_stim_latemove_hit,batch_vars.n_im_stim_latemove_hit);
    
    im_stim_earlymove_miss_avg = bsxfun(@rdivide,batch_vars.im_stim_earlymove_miss,batch_vars.n_im_stim_earlymove_miss);
    im_stim_latemove_miss_avg = bsxfun(@rdivide,batch_vars.im_stim_latemove_miss,batch_vars.n_im_stim_latemove_miss);
    
    % Save
    save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
    save([save_path filesep animal 'im_stim_earlymove_hit_avg'],'im_stim_earlymove_hit_avg','-v7.3');
    save([save_path filesep animal '_im_stim_latemove_hit_avg'],'im_stim_latemove_hit_avg','-v7.3');
    save([save_path filesep animal 'im_stim_earlymove_miss_avg'],'im_stim_earlymove_miss_avg','-v7.3');
    save([save_path filesep animal '_im_stim_latemove_miss_avg'],'im_stim_latemove_miss_avg','-v7.3');
    
    disp(['Finished ' animal]);
    
end

disp('Finished batch.')
warning('This uses -v7.3 and therefore compresses data, switch to dat in the future');


%% Batch widefield > striatum maps

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% protocol = 'vanillaChoiceworld';
% protocol = 'stimSparseNoiseUncorrAsync';
% protocol = 'stimKalatsky';
protocol = 'AP_choiceWorldStimPassive';

batch_vars = struct;
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    % Skip if this animal doesn't have this experiment
    if isempty(experiments)
        continue
    end
    
    disp(animal);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        
        sample_rate = (1/median(diff(frame_t)))*1;
        
        % Skip the first/last n seconds to do this
        skip_seconds = 60;
        
        time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Group multiunit by depth
        n_depth_groups = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
        depth_group_edges_use = depth_group_edges;
        
        [depth_group_n,depth_group] = histcounts(spikeDepths,depth_group_edges_use);
        depth_groups_used = unique(depth_group);
        depth_group_centers = depth_group_edges_use(1:end-1)+(diff(depth_group_edges_use)/2);
        
        binned_spikes = zeros(length(depth_group_edges_use)-1,length(time_bins)-1);
        for curr_depth = 1:length(depth_group_edges_use)-1           
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);           
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
        
        batch_vars(curr_animal).r_px_com{curr_day} = r_px_com;
        batch_vars(curr_animal).r_px_weight{curr_day} = r_px_weight;
        batch_vars(curr_animal).explained_var{curr_day} = explained_var.total;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
      
    disp(['Finished ' animal]);
    
end

% Save
save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_ephys'];
save([save_path filesep 'wf_ephys_maps_' protocol],'batch_vars');

disp('Finished batch');

%% Batch cortex > striatum prediction by condition and depth

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';
batch_vars = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    disp(animal);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        
        sample_rate_factor = 1;
        
        % Group multiunit by depth
        n_depth_groups = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
        
        depth_group = discretize(spikeDepths,depth_group_edges);
        depth_groups_used = unique(depth_group);
        depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);
        
        % Skip the first/last n seconds for prediction
        skip_seconds = 60*1;
        
        sample_rate = (1/median(diff(frame_t)))*sample_rate_factor;
        
        time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        binned_spikes = zeros(length(depth_group_edges)-1,length(time_bins)-1);
        for curr_depth = 1:length(depth_group_edges)-1
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        
        stim_conditions = signals_events.trialContrastValues.*signals_events.trialSideValues;
        conditions = unique(block.events.sessionPerformanceValues(1,:));
        
        % (only use trials that were completed and not repeat trials)
        use_trials = any([signals_events.hitValues;signals_events.missValues],1) & ~signals_events.repeatTrialValues;
        align_times = reshape(stimOn_times(use_trials),[],1);
        
        interval_surround = [-0.5,1.5];
        t_surround = interval_surround(1):1/sample_rate:interval_surround(2);
        t_peri_event = bsxfun(@plus,align_times,t_surround);
        
        % Regress spikes from cortex
        use_svs = 1:50;
        kernel_frames = -35:17;
        lambda = 2e5;
        zs = [false,true];
        cvfold = 5;
        
        fVdf_resample = interp1(frame_t,fVdf(use_svs,:)',time_bin_centers)';
        dfVdf_resample = interp1(conv(frame_t,[1,1]/2,'valid'),diff(fVdf(use_svs,:),[],2)',time_bin_centers)';
        
        % [k,predicted_spikes,explained_var] = ...
        %     AP_regresskernel(fVdf_resample, ...
        %     binned_spikes,kernel_frames,lambda,zs,cvfold);
        [k,predicted_spikes,explained_var] = ...
            AP_regresskernel(dfVdf_resample, ...
            binned_spikes,kernel_frames,lambda,zs,cvfold);
        
        binned_spikes_std = std(binned_spikes,[],2);
        binned_spikes_mean = mean(binned_spikes,2);
        predicted_spikes_reranged = bsxfun(@plus,bsxfun(@times,predicted_spikes, ...
            binned_spikes_std),binned_spikes_mean);
        
        % Find nonlinearity for predicted spikes
        nlin_fit = @(b,x) (x+b(1)).^b(2);
        nlin_params = nan(n_depth_groups,2);
        predicted_spikes_nlin = nan(size(predicted_spikes_reranged));
        for curr_depth = 1:n_depth_groups
            
            % Get the median predicted spikes per binned spikes
            [grp_binned_spikes,grp_predicted_spikes]= ...
                grpstats(predicted_spikes_reranged(curr_depth,:),binned_spikes(curr_depth,:),{'gname','median'});
            grp_binned_spikes = cellfun(@str2num,grp_binned_spikes);
            
            % Pick spikes to fit nonlinearity
            fit_spikes_thresh = prctile(binned_spikes(curr_depth,:),99);
            fit_spikes = grp_binned_spikes > 0 & ...
                grp_binned_spikes < fit_spikes_thresh;
            
            % Less than 5 points to fit, don't bother
            if sum(fit_spikes) < 5
                continue
            end
            
            % Fit the nonlinearity
            fit_options = statset('MaxIter',1000);
            beta0 = [grp_predicted_spikes(1),1];
            beta = nlinfit(grp_predicted_spikes(fit_spikes),grp_binned_spikes(fit_spikes),nlin_fit,beta0,fit_options);
            nlin_params(curr_depth,:) = beta;
            
            % Fit model and rectify
            predicted_spikes_rectif = nlin_fit(beta,predicted_spikes_reranged(curr_depth,:));
            predicted_spikes_rectif(imag(predicted_spikes_rectif) ~= 0) = 0;
            predicted_spikes_nlin(curr_depth,:) = predicted_spikes_rectif;
            
%             % Plot model fit
%             figure; hold on;
%             plot(grp_predicted_spikes,grp_binned_spikes,'k')
%             grp_predicted_spikes_nlin = nlin_fit(beta,grp_predicted_spikes);
%             grp_predicted_spikes_nlin(imag(grp_predicted_spikes_nlin) ~= 0) = 0;
%             plot(grp_predicted_spikes_nlin,grp_binned_spikes,'r')
%             line([0,fit_spikes_thresh],[0,fit_spikes_thresh]);
%             xlabel('Predicted spikes')
%             ylabel('Real spikes')
%             legend({'Raw','Nlin'})
        end
        
        if ~isreal(nlin_params)
            warning('Imaginary parameter fits')
        end
        
        % % Get new explained variance (this can't be right...)
        % sse_real_spikes = sum(bsxfun(@minus,binned_spikes,nanmean(binned_spikes,2)).^2,2);
        % sse_total_residual = sum((predicted_spikes_nlin-binned_spikes).^2,2);
        % explained_var_nlin = (sse_real_spikes - sse_total_residual)./sse_real_spikes;
        
        % Align real and predicted spikes to event
        mua_stim_hit = nan(n_depth_groups,length(t_surround),length(conditions));
        mua_stim_hit_pred = nan(n_depth_groups,length(t_surround),length(conditions));
        
        mua_stim_miss = nan(n_depth_groups,length(t_surround),length(conditions));
        mua_stim_miss_pred = nan(n_depth_groups,length(t_surround),length(conditions));
        
        for curr_depth = 1:n_depth_groups
            for curr_condition_idx = 1:length(conditions)
                mua_stim = interp1(time_bin_centers,binned_spikes(curr_depth,:),t_peri_event);
                mua_stim_pred = interp1(time_bin_centers,predicted_spikes_nlin(curr_depth,:),t_peri_event);
                
                curr_trials = signals_events.hitValues(use_trials) == 1 & stim_conditions(use_trials) == conditions(curr_condition_idx);
                if sum(curr_trials) > 5
                    mua_stim_hit(curr_depth,:,curr_condition_idx) = nanmean(mua_stim(curr_trials,:),1);
                    mua_stim_hit_pred(curr_depth,:,curr_condition_idx) = nanmean(mua_stim_pred(curr_trials,:),1);
                end
                
                curr_trials = signals_events.missValues(use_trials) == 1 & stim_conditions(use_trials) == conditions(curr_condition_idx);
                if sum(curr_trials) > 5
                    mua_stim_miss(curr_depth,:,curr_condition_idx) = nanmean(mua_stim(curr_trials,:),1);
                    mua_stim_miss_pred(curr_depth,:,curr_condition_idx) = nanmean(mua_stim_pred(curr_trials,:),1);
                end
            end
        end
        
        % Store variables
        batch_vars(curr_animal).mua_stim_hit(:,:,:,curr_day) = mua_stim_hit;
        batch_vars(curr_animal).mua_stim_hit_pred(:,:,:,curr_day) = mua_stim_hit_pred;
        
        batch_vars(curr_animal).mua_stim_miss(:,:,:,curr_day) = mua_stim_miss;
        batch_vars(curr_animal).mua_stim_miss_pred(:,:,:,curr_day) = mua_stim_miss_pred;
        
        batch_vars(curr_animal).nlin_params(:,:,curr_day) = nlin_params;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal])
    
end


%% Batch striatum responses to choiceworld

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';
batch_vars = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    disp(animal);
    
    experiments = experiments([experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = false;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment        
        conditions = unique(block.events.sessionPerformanceValues(1,:));
        
        % Define trials to use
        n_trials = length(block.paramsValues);
        trial_conditions = signals_events.trialSideValues(1:n_trials).*signals_events.trialContrastValues(1:n_trials);
        trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
        stim_to_move = padarray(wheel_move_time - stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
        stim_to_feedback = padarray(signals_events.responseTimes, ...
            [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
            padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');    
      
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Group multiunit by depth
        n_depth_groups = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
        depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
        
        depth_group = discretize(spikeDepths,depth_group_edges);
        
        raster_window = [-0.5,3];
        psth_bin_size = 0.001;
        t = raster_window(1):psth_bin_size:raster_window(2);
        t_bins = t(1:end-1) + diff(t);
        
        mua_stim_earlymove_hit = nan(6,length(t_bins),length(conditions));
        mua_stim_latemove_hit = nan(6,length(t_bins),length(conditions));
        
        mua_stim_earlymove_miss = nan(6,length(t_bins),length(conditions));
        mua_stim_latemove_miss = nan(6,length(t_bins),length(conditions));
        
        for curr_depth = 1:n_depth_groups
            
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            
            for curr_condition_idx = 1:length(conditions)
                curr_condition = conditions(curr_condition_idx);                
                
                use_stims = use_trials & trial_conditions == curr_condition & trial_outcome == 1 & stim_to_move < 0.5;
                use_stim_onsets = stimOn_times(use_stims);
                if length(use_stim_onsets) > 5
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_stim_onsets,raster_window,psth_bin_size);
                    mua_stim_earlymove_hit(curr_depth,:,curr_condition_idx) = psth;
                end
                
                use_stims = use_trials & trial_conditions == curr_condition & trial_outcome == 1 & stim_to_move >= 0.5;
                use_stim_onsets = stimOn_times(use_stims);
                if length(use_stim_onsets) > 5
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_stim_onsets,raster_window,psth_bin_size);
                    mua_stim_latemove_hit(curr_depth,:,curr_condition_idx) = psth;
                end
                
                use_stims = use_trials & trial_conditions == curr_condition & trial_outcome == -1 & stim_to_move < 0.5;
                use_stim_onsets = stimOn_times(use_stims);
                if length(use_stim_onsets) > 5
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_stim_onsets,raster_window,psth_bin_size);
                    mua_stim_earlymove_miss(curr_depth,:,curr_condition_idx) = psth;
                end
                
                use_stims = use_trials & trial_conditions == curr_condition & trial_outcome == -1 & stim_to_move >= 0.5;
                use_stim_onsets = stimOn_times(use_stims);
                if length(use_stim_onsets) > 5
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_stim_onsets,raster_window,psth_bin_size);
                    mua_stim_latemove_miss(curr_depth,:,curr_condition_idx) = psth;
                end
                
            end
        end
        
        batch_vars(curr_animal).mua_stim_earlymove_hit(:,:,:,curr_day) = mua_stim_earlymove_hit;
        batch_vars(curr_animal).mua_stim_latemove_hit(:,:,:,curr_day) = mua_stim_latemove_hit;
        
        batch_vars(curr_animal).mua_stim_earlymove_miss(:,:,:,curr_day) = mua_stim_earlymove_miss;
        batch_vars(curr_animal).mua_stim_latemove_miss(:,:,:,curr_day) = mua_stim_latemove_miss;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal])
    
end

% Save
save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
save([save_path filesep 'mua_stim_choiceworld'],'batch_vars');

%% Batch striatum responses to passive

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% protocol = 'stimKalatsky';
protocol = 'AP_choiceWorldStimPassive';

batch_vars = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    disp(animal);
    
    experiments = experiments([experiments.ephys]);
    
    % Skip if this animal doesn't have this experiment
    if isempty(experiments)
        continue
    end
    
    load_parts.cam = false;
    load_parts.imaging = false;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        
        conditions = unique(stimIDs);

        % Group multiunit by depth
        n_depth_groups = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
        depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
        
        depth_group = discretize(spikeDepths,depth_group_edges);
        
        raster_window = [-0.5,5];
        psth_bin_size = 0.001;
        t = raster_window(1):psth_bin_size:raster_window(2);
        t_bins = t(1:end-1) + diff(t);
        
        mua_stim = nan(6,length(t_bins),length(conditions));
        for curr_depth = 1:n_depth_groups
            
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            
            for curr_condition_idx = 1:length(conditions)
                curr_condition = conditions(curr_condition_idx);
                
                use_stims = find(stimIDs == curr_condition);
                use_stim_onsets = stimOn_times(use_stims(2:end));
                use_stim_onsets([1,end]) = [];
                
                if length(use_stim_onsets) > 5
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_stim_onsets,raster_window,psth_bin_size);
                    mua_stim(curr_depth,:,curr_condition_idx) = psth;
                end
                
            end
        end
        
        batch_vars(curr_animal).mua_stim(:,:,:,curr_day) = mua_stim;

        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal])
    
end

% Save
save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\passive'];
save([save_path filesep 'mua_stim_' protocol],'batch_vars');


























