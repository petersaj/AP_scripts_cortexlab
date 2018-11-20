%% Batch processes for naive animals
% (these should ultimately be combined with regular animals into a master
% script that prepares all data for figures)


%% !!!!                 DEFINE CORTEX-STRIATUM PARAMETERS               !!!!

%% 1) Get boundaries of striatum across experiments
disp('Getting boundaries of striatum');

animals = {'AP032','AP033','AP034','AP035','AP036'};
protocol = 'AP_choiceWorldStimPassive';

ephys_depth_align = struct;

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
    load_parts.imaging = false;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        AP_load_experiment
       
        % Store MUA corr, template by depth, str depths
        % (this is all calculated during load)
        ephys_depth_align(curr_animal).animal = animal;
        ephys_depth_align(curr_animal).day{curr_day} = day;
        ephys_depth_align(curr_animal).mua_corr{curr_day} = mua_corr;
        ephys_depth_align(curr_animal).templateDepths{curr_day} = templateDepths;
        ephys_depth_align(curr_animal).str_depth(curr_day,:) = str_depth;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal load_parts ephys_depth_align
        
    end
    
    disp(['Finished ' animal])
    
end

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
save([save_path filesep 'ephys_depth_align_naive'],'ephys_depth_align');
disp('Finished batch');


%% 2) Estimate imaging-ephys lambda
disp('Estimating lambda values');

animals = {'AP032','AP033','AP034','AP035','AP036'};

protocol = 'AP_choiceWorldStimPassive';

ctx_str_lambda = struct;

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
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        str_align = 'none';
        AP_load_experiment
                
        % Set upsample factor and sample rate
        upsample_factor = 2;
        sample_rate = (1/median(diff(frame_t)))*upsample_factor;
        
        % Skip the first n seconds to do this
        skip_seconds = 60;
        time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Use all spikes in striatum
        use_spikes = spike_times_timeline(ismember(spike_templates, ...
            find(templateDepths > str_depth(1) & templateDepths <= str_depth(2))));
        binned_spikes = single(histcounts(use_spikes,time_bins));
        
        use_svs = 1:50;
        kernel_t = [-0.3,0.3];
        kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
        zs = [false,true]; % z-score MUA, not V's
        cvfold = 2;
        
        % Resample and get derivative of V
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf(use_svs,:),[],2)',time_bin_centers)';
        
        n_update_lambda = 1;
        lambda_range = [0,1.5]; % ^10 (used to be [3,7] before df/f change
        n_lambdas = 50;
        
        for curr_update_lambda = 1:n_update_lambda
            
            lambdas = logspace(lambda_range(1),lambda_range(2),n_lambdas)';
            explained_var_lambdas = nan(n_lambdas,1);
            
            for curr_lambda_idx = 1:length(lambdas)
                
                curr_lambda = lambdas(curr_lambda_idx);
                
                [~,~,explained_var] = ...
                    AP_regresskernel(dfVdf_resample, ...
                    binned_spikes,kernel_frames,curr_lambda,zs,cvfold);
                
                explained_var_lambdas(curr_lambda_idx) = explained_var.total;
                
            end
            
            lambda_bin_size = diff(lambda_range)/n_lambdas;
            explained_var_lambdas_smoothed = smooth(explained_var_lambdas,10);
            [best_lambda_explained_var,best_lambda_idx] = max(explained_var_lambdas_smoothed);
            best_lambda = lambdas(best_lambda_idx);
            lambda_range = log10([best_lambda,best_lambda]) + ...
                [-lambda_bin_size,lambda_bin_size];
            
        end
                
        ctx_str_lambda(curr_animal).animal = animal;
        ctx_str_lambda(curr_animal).day{curr_day} = day;
        ctx_str_lambda(curr_animal).best_lambda(curr_day) = best_lambda;
        ctx_str_lambda(curr_animal).lambdas{curr_day} = lambdas;        
        ctx_str_lambda(curr_animal).explained_var_lambdas{curr_day} = explained_var_lambdas;        
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except load_parts animals animal curr_animal protocol experiments curr_day animal ctx_str_lambda 
        
    end
end

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
save_fn = 'ctx-str_lambda';
save([save_path filesep save_fn '_naive'],'ctx_str_lambda');
disp('Saved cortex-striatum lambda values');



%% 3) Get kernels at regular depths along striatum
% in the future maybe this should concatenate all imaging for the day to
% have cleaner data (this was done in the code for EJ's paper)
disp('Getting kernels for regular depths in striatum');

n_aligned_depths = 4;

animals = {'AP032','AP033','AP034','AP035','AP036'};

protocol = 'AP_choiceWorldStimPassive';

ephys_kernel_depth = struct;

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
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
                
        % Load data and align striatum by depth
        str_align = 'none';  
        AP_load_experiment         
        
        %%% Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda_naive';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end       
        
        % Try not doing this anymore
%         % Bump it up - no this doesn't make sense, but it makes it better
%         lambda = lambda*10;
        
        % Set upsample value for regression
        upsample_factor = 2;
        sample_rate = (1/median(diff(frame_t)))*upsample_factor;
        
        % Skip the first/last n seconds to do this
        skip_seconds = 60;
        time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get multiunit in ~200 um depths
        n_depths = round(diff(str_depth)/200);
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        [depth_group_n,depth_group] = histc(spikeDepths,depth_group_edges);
        depth_groups_used = unique(depth_group);
        depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);
        
        binned_spikes = zeros(n_depths,length(time_bins)-1);
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        
        % Get rid of NaNs (if no data?)
        binned_spikes(isnan(binned_spikes)) = 0;
        
        use_svs = 1:50;
        kernel_t = [-0.3,0.3];
        kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
        % (lambda taken from estimate above)
        zs = [false,true];
        cvfold = 10;
        
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf(use_svs,:),[],2)',time_bin_centers)';
        
%         % Regress from smoothed diff df/f
%         dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
%             diff(conv2(fVdf(use_svs,:),ones(1,5)/5,'same'),[],2)',time_bin_centers)';
        
        % Regress spikes from cortical fluorescence
        [k,predicted_spikes,explained_var] = ...
            AP_regresskernel(dfVdf_resample, ...
            binned_spikes,kernel_frames,lambda,zs,cvfold);
        
        % Reshape kernel and convert to pixel space
        r = reshape(k,length(use_svs),length(kernel_frames),size(binned_spikes,1));
        
        r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
        for curr_spikes = 1:size(r,3)
            r_px(:,:,:,curr_spikes) = svdFrameReconstruct(Udf(:,:,use_svs),r(:,:,curr_spikes));
        end
        
        % Make single-frames from kernels
        t = kernel_frames/sample_rate;
        use_t = t >= 0 & t <= 0;
        r_px_use = squeeze(nanmean(r_px(:,:,use_t,:),3));
        r_px_use_norm = bsxfun(@rdivide,r_px_use, ...
            permute(max(reshape(r_px_use,[],n_depths),[],1),[1,3,2]));
        r_px_use_norm(isnan(r_px_use_norm)) = 0;
        
        r_px_max_aligned = AP_align_widefield(animal,day,r_px_use_norm);
    
        % Package in structure
        ephys_kernel_depth(curr_animal).animal = animal;
        ephys_kernel_depth(curr_animal).days{curr_day} = day;       
        ephys_kernel_depth(curr_animal).r_px{curr_day} = r_px_max_aligned;
               
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except n_aligned_depths animals animal curr_animal protocol experiments curr_day animal ephys_kernel_depth load_parts
        
    end
    disp(['Finished ' animal]);
end

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
save_fn = ['ephys_kernel_depth'];
save([save_path filesep save_fn '_naive'],'ephys_kernel_depth');
disp('Saved ephys depth kernels');


%% 4) Get template kernels by K-means of depth kernels

% Don't do this for naive mice: use the templates from the trained animals


%% 5) Align striatum recordings with template kernels
disp('Aligning striatum recordings with template kernels');

warning('If done here, needs to be concatenated to original manually')

n_aligned_depths = 4;
animals = {'AP032','AP033','AP034','AP035','AP036'};
protocol = 'AP_choiceWorldStimPassive';

% Load kernels by depths
kernel_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
kernel_fn = ['ephys_kernel_depth_naive'];
load([kernel_path filesep kernel_fn])

ephys_kernel_align_new = struct;

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
    load_parts.imaging = false;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
                
        % Load data
        str_align = 'none';  
        AP_load_experiment 
        
        % Get depth groups (correspond to depth kernels - save in future?)
        n_depths = round(diff(str_depth)/200);
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        [depth_group_n,depth_group] = histc(spikeDepths,depth_group_edges);
        depth_groups_used = unique(depth_group);
        depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);
        
        % Load the template kernels
        kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template'  '_' num2str(n_aligned_depths) '_depths.mat'];
        load(kernel_template_fn);
        
        % Pull out the relevant kernels and normalize
        curr_animal_idx = strcmp(animal,{ephys_kernel_depth.animal});
        curr_day_idx = strcmp(day,ephys_kernel_depth(curr_animal_idx).days);
       
        r_px = ephys_kernel_depth(curr_animal_idx).r_px{curr_day_idx};        
%         kernel_align = mat2gray(bsxfun(@rdivide,r_px,permute(prctile(reshape( ...
%             r_px,[],size(r_px,3)),99,1),[1,3,2])),[0,1]);
        % (don't normalize now)
        kernel_align = r_px;
    
        % Get best match from kernels to kernel templates
        kernel_align_reshape = reshape(kernel_align,[],size(kernel_align,3));
        kernel_align_medfilt = medfilt2(kernel_align_reshape,[1,1]);
        kernel_corr = (zscore(kernel_align_medfilt,[],1)'* ...
            zscore(reshape(kernel_template,[],size(kernel_template,3)),[],1))./ ...
            (size(kernel_align_medfilt,1)-1);
        [kernel_match_corr,kernel_match_raw] = max(kernel_corr,[],2);
        
        % CLEAN UP KERNEL MATCH, NOT SURE HOW TO DO THIS YET
        % Median filter kernel match to get rid of blips
        kernel_match = medfilt1(kernel_match_raw,3);
        
        % If kernel match goes backwards, then set it to nearest neighbor
        replace_kernel_match = kernel_match < cummax(kernel_match);
        kernel_match(replace_kernel_match) = ...
            interp1(find(~replace_kernel_match),kernel_match(~replace_kernel_match), ...
            find(replace_kernel_match),'nearest');
        
        % Assign depth edges and kernel template index numbers
        kernel_match_boundary_idx = unique([1;find(diff(kernel_match) ~= 0)+1;n_depths+1]);
        
        kernel_match_depth_edges = depth_group_edges(kernel_match_boundary_idx);
        kernel_match_idx = kernel_match(kernel_match_boundary_idx(1:end-1));
        
        % Categorize template by kernel match
        aligned_str_depth_group = discretize(spikeDepths,kernel_match_depth_edges,kernel_match_idx);
        
        % Package in structure
        ephys_kernel_align_new(curr_animal).animal = animal;
        ephys_kernel_align_new(curr_animal).days{curr_day} = day;       

        ephys_kernel_align_new(curr_animal).kernel_corr{curr_day} = kernel_corr;
        ephys_kernel_align_new(curr_animal).kernel_match_raw{curr_day} = kernel_match_raw;
        ephys_kernel_align_new(curr_animal).kernel_match{curr_day} = kernel_match;
        
        ephys_kernel_align_new(curr_animal).aligned_str_depth_group{curr_day} = aligned_str_depth_group;
        ephys_kernel_align_new(curr_animal).n_aligned_depths(curr_day) = n_aligned_depths;
                
        AP_print_progress_fraction(curr_day,length(experiments));
        
        clearvars -except n_aligned_depths animals ephys_kernel_depth ...
            animal curr_animal protocol experiments ...
            curr_day animal ephys_kernel_align_new load_parts
        
    end
    disp(['Finished ' animal]);
end

% Overwrite old alignment
ephys_kernel_align = ephys_kernel_align_new;

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
save_fn = ['ephys_kernel_align_' num2str(n_aligned_depths) '_depths'];
save([save_path filesep save_fn '_naive'],'ephys_kernel_align');
disp('Saved ephys kernel alignment');


%% 6) Cortex -> striatum maps with kernel alignment

% ADD THIS IN LATER: important to see if different striatal domains

%% !!!!                                                                                                       !!!!

%% Batch load and save activity from all passive fullscreen (common U)

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

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
save_fn = ['all_trial_activity_Udf_kernel-str_passive_fullscreen_' num2str(n_aligned_depths) '_depths'];
save([save_path filesep save_fn '_naive']);

%% Batch load and save activity from all passive choiceworld (common U)

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

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
save_fn = ['all_trial_activity_Udf_kernel-str_passive_choiceworld_' num2str(n_aligned_depths) '_depths'];
save([save_path filesep save_fn '_naive'],'-v7.3');








