% Preprocessing to align striatal electrophysiology by relationship to cortex

%% 1) Get boundaries of striatum across experiments
disp('Getting boundaries of striatum across experiments');

animals = {'AP024','AP025','AP026','AP027','AP028','AP029', ...
    'AP032','AP033','AP034','AP035','AP036'};

ephys_depth_align = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    % Find experiments
    % (use only behavior days because cortical recordings afterwards)
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    experiments(~([experiments.imaging] & [experiments.ephys])) = [];
    if isempty(experiments)
        % (if no behavior days then it was a naive mouse - use passive expt)
        protocol = 'AP_choiceWorldStimPassive';
        experiments = AP_find_experiments(animal,protocol);
        experiments(~([experiments.imaging] & [experiments.ephys])) = [];
    end
    
    load_parts.cam = false;
    load_parts.imaging = false;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
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
save([save_path filesep 'ephys_depth_align'],'ephys_depth_align');
disp('Finished batch');

%% 2) Estimate imaging-ephys lambda (concat experiments)
disp('Estimating imaging-ephys lambda');

% Parameters for regression
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 60;
regression_params.upsample_factor = 2;
regression_params.kernel_t = [-0.3,0.3];

animals = {'AP024','AP025','AP026','AP027','AP028','AP029', ...
    'AP032','AP033','AP034','AP035','AP036'};

ctx_str_lambda = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    % Find experiments
    % (use only behavior days because cortical recordings afterwards)
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    days = {experiments([experiments.imaging] & [experiments.ephys]).day};
    if isempty(experiments)
        % (if no behavior days then it was a naive mouse - use passive expt)
        protocol = 'AP_choiceWorldStimPassive';
        experiments = AP_find_experiments(animal,protocol);
        days = {experiments([experiments.imaging] & [experiments.ephys]).day};
    end

    disp(animal); 
    
    for curr_day = 1:length(days)

        day = days{curr_day};
        
        % Find all experiments for that day
        % (get folders with only a number - those're the experiment folders)
        curr_experiments_dir = dir(AP_cortexlab_filename(animal,day,[],'expInfo'));
        curr_experiments_idx = cellfun(@(x) ~isempty(x), regexp({curr_experiments_dir.name},'^\d*$'));
        curr_experiments = cellfun(@str2num,{experiments_dir(curr_experiments_idx).name});
        
        % Loop through experiments, collate data
        time_bin_centers_all = cell(size(curr_experiments));
        dfVdf_all = cell(size(curr_experiments));
        binned_spikes_all = cell(size(curr_experiments));     
        
        disp('Loading and concatenating experiments...');
        for curr_exp = 1:length(curr_experiments)
            experiment = curr_experiments(curr_exp);
            AP_load_experiment;
            
            % Get time points to query
            sample_rate = framerate*regression_params.upsample_factor;
            time_bins = frame_t(find(frame_t > ...
                regression_params.skip_seconds,1)):1/sample_rate: ...
                frame_t(find(frame_t-frame_t(end) < ...
                -regression_params.skip_seconds,1,'last'));
            time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
            time_bin_centers_all{curr_exp} = time_bin_centers;
            
            % Get upsampled dVdf's
            dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
                diff(fVdf(regression_params.use_svs,:),[],2)',time_bin_centers)';
            dfVdf_all{curr_exp} = dfVdf_resample;
            
            % Get striatum depth group by across-experiment alignment
            n_depths = n_aligned_depths;
            depth_group = aligned_str_depth_group;
            
            binned_spikes = zeros(n_depths,length(time_bin_centers));
            for curr_depth = 1:n_depths
                curr_spike_times = spike_times_timeline(depth_group == curr_depth);
                binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
            end
            
            binned_spikes_all{curr_exp} = binned_spikes;
            
            AP_print_progress_fraction(curr_exp,length(curr_experiments));
        end
        
        % Concatenate all data
        time_bin_centers = cat(2,time_bin_centers_all{:});
        dfVdf_resample = cat(2,dfVdf_all{:});
        binned_spikes = cat(2,binned_spikes_all{:});       
        
        % Do regression over range of lambdas, get explained variance
        n_update_lambda = 1;
        lambda_range = [0,1.5]; % ^10 (used to be [3,7] before df/f change
        n_lambdas = 50;
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        zs = [false,true];
        cvfold = 10;
        
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
        
        AP_print_progress_fraction(curr_day,length(days));
        clearvars -except regression_params animals animal curr_animal protocol days curr_day animal ctx_str_lambda 
        
    end
end

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
save_fn = 'ctx-str_lambda';
save([save_path filesep save_fn],'ctx_str_lambda');
disp('Saved cortex-striatum lambda values');

%% 3) Get kernels at regular depths along striatum (concat experiments)
disp('Getting kernels at regular depths along striatum');

% Parameters for regression
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 60;
regression_params.upsample_factor = 2;
regression_params.kernel_t = [-0.3,0.3];

animals = {'AP024','AP025','AP026','AP027','AP028','AP029', ...
    'AP032','AP033','AP034','AP035','AP036'};

ephys_kernel_depth = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    % Find experiments
    % (use only behavior days because cortical recordings afterwards)
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    days = {experiments([experiments.imaging] & [experiments.ephys]).day};
    if isempty(experiments)
        % (if no behavior days then it was a naive mouse - use passive expt)
        protocol = 'AP_choiceWorldStimPassive';
        experiments = AP_find_experiments(animal,protocol);
        days = {experiments([experiments.imaging] & [experiments.ephys]).day};
    end

    disp(animal); 
    
    for curr_day = 1:length(days)
        
        day = days{curr_day};
        
        % Find all experiments for that day
        % (get folders with only a number - those're the experiment folders)
        curr_experiments_dir = dir(AP_cortexlab_filename(animal,day,[],'expInfo'));
        curr_experiments_idx = cellfun(@(x) ~isempty(x), regexp({curr_experiments_dir.name},'^\d*$'));
        curr_experiments = cellfun(@str2num,{experiments_dir(curr_experiments_idx).name});
        
        % Loop through experiments, collate data
        time_bin_centers_all = cell(size(curr_experiments));
        dfVdf_resample_all = cell(size(curr_experiments));
        binned_spikes_all = cell(size(curr_experiments));     
        
        disp('Loading and concatenating experiments...');
        for curr_exp = 1:length(curr_experiments)
            experiment = curr_experiments(curr_exp);
            str_align = 'none'; 
            AP_load_experiment;
            
            % Get time points to query
            sample_rate = framerate*regression_params.upsample_factor;
            time_bins = frame_t(find(frame_t > ...
                regression_params.skip_seconds,1)):1/sample_rate: ...
                frame_t(find(frame_t-frame_t(end) < ...
                -regression_params.skip_seconds,1,'last'));
            time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
            time_bin_centers_all{curr_exp} = time_bin_centers;
            
            % Get upsampled dVdf's
            dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
                diff(fVdf(regression_params.use_svs,:),[],2)',time_bin_centers)';
            dfVdf_resample_all{curr_exp} = dfVdf_resample;
            
            % Get striatum multiunit in ~200 um chunks
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
            
            binned_spikes_all{curr_exp} = binned_spikes;
            
            AP_print_progress_fraction(curr_exp,length(curr_experiments));
        end
        
        % Concatenate all data
        time_bin_centers = cat(2,time_bin_centers_all{:});
        dfVdf_resample = cat(2,dfVdf_resample_all{:});
        binned_spikes = cat(2,binned_spikes_all{:});

        % Regress fluorescence to spikes
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        zs = [false,true];
        cvfold = 10;
        
        [k,predicted_spikes,explained_var] = ...
            AP_regresskernel(dfVdf_resample, ...
            binned_spikes,kernel_frames,lambda,zs,cvfold);
        
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        k_px = zeros(size(Udf_aligned,1),size(Udf_aligned,2),size(k,2),size(k,3),'single');
        for curr_spikes = 1:size(k,3)
            k_px(:,:,:,curr_spikes) = ...
                svdFrameReconstruct(Udf_aligned(:,:,regression_params.use_svs),k(:,:,curr_spikes));
        end
        
        % NaN-out depths with no spikes
        k_px(:,:,:,~any(binned_spikes,2)) = NaN;
        
        % Keep kernel at t = 0
        k_px_frame = k_px(:,:,kernel_frames == 0,:);
    
        % Package in structure
        ephys_kernel_depth(curr_animal).animal = animal;
        ephys_kernel_depth(curr_animal).days{curr_day} = day;       
        ephys_kernel_depth(curr_animal).k_px{curr_day} = k_px_frame;
               
        AP_print_progress_fraction(curr_day,length(days));
        clearvars -except regression_params animals animal curr_animal days curr_day ephys_kernel_depth
        
    end
    disp(['Finished ' animal]);
end

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
save_fn = ['ephys_kernel_depth'];
save([save_path filesep save_fn],'ephys_kernel_depth');
disp('Saved ephys depth kernels');

%% 4) ** Get template kernels by K-means of depth kernels **
disp('Getting template kernels');

warning('Change this to deterministic, not K-means');

% Load kernels by depths
kernel_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
kernel_fn = ['ephys_kernel_depth'];
load([kernel_path filesep kernel_fn])

% Concatenate all kernels, do K-means
k_px_cat = [ephys_kernel_depth(:).k_px];
k_px_cat = cat(3,k_px_cat{:});

% Not normalizing for now
% k_px_cat_norm = mat2gray(bsxfun(@rdivide,k_px_cat,permute(prctile(reshape( ...
%     k_px_cat,[],size(k_px_cat,3)),99,1),[1,3,2])),[0,1]);
k_px_cat_norm = k_px_cat;

k_px_cat_norm_reshape = reshape(k_px_cat_norm,[],size(k_px_cat_norm,3));

use_k_px = find(std(k_px_cat_norm_reshape,[],1) ~= 0);

n_aligned_depths = 4;
kidx = kmeans(k_px_cat_norm_reshape(:,use_k_px)',n_aligned_depths,'Distance','correlation');

% Get average depth for each group
total_depths = 1:max(cellfun(@(x) size(x,3),[ephys_kernel_depth.k_px]));
k_px_depths = cellfun(@(x) total_depths(end-size(x,3)+1:end),[ephys_kernel_depth.k_px],'uni',false);
k_px_depth_cat = horzcat(k_px_depths{:});
k_px_depth_grp = grpstats(k_px_depth_cat(use_k_px),kidx);

% Sort by k-means and group depth
[~,all_depth_sort_idx] = sort(k_px_depth_cat(use_k_px));
[~,k_sort_idx] = sort(kidx);
[~,depth_sort_idx] = sort(k_px_depth_grp);

% Plot k-means groups by depth
k_grp = reshape(grpstats(k_px_cat_norm_reshape(:,use_k_px)',kidx)', ...
    size(k_px_cat,1),size(k_px_cat,2),[]);
k_grp_ordered = k_grp(:,:,depth_sort_idx);

figure;
for i = 1:n_aligned_depths
    subplot(1,n_aligned_depths,i);
    imagesc(reshape(k_grp_ordered(:,:,i),size(k_px_cat,1),[]))
    caxis([-max(abs(k_grp(:))),max(abs(k_grp(:)))]);
    axis image off;
    colormap(brewermap([],'*RdBu'));
    AP_reference_outline('ccf_aligned','k');
end

% Plot correlations of kernels by kidx (renumber by depth)
kidx_depth = depth_sort_idx(kidx);
[~,k_depth_sort_idx] = sort(kidx_depth);
figure;
imagesc(corrcoef(k_px_cat_norm_reshape(:,use_k_px(k_depth_sort_idx))))
title('Sorted kernel correlations');
caxis([-0.5,0.5]); colormap(brewermap([],'*RdBu'));
axis square;
for i = 2:n_aligned_depths
    line(xlim,repmat(sum(kidx_depth < i),2,1),'color','k','linewidth',2);
    line(repmat(sum(kidx_depth < i),2,1),ylim,'color','k','linewidth',2);
end

% % Save template kernels
% kernel_template = k_grp_ordered;
% kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template'  '_' num2str(n_aligned_depths) '_depths'];
% save(kernel_template_fn,'kernel_template');
% disp('Saved kernel template');


%% 5) Align striatum recordings from template kernels
disp('Aligning striatum recordings from template kernels');

n_aligned_depths = 4;
animals = {'AP024','AP025','AP026','AP027','AP028','AP029', ...
    'AP032','AP033','AP034','AP035','AP036'};

% Load kernels by depths
kernel_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
kernel_fn = ['ephys_kernel_depth'];
load([kernel_path filesep kernel_fn])

ephys_kernel_align_new = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    % Find experiments
    % (use only behavior days because cortical recordings afterwards)
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    days = {experiments([experiments.imaging] & [experiments.ephys]).day};
    if isempty(experiments)
        % (if no behavior days then it was a naive mouse - use passive expt)
        protocol = 'AP_choiceWorldStimPassive';
        experiments = AP_find_experiments(animal,protocol);
        days = {experiments([experiments.imaging] & [experiments.ephys]).day};
    end

    disp(animal);  
    
    load_parts.cam = false;
    load_parts.imaging = false;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
                
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
       
        k_px = ephys_kernel_depth(curr_animal_idx).k_px{curr_day_idx};        
%         kernel_align = mat2gray(bsxfun(@rdivide,k_px,permute(prctile(reshape( ...
%             k_px,[],size(k_px,3)),99,1),[1,3,2])),[0,1]);
        % (don't normalize now)
        kernel_align = k_px;
    
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
                
        AP_print_progress_fraction(curr_day,length(days));
        
        clearvars -except n_aligned_depths ...
            animals animal curr_animal days curr_day...
            ephys_kernel_depth ephys_kernel_align_new load_parts
        
    end
    disp(['Finished ' animal]);
end

% Overwrite old alignment
ephys_kernel_align = ephys_kernel_align_new;

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
save_fn = ['ephys_kernel_align_' num2str(n_aligned_depths) '_depths'];
save([save_path filesep save_fn],'ephys_kernel_align');
disp('Saved ephys kernel alignment');

%% 6) Cortex -> kernel-aligned striatum regression (concat experiments)
disp('Cortex -> kernel-aligned striatum regression');

% Parameters for regression
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 60;
regression_params.upsample_factor = 2;
regression_params.kernel_t = [-0.3,0.3];

animals = {'AP024','AP025','AP026','AP027','AP028','AP029', ...
    'AP032','AP033','AP034','AP035','AP036'};

batch_vars = struct;
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    % Find experiments
    % (use only behavior days because cortical recordings afterwards)
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    days = {experiments([experiments.imaging] & [experiments.ephys]).day};
    if isempty(experiments)
        % (if no behavior days then it was a naive mouse - use passive expt)
        protocol = 'AP_choiceWorldStimPassive';
        experiments = AP_find_experiments(animal,protocol);
        days = {experiments([experiments.imaging] & [experiments.ephys]).day};
    end

    disp(animal);    
    
    for curr_day = 1:length(days)
        
        day = days{curr_day};
        
        % Find all experiments for that day
        % (get folders with only a number - those're the experiment folders)
        curr_experiments_dir = dir(AP_cortexlab_filename(animal,day,[],'expInfo'));
        curr_experiments_idx = cellfun(@(x) ~isempty(x), regexp({curr_experiments_dir.name},'^\d*$'));
        curr_experiments = cellfun(@str2num,{experiments_dir(curr_experiments_idx).name});
        
        % Loop through experiments, collate data
        time_bin_centers_all = cell(size(curr_experiments));
        dfVdf_resample_all = cell(size(curr_experiments));
        binned_spikes_all = cell(size(curr_experiments));     
        
        disp('Loading and concatenating experiments...');
        for curr_exp = 1:length(curr_experiments)
            experiment = curr_experiments(curr_exp);
            AP_load_experiment;
            
            % Get time points to query
            sample_rate = framerate*regression_params.upsample_factor;
            time_bins = frame_t(find(frame_t > ...
                regression_params.skip_seconds,1)):1/sample_rate: ...
                frame_t(find(frame_t-frame_t(end) < ...
                -regression_params.skip_seconds,1,'last'));
            time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
            time_bin_centers_all{curr_exp} = time_bin_centers;
            
            % Get upsampled dVdf's
            dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
                diff(fVdf(regression_params.use_svs,:),[],2)',time_bin_centers)';
            dfVdf_resample_all{curr_exp} = dfVdf_resample;
            
            % Get striatum depth group by across-experiment alignment
            n_depths = n_aligned_depths;
            depth_group = aligned_str_depth_group;
            
            binned_spikes = zeros(n_depths,length(time_bin_centers));
            for curr_depth = 1:n_depths
                curr_spike_times = spike_times_timeline(depth_group == curr_depth);
                binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
            end
            
            binned_spikes_all{curr_exp} = binned_spikes;
            
            AP_print_progress_fraction(curr_exp,length(curr_experiments));
        end
        
        % Concatenate all data
        time_bin_centers = cat(2,time_bin_centers_all{:});
        dfVdf_resample = cat(2,dfVdf_resample_all{:});
        binned_spikes = cat(2,binned_spikes_all{:});
        
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
        
        % Regress fluorescence to spikes
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        zs = [false,true];
        cvfold = 10;
        
        [k,predicted_spikes,explained_var] = ...
            AP_regresskernel(dfVdf_resample, ...
            binned_spikes,kernel_frames,lambda,zs,cvfold);
        
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        k_px = zeros(size(Udf_aligned,1),size(Udf_aligned,2),size(k,2),size(k,3),'single');
        for curr_spikes = 1:size(k,3)
            k_px(:,:,:,curr_spikes) = ...
                svdFrameReconstruct(Udf_aligned(:,:,regression_params.use_svs),k(:,:,curr_spikes));
        end
        
        % NaN-out depths with no spikes
        k_px(:,:,:,~any(binned_spikes,2)) = NaN;
        
        % Store kernel pixels
        batch_vars(curr_animal).animal = animal;
        batch_vars(curr_animal).regression_params = regression_params;
        batch_vars(curr_animal).day{curr_day} = day;
        batch_vars(curr_animal).t{curr_day} = kernel_frames/sample_rate;
        batch_vars(curr_animal).k_px{curr_day} = k_px;
        batch_vars(curr_animal).explained_var{curr_day} = explained_var.total;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except regression_params animals animal curr_animal days curr_day batch_vars        
    end
    disp(['Finished ' animal]);
end

% Save
save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_ephys'];
save([save_path filesep 'wf_ephys_maps_concat_' num2str(n_aligned_depths) '_depths_kernel'],'batch_vars','-v7.3');
warning('saving -v7.3');
disp('Finished batch');



