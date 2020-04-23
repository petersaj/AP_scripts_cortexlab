%% Playground for assorted paper revision code tests

%% Testing: "baseline" cortical prediction: remove cortical regions

% Set example experiment to use
animal = 'AP025';
day = '2017-10-04';
experiment = 1;
AP_load_experiment;


% Align U, recast V's excluding high kernel weight regions
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
Udf_aligned = AP_align_widefield(Udf,animal,day);

kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
load(kernel_roi_fn);

zero_region_idx = 1;
wf_region_zero = kernel_roi.bw(:,:,zero_region_idx) + ...
    AP_reflect_widefield(kernel_roi.bw(:,:,zero_region_idx));
Udf_aligned_regionzero = Udf_aligned.*~wf_region_zero;
fVdf_regionzero_altU = ChangeU(Udf_aligned,fVdf,Udf_aligned_regionzero);

% (those U's aren't orthonormal, recast back to original Udf_aligned)
fVdf_regionzero = ChangeU(Udf_aligned_regionzero,fVdf_regionzero_altU,Udf_aligned);



fVdf_deconv = AP_deconv_wf(fVdf_regionzero);






% Parameters for regression
regression_params.use_svs = 1:100;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

% Get time points to bin
sample_rate = framerate*regression_params.upsample_factor;
time_bins = frame_t(find(frame_t > ...
    regression_params.skip_seconds,1)):1/sample_rate: ...
    frame_t(find(frame_t-frame_t(end) < ...
    -regression_params.skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Resample deconvolved fluorescence
fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';

% Bin spikes to match widefield frames
binned_spikes = nan(n_depths,length(time_bin_centers));
for curr_depth = 1:n_depths
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    
    % (skip if no spikes at this depth)
    if isempty(curr_spike_times)
        continue
    end
    
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
end
binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);

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

kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
    round(regression_params.kernel_t(2)*sample_rate);

% Regress cortex to striatum (region excluded)
[ctx_str_k,ctxpred_spikes_std,explained_var] = ...
    AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
    binned_spikes_std,kernel_frames,lambda, ...
    regression_params.zs,regression_params.cvfold, ...
    true,regression_params.use_constant);


% Reshape kernel and convert to pixel space
k_px = zeros(size(Udf_aligned_regionzero,1),size(Udf_aligned_regionzero,2), ...
    size(ctx_str_k{1},2),size(ctx_str_k{1},3));
for curr_spikes = 1:size(ctx_str_k{1},3)
    k_px(:,:,:,curr_spikes) = ...
        svdFrameReconstruct(Udf_aligned_regionzero(:,:,regression_params.use_svs), ...
        ctx_str_k{1}(:,:,curr_spikes));
end

          



%% QUICK DIRTY VERSION: regress on trial data
% I doubt this is even usable because it's not comparable to doing it on
% each experiment basis so this whole thing is just a goddamn proof of
% concept



% (load in a dataset first)



% Regress kernel ROI activity to striatum domain activity (per recording)
regression_params.use_svs = 1:100;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;
lambda = 20;
kernel_frames = floor(regression_params.kernel_t(1)*sample_rate): ...
    ceil(regression_params.kernel_t(2)*sample_rate);

mua_allcat_exp = mat2cell(mua_allcat,trials_recording,length(t),n_depths);
fluor_allcat_exp = mat2cell(fluor_allcat,use_split,length(t),n_vs);

mua_ctxtrialpred_exp = cellfun(@(x) nan(size(x)),mua_allcat_exp,'uni',false);
mua_ctxtrialpred_k = nan(length(regression_params.use_svs),length(kernel_frames),n_depths,length(use_split));
mua_ctxtrialpred_regionzero_exp = cellfun(@(x) nan(size(x)),mua_allcat_exp,'uni',false);
mua_ctxtrialpred_regionzero_k = nan(n_vs,37,n_depths,length(use_split));

for curr_exp = 1:length(trials_recording)
    for curr_depth = 1:n_depths
        
        curr_mua = reshape(mua_allcat_exp{curr_exp}(:,:,curr_depth)',1,[]);
        curr_fluor = reshape(permute( ...
            fluor_allcat_exp{curr_exp},[2,1,3]),[],n_vs)';
        
        %         % Zero out region of cortex related to each depth
        %         curr_fluor_regionzero = nan(n_vs,size(curr_fluor,2),n_depths);
        %         for curr_regionzero = 1:n_depths
        %             wf_region_zero = kernel_roi.bw(:,:,curr_regionzero) + ...
        %                 AP_reflect_widefield(kernel_roi.bw(:,:,curr_regionzero));
        %             U_master_regionzero = U_master.*~wf_region_zero;
        %             fVdf_regionzero_altU = ChangeU(U_master(:,:,1:n_vs),curr_fluor,U_master_regionzero(:,:,1:n_vs));
        %             % (those U's aren't orthonormal, recast back to original Udf_aligned)
        %             fVdf_regionzero = ChangeU(U_master_regionzero(:,:,1:n_vs),fVdf_regionzero_altU,U_master(:,:,1:n_vs));
        %             curr_fluor_regionzero(:,:,curr_regionzero) = fVdf_regionzero;
        %         end
        
        % Load kernel templates
        n_aligned_depths = 4;
        kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template_' num2str(n_aligned_depths) '_depths.mat'];
        load(kernel_template_fn);
        
        % Get regions to zero based on kernel template weight
        frac_max_weight = 0.1; % zero pixels > max weight * this
        wf_region_zero_depths = kernel_template > ...
            max(reshape(kernel_template,[],1,n_depths),[],1)*frac_max_weight;
        wf_region_zero = wf_region_zero_depths(:,:,curr_depth);
        
        % Zero out only current depth region      
        U_master_regionzero = U_master.*~wf_region_zero;
        fVdf_regionzero_altU = ChangeU(U_master(:,:,1:n_vs),curr_fluor,U_master_regionzero(:,:,1:n_vs));
        % (those U's aren't orthonormal, recast back to original Udf_aligned)
        fVdf_regionzero = ChangeU(U_master_regionzero(:,:,1:n_vs),fVdf_regionzero_altU,U_master(:,:,1:n_vs));
        curr_fluor_regionzero = fVdf_regionzero;
        
        % Skip if no data
        if all(isnan(curr_mua))
            continue
        end
        
        % Set discontinuities in trial data
        trial_discontinuities = false(size(mua_allcat_exp{curr_exp}(:,:,curr_depth)));
        trial_discontinuities(:,1) = true;
        trial_discontinuities = reshape(trial_discontinuities',[],1)';
        
        % Full cortex regression
        [k_fluor,curr_mua_fluorpred,explained_var] = ...
            AP_regresskernel(curr_fluor(regression_params.use_svs,:), ...
            curr_mua,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant,trial_discontinuities);
        
        mua_ctxtrialpred_k(:,:,curr_depth,curr_exp) = k_fluor;
        mua_ctxtrialpred_exp{curr_exp}(:,:,curr_depth) = ...
            reshape(curr_mua_fluorpred,length(t),[])';
        
        % Region-zeroed cortex regression
        [k_fluorregionzero,curr_mua_fluorregionzeropred,explained_var] = ...
            AP_regresskernel(curr_fluor_regionzero, ...
            curr_mua,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant,trial_discontinuities);
        
        mua_ctxtrialpred_regionzero_k(:,:,curr_depth,curr_exp) = k_fluorregionzero;
        mua_ctxtrialpred_regionzero_exp{curr_exp}(:,:,curr_depth) = ...
            reshape(curr_mua_fluorregionzeropred,length(t),[])';
        
    end   
    AP_print_progress_fraction(curr_exp,length(trials_recording));   
end

% after this: make fig 3 and S8 with the ctxpred, the ctxtrialpred, and the
% ctxtrialpred_regionzero (ctxpred should ~ ctxtrialpred, regionzero should
% hopefully be much worse)








%% CTX EPHYS: passive trial activity WITH NEW MUA FILTER
clear all

animals = {'AP043','AP060','AP061'};
protocol = 'AP_lcrGratingPassive';
recording_site = {'ctxstrephys_str_muafilt','ctxstrephys_ctx_muafilt'};

% Load ephys alignment
vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
load(vis_ctx_ephys_fn);

for curr_site = 1:2
    
    % Initialize save variable
    trial_data_all = struct;
    
    for curr_animal = 1:length(animals)
        
        animal = animals{curr_animal};
        
        experiments = AP_find_experiments(animal,protocol);
        experiments = experiments([experiments.imaging] & [experiments.ephys]);
        
        disp(['Loading ' animal]);
        
        for curr_day = 1:length(experiments)
            
            day = experiments(curr_day).day;
            experiment = experiments(curr_day).experiment(1);
            
            % Load experiment
            site = curr_site;
            switch site
                case 1
                    % Site 1 = striautm
                    str_align = 'kernel';
                case 2
                    % Site 2 = cortex
                    str_align = 'none';
            end
            AP_load_experiment;
            
            % If cortical experiment - depth-align cortical MUA
            if site == 2
                % Get depth group boundaries
                % (these are just manual based on the CSD)
                n_aligned_depths = 2;
                ctx_depth_edges = linspace(0,1200,n_aligned_depths+1);
                
                % Depth-align templates
                curr_animal_idx = strcmp(animal,{vis_ctx_ephys.animal});
                curr_day_idx = strcmp(day,vis_ctx_ephys(curr_animal_idx).day);
                curr_csd_depth = vis_ctx_ephys(curr_animal_idx).stim_csd_depth{curr_day_idx};
                curr_csd_depth_aligned = vis_ctx_ephys(curr_animal_idx).stim_csd_depth_aligned{curr_day_idx};
                
                template_depths_aligned = interp1(curr_csd_depth,curr_csd_depth_aligned,template_depths);
                spike_depths_aligned = template_depths_aligned(spike_templates);
                
                % Find cortex end by largest gap between templates
                sorted_template_depths = sort([template_depths_aligned]);
                [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
                ctx_end = sorted_template_depths(max_gap_idx)+1;               
                ctx_depth = [sorted_template_depths(1),ctx_end];
                ctx_units = template_depths_aligned <= ctx_depth(2);
                
                % Assign templates to depth groups
                % (NOTE: variable still called 'str' to make it easier to
                % fit into other code)
                ctx_spike_depths = spike_depths_aligned;
                ctx_spike_depths(spike_depths_aligned < ctx_depth(1) | spike_depths_aligned > ctx_depth(2)) = NaN;                        
                aligned_str_depth_group = discretize(ctx_spike_depths,ctx_depth_edges);               
            end
            
            % Pull out trial data
            filter_mua = true; % Filter MUA to match widefield
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
                        
            % Clear for next loop
            clearvars -except ...
                vis_ctx_ephys ...
                recording_site curr_site ...
                animals curr_animal animal protocol experiments curr_day ...
                trial_data_all
            
            AP_print_progress_fraction(curr_day,length(experiments));
        end
    end
    
    clearvars -except ...
        vis_ctx_ephys ...
        recording_site curr_site animals protocol ...
        trial_data_all
    
    disp('Finished loading all')
    
    save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
    save_fn = ['trial_activity_' protocol '_' recording_site{curr_site}];
    save([save_path filesep save_fn],'-v7.3');
    disp(['Saved ' save_fn]);
    
end







