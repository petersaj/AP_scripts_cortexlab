%% post to alyx

onLoad; % initializes missing-http

myAlyx = alyx.getToken([], 'andy', 'xpr1mnt');

use_date = datestr(now,'yyyy-mm-dd');
use_time = datestr(now,'HH:MM');

clear d
d.user = 'andy';

animals = {'AP030','AP031'};

water = [1.2,1.2];
for curr_animal = 1:length(animals)
    d.subject = animals{curr_animal};
    d.water_administered = water(curr_animal);
    d.date_time = sprintf('%sT%s',use_date,use_time);
    newWater = alyx.postData(myAlyx, 'water-administrations', d);
end

weight = [26.7,26.3];
for curr_animal = 1:length(animals)
    d.subject = animals{curr_animal};
    d.weight = weight(curr_animal);
    d.date_time = sprintf('%sT%s',use_date,use_time);
    newWeight = alyx.postData(myAlyx, 'weighings', d);
end


%% Hemidiff in V (V-V_mirror) - keep for future reference
% looks good 

U_r = reshape(U_master(:,:,1:n_vs),[],n_vs);
U_mirror_r = reshape(AP_reflect_widefield(U_master(:,:,1:n_vs)),[],n_vs);
mirror_matrix = U_r'*U_mirror_r;
fluor_allcat_mirror = reshape(transpose( ...
    mirror_matrix*reshape(fluor_allcat,[],n_vs)'),size(fluor_allcat));

% sanity check: plot regular and mirrored
use_trials = trial_contrast_allcat > 0 & trial_side_allcat == 1 & ...
    trial_choice_allcat == -1 & move_t > 0.1 & move_t < 0.3;

V_mean = squeeze(nanmean(fluor_allcat(use_trials,:,:),1))';
V_mirror_mean = squeeze(nanmean(fluor_allcat_mirror(use_trials,:,:),1))';

px_mean = svdFrameReconstruct(U_master(:,:,1:n_vs),V_mean);
px_mirror_mean = svdFrameReconstruct(U_master(:,:,1:n_vs),V_mirror_mean);

AP_image_scroll([px_mean,px_mirror_mean],t);
caxis([-prctile(abs(px_mean(:)),99),prctile(abs(px_mean(:)),99)]);
colormap(brewermap([],'*RdBu'));
axis image;

V_hemidiff = V_mean - V_mirror_mean;
px_hemidiff = svdFrameReconstruct(U_master(:,:,1:n_vs),V_hemidiff);
AP_image_scroll(px_hemidiff,t);
caxis([-prctile(abs(px_hemidiff(:)),99),prctile(abs(px_hemidiff(:)),99)]);
colormap(brewermap([],'*RdBu'));
axis image;

% sanity check
use_trials = trial_contrast_allcat > 0 & trial_side_allcat == 1 & ...
    trial_choice_allcat == -1 & move_t > 0.1 & move_t < 0.3;

V_check = squeeze(nanmean(fluor_allcat_downsamp_filt(use_trials,:,:),1))';
px_check = svdFrameReconstruct(U_master(:,:,1:n_vs),V_check);

AP_image_scroll(px_check,t);
caxis([-prctile(abs(px_check(:)),100),prctile(abs(px_check(:)),100)]);
colormap(brewermap([],'*RdBu'));
axis image;
AP_reference_outline('ccf_aligned','k');


% compatible with code: 

mirror_matrix = reshape(U_master(:,:,1:n_vs),[],n_vs)'* ...
    reshape(AP_reflect_widefield(U_master(:,:,1:n_vs)),[],n_vs);
fluor_allcat_downsamp_mirror = reshape(transpose( ...
    mirror_matrix*reshape(fluor_allcat_downsamp,[],n_vs)'),size(fluor_allcat_downsamp));
fluor_allcat_downsamp = fluor_allcat_downsamp - fluor_allcat_downsamp_mirror;


%% Checking waveform stuff for eventual spike sorting?

% NEED TO ADD IN WINV TO GET THE RIGHT SCALING

% plot templates by channel

figure; hold on;
p = arrayfun(@(x) plot(0,0,'k'),1:size(templates,3));
for curr_template = 1:size(templates,1)
    
    y = permute(templates(curr_template,:,:),[3,2,1]);
    y = y + (max(channel_positions(:,2)) - channel_positions(:,2))/200;   
    x = (1:size(templates,2)) + channel_positions(:,1)*7;
    
    nonzero_channels = squeeze(any(templates(curr_template,:,:),2));
    [~,max_channel] = max(max(abs(templates(curr_template,:,:)),[],2),[],3);
    
    arrayfun(@(ch) set(p(ch),'XData',x(ch,:),'YData',y(ch,:)),1:size(templates,3));
    arrayfun(@(ch) set(p(ch),'Color','r'),find(nonzero_channels));
    arrayfun(@(ch) set(p(ch),'Color','k'),find(~nonzero_channels));
    set(p(max_channel),'Color','b');
    
%     ylim([min(reshape(y(nonzero_channels,:),[],1)), ...
%         max(reshape(y(nonzero_channels,:),[],1))]);
    title(curr_template);
    waitforbuttonpress;
    
end

% Get rid of zero-weight template channels
templates_permute = permute(templates,[3,2,1]);
used_template_channels = any(templates_permute,2);
if length(unique(sum(used_template_channels,1))) ~= 1
    error('Different number of unused channels');
end
templates_used = reshape( ...
    templates_permute(repmat(used_template_channels,1,size(templates_permute,2),1)), ...
    [],size(templates_permute,2),size(templates_permute,3));

template_corr = arrayfun(@(x) nanmean(AP_itril(corrcoef(templates_used(:,:,x)'),-1)),1:size(templates_used,3));
template_channel_skewness = squeeze(skewness(sum(templates_used.^2,2),[],1));
figure;plot(template_channel_skewness,template_corr,'.k');
xlabel('Channel skewness');
ylabel('Channel correlation');

figure;plot3(template_channel_skewness,template_corr,1:size(templates,1),'.k');


% 
% % Plot template skewness by depth
% skewness_time = max(skewness(templates.^2,[],2),[],3);
% skewness_channel = skewness(sum(templates.^2,2),[],3);
% 
% figure;plot3(skewness_time,skewness_channel,templateDepths,'.k')
% axis vis3d
% set(gca,'ZDir','reverse');
% xlabel('Skewness time');
% ylabel('Skewness channel');
% zlabel('Template depth');



%% Trying regression with deconv

clear all
disp('Passive fullscreen trial activity (trained)')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 60;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.2,0.2];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

fluor_all = cell(length(animals),1);
mua_all = cell(length(animals),1);
predicted_mua_std_all = cell(length(animals),1);
wheel_all = cell(length(animals),1);
D_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    % Only use days with choiceworld (sometimes recorded in cortex, no bhv)
    protocol = 'vanillaChoiceworld';
    behavior_experiments = AP_find_experiments(animal,protocol);
    
    protocol = 'stimKalatsky';
    passive_experiments = AP_find_experiments(animal,protocol);
    
    behavior_day = ismember({passive_experiments.day},{behavior_experiments.day});
    
    experiments = passive_experiments([passive_experiments.imaging] & [passive_experiments.ephys] & behavior_day);

    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    fluor_all{curr_animal} = cell(length(experiments),1);
    mua_all{curr_animal} = cell(length(experiments),1);
    predicted_mua_std_all{curr_animal} = cell(length(experiments),1);
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
        raster_sample_rate = 1/(framerate*regression_params.upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);      
        
        use_align = stimOn_times;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% MUA
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
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;       
        
        %%%%%% TESTING %%%%%%%
%         % Get upsampled dVdf's
%         dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
%             diff(fVdf(regression_params.use_svs,:),[],2)',time_bin_centers)';
        
%         % Use smoothed-trace derivative
%         deriv_smooth = 3;       
%         frame_t_smooth_diff = conv(conv(frame_t,ones(1,deriv_smooth)/deriv_smooth,'valid'),[1,1]/2,'valid');
%         dfVdf_resample = interp1(frame_t_smooth_diff,diff( ...
%             convn(fVdf(regression_params.use_svs,:), ...
%             ones(1,deriv_smooth)/deriv_smooth,'valid'),[],2)', ...
%             time_bin_centers,'linear','extrap')';       
        
        % Use deconvolved fluorescence
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\gcamp_kernel\gcamp6s_kernel.mat');        
        gcamp6s_kernel_cat = vertcat(gcamp6s_kernel.regression{:});
        gcamp6s_kernel = nanmean(gcamp6s_kernel_cat./max(gcamp6s_kernel_cat,[],2),1);
        fluor_allcat_deriv = convn(fVdf,gcamp6s_kernel,'same');
        dfVdf_resample = interp1(frame_t, ...
           fluor_allcat_deriv(regression_params.use_svs,:)',time_bin_centers)';
        
        %%%%%% TESTING %%%%%%%
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        
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
        
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;
        
        [~,predicted_spikes_std,explained_var] = ...
            AP_regresskernel(dfVdf_resample, ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant);
        
        event_aligned_predicted_mua_std = ...
            interp1(time_bin_centers,predicted_spikes_std',t_peri_event);
        
        %%% Wheel
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
        predicted_mua_std_all{curr_animal}{curr_day} = event_aligned_predicted_mua_std;
        wheel_all{curr_animal}{curr_day} = event_aligned_wheel;
        D_all{curr_animal}{curr_day} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
    
    clearvars -except n_aligned_depths regression_params animals curr_animal ...
        t fluor_all mua_all predicted_mua_std_all wheel_all D_all
    
end
clearvars -except n_aligned_depths t fluor_all mua_all predicted_mua_std_all wheel_all D_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_passive_fullscreen_DECONVTEST'];
save([save_path filesep save_fn],'-v7.3');

%% DECONV TEST: Choiceworld trial activity
clear all
disp('Choiceworld trial activity')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 60;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.2,0.2];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

fluor_all = cell(length(animals),1);
mua_all = cell(length(animals),1);
predicted_mua_std_all = cell(length(animals),1);
wheel_all = cell(length(animals),1);
reward_all = cell(length(animals),1);
D_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    fluor_all{curr_animal} = cell(length(experiments),1);
    mua_all{curr_animal} = cell(length(experiments),1);
    predicted_mua_std_all{curr_animal} = cell(length(experiments),1);
    wheel_all{curr_animal} = cell(length(experiments),1);
    D_all{curr_animal} = cell(length(experiments),1);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,3];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Align (only to stim)
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf(regression_params.use_svs,:),[],2)',time_bin_centers)';
        
        %%%% TESTING: DECONV
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv(regression_params.use_svs,:)',time_bin_centers)';
        %%%%
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        
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
        
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;
        
%         [~,predicted_spikes_std,explained_var] = ...
%             AP_regresskernel(dfVdf_resample, ...
%             binned_spikes_std,kernel_frames,lambda, ...
%             regression_params.zs,regression_params.cvfold, ...
%             false,regression_params.use_constant);

        %%%% TESTING DECONV
        [~,predicted_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample, ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant);
        %%%%
        
        event_aligned_predicted_mua_std = ...
            interp1(time_bin_centers,predicted_spikes_std',t_peri_event);
        
        %%% Wheel
        event_aligned_wheel_raw = interp1(Timeline.rawDAQTimestamps, ...
            wheel_position,t_peri_event);
        event_aligned_wheel = bsxfun(@minus,event_aligned_wheel_raw, ...
            nanmedian(event_aligned_wheel_raw(:,t < 0),2));
        
        %%% Reward
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];      
        event_aligned_reward = (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_bins(x,:)), ...
            [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate) > 0;
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        % D.response = (trial_choice(use_trials)'+1)/2+1;
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);
        
        % Store activity and stim
        fluor_all{curr_animal}{curr_day} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal}{curr_day} = event_aligned_mua(use_trials,:,:,:);
        predicted_mua_std_all{curr_animal}{curr_day} = event_aligned_predicted_mua_std(use_trials,:,:,:);
        wheel_all{curr_animal}{curr_day} = event_aligned_wheel(use_trials,:,:);
        reward_all{curr_animal}{curr_day} = event_aligned_reward(use_trials,:,:);
        D_all{curr_animal}{curr_day} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except n_aligned_depths regression_params animals curr_animal ...
        t fluor_all mua_all predicted_mua_std_all wheel_all reward_all D_all
    
end
clearvars -except n_aligned_depths t fluor_all mua_all predicted_mua_std_all wheel_all reward_all D_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_DECONVTEST'];
save([save_path filesep save_fn],'-v7.3');

%% DECONVTEST: Cortex -> kernel-aligned striatum map (protocols separately)

% Parameters for regression
n_aligned_depths = 4;
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 60;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.2,0.2];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;
 
%  protocols = {'vanillaChoiceworld', ...
%      'stimSparseNoiseUncorrAsync', ...
%      'stimKalatsky'};

protocols = {'vanillaChoiceworld'};

for protocol = protocols

    n_aligned_depths = 4;
    
    animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
    
    batch_vars = struct;
    for curr_animal = 1:length(animals)
        
        animal = animals{curr_animal};
        
        % Only use days with choiceworld (sometimes recorded in cortex, no bhv)
        behavior_protocol = 'vanillaChoiceworld';
        behavior_experiments = AP_find_experiments(animal,behavior_protocol);
        
        curr_protocol = cell2mat(protocol);
        curr_experiments = AP_find_experiments(animal,curr_protocol);
        
        behavior_day = ismember({curr_experiments.day},{behavior_experiments.day});
        
        experiments = curr_experiments([curr_experiments.imaging] & [curr_experiments.ephys] & behavior_day);
        
        % Skip if this animal doesn't have this experiment
        if isempty(experiments)
            continue
        end
        
        disp(animal);
        
        load_parts.cam = false;
        load_parts.imaging = true;
        load_parts.ephys = true;
        
        for curr_day = 1:length(experiments)
            
            day = experiments(curr_day).day;
            experiment = experiments(curr_day).experiment(end);
            
            % Load data and align striatum by depth
            str_align = 'kernel';
            AP_load_experiment;
            
            %%% Load lambda from previously estimated and saved
            lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
            load(lambda_fn);
            curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
            if any(curr_animal_idx)
                curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
                if any(curr_day_idx)
                    lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
                end
            end
            
            %%% Prepare data for regression
            
            % Get time points to bin
            sample_rate = framerate*regression_params.upsample_factor;
            time_bins = frame_t(find(frame_t > ...
                regression_params.skip_seconds,1)):1/sample_rate: ...
                frame_t(find(frame_t-frame_t(end) < ...
                -regression_params.skip_seconds,1,'last'));
            time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
            
            % Get upsampled dVdf's
            dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
                diff(fVdf(regression_params.use_svs,:),[],2)',time_bin_centers)';
            
            %%%% TESTING: DECONV
            fVdf_deconv = AP_deconv_wf(fVdf);
            fVdf_deconv(isnan(fVdf_deconv)) = 0;
            fVdf_deconv_resample = interp1(frame_t,fVdf_deconv(regression_params.use_svs,:)',time_bin_centers)';
            %%%%
            
            % Get striatum depth group by across-experiment alignment
            n_depths = n_aligned_depths;
            depth_group = aligned_str_depth_group;
            
            binned_spikes = zeros(n_depths,length(time_bin_centers));
            for curr_depth = 1:n_depths
                curr_spike_times = spike_times_timeline(depth_group == curr_depth);
                binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
            end
            
            %%% Regress MUA from cortex
            kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
                round(regression_params.kernel_t(2)*sample_rate);
            
            binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
            binned_spikes_std(isnan(binned_spikes_std)) = 0;
            
%             [k,predicted_spikes,explained_var] = ...
%                 AP_regresskernel(dfVdf_resample, ...
%                 binned_spikes_std,kernel_frames,lambda, ...
%                 regression_params.zs,regression_params.cvfold, ...
%                 false,regression_params.use_constant);
            
            %%%% TESTING DECONV
            [k,predicted_spikes_std,explained_var] = ...
                AP_regresskernel(fVdf_deconv_resample, ...
                binned_spikes_std,kernel_frames,lambda, ...
                regression_params.zs,regression_params.cvfold, ...
                false,regression_params.use_constant);
            %%%%

             
             % Reshape kernel and convert to pixel space
             r = reshape(k,length(regression_params.use_svs),length(kernel_frames),size(binned_spikes,1));
             
             aUdf = single(AP_align_widefield(animal,day,Udf));
             r_px = zeros(size(aUdf,1),size(aUdf,2),size(r,2),size(r,3),'single');
             for curr_spikes = 1:size(r,3)
                 r_px(:,:,:,curr_spikes) = svdFrameReconstruct(aUdf(:,:,regression_params.use_svs),r(:,:,curr_spikes));
             end
             
             % Get center of mass for each pixel
             t = kernel_frames/sample_rate;
             use_t = t >= 0 & t <= 0;
             r_px_max = squeeze(nanmean(r_px(:,:,use_t,:),3));
             r_px_max_norm = bsxfun(@rdivide,r_px_max, ...
                 permute(max(reshape(r_px_max,[],n_depths),[],1),[1,3,2]));
             r_px_max_norm(isnan(r_px_max_norm)) = 0;
             
             r_px_com = sum(bsxfun(@times,r_px_max_norm,permute(1:n_depths,[1,3,2])),3)./sum(r_px_max_norm,3);
             
             r_px_weight = max(abs(r_px_max),[],3);
             
             % Store all variables to save
             batch_vars(curr_animal).r_px{curr_day} = r_px;
             batch_vars(curr_animal).r_px_com{curr_day} = r_px_com;
             batch_vars(curr_animal).r_px_weight{curr_day} = r_px_weight;
             batch_vars(curr_animal).explained_var{curr_day} = explained_var.total;
             
             AP_print_progress_fraction(curr_day,length(experiments));
             clearvars -except regression_params n_aligned_depths ...
                 animals animal curr_animal protocol ...
                 experiments curr_day animal batch_vars load_parts
             
         end
         
         disp(['Finished ' animal]);
         
     end
     
     % Save
     curr_protocol = cell2mat(protocol);
     save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_ephys'];
     save([save_path filesep 'wf_ephys_maps_' curr_protocol '_' num2str(n_aligned_depths) '_depths_kernel_DECONVTEST'],'batch_vars','-v7.3');
     warning('saving -v7.3');
     disp(['Finished batch ' curr_protocol]);
     
 end


%% TESTING: Plot mean waveforms from save_phy

spike_clusters = readNPY(['C:\data_temp\phy\spike_clusters.npy']);
unique_clusters = unique(spike_clusters);
channel_positions = readNPY(['C:\data_temp\phy\channel_positions.npy']);

cluster_mean_waveforms_norm = ...
    (cluster_mean_waveforms-cluster_mean_waveforms(:,1,:))./ ...
    abs(max(max(cluster_mean_waveforms,[],2),[],3));

figure; hold on;
p = arrayfun(@(x) plot(0,0,'k'),1:size(cluster_mean_waveforms_norm,3));
for curr_template = 1:size(cluster_mean_waveforms_norm,1)
    
    y = permute(cluster_mean_waveforms_norm(curr_template,:,:),[3,2,1]);
    y = y + (max(channel_positions(:,2)) - channel_positions(:,2))/100;   
    x = (1:size(cluster_mean_waveforms_norm,2)) + channel_positions(:,1)*8;
    
    nonzero_channels = squeeze(any(cluster_mean_waveforms_norm(curr_template,:,:),2));
    [~,max_channel] = max(max(abs(cluster_mean_waveforms_norm(curr_template,:,:)),[],2),[],3);
    
    arrayfun(@(ch) set(p(ch),'XData',x(ch,:),'YData',y(ch,:)),1:size(cluster_mean_waveforms_norm,3));
    arrayfun(@(ch) set(p(ch),'Color','r'),find(nonzero_channels));
    arrayfun(@(ch) set(p(ch),'Color','k'),find(~nonzero_channels));
    set(p(max_channel),'Color','b');
    
%     ylim([min(reshape(y(nonzero_channels,:),[],1)), ...
%         max(reshape(y(nonzero_channels,:),[],1))]);
    
    curr_cluster = unique_clusters(curr_template);
    title(curr_cluster);
    waitforbuttonpress;
    
end

 

