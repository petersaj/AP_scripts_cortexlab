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



%% DECONV TEST: trained passive fullscreen

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
        
        %%%% TESTING: DECONV
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv(regression_params.use_svs,:)',time_bin_centers)';
        %%%%
        
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

%% DECONV TEST: Passive fullscreen trial activity (naive)

clear all
disp('Passive fullscreen trial activity (naive)')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 60;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.2,0.2];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP032','AP033','AP034','AP035','AP036'};

fluor_all = cell(length(animals),1);
mua_all = cell(length(animals),1);
predicted_mua_std_all = cell(length(animals),1);
wheel_all = cell(length(animals),1);
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
save_fn = ['trial_activity_passive_fullscreen_naive_DECONVTEST'];
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
mua_ctxpred_all = cell(length(animals),1);
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
        
        % Filter and std-normalize spikes
        lowpassCutoff = 6; % Hz
        [b100s, a100s] = butter(2, lowpassCutoff/((sample_rate)/2), 'low');
        binned_spikes_filt = filter(b100s,a100s,binned_spikes,[],2);
        binned_spikes_filt_std = binned_spikes_filt./nanstd(binned_spikes_filt,[],2);
        binned_spikes_filt_std(isnan(binned_spikes_filt_std)) = 0;
        
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
        
%         [~,predicted_spikes_std,explained_var] = ...
%             AP_regresskernel(dfVdf_resample, ...
%             binned_spikes_std,kernel_frames,lambda, ...
%             regression_params.zs,regression_params.cvfold, ...
%             false,regression_params.use_constant);

        %%%% TESTING DECONV
        [k,predicted_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample, ...
            binned_spikes_filt_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        %%%%
        
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        predicted_spikes = (predicted_spikes_std - squeeze(k{end})).* ...
            nanstd(binned_spikes_filt,[],2) + ...
            nanstd(binned_spikes_filt,[],2).*squeeze(k{end});
        
        event_aligned_predicted_mua = ...
            interp1(time_bin_centers,predicted_spikes',t_peri_event)./raster_sample_rate;
        
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
        mua_ctxpred_all{curr_animal}{curr_day} = event_aligned_predicted_mua(use_trials,:,:,:);
        wheel_all{curr_animal}{curr_day} = event_aligned_wheel(use_trials,:,:);
        reward_all{curr_animal}{curr_day} = event_aligned_reward(use_trials,:,:);
        D_all{curr_animal}{curr_day} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except n_aligned_depths regression_params animals curr_animal ...
        t fluor_all mua_all mua_ctxpred_all wheel_all reward_all D_all
    
end
clearvars -except n_aligned_depths t fluor_all mua_all mua_ctxpred_all wheel_all reward_all D_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_DECONVTEST'];
save([save_path filesep save_fn],'-v7.3');




%% DECONV TEST: Choiceworld trial activity (split pre/post spike)
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

        %%%% TESTING DECONV SPLIT FORWARD/BACK SPIKE
        [~,predicted_spikes_std_backward,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample, ...
            binned_spikes_std,kernel_frames(kernel_frames < 0),lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant);
        
        [~,predicted_spikes_std_forward,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample, ...
            binned_spikes_std,kernel_frames(kernel_frames > 0),lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant);
        %%%%
        
        event_aligned_predicted_mua_std_backward = ...
            interp1(time_bin_centers,predicted_spikes_std_backward',t_peri_event);
        event_aligned_predicted_mua_std_forward = ...
            interp1(time_bin_centers,predicted_spikes_std_forward',t_peri_event);
        
        event_aligned_predicted_mua_std = cat(4, ...
            event_aligned_predicted_mua_std_forward, ...
            event_aligned_predicted_mua_std_backward);
        
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
save_fn = ['trial_activity_choiceworld_FORWARDBACKWARDTEST'];
save([save_path filesep save_fn],'-v7.3');


%% DECONVTEST: Cortex -> kernel-aligned striatum map (protocols separately)

n_aligned_depths = 4;

% Parameters for regression
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;
 
%  protocols = {'vanillaChoiceworld', ...
%      'stimSparseNoiseUncorrAsync', ...
%      'stimKalatsky'};

protocols = {'vanillaChoiceworld'};

for protocol = protocols
    
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
            warning('Overriding lambda');
            lambda = 5;
            
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
            
            binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
            binned_spikes_std(isnan(binned_spikes_std)) = 0;         
                                 
            %%% Regress MUA from cortex
            kernel_frames = floor(regression_params.kernel_t(1)*sample_rate): ...
                ceil(regression_params.kernel_t(2)*sample_rate);
            
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

 
%% (to plot fluorescence by trial group)

% Get conditions for each trial, plot selected
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];

contrast_side_col = colormap_BlueWhiteRed(5);
contrast_side_col(6,:) = 0;
contrast_side_val = unique(sort([-contrasts,contrasts]))';

% % contrast, side, choice
plot_conditions = ...
    [contrasts(2:end),contrasts(2:end); ...
    -ones(1,5),ones(1,5); ...
    ones(1,5),-ones(1,5)]';

% plot_conditions = ...
%     [1,0,1,1,0,1; ...
%     1,-1,-1,-1,-1,1; ...
%     -1,-1,-1,1,1,1]';

% plot_conditions = ...
%     [1,0,1,1; ...
%     -1,-1,1,-1; ...
%     -1,-1,-1,1]';

% plot_conditions = ...
%     [1; ...
%     1,; ...
%     -1]';

% plot_conditions = ...
%     [0.06,0.06,0.06; ...
%     1,1,-1; ...
%     -1,1,-1]';

use_rxn = move_t > 0 & move_t < 0.5;
    
[~,plot_id] = ismember( ...
    [trial_contrast_allcat,trial_side_allcat,trial_choice_allcat], ...
    plot_conditions,'rows');

move_align = true;

spacing = nanstd(fluor_roi_deconv(:))*5;

fluor_roi_deconv_hemidiff = fluor_roi_deconv(:,:,1:size(wf_roi,1)) - fluor_roi_deconv(:,:,size(wf_roi,1)+1:end);

figure; hold on;
for curr_plot_condition = 1:size(plot_conditions,1)
    
    curr_trials = plot_id == curr_plot_condition & use_rxn;    
    curr_data = fluor_roi_deconv_hemidiff(curr_trials,:,:);
    
    if move_align
        % Re-align to movement onset
        t_leeway = 0.5;
        leeway_samples = round(t_leeway*sample_rate);
        curr_move_idx = move_idx(curr_trials);
        for i = 1:size(curr_data,1)
            curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
        end
    end
    
    curr_data_mean = squeeze(nanmean(curr_data,1));

    curr_contrast_side = max(trial_contrast_allcat(curr_trials))*sign(max(trial_side_allcat(curr_trials)));
    contrast_side_idx = find(curr_contrast_side == contrast_side_val);
    curr_col = contrast_side_col(contrast_side_idx,:);   
    
    curr_linewidth = mean(trial_choice_allcat(curr_trials))+2;
    
    if curr_contrast_side == 0
        switch max(trial_choice_allcat(curr_trials))
            case -1
                curr_col = 'm';
            case 1
                curr_col = 'c';
        end
    end
    
    AP_stackplot(curr_data_mean,t,spacing,false,curr_col,{wf_roi(:,1).area});
    
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');
title('Widefield ROI');



% (rank difference)
trial_groups = {'Stim','Move onset','Reward'};
figure;
for curr_group = 1:length(trial_groups)
    
    curr_regressor_idx = strcmp(trial_groups{curr_group},regressor_labels);
    
    fluor_predicted_roi_reduced = permute(reshape( ...
        AP_svd_roi(U_master(:,:,1:n_vs), ...
        reshape(permute(fluor_allcat_predicted_reduced(:,:,:,curr_regressor_idx),[3,2,1]), ...
        n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
        n_rois,[],size(fluor_allcat_predicted,1)),[3,2,1]);
    
    curr_fluor =  fluor_roi_deriv - fluor_predicted_roi_reduced;
    curr_fluor = curr_fluor(:,:,1:10) - curr_fluor(:,:,11:end);
    
    % (if movement, align to move onset)
    if strcmp(trial_groups{curr_group},'Move onset')
        t_leeway = -t(1);
        leeway_samples = round(t_leeway*(sample_rate));
        for i = 1:size(fluor_allcat_deriv,1)
            curr_fluor(i,:,:) = circshift(curr_fluor(i,:,:),-move_idx(i)+leeway_samples,2);
        end
    end
    
    % Set trials and grouping (-1,1) to use
    switch trial_groups{curr_group}
        case 'Stim'
            use_trials = move_t < 1;
            trial_group = trial_side_allcat;
            group1 = 1;
            group2 = -1;
        case 'Move onset'
            use_trials = move_t < 1 & trial_side_allcat == 1;
            trial_group = trial_choice_allcat;
            group1 = -1;
            group2 = 1;
        case 'Reward'
            use_trials = move_t < 1;
            trial_group = trial_choice_allcat == -trial_side_allcat;
            group1 = 1;
            group2 = 0;
    end
    
    act_rank = tiedrank(curr_fluor(use_trials,:,:));
    act_rank_difference = squeeze(( ...
        nanmean(act_rank(trial_group(use_trials) == group1,:,:),1) - ...
        nanmean(act_rank(trial_group(use_trials) == group2,:,:),1))./sum(use_trials));   
    
    p1 = subplot(1,length(trial_groups),curr_group); 
    hold on; set(gca,'ColorOrder',jet(size(curr_fluor,3)));
    plot(t,act_rank_difference,'linewidth',2);
    axis tight; line([0,0],ylim,'color','k');
    xlabel('Time');
    ylabel('Mean rank difference');
    title([trial_groups{curr_group}]);
   
    
    
%     %%% (testing shuffle)    
%     n_shuff = 1000;
%     rank_diff_shuff = nan(size(act_rank,2),size(act_rank,3),n_shuff);
%     for curr_shuff = 1:n_shuff
%         trial_group_shuff = AP_shake(trial_group(use_trials));
%         rank_diff_shuff(:,:,curr_shuff) = squeeze(( ...
%             nanmean(act_rank(trial_group_shuff == group1,:,:),1) - ...
%             nanmean(act_rank(trial_group_shuff == group2,:,:),1))./sum(use_trials));  
%         AP_print_progress_fraction(curr_shuff,n_shuff);
%     end
%     
%     rank_diff_rank = permute(tiedrank(permute(cat(3,act_rank_difference,rank_diff_shuff),[3,1,2])),[2,3,1]);
%     rank_diff_p = rank_diff_rank(:,:,1)./size(rank_diff_rank,3);
%     
%     max_plot = max(act_rank_difference(:));
%     spacing = max_plot/10;
%     
%     rank_diff_p_plot = +(rank_diff_p > 0.95);
%     rank_diff_p_plot(~rank_diff_p_plot) = NaN;
%     rank_diff_p_plot = +rank_diff_p_plot.*(1:n_rois)*spacing+max_plot;  
%     plot(p1,t,rank_diff_p_plot,'.','MarkerSize',5);
%     axis tight;
    
   
    
    
end





% Get reduced-model fluorescence
% fluor_allcat_deconv_reduced_stim = fluor_allcat_deconv - fluor_allcat_predicted_reduced(:,:,:,1);

fluor_allcat_deconv_move = fluor_allcat_deconv;
% fluor_allcat_deconv_reduced_move = fluor_allcat_deconv - fluor_allcat_predicted_reduced(:,:,:,2);
t_leeway = 0.5;
leeway_samples = round(t_leeway*(sample_rate));
for i = 1:size(fluor_allcat_deconv,1)
    fluor_allcat_deconv_move(i,:,:) = circshift(fluor_allcat_deconv_move(i,:,:),-move_idx(i)+leeway_samples,2);
%     fluor_allcat_deconv_reduced_move(i,:,:) = circshift(fluor_allcat_deconv_reduced_move(i,:,:),-move_idx(i)+leeway_samples,2);
end

plot_rxn_time = [0,0.5];

plot_trials = trial_side_allcat == 1 & ...
    trial_contrast_allcat > 0 & ...
    trial_choice_allcat == -1 & ...
    move_t >= plot_rxn_time(1) & ...
    move_t <= plot_rxn_time(2);

plot_px = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    squeeze(nanmean(fluor_allcat_deconv(plot_trials,:,:),1))');

AP_image_scroll(plot_px,t);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k');




% Plot fluorescence mirroring and combining sides (side = ipsi/contra now)
mirror_matrix = reshape(U_master(:,:,1:n_vs),[],n_vs)'* ...
    reshape(AP_reflect_widefield(U_master(:,:,1:n_vs)),[],n_vs);
fluor_allcat_deconv_mirror = reshape(transpose( ...
    mirror_matrix*reshape(fluor_allcat_deconv,[],n_vs)'),size(fluor_allcat_deconv));
fluor_allcat_deconv_move_mirror = reshape(transpose( ...
    mirror_matrix*reshape(fluor_allcat_deconv_move,[],n_vs)'),size(fluor_allcat_deconv_move));

plot_rxn_time = [0,0.5];

plot_side = 1;
plot_choice = -1;
plot_contrast = trial_contrast_allcat > 0;
plot_move_t = move_t >= plot_rxn_time(1) & move_t <= plot_rxn_time(2);

curr_trials_standard = ...%trial_side_allcat == plot_side & ...
    trial_choice_allcat == plot_choice & ...
    plot_contrast & plot_move_t; 
curr_trials_mirror = ...%trial_side_allcat == -plot_side & ...
    trial_choice_allcat == -plot_choice & ...
    plot_contrast & plot_move_t; 

curr_px = svdFrameReconstruct(U_master(:,:,1:n_vs),squeeze(nanmean( ...
    [fluor_allcat_deconv(curr_trials_standard,:,:); ...
    fluor_allcat_deconv_mirror(curr_trials_mirror,:,:)],1))');

AP_image_scroll(curr_px,t);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k');



%% testing ctx->wheel with mirror symmetry

clear all

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);
mua_taskpred_k_all = cell(1,1);
mua_taskpred_all = cell(1,1);
mua_taskpred_reduced_all = cell(1,1);
mua_ctxpred_taskpred_k_all = cell(1,1);
mua_ctxpred_taskpred_all = cell(1,1);
mua_ctxpred_taskpred_reduced_all = cell(1,1);
wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;    
    
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
        raster_window = [-0.5,2];
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
            diff(fVdf,[],2)',time_bin_centers)';
        
        %%%% TESTING: DECONV
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        %%%%
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;        
        
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
        
        %%%%% TEST
        warning('OVERRIDING LAMBDA')
        lambda = 5;
        %%%%%%%
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        %         [~,predicted_spikes_std,explained_var] = ...
        %             AP_regresskernel(dfVdf_resample, ...
        %             binned_spikes_std,kernel_frames,lambda, ...
        %             regression_params.zs,regression_params.cvfold, ...
        %             false,regression_params.use_constant);
        
        %%%% TESTING DECONV
        [ctx_str_k,predicted_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Shove the k's into the master to save small V's 
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
                     
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        predicted_spikes = (predicted_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,predicted_spikes',t_peri_event)./raster_sample_rate;
        
        %%%%
        
        %%%% TESTING CTX->WHEEL VELOCITY
        
        % TEST: add mirror to force symmetric kernel
        mirror_matrix = reshape(Udf,[],size(Udf,3))'* ...
            reshape(AP_reflect_widefield(Udf),[],size(Udf,3));
        fVdf_deconv_resample_mirror = mirror_matrix*fVdf_deconv_resample;
            
        % Resample wheel, predict both velocity and speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample_mirror(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Shove the k's into the master to save small V's 
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%%%
        
        %%% Wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_bins(x,:)), ...
            find(trial_outcome == 1)','uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_bins(x,:)), ...
            find(trial_outcome == -1)','uni',false))) > 0;
        
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
        
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);    
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING: TASK REGRESSION (LONG TRACE)
                
        %%% Regression (separate regressors to get partial explained variance)
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times)).* ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end    
        
        % Stim move regressors (one for each stim when it starts to move)       
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));            
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
                        
        end      
        
        % Stim center regressors (one for each stim when it's stopped during reward)    
        unique_contrasts = unique(contrasts(contrasts > 0));
        stim_contrasts = ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ... 
                stim_contrasts == unique_contrasts(curr_contrast));
                        
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));                  
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
                        
        end    

        % Move onset regressors (L/R)          
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move ongoing regressors (L/R)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,wheel_velocity_interp < 0) = 1;
        move_ongoing_regressors(2,wheel_velocity_interp > 0) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        if length(signals_events.interactiveOnTimes) ~= length(move_t)
            error('Different number of interactive ons and move times')
        end
        
        go_cue_regressors = zeros(2,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        outcome_regressors = zeros(2,length(time_bin_centers));
        
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate regressors, set parameters
        regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        regressor_labels = {'Stim','Move onset','Go cue','Outcome'};     
        
        t_shifts = {[0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [-0.1,0.5]; ... % go cue
            [-0.5,1]}; % outcome
        
        sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;       
        
        % Regression task->MUA
        mua_taskpred_k = cell(length(regressors)+return_constant,n_depths);
        mua_taskpred = nan(size(event_aligned_mua));
        mua_taskpred_reduced = ...
            repmat(nan(size(event_aligned_mua)),1,1,1,length(regressors));       
        for curr_depth = 1:n_depths
                        
            baseline = nanmean(reshape(event_aligned_mua(:,t < 0,curr_depth),[],1)*raster_sample_rate);
            activity = single(binned_spikes(curr_depth,:)) - baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
                   
            [task_kernel,activity_predicted,expl_var,activity_predicted_reduced] = ...
                AP_regresskernel(regressors,activity,sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            mua_taskpred_k(:,curr_depth) = task_kernel;
            
            mua_taskpred(:,:,curr_depth) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event)./raster_sample_rate;
            
            mua_taskpred_reduced(:,:,curr_depth,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event)./raster_sample_rate,permute(1:length(regressors),[1,3,4,2]),'uni',false));       
            
        end
        
        % Regression task->MUA-ctxpred
        mua_ctxpred_taskpred_k = cell(length(regressors)+return_constant,n_depths);
        mua_ctxpred_taskpred = nan(size(event_aligned_mua));
        mua_ctxpred_taskpred_reduced = ...
            repmat(nan(size(event_aligned_mua)),1,1,1,length(regressors));       
        for curr_depth = 1:n_depths
                      
            baseline = nanmean(reshape(event_aligned_mua_ctxpred(:,t < 0,curr_depth),[],1)*raster_sample_rate);
            activity = single(predicted_spikes(curr_depth,:))-baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
            
            % (there are NaNs because of skipped edges)
            nan_samples = isnan(activity);
            activity_predicted = nan(size(activity));
            activity_predicted_reduced = nan(size(activity,1), ...
                size(activity,2),length(regressors));
            
            [task_kernel,activity_predicted(~nan_samples), ...
                expl_var,activity_predicted_reduced(:,~nan_samples,:)] = ...
                AP_regresskernel( ...
                cellfun(@(x) x(:,~nan_samples),regressors,'uni',false), ...
                activity(~nan_samples),sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            mua_ctxpred_taskpred_k(:,curr_depth) = task_kernel;
            
            mua_ctxpred_taskpred(:,:,curr_depth) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event)./raster_sample_rate;
            
            mua_ctxpred_taskpred_reduced(:,:,curr_depth,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event)./raster_sample_rate,permute(1:length(regressors),[1,3,4,2]),'uni',false));       
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        mua_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_taskpred_k;
        mua_taskpred_all{curr_animal,1}{curr_day,1} = mua_taskpred(use_trials,:,:,:);
        mua_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_taskpred_reduced(use_trials,:,:,:);
        
        mua_ctxpred_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_k;
        mua_ctxpred_taskpred_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred(use_trials,:,:,:);
        mua_ctxpred_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_reduced(use_trials,:,:,:);
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);   
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except curr_animal ...
        n_aligned_depths regression_params animals t ...
        fluor_all ...
        mua_all ...
        mua_ctxpred_all ...
        mua_taskpred_k_all ...
        mua_taskpred_all ...
        mua_taskpred_reduced_all ...
        mua_ctxpred_taskpred_k_all ...
        mua_ctxpred_taskpred_all ...
        mua_ctxpred_taskpred_reduced_all ...
        wheel_ctxpred_all ...
        ctx_str_k_all ...
        ctx_wheel_k_all ...
        wheel_all ...
        outcome_all ...
        D_all
    
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    outcome_all ...
    D_all

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_SYMWHEELTEST'];
save([save_path filesep save_fn],'-v7.3');






%% temp plot for bilab

% figure; hold on
% 
% AP_stackplot(squeeze(nanmean(mua_allcat,1)),t,5,false,'k',1:4);
% 
% AP_stackplot(squeeze(nanmean(mua_allcat(trial_contrastside_allcat > 0,:,:),1)),t,5,false,'r',1:4);
% AP_stackplot(squeeze(nanmean(mua_allcat(trial_contrastside_allcat < 0,:,:),1)),t,5,false,'b',1:4);
% 
% AP_stackplot(squeeze(nanmean(mua_allcat(trial_choice_allcat == -1,:,:),1)),t,5,false,[0.6,0,0.6],1:4);
% AP_stackplot(squeeze(nanmean(mua_allcat(trial_choice_allcat == 1,:,:),1)),t,5,false,[0,0.7,0],1:4);
% 
% AP_stackplot(squeeze(nanmean(mua_allcat(move_t > 0 & move_t < 0.2,:,:),1)),t,5,false,[0,0.3,0.3],1:4);
% AP_stackplot(squeeze(nanmean(mua_allcat(move_t > 0.2 & move_t < 0.4,:,:),1)),t,5,false,[0,0.5,0.5],1:4);
% AP_stackplot(squeeze(nanmean(mua_allcat(move_t > 0.4 & move_t < 0.6,:,:),1)),t,5,false,[0,0.7,0.7],1:4);
% 
% AP_stackplot(squeeze(nanmean(mua_allcat(trial_outcome_allcat == -1,:,:),1)),t,5,false,'m',1:4);
% AP_stackplot(squeeze(nanmean(mua_allcat(trial_outcome_allcat == 1,:,:),1)),t,5,false,'c',1:4);



% figure; hold on
% AP_stackplot(squeeze(nanmean(mua_allcat,1)),t,5,false,'k',1:4);
% AP_stackplot(squeeze(nanmean(mua_ctxpred_allcat,1)),t,5,false,[0,0.7,0],1:4);



t_shifts = {[0,0.5]; ... % stim
    [-0.5,1]; ... % move
    [-0.1,0.5]; ... % go cue
    [-0.5,1]}; % outcome

sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
    round(x(2)*(sample_rate)),t_shifts,'uni',false);

curr_regressor = 4;
curr_k_v = permute(cell2mat(permute(fluor_kernel(curr_regressor,1:n_vs),[1,3,2])),[3,2,1]);
curr_k_px = reshape(svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    reshape(curr_k_v,n_vs,[])),size(U_master,1),size(U_master,2),[],size(curr_k_v,3));



% figure; 
% plot_depth = 1;
% plot_reduction = 1;
% 
% p1 = subplot(1,3,1); hold on;
% plot(t,nanmean(mua_allcat(move_t < 0.5 & trial_side_allcat == 1,:,plot_depth),1));
% plot(t,nanmean(mua_allcat(move_t < 0.5 & trial_side_allcat == -1,:,plot_depth),1));
% 
% p2 = subplot(1,3,2); hold on;
% plot(t,nanmean(mua_taskpred_reduced_allcat(move_t < 0.5 & trial_side_allcat == 1,:,plot_depth,plot_reduction),1));
% plot(t,nanmean(mua_taskpred_reduced_allcat(move_t < 0.5 & trial_side_allcat == -1,:,plot_depth,plot_reduction),1));
% 
% p3 = subplot(1,3,3); hold on;
% plot(t,nanmean(mua_allcat(move_t < 0.5 & trial_side_allcat == 1,:,plot_depth),1) - ...
%     nanmean(mua_taskpred_reduced_allcat(move_t < 0.5 & trial_side_allcat == 1,:,plot_depth,plot_reduction),1));
% plot(t,nanmean(mua_allcat(move_t < 0.5 & trial_side_allcat == -1,:,1),1) - ...
%     nanmean(mua_taskpred_reduced_allcat(move_t < 0.5 & trial_side_allcat == -1,:,plot_depth,plot_reduction),1));
% 
% linkaxes([p1,p2,p3]);



%% quick test - apparently the ROIs can distinguish L/R too?? ugh


% Load widefield kernel ROIs
kernel_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\kernel_roi';
load(kernel_roi_fn);

kernel_roi_bw_hemidiff = kernel_roi.bw - AP_reflect_widefield(kernel_roi.bw);

% Get predicted fluorescence in ROIs
fluor_kernel_roi_bw = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],kernel_roi.bw), ...
    size(kernel_roi.bw,3),[],size(fluor_allcat_deconv,1)),[3,2,1]);

fluor_kernel_roi_bw_hemidiff = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],kernel_roi_bw_hemidiff), ...
    size(kernel_roi_bw_hemidiff,3),[],size(fluor_allcat_deconv,1)),[3,2,1]);

fluor_kernel_roi_weighted = permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],kernel_roi.max_weighted), ...
    size(kernel_roi.max_weighted,3),[],size(fluor_allcat_deconv,1)),[3,2,1]);




fluor_kernel_roi_bw_reduced = arrayfun(@(x) ...
    permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_predicted_reduced(:,:,:,x),[3,2,1]), ...
    n_vs,[]),[],[],kernel_roi.bw), ...
    size(kernel_roi.bw,3),[],size(fluor_allcat_predicted,1)),[3,2,1]), ...
    1:size(fluor_allcat_predicted_reduced,4),'uni',false);
fluor_kernel_roi_bw_reduced = cat(4,fluor_kernel_roi_bw_reduced{:});


fluor_roi_deconv_reduced = arrayfun(@(x) ...
    permute(reshape( ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_predicted_reduced(:,:,:,x),[3,2,1]), ...
    n_vs,[]),[],[],cat(3,wf_roi.mask)), ...
    size(cat(3,wf_roi.mask),3),[],size(fluor_allcat_predicted,1)),[3,2,1]), ...
    1:size(fluor_allcat_predicted_reduced,4),'uni',false);
fluor_roi_deconv_reduced = cat(4,fluor_roi_deconv_reduced{:});



% Rank differences by experiment
trials_animal = arrayfun(@(x) size(vertcat(mua_all{x}{:}),1),1:size(mua_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(mua_all{:}));
use_split = trials_recording;

mua_exp = mat2cell(fluor_roi_deconv,use_split,length(t),n_depths);
mua_ctxpred_exp = mat2cell(fluor_roi_deconv,use_split,length(t),n_depths);

mua_taskpred_reduced_exp = mat2cell(fluor_roi_deconv_reduced,use_split,length(t),n_depths, ...
    size(mua_taskpred_reduced_allcat,4));
mua_ctxpred_taskpred_reduced_exp = mat2cell(fluor_roi_deconv_reduced,use_split,length(t),n_depths, ...
    size(mua_ctxpred_taskpred_reduced_allcat,4));

move_t_exp = mat2cell(move_t,use_split,1);
move_idx_exp = mat2cell(move_idx,use_split,1);
outcome_idx_exp = mat2cell(outcome_idx,use_split,1);
trial_side_allcat_exp = mat2cell(trial_side_allcat,use_split,1);
trial_contrast_allcat_exp = mat2cell(trial_contrast_allcat,use_split,1);
trial_choice_allcat_exp = mat2cell(trial_choice_allcat,use_split,1);

regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
trial_groups = {'Stim','Move onset','Outcome'};
t_groups = {[0,0.3],[-0.2,0.2],[0,0.2]};

act_rank_difference = nan(length(t),n_depths,length(use_split),length(trial_groups));
predicted_act_rank_difference = nan(length(t),n_depths,length(use_split),length(trial_groups));

act_rank_difference_trial = nan(n_depths,length(use_split),length(trial_groups));
predicted_act_rank_difference_trial = nan(n_depths,length(use_split),length(trial_groups));

for curr_exp = 1:length(mua_exp)
    for curr_group = 1:length(trial_groups)
        
        % Get MUA/predicted using reduced model
        curr_regressor_idx = strcmp(trial_groups{curr_group},regressor_labels);
        curr_mua =  mua_exp{curr_exp} - mua_taskpred_reduced_exp{curr_exp}(:,:,:,curr_regressor_idx);
        curr_mua_ctx = mua_ctxpred_exp{curr_exp} - mua_ctxpred_taskpred_reduced_exp{curr_exp}(:,:,:,curr_regressor_idx);
        
        % Set common NaNs
        nan_samples = isnan(curr_mua) | isnan(curr_mua_ctx);
        curr_mua(nan_samples) = NaN;
        curr_mua_ctx(nan_samples) = NaN;
        
        % (movement: align to move onset)
        if any(strfind(lower(trial_groups{curr_group}),'move'))
            t_leeway = -t(1);
            leeway_samples = round(t_leeway*(sample_rate));
            for i = 1:size(curr_mua,1)
                curr_mua(i,:,:) = circshift(curr_mua(i,:,:),-move_idx_exp{curr_exp}(i)+leeway_samples,2);
                curr_mua_ctx(i,:,:) = circshift(curr_mua_ctx(i,:,:),-move_idx_exp{curr_exp}(i)+leeway_samples,2);
            end
        end
        
        % (outcome: align to outcome)
        if any(strfind(lower(trial_groups{curr_group}),'outcome'))
            t_leeway = -t(1);
            leeway_samples = round(t_leeway*(sample_rate));
            for i = 1:size(curr_mua,1)
                curr_mua(i,:,:) = circshift(curr_mua(i,:,:),-outcome_idx_exp{curr_exp}(i)+leeway_samples,2);
                curr_mua_ctx(i,:,:) = circshift(curr_mua_ctx(i,:,:),-outcome_idx_exp{curr_exp}(i)+leeway_samples,2);
            end
        end
        
        % Set trials and grouping to use
        switch trial_groups{curr_group}
            case 'Stim'
                use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} < 0.5 & trial_contrast_allcat_exp{curr_exp} > 0;
                trial_group = trial_side_allcat_exp{curr_exp};
                group1 = 1;
                group2 = -1;
            case 'Move onset'
                use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} < 0.5;
                trial_group = trial_choice_allcat_exp{curr_exp};
                group1 = -1;
                group2 = 1;                
            case 'Outcome'
                use_trials = move_t_exp{curr_exp} > 0 & move_t_exp{curr_exp} < 0.5;
                trial_group = trial_choice_allcat_exp{curr_exp} == -trial_side_allcat_exp{curr_exp};
                group1 = 1;
                group2 = 0;
        end
        
        act_rank = tiedrank(curr_mua(use_trials,:,:));
        predicted_act_rank = tiedrank(curr_mua_ctx(use_trials,:,:));
        
        act_rank_difference(:,:,curr_exp,curr_group) = squeeze(( ...
            nanmean(act_rank(trial_group(use_trials) == group1,:,:),1) - ...
            nanmean(act_rank(trial_group(use_trials) == group2,:,:),1))./max(act_rank,[],1));
        predicted_act_rank_difference(:,:,curr_exp,curr_group) = squeeze(( ...
            nanmean(predicted_act_rank(trial_group(use_trials) == group1,:,:),1) - ...
            nanmean(predicted_act_rank(trial_group(use_trials) == group2,:,:),1))./max(predicted_act_rank,[],1));
        
        use_t = t >= t_groups{curr_group}(1) & t <= t_groups{curr_group}(2);
        act_rank_trial = tiedrank(squeeze(nanmean(curr_mua(use_trials,use_t,:),2)));
        predicted_act_rank_trial = tiedrank(squeeze(nanmean(curr_mua_ctx(use_trials,use_t,:),2)));
        
        act_rank_difference_trial(:,curr_exp,curr_group) = ...
            (nanmean(act_rank_trial(trial_group(use_trials) == group1,:),1) - ...
            nanmean(act_rank_trial(trial_group(use_trials) == group2,:),1))./max(act_rank_trial,[],1);
        predicted_act_rank_difference_trial(:,curr_exp,curr_group) = ...
            (nanmean(predicted_act_rank_trial(trial_group(use_trials) == group1,:),1) - ...
            nanmean(predicted_act_rank_trial(trial_group(use_trials) == group2,:),1))./max(predicted_act_rank_trial,[],1);
        
%         n_shuff = 100;
%         shuff_rank_difference = nan(length(t),n_depths,n_shuff);
%         shuff_predicted_rank_difference = nan(length(t),n_depths,n_shuff);
%         for curr_shuff = 1:n_shuff
%             curr_shake_group = AP_shake(trial_group(use_trials));
%             shuff_rank_difference(:,:,curr_shuff) = squeeze(( ...
%                 nanmean(act_rank(curr_shake_group == group1,:,:),1) - ...
%                 nanmean(act_rank(curr_shake_group == group2,:,:),1))./max(act_rank,[],1));
%             shuff_predicted_rank_difference(:,:,curr_shuff) = squeeze(( ...
%                 nanmean(predicted_act_rank(curr_shake_group == group1,:,:),1) - ...
%                 nanmean(predicted_act_rank(curr_shake_group == group2,:,:),1))./max(predicted_act_rank,[],1));
%             AP_print_progress_fraction(curr_shuff,n_shuff);
%         end
        
    end
end

act_rank_difference_mean = squeeze(nanmean(act_rank_difference,3));
predicted_act_rank_difference_mean = squeeze(nanmean(predicted_act_rank_difference,3));

figure;
for curr_group = 1:length(trial_groups)
    p1 = subplot(2,length(trial_groups),curr_group);
    hold on; set(gca,'ColorOrder',copper(n_depths));
    plot(t,act_rank_difference_mean(:,:,curr_group),'linewidth',2);
    axis tight; line([0,0],ylim,'color','k');
    xlabel('Time');
    ylabel('Mean rank difference');
    title([trial_groups{curr_group} ': Measured']);
    
    p2 = subplot(2,length(trial_groups),length(trial_groups)+curr_group);
    hold on; set(gca,'ColorOrder',copper(n_depths));
    plot(t,predicted_act_rank_difference_mean(:,:,curr_group),'linewidth',2);
    axis tight; line([0,0],ylim,'color','k');
    xlabel('Time');
    ylabel('Mean rank difference');
    title([trial_groups{curr_group} ': Predicted']);
    linkaxes([p1,p2],'xy');
end

figure;
p1 = subplot(1,2,1);
errorbar(squeeze(nanmean(act_rank_difference_trial,2)), ...
    squeeze(nanstd(act_rank_difference_trial,[],2)./ ...
    sqrt(sum(~isnan(act_rank_difference_trial),2))),'linewidth',2);
xlabel('Striatum depth');
ylabel('Rank difference');
legend(trial_groups);
title('Striatum');

p2 = subplot(1,2,2);
errorbar(squeeze(nanmean(predicted_act_rank_difference_trial,2)), ...
    squeeze(nanstd(predicted_act_rank_difference_trial,[],2)./ ...
    sqrt(sum(~isnan(predicted_act_rank_difference_trial),2))),'linewidth',2);
xlabel('Striatum depth');
ylabel('Rank difference');
legend(trial_groups);
title('Cortex-predicted striatum');

linkaxes([p1,p2]);


























