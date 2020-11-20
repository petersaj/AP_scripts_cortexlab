%% Testing for corticostriatal-specific imaging

%% ~~~~~~~~~ TESTS ~~~~~~~~~

%% Get explained variance in sliding window

% Load data
animal = 'AP068';
day = '2020-03-17';
experiments = AP_list_experiments(animal,day);
% experiment = experiments(find(contains({experiments.protocol},'vanillaChoiceworld'),1,'first')).experiment;
experiment = experiments(find(contains({experiments.protocol},'AP_lcrGratingPassive'),1,'last')).experiment;
verbose = true;
AP_load_experiment;


% Get striatum MUA in sliding window
skip_seconds = 60;
spike_binning_t = 1/framerate; % seconds
spike_binning_t_edges = frame_t(1)+skip_seconds:spike_binning_t:frame_t(end)-skip_seconds;
spike_binning_t_centers = spike_binning_t_edges(1:end-1) + diff(spike_binning_t_edges)/2;

mua_window = 500; % MUA window in microns
mua_window_spacing = 250; % MUA window spacing in microns

mua_bins = [str_depth(1):mua_window_spacing:(str_depth(2)-mua_window); ...
    (str_depth(1):mua_window_spacing:(str_depth(2)-mua_window))+mua_window];
mua_bin_centers = mua_bins(1,:) + diff(mua_bins,[],1)/2;
n_depths = length(mua_bin_centers);

striatum_mua = zeros(size(mua_bins,2),length(spike_binning_t_edges)-1);
for curr_depth = 1:size(mua_bins,2)
    curr_depth_templates_idx = ...
        find(template_depths >= mua_bins(1,curr_depth) & ...
        template_depths < mua_bins(2,curr_depth));
    
    striatum_mua(curr_depth,:) = histcounts(spike_times_timeline( ...
        ismember(spike_templates,curr_depth_templates_idx)),spike_binning_t_edges);
end

striatum_mua_std = striatum_mua./nanstd(striatum_mua,[],2);
striatum_mua_std(isnan(striatum_mua_std)) = 0;



use_svs = 1:200;
kernel_t = [-0.5,0.5];
kernel_frames = round(kernel_t(1)*framerate):round(kernel_t(2)*framerate);
lambda = 10;
zs = [false,false];
cvfold = 5;
return_constant = false;
use_constant = true;
        
fVdf_deconv = AP_deconv_wf(fVdf);
fVdf_deconv(isnan(fVdf_deconv)) = 0;
fVdf_deconv_resample = interp1(frame_t,fVdf_deconv(use_svs,:)',spike_binning_t_centers)';
        
[k,striatum_mua_std_predicted,explained_var] = ...
    AP_regresskernel(fVdf_deconv_resample, ...
    striatum_mua_std,kernel_frames,lambda,zs,cvfold,return_constant,use_constant);

% Convert kernel to pixel space
r_px = zeros(size(Udf,1),size(Udf,2),size(k,2),size(k,3),'single');
for curr_spikes = 1:size(k,3)
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(Udf(:,:,use_svs),k(:,:,curr_spikes));
end

AP_image_scroll(r_px,kernel_frames/framerate);
caxis([-prctile(r_px(:),99.9),prctile(r_px(:),99.9)])
colormap(colormap_BlueWhiteRed);
axis image;


% Get center of mass for each pixel
% (get max r for each pixel, filter out big ones)
r_px_max = squeeze(max(r_px,[],3));
r_px_max(isnan(r_px_max)) = 0;
r_px_max_norm = bsxfun(@rdivide,r_px_max, ...
    permute(max(reshape(r_px_max,[],n_depths),[],1),[1,3,2]));
r_px_max_norm(isnan(r_px_max_norm)) = 0;
r_px_com = sum(bsxfun(@times,r_px_max_norm,permute(1:n_depths,[1,3,2])),3)./sum(r_px_max_norm,3);

% Plot map of cortical pixel by preferred depth of probe
r_px_com_col = ind2rgb(round(mat2gray(r_px_com,[1,n_depths])*255),jet(255));
figure;
a1 = axes('YDir','reverse');
imagesc(avg_im); colormap(gray); caxis([0,prctile(avg_im(:),99)]);
axis off; axis image;
a2 = axes('Visible','off'); 
p = imagesc(r_px_com_col);
axis off; axis image;
set(p,'AlphaData',mat2gray(max(r_px_max_norm,[],3), ...
     [0,double(prctile(reshape(max(r_px_max_norm,[],3),[],1),95))]));
set(gcf,'color','w');

c1 = colorbar('peer',a1,'Visible','off');
c2 = colorbar('peer',a2);
ylabel(c2,'Depth (\mum)');
colormap(c2,jet);
set(c2,'YDir','reverse');
set(c2,'YTick',linspace(0,1,n_depths));
set(c2,'YTickLabel',round(linspace(mua_bins(1),mua_bins(end),n_depths)));

% Plot MUA correlation and explained variance by cortex
figure;
subplot(4,1,1); 
plot(mua_bin_centers,explained_var.total,'k','linewidth',2);
line(xlim,[0,0],'color','r');

subplot(4,1,2:4);
imagesc(mua_bin_centers,mua_bin_centers,corrcoef(striatum_mua_std'));
colormap(hot);
caxis([0,1]);


% Plot measured and predicted stim-aligned activity

% Set options
surround_window = [-0.5,3];
baseline_window = [-0.1,0];

surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);
baseline_surround_time = baseline_window(1):surround_samplerate:baseline_window(2);

stimIDs = trial_conditions(:,1).*trial_conditions(:,2);
use_stims = find(stimIDs > 0);

use_stimOn_times = stimOn_times(use_stims);

stim_surround_times = bsxfun(@plus, use_stimOn_times(:), surround_time);
stim_baseline_surround_times = bsxfun(@plus, use_stimOn_times(:), baseline_surround_time);

psth_measured_stim = permute(interp1(spike_binning_t_centers, ...
    striatum_mua_std',stim_surround_times),[3,2,1]);
psth_predicted_stim = permute(interp1(spike_binning_t_centers, ...
    striatum_mua_std_predicted',stim_surround_times),[3,2,1]);

psth_measured_baseline = permute(nanmean(interp1(spike_binning_t_centers, ...
    striatum_mua_std',stim_baseline_surround_times),2),[3,2,1]);
psth_predicted_baseline = permute(nanmean(interp1(spike_binning_t_centers, ...
    striatum_mua_std_predicted',stim_baseline_surround_times),2),[3,2,1]);

psth_measured = nanmean(psth_measured_stim - psth_measured_baseline,3);
psth_predicted = nanmean(psth_predicted_stim - psth_predicted_baseline,3);

figure; hold on;
AP_stackplot(psth_measured',surround_time,2,false,'k',mua_bin_centers);
AP_stackplot(psth_predicted',surround_time,2,false,[0,0.7,0]);
line([0,0],ylim,'color','k','linestyle','--');

%% ~~~~~~~~~ ALIGNMENT ~~~~~~~~~

%% Get and save average response to right grating (for animal alignment)

animals = {'AP085','AP086','AP087'};

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    % (get all passive gratings)
    protocol = 'AP_lcrGratingPassive';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    experiments = experiments([experiments.imaging]);
    
%     % (get all imaging days with behavior)
%     bhv_protocol = 'choiceworld';
%     bhv_experiments = AP_find_experiments(animal,bhv_protocol,true);
%     bhv_experiments = bhv_experiments([bhv_experiments.imaging]);
%     
%     % Use only post-trained (larger response, match to trained passive)
%     trained_experiments = ismember({experiments.day},{bhv_experiments.day});
%     experiments = experiments(trained_experiments);
    
    im_stim_all = cell(length(experiments),1);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        AP_load_experiment;
        
        Udf_aligned = AP_align_widefield(Udf,animal,day,'day_only');
        
        % Set options
        surround_window = [-0.5,3];
        baseline_window = [-0.1,0];
        
        surround_samplerate = 1/(framerate*1);
        surround_time = surround_window(1):surround_samplerate:surround_window(2);
        baseline_surround_time = baseline_window(1):surround_samplerate:baseline_window(2);
        
        % Average (time course) responses
        use_vs = 1:size(U,3);
        
        conditions = unique(stimIDs);
        im_stim = nan(size(Udf_aligned,1),size(Udf_aligned,2),length(surround_time),length(conditions));
        for curr_condition_idx = 1:length(conditions)
            curr_condition = conditions(curr_condition_idx);
            
            use_stims = find(stimIDs == curr_condition);
            use_stimOn_times = stimOn_times(use_stims(2:end));
            use_stimOn_times([1,end]) = [];
            
            stim_surround_times = bsxfun(@plus, use_stimOn_times(:), surround_time);
            stim_baseline_surround_times = bsxfun(@plus, use_stimOn_times(:), baseline_surround_time);
            
            peri_stim_v = permute(interp1(frame_t,fVdf',stim_surround_times),[3,2,1]);
            baseline_v = permute(nanmean(interp1(frame_t,fVdf',stim_baseline_surround_times),2),[3,2,1]);
            
            stim_v_mean = nanmean(peri_stim_v - baseline_v,3);
            
            im_stim(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf_aligned(:,:,use_vs), ...
                stim_v_mean(use_vs,:));
        end        
        
        use_stim = 3;
        im_stim_all{curr_day} = im_stim(:,:,:,use_stim);
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
    
    % Get average stim response across days
    stim_time = [0.2,0.4];
    use_time = surround_time >= stim_time(1) & surround_time <= stim_time(2);
    right_grating = nanmean(cell2mat(reshape(cellfun(@(x) ...
        nanmean(x(:,:,use_time),3),im_stim_all,'uni',false),1,1,[])),3);
   
    % Save average stim response
    save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment\right_grating';
    save_fn = [save_path filesep animal '_right_grating.mat'];
    save(save_fn,'right_grating');
    
end


%% Align corticostriatal across animals (using right grating response)

% % Make master from passive grating response
% % (only run once)
% data_fn = 'trial_activity_AP_choiceWorldStimPassive_trained';
% AP_load_concat_normalize_ctx_str;
% 
% quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > 0.025,2);
% 
% stim_time = [0.2,0.4];
% use_time = t >= stim_time(1) & t <= stim_time(2);
% right_grating_master = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
%     squeeze(nanmean(nanmean(fluor_allcat(trial_stim_allcat == 1 & quiescent_trials,use_time,:),2),1)));
% 
% right_grating_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment\right_grating';
% right_grating_master_fn = [right_grating_path filesep 'right_grating_master.mat'];
% save(right_grating_master_fn,'right_grating_master');

% Load right grating master
right_grating_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment\right_grating';
right_grating_master_fn = [right_grating_path filesep 'right_grating_master.mat'];
load(right_grating_master_fn);

% Align animal to master
animal = 'AP087';

right_grating_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment\right_grating';
right_grating_fn = [right_grating_path filesep animal '_right_grating.mat'];
load(right_grating_fn);

AP_align_widefield(right_grating,animal,[],'new_animal',right_grating_master);


%% ~~~~~~~~~ BATCH ANALYSIS ~~~~~~~~~

%% Get passive responses before/after training

animals = {'AP063','AP064','AP066','AP068','AP071'};

im_stim_all = cell(length(animals),2);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    % (get all passive gratings)
    protocol = 'AP_lcrGratingPassive';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    experiments = experiments([experiments.imaging]);
    
    % (get all imaging days with behavior)
    bhv_protocol = 'choiceworld';
    bhv_experiments = AP_find_experiments(animal,bhv_protocol,true);
    bhv_experiments = bhv_experiments([bhv_experiments.imaging]);
    
    % Split days by naive/trained
    naive_experiments = ~ismember({experiments.day},{bhv_experiments.day});
    trained_experiments = ismember({experiments.day},{bhv_experiments.day});
    
    disp(animal);
    for curr_training = 1:2
        
        switch curr_training
            case 1
                % Naive
                curr_experiments = experiments(naive_experiments);
            case 2
                % Trained
                curr_experiments = experiments(trained_experiments);
        end
        
        for curr_day = 1:length(curr_experiments)
            
            day = curr_experiments(curr_day).day;
            experiment = curr_experiments(curr_day).experiment(end);
            load_parts.imaging = true;
            AP_load_experiment;
            
            Udf_aligned = AP_align_widefield(Udf,animal,day);
            
            % Set options
            surround_window = [-0.5,3];
            baseline_window = [-0.1,0];
            
            surround_samplerate = 1/(framerate*1);
            surround_time = surround_window(1):surround_samplerate:surround_window(2);
            baseline_surround_time = baseline_window(1):surround_samplerate:baseline_window(2);
            
            % Get wheel movements during stim, only use quiescent trials
            wheel_window = [0,0.5];
            wheel_window_t = wheel_window(1):1/framerate:wheel_window(2);
            wheel_window_t_peri_event = bsxfun(@plus,stimOn_times,wheel_window_t);
            event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
                wheel_velocity,wheel_window_t_peri_event);
            wheel_thresh = 0.025;
            quiescent_trials = ~any(abs(event_aligned_wheel) > wheel_thresh,2);
                     
            % Average (time course) responses
            use_vs = 1:size(U,3);
            
            conditions = unique(stimIDs);
            im_stim = nan(size(Udf_aligned,1),size(Udf_aligned,2),length(surround_time),length(conditions));
            for curr_condition_idx = 1:length(conditions)
                curr_condition = conditions(curr_condition_idx);
                
                use_stims = find(stimIDs == curr_condition & quiescent_trials);
                use_stimOn_times = stimOn_times(use_stims(2:end));
                use_stimOn_times([1,end]) = [];
                
                stim_surround_times = bsxfun(@plus, use_stimOn_times(:), surround_time);
                stim_baseline_surround_times = bsxfun(@plus, use_stimOn_times(:), baseline_surround_time);
                
                peri_stim_v = permute(interp1(frame_t,fVdf',stim_surround_times),[3,2,1]);
                baseline_v = permute(nanmean(interp1(frame_t,fVdf',stim_baseline_surround_times),2),[3,2,1]);
                
                stim_v_mean = nanmean(peri_stim_v - baseline_v,3);
                
                im_stim(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf_aligned(:,:,use_vs), ...
                    stim_v_mean(use_vs,:));
            end
            
            im_stim_all{curr_animal,curr_training}{curr_day} = im_stim;
            AP_print_progress_fraction(curr_day,length(curr_experiments));
        end       
    end        
end

% Get average pre/post for each animal
im_stim_avg = cellfun(@(x) nanmean(cat(5,x{:}),5),im_stim_all,'uni',false);

AP_image_scroll(im_stim_avg{1,2});
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));

% Get average stim response in time and plot
use_stim = 3;
use_t = [0.2,0.4];
use_frames = surround_time >= use_t(1) & surround_time <= use_t(2);
im_stim_avg_t = cellfun(@(x) nanmean(x(:,:,use_frames,use_stim),3),im_stim_avg,'uni',false);
figure;
for curr_animal = 1:size(im_stim_avg_t,1)
    for curr_training = 1:2
        subplot(size(im_stim_avg_t,1),3,(curr_animal-1)*3+curr_training);
        imagesc(im_stim_avg_t{curr_animal,curr_training});
        axis image off;
        caxis([-max(abs(caxis)),max(abs(caxis))]);
        colormap(brewermap([],'*RdBu'));
    end
    subplot(size(im_stim_avg_t,1),3,(curr_animal-1)*3+3);
    imagesc(im_stim_avg_t{curr_animal,2} - ...
        im_stim_avg_t{curr_animal,1});
    axis image off;
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'*RdBu'));
end

stim_naive_avg = nanmean(cat(3,im_stim_avg_t{:,1}),3);
stim_trained_avg = nanmean(cat(3,im_stim_avg_t{:,2}),3);
figure;
subplot(1,3,1);
imagesc(stim_naive_avg);
axis image off;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k');
title('Naive');

subplot(1,3,2);
imagesc(stim_trained_avg);
axis image off;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k');
title('Trained');

subplot(1,3,3);
imagesc(stim_trained_avg - stim_naive_avg);
axis image off;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k');
title('Difference');


%% Grab and save corticostriatal choiceworld trial data (widefield-only days)

clear all
disp('Choiceworld trial activity (widefield-only days)')

animals = {'AP063','AP064','AP066','AP068','AP071'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & ~[experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        AP_load_experiment;
        
        % Pull out trial data
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
        trial_data_all.task_regressor_labels = task_regressor_labels;
        trial_data_all.task_regressor_sample_shifts = task_regressor_sample_shifts;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\corticostriatal\data';
save_fn = ['trial_activity_choiceworld_corticostriatal_wfonly'];
save([save_path filesep save_fn],'-v7.3');

%% Grab and save corticostriatal passive trial data (new mice, widefield-only days)

clear all
disp('Choiceworld trial activity (widefield-only days)')

animals = {'AP077','AP079','AP085','AP086','AP087'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & ~[experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        AP_load_experiment;
        
        % Pull out trial data
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
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\corticostriatal\data';
save_fn = ['trial_activity_AP_lcrGratingPassive_corticostriatal_wfonly'];
save([save_path filesep save_fn],'-v7.3');


%% Load corticostriatal trial data

% Task
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\corticostriatal\data';

% data_fn = 'trial_activity_choiceworld_corticostriatal_wfonly';
data_fn = 'trial_activity_AP_lcrGratingPassive_corticostriatal_wfonly';

AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_recording;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));


















