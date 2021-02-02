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


%% Get widefield correlation boundaries
% (doesn't work - not enough fluorescence to have strong boundaries)

clear all
animals = {'AP063','AP064','AP066','AP068','AP071','AP085','AP086','AP087'};
protocol = 'vanillaChoiceworld';

wf_corr_borders = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    disp(animal);
    
    experiments = AP_find_experiments(animal,protocol);
    % (wf-only days: no craniotomy)
    experiments = experiments([experiments.imaging] & ~[experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = false;
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment

        % Get V covariance
        Ur = reshape(Udf, size(Udf,1)*size(Udf,2),[]); % P x S
        covV = cov(fVdf'); % S x S % this is the only one that takes some time really
        varP = dot((Ur*covV)', Ur'); % 1 x P
        
        ySize = size(Udf,1); xSize = size(Udf,2);

        px_spacing = 20;
        use_y = 1:px_spacing:size(Udf,1);
        use_x = 1:px_spacing:size(Udf,2);
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
        
        % Get downsampled correlation maps
        downsample_factor = 10;
        corr_map_downsamp = cellfun(@(x) ...
            imresize(x,1/downsample_factor,'bilinear'),corr_map,'uni',false);
                
%         % OLD: Correlation map edge detection
%         corr_map_cat = cat(3,corr_map{:});
%         corr_map_edge = corr_map_cat-imgaussfilt(corr_map_cat,20);
%         corr_map_edge_norm = reshape(zscore(reshape(corr_map_edge,[], ...
%             size(corr_map_edge,3)),[],1),size(corr_map_edge));      
%         corr_edges = nanmean(corr_map_edge,3);
        
        % NEW: Correlation map edge detection
        corr_map_cat = cat(3,corr_map{:});
        corr_map_edge = corr_map_cat-imfilter(corr_map_cat,fspecial('disk',20));
        corr_map_edge_norm = reshape(zscore(reshape(corr_map_edge,[], ...
            size(corr_map_edge,3)),[],1),size(corr_map_edge));      
        corr_edges = nanmean(corr_map_edge,3);
        % (get rid of high frequency of edges - vasculature/bone)
        corr_edges_highpass = corr_edges - imfilter(corr_edges,fspecial('disk',10));
        corr_edges_lowpass = corr_edges - corr_edges_highpass;
        
        wf_corr_borders(curr_animal).corr_map_downsamp{curr_day} = corr_map_downsamp;
        wf_corr_borders(curr_animal).corr_edges{curr_day} = corr_edges_lowpass;
        
        clearvars -except animals animal curr_animal protocol experiments curr_day animal wf_corr_borders load_parts
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
end

% (OLD - DON'T OVERWRITE, CHOOSE NEW PLACE?)
% % Save widefield correlation borders
% save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_borders'];
% save_fn = [save_path filesep 'wf_corr_borders'];
% save(save_fn,'wf_corr_borders');
% disp(['Saved ' save_fn]);



%% ~~~~~~~~~ ALIGNMENT ~~~~~~~~~

%% Get and save average response to right grating (for animal alignment)

animals = {'AP075'};

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    % (get all passive gratings)
    protocol = 'AP_lcrGratingPassive';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    experiments = experiments([experiments.imaging]);
    
    % Use only pre-trained
    bhv_protocol = 'choiceworld';
    bhv_experiments = AP_find_experiments(animal,bhv_protocol,true);
    bhv_experiments = bhv_experiments([bhv_experiments.imaging]);
    trained_experiments = ismember({experiments.day},{bhv_experiments.day});
    
    experiments = experiments(~trained_experiments);
    
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
animal = 'AP075';

right_grating_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment\right_grating';
right_grating_fn = [right_grating_path filesep animal '_right_grating.mat'];
load(right_grating_fn);

AP_align_widefield(right_grating,animal,[],'new_animal',right_grating_master);


%% Load and align grating responses (to check)

% animals = {'AP063','AP064','AP066','AP068','AP071'};
% animals = {'AP085','AP086','AP087'};

animals = {'AP063','AP064','AP066','AP068','AP071','AP075','AP077','AP079','AP085','AP086','AP087'};

% Load right grating master
right_grating_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment\right_grating';
right_grating_master_fn = [right_grating_path filesep 'right_grating_master.mat'];
load(right_grating_master_fn);

% Load and align right gratings
right_grating_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment\right_grating';
right_gratings_aligned = nan(size(right_grating_master));
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    right_grating_fn = [right_grating_path filesep animal '_right_grating.mat'];
    load(right_grating_fn);
    right_grating_aligned(:,:,curr_animal) = ...
        AP_align_widefield(right_grating,animal,[],'animal_only');
end

AP_image_scroll(cat(3,right_grating_master,right_grating_aligned),['Master',animals]);
axis image;
colormap(brewermap([],'PRGn'));
caxis([-max(abs(caxis)),max(abs(caxis))]);
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);




%% ~~~~~~~~~ SINGLE-RECORDING ANALYSIS ~~~~~~~~~



%% ~~~~~~~~~ GRAB & SAVE BATCH  ~~~~~~~~~

%% Get passive responses before/after training

clear all
disp('Getting pre/post training passive stim')

animals = {'AP063','AP064','AP066','AP071','AP068','AP085','AP086','AP087'};

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

AP_image_scroll([im_stim_avg{3,1},im_stim_avg{3,2}]);
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
    title(animals{curr_animal});
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



% Deconvolve average stim response and plot time-average
im_stim_avg_deconv = cellfun(@(x) ...
    reshape(AP_deconv_wf(reshape(x,size(x,1)*size(x,2),[])),size(x)), ...
    im_stim_avg,'uni',false);

use_stim = 3;
use_t = [0.07,0.12];
use_frames = surround_time >= use_t(1) & surround_time <= use_t(2);
im_stim_avg_deconv_t = cellfun(@(x) nanmean(x(:,:,use_frames,use_stim),3),im_stim_avg_deconv,'uni',false);
figure;
for curr_animal = 1:size(im_stim_avg_deconv_t,1)
    for curr_training = 1:2
        subplot(3,length(animals),length(animals)*(curr_training-1)+curr_animal);
        imagesc(im_stim_avg_deconv_t{curr_animal,curr_training});
        axis image off;
        caxis([-max(abs(caxis)),max(abs(caxis))]);
        colormap(brewermap([],'*RdBu'));
    end
    subplot(3,length(animals),length(animals)*2+curr_animal);
    imagesc(im_stim_avg_deconv_t{curr_animal,2} - ...
        im_stim_avg_deconv_t{curr_animal,1});
    axis image off;
    caxis([-max(abs(caxis)),max(abs(caxis))]);
    colormap(brewermap([],'*RdBu'));
    title(animals{curr_animal});
end




%% Choiceworld trial data (DCS)

clear all
disp('Choiceworld trial activity (DCS)')

animals = {'AP063','AP064','AP066','AP071'};

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
save_fn = ['trial_activity_choiceworld_cstr_dcs'];
save([save_path filesep save_fn],'-v7.3');

%% Choiceworld trial data (DMS)

clear all
disp('Choiceworld trial activity (DMS)')

animals = {'AP068','AP085','AP086','AP087'};

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
save_fn = ['trial_activity_choiceworld_cstr_dms'];
save([save_path filesep save_fn],'-v7.3');



%% ~~~~~~~~~ BATCH ANALYSIS ~~~~~~~~~


%% [[Load corticostriatal trial data]]

% Task
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\corticostriatal\data';

% data_fn = 'trial_activity_choiceworld_cstr_dcs';
data_fn = 'trial_activity_choiceworld_cstr_dms';

AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_recording;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));



%% Task > cortex kernels

% Get task>cortex parameters
n_regressors = length(task_regressor_labels);
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

% Get average task > cortex kernels (V's and ROIs)
regressor_v = cell(n_regressors,1);
for curr_regressor = 1:n_regressors
    regressor_v{curr_regressor} = nanmean(cell2mat(cellfun(@(x) x{curr_regressor}, ...
        permute(vertcat(fluor_taskpred_k_all{:}),[2,3,4,1]),'uni',false)),4);   
    AP_print_progress_fraction(curr_regressor,n_regressors);
end

% Get regressor pixels
regressor_px = cellfun(@(v) cell2mat(arrayfun(@(subregressor) ...
    svdFrameReconstruct(U_master(:,:,1:n_vs),permute(v(subregressor,:,:),[3,2,1])), ...
    permute(1:size(v,1),[1,3,4,2]),'uni',false)),regressor_v,'uni',false);

plot_regressor = 1;
AP_image_scroll(regressor_px{plot_regressor})
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
axis image;


%% Choice regression (multiple timepoints)

use_svs = 1:100;
lambda = 2;
zs = [false,false];
cvfold = 5;

use_trials = trial_stim_allcat == 0;
use_t = t < 0.5;

predict_choice = trial_choice_allcat(use_trials)';

[k,predicted_choice] = ...
    AP_regresskernel( ...
    reshape(permute(fluor_allcat_deconv(use_trials,use_t,use_svs),[2,3,1]),[],sum(use_trials)), ...
    predict_choice,0,lambda,zs,cvfold);

k_px = svdFrameReconstruct(U_master(:,:,use_svs), ...
    reshape(k,sum(use_t),length(use_svs))');

AP_image_scroll(k_px,t(use_t));
axis image;
colormap(brewermap([],'PrGn'));
caxis([-max(abs(caxis)),max(abs(caxis))]);
AP_reference_outline('ccf_aligned','k');

sign_accuracy = nanmean(sign(predicted_choice) == sign(predict_choice));
disp(['Sign accuracy: ' num2str(sign_accuracy)]);


%% Choice regression (all single timepoints)

% % Make move-aligned fluorescence
% fluor_allcat_deconv_move = fluor_allcat_deconv;
% t_leeway = -t(1);
% leeway_samples = round(t_leeway*(sample_rate));
% for i = 1:size(fluor_allcat_deconv,1)
%     fluor_allcat_deconv_move(i,:,:,:) = circshift(fluor_allcat_deconv_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
% end

use_svs = 1:100;
lambda = 2;
zs = [false,false];
cvfold = 5;

use_trials = trial_stim_allcat == 0 & move_t > 0.5;
predict_choice = AP_shake(trial_choice_allcat(use_trials)');
    
k_choice = nan(length(use_svs),length(t));
sign_accuracy = nan(size(t));
for curr_t = 1:length(t)
      
    [k_choice(:,curr_t),predicted_choice] = ...
        AP_regresskernel( ...
        reshape(permute(fluor_allcat_deconv(use_trials,curr_t,use_svs),[2,3,1]),[],sum(use_trials)), ...
        predict_choice,0,lambda,zs,cvfold);
    
    sign_accuracy(curr_t) = nanmean(sign(predicted_choice) == sign(predict_choice));
    AP_print_progress_fraction(curr_t,length(t));
    
end

k_px = svdFrameReconstruct(U_master(:,:,use_svs),k_choice);
AP_image_scroll(k_px,t);
axis image;
colormap(brewermap([],'PrGn'));
caxis([-max(abs(caxis)),max(abs(caxis))]);
AP_reference_outline('ccf_aligned','k');

figure;
plot(t,sign_accuracy,'k');
xlabel('Time');
ylabel('Sign accuracy');


%% Choice regression (ROI: all single timepoints)

guide_im = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    squeeze(nanmean(nanmean(fluor_allcat_deconv(trial_stim_allcat == 1,21:23,:),1),2)));
[roi_trace_long,roi_mask] = AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),guide_im);

roi_trace = reshape(roi_trace_long,length(t),[])';

lambda = 0;
zs = [false,false];
cvfold = 5;

use_trials = trial_stim_allcat == 0 & move_t > 0.5;
predict_choice = trial_choice_allcat(use_trials)';
    
k_choice = nan(size(t));
sign_accuracy = nan(size(t));
for curr_t = 1:length(t)
      
    [k_choice(:,curr_t),predicted_choice] = ...
        AP_regresskernel( ...
        roi_trace(use_trials,curr_t)', ...
        predict_choice,0,lambda,zs,cvfold);
    
    sign_accuracy(curr_t) = nanmean(sign(predicted_choice) == sign(predict_choice));
    AP_print_progress_fraction(curr_t,length(t));
    
end

figure;
plot(t,sign_accuracy,'k');
xlabel('Time');
ylabel('Sign accuracy');

%% ROI stim activity by choice

guide_im = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    squeeze(nanmean(nanmean(fluor_allcat_deconv(trial_stim_allcat == 1,21:23,:),1),2)));
[roi_trace_long,roi_mask] = AP_svd_roi(U_master(:,:,1:n_vs), ...
    reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),guide_im);
roi_trace = reshape(roi_trace_long,length(t),[])';

[roi_trace_grp,grp_name] = grpstats(roi_trace, ...
    [trial_stim_allcat,trial_choice_allcat],{'mean','gname'});

grp = cellfun(@str2num,grp_name);

figure;
stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = 0;
subplot(1,2,1); hold on;
set(gca,'ColorOrder',stim_col);
plot(t,roi_trace_grp(grp(:,2) == 1,:)');
title('Orient left');
subplot(1,2,2); hold on;
set(gca,'ColorOrder',stim_col);
plot(t,roi_trace_grp(grp(:,2) == -1,:)');
title('Orient right');
linkaxes(get(gcf,'Children'),'xy')





move_col = [0.6,0,0.6;1,0.6,0];



%% (testing Peter model)


use_trials = true(size(trial_stim_allcat));

% Get behavioural data
D = struct;
D.stimulus = zeros(sum(use_trials),2);

L_trials = trial_choice_allcat == -1;
R_trials = trial_choice_allcat == 1;

D.stimulus(L_trials(use_trials),1) = trial_stim_allcat(L_trials & use_trials);
D.stimulus(R_trials(use_trials),2) = trial_stim_allcat(R_trials & use_trials);

D.response = ((trial_choice_allcat(use_trials)+1)/2)+1;
D.repeatNum = ones(sum(use_trials),1);

% Fit fluor model (across all timepoints)
cvfold = 5;
use_Vs = 1:100;

pV = nan(size(U_master,1),size(U_master,2),length(t),2);
for curr_t = 1:length(t)
    
    %                 D.offset_ZL = g.ZL(behavParameterFit, g.Zinput(g.data));
    D.neur = double(squeeze(fluor_allcat_deconv(use_trials,curr_t,use_Vs)));
    g_neur = GLM(D).setModel('AP_test_V').fitCV(cvfold);
    %                 pL = g_neur.p_hat(:,1);
    %                 pR = g_neur.p_hat(:,2);
    %                 likelihood = pL.*(g_neur.data.response==1) + pR.*(g_neur.data.response==2);
    %                 loglik_bpt_fluor = mean(log2(likelihood));
    
    pV_mean_cv = nanmean(g_neur.parameterFits(:,2:end),1);
    pV(:,:,curr_t) = svdFrameReconstruct(U_master(:,:,use_Vs),pV_mean_cv');
    AP_print_progress_fraction(curr_t,length(t));
end
















