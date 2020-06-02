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

animals = {'AP063','AP068'};

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    experiments = experiments([experiments.imaging]);
    
    im_stim_all = cell(length(experiments),1);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        AP_load_experiment;
        
        Udf_aligned = AP_align_widefield(Udf,animal,day,'day_only');
        fVdf_deconv = AP_deconv_wf(fVdf);
        
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
            
            peri_stim_v = permute(interp1(frame_t,fVdf_deconv',stim_surround_times),[3,2,1]);
            baseline_v = permute(nanmean(interp1(frame_t,fVdf_deconv',stim_baseline_surround_times),2),[3,2,1]);
            
            stim_v_mean = nanmean(peri_stim_v - baseline_v,3);
            
            im_stim(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf_aligned(:,:,use_vs), ...
                stim_v_mean(use_vs,:));
        end        
        
        use_stim = 3;
        im_stim_all{curr_day} = im_stim(:,:,:,use_stim);
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
    
    % Get average stim response across days
    stim_time = [0.05,0.15];
    use_time = surround_time >= stim_time(1) & surround_time <= stim_time(2);
    right_grating = nanmean(cell2mat(reshape(cellfun(@(x) ...
        nanmean(x(:,:,use_time),3),im_stim_all,'uni',false),1,1,[])),3);
   
    % Save average stim response
    save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment\right_grating';
    save_fn = [save_path filesep animal '_right_grating.mat'];
    save(save_fn,'right_grating');
    
end


%% Align corticostriatal across animals (using right grating response)

animal = 'AP068';

right_grating_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment\right_grating';
right_grating_fn = [right_grating_path filesep animal '_right_grating.mat'];
load(right_grating_fn);

% Make master from passive grating response
data_fn = 'trial_activity_AP_choiceWorldStimPassive_trained';
AP_load_concat_normalize_ctx_str;

stim_time = [0.05,0.15];
use_time = t >= stim_time(1) & t <= stim_time(2);
right_grating_master = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    squeeze(nanmean(nanmean(fluor_allcat_deconv(trial_stim_allcat == 1,use_time,:),2),1)));

% Align animal to master
AP_align_widefield(right_grating,animal,[],'new_animal',right_grating_master);


%% ~~~~~~~~~ BATCH ANALYSIS ~~~~~~~~~

%% Get passive responses before/after training

animals = {'AP063','AP068'};

im_stim_all = cell(size(animals));

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    experiments = experiments([experiments.imaging]);
        
    disp(animal);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        load_parts.imaging = true;
        AP_load_experiment;
        
        Udf_aligned = AP_align_widefield(Udf,animal,day);
        fVdf_deconv = AP_deconv_wf(fVdf);
        
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
            
            peri_stim_v = permute(interp1(frame_t,fVdf_deconv',stim_surround_times),[3,2,1]);
            baseline_v = permute(nanmean(interp1(frame_t,fVdf_deconv',stim_baseline_surround_times),2),[3,2,1]);
            
            stim_v_mean = nanmean(peri_stim_v - baseline_v,3);
            
            im_stim(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf_aligned(:,:,use_vs), ...
                stim_v_mean(use_vs,:));
        end        
        
        im_stim_all{curr_animal}{curr_day} = im_stim;       
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end    
    
end

use_stim = 3;
use_animal = 2;
curr_im = cell2mat(reshape(cellfun(@(x) x(:,:,:,use_stim), ...
    im_stim_all{use_animal},'uni',false),1,1,1,[]));

AP_image_scroll(curr_im);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));





%% Load corticostriatal trial data

% Task
data_fn = 'trial_activity_choiceworld_corticostriatal';
AP_load_concat_normalize_ctx_str;

% Passive
data_fn = 'trial_activity_AP_lcrGratingPassive_corticostriatal';
AP_load_concat_normalize_ctx_str;

% Passive naive
data_fn = 'trial_activity_AP_lcrGratingPassive_corticostriatal_naive';
AP_load_concat_normalize_ctx_str;


% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_recording;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));



% (testing: look at cortex>striatum kernel - garbage at the moment)
% Get average ctx>str kernel in pixels (flip time to be lag-to-MUA)
ctx_str_k_mean = fliplr(nanmean(cell2mat(permute(vertcat(ctx_str_k_all{:}),[2,3,4,5,1])),5));
ctx_str_k_mean_px = cell2mat(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:100), ...
    ctx_str_k_mean(:,:,x)),permute(1:n_depths,[1,3,4,2]),'uni',false));

ctx_str_kernel_frames_t = [-0.5,0.5];
ctx_str_kernel_frames = round(ctx_str_kernel_frames_t(1)*sample_rate): ...
    round(ctx_str_kernel_frames_t(2)*sample_rate);
ctx_str_kernel_t = ctx_str_kernel_frames./sample_rate;

AP_image_scroll(ctx_str_k_mean_px,ctx_str_kernel_t)
axis image
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));



% Average fluorescence
quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > 0.025,2);

unique_stim = unique(trial_stim_allcat);
stim_px = nan(size(U_master,1),size(U_master,2),length(t),length(unique_stim),length(use_split));
for curr_exp = 1:length(use_split)
    for curr_stim = 1:length(unique_stim)
        curr_v = permute(nanmean(fluor_allcat_deconv(split_idx == curr_exp & ...
            trial_stim_allcat == unique_stim(curr_stim),:,:),1),[3,2,1]);
        curr_px = svdFrameReconstruct(U_master(:,:,1:n_vs),curr_v);
        stim_px(:,:,:,curr_stim,curr_exp) = curr_px;
    end
end

trial_stim_allcat_exp = mat2cell(trial_stim_allcat,trials_recording,1);
fluor_allcat_deconv_exp = mat2cell(fluor_allcat_deconv,use_split,length(t),n_vs);
quiescent_trials_exp = mat2cell(quiescent_trials,trials_recording,1);

curr_v = cellfun(@(act,stim,quiesc) ...
    permute(nanmean(act(stim == 1 & quiesc,:,:),1),[3,2,1]), ...
    fluor_allcat_deconv_exp,trial_stim_allcat_exp,quiescent_trials_exp,'uni',false);

curr_v_mean = nanmean(cat(3,curr_v{:}),3);

curr_px = svdFrameReconstruct(U_master(:,:,1:n_vs),curr_v_mean);

AP_image_scroll(curr_px,t)
axis image;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));

















