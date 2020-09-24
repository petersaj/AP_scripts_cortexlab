%% Test analysis for recording across learning

%% Plot psychometrics and reaction times

animal = 'AP041';
protocol = 'vanillaChoiceworld';
experiments = AP_find_experiments(animal,protocol);

% experiments = experiments([experiments.imaging] & [experiments.ephys]);
experiments = experiments([experiments.imaging] & ~[experiments.ephys]);

init_array = cell(size(experiments));
bhv = struct('move_t',init_array,'stim_contrastsides',init_array,'trial_choice',init_array);

for curr_day = 1:length(experiments)
    
    load_parts.cam = false;
    load_parts.imaging = false;
    load_parts.ephys = false;
    
    day = experiments(curr_day).day;
    experiment = experiments(curr_day).experiment(end);
    
    AP_load_experiment
    framerate = 35;
    
    % Get event-aligned activity
    raster_window = [-0.5,2];
    upsample_factor = 1;
    raster_sample_rate = 1/(framerate*upsample_factor);
    t = raster_window(1):raster_sample_rate:raster_window(2);
    
    % Get align times
    use_align = stimOn_times;
    use_align(isnan(use_align)) = 0;
    
    t_peri_event = bsxfun(@plus,use_align,t);
    t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
    
    %%% Trial-align wheel velocity
    event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
        wheel_velocity,t_peri_event);
    
    % Get reaction time for building regressors
    [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.03,[],2);
    move_idx(~move_trial) = NaN;
    move_t = nan(size(move_idx));
    move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
    
    % Get stim
    unique_stim = unique(contrasts(contrasts > 0).*sides');
    stim_contrastsides = ...
        signals_events.trialSideValues(1:length(stimOn_times))'.* ...
        signals_events.trialContrastValues(1:length(stimOn_times))';
    
    % Package
    bhv(curr_day).day = day;
    bhv(curr_day).move_t = move_t;
    bhv(curr_day).stim_contrastsides = stim_contrastsides;
    bhv(curr_day).trial_choice = trial_choice;
    
    AP_print_progress_fraction(curr_day,length(experiments));
    
end

% Plot all psychometrics and reaction times overlaid
figure; hold on
psychometric_ax = subplot(1,2,1); hold on;
line(psychometric_ax,[-1,1],[0.5,0.5],'linestyle','--','color','k');
line(psychometric_ax,[0,0],[0,1],'linestyle','--','color','k');

rxn_ax = subplot(1,2,2); hold on;
line(rxn_ax,[-1,1],[0.5,0.5],'linestyle','--','color','k');
line(rxn_ax,[0,0],[0,1],'linestyle','--','color','k');

% Plot psychometric and number of trials over days
unique_stim_contrastsides = unique(vertcat(bhv.stim_contrastsides));
frac_left_all = nan(length(bhv),11);
rxn_all = nan(length(bhv),11);
for i = 1:length(bhv)
    n_trials = length(bhv(i).move_t);
    frac_left = grpstats(bhv(i).trial_choice(1:n_trials) == -1,bhv(i).stim_contrastsides(1:n_trials));
    plot(psychometric_ax,unique(bhv(i).stim_contrastsides),frac_left,'k','linewidth',2);
    
    rxn = grpstats(bhv(i).move_t(1:n_trials),bhv(i).stim_contrastsides(1:n_trials),{'nanmedian'});
    plot(rxn_ax,unique(bhv(i).stim_contrastsides),rxn,'k','linewidth',2);
    
    curr_stim_contrastsides = unique(bhv(i).stim_contrastsides(1:n_trials));
    [~,curr_stim_contrastsides_idx] = ismember(curr_stim_contrastsides,unique_stim_contrastsides);
    frac_left_all(i,curr_stim_contrastsides_idx) = frac_left;
    rxn_all(i,curr_stim_contrastsides_idx) = rxn;
end

n_trials_all = cellfun(@length,{bhv.move_t});

figure('name',animal); 

subplot(3,1,1);
plot(n_trials_all,'k','linewidth',2);
ylabel('Number of trials');
c = colorbar;
set(c,'visible','off');
set(gca,'XTick',1:length(experiments))
set(gca,'XTickLabel',{experiments.day})
set(gca,'XTickLabelRotation',30)

subplot(3,1,2);
imagesc([],unique_stim_contrastsides,frac_left_all','AlphaData',~isnan(frac_left_all'));
colormap(brewermap([],'*RdBu'));
set(gca,'color',[0.5,0.5,0.5]);
c = colorbar; ylabel(c,'Fraction orient left');
xlabel('Session');
ylabel('Stim contrast*side');

subplot(3,1,3);
imagesc([],unique_stim_contrastsides,rxn_all','AlphaData',~isnan(rxn_all'));
colormap(brewermap([],'*RdBu'));
set(gca,'color',[0.5,0.5,0.5]);
caxis([0.2,0.8]);
c = colorbar; ylabel(c,'Time to first move');
xlabel('Session');
ylabel('Stim contrast*side');


%% Passive stim (pixels)

animal = 'AP077';

protocol = 'AP_lcrGratingPassive';
experiments = AP_find_experiments(animal,protocol);

% Use only days with imaging
experiments = experiments([experiments.imaging]);

im_stim_all = cell(size(experiments));
for curr_day = 1:length(experiments)
    
    % Load data
    day = experiments(curr_day).day;
    experiment = experiments(curr_day).experiment(end);
    load_parts.imaging = true;
    AP_load_experiment;
    
%     fVdf_deconv = AP_deconv_wf(fVdf);
    %%%% TEMP
    fVdf_deconv = fVdf;
    %%%%
    
    % Get wheel movements during stim
    wheel_window = [0,0.5];
    wheel_window_t = wheel_window(1):1/framerate:wheel_window(2);
    wheel_window_t_peri_event = bsxfun(@plus,stimOn_times,wheel_window_t);
    event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
        wheel_velocity,wheel_window_t_peri_event);
    wheel_thresh = 0.025;
    quiescent_trials = ~any(abs(event_aligned_wheel) > wheel_thresh,2);
    
    % Set options
    surround_window = [-0.5,1];
    baseline_window = [-0.1,0];
    
    surround_samplerate = 1/(framerate*1);
    surround_time = surround_window(1):surround_samplerate:surround_window(2);
    baseline_surround_time = baseline_window(1):surround_samplerate:baseline_window(2);
    
    % Average (time course) responses
    use_vs = 1:size(U,3);
    
    conditions = unique(stimIDs);
    im_stim = nan(size(U,1),size(U,2),length(surround_time),length(conditions));
    for curr_condition_idx = 1:length(conditions)
        curr_condition = conditions(curr_condition_idx);
        
        use_stims = find(stimIDs == curr_condition & quiescent_trials);
        use_stimOn_times = stimOn_times(use_stims);
        
        stim_surround_times = bsxfun(@plus, use_stimOn_times(:), surround_time);
        stim_baseline_surround_times = bsxfun(@plus, use_stimOn_times(:), baseline_surround_time);
        
        peri_stim_v = permute(interp1(frame_t,fVdf_deconv',stim_surround_times),[3,2,1]);
        baseline_v = permute(nanmean(interp1(frame_t,fVdf_deconv',stim_baseline_surround_times),2),[3,2,1]);
        
        stim_v_mean = nanmean(peri_stim_v - baseline_v,3);
        
        im_stim(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf(:,:,use_vs), ...
            stim_v_mean(use_vs,:));
    end
    
    im_stim_aligned = AP_align_widefield(im_stim,animal,day);
    
    im_stim_all{curr_day} = im_stim_aligned;

    AP_print_progress_fraction(curr_day,length(experiments));
end

im_stim_cat = cell2mat(permute(cellfun(@(im) ...
    reshape(permute(im,[1,2,4,3]),size(im,1),[],size(im,3)), ...
    im_stim_all,'uni',false),[2,3,4,1]));

AP_image_scroll(im_stim_cat);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));

% t_stim = 21:23;
t_stim = 26:29;
im_stim_avg = squeeze(max(im_stim_cat(:,:,t_stim,:),[],3));
AP_image_scroll(im_stim_avg);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));


%% Passive stim (U master)

animal = 'AP077';

protocol = 'AP_lcrGratingPassive';
experiments = AP_find_experiments(animal,protocol);

% Remove ephys days - those had muscimol
experiments = experiments([experiments.imaging] & ~[experiments.ephys]);

stim_v_all = cell(length(experiments),3);
for curr_day = 1:length(experiments)
    
    % Load data
    day = experiments(curr_day).day;
    experiment = experiments(curr_day).experiment(end);
    load_parts.imaging = true;
    AP_load_experiment;
    
    % Convert U to master U
    load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
    Udf_aligned = AP_align_widefield(Udf,animal,day);
    fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);

    % Deconvolve
    fVdf_recast_deconv = AP_deconv_wf(fVdf_recast);
    
    % Get wheel movements during stim
    wheel_window = [0,0.5];
    wheel_window_t = wheel_window(1):1/framerate:wheel_window(2);
    wheel_window_t_peri_event = bsxfun(@plus,stimOn_times,wheel_window_t);
    event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
        wheel_velocity,wheel_window_t_peri_event);
    wheel_thresh = 0.025;
    quiescent_trials = ~any(abs(event_aligned_wheel) > wheel_thresh,2);
    
    % Set options
    surround_window = [-0.5,1];
    baseline_window = [-0.1,0];
    
    surround_samplerate = 1/(framerate*1);
    surround_time = surround_window(1):surround_samplerate:surround_window(2);
    baseline_surround_time = baseline_window(1):surround_samplerate:baseline_window(2);
    
    % Average (time course) responses
    use_vs = 1:size(U,3);
    
    conditions = unique(stimIDs);
    stim_v = nan(size(U,3),length(surround_time),length(conditions));
    for curr_condition_idx = 1:length(conditions)
        curr_condition = conditions(curr_condition_idx);
        
        use_stims = find(stimIDs == curr_condition & quiescent_trials);
        use_stimOn_times = stimOn_times(use_stims);
        
        stim_surround_times = bsxfun(@plus, use_stimOn_times(:), surround_time);
        stim_baseline_surround_times = bsxfun(@plus, use_stimOn_times(:), baseline_surround_time);
        
        peri_stim_v = permute(interp1(frame_t,fVdf_recast_deconv',stim_surround_times),[3,2,1]);
        baseline_v = permute(nanmean(interp1(frame_t,fVdf_recast_deconv',stim_baseline_surround_times),2),[3,2,1]);
        
        stim_v = peri_stim_v - baseline_v;
       
        stim_v_all{curr_day,curr_condition_idx} = stim_v;
        
    end
        
    AP_print_progress_fraction(curr_day,length(experiments));
end

% Plot average from one stim
curr_px = svdFrameReconstruct(U_master,nanmean(cat(3,stim_v_all{:,1}),3));
AP_image_scroll(curr_px);
axis image off;
colormap(brewermap([],'*RdBu'));
caxis([-max(abs(caxis)),max(abs(caxis))]);

% Plot mean or median within each day
a = cell2mat(permute(cellfun(@(x) svdFrameReconstruct(U_master,nanmean(x,3)),stim_v_all(:,1),'uni',false),[2,3,4,1]));
AP_image_scroll(a);
axis image off;
colormap(brewermap([],'*RdBu'));
caxis([-max(abs(caxis)),max(abs(caxis))]);

b = squeeze(max(a(:,:,20:23,:),[],3));
AP_image_scroll(b);
axis image off;
colormap(brewermap([],'*RdBu'));
caxis([-max(abs(caxis)),max(abs(caxis))]);

% Plot smoothed average
a = cat(3,stim_v_all{:,1});
b = squeeze(nanmean(a(:,20:23,:),2));
c = convn(b,ones(1,150)./150,'same');
curr_px = svdFrameReconstruct(U_master,c);
AP_image_scroll(curr_px);
axis image off;
colormap(brewermap([],'*RdBu'));
caxis([-max(abs(caxis)),max(abs(caxis))]);




%% Get task -> widefield kernel

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

% Prepare fluorescence

% % (if aligned)
% % Convert U to master U
% load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
% Udf_aligned = AP_align_widefield(Udf,animal,day);
% fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);

% (if not aligned)
U_master = Udf;
Udf_aligned = Udf;
fVdf_recast = fVdf;

% Set components to keep
use_components = 1:200;

% Get event-aligned activity
raster_window = [-0.5,2];
upsample_factor = 1;
raster_sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):raster_sample_rate:raster_window(2);

% Get align times
use_align = stimOn_times;
use_align(isnan(use_align)) = 0;

t_peri_event = bsxfun(@plus,use_align,t);
t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];

%%% Trial-align cortex
event_aligned_V = ...
    interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);

% Get time points to bin
sample_rate = framerate*regression_params.upsample_factor;
time_bins = frame_t(find(frame_t > ...
    regression_params.skip_seconds,1)):1/sample_rate: ...
    frame_t(find(frame_t-frame_t(end) < ...
    -regression_params.skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Deconvolve fluoresence
fVdf_deconv = AP_deconv_wf(fVdf);
fVdf_deconv(isnan(fVdf_deconv)) = 0;
fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';

%%% Trial-align wheel velocity
event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
    wheel_velocity,t_peri_event);

%%% Regress task to cortex/striatum/cortex-predicted striatum

% Get reaction time for building regressors
% (wheel thresh is the max in the quiescent period)
wheel_thresh = max(reshape(abs(event_aligned_wheel(:,t < -0.2)),[],1));
[move_trial,move_idx] = max(abs(event_aligned_wheel) > wheel_thresh,[],2);
move_idx(~move_trial) = NaN;
move_t = nan(size(move_idx));
move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';

% Build regressors (only a subset of these are used)

% Stim regressors
unique_stim = unique(contrasts(contrasts > 0).*sides');
stim_contrastsides = ...
    signals_events.trialSideValues(1:length(stimOn_times))'.* ...
    signals_events.trialContrastValues(1:length(stimOn_times))';

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

stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
for curr_contrast = 1:length(unique_contrasts)
    
    % (find the last photodiode flip before the reward)
    curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ...
        abs(stim_contrastsides) == unique_contrasts(curr_contrast));
    
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
    find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == -1));
move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
    find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == 1));

move_onset_regressors = zeros(2,length(time_bin_centers));
move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);

% Move onset x stim regressors (one for each contrast/side)
move_onset_stim_time_absolute = arrayfun(@(curr_stim) ...
    arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
    find(~isnan(move_idx) & stim_contrastsides == unique_stim(curr_stim))), ...
    1:length(unique_stim),'uni',false);

move_onset_stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
for curr_stim = 1:length(unique_stim)
    move_onset_stim_regressors(curr_stim,:) = ...
        histcounts(move_onset_stim_time_absolute{curr_stim},time_bins);
end

% Move ongoing regressors (L/R choice for duration of movement)
wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);

move_stopped_t = 0.5;
move_stopped_samples = round(sample_rate*move_stopped_t);
wheel_thresh = 0.025;
wheel_moving_conv = convn((abs(wheel_velocity_interp) > wheel_thresh), ...
    ones(1,move_stopped_samples),'full') > 0;
wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);

move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
    find(time_bin_centers > x): ...
    find(time_bin_centers > x & ~wheel_moving,1), ...
    move_time_L_absolute','uni',false));
move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
    find(time_bin_centers > x): ...
    find(time_bin_centers > x & ~wheel_moving,1), ...
    move_time_R_absolute','uni',false));

move_ongoing_regressors = zeros(2,length(time_bin_centers));
move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
move_ongoing_regressors(2,move_ongoing_R_samples) = 1;

% Go cue regressors - separate for early/late move
% (using signals timing - not precise but looks good)
% (for go cue only on late move trials)
%         go_cue_regressors = histcounts( ...
%             signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
% (for go cue with early/late move trials)
go_cue_regressors = zeros(1,length(time_bin_centers));
go_cue_regressors(1,:) = histcounts( ...
    signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
go_cue_regressors(2,:) = histcounts( ...
    signals_events.interactiveOnTimes(move_t > 0.5),time_bins);

% Outcome regressors
% (using signals timing - not precise but looks good)
% (regressors for hit only)
%         outcome_regressors = histcounts(reward_t_timeline,time_bins);
% (regressors for both hit and miss)
outcome_regressors = zeros(2,length(time_bin_centers));
outcome_regressors(1,:) = histcounts( ...
    reward_t_timeline,time_bins);
outcome_regressors(2,:) = histcounts( ...
    signals_events.responseTimes(trial_outcome == -1),time_bins);

% Concatenate selected regressors, set parameters

task_regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
task_regressor_labels = {'Stim','Move onset','Go cue','Outcome'};

task_t_shifts = { ...
    [0,0.5]; ... % stim
    [-0.5,1]; ... % move
    [0,0.5]; ... % go cue
    [0,0.5]}; % outcome

task_regressor_sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
    round(x(2)*(sample_rate)),task_t_shifts,'uni',false);

lambda = 0;
zs = [false,false];
cvfold = 5;
use_constant = false;
return_constant = false;

% Regression task -> (master U, deconvolved) fluor
event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,Udf_aligned);

baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,:),[],size(event_aligned_V_deconv,3)))';
activity = single(fVdf_deconv_resample_recast(use_components,:))-baseline;

[fluor_taskpred_k,fluor_taskpred_long,fluor_taskpred_expl_var,fluor_taskpred_reduced_long] = ...
    AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
    lambda,zs,cvfold,return_constant,use_constant);


% Get regressor pixels
regressor_px = cellfun(@(v) cell2mat(arrayfun(@(subregressor) ...
    svdFrameReconstruct(U_master(:,:,use_components),permute(v(subregressor,:,:),[3,2,1])), ...
    permute(1:size(v,1),[1,3,4,2]),'uni',false)),fluor_taskpred_k,'uni',false);


AP_image_scroll(regressor_px{4});
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));


%% Previous passive data: get stimulus responses across days
% Q: does stimulus response change in passive over days?

% Load passive data
data_fn = 'trial_activity_AP_choiceWorldStimPassive_naive';
exclude_data = false;
AP_load_concat_normalize_ctx_str;

% Split data into experiment
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));

% Reshape data into day x animal
animal_n_days = cellfun(@length,fluor_all);
day_animal_idx = false(max(animal_n_days),length(fluor_all));
for curr_animal = 1:length(fluor_all)
   day_animal_idx(1:animal_n_days(curr_animal),curr_animal) = true; 
end

fluor_allcat_deconv_exp = mat2cell(fluor_allcat_deconv,trials_recording,length(t),n_vs);
fluor_allcat_deconv_exp_reshape = cell(max(animal_n_days),length(fluor_all));
fluor_allcat_deconv_exp_reshape(day_animal_idx) = fluor_allcat_deconv_exp;

stim_exp = mat2cell(D_allcat.stimulus,trials_recording,1);
stim_exp_reshape = cell(max(animal_n_days),length(fluor_all));
stim_exp_reshape(day_animal_idx) = stim_exp;

wheel_thresh = 0.025;
quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > wheel_thresh,2);
quiescent_trials_exp = mat2cell(quiescent_trials,trials_recording,1);
quiescent_trials_exp_reshape = cell(max(animal_n_days),length(fluor_all));
quiescent_trials_exp_reshape(day_animal_idx) = quiescent_trials_exp;

% Get averages of single stimuli during quiescent trials
fluor_stim = cellfun(@(fluor,stim,quiescent) ...
    permute(nanmean(fluor(stim == 1 & quiescent,:,:),1),[3,2,1]), ...
    fluor_allcat_deconv_exp_reshape,stim_exp_reshape,quiescent_trials_exp_reshape,'uni',false);

fluor_stim_day = cell(max(animal_n_days),1);
for curr_day = 1:max(animal_n_days)
    use_animals = cellfun(@(x) ~isempty(x),fluor_stim(curr_day,:));
    use_animals(2) = false;
    curr_fluor = nanmean(cat(3,fluor_stim{curr_day,use_animals}),3);
    fluor_stim_day{curr_day} = svdFrameReconstruct(U_master(:,:,1:n_vs),curr_fluor);
end

fluor_stim_day = cat(4,fluor_stim_day{:});

AP_image_scroll(fluor_stim_day,t);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k');































