%% ~~~~~~ Test analysis for muscimol during recording ~~~~~~
% (AP040 and AP041 at the moment)

%% Plot psychometrics and reaction times (only one choiceworld)

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

figure; 

subplot(2,1,1);
plot(n_trials_all,'k','linewidth',2);
xlabel('Session');
ylabel('Number of trials');

subplot(2,1,2);
imagesc(frac_left_all');
colormap(brewermap([],'*RdBu'));

%% Plot psychometrics and reaction times (pre/post muscimol choiceworld)

animal = 'AP047';
protocol = 'vanillaChoiceworld';
flexible_name = true;
experiments = AP_find_experiments(animal,protocol,flexible_name);

% experiments = experiments([experiments.imaging] & [experiments.ephys]);

init_array = cell(size(experiments));

curr_day = length(experiments);
day = experiments(curr_day).day;

if length(experiments(curr_day).experiment) ~= 2
    error('~= 2 experiments')
end

bhv = struct('frac_orient_right',init_array,'rxn_time',init_array);
for curr_exp = 1:2
    
    experiment = experiments(curr_day).experiment(curr_exp);
    
    [block_filename, block_exists] = AP_cortexlab_filename(animal,day,experiment,'block');
    
    % Load the block file
    load(block_filename)
    
    % Get protocol name
    [~,curr_protocol] = fileparts(block.expDef);
    
    % Time of session (in minutes)
    session_duration = block.duration/60;
    
    % Trial counts
    n_trials = length(block.paramsValues);
    total_water = sum(block.outputs.rewardValues);
    
    % Estimate reaction time
    % (evenly resample velocity - not even natively);
    wheel_resample_rate = 1000;
    wheel_t_resample = block.inputs.wheelTimes(1):1/wheel_resample_rate:block.inputs.wheelTimes(end);
    wheel_values_resample = interp1(block.inputs.wheelTimes,block.inputs.wheelValues,wheel_t_resample);
    
    wheel_smooth_t = 0.05; % seconds
    wheel_smooth_samples = round(wheel_smooth_t*wheel_resample_rate);
    wheel_velocity = interp1(conv(wheel_t_resample,[1,1]/2,'valid'), ...
        diff(smooth(wheel_values_resample,wheel_smooth_samples)),wheel_t_resample)';
    
    wheel_thresh = 0.025;
    wheel_starts = wheel_t_resample(abs(wheel_velocity(1:end-1)) < wheel_thresh & ...
        abs(wheel_velocity(2:end)) > wheel_thresh);
    
    response_trials = 1:length(block.events.responseValues);
    trial_wheel_starts = arrayfun(@(x) ...
        wheel_starts(find(wheel_starts > block.events.stimOnTimes(x),1)), ...
        response_trials);
    trial_move_t = trial_wheel_starts - block.events.stimOnTimes(response_trials);
    
    % Wheel movements/biases
    left_wheel_velocity = abs(wheel_velocity.*(wheel_velocity < 0));
    right_wheel_velocity = abs(wheel_velocity.*(wheel_velocity > 0));
    wheel_bias = (nansum(right_wheel_velocity)-nansum(left_wheel_velocity))/ ...
        (nansum(right_wheel_velocity)+nansum(left_wheel_velocity));
    
    % Get reaction times for each stimulus
    trial_stim = block.events.trialContrastValues(response_trials).*block.events.trialSideValues(response_trials);
    stim_list = unique(reshape(unique(block.events.contrastsValues).*[-1;1],[],1));
    [~,trial_stim_idx] = ismember(trial_stim,stim_list);
    stim_rxn_time = accumarray(trial_stim_idx(response_trials)',trial_move_t',[11,1],@nanmedian,nan);
    
    % Performance (note that this excludes repeat on incorrect trials)
    performance = block.events.sessionPerformanceValues(:,end-10:end);
    
    % Get whether all contrasts were used
    use_all_contrasts = all(block.events.useContrastsValues(end-10:end));
    
    
    % Package final 
    bhv(curr_exp).frac_orient_right = performance(3,:)./performance(2,:);
    bhv(curr_exp).rxn_time = stim_rxn_time;
    
end

figure; 
subplot(1,2,1); hold on;
plot(performance(1,:),bhv(1).frac_orient_right,'k','linewidth',2);
plot(performance(1,:),bhv(2).frac_orient_right,'r','linewidth',2);
line(xlim,[0.5,0.5],'color','k','linestyle','--');
line([0,0],[0,1],'color','k','linestyle','--');

subplot(1,2,2); hold on;
plot(performance(1,:),bhv(1).rxn_time,'k','linewidth',2);
plot(performance(1,:),bhv(2).rxn_time,'r','linewidth',2);
line([-1,1],[0.5,0.5],'color','k','linestyle','--');

%% Plot psychometrics and reaction times (pre/post muscimol choiceworld, BATCH)

animals = {'AP045','AP054','AP055','AP053','AP047','AP048'};

bhv = struct('frac_orient_right',cell(length(animals),1), ...
    'rxn_time',cell(length(animals),1));

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);   
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        
        if length(experiments(curr_day).experiment) ~= 2
            error('~= 2 experiments')
        end
        
        for curr_exp = 1:2
            
            experiment = experiments(curr_day).experiment(curr_exp);
            
            [block_filename, block_exists] = AP_cortexlab_filename(animal,day,experiment,'block');
            
            % Load the block file
            load(block_filename)
            
            % Get protocol name
            [~,curr_protocol] = fileparts(block.expDef);
            
            % Time of session (in minutes)
            session_duration = block.duration/60;
            
            % Trial counts
            n_trials = length(block.paramsValues);
            total_water = sum(block.outputs.rewardValues);
            
            % Estimate reaction time
            % (evenly resample velocity - not even natively);
            wheel_resample_rate = 1000;
            wheel_t_resample = block.inputs.wheelTimes(1):1/wheel_resample_rate:block.inputs.wheelTimes(end);
            wheel_values_resample = interp1(block.inputs.wheelTimes,block.inputs.wheelValues,wheel_t_resample);
            
            wheel_smooth_t = 0.05; % seconds
            wheel_smooth_samples = round(wheel_smooth_t*wheel_resample_rate);
            wheel_velocity = interp1(conv(wheel_t_resample,[1,1]/2,'valid'), ...
                diff(smooth(wheel_values_resample,wheel_smooth_samples)),wheel_t_resample)';
            
            wheel_thresh = 0.025;
            wheel_starts = wheel_t_resample(abs(wheel_velocity(1:end-1)) < wheel_thresh & ...
                abs(wheel_velocity(2:end)) > wheel_thresh);
            
            response_trials = 1:length(block.events.responseValues);
            trial_wheel_starts = arrayfun(@(x) ...
                wheel_starts(find(wheel_starts > block.events.stimOnTimes(x),1)), ...
                response_trials);
            trial_move_t = trial_wheel_starts - block.events.stimOnTimes(response_trials);
            
            % Wheel movements/biases
            left_wheel_velocity = abs(wheel_velocity.*(wheel_velocity < 0));
            right_wheel_velocity = abs(wheel_velocity.*(wheel_velocity > 0));
            wheel_bias = (nansum(right_wheel_velocity)-nansum(left_wheel_velocity))/ ...
                (nansum(right_wheel_velocity)+nansum(left_wheel_velocity));
            
            % Get reaction times for each stimulus
            trial_stim = block.events.trialContrastValues(response_trials).*block.events.trialSideValues(response_trials);
            stim_list = unique(reshape(unique(block.events.contrastsValues).*[-1;1],[],1));
            [~,trial_stim_idx] = ismember(trial_stim,stim_list);
            stim_rxn_time = accumarray(trial_stim_idx(response_trials)',trial_move_t',[11,1],@nanmedian,nan);
            
            % Performance (note that this excludes repeat on incorrect trials)
            performance = block.events.sessionPerformanceValues(:,end-10:end);
            
            % Get whether all contrasts were used
            use_all_contrasts = all(block.events.useContrastsValues(end-10:end));
            
            
            % Package final
            bhv(curr_animal).frac_orient_right{curr_day}(curr_exp,:) = performance(3,:)./performance(2,:);
            bhv(curr_animal).rxn_time{curr_day}(curr_exp,:) = stim_rxn_time;
            
        end  
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
end

frac_orient_right_cat = cell2mat(permute(cat(2,bhv.frac_orient_right),[1,3,2]));
rxn_time_cat = cell2mat(permute(cat(2,bhv.rxn_time),[1,3,2]));

figure; 
subplot(2,2,1); hold on;
plot(performance(1,:),squeeze(frac_orient_right_cat(1,:,:)),'color',[0.5,0.5,0.5],'linewidth',1);
plot(performance(1,:),squeeze(frac_orient_right_cat(2,:,:)),'color',[1,0.5,0.5],'linewidth',1);
plot(performance(1,:),nanmean(frac_orient_right_cat(1,:,:),3),'color','k','linewidth',3);
plot(performance(1,:),nanmean(frac_orient_right_cat(2,:,:),3),'color','r','linewidth',3);
line(xlim,[0.5,0.5],'color','k','linestyle','--');
line([0,0],[0,1],'color','k','linestyle','--');
xlabel('Side*Contrast');
ylabel('Fraction orient right')

subplot(2,2,2); hold on;
plot(performance(1,:),squeeze(rxn_time_cat(1,:,:)),'color',[0.5,0.5,0.5],'linewidth',1);
plot(performance(1,:),squeeze(rxn_time_cat(2,:,:)),'color',[1,0.5,0.5],'linewidth',1);
plot(performance(1,:),nanmean(rxn_time_cat(1,:,:),3),'color','k','linewidth',3);
plot(performance(1,:),nanmean(rxn_time_cat(2,:,:),3),'color','r','linewidth',3);
line([-1,1],[0.5,0.5],'color','k','linestyle','--');
xlabel('Side*Contrast');
ylabel('Reaction time');

subplot(2,2,3); hold on;
plot(performance(1,:),squeeze(diff(frac_orient_right_cat,[],1)),'color',[0.5,0.5,0.5],'linewidth',1);
plot(performance(1,:),nanmean(diff(frac_orient_right_cat,[],1),3),'color','k','linewidth',3);
line(xlim,[0,0],'color','k','linestyle','--');
line([0,0],[-1,1],'color','k','linestyle','--');
xlabel('Side*Contrast');
ylabel('\Delta Fraction orient right')

subplot(2,2,4); hold on;
plot(performance(1,:),squeeze(diff(rxn_time_cat,[],1)),'color',[0.5,0.5,0.5],'linewidth',1);
plot(performance(1,:),nanmean(diff(rxn_time_cat,[],1),3),'color','k','linewidth',3);
line([-1,1],[0,0],'color','k','linestyle','--');
xlabel('Side*Contrast');
ylabel('\Delta Reaction time');

%% Get VFS pre/post musicmol

animals = {'AP045','AP054','AP055','AP053','AP047','AP048'};

vfs = cell(length(animals),1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_sparseNoise';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);   
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        
        if length(experiments(curr_day).experiment) ~= 2
            warning([animal ' ' day '~= 2 experiments, using first/last']);
            experiments(curr_day).experiment = experiments(curr_day).experiment([1,end]);
        end
        
        for curr_exp = 1:2
            
            experiment = experiments(curr_day).experiment(curr_exp);
            
            AP_load_experiment;
            lilrig_retinotopy;
            vfs_aligned = AP_align_widefield(vfs_median,animal,day);
            
            % Package final
            vfs{curr_animal}{curr_day,curr_exp} = vfs_aligned;
            
        end         
        AP_print_progress_fraction(curr_day,length(experiments));
    end
end

% Plot pre/post for each animal and day
for curr_animal = 1:length(vfs)
    figure('Name',animals{curr_animal});
    for curr_day = 1:size(vfs{curr_animal},1)
        subplot(size(vfs{curr_animal},1),2,(curr_day-1)*2+1);
        imagesc(vfs{curr_animal}{curr_day,1});
        axis image off;
        colormap(brewermap([],'*RdBu'));
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        
        subplot(size(vfs{curr_animal},1),2,(curr_day-1)*2+2);
        imagesc(vfs{curr_animal}{curr_day,2});
        axis image off;
        colormap(brewermap([],'*RdBu'));
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    end
end

vfs_cat = vertcat(vfs{:});

% Plot pre/post mean
vfs_pre_mean = nanmean(cat(3,vfs_cat{:,1}),3);
vfs_post_mean = nanmean(cat(3,vfs_cat{:,2}),3);

figure; 
subplot(1,2,1);
imagesc(vfs_pre_mean);
axis image off
colormap(brewermap([],'*RdBu'));
caxis([-1,1])
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('Pre-muscimol');

subplot(1,2,2);
imagesc(vfs_post_mean);
axis image off
colormap(brewermap([],'*RdBu'));
caxis([-1,1])
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('Post-muscimol');

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\tests\muscimol_test';
save_fn = 'muscimol_vfs.mat';
save([save_path filesep save_fn],'vfs');
disp(['Saved ' [save_path filesep save_fn]]);


%% Quantify/plot VFS pre/post musicmol

muscimiol_vfs_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\tests\muscimol_test\muscimol_vfs.mat';
load(muscimiol_vfs_fn);

vfs_cat = cell2mat(permute(vertcat(vfs{:}),[3,2,1]));
AP_image_scroll(vfs_cat)
axis image
colormap(brewermap([],'*RdBu'));
caxis([-1,1]);


vfs_diff = diff(abs(cell2mat(permute(vertcat(vfs{:}),[3,4,1,2]))),[],4);





%% Plot pre/post of cortical muscimol recordings

spike_rate_exp = nan(max(spike_templates),4);
n_exp = length(sync(2).timestamps)/2;

for curr_exp = 1:n_exp
    
    curr_exp_start = sync(2).timestamps(curr_exp*2 - 1);
    curr_exp_stop = sync(2).timestamps(curr_exp*2);
    
    curr_use_spikes = spike_times >= curr_exp_start & ...
        spike_times <= curr_exp_stop;
    
    spike_rate_exp(:,curr_exp) = ...
        accumarray(spike_templates(curr_use_spikes),1,[max(spike_templates),1])./ ...
        (curr_exp_stop - curr_exp_start);
    
end

figure; 
p = gobjects(n_exp,1);
for i = 1:n_exp
   p(i) = subplot(1,n_exp,i);
   plot(log10(spike_rate_exp(:,i)),template_depths,'.k');
   set(gca,'YDir','reverse');
   xlabel('Log_{10} norm spikes/s')
end
linkaxes(p,'x');

use_baseline = 1:2;
use_muscimol = 3:4;

baseline_spikes = nanmean(spike_rate_exp(:,use_baseline),2);
muscimol_spikes = nanmean(spike_rate_exp(:,use_muscimol),2);

figure; hold on; colormap(jet);
scatter(baseline_spikes,muscimol_spikes,50,template_depths,'filled');
axis square; ylim(xlim);
xlabel('Baseline spikes/s');
ylabel('Muscimol spikes/s');
line(xlim,xlim,'color','k');



%% ~~~~~~~~ After preprocessing/saving trial activity ~~~~~~~~

%% Get cortex/striatum activity pre/post muscimol

data_fns = { ...
    'trial_activity_AP_lcrGratingPassive_pre_muscimol', ...
    'trial_activity_AP_lcrGratingPassive_post_muscimol'};

mua_stim_act_exp = [];

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Split data by experiment
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    
    trial_stim_allcat_exp = mat2cell(trial_stim_allcat,trials_recording,1);
    mua_allcat_exp = mat2cell(mua_allcat, ...
        trials_recording,size(mua_allcat,2),size(mua_allcat,3));
    mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat, ...
        trials_recording,size(mua_ctxpred_allcat,2),size(mua_ctxpred_allcat,3));
    fluor_roi_deconv_exp = mat2cell(fluor_roi_deconv, ...
        trials_recording,size(fluor_roi_deconv,2),size(fluor_roi_deconv,3));
    
    % Get trials with movement during stim to exclude
    wheel_thresh = 0.025;
    quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > wheel_thresh,2);
    quiescent_trials_exp = mat2cell(quiescent_trials,trials_recording,1);
    
    % Get stimulus activity for each recording
    use_stim_time = t >= 0.05 & t <= 0.15;
    use_stim = 3;
    
    mua_stim_act_exp(:,:,curr_data) = cell2mat(cellfun(@(act,stim,quiescent) ...
        squeeze(nanmean(nanmean(act(stim == use_stim & quiescent,use_stim_time,:,:),2),1)), ...
        mua_allcat_exp,trial_stim_allcat_exp,quiescent_trials_exp,'uni',false)');
    
    mua_ctxpred_stim_act_exp(:,:,curr_data) = cell2mat(cellfun(@(act,stim,quiescent) ...
        squeeze(nanmean(nanmean(act(stim == use_stim & quiescent,use_stim_time,:,:),2),1)), ...
        mua_ctxpred_allcat_exp,trial_stim_allcat_exp,quiescent_trials_exp,'uni',false)');
    
    fluor_roi_stim_act_exp(:,:,curr_data) = cell2mat(cellfun(@(act,stim,quiescent) ...
        squeeze(nanmean(nanmean(act(stim == use_stim & quiescent,use_stim_time,:,:),2),1)), ...
        fluor_roi_deconv_exp,trial_stim_allcat_exp,quiescent_trials_exp,'uni',false)');
    
    clearvars -except data_fns curr_data mua_stim_act_exp mua_ctxpred_stim_act_exp fluor_roi_stim_act_exp
    
    AP_print_progress_fraction(curr_data,length(data_fns));
    
end

% Load WF ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);

% Plot (pre-post)/(pre+post) for all areas
mua_change_idx = (mua_stim_act_exp(:,:,1)-mua_stim_act_exp(:,:,2))./(mua_stim_act_exp(:,:,1)+mua_stim_act_exp(:,:,2));
fluor_change_idx = (fluor_roi_stim_act_exp(:,:,1)-fluor_roi_stim_act_exp(:,:,2))./(fluor_roi_stim_act_exp(:,:,1)+fluor_roi_stim_act_exp(:,:,2));

figure; 
subplot(1,5,1:4);
plot(nanmedian(fluor_change_idx,2),'k','linewidth',2)
set(gca,'XTick',1:numel(wf_roi),'XTickLabel',{wf_roi.area});
ylim([-1,1]);
line(xlim,[0,0]);
ylabel('(Pre-Post)/(Pre+Post)')
subplot(1,5,5);
plot(nanmedian(mua_change_idx,2),'k','linewidth',2)
set(gca,'XTick',1:4,'XTickLabel',cellfun(@(x) ['Str ' num2str(x)],num2cell(1:4),'uni',false));
ylim([-1,1]);
line(xlim,[0,0]);
ylabel('(Pre-Post)/(Pre+Post)')

%% Get stim responses during task with muscimol

data_fns = { ...
    'trial_activity_vanillaChoiceworldNoRepeats_pre_muscimol', ...
    'trial_activity_vanillaChoiceworldNoRepeats_post_muscimol'};

task_str_kernel = cell(2,1);
task_str_ctxpred_kernel = cell(2,1);
mua_norm_exp = cell(2,1);

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;

    % Keep task > str kernels, task > ctx-pred str norm factor
    task_str_kernel{curr_data} = vertcat(mua_taskpred_k_all{:});
    task_str_ctxpred_kernel{curr_data} = vertcat(mua_ctxpred_taskpred_k_all{:});
    mua_norm_exp{curr_data} = vertcat(mua_norm{:});
    
end

mua_task_k_premuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_kernel{1},mua_norm_exp{1},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

mua_ctxpred_task_k_premuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_ctxpred_kernel{1},mua_norm_exp{1},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

mua_task_k_postmuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_kernel{2},mua_norm_exp{1},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

mua_ctxpred_task_k_postmuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_ctxpred_kernel{2},mua_norm_exp{1},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

% Get task>striatum parameters
n_regressors = length(task_regressor_labels);

% Plot task>striatum kernels
stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = [];
move_col = [0.6,0,0.6;0,0.6,0];
go_col = [0.8,0.8,0.2;0.5,0.5,0.5];
outcome_col = [0,0,0.8;0.8,0,0];
task_regressor_cols = {stim_col,move_col,go_col,outcome_col};
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

figure;
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors  
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);  
        
        curr_kernels = mua_ctxpred_task_k_premuscimol{curr_regressor}(:,:,curr_depth,:);
        n_subregressors = size(mua_ctxpred_task_k_premuscimol{curr_regressor},1);
        col = task_regressor_cols{curr_regressor};
        for curr_subregressor = 1:n_subregressors
            AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
                nanmean(curr_kernels(curr_subregressor,:,:,:),4), ...
                AP_sem(curr_kernels(curr_subregressor,:,:,:),4), ...
                col(curr_subregressor,:),0.5);
        end
        
        xlabel('Time (s)');
        ylabel('Weight');
        title(task_regressor_labels{curr_regressor});
        line([0,0],ylim,'color','k');
        
    end
end
linkaxes(p);
y_scale = 1;
t_scale = 0.5;
line(min(xlim) + [0,t_scale],repmat(min(ylim),2,1),'color','k','linewidth',3);
line(repmat(min(xlim),2,1),min(ylim) + [0,y_scale],'color','k','linewidth',3);


figure;
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors  
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);  
        
        curr_kernels = mua_task_k_postmuscimol{curr_regressor}(:,:,curr_depth,:);
        n_subregressors = size(mua_task_k_postmuscimol{curr_regressor},1);
        col = task_regressor_cols{curr_regressor};
        for curr_subregressor = 1:n_subregressors
            AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
                nanmean(curr_kernels(curr_subregressor,:,:,:),4), ...
                AP_sem(curr_kernels(curr_subregressor,:,:,:),4), ...
                col(curr_subregressor,:),0.5);
        end
        
        xlabel('Time (s)');
        ylabel('Weight');
        title(task_regressor_labels{curr_regressor});
        line([0,0],ylim,'color','k');
        
    end
end
linkaxes(p);
y_scale = 1;
t_scale = 0.5;
line(min(xlim) + [0,t_scale],repmat(min(ylim),2,1),'color','k','linewidth',3);
line(repmat(min(xlim),2,1),min(ylim) + [0,y_scale],'color','k','linewidth',3);

use_stim_time = [0.05,0.15];
use_regressor_time = task_regressor_t_shifts{1} >= use_stim_time(1) & ...
    task_regressor_t_shifts{1} <= use_stim_time(2);

a = permute(nanmean(mua_task_k_premuscimol{1}(:,use_regressor_time,:,:),2),[1,3,4,2]);
b = permute(nanmean(mua_task_k_postmuscimol{1}(:,use_regressor_time,:,:),2),[1,3,4,2]);
c = permute(nanmean(mua_ctxpred_task_k_premuscimol{1}(:,use_regressor_time,:,:),2),[1,3,4,2]);
d = permute(nanmean(mua_ctxpred_task_k_postmuscimol{1}(:,use_regressor_time,:,:),2),[1,3,4,2]);

stim = unique([0.06,0.125,0.25,0.5,1].*[-1;1]);

figure; hold on;
plot(stim,nanmean(a(:,1,:),3))
plot(stim,nanmean(b(:,1,:),3))
plot(stim,nanmean(c(:,1,:),3))
plot(stim,nanmean(d(:,1,:),3))

% (do something here to see if the change in striatum is proportional to
% the change in associated cortex)








