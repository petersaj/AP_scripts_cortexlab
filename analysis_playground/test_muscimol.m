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

animal = 'AP053';
protocol = 'vanillaChoiceworld';
flexible_name = true;
experiments = AP_find_experiments(animal,protocol,flexible_name);

experiments = experiments([experiments.imaging] & [experiments.ephys]);

init_array = cell(size(experiments));

curr_day = 1;%length(experiments);
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
            stim_rxn_time_gocue = accumarray(trial_stim_idx(response_trials)',trial_move_t' > 0.5,[11,1],@nanmean,nan);
            
            % Performance (note that this excludes repeat on incorrect trials)
            performance = block.events.sessionPerformanceValues(:,end-10:end);
            
            % Get whether all contrasts were used
            use_all_contrasts = all(block.events.useContrastsValues(end-10:end));
            
            
            % Package final
            bhv(curr_animal).frac_orient_right{curr_day}(curr_exp,:) = performance(3,:)./performance(2,:);
            bhv(curr_animal).rxn_time{curr_day}(curr_exp,:) = stim_rxn_time;
            bhv(curr_animal).rxn_time_gocue{curr_day}(curr_exp,:) = stim_rxn_time_gocue;
            
        end
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
end

% Save for quick loading later
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\tests\muscimol_test';
save_fn = [save_path filesep 'muscimol_bhv.mat'];
save(save_fn,'bhv');
disp(['Saved ' save_fn]);

% Plot
frac_orient_right_cat = cell2mat(permute(cat(2,bhv.frac_orient_right),[1,3,2]));
rxn_time_cat = cell2mat(permute(cat(2,bhv.rxn_time),[1,3,2]));
rxn_time_gocue_cat = cell2mat(permute(cat(2,bhv.rxn_time_gocue),[1,3,2]));

contrastsides = unique([1,0.5,0.25,0.125,0.06,0].*[-1;1]);

figure;
subplot(2,2,1); hold on;
plot(contrastsides,squeeze(frac_orient_right_cat(1,:,:)),'color',[0.5,0.5,0.5],'linewidth',1);
plot(contrastsides,squeeze(frac_orient_right_cat(2,:,:)),'color',[1,0.5,0.5],'linewidth',1);
plot(contrastsides,nanmean(frac_orient_right_cat(1,:,:),3),'color','k','linewidth',3);
plot(contrastsides,nanmean(frac_orient_right_cat(2,:,:),3),'color','r','linewidth',3);
line(xlim,[0.5,0.5],'color','k','linestyle','--');
line([0,0],[0,1],'color','k','linestyle','--');
xlabel('Side*Contrast');
ylabel('Fraction orient right')

subplot(2,2,2); hold on;
plot(contrastsides,squeeze(rxn_time_cat(1,:,:)),'color',[0.5,0.5,0.5],'linewidth',1);
plot(contrastsides,squeeze(rxn_time_cat(2,:,:)),'color',[1,0.5,0.5],'linewidth',1);
plot(contrastsides,nanmean(rxn_time_cat(1,:,:),3),'color','k','linewidth',3);
plot(contrastsides,nanmean(rxn_time_cat(2,:,:),3),'color','r','linewidth',3);
line([-1,1],[0.5,0.5],'color','k','linestyle','--');
xlabel('Side*Contrast');
ylabel('Reaction time');

subplot(2,2,3); hold on;
plot(contrastsides,squeeze(diff(frac_orient_right_cat,[],1)),'color',[0.5,0.5,0.5],'linewidth',1);
plot(contrastsides,nanmean(diff(frac_orient_right_cat,[],1),3),'color','k','linewidth',3);
line(xlim,[0,0],'color','k','linestyle','--');
line([0,0],[-1,1],'color','k','linestyle','--');
xlabel('Side*Contrast');
ylabel('\Delta Fraction orient right')

subplot(2,2,4); hold on;
plot(contrastsides,squeeze(diff(rxn_time_cat,[],1)),'color',[0.5,0.5,0.5],'linewidth',1);
plot(contrastsides,nanmean(diff(rxn_time_cat,[],1),3),'color','k','linewidth',3);
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
AP_image_scroll(vfs_diff)
axis image
colormap(brewermap([],'*RdBu'));
caxis([-1,1]);





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

%% Get muscimol change in passive

data_fns = { ...
    'trial_activity_AP_lcrGratingPassive_pre_muscimol', ...
    'trial_activity_AP_lcrGratingPassive_post_muscimol'};

stimIDs = cell(2,1);
mua_muscimol = cell(2,1);
mua_ctxpred_muscimol = cell(2,1);
fluor_muscimol = cell(2,1);
fluor_roi_muscimol = cell(2,1);

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Split by experiment
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    
    % Get trials with movement during stim to exclude
    wheel_thresh = 0.025;
    quiescent_trials = ~any(abs(wheel_allcat(:,t >= -0.5 & t <= 0.5)) > wheel_thresh,2);
    quiescent_trials_exp = mat2cell(quiescent_trials,trials_recording,1);
    
    % Get stim by experiment
    trial_stim_allcat_exp = mat2cell(trial_stim_allcat,trials_recording,1);
    
    % Get fluorescence in ROIs deconv but not baseline-subtracted
    fluor_deconv_exp = mat2cell(fluor_allcat_deconv,trials_recording,length(t),n_vs);
    fluor_roi_deconv_exp = mat2cell(fluor_roi_deconv,trials_recording,length(t),numel(wf_roi));
    mua_allcat_exp = mat2cell(mua_allcat,trials_recording,length(t),n_depths);
    mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,trials_recording,length(t),n_depths);
    
    % Grab data
    stimIDs{curr_data} = cellfun(@(data,trials) data(trials,:,:), ...
        trial_stim_allcat_exp,quiescent_trials_exp,'uni',false);
    
    mua_muscimol{curr_data} = cellfun(@(data,trials) data(trials,:,:), ...
        mua_allcat_exp,quiescent_trials_exp,'uni',false);
    
    mua_ctxpred_muscimol{curr_data} = cellfun(@(data,trials) data(trials,:,:), ...
        mua_ctxpred_allcat_exp,quiescent_trials_exp,'uni',false);
    
    fluor_muscimol{curr_data} = cellfun(@(data,trials) data(trials,:,:), ...
        fluor_deconv_exp,quiescent_trials_exp,'uni',false);
    
    fluor_roi_muscimol{curr_data} = cellfun(@(data,trials) data(trials,:,:), ...
        fluor_roi_deconv_exp,quiescent_trials_exp,'uni',false);
    
end

% Plot average fluorescence
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');

use_stim = 1;

fluor_premuscimol_mean = ...
    svdFrameReconstruct(U_master(:,:,1:200), ...
    permute(nanmean(cell2mat(cellfun(@(x,stim) ...
    nanmean(x(stim == use_stim,:,:),1),fluor_muscimol{1},stimIDs{1},'uni',false)),1),[3,2,1]));
fluor_postmuscimol_mean = ...
    svdFrameReconstruct(U_master(:,:,1:200), ...
    permute(nanmean(cell2mat(cellfun(@(x,stim) ...
    nanmean(x(stim == use_stim,:,:),1),fluor_muscimol{2},stimIDs{2},'uni',false)),1),[3,2,1]));

AP_image_scroll([fluor_premuscimol_mean,fluor_postmuscimol_mean]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
axis image;
colormap(brewermap([],'*RdBu'));

% Get pre/post stim response
use_stim = 3;

mua_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_muscimol{1},stimIDs{1},'uni',false));
mua_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_muscimol{2},stimIDs{2},'uni',false));

mua_ctxpred_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_ctxpred_muscimol{1},stimIDs{1},'uni',false));
mua_ctxpred_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_ctxpred_muscimol{2},stimIDs{2},'uni',false));

fluor_roi_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_roi_muscimol{1},stimIDs{1},'uni',false));
fluor_roi_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_roi_muscimol{2},stimIDs{2},'uni',false));

% Plot all str responses
figure;
for curr_str = 1:n_depths
    subplot(n_depths,1,curr_str);
    AP_errorfill(t,nanmean(mua_premuscimol_mean(:,:,curr_str),1), ...
    AP_sem(mua_premuscimol_mean(:,:,curr_str),1),'k',1,false);
    AP_errorfill(t,nanmean(mua_postmuscimol_mean(:,:,curr_str),1), ...
    AP_sem(mua_postmuscimol_mean(:,:,curr_str),1),'r',1,false);
end
linkaxes(get(gcf,'Children'))

% Plot ctx responses
figure;
plot_ctx = [1,3,7];
for curr_ctx_idx = 1:length(plot_ctx)
    curr_ctx = plot_ctx(curr_ctx_idx);
    subplot(length(plot_ctx),1,curr_ctx_idx);
    AP_errorfill(t,nanmean(fluor_roi_premuscimol_mean(:,:,curr_ctx),1), ...
        AP_sem(fluor_roi_premuscimol_mean(:,:,curr_ctx),1),'k',1,false);
    AP_errorfill(t,nanmean(fluor_roi_postmuscimol_mean(:,:,curr_ctx),1), ...
        AP_sem(fluor_roi_postmuscimol_mean(:,:,curr_ctx),1),'r',1,false);
    ylabel(wf_roi(curr_ctx).area);
end
linkaxes(get(gcf,'Children'))

% Plot pre/post muscimol for pair of str/ctx
plot_str = 1;
plot_ctx = 3;

figure;
subplot(2,3,1);hold on;
AP_errorfill(t,nanmean(mua_premuscimol_mean(:,:,plot_str),1), ...
    AP_sem(mua_premuscimol_mean(:,:,plot_str),1),'k');
AP_errorfill(t,nanmean(mua_postmuscimol_mean(:,:,plot_str),1), ...
    AP_sem(mua_postmuscimol_mean(:,:,plot_str),1),'r');
xlim([-0.2,1])
xlabel('Time from stim (s)')
ylabel(['Str ' num2str(plot_str)]);
axis square

subplot(2,3,2);hold on;
AP_errorfill(t,nanmean(fluor_roi_premuscimol_mean(:,:,plot_ctx),1), ...
    AP_sem(fluor_roi_premuscimol_mean(:,:,plot_ctx),1),'k');
AP_errorfill(t,nanmean(fluor_roi_postmuscimol_mean(:,:,plot_ctx),1), ...
    AP_sem(fluor_roi_postmuscimol_mean(:,:,plot_ctx),1),'r');
xlim([-0.2,1])
xlabel('Time from stim (s)')
ylabel(wf_roi(plot_ctx).area);
axis square

subplot(2,3,5);hold on;
AP_errorfill(t,nanmean(mua_ctxpred_premuscimol_mean(:,:,plot_str),1), ...
    AP_sem(mua_ctxpred_premuscimol_mean(:,:,plot_str),1),'k');
AP_errorfill(t,nanmean(mua_ctxpred_postmuscimol_mean(:,:,plot_str),1), ...
    AP_sem(mua_ctxpred_postmuscimol_mean(:,:,plot_str),1),'r');
xlim([-0.2,1])
xlabel('Time from stim (s)')
ylabel(['Ctx-pred: Str ' num2str(plot_str)]);
axis square

% Get average responses
t_stim = t >= 0.05 & t <= 0.15;
mua_avg_premuscimol = permute(nanmean(mua_premuscimol_mean(:,t_stim,:),2),[1,3,2]);
mua_avg_postmuscimol = permute(nanmean(mua_postmuscimol_mean(:,t_stim,:),2),[1,3,2]);
mua_avg_postpre_change = (mua_avg_postmuscimol-mua_avg_premuscimol)./(mua_avg_premuscimol);

mua_ctxpred_avg_premuscimol = permute(nanmean(mua_ctxpred_premuscimol_mean(:,t_stim,:),2),[1,3,2]);
mua_ctxpred_avg_postmuscimol = permute(nanmean(mua_ctxpred_postmuscimol_mean(:,t_stim,:),2),[1,3,2]);
mua_ctxpred_avg_postpre_change = (mua_ctxpred_avg_postmuscimol-mua_ctxpred_avg_premuscimol)./(mua_ctxpred_avg_premuscimol);

fluor_avg_premuscimol = permute(nanmean(fluor_roi_premuscimol_mean(:,t_stim,:),2),[1,3,2]);
fluor_avg_postmuscimol = permute(nanmean(fluor_roi_postmuscimol_mean(:,t_stim,:),2),[1,3,2]);
fluor_avg_postpre_change = (fluor_avg_postmuscimol-fluor_avg_premuscimol)./(fluor_avg_premuscimol);

subplot(2,3,3);
plot(fluor_avg_postpre_change(:,plot_ctx),mua_avg_postpre_change(:,plot_str),'.k','MarkerSize',20)
xlabel(wf_roi(plot_ctx).area);
ylabel(['Str ' num2str(plot_str)]);
line([-1,1],[-1,1],'color','k');
line([-1,1],[0,0],'color','k','linestyle','--');
line([0,0],[-1,1],'color','k','linestyle','--');
axis square;

nonan_points = ~isnan(mua_avg_postpre_change(:,plot_str)) & ...
    ~isnan(fluor_avg_postpre_change(:,plot_ctx));
[r,p] = corrcoef(fluor_avg_postpre_change(nonan_points,plot_ctx), ...
    mua_avg_postpre_change(nonan_points,plot_str));

subplot(2,3,6);
plot(mua_ctxpred_avg_postpre_change(:,plot_str),mua_avg_postpre_change(:,plot_str),'.k','MarkerSize',20)
ylabel(['Str ' num2str(plot_str)]);
xlabel(['Ctx-pred: Str ' num2str(plot_str)]);
line([-1,1],[-1,1],'color','k');
line([-1,1],[0,0],'color','k','linestyle','--');
line([0,0],[-1,1],'color','k','linestyle','--');
axis square;

nonan_points = ~isnan(mua_ctxpred_avg_postpre_change(:,plot_str)) & ...
    ~isnan(mua_avg_postpre_change(:,plot_str));
[r,p] = corrcoef(mua_ctxpred_avg_postpre_change(nonan_points,plot_str), ...
    mua_avg_postpre_change(nonan_points,plot_str));



%%%%%%% TESTING

% fluor_premuscimol_mean = cell2mat(permute(cellfun(@(x) ...
%     svdFrameReconstruct(U_master(:,:,1:200),permute(x,[3,2,1])), ...
%     cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1), ...
%     fluor_muscimol{1},stimIDs{1},'uni',false),'uni',false),[2,3,4,1]));
% 
% fluor_postmuscimol_mean = cell2mat(permute(cellfun(@(x) ...
%     svdFrameReconstruct(U_master(:,:,1:200),permute(x,[3,2,1])), ...
%     cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1), ...
%     fluor_muscimol{2},stimIDs{2},'uni',false),'uni',false),[2,3,4,1]));
% 
% fluor_avg_premuscimol = squeeze(nanmean(fluor_premuscimol_mean(:,:,t_stim,:),3));
% fluor_avg_postmuscimol = squeeze(nanmean(fluor_postmuscimol_mean(:,:,t_stim,:),3));



%% Get task stim kernels with muscimol

data_fns = { ...
    'trial_activity_vanillaChoiceworldNoRepeats_pre_muscimol', ...
    'trial_activity_vanillaChoiceworldNoRepeats_post_muscimol'};

mua_norm_exp = cell(2,1);
task_str_kernel = cell(2,1);
task_str_ctxpred_kernel = cell(2,1);
task_ctx_roi_kernel = cell(2,1);

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Get average task > cortex kernels (ROIs)
    n_regressors = length(task_regressor_labels);
    regressor_roi = cell(n_regressors,1);
    for curr_regressor = 1:n_regressors
        
        curr_k = cell2mat(cellfun(@(x) x{curr_regressor}, ...
            permute(vertcat(fluor_taskpred_k_all{:}),[2,3,4,1]),'uni',false));
        
        curr_k_roi = nan(size(curr_k,1),size(curr_k,2),n_rois,size(curr_k,4));
        for curr_subregressor = 1:size(curr_k,1)
            for curr_exp = 1:size(curr_k,4)
                curr_k_roi(curr_subregressor,:,:,curr_exp) = ...
                    permute(AP_svd_roi(U_master(:,:,1:n_vs), ...
                    permute(curr_k(curr_subregressor,:,:,curr_exp),[3,2,1]), ...
                    [],[],cat(3,wf_roi.mask)),[3,2,1]);
            end
        end
        
        regressor_roi{curr_regressor} = curr_k_roi;
        
        AP_print_progress_fraction(curr_regressor,n_regressors);
    end
    
    % Keep task > str kernels, task > ctx-pred str norm factor
    mua_norm_exp{curr_data} = vertcat(mua_norm{:});
    task_str_kernel{curr_data} = vertcat(mua_taskpred_k_all{:});
    task_str_ctxpred_kernel{curr_data} = vertcat(mua_ctxpred_taskpred_k_all{:});
    task_ctx_roi_kernel{curr_data} = regressor_roi;
    
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
    task_str_kernel{2},mua_norm_exp{2},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

mua_ctxpred_task_k_postmuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_ctxpred_kernel{2},mua_norm_exp{2},'uni',false), ...
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
        
        curr_kernels = mua_ctxpred_task_k_postmuscimol{curr_regressor}(:,:,curr_depth,:);
        n_subregressors = size(mua_ctxpred_task_k_postmuscimol{curr_regressor},1);
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
        
        curr_kernels = mua_task_k_premuscimol{curr_regressor}(:,:,curr_depth,:);
        n_subregressors = size(mua_task_k_postmuscimol{curr_regressor},1);
        col = task_regressor_cols{curr_regressor};
        for curr_subregressor = 1:n_subregressors
            AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
                nanmean(curr_kernels(curr_subregressor,:,:,:),4), ...
                AP_sem(curr_kernels(curr_subregressor,:,:,:),4), ...
                col(curr_subregressor,:));
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

task_ctx_stim_premuscimol = permute(sum(task_ctx_roi_kernel{1}{1},2),[1,3,4,2]);
task_ctx_stim_postmuscimol = permute(sum(task_ctx_roi_kernel{2}{1},2),[1,3,4,2]);
task_ctx_stim_change = (task_ctx_stim_postmuscimol-task_ctx_stim_premuscimol)./ ...
    (abs(task_ctx_stim_postmuscimol)+abs(task_ctx_stim_premuscimol));

task_mua_stim_premuscimol = permute(sum(mua_task_k_premuscimol{1},2),[1,3,4,2]);
task_mua_stim_postmuscimol = permute(sum(mua_task_k_postmuscimol{1},2),[1,3,4,2]);
task_mua_stim_change = (task_mua_stim_postmuscimol-task_mua_stim_premuscimol)./ ...
    (abs(task_mua_stim_postmuscimol)+abs(task_mua_stim_premuscimol));


plot_str = 1;
plot_ctx = 3;

stim = unique([0.06,0.125,0.25,0.5,1].*[-1;1]);
stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = [];

figure;
p = gobjects(2,3);

p(1,1) = subplot(2,3,1);
curr_regressor = 1;
curr_kernels = mua_task_k_premuscimol{curr_regressor}(:,:,plot_str,:);
n_subregressors = size(mua_task_k_premuscimol{curr_regressor},1);
for curr_subregressor = 1:length(stim)
    AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
        nanmean(curr_kernels(curr_subregressor,:,:,:),4), ...
        AP_sem(curr_kernels(curr_subregressor,:,:,:),4), ...
        stim_col(curr_subregressor,:),0.5);
end
xlabel('Time (s)');
ylabel('Weight');
ylabel(['Str ' num2str(plot_str) ' weight']);
title('Pre-muscimol');

p(2,1) = subplot(2,3,4);
curr_regressor = 1;
curr_kernels = mua_task_k_postmuscimol{curr_regressor}(:,:,plot_str,:);
n_subregressors = size(mua_task_k_postmuscimol{curr_regressor},1);
for curr_subregressor = 1:length(stim)
    AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
        nanmean(curr_kernels(curr_subregressor,:,:,:),4), ...
        AP_sem(curr_kernels(curr_subregressor,:,:,:),4), ...
        stim_col(curr_subregressor,:),0.5);
end
xlabel('Time (s)');
ylabel(['Str ' num2str(plot_str) ' weight']);
title('Post-muscimol');

p(1,2) = subplot(2,3,2);
curr_regressor = 1;
curr_kernels = task_ctx_roi_kernel{1}{curr_regressor}(:,:,plot_ctx,:);
n_subregressors = size(task_ctx_roi_kernel{1}{curr_regressor},1);
for curr_subregressor = 1:length(stim)
    AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
        nanmean(curr_kernels(curr_subregressor,:,:,:),4), ...
        AP_sem(curr_kernels(curr_subregressor,:,:,:),4), ...
        stim_col(curr_subregressor,:),0.5);
end
xlabel('Time (s)');
ylabel([wf_roi(plot_ctx).area 'weight']);
title('Pre-muscimol');

p(2,2) = subplot(2,3,5);
curr_regressor = 1;
curr_kernels = task_ctx_roi_kernel{2}{curr_regressor}(:,:,plot_ctx,:);
n_subregressors = size(task_ctx_roi_kernel{2}{curr_regressor},1);
for curr_subregressor = 1:length(stim)
    AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
        nanmean(curr_kernels(curr_subregressor,:,:,:),4), ...
        AP_sem(curr_kernels(curr_subregressor,:,:,:),4), ...
        stim_col(curr_subregressor,:),0.5);
end
xlabel('Time (s)');
ylabel([wf_roi(plot_ctx).area 'weight']);
title('Post-muscimol');

linkaxes([p(1,1),p(2,1)],'y');
linkaxes([p(1,2),p(2,2)],'y');

subplot(2,3,3);
plot(squeeze(nanmean(task_ctx_stim_change(stim > 0,plot_ctx,:),1)), ...
    squeeze(nanmean(task_mua_stim_change(stim > 0,plot_str,:),1)),'.k','MarkerSize',10);
line([0,0],[-1,1],'color','k','linestyle','--');
line([-1,1],[0,0],'color','k','linestyle','--');
line([-1,1],[-1,1],'color','k');





a = cellfun(@(x) permute(sum(x,2),[1,3,4,2]),mua_task_k_premuscimol,'uni',false);
b = cellfun(@(x) permute(sum(x,2),[1,3,4,2]),mua_task_k_postmuscimol,'uni',false);

figure;
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);
        hold on
        
        AP_errorfill([],nanmean(a{curr_regressor}(:,curr_depth,:),3), ...
            AP_sem(a{curr_regressor}(:,curr_depth,:),3),'k');
        AP_errorfill([],nanmean(b{curr_regressor}(:,curr_depth,:),3), ...
            AP_sem(b{curr_regressor}(:,curr_depth,:),3),'r');
        
    end
end

linkaxes(p,'y');




%% Cortex-explained striatal variance pre/post muscimol

data_fns = { ...
    'trial_activity_vanillaChoiceworldNoRepeats_pre_muscimol', ...
    'trial_activity_vanillaChoiceworldNoRepeats_post_muscimol'};

figure;
for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Split data by recording
    trials_allcat = size(wheel_allcat,1);
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    use_split = trials_recording;
    split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
        [1:length(use_split)]',reshape(use_split,[],1),'uni',false));
    
    % Get R^2 for task, cortex full, and cortex ROI predictions
    taskpred_r2 = nan(max(split_idx),n_depths);
    ctxpred_r2 = nan(max(split_idx),n_depths);
    ctxroipred_r2 = nan(max(split_idx),n_depths);
    for curr_exp = 1:max(split_idx)
        
        curr_data = reshape(permute(mua_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
        curr_taskpred_data = reshape(permute(mua_taskpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
        curr_ctxpred_data = reshape(permute(mua_ctxpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
        
        % Set common NaNs
        nan_samples = isnan(curr_data) | isnan(curr_taskpred_data) | isnan(curr_ctxpred_data);
        curr_data(nan_samples) = NaN;
        curr_taskpred_data(nan_samples) = NaN;
        curr_ctxpred_data(nan_samples) = NaN;
        
        taskpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_taskpred_data).^2,1)./ ...
            nansum((curr_data-nanmean(curr_data,1)).^2,1));
        ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
            nansum((curr_data-nanmean(curr_data,1)).^2,1));
    end
    
    subplot(2,1,curr_data);
    errorbar(nanmean(taskpred_r2,1),AP_sem(taskpred_r2,1),'b','linewidth',2,'CapSize',0);
    errorbar(nanmean(ctxpred_r2,1),AP_sem(ctxpred_r2,1),'color',[0,0.6,0],'linewidth',2,'CapSize',0);
    xlabel('Striatum depth');
    ylabel('Task explained variance');
    legend({'Task','Cortex'});
    
    
    
end


%% Offset between measured and predicted pre/post muscimol

data_fns = { ...
    'trial_activity_vanillaChoiceworldNoRepeats_pre_muscimol', ...
    'trial_activity_vanillaChoiceworldNoRepeats_post_muscimol'};

ctx_v_error = cell(2,2);
measured_pred_error_fig = figure;
measured_v_pred_fig = figure;
for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Split data by recording
    trials_allcat = size(wheel_allcat,1);
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    use_split = trials_recording;
    split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
        [1:length(use_split)]',reshape(use_split,[],1),'uni',false));
    
    % Set alignment shifts
    t_leeway = -t(1);
    leeway_samples = round(t_leeway*(sample_rate));
    stim_align = zeros(size(mua_allcat,1),1);
    move_align = -move_idx + leeway_samples;
    outcome_align = -outcome_idx + leeway_samples;
    
    % Set windows to average activity
    timeavg_labels = {'Stim','Move early','Move late','Outcome'};
    timeavg_t = {[0.05,0.15],[-0.05,0.05],[0.05,0.15],[0.05,0.15]};
    timeavg_align = {stim_align,move_align,move_align,outcome_align};
    timeavg_trial_conditions = ...
        {[trial_contrastside_allcat > 0,trial_contrastside_allcat == 0,trial_contrastside_allcat < 0], ...
        [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
        [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
        [trial_outcome_allcat == 1, trial_outcome_allcat == -1]};
    timeavg_condition_colors = ...
        {[1,0,0;0,0,0;0,0,1], ...
        [0.6,0,0.6;0,0.6,0], ...
        [0.6,0,0.6;0,0.6,0], ...
        [0,0,0.7;0,0,0]};
    timeavg_event_offset = {'Stim','Move onset','Move onset','Outcome'};
    
    % Set activity percentiles and bins
    act_prctile = [5,95];
    n_act_bins = 5;
    
    % Set areas and conditions
    plot_areas = [1,2,3,4];
    
    % Loop across area pairs, plot binned predicted v measured activity
    curr_act_allcat = mua_allcat;
    curr_act_pred_allcat = mua_ctxpred_allcat;

    % (old: to use all events)
    % task_fix = mua_taskpred_allcat - mua_ctxpred_taskpred_allcat;
    % curr_act_pred_fix_allcat = curr_act_pred_allcat + task_fix;
    
    % (to use individual events)
    task_fix = (mua_taskpred_allcat - mua_taskpred_reduced_allcat) - ...
        (mua_ctxpred_taskpred_allcat - mua_ctxpred_taskpred_reduced_allcat);
    curr_act_pred_fix_allcat = curr_act_pred_allcat + task_fix;
    
    for curr_area_idx = 1:length(plot_areas)
        
        plot_area = plot_areas(curr_area_idx);
        
        for curr_timeavg = 1:length(timeavg_labels)
            
            trial_conditions = timeavg_trial_conditions{curr_timeavg};
            curr_event_offset_idx = find(strcmp(timeavg_event_offset{curr_timeavg},task_regressor_labels));
            curr_act_pred_fix_allcat_event = ...
                curr_act_pred_fix_allcat(:,:,:,curr_event_offset_idx);
            
            % (re-align and split activity)
            act_title = timeavg_labels{curr_timeavg};
            
            curr_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_act_allcat(trial,:,plot_area), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_allcat,1)),'uni',false)), ...
                use_split,length(t));
            
            curr_act_pred = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_act_pred_allcat(trial,:,plot_area), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
                use_split,length(t));
            
            curr_act_predfix = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_act_pred_fix_allcat_event(trial,:,plot_area), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_fix_allcat_event,1)),'uni',false)), ...
                use_split,length(t));
            
            trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
            
            % (get average activity within window)
            curr_event_t = t >= timeavg_t{curr_timeavg}(1) & t <= timeavg_t{curr_timeavg}(2);
            curr_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act,'uni',false);
            curr_act_pred_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_pred,'uni',false);
            curr_act_predfix_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_predfix,'uni',false);

            %%% Plot total error in prediction           
            curr_meas_pred_error = cellfun(@(x,y) x-y,curr_act_avg,curr_act_pred_avg,'uni',false);
            curr_meas_pred_error_cond = cell2mat(cellfun(@(x,cond) arrayfun(@(curr_cond) ...
                nanmean(x(cond(:,curr_cond))),1:size(cond,2)),curr_meas_pred_error,trial_conditions_exp,'uni',false));
                   
            figure(measured_pred_error_fig);
            
            subplot(length(plot_areas),length(timeavg_labels), ...
                sub2ind(fliplr([length(plot_areas),length(timeavg_labels)]),curr_timeavg,curr_area_idx));
            hold on;
            
            switch curr_data
                case 1
                    col = 'k';
                case 2
                    col = 'r';
            end
            
            errorbar(nanmean(curr_meas_pred_error_cond,1),AP_sem(curr_meas_pred_error_cond,1),col,'linewidth',2);
            ylabel('Measured-Predicted');
            xlabel('Condition');
            
            drawnow;
            
            %%% Get and plot offset for Str 1 stim only
            
            if curr_area_idx == 1 &&  curr_timeavg == 1
                                
                % (bin predicted data across percentile range)
                pred_bin_edges = prctile(cell2mat(curr_act_pred_avg),linspace(act_prctile(1),act_prctile(2),n_act_bins+1));
                pred_bin_centers = pred_bin_edges(1:end-1) + diff(pred_bin_edges)./2;
                pred_trial_bins = cellfun(@(x) discretize(x,pred_bin_edges),curr_act_pred_avg,'uni',false);
                
                % (get activity binned by predicted)
                pred_use_trials = cellfun(@(act,act_pred,trial_bins) ...
                    squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
                    squeeze(~any(isnan(act_pred(:,curr_event_t)),2)) & ...
                    ~isnan(trial_bins), ...
                    curr_act,curr_act_pred,pred_trial_bins,'uni',false);
                
                act_predbinmean = cell2mat(arrayfun(@(condition) ...
                    cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                    accumarray(bins(use_trials & trial_cond(:,condition)), ...
                    act(use_trials & trial_cond(:,condition)), ...
                    [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                    curr_act_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
                    permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
                
                act_pred_predbinmean = cell2mat(arrayfun(@(condition) ...
                    cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                    accumarray(bins(use_trials & trial_cond(:,condition)), ...
                    act(use_trials & trial_cond(:,condition)), ...
                    [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                    curr_act_pred_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
                    permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
                
                % (get "fixed" predicted activity binned by predicted)
                act_predfix_predbinmean = cell2mat(arrayfun(@(condition) ...
                    cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                    accumarray(bins(use_trials & trial_cond(:,condition)), ...
                    act(use_trials & trial_cond(:,condition)), ...
                    [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                    curr_act_predfix_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
                    permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
                
                % Get condition difference significance from shuffle
                n_shuff = 1000;
                trial_conditions_shuff = cellfun(@(x) AP_shake(x,1), ...
                    repmat(trial_conditions_exp,1,n_shuff),'uni',false);
                
                act_predbinmean_condshuff = cell2mat(arrayfun(@(curr_shuff) ...
                    cell2mat(arrayfun(@(condition) ...
                    cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                    accumarray(bins(use_trials & trial_cond(:,condition)), ...
                    act(use_trials & trial_cond(:,condition)), ...
                    [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                    curr_act_avg,pred_trial_bins,trial_conditions_shuff(:,curr_shuff), ...
                    pred_use_trials,'uni',false)'), ...
                    permute(1:size(trial_conditions,2),[1,3,2]),'uni',false)), ...
                    permute(1:n_shuff,[1,3,4,2]),'uni',false));
                
                act_predfix_predbinmean_condshuff = cell2mat(arrayfun(@(curr_shuff) ...
                    cell2mat(arrayfun(@(condition) ...
                    cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                    accumarray(bins(use_trials & trial_cond(:,condition)), ...
                    act(use_trials & trial_cond(:,condition)), ...
                    [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                    curr_act_predfix_avg,pred_trial_bins,trial_conditions_shuff(:,curr_shuff), ...
                    pred_use_trials,'uni',false)'), ...
                    permute(1:size(trial_conditions,2),[1,3,2]),'uni',false)), ...
                    permute(1:n_shuff,[1,3,4,2]),'uni',false));
                
                % (measured data: null = no difference between conditions)
                cond_combos = nchoosek(1:size(trial_conditions,2),2);
                cond_sig_diff = false(n_act_bins,size(cond_combos,1));
                predfix_cond_sig_diff = false(n_act_bins,size(cond_combos,1));
                for curr_cond_combo = 1:size(cond_combos,1)
                    curr_combo_diff = nanmean(act_predbinmean(:,:,cond_combos(curr_cond_combo,1)) - ...
                        act_predbinmean(:,:,cond_combos(curr_cond_combo,2)),2);
                    shuff_prctile = squeeze(prctile(nanmean(act_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,1),:) - ...
                        act_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,2),:),2),[2.5,97.5],4));
                    cond_sig_diff(:,curr_cond_combo) = ...
                        curr_combo_diff < shuff_prctile(:,1) | ...
                        curr_combo_diff > shuff_prctile(:,2);
                    
                    predfix_curr_combo_diff = nanmean(act_predfix_predbinmean(:,:,cond_combos(curr_cond_combo,1)) - ...
                        act_predfix_predbinmean(:,:,cond_combos(curr_cond_combo,2)),2);
                    predfix_shuff_prctile = squeeze(prctile(nanmean(act_predfix_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,1),:) - ...
                        act_predfix_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,2),:),2),[2.5,97.5],4));
                    predfix_cond_sig_diff(:,curr_cond_combo) = ...
                        predfix_curr_combo_diff < predfix_shuff_prctile(:,1) | ...
                        predfix_curr_combo_diff > predfix_shuff_prctile(:,2);
                end
                
                % Plot binned measured, predicted, and error (by predicted bins)
                measured_pred_fig = figure('color','w','Name', ...
                    ['Str ' num2str(plot_area)' ', ' timeavg_labels{curr_timeavg}]);
                n_col_bins = n_act_bins + 2;
                
                [binned_act_pred_t,binned_act_pred_grp] = grpstats(cell2mat(curr_act_pred), ...
                    [cell2mat(pred_trial_bins),trial_conditions],{'nanmean','gname'});
                binned_act_pred_grp = cellfun(@str2num,binned_act_pred_grp);
                
                [binned_act_t,binned_act_grp] = grpstats(cell2mat(curr_act), ...
                    [cell2mat(pred_trial_bins),trial_conditions],{'nanmean','gname'});
                binned_act_grp = cellfun(@str2num,binned_act_grp);
                
                % (plot predicted data timecourse)
                for curr_cond = 1:size(trial_conditions,2)
                    subplot(2,size(trial_conditions,2), ...
                        sub2ind(fliplr([2,size(trial_conditions,2)]),curr_cond,1)); hold on;
                    
                    curr_mean = nanmean(cell2mat(cellfun(@(act,use_trials,cond) ...
                        act(use_trials & cond(:,curr_cond),:), ...
                        curr_act_pred,pred_use_trials,trial_conditions_exp,'uni',false)),1);
                    curr_std = nanstd(cell2mat(cellfun(@(act,use_trials,cond) ...
                        act(use_trials & cond(:,curr_cond),:), ...
                        curr_act_pred,pred_use_trials,trial_conditions_exp,'uni',false)),[],1);
                    
                    AP_errorfill(t,curr_mean,curr_std,'k',0.5,true);
                    
                    set(gca,'ColorOrder',brewermap(n_col_bins,'*Greens'));
                    plot(t,binned_act_pred_t(binned_act_pred_grp(:,curr_cond + 1) == 1,:)','linewidth',2);
                    xlabel('Time'); ylabel('Predicted data');
                    title(['Condition ' num2str(curr_cond)]);
                    line(repmat(timeavg_t{curr_timeavg}(1),2,1),ylim,'color','k');
                    line(repmat(timeavg_t{curr_timeavg}(2),2,1),ylim,'color','k');
                end
                
                % (plot measured data timecourse)
                for curr_cond = 1:size(trial_conditions,2)
                    subplot(2,size(trial_conditions,2), ...
                        sub2ind(fliplr([2,size(trial_conditions,2)]),curr_cond,2)); hold on;
                    set(gca,'ColorOrder',[brewermap(n_col_bins,'*Greys')]);
                    plot(t,binned_act_t(binned_act_grp(:,curr_cond + 1) == 1,:)','linewidth',2);
                    xlabel('Time'); ylabel('Measured data');
                    title(['Condition ' num2str(curr_cond)]);
                    line(repmat(timeavg_t{curr_timeavg}(1),2,1),ylim,'color','k');
                    line(repmat(timeavg_t{curr_timeavg}(2),2,1),ylim,'color','k');
                end
                
                linkaxes(get(measured_pred_fig,'Children'),'xy');
                
                % Plot measured v binned predicted
                figure(measured_v_pred_fig);
                subplot(1,2,curr_data); hold on;
                set(gca,'ColorOrder',timeavg_condition_colors{curr_timeavg});
                
                fill_cols = min(timeavg_condition_colors{curr_timeavg} + 0.5,1);
                for curr_cond = 1:size(trial_conditions,2)
                    AP_errorfill( ...
                        squeeze(nanmean(act_pred_predbinmean(:,:,curr_cond),2)), ...
                        squeeze(nanmean(act_predfix_predbinmean(:,:,curr_cond),2)), ...
                        squeeze(AP_sem(act_predfix_predbinmean(:,:,curr_cond),2)), ...
                        fill_cols(curr_cond,:),1,false);
                end
                
                errorbar( ...
                    squeeze(nanmean(act_pred_predbinmean,2)), ...
                    squeeze(nanmean(act_predbinmean,2)), ...
                    squeeze(AP_sem(act_predbinmean,2)), ...
                    'linestyle','none','linewidth',2,'CapSize',10);
                
                ctx_col = brewermap(n_col_bins,'*Greens');
                scatter( ...
                    reshape(squeeze(nanmean(act_pred_predbinmean,2)),[],1), ...
                    reshape(squeeze(nanmean(act_predbinmean,2)),[],1), ...
                    60,repmat(ctx_col(1:n_act_bins,:),size(trial_conditions,2),1),'filled');
                
                xlabel(['Predicted (' num2str(plot_area) ')']);
                ylabel(['Measured (' num2str(plot_area) ')'])
                axis square tight;
                xlim(xlim + [-0.1,0.1]);
                title(act_title);
                
                % (plot significant measured condition differences)
                % (* and o = significant in both measured and "fixed" predicted)
                curr_ylim = max(ylim);
                for curr_cond_combo = 1:size(cond_combos,1)
                    % (plot * for measured condition differences)
                    sig_y = curr_ylim + 0.1*curr_cond_combo;
                    sig_x = pred_bin_centers;
                    if any(cond_sig_diff(:,curr_cond_combo))
                        plot(sig_x(cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                            '*','MarkerSize',10,'color', ...
                            timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,1),:));
                        plot(sig_x(cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                            '*','MarkerSize',5,'color', ...
                            timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,2),:));
                    end
                    % (plot o for [predicted condition differences)
                    if any(predfix_cond_sig_diff(:,curr_cond_combo))
                        plot(sig_x(predfix_cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                            'o','MarkerSize',15,'color', ...
                            timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,1),:));
                        plot(sig_x(predfix_cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                            'o','MarkerSize',10,'color', ...
                            timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,2),:));
                    end
                end   
                drawnow;
                
                %%%%%%%%%% TESTING %%%%%%%%%%%%%
                % Store the cortical activity and error
%                 use_cond = 1;                
%                 ctx_v_error{curr_data,1} = ...
%                     cellfun(@(x,cond) nanmean(x(cond(:,use_cond)),1),curr_act_ctx_avg,trial_conditions_exp);
%                 ctx_v_error{curr_data,2} = ...
%                     cellfun(@(x,cond) nanmean(x(cond(:,use_cond)),1),curr_meas_pred_error,trial_conditions_exp);
                % Store error but in bins
                use_cond = 1;     
                curr_offset = nanmean(act_predbinmean(:,:,1) - act_predbinmean(:,:,3),1);
                ctx_v_error{curr_data,1} = ...
                    cellfun(@(x,cond,pred_use_trials) ...
                    nanmean(x(pred_use_trials & cond(:,use_cond)),1), ...
                    curr_act_pred_avg,trial_conditions_exp,pred_use_trials);
                ctx_v_error{curr_data,2} = ...
                    curr_offset';
                
            end
        end
    end
end


% Plot striatum/cortex offset relative to cortical activity
ctx_change = (ctx_v_error{2,1} - ctx_v_error{1,1})./(abs(ctx_v_error{2,1}) + abs(ctx_v_error{1,1}));
error_change = (ctx_v_error{2,2} - ctx_v_error{1,2})./(abs(ctx_v_error{2,2}) + abs(ctx_v_error{1,2}));
figure; 
plot(ctx_change,error_change,'.k','MarkerSize',20);
line([-1,1],[0,0],'color','k','linestyle','--');
line([0,0],[-1,1],'color','k','linestyle','--');
line([-1,1],[-1,1],'color','k');
xlabel('Cortex (pred striatum) change');
ylabel('Error change');

nonan_points = ~isnan(ctx_change) & ~isnan(error_change);
[r,p] = corrcoef(ctx_change(nonan_points),error_change(nonan_points));
disp(r(2));
disp(p(2));

% HAVEN'T GOTTEN ABOVE TO LOOK RIGHT - NOT SURE WHAT'S RIGHT TO DO HERE
% above is biased based on which conditions mouse was shown, maybe just use
% std or something?


%% Cortex/striatum offset pre/post muscimol (passive)

data_fns = { ...
    'trial_activity_AP_lcrGratingPassive_pre_muscimol', ...
    'trial_activity_AP_lcrGratingPassive_post_muscimol'};

ctx_v_error = cell(2,2);

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Split data by recording
    trials_allcat = size(wheel_allcat,1);
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    use_split = trials_recording;
    split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
        [1:length(use_split)]',reshape(use_split,[],1),'uni',false));
       
    % Set alignment shifts
    t_leeway = -t(1);
    leeway_samples = round(t_leeway*(sample_rate));
    stim_align = zeros(size(mua_allcat,1),1);
    
    % Get trials with movement during stim to exclude
    wheel_thresh = 0.025;
    quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > wheel_thresh,2);
    
    % Set windows to average activity
    timeavg_labels = {'Stim'};
    timeavg_t = {[0.05,0.15]};
    timeavg_align = {stim_align};
    timeavg_trial_conditions = ...
        {[trial_contrastside_allcat > 0 & quiescent_trials, ...
        trial_contrastside_allcat < 0 & quiescent_trials]};
    timeavg_condition_colors = ...
        {[1,0,0;0,0,1]};
    
    % Set activity percentiles and bins
    act_prctile = [5,95];
    n_act_bins = 5;
    
    % Set areas and conditions
    plot_areas = [1];
    
    % Loop across area pairs, plot binned predicted v measured activity
    curr_act_allcat = mua_allcat;
    
    % (ctx-predicted)
    curr_act_pred_allcat = mua_ctxpred_allcat;
    
    % (fix by average stim response within experiment)
    mua_allcat_exp = mat2cell(mua_allcat,use_split,length(t),n_depths);
    mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,use_split,length(t),n_depths);
    trial_contrastside_allcat_exp = mat2cell(trial_contrastside_allcat,use_split,1);
    trial_conditions_exp = mat2cell(timeavg_trial_conditions{1},use_split, ...
        size(timeavg_trial_conditions{1},2));
    
    str_stim_exp = nan(length(unique(trial_contrastside_allcat)),length(t),n_depths,length(use_split));
    ctx_stim_exp = nan(length(unique(trial_contrastside_allcat)),length(t),n_depths,length(use_split));
    curr_act_pred_fix_allcat_exp = mua_ctxpred_allcat_exp;
    for curr_exp = 1:length(use_split)
        
        curr_stim = unique(trial_contrastside_allcat_exp{curr_exp});
        for curr_stim_idx = 1:length(curr_stim)
            curr_trials = any(trial_conditions_exp{curr_exp},2) & ...
                trial_contrastside_allcat_exp{curr_exp} == curr_stim(curr_stim_idx);
            curr_act_stim_avg = nanmean(mua_allcat_exp{curr_exp}(curr_trials,:,:),1);
            curr_act_pred_stim_avg = nanmean(mua_ctxpred_allcat_exp{curr_exp}(curr_trials,:,:),1);
            
            str_stim_exp(curr_stim_idx,:,:,curr_exp) = curr_act_stim_avg;
            ctx_stim_exp(curr_stim_idx,:,:,curr_exp) = curr_act_pred_stim_avg;
            
            curr_stim_fix = curr_act_stim_avg - curr_act_pred_stim_avg;
            curr_act_pred_fix_allcat_exp{curr_exp}(curr_trials,:,:) = ...
                curr_act_pred_fix_allcat_exp{curr_exp}(curr_trials,:,:) + curr_stim_fix;
        end
    end
    curr_act_pred_fix_allcat = cell2mat(curr_act_pred_fix_allcat_exp);
    
    % Plot stim response
    stim_col = colormap_BlueWhiteRed(5);
    stim_col(6,:) = 0.5;
    stim_col_contrastside = unique([0,0.06,0.125,0.25,0.5,1].*[-1;1]);
    [~,used_stim_idx] = ismember(unique(trial_contrastside_allcat), ...
        stim_col_contrastside,'rows');
    
    str_stim_avg = nanmean(str_stim_exp,4);
    figure('color','w');
    for curr_area_idx = 1:length(plot_areas)
        curr_area = plot_areas(curr_area_idx);
        subplot(length(plot_areas),2,length(plot_areas)*(curr_area_idx-1)+1);
        hold on;
        for curr_stim_idx = 1:size(str_stim_avg,1)
            AP_errorfill(t, ...
                nanmean(str_stim_exp(curr_stim_idx,:,curr_area,:),4), ...
                AP_sem(str_stim_exp(curr_stim_idx,:,curr_area,:),4), ...
                stim_col(used_stim_idx(curr_stim_idx),:),0.5,true);
        end
        line([0,0],ylim,'color','k');
        ylabel(['Measured (' num2str(curr_area) ')']);
        xlabel('Time from stim');
        
        subplot(length(plot_areas),2,length(plot_areas)*(curr_area_idx-1)+2);
        hold on;
        for curr_stim_idx = 1:size(str_stim_avg,1)
            AP_errorfill(t, ...
                nanmean(ctx_stim_exp(curr_stim_idx,:,curr_area,:),4), ...
                AP_sem(ctx_stim_exp(curr_stim_idx,:,curr_area,:),4), ...
                stim_col(used_stim_idx(curr_stim_idx),:),0.5,true);
        end
        line([0,0],ylim,'color','k');
        ylabel(['Predicted (' num2str(curr_area) ')']);
        xlabel('Time from stim');
    end
    linkaxes(get(gcf,'Children'),'xy');
    
    measured_v_pred_fig = figure('color','w');
    for curr_area_idx = 1:length(plot_areas)
        
        plot_area = plot_areas(curr_area_idx);
        
        for curr_timeavg = 1:length(timeavg_labels)
            
            trial_conditions = timeavg_trial_conditions{curr_timeavg};
            
            % (re-align and split activity)
            act_title = timeavg_labels{curr_timeavg};
            
            curr_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_act_allcat(trial,:,plot_area), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_allcat,1)),'uni',false)), ...
                use_split,length(t));
            
            curr_act_pred = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_act_pred_allcat(trial,:,plot_area), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
                use_split,length(t));
            
            curr_act_predfix = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_act_pred_fix_allcat(trial,:,plot_area), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_fix_allcat,1)),'uni',false)), ...
                use_split,length(t));
            
            trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
            
            % (get average activity within window)
            curr_event_t = t >= timeavg_t{curr_timeavg}(1) & t <= timeavg_t{curr_timeavg}(2);
            curr_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act,'uni',false);
            curr_act_pred_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_pred,'uni',false);
            curr_act_predfix_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_predfix,'uni',false);
            
            % (bin predicted data across percentile range)
            pred_bin_edges = prctile(cell2mat(curr_act_pred_avg),linspace(act_prctile(1),act_prctile(2),n_act_bins+1));
            pred_bin_centers = pred_bin_edges(1:end-1) + diff(pred_bin_edges)./2;
            pred_trial_bins = cellfun(@(x) discretize(x,pred_bin_edges),curr_act_pred_avg,'uni',false);
            
            % (get activity binned by predicted)
            pred_use_trials = cellfun(@(act,act_pred,trial_bins) ...
                squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
                squeeze(~any(isnan(act_pred(:,curr_event_t)),2)) & ...
                ~isnan(trial_bins), ...
                curr_act,curr_act_pred,pred_trial_bins,'uni',false);
            
            act_predbinmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                accumarray(bins(use_trials & trial_cond(:,condition)), ...
                act(use_trials & trial_cond(:,condition)), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                curr_act_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
                permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
            
            act_pred_predbinmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                accumarray(bins(use_trials & trial_cond(:,condition)), ...
                act(use_trials & trial_cond(:,condition)), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                curr_act_pred_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
                permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
            
            % (get "fixed" predicted activity binned by predicted)
            act_predfix_predbinmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                accumarray(bins(use_trials & trial_cond(:,condition)), ...
                act(use_trials & trial_cond(:,condition)), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                curr_act_predfix_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
                permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
            
            % Get condition difference significance from shuffle
            n_shuff = 1000;
            trial_conditions_shuff = cellfun(@(x) AP_shake(x,1), ...
                repmat(trial_conditions_exp,1,n_shuff),'uni',false);
            
            act_predbinmean_condshuff = cell2mat(arrayfun(@(curr_shuff) ...
                cell2mat(arrayfun(@(condition) ...
                cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                accumarray(bins(use_trials & trial_cond(:,condition)), ...
                act(use_trials & trial_cond(:,condition)), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                curr_act_avg,pred_trial_bins,trial_conditions_shuff(:,curr_shuff), ...
                pred_use_trials,'uni',false)'), ...
                permute(1:size(trial_conditions,2),[1,3,2]),'uni',false)), ...
                permute(1:n_shuff,[1,3,4,2]),'uni',false));
            
            act_predfix_predbinmean_condshuff = cell2mat(arrayfun(@(curr_shuff) ...
                cell2mat(arrayfun(@(condition) ...
                cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
                accumarray(bins(use_trials & trial_cond(:,condition)), ...
                act(use_trials & trial_cond(:,condition)), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
                curr_act_predfix_avg,pred_trial_bins,trial_conditions_shuff(:,curr_shuff), ...
                pred_use_trials,'uni',false)'), ...
                permute(1:size(trial_conditions,2),[1,3,2]),'uni',false)), ...
                permute(1:n_shuff,[1,3,4,2]),'uni',false));
            
            % (measured data: null = no difference between conditions)
            cond_combos = nchoosek(1:size(trial_conditions,2),2);
            cond_sig_diff = false(n_act_bins,size(cond_combos,1));
            predfix_cond_sig_diff = false(n_act_bins,size(cond_combos,1));
            for curr_cond_combo = 1:size(cond_combos,1)
                curr_combo_diff = nanmean(act_predbinmean(:,:,cond_combos(curr_cond_combo,1)) - ...
                    act_predbinmean(:,:,cond_combos(curr_cond_combo,2)),2);
                shuff_prctile = squeeze(prctile(nanmean(act_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,1),:) - ...
                    act_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,2),:),2),[2.5,97.5],4));
                cond_sig_diff(:,curr_cond_combo) = ...
                    curr_combo_diff < shuff_prctile(:,1) | ...
                    curr_combo_diff > shuff_prctile(:,2);
                
                predfix_curr_combo_diff = nanmean(act_predfix_predbinmean(:,:,cond_combos(curr_cond_combo,1)) - ...
                    act_predfix_predbinmean(:,:,cond_combos(curr_cond_combo,2)),2);
                predfix_shuff_prctile = squeeze(prctile(nanmean(act_predfix_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,1),:) - ...
                    act_predfix_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,2),:),2),[2.5,97.5],4));
                predfix_cond_sig_diff(:,curr_cond_combo) = ...
                    predfix_curr_combo_diff < predfix_shuff_prctile(:,1) | ...
                    predfix_curr_combo_diff > predfix_shuff_prctile(:,2);
            end
            
            % Plot measured v predicted in bins
            figure(measured_v_pred_fig)
            
            % (measured vs binned predicted)
            subplot(length(plot_areas),length(timeavg_labels), ...
                sub2ind(fliplr([length(plot_areas),length(timeavg_labels)]),curr_timeavg,curr_area_idx));
            hold on;
            set(gca,'ColorOrder',timeavg_condition_colors{curr_timeavg});
            
            fill_cols = min(timeavg_condition_colors{curr_timeavg} + 0.5,1);
            for curr_cond = 1:size(trial_conditions,2)
                AP_errorfill( ...
                    squeeze(nanmean(act_pred_predbinmean(:,:,curr_cond),2)), ...
                    squeeze(nanmean(act_predfix_predbinmean(:,:,curr_cond),2)), ...
                    squeeze(AP_sem(act_predfix_predbinmean(:,:,curr_cond),2)), ...
                    fill_cols(curr_cond,:),1,false);
            end
            
            errorbar( ...
                squeeze(nanmean(act_pred_predbinmean,2)), ...
                squeeze(nanmean(act_predbinmean,2)), ...
                squeeze(AP_sem(act_predbinmean,2)), ...
                'linestyle','none','linewidth',2,'CapSize',10);
            
            n_col_bins = n_act_bins + 2;
            ctx_col = brewermap(n_col_bins,'*Greens');
            scatter( ...
                reshape(squeeze(nanmean(act_pred_predbinmean,2)),[],1), ...
                reshape(squeeze(nanmean(act_predbinmean,2)),[],1), ...
                60,repmat(ctx_col(1:n_act_bins,:),size(trial_conditions,2),1),'filled');
            
            xlabel(['Predicted (' num2str(plot_area) ')']);
            ylabel(['Measured (' num2str(plot_area) ')'])
            axis square;
            title(act_title);
            
            % (plot significant measured condition differences)
            % (* and o = significant in both measured and "fixed" predicted)
            curr_ylim = max(ylim);
            for curr_cond_combo = 1:size(cond_combos,1)
                % (plot * for measured condition differences)
                sig_y = curr_ylim + 0.1*curr_cond_combo;
                sig_x = pred_bin_centers;
                if any(cond_sig_diff(:,curr_cond_combo))
                    plot(sig_x(cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                        '*','MarkerSize',10,'color', ...
                        timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,1),:));
                    plot(sig_x(cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                        '*','MarkerSize',5,'color', ...
                        timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,2),:));
                end
                % (plot o for [predicted condition differences)
                if any(predfix_cond_sig_diff(:,curr_cond_combo))
                    plot(sig_x(predfix_cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                        'o','MarkerSize',15,'color', ...
                        timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,1),:));
                    plot(sig_x(predfix_cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                        'o','MarkerSize',10,'color', ...
                        timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,2),:));
                end
            end
            
            drawnow;
            
        end
    end
    
    %%%%%%%%%% TESTING %%%%%%%%%%%%%
    
    predfit = [];
    for i = 1:size(act_predbinmean,2)
        use_points = ~isnan(act_predbinmean(:,i,1)) & ~isnan(act_predbinmean(:,i,1));
        predfit(i,:) = polyfit(act_predbinmean(use_points,i,1),act_predbinmean(use_points,i,2),1);
    end
    
    % Store error in bins
    use_cond = 1;
    curr_offset = nanmean(act_predbinmean(:,:,1) - act_predbinmean(:,:,2),1);
    ctx_v_error{curr_data,1} = ...
        cellfun(@(x,cond,pred_use_trials) ...
        nanmean(x(pred_use_trials & cond(:,use_cond)),1), ...
        curr_act_pred_avg,trial_conditions_exp,pred_use_trials);
    ctx_v_error{curr_data,2} = ...
        curr_offset';
    
    linkaxes(get(measured_v_pred_fig,'Children'),'xy');
    
end

% Plot striatum/cortex offset relative to cortical activity
ctx_change = (ctx_v_error{2,1} - ctx_v_error{1,1})./(abs(ctx_v_error{2,1}) + abs(ctx_v_error{1,1}));
error_change = (ctx_v_error{2,2} - ctx_v_error{1,2})./(abs(ctx_v_error{2,2}) + abs(ctx_v_error{1,2}));
figure; 
plot(ctx_change,error_change,'.k','MarkerSize',20);
line([-1,1],[0,0],'color','k','linestyle','--');
line([0,0],[-1,1],'color','k','linestyle','--');
line([-1,1],[-1,1],'color','k');
xlabel('Cortex (pred striatum) change');
ylabel('Error change');

nonan_points = ~isnan(ctx_change) & ~isnan(error_change);
[r,p] = corrcoef(ctx_change(nonan_points),error_change(nonan_points));
disp(r(2));
disp(p(2));


a = (act_predbinmean(:,:,1) - act_predbinmean(:,:,2))./(abs(act_predbinmean(:,:,1)) + abs(act_predbinmean(:,:,2)));











