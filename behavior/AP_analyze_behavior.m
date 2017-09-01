%% Plot licking piezo aligned to water

% Get licks
lick_piezo_name = 'piezoLickDetector';
lick_pizeo_idx = strcmp({Timeline.hw.inputs.name}, lick_piezo_name);
licking_energy = abs(hilbert(Timeline.rawDAQData(:,lick_pizeo_idx)));

% Get reward delivery times
water_name = 'rewardEcho';
water_idx = strcmp({Timeline.hw.inputs.name}, water_name);
water_samples = find(Timeline.rawDAQData(1:end-1,water_idx) <= 2 & ...
    Timeline.rawDAQData(2:end,water_idx) > 2) + 1;

% How much licking to plot
surround_time = [-0.5,2];
surround_samples = surround_time/Timeline.hw.samplingInterval;

% Get lick aligned to water delivery
surround_samples = surround_time/Timeline.hw.samplingInterval;
water_surround_lick = bsxfun(@plus,water_samples,surround_samples(1):surround_samples(2));
water_aligned_lick = reshape(Timeline.rawDAQData(water_surround_lick,lick_pizeo_idx),size(water_surround_lick));
water_aligned_lick_energy = reshape(licking_energy(water_surround_lick),size(water_surround_lick));

% Plot
figure;
t = (surround_samples(1):surround_samples(2))*Timeline.hw.samplingInterval;

subplot(1,2,1);
imagesc(t,1:size(water_aligned_lick,1),water_aligned_lick)
line([0,0],ylim,'color','r');
caxis([-0.2 0.2])
colormap(colormap_BlueWhiteRed);
ylabel('Reward');
xlabel('Time from water delivery (s)');

subplot(1,2,2);
plot(t,nanmean(water_aligned_lick_energy,1),'k');
y = ylim; 
line([0,0],ylim,'color','r'); 
ylim(y);
xlabel('Time from water delivery (s)');
ylabel('Lick power');

%% Plot licking beam aligned to water

% Get licks
lick_piezo_name = 'beamLickDetector';
lick_pizeo_idx = strcmp({Timeline.hw.inputs.name}, lick_piezo_name);
lick_times = Timeline.rawDAQTimestamps(find(Timeline.rawDAQData(1:end-1,lick_pizeo_idx) <= 2 & ...
    Timeline.rawDAQData(2:end,lick_pizeo_idx) > 2) + 1);

% Get reward delivery times
water_name = 'rewardEcho';
water_idx = strcmp({Timeline.hw.inputs.name}, water_name);
water_samples = find(Timeline.rawDAQData(1:end-1,water_idx) <= 2 & ...
    Timeline.rawDAQData(2:end,water_idx) > 2) + 1;

% How much licking to plot
surround_time = [-1,10];
surround_samples = surround_time/Timeline.hw.samplingInterval;

% Get lick aligned to water delivery
surround_samples = surround_time/Timeline.hw.samplingInterval;
water_surround_lick = bsxfun(@plus,water_samples,surround_samples(1):surround_samples(2));
water_surround_lick(water_surround_lick < 1 | water_surround_lick > size(Timeline.rawDAQData,1)) = 1;
water_aligned_lick = reshape(Timeline.rawDAQData(water_surround_lick,lick_pizeo_idx),size(water_surround_lick)) < 2.5;

% Plot
figure;
t = (surround_samples(1):surround_samples(2))*Timeline.hw.samplingInterval;

imagesc(t,1:size(water_aligned_lick,1),water_aligned_lick)
line([0,0],ylim,'color','r');
colormap(gray);
ylabel('Reward');
xlabel('Time from water delivery (s)');

%% Plot lever presses aligned to stim

% Get photodiode flips closest to stim presentations
photodiode_name = 'photoDiode';
photodiode_idx = strcmp({Timeline.hw.inputs.name}, photodiode_name);
photodiode_flip_samples = find((Timeline.rawDAQData(1:end-1,photodiode_idx) <= 2) ~= ...
    (Timeline.rawDAQData(2:end,photodiode_idx) <= 2)) + 1;

[~,closest_stimOn_photodiode] = ...
    arrayfun(@(x) min(abs(signals_events.stimOnTimes(x) - ...
    Timeline.rawDAQTimestamps(photodiode_flip_samples))), ...
    1:length(signals_events.stimOnTimes));
stimOn_samples = photodiode_flip_samples(closest_stimOn_photodiode);

% How much lever to plot
surround_time = [-4,4];
surround_samples = surround_time/Timeline.hw.samplingInterval;

% Get lever around stim presentations
lever_r_name = 'lever_r';
lever_r_idx = strcmp({Timeline.hw.inputs.name}, lever_r_name);

stim_surround_lever = bsxfun(@plus,stimOn_samples,surround_samples(1):surround_samples(2));
stim_surround_lever(any(stim_surround_lever <= 0,2) | ...
    any(stim_surround_lever > size(Timeline.rawDAQData,1),2),:) = [];

stim_aligned_lever = reshape(Timeline.rawDAQData(stim_surround_lever,lever_r_idx),size(stim_surround_lever)) > 2.5;

figure;
t = (surround_samples(1):surround_samples(2))*Timeline.hw.samplingInterval;

subplot(1,2,1);
imagesc(t,1:size(stim_aligned_lever,1),stim_aligned_lever);
colormap(gray);
line([0,0],ylim,'color','r');
ylabel('Stim presentation');
xlabel('Time from stim onset (s)');
title('Lever pressing');

subplot(1,2,2);
plot(t,nanmean(stim_aligned_lever,1));
line([0,0],ylim,'color','r');
xlabel('Time from stim onset (s)');
ylabel('Fraction of trials with press');


%% Plot licking (energy) aligned to all lever presses/water

licking_energy = abs(hilbert(Timeline.rawDAQData(:,lick_pizeo_idx)));

% Get lick aligned to lever press times
lick_piezo_name = 'piezoLickDetector';
lick_pizeo_idx = strcmp({Timeline.hw.inputs.name}, lick_piezo_name);
lever_surround_samples = bsxfun(@plus,lever_r_push_samples,surround_samples(1):surround_samples(2));
lever_aligned_lick = reshape(licking_energy(lever_surround_samples),size(lever_surround_samples));

% Get lick aligned to water delivery
surround_time = [-4,4];
surround_samples = surround_time/Timeline.hw.samplingInterval;
water_surround_lick = bsxfun(@plus,water_samples,surround_samples(1):surround_samples(2));
water_aligned_lick = reshape(licking_energy(water_surround_lick),size(water_surround_lick));

% Plot
t = (surround_samples(1):surround_samples(2))/Timeline.hw.samplingInterval;
t0 = find(surround_samples(1):surround_samples(2) == 0);
figure;
plot(t,nanmean(water_aligned_lick,1),'k');
line([t0,t0],ylim,'color','r');
xlabel('Time from water delivery');
ylabel('Lick energy');

%% Plot piezo licking aligned to stimuli by azimuths

% Get licks
lick_piezo_name = 'piezoLickDetector';
lick_pizeo_idx = strcmp({Timeline.hw.inputs.name}, lick_piezo_name);
licking_energy = abs(hilbert(Timeline.rawDAQData(:,lick_pizeo_idx)));

% Get photodiode flips closest to stim presentations
photodiode_name = 'photoDiode';
photodiode_idx = strcmp({Timeline.hw.inputs.name}, photodiode_name);
photodiode_flip_samples = find(abs(Timeline.rawDAQData(1:end-1,photodiode_idx) - ...
    Timeline.rawDAQData(2:end,photodiode_idx)) > 2) + 1;

[~,closest_stimOn_photodiode] = ...
    arrayfun(@(x) min(abs(signals_events.stimOnTimes(x) - ...
    Timeline.rawDAQTimestamps(photodiode_flip_samples))), ...
    1:length(signals_events.stimOnTimes));
stimOn_samples = photodiode_flip_samples(closest_stimOn_photodiode);

% Define surround time
surround_time = [-0.2,3];
surround_samples = surround_time/Timeline.hw.samplingInterval;

% Get licking aligned to each azimuth
[azimuths,~,azimuths_idx] = unique(signals_events.trialAzimuthValues);
azimuths_n = accumarray(azimuths_idx,1);

azimuth_aligned_lick = cellfun(@(x) nan(x, ...
    length(surround_samples(1):surround_samples(2))),num2cell(azimuths_n),'uni',false);
azimuth_aligned_lick_energy = cellfun(@(x) nan(x, ...
    length(surround_samples(1):surround_samples(2))),num2cell(azimuths_n),'uni',false);
for trialAzimuth_idx = 1:length(azimuths)
        
        curr_azimuth = azimuths(trialAzimuth_idx);       
        use_trials = signals_events.trialAzimuthValues == curr_azimuth;
        curr_stim_samples = stimOn_samples(use_trials);
        stim_surround_lick = bsxfun(@plus,curr_stim_samples,surround_samples(1):surround_samples(2));
        
        stim_surround_lick(any(stim_surround_lick <= 0,2) | ...
            any(stim_surround_lick > size(Timeline.rawDAQData,1),2),:) = [];
        
        azimuth_aligned_lick{trialAzimuth_idx} = reshape(Timeline.rawDAQData( ...
            stim_surround_lick,lick_pizeo_idx),size(stim_surround_lick));
        
        azimuth_aligned_lick_energy{trialAzimuth_idx} = reshape(licking_energy(stim_surround_lick), ...
            size(stim_surround_lick));
        
end

% Plot
figure;

t = (surround_samples(1):surround_samples(2))*Timeline.hw.samplingInterval;

subplot(1,2,1);
imagesc(t,1:size(vertcat(azimuth_aligned_lick{:}),1),vertcat(azimuth_aligned_lick{:}))
line([0,0],ylim,'color','k');
line([1,1],ylim,'color','k');
for i = 1:length(azimuths)
    line(xlim,repmat(sum(azimuths_n(1:i))+0.5,2,1),'color','k')
end
caxis([-0.2 0.2])
colormap(colormap_BlueWhiteRed);
ylabel('Stim presentation');
xlabel('Time from stim onset (s)')

subplot(1,2,2); hold on;
set(gca,'ColorOrder',copper(length(azimuths)));
plot(t,cell2mat(cellfun(@(x) nanmean(x,1),azimuth_aligned_lick_energy,'uni',false))');
line([0,0],ylim,'color','k');
line([1,1],ylim,'color','k');
xlabel('Time from stim onset (s)');
ylabel('Lick energy');
legend(cellfun(@(x) num2str(x),num2cell(azimuths),'uni',false));

%% Plot beam licking aligned to stimuli by condition

% Define conditions
%signals_conditions = signals_events.trialAzimuthValues;
signals_conditions = signals_events.trialOrientationValues;

% Get licks
lick_piezo_name = 'beamLickDetector';
lick_pizeo_idx = strcmp({Timeline.hw.inputs.name}, lick_piezo_name);

% Get reward delivery times
water_name = 'rewardEcho';
water_idx = strcmp({Timeline.hw.inputs.name}, water_name);
water_samples = find(Timeline.rawDAQData(1:end-1,water_idx) <= 2 & ...
    Timeline.rawDAQData(2:end,water_idx) > 2) + 1;

% Define surround time
surround_time = [-5,5];

% Get licking aligned to each condition
[conditions,~,conditions_idx] = unique(signals_conditions);
conditions_n = accumarray(conditions_idx,1);

condition_aligned_lick = cellfun(@(x) nan(x, ...
    ceil(diff(surround_time)/Timeline.hw.samplingInterval)),num2cell(conditions_n),'uni',false);
condition_aligned_water = cellfun(@(x) nan(x, ...
    ceil(diff(surround_time)/Timeline.hw.samplingInterval)),num2cell(conditions_n),'uni',false);
for trialCondition_idx = 1:length(conditions)
        
        curr_condition = conditions(trialCondition_idx);       
        use_trials = signals_conditions == curr_condition;
      
        curr_stim_times = stimOn_times(use_trials);
        stim_surround_lick = bsxfun(@plus,curr_stim_times, ...
            surround_time(1):Timeline.hw.samplingInterval:surround_time(2));
        
        stim_surround_lick(any(stim_surround_lick <= min(Timeline.rawDAQTimestamps),2) | ...
            any(stim_surround_lick > max(Timeline.rawDAQTimestamps),2),:) = [];
        
        condition_aligned_lick{trialCondition_idx} = interp1(Timeline.rawDAQTimestamps, ...
            Timeline.rawDAQData(:,lick_pizeo_idx),stim_surround_lick) < 2.5;
        condition_aligned_water{trialCondition_idx} = interp1(Timeline.rawDAQTimestamps, ...
            Timeline.rawDAQData(:,water_idx),stim_surround_lick) < 2.5;
                
end

% Plot
figure;
t = surround_time(1):Timeline.hw.samplingInterval:surround_time(2);
stim_time = median(stimOff_times - stimOn_times);

% Licks
subplot(1,2,1);
imagesc(t,1:size(vertcat(condition_aligned_lick{:}),1),vertcat(condition_aligned_lick{:}))
ylim([0 size(vertcat(condition_aligned_lick{:}),1)]);
xlim(surround_time)
line([0,0],ylim,'color','r');
line([stim_time,stim_time],ylim,'color','r');
for i = 1:length(conditions)
    line(xlim,repmat(sum(conditions_n(1:i))+0.5,2,1),'color','r')
end
colormap(gray);
ylabel('Stim presentation');
xlabel('Time from stim onset (s)')
title('Licks')

% Reward
subplot(1,2,2);
imagesc(t,1:size(vertcat(condition_aligned_water{:}),1),vertcat(condition_aligned_water{:}))
ylim([0 size(vertcat(condition_aligned_water{:}),1)]);
xlim(surround_time)
line([0,0],ylim,'color','r');
line([stim_time,stim_time],ylim,'color','r');
for i = 1:length(conditions)
    line(xlim,repmat(sum(conditions_n(1:i))+0.5,2,1),'color','r')
end
colormap(gray);
ylabel('Stim presentation');
xlabel('Time from stim onset (s)')
title('Water')



%% Choiceworld: wheel movement around stim

surround_time = [-0.5,2];
surround_samples = surround_time/Timeline.hw.samplingInterval;

% Get wheel aligned to stim onset
rotaryEncoder_idx = strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder');
surround_time = surround_time(1):Timeline.hw.samplingInterval:surround_time(2);
pull_times = bsxfun(@plus,signals_events.stimOnTimes',surround_time);

stim_aligned_wheel_raw = interp1(Timeline.rawDAQTimestamps, ...
    Timeline.rawDAQData(:,rotaryEncoder_idx),pull_times);
stim_aligned_wheel = bsxfun(@minus,stim_aligned_wheel_raw, ...
    nanmedian(stim_aligned_wheel_raw(:,surround_time < 0),2));

figure; hold on;

use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues > 0 & signals_events.hitValues;
plot(surround_time,nanmedian(stim_aligned_wheel(use_trials,:),1),'k','linewidth',2)

use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues > 0 & ~signals_events.hitValues;
plot(surround_time,nanmedian(stim_aligned_wheel(use_trials,:),1),'r','linewidth',2)

use_trials = signals_events.trialSideValues == -1 & signals_events.trialContrastValues > 0 & signals_events.hitValues;
plot(surround_time,nanmedian(stim_aligned_wheel(use_trials,:),1),'b','linewidth',2)

use_trials = signals_events.trialSideValues == -1 & signals_events.trialContrastValues > 0 & ~signals_events.hitValues;
plot(surround_time,nanmedian(stim_aligned_wheel(use_trials,:),1),'m','linewidth',2)

legend({'Stim right hit','Stim right miss','Stim left hit','Stim left miss'})
xlabel('Time from stim');
ylabel('Wheel displacement');

% Define time to first wheel movement
thresh_displacement = 2;
[~,wheel_move_sample] = max(abs(stim_aligned_wheel) > thresh_displacement,[],2);
wheel_move_time = arrayfun(@(x) pull_times(x,wheel_move_sample(x)),1:size(pull_times,1));
wheel_move_time(wheel_move_sample == 1) = NaN;

%% Choiceworld: wheel movement around wheel movement start (sanity check)

% Get wheel movement time for each trial
surround_time = [-0.5,2];
surround_samples = surround_time/Timeline.hw.samplingInterval;

rotaryEncoder_idx = strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder');
surround_time = surround_time(1):Timeline.hw.samplingInterval:surround_time(2);
pull_times = bsxfun(@plus,signals_events.stimOnTimes',surround_time);

stim_aligned_wheel_raw = interp1(Timeline.rawDAQTimestamps, ...
    Timeline.rawDAQData(:,rotaryEncoder_idx),pull_times);
stim_aligned_wheel = bsxfun(@minus,stim_aligned_wheel_raw, ...
    nanmedian(stim_aligned_wheel_raw(:,surround_time < 0),2));

thresh_displacement = 10;
[~,wheel_move_sample] = max(abs(stim_aligned_wheel) > thresh_displacement,[],2);
wheel_move_time = arrayfun(@(x) pull_times(x,wheel_move_sample(x)),1:size(pull_times,1));
wheel_move_time(wheel_move_sample == 1) = NaN;

% Get wheel aligned to stim onset
surround_time = [-0.5,2];
surround_samples = surround_time/Timeline.hw.samplingInterval;

rotaryEncoder_idx = strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder');
surround_time = surround_time(1):Timeline.hw.samplingInterval:surround_time(2);
pull_times = bsxfun(@plus,wheel_move_time',surround_time);

stim_aligned_wheel_raw = interp1(Timeline.rawDAQTimestamps, ...
    Timeline.rawDAQData(:,rotaryEncoder_idx),pull_times);
stim_aligned_wheel = bsxfun(@minus,stim_aligned_wheel_raw, ...
    nanmedian(stim_aligned_wheel_raw(:,surround_time < 0),2));

figure; hold on;

use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues > 0 & signals_events.hitValues & ~signals_events.repeatTrialValues;
plot(surround_time,nanmedian(stim_aligned_wheel(use_trials,:),1),'k','linewidth',2)

use_trials = signals_events.trialSideValues == 1 & signals_events.trialContrastValues > 0 & ~signals_events.hitValues & ~signals_events.repeatTrialValues;
plot(surround_time,nanmedian(stim_aligned_wheel(use_trials,:),1),'r','linewidth',2)

use_trials = signals_events.trialSideValues == -1 & signals_events.trialContrastValues > 0 & signals_events.hitValues & ~signals_events.repeatTrialValues;
plot(surround_time,nanmedian(stim_aligned_wheel(use_trials,:),1),'b','linewidth',2)

use_trials = signals_events.trialSideValues == -1 & signals_events.trialContrastValues > 0 & ~signals_events.hitValues & ~signals_events.repeatTrialValues;
plot(surround_time,nanmedian(stim_aligned_wheel(use_trials,:),1),'m','linewidth',2)

legend({'Stim right hit','Stim right miss','Stim left hit','Stim left miss'})
xlabel('Time from stim');
ylabel('Wheel displacement');

% Define time to first wheel movement
thresh_displacement = 2;
[~,wheel_move_sample] = max(abs(stim_aligned_wheel) > thresh_displacement,[],2);
wheel_move_time = arrayfun(@(x) pull_times(x,wheel_move_sample(x)),1:size(pull_times,1));
wheel_move_time(wheel_move_sample == 1) = NaN;


%% Get licking bout statistics (from beam lick)
% ALSO KEEP IN MIND THIS PROBABLY CHANGES! mice start licking more
% specifically to stimuli at the end maybe?

lick_bout_cutoff = 1; % seconds between what's considered a lick bout
lick_bout_times = signals_events.n_licksTimes(find(diff( ...
    signals_events.n_licksTimes) >= lick_bout_cutoff)+1);

min_lick_to_reward_time = 0.01; % time from lick to reward to define rewarded lick
lick_to_reward_time = min(abs(bsxfun(@minus,lick_bout_times',signals_events.hitTimes)),[],2);
rewarded_licks = lick_to_reward_time < min_lick_to_reward_time;

% I don't think this makes sense... but get some kind of statistics for
% probability of licking during stim vs. not
% use_stim = min(length(signals_events.stimOnTimes),length(signals_events.stimOffTimes));
% epoch_times = reshape([signals_events.stimOnTimes(1:use_stim);signals_events.stimOffTimes(1:use_stim)],[],1);
% epoch_licks = histcounts(lick_bout_times, epoch_times);
% 
% epoch_lick_rate = epoch_licks'./diff(epoch_times);

















