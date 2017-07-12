%% Plot licking aligned to water

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

%% Plot licking aligned to stimuli by azimuths

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
 


