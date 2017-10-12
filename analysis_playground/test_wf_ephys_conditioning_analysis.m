%% Loop through days, get licking piezo aligned to azimuths

animal = 'AP015';
expInfo_path = ['\\zserver.cortexlab.net\Data\expInfo\' animal];
expInfo_dir = dir(expInfo_path);
days = {expInfo_dir(find([expInfo_dir(3:end).isdir])+2).name};
days = days(end-5:end-1);

for curr_day = 1:length(days)
    day = days{curr_day};
    experiment =1;
    load_parts = struct;
    AP_load_experiment
    
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
    
end





