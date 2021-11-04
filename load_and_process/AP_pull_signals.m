% AP_pull_signals
%
% Pull out the relevant signals for specific expDef
% (not function - to be used in context of variables from
% AP_load_experiment)

[~,expDef] = fileparts(block.expDef);
switch expDef
    case {'vanillaChoiceworld','vanillaChoiceworldBias', ...
            'vanillaChoiceworldNoRepeats','vanillaChoiceworldFastwheel', ...
            'vanillaChoiceworldNoCue','vanillaChoiceworldBiasNoCue'}
        % Hit/miss recorded for previous trial, circshift to align
        signals_events.hitValues = circshift(signals_events.hitValues,[0,-1]);
        signals_events.missValues = circshift(signals_events.missValues,[0,-1]);
        
        % Get number of completed trials (if uncompleted last trial)
        n_trials = length(signals_events.endTrialTimes);
        
        % Get stim on times by closest photodiode flip
        [~,closest_stimOn_photodiode] = ...
            arrayfun(@(x) min(abs(signals_events.stimOnTimes(x) - ...
            photodiode_flip_times)), ...
            1:n_trials);
        stimOn_times = photodiode_flip_times(closest_stimOn_photodiode);
        
        % Check that the stim times aren't off by a certain threshold
        % (skip the first one - that's usually delayed a little)
        stim_time_offset_thresh = 0.05;
        if any(abs(stimOn_times(2:end) - signals_events.stimOnTimes(2:n_trials)') >= ...
                stim_time_offset_thresh)
            figure;
            plot(stimOn_times - signals_events.stimOnTimes(1:n_trials)','.k')
            line(xlim,repmat(stim_time_offset_thresh,2,1),'color','r');
            line(xlim,repmat(-stim_time_offset_thresh,2,1),'color','r');
            warning('Stim signals/photodiode offset over threshold');
            xlabel('Stim number');
            ylabel('Photodiode - signals stim time');
            title([animal ' ' day ' ' num2str(experiment)]);
        end
        
        % Get first movement time after stim onset
        surround_time = [-0.5,2];
        surround_sample_rate = 1/Timeline.hw.samplingInterval; % (match this to framerate)
        surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);
        pull_times = bsxfun(@plus,stimOn_times,surround_time_points);
        
        stim_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,pull_times);
        
        % (set a threshold in speed and time for wheel movement)
        thresh_displacement = 0.025;
        time_over_thresh = 0.05; % ms over velocity threshold to count
        samples_over_thresh = time_over_thresh.*surround_sample_rate;
        wheel_over_thresh_fullconv = convn( ...
            abs(stim_aligned_wheel) > thresh_displacement, ...
            ones(1,samples_over_thresh)) >= samples_over_thresh;
        wheel_over_thresh = wheel_over_thresh_fullconv(:,end-size(stim_aligned_wheel,2)+1:end);
        
        [move_trial,wheel_move_sample] = max(wheel_over_thresh,[],2);
        wheel_move_time = arrayfun(@(x) pull_times(x,wheel_move_sample(x)),1:size(pull_times,1))';
        wheel_move_time(~move_trial) = NaN;
        
        % Get conditions for all trials
        
        % (trial_timing)
        stim_to_move = padarray(wheel_move_time - stimOn_times,[n_trials-length(stimOn_times),0],NaN,'post');
        stim_to_feedback = signals_events.responseTimes(1:n_trials)' - stimOn_times(1:n_trials);
        
        % (early vs late move)
        trial_timing = 1 + (stim_to_move > 0.5);
        
        % (choice and outcome)
        go_left = (signals_events.trialSideValues == 1 & signals_events.hitValues == 1) | ...
            (signals_events.trialSideValues == -1 & signals_events.missValues == 1);
        go_right = (signals_events.trialSideValues == -1 & signals_events.hitValues == 1) | ...
            (signals_events.trialSideValues == 1 & signals_events.missValues == 1);
        trial_choice = go_right(1:n_trials)' - go_left(1:n_trials)';
        trial_outcome = signals_events.hitValues(1:n_trials)'-signals_events.missValues(1:n_trials)';
        
        % (trial conditions: [contrast,side,choice,timing])
        contrasts = [0,0.06,0.125,0.25,0.5,1];
        sides = [-1,1];
        choices = [-1,1];
        timings = [1,2];
        
        conditions = combvec(contrasts,sides,choices,timings)';
        n_conditions = size(conditions,1);
        
        trial_conditions = ...
            [signals_events.trialContrastValues(1:n_trials)', signals_events.trialSideValues(1:n_trials)', ...
            trial_choice(1:n_trials), trial_timing(1:n_trials)];
        [~,trial_id] = ismember(trial_conditions,conditions,'rows');
        
    case {'AP_stimWheelRight','AP_stimWheelLeft'}
        % Hit/miss recorded for previous trial, circshift to align
        signals_events.hitValues = circshift(signals_events.hitValues,[0,-1]);
        signals_events.missValues = circshift(signals_events.missValues,[0,-1]);
        
        % Get number of completed trials (if uncompleted last trial)
        n_trials = length(signals_events.endTrialTimes);
        
        % Get stim on times by closest photodiode flip
        [~,closest_stimOn_photodiode] = ...
            arrayfun(@(x) min(abs(signals_events.stimOnTimes(x) - ...
            photodiode_flip_times)), ...
            1:n_trials);
        stimOn_times = photodiode_flip_times(closest_stimOn_photodiode);
        
        % Check that the stim times aren't off by a certain threshold
        % (skip the first one - that's usually delayed a little)
        stim_time_offset_thresh = 0.05;
        if any(abs(stimOn_times(2:end) - signals_events.stimOnTimes(2:n_trials)') >= ...
                stim_time_offset_thresh)
            figure;
            plot(stimOn_times - signals_events.stimOnTimes(1:n_trials)','.k')
            line(xlim,repmat(stim_time_offset_thresh,2,1),'color','r');
            line(xlim,repmat(-stim_time_offset_thresh,2,1),'color','r');
            warning('Stim signals/photodiode offset over threshold');
            xlabel('Stim number');
            ylabel('Photodiode - signals stim time');
            title([animal ' ' day ' ' num2str(experiment)]);
        end
                
        % Get wheel movement on/offsets
        wheel_starts = Timeline.rawDAQTimestamps(diff([0;wheel_move]) == 1)';
        wheel_stops = Timeline.rawDAQTimestamps(diff([wheel_move;0]) == -1)';
        
        % (stim move: first move after stim)
        % (give this a little leeway, sometimes movement starts early but
        % stim comes on anyway)
        stim_leeway = 0.1;
        wheel_move_stim_idx = ...
            arrayfun(@(stim) find(wheel_starts > stim-stim_leeway,1,'first'), ...
            stimOn_times);
        
        % (response move: last move start before response signal)
        wheel_move_response_idx = ...
            arrayfun(@(response) find(wheel_starts <= response,1,'last'), ...
            signals_events.responseTimes(1:n_trials)');
        
        % (iti move: move start with no stim on screen)
        stimOff_times = signals_events.stimOffTimes';
        stimOn_epochs = logical(interp1([0;stimOn_times;stimOff_times], ...
            [0;ones(size(stimOn_times));zeros(size(stimOff_times))], ...
            Timeline.rawDAQTimestamps','previous','extrap'));
        wheel_move_iti_idx = find(ismember(wheel_starts, ...
            Timeline.rawDAQTimestamps(~stimOn_epochs)));
        
        % Get time from stim to rewarded movement onset and feedback
        stim_to_move = wheel_starts(wheel_move_stim_idx) - stimOn_times(1:n_trials);
        stim_to_feedback = signals_events.responseTimes(1:n_trials)' - stimOn_times(1:n_trials);

        % (choice and outcome)
        go_left = (signals_events.trialSideValues == 1 & signals_events.hitValues == 1) | ...
            (signals_events.trialSideValues == -1 & signals_events.missValues == 1);
        go_right = (signals_events.trialSideValues == -1 & signals_events.hitValues == 1) | ...
            (signals_events.trialSideValues == 1 & signals_events.missValues == 1);
        trial_choice = go_right(1:n_trials)' - go_left(1:n_trials)';
        trial_outcome = signals_events.hitValues(1:n_trials)'-signals_events.missValues(1:n_trials)';
        
    case {'AP_sparseNoise'}
        % Don't do anything: stim info is pulled out in
        % lilrig_retinotopy
        
    case {'AP_visualAuditoryPassive','AP_visualAuditoryPairing','AP_visualAuditoryPairingHalf'}
        % Get stim times (first flip is initializing gray to black)
        stimOn_times = photodiode_flip_times(2:2:end);
        vis_azimuth = signals_events.visAzimuthValues;
        aud_freq = signals_events.auditoryFrequencyValues;
        if isfield(signals_events,'visContrastValues')
            vis_contrast = signals_events.visContrastValues;
        else
            vis_contrast = ones(size(vis_azimuth));
        end
        
        % (temporary - set stim IDs)
        aud_freq_nonan = aud_freq;
        aud_freq_nonan(isnan(aud_freq)) = -1;
        trial_conditions = [vis_azimuth;aud_freq_nonan;vis_contrast];
        conds = unique(trial_conditions','rows');
        [~,stimIDs] = ismember(trial_conditions',conds,'rows');
        
    case 'AP_choiceWorldStimPassive'
        % This is kind of a dumb hack to get the stimOn times, maybe not
        % permanent unless it works fine: get stim times by checking for
        % close to the median photodiode flip difference
        block_stim_iti = mean(diff(block.stimWindowUpdateTimes));
        
        photodiode_flip_diff = diff(stimScreen_on_t(photodiode_flip));
        median_photodiode_flip_diff = mode(round(photodiode_flip_diff*10)/10);
        
        stimOn_idx = find(abs(photodiode_flip_diff-median_photodiode_flip_diff) < 0.1);
        
        stimOn_times = stimScreen_on_t(photodiode_flip(stimOn_idx))';
        
        % Set stimID as the contrast*side
        % (use last n values - sometimes short buffer times means some
        % stimuli in the beginning could be missed)
        use_signals_stim = size(signals_events.visualParamsValues,2)-length(stimOn_times)+1: ...
            size(signals_events.visualParamsValues,2);
        stimIDs = sign(signals_events.visualParamsValues(1,use_signals_stim))'.* ...
            signals_events.visualParamsValues(2,use_signals_stim)';
        
    case {'AP_lcrGratingPassive','AP_contrastGratingPassiveRight'}
        % Get stim times (first flip is initializing gray to black)
        stimOn_times = photodiode_flip_times(2:2:end);
        
        % Check number of stim matches photodiode
        % (once I saw some weird extra flip at the end of an
        % experiment? so added case to only use first n flips)
        if length(signals_events.stimAzimuthValues) ~= length(stimOn_times)
            warning([animal ' ' day ': different stim number signals and photodiode']);
            stimOn_times = stimOn_times(1:length(signals_events.stimAzimuthValues));
        end
        
        % Get stim ID and conditions
        contrasts = unique(signals_events.stimContrastValues);
        azimuths = unique(signals_events.stimAzimuthValues);
        
        conditions = combvec(contrasts,azimuths)';
        n_conditions = size(conditions,1);
        
        trial_conditions = ...
            [signals_events.stimContrastValues; signals_events.stimAzimuthValues]';
        [~,stimIDs] = ismember(trial_conditions,conditions,'rows');
        
    case 'AP_lcrGratingPassiveFlicker'
        % Flickering stim: get first photodiode after long gap
        % (ignore the first flip because that's initializing)
        iti_min = 0.5;
        stimOn_times = photodiode_flip_times([find(diff(photodiode_flip_times) > iti_min)+1]);
        
        % Check number of stim matches photodiode
        if length(signals_events.stimAzimuthValues) ~= length(stimOn_times)
            error('Different stim number signals and photodiode')
        end
        
        % Get stim ID and conditions
        contrasts = unique(signals_events.stimContrastValues);
        azimuths = unique(signals_events.stimAzimuthValues);
        
        conditions = combvec(contrasts,azimuths)';
        n_conditions = size(conditions,1);
        
        trial_conditions = ...
            [signals_events.stimContrastValues; signals_events.stimAzimuthValues]';
        [~,stimIDs] = ismember(trial_conditions,conditions,'rows');
        
    case 'AP_localize_choiceWorldStimPassive'
        % get stim times - first stim photodiode is messed up so throw it out
        stimOn_times = photodiode_flip_times(2:2:end);
        
        % sanity check: times between stim on times in signals
        signals_photodiode_iti_diff = diff(signals_events.stimOnTimes(2:end)) - diff(stimOn_times)';
        if any(signals_photodiode_iti_diff > 0.1)
            error('mismatching signals/photodiode stim ITIs')
        end
        
        % Get stim ID and conditions
        azimuths = unique(signals_events.stimAzimuthValues);
        altitudes = unique(signals_events.stimAltitudeValues);
        
        trial_conditions = reshape(signals_events.visualParamsValues,2,[])';
        
        conditions = unique(trial_conditions,'rows');
        n_conditions = size(conditions,1);
        
        [~,stimIDs] = ismember(trial_conditions,conditions,'rows');
        
        % Get rid of the first one for now
        trial_conditions = trial_conditions(2:end);
        stimIDs = stimIDs(2:end);
        
    case 'AP_auditoryStim'
        % Auditory stim only, use audioOut to get times
        speaker_idx = strcmp({Timeline.hw.inputs.name}, 'audioOut');
        speaker_threshold = 0.04; % eyeballed this
        speaker_flip_times = ...
            Timeline.rawDAQTimestamps( ...
            find(abs(Timeline.rawDAQData(1:end-1,speaker_idx)) < speaker_threshold & ...
            abs(Timeline.rawDAQData(2:end,speaker_idx)) > speaker_threshold)+1);
        
        iti_min = 1.9;
        stimOn_times = speaker_flip_times([1,find(diff(speaker_flip_times) > iti_min)+1]);
        
        % TEMPORARY: use first and last to interpolate Signals
        first_last_stim_tl = stimOn_times([1,end]);
        first_last_stim_block = block.events.stimOnTimes([1,end]);
        
        block_fieldnames = fieldnames(block.events);
        block_values_idx = cellfun(@(x) ~isempty(x),strfind(block_fieldnames,'Values'));
        block_times_idx = cellfun(@(x) ~isempty(x),strfind(block_fieldnames,'Times'));
        for curr_times = find(block_times_idx)'
            if isempty(signals_events.(block_fieldnames{curr_times}))
                % skip if empty
                continue
            end
            signals_events.(block_fieldnames{curr_times}) = ...
                interp1(first_last_stim_block,first_last_stim_tl, ...
                block.events.(block_fieldnames{curr_times}),'linear','extrap');
        end
        stimOn_times = signals_events.stimOnTimes;
        stimIDs = signals_events.stimFrequencyValues;
        
    otherwise
        warning(['Signals protocol with no analysis script:' expDef]);
end


