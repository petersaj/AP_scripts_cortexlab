function AP_visualAuditoryPairing(t, events, parameters, visStim, inputs, outputs, audio)

% Pair one visual stimulus with noise

%% Stimuli

% Duration and ITI
stim_time = 0.5;
min_iti = 2;
max_iti = 3;
step_iti = 0.1;

% Visual stimuli
azimuths = [35,90];
sigma = [20,20];
contrast = 1;
spatialFreq = 1/15;

% Auditory stim (one for each azimuth)
noise_volume = 0.05;

audDev = audio.Devices('Strix');
audioSampleRate = audDev.DefaultSampleRate;
audioChannels = audDev.NrOutputChannels;

noise_samples = {randn(2,audioSampleRate*stim_time)*noise_volume};
auditory_samples = [cell(1,1),noise_samples];


%% Set trial data

% Signals garbage: things can't happen at the exact same time as newTrial
new_trial_set = events.newTrial.delay(0);

% Start clock for trial
trial_t = t - t.at(new_trial_set);

% Set the stim order and ITIs for this trial
stimOrder = new_trial_set.map(@(x) randperm(length(azimuths));
stimITIs = new_trial_set.map(@(x) randsample(min_iti:step_iti:max_iti,length(azimuths),true));

% Get the stim on times and the trial end time
trial_stimOn_times = stimITIs.map(@(x) [0,cumsum(x(1:end-1) + stim_time)]);
trial_end_time = stimITIs.map(@(x) sum(x) + stim_time*length(azimuths));

%% Present stim

% Get which condition should be presented
condition_num = trial_t.ge(trial_stimOn_times).sum.skipRepeats;
condition_id = map2(stimOrder,condition_num,@(stim_order,stim_num) stim_order(stim_num));

% Visual
vis_azimuth = stim_id.map(@(x) azimuths(x));

vis_stim = vis.grating(t, 'square', 'gaussian');
vis_stim.azimuth = vis_azimuth;
vis_stim.contrast = 1;
vis_stimOn = condition_id.to(condition_id.delay(stim_time));
vis_stim.show = vis_stimOn;
visStim.stim = stim;

% Auditory
audio.Strix = condition_id.map(@(x) auditory_samples{x});

% End trial after all stim are presented
endTrial = events.newTrial.setTrigger(trial_t.gt(trial_end_time));






















