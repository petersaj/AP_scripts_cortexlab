function AP_visualAuditoryPassive(t, events, parameters, visStim, inputs, outputs, audio)

% Present visual and auditory stimuli separately and randomly within trial

%% Set up stimuli

% Duration and ITI
stim_time = 0.5;
min_iti = 2;
max_iti = 3;
step_iti = 0.1;

% Visual stim
sigma = [20,20];
azimuths = [-90,0,90];
contrast = 1;
spatialFreq = 1/15;

% Auditory stim
audio_freqs = [8000,14000,20000];
audioRampDuration = 0.01;
tone_volume = 0.05;
noise_volume = 0.05;

audDev = audio.Devices('Strix');
audioSampleRate = audDev.DefaultSampleRate;
audioChannels = audDev.NrOutputChannels;

tone_samples = arrayfun(@(freq) ...
    aud.pureTone(freq,stim_time,audioSampleRate, ...
    audioRampDuration,audioChannels)*tone_volume,audio_freqs,'uni',false);
noise_samples = {randn(2,audioSampleRate*stim_time)*noise_volume};

audio_samples = [tone_samples,noise_samples];
audio_labels = [audio_freqs,0];

% Concatenate conditions
n_conditions = length(azimuths) + length(audio_samples);
condition_vis_azimuth = [azimuths,zeros(1,length(audio_samples))];
condition_vis_contrast = [ones(1,length(azimuths)),zeros(1,length(audio_samples))];
condition_auditory_samples = [cell(1,length(azimuths)),audio_samples];
condition_auditory_labels = [NaN(1,length(azimuths)),audio_labels];

%% Set trial data

% Signals garbage: things can't happen at the exact same time as newTrial
new_trial_set = events.newTrial.delay(0);

% Start clock for trial
trial_t = t - t.at(new_trial_set);

% Set the stim order and ITIs for this trial
stimOrder = new_trial_set.map(@(x) randperm(n_conditions));
stimITIs = new_trial_set.map(@(x) randsample(min_iti:step_iti:max_iti,n_conditions,true));

% Get the stim on times and the trial end time
trial_stimOn_times = stimITIs.map(@(x) [0,cumsum(x(1:end-1) + stim_time)]);
trial_end_time = stimITIs.map(@(x) sum(x) + stim_time*n_conditions);

%% Present stim

% Get which condition should be presented
condition_num = trial_t.ge(trial_stimOn_times).sum.skipRepeats;
condition_id = map2(stimOrder,condition_num,@(stim_order,stim_num) stim_order(stim_num));

% Visual
vis_azimuth = stim_id.map(@(x) condition_vis_azimuth(x));
vis_contrast = stim_id.map(@(x) condition_vis_contrast(x));

vis_stim = vis.grating(t, 'square', 'gaussian');
vis_stim.azimuth = vis_azimuth;
vis_stim.contrast = vis_contrast;
vis_stimOn = condition_id.to(condition_id.delay(stim_time));
vis_stim.show = vis_stimOn;
visStim.stim = stim;

% Auditory
audio.Strix = condition_id.map(@(x) condition_auditory_samples{x});

% End trial after all stim are presented
endTrial = events.newTrial.setTrigger(trial_t.gt(trial_end_time));


%% Events

events.stimITIs = stimITIs;
events.stimOn = vis_stimOn;

events.visAzimuth = vis_azimuth;
events.visContrast = vis_contrast;
events.auditoryFrequency = condition_id.map(@(x) condition_auditory_labels(x));

events.endTrial = endTrial;





