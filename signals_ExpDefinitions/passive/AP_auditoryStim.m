function AP_auditoryStim(t, events, parameters, visStim, inputs, outputs, audio)
% Present tones and white noise
% Pseudorandom (each stim presented once for each trial)
% Number of trials = number of repeats
%
% (invisible stim comes on with sound so that photodiode flip can be used
% as a sync)


%% Set up stimuli

% Duration and ITI
stim_time = 0.5;
min_iti = 2;
max_iti = 3;
step_iti = 0.1;

% Sounds
audio_freqs = linspace(8000,20000,4);
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

%% Set trial data

% Signals garbage: things can't happen at the exact same time as newTrial
new_trial_set = events.newTrial.delay(0);

% Start clock for trial
trial_t = t - t.at(new_trial_set);

% Set the stim order and ITIs for this trial
stimOrder = new_trial_set.map(@(x) randperm(length(audio_samples)));
stimITIs = new_trial_set.map(@(x) randsample(min_iti:step_iti:max_iti,length(audio_samples),true));

% Get the stim on times and the trial end time
trial_stimOn_times = stimITIs.map(@(x) [0,cumsum(x(1:end-1) + stim_time)]);
trial_end_time = stimITIs.map(@(x) sum(x) + stim_time*length(audio_samples));


%% Present stim

% Auditory
stim_num = trial_t.ge(trial_stimOn_times).sum.skipRepeats;
stim_id = map2(stimOrder,stim_num,@(stim_order,stim_num) stim_order(stim_num));
audio.Strix = stim_id.map(@(x) audio_samples{x});

% Visual
% (this is a quick and bad hack: have a 0% contrast come up on the screen
% so that there's something to synchronize the auditory with)
vis_stim = vis.grating(t, 'square', 'gaussian');
vis_stim.contrast = 0;
vis_stimOn = stim_id.to(stim_id.delay(stim_time));
vis_stim.show = vis_stimOn;
visStim.stim = stim;

endTrial = events.newTrial.setTrigger(trial_t.gt(trial_end_time));

%% Events

events.stimITIs = stimITIs;
events.stimOn = t.at(stim_id);
events.vis_stimOn = vis_stimOn; % hack to get photodiode flip with sound

events.stimFrequency = stim_id.map(@(x) audio_labels(x));

events.endTrial = endTrial;

end

















