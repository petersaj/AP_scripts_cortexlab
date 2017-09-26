function AP_visAudioPassive(t, events, parameters, visStim, inputs, outputs, audio)
% AP 2017-09-22: present random tones and visual stimuli

%% Set up stimuli

stimTime = 0.5;
itiTimeRange= [1:0.1:2];

% Sound stim
audioFreqRange = linspace(5000,20000,5);
audioSampleRate = 192e3;

audioRampDuration = 0.01;
audioChannels = 2;

% Visual stim
vis_az_el_sig = {[-60,90,45,360],[0,90,45,360],[60,90,45,360]};
spatialFrequency = 0.01;
contrast = 1;
stimFlickerFrequency = 5;
orientation = 45; 


%% Stim times

n_audio_reps = 5;
n_visual_reps = 5;

% Audio
audio_stim_uniform = repmat(audioFreqRange,n_audio_reps,1);
audio_stim_shuffle = audio_stim_uniform(randperm(numel(audio_stim_uniform)));

audio_itiTimes = randsample(itiTimeRange,length(audio_stim_shuffle),true);
audio_startTimes = cumsum(audio_itiTimes + stimTime);

% Visual
visual_stim_uniform = repmat(vis_az_el_sig,n_visual_reps,1);
visual_stim_shuffle = visual_stim_uniform(randperm(numel(visual_stim_uniform)));

visual_itiTimes = randsample(itiTimeRange,length(visual_stim_shuffle),true);
visual_startTimes = cumsum(visual_itiTimes + stimTime);


%% Present stim

% Audio 
audioOnset = t.map(@(t) sum(t > audio_startTimes)).skipRepeats;

toneSamples = audioOnset.at(audioOnset).map(@(x) ...
    aud.pureTone(audio_stim_shuffle(x),stimTime,audioSampleRate, ...
    audioRampDuration,audioChannels));
audio.playTone = toneSamples;

% Visual
visualOnset = t.map(@(t) sum(t > visual_startTimes)).skipRepeats;

stimFlicker = mod(skipRepeats(floor((t - t.at(visualOnset))/(1/stimFlickerFrequency))),2);
stim = vis.grating(t, 'square', 'gaussian');
stim.spatialFrequency = spatialFrequency;
stim.contrast = contrast;
stim.orientation = orientation;
stim.phase = pi*stimFlicker;

stim.azimuth = visualOnset.at(visualOnset).map(@(x) visual_stim_shuffle{x}(1));
stim.elevation = visualOnset.at(visualOnset).map(@(x) visual_stim_shuffle{x}(2));
stim.sigma = visualOnset.at(visualOnset).map(@(x) visual_stim_shuffle{x}(3:4));
stim.show = visualOnset.to(visualOnset.delay(stimTime));
visStim.stim = stim;


endTrial = t.ge(max(audio_startTimes(end),visual_startTimes(end))).skipRepeats;

%% Events

events.audioOnset = audioOnset;
events.audioParams = events.expStart.map(@(x) audio_stim_shuffle);

events.visualOnset = visualOnset;
events.visualParams = events.expStart.map(@(x) visual_stim_shuffle);

events.endTrial = endTrial;

end

















