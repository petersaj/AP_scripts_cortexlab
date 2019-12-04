function AP_visAudioPassive(t, events, parameters, visStim, inputs, outputs, audio)
% AP 2017-09-22: present random tones and visual stimuli

%% Set up stimuli

% Make these static for now (dumb signals problems makes this difficult)
staticParameters.numRepeats = 100;
staticParameters.stimTime = 0.5;
staticParameters.minITI = 1;
staticParameters.maxITI = 2;
staticParameters.stepITI = 0.1;
staticParameters.volume = 0.3;
staticParameters.bufferTime = 10;

% Sound stim
audioFreqRange = linspace(5000,20000,3);
audioSampleRate = 192e3;

audioRampDuration = 0.01;
audioChannels = 2;

% Visual stim
vis_az_el_sig = {[-90,90,20,360],[0,90,20,360],[90,90,20,360]};
spatialFreq = 1/15;
contrast = 1;
stimFlickerFrequency = 5;
orientation = 45; 


%% Stim times

% Audio
audio_stim_uniform = repmat(audioFreqRange,staticParameters.numRepeats,1);
audio_stim_shuffle = audio_stim_uniform(randperm(numel(audio_stim_uniform)));

audio_itiTimes = randsample(staticParameters.minITI:staticParameters.stepITI:staticParameters.maxITI,length(audio_stim_shuffle),true);
audio_startTimes = staticParameters.bufferTime + cumsum(audio_itiTimes + staticParameters.stimTime);

% Visual
visual_stim_uniform = repmat(vis_az_el_sig,staticParameters.numRepeats,1);
visual_stim_shuffle = visual_stim_uniform(randperm(numel(visual_stim_uniform)));

visual_itiTimes = randsample(staticParameters.minITI:staticParameters.stepITI:staticParameters.maxITI,length(visual_stim_shuffle),true);
visual_startTimes = staticParameters.bufferTime + cumsum(visual_itiTimes + staticParameters.stimTime);


%% Present stim

% Audio 
audioOnset = t.map(@(t) sum(t > audio_startTimes)).skipRepeats;

toneSamples = audioOnset.at(audioOnset).map(@(x) ...
    aud.pureTone(audio_stim_shuffle(x),staticParameters.stimTime,audioSampleRate, ...
    audioRampDuration,audioChannels))*staticParameters.volume;
audio.playTone = toneSamples;

% Visual
visualOnset = t.map(@(t) sum(t > visual_startTimes)).skipRepeats;

stimFlicker = mod(skipRepeats(floor((t - t.at(visualOnset))/(1/stimFlickerFrequency))),2);
stim = vis.grating(t, 'square', 'gaussian');
stim.spatialFreq = spatialFreq;
stim.contrast = contrast;
stim.orientation = orientation;
stim.phase = pi*stimFlicker;

stim.azimuth = visualOnset.at(visualOnset).map(@(x) visual_stim_shuffle{x}(1));
stim.elevation = visualOnset.at(visualOnset).map(@(x) visual_stim_shuffle{x}(2));
stim.sigma = visualOnset.at(visualOnset).map(@(x) visual_stim_shuffle{x}(3:4));
stim.show = visualOnset.to(visualOnset.delay(staticParameters.stimTime));
visStim.stim = stim;

endTrial = t.ge(staticParameters.bufferTime + max(audio_startTimes(end),visual_startTimes(end))).skipRepeats;

%% Events

events.audioOnset = audioOnset;
events.audioParams = events.expStart.map(@(x) audio_stim_shuffle);

events.visualOnset = visualOnset;
events.visualParams = events.expStart.map(@(x) visual_stim_shuffle);

events.endTrial = endTrial;

end

















