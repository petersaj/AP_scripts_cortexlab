function AP_stimPassive(t, events, parameters, visStim, inputs, outputs, audio)
% AP 2017-06-23: show passive gratings

%% Fixed parameters

% Stimuli
stimAzimuths = [-90,0,90];
spatialFreq = 1/15;
contrast = 0;
sigma = [20,20];

% Timing
stimTime = 0.5; % time stimulus is on the screen
itiTimes = 1:2;

%% Set up trial parameters and events

% Choose random stimulus on each trial
trialAzimuth = events.newTrial.mapn(@(x) randsample(stimAzimuths,1));

% Turn stim on for fixed interval
stimOn = at(true,events.newTrial); 
stimOff = stimOn.delay(stimTime);

% Start the ITI after the stimulus turns off
trialITI = events.newTrial.mapn(@(x) randsample(itiTimes,1));
endTrial = stimOff.delay(trialITI);

%% Visual stimuli

stim = vis.grating(t, 'square', 'gaussian');
stim.sigma = sigma;
stim.spatialFreq = spatialFreq;
stim.phase = 2*pi*events.newTrial.map(@(v)rand);
stim.azimuth = trialAzimuth.at(stimOn);
stim.contrast = contrast;
stim.show = stimOn.to(stimOff);

visStim.stim = stim;

%% Events

events.stimOn = stimOn;
events.stimOff = stimOff;
events.trialAzimuth = trialAzimuth;

events.endTrial = endTrial;























