function AP_opticFlow(t, events, parameters, visStim, inputs, outputs, audio)
% AP 2017-06-23: pair visual stimuli with rewards

%% Fixed parameters

% Stimuli
stimAzimuth = 180;
spatialFreq = 1/15;
contrast = 1;
sigma = [100,100];

temporalFreqs = [1,2,5];

% Timing
stimTime = 2; % time stimulus is on the screen
itiTimes = 2;

%% Set up trial parameters and events

% Choose random stimulus on each trial
trialTemporalFreq = events.newTrial.mapn(@(x) randsample(temporalFreqs,1));

% Turn stim on for fixed interval
stimOn = at(true,events.newTrial); 
stimOff = stimOn.delay(stimTime);

% Start the ITI after the stimulus turns off
trialITI = events.newTrial.mapn(@(x) randsample(itiTimes,1));
endTrial = stimOff.delay(trialITI);

%% Visual stimuli

stimPhase = t*2*pi*trialTemporalFreq.at(events.newTrial);

% Set up L/R stim
stimR = vis.grating(t, 'square', 'gaussian');
stimR.sigma = sigma;
stimR.spatialFreq = spatialFreq;
stimR.phase = stimPhase;
stimR.contrast = contrast;
stimR.azimuth = stimAzimuth;
stimR.show = stimOn.to(stimOff);

stimL = vis.grating(t, 'square', 'gaussian');
stimL.sigma = sigma;
stimL.spatialFreq = spatialFreq;
stimL.phase = -stimPhase;
stimL.contrast = contrast;
stimL.azimuth = -stimAzimuth;
stimL.show = stimOn.to(stimOff);

% Black out the center screen
stimC = vis.patch(t, 'rect');
stimC.azimuth = 0;
stimC.altitude = 0;
stimC.dims = [90,180];
stimC.colour = [0,0,0];
stimC.show = events.expStart;

% Package into visual stim
visStim.stimL = stimL;
visStim.stimR = stimR;
% visStim.stimC = stimC; 

%% Events

events.stimOn = stimOn;
events.stimOff = stimOff;
events.trialTemporalFreq = trialTemporalFreq;

events.endTrial = endTrial;























