function AP_stimReward(t, events, parameters, visStim, inputs, outputs, audio)
% AP 2017-06-23: pair visual stimuli with rewards

%% Fixed parameters

% Reward
rewardSize = 3;

% Stimuli
stimAzimuths = [-90,0,90];
spatialFreq = 1/15;
contrast = 1;
sigma = [90,90];
stimFlickerFrequency = 5;
orientation = 45;

% Timing
stimTime = 1.5; % time stimulus is on the screen
stimRewardTime = 1; % time after stimulus onset when reward is given
itiTimes = 5:7;

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

stimFlicker = mod(skipRepeats(floor((t - t.at(stimOn))/(1/stimFlickerFrequency))),2);
%stim.contrast = trialContrast.at(stimOn)*stimFlicker;

stim = vis.grating(t, 'square', 'gaussian');
stim.sigma = sigma;
stim.spatialFreq = spatialFreq;
stim.phase = pi*stimFlicker;
stim.azimuth = trialAzimuth.at(stimOn);
stim.contrast = contrast;
stim.orientation = orientation;
stim.show = stimOn.to(stimOff);

visStim.stim = stim;

%% Events

events.stimOn = stimOn;
events.stimOff = stimOff;
events.stimFlicker = stimFlicker;
events.trialAzimuth = trialAzimuth;

events.endTrial = endTrial;























