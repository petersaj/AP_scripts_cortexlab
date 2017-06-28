function AP_stimReward(t, events, parameters, visStim, inputs, outputs, audio)
% AP 2017-06-23: pair visual stimuli with rewards

%% Fixed parameters

% Reward (chosen to give ~400 ul in 30 min - otherwise they stop early)
rewardSize = 3;

% Stimuli
stimAzimuths = [-30,0,30];
rewardedStim = 30; % which stimulus is paired with reward
spatialFrequency = 0.01;
contrast = 1;
sigma = [5,90];
stimFlickerFrequency = 4;

% Timing
stimTime = 1.5; % time stimulus is on the screen
stimRewardTime = 1; % time after stimulus onset when reward is given
itiTimes = 2:4;

%% Set up trial parameters and events

% Choose trial parameters (random: azimuth and ITI)
trialAzimuth = events.newTrial.mapn(@(x) randsample(stimAzimuths,1));
trialITI = events.newTrial.mapn(@(x) randsample(itiTimes,1));

% Turn stim on for fixed interval 
% (delay this by 10 ms because concurrent events causes a crash)
stimOn = delay(at(true,events.newTrial),0.01); 
stimOff = stimOn.delay(stimTime);

% Reward on rewarded stimulus after set time
water = at(rewardSize*trialAzimuth.eq(rewardedStim),stimOn.delay(stimRewardTime));  
outputs.reward = water;
totalWater = water.scan(@plus,0);

% Start the ITI after the stimulus turns off
endTrial = stimOff.delay(trialITI.at(stimOff));

%% Visual stimuli

stimFlicker = mod(skipRepeats(floor((t - t.at(stimOn))/(1/stimFlickerFrequency))),2);
%stim.contrast = trialContrast.at(stimOn)*stimFlicker;

stim = vis.grating(t, 'square', 'gaussian');
stim.sigma = sigma;
stim.spatialFrequency = spatialFrequency;
stim.phase = pi*stimFlicker;
stim.azimuth = trialAzimuth.at(stimOn);
stim.contrast = contrast;
stim.show = stimOn.to(stimOff);

visStim.stim = stim;

%% Events

events.stimOn = stimOn;
events.stimOff = stimOff;
events.stimFlicker = stimFlicker;
events.trialAzimuth = trialAzimuth;
events.totalWater = totalWater;

events.endTrial = endTrial;























