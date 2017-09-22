function AP_singleStimLickReward_fullfield(t, events, parameters, visStim, inputs, outputs, audio)
% AP 2017-06-23: pair visual stimuli with rewards

%% Fixed parameters

% Reward
rewardSize = 3;

% Stimuli
stimAzimuths = [0];
rewardedStim = 0; % which stimulus is paired with reward
spatialFrequency = 0.01;
contrast = 1;
orientation = 45;
sigma = [360,360];
stimFlickerFrequency = 5;

% Timing
stimTime = 30; % time stimulus is on the screen
itiTimes = 2:4;

% Lick detector
lick = inputs.lickDetector.delta.skipRepeats;

%% Set up trial parameters and events

% Choose trial parameters (random: azimuth and ITI)
trialAzimuth = events.newTrial.mapn(@(x) stimAzimuths(randperm(length(stimAzimuths),1)));
trialITI = events.newTrial.mapn(@(x) itiTimes(randperm(length(itiTimes),1)));

% Turn stim on for fixed interval 
% (delay this by 10 ms because concurrent events causes a crash)
stimOn = delay(at(true,events.newTrial),0.01); 

% Reward on rewarded stimulus on lever press
timeout = skipRepeats(ge((t-t.at(stimOn)),stimTime));
hit = stimOn.setTrigger(lick.eq(1) & trialAzimuth.eq(rewardedStim)).keepWhen(stimOn.to(timeout));
miss = stimOn.setTrigger(lick.eq(1) & ~trialAzimuth.eq(rewardedStim)).keepWhen(stimOn.to(timeout));

% Turn stim off after a timeout or a hit
stimOff = stimOn.setTrigger(merge(hit,timeout));

water = at(rewardSize,hit);  
outputs.reward = water;
totalWater = water.scan(@plus,0);

% Start the ITI after hit or the stimulus turns off
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
stim.orientation = orientation;
stim.show = stimOn.to(stimOff);

visStim.stim = stim;

%% Events

events.stimOn = stimOn;
events.stimOff = stimOff;
events.stimFlicker = stimFlicker;
events.trialAzimuth = trialAzimuth;
events.totalWater = totalWater;
events.lick = lick;
events.hit = hit;
events.miss = miss;
events.timeout = timeout;

events.endTrial = endTrial;























