function AP_stimLeverReward(t, events, parameters, visStim, inputs, outputs, audio)
% AP 2017-06-23: pair visual stimuli with rewards

%% Fixed parameters

% Reward
rewardSize = 3;

% Stimuli
stimAzimuths = [30];
rewardedStim = 30; % which stimulus is paired with reward
spatialFrequency = 0.01;
contrast = 1;
sigma = [5,90];
stimFlickerFrequency = 4;

% Timing
stimTime = 2; % time stimulus is on the screen
itiTimes = 5:7;

% Lever
lever_r_flip = inputs.lever_r.delta.skipRepeats;

%% Set up trial parameters and events

% Choose trial parameters (random: azimuth and ITI)
trialAzimuth = events.newTrial.mapn(@(x) stimAzimuths(randperm(length(stimAzimuths),1)));
trialITI = events.newTrial.mapn(@(x) itiTimes(randperm(length(itiTimes),1)));

% Turn stim on for fixed interval 
% (delay this by 10 ms because concurrent events causes a crash)
stimOn = delay(at(true,events.newTrial),0.01); 

% Reward on rewarded stimulus on lever press
timeout = skipRepeats(ge((t-t.at(stimOn)),stimTime));
hit = stimOn.setTrigger(lever_r_flip.eq(1) & trialAzimuth.eq(rewardedStim)).keepWhen(stimOn.to(timeout));
miss = stimOn.setTrigger(lever_r_flip.eq(1) & ~trialAzimuth.eq(rewardedStim)).keepWhen(stimOn.to(timeout));

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
stim.show = stimOn.to(stimOff);

visStim.stim = stim;

%% Events

events.stimOn = stimOn;
events.stimOff = stimOff;
events.stimFlicker = stimFlicker;
events.trialAzimuth = trialAzimuth;
events.totalWater = totalWater;
events.lever_r_flip = lever_r_flip;
events.hit = hit;
events.miss = miss;
events.timeout = timeout;

events.endTrial = endTrial;























