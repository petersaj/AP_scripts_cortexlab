function AP_azimuthStimLickReward(t, events, parameters, visStim, inputs, outputs, audio)
% AP 2017-06-23: pair visual stimuli with rewards

%% Fixed parameters

% Reward
rewardSize = 3;

% Stimuli
stimAzimuths = [0,90];
rewardedStim = [90]; % which stimuli is paired with reward
spatialFrequency = 0.01;
contrast = 1;
orientation = 45;
sigma = [20,20];
stimFlickerFrequency = 5;

% Timing
minLickStimTime = 2; % minimum time from last lick to show stim
stimTime = 3; % time stimulus is on the screen
itiTimes = 10:20;

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
hit = stimOn.setTrigger(lick.eq(1) & trialAzimuth.map(@(x) ismember(x,rewardedStim))).keepWhen(stimOn.to(timeout));
miss = stimOn.setTrigger(lick.eq(1) & ~trialAzimuth.map(@(x) ismember(x,rewardedStim))).keepWhen(stimOn.to(timeout));

% Turn stim off after a timeout or a hit
%stimOff = stimOn.setTrigger(merge(hit,timeout));
% Turn stim off after timeout
stimOff = stimOn.setTrigger(timeout);

water = at(rewardSize,hit);  
outputs.reward = water;
totalWater = water.scan(@plus,0);

% ITI is time after stimulus off only
%endTrial = stimOff.delay(trialITI.at(stimOff));
% ITI is time after stim off but resets if lick during ITI
% define the last lick time as infinite if no licks yet
n_licks = lick.ge(1).scan(@plus,0);
lastLickTime = cond(...
    n_licks.gt(0), t-t.at(lick),...
    true, Inf);
endTrial = stimOff.setTrigger(ge(t-t.at(stimOff),trialITI.at(stimOff)) & ge(lastLickTime,minLickStimTime));

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
events.trialITI = trialITI;
events.totalWater = totalWater;
events.lick = lick;
events.hit = hit;
events.miss = miss;
events.timeout = timeout;

events.n_licks = n_licks;

events.endTrial = endTrial;























