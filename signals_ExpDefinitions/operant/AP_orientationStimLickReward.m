function AP_orientationStimLickReward(t, events, parameters, visStim, inputs, outputs, audio)

%% Parameters

% Reward
rewardSize = 4;

% Stimuli
stimOrientations = [-45,45];
rewardedStim = stimOrientations; % which stimuli is paired with reward
spatialFrequency = 0.04;
contrast = 1;
azimuth = 0;
sigma = [20,20];
stimFlickerFrequency = 5;

% Timing (these are MC-parameterized now)
parameters.minLickStimTime; % = 2; % minimum time from last lick to show stim
parameters.stimTime; % =  3; % time stimulus is on the screen
parameters.minITI; % = 3;
parameters.maxITI; % = 5;

% Lick detector
lick = inputs.lickDetector.delta.skipRepeats;

%% Set up trial parameters and events

% Choose trial parameters (random: condition and ITI)
trialOrientation = events.newTrial.mapn(@(x) stimOrientations(randperm(length(stimOrientations),1)));
%trialITI = events.newTrial.mapn(@(x) randi(currentMinITI,currentMaxITI,1));
trialITI = at(mapn(parameters.minITI,parameters.maxITI,@(minITI,maxITI) randi([minITI,maxITI])),events.newTrial);

% Turn stim on for fixed interval 
% (delay this by 10 ms because concurrent events causes a crash)
stimOn = delay(at(true,events.newTrial),0.01); 

% Reward on rewarded stimulus on lever press
stimTimeClock = skipRepeats(ge((t-t.at(stimOn)),parameters.stimTime));
hit = stimOn.setTrigger(lick.eq(1) & trialOrientation.map(@(x) ismember(x,rewardedStim))).keepWhen(stimOn.to(stimTimeClock));
falseAlarm = stimOn.setTrigger(lick.eq(1) & ...
    ~trialOrientation.map(@(x) ismember(x,rewardedStim))).keepWhen(stimOn.to(stimTimeClock));

% Turn stim off after a timeout or a hit
%stimOff = stimOn.setTrigger(merge(hit,timeout));
% Turn stim off after timeout
stimOff = stimOn.setTrigger(stimTimeClock);

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
endTrial = stimOff.setTrigger(ge(t-t.at(stimOff),trialITI.at(stimOff)) & ge(lastLickTime,parameters.minLickStimTime));

%% Visual stimuli

stimFlicker = mod(skipRepeats(floor((t - t.at(stimOn))/(1/stimFlickerFrequency))),2);
%stim.contrast = trialContrast.at(stimOn)*stimFlicker;

stim = vis.grating(t, 'square', 'gaussian');
stim.sigma = sigma;
stim.spatialFrequency = spatialFrequency;
stim.phase = pi*stimFlicker;
stim.orientation = trialOrientation.at(stimOn);
stim.contrast = contrast;
stim.azimuth = azimuth;
stim.show = stimOn.to(stimOff);

visStim.stim = stim;

%% Events

events.stimOn = stimOn;
events.stimOff = stimOff;
events.stimFlicker = stimFlicker;
events.trialOrientation = trialOrientation;
events.trialITI = trialITI;
events.totalWater = totalWater;
events.lick = lick;
events.hit = hit;
events.falseAlarm = falseAlarm;

events.n_licks = n_licks;

events.endTrial = endTrial;























