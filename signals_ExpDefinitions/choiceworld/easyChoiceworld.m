function easyChoiceworld(t, events, parameters, visStim, inputs, outputs, audio)
% easyChoiceworld(t, events, parameters, visStim, inputs, outputs, audio)
% 170510 - AP
%
% Easy choiceworld to build wheel-stim to center-water association first
%
% - prestim quiescent period
% - 2 visual stimuli, same high contrast but varying across trials
% - going left or right always gets reward

%% Fixed parameters

% Reward
rewardSize = 3;

% Stimulus/target
% (which contrasts to use)
contrasts = [1,0.5];
% (stim parameters)
sigma = [15,15];
spatialFreq = 0.01;
startingAzimuth = 90;
responseDisplacement = 90;

% Timing
prestimQuiescentTime = 0.5;
cueInteractiveDelay = 0;
iti = 5;

% Wheel parameters
quiescThreshold = 1;
wheelGain = 2;

%% Initialize trial data

trialDataInit = events.expStart.mapn(contrasts,@initializeTrialData).subscriptable;

%% Set up wheel 

wheel = inputs.wheel.skipRepeats();

%% Trial event times
% (this is set up to be independent of trial conditon, that way the trial
% condition can be chosen in a performance-dependent manner)

% Resetting pre-stim quiescent period
prestimQuiescentPeriod = at(prestimQuiescentTime,events.newTrial.delay(0)); 
preStimQuiescence = sig.quiescenceWatch(prestimQuiescentPeriod, t, wheel, quiescThreshold); 

% Stimulus onset
%stimOn = sig.quiescenceWatch(preStimQuiescPeriod, t, wheel, quiescThreshold); 
stimOn = at(true,preStimQuiescence); 

% Fixed cue interactive delay
interactiveOn = stimOn.delay(cueInteractiveDelay); 

% Response
% (wheel displacement zeroed at interactiveOn)
stimDisplacement = wheelGain*(wheel - wheel.at(interactiveOn));

response = keepWhen(interactiveOn.setTrigger(abs(stimDisplacement) ...
    >= responseDisplacement),interactiveOn.to(events.newTrial));

%% Update performance at response
% (NOTE: this cannot be done at endTrial: concurrent events break things)

% Update performance
trialData = stimDisplacement.at(response).scan(@updateTrialData,trialDataInit).subscriptable;

% Set trial contrast (chosen when updating performance)
trialContrast = trialData.trialContrast;


%% Give feedback and end trial

% Give reward on hit
water = at(rewardSize,response);  
outputs.reward = water;
totalWater = water.scan(@plus,0);

% ITI defined by outcome
iti = response.delay(iti);

% Stim stays on until the end of the ITI
stimOff = iti;

% End trial at the end of the ITI
endTrial = iti;

%% Visual stimulus

% Azimuth control
% Stim fixed in place before interactive and after response, wheel-conditional otherwise
stimAzimuth = cond( ...
    events.newTrial.to(interactiveOn), startingAzimuth, ...
    interactiveOn.to(response), startingAzimuth + stimDisplacement);

stim_right = vis.grating(t, 'square', 'gaussian');
stim_right.sigma = sigma;
stim_right.spatialFreq = spatialFreq;
stim_right.phase = 2*pi*events.newTrial.map(@(v)rand);
stim_right.azimuth = stimAzimuth;
stim_right.contrast = trialContrast.at(stimOn);
stim_right.show = stimOn.to(stimOff);

stim_left = vis.grating(t, 'square', 'gaussian');
stim_left.sigma = sigma;
stim_left.spatialFreq = spatialFreq;
stim_left.phase = 2*pi*events.newTrial.map(@(v)rand);
stim_left.azimuth = -startingAzimuth*2 + stimAzimuth;
stim_left.contrast = trialContrast.at(stimOn);
stim_left.show = stimOn.to(stimOff);

visStim.stim_left = stim_left;
visStim.stim_right = stim_right;

%% Display and save

% Wheel and stim
events.stimAzimuth = stimAzimuth;

% Trial times
events.stimOn = stimOn;
events.stimOff = stimOff;
events.interactiveOn = interactiveOn;
events.response = response;
events.endTrial = endTrial;

% Performance
events.contrasts = trialData.contrasts;
events.trialContrast = trialData.trialContrast;
events.totalWater = totalWater;

end

function trialDataInit = initializeTrialData(subject_info,contrasts)

% Store the contrasts which are used
trialDataInit.contrasts = contrasts;

% Set the first contrast to 1
trialDataInit.trialContrast = 1;

end

function trialData = updateTrialData(trialData,stimDisplacement)
% Pick the next contrast

% (randomly)
trialData.trialContrast = randsample(trialData.contrasts,1);

end


