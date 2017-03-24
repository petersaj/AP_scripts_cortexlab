function adaptiveChoiceWorld3(t, events, pars, visStim, inputs, outputs, audio)
% adaptiveChoiceWorld3(t, evts, pars, vs, in, out, audio)
% Choice world that adapts with behavior
% 170309 - MW/AP

% CURRENT PLACE NOTES
% RESPONSE/ITI is currently broken - this needs to change so that the scan
% function doesn't magically update on the new trial for no known reason

%% Fixed parameters

% Stimulus/target
sigma = [9,9];
spatialFrequency = 0.01;
startingAzimuth = 90;
hitDisplacement = 90;
missDisplacement = 90;

% Timing
prestimQuiescentTime = 0.5;
cueInteractiveDelay = 0.5;
itiReward = 1;
itiPunish = 2;

% Wheel parameters
quiescThreshold = 1; % what's a reasonable value for this?
wheelGain = 1; % I guess ultimately defined per rig...

%% Conditional parameters

% Initialize performance structure
contrasts = [1,0.5,0.25,0.125,0.06,0];
performanceInit = struct;
performanceInit.contrasts = contrasts;
% (use only 100 and 50 from the start)
performanceInit.use_contrasts = [true,true,false,false,false,false];
performanceInit.conditions = unique(sort([contrasts,-contrasts]));
performanceInit.n_trials = zeros(size(performanceInit.conditions));
performanceInit.n_correct = zeros(size(performanceInit.conditions));
% (convert the performance structure to a signal on experiment start)
performanceInit = events.expStart.map(@(x) performanceInit).subscriptable;

% Choose trial condition
trialCondition = events.newTrial.mapn(performanceInit,@chooseTrialCondition);

%% Set up wheel 

wheel = inputs.wheel.skipRepeats();

%% Trial event times

% Task structure: 
% Start trial
% Resetting pre-stim quiescent period
% Stimulus onset
% Fixed cue interactive delay
% Infinite time for response, fix stim azimuth on response
% Short ITI on reward, long ITI on punish, then turn stim off
% End trial

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

hit = keepWhen(interactiveOn.setTrigger(stimDisplacement*sign(trialCondition) ...
    <= -hitDisplacement),interactiveOn.to(events.newTrial));
miss = keepWhen(interactiveOn.setTrigger(stimDisplacement*sign(trialCondition) ...
    >= missDisplacement),interactiveOn.to(events.newTrial));

% ITI defined by outcome
iti = merge(hit.delay(itiReward),miss.delay(itiPunish));

% Stim fixed in place before interactive and after response, wheel-conditional otherwise
stimAzimuth = cond( ...
    events.newTrial.to(interactiveOn), startingAzimuth*sign(trialCondition), ...
    interactiveOn.to(merge(hit,miss)), startingAzimuth*sign(trialCondition) + stimDisplacement, ...
    hit.to(events.newTrial), (startingAzimuth-hitDisplacement)*sign(trialCondition), ...
    miss.to(events.newTrial), (startingAzimuth+missDisplacement)*sign(trialCondition));

% Stim stays on until the end of the ITI
stimOff = iti;

% End trial at the end of the ITI
endTrial = iti;

%% Update performance during ITI
% NOTE: this cannot be done at endTrial: concurrent events break things

outcome = merge( ...
    at(true,hit), ...
    at(false,miss));

trialInfo = ...
    [trialCondition.at(merge(hit,miss)), ...
    outcome.at(merge(hit,miss))];

performance = trialInfo.scan(@updatePerformance,performanceInit).subscriptable;


%% Visual stimulus

% Trial-independent parameters
stim = vis.grating(t, 'square', 'gaussian');
stim.sigma = sigma;
stim.spatialFrequency = spatialFrequency;
stim.phase = 2*pi*events.newTrial.map(@(v)rand);

% Position
stim.azimuth = stimAzimuth;

% Contrast
stim.contrast = abs(trialCondition);

% Stimulus on/off times
stim.show = stimOn.to(stimOff);

% Send stimulus object to the visual stimulus handler
visStim.target = stim;

%% Display and save
events.preStimQuiescence = preStimQuiescence;
events.stimOn = stimOn;
events.stimOff = stimOff;
events.interactiveOn = interactiveOn;
events.iti = iti;
events.stimAzimuth = stimAzimuth;
events.stimDisplacement = stimDisplacement;
events.endTrial = endTrial;

%at(true,interactiveOn.setTrigger(stimDisplacement*sign(trialAzimuth) <= hitAzimuth));
%at(false,interactiveOn.setTrigger(stimDisplacement*sign(trialAzimuth) >= missAzimuth));
end

function performance = updatePerformance(performance,trialInfo)
% Update the performance structure at the end of each trial

trialCondition = trialInfo(1);
hit = trialInfo(2);

performance.n_trials(performance.conditions == trialCondition) = ...
    performance.n_trials(performance.conditions == trialCondition) + 1;

performance.n_correct(performance.conditions == trialCondition) = ...
    performance.n_correct(performance.conditions == trialCondition) + hit;

end

function trialCondition = chooseTrialCondition(~,performance)
% Choose the next trial condition at the start of each trial

trialSide = randsample([-1,1],1);
trialContrast = randsample(performance.contrasts(performance.use_contrasts),1);

trialCondition = trialSide*trialContrast;

end







