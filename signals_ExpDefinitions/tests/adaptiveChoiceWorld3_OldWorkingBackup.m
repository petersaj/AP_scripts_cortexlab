function adaptiveChoiceWorld3(t, events, pars, visStim, inputs, outputs, audio)
% adaptiveChoiceWorld3(t, evts, pars, vs, in, out, audio)
% Choice world that adapts with behavior
% 170309 - MW/AP


%% Fixed parameters

% Trial choice parameters
% (staircase: every staircaseTrials, engage an staircaseHit forward / staircaseMiss back
% choice of contrast)
staircaseTrials = 2; 
staircaseHit = 3;
staircaseMiss = 1;

% Stimulus/target
contrasts = [1,0.5,0.25,0.125,0.06,0];
startingContrasts = [true,true,false,false,false,false];
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

% Initialize new or load old performance
performanceInit = events.expStart.mapn(contrasts,startingContrasts,@initializePerformance).subscriptable;

% Initialize staircase
staircase = events.trialNum.mapn(@(trialNum) mod(trialNum,staircaseTrials) == 0);

% Choose trial condition
trialCondition = events.trialNum.mapn( ...
    performanceInit,staircase.at(events.trialNum),@chooseTrialCondition);

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

% This is a dirty way to package things, but scan only accepts one argument
trialInfo = ...
    [trialCondition.at(merge(hit,miss)), ...
    outcome.at(merge(hit,miss)), ...
    staircase.at(merge(hit,miss)), ...
    staircaseHit,staircaseMiss];

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
events.staircase = staircase;

%at(true,interactiveOn.setTrigger(stimDisplacement*sign(trialAzimuth) <= hitAzimuth));
%at(false,interactiveOn.setTrigger(stimDisplacement*sign(trialAzimuth) >= missAzimuth));
end

function performanceInit = initializePerformance(~,contrasts,startingContrasts)

% Initialize performance structure
performanceInit = struct;
performanceInit.contrasts = contrasts;
performanceInit.use_contrasts = startingContrasts;
performanceInit.conditions = unique(sort([contrasts,-contrasts]));
performanceInit.n_trials = zeros(size(performanceInit.conditions));
performanceInit.n_correct = zeros(size(performanceInit.conditions));

% Initialize the staircase: [current contrast, hits, misses]
performanceInit.staircase = [contrasts(1),0,0];

end

function performance = updatePerformance(performance,trialInfo)
% Update the performance structure at the end of each trial

% Unpackage (have to be packaged: only one allowable input argument)
trialCondition = trialInfo(1);
hit = trialInfo(2);
staircaseTrial = trialInfo(3);
staircaseHit = trialInfo(4);
staircaseMiss = trialInfo(5);

% Number of trials in each condition
performance.n_trials(performance.conditions == trialCondition) = ...
    performance.n_trials(performance.conditions == trialCondition) + 1;

% Number of hits in each condition
performance.n_correct(performance.conditions == trialCondition) = ...
    performance.n_correct(performance.conditions == trialCondition) + hit;

% Update staircase on staircase trials
if staircaseTrial
    
    % Update hit/miss counter
    performance.staircase(2) = performance.staircase(2) + hit;
    performance.staircase(3) = performance.staircase(3) + ~hit;
    
    % Move staircase on hit/miss counter threshold
    if performance.staircase(2) >= staircaseHit
        % On hit threshold, move the staircase forward and reset hit/miss
        newStaircaseContrast = performance.contrasts(...
            min(find(performance.staircase(1) == performance.contrasts)+1, ...
            sum(performance.use_contrasts)));
        performance.staircase(1) = newStaircaseContrast;
        performance.staircase(2:3) = 0;
    elseif performance.staircase(3) >= staircaseMiss
        % On miss threshold, move staircase backward and reset hit/miss
        newStaircaseContrast = performance.contrasts(...
            max(find(performance.staircase(1) == performance.contrasts)-1,1));
        performance.staircase(1) = newStaircaseContrast;
        performance.staircase(2:3) = 0;
    end
    
end

end

function trialCondition = chooseTrialCondition(~,performance,staircaseTrial)
% Choose the next trial condition at the start of each trial

% Choose a side randomly
trialSide = randsample([-1,1],1);

% Alternate random and staircase contrast choice (to more often sample the
% harder trials given the performance)
if ~staircaseTrial
    trialContrast = randsample(performance.contrasts(performance.use_contrasts),1);      
elseif staircaseTrial
    trialContrast = performance.staircase(1);    
end

trialCondition = trialSide*trialContrast;

end







