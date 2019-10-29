function vanillaChoiceworldBias(t, events, parameters, visStim, inputs, outputs, audio)
% vanillaChoiceworldBias(t, events, parameters, visStim, inputs, outputs, audio)
%
% Choice world that only uses 2 contrasts and tries to balance bias
%
% Task structure: 
% Start trial
% Resetting pre-stim quiescent period
% Stimulus onset
% Fixed cue interactive delay
% Infinite time for response, fix stim azimuth on response
% Short ITI on reward, long ITI on punish, then turn stim off
% End trial


%% Fixed parameters

% Reward
rewardSize = 2;

% Stimulus/target
% (which contrasts to use)
contrasts = [1,0.5,0.25,0.125,0.06,0];
% (which conrasts to use at the beginning of training)
startingContrasts = [true,true,false,false,false,false];
% (which contrasts to repeat on miss)
repeatOnMiss = [true,true,false,false,false,false];
% (number of trials to judge rolling performance)
trialsToBuffer = 5;

sigma = [20,20];
spatialFreq = 1/15;
stimFlickerFrequency = 5; % DISABLED BELOW
startingAzimuth = 90;
responseDisplacement = 90;

% Timing
prestimQuiescentTime = 0.5;
cueInteractiveDelay = 0.5;
itiHit = 1;
itiMiss = 2;

% Sounds
% audioSampleRate = 192e3; % (Not used after new audio system 20181107)
audDev = audio.Devices('Strix');
audioSampleRate = audDev.DefaultSampleRate;
audioChannels = audDev.NrOutputChannels;

onsetToneAmplitude = 0.3; % Changed from 1 when moving to Strix
onsetToneFreq = 12000;
onsetToneDuration = 0.1;
onsetToneRampDuration = 0.01;
toneSamples = onsetToneAmplitude*events.expStart.map(@(x) ...
    aud.pureTone(onsetToneFreq,onsetToneDuration,audioSampleRate, ...
    onsetToneRampDuration,audioChannels));

missNoiseDuration = itiMiss;
missNoiseAmplitude = 0.03;
missNoiseSamples = missNoiseAmplitude*events.expStart.map(@(x) ...
    randn(2, audioSampleRate*missNoiseDuration));

% Wheel parameters
quiescThreshold = 1;
wheelGain = 2;

% Key press for manual reward
rewardKeyPressed = inputs.keyboard.strcmp('w');

%% Initialize trial data

trialDataInit = events.expStart.mapn( ...
    contrasts,startingContrasts,repeatOnMiss, ...
    trialsToBuffer,@initializeTrialData).subscriptable;

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

% Play tone at interactive onset
audio.Strix = toneSamples.at(interactiveOn);
% (20180711 changed from audio.onsetTone for new audio handling)

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
% NOTE: there is a 10ms delay for water output, because otherwise water and
% stim output compete and stim is delayed
water = at(rewardSize,merge(rewardKeyPressed,trialData.hit.delay(0.01)));  
outputs.reward = water;
totalWater = water.scan(@plus,0);

% Play noise on miss
audio.Strix = missNoiseSamples.at(trialData.miss.delay(0.01));
% (20180711 changed from audio.missNoise for new audio handling)

% ITI defined by outcome
iti = response.delay(trialData.hit.at(response)*itiHit + trialData.miss.at(response)*itiMiss);

% Stim stays on until the end of the ITI
stimOff = iti;

% End trial at the end of the ITI
endTrial = iti;

%% Visual stimulus

% Azimuth control
% 1) stim fixed in place until interactive on
% 2) wheel-conditional during interactive
% 3) fixed at response displacement azimuth after response
stimAzimuth = cond( ...
    events.newTrial.to(interactiveOn), startingAzimuth*trialData.trialSide, ...
    interactiveOn.to(response), startingAzimuth*trialData.trialSide + stimDisplacement, ...
    response.to(events.newTrial), ...
    startingAzimuth*trialData.trialSide.at(interactiveOn) + sign(stimDisplacement.at(response))*responseDisplacement);

% Stim flicker
stimFlicker = sin((t - t.at(stimOn))*stimFlickerFrequency*2*pi) > 0;

stim = vis.grating(t, 'square', 'gaussian');
stim.sigma = sigma;
stim.spatialFreq = spatialFreq;
stim.phase = 2*pi*events.newTrial.map(@(v)rand);
stim.azimuth = stimAzimuth;
%stim.contrast = trialContrast.at(stimOn)*stimFlicker;
stim.contrast = trialContrast.at(stimOn);
stim.show = stimOn.to(stimOff);

visStim.stim = stim;

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
events.repeatOnMiss = trialData.repeatOnMiss;
events.trialContrast = trialData.trialContrast;
events.trialSide = trialData.trialSide;
events.repeatTrial = trialData.repeatTrial;
events.hit = trialData.hit;
events.miss = trialData.miss;
events.useContrasts = trialData.useContrasts;
events.hitBuffer = trialData.hitBuffer;
events.bias = trialData.bias;
events.sessionPerformance = trialData.sessionPerformance;
events.totalWater = totalWater;

end

function trialDataInit = initializeTrialData(subject_info, ...
    contrasts,startingContrasts,repeatOnMiss,trialsToBuffer)

%%%% Get the subject
% (from events.expStart - MC gives subject after last underscore)
subject_info_underscore_idx = strfind(subject_info,'_');
if ~isempty(subject_info_underscore_idx)
    subject = subject_info(subject_info_underscore_idx(end)+1:end);
else
    % (if there are no underscores, set subject to nothing)
    subject = '';
end

%%%% Initialize all of the session-independent performance values
trialDataInit = struct;

% Store the contrasts which are used
trialDataInit.contrasts = contrasts;
% Store which trials are repeated on miss
trialDataInit.repeatOnMiss = repeatOnMiss;
% Set the first contrast to 1
trialDataInit.trialContrast = 1;
% Set the first trial side randomly
trialDataInit.trialSide = randsample([-1,1],1);
% Set up the flag for repeating incorrect
trialDataInit.repeatTrial = false;
% Initialize hit/miss
%%%%%%%%%%%% NOTE THAT THIS IS WEIRD (FIGURED THIS OUT LATE): 
%%%%%%%%%%%% This sets up the hit/miss signal AFTER a trial (i.e. the first
%%%%%%%%%%%% value is always undefined, the second value represents the
%%%%%%%%%%%% first trial). But this is confusing because contrasts are set
%%%%%%%%%%%% up BEFORE a trial (i.e. the first value is the first trial).
%%%%%%%%%%%% This means that in order to align them, in analysis they need
%%%%%%%%%%%% to be shifted one spot to the left. Maybe this could be done
%%%%%%%%%%%% here, but would be non-trivial, so doing in anaylsis.
trialDataInit.hit = false;
trialDataInit.miss = false;
% Initialize the day's performance to plot
% [conditions, number of trials, number of move left choices]
trialDataInit.sessionPerformance = ...
    [sort(unique([-contrasts,contrasts])); ...
    zeros(size(unique([-contrasts,contrasts]))); ...
    zeros(size(unique([-contrasts,contrasts])))];

% Initialize performance
trialDataInit.useContrasts = startingContrasts;
trialDataInit.hitBuffer = nan(trialsToBuffer,2);
trialDataInit.bias = 0;

end

function trialData = updateTrialData(trialData,stimDisplacement)
% Update the performance and pick the next contrast

%%%% Get current side
currentSideIdx = (trialData.trialSide == 1) + 1;

%%%% Define response type based on trial condition
trialData.hit = stimDisplacement*trialData.trialSide < 0;
trialData.miss = stimDisplacement*trialData.trialSide > 0;

%%%% Update buffers and counters if not a repeat trial
if ~trialData.repeatTrial
    
    %%%% Contrast-adding performance buffer
    % Update hit buffer for running performance
    trialData.hitBuffer(:,currentSideIdx) = ...
        [trialData.hit;trialData.hitBuffer(1:end-1,currentSideIdx)];       
    
    %%% Session performance for plotting
    currCondition = trialData.trialContrast*trialData.trialSide;
    currConditionIdx = trialData.sessionPerformance(1,:) == currCondition;
    trialData.sessionPerformance(2,currConditionIdx) = ...
        trialData.sessionPerformance(2,currConditionIdx) + 1;
    trialData.sessionPerformance(3,currConditionIdx) = ...
        trialData.sessionPerformance(3,currConditionIdx) + (stimDisplacement < 0);
   
    %%% Get bias
    curr_bias = diff(nanmean(trialData.hitBuffer,1));
    trialData.bias = curr_bias;

end

%%%% Set flag to repeat - skip trial choice if so
if trialData.miss && ...
        ismember(trialData.trialContrast,trialData.contrasts(trialData.repeatOnMiss))
    trialData.repeatTrial = true;
    return
else
    trialData.repeatTrial = false;
end

%%%% Pick next contrast at random
trialData.trialContrast = randsample(trialData.contrasts(trialData.useContrasts),1);

%%%% Next side is probabilistically drawn based on bias
curr_bias = diff(nanmean(trialData.hitBuffer,1));
if ~isnan(curr_bias)
    trialData.trialSide = sign(rand-0.5-curr_bias/2);
else
    trialData.trialSide = sign(rand-0.5);
end

end


