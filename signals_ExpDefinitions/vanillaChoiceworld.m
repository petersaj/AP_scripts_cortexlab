function vanillaChoiceworld(t, events, parameters, visStim, inputs, outputs, audio)
% vanillaChoiceworld(t, events, parameters, visStim, inputs, outputs, audio)
% 170309 - AP
%
% Choice world that adapts with behavior
%
% Task structure: 
% Start trial
% Resetting pre-stim quiescent period
% Stimulus onset
% Fixed cue interactive delay
% Infinite time for response, fix stim azimuth on response
% Short ITI on reward, long ITI on punish, then turn stim off
% End trial
% 
% TO DO: finish writing properties

%% Fixed parameters

% Reward
rewardSize = 0.2; %%%% TO DO: CHECK OTHER PROTOCOLS

% Trial choice parameters
% Staircase trial choice
% (how often staircase trials appear - every staircaseTrials trials)
staircaseTrials = 2; 
% (how many hits to move forward on the staircase)
staircaseHit = 3;
% (how many misses to move backward on the staircase)
staircaseMiss = 1;

% Stimulus/target
% (which contrasts to use)
contrasts = [1,0.5,0.25,0.125,0.06,0];
% (which conrasts to use at the beginning of training)
startingContrasts = [true,true,false,false,false,false];
% (which contrasts to repeat on miss)
repeatOnMiss = [true,true,false,false,false,false];
% (number of trials to judge rolling performance)
trialsToBuffer = 10; %%%% TO DO: make a reasonable number here (200?)
% (number of trials after introducing 12.5% contrast to introduce 0%)
trialsToZeroContrast = 500;
sigma = [9,9];
spatialFrequency = 0.01;
startingAzimuth = 90;
responseDisplacement = 90;

% Timing
prestimQuiescentTime = 0.5;
cueInteractiveDelay = 0.5;
itiHit = 1;
itiMiss = 2;

% Sounds
audio_sample_rate = 192e3; %%%% TO DO: CHECK OTHER PROTOCOLS

onsetToneAmplitude = 1; %%%% TO DO: CHECK OTHER PROTOCOLS
onsetToneFreq = 500; %%%% TO DO: CHECK OTHER PROTOCOLS
onsetToneDuration = 0.1; %%%% TO DO: CHECK OTHER PROTOCOLS
toneSamples = onsetToneAmplitude*events.expStart.map(@(x) ...
    aud.pureTone(onsetToneFreq, onsetToneDuration, audio_sample_rate));

missNoiseDuration = itiMiss;
missNoiseAmplitude = 0.2; %%%% TO DO: CHECK OTHER PROTOCOLS
missNoiseSamples = missNoiseAmplitude*events.expStart.map(@(x) ...
    randn(2, audio_sample_rate*missNoiseDuration));

% Wheel parameters
quiescThreshold = 1; % what's a reasonable value for this?
wheelGain = 1; % I guess ultimately defined per rig...

%% Conditional parameters

subject = 'AP001'; %%%% TO DO: placeholder, get this from MC

% Initialize new or load old performance
trialDataInit = events.expStart.mapn( ...
    subject,contrasts,startingContrasts,repeatOnMiss, ...
    trialsToBuffer,trialsToZeroContrast,staircaseTrials,staircaseHit,staircaseMiss, ...
    @initializeTrialData).subscriptable;

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
audio.onsetTone = toneSamples.at(interactiveOn);

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
outputs.reward = at(rewardSize,trialData.hit);  

% Play noise on miss
audio.missNoise = missNoiseSamples.at(trialData.miss);

% ITI defined by outcome
iti = response.delay(trialData.hit.at(response)*itiHit + trialData.miss.at(response)*itiMiss);

% Stim stays on until the end of the ITI
stimOff = iti;

% End trial at the end of the ITI
endTrial = iti;

%% Visual stimulus

% Azimuth control
% Stim fixed in place before interactive and after response, wheel-conditional otherwise
stimAzimuth = cond( ...
    events.newTrial.to(interactiveOn), startingAzimuth*trialData.trialSide, ...
    interactiveOn.to(response), startingAzimuth*trialData.trialSide + stimDisplacement);

stim = vis.grating(t, 'square', 'gaussian');
stim.sigma = sigma;
stim.spatialFrequency = spatialFrequency;
stim.phase = 2*pi*events.newTrial.map(@(v)rand);
stim.azimuth = stimAzimuth;
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
events.conditions = trialData.conditions;
events.trialContrast = trialData.trialContrast;
events.trialSide = trialData.trialSide;
events.repeatTrial = trialData.repeatTrial;
events.hit = trialData.hit;
events.miss = trialData.miss;
events.staircase = trialData.staircase;
events.use_contrasts = trialData.use_contrasts;
events.hit_buffer = trialData.hit_buffer;
events.trialsToZeroContrast = trialData.trialsToZeroContrast;

end

function trialDataInit = initializeTrialData(~, ...
    subject,contrasts,startingContrasts,repeatOnMiss,trialsToBuffer, ...
    trialsToZeroContrast,staircaseTrials,staircaseHit,staircaseMiss)

%%%% Initialize all of the session-independent performance values
trialDataInit = struct;

% Store the contrasts which are used
trialDataInit.contrasts = contrasts;
% Store which trials are repeated on miss
trialDataInit.repeatOnMiss = repeatOnMiss;
% Define conditions as side*contrast
trialDataInit.conditions = unique(sort([contrasts,-contrasts]));
% Set the first contrast to 1
trialDataInit.trialContrast = 1;
% Set the first trial side randomly
trialDataInit.trialSide = randsample([-1,1],1);
% Set up the flag for repeating incorrect
trialDataInit.repeatTrial = false;
% Initialize hit/miss
trialDataInit.hit = false;
trialDataInit.miss = false;
% Initialize the staircase: 
% [current contrast, hits, misses, staircase trial counter, 
% staircase every n trials, hit requirement, miss requirement]
trialDataInit.staircase = ...
    [contrasts(1),0,0,0, ...
    staircaseTrials,staircaseHit,staircaseMiss];

n_conditions = length(trialDataInit.conditions);

%%%% Load the last experiment for the subject if it exists
% (note: MC creates folder on initilization, so look for > 1)
expRef = dat.listExps(subject);
if length(expRef) > 1
    previousBlockFilename = dat.expFilePath(expRef{end-1}, 'block', 'master');
    if exist(previousBlockFilename,'file')
        previousBlock = load(previousBlockFilename);
    end
end

if exist('previousBlock','var') && all(isfield(previousBlock.block, ...
        {'use_contrasts','hit_buffer','trialsToZeroContrast'}))
    % If the last experiment file has the relevant fields, set up performance
    
    % Which contrasts are currently in use
    trialDataInit.use_contrasts = previousBlock.block.use_contrasts(end-n_conditions+1:end);
    % The buffer to judge recent performance for adding contrasts
    trialDataInit.hit_buffer = reshape( ...
        previousBlock.block.hit_buffer(end-trialsToBuffer*n_conditions+1:end), ...
        trialsToBuffer,n_conditions);
    % The countdown to adding 0% contrast
    trialDataInit.trialsToZeroContrast = previousBlock.block.trialsToZeroContrast(end);      
    
else    
    % If this animal has no previous experiments, initialize performance
    trialDataInit.use_contrasts = startingContrasts;
    trialDataInit.hit_buffer = nan(trialsToBuffer,n_conditions);
    trialDataInit.trialsToZeroContrast = trialsToZeroContrast;  
end

end

function trialData = updateTrialData(trialData,stimDisplacement)
% Update the performance and pick the next contrast

%%%% Define the current trial condition
currentTrialCondition = trialData.trialSide*trialData.trialContrast;
conditionIdx = trialData.conditions == currentTrialCondition;

%%%% Define response type based on trial condition
trialData.hit = stimDisplacement*trialData.trialSide < 0;
trialData.miss = stimDisplacement*trialData.trialSide > 0;

%%%% Update hit/miss buffer if not a repeat trial
if ~trialData.repeatTrial
    trialData.hit_buffer(:,conditionIdx) = ...
        [trialData.hit;trialData.hit_buffer(1:end-1,conditionIdx)];
    trialsToBuffer = size(trialData.hit_buffer,1);
end

%%%% Set flag to repeat - skip trial choice if so, choose trial if not
if trialData.miss && ...
        ismember(trialData.trialContrast,trialData.contrasts(trialData.repeatOnMiss))
    trialData.repeatTrial = true;
    return
else
    trialData.repeatTrial = false;
end

%%%% Add new contrasts as necessary given performance
% This is based on the last trialsToBuffer trials for rolling performance
% (these parameters are hard-coded because too specific)
% (these are side-independent)
current_min_contrast = min(trialData.contrasts(trialData.use_contrasts & trialData.contrasts ~= 0));
switch current_min_contrast
    
    case 0.5
        % Lower from 0.5 contrast after > 75% correct
        min_hit_percentage = 0.75;
        
        curr_condition = ismember(abs(trialData.conditions),[0.5,1]);
        condition_total_trials = sum(sum(~isnan(trialData.hit_buffer(:,curr_condition))));
        % If there have been enough buffer trials, check performance
        if condition_total_trials >= size(trialData.hit_buffer,1)
            % Sample as evenly as possible across pooled conditions
            pooled_hits = reshape(trialData.hit_buffer(:,curr_condition)',[],1);
            use_hits = sum(pooled_hits(find(~isnan(pooled_hits),trialsToBuffer)));
            min_hits = find(1 - binocdf(1:trialsToBuffer,trialsToBuffer,min_hit_percentage) < 0.05,1);
            if use_hits >= min_hits
                trialData.use_contrasts(find(~trialData.use_contrasts,1)) = true;
            end
        end

    case 0.25
        % Lower from 0.25 contrast after > 50% correct
        min_hit_percentage = 0.5;
        
        curr_condition = ismember(abs(trialData.conditions),current_min_contrast);
        condition_total_trials = sum(sum(~isnan(trialData.hit_buffer(:,curr_condition))));
        % If there have been enough buffer trials, check performance
        if condition_total_trials >= size(trialData.hit_buffer,1)
            % Sample as evenly as possible across pooled conditions
            pooled_hits = reshape(trialData.hit_buffer(:,curr_condition)',[],1);
            use_hits = sum(pooled_hits(find(~isnan(pooled_hits),trialsToBuffer)));
            min_hits = find(1 - binocdf(1:trialsToBuffer,trialsToBuffer,min_hit_percentage) < 0.05,1);
            if use_hits >= min_hits
                trialData.use_contrasts(find(~trialData.use_contrasts,1)) = true;
            end
        end
        
    case 0.125
        % Lower from 0.25 contrast after > 50% correct
        min_hit_percentage = 0.5;
        
        curr_condition = ismember(abs(trialData.conditions),current_min_contrast);
        condition_total_trials = sum(sum(~isnan(trialData.hit_buffer(:,curr_condition))));
        % If there have been enough buffer trials, check performance
        if condition_total_trials >= size(trialData.hit_buffer,1)
            % Sample as evenly as possible across pooled conditions
            pooled_hits = reshape(trialData.hit_buffer(:,curr_condition)',[],1);
            use_hits = sum(pooled_hits(find(~isnan(pooled_hits),trialsToBuffer)));
            min_hits = find(1 - binocdf(1:trialsToBuffer,trialsToBuffer,min_hit_percentage) < 0.05,1);
            if use_hits >= min_hits
                trialData.use_contrasts(find(~trialData.use_contrasts,1)) = true;
            end
        end          
        
end

% Add 0 contrast after trialsToZeroContrast trials with 0.125 contrast
if min(trialData.contrasts(trialData.use_contrasts)) <= 0.125 && ...
        trialData.trialsToZeroContrast > 0
    % Subtract one from the countdown
    trialData.trialsToZeroContrast = trialData.trialsToZeroContrast-1;
    % If at zero, add the 0 contrast condition
    if trialData.trialsToZeroContrast == 0
        trialData.use_contrasts(trialData.contrasts == 0) = true;
    end    
end

%%%% Pick next contrast (and update staircase on staircase trials)
staircaseTrial = trialData.staircase(4) == 0;

if ~staircaseTrial
    
    % Next contrast is random from current contrast set
    trialData.trialContrast = randsample(trialData.contrasts(trialData.use_contrasts),1);
    
elseif staircaseTrial
        
    % Update hit/miss counter
    trialData.staircase(2) = trialData.staircase(2) + trialData.hit;
    trialData.staircase(3) = trialData.staircase(3) + trialData.miss;
    
    % Update trial counter
    trialData.staircase(4) = trialData.staircase(4) + 1;
    if trialData.staircase(4) >= trialData.staircase(5)
        trialData.staircase(4) = 0;
    end
    
    % Move staircase on hit/miss counter threshold
    if trialData.staircase(2) >= trialData.staircase(6)
        % On hit threshold, move the staircase forward and reset hit/miss
        newStaircaseContrast = trialData.contrasts(...
            min(find(trialData.staircase(1) == trialData.contrasts)+1, ...
            sum(trialData.use_contrasts)));
        trialData.staircase(1) = newStaircaseContrast;
        trialData.staircase(2:3) = 0;
    elseif trialData.staircase(3) >= trialData.staircase(7)
        % On miss threshold, move staircase backward and reset hit/miss
        newStaircaseContrast = trialData.contrasts(...
            max(find(trialData.staircase(1) == trialData.contrasts)-1,1));
        trialData.staircase(1) = newStaircaseContrast;
        trialData.staircase(2:3) = 0;
    end
    
    % Next contrast is defined by the staircase
    trialData.trialContrast = trialData.staircase(1);
    
end

%%%% Pick next side (this is done at random)
trialData.trialSide = randsample([-1,1],1);

end


