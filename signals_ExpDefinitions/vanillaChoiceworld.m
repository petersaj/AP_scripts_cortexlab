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
performanceInit = events.expStart.mapn( ...
    subject,contrasts,startingContrasts,repeatOnMiss, ...
    trialsToBuffer,trialsToZeroContrast,staircaseTrials,staircaseHit,staircaseMiss, ...
    @initializePerformance).subscriptable;

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
performance = stimDisplacement.at(response).scan(@updatePerformance,performanceInit).subscriptable;

% Set trial contrast (chosen when updating performance)
trialContrast = performance.trialContrast;


%% Give feedback and end trial

% Give reward on hit
outputs.reward = at(rewardSize,performance.hit);  

% Play noise on miss
audio.missNoise = missNoiseSamples.at(performance.miss);

% ITI defined by outcome
iti = response.delay(performance.hit.at(response)*itiHit + performance.miss.at(response)*itiMiss);

% Stim stays on until the end of the ITI
stimOff = iti;

% End trial at the end of the ITI
endTrial = iti;

%% Visual stimulus

% Azimuth control
% Stim fixed in place before interactive and after response, wheel-conditional otherwise
stimAzimuth = cond( ...
    events.newTrial.to(interactiveOn), startingAzimuth*performance.trialSide, ...
    interactiveOn.to(response), startingAzimuth*performance.trialSide + stimDisplacement);

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
events.wheel = wheel;
events.stimAzimuth = stimAzimuth;

% Trial times
events.stimOn = stimOn;
events.stimOff = stimOff;
events.interactiveOn = interactiveOn;
events.response = response;
events.endTrial = endTrial;

% Performance
events.trialSide = performance.trialSide;
events.trialContrast = performance.trialContrast;
events.hit = performance.hit;
events.miss = performance.miss;
events.use_contrasts = performance.use_contrasts;
events.hit_buffer = performance.hit_buffer;
events.repeatTrial = performance.repeatTrial;
events.trialsToZeroContrast = performance.trialsToZeroContrast;

end

function performanceInit = initializePerformance(~, ...
    subject,contrasts,startingContrasts,repeatOnMiss,trialsToBuffer, ...
    trialsToZeroContrast,staircaseTrials,staircaseHit,staircaseMiss)

%%%% Initialize all of the session-independent performance values
performanceInit = struct;

% Store the contrasts which are used
performanceInit.contrasts = contrasts;
% Store which trials are repeated on miss
performanceInit.repeatOnMiss = repeatOnMiss;
% Define conditions as side*contrast
performanceInit.conditions = unique(sort([contrasts,-contrasts]));
% Set the first contrast to 1
performanceInit.trialContrast = 1;
% Set the first trial side randomly
performanceInit.trialSide = randsample([-1,1],1);
% Set up the flag for repeating incorrect
performanceInit.repeatTrial = false;
% Initialize hit/miss
performanceInit.hit = false;
performanceInit.miss = false;
% Initialize the staircase: 
% [current contrast, hits, misses, staircase trial counter, 
% staircase every n trials, hit requirement, miss requirement]
performanceInit.staircase = ...
    [contrasts(1),0,0,0, ...
    staircaseTrials,staircaseHit,staircaseMiss];

n_conditions = length(performanceInit.conditions);

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
    performanceInit.use_contrasts = previousBlock.block.use_contrasts(end-n_conditions+1:end);
    % The buffer to judge recent performance for adding contrasts
    performanceInit.hit_buffer = reshape( ...
        previousBlock.block.hit_buffer(end-trialsToBuffer*n_conditions+1:end), ...
        trialsToBuffer,n_conditions);
    % The countdown to adding 0% contrast
    performanceInit.trialsToZeroContrast = previousBlock.block.trialsToZeroContrast(end);      
    
else    
    % If this animal has no previous experiments, initialize performance
    performanceInit.use_contrasts = startingContrasts;
    performanceInit.hit_buffer = nan(trialsToBuffer,n_conditions);
    performanceInit.trialsToZeroContrast = trialsToZeroContrast;  
end

end

function performance = updatePerformance(performance,stimDisplacement)
% Update the performance and pick the next contrast

%%%% Define the current trial condition
currentTrialCondition = performance.trialSide*performance.trialContrast;
conditionIdx = performance.conditions == currentTrialCondition;

%%%% Define response type based on trial condition
performance.hit = stimDisplacement*performance.trialSide < 0;
performance.miss = stimDisplacement*performance.trialSide > 0;

%%%% Update hit/miss buffer if not a repeat trial
if ~performance.repeatTrial
    performance.hit_buffer(:,conditionIdx) = ...
        [performance.hit;performance.hit_buffer(1:end-1,conditionIdx)];
    trialsToBuffer = size(performance.hit_buffer,1);
end

%%%% Set flag to repeat - skip trial choice if so, choose trial if not
if performance.miss && ...
        ismember(performance.trialContrast,performance.contrasts(performance.repeatOnMiss))
    performance.repeatTrial = true;
    return
else
    performance.repeatTrial = false;
end

%%%% Add new contrasts as necessary given performance
% This is based on the last trialsToBuffer trials for rolling performance
% (these parameters are hard-coded because too specific)
% (these are side-independent)
current_min_contrast = min(performance.contrasts(performance.use_contrasts & performance.contrasts ~= 0));
switch current_min_contrast
    
    case 0.5
        % Lower from 0.5 contrast after > 75% correct
        min_hit_percentage = 0.75;
        
        curr_condition = ismember(abs(performance.conditions),[0.5,1]);
        condition_total_trials = sum(sum(~isnan(performance.hit_buffer(:,curr_condition))));
        % If there have been enough buffer trials, check performance
        if condition_total_trials >= size(performance.hit_buffer,1)
            % Sample as evenly as possible across pooled conditions
            pooled_hits = reshape(performance.hit_buffer(:,curr_condition)',[],1);
            use_hits = sum(pooled_hits(find(~isnan(pooled_hits),trialsToBuffer)));
            min_hits = find(1 - binocdf(1:trialsToBuffer,trialsToBuffer,min_hit_percentage) < 0.05,1);
            if use_hits >= min_hits
                performance.use_contrasts(find(~performance.use_contrasts,1)) = true;
            end
        end

    case 0.25
        % Lower from 0.25 contrast after > 50% correct
        min_hit_percentage = 0.5;
        
        curr_condition = ismember(abs(performance.conditions),current_min_contrast);
        condition_total_trials = sum(sum(~isnan(performance.hit_buffer(:,curr_condition))));
        % If there have been enough buffer trials, check performance
        if condition_total_trials >= size(performance.hit_buffer,1)
            % Sample as evenly as possible across pooled conditions
            pooled_hits = reshape(performance.hit_buffer(:,curr_condition)',[],1);
            use_hits = sum(pooled_hits(find(~isnan(pooled_hits),trialsToBuffer)));
            min_hits = find(1 - binocdf(1:trialsToBuffer,trialsToBuffer,min_hit_percentage) < 0.05,1);
            if use_hits >= min_hits
                performance.use_contrasts(find(~performance.use_contrasts,1)) = true;
            end
        end
        
    case 0.125
        % Lower from 0.25 contrast after > 50% correct
        min_hit_percentage = 0.5;
        
        curr_condition = ismember(abs(performance.conditions),current_min_contrast);
        condition_total_trials = sum(sum(~isnan(performance.hit_buffer(:,curr_condition))));
        % If there have been enough buffer trials, check performance
        if condition_total_trials >= size(performance.hit_buffer,1)
            % Sample as evenly as possible across pooled conditions
            pooled_hits = reshape(performance.hit_buffer(:,curr_condition)',[],1);
            use_hits = sum(pooled_hits(find(~isnan(pooled_hits),trialsToBuffer)));
            min_hits = find(1 - binocdf(1:trialsToBuffer,trialsToBuffer,min_hit_percentage) < 0.05,1);
            if use_hits >= min_hits
                performance.use_contrasts(find(~performance.use_contrasts,1)) = true;
            end
        end          
        
end

% Add 0 contrast after trialsToZeroContrast trials with 0.125 contrast
if min(performance.contrasts(performance.use_contrasts)) <= 0.125 && ...
        performance.trialsToZeroContrast > 0
    % Subtract one from the countdown
    performance.trialsToZeroContrast = performance.trialsToZeroContrast-1;
    % If at zero, add the 0 contrast condition
    if performance.trialsToZeroContrast == 0
        performance.use_contrasts(performance.contrasts == 0) = true;
    end    
end

%%%% Pick next contrast (and update staircase on staircase trials)
staircaseTrial = performance.staircase(4) == 0;

if ~staircaseTrial
    
    % Next contrast is random from current contrast set
    performance.trialContrast = randsample(performance.contrasts(performance.use_contrasts),1);
    
elseif staircaseTrial
        
    % Update hit/miss counter
    performance.staircase(2) = performance.staircase(2) + performance.hit;
    performance.staircase(3) = performance.staircase(3) + performance.miss;
    
    % Update trial counter
    performance.staircase(4) = performance.staircase(4) + 1;
    if performance.staircase(4) >= performance.staircase(5)
        performance.staircase(4) = 0;
    end
    
    % Move staircase on hit/miss counter threshold
    if performance.staircase(2) >= performance.staircase(6)
        % On hit threshold, move the staircase forward and reset hit/miss
        newStaircaseContrast = performance.contrasts(...
            min(find(performance.staircase(1) == performance.contrasts)+1, ...
            sum(performance.use_contrasts)));
        performance.staircase(1) = newStaircaseContrast;
        performance.staircase(2:3) = 0;
    elseif performance.staircase(3) >= performance.staircase(7)
        % On miss threshold, move staircase backward and reset hit/miss
        newStaircaseContrast = performance.contrasts(...
            max(find(performance.staircase(1) == performance.contrasts)-1,1));
        performance.staircase(1) = newStaircaseContrast;
        performance.staircase(2:3) = 0;
    end
    
    % Next contrast is defined by the staircase
    performance.trialContrast = performance.staircase(1);
    
end

%%%% Pick next side (this is done at random)
performance.trialSide = randsample([-1,1],1);

end


