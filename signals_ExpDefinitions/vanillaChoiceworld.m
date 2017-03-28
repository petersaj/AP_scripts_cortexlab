function vanillaChoiceworld(t, events, parameters, visStim, inputs, outputs, audio)
% vanillaChoiceworld(t, events, parameters, visStim, inputs, outputs, audio)
% Choice world that adapts with behavior
% 170309 - AP

%% Fixed parameters

% Reward
rewardSize = 0.2; %%%% TO DO: CHECK OTHER PROTOCOLS

% Trial choice parameters
% (staircase: every staircaseTrials, engage an staircaseHit forward / staircaseMiss back
% choice of contrast)
staircaseTrials = 2; 
staircaseHit = 3;
staircaseMiss = 1;

% Stimulus/target
contrasts = [1,0.5,0.25,0.125,0.06,0];
startingContrasts = [true,true,false,false,false,false];
trialsToBuffer = 10; % number of trials to judge current performance
trialsToZeroContrast = 500; % number of trials after introducing 0.125
sigma = [9,9];
spatialFrequency = 0.01;
startingAzimuth = 90;
hitDisplacement = 90;
missDisplacement = 90;

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
    subject,contrasts,startingContrasts,trialsToBuffer,trialsToZeroContrast,@initializePerformance).subscriptable;

% Initialize staircase
staircase = events.trialNum.mapn(@(trialNum) mod(trialNum,staircaseTrials) == 0);

% Choose trial condition
trialSide = events.newTrial.mapn(@(x) randsample([-1,1],1));

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

% Play tone at interactive onset
audio.onsetTone = toneSamples.at(interactiveOn);

% Response
% (wheel displacement zeroed at interactiveOn)
stimDisplacement = wheelGain*(wheel - wheel.at(interactiveOn));

hit = keepWhen(interactiveOn.setTrigger(stimDisplacement*trialSide ...
    <= -hitDisplacement),interactiveOn.to(events.newTrial));
miss = keepWhen(interactiveOn.setTrigger(stimDisplacement*trialSide ...
    >= missDisplacement),interactiveOn.to(events.newTrial));

% Give reward on hit
outputs.reward = at(rewardSize,hit);  

% Play noise on miss
audio.missNoise = missNoiseSamples.at(miss);

% ITI defined by outcome
iti = merge(hit.delay(itiHit),miss.delay(itiMiss));

% Stim fixed in place before interactive and after response, wheel-conditional otherwise
stimAzimuth = cond( ...
    events.newTrial.to(interactiveOn), startingAzimuth*trialSide, ...
    interactiveOn.to(merge(hit,miss)), startingAzimuth*trialSide + stimDisplacement, ...
    hit.to(events.newTrial), (startingAzimuth-hitDisplacement)*trialSide, ...
    miss.to(events.newTrial), (startingAzimuth+missDisplacement)*trialSide);

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
    [trialSide.at(merge(hit,miss)), ...
    outcome.at(merge(hit,miss)), ...
    staircase.at(merge(hit,miss)), ...
    staircaseHit,staircaseMiss];

performance = trialInfo.scan(@updatePerformance,performanceInit).subscriptable;
trialContrast = performance.nextTrialContrast;

%% Visual stimulus

% Trial-independent parameters
stim = vis.grating(t, 'square', 'gaussian');
stim.sigma = sigma;
stim.spatialFrequency = spatialFrequency;
stim.phase = 2*pi*events.newTrial.map(@(v)rand);
stim.azimuth = stimAzimuth;
stim.contrast = trialContrast.at(stimOn);
stim.show = stimOn.to(stimOff);

visStim.stim = stim;

%% Display and save

% Trial conditions and outcome
events.trialSide = trialSide;
events.trialContrast = trialContrast;

% Animal-dependent events
events.wheel = wheel;
events.hit = hit;
events.miss = miss;

% Trial times
events.stimOn = stimOn;
events.stimOff = stimOff;
events.interactiveOn = interactiveOn;
events.stimAzimuth = stimAzimuth;
events.endTrial = endTrial;

% Performance
events.use_contrasts = performance.use_contrasts;
events.hit_buffer = performance.hit_buffer;
events.trialsToZeroContrast = performance.trialsToZeroContrast;

end

function performanceInit = initializePerformance(~,subject,contrasts,startingContrasts,trialsToBuffer,trialsToZeroContrast)

performanceInit = struct;
performanceInit.contrasts = contrasts;
performanceInit.conditions = unique(sort([contrasts,-contrasts]));
performanceInit.nextTrialContrast = 1;
% Initialize the staircase: [current contrast, hits, misses]
performanceInit.staircase = [contrasts(1),0,0];

n_conditions = length(performanceInit.conditions);

% Load the last experiment for the subject if it exists
% (note: MC creates folder on initilization, so look for > 1)
expRef = dat.listExps(subject);
if length(expRef) > 1
    previousBlockFilename = dat.expFilePath(expRef{end-1}, 'block', 'master');
    previousBlock = load(previousBlockFilename);
end

if exist('previousBlock','var') && all(isfield(previousBlock.block, ...
        {'use_contrasts','hit_buffer','trialsToZeroContrast'}))
    % If the last experiment file has the relevant fields, set up performance   
    performanceInit.use_contrasts = previousBlock.block.use_contrasts(end-n_conditions+1:end);
    performanceInit.hit_buffer = reshape( ...
        previousBlock.block.hit_buffer(end-trialsToBuffer*n_conditions+1:end), ...
        trialsToBuffer,n_conditions);
    performanceInit.trialsToZeroContrast = previousBlock.block.trialsToZeroContrast(end);        
else    
    % If this animal has no previous experiments, initialize performance
    performanceInit.use_contrasts = startingContrasts;
    performanceInit.hit_buffer = nan(trialsToBuffer,n_conditions);
    performanceInit.trialsToZeroContrast = trialsToZeroContrast;  
end

end

function performance = updatePerformance(performance,trialInfo)
% Update the performance and pick the next contrast

%%%% Unpackage (have to be packaged: only one allowable input argument)
trialSide = trialInfo(1);
hit = trialInfo(2);
staircaseTrial = trialInfo(3);
staircaseHit = trialInfo(4);
staircaseMiss = trialInfo(5);

thisTrialCondition = trialSide*performance.nextTrialContrast;

%%%% Update performance
conditionIdx = performance.conditions == thisTrialCondition;

% Hit/miss in buffer
performance.hit_buffer(:,conditionIdx) = ...
    [hit;performance.hit_buffer(1:end-1,conditionIdx)];
trialsToBuffer = size(performance.hit_buffer,1);

%%%% Add new contrasts as necessary given performance
% This is based on the last trialsToBuffer trials for rolling performance
% (these parameters are hard-coded because too specific)
% (these are side-independent)
current_min_contrast = min(performance.contrasts(performance.use_contrasts & performance.contrasts ~= 0));
switch current_min_contrast
    
    case 0.5
        % Lower from 0.5 contrast after > 75% correct
        min_hit_percentage = 0.75;
        
        curr_condition = ismember(abs(performance.conditions),[current_min_contrast,1]);
        condition_total_trials = sum(sum(~isnan(performance.hit_buffer(:,curr_condition))));
        % If there have been enough buffer trials, check performance
        if condition_total_trials >= size(performance.hit_buffer,1)
            % Sample as evenly as possible across pooled conditions
            pooled_hits = reshape(performance.hit_buffer(:,curr_condition)',[],1);
            use_hits = sum(pooled_hits(find(~isnan(pooled_hits),trialsToBuffer)));
            min_hits = find(1 - binocdf(1:trialsToBuffer,trialsToBuffer,min_hit_percentage) < 0.05,1);
            if use_hits > min_hits
                performance.use_contrasts(find(~performance.use_contrasts,1)) = true;
            end
        end

    case 0.25
        % Lower from 0.25 contrast after > 50% correct
        min_hit_percentage = 0.5;
        
        curr_condition = ismember(abs(performance.conditions),[current_min_contrast,1]);
        condition_total_trials = sum(sum(~isnan(performance.hit_buffer(:,curr_condition))));
        % If there have been enough buffer trials, check performance
        if condition_total_trials >= size(performance.hit_buffer,1)
            % Sample as evenly as possible across pooled conditions
            pooled_hits = reshape(performance.hit_buffer(:,curr_condition)',[],1);
            use_hits = sum(pooled_hits(find(~isnan(pooled_hits),trialsToBuffer)));
            min_hits = find(1 - binocdf(1:trialsToBuffer,trialsToBuffer,min_hit_percentage) < 0.05,1);
            if use_hits > min_hits
                performance.use_contrasts(find(~performance.use_contrasts,1)) = true;
            end
        end
        
    case 0.125
        % Lower from 0.25 contrast after > 50% correct
        min_hit_percentage = 0.5;
        
        curr_condition = ismember(abs(performance.conditions),[current_min_contrast,1]);
        condition_total_trials = sum(sum(~isnan(performance.hit_buffer(:,curr_condition))));
        % If there have been enough buffer trials, check performance
        if condition_total_trials >= size(performance.hit_buffer,1)
            % Sample as evenly as possible across pooled conditions
            pooled_hits = reshape(performance.hit_buffer(:,curr_condition)',[],1);
            use_hits = sum(pooled_hits(find(~isnan(pooled_hits),trialsToBuffer)));
            min_hits = find(1 - binocdf(1:trialsToBuffer,trialsToBuffer,min_hit_percentage) < 0.05,1);
            if use_hits > min_hits
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
if ~staircaseTrial
    
    % Next contrast is random from current contrast set
    performance.nextTrialContrast = randsample(performance.contrasts(performance.use_contrasts),1);
    
elseif staircaseTrial
    
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
    
    % Next contrast is defined by the staircase
    performance.nextTrialContrast = performance.staircase(1);
    
end


end


