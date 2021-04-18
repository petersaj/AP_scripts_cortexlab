function AP_stimWheelLeft(t, events, parameters, visStim, inputs, outputs, audio)
% 2021-03-09 - AP
%
% Modified from vanillaChoiceworld:
% - stim only on one side
% - 100% contrast only
% - no interactive delay
% - no audio cue
% - longer ITIs
% - variable quiescence period

%% User parameters

% Reward
rewardSize = parameters.rewardSize.at(events.expStart);

% ITI
itiMin = parameters.itiMin;
itiMax = parameters.itiMax;

% Quiescence
quiescenceMin = parameters.quiescenceMin;
quiescenceMax = parameters.quiescenceMax;


%% Fixed parameters

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
contrasts = [1];
% (which conrasts to use at the beginning of training)
startingContrasts = [true];
% (which contrasts to repeat on miss)
repeatOnMiss = [false];
% (number of trials to judge rolling performance)
trialsToBuffer = 50;
% (number of trials after introducing 12.5% contrast to introduce 0%)
trialsToZeroContrast = 500;
sigma = [20,20];
spatialFreq = 1/15;
stimFlickerFrequency = 5; % DISABLED BELOW
startingAzimuth = 90;
responseDisplacement = 90;

% Timing
cueInteractiveDelay = 0;
stimOffDelay = 1; % (time stim on screen after response)

% Sounds
% audioSampleRate = 192e3; % (Not used after new audio system 20181107)
audDev = audio.Devices('Strix');
audioSampleRate = audDev.DefaultSampleRate;
audioChannels = audDev.NrOutputChannels;

onsetToneAmplitude = 0; % Changed from 1 when moving to Strix
onsetToneFreq = 12000;
onsetToneDuration = 0.1;
onsetToneRampDuration = 0.01;
toneSamples = onsetToneAmplitude*events.expStart.map(@(x) ...
    aud.pureTone(onsetToneFreq,onsetToneDuration,audioSampleRate, ...
    onsetToneRampDuration,audioChannels));

missNoiseDuration = 0.5;
missNoiseAmplitude = 0.03;
missNoiseSamples = missNoiseAmplitude*events.expStart.map(@(x) ...
    randn(audioChannels,audioSampleRate*missNoiseDuration));

% Wheel parameters
quiescThreshold = 1;
wheelGain = 8; % deg/mm

% Key press for manual reward
rewardKeyPressed = inputs.keyboard.strcmp('w');

%% Default parameters

try
    parameters.rewardSize = 6;
    parameters.itiMin = 1;
    parameters.itiMax = 3;
    parameters.quiescenceMin = 0.5;
    parameters.quiescenceMax = 1;
catch
end

%% Initialize trial data

trialDataInit = events.expStart.mapn( ...
    contrasts,startingContrasts,repeatOnMiss, ...
    trialsToBuffer,trialsToZeroContrast,staircaseTrials,staircaseHit,staircaseMiss, ...
    @initializeTrialData).subscriptable;

%% Set up wheel 

wheel = inputs.wheelMM.skipRepeats();

%% Trial event times
% (this is set up to be independent of trial conditon, that way the trial
% condition can be chosen in a performance-dependent manner)

% Resetting pre-stim quiescent period
trialQuiescence = events.newTrial.mapn(quiescenceMin,quiescenceMax,@chooseQuiescence);
preStimQuiescence = sig.quiescenceWatch(trialQuiescence,t,wheel,quiescThreshold); 

% Stimulus onset
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
trialITI = events.newTrial.mapn(itiMin,itiMax,@chooseITI);
iti = response.delay(trialITI);

% Stim stays on for delay period
stimOff = response.delay(stimOffDelay);

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

% Trial times (set by trial)
events.trialQuiescence = trialQuiescence;
events.trialITI = trialITI;

% Trial times (fixed)
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
events.staircase = trialData.staircase;
events.useContrasts = trialData.useContrasts;
events.trialsToZeroContrast = trialData.trialsToZeroContrast;
events.hitBuffer = trialData.hitBuffer;
events.sessionPerformance = trialData.sessionPerformance;
events.totalWater = totalWater;

end

function trialDataInit = initializeTrialData(subject_info, ...
    contrasts,startingContrasts,repeatOnMiss,trialsToBuffer, ...
    trialsToZeroContrast,staircaseTrials,staircaseHit,staircaseMiss)

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
% Set the first trial side (always left)
trialDataInit.trialSide = -1;
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
% Initialize the staircase: 
% [current contrast, hits, misses, staircase trial counter, 
% staircase every n trials, hit requirement, miss requirement]
trialDataInit.staircase = ...
    [contrasts(1),0,0,0, ...
    staircaseTrials,staircaseHit,staircaseMiss];
% Initialize the day's performance to plot
% [conditions, number of trials, number of move left choices]
trialDataInit.sessionPerformance = ...
    [sort(unique([-contrasts,contrasts])); ...
    zeros(size(unique([-contrasts,contrasts]))); ...
    zeros(size(unique([-contrasts,contrasts])))];

%%%% Load the last experiment for the subject if it exists
% (note: MC creates folder on initilization, so start search at 1-back)
expRef = dat.listExps(subject);
useOldParams = false;
if length(expRef) > 1
    % Get name of current expDef
    [~,curr_expDef,~] = fileparts(mfilename);

    % Loop through blocks from latest to oldest, if any have the relevant
    % parameters then carry them over
    for check_expt = length(expRef)-1:-1:1
        previousBlockFilename = dat.expFilePath(expRef{check_expt}, 'block', 'master');
        if exist(previousBlockFilename,'file')
            previousBlock = load(previousBlockFilename);
        else
            continue
        end
        % Check if it was the same expDef and fields are populated
        if isfield(previousBlock.block,'expDef')
            [~,old_expDef,~] = fileparts(previousBlock.block.expDef);
        else 
            continue
        end
        if strcmp(curr_expDef,old_expDef) && ~( ...
                isempty(previousBlock.block.events.useContrastsValues) || ...
                isempty(previousBlock.block.events.hitBufferValues) || ...
                isempty(previousBlock.block.events.trialsToZeroContrastValues))
            % Break the loop and use these parameters
            useOldParams = true;
            break
        end
    end
end

if useOldParams
    % If the last experiment file has the relevant fields, set up performance
    
    % Which contrasts are currently in use
    trialDataInit.useContrasts = previousBlock.block. ...
        events.useContrastsValues(end-length(contrasts)+1:end);
    
    % The buffer to judge recent performance for adding contrasts
    trialDataInit.hitBuffer = ...
        previousBlock.block. ...
        events.hitBufferValues(:,end-length(contrasts)+1:end);
    
    % The countdown to adding 0% contrast
    trialDataInit.trialsToZeroContrast = previousBlock.block. ...
        events.trialsToZeroContrastValues(end);
    
else
    % If this animal has no previous experiments, initialize performance
    trialDataInit.useContrasts = startingContrasts;
    trialDataInit.hitBuffer = nan(trialsToBuffer,length(contrasts));
    trialDataInit.trialsToZeroContrast = trialsToZeroContrast;  
end

end

function trialData = updateTrialData(trialData,stimDisplacement)
% Update the performance and pick the next contrast

%%%% Get index of current trial contrast
currentContrastIdx = trialData.trialContrast == trialData.contrasts;

%%%% Define response type based on trial condition
trialData.hit = stimDisplacement*trialData.trialSide < 0;
trialData.miss = stimDisplacement*trialData.trialSide > 0;

%%%% Update buffers and counters if not a repeat trial
if ~trialData.repeatTrial
    
    %%%% Contrast-adding performance buffer
    % Update hit buffer for running performance
    trialData.hitBuffer(:,currentContrastIdx) = ...
        [trialData.hit;trialData.hitBuffer(1:end-1,currentContrastIdx)];
    
    %%%% Staircase
    % Update staircase trial counter
    trialData.staircase(4) = trialData.staircase(4) + 1;
    if trialData.staircase(4) >= trialData.staircase(5)
        trialData.staircase(4) = 0;
    end
    
    % Update hit/miss counter
    trialData.staircase(2) = trialData.staircase(2) + trialData.hit;
    trialData.staircase(3) = trialData.staircase(3) + trialData.miss;
        
    % Move staircase on hit/miss counter threshold
    if trialData.staircase(2) >= trialData.staircase(6)
        % On hit threshold, move the staircase forward and reset hit/miss
        newStaircaseContrast = trialData.contrasts(...
            min(find(trialData.staircase(1) == trialData.contrasts)+1, ...
            sum(trialData.useContrasts)));
        trialData.staircase(1) = newStaircaseContrast;
        trialData.staircase(2:3) = 0;
    elseif trialData.staircase(3) >= trialData.staircase(7)
        % On miss threshold, move staircase backward and reset hit/miss
        newStaircaseContrast = trialData.contrasts(...
            max(find(trialData.staircase(1) == trialData.contrasts)-1,1));
        trialData.staircase(1) = newStaircaseContrast;
        trialData.staircase(2:3) = 0;
    end
    
    %%% Session performance for plotting
    currCondition = trialData.trialContrast*trialData.trialSide;
    currConditionIdx = trialData.sessionPerformance(1,:) == currCondition;
    trialData.sessionPerformance(2,currConditionIdx) = ...
        trialData.sessionPerformance(2,currConditionIdx) + 1;
    trialData.sessionPerformance(3,currConditionIdx) = ...
        trialData.sessionPerformance(3,currConditionIdx) + (stimDisplacement < 0);
    
end

%%%% Add new contrasts as necessary given performance
% This is based on the last trialsToBuffer trials for rolling performance
% (these parameters are hard-coded because too specific)
% (these are side-independent)
current_min_contrast = min(trialData.contrasts(trialData.useContrasts & trialData.contrasts ~= 0));
trialsToBuffer = size(trialData.hitBuffer,1);
switch current_min_contrast
    
    case 0.5
        % Lower from 0.5 contrast after > 70% correct
        min_hit_percentage = 0.70;
        
        contrast_buffer_idx = ismember(trialData.contrasts,[0.5,1]);
        contrast_total_trials = sum(sum(~isnan(trialData.hitBuffer(:,contrast_buffer_idx))));
        % If there have been enough buffer trials, check performance
        if contrast_total_trials >= size(trialData.hitBuffer,1)
            % Sample as evenly as possible across pooled contrasts
            pooled_hits = reshape(trialData.hitBuffer(:,contrast_buffer_idx)',[],1);
            use_hits = sum(pooled_hits(find(~isnan(pooled_hits),trialsToBuffer)));
            min_hits = find(1 - binocdf(1:trialsToBuffer,trialsToBuffer,min_hit_percentage) < 0.05,1);
            if use_hits >= min_hits
                trialData.useContrasts(find(~trialData.useContrasts,1)) = true;
            end
        end

    case 0.25
        % Lower from 0.25 contrast after > 50% correct
        min_hit_percentage = 0.5;
        
        contrast_buffer_idx = ismember(trialData.contrasts,current_min_contrast);
        contrast_total_trials = sum(sum(~isnan(trialData.hitBuffer(:,contrast_buffer_idx))));
        % If there have been enough buffer trials, check performance
        if contrast_total_trials >= size(trialData.hitBuffer,1)
            % Sample as evenly as possible across pooled contrasts
            pooled_hits = reshape(trialData.hitBuffer(:,contrast_buffer_idx)',[],1);
            use_hits = sum(pooled_hits(find(~isnan(pooled_hits),trialsToBuffer)));
            min_hits = find(1 - binocdf(1:trialsToBuffer,trialsToBuffer,min_hit_percentage) < 0.05,1);
            if use_hits >= min_hits
                trialData.useContrasts(find(~trialData.useContrasts,1)) = true;
            end
        end
        
    case 0.125
        % Lower from 0.25 contrast after > 50% correct
        min_hit_percentage = 0.5;
        
        contrast_buffer_idx = ismember(trialData.contrasts,current_min_contrast);
        contrast_total_trials = sum(sum(~isnan(trialData.hitBuffer(:,contrast_buffer_idx))));
        % If there have been enough buffer trials, check performance
        if contrast_total_trials >= size(trialData.hitBuffer,1)
            % Sample as evenly as possible across pooled contrasts
            pooled_hits = reshape(trialData.hitBuffer(:,contrast_buffer_idx)',[],1);
            use_hits = sum(pooled_hits(find(~isnan(pooled_hits),trialsToBuffer)));
            min_hits = find(1 - binocdf(1:trialsToBuffer,trialsToBuffer,min_hit_percentage) < 0.05,1);
            if use_hits >= min_hits
                trialData.useContrasts(find(~trialData.useContrasts,1)) = true;
            end
        end          
        
end

% Add 0 contrast after trialsToZeroContrast trials with 0.125 contrast
if min(trialData.contrasts(trialData.useContrasts)) <= 0.125 && ...
        trialData.trialsToZeroContrast > 0
    % Subtract one from the countdown
    trialData.trialsToZeroContrast = trialData.trialsToZeroContrast-1;
    % If at zero, add 0 contrast
    if trialData.trialsToZeroContrast == 0
        trialData.useContrasts(trialData.contrasts == 0) = true;
    end    
end

%%%% Set flag to repeat - skip trial choice if so
if trialData.miss && ...
        ismember(trialData.trialContrast,trialData.contrasts(trialData.repeatOnMiss))
    trialData.repeatTrial = true;
    return
else
    trialData.repeatTrial = false;
end

%%%% Pick next contrast

% Define whether this is a staircase trial
staircaseTrial = trialData.staircase(4) == 0;

if ~staircaseTrial
    % Next contrast is random from current contrast set
    trialData.trialContrast = randsample(trialData.contrasts(trialData.useContrasts),1);    
elseif staircaseTrial  
    % Next contrast is defined by the staircase
    trialData.trialContrast = trialData.staircase(1);    
end

%%%% Pick next side (always left side)
trialData.trialSide = -1;

end

function trialQuiescence = chooseQuiescence(~,quiescenceMin,quiescenceMax)
quiescence_interval = 0.1;
trialQuiescence = randsample(quiescenceMin:quiescence_interval:quiescenceMax,1);
end

function trialITI = chooseITI(~,itiMin,itiMax)
iti_interval = 0.1;
trialITI = randsample(itiMin:iti_interval:itiMax,1);
end



