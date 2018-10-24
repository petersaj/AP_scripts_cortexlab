function AP_localize_choiceWorldStimPassive(t, events, parameters, visStim, inputs, outputs, audio)
% AP 2017-12-13: present choiceworld stimuli passively

%% Set up stimuli

numRepeats = 30;
stimTime = 0.5;
minITI = 2;
maxITI = 4;
stepITI = 1;
prestimQuiescentTime = 0.5;

% Visual stim
position_offsets = 10;
contrast = 1;
spatialFreq = 0.1;

sigma = [8,8];
az_el = [65,0]; % choiceworld stim (will be mirrored L/R)

stim_positions = bsxfun(@plus,az_el, ...
    [0,0; ...
    -position_offsets,0; ...
    position_offsets,0; ...
    0,position_offsets; ...
    0,-position_offsets;]);

vis_params = mat2cell(stim_positions,ones(1,size(stim_positions,1)),2);

%% Set up wheel (for pre-stim quiescence)

wheel = inputs.wheel.skipRepeats();

%% Set up stim parameters for all trials

% Visual
visual_stim_uniform = repmat(vis_params,numRepeats,1);
visual_stim_shuffle = visual_stim_uniform(randperm(numel(visual_stim_uniform)));

visual_itiTimes = randsample(minITI:stepITI:maxITI,length(visual_stim_shuffle),true);


%% Present stim

% Present stim with resetting pre-stim quiescent period
quiescThreshold = 1;
prestimQuiescentPeriod = at(prestimQuiescentTime,events.newTrial.delay(0)); 
preStimQuiescence = sig.quiescenceWatch(prestimQuiescentPeriod, t, wheel, quiescThreshold); 
stimOn = at(true,preStimQuiescence); 
stimOff = stimOn.delay(stimTime);

% Define common stim features
stim1 = vis.grating(t, 'square', 'gaussian');
stim1.spatialFreq = spatialFreq;
stim1.sigma = sigma;
stim1.phase = 2*pi*events.newTrial.map(@(v)rand);
stim1.contrast = contrast;

stim2 = vis.grating(t, 'square', 'gaussian');
stim2.spatialFreq = spatialFreq;
stim2.sigma = sigma;
stim2.phase = 2*pi*events.newTrial.map(@(v)rand);
stim1.contrast = contrast;

% Get trial stim features
stim_azimuth = events.trialNum.map(@(x) visual_stim_shuffle{x}(1));
stim_altitude = events.trialNum.map(@(x) visual_stim_shuffle{x}(2));

stim1.azimuth = stim_azimuth.at(stim_azimuth);
stim1.altitude = stim_altitude.at(stim_altitude);

stim2.azimuth = -stim_azimuth.at(stim_azimuth);
stim2.altitude = stim_altitude.at(stim_altitude);

stim1.show = stimOn.to(stimOff);
stim2.show = stimOn.to(stimOff);

visStim.stim1 = stim1;
visStim.stim2 = stim2;

% End the trial after stim onset + stim time + trial ITI
endTrial = stimOn.delay(stimTime + events.trialNum.map(@(x) visual_itiTimes(x)));

%% Events

events.visualParams = events.expStart.map(@(x) visual_stim_shuffle);
events.stimOn = stimOn;
events.stimOff = stimOff;
events.stimAzimuth = stim_azimuth;
events.stimAltitude = stim_altitude;

events.endTrial = endTrial;

end

















