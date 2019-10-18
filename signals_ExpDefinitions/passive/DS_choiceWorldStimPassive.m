function DS_choiceWorldStimPassive(t, events, parameters, visStim, inputs, outputs, audio)
% AP 2017-12-13: present choiceworld stimuli passively
error('Zero stim doesn''t work')

%% Set up stimuli

numRepeats = 30;
stimTime = 0.5;
minITI = 2;
maxITI = 4;
stepITI = 1;
prestimQuiescentTime = 0.5;

% Visual stim
sigma = [8,8];
azimuths = [-80,80];
contrasts = [1,0.5,0.25,0.12,0];
vis_params = mat2cell(CombVec(azimuths,contrasts),2,ones(1,length(azimuths)*length(contrasts)));
spatialFreq = 0.1;

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
stim = vis.grating(t, 'square', 'gaussian');
stim.spatialFreq = spatialFreq;
stim.sigma = sigma;
stim.phase = 2*pi*events.newTrial.map(@(v)rand);

% Get trial stim features
stim_azimuth = events.trialNum.map(@(x) visual_stim_shuffle{x}(1));
stim_contrast = events.trialNum.map(@(x) visual_stim_shuffle{x}(2));

stim.azimuth = stim_azimuth.at(stim_azimuth);
stim.contrast = stim_contrast.at(stim_contrast);

stim.show = stimOn.to(stimOff);
visStim.stim = stim;

% End the trial after stim onset + stim time + trial ITI
endTrial = stimOn.delay(stimTime + events.trialNum.map(@(x) visual_itiTimes(x)));

%% Events

events.visualParams = events.expStart.map(@(x) visual_stim_shuffle);
events.stimOn = stimOn;
events.stimOff = stimOff;
events.stimAzimuth = stim_azimuth;
events.stimContrast = stim_contrast;

events.endTrial = endTrial;

end

















