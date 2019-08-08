function AP_lcrGratingPassive(t, events, parameters, visStim, inputs, outputs, audio)
% Present static grating (as in choiceworld) passively left/center/right 100% contrast
% NOTE: SET ONLY ONE TRIAL IN MC, RUNS THROUGH ALL STIM

%% Set up stimuli

staticParameters.numRepeats = 100;
staticParameters.stimTime = 0.5;
staticParameters.minITI = 2;
staticParameters.maxITI = 3;
staticParameters.stepITI = 0.1;
staticParameters.bufferTime = 10; % time to delay experiment start and end to wait for timeline etc

% Visual stim
sigma = [20,20];
azimuths = [-90,0,90];
contrasts = [1];
vis_params = mat2cell(CombVec(azimuths,contrasts),2,ones(1,length(azimuths)*length(contrasts)));
spatialFreq = 1/15;


%% Stim times

% Visual
visual_stim_uniform = repmat(vis_params,staticParameters.numRepeats,1);
visual_stim_shuffle = visual_stim_uniform(randperm(numel(visual_stim_uniform)));

visual_itiTimes = randsample(staticParameters.minITI:staticParameters.stepITI:staticParameters.maxITI,length(visual_stim_shuffle),true);
visual_startTimes = staticParameters.bufferTime + cumsum(visual_itiTimes + staticParameters.stimTime);


%% Present stim

% Visual
visualOnset = t.map(@(t) sum(t > visual_startTimes)).skipRepeats;

stim = vis.grating(t, 'square', 'gaussian');
stim.spatialFreq = spatialFreq;
stim.sigma = sigma;
stim.phase = 2*pi*events.newTrial.map(@(v)rand);

stim.azimuth = visualOnset.at(visualOnset).map(@(x) visual_stim_shuffle{x}(1));
stim.contrast = visualOnset.at(visualOnset).map(@(x) visual_stim_shuffle{x}(2));

stim.show = visualOnset.to(visualOnset.delay(staticParameters.stimTime));
visStim.stim = stim;

endTrial = t.ge(staticParameters.bufferTime + visual_startTimes(end)).skipRepeats;

%% Events

events.visualOnset = visualOnset;
events.visualParams = events.expStart.map(@(x) visual_stim_shuffle);

events.endTrial = endTrial;

end

















