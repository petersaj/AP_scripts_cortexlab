function AP_lcrGratingPassive(t, events, parameters, visStim, inputs, outputs, audio)
% Present static grating (as in choiceworld) passively left/center/right 100% contrast
% Pseudorandom (each stim presented once for each trial)
% Number of trials = number of repeats

% this doesn't work
% why the actual fuck doesn't this work and why is it so hard


%% Set up stimuli

stim_time = 0.5;
min_iti = 2;
max_iti = 3;
step_iti = 0.1;

% Visual stim
sigma = [20,20];
azimuths = [-90,0,90];
contrasts = [1];
vis_params = CombVec(azimuths,contrasts);
spatialFreq = 1/15;


%% Set trial data

% Signals garbage: things can't happen at the exact same time as newTrial
new_trial_set = events.newTrial.delay(0);

% Start clock for trial
trial_t = t - t.at(new_trial_set);

% Set the stim order and ITIs for this trial
stim_order = new_trial_set.map(@(x) randperm(size(vis_params,2)));
stim_itis = new_trial_set.map(@(x) randsample(min_iti:step_iti:max_iti,size(vis_params,2),true));

% Get the stim on times and the trial end time
stimOn_times = stim_itis.map(@(x) [0,cumsum(x(1:end-1) + stim_time)]);
trial_end_time = stim_itis.map(@(x) sum(x) + stim_time);


%% Present stim

% % Visual

stim_num = trial_t.ge(stimOn_times).sum.skipRepeats;
stim_id = map2(stim_order,stim_num,@(stim_order,stim_num) stim_order(stim_num));
stim_azimuth = stim_id.map(@(x) vis_params(1,x));
stim_contrast = stim_id.map(@(x) vis_params(2,x));

stim = vis.grating(t, 'square', 'gaussian');
stim.spatialFreq = spatialFreq;
stim.sigma = sigma;
stim.phase = 2*pi*events.newTrial.map(@(x)rand);

stim.azimuth = stim_azimuth;
stim.contrast = stim_contrast;

stim.show = stim_id.to(stim_id.delay(stim_time));
visStim.stim = stim;

endTrial = events.newTrial.setTrigger(trial_t.gt(trial_end_time));

%% Events

% events.trial_t = trial_t;
% events.stim_order = stim_order;
% events.stim_itis = stim_itis;
% events.trial_end_time = trial_end_time;

events.stim_num = stim_num;
events.stim_id = stim_id;

events.stim_azimuth = stim_azimuth;
events.stim_contrast = stim_contrast;

events.endTrial = endTrial;

end

















