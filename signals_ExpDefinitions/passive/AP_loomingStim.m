function AP_loomingStim(t, events, parameters, visStim, inputs, outputs, audio)
% Present looming stimulus on left/center/right
% (looming rectangle: vis.patch circles don't support signals radii)
% Pseudorandom (each stim presented once for each trial)
% Number of trials = number of repeats


%% Set up stimuli

stim_time = 0.5;
min_iti = 2;
max_iti = 3;
step_iti = 0.1;

% Visual stim
azimuths = [-90,0,90];
vis_params = CombVec(azimuths);


%% Set trial data

% Signals garbage: things can't happen at the exact same time as newTrial
new_trial_set = events.newTrial.delay(0);

% Start clock for trial
trial_t = t - t.at(new_trial_set);

% Set the stim order and ITIs for this trial
stimOrder = new_trial_set.map(@(x) randperm(size(vis_params,2)));
stimITIs = new_trial_set.map(@(x) randsample(min_iti:step_iti:max_iti,size(vis_params,2),true));

% Get the stim on times and the trial end time
trial_stimOn_times = stimITIs.map(@(x) [0,cumsum(x(1:end-1) + stim_time)]);
trial_end_time = stimITIs.map(@(x) sum(x) + stim_time);


%% Present stim

% Visual
stim_num = trial_t.ge(trial_stimOn_times).sum.skipRepeats;
stim_id = map2(stimOrder,stim_num,@(stim_order,stim_num) stim_order(stim_num));
stimAzimuth = stim_id.map(@(x) vis_params(1,x));

stim = vis.patch(t, 'rectangle');
stim.colour = [0,0,0];

stimOn = stim_id.to(stim_id.delay(stim_time));
stim_t = t - t.at(stimOn);

stim.dims = stim_t.map(@(x) repmat(x,1,2))*(90/stim_time);
stim.azimuth = stimAzimuth;

stim.show = stimOn;
visStim.stim = stim;

endTrial = events.newTrial.setTrigger(trial_t.gt(trial_end_time));

%% Events

events.stimITIs = stimITIs;
events.stimOn = stimOn;

events.stimAzimuth = stimAzimuth;

%%% TEMP
events.stimDims = stim.dims;

events.endTrial = endTrial;

end

















