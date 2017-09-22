function AP_brightLever_test(t, events, parameters, visStim, inputs, outputs, audio)
% AP 2017-06-23: pair visual stimuli with rewards

%% Fixed parameters

% Trial params
iti = 5;
colorRate = 1; % number of seconds to reach maximum brightness

% Lever
lever_r_flip = inputs.lever_r.delta.skipRepeats;

% Stimulus
stimOnset = events.newTrial;
stimOffset = stimOnset.setTrigger(lever_r_flip.ge(1));

rectGray = skipRepeats(map((t - t.at(stimOnset))/colorRate,@(x) min(x,1)));
rectColor = [rectGray,rectGray,rectGray];

rectStim = vis.patch(t, 'rect');
rectStim.azimuth = 90;
rectStim.colour = rectColor;
rectStim.show = stimOnset.to(stimOffset);
rectStim.dims = [90,90];

visStim.rectStim = rectStim;

% Saved events
events.stimOnset = stimOnset;
events.stimOffset = stimOffset;
events.rectGray = rectGray;
events.endTrial = stimOffset.delay(iti);
























