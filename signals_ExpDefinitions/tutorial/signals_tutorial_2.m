function signals_tutorial_basic(t, events, pars, visStim, inputs, outputs, audio)
% Signals tutorial 2
%
% TO DO HERE
% give the order of later tutorials:
% - stimMoveTrigger (includes keepWhen and setTrigger)?
% - scan
%
%% Trial contingencies %%

% -- UNCOMMENT --
wheel = inputs.wheel.skipRepeats;
stim = vis.grating(t, 'sinusoid', 'gaussian');

stim_azimuth = wheel - wheel.at(
stim.azimuth = cond( ...
    events.expStart,wheel, ...
    true, 0);

stim.show = true;
visStim.gaborStim = stim;
% ---------------
