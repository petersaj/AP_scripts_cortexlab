function signals_vistest(t, events, pars, visStim, inputs, outputs, audio)

stimOnset = events.newTrial;
stimOffset = events.newTrial.delay(1);
stimOnOff = stimOnset.to(stimOffset);

colorRate = 1; % number of seconds to reach maximum
rectGray = map((t - t.at(stimOnset))/colorRate,@(x) min(x,1));
rectColor = [rectGray,rectGray,rectGray];

rectStim = vis.patch(t, 'rect');
rectStim.azimuth = 90;
rectStim.colour = rectColor;
rectStim.show = stimOnOff;
rectStim.dims = [90,90];

visStim.rectStim = rectStim;

events.rectGray = rectGray;
events.endTrial = stimOffset.delay(0.5);








