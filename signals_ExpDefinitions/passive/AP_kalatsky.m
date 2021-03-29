function AP_kalatsky(t, events, parameters, visStim, inputs, outputs, audio)

% Kalatsky-style sweeping bar stimulus for retinotopic mapping
% (modified from example in the signals docs)

%% Timing
iti = 2; % time between periodic stimuli

sweepFreq = 0.1;% sweep frequency in sweeps/sec
n_sweeps = 10; % sweeps per stim

%% Trial structure

stimTime = n_sweeps/sweepFreq;

% Alternate vertical/horizontal every 2 trials, direction every 4 trials
% (NOTE: to have equal numbers of all, total trial number must be x4)
vertHorzTrial = mod(events.trialNum,2);
forwardBackwardTrial = 1-mod(floor((events.trialNum+1)/2),2);

% Turn stim on for fixed time based on repeats and frequency
stimOn_vert = at(true,vertHorzTrial); 
stimOff_vert = stimOn_vert.delay(stimTime);

stimOn_horz = at(true,~vertHorzTrial); 
stimOff_horz = stimOn_horz.delay(stimTime);

endTrial = merge(stimOff_vert.delay(iti),stimOff_horz.delay(iti));

%% Vertical stimuli

% Define some parameters used by both bars.  These are constants but could
% be made into parameter signals or be derived from other signals
aziRange = [135,0]; % The horizontal sweep bounds of the bars
w = 10; % bar width in visual degrees
flipFreq = 2; % Frequency of colour reversal in Hz

% Compute colour flips and position
white = skipRepeats(t.mod(1) < 1/flipFreq); % periodically flips between true and false

trial_t = t - t.at(events.newTrial);
azimuth = range(aziRange)*forwardBackwardTrial - ...
    (aziRange(1) - mod(trial_t*(range(aziRange)*sweepFreq),range(aziRange)));

% Define the two bars (left hemifield)
barLeft1 = vis.patch(t, 'rect');
% Current azimuth shifted left by half the total width
barLeft1.azimuth = -(azimuth - w/2);
barLeft1.colour = [1 1 1] * white;
barLeft1.dims = [w 150]; % [azimuth altitude]

barLeft2 = vis.patch(t, 'rect');
% Current azimuth shifted right by half the total width
barLeft2.azimuth = -(azimuth + w/2);
barLeft2.colour = [1 1 1] * ~white;
barLeft2.dims = [w 150]; % [azimuth altitude]

% Define two bars (right hemifield)
barRight1 = vis.patch(t, 'rect');
barRight1.azimuth = azimuth - w/2;
barRight1.colour = [1 1 1] * white;
barRight1.dims = [w 150]; % [azimuth altitude]

barRight2 = vis.patch(t, 'rect');
barRight2.azimuth = azimuth + w/2;
barRight2.colour = [1 1 1] * ~white;
barRight2.dims = [w 150]; % [azimuth altitude]

% Show both bar stimuli
barLeft1.show = stimOn_vert.to(stimOff_vert);
barLeft2.show = stimOn_vert.to(stimOff_vert);
barRight1.show = stimOn_vert.to(stimOff_vert);
barRight2.show = stimOn_vert.to(stimOff_vert);

% Render stimulus
visStim.barLeft1 = barLeft1;
visStim.barLeft2 = barLeft2;
visStim.barRight1 = barRight1;
visStim.barRight2 = barRight2;


%% Horizontal stimuli

% Define some parameters used by both bars.  These are constants but could
% be made into parameter signals or be derived from other signals
altRange = [-45,45]; % The horizontal sweep bounds of the bars
w = 10; % bar width in visual degrees
flipFreq = 2; % Frequency of colour reversal in Hz

% Compute colour flips and position
white = skipRepeats(t.mod(1) < 1/flipFreq); % periodically flips between true and false

trial_t = t - t.at(events.newTrial);
altitude = sign(0.5-forwardBackwardTrial)* ...
    (altRange(1) + mod(trial_t*(range(altRange)*sweepFreq),range(altRange)));

% Define the two bars (left hemifield)
barTop = vis.patch(t, 'rect');
barTop.altitude = altitude - w/2;
barTop.colour = [1 1 1] * white;
barTop.dims = [270 w]; % [azimuth altitude]

barBottom = vis.patch(t, 'rect');
barBottom.altitude = altitude + w/2;
barBottom.colour = [1 1 1] * ~white;
barBottom.dims = [270 w]; % [azimuth altitude]

% Show both bar stimuli
barTop.show = stimOn_horz.to(stimOff_horz);
barBottom.show = stimOn_horz.to(stimOff_horz);

% Render stimulus
visStim.barTop = barTop;
visStim.barBottom = barBottom;


%% TEST

% % Create a sampler at 60 Hz
% samplerFs = 60; % Hz
% sampler = skipRepeats(floor(t*samplerFs)); % updates at our sampling rate
% 
% % Create the checker stimulus
% noise = vis.checker(t);
% noise.show = true; % Turn on
% noise.colour = [1 1 1]; % Normalized RGB values
% 
% % A grid of 60 rectangles along the altitude, 60 along the azimuth
% gridSize = [60 60];
% % Each time the sampler updates, call randi, which will produce a [10 20]
% % array of values between 0 and 3.  Subtract 2 from this to make values
% % between -1 and 1.
% noise.pattern = sampler.map(@(~) randi(3, gridSize)-2);
% % Each rectangle takes up all of its space, leaving no gap between
% % neighbouring rectangles.
% noise.rectSizeFrac = [1 1];
% % Set the range of the checkerboard in visual degrees.  The pulasation in
% % each dimension will happen in anti-phase and the maximum range will be
% % [-60 60] along the azimuth, and [-30 30] along the altitude:
% noise.azimuthRange = (0.75 + 0.25*cos(2*pi*t))*[-60 60];
% noise.altitudeRange = (0.75 + 0.25*sin(2*pi*t))*[-30 30];
% % Assign to our visual stimulus object for rendering
% visStim.noise = noise;

%% Events

events.stimOn_vert = stimOn_vert;
events.stimOff_vert = stimOff_vert;
events.stimOn_horz = stimOn_horz;
events.stimOff_horz = stimOff_horz;

events.azimuth = azimuth;
events.altitude = altitude;

events.endTrial = endTrial;






