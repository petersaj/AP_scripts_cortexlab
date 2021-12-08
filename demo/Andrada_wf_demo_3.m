%% Load example dataset 
% (for the next things, let's use the same example dataset as before)

% Here's an example animal/day
animal = 'AP025';
day = '2017-09-28';
experiment = 1;
verbose = true;

AP_load_experiment;

%% Exploring data

% Here are two functions which could be useful for exploring data.

% Pixel correlation viewer: this plots the correlation of a particular
% pixel with all other pixels. This doesn't tell you much about an
% experiment, but I think it's useful for getting an intuition for brain
% areas. 
pixelCorrelationViewerSVD(Udf,fVdf);

% You can now click on any pixel, and it will show you the correlation
% between fluorescence in that pixel and all other pixels (yellow areas are
% highly correlated, black is 0 correlation, blue is negative
% correlation). Also - you can click on the figure and hit 'h' which puts
% it into 'hover mode', so now you can just move your mouse around to
% change which pixel is the reference pixel.

% EXERCISE: align the U's, deconvolve the V's, then feed those into the
% pixel correlation viewer. Then overlay the atlas regions on the image and
% look around at the data:
%
% - Based on correlated areas, if you were to divide the brain into a few
% correlated chunks, where would the be? What areas/functions do those
% relate to?
%
% - Sometimes you'll see a pattern where you'll have separated spots of
% correlated pixels that merge together as you move the mouse. Where are
% these? Why do they show up where they do?

% This next tool is a little more useful for exploring data as it relates
% to the experiment: it plots the face camera, eye camera, and widefield
% images (after you open, you can scroll with the scroll bar or the scroll
% wheel on the mouse).
AP_wfmovies(Udf,fVdf,frame_t,eyecam_fn,eyecam_t,facecam_fn,facecam_t)

% EXERCISE: look through the data and come up with some observations about
% the relationship between behavior and activity. For example: what areas
% are active while the mouse is moving? What's going on in the brain when
% the mouse isn't moving? Are there types of activity that you see often?
% Are there types of activity you see rarely? 
%
% Then: run that tool using the deconvolved V's instead. You'll probably
% need to adjust the color balance as we did before - click on the
% widefield image, and in the command line type caxis(caxis/2) and repeat
% that command until you can see activity well. What's the difference
% between this and the non-deconvolved data? Can you make any observations
% about activity/behavior in this data that you couldn't see in the raw
% data?

%% Wheel movement

% The last component of our experiments is the wheel that the mice can
% turn, which is recorded into Timeline.

% EXERCISE: the wheel is recorded into Timeline under the name
% 'rotaryEncoder'. Plot the raw wheel trace in Timeline with time in
% seconds on the x-axis. This is the cumulative wheel position over the
% course of the whole experiment (i.e. each point is relative to the very
% start of the experiment, which is set at 0). Positive values means the
% wheel is turning right, negative values means the wheel is turning left.

% EXERCISE: using the interp1 function, get and plot the wheel position
% aligned to each stimulus side/contrast. Note - you'll want "0" in your
% plot to be the wheel position before the stimulus comes on, so baseline
% subtraction will be important here to make sense of the data. Do these
% traces make sense given what the task is?

% EXERCISE: instead of using the wheel position, we can also look at wheel
% velocity. Plot the derivative of the wheel position using the 'diff'
% function, and zoom into the plot: you should be able to see discrete
% positive or negative pulses, which is what the wheel actually sends out
% to the computer.

% Those blips aren't super useful to look at, so I have a function to
% smooth these and break them into "movement" and "quiescence" periods.
% This is done in the loading script, so it's already available. Let's plot
% this processed velocity ('wheel_velocity'):
t = Timeline.rawDAQTimestamps;
figure; hold on;
plot(t,wheel_velocity,'k');
% The movement/quiescence periods are stored in 'wheel_move': this is 0
% when the mouse is quiescent and 1 when the mouse is moving. Let's create
% a new trace which is the wheel velocity during movement and NaN during
% quiescence. NaNs (not-a-number) are placeholder values. In this case
% they're useful because NaNs aren't plotted, so we plot the wheel velocity
% during movement in a different color.
% (copy the wheel velocity trace)
wheel_velocity_movement = wheel_velocity;
% (whenever the mouse is quiescent, set those values to NaN)
wheel_velocity_movement(~wheel_move) = NaN;
% (plot the wheel velocity during movement in red)
plot(t,wheel_velocity_movement,'r');
% Zoom into this plot to see how the trace is either black+flat or
% red+wobbly.

% EXERCISE: using the 'wheel_move' trace, pull out times of movement
% onsets. Align the widefield data to these movement onsets and plot the
% average widefield fluorescence aligned to movement. What's going on in
% this movie?

% EXERCISE: is the movement-aligned activity different depending on how
% fast and in which direction the mouse is moving the wheel? This will be
% more complicated and will need a few steps: 
% - for each movement onset, get a number summarizing that movement (maybe
% the amplitude and sign of the velocity: remember positive values means
% the wheel is turning right, and negative values means it's turning left).
% One option for this is to split the velocity trace into a cell array,
% where each cell contains a full movement, which can be done with the
% mat2cell function. 
% - after you get a value for each movement representing how fast and in
% which direction the mouse is turning, you'll want to discretize those
% values into categories. Try using the 'discretize' function. 
% - once you have discrete values of movement (e.g. very fast left turn,
% small right turn, etc.), align the widefield images to those movement
% onsets and see how the activity varies with movement speed and direction.

%% Combining stimulus, movements, and imaging

% Now you've got the general tools to work with all of the data: let's try
% using them to get an example of separating stimulus- and movement-related
% activity. 

% EXERCISE: find the "reaction time" for each stimulus as the first time
% the mouse moves the wheel following the onset of each stimulus. As a
% sanity check, plot the median reaction time vs. stimulus and make sure
% this plot looks like what you expect it to look like. 

% EXERCISE: align widefield activity to the onset of all stimuli on the
% righthand screen (use all contrasts except zero - just lump them all
% together). Split these trials by reaction times into a few groups, so
% that you'll have trials with slow/medium/fast etc. reaction times. Make an
% average stimulus-aligned activity for each one of these reaction time
% groups. What changes in the activity across these reaction time groups?


















