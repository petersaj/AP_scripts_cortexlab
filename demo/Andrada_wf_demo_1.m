%% Repos to download

% Clone these repos: 
% https://github.com/petersaj/AP_scripts_cortexlab
% https://github.com/cortex-lab/widefield
% https://github.com/kwikteam/npy-matlab

%% Load example dataset 

% Here's an example animal/day
animal = 'AP025';
day = '2017-09-28';
experiment = 1;
verbose = true;

% This is my catch-all code for loading data
% (finds the data on the server, loads, processes)
AP_load_experiment;

%% Timeline introduction
% Timeline is our code environment for input/output

% Anything in the experiment with a signal is recorded into timeline, which
% is saved into this structure:
Timeline

% The names of the inputs into timeline are stored here: 
{Timeline.hw.inputs.name};

% The timestamps for recorded signals are stored here: 
Timeline.rawDAQTimestamps;

% The recorded data is stored here (N timestamps x N signals)
Timeline.rawDAQData;

% The main inputs we'll use are: 
% photoDiode - this detects changes on the screen (e.g. stim presentation)
% rotaryEncoder - this is the position of the mouse's wheel
% rewardEcho - this is when a reward was given
% pcoExposure - this is when the widefield camera took a picture

% Example: if we want to plot the photodiode signal, we do that with this:
% 1) figure out which input belongs to the reward
photodiode_index = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
% 2) use that to index which raw data to plot
figure;
plot(Timeline.rawDAQTimestamps,Timeline.rawDAQData(:,photodiode_index));
xlabel('Time (s)');
ylabel('Photodiode (volts)');

% EXERCISE: find the times that a reward was given by determining when the
% 'rewardEcho' signal turned on. That's been done in the load script, so
% compare your answer with 'reward_t_timeline'


%% Signals introduction
% Signals is the code environment for our experiment protocols

% Signals saves experiment information for each trial, including things
% that can't be recorded by timeline (e.g. which stimulus was shown?)

% Each of these signals has both a value and a time, we usually only care
% about the values (because more accurate timings are in Timeline). These
% are stored in signals_events.xValues and signals_events.xTimes, for
% example trial numbers are stored in signals_events.trialNumValues and
% signals_events.trialNumTimes.

% Signals most important for this task are
% The side of the stim on each trial:
signals_events.trialSideValues; % (-1 is on the left, 1 is on the right)
% The contrast of the stim on each trial
signals_events.trialContrastValues;
% Correct (hit) or incorrect (miss) for each trial:
signals_events.hitValues;

% EXERCISE: using the above signals, make a plot showing the fraction of
% correct trials for each unique stimulus (side and contrast).


%% Widefield data introduction
% Widefield data is saved in SVD format rather than in pixels because it
% takes up way less space but contains most of the information. 

% This has two parts to it:
Udf; % U's: the spatial components (Y pixels x X pixels x N components)
fVdf; % V's: the temporal components (N components x N timepoints)
frame_t; % these are the timestamps for the temporal components (seconds)

% The spatial components represent modes that all the pixels can vary in.
% They decrease in how much variance they explain, so the first one looks
% like a brain and explains most of the variance, and the last one looks
% like noise and explains very little: 
figure; 
subplot(1,2,1);
imagesc(Udf(:,:,1));
axis image off;
title('Spatial component 1')

subplot(1,2,2);
imagesc(Udf(:,:,2000));
axis image off;
title('Spatial component 2000')

% The temporal components represent how strong each spatial component is at
% each timepoint. The first one is usually big (that component explains
% lots of variance and the last one is usually very small (that component
% explains very little variance) - try zooming in to see the variance in
% the last component:
figure; hold on
plot(frame_t,fVdf([1,end],:)');
xlabel('Time (s)');
ylabel('Component weight');
legend({'Temporal component 1','Temporal component 2'});

% Pixel values for each frame ('reconstructing') can then be gotten by 
% matrix multiplcation of the spatial and temporal components 
% (this has some reshapes: the spatial components need to be changed from 
% 3D to 2D for matrix multiplcation then the resulting frame has to be 
% changed from one long vector into a Y pixels x X pixels size)
example_frame = 500;
example_fluorescence_long = ...
    reshape(Udf,[],size(Udf,3))*fVdf(:,example_frame);
example_fluorescence = ...
    reshape(example_fluorescence_long,size(Udf,1),size(Udf,2));
figure;
imagesc(example_fluorescence);
axis image off;
title('Example frame fluorescence');

% We've got this function to do the above quickly, so for example we can
% reconstruct a few frames together, which makes a 3D matrix of Y pixels x
% X pixels x N frames
example_frames = 500:600;
example_fluorescence = AP_svdFrameReconstruct(Udf,fVdf(:,example_frames));
% and then we can view that with a function I have for scrolling through 3D
% matricies: 
AP_imscroll(example_fluorescence);
axis image;

% The other nice thing about working with SVD data is that any linear
% operation can be done on the temporal component ('V space') before
% reconstruction ('pixel space'). For example, we can get the average
% activity across the whole day simply by taking the average V, then
% reconstructing: 
avg_V = nanmean(fVdf,2);
avg_fluorescence = AP_svdFrameReconstruct(Udf,avg_V);
figure;imagesc(avg_fluorescence);
axis image off
title('Average fluoresence');

% Most of the things we'd want to do are linear
% (adding/subtracting/multiplying/dividing: things that can be done in any
% order). Some things we might be interested in aren't linear so we'd have
% to do in pixel space (e.g. maximum fluorescence in each pixel, standard
% devation - because that involves a square root).

% EXERCISE: make a reward-triggered average movie, showing the average
% fluorescence -1:1 second around rewards. You found the reward times
% above, and since this is averaging you can do it in V-space before
% reconstructing into pixel space.












