%% Repos to download 

% Clone these repos and add them to your path:
% https://github.com/cortex-lab/Lilrig

% I've copied some files to
% \\zserver.cortexlab.net\Lab\Share\ajpeters\for_Andrada\widefield_alignment,
% copy that folder to your computer

% You'll have to edit 2 functions to use the path where you copied the
% above files: 
% edit AP_align_widefield, enter the path in line 30
% edit AP_reference_outline, enter the path in line 30


%% Retinotopy, alignment, and atlas overlay
% Now that we can work with widefield data, we want to know where different
% brain areas are in our image. The atlas we use is from the Allen
% Institute, sometimes called the common coordinate framework (CCF) or the
% Allen reference atlas (ARA). Here's the paper describing it that has some
% nice figures listing all the top-down areas we see in widefield:
% https://www.sciencedirect.com/science/article/pii/S0092867420304025?via%3Dihub
%
% We need to align the CCF atlas to our image to figure out which areas are
% which, and we also need something to align from one animal to another so
% we can combine our data. We do this both of these with retinotopy.
%
% 'Retinotopy' means there's a 1:1 map of visual space onto the brain,
% which happens in all visual cortical areas - this gives us a few maps,
% one for primary visual cortex and one for all secondary visual areas. The
% most recent mapping we use for reference is in this paper: 
% https://elifesciences.org/articles/18372
%
% Each visual area in that paper has a color, red or blue, which is the
% direction of the map. In V1 the map is flipped from the screen: 
% leftward on the screen gives rightward on the brain, which is colored 
% as blue. This reverses at each visual step, so in visual area AM for 
% example leftward on the screen is leftward on the brain, which is 
% colored red. 
%
% We get our retinotopic maps from a "sparse noise" stimulus: show a bunch
% of dots on the screen at random, get a screen-to-brain mapping, then
% color each map by whether it is oriented the same way as the screen or
% flipped. Here's an example recording and retinotopic map:

% Load in a sparse noise experiment
animal = 'AP025';
day = '2017-09-29';
experiment = 2;
load_parts.imaging = true;
verbose = true;
AP_load_experiment;

% Run code to pull out the retinotopy (takes a few minutes):
lilrig_retinotopy

% That gives you the retinotopic map (blue/red) by itself, overlaid on the
% average image, and with the CCF atlas areas aligned and overlaid.

% I align mice by aligning the retinotopic map for each mouse to a "master"
% retinotopic map (which is the aligned average of a bunch of mice). I
% then store the alignment transformation for each mouse/day, which I load
% to do the alignment. In the first part of this script, you copied the
% alignment file and changed the paths in the relevant files. 

% Let's align the example day's retinotopic map, which is in the variable
% 'vfs_boot_mean': 
vfs_boot_mean_aligned = AP_align_widefield(vfs_boot_mean,animal,day);

% We can plot the aligned retinotopic map, and then use the function
% AP_reference_outline to draw the aligned CCF regions on top
figure; imagesc(vfs_boot_mean_aligned); axis image off;
AP_reference_outline('ccf_aligned','k');

% The big yellow blob is V1, which should fit into the top part of V1 in
% the CCF. The little yellow blob sticking out of the top of V1 is visual
% area AM - that should be in the corner of the triangle that's VisAM.

% For reference, let's do that plot for the unaligned data: you can see it
% doesn't match up as well (it's not super far off, because each animal is
% imaged in roughly the same orientation).
figure; imagesc(vfs_boot_mean); axis image off;
AP_reference_outline('ccf_aligned','k');

%% Load example dataset 
% (for the next things, let's use the same example dataset as before)

% Here's an example animal/day
animal = 'AP025';
day = '2017-09-28';
experiment = 1;
verbose = true;

% This is my catch-all code for loading data
% (finds the data on the server, loads, processes)
AP_load_experiment;

%% Align stimulus responses

% EXERCISE: you previously got the average widefield response to each
% stimulus - do that again, but with master-aligned data. Our widefield
% data comes in U's (spatial components) and V's (spatial components), so
% all you have to do is align the U's once (feed them into
% AP_align_widefield) to get aligned U's and then do your analysis on that
% data. View those movies with AP_image_scroll, and overlay the CCF areas
% with AP_reference_outline. 

%% Regions of interest (ROIs)

% When we want to plot a line trace for fluorescence in one area, we can
% draw a region of interest (ROI), average all the pixels in the ROI, and
% plot the ROI fluorescence across time. 



%% Deconvolution




%% Hemodynamics

% Our imaging is done with fluorescence, excitation light of one wavelength
% (blue or violet) triggers light emission from GCaMP in a different
% wavelength (green). We use two different excitation colors because they
% each excite GCaMP in a different way: 






















