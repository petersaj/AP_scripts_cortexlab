%% Files to download

% I've copied some more files to
% \\zserver.cortexlab.net\Lab\Share\ajpeters\for_Andrada\widefield_alignment,
% copy that folder to your computer


%% Standardizing data for combining

% Each widefield recording has a unique set of spatial components (U),
% so to combine data across recordings we need to put our data into a set
% of common spatial components. My common spatial components, which I call
% my "master U", were made from SVDing U's across a bunch of recordings to
% capture consistent spatial patterns across recordings. I think that math
% is dubious, but it doesn't really matter: your U just needs to capture
% the variance in the data, so it can be a little arbitrary. 

% Let's load in an example data set
animal = 'AP025';
day = '2017-09-28';
experiment = 1;
verbose = true;
AP_load_experiment;

% I added the master U into your share folder, so just point this variable
% to wherever you've stored it on your local computer: 
master_u_fn = 'put your master U file path here';
load(master_u_fn);

% After loading above, you'll have a variable in your workspace U_master

% Converting widefield data from an experiment into this master U space
% needs 2 steps:

% 1) align the data (the master U is aligned to the master retinotopy, so
% the current day's data needs to be aligned to the master retinotopy)
Udf_aligned = AP_align_widefield(Udf,animal,day);

% 2) we use the ChangeU function to give us a new set of temporal
% components (V) corresponding to our master U
fVdf_Umaster = ChangeU(Udf_aligned,fVdf,U_master);

% EXERCISE: convince yourself (hopefully) that reconstructed data with the
% experiment U's and the master U's are similar. Try reconstructing a bunch
% of frames using Udf/fVdf and U_master/fVdf_Umaster and scroll through
% them side-by-side, do they look similar?

% EXERCISE: try comparing the ROIs from the experiment/master U's to check
% how similar they are. In order to do this, we'll need to define an ROI,
% and use the same ROI for both reconstructions.
%
% First draw an ROI on top of this figure using the 'roipoly' function
avg_im_aligned = AP_align_widefield(avg_im,animal,day);
figure;
imagesc(avg_im_aligned);
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
axis image off
% roi = (use roipoly here to get an ROI)
%
% You should have an ROI mask now (mask is true inside ROI, false
% outside),you can now put this into AP_svd_roi - so far you've only used
% this function to draw ROIs, but it can also take pre-drawn ROIs
experiment_U_trace = AP_svd_roi(Udf_aligned,fVdf,[],[],roi);
master_U_trace = AP_svd_roi(U_master,fVdf_Umaster,[],[],roi);

% Plot those traces on top of each other - are they similar? If so, it
% means we can swap out the experiment U for the master U and it gives us
% functionally the same data.

% Let's quantify how similar they are: calculate the fraction of explained
% variance when moving from the master U to the experiment U.

% When we're combining data, we don't necessarily want to keep all possible
% components since it takes up more room and is slower to work with. How
% many components should we keep? Get ROI traces using different number of
% components (3rd dimension in the U's, 1st dimension in the V's) and
% calculating the explained variance. Make a plot with explained variance
% vs. number of components.

% The raw data is dominated by large slow events, but we care about the
% faster deconvolved data. Make the explained variance vs. number of
% components plot above but this time using deconvolved data - is there a
% difference?

%% Combining data

% When we do our full analyses, we'll want to work on all of our data in a
% format that makes it easy to load and combine across recordings.

% EXERCISE: for a few of your mice (AP107/8/9 are probably a good bet),
% grab and save their passive data so you can analyze it all together. It
% will need these steps:
% 1) loop through each animal/recording
% 2) load that day's data
% 3) convert that day's data into the master U format
% 4) grab activity aligned to each stimulus
% 5) store the the activity and trial information, move to next loop
% 6) save
%
% Then write some code to load in the data and plot the average stimulus
% response across mice for each day (e.g. day 1 stimulus response averaged
% across all mice). What changes in the activity?














