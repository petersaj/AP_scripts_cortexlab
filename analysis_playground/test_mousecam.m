%% Tests for processing the mouse cameras (face/eye)


%% Get and save movement energy in movie

animal = 'AP100';
day = '2021-05-10';
experiments = 2;

AP_mouse_movie_movement(animal,day,experiments);


%% Get average movie around event
% (ADD INTO THIS: average abs movie difference)

% animal = 'AP100';
% day = '2021-05-14';
% experiment = 2;
% load_parts.cam = true;
% verbose = true;
% AP_load_experiment

use_cam = facecam_fn; % facecam_fn, eyecam_fn
use_t = facecam_t; % facecam_t, eyecam_t

% Get wheel movements during stim, only use quiescent trials
framerate = 30;
wheel_window = [0,0.5];
wheel_window_t = wheel_window(1):1/framerate:wheel_window(2);
wheel_window_t_peri_event = bsxfun(@plus,stimOn_times,wheel_window_t);
event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
    wheel_velocity,wheel_window_t_peri_event);
wheel_thresh = 0;
quiescent_trials = ~any(abs(event_aligned_wheel) > 0,2);

use_stim = 3;
use_align = stimOn_times(stimIDs == use_stim & quiescent_trials);
surround_frames = 30;

% Initialize video reader, get average and average difference
vr = VideoReader(use_cam);
cam_im1 = read(vr,1);

cam_align_avg = zeros(size(cam_im1,1),size(cam_im1,2), ...
    surround_frames*2+1);
cam_align_diff_avg = zeros(size(cam_im1,1),size(cam_im1,2), ...
    surround_frames*2);

frame_t_offset = nan(size(use_align));
for curr_align = 1:length(use_align)
    
    % Find closest camera frame to timepoint
    [frame_t_offset(curr_align),curr_frame] = ...
        min(abs(use_align(curr_align) - use_t));
    
    % Pull surrounding frames
    curr_surround_frames = curr_frame + [-surround_frames,surround_frames];
    curr_clip = double(squeeze(read(vr,curr_surround_frames)));
    curr_clip_diff = abs(diff(curr_clip,[],3));

    cam_align_avg = cam_align_avg + curr_clip./length(use_align);
    cam_align_diff_avg = cam_align_diff_avg + curr_clip_diff./length(use_align);

    AP_print_progress_fraction(curr_align,length(use_align));
end

surround_t = [-surround_frames:surround_frames]./vr.FrameRate;
AP_image_scroll(cam_align_avg,surround_t)
axis image;

surround_t = [-surround_frames:surround_frames]./vr.FrameRate;
AP_image_scroll(cam_align_diff_avg,surround_t(2:end))
axis image;

% Plot difference within window
use_t = [0,0.2];
use_t_idx = surround_t >= use_t(1) & surround_t <= use_t(2);
figure;
imagesc(nanmean(cam_align_diff_avg(:,:,use_t_idx(2:end)),3));
axis image off;
title(use_stim);


%% Get line ROI for whiskers

AP_mousemovie(facecam_fn,facecam_t)
whisker_line = drawline;
whisker_mask = createMask(whisker_line);

vr = VideoReader(facecam_fn);
n_frames = vr.NumFrames;

whisker_px = nan(n_frames,sum(whisker_mask(:)));
for curr_frame = 1:n_frames
    curr_im = read(vr,curr_frame);
    whisker_px(curr_frame,:) = curr_im(whisker_mask);
    AP_print_progress_fraction(curr_frame,n_frames);
end

whisker_move = sum(abs(diff(whisker_px,[],1)),2);
whisker_move_t = facecam_t(2:end);

figure;
plot(whisker_move_t,whisker_move);
xlabel('Time (s)');
ylabel('Whisker movement');

% Get wheel movements during stim, only use quiescent trials
framerate = 30;
wheel_window = [0,0.5];
wheel_window_t = wheel_window(1):1/framerate:wheel_window(2);
wheel_window_t_peri_event = bsxfun(@plus,stimOn_times,wheel_window_t);
event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
    wheel_velocity,wheel_window_t_peri_event);
wheel_thresh = 0;
quiescent_trials = ~any(abs(event_aligned_wheel) > 0,2);

% Stim-align and plot whisker times
align_times = stimOn_times;
surround_window = [-0.2,1];

surround_samplerate = 50;
t = surround_window(1):1/surround_samplerate:surround_window(2);
peri_event_t = reshape(align_times,[],1) + reshape(t,1,[]);

stim_aligned_whisker = interp1(whisker_move_t(~isnan(whisker_move_t)), ...
    whisker_move(~isnan(whisker_move_t)),peri_event_t,'previous');

[stim_grp_mean,stim_group_sem] = grpstats( ...
    stim_aligned_whisker(quiescent_trials,:), ...
    stimIDs(quiescent_trials),{'nanmean','sem'});
figure; hold on
AP_errorfill(t,stim_grp_mean',stim_group_sem')
xline(0);
xlabel('Time from stim (s)');
ylabel('Whisker movement');



%% Get motion energy of whole frame - exclude wf-illumination pixels

vr = VideoReader(facecam_fn);

% (frame before sync = infrared + widefield)
ir_wf_faceframe = double(read(vr,find(facecam_t < ...
    Timeline.rawDAQTimestamps(camSync_flip(1)),1,'last')));
% (frame during sync = just widefield)
wf_faceframe = double(read(vr,find(facecam_t > ...
    Timeline.rawDAQTimestamps(camSync_flip(1)),1,'first')));

% (get IR mask: remove bright pixels that don't drop with IR strobe)
faceframe_wf_mask = (ir_wf_faceframe > prctile(ir_wf_faceframe,50)) & ...
    (wf_faceframe./ir_wf_faceframe > 0.5);
ir_mask_leeway = 10; % border around wf pixels to exclude
faceframe_ir_mask = ~(imdilate(faceframe_wf_mask,ones(ir_mask_leeway)));

sync_frame_diff = ir_wf_faceframe - wf_faceframe;

figure;
subplot(1,3,1);
imagesc(ir_wf_faceframe);
axis image off;
subplot(1,3,2);
imagesc(wf_faceframe);
axis image off;
subplot(1,3,3);
imagesc(imoverlay(sync_frame_diff,~faceframe_ir_mask,'r'));
axis image off;



% Open experiment movie
n_frames = vr.NumberOfFrames;
n_frames_surround = 1000;

% Set chunks to load in (for memory)
chunk_frames_n = 5000;
chunk_frames = round(linspace(1,n_frames,ceil(n_frames/chunk_frames_n)+1));

frame_movement = nan(n_frames,1);
for curr_chunk = 1:length(chunk_frames)-1
    % (back up one frame if not at the start to fill in all diffs)
    curr_chunk_frames = max(chunk_frames(curr_chunk)-1,1):chunk_frames(curr_chunk+1);

    % (load movie, get abs diff)
    curr_chunk_movie = permute((read(vr, ...
        [curr_chunk_frames(1),curr_chunk_frames(end)])),[1,2,4,3]);
    curr_chunk_movie_diff = abs(diff(double(curr_chunk_movie),[],3));

    % (get sum movie difference in IR pixels)
    curr_movement = ...
        reshape(curr_chunk_movie_diff,[],size(curr_chunk_movie_diff,3))'* ...
        faceframe_ir_mask(:); 
    frame_movement(curr_chunk_frames(2:end)) = curr_movement;

    AP_print_progress_fraction(curr_chunk,length(chunk_frames)-1);
end



% Get wheel movements during stim, only use quiescent trials
framerate = 30;
wheel_window = [0,0.5];
wheel_window_t = wheel_window(1):1/framerate:wheel_window(2);
wheel_window_t_peri_event = bsxfun(@plus,stimOn_times,wheel_window_t);
event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
    wheel_velocity,wheel_window_t_peri_event);
wheel_thresh = 0;
quiescent_trials = ~any(abs(event_aligned_wheel) > 0,2);

% Stim-align and plot whisker times
align_times = stimOn_times;
surround_window = [-0.2,1];

surround_samplerate = 100;
t = surround_window(1):1/surround_samplerate:surround_window(2);
peri_event_t = reshape(align_times,[],1) + reshape(t,1,[]);

stim_aligned_whisker = interp1(whisker_move_t(~isnan(whisker_move_t)), ...
    frame_movement(~isnan(whisker_move_t)),peri_event_t,'previous');

[stim_grp_mean,stim_group_sem] = grpstats( ...
    stim_aligned_whisker(quiescent_trials,:), ...
    stimIDs(quiescent_trials),{'nanmean','sem'});
figure; hold on
AP_errorfill(t,stim_grp_mean',stim_group_sem')
xline(0);


%% ~~~~~~~~~~~~~ TESTING BATCH
% (moved to AP_operant_learning_preprocessing and developing from there)

%% Collect sample facecam frame from all experiments

animals = {'AP100','AP101','AP103','AP104','AP105','AP106','AP107','AP108','AP109','AP111','AP112'};

% Initialize save variable
facecam_sampleframe_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);
    
    % Get days with muscimol
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    if any(muscimol_animal_idx)
        muscimol_start_day = muscimol(muscimol_animal_idx).day{1};
        muscimol_experiments = datenum({experiments.day})' >= datenum(muscimol_start_day);
    else
        muscimol_experiments = false(size({experiments.day}));
    end
    
    % Set experiments to use (imaging, not muscimol)
    experiments = experiments([experiments.imaging] & ~muscimol_experiments);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment (only cam)
        load_parts.cam = true;
        AP_load_experiment       

        if ~facecam_exists
            continue
        end

        % Grab sample frame (last frame)
        vr = VideoReader(facecam_fn);
        grab_frame = vr.NumFrames;
        facecam_sample_frame = read(vr,grab_frame);
        
        facecam_sampleframe_all(curr_animal).animal = animal;
        facecam_sampleframe_all(curr_animal).day{curr_day} = day;
        facecam_sampleframe_all(curr_animal).im{curr_day} = facecam_sample_frame;

        % Prep for next loop
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars('-except',preload_vars{:});

    end
end
disp('Done loading all');

% % Save
% save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\facecam_processing';
% save_fn = ['facecam_sampleframe_all'];
% save([save_path filesep save_fn],'facecam_sampleframe_all','-v7.3');
% disp(['Saved: ' save_path filesep save_fn])

%% Align facecam frames (auto)
% Align facecam sample frames across mice/days
% (copied from AP_align_widefield)

% Load facecam sample frames
facecam_processing_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\facecam_processing';
facecam_sampleframe_all_fn = fullfile(facecam_processing_path,'facecam_sampleframe_all.mat');
load(facecam_sampleframe_all_fn);

im_unaligned = cellfun(@double,[facecam_sampleframe_all.im],'uni',false);

% Pad the images for checking
% (set output size as the largest image)
[im_y,im_x] = cellfun(@size,im_unaligned);

im_y_max = max(im_y);
im_x_max = max(im_x);
ref_size = [im_y_max,im_x_max];

im_unaligned_pad = ...
    cell2mat(reshape(cellfun(@(im) padarray(im, ...
    [im_y_max - size(im,1),im_x_max - size(im,2)],0,'post'), ...
    im_unaligned,'uni',false),1,1,[]));

AP_image_scroll(im_unaligned_pad);
axis image;

% (for RegularStepGradientDescent)
[optimizer, metric] = imregconfig('monomodal');
optimizer.GradientMagnitudeTolerance = 1e-6;
optimizer.MaximumIterations = 100;
optimizer.MaximumStepLength = 1e-4;
optimizer.MinimumStepLength = 1e-6;
optimizer.RelaxationFactor = 0.6;

im_unaligned_edge = cellfun(@(x) x-imgaussfilt(x,10),im_unaligned,'uni',false);

% (rigid align to first image)
disp('Rigid aligning images...')
% (use only top half of image - excludes hands which are random)
im_ref = im_unaligned_edge{1}(1:round(size(im_unaligned{1},1)/2),:);

im_rigid_aligned = nan(ref_size(1),ref_size(2),length(im_unaligned));
rigid_tform = cell(size(im_unaligned));
for curr_im = 1:length(im_unaligned)
    if isempty(im_unaligned{curr_im})
        continue
    end
    
    im_source = im_unaligned_edge{curr_im}(1:round(size(im_unaligned{curr_im},1)/2),:);

    tformEstimate_affine = imregtform(im_source, ...
        im_ref,'rigid',optimizer,metric,'PyramidLevels',4);
    curr_im_reg = imwarp(im_unaligned{curr_im}, ...
        tformEstimate_affine,'Outputview',imref2d(ref_size));
    rigid_tform{curr_im} = tformEstimate_affine.T;
    im_rigid_aligned(:,:,curr_im) = curr_im_reg;
    
    if exist('h','var');delete(h);end
    h = figure;
    imshowpair(im_ref,im_rigid_aligned(:,:,curr_im));
    drawnow

    AP_print_progress_fraction(curr_im,length(im_unaligned));
end

AP_image_scroll(im_rigid_aligned)
axis image;

% Save
facecam_align_refsize = ref_size;
facecam_align_tform = rigid_tform;

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\facecam_processing';
save_fn = ['facecam_align_tform'];
save([save_path filesep save_fn],'facecam_align_refsize','facecam_align_tform','-v7.3');
disp(['Saved: ' save_path filesep save_fn])

%% Align facecam frames (control point)

% Load facecam sample frames
facecam_processing_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\facecam_processing';
facecam_sampleframe_all_fn = fullfile(facecam_processing_path,'facecam_sampleframe_all.mat');
load(facecam_sampleframe_all_fn);

% Plot images and select control points
im_unaligned = cellfun(@double,[facecam_sampleframe_all.im],'uni',false);

target_im = im_unaligned{1};
source_im = im_unaligned{20};

figure;
target_ax = subplot(1,2,1);
imagesc(target_im);
axis image off; hold on;
source_ax = subplot(1,2,2);
imagesc(source_im);
axis image off; hold on;

target_ctrl_points = ginput;
plot(target_ax,target_ctrl_points(:,1),target_ctrl_points(:,2),'.r','MarkerSize',20);

source_ctrl_points = ginput;
plot(source_ax,source_ctrl_points(:,1),source_ctrl_points(:,2),'.b','MarkerSize',20);

% Align images
tform = fitgeotrans(source_ctrl_points,target_ctrl_points,'nonreflectivesimilarity');
tform_size = imref2d(size(target_im));
source_im_aligned = imwarp(source_im,tform,'OutputView',tform_size);

target_im_edge = target_im - imgaussfilt(target_im,10);



%% Draw facecam ROI and align to each recording

% (TO DO)














