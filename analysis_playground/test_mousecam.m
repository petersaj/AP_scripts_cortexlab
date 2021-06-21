%% Tests for processing the mouse cameras (face/eye)


%% Get and save movement energy in movie

animal = 'AP100';
day = '2021-05-10';
experiments = 2;

AP_mouse_movie_movement(animal,day,experiments);


%% Get average movie around event

% animal = 'AP100';
% day = '2021-05-10';
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
wheel_thresh = 0.025;
quiescent_trials = ~any(abs(event_aligned_wheel) > 0,2);

use_align = stimOn_times(stimIDs == 3 & quiescent_trials);
surround_frames = 30;

% Initialize video reader and average
vr = VideoReader(use_cam);
cam_im1 = read(vr,1);
cam_align_avg = zeros(size(cam_im1,1),size(cam_im1,2), ...
    surround_frames*2+1);

frame_t_offset = nan(size(use_align));
for curr_align = 1:length(use_align)
    
    % Find closest camera frame to timepoint
    [frame_t_offset(curr_align),curr_frame] = ...
        min(abs(use_align(curr_align) - use_t));
    
    % Pull surrounding frames
    curr_surround_frames = curr_frame + [-surround_frames,surround_frames];
    curr_clip = double(squeeze(read(vr,curr_surround_frames)));
    cam_align_avg = cam_align_avg + curr_clip./length(use_align);
    
    AP_print_progress_fraction(curr_align,length(use_align));
end

surround_t = [-surround_frames:surround_frames]./vr.FrameRate;
AP_image_scroll(cam_align_avg,surround_t)
axis image;














