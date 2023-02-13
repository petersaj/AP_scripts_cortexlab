%% Tests for processing the mouse cameras (face/eye)


%% Get and save movement energy in movie

animal = 'AP100';
day = '2021-05-10';
experiments = 2;

AP_mouse_movie_movement(animal,day,experiments);


%% Get average movie around event

% animal = 'AP100';
% day = '2021-05-14';
% experiment = 2;
% load_parts.cam = true;
% verbose = true;
% AP_load_experiment

use_cam = facecam_fn;
use_t = facecam_t;

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
AP_imscroll(cam_align_avg,surround_t)
axis image;

surround_t = [-surround_frames:surround_frames]./vr.FrameRate;
AP_imscroll(cam_align_diff_avg,surround_t(2:end))
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
stim_col = [0,0,1;0,0,0;1,0,0];
AP_errorfill(t,stim_grp_mean',stim_group_sem',stim_col)
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

facecam_align = struct;
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
% save_fn = ['facecam_align'];
% save([save_path filesep save_fn],'facecam_align','-v7.3');
% disp(['Saved: ' save_path filesep save_fn])


%% Align facecam frames (nose/eye control point)

% Load facecam sample frames
facecam_processing_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\facecam_processing';
facecam_align_fn = fullfile(facecam_processing_path,'facecam_align.mat');
load(facecam_align_fn);

% Plot images and select control points
im_unaligned = cellfun(@double,[facecam_align.im],'uni',false);

target_im = im_unaligned{1};
target_size = size(target_im);

figure;
target_ax = subplot(1,2,1);
imagesc(target_im);
axis image off; hold on;
source_ax = subplot(1,2,2);
source_h = imagesc([]);
axis image off; hold on;

title(target_ax,'Click: nose eye');
target_ctrl_points = ginput(2);
plot(target_ax,target_ctrl_points(:,1),target_ctrl_points(:,2),'.r','MarkerSize',20);
title(target_ax,'');

im_aligned = nan(target_size(1),target_size(2),length(im_unaligned));
source_ctrl_points = cell(length(im_unaligned),1);
cam_tform = cell(length(im_unaligned),1);
for curr_im = 1:length(im_unaligned)
    source_im = im_unaligned{curr_im};
    if isempty(source_im)
        continue
    end

    % Click control points
    title(source_ax,{['Click: nose eye'], ...
        [sprintf('%d/%d',curr_im,length(im_unaligned))]});
    set(source_h,'CData',source_im);
    source_ctrl_points{curr_im} = ginput(2);

    % Store tform
    cam_tform{curr_im} = fitgeotrans(source_ctrl_points{curr_im}, ...
        target_ctrl_points,'nonreflectivesimilarity');
    tform_size = imref2d(target_size);
    im_aligned(:,:,curr_im) = ...
        imwarp(source_im,cam_tform{curr_im},'OutputView',tform_size);

end

% Plot aligned
AP_imscroll(im_aligned); axis image

% % Save transform (into original struct)
% % (package back into animals)
% n_days_animal = cellfun(@length,{facecam_align.day});
% cam_tform_animal = mat2cell(cam_tform,n_days_animal);
% [facecam_align.tform] = cam_tform_animal{:};
% 
% save(facecam_align_fn,'facecam_align');
% disp(['Saved: ' facecam_align_fn]);


%% Load and align (sanity check)

% Load facecam align
facecam_processing_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\facecam_processing';
facecam_align_fn = fullfile(facecam_processing_path,'facecam_align.mat');
load(facecam_align_fn);

% Align facecams
im_unaligned_cat = cellfun(@double,[facecam_align.im],'uni',false);
im_tform_cat = cat(1,facecam_align.tform);
im_ref = im_unaligned_cat{1};

im_aligned = nan(size(im_ref,1),size(im_ref,2),length(im_unaligned_cat));
for curr_im =1:length(im_unaligned_cat)
    if isempty(im_unaligned_cat{curr_im})
        continue
    end
    tform_size = imref2d(size(im_ref));
    im_aligned(:,:,curr_im) = ...
        imwarp(im_unaligned_cat{curr_im},im_tform_cat{curr_im}, ...
        'OutputView',tform_size);
end

AP_imscroll(im_aligned); axis image;

%% Draw facecam whisker ROI and align to each recording

% Load facecam align
facecam_processing_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\facecam_processing';
facecam_align_fn = fullfile(facecam_processing_path,'facecam_align.mat');
load(facecam_align_fn);

% Align facecams
im_unaligned_cat = cellfun(@double,[facecam_align.im],'uni',false);
im_tform_cat = cat(1,facecam_align.tform);
im_ref = im_unaligned_cat{1};

im_aligned = nan(size(im_ref,1),size(im_ref,2),length(im_unaligned_cat));
for curr_im =1:length(im_unaligned_cat)
    if isempty(im_unaligned_cat{curr_im})
        continue
    end
    tform_size = imref2d(size(im_ref));
    im_aligned(:,:,curr_im) = ...
        imwarp(im_unaligned_cat{curr_im},im_tform_cat{curr_im}, ...
        'OutputView',tform_size);
end

% Plot average image, draw whisker ROI
figure;imagesc(nanmean(im_aligned,3));axis image off;
title('Draw whisker line ROI');
whisker_line = drawline;
master_whisker_mask = createMask(whisker_line);
close(gcf);

figure;
image(imoverlay(mat2gray(nanmean(im_aligned,3)),master_whisker_mask,'r'));
axis image off
title('Master whisker mask');

% Align whisker mask to individual day
whisker_mask = cell(size(im_unaligned_cat));
for curr_im =1:length(im_unaligned_cat)
    if isempty(im_unaligned_cat{curr_im})
        continue
    end

    tform_size = imref2d(size(im_unaligned_cat{curr_im}));
    whisker_mask{curr_im} = ...
        imwarp(master_whisker_mask,invert(im_tform_cat{curr_im}), ...
        'OutputView',tform_size); 
end

% Package into original structure
n_days_animal = cellfun(@length,{facecam_align.day});
whisker_mask_animal = mat2cell(whisker_mask,1,n_days_animal);
[facecam_align.whisker_mask] = whisker_mask_animal{:};

% Plot through all days to check the aligned whisker ROIs
whisker_mask_cat = [facecam_align.whisker_mask];
figure; h = image(im_unaligned_cat{1}); axis image off;
for curr_im =1:length(im_unaligned_cat)
    if isempty(im_unaligned_cat{curr_im})
        continue
    end
    set(h,'CData', ...
        imoverlay(mat2gray(im_unaligned_cat{curr_im}), ...
        whisker_mask_cat{curr_im},'r'));
    pause(0.1);
end

% % Save whisker mask (into original struct)
% save(facecam_align_fn,'facecam_align');
% disp(['Saved: ' facecam_align_fn]);











