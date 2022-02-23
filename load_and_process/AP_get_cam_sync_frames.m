function [cam_sync_frames,n_frames] = AP_get_cam_sync_frames(fn,draw_roi)
% [cam_strobe_frames,n_frames) = AP_get_cam_strobe_frames(fn,draw_roi)
%
% Get the frames for the onset of the sync strobes
% Input:
% fn - video filename
% draw_roi - true/false(default), manually draw roi or use lower half
% Output: 
% cam_sync_frames - video frames at strobe points
% n_frames - number of total frames in the video
%
% This assumes the minimum point = the strobe signal
% The sync is set as the first/last time the signal dips past thresh


% Known issue: disable erroneous matlab warning
warning off MATLAB:subscripting:noSubscriptsSpecified


if ~exist('draw_roi','var') || isempty(draw_roi)
    draw_roi = false;
end

% Find video, set up video reader
vr = VideoReader(fn);

% Get first frame, draw ROI for strobe detection
 frame1 = readFrame(vr);
if draw_roi
    f = figure;
    imagesc(frame1);colormap(gray);
    axis off
    title('Draw ROI to find cam sync');
    roi_mask = roipoly;
    close(f);
    drawnow;
else
    % If no ROI is specified, use a mask of the brightest pixels
    roi_mask_thresh = prctile(frame1(:),80);
    roi_mask = frame1 >= roi_mask_thresh;
end

% Assume that the strobe happens within first and last n frames
n_frames_check = 1000;

clear vr
vr = VideoReader(fn);

disp('Getting sync frames...')

% Find the start strobe
n_frames = vr.NumberOfFrames;

read_frames_start = 1:min(n_frames_check,n_frames);
strobe_roi_start = zeros(length(read_frames_start),1);
for curr_frame_idx = 1:length(read_frames_start)
    curr_frame = read_frames_start(curr_frame_idx);
    curr_im = read(vr,curr_frame);
    strobe_roi_start(curr_frame_idx) = nanmean(curr_im(roi_mask));
end
strobe_roi_thresh = min(strobe_roi_start) + ...
    (median(strobe_roi_start(n_frames_check/2:end)) - min(strobe_roi_start))*0.5;
strobe_start_frame_idx = find(diff(strobe_roi_start < strobe_roi_thresh) == 1,1,'first')+1;
strobe_start_frame = read_frames_start(strobe_start_frame_idx);

% Find the end strobe
% (refine ROI based on first strobe to ignore imaging lights)
strobe_start_framediff = read(vr,strobe_start_frame-1) - ...
    read(vr,strobe_start_frame);
roi_mask_thresh = prctile(strobe_start_framediff(:),80);
roi_mask = strobe_start_framediff >= roi_mask_thresh;

read_frames_end = max(1,(n_frames-n_frames_check)):n_frames;
strobe_roi_end = zeros(length(read_frames_end),1);
for curr_frame_idx = 1:length(read_frames_end)
    curr_frame = read_frames_end(curr_frame_idx);
    curr_im = read(vr,curr_frame);
    strobe_roi_end(curr_frame_idx) = nanmean(curr_im(roi_mask));
end
strobe_roi_thresh = min(strobe_roi_end) + ...
    (median(strobe_roi_end(1:n_frames_check/2)) - min(strobe_roi_end))*0.5;
strobe_end_frame_idx = find(diff(strobe_roi_end < strobe_roi_thresh) == 1,1,'last')+1;
strobe_end_frame = read_frames_end(strobe_end_frame_idx);

cam_sync_frames = [strobe_start_frame,strobe_end_frame];

if length(cam_sync_frames) ~= 2
    disp('Cam sync not detected')
    cam_sync_frames = [];
    return
end

% Plot the thresholds and crossings
figure;
subplot(2,2,1);
plot(strobe_roi_start,'k');
line(repmat(strobe_start_frame,1,2),ylim,'color','r','linestyle','--');
line(xlim,repmat(strobe_roi_start(strobe_start_frame_idx),2,1),'color','r','linestyle','--');
xlabel('Frame');
ylabel('Intensity');
title({'Recording start',fn},'interpreter','none')

subplot(2,4,5);
imagesc(read(vr,strobe_start_frame-1));
colormap(gray); axis image off; caxis([0,255]);
title('Pre-sync frame (start)');
subplot(2,4,6);
imagesc(read(vr,strobe_start_frame));
colormap(gray); axis image off; caxis([0,255]);
title('Sync frame (start)');

subplot(2,2,2);
plot(strobe_roi_end,'k');
line(repmat(find(read_frames_end == strobe_end_frame),1,2),ylim,'color','r','linestyle','--');
line(xlim,repmat(strobe_roi_end(strobe_end_frame_idx),2,1),'color','r','linestyle','--');
xlabel('Frame');
ylabel('Intensity');
title({'Recording stop',fn},'interpreter','none')
disp('Done')

subplot(2,4,7);
imagesc(read(vr,strobe_end_frame-1));
colormap(gray); axis image off; caxis([0,255]);
title('Pre-sync frame (end)');
subplot(2,4,8);
imagesc(read(vr,strobe_end_frame));
colormap(gray); axis image off; caxis([0,255]);
title('Sync frame (end)');

drawnow;

% Re-enable matlab warning
warning on MATLAB:subscripting:noSubscriptsSpecified

