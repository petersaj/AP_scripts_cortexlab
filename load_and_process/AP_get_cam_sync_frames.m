function [cam_sync_frames,n_frames] = AP_get_cam_sync_frames(fn,draw_roi)
% cam_strobe_frames = AP_get_cam_strobe_frames(fn)
%
% Get the frames for the onset of the sync strobes

if ~exist('draw_roi','var') || isempty(draw_roi)
    draw_roi = false;
end

% Find video, set up video reader
vr = VideoReader(fn);

% Get first frame, draw ROI for strobe detection
if draw_roi
    frame1 = readFrame(vr);
    f = figure;
    imagesc(frame1);colormap(gray);
    axis off
    title('Draw ROI to find cam sync');
    roiMask = roipoly;
    close(f);
    drawnow;
else
    % If no ROI is specified, use the lower half of the image (less likely
    % to have the scope illumination included)
    frame1 = readFrame(vr);
    roiMask = false(size(frame1));
    roiMask(round(size(roiMask,1)/2):end,:) = true;
end

% Assume that the strobe happens within first and last 1000 frames
clear vr
vr = VideoReader(fn);

disp('Getting sync frames...')

% Find the start strobe (within first n frames);
n_frames = vr.NumberOfFrames;

read_frames_start = 1:min(6000,n_frames);
strobe_roi_start = zeros(length(read_frames_start),1);
for curr_frame_idx = 1:length(read_frames_start)
    curr_frame = read_frames_start(curr_frame_idx);
    curr_im = read(vr,curr_frame);
    strobe_roi_start(curr_frame_idx) = nanmean(curr_im(roiMask));
end
strobe_roi_thresh = prctile(strobe_roi_start,50)/1.5;
strobe_start_frame_idx = find(strobe_roi_start < strobe_roi_thresh,1);
strobe_start_frame = read_frames_start(strobe_start_frame_idx);

% Find the end strobe (longer, within n frames?)
read_frames_end = max(1,(n_frames-10000)):n_frames;
strobe_roi_end = zeros(length(read_frames_end),1);
for curr_frame_idx = 1:length(read_frames_end)
    curr_frame = read_frames_end(curr_frame_idx);
    curr_im = read(vr,curr_frame);
    strobe_roi_end(curr_frame_idx) = nanmean(curr_im(roiMask));
end
strobe_roi_thresh = prctile(strobe_roi_end,50)/1.5;
strobe_end_frame_idx = find(strobe_roi_end < strobe_roi_thresh,1);
strobe_end_frame = read_frames_end(strobe_end_frame_idx);

cam_sync_frames = [strobe_start_frame,strobe_end_frame];

if length(cam_sync_frames) ~= 2
    display('Cam sync not detected')
    cam_sync_frames = [];
    return
end

% Plot the thresholds and crossings
figure;
subplot(2,2,1);
plot(strobe_roi_start,'k');
line(repmat(strobe_start_frame,1,2),ylim,'color','r','linewidth',2);
line(xlim,repmat(strobe_roi_start(strobe_start_frame_idx),2,1),'color','r','linewidth',2);
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
line(repmat(find(read_frames_end == strobe_end_frame),1,2),ylim,'color','r','linewidth',2);
line(xlim,repmat(strobe_roi_end(strobe_end_frame_idx),2,1),'color','r','linewidth',2);
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



