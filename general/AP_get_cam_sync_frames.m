function [cam_sync_frames,n_frames] = AP_get_cam_sync_frames(fn)
% cam_strobe_frames = AP_get_cam_strobe_frames(fn)
%
% Get the frames for the onset of the sync strobes

% Find video, set up video reader
vr = VideoReader(fn);

% Get first frame, draw ROI for strobe detection
frame1 = readFrame(vr);
f = figure;
imagesc(frame1);colormap(gray);
axis off
title('Draw ROI to find cam sync');
roiMask = roipoly;
close(f);
drawnow;

% Assume that the strobe happens within first and last 1000 frames
clear vr
vr = VideoReader(fn);

% Find the start strobe (within first 1000 frames);
n_frames = vr.NumberOfFrames;

read_frames = 1:min(6000,n_frames);
strobe_roi = zeros(length(read_frames),1);
for curr_frame_idx = 1:length(read_frames);
    curr_frame = read_frames(curr_frame_idx);
    curr_im = read(vr,curr_frame);
    strobe_roi(curr_frame_idx) = nanmean(curr_im(roiMask));
end
strobe_roi_thresh = prctile(strobe_roi,90)/2;
strobe_start_frame = read_frames(find(strobe_roi < strobe_roi_thresh,1));

% Find the end strobe (longer, within 4000 frames?)
read_frames = max(1,(n_frames-6000)):n_frames;
strobe_roi = zeros(length(read_frames),1);
for curr_frame_idx = 1:length(read_frames);
    curr_frame = read_frames(curr_frame_idx);
    curr_im = read(vr,curr_frame);
    strobe_roi(curr_frame_idx) = nanmean(curr_im(roiMask));
end
strobe_roi_thresh = prctile(strobe_roi,90)/2;
strobe_stop_frame = read_frames(find(strobe_roi < strobe_roi_thresh,1));

cam_sync_frames = [strobe_start_frame,strobe_stop_frame];
