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

read_frames = 1:min(4000,n_frames);
strobe_roi = zeros(length(read_frames),1);
for curr_frame_idx = 1:length(read_frames);
    curr_frame = read_frames(curr_frame_idx);
    curr_im = read(vr,curr_frame);
    strobe_roi(curr_frame_idx) = nanmean(curr_im(roiMask));
end
strobe_roi_thresh = prctile(strobe_roi,90)/2;
strobe_start_frame = read_frames(find(strobe_roi < strobe_roi_thresh,1));

% Find the end strobe (longer, within n frames?)
read_frames = max(1,(n_frames-4000)):n_frames;
strobe_roi = zeros(length(read_frames),1);
for curr_frame_idx = 1:length(read_frames);
    curr_frame = read_frames(curr_frame_idx);
    curr_im = read(vr,curr_frame);
    strobe_roi(curr_frame_idx) = nanmean(curr_im(roiMask));
end
strobe_roi_thresh = prctile(strobe_roi,90)/2;
strobe_stop_frame = read_frames(find(strobe_roi < strobe_roi_thresh,1));

cam_sync_frames = [strobe_start_frame,strobe_stop_frame];

if length(cam_sync_frames) ~= 2
    display('Cam sync not detected')
    cam_sync_frames = [];
    return
end

disp('Done')
