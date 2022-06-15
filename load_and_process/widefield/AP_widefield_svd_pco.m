function [U,Vrec,im_color_avg,frame_info] = AP_widefield_svd_pco(im_path)
% [U,Vrec,im_color_avg,frame_info] = AP_widefield_svd_pco(im_path)
%
% SVD-compress widefield imaging from PCO Edge 5.5 camera
% Assumes: 
% - alternating 2-color imaging that resets order on recording start
% - recordings are defined by timestamp gaps of >2s
% - binary and ASCII timestamps are turned on


%% Get image filenames

im_files = dir(fullfile(im_path,'*.tif'));


%% Get header information from all frames
% PCO timestamps in binary coded decimal (BCD) as:
% Pixels code: 
% 1:4          5 6 7 8 9 10 11 12    13    14
% framenum(x4) Y Y M D H M  S  10kus 100us us
% BCD converts to decimal as: a decimal number that codes for one byte, 
% with two nibbles that are one decimal each
% (e.g. 10dec = 16bcd: 16 -> 0001 0000 = "1,0" = 10)

tic;
disp('Getting image headers...');

% Set header position (first row, first 14 pixels)
header_px_loc = {[1,1],[1,14],[1,Inf]};

% Loop through all files, get frame number and timestamp
frame_info = struct('frame_num',cell(size(im_files)),'timestamp',cell(size(im_files)));
for curr_im = 1:length(im_files)

    % Grab pixels with header information
    curr_im_fn = fullfile(im_path,im_files(curr_im).name);
    pco_timestamps_bcd = permute(tiffreadVolume(curr_im_fn, ...
        'PixelRegion',header_px_loc),[2,3,1]);

    % Convert from BCD to decimal
    % (decimal to binary, split into 2 nibbles, binary to decimal)
    pco_timestamps_decimal = cellfun(@(x) num2str(bin2dec(reshape(dec2bin(x,8),4,[])')'), ...
        num2cell(pco_timestamps_bcd),'uni',false);

    % (frame number: join digits without spaces, convert to numbers)
    pco_frame_num = cellfun(@(x) str2num(strrep(x,' ','')), ...
        join(pco_timestamps_decimal(1:4,:)'));
    % (time stamp: get joined digits without spaces, then convert to datenum)
    pco_timestamp_str = cellfun(@(x) strrep(x,' ',''), ...
        join(pco_timestamps_decimal(5:13,:)'),'uni',false);
    pco_timestamp = cellfun(@(x) ...
        datenum(x(1:end-1),'yyyymmddHHMMSSFFF'), ...
        pco_timestamp_str);

    % Store header info for file
    frame_info(curr_im).frame_num = pco_frame_num;
    frame_info(curr_im).timestamp = pco_timestamp;

    AP_print_progress_fraction(curr_im,length(im_files));

end

toc;


%% Get illumination color for each frame
% Assume: blue/violet alternating, with blue starting whenever there is a >
% 2s gap (= new recording, which starts on blue)

% Sanity check: make sure frame numbers are consecutive
if any(diff(vertcat(frame_info.frame_num)) ~= 1)
    error('Frame numbers not consecutive');
end

% Get time difference between frames (convert datenum days to seconds)
timestamp_diff = diff(vertcat(frame_info.timestamp)*24*60*60);

% Find recording boundaries (time between frames > threshold)
recording_boundary_thresh = 2; % seconds between frames to define recording
recording_frame_boundary = [0,find(timestamp_diff > recording_boundary_thresh),frame_info(end).frame_num(end)]+1;

% Get recording for each frame
im_rec_idx = mat2cell( ....
    interp1(recording_frame_boundary,1:length(recording_frame_boundary), ...
    [1:frame_info(end).frame_num(end)],'previous'), ...
    1,cellfun(@length,{frame_info.frame_num}));
[frame_info.rec_idx] = im_rec_idx{:};

% Get illumination color of frame (alternate starting at each recording)
n_colors = 2;
im_frame_color = mat2cell( ...
    1+mod([1:frame_info(end).frame_num(end)] - ...
    interp1(recording_frame_boundary,recording_frame_boundary,...
    1:frame_info(end).frame_num(end),'previous'),n_colors), ...
    1,cellfun(@length,{frame_info.frame_num}));
[frame_info.color] = im_frame_color{:};


%% Set pixels to keep 
% (remove timestamp - first 8 pixel rows)

im_info = imfinfo(fullfile(im_path,im_files(1).name));
im_size = [im_info(1).Height,im_info(1).Width];

exclude_stamp_height = 8;

im_px_loc = {[exclude_stamp_height+1,im_size(1)],[1,im_size(2)],[1,Inf]};
im_grab_size = cellfun(@(x) diff(x)+1,im_px_loc(1:2));


%% Create moving- and total-average images

% Set moving average number of frames
n_frame_avg = 15;

tic;
disp('Building image moving average by color...');

% Loop through images, cumulatively build averages by illumination color
im_color_avg = zeros(im_grab_size(1),im_grab_size(2),n_colors);
im_color_mov_avg = cell(length(im_files),n_colors);
for curr_im = 1:length(im_files)
   
    curr_im_fn = fullfile(im_path,im_files(curr_im).name);
    im = single(tiffreadVolume(curr_im_fn,'PixelRegion',im_px_loc));

    % Loop through illumination colors
    for curr_color = 1:n_colors
        % Get color index for frames in current image
        curr_frame_color_idx = frame_info(curr_im).color == curr_color;

        % Cumulatively add average image
        curr_color_partial = ...
            sum(im(:,:,curr_frame_color_idx)./sum([frame_info.color] == curr_color),3);
        im_color_avg(:,:,curr_color) = im_color_avg(:,:,curr_color) + ...
            curr_color_partial;

        % Get moving average (truncate based on moving avg modulus)
        curr_n_frames = sum(curr_frame_color_idx);
        curr_frame_avg_idx = find(curr_frame_color_idx, ...
            curr_n_frames - mod(curr_n_frames,n_frame_avg));
        im_color_mov_avg{curr_im,curr_color} = ...
            permute(mean(reshape(im(:,:,curr_frame_avg_idx), ...
            size(im,1),size(im,2),n_frame_avg,[]),3),[1,2,4,3]);
    end

    AP_print_progress_fraction(curr_im,length(im_files));

end

toc;


%% Do SVD on moving-average images
% (keep U and S, don't keep V since U will be applied to full dataset next)

tic;
disp('Running SVD on moving average images by color...');

[U,~,~] = arrayfun(@(color) ...
    svd(reshape(cat(3,im_color_mov_avg{:,color}),prod(im_grab_size),[]),'econ'), ...
    1:n_colors,'uni',false);

% Reshape U into pixels row x pixels column x components
U = cellfun(@(x) reshape(x,im_grab_size(1),im_grab_size(2),[]),U,'uni',false);

toc;


%% Apply spatial components (U's) from SVD to full data
% (note: spatial components are applied as U' * mean-subtracted data, so
% the resulting matrix is equivalent to S*V but just called 'V')

tic;
disp('Applying SVD spatial components to full data...');

V = cell(length(im_files),n_colors);
for curr_im = 1:length(im_files)
   
    curr_im_fn = fullfile(im_path,im_files(curr_im).name);
    im = single(tiffreadVolume(curr_im_fn,'PixelRegion',im_px_loc));

    % Loop through illumination colors
    for curr_color = 1:n_colors
        % Get color index for frames in current image
        curr_frame_color_idx = frame_info(curr_im).color == curr_color;

        % Apply spatial components to mean-subtracted data
        V{curr_im,curr_color} = ...
            reshape(U{curr_color},[],size(U{curr_color},3))' * ...
            (reshape(im(:,:,curr_frame_color_idx),[],sum(curr_frame_color_idx)) - ...
            reshape(im_color_avg(:,:,curr_color),[],1));
    end

    AP_print_progress_fraction(curr_im,length(im_files));

end

toc;


%% Split V's by recording (instead of by file)

disp('Applying SVD spatial components to full data...');

% Store V's as recordings x color
Vrec = cell(length(recording_frame_boundary)-1,n_colors);

frame_num_cat = vertcat(frame_info.frame_num);

frame_color_cat = horzcat(frame_info.color);
frame_rec_idx_cat = horzcat(frame_info.rec_idx);
for curr_color = 1:n_colors

    % Split concatenated V by recording index
    curr_V_cat = horzcat(V{:,curr_color});
    Vrec(:,curr_color) = ...
        mat2cell(curr_V_cat,size(curr_V_cat,1), ...
        accumarray(frame_rec_idx_cat(frame_color_cat == curr_color),1));
    
end

disp('Finished SVD.')













































