%% TESTING: replacing pipelineHere

%% Get image filenames
im_path = 'G:\test_widefield';
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
im_header = struct('frame_num',cell(size(im_files)),'timestamp',cell(size(im_files)));
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
    im_header(curr_im).frame_num = pco_frame_num;
    im_header(curr_im).timestamp = pco_timestamp;

    AP_print_progress_fraction(curr_im,length(im_files));

end
toc;

%% Get illumination color for each frame
% Assume: blue/violet alternating, with blue starting whenever there is a >
% 2s gap (= new experiment, which starts on blue)

% Sanity check: make sure frame numbers are consecutive
if any(diff(vertcat(im_header.frame_num)) ~= 1)
    error('Frame numbers not consecutive');
end

% Get time difference between frames (convert datenum days to seconds)
timestamp_diff = diff(vertcat(im_header.timestamp)*24*60*60);

% Find experiment boundaries (time between frames > threshold)
exp_boundary_thresh = 2; % seconds between frames to define experiment
exp_frame_boundary = [0,find(timestamp_diff > exp_boundary_thresh),im_header(end).frame_num(end)]+1;

% Get illumination color of frame (alternate from experiment starts)
frame_color = 1+mod([1:im_header(end).frame_num(end)] - ...
    interp1(exp_frame_boundary,exp_frame_boundary,...
    1:im_header(end).frame_num(end),'previous'),2);

%% For image processing: ignore pixels with timestamps
% (first 8 pixel rows)

im_info = imfinfo(fullfile(im_path,im_files(1).name));
im_size = [im_info(1).Height,im_info(1).Width];

exclude_stamp_height = 8;

im_px_loc = {[exclude_stamp_height+1,im_size(1)],[1,im_size(2)],[1,Inf]};
im_grab_size = cellfun(@(x) diff(x)+1,im_px_loc(1:2));

%% Get average image for each color

tic;
disp('Building average images by color...');

% Loop through images, build average image for each color
im_color_avg = zeros(im_grab_size(1),im_grab_size(2),max(frame_color));
for curr_im = 1:length(im_files)
   
    curr_im_fn = fullfile(im_path,im_files(curr_im).name);
    im = tiffreadVolume(curr_im_fn,'PixelRegion',im_px_loc);

    % Cumulatively add partial color average
    for curr_color = 1:max(frame_color)
        curr_frame_color_idx = frame_color(im_header(curr_im).frame_num) == curr_color;
        curr_color_partial = ...
            sum(double(im(:,:,curr_frame_color_idx))./sum(frame_color == curr_color),3);
        im_color_avg(:,:,curr_color) = im_color_avg(:,:,curr_color) + ...
            curr_color_partial;
    end

    AP_print_progress_fraction(curr_im,length(im_files));

end

toc;

%% Build covariance matrix for each color 

tic;
disp('Building covariance matricies by color...');

im_color_cov = zeros(prod(im_grab_size),prod(im_grab_size),max(frame_color));

for curr_im = 1:length(im_files)
   
    curr_im_fn = fullfile(im_path,im_files(curr_im).name);
    im = tiffreadVolume(curr_im_fn,'PixelRegion',im_px_loc);

    % Cumulatively add partial color covariance
    for curr_color = 1:max(frame_color)
        curr_frame_color_idx = frame_color(im_header(curr_im).frame_num) == curr_color;
        curr_color_partial = ...
            sum(double(im(:,:,curr_frame_color_idx))./sum(frame_color == curr_color),3);
        im_color_avg(:,:,curr_color) = im_color_avg(:,:,curr_color) + ...
            curr_color_partial;
    end

    AP_print_progress_fraction(curr_im,length(im_files));

end

toc;





%% (faster than imread)


tic; V = tiffreadVolume(curr_im_fn); toc;

%% Read section

header_px_loc = {[1,1],[1,14],[1,Inf]};
tic; V = tiffreadVolume(curr_im_fn,'PixelRegion',header_px_loc); toc;

%% Read alternate frames

header_px_loc = {[1,Inf],[1,Inf],[1,2,Inf]};
tic; V = tiffreadVolume(curr_im_fn,'PixelRegion',header_px_loc); toc;



%% Do datastore to keep all files accessible without loading together?
% https://uk.mathworks.com/matlabcentral/answers/350406-how-do-i-read-large-set-of-multi-page-tiff-files-into-tall-array-using-a-datastore

ds = fileDatastore({im_path},'ReadFcn', @tiffreadVolume);


ds = imageDatastore(im_path,'FileExtensions','.tif');

%% 

tic;read(ds);toc;

%% Load current lilrig options 

load('bluePurpleOps.mat');


% averages 7500 frames together??
ops.Nframes = size(V,3)*17;


%% Build up covariance matrix?

% Nicks loaded in ~1000 frames and averaged groups of ~15 frames (~500ms)

% grab one channel
V2 = single(V(:,:,1:2:6540));

% option 1: average groups of N frames? 
V3 = squeeze(nanmean(reshape(V2,size(V2,1),size(V2,2),15,[]),3));

% option 2: just get covariance of these pixels?
V2s = single(V2);
x = reshape(V2s,[],size(V2,3))'*reshape(V2s,[],size(V2,3));




