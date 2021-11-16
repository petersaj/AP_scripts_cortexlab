%% Histology processing helpers

%% Concatenate split channel images into single files
% (images were saved in folders as separate files for each channel, combine
% into one file per image)

% Directory with images
im_split_path = '\\znas.cortexlab.net\Subjects\AP087\histology';
im_subfolder = 'Default'; % (Olympus subfolder default name)

% Directory to put channel-concatenated images
im_cat_path = [im_split_path filesep 'channelconcat_images'];

% Loop through subfolders, load and concatenate channels, save
if ~exist(im_cat_path,'dir')
    mkdir(im_cat_path)
end

im_dir = dir(im_split_path);
im_folders = natsortfiles(cellfun(@(x) ...
    [im_split_path filesep x], ...
    {im_dir(setdiff(find([im_dir.isdir]),[1,2])).name},'uni',false));

disp('Concatenating channels and saving...');
for curr_im_idx = 1:length(im_folders)
    % Get all tiff images in folder
    curr_im_dir = dir([im_folders{curr_im_idx} filesep im_subfolder filesep '*.tif*']);
    if isempty(curr_im_dir)
        continue
    end
    
    % First image: copy, subsequent images: append
    curr_im_channelcat_fn = [im_cat_path filesep 'image_' num2str(curr_im_idx) '.tif'];
    for curr_channel_idx = 1:length(curr_im_dir)
        curr_channel_fn = [curr_im_dir(curr_channel_idx).folder filesep ...
                curr_im_dir(curr_channel_idx).name];
        if curr_channel_idx == 1
            copyfile(curr_channel_fn,curr_im_channelcat_fn)
        else
            curr_channel_im = imread(curr_channel_fn);
            imwrite(curr_channel_im,curr_im_channelcat_fn,'tif', ...
                'Compression','none','WriteMode','append');
        end
    end
    
    AP_print_progress_fraction(curr_im_idx,length(im_dir));
    
end
disp('Done.');

%% Flip hemispheres in aligned histology
% for when I did it the wrong way

% Set location
im_path = '\\znas.cortexlab.net\Subjects\AP101\histology\all_images';
slice_path = [im_path filesep 'slices'];


mirror_matrix = eye(3).*[-1;1;1];


% Load in slice images
gui_data.slice_im_path = slice_path;
slice_im_dir = dir([slice_path filesep '*.tif']);
slice_im_fn = natsortfiles(cellfun(@(path,fn) [path filesep fn], ...
    {slice_im_dir.folder},{slice_im_dir.name},'uni',false));
gui_data.slice_im = cell(length(slice_im_fn),1);
for curr_slice = 1:length(slice_im_fn)
    gui_data.slice_im{curr_slice} = imread(slice_im_fn{curr_slice});
end

curr_slice = 10;
curr_slice_im = gui_data.slice_im{curr_slice};

curr_av_slice = histology_ccf(curr_slice).av_slices;
curr_av_slice(isnan(curr_av_slice)) = 1;

tform = affine2d;
tform.T = atlas2histology_tform{curr_slice}*mirror_matrix;
tform_size = imref2d([size(curr_slice_im,1),size(curr_slice_im,2)]);
curr_av_slice_warp = imwarp(curr_av_slice,tform,'nearest','OutputView',tform_size);

figure;imagesc(curr_av_slice_warp)
figure;imagesc(curr_slice_im)
















