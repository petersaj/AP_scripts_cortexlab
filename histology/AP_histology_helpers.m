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












