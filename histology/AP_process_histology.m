function AP_process_histology(im_path)
% AP_process_histology(im_path)
%
% Resize and white balance histology images and extract images of each slice
% Andy Peters (peters.andrew.j@gmail.com)

% Get and sort image files
im_path_dir = dir([im_path filesep '*.tif*']);
im_fn = natsortfiles(cellfun(@(path,fn) [path filesep fn], ...
    {im_path_dir.folder},{im_path_dir.name},'uni',false));

% Get microns/pixel from metadata (if ome.tiff) 
im_info = imfinfo(im_fn{1});
im_description = im_info(1).ImageDescription;

im_um = regexp(im_description,'PhysicalSizeX="(\S*)".*PhysicalSizeY="(\S*)"','tokens');
im_um_x = str2num(im_um{1}{1});
im_um_y = str2num(im_um{1}{2});

if im_um_x ~= im_um_y
   error('Pixel X/Y values different') 
end

% Set resize factor to match to Allen CCF
allen_um2px = 10;
im_rescale_factor = im_um_x/allen_um2px;

% Load and resize images
n_im = length(im_fn);
im_resized = cell(n_im,3);

h = waitbar(0,'Loading and resizing images...');
for curr_im = 1:n_im
    for curr_channel = 1:3
        im_resized{curr_im,curr_channel} = imresize(imread(im_fn{curr_im},curr_channel),im_rescale_factor);
    end
    waitbar(curr_im/n_im,h,['Loading and resizing images (' num2str(curr_im) '/' num2str(n_im) ')...']);
end
close(h);

% Estimate white balance within each channel
% (dirty: assume one peak for background, one for signal)
h = figure;
im_montage = cell(3,1);
channel_caxis = nan(3,2);
for curr_channel = 1:3
    
   curr_montage = montage(im_resized(:,curr_channel));
   
   im_montage{curr_channel} = curr_montage.CData;

   im_hist = histcounts(im_montage{curr_channel}(im_montage{curr_channel} > 0),0:max(im_montage{curr_channel}(:)));
   im_hist_deriv = diff(smooth(im_hist,10));
   
   [~,bg_median] = max(im_hist);
   bg_signal_min = find(im_hist_deriv(bg_median:end) > 0,1) + bg_median;
   [~,bg_median_rel] = max(im_hist(bg_signal_min:end));
   signal_median = bg_median_rel + bg_signal_min;
   
   cmin = bg_signal_min;
   cmax = signal_median*3;
   caxis([cmin,cmax]);
   
   channel_caxis(curr_channel,:) = [cmin,cmax];
   
   white_confirm = questdlg('White balance ok?');
   
end
close(h)

% Display montage of final balanced image, sort color channels by RGB
color_order = [2,3,1]; % (which channel is R,G,B)
im_montage_rgb = cell2mat(arrayfun(@(ch) rescale(im_montage{ch}, ...
    'InputMin',channel_caxis(ch,1),'InputMax',channel_caxis(ch,2)), ...
    permute(color_order,[1,3,2]),'uni',false));
figure;imshow(im_montage_rgb);

% Store RGB for each slide
im_rgb = cellfun(@(x) zeros(size(x,1),size(x,2),size(im_resized,3),class(x)),im_resized(:,1),'uni',false);
for curr_im = 1:n_im
    im_rgb{curr_im} = cell2mat(arrayfun(@(ch) rescale(im_resized{curr_im,ch}, ...
        'InputMin',channel_caxis(ch,1),'InputMax',channel_caxis(ch,2)), ...
        permute(color_order,[1,3,2]),'uni',false));  
end

% Set up GUI to pick slices on slide to extract
slice_fig = figure('KeyPressFcn',@slice_keypress);

% Initialize data
slice_data = struct;
slice_data.im_path = im_path;
slice_data.im_fn = im_fn;
slice_data.im_rescale_factor = im_rescale_factor;
slice_data.im_rgb = im_rgb;
slice_data.curr_slide = 0;
slice_data.slice_mask = cell(0,0);
slice_data.slice_rgb = cell(0,0);

% Update gui data
guidata(slice_fig, slice_data);

% Update slide
update_slide(slice_fig);

end

function slice_click(slice_fig,eventdata)
% On slice click, mark to extract

slice_data = guidata(slice_fig);

selected_slice_bw = bwselect(slice_data.mask,eventdata.IntersectionPoint(1),eventdata.IntersectionPoint(2));

if eventdata.Button == 1
    % If left button pressed, create new slice ROI
    roi_num = size(slice_data.user_masks,3) + 1;
    
    % Make new mask with object
    slice_data.user_masks(:,:,roi_num) = selected_slice_bw;
    
    % Draw bounding box around object
    box_x = find(any(slice_data.user_masks(:,:,roi_num),1),1);
    box_y = find(any(slice_data.user_masks(:,:,roi_num),2),1);
    box_w = find(any(slice_data.user_masks(:,:,roi_num),1),1,'last') - box_x;
    box_h = find(any(slice_data.user_masks(:,:,roi_num),2),1,'last') - box_y;
    slice_data.user_rectangles(roi_num) = ...
        rectangle('Position',[box_x,box_y,box_w,box_h],'EdgeColor','w');
       
elseif eventdata.Button == 3
    % if right button pressed, join to last slice ROI
    roi_num = size(slice_data.user_masks,3);
    
    % Join old and new objects in mask
    slice_data.user_masks(:,:,roi_num) = ...
        slice_data.user_masks(:,:,roi_num) | selected_slice_bw;
    
    % Draw bounding box around object
    box_x = find(any(slice_data.user_masks(:,:,roi_num),1),1);
    box_y = find(any(slice_data.user_masks(:,:,roi_num),2),1);
    box_w = find(any(slice_data.user_masks(:,:,roi_num),1),1,'last') - box_x;
    box_h = find(any(slice_data.user_masks(:,:,roi_num),2),1,'last') - box_y;
    set(slice_data.user_rectangles(roi_num),'Position', ...
        [box_x,box_y,box_w,box_h]);
    
end

% Update gui data
guidata(slice_fig, slice_data);

end

function slice_keypress(slice_fig,eventdata)
% Move to next slide with spacebar

if strcmp(eventdata.Key,'space')
    update_slide(slice_fig) 
end

end

function update_slide(slice_fig)
% Find slices on slide by over-threshold objects of a large enough size

slice_data = guidata(slice_fig);

% Pull the images from selected slices (not during initialization)
if slice_data.curr_slide > 0
    extract_slice_rgb(slice_fig);
    slice_data = guidata(slice_fig);
end

% After the last slice, save the images and close out
if slice_data.curr_slide == length(slice_data.im_rgb)
    save_slice_rgb(slice_fig);
    close(slice_fig);
    return
end

slice_data.curr_slide = slice_data.curr_slide + 1;

min_slice = (1000/10)^2; % (um/10(CCF units))^2
slice_mask = bwareaopen(mean( ...
    slice_data.im_rgb{slice_data.curr_slide},3) > 0.01,min_slice);
slice_conncomp = bwconncomp(slice_mask);

im_handle = imshow(slice_data.im_rgb{slice_data.curr_slide});
set(im_handle,'ButtonDownFcn',@slice_click);
title('Click slices to extract (left = new, right = add to last), spacebar to finish slide');

slice_boundaries = bwboundaries(slice_mask);
slice_lines = gobjects(length(slice_boundaries),1);
for curr_slice = 1:length(slice_boundaries)
    slice_lines(curr_slice) = line(slice_boundaries{curr_slice}(:,2), ...
        slice_boundaries{curr_slice}(:,1),'color','w','linewidth',2,'LineSmoothing','on','linestyle','--');
end

slice_data.im_h = im_handle;
slice_data.mask = slice_mask;
slice_data.lines = slice_lines;
slice_data.user_masks = zeros(size(slice_mask,1),size(slice_mask,2),0,'logical');
slice_data.user_rectangles = gobjects(0);

% Update gui data
guidata(slice_fig, slice_data);

end


function extract_slice_rgb(slice_fig)
% When changing slide, extract the selected slice images

slice_data = guidata(slice_fig);

n_slices = size(slice_data.user_masks,3);
curr_slice_mask = cell(n_slices,1);
curr_slice_rgb = cell(n_slices,1);
for curr_slice = 1:n_slices
    % Pull a rectangular area, exclude spaces (e.g. between torn piece)
    dilate_size = 30;
    curr_mask = imdilate(logical(any(slice_data.user_masks(:,:,curr_slice),2).* ...
        any(slice_data.user_masks(:,:,curr_slice),1)),ones(dilate_size));
    
    curr_rgb = reshape(slice_data.im_rgb{slice_data.curr_slide}( ...
        repmat(curr_mask,1,3)),sum(any(curr_mask,2)),sum(any(curr_mask,1)),3);
    
    curr_slice_mask{curr_slice} = curr_mask;
    curr_slice_rgb{curr_slice} = curr_rgb;
    
end

% Store the image and mask for each slice
slice_data.slice_mask{slice_data.curr_slide} = curr_slice_mask;
slice_data.slice_rgb{slice_data.curr_slide} = curr_slice_rgb;

% Update gui data
guidata(slice_fig, slice_data);

end


function save_slice_rgb(slice_fig)
% After the last slide, save the slice images

slice_data = guidata(slice_fig);

% Set save directory as subdirectory within original
save_dir = [slice_data.im_path filesep 'slices'];
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end

% Concatenate all slice images
slice_rgb_cat = vertcat(slice_data.slice_rgb{:});

% Write all slice images to separate files
for curr_im = 1:length(slice_rgb_cat)
    curr_fn = [save_dir filesep num2str(curr_im) '.tif'];
    imwrite(slice_rgb_cat{curr_im},curr_fn,'tif');
end

% Get rows and columns for each slice corresponding to full size image
slice_slide_locations = cell(size(slice_data.slice_mask));
for curr_slide = 1:length(slice_data.slice_mask)
    for curr_slice = 1:length(slice_data.slice_mask{curr_slide})
        
        curr_mask = slice_data.slice_mask{curr_slide}{curr_slice};
        
        mask_x = find(interp1(1:size(curr_mask,2),+any(curr_mask,1), ...
            linspace(1,size(curr_mask,2), ...
            round(size(curr_mask,2)/slice_data.im_rescale_factor)),'nearest'));
        mask_y = find(interp1(1:size(curr_mask,1),+any(curr_mask,2), ...
            linspace(1,size(curr_mask,1), ...
            round(size(curr_mask,1)/slice_data.im_rescale_factor)),'nearest'));
        
        slice_slide_locations{curr_slide}{curr_slice} = ...
            {mask_y,mask_x};
        
    end
end

slice_slide_locations_fn = [save_dir filesep 'slice_slide_locations.mat'];
save(slice_slide_locations_fn,'slice_slide_locations');

disp(['Slices saved in ' save_dir]);

end



