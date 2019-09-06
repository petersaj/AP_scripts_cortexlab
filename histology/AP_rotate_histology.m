function AP_rotate_histology(im_path)
% AP_rotate_histology(im_path)
%
% Pad, center, and rotate images of histological slices
% Andy Peters (peters.andrew.j@gmail.com)

slice_dir = dir([im_path filesep '*.tif']);
slice_fn = natsortfiles(cellfun(@(path,fn) [path filesep fn], ...
    {slice_dir.folder},{slice_dir.name},'uni',false));

slice_im = cell(length(slice_fn),1);
for curr_slice = 1:length(slice_fn)
   slice_im{curr_slice} = imread(slice_fn{curr_slice});  
end

% Pad all slices centrally to the largest slice and make matrix
slice_size_max = max(cell2mat(cellfun(@size,slice_im,'uni',false)),[],1);
slice_im_pad = ...
    cell2mat(cellfun(@(x) x(1:slice_size_max(1),1:slice_size_max(2),:), ...
    reshape(cellfun(@(im) padarray(im, ...
    [ceil((slice_size_max(1) - size(im,1))./2), ...
    ceil((slice_size_max(2) - size(im,2))./2)],0,'both'), ...
    slice_im,'uni',false),1,1,1,[]),'uni',false));

% Draw line to indicate midline for rotation
rotation_fig = figure;

align_axis = nan(2,2,length(slice_im));
for curr_im = 1:length(slice_im)
    imshow(slice_im_pad(:,:,:,curr_im));
    title('Click and drag reference line (e.g. midline)')
    curr_line = imline;
    align_axis(:,:,curr_im) = curr_line.getPosition;  
end
close(rotation_fig);

% NOTE FOR FUTURE: make target angle either vert or horiz for whatever's
% closest to the average reference? then switch dx vs dy below
target_angle = 90;
target_position = round(slice_size_max(2)./2);

% Get angle for all axes
align_angle = squeeze(acosd(diff(align_axis(:,1,:),[],1)./diff(align_axis(:,2,:),[],1)));
align_center = squeeze(nanmean(align_axis,1));

im_aligned = zeros(size(slice_im_pad),class(slice_im_pad));

for curr_im = 1:length(slice_im)
    
    angle_diff = align_angle(curr_im) - target_angle;
    x_diff = target_position - align_center(1,curr_im);
    y_diff = 0;
    
    im_aligned(:,:,:,curr_im) = ...
        imrotate(imtranslate(slice_im_pad(:,:,:,curr_im), ...
        [x_diff,y_diff]),angle_diff,'bilinear','crop');
    
end

% Overwrite old images with new ones
for curr_im = 1:length(slice_im)
    imwrite(im_aligned(:,:,:,curr_im),slice_fn{curr_im},'tif');
end

disp(['Saved padded/centered/rotated slices']);










