function AP_process_histology
% AP_process_histology
%
% This is meant to be a non-shitty version of Process_Histology


%% Testing histology image preprocessing

% Get and sort image files
im_path = 'C:\Users\Andrew\Desktop\test_histology';
im_path_dir = dir([im_path filesep '*.ome.tiff']);
im_fn = natsortfiles(cellfun(@(path,fn) [path filesep fn], ...
    {im_path_dir.folder},{im_path_dir.name},'uni',false));


% Get microns/pixel from metadata (ome.tiff) 
% NOTE: apparently PhysicalSizeX isn't necessarily real??

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


n_im = length(im_fn);

im_resized = cell(n_im,3);

h = waitbar(0,'Loading and resizing images...');
for curr_im = 1:n_im
    for curr_channel = 1:3
        im_resized{curr_im,curr_channel} = imresize(imread(im_fn{curr_im},curr_channel),im_rescale_factor);
    end
    waitbar(curr_im/n_im,h,'Loading and resizing images...');
end
close(h);

h = figure;
im_montage = cell(3,1);
channel_caxis = nan(3,2);
for curr_channel = 1:3
    
   curr_montage = montage(im_resized(:,curr_channel));
   
   im_montage{curr_channel} = curr_montage.CData;
   
   % Estimate white balance
   % (very dirty: assume one peak for background, one for signal)
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

color_order = [2,3,1]; % (which channel is R,G,B)
im_montage_rgb = cell2mat(arrayfun(@(ch) rescale(im_montage{ch}, ...
    'InputMin',channel_caxis(ch,1),'InputMax',channel_caxis(ch,2)), ...
    permute(color_order,[1,3,2]),'uni',false));

figure;imshow(im_montage_rgb);


im_rgb = cellfun(@(x) zeros(size(x,1),size(x,2),size(im_resized,3),class(x)),im_resized(:,1),'uni',false);
for curr_im = 1:n_im
    im_rgb{curr_im} = cell2mat(arrayfun(@(ch) rescale(im_resized{curr_im,ch}, ...
        'InputMin',channel_caxis(ch,1),'InputMax',channel_caxis(ch,2)), ...
        permute(color_order,[1,3,2]),'uni',false));  
end



curr_im = 4;

min_blob = 1000;
curr_mask = bwareaopen(mean(im_rgb{curr_im},3) > 0.01,min_blob);
cc = bwconncomp(curr_mask);

blob_boundaries = bwboundaries(curr_mask);

figure;
imshow(im_rgb{curr_im}); hold on;
blob_line = gobjects(length(blob_boundaries),1);
for curr_blob = 1:length(blob_boundaries)
   blob_line(curr_blob) = line(blob_boundaries{curr_blob}(:,2), ...
       blob_boundaries{curr_blob}(:,1),'color',[0.5,0.5,0.5],'linewidth',3,'LineSmoothing','on');
end

title('Select slices');
[blob_x,blob_y] = ginput;
blob_ind = sub2ind(size(curr_mask),round(blob_y),round(blob_x));
selected_blob = cellfun(@(x) any(ismember(x,blob_ind)),cc.PixelIdxList);
set(blob_line(selected_blob),'color','w');


% dilate the mask, save that as the ROI


















