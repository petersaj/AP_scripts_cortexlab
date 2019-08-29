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
        im_resized{curr_im} = imresize(imread(im_fn{curr_im},curr_channel),im_rescale_factor);
    end
    waitbar(curr_im/n_im,h,'Loading and resizing images...');
end



figure;
im_montage = montage(im_resized);
caxis([0,max(im_montage.CData(:))])

% ROI method
% bg_roi = roipoly;
% signal_roi = roipoly;
% 
% figure; hold on
% histogram(im_montage.CData(bg_roi))
% histogram(im_montage.CData(signal_roi))
% 
% cmin = max(im_montage.CData(bg_roi));
% cmax = max(im_montage.CData(signal_roi));
% caxis([cmin,cmax]);


% Fit 2 mixed gaussians
gauss_fit = fitgmdist(single(im_montage.CData(im_montage.CData > 0)),2,'RegularizationValue',0.01);
cmin = norminv(0.99,gauss_fit.mu(1),sqrt(gauss_fit.Sigma(1)));
cmax = 2*norminv(0.99,gauss_fit.mu(2),sqrt(gauss_fit.Sigma(2)));
caxis([cmin,cmax]);


white_confirm = questdlg('White balance ok?');






[~,threshold] = edge(im_montage.CData(im_montage.CData > 0),'sobel');
fudgeFactor = 0.5;
BWs = edge(im_montage.CData,'sobel',threshold * fudgeFactor);
figure;imshow(BWs);




















