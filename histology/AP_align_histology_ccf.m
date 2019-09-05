function AP_align_histology_ccf(tv,av,st,slice_im_path)
% Grab CCF slices corresponding to histology slices

% Initialize guidata
gui_data = struct;
gui_data.tv = tv;
gui_data.av = av;
gui_data.st = st;

% Load in slice images
gui_data.slice_im_path = slice_im_path;
slice_im_dir = dir([slice_im_path filesep '*.tif']);
slice_im_fn = natsortfiles(cellfun(@(path,fn) [path filesep fn], ...
    {slice_im_dir.folder},{slice_im_dir.name},'uni',false));
gui_data.slice_im = cell(length(slice_im_fn),1);
for curr_slice = 1:length(slice_im_fn)
    gui_data.slice_im{curr_slice} = imread(slice_im_fn{curr_slice});
end



end


%% Storing here for now - transform
function asdf

example_slice_fn = 'C:\Users\Andrew\Desktop\test_histology\processed\ps_slices\12.tif';
example_slice = imread(example_slice_fn);
figure;imshow(example_slice);


allen_atlas_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']); % grey-scale "background signal intensity"
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
ccf_fig = allenAtlasBrowser(tv,av,st);


% tv_slice/av_slice = get ccf fig
[moving_points,fixed_points] = cpselect(example_slice,rescale(tv_slice),'Wait',true);



% this gui isn't good and need to have something that updates on every
% click, because otherwise no idea how it'll turn out

tform = fitgeotrans(moving_points,fixed_points,'affine');
R = imref2d(size(tv_slice));
example_slice_warp = imwarp(example_slice, tform, 'OutputView',R);


tform = fitgeotrans(fixed_points,moving_points,'affine');
R = imref2d(size(example_slice));
av_slice_warp = imwarp(av_slice, tform,'nearest','OutputView',R);


figure;

% slice -> ccf
slice_boundaries = bwboundaries(round(conv2(av_slice,ones(3)./9,'same')) ~= av_slice,4);

subplot(1,2,1);
imshow(example_slice_warp);

for i = 1:length(slice_boundaries)
   line(slice_boundaries{i}(:,2),slice_boundaries{i}(:,1), ...
       'color',[0.5,0.5,0.5],'linewidth',1); 
end
title('Slice -> CCF');

% ccf -> slice
slice_boundaries = bwboundaries(round(conv2(av_slice_warp,ones(3)./9,'same')) ~= av_slice_warp,4);

subplot(1,2,2);
imshow(example_slice);

for i = 1:length(slice_boundaries)
   line(slice_boundaries{i}(:,2),slice_boundaries{i}(:,1), ...
       'color',[0.5,0.5,0.5],'linewidth',1); 
end
title('CCF -> Slice');

% caxis([1,1305]);
% colormap(allen_ccf_colormap)
% axis image off





end




















