%% Get fluorescence in ROI

% Set path with mouse
im_path = '\\znas.cortexlab.net\Subjects\AP077\histology';
slice_path = [im_path filesep 'slices'];

% Load CCF atlas structure tree
allen_atlas_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

% Load CCF slices
ccf_slice_fn = [slice_path filesep 'histology_ccf.mat'];
load(ccf_slice_fn);

% Load histology/CCF alignment
ccf_alignment_fn = [slice_path filesep 'atlas2histology_tform.mat'];
load(ccf_alignment_fn);

% Set ROI, get areas within that area hierarchy
roi_area = 'Anterior cingulate area';

roi_id_path = st.structure_id_path(find(strcmp(st.safe_name,roi_area)));
roi_idx = find(contains(st.structure_id_path,roi_id_path));

% Loop through images, get fluorescence in ROI
slice_im_dir = dir([slice_path filesep '*.tif']);
slice_im_fn = natsortfiles(cellfun(@(path,fn) [path filesep fn], ...
    {slice_im_dir.folder},{slice_im_dir.name},'uni',false));

slice_fluor = nan(length(slice_im_fn),1);
for curr_slice = 1:length(slice_im_fn)
    
    % Load slice image
    curr_slice_im = imread(slice_im_fn{curr_slice});
    
    % Warp AV slice to align
    curr_slice_av_unaligned = histology_ccf(curr_slice).av_slices;
    curr_slice_av_unaligned(isnan(curr_slice_av_unaligned)) = 1;
    
    tform = affine2d;
    tform.T = atlas2histology_tform{curr_slice};
    tform_size = imref2d([size(curr_slice_im,1),size(curr_slice_im,2)]);
    curr_slice_av = imwarp(curr_slice_av_unaligned, ...
        tform,'nearest','OutputView',tform_size);
    
    % Get summed fluorescence of channel within ROI
    curr_slice_av_roimask = ismember(curr_slice_av,roi_idx);
    
    use_channel = 2;
    slice_fluor(curr_slice) = ...
        reshape(double(curr_slice_im(:,:,use_channel)),1,[])* ...
        reshape(curr_slice_av_roimask,[],1);
    
end

figure;
plot(slice_fluor)
ylabel(['Total fluorescence in ' roi_area]);
xlabel('Slice');









