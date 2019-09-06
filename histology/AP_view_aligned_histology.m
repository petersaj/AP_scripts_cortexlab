function AP_view_aligned_histology(tv,av,st,slice_im_path)
% AP_view_aligned_histology(tv,av,st,slice_im_path)
%
% View histology slices with overlaid aligned CCF areas
% Andy Peters (peters.andrew.j@gmail.com)

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

% Load corresponding CCF slices
ccf_slice_fn = [slice_im_path filesep 'histology_ccf.mat'];
load(ccf_slice_fn);
gui_data.histology_ccf = histology_ccf;

% Load histology/CCF alignment
ccf_alignment_fn = [slice_im_path filesep 'histology_ccf_alignment.mat'];
load(ccf_alignment_fn);
gui_data.histology_ccf_alignment = histology_ccf_alignment;

% Create figure, set button functions
gui_fig = figure('KeyPressFcn',@keypress,'WindowButtonMotionFcn',@mousehover);
gui_data.curr_slice = 1;

% Set up axis for histology image
gui_data.histology_ax = axes('YDir','reverse'); 
hold on; colormap(gray); axis image off;
gui_data.histology_im_h = image(gui_data.slice_im{1}, ...
    'Parent',gui_data.histology_ax);

% Set up histology-aligned atlas overlay
% (and make it invisible to mouse clicks)
histology_aligned_atlas_boundaries_init = ...
    zeros(size(gui_data.slice_im{1},1),size(gui_data.slice_im{1},2));
gui_data.histology_aligned_atlas_boundaries = ...
    imagesc(histology_aligned_atlas_boundaries_init,'Parent',gui_data.histology_ax, ...
    'AlphaData',histology_aligned_atlas_boundaries_init,'PickableParts','none');

% Upload gui data
guidata(gui_fig,gui_data);

% Update the slice
update_slice(gui_fig);

end


function keypress(gui_fig,eventdata)

% Get guidata
gui_data = guidata(gui_fig);

switch eventdata.Key
    
    % Left/right: move slice
    case 'leftarrow'
        gui_data.curr_slice = max(gui_data.curr_slice - 1,1);
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);
        
    case 'rightarrow'
        gui_data.curr_slice = ...
            min(gui_data.curr_slice + 1,length(gui_data.slice_im));
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);
        
end

end

function mousehover(gui_fig,eventdata)
% Display area of atlas on mouse hover
keyboard
end

function align_ccf_to_histology(gui_fig)

% Get guidata
gui_data = guidata(gui_fig);

curr_av_slice = gui_data.histology_ccf(gui_data.curr_slice).av_slices;
curr_av_slice(isnan(curr_av_slice)) = 1;
curr_slice_im = gui_data.slice_im{gui_data.curr_slice};

tform = fitgeotrans(gui_data.histology_ccf_alignment.atlas_control_points{gui_data.curr_slice}, ...
    gui_data.histology_ccf_alignment.histology_control_points{gui_data.curr_slice},'affine');
tform_size = imref2d([size(curr_slice_im,1),size(curr_slice_im,2)]);
curr_av_slice_warp = imwarp(curr_av_slice, tform, 'OutputView',tform_size);

av_warp_boundaries = round(conv2(curr_av_slice_warp,ones(3)./9,'same')) ~= curr_av_slice_warp;

set(gui_data.histology_aligned_atlas_boundaries, ...
    'CData',av_warp_boundaries, ...
    'AlphaData',av_warp_boundaries*0.3);

% Upload gui data
guidata(gui_fig, gui_data);

end

function update_slice(gui_fig)
% Draw histology and CCF slice

% Get guidata
gui_data = guidata(gui_fig);

% Set next histology slice
set(gui_data.histology_im_h,'CData',gui_data.slice_im{gui_data.curr_slice})

% Upload gui data
guidata(gui_fig, gui_data);

% Update atlas boundaries
align_ccf_to_histology(gui_fig)

end




















