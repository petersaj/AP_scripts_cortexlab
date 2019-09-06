function AP_get_probe_histology(slice_im_path)
% AP_get_probe_histology(tv,av,st,slice_im_path)
%
% Get probe trajectory in hostology and convert to ccf
% Andy Peters (peters.andrew.j@gmail.com)

% Load in slice images
slice_im_path = slice_im_path;
slice_im_dir = dir([slice_im_path filesep '*.tif']);
slice_im_fn = natsortfiles(cellfun(@(path,fn) [path filesep fn], ...
    {slice_im_dir.folder},{slice_im_dir.name},'uni',false));
slice_im = cell(length(slice_im_fn),1);
for curr_slice = 1:length(slice_im_fn)
    slice_im{curr_slice} = imread(slice_im_fn{curr_slice});
end

% Load corresponding CCF slices
ccf_slice_fn = [slice_im_path filesep 'histology_ccf.mat'];
load(ccf_slice_fn);

% Load histology/CCF alignment
ccf_alignment_fn = [slice_im_path filesep 'histology_ccf_alignment.mat'];
load(ccf_alignment_fn);

% Draw line along probe trajectory
probe_fig = figure;

probe_points_histology = nan(2,2,length(slice_im));
for curr_slice = 1:length(slice_im)
    imshow(slice_im{curr_slice});
    title('Click and drag probe trajectory (or just click if no probe)')
    curr_line = imline;
    % If the line is just a click, don't include
    curr_line_length = sqrt(sum(abs(diff(curr_line.getPosition,[],1)).^2));
    if curr_line_length == 0
        continue
    end   
    probe_points_histology(:,:,curr_slice) = curr_line.getPosition;  
end
close(probe_fig);

% Convert probe points to CCF points by alignment
probe_ccf = zeros(0,3);
for curr_slice = 1:length(slice_im)
    
    % Transform histology to atlas slice
    tform = fitgeotrans(histology_ccf_alignment.histology_control_points{curr_slice}, ...
        histology_ccf_alignment.atlas_control_points{curr_slice},'affine');    
    probe_points_atlas_slice = round([probe_points_histology(:,:,curr_slice),zeros(2,1)]*tform.T);
    
    % Get CCF coordinates corresponding to atlas slice points
    use_points = find(~all(isnan(probe_points_atlas_slice),2));
    for curr_point = 1:length(use_points)
        ccf_x = histology_ccf(curr_slice). ...
            plane_x(probe_points_atlas_slice(curr_point,1), ...
            probe_points_atlas_slice(curr_point,2));
        ccf_y = histology_ccf(curr_slice). ...
            plane_y(probe_points_atlas_slice(curr_point,1), ...
            probe_points_atlas_slice(curr_point,2));
        ccf_z = histology_ccf(curr_slice). ...
            plane_z(probe_points_atlas_slice(curr_point,1), ...
            probe_points_atlas_slice(curr_point,2));
        probe_ccf = vertcat(probe_ccf,[ccf_x,ccf_y,ccf_z]);
    end 
    
end

% Create 3D axis with CCF outlne
figure;
axes_atlas = axes;
[~, brain_outline] = plotBrainGrid([],axes_atlas);
hold(axes_atlas,'on');
axis vis3d equal off manual
view([-30,25]);
caxis([0 300]);
[ap_max,dv_max,ml_max] = deal(1320,800,1140);
xlim([-10,ap_max+10])
ylim([-10,ml_max+10])
zlim([-10,dv_max+10])

% Plot points and line of best fit
r0 = mean(probe_ccf,1);
xyz = bsxfun(@minus,probe_ccf,r0);
[~,~,V] = svd(xyz,0);
histology_probe_direction = V(:,1);

line_eval = [-1000,1000];
probe_fit_line = bsxfun(@plus,bsxfun(@times,line_eval',histology_probe_direction'),r0);
plot3(probe_ccf(:,1),probe_ccf(:,2),probe_ccf(:,3),'.r','MarkerSize',20);
line(probe_fit_line(:,1),probe_fit_line(:,2),probe_fit_line(:,3),'color','k','linewidth',2)

% Make 3D rotation the default state (toggle on/off with 'r')
h = rotate3d(axes_atlas);
h.Enable = 'on';

% Save probe CCF points
save_fn = [slice_im_path filesep 'probe_ccf.mat'];
save(save_fn,'probe_ccf');











