% Generate revision figures for ctx-str paper
% (trials data is prepared in AP_ctx_str_trial_preprocessing)

% NOTE: these are just in order that I wrote them at the moment

%% [[LOAD DATASETS]]

% Load data

% (task)
% data_fn = 'trial_activity_choiceworld'; % Primary dataset
% data_fn = 'trial_activity_choiceworld_4strdepth'; % Depth-aligned striatum
% exclude_data = true;

% (passive)
data_fn = 'trial_activity_AP_choiceWorldStimPassive_trained';
% data_fn = 'trial_activity_AP_choiceWorldStimPassive_naive';
% data_fn = 'trial_activity_stimKalatsky_naive';
% data_fn = 'trial_activity_stimKalatsky_trained';
exclude_data = false;

% (unused at the moment)
% data_fn = 'trial_activity_choiceworld_wfonly'; % Widefield-only days (no craniotomy, so cleaner)
% exclude_data = true;

AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_recording;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));


%% CCF domain location / Allen corticostriatal projections

warning('THIS IS NOW SUBBED INTO MAIN FIGURES')

% Load the kernel template matches
n_aligned_depths = 3;
kernel_match_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing';
kernel_match_fn = ['ephys_kernel_align_' num2str(n_aligned_depths) '_depths.mat'];
load([kernel_match_path filesep kernel_match_fn]);

% Pad and concatenate kernel matches
max_n_kernel_match = max(cell2mat(cellfun(@(x) ...
    cellfun(@length,x),{ephys_kernel_align.kernel_match},'uni',false)));

kernel_match_all = cellfun(@(x) cell2mat(cellfun(@(x) ...
    padarray(x,max_n_kernel_match-length(x),NaN,'pre'),x,'uni',false)), ...
    {ephys_kernel_align.kernel_match},'uni',false);

kernel_match_cat = cell2mat(kernel_match_all);

% (use only depths in at least x% of recordings)
frac_depths = nanmean(~isnan(kernel_match_cat),2);
use_depths = frac_depths > 0.5;
kernel_match_cat_use = kernel_match_cat(use_depths,:);

% Get fraction of each domain for each depth
kernel_match_frac = nan(sum(use_depths),n_aligned_depths);
for curr_depth = 1:n_aligned_depths
    curr_kernel_match = +(kernel_match_cat_use == curr_depth);
    curr_kernel_match(isnan(kernel_match_cat_use)) = NaN;
    kernel_match_frac(:,curr_depth) = nanmean(curr_kernel_match,2);
end

% Get center-of-mass for each domain
kernel_match_com = nan(n_aligned_depths,1);
for curr_depth = 1:n_aligned_depths
    curr_kernel_match = +(kernel_match_cat_use == curr_depth);
    curr_kernel_match(isnan(kernel_match_cat_use)) = NaN;
        kernel_match_com(curr_depth) = ...
            nansum(nanmean(curr_kernel_match,2).*(1:sum(use_depths))') ...
            ./nansum(nanmean(curr_kernel_match,2));
end

% Plot domain locations
str_col = copper(n_aligned_depths);
figure; hold on;
colormap(str_col);
area(kernel_match_frac,'FaceColor','flat');
for curr_depth = 1:n_aligned_depths
    line(repmat(kernel_match_com(curr_depth),2,1),ylim, ...
        'color',min(str_col(curr_depth,:)+0.2,1),'linewidth',3);
end
axis tight;
xlabel('Estimated depth');
ylabel('Fraction domain match');

% Get relative domain COM along trajectory and plot in CCF
kernel_match_com_relative = kernel_match_com./sum(use_depths);

% Define the probe vector manually according to the targeted trajectory
probe_vector_ccf = [520,240,510;520,511,239];

% Load in the annotated Allen volume and names
allen_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
av = readNPY([allen_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean

% Get probe location per micron
probe_size = pdist2(probe_vector_ccf(1,:),probe_vector_ccf(2,:))*10;
probe_depths = ...
    round([linspace(probe_vector_ccf(1,1)',probe_vector_ccf(2,1)',probe_size); ...
    linspace(probe_vector_ccf(1,2)',probe_vector_ccf(2,2)',probe_size); ...
    linspace(probe_vector_ccf(1,3)',probe_vector_ccf(2,3)',probe_size)]');

% Eliminiate trajectory points that are off the atlas
eliminate_depths = ...
    probe_depths(:,1) < 1 | probe_depths(:,1) > size(av,1) | ...
    probe_depths(:,2) < 1 | probe_depths(:,2) > size(av,2) | ...
    probe_depths(:,3) < 1 | probe_depths(:,3) > size(av,3);
probe_depths(eliminate_depths,:) = [];

% Convert probe depths subscripts to indicies
probe_depths_ind = sub2ind(size(av),probe_depths(:,1),probe_depths(:,2),probe_depths(:,3));

% Get structures that the probe is in
probe_structures = av(probe_depths_ind);

% Get target relative depths through striatum
str_id = find(strcmp(st.safe_name,'Caudoputamen'));
probe_structures_str = probe_structures == str_id;
probe_str = [probe_depths(find(probe_structures_str,1,'first'),:); ...
    probe_depths(find(probe_structures_str,1,'last'),:)];
kernel_depth_ccf = interp1([0,1],probe_str,kernel_match_com_relative);

% Plot brain to overlay probes
% (note the CCF is rotated to allow for dim 1 = x)
h = figure; ccf_axes = axes; hold on
slice_spacing = 10;
target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) > 1,[2,1,3]);
structure_patch = isosurface(target_volume,0);
structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
target_structure_color = [0.7,0.7,0.7];
brain_outline = patch('Vertices',structure_wire.vertices*slice_spacing, ...
    'Faces',structure_wire.faces, ...
    'FaceColor','none','EdgeColor',target_structure_color);

str_id = find(strcmp(st.safe_name,'Caudoputamen'));
target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) == str_id,[2,1,3]);
structure_patch = isosurface(target_volume,0);
structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
target_structure_color = [0.7,0,0.7];
striatum_outline = patch('Vertices',structure_wire.vertices*slice_spacing, ...
    'Faces',structure_wire.faces, ...
    'FaceColor','none','EdgeColor',target_structure_color);

axis image vis3d off;
view([-30,25]);
cameratoolbar(h,'SetCoordSys','y');
cameratoolbar(h,'SetMode','orbit');

scatter3(kernel_depth_ccf(:,1),kernel_depth_ccf(:,2),kernel_depth_ccf(:,3), ...
    100,copper(n_aligned_depths),'filled');

% Plot brain area outlines at slice
callosum_id = find(strcmp(st.safe_name,'corpus callosum body'));
ventricle_id = find(strcmp(st.safe_name,'lateral ventricle'));

av_slice = permute(av(probe_vector_ccf(1),:,:),[2,3,1]);
slice_brain_outline = bwboundaries(av_slice > 1,'noholes');
slice_str_outline = bwboundaries(av_slice == str_id,'noholes');
slice_callosum_outline = bwboundaries(av_slice == callosum_id,'noholes');
slice_ventricle_outline = bwboundaries(av_slice == ventricle_id,'noholes');

figure; hold on; axis image on; box on; grid on; set(gca,'YDir','reverse');
plot(slice_brain_outline{1}(:,2),slice_brain_outline{1}(:,1),'k','linewidth',2);
cellfun(@(x) plot(x(:,2),x(:,1),'b','linewidth',2),slice_str_outline);
cellfun(@(x) fill(x(:,2),x(:,1),'k','linewidth',2),slice_callosum_outline);
cellfun(@(x) fill(x(:,2),x(:,1),'k','linewidth',2),slice_ventricle_outline);

line(probe_vector_ccf(:,3),probe_vector_ccf(:,2),'linewidth',2,'color','r');
scatter(kernel_depth_ccf(:,3),kernel_depth_ccf(:,2), ...
    100,copper(n_aligned_depths),'filled');

% Mirror all of the locations across the midline and get projections from
% both (because the cortical injections in the database aren't evenly
% distributed, so this is more comprehensive)
kernel_depth_um = round(kernel_depth_ccf*10);
bregma_um = allenCCFbregma*10;
ccf_midline = bregma_um(3);
hemisphere = sign(kernel_depth_um(1,3) - ccf_midline);
str_depths_mirror = [kernel_depth_um(:,1:2),ccf_midline - hemisphere*abs(kernel_depth_um(:,3)-ccf_midline)];

str_depths_query = [kernel_depth_um;str_depths_mirror];
injection_parameters = get_allen_projection(str_depths_query);
injection_coordinates = {injection_parameters.coordinates};

% Standardize injection coordinates by hemisphere (left = contra, right =
% ipsi)
injection_coordinates_standardized = injection_coordinates;
for curr_coord = 1:length(injection_coordinates)
    
    target_hemisphere = sign(ccf_midline - str_depths_query(curr_coord,3));
    injection_coords_ml_offset = abs(injection_coordinates{curr_coord}(:,3) - ccf_midline);
    injection_coordinates_hemisphere = sign(injection_coordinates{curr_coord}(:,3) - ccf_midline);
    
    injection_coords_ipsi = injection_coordinates_hemisphere == target_hemisphere;
    injection_coords_contra = injection_coordinates_hemisphere == -target_hemisphere;
    
    injection_coordinates_standardized{curr_coord}(injection_coords_ipsi,3) = ...
        ccf_midline + injection_coords_ml_offset(injection_coords_ipsi);
    injection_coordinates_standardized{curr_coord}(injection_coords_contra,3) = ...
        ccf_midline - injection_coords_ml_offset(injection_coords_contra);
    
end

% Combine side-standardized injection data across hemispheres
injection_coordinates_bilateral = arrayfun(@(x) ...
    vertcat(injection_coordinates_standardized{x}, ...
    injection_coordinates_standardized{n_aligned_depths+x}), ...
    1:n_aligned_depths,'uni',false);

% Convert points from CCF to widefield
alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment';
ccf_tform_fn = [alignment_path filesep 'ccf_tform.mat'];
load(ccf_tform_fn);

um2pixel = 20.6;
injection_coordinates_bilateral_wf = cellfun(@(x) ...
    [x(:,[3,1]).*(1/(um2pixel)),ones(size(x,1),1)]*ccf_tform.T, ...
    injection_coordinates_bilateral,'uni',false);

% Get projection density
projection_strength = {injection_parameters.density};

projection_strength_bilateral = arrayfun(@(x) ...
    vertcat(projection_strength{x}', ...
    projection_strength{n_aligned_depths+x}'), ...
    1:n_aligned_depths,'uni',false);

% Load kernel templates for overlay
n_aligned_depths = 3;
kernel_template_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\kernel_template_' num2str(n_aligned_depths) '_depths.mat'];
load(kernel_template_fn);

figure; 
colormap(brewermap([],'*RdBu'));
for curr_depth = 1:n_aligned_depths
    subplot(n_aligned_depths,1,curr_depth); hold on; axis image off;
    set(gca,'YDir','reverse');
    
    scatter(injection_coordinates_bilateral_wf{curr_depth}(:,1), ...
        injection_coordinates_bilateral_wf{curr_depth}(:,2), ...
        projection_strength_bilateral{curr_depth}*50 + 10, ...
        'k','filled');
    
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
end


    
    

%% Goodness-of-fit (cortex and striatum)


% Regress kernel ROI activity to striatum domain activity (per recording)
regression_params.kernel_t = [0,0];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;
lambda = 0;
kernel_frames = floor(regression_params.kernel_t(1)*sample_rate): ...
    ceil(regression_params.kernel_t(2)*sample_rate);

mua_ctxroi_k = arrayfun(@(x) nan(length(kernel_frames),n_depths),1:length(use_split),'uni',false);
mua_ctxroipred_allcat = nan(size(mua_allcat));
for curr_exp = 1:length(trials_recording)
    for curr_depth = 1:n_depths
        
        curr_mua = reshape(mua_allcat(split_idx == curr_exp,:,curr_depth)',[],1)';
        curr_fluor_kernelroi = reshape(fluor_kernelroi_deconv(split_idx == curr_exp,:,curr_depth)',[],1)';
        
        % Skip if no data
        if all(isnan(curr_mua))
            continue
        end
        
        % Set discontinuities in trial data
        trial_discontinuities = false(size(mua_allcat(split_idx == curr_exp,:,curr_depth)));
        trial_discontinuities(:,1) = true;
        trial_discontinuities = reshape(trial_discontinuities',[],1)';
        
        % Do regression
        [k,curr_mua_kernelroipred,explained_var] = ...
            AP_regresskernel(curr_fluor_kernelroi, ...
            curr_mua,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,regression_params.use_constant,trial_discontinuities);
              
        mua_ctxroi_k{curr_exp}(:,curr_depth) = k;
        mua_ctxroipred_allcat(split_idx == curr_exp,:,curr_depth) = ...
            reshape(curr_mua_kernelroipred,length(t),[])';

    end
end






% Get R^2 for task, cortex full, and cortex ROI predictions
taskpred_r2 = nan(max(split_idx),n_depths);
ctxpred_r2 = nan(max(split_idx),n_depths);
ctxroipred_r2 = nan(max(split_idx),n_depths);
for curr_exp = 1:max(split_idx)
       
    curr_data = reshape(permute(mua_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_taskpred_data = reshape(permute(mua_taskpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_ctxpred_data = reshape(permute(mua_ctxpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_ctxroipred_data = reshape(permute(mua_ctxroipred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_taskpred_data) | isnan(curr_ctxpred_data) | isnan(curr_ctxroipred_data);
    curr_data(nan_samples) = NaN;
    curr_taskpred_data(nan_samples) = NaN;
    curr_ctxpred_data(nan_samples) = NaN;
    curr_ctxroipred_data(nan_samples) = NaN;

    taskpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_taskpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    ctxroipred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxroipred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
end
figure; hold on;
errorbar(nanmean(taskpred_r2,1),AP_sem(taskpred_r2,1),'b','linewidth',2,'CapSize',0);
errorbar(nanmean(ctxpred_r2,1),AP_sem(ctxpred_r2,1),'color',[0,0.6,0],'linewidth',2,'CapSize',0);
errorbar(nanmean(ctxroipred_r2,1),AP_sem(ctxroipred_r2,1),'color',[1,0.5,0],'linewidth',2,'CapSize',0);
xlabel('Striatum depth');
ylabel('Explained variance');
legend({'Task','Cortex (Full)','Cortex (ROI)'});

% Get significance between cortex kernel and ROI
ctx_kernel_roi_p = nan(n_depths,1);
for curr_depth = 1:n_depths
   ctx_kernel_roi_p(curr_depth) = signrank(ctxroipred_r2(:,curr_depth), ...
       ctxpred_r2(:,curr_depth));
   disp(['Str ' num2str(curr_depth) ' kernel vs ROI: p = ' ...
       num2str(ctx_kernel_roi_p(curr_depth))]);
end



% Get R^2 for task in cortex ROI
taskpred_r2 = nan(max(split_idx),n_depths);
for curr_exp = 1:max(split_idx)
       
    curr_data = reshape(permute(fluor_kernelroi_deconv(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    curr_taskpred_data = reshape(permute(fluor_kernelroi_taskpred(split_idx == curr_exp,:,:),[2,1,3]),[],n_depths);
    
    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_taskpred_data);
    curr_data(nan_samples) = NaN;
    curr_taskpred_data(nan_samples) = NaN;

    taskpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_taskpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));

end
figure; hold on;
errorbar(nanmean(taskpred_r2,1),AP_sem(taskpred_r2,1),'b','linewidth',2,'CapSize',0);
xlabel('Cortex ROI');
ylabel('Task explained variance');


% Cortex explained variance

% (spatial explained variance in pixels)
px_taskpred_r2 = nan(size(U_master,1),size(U_master,2),max(split_idx));
for curr_exp = 1:max(split_idx)  
    px_taskpred_r2(:,:,curr_exp) = AP_spatial_explained_var(U_master(:,:,1:n_vs), ...
        reshape(permute(fluor_allcat_deconv(split_idx == curr_exp,:,:),[2,1,3]),[],n_vs)', ...
        reshape(permute(fluor_taskpred_allcat(split_idx == curr_exp,:,:),[2,1,3]),[],n_vs)',10);
    AP_print_progress_fraction(curr_exp,max(split_idx));
end

figure;imagesc(nanmedian(px_taskpred_r2,3));
axis image off; 
colormap(brewermap([],'Reds'));
caxis([0,1]); 
c = colorbar;
ylabel(c,'Task R^2');
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);


%% Widefield correlation borders

wf_corr_borders_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_borders\wf_corr_borders.mat';
load(wf_corr_borders_fn);

wf_corr_borders_cat = cell2mat(reshape([wf_corr_borders(:).corr_edges],1,1,[]));
figure;
imagesc(nanmean(wf_corr_borders_cat,3));
axis image off;
caxis([0,max(caxis)])
colormap(brewermap([],'Greys'))
ccf_outline = AP_reference_outline('ccf_aligned',[1,0,0]);
cellfun(@(x) set(x,'linewidth',2),vertcat(ccf_outline{:}));


%% Probe location variation
% (plot widefield-estimated/histology-aligned probe location)

animal_days = {...
    'AP032','2018-10-26';
    'AP033','2018-10-26';
    'AP034','2018-10-26';
    'AP036','2018-11-14';
    'AP043','2019-12-09';
    'AP060','2019-12-06';
    'AP061','2019-12-09'};

ccf_fig = figure;
plot_x = ceil(sqrt(length(animal_days)));
plot_y = plot_x;

wf_fig = figure;

for curr_animal_day = 1:length(animal_days)
    
    animal = animal_days{curr_animal_day,1};
    day = animal_days{curr_animal_day,2};
    
    probe_angle = 45; % from horizontal
    
    [img_path,img_exists] = AP_cortexlab_filename(animal,day,[],'imaging');
    avg_im = readNPY([img_path filesep 'meanImage_blue.npy']);
    
    avg_im_aligned = AP_align_widefield(avg_im,animal,day);
    
    h = figure('units','normalized','outerposition',[0 0 1 1]);imagesc(avg_im_aligned);
    axis image off;
    colormap(gray);caxis([0,prctile(avg_im_aligned(:),95)]);
    title('Click probe start/end');
    probe_wf = ginput(2);
    close(h)
    
    % Load and invert master CCF tform
    alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment';
    ccf_tform_fn = [alignment_path filese 'ccf_tform.mat'];
    load(ccf_tform_fn);
    ccf_tform_inverse = invert(ccf_tform);
    
    % Convert aligned widefield probe points to CCF
    % (master alignment to downsampled CCF, need to upsample)
    um2pixel = 20.6;
    probe_ccf_tformed = [probe_wf,[1;1]]*ccf_tform_inverse.T;
    probe_ccf = round(probe_ccf_tformed(:,1:2)*(um2pixel/10));
    
    % Plot probe points on CCF next to average image
    figure(wf_fig);
    
    subplot(length(animal_days),2,((curr_animal_day-1)*2)+1);
    imagesc(avg_im_aligned);
    axis image off;
    colormap(gray);caxis([0,prctile(avg_im_aligned(:),95)]);
    AP_reference_outline('ccf_aligned','r');
    line(probe_wf(:,1),probe_wf(:,2),'color','b','linewidth',1,'linestyle','--');
    title('Widefield');
    
    subplot(length(animal_days),2,((curr_animal_day-1)*2)+2);
    hold on; set(gca,'YDir','reverse');axis image;
    load('C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\cortical_area_boundaries.mat');
    for curr_area_idx =1:length(cortical_area_boundaries)
        p = cellfun(@(outline) plot(outline(:,2),outline(:,1),'color','k'), ...
            cortical_area_boundaries{curr_area_idx},'uni',false);
    end
    line(probe_ccf(:,1),probe_ccf(:,2),'color','r','linewidth',2);
    title('CCF');
    
    % Load in the annotated Allen volume and names
    allen_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
    av = readNPY([allen_path filesep 'annotation_volume_10um_by_index.npy']);
    st = loadStructureTree([allen_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
    
    % Estimate probe entry point by clicked probe start
    depth_start = find(av(probe_ccf(1,2),:,probe_ccf(1,1)) > 1,1);
    probe_entry_ccf = [probe_ccf(1,2),depth_start,probe_ccf(1,1)];
    
    % Estimate height of clicked point from user-set probe angle
    probe_sample_length = pdist2(probe_ccf(1,:),probe_ccf(2,:));
    probe_sample_height = round(probe_sample_length/tand(90-probe_angle));
    probe_air_ccf = [probe_ccf(2,2),depth_start-probe_sample_height,probe_ccf(2,1)];
    
    % Concatenate probe CCF coordinates (in up-down direction);
    probe_ccf = [probe_air_ccf;probe_entry_ccf];
    
    % Get estimated probe vector (widefield)
    r0 = mean(probe_ccf,1);
    xyz = bsxfun(@minus,probe_ccf,r0);
    [~,~,V] = svd(xyz,0);
    probe_direction = V(:,1);
    
    probe_vector_evaluate = [0,sign(probe_direction(2))*1000];
    probe_vector_ccf = round(bsxfun(@plus,bsxfun(@times,probe_vector_evaluate',probe_direction'),r0));
    
    % Get estimated probe vector (histology)
    [probe_filename,probe_filename_exists] = AP_cortexlab_filename(animal,[],[],'probe_ccf');
    load(probe_filename);
    
    use_probe = 1;
    histology_points = probe_ccf(1).points;
    r0 = mean(histology_points,1);
    xyz = bsxfun(@minus,histology_points,r0);
    [~,~,V] = svd(xyz,0);
    histology_probe_direction = V(:,1);
    
    probe_vector_evaluate = [-sign(probe_direction(2))*500,sign(probe_direction(2))*500];
    probe_vector = bsxfun(@plus,bsxfun(@times,probe_vector_evaluate',histology_probe_direction'),r0);
    probe_vector_ccf_histology = round(probe_vector);
    
    % histology done looking A->P so flipped, mirror about bregma
    bregma = allenCCFbregma;
    probe_vector_ccf_histology_mirrored = probe_vector_ccf_histology;
    probe_vector_ccf_histology_mirrored(:,3) = ...
        bregma(3) - (probe_vector_ccf_histology(:,3) - bregma(3));
    
    % Plot estimated probe location on CCF
    % (note the CCF is rotated to allow for dim 1 = x)
    figure(ccf_fig);
    ccf_axes = subplot(plot_y,plot_x,curr_animal_day); hold on
    slice_spacing = 10;
    target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) > 1,[2,1,3]);
    structure_patch = isosurface(target_volume,0);
    structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
    target_structure_color = [0.7,0.7,0.7];
    brain_outline = patch('Vertices',structure_wire.vertices*slice_spacing, ...
        'Faces',structure_wire.faces, ...
        'FaceColor','none','EdgeColor',target_structure_color);
    
    str_id = find(strcmp(st.safe_name,'Caudoputamen'));
    target_volume = permute(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end) == str_id,[2,1,3]);
    structure_patch = isosurface(target_volume,0);
    structure_wire = reducepatch(structure_patch.faces,structure_patch.vertices,0.01);
    target_structure_color = [0.7,0,0.7];
    striatum_outline = patch('Vertices',structure_wire.vertices*slice_spacing, ...
        'Faces',structure_wire.faces, ...
        'FaceColor','none','EdgeColor',target_structure_color);
    
    axis image vis3d off;
    view([-30,25]);
    cameratoolbar(ccf_fig,'SetCoordSys','y');
    cameratoolbar(ccf_fig,'SetMode','orbit');
    title([animal ' ' day]);
    
    % Plot current probe vector
    probe_wf_line = line(ccf_axes,probe_vector_ccf(:,1),probe_vector_ccf(:,2),probe_vector_ccf(:,3),'color','k','linewidth',3);
    % (make sure histology probe enters on the left)
    if  histology_probe_direction(3) < 0
        probe_histology_line = line(ccf_axes, ...
            probe_vector_ccf_histology(:,1), ...
            probe_vector_ccf_histology(:,2), ...
            probe_vector_ccf_histology(:,3),'color','r','linewidth',3);
    else
        probe_histology_line = line(ccf_axes, ...
            probe_vector_ccf_histology_mirrored(:,1), ...
            probe_vector_ccf_histology_mirrored(:,2), ...
            probe_vector_ccf_histology_mirrored(:,3),'color','r','linewidth',3);
    end
    
    drawnow;
    
    if curr_animal_day == 1
        legend([probe_wf_line,probe_histology_line],{'Widefield-estimated','Histology'},'location','nw');
    end
    
end


%% Widefield + Striatum + Cortex ephys


%% ^^^ Correlation between wf/ctx-mua and ctx-mua/str-mua

%%% Load correlation data
use_protocol = 'vanillaChoiceworld';
% use_protocol = 'AP_sparseNoise';

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
data_fn = [data_path filesep 'ctx_fluor_mua_corr_' use_protocol];
load(data_fn);

%%% Get average aligned CSD

% Load CSD 
vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
load(vis_ctx_ephys_fn);

% Concatenate for plotting
animal_cat = cellfun(@(x,y) repmat({x},1,length(y)),{vis_ctx_ephys.animal},{vis_ctx_ephys.day},'uni',false);
animal_cat = horzcat(animal_cat{:});
day_cat = [vis_ctx_ephys(:).day];
n_recordings = length(day_cat);

stim_lfp_t = vis_ctx_ephys(1).stim_lfp_t{1};
stim_lfp_depth = vis_ctx_ephys(1).stim_lfp_depth{1};
stim_csd_depth = vis_ctx_ephys(1).stim_csd_depth{1};

stim_lfp_cat = cell2mat(permute([vis_ctx_ephys(:).stim_lfp],[1,3,2]));
stim_csd_cat = cell2mat(permute([vis_ctx_ephys(:).stim_csd],[1,3,2]));
csd_slice_cat = cell2mat([vis_ctx_ephys(:).csd_slice]);

% Get aligned depth for each recording, plot, save
figure;
curr_recording = 1;
stim_csd_aligned_scaled_cat = nan(0);
for curr_animal = 1:length(vis_ctx_ephys)
    for curr_day = 1:length(vis_ctx_ephys(curr_animal).day)
        
        stim_csd = vis_ctx_ephys(curr_animal).stim_csd{curr_day};
        csd_slice = vis_ctx_ephys(curr_animal).csd_slice{curr_day};
        stim_csd_depth = vis_ctx_ephys(curr_animal).stim_csd_depth{curr_day};
        stim_csd_depth_aligned = vis_ctx_ephys(curr_animal).stim_csd_depth_aligned{curr_day};
                
        % Plot aligned CSD
        depth_align_interp = [-500:20:2000];
        stim_csd_aligned = interp1(...
            stim_csd_depth_aligned,stim_csd_cat(:,:,curr_recording),depth_align_interp, ...
            'linear','extrap');
        
        animal = vis_ctx_ephys(curr_animal).animal;
        day = vis_ctx_ephys(curr_animal).day{curr_day};
        stim_lfp_t = vis_ctx_ephys(curr_animal).stim_lfp_t{curr_day};
        stim_lfp_depth = vis_ctx_ephys(curr_animal).stim_lfp_depth{curr_day};
        stim_lfp = vis_ctx_ephys(curr_animal).stim_lfp{curr_day};
        
        subplot(3,n_recordings,curr_recording);
        imagesc(stim_lfp_t,stim_lfp_depth,stim_lfp);
        caxis([-max(abs(caxis))/2,max(abs(caxis))/2]);
        colormap(brewermap([],'*RdBu'));
        title({animal,day,'LFP'});
        
        subplot(3,n_recordings,size(stim_csd_cat,3)+curr_recording);
        imagesc(stim_lfp_t,stim_csd_depth,stim_csd);
        caxis([-max(abs(caxis))/2,max(abs(caxis))/2]);
        colormap(brewermap([],'*RdBu'));
        title('CSD');
        
        subplot(3,n_recordings,size(stim_csd_cat,3)*2+curr_recording);
        imagesc(stim_lfp_t,stim_csd_depth_aligned,stim_csd_aligned);
        caxis([-max(abs(caxis))/2,max(abs(caxis))/2]);
        colormap(brewermap([],'*RdBu'));
        title('Aligned CSD');
                
        % Keep aligned CSD for averaging
        stim_csd_aligned_scaled_cat(:,:,curr_recording) = stim_csd_aligned./min(csd_slice);       
        
        % Iterate recording index (only used for plotting)
        curr_recording = curr_recording + 1;

    end
end

stim_csd_aligned_mean = nanmean(stim_csd_aligned_scaled_cat,3);

%%% Plot correlations

mua_depth = data(1).cortex_mua_depth{1}; % they're all the same, use 1st
cortex_fluor_corr_cat = cell2mat(horzcat(data.cortex_fluor_corr));
cortex_striatum_corr_cat = cell2mat(permute(horzcat(data.cortex_striatum_corr),[1,3,2]));

figure;

subplot(1,4,1,'YDir','reverse');
imagesc(stim_lfp_t,depth_align_interp,stim_csd_aligned_mean);
caxis([-1,1]);
line([0,0],ylim,'color','k');
colormap(brewermap([],'*RdBu'));
xlabel('Time from stim');
ylabel('Aligned depth');
title('Average aligned CSD')

subplot(1,4,2,'YDir','reverse'); hold on;
plot(cortex_fluor_corr_cat',mua_depth,'color',[0.2,0.8,0.2]);
errorbar(nanmean(cortex_fluor_corr_cat,2), ...
    mua_depth,AP_sem(cortex_fluor_corr_cat,2), ...
    'horizontal','color',[0,0.6,0],'linewidth',2)
xlabel('Correlation');
ylabel('Cortical MUA aligned depth');
title('Cortical fluorescence');

subplot(1,4,3,'YDir','reverse'); hold on;
set(gca,'ColorOrder',copper(4));
plot_str = 1;
plot(permute(cortex_striatum_corr_cat(:,plot_str,:),[1,3,2]),mua_depth,'color',[0.5,0.5,0.5]);
errorbar(nanmean(cortex_striatum_corr_cat(:,plot_str,:),3), ...
    mua_depth,AP_sem(cortex_striatum_corr_cat(:,plot_str,:),3),'k','horizontal','linewidth',2)
xlabel('Correlation');
ylabel('Cortical MUA aligned depth');
title(['Str ' num2str(plot_str) ' multiunit']);

subplot(2,4,4); hold on;
for i = 1:size(cortex_fluor_corr_cat,2)
    plot(cortex_fluor_corr_cat(:,i), ...
        cortex_striatum_corr_cat(:,plot_str,i), ...
        'color',[0.5,0.5,0.5]);
end
plot(nanmean(cortex_fluor_corr_cat,2), ...
    nanmean(cortex_striatum_corr_cat(:,plot_str,:),3),'k','linewidth',2);
xlabel('Fluorescence - cortical MUA correlation');
ylabel('Cortical MUA - striatal MUA correlation')

subplot(2,4,8); hold on;
plot(permute(max(cortex_striatum_corr_cat,[],1),[2,3,1]),'color',[0.5,0.5,0.5]);
errorbar(squeeze(nanmean(max(cortex_striatum_corr_cat,[],1),3)), ...
    squeeze(AP_sem(max(cortex_striatum_corr_cat,[],1),3)),'k','linewidth',2);
xlim([0.5,3.5]);
xlabel('Striatal domain');
ylabel('Cortical MUA max corr');


%% ^^^ Example recording

% AP060 2019-12-06 looks like the best?

animal = 'AP060';
day = '2019-12-06';
experiment = 1;
site = 2; % (cortex)
str_align= 'none'; % (cortex)
verbose = true;
AP_load_experiment;

mua_fig = figure;
raster_fig = figure;
plot_t = [100,300];

% Plot CSD
vis_ctx_ephys_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\vis_ctx_ephys.mat';
load(vis_ctx_ephys_fn);

curr_animal_idx = strcmp(animal,{vis_ctx_ephys.animal});
curr_day_idx = strcmp(day,vis_ctx_ephys(curr_animal_idx).day);

figure;
imagesc(vis_ctx_ephys(curr_animal_idx).stim_lfp_t{curr_day_idx}, ...
    vis_ctx_ephys(curr_animal_idx).stim_csd_depth{curr_day_idx}, ...
    vis_ctx_ephys(curr_animal_idx).stim_csd{curr_day_idx});
caxis([-max(abs(caxis))/2,max(abs(caxis))/2]);
colormap(brewermap([],'*RdBu'));
ylabel('Depth (\mum)');
xlabel('Time from stim');
colorbar;

% (COPIED FROM ABOVE: PLOT CORTEX MULTIUNIT AND FLUORESCENCE)

%%% DEPTH-ALIGN TEMPLATES, FIND CORTEX BOUNDARY
curr_animal_idx = strcmp(animal,{vis_ctx_ephys.animal});
curr_day_idx = strcmp(day,vis_ctx_ephys(curr_animal_idx).day);
curr_csd_depth = vis_ctx_ephys(curr_animal_idx).stim_csd_depth{curr_day_idx};
curr_csd_depth_aligned = vis_ctx_ephys(curr_animal_idx).stim_csd_depth_aligned{curr_day_idx};
template_depths_aligned = interp1(curr_csd_depth,curr_csd_depth_aligned,template_depths);
spike_depths_aligned = interp1(curr_csd_depth,curr_csd_depth_aligned,spike_depths);

% Find cortex end by largest gap between templates
sorted_template_depths = sort([template_depths_aligned]);
[max_gap,max_gap_idx] = max(diff(sorted_template_depths));
ctx_end = sorted_template_depths(max_gap_idx)+1;

ctx_depth = [sorted_template_depths(1),ctx_end];
ctx_units = template_depths_aligned <= ctx_depth(2);

%%% GET FLUORESCENCE AND SPIKES BY DEPTH

% Set binning time
skip_seconds = 60;
spike_binning_t = 1/framerate; % seconds
spike_binning_t_edges = frame_t(1)+skip_seconds:spike_binning_t:frame_t(end)-skip_seconds;
spike_binning_t_centers = spike_binning_t_edges(1:end-1) + diff(spike_binning_t_edges)/2;

% Get fluorescence in pre-drawn ROI
curr_ctx_roi = vis_ctx_ephys(curr_animal_idx).ctx_roi{curr_day_idx};

fVdf_deconv = AP_deconv_wf(fVdf);
fluor_roi = AP_svd_roi(Udf,fVdf_deconv,avg_im,[],curr_ctx_roi);
fluor_roi_interp = interp1(frame_t,fluor_roi,spike_binning_t_centers);

% Set sliding depth window of MUA
depth_corr_range = [-200,1500];
depth_corr_window = 200; % MUA window in microns
depth_corr_window_spacing = 50; % MUA window spacing in microns

depth_corr_bins = [depth_corr_range(1):depth_corr_window_spacing:(depth_corr_range(2)-depth_corr_window); ...
    (depth_corr_range(1):depth_corr_window_spacing:(depth_corr_range(2)-depth_corr_window))+depth_corr_window];
depth_corr_bin_centers = depth_corr_bins(1,:) + diff(depth_corr_bins,[],1)/2;

cortex_mua = zeros(size(depth_corr_bins,2),length(spike_binning_t_centers));
for curr_depth = 1:size(depth_corr_bins,2)
    curr_depth_templates_idx = ...
        find(ctx_units & ...
        template_depths_aligned >= depth_corr_bins(1,curr_depth) & ...
        template_depths_aligned < depth_corr_bins(2,curr_depth));
    
    cortex_mua(curr_depth,:) = histcounts(spike_times_timeline( ...
        ismember(spike_templates,curr_depth_templates_idx)),spike_binning_t_edges);
end

% Plot cortex raster
figure(raster_fig);
subplot(3,1,1); hold on;
plot_t_idx = spike_binning_t_centers >= plot_t(1) & ...
    spike_binning_t_centers <= plot_t(2);
plot(spike_binning_t_centers(plot_t_idx), ...
    fluor_roi_interp(plot_t_idx),'linewidth',2,'color',[0,0.7,0]);

subplot(3,1,2,'YDir','reverse'); hold on;
plot_spikes = spike_times_timeline >= plot_t(1) & ...
    spike_times_timeline <= plot_t(2) & ...
    spike_depths_aligned <= ctx_depth(2);
plot(spike_times_timeline(plot_spikes),spike_depths(plot_spikes),'.k');
ylabel('Depth (\mum)');
xlabel('Time (s)');

%%% LOAD STRIATUM EPHYS AND GET MUA BY DEPTH
clear load_parts
load_parts.ephys = true;
site = 1; % (striatum is always on probe 1)
str_align = 'kernel';
AP_load_experiment;

striatum_mua = nan(n_aligned_depths,length(spike_binning_t_centers));
for curr_depth = 1:n_aligned_depths
    curr_spike_times = spike_times_timeline(aligned_str_depth_group == curr_depth);
    % Skip if no spikes at this depth
    if isempty(curr_spike_times)
        continue
    end
    striatum_mua(curr_depth,:) = histcounts(curr_spike_times,spike_binning_t_edges);
end

% Plot striatum raster
figure(raster_fig);
subplot(3,1,3,'YDir','reverse'); hold on;
plot_spikes = spike_times_timeline >= plot_t(1) & ...
    spike_times_timeline <= plot_t(2) & ...
    spike_depths >= str_depth(1) & spike_depths <= str_depth(2);
plot(spike_times_timeline(plot_spikes),spike_depths(plot_spikes),'.k');
ylabel('Depth (\mum)');
xlabel('Time (s)');
linkaxes(get(raster_fig,'Children'),'x');


% Plot multiunit
figure(mua_fig);
subplot(5,1,1);
plot(spike_binning_t_centers,fluor_roi_interp,'linewidth',2,'color',[0,0.7,0]);
title('Fluorescence');
subplot(5,1,2:4);
imagesc(spike_binning_t_centers,[],cortex_mua)
caxis([0,10]);
colormap(brewermap([],'Greys'));
title('Cortex MUA');
subplot(5,1,5);
imagesc(spike_binning_t_centers,[],striatum_mua);
caxis([0,10]);
title('Striatum MUA');

linkaxes(get(gcf,'Children'),'x');

xlim(plot_t);




% %%%%%% TESTING
% 
% 
% % (wheel velocity)
% wheel_axes = subplot(6,1,6);
% plot_wheel_idx = Timeline.rawDAQTimestamps >= plot_t(1) & ...
%     Timeline.rawDAQTimestamps <= plot_t(2);
% plot(wheel_axes,Timeline.rawDAQTimestamps(plot_wheel_idx), ...
%     wheel_velocity(plot_wheel_idx),'k','linewidth',2);
% ylabel('Wheel velocity');
% axis off
% 
% % (stimuli)
% stim_col = colormap_BlueWhiteRed(5);
% [~,trial_contrast_idx] = ...
%     ismember(trial_conditions(:,1).*trial_conditions(:,2),unique(contrasts'.*sides),'rows');
% stim_lines = arrayfun(@(x) line(wheel_axes,repmat(stimOn_times(x),1,2),ylim(wheel_axes),'color', ...
%     stim_col(trial_contrast_idx(x),:),'linewidth',2), ...
%     find(stimOn_times >= plot_t(1) & stimOn_times <= plot_t(2)));
% 
% % (movement starts)
% move_col = [0.6,0,0.6;0,0.6,0];
% [~,trial_choice_idx] = ismember(trial_conditions(:,3),[-1;1],'rows');
% move_lines = arrayfun(@(x) line(wheel_axes,repmat(wheel_move_time(x),1,2),ylim(wheel_axes),'color', ...
%     move_col(trial_choice_idx(x),:),'linewidth',2), ...
%     find(wheel_move_time >= plot_t(1) & wheel_move_time <= plot_t(2)));
% 
% % (go cues)
% go_col = [0.8,0.8,0.2];
% go_cue_times = signals_events.interactiveOnTimes(1:n_trials);
% go_cue_lines = arrayfun(@(x) line(wheel_axes,repmat(go_cue_times(x),1,2),ylim(wheel_axes),'color', ...
%     go_col,'linewidth',2), ...
%     find(go_cue_times >= plot_t(1) & go_cue_times <= plot_t(2)));
% 
% % (outcomes)
% outcome_col = [0,0,0.8;0.5,0.5,0.5];
% reward_lines = arrayfun(@(x) line(wheel_axes,repmat(reward_t_timeline(x),1,2),ylim(wheel_axes),'color', ...
%     outcome_col(1,:),'linewidth',2), ...
%     find(reward_t_timeline >= plot_t(1) & reward_t_timeline <= plot_t(2)));
% punish_times = signals_events.responseTimes(trial_outcome == -1);
% punish_lines = arrayfun(@(x) line(wheel_axes,repmat(punish_times(x),1,2),ylim(wheel_axes),'color', ...
%     outcome_col(2,:),'linewidth',2), ...
%     find(punish_times >= plot_t(1) & punish_times <= plot_t(2)));



%%  Cortical muscimol


%% ^^^ Cortical spike rate pre/post muscimol
% (use spike rate over all experiments pre/post muscimol)

animal_days = { ...
    'AP052','2019-09-20';
    'AP058','2019-12-06'};

figure;

for curr_animalday = 1:length(animal_days)
    
    animal = animal_days{curr_animalday,1};
    day = animal_days{curr_animalday,2};
    
    % Load data (first experiment - but spikes throughout used)
    experiment = 1;
    AP_load_experiment
    
    % Estimate boundaries of cortex (the dumb way: first template/gap)
    sorted_template_depths = sort(template_depths);
    ctx_start = sorted_template_depths(1) - 1;
    [max_gap,max_gap_idx] = max(diff(sorted_template_depths));
    ctx_end = sorted_template_depths(max_gap_idx)+1;
    ctx_depth = [ctx_start,ctx_end];
    ctx_templates = template_depths <= ctx_depth(2);
    
    % Set experiments in conditions (1-2 = pre-muscimol, 3-4 = post-muscimol)
    cond_expts = {[1,2],[3,4]};
    
    spike_rate_cond = nan(max(spike_templates),2);
    for curr_cond = 1:2
        
        exp_starts = sync(2).timestamps(sync(2).values == 1);
        exp_stops = sync(2).timestamps(sync(2).values == 0);
        
        curr_exp_start = exp_starts(cond_expts{curr_cond}(1));
        curr_exp_stop = exp_stops(cond_expts{curr_cond}(end));
        
        curr_use_spikes = spike_times >= curr_exp_start & ...
            spike_times <= curr_exp_stop;
        
        spike_rate_cond(:,curr_cond) = ...
            accumarray(spike_templates(curr_use_spikes),1,[max(spike_templates),1])./ ...
            (curr_exp_stop - curr_exp_start);
        
    end
    
    spike_rate_change = (spike_rate_cond(:,2) - spike_rate_cond(:,1))./(spike_rate_cond(:,1)+spike_rate_cond(:,2));
    
    subplot(length(animal_days),2,(curr_animalday-1)*2+1,'YDir','reverse'); hold on;
    plot(reshape([spike_rate_cond,nan(size(template_depths))]',[],1), ...
        reshape(repmat(template_depths,1,3)',[],1),'r');
    p1 = plot(spike_rate_cond(:,1),template_depths,'.k','MarkerSize',10);
    p2 = plot(spike_rate_cond(:,2),template_depths,'.r','MarkerSize',10);
    xlabel('Spikes/s')
    ylabel('Depth (\mum)');
    legend([p1,p2],{'Pre-muscimol','Post-muscimol'});
    axis tight;
    xlim([-1,prctile(spike_rate_cond(:,1),95)])
   
    subplot(length(animal_days),2,(curr_animalday-1)*2+2); hold on;
    plot(spike_rate_change,template_depths,'.k','MarkerSize',10);   
    axis tight;
    xlim([-1.1,1.1]);
    line([0,0],ylim);
    set(gca,'YDir','reverse');
    xlabel('(Post-pre)/(pre+post)');
    ylabel('Depth (\mum)');
    title({animal,day,'Change'});
    
    
end






%% ^^^ VFS pre/post musicmol

muscimol_wf_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\muscimol_wf.mat';
load(muscimol_wf_fn);

std_cat = vertcat(muscimol_wf.std);
vfs_cat = vertcat(muscimol_wf.vfs);

std_change = (cat(3,std_cat{:,2})-cat(3,std_cat{:,1}))./(cat(3,std_cat{:,2})+cat(3,std_cat{:,1}));
vfs_softnorm = 0.2;
vfs_change = (abs(cat(3,vfs_cat{:,2}))-abs(cat(3,vfs_cat{:,1})))./(vfs_softnorm+abs(cat(3,vfs_cat{:,2}))+abs(cat(3,vfs_cat{:,1})));

figure; 

% Plot std
subplot(2,3,1);
imagesc(nanmean(cat(3,std_cat{:,1}),3));
axis image off;
caxis([0,0.03]);
colormap(gca,brewermap([],'*Greys'));
AP_reference_outline('ccf_aligned','r');
title('Std (pre-muscimol');
colorbar

subplot(2,3,2);
imagesc(nanmean(cat(3,std_cat{:,2}),3));
axis image off;
caxis([0,0.03]);
colormap(gca,brewermap([],'*Greys'));
AP_reference_outline('ccf_aligned','r');
title('Std (post-muscimol');
colorbar

subplot(2,3,3);
imagesc(nanmean(std_change,3));
axis image off;
caxis([-0.5,0.5]);
colormap(gca,brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('(post-pre)/(post+pre)');
colorbar

% Plot VFS
subplot(2,3,4);
imagesc(nanmean(cat(3,vfs_cat{:,1}),3));
axis image off;
caxis([-1,1]);
colormap(gca,brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('VFS pre-muscimol');
colorbar

subplot(2,3,5);
imagesc(nanmean(cat(3,vfs_cat{:,2}),3));
axis image off;
caxis([-1,1]);
colormap(gca,brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('VFS post-muscimol');
colorbar

subplot(2,3,6);
imagesc(nanmean(vfs_change,3));
axis image off;
caxis([-0.5,0.5]);
colormap(gca,brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('(post-pre)/(post+pre)');
colorbar


%% ^^^ Striatal spike rate pre/post muscimol

animals = {'AP045','AP054','AP055','AP053','AP047','AP048'};

spike_rate_cond = cell(size(animals));

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworldNoRepeats';
    experiments = AP_find_experiments(animal,protocol);
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
        
    disp(animal);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        
        % Load data (first experiment - but spikes throughout used)
        experiment = experiments(curr_day).experiment(1);
        load_parts.ephys = true;
        AP_load_experiment
        
        % Set experiments in conditions 
        % (1-3 pre-muscimol, 4+ post-muscimol)
        % (assumes all repeated expts/failures were post-muscimol)
        curr_experiments = AP_list_experiments(animal,day);
        cond_expts = {[1:3],[4:length(curr_experiments)]};
        
        for curr_cond = 1:2
                exp_starts = sync(2).timestamps(sync(2).values == 1);
                exp_stops = sync(2).timestamps(sync(2).values == 0);
                
                curr_exp_start = exp_starts(cond_expts{curr_cond}(1));
                curr_exp_stop = exp_stops(cond_expts{curr_cond}(end));
                
                curr_use_spikes = spike_times >= curr_exp_start & ...
                    spike_times <= curr_exp_stop & ~isnan(aligned_str_depth_group);
                
                spike_rate_cond{curr_animal}{curr_day}(:,curr_cond) = ...
                    accumarray(aligned_str_depth_group(curr_use_spikes),1, ...
                    [n_aligned_depths,1])./(curr_exp_stop - curr_exp_start);
        end
         
        clearvars -except animals curr_animal animal experiments curr_day ...
            spike_rate_cond
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
    
end

% Concatenate across recordings and get change
spike_rate_cond_cat = cell2mat(reshape([spike_rate_cond{:}],1,1,[]));
spike_rate_cond_cat_change = ...
    squeeze((spike_rate_cond_cat(:,1,:) - spike_rate_cond_cat(:,2,:))./ ...
    (spike_rate_cond_cat(:,1,:) + spike_rate_cond_cat(:,2,:)));
n_depths = size(spike_rate_cond_cat,1);

figure;
set(gca,'YDir','reverse'); hold on;
plot(spike_rate_cond_cat_change,1:n_depths,'color',[0.5,0.5,0.5]);
errorbar(nanmean(spike_rate_cond_cat_change,2),1:n_depths, ...
    AP_sem(spike_rate_cond_cat_change,2),'k','horizontal','linewidth',2);
line([0,0],ylim,'color','r')
xlim([-1.1,1.1]);
ylim([1-0.2,n_depths+0.2])
xlabel('(post-pre)/(post+pre)');
ylabel('Striatal depth');
title('Muscimol change');


%% ^^^ Striatum cortical kernels pre/post muscimol
error('I don''t think this should be included');

protocols = {'vanillaChoiceworldNoRepeats_pre_muscimol','vanillaChoiceworldNoRepeats_post_muscimol'};

for protocol = protocols 
    protocol = cell2mat(protocol);
    
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
    k_fn = [data_path filesep 'ctx_str_kernels_' protocol];
    load(k_fn);
    
    framerate = 35;
    upsample_factor = 1;
    sample_rate = framerate*upsample_factor;
    kernel_t = [-0.1,0.1];
    kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
    t = kernel_frames/sample_rate;
    
    % Concatenate explained variance
    expl_var_animal = cell2mat(cellfun(@(x) nanmean(horzcat(x{:}),2),ctx_str_expl_var','uni',false));
    figure('Name',protocol);
    errorbar(nanmean(expl_var_animal,2),AP_sem(expl_var_animal,2),'k','linewidth',2);
    xlabel('Striatal depth');
    ylabel('Fraction explained variance');
    
    % Concatenate and mean
    % (kernel is -:+ fluorescence lag, flip to be spike-oriented)
    k_px_timeflipped = cellfun(@(x) cellfun(@(x) x(:,:,end:-1:1,:),x,'uni',false),ctx_str_kernel,'uni',false);
    k_px_animal = cellfun(@(x) nanmean(cat(5,x{:}),5),k_px_timeflipped,'uni',false);
    k_px = nanmean(double(cat(5,k_px_animal{:})),5);
    
    % Get center-of-mass maps
    k_px_positive = k_px;
    k_px_positive(k_px_positive < 0) = 0;
    k_px_com = sum(k_px_positive.*permute(1:n_aligned_depths,[1,3,4,2]),4)./sum(k_px_positive,4);
    k_px_com_colored = nan(size(k_px_com,1),size(k_px_com,2),3,size(k_px_com,3));
    
    use_colormap = min(jet(255)-0.2,1);
    for curr_frame = 1:size(k_px_com,3)
        k_px_com_colored(:,:,:,curr_frame) = ...
            ind2rgb(round(mat2gray(k_px_com(:,:,curr_frame),...
            [1,n_aligned_depths])*size(use_colormap,1)),use_colormap);
    end
    
    % Plot center kernel frames independently at t = 0
    figure('Name',protocol);
    plot_frame = kernel_frames == 0;
    for curr_depth = 1:n_aligned_depths
       subplot(n_aligned_depths,1,curr_depth);
       imagesc(k_px(:,:,plot_frame,curr_depth));
       AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
       axis image off;
       colormap(brewermap([],'*RdBu'));
       caxis([-0.01,0.01]);
    end
    
    % Plot center-of-mass color at select time points 
    plot_t = [-0.05:0.025:0.05];
    
    k_px_com_colored_t = ...
        permute(reshape(interp1(t,permute(reshape(k_px_com_colored,[],3,length(t)), ...
        [3,1,2]),plot_t),length(plot_t),size(k_px_com_colored,1), ...
        size(k_px_com_colored,2),3),[2,3,4,1]);
    
    k_px_max = squeeze(max(k_px,[],4));
    k_px_max_t = ...
        permute(reshape(interp1(t,reshape(k_px_max,[],length(t))', ...
        plot_t),length(plot_t),size(k_px_max,1), ...
        size(k_px_max,2)),[2,3,1]);
    
    weight_max = 0.005;
    figure('Name',protocol);
    for t_idx = 1:length(plot_t)
        subplot(1,length(plot_t),t_idx);
        p = image(k_px_com_colored_t(:,:,:,t_idx));
        set(p,'AlphaData', ...
            mat2gray(k_px_max_t(:,:,t_idx),[0,weight_max]));
        axis image off;
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        title([num2str(plot_t(t_idx)),' s']);
    end    
        
    % Plot movie of kernels
    AP_image_scroll(reshape(permute(k_px,[1,4,2,3]),size(k_px,1)*size(k_px,4),size(k_px,2),length(t)),t);
    colormap(brewermap([],'*RdBu'));
    caxis([-max(caxis),max(caxis)]);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5],[],[size(k_px,1),size(k_px,2),size(k_px,4),1]);
    axis image off
    
    drawnow;
    
end


%% ^^^ Cortex/striatum passive gratings pre/post muscimol

data_fns = { ...
    'trial_activity_AP_lcrGratingPassive_pre_muscimol', ...
    'trial_activity_AP_lcrGratingPassive_post_muscimol'};

mua_prepost_norm = cell(2,1);
fluor_kernelroi_prepost_norm = cell(2,1);

stimIDs = cell(2,1);
mua_muscimol = cell(2,1);
mua_ctxpred_muscimol = cell(2,1);
fluor_muscimol = cell(2,1);
fluor_roi_muscimol = cell(2,1);
fluor_kernelroi_muscimol = cell(2,1);

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;

    % Split by experiment
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    
    % Get trials with movement during stim to exclude
    wheel_thresh = 0.025;
    quiescent_trials = ~any(abs(wheel_allcat(:,t >= -0.5 & t <= 0.5)) > wheel_thresh,2);
    quiescent_trials_exp = mat2cell(quiescent_trials,trials_recording,1);
    
    % Get stim and activity by experiment
    trial_stim_allcat_exp = mat2cell(trial_stim_allcat,trials_recording,1);  
    fluor_deconv_exp = mat2cell(fluor_allcat_deconv,trials_recording,length(t),n_vs);
    fluor_roi_deconv_exp = mat2cell(fluor_roi_deconv,trials_recording,length(t),numel(wf_roi));
    fluor_kernelroi_deconv_exp = mat2cell(fluor_kernelroi_deconv,trials_recording,length(t),n_depths);
    mua_allcat_exp = mat2cell(mua_allcat,trials_recording,length(t),n_depths);
    mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,trials_recording,length(t),n_depths);
            
    % (use the norm from pre-muscimol)     
    mua_prepost_norm{curr_data} = vertcat(mua_norm{:});
    mua_allcat_exp = cellfun(@(act,curr_norm,use_norm) ...
        (act.*curr_norm)./use_norm,mua_allcat_exp, ...
        vertcat(mua_norm{:}),mua_prepost_norm{1},'uni',false);
    mua_ctxpred_allcat_exp = cellfun(@(act,curr_norm,use_norm) ...
        (act.*curr_norm)./use_norm,mua_ctxpred_allcat_exp, ...
        vertcat(mua_norm{:}),mua_prepost_norm{1},'uni',false);
    
    fluor_kernelroi_prepost_norm{curr_data} = fluor_kernelroi_norm;
    fluor_kernelroi_deconv_exp = cellfun(@(act,curr_norm,use_norm) ...
        (act.*curr_norm)./use_norm,fluor_kernelroi_deconv_exp, ...
        fluor_kernelroi_norm,fluor_kernelroi_prepost_norm{1},'uni',false);
    
    % Exclude trials with fluorescence spikes
    % (this is a dirty way to do this but don't have a better alt)
    fluor_spike_thresh = 100;
    fluor_spike_trial = cellfun(@(x) any(any(x > fluor_spike_thresh,2),3), ...
        fluor_kernelroi_deconv_exp,'uni',false);
    
    % Grab data
    stimIDs{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        trial_stim_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    mua_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        mua_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    mua_ctxpred_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        mua_ctxpred_allcat_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_roi_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_roi_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
    fluor_kernelroi_muscimol{curr_data} = cellfun(@(data,quiesc,fluor_spike) data(quiesc & ~fluor_spike,:,:), ...
        fluor_kernelroi_deconv_exp,quiescent_trials_exp,fluor_spike_trial,'uni',false);
    
end

% Plot average fluorescence
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');

use_stim = 1;

fluor_premuscimol_mean = ...
    svdFrameReconstruct(U_master(:,:,1:200), ...
    permute(nanmean(cell2mat(cellfun(@(x,stim) ...
    nanmean(x(stim == use_stim,:,:),1),fluor_muscimol{1},stimIDs{1},'uni',false)),1),[3,2,1]));
fluor_postmuscimol_mean = ...
    svdFrameReconstruct(U_master(:,:,1:200), ...
    permute(nanmean(cell2mat(cellfun(@(x,stim) ...
    nanmean(x(stim == use_stim,:,:),1),fluor_muscimol{2},stimIDs{2},'uni',false)),1),[3,2,1]));

AP_image_scroll([fluor_premuscimol_mean,fluor_postmuscimol_mean]);
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5],[],[size(U_master,1),size(U_master,2),1,2]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
axis image;
colormap(brewermap([],'*RdBu'));


figure;
t_stim = t >= 0.05 & t <= 0.15;

subplot(1,2,1)
imagesc(nanmean(fluor_premuscimol_mean(:,:,t_stim),3));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
c = caxis;
axis image off;
colormap(brewermap([],'*RdBu'));
title('Pre-muscimol');

subplot(1,2,2)
imagesc(nanmean(fluor_postmuscimol_mean(:,:,t_stim),3));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
caxis(c);
axis image off;
colormap(brewermap([],'*RdBu'));
title('Post-muscimol');

% Get pre/post stim response
use_stim = 1;

mua_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_muscimol{1},stimIDs{1},'uni',false));
mua_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_muscimol{2},stimIDs{2},'uni',false));

mua_ctxpred_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_ctxpred_muscimol{1},stimIDs{1},'uni',false));
mua_ctxpred_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),mua_ctxpred_muscimol{2},stimIDs{2},'uni',false));

fluor_roi_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_roi_muscimol{1},stimIDs{1},'uni',false));
fluor_roi_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_roi_muscimol{2},stimIDs{2},'uni',false));

fluor_kernelroi_premuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_kernelroi_muscimol{1},stimIDs{1},'uni',false));
fluor_kernelroi_postmuscimol_mean = ...
    cell2mat(cellfun(@(x,stim) nanmean(x(stim == use_stim,:,:),1),fluor_kernelroi_muscimol{2},stimIDs{2},'uni',false));

% Plot all str responses
figure;
for curr_str = 1:n_depths
    subplot(n_depths,1,curr_str);
    AP_errorfill(t,nanmean(mua_premuscimol_mean(:,:,curr_str),1)', ...
        AP_sem(mua_premuscimol_mean(:,:,curr_str),1)','k',1,false);
    AP_errorfill(t,nanmean(mua_postmuscimol_mean(:,:,curr_str),1)', ...
        AP_sem(mua_postmuscimol_mean(:,:,curr_str),1)','r',1,false);
    xlabel('Time from stim (s)');
    ylabel('Spikes (std)');
    title(['Str ' num2str(curr_str)]);
end
linkaxes(get(gcf,'Children'))

% Plot ctx responses
figure;
plot_ctx = [1,3,7];
for curr_ctx_idx = 1:length(plot_ctx)
    curr_ctx = plot_ctx(curr_ctx_idx);
    subplot(length(plot_ctx),1,curr_ctx_idx);
    AP_errorfill(t,nanmean(fluor_roi_premuscimol_mean(:,:,curr_ctx),1)', ...
        AP_sem(fluor_roi_premuscimol_mean(:,:,curr_ctx),1)','k',1,false);
    AP_errorfill(t,nanmean(fluor_roi_postmuscimol_mean(:,:,curr_ctx),1)', ...
        AP_sem(fluor_roi_postmuscimol_mean(:,:,curr_ctx),1),'r',1,false);
    ylabel(wf_roi(curr_ctx).area);
end
linkaxes(get(gcf,'Children'))

% Plot pre/post muscimol and repsonse change for pair of str/ctx
figure;
for plot_str = 1:n_depths
    
    t_stim = t >= 0.05 & t <= 0.15;
    mua_avg_premuscimol = permute(nanmean(mua_premuscimol_mean(:,t_stim,:),2),[1,3,2]);
    mua_avg_postmuscimol = permute(nanmean(mua_postmuscimol_mean(:,t_stim,:),2),[1,3,2]);
    mua_avg_postpre_change = (mua_avg_postmuscimol-mua_avg_premuscimol)./(mua_avg_premuscimol+mua_avg_postmuscimol);
    
    fluor_avg_premuscimol = permute(nanmean(fluor_kernelroi_premuscimol_mean(:,t_stim,:),2),[1,3,2]);
    fluor_avg_postmuscimol = permute(nanmean(fluor_kernelroi_postmuscimol_mean(:,t_stim,:),2),[1,3,2]);
    fluor_avg_postpre_change = (fluor_avg_postmuscimol-fluor_avg_premuscimol)./(fluor_avg_premuscimol+fluor_avg_postmuscimol);
    
    subplot(n_depths,3,(plot_str-1)*n_depths+1);hold on;
    AP_errorfill(t,nanmean(mua_premuscimol_mean(:,:,plot_str),1)', ...
        AP_sem(mua_premuscimol_mean(:,:,plot_str),1)','k');
    AP_errorfill(t,nanmean(mua_postmuscimol_mean(:,:,plot_str),1)', ...
        AP_sem(mua_postmuscimol_mean(:,:,plot_str),1)','r');
    xlim([-0.2,1])
    xlabel('Time from stim (s)')
    ylabel(['Str ' num2str(plot_str)]);
    axis square
    
    subplot(n_depths,3,(plot_str-1)*n_depths+2);hold on;
    AP_errorfill(t,nanmean(fluor_kernelroi_premuscimol_mean(:,:,plot_str),1)', ...
        AP_sem(fluor_kernelroi_premuscimol_mean(:,:,plot_str),1)','k');
    AP_errorfill(t,nanmean(fluor_kernelroi_postmuscimol_mean(:,:,plot_str),1)', ...
        AP_sem(fluor_kernelroi_postmuscimol_mean(:,:,plot_str),1)','r');
    xlim([-0.2,1])
    xlabel('Time from stim (s)')
    ylabel('Cortex ROI');
    axis square
    
    subplot(n_depths,3,(plot_str-1)*n_depths+3);
    plot(fluor_avg_postpre_change(:,plot_str),mua_avg_postpre_change(:,plot_str),'.k','MarkerSize',20)
    xlabel(['Cortex ROI (post-pre)']);
    ylabel(['Str ' num2str(plot_str) ' (post-pre)/(pre+post)']);
    line([-1,1],[0,0],'color','k','linestyle','--');
    line([0,0],[-1,1],'color','k','linestyle','--');
    line([-1,1],[-1,1],'color','k');
    % xlim([-1,1]);
    % ylim([-1,1]);
    axis square;
    
    nonan_points = ~isnan(mua_avg_postpre_change(:,plot_str)) & ...
        ~isnan(fluor_avg_postpre_change(:,plot_str));
    [r,p] = corrcoef(fluor_avg_postpre_change(nonan_points,plot_str), ...
        mua_avg_postpre_change(nonan_points,plot_str));
    
end

error('Need to fix this metric - also plot for 2 (and 3?)');


%% ^^^ Striatal task trial activity pre/post muscimol

data_fns = { ...
    'trial_activity_vanillaChoiceworldNoRepeats_pre_muscimol', ...
    'trial_activity_vanillaChoiceworldNoRepeats_post_muscimol'};

mua_prepost_norm = cell(2,1);

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Split data by recording
    trials_allcat = size(wheel_allcat,1);
    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
    use_split = trials_recording;
    split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
        [1:length(use_split)]',reshape(use_split,[],1),'uni',false));
    
    % (use the norm from pre-muscimol)
    mua_prepost_norm{curr_data} = vertcat(mua_norm{:});
    mua_allcat = cell2mat(cellfun(@(act,curr_norm,use_norm) ...
        (act.*curr_norm)./use_norm, ...
        mat2cell(mua_allcat,trials_recording,length(t),n_depths), ...
        vertcat(mua_norm{:}),mua_prepost_norm{1},'uni',false));   
    mua_ctxpred_allcat = cell2mat(cellfun(@(act,curr_norm,use_norm) ...
        (act.*curr_norm)./use_norm, ...
        mat2cell(mua_ctxpred_allcat,trials_recording,length(t),n_depths), ...
        vertcat(mua_norm{:}),mua_prepost_norm{1},'uni',false));    
    mua_taskpred_allcat = cell2mat(cellfun(@(act,curr_norm,use_norm) ...
        (act.*curr_norm)./use_norm, ...
        mat2cell(mua_taskpred_allcat,trials_recording,length(t),n_depths), ...
        vertcat(mua_norm{:}),mua_prepost_norm{1},'uni',false));
        
    % Plot stim-aligned/sorted measured and predicted striatum activity
    % (correct contra trials)      
    switch curr_data
        case 1
            cond_name = 'Pre-muscimol';
        case 2
            cond_name = 'Post-muscimol';
    end
    
    for curr_trial_set = 1:2
        switch curr_trial_set
            case 1
                plot_trials = move_t < Inf & trial_stim_allcat > 0 & trial_choice_allcat == -1;
                figure('Name',[cond_name ': Correct contra trials']);
            case 2
                plot_trials = move_t < Inf & trial_stim_allcat < 0 & trial_choice_allcat == 1;
                figure('Name',[cond_name ': Correct ipsi trials']);
        end
        
        p = gobjects(n_depths,4);
        colormap(brewermap([],'Greys'));
        for curr_depth = 1:n_depths
            
            % Get trials to plot, sort by reaction time
            curr_trials = plot_trials & ~all(isnan(mua_allcat(:,:,curr_depth)),2);
            curr_trials_idx = find(curr_trials);
            [~,rxn_sort_idx] = sort(move_t(curr_trials_idx));
            
            sorted_plot_trials = curr_trials_idx(rxn_sort_idx);
            
            curr_plot = mua_allcat(sorted_plot_trials,:,curr_depth);
            curr_taskpred_plot = mua_taskpred_allcat(sorted_plot_trials,:,curr_depth);
            curr_ctxpred_plot = mua_ctxpred_allcat(sorted_plot_trials,:,curr_depth);
            
            % Smooth and plot with stim/move/reward times
            % (as conv(nans-zeroed)./conv(non-nan) to ignore in nans in conv)
            smooth_filt = [50,1]; % (trials x frames)
            
            curr_plot_smooth = conv2(curr_plot,ones(smooth_filt),'same')./ ...
                conv2(~isnan(curr_plot),ones(smooth_filt),'same');
            
            curr_taskpred_plot_smooth = curr_taskpred_plot;
            curr_taskpred_plot_smooth(isnan(curr_taskpred_plot_smooth)) = 0;
            curr_taskpred_plot_smooth = conv2(curr_taskpred_plot_smooth,ones(smooth_filt),'same')./ ...
                conv2(~isnan(curr_taskpred_plot),ones(smooth_filt),'same');
            
            curr_ctxpred_plot_smooth = curr_ctxpred_plot;
            curr_ctxpred_plot_smooth(isnan(curr_ctxpred_plot_smooth)) = 0;
            curr_ctxpred_plot_smooth = conv2(curr_ctxpred_plot_smooth,ones(smooth_filt),'same')./ ...
                conv2(~isnan(curr_ctxpred_plot),ones(smooth_filt),'same');
            
            p(curr_depth,1) = subplot(n_depths,4,1+(curr_depth-1)*4,'YDir','reverse'); hold on
            imagesc(t,[],curr_plot_smooth);
            plot(zeros(size(curr_trials_idx)),1:length(curr_trials_idx),'MarkerSize',1,'color','r');
            plot(move_t(sorted_plot_trials),1:length(curr_trials_idx),'MarkerSize',1,'color',[0.8,0,0.8]);
            %         plot(outcome_t(sorted_plot_trials),1:length(curr_trials_idx),'.','MarkerSize',1,'color','b');
            axis tight;
            xlim([-0.2,1]);
            xlabel('Time from stim');
            ylabel('Trials (rxn sorted)');
            title('Measured');
            
            p(curr_depth,2) = subplot(n_depths,4,2+(curr_depth-1)*4,'YDir','reverse'); hold on
            imagesc(t,[],curr_taskpred_plot_smooth);
            plot(zeros(size(curr_trials_idx)),1:length(curr_trials_idx),'MarkerSize',1,'color','r');
            plot(move_t(sorted_plot_trials),1:length(curr_trials_idx),'MarkerSize',1,'color',[0.8,0,0.8]);
            %         plot(outcome_t(sorted_plot_trials),1:length(curr_trials_idx),'.','MarkerSize',1,'color','b');
            axis tight;
            xlim([-0.2,1]);
            xlabel('Time from stim');
            ylabel('Trials (rxn sorted)');
            title('Task-predicted');
            
            p(curr_depth,3) = subplot(n_depths,4,3+(curr_depth-1)*4,'YDir','reverse'); hold on
            imagesc(t,[],curr_ctxpred_plot_smooth);
            plot(zeros(size(curr_trials_idx)),1:length(curr_trials_idx),'MarkerSize',1,'color','r');
            plot(move_t(sorted_plot_trials),1:length(curr_trials_idx),'MarkerSize',1,'color',[0.8,0,0.8]);
            %         plot(outcome_t(sorted_plot_trials),1:length(curr_trials_idx),'.','MarkerSize',1,'color','b');
            axis tight;
            xlim([-0.2,1]);
            xlabel('Time from stim');
            ylabel('Trials (rxn sorted)');
            title('Cortex-predicted');
            
            % Split and average trials by animal
            curr_trials_exp = mat2cell(curr_trials,use_split,1);
            
            curr_mua_exp = mat2cell(mua_allcat(:,:,curr_depth),use_split,length(t));
            curr_mua_exp_mean = cell2mat(cellfun(@(data,trials) nanmean(data(trials,:),1), ...
                curr_mua_exp,curr_trials_exp,'uni',false));
            
            curr_taskpred_mua_exp = mat2cell(mua_taskpred_allcat(:,:,curr_depth),use_split,length(t));
            curr_taskpred_mua_exp_mean = cell2mat(cellfun(@(data,trials) nanmean(data(trials,:),1), ...
                curr_taskpred_mua_exp,curr_trials_exp,'uni',false));
            
            curr_ctxpred_mua_exp = mat2cell(mua_ctxpred_allcat(:,:,curr_depth),use_split,length(t));
            curr_ctxpred_mua_exp_mean = cell2mat(cellfun(@(data,trials) nanmean(data(trials,:),1), ...
                curr_ctxpred_mua_exp,curr_trials_exp,'uni',false));
            
            % Plot PSTH (measured, task-predicted, cortex-predicted);
            p(curr_depth,4) = subplot(n_depths,4,4+(curr_depth-1)*4); hold on
            p1 = AP_errorfill(t,nanmean(curr_mua_exp_mean,1)', ...
                AP_sem(curr_mua_exp_mean,1)','k',0.5);
            p2 = AP_errorfill(t,nanmean(curr_taskpred_mua_exp_mean,1)', ...
                AP_sem(curr_taskpred_mua_exp_mean,1)',[0,0,0.7],0.5);
            p3 = AP_errorfill(t,nanmean(curr_ctxpred_mua_exp_mean,1)', ...
                AP_sem(curr_ctxpred_mua_exp_mean,1)',[0,0.7,0],0.5);
            xlim([-0.2,1])
            line([0,0],ylim,'color','r');
            line(repmat(median(move_t(sorted_plot_trials)),1,2),ylim,'color',[0.8,0,0.8],'linestyle','--');
            line(repmat(median(outcome_t(sorted_plot_trials)),1,2),ylim,'color','b','linestyle','--');
            xlabel('Time from stim');
            ylabel('Spikes (std)');
            legend([p1,p2,p3],{'Measured','Task-predicted','Cortex-predicted'});
            
        end
        % Link the x-axes, set the c/y-axes same within a row
        linkaxes(p(:),'x');
        
        for curr_row = 1:size(p,1)
            curr_ylim = ylim(p(curr_row,4));
            caxis(p(curr_row,1),[0,curr_ylim(2)]);
            caxis(p(curr_row,2),[0,curr_ylim(2)]);
            caxis(p(curr_row,3),[0,curr_ylim(2)]);
        end
        
        trial_scale = 500;
        t_scale = 0.5;
        y_scale = 1;
        line(p(1,1),min(xlim(p(1,1))) + [0,t_scale],repmat(min(ylim(p(1,1))),2,1),'color','k','linewidth',3);
        line(p(1,4),min(xlim(p(1,4))) + [0,t_scale],repmat(min(ylim(p(1,4))),2,1),'color','k','linewidth',3);
        line(p(1,1),repmat(min(xlim(p(1,1))),2,1),min(ylim(p(1,1))) + [0,trial_scale],'color','k','linewidth',3);
        line(p(1,4),repmat(min(xlim(p(1,4))),2,1),min(ylim(p(1,4))) + [0,y_scale],'color','k','linewidth',3);
        drawnow;
        
    end
    
end


%% ^^^ Striatal task kernels pre/post muscimol

data_fns = { ...
    'trial_activity_vanillaChoiceworldNoRepeats_pre_muscimol', ...
    'trial_activity_vanillaChoiceworldNoRepeats_post_muscimol'};

mua_norm_exp = cell(2,1);
task_str_kernel = cell(2,1);
task_str_ctxpred_kernel = cell(2,1);
task_ctx_roi_kernel = cell(2,1);

for curr_data = 1:length(data_fns)
    
    % Load data
    data_fn = data_fns{curr_data};
    AP_load_concat_normalize_ctx_str;
    
    % Get average task > cortex kernels (ROIs)
    n_regressors = length(task_regressor_labels);
    regressor_roi = cell(n_regressors,1);
    for curr_regressor = 1:n_regressors
        
        curr_k = cell2mat(cellfun(@(x) x{curr_regressor}, ...
            permute(vertcat(fluor_taskpred_k_all{:}),[2,3,4,1]),'uni',false));
        
        curr_k_roi = nan(size(curr_k,1),size(curr_k,2),n_rois,size(curr_k,4));
        for curr_subregressor = 1:size(curr_k,1)
            for curr_exp = 1:size(curr_k,4)
                curr_k_roi(curr_subregressor,:,:,curr_exp) = ...
                    permute(AP_svd_roi(U_master(:,:,1:n_vs), ...
                    permute(curr_k(curr_subregressor,:,:,curr_exp),[3,2,1]), ...
                    [],[],cat(3,wf_roi.mask)),[3,2,1]);
            end
        end
        
        regressor_roi{curr_regressor} = curr_k_roi;
        
        AP_print_progress_fraction(curr_regressor,n_regressors);
    end
    
    % Keep task > str kernels, task > ctx-pred str norm factor
    mua_norm_exp{curr_data} = vertcat(mua_norm{:});
    task_str_kernel{curr_data} = vertcat(mua_taskpred_k_all{:});
    task_str_ctxpred_kernel{curr_data} = vertcat(mua_ctxpred_taskpred_k_all{:});
    task_ctx_roi_kernel{curr_data} = regressor_roi;
    
end

% Normalize and concatenate task kernels
mua_task_k_premuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_kernel{1},mua_norm_exp{1},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

mua_ctxpred_task_k_premuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_ctxpred_kernel{1},mua_norm_exp{1},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

mua_task_k_postmuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_kernel{2},mua_norm_exp{2},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

mua_ctxpred_task_k_postmuscimol = arrayfun(@(regressor) ...
    cell2mat(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    task_str_ctxpred_kernel{2},mua_norm_exp{2},'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false)';

% Get task>striatum parameters
n_regressors = length(task_regressor_labels);

% Plot task>striatum kernels
stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = [];
move_col = [0.6,0,0.6;0,0.6,0];
go_col = [0.8,0.8,0.2;0.5,0.5,0.5];
outcome_col = [0,0,0.8;0.8,0,0];
task_regressor_cols = {stim_col,move_col,go_col,outcome_col};
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

figure('Name','Pre-muscimol');
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);
        
        curr_kernels = mua_task_k_premuscimol{curr_regressor}(:,:,curr_depth,:);
        n_subregressors = size(mua_task_k_premuscimol{curr_regressor},1);
        col = task_regressor_cols{curr_regressor};
        for curr_subregressor = 1:n_subregressors
            AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
                nanmean(curr_kernels(curr_subregressor,:,:,:),4)', ...
                AP_sem(curr_kernels(curr_subregressor,:,:,:),4)', ...
                col(curr_subregressor,:));
        end
        
        xlabel('Time (s)');
        ylabel('Weight');
        title(task_regressor_labels{curr_regressor});
        line([0,0],ylim,'color','k');
        
    end
end
linkaxes(p);
y_scale = 1;
t_scale = 0.5;
line(min(xlim) + [0,t_scale],repmat(min(ylim),2,1),'color','k','linewidth',3);
line(repmat(min(xlim),2,1),min(ylim) + [0,y_scale],'color','k','linewidth',3);

figure('Name','Post-muscimol');
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);
        
        curr_kernels = mua_task_k_postmuscimol{curr_regressor}(:,:,curr_depth,:);
        n_subregressors = size(mua_task_k_postmuscimol{curr_regressor},1);
        col = task_regressor_cols{curr_regressor};
        for curr_subregressor = 1:n_subregressors
            AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
                nanmean(curr_kernels(curr_subregressor,:,:,:),4)', ...
                AP_sem(curr_kernels(curr_subregressor,:,:,:),4)', ...
                col(curr_subregressor,:),0.5);
        end
        
        xlabel('Time (s)');
        ylabel('Weight');
        title(task_regressor_labels{curr_regressor});
        line([0,0],ylim,'color','k');
        
    end
end
linkaxes(p);
y_scale = 1;
t_scale = 0.5;
line(min(xlim) + [0,t_scale],repmat(min(ylim),2,1),'color','k','linewidth',3);
line(repmat(min(xlim),2,1),min(ylim) + [0,y_scale],'color','k','linewidth',3);

% Plot kernel sums pre/post muscimol
str_k_sum_premuscimol = cellfun(@(x) permute(sum(x,2),[1,3,4,2]),mua_task_k_premuscimol,'uni',false);
str_k_sum_postmuscimol = cellfun(@(x) permute(sum(x,2),[1,3,4,2]),mua_task_k_postmuscimol,'uni',false);

figure;
p = nan(n_depths,n_regressors);
for curr_depth = 1:n_depths
    for curr_regressor = 1:n_regressors
        if curr_regressor == 1
            x = unique([0.06,0.125,0.25,0.5,1].*[-1;1]);
        else
            x = 1:size(str_k_sum_premuscimol{curr_regressor},1);
        end
        
        p(curr_depth,curr_regressor) = ...
            subplot(n_depths,n_regressors,curr_regressor+(curr_depth-1)*n_regressors);
        hold on
        
        errorbar(x,nanmean(str_k_sum_premuscimol{curr_regressor}(:,curr_depth,:),3), ...
            AP_sem(str_k_sum_premuscimol{curr_regressor}(:,curr_depth,:),3),'k','linewidth',2);
        
        errorbar(x,nanmean(str_k_sum_postmuscimol{curr_regressor}(:,curr_depth,:),3), ...
            AP_sem(str_k_sum_postmuscimol{curr_regressor}(:,curr_depth,:),3),'r','linewidth',2);
        
        axis tight
        xlim(xlim + [-0.2,0.2]);
        
        xlabel('Condition');
        ylabel('Weight sum');
        title(task_regressor_labels{curr_regressor});
        
    end
end

linkaxes(p,'y');















