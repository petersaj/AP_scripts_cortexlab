%% Notes
%
%
% Batch scripts to save preprocessed data here, saved to:
% C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning


%% ~~~~~~~~~~~~~ Facecam alignment ~~~~~~~~~~~~~

%% Collect sample facecam frame from all experiments

animals = {'AP100','AP101','AP103','AP104','AP105','AP106','AP107','AP108','AP109','AP111','AP112'};

% Initialize save variable
facecam_sampleframe_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);
    
    % Get days with muscimol
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    if any(muscimol_animal_idx)
        muscimol_start_day = muscimol(muscimol_animal_idx).day{1};
        muscimol_experiments = datenum({experiments.day})' >= datenum(muscimol_start_day);
    else
        muscimol_experiments = false(size({experiments.day}));
    end
    
    % Set experiments to use (imaging, not muscimol)
    experiments = experiments([experiments.imaging] & ~muscimol_experiments);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment (only cam)
        load_parts.cam = true;
        AP_load_experiment       

        if ~facecam_exists
            continue
        end

        % Grab sample frame (last frame)
        vr = VideoReader(facecam_fn);
        grab_frame = vr.NumFrames;
        facecam_sample_frame = read(vr,grab_frame);
        
        facecam_sampleframe_all(curr_animal).animal = animal;
        facecam_sampleframe_all(curr_animal).day{curr_day} = day;
        facecam_sampleframe_all(curr_animal).im{curr_day} = facecam_sample_frame;

        % Prep for next loop
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars('-except',preload_vars{:});

    end
end
disp('Done loading all');

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\facecam_processing';
save_fn = ['facecam_sampleframe_all'];
save([save_path filesep save_fn],'facecam_sampleframe_all','-v7.3');
disp(['Saved: ' save_path filesep save_fn])


%% Align facecam frames (auto)
% Align facecam sample frames across mice/days
% (copied from AP_align_widefield)

% Load facecam sample frames
facecam_processing_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\facecam_processing';
facecam_sampleframe_all_fn = fullfile(facecam_processing_path,'facecam_sampleframe_all.mat');
load(facecam_sampleframe_all_fn);

im_unaligned = cellfun(@double,[facecam_sampleframe_all.im],'uni',false);

% Pad the images for checking
% (set output size as the largest image)
[im_y,im_x] = cellfun(@size,im_unaligned);

im_y_max = max(im_y);
im_x_max = max(im_x);
ref_size = [im_y_max,im_x_max];

im_unaligned_pad = ...
    cell2mat(reshape(cellfun(@(im) padarray(im, ...
    [im_y_max - size(im,1),im_x_max - size(im,2)],0,'post'), ...
    im_unaligned,'uni',false),1,1,[]));

AP_image_scroll(im_unaligned_pad);
axis image;

% (for RegularStepGradientDescent)
[optimizer, metric] = imregconfig('monomodal');
optimizer.GradientMagnitudeTolerance = 1e-6;
optimizer.MaximumIterations = 100;
optimizer.MaximumStepLength = 1e-4;
optimizer.MinimumStepLength = 1e-6;
optimizer.RelaxationFactor = 0.6;

im_unaligned_edge = cellfun(@(x) x-imgaussfilt(x,10),im_unaligned,'uni',false);

% (rigid align to first image)
disp('Rigid aligning images...')
% (use only top half of image - excludes hands which are random)
im_ref = im_unaligned_edge{1}(1:round(size(im_unaligned{1},1)/2),:);

im_rigid_aligned = nan(ref_size(1),ref_size(2),length(im_unaligned));
rigid_tform = cell(size(im_unaligned));
for curr_im = 1:length(im_unaligned)
    if isempty(im_unaligned{curr_im})
        continue
    end
    
    im_source = im_unaligned_edge{curr_im}(1:round(size(im_unaligned{curr_im},1)/2),:);

    tformEstimate_affine = imregtform(im_source, ...
        im_ref,'rigid',optimizer,metric,'PyramidLevels',4);
    curr_im_reg = imwarp(im_unaligned{curr_im}, ...
        tformEstimate_affine,'Outputview',imref2d(ref_size));
    rigid_tform{curr_im} = tformEstimate_affine.T;
    im_rigid_aligned(:,:,curr_im) = curr_im_reg;
    
    if exist('h','var');delete(h);end
    h = figure;
    imshowpair(im_ref,im_rigid_aligned(:,:,curr_im));
    drawnow

    AP_print_progress_fraction(curr_im,length(im_unaligned));
end

AP_image_scroll(im_rigid_aligned)
axis image;

% Save
facecam_align_refsize = ref_size;
facecam_align_tform = rigid_tform;

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\facecam_processing';
save_fn = ['facecam_align_tform'];
save([save_path filesep save_fn],'facecam_align_refsize','facecam_align_tform','-v7.3');
disp(['Saved: ' save_path filesep save_fn])

%% Align facecam frames (control point)

% Load facecam sample frames
facecam_processing_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\facecam_processing';
facecam_sampleframe_all_fn = fullfile(facecam_processing_path,'facecam_sampleframe_all.mat');
load(facecam_sampleframe_all_fn);

% Plot images and select control points
im_unaligned = cellfun(@double,[facecam_sampleframe_all.im],'uni',false);

target_im = im_unaligned{1};
source_im = im_unaligned{20};

figure;
target_ax = subplot(1,2,1);
imagesc(target_im);
axis image off; hold on;
source_ax = subplot(1,2,2);
imagesc(source_im);
axis image off; hold on;

target_ctrl_points = ginput;
plot(target_ax,target_ctrl_points(:,1),target_ctrl_points(:,2),'.r','MarkerSize',20);

source_ctrl_points = ginput;
plot(source_ax,source_ctrl_points(:,1),source_ctrl_points(:,2),'.b','MarkerSize',20);

% Align images
tform = fitgeotrans(source_ctrl_points,target_ctrl_points,'nonreflectivesimilarity');
tform_size = imref2d(size(target_im));
source_im_aligned = imwarp(source_im,tform,'OutputView',tform_size);

target_im_edge = target_im - imgaussfilt(target_im,10);









