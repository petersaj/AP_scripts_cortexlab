function im_aligned = AP_align_widefield(animal,days,im_unaligned)
% im_aligned = AP_align_widefield(animals,days,im_unaligned)
%
% Get transform matrix to align widefield images across days, combine
% within- and across- animal transforms if available
%
% animals - only more than 1 if across-animal (rare, run manually)
% days - 'yyyy-mm-dd',multiple if unaligned images across days
% im_unaligned - unaligned images (if 'new', then create and save tform)
%
% NOTE: some manual one-time-only code lives at the bottom 

%% Check arguments

% Set flag to create new alignment or load old
if ischar(im_unaligned) && strcmp(im_unaligned,'new')
    new_alignment = true;
else
    new_alignment = false;
    if ~iscell(im_unaligned)
        im_unaligned = {im_unaligned};
    end
end

% Standardize variable types and do checks
if ~new_alignment
    if isempty(days)
        if ~iscell(animal)
            animal = {animal};
        end
    else
        if ~iscell(days)
            days = {days};
        end
        if length(days) ~= length(im_unaligned)
            error('Different number of days and images')
        end
    end
end

% Set alignment save path
alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment';

%% Apply alignments to data

% Load across animal alignments
load([alignment_path filesep 'animal_wf_tform']);

if ~new_alignment
    
    % Load within-animal alignments
    load([alignment_path filesep animal '_wf_tform']);
    
    % Load across-animal alignments
    curr_animal_idx = strcmp(animal,{animal_wf_tform.animal});
    globally_aligned = any(curr_animal_idx);
    if ~globally_aligned
        warning(['No global alignment: ' animal]);
    end
    
    % Apply to each day
    n_frames = unique(cellfun(@(x) size(x,3),im_unaligned));
    if length(n_frames) > 1
        error('Different number of frames across animals');
    end
    n_conditions = unique(cellfun(@(x) size(x,4),im_unaligned));
    if length(n_conditions) > 1
        error('Different number of frames across animals');
    end
    
    % Get target size (master if across animals, first if within)
    if globally_aligned
        n_px_y = animal_wf_tform(curr_animal_idx).im_size(1);
        n_px_x = animal_wf_tform(curr_animal_idx).im_size(2);
    else
        n_px_y = size(im_unaligned{1},1);
        n_px_x = size(im_unaligned{1},2);
    end
    
    im_aligned = nan(n_px_y,n_px_x,n_frames,n_conditions,length(days));
    
    for curr_day = 1:length(days)
        
        curr_day_idx = strcmp(days{curr_day},{wf_tform.day});
        if ~any(curr_day_idx)
            error(['No transform matrix for ' animal ' ' days{curr_day}]);
        end
        
        tform = affine2d;
        if globally_aligned
            tform.T = wf_tform(curr_day_idx).t*animal_wf_tform(curr_animal_idx).t;
        else
            tform.T = wf_tform(curr_day_idx).t;
        end
        
        % (this might be bad for the future, keep it in mind)
        curr_im = im_unaligned{curr_day};
        curr_im(isnan(curr_im)) = 0;
        
        im_aligned(:,:,:,:,curr_day) = imwarp(curr_im,tform, ...
            'Outputview',imref2d([n_px_y,n_px_x]));
    end
    % Squeeze any extra dimensions out
    im_aligned = squeeze(im_aligned);
    
end

%% If new tform, register average images and save tform
% (this should only be run once per animal)

if new_alignment
    disp('Finding alignment transform...')
    
    avg_im_blue = cell(length(days),1);
    for curr_day = 1:length(days)
        [img_path,img_exists] = AP_cortexlab_filename(animal,days{curr_day},[],'imaging');
        avg_im_blue{curr_day} = readNPY([img_path filesep 'meanImage_blue.npy']);
    end
    
    avg_im_purple = cell(length(days),1);
    for curr_day = 1:length(days)
        [img_path,img_exists] = AP_cortexlab_filename(animal,days{curr_day},[],'imaging');
        avg_im_purple{curr_day} = readNPY([img_path filesep 'meanImage_purple.npy']);
    end
    
    % Choose reference image
    % (currently: align to first day, align with purple)
    ref_im_full = avg_im_purple{1};
    ref_size = size(ref_im_full);
    
    border_pixels = 30; % default: 50
    
    % Remove border for aligning
    im_align = cellfun(@(x) x(border_pixels:end-border_pixels+1,border_pixels:end-border_pixels+1),avg_im_purple,'uni',false);
    ref_im = ref_im_full(border_pixels:end-border_pixels+1,border_pixels:end-border_pixels+1);
    
    fprintf('Registering average images:')
    fprintf('\n');
    
    tform_matrix = cell(length(avg_im_blue),1);
    im_aligned = nan(ref_size(1),ref_size(2),length(days));
    for curr_day = 1:length(days)
        
        [optimizer, metric] = imregconfig('monomodal');
        optimizer = registration.optimizer.OnePlusOneEvolutionary();
        optimizer.MaximumIterations = 200; % default: 200
        optimizer.GrowthFactor = 1+1e-6; % default: 1+1e-6
        optimizer.InitialRadius = 1e-4; % default: 1e-4
        
        tformEstimate_affine = imregtform(im_align{curr_day},ref_im,'affine',optimizer,metric);
        curr_im_reg = imwarp(avg_im_blue{curr_day},tformEstimate_affine,'Outputview',imref2d(ref_size));
        tform_matrix{curr_day} = tformEstimate_affine.T;
        
        im_aligned(:,:,curr_day) = curr_im_reg;
        
        AP_print_progress_fraction(curr_day,length(days))
        
    end
    
    % Show aligned images
    AP_image_scroll(im_aligned,days);
    axis image off;
    colormap(gray);
    caxis([0,prctile(im_aligned(:),95)]);
    set(gcf,'Name',['Aligned images: ' animal]);
    
    % Save alignments
    wf_tform = struct('day',reshape(days,[],1),'t',tform_matrix,'im_size',size(ref_im_full));
    
    alignment_filename = [alignment_path filesep animal '_wf_tform'];
    save(alignment_filename,'wf_tform');
    
end


%% STORING HERE/RUN ONCE: Align animals by retinotopy, save U_master_animal
% (TO DO: INTEGRATE INTO ABOVE FUNCTION 'ANIMALS' OPTION)
% (this creates a new combined retinotopy - run AP_vfs_ccf_align afterwards)

% Secure this to be manual
if false
    
    % Set reference animal (never change)
    ref_animal = 'AP026';
    
    % Get retinotopy from all animals
    % (these are kept as figures because I never updated to save as mats)
    retinotopy_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\retinotopy';
    retinotopy_dir = dir(retinotopy_path);
    
    animal_retinotopy_idx = cellfun(@(x) ~isempty(x), regexp({retinotopy_dir.name},'AP\d*_retinotopy'));
    animals_tokens = cellfun(@(x) regexp({x},'(AP\d*)_retinotopy','tokens'),{retinotopy_dir.name});
    animals = cellfun(@(x) cell2mat(x{:}),animals_tokens(cellfun(@(x) ~isempty(x),animals_tokens)),'uni',false);
    
    retinotopy_unaligned = cell(size(animals));
    for curr_animal = 1:length(animals)
        h = open([retinotopy_path filesep animals{curr_animal} '_retinotopy.fig']);
        retinotopy_unaligned{curr_animal} = get(get(subplot(1,2,1),'Children'),'CData');
        close(h);
    end
    
    ref_animal_idx = strcmp(ref_animal,animals);
    ref_im = retinotopy_unaligned{ref_animal_idx};
    ref_size = size(ref_im);
    
    [optimizer, metric] = imregconfig('monomodal');
    optimizer = registration.optimizer.OnePlusOneEvolutionary();
    optimizer.MaximumIterations = 200;
    optimizer.GrowthFactor = 1+1e-6;
    optimizer.InitialRadius = 1e-4;
    
    disp('Aligning retinotopy across animals...');
    im_aligned = nan(ref_size(1),ref_size(2),length(animals));
    for curr_animal = 1:length(animals)
        tformEstimate_affine = imregtform(retinotopy_unaligned{curr_animal},ref_im,'affine',optimizer,metric);
        curr_im_reg = imwarp(retinotopy_unaligned{curr_animal},tformEstimate_affine,'Outputview',imref2d(ref_size));
        tform_matrix{curr_animal} = tformEstimate_affine.T;
        
        im_aligned(:,:,curr_animal) = curr_im_reg;
        AP_print_progress_fraction(curr_animal,length(animals));
    end
    
    AP_image_scroll(im_aligned,animals);axis image
    
    h = figure;imagesc(nanmean(im_aligned,3));
    axis image off;
    colormap(brewermap([],'*RdBu'));
    caxis([-0.5,0.5]);
    title('Combined retinotopy');
    savefig(h,'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\retinotopy\combined_retinotopy')
    
    animal_wf_tform = struct('animal',animals','t',tform_matrix','im_size',ref_size);
    
    alignment_filename = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\animal_wf_tform.mat'];
    save(alignment_filename,'animal_wf_tform');
    
end

%% STORING HERE/RUN ONCE: Save U_master for each animal
% (just keeping this here, independent from above function, run manually)
% (using last day before ephys as master day since probably best
% performance and least bugs)

% Secure this to be manual
if false
    
    warning('Alignment produces non-orthonormal Us!!')
    
%     animals = {'AP024','AP025','AP026','AP027','AP028','AP029', ...
%         'AP032','AP033','AP034','AP035','AP036'};

    % (only trained animals at the moment: naive doesn't have wf only days)
    animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

    for curr_animal = 1:length(animals)
        
        clearvars -except animals curr_animal
        
        animal = animals{curr_animal};
        % Find imaging experiments without ephys (unobstructed cortex)
        % (use only behavior days because cortical recordings afterwards)
        protocol = 'vanillaChoiceworld';
        experiments = AP_find_experiments(animal,protocol);
        experiments(~[experiments.imaging] | [experiments.ephys]) = [];
        if isempty(experiments)
            % (if no behavior days then it was a naive mouse - use passive expt)
            protocol = 'AP_choiceWorldStimPassive';
            experiments = AP_find_experiments(animal,protocol);
            experiments(~[experiments.imaging] | [experiments.ephys]) = [];
        end
        
        load_parts.cam = false;
        load_parts.imaging = true;
        load_parts.ephys = false;
        
        % Use the last day before ephys as the master day
        day = experiments(end).day;
        experiment = experiments(end).experiment(end);
        
        % Align U's of the day to master alignment, save as animal master U
        AP_load_experiment
        U_master_animal = single(AP_align_widefield(animal,day,Udf));
        
        wf_alignment_dir = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment';
        U_master_animal_filename = [wf_alignment_dir filesep 'U_master_' animal];
        save(U_master_animal_filename,'U_master_animal');
        
        AP_print_progress_fraction(curr_animal,length(animals));
        
    end
    
end


%% STORING HERE/RUN ONCE: Create U_master from SVD across animals

% Secure this to be manual
if false
    
    % Load the animal master Us
    animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
    U_master_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment';
    
    U_master_animals = cell(size(animals));
    for curr_animal = 1:length(animals)
        load([U_master_path filesep 'U_master_' animals{curr_animal}]);
        U_master_animals{curr_animal} = U_master_animal;
        clear U_master_animal
        AP_print_progress_fraction(curr_animal,length(animals))
    end    
    
    % Compute SVD across animal master Us 
    % (2000 components is over the GPU limit, slow but only run once)
    n_vs = 2000;
    U_master_animal_cat = cell2mat(cellfun(@(x) ...
            reshape(x(:,:,1:n_vs),[],n_vs),U_master_animals,'uni',false));
        
    [U_U,~,~] = svd(U_master_animal_cat,'econ');
      
    U_master = reshape(U_U(:,1:n_vs),size(U_master_animals{1},1),size(U_master_animals{1},2),n_vs);
    save([U_master_path filesep 'U_master.mat'],'U_master');    
    
end












