function im_aligned = AP_align_widefield(animals,days,im_unaligned)
% im_aligned = AP_align_widefield(animals,days,im_unaligned)
%
% Get transform matrix to align widefield images across days
% im_unaligned - unaligned images (if 'new', then create and save tform)
%
% If no days, assume already aligned to within-animal master

%% Check arguments

% Whether to create new alignment or load old
if ischar(im_unaligned) && strcmp(im_unaligned,'new')
    new_alignment = true;
else
    new_alignment = false;
    if ~iscell(im_unaligned)
        im_unaligned = {im_unaligned};
    end
end

% Whether to align within animal across days or across animals
if isempty(days)
    animal_aligned = true;
    if ~iscell(animals)
        animals = {animals};
    end
else
    if ~iscell(days)
        days = {days};
    end
    animal_aligned = false;
end

if ~animal_aligned && iscell(animals) && length(animals) > 1
    error('Only within or across animal alignment possible at a time')
end

% Set alignment save path
alignment_path = '\\basket.cortexlab.net\data\ajpeters\wf_alignment';

%% Apply alignments to data

% Load across animal alignments
load([alignment_path filesep 'animal_wf_tform']);

if ~new_alignment
    
    if ~animal_aligned
        
        % Load within animal alignments
        load([alignment_path filesep animals '_wf_tform']);
        
        curr_animal_idx = strcmp(animals,{animal_wf_tform.animal});
        globally_aligned = any(curr_animal_idx);
        if ~globally_aligned
            warning(['No global alignment: ' animals]);
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
        
        if globally_aligned
            n_px_y = animal_wf_tform(curr_animal_idx).im_size(1);
            n_px_x = animal_wf_tform(curr_animal_idx).im_size(2);
        else
            n_px_y = size(im_unaligned{1},1);
            n_px_x = size(im_unaligned{1},2);
        end
        
        im_aligned = nan(n_px_y,n_px_x,n_frames,n_conditions,length(days));
        
        for curr_day = 1:length(days);
            
            curr_day_idx = strcmp(days{curr_day},{wf_tform.day});
            if ~any(curr_day_idx)
                error(['No transform matrix for ' days{curr_day}]);
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
        
    elseif animal_aligned
        
        % Apply to each day
        n_frames = unique(cellfun(@(x) size(x,3),im_unaligned));
        if length(n_frames) > 1
           error('Different number of frames across animals'); 
        end
        n_conditions = unique(cellfun(@(x) size(x,4),im_unaligned));
        if length(n_conditions) > 1
           error('Different number of frames across animals'); 
        end
        
        im_aligned = nan(animal_wf_tform(1).im_size(1),animal_wf_tform(1).im_size(2),n_frames,n_conditions,length(animals));
        for curr_animal = 1:length(animals);
            
            curr_animal_idx = strcmp(animals{curr_animal},{animal_wf_tform.animal});
            if ~any(curr_animal_idx)
                error(['No transform matrix for ' animals{curr_animal}]);
            end
            
            tform = affine2d;
            tform.T = animal_wf_tform(curr_animal_idx).t;
            
            % (this might be bad for the future, keep it in mind)
            curr_im = im_unaligned{curr_animal};
            curr_im(isnan(curr_im)) = 0;
            
            im_aligned(:,:,:,:,curr_animal) = imwarp(curr_im,tform, ...
                'Outputview',imref2d(animal_wf_tform(curr_animal_idx).im_size));
        end
        % Squeeze any extra dimensions out
        im_aligned = squeeze(im_aligned);
        
    end
    
end

%% If new tform, register average images and save tform
% (this should only be run once per animal)

if new_alignment
    disp('Finding alignment transform...')
    
    avg_im_blue = cell(length(days),1);
    for curr_day = 1:length(days)
        [img_path,img_exists] = AP_cortexlab_filename(animals,days{curr_day},[],'imaging');
        avg_im_blue{curr_day} = readNPY([img_path filesep 'meanImage_blue.npy']);
    end
    
    avg_im_purple = cell(length(days),1);
    for curr_day = 1:length(days)
        [img_path,img_exists] = AP_cortexlab_filename(animals,days{curr_day},[],'imaging');
        avg_im_purple{curr_day} = readNPY([img_path filesep 'meanImage_purple.npy']);
    end
    
    % Choose reference image
    ref_im_full = avg_im_purple{1};
    ref_size = size(ref_im_full);
    
    border_pixels = 50;
    
    % Remove border for aligning
    im_align = cellfun(@(x) x(border_pixels:end-border_pixels+1,border_pixels:end-border_pixels+1),avg_im_purple,'uni',false);
    ref_im = ref_im_full(border_pixels:end-border_pixels+1,border_pixels:end-border_pixels+1);
    
    fprintf('Registering average images:')
    fprintf('\n');
    
    tform_matrix = cell(length(avg_im_blue),1);
    im_aligned = nan(ref_size(1),ref_size(2),length(days));
    for curr_day = 1:length(days);
        
        [optimizer, metric] = imregconfig('monomodal');
        optimizer = registration.optimizer.OnePlusOneEvolutionary();
        optimizer.MaximumIterations = 200;
        optimizer.GrowthFactor = 1+1e-6;
        optimizer.InitialRadius = 1e-4;
        
        tformEstimate_affine = imregtform(im_align{curr_day},ref_im,'affine',optimizer,metric);
        curr_im_reg = imwarp(avg_im_blue{curr_day},tformEstimate_affine,'Outputview',imref2d(ref_size));
        tform_matrix{curr_day} = tformEstimate_affine.T;
        
        im_aligned(:,:,curr_day) = curr_im_reg;
        
        AP_print_progress_fraction(curr_day,length(days))
        
    end
    
    % Save alignments
    wf_tform = struct('day',reshape(days,[],1),'t',tform_matrix,'im_size',size(ref_im_full));
    
    alignment_filename = [alignment_path filesep animals '_wf_tform'];
    save(alignment_filename,'wf_tform');
    
end




%% Align animals by retiotopy 
% (just storing this here for now, already done)

% % Load all retinotopic maps in structure 'retinotopy'
% 
% ref_im = retinotopy{3};
% ref_size = size(ref_im);
% 
% im_align = retinotopy;
% 
% [optimizer, metric] = imregconfig('monomodal');
% optimizer = registration.optimizer.OnePlusOneEvolutionary();
% optimizer.MaximumIterations = 200;
% optimizer.GrowthFactor = 1+1e-6;
% optimizer.InitialRadius = 1e-4;
% 
% im_aligned = nan(ref_size(1),ref_size(2),6);
% for curr_animal = 1:6
%     tformEstimate_affine = imregtform(im_align{curr_animal},ref_im,'affine',optimizer,metric);
%     curr_im_reg = imwarp(retinotopy{curr_animal},tformEstimate_affine,'Outputview',imref2d(ref_size));
%     tform_matrix{curr_animal} = tformEstimate_affine.T;
%     
%     im_aligned(:,:,curr_animal) = curr_im_reg;
% end
% 
% AP_image_scroll(im_aligned);axis image
% 
% animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
% wf_tform = struct('animal',animals','t',tform_matrix','im_size',ref_size);    
% 
% alignment_filename = ['\\basket.cortexlab.net\data\ajpeters\wf_alignment\mice_wf_tform.mat'];
% % save(alignment_filename,'wf_tform');












