function im_aligned = AP_align_widefield(animal,days,im_unaligned)
% im_aligned = AP_align_widefield(animal,days,im_unaligned)
% 
% Get transform matrix to align widefield images across days
% im_unaligned - unaligned images (if 'new', then create and save tform)

%% Check arguments
if iscell(im_unaligned)
    new_alignment = false;
elseif ischar(im_unaligned) && strcmp(im_unaligned,'new')
    new_alignment = true;
else
    error('Alignment option invalid');
end

%% Set saved alignment filename

alignment_path = '\\basket.cortexlab.net\data\ajpeters\wf_alignment';
alignment_filename = [alignment_path filesep animal '_wf_tform'];

%% If new tform, register average images and save tform

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
    save(alignment_filename,'wf_tform');
    
end

%% If not, apply old alignment to data
if ~new_alignment
    
    % Load alignment tform
    load(alignment_filename,'wf_tform');
    
    % Apply to each day
    im_aligned = nan(wf_tform(1).im_size(1),wf_tform(1).im_size(2),length(days));
    for curr_day = 1:length(days);
        
        curr_day_idx = strcmp(days{curr_day},{wf_tform.day});
        if ~any(curr_day_idx)
            error(['No transform matrix for ' days{curr_day}]);
        end
        
        tform = affine2d;
        tform.T = wf_tform(curr_day_idx).t;
        
        % (this might be bad for the future, keep it in mind)
        curr_im = im_unaligned{curr_day};
        curr_im(isnan(curr_im)) = 0;
        
        im_aligned(:,:,curr_day) = imwarp(curr_im,tform, ...
            'Outputview',imref2d(wf_tform(curr_day_idx).im_size));
        
    end

end





















