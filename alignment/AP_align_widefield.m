function im_aligned = AP_align_widefield(im_unaligned,animal,day,align_type)
% im_aligned = AP_align_widefield(im_unaligned,animal,day,align_type)
%
% Align widefield images across days and animals
%
% im_unaligned - unaligned images
% animal - animal corresponding to im_unaligned, e.g. 'XX001'
% day - 'yyyy-mm-dd' corresponding to im_unaligned (cell array if multiple)
% align_type - type of alignment:
% > default (otherwise/empty) - apply saved day/animal transform to image(s)
% > 'new_days' - make new alignment within animal across days
% > 'new_animal' - make new alignment from animal to master 
% > 'create_master' - create new master alignment reference
%
% USAGE: 
% 1) Create master alignment reference ('create_master') from average aligned VFSs from multiple animals
% 2) Across-day align ('new_days') using vasulature from all days of one animal
% 3) Across-animal align ('new_animal) using the aligned average VFS from new animal
% 4) Apply saved alignments to any images (specifying animal and day source)


%% Initialize

% Set path with saved alignments
alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment';

% Load transform structure
wf_tform_fn = [alignment_path filesep 'wf_tform.mat'];
% (if it doesn't exist yet, create it)
if ~exist(wf_tform_fn,'file')
    wf_tform = struct('animal',cell(0),'day',cell(0),'day_tform',cell(0),'ref_size',cell(0),'animal_tform',cell(0),'master_ref_size',cell(0));
    save(wf_tform_fn,'wf_tform');
else
    load(wf_tform_fn);
end

% Set empty align_type if unassigned
if ~exist('align_type','var')
    align_type = '';
end

%% Apply or create alignment

switch align_type
    
    case 'new_days'
        %% Align all days from one animal
        % (only run once for each animal - aligns to iterative average)
        
        % Intensity-normalize the unaligned images
        % (if there are big differences, this greatly improves alignment)
        im_unaligned_norm = cellfun(@(x) x./mad(x(:),true),im_unaligned,'uni',false);
        
        % Set output size as the largest image
        [im_y,im_x] = cellfun(@size,im_unaligned_norm);
        
        im_y_max = max(im_y);
        im_x_max = max(im_x);
        ref_size = [im_y_max,im_x_max];
        
        im_unaligned_pad = ...
            cell2mat(reshape(cellfun(@(im) padarray(im, ...
            [im_y_max - size(im,1),im_x_max - size(im,2)],0,'post'), ...
            im_unaligned_norm,'uni',false),1,1,[]));
        
        % (set transform optimizer)
        [optimizer, metric] = imregconfig('monomodal');
        optimizer = registration.optimizer.OnePlusOneEvolutionary();
        optimizer.MaximumIterations = 200;
        optimizer.GrowthFactor = 1+1e-6;
        optimizer.InitialRadius = 1e-4;      
        
        % (first pass: rigid transform to day 1)  
        disp('Rigid aligning images...')
        im_ref = im_unaligned_norm{1};
        im_rigid_aligned = nan(ref_size(1),ref_size(2),length(im_unaligned_norm));       
        rigid_tform = cell(size(im_unaligned_norm));        
        for curr_im = 1:length(im_unaligned_norm)
            tformEstimate_affine = imregtform(im_unaligned_norm{curr_im},im_ref,'rigid',optimizer,metric);
            curr_im_reg = imwarp(im_unaligned_norm{curr_im},tformEstimate_affine,'Outputview',imref2d(ref_size));
            rigid_tform{curr_im} = tformEstimate_affine.T;
            im_rigid_aligned(:,:,curr_im) = curr_im_reg;
        end
        
        % (user draw ROI around most stable area for alignment reference)
        draw_roi = true;
        while draw_roi
            f = AP_image_scroll(im_rigid_aligned, ...
                repmat({'Draw ROI over largest stable area'},length(im_unaligned_norm),1));
            caxis([0,prctile(im_unaligned_pad(:),99)]);
            axis image;
            h = imrect;
            position = round(wait(h));
            close(f);
            
            im_rigid_aligned_roi = im_rigid_aligned( ...
                position(2):position(2)+position(4), ...
                position(1):position(1)+position(3),:);
                                    
            disp('Rigid aligning selection...')
            im_ref = im_rigid_aligned_roi(:,:,1);
            im_aligned = nan(ref_size(1),ref_size(2),length(im_unaligned_norm));
            tform_matrix = cell(length(im_unaligned_norm),1);
            for curr_im = 1:length(im_unaligned_norm)
                tformEstimate_affine = imregtform(im_rigid_aligned_roi(:,:,curr_im),im_ref,'rigid',optimizer,metric);
                
                tform_combine = rigid_tform{curr_im}*tformEstimate_affine.T;
                tform_matrix{curr_im} = tform_combine;
                
                curr_tform = affine2d;
                curr_tform.T = tform_combine;
                curr_im_reg = imwarp(im_unaligned_norm{curr_im},curr_tform,'Outputview',imref2d(ref_size));
                
                tform_matrix{curr_im} = tform_combine;
                im_aligned(:,:,curr_im) = curr_im_reg;
            end
            
            f = AP_image_scroll(im_aligned, ...
                cellfun(@(x) ['Final alignment: ' animal ' ' x],day,'uni',false));
            caxis([0,prctile(im_aligned(:),99)]);
            axis image;
            confirm_save = input(['Save day transforms for ' animal '? (y/n/redo): '],'s');
            close(f)
            
            if ~strcmp(confirm_save,'redo')
                draw_roi = false;
            end
        end
        
        % Save transform matrix into structure
        if strcmp(confirm_save,'y')
            curr_animal_idx = strcmp(animal,{wf_tform.animal});
            if isempty(curr_animal_idx) || ~any(curr_animal_idx)
                % (if not already extant, make new)
                curr_animal_idx = length(wf_tform) + 1;
                wf_tform(curr_animal_idx).animal = animal;
            end
            
            wf_tform(curr_animal_idx).day = reshape(day,[],1);
            wf_tform(curr_animal_idx).day_tform = tform_matrix;
            wf_tform(curr_animal_idx).ref_size = ref_size;
            
            save(wf_tform_fn,'wf_tform');
            disp(['Saved day transforms for ' animal '.'])
        else
            disp('Not saved');
        end
        
        
    case 'create_master'
        %% Create master VFS
        % NOTE: changing this requires re-aligning animals, ideally only done once
        
        disp('Creating new master VFS...')
        
        % Check inputs
        if ~iscell(im_unaligned)
            error('Expected cell array of unaligned VFS''s');
        end
        
        % Pad and concatenate unaligned VFS's
        [vfs_y,vfs_x] = cellfun(@size,im_unaligned);
        
        vfs_y_max = max(vfs_y);
        vfs_x_max = max(vfs_x);
        
        vfs_unaligned_pad = ...
            cell2mat(reshape(cellfun(@(vfs) padarray(vfs, ...
            [vfs_y_max - size(vfs,1),vfs_y_max - size(vfs,2)],0,'post'), ...
            im_unaligned,'uni',false),1,1,[]));
        
        % Iterate averaging the VFS and the aligning
        disp('Iterating aligning average VFS...')
        
        vfs_ref = nanmean(vfs_unaligned_pad,3);
        ref_size = size(vfs_ref);
        
        [optimizer, metric] = imregconfig('monomodal');
        optimizer = registration.optimizer.OnePlusOneEvolutionary();
        optimizer.MaximumIterations = 200;
        optimizer.GrowthFactor = 1+1e-6;
        optimizer.InitialRadius = 1e-4;
        
        n_loops = 5;
        for curr_loop = 1:n_loops
            
            vfs_aligned = nan(ref_size(1),ref_size(2),length(im_unaligned));
            for curr_animal = 1:length(im_unaligned)
                tformEstimate_affine = imregtform(im_unaligned{curr_animal},vfs_ref,'affine',optimizer,metric);
                curr_im_reg = imwarp(im_unaligned{curr_animal},tformEstimate_affine,'Outputview',imref2d(ref_size));
                tform_matrix{curr_animal} = tformEstimate_affine.T;
                vfs_aligned(:,:,curr_animal) = curr_im_reg;
            end
            
            vfs_ref = nanmean(vfs_aligned,3);
            AP_print_progress_fraction(curr_loop,n_loops);
        end
        
        % Symmetrize the average aligned VFS
        disp('Symmetrizing aligned average VFS...')
        
        n_loops = 10;
        ref_size = size(vfs_ref);
        
        ref_reg = nan(size(vfs_ref,1),size(vfs_ref,2),n_loops);
        for curr_loop = 1:n_loops
            ref_im_symm = (vfs_ref + fliplr(vfs_ref))./2;
            
            [optimizer, metric] = imregconfig('monomodal');
            optimizer = registration.optimizer.OnePlusOneEvolutionary();
            optimizer.MaximumIterations = 200;
            optimizer.GrowthFactor = 1+1e-6;
            optimizer.InitialRadius = 1e-4;
            
            tformEstimate_affine = imregtform(vfs_ref,ref_im_symm,'rigid',optimizer,metric);
            curr_im_reg = imwarp(vfs_ref,tformEstimate_affine,'Outputview',imref2d(ref_size));
            tform_matrix = tformEstimate_affine.T;
            vfs_ref = curr_im_reg;
            
            ref_reg(:,:,curr_loop) = vfs_ref;
            AP_print_progress_fraction(curr_loop,n_loops);
        end
        
        % Set make symmetric-average map the master
        symm_vfs = ref_reg(:,:,end);
        symm_vfs_mirror_avg = (symm_vfs + fliplr(symm_vfs))./2;
        
        master_vfs = symm_vfs_mirror_avg;
        
        figure;
        imagesc(master_vfs);
        axis image off
        colormap(brewermap([],'*RdBu'));
        caxis([-1,1]);
        title('Master VFS');
        
        % Save the master
        confirm_save = input(['Overwrite master VFS? (y/n): '],'s');
        if strcmp(confirm_save,'y')
            master_vfs_fn = [alignment_path filesep 'master_vfs.mat'];
            save(master_vfs_fn,'master_vfs');
            disp('Saved new master VFS.');
            confirm_align_ccf = input(['Align CCF to new master VFS? (y/n): '],'s');
            if strcmp(confirm_save,'y')
                AP_vfs_ccf_align;
            end
        else
            disp('Not saved.')
        end
        

    case 'new_animal'
        %% Align animal to master VFS
        
        disp('Aligning animal VFS to master VFS...')
        
        % Load master VFS
        master_vfs_fn = [alignment_path filesep 'master_vfs.mat'];
        load(master_vfs_fn);
        
        % Align animal VFS to master VFS
        ref_size = size(master_vfs);
        
        [optimizer, metric] = imregconfig('monomodal');
        optimizer = registration.optimizer.OnePlusOneEvolutionary();
        optimizer.MaximumIterations = 200;
        optimizer.GrowthFactor = 1+1e-6;
        optimizer.InitialRadius = 1e-4;
        
        tformEstimate_affine = imregtform(im_unaligned,master_vfs,'affine',optimizer,metric);
        im_aligned = imwarp(im_unaligned,tformEstimate_affine,'Outputview',imref2d(ref_size));
        tform_matrix = tformEstimate_affine.T;
        
        figure;
        subplot(1,2,1);
        imshowpair(master_vfs,im_unaligned);
        title('Unaligned');
        subplot(1,2,2);
        imshowpair(master_vfs,im_aligned);
        title('Aligned');
        
        % Save transform matrix into structure
        curr_animal_idx = strcmp(animal,{wf_tform.animal});
        if isempty(curr_animal_idx) || ~any(curr_animal_idx)
            % (if not already extant, make new)
            new_animal_idx = length(wf_tform) + 1;
            wf_tform(new_animal_idx).animal = animal;
            wf_tform(new_animal_idx).animal_tform = tform_matrix;
            save(wf_tform_fn,'wf_tform');
            disp(['Saved new animal alignment for ' animal '.'])
        else
            % (if extant, prompt to overwrite)
            confirm_save = input(['Overwrite animal transform for ' animal '? (y/n): '],'s');
            if strcmp(confirm_save,'y')
                wf_tform(curr_animal_idx).animal_tform = tform_matrix;
                save(wf_tform_fn,'wf_tform');
                disp(['Saved new animal alignment for ' animal '.'])
            else
                disp('Not overwriting.')
            end
        end        
        
        
    otherwise
        %% Apply alignments to data
        
        % Find animal and day index within wf_tform structure
        curr_animal_idx = strcmp(animal,{wf_tform.animal});
        curr_day_idx = strcmp(day,wf_tform(curr_animal_idx).day);
        
        if ~any(curr_animal_idx)
            disp(['No alignments found for ' animal]);
        elseif ~any(curr_day_idx)
            error(['No ' animal ' alignment found for ' day]);
        end
        
        % Get transform for day
        curr_day_tform = wf_tform(curr_animal_idx).day_tform{curr_day_idx};
        
        % Get transform for animal (if it exists) and set image size
        if ~isempty(wf_tform(curr_animal_idx).animal_tform)
            curr_animal_tform = wf_tform(curr_animal_idx).animal_tform;
            curr_tform = curr_day_tform*curr_animal_tform;
            ref_size = wf_tform.master_ref_size;
        else
            warning(['No master animal alignment for ' animal]);
            curr_tform = curr_day_tform;
            ref_size = wf_tform(curr_animal_idx).ref_size;
        end
        
        % Transform image, cast to input data type
        tform = affine2d;
        tform.T = curr_tform;
        im_aligned = cast(imwarp(im_unaligned,tform,'Outputview',imref2d(ref_size)), ...
            class(im_unaligned));
        

end











