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
% 0) Create master alignment reference ('create_master') from average aligned VFSs from multiple animals
% 1) Across-day align ('new_days') using vasulature from all days of one animal
% 2) Across-animal align ('new_animal) using the aligned average VFS from new animal
% 3) Apply saved alignments to any images (specifying animal and day source)


%% Initialize

% Set path with saved alignments
alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment';

% Load transform structure
wf_tform_fn = [alignment_path filesep 'wf_tform.mat'];
% (if it doesn't exist yet, create it)
if ~exist(wf_tform_fn,'file')
    wf_tform = struct('animal',cell(0),'day',cell(0),'day_tform',cell(0),'animal_tform',cell(0));
    save(wf_tform_fn,'wf_tform');
else
    load(wf_tform_fn);
end

% Set empty align_type if unassigned
if ~exist('align_type','var')
    align_type = [];
end

%% Apply or create alignment

switch align_type
    
    case 'new_days'
        %% Align all days from one animal
        % (only run once for each animal - aligns to iterative average)
        
        % Set output size as the largest image
        [im_y,im_x] = cellfun(@size,im_unaligned);
        
        im_y_max = max(im_y);
        im_x_max = max(im_x);
        ref_size = [im_y_max,im_x_max];
        
        im_unaligned_pad = ...
            cell2mat(reshape(cellfun(@(im) padarray(im, ...
            [im_y_max - size(im,1),im_y_max - size(im,2)],0,'post'), ...
            im_unaligned,'uni',false),1,1,[]));
        
        % (set transform optimizer)
        [optimizer, metric] = imregconfig('monomodal');
        optimizer = registration.optimizer.OnePlusOneEvolutionary();
        optimizer.MaximumIterations = 200;
        optimizer.GrowthFactor = 1+1e-6;
        optimizer.InitialRadius = 1e-4;
        
        % (first pass: rigid transform to day 1)             
        im_ref = im_unaligned{1};
        im_rigid_aligned = nan(ref_size(1),ref_size(2),length(im_unaligned));       
        rigid_tform = cell(size(im_unaligned));        
        for curr_im = 1:length(im_unaligned)
            tformEstimate_affine = imregtform(im_unaligned{curr_im},im_ref,'rigid',optimizer,metric);
            curr_im_reg = imwarp(im_unaligned{curr_im},tformEstimate_affine,'Outputview',imref2d(ref_size));
            rigid_tform{curr_im} = tformEstimate_affine.T;
            im_rigid_aligned(:,:,curr_im) = curr_im_reg;
        end
        
        AP_image_scroll(im_rigid_aligned); axis image    
        
        % (user draw ROI around most stable area for alignment reference)               
        f = AP_image_scroll(im_rigid_aligned, ...
            repmat({'Draw ROI over largest stable area'},length(im_unaligned),1));
        axis image;
        h = imrect;
        position = round(wait(h));
        close(f);

        im_rigid_aligned_roi = im_rigid_aligned( ...
            position(2):position(2)+position(4), ...
            position(1):position(1)+position(3),:);
        
        im_ref = im_rigid_aligned_roi(:,:,1);
        im_aligned = nan(ref_size(1),ref_size(2),length(im_unaligned));
        tform_matrix = cell(length(im_unaligned),1);
        for curr_im = 1:length(im_unaligned)
            tformEstimate_affine = imregtform(im_rigid_aligned_roi(:,:,curr_im),im_ref,'affine',optimizer,metric);
            
            tform_combine = rigid_tform{curr_im}*tformEstimate_affine.T;
            tform_matrix{curr_im} = tform_combine;
            
            curr_tform = affine2d;
            curr_tform.T = tform_combine;
            curr_im_reg = imwarp(im_unaligned{curr_im},curr_tform,'Outputview',imref2d(ref_size));
            
            tform_matrix{curr_im} = tformEstimate_affine.T;
            im_aligned(:,:,curr_im) = curr_im_reg;
        end
        
        f = AP_image_scroll(im_aligned, ...
            cellfun(@(x) ['Final alignment: ' animal ' ' x],day,'uni',false));
        axis image;
        confirm_save = input(['Save day transforms for ' animal '? (y/n): '],'s');
        close(f)
        
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
        
        im_aligned = nan(n_px_y,n_px_x,n_frames,n_conditions,length(day));
        
        for curr_day = 1:length(day)
            
            curr_day_idx = strcmp(day{curr_day},{wf_tform.day});
            if ~any(curr_day_idx)
                error(['No transform matrix for ' animal ' ' day{curr_day}]);
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
        % Squeeze any extra dimensions out, cast to original data type
        im_aligned = cast(squeeze(im_aligned),class(im_unaligned{curr_day}));
        
end











