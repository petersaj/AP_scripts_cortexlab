function im_aligned = AP_align_widefield(im_unaligned,animal,day,align_type,master_align)
% im_aligned = AP_align_widefield(im_unaligned,animal,day,align_type,master_align)
%
% Align widefield images across days and animals
%
% im_unaligned - unaligned images
% animal - animal corresponding to im_unaligned, e.g. 'XX001'
% day - 'yyyy-mm-dd' corresponding to im_unaligned (cell array if multiple)
% align_type - type of alignment:
% master_align - used with 'new_animal': master template for alignment
% > default (otherwise/empty) - apply saved day/animal transform to image(s)
% > 'day_only' - apply only day (only) transform to image(s)
% > 'animal_only' - apply only animal transform (e.g. already day-aligned)
% > 'new_days' - make new alignment within animal across days
% > 'new_animal' - make new alignment from animal to master 
% > 'create_master' - create new master alignment reference (ONCE!!)
% > 'create_submaster' - create sub-master (align to master)
% master_align - if 'new_animal', use this instead of master_vfs
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
        
        % Set output size as the largest image
        [im_y,im_x] = cellfun(@size,im_unaligned);
        
        im_y_max = max(im_y);
        im_x_max = max(im_x);
        ref_size = [im_y_max,im_x_max];
        
        im_unaligned_pad = ...
            cell2mat(reshape(cellfun(@(im) padarray(im, ...
            [im_y_max - size(im,1),im_x_max - size(im,2)],0,'post'), ...
            im_unaligned,'uni',false),1,1,[]));
        
        % (set transform optimizer)
%         % (for OnePlusOneEvolutionary)
%         [optimizer, metric] = imregconfig('monomodal');
%         optimizer = registration.optimizer.OnePlusOneEvolutionary();
%         optimizer.MaximumIterations = 200; % 200
%         optimizer.GrowthFactor = 1+1e-6; % 1+1e-6
%         optimizer.InitialRadius = 5e-5; % 1e-5
        
        % (for RegularStepGradientDescent)
        [optimizer, metric] = imregconfig('monomodal');
        optimizer.GradientMagnitudeTolerance = 1e-7;
        optimizer.MaximumIterations = 300;
        optimizer.MaximumStepLength = 1e-3;
        optimizer.MinimumStepLength = 1e-7;
        optimizer.RelaxationFactor = 0.6;
        
        % (first pass: rigid transform to day 1)  
        % (increasing the PyramidLevels made the biggest improvement)
        disp('Rigid aligning images...')
        im_ref = im_unaligned{1};
        im_rigid_aligned = nan(ref_size(1),ref_size(2),length(im_unaligned));       
        rigid_tform = cell(size(im_unaligned));        
        for curr_im = 1:length(im_unaligned)
            tformEstimate_affine = imregtform(im_unaligned{curr_im}, ...
                im_ref,'rigid',optimizer,metric,'PyramidLevels',5);
            curr_im_reg = imwarp(im_unaligned{curr_im}, ...
                tformEstimate_affine,'Outputview',imref2d(ref_size));
            rigid_tform{curr_im} = tformEstimate_affine.T;
            im_rigid_aligned(:,:,curr_im) = curr_im_reg;
        end
        
        % (user draw ROI around most stable area for alignment reference)
        % (and lower the search radius)
        draw_roi = true;
        while draw_roi
            f = AP_imscroll(im_rigid_aligned, ...
                repmat({'Draw ROI or double-click to accept'},length(im_unaligned),1));
            caxis([prctile(im_unaligned_pad(:),1),prctile(im_unaligned_pad(:),99)]);
            axis image;
            h = imrect;
            position = round(wait(h));
            close(f);
            
            % If no ROI drawn, accept and save transform
            if sum(position(3:4)) == 0
                tform_matrix = rigid_tform;
                draw_roi = false; 
                confirm_save = 'y';
                continue
            end
            
            im_rigid_aligned_roi = im_rigid_aligned( ...
                position(2):position(2)+position(4), ...
                position(1):position(1)+position(3),:);
                                    
            disp('Rigid aligning selection...')
            im_ref = im_rigid_aligned_roi(:,:,1);
            im_aligned = nan(ref_size(1),ref_size(2),length(im_unaligned));
            tform_matrix = cell(length(im_unaligned),1);
            for curr_im = 1:length(im_unaligned)
                tformEstimate_affine = imregtform( ...
                    im_rigid_aligned_roi(:,:,curr_im), ...
                    im_ref,'rigid',optimizer,metric,'PyramidLevels',4);
                
                tform_combine = rigid_tform{curr_im}*tformEstimate_affine.T;
                tform_matrix{curr_im} = tform_combine;
                
                curr_tform = affine2d;
                curr_tform.T = tform_combine;
                curr_im_reg = imwarp(im_unaligned{curr_im},curr_tform, ...
                    'Outputview',imref2d(ref_size));
                
                tform_matrix{curr_im} = tform_combine;
                im_aligned(:,:,curr_im) = curr_im_reg;
            end
            
            f = AP_imscroll(im_aligned, ...
                cellfun(@(x) ['Final alignment: ' animal ' ' x],day,'uni',false));
            caxis([prctile(im_aligned(:),1),prctile(im_aligned(:),99)]);
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
        % Whenever this is run, any dependent analysis must be re-run
        error('Really create new master? Manually override');
        
        disp('Creating new master VFS...')
        
        % Check inputs
        if ~iscell(im_unaligned)
            error('Expected cell array of unaligned VFS''s');
        end
        
        % Check with user: are L and R hemi same or opposite chirality?
        LR_symmetry = questdlg('Are left/right hemisphere colors same or opposite?','Choose VFS type','symmetric','opposite','symmetric');
        switch LR_symmetry
            case 'symmetric'
                side_sign = 1;
            case 'opposite'
                side_sign = -1;
        end
        % Pad and concatenate unaligned VFS's
        [vfs_y,vfs_x] = cellfun(@size,im_unaligned);
        
        vfs_y_max = max(vfs_y);
        vfs_x_max = max(vfs_x);
        
        vfs_unaligned_pad = ...
            cell2mat(reshape(cellfun(@(vfs) padarray(vfs, ...
            [vfs_y_max - size(vfs,1),vfs_x_max - size(vfs,2)],0,'post'), ...
            im_unaligned,'uni',false),1,1,[]));
        
        % Set alignment parameters     
        [optimizer, metric] = imregconfig('monomodal');
        optimizer = registration.optimizer.OnePlusOneEvolutionary();
        optimizer.MaximumIterations = 200;
        optimizer.GrowthFactor = 1+1e-6;
        optimizer.InitialRadius = 1e-4;
        
        % Iterate averaging the VFS and the aligning
        disp('Iterating aligning average VFS...')
        
        vfs_ref = nanmean(vfs_unaligned_pad,3);
        ref_size = size(vfs_ref);       
        
        n_loops = 5;
        for curr_loop = 1:n_loops
            
            vfs_aligned = nan(ref_size(1),ref_size(2),length(im_unaligned));
            for curr_animal = 1:length(im_unaligned)
                tformEstimate_affine = imregtform(im_unaligned{curr_animal}, ...
                    vfs_ref,'affine',optimizer,metric,'PyramidLevels',5);
                curr_im_reg = imwarp(im_unaligned{curr_animal}, ...
                    tformEstimate_affine,'Outputview',imref2d(ref_size));
                tform_matrix{curr_animal} = tformEstimate_affine.T;
                vfs_aligned(:,:,curr_animal) = curr_im_reg;
            end
            
            vfs_ref = nanmean(vfs_aligned,3);
            AP_print_progress_fraction(curr_loop,n_loops);
        end
        
        % Symmetrize the average aligned VFS
        disp('Symmetrizing aligned average VFS...')
        
        n_loops = 15;
        ref_size = size(vfs_ref);
        
        ref_reg = nan(size(vfs_ref,1),size(vfs_ref,2),n_loops);
        for curr_loop = 1:n_loops
            
            ref_im_symm = (vfs_ref + side_sign*fliplr(vfs_ref))./2;            
            
            tformEstimate_affine = imregtform(vfs_ref,ref_im_symm,'rigid', ...
                optimizer,metric,'PyramidLevels',5);
            curr_im_reg = imwarp(vfs_ref,tformEstimate_affine,'Outputview',imref2d(ref_size));
            tform_matrix = tformEstimate_affine.T;
            vfs_ref = curr_im_reg;
            
            ref_reg(:,:,curr_loop) = vfs_ref;
            AP_print_progress_fraction(curr_loop,n_loops);
        end
        
        % Set make symmetric-average map the master
        symm_vfs = ref_reg(:,:,end);
        symm_vfs_mirror_avg = (symm_vfs + side_sign*fliplr(symm_vfs))./2;
        
        master_vfs = symm_vfs_mirror_avg;
        
        figure;
        imagesc(master_vfs);
        axis image off
        colormap(brewermap([],'*RdBu'));
        caxis([-1,1]);
        title('Master VFS');
        
        % Save the master
        confirm_save = input(['Save new master VFS? (y/n): '],'s');
        if strcmp(confirm_save,'y')
            [master_vfs_filename,master_vfs_path] = ...
                uiputfile([alignment_path filesep '*.mat'],'Save master VFS');
            master_vfs_fn = [master_vfs_path master_vfs_filename];
            save(master_vfs_fn,'master_vfs');
            disp(['Saved new master VFS: ' master_vfs_fn]);
            confirm_align_ccf = strcmp(input(['Align CCF to new master VFS? (y/n): '],'s'),'y');
            if confirm_align_ccf
                AP_vfs_ccf_align;
            end
        else
            disp('Not saved.')
        end
        
    case 'create_submaster'
        %% Make submaster for new type, aligned to master
        
        disp('Creating new submaster VFS...')
        
        % Check inputs
        if ~iscell(im_unaligned)
            error('Expected cell array of unaligned VFS''s');
        end
        
        % Check with user: are L and R hemi same or opposite chirality?
        LR_symmetry = questdlg('Are left/right hemisphere colors same or opposite?','Choose VFS type','symmetric','opposite','symmetric');
        switch LR_symmetry
            case 'symmetric'
                side_sign = 1;
            case 'opposite'
                side_sign = -1;
        end
        
        % Pad and concatenate unaligned VFS's
        [vfs_y,vfs_x] = cellfun(@size,im_unaligned);
        
        vfs_y_max = max(vfs_y);
        vfs_x_max = max(vfs_x);
        
        vfs_unaligned_pad = ...
            cell2mat(reshape(cellfun(@(vfs) padarray(vfs, ...
            [vfs_y_max - size(vfs,1),vfs_x_max - size(vfs,2)],0,'post'), ...
            im_unaligned,'uni',false),1,1,[]));
        
        % Set alignment parameters     
        [optimizer, metric] = imregconfig('monomodal');
        optimizer = registration.optimizer.OnePlusOneEvolutionary();
        optimizer.MaximumIterations = 200;
        optimizer.GrowthFactor = 1+1e-6;
        optimizer.InitialRadius = 1e-4;
        
        % Iterate averaging the VFS and the aligning
        disp('Iterating aligning average VFS...')
        
        vfs_ref = nanmean(vfs_unaligned_pad,3);
        ref_size = size(vfs_ref);
        
        n_loops = 5;
        for curr_loop = 1:n_loops
            
            vfs_aligned = nan(ref_size(1),ref_size(2),length(im_unaligned));
            for curr_animal = 1:length(im_unaligned)
                tformEstimate_affine = imregtform(im_unaligned{curr_animal}, ...
                    vfs_ref,'affine',optimizer,metric,'PyramidLevels',5);
                curr_im_reg = imwarp(im_unaligned{curr_animal}, ...
                    tformEstimate_affine,'Outputview',imref2d(ref_size));
                tform_matrix{curr_animal} = tformEstimate_affine.T;
                vfs_aligned(:,:,curr_animal) = curr_im_reg;
            end
            
            vfs_ref = nanmean(vfs_aligned,3);
            AP_print_progress_fraction(curr_loop,n_loops);
        end
        
        % Symmetrize the average aligned VFS
        disp('Symmetrizing aligned average VFS...')
        
        n_loops = 15;
        ref_size = size(vfs_ref);
        
        ref_reg = nan(size(vfs_ref,1),size(vfs_ref,2),n_loops);
        for curr_loop = 1:n_loops
            
            ref_im_symm = (vfs_ref + side_sign*fliplr(vfs_ref))./2;
   
            tformEstimate_affine = imregtform(vfs_ref,ref_im_symm,'rigid', ...
                optimizer,metric,'PyramidLevels',5);
            curr_im_reg = imwarp(vfs_ref,tformEstimate_affine,'Outputview',imref2d(ref_size));
            tform_matrix = tformEstimate_affine.T;
            vfs_ref = curr_im_reg;
            
            ref_reg(:,:,curr_loop) = vfs_ref;
            AP_print_progress_fraction(curr_loop,n_loops);
        end
        
        % Set make symmetric-average map the master
        symm_vfs = ref_reg(:,:,end);
        symm_vfs_mirror_avg = (symm_vfs + side_sign*fliplr(symm_vfs))./2;
                
        % Load the master, align the submaster
        alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment';
        master_vfs_fn = [alignment_path filesep 'master_vfs.mat'];
        load(master_vfs_fn);
        % (if hemispheres opposite, flip left side)
        if strcmp(LR_symmetry,'opposite')
           master_vfs(:,1:round(size(master_vfs,2)/2)) = ...
               master_vfs(:,1:round(size(master_vfs,2)/2))*-1;
        end
                
        tformEstimate_affine = imregtform(symm_vfs_mirror_avg, ...
            master_vfs,'affine',optimizer,metric,'PyramidLevels',5);
        submaster_vfs = imwarp(symm_vfs_mirror_avg, ...
            tformEstimate_affine,'Outputview',imref2d(size(master_vfs)));
        
        % Plot submaster and alignment
        figure; 
        subplot(1,3,1);
        imagesc(master_vfs);
        axis image off
        colormap(brewermap([],'*RdBu'));
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        title('Master');
        subplot(1,3,2);
        imagesc(submaster_vfs);
        axis image off
        colormap(brewermap([],'*RdBu'));
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        title('New aligned submaster');
        subplot(1,3,3);
        imshowpair(abs(master_vfs),abs(submaster_vfs)); 
        title('Master/submaster alignment');
        
        % Save the submaster
        confirm_save = input(['Save new submaster VFS? (y/n): '],'s');
        if strcmp(confirm_save,'y')
            [master_vfs_filename,master_vfs_path] = ...
                uiputfile([alignment_path filesep '*.mat'],'Save submaster VFS');
            master_vfs_fn = [master_vfs_path master_vfs_filename];
            save(master_vfs_fn,'submaster_vfs');
            disp(['Saved new submaster VFS: ' master_vfs_fn]);
        else
            disp('Not saved.')
        end
        

    case 'new_animal'
        %% Align animal to master VFS
        
        if ~exist('master_align') || isempty(master_align)
            % Align master VFS by default          
            disp('Aligning animal VFS to master VFS...')
            % Load master VFS
            master_vfs_fn = [alignment_path filesep 'master_vfs.mat'];
            load(master_vfs_fn);
            master_align = master_vfs;
        else
            disp('Aligning animal image to master image')
        end
  
        % Align animal image to master image
        ref_size = size(master_align);
        
        [optimizer, metric] = imregconfig('monomodal');
        optimizer = registration.optimizer.OnePlusOneEvolutionary();
        optimizer.MaximumIterations = 200;
        optimizer.GrowthFactor = 1+1e-6;
        optimizer.InitialRadius = 1e-4;
        
        tformEstimate_affine = imregtform(im_unaligned,master_align,'affine',optimizer,metric);
        im_aligned = imwarp(im_unaligned,tformEstimate_affine,'Outputview',imref2d(ref_size));
        tform_matrix = tformEstimate_affine.T;
        
        figure;
        subplot(2,2,1);
        imshowpair(master_align,im_unaligned);
        title('Unaligned');
        subplot(2,2,2);
        imshowpair(master_align,im_aligned);
        title('Aligned');
        subplot(2,2,3);
        imagesc(master_align);
        axis image off
        caxis([-1,1]);
        colormap(brewermap([],'*RdBu'));
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        title('Master')
        subplot(2,2,4);
        imagesc(im_aligned);
        axis image off
        caxis([-1,1]);
        colormap(brewermap([],'*RdBu'));
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        title('Aligned')
        
        % Save transform matrix into structure
        curr_animal_idx = strcmp(animal,{wf_tform.animal});
        if isempty(curr_animal_idx) || ~any(curr_animal_idx)
            % (if not already extant, make new)
            curr_animal_idx = length(wf_tform) + 1;
            confirm_save = true;
        else
            % (if extant, prompt to overwrite)
            confirm_save = strcmp(input(['Overwrite animal transform for ' animal '? (y/n): '],'s'),'y');
        end
        if confirm_save
            wf_tform(curr_animal_idx).animal = animal;
            wf_tform(curr_animal_idx).animal_tform = tform_matrix;
            wf_tform(curr_animal_idx).master_ref_size = ref_size;
            save(wf_tform_fn,'wf_tform');
            disp(['Saved new animal alignment for ' animal '.'])
        else
            disp('Not overwriting.')
        end
        

        
        
    otherwise
        %% Apply alignments to data
        
        % Find animal and day index within wf_tform structure
        curr_animal_idx = strcmp(animal,{wf_tform.animal});
        if ~any(curr_animal_idx)
           error(['No alignments found for ' animal]);
        end
        
        curr_day_idx = strcmp(day,wf_tform(curr_animal_idx).day);
        if ~any(curr_day_idx) && ~strcmp(align_type,'animal_only')
           error(['No ' animal ' alignment found for ' day]);
        end
             
        switch align_type
            case 'day_only'
                % Align day only (ignore animal, even if present)
                curr_tform = wf_tform(curr_animal_idx).day_tform{curr_day_idx};
                ref_size = wf_tform(curr_animal_idx).ref_size;
            case 'animal_only'
                % Align animal only (used if already day-aligned)
                curr_tform = wf_tform(curr_animal_idx).animal_tform;
                ref_size = wf_tform(curr_animal_idx).master_ref_size;
            otherwise
                % Apply both day and animal alignments (if available)
                curr_day_tform = wf_tform(curr_animal_idx).day_tform{curr_day_idx};
                if ~isempty(wf_tform(curr_animal_idx).animal_tform)
                    curr_animal_tform = wf_tform(curr_animal_idx).animal_tform;
                    curr_tform = curr_day_tform*curr_animal_tform;
                    ref_size = wf_tform(curr_animal_idx).master_ref_size;
                else
                    warning([animal ' ' day ': No animal alignment']);
                    curr_tform = curr_day_tform;
                    ref_size = wf_tform(curr_animal_idx).ref_size;
                end               
        end
            
        % Transform image, cast to input data type
        tform = affine2d;
        tform.T = curr_tform;
        im_aligned = cast(imwarp(im_unaligned,tform,'Outputview',imref2d(ref_size)), ...
            class(im_unaligned));
        

end











