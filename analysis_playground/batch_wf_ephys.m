%% Create widefield ROIs

alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment';
load([alignment_path filesep 'animal_wf_tform']);
im_size = animal_wf_tform(1).im_size;

roi_areas = {'V1','AM','RSPp','FRm','FRl','SMl','SMf'};

wf_roi = struct('area',cell(length(roi_areas),2),'mask',cell(length(roi_areas),2));

% (ideally this should mirror, but that'd involve a transform etc)
curr_roi = 1;
for curr_hemi = {'L','R'};
    for curr_area = roi_areas;
        curr_roi_name = [curr_area{:} '_' curr_hemi{:}];
        disp(curr_roi_name);
        wf_roi(curr_roi).area = curr_roi_name;
        [~,wf_roi(curr_roi).mask] = AP_svd_roi(nan(im_size),[],'master');
        curr_roi = curr_roi + 1;
    end
end

wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
save(wf_roi_fn,'wf_roi');
disp('Saved new widefield ROIs');

%% Batch widefield area boundaries

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = false;
    
    batch_vars = struct;
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        
        % Get V covariance
        Ur = reshape(U, size(U,1)*size(U,2),[]); % P x S
        covV = cov(fV'); % S x S % this is the only one that takes some time really
        varP = dot((Ur*covV)', Ur'); % 1 x P
        
        ySize = size(U,1); xSize = size(U,2);
        
        px_spacing = 20;
        use_y = 1:px_spacing:size(U,1);
        use_x = 1:px_spacing:size(U,2);
        corr_map = cell(length(use_y),length(use_x));
        for curr_x_idx = 1:length(use_x)
            curr_x = use_x(curr_x_idx);
            for curr_y_idx = 1:length(use_y)
                curr_y = use_y(curr_y_idx);
                
                pixel = [curr_y,curr_x];
                pixelInd = sub2ind([ySize, xSize], pixel(1), pixel(2));
                
                covP = Ur(pixelInd,:)*covV*Ur'; % 1 x P
                stdPxPy = varP(pixelInd).^0.5 * varP.^0.5; % 1 x P
                corrMat = reshape(covP./stdPxPy,ySize,xSize); % 1 x P
                
                corr_map{curr_y_idx,curr_x_idx} = corrMat;
            end
        end
        
        % Correlation map edge detection
        corr_map_norm = mat2gray(cat(3,corr_map{:}),[0.5,1]);
        corr_map_edge = imgaussfilt(corr_map_norm,5)-imgaussfilt(corr_map_norm,20);
        corr_edges = nanmean(corr_map_edge,3);
        
        batch_vars.corr_edges{curr_day} = corr_edges;
                     
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
        
    % Align images from batch processing
    days = {experiments.day};
    corr_edges_aligned = AP_align_widefield(animal,days,batch_vars.corr_edges);
    
    wf_borders_fig = figure('Name',animal);
    imagesc(nanmean(corr_edges_aligned,3));
    axis image off; colormap(gray); caxis([0,0.05])
    
    fn = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_borders' filesep animal '_wf_borders'];
    saveas(wf_borders_fig,fn);
   
    AP_print_progress_fraction(curr_animal,length(animals))
    
end

disp('Finished batch.')


%% Batch widefield responses to passive stim

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% protocol = 'stimKalatsky';
protocol = 'AP_choiceWorldStimPassive';

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    % Skip if this animal doesn't have this experiment
    if isempty(experiments)
        continue
    end
    
    disp(animal);
    
    experiments = experiments([experiments.imaging]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = false;
    
    batch_vars = struct;
    for curr_day = 1:length(experiments);
        
        % Load the experiment
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;        
        AP_load_experiment
             
        % Set options
        surround_window = [-0.5,5];
        framerate = 1./median(diff(frame_t));
        surround_samplerate = 1/(framerate*1);
        t_surround = surround_window(1):surround_samplerate:surround_window(2);
        
        % Average (time course) responses
        conditions = unique(stimIDs);
        im_stim = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        for curr_condition_idx = 1:length(conditions)
            curr_condition = conditions(curr_condition_idx);
            
            use_stims = find(stimIDs == curr_condition);
            use_stim_onsets = stimOn_times(use_stims(2:end));
            use_stim_onsets([1,end]) = [];
            
            stim_surround_times = bsxfun(@plus, use_stim_onsets(:), t_surround);
            peri_stim_v = permute(mean(interp1(frame_t,fVdf',stim_surround_times),1),[3,2,1]);
            
            im_stim(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
        end
        
        batch_vars.im_stim{curr_day} = im_stim;        
        
        % Prepare for next loop
        AP_print_progress_fraction(curr_day,length(experiments))
        clearvars -except animals protocol curr_animal experiments curr_day animal batch_vars load_parts
        
    end
    
    %     % Get ddf
    %     ddf_im = batch_vars.im_stim;
    %     for curr_day = 1:length(experiments)
    %         curr_im = ddf_im{curr_day};
    %         curr_im(isnan(curr_im)) = 0;
    %         curr_im = imgaussfilt(diff(curr_im,[],3),2);
    %         curr_im(curr_im < 0) = 0;
    %         ddf_im{curr_day} = curr_im;
    %     end
    
    % Align
    days = {experiments.day};
    im_aligned = AP_align_widefield(animal,days,batch_vars.im_stim);
    im_aligned_average = nanmean(im_aligned,5);
    
    %     % Plot
    %     surround_window = [-0.5,5];
    %     t_surround = linspace(surround_window(1),surround_window(2),size(im_aligned_average,3));
    %     AP_image_scroll(im_aligned_average,t_surround);
    %     colormap(gray); axis image off;
    %     title([animal ': passive stimuli']);
    
    % Save
    save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\passive'];
    save([save_path filesep animal '_' protocol],'im_aligned_average');    
    
    disp(['Finished ' animal]);
    
end

disp('Finished batch');

%% Batch widefield choiceworld (stim onset)

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);    
    
    disp(animal);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = false;
    
    batch_vars = struct;
    
    % (initialize these because iteratively added below)
    batch_vars.im_stim_earlymove_hit = [];
    batch_vars.im_stim_latemove_hit = [];
    
    batch_vars.im_stim_earlymove_miss = [];
    batch_vars.im_stim_latemove_miss = [];
    
    batch_vars.n_im_stim_earlymove_hit = [];
    batch_vars.n_im_stim_latemove_hit = [];
    
    batch_vars.n_im_stim_earlymove_miss = [];
    batch_vars.n_im_stim_latemove_miss = [];

    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        conditions = unique(block.events.sessionPerformanceValues(1,:));

        % Define trials to use
        n_trials = length(block.paramsValues);
        trial_conditions = signals_events.trialSideValues(1:n_trials).*signals_events.trialContrastValues(1:n_trials);
        trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
        stim_to_move = padarray(wheel_move_time - stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
        stim_to_feedback = padarray(signals_events.responseTimes, ...
            [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
            padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');    
      
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Set options
        surround_window = [-0.5,5];
        framerate = 1./median(diff(frame_t));
        surround_samplerate = 1/(framerate*1);
        t_surround = surround_window(1):surround_samplerate:surround_window(2);
        
        % Average (time course) responses
        im_stim_earlymove_hit = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        im_stim_latemove_hit = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        
        im_stim_earlymove_miss = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        im_stim_latemove_miss = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        
        for curr_condition_idx = 1:length(conditions)
            curr_condition = conditions(curr_condition_idx);
            
            curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == 1 & stim_to_move < 0.5;
            curr_align = stimOn_times(curr_trials);           
            if length(curr_align) > 5
                curr_surround_times = bsxfun(@plus, curr_align(:), t_surround);
                peri_stim_v = permute(mean(interp1(frame_t,fVdf',curr_surround_times),1),[3,2,1]);
                im_stim_earlymove_hit(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
            end      
            
            curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == 1 & stim_to_move >= 0.5;
            curr_align = stimOn_times(curr_trials);           
            if length(curr_align) > 5
                curr_surround_times = bsxfun(@plus, curr_align(:), t_surround);
                peri_stim_v = permute(mean(interp1(frame_t,fVdf',curr_surround_times),1),[3,2,1]);
                im_stim_latemove_hit(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
            end      
                        
            curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == -1 & stim_to_move < 0.5;
            curr_align = stimOn_times(curr_trials);  
            if length(curr_align) > 5
                curr_surround_times = bsxfun(@plus, curr_align(:), t_surround);
                peri_stim_v = permute(mean(interp1(frame_t,fVdf',curr_surround_times),1),[3,2,1]);
                im_stim_earlymove_miss(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
            end
            
            curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == -1 & stim_to_move >= 0.5;
            curr_align = stimOn_times(curr_trials);  
            if length(curr_align) > 5
                curr_surround_times = bsxfun(@plus, curr_align(:), t_surround);
                peri_stim_v = permute(mean(interp1(frame_t,fVdf',curr_surround_times),1),[3,2,1]);
                im_stim_latemove_miss(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
            end
            
        end
        
        % Align and average as it goes, otherwise the variable's too big
        im_stim_earlymove_hit_align = AP_align_widefield(animal,day,im_stim_earlymove_hit);
        im_stim_latemove_hit_align = AP_align_widefield(animal,day,im_stim_latemove_hit);
        
        im_stim_earlymove_miss_align = AP_align_widefield(animal,day,im_stim_earlymove_miss);
        im_stim_latemove_miss_align = AP_align_widefield(animal,day,im_stim_latemove_miss);
        
        batch_vars.im_stim_earlymove_hit = nansum(cat(5,batch_vars.im_stim_earlymove_hit,im_stim_earlymove_hit_align),5);
        batch_vars.im_stim_latemove_hit = nansum(cat(5,batch_vars.im_stim_latemove_hit,im_stim_latemove_hit_align),5);
        
        batch_vars.im_stim_earlymove_miss = nansum(cat(5,batch_vars.im_stim_earlymove_miss,im_stim_earlymove_miss_align),5);
        batch_vars.im_stim_latemove_miss = nansum(cat(5,batch_vars.im_stim_latemove_miss,im_stim_latemove_miss_align),5);
        
        % Count conditions to divide at end
        batch_vars.n_im_stim_earlymove_hit = ...
            sum(cat(5,batch_vars.n_im_stim_earlymove_hit,any(any(any(im_stim_earlymove_hit_align,1),2),3)),5);
        batch_vars.n_im_stim_latemove_hit = ...
            sum(cat(5,batch_vars.n_im_stim_latemove_hit,any(any(any(im_stim_latemove_hit_align,1),2),3)),5);
        
        batch_vars.n_im_stim_earlymove_miss = ...
            sum(cat(5,batch_vars.n_im_stim_earlymove_miss,any(any(any(im_stim_earlymove_miss_align,1),2),3)),5);
        batch_vars.n_im_stim_latemove_miss = ...
            sum(cat(5,batch_vars.n_im_stim_latemove_miss,any(any(any(im_stim_latemove_miss_align,1),2),3)),5);
        
        % Prepare for next loop
        AP_print_progress_fraction(curr_day,length(experiments))
        clearvars -except animals protocol curr_animal experiments curr_day animal batch_vars load_parts
        
    end
    
    % Divide sum to get average
    im_stim_earlymove_hit_avg = bsxfun(@rdivide,batch_vars.im_stim_earlymove_hit,batch_vars.n_im_stim_earlymove_hit);
    im_stim_latemove_hit_avg = bsxfun(@rdivide,batch_vars.im_stim_latemove_hit,batch_vars.n_im_stim_latemove_hit);
    
    im_stim_earlymove_miss_avg = bsxfun(@rdivide,batch_vars.im_stim_earlymove_miss,batch_vars.n_im_stim_earlymove_miss);
    im_stim_latemove_miss_avg = bsxfun(@rdivide,batch_vars.im_stim_latemove_miss,batch_vars.n_im_stim_latemove_miss);
    
    % Save
    save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
    save([save_path filesep animal '_im_stim_earlymove_hit_avg'],'im_stim_earlymove_hit_avg','-v7.3');
    save([save_path filesep animal '_im_stim_latemove_hit_avg'],'im_stim_latemove_hit_avg','-v7.3');
    save([save_path filesep animal '_im_stim_earlymove_miss_avg'],'im_stim_earlymove_miss_avg','-v7.3');
    save([save_path filesep animal '_im_stim_latemove_miss_avg'],'im_stim_latemove_miss_avg','-v7.3');
    
    disp(['Finished ' animal]);
    
end

disp('Finished batch.')
warning('This uses -v7.3 and therefore compresses data, switch to dat in the future');

%% Batch widefield choiceworld (move onset)

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);    
    
    disp(animal);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = false;
    
    batch_vars = struct;
    
    % (initialize these because iteratively added below)
    batch_vars.im_move_earlymove_hit = [];
    batch_vars.im_move_latemove_hit = [];
    
    batch_vars.im_move_earlymove_miss = [];
    batch_vars.im_move_latemove_miss = [];
    
    batch_vars.n_im_move_earlymove_hit = [];
    batch_vars.n_im_move_latemove_hit = [];
    
    batch_vars.n_im_move_earlymove_miss = [];
    batch_vars.n_im_move_latemove_miss = [];

    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        conditions = unique(block.events.sessionPerformanceValues(1,:));

        % Define trials to use
        n_trials = length(block.paramsValues);
        trial_conditions = signals_events.trialSideValues(1:n_trials).*signals_events.trialContrastValues(1:n_trials);
        trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
        stim_to_move = padarray(wheel_move_time - stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
        stim_to_feedback = padarray(signals_events.responseTimes, ...
            [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
            padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');    
      
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Set options
        surround_window = [-0.5,2];
        upsample_factor = 3;
        
        framerate = 1./median(diff(frame_t));
        surround_samplerate = 1/(framerate*upsample_factor);
        t_surround = surround_window(1):surround_samplerate:surround_window(2);
        
        % Average (time course) responses
        im_move_earlymove_hit = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        im_move_latemove_hit = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        
        im_move_earlymove_miss = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        im_move_latemove_miss = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        
        for curr_condition_idx = 1:length(conditions)
            curr_condition = conditions(curr_condition_idx);
            
            curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == 1 & stim_to_move < 0.5;
            curr_align = wheel_move_time(curr_trials);                 
            if length(curr_align) > 5
                curr_surround_times = bsxfun(@plus, curr_align(:), t_surround);
                peri_stim_v = permute(mean(interp1(frame_t,fVdf',curr_surround_times),1),[3,2,1]);
                im_move_earlymove_hit(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
            end      
            
            curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == 1 & stim_to_move >= 0.5;
            curr_align = wheel_move_time(curr_trials);                 
            if length(curr_align) > 5
                curr_surround_times = bsxfun(@plus, curr_align(:), t_surround);
                peri_stim_v = permute(mean(interp1(frame_t,fVdf',curr_surround_times),1),[3,2,1]);
                im_move_latemove_hit(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
            end      
                        
            curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == -1 & stim_to_move < 0.5;
            curr_align = wheel_move_time(curr_trials);                 
            if length(curr_align) > 5
                curr_surround_times = bsxfun(@plus, curr_align(:), t_surround);
                peri_stim_v = permute(mean(interp1(frame_t,fVdf',curr_surround_times),1),[3,2,1]);
                im_move_earlymove_miss(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
            end
            
            curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == -1 & stim_to_move >= 0.5;
            curr_align = wheel_move_time(curr_trials);                 
            if length(curr_align) > 5
                curr_surround_times = bsxfun(@plus, curr_align(:), t_surround);
                peri_stim_v = permute(mean(interp1(frame_t,fVdf',curr_surround_times),1),[3,2,1]);
                im_move_latemove_miss(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
            end
            
        end
        
        % Align and average as it goes, otherwise the variable's too big
        im_move_earlymove_hit_align = AP_align_widefield(animal,day,im_move_earlymove_hit);
        im_move_latemove_hit_align = AP_align_widefield(animal,day,im_move_latemove_hit);
        
        im_move_earlymove_miss_align = AP_align_widefield(animal,day,im_move_earlymove_miss);
        im_move_latemove_miss_align = AP_align_widefield(animal,day,im_move_latemove_miss);
        
        batch_vars.im_move_earlymove_hit = nansum(cat(5,batch_vars.im_move_earlymove_hit,im_move_earlymove_hit_align),5);
        batch_vars.im_move_latemove_hit = nansum(cat(5,batch_vars.im_move_latemove_hit,im_move_latemove_hit_align),5);
        
        batch_vars.im_move_earlymove_miss = nansum(cat(5,batch_vars.im_move_earlymove_miss,im_move_earlymove_miss_align),5);
        batch_vars.im_move_latemove_miss = nansum(cat(5,batch_vars.im_move_latemove_miss,im_move_latemove_miss_align),5);
        
        % Count conditions to divide at end
        batch_vars.n_im_move_earlymove_hit = ...
            sum(cat(5,batch_vars.n_im_move_earlymove_hit,any(any(any(im_move_earlymove_hit_align,1),2),3)),5);
        batch_vars.n_im_move_latemove_hit = ...
            sum(cat(5,batch_vars.n_im_move_latemove_hit,any(any(any(im_move_latemove_hit_align,1),2),3)),5);
        
        batch_vars.n_im_move_earlymove_miss = ...
            sum(cat(5,batch_vars.n_im_move_earlymove_miss,any(any(any(im_move_earlymove_miss_align,1),2),3)),5);
        batch_vars.n_im_move_latemove_miss = ...
            sum(cat(5,batch_vars.n_im_move_latemove_miss,any(any(any(im_move_latemove_miss_align,1),2),3)),5);
        
        % Prepare for next loop
        AP_print_progress_fraction(curr_day,length(experiments))
        clearvars -except animals protocol curr_animal experiments curr_day animal batch_vars load_parts
        
    end
    
    % Divide sum to get average
    im_move_earlymove_hit_avg = bsxfun(@rdivide,batch_vars.im_move_earlymove_hit,batch_vars.n_im_move_earlymove_hit);
    im_move_latemove_hit_avg = bsxfun(@rdivide,batch_vars.im_move_latemove_hit,batch_vars.n_im_move_latemove_hit);
    
    im_move_earlymove_miss_avg = bsxfun(@rdivide,batch_vars.im_move_earlymove_miss,batch_vars.n_im_move_earlymove_miss);
    im_move_latemove_miss_avg = bsxfun(@rdivide,batch_vars.im_move_latemove_miss,batch_vars.n_im_move_latemove_miss);
    
    % Save
    save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
    save([save_path filesep animal '_im_move_earlymove_hit_avg'],'im_move_earlymove_hit_avg','-v7.3');
    save([save_path filesep animal '_im_move_latemove_hit_avg'],'im_move_latemove_hit_avg','-v7.3');
    save([save_path filesep animal '_im_move_earlymove_miss_avg'],'im_move_earlymove_miss_avg','-v7.3');
    save([save_path filesep animal '_im_move_latemove_miss_avg'],'im_move_latemove_miss_avg','-v7.3');
    
    disp(['Finished ' animal]);
    
end

disp('Finished batch.')
warning('This uses -v7.3 and therefore compresses data, switch to dat in the future');

%% Batch widefield choiceworld (stim onset) !UPSAMPLED: CAN'T DO ALL AT ONCE!

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);    
    
    disp(animal);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = false;
    
    batch_vars = struct;
    
    % (initialize these because iteratively added below)
%     batch_vars.im_stim_earlymove_hit = [];
    batch_vars.im_stim_latemove_hit = [];
    
%     batch_vars.im_stim_earlymove_miss = [];
%     batch_vars.im_stim_latemove_miss = [];
    
%     batch_vars.n_im_stim_earlymove_hit = [];
    batch_vars.n_im_stim_latemove_hit = [];
    
%     batch_vars.n_im_stim_earlymove_miss = [];
%     batch_vars.n_im_stim_latemove_miss = [];

    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        conditions = unique(block.events.sessionPerformanceValues(1,:));

        % Define trials to use
        n_trials = length(block.paramsValues);
        trial_conditions = signals_events.trialSideValues(1:n_trials).*signals_events.trialContrastValues(1:n_trials);
        trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
        stim_to_move = padarray(wheel_move_time - stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
        stim_to_feedback = padarray(signals_events.responseTimes, ...
            [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
            padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');    
      
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Set options
        surround_window = [-0.5,2];
        upsample_factor = 3;
        
        framerate = 1./median(diff(frame_t));
        surround_samplerate = 1/(framerate*upsample_factor);
        t_surround = surround_window(1):surround_samplerate:surround_window(2);
        
        % Average (time course) responses
%         im_stim_earlymove_hit = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        im_stim_latemove_hit = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        
%         im_stim_earlymove_miss = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
%         im_stim_latemove_miss = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        
        for curr_condition_idx = 1:length(conditions)
            curr_condition = conditions(curr_condition_idx);
            
%             curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == 1 & stim_to_move < 0.5;
%             curr_align = stimOn_times(curr_trials);           
%             if length(curr_align) > 5
%                 curr_surround_times = bsxfun(@plus, curr_align(:), t_surround);
%                 peri_stim_v = permute(mean(interp1(frame_t,fVdf',curr_surround_times),1),[3,2,1]);
%                 im_stim_earlymove_hit(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
%             end      
            
            curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == 1 & stim_to_move >= 0.5;
            curr_align = stimOn_times(curr_trials);           
            if length(curr_align) > 5
                curr_surround_times = bsxfun(@plus, curr_align(:), t_surround);
                peri_stim_v = permute(mean(interp1(frame_t,fVdf',curr_surround_times),1),[3,2,1]);
                im_stim_latemove_hit(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
            end      
%                         
%             curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == -1 & stim_to_move < 0.5;
%             curr_align = stimOn_times(curr_trials);  
%             if length(curr_align) > 5
%                 curr_surround_times = bsxfun(@plus, curr_align(:), t_surround);
%                 peri_stim_v = permute(mean(interp1(frame_t,fVdf',curr_surround_times),1),[3,2,1]);
%                 im_stim_earlymove_miss(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
%             end
%             
%             curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == -1 & stim_to_move >= 0.5;
%             curr_align = stimOn_times(curr_trials);  
%             if length(curr_align) > 5
%                 curr_surround_times = bsxfun(@plus, curr_align(:), t_surround);
%                 peri_stim_v = permute(mean(interp1(frame_t,fVdf',curr_surround_times),1),[3,2,1]);
%                 im_stim_latemove_miss(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
%             end
            
        end
        
        % Align and average as it goes, otherwise the variable's too big
%         im_stim_earlymove_hit_align = AP_align_widefield(animal,day,im_stim_earlymove_hit);
        im_stim_latemove_hit_align = AP_align_widefield(animal,day,im_stim_latemove_hit);
%         
%         im_stim_earlymove_miss_align = AP_align_widefield(animal,day,im_stim_earlymove_miss);
%         im_stim_latemove_miss_align = AP_align_widefield(animal,day,im_stim_latemove_miss);
        
%         batch_vars.im_stim_earlymove_hit = nansum(cat(5,batch_vars.im_stim_earlymove_hit,im_stim_earlymove_hit_align),5);
        batch_vars.im_stim_latemove_hit = nansum(cat(5,batch_vars.im_stim_latemove_hit,im_stim_latemove_hit_align),5);
%         
%         batch_vars.im_stim_earlymove_miss = nansum(cat(5,batch_vars.im_stim_earlymove_miss,im_stim_earlymove_miss_align),5);
%         batch_vars.im_stim_latemove_miss = nansum(cat(5,batch_vars.im_stim_latemove_miss,im_stim_latemove_miss_align),5);
        
        % Count conditions to divide at end
%         batch_vars.n_im_stim_earlymove_hit = ...
%             sum(cat(5,batch_vars.n_im_stim_earlymove_hit,any(any(any(im_stim_earlymove_hit_align,1),2),3)),5);
        batch_vars.n_im_stim_latemove_hit = ...
            sum(cat(5,batch_vars.n_im_stim_latemove_hit,any(any(any(im_stim_latemove_hit_align,1),2),3)),5);
%         
%         batch_vars.n_im_stim_earlymove_miss = ...
%             sum(cat(5,batch_vars.n_im_stim_earlymove_miss,any(any(any(im_stim_earlymove_miss_align,1),2),3)),5);
%         batch_vars.n_im_stim_latemove_miss = ...
%             sum(cat(5,batch_vars.n_im_stim_latemove_miss,any(any(any(im_stim_latemove_miss_align,1),2),3)),5);
        
        % Prepare for next loop
        AP_print_progress_fraction(curr_day,length(experiments))
        clearvars -except animals protocol curr_animal experiments curr_day animal batch_vars load_parts
        
    end
    
    % Divide sum to get average
%     im_stim_earlymove_hit_avg = bsxfun(@rdivide,batch_vars.im_stim_earlymove_hit,batch_vars.n_im_stim_earlymove_hit);
    im_stim_latemove_hit_avg = bsxfun(@rdivide,batch_vars.im_stim_latemove_hit,batch_vars.n_im_stim_latemove_hit);
%     
%     im_stim_earlymove_miss_avg = bsxfun(@rdivide,batch_vars.im_stim_earlymove_miss,batch_vars.n_im_stim_earlymove_miss);
%     im_stim_latemove_miss_avg = bsxfun(@rdivide,batch_vars.im_stim_latemove_miss,batch_vars.n_im_stim_latemove_miss);
    
    % Save
    save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
%     save([save_path filesep animal '_im_stim_earlymove_hit_avg'],'im_stim_earlymove_hit_avg','-v7.3');
    save([save_path filesep animal '_im_stim_latemove_hit_avg'],'im_stim_latemove_hit_avg','-v7.3');
%     save([save_path filesep animal '_im_stim_earlymove_miss_avg'],'im_stim_earlymove_miss_avg','-v7.3');
%     save([save_path filesep animal '_im_stim_latemove_miss_avg'],'im_stim_latemove_miss_avg','-v7.3');
    
    disp(['Finished ' animal]);
    
end

disp('Finished batch.')
warning('This uses -v7.3 and therefore compresses data, switch to dat in the future');

%% Batch widefield choiceworld (move onset) !UPSAMPLED: CAN'T DO ALL AT ONCE!

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);    
    
    disp(animal);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = false;
    
    batch_vars = struct;
    
    % (initialize these because iteratively added below)
    batch_vars.im_move_earlymove_hit = [];
%     batch_vars.im_move_latemove_hit = [];
%     
%     batch_vars.im_move_earlymove_miss = [];
%     batch_vars.im_move_latemove_miss = [];
    
    batch_vars.n_im_move_earlymove_hit = [];
%     batch_vars.n_im_move_latemove_hit = [];
%     
%     batch_vars.n_im_move_earlymove_miss = [];
%     batch_vars.n_im_move_latemove_miss = [];

    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        conditions = unique(block.events.sessionPerformanceValues(1,:));

        % Define trials to use
        n_trials = length(block.paramsValues);
        trial_conditions = signals_events.trialSideValues(1:n_trials).*signals_events.trialContrastValues(1:n_trials);
        trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
        stim_to_move = padarray(wheel_move_time - stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
        stim_to_feedback = padarray(signals_events.responseTimes, ...
            [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
            padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');    
      
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Set options
        surround_window = [-0.5,2];
        upsample_factor = 3;
        
        framerate = 1./median(diff(frame_t));
        surround_samplerate = 1/(framerate*upsample_factor);
        t_surround = surround_window(1):surround_samplerate:surround_window(2);
        
        % Average (time course) responses
        im_move_earlymove_hit = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
%         im_move_latemove_hit = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
%         
%         im_move_earlymove_miss = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
%         im_move_latemove_miss = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
        
        for curr_condition_idx = 1:length(conditions)
            curr_condition = conditions(curr_condition_idx);
            
            curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == 1 & stim_to_move < 0.5;
            curr_align = wheel_move_time(curr_trials);                 
            if length(curr_align) > 5
                curr_surround_times = bsxfun(@plus, curr_align(:), t_surround);
                peri_stim_v = permute(mean(interp1(frame_t,fVdf',curr_surround_times),1),[3,2,1]);
                im_move_earlymove_hit(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
            end      
            
%             curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == 1 & stim_to_move >= 0.5;
%             curr_align = wheel_move_time(curr_trials);                 
%             if length(curr_align) > 5
%                 curr_surround_times = bsxfun(@plus, curr_align(:), t_surround);
%                 peri_stim_v = permute(mean(interp1(frame_t,fVdf',curr_surround_times),1),[3,2,1]);
%                 im_move_latemove_hit(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
%             end      
%                         
%             curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == -1 & stim_to_move < 0.5;
%             curr_align = wheel_move_time(curr_trials);                 
%             if length(curr_align) > 5
%                 curr_surround_times = bsxfun(@plus, curr_align(:), t_surround);
%                 peri_stim_v = permute(mean(interp1(frame_t,fVdf',curr_surround_times),1),[3,2,1]);
%                 im_move_earlymove_miss(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
%             end
%             
%             curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == -1 & stim_to_move >= 0.5;
%             curr_align = wheel_move_time(curr_trials);                 
%             if length(curr_align) > 5
%                 curr_surround_times = bsxfun(@plus, curr_align(:), t_surround);
%                 peri_stim_v = permute(mean(interp1(frame_t,fVdf',curr_surround_times),1),[3,2,1]);
%                 im_move_latemove_miss(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
%             end
            
        end
        
        % Align and average as it goes, otherwise the variable's too big
        im_move_earlymove_hit_align = AP_align_widefield(animal,day,im_move_earlymove_hit);
%         im_move_latemove_hit_align = AP_align_widefield(animal,day,im_move_latemove_hit);
%         
%         im_move_earlymove_miss_align = AP_align_widefield(animal,day,im_move_earlymove_miss);
%         im_move_latemove_miss_align = AP_align_widefield(animal,day,im_move_latemove_miss);
        
        batch_vars.im_move_earlymove_hit = nansum(cat(5,batch_vars.im_move_earlymove_hit,im_move_earlymove_hit_align),5);
%         batch_vars.im_move_latemove_hit = nansum(cat(5,batch_vars.im_move_latemove_hit,im_move_latemove_hit_align),5);
%         
%         batch_vars.im_move_earlymove_miss = nansum(cat(5,batch_vars.im_move_earlymove_miss,im_move_earlymove_miss_align),5);
%         batch_vars.im_move_latemove_miss = nansum(cat(5,batch_vars.im_move_latemove_miss,im_move_latemove_miss_align),5);
        
        % Count conditions to divide at end
        batch_vars.n_im_move_earlymove_hit = ...
            sum(cat(5,batch_vars.n_im_move_earlymove_hit,any(any(any(im_move_earlymove_hit_align,1),2),3)),5);
%         batch_vars.n_im_move_latemove_hit = ...
%             sum(cat(5,batch_vars.n_im_move_latemove_hit,any(any(any(im_move_latemove_hit_align,1),2),3)),5);
%         
%         batch_vars.n_im_move_earlymove_miss = ...
%             sum(cat(5,batch_vars.n_im_move_earlymove_miss,any(any(any(im_move_earlymove_miss_align,1),2),3)),5);
%         batch_vars.n_im_move_latemove_miss = ...
%             sum(cat(5,batch_vars.n_im_move_latemove_miss,any(any(any(im_move_latemove_miss_align,1),2),3)),5);
        
        % Prepare for next loop
        AP_print_progress_fraction(curr_day,length(experiments))
        clearvars -except animals protocol curr_animal experiments curr_day animal batch_vars load_parts
        
    end
    
    % Divide sum to get average
    im_move_earlymove_hit_avg = bsxfun(@rdivide,batch_vars.im_move_earlymove_hit,batch_vars.n_im_move_earlymove_hit);
%     im_move_latemove_hit_avg = bsxfun(@rdivide,batch_vars.im_move_latemove_hit,batch_vars.n_im_move_latemove_hit);
%     
%     im_move_earlymove_miss_avg = bsxfun(@rdivide,batch_vars.im_move_earlymove_miss,batch_vars.n_im_move_earlymove_miss);
%     im_move_latemove_miss_avg = bsxfun(@rdivide,batch_vars.im_move_latemove_miss,batch_vars.n_im_move_latemove_miss);
    
    % Save
    save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
    save([save_path filesep animal '_im_move_earlymove_hit_avg'],'im_move_earlymove_hit_avg','-v7.3');
%     save([save_path filesep animal '_im_move_latemove_hit_avg'],'im_move_latemove_hit_avg','-v7.3');
%     save([save_path filesep animal '_im_move_earlymove_miss_avg'],'im_move_earlymove_miss_avg','-v7.3');
%     save([save_path filesep animal '_im_move_latemove_miss_avg'],'im_move_latemove_miss_avg','-v7.3');
    
    disp(['Finished ' animal]);
    
end

disp('Finished batch.')
warning('This uses -v7.3 and therefore compresses data, switch to dat in the future');

%% Batch widefield > striatum maps (upsampled, ddf/f, lambda estimate)

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

protocol = 'vanillaChoiceworld';
% protocol = 'stimSparseNoiseUncorrAsync';
% protocol = 'stimKalatsky';
% protocol = 'AP_choiceWorldStimPassive';

batch_vars = struct;
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    % Skip if this animal doesn't have this experiment
    if isempty(experiments)
        continue
    end
    
    disp(animal);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment       
        
        upsample_factor = 2;
        sample_rate = (1/median(diff(frame_t)))*upsample_factor;
        
        skip_seconds = 60;
        time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        use_svs = 1:50;
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf(use_svs,:),[],2)',time_bin_centers)';
        
        % Estimate lambda using entire striatum       
        use_frames = (frame_t > skip_seconds & frame_t < frame_t(end)-skip_seconds);        
        use_spikes = spike_times_timeline(ismember(spike_templates, ...
            find(templateDepths > str_depth(1) & templateDepths <= str_depth(2))));
        
        [binned_spikes,~,spike_frames] = histcounts(use_spikes,time_bins);
        binned_spikes = single(binned_spikes);
        
        lambda_start = 1e3;
        n_reduce_lambda = 3;
        
        kernel_t = [-0.1,0.1];
        kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
        zs = [false,true];
        cvfold = 5;
        use_frames_idx = find(use_frames);
        
        update_lambda = n_reduce_lambda;
        curr_lambda = lambda_start;
        lambdas = 0;
        explained_var_lambdas = 0;
        
        while update_lambda;
            
            % TO USE dfV
            [~,~,explained_var] = ...
                AP_regresskernel(dfVdf_resample, ...
                binned_spikes,kernel_frames,curr_lambda,zs,cvfold);
            
            lambdas(end+1) = curr_lambda;
            explained_var_lambdas(end+1) = explained_var.total;
            
            if explained_var_lambdas(end) > explained_var_lambdas(end-1)
                curr_lambda = curr_lambda*10;
            else
                lambdas(end) = [];
                explained_var_lambdas(end) = [];
                curr_lambda = curr_lambda/2;
                update_lambda = update_lambda-1;
            end
            
        end
        
        lambda = lambdas(end);         
        
        % Get cortex/striatum regression by depth
               
        % Group striatum depths
        n_depth_groups = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
        [depth_group_n,~,depth_group] = histcounts(spikeDepths,depth_group_edges);
          
        binned_spikes = zeros(n_depth_groups,length(time_bin_centers));
        for curr_depth = 1:n_depth_groups           
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);           
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);            
        end
        
        kernel_t = [-0.3,0.3];
        kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
        zs = [false,true];
        cvfold = 5;      

        % TO USE dfV
        [k,predicted_spikes,explained_var] = ...
            AP_regresskernel(dfVdf_resample, ...
            binned_spikes,kernel_frames,lambda,zs,cvfold);
        
        % Reshape kernel and convert to pixel space
        r = reshape(k,length(use_svs),length(kernel_frames),size(binned_spikes,1));
        
        r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
        for curr_spikes = 1:size(r,3);
            r_px(:,:,:,curr_spikes) = svdFrameReconstruct(Udf(:,:,use_svs),r(:,:,curr_spikes));
        end
            
        % Get center of mass for each pixel
        r_px_max = squeeze(max(r_px,[],3));
        r_px_max_zeronan = r_px_max;
        r_px_max_zeronan(isnan(r_px_max_zeronan)) = 0;
        r_px_max_norm = bsxfun(@rdivide,bsxfun(@minus,r_px_max_zeronan,min(r_px_max_zeronan,[],3)), ...
            max(bsxfun(@minus,r_px_max_zeronan,min(r_px_max_zeronan,[],3)),[],3));
        r_px_com = sum(bsxfun(@times,r_px_max_norm,permute(1:n_depth_groups,[1,3,2])),3)./sum(r_px_max_norm,3);
        
        r_px_weight = max(r_px_max,[],3);
        
        % Store all variables to save
        batch_vars(curr_animal).r_px{curr_day} = r_px;
        batch_vars(curr_animal).r_px_com{curr_day} = r_px_com;
        batch_vars(curr_animal).r_px_weight{curr_day} = r_px_weight;
        batch_vars(curr_animal).explained_var{curr_day} = explained_var.total;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
      
    disp(['Finished ' animal]);
    
end

% Save
save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_ephys'];
save([save_path filesep 'wf_ephys_maps_' protocol],'batch_vars','-v7.3');
warning('saving -v7.3');
disp('Finished batch');

%% Batch cortex > striatum prediction by condition and depth (NL fit)

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';
batch_vars = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    disp(animal);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        
        sample_rate_factor = 1;
        
        % Group multiunit by depth
        n_depth_groups = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
        
        depth_group = discretize(spikeDepths,depth_group_edges);
        depth_groups_used = unique(depth_group);
        depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);
        
        % Skip the first/last n seconds for prediction
        skip_seconds = 60*1;
        
        sample_rate = (1/median(diff(frame_t)))*sample_rate_factor;
        
        time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        binned_spikes = zeros(length(depth_group_edges)-1,length(time_bins)-1);
        for curr_depth = 1:length(depth_group_edges)-1
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        
        stim_conditions = signals_events.trialContrastValues.*signals_events.trialSideValues;
        conditions = unique(block.events.sessionPerformanceValues(1,:));
        
        % (only use trials that were completed and not repeat trials)
        use_trials = any([signals_events.hitValues;signals_events.missValues],1) & ~signals_events.repeatTrialValues;
        align_times = reshape(stimOn_times(use_trials),[],1);
        
        interval_surround = [-0.5,1.5];
        t_surround = interval_surround(1):1/sample_rate:interval_surround(2);
        t_peri_event = bsxfun(@plus,align_times,t_surround);
        
        % Regress spikes from cortex
        use_svs = 1:50;
        kernel_frames = -35:17;
        lambda = 2e5;
        zs = [false,true];
        cvfold = 5;
        
        fVdf_resample = interp1(frame_t,fVdf(use_svs,:)',time_bin_centers)';
        dfVdf_resample = interp1(conv(frame_t,[1,1]/2,'valid'),diff(fVdf(use_svs,:),[],2)',time_bin_centers)';
        
        % [k,predicted_spikes,explained_var] = ...
        %     AP_regresskernel(fVdf_resample, ...
        %     binned_spikes,kernel_frames,lambda,zs,cvfold);
        [k,predicted_spikes,explained_var] = ...
            AP_regresskernel(dfVdf_resample, ...
            binned_spikes,kernel_frames,lambda,zs,cvfold);
        
        binned_spikes_std = std(binned_spikes,[],2);
        binned_spikes_mean = mean(binned_spikes,2);
        predicted_spikes_reranged = bsxfun(@plus,bsxfun(@times,predicted_spikes, ...
            binned_spikes_std),binned_spikes_mean);
        
        % Find nonlinearity for predicted spikes
        nlin_fit = @(b,x) (x+b(1)).^b(2);
        nlin_params = nan(n_depth_groups,2);
        predicted_spikes_nlin = nan(size(predicted_spikes_reranged));
        for curr_depth = 1:n_depth_groups
            
            % Get the median predicted spikes per binned spikes
            [grp_binned_spikes,grp_predicted_spikes]= ...
                grpstats(predicted_spikes_reranged(curr_depth,:),binned_spikes(curr_depth,:),{'gname','median'});
            grp_binned_spikes = cellfun(@str2num,grp_binned_spikes);
            
            % Pick spikes to fit nonlinearity
            fit_spikes_thresh = prctile(binned_spikes(curr_depth,:),99);
            fit_spikes = grp_binned_spikes > 0 & ...
                grp_binned_spikes < fit_spikes_thresh;
            
            % Less than 5 points to fit, don't bother
            if sum(fit_spikes) < 5
                continue
            end
            
            % Fit the nonlinearity
            fit_options = statset('MaxIter',1000);
            beta0 = [grp_predicted_spikes(1),1];
            beta = nlinfit(grp_predicted_spikes(fit_spikes),grp_binned_spikes(fit_spikes),nlin_fit,beta0,fit_options);
            nlin_params(curr_depth,:) = beta;
            
            % Fit model and rectify
            predicted_spikes_rectif = nlin_fit(beta,predicted_spikes_reranged(curr_depth,:));
            predicted_spikes_rectif(imag(predicted_spikes_rectif) ~= 0) = 0;
            predicted_spikes_nlin(curr_depth,:) = predicted_spikes_rectif;
            
%             % Plot model fit
%             figure; hold on;
%             plot(grp_predicted_spikes,grp_binned_spikes,'k')
%             grp_predicted_spikes_nlin = nlin_fit(beta,grp_predicted_spikes);
%             grp_predicted_spikes_nlin(imag(grp_predicted_spikes_nlin) ~= 0) = 0;
%             plot(grp_predicted_spikes_nlin,grp_binned_spikes,'r')
%             line([0,fit_spikes_thresh],[0,fit_spikes_thresh]);
%             xlabel('Predicted spikes')
%             ylabel('Real spikes')
%             legend({'Raw','Nlin'})
        end
        
        if ~isreal(nlin_params)
            warning('Imaginary parameter fits')
        end
        
        % % Get new explained variance (this can't be right...)
        % sse_real_spikes = sum(bsxfun(@minus,binned_spikes,nanmean(binned_spikes,2)).^2,2);
        % sse_total_residual = sum((predicted_spikes_nlin-binned_spikes).^2,2);
        % explained_var_nlin = (sse_real_spikes - sse_total_residual)./sse_real_spikes;
        
        % Align real and predicted spikes to event
        mua_stim_hit = nan(n_depth_groups,length(t_surround),length(conditions));
        mua_stim_hit_pred = nan(n_depth_groups,length(t_surround),length(conditions));
        
        mua_stim_miss = nan(n_depth_groups,length(t_surround),length(conditions));
        mua_stim_miss_pred = nan(n_depth_groups,length(t_surround),length(conditions));
        
        for curr_depth = 1:n_depth_groups
            for curr_condition_idx = 1:length(conditions)
                mua_stim = interp1(time_bin_centers,binned_spikes(curr_depth,:),t_peri_event);
                mua_stim_pred = interp1(time_bin_centers,predicted_spikes_nlin(curr_depth,:),t_peri_event);
                
                curr_trials = signals_events.hitValues(use_trials) == 1 & stim_conditions(use_trials) == conditions(curr_condition_idx);
                if sum(curr_trials) > 5
                    mua_stim_hit(curr_depth,:,curr_condition_idx) = nanmean(mua_stim(curr_trials,:),1);
                    mua_stim_hit_pred(curr_depth,:,curr_condition_idx) = nanmean(mua_stim_pred(curr_trials,:),1);
                end
                
                curr_trials = signals_events.missValues(use_trials) == 1 & stim_conditions(use_trials) == conditions(curr_condition_idx);
                if sum(curr_trials) > 5
                    mua_stim_miss(curr_depth,:,curr_condition_idx) = nanmean(mua_stim(curr_trials,:),1);
                    mua_stim_miss_pred(curr_depth,:,curr_condition_idx) = nanmean(mua_stim_pred(curr_trials,:),1);
                end
            end
        end
        
        % Store variables
        batch_vars(curr_animal).mua_stim_hit(:,:,:,curr_day) = mua_stim_hit;
        batch_vars(curr_animal).mua_stim_hit_pred(:,:,:,curr_day) = mua_stim_hit_pred;
        
        batch_vars(curr_animal).mua_stim_miss(:,:,:,curr_day) = mua_stim_miss;
        batch_vars(curr_animal).mua_stim_miss_pred(:,:,:,curr_day) = mua_stim_miss_pred;
        
        batch_vars(curr_animal).nlin_params(:,:,curr_day) = nlin_params;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal])
    
end

%% Batch cortex > striatum prediction by condition and depth (raw)

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';
batch_vars = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    disp(animal);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        conditions = unique(block.events.sessionPerformanceValues(1,:));

        % Define trials to use
        n_trials = length(block.paramsValues);
        trial_conditions = signals_events.trialSideValues(1:n_trials).*signals_events.trialContrastValues(1:n_trials);
        trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
        stim_to_move = padarray(wheel_move_time - stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
        stim_to_feedback = padarray(signals_events.responseTimes, ...
            [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
            padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');    
      
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
               
        % Group multiunit by depth
        n_depth_groups = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
        
        depth_group = discretize(spikeDepths,depth_group_edges);
        depth_groups_used = unique(depth_group);
        depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);
        
        % Skip the first/last n seconds for prediction
        skip_seconds = 60*1;
        
        sample_rate_factor = 1;
        sample_rate = (1/median(diff(frame_t)))*sample_rate_factor;
        
        time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        binned_spikes = zeros(length(depth_group_edges)-1,length(time_bins)-1);
        for curr_depth = 1:length(depth_group_edges)-1
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        
        stim_conditions = signals_events.trialContrastValues.*signals_events.trialSideValues;
        align_times = reshape(stimOn_times,[],1);
        
        interval_surround = [-0.5,3];
        t_surround = interval_surround(1):1/sample_rate:interval_surround(2);
        t_peri_event = bsxfun(@plus,align_times,t_surround);
        
        % Regress spikes from cortex
        use_svs = 1:50;
        kernel_frames = -35:17;
        lambda = 2e6;
        zs = [false,true];
        cvfold = 5;
        
        fVdf_resample = interp1(frame_t,fVdf(use_svs,:)',time_bin_centers)';
        dfVdf_resample = interp1(conv(frame_t,[1,1]/2,'valid'),diff(fVdf(use_svs,:),[],2)',time_bin_centers)';
        
        % (to use df/f)
        [k,predicted_spikes,explained_var] = ...
            AP_regresskernel(fVdf_resample, ...
            binned_spikes,kernel_frames,lambda,zs,cvfold);
        % (to use ddf/f)
%         [k,predicted_spikes,explained_var] = ...
%             AP_regresskernel(dfVdf_resample, ...
%             binned_spikes,kernel_frames,lambda,zs,cvfold);
        
        binned_spikes_std = std(binned_spikes,[],2);
        binned_spikes_mean = mean(binned_spikes,2);
        predicted_spikes_reranged = bsxfun(@plus,bsxfun(@times,predicted_spikes, ...
            binned_spikes_std),binned_spikes_mean);        
        
        % Align real and predicted spikes to event
        mua_stim_earlymove_hit = nan(n_depth_groups,length(t_surround),length(conditions));
        mua_stim_earlymove_hit_pred = nan(n_depth_groups,length(t_surround),length(conditions));
        
        mua_stim_latemove_hit = nan(n_depth_groups,length(t_surround),length(conditions));
        mua_stim_latemove_hit_pred = nan(n_depth_groups,length(t_surround),length(conditions));
        
        mua_stim_earlymove_miss = nan(n_depth_groups,length(t_surround),length(conditions));
        mua_stim_earlymove_miss_pred = nan(n_depth_groups,length(t_surround),length(conditions));
        
        mua_stim_latemove_miss = nan(n_depth_groups,length(t_surround),length(conditions));
        mua_stim_latemove_miss_pred = nan(n_depth_groups,length(t_surround),length(conditions));
        
        for curr_depth = 1:n_depth_groups
            for curr_condition_idx = 1:length(conditions)
                
                curr_condition = conditions(curr_condition_idx);
                
                mua_stim = interp1(time_bin_centers,binned_spikes(curr_depth,:),t_peri_event);
                mua_stim_pred = interp1(time_bin_centers,predicted_spikes_reranged(curr_depth,:),t_peri_event);
                
                curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == 1 & stim_to_move < 0.5;
                if sum(curr_trials) > 5
                    mua_stim_earlymove_hit(curr_depth,:,curr_condition_idx) = nanmean(mua_stim(curr_trials,:),1);
                    mua_stim_earlymove_hit_pred(curr_depth,:,curr_condition_idx) = nanmean(mua_stim_pred(curr_trials,:),1);
                end
                
                curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == 1 & stim_to_move >= 0.5;
                if sum(curr_trials) > 5
                    mua_stim_latemove_hit(curr_depth,:,curr_condition_idx) = nanmean(mua_stim(curr_trials,:),1);
                    mua_stim_latemove_hit_pred(curr_depth,:,curr_condition_idx) = nanmean(mua_stim_pred(curr_trials,:),1);
                end
                
                curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == -1 & stim_to_move < 0.5;
                if sum(curr_trials) > 5
                    mua_stim_earlymove_miss(curr_depth,:,curr_condition_idx) = nanmean(mua_stim(curr_trials,:),1);
                    mua_stim_earlymove_miss_pred(curr_depth,:,curr_condition_idx) = nanmean(mua_stim_pred(curr_trials,:),1);
                end
                
                curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == -1 & stim_to_move >= 0.5;
                if sum(curr_trials) > 5
                    mua_stim_latemove_miss(curr_depth,:,curr_condition_idx) = nanmean(mua_stim(curr_trials,:),1);
                    mua_stim_latemove_miss_pred(curr_depth,:,curr_condition_idx) = nanmean(mua_stim_pred(curr_trials,:),1);
                end
            end
        end
        
        % Store variables
        batch_vars(curr_animal).mua_stim_earlymove_hit(:,:,:,curr_day) = mua_stim_earlymove_hit;
        batch_vars(curr_animal).mua_stim_earlymove_hit_pred(:,:,:,curr_day) = mua_stim_earlymove_hit_pred;
        
        batch_vars(curr_animal).mua_stim_latemove_hit(:,:,:,curr_day) = mua_stim_latemove_hit;
        batch_vars(curr_animal).mua_stim_latemove_hit_pred(:,:,:,curr_day) = mua_stim_latemove_hit_pred;
        
        batch_vars(curr_animal).mua_stim_earlymove_miss(:,:,:,curr_day) = mua_stim_earlymove_miss;
        batch_vars(curr_animal).mua_stim_earlymove_miss_pred(:,:,:,curr_day) = mua_stim_earlymove_miss_pred;
        
        batch_vars(curr_animal).mua_stim_latemove_miss(:,:,:,curr_day) = mua_stim_latemove_miss;
        batch_vars(curr_animal).mua_stim_latemove_miss_pred(:,:,:,curr_day) = mua_stim_latemove_miss_pred;
                
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal])
    
end

save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
save([save_path filesep 'mua_stim_choiceworld_pred'],'batch_vars');

%% Batch striatum responses to choiceworld

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';
batch_vars = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    disp(animal);
    
    experiments = experiments([experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = false;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment        
        conditions = unique(block.events.sessionPerformanceValues(1,:));
        
        % Define trials to use
        n_trials = length(block.paramsValues);
        trial_conditions = signals_events.trialSideValues(1:n_trials).*signals_events.trialContrastValues(1:n_trials);
        trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
        stim_to_move = padarray(wheel_move_time - stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
        stim_to_feedback = padarray(signals_events.responseTimes, ...
            [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
            padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');    
      
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Group multiunit by depth
        n_depth_groups = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
        depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
        
        depth_group = discretize(spikeDepths,depth_group_edges);
        
        raster_window = [-0.5,3];
        psth_bin_size = 0.001;
        t = raster_window(1):psth_bin_size:raster_window(2);
        t_bins = t(1:end-1) + diff(t);
        
        mua_stim_earlymove_hit = nan(6,length(t_bins),length(conditions));
        mua_stim_latemove_hit = nan(6,length(t_bins),length(conditions));        
        mua_stim_earlymove_miss = nan(6,length(t_bins),length(conditions));
        mua_stim_latemove_miss = nan(6,length(t_bins),length(conditions));
        
        mua_move_earlymove_hit = nan(6,length(t_bins),length(conditions));
        mua_move_latemove_hit = nan(6,length(t_bins),length(conditions));        
        mua_move_earlymove_miss = nan(6,length(t_bins),length(conditions));
        mua_move_latemove_miss = nan(6,length(t_bins),length(conditions));
        
        mua_feedback_earlymove_hit = nan(6,length(t_bins),length(conditions));
        mua_feedback_latemove_hit = nan(6,length(t_bins),length(conditions));        
        mua_feedback_earlymove_miss = nan(6,length(t_bins),length(conditions));
        mua_feedback_latemove_miss = nan(6,length(t_bins),length(conditions));
        
        for curr_depth = 1:n_depth_groups
            
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            
            for curr_condition_idx = 1:length(conditions)
                curr_condition = conditions(curr_condition_idx);                
                
                curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == 1 & stim_to_move < 0.5;
                use_stim_onsets = stimOn_times(curr_trials);
                use_move_onsets = wheel_move_time(curr_trials);
                use_feedback_onsets = signals_events.responseTimes(curr_trials);
                if length(use_stim_onsets) > 5
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_stim_onsets,raster_window,psth_bin_size);
                    mua_stim_earlymove_hit(curr_depth,:,curr_condition_idx) = psth;
                    
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_move_onsets,raster_window,psth_bin_size);
                    mua_move_earlymove_hit(curr_depth,:,curr_condition_idx) = psth;
                    
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_feedback_onsets,raster_window,psth_bin_size);
                    mua_feedback_earlymove_hit(curr_depth,:,curr_condition_idx) = psth;
                end
                
                curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == 1 & stim_to_move >= 0.5;
                use_stim_onsets = stimOn_times(curr_trials);
                use_move_onsets = wheel_move_time(curr_trials);
                use_feedback_onsets = signals_events.responseTimes(curr_trials);
                if length(use_stim_onsets) > 5
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_stim_onsets,raster_window,psth_bin_size);
                    mua_stim_latemove_hit(curr_depth,:,curr_condition_idx) = psth;
                    
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_move_onsets,raster_window,psth_bin_size);
                    mua_move_latemove_hit(curr_depth,:,curr_condition_idx) = psth;
                    
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_feedback_onsets,raster_window,psth_bin_size);
                    mua_feedback_latemove_hit(curr_depth,:,curr_condition_idx) = psth;
                end
                
                curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == -1 & stim_to_move < 0.5;
                use_stim_onsets = stimOn_times(curr_trials);
                use_move_onsets = wheel_move_time(curr_trials);
                use_feedback_onsets = signals_events.responseTimes(curr_trials);
                if length(use_stim_onsets) > 5
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_stim_onsets,raster_window,psth_bin_size);
                    mua_stim_earlymove_miss(curr_depth,:,curr_condition_idx) = psth;
                    
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_move_onsets,raster_window,psth_bin_size);
                    mua_move_earlymove_miss(curr_depth,:,curr_condition_idx) = psth;
                    
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_feedback_onsets,raster_window,psth_bin_size);
                    mua_feedback_earlymove_miss(curr_depth,:,curr_condition_idx) = psth;
                end
                
                curr_trials = use_trials & trial_conditions == curr_condition & trial_outcome == -1 & stim_to_move >= 0.5;
                use_stim_onsets = stimOn_times(curr_trials);
                use_move_onsets = wheel_move_time(curr_trials);
                use_feedback_onsets = signals_events.responseTimes(curr_trials);
                if length(use_stim_onsets) > 5
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_stim_onsets,raster_window,psth_bin_size);
                    mua_stim_latemove_miss(curr_depth,:,curr_condition_idx) = psth;
                    
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_move_onsets,raster_window,psth_bin_size);
                    mua_move_latemove_miss(curr_depth,:,curr_condition_idx) = psth;
                    
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_feedback_onsets,raster_window,psth_bin_size);
                    mua_feedback_latemove_miss(curr_depth,:,curr_condition_idx) = psth;
                end
                
            end
        end
        
        batch_vars(curr_animal).mua_stim_earlymove_hit(:,:,:,curr_day) = mua_stim_earlymove_hit;
        batch_vars(curr_animal).mua_stim_latemove_hit(:,:,:,curr_day) = mua_stim_latemove_hit;       
        batch_vars(curr_animal).mua_stim_earlymove_miss(:,:,:,curr_day) = mua_stim_earlymove_miss;
        batch_vars(curr_animal).mua_stim_latemove_miss(:,:,:,curr_day) = mua_stim_latemove_miss;
        
        batch_vars(curr_animal).mua_move_earlymove_hit(:,:,:,curr_day) = mua_move_earlymove_hit;
        batch_vars(curr_animal).mua_move_latemove_hit(:,:,:,curr_day) = mua_move_latemove_hit;       
        batch_vars(curr_animal).mua_move_earlymove_miss(:,:,:,curr_day) = mua_move_earlymove_miss;
        batch_vars(curr_animal).mua_move_latemove_miss(:,:,:,curr_day) = mua_move_latemove_miss;
        
        batch_vars(curr_animal).mua_feedback_earlymove_hit(:,:,:,curr_day) = mua_feedback_earlymove_hit;
        batch_vars(curr_animal).mua_feedback_latemove_hit(:,:,:,curr_day) = mua_feedback_latemove_hit;       
        batch_vars(curr_animal).mua_feedback_earlymove_miss(:,:,:,curr_day) = mua_feedback_earlymove_miss;
        batch_vars(curr_animal).mua_feedback_latemove_miss(:,:,:,curr_day) = mua_feedback_latemove_miss;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal])
    
end

% Save
save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
save([save_path filesep 'mua_choiceworld'],'batch_vars');

%% Batch striatum responses to passive

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% protocol = 'stimKalatsky';
protocol = 'AP_choiceWorldStimPassive';

batch_vars = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    disp(animal);
    
    experiments = experiments([experiments.ephys]);
    
    % Skip if this animal doesn't have this experiment
    if isempty(experiments)
        continue
    end
    
    load_parts.cam = false;
    load_parts.imaging = false;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        
        conditions = unique(stimIDs);

        % Group multiunit by depth
        n_depth_groups = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
        depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
        
        depth_group = discretize(spikeDepths,depth_group_edges);
        
        raster_window = [-0.5,5];
        psth_bin_size = 0.001;
        t = raster_window(1):psth_bin_size:raster_window(2);
        t_bins = t(1:end-1) + diff(t);
        
        mua_stim = nan(6,length(t_bins),length(conditions));
        for curr_depth = 1:n_depth_groups
            
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            
            for curr_condition_idx = 1:length(conditions)
                curr_condition = conditions(curr_condition_idx);
                
                use_stims = find(stimIDs == curr_condition);
                use_stim_onsets = stimOn_times(use_stims(2:end));
                use_stim_onsets([1,end]) = [];
                
                if length(use_stim_onsets) > 5
                    [psth, bins, rasterX, rasterY] = psthAndBA(curr_spike_times,use_stim_onsets,raster_window,psth_bin_size);
                    mua_stim(curr_depth,:,curr_condition_idx) = psth;
                end
                
            end
        end
        
        batch_vars(curr_animal).mua_stim(:,:,:,curr_day) = mua_stim;

        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal])
    
end

% Save
save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\passive'];
save([save_path filesep 'mua_stim_' protocol],'batch_vars');


%% Batch trial-trial MUA-fluorescence correlation (stim-aligned)
% NOTE: I put a bunch of work into doing this on the GPU instead, but the
% memory-limiting steps required FOR loops and the resulting code was
% barely faster

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

protocol = 'vanillaChoiceworld';
% protocol = 'stimSparseNoiseUncorrAsync';
% protocol = 'stimKalatsky';
% protocol = 'AP_choiceWorldStimPassive';

batch_vars = struct;
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    % Skip if this animal doesn't have this experiment
    if isempty(experiments)
        continue
    end
    
    disp(animal);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        
        % Load widefield ROIs
        wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
        load(wf_roi_fn);
        
        % Group striatum depths
        n_depth_groups = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
        [depth_group_n,~,depth_group] = histcounts(spikeDepths,depth_group_edges);
     
        % Define times to align
        n_trials = length(block.paramsValues);
        trial_conditions = signals_events.trialSideValues(1:n_trials).*signals_events.trialContrastValues(1:n_trials);
        trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
        stim_to_move = padarray(wheel_move_time - stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
        stim_to_feedback = padarray(signals_events.responseTimes, ...
            [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
            padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
        
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_move < 0.5 & ...
            stim_to_feedback < 1.5;
        
        align_times = reshape(stimOn_times(use_trials),[],1);
        
        sample_rate_factor = 3;
        interval_surround = [-0.5,1.5];
        sample_rate = framerate*sample_rate_factor;
        t_surround = interval_surround(1):1/sample_rate:interval_surround(2);
        t_peri_event = bsxfun(@plus,align_times,t_surround);
        
        % Get fluorescence within saved ROIs
        Udf_aligned = AP_align_widefield(animal,day,Udf);
        roi_traces = AP_svd_roi(Udf_aligned,fVdf,[],[],cat(3,wf_roi.mask));
        event_aligned_f = interp1(frame_t,roi_traces',t_peri_event);
        event_aligned_df = interp1(conv(frame_t,[1,1]/2,'valid'),diff(roi_traces,[],2)',t_peri_event);
        event_aligned_df(event_aligned_df < 0) = 0;
        
        % Pull out MUA across depths
        t_peri_event_bins = [t_peri_event - 1/(sample_rate*2), ...
            t_peri_event(:,end) + 1/(sample_rate*2)];
        event_aligned_spikes = nan(size(t_peri_event_bins,1),size(t_peri_event_bins,2)-1,n_depth_groups);
        for curr_depth = 1:n_depth_groups
            use_spikes = spike_times_timeline(depth_group == curr_depth);
            event_aligned_spikes(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(use_spikes,t_peri_event_bins(x,:)),[1:length(align_times)]','uni',false));
        end
        
        % Get wheel movements around all aligned events
        event_aligned_wheel = interp1(conv2(Timeline.rawDAQTimestamps,[0.5,0.5],'valid'), ...
            wheel_velocity,t_peri_event);
        
        % Choice of all events (L/R movement)
        choice = ...
            ((signals_events.trialSideValues(use_trials) == 1 & trial_outcome(use_trials) == -1 | ...
            signals_events.trialSideValues(use_trials) == -1 & trial_outcome(use_trials) == 1) - ...
            (signals_events.trialSideValues(use_trials) == -1 & trial_outcome(use_trials) == -1 | ...
            signals_events.trialSideValues(use_trials) == 1 & trial_outcome(use_trials) == 1))';
        if any(~ismember(choice,[-1,1]))
            error('Non-valid choice')
        end
        
        % Get sig correlations across all time points across modalities
        n_shuff = 1000;
        warning off;
        shuff_idx = ...
            shake(repmat(reshape(1:length(align_times)*length(t_surround), ...
            length(align_times),length(t_surround)),1,1,n_shuff),1);
        warning on ;
        
        % MUA-MUA
        corr_mua_mua = cell(size(event_aligned_spikes,3),size(event_aligned_spikes,3));
        for curr_mua1 = 1:size(event_aligned_spikes,3)
            for curr_mua2 = 1:size(event_aligned_spikes,3)                            
                
                curr_data1 = zscore(event_aligned_spikes(:,:,curr_mua1),[],1);
                curr_data2 = zscore(event_aligned_spikes(:,:,curr_mua2),[],1);               
                curr_data2_shuff = curr_data2(shuff_idx);
                              
                corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
                
                corr_shuff = ...
                    gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                    gpuArray(curr_data2_shuff)))./(length(align_times)-1);
                           
                corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
                corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
                
                corr_mua_mua{curr_mua1,curr_mua2} = corr_real;
                corr_mua_mua{curr_mua1,curr_mua2}(~corr_sig) = 0;
                
            end
        end      
        
        % Fluor-Fluor
        corr_fluor_fluor = cell(size(event_aligned_df,3),size(event_aligned_df,3));
        for curr_fluor1 = 1:size(event_aligned_df,3)
            for curr_fluor2 = 1:size(event_aligned_df,3)
                
                curr_data1 = zscore(event_aligned_df(:,:,curr_fluor1),[],1);
                curr_data2 = zscore(event_aligned_df(:,:,curr_fluor2),[],1);               
                curr_data2_shuff = curr_data2(shuff_idx);
                              
                corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
                
                corr_shuff = ...
                    gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                    gpuArray(curr_data2_shuff)))./(length(align_times)-1);
                
                corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
                corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
                
                corr_fluor_fluor{curr_fluor1,curr_fluor2} = corr_real;
                corr_fluor_fluor{curr_fluor1,curr_fluor2}(~corr_sig) = 0;
                
            end
        end     
        
        % Fluor-MUA
        corr_fluor_mua = cell(size(event_aligned_df,3),size(event_aligned_spikes,3));
        for curr_fluor = 1:size(event_aligned_df,3)
            for curr_mua = 1:size(event_aligned_spikes,3)
                
                curr_data1 = zscore(event_aligned_df(:,:,curr_fluor),[],1);
                curr_data2 = zscore(event_aligned_spikes(:,:,curr_mua),[],1);               
                curr_data2_shuff = curr_data2(shuff_idx);
                              
                corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
                
                corr_shuff = ...
                    gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                    gpuArray(curr_data2_shuff)))./(length(align_times)-1);
                
                corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
                corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
                
                corr_fluor_mua{curr_fluor,curr_mua} = corr_real;
                corr_fluor_mua{curr_fluor,curr_mua}(~corr_sig) = 0;
                
            end
        end      
        
        % MUA-wheel/choice
        corr_mua_wheel = cell(size(event_aligned_spikes,3),1);
        corr_mua_choice = cell(size(event_aligned_spikes,3),1);
        for curr_mua = 1:size(event_aligned_spikes,3)
            
            % Wheel
            curr_data1 = zscore(event_aligned_spikes(:,:,curr_mua),[],1);
            curr_data2 = zscore(event_aligned_wheel,[],1);
            curr_data2_shuff = curr_data2(shuff_idx);
            
            corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
            
            corr_shuff = ...
                gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                gpuArray(curr_data2_shuff)))./(length(align_times)-1);
            
            corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
            corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
            
            corr_mua_wheel{curr_mua} = corr_real;
            corr_mua_wheel{curr_mua}(~corr_sig) = 0;
            
            % Choice
            curr_data1 = zscore(event_aligned_spikes(:,:,curr_mua),[],1);
            curr_data2 = zscore(choice,[],1);
            curr_data2_shuff = curr_data2(shuff_idx(:,1,:));
            
            corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
            
            corr_shuff = ...
                gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                gpuArray(curr_data2_shuff)))./(length(align_times)-1);
            
            corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
            corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
            
            corr_mua_choice{curr_mua} = corr_real;
            corr_mua_choice{curr_mua}(~corr_sig) = 0;
            
        end
        
        % Fluor-wheel/choice
        corr_fluor_wheel = cell(size(event_aligned_df,3),1);
        corr_fluor_choice = cell(size(event_aligned_df,3),1);
        for curr_fluor = 1:size(event_aligned_df,3)
            
            % Wheel
            curr_data1 = zscore(event_aligned_df(:,:,curr_fluor),[],1);
            curr_data2 = zscore(event_aligned_wheel,[],1);
            curr_data2_shuff = curr_data2(shuff_idx);
            
            corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
            
            corr_shuff = ...
                gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                gpuArray(curr_data2_shuff)))./(length(align_times)-1);
            
            corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
            corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
            
            corr_fluor_wheel{curr_fluor} = corr_real;
            corr_fluor_wheel{curr_fluor}(~corr_sig) = 0;
            
            % Choice
            curr_data1 = zscore(event_aligned_df(:,:,curr_fluor),[],1);
            curr_data2 = zscore(choice,[],1);
            curr_data2_shuff = curr_data2(shuff_idx(:,1,:));
            
            corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
            
            corr_shuff = ...
                gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                gpuArray(curr_data2_shuff)))./(length(align_times)-1);
            
            corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
            corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
            
            corr_fluor_choice{curr_fluor} = corr_real;
            corr_fluor_choice{curr_fluor}(~corr_sig) = 0;
            
        end
        
        batch_vars(curr_animal).corr_mua_mua(:,:,curr_day) = corr_mua_mua;
        batch_vars(curr_animal).corr_fluor_fluor(:,:,curr_day) = corr_fluor_fluor;
        batch_vars(curr_animal).corr_fluor_mua(:,:,curr_day) = corr_fluor_mua;
        
        batch_vars(curr_animal).corr_mua_wheel(:,:,curr_day) = corr_mua_wheel;
        batch_vars(curr_animal).corr_mua_choice(:,:,curr_day) = corr_mua_choice;
        
        batch_vars(curr_animal).corr_fluor_wheel(:,:,curr_day) = corr_fluor_wheel;
        batch_vars(curr_animal).corr_fluor_choice(:,:,curr_day) = corr_fluor_choice;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal]);
    
end

fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\corr_mua_fluor_stim';
save(fn,'batch_vars','-v7.3');

%% Batch trial-trial MUA-fluorescence correlation (move-aligned)
% NOTE: I put a bunch of work into doing this on the GPU instead, but the
% memory-limiting steps required FOR loops and the resulting code was
% barely faster

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

protocol = 'vanillaChoiceworld';
% protocol = 'stimSparseNoiseUncorrAsync';
% protocol = 'stimKalatsky';
% protocol = 'AP_choiceWorldStimPassive';

batch_vars = struct;
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    % Skip if this animal doesn't have this experiment
    if isempty(experiments)
        continue
    end
    
    disp(animal);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        
        % Load widefield ROIs
        wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
        load(wf_roi_fn);
        
        % Group striatum depths
        n_depth_groups = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
        [depth_group_n,~,depth_group] = histcounts(spikeDepths,depth_group_edges);
     
        % Define times to align
        n_trials = length(block.paramsValues);
        trial_conditions = signals_events.trialSideValues(1:n_trials).*signals_events.trialContrastValues(1:n_trials);
        trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
        stim_to_move = padarray(wheel_move_time - stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
        stim_to_feedback = padarray(signals_events.responseTimes, ...
            [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
            padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
        
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_move < 0.5 & ...
            stim_to_feedback < 1.5;
        
        align_times = reshape(wheel_move_time(use_trials),[],1);
        
        sample_rate_factor = 3;
        interval_surround = [-0.5,1.5];
        sample_rate = framerate*sample_rate_factor;
        t_surround = interval_surround(1):1/sample_rate:interval_surround(2);
        t_peri_event = bsxfun(@plus,align_times,t_surround);
        
        % Get fluorescence within saved ROIs
        Udf_aligned = AP_align_widefield(animal,day,Udf);
        roi_traces = AP_svd_roi(Udf_aligned,fVdf,[],[],cat(3,wf_roi.mask));
        event_aligned_f = interp1(frame_t,roi_traces',t_peri_event);
        event_aligned_df = interp1(conv(frame_t,[1,1]/2,'valid'),diff(roi_traces,[],2)',t_peri_event);
        event_aligned_df(event_aligned_df < 0) = 0;
        
        % Pull out MUA across depths
        t_peri_event_bins = [t_peri_event - 1/(sample_rate*2), ...
            t_peri_event(:,end) + 1/(sample_rate*2)];
        event_aligned_spikes = nan(size(t_peri_event_bins,1),size(t_peri_event_bins,2)-1,n_depth_groups);
        for curr_depth = 1:n_depth_groups
            use_spikes = spike_times_timeline(depth_group == curr_depth);
            event_aligned_spikes(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(use_spikes,t_peri_event_bins(x,:)),[1:length(align_times)]','uni',false));
        end
        
        % Get wheel movements around all aligned events
        event_aligned_wheel = interp1(conv2(Timeline.rawDAQTimestamps,[0.5,0.5],'valid'), ...
            wheel_velocity,t_peri_event);
        
        % Choice of all events (L/R movement)
        choice = ...
            ((signals_events.trialSideValues(use_trials) == 1 & trial_outcome(use_trials) == -1 | ...
            signals_events.trialSideValues(use_trials) == -1 & trial_outcome(use_trials) == 1) - ...
            (signals_events.trialSideValues(use_trials) == -1 & trial_outcome(use_trials) == -1 | ...
            signals_events.trialSideValues(use_trials) == 1 & trial_outcome(use_trials) == 1))';
        if any(~ismember(choice,[-1,1]))
            error('Non-valid choice')
        end
        
        % Get sig correlations across all time points across modalities
        n_shuff = 1000;
        warning off;
        shuff_idx = ...
            shake(repmat(reshape(1:length(align_times)*length(t_surround), ...
            length(align_times),length(t_surround)),1,1,n_shuff),1);
        warning on ;
        
        % MUA-MUA
        corr_mua_mua = cell(size(event_aligned_spikes,3),size(event_aligned_spikes,3));
        for curr_mua1 = 1:size(event_aligned_spikes,3)
            for curr_mua2 = 1:size(event_aligned_spikes,3)                            
                
                curr_data1 = zscore(event_aligned_spikes(:,:,curr_mua1),[],1);
                curr_data2 = zscore(event_aligned_spikes(:,:,curr_mua2),[],1);               
                curr_data2_shuff = curr_data2(shuff_idx);
                              
                corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
                
                corr_shuff = ...
                    gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                    gpuArray(curr_data2_shuff)))./(length(align_times)-1);
                           
                corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
                corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
                
                corr_mua_mua{curr_mua1,curr_mua2} = corr_real;
                corr_mua_mua{curr_mua1,curr_mua2}(~corr_sig) = 0;
                
            end
        end      
        
        % Fluor-Fluor
        corr_fluor_fluor = cell(size(event_aligned_df,3),size(event_aligned_df,3));
        for curr_fluor1 = 1:size(event_aligned_df,3)
            for curr_fluor2 = 1:size(event_aligned_df,3)
                
                curr_data1 = zscore(event_aligned_df(:,:,curr_fluor1),[],1);
                curr_data2 = zscore(event_aligned_df(:,:,curr_fluor2),[],1);               
                curr_data2_shuff = curr_data2(shuff_idx);
                              
                corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
                
                corr_shuff = ...
                    gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                    gpuArray(curr_data2_shuff)))./(length(align_times)-1);
                
                corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
                corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
                
                corr_fluor_fluor{curr_fluor1,curr_fluor2} = corr_real;
                corr_fluor_fluor{curr_fluor1,curr_fluor2}(~corr_sig) = 0;
                
            end
        end     
        
        % Fluor-MUA
        corr_fluor_mua = cell(size(event_aligned_df,3),size(event_aligned_spikes,3));
        for curr_fluor = 1:size(event_aligned_df,3)
            for curr_mua = 1:size(event_aligned_spikes,3)
                
                curr_data1 = zscore(event_aligned_df(:,:,curr_fluor),[],1);
                curr_data2 = zscore(event_aligned_spikes(:,:,curr_mua),[],1);               
                curr_data2_shuff = curr_data2(shuff_idx);
                              
                corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
                
                corr_shuff = ...
                    gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                    gpuArray(curr_data2_shuff)))./(length(align_times)-1);
                
                corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
                corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
                
                corr_fluor_mua{curr_fluor,curr_mua} = corr_real;
                corr_fluor_mua{curr_fluor,curr_mua}(~corr_sig) = 0;
                
            end
        end      
        
        % MUA-wheel/choice
        corr_mua_wheel = cell(size(event_aligned_spikes,3),1);
        corr_mua_choice = cell(size(event_aligned_spikes,3),1);
        for curr_mua = 1:size(event_aligned_spikes,3)
            
            % Wheel
            curr_data1 = zscore(event_aligned_spikes(:,:,curr_mua),[],1);
            curr_data2 = zscore(event_aligned_wheel,[],1);
            curr_data2_shuff = curr_data2(shuff_idx);
            
            corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
            
            corr_shuff = ...
                gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                gpuArray(curr_data2_shuff)))./(length(align_times)-1);
            
            corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
            corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
            
            corr_mua_wheel{curr_mua} = corr_real;
            corr_mua_wheel{curr_mua}(~corr_sig) = 0;
            
            % Choice
            curr_data1 = zscore(event_aligned_spikes(:,:,curr_mua),[],1);
            curr_data2 = zscore(choice,[],1);
            curr_data2_shuff = curr_data2(shuff_idx(:,1,:));
            
            corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
            
            corr_shuff = ...
                gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                gpuArray(curr_data2_shuff)))./(length(align_times)-1);
            
            corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
            corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
            
            corr_mua_choice{curr_mua} = corr_real;
            corr_mua_choice{curr_mua}(~corr_sig) = 0;
            
        end
        
        % Fluor-wheel/choice
        corr_fluor_wheel = cell(size(event_aligned_df,3),1);
        corr_fluor_choice = cell(size(event_aligned_df,3),1);
        for curr_fluor = 1:size(event_aligned_df,3)
            
            % Wheel
            curr_data1 = zscore(event_aligned_df(:,:,curr_fluor),[],1);
            curr_data2 = zscore(event_aligned_wheel,[],1);
            curr_data2_shuff = curr_data2(shuff_idx);
            
            corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
            
            corr_shuff = ...
                gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                gpuArray(curr_data2_shuff)))./(length(align_times)-1);
            
            corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
            corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
            
            corr_fluor_wheel{curr_fluor} = corr_real;
            corr_fluor_wheel{curr_fluor}(~corr_sig) = 0;
            
            % Choice
            curr_data1 = zscore(event_aligned_df(:,:,curr_fluor),[],1);
            curr_data2 = zscore(choice,[],1);
            curr_data2_shuff = curr_data2(shuff_idx(:,1,:));
            
            corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
            
            corr_shuff = ...
                gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                gpuArray(curr_data2_shuff)))./(length(align_times)-1);
            
            corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
            corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
            
            corr_fluor_choice{curr_fluor} = corr_real;
            corr_fluor_choice{curr_fluor}(~corr_sig) = 0;
            
        end
        
        batch_vars(curr_animal).corr_mua_mua(:,:,curr_day) = corr_mua_mua;
        batch_vars(curr_animal).corr_fluor_fluor(:,:,curr_day) = corr_fluor_fluor;
        batch_vars(curr_animal).corr_fluor_mua(:,:,curr_day) = corr_fluor_mua;
        
        batch_vars(curr_animal).corr_mua_wheel(:,:,curr_day) = corr_mua_wheel;
        batch_vars(curr_animal).corr_mua_choice(:,:,curr_day) = corr_mua_choice;
        
        batch_vars(curr_animal).corr_fluor_wheel(:,:,curr_day) = corr_fluor_wheel;
        batch_vars(curr_animal).corr_fluor_choice(:,:,curr_day) = corr_fluor_choice;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal]);
    
end

fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\corr_mua_fluor_move';
save(fn,'batch_vars','-v7.3');


%% Batch trial-trial correlation (stim-aligned, within-condition shuffle)
% 
% NOTE: I put a bunch of work into doing this on the GPU instead, but the
% memory-limiting steps required FOR loops and the resulting code was
% barely faster

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

protocol = 'vanillaChoiceworld';
% protocol = 'stimSparseNoiseUncorrAsync';
% protocol = 'stimKalatsky';
% protocol = 'AP_choiceWorldStimPassive';

batch_vars = struct;
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    % Skip if this animal doesn't have this experiment
    if isempty(experiments)
        continue
    end
    
    disp(animal);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        
        % Load widefield ROIs
        wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
        load(wf_roi_fn);
        
        % Group striatum depths
        n_depth_groups = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
        [depth_group_n,~,depth_group] = histcounts(spikeDepths,depth_group_edges);
     
        % Define times to align
        n_trials = length(block.paramsValues);
        trial_conditions = signals_events.trialSideValues(1:n_trials).*signals_events.trialContrastValues(1:n_trials);
        trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
        stim_to_move = padarray(wheel_move_time - stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
        stim_to_feedback = padarray(signals_events.responseTimes, ...
            [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
            padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
        
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_move > 0.5 & ...
            stim_to_feedback < 1.5;
        
        align_times = reshape(stimOn_times(use_trials),[],1);
        
        sample_rate_factor = 3;
        interval_surround = [-0.5,1.5];
        sample_rate = framerate*sample_rate_factor;
        t_surround = interval_surround(1):1/sample_rate:interval_surround(2);
        t_peri_event = bsxfun(@plus,align_times,t_surround);
        
        % Get fluorescence within saved ROIs
        Udf_aligned = AP_align_widefield(animal,day,Udf);
        roi_traces = AP_svd_roi(Udf_aligned,fVdf,[],[],cat(3,wf_roi.mask));
        event_aligned_f = interp1(frame_t,roi_traces',t_peri_event);
        event_aligned_df = interp1(conv(frame_t,[1,1]/2,'valid'),diff(roi_traces,[],2)',t_peri_event);
        event_aligned_df(event_aligned_df < 0) = 0;
        
        % Make the R hemisphere ROIs actually be L-R
        event_aligned_df(:,:,size(wf_roi,1)+1:end) = ...
            event_aligned_df(:,:,1:size(wf_roi,1)) - ...
            event_aligned_df(:,:,size(wf_roi,1)+1:end);
        
        % Pull out MUA across depths
        t_peri_event_bins = [t_peri_event - 1/(sample_rate*2), ...
            t_peri_event(:,end) + 1/(sample_rate*2)];
        event_aligned_spikes = nan(size(t_peri_event_bins,1),size(t_peri_event_bins,2)-1,n_depth_groups);
        for curr_depth = 1:n_depth_groups
            use_spikes = spike_times_timeline(depth_group == curr_depth);
            event_aligned_spikes(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(use_spikes,t_peri_event_bins(x,:)),[1:length(align_times)]','uni',false));
        end
        
        % Get wheel movements around all aligned events
        event_aligned_wheel = interp1(conv2(Timeline.rawDAQTimestamps,[0.5,0.5],'valid'), ...
            wheel_velocity,t_peri_event);
        
        % Choice of all events (L/R movement)
        choice = ...
            ((signals_events.trialSideValues(use_trials) == 1 & trial_outcome(use_trials) == -1 | ...
            signals_events.trialSideValues(use_trials) == -1 & trial_outcome(use_trials) == 1) - ...
            (signals_events.trialSideValues(use_trials) == -1 & trial_outcome(use_trials) == -1 | ...
            signals_events.trialSideValues(use_trials) == 1 & trial_outcome(use_trials) == 1))';
        if any(~ismember(choice,[-1,1]))
            error('Non-valid choice')
        end
        
        % Get sig correlations across all time points across modalities
        % (SHUFFLE ONLY WITHIN CONDITION!)
        n_shuff = 1000;
        warning off;
        use_conditions = unique(trial_conditions(use_trials));
        data_idx = reshape(1:length(align_times)*length(t_surround), ...
            length(align_times),length(t_surround));
        shuff_idx = nan(length(align_times),length(t_surround),n_shuff);
        for curr_condition_idx = 1:length(use_conditions)
            curr_condition = use_conditions(curr_condition_idx);
            curr_trials = trial_conditions(use_trials) == curr_condition;
            shuff_idx(curr_trials,:,:) = ...
                shake(repmat(data_idx(curr_trials,:,:),1,1,n_shuff),1);
        end
        warning on ;
        
        % MUA-MUA
        corr_mua_mua = cell(size(event_aligned_spikes,3),size(event_aligned_spikes,3));
        for curr_mua1 = 1:size(event_aligned_spikes,3)
            for curr_mua2 = 1:size(event_aligned_spikes,3)                            
                
                curr_data1 = zscore(event_aligned_spikes(:,:,curr_mua1),[],1);
                curr_data2 = zscore(event_aligned_spikes(:,:,curr_mua2),[],1);               
                curr_data2_shuff = curr_data2(shuff_idx);
                              
                corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
                
                corr_shuff = ...
                    gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                    gpuArray(curr_data2_shuff)))./(length(align_times)-1);
                           
                corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
                corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
                
                corr_mua_mua{curr_mua1,curr_mua2} = corr_real;
                corr_mua_mua{curr_mua1,curr_mua2}(~corr_sig) = 0;
                
            end
        end      
        
        % Fluor-Fluor
        corr_fluor_fluor = cell(size(event_aligned_df,3),size(event_aligned_df,3));
        for curr_fluor1 = 1:size(event_aligned_df,3)
            for curr_fluor2 = 1:size(event_aligned_df,3)
                
                curr_data1 = zscore(event_aligned_df(:,:,curr_fluor1),[],1);
                curr_data2 = zscore(event_aligned_df(:,:,curr_fluor2),[],1);               
                curr_data2_shuff = curr_data2(shuff_idx);
                              
                corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
                
                corr_shuff = ...
                    gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                    gpuArray(curr_data2_shuff)))./(length(align_times)-1);
                
                corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
                corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
                
                corr_fluor_fluor{curr_fluor1,curr_fluor2} = corr_real;
                corr_fluor_fluor{curr_fluor1,curr_fluor2}(~corr_sig) = 0;
                
            end
        end     
        
        % Fluor-MUA
        corr_fluor_mua = cell(size(event_aligned_df,3),size(event_aligned_spikes,3));
        for curr_fluor = 1:size(event_aligned_df,3)
            for curr_mua = 1:size(event_aligned_spikes,3)
                
                curr_data1 = zscore(event_aligned_df(:,:,curr_fluor),[],1);
                curr_data2 = zscore(event_aligned_spikes(:,:,curr_mua),[],1);               
                curr_data2_shuff = curr_data2(shuff_idx);
                              
                corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
                
                corr_shuff = ...
                    gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                    gpuArray(curr_data2_shuff)))./(length(align_times)-1);
                
                corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
                corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
                
                corr_fluor_mua{curr_fluor,curr_mua} = corr_real;
                corr_fluor_mua{curr_fluor,curr_mua}(~corr_sig) = 0;
                
            end
        end      
        
        % MUA-wheel/choice
        corr_mua_wheel = cell(size(event_aligned_spikes,3),1);
        corr_mua_choice = cell(size(event_aligned_spikes,3),1);
        for curr_mua = 1:size(event_aligned_spikes,3)
            
            % Wheel
            curr_data1 = zscore(event_aligned_spikes(:,:,curr_mua),[],1);
            curr_data2 = zscore(event_aligned_wheel,[],1);
            curr_data2_shuff = curr_data2(shuff_idx);
            
            corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
            
            corr_shuff = ...
                gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                gpuArray(curr_data2_shuff)))./(length(align_times)-1);
            
            corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
            corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
            
            corr_mua_wheel{curr_mua} = corr_real;
            corr_mua_wheel{curr_mua}(~corr_sig) = 0;
            
            % Choice
            curr_data1 = zscore(event_aligned_spikes(:,:,curr_mua),[],1);
            curr_data2 = zscore(choice,[],1);
            curr_data2_shuff = curr_data2(shuff_idx(:,1,:));
            
            corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
            
            corr_shuff = ...
                gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                gpuArray(curr_data2_shuff)))./(length(align_times)-1);
            
            corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
            corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
            
            corr_mua_choice{curr_mua} = corr_real;
            corr_mua_choice{curr_mua}(~corr_sig) = 0;
            
        end
        
        % Fluor-wheel/choice
        corr_fluor_wheel = cell(size(event_aligned_df,3),1);
        corr_fluor_choice = cell(size(event_aligned_df,3),1);
        for curr_fluor = 1:size(event_aligned_df,3)
            
            % Wheel
            curr_data1 = zscore(event_aligned_df(:,:,curr_fluor),[],1);
            curr_data2 = zscore(event_aligned_wheel,[],1);
            curr_data2_shuff = curr_data2(shuff_idx);
            
            corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
            
            corr_shuff = ...
                gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                gpuArray(curr_data2_shuff)))./(length(align_times)-1);
            
            corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
            corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
            
            corr_fluor_wheel{curr_fluor} = corr_real;
            corr_fluor_wheel{curr_fluor}(~corr_sig) = 0;
            
            % Choice
            curr_data1 = zscore(event_aligned_df(:,:,curr_fluor),[],1);
            curr_data2 = zscore(choice,[],1);
            curr_data2_shuff = curr_data2(shuff_idx(:,1,:));
            
            corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
            
            corr_shuff = ...
                gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                gpuArray(curr_data2_shuff)))./(length(align_times)-1);
            
            corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
            corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
            
            corr_fluor_choice{curr_fluor} = corr_real;
            corr_fluor_choice{curr_fluor}(~corr_sig) = 0;
            
        end
        
        batch_vars(curr_animal).corr_mua_mua(:,:,curr_day) = corr_mua_mua;
        batch_vars(curr_animal).corr_fluor_fluor(:,:,curr_day) = corr_fluor_fluor;
        batch_vars(curr_animal).corr_fluor_mua(:,:,curr_day) = corr_fluor_mua;
        
        batch_vars(curr_animal).corr_mua_wheel(:,:,curr_day) = corr_mua_wheel;
        batch_vars(curr_animal).corr_mua_choice(:,:,curr_day) = corr_mua_choice;
        
        batch_vars(curr_animal).corr_fluor_wheel(:,:,curr_day) = corr_fluor_wheel;
        batch_vars(curr_animal).corr_fluor_choice(:,:,curr_day) = corr_fluor_choice;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal]);
    
end

fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\corr_mua_fluor_stim_latemove_conditionshuff';
save(fn,'batch_vars','-v7.3');


%% Batch trial-trial correlation (move-aligned, within-condition shuffle)
% 
% NOTE: I put a bunch of work into doing this on the GPU instead, but the
% memory-limiting steps required FOR loops and the resulting code was
% barely faster

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

protocol = 'vanillaChoiceworld';
% protocol = 'stimSparseNoiseUncorrAsync';
% protocol = 'stimKalatsky';
% protocol = 'AP_choiceWorldStimPassive';

batch_vars = struct;
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    % Skip if this animal doesn't have this experiment
    if isempty(experiments)
        continue
    end
    
    disp(animal);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        
        % Load widefield ROIs
        wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
        load(wf_roi_fn);
        
        % Group striatum depths
        n_depth_groups = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depth_groups+1));
        [depth_group_n,~,depth_group] = histcounts(spikeDepths,depth_group_edges);
     
        % Define times to align
        n_trials = length(block.paramsValues);
        trial_conditions = signals_events.trialSideValues(1:n_trials).*signals_events.trialContrastValues(1:n_trials);
        trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
        stim_to_move = padarray(wheel_move_time - stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
        stim_to_feedback = padarray(signals_events.responseTimes, ...
            [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
            padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
        
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_move > 0.5 & ...
            stim_to_feedback < 1.5;
        
        align_times = reshape(wheel_move_time(use_trials),[],1);
        
        sample_rate_factor = 3;
        interval_surround = [-0.5,1.5];
        sample_rate = framerate*sample_rate_factor;
        t_surround = interval_surround(1):1/sample_rate:interval_surround(2);
        t_peri_event = bsxfun(@plus,align_times,t_surround);
        
        % Get fluorescence within saved ROIs
        Udf_aligned = AP_align_widefield(animal,day,Udf);
        roi_traces = AP_svd_roi(Udf_aligned,fVdf,[],[],cat(3,wf_roi.mask));
        event_aligned_f = interp1(frame_t,roi_traces',t_peri_event);
        event_aligned_df = interp1(conv(frame_t,[1,1]/2,'valid'),diff(roi_traces,[],2)',t_peri_event);
        event_aligned_df(event_aligned_df < 0) = 0;
        
        % Make the R hemisphere ROIs actually be L-R
        event_aligned_df(:,:,size(wf_roi,1)+1:end) = ...
            event_aligned_df(:,:,1:size(wf_roi,1)) - ...
            event_aligned_df(:,:,size(wf_roi,1)+1:end);
        
        % Pull out MUA across depths
        t_peri_event_bins = [t_peri_event - 1/(sample_rate*2), ...
            t_peri_event(:,end) + 1/(sample_rate*2)];
        event_aligned_spikes = nan(size(t_peri_event_bins,1),size(t_peri_event_bins,2)-1,n_depth_groups);
        for curr_depth = 1:n_depth_groups
            use_spikes = spike_times_timeline(depth_group == curr_depth);
            event_aligned_spikes(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(use_spikes,t_peri_event_bins(x,:)),[1:length(align_times)]','uni',false));
        end
        
        % Get wheel movements around all aligned events
        event_aligned_wheel = interp1(conv2(Timeline.rawDAQTimestamps,[0.5,0.5],'valid'), ...
            wheel_velocity,t_peri_event);
        
        % Choice of all events (L/R movement)
        choice = ...
            ((signals_events.trialSideValues(use_trials) == 1 & trial_outcome(use_trials) == -1 | ...
            signals_events.trialSideValues(use_trials) == -1 & trial_outcome(use_trials) == 1) - ...
            (signals_events.trialSideValues(use_trials) == -1 & trial_outcome(use_trials) == -1 | ...
            signals_events.trialSideValues(use_trials) == 1 & trial_outcome(use_trials) == 1))';
        if any(~ismember(choice,[-1,1]))
            error('Non-valid choice')
        end
        
        % Get sig correlations across all time points across modalities
        % (SHUFFLE ONLY WITHIN CONDITION!)
        n_shuff = 1000;
        warning off;
        use_conditions = unique(trial_conditions(use_trials));
        data_idx = reshape(1:length(align_times)*length(t_surround), ...
            length(align_times),length(t_surround));
        shuff_idx = nan(length(align_times),length(t_surround),n_shuff);
        for curr_condition_idx = 1:length(use_conditions)
            curr_condition = use_conditions(curr_condition_idx);
            curr_trials = trial_conditions(use_trials) == curr_condition;
            shuff_idx(curr_trials,:,:) = ...
                shake(repmat(data_idx(curr_trials,:,:),1,1,n_shuff),1);
        end
        warning on ;
        
        % MUA-MUA
        corr_mua_mua = cell(size(event_aligned_spikes,3),size(event_aligned_spikes,3));
        for curr_mua1 = 1:size(event_aligned_spikes,3)
            for curr_mua2 = 1:size(event_aligned_spikes,3)                            
                
                curr_data1 = zscore(event_aligned_spikes(:,:,curr_mua1),[],1);
                curr_data2 = zscore(event_aligned_spikes(:,:,curr_mua2),[],1);               
                curr_data2_shuff = curr_data2(shuff_idx);
                              
                corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
                
                corr_shuff = ...
                    gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                    gpuArray(curr_data2_shuff)))./(length(align_times)-1);
                           
                corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
                corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
                
                corr_mua_mua{curr_mua1,curr_mua2} = corr_real;
                corr_mua_mua{curr_mua1,curr_mua2}(~corr_sig) = 0;
                
            end
        end      
        
        % Fluor-Fluor
        corr_fluor_fluor = cell(size(event_aligned_df,3),size(event_aligned_df,3));
        for curr_fluor1 = 1:size(event_aligned_df,3)
            for curr_fluor2 = 1:size(event_aligned_df,3)
                
                curr_data1 = zscore(event_aligned_df(:,:,curr_fluor1),[],1);
                curr_data2 = zscore(event_aligned_df(:,:,curr_fluor2),[],1);               
                curr_data2_shuff = curr_data2(shuff_idx);
                              
                corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
                
                corr_shuff = ...
                    gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                    gpuArray(curr_data2_shuff)))./(length(align_times)-1);
                
                corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
                corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
                
                corr_fluor_fluor{curr_fluor1,curr_fluor2} = corr_real;
                corr_fluor_fluor{curr_fluor1,curr_fluor2}(~corr_sig) = 0;
                
            end
        end     
        
        % Fluor-MUA
        corr_fluor_mua = cell(size(event_aligned_df,3),size(event_aligned_spikes,3));
        for curr_fluor = 1:size(event_aligned_df,3)
            for curr_mua = 1:size(event_aligned_spikes,3)
                
                curr_data1 = zscore(event_aligned_df(:,:,curr_fluor),[],1);
                curr_data2 = zscore(event_aligned_spikes(:,:,curr_mua),[],1);               
                curr_data2_shuff = curr_data2(shuff_idx);
                              
                corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
                
                corr_shuff = ...
                    gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                    gpuArray(curr_data2_shuff)))./(length(align_times)-1);
                
                corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
                corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
                
                corr_fluor_mua{curr_fluor,curr_mua} = corr_real;
                corr_fluor_mua{curr_fluor,curr_mua}(~corr_sig) = 0;
                
            end
        end      
        
        % MUA-wheel/choice
        corr_mua_wheel = cell(size(event_aligned_spikes,3),1);
        corr_mua_choice = cell(size(event_aligned_spikes,3),1);
        for curr_mua = 1:size(event_aligned_spikes,3)
            
            % Wheel
            curr_data1 = zscore(event_aligned_spikes(:,:,curr_mua),[],1);
            curr_data2 = zscore(event_aligned_wheel,[],1);
            curr_data2_shuff = curr_data2(shuff_idx);
            
            corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
            
            corr_shuff = ...
                gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                gpuArray(curr_data2_shuff)))./(length(align_times)-1);
            
            corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
            corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
            
            corr_mua_wheel{curr_mua} = corr_real;
            corr_mua_wheel{curr_mua}(~corr_sig) = 0;
            
            % Choice
            curr_data1 = zscore(event_aligned_spikes(:,:,curr_mua),[],1);
            curr_data2 = zscore(choice,[],1);
            curr_data2_shuff = curr_data2(shuff_idx(:,1,:));
            
            corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
            
            corr_shuff = ...
                gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                gpuArray(curr_data2_shuff)))./(length(align_times)-1);
            
            corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
            corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
            
            corr_mua_choice{curr_mua} = corr_real;
            corr_mua_choice{curr_mua}(~corr_sig) = 0;
            
        end
        
        % Fluor-wheel/choice
        corr_fluor_wheel = cell(size(event_aligned_df,3),1);
        corr_fluor_choice = cell(size(event_aligned_df,3),1);
        for curr_fluor = 1:size(event_aligned_df,3)
            
            % Wheel
            curr_data1 = zscore(event_aligned_df(:,:,curr_fluor),[],1);
            curr_data2 = zscore(event_aligned_wheel,[],1);
            curr_data2_shuff = curr_data2(shuff_idx);
            
            corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
            
            corr_shuff = ...
                gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                gpuArray(curr_data2_shuff)))./(length(align_times)-1);
            
            corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
            corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
            
            corr_fluor_wheel{curr_fluor} = corr_real;
            corr_fluor_wheel{curr_fluor}(~corr_sig) = 0;
            
            % Choice
            curr_data1 = zscore(event_aligned_df(:,:,curr_fluor),[],1);
            curr_data2 = zscore(choice,[],1);
            curr_data2_shuff = curr_data2(shuff_idx(:,1,:));
            
            corr_real = (curr_data1'*curr_data2)./(length(align_times)-1);
            
            corr_shuff = ...
                gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                gpuArray(curr_data2_shuff)))./(length(align_times)-1);
            
            corr_shuff_cutoff = prctile(corr_shuff,[2.5,97.5],3);
            corr_sig = (corr_real > corr_shuff_cutoff(:,:,2)) | (corr_real < corr_shuff_cutoff(:,:,1));
            
            corr_fluor_choice{curr_fluor} = corr_real;
            corr_fluor_choice{curr_fluor}(~corr_sig) = 0;
            
        end
        
        batch_vars(curr_animal).corr_mua_mua(:,:,curr_day) = corr_mua_mua;
        batch_vars(curr_animal).corr_fluor_fluor(:,:,curr_day) = corr_fluor_fluor;
        batch_vars(curr_animal).corr_fluor_mua(:,:,curr_day) = corr_fluor_mua;
        
        batch_vars(curr_animal).corr_mua_wheel(:,:,curr_day) = corr_mua_wheel;
        batch_vars(curr_animal).corr_mua_choice(:,:,curr_day) = corr_mua_choice;
        
        batch_vars(curr_animal).corr_fluor_wheel(:,:,curr_day) = corr_fluor_wheel;
        batch_vars(curr_animal).corr_fluor_choice(:,:,curr_day) = corr_fluor_choice;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal]);
    
end

fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\corr_mua_fluor_move_latemove_conditionshuff';
save(fn,'batch_vars','-v7.3');






