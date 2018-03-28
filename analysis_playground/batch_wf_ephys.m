%% Create widefield ROIs

alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment';
load([alignment_path filesep 'animal_wf_tform']);
im_size = animal_wf_tform(1).im_size;

roi_areas = {'V1','AM','RSPp','PPC','FRm','FRl','SMl','SMf'};

wf_roi = struct('area',cell(length(roi_areas),2),'mask',cell(length(roi_areas),2));

% Draw all ROIs on left hemisphere
for curr_area = 1:length(roi_areas);
    curr_roi_name = [roi_areas{curr_area} '_L'];
    disp(curr_roi_name);
    wf_roi(curr_area).area = curr_roi_name;
    [~,wf_roi(curr_area).mask] = AP_svd_roi(nan(im_size),[],'master');
end

% Reflect all ROIs to right hemisphere
for curr_area = 1:length(roi_areas);
    curr_roi_name = [roi_areas{curr_area} '_R'];
    wf_roi(curr_area,2).area = curr_roi_name;
    
    L_roi = wf_roi(curr_area,1).mask;
    R_roi = AP_reflect_widefield(L_roi) > 0;
    
    wf_roi(curr_area,2).mask = R_roi;
end

wf_roi_plot = sum(cat(3,wf_roi(:).mask),3);
figure;imagesc(wf_roi_plot);
colormap(flipud(gray));
AP_reference_outline('ccf_aligned','r');AP_reference_outline('retinotopy','b');
axis image off;
title('New ROIs')

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

%% Batch widefield choiceworld (GENERAL - UPSAMPLED)

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';

for trialtype_align = {'stim','move'};
    trialtype_align = cell2mat(trialtype_align);
    for trialtype_timing = {'earlymove','latemove'};
        trialtype_timing = cell2mat(trialtype_timing);
        for trialtype_success = {'hit','miss'};
            trialtype_success = cell2mat(trialtype_success);
            
            for curr_animal = 1:length(animals)
                
                animal = animals{curr_animal};
                experiments = AP_find_experiments(animal,protocol);
                
                disp(animal);
                
                experiments = experiments([experiments.imaging] & [experiments.ephys]);
                
                load_parts.cam = false;
                load_parts.imaging = true;
                load_parts.ephys = false;
                
                batch_vars = struct;
                
                batch_vars.im_aligned = [];
                batch_vars.n_im_aligned = [];
                
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
                    
                    switch trialtype_align
                        case 'stim'
                            use_align = stimOn_times;
                        case 'move'
                            use_align = wheel_move_time;
                    end
                    
                    switch trialtype_timing
                        case 'earlymove'
                            use_stim_to_move = stim_to_move < 0.5;
                            use_stim_to_feedback = stim_to_feedback < 1.5;
                        case 'latemove'
                            use_stim_to_move = stim_to_move > 0.5;
                            use_stim_to_feedback = stim_to_feedback < 1.5;
                        case 'nomove'
                            use_stim_to_move = stim_to_move > 0.5;
                            use_stim_to_feedback = stim_to_feedback >= 1.5;
                    end
                    
                    switch trialtype_success
                        case 'hit'
                            use_outcome = trial_outcome == 1;
                        case 'miss'
                            use_outcome = trial_outcome == -1;
                    end
                    
                    use_trials = ...
                        trial_outcome ~= 0 & ...
                        ~signals_events.repeatTrialValues(1:n_trials);
                    
                    % Set options
                    surround_window = [-0.5,2];
                    upsample_factor = 3;
                    
                    framerate = 1./median(diff(frame_t));
                    surround_samplerate = 1/(framerate*upsample_factor);
                    t_surround = surround_window(1):surround_samplerate:surround_window(2);
                    
                    t_baseline = [-0.5,0]; % (from stim onset)
                    t_baseline_surround = t_baseline(1):surround_samplerate:t_baseline(2);
                    
                    % Average (time course) responses
                    im = nan(size(U,1),size(U,2),length(t_surround),length(conditions));
                    
                    for curr_condition_idx = 1:length(conditions)
                        curr_condition = conditions(curr_condition_idx);
                        
                        curr_trials = use_trials & trial_conditions == curr_condition & ...
                            use_outcome & use_stim_to_move & use_stim_to_feedback;
                        curr_align = use_align(curr_trials);
                        curr_align_baseline = stimOn_times(curr_trials);
                        
                        if length(curr_align) > 5
                            curr_surround_times = bsxfun(@plus, curr_align(:), t_surround);
                            curr_surround_baseline = bsxfun(@plus, curr_align_baseline(:), t_baseline_surround);
                            
                            peri_stim_v = permute(mean(interp1(frame_t,fVdf',curr_surround_times),1),[3,2,1]);
                            peri_stim_v_baseline = permute(mean(interp1(frame_t,fVdf',curr_surround_baseline),1),[3,2,1]);
                            
                            peri_stim_v_baselinesub = bsxfun(@minus, ...
                                peri_stim_v,nanmean(peri_stim_v_baseline,2));
                            
                            im(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf,peri_stim_v);
                        end
                        
                    end
                    
                    % Align and average as it goes, otherwise the variable's too big
                    im_aligned = AP_align_widefield(animal,day,im);
                    batch_vars.im_aligned = nansum(cat(5,batch_vars.im_aligned,im_aligned),5);
                    
                    % Count conditions to divide at end
                    batch_vars.n_im_aligned = ...
                        sum(cat(5,batch_vars.n_im_aligned,any(any(any(im_aligned,1),2),3)),5);
                    
                    % Prepare for next loop
                    AP_print_progress_fraction(curr_day,length(experiments))
                    clearvars -except animals protocol trialtype_align ...
                        trialtype_timing trialtype_success curr_animal  ...
                        experiments curr_day animal batch_vars load_parts
                    
                end
                
                % Divide sum to get average
                im_aligned_avg = bsxfun(@rdivide,batch_vars.im_aligned,batch_vars.n_im_aligned);
                
                % Save
                save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
                save([save_path filesep animal '_im_' trialtype_align '_' trialtype_timing '_' trialtype_success],'im_aligned_avg','-v7.3');
                
                disp(['Finished ' animal]);
                
            end
            
            disp('Finished batch.')
            warning('This uses -v7.3 and therefore compresses data, switch to dat in the future');
            
        end
    end
end

%% Batch widefield choiceworld (TRIAL ID - ROIS ONLY, also regression)

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';
batch_vars = struct;

% Load widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    disp(animal);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = false;
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        AP_load_experiment
        
        % Get ROI traces
        aUdf = single(AP_align_widefield(animal,day,Udf));
        roi_traces = AP_svd_roi(aUdf,fVdf,[],[],cat(3,wf_roi.mask));
        
        % Make DDF
        roi_traces_derivative = diff(roi_traces,[],2);
        roi_traces_derivative(roi_traces_derivative < 0) = 0;
        
        % Define trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Get event-aligned fluorescence
        raster_window = [-0.5,3];
        upsample_factor = 5;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        baseline_window = [-0.2,0]; % (from stim onset)
        t_baseline = baseline_window(1):raster_sample_rate:baseline_window(2);
        
        % [condition, time, roi, alignment, timing, ddf/df]
        roi_psth = nan(n_conditions,length(t),numel(wf_roi),2,2);
        for curr_align = 1:2
            switch curr_align
                case 1
                    use_align = stimOn_times;
                case 2
                    use_align = wheel_move_time';
                    use_align(isnan(use_align)) = 0;
            end
            
            t_peri_event = bsxfun(@plus,use_align,t);
                        
            % DDF/F
            event_aligned_ddf = interp1(conv2(frame_t,[1,1]/2,'valid'),roi_traces_derivative',t_peri_event);           
            
            % DF/F (subtract baseline)
            use_align_baseline = stimOn_times;
            t_peri_event_baseline = bsxfun(@plus,use_align_baseline,t_baseline);
            event_aligned_df_raw = interp1(frame_t,roi_traces',t_peri_event);
            event_aligned_df_baseline = interp1(frame_t,roi_traces',t_peri_event_baseline);
            event_aligned_df = bsxfun(@minus,event_aligned_df_raw,nanmean(event_aligned_df_baseline,2));
                       
            for curr_roi = 1:numel(wf_roi)
                
                % DDF/F
                [curr_ids,curr_mean_psth_ddf] = ...
                    grpstats(event_aligned_ddf(use_trials,:,curr_roi), ...
                    trial_id(use_trials),{'gname',@(x) mean(x,1)});
                curr_ids = cellfun(@str2num,curr_ids);
                
                % DF/F
                [curr_ids,curr_mean_psth_df] = ...
                    grpstats(event_aligned_df(use_trials,:,curr_roi), ...
                    trial_id(use_trials),{'gname',@(x) mean(x,1)});
                curr_ids = cellfun(@str2num,curr_ids);
                
                roi_psth(curr_ids,:,curr_roi,curr_align,1) = curr_mean_psth_ddf;
                roi_psth(curr_ids,:,curr_roi,curr_align,2) = curr_mean_psth_df;
                
            end
        end
        
        % Count trials per condition
        condition_counts = histcounts(trial_id(use_trials), ...
            'BinLimits',[1,n_conditions],'BinMethod','integers')';
          
        % Regress out contrast/choice from both real and predicted
        
        % (make R ROIs into L-R)
        roi_traces(size(wf_roi,1)+1:end,:) = ...
            roi_traces(1:size(wf_roi,1),:) - roi_traces(size(wf_roi,1)+1:end,:);
        
        roi_traces_derivative(size(wf_roi,1)+1:end,:) = ...
            roi_traces_derivative(1:size(wf_roi,1),:) - roi_traces_derivative(size(wf_roi,1)+1:end,:);
       
        contrast_exp = 1;
        
        % [roi,time,align,timing,param,ddf/df]
        n_rois = numel(wf_roi);
        activity_model_params = nan(n_rois,length(t),2,2,4,2);
        for curr_timing = 1:2
            for curr_align = 1:2
                switch curr_align
                    case 1
                        use_align = stimOn_times;
                    case 2
                        use_align = wheel_move_time';
                        use_align(isnan(use_align)) = 0;
                end
                
                use_timing_trials = use_trials' & trial_conditions(:,4) == curr_timing;
                               
                t_peri_event = bsxfun(@plus,use_align(use_timing_trials),t);
                
                % DDF/F
                curr_aligned_ddf = interp1(conv2(frame_t,[1,1]/2,'valid'),roi_traces_derivative',t_peri_event);
                
                % DF/F (subtract baseline)
                use_align_baseline = stimOn_times;
                t_peri_event_baseline = bsxfun(@plus,use_align_baseline(use_timing_trials),t_baseline);
                curr_aligned_df_raw = interp1(frame_t,roi_traces',t_peri_event);
                curr_aligned_df_baseline = interp1(frame_t,roi_traces',t_peri_event_baseline);
                curr_aligned_df = bsxfun(@minus,curr_aligned_df_raw,nanmean(curr_aligned_df_baseline,2));                
                
                for curr_roi = 1:n_rois
            
                    curr_contrast_sides = zeros(sum(use_timing_trials),2);
                    curr_contrast_sides(trial_conditions(use_timing_trials,2) == -1,1) = ...
                       trial_conditions(trial_conditions(use_timing_trials,2) == -1,1);
                    curr_contrast_sides(trial_conditions(use_timing_trials,2) == 1,2) = ...
                       trial_conditions(trial_conditions(use_timing_trials,2) == 1,1);
                    
                    curr_choices = trial_conditions(use_timing_trials,3);                                       
                    
                    valid_trials = all(~isnan([curr_aligned_ddf(:,:,curr_roi)]),2);
                    
                    for curr_t = 1:length(t)            
                        % DDF/F
                        curr_act = curr_aligned_ddf(:,curr_t,curr_roi);
                        params_fit = [ones(sum(valid_trials),1), ...
                            curr_contrast_sides(valid_trials,:).^contrast_exp, ...
                            curr_choices(valid_trials)]\curr_act(valid_trials);
                        activity_model_params(curr_roi,curr_t,curr_align,curr_timing,:,1) = ...
                            params_fit;
                        
                        % DF/F
                        curr_act = curr_aligned_df(:,curr_t,curr_roi);
                        params_fit = [ones(sum(valid_trials),1), ...
                            curr_contrast_sides(valid_trials,:).^contrast_exp, ...
                            curr_choices(valid_trials)]\curr_act(valid_trials);
                        activity_model_params(curr_roi,curr_t,curr_align,curr_timing,:,2) = ...
                            params_fit;
                        
                    end
                end
            end
        end     
        
        % Store
        batch_vars(curr_animal).roi_psth(:,:,:,:,:,curr_day) = roi_psth;
        batch_vars(curr_animal).condition_counts(:,curr_day) = condition_counts;
        batch_vars(curr_animal).activity_model_params(:,:,:,:,:,:,curr_day) = activity_model_params;
        
        % Prep for next loop
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts wf_roi
        
    end
end

% Save
save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
save([save_path filesep 'roi_choiceworld'],'batch_vars');

disp(['Finished batch'])

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
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        [depth_group_n,~,depth_group] = histcounts(spikeDepths,depth_group_edges);
          
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths           
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
        r_px_com = sum(bsxfun(@times,r_px_max_norm,permute(1:n_depths,[1,3,2])),3)./sum(r_px_max_norm,3);
        
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
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        
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
        nlin_params = nan(n_depths,2);
        predicted_spikes_nlin = nan(size(predicted_spikes_reranged));
        for curr_depth = 1:n_depths
            
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
        mua_stim_hit = nan(n_depths,length(t_surround),length(conditions));
        mua_stim_hit_pred = nan(n_depths,length(t_surround),length(conditions));
        
        mua_stim_miss = nan(n_depths,length(t_surround),length(conditions));
        mua_stim_miss_pred = nan(n_depths,length(t_surround),length(conditions));
        
        for curr_depth = 1:n_depths
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
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        
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
        mua_stim_earlymove_hit = nan(n_depths,length(t_surround),length(conditions));
        mua_stim_earlymove_hit_pred = nan(n_depths,length(t_surround),length(conditions));
        
        mua_stim_latemove_hit = nan(n_depths,length(t_surround),length(conditions));
        mua_stim_latemove_hit_pred = nan(n_depths,length(t_surround),length(conditions));
        
        mua_stim_earlymove_miss = nan(n_depths,length(t_surround),length(conditions));
        mua_stim_earlymove_miss_pred = nan(n_depths,length(t_surround),length(conditions));
        
        mua_stim_latemove_miss = nan(n_depths,length(t_surround),length(conditions));
        mua_stim_latemove_miss_pred = nan(n_depths,length(t_surround),length(conditions));
        
        for curr_depth = 1:n_depths
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

%% Batch cortex > striatum prediction by condition and depth (UPDATED)
% still uses raw output, but now estimates lambda and uses trial ID
% also regression from task parameters
% NOTE: ran once with -0.3:0 kernel, now doing 0:0.3 (ctx->str or vv)

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
        
        upsample_factor = 2;
        sample_rate = (1/median(diff(frame_t)))*upsample_factor;
        
        skip_seconds = 60;
        time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        use_svs = 1:100;
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
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        [depth_group_n,~,depth_group] = histcounts(spikeDepths,depth_group_edges);
          
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths           
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);           
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);            
        end
        
%         kernel_t = [-0.3,0]; % TO PREDICT CORTEX -> STRIATUM
        kernel_t = [0,0.3]; % TO PREDICT STRIATUM -> CORTEX
        kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
        zs = [false,true];
        cvfold = 5;      

        % Predict striatal spikes from cortical fluorescence (dfV)
        [k,predicted_spikes,explained_var] = ...
            AP_regresskernel(dfVdf_resample, ...
            binned_spikes,kernel_frames,lambda,zs,cvfold);          
        
        binned_spikes_std = std(binned_spikes,[],2);
        binned_spikes_mean = mean(binned_spikes,2);
        predicted_spikes_reranged = bsxfun(@plus,bsxfun(@times,predicted_spikes, ...
            binned_spikes_std),binned_spikes_mean);     
                        
        % Set times for PSTH
        interval_surround = [-0.5,3];
        t = interval_surround(1):1/sample_rate:interval_surround(2);
        
        % Get trial properties
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Loop through all trial conditions, PSTH        
        depth_psth = nan(n_conditions,length(t),n_depths,2);
        predicted_depth_psth = nan(n_conditions,length(t),n_depths,2);
        for curr_align = 1:2
            switch curr_align
                case 1
                    use_align = stimOn_times;
                case 2
                    use_align = wheel_move_time';
                    use_align(isnan(use_align)) = 0;
            end
            t_peri_event = bsxfun(@plus,use_align,t);
            for curr_depth = 1:n_depths
                
                curr_spikes_binned = interp1(time_bin_centers,binned_spikes(curr_depth,:),t_peri_event)*sample_rate;
                curr_spikes_binned_pred = interp1(time_bin_centers,predicted_spikes_reranged(curr_depth,:),t_peri_event)*sample_rate;
      
                [curr_ids_str,curr_mean_psth] = ...
                    grpstats(curr_spikes_binned(use_trials,:), ...
                    trial_id(use_trials),{'gname','mean'});
                curr_ids = cellfun(@str2num,curr_ids_str);
                
                [curr_ids_str,curr_mean_pred_psth] = ...
                    grpstats(curr_spikes_binned_pred(use_trials,:), ...
                    trial_id(use_trials),{'gname','mean'});
                curr_ids = cellfun(@str2num,curr_ids_str);
                
                depth_psth(curr_ids,:,curr_depth,curr_align) = curr_mean_psth;
                predicted_depth_psth(curr_ids,:,curr_depth,curr_align) = curr_mean_pred_psth;
            end
        end
        
        % Count trials per condition
        condition_counts = histcounts(trial_id(use_trials), ...
            'BinLimits',[1,n_conditions],'BinMethod','integers')';
        
        % Regress out contrast/choice from both real and predicted
        contrast_exp = 0.3;
        
        % [depth,time,align,timing,param]
        activity_model_params = nan(n_depths,length(t),2,2,4);
        predicted_activity_model_params = nan(n_depths,length(t),2,2,4);
        for curr_timing = 1:2
            for curr_align = 1:2
                switch curr_align
                    case 1
                        use_align = stimOn_times;
                    case 2
                        use_align = wheel_move_time';
                        use_align(isnan(use_align)) = 0;
                end
                
                use_timing_trials = use_trials' & trial_conditions(:,4) == curr_timing;
                
                t_peri_event = bsxfun(@plus,use_align(use_timing_trials),t);
                for curr_depth = 1:n_depths
                    
                    curr_spikes_binned = interp1(time_bin_centers,binned_spikes(curr_depth,:),t_peri_event)*sample_rate;
                    curr_spikes_binned_pred = interp1(time_bin_centers,predicted_spikes_reranged(curr_depth,:),t_peri_event)*sample_rate;                    
                    
                    curr_contrast_sides = zeros(sum(use_timing_trials),2);
                    curr_contrast_sides(trial_conditions(use_timing_trials,2) == -1,1) = ...
                       trial_conditions(trial_conditions(use_timing_trials,2) == -1,1);
                    curr_contrast_sides(trial_conditions(use_timing_trials,2) == 1,2) = ...
                       trial_conditions(trial_conditions(use_timing_trials,2) == 1,1);
                    
                    curr_choices = trial_conditions(use_timing_trials,3);                                       
                    
                    valid_trials = all(~isnan([curr_spikes_binned,curr_spikes_binned_pred]),2);
                    
                    for curr_t = 1:length(t)                       
                        curr_act = curr_spikes_binned(:,curr_t);
                        params_fit = [ones(sum(valid_trials),1), ...
                            curr_contrast_sides(valid_trials,:).^contrast_exp, ...
                            curr_choices(valid_trials)]\curr_act(valid_trials);
                        activity_model_params(curr_depth,curr_t,curr_align,curr_timing,:) = ...
                            params_fit;
                        
                        curr_act = curr_spikes_binned_pred(:,curr_t);
                        params_fit = [ones(sum(valid_trials),1), ...
                            curr_contrast_sides(valid_trials,:).^contrast_exp, ...
                            curr_choices(valid_trials)]\curr_act(valid_trials);
                        predicted_activity_model_params(curr_depth,curr_t,curr_align,curr_timing,:) = ...
                            params_fit;
                    end
                end
            end
        end
              
        % Store
        batch_vars(curr_animal).depth_psth(:,:,:,:,curr_day) = depth_psth; 
        batch_vars(curr_animal).predicted_depth_psth(:,:,:,:,curr_day) = predicted_depth_psth; 
        batch_vars(curr_animal).condition_counts(:,curr_day) = condition_counts;
        
        batch_vars(curr_animal).activity_model_params(:,:,:,:,:,curr_day) = activity_model_params;
        batch_vars(curr_animal).predicted_activity_model_params(:,:,:,:,:,curr_day) = predicted_activity_model_params;
        
        % Prepare for next loop       
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
     
    disp(['Finished ' animal])
    
end

save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
save([save_path filesep 'mua_choiceworld_predicted_str-ctx'],'batch_vars');

disp('Finished batch');

%% Batch striatum responses to choiceworld (OLD)

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
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
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
        
        for curr_depth = 1:n_depths
            
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

%% Batch striatum responses to choiceworld (OLD - TRIAL TYPE BY STRUCT)

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
        n_conditions = length(conditions);
        
        % Group multiunit by depth
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
        
        depth_group = discretize(spikeDepths,depth_group_edges);
        
        % Set times for PSTH
        raster_window = [-0.5,3];
        psth_bin_size = 0.001;
        t = raster_window(1):psth_bin_size:raster_window(2);
        
        % Get trial properties 
        n_trials = length(block.paramsValues);
        trial_conditions = signals_events.trialSideValues(1:n_trials).*signals_events.trialContrastValues(1:n_trials);
        trial_outcome = signals_events.hitValues(1:n_trials)-signals_events.missValues(1:n_trials);
        stim_to_move = padarray(wheel_move_time - stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');
        stim_to_feedback = padarray(signals_events.responseTimes, ...
            [0,n_trials-length(signals_events.responseTimes)],NaN,'post') - ...
            padarray(stimOn_times',[0,n_trials-length(stimOn_times)],NaN,'post');    
        
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials);
        
        % Loop through all trial conditions, PSTH
        mua = struct;
        
        for trialtype_align = {'stim','move'};
            trialtype_align = cell2mat(trialtype_align);
            for trialtype_timing = {'earlymove','latemove','nomove'};
                trialtype_timing = cell2mat(trialtype_timing);
                for trialtype_success = {'hit','miss'};
                    trialtype_success = cell2mat(trialtype_success);
                    
                    trialtype = [trialtype_align '_' trialtype_timing '_' trialtype_success];
                    
                    switch trialtype_align
                        case 'stim'
                            use_align = stimOn_times;
                        case 'move'
                            use_align = wheel_move_time;
                    end
                    
                    switch trialtype_timing
                        case 'earlymove'
                            use_stim_to_move = stim_to_move < 0.5;
                            use_stim_to_feedback = stim_to_feedback < 1.5;
                        case 'latemove'
                            use_stim_to_move = stim_to_move > 0.5;
                            use_stim_to_feedback = stim_to_feedback < 1.5;
                        case 'nomove'
                            use_stim_to_move = stim_to_move > 0.5;
                            use_stim_to_feedback = stim_to_feedback >= 1.5;
                    end
                    
                    switch trialtype_success
                        case 'hit'
                            use_outcome = trial_outcome == 1;
                        case 'miss'
                            use_outcome = trial_outcome == -1;
                    end
                    
                    
                    mua.(trialtype) = nan(n_depths,length(t)-1,n_conditions);
                    for curr_depth = 1:n_depths                       
                        for curr_condition_idx = 1:length(conditions)
                            
                            curr_condition = conditions(curr_condition_idx);                            
                            curr_trials = use_trials & trial_conditions == curr_condition & ...
                                use_outcome & use_stim_to_move & use_stim_to_feedback;
                            curr_align = reshape(use_align(curr_trials),[],1);
                            
                            if length(curr_align) > 5
                                t_peri_event = bsxfun(@plus,curr_align,t);
                                use_spikes = spike_times_timeline(depth_group == curr_depth);
                                psth = cell2mat(arrayfun(@(x) ...
                                    histcounts(use_spikes,t_peri_event(x,:)), ...
                                    [1:size(t_peri_event,1)]','uni',false));
                                
                                mua.(trialtype)(curr_depth,:,curr_condition_idx) = ...
                                    nanmean(psth./psth_bin_size,1);                            
                            end                           
                        end                       
                    end                                      
                end
            end
        end    
        
        batch_vars(curr_animal).mua(curr_day) = mua;      
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal])
    
end

% Save
save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
save([save_path filesep 'mua_choiceworld'],'batch_vars');

%% Batch striatum responses to choiceworld (NEW - TRIAL ID)

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
        
        % Group multiunit by depth
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
        
        depth_group = discretize(spikeDepths,depth_group_edges);
        
        % Set times for PSTH
        raster_window = [-0.5,3];
        psth_bin_size = 0.001;
        t = raster_window(1):psth_bin_size:raster_window(2);
        
        % Get trial properties
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Loop through all trial conditions, PSTH        
        depth_psth = nan(n_conditions,length(t)-1,n_depths,2);
        for curr_align = 1:2
            switch curr_align
                case 1
                    use_align = stimOn_times;
                case 2
                    use_align = wheel_move_time';
                    use_align(isnan(use_align)) = 0;
            end
            t_peri_event = bsxfun(@plus,use_align,t);
            for curr_depth = 1:n_depths
                
                curr_spikes = spike_times_timeline(depth_group == curr_depth);
                
                curr_spikes_binned = cell2mat(arrayfun(@(x) ...
                    histcounts(curr_spikes,t_peri_event(x,:)), ...
                    [1:size(t_peri_event,1)]','uni',false))./psth_bin_size;
                
                [curr_ids_str,curr_mean_psth] = ...
                    grpstats(curr_spikes_binned(use_trials,:), ...
                    trial_id(use_trials),{'gname','mean'});
                curr_ids = cellfun(@str2num,curr_ids_str);
                
                depth_psth(curr_ids,:,curr_depth,curr_align) = curr_mean_psth;
            end
        end
        
        % Count trials per condition
        condition_counts = histcounts(trial_id(use_trials), ...
            'BinLimits',[1,n_conditions],'BinMethod','integers')';
       
        % Store
        batch_vars(curr_animal).depth_psth(:,:,:,:,curr_day) = depth_psth; 
        batch_vars(curr_animal).condition_counts(:,curr_day) = condition_counts; 
        
        % Prep for next loop
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal])
    
end

% Save
save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
save([save_path filesep 'mua_choiceworld'],'batch_vars');

disp(['Finished batch'])

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
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
        
        depth_group = discretize(spikeDepths,depth_group_edges);
        
        raster_window = [-0.5,5];
        psth_bin_size = 0.001;
        t = raster_window(1):psth_bin_size:raster_window(2);
        t_bins = t(1:end-1) + diff(t);
        
        mua_stim = nan(6,length(t_bins),length(conditions));
        for curr_depth = 1:n_depths
            
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
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
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
        event_aligned_spikes = nan(size(t_peri_event_bins,1),size(t_peri_event_bins,2)-1,n_depths);
        for curr_depth = 1:n_depths
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
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
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
        event_aligned_spikes = nan(size(t_peri_event_bins,1),size(t_peri_event_bins,2)-1,n_depths);
        for curr_depth = 1:n_depths
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


%% Batch trial-trial correlations (UPDATED: within-condition shuffle, across alignment and timing conditions)

for trialtype_align = {'stim','move'};
    trialtype_align = cell2mat(trialtype_align);
    for trialtype_timing = {'earlymove','latemove'};
        trialtype_timing = cell2mat(trialtype_timing);
                
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
                
                switch trialtype_align
                    case 'stim'
                        use_align = stimOn_times;
                    case 'move'
                        use_align = wheel_move_time;
                end
                
                switch trialtype_timing
                    case 'earlymove'
                        use_stim_to_move = stim_to_move < 0.5;
                        use_stim_to_feedback = stim_to_feedback < 1.5;
                    case 'latemove'
                        use_stim_to_move = stim_to_move > 0.5;
                        use_stim_to_feedback = stim_to_feedback < 1.5;
                    case 'nomove'
                        use_stim_to_move = stim_to_move > 0.5;
                        use_stim_to_feedback = stim_to_feedback >= 1.5;
                end
                
                
                % Load widefield ROIs
                wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
                load(wf_roi_fn);
                
                % Group striatum depths
                n_depths = 6;
                depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
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
                    use_stim_to_move & ...
                    use_stim_to_feedback < 1.5;
                
                align_times = reshape(use_align(use_trials),[],1);
                
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
                event_aligned_spikes = nan(size(t_peri_event_bins,1),size(t_peri_event_bins,2)-1,n_depths);
                for curr_depth = 1:n_depths
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
                clearvars -except animals animal curr_animal  ...
                    protocol experiments curr_day animal batch_vars load_parts ...
                    trialtype_align trialtype_timing
                
            end
            
            disp(['Finished ' animal]);
            
        end
        
        save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\';
        save_fn = [save_path filesep 'corr_mua_fluor_' trialtype_align '_' trialtype_timing '_conditionshuff'];
        save(save_fn,'batch_vars','-v7.3');
        disp(['Saved ' save_fn]);      
        
    end
end


%% Batch trial-trial CHOICE DIFF ONLY (within-condition shuffle)

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
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        [depth_group_n,~,depth_group] = histcounts(spikeDepths,depth_group_edges);
     
        % Define times to align   
        sample_rate_factor = 3;
        interval_surround = [-0.5,1.5];
        sample_rate = framerate*sample_rate_factor;
        t_surround = interval_surround(1):1/sample_rate:interval_surround(2);
        
        % Get p-values for choice-specific differences
        % [region,t,alignment,timing]
        mua_choice_p = nan(n_depths,length(t_surround),2,2);
        fluor_choice_p = nan(numel(wf_roi),length(t_surround),2,2);
        
        for curr_align = 1:2
            
            switch curr_align
                case 1
                    align_times = stimOn_times;
                case 2
                    align_times = wheel_move_time';
                    align_times(isnan(align_times)) = 0;
            end
            
            for curr_timing = 1:2
                
                switch curr_timing
                    case 1
                        use_stim_to_move = stim_to_move < 0.5;
                    case 2
                        use_stim_to_move = stim_to_move > 0.5;
                end
                
                % Get trials to use and surrounding time points
                use_trials = ...
                    trial_outcome ~= 0 & ...
                    ~signals_events.repeatTrialValues(1:n_trials) & ...
                    stim_to_feedback < 1.5 & ...
                    use_stim_to_move;
                
                t_peri_event = bsxfun(@plus,align_times(use_trials),t_surround);
                
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
                event_aligned_spikes = nan(size(t_peri_event_bins,1),size(t_peri_event_bins,2)-1,n_depths);
                for curr_depth = 1:n_depths
                    use_spikes = spike_times_timeline(depth_group == curr_depth);
                    event_aligned_spikes(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                        histcounts(use_spikes,t_peri_event_bins(x,:)),[1:size(t_peri_event_bins,1)]','uni',false));
                end
                
                % Get sig correlations across all time points across modalities
                % (SHUFFLE ONLY WITHIN CONTRAST/SIDE)
                n_shuff = 1000;
                warning off;
                use_comb = combvec(contrasts,sides)';
                data_idx = reshape(1:sum(use_trials)*length(t_surround), ...
                    sum(use_trials),length(t_surround));
                shuff_idx = nan(sum(use_trials),length(t_surround),n_shuff);
                for curr_condition = 1:size(use_comb,1)
                    curr_trials = ismember(trial_conditions(use_trials,1:2),use_comb(curr_condition,:),'rows');
                    shuff_idx(curr_trials,:,:) = ...
                        shake(repmat(data_idx(curr_trials,:,:),1,1,n_shuff),1);
                end
                warning on;                
                
                % MUA-choice               
                for curr_mua = 1:size(event_aligned_spikes,3)
                    
                    curr_data1 = event_aligned_spikes(:,:,curr_mua);
                    curr_data2 = trial_choice(use_trials)';
                    curr_data2_shuff = curr_data2(shuff_idx(:,1,:));
                    
                    corr_real = curr_data1'*curr_data2;
                    
                    corr_shuff = ...
                        gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                        gpuArray(curr_data2_shuff)));
                    
                    corr_rank = tiedrank([corr_real';permute(corr_shuff,[3,1,2])]);
                    corr_p = corr_rank(1,:)/(n_shuff+1);
                    
                    mua_choice_p(curr_mua,:,curr_align,curr_timing) = corr_p;
                    
                end
                
                % Fluor-choice
                for curr_fluor = 1:size(event_aligned_df,3)
                    
                    curr_data1 = event_aligned_df(:,:,curr_fluor);
                    curr_data2 = trial_choice(use_trials)';
                    curr_data2_shuff = curr_data2(shuff_idx(:,1,:));
                    
                    corr_real = curr_data1'*curr_data2;
                    
                    corr_shuff = ...
                        gather(pagefun(@mtimes,gpuArray(curr_data1'), ...
                        gpuArray(curr_data2_shuff)));
                    
                    corr_rank = tiedrank([corr_real';permute(corr_shuff,[3,1,2])]);
                    corr_p = corr_rank(1,:)/(n_shuff+1);
                    
                    fluor_choice_p(curr_fluor,:,curr_align,curr_timing) = corr_p;
                    
                end               
            end
        end
        
        batch_vars(curr_animal).mua_choice_p{curr_day} = mua_choice_p;       
        batch_vars(curr_animal).fluor_choice_p{curr_day} = fluor_choice_p;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal]);
    
end

fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\choice_p';
save(fn,'batch_vars','-v7.3');


%% Batch Peter-modeling

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

protocol = 'vanillaChoiceworld';

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
               
        % Prepare fluorescence
        % (load widefield ROIs)
        wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
        load(wf_roi_fn);
        n_rois = numel(wf_roi);
        
        aUdf = single(AP_align_widefield(animal,day,Udf));
        roi_trace = AP_svd_roi(aUdf,fVdf,[],[],cat(3,wf_roi.mask));
        
        % (get ddf)
        roi_traces_derivative = diff(roi_trace,[],2);
        roi_traces_derivative(roi_traces_derivative < 0) = 0;
        
        % (make L-R traces)
        roi_traces_derivative(size(wf_roi,1)+1:end,:) = ...
            roi_traces_derivative(1:size(wf_roi,1),:) - roi_traces_derivative(size(wf_roi,1)+1:end,:);
        
        % Prepare MUA
        % (group striatum depths)
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        depth_group = discretize(spikeDepths,depth_group_edges);
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5 & ...
            stim_to_move < 0.5;
        
        % Get event-aligned activity        
        raster_window = [-0.5,1];
        upsample_factor = 3;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        event_aligned_ddf = nan(length(stimOn_times),length(t),n_rois,2);
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths,2);
        for curr_align = 1:2
            switch curr_align
                case 1
                    use_align = stimOn_times;
                    use_align(isnan(use_align)) = 0;
                case 2
                    use_align = wheel_move_time';
                    use_align(isnan(use_align)) = 0;
            end
            
            t_peri_event = bsxfun(@plus,use_align,t);
            
            % Fluorescence
            event_aligned_ddf(:,:,:,curr_align) = ...
                interp1(conv2(frame_t,[1,1]/2,'valid'),roi_traces_derivative',t_peri_event);
            
            % MUA
            % (raster times)
            t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
            for curr_depth = 1:n_depths
                curr_spikes = spike_times_timeline(depth_group == curr_depth);
                event_aligned_mua(:,:,curr_depth,curr_align) = cell2mat(arrayfun(@(x) ...
                    histcounts(curr_spikes,t_bins(x,:)), ...
                    [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
            end
        end
        
        % Smooth MUA
        smooth_size = 20;
        gw = gausswin(smooth_size,3)';
        smWin = gw./sum(gw);
        event_aligned_mua = convn(event_aligned_mua,smWin,'same');
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = ((trial_choice(use_trials)'+1)/2)+1;
        D.repeatNum = ones(sum(use_trials),1);
        
        % 2) Fit behavioural model
        use_model = 'AP_test_noneur';
        g_stim_all = GLM(D).setModel(use_model).fit;
        behavParameterFit = g_stim_all.parameterFits;
        
        % 3) Do a cross-validated fit for the behavioural model, to provide a
        % measure of the baseline behavioural model log likelihood
        g_stim_cv = GLM(D).setModel(use_model).fitCV(10);
        pL = g_stim_cv.p_hat(:,1);    pR = g_stim_cv.p_hat(:,2);
        likelihood = pL.*(g_stim_cv.data.response==1) + pR.*(g_stim_cv.data.response==2);
        loglik_bpt_stim = nanmean(log2(likelihood(likelihood ~= 0)));
        
        % 3.5) Get a cross-validated probability in an "empirical model" to
        % set a quasi-upper bound on predictability
        cvfold = 10;
        trial_partition_ordered = round(linspace(1,cvfold,sum(use_trials)))';
        trial_partition = trial_partition_ordered(randperm(length(trial_partition_ordered)));
        
        contrasts = [0,0.06,0.125,0.25,0.5,1];
        contrast_sides = unique([-1*contrasts,contrasts]);
        contrast_side_trial = diff(D.stimulus,[],2);
        [~,contrast_sides_id] = ismember(contrast_side_trial,contrast_sides','rows');
        
        pL_hat = nan(sum(use_trials),1);
        for curr_cv = 1:cvfold   
            
            train_trials = trial_partition ~= curr_cv;
            test_trials = trial_partition == curr_cv;
            
            go_left_n_train = accumarray(contrast_sides_id(train_trials), ...
                D.response(train_trials) == 1,[length(contrast_sides),1]);
            contrast_sides_n_train = accumarray(contrast_sides_id(train_trials),1, ...
                [length(contrast_sides),1]);
            go_left_frac_train = go_left_n_train./contrast_sides_n_train;
            
            pL_hat(test_trials) = go_left_frac_train(contrast_sides_id(test_trials));
            
        end
        
        likelihood = pL_hat.*(D.response==1) + (1-pL_hat).*(D.response==2);
        loglik_bpt_empirical = mean(log2(likelihood(likelihood ~= 0)));
        loglik_increase_empirical = loglik_bpt_empirical - loglik_bpt_stim;
        
        % 4) Now do a cross-validated fit for the neural activity. Use the behavioural model
        % parameter (non-crossval fit) to provide an 'offset' for each trial, which
        % reflects the contribution of the behavioural model
        warning off
        loglik_bpt_fluor = nan(n_rois,length(t),2);
        for curr_align = 1:2
            for curr_roi = 1:n_rois
                for curr_t = 1:length(t)
                    D.offset_ZL = g_stim_all.ZL(behavParameterFit, g_stim_all.Zinput(g_stim_all.data));
                    D.neur = reshape(event_aligned_ddf(use_trials,curr_t,curr_roi,curr_align),[],1);
                    g_neur = GLM(D).setModel('AP_test_neur_stimoffset').fitCV(10);
                    pL = g_neur.p_hat(:,1);
                    pR = g_neur.p_hat(:,2);
                    likelihood = pL.*(g_neur.data.response==1) + pR.*(g_neur.data.response==2);
                    loglik_bpt_fluor(curr_roi,curr_t,curr_align) = mean(log2(likelihood(likelihood ~= 0)));
                end
            end
        end
        loglik_increase_fluor = loglik_bpt_fluor - loglik_bpt_stim;
        
        loglik_bpt_mua = nan(n_depths,length(t),2);
        for curr_align = 1:2
            for curr_depth = 1:n_depths
                for curr_t = 1:length(t)
                    D.offset_ZL = g_stim_all.ZL(behavParameterFit, g_stim_all.Zinput(g_stim_all.data));
                    D.neur = reshape(event_aligned_mua(use_trials,curr_t,curr_depth,curr_align),[],1);
                    g_neur = GLM(D).setModel('AP_test_neur_stimoffset').fitCV(10);
                    pL = g_neur.p_hat(:,1);
                    pR = g_neur.p_hat(:,2);
                    likelihood = pL.*(g_neur.data.response==1) + pR.*(g_neur.data.response==2);
                    loglik_bpt_mua(curr_depth,curr_t,curr_align) = mean(log2(likelihood(likelihood ~= 0)));
                end
            end
        end
        loglik_increase_mua = loglik_bpt_mua - loglik_bpt_stim;
        warning on;
        
        % Store
        batch_vars(curr_animal).loglik_bpt_stim(curr_day) = loglik_bpt_stim;
        batch_vars(curr_animal).loglik_increase_empirical(curr_day) = loglik_increase_empirical;
        batch_vars(curr_animal).loglik_increase_fluor(:,:,:,curr_day) = loglik_increase_fluor;
        batch_vars(curr_animal).loglik_increase_mua(:,:,:,curr_day) = loglik_increase_mua;       
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal]);
    
end

fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\loglik_increase_early';
save(fn,'batch_vars','-v7.3');
disp('Finished batch');

%% Batch modeling: fit C, fit C + N, get LL difference (MEAN TIME)

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

protocol = 'vanillaChoiceworld';

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
        
        % Prepare fluorescence
        % (load widefield ROIs)
        wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
        load(wf_roi_fn);
        n_rois = numel(wf_roi);
        
        aUdf = single(AP_align_widefield(animal,day,Udf));
        roi_trace = AP_svd_roi(aUdf,fVdf,[],[],cat(3,wf_roi.mask));
        
        % (get ddf)
        roi_traces_derivative = diff(roi_trace,[],2);
        roi_traces_derivative(roi_traces_derivative < 0) = 0;
        
        % (make L-R traces)
        roi_traces_derivative(size(wf_roi,1)+1:end,:) = ...
            roi_traces_derivative(1:size(wf_roi,1),:) - roi_traces_derivative(size(wf_roi,1)+1:end,:);
        
        % Prepare MUA
        % (group striatum depths)
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        depth_group = discretize(spikeDepths,depth_group_edges);
        
        % Get event-aligned activity
        raster_window = [-0.5,1];
        upsample_factor = 3;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        event_aligned_ddf = nan(length(stimOn_times),length(t),n_rois,2);
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths,2);
        for curr_align = 1:2
            switch curr_align
                case 1
                    use_align = stimOn_times;
                    use_align(isnan(use_align)) = 0;
                case 2
                    use_align = wheel_move_time';
                    use_align(isnan(use_align)) = 0;
            end
            
            t_peri_event = bsxfun(@plus,use_align,t);
            
            % Fluorescence
            event_aligned_ddf(:,:,:,curr_align) = ...
                interp1(conv2(frame_t,[1,1]/2,'valid'),roi_traces_derivative',t_peri_event);
            
            % MUA
            % (raster times)
            t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
            for curr_depth = 1:n_depths
                curr_spikes = spike_times_timeline(depth_group == curr_depth);
                event_aligned_mua(:,:,curr_depth,curr_align) = cell2mat(arrayfun(@(x) ...
                    histcounts(curr_spikes,t_bins(x,:)), ...
                    [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
            end
        end
        
        % Smooth MUA
        smooth_size = 20;
        gw = gausswin(smooth_size,3)';
        smWin = gw./sum(gw);
        event_aligned_mua = convn(event_aligned_mua,smWin,'same');
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5 & ...
            stim_to_move < 0.5;
                
        cvfold = 10;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = ((trial_choice(use_trials)'+1)/2)+1;
        D.repeatNum = ones(sum(use_trials),1);
        
        % Fit stim crossvalidated
        use_model = 'AP_test_noneur';
        
        g_stim = GLM(D).setModel(use_model).fitCV(cvfold);
        pL = g_stim.p_hat(:,1);    pR = g_stim.p_hat(:,2);
        likelihood = pL.*(g_stim.data.response==1) + pR.*(g_stim.data.response==2);
        
        stim_params = nanmean(g_stim.parameterFits,1);
        loglik_bpt_stim = nanmean(log2(likelihood(likelihood ~= 0)));
        
        % Fit stim all
        use_model = 'AP_test_noneur';
        g_stim_all = GLM(D).setModel(use_model).fit;
        behavParameterFit = g_stim_all.parameterFits;
        
        % Fit stim w/ neural
        use_model = 'AP_test_neur_stimoffset';
        
        use_t_stim = t > 0.05 & t < 0.15;
        use_t_move = t > -0.15 & t < -0.02;
        
        warning off
        fluor_params = nan(n_rois,2,2);
        loglik_bpt_fluor = nan(n_rois,2);
        for curr_align = 1:2
            switch curr_align
                case 1
                    use_t = use_t_stim;
                case 2
                    use_t = use_t_move;
            end
            for curr_roi = 1:n_rois
                D.offset_ZL = g_stim_all.ZL(behavParameterFit, g_stim_all.Zinput(g_stim_all.data));
                D.neur = reshape(nanmean(event_aligned_ddf(use_trials,use_t,curr_roi,curr_align),2),[],1);
                clear g_fluor
                g_fluor = GLM(D).setModel(use_model).fitCV(cvfold);
                pL = g_fluor.p_hat(:,1);
                pR = g_fluor.p_hat(:,2);
                likelihood = pL.*(g_fluor.data.response==1) + pR.*(g_fluor.data.response==2);
                
                fluor_params(curr_roi,:,curr_align) = nanmean(g_fluor.parameterFits,1);
                loglik_bpt_fluor(curr_roi,curr_align) = nanmean(log2(likelihood(likelihood ~= 0)));
            end
        end
        loglik_increase_fluor = (loglik_bpt_fluor - loglik_bpt_stim);
        
        mua_params = nan(n_depths,2,2);
        loglik_bpt_mua = nan(n_depths,2);
        for curr_align = 1:2
            switch curr_align
                case 1
                    use_t = use_t_stim;
                case 2
                    use_t = use_t_move;
            end
            for curr_depth = 1:n_depths
                D.offset_ZL = g_stim_all.ZL(behavParameterFit, g_stim_all.Zinput(g_stim_all.data));
                D.neur = reshape(nanmean(event_aligned_mua(use_trials,use_t,curr_depth,curr_align),2),[],1);
                clear g_mua
                g_mua = GLM(D).setModel(use_model).fitCV(cvfold);
                pL = g_mua.p_hat(:,1);
                pR = g_mua.p_hat(:,2);
                likelihood = pL.*(g_mua.data.response==1) + pR.*(g_mua.data.response==2);
                
                mua_params(curr_depth,:,curr_align) = nanmean(g_mua.parameterFits,1);
                loglik_bpt_mua(curr_depth,curr_align) = nanmean(log2(likelihood(likelihood ~= 0)));
            end
        end
        loglik_increase_mua = (loglik_bpt_mua - loglik_bpt_stim);
        warning on
        
        % Store
        batch_vars(curr_animal).stim_params(curr_day,:) = stim_params;
        batch_vars(curr_animal).fluor_params(:,:,:,curr_day) = fluor_params;
        batch_vars(curr_animal).mua_params(:,:,:,curr_day) = mua_params;
        batch_vars(curr_animal).loglik_increase_fluor(:,:,curr_day) = loglik_increase_fluor;
        batch_vars(curr_animal).loglik_increase_mua(:,:,curr_day) = loglik_increase_mua;       
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal]);
    
end

fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\c_cn_modeling_ll_earlymove';
save(fn,'batch_vars','-v7.3');

disp('Finished batch');

%% Batch modeling: fit C, fit C + N, get LL difference (MEAN TIME, L/R HEMISPHERE)
% using all trial timings at the moment

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

protocol = 'vanillaChoiceworld';

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
        
        % Prepare fluorescence
        % (load widefield ROIs)
        wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
        load(wf_roi_fn);
        n_rois = numel(wf_roi);
        
        aUdf = single(AP_align_widefield(animal,day,Udf));
        roi_trace = AP_svd_roi(aUdf,fVdf,[],[],cat(3,wf_roi.mask));
        
        % (get ddf)
        roi_traces_derivative = diff(roi_trace,[],2);
        roi_traces_derivative(roi_traces_derivative < 0) = 0;
                
        % Prepare MUA
        % (group striatum depths)
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        depth_group = discretize(spikeDepths,depth_group_edges);
        
        % Get event-aligned activity
        raster_window = [-0.5,1];
        upsample_factor = 3;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        event_aligned_ddf = nan(length(stimOn_times),length(t),n_rois,2);
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths,2);
        for curr_align = 1:2
            switch curr_align
                case 1
                    use_align = stimOn_times;
                    use_align(isnan(use_align)) = 0;
                case 2
                    use_align = wheel_move_time';
                    use_align(isnan(use_align)) = 0;
            end
            
            t_peri_event = bsxfun(@plus,use_align,t);
            
            % Fluorescence
            event_aligned_ddf(:,:,:,curr_align) = ...
                interp1(conv2(frame_t,[1,1]/2,'valid'),roi_traces_derivative',t_peri_event);
            
            % MUA
            % (raster times)
            t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
            for curr_depth = 1:n_depths
                curr_spikes = spike_times_timeline(depth_group == curr_depth);
                event_aligned_mua(:,:,curr_depth,curr_align) = cell2mat(arrayfun(@(x) ...
                    histcounts(curr_spikes,t_bins(x,:)), ...
                    [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
            end
        end
        
        % Smooth MUA
        smooth_size = 20;
        gw = gausswin(smooth_size,3)';
        smWin = gw./sum(gw);
        event_aligned_mua = convn(event_aligned_mua,smWin,'same');
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5 & ...
            stim_to_move < Inf;
        
        cvfold = 10;
        
        % 1) Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues == -1;
        R_trials = signals_events.trialSideValues == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = ((trial_choice(use_trials)'+1)/2)+1;
        D.repeatNum = ones(sum(use_trials),1);
        
        % 2) Fit without neural
        use_model = 'AP_test_noneur';
        
        g_stim = GLM(D).setModel(use_model).fitCV(cvfold);
        pL = g_stim.p_hat(:,1);    pR = g_stim.p_hat(:,2);
        likelihood = pL.*(g_stim.data.response==1) + pR.*(g_stim.data.response==2);
        
        stim_params = nanmean(g_stim.parameterFits,1);
        loglik_bpt_stim = nanmean(log2(likelihood(likelihood ~= 0)));
        
        % 3) Fit with neural        
        use_t_stim = t > 0.05 & t < 0.15;
        use_t_move = t > -0.15 & t < -0.02;
        
        fluor_params = nan(n_rois/2,6,2);
        loglik_bpt_fluor = nan(n_rois/2,2);
        for curr_align = 1:2
            switch curr_align
                case 1
                    use_t = use_t_stim;
                case 2
                    use_t = use_t_move;
            end
            for curr_roi = 1:n_rois/2
                D.neur = [reshape(nanmean(event_aligned_ddf(use_trials,use_t,curr_roi,curr_align),2),[],1), ...
                    reshape(nanmean(event_aligned_ddf(use_trials,use_t,curr_roi+size(wf_roi,1),curr_align),2),[],1)];
                clear g_fluor
                g_fluor = GLM(D).setModel('AP_test_neur_hemisphere').fitCV(cvfold);
                pL = g_fluor.p_hat(:,1);
                pR = g_fluor.p_hat(:,2);
                likelihood = pL.*(g_fluor.data.response==1) + pR.*(g_fluor.data.response==2);
                
                fluor_params(curr_roi,:,curr_align) = nanmean(g_fluor.parameterFits,1);
                loglik_bpt_fluor(curr_roi,curr_align) = nanmean(log2(likelihood(likelihood ~= 0)));
            end
        end
        loglik_increase_fluor = (loglik_bpt_fluor - loglik_bpt_stim);
        
        mua_params = nan(n_depths,5,2);
        loglik_bpt_mua = nan(n_depths,2);
        for curr_align = 1:2
            switch curr_align
                case 1
                    use_t = use_t_stim;
                case 2
                    use_t = use_t_move;
            end
            for curr_depth = 1:n_depths
                D.neur = reshape(nanmean(event_aligned_mua(use_trials,use_t,curr_depth,curr_align),2),[],1);
                clear g_mua
                g_mua = GLM(D).setModel('AP_test_neur').fitCV(cvfold);
                pL = g_mua.p_hat(:,1);
                pR = g_mua.p_hat(:,2);
                likelihood = pL.*(g_mua.data.response==1) + pR.*(g_mua.data.response==2);
                
                mua_params(curr_depth,:,curr_align) = nanmean(g_mua.parameterFits,1);
                loglik_bpt_mua(curr_depth,curr_align) = nanmean(log2(likelihood(likelihood ~= 0)));
            end
        end
        loglik_increase_mua = (loglik_bpt_mua - loglik_bpt_stim);
        
        % Store
        batch_vars(curr_animal).stim_params(curr_day,:) = stim_params;
        batch_vars(curr_animal).fluor_params(:,:,:,curr_day) = fluor_params;
        batch_vars(curr_animal).mua_params(:,:,:,curr_day) = mua_params;
        batch_vars(curr_animal).loglik_bpt_stim(:,curr_day) = loglik_bpt_stim;      
        batch_vars(curr_animal).loglik_increase_fluor(:,:,curr_day) = loglik_increase_fluor;
        batch_vars(curr_animal).loglik_increase_mua(:,:,curr_day) = loglik_increase_mua;       
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal]);
    
end

fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\c_cn_modeling_ll_hemisphere';
save(fn,'batch_vars','-v7.3');


%% Batch logistic choice modeling - get parameters from Vs

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

protocol = 'vanillaChoiceworld';

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
    load_parts.ephys = false;
    
    for curr_day = 1:length(experiments);
        
        % Load experiment
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        AP_load_experiment
        
        % Align U's
        aUdf = single(AP_align_widefield(animal,day,Udf));
        
        % Aligned V's to trial events
        raster_window = [-0.5,1];
        upsample_factor = 3;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        event_aligned_V = nan(length(stimOn_times),length(t),size(U,3),2);
        for curr_align = 1:2
            
            % Align by stim or movement
            switch curr_align
                case 1
                    use_align = stimOn_times;
                case 2
                    use_align = wheel_move_time';
                    use_align(isnan(use_align)) = 0;
            end
            
            % Get event-aligned Vs
            t_peri_event = bsxfun(@plus,use_align,t);
            event_aligned_V(:,:,:,curr_align) = ...
                interp1(conv2(frame_t,[1,1]/2,'valid'),diff(fVdf,[],2)',t_peri_event);
            
        end
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5 & ...
            stim_to_move < 0.5;
                
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues == -1;
        R_trials = signals_events.trialSideValues == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = ((trial_choice(use_trials)'+1)/2)+1;
        D.repeatNum = ones(sum(use_trials),1);
        
        % Fit fluor model (across all timepoints)
        cvfold = 5;
        use_Vs = 1:100;
        
        pV = nan(size(aUdf,1),size(aUdf,2),length(t),2);
        for curr_align = 1:2
            for curr_t = 1:length(t)
                
%                 D.offset_ZL = g.ZL(behavParameterFit, g.Zinput(g.data));
                D.neur = double(squeeze(event_aligned_V(use_trials,curr_t,use_Vs,curr_align)));
                g_neur = GLM(D).setModel('AP_test_V').fitCV(cvfold);
%                 pL = g_neur.p_hat(:,1);
%                 pR = g_neur.p_hat(:,2);
%                 likelihood = pL.*(g_neur.data.response==1) + pR.*(g_neur.data.response==2);
%                 loglik_bpt_fluor = mean(log2(likelihood));
                
                pV_mean_cv = nanmean(g_neur.parameterFits(:,2:end),1);
                pV(:,:,curr_t,curr_align) = svdFrameReconstruct(aUdf(:,:,use_Vs),pV_mean_cv');
            end
        end
        
        % Store
        batch_vars(curr_animal).pV(:,:,:,:,curr_day) = pV;       
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal]);
    
end

fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\logistic_regression_pV';
save(fn,'batch_vars','-v7.3');

disp('Finished batch')

%% Batch average activity within times of interest for all trial types
% using all trial timings at the moment

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

protocol = 'vanillaChoiceworld';

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
        
        % Prepare fluorescence
        % (load widefield ROIs)
        wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
        load(wf_roi_fn);
        n_rois = numel(wf_roi);
        
        aUdf = single(AP_align_widefield(animal,day,Udf));
        roi_trace = AP_svd_roi(aUdf,fVdf,[],[],cat(3,wf_roi.mask));
        
        % (get ddf)
        roi_traces_derivative = diff(roi_trace,[],2);
        roi_traces_derivative(roi_traces_derivative < 0) = 0;
        
        % (make L-R traces)
        roi_traces_derivative(size(wf_roi,1)+1:end,:) = ...
            roi_traces_derivative(1:size(wf_roi,1),:) - roi_traces_derivative(size(wf_roi,1)+1:end,:);
        
        % Prepare MUA
        % (group striatum depths)
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        depth_group = discretize(spikeDepths,depth_group_edges);
        
        % Get event-aligned activity
        raster_window = [-0.5,1];
        upsample_factor = 3;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        event_aligned_ddf = nan(length(stimOn_times),length(t),n_rois,2);
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths,2);
        for curr_align = 1:2
            switch curr_align
                case 1
                    use_align = stimOn_times;
                    use_align(isnan(use_align)) = 0;
                case 2
                    use_align = wheel_move_time';
                    use_align(isnan(use_align)) = 0;
            end
            
            t_peri_event = bsxfun(@plus,use_align,t);
            
            % Fluorescence
            event_aligned_ddf(:,:,:,curr_align) = ...
                interp1(conv2(frame_t,[1,1]/2,'valid'),roi_traces_derivative',t_peri_event);
            
            % MUA
            % (raster times)
            t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
            for curr_depth = 1:n_depths
                curr_spikes = spike_times_timeline(depth_group == curr_depth);
                event_aligned_mua(:,:,curr_depth,curr_align) = cell2mat(arrayfun(@(x) ...
                    histcounts(curr_spikes,t_bins(x,:)), ...
                    [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
            end
        end
        
        % Smooth MUA
        smooth_size = 20;
        gw = gausswin(smooth_size,3)';
        smWin = gw./sum(gw);
        event_aligned_mua = convn(event_aligned_mua,smWin,'same');
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Pick times to average across
        use_t_stim = t > 0.05 & t < 0.15;
        use_t_move = t > -0.15 & t < -0.02;
        use_t_beep = t > 0.55 & t < 0.65;
        
        % Get average activity for both alignments
        event_aligned_ddf_avg = nan(sum(use_trials),n_rois,3);
        event_aligned_mua_avg = nan(sum(use_trials),n_depths,3);
        
        for curr_align = 1:3
            switch curr_align
                case 1
                    use_align = 1;
                    use_t = use_t_stim;
                case 2
                    use_align = 2;
                    use_t = use_t_move;
                case 3
                    use_align = 1;
                    use_t = use_t_beep;
            end
            event_aligned_ddf_avg(:,:,curr_align) = ...
                squeeze(nanmean(event_aligned_ddf(use_trials,use_t,:,use_align),2));
            event_aligned_mua_avg(:,:,curr_align) = ...
                squeeze(nanmean(event_aligned_mua(use_trials,use_t,:,use_align),2));
        end
        
        % Package all activity by trial ID (must be more elegant way but whatever)
        fluor_trial_act = cell(n_rois,n_conditions,3);
        for curr_roi = 1:n_rois
            for curr_align = 1:3
                for curr_condition = 1:n_conditions
                    fluor_trial_act{curr_roi,curr_condition,curr_align} = ...
                        event_aligned_ddf_avg(trial_id(use_trials) == curr_condition,curr_roi,curr_align);
                end
            end
        end
        
        mua_trial_act = cell(n_depths,n_conditions,3);
        for curr_depth = 1:n_depths
            for curr_align = 1:3
                for curr_condition = 1:n_conditions
                    mua_trial_act{curr_depth,curr_condition,curr_align} = ...
                        event_aligned_mua_avg(trial_id(use_trials) == curr_condition,curr_depth,curr_align);
                end
            end
        end
        
        % Store
        batch_vars(curr_animal).fluor_trial_act(:,:,:,curr_day) = fluor_trial_act;
        batch_vars(curr_animal).mua_trial_act(:,:,:,curr_day) = mua_trial_act;
                
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal]);
    
end

fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\trial_activity';
save(fn,'batch_vars','-v7.3');

disp('Finished batch');


%% Batch kernel regression of activity from task parameters
%%%%%% THIS IS A WORK IN PROGRESS
% I started this because I was thinking that kernel regression would be
% cleaner because time points wouldn't be independent, but they still would
% be, so this doesn't add anything. On the other hand, it's cleaner and can
% give partial variance explained, so it's worth eventually moving the
% stuff above to be stand-alone down here.
%%%%%%

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

protocol = 'vanillaChoiceworld';

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
        
        % Prepare fluorescence
        % (load widefield ROIs)
        wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
        load(wf_roi_fn);
        n_rois = numel(wf_roi);
        
        aUdf = single(AP_align_widefield(animal,day,Udf));
        roi_trace = AP_svd_roi(aUdf,fVdf,[],[],cat(3,wf_roi.mask));
        
        % (get ddf)
        roi_traces_derivative = diff(roi_trace,[],2);
        roi_traces_derivative(roi_traces_derivative < 0) = 0;
        
        % (make L-R traces)
        roi_traces_derivative(size(wf_roi,1)+1:end,:) = ...
            roi_traces_derivative(1:size(wf_roi,1),:) - roi_traces_derivative(size(wf_roi,1)+1:end,:);
        
        % Prepare MUA
        % (group striatum depths)
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        depth_group = discretize(spikeDepths,depth_group_edges);
        
        
        
        
        
        
        
        %%%%%%% Task regressors
        contrast_sides = zeros(length(stimOn_times),2);
        contrast_sides(trial_conditions(:,2) == -1,1) = ...
            trial_conditions(trial_conditions(:,2) == -1,1);
        contrast_sides(trial_conditions(:,2) == 1,2) = ...
            trial_conditions(trial_conditions(:,2) == 1,1);
        
        contrast_exp = 0.3;
        task_regressors = [ones(length(stimOn_times),1), ...
                            contrast_sides.^contrast_exp, ...
                            trial_choice(1:n_trials)'];
        
        
        
        %%%%%%% MUA
        
        % Set times for PSTH
        raster_window = [-0.5,3];
        psth_bin_size = 0.001;
        t = raster_window(1):psth_bin_size:raster_window(2);
        
        % Loop through all trial conditions, PSTH        
        depth_psth = nan(length(stimOn_times),length(t)-1,n_depths,2);
        for curr_align = 1:2
            switch curr_align
                case 1
                    use_align = stimOn_times;
                case 2
                    use_align = wheel_move_time';
                    use_align(isnan(use_align)) = 0;
            end
            t_peri_event = bsxfun(@plus,use_align,t);
            for curr_depth = 1:n_depths
                
                curr_spikes = spike_times_timeline(depth_group == curr_depth);
                
                curr_spikes_binned = cell2mat(arrayfun(@(x) ...
                    histcounts(curr_spikes,t_peri_event(x,:)), ...
                    [1:size(t_peri_event,1)]','uni',false))./psth_bin_size;
                
                depth_psth(:,:,curr_depth,curr_align) = curr_spikes_binned;
                
            end
        end
        
        smooth_size = 50;
        gw = gausswin(smooth_size,3)';
        smWin = gw./sum(gw);
        depth_psth_smooth = convn(depth_psth,smWin,'same');
        
        % Regress out contrast/choice from both real and predicted
        contrast_exp = 0.3;
        
        % [depth,time,align,timing,param]
        activity_model_params = nan(n_depths,length(t),2,2,4);
        predicted_activity_model_params = nan(n_depths,length(t),2,2,4);
        for curr_timing = 1:2
            for curr_align = 1:2
                switch curr_align
                    case 1
                        use_align = stimOn_times;
                    case 2
                        use_align = wheel_move_time';
                        use_align(isnan(use_align)) = 0;
                end
                
                use_timing_trials = use_trials' & trial_conditions(:,4) == curr_timing;
                
                t_peri_event = bsxfun(@plus,use_align(use_timing_trials),t);
                for curr_depth = 1:n_depths
                    
                    curr_spikes_binned = interp1(time_bin_centers,binned_spikes(curr_depth,:),t_peri_event)*sample_rate;
                    curr_spikes_binned_pred = interp1(time_bin_centers,predicted_spikes_reranged(curr_depth,:),t_peri_event)*sample_rate;                    
                    
                    curr_contrast_sides = zeros(sum(use_timing_trials),2);
                    curr_contrast_sides(trial_conditions(use_timing_trials,2) == -1,1) = ...
                       trial_conditions(trial_conditions(use_timing_trials,2) == -1,1);
                    curr_contrast_sides(trial_conditions(use_timing_trials,2) == 1,2) = ...
                       trial_conditions(trial_conditions(use_timing_trials,2) == 1,1);
                    
                    curr_choices = trial_conditions(use_timing_trials,3);                                       
                    
                    valid_trials = all(~isnan([curr_spikes_binned,curr_spikes_binned_pred]),2);
                    
                    for curr_t = 1:length(t)                       
                        curr_act = curr_spikes_binned(:,curr_t);
                        params_fit = [ones(sum(valid_trials),1), ...
                            curr_contrast_sides(valid_trials,:).^contrast_exp, ...
                            curr_choices(valid_trials)]\curr_act(valid_trials);
                        activity_model_params(curr_depth,curr_t,curr_align,curr_timing,:) = ...
                            params_fit;
                        
                        curr_act = curr_spikes_binned_pred(:,curr_t);
                        params_fit = [ones(sum(valid_trials),1), ...
                            curr_contrast_sides(valid_trials,:).^contrast_exp, ...
                            curr_choices(valid_trials)]\curr_act(valid_trials);
                        predicted_activity_model_params(curr_depth,curr_t,curr_align,curr_timing,:) = ...
                            params_fit;
                    end
                end
            end
        end
        
        %%%%%%% HERE
        
        
        
        
        
        
        
        
        % Store
        batch_vars(curr_animal).fluor_trial_act(:,:,:,curr_day) = fluor_trial_act;
        batch_vars(curr_animal).mua_trial_act(:,:,:,curr_day) = mua_trial_act;
                
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts
        
    end
    
    disp(['Finished ' animal]);
    
end

fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\trial_activity';
save(fn,'batch_vars','-v7.3');

disp('Finished batch');


%% Batch cortex ROI -> striatum regression

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';
batch_vars = struct;

% Load widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);

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
        
        upsample_factor = 2;
        sample_rate = (1/median(diff(frame_t)))*upsample_factor;
        
        skip_seconds = 60;
        time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
        t = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled roi traces
        aUdf = single(AP_align_widefield(animal,day,Udf));
        roi_traces = AP_svd_roi(aUdf,fVdf,[],[],cat(3,wf_roi.mask));
        n_rois = numel(wf_roi);        

        % Make DDF, upsample
        roi_traces_derivative = diff(roi_traces,[],2);
        roi_traces_derivative(roi_traces_derivative < 0) = 0;
        
        roi_traces_derivative_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            roi_traces_derivative',t)';

        % Group striatum depths
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        [depth_group_n,~,depth_group] = histcounts(spikeDepths,depth_group_edges);
          
        binned_spikes = zeros(n_depths,length(t));
        for curr_depth = 1:n_depths           
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);           
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins)*sample_rate;            
        end
        
        % Align activity, regress from fluor to MUA
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        interval_surround = [-0.5,3];
        t_surround = interval_surround(1):1/sample_rate:interval_surround(2);
        t_surround_bins = [t_surround-(1/sample_rate)/2,t_surround(:,end)+(1/sample_rate)/2];
        
        kernel_t = [-0.3,0.3];
        kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
        zs = [true,true];
        cvfold = 1;
        lambda = 1e2;

        roi_depth_kernel = nan(n_rois,length(kernel_frames),n_depths,2,2);
        activity_model_params = nan(n_depths,length(t_surround),2,2,4);
        predicted_activity_model_params = nan(n_depths,length(t_surround),2,2,4);                       
        for curr_timing = 1:2
            
            use_timing_trials = use_trials' & trial_conditions(:,4) == curr_timing;
            
            for curr_align = 1:2
                switch curr_align
                    case 1
                        use_align = stimOn_times(use_timing_trials);
                    case 2
                        use_align = wheel_move_time(use_timing_trials)';
                        use_align(isnan(use_align)) = 0;
                end
                t_peri_event = bsxfun(@plus,use_align,t_surround);
                t_peri_event_bins = bsxfun(@plus,use_align,t_surround_bins);
                
                % Fluor
                event_aligned_ddf = interp1(conv2(frame_t,[1,1]/2,'valid'),roi_traces_derivative',t_peri_event);
                event_aligned_ddf_reshape = reshape(permute(event_aligned_ddf,[2,1,3]),[],n_rois);
                
                % MUA
                event_aligned_mua = nan(length(use_align),length(t_surround),n_depths);
                for curr_depth = 1:n_depths
                    curr_spikes = spike_times_timeline(depth_group == curr_depth);
                    event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                        histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                        [1:size(t_peri_event,1)]','uni',false))*sample_rate;
                end
                event_aligned_mua_reshape = reshape(permute(event_aligned_mua,[2,1,3]),[],n_depths);
           
                % Regress from fluor to MUA
                [k,predicted_spikes,explained_var] = ...
                    AP_regresskernel(event_aligned_ddf_reshape', ...
                    event_aligned_mua_reshape',kernel_frames,lambda,zs,cvfold);
                
                k_reshape = reshape(k,n_rois,length(kernel_frames),n_depths);      
                
                % Rescale predictions to un-zscore
                mua_std = std(event_aligned_mua_reshape,[],1)';
                mua_mean = mean(event_aligned_mua_reshape,1)';                
                predicted_spikes_reranged = bsxfun(@plus,bsxfun(@times,predicted_spikes, ...
                    mua_std),mua_mean);
                
                event_aligned_mua_predicted = permute(reshape(predicted_spikes_reranged',length(t_surround),length(use_align),n_depths),[2,1,3]);
                
                roi_depth_kernel(:,:,:,curr_align,curr_timing) = fliplr(k_reshape);
                
                % Regress from task parameters to MUA/ fluor-predicted MUA
                contrast_exp = 0.3;
                for curr_depth = 1:n_depths                    
                 
                    curr_contrast_sides = zeros(sum(use_timing_trials),2);
                    curr_contrast_sides(trial_conditions(use_timing_trials,2) == -1,1) = ...
                        trial_conditions(trial_conditions(use_timing_trials,2) == -1,1);
                    curr_contrast_sides(trial_conditions(use_timing_trials,2) == 1,2) = ...
                        trial_conditions(trial_conditions(use_timing_trials,2) == 1,1);
                    
                    curr_choices = trial_conditions(use_timing_trials,3);
                    
                    valid_trials = all(~isnan([event_aligned_mua(:,:,curr_depth), ...
                        event_aligned_mua_predicted(:,:,curr_depth)]),2);
                    
                    for curr_t = 1:length(t_surround)
                        curr_act = event_aligned_mua(:,curr_t,curr_depth);
                        params_fit = [ones(sum(valid_trials),1), ...
                            curr_contrast_sides(valid_trials,:).^contrast_exp, ...
                            curr_choices(valid_trials)]\curr_act(valid_trials);
                        activity_model_params(curr_depth,curr_t,curr_align,curr_timing,:) = ...
                            params_fit;
                        
                        curr_act = event_aligned_mua_predicted(:,curr_t,curr_depth);
                        params_fit = [ones(sum(valid_trials),1), ...
                            curr_contrast_sides(valid_trials,:).^contrast_exp, ...
                            curr_choices(valid_trials)]\curr_act(valid_trials);
                        predicted_activity_model_params(curr_depth,curr_t,curr_align,curr_timing,:) = ...
                            params_fit;
                    end
                end
                
            end
        end
          
        % Store
        batch_vars(curr_animal).roi_depth_kernel(:,:,:,:,:,curr_day) = roi_depth_kernel;        
        batch_vars(curr_animal).activity_model_params(:,:,:,:,:,curr_day) = activity_model_params;
        batch_vars(curr_animal).predicted_activity_model_params(:,:,:,:,:,curr_day) = predicted_activity_model_params;
        
        % Prepare for next loop       
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animals animal curr_animal protocol experiments curr_day animal batch_vars load_parts wf_roi
        
    end
     
    disp(['Finished ' animal])
    
end

save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
save([save_path filesep 'roi-mua_choiceworld_predicted'],'batch_vars');

disp('Finished batch');

%% Batch logistic regression from activity (mean time) - concatenate within animal
% MUA normalization could use some work

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';
batch_vars = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    fluor_all = cell(length(experiments),1);
    mua_all = cell(length(experiments),1);
    D_all = cell(length(experiments),1);
    
    disp(['Loading ' animal]);
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        AP_load_experiment;
        
        % Prepare fluorescence
        % (load widefield ROIs)
        wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
        load(wf_roi_fn);
        n_rois = numel(wf_roi);
        
        aUdf = single(AP_align_widefield(animal,day,Udf));
        roi_trace = AP_svd_roi(aUdf,fVdf,[],[],cat(3,wf_roi.mask));
        
        % (get ddf)
        roi_traces_derivative = diff(roi_trace,[],2);
        roi_traces_derivative(roi_traces_derivative < 0) = 0;
        
        % Prepare MUA
        % (group striatum depths)
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        depth_group = discretize(spikeDepths,depth_group_edges);
        
        % Get event-aligned activity
        raster_window = [-0.5,1];
        upsample_factor = 3;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        event_aligned_ddf = nan(length(stimOn_times),length(t),n_rois,2);
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths,2);
        for curr_align = 1:2
            switch curr_align
                case 1
                    use_align = stimOn_times;
                    use_align(isnan(use_align)) = 0;
                case 2
                    use_align = wheel_move_time';
                    use_align(isnan(use_align)) = 0;
            end
            
            t_peri_event = bsxfun(@plus,use_align,t);
            
            % Fluorescence
            event_aligned_ddf(:,:,:,curr_align) = ...
                interp1(conv2(frame_t,[1,1]/2,'valid'),roi_traces_derivative',t_peri_event);
            
            % MUA
            t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
            for curr_depth = 1:n_depths
                curr_spikes = spike_times_timeline(depth_group == curr_depth);
                event_aligned_mua(:,:,curr_depth,curr_align) = cell2mat(arrayfun(@(x) ...
                    histcounts(curr_spikes,t_bins(x,:)), ...
                    [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
            end
            
        end
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5 & ...
            stim_to_move < 0.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        % D.response = (trial_choice(use_trials)'+1)/2+1;
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);
        
        % Store activity and stim
        fluor_all{curr_day} = event_aligned_ddf(use_trials,:,:,:);
        mua_all{curr_day} = event_aligned_mua(use_trials,:,:,:);
        D_all{curr_day} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
    end
    
    % Get rid of everything that's not the data
    clearvars -except animals protocol batch_vars curr_animal fluor_all mua_all D_all
    
    % Get time
    framerate = 35.2;
    raster_window = [-0.5,1];
    upsample_factor = 3;
    raster_sample_rate = 1/(framerate*upsample_factor);
    t = raster_window(1):raster_sample_rate:raster_window(2);
    
    use_t_stim = t > 0.05 & t < 0.15;
    use_t_move = t > -0.15 & t < -0.02;
    use_t_align = [use_t_stim;use_t_move];
    
    % Get widefield ROIs
    wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
    load(wf_roi_fn);
    n_rois = numel(wf_roi);
    
    n_depths = 6;
    
    cvfold = 20;
    
    % Concatenate activity (and normalize MUA)
    fluor_cat = cat(1,fluor_all{:});
    
    mua_cat_raw = cat(1,mua_all{:});
    t_baseline = t < 0;
    softnorm = 5;
    mua_baseline = nanmean(mua_cat_raw(:,t_baseline,:,1),2);
    mua_cat = bsxfun(@rdivide,mua_cat_raw,mua_baseline + softnorm);
    
    % Get L-R activity
    fluor_cat_hemidiff = cat(3,fluor_cat(:,:,1:n_rois/2,:),fluor_cat(:,:,1:n_rois/2,:) - ...
        fluor_cat(:,:,n_rois/2+1:end,:));
    
    % Concatenate behavioural data
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all,'uni',false));
    D.response = cell2mat(cellfun(@(x) x.response,D_all,'uni',false));
    D.repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all,'uni',false));
    
    % Get zero-contrasts as a subset
    zero_contrasts = D.stimulus(:,1) == 0 & D.stimulus(:,2) == 0;
    D_zero = struct;
    D_zero.stimulus = D.stimulus(zero_contrasts,:);
    D_zero.response = D.response(zero_contrasts);
    D_zero.repeatNum = D.repeatNum(zero_contrasts);
    
    % Fit stim all
    use_model = 'AP_test_stim';
    g_stim_all = GLM(D).setModel(use_model).fit;
    behavParameterFit = g_stim_all.parameterFits;
    
    % Fit stim cross-validated
    use_model = 'AP_test_stim';
    
    g_stim = GLM(D).setModel(use_model).fitCV(cvfold);
    pL = g_stim.p_hat(:,1);    pR = g_stim.p_hat(:,2);
    likelihood = pL.*(g_stim.data.response == 1) + pR.*(g_stim.data.response == 2);
    
    loglik_bpt_stim = nanmean(log2(likelihood));
    
    % Fit stim + activity
    use_model = 'AP_test_neur_stim';
    
    loglik_bpt_fluor = nan(n_rois,2);
    for curr_align = 1:2
        for curr_area = 1:n_rois
            
            D.neur = reshape(nanmean(fluor_cat_hemidiff(:,use_t_align(curr_align,:),curr_area,curr_align),2),[],1);
            
            clear g_act
            g_act = GLM(D).setModel(use_model).fitCV(cvfold);
            pL = g_act.p_hat(:,1);
            pR = g_act.p_hat(:,2);
            likelihood = pL.*(g_act.data.response==1) + pR.*(g_act.data.response==2);
            
            loglik_bpt_fluor(curr_area,curr_align) = nanmean(log2(likelihood));
            
        end
    end
    loglik_increase_fluor = loglik_bpt_fluor - loglik_bpt_stim;
    
    loglik_bpt_mua = nan(n_depths,2);
    for curr_align = 1:2
        for curr_area = 1:n_depths
            
            D.neur = reshape(nanmean(mua_cat(:,use_t_align(curr_align,:),curr_area,curr_align),2),[],1);
            
            clear g_act
            g_act = GLM(D).setModel(use_model).fitCV(cvfold);
            pL = g_act.p_hat(:,1);
            pR = g_act.p_hat(:,2);
            likelihood = pL.*(g_act.data.response==1) + pR.*(g_act.data.response==2);
            
            loglik_bpt_mua(curr_area,curr_align) = nanmean(log2(likelihood));
            
        end
    end
    loglik_increase_mua = loglik_bpt_mua - loglik_bpt_stim;
    
    % Fit stim + activity (ZERO CONTRASTS)
    use_model = 'AP_test_neur_nostim';
    
    loglik_bpt_fluor_zerocontrast = nan(n_rois,2);
    loglik_bpt_guess = nan(n_rois,2);
    for curr_align = 1:2
        for curr_area = 1:n_rois
            
            D_zero.neur = reshape(nanmean(fluor_cat_hemidiff(zero_contrasts,use_t_align(curr_align,:),curr_area,curr_align),2),[],1);
            
            clear g_act
            g_act = GLM(D_zero).setModel(use_model).fitCV(cvfold);
            pL = g_act.p_hat(:,1);
            pR = g_act.p_hat(:,2);
            likelihood = pL.*(g_act.data.response==1) + pR.*(g_act.data.response==2);
            
            loglik_bpt_fluor_zerocontrast(curr_area,curr_align) = nanmean(log2(likelihood));
            loglik_bpt_guess(curr_area,curr_align) = g_act.guess_bpt;
            
        end
    end
    loglik_increase_fluor_zerocontrast = bsxfun(@minus,loglik_bpt_fluor_zerocontrast,nanmean(loglik_bpt_guess,1));
    
    loglik_bpt_mua_zerocontrast = nan(n_depths,2);
    loglik_bpt_guess = nan(n_depths,2);
    for curr_align = 1:2
        for curr_area = 1:n_depths
            
            D_zero.neur = reshape(nanmean(mua_cat(zero_contrasts,use_t_align(curr_align,:),curr_area,curr_align),2),[],1);
            
            clear g_act
            g_act = GLM(D_zero).setModel(use_model).fitCV(cvfold);
            pL = g_act.p_hat(:,1);
            pR = g_act.p_hat(:,2);
            likelihood = pL.*(g_act.data.response==1) + pR.*(g_act.data.response==2);
            
            loglik_bpt_mua_zerocontrast(curr_area,curr_align) = nanmean(log2(likelihood));
            loglik_bpt_guess(curr_area,curr_align) = g_act.guess_bpt;
            
        end
    end
    loglik_increase_mua_zerocontrast = bsxfun(@minus,loglik_bpt_mua_zerocontrast,nanmean(loglik_bpt_guess,1));
    
    % Fit stim + activity (all L/R ROIs together)
    use_model = 'AP_test_roi_stimoffset';
    roi_params = nan(n_rois,2);
    for curr_align = 1:2
        
        D.offset_ZL = g_stim_all.ZL(behavParameterFit, g_stim_all.Zinput(g_stim_all.data));
        D.neur = squeeze(nanmean(fluor_cat(:,use_t_align(curr_align,:),:,curr_align),2));
        
        clear g_act
        g_act = GLM(D).setModel(use_model).fit;
        roi_params(:,curr_align) = g_act.parameterFits(2:end);
        
    end
   
    % Store
    batch_vars.loglik_increase_fluor(:,:,curr_animal) = loglik_increase_fluor;
    batch_vars.loglik_increase_mua(:,:,curr_animal) = loglik_increase_mua;
    batch_vars.loglik_increase_fluor_zerocontrast(:,:,curr_animal) = loglik_increase_fluor_zerocontrast;
    batch_vars.loglik_increase_mua_zerocontrast(:,:,curr_animal) = loglik_increase_mua_zerocontrast;   
    batch_vars.roi_params(:,:,curr_animal) = roi_params;
 
    % Clear out for next animal
    disp(['Finished ' animals{curr_animal}]);
    clearvars -except animals protocol batch_vars curr_animal
    
end

save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
save_fn = ['activity_choice_logistic_regression_meantime'];
save([save_path filesep save_fn],'batch_vars');
disp('Finished batch');

%% Batch logistic regression from activity (all time) - concatenate within animal
% MUA normalization could use some work

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};
protocol = 'vanillaChoiceworld';
batch_vars = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;   
    
    fluor_all = cell(length(experiments),1);
    mua_all = cell(length(experiments),1);
    D_all = cell(length(experiments),1);
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        AP_load_experiment;
        
        % Prepare fluorescence
        % (load widefield ROIs)
        wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
        load(wf_roi_fn);
        n_rois = numel(wf_roi);
        
        aUdf = single(AP_align_widefield(animal,day,Udf));
        roi_trace = AP_svd_roi(aUdf,fVdf,[],[],cat(3,wf_roi.mask));
        
        % (get ddf)
        roi_traces_derivative = diff(roi_trace,[],2);
        roi_traces_derivative(roi_traces_derivative < 0) = 0;
        
        % Prepare MUA
        % (group striatum depths)
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        depth_group = discretize(spikeDepths,depth_group_edges);
        
        % Get event-aligned activity
        raster_window = [-0.5,1];
        upsample_factor = 3;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        event_aligned_ddf = nan(length(stimOn_times),length(t),n_rois,2);
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths,2);
        for curr_align = 1:2
            switch curr_align
                case 1
                    use_align = stimOn_times;
                    use_align(isnan(use_align)) = 0;
                case 2
                    use_align = wheel_move_time';
                    use_align(isnan(use_align)) = 0;
            end
            
            t_peri_event = bsxfun(@plus,use_align,t);
            
            % Fluorescence
            event_aligned_ddf(:,:,:,curr_align) = ...
                interp1(conv2(frame_t,[1,1]/2,'valid'),roi_traces_derivative',t_peri_event);
            
            % MUA
            t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
            for curr_depth = 1:n_depths
                curr_spikes = spike_times_timeline(depth_group == curr_depth);
                event_aligned_mua(:,:,curr_depth,curr_align) = cell2mat(arrayfun(@(x) ...
                    histcounts(curr_spikes,t_bins(x,:)), ...
                    [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
            end
            
        end
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5 & ...
            stim_to_move < 0.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        % D.response = (trial_choice(use_trials)'+1)/2+1;
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);
        
        % Store activity and stim
        fluor_all{curr_day} = event_aligned_ddf(use_trials,:,:,:);
        mua_all{curr_day} = event_aligned_mua(use_trials,:,:,:);
        D_all{curr_day} = D;
        
    end
    
    % Get rid of everything that's not the data
    clearvars -except animals protocol batch_vars curr_animal fluor_all mua_all D_all
    
    % Get time
    framerate = 35.2;
    raster_window = [-0.5,1];
    upsample_factor = 3;
    raster_sample_rate = 1/(framerate*upsample_factor);
    t = raster_window(1):raster_sample_rate:raster_window(2);
    
    % Get widefield ROIs
    wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
    load(wf_roi_fn);
    n_rois = numel(wf_roi);
    
    n_depths = 6;
    
    cvfold = 20;
    
    % Concatenate activity (and normalize MUA)
    % (this normalization doesn't look very good right now...)
    fluor_cat = cat(1,fluor_all{:});
    
    mua_cat_raw = cat(1,mua_all{:});
    t_baseline = t < 0;
    softnorm = 5;
    mua_baseline = nanmean(mua_cat_raw(:,t_baseline,:,1),2);
    mua_cat = bsxfun(@rdivide,mua_cat_raw,mua_baseline + softnorm);
    
    % Get L-R activity
    fluor_cat_hemidiff = cat(3,fluor_cat(:,:,1:n_rois/2,:),fluor_cat(:,:,1:n_rois/2,:) - ...
        fluor_cat(:,:,n_rois/2+1:end,:));
    
    % Concatenate behavioural data
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all,'uni',false));
    D.response = cell2mat(cellfun(@(x) x.response,D_all,'uni',false));
    D.repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all,'uni',false));
    
    % Get zero-contrasts as a subset
    zero_contrasts = D.stimulus(:,1) == 0 & D.stimulus(:,2) == 0;
    D_zero = struct;
    D_zero.stimulus = D.stimulus(zero_contrasts,:);
    D_zero.response = D.response(zero_contrasts);
    D_zero.repeatNum = D.repeatNum(zero_contrasts);
    
    % Fit stim all
    use_model = 'AP_test_stim';
    g_stim_all = GLM(D).setModel(use_model).fit;
    behavParameterFit = g_stim_all.parameterFits;
    
    D.offset_ZL = g_stim_all.ZL(behavParameterFit, g_stim_all.Zinput(g_stim_all.data));
    
    % Fit stim cross-validated
    use_model = 'AP_test_stim';
    
    g_stim = GLM(D).setModel(use_model).fitCV(cvfold);
    pL = g_stim.p_hat(:,1);    pR = g_stim.p_hat(:,2);
    likelihood = pL.*(g_stim.data.response == 1) + pR.*(g_stim.data.response == 2);
    
    loglik_bpt_stim = nanmean(log2(likelihood));
    
    % Fit stim + activity
    use_model = 'AP_test_neur_stimoffset';
    
    loglik_bpt_fluor = nan(length(t),n_rois,2);
    for curr_align = 1:2
        for curr_area = 1:n_rois
            for curr_t = 1:length(t)
                
                D.neur = reshape(fluor_cat_hemidiff(:,curr_t,curr_area,curr_align),[],1);
                
                clear g_act
                g_act = GLM(D).setModel(use_model).fitCV(cvfold);
                pL = g_act.p_hat(:,1);
                pR = g_act.p_hat(:,2);
                likelihood = pL.*(g_act.data.response==1) + pR.*(g_act.data.response==2);
                
                loglik_bpt_fluor(curr_t,curr_area,curr_align) = nanmean(log2(likelihood));
                
            end
        end
    end
    loglik_increase_fluor = loglik_bpt_fluor - loglik_bpt_stim;
    
    loglik_bpt_mua = nan(length(t),n_depths,2);
    for curr_align = 1:2
        for curr_area = 1:n_depths
            for curr_t = 1:length(t)
                
                D.neur = reshape(mua_cat(:,curr_t,curr_area,curr_align),[],1);
                
                clear g_act
                g_act = GLM(D).setModel(use_model).fitCV(cvfold);
                pL = g_act.p_hat(:,1);
                pR = g_act.p_hat(:,2);
                likelihood = pL.*(g_act.data.response==1) + pR.*(g_act.data.response==2);
                
                loglik_bpt_mua(curr_t,curr_area,curr_align) = nanmean(log2(likelihood));
                
            end
        end
    end
    loglik_increase_mua = loglik_bpt_mua - loglik_bpt_stim;
    
    % Fit stim + activity (ZERO CONTRASTS)
    use_model = 'AP_test_neur_nostim';
    
    loglik_bpt_fluor_zerocontrast = nan(length(t),n_rois,2);
    loglik_bpt_guess = nan(length(t),n_rois,2);
    for curr_align = 1:2
        for curr_area = 1:n_rois
            for curr_t = 1:length(t)
                
                D_zero.neur = reshape(fluor_cat_hemidiff(zero_contrasts,curr_t,curr_area,curr_align),[],1);
                
                clear g_act
                g_act = GLM(D_zero).setModel(use_model).fitCV(cvfold);
                pL = g_act.p_hat(:,1);
                pR = g_act.p_hat(:,2);
                likelihood = pL.*(g_act.data.response==1) + pR.*(g_act.data.response==2);
                
                loglik_bpt_fluor_zerocontrast(curr_t,curr_area,curr_align) = nanmean(log2(likelihood));
                loglik_bpt_guess(curr_t,curr_area,curr_align) = g_act.guess_bpt;
                
            end
        end
    end
    loglik_increase_fluor_zerocontrast = bsxfun(@minus,loglik_bpt_fluor_zerocontrast,nanmean(loglik_bpt_guess,2));
    
    loglik_bpt_mua_zerocontrast = nan(length(t),n_depths,2);
    loglik_bpt_guess = nan(length(t),n_depths,2);
    for curr_align = 1:2
        for curr_area = 1:n_depths
            for curr_t = 1:length(t)
                
                D_zero.neur = reshape(mua_cat(zero_contrasts,curr_t,curr_area,curr_align),[],1);
                
                clear g_act
                g_act = GLM(D_zero).setModel(use_model).fitCV(cvfold);
                pL = g_act.p_hat(:,1);
                pR = g_act.p_hat(:,2);
                likelihood = pL.*(g_act.data.response==1) + pR.*(g_act.data.response==2);
                
                loglik_bpt_mua_zerocontrast(curr_t,curr_area,curr_align) = nanmean(log2(likelihood));
                loglik_bpt_guess(curr_t,curr_area,curr_align) = g_act.guess_bpt;
                
            end
        end
    end
    loglik_increase_mua_zerocontrast = bsxfun(@minus,loglik_bpt_mua_zerocontrast,nanmean(loglik_bpt_guess,2));
    
    % Fit stim + activity (all L/R ROIs together)
    use_model = 'AP_test_roi_stimoffset';
    roi_params = nan(length(t),n_rois,2);
    for curr_align = 1:2
        for curr_t = 1:length(t)
                       
            D.neur = squeeze(fluor_cat(:,curr_t,:,curr_align));
            
            clear g_act
            g_act = GLM(D).setModel(use_model).fit;
            roi_params(curr_t,:,curr_align) = g_act.parameterFits(2:end);
            
        end
    end
   
    % Store
    batch_vars.loglik_increase_fluor(:,:,:,curr_animal) = loglik_increase_fluor;
    batch_vars.loglik_increase_mua(:,:,:,curr_animal) = loglik_increase_mua;
    batch_vars.loglik_increase_fluor_zerocontrast(:,:,:,curr_animal) = loglik_increase_fluor_zerocontrast;
    batch_vars.loglik_increase_mua_zerocontrast(:,:,:,curr_animal) = loglik_increase_mua_zerocontrast;   
    batch_vars.roi_params(:,:,:,curr_animal) = roi_params;
    
    % Clear out for next animal
    animal = animals{curr_animal};
    disp(['Finished ' animal]);
    clearvars -except animals protocol batch_vars curr_animal
    
end

save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
save_fn = ['activity_choice_logistic_regression_alltime' ];
save([save_path filesep save_fn],'batch_vars');

disp('Finished batch');

%% Batch load and save activity from all chosen trials for all animals

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

fluor_all = cell(length(animals),1);
mua_all = cell(length(animals),1);
D_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;
    
    fluor_all{curr_animal} = cell(length(experiments),1);
    mua_all{curr_animal} = cell(length(experiments),1);
    D_all{curr_animal} = cell(length(experiments),1);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments);
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        verbose = false; AP_load_experiment;
        
        % Prepare fluorescence
        % (load widefield ROIs)
        wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
        load(wf_roi_fn);
        n_rois = numel(wf_roi);
        
        % OLD DDFF
        %     aUdf = single(AP_align_widefield(animal,day,Udf));
        %     roi_traces_df = AP_svd_roi(aUdf,fVdf,[],[],cat(3,wf_roi.mask));
        %
        %     roi_traces_derivative = diff(roi_traces_df,[],2);
        %     roi_traces_derivative(roi_traces_derivative < 0) = 0;
        %     frame_t_derivative = conv2(frame_t,[1,1]/2,'valid');
        
        % NEW DFF
        aU = single(AP_align_widefield(animal,day,U));
        roi_traces_f = AP_svd_roi(aU,fV,[],[],cat(3,wf_roi.mask));
        
        avg_im_aligned = AP_align_widefield(animal,day,avg_im);
        f_avg = cellfun(@(x) nanmedian(avg_im_aligned(x)),{wf_roi.mask})';
        
        baseline_window = [-0.5,0];
        t_baseline = baseline_window(1):1/(framerate):baseline_window(2);
        t_peri_baseline = bsxfun(@plus,stimOn_times,t_baseline);
        roi_f_baseline = ...
            squeeze(nanmedian(nanmedian(interp1(frame_t,roi_traces_f',t_peri_baseline),2),1));
        
        roi_traces_df = bsxfun(@minus,roi_traces_f,roi_f_baseline);
        roi_traces_dff = bsxfun(@rdivide,roi_traces_df,f_avg);
                
        % Just fluorescence, no derivative
        roi_traces_derivative = roi_traces_dff;
        frame_t_derivative = frame_t;
        
%         % Derivative via diff
%         roi_traces_derivative = diff(roi_traces_dff,[],2);
%         roi_traces_derivative(roi_traces_derivative < 0) = 0;
%         frame_t_derivative = conv2(frame_t,[1,1]/2,'valid');
        
%         % Derivative via low-pass filter, then derivative (from Nick)       
%         Nf = 50;
%         Fpass = 8.5;
%         Fstop = 10;
%         Fs = 1/mean(diff(frame_t)); %sampling frequency
%         
%         d = designfilt('differentiatorfir','FilterOrder',Nf, ...
%             'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
%             'SampleRate',Fs);
%         roi_traces_derivative = filter(d,roi_traces_dff')';
%         delay = mean(grpdelay(d));
%         roi_traces_derivative = circshift(roi_traces_derivative,[0,-delay]);
%         roi_traces_derivative(roi_traces_derivative < 0) = 0;
%         
%         frame_t_derivative = frame_t;
                
        % Prepare MUA
        % (group striatum depths)
        n_depths = 6;
        depth_group_edges = round(linspace(str_depth(1),str_depth(2),n_depths+1));
        depth_group = discretize(spikeDepths,depth_group_edges);
        
        % Get event-aligned activity
        raster_window = [-0.5,1];
        upsample_factor = 3;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        event_aligned_ddf = nan(length(stimOn_times),length(t),n_rois,2);
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths,2);
        for curr_align = 1:2
            switch curr_align
                case 1
                    use_align = stimOn_times;
                    use_align(isnan(use_align)) = 0;
                case 2
                    use_align = wheel_move_time';
                    use_align(isnan(use_align)) = 0;
            end
            
            t_peri_event = bsxfun(@plus,use_align,t);
            
            % Fluorescence
            event_aligned_ddf(:,:,:,curr_align) = ...
                interp1(frame_t_derivative,roi_traces_derivative',t_peri_event);
            
            % MUA
            t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
            for curr_depth = 1:n_depths
                curr_spikes = spike_times_timeline(depth_group == curr_depth);
                event_aligned_mua(:,:,curr_depth,curr_align) = cell2mat(arrayfun(@(x) ...
                    histcounts(curr_spikes,t_bins(x,:)), ...
                    [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
            end
            
        end
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5 & ...
            stim_to_move > 0.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        % D.response = (trial_choice(use_trials)'+1)/2+1;
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);
        
        % Store activity and stim
        fluor_all{curr_animal}{curr_day} = event_aligned_ddf(use_trials,:,:,:);
        mua_all{curr_animal}{curr_day} = event_aligned_mua(use_trials,:,:,:);
        D_all{curr_animal}{curr_day} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except animals curr_animal fluor_all mua_all D_all
    
end
clearvars -except fluor_all mua_all D_all
disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
save_fn = 'all_trial_activity_f_latemove.mat';
save([save_path filesep save_fn]);

%% Batch modeling on saved day-concatenated activity

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = 'all_trial_activity_diff_earlymove.mat';
load([data_path filesep data_fn]);

% Get time
framerate = 35.2;
raster_window = [-0.5,1];
upsample_factor = 3;
raster_sample_rate = 1/(framerate*upsample_factor);
t = raster_window(1):raster_sample_rate:raster_window(2);

use_t_stim = t > 0.05 & t < 0.15;
use_t_move = t > -0.15 & t < -0.02;
use_t_align = [use_t_stim;use_t_move];

% Get widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

n_depths = 6;

n_animals = length(D_all);

loglik_increase_fluor = nan(length(t),n_rois,2,n_animals);
loglik_increase_mua = nan(length(t),n_depths,2,n_animals);
for curr_animal = 1:n_animals
    
    % Concatenate fluorescence
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    smooth_size = 8;
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    
    mua_cat_raw_smoothed = convn(mua_cat_raw,smWin,'same');
    
    t_baseline = t < 0;
    softnorm = 10;
    mua_baseline = nanmean(mua_cat_raw_smoothed(:,t_baseline,:,1),2);
    mua_cat = bsxfun(@rdivide,mua_cat_raw_smoothed,mua_baseline + softnorm)-1;
    
    % Normalize all activity by std
    fluor_cat = bsxfun(@rdivide,fluor_cat,permute(std(reshape(permute(fluor_cat,[1,2,4,3]),[],n_rois),[],1),[1,3,2,4]));
    mua_cat = bsxfun(@rdivide,mua_cat,permute(std(reshape(permute(mua_cat,[1,2,4,3]),[],n_depths),[],1),[1,3,2,4]));
    
    % Get L-R activity
    fluor_cat_hemidiff = cat(3,fluor_cat(:,:,1:n_rois/2,:),fluor_cat(:,:,1:n_rois/2,:) - ...
        fluor_cat(:,:,n_rois/2+1:end,:));
    
    % Concatenate behavioural data
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
    D.response = cell2mat(cellfun(@(x) x.response,D_all{curr_animal},'uni',false));
    D.repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all{curr_animal},'uni',false));
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    D.day = trial_day;
    
    max_contrast = max(D.stimulus,[],2);
    use_trials = max_contrast >= 0 & max_contrast < 1;   
    
    %%% Modeling
    cvfold = 20;
    
    warning off;
    % Fit stim all
    use_model = 'AP_test_stim';
    g_stim_all = GLM(D).setModel(use_model).fit;
    behavParameterFit = g_stim_all.parameterFits;
    
    D.offset_ZL = g_stim_all.ZL(behavParameterFit, g_stim_all.Zinput(g_stim_all.data));
    
    % Fit stim cross-validated
    use_model = 'AP_test_stim';
    
    D_use = structfun(@(x) x(use_trials,:),D,'uni',false);
    
    g_stim = GLM(D_use).setModel(use_model).fitCV(cvfold);
    pL = g_stim.p_hat(:,1);    pR = g_stim.p_hat(:,2);
    likelihood = pL.*(g_stim.data.response == 1) + pR.*(g_stim.data.response == 2);
    
    loglik_bpt_stim = nanmean(log2(likelihood));
    
    % Fit stim + activity (all time)
    use_model = 'AP_test_neur_stimoffset';
    
    loglik_bpt_fluor = nan(length(t),n_rois,2);
    for curr_align = 1:2
        for curr_area = 1:n_rois
            for curr_t = 1:length(t)
                
                % Set the activity
                D.neur = reshape(fluor_cat_hemidiff(:,curr_t,curr_area,curr_align),[],1);
                
                % Pick subset of trials
                D_use = structfun(@(x) x(use_trials,:),D,'uni',false);
                
                clear g_act
                g_act = GLM(D_use).setModel(use_model).fitCV(cvfold);
                pL = g_act.p_hat(:,1);
                pR = g_act.p_hat(:,2);
                likelihood = pL.*(g_act.data.response==1) + pR.*(g_act.data.response==2);
                
                loglik_bpt_fluor(curr_t,curr_area,curr_align) = nanmean(log2(likelihood));
            end
        end
    end
    loglik_increase_fluor(:,:,:,curr_animal) = loglik_bpt_fluor - loglik_bpt_stim;
    
    loglik_bpt_mua = nan(length(t),n_depths,2);
    for curr_align = 1:2
        for curr_area = 1:n_depths
            for curr_t = 1:length(t)
                
                % Set the activity
                D.neur = reshape(mua_cat(:,curr_t,curr_area,curr_align),[],1);

                % Pick subset of trials
                D_use = structfun(@(x) x(use_trials,:),D,'uni',false);
                
                clear g_act
                g_act = GLM(D_use).setModel(use_model).fitCV(cvfold);
                pL = g_act.p_hat(:,1);
                pR = g_act.p_hat(:,2);
                likelihood = pL.*(g_act.data.response==1) + pR.*(g_act.data.response==2);
                
                loglik_bpt_mua(curr_t,curr_area,curr_align) = nanmean(log2(likelihood));
            end
        end
    end
    loglik_increase_mua(:,:,:,curr_animal) = loglik_bpt_mua - loglik_bpt_stim;
    
    warning on;
    
    AP_print_progress_fraction(curr_animal,n_animals);
end

clearvars -except loglik_increase_fluor loglik_increase_mua
save_path = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld'];
save_fn = ['activity_sessioncat_logistic_regression_earlymove'];
save([save_path filesep save_fn])

disp('Finished');

