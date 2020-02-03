%% Testing for corticostriatal-specific imaging


%% Load and average stim-aligned image

animals = {'AP063','AP064','AP066','AP068'};
day = '2020-02-03';
experiment = 1;
verbose = true; 

im_stim_all = cell(length(animals),1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    AP_load_experiment;
    
    % Set options
    surround_window = [-0.5,3];
    baseline_window = [-0.1,0];
    
    surround_samplerate = 1/(framerate*1);
    surround_time = surround_window(1):surround_samplerate:surround_window(2);
    baseline_surround_time = baseline_window(1):surround_samplerate:baseline_window(2);
    
    % Average (time course) responses
    use_vs = 1:size(U,3);
    
    conditions = unique(stimIDs);
    im_stim = nan(size(U,1),size(U,2),length(surround_time),length(conditions));
    for curr_condition_idx = 1:length(conditions)
        curr_condition = conditions(curr_condition_idx);
        
        use_stims = find(stimIDs == curr_condition);
        use_stimOn_times = stimOn_times(use_stims(2:end));
        use_stimOn_times([1,end]) = [];
        
        stim_surround_times = bsxfun(@plus, use_stimOn_times(:), surround_time);
        stim_baseline_surround_times = bsxfun(@plus, use_stimOn_times(:), baseline_surround_time);
        
        peri_stim_v = permute(interp1(frame_t,fVdf',stim_surround_times),[3,2,1]);
        baseline_v = permute(nanmean(interp1(frame_t,fVdf',stim_baseline_surround_times),2),[3,2,1]);
        
        stim_v_mean = nanmean(peri_stim_v - baseline_v,3);
        
        im_stim(:,:,:,curr_condition_idx) = svdFrameReconstruct(Udf(:,:,use_vs), ...
            stim_v_mean(use_vs,:));
    end
    
    im_stim_all{curr_animal} = im_stim;
    
end

im_stim_cat = cat(5,im_stim_all{:});
AP_image_scroll(nanmean(im_stim_cat(:,:,:,2,:),5));
axis image










