%% This is to load and process animal 65, voltage imaging by DS

%% StackSet: Load image from single trial

animal = '65';
day = '20151012';
seq = 5; % this is experiment within day

% Stackset class written by DS and MC
addpath('\\zserver\Code\stacks'); %stackset

% Path for raw data
info_data_path = '\\zserver3.cortexlab.net\Data\Stacks';

% If day entered as number, convert to string
if isnumeric(day)
    day = num2str(day);
end

% Find filenames for behavior/input
day_dash = datestr(datenum(day,'yyyymmdd'),'yyyy-mm-dd');
timeline_filename = dat.expFilePath(animal,day_dash,seq,'timeline','master');
parameters_filename = dat.expFilePath(animal,day_dash,seq,'parameters','master');
block_filename = dat.expFilePath(animal,day_dash,seq,'block','master');

% Load behavior/input
load(timeline_filename);
load(parameters_filename);
load(block_filename);

% Make protocol structure required for StackSet loading
protocol.animal = animal;
protocol.iseries = str2num(day);
protocol.iexp = seq;
protocol.nstim = 1;
protocol.seqnums = (1:block.numCompletedTrials)-1; % 1st trial is not imaged

% Load StackSet of single trial
resize_fac = 0.25; % hard coded for now - this must be in a header somewhere
iStim = 1; % "averaged repeat for stimulus iStim" ??
tr = 20; % trial to load
suffix = 'ar gd'; % naming convention from DS
cam = tools.Imager('PCO', [], '_ratio', 1.6, 4); % camera settings
singleStack = StackSet.LoadStacks(info_data_path,protocol,resize_fac,iStim,tr,cam,'suffix',suffix);


%% StackSet: Load many trials and concatenate

% This uses a bunch of lab code
addpath(genpath('\\zserver.cortexlab.net\Code\stacks')); % stackset
addpath(genpath('\\zserver.cortexlab.net\Code\Rigging\main')); % CB stuff
addpath(genpath('\\zserver.cortexlab.net\Code\Rigging\cb-tools')) % more CB stuff

animal = 'M160128_SF';
day = '20160526';
seq = 2; % this is experiment within day

% Path for raw data
info_data_path = '\\zserver3.cortexlab.net\Data\Stacks';

% If day entered as number, convert to string
if isnumeric(day)
    day = num2str(day);
end

% Find filenames for behavior/input
day_dash = datestr(datenum(day,'yyyymmdd'),'yyyy-mm-dd');
timeline_filename = dat.expFilePath(animal,day_dash,seq,'timeline','master');
parameters_filename = dat.expFilePath(animal,day_dash,seq,'parameters','master');
block_filename = dat.expFilePath(animal,day_dash,seq,'block','master');
protocol_filename = get_cortexlab_filename(animal,day_dash,seq,'protocol','8digit');

% Load behavior/input
try
    load(timeline_filename);
catch me
    disp('No timeline')
end
try
    load(parameters_filename);
catch me
    disp('No parameters')
end
try
    load(block_filename);
catch me
    disp('No block')
end
try
    load(protocol_filename);
catch me
    disp('No protocol')
end

% Make protocol structure required for StackSet loading
protocol.animal = animal;
protocol.iseries = str2num(day);
protocol.iexp = seq;
protocol.nstim = 1;
protocol.seqnums = (1:block.numCompletedTrials)-1; % 1st trial is not imaged

% Load StackSet of selected trials (first trial doesn't work?)
load_trials = 2:40;
resize_fac = 0.25; % hard coded for now - this must be in a header somewhere
iStim = 1; % "averaged repeat for stimulus iStim" ??

% Load 2nd trial to get size
tr = 8;
suffix = 'ar gd'; % naming convention from DS
suffix = ''; % if not DS
cam = tools.Imager('PCO', [], '_ratio', 1.6, 4); % camera settings (voltage)
cam = tools.Imager('PCO', [], '', 1.6, 4); % camera settings (gcamp)

singleStack = StackSet.LoadStacks(info_data_path,protocol,resize_fac,iStim,tr,cam,'suffix',suffix);
[im_y,im_x,numframes] = size(singleStack.Values);

% Load in all selected trials
im_concat = [];
for curr_trial_idx = 1:length(load_trials)
    curr_trial = load_trials(curr_trial_idx);
    suffix = 'ar gd'; % naming convention from DS
    suffix = ''; % if not DS
    cam = tools.Imager('PCO', [], '_ratio', 1.6, 4); % camera settings
    singleStack = StackSet.LoadStacks(info_data_path,protocol,resize_fac,iStim,curr_trial,cam,'suffix',suffix);
    
    im_concat = [im_concat,reshape(singleStack.Values,im_x*im_y,[])];
    disp(['Loaded ' num2str(curr_trial)]);
end

%% StackSet: load others

% This uses a bunch of lab code
addpath(genpath('\\zserver.cortexlab.net\Code\stacks')); % stackset
addpath(genpath('\\zserver.cortexlab.net\Code\Rigging\main')); % CB stuff
addpath(genpath('\\zserver.cortexlab.net\Code\Rigging\cb-tools')) % more CB stuff

% Animal info
Exps.animal = 'M150723_SF_cam2';
Exps.iseries = '001';
Exps.iexp = '002';
suffix = '';

stack_dir = '\\zserver3.cortexlab.net\Data\Stacks';
stim_dir = 'Resize 50 Stim 001 Repeat All';

stackDir = [stack_dir filesep Exps.animal filesep Exps.iseries filesep Exps.iexp filesep stim_dir];

% Get stack filenames
stack_dir = dir([stackDir filesep '*.bin']);
stack_filenames = {stack_dir.name};

n_end_frames = 3; % frames to cut off at the end
% use_stacks = length(stack_filenames)
use_stacks = 10;

clear all_stacks;
all_stacks(use_stacks) = StackSet;
for curr_stack = 1:use_stacks
    all_stacks(curr_stack) = tools.LoadMyStacks(stackDir, stack_filenames{curr_stack}(1:end-4));
    all_stacks(curr_stack).TimeVec(end-n_end_frames+1:end) = [];
    all_stacks(curr_stack).Values(:,:,end-n_end_frames+1:end) = [];
end

t_all = {all_stacks.TimeVec};
t_ends = num2cell([0,cumsum(cellfun(@(x) x(end),t_all(1:end-1)))]);
t_cat = cell2mat(cellfun(@(t_all,t_ends) t_all+t_ends,t_all,t_ends,'uni',false));
stack_cat = cat(3,all_stacks.Values);

mean_im = nanmean(stack_cat,3);
dff = bsxfun(@rdivide,stack_cat,mean_im) - 1;

figure;imagesc(mean_im);
title(Exps.animal);colormap(gray);
set(gca,'YDir','normal');
axis off;

AP_image_scroll(dff,t_cat);
set(gca,'YDir','normal');

%% Get average movie around stimulus of all selected trials

frames_back = 50;
frames_forward = 100;

animal = '65';
days = {'20151012','20151013','20151014'};

% Stackset class written by DS and MC
addpath('\\zserver.cortexlab.net\Code\stacks'); %stackset

% Path for raw data
info_data_path = '\\zserver3.cortexlab.net\Data\Stacks';

% Loop through days, get average aligned image

im_y = 120; im_x = 92; % just hardcode for now
% Sum image: y,x,frames,stim side,contrast,correct
im_sum = zeros(im_y,im_x,frames_back+frames_forward+1,2,4,2);
% Trials per image: stim side,contrast,correct
n_trials_condition = zeros(2,4,2);

for curr_day = 1:length(days)
    
    day = days{curr_day};
    
    % Get experiment number within day (should only be one folder)
    seq_dir = dir([info_data_path filesep animal '_ratio' filesep day]);
    if length(seq_dir) > 3
        error('Day folder has more than one entry')
    end
    seq = str2num(seq_dir(3).name);
    
    
    % If day entered as number, convert to string
    if isnumeric(day)
        day = num2str(day);
    end
    
    % Find filenames for behavior/input
    day_dash = datestr(datenum(day,'yyyymmdd'),'yyyy-mm-dd');
    timeline_filename = dat.expFilePath(animal,day_dash,seq,'timeline','master');
    parameters_filename = dat.expFilePath(animal,day_dash,seq,'parameters','master');
    block_filename = dat.expFilePath(animal,day_dash,seq,'block','master');
    
    % Load behavior/input
    load(timeline_filename);
    load(parameters_filename);
    load(block_filename);
    
    timeline_sample_rate = Timeline.hw.daqSampleRate;
    n_trials = block.numCompletedTrials;
    
    % Imaging epoch signal from timeline
    timeline_syncEcho_idx = strcmp({Timeline.hw.inputs.name}, 'syncEcho');
    syncEcho_start_samples = find(Timeline.rawDAQData(1:end-1,timeline_syncEcho_idx) <= 2 & ...
        Timeline.rawDAQData(2:end,timeline_syncEcho_idx) > 2);
    syncEcho_stop_samples = find(Timeline.rawDAQData(1:end-1,timeline_syncEcho_idx) > 2 & ...
        Timeline.rawDAQData(2:end,timeline_syncEcho_idx) <= 2);
    
    syncEcho_start_time = syncEcho_start_samples./timeline_sample_rate;
    syncEcho_stop_time = syncEcho_stop_samples./timeline_sample_rate;
    
    % Camera signal from timeline
    timeline_cam2_idx = strcmp({Timeline.hw.inputs.name}, 'cam2');
    cam2_samples = find(Timeline.rawDAQData(1:end-1,timeline_cam2_idx) <= 2 & ...
        Timeline.rawDAQData(2:end,timeline_cam2_idx) > 2);
    cam2_time = cam2_samples./timeline_sample_rate;
    
    % Get the imaging epoch start/stop time for each trial (first trial not imaged)
    trial_start_time = [block.trial.trialStartedTime];
    imaging_epoch_start_time = syncEcho_start_time(arrayfun(@(x) find(syncEcho_start_time <= ...
        trial_start_time(x),1,'last'),2:n_trials));
    imaging_epoch_stop_time = syncEcho_stop_time(arrayfun(@(x) find(syncEcho_stop_time > ...
        imaging_epoch_start_time(x),1),1:n_trials-1));
    
    % Get the frame time indicies for each imaging epoch
    imaging_epoch_frame_time = arrayfun(@(x) cam2_time( ...
        cam2_time > imaging_epoch_start_time(x) & ...
        cam2_time < imaging_epoch_stop_time(x)),1:n_trials-1,'uni',false);
    
    % Photodiode signal from timeline
    photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
    photodiode = Timeline.rawDAQData(:,photodiode_idx);
    photodiode_flip = find((photodiode(1:end-1) <= 1.5 & photodiode(2:end) > 1.5) | ...
        (photodiode(1:end-1) >= 1.5 & photodiode(2:end) < 1.5))+1;
    photodiode_sec = photodiode_flip/timeline_sample_rate;
    
    % Get stimulus presentation times relative to imaging frames
    stim_times_block = [block.trial.stimulusCueStartedTime];
    stim_samples_photodiode = photodiode_flip(arrayfun(@(x) find(photodiode_sec >= ...
        stim_times_block(x),1),2:n_trials));
    stim_times_photodiode = stim_samples_photodiode./timeline_sample_rate;
    [~,stim_frames_photodiode] = arrayfun(@(x) ...
        min(abs(imaging_epoch_frame_time{x} - stim_times_photodiode(x))), ...
        1:n_trials-1);
    
    % Make protocol structure required for StackSet loading
    protocol.animal = animal;
    protocol.iseries = str2num(day);
    protocol.iexp = seq;
    protocol.nstim = 1;
    protocol.seqnums = (1:block.numCompletedTrials)-1; % 1st trial is not imaged
    
    % Load StackSet of correct trials across conditions
    contrasts = [0,0.12,0.25,0.5];
    
    % Load 2nd trial to get size
    tr = 2;
    resize_fac = 0.25; % hard coded for now - this must be in a header somewhere
    iStim = 1; % "averaged repeat for stimulus iStim" ??
    suffix = 'ar gd'; % naming convention from DS
    cam = tools.Imager('PCO', [], '_ratio', 1.6, 4); % camera settings
    singleStack = StackSet.LoadStacks(info_data_path,protocol,resize_fac,iStim,tr,cam,'suffix',suffix);
    [im_y,im_x,numframes] = size(singleStack.Values);
    %stim side,contrast,correct
    for curr_side = 1:2       
        for curr_contrast_idx = 1:length(contrasts)            
            curr_contrast = contrasts(curr_contrast_idx);           
            for curr_correct = 1:2
                
                switch curr_correct
                    case 1
                        curr_correct_vector = [block.trial(1:n_trials).feedbackType]' == 1;
                    case 2
                        curr_correct_vector = [block.trial(1:n_trials).feedbackType]' == -1;                       
                end
                
                trial_condition = cell2mat(cellfun(@(x) x.visCueContrast,{block.trial(1:n_trials).condition},'uni',false))';
                trial_correct = [block.trial(1:n_trials).feedbackType]' == 1;
                
                load_trials = find(trial_condition(:,curr_side) == curr_contrast ...
                    & trial_condition(:,setdiff(1:2,curr_side)) == 0 & curr_correct_vector);
                %load_trials = randperm(n_trials,length(load_trials)); % load random trials: sanity check
                load_trials(load_trials == 1) = []; % (first trial not imaged)
                
                % Load in and average all selected trials
                curr_im_sum = zeros(im_y,im_x,frames_back+frames_forward+1);
                sum_trials = 0;
                for curr_trial_idx = 1:length(load_trials)
                    curr_trial = load_trials(curr_trial_idx);
                    suffix = 'ar gd'; % naming convention from DS
                    cam = tools.Imager('PCO', [], '_ratio', 1.6, 4); % camera settings
                    try
                        singleStack = StackSet.LoadStacks(info_data_path,protocol,resize_fac,iStim,curr_trial,cam,'suffix',suffix);
                        curr_im = double(singleStack.Values);
                    catch me
                        continue
                    end
                    
                    curr_stim_frame = stim_frames_photodiode(curr_trial - 1);
                    if curr_stim_frame + frames_forward > size(curr_im,3) || ...
                            curr_stim_frame - frames_back < 0;
                        continue
                    end
                    
                    curr_im_sum = curr_im_sum + curr_im(:,:,curr_stim_frame - frames_back: ...
                        curr_stim_frame + frames_forward);
                    
                    sum_trials = sum_trials + 1;
                    disp(['Loaded ' num2str(curr_trial_idx) '/' num2str(length(load_trials))]);
                end
                
                % Condition order: stim side,contrast,correct
                n_trials_condition(curr_side,curr_contrast_idx,curr_correct) = ...
                    n_trials_condition(curr_side,curr_contrast_idx,curr_correct) + sum_trials;
                im_sum(:,:,:,curr_side,curr_contrast_idx,curr_correct) = ...
                    im_sum(:,:,:,curr_side,curr_contrast_idx,curr_correct) + curr_im_sum;
                
            end
        end
    end
    
end

% Get average across each condition
im_mean = nan(size(im_sum));
for curr_side = 1:2
    for curr_contrast = 1:4
        for curr_correct = 1:2
            im_mean(:,:,:,curr_side,curr_contrast,curr_correct) = ...
                im_sum(:,:,:,curr_side,curr_contrast,curr_correct)./ ...
                n_trials_condition(curr_side,curr_contrast,curr_correct);
        end
    end
end

% Denoise by PCA
[coeff,score,latent] = princomp(reshape(im_mean,size(im_mean,1)*size(im_mean,2),[])');
im_denoised = reshape(coeff(:,1:50)*score(:,1:50)',size(im_mean,1), ...
    size(im_mean,2),size(im_mean,3),size(im_mean,4),size(im_mean,5),size(im_mean,6));

% ROI positions (these were chosen beforehand)
roi = struct;

roi.x{1} = 57:64;
roi.y{1} = 83:90;
roi.area{1} = 'V1 (loc 1)';

roi.x{2} = 38:45;
roi.y{2} = 63:70;
roi.area{2} = 'RL';

roi.x{3} = 30:37;
roi.y{3} = 75:82;
roi.area{3} = 'AL';

roi.x{4} = 47:54;
roi.y{4} = 83:90;
roi.area{4} = 'V1 (loc 2)';

roi.x{5} = 62:69;
roi.y{5} = 77:84;
roi.area{5} = 'PM';

roi.x{6} = 57:64;
roi.y{6} = 83:90;
roi.area{6} = 'V1 (loc 1, RF match)';

roi.x{7} = 70:77;
roi.y{7} = 36:43;
roi.area{7} = 'Cg';

% Plot ROIs over average image
figure;
imagesc(nanmean(sum(sum(sum(im_sum,4),5),6),3)./sum(n_trials_condition(:)));
colormap(gray);
caxis([0.999,1.002])
for curr_roi = 1:length(roi.x)
    rectangle('Position',[roi.x{curr_roi}(1),roi.y{curr_roi}(1), ...
        roi.x{curr_roi}(end) - roi.x{curr_roi}(1), ...
        roi.y{curr_roi}(end) - roi.y{curr_roi}(1)],'EdgeColor','r', ...
        'LineWidth',2);
    text_offset_x = 2;
    if curr_roi == 4
        text_offset_y = 10;
    elseif curr_roi == 6
        text_offset_y = 5;
    else 
        text_offset_y = 2;
    end
    text(roi.x{curr_roi}(1)+text_offset_x,roi.y{curr_roi}(1)+text_offset_y, ...
        [roi.area{curr_roi} ' (' num2str(curr_roi) ')'],'color','y');
end
axis off;

% Condition order: stim side,contrast,correct

% For now: plot all contrasts, contralateral correct

im_avg = bsxfun(@rdivide,squeeze(im_sum(:,:,:,2,:,1)),permute(n_trials_condition(2,:,1),[4,5,1,2,3]));

% Save movie
im_avg_diff = im_avg(:,:,:,4)./im_avg(:,:,:,1);
plot_frames = frames_back-5:frames_back+51;
framerate = 10;
color_axis = [0.998,1.003];
figure_position =  [76,250,580,716];
savefile = 'C:\Users\Andy\Documents\CarandiniHarrisLab\65_analysis\65_50contrast-0contrast.avi';
 AP_movie2avi(im_avg_diff(:,:,plot_frames),framerate,color_axis,figure_position,savefile)

figure;
for curr_roi = 1:length(roi.x);
    % Plot timecourse
    subplot(length(roi.x),2,(curr_roi-1)*2+1); hold on;
    roi_timecourse = ...
        squeeze(nanmean(reshape(im_avg(roi.y{curr_roi},roi.x{curr_roi},:,:), ...
        length(roi.y{curr_roi})*length(roi.x{curr_roi}),size(im_avg,3),size(im_avg,4)),1));
    roi_timecourse_norm = bsxfun(@rdivide,roi_timecourse,nanmean(roi_timecourse(frames_back-20:frames_back-5,:),1))-1;
    set(gca,'ColorOrder',copper(size(roi_timecourse,2)));
    
    t = (plot_frames-frames_back+1)/50;
    
    plot(t,roi_timecourse_norm(plot_frames,:),'linewidth',2);
    xlim([min(t),max(t)]);
    ylim([min(reshape(roi_timecourse_norm(plot_frames,:),[],1)), ...
        max(reshape(roi_timecourse_norm(plot_frames,:),[],1))]);
    
    line([0,0],ylim,'color','r','linewidth',2);
    
    set(gca,'YTick',[]);
    ylabel(roi.area{curr_roi});
    xlabel('Time from stim onset (s)');
    
    % Plot mean activity as a function of contrast
    avg_times = frames_back+6:frames_back+11; % 100-200 ms, used by DS
    contrast_act = nanmean(roi_timecourse_norm(avg_times,:),1);
    subplot(length(roi.x),2,(curr_roi-1)*2+2);
    plot(contrasts,contrast_act,'k','linewidth',2);
    ylabel('Avg act');
    xlabel('Contrast');
    xlim([min(contrasts),max(contrasts)]);
    set(gca,'YTick',[]);

end

% Plot correct vs incorrect activity timecourse
figure;
for curr_roi = 1:7
    
    im_avg_correct = bsxfun(@rdivide,squeeze(im_sum(:,:,:,2,:,1)),permute(n_trials_condition(2,:,1),[4,5,1,2,3]));
    roi_timecourse_correct = ...
        squeeze(nanmean(reshape(im_avg_correct(roi.y{curr_roi},roi.x{curr_roi},:,:), ...
        length(roi.y{curr_roi})*length(roi.x{curr_roi}),size(im_avg_correct,3),size(im_avg_correct,4)),1));
    roi_timecourse_correct_norm = bsxfun(@rdivide,roi_timecourse_correct,nanmean(roi_timecourse_correct(frames_back-20:frames_back-5,:),1))-1;
    
    im_avg_incorrect = bsxfun(@rdivide,squeeze(im_sum(:,:,:,2,:,2)),permute(n_trials_condition(2,:,2),[4,5,1,2,3]));
    roi_timecourse_incorrect = ...
        squeeze(nanmean(reshape(im_avg_incorrect(roi.y{curr_roi},roi.x{curr_roi},:,:), ...
        length(roi.y{curr_roi})*length(roi.x{curr_roi}),size(im_avg_incorrect,3),size(im_avg_incorrect,4)),1));
    roi_timecourse_incorrect_norm = bsxfun(@rdivide,roi_timecourse_incorrect,nanmean(roi_timecourse_incorrect(frames_back-20:frames_back-5,:),1))-1;
    
    for curr_contrast = 1:4
        subplot(4,7,(curr_contrast-1)*7+curr_roi); hold on;
        plot(t,roi_timecourse_correct_norm(plot_frames,curr_contrast),'k','linewidth',2);
        plot(t,roi_timecourse_incorrect_norm(plot_frames,curr_contrast),'r','linewidth',2);
        xlim([min(t),max(t)]);
        
        if curr_contrast == 1
            title(roi.area{curr_roi})
        end
        if curr_roi == 1
            ylabel(['Contrast ' num2str(contrasts(curr_contrast))]);
        end
        if curr_roi == 7 && curr_contrast == 4;
            legend({'Correct','Incorrect'});
        end
    end
    
end
 
%% TO DO: wheel turning + motor cortex?
% keep in mind that often the take their arms off the wheel as it's
% spinning, so the wheel is definitely not a good proxy for forelimb
% movement

%% ====== EVERYTHING PAST THIS POINT IS GENERAL-PURPOSE =========

%% Align widefield to events

% Set options
surround_window = [-1,5];
baseline_surround_window = [-1,0];
framerate = 1./median(diff(frame_t));
surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);
baseline_surround_time = baseline_surround_window(1):surround_samplerate:baseline_surround_window(2);

% Gui for plotting responses
pixelTuningCurveViewerSVD(U,fV,frame_t,stim_onsets,stimIDs,surround_window);

% Average (time course) responses
im_stim = nan(size(U,1),size(U,2),length(surround_time),max(unique(stimIDs)));
for curr_condition = unique(stimIDs)'

    use_stims = find(stimIDs == curr_condition);
    use_stim_onsets = stim_onsets(use_stims(2:end));
    use_stim_onsets([1,end]) = [];
        
    stim_surround_times = bsxfun(@plus, use_stim_onsets(:), surround_time);
    peri_stim_v = permute(mean(interp1(frame_t,fV',stim_surround_times),1),[3,2,1]);
    
    im_stim(:,:,:,curr_condition) = svdFrameReconstruct(U,peri_stim_v);   
end

AP_image_scroll(im_stim,surround_time);

% Average (single frame) responses to stimuli
surround_window = [0.1,0.2];
framerate = 1./median(diff(frame_t));
surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);

peri_stim_v = nan(size(fV,1),max(stimIDs));
for curr_stim = 1:max(stimIDs)
    align_times = stim_onsets(stimIDs == curr_stim);
    align_times(align_times + surround_time(1) < frame_t(2) | ...
        align_times + surround_time(2) > frame_t(end)) = [];
    
    align_surround_times = bsxfun(@plus, align_times, surround_time);
    peri_stim_v(:,curr_stim) = nanmean(permute(nanmean(interp1(frame_t,fV',align_surround_times),1),[3,2,1]),2);
end

surround_im = svdFrameReconstruct(U,peri_stim_v);




%% Autocorrelation of each pixel
% (i.e. do medial areas have wider autocorrelation shapes?)

% Lifted from pixelCorrelationViewerSVD
Ur = reshape(U, size(U,1)*size(U,2),[]); % P x S
covV = cov(fV'); % S x S % this is the only one that takes some time really
varP = dot((Ur*covV)', Ur'); % 1 x P

pixelInd = sub2ind([size(U,1), size(U,2)], 250, 350);
covP = Ur(pixelInd,:)*covV*Ur'; % 1 x P
stdPxPy = varP(pixelInd).^0.5 * varP.^0.5; % 1 x P
covP = Ur(pixelInd,:)*covV*Ur'; % 1 x P
stdPxPy = varP(pixelInd).^0.5 * varP.^0.5; % 1 x P
corrMat = covP./stdPxPy; % 1 x P


% Autocorrelation of all V's;
maxlag = 100;
V_autocorr = nan(size(V,1),maxlag*2+1);
for i = 1:size(V,1)
    V_autocorr(i,:) = xcorr(fV(i,:),maxlag);
end

svd_xcorr = abs(Ur*V_autocorr);
svd_xcorr_norm = bsxfun(@minus,svd_xcorr,min(svd_xcorr,[],2));
fwhm = reshape(arrayfun(@(x) 2*((maxlag+1)-find(svd_xcorr_norm(x,:) > ...
    max(svd_xcorr_norm(x,:))/2,1)),1:size(svd_xcorr_norm,1)),size(U,1),size(U,2));
figure;imagesc(fwhm);colormap(gray);


U_px = squeeze(U(218,167,:));




% Choose ROI
h = figure;
imagesc(avg_im);
set(gca,'YDir','normal');
colormap(gray);
caxis([0 10000]);
roiMask = roipoly;
close(h);

% Get fluorescence across session in ROI
U_roi = reshape(U(repmat(roiMask,1,1,size(U,3))),[],size(U,3));
roi_trace = nanmean(U_roi*fV);


%% Get fluorescence trace of ROI

% % Choose ROI
% h = figure;
% imagesc(avg_im);
% set(gca,'YDir','reverse');
% colormap(gray);
% caxis([0 prctile(avg_im(:),90)]);
% roiMask = roipoly;
% close(h);
% 
% % Get fluorescence across session in ROI
% U_roi = reshape(U(repmat(roiMask,1,1,size(U,3))),[],size(U,3));
% roi_trace = nanmean(U_roi*fV);
% 
% figure;
% plot(frame_t,roi_trace,'k');
% xlabel('Time (s)')
% ylabel('ROI Fluorescence')

% Now have a function to do this:
roi_trace = AP_svd_roi(U,fV);
figure;plot(frame_t,roi_trace,'k');

%% Correlate fluorescence with trace

% Set the trace to use
use_trace = interp1(facecam_t(~isnan(facecam_t)),facecam.proc.data.whisker.motion(~isnan(facecam_t)),frame_t);
%use_trace = frame_spikes_conv;
%use_trace = wheel_speed;
%use_trace = roi_trace;
%use_trace = interp1(t,beta_power,frame_t);

use_trace_meansub = use_trace - mean(use_trace);

% Get cross correlation with all pixels
corr_lags = 35*30;

framerate = 1./nanmedian(diff(frame_t));
v_xcorr = nan(size(fV,1),corr_lags*2+1);
for curr_u = 1:size(U,3)  
    v_xcorr(curr_u,:) = xcorr(fV(curr_u,:)-mean(fV(curr_u,:)),use_trace_meansub,corr_lags,'biased');
end
lags_t = (-corr_lags:corr_lags)/framerate;
svd_xcorr = svdFrameReconstruct(U,v_xcorr);

% Normalize
svd_xcorr_norm = bsxfun(@rdivide,svd_xcorr,px_std*std(use_trace));

% Draw the movie
AP_image_scroll(svd_xcorr_norm,lags_t);


% % Correlation in pixel space
% U_downsample_factor = 10;
% maxlags = 35;
% 
% % Make mask
% h = figure;
% imagesc(avg_im);
% set(gca,'YDir','reverse');
% colormap(gray);
% caxis([0 prctile(avg_im(:),90)]);
% title('Draw mask to use pixels')
% roiMask = roipoly;
% delete(h);
% roiMaskd = imresize(roiMask,1/U_downsample_factor,'bilinear') > 0;
% 
% Ud = imresize(U,1/U_downsample_factor,'bilinear');
% Ud_flat = reshape(Ud,[],size(U,3));
% svd_corr = nan(size(Ud,1),size(Ud,2));
% for curr_px = find(roiMaskd)';
%     curr_trace = Ud_flat(curr_px,:)*fV;
%     %curr_corr = corrcoef(curr_trace,use_trace);
%     %svd_corr(curr_px) = curr_corr(2);
%     curr_corr = xcorr(curr_trace,use_trace,maxlags,'coeff');
%     svd_corr(curr_px) = max(curr_corr);
%     disp(curr_px/find(roiMaskd,1,'last'));
% end
% figure;imagesc(svd_corr);colormap(gray)



%% Cluster pixel correlations, get fluorescence in clusters

U_downsample_factor = 10;

% Make mask
h = figure;
imagesc(avg_im);
colormap(gray);
caxis([0 prctile(avg_im(:),90)]);
title('Draw mask to use pixels');
roiMask = roipoly;
delete(h);
roiMaskd = imresize(roiMask,1/U_downsample_factor,'bilinear') > 0;

Ud = imresize(U,1/U_downsample_factor,'bilinear');
Ud = bsxfun(@times,Ud,roiMaskd);
px_traces = reshape(Ud,[],size(Ud,3))*fV;
brain_px = std(px_traces,[],2) > 0.0001;

% Get fluorescence traces by grouped pixels
use_kgrps = 6;

kidx = kmeans(px_traces(brain_px,:),use_kgrps,'Distance','correlation');

kgrp_downsample = zeros(size(Ud,1),size(Ud,2));
kgrp_downsample(brain_px) = kidx;

% Sample back up to original size, fix edges with brain mask
kgrp = imresize(kgrp_downsample,[size(U,1),size(U,2)],'nearest');
kgrp(~any(U,3)) = 0;

% Get borders (for other stuff)
kgrp = imresize(kgrp_downsample,[size(U,1),size(U,2)],'nearest');
kgrp_borders = arrayfun(@(x) bwboundaries(kgrp == x),1:use_kgrps,'uni',false);

figure;
kgrp_im = imagesc(kgrp);
colormap(hot);
drawnow;
% If desired, cut areas bilaterally
bilateral_separate = true;
if bilateral_separate
    % Select midline
    [y,x] = ginput(1);
    y = round(y);
    % Everything on one side of the midline, add max group number
    kgrp(:,y:end) = kgrp(:,y:end) + ...
        use_kgrps*(kgrp(:,y:end) ~= 0);
    % Relabel with consecutive numbers (in case something skipped)
    unique_groups = unique(kgrp(kgrp ~= 0));
    for curr_grp = 1:length(unique_groups)
       kgrp(kgrp == unique_groups(curr_grp)) = curr_grp; 
    end
    set(kgrp_im,'CData',kgrp);    
    drawnow;
end

% Get fluorescence across session in ROI
% (there must be a faster way to do this...)
disp('Getting cluster traces...')
cluster_trace = zeros(max(kgrp(:)),size(fV,2));
U_reshape = reshape(U,[],size(U,3));
for i = 1:max(kgrp(:))
    roiMask = kgrp == i;
    cluster_trace(i,:) = nanmean(U_reshape(roiMask(:),:)*fV);
end
disp('Done')


%% Topography by correlation?
% not exactly sure where I was going with this... reversals in map, i.e.
% can pick up both visual areas and somatomotor areas based on correlation?

U_downsample_factor = 10;
Ud = imresize(U,1/U_downsample_factor,'bilinear');

% Lifted from pixelCorrelationViewerSVD

% to compute just the correlation with one pixel and the rest:
% 1) ahead of time:
fprintf(1, 'pre-computation...\n');
Ur = reshape(U, size(U,1)*size(U,2),[]); % P x S
covV = cov(fV'); % S x S % this is the only one that takes some time really
varP = dot((Ur*covV)', Ur'); % 1 x P
fprintf(1, 'done.\n');

ySize = size(U,1); xSize = size(U,2);

pixel = [349,360];
pixelInd = sub2ind([ySize, xSize], pixel(1), pixel(2));

covP = Ur(pixelInd,:)*covV*Ur'; % 1 x P
stdPxPy = varP(pixelInd).^0.5 * varP.^0.5; % 1 x P
corrMat = reshape(covP./stdPxPy,ySize,xSize); % 1 x P


corr_max_thresh = 0.005;
curr_local_max = imextendedmax(corrMat,corr_max_thresh);

% iterate below

U_downsample_factor = 10;
Ud = imresize(U,1/U_downsample_factor,'bilinear');

ySize = size(Ud,1);
xSize = size(Ud,2);

fprintf(1, 'pre-computation...\n');
Ur = reshape(Ud, size(Ud,1)*size(Ud,2),[]); % P x S
covV = cov(fV'); % S x S % this is the only one that takes some time really
varP = dot((Ur*covV)', Ur'); % 1 x P
fprintf(1, 'done.\n');

use_y = 1:ySize;
use_x = 1:xSize;

corr_max_thresh = 0.001;

max_pix = false(ySize,xSize,length(use_x));
max_corr_map = zeros(ySize,xSize);
max_corr_x_idx = zeros(ySize,xSize);
for y_idx = 1:length(use_y)
    y = use_y(y_idx);
    for x_idx = 1:length(use_x)
        x = use_x(x_idx);
        pixel = [y,x];
        pixelInd = sub2ind([ySize, xSize], pixel(1), pixel(2));
        
        covP = Ur(pixelInd,:)*covV*Ur'; % 1 x P
        stdPxPy = varP(pixelInd).^0.5 * varP.^0.5; % 1 x P
        corrMat = reshape(covP./stdPxPy,ySize,xSize); % 1 x P
        
        curr_local_max = imextendedmax(corrMat,corr_max_thresh);
        
        % get rid of the local max that's the selected pixel
        curr_selected = bwselect(curr_local_max,pixel(2),pixel(1));
        curr_local_max_external = curr_local_max;
        curr_local_max_external(curr_selected) = false;
        
        curr_local_max_corr = zeros(ySize,xSize);
        curr_local_max_corr(curr_local_max_external) = corrMat(curr_local_max_external);
        
        replace_idx = curr_local_max_corr > max_corr_map;
        max_corr_x_idx(replace_idx) = x;
        
        max_pix(:,:,x_idx) = curr_local_max;       
    end
    disp(y_idx/length(use_y));
end

max_pix = false(ySize,xSize,length(use_y));
corr_pix = false(ySize,xSize,length(use_y));
max_corr_map = zeros(ySize,xSize);
max_corr_y_idx = zeros(ySize,xSize);
for x_idx = 1:length(use_x)
    x = use_x(x_idx);
    for y_idx = 1:length(use_y)
        y = use_y(y_idx);
        pixel = [y,x];
        pixelInd = sub2ind([ySize, xSize], pixel(1), pixel(2));
        
        covP = Ur(pixelInd,:)*covV*Ur'; % 1 x P
        stdPxPy = varP(pixelInd).^0.5 * varP.^0.5; % 1 x P
        corrMat = reshape(covP./stdPxPy,ySize,xSize); % 1 x P
        
        curr_local_max = imextendedmax(corrMat,corr_max_thresh);
        
        % get rid of the local max that's the selected pixel
        curr_selected = bwselect(curr_local_max,pixel(2),pixel(1));
        curr_local_max_external = curr_local_max;
        curr_local_max_external(curr_selected) = false;
        
        curr_local_max_corr = zeros(ySize,xSize);
        curr_local_max_corr(curr_local_max_external) = corrMat(curr_local_max_external);
        
        replace_idx = curr_local_max_corr > max_corr_map;
        max_corr_y_idx(replace_idx) = y;
        
        corr_pix(:,:,y_idx) = corrMat;
        max_pix(:,:,y_idx) = curr_local_max;
    end
    disp(x_idx/length(use_x));
end

% get gradient direction
[Xmag,Xdir] = imgradient(imgaussfilt(max_corr_x_idx,1));
[Ymag,Ydir] = imgradient(imgaussfilt(max_corr_y_idx,1));

angle_diff = sind(Xdir-Ydir);

figure;imagesc(angle_diff);









max_pix = false(ySize,xSize,length(use_y));
corr_pix = false(ySize,xSize,length(use_y));
max_corr_map = zeros(ySize,xSize);
max_corr_y_idx = zeros(ySize,xSize);
for x_idx = 1:length(use_x)
    x = use_x(x_idx);
    for y_idx = 1:length(use_y)
        y = use_y(y_idx);
        pixel = [y,x];
        pixelInd = sub2ind([ySize, xSize], pixel(1), pixel(2));
        
        covP = Ur(pixelInd,:)*covV*Ur'; % 1 x P
        stdPxPy = varP(pixelInd).^0.5 * varP.^0.5; % 1 x P
        corrMat = reshape(covP./stdPxPy,ySize,xSize); % 1 x P
        
        curr_local_max = imextendedmax(corrMat,corr_max_thresh);
        
        % get rid of the local max that's the selected pixel
        curr_selected = bwselect(curr_local_max,pixel(2),pixel(1));
        curr_local_max_external = curr_local_max;
        curr_local_max_external(curr_selected) = false;
        
        curr_local_max_corr = zeros(ySize,xSize);
        curr_local_max_corr(curr_local_max_external) = corrMat(curr_local_max_external);
        
        replace_idx = curr_local_max_corr > max_corr_map;
        max_corr_y_idx(replace_idx) = y;
        
        corr_pix(:,:,y_idx) = corrMat;
        max_pix(:,:,y_idx) = curr_local_max;
    end
    disp(x_idx/length(use_x));
end






%% xcorr trace/facecam, this probably doesn't make any sense... 

surround_frames = 60;

facecam_t_idx = facecam_t > min(frame_t) & facecam_t < max(frame_t);

% interpolate ROI trace to facecam trace
roi_trace_interp = zscore(interp1(frame_t,roi_trace,facecam_t(facecam_t_idx)));

use_frame_idx = find(facecam_t_idx);

facecam_vr = VideoReader(facecam_filename);

facecam_first_frame = read(facecam_vr,1);
facecam_surround = double(repmat(facecam_first_frame,[1,1,surround_frames*2+1]));

for curr_frame_idx = surround_frames+1:length(use_frame_idx)-surround_frames
    
    curr_facecam_frame = use_frame_idx(curr_frame_idx);
    
    curr_trace = roi_trace_interp(curr_frame_idx-surround_frames:curr_frame_idx+surround_frames);
    
    facecam_im_pre = read(facecam_vr,curr_facecam_frame-1);
    facecam_im = read(facecam_vr,curr_facecam_frame);
    
    facecam_im_diff = abs(facecam_im - facecam_im_pre);
    
    curr_facecam_surround = double(repmat(facecam_im_diff,[1,1,surround_frames*2+1]));
    curr_facecam_surround = bsxfun(@times,curr_facecam_surround,permute(curr_trace,[2,3,1]));
    
    facecam_surround = facecam_surround+curr_facecam_surround;
    
    disp(curr_frame_idx./length(use_frame_idx));
    
end


%% Get STD of pixels (for use in getting normalized xcorr)

px_std_sq = zeros(size(U,1),size(U,2));

px_mean = svdFrameReconstruct(U,nanmean(fV,2));

% Do this in chunks of 1000 frames
% Don't use the first n frames (can be weird)
skip_start_frames = 35*10;
n_frames = size(fV,2) - skip_start_frames + 1;
chunk_size = 1000;
frame_chunks = unique([skip_start_frames:chunk_size:size(fV,2),size(fV,2)]);

for curr_chunk = 1:length(frame_chunks)-1
    curr_im = svdFrameReconstruct(U,fV(:,frame_chunks(curr_chunk): ...
        frame_chunks(curr_chunk+1)));
    px_std_sq = px_std_sq + sum((bsxfun(@minus,curr_im,px_mean).^2./n_frames),3);
    disp(curr_chunk/(length(frame_chunks)-1));
end
px_std = sqrt(px_std_sq);

% don't know if this is legit
px_10prct = svdFrameReconstruct(U,prctile(fV(:,skip_start_frames:end),10,2));


%% Get average fluorescence to ChoiceWorld event

% Define the window to get an aligned response to
surround_window = [-0.2,1];

% Define the times to align to
use_trials = ismember(choiceworld.trialContrastValues,[1]) &  ...
    ismember(choiceworld.trialSideValues,[-1]) & ...
    ismember(choiceworld.hitValues,[1]);
align_times = choiceworld.stimOnTimes(use_trials(1:length(choiceworld.stimOnTimes)))';

% Get the surround time
framerate = 1./nanmedian(diff(frame_t));
surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);

% Don't use times that fall outside of imaging
align_times(align_times + surround_time(1) < frame_t(2) | ...
    align_times + surround_time(2) > frame_t(end)) = [];

% Use closest frames to times
align_surround_times = bsxfun(@plus, align_times, surround_time);
frame_edges = [frame_t,frame_t(end)+1/framerate];
align_frames = discretize(align_surround_times,frame_edges);

% If any aligned V's are NaNs (when does this happen?), don't use
align_frames(any(isnan(align_frames),2),:) = [];

aligned_V = reshape(fV(:,align_frames'), ...
    size(fV,1),size(align_frames,2),size(align_frames,1));

% (to bootstrap) this doesn't help
% n_boot = 1000;
% mean_aligned_V = ....
%     squeeze(nanmean(reshape(bootstrp(n_boot,@mean, ...
%     permute(aligned_V,[3,1,2])),n_boot,size(aligned_V,1),[]),1));

mean_aligned_V = nanmean(aligned_V,3);

% Get and plot the average fluorescence around event
mean_aligned_px = svdFrameReconstruct(U,mean_aligned_V);

AP_image_scroll(mean_aligned_px,surround_time);
warning off; truesize; warning on;


%% Regress fluorescence to choiceworld event

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = (frame_t > skip_seconds);

% Make choiceworld event trace
use_trials = ismember(choiceworld.trialContrastValues,[1,0.5]) &  ...
    ismember(choiceworld.trialSideValues,[-1]) & ...
    ismember(choiceworld.hitValues,[1]);
align_times = choiceworld.stimOnTimes(use_trials(1:length(choiceworld.stimOnTimes)))';

frame_edges = [frame_t,frame_t(end)+1/framerate];
choiceworld_event_trace = histcounts(align_times,frame_edges);
 
% for stim position, this didn't really work though...
% choiceworld_event_trace = interp1(choiceworld.stimAzimuthTimes,choiceworld.stimAzimuthValues,frame_t);
% choiceworld_event_trace(isnan(choiceworld_event_trace)) = 0;

% for timeline inputs
%choiceworld_event_trace = licking_trace;
%choiceworld_event_trace = wheel_speed;

use_svs = 1:50;
kernel_frames = -35:7;
downsample_factor = 1;
lambda = 1e6;
zs = false;
cvfold = 5;

kernel_frames_downsample = round(downsample(kernel_frames,downsample_factor)/downsample_factor);

[k,predicted_choiceworld_events,explained_var] = ...
    AP_regresskernel(downsample(fV(use_svs,use_frames)',downsample_factor)', ...
    downsample(choiceworld_event_trace(:,use_frames)',downsample_factor)',kernel_frames_downsample,lambda,zs,cvfold);

% Reshape kernel and convert to pixel space
r = reshape(k,length(use_svs),length(kernel_frames_downsample),size(choiceworld_event_trace,1));

r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3);
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(U(:,:,use_svs),r(:,:,curr_spikes));
end

AP_image_scroll(r_px,(kernel_frames_downsample*downsample_factor)/framerate);
caxis([prctile(r_px(:),[1,99])]*4);
truesize


%% Regress fluorescence to 2 opposing choiceworld events

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = (frame_t > skip_seconds);

% Make choiceworld event trace

use_trials_1 = ismember(choiceworld.trialContrastValues,[1,0.5]) &  ...
    ismember(choiceworld.trialSideValues,[1]) & ...
    ismember(choiceworld.hitValues,[1]);
align_times_1 = choiceworld.interactiveOnTimes(use_trials_1(1:length(choiceworld.interactiveOnTimes)))';

use_trials_2 = ismember(choiceworld.trialContrastValues,[1,0.5]) &  ...
    ismember(choiceworld.trialSideValues,[1]) & ...
    ismember(choiceworld.hitValues,[0]);
align_times_2 = choiceworld.interactiveOnTimes(use_trials_2(1:length(choiceworld.interactiveOnTimes)))';

frame_edges = [frame_t,frame_t(end)+1/framerate];
choiceworld_event_trace_1 = histcounts(align_times_1,frame_edges);
choiceworld_event_trace_2 = histcounts(align_times_2,frame_edges);

choiceworld_event_trace_combined = choiceworld_event_trace_1 - choiceworld_event_trace_2;
 
use_svs = 1:50;
kernel_frames = -35:7;
downsample_factor = 1;
lambda = 1e7;
zs = false;
cvfold = 5;

kernel_frames_downsample = round(downsample(kernel_frames,downsample_factor)/downsample_factor);

[k,predicted_choiceworld_events,explained_var] = ...
    AP_regresskernel(downsample(fV(use_svs,use_frames)',downsample_factor)', ...
    downsample(choiceworld_event_trace_combined(:,use_frames)',downsample_factor)',kernel_frames_downsample,lambda,zs,cvfold);

% Reshape kernel and convert to pixel space
r = reshape(k,length(use_svs),length(kernel_frames_downsample),size(choiceworld_event_trace_combined,1));

r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3);
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(U(:,:,use_svs),r(:,:,curr_spikes));
end

AP_image_scroll(r_px,(kernel_frames_downsample*downsample_factor)/framerate);
caxis([prctile(r_px(:),[1,99])]*4);
truesize

%% Regress fluorescence to multiple choiceworld events

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = (frame_t > skip_seconds);

% Make choiceworld event trace
frame_edges = [frame_t,frame_t(end)+1/framerate];
choiceworld_event_trace = [];

use_trials = ismember(choiceworld.trialContrastValues,[1,0.5]) &  ...
    ismember(choiceworld.trialSideValues,[-1]) & ...
    ismember(choiceworld.hitValues,[1]);
align_times = choiceworld.stimOnTimes(use_trials(1:length(choiceworld.stimOnTimes)))';
choiceworld_event_trace = [choiceworld_event_trace;histcounts(align_times,frame_edges)];

use_trials = ismember(choiceworld.trialContrastValues,[1,0.5]) &  ...
    ismember(choiceworld.trialSideValues,[1]) & ...
    ismember(choiceworld.hitValues,[1]);
align_times = choiceworld.stimOnTimes(use_trials(1:length(choiceworld.stimOnTimes)))';
choiceworld_event_trace = [choiceworld_event_trace;histcounts(align_times,frame_edges)];

% for stim position, this didn't really work though...
% choiceworld_event_trace = interp1(choiceworld.stimAzimuthTimes,choiceworld.stimAzimuthValues,frame_t);
% choiceworld_event_trace(isnan(choiceworld_event_trace)) = 0;

% for timeline inputs
%choiceworld_event_trace = licking_trace;
%choiceworld_event_trace = wheel_speed;

use_svs = 1:50;
kernel_frames = -35:35;
downsample_factor = 1;
lambda = 1e7;
zs = false;
cvfold = 5;

kernel_frames_downsample = round(downsample(kernel_frames,downsample_factor)/downsample_factor);

[k,predicted_choiceworld_events,explained_var] = ...
    AP_regresskernel(downsample(fV(use_svs,use_frames)',downsample_factor)', ...
    downsample(choiceworld_event_trace(:,use_frames)',downsample_factor)',kernel_frames_downsample,lambda,zs,cvfold);

% Reshape kernel and convert to pixel space
r = reshape(k,length(use_svs),length(kernel_frames_downsample),size(choiceworld_event_trace,1));

r_px = zeros(size(U,1),size(U,2),size(r,2),size(r,3),'single');
for curr_spikes = 1:size(r,3);
    r_px(:,:,:,curr_spikes) = svdFrameReconstruct(U(:,:,use_svs),r(:,:,curr_spikes));
end

AP_image_scroll(r_px,(kernel_frames_downsample*downsample_factor)/framerate);
caxis([prctile(r_px(:),[1,99])]*4);
truesize

%% Regress choiceworld event(s) to fluorescence

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = (frame_t > skip_seconds);

% Make choiceworld event traces
frame_edges = [frame_t,frame_t(end)+1/framerate];
choiceworld_event_trace = [];

use_contrasts = [0,0.125,0.25,0.5,1];
use_trialSides = [-1,1];

for trialSide_idx = 1:length(use_trialSides)
    for contrast_idx = 1:length(use_contrasts)
        
        curr_trialSide = use_trialSides(trialSide_idx);
        curr_contrast = use_contrasts(contrast_idx);
        
        use_trials = ismember(choiceworld.trialContrastValues,[curr_contrast]) &  ...
            choiceworld.trialSideValues == curr_trialSide & ...
            choiceworld.hitValues == 1;
        align_times = choiceworld.stimOnTimes(use_trials(1:length(choiceworld.stimOnTimes)))';
        choiceworld_event_trace = [choiceworld_event_trace;histcounts(align_times,frame_edges)];
        
    end
end

%choiceworld_event_trace = [choiceworld_event_trace;licking_trace];
%choiceworld_event_trace = [choiceworld_event_trace;wheel_speed];

%choiceworld_event_trace = zscore(choiceworld_event_trace,[],2);

use_svs = 1:50;
kernel_frames = -7:35;
lambda = 0;
zs = false;
cvfold = 5;

[k,predicted_fluor,explained_var] = ...
    AP_regresskernel(choiceworld_event_trace(:,use_frames), ...
    fV(use_svs,use_frames), ...
    kernel_frames,lambda,zs,cvfold);

% Reshape kernel and convert to pixel space
k_r = permute(reshape(k,size(choiceworld_event_trace,1),length(kernel_frames),length(use_svs)),[3,2,1]);

r_px = zeros(size(U,1),size(U,2),size(k_r,2),size(k_r,3),'single');
for curr_event = 1:size(k_r,3);
    r_px(:,:,:,curr_event) = svdFrameReconstruct(U(:,:,use_svs),k_r(:,:,curr_event));
end

AP_image_scroll(r_px,kernel_frames/framerate);
caxis([prctile(r_px(:),[1,99])]*4);
truesize

%% Align widefield images across days (/ get transform matricies)

animal = 'AP017';
expInfo_path = ['\\zserver.cortexlab.net\Data\expInfo\' animal];
expInfo_dir = dir(expInfo_path);
days = {expInfo_dir(find([expInfo_dir(3:end).isdir])+2).name};

avg_im = cell(length(days),1);
for curr_day = 1:length(days)
    data_path = ['\\zserver.cortexlab.net\Data\Subjects\' animal filesep days{curr_day}];
    avg_im{curr_day} = readNPY([data_path filesep 'meanImage_blue.npy']);
end

% border_pixels = 20;
% avg_im = cellfun(@(x) imgaussfilt(x(border_pixels:end-border_pixels+1,border_pixels:end-border_pixels+1),3),avg_im,'uni',false);

% Register days to each other
ref_im_num = length(avg_im);

disp('Registering average images')
tform_matrix = cell(length(avg_im),1);
tform_matrix{1} = eye(3);

avg_im_reg = nan(size(avg_im{ref_im_num},1),size(avg_im{ref_im_num},2),length(avg_im));
avg_im_reg(:,:,ref_im_num) = avg_im{ref_im_num};

for curr_session = setdiff(1:length(avg_im),ref_im_num)
    
    % OLD
    
    %     [optimizer, metric] = imregconfig('monomodal');
    %     optimizer.MaximumStepLength = 0.02;
    %     %optimizer.MaximumStepLength = 0.0001;
    %     optimizer.RelaxationFactor = 0.1;
    %     optimizer.GradientMagnitudeTolerance = 1e-5;
    %     optimizer.MaximumIterations = 300;
    %     %         optimizer = registration.optimizer.OnePlusOneEvolutionary();
    %     %         optimizer.MaximumIterations = 500;
    %     %         optimizer.GrowthFactor = 1.00001;
    %     %         optimizer.InitialRadius = 1e-5;
    %
    %     % Perform the registration on the maximum image
    %     tformEstimate = imregtform(avg_im{curr_session},avg_im{1},'affine',optimizer,metric);
    %     %tformEstimate = imregcorr(summed_max(:,:,curr_session),summed_max(:,:,1),'similarity');
    %
    %     curr_im_reg = imwarp(avg_im{curr_session},tformEstimate,'Outputview',imref2d(size(avg_im{1})));
    %
    %     tform_matrix{curr_session} = tformEstimate.T;
    %     summed_max_reg(:,:,curr_session) = max_im_reg;
    %     summed_mean_reg(:,:,curr_session) = mean_im_reg;
     
    
    %%%%%%%%%%%%%
    
    % This is to do correlation, then affine (if above doesn't work)
    [optimizer, metric] = imregconfig('monomodal');
    optimizer = registration.optimizer.OnePlusOneEvolutionary();
    optimizer.MaximumIterations = 1000;
    optimizer.GrowthFactor = 1+1e-6;
    optimizer.InitialRadius = 1e-4;
    
%     % Register 1) correlation
%     tformEstimate_corr = imregcorr(avg_im{curr_session},avg_im{ref_im_num},'similarity');
%     curr_im_reg_corr = imwarp(avg_im{curr_session},tformEstimate_corr,'Outputview',imref2d(size(avg_im{ref_im_num})));
%     
%     % Register 2) affine
%     tformEstimate_affine = imregtform(curr_im_reg_corr,avg_im{ref_im_num},'affine',optimizer,metric);
%     
%     tformEstimate_combined = tformEstimate_corr;
%     tformEstimate_combined.T = tformEstimate_affine.T*tformEstimate_corr.T;
%     
%     curr_im_reg = imwarp(avg_im{curr_session},tformEstimate_combined,'Outputview',imref2d(size(avg_im{ref_im_num})));
    
%     tform_matrix{curr_session} = tformEstimate_combined.T;

    %%% TEST: just affine
    tformEstimate_affine = imregtform(avg_im{curr_session},avg_im{ref_im_num},'affine',optimizer,metric);
    curr_im_reg = imwarp(avg_im{curr_session},tformEstimate_affine,'Outputview',imref2d(size(avg_im{ref_im_num})));
    tform_matrix{curr_session} = tformEstimate_affine.T;
    %%%%
    
    avg_im_reg(:,:,curr_session) = curr_im_reg;
    
    disp(curr_session);
    
end

AP_image_scroll(avg_im_reg)














