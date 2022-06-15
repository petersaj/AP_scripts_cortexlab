% Batch preprocess local widefield data on Lilrig

% Grab this by just looking for populated folders
animal = 'AP113';
day = '2022-03-21';



%% Run SVD on day's data
% Single U for day, split V by recording

im_path = 'G:\test_widefield';
[U,Vrec,im_color_avg,frame_info] = AP_widefield_svd_pco(im_path);


%% Match frames to recording number
% Determine correct recording folder by closest folder and frame time

% Get recordings (numbered folders) in experiment path
experiment_path =  AP_cortexlab_filename(animal,day,1,'expInfo');
experiment_dir =  dir(experiment_path);
recording_dir_idx = cellfun(@(x) ~isempty(x), regexp({experiment_dir.name},'^\d*$'));
recording_dir = experiment_dir(recording_dir_idx);

% Find creation time of recording folders
rec_dir_starttime = nan(size(recording_dir));
for curr_recording = 1:length(recording_dir)
    curr_t = System.IO.File.GetCreationTime( ...
        fullfile(recording_dir(curr_recording).folder, ...
        recording_dir(curr_recording).name));
    
    rec_dir_starttime(curr_recording) = ...
        datenum(double([curr_t.Year,curr_t.Month,curr_t.Day, ...
        curr_t.Hour,curr_t.Minute,curr_t.Second]));
end

% Get nearest recording folder time for each imaging recording start
[~,rec_frame_start_idx] = unique([frame_info.rec_idx]);
frame_timestamp_cat = vertcat(frame_info.timestamp);
rec_im_time = frame_timestamp_cat(rec_frame_start_idx);

[~,im_rec_idx] = min(abs(rec_im_time - rec_dir_starttime'),[],2);


%% Save preprocessed widefield data on server

% Set number of components to save
n_components_save = 2000;

% Assume 2 colors in order of blue/purple
color_names = {'blue','purple'};

% Save frame information in experiment folder
frame_info_fn = fullfile(experiment_path,'widefield_frame_info');
save(frame_info_fn,'frame_info','-v7.3');

% Save mean images in experiment folder by color
for curr_color = 1:lenth(color_names)
    curr_mean_im_fn = fullfile(experiment_path, ...
        sprintf('meanImage_%s.npy',color_names{curr_color}));
    writeNPY(im_avg_color(:,:,curr_color),curr_mean_im_fn);
end

% Save spatial components in experiment (animal/day) folder by color
for curr_color = 1:lenth(color_names)


end





% Save temporal components in associated recording folders


fn = fullfile(Upath, ['svdSpatialComponents_' ops.vidName]);
fnMeanImage = fullfile(Upath, ['meanImage_' ops.vidName]);

if isfield(ops, 'saveAsNPY') && ops.saveAsNPY
    writeUVtoNPY(svdSpatialComponents, [], fn, []);
    writeNPY(meanImage, [fnMeanImage '.npy']);
else
    save(fn, '-v7.3', 'svdSpatialComponents');
    save(fnMeanImage, 'meanImage');
end




% Expects U to be yPix x xPix x nSV
% Expects V to be nSV x nTimePoints
% For U this is uncomplicated but for V you want to transpose first so that
% each component's time course will be written in order. This facilitates
% loading just a limited number of components later. 

if ~isempty(uFilePath)
    writeNPY(U, [uFilePath '.npy']);
end
if ~isempty(vFilePath)
    writeNPY(V', [vFilePath '.npy']);
end













