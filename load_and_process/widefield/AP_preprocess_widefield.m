% Preprocess widefield data

% Grab this by just looking for populated folders
animal = 'AP113';
day = '2022-03-21';



%% Run SVD on day's data
% Single U for day, split V by recording

im_path = 'G:\test_widefield';
[U,Vrec,im_color_avg,frame_info] = AP_widefield_svd_pco(im_path);


%%

% Determine the recording number by the last created folder

% Get recording paths (numbered folders) in experiment path
experiment_dir =  dir(AP_cortexlab_filename(animal,day,1,'expInfo'));
recording_dir_idx = cellfun(@(x) ~isempty(x), regexp({experiments_dir.name},'^\d*$'));
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






















