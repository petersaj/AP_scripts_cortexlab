%% If the V's were parsed but the timestamps weren't saved

animal = 'AP032';
day = '2018-10-27';
experiments = [1,2,3];

[data_path,file_exists] = AP_cortexlab_filename(animal,day,[],'imaging');
dataSummary_fn = [data_path filesep 'dataSummary_blue'];
load(dataSummary_fn);

cam_times_blue = cell(length(experiments),1);
cam_times_purple = cell(length(experiments),1);
exp_start_frames = [1,find(diff(dataSummary.timeStampsFromStamp) > 2)+1,length(dataSummary.timeStampsFromStamp)+1];

for curr_exp_idx = 1:length(experiments)   
    curr_exp = experiments(curr_exp_idx);
    
    % Load timeline
    timeline_filename = AP_cortexlab_filename(animal,day,curr_exp,'timeline');
    load(timeline_filename);
    
    timeline_cam_idx = strcmp({Timeline.hw.inputs.name}, 'pcoExposure');
    cam_samples = find(Timeline.rawDAQData(1:end-1,timeline_cam_idx) <= 2 & ...
        Timeline.rawDAQData(2:end,timeline_cam_idx) > 2) + 1;

    cam_time = Timeline.rawDAQTimestamps(cam_samples);
    
    % Load V's to check
    V_blue_fn = [data_path filesep num2str(curr_exp) filesep 'svdTemporalComponents_blue.npy'];
    V_purple_fn = [data_path filesep num2str(curr_exp) filesep 'svdTemporalComponents_purple.npy'];
    V_blue = readNPY(V_blue_fn);
    V_purple = readNPY(V_purple_fn);
    
    if length(cam_samples) ~= size(V_blue,1) + size(V_purple,1)
        error('Wrong number of recorded frames');
    end
    
    cam_times_blue{curr_exp_idx} = cam_time(1:2:end);
    cam_times_purple{curr_exp_idx} = cam_time(2:2:end);
end

frame_times_fn = 'svdTemporalComponents_blue.timestamps.npy';
for curr_exp_idx = 1:length(experiments)   
    curr_exp = experiments(curr_exp_idx);
    curr_exp_dir = [data_path filesep num2str(curr_exp)];
    
    if ~exist(curr_exp_dir,'dir')
        mkdir(curr_exp_dir);
    end
    writeNPY(cam_times_blue{curr_exp_idx},[curr_exp_dir filesep frame_times_fn]);
    disp('Saved blue timestamps')
end


frame_times_fn = 'svdTemporalComponents_purple.timestamps.npy';
for curr_exp_idx = 1:length(experiments)   
    curr_exp = experiments(curr_exp_idx);
    curr_exp_dir = [data_path filesep num2str(curr_exp)];
    
    if ~exist(curr_exp_dir,'dir')
        mkdir(curr_exp_dir);
    end
    writeNPY(cam_times_purple{curr_exp_idx},[curr_exp_dir filesep frame_times_fn]);
    disp('Saved purple timestamps')
end


%% If the V's were parsed incorrectly
% (doesn't save anything, just loads in to look)

animal = 'AP027';
day = '2017-11-25';
experiments = [1,5,6];

[data_path,file_exists] = AP_cortexlab_filename(animal,day,[],'imaging');

dataSummary_blue_fn = [data_path filesep 'dataSummary_blue'];
dataSummary_blue = load(dataSummary_blue_fn);

dataSummary_purple_fn = [data_path filesep 'dataSummary_purple'];
dataSummary_purple = load(dataSummary_purple_fn);

cam_times_blue = cell(size(experiments));
cam_times_purple = cell(size(experiments));

V_blue = cell(size(experiments));
V_purple = cell(size(experiments));
for curr_exp_idx = 1:length(experiments)
    
    curr_exp = experiments(curr_exp_idx);
    
    % Load frame times
    timeline_filename = AP_cortexlab_filename(animal,day,curr_exp,'timeline');
    load(timeline_filename);
    
    timeline_cam_idx = strcmp({Timeline.hw.inputs.name}, 'pcoExposure');
    cam_samples = find(Timeline.rawDAQData(1:end-1,timeline_cam_idx) <= 2 & ...
        Timeline.rawDAQData(2:end,timeline_cam_idx) > 2) + 1;

    cam_time = Timeline.rawDAQTimestamps(cam_samples);    
    cam_times_blue{curr_exp_idx} = cam_time(1:2:end);
    cam_times_purple{curr_exp_idx} = cam_time(2:2:end);
    
    % Load V's
    V_blue_fn = [data_path filesep num2str(curr_exp) filesep 'svdTemporalComponents_blue.npy'];
    V_purple_fn = [data_path filesep num2str(curr_exp) filesep 'svdTemporalComponents_purple.npy'];
    V_blue{curr_exp_idx} = readNPY(V_blue_fn);
    V_purple{curr_exp_idx} = readNPY(V_purple_fn);

    AP_print_progress_fraction(curr_exp_idx,length(experiments));
end

% Checks
if any(diff(sort([dataSummary_blue.dataSummary.frameNumbersFromStamp, ...
        dataSummary_purple.dataSummary.frameNumbersFromStamp])) ~= 1)
    error('Frame numbers not continuous');
end

time_stamp_cat = sort([dataSummary_blue.dataSummary.timeStampsFromStamp, ...
    dataSummary_purple.dataSummary.timeStampsFromStamp]);
exp_borders = [0,find(diff(time_stamp_cat) > 2),length(time_stamp_cat)];
    
if ~all(diff(exp_borders) == ...
        cellfun(@(x) size(x,1),V_blue) + cellfun(@(x) size(x,1),V_purple))
    error('V''s parsed incorrectly');
end

frame_diff = ...
    (cellfun(@(x) length(x),cam_times_blue) + cellfun(@(x) length(x),cam_times_purple)) - ...
    (cellfun(@(x) size(x,1),V_blue) + cellfun(@(x) size(x,1),V_purple));












