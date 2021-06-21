%% If the V's were parsed but the timestamps weren't saved
% This loads timeline, checks that the number of frames logged in timeline
% matches the frames recorded in the V's, and if so it saves the frame
% times in .timestamps file
%
% This is only used if something unusual happened in the preprocessing
% (most commonly because of a mismatch in actual and preprocessed
% experiments, e.g. experiments 1-4 were imaged but in the 'experiments'
% field of the preprocessing option there was only 2-4)
%
% If there is a larger issue e.g. timeline logged less frames than were
% actually recorded, manual inspection is necessary. In that case, run the
% next cell to load and check frame information.

animal = 'AP045'; % animal name
day = '2019-10-09'; % yyyy-mm-dd
experiments = [2,3,4]; % all experiments run that day (e.g. [1,2,3,4])

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

%% If the V's were not parsed into experiments

animal = 'AP104'; % animal name
day = '2021-06-05'; % yyyy-mm-dd
experiments = [1:2]; % all experiments run that day (e.g. [1,2,3,4])

[data_path,file_exists] = AP_cortexlab_filename(animal,day,[],'imaging');

dataSummary_blue_fn = [data_path filesep 'dataSummary_blue'];
dataSummary_blue = load(dataSummary_blue_fn);

dataSummary_purple_fn = [data_path filesep 'dataSummary_purple'];
dataSummary_purple = load(dataSummary_purple_fn);

cam_times_blue = cell(size(experiments));
cam_times_purple = cell(size(experiments));

[data_path,file_exists] = AP_cortexlab_filename(animal,day,[],'imaging');
dataSummary_fn = [data_path filesep 'dataSummary_blue'];
load(dataSummary_fn);

% Get cam times from each experiment
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
    
    cam_times_blue{curr_exp_idx} = cam_time(1:2:end);
    cam_times_purple{curr_exp_idx} = cam_time(2:2:end);
end

% Load unparsed V's
V_blue_fn = [data_path filesep 'svdTemporalComponents_blue.npy'];
V_purple_fn = [data_path filesep 'svdTemporalComponents_purple.npy'];
V_blue = readNPY(V_blue_fn);
V_purple = readNPY(V_purple_fn);

% Checks
if any(diff(sort([dataSummary_blue.dataSummary.frameNumbersFromStamp, ...
        dataSummary_purple.dataSummary.frameNumbersFromStamp])) ~= 1)
    error('Frame numbers not continuous');
end

time_stamp_cat = sort([dataSummary_blue.dataSummary.timeStampsFromStamp, ...
    dataSummary_purple.dataSummary.timeStampsFromStamp]);
exp_borders = [0,find(diff(time_stamp_cat) > 2),length(time_stamp_cat)];

if length(exp_borders) ~= length(experiments)+1
    error('Mismatching number of experiments and long frame gaps')
end

% Get number of timeline frames in each experiment
n_cam_blue = cellfun(@length,cam_times_blue);
n_cam_purple = cellfun(@length,cam_times_purple);

% To parse by number of timeline frames (if fine, just not parsed)
V_blue_exp = mat2cell(V_blue,n_cam_blue,size(V_blue,2));
V_purple_exp = mat2cell(V_purple,n_cam_purple,size(V_purple,2));

% % To manually parse (if something lost or messed up)
% V_blue_exp = mat2cell(V_blue,n_cam_blue + [0;0;size(V_blue,1)-sum(n_cam_blue);0],size(V_blue,1));
% V_purple_exp = mat2cell(V_purple,n_cam_purple + [0;0;size(V_purple,1)-sum(n_cam_purple);0],size(V_purple,1));

for curr_exp_idx = 1:length(experiments)
    curr_exp = experiments(curr_exp_idx);
    curr_exp_dir = [data_path filesep num2str(curr_exp)];
    
    % Save parsed V's
    V_blue_fn = [curr_exp_dir filesep 'svdTemporalComponents_blue.npy'];
    V_purple_fn = [curr_exp_dir filesep 'svdTemporalComponents_purple.npy'];
    
    writeNPY(V_blue_exp{curr_exp_idx},V_blue_fn);
    writeNPY(V_purple_exp{curr_exp_idx},V_purple_fn);
    
    % Save timestamps
    frame_times_fn = 'svdTemporalComponents_blue.timestamps.npy';
    writeNPY(cam_times_blue{curr_exp_idx},[curr_exp_dir filesep frame_times_fn]);
    
    frame_times_fn = 'svdTemporalComponents_purple.timestamps.npy';
    writeNPY(cam_times_purple{curr_exp_idx},[curr_exp_dir filesep frame_times_fn]);
    
    disp('Saved manually parsed V''s and timestamps');

end



%% If the V's were parsed incorrectly into experiments
% Loads in frame information to figure out what went wrong
% (no automatic fix)

animal = 'AP047'; % animal name
day = '2019-10-11'; % yyyy-mm-dd
experiments = [2:4]; % all experiments run that day (e.g. [1,2,3,4])

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

% Check that the total frames match the camera times
sum(cellfun(@length,cam_times_blue)) == sum(cellfun(@(x) size(x,1),V_blue));











