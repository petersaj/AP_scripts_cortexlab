
animal = 'AP014';
day = '2017-02-27';
experiments = 1:2;

data_path = ['\\zserver.cortexlab.net\Data\Subjects\' animal filesep day];
V_fn = [data_path filesep 'svdTemporalComponents_cam2.npy'];

V = readNPY(V_fn);

dataSummary_fn = [data_path filesep 'dataSummary_cam2'];
load(dataSummary_fn);

exp_start_frames = [1,find(diff(dataSummary.timeStampsFromStamp) > 2)+1,size(V,1)+1];

cam_times = cell(length(experiments),1);
for curr_exp_idx = 1:length(experiments)   
    curr_exp = experiments(curr_exp_idx);
    
    % Load timeline
    timeline_filename = get_cortexlab_filename(animal,day,curr_exp,'timeline');
    load(timeline_filename);
    
    timeline_cam_idx = strcmp({Timeline.hw.inputs.name}, 'cam2');
    cam_samples = find(Timeline.rawDAQData(1:end-1,timeline_cam_idx) <= 2 & ...
        Timeline.rawDAQData(2:end,timeline_cam_idx) > 2) + 1;
    cam_time = Timeline.rawDAQTimestamps(cam_samples);
    
    cam_times{curr_exp_idx} = cam_time;
end

% Sanity check: total cam_times = size V
if ~sum(cellfun(@length,cam_times)) == size(V,1);
    error('Wrong number of recorded frames')
end

temporal_components_fn = 'svdTemporalComponents_cam2.npy';
frame_times_fn = 'svdTemporalComponents_cam2.timestamps.npy';
for curr_exp_idx = 1:length(experiments)   
    curr_exp = experiments(curr_exp_idx);
    curr_exp_dir = [data_path filesep num2str(curr_exp)];
    
    mkdir(curr_exp_dir);
    writeNPY(cam_times{curr_exp_idx},[curr_exp_dir filesep frame_times_fn]);
    writeNPY(V(exp_start_frames(curr_exp_idx):exp_start_frames(curr_exp_idx+1)-1,:),[curr_exp_dir filesep temporal_components_fn]);
end



%% If the V's were parsed but the timestamps weren't saved

animal = 'AP014';
day = '2017-02-27';
experiments = 1:2;

data_path = ['\\zserver.cortexlab.net\Data\Subjects\' animal filesep day];
dataSummary_fn = [data_path filesep 'dataSummary_blue'];
load(dataSummary_fn);

cam_times = cell(length(experiments),1);
exp_start_frames = [1,find(diff(dataSummary.timeStampsFromStamp) > 2)+1,length(dataSummary.timeStampsFromStamp)-1];

for curr_exp_idx = 1:length(experiments)   
    curr_exp = experiments(curr_exp_idx);
    
    % Load timeline
    timeline_filename = get_cortexlab_filename(animal,day,curr_exp,'timeline');
    load(timeline_filename);
    
    timeline_cam_idx = strcmp({Timeline.hw.inputs.name}, 'pcoExposure');
    cam_samples = find(Timeline.rawDAQData(1:end-1,timeline_cam_idx) <= 2 & ...
        Timeline.rawDAQData(2:end,timeline_cam_idx) > 2) + 1;
    cam_time = Timeline.rawDAQTimestamps(cam_samples( ...
        exp_start_frames(curr_exp):exp_start_frames(curr_exp+1)-1));
    
    % Load V's to check
    V_blue_fn = [data_path filesep num2str(curr_exp) filesep 'svdTemporalComponents_blue.npy'];
    V_purple_fn = [data_path filesep num2str(curr_exp) filesep 'svdTemporalComponents_purple.npy'];
    V_blue = readNPY(V_blue_fn);
    V_purple = readNPY(V_purple_fn);
    
    if length(cam_time) ~= size(V,1)
        error('Wrong number of recorded frames');
    end
    
    cam_times{curr_exp_idx} = cam_time;
end

frame_times_fn = 'svdTemporalComponents_blue.timestamps.npy';
for curr_exp_idx = 1:length(experiments)   
    curr_exp = experiments(curr_exp_idx);
    curr_exp_dir = [data_path filesep num2str(curr_exp)];
    
    mkdir(curr_exp_dir);
    writeNPY(cam_times{curr_exp_idx},[curr_exp_dir filesep frame_times_fn]);
end




