function [filename,file_exists] = AP_cortexlab_filename(animal,day,experiment,file,site)
% [filename,file_exists] = AP_cortexlab_filename(animal,day,experiment,file,site)
%
% file - can include:
% expInfo
% timeline
% block
% parameters
% protocol
% eyecam
% eyecam_processed
% facecam
% facecam_processed
% facecam_movement
% hardware
% imaging
% ephys
% ephys_ks1
% ephys_dir
% ephys_ap
% probe_histology

% Make inputs strings if they're numbers
if isnumeric(experiment)
    experiment = num2str(experiment);
end

if isnumeric(day)
    experiment = num2str(day);
end

if exist('site','var') && ~isempty(site)
    if isnumeric(site)
        site = num2str(site);
    end
    site_dir = [filesep 'site' site];
else
    site = [];
    site_dir = [];
end

%%%% DATE FORMAT: NOT USED ANY MORE (the standard is yyyy-mm-dd)
% % There is no lab convention for dates so it's all mixed up, make both
% % versions and use as necessary
% isdaydash = any(strfind(day,'-'));
% if isdaydash
%     day_dash = day;
%     day_8digit = datestr(datenum(day,'yyyy-mm-dd'),'yyyymmdd');
% elseif ~isdaydash
%     day_dash = datestr(datenum(day,'yyyymmdd'),'yyyy-mm-dd');
%     day_8digit = day;
% end
%
% % Switch the used day format if necessary
% if exist('dayformat','var') && ~isempty(dayformat)
%     switch dayformat
%         case 'dash'
%             day = day_dash;
%         case '8digit'
%             day = day_8digit;
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% server1 = '\\zclone.cortexlab.net';
server1 = '\\zserver.cortexlab.net';
server2 = '\\zubjects.cortexlab.net';

% Check that servers are accessible (login needed on restart)
if ~exist([server1 filesep 'Data'])
    error('Zserver not available');
end
if ~exist([server2 filesep 'Subjects'])
    error('Zubjects not available');
end

switch file
    
    case 'expInfo'
        filepath = [server1 filesep 'Data\expInfo'];
        filename = [filepath filesep animal filesep day filesep];
        % CHECK SERVER2 IF IT DOESN'T EXIST
        if ~exist(filename,'file')
            filepath = [server2 filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep];
        end
        
    case 'timeline'
        filepath = [server1 filesep 'Data\expInfo'];
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep day '_' experiment '_' animal '_Timeline.mat'];
        % CHECK ANIMAL FOLDER IF IT DOESN'T EXIST
        if ~exist(filename,'file')
            filepath = [server1 filesep 'Data' filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep experiment ...
                filesep day '_' experiment '_' animal '_Timeline.mat'];
        end
        % CHECK SERVER2 IF IT DOESN'T EXIST
        if ~exist(filename,'file')
            filepath = [server2 filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep experiment ...
                filesep day '_' experiment '_' animal '_Timeline.mat'];
        end
        
    case 'block'
        filepath = [server1 filesep 'Data\expInfo'];
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep day '_' experiment '_' animal '_Block.mat'];
        % CHECK SERVER2 IF IT DOESN'T EXIST
        if ~exist(filename,'file')
            filepath = [server2 filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep experiment ...
                filesep day '_' experiment '_' animal '_Block.mat'];
        end
        
    case 'parameters'
        filepath = [server1 filesep 'Data\expInfo'];
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep day '_' experiment '_' animal '_parameters.mat'];
        % CHECK SERVER2 IF IT DOESN'T EXIST
        if ~exist(filename,'file')
            filepath = [server2 filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep experiment ...
                filesep day '_' experiment '_' animal '_parameters.mat'];
        end
        
    case 'protocol'
        filepath = [server1 filesep 'Data\trodes'];
        filename = [filepath filesep animal filesep day filesep experiment filesep 'Protocol.mat'];
        % CHECK ANIMAL FOLDER IF IT DOESN'T EXIST
        if ~exist(filename,'file')
            filepath = [server1 filesep 'Data' filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep experiment filesep 'Protocol.mat'];
        end
        % CHECK SERVER2 IF IT DOESN'T EXIST
        if ~exist(filename,'file')
            filepath = [server2 filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep experiment filesep 'Protocol.mat'];
        end
        
    case 'eyecam'
        filepath = [server1 filesep 'Data\EyeCamera'];
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep 'eye.mj2'];
        % CHECK ANIMAL FOLDER IF IT DOESN'T EXIST
        if ~exist(filename,'file')
            filepath = [server1 filesep 'Data' filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep experiment filesep 'eye.mj2'];
        end
        % CHECK SERVER2 IF IT DOESN'T EXIST
        if ~exist(filename,'file')
            filepath = [server2 filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep experiment filesep 'eye.mj2'];
        end
        
    case 'facecam'
        filepath = [server1 filesep 'Data\EyeCamera'];
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep 'face.mj2'];
        % CHECK ANIMAL FOLDER IF IT DOESN'T EXIST
        if ~exist(filename,'file')
            filepath = [server1 filesep 'Data' filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep experiment filesep 'face.mj2'];
        end
        % CHECK SERVER2 IF IT DOESN'T EXIST
        if ~exist(filename,'file')
            filepath = [server2 filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep experiment filesep 'face.mj2'];
        end
        
        % UNUSED: processing with etGUI/eyeGUI/facemap
        
        %     % Old (with etGUI)
        %     case 'eyecam_processed'
        %         filepath = '\\zserver1.cortexlab.net\Data\EyeCamera';
        %         filename = [filepath filesep mouse filesep day filesep experiment ...
        %             filesep day '_' experiment '_' mouse '_eye_processed.mat'];
        
        % New (with eyeGUI)
    case 'eyecam_processed'
        filepath = [server1 filesep 'Data\EyeCamera'];
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep 'eye_proc.mat'];
        
    case 'facecam_processed'
        filepath = [server1 filesep 'Data\EyeCamera'];
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep 'face_proc.mat'];
        
        % Output from AP_mouse_movie_movement
    case 'facecam_movement'
        filepath = [server1 filesep 'Data\EyeCamera'];
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep 'facecam_movement.mat'];
        % CHECK SERVER2 IF IT DOESN'T EXIST
        if ~exist(filename,'file')
            filepath = [server2 filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep experiment filesep 'facecam_movement.mat'];
        end
        
    case 'hardware'
        filepath = [server1 filesep 'Data\expInfo'];
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep day '_' experiment '_' animal '_hardwareInfo.mat'];
        % CHECK ANIMAL FOLDER IF IT DOESN'T EXIST
        if ~exist(filename,'file')
            filepath = [server1 filesep 'Data' filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep experiment ...
            filesep day '_' experiment '_' animal '_hardwareInfo.mat'];
        end
        % CHECK SERVER2 IF IT DOESN'T EXIST
        if ~exist(filename,'file')
            filepath = [server2 filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep experiment ...
                filesep day '_' experiment '_' animal '_hardwareInfo.mat'];
        end
        
    case 'imaging'
        filepath = [server1 filesep 'Data\Subjects'];
        filename = [filepath filesep animal filesep day site_dir filesep];
        check_file = [filename filesep 'svd*'];
        
        % CHECK SERVER2 IF IT DOESN'T EXIST
        if isempty(dir(check_file))
            filepath = [server2 filesep 'Subjects'];
            filename = [filepath filesep animal filesep day site_dir filesep];
            check_file = [filename filesep 'svd*'];
        end
        
        file_exists = ~isempty(dir(check_file));
        
    case 'ephys_dir'
        % the path where the ephys data is kept
        filename = [server1 filesep 'Data\Subjects' filesep animal filesep day filesep 'ephys' site_dir];
        
        % CHECK SERVER2 IF IT DOESN'T EXIST
        if ~exist(filename,'dir')
            filename = [server2 filesep 'Subjects' filesep animal filesep day filesep 'ephys' site_dir];
        end
        
    case 'ephys_ap'
        % the raw action potential band data file
        
        filepath = [server1 filesep 'Data\Subjects' filesep animal filesep day filesep 'ephys' site_dir];
        filepattern = [filepath filesep 'experiment1_10*-0_0.dat'];
        filepattern_dir = dir(filepattern);
        filename = [filepath filesep filepattern_dir.name];
        
        % CHECK SERVER2 IF IT DOESN'T EXIST
        if isempty(filepattern_dir)
            filepath = [server2 filesep 'Subjects' filesep animal filesep day filesep 'ephys' site_dir];
            filepattern = [filepath filesep 'experiment1_10*-0_0.dat'];
            filepattern_dir = dir(filepattern);
            filename = [filepath filesep filepattern_dir.name];
        end
        
        % CHECK NEW OPEN EPHYS ON SERVER2 IF IT DOESN'T EXIST
        % (after server switch, should never be on zserver)
        if isempty(filepattern_dir)
            filepath = [server2 filesep 'Subjects' filesep animal filesep day filesep ...
                'ephys' site_dir filesep 'experiment1' filesep 'recording1'];
            filename = [filepath filesep 'continuous'  filesep 'Neuropix-3a-100.0' filesep 'continuous.dat'];
        end
        
    case 'ephys'
        % folder with kilosort/phy outputs
        
        kilosort_version = 2; % (kilosort 2 by default)
        
        % Self-call to match ephys_dir above
        [ephys_dir_site,~] = AP_cortexlab_filename(animal,day,experiment,'ephys_dir',site);
        ephys_dir_parent = ephys_dir_site(1:end-length(site_dir));
   
        % Drop the kilosort version in the base workspace
        assignin('base','kilosort_version',kilosort_version); 
        
        if kilosort_version == 1
            filename = [ephys_dir_parent filesep 'kilosort' filesep site_dir];
        elseif kilosort_version == 2
            filename = [ephys_dir_parent filesep 'kilosort2' filesep site_dir];
        end
        
    case 'ephys_ks1'
        % folder with kilosort/phy outputs
        
        kilosort_version = 1; % (kilosort 2 by default)
        
        % Self-call to match ephys_dir above
        [ephys_dir_site,~] = AP_cortexlab_filename(animal,day,experiment,'ephys_dir',site);
        ephys_dir_parent = ephys_dir_site(1:end-length(site_dir));
   
        % Drop the kilosort version in the base workspace
        assignin('base','kilosort_version',kilosort_version); 
        
        if kilosort_version == 1
            filename = [ephys_dir_parent filesep 'kilosort' filesep site_dir];
        elseif kilosort_version == 2
            filename = [ephys_dir_parent filesep 'kilosort2' filesep site_dir];
        end
        
    case 'probe_histology'
        % (the output from P Shamash's program for histology to probe)
        filepath = [server2 filesep 'Subjects'];
        filepattern = [filepath filesep animal filesep 'histology' filesep 'probe' filesep 'processed' filesep 'probe_points*.mat'];
        filepattern_dir = dir(filepattern);
        if ~isempty(filepattern_dir) && length(filepattern_dir) == 1
            filename = [filepattern_dir.folder filesep filepattern_dir.name];
        else
            filename = filepattern;
        end
        
end

% Check whether the file exists (unless that check was performed above)
if ~exist('file_exists','var')
    file_exists = any(exist(filename));
end







