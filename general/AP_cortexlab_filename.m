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
% hardware
% imaging
% ephys
% ephysraw
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
        % CHECK SERVER2 IF IT DOESN'T EXIST
        if ~exist(filename,'file')
            filepath = [server2 filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep experiment filesep 'Protocol.mat'];
        end
        
    case 'eyecam'
        filepath = [server1 filesep 'Data\EyeCamera'];
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep 'eye.mj2'];
        % CHECK SERVER2 IF IT DOESN'T EXIST
        if ~exist(filename,'file')
            filepath = [server2 filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep experiment filesep 'eye.mj2'];
        end
        
    case 'facecam'
        filepath = [server1 filesep 'Data\EyeCamera'];
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep 'face.mj2'];
        % CHECK SERVER2 IF IT DOESN'T EXIST
        if ~exist(filename,'file')
            filepath = [server2 filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep experiment filesep 'face.mj2'];
        end
        
    % Old (with etGUI)
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
        
    case 'hardware'
        filepath = [server1 filesep 'Data\expInfo'];
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep day '_' experiment '_' animal '_hardwareInfo.mat'];
        % CHECK SERVER2 IF IT DOESN'T EXIST
        if ~exist(filename,'file')
            filepath = [server2 filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep experiment ...
            filesep day '_' experiment '_' animal '_hardwareInfo.mat'];
        end
        
    case 'imaging'
        filepath = [server1 filesep 'Data\Subjects'];
        filename = [filepath filesep animal filesep day site_dir filesep];
        % CHECK SERVER2 IF IT DOESN'T EXIST
        if ~exist(filename,'dir')
            filepath = [server2 filesep 'Subjects'];
            filename = [filepath filesep animal filesep day site_dir filesep];
        end
        
    case 'ephys'
        filepath = [server1 filesep 'Data\Subjects'];
        filename = [filepath filesep animal filesep day filesep 'ephys' filesep 'kilosort' site_dir filesep];
        % CHECK SERVER2 IF IT DOESN'T EXIST
        if ~exist(filename,'dir')
            filepath = [server2 filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep 'ephys' filesep 'kilosort' site_dir filesep];
        end
        
    case 'ephysraw'
        filepath = [server1 filesep 'Data\Subjects'];
        filename = [filepath filesep animal filesep day filesep 'ephys' site_dir filesep 'experiment1_10*-0_0.dat'];
        % CHECK SERVER2 IF IT DOESN'T EXIST
        if isempty(dir(filename))
            filepath = [server2 filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep 'ephys' site_dir filesep 'experiment1_10*-0_0.dat'];
        end
        % CHECK NEW OPEN EPHYS ON SERVER2 IF IT DOESN'T EXIST
        if isempty(dir(filename))
            filepath = [server2 filesep 'Subjects'];
            filename = [filepath filesep animal filesep day filesep ...
                'ephys' site_dir filesep 'experiment1' filesep 'recording1' ...
                filesep 'continuous' filesep 'Neuropix-3a-100.0' filesep 'continuous.dat'];
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

file_exists = exist(filename) > 0;



