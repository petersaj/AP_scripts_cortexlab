function [filename,file_exists] = AP_cortexlab_filename(animal,day,experiment,file,site)
% filename = AP_cortexlab_filename(animal,day,experiment,file,site);
%
% file - can include:
% timeline
% block
% parameters
% protocol
% eyecam
% eyecam_processed
% facecam
% facecam_processed
% hardware
% datapath
%
% dayformat - dash or 8digit (default is to use whatever day is used as)

% I know this is redundant with dat.expFilePath, but I don't like that it
% uses a lot of server/nested functions to do something so simple. I also
% don't like reservation of 'dat'. So I just made this one.
%

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

switch file
    
    case 'timeline'
        filepath = '\\zserver.cortexlab.net\Data\expInfo';
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep day '_' experiment '_' animal '_Timeline.mat'];
        
    case 'block'
        filepath = '\\zserver.cortexlab.net\Data\expInfo';
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep day '_' experiment '_' animal '_Block.mat'];
        
    case 'parameters' 
        filepath = '\\zserver.cortexlab.net\Data\expInfo';
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep day '_' experiment '_' animal '_parameters.mat'];
        
    case 'protocol'
        filepath = '\\zserver.cortexlab.net\Data\trodes';
        filename = [filepath filesep animal filesep day filesep experiment filesep 'Protocol.mat'];
        
    case 'eyecam'
        filepath = '\\zserver.cortexlab.net\Data\EyeCamera';
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep 'eye.mj2'];
        
    case 'facecam'
        filepath = '\\zserver.cortexlab.net\Data\EyeCamera';
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep 'face.mj2'];
        
    % Old (with etGUI)
%     case 'eyecam_processed'
%         filepath = '\\zserver.cortexlab.net\Data\EyeCamera';
%         filename = [filepath filesep mouse filesep day filesep experiment ...
%             filesep day '_' experiment '_' mouse '_eye_processed.mat'];

    % New (with eyeGUI)
    case 'eyecam_processed'
        filepath = '\\zserver.cortexlab.net\Data\EyeCamera';
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep 'eye_proc.mat'];
        
    case 'facecam_processed'
        filepath = '\\zserver.cortexlab.net\Data\EyeCamera';
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep 'face_proc.mat'];
        
    case 'hardware'
        filepath = '\\zserver.cortexlab.net\Data\expInfo';
        filename = [filepath filesep animal filesep day filesep experiment ...
            filesep day '_' experiment '_' animal '_hardwareInfo.mat'];
        
    case 'imaging'
        filepath = '\\zserver.cortexlab.net\Data\Subjects';
        filename = [filepath filesep animal filesep day site_dir filesep];
        
    case 'ephys'
        filepath = '\\basket.cortexlab.net\data\ajpeters';
        filename = [filepath filesep animal filesep day filesep 'ephys' site_dir filesep];
        
    case 'ephysraw'
        filepath = '\\zserver.cortexlab.net\Data\Subjects';
        filename = [filepath filesep animal filesep day filesep 'ephys' site_dir filesep];

end

file_exists = exist(filename) > 0;



