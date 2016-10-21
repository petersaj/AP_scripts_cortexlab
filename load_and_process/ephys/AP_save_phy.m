function AP_save_phy
% AP_save_phy
%
% Move clustering data from current phy folder to basket
% NOTE: annoying because mpep means different animal names, can use this
% later I guess, have to manually do for old data...

% Set local phy data directory
local_phy_dir = 'C:\Users\Andrew\Documents\CarandiniHarrisLab\data\phy_temp';

% Get header info
header_path = [local_phy_dir filesep 'dat_params.txt'];

header_fid = fopen(header_path);
header_info = textscan(header_fid,'%s %s', 'delimiter',{' = '});
fclose(header_fid);

raw_path_idx = strcmp(header_info{1},'raw_path');
raw_path = header_info{2}{raw_path_idx};

% Get animal-specific folder
zserver_path = '\\zserver.cortexlab.net\Data\Subjects\';
basket_path = ['\\basket.cortexlab.net\data\ajpeters' raw_path(length(zserver_path)+1:end)];
