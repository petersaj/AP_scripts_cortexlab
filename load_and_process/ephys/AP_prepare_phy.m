function AP_prepare_phy(data_filename)
% AP_prepare_phy(data_filename)
%
% Clears whatever is currently in local phy folder
% Sets up params file
% Loads raw and kilosorted data to local file
%
% data_filename = raw data filename

%% Get parts of raw filename

[data_path,data_file,data_ext] = fileparts(data_filename);


%% Clear out local phy folder

local_phy_dir = 'C:\data_temp\phy';
if ~exist(local_phy_dir)
    mkdir(local_phy_dir)
elseif length(dir(local_phy_dir)) > 2
    rmdir(local_phy_dir,'s');
    mkdir(local_phy_dir);
end

%% Write parameters

% Read in the new format header file, if it exists. If it doesn''t, write
% the defaults of the old format
header_path = [data_path filesep 'dat_params.txt'];
if exist(header_path,'file');
    
    header_fid = fopen(header_path);
    header_info = textscan(header_fid,'%s %s', 'delimiter',{' = '});
    fclose(header_fid);
    
    header = struct;
    for i = 1:length(header_info{1})
        header.(header_info{1}{i}) = header_info{2}{i};
    end
    
    params = {'dat_path',['''' data_file,data_ext '''']; ...
        'n_channels_dat',header.n_channels; ...
        'dtype',['''int16''']; ...
        'offset','0'; ...
        'sample_rate',header.sample_rate; ...
        'hp_filtered','False'};
    
else
    params = {'dat_path',['''' data_file,data_ext '''']; ...
        'n_channels_dat','129'; ...
        'dtype',['''int16''']; ...
        'offset','0'; ...
        'sample_rate','25000'; ...
        'hp_filtered','False'};
end

param_filename = [local_phy_dir filesep 'params.py'];

formatSpec = '%s = %s \r\n';
fid = fopen(param_filename,'w');
for curr_param = 1:size(params,1)
    fprintf(fid,formatSpec,params{curr_param,:});
end
fclose(fid);


%% Copy raw data and kilosorted data to local

disp('Copying data to local...')
copyfile(data_path,local_phy_dir);
disp('Done.');





































