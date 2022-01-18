function experiments = AP_find_experiments(animal,protocol,flexible_name)
% experiments = AP_find_experiments(animal,protocol,flexible_name)
%
% flexible_name - if true, uses strfind(lower)
%
% Find all experiments from an animal with a given protocol name
% 
% Note: file locations changed many times so this is necessarily totally
% disorganized and ridiculous

% If no protocol specified, return all experiments
if ~exist('protocol','var') || isempty(protocol)
    protocol = [];
end

% If no flexible name specified, use exact
if ~exist('flexible_name','var') || isempty(flexible_name)
    flexible_name = false;
end


%% Server locations

% Initialize pathname, add to it with each server location
days_combined = {};
days_pathnames_combined = {};

% (server 1 expInfo - old)
expInfo_path = ['\\zserver.cortexlab.net\Data\expInfo\' animal];
expInfo_dir = dir(expInfo_path);
day_paths = cellfun(@(x) ~isempty(regexp(x,'\d\d\d\d-\d\d-\d\d')),{expInfo_dir.name}) &...
    [expInfo_dir.isdir];
curr_days = {expInfo_dir(day_paths).name};
curr_days_pathname = cellfun(@(x) [expInfo_path filesep x],curr_days,'uni',false);

days_combined = [days_combined,curr_days];
days_pathnames_combined = [days_pathnames_combined,curr_days_pathname];

% (server 1 subjects - new)
expInfo_path = ['\\zserver.cortexlab.net\Data\Subjects\' animal];
expInfo_dir = dir(expInfo_path);
day_paths = cellfun(@(x) ~isempty(regexp(x,'\d\d\d\d-\d\d-\d\d')),{expInfo_dir.name}) &...
    [expInfo_dir.isdir];
curr_days = {expInfo_dir(day_paths).name};
curr_days_pathname = cellfun(@(x) [expInfo_path filesep x],curr_days,'uni',false);

days_combined = [days_combined,curr_days];
days_pathnames_combined = [days_pathnames_combined,curr_days_pathname];

% (server 2)
expInfo_path = ['\\zubjects.cortexlab.net\Subjects\' animal];
expInfo_dir = dir(expInfo_path);
day_paths = cellfun(@(x) ~isempty(regexp(x,'\d\d\d\d-\d\d-\d\d')),{expInfo_dir.name}) &...
    [expInfo_dir.isdir];
curr_days = {expInfo_dir(day_paths).name};
curr_days_pathname = cellfun(@(x) [expInfo_path filesep x],curr_days,'uni',false);

days_combined = [days_combined,curr_days];
days_pathnames_combined = [days_pathnames_combined,curr_days_pathname];

% (server 3) 
expInfo_path = ['\\znas.cortexlab.net\Subjects\' animal];
expInfo_dir = dir(expInfo_path);
day_paths = cellfun(@(x) ~isempty(regexp(x,'\d\d\d\d-\d\d-\d\d')),{expInfo_dir.name}) &...
    [expInfo_dir.isdir];
curr_days = {expInfo_dir(day_paths).name};
curr_days_pathname = cellfun(@(x) [expInfo_path filesep x],curr_days,'uni',false);

days_combined = [days_combined,curr_days];
days_pathnames_combined = [days_pathnames_combined,curr_days_pathname];

% (server 4)
expInfo_path = ['\\zinu.cortexlab.net\Subjects\' animal];
expInfo_dir = dir(expInfo_path);
day_paths = cellfun(@(x) ~isempty(regexp(x,'\d\d\d\d-\d\d-\d\d')),{expInfo_dir.name}) &...
    [expInfo_dir.isdir];
curr_days = {expInfo_dir(day_paths).name};
curr_days_pathname = cellfun(@(x) [expInfo_path filesep x],curr_days,'uni',false);

days_combined = [days_combined,curr_days];
days_pathnames_combined = [days_pathnames_combined,curr_days_pathname];

% If multiple days: experiment number folders should be preserved across
% all folders so just pick one
[days,unique_day_idx] = unique(days_combined);
days_pathnames = days_pathnames_combined(unique_day_idx);


%% Find experiments
% Find experiments with chosen protocol and which modalities were recorded

protocol_expts = cell(size(days));
imaging_expts = cell(size(days));
ephys_expts = cell(size(days));

for curr_day = 1:length(days)  
    
    day = days{curr_day};
    % Find all experiments of that day (defined as number-only folders)
    expDay_dir = dir(days_pathnames{curr_day});
    exp_folders = cellfun(@any,regexp({expDay_dir.name},'^\d*$'));
    exp_nums = cellfun(@str2num,{expDay_dir(exp_folders).name});
    
    % If looking for specific protocol, find amongst days's experiments
    if ~isempty(protocol)
        use_exp = false(size(exp_nums));
        for curr_exp = 1:length(exp_nums)
            % Check for signals or MPEP
            [block_filename, block_exists] = AP_cortexlab_filename(animal,day,exp_nums(curr_exp),'block');
            [protocol_filename,protocol_exists] = AP_cortexlab_filename(animal,day,exp_nums(curr_exp),'protocol');
            
            if block_exists
                % If signals
                load(block_filename)
                if isfield(block,'expType') % old name
                    [~,expDef] = fileparts(block.expType);
                    if flexible_name
                        use_exp(curr_exp) = contains(lower(expDef),lower(protocol));
                    elseif ~flexible_name
                        use_exp(curr_exp) = strcmp(expDef,protocol);
                    end
                end
                if isfield(block,'expDef') % new name
                    [~,expDef] = fileparts(block.expDef);
                    if flexible_name
                        use_exp(curr_exp) = contains(lower(expDef),lower(protocol));
                    elseif ~flexible_name
                        use_exp(curr_exp) = strcmp(expDef,protocol);
                    end
                end
            elseif protocol_exists
                % If MPEP
                load(protocol_filename)
                [~,expDef] = fileparts(Protocol.xfile);
                if flexible_name
                    use_exp(curr_exp) = contains(lower(expDef),lower(protocol));
                elseif ~flexible_name
                    use_exp(curr_exp) = strcmp(expDef,protocol);
                end
            else
                continue
            end
        end
    else
        use_exp = true(size(exp_nums));
    end
    
    % Find days with imaging/electrophysiology
    if any(use_exp)
        protocol_expts{curr_day} = exp_nums(use_exp);
        imaging_path = AP_cortexlab_filename(animal,day,exp_nums(use_exp),'imaging');
        imaging_expts{curr_day} = exist([imaging_path filesep 'meanImage_blue.npy'],'file') > 0;
        [~,ephys_expts{curr_day}] = AP_cortexlab_filename(animal,day,exp_nums(use_exp),'ephys_dir');
    end
    
end

%% Package experiment info
use_days = ~cellfun(@isempty,protocol_expts);
experiments = struct('day',cell(sum(use_days),1),'experiment',cell(sum(use_days),1),...
    'imaging',cell(sum(use_days),1),'ephys',cell(sum(use_days),1));

[experiments.day] = deal(days{use_days});
[experiments.experiment] = deal(protocol_expts{use_days});
[experiments.imaging] = deal(imaging_expts{use_days});
[experiments.ephys] = deal(ephys_expts{use_days});









