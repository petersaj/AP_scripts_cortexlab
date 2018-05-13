function experiments = AP_find_experiments(animal,protocol)
% experiments = AP_find_experiments(animal,protocol)
%
% Find all experiments from an animal with a given protocol name

% Days can either be on Data\expInfo or Data2 (to be replaced, eventually)

% (look in server 1)
expInfo_path = ['\\zserver.cortexlab.net\Data\expInfo\' animal];
expInfo_dir = dir(expInfo_path);
day_paths = cellfun(@(x) ~isempty(regexp(x,'\d\d\d\d-\d\d-\d\d')),{expInfo_dir.name}) &...
    [expInfo_dir.isdir];
days_Data = {expInfo_dir(day_paths).name};
days_Data_pathname = cellfun(@(x) [expInfo_path filesep x],days_Data,'uni',false);

% (look in server 2)
expInfo_path = ['\\zubjects.cortexlab.net\Subjects\' animal];
expInfo_dir = dir(expInfo_path);
day_paths = cellfun(@(x) ~isempty(regexp(x,'\d\d\d\d-\d\d-\d\d')),{expInfo_dir.name}) &...
    [expInfo_dir.isdir];
days_Data2 = {expInfo_dir(day_paths).name};
days_Data2_pathname = cellfun(@(x) [expInfo_path filesep x],days_Data2,'uni',false);

% (take unique days from both)
days_combined = [days_Data,days_Data2];
days_pathname_combined = [days_Data_pathname,days_Data2_pathname];
[days,unique_day_idx] = unique(days_combined);
days_pathnames = days_pathname_combined(unique_day_idx);

% Find experiments with chosen protocol and which modalities were recorded
protocol_expts = cell(size(days));
imaging_expts = cell(size(days));
ephys_expts = cell(size(days));

for curr_day = 1:length(days)  
    day = days{curr_day};
    % Check all experiments of that day
    expDay_dir = dir(days_pathnames{curr_day});
    exp_nums = cellfun(@str2num,{expDay_dir(3:end).name});
    use_exp = false(size(exp_nums));
    for curr_exp = 1:length(exp_nums);
        % Check for signals or MPEP
        [block_filename, block_exists] = AP_cortexlab_filename(animal,day,exp_nums(curr_exp),'block');
        [protocol_filename,protocol_exists] = AP_cortexlab_filename(animal,day,exp_nums(curr_exp),'protocol');
        
        if block_exists
            % If signals
            load(block_filename)
            if isfield(block,'expDef');
                [~,expDef] = fileparts(block.expDef);
                use_exp(curr_exp) = ~isempty(strfind(expDef,protocol));
            end
        elseif protocol_exists
            % If MPEP
            load(protocol_filename)
            [~,expDef] = fileparts(Protocol.xfile);
            use_exp(curr_exp) = strcmp(expDef,protocol);
        else
            continue
        end       
    end
    
    if any(use_exp)
        protocol_expts{curr_day} = exp_nums(use_exp);
        imaging_path = AP_cortexlab_filename(animal,day,exp_nums(use_exp),'imaging');
        imaging_expts{curr_day} = exist([imaging_path filesep 'meanImage_blue.npy'],'file') > 0;
        [~,ephys_expts{curr_day}] = AP_cortexlab_filename(animal,day,exp_nums(use_exp),'ephys');
    end
end

% Package experiment info
use_days = ~cellfun(@isempty,protocol_expts);
experiments = struct('day',cell(sum(use_days),1),'experiment',cell(sum(use_days),1),...
    'imaging',cell(sum(use_days),1),'ephys',cell(sum(use_days),1));

[experiments.day] = deal(days{use_days});
[experiments.experiment] = deal(protocol_expts{use_days});
[experiments.imaging] = deal(imaging_expts{use_days});
[experiments.ephys] = deal(ephys_expts{use_days});









