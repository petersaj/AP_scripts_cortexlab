function experiments = AP_find_experiments(animal,protocol)
% experiments = AP_find_experiments(animal,protocol)
%
% Find all experiments from an animal with a given protocol name

expInfo_path = ['\\zserver.cortexlab.net\Data\expInfo\' animal];
expInfo_dir = dir(expInfo_path);
days = {expInfo_dir(find([expInfo_dir(3:end).isdir])+2).name};

% Find experiments with chosen protocol and which modalities were recorded
protocol_expts = cell(size(days));
imaging_expts = cell(size(days));
ephys_expts = cell(size(days));

for curr_day = 1:length(days)  
    day = days{curr_day};
    % Check all experiments of that day
    expDay_dir = dir([expInfo_path filesep days{curr_day}]);
    exp_nums = cellfun(@str2num,{expDay_dir(3:end).name});
    use_exp = false(size(exp_nums));
    for curr_exp = 1:length(exp_nums);
        % Check for signals or MPEP
        [block_filename, block_exists] = AP_cortexlab_filename(animal,day,exp_nums(curr_exp),'block');
        [protocol_filename,protocol_exists] = AP_cortexlab_filename(animal,day,exp_nums(curr_exp),'protocol');
        
        if block_exists
            % If signals
            load(block_filename)
            [~,expDef] = fileparts(block.expDef);
            use_exp(curr_exp) = ~isempty(strfind(expDef,protocol));
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
        [~,imaging_expts{curr_day}] = AP_cortexlab_filename(animal,day,exp_nums(use_exp),'imaging');
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









