function protocols = AP_list_experiments(animal,day)
% experiments = AP_list_experiments(animal,day)
%
% List experiment numbers and protocols for given animal/day
% (indicate if anything is doubled)

% Find all experiments for that day (numbered folder)
experiments_dir = dir(AP_cortexlab_filename(animal,day,1,'expInfo'));
experiments_num_idx = cellfun(@(x) ~isempty(x), regexp({experiments_dir.name},'^\d*$'));
experiment_num = cellfun(@str2num,{experiments_dir(experiments_num_idx).name});

% Get protocol for each experiment
protocols = struct('experiment',num2cell(experiment_num), ...
    'protocol',cell(size(experiment_num)), ...
    'multiple',false);

for curr_exp = 1:length(experiment_num)
    % Check for signals or MPEP
    [block_filename, block_exists] = AP_cortexlab_filename(animal,day,experiment_num(curr_exp),'block');
    [protocol_filename,protocol_exists] = AP_cortexlab_filename(animal,day,experiment_num(curr_exp),'protocol');
    
    if block_exists
        % If signals
        load(block_filename)
        if isfield(block,'expDef')
            [~,protocols(curr_exp).protocol] = fileparts(block.expDef);
        end
    elseif protocol_exists
        % If MPEP
        load(protocol_filename)
        [~,protocols(curr_exp).protocol] = fileparts(Protocol.xfile);
    else
        protocols(curr_exp).protocol = 'manual';
    end
end

% Indicate multiple experiments (usually because of error, last is good)
[~,last_unique_protocols] = unique({protocols.protocol},'last');
muliple_protocols = setdiff(1:length(protocols),last_unique_protocols);
[protocols(muliple_protocols).multiple] = deal(true);







