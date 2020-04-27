function matrix_idx = AP_index_ans(matrix,varargin)
% matrix_idx = AP_index_ans(matrix,idx)
% This function is janky and bad practice - just meant for quick checks
%
% Return 'idx' values of 'matrix'
% (used for one-line operation and indexes, e.g.
% first_val = AP_index_ans(rand(1,10)*rand(10,2),1)
%
% idx can be one argument for each dimension or skip for full (e.g.
% use_trials,[],1)
% Or can be string (e.g. '(:,3:end)')

if length(varargin) == 1 && ischar(varargin{1})
    matrix_idx = eval(['matrix' idx '']);
elseif length(varargin) == ndims(matrix)
    idx = varargin;
    % empty dimensions: use full
    empty_idx = cellfun(@isempty,idx);
    idx(empty_idx) = cellfun(@(dim) 1:size(matrix,dim),num2cell(find(empty_idx)),'uni',false);
    % logical dimensions: turn into index
    logical_idx = cellfun(@islogical,idx);
    idx(logical_idx) = cellfun(@find,idx(logical_idx),'uni',false);
    % turn into string
    idx_char = cell2mat(cellfun(@(x) ['[' num2str(reshape(x,1,[])) '],'],idx,'uni',false));
    matrix_idx = eval(['matrix(' idx_char(1:end-1) ')']);
else
    error('Invalid index');
end

