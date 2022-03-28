function grp_data = AP_groupfun(data,grp,fun)
% grp_data = AP_grpfun(data,grp,fun)
%
% Apply function to data grouped by grouping variables
%
% data - ND array
% grp - cell array of grouping variables (oriented in the demension they
% index, e.g. Nx1 = first dim, 1x1xN = third dim)
% fun - function to apply to grouped data
%
% I made this because 1) grpstats only works with 2D arrays and 1)
% accumarray requires ND index matricies which can be ridiculously huge,
% I'm surprised there's no better built-in function to do this...

error('ABANDONED FOR NOW - SUPER SLOW')

% Check that each grouping variable is only 1D
grp_isvec = cellfun(@(x) isvector(squeeze(x)),grp);
if ~all(grp_isvec)
    error('Not all groups are vectors');
end

% Check that the size of each grouping variable matches it's data
grp_match_datasize = cellfun(@(x) size(x,find(size(x)~=1)) == ...
    size(data,find(size(x)~=1)),grp);
if ~all(grp_match_datasize)
    error('Group vectors don''t match data size');
end

% Turn grouping variables into unique indicies
[grp_unique,~,grp_idx] = cellfun(@unique,grp,'uni',false);
grp_idx = cellfun(@(idx,grp) reshape(idx,size(grp)),grp_idx,grp,'uni',false);

grp_unique_n = cellfun(@length,grp_unique);
n_combinations = prod(cellfun(@(x) length(x),grp_unique));

% Initialize group data (as NaNs)
grp_data = nan(grp_unique_n,class(data));

data_size = size(data);

% Loop through all possible group combinations
% (set combination index as subscript)
curr_grp_sub = ones(size(grp));
curr_grp_subchange = 1;
for curr_comb = 1:n_combinations
    % Get current group combination
    curr_grp_idx = cellfun(@(x,idx) x == idx,grp_idx,num2cell(curr_grp_sub),'uni',false);
    % Loop through current combination indicies, get intersection
    curr_grp_idx_block = curr_grp_idx{1};
    for curr_grp_comb = 2:length(grp)
        curr_grp_idx_block = curr_grp_idx_block & curr_grp_idx{curr_grp_comb};
    end
    curr_grp_data = feval(fun,data(curr_grp_idx_block));

    % Convert current group subscripts to index
    curr_grp_idx = (curr_grp_sub-1)*cumprod([1 grp_unique_n(1:end-1)]')+1;

    % Store group data in output matrix
    grp_data(curr_grp_idx) = curr_grp_data;

    % Update current group subscripts
    curr_grp_sub(curr_grp_subchange) = curr_grp_sub(curr_grp_subchange)+1;

    %%%% UNFINISHED; NEED TO TICK OVER CURR_GRP_SUB
    AP_print_progress_fraction(curr_comb,n_combinations);

end






