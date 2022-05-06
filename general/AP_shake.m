function x_shaken = AP_shake(x,dim,grp_idx)
% data_shaken = AP_shake(data,dim,grp_idx)
%
% Randomly permute data along the specified direction
% dim - dimension to shake (can be multiple dimensions)
% grp_idx - grouping to shake within (only shake elements within group)

% If no dim input, major axis if vector and break if not
if ~exist('dim','var') || isempty(dim)
    if sum(size(x) > 1) == 1
        dim = find(size(x) > 1);
    else
        error('No dimension specified on non-vector data');
    end
end

% Grouping index: 
if exist('grp_idx','var') && ~isempty(grp_idx)
    % Ensure column group matrix
    grp_idx = reshape(grp_idx,[],1);

    % Check that group index matches size of shaken dimensions
    if length(grp_idx) ~= size(x,dim)
        error('Grouping variable size doesn''t match shaken elements')
    end
else
    % If no group index, put everything in one group
    grp_idx = ones(prod(size(x,dim)),1);
end

x_size = size(x);

% Put the relevant dimension(s) first
dim_reorder = [dim,setxor(dim,1:ndims(x))];
x_shaken = permute(x,dim_reorder);

% Reshape the data to be shakable x non-shakable
x_shaken = reshape(x_shaken,prod(x_size(dim)),[]);

% Shake within each group separately
for curr_grp = reshape(unique(grp_idx),1,[])

    curr_grp_idx = find(grp_idx == curr_grp);
    
    % Run through each column and shake the data
    [~,shaken_column_idx] = sort(rand(size(x_shaken(curr_grp_idx,:))),1);
    shaken_element_idx = curr_grp_idx(shaken_column_idx) + size(x_shaken,1)*(0:size(x_shaken,2)-1);
    x_shaken(curr_grp_idx,:) = x_shaken(shaken_element_idx);

end

% Reshape and permute the shaken data back to normal
[~,dim_reorder_reverse] = sort(dim_reorder);
x_shaken = permute(reshape(x_shaken,x_size(dim_reorder)),dim_reorder_reverse);















