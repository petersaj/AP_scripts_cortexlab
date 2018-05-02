function x_shaken = AP_shake(x,dim)
% data_shaken = AP_shake(data,dim)
%
% Randomly permute data along the specified direction
% dim can be multiple dimensions

x_size = size(x);

% Put the relevant dimension(s) first
dim_reorder = [dim,setxor(dim,1:ndims(x))];
x_shaken = permute(x,dim_reorder);

% Reshape the data to be shakable x non-shakable
x_shaken = reshape(x_shaken,prod(x_size(dim)),[]);

% Run through each column and shake the data
[~,shaken_column_idx] = sort(rand(size(x_shaken)),1);
shaken_element_idx = bsxfun(@plus,shaken_column_idx,size(x_shaken,1)*(0:size(x_shaken,2)-1));
x_shaken = x_shaken(shaken_element_idx);

% Reshape and permute the shaken data back to normal
[~,dim_reorder_reverse] = sort(dim_reorder);
x_shaken = permute(reshape(x_shaken,x_size(dim_reorder)),dim_reorder_reverse);


















