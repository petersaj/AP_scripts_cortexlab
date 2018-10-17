function signed_max = AP_signed_max(x,dim)
% Maximum amplitude (signed maximum absolute value)

% Maximum value and dimension subscript
[max_abs,max_abs_sub] = max(abs(x),[],dim);

% (this is a ridiculous solution but it's the best I could come up with)
[~,dim_reorder_reverse] = sort([dim,setxor(dim,1:ndims(x))]);
dim_idx = permute([1:size(x,dim)]',dim_reorder_reverse);

max_abs_idx = bsxfun(@eq,max_abs_sub,dim_idx);
signed_max = sum(x.*max_abs_idx,dim);

% sanity check that max values are indexed correctly
if ~all(max_abs(:) == abs(signed_max(:)))
    error('indexing error')
end




