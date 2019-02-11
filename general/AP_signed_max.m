function signed_max = AP_signed_max(x,dim)
% signed_max = AP_signed_max(x,dim)
% Maximum amplitude (signed maximum absolute value)

% Maximum value and dimension subscript
[max_abs,max_abs_sub] = nanmax(abs(x),[],dim);

% (this is a ridiculous solution but it's the best I could come up with)
[~,dim_reorder_reverse] = sort([dim,setxor(dim,1:ndims(x))]);
dim_idx = permute([1:size(x,dim)]',dim_reorder_reverse);

max_abs_idx = bsxfun(@eq,max_abs_sub,dim_idx);
signed_max = nansum(x.*max_abs_idx,dim);

% sanity check that max values are indexed correctly (ignore NaNs)
nan_vals = isnan(signed_max) | isnan(max_abs);
if ~all(max_abs(~nan_vals) == abs(signed_max(~nan_vals)))
    error('AP_signed_max indexing error')
end




