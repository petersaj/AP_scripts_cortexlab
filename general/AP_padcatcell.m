function x_padcat = AP_padcatcell(x)
% x_padcat = AP_padcatcell(x)
%
% Pad (with NaNs) and concatenate cell array

% Concatenate along the first consistently singleton dimension
x_size = cell2mat(cellfun(@size,reshape(x,[],1),'uni',false));
x_size_max = max(x_size,[],1);
cat_dim = find(~all(x_size == 1,1),1,'last') + 1;

x_pad = cellfun(@(x) padarray(x,x_size_max-size(x),NaN,'post'),x,'uni',false);
x_padcat = cat(cat_dim,x_pad{:});

