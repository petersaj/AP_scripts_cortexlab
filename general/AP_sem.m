function sem = AP_sem(data,dim)
% sem = AP_sem(data,dim);
% Standard error (nanstd/sqrt(sum(~isnan(x)))

% If dim unspecified and data is 1D, vectorize
if ~exist('dim','var') && sum(size(data) ~= 1)
    data = reshape(data,[],1);
    dim = 1;    
end

sem = nanstd(data,[],dim)./sqrt(sum(~isnan(data),dim));
