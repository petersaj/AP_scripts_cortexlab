function sem = AP_sem(data,dim)
% sem = AP_sem(data,dim);
% Standard error (nanstd/sqrt(sum(~isnan(x)))

sem = nanstd(data,[],dim)./sqrt(sum(~isnan(data),dim));
