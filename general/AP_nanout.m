function idx_nanout = AP_nanout(idx)
% idx_nanout = AP_nanout(idx)
% 
% Turn logical array (0/1) into (1/NaN), used to NaN-out values in an nd
% array with multiplication
%
% e.g. x = [2,5;-2,5], want to get the mean of non-negative values:
% nanmean(x.*AP_nanout(x < 0))
%
% Input: idx - logical matrix (1/0)
% Output: idx_nanout - 1/NaN


idx_nanout = +~idx;
idx_nanout(idx) = NaN;
