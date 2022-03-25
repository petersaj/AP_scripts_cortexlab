function group_data = AP_groupfun(data,idx,fun)
% group_data = AP_groupfun(data,idx,fun)
%
% Apply function to data grouped by index
%
% data - ND array
% idx - cell array of vector indicies (oriented in the demension they
% index, e.g. Nx1 = first dim, 1x1xN = third dim)
% fun - function to apply to grouped data
%
% I made this because 1) grpstats only works with 2D arrays and 1)
% accumarray requires ND index matricies which can be ridiculously huge,
% I'm surprised there's no better built-in function to do this...

keyboard






