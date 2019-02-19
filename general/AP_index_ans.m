function matrix_idx = AP_index_ans(matrix,idx)
% matrix_idx = AP_index_ans(matrix,idx)
%
% Return 'idx' values of 'matrix'
% (used for one-line operation and indexes, e.g.
% first_val = AP_index_ans(rand(1,10)*rand(10,2),1)
%
% idx can be a number (e.g. [1,0,0], 3) or string (e.g. ':,3:end')

if ischar(idx)
    matrix_idx = eval(['matrix(' idx ')']);
else
    matrix_idx = matrix(idx);
end

