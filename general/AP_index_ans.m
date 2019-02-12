function matrix_idx = AP_index_ans(matrix,idx)
% matrix_idx = AP_index_ans(matrix,idx)
%
% Return 'idx' values of 'matrix'
% (used for one-line operation and indexes, e.g.
% first_val = AP_index_ans(rand(1,10)*rand(10,2),1)

matrix_idx = matrix(idx);