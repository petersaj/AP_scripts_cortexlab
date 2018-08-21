function x_tril = AP_itril(x,k)
% x_tril = AP_itril(x,k)
% Return the vectorized lower triagular part of a matrix from kth diagonal


if nargin == 1
    k = 0;
end

tril_idx = tril(true(size(x)),k);
x_tril = x(tril_idx);
