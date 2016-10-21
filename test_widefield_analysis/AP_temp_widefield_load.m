%% Convert StackSet to dat file

test_savefile = 'C:\Users\Andy\Documents\CarandiniHarrisLab\data\widefield_test\test.dat';

[px_y,px_x,num_frames] = size(singleStack.Values);

int16_max = 2^15-1;
im_reshape = reshape(int16(mat2gray(singleStack.Values)*int16_max),[],size(singleStack.Values,3));



%% Try SVDs on test dat file

ops.RegFile = test_savefile;
[ops.Ly,ops.Lx,ops.Nframes] = size(singleStack.Values);


[ops, U, Sv, V, totalVar] = get_svdcomps(ops);


%% Try SVDs on just the data from the matlab image

im_minsub = bsxfun(@minus, double(im_reshape), mean(im_reshape,2));
COV = im_minsub' * im_minsub/size(im_minsub,1);

% total variance of data. If you ask for all Svs back then you will see
% this is equal to sum(Sv). In this case Sv are the singular values *of the
% covariance matrix* not of the original data - they are equal to the Sv of
% the original data squared (the variances per dimension).
totalVar = sum(diag(COV));

%[V, Sv] = svd(gpuArray(double(COV)));
[V, Sv] = svd(double(COV));
V = gather(V);
Sv = gather(Sv);
U = normc(im_minsub * V);
U = single(U);
Sv = single(diag(Sv));


%% SVD movie?

U_r = reshape(U,px_y,px_x,[]);




















