function px = AP_svdFrameReconstruct(U,V)
% frameRecon = AP_svdFrameReconstruct(U,V)
%
% U is Y x X x nSVD
% V is nSVD x nFrames x ...[nConditions]

U_size = size(U);
V_size = size(V);

px = reshape(reshape(U,[],size(U,3))* ...
    reshape(V,size(V,1),[]),[U_size(1:2),V_size(2:end)]);