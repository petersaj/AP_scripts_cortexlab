function frameRecon = AP_svdFrameReconstruct(U,V)
% frameRecon = AP_svdFrameReconstruct(U,V)
%
% U is Y x X x nSVD
% V is nSVD x nFrames x nConditions
% 
% Copied from NS svdFrameReconstruct: added 3rd V dimension

frameRecon = reshape(reshape(U,[],size(U,3))* ...
    reshape(V,size(V,1),[]),size(U,1),size(U,2),size(V,2),size(V,3));