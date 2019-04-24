function spatial_explained_var = AP_spatial_explained_var(U,V,V_pred,downsample_factor)
% spatial_explained_var = AP_spatial_explained_var(U,V,V_pred,downsample_factor)
%
% Reconstructs traces and gets variance explained
% downsample_factor - amount to spatially downsample U

% Get measure and predicted fluorescence trace for each pixel
Ud = imresize(U,1/downsample_factor,'bilinear');
fluor_px = svdFrameReconstruct(Ud,V);
fluor_px_pred = svdFrameReconstruct(Ud,V_pred);

% Find total explained variance (R^2)
nan_frames = any(isnan(reshape(fluor_px,[],size(fluor_px,3))),1) | ...
    any(isnan(reshape(fluor_px_pred,[],size(fluor_px_pred,3))),1);
sse_residual = sum((fluor_px(:,:,~nan_frames)-fluor_px_pred(:,:,~nan_frames)).^2,3);
sse_total = sum((fluor_px(:,:,~nan_frames)-nanmean(fluor_px(:,:,~nan_frames),3)).^2,3);
spatial_explained_var_d = 1 - (sse_residual./sse_total);

% Upsample explained variance map to original size
spatial_explained_var = imresize(spatial_explained_var_d,size(U(:,:,1)),'bilinear');
