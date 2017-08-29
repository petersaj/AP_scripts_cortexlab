function spatial_explained_var = AP_spatial_explained_var(U,V,V_pred,downsample_factor)

Ud = imresize(U,1/downsample_factor,'bilinear');
f = svdFrameReconstruct(Ud,V);
f_pred = svdFrameReconstruct(Ud,V_pred);
spatial_explained_var_d = (nansum(f.^2,3)-nansum((f-f_pred).^2,3))./nansum(f.^2,3);
spatial_explained_var = imresize(spatial_explained_var_d,size(U(:,:,1)),'bilinear');
