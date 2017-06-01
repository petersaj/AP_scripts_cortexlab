function roi_trace = AP_svd_roi(U,V)
% roi_trace = AP_svd_roi(U,V)
%
% Draw roi on average SVD'd image, return reconstructed trace in time

% Choose ROI
h = figure;
avg_svd = reshape(reshape(U,[],size(U,3))*mean(V,2),size(U,1),size(U,2));
imagesc(avg_svd);
set(gca,'YDir','reverse');
colormap(gray);
caxis([prctile(avg_svd(:),0) prctile(avg_svd(:),100)]);
roiMask = roipoly;
close(h);

% Get fluorescence across session in ROI
U_roi = reshape(U(repmat(roiMask,1,1,size(U,3))),[],size(U,3));
roi_trace = nanmean(U_roi*V);




