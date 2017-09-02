function roi_trace = AP_svd_roi(U,V,avg_im)
% roi_trace = AP_svd_roi(U,V,avg_im)
%
% Draw roi on average SVD'd image, return reconstructed trace in time

% Choose ROI
h = figure;
imagesc(avg_im);
set(gca,'YDir','reverse');
colormap(gray);
caxis([prctile(avg_im(:),0) prctile(avg_im(:),95)]);
roiMask = roipoly;
close(h);

% Get fluorescence across session in ROI
U_roi = reshape(U(repmat(roiMask,1,1,size(U,3))),[],size(U,3));
roi_trace = nanmean(U_roi*V);




