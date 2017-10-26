function [roi_trace,roi_mask] = AP_svd_roi(U,V,avg_im,overlay,roi_mask)
% [roi_trace,roi_mask] = AP_svd_roi(U,V,avg_im,overlay,roi_mask)
%
% Draw roi on average SVD'd image, return reconstructed trace in time
% (or use already drawn mask)

% Choose ROI
if ~exist('roi_mask','var') || isempty(roi_mask)
    h = figure;
    imagesc(avg_im);
    set(gca,'YDir','reverse');
    axis image off;
    colormap(gray);
    caxis([prctile(avg_im(:),0) prctile(avg_im(:),99)]);
    roi_mask = roipoly;
    close(h);
end

% Get fluorescence across session in ROI
U_roi = reshape(U(repmat(roi_mask,1,1,size(U,3))),[],size(U,3));
roi_trace = nanmean(U_roi*V);




