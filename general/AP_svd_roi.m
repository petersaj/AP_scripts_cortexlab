function [roi_trace,roi_mask] = AP_svd_roi(U,V,guide_im,overlay,roi_mask)
% [roi_trace,roi_mask] = AP_svd_roi(U,V,guide_im,overlay,roi_mask)
%
% Draw roi on average SVD'd image, return reconstructed trace in time
% guide_im - can be an image or 'master' if master-aligned CCF/retinotopy
% overlay - optional, if drawing ROI mask, overlays this image
% roi_mask - optional, used if provided and prompted to draw otherwise (can
% also be Y x X x N rois)

% Choose ROI
if ~exist('roi_mask','var') || isempty(roi_mask)
    h = figure;
    
    ax_guide_im = axes;
    if isnumeric(guide_im)
        imagesc(ax_guide_im,guide_im);
        axis image off;
        colormap(ax_guide_im, gray);
        caxis(ax_guide_im, [prctile(guide_im(:),5) prctile(guide_im(:),95)]);
    elseif ischar(guide_im) && strcmp(guide_im,'master')
        imagesc(ax_guide_im,ones(size(U,1),size(U,2)));
        axis image off;
        caxis([0,1])
        colormap(ax_guide_im, gray);
        set(gca,'YDir','Reverse');
        AP_reference_outline('retinotopy','r');
        AP_reference_outline('ccf_aligned','k');
    end
    if exist('overlay','var') && ~isempty(overlay) && all(size(overlay) == size(guide_im))
        ax_overlay = axes;
        overlay_im = imagesc(ax_overlay,overlay);
        axis image off;
        colormap(ax_overlay, colormap_BlueWhiteRed);
        overlay_thresh = double(prctile(abs(overlay(:)),95));
        caxis(ax_overlay, [-overlay_thresh,overlay_thresh]);
        set(overlay_im,'AlphaData',0.5*mat2gray(abs(overlay),[0,overlay_thresh]));
    end
    
    roi_mask = roipoly;
    close(h);
end

% Get fluorescence across session in ROIs
roi_trace = nan(size(roi_mask,3),size(V,2));
for curr_roi = 1:size(roi_mask,3)
    U_roi = reshape(U(repmat(roi_mask(:,:,curr_roi),1,1,size(U,3))),[],size(U,3));
    roi_trace(curr_roi,:) = nanmean(U_roi*V);
end




