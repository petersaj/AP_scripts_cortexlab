function h = AP_bregma_grid(reference_im)
% h = AP_bregma_grid(reference_im)
%
% If a reference image is provided, bregma is defined there

if exist('reference_im','var') && ~isempty(reference_im)
    temp_f = figure;
    imagesc(reference_im);
    colormap(gray);
    axis off image;
    [bregma_x,bregma_y] = ginput(1);
    close(temp_f);
else
    [bregma_x,bregma_y] = ginput(1);
end

hold on;

pixel2um = 20.6;
spacing_um = 1000;
spacing_pixels = spacing_um/pixel2um;

xlines_pos = bregma_y + spacing_pixels*(ceil((min(ylim)-bregma_y)./spacing_pixels):floor((max(ylim)-bregma_y)./spacing_pixels));
ylines_pos = bregma_x + spacing_pixels*(ceil((min(xlim)-bregma_x)./spacing_pixels):floor((max(xlim)-bregma_x)./spacing_pixels));

h = struct;

for curr_xline = 1:length(xlines_pos)
    h.xlines(curr_xline) = line(xlim,repmat(xlines_pos(curr_xline),1,2),'color','w','linestyle','--');
end

for curr_yline = 1:length(ylines_pos)
    h.ylines(curr_yline) = line(repmat(ylines_pos(curr_yline),1,2),ylim,'color','w','linestyle','--');
end

h.bregma = plot(bregma_x,bregma_y,'.r','MarkerSize',30);


