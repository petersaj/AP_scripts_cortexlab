function h = AP_reference_outline(reference_im)
% h = AP_reference_outline(reference_im)
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



%% for later? 

if false
    
    % (from Nick)
    % Get first brain pixel from top-down, get annotation at that point
    [~,top_down_depth] = max(av>1, [], 2);
    top_down_depth = squeeze(top_down_depth);
        
    [xx,yy] = meshgrid(1:size(top_down_depth,2), 1:size(top_down_depth,1));    
    top_down_annotation = reshape(av(sub2ind(size(av),yy(:),top_down_depth(:),xx(:))), size(av,1), size(av,3));
    
    % Get all labelled areas
    used_areas = unique(top_down_annotation(:));

    % Restrict to only cortical areas
    structure_id_path = cellfun(@(x) textscan(x(2:end),'%d', 'delimiter',{'/'}),st.structure_id_path);
    
    ctx_path = [997,8,567,688,695,315];
    ctx_idx = find(cellfun(@(id) length(id) > length(ctx_path) & ...
        all(id(min(length(id),length(ctx_path))) == ctx_path(min(length(id),length(ctx_path)))),structure_id_path));
    
    plot_areas = intersect(used_areas,ctx_idx);
    
    bregma = allenCCFbregma;
    
    figure; hold on;
    for curr_area = plot_areas'       
        curr_boundaries = bwboundaries(top_down_annotation == curr_area);
        cellfun(@(x) plot((x(:,2)-bregma(3))/100,(x(:,1)-bregma(1))/100,'k'),curr_boundaries,'uni',false);
    end
  
    
end

















