function h = AP_reference_outline(type,color,reference_im)
% h = AP_reference_outline(type,color,reference_im)
%
% If a reference image is provided, bregma is defined there
% type - grid, ccf, ccf_wf, ccf_aligned, retinotopy (master only at the moment)

if ~exist('color','var')
    color = 'k';
end

%% Define bregma

% Hard-code microns per pixel
um2pixel = 20.6;

if exist('reference_im','var') && ~isempty(reference_im)
    temp_f = figure;
    imagesc(reference_im);
    colormap(gray);
    axis off image;
    [bregma_offset_x,bregma_offset_y] = ginput(1);
    close(temp_f);
end

hold on;

switch type   
    
    case 'grid'
        % Plot mm grid for widefield
        spacing_um = 1000;
        spacing_pixels = spacing_um/um2pixel;
        
        xlines_pos = bregma_offset_y + spacing_pixels*(ceil((min(ylim)-bregma_offset_y)./spacing_pixels):floor((max(ylim)-bregma_offset_y)./spacing_pixels));
        ylines_pos = bregma_offset_x + spacing_pixels*(ceil((min(xlim)-bregma_offset_x)./spacing_pixels):floor((max(xlim)-bregma_offset_x)./spacing_pixels));
        
        h = struct;
        
        for curr_xline = 1:length(xlines_pos)
            h.xlines(curr_xline) = line(xlim,repmat(xlines_pos(curr_xline),1,2),'color',color,'linestyle','--');
        end
        
        for curr_yline = 1:length(ylines_pos)
            h.ylines(curr_yline) = line(repmat(ylines_pos(curr_yline),1,2),ylim,'color',color,'linestyle','--');
        end
        
        h.bregma = plot(bregma_offset_x,bregma_offset_y,'.r','MarkerSize',30);
                
    case 'ccf'
        % Plot allen CCF outlines       
        load('C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\cortical_area_boundaries.mat');
        bregma = allenCCFbregma;
        for curr_area_idx =1:length(cortical_area_boundaries)
            h = cellfun(@(outline) plot((outline(:,2)-bregma(3))*10, ...
                (bregma(1)-outline(:,1))*10,color),cortical_area_boundaries{curr_area_idx},'uni',false);
        end
        
    case 'ccf_wf'        
        % Plot allen CCF outlines scaled for widefield
        load('C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\cortical_area_boundaries.mat');
        bregma = allenCCFbregma;
        for curr_area_idx =1:length(cortical_area_boundaries)
            h = cellfun(@(outline) plot(bregma_offset_x + ((outline(:,2)-bregma(3))*10)/um2pixel, ...
                bregma_offset_y + ((outline(:,1)-bregma(1))*10)/um2pixel,color),cortical_area_boundaries{curr_area_idx},'uni',false);
        end
        
    case 'ccf_aligned'
        % Plot CCF borders aligned to master retinotopy
        load('C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\cortical_area_boundaries_aligned.mat');
        h = cellfun(@(areas) cellfun(@(outline) plot(outline(:,2),outline(:,1),color),areas,'uni',false), ...
            cortical_area_boundaries_aligned,'uni',false);
        
    case 'retinotopy'
        % Plot master retinotopic borders        
        load('\\basket.cortexlab.net\data\ajpeters\retinotopy\retinotopic_boundaries.mat');
        h = cellfun(@(outline) plot(outline(:,2),outline(:,1),color),retinotopic_boundaries,'uni',false);
        
end


% ONLY NEEDED ONCE: create top-down cortical boundaries from CCF
% % Get first brain pixel from top-down, get annotation at that point
% [~,top_down_depth] = max(av>1, [], 2);
% top_down_depth = squeeze(top_down_depth);
%
% [xx,yy] = meshgrid(1:size(top_down_depth,2), 1:size(top_down_depth,1));
% top_down_annotation = reshape(av(sub2ind(size(av),yy(:),top_down_depth(:),xx(:))), size(av,1), size(av,3));
%
% % Get all labelled areas
% used_areas = unique(top_down_annotation(:));
%
% % Restrict to only cortical areas
% structure_id_path = cellfun(@(x) textscan(x(2:end),'%d', 'delimiter',{'/'}),st.structure_id_path);
%
% ctx_path = [997,8,567,688,695,315];
% ctx_idx = find(cellfun(@(id) length(id) > length(ctx_path) & ...
%     all(id(min(length(id),length(ctx_path))) == ctx_path(min(length(id),length(ctx_path)))),structure_id_path));
%
% plot_areas = intersect(used_areas,ctx_idx);
%
% bregma = allenCCFbregma;
%
% % Get outlines of all areas
% cortical_area_boundaries = cell(size(plot_areas));
% for curr_area_idx = 1:length(plot_areas)
%     cortical_area_boundaries{curr_area_idx} = bwboundaries(top_down_annotation == plot_areas(curr_area_idx));
% end













