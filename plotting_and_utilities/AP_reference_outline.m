function h = AP_reference_outline(type,color,reference_im,repeats)
% h = AP_reference_outline(type,color,reference_im,repeats)
%
% If a reference image is provided, bregma is defined there
% type - grid_aligned, ccf, ccf_wf, ccf_aligned, retinotopy (master only at the moment)
%
% repeats - [pixels_y,pixels_x,repeats_y,repeats_x] (only for ccf_aligned)

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

% alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment';

switch type   
    
    case 'grid_aligned'
        % Plot mm grid for widefield
        
        % Get bregma from aligned (assume aligned atm)
        bregma = allenCCFbregma;
        bregma(3) = bregma(3) + 0.5;
        ccf_tform_fn = [alignment_path filesep 'ccf_tform.mat'];
        load(ccf_tform_fn);
        bregma_resize = bregma*(10/um2pixel);
        bregma_align = [bregma_resize([3,1]),1]*ccf_tform.T;
        bregma_offset_x = bregma_align(1);
        bregma_offset_y = bregma_align(2);
        
        spacing_um = 500;
        spacing_pixels = spacing_um/um2pixel;
        
        xlines_pos = bregma_offset_y + spacing_pixels*(ceil((min(ylim)-bregma_offset_y)./spacing_pixels):floor((max(ylim)-bregma_offset_y)./spacing_pixels));
        ylines_pos = bregma_offset_x + spacing_pixels*(ceil((min(xlim)-bregma_offset_x)./spacing_pixels):floor((max(xlim)-bregma_offset_x)./spacing_pixels));
        
        h = struct;
        
        for curr_xline = 1:length(xlines_pos)
            h.xlines(curr_xline) = line(xlim,repmat(xlines_pos(curr_xline),1,2),'color',color,'linestyle','-');
        end
        
        for curr_yline = 1:length(ylines_pos)
            h.ylines(curr_yline) = line(repmat(ylines_pos(curr_yline),1,2),ylim,'color',color,'linestyle','-');
        end
        
        h.bregma = plot(bregma_offset_x,bregma_offset_y,'xr','MarkerSize',10);
                
    case 'ccf'
        % Plot allen CCF outlines       
        load([alignment_path filesep 'cortical_area_boundaries.mat']);
        bregma = allenCCFbregma;
        for curr_area_idx =1:length(cortical_area_boundaries)
            h = cellfun(@(outline) plot((outline(:,2)-bregma(3))*10, ...
                (bregma(1)-outline(:,1))*10,'color',color),cortical_area_boundaries{curr_area_idx},'uni',false);
        end
        
    case 'ccf_wf'        
        % Plot allen CCF outlines scaled for widefield
        load([alignment_path filesep 'cortical_area_boundaries.mat']);
        bregma = allenCCFbregma;
        for curr_area_idx =1:length(cortical_area_boundaries)
            h = cellfun(@(outline) plot(bregma_offset_x + ((outline(:,2)-bregma(3))*10)/um2pixel, ...
                bregma_offset_y + ((outline(:,1)-bregma(1))*10)/um2pixel,'color',color), ...
                cortical_area_boundaries{curr_area_idx},'uni',false);
        end
        
    case 'ccf_aligned'
        % Plot CCF borders aligned to master retinotopy     
        load([alignment_path filesep 'cortical_area_boundaries_aligned.mat']);

        if exist('repeats','var') && ~isempty(repeats) && numel(repeats) == 4
            % Plot multiple in grid
            h = cell(repeats(3),repeats(4));
            for curr_y_rep = 1:repeats(3)
                y_add = repeats(1)*(curr_y_rep-1);
                for curr_x_rep = 1:repeats(4)
                    x_add = repeats(2)*(curr_x_rep-1);
                    h{repeats(3)*repeats(4)} = ...
                        cellfun(@(areas) cellfun(@(outline) ...
                        plot(outline(:,2)+x_add,outline(:,1)+y_add,'color',color),areas,'uni',false), ...
                        cortical_area_boundaries_aligned,'uni',false);
                end
            end
        else
            % Just plot one
            h = cellfun(@(areas) cellfun(@(outline) plot(outline(:,2),outline(:,1),'color',color),areas,'uni',false), ...
                cortical_area_boundaries_aligned,'uni',false);
        end
        
    case 'retinotopy'
        % Plot master retinotopic borders        
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\retinotopy\retinotopic_boundaries.mat');
        h = cellfun(@(outline) plot(outline(:,2),outline(:,1),'color',color),retinotopic_boundaries,'uni',false);
        
    otherwise
        warning(['Invalid reference: ' type]);
        
end


% % ONLY NEEDED ONCE: create top-down cortical boundaries from CCF
%
% Load in the annotated Allen volume and names
% allen_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
% av = readNPY([allen_path filesep 'annotation_volume_10um_by_index.npy']);
% st = loadStructureTree([allen_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean
%
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
%
% % Save borders
% save_fn = [alignment_path filesep 'cortical_area_boundaries'];
% save(save_fn,'cortical_area_boundaries');
% disp(['Saved ' save_fn]);










