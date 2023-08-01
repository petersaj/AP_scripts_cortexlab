function h = AP_errorfill(x,y,ye,color,alpha,plot_mean,linewidth)
% AP_errorfill(x,y,ye,color,alpha,plot_mean,linewidth);
%
% Draw filled polygon as errorbars
%
% y - time x lines
% ye - time x lines x lower/upper (+/- error if length 1, values if length 2)
% color - lines x 3 (RGB)
% alpha - alpha of error fill
% plot_mean = plot the mean as a line

% If no y (just plotting error bars), default to NaNs size ye
if isempty(y)
    y = nan(size(ye,1),1);
end

% Define defaults, reshape to column vectors
if isempty(x)
    x = 1:length(y);
end
x = reshape(x,[],1);

% (allow single dimension y's in any orientation)
if sum(size(y) > 1) == 1
    y = reshape(y,[],1);
    if size(ye,2) == length(y)
        ye = ye';
    end
    % (if ye as +/- on dim 2, put on dim 3)
    if size(ye,2) > 1
        ye = permute(ye,[1,3,2]);
    end
end

% If no color set, continue through axis color order
% (after max number of any object type already on plot)
if ~exist('color','var') || isempty(color)
    ax_color_order = get(gca,'ColorOrder');
    ax_children = get(gca,'Children');
    if ~isempty(ax_children)
        ax_child_type = get(ax_children,'type');
        child_type_max_n = max(accumarray(findgroups(ax_child_type),1));
    else
        child_type_max_n = 0;
    end
    color_idx = mod(child_type_max_n+[1:size(y,2)]-1,size(ax_color_order,1))+1;
    color = ax_color_order(color_idx,:);
elseif size(color,1) == 1 && size(y,2) > 1
    % If only one color specified, replicate for every line
    color = repmat(color,size(y,2),1);
end

if ~exist('alpha','var') || isempty(alpha)
    alpha = 0.5;
end

if ~exist('plot_mean','var') || isempty(plot_mean)
    plot_mean = true;
end

if ~exist('linewidth','var') || isempty(linewidth)
    linewidth = 2;
end

% Set error values
if size(ye,3) == 1
    ye_pos = y + ye;
    ye_neg = y - ye;
elseif size(ye,3) == 2
    ye_pos = ye(:,:,1);
    ye_neg = ye(:,:,2);
elseif size(ye,3) > 2
    error('Y error dim 3 needs to be 1 (symm) or 2 (upper/lower)')
end

% Check no mismatching NaNs in error values
if any(any(isnan(ye_pos) ~= isnan(ye_neg)))
    error('Mismatching NaNs in +/- error values')
end

% Draw plot
hold on;
h = gobjects(size(y,2),1);
for curr_line = 1:size(y,2)
    % Error polygon (separate shapes if NaN breaks)
    if ~any(isnan(ye_pos(:,curr_line)))
        h_fill = fill([x;flipud(x)], ...
            [ye_pos(:,curr_line);flipud(ye_neg(:,curr_line))], ...
            color(curr_line,:),'FaceAlpha',alpha,'EdgeColor','none');
    else
        fill_starts = [find(diff([true;isnan(ye_pos(:,curr_line))]) == -1)];
        fill_stops = [find(diff(isnan(ye_pos(:,curr_line))) == 1); ...
            find(~isnan(ye_pos(:,curr_line)),1,'last')];
        for curr_fill = 1:length(fill_starts)
            curr_idx = fill_starts(curr_fill):fill_stops(curr_fill);
            h_fill = fill([x(curr_idx);flipud(x(curr_idx))], ...
                [ye_pos(curr_idx,curr_line);flipud(ye_neg(curr_idx,curr_line))], ...
                color(curr_line,:),'FaceAlpha',alpha,'EdgeColor','none');
        end
    end

    % Plot central line
    if plot_mean
        h(curr_line) = plot(x,y(:,curr_line), ...
            'color',color(curr_line,:),'linewidth',linewidth);
    else
        % (if no central line, set the fill as the handle)
        h(curr_line) = h_fill(1);
    end
end





