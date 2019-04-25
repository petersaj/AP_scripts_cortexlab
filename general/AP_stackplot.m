function h = AP_stackplot(x,t,spacing,zs,color,ylabels,zero_lines)
% h = AP_stackplot(x,t,spacing,zs,color,ylabels,zero_lines)
%
% Plot lines stacked on each other
% x - matrix to plot, N points x M lines
% spacing - spacing between lines
% zscore - true/false, zscore traces
% ylabels - can set the ylabels for each line being plotted
% zero_lines - plot horizontal lines at zero for each line

% Set defaults
if ~exist('t','var') || isempty(t)
    t = 1:size(x,1);
end
if ~exist('spacing','var') || isempty(spacing)
    % default spacing is a little larger than the range of the data
    spacing = range(x(:))*1.2;
end
if ~exist('zs','var') || isempty(zs)
   zs = false; 
end
if ~exist('zero_lines','var') || isempty(zero_lines)
    zero_lines = false;
end

% Stack from top to bottom
spacing_add = spacing*[size(x,2):-1:1];
if ~zs
    x_spaced = x + spacing_add;
elseif zs
    % (zscore ignoring nans)
    x_spaced = (x-nanmean(x,1))./nanstd(x,[],1) + spacing_add;
end

if exist('color','var') && ~isempty(color)
    if length(color) == 1
        h = plot(t,x_spaced,'linewidth',2,'color',color);
    elseif size(color,1) == 1 && size(color,2) == 3
        h = plot(t,x_spaced,'linewidth',2,'color',color);
    elseif size(color,1) == size(x_spaced,2)
        axis; hold on;
        set(gca,'ColorOrder',color);
        h = plot(t,x_spaced,'linewidth',2);
    end
else
    h = plot(t,x_spaced,'linewidth',2);
end

% Set ylabels
if exist('ylabels','var')
    set(gca,'YTick',sort(spacing_add));
    set(gca,'YTickLabel',flipud(ylabels(:)));
end

% Plot horizontal zero-lines 
if zero_lines
    for i = 1:size(x,2)
        axis tight
        line(xlim,repmat(spacing_add(i),1,2),'color','k');
    end
end






