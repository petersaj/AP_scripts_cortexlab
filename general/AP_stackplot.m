function h = AP_stackplot(x,t,spacing,zs,color,ylabels)
% h = AP_stackplot(x,t,spacing,zs,color,ylabels)
%
% Plot lines stacked on each other
% x - 2d matrix of lines to plot
% spacing - spacing between lines
% zscore - true/false, zscore traces
% ylabels - can set the ylabels for each line being plotted

if ~exist('t','var') || isempty(t)
    t = 1:size(x,1);
end
if ~exist('spacing','var') || isempty(spacing)
    spacing = 5;
end
if ~exist('zs','var') || isempty(zs)
   zs = false; 
end

% Stack from top to bottom
spacing_add = spacing*[size(x,2):-1:1];
if ~zs
    x_spaced = bsxfun(@plus,x,spacing_add);
elseif zs
    x_spaced = bsxfun(@plus,zscore(x,[],1),spacing_add);
end

if exist('color','var') && ~isempty(color);
    if length(color) == 1
        h = plot(t,x_spaced,'linewidth',2,'color',color);
    elseif size(color,1) == size(x_spaced,2)
        axis; hold on;
        set(gca,'ColorOrder',color);
        h = plot(t,x_spaced,'linewidth',2);
    end
else
    h = plot(t,x_spaced,'linewidth',2);
end

% Set ylabels if desired
if exist('ylabels','var')
    set(gca,'YTick',sort(spacing_add));
    set(gca,'YTickLabel',flipud(ylabels(:)));
end


