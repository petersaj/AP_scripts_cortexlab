function h = AP_stackplot(x,t,spacing,yscale,color,ylabels,zero_lines)
% h = AP_stackplot(x,t,spacing,yscale,color,ylabels,zero_lines)
%
% Plot lines stacked on each other
% x - matrix to plot, N points x M lines
% t - time (x-axis)
% spacing - number: spacing between traces, vector: y-position of traces
% yscale - value to rescale y-values
% ylabels - can set the ylabels for each line being plotted
% zero_lines - plot horizontal lines at zero for each line

% Set defaults
if ~exist('t','var') || isempty(t)
    t = 1:size(x,1);
end
if ~exist('spacing','var') || isempty(spacing)
    % default spacing is a little larger than the range of the data
    spacing = range(x(:))*1.2;
elseif length(spacing) > 1 && length(spacing) ~= size(x,2)
    error('Different number of lines and y-values');
else
    % Force row orientation of spacing vector
    spacing = reshape(spacing,1,[]);
end
if ~exist('yscale','var') || isempty(yscale)
   yscale = 1; 
end
if ~exist('zero_lines','var') || isempty(zero_lines)
    zero_lines = false;
end

% Hold the current plot to retain properties
hold on

% Stack from top to bottom
if length(spacing) == 1
    spacing_add = spacing*cast([size(x,2):-1:1],class(spacing));
else
    spacing_add = spacing;
end
x_spaced = x*yscale + spacing_add;

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






