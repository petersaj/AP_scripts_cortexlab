function AP_stackplot(x,t,spacing,zs,color)
% fig = AP_stackplot(x,t,spacing,zs)
%
% Plot lines stacked on each other
% x - 2d matrix of lines to plot
% spacing - spacing between lines
% zscore - true/false, zscore traces

if ~exist('t','var') || isempty(t)
    t = 1:size(x,1);
end
if ~exist('spacing','var') || isempty(spacing)
    spacing = 5;
end
if ~exist('zs','var') || isempty(zs)
   zs = false; 
end

spacing_add = spacing*[1:size(x,2)];
if ~zs
    x_spaced = bsxfun(@plus,x,spacing_add);
elseif zs
    x_spaced = bsxfun(@plus,zscore(x,[],1),spacing_add);
end

if exist('color','var') && ~isempty(color);
    plot(t,x_spaced,'linewidth',2,'color',color);
else
    plot(t,x_spaced,'linewidth',2);
end


