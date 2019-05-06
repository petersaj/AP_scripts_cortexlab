function h = AP_errorfill(x,y,ye,color,alpha,plot_mean)
% AP_errorfill(x,y,ye,color,alpha,plot_mean);
%
% Draw filled polygon as errorbars
% ye = y error, can be 1 or 2 vectors (2 elements is ambiguous)
% plot_mean = plot the mean as a line

% Define defaults, reshape to column vectors
if isempty(x)
    x = 1:length(y);
end

if ~exist('color','var') || isempty(color)
    color = 'k';
end

if ~exist('alpha','var') || isempty(alpha)
    alpha = 0.5;
end

if ~exist('plot_mean','var') || isempty(plot_mean)
    plot_mean = true;
end

x = reshape(x,[],1);
y = reshape(y,[],1);
[d1,d2] = size(ye);
if d2 > d1
    ye = ye';
end

% Set error values
if size(ye,2) == 1
    ye_pos = y + ye;
    ye_neg = y - ye;
elseif size(ye,2) == 2
    ye_pos = y + ye(:,1);
    ye_neg = y + ye(:,2);
else
    error('Y error size wrong');
end

hold on;

% Draw 
fill([x;flipud(x)],[ye_pos;flipud(ye_neg)], ...
    color,'FaceAlpha',alpha,'EdgeColor','none')

% Plot central line
if plot_mean
    h = plot(x,y,'color',color,'linewidth',2);
end






