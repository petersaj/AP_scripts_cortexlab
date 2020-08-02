function h = AP_errorfill(x,y,ye,color,alpha,plot_mean,linewidth)
% AP_errorfill(x,y,ye,color,alpha,plot_mean,linewidth);
%
% Draw filled polygon as errorbars
%
% y - time x lines
% ye - time x lines x lower/upper (if only 1, symmetrical error)
% color - lines x 3 (RGB)
% alpha - alpha of error fill
% plot_mean = plot the mean as a line

% Define defaults, reshape to column vectors
if isempty(x)
    x = 1:length(y);
end
x = reshape(x,[],1);

if ~exist('color','var') || isempty(color)
    color = lines(size(y,2));
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
    ye_pos = y + ye(:,:,1);
    ye_neg = y + ye(:,:,2);
elseif size(ye,3) > 2
    error('Y error dim 3 needs to be 1 (symm) or 2 (upper/lower)')
end

% Draw plot
hold on;
for curr_line = 1:size(y,2)
    % Error polygon
    fill([x;flipud(x)],[ye_pos(:,curr_line);flipud(ye_neg(:,curr_line))], ...
        color(curr_line,:),'FaceAlpha',alpha,'EdgeColor','none')
    
    % Plot central line
    if plot_mean
        h = plot(x,y(:,curr_line),'color',color(curr_line,:),'linewidth',linewidth);
    end
end





