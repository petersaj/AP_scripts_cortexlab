function AP_errorfill(x,y,ye,color,alpha)
% AP_errorfill(x,y,ye,color,alpha);
%
% Draw filled polygon as errorbars
% ye = y error, can be vector or two-row matrix (first +, second -)

% Define defaults
if isempty(x)
    x = 1:length(y);
end

if nargin < 4 || isempty(color)
    color = 'k';
end

if nargin < 5 || isempty(alpha)
    alpha = 0.5;
end

% Set error values
if size(ye,1) == 1;
    ye_pos = y + ye;
    ye_neg = y - ye;
elseif size(ye,1) == 2
    ye_pos = y + ye(1,:);
    ye_neg = y - ye(2,:);
else
    error('Y error size wrong');
end

hold on;

% Draw 
fill([x,fliplr(x)],[ye_pos,fliplr(ye_neg)], ...
    color,'FaceAlpha',alpha,'EdgeColor','none')

% Plot central line
plot(x,y,'color','k','linewidth',2);






