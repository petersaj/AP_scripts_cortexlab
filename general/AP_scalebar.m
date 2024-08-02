function AP_scalebar(x_scale,y_scale)
% AP_scalebar(x_scale,y_scale)
% 
% Draws scalebars in the bottom left corner of a plot 

% Place on current axis
draw_ax = gca;

% Set properties
% (draw magenta for now so they're obvious when transferring)
scalebar_col = 'm';
scalebar_width = 2;

if ~isempty(x_scale)
    line(draw_ax,min(xlim(draw_ax)) + [0,x_scale],repmat(min(ylim(draw_ax)),2,1),'color',scalebar_col,'linewidth',scalebar_width)
end
if ~isempty(y_scale)
    line(draw_ax,repmat(min(xlim(draw_ax)),2,1),min(ylim(draw_ax)) + [0,y_scale],'color',scalebar_col,'linewidth',scalebar_width)
end

