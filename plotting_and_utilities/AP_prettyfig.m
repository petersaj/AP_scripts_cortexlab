function AP_prettyfig(format_type)
% AP_prettyfig(format_type)
% 
% Set figure properties for presentations and documents
% format_type - ppt/eps (ppt default)

if ~exist('format_type') || isempty(format_type)
    format_type = 'ppt';
end

% Get all axes in figure
% (from axes or tiledlayout creation)
fig_child = allchild(gcf);

fig_child_ax_idx = strcmp(get(fig_child,'type'),'axes');
fig_child_tiledlayout_idx = strcmp(get(fig_child,'type'),'tiledlayout');

fig_ax = [fig_child(fig_child_ax_idx); ...
    allchild(fig_child(fig_child_tiledlayout_idx))];

% Loop through axes
for curr_ax = 1:length(fig_ax)
    try
    % Set axis properties
    set(fig_ax(curr_ax), ...
        'TickLength',[0.015,0.015], ... % larger ticks
        'TickDir','out', ... % ticks outside plot
        'box','off', ... % no bounding box
        'FontSize',14,'FontName','Calibri'); % larger/presentation font
    catch me
        continue
    end
    
    % Remove axis exponent labels
    fig_ax(curr_ax).XAxis.Exponent = 0;
    fig_ax(curr_ax).YAxis.Exponent = 0;

end

% Powerpoint: transparency doesn't copy, so keep color and make opaque
if strcmp(format_type,'ppt')
    for curr_ax = 1:length(fig_ax)
        ax_objs = get(fig_ax(curr_ax),'Children');
        for curr_obj = 1:length(ax_objs)
            try
                old_transparency = get(ax_objs(curr_obj),'FaceAlpha');
                old_color = get(ax_objs(curr_obj),'FaceColor');
                new_color = old_color+((1-old_color)*old_transparency);
                set(ax_objs(curr_obj),'FaceColor',new_color);
                set(ax_objs(curr_obj),'FaceAlpha',1);
            catch me
                continue
            end
        end
    end
end



