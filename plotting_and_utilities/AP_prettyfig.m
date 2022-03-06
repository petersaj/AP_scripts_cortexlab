function AP_prettyfig
% Set figure properties for presentations and documents

% Get all axes in figure
fig_ax = get(gcf,'Children');

% Loop through axes
for curr_ax = 1:length(fig_ax)
    try
    % Set larger tick length
    set(fig_ax(curr_ax),'TickLength',[0.03,0.03])
    
    % Turn off bounding box
    set(fig_ax(curr_ax),'box','off')
    
    % Set font and font size
    set(fig_ax(curr_ax),'FontSize',14,'FontName','Calibri');
    catch me
        continue
    end        
end

% Keeping this here for now? if alpha < 1 then matlab can't handle, so make
% opaque but keep relative color
for curr_ax = 1:length(fig_ax)
    ax_objs = get(fig_ax(curr_ax),'Children');
    for curr_obj = 1:length(ax_objs)
        try
            old_transparency = get(ax_objs(curr_obj),'FaceAlpha');
            old_color = get(ax_objs(curr_obj),'FaceColor');
            new_color = old_color+((1-old_color)*old_transparency);
            set(ax_objs(curr_obj),'FaceColor',new_color);
            set(ax_objs(curr_obj),'FaceAlpha',1);
            disp('Removed transparency');
        catch me
            continue
        end
    end
end



