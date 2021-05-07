function AP_prettyfig
% Set figure properties for presentations


% Keeping this here for now? if alpha < 1 then matlab can't handle, so make
% opaque but keep relative color


fig_ax = get(gcf,'Children');
for curr_ax = 1:length(fig_ax)
    ax_objs = get(fig_ax(curr_ax),'Children');
    for curr_obj = 1:length(ax_objs)
        try
            old_transparency = get(ax_objs(curr_obj),'FaceAlpha');
            set(ax_objs(curr_obj),'FaceColor', ...
                min(get(ax_objs(curr_obj),'FaceColor') + old_transparency,1));
            set(ax_objs(curr_obj),'FaceAlpha',1);           
        catch me
            continue
        end
    end
end
disp('Removed transparency');