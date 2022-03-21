function AP_prettyfig(format_type)
% AP_prettyfig(format_type)
% 
% Set figure properties for presentations and documents
% format_type - ppt/eps (ppt default)

if ~exist('format_type') || isempty(format_type)
    format_type = 'ppt';
end

% Grab current figure
curr_fig = gcf;

% Set figure background white
curr_fig.Color = 'w';

% Get all axes in figure
% (from axes or tiledlayout creation)
fig_child = allchild(curr_fig);
fig_child_type = get(fig_child,'type');

fig_child_ax_idx = strcmp(fig_child_type,'axes');
fig_child_tiledlayout_idx = strcmp(fig_child_type,'tiledlayout');

fig_ax = [fig_child(fig_child_ax_idx); ...
    allchild(fig_child(fig_child_tiledlayout_idx))];

% Loop through axes
for curr_ax = 1:length(fig_ax)
    
    % Set axis properties
    try      
        fig_ax(curr_ax).TickLength = [0.015,0.015]; % larger ticks
        fig_ax(curr_ax).TickDir = 'out'; % ticks outside plot
        fig_ax(curr_ax).Box = 'off'; % no bounding box
        fig_ax(curr_ax).FontName = 'MyriadPro'; % font
        fig_ax(curr_ax).TitleFontSizeMultiplier = 1; % relative title size
        
        tick_fontsize = 10; % axis tick font size
        fig_ax(curr_ax).XAxis.FontSize = tick_fontsize;
        fig_ax(curr_ax).YAxis.FontSize = tick_fontsize;
        
        axis_fontsize = 14; % axis label font size
        fig_ax(curr_ax).FontSize = axis_fontsize;
        fig_ax(curr_ax).XLabel.FontSize = axis_fontsize;
        fig_ax(curr_ax).YLabel.FontSize = axis_fontsize;
        
        % remove axis exponent labels
        fig_ax(curr_ax).XAxis.Exponent = 0; 
        fig_ax(curr_ax).YAxis.Exponent = 0;      
    catch me
    end
    
    % Turn off legend boxes
    if strcmp(get(fig_ax(curr_ax),'type'),'legend')
        set(fig_ax(curr_ax),'box','off');
    end

    % Powerpoint: transparency doesn't copy, so keep color and make opaque
    if strcmp(format_type,'ppt')
        ax_objs = get(fig_ax(curr_ax),'Children');
        for curr_obj = 1:length(ax_objs)
            try
                old_transparency = get(ax_objs(curr_obj),'FaceAlpha');
                old_color = get(ax_objs(curr_obj),'FaceColor');
                new_color = old_color+((1-old_color)*old_transparency);
                set(ax_objs(curr_obj),'FaceColor',new_color);
                set(ax_objs(curr_obj),'FaceAlpha',1);
            catch me
            end
        end
    end
end



