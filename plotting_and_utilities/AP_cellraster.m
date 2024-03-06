function AP_cellraster(align_times,align_groups)
% AP_cellraster(align_times,align_groups)
%
% Raster viewer for Neuropixels
%
% align_times - vector or cell array of times
% align_groups - categorical vectors of trial types (optional)
%
% These variables are required in the base workspace: 
% templates (output from kilosort)
% channel_positions (output from kilosort)
% template_depths (calculated from center-of-mass or otherwise)
% spike_times_timeline (spike_times aligned to timeline)
% spike_templates (output from kilosort+1)
% template_amplitudes (output from kilosort)
%
% Controls: 
% up/down - switch between units (clicking on unit also selects)
% left/right - switch between alignments
% pageup/pagedown - switch between trial groupings
% m - select depth range to plot multiunit
% u - go to unit number


% Initiate align_times

% (align times required as input)
if ~exist('align_times','var')
    error('No align times');
end

% (put align_times into cell array if it isn't already, standardize dim)
if ~iscell(align_times)
    align_times = {align_times};
end
align_times = reshape(cellfun(@(x) reshape(x,[],1),align_times,'uni',false),1,[]);

% Initiate align_groups

% (if no align groups specified, create one group for each alignment)
if ~exist('align_groups','var') || isempty(align_groups)
   align_groups =  cellfun(@(x) ones(size(x,1),1),align_times,'uni',false);
end

% (put groups into cell array if it isn't already)
if ~iscell(align_groups)
    align_groups = {align_groups};
end
align_groups = reshape(align_groups,1,[]);

% (replicate for each align_times if only one)
if length(align_times) > 1 && length(align_groups) == 1
   align_groups = repmat(align_groups,size(align_times));
elseif length(align_times) > 1 && length(align_groups) < length(align_times)
    error('Mismatching align time/group sets')
end

% (check group against time dimensions, orient align times x groups)
group_dim = cellfun(@(align,group) find(ismember(size(group),length(align))),align_times,align_groups,'uni',false);
if any(cellfun(@isempty,group_dim))
    error('Mismatching times/groups within align set')
end
align_groups = cellfun(@(groups,dim) shiftdim(groups,dim-1),align_groups,group_dim,'uni',false);

% (if there isn't an all ones category first, make one)
align_groups = cellfun(@(x) padarray(x,[0,1-all(x(:,1) == 1)],1,'pre'),align_groups,'uni',false);

% Pull standard ephys variables from base workspace
try
templates = evalin('base','templates');
channel_positions = evalin('base','channel_positions');
template_depths = evalin('base','template_depths');
spike_times = evalin('base','spike_times_timeline');
spike_templates = evalin('base','spike_templates');
template_amplitudes = evalin('base','template_amplitudes');
catch me
    missing_var = textscan(me.message,'Undefined function or variable '' %s');
    error(['Ephys variable missing from base workspace: ' cell2mat(missing_var{:})]);
end

% Pull regions from base workspace, if available
try
    probe_areas = evalin('base','probe_areas{1}');
catch me
end

% Initialize figure and axes
cellraster_gui = figure('color','w','units','normalized','position',[0.2,0.1,0.5,0.8]);
tiledlayout(cellraster_gui,6,4,'TileSpacing','tight');

% (unit dots: plot depths vs spike number, color background by area)
unit_axes = nexttile([5,1]);
set(unit_axes,'YDir','reverse');
hold on;

if exist('probe_areas','var')
    try
    probe_areas_rgb = permute(cell2mat(cellfun(@(x) hex2dec({x(1:2),x(3:4),x(5:6)})'./255, ...
        probe_areas.color_hex_triplet,'uni',false)),[1,3,2]);

    % Get area depth on probe (update this eventually)
    if any(strcmp('probe_depth',probe_areas.Properties.VariableNames))
        % Old: "probe_depth" relative to top to bank 0
    probe_areas_boundaries = probe_areas.probe_depth;
    elseif any(strcmp('probe_tip_distance',probe_areas.Properties.VariableNames))
        % New: "probe_tip_distance" relative to probe tip
        probe_areas_boundaries = 3.84-probe_areas.probe_tip_distance;
    end
    probe_areas_centers = mean(probe_areas_boundaries,2);

    probe_areas_image_depth = 0:1:max(probe_areas_boundaries,[],'all');
    probe_areas_image_idx = interp1(probe_areas_boundaries(:,1), ...
        1:height(probe_areas),probe_areas_image_depth, ...
        'previous','extrap');
    probe_areas_image = probe_areas_rgb(probe_areas_image_idx,:,:);

    image(unit_axes,[0,1],probe_areas_image_depth,probe_areas_image);
    yline(unique(probe_areas_boundaries(:)),'color','k','linewidth',1);
    set(unit_axes,'YTick',probe_areas_centers,'YTickLabels',probe_areas.acronym);
    catch me
        warning('Probe areas not found: old format?');
    end
end

% (plot unit depths by depth and relative number of spikes)
norm_spike_n = mat2gray(log10(accumarray(findgroups(spike_templates),1)+1));
unit_dots = scatter3( ...
    norm_spike_n,template_depths(unique(spike_templates)), ...
    unique(spike_templates),20,'k','filled','ButtonDownFcn',@unit_click);
xlim(unit_axes,[-0.1,1]);
ylim([-50, max(channel_positions(:,2))+50]);
ylabel('Depth (\mum)')
xlabel('Normalized log rate')

% (plot of waveform across the probe)
waveform_axes =  nexttile([5,1]);
set(waveform_axes,'visible','off','YDir','reverse');
hold on;
ylim([-50, max(channel_positions(:,2))+50]);
waveform_lines = plot(waveform_axes,0,0,'k','linewidth',1);

% (smoothed psth)
psth_axes = nexttile([1,2]);
hold on;
set(psth_axes,'YAxisLocation','right');
xlabel('Time from event (s)');
ylabel('Spikes/s/trial');

% (raster)
raster_axes = nexttile([4,2]);
set(raster_axes,'YDir','reverse','YAxisLocation','right');
hold on;
raster_dots = scatter(NaN,NaN,5,'k','filled');
raster_image = imagesc(NaN,'visible','off'); colormap(raster_axes,hot);
xlabel('Time from event (s)');
ylabel('Trial');

% (spike amplitude across the recording)
amplitude_axes = nexttile([1,4]); hold on;
amplitude_plot = plot(NaN,NaN,'.k');
amplitude_lines = xline(amplitude_axes,[0,0],'linewidth',2,'color','r');
xlabel('Experiment time (s)');
ylabel('Template amplitude');
axis tight

% Set default raster times
raster_window = [-0.5,2];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t = t_bins(1:end-1) + diff(t_bins)./2;
use_align = reshape(align_times{1},[],1);
t_peri_event = use_align + t_bins;
% (handle NaNs by setting rows with NaN times to 0)
t_peri_event(any(isnan(t_peri_event),2),:) = 0;

% Set functions for key presses
set(cellraster_gui,'KeyPressFcn',@key_press);

% Package gui data
gui_data = struct;

% (plots)
gui_data.unit_dots = unit_dots;
gui_data.waveform_lines = waveform_lines;
gui_data.psth_axes = psth_axes;
gui_data.raster_axes = raster_axes;
gui_data.raster_dots = raster_dots;
gui_data.raster_image = raster_image;
gui_data.amplitude_plot = amplitude_plot;
gui_data.amplitude_lines = amplitude_lines; 

% (raster times)
gui_data.t = t;
gui_data.t_bins = t_bins;
gui_data.t_peri_event = t_peri_event;

% (user inputs)
gui_data.align_times = align_times;
gui_data.align_groups = align_groups;

% (spike data)
gui_data.templates = templates;
gui_data.channel_positions = channel_positions;
gui_data.spike_times = spike_times;
gui_data.spike_templates = spike_templates;
gui_data.template_amplitudes = template_amplitudes;

% (current settings)
gui_data.curr_unit = spike_templates(1);
gui_data.curr_align = 1;
gui_data.curr_group = 1;

% Upload gui data and draw
guidata(cellraster_gui, gui_data);
update_plot(cellraster_gui);

end


function update_plot(cellraster_gui,eventdata)

% Get guidata
gui_data = guidata(cellraster_gui);

% Turn on/off the appropriate graphics
if length(gui_data.curr_unit) == 1
    set(gui_data.raster_dots,'visible','on');
    set(gui_data.raster_image,'visible','off');
elseif length(gui_data.curr_unit) > 1
    set(gui_data.raster_dots,'visible','off');
    set(gui_data.raster_image,'visible','on');
end

% Set color and size on unit dots
curr_unit_unique_idx = ismember(unique(gui_data.spike_templates),gui_data.curr_unit);
gui_data.unit_dots.CData = zeros(length(curr_unit_unique_idx),3) + [1,0,0].*curr_unit_unique_idx;
gui_data.unit_dots.SizeData = 20*ones(size(curr_unit_unique_idx)) + 50*curr_unit_unique_idx;

% Plot waveform across probe (reversed YDir, flip Y axis and plot depth)
template_xscale = 3;
template_yscale = 0.01;
template_n_surround_plot = 5; % plot N channels around max

[~,max_channel] = max(max(abs(gui_data.templates(gui_data.curr_unit,:,:)),[],2),[],3);
center_channel = round(median(max_channel));
template_plot_channels = [-template_n_surround_plot:template_n_surround_plot] + round(median(center_channel));
template_plot_channels = template_plot_channels(template_plot_channels > 0 & ...
    template_plot_channels <= size(gui_data.channel_positions,1));

template_y = permute(-gui_data.templates(gui_data.curr_unit,:,template_plot_channels) ...
    *template_yscale + permute(gui_data.channel_positions(template_plot_channels,2),[3,2,1]),[2,3,1]);
template_x = repmat((1:size(gui_data.templates,2))' + ...
    gui_data.channel_positions(template_plot_channels,1)'*template_xscale,[1,1,length(gui_data.curr_unit)]);
set(gui_data.waveform_lines, ...
    'XData',reshape(padarray(template_x,[1,0],NaN,'post'),[],1), ...
    'YData',reshape(padarray(template_y,[1,0],NaN,'post'),[],1));

yrange = range(gui_data.channel_positions(:,2))*0.03.*[-1,1];
ylim(get(gui_data.waveform_lines,'Parent'),[gui_data.channel_positions(center_channel,2) + yrange]);

set(gui_data.waveform_lines, ...
    'XData',reshape(padarray(template_x,[1,0],NaN,'post'),[],1), ...
    'YData',reshape(padarray(template_y,[1,0],NaN,'post'),[],1));
axis(get(gui_data.waveform_lines(1),'Parent'),'tight')

% Bin spikes (use only spikes within time range, big speed-up)
curr_spikes_idx = ismember(gui_data.spike_templates,gui_data.curr_unit);
curr_raster_spike_times = gui_data.spike_times(curr_spikes_idx);
curr_raster_spike_times(curr_raster_spike_times < min(gui_data.t_peri_event(:)) | ...
    curr_raster_spike_times > max(gui_data.t_peri_event(:))) = [];

if ~any(diff(reshape(gui_data.t_peri_event',[],1)) < 0)
    % (if no backward time jumps, can do long bin and cut out in-between, faster)
    curr_raster_continuous = reshape([histcounts(curr_raster_spike_times, ...
        reshape(gui_data.t_peri_event',[],1)),NaN],size(gui_data.t_peri_event'))';
    curr_raster = curr_raster_continuous(:,1:end-1);   
else
    % (otherwise, bin trial-by-trial)
    curr_raster = cell2mat(arrayfun(@(x) ...
        histcounts(curr_raster_spike_times,gui_data.t_peri_event(x,:)), ...
        [1:size(gui_data.t_peri_event,1)]','uni',false));
end

% Set color scheme
curr_group = gui_data.align_groups{gui_data.curr_align}(:,gui_data.curr_group);
if length(unique(curr_group)) == 1
    % Black if one group or number
    group_colors = [0,0,0];
elseif length(unique(curr_group)) == length(curr_group)
    % Black if number of groups = trials (e.g. sort index)
    % (and set grouping  as one's and 
    group_colors = [0,0,0];
    trial_sort = gui_data.align_groups{gui_data.curr_align}(:,gui_data.curr_group);
    curr_group = ones(size(curr_group));
elseif length(unique(sign(curr_group(curr_group ~= 0)))) == 1
    % 'Lines' colors if all groups positive
    n_groups = length(unique(curr_group));
    group_colors = lines(n_groups);
elseif length(unique(sign(curr_group(curr_group ~= 0)))) == 2
    % Symmetrical blue-black-red if negative and positive groups
    n_groups_pos = length(unique(curr_group(curr_group > 0)));
    group_colors_pos = [linspace(0.3,1,n_groups_pos)',zeros(n_groups_pos,1),zeros(n_groups_pos,1)];
    
    n_groups_neg = length(unique(curr_group(curr_group < 0)));
    group_colors_neg = [zeros(n_groups_neg,1),zeros(n_groups_neg,1),linspace(0.3,1,n_groups_neg)'];
    
    n_groups_zero = length(unique(curr_group(curr_group == 0)));
    group_colors_zero = [zeros(n_groups_zero,1),zeros(n_groups_zero,1),zeros(n_groups_zero,1)];
    
    group_colors = [flipud(group_colors_neg);group_colors_zero;group_colors_pos];    
end

% Plot smoothed PSTH
smooth_size = 51;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
bin_t = mean(diff(gui_data.t));

curr_psth = grpstats(curr_raster,curr_group,@(x) mean(x,1));
curr_smoothed_psth = conv2(padarray(curr_psth, ...
    [0,floor(length(smWin)/2)],'replicate','both'), ...
    smWin,'valid')./bin_t;

cla(gui_data.psth_axes);
set(gui_data.psth_axes,'ColorOrder',group_colors);
plot(gui_data.psth_axes,gui_data.t,curr_smoothed_psth','linewidth',2);

ylim(gui_data.psth_axes,[min(curr_smoothed_psth(:)), ...
    max(max(curr_smoothed_psth(:),min(curr_smoothed_psth(:))+1))]);
if length(gui_data.curr_unit) == 1
title(gui_data.psth_axes, ...
    ['Unit ' num2str(gui_data.curr_unit) ...
    ', Align ' num2str(gui_data.curr_align) ...
    ', Group ' num2str(gui_data.curr_group)],'FontSize',14);
elseif length(gui_data.curr_unit) > 1
title(gui_data.psth_axes, ...
    ['Multiunit, Align ' num2str(gui_data.curr_align)],'FontSize',14);
end

% Plot raster
% (sort raster by group)
if ~exist('trial_sort','var')
    [~,trial_sort] = sort(curr_group);
end
curr_raster_sorted = curr_raster(trial_sort,:);

if all(trial_sort == [1:length(trial_sort)]')
    ylabel(gui_data.raster_axes,'Trial')
else
    ylabel(gui_data.raster_axes,'Sorted trial')
end
if length(gui_data.curr_unit) == 1
    % (single unit mode)

    % (plot raster matrix as x,y)
    [raster_y,raster_x] = find(curr_raster_sorted);
    set(gui_data.raster_dots,'XData',gui_data.t(raster_x),'YData',raster_y);
    xlim(get(gui_data.raster_dots,'Parent'),[gui_data.t_bins(1),gui_data.t_bins(end)]);
    ylim(get(gui_data.raster_dots,'Parent'),[0,size(gui_data.t_peri_event,1)]);
    
    % (set dot color by group)
    [~,~,row_group] = unique(curr_group(trial_sort),'sorted');
    raster_dot_color = group_colors(row_group(raster_y),:);
    set(gui_data.raster_dots,'CData',raster_dot_color);
    
elseif length(gui_data.curr_unit) > 1
    % (multiunit mode)
     
    % (plot raster matrix as smoothed heatmap)
    raster_heatmap = imgaussfilt(curr_raster_sorted,[3,5]);
    set(gui_data.raster_image,'XData',gui_data.t,'YData', ...
        1:size(gui_data.t_peri_event,1),'CData',raster_heatmap);
    caxis(get(gui_data.raster_image,'Parent'),prctile(raster_heatmap(:),[0,100]));
    axis(get(gui_data.raster_image,'Parent'),'tight');
end

% Plot template amplitude over whole experiment
if length(gui_data.curr_unit) == 1
    set(gui_data.amplitude_plot,'XData', ...
        gui_data.spike_times(curr_spikes_idx), ...
        'YData',gui_data.template_amplitudes(curr_spikes_idx),'linestyle','none');
elseif length(gui_data.curr_unit) > 1
    long_bin_size = 60;
    long_bins = gui_data.spike_times(1):long_bin_size:gui_data.spike_times(end);
    long_bins_t = long_bins(1:end-1) + diff(long_bins)/2;
    long_spikes_binned = discretize(gui_data.spike_times,long_bins);
    amplitude_binned = accumarray(long_spikes_binned(curr_spikes_idx & ~isnan(long_spikes_binned)), ...
        gui_data.template_amplitudes(curr_spikes_idx & ~isnan(long_spikes_binned)),size(long_bins_t'),@nansum,NaN);
    set(gui_data.amplitude_plot,'XData',long_bins_t,'YData',amplitude_binned,'linestyle','-');
end

[gui_data.amplitude_lines.Value] = deal(min(gui_data.t_peri_event(:)),max(gui_data.t_peri_event(:)));

end


function key_press(cellraster_gui,eventdata)

% Get guidata
gui_data = guidata(cellraster_gui);

switch eventdata.Key
    case 'downarrow'
        % One unit down      
        template_depths = get(gui_data.unit_dots,'YData');
        template_id = get(gui_data.unit_dots,'ZData');
        depth_sort = sortrows([template_depths',template_id'],1);
        new_unit = depth_sort(circshift(ismember(depth_sort(:,2),gui_data.curr_unit),1),2);
        gui_data.curr_unit = new_unit;
        
    case 'uparrow'
        % One unit up
        template_depths = get(gui_data.unit_dots,'YData');
        template_id = get(gui_data.unit_dots,'ZData');
        depth_sort = sortrows([template_depths',template_id'],1);
        new_unit = depth_sort(circshift(ismember(depth_sort(:,2),gui_data.curr_unit),-1),2);
        gui_data.curr_unit = new_unit;
        
    case 'rightarrow'
        % Next alignment
        new_align = gui_data.curr_align + 1;
        if new_align > length(gui_data.align_times)
            new_align = 1;
        end
        use_align = reshape(gui_data.align_times{new_align},[],1);
        t_peri_event = use_align + gui_data.t_bins;
        
        % (handle NaNs by setting rows with NaN times to 0)
        t_peri_event(any(isnan(t_peri_event),2),:) = 0;
        
        gui_data.curr_align = new_align;
        gui_data.t_peri_event = t_peri_event;
        gui_data.curr_group = 1;
        
    case 'leftarrow'
        % Previous alignment
        new_align = gui_data.curr_align - 1;
        if new_align < 1
            new_align = length(gui_data.align_times);
        end
        use_align = reshape(gui_data.align_times{new_align},[],1);
        t_peri_event = use_align + gui_data.t_bins;
        
        % (handle NaNs by setting rows with NaN times to 0)
        t_peri_event(any(isnan(t_peri_event),2),:) = 0;
        
        gui_data.curr_align = new_align;
        gui_data.t_peri_event = t_peri_event;
        gui_data.curr_group = 1;

    case {'pagedown','space'}
        % Next group
        next_group = gui_data.curr_group + 1;
        if next_group > size(gui_data.align_groups{gui_data.curr_align},2)
            next_group = 1;
        end
        gui_data.curr_group = next_group;
        
    case 'pageup'
        % Previous group
        next_group = gui_data.curr_group - 1;
        if next_group < 1
            next_group = size(gui_data.align_groups{gui_data.curr_align},2);
        end
        gui_data.curr_group = next_group;
                
    case 'm'
        % Multiunit (select on unit depth plot)
        [~,multiunit_top] = ginput(1);
        [~,multiunit_bottom] = ginput(1);
        
        template_depths = get(gui_data.unit_dots,'YData');
        template_id = get(gui_data.unit_dots,'ZData');
        gui_data.curr_unit = template_id(template_depths >= multiunit_top & ...
            template_depths <= multiunit_bottom);       
        
    case 'u'
        % Enter and go to unit
        new_unit = str2num(cell2mat(inputdlg('Go to unit:')));
        if ~ismember(new_unit,unique(gui_data.spike_templates))
            error(['Unit ' num2str(new_unit) ' not present'])
        else
            gui_data.curr_unit = new_unit;
        end

    case 't'
        % Change time
        raster_window = str2num(cell2mat(inputdlg('Peri-event times:')));

        % Create new bin times
        psth_bin_size = 0.001;
        t_bins = raster_window(1):psth_bin_size:raster_window(2);
        t = t_bins(1:end-1) + diff(t_bins)./2;

        use_align = reshape(gui_data.align_times{gui_data.curr_align},[],1);
        t_peri_event = use_align + t_bins;

        % (handle NaNs by setting rows with NaN times to 0)
        t_peri_event(any(isnan(t_peri_event),2),:) = 0;

        % Set new bin times
        gui_data.t = t;
        gui_data.t_bins = t_bins;
        gui_data.t_peri_event = t_peri_event;

end

% Upload gui data and draw
guidata(cellraster_gui,gui_data);
update_plot(cellraster_gui);
drawnow;
        
end

function unit_click(cellraster_gui,eventdata)

% Get guidata
gui_data = guidata(cellraster_gui);

% Get the clicked unit, update current unit
unit_x = get(gui_data.unit_dots,'XData');
unit_y = get(gui_data.unit_dots,'YData');

[~,clicked_dot] = min(sqrt(sum(([unit_x;unit_y] - ...
    eventdata.IntersectionPoint(1:2)').^2,1)));

template_id = get(gui_data.unit_dots,'ZData');
clicked_unit = template_id(clicked_dot);

gui_data.curr_unit = clicked_unit;

% Upload gui data and draw
guidata(cellraster_gui,gui_data);
update_plot(cellraster_gui);

end




















