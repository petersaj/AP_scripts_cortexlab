function AP_cellraster(align_times,unit_sort,templates,channel_positions,template_depths,spike_times,spike_templates,template_amplitudes)
% Working on: a custom equivalent of psthViewer

% Align times required as input
if ~exist('align_times','var')
    error('Align times required');
end

% Put align_times into cell array if it isn't already
if ~iscell(align_times)
    align_times = {align_times};
end

% For standard things: if not specified, pull from base workspace
if ~exist('templates','var')
    templates = evalin('base','templates');
end

if ~exist('channel_positions','var')
    channel_positions = evalin('base','channel_positions');
end

if ~exist('template_depths','var')
    template_depths = evalin('base','template_depths');
end

if ~exist('spike_times','var')
    spike_times = evalin('base','spike_times_timeline');
end

if ~exist('spike_templates','var')
    spike_templates = evalin('base','spike_templates');
end

if ~exist('template_amplitudes','var')
    template_amplitudes = evalin('base','template_amplitudes');
end

% Sort the units by depth if not specified
if ~exist('unit_sort','var')
   [~,unit_sort] = sort(template_depths); 
end

% Initialize figure and axes
cellraster_gui = figure;

unit_axes = subplot(5,5,[1:5:20],'visible','off','YDir','reverse');
hold on;
unit_dots = plot(rand(size(template_depths)),template_depths,'.k','MarkerSize',20,'ButtonDownFcn',@unit_click);
curr_unit_dots = plot(0,0,'.r','MarkerSize',20);
multiunit_lines = arrayfun(@(x) line(xlim,[0,0],'linewidth',2,'visible','off'),1:2);
xlim(unit_axes,[0,1]);
ylim([-50, max(channel_positions(:,2))+50]);

waveform_axes = subplot(5,5,[2:5:20],'visible','off','YDir','reverse');
hold on;
ylim([-50, max(channel_positions(:,2))+50]);
waveform_lines = arrayfun(@(x) plot(waveform_axes,0,0,'k','linewidth',1),1:size(templates,3));

linkaxes([unit_axes,waveform_axes],'y');

psth_axes = subplot(5,5,[3,4,5],'YAxisLocation','right');
hold on;
psth_line = plot(0,0,'k','linewidth',2);
xlabel('Time from event (s)');
ylabel('Spikes/s/trial');

raster_axes = subplot(5,5,[8,9,10,13,14,15,18,19,20],'YDir','reverse','YAxisLocation','right');
hold on;
raster_dots = plot(0,0,'.k','linewidth',2);
raster_image = imagesc(0,'visible','off'); colormap(raster_axes,hot);
xlabel('Time from event (s)');
ylabel('Trial');

amplitude_axes = subplot(5,5,21:25); hold on;
amplitude_plot = plot(0,0,'.k');
amplitude_lines = arrayfun(@(x) line([0,0],ylim,'linewidth',2),1:2);
xlabel('Experiment time (s)');
ylabel('Template amplitude');
axis tight

% Set default raster times
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t = t_bins(1:end-1) + diff(t_bins)./2;
use_align = reshape(align_times{1}(~isnan(align_times{1})),[],1);
t_peri_event = use_align + t_bins;

% Set functions for key presses
set(cellraster_gui,'KeyPressFcn',@key_press);

% Package gui data
gui_data = struct;

% (plots)
gui_data.unit_dots = unit_dots;
gui_data.curr_unit_dots = curr_unit_dots;
gui_data.multiunit_lines = multiunit_lines;
gui_data.waveform_lines = waveform_lines;
gui_data.psth_line = psth_line;
gui_data.raster_dots = raster_dots;
gui_data.raster_image = raster_image;
gui_data.amplitude_plot = amplitude_plot;
gui_data.amplitude_lines = amplitude_lines; 

% (raster times)
gui_data.t = t;
gui_data.t_bins = t_bins;
gui_data.t_peri_event = t_peri_event;

% (spike data)
gui_data.align_times = align_times;
gui_data.unit_sort = unit_sort;
gui_data.templates = templates;
gui_data.channel_positions = channel_positions;
gui_data.spike_times = spike_times;
gui_data.spike_templates = spike_templates;
gui_data.template_amplitudes = template_amplitudes;

% (current settings)
gui_data.curr_unit = 1;
gui_data.curr_align = 1;

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
    set(gui_data.multiunit_lines,'visible','off');
    set(gui_data.raster_image,'visible','off');
elseif length(gui_data.curr_unit) > 1
    set(gui_data.raster_dots,'visible','off');
    set(gui_data.multiunit_lines,'visible','on');
    set(gui_data.raster_image,'visible','on');
end


% Plot depth location on probe
unit_x = get(gui_data.unit_dots,'XData');
unit_y = get(gui_data.unit_dots,'YData');
set(gui_data.curr_unit_dots,'XData',unit_x(gui_data.curr_unit), ...
    'YData',unit_y(gui_data.curr_unit));

% Plot waveform across probe (reversed YDir, flip Y axis and plot depth)
template_xscale = 7;
template_yscale = 5;

template_y = permute(mean(gui_data.templates(gui_data.curr_unit,:,:),1),[3,2,1]);
template_y = -template_y*template_yscale + gui_data.channel_positions(:,2);
template_x = (1:size(gui_data.templates,2)) + gui_data.channel_positions(:,1)*template_xscale;

template_channel_amp = range(gui_data.templates(gui_data.curr_unit,:,:),2);
template_thresh = max(template_channel_amp,[],3)*0.2;
template_use_channels = any(template_channel_amp > template_thresh,1);
[~,max_channel] = max(max(abs(gui_data.templates(gui_data.curr_unit,:,:)),[],2),[],3);

arrayfun(@(ch) set(gui_data.waveform_lines(ch),'XData',template_x(ch,:),'YData',template_y(ch,:)),1:size(gui_data.templates,3));
arrayfun(@(ch) set(gui_data.waveform_lines(ch),'Color','r'),find(template_use_channels));
arrayfun(@(ch) set(gui_data.waveform_lines(ch),'Color','k'),find(~template_use_channels));
set(gui_data.waveform_lines(max_channel),'Color','b');

% Bin spikes (use only spikes within time range, big speed-up)
curr_spikes_idx = ismember(gui_data.spike_templates,gui_data.curr_unit);
curr_raster_spike_times = gui_data.spike_times(curr_spikes_idx);
curr_raster_spike_times(curr_raster_spike_times < min(gui_data.t_peri_event(:)) | ...
    curr_raster_spike_times > max(gui_data.t_peri_event(:))) = [];

curr_raster = cell2mat(arrayfun(@(x) ...
    histcounts(curr_raster_spike_times,gui_data.t_peri_event(x,:)), ...
    [1:size(gui_data.t_peri_event,1)]','uni',false));

% Plot smoothed PSTH
smooth_size = 51;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
bin_trial_t = mean(diff(gui_data.t)).*size(gui_data.t_peri_event,1);
curr_psth = conv(padarray(sum(curr_raster,1), ...
    [0,floor(length(smWin)/2)],'replicate','both'), ...
    smWin,'valid')./bin_trial_t;

set(gui_data.psth_line,'XData',gui_data.t,'YData',curr_psth);
ylim(get(gui_data.psth_line,'Parent'),[min(curr_psth),max(max(curr_psth,min(curr_psth)+1))]);
if length(gui_data.curr_unit) == 1
title(get(gui_data.psth_line,'Parent'), ...
    ['Unit ' num2str(gui_data.curr_unit) ', Align ' num2str(gui_data.curr_align)],'FontSize',14);
elseif length(gui_data.curr_unit) > 1
title(get(gui_data.psth_line,'Parent'), ...
    ['Multiunit, Align ' num2str(gui_data.curr_align)],'FontSize',14);
end

% Plot raster
if length(gui_data.curr_unit) == 1
    [raster_y,raster_x] = find(curr_raster);
    set(gui_data.raster_dots,'XData',gui_data.t(raster_x),'YData',raster_y);
    xlim(get(gui_data.raster_dots,'Parent'),[gui_data.t_bins(1),gui_data.t_bins(end)]);
    ylim(get(gui_data.raster_dots,'Parent'),[0,size(gui_data.t_peri_event,1)]);
elseif length(gui_data.curr_unit) > 1
    raster_heatmap = imgaussfilt(curr_raster,[5,10]);
    set(gui_data.raster_image,'XData',gui_data.t,'YData', ...
        1:size(gui_data.t_peri_event,1),'CData',raster_heatmap);
    caxis(get(gui_data.raster_image,'Parent'),prctile(raster_heatmap(:),[0.05,99.5]));
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

[ymin,ymax] = bounds(get(gui_data.amplitude_plot,'YData'));
set(gui_data.amplitude_lines(1),'XData',repmat(min(gui_data.t_peri_event(:)),2,1),'YData',[ymin,ymax]);
set(gui_data.amplitude_lines(2),'XData',repmat(max(gui_data.t_peri_event(:)),2,1),'YData',[ymin,ymax]);

end


function key_press(cellraster_gui,eventdata)

% Get guidata
gui_data = guidata(cellraster_gui);

% Get unique units
unique_templates = unique(gui_data.spike_templates);

switch eventdata.Key
    case 'rightarrow'
        % Right = next unit
        curr_unit_idx = gui_data.curr_unit(1) == gui_data.unit_sort;
        new_unit = gui_data.unit_sort(circshift(curr_unit_idx,1));
        gui_data.curr_unit = new_unit;
    case 'leftarrow'
        % Left = previous unit
        curr_unit_idx = gui_data.curr_unit(end) == gui_data.unit_sort;
        new_unit = gui_data.unit_sort(circshift(curr_unit_idx,-1));
        gui_data.curr_unit = new_unit;
    case 'downarrow'
        % Down = next alignment
        next_align = gui_data.curr_align + 1;
        if next_align > length(gui_data.align_times)
            next_align = 1;
        end
        use_align = reshape(gui_data.align_times{next_align}(~isnan(gui_data.align_times{next_align})),[],1);
        t_peri_event = use_align + gui_data.t_bins;
        gui_data.curr_align = next_align;
        gui_data.t_peri_event = t_peri_event;
    case 'uparrow'
        % Up = next alignment
        next_align = gui_data.curr_align - 1;
        if next_align < 1
            next_align = length(gui_data.align_times);
        end
        use_align = reshape(gui_data.align_times{next_align}(~isnan(gui_data.align_times{next_align})),[],1);
        t_peri_event = use_align + gui_data.t_bins;
        gui_data.curr_align = next_align;
        gui_data.t_peri_event = t_peri_event;
    case 'm'
        % m = multiunit (select on unit depth plot)
        [~,multiunit_top] = ginput(1);
        set(gui_data.multiunit_lines(1),'visible','on','YData',repmat(multiunit_top,1,2));
        [~,multiunit_bottom] = ginput(1);
        set(gui_data.multiunit_lines(2),'visible','on','YData',repmat(multiunit_bottom,1,2));
        
        template_depths = get(gui_data.unit_dots,'YData');
        gui_data.curr_unit = find(template_depths >= multiunit_top & ...
            template_depths <= multiunit_bottom);        
end

% Upload gui data and draw
guidata(cellraster_gui,gui_data);
update_plot(cellraster_gui);
        
end

function unit_click(cellraster_gui,eventdata)

% Get guidata
gui_data = guidata(cellraster_gui);

% Get the clicked unit, update current unit
unit_x = get(gui_data.unit_dots,'XData');
unit_y = get(gui_data.unit_dots,'YData');

[~,clicked_unit] = min(sqrt(sum(([unit_x;unit_y] - ...
    eventdata.IntersectionPoint(1:2)').^2,1)));

gui_data.curr_unit = clicked_unit;

% Upload gui data and draw
guidata(cellraster_gui,gui_data);
update_plot(cellraster_gui);

end




















