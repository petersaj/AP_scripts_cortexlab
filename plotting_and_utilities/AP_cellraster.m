function AP_cellraster(align_times,align_groups,unit_sort)
% AP_cellraster(align_times,align_groups,unit_sort)
%
% Raster viewer for Neuropixels
%
% align_times - vector or cell array of times
% align_groups - categorical vectors of trial types (optional)
% unit_sort - sorting for scrolling through units (default = depth)
%
% These variables are required in the base workspace: 
% templates (output from kilosort)
% channel_positions (output from kilosort)
% template_depths (calculated from center-of-mass or otherwise)
% spike_times_timeline (spike_times aligned to timeline)
% spike_templates (output from kilosort)
% template_amplitudes (output from kilosort)
%
% Controls: 
% up/down - switch between units (clicking on unit also selects)
% left/right - switch between alignments (if multiple)
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
if ~exist('align_groups','var') || isempty('align_groups')
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

% Sort the units by depth if not specified
if ~exist('unit_sort','var')
   [~,unit_sort] = sort(template_depths); 
end

% Initialize figure and axes
cellraster_gui = figure('color','w');

% (plot unit depths by depth and relative number of spikes)
unit_axes = subplot(5,5,[1:5:20],'YDir','reverse');
hold on;

norm_spike_n = mat2gray(log(accumarray(spike_templates,1)+1));
unit_dots = plot(norm_spike_n,template_depths,'.k','MarkerSize',20,'ButtonDownFcn',@unit_click);
curr_unit_dots = plot(0,0,'.r','MarkerSize',20);
multiunit_lines = arrayfun(@(x) line(xlim,[0,0],'linewidth',2,'visible','off'),1:2);
xlim(unit_axes,[-0.1,1]);
ylim([-50, max(channel_positions(:,2))+50]);
ylabel('Depth (\mum)')
xlabel('Normalized log rate')

% (plot of waveform across the probe)
waveform_axes = subplot(5,5,[2:5:20],'visible','off','YDir','reverse');
hold on;
ylim([-50, max(channel_positions(:,2))+50]);
waveform_lines = arrayfun(@(x) plot(waveform_axes,0,0,'k','linewidth',1),1:size(templates,3));

linkaxes([unit_axes,waveform_axes],'y');

% (smoothed psth)
psth_axes = subplot(5,5,[3,4,5],'YAxisLocation','right');
hold on;
max_n_groups = max(cell2mat(cellfun(@(x) 1+sum(diff(sort(x,1),[],1) ~= 0),align_groups,'uni',false)));
if max_n_groups == 1
    psth_colors = [0,0,0];
else
    psth_colors = [linspace(0,0.8,max_n_groups)',zeros(max_n_groups,1),zeros(max_n_groups,1)];
end

psth_lines = arrayfun(@(x) plot(NaN,NaN,'linewidth',2,'color',psth_colors(x,:)),1:max_n_groups);
xlabel('Time from event (s)');
ylabel('Spikes/s/trial');

% (raster)
raster_axes = subplot(5,5,[8,9,10,13,14,15,18,19,20],'YDir','reverse','YAxisLocation','right');
hold on;
raster_dots = scatter(NaN,NaN,5,'k','filled');
raster_image = imagesc(NaN,'visible','off'); colormap(raster_axes,hot);
xlabel('Time from event (s)');
ylabel('Trial');

% (spike amplitude across the recording)
amplitude_axes = subplot(5,5,21:25); hold on;
amplitude_plot = plot(NaN,NaN,'.k');
amplitude_lines = arrayfun(@(x) line([0,0],ylim,'linewidth',2),1:2);
xlabel('Experiment time (s)');
ylabel('Template amplitude');
axis tight

% Set default raster times
raster_window = [-0.5,1];
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
gui_data.curr_unit_dots = curr_unit_dots;
gui_data.multiunit_lines = multiunit_lines;
gui_data.waveform_lines = waveform_lines;
gui_data.psth_lines = psth_lines;
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
gui_data.unit_sort = unit_sort;

% (spike data)
gui_data.templates = templates;
gui_data.channel_positions = channel_positions;
gui_data.spike_times = spike_times;
gui_data.spike_templates = spike_templates;
gui_data.template_amplitudes = template_amplitudes;

% (current settings)
gui_data.curr_unit = 1;
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
bin_t = mean(diff(gui_data.t));

curr_psth = grpstats(curr_raster, ...
    gui_data.align_groups{gui_data.curr_align}(:,gui_data.curr_group),@(x) mean(x,1));
curr_smoothed_psth = conv2(padarray(curr_psth, ...
    [0,floor(length(smWin)/2)],'replicate','both'), ...
    smWin,'valid')./bin_t;

% (set the first n lines by group, set all others to NaN)
arrayfun(@(align_groups) set(gui_data.psth_lines(align_groups), ...
    'XData',gui_data.t,'YData',curr_smoothed_psth(align_groups,:)), ...
    1:size(curr_psth,1));
arrayfun(@(align_groups) set(gui_data.psth_lines(align_groups), ...
    'XData',NaN,'YData',NaN), ...
    size(curr_psth,1)+1:length(gui_data.psth_lines));

ylim(get(gui_data.psth_lines(1),'Parent'),[min(curr_smoothed_psth(:)), ...
    max(max(curr_smoothed_psth(:),min(curr_smoothed_psth(:))+1))]);
if length(gui_data.curr_unit) == 1
title(get(gui_data.psth_lines(1),'Parent'), ...
    ['Unit ' num2str(gui_data.curr_unit) ...
    ', Align ' num2str(gui_data.curr_align) ...
    ', Group ' num2str(gui_data.curr_group)],'FontSize',14);
elseif length(gui_data.curr_unit) > 1
title(get(gui_data.psth_lines(1),'Parent'), ...
    ['Multiunit, Align ' num2str(gui_data.curr_align)],'FontSize',14);
end

% Plot raster
if length(gui_data.curr_unit) == 1
    % (single unit mode)
    [raster_y,raster_x] = find(curr_raster);
    set(gui_data.raster_dots,'XData',gui_data.t(raster_x),'YData',raster_y);
    xlim(get(gui_data.raster_dots,'Parent'),[gui_data.t_bins(1),gui_data.t_bins(end)]);
    ylim(get(gui_data.raster_dots,'Parent'),[0,size(gui_data.t_peri_event,1)]);
    % (set dot color by group)
    [~,~,row_group] = unique(gui_data.align_groups{gui_data.curr_align}(:,gui_data.curr_group),'sorted');
    psth_colors = get(gui_data.psth_lines,'color');
    if iscell(psth_colors); psth_colors = cell2mat(psth_colors); end
    raster_dot_color = psth_colors(row_group(raster_y),:);
    set(gui_data.raster_dots,'CData',raster_dot_color);
    
elseif length(gui_data.curr_unit) > 1
    % (multiunit mode)
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

switch eventdata.Key
    case 'downarrow'
        % Next unit
        curr_unit_idx = gui_data.curr_unit(1) == gui_data.unit_sort;
        new_unit = gui_data.unit_sort(circshift(curr_unit_idx,1));
        gui_data.curr_unit = new_unit;
        
    case 'uparrow'
        % Previous unit
        curr_unit_idx = gui_data.curr_unit(end) == gui_data.unit_sort;
        new_unit = gui_data.unit_sort(circshift(curr_unit_idx,-1));
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

    case 'pagedown'
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
        set(gui_data.multiunit_lines(1),'visible','on','YData',repmat(multiunit_top,1,2));
        [~,multiunit_bottom] = ginput(1);
        set(gui_data.multiunit_lines(2),'visible','on','YData',repmat(multiunit_bottom,1,2));
        
        template_depths = get(gui_data.unit_dots,'YData');
        gui_data.curr_unit = find(template_depths >= multiunit_top & ...
            template_depths <= multiunit_bottom);       
        
    case 'u'
        % Enter and go to unit
        new_unit = str2num(cell2mat(inputdlg('Go to unit:')));
        if ~ismember(new_unit,unique(gui_data.spike_templates))
            error(['Unit ' num2str(new_unit) ' not present'])
        end
        gui_data.curr_unit = new_unit;
        
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




















