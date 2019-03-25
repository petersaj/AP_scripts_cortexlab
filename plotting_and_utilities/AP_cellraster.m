function AP_cellraster(align_times,templates,channel_positions,spike_times,spike_templates,template_amplitudes)
% Working on: a custom equivalent of psthViewer

% For standard things: if not specified, pull from base workspace
if ~exist('templates','var')
    templates = evalin('base','templates');
end

if ~exist('channel_positions','var')
    channel_positions = evalin('base','channel_positions');
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

% Initialize figure and axes
cellraster_gui = figure;

waveform_axes = subplot(4,3,[1,4,7],'visible','off'); hold on;
waveform_lines = arrayfun(@(x) plot(waveform_axes,0,0,'k','linewidth',1),1:size(templates,3));

psth_axes = subplot(4,3,[2,3]);
psth_line = plot(0,0,'k','linewidth',2);
xlabel('Time from event (s)');
ylabel('Spikes');

raster_axes = subplot(4,3,[5,6,8,9],'YDir','reverse');
raster_dots = plot(0,0,'.k','linewidth',2);
xlabel('Time from event (s)');
ylabel('Trial');
axis tight;

amplitude_axes = subplot(4,1,4);
amplitude_dots = plot(0,0,'.k');
amplitude_lines = arrayfun(@(x) line([0,0],ylim),1:2);
xlabel('Experiment time (min)');
ylabel('Template amplitude');

% Set default raster times
raster_window = [-0.5,1];
psth_bin_size = 0.001;
t_bins = raster_window(1):psth_bin_size:raster_window(2);
t = t_bins(1:end-1) + diff(t_bins)./2;
use_align = reshape(align_times(~isnan(align_times)),[],1);
t_peri_event = use_align + t_bins;

% Set functions for key presses
set(cellraster_gui,'KeyPressFcn',@key_press);

% Package gui data
gui_data = struct;

% (plots)
gui_data.waveform_lines = waveform_lines;
gui_data.psth_line = psth_line;
gui_data.raster_dots = raster_dots;
gui_data.amplitude_dots = amplitude_dots;
gui_data.amplitude_lines = amplitude_lines; 

% (raster times)
gui_data.t = t;
gui_data.t_bins = t_bins;
gui_data.t_peri_event = t_peri_event;

% (spike data)
gui_data.templates = templates;
gui_data.channel_positions = channel_positions;
gui_data.curr_unit = 1;
gui_data.spike_times = spike_times;
gui_data.spike_templates = spike_templates;
gui_data.template_amplitudes = template_amplitudes;
gui_data.align_times = align_times;

% Upload gui data and draw
guidata(cellraster_gui, gui_data);
update_plot(cellraster_gui);

end


function update_plot(cellraster_gui,eventdata)

% Get guidata
gui_data = guidata(cellraster_gui);

% Index current spikes
curr_spikes_idx = gui_data.spike_templates == gui_data.curr_unit;

% Plot waveform across probe
yscale = 0.4;
xscale = 7;

y = permute(gui_data.templates(gui_data.curr_unit,:,:),[3,2,1]);
y = y - gui_data.channel_positions(:,2)*yscale;
x = (1:size(gui_data.templates,2)) + gui_data.channel_positions(:,1)*xscale;

template_channel_amp = squeeze(range(gui_data.templates(gui_data.curr_unit,:,:),2));
template_thresh = max(template_channel_amp)*0.2;
template_use_channels = template_channel_amp > template_thresh;
[~,max_channel] = max(max(abs(gui_data.templates(gui_data.curr_unit,:,:)),[],2),[],3);

arrayfun(@(ch) set(gui_data.waveform_lines(ch),'XData',x(ch,:),'YData',y(ch,:)),1:size(gui_data.templates,3));
arrayfun(@(ch) set(gui_data.waveform_lines(ch),'Color','r'),find(template_use_channels));
arrayfun(@(ch) set(gui_data.waveform_lines(ch),'Color','k'),find(~template_use_channels));
set(gui_data.waveform_lines(max_channel),'Color','b');

% Bin spikes, plot PSTH and raster


% (use only spikes within time range, big speed-up)
curr_raster_spike_times = gui_data.spike_times(curr_spikes_idx);
curr_raster_spike_times(curr_raster_spike_times < min(gui_data.t_peri_event(:)) | curr_raster_spike_times > max(gui_data.t_peri_event(:))) = [];

curr_raster = cell2mat(arrayfun(@(x) ...
    histcounts(curr_raster_spike_times,gui_data.t_peri_event(x,:)), ...
    [1:size(gui_data.t_peri_event,1)]','uni',false));

smooth_size = 50;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);
curr_psth = conv(sum(curr_raster,1),smWin,'same');

set(gui_data.psth_line,'XData',gui_data.t,'YData',curr_psth);
title(get(gui_data.psth_line,'Parent'),['Unit ' num2str(gui_data.curr_unit)]);

[raster_y,raster_x] = find(curr_raster);
set(gui_data.raster_dots,'XData',gui_data.t(raster_x),'YData',raster_y);
axis tight;

% Plot template amplitude over whole experiment
set(gui_data.amplitude_dots,'XData', ...
    gui_data.spike_times(curr_spikes_idx), ...
    'YData',gui_data.template_amplitudes(curr_spikes_idx));
set(gui_data.amplitude_lines(1),'XData',repmat(min(gui_data.t_peri_event(:)),2,1),'YData',ylim);
set(gui_data.amplitude_lines(2),'XData',repmat(max(gui_data.t_peri_event(:)),2,1),'YData',ylim);

axis tight

end


function key_press(cellraster_gui,eventdata)

% Get guidata
gui_data = guidata(cellraster_gui);

% Get unique units
unique_templates = unique(gui_data.spike_templates);

switch eventdata.Key
    case 'leftarrow'
        new_unit = unique_templates(find(unique_templates < gui_data.curr_unit,1,'last'));
        if isempty(new_unit)
            new_unit = unique_templates(end); 
        end
        gui_data.curr_unit = new_unit;       
    case 'rightarrow'
        new_unit = unique_templates(find(unique_templates > gui_data.curr_unit,1));
        if isempty(new_unit)
            new_unit = unique_templates(1); 
        end
        gui_data.curr_unit = new_unit;
end

% Upload gui data and draw
guidata(cellraster_gui, gui_data);
update_plot(cellraster_gui);
        
end






















