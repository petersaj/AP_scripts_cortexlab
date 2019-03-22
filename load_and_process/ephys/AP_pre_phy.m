function AP_pre_phy(animal,day)
% First-pass blunt good/bad
% WORKING ON THIS AT THE MOMENT - was unused previously

[ephys_path,ephys_exists] = AP_cortexlab_filename(animal,day,experiment,'ephys',site);

% Load spike data
if isfield(header,'sample_rate')
    ephys_sample_rate = str2num(header.sample_rate);
elseif isfield(header,'ap_sample_rate')
    ephys_sample_rate = str2num(header.ap_sample_rate);
end
spike_times = double(readNPY([ephys_path filesep 'spike_times.npy']))./ephys_sample_rate;
spike_templates = readNPY([ephys_path filesep 'spike_templates.npy']);
templates_whitened = readNPY([ephys_path filesep 'templates.npy']);
channel_positions = readNPY([ephys_path filesep 'channel_positions.npy']);
channel_map = readNPY([ephys_path filesep 'channel_map.npy']);
winv = readNPY([ephys_path filesep 'whitening_mat_inv.npy']);
template_amplitudes = readNPY([ephys_path filesep 'amplitudes.npy']);

% Default channel map/positions are from end: make from surface
channel_positions(:,2) = max(channel_positions(:,2)) - channel_positions(:,2);

% Unwhiten templates
templates = zeros(size(templates_whitened));
for t = 1:size(templates_whitened,1)
    templates(t,:,:) = squeeze(templates_whitened(t,:,:))*winv;
end

% Get the waveform of all templates (channel with largest amplitude)
[~,max_site] = max(max(abs(templates),[],2),[],3);
templates_max = nan(size(templates,1),size(templates,2));
for curr_template = 1:size(templates,1)
    templates_max(curr_template,:) = ...
        templates(curr_template,:,max_site(curr_template));
end
waveforms = templates_max;

% Get depth of each template (by center-of-mass)
template_chan_amp = squeeze(range(templates,2));
templateDepths = sum(template_chan_amp.*channel_positions(:,2)',2)./sum(template_chan_amp,2);


%%%%%%%%%%%%%%%%%%%%%%% (old, not dealt with)


local_phy_path = 'C:\data_temp\phy';

% Get the maximum amplitude waveform for each templates
templates = readNPY([ephys_path filesep 'templates.npy']);
channel_positions = readNPY([ephys_path filesep 'channel_positions.npy']);
channel_map = readNPY([ephys_path filesep 'channel_map.npy']);
winv = readNPY([ephys_path filesep 'whitening_mat_inv.npy']);

% Default channel map/positions are from end: make from surface
channel_positions(:,2) = max(channel_positions(:,2)) - channel_positions(:,2);

% Get the depths of each template (REMOVE DEPENDENCY EVENTUALLY)
% (by COM - this used to not work but now looks ok)
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = ...
    templatePositionsAmplitudes(templates,winv,channel_positions(:,2),spike_templates,template_amplitudes);

% Set the templates to be the unwhitened templates
templates = tempsUnW;


templates = readNPY([local_phy_path filesep 'templates.npy']);
template_abs = permute(max(abs(templates),[],2),[3,1,2]);
[~,max_channel_idx] =  max(template_abs,[],1);

templates_maxchan = cell2mat(arrayfun(@(x) permute(templates(x,:,max_channel_idx(x)),[2,1,3]),1:size(templates,1),'uni',false));

% Load clusters, if they exist
waveform_classification = zeros(size(templates,1),1);

cluster_filepattern = [ephys_path 'cluster_group*'];
cluster_filedir = dir(cluster_filepattern);
if ~isempty(cluster_filedir)
    cluster_filename = [ephys_path cluster_filedir.name];
    fid = fopen(cluster_filename);
    cluster_groups = textscan(fid,'%d%s','HeaderLines',1);
    fclose(fid);
end

% Set up the figure and GUI data
gui_fig = figure;
set(gui_fig,'KeyPressFcn',@key_press);

waveform_plot = plot(templates_maxchan(:,1),'linewidth',3,'color','k');

gui_data = struct;
gui_data.waveform_plot = waveform_plot;
gui_data.current_waveform = 1;
gui_data.waveforms = templates_maxchan;
gui_data.waveform_classification = waveform_classification;

% Set figure color based on classification
if gui_data.waveform_classification(gui_data.current_waveform) == 1
    set(gui_fig,'color','g');
elseif gui_data.waveform_classification(gui_data.current_waveform) == -1
    set(gui_fig,'color','r');
else
    set(gui_fig,'color','w');
end

% Update GUI data
guidata(gui_fig, gui_data);

end

function key_press(gui_fig,eventdata)

% Get GUI data
gui_data = guidata(gui_fig);

switch eventdata.Key
    
    case 'uparrow'
        % Classify as good, move to next
        gui_data.waveform_classification(gui_data.current_waveform) = 1;
        
        if gui_data.current_waveform < size(gui_data.waveforms,2)
            gui_data.current_waveform = gui_data.current_waveform + 1;
        end
        set(gui_data.waveform_plot,'YData',gui_data.waveforms(:,gui_data.current_waveform));
        
        if gui_data.waveform_classification(gui_data.current_waveform) == 1
            set(gui_fig,'color','g');
        elseif gui_data.waveform_classification(gui_data.current_waveform) == -1
            set(gui_fig,'color','r');
        else
            set(gui_fig,'color','w');
        end
        disp(gui_data.current_waveform);
                
    case 'downarrow'
        % Classify as bad, move to next
        gui_data.waveform_classification(gui_data.current_waveform) = -1;
        
        if gui_data.current_waveform < size(gui_data.waveforms,2)
            gui_data.current_waveform = gui_data.current_waveform + 1;
        end
        set(gui_data.waveform_plot,'YData',gui_data.waveforms(:,gui_data.current_waveform));
        
        if gui_data.waveform_classification(gui_data.current_waveform) == 1
            set(gui_fig,'color','g');
        elseif gui_data.waveform_classification(gui_data.current_waveform) == -1
            set(gui_fig,'color','r');
        else
            set(gui_fig,'color','w');
        end
        disp(gui_data.current_waveform);
        
    case 'rightarrow'
        % Move to next waveform
        if gui_data.current_waveform < size(gui_data.waveforms,2)
            gui_data.current_waveform = gui_data.current_waveform + 1;
        end
        set(gui_data.waveform_plot,'YData',gui_data.waveforms(:,gui_data.current_waveform));
        
        if gui_data.waveform_classification(gui_data.current_waveform) == 1
            set(gui_fig,'color','g');
        elseif gui_data.waveform_classification(gui_data.current_waveform) == -1
            set(gui_fig,'color','r');
        else
            set(gui_fig,'color','w');
        end
        disp(gui_data.current_waveform);
        
    case 'leftarrow'
        % Move to last waveform
        if gui_data.current_waveform > 1
            gui_data.current_waveform = gui_data.current_waveform - 1;
        end
        set(gui_data.waveform_plot,'YData',gui_data.waveforms(:,gui_data.current_waveform));
        
        if gui_data.waveform_classification(gui_data.current_waveform) == 1
            set(gui_fig,'color','g');
        elseif gui_data.waveform_classification(gui_data.current_waveform) == -1
            set(gui_fig,'color','r');
        else
            set(gui_fig,'color','w');
        end
        disp(gui_data.current_waveform);
        
    case 's'
        % Save classifications in CSV
        cluster_groups_text = cell(size(gui_data.waveforms,2),2);        
        cluster_groups_text(1:end,1) = num2cell(0:size(gui_data.waveforms,2)-1);       
        cluster_groups_text(gui_data.waveform_classification == 0,2) = {'unsorted'};
        cluster_groups_text(gui_data.waveform_classification == 1,2) = {'good'};
        cluster_groups_text(gui_data.waveform_classification == -1,2) = {'noise'};
        
%         % Write
%         fid = fopen(gui_data.local_cluster_groups_filename,'w');
%         fprintf(fid,'%s\t%s\n','cluster_id','group');
%         for curr_template = 1:size(gui_data.waveforms,2)
%             fprintf(fid,'%d\t%s\n',cluster_groups_text{curr_template,:});
%         end
%         fclose(fid);
%         
%         disp('Saved classifications')
        
end

guidata(gui_fig, gui_data);

end




