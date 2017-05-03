function AP_pre_phy
% First-pass blunt good/bad

local_phy_path = 'C:\data_temp\phy';
local_cluster_groups_filename = [local_phy_path filesep 'cluster_groups.csv'];

% Get the maximum amplitude waveform for each templates
templates = readNPY([local_phy_path filesep 'templates.npy']);
template_abs = permute(max(abs(templates),[],2),[3,1,2]);
[~,max_channel_idx] =  max(template_abs,[],1);

templates_maxchan = cell2mat(arrayfun(@(x) permute(templates(x,:,max_channel_idx(x)),[2,1,3]),1:size(templates,1),'uni',false));

% Load clusters, if they exist
waveform_classification = zeros(size(templates,1),1);
if exist(local_cluster_groups_filename,'file')
    fid = fopen(local_cluster_groups_filename);
    cluster_groups = textscan(fid,'%d%s','HeaderLines',1);
    
    waveform_classification(strcmp(cluster_groups{2},'good')) = 1;
    waveform_classification(strcmp(cluster_groups{2},'noise')) = -1;

    fclose(fid);   
end

% Set up the figure and GUI data
gui_fig = figure;
set(gui_fig,'KeyPressFcn',@key_press);

waveform_plot = plot(templates_maxchan(:,1),'linewidth',3,'color','k');

gui_data = struct;
gui_data.local_cluster_groups_filename = local_cluster_groups_filename;
gui_data.waveform_plot = waveform_plot;
gui_data.current_waveform = 1;
gui_data.waveforms = templates_maxchan;
gui_data.waveform_classification = waveform_classification;

% Set figure color based on classification
if gui_data.waveform_classification(gui_data.current_waveform) == 1;
    set(gui_fig,'color','g');
elseif gui_data.waveform_classification(gui_data.current_waveform) == -1;
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
        
        if gui_data.current_waveform < size(gui_data.waveforms,2);
            gui_data.current_waveform = gui_data.current_waveform + 1;
        end
        set(gui_data.waveform_plot,'YData',gui_data.waveforms(:,gui_data.current_waveform));
        
        if gui_data.waveform_classification(gui_data.current_waveform) == 1;
            set(gui_fig,'color','g');
        elseif gui_data.waveform_classification(gui_data.current_waveform) == -1;
            set(gui_fig,'color','r');
        else
            set(gui_fig,'color','w');
        end
        disp(gui_data.current_waveform);
                
    case 'downarrow'
        % Classify as bad, move to next
        gui_data.waveform_classification(gui_data.current_waveform) = -1;
        
        if gui_data.current_waveform < size(gui_data.waveforms,2);
            gui_data.current_waveform = gui_data.current_waveform + 1;
        end
        set(gui_data.waveform_plot,'YData',gui_data.waveforms(:,gui_data.current_waveform));
        
        if gui_data.waveform_classification(gui_data.current_waveform) == 1;
            set(gui_fig,'color','g');
        elseif gui_data.waveform_classification(gui_data.current_waveform) == -1;
            set(gui_fig,'color','r');
        else
            set(gui_fig,'color','w');
        end
        disp(gui_data.current_waveform);
        
    case 'rightarrow'
        % Move to next waveform
        if gui_data.current_waveform < size(gui_data.waveforms,2);
            gui_data.current_waveform = gui_data.current_waveform + 1;
        end
        set(gui_data.waveform_plot,'YData',gui_data.waveforms(:,gui_data.current_waveform));
        
        if gui_data.waveform_classification(gui_data.current_waveform) == 1;
            set(gui_fig,'color','g');
        elseif gui_data.waveform_classification(gui_data.current_waveform) == -1;
            set(gui_fig,'color','r');
        else
            set(gui_fig,'color','w');
        end
        disp(gui_data.current_waveform);
        
    case 'leftarrow'
        % Move to last waveform
        if gui_data.current_waveform > 1;
            gui_data.current_waveform = gui_data.current_waveform - 1;
        end
        set(gui_data.waveform_plot,'YData',gui_data.waveforms(:,gui_data.current_waveform));
        
        if gui_data.waveform_classification(gui_data.current_waveform) == 1;
            set(gui_fig,'color','g');
        elseif gui_data.waveform_classification(gui_data.current_waveform) == -1;
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
        
        % Write
        fid = fopen(gui_data.local_cluster_groups_filename,'w');
        fprintf(fid,'%s\t%s\n','cluster_id','group');
        for curr_template = 1:size(gui_data.waveforms,2)
            fprintf(fid,'%d\t%s\n',cluster_groups_text{curr_template,:});
        end
        fclose(fid);
        
        disp('Saved classifications')
        
end

guidata(gui_fig, gui_data);

end




