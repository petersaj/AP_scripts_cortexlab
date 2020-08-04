function gui_fig = AP_image_scroll(images,t,colored)
% gui_fig = AP_image_scroll(images,t,colored)
% Make scrollable window from set of images
% images:
% 3D: grayscale (M,N,frames)
% 4D: can either be color (M,N,color,frames) 
%     or condition (M,N,frames,condition)
%     NOTE: in cases where the 3rd dimension is 3 this is ambiguous and
%     requires 'colored' to be true to display color
% t(optional): timepoints
%
% 'r' key: draw ROI and return temporal trace to workspace as 'roi'
% 4D condition images:
% up/down arrows switch between conditions
% 'r' returns one trace for each condition

dims = ndims(images);
if dims < 2
    error('Images matrix < 2D')
end

if ~exist('colored','var')
    colored = false;
end

if ~exist('t','var')
    t = [];
elseif ~iscell(t) && ~colored
    t = cellfun(@num2str,num2cell(t),'uni',false);
    if length(t) ~= size(images,3)
        error('t not same size as image dim 3')
    end
elseif colored && ~isempty(t) && ~iscell(t)
    t = cellfun(@num2str,num2cell(t),'uni',false);
    if length(t) ~= size(images,4)
        error('t not same size as image dim 4 (colored)')
    end
end

% Set up handles structure
handles = struct;
handles.data = images;
handles.t = t;
handles.curr_condition = 1;
handles.colored = colored;

% Store dimensions of matrix
handles.dims = dims;

% Create figure for scrolling and ROIs
gui_fig = figure;
set(gui_fig,'WindowScrollWheelFcn',{@imgSlider_MouseWheel, gui_fig});
set(gui_fig, 'KeyPressFcn', {@im_keypress, gui_fig});

% Set up scrollbar
ypos = [0 0 1 0.05];
handles.imgSlider = uicontrol('style','slider','units','normalized','position',ypos,'min',0,'max',1,'value',0);
set(handles.imgSlider,'Callback',{@imgSlider_Listener, gui_fig});

set(handles.imgSlider,'Enable','on');
if dims < 4
    n_images = size(images,3);
elseif dims == 4 && ~handles.colored
    n_images = size(images,3);
elseif dims == 4 && handles.colored
    n_images = size(images,4);
end
set(handles.imgSlider,'Min',1);
set(handles.imgSlider,'Max',n_images);
set(handles.imgSlider,'Value',1);
set(handles.imgSlider,'SliderStep',[1/n_images, 1/n_images]);

handles.n_images = n_images;

% Create axes, plot first image
handles.curr_im = 1;
if dims < 4
    handles.im = imagesc(handles.data(:,:,handles.curr_im));
elseif dims == 4 && ~handles.colored
    handles.im = imagesc(handles.data(:,:,handles.curr_im,handles.curr_condition));
elseif dims == 4 && handles.colored
    handles.im = image(handles.data(:,:,:,handles.curr_im));
end
axis off;
colormap(gray);

% Set the title with the frame number
if dims == 4 && ~handles.colored
    cond_text = [', Condition: ' num2str(handles.curr_condition)];
else
    cond_text = [];
end

if ~isempty(handles.t)
    handles.frame_num = title([handles.t{handles.curr_im} cond_text]);
else
    handles.frame_num = title(['Frame: ' num2str(handles.curr_im) cond_text]);
end

% Scale to max/min of entire image set
caxis(double([min(images(:)) max(images(:))]))

% Make figure toolbar available
set(gui_fig,'toolbar','figure');

% Update gui data
guidata(gui_fig,handles);



function imgSlider_Listener(currentObject, eventdata, gui_fig)

% Get guidata
handles = guidata(gui_fig);

% Get frame number from slider, round appropriately
new_im = get(handles.imgSlider,'Value');
new_im = round(new_im);
set(handles.imgSlider,'Value',new_im);

% Update the image
if handles.dims == 3
    set(handles.im,'Cdata',handles.data(:,:,new_im));
elseif handles.dims == 4 && ~handles.colored
    set(handles.im,'Cdata',handles.data(:,:,new_im,handles.curr_condition));
elseif handles.dims == 4 && handles.colored
    set(handles.im,'Cdata',handles.data(:,:,:,new_im));
end

% Update the title
if handles.dims == 4 && ~handles.colored
    cond_text = [', Condition: ' num2str(handles.curr_condition)];
else
    cond_text = [];
end

if ~isempty(handles.t)
    set(handles.frame_num,'String',[handles.t{new_im} cond_text]);
else
    set(handles.frame_num,'String',['Frame: ' num2str(new_im) cond_text]);
end

% Update guidata
handles.curr_im = new_im;
guidata(gui_fig, handles);



function imgSlider_MouseWheel(currentObject, eventdata, gui_fig)
% Executes when mouse wheel is scrolled in figure

% Get guidata
handles = guidata(gui_fig);

% Get current image
curr_im = handles.curr_im;

% Update current frame based on mouse wheel
mouse_wheel_count = eventdata.VerticalScrollCount;
new_im = curr_im + mouse_wheel_count;
if new_im > handles.n_images
    new_im = handles.n_images;
elseif new_im < 1
    new_im = 1;
end

% Update the image
if handles.dims == 3
    set(handles.im,'Cdata',handles.data(:,:,new_im));
elseif handles.dims == 4 && ~handles.colored
    set(handles.im,'Cdata',handles.data(:,:,new_im,handles.curr_condition));
elseif handles.dims == 4 && handles.colored
    set(handles.im,'Cdata',handles.data(:,:,:,new_im));
end

% Update slider
set(handles.imgSlider,'Value',new_im);

% Update the title
if handles.dims == 4 && ~handles.colored
    cond_text = [', Condition: ' num2str(handles.curr_condition)];
else
    cond_text = [];
end

if ~isempty(handles.t)
    set(handles.frame_num,'String',[handles.t{new_im} cond_text]);
else
    set(handles.frame_num,'String',['Frame: ' num2str(new_im) cond_text]);
end

% Update guidata
handles.curr_im = new_im;
guidata(gui_fig, handles);



function im_keypress(currentObject, eventdata, gui_fig)
% Executes when a key is pressed

% Get guidata
handles = guidata(gui_fig);

switch eventdata.Key
    case 'r'
        % draw ROI, output average traces and mask to workspace
        roiMask = roipoly;
        if handles.dims == 3
            roi.trace = nanmean(reshape(handles.data(repmat(roiMask,1,1,size(handles.data,3))),sum(roiMask(:)),[]),1);
        elseif handles.dims == 4 && ~handles.colored
            roi.trace = permute(nanmean(reshape(handles.data(repmat(roiMask,1,1,size(handles.data,3),size(handles.data,4))), ...
                sum(roiMask(:)),[],size(handles.data,4)),1),[3,2,1]);
        end
        roi.mask = roiMask;
        % push to base workspace
        assignin('base','roi',roi);
        
    case 'uparrow'
        if handles.dims == 4 && ~handles.colored
            new_condition = handles.curr_condition + 1;
            if new_condition > size(handles.data,4)
                return
            end
            
            set(handles.im,'Cdata',handles.data(:,:,handles.curr_im,new_condition));
            handles.curr_condition = new_condition;
            
            % Update the title
            cond_text = [', Condition: ' num2str(handles.curr_condition)];           
            if ~isempty(handles.t)
                set(handles.frame_num,'String',[handles.t{handles.curr_im} cond_text]);
            else
                set(handles.frame_num,'String',['Frame: ' num2str(handles.curr_im) cond_text]);
            end

        end
        
    case 'downarrow'
        if handles.dims == 4 && ~handles.colored
            new_condition = handles.curr_condition - 1;
            if new_condition < 1
                return
            end
            
            set(handles.im,'Cdata',handles.data(:,:,handles.curr_im,new_condition));
            handles.curr_condition = new_condition;
            
            % Update the title
            cond_text = [', Condition: ' num2str(handles.curr_condition)];           
            if ~isempty(handles.t)
                set(handles.frame_num,'String',[handles.t{handles.curr_im} cond_text]);
            else
                set(handles.frame_num,'String',['Frame: ' num2str(handles.curr_im) cond_text]);
            end
            
        end
        
end

% Update guidata
guidata(gui_fig, handles);




