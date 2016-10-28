function AP_wfmovies(U,V,frame_t,eyecam_fn,eyecam_t,facecam_fn,facecam_t,trace,trace_t)
% AP_wfmovies(U,V,frame_t,eyecam_fn,eyecam_t,facecam_fn,facecam_t)
%
% U = wf U
% V = wf V
% frame_t = wf frames in timeline time
%
% eyecam_fn = filename of eyecam movie
% eyecam_t = eyecam frames in timeline time
%
% facecam_fn = filename of facecam movie
% facecam_t = facecam frames in timeline time
%
% trace = trace to plot
% trace_t = trace times

if exist('eyecam_fn','var') && exist('eyecam_t','var')
   handles.plot_eyecam = true; 
else
    handles.plot_eyecam = false;
end

if exist('facecam_fn','var') && exist('facecam_t','var')
   handles.plot_facecam = true; 
else
    handles.plot_facecam = false;
end

if exist('trace','var') && exist('trace_t','var')
   handles.plot_trace = true; 
else
    handles.plot_trace = false;
end

% Set up subplots (4: widefield, two behavior cams, one for traces)
gui_fig = figure; colormap(gray)
set(gui_fig,'WindowScrollWheelFcn',{@imgSlider_MouseWheel, gui_fig});
set(gui_fig, 'KeyPressFcn', {@im_keypress, gui_fig});

handles.eyecam_axis = subplot(2,2,1); 
handles.facecam_axis = subplot(2,2,2);
handles.wf_axis = subplot(2,2,3);
handles.trace_axis = subplot(2,2,4);

% Set up widefield data and image
handles.V = V; 
handles.U = U;
n_frames = size(V,2);
framerate = median(diff(frame_t));

handles.frame_t = frame_t;

wf_frame = 1;
handles.t = frame_t(wf_frame);

wf_im = svdFrameReconstruct(U,V(:,wf_frame));
handles.wf_im = imagesc(handles.wf_axis,wf_im); axis(handles.wf_axis,'off');

% Set up videoreaders, relative times, and images of cameras
discretize_times = [frame_t,frame_t(end)+framerate];
wf_frame_idx = frame_t;

if handles.plot_eyecam
    handles.eyecam_vr = VideoReader(eyecam_fn);
    handles.eyecam_frame_idx = discretize(eyecam_t,discretize_times);
    
    eyecam_frame = find(handles.eyecam_frame_idx == wf_frame,1);
    eyecam_im = read(handles.eyecam_vr,eyecam_frame);
    handles.eyecam_im = imagesc(handles.eyecam_axis,eyecam_im); axis(handles.eyecam_axis,'off');   
end

if handles.plot_facecam
    handles.facecam_vr = VideoReader(facecam_fn);
    handles.facecam_frame_idx = discretize(facecam_t,discretize_times);
        
    facecam_frame = find(handles.facecam_frame_idx == wf_frame,1);
    facecam_im = read(handles.facecam_vr,facecam_frame);
    handles.facecam_im = imagesc(handles.facecam_axis,facecam_im); axis(handles.facecam_axis,'off');
end

% Set up trace
if handles.plot_trace
    handles.trace_plot = plot(handles.trace_axis,trace_t,trace);
    handles.trace_xlim = [handles.t-5,handles.t+5];
    handles.trace_tmark = line([handles.t,handles.t],ylim,'color','k');
    xlim(handles.trace_axis,handles.trace_xlim);
end

% Set up scrollbar (use timer function to prevent lag)
ypos = [0 0 1 0.05];
handles.imgSlider = uicontrol('style','slider','units','normalized','position',ypos,'min',0,'max',1,'value',0);
set(handles.imgSlider,'Callback',{@imgSlider_Listener_timer,gui_fig});
handles.scroll_timer = timer('ExecutionMode','singleShot','TimerFcn',{@imgSlider_Listener,gui_fig});

set(handles.imgSlider,'Min',1);
set(handles.imgSlider,'Max',n_frames);
set(handles.imgSlider,'Value',1);
set(handles.imgSlider,'SliderStep',[10/n_frames, 100/n_frames]);

% Set up time title
handles.time_text = uicontrol('Style','text','String', ...
    ['Time: ' num2str(handles.t) ,'s'],'FontSize',14,'Units', ...
    'Normalized','Position',[0.3,0.93,0.4,0.07]);

% Update guidata
handles.wf_frame = wf_frame;
guidata(gui_fig, handles);

drawnow;


function imgSlider_Listener_timer(currentObject, eventdata, gui_fig)
% Get guidata
handles = guidata(gui_fig);

% Only run if not already running
if strcmp(handles.scroll_timer.Running,'off')
    start(handles.scroll_timer);
end


function imgSlider_Listener(currentObject, eventdata, gui_fig)
% Executes whenever the slider is pressed

% Get guidata
handles = guidata(gui_fig);

% Get frame number from slider, round appropriately
wf_frame = get(handles.imgSlider,'Value');
wf_frame = round(wf_frame);
set(handles.imgSlider,'Value',wf_frame);

% Update the images
update_im(handles,gui_fig,wf_frame);


function imgSlider_MouseWheel(currentObject, eventdata, gui_fig)
% Executes when mouse wheel is scrolled in figure

% Get guidata
handles = guidata(gui_fig);

% Get current frame
wf_frame = handles.wf_frame;

% Update current frame based on mouse wheel
mouse_wheel_count = eventdata.VerticalScrollCount;
wf_frame = wf_frame + mouse_wheel_count;
if wf_frame > length(handles.frame_t)
    wf_frame = length(handles.frame_t);
elseif wf_frame < 1
    wf_frame = 1;
end

% Set the slider
set(handles.imgSlider,'Value',wf_frame);

% Update the images
update_im(handles,gui_fig,wf_frame);


function im_keypress(currentObject, eventdata, gui_fig)
% Executes when a key is pressed

% Get guidata
handles = guidata(gui_fig);

switch eventdata.Key
    
    % Save section of images as movie
    case 's'
        
        % Get options
        disp('Preparing to make movie:')
        movie_t = input('Start/stop time (e.g. [0 5]): ');
        movie_framerate = input('Framerate: ');
        save_filename = uiputfile('.avi','Choose save location');
        
        movie_wf_frames = handles.frame_t > movie_t(1) & ...
            handles.frame_t < movie_t(2);
        n_movie_frames = sum(movie_wf_frames);
        
        % Run through selected frames and save
        disp('Recording...')
        movie_frames(n_movie_frames) = struct('cdata',[],'colormap',[]);
        for curr_movie_frame = find(movie_wf_frames);
            % Update images
            update_im(handles,gui_fig,curr_movie_frame);
            movie_frames(curr_movie_frame) = getframe(gui_fig);
        end
        
        % Write movie
        disp('Saving...')
        writerObj = VideoWriter(save_filename);
        writerObj.FrameRate = movie_framerate;
        open(writerObj);
        writeVideo(writerObj,movie_frames);
        close(writerObj);
        
        disp('Done.')
        
end

% Update guidata
guidata(gui_fig, handles);


function update_im(handles, gui_fig, wf_frame)

handles.t = handles.frame_t(wf_frame);

% Update the images
if handles.plot_eyecam
    eyecam_frame = find(handles.eyecam_frame_idx == wf_frame,1);
    eyecam_im = read(handles.eyecam_vr,eyecam_frame);
    set(handles.eyecam_im,'Cdata',eyecam_im);
end

if handles.plot_facecam
    facecam_frame = find(handles.facecam_frame_idx == wf_frame,1);
    facecam_im = read(handles.facecam_vr,facecam_frame);
    set(handles.facecam_im,'Cdata',facecam_im);
end

wf_im = svdFrameReconstruct(handles.U,handles.V(:,wf_frame));

set(handles.wf_im,'Cdata',wf_im);

% Update trace
if handles.plot_trace
    handles.trace_xlim = [handles.t-5,handles.t+5];
    xlim(handles.trace_axis,handles.trace_xlim);
    set(handles.trace_tmark,'XData',[handles.t,handles.t]);
end

% Update the time text
set(handles.time_text,'String',['Time: ' num2str(handles.t) ,'s']);

% Update guidata
handles.wf_frame = wf_frame;
guidata(gui_fig, handles);

drawnow;


