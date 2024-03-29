function AP_mousemovie(movie_fn,movie_t,dlc)
% AP_mousemovie(movie_fn,movie_t,dlc)
%
% Scroll through a recorded movie of the mouse (face/eye)
% movie_fn = filename of movie
% movie_t = timestamps of movie (optional)
% dlc = DeepLabCut outputs (optional)

if ~exist(movie_fn,'file')
   error(['No file: ' movie_fn]); 
end

disp('Press s to save section of movies');

% Set up figure
gui_fig = figure; colormap(gray)
set(gui_fig,'WindowScrollWheelFcn',{@imgSlider_MouseWheel, gui_fig});
set(gui_fig, 'KeyPressFcn', {@im_keypress, gui_fig});

handles.cam_axis = axes('Ydir','reverse'); hold on; caxis([0,255]);
axis(handles.cam_axis,'image','off');

% Set up videoreader
handles.movie_vr = VideoReader(movie_fn);
movie_frame = 1;
movie_im = read(handles.movie_vr,movie_frame);
handles.movie_im = imagesc(handles.cam_axis,movie_im);
handles.n_frames = handles.movie_vr.NumberOfFrames;
handles.framerate = handles.movie_vr.FrameRate;

if ~exist('movie_t') || isempty(movie_t)
    handles.movie_t = (1:handles.n_frames)/handles.framerate;
elseif length(movie_t) ~= handles.n_frames
    error('Different number timestamps and frames');
else
    handles.movie_t = movie_t;
end

% Check DLC size if loaded
if exist('dlc','var')
    % Get marker names and x/y points in grid
    dlc_markers = fieldnames(dlc);
    dlc_points_x = cell2mat(struct2cell(structfun(@(s) s.x,dlc,'UniformOutput',false))');
    dlc_points_y = cell2mat(struct2cell(structfun(@(s) s.y,dlc,'UniformOutput',false))');

    if size(dlc_points_x,1) ~= handles.n_frames
        error('Different number frames and DLC points')
    end

    % Draw/store in handles
    handles.dlc.markers = dlc_markers;
    handles.dlc.x = dlc_points_x;
    handles.dlc.y = dlc_points_y;

    handles.dlc.dots = ...
        scatter(handles.dlc.x(movie_frame,:),handles.dlc.y(movie_frame,:), ...
        50,jet(length(dlc_markers)),'filled');
end


% Set up scrollbar (use timer function to prevent lag)
ypos = [0 0 1 0.05];
handles.imgSlider = uicontrol('style','slider','units','normalized','position',ypos,'min',0,'max',1,'value',0);
set(handles.imgSlider,'Callback',{@imgSlider_Listener_timer,gui_fig});
handles.scroll_timer = timer('ExecutionMode','singleShot','TimerFcn',{@imgSlider_Listener,gui_fig});

set(handles.imgSlider,'Min',1);
set(handles.imgSlider,'Max',handles.n_frames);
set(handles.imgSlider,'Value',1);
set(handles.imgSlider,'SliderStep',[1/handles.n_frames, 100/handles.n_frames]);

% Set up frame title
handles.frame_text = uicontrol('Style','text','String', ...
    sprintf('Time: %0.2f s',1/handles.framerate),'FontSize',14,'Units', ...
    'Normalized','Position',[0.3,0.93,0.4,0.07]);

% Update guidata
handles.movie_frame = movie_frame;
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
movie_frame = get(handles.imgSlider,'Value');
movie_frame = round(movie_frame);
set(handles.imgSlider,'Value',movie_frame);

% Update the images
update_im(handles,gui_fig,movie_frame);


function imgSlider_MouseWheel(currentObject, eventdata, gui_fig)
% Executes when mouse wheel is scrolled in figure

% Get guidata
handles = guidata(gui_fig);

% Get current frame
movie_frame = handles.movie_frame;

% Update current frame based on mouse wheel
mouse_wheel_count = eventdata.VerticalScrollCount;
movie_frame = movie_frame + mouse_wheel_count;
if movie_frame > handles.n_frames
    movie_frame = handles.n_frames;
elseif movie_frame < 1
    movie_frame = 1;
end

% Set the slider
set(handles.imgSlider,'Value',movie_frame);

% Update the images
update_im(handles,gui_fig,movie_frame);


function im_keypress(currentObject, eventdata, gui_fig)
% Executes when a key is pressed
% HAVEN'T UPDATED THIS FROM AP_wfmovies YET

% Get guidata
handles = guidata(gui_fig);

switch eventdata.Key
    
    % Save section of images as movie
    case 's'
        
        % Get options
        disp('Preparing to make movie:')
        record_t = input('Start/stop time (e.g. [0 5]): ');
        record_framerate = input('Framerate: ');
        [save_file,save_path] = uiputfile('.avi','Choose save location');
        save_filename = [save_path save_file];
        
        frame_times = (1:handles.n_frames)/handles.framerate;
        record_frames = find(frame_times >= record_t(1) & ...
            frame_times <= record_t(2));
        n_record_frames = length(record_frames);       
        
        % Run through selected frames and save
        disp('Recording...')
        movie_frames(n_record_frames) = struct('cdata',[],'colormap',[]);
        for curr_movie_frame_idx = 1:n_record_frames
            curr_movie_frame = record_frames(curr_movie_frame_idx);
            % Update images
            update_im(handles,gui_fig,curr_movie_frame);
            movie_frames(curr_movie_frame_idx) = getframe(gui_fig);
        end
        
        % Write movie
        disp('Saving...')
        writerObj = VideoWriter(save_filename);
        writerObj.FrameRate = record_framerate;
        open(writerObj);
        writeVideo(writerObj,movie_frames);
        close(writerObj);
        
        disp('Done.')
        
end

% Update guidata
guidata(gui_fig, handles);


function update_im(handles, gui_fig, movie_frame)

% Update the images
movie_im = read(handles.movie_vr,movie_frame);
set(handles.movie_im,'Cdata',movie_im);

% Update the frame text
set(handles.frame_text,'String',sprintf('Time: %0.2f s',handles.movie_t(movie_frame)));

% Update DLC (if loaded)
if isfield(handles,'dlc')
    set(handles.dlc.dots,'XData',handles.dlc.x(movie_frame,:), ...
        'YData',handles.dlc.y(movie_frame,:));
end

% Update guidata
handles.movie_frame = movie_frame;
guidata(gui_fig, handles);

drawnow;



