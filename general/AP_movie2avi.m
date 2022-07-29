function AP_movie2avi(im,framerate,color_map,color_axis,figure_position,savefile,t_annotation,movie_annotation)
% AP_movie2avi(im,framerate,color_map,color_axis,figure_position,savefile,t_annotation,movie_annotation)
% 
% Makes AVI movie from 3D matrix
% im - 3D (time in dim 3) or 4D (multiple movies in dim 4)
% t_annotation - annotation for text (upper left corner)
% movie_annotation - title of each movie (if 4D)

% Prepare annotation text
if exist('t_annotation','var') && ~isempty(t_annotation)
    if length(t_annotation) ~= size(im,3)
        error('Annotation text different size than frame numbers')
    end
    if ~iscell(t_annotation)
        t_annotation = num2cell(t_annotation);
    end
end

% Create figure
if ~isempty(figure_position)
    f = figure('Position',figure_position,'color','w');
else
    f = figure('color','w');
end

% Make one axis per movie (if > 4 movies, split into 2 rows)
n_movies = size(im,4);
ax = gobjects(1,n_movies);
im_plot = gobjects(1,n_movies);
for curr_movie = 1:n_movies
    if n_movies <= 4
        subplot_y = 1;
        subplot_x = n_movies;
    elseif n_movies > 4
        subplot_y = 2;
        subplot_x = ceil(n_movies/subplot_y);
    end
    
    ax(curr_movie) = subplot(subplot_y,subplot_x,curr_movie);
    im_plot(curr_movie) = imagesc(ax(curr_movie),im(:,:,1,curr_movie));
    axis image off;
    
    if ~isempty(color_axis)
        caxis(ax(curr_movie),color_axis);
    else
        caxis(ax(curr_movie),[-max(abs(caxis)),max(abs(caxis))]);
    end
    colormap(color_map);
    
    if exist('movie_annotation','var') && ~isempty(movie_annotation)
        title(movie_annotation{curr_movie},'FontSize',16);
    end

end

% Add dorsal cortex CCF outline
for curr_movie = 1:n_movies
    axes(ax(curr_movie));
    AP_reference_outline('ccf_aligned_lefthemi',[0.5,0.5,0.5]);
end

for curr_t = 1:size(im,3)
    for curr_movie = 1:n_movies
        set(im_plot(curr_movie),'CData',im(:,:,curr_t,curr_movie));             
        if curr_movie == 1 && exist('t_annotation','var') && ~isempty(t_annotation)
            annotation('textbox','Position',[0,0.93,0,0.07], ...
                'FitBoxToText','on','String',t_annotation{curr_t},'Color','w','BackgroundColor','k','FontSize',16)
        end
    end
    drawnow;
    frames(curr_t) = getframe(f);
end

close(f);
drawnow;

writerObj = VideoWriter(savefile);
writerObj.FrameRate = framerate;
open(writerObj);
writeVideo(writerObj,frames);
close(writerObj);





