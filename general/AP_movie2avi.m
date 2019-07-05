function AP_movie2avi(im,framerate,color_map,color_axis,figure_position,savefile,annotation_text,ccf_overlay)
% AP_movie2avi(im,framerate,color_map,color_axis,figure_position,savefile,annotation_text,ccf_overlay)
% 
%
% Makes AVI movie from 3D matrix
% annotation_text - optional: string or cell array of text to put on each frame

if exist('annotation_text','var') && ~isempty(annotation_text)
   if ~iscell(annotation_text)
       annotation_text = repmat({annotation_text},size(im,3));
   end
   if length(annotation_text) ~= size(im,3)
       error('Annotation text different size than frame numbers')
   end
end

if ~isempty(figure_position)
    f = figure('Position',figure_position,'color','w');
else
    f = figure('color','w');
end
a = axes('Position',[0,0,1,1]);

im_plot = imagesc(a,im(:,:,1));
axis image off;

% Add dorsal cortex CCF outline
if exist('ccf_overlay','var') && ~isempty(ccf_overlay)
    AP_reference_outline('ccf_aligned','k',[],[size(im,1),size(im,2)/ccf_overlay,1,ccf_overlay]);
end

for i = 1:size(im,3)
    set(im_plot,'CData',im(:,:,i));
    if ~isempty(color_axis)
        caxis(color_axis);
    else
        caxis([-max(abs(caxis)),max(abs(caxis))]);
    end
    
    colormap(color_map);
   
   if exist('annotation_text','var') && ~isempty(annotation_text)
       annotation('textbox','Position',[0,0.93,0,0.07], ...
           'FitBoxToText','on','String',annotation_text{i},'Color','w','BackgroundColor','k','FontSize',20)
   end
   
   frames(i) = getframe(f);
end
close(f);
drawnow;

writerObj = VideoWriter(savefile);
writerObj.FrameRate = framerate;
open(writerObj);
writeVideo(writerObj,frames);
close(writerObj);





