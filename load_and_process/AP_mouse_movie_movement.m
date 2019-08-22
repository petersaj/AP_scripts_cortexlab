function AP_mouse_movie_movement(animal,day,experiment)
% AP_mouse_movie_movement(animal,day,experiment)
%
% Simple process for getting movement from movie: threshold, sum(abs(diff))
% across frames

[facecam_fn,facecam_exists] = AP_cortexlab_filename(animal,day,experiment,'facecam');

if facecam_exists
    disp(['Getting movement for ' animal ' ' day ' ' num2str(experiment)]);
else
    error(['No facecam for ' animal ' ' day ' ' num2str(experiment)])
end

% Load in set of frames from middle of experiment
vr = VideoReader(facecam_fn);
n_frames = vr.NumberOfFrames;
n_frames_surround = 1000;

time_start = 2; % minutes in to start
frame_start = round(vr.FrameRate*time_start*60);
read_frames = [frame_start - round(n_frames_surround/2),frame_start + round(n_frames_surround/2)];
sample_frames = read(vr,read_frames);

draw_sample_frame = sample_frames(:,:,1);

f = figure('units','normalized','outerposition',[0 0 1 1]);
h = axes; colormap(gray); caxis([0,255]);
imagesc(draw_sample_frame);
axis image off
title(h,'Set vertical cutoff');
[~,vertical_cutoff] = ginput(1);
vertical_cutoff = round(vertical_cutoff); 

title(h,'ROI over hand');
hand_roi = roipoly;

close(f)

hand_px = draw_sample_frame(hand_roi);
hand_thresh = floor(nanmedian(hand_px)/2);

sample_frames_handmove = padarray(diff( ...
    sample_frames(vertical_cutoff:end,:,:) > hand_thresh,[],3) > 0,[0,0,1],0,'pre');

sample_frames_move = ...
    nanmean(reshape(sample_frames_handmove,[],size(sample_frames_handmove,3)),1);

% Summary plot
f = figure('units','normalized','outerposition',[0 0 1 1]);

% Summary plot: mean, std, motion
subplot(4,2,1); colormap(gray);
imagesc(draw_sample_frame(vertical_cutoff:end,:) > hand_thresh);
axis image off
title('Sample frame thresh');

subplot(4,2,3); colormap(gray);
imagesc(nansum(sample_frames_handmove,3));
axis image off
title('Mean');

subplot(4,2,5); colormap(gray);
imagesc(nanstd(sample_frames_handmove,[],3));
axis image off
title('Std');

subplot(4,2,7);
plot(sample_frames_move,'k');
ylabel('Movement');
xlabel('Frame');

% Summary plot: get and plot first 10 components of svd
[u,s,v] = svd(+reshape(sample_frames_handmove,[],size(sample_frames_handmove,3)),'econ');
u = reshape(u,size(sample_frames_handmove,1),size(sample_frames_handmove,2),[]);
p4 = subplot(1,2,2);
imagesc(reshape(permute(u(:,:,1:10),[1,3,2]),[],size(u,2)));
axis image off
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(p4,brewermap([],'*RdBu'));
caxis(caxis/2)
title('SVD');

% Check with user that it looks good, then load and process all frames
user_check = input('Get movement on all frames? (y/n) ','s');
close(f);

if strcmp(lower(user_check),'y')
    
    chunk_frames_n = 5000;
    chunk_frames = round(linspace(1,n_frames,round(n_frames/chunk_frames_n)));
    
    frame_movement = nan(n_frames,1);
    for curr_chunk = 1:length(chunk_frames)-1        
        if curr_chunk == 1
            curr_chunk_frames = chunk_frames(curr_chunk):chunk_frames(curr_chunk+1);
        else
            % (fill chunk gaps by using one overlapping frame)
            curr_chunk_frames = chunk_frames(curr_chunk)-1:chunk_frames(curr_chunk+1);
        end
               
        curr_movie = read(vr,[curr_chunk_frames(1),curr_chunk_frames(end)]);

        curr_movie_movement = diff(curr_movie(vertical_cutoff:end,:,:) > hand_thresh,[],3) > 0;
        frame_movement(curr_chunk_frames(2:end)) = ...
            nanmean(reshape(curr_movie_movement,[],size(curr_movie_movement,3)),1);
        
        AP_print_progress_fraction(curr_chunk,length(chunk_frames)-1);       
    end
    
    facecam_dir = fileparts(facecam_fn);
    facecam_movement_fn = [facecam_dir filesep 'facecam_movement.mat'];
    save(facecam_movement_fn,'frame_movement');
    
else
    disp('User cancelled')
end













