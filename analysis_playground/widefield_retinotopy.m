%% Retinotopy via sweep fourier (set to use for hemodynamic)

% Downsample U
U_downsample_factor = 4; %2 if max method
Ud = imresize(Uh,1/U_downsample_factor,'bilinear');

framerate = 1./nanmedian(diff(frame_t));

% Get average response to stim

% Get parameters of stim
stim_duration = unique(Protocol.pars(strcmp('dur',Protocol.parnames),:)/10);
if length(stim_duration) > 1
    error('Different stim durations')
end
stim_freq = Protocol.pars(strcmp('tf',Protocol.parnames),:)/100;
stim_direction = Protocol.pars(strcmp('dir',Protocol.parnames),:);
stim_orientation = Protocol.pars(strcmp('ori',Protocol.parnames),:);
stim_range = [Protocol.pars(strcmp('start',Protocol.parnames),:); ...
    Protocol.pars(strcmp('end',Protocol.parnames),:)];

surround_window = [0,stim_duration];
surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);

% Loop through conditions, get average response
im_stim = nan(size(Ud,1),size(Ud,2),length(surround_time),max(unique(stimIDs)));
for curr_condition = unique(stimIDs)'
    
    use_stims = find(stimIDs == curr_condition);
    use_stim_onsets = stimOn_times(use_stims);
    
    stim_surround_times = bsxfun(@plus, use_stim_onsets(:), surround_time);
    peri_stim_v = permute(nanmean(interp1(frame_t,fVh',stim_surround_times),1),[3,2,1]);
    
    im_stim(:,:,:,curr_condition) = svdFrameReconstruct(Ud,peri_stim_v);
end

% Get fourier component at stimulus frequency
ComplexMaps = zeros(size(im_stim,1),size(im_stim,2),size(im_stim,4));
yy = permute(2*exp(-surround_time*2*pi*1i*stim_freq(1)),[1,3,2]);
aaa = repmat(yy,[size(im_stim,1),size(im_stim,2)]);

for curr_condition = unique(stimIDs)'
    ComplexMaps(:,:,curr_condition) = double(mean(im_stim(:,:,:,curr_condition).*aaa, 3));
end

AbsMaps = abs(ComplexMaps);
AngleMaps = angle(ComplexMaps);

% Combine maps of same orientation and opposite directions (just hard coded now)
angle_maps = nan(size(im_stim,1),size(im_stim,2),2);
retinotopy_maps = nan(size(im_stim,1),size(im_stim,2),2);
for curr_orientation = 1:2
    
    curr_stims = find(stim_orientation == curr_orientation');
    
    AbsolutePhaseS = sum(bsxfun(@times,AngleMaps(:,:,curr_stims),permute(stim_direction(curr_stims),[1,3,2])),3);
    DoubleDelayMap = sum(AngleMaps(:,:,curr_stims),3);
    DoubleDelayMap(DoubleDelayMap<0)= DoubleDelayMap(DoubleDelayMap<0) + 2*pi;
    DelayMap = DoubleDelayMap/2;
    
    AbsPhase1 = AngleMaps(:,:,curr_stims(1))-DelayMap;
    AbsPhase2 = AngleMaps(:,:,curr_stims(2))-DelayMap;
    
    AbsPhase1(sign(AbsPhase1) == stim_direction(curr_stims(1))) = AbsPhase1(sign(AbsPhase1) == stim_direction(curr_stims(1))) + 2*pi*-stim_direction(curr_stims(1)); %range=[-2*pi;0]
    AbsPhase1 = AbsPhase1*-stim_direction(curr_stims(1));
    AbsPhase2(sign(AbsPhase2) == stim_direction(curr_stims(2))) = AbsPhase2(sign(AbsPhase2) == stim_direction(curr_stims(2))) + 2*pi*-stim_direction(curr_stims(2)); %range=[-2*pi;0]
    AbsPhase2 = AbsPhase2*-stim_direction(curr_stims(2));
    
    meanAngleMaps = (AbsPhase1 + AbsPhase2)/2;
    meanAmpMaps = (AbsMaps(:,:,curr_stims(1)) + AbsMaps(:,:,curr_stims(2)))/2;
    
    angle_maps(:,:,curr_orientation) = meanAngleMaps;
    retinotopy_maps(:,:,curr_orientation) = meanAmpMaps.*exp(meanAngleMaps*sqrt(-1));
    
end

% old - who knows where this code is
% [f,jointImage] = bpViewComplexMaps(retinotopy_maps,[],[],[],[],[],'complex');
% close(f);
% % Get amplitude to plot alpha... does this make sense?
% alpha_maps = cellfun(@(x) mat2gray(sqrt(sum((x-1).^2,3))),jointImage,'uni',false);
% 
% % Plot retinotopy
% figure;
% for curr_orientation = 1:2
%     subplot(1,2,curr_orientation); hold on;
%     imagesc(avg_im); colormap(gray); caxis([0 6000]);
%     h = imshow(jointImage{curr_orientation});
%     set(h,'AlphaData',alpha_maps{curr_orientation});
% end

% Visual sign map

% 1) get gradient
[dhdx,dhdy] = imgradientxy(angle_maps(:,:,1));
[dvdx,dvdy] = imgradientxy(angle_maps(:,:,2));

% 2) get direction of gradient
[~,Vdir] = imgradient(dvdx,dvdy);
[~,Hdir] = imgradient(dhdx,dhdy);

% 3) get sin(difference in direction) if retinotopic, H/V should be
% orthogonal, so the closer the orthogonal the better (and get sign)
vfs = sind(Vdir-Hdir);

figure('Name',animal);
ax1 = axes;
subplot(1,2,1,ax1);
imagesc(vfs);
caxis([-1,1]);
axes(ax1); axis image off;
colormap(colormap_BlueWhiteRed)

ax2 = axes;
ax3 = axes;
subplot(1,2,2,ax2);
subplot(1,2,2,ax3);
h1 = imagesc(ax2,avg_im);
colormap(ax2,gray);
caxis(ax2,[0 prctile(avg_im(:),99.9)]);
h2 = imagesc(ax3,vfs);
colormap(ax3,colormap_BlueWhiteRed);
caxis([-1,1]);
set(ax2,'Visible','off');
axes(ax2); axis image off;
set(ax3,'Visible','off');
axes(ax3); axis image off;
set(h2,'AlphaData',mat2gray(abs(vfs))*0.5);
colormap(colormap_BlueWhiteRed)

%% Sparse noise retinotopy

% TO DO HERE: use pagefun or bsxfun on GPU for all pixels quickly?

[Uy,Ux,nSV] = size(U);

myScreenInfo.windowPtr = NaN; % so we can call the stimulus generation and it won't try to display anything
stimNum = 1;
ss = eval([Protocol.xfile(1:end-2) '(myScreenInfo, Protocol.pars(:,stimNum));']);
stim_screen = cat(3,ss.ImageTextures{:});
ny = size(stim_screen,1);
nx = size(stim_screen,2);

switch photodiode_type
    case 'flicker'
        % Check for case of mismatch between photodiode and stimuli:
        % odd number of stimuli, but one extra photodiode flip to come back down
        if mod(size(stim_screen,3),2) == 1 && ...
                length(photodiode.timestamps) == size(stim_screen,3) + 1;
            photodiode.timestamps(end) = [];
            photodiode.values(end) = [];
            warning('Odd number of stimuli, removed last photodiode');
        end
        
        % If there's still a mismatch, break
        if size(stim_screen,3) ~= length(photodiode.timestamps);
            warning([num2str(size(stim_screen,3)) ' stimuli, ', ...
                num2str(length(photodiode.timestamps)) ' photodiode pulses']);
            
            % Try to estimate which stim were missed by time difference
            photodiode_diff = diff(photodiode.timestamps);
            max_regular_diff_time = prctile(diff(photodiode.timestamps),99);
            skip_cutoff = max_regular_diff_time*2;
            photodiode_skip = find(photodiode_diff > skip_cutoff);
            est_n_pulse_skip = ceil(photodiode_diff(photodiode_skip)/max_regular_diff_time)-1;
            stim_skip = cell2mat(arrayfun(@(x) photodiode_skip(x):photodiode_skip(x)+est_n_pulse_skip(x)-1, ...
                1:length(photodiode_skip),'uni',false));
            
            if isempty(est_n_pulse_skip) || length(photodiode.timestamps) + sum(est_n_pulse_skip) ~= size(stim_screen,3)
                error('Can''t match photodiode events to stimuli')
            end
        end
        
        stim_times = photodiode.timestamps;
        
    case 'Steady'
        % If the photodiode is on steady: extrapolate the stim times
        if length(photodiode.timestamps) ~= 2
            error('Steady photodiode, but not 2 flips')
        end
        stim_duration = diff(photodiode.timestamps)/size(stim_screen,3);
        stim_times = linspace(photodiode.timestamps(1), ...
            photodiode.timestamps(2)-stim_duration,size(stim_screen,3))';
        
end

% Get average response to each stimulus
surround_window = [0.3,0.6];
framerate = 1./nanmedian(diff(frame_t));
surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);
response_n = nan(ny,nx);
response_grid = cell(ny,nx);
for px_y = 1:ny;
    for px_x = 1:nx;
        
        % Use first frame of dark or light stim
        align_stims = (stim_screen(px_y,px_x,2:end)~= 0) & ...
            (diff(stim_screen(px_y,px_x,:),[],3) ~= 0);
        align_times = stim_times(find(align_stims)+1);
        
        align_times = align_times(round(length(align_times)/2):end);
        
        response_n(px_y,px_x) = length(align_times);
        
        % Don't use times that fall outside of imaging
        align_times(align_times + surround_time(1) < frame_t(2) | ...
            align_times + surround_time(2) > frame_t(end)) = [];
        
        % Get stim-aligned responses, 2 choices:
        
        % 1) Interpolate times (slow - but supersamples so better)
        %         align_surround_times = bsxfun(@plus, align_times, surround_time);
        %         peri_stim_v = permute(mean(interp1(frame_t,fV',align_surround_times),1),[3,2,1]);
        
        % 2) Use closest frames to times (much faster - not different)
        align_surround_times = bsxfun(@plus, align_times, surround_time);
        frame_edges = [frame_t,frame_t(end)+1/framerate];
        align_frames = discretize(align_surround_times,frame_edges);
        
        align_frames(any(isnan(align_frames),2),:) = [];
        
        peri_stim_v = permute(nanmean(reshape(fV(:,align_frames)', ...
            size(align_frames,1),size(align_frames,2),[]),1),[3,2,1]);
        
        % Save V's
        response_grid{px_y,px_x} = peri_stim_v;
        
    end
end

% Get position preference for every pixel
U_downsample_factor = 3;
resize_scale = 3;
filter_sigma = (resize_scale*1.5);

% Downsample U
use_u_y = 1:Uy;
Ud = imresize(U(use_u_y,:,:),1/U_downsample_factor,'bilinear');

% Convert V responses to pixel responses
use_svs = 1:size(U,3);
response_mean = cell2mat(cellfun(@(x) nanmean(x,2),response_grid(:),'uni',false)');
stim_im_px = reshape(permute(svdFrameReconstruct(Ud(:,:,use_svs),response_mean(use_svs,:)),[3,1,2]),ny,nx,[]);

% Upsample each pixel's response map and find maximum
gauss_filt = fspecial('gaussian',[ny,nx],filter_sigma);
stim_im_smoothed = imfilter(imresize(stim_im_px,resize_scale,'bilinear'),gauss_filt);
[mc,mi] = max(reshape(stim_im_smoothed,[],size(stim_im_px,3)),[],1);
[m_y,m_x] = ind2sub(size(stim_im_smoothed),mi);
m_yr = reshape(m_y,size(Ud,1),size(Ud,2));
m_xr = reshape(m_x,size(Ud,1),size(Ud,2));

% Plot retinotopy in different ways:
figure;

% Plot retinotopy by color
screen_pos_col = nan(ny,nx,3);
[screen_pos_col(:,:,1),screen_pos_col(:,:,2)] = ...
    meshgrid(mat2gray(1:nx),mat2gray(1:ny));
screen_pos_col(:,:,3) = 1- ...
    (screen_pos_col(:,:,1)+screen_pos_col(:,:,2))/2;
screen_pos_col_upsample = reshape(imresize(screen_pos_col,resize_scale),[],3);

retinotopy_colormap = reshape(screen_pos_col_upsample(mi,:),size(Ud,1),size(Ud,2),3);

subplot(3,1,1);
h = imagesc(retinotopy_colormap);
set(gca,'YDir','Normal');
axis off;
title('Positional retinotopy');

% Plot vertical/horizontal lines
stripe_spacing = 5;
vert_stripes = imgaussfilt(double(mod(m_xr-1,stripe_spacing) == 0),1);
horz_stripes = imgaussfilt(double(mod(m_yr-1,stripe_spacing) == 0),1);
subplot(3,1,2);
imagesc(padarray(cat(3,mat2gray(vert_stripes),mat2gray(horz_stripes)),[0,0,1],0,'post'));
set(gca,'YDir','normal')
axis off;
title('Iso-horizontal/vertical retinotopy');

% (get all vertical/horizontal stripes for plotting manually)
vert_response = nan(size(Ud,1),size(Ud,2),size(stim_im_smoothed,2));
for i = 1:size(stim_im_smoothed,2)
    vert_response(:,:,i) = m_xr == i;
end
vert_response = imgaussfilt(vert_response,1);

horz_response = nan(size(Ud,1),size(Ud,2),size(stim_im_smoothed,1));
for i = 1:size(stim_im_smoothed,1)
    horz_response(:,:,i) = m_yr == i;
end
horz_response = imgaussfilt(horz_response,1);

% Calculate and plot sign map (do this just with dot product between horz / vert grad?)

% 1) get gradient direction
[Vmag,Vdir] = imgradient(imgaussfilt(m_yr,1));
[Hmag,Hdir] = imgradient(imgaussfilt(m_xr,1));

% 3) get sin(difference in direction) if retinotopic, H/V should be
% orthogonal, so the closer the orthogonal the better (and get sign)
angle_diff = sind(Vdir-Hdir);

subplot(3,1,3);
imagesc(imgaussfilt(angle_diff,1));
set(gca,'YDir','normal')
axis off;
title('Visual sign field');

%%% TO DO: 2 dimensions, sign and orientation?
% average direction between gradients
retinotopy_dir = (Vdir+Hdir)/2;

%%% TO DO NEXT: use MakeSeparable function instead (just SVD, but gets the
%%% temporal and spatial kernel given first temporal without other kinds of stuff)
%%% (this doesn't really look like it's worth it, just use 200-300 ms)
%%% potentially get rid of not visually-related signals by something other
%%% than ignoring early SVs (maybe visually responsive SVs?)
%%% Also try processing all pixels w/ gpu to get reasonable speed?






%%%% FOR TESTING PURPOSES


% % (After concatenating the response images to make response map)
% %im = reshape(reshape(stim_im_px,[],size(stim_im_px,3))',size(Ud,1),size(Ud,2),[]);
% a = reshape(roi.trace,ny,nx);
% % upsample and filter
% upsample_factor = 10;
% b = imgaussfilt(imresize(a,upsample_factor,'bilinear'),upsample_factor*filter_sigma);
% figure;imagesc(b);colormap(redblue);
% c = caxis;
% caxis([-max(abs(c)),max(abs(c))]);


%% Sparse noise subtracting surround responses

[Uy,Ux,nSV] = size(U);

myScreenInfo.windowPtr = NaN; % so we can call the stimulus generation and it won't try to display anything
stimNum = 1;
ss = eval([Protocol.xfile(1:end-2) '(myScreenInfo, Protocol.pars(:,stimNum));']);
stim_screen = cat(3,ss.ImageTextures{:});
ny = size(stim_screen,1);
nx = size(stim_screen,2);

switch photodiode_type
    case 'flicker'
        % Check for case of mismatch between photodiode and stimuli:
        % odd number of stimuli, but one extra photodiode flip to come back down
        if mod(size(stim_screen,3),2) == 1 && ...
                length(photodiode.timestamps) == size(stim_screen,3) + 1;
            photodiode.timestamps(end) = [];
            photodiode.values(end) = [];
            warning('Odd number of stimuli, removed last photodiode');
        end
        
        % If there's still a mismatch, break
        if size(stim_screen,3) ~= length(photodiode.timestamps);
            warning([num2str(size(stim_screen,3)) ' stimuli, ', ...
                num2str(length(photodiode.timestamps)) ' photodiode pulses']);
            
            % Try to estimate which stim were missed by time difference
            photodiode_diff = diff(photodiode.timestamps);
            max_regular_diff_time = prctile(diff(photodiode.timestamps),99);
            skip_cutoff = max_regular_diff_time*2;
            photodiode_skip = find(photodiode_diff > skip_cutoff);
            est_n_pulse_skip = ceil(photodiode_diff(photodiode_skip)/max_regular_diff_time)-1;
            stim_skip = cell2mat(arrayfun(@(x) photodiode_skip(x):photodiode_skip(x)+est_n_pulse_skip(x)-1, ...
                1:length(photodiode_skip),'uni',false));
            
            if isempty(est_n_pulse_skip) || length(photodiode.timestamps) + sum(est_n_pulse_skip) ~= size(stim_screen,3)
                error('Can''t match photodiode events to stimuli')
            end
        end
        
        stim_times = photodiode.timestamps;
        
    case 'Steady'
        % If the photodiode is on steady: extrapolate the stim times
        if length(photodiode.timestamps) ~= 2
            error('Steady photodiode, but not 2 flips')
        end
        stim_duration = diff(photodiode.timestamps)/size(stim_screen,3);
        stim_times = linspace(photodiode.timestamps(1), ...
            photodiode.timestamps(2)-stim_duration,size(stim_screen,3));
        
end

% Get average response to each stimulus
surround_window = [0.2,0.3];
framerate = 1./nanmedian(diff(frame_t));
surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);
response_n = nan(ny,nx);
response_grid = cell(ny,nx);
for px_y = 1:ny;
    for px_x = 1:nx;
        
        % Use first frame of dark or light stim
        align_stims = (stim_screen(px_y,px_x,2:end)~= 0) & ...
            (diff(stim_screen(px_y,px_x,:),[],3) ~= 0);
        align_times = photodiode.timestamps(find(align_stims)+1);
        response_n(px_y,px_x) = length(align_times);
        
        % Don't use times that fall outside of imaging
        align_times(align_times + surround_time(1) < frame_t(2) | ...
            align_times + surround_time(2) > frame_t(end)) = [];
        
        % Get stim-aligned responses, 2 choices:
        
        % 1) Interpolate times (slow - but supersamples so better)
        %         align_surround_times = bsxfun(@plus, align_times, surround_time);
        %         peri_stim_v = permute(mean(interp1(frame_t,fV',align_surround_times),1),[3,2,1]);
        
        % 2) Use closest frames to times (much faster - not different)
        align_surround_times = bsxfun(@plus, align_times, surround_time);
        frame_edges = [frame_t,frame_t(end)+1/framerate];
        align_frames = discretize(align_surround_times,frame_edges);
        
        align_frames(any(isnan(align_frames),2),:) = [];
        
        peri_stim_v = permute(nanmean(reshape(fV(:,align_frames)', ...
            size(align_frames,1),size(align_frames,2),[]),1),[3,2,1]);
        
        % Save V's
        response_grid{px_y,px_x} = peri_stim_v;
        
    end
end

% Get position preference for every pixel

U_downsample_factor = 1;
resize_scale = 5;
filter_sigma = (resize_scale*1.5);

% Downsample U
%use_u_y = 1:Uy;
use_u_y = 20:200;
Ud = imresize(U(use_u_y,:,:),1/U_downsample_factor,'bilinear');

% Convert V responses to pixel responses
%use_svs = 30:size(U,3);
use_svs = 10:size(U,3);
response_mean = cell2mat(cellfun(@(x) nanmean(x,2),response_grid(:),'uni',false)');

response_grid_mean = cellfun(@(x) nanmean(x,2),response_grid,'uni',false);
response_mean_subsurr = nan(size(response_mean));

for i = 1:numel(response_grid)
    
    curr_px = false(size(response_grid));
    curr_px(i) = true;
    
    surround_radius = 3;
    se_surround = strel('disk',surround_radius,0);
    se_inside = strel('disk',surround_radius-1,0);
    
    surround_px = imdilate(curr_px,se_surround) & ~imdilate(curr_px,se_inside);
    
    surround_response = nanmean(horzcat(response_grid_mean{surround_px}),2);
    
    response_mean_subsurr(:,i) = response_grid_mean{i} - surround_response;
    
end

stim_im_px = reshape(permute(svdFrameReconstruct(Ud(:,:,use_svs),response_mean_subsurr(use_svs,:)),[3,1,2]),ny,nx,[]);

U_downsample_factor = 1;
resize_scale = 5;
filter_sigma = (resize_scale*1.5);

% Upsample each pixel's response map and find maximum
gauss_filt = fspecial('gaussian',[ny,nx],filter_sigma);
stim_im_smoothed = imfilter(imresize(stim_im_px,resize_scale,'bilinear'),gauss_filt);
stim_im_skew = reshape(skewness(abs(reshape(stim_im_smoothed,[],size(stim_im_px,3))),[],1),size(Ud,1),size(Ud,2));
[mc,mi] = max(reshape(stim_im_smoothed,[],size(stim_im_px,3)),[],1);
[m_y,m_x] = ind2sub(size(stim_im_smoothed),mi);
m_yr = reshape(m_y,size(Ud,1),size(Ud,2));
m_xr = reshape(m_x,size(Ud,1),size(Ud,2));

% Plot retinotopy in different ways:
figure;

% Plot retinotopy by color
screen_pos_col = nan(ny,nx,3);
[screen_pos_col(:,:,1),screen_pos_col(:,:,2)] = ...
    meshgrid(mat2gray(1:nx),mat2gray(1:ny));
screen_pos_col(:,:,3) = 1- ...
    (screen_pos_col(:,:,1)+screen_pos_col(:,:,2))/2;
screen_pos_col_upsample = reshape(imresize(screen_pos_col,resize_scale),[],3);

retinotopy_colormap = reshape(screen_pos_col_upsample(mi,:),size(Ud,1),size(Ud,2),3);

subplot(3,1,1);
h = imagesc(retinotopy_colormap);
set(h,'AlphaData',mat2gray(stim_im_skew));
set(gca,'YDir','Normal');
axis off;
title('Positional retinotopy');

% Plot isovertical/horizontal lines
stripe_spacing = 5;
vert_stripes = imgaussfilt(double(mod(m_xr-1,stripe_spacing) == 0),1);
horz_stripes = imgaussfilt(double(mod(m_yr-1,stripe_spacing) == 0),1);
subplot(3,1,2);
imagesc(padarray(cat(3,mat2gray(vert_stripes),mat2gray(horz_stripes)),[0,0,1],0,'post'));
set(gca,'YDir','normal')
axis off;
title('Iso-horizontal/vertical retinotopy');

% (get all vertical/horizontal stripes for plotting manually)
vert_response = nan(size(Ud,1),size(Ud,2),size(stim_im_smoothed,2));
for i = 1:size(stim_im_smoothed,2)
    vert_response(:,:,i) = m_xr == i;
end
vert_response = imgaussfilt(vert_response,1);

horz_response = nan(size(Ud,1),size(Ud,2),size(stim_im_smoothed,1));
for i = 1:size(stim_im_smoothed,1)
    horz_response(:,:,i) = m_yr == i;
end
horz_response = imgaussfilt(horz_response,1);

% Calculate and plot sign map (do this just with dot product between horz / vert grad?)

% 1) get gradient
[dhdx,dhdy] = imgradientxy(m_xr);
[dvdx,dvdy] = imgradientxy(m_yr);

% 2) get direction of gradient
[Vmag,Vdir] = imgradient(dvdx,dvdy);
[Hmag,Hdir] = imgradient(dhdx,dhdy);

% 3) get sin(difference in direction) if retinotopic, H/V should be
% orthogonal, so the closer the orthogonal the better (and get sign)
angle_diff = sind(Vdir-Hdir);

subplot(3,1,3);
imagesc(imgaussfilt(angle_diff,1));
set(gca,'YDir','normal')
axis off;
title('Visual sign field');

%% Sparse noise fitting with gaussian / forward model

[Uy,Ux,nSV] = size(U);

myScreenInfo.windowPtr = NaN; % so we can call the stimulus generation and it won't try to display anything
stimNum = 1;
ss = eval([Protocol.xfile(1:end-2) '(myScreenInfo, Protocol.pars(:,stimNum));']);
stim_screen = cat(3,ss.ImageTextures{:});
ny = size(stim_screen,1);
nx = size(stim_screen,2);

switch photodiode_type
    case 'flicker'
        % Check for case of mismatch between photodiode and stimuli:
        % odd number of stimuli, but one extra photodiode flip to come back down
        if mod(size(stim_screen,3),2) == 1 && ...
                length(photodiode.timestamps) == size(stim_screen,3) + 1;
            photodiode.timestamps(end) = [];
            photodiode.values(end) = [];
            warning('Odd number of stimuli, removed last photodiode');
        end
        
        % If there's still a mismatch, break
        if size(stim_screen,3) ~= length(photodiode.timestamps);
            warning([num2str(size(stim_screen,3)) ' stimuli, ', ...
                num2str(length(photodiode.timestamps)) ' photodiode pulses']);
            
            % Try to estimate which stim were missed by time difference
            photodiode_diff = diff(photodiode.timestamps);
            max_regular_diff_time = prctile(diff(photodiode.timestamps),99);
            skip_cutoff = max_regular_diff_time*2;
            photodiode_skip = find(photodiode_diff > skip_cutoff);
            est_n_pulse_skip = ceil(photodiode_diff(photodiode_skip)/max_regular_diff_time)-1;
            stim_skip = cell2mat(arrayfun(@(x) photodiode_skip(x):photodiode_skip(x)+est_n_pulse_skip(x)-1, ...
                1:length(photodiode_skip),'uni',false));
            
            if isempty(est_n_pulse_skip) || length(photodiode.timestamps) + sum(est_n_pulse_skip) ~= size(stim_screen,3)
                error('Can''t match photodiode events to stimuli')
            end
        end
        
        stim_times = photodiode.timestamps;
        
    case 'Steady'
        % If the photodiode is on steady: extrapolate the stim times
        if length(photodiode.timestamps) ~= 2
            error('Steady photodiode, but not 2 flips')
        end
        stim_duration = diff(photodiode.timestamps)/size(stim_screen,3);
        stim_times = linspace(photodiode.timestamps(1), ...
            photodiode.timestamps(2)-stim_duration,size(stim_screen,3));
        
end


% Get average response to each stimulus
surround_window = [0,3];
framerate = 1./nanmedian(diff(frame_t));
surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);
response_n = nan(ny,nx);
response_grid = cell(ny,nx);
for px_y = 1:ny;
    for px_x = 1:nx;
        
        % Use first frame of dark or light stim
        align_stims = (stim_screen(px_y,px_x,2:end)~= 0) & ...
            (diff(stim_screen(px_y,px_x,:),[],3) ~= 0);
        align_times = stim_times(find(align_stims)+1);
        
        align_times = align_times(round(length(align_times)/2):end);
        
        response_n(px_y,px_x) = length(align_times);
        
        % Don't use times that fall outside of imaging
        align_times(align_times + surround_time(1) < frame_t(2) | ...
            align_times + surround_time(2) > frame_t(end)) = [];
        
        % Get stim-aligned responses, 2 choices:
        
        % 1) Interpolate times (slow - but supersamples so better)
        %         align_surround_times = bsxfun(@plus, align_times, surround_time);
        %         peri_stim_v = permute(mean(interp1(frame_t,fV',align_surround_times),1),[3,2,1]);
        
        % 2) Use closest frames to times (much faster - not different)
        align_surround_times = bsxfun(@plus, align_times, surround_time);
        frame_edges = [frame_t,frame_t(end)+1/framerate];
        align_frames = discretize(align_surround_times,frame_edges);
        
        align_frames(any(isnan(align_frames),2),:) = [];
        
        peri_stim_v = permute(nanmean(reshape(fV(:,align_frames)', ...
            size(align_frames,1),size(align_frames,2),[]),1),[3,2,1]);
        
        % Save V's
        response_grid{px_y,px_x} = peri_stim_v;
        
    end
end

% Get position preference for every pixel

U_downsample_factor = 3;
resize_scale = 1;
filter_sigma = (resize_scale*1.5);

% Downsample U
%use_u_y = 1:Uy;
use_u_y = 20:200;
Ud = imresize(U(use_u_y,:,:),1/U_downsample_factor,'bilinear');

% Convert V responses to pixel responses
use_svs = 30:size(U,3);
response_mean = cell2mat(cellfun(@(x) nanmean(x,2),response_grid(:),'uni',false)');
stim_im_px = reshape(permute(svdFrameReconstruct(Ud(:,:,use_svs),response_mean(use_svs,:)),[3,1,2]),ny,nx,[]);

% Upsample each pixel's response map and find maximum
gauss_filt = fspecial('gaussian',[ny,nx],filter_sigma);
stim_im_smoothed = imfilter(imresize(stim_im_px,resize_scale,'bilinear'),gauss_filt);
stim_im_skew = reshape(skewness(abs(reshape(stim_im_smoothed,[],size(stim_im_px,3))),[],1),size(Ud,1),size(Ud,2));
[mc,mi] = max(reshape(stim_im_smoothed,[],size(stim_im_px,3)),[],1);
[m_y,m_x] = ind2sub(size(stim_im_smoothed),mi);
m_yr = reshape(m_y,size(Ud,1),size(Ud,2));
m_xr = reshape(m_x,size(Ud,1),size(Ud,2));


% Make temporal kernel for response

% get grid of frames when stim is on screen
stim_frames = zeros(ny,nx,size(V,2));
for px_y = 1:ny
    for px_x = 1:nx
        use_stims = (stim_screen(px_y,px_x,2:end)~= 0) & ...
            (diff(stim_screen(px_y,px_x,:),[],3) ~= 0);
        stim_times = photodiode.timestamps(find(use_stims)+1);
        
        frame_edges = [frame_t,frame_t(end)+1/framerate];
        stim_frames(px_y,px_x,discretize(stim_times,frame_edges)) = 1;
    end
end
stim_frames = reshape(stim_frames,ny*nx,[]);


% at the moment, this is just the dumbest one there is
temporal_kernel = [zeros(1,5),ones(1,10)];

stim_trace = conv2(stim_frames,temporal_kernel);
stim_trace = stim_trace(:,1:size(stim_frames,2));

Ud_reshape = reshape(Ud,[],size(Ud,3));

gauss_fit_x = nan(size(Ud,1),size(Ud,2));
gauss_fit_y = nan(size(Ud,1),size(Ud,2));
gauss_fit_sigma = nan(size(Ud,1),size(Ud,2));

for curr_px = 1:size(Ud_reshape,1)
    
    px_trace = permute(Ud_reshape(curr_px,:)*fV,[2,3,1]);
    
    x_start = m_x(curr_px);
    y_start = m_y(curr_px);
    
    gauss_fit_params = fminsearch(@(gauss_fit_params) AP_2d_gaussian_fit_cost ...
        (gauss_fit_params,[ny,nx],stim_trace,px_trace),[x_start,y_start,1]);
    
    gauss_fit_x(curr_px) = gauss_fit_params(1);
    gauss_fit_y(curr_px) = gauss_fit_params(2);
    gauss_fit_sigma(curr_px) = gauss_fit_params(3);
    
    disp(curr_px/size(Ud_reshape,1));
    
end

% Sign map
% 1) get gradient
[dhdx,dhdy] = imgradientxy(gauss_fit_x);
[dvdx,dvdy] = imgradientxy(gauss_fit_y);

% 2) get direction of gradient
[Vmag,Vdir] = imgradient(dvdx,dvdy);
[Hmag,Hdir] = imgradient(dhdx,dhdy);

% 3) get sin(difference in direction) if retinotopic, H/V should be
% orthogonal, so the closer the orthogonal the better (and get sign)
angle_diff = sind(Vdir-Hdir);

figure;
imagesc(imgaussfilt(angle_diff,1));
set(gca,'YDir','normal')
axis off;
title('Visual sign field');

%% (OTHER TEST: FIT GAUSSIAN WITHOUT RESHAPE)

px_gauss = nan(size(stim_im_px,3),6);
for curr_px = 1:size(stim_im_px,3)
    px_gauss(curr_px,:) = fit2dGaussRF(double(stim_im_px(:,:,curr_px)));
    disp(curr_px/size(stim_im_px,3));
end

x_center = reshape(px_gauss(:,2),size(Ud,1),size(Ud,2));
y_center = reshape(px_gauss(:,4),size(Ud,1),size(Ud,2));

% 1) get gradient
[dhdx,dhdy] = imgradientxy(x_center);
[dvdx,dvdy] = imgradientxy(y_center);

% 2) get direction of gradient
[Vmag,Vdir] = imgradient(dvdx,dvdy);
[Hmag,Hdir] = imgradient(dhdx,dhdy);

% 3) get sin(difference in direction) if retinotopic, H/V should be
% orthogonal, so the closer the orthogonal the better (and get sign)
angle_diff = sind(Vdir-Hdir);

imagesc(imgaussfilt(angle_diff,1));
set(gca,'YDir','normal')
axis off;
title('Visual sign field');

% try on gpu?
px_gauss = pagefun(@fit2dGaussRF,gpuArray(double(stim_im_px)));


%% XY Pos localize spot

refresh_rate_cutoff = 1/10;
stim_onsets = photodiode_onsets( ...
    [1;find(diff(photodiode_onsets) > refresh_rate_cutoff) + 1]);
stim_offsets = photodiode_offsets([arrayfun(@(x) find(photodiode_offsets < x,1,'last'),stim_onsets(2:end));end]);

if length(stim_onsets) ~= numel(Protocol.seqnums)
    error('Mismatching number of stims and photodiode events')
end

stimIDs = zeros(size(stim_onsets));
for q = 1:size(Protocol.seqnums,1)
    stimIDs(Protocol.seqnums(q,:)) = q;
end

stim_duration = Protocol.pars(strcmp('dur',Protocol.parnames),:)/10;
stim_x = Protocol.pars(strcmp('xfocus',Protocol.parnames),:)/10;
stim_y = Protocol.pars(strcmp('yfocus',Protocol.parnames),:)/10;

n_stim = size(Protocol.pars,2);

% Get average response to each stim position
framerate = 1./nanmedian(diff(frame_t));

surround_window = [0,max(stim_duration)];
surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);

peri_stim_v = nan(size(U,3),length(surround_time),n_stim);
for curr_condition = 1:n_stim
    
    use_stims = find(stimIDs == curr_condition);
    use_stim_onsets = stim_onsets(use_stims);
    
    use_stim_onsets_periods = use_stim_onsets;
    
    stim_surround_times = bsxfun(@plus, use_stim_onsets_periods(:), surround_time);
    peri_stim_v(:,:,curr_condition) = permute(mean(interp1(frame_t,fV',stim_surround_times),1),[3,2,1]);
    
end

% Define "baseline" as non-stim times (in this convention, odd = during
% stim and even = during no-stim)
% give leeway to stim to allow for hemo bounceback
stim_offsets_leeway = stim_offsets + 0.5;
no_stim_frames = mod(discretize(frame_t,reshape([stim_onsets';stim_offsets_leeway';],[],1)),2) == 0;
baseline_mean = mean(svdFrameReconstruct(U,fV(:,no_stim_frames)),3);
baseline_std = std(svdFrameReconstruct(U,fV(:,no_stim_frames)),[],3);

stim_im = permute(cell2mat(permute(arrayfun(@(x) ...
    svdFrameReconstruct(U,nanmean(peri_stim_v(:,:,x),2)),1:n_stim,'uni',false),[1,3,4,2])),[1,2,4,3]);
stim_im_norm = bsxfun(@rdivide,stim_im,baseline_std);

% Subtract 1/4 of 4 surrounding stim to get localized spot
% (just hardcode which stimuli these are for now)
r_surround_stim = [1,3,4,5];
r_center_stim = 2;

l_surround_stim = [6,8,9,10];
l_center_stim = 7;

r_localize_response = stim_im_norm(:,:,r_center_stim) - nanmean(stim_im_norm(:,:,r_surround_stim),3);
l_localize_response = stim_im_norm(:,:,l_center_stim) - nanmean(stim_im_norm(:,:,l_surround_stim),3);

figure;
subplot(1,2,1);
imagesc(r_localize_response);colormap(gray);
set(gca,'YDir','normal');
axis off
subplot(1,2,2);
imagesc(l_localize_response);colormap(gray);
set(gca,'YDir','normal');
axis off

%
%
% Ur = gpuArray(reshape(U, size(U,1)*size(U,2), size(U,3)));
% fVg = gpuArray(fV);
%
% a = pagefun(@(x) x*fVg,permute(Ur,[3,2,1]));


%% Sparse noise retinotopy (bootstrap across trials to get less noise)

[Uy,Ux,nSV] = size(U);

myScreenInfo.windowPtr = NaN; % so we can call the stimulus generation and it won't try to display anything
stimNum = 1;
ss = eval([Protocol.xfile(1:end-2) '(myScreenInfo, Protocol.pars(:,stimNum));']);
stim_screen = cat(3,ss.ImageTextures{:});
ny = size(stim_screen,1);
nx = size(stim_screen,2);

% This is garbage: higher threshold for this photodiode flip on flicker
photodiode_thresh = 4;
photodiode_trace = Timeline.rawDAQData(stimScreen_on,photodiode_idx) > photodiode_thresh;
% (medfilt because photodiode can be intermediate value when backlight
% coming on)
photodiode_trace_medfilt = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
    photodiode_idx),1) > photodiode_thresh;
photodiode_flip = find((~photodiode_trace_medfilt(1:end-1) & photodiode_trace_medfilt(2:end)) | ...
    (photodiode_trace_medfilt(1:end-1) & ~photodiode_trace_medfilt(2:end)))+1;
photodiode_flip_times = stimScreen_on_t(photodiode_flip)';

if strcmp(photodiode_type,'flicker') && length(photodiode_flip_times) ~= size(stim_screen,3)
    warning('Mismatching stim times: interpolating start/end');
    photodiode_flip_times = photodiode_flip_times([1,end]);
    photodiode_type = 'steady';
end

switch lower(photodiode_type)
    case 'flicker'
        % Check for case of mismatch between photodiode and stimuli:
        % odd number of stimuli, but one extra photodiode flip to come back down
        if mod(size(stim_screen,3),2) == 1 && ...
                length(photodiode_flip_times) == size(stim_screen,3) + 1
            photodiode_flip_times(end) = [];
            warning('Odd number of stimuli, removed last photodiode');
        end
        
        % If there's still a mismatch, attempt to fix
        if size(stim_screen,3) ~= length(photodiode_flip_times)
            warning([num2str(size(stim_screen,3)) ' stimuli, ', ...
                num2str(length(photodiode_flip_times)) ' photodiode pulses']);
            
            % Try to estimate which stim were missed by time difference
            photodiode_diff = diff(photodiode_flip_times);
            max_regular_diff_time = prctile(diff(photodiode_flip_times),99);
            skip_cutoff = max_regular_diff_time*2;
            photodiode_skip = find(photodiode_diff > skip_cutoff);
            est_n_pulse_skip = ceil(photodiode_diff(photodiode_skip)/max_regular_diff_time)-1;
            stim_skip = cell2mat(arrayfun(@(x) photodiode_skip(x):photodiode_skip(x)+est_n_pulse_skip(x)-1, ...
                1:length(photodiode_skip),'uni',false));
            
            if isempty(est_n_pulse_skip) || length(photodiode_flip_times) + sum(est_n_pulse_skip) ~= size(stim_screen,3)
                error('Can''t match photodiode events to stimuli')
            end
        end
        
        stim_times = photodiode_flip_times;
        
    case 'steady'
        % If the photodiode is on steady: extrapolate the stim times
        if length(photodiode_flip_times) ~= 2
            error('Steady photodiode, but not 2 flips')
        end
        stim_duration = diff(photodiode_flip_times)/size(stim_screen,3);
        stim_times = linspace(photodiode_flip_times(1), ...
            photodiode_flip_times(2)-stim_duration,size(stim_screen,3))';
        
end

% Get average response to each stimulus
surround_window = [0.1,0.5]; % 6s = [0.1,0.5], 6f = [0.05,0.2]
framerate = 1./nanmedian(diff(frame_t));
surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);
response_n = nan(ny,nx);
response_grid = cell(ny,nx);
for px_y = 1:ny
    for px_x = 1:nx
        
        % Use first frame of dark or light stim
        align_stims = (stim_screen(px_y,px_x,2:end)~= 0) & ...
            (diff(stim_screen(px_y,px_x,:),[],3) ~= 0);
        align_times = stim_times(find(align_stims)+1);
        
%         error('What the hell is this? throw away half of the data?')
%         align_times = align_times(round(length(align_times)/2):end);
        
        response_n(px_y,px_x) = length(align_times);
        
        % Don't use times that fall outside of imaging
        align_times(align_times + surround_time(1) < frame_t(2) | ...
            align_times + surround_time(2) > frame_t(end)) = [];
        
        % Get stim-aligned responses, 2 choices:
        
        % 1) Interpolate times (slow - but supersamples so better)
        %         align_surround_times = bsxfun(@plus, align_times, surround_time);
        %         peri_stim_v = permute(mean(interp1(frame_t,fV',align_surround_times),1),[3,2,1]);
        
        % 2) Use closest frames to times (much faster - not different)
        align_surround_times = bsxfun(@plus, align_times, surround_time);
        frame_edges = [frame_t,frame_t(end)+1/framerate];
        align_frames = discretize(align_surround_times,frame_edges);
        
        align_frames(any(isnan(align_frames),2),:) = [];
        
        % Define the peri-stim V's as subtracting first frame (baseline)
        peri_stim_v = bsxfun(@minus, ...
            reshape(fV(:,align_frames)',size(align_frames,1),size(align_frames,2),[]), ...
            reshape(fV(:,align_frames(:,1))',size(align_frames(:,1),1),size(align_frames(:,1),2),[]));

        mean_peri_stim_v = permute(mean(peri_stim_v,2),[3,1,2]);
        
        % Save V's
        response_grid{px_y,px_x} = mean_peri_stim_v;
        
    end
end

% Get position preference for every pixel
U_downsample_factor = 1; %2 if max method
screen_resize_scale = 1; %3 if max method
filter_sigma = (screen_resize_scale*2);

% Downsample U
use_u_y = 1:Uy;
Ud = imresize(U(use_u_y,:,:),1/U_downsample_factor,'bilinear');

% Convert V responses to pixel responses
use_svs = 1:size(U,3);
n_boot = 10;

response_mean_boostrap = cellfun(@(x) bootstrp(n_boot,@mean,x')',response_grid,'uni',false);

% (to split trials instead of bootstrap)
%split_trials = cellfun(@(x) shake(discretize(1:size(x,2),round(linspace(1,size(x,2),n_boot+1)))),response_grid,'uni',false);
%response_mean_boostrap = cellfun(@(x,y) grpstats(x',y','mean')',response_grid,split_trials,'uni',false);
use_method = 'com'; % max or com
vfs_boot = nan(size(Ud,1),size(Ud,2),n_boot);
for curr_boot = 1:n_boot
    
    response_mean = cell2mat(cellfun(@(x) x(:,curr_boot),response_mean_boostrap(:),'uni',false)');
    stim_im_px = reshape(permute(svdFrameReconstruct(Ud(:,:,use_svs),response_mean(use_svs,:)),[3,1,2]),ny,nx,[]);
    gauss_filt = fspecial('gaussian',[ny,nx],filter_sigma);
    stim_im_smoothed = imfilter(imresize(stim_im_px,screen_resize_scale,'bilinear'),gauss_filt);

    switch use_method
        case 'max'
            % Upsample each pixel's response map and find maximum            
            [~,mi] = max(reshape(stim_im_smoothed,[],size(stim_im_px,3)),[],1);
            [m_y,m_x] = ind2sub(size(stim_im_smoothed),mi);
            m_yr = reshape(m_y,size(Ud,1),size(Ud,2));
            m_xr = reshape(m_x,size(Ud,1),size(Ud,2));
            
        case 'com'
            % Conversely, do COM on original^2
            [xx,yy] = meshgrid(1:size(stim_im_smoothed,2),1:size(stim_im_smoothed,1));
            m_xr = reshape(sum(sum(bsxfun(@times,stim_im_smoothed.^2,xx),1),2)./sum(sum(stim_im_smoothed.^2,1),2),size(Ud,1),size(Ud,2));
            m_yr = reshape(sum(sum(bsxfun(@times,stim_im_smoothed.^2,yy),1),2)./sum(sum(stim_im_smoothed.^2,1),2),size(Ud,1),size(Ud,2));
    end
    
    % Calculate and plot sign map (dot product between horz & vert gradient)
    
    % 1) get gradient direction
    [~,Vdir] = imgradient(imgaussfilt(m_yr,1));
    [~,Hdir] = imgradient(imgaussfilt(m_xr,1));
    
    % 3) get sin(difference in direction) if retinotopic, H/V should be
    % orthogonal, so the closer the orthogonal the better (and get sign)
    angle_diff = sind(Vdir-Hdir);
    angle_diff(isnan(angle_diff)) = 0;
    
    vfs_boot(:,:,curr_boot) = angle_diff;
    
    disp(curr_boot);
    
end

vfs_median = imgaussfilt(nanmedian(vfs_boot,3),2);

figure('Name',animal);
ax1 = axes;
subplot(1,2,1,ax1);
imagesc(vfs_median);
caxis([-1,1]);
axes(ax1); axis image off;
colormap(colormap_BlueWhiteRed)

ax2 = axes;
ax3 = axes;
subplot(1,2,2,ax2);
subplot(1,2,2,ax3);
h1 = imagesc(ax2,avg_im(use_u_y,:));
colormap(ax2,gray);
caxis(ax2,[0 prctile(avg_im(:),99.9)]);
h2 = imagesc(ax3,vfs_median);
colormap(ax3,colormap_BlueWhiteRed);
caxis([-1,1]);
set(ax2,'Visible','off');
axes(ax2); axis image off;
set(ax3,'Visible','off');
axes(ax3); axis image off;
set(h2,'AlphaData',mat2gray(abs(vfs_median))*0.5);
colormap(colormap_BlueWhiteRed)
colormap(ax2,gray);

% vfs_cutoff = 1.5*std(vfs_mean(:));
% vis_thresh = abs(vfs_mean) > vfs_cutoff;
% vis_boundaries = cellfun(@(x) x*U_downsample_factor,bwboundaries(vis_thresh),'uni',false);
% 
% figure; hold on;
% set(gca,'YDir','reverse');
% imagesc(avg_im); colormap(gray);
% caxis([0 prctile(avg_im(:),95)]);
% for i = 1:length(vis_boundaries)
%     plot(vis_boundaries{i}(:,2),vis_boundaries{i}(:,1),'m','linewidth',1);
% end
% ylim([0,size(avg_im,1)]);
% xlim([0,size(avg_im,2)]);
% axis off;


%% Signals: sparse noise retinotopy (bootstrap across trials to get less noise)

[Uy,Ux,nSV] = size(U);

ny = size(block.events.stimuliOnValues,1);
nx = size(block.events.stimuliOnValues,2)/ ...
    size(block.events.stimuliOnTimes,2);
stim_screen = reshape(block.events.stimuliOnValues,ny,nx,[]);

% Get average response to each stimulus
surround_window = [0,0.5]; % 6s = [0.1,0.5], 6f = [0.05,0.2]
framerate = 1./nanmean(diff(frame_t));
surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);
response_n = nan(ny,nx);
response_grid = cell(ny,nx);
for px_y = 1:ny
    for px_x = 1:nx
        
        % Get all square flips to white only
        align_stims = (stim_screen(px_y,px_x,2:end) - ...
            stim_screen(px_y,px_x,1:end-1)) > 0;
        
        % Get corresponding stim times (don't add 1: first screen blank
        % with no photodiode flip)
        align_times = stimOn_times(find(align_stims));
                
        response_n(px_y,px_x) = length(align_times);
        
        % Don't use times that fall outside of imaging
        align_times(align_times + surround_time(1) < frame_t(2) | ...
            align_times + surround_time(2) > frame_t(end)) = [];
        
        % Get stim-aligned responses, 2 choices:
        
        % 1) Interpolate times (slow - but supersamples so better)
        %         align_surround_times = bsxfun(@plus, align_times, surround_time);
        %         peri_stim_v = permute(mean(interp1(frame_t,fV',align_surround_times),1),[3,2,1]);
        
        % 2) Use closest frames to times (much faster - not different)
        align_surround_times = bsxfun(@plus, align_times, surround_time);
        frame_edges = [frame_t,frame_t(end)+1/framerate];
        align_frames = discretize(align_surround_times,frame_edges);
        
        align_frames(any(isnan(align_frames),2),:) = [];
        
        % Define the peri-stim V's as subtracting first frame (baseline)
        peri_stim_v = bsxfun(@minus, ...
            reshape(fV(:,align_frames)',size(align_frames,1),size(align_frames,2),[]), ...
            reshape(fV(:,align_frames(:,1))',size(align_frames(:,1),1),size(align_frames(:,1),2),[]));

        mean_peri_stim_v = permute(mean(peri_stim_v,2),[3,1,2]);
        
        % Save V's
        response_grid{px_y,px_x} = mean_peri_stim_v;
        
    end
end

% Get position preference for every pixel
U_downsample_factor = 1; %2 if max method
screen_resize_scale = 1; %3 if max method
filter_sigma = (screen_resize_scale*2);

% Downsample U
use_u_y = 1:Uy;
Ud = imresize(U(use_u_y,:,:),1/U_downsample_factor,'bilinear');

% Convert V responses to pixel responses
use_svs = 1:size(U,3);
n_boot = 10;

response_mean_boostrap = cellfun(@(x) bootstrp(n_boot,@mean,x')',response_grid,'uni',false);

% (to split trials instead of bootstrap)
%split_trials = cellfun(@(x) shake(discretize(1:size(x,2),round(linspace(1,size(x,2),n_boot+1)))),response_grid,'uni',false);
%response_mean_boostrap = cellfun(@(x,y) grpstats(x',y','mean')',response_grid,split_trials,'uni',false);
use_method = 'max'; % max or com
vfs_boot = nan(size(Ud,1),size(Ud,2),n_boot);
for curr_boot = 1:n_boot
    
    response_mean = cell2mat(cellfun(@(x) x(:,curr_boot),response_mean_boostrap(:),'uni',false)');
    stim_im_px = reshape(permute(svdFrameReconstruct(Ud(:,:,use_svs),response_mean(use_svs,:)),[3,1,2]),ny,nx,[]);
    gauss_filt = fspecial('gaussian',[ny,nx],filter_sigma);
    stim_im_smoothed = imfilter(imresize(stim_im_px,screen_resize_scale,'bilinear'),gauss_filt);

    switch use_method
        case 'max'
            % Upsample each pixel's response map and find maximum            
            [~,mi] = max(reshape(stim_im_smoothed,[],size(stim_im_px,3)),[],1);
            [m_y,m_x] = ind2sub(size(stim_im_smoothed),mi);
            m_yr = reshape(m_y,size(Ud,1),size(Ud,2));
            m_xr = reshape(m_x,size(Ud,1),size(Ud,2));
            
        case 'com'
            % Conversely, do COM on original^2
            [xx,yy] = meshgrid(1:size(stim_im_smoothed,2),1:size(stim_im_smoothed,1));
            m_xr = reshape(sum(sum(bsxfun(@times,stim_im_smoothed.^2,xx),1),2)./sum(sum(stim_im_smoothed.^2,1),2),size(Ud,1),size(Ud,2));
            m_yr = reshape(sum(sum(bsxfun(@times,stim_im_smoothed.^2,yy),1),2)./sum(sum(stim_im_smoothed.^2,1),2),size(Ud,1),size(Ud,2));
    end
    
    % Calculate and plot sign map (dot product between horz & vert gradient)
    
    % 1) get gradient direction
    [~,Vdir] = imgradient(imgaussfilt(m_yr,1));
    [~,Hdir] = imgradient(imgaussfilt(m_xr,1));
    
    % 3) get sin(difference in direction) if retinotopic, H/V should be
    % orthogonal, so the closer the orthogonal the better (and get sign)
    angle_diff = sind(Vdir-Hdir);
    angle_diff(isnan(angle_diff)) = 0;
    
    vfs_boot(:,:,curr_boot) = angle_diff;
    
    disp(curr_boot);
    
end

vfs_median = imgaussfilt(nanmedian(vfs_boot,3),2);

figure('Name',animal);
ax1 = axes;
subplot(1,2,1,ax1);
imagesc(vfs_median);
caxis([-1,1]);
axes(ax1); axis image off;
colormap(colormap_BlueWhiteRed)

ax2 = axes;
ax3 = axes;
subplot(1,2,2,ax2);
subplot(1,2,2,ax3);
h1 = imagesc(ax2,avg_im(use_u_y,:));
colormap(ax2,gray);
caxis(ax2,[0 prctile(avg_im(:),99.9)]);
h2 = imagesc(ax3,vfs_median);
colormap(ax3,colormap_BlueWhiteRed);
caxis([-1,1]);
set(ax2,'Visible','off');
axes(ax2); axis image off;
set(ax3,'Visible','off');
axes(ax3); axis image off;
set(h2,'AlphaData',mat2gray(abs(vfs_median))*0.3);
colormap(ax2,gray);

% vfs_cutoff = 1.5*std(vfs_mean(:));
% vis_thresh = abs(vfs_mean) > vfs_cutoff;
% vis_boundaries = cellfun(@(x) x*U_downsample_factor,bwboundaries(vis_thresh),'uni',false);
% 
% figure; hold on;
% set(gca,'YDir','reverse');
% imagesc(avg_im); colormap(gray);
% caxis([0 prctile(avg_im(:),95)]);
% for i = 1:length(vis_boundaries)
%     plot(vis_boundaries{i}(:,2),vis_boundaries{i}(:,1),'m','linewidth',1);
% end
% ylim([0,size(avg_im,1)]);
% xlim([0,size(avg_im,2)]);
% axis off;



%% BLACK BACKGROUND: Sparse noise retinotopy (bootstrap across trials to get less noise)

[Uy,Ux,nSV] = size(U);

myScreenInfo.windowPtr = NaN; % so we can call the stimulus generation and it won't try to display anything
stimNum = 1;
ss = eval([Protocol.xfile(1:end-2) '(myScreenInfo, Protocol.pars(:,stimNum));']);
stim_screen = cat(3,ss.ImageTextures{:});
ny = size(stim_screen,1);
nx = size(stim_screen,2);

switch lower(photodiode_type)
    case 'flicker'
        % Check for case of mismatch between photodiode and stimuli:
        % odd number of stimuli, but one extra photodiode flip to come back down
        if mod(size(stim_screen,3),2) == 1 && ...
                length(photodiode.timestamps) == size(stim_screen,3) + 1;
            photodiode.timestamps(end) = [];
            photodiode.values(end) = [];
            warning('Odd number of stimuli, removed last photodiode');
        end
        
        % If there's still a mismatch, break
        if size(stim_screen,3) ~= length(photodiode.timestamps);
            warning([num2str(size(stim_screen,3)) ' stimuli, ', ...
                num2str(length(photodiode.timestamps)) ' photodiode pulses']);
            
            % Try to estimate which stim were missed by time difference
            photodiode_diff = diff(photodiode.timestamps);
            max_regular_diff_time = prctile(diff(photodiode.timestamps),99);
            skip_cutoff = max_regular_diff_time*2;
            photodiode_skip = find(photodiode_diff > skip_cutoff);
            est_n_pulse_skip = ceil(photodiode_diff(photodiode_skip)/max_regular_diff_time)-1;
            stim_skip = cell2mat(arrayfun(@(x) photodiode_skip(x):photodiode_skip(x)+est_n_pulse_skip(x)-1, ...
                1:length(photodiode_skip),'uni',false));
            
            if isempty(est_n_pulse_skip) || length(photodiode.timestamps) + sum(est_n_pulse_skip) ~= size(stim_screen,3)
                error('Can''t match photodiode events to stimuli')
            end
        end
        
        stim_times = photodiode.timestamps;
        
    case 'steady'
        % If the photodiode is on steady: extrapolate the stim times
        if length(photodiode.timestamps) ~= 2
            error('Steady photodiode, but not 2 flips')
        end
        stim_duration = diff(photodiode.timestamps)/size(stim_screen,3);
        stim_times = linspace(photodiode.timestamps(1), ...
            photodiode.timestamps(2)-stim_duration,size(stim_screen,3))';
        
end

% Get average response to each stimulus
surround_window = [0.2,0.3];
framerate = 1./nanmedian(diff(frame_t));
surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);
response_n = nan(ny,nx);
response_grid = cell(ny,nx);
for px_y = 1:ny;
    for px_x = 1:nx;
        
        % Use first frame of light stim
        align_stims = (stim_screen(px_y,px_x,2:end) == 1) & ...
            (diff(stim_screen(px_y,px_x,:),[],3) == 1);
        align_times = stim_times(find(align_stims)+1);
        
        align_times = align_times(round(length(align_times)/2):end);
        
        response_n(px_y,px_x) = length(align_times);
        
        % Don't use times that fall outside of imaging
        align_times(align_times + surround_time(1) < frame_t(2) | ...
            align_times + surround_time(2) > frame_t(end)) = [];
        
        % Get stim-aligned responses, 2 choices:
        
        % 1) Interpolate times (slow - but supersamples so better)
        %         align_surround_times = bsxfun(@plus, align_times, surround_time);
        %         peri_stim_v = permute(mean(interp1(frame_t,fV',align_surround_times),1),[3,2,1]);
        
        % 2) Use closest frames to times (much faster - not different)
        align_surround_times = bsxfun(@plus, align_times, surround_time);
        frame_edges = [frame_t,frame_t(end)+1/framerate];
        align_frames = discretize(align_surround_times,frame_edges);
        
        align_frames(any(isnan(align_frames),2),:) = [];
        
        peri_stim_v = squeeze(nanmean(reshape(fV(:,align_frames)', ...
            size(align_frames,1),size(align_frames,2),[]),2))';
        
        % Save V's
        response_grid{px_y,px_x} = peri_stim_v;
        
    end
end

% Get position preference for every pixel
U_downsample_factor = 3;
resize_scale = 3;
filter_sigma = (resize_scale*1.5);

% Downsample U
use_u_y = 1:Uy;
Ud = imresize(U(use_u_y,:,:),1/U_downsample_factor,'bilinear');

% Convert V responses to pixel responses
use_svs = 1:size(U,3);
n_boot = 10;

response_mean_boostrap = cellfun(@(x) bootstrp(n_boot,@mean,x')',response_grid,'uni',false);

% (to split trials instead of bootstrap)
%split_trials = cellfun(@(x) shake(discretize(1:size(x,2),round(linspace(1,size(x,2),n_boot+1)))),response_grid,'uni',false);
%response_mean_boostrap = cellfun(@(x,y) grpstats(x',y','mean')',response_grid,split_trials,'uni',false);

vfs_boot = nan(size(Ud,1),size(Ud,2),n_boot);
for curr_boot = 1:n_boot
    
    response_mean = cell2mat(cellfun(@(x) x(:,curr_boot),response_mean_boostrap(:),'uni',false)');
    stim_im_px = reshape(permute(svdFrameReconstruct(Ud(:,:,use_svs),response_mean(use_svs,:)),[3,1,2]),ny,nx,[]);
    
    % Upsample each pixel's response map and find maximum
    gauss_filt = fspecial('gaussian',[ny,nx],filter_sigma);
    stim_im_smoothed = imfilter(imresize(stim_im_px,resize_scale,'bilinear'),gauss_filt);
    [mc,mi] = max(reshape(stim_im_smoothed,[],size(stim_im_px,3)),[],1);
    [m_y,m_x] = ind2sub(size(stim_im_smoothed),mi);
    m_yr = reshape(m_y,size(Ud,1),size(Ud,2));
    m_xr = reshape(m_x,size(Ud,1),size(Ud,2));
    
    % Calculate and plot sign map (do this just with dot product between horz / vert grad?)
    
    % 1) get gradient direction
    [Vmag,Vdir] = imgradient(imgaussfilt(m_yr,1));
    [Hmag,Hdir] = imgradient(imgaussfilt(m_xr,1));
    
    % 3) get sin(difference in direction) if retinotopic, H/V should be
    % orthogonal, so the closer the orthogonal the better (and get sign)
    angle_diff = sind(Vdir-Hdir);
    
    vfs_boot(:,:,curr_boot) = angle_diff;
    
    disp(curr_boot);
    
end

vfs_mean = nanmean(imgaussfilt(vfs_boot,1),3);

figure;
imagesc(vfs_mean);
caxis([-1,1]);
axis off;
title('Visual field sign')

figure;
subplot(1,3,1);
imagesc(vfs_mean);
caxis([-1,1]);
axis off;
title('All')

subplot(1,3,2);
p = reshape(AP_signrank_matrix(reshape(imgaussfilt(vfs_boot,1),[],n_boot)'), ...
    size(vfs_mean,1),size(vfs_mean,2));
h = imagesc(vfs_mean);
caxis([-1,1]);
axis off;
set(h,'AlphaData',p < 0.05);
title('Sign rank p < 0.05');

subplot(1,3,3);
h = imagesc(vfs_mean);
caxis([-1,1]);
axis off;
vfs_cutoff = 1.5*std(vfs_mean(:));
set(h,'AlphaData',abs(vfs_mean) > vfs_cutoff);
title('Abs > 1.5*std');

vis_thresh = abs(vfs_mean) > vfs_cutoff;
vis_boundaries = cellfun(@(x) x*U_downsample_factor,bwboundaries(vis_thresh),'uni',false);

figure; hold on;
set(gca,'YDir','reverse');
imagesc(avg_im); colormap(gray);
caxis([0 prctile(avg_im(:),95)]);
for i = 1:length(vis_boundaries)
    plot(vis_boundaries{i}(:,2),vis_boundaries{i}(:,1),'m','linewidth',1);
end
ylim([0,size(avg_im,1)]);
xlim([0,size(avg_im,2)]);
axis off;


%% Sparse noise retinotopy (Kenneth's suggestions for reducing non-vis)

% "Whiten" noise in V
% 1) filter V so all frequencies are equal (I forget why?)
% highpassCutoff = 0.2; % Hz
% [b100s, a100s] = butter(2, highpassCutoff/(Fs/2), 'high');
% fV_white = filter(b100s,a100s,fV,[],2);
fV_white = fV./conv2(fV,[1,-0.9],'same');

% 2) normalize values of V (in V/(std(V)+epsilon) where epsilon is the estimated noise)
epsilon = 1*std(fV(end,:));
fV_white = bsxfun(@rdivide,fV_white,std(fV_white,[],2)+epsilon);

[Uy,Ux,nSV] = size(U);

myScreenInfo.windowPtr = NaN; % so we can call the stimulus generation and it won't try to display anything
stimNum = 1;
ss = eval([Protocol.xfile(1:end-2) '(myScreenInfo, Protocol.pars(:,stimNum));']);
stim_screen = cat(3,ss.ImageTextures{:});
ny = size(stim_screen,1);
nx = size(stim_screen,2);

% Check number of photodiode pulses vs. number of stimuli
% (this is really not a good way to do this, why does it skip photodiode
% events??)
if size(stim_screen,3) ~= length(photodiode.timestamps);
    warning([num2str(size(stim_screen,3)) ' stimuli, ', ...
        num2str(length(photodiode.timestamps)) ' photodiode pulses']);
    
    % Try to estimate which stim were missed by time difference
    photodiode_diff = diff(photodiode.timestamps);
    max_regular_diff_time = prctile(diff(photodiode.timestamps),99);
    skip_cutoff = max_regular_diff_time*2;
    photodiode_skip = find(photodiode_diff > skip_cutoff);
    est_n_pulse_skip = ceil(photodiode_diff(photodiode_skip)/max_regular_diff_time)-1;
    stim_skip = cell2mat(arrayfun(@(x) photodiode_skip(x):photodiode_skip(x)+est_n_pulse_skip(x)-1, ...
        1:length(photodiode_skip),'uni',false));
    
    if isempty(est_n_pulse_skip) || length(photodiode.timestamps) + sum(est_n_pulse_skip) ~= size(stim_screen,3)
        error('Can''t match photodiode events to stimuli')
    else
        warning('Estimated skipped stimuli, removing those and continuing');
        stim_screen(:,:,stim_skip) = [];
    end
end

% Get average response to each stimulus
surround_window = [0.2,0.3];
framerate = 1./nanmedian(diff(frame_t));
surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);
response_n = nan(ny,nx);
response_grid = cell(ny,nx);
for px_y = 1:ny;
    for px_x = 1:nx;
        
        % Use first frame of dark or light stim
        align_stims = (stim_screen(px_y,px_x,2:end)~= 0) & ...
            (diff(stim_screen(px_y,px_x,:),[],3) ~= 0);
        align_times = photodiode.timestamps(find(align_stims)+1);
        
        align_times = align_times(round(length(align_times)/2):end);
        
        response_n(px_y,px_x) = length(align_times);
        
        % Don't use times that fall outside of imaging
        align_times(align_times + surround_time(1) < frame_t(2) | ...
            align_times + surround_time(2) > frame_t(end)) = [];
        
        % Get stim-aligned responses, 2 choices:
        
        % 1) Interpolate times (slow - but supersamples so better)
        %         align_surround_times = bsxfun(@plus, align_times, surround_time);
        %         peri_stim_v = permute(mean(interp1(frame_t,fV',align_surround_times),1),[3,2,1]);
        
        % 2) Use closest frames to times (much faster - not different)
        align_surround_times = bsxfun(@plus, align_times, surround_time);
        frame_edges = [frame_t,frame_t(end)+1/framerate];
        align_frames = discretize(align_surround_times,frame_edges);
        
        align_frames(any(isnan(align_frames),2),:) = [];
        
        peri_stim_v = permute(nanmean(reshape(fV_white(:,align_frames)', ...
            size(align_frames,1),size(align_frames,2),[]),1),[3,2,1]);
        
        % Save V's
        response_grid{px_y,px_x} = peri_stim_v;
        
    end
end

% Get position preference for every pixel
U_downsample_factor = 3;
resize_scale = 3;
filter_sigma = (resize_scale*1.5);

% Downsample U
use_u_y = 1:Uy;
Ud = imresize(U(use_u_y,:,:),1/U_downsample_factor,'bilinear');

% Convert V responses to pixel responses
use_svs = 1:size(U,3);
response_mean = cell2mat(cellfun(@(x) nanmean(x,2),response_grid(:),'uni',false)');
stim_im_px = reshape(permute(svdFrameReconstruct(Ud(:,:,use_svs),response_mean(use_svs,:)),[3,1,2]),ny,nx,[]);

% Upsample each pixel's response map and find maximum
gauss_filt = fspecial('gaussian',[ny,nx],filter_sigma);
stim_im_smoothed = imfilter(imresize(stim_im_px,resize_scale,'bilinear'),gauss_filt);
[mc,mi] = max(reshape(stim_im_smoothed,[],size(stim_im_px,3)),[],1);
[m_y,m_x] = ind2sub(size(stim_im_smoothed),mi);
m_yr = reshape(m_y,size(Ud,1),size(Ud,2));
m_xr = reshape(m_x,size(Ud,1),size(Ud,2));

% Plot retinotopy in different ways:
figure;

% Plot retinotopy by color
screen_pos_col = nan(ny,nx,3);
[screen_pos_col(:,:,1),screen_pos_col(:,:,2)] = ...
    meshgrid(mat2gray(1:nx),mat2gray(1:ny));
screen_pos_col(:,:,3) = 1- ...
    (screen_pos_col(:,:,1)+screen_pos_col(:,:,2))/2;
screen_pos_col_upsample = reshape(imresize(screen_pos_col,resize_scale),[],3);

retinotopy_colormap = reshape(screen_pos_col_upsample(mi,:),size(Ud,1),size(Ud,2),3);

subplot(3,1,1);
h = imagesc(retinotopy_colormap);
set(gca,'YDir','Normal');
axis off;
title('Positional retinotopy');

% Plot isovertical/horizontal lines
stripe_spacing = 5;
vert_stripes = imgaussfilt(double(mod(m_xr-1,stripe_spacing) == 0),1);
horz_stripes = imgaussfilt(double(mod(m_yr-1,stripe_spacing) == 0),1);
subplot(3,1,2);
imagesc(padarray(cat(3,mat2gray(vert_stripes),mat2gray(horz_stripes)),[0,0,1],0,'post'));
set(gca,'YDir','normal')
axis off;
title('Iso-horizontal/vertical retinotopy');

% (get all vertical/horizontal stripes for plotting manually)
vert_response = nan(size(Ud,1),size(Ud,2),size(stim_im_smoothed,2));
for i = 1:size(stim_im_smoothed,2)
    vert_response(:,:,i) = m_xr == i;
end
vert_response = imgaussfilt(vert_response,1);

horz_response = nan(size(Ud,1),size(Ud,2),size(stim_im_smoothed,1));
for i = 1:size(stim_im_smoothed,1)
    horz_response(:,:,i) = m_yr == i;
end
horz_response = imgaussfilt(horz_response,1);

% Calculate and plot sign map (do this just with dot product between horz / vert grad?)

% 1) get gradient direction
[Vmag,Vdir] = imgradient(imgaussfilt(m_yr,1));
[Hmag,Hdir] = imgradient(imgaussfilt(m_xr,1));

% 3) get sin(difference in direction) if retinotopic, H/V should be
% orthogonal, so the closer the orthogonal the better (and get sign)
angle_diff = sind(Vdir-Hdir);

subplot(3,1,3);
imagesc(imgaussfilt(angle_diff,1));
set(gca,'YDir','normal')
axis off;
title('Visual sign field');


%% Sparse noise retinotopy by regression

[Uy,Ux,nSV] = size(U);

myScreenInfo.windowPtr = NaN; % so we can call the stimulus generation and it won't try to display anything
stimNum = 1;
ss = eval([Protocol.xfile(1:end-2) '(myScreenInfo, Protocol.pars(:,stimNum));']);
stim_screen = cat(3,ss.ImageTextures{:});
ny = size(stim_screen,1);
nx = size(stim_screen,2);


switch photodiode_type
    case 'flicker'
        % Check for case of mismatch between photodiode and stimuli:
        % odd number of stimuli, but one extra photodiode flip to come back down
        if mod(size(stim_screen,3),2) == 1 && ...
                length(photodiode.timestamps) == size(stim_screen,3) + 1;
            photodiode.timestamps(end) = [];
            photodiode.values(end) = [];
            warning('Odd number of stimuli, removed last photodiode');
        end
        
        % If there's still a mismatch, break
        if size(stim_screen,3) ~= length(photodiode.timestamps);
            warning([num2str(size(stim_screen,3)) ' stimuli, ', ...
                num2str(length(photodiode.timestamps)) ' photodiode pulses']);
            
            % Try to estimate which stim were missed by time difference
            photodiode_diff = diff(photodiode.timestamps);
            max_regular_diff_time = prctile(diff(photodiode.timestamps),99);
            skip_cutoff = max_regular_diff_time*2;
            photodiode_skip = find(photodiode_diff > skip_cutoff);
            est_n_pulse_skip = ceil(photodiode_diff(photodiode_skip)/max_regular_diff_time)-1;
            stim_skip = cell2mat(arrayfun(@(x) photodiode_skip(x):photodiode_skip(x)+est_n_pulse_skip(x)-1, ...
                1:length(photodiode_skip),'uni',false));
            
            if isempty(est_n_pulse_skip) || length(photodiode.timestamps) + sum(est_n_pulse_skip) ~= size(stim_screen,3)
                error('Can''t match photodiode events to stimuli')
            end
        end
        
        stim_times = photodiode.timestamps;
        
    case 'Steady'
        % If the photodiode is on steady: extrapolate the stim times
        if length(photodiode.timestamps) ~= 2
            error('Steady photodiode, but not 2 flips')
        end
        stim_duration = diff(photodiode.timestamps)/size(stim_screen,3);
        stim_times = linspace(photodiode.timestamps(1), ...
            photodiode.timestamps(2)-stim_duration,size(stim_screen,3))';
        
end

% Interpolate stim screen for all times
stim_screen_interp = abs(single(interp1(stim_times,reshape(stim_screen,[],length(stim_times))',frame_t,'nearest')'));
stim_screen_interp(isnan(stim_screen_interp)) = 0;

% Use only the onsets of stim
stim_screen_interp = single([zeros(size(stim_screen_interp,1),1),diff(stim_screen_interp,[],2)] == 1);

% Skip the first n seconds to do this
skip_seconds = 10;
use_frames = (frame_t > skip_seconds);

use_svs = 100:500;
fluor_kernel_frames = -17:-6;
lambda = 1e6;
zs = [false,true];

k = AP_regresskernel(fV(use_svs,use_frames), ...
    stim_screen_interp(:,use_frames),fluor_kernel_frames,lambda,zs,5);

% Get maximum kernel response for all pixels
r = permute(mean(reshape(k,length(use_svs), ...
    length(fluor_kernel_frames),size(stim_screen_interp,1)),2),[1,3,2]);

% Get position preference for every pixel
U_downsample_factor = 1;
resize_scale = 3;
filter_sigma = (resize_scale*1.5);

% Downsample U
use_u_y = 1:Uy;
Ud = imresize(U(use_u_y,:,:),1/U_downsample_factor,'bilinear');

% Convert V responses to pixel responses
stim_im_px = reshape(permute(svdFrameReconstruct(Ud(:,:,use_svs),r),[3,1,2]),ny,nx,[]);

% Upsample each pixel's response map and find maximum
gauss_filt = fspecial('gaussian',[ny,nx],filter_sigma);
stim_im_smoothed = imfilter(imresize(stim_im_px,resize_scale,'bilinear'),gauss_filt);
[mc,mi] = max(reshape(stim_im_smoothed,[],size(stim_im_px,3)),[],1);
[m_y,m_x] = ind2sub(size(stim_im_smoothed),mi);
m_yr = reshape(m_y,size(Ud,1),size(Ud,2));
m_xr = reshape(m_x,size(Ud,1),size(Ud,2));

% 1) get gradient direction
[Vmag,Vdir] = imgradient(imgaussfilt(m_yr,1));
[Hmag,Hdir] = imgradient(imgaussfilt(m_xr,1));

% 2) get sin(difference in direction) if retinotopic, H/V should be
% orthogonal, so the closer the orthogonal the better (and get sign)
angle_diff = sind(Vdir-Hdir);

vfs_mean = nanmean(imgaussfilt(angle_diff,1),3);

figure;
imagesc(vfs_mean);
caxis([-1,1]);
axis off;
title('All')

vfs_cutoff = 1.5*std(vfs_mean(:));
vis_thresh = abs(vfs_mean) > vfs_cutoff;
vis_boundaries = cellfun(@(x) x*U_downsample_factor,bwboundaries(vis_thresh),'uni',false);

figure; hold on;
set(gca,'YDir','reverse');
imagesc(avg_im); colormap(gray);
caxis([0 prctile(avg_im(:),95)]);
for i = 1:length(vis_boundaries)
    plot(vis_boundaries{i}(:,2),vis_boundaries{i}(:,1),'m','linewidth',1);
end
ylim([0,size(avg_im,1)]);
xlim([0,size(avg_im,2)]);
axis off;


%% Sparse noise retinotopy with Nick's code

[Uy,Ux,nSV] = size(U);

myScreenInfo.windowPtr = NaN; % so we can call the stimulus generation and it won't try to display anything
stimNum = 1;
ss = eval([Protocol.xfile(1:end-2) '(myScreenInfo, Protocol.pars(:,stimNum));']);
stim_screen = cat(3,ss.ImageTextures{:});
ny = size(stim_screen,1);
nx = size(stim_screen,2);

switch lower(photodiode_type)
    case 'flicker'
        % Check for case of mismatch between photodiode and stimuli:
        % odd number of stimuli, but one extra photodiode flip to come back down
        if mod(size(stim_screen,3),2) == 1 && ...
                length(photodiode.timestamps) == size(stim_screen,3) + 1;
            photodiode.timestamps(end) = [];
            photodiode.values(end) = [];
            warning('Odd number of stimuli, removed last photodiode');
        end
        
        % If there's still a mismatch, break
        if size(stim_screen,3) ~= length(photodiode.timestamps);
            warning([num2str(size(stim_screen,3)) ' stimuli, ', ...
                num2str(length(photodiode.timestamps)) ' photodiode pulses']);
            
            % Try to estimate which stim were missed by time difference
            photodiode_diff = diff(photodiode.timestamps);
            max_regular_diff_time = prctile(diff(photodiode.timestamps),99);
            skip_cutoff = max_regular_diff_time*2;
            photodiode_skip = find(photodiode_diff > skip_cutoff);
            est_n_pulse_skip = ceil(photodiode_diff(photodiode_skip)/max_regular_diff_time)-1;
            stim_skip = cell2mat(arrayfun(@(x) photodiode_skip(x):photodiode_skip(x)+est_n_pulse_skip(x)-1, ...
                1:length(photodiode_skip),'uni',false));
            
            if isempty(est_n_pulse_skip) || length(photodiode.timestamps) + sum(est_n_pulse_skip) ~= size(stim_screen,3)
                error('Can''t match photodiode events to stimuli')
            end
        end
        
        stim_times = photodiode.timestamps;
        
    case 'steady'
        % If the photodiode is on steady: extrapolate the stim times
        if length(photodiode.timestamps) ~= 2
            error('Steady photodiode, but not 2 flips')
        end
        stim_duration = diff(photodiode.timestamps)/size(stim_screen,3);
        stim_times = linspace(photodiode.timestamps(1), ...
            photodiode.timestamps(2)-stim_duration,size(stim_screen,3))';
        
end

% Get times and positions of each stimulus
sparse_noise_position_times = cell(ny,nx);
for px_y = 1:ny;
    for px_x = 1:nx;
        
        % Use first frame of dark or light stim
        align_stims = (stim_screen(px_y,px_x,2:end)~= 0) & ...
            (diff(stim_screen(px_y,px_x,:),[],3) ~= 0);
        align_times = stim_times(find(align_stims)+1);
        
        align_times = align_times(round(length(align_times)/2):end);
        sparse_noise_position_times{px_y,px_x} = align_times;
        
    end
end
[x,y] = meshgrid(1:nx,1:ny);
sparse_noise_positions = cellfun(@(times,px_x,px_y) repmat([px_x,px_y], ...
    length(times),1),sparse_noise_position_times,num2cell(x),num2cell(y),'uni',false);

[signMap, xMap, yMap] = sparseRetinotopy(Udf,fVdf,frame_t, ...
    vertcat(sparse_noise_positions{:}),vertcat(sparse_noise_position_times{:}),avg_im);

signMap_blur = imgaussfilt(signMap,2);
vfs_cutoff = 0.8*nanstd(signMap_blur(:));
vis_thresh = abs(signMap_blur) > vfs_cutoff;
vis_boundaries = cellfun(@(x) x*U_downsample_factor,bwboundaries(vis_thresh),'uni',false);

figure; hold on;
set(gca,'YDir','reverse');
imagesc(avg_im); colormap(gray);
caxis([0 prctile(avg_im(:),95)]);
for i = 1:length(vis_boundaries)
    plot(vis_boundaries{i}(:,2),vis_boundaries{i}(:,1),'m','linewidth',1);
end
ylim([0,size(avg_im,1)]);
xlim([0,size(avg_im,2)]);
axis off;

%% Sparse noise with Nick's code in batch

animal = 'AP025';
protocol = 'stimSparseNoiseUncorrAsync';
experiments = AP_find_experiments(animal,protocol);

% only use experiments with ephys + imaging
experiments = experiments([experiments.imaging] & [experiments.ephys]);

load_parts.cam = false;
load_parts.imaging = true;
load_parts.ephys = false;

batch_vars = struct;
for curr_day = 1:length(experiments);
    
    day = experiments(curr_day).day;
    experiment = experiments(curr_day).experiment;
    
    AP_load_experiment
    
    %%%%%%%%%%%%%%%
    % DO THE STUFF
    %%%%%%%%%%%%%%%
    
    [Uy,Ux,nSV] = size(U);
    
    myScreenInfo.windowPtr = NaN; % so we can call the stimulus generation and it won't try to display anything
    stimNum = 1;
    ss = eval([Protocol.xfile(1:end-2) '(myScreenInfo, Protocol.pars(:,stimNum));']);
    stim_screen = cat(3,ss.ImageTextures{:});
    ny = size(stim_screen,1);
    nx = size(stim_screen,2);
    
    switch lower(photodiode_type)
        case 'flicker'
            % Check for case of mismatch between photodiode and stimuli:
            % odd number of stimuli, but one extra photodiode flip to come back down
            if mod(size(stim_screen,3),2) == 1 && ...
                    length(photodiode.timestamps) == size(stim_screen,3) + 1;
                photodiode.timestamps(end) = [];
                photodiode.values(end) = [];
                warning('Odd number of stimuli, removed last photodiode');
            end
            
            % If there's still a mismatch, break
            if size(stim_screen,3) ~= length(photodiode.timestamps);
                warning([num2str(size(stim_screen,3)) ' stimuli, ', ...
                    num2str(length(photodiode.timestamps)) ' photodiode pulses']);
                
                % Try to estimate which stim were missed by time difference
                photodiode_diff = diff(photodiode.timestamps);
                max_regular_diff_time = prctile(diff(photodiode.timestamps),99);
                skip_cutoff = max_regular_diff_time*2;
                photodiode_skip = find(photodiode_diff > skip_cutoff);
                est_n_pulse_skip = ceil(photodiode_diff(photodiode_skip)/max_regular_diff_time)-1;
                stim_skip = cell2mat(arrayfun(@(x) photodiode_skip(x):photodiode_skip(x)+est_n_pulse_skip(x)-1, ...
                    1:length(photodiode_skip),'uni',false));
                
                if isempty(est_n_pulse_skip) || length(photodiode.timestamps) + sum(est_n_pulse_skip) ~= size(stim_screen,3)
                    error('Can''t match photodiode events to stimuli')
                end
            end
            
            stim_times = photodiode.timestamps;
            
        case 'steady'
            % If the photodiode is on steady: extrapolate the stim times
            if length(photodiode.timestamps) ~= 2
                error('Steady photodiode, but not 2 flips')
            end
            stim_duration = diff(photodiode.timestamps)/size(stim_screen,3);
            stim_times = linspace(photodiode.timestamps(1), ...
                photodiode.timestamps(2)-stim_duration,size(stim_screen,3))';
            
    end
    
    % Get times and positions of each stimulus
    sparse_noise_position_times = cell(ny,nx);
    for px_y = 1:ny;
        for px_x = 1:nx;
            
            % Use first frame of dark or light stim
            align_stims = (stim_screen(px_y,px_x,2:end)~= 0) & ...
                (diff(stim_screen(px_y,px_x,:),[],3) ~= 0);
            align_times = stim_times(find(align_stims)+1);
            
            align_times = align_times(round(length(align_times)/2):end);
            sparse_noise_position_times{px_y,px_x} = align_times;
            
        end
    end
    [x,y] = meshgrid(1:nx,1:ny);
    sparse_noise_positions = cellfun(@(times,px_x,px_y) repmat([px_x,px_y], ...
        length(times),1),sparse_noise_position_times,num2cell(x),num2cell(y),'uni',false);
    
    [signMap, xMap, yMap] = sparseRetinotopy(Udf,fVdf,frame_t, ...
        vertcat(sparse_noise_positions{:}),vertcat(sparse_noise_position_times{:}),avg_im);
    
    batch_vars.signMap{curr_day} = signMap;
    
    %%%%%%%%%%%%%%%%%%%%
    % THE STUFF IS DONE
    %%%%%%%%%%%%%%%%%%%%
    
    drawnow
    disp(curr_day)
    clearvars -except experiments curr_day animal batch_vars load_parts
    
end

disp('Finished batch.')

% Align
days = {experiments.day};
tform_matrix = AP_align_widefield(animal,days);

batch_vars_reg = batch_vars;

% Align across days and replace nans
for curr_day = setdiff(1:length(experiments),ref_im_num);
    
    tform = affine2d;
    tform.T = tform_matrix{curr_day};
    
    curr_im = batch_vars_reg.signMap{curr_day};
    nan_cond = squeeze(all(all(all(isnan(curr_im),1),2),3));
    curr_im(isnan(curr_im)) = 0;
    curr_im = imwarp(curr_im,tform, ...
        'Outputview',imref2d(size(avg_im{ref_im_num})));
    curr_im(:,:,:,nan_cond) = NaN;
    batch_vars_reg.signMap{curr_day} = curr_im;

end

signMap_cat = cat(3,batch_vars_reg.signMap{:});

signMap_mean = nanmean(signMap_cat,3);
avg_im_mean = nanmean(avg_im_reg,3);

% Get boundaries
% vfs_cutoff = 1.5*std(signMap_mean(:));
% vis_thresh = abs(signMap_mean) > vfs_cutoff;
% vis_boundaries = bwboundaries(vis_thresh);
% 
% figure; hold on;
% set(gca,'YDir','reverse');
% imagesc(avg_im_mean); colormap(gray);
% caxis([0 prctile(avg_im_mean(:),99)]);
% for i = 1:length(vis_boundaries)
%     plot(vis_boundaries{i}(:,2),vis_boundaries{i}(:,1),'m','linewidth',1);
% end
% axis off;

% Plot
figure;
ax1 = axes;
p1 = imagesc(avg_im_mean);
colormap(ax1,gray);
caxis([0 prctile(avg_im_mean(:),99.5)]);
axis image off

ax2 = axes;
p2 = imagesc(signMap_mean);
colormap(ax2,colormap_BlueWhiteRed);
set(p2,'AlphaData',mat2gray(abs(signMap_mean))*0.3);
axis image off

%% Sparse noise with my code in batch

animals = {'AP055'};

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    disp(animal);
    
%     protocol = 'stimSparseNoiseUncorrAsync';
    protocol = 'AP_sparseNoise';
    experiments = AP_find_experiments(animal,protocol);
%     experiments = experiments([experiments.imaging]);
    experiments = experiments([experiments.imaging] & ~[experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = false;
    
    vfs = cell(length(experiments),1);
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        AP_load_experiment;
        
        lilrig_retinotopy;
        vfs{curr_day} = vfs_median;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        clearvars -except animal experiments load_parts curr_day  vfs animals curr_animal
        
    end
    
    disp('Finished batch.')
    
    % Save retinotopy from all days
    retinotopy = struct('day',reshape({experiments.day},[],1),'vfs',vfs);
    
    alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment';
    retinotopy_path = [alignment_path filesep 'retinotopy'];
    retinotopy_fn = [retinotopy_path filesep animal '_retinotopy.mat'];
    save(retinotopy_fn,'retinotopy');
    
    clearvars -except animals curr_animal
    
end










