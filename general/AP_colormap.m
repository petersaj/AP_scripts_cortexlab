function cmap = AP_colormap(cmap_type,n_colors,cmap_gamma)
% cmap = AP_colormap(cmap_type,n_colors,cmap_gamma)
%
% cmap_type: gradient/diverging between black/white and color
% (currently defined: W,K,G,B,P,R - base must be white/black)%
% e.g. 'WR' = white-to-red, = 'BKR' blue-black-red
%
% n_colors: number of colors along the map (default = 2^8 = 256)
%
% cmap_gamma: value to gamma correct (resamples colormap values according to
% x^gamma, which emphasizes higher/lower values to give better contrast in
% different places) (default = 1: no gamma change)
%
% Interpolate in CIELAB colorspace to ensure linear luminance change, see:
% https://blogs.mathworks.com/steve/2006/05/09/a-lab-based-uniform-color-scale/

% Force cmap_type upper-case
cmap_type = upper(cmap_type);

% Set default number of colors
if ~exist('n_colors','var') || isempty(n_colors)
    n_colors = 2^8; % default 8-bit color
end

% If diverging colormap, force number of colors to be odd
if length(cmap_type) == 3
    n_colors = n_colors - (1-mod(n_colors,2));
end

% Set default gamma
if ~exist('cmap_gamma','var') || isempty(cmap_gamma)
    cmap_gamma = 1; % default: no gamma change
end

% If gamma, supersample colors to resample finer gradients
% (note: randomly 100-fold now, should base on gamma...)
% (force odd number of colors if diverging colormap)
if cmap_gamma == 1
    n_colors_generate = n_colors;
elseif    cmap_gamma ~= 1
    n_colors_generate = n_colors*100 + (length(cmap_type) == 3);
end

% Set luminance/chroma limits for colors
min_lum = 20;
max_lum = 50;

min_chroma = 10;
max_chroma = 40;

% Set black/white luminance
K_lum = 1;
W_lum = 99;

% Set colors in CIELAB space (angle, min chroma, max chroma)
col = struct;
col.W = [0,0,0];
col.K = [0,0,0];
col.G = [139,min_chroma,max_chroma];
col.P = [315,min_chroma,max_chroma];
col.B = [250,min_chroma,max_chroma];
col.R = [30,min_chroma,max_chroma];

% (to plot colors while writing/troubleshooting)
plot_col = false;
if plot_col
    all_col = cell2mat(cellfun(@(x) permute(lab2rgb(interp1([min_lum,max_lum], ...
        [min_lum,x(2)*cosd(x(1)),x(2)*sind(x(1)); ...
        max_lum,x(3)*cosd(x(1)),x(3)*sind(x(1))], ...
        linspace(min_lum,max_lum,256))),[1,3,2]),struct2cell(col),'uni',false)');
    figure;imagesc(all_col);
    set(gca,'XTickLabels',fieldnames(col),'YTick','');
end

% Check base color is black or white (others not supported now)
if ~ismember(cmap_type(ceil(length(cmap_type)/2)),{'W','K'})
    error('Base color must be white or black')
end

% Set luminance/chroma interpolation direction from base color
switch cmap_type(ceil(length(cmap_type)/2))
    case 'W'
        lum_interp = [W_lum,max_lum,min_lum];
        chroma_interp = [3,2];
    case 'K'
        lum_interp = [K_lum,min_lum,max_lum];
        chroma_interp = [2,3];
end

% Define gradients by interpolating in CIELAB space and converting to RGB
if length(cmap_type) == 2
    % 2 colors = gradient

    cmap_full = lab2rgb(interp1(lum_interp, ...
        [lum_interp(1),0,0;
        lum_interp(2), ...
        col.(cmap_type(2))(chroma_interp(1))*cosd(col.(cmap_type(2))(1)), ...
        col.(cmap_type(2))(chroma_interp(1))*sind(col.(cmap_type(2))(1)); ...
        lum_interp(3), ...
        col.(cmap_type(2))(chroma_interp(2))*cosd(col.(cmap_type(2))(1)), ...
        col.(cmap_type(2))(chroma_interp(2))*sind(col.(cmap_type(2))(1))], ...
        linspace(lum_interp(1),lum_interp(end),ceil(n_colors_generate/2)),'spline'));

elseif length(cmap_type) == 3
    % 3 colors = diverging 

    cmap_bot = lab2rgb(interp1(lum_interp, ...
        [lum_interp(1),0,0;
        lum_interp(2), ...
        col.(cmap_type(1))(chroma_interp(1))*cosd(col.(cmap_type(1))(1)), ...
        col.(cmap_type(1))(chroma_interp(1))*sind(col.(cmap_type(1))(1)); ...
        lum_interp(3), ...
        col.(cmap_type(1))(chroma_interp(2))*cosd(col.(cmap_type(1))(1)), ...
        col.(cmap_type(1))(chroma_interp(2))*sind(col.(cmap_type(1))(1))], ...
        linspace(lum_interp(1),lum_interp(end),ceil(n_colors_generate/2)),'spline'));

    cmap_top = lab2rgb(interp1(lum_interp, ...
        [lum_interp(1),0,0;
        lum_interp(2), ...
        col.(cmap_type(3))(chroma_interp(1))*cosd(col.(cmap_type(3))(1)), ...
        col.(cmap_type(3))(chroma_interp(1))*sind(col.(cmap_type(3))(1)); ...
        lum_interp(3), ...
        col.(cmap_type(3))(chroma_interp(2))*cosd(col.(cmap_type(3))(1)), ...
        col.(cmap_type(3))(chroma_interp(2))*sind(col.(cmap_type(3))(1))], ...
        linspace(lum_interp(1),lum_interp(end),ceil(n_colors_generate/2)),'spline'));
    
    % Combine bottom/top (center color is replicated - remove one)
    cmap_full = vertcat(flipud(cmap_bot),cmap_top(2:end,:));

end

% Force colormap into RGB gamut (probably not valid...)
cmap_full(cmap_full > 1) = 1;
cmap_full(cmap_full < 0) = 0;

% If gamma, resample colormap (created as linear x, sampled at x^gamma)
if cmap_gamma == 1
    cmap = cmap_full;
elseif    cmap_gamma ~= 1
    if length(cmap_type) == 2
        cmap_x = linspace(0,1,n_colors);
        cmap_x_gamma = sign(cmap_x).*abs(cmap_x).^cmap_gamma;
        cmap_gamma_idx = 1+floor((n_colors_generate-1)* ...
            (cmap_x_gamma-min(cmap_x_gamma))./range(cmap_x_gamma));
        cmap = cmap_full(cmap_gamma_idx,:);
    elseif length(cmap_type) == 3
        cmap_x = linspace(-1,1,n_colors);
        cmap_x_gamma = sign(cmap_x).*abs(cmap_x).^cmap_gamma;
        cmap_gamma_idx = 1+floor((n_colors_generate-1)* ...
            (cmap_x_gamma-min(cmap_x_gamma))./range(cmap_x_gamma));
        cmap = cmap_full(cmap_gamma_idx,:);
    end
end












