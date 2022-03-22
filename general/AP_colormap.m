function cmap = AP_colormap(cmap_type,n_colors)
% cmap = AP_colormap(cmap_type,n_colors)
%
% cmap_type: 2-3 colors to transition between
% (currently defined: W,K,R,B,G,P)
%
% Interpolate in LAB colorspace, as in Brewermap and from blog post: 
% https://blogs.mathworks.com/steve/2006/05/09/a-lab-based-uniform-color-scale/


% Force cmap_type upper-case
cmap_type = upper(cmap_type);

% Set default number of colors
if ~exist('n_colors','var') || isempty(n_colors)
    n_colors = 2^8; % default 8-bit color
end

% Define single-letter colors (mostly approx from brewermap)
col = struct;
col.W = [1,1,1]; % white
col.K = [0,0,0]; % black
col.R = [0.5,0,0]; % red
col.B = [0,0,0.5]; % blue
col.G = [0,0.5,0]; % green
col.P = [0.5,0,0.5]; % purple

%%%%% WORKING HERE: GET ISOLUMINANT BASE COLORS?
% (colors, equal luminance, defined as angle in L*a*b* space)
% (olors picked from: https://rufflewind.com/_urandom/colorpicker/#00464b)
lum = 40;
col.R = lab2rgb([lum,58,27]);
col.G = lab2rgb([lum,-41,35]);
col.P = lab2rgb([lum,51,-40]);
col.B = lab2rgb([lum,0,30 ]);

% Check colors are valid RGB values
if(any(cell2mat(struct2cell(col)) < 0))
    error('Not all colors in RGB gamut')
end

% % (to view the colors)
% figure;
% imagesc(permute(cell2mat(struct2cell(col)),[1,3,2]));
% set(gca,'XTick','','YTickLabel',fieldnames(col));

% Make sure chosen colors are defined
defined_colors = ismember(cmap_type,cell2mat(fieldnames(col)));
if ~all(defined_colors)
    error('Color(s) %s not defined',cmap_type(~defined_colors));
end

if length(cmap_type) == 2
    % 2 colors = gradient

    cmap = lab2rgb(interp1([1,n_colors], ...
        rgb2lab(cat(1,col.(cmap_type(1)),col.(cmap_type(2)))),1:n_colors));

elseif length(cmap_type) == 3
    % 3 colors = diverging

    cmap_bot = lab2rgb(interp1([1,n_colors], ...
        rgb2lab(cat(1,col.(cmap_type(1)),col.(cmap_type(2)))),1:n_colors));

    cmap_top = lab2rgb(interp1([1,n_colors], ...
        rgb2lab(cat(1,col.(cmap_type(2)),col.(cmap_type(3)))),1:n_colors));
    
    if mod(n_colors,2)
        % (if odd number of colors insert middle color at center)
        cmap = [cmap_bot;col,cmap_type(2);cmap_top];
    else
        % (if even number of colors, just concatenate)
        cmap = [cmap_bot;cmap_top];
    end

end

% Apply colormap to current axes
if nargout == 0
    % If no output argument specified, apply to current axes
    colormap(gca,cmap);
    clear cmap
else
    % If output specified, return map without applying to anything
    return
end


