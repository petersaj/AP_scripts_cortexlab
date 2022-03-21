function cmap = AP_colormap(cmap_type,n_colors)
% cmap = AP_colormap(cmap_type,n_colors)
%
% cmap_type: 2-3 colors to transition between
% (currently defined: W,K,R,B,G,P)

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
col.R = [0.7,0,0]; % red
col.B = [0,0,0.7]; % blue
col.G = [0,0.7,0]; % green
col.P = [0.7,0,0.7]; % purple

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

    cmap = col.(cmap_type(1)) + ...
        linspace(0,1,n_colors)'.* ...
        ((col.(cmap_type(2)) - col.(cmap_type(1))));

elseif length(cmap_type) == 3
    % 3 colors = diverging

    cmap_bot = col.(cmap_type(1)) + ...
        linspace(0,1,floor(n_colors/2))'.* ...
        ((col.(cmap_type(2)) - col.(cmap_type(1))));

    cmap_top = col.(cmap_type(2)) + ...
        linspace(0,1,floor(n_colors/2))'.* ...
        ((col.(cmap_type(3)) - col.(cmap_type(2))));
    
    if mod(n_colors,2)
        % (if odd number of colors insert middle color at center)
        cmap = [cmap_bot;col,cmap_type(2);cmap_top];
    else
        % (if even number of colors, just concatenate)
        cmap = [cmap_bot;cmap_top];
    end

end

% If no output argument specified, apply colormap to current axes
if nargout == 0
    colormap(gca,cmap);
end


