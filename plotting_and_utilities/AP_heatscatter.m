function AP_heatscatter(x,y,n_bins,point_thresh)
% AP_heatscatter(x,y,n_bins,point_thresh)
%
% Scatter plot turning into heatmap when dense
% x = x data
% y = y data
% n_bins = number of histogram bins for x and y (1 or 2 elements)
% point_thresh = # points threshold for scatter vs heatmap

% Force column orientation
x = reshape(x,[],1);
y = reshape(y,[],1);

% Cut out NaNs
nan_points = isnan(x) | isnan(y);
x(nan_points) = [];
y(nan_points) = [];

% If one n_bins provided, double
if length(n_bins) == 1
    n_bins = repmat(n_bins,1,2);
end

% Set default heat/scatter threshold
if ~exist('point_thresh','var') || isempty(point_thresh)
    point_thresh = 10; % threshold for heatmap vs scatter
end

% Set bins and bin data
x_bins = linspace(min(x),max(x),n_bins(1)+1);
y_bins = linspace(min(y),max(y),n_bins(2)+1);

x_bin_centers = x_bins(1:end-1) + diff(x_bins)./2;
y_bin_centers = y_bins(1:end-1) + diff(y_bins)./2;

[bin_n,~,~,bin_x,bin_y] = histcounts2(x,y,x_bins,y_bins);
bin_n = bin_n'; % histcounts2 returns a transposed x/y?! that's dumb
bin_idx = sub2ind([size(bin_n,1),size(bin_n,2)],bin_y,bin_x);

% Set heat bins to image and scatter points to plot
heat_bins = bin_n > point_thresh;
scatter_bins = bin_n <= point_thresh; 
scatter_points = scatter_bins(bin_idx);

% TEMP: sqrt bin_n for skewed data?
bin_n = sqrt(bin_n);

% Convert bins to 255-scale color index, get color by density bin
density_col = hot(255);
bin_n_col_idx = round((bin_n - min(bin_n(:)))./(max(bin_n(:)) - min(bin_n(:)))*254)+1;
point_col_idx = bin_n_col_idx(bin_idx);
point_col = density_col(point_col_idx,:);

% Plot scatter plot colored by density
hold on; colormap(hot)

heat_points = imagesc(x_bin_centers,y_bin_centers,bin_n);
set(heat_points,'AlphaData',heat_bins);

scatter(x(scatter_points),y(scatter_points),5,point_col(scatter_points,:),'filled');




