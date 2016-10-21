function corr_diff = AP_2d_gaussian_fit_cost(gauss_fit_params,field_size,stim_trace,px_trace)

% gauss_fit_params = [x,y,sigma];
x = gauss_fit_params(1);
y = gauss_fit_params(2);
sigma = abs(gauss_fit_params(3));

% Make a gaussian filter that can be centered off the screen
[field_x,field_y] = meshgrid(1:field_size(2),1:field_size(1));
x_gauss = normpdf(field_x,x,sigma);
y_gauss = normpdf(field_y,y,sigma);
field_gaussian = mvnpdf([field_x(:),field_y(:)],[x,y],repmat(sigma,1,2));

predicted_px_trace = field_gaussian'*stim_trace;

curr_corr = corrcoef(predicted_px_trace,px_trace);
corr_diff = 1-curr_corr(2);


























