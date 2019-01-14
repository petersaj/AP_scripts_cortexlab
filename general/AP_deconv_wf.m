function deconv_fluorescence = AP_deconv_wf(fluorescence)
% deconv_fluorescence = AP_deconv_wf(fluorescence)
%
% Deconvolve widefield from tetO-GC6s mice
% fluorescence = ND array, time in 2nd dim

% Load GCaMP6s widefield kernel
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\gcamp_kernel\gcamp6s_kernel.mat');
gcamp6s_kernel_cat = vertcat(gcamp6s_kernel.regression{:});
gcamp6s_kernel_mean = nanmean(gcamp6s_kernel_cat./max(abs(gcamp6s_kernel_cat),[],2),1);

% Use -300:300ms of the kernel
use_kernel_t = abs(gcamp6s_kernel.regression_t) < 0.3;
gcamp6s_kernel_mean = gcamp6s_kernel_mean(use_kernel_t);

% Get valid time from kernel (assume time-symmetric)
t_rel = 1:size(fluorescence,2);
t_rel_conv_valid = round(conv(t_rel,ones(1,length(gcamp6s_kernel_mean))./length(gcamp6s_kernel_mean),'valid'));

% Remove NaNs to make convolution faster (put back in later)
fluorescence_nonan = fluorescence;
fluorescence_nonan(isnan(fluorescence_nonan)) = 0;

% Deconvolve widefield, set invalid points to NaN
dims = 1:ndims(fluorescence);
dims_permute = circshift(dims,-1);
[~,dims_unpermute] = sort(dims_permute);

deconv_fluorescence = permute(interp1(t_rel_conv_valid, ...
    permute(convn(fluorescence_nonan,gcamp6s_kernel_mean,'valid'), ...
    dims_permute),t_rel),dims_unpermute);

% Put original NaNs back in
deconv_fluorescence(isnan(fluorescence)) = NaN;












