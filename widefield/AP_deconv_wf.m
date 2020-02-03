function deconv_fluorescence = AP_deconv_wf(fluorescence)
% deconv_fluorescence = AP_deconv_wf(fluorescence)
%
% Deconvolve widefield from tetO-GC6s mice
% fluorescence = ND array, time in 2nd dim

%% Load GCaMP6s widefield kernel

% Use the kernel within the directory with the function
kernel_path = fileparts(mfilename('fullpath'));
kernel_fn = [kernel_path filesep 'gcamp6s_kernel.mat'];
load(kernel_fn);

gcamp6s_kernel_cat = vertcat(gcamp6s_kernel.regression{:});
gcamp6s_kernel_mean = nanmean(gcamp6s_kernel_cat./max(abs(gcamp6s_kernel_cat),[],2),1);

%% Not fancy way: 

% Remove NaNs for convolution (put them back after)
fluorescence_nonan = fluorescence;
fluorescence_nonan(isnan(fluorescence_nonan)) = 0;

% Deconvolve to same size with replication padding
deconv_fluorescence = convn(padarray(fluorescence_nonan, ...
    [0,floor(length(gcamp6s_kernel_mean)/2)],'replicate','both'), ...
    gcamp6s_kernel_mean,'valid');

% Put original NaNs back in
deconv_fluorescence(isnan(fluorescence)) = NaN;

%% Fancy way: filter kernel, use only valid convolution points

% % Filter kernel
% sample_rate = 1/mean(diff(gcamp6s_kernel.regression_t));
% lowpassCutoff = 10;
% [b100s, a100s] = butter(2, lowpassCutoff/((sample_rate)/2),'low');
% gcamp6s_kernel_mean_filt = single(filtfilt(b100s,a100s,double(gcamp6s_kernel_mean)));
% 
% % Use -300:300ms of the kernel
% use_kernel_t = abs(gcamp6s_kernel.regression_t) < 0.3;
% gcamp6s_kernel_mean_use = gcamp6s_kernel_mean_filt(use_kernel_t);
% 
% % Get valid time from kernel (assume time-symmetric)
% t_rel = 1:size(fluorescence,2);
% t_rel_conv_valid = round(conv(t_rel,ones(1,length(gcamp6s_kernel_mean_use))./length(gcamp6s_kernel_mean_use),'valid'));
% 
% % Remove NaNs to make convolution faster (put back in later)
% fluorescence_nonan = fluorescence;
% fluorescence_nonan(isnan(fluorescence_nonan)) = 0;
% 
% % Deconvolve widefield, set invalid points to NaN
% dims = 1:ndims(fluorescence);
% dims_permute = circshift(dims,-1);
% [~,dims_unpermute] = sort(dims_permute);
% 
% deconv_fluorescence = permute(interp1(t_rel_conv_valid, ...
%     permute(convn(fluorescence_nonan,gcamp6s_kernel_mean_use,'valid'), ...
%     dims_permute),t_rel),dims_unpermute);
% 
% % Put original NaNs back in
% deconv_fluorescence(isnan(fluorescence)) = NaN;





















