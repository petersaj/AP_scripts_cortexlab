function deconvolved_activity = AP_deconv_wf(activity,spikes_flag,sample_rate)
% deconvolved_activity = AP_deconv_wf(activity,spikes_flag)
%
% Deconvolve fluorescence from tetO-GC6s widefield to match spikes
% OR convolve spikes with kernel to match deconvolved widefield signal
%
% Kernel created in AP_deconv_wf_kernelfit
%
% activity = ND array, time in 2nd dim
% spikes_flag = 
% false (default), apply widefield > spikes kernel
% true, apply spike > deconved-widefield kernel
% sample_rate = sample rate of signal to deconvolve (default = 35Hz)

%% Set defaults

if ~exist('spikes_flag','var') || isempty(spikes_flag)
    spikes_flag = false;
elseif ~islogical(spikes_flag)
    error('spikes_flag isn''t a logical');
end

if ~exist('sample_rate','var')
    sample_rate = 35;
end

%% Load GCaMP6s widefield kernel

% Use the kernel within the directory with the function
% (tried regression from supersample - didn't change shape)
kernel_path = fileparts(mfilename('fullpath'));
kernel_fn = [kernel_path filesep 'gcamp6s_kernel.mat'];
load(kernel_fn);

kernel_sample_rate = 1/mean(diff(gcamp6s_kernel.regression_t));

% Choose kernel (fluor->spikes or spikes->fluor)
if ~spikes_flag
    kernel_cat = vertcat(gcamp6s_kernel.regression{:});
elseif spikes_flag
    kernel_cat = vertcat(gcamp6s_kernel.spikes_regression{:});
end

% Max-normalize, filter at nyquist, average, and squared-sum normalize (arbitrary)
kernel_cat_filt = lowpass(kernel_cat',kernel_sample_rate/2-1,kernel_sample_rate)';
kernel_mean = nanmean(kernel_cat_filt./vecnorm(kernel_cat_filt,2,2),1);

% Resample kernel to given sampling rate
% (native/default sampling rate = 35Hz)
% (get range as nearest integer sample)
kernel_resample_t_range = ...
    floor(sample_rate*max(abs(gcamp6s_kernel.regression_t)))/sample_rate;
kernel_resample_t = -kernel_resample_t_range:1/sample_rate:kernel_resample_t_range;
kernel_resample = interp1(gcamp6s_kernel.regression_t,kernel_mean,kernel_resample_t);
kernel = kernel_resample./norm(kernel_resample);


%% Not fancy way: 

% Remove NaNs for convolution (put them back after)
activity_nonan = activity;
activity_nonan(isnan(activity_nonan)) = 0;

% Deconvolve to same size with replication padding
deconvolved_activity = convn(padarray(activity_nonan, ...
    [0,floor(length(kernel)/2)],'replicate','both'), ...
    kernel,'valid');

% Put original NaNs back in
deconvolved_activity(isnan(activity)) = NaN;























