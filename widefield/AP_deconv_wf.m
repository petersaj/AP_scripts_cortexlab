function convolved_activity = AP_deconv_wf(activity,spikes_flag,sample_rate)
% filtered_activity = AP_deconv_wf(activity,spikes_flag)
%
% Convolved fluorescence from tetO-GC6s widefield to match spikes
% OR convolve spikes with kernel to match deconvolved widefield signal
%
% Kernel created in AP_deconv_wf_kernelfit
%
% activity = ND array, time in 2nd dim
% spikes_flag = 
% false (default), apply widefield > spikes kernel
% true, apply spike > deconved-widefield kernel

%% Set defaults

if ~exist('spikes_flag','var')
    spikes_flag = false;
elseif ~islogical(spikes_flag)
    error('spikes_flag isn''t a logical');
end

if ~exist('sample_rate','var')
    sample_rate = []; 
end

%% Load GCaMP6s widefield kernel

% Use the kernel within the directory with the function
% (tried regression from supersample - didn't change shape)
kernel_path = fileparts(mfilename('fullpath'));
if isempty(sample_rate)
    kernel_fn = [kernel_path filesep 'gcamp6s_kernel.mat'];
else
    kernel_fn = [kernel_path filesep 'gcamp6s_kernel_' num2str(sample_rate) 'Hz.mat'];
end
load(kernel_fn);

kernel_sample_rate = 1/mean(diff(gcamp6s_kernel.regression_t));

% Max-normalize, filter at nyquist, average, and squared-sum normalize (arbitrary)
if ~spikes_flag
    kernel_cat = vertcat(gcamp6s_kernel.regression{:});
    kernel_cat_filt = lowpass(kernel_cat',kernel_sample_rate/2-1,kernel_sample_rate)';
    kernel_mean = nanmean(kernel_cat_filt./max(abs(kernel_cat_filt),[],2),1);
    kernel = kernel_mean./norm(kernel_mean);
elseif spikes_flag
    kernel_cat = vertcat(gcamp6s_kernel.spikes_regression{:});
    kernel_cat_filt = lowpass(kernel_cat',kernel_sample_rate/2-1,kernel_sample_rate)';
    kernel_mean = nanmean(kernel_cat_filt./max(abs(kernel_cat_filt),[],2),1);
    kernel = kernel_mean./norm(kernel_mean);
end


%% Not fancy way: 

% Remove NaNs for convolution (put them back after)
activity_nonan = activity;
activity_nonan(isnan(activity_nonan)) = 0;

% Deconvolve to same size with replication padding
convolved_activity = convn(padarray(activity_nonan, ...
    [0,floor(length(kernel)/2)],'replicate','both'), ...
    kernel,'valid');

% Put original NaNs back in
convolved_activity(isnan(activity)) = NaN;























