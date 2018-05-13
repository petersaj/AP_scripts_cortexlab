%% Estimate the best upsample factor 
% This is unused at the moment: the question was to look at what binning
% scale the cortex was best predictive of the striatum. In one dataset, it
% looked like this was ~ 0.7x framerate (~400 ms bins). I think I still
% want higher temporal resolution of the cortex so upsampling is still
% legal even if the explained variance goes down, since it's a slightly
% different question

% upsample_factors = linspace(0.01,2,20);
% 
% explained_var_upsamples = nan(length(upsample_factors),1);
% 
% if plot_lambda_estimate
%     figure; hold on;
%     xlabel('Upsample factor');
%     ylabel('Explained variance');
%     drawnow;
%     curr_plot = plot(upsample_factors,explained_var_upsamples,'linewidth',2);
% end
% 
% for curr_upsample_idx = 1:length(upsample_factors)
%     
%     curr_upsample_factor = upsample_factors(curr_upsample_idx);
%     
%     sample_rate = (1/median(diff(frame_t)))*curr_upsample_factor;
%     
%     % Skip the first n seconds to do this
%     skip_seconds = 60;
%     time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
%     time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
%     
%     % Use all spikes in striatum
%     use_spikes = spike_times_timeline(ismember(spike_templates, ...
%         find(templateDepths > str_depth(1) & templateDepths <= str_depth(2))));
%     binned_spikes = single(histcounts(use_spikes,time_bins));
%     
%     use_svs = 1:50;
%     kernel_t = [-0.3,0.3];
%     kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
%     zs = [false,true];
%     lambda = 0;
%     cvfold = 2;
%     
%     % Resample and get derivative of V
%     dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
%         diff(fVdf(use_svs,:),[],2)',time_bin_centers)';
%     
%     [~,predicted_spikes,explained_var] = ...
%         AP_regresskernel(dfVdf_resample, ...
%         binned_spikes,kernel_frames,lambda,zs,cvfold);
%     
%     explained_var_upsamples(curr_upsample_idx) = explained_var.total;
%     
%     if plot_lambda_estimate
%         set(curr_plot,'YData',explained_var_upsamples);
%         drawnow;
%     end
%     
% end


%% Estimate lambda between widefield and spikes via cross-validation
% (good to make general in the future)

upsample_factor = 2;
sample_rate = (1/median(diff(frame_t)))*upsample_factor;

% Skip the first n seconds to do this
skip_seconds = 60;
time_bins = frame_t(find(frame_t > skip_seconds,1)):1/sample_rate:frame_t(find(frame_t-frame_t(end) < -skip_seconds,1,'last'));
time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

% Use all spikes in striatum
use_spikes = spike_times_timeline(ismember(spike_templates, ...
    find(templateDepths > str_depth(1) & templateDepths <= str_depth(2))));
binned_spikes = single(histcounts(use_spikes,time_bins));

use_svs = 1:50;
kernel_t = [-0.3,0.3];
kernel_frames = round(kernel_t(1)*sample_rate):round(kernel_t(2)*sample_rate);
zs = [false,true]; % MUA has to be z-scored if large lambda
cvfold = 2;

% Resample and get derivative of V
dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
    diff(fVdf(use_svs,:),[],2)',time_bin_centers)';

n_update_lambda = 1;
lambda_range = [3,7]; % ^10
n_lambdas = 50;

if plot_lambda_estimate
    figure; hold on;
    set(gca,'XScale','log');
    xlabel('\lambda');
    ylabel('Explained variance');
    drawnow;
end

for curr_update_lambda = 1:n_update_lambda
    
    lambdas = logspace(lambda_range(1),lambda_range(2),n_lambdas);
    explained_var_lambdas = nan(n_lambdas,1);
    
    if plot_lambda_estimate
        curr_plot = plot(lambdas,explained_var_lambdas,'linewidth',2);
    end
    
    for curr_lambda_idx = 1:length(lambdas)
        
        curr_lambda = lambdas(curr_lambda_idx);
        
        [~,predicted_spikes,explained_var] = ...
            AP_regresskernel(dfVdf_resample, ...
            binned_spikes,kernel_frames,curr_lambda,zs,cvfold);
        
        explained_var_lambdas(curr_lambda_idx) = explained_var.total;
        
        if plot_lambda_estimate
            set(curr_plot,'YData',explained_var_lambdas);
            drawnow;
        end
        
    end
    
    lambda_bin_size = diff(lambda_range)/n_lambdas;
    explained_var_lambdas_smoothed = smooth(explained_var_lambdas,3);
    [best_lambda_explained_var,best_lambda_idx] = max(explained_var_lambdas_smoothed);
    best_lambda = lambdas(best_lambda_idx);
    lambda_range = log10([best_lambda,best_lambda]) + ...
        [-lambda_bin_size,lambda_bin_size];
    
end

if plot_lambda_estimate
    plot(lambdas,explained_var_lambdas_smoothed,'r');
    line(xlim,repmat(best_lambda_explained_var,1,2),'color','k');
    line(repmat(best_lambda,1,2),ylim,'color','k');
    title(['\lambda = ' num2str(best_lambda)]);
    disp(['Best lambda = ' num2str(best_lambda) ', Frac var = ' num2str(explained_var_lambdas(best_lambda_idx))]);
    drawnow;
end

lambda = best_lambda;












