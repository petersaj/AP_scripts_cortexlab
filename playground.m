%% post to alyx

onLoad; % initializes missing-http

myAlyx = alyx.getToken([], 'andy', 'xpr1mnt');

use_date = datestr(now,'yyyy-mm-dd');
use_time = datestr(now,'HH:MM');

clear d
d.user = 'andy';

animals = {'AP030','AP031'};

water = [1.2,1.2];
for curr_animal = 1:length(animals)
    d.subject = animals{curr_animal};
    d.water_administered = water(curr_animal);
    d.date_time = sprintf('%sT%s',use_date,use_time);
    newWater = alyx.postData(myAlyx, 'water-administrations', d);
end

weight = [26.7,26.3];
for curr_animal = 1:length(animals)
    d.subject = animals{curr_animal};
    d.weight = weight(curr_animal);
    d.date_time = sprintf('%sT%s',use_date,use_time);
    newWeight = alyx.postData(myAlyx, 'weighings', d);
end


%% Deconvolve fluorescence
% This looks like it works - it just smooths it out too much

load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\gcamp_kernel\gcamp6s_kernel.mat');

gcamp6s_kernel_resample = interp1(gcamp6s_kernel_t,gcamp6s_kernel,gcamp6s_kernel_t(1):1/framerate:gcamp6s_kernel_t(end));
gcamp6s_kernel_long = padarray(gcamp6s_kernel,[0,size(fVdf,2)-length(gcamp6s_kernel)],0,'post');

% lowpassCutoff = 100; % Hz
% [b100s, a100s] = butter(2, lowpassCutoff/(framerate/2), 'low');
% gcamp6s_kernel_long_filt = single(filtfilt(b100s,a100s,double(gcamp6s_kernel_long)')');
       
softnorm = 0.1;
dfVdf = real(ifft(fft(fVdf').*conj(fft(gcamp6s_kernel_long'))./(abs(fft(gcamp6s_kernel_long')).^2 + softnorm)))';


%% Get px combined from predicted (for below)

% Settings to plot
% rxn_time_bins = {[0,0.1],[0.1,0.2],[0.2,0.3],[0.3,0.4]};
% rxn_time_bins = {[0,0.15],[0.15,0.4],[0.6,0.7]};
% rxn_time_bins = {[0,0.5],[0.5,1]};
rxn_time_bins = {[0.1,0.3],[0.6,0.7]};
% rxn_time_bins = {[0.1,0.3]};
% rxn_time_bins = {[0.6,0.7]};

normalize_px = false;

% Get major trial types
vis_L_trials_hit = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == -1;

vis_R_trials_hit = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == -1 & ...
    trial_choice_allcat == 1;

vis_L_trials_miss = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == -1 & ...
    trial_choice_allcat == -1;

vis_R_trials_miss = ...
    trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & ...
    trial_choice_allcat == 1;

zero_L_trials = ...
    trial_contrast_allcat == 0 & ...
    trial_choice_allcat == -1;

zero_R_trials = ...
    trial_contrast_allcat == 0 & ...
    trial_choice_allcat == 1;

trial_types = ...
    [vis_L_trials_hit, vis_R_trials_hit, ...
    vis_L_trials_miss, vis_R_trials_miss, ...
    zero_L_trials, zero_R_trials];

%%% Get and plot fluorescence
px_trial_types = nan(size(U_master,1),size(U_master,2), ...
    length(t_downsample_diff),size(trial_types,2),length(rxn_time_bins));
for curr_rxn = 1:length(rxn_time_bins)
    for curr_trial_type = 1:size(trial_types,2)
        
        curr_trials = trial_types(:,curr_trial_type) & ...
            move_t > rxn_time_bins{curr_rxn}(1) & ...
            move_t < rxn_time_bins{curr_rxn}(2);
        
        curr_data = fluor_allcat_residual(curr_trials,:,1:use_svs);

        % re-align to movement onset
        t_leeway = 0.5;
        leeway_samples = round(t_leeway*(sample_rate/downsample_factor));
        curr_move_idx = move_idx(curr_trials);
        for i = 1:size(curr_data,1)
            curr_data(i,:) = circshift(curr_data(i,:),-curr_move_idx(i)+leeway_samples,2);
        end
        
        curr_data_mean = squeeze(nanmean(curr_data,1))';
        
        curr_px = svdFrameReconstruct(U_master(:,:,1:use_svs),curr_data_mean);
        
        px_trial_types(:,:,:,curr_trial_type,curr_rxn) = curr_px;
        
    end
    AP_print_progress_fraction(curr_rxn,length(rxn_time_bins));
end


% Flip and combine trial types
% 1) visual hit (contra), 2) visual miss (ipsi), 3) zero (L)
px_combined = cat(4, ...
    (px_trial_types(:,:,:,1,:) + AP_reflect_widefield(px_trial_types(:,:,:,2,:)))./2, ...
    (px_trial_types(:,:,:,3,:) + AP_reflect_widefield(px_trial_types(:,:,:,4,:)))./2, ...
    (px_trial_types(:,:,:,5,:) + AP_reflect_widefield(px_trial_types(:,:,:,6,:)))./2);

px_combined_hemidiff = px_combined - AP_reflect_widefield(px_combined);

if normalize_px
    % Normalize by dividing by max of each frame
    px_dims = size(px_combined);
    px_combined = bsxfun(@rdivide,px_combined,max(abs(reshape( ...
         px_combined,[px_dims(1)*px_dims(2),1,px_dims(3:end)])),[],1));
     
    px_combined_hemidiff = bsxfun(@rdivide,px_combined_hemidiff,max(abs(reshape( ...
         px_combined_hemidiff,[px_dims(1)*px_dims(2),1,px_dims(3:end)])),[],1));
end


disp('done')



%% Comparing movement-aligned activity by reaction time
% PUT THIS INTO THE MAIN SCRIPT?

%%% Get differences between reaction times
compare_rxn = [1,2];
px_condition_diff = px_combined(:,:,:,:,compare_rxn(1)) - px_combined(:,:,:,:,compare_rxn(2));

% Get similarities: subtract off differences
px_condition_diff_1 = px_condition_diff;
px_condition_diff_1(px_condition_diff_1 < 0) = 0;
px_condition_similarity_1 = px_combined(:,:,:,:,compare_rxn(2)) - px_condition_diff_1;

px_condition_diff_2 = -px_condition_diff;
px_condition_diff_2(px_condition_diff_2 < 0) = 0;
px_condition_similarity_2 = px_combined(:,:,:,:,compare_rxn(2)) - px_condition_diff_2;

% Plot single condition
plot_condition = 1;
AP_image_scroll([px_condition_diff(:,:,:,plot_condition), ...
    px_condition_similarity_1(:,:,:,plot_condition), ...
    px_condition_similarity_2(:,:,:,plot_condition)],t_diff);
axis image; caxis([-1e-3,1e-3]); 
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k');

% Plot all concatenated
AP_image_scroll([reshape(permute(px_condition_diff,[1,4,2,3]),[],size(px_combined,2),size(px_combined,3)), ...
    reshape(permute(px_condition_similarity_1,[1,4,2,3]),[],size(px_combined,2),size(px_combined,3)), ...
    reshape(permute(px_condition_similarity_2,[1,4,2,3]),[],size(px_combined,2),size(px_combined,3))],t_diff);
axis image; caxis([-1e-3,1e-3]); 
colormap(brewermap([],'*RdBu'))


%%% Get differences between conditions
compare_conditions = [1,3];
use_rxn_bin = 1;

px_condition_diff = px_combined(:,:,:,compare_conditions(1),use_rxn_bin) - ...
    px_combined(:,:,:,compare_conditions(2),use_rxn_bin);

% Get similarities: subtract off differences
px_condition_diff_1 = px_condition_diff;
px_condition_diff_1(px_condition_diff_1 < 0) = 0;
px_condition_similarity_1 = px_combined(:,:,:,compare_conditions(1),use_rxn_bin) - px_condition_diff_1;

px_condition_diff_2 = -px_condition_diff;
px_condition_diff_2(px_condition_diff_2 < 0) = 0;
px_condition_similarity_2 = px_combined(:,:,:,compare_conditions(2),use_rxn_bin) - px_condition_diff_2;

AP_image_scroll([px_condition_diff, ...
    px_condition_similarity_1, ...
    px_condition_similarity_2],t_diff);
axis image; caxis([-1e-3,1e-3]); 
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k');










%% Ctx->Str using specific time window
% the idea is if you can predict the bump around movement if you only use
% around movement
% (after v->str load code)

n_vs = size(fluor_allcat,3);
n_rois = n_vs;

% Downsample (otherwise it's too much for regression) and d(smooth(fluor))
downsample_factor = 2;
t_downsample = linspace(t(1),t(end),round(length(t)/downsample_factor));

t_diff =  conv(t,[1,1]/2,'valid');
t_downsample_diff = conv(t_downsample,[1,1]/2,'valid');

wheel = interp1(t,wheel_allcat(:,:,1)',t_downsample_diff)';
[~,move_idx] = max(abs(wheel(:,:,1)) > 2,[],2);


% To use straight fluorescence
% smooth_factor = 1;
% fluor_allcat_downsamp = permute(interp1(t,permute(convn( ...
%     fluor_allcat,ones(1,smooth_factor)/smooth_factor,'same'),[2,1,3,4]),t_downsample_diff),[2,1,3,4]);

% % To use derivative
smooth_factor = 3;
fluor_allcat_downsamp = permute(interp1(t_diff,permute(diff(convn( ...
    fluor_allcat,ones(1,smooth_factor)/smooth_factor,'same'),[],2),[2,1,3,4]),t_downsample_diff),[2,1,3,4]);

% (there's a nan in one trial??)
fluor_allcat_downsamp(isnan(fluor_allcat_downsamp)) = 0;

mua_allcat_downsamp = permute(interp1(t,permute(mua_allcat,[2,1,3,4]),t_downsample_diff),[2,1,3,4]);

% Low-pass filter activity (where does this ~10 Hz crap come from?)
lowpassCutoff = 6; % Hz
[b100s, a100s] = butter(2, lowpassCutoff/((sample_rate/downsample_factor)/2), 'low');
fluor_allcat_downsamp_filt = filter(b100s,a100s,fluor_allcat_downsamp,[],2);
mua_allcat_downsamp_filt = filter(b100s,a100s,mua_allcat_downsamp,[],2);

% re-align to movement onset
fluor_allcat_downsamp_filt_move = fluor_allcat_downsamp_filt;
mua_allcat_downsamp_filt_move = mua_allcat_downsamp_filt;
t_leeway = 0.5;
leeway_samples = round(t_leeway*(sample_rate/downsample_factor));
for i = 1:size(fluor_allcat_downsamp_filt_move,1)
    fluor_allcat_downsamp_filt_move(i,:,:) = circshift(fluor_allcat_downsamp_filt_move(i,:,:),-move_idx(i)+leeway_samples,2);
    mua_allcat_downsamp_filt_move(i,:,:) = circshift(mua_allcat_downsamp_filt_move(i,:,:),-move_idx(i)+leeway_samples,2);
end

% Regress from fluor to MUA
kernel_t = [0,0];
kernel_frames = round(kernel_t(1)*sample_rate/downsample_factor):round(kernel_t(2)*sample_rate/downsample_factor);
zs = [false,false];
cvfold = 1;
lambda = 6;
return_constant = true;

% Set trials to do regression on (kernel then applied to all trials)
% (to use all - for cv)
% kernel_trials = true(size(fluor_allcat,1),1);
% (to use only correct leftward trials)
kernel_trials = trial_side_allcat == 1 & ...
    trial_choice_allcat == -1 & trial_contrast_allcat > 0;
% (to use only correct rightward trials)
% kernel_trials = trial_side_allcat == -1 & ...
%     trial_choice_allcat == 1 & trial_contrast_allcat > 0;
% (to only use incorrect/zero trials)
% kernel_trials = trial_contrast_allcat == 0 | ...
%     trial_side_allcat == trial_choice_allcat;
% (to only use particular choice trials)
% kernel_trials = trial_choice_allcat == 1;
% (to only use timing)
% kernel_trials = move_t > 0.6 & move_t < 0.7;

kernel_time = find(t_downsample_diff > 0,1);

mua_nonan_trials = ~squeeze(any(any(isnan(mua_allcat),2),4));
mua_allcat_predicted = nan(size(mua_allcat_downsamp_filt_move));
fluor_k_reshape = nan(n_rois,length(kernel_frames),n_depths);
for curr_depth = 1:n_depths
    
    curr_nonan_trials = mua_nonan_trials(:,curr_depth);
    curr_kernel_trials = kernel_trials & mua_nonan_trials(:,curr_depth);
    
    [fluor_k,curr_mua_predicted,explained_var] = ...
        AP_regresskernel(reshape(permute( ...
        fluor_allcat_downsamp_filt_move(curr_kernel_trials,kernel_time,:,:), ...
        [2,1,4,3]),[],n_rois)', ...
        reshape(permute(mua_allcat_downsamp_filt_move(curr_kernel_trials,kernel_time,curr_depth,:), ...
        [2,1,4,3]),[],1)',kernel_frames,lambda,zs,cvfold,return_constant);
    
    fluor_k_reshape(:,:,curr_depth) = reshape(fluor_k(1:end-1),n_rois,[]);
       
%     % To use predicted straight from regression (cross-validated)
%     mua_allcat_predicted(curr_kernel_trials,:,curr_depth,:) = ...
%         permute(reshape(curr_mua_predicted',length(t_downsample_diff),sum(curr_kernel_trials)),[2,1,4,3]);
    
    % To apply kernel from regressed subset to all (CV NOT APPLICABLE)
    % Create design matrix of all time-shifted regressors
    regressor_design = repmat(reshape(permute( ...
        fluor_allcat_downsamp_filt_move(curr_nonan_trials,:,:,:), ...
        [2,1,4,3]),[],n_rois), ...
        [1,1,length(kernel_frames)]);
    
    % Temporally shift each page
    for curr_kernel_frame = 1:length(kernel_frames)
        regressor_design(:,:,curr_kernel_frame) = ...
            circshift(regressor_design(:,:,curr_kernel_frame), ...
            [kernel_frames(curr_kernel_frame),0,0]);
    end
    
    regressor_design = ...
        [reshape(regressor_design,[],size(regressor_design,2)*size(regressor_design,3)) ...
        ones(size(regressor_design,1),1)]; 
    
    predicted_spikes = regressor_design*fluor_k;
    
    mua_allcat_predicted(curr_nonan_trials,:,curr_depth,:) = ...
        permute(reshape(predicted_spikes',length(t_downsample_diff),sum(curr_nonan_trials)),[2,1,4,3]);
    
    AP_print_progress_fraction(curr_depth,n_depths);
    
end

fluor_k_px = reshape(svdFrameReconstruct(U_master(:,:,1:200),reshape(fluor_k_reshape,200,[])), ...
    size(U_master,1),size(U_master,2),length(kernel_frames),[]);

AP_image_scroll(fluor_k_px);
caxis([-prctile(abs(fluor_k_px(:)),100),prctile(abs(fluor_k_px(:)),100)]);
colormap(brewermap([],'*RdBu'));
axis image;
AP_reference_outline('ccf_aligned','k');

fluor_k_px_max = squeeze(max(abs(fluor_k_px),[],3));
figure;imagesc(reshape(fluor_k_px_max,size(U_master,1),[]));
colormap(brewermap([],'*RdBu'));
caxis([-prctile(abs(fluor_k_px_max(:)),99),prctile(abs(fluor_k_px_max(:)),99)]);
axis image off;

% % Apply empirical static nonlinearity (just using the regressed bit)
% figure;
% mua_allcat_predicted_nlin = mua_allcat_predicted;
% for curr_depth = 1:n_depths       
%     
%     measured_data = reshape(mua_allcat_downsamp_filt_move(kernel_trials,kernel_time,curr_depth),[],1);
%     predicted_data = reshape(mua_allcat_predicted(kernel_trials,kernel_time,curr_depth),[],1);
%     predicted_data_total = reshape(mua_allcat_predicted(:,:,curr_depth),[],1);
%     
%     n_bins = 300;
%     activity_bounds = linspace(-1,6,n_bins+1);
%     activity_bin_centers = conv2(activity_bounds,[1,1]/2,'valid');
%   
%     measured_bins = discretize(measured_data,activity_bounds);
%     predicted_bins = discretize(predicted_data,activity_bounds);
%     predicted_total_bins = discretize(predicted_data_total,activity_bounds);
%     
%     measured_data_binmean = accumarray(predicted_bins(~isnan(predicted_bins)), ...
%         measured_data(~isnan(predicted_bins)),[n_bins,1],@mean,nan);
%     predicted_data_binmean = accumarray(predicted_bins(~isnan(predicted_bins)),...
%         predicted_data(~isnan(predicted_bins)),[n_bins,1],@mean,nan);
%     
%     predicted_data_nlin = nan(size(predicted_data_total));
%     predicted_data_nlin(~isnan(predicted_total_bins)) = measured_data_binmean(predicted_total_bins(~isnan(predicted_total_bins)));
%     
%     predicted_data_nlin_bins = discretize(predicted_data_nlin,activity_bounds);
%     predicted_data_nlin_binmean = accumarray( ...
%         predicted_bins(~isnan(predicted_bins)), ...
%         predicted_data_nlin(~isnan(predicted_bins)),[n_bins,1],@mean,nan);
%     
%     mua_allcat_predicted_nlin(:,:,curr_depth) = ...
%         reshape(predicted_data_nlin, ...
%         size(mua_allcat_predicted,1),size(mua_allcat_predicted,2));
%     
%     subplot(1,n_depths,curr_depth); hold on;
%     plot(measured_data_binmean,predicted_data_binmean,'linewidth',2);
%     plot(measured_data_binmean,predicted_data_nlin_binmean,'linewidth',2);
%     xlim([-2,5]);ylim([-2,5]);
%     line(xlim,ylim,'color','k');
%     xlabel('Measured')
%     ylabel('Predicted')
%     axis square;
% end

% Apply empirical static nonlinearity (whole dataset)
figure;
mua_allcat_predicted_nlin = mua_allcat_predicted;
for curr_depth = 1:n_depths       
    measured_data = reshape(mua_allcat_downsamp_filt_move(:,:,curr_depth),[],1);
    predicted_data = reshape(mua_allcat_predicted(:,:,curr_depth),[],1);
    
    n_bins = 300;
    activity_bounds = linspace(-1,6,n_bins+1);
    activity_bin_centers = conv2(activity_bounds,[1,1]/2,'valid');
  
    measured_bins = discretize(measured_data,activity_bounds);
    predicted_bins = discretize(predicted_data,activity_bounds);
    
    measured_data_binmean = accumarray(predicted_bins(~isnan(predicted_bins)), ...
        measured_data(~isnan(predicted_bins)),[n_bins,1],@mean,nan);
    predicted_data_binmean = accumarray(predicted_bins(~isnan(predicted_bins)),...
        predicted_data(~isnan(predicted_bins)),[n_bins,1],@mean,nan);
    
    predicted_data_nlin = nan(size(predicted_data));
    predicted_data_nlin(~isnan(predicted_bins)) = measured_data_binmean(predicted_bins(~isnan(predicted_bins)));
    
    predicted_data_nlin_bins = discretize(predicted_data_nlin,activity_bounds);
    predicted_data_nlin_binmean = accumarray( ...
        predicted_bins(~isnan(predicted_bins)), ...
        predicted_data_nlin(~isnan(predicted_bins)),[n_bins,1],@mean,nan);
    
    mua_allcat_predicted_nlin(:,:,curr_depth) = ...
        reshape(predicted_data_nlin, ...
        size(mua_allcat_predicted,1),size(mua_allcat_predicted,2));
    
    subplot(1,n_depths,curr_depth); hold on;
    plot(measured_data_binmean,predicted_data_binmean,'linewidth',2);
    plot(measured_data_binmean,predicted_data_nlin_binmean,'linewidth',2);
    xlim([-2,5]);ylim([-2,5]);
    line(xlim,ylim,'color','k');
    xlabel('Measured')
    ylabel('Predicted')
    axis square;
end

% Get explained variance within time
mua_allcat_residual = mua_allcat_downsamp_filt_move - mua_allcat_predicted;

t_var = t_downsample_diff > -inf & t_downsample_diff < inf;
mua_sse_measured = squeeze(nansum(mua_allcat_downsamp_filt_move(:,t_var,:).^2,1));
mua_sse_residual = squeeze(nansum(mua_allcat_residual(:,t_var,:).^2,1));
mua_expl_var = (mua_sse_measured - mua_sse_residual)./mua_sse_measured;
figure;
plot((sum(mua_sse_measured,1)-sum(mua_sse_residual,1))./ ...
    sum(mua_sse_measured,1),'k','linewidth',2)
ylabel('Explained variance');
xlabel('Striatum depth');

% Get conditions for each trial, plot selected
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];

contrast_side_col = colormap_BlueWhiteRed(5);
contrast_side_col(6,:) = 0;
contrast_side_val = unique(sort([-contrasts,contrasts]))';

% contrast, side, choice
plot_conditions = ...
    [contrasts,contrasts; ...
    -ones(1,6),-1,ones(1,5); ...
    ones(1,6),-ones(1,6)]';
%     plot_conditions = ...
%         [0.125,1,0.125,1; ...
%         -1,-1,1,1; ...
%         1,1,-1,-1]';
%     % plot_conditions = ...
%     [contrasts(2:end),contrasts(2:end); ...
%     ones(1,10); ...
%     ones(1,5),-ones(1,5)]';
% plot_conditions = ...
%     [contrasts(2:end),contrasts(2:end); ...
%     -ones(1,5),ones(1,5); ...
%     ones(1,5),-ones(1,5)]';
% plot_conditions = ...
%     [contrasts(2:end),contrasts(2:end); ...
%     -ones(1,5),ones(1,5); ...
%     ones(1,5),ones(1,5)]';
% plot_conditions = ...
%     [0,0; ...
%     -1,-1; ...
%     -1,1]';

use_rxn = move_t > 0 & move_t < 0.5;

[~,plot_id] = ismember( ...
    [trial_contrast_allcat,trial_side_allcat,trial_choice_allcat], ...
    plot_conditions,'rows');

% Plot striatum
figure; hold on;
for curr_plot = 1:3
    
    switch curr_plot
        case 1
            plot_data = mua_allcat_downsamp_filt_move;
            plot_title = 'Measured';
        case 2
            plot_data = mua_allcat_predicted_nlin;
            plot_title = 'Predicted';
        case 3
            plot_data = mua_allcat_downsamp_filt_move - mua_allcat_predicted_nlin;
            plot_title = 'Residual';
    end
    
    p_str(curr_plot) = subplot(1,3,curr_plot); hold on;
    for curr_plot_condition = 1:size(plot_conditions,1)
        
        curr_trials = plot_id == curr_plot_condition & use_rxn;
        curr_data = plot_data(curr_trials,:,:);
        
%         % re-align to movement onset
%         t_leeway = 0.5;
%         leeway_samples = round(t_leeway*(sample_rate/downsample_factor));
%         curr_move_idx = move_idx(curr_trials);
%         for i = 1:size(curr_data,1)
%             curr_data(i,:) = circshift(curr_data(i,:),-curr_move_idx(i)+leeway_samples,2);
%         end
        
        curr_data_mean = squeeze(nanmean(curr_data,1));
        
        curr_contrast_side = max(trial_contrast_allcat(curr_trials))*sign(max(trial_side_allcat(curr_trials)));
        contrast_side_idx = find(curr_contrast_side == contrast_side_val);
        curr_col = contrast_side_col(contrast_side_idx,:);
        
        if curr_contrast_side == 0
            switch max(trial_choice_allcat(curr_trials))
                case -1
                    curr_col = 'm';
                case 1
                    curr_col = 'c';
            end
        end
        
        AP_stackplot(curr_data_mean,t_downsample_diff,2,false,curr_col,1:n_depths);
        
    end
    line([0,0],ylim,'color','k');
    xlabel('Time from stim');
    title(plot_title);
end
linkaxes(p_str)

%% FROM ABOVE: do kernels at each time point

% Regress from fluor to MUA
kernel_t = [0,0];
kernel_frames = round(kernel_t(1)*sample_rate/downsample_factor):round(kernel_t(2)*sample_rate/downsample_factor);
zs = [false,false];
cvfold = 5;
lambda = 2;
return_constant = false;
use_constant = false;

% Set trials to do regression on (kernel then applied to all trials)
% (to use all - for cv)
% kernel_trials = true(size(fluor_allcat,1),1);
% (to use only correct leftward trials)
% kernel_trials = trial_side_allcat == 1 & ...
%     trial_choice_allcat == 1 & trial_contrast_allcat > 0;
% (to use only correct rightward trials)
% kernel_trials = trial_side_allcat == -1 & ...
%     trial_choice_allcat == 1 & trial_contrast_allcat > 0;
% (to only use incorrect/zero trials)
% kernel_trials = trial_contrast_allcat == 0 | ...
%     trial_side_allcat == trial_choice_allcat;
% (to only use particular choice trials)
% kernel_trials = trial_contrast_allcat == 0;
% (to only use timing)
% kernel_trials = move_t > 0.6 & move_t < 0.7;

kernel_trials = move_t > 0 & move_t < 0.5 & ...
    trial_choice_allcat == -1 & trial_contrast_allcat == 0;

mua_trials_predicted = nan(size(mua_allcat_downsamp));

fluor_k_px_max_all = nan(size(U_master,1),size(U_master,2)*n_depths,length(t_downsample_diff));
k_constant_all = nan(length(t_downsample_diff),n_depths);
for kernel_time = 1:length(t_downsample_diff)
    
    mua_nonan_trials = ~squeeze(any(any(isnan(mua_allcat),2),4));
    mua_allcat_predicted = nan(size(mua_allcat_downsamp_filt_move));
    fluor_k_reshape = nan(n_rois,length(kernel_frames),n_depths);
    for curr_depth = 1:n_depths
        
        curr_nonan_trials = mua_nonan_trials(:,curr_depth);
        curr_kernel_trials = kernel_trials & mua_nonan_trials(:,curr_depth);
        
        [fluor_k,curr_mua_predicted,explained_var] = ...
            AP_regresskernel(reshape(permute( ...
            fluor_allcat_downsamp_filt_move(curr_kernel_trials,kernel_time,:,:), ...
            [2,1,4,3]),[],n_rois)', ...
            reshape(permute(mua_allcat_downsamp_filt_move(curr_kernel_trials,kernel_time,curr_depth,:), ...
            [2,1,4,3]),[],1)',kernel_frames,lambda,zs,cvfold,return_constant,use_constant);
        
        fluor_k_reshape(:,:,curr_depth) = reshape(fluor_k(1:end-use_constant),n_rois,[]);
        k_constant_all(kernel_time,curr_depth) = fluor_k(end);
        
        % To use predicted straight from regression (cross-validated)
        mua_trials_predicted(curr_kernel_trials,kernel_time,curr_depth) = curr_mua_predicted;
        
%         % To apply kernel from regressed subset to all (CV NOT APPLICABLE)
%         % Create design matrix of all time-shifted regressors
%         regressor_design = repmat(reshape(permute( ...
%             fluor_allcat_downsamp_filt_move(curr_nonan_trials,:,:,:), ...
%             [2,1,4,3]),[],n_rois), ...
%             [1,1,length(kernel_frames)]);
%         
%         % Temporally shift each page
%         for curr_kernel_frame = 1:length(kernel_frames)
%             regressor_design(:,:,curr_kernel_frame) = ...
%                 circshift(regressor_design(:,:,curr_kernel_frame), ...
%                 [kernel_frames(curr_kernel_frame),0,0]);
%         end
%         
%         regressor_design = ...
%             [reshape(regressor_design,[],size(regressor_design,2)*size(regressor_design,3)) ...
%             ones(size(regressor_design,1),1)];
%         
%         predicted_spikes = regressor_design*fluor_k;
%         
%         mua_allcat_predicted(curr_nonan_trials,:,curr_depth,:) = ...
%             permute(reshape(predicted_spikes',length(t_downsample_diff),sum(curr_nonan_trials)),[2,1,4,3]);
        
%         AP_print_progress_fraction(curr_depth,n_depths);
        
    end
    
    fluor_k_px = reshape(svdFrameReconstruct(U_master(:,:,1:200),reshape(fluor_k_reshape,200,[])), ...
        size(U_master,1),size(U_master,2),length(kernel_frames),[]);  
    fluor_k_px_max = squeeze(max(abs(fluor_k_px),[],3));

    fluor_k_px_max_all(:,:,kernel_time) = reshape(fluor_k_px_max,size(U_master,1),[]);
    
    AP_print_progress_fraction(kernel_time,length(t_downsample_diff));
end

% Plot cortex kernels
AP_image_scroll(fluor_k_px_max_all,t_downsample_diff);
axis image off
caxis([-prctile(abs(fluor_k_px_max_all(:)),100),prctile(abs(fluor_k_px_max_all(:)),100)]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k',[],[size(U_master,1),size(U_master,2),1,4])

% Plot constant if used
if return_constant
    figure; hold on; set(gca,'ColorOrder',copper(n_depths));
    plot(t_downsample_diff,k_constant_all,'linewidth',2);
    ylabel('Constant term');
end

% Plot predicted over real
avg_measured = squeeze(nanmean(mua_allcat_downsamp_filt_move(kernel_trials,:,:),1));
avg_predicted = squeeze(nanmean(mua_trials_predicted,1));
figure; hold on; set(gca,'ColorOrder',copper(n_depths));
plot(t_downsample_diff,avg_measured,'linewidth',2);
plot(t_downsample_diff,avg_predicted,'--','linewidth',2);
ylabel('Average MUA');
xlabel('Time from move');



%% KEEP THIS: fixing MUA downsample to be summed and not interp

%%%%%%%%%%%%%%%%%%%%

downsample_factor = 4;
t_downsample = mean(reshape(t(1:end-mod(length(t), ...
    downsample_factor)),downsample_factor,[]),1);
t_downsample_diff = conv(t_downsample,[1,1]/2,'valid');

if load_task
    wheel = interp1(t,wheel_allcat(:,:,1)',t_downsample_diff)';
    [~,move_idx] = max(abs(wheel(:,:,1)) > 2,[],2);
end

mua_allcat_downsamp = permute(reshape(squeeze(sum(reshape(reshape(permute( ...
    mua_allcat(:,1:end-mod(size(mua_allcat,2),downsample_factor),:), ...
    [2,1,3]),[],n_depths),downsample_factor,[],n_depths),1)), ...
    [],size(mua_allcat,1),n_depths),[2,1,3]);

fluor_allcat_downsamp = permute(reshape(squeeze(sum(reshape(reshape(permute( ...
    fluor_allcat(:,1:end-mod(size(fluor_allcat,2),downsample_factor),:), ...
    [2,1,3]),[],n_rois),downsample_factor,[],n_rois),1)), ...
    [],size(fluor_allcat,1),n_rois),[2,1,3]);

smooth_factor = 3;
fluor_allcat_downsamp = diff(convn( ...
    fluor_allcat_downsamp,ones(1,smooth_factor)/smooth_factor,'same'),[],2);
mua_allcat_downsamp = permute(interp1(t_downsample, ...
    permute(mua_allcat_downsamp,[2,1,3,4]),t_downsample_diff),[2,1,3,4]);

% (there's a nan in one trial??)
fluor_allcat_downsamp(isnan(fluor_allcat_downsamp)) = 0;

%%%%%%%%%%%%%%%%%%%%


%% Hemidiff in V (V-V_mirror) - keep for future reference
% looks good 

U_r = reshape(U_master(:,:,1:n_vs),[],n_vs);
U_mirror_r = reshape(AP_reflect_widefield(U_master(:,:,1:n_vs)),[],n_vs);
mirror_matrix = U_r'*U_mirror_r;
fluor_allcat_mirror = reshape(transpose( ...
    mirror_matrix*reshape(fluor_allcat,[],n_vs)'),size(fluor_allcat));

% sanity check: plot regular and mirrored
use_trials = trial_contrast_allcat > 0 & trial_side_allcat == 1 & ...
    trial_choice_allcat == -1 & move_t > 0.1 & move_t < 0.3;

V_mean = squeeze(nanmean(fluor_allcat(use_trials,:,:),1))';
V_mirror_mean = squeeze(nanmean(fluor_allcat_mirror(use_trials,:,:),1))';

px_mean = svdFrameReconstruct(U_master(:,:,1:n_vs),V_mean);
px_mirror_mean = svdFrameReconstruct(U_master(:,:,1:n_vs),V_mirror_mean);

AP_image_scroll([px_mean,px_mirror_mean],t);
caxis([-prctile(abs(px_mean(:)),99),prctile(abs(px_mean(:)),99)]);
colormap(brewermap([],'*RdBu'));
axis image;

V_hemidiff = V_mean - V_mirror_mean;
px_hemidiff = svdFrameReconstruct(U_master(:,:,1:n_vs),V_hemidiff);
AP_image_scroll(px_hemidiff,t);
caxis([-prctile(abs(px_hemidiff(:)),99),prctile(abs(px_hemidiff(:)),99)]);
colormap(brewermap([],'*RdBu'));
axis image;

% sanity check
use_trials = trial_contrast_allcat > 0 & trial_side_allcat == 1 & ...
    trial_choice_allcat == -1 & move_t > 0.1 & move_t < 0.3;

V_check = squeeze(nanmean(fluor_allcat_downsamp_filt(use_trials,:,:),1))';
px_check = svdFrameReconstruct(U_master(:,:,1:n_vs),V_check);

AP_image_scroll(px_check,t);
caxis([-prctile(abs(px_check(:)),100),prctile(abs(px_check(:)),100)]);
colormap(brewermap([],'*RdBu'));
axis image;
AP_reference_outline('ccf_aligned','k');


% compatible with code: 

mirror_matrix = reshape(U_master(:,:,1:n_vs),[],n_vs)'* ...
    reshape(AP_reflect_widefield(U_master(:,:,1:n_vs)),[],n_vs);
fluor_allcat_downsamp_mirror = reshape(transpose( ...
    mirror_matrix*reshape(fluor_allcat_downsamp,[],n_vs)'),size(fluor_allcat_downsamp));
fluor_allcat_downsamp = fluor_allcat_downsamp - fluor_allcat_downsamp_mirror;

















