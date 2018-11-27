% The original scripts here were in test_wf_ephys_choiceworld_analysis

%% TO DO:
% - ctx->str alignment now in separate script, running
% - put other batch preprocessing into separate script

%% Fig 1b: Example traces

warning('Probably not best example, check others');

% Load and align
str_align = 'kernel';
animal = 'AP028'; 
day = '2017-12-16'; 
experiment = 1; 
verbose = false; 
AP_load_experiment;

avg_im_aligned = AP_align_widefield(animal,day,avg_im);
Udf_aligned = single(AP_align_widefield(animal,day,Udf));

% Define ROIs and get fluorescence traces
roi_circle_size = 20;
roi_x = [131,174,110,51];
roi_y = [297,96,71,144];
[x,y] = meshgrid(1:size(avg_im_aligned,1),1:size(avg_im_aligned,2));
roi_mask = cell2mat(arrayfun(@(roi) sqrt((x-roi_x(roi)).^2 + (y-roi_y(roi)).^2) <= ...
    roi_circle_size,permute(1:length(roi_x),[1,3,2]),'uni',false));
roi_trace = AP_svd_roi(Udf_aligned,fVdf,[],[],roi_mask);

roi_trace_deriv = diff(roi_trace,[],2);
roi_trace_deriv(roi_trace_deriv < 0) = 0;
frame_t_deriv = conv(frame_t,[1,1]/2,'valid');

% Bin spikes by aligned depth
n_depths = n_aligned_depths;
depth_group = aligned_str_depth_group;

time_bins = [frame_t_deriv,frame_t_deriv(end)+1/framerate];
binned_spikes = zeros(n_depths,length(time_bins)-1);
for curr_depth = 1:n_depths
    curr_spike_times = spike_times_timeline(depth_group == curr_depth);
    binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
end

% Plot ROIs and traces
figure;
subplot(1,6,1);

roi_boundaries = bwboundaries(sum(roi_mask,3));
imagesc(avg_im_aligned);colormap(gray);
caxis([0,prctile(avg_im_aligned(:),99)]);
axis image;
AP_reference_outline('ccf_aligned','r');
p = cellfun(@(x) plot(x(:,2),x(:,1),'b','linewidth',2),roi_boundaries);

subplot(1,6,2:6); hold on;
p1 = AP_stackplot(bsxfun(@rdivide,binned_spikes,std(binned_spikes,[],2))', ...
    frame_t_deriv,10,false,'k');
p2 = AP_stackplot(bsxfun(@rdivide,roi_trace_deriv,std(roi_trace_deriv,[],2))', ...
    frame_t_deriv,10,false,[0,0.7,0]);
xlabel('Time (seconds)');
ylabel('Activity (std)');
legend([p1(1),p2(1)],{'MUA','\DeltaFluorescence'});
xlim([177,200]);

%% Fig 1b: Average regression maps
% (half-done: variable names etc will change when new code finished)

n_aligned_depths = 4;

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_ephys';

k_fn = [data_path filesep 'wf_ephys_maps_concat_' num2str(n_aligned_depths) '_depths_kernel'];
load(k_fn);

% (scale r_px's because different lambdas give different weights)
% (do in a loop because memory can't handle a cellfun??)
k_px = nan(437,416,43,n_aligned_depths,length(batch_vars));
for curr_animal = 1:length(batch_vars)
    curr_animal_k_px = nan(437,416,43,n_aligned_depths,length(days));
    for curr_day = 1:length(batch_vars(curr_animal).k_px)        
        
        curr_k_px = batch_vars(curr_animal).k_px{curr_day};
        curr_scaled_k_px = mat2gray(bsxfun(@rdivide, ...
            curr_k_px,permute(prctile(abs( ...
            reshape(curr_k_px,[],size(curr_k_px,5))),95)',[2,3,4,5,1])),[-1,1])*2-1;     
        
        % Set any NaN explained (no MUA data probably) to NaN
        curr_scaled_k_px(:,:,:,isnan(batch_vars(curr_animal).explained_var{curr_day})) = NaN;
        
        curr_animal_k_px(:,:,:,:,curr_day) = curr_scaled_k_px;
        
    end
    k_px(:,:,:,:,curr_animal) = nanmean(curr_animal_k_px,5);
    disp(curr_animal);
end

k_px_mean = nanmean(k_px,5);
AP_image_scroll(k_px_mean);
caxis([-1,1]); axis image;
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k');


%% Fig 1b: Allen projection maps vs regression maps?
% (maybe get average centroid of str 1/2/3/4 then get one map?)
% (copy code from AP_ctx2str_probe)

probe_vector_ccf = [520,240,510;520,511,239];


%% Fig 2a: Behavior psychometric 
% (from AP_vanillaChoiceworld_behavior - currently no eliminations)

% Load behavior
bhv_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\bhv_processing\bhv.mat';
load(bhv_fn);

animals = {bhv.animal};

conditions = unique(vertcat(bhv.conditions),'rows');
trial_choice_cat = arrayfun(@(x) horzcat(bhv(x).trial_choice{:}),1:length(bhv),'uni',false);
trial_outcome_cat = arrayfun(@(x) horzcat(bhv(x).trial_outcome{:}),1:length(bhv),'uni',false);
trial_side_cat = arrayfun(@(x) horzcat(bhv(x).trial_side{:}),1:length(bhv),'uni',false);
trial_contrast_cat = arrayfun(@(x) horzcat(bhv(x).trial_contrast{:}),1:length(bhv),'uni',false);
trial_condition_cat = cellfun(@(side,contrast) side.*contrast,trial_side_cat,trial_contrast_cat,'uni',false);
trial_wheel_velocity_cat = arrayfun(@(x) vertcat(bhv(x).trial_wheel_velocity{:})',1:length(bhv),'uni',false);
stim_to_move_cat = arrayfun(@(x) horzcat(bhv(x).stim_to_move{:}),1:length(bhv),'uni',false);
stim_to_feedback_cat = arrayfun(@(x) horzcat(bhv(x).stim_to_feedback{:}),1:length(bhv),'uni',false);

% Distinguish early/late movements
go_time = 0.5;
trial_timing = arrayfun(@(animal) cellfun(@(x) 1+(x > go_time), ...
    bhv(animal).stim_to_move,'uni',false),1:length(bhv),'uni',false);
trial_timing_cat = arrayfun(@(animal) ...
    horzcat(trial_timing{animal}{:}),1:length(bhv),'uni',false);

% Plot psychometric 
frac_left = cell2mat(cellfun(@(choice,condition) ...
    grpstats(choice == -1,condition),trial_choice_cat,trial_condition_cat,'uni',false));

frac_left_earlymove = cell2mat(cellfun(@(choice,condition,stim_to_move) ...
    grpstats(choice(stim_to_move < 0.5) == -1,condition(stim_to_move < 0.5)), ...
    trial_choice_cat,trial_condition_cat,stim_to_move_cat,'uni',false));

frac_left_latemove = cell2mat(cellfun(@(choice,condition,stim_to_move) ...
    grpstats(choice(stim_to_move >= 0.5) == -1,condition(stim_to_move >= 0.5)), ...
    trial_choice_cat,trial_condition_cat,stim_to_move_cat,'uni',false));

figure;

subplot(1,3,1); hold on; axis square;
plot(conditions,frac_left,'linewidth',2,'color',[0.5,0.5,0.5]);
plot(conditions,nanmean(frac_left,2),'linewidth',5,'color','k');
xlim([-1,1]);
ylim([0,1]);
line([0,0],ylim,'linestyle','--','color','k');
line(xlim,[0.5,0.5],'linestyle','--','color','k');
xlabel('Condition');
ylabel('Fraction go left');
title('All trials');

subplot(1,3,2); hold on; axis square;
plot(conditions,frac_left_earlymove,'linewidth',2,'color',[0.5,0.5,0.5]);
plot(conditions,nanmean(frac_left_earlymove,2),'linewidth',5,'color','k');
xlim([-1,1]);
ylim([0,1]);
line([0,0],ylim,'linestyle','--','color','k');
line(xlim,[0.5,0.5],'linestyle','--','color','k');
xlabel('Condition');
ylabel('Fraction go left');
title('Early move');

subplot(1,3,3); hold on; axis square;
plot(conditions,frac_left_latemove,'linewidth',2,'color',[0.5,0.5,0.5]);
plot(conditions,nanmean(frac_left_latemove,2),'linewidth',5,'color','k');
xlim([-1,1]);
ylim([0,1]);
line([0,0],ylim,'linestyle','--','color','k');
line(xlim,[0.5,0.5],'linestyle','--','color','k');
xlabel('Condition');
ylabel('Fraction go left');
title('Late move');

%% Fig 2b/c: Cortical/striatal activity during task

% Load data
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
data_fn = ['all_trial_activity_Udf_kernel-str_4_depths_long'];

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Load use experiments and cut out bad ones
exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
exclude_fn{1} = 'bhv_use_experiments';
% exclude_fn{2} = 'expl_var_use_experiments';
use_experiments_all = {};
for curr_exclude = 1:length(exclude_fn)
    curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
    use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
end
use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
    1:size(use_experiments_all,2),'uni',false);

D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);
reward_all = cellfun(@(data,use_expts) data(use_expts),reward_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);
reward_all = reward_all(use_animals);

% Get number of widefield components and MUA depths
n_vs = size(fluor_all{1}{1},3);
n_depths = size(mua_all{1}{1},3);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_vs,1);
mua_allcat = nan(0,length(t),n_depths,1);
wheel_allcat = nan(0,length(t),1);
reward_allcat = nan(0,length(t),1);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

D_cat = struct('stimulus',cell(1,1),'response',cell(1,1),'repeatNum',cell(1,1));

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
        
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths),'uni',false));
    
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x,[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
 
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
%     wheel_cat_norm = wheel_cat./prctile(max(max(abs(wheel_cat),[],2),[],3),95);
    
    % Concatenate behavioural data
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
    D.response = cell2mat(cellfun(@(x) x.response,D_all{curr_animal},'uni',false));
    D.repeatNum = cell2mat(cellfun(@(x) x.repeatNum,D_all{curr_animal},'uni',false));
    
    D.day = trial_day;
    
    D_cat.stimulus = vertcat(D_cat.stimulus,D.stimulus);
    D_cat.response = vertcat(D_cat.response,D.response);
    D_cat.repeatNum = vertcat(D_cat.repeatNum,D.repeatNum);
    
    % Get trial ID   
    trial_contrast = max(D.stimulus,[],2);
    [~,side_idx] = max(D.stimulus > 0,[],2);
    trial_side = (side_idx-1.5)*2;
    trial_choice = -(D.response-1.5)*2;
    
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    reward_allcat = [reward_allcat;vertcat(reward_all{curr_animal}{:})];
    
    trial_contrast_allcat = [trial_contrast_allcat;trial_contrast];
    trial_side_allcat = [trial_side_allcat;trial_side];
    trial_choice_allcat = [trial_choice_allcat;trial_choice];
   
end

% Get max velocity in chosen direction
wheel_velocity_allcat = [zeros(size(wheel_allcat,1),1),diff(wheel_allcat,[],2)];
[max_speed,max_vel_idx] = max(abs(wheel_velocity_allcat(:,:,1).* ...
    (bsxfun(@times,wheel_velocity_allcat(:,:,1),trial_choice_allcat) > 0)),[],2);
vel_t = t(max_vel_idx);
max_vel = max_speed.*trial_choice_allcat;

% Get reaction time
[~,move_idx] = max(abs(wheel_allcat(:,:,1)) > 2,[],2);
move_t = t(move_idx)';

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Reaction times to plot
rxn_time_bins = {[0,0.5]};

move_align = false;
normalize_px = true;

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
    size(fluor_allcat,2)-1,size(trial_types,2),length(rxn_time_bins));
for curr_rxn = 1:length(rxn_time_bins)
    for curr_trial_type = 1:size(trial_types,2)
        
        curr_trials = trial_types(:,curr_trial_type) & ...
            move_t > rxn_time_bins{curr_rxn}(1) & ...
            move_t < rxn_time_bins{curr_rxn}(2);
        
        % Fluorescence derivative
        curr_data = fluor_allcat(curr_trials,:,:);
        
%         % Smoothed Fluorescence derivative
%         smooth_factor = 3;
%         curr_data = convn(fluor_allcat(curr_trials,:,:), ...
%             ones(1,smooth_factor)/smooth_factor,'same');
        
        if move_align
            % Re-align to movement onset
            t_leeway = 0.5;
            leeway_samples = round(t_leeway*sample_rate);
            curr_move_idx = move_idx(curr_trials);
            for i = 1:size(curr_data,1)
                curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
            end
        end
        
        curr_data_mean = squeeze(nanmean(curr_data,1))';       
        curr_px = diff(svdFrameReconstruct(U_master(:,:,1:n_vs),curr_data_mean),[],3);
%         curr_px(curr_px < 0) = 0;
        
        px_trial_types(:,:,:,curr_trial_type,curr_rxn) = curr_px;
    end
end

t_diff = conv(t,[1,1]/2,'valid');

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

% Plot
plot_rxn = 1;
AP_image_scroll(cat(4,px_combined(:,:,:,:,plot_rxn),px_combined_hemidiff(:,:,:,:,plot_rxn)),t_diff);
axis image; caxis([-1,1]); 
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k');

% Plot all concatenated
px_dims = size(px_combined);
AP_image_scroll([reshape(permute(px_combined,[1,4,2,5,3]), ...
    [],size(px_combined,2)*size(px_combined,5),size(px_combined,3)), ...
    reshape(permute(px_combined_hemidiff,[1,4,2,5,3]), ...
    [],size(px_combined_hemidiff,2)*size(px_combined_hemidiff,5), ...
    size(px_combined_hemidiff,3))],t_diff);
axis image;
caxis([-1,1]);
colormap(brewermap([],'*RdBu'))

%%% Get and plot MUA 
mua_trial_types = nan(n_depths, ...
    size(mua_allcat,2),size(trial_types,2),length(rxn_time_bins));
for curr_rxn = 1:length(rxn_time_bins)
    for curr_trial_type = 1:size(trial_types,2)
        
        curr_trials = trial_types(:,curr_trial_type) & ...
            move_t > rxn_time_bins{curr_rxn}(1) & ...
            move_t < rxn_time_bins{curr_rxn}(2);
        
        curr_data = mua_allcat(curr_trials,:,:);
        
        if move_align
            % Re-align to movement onset
            t_leeway = 0.5;
            leeway_samples = round(t_leeway*sample_rate);
            curr_move_idx = move_idx(curr_trials);
            for i = 1:size(curr_data,1)
                curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
            end
        end
        
        mua_trial_types(:,:,curr_trial_type,curr_rxn) = squeeze(nanmean(curr_data,1))';              
    end
end

% mua_trial_types_pseudodiff = cat(3,mua_trial_types(:,:,[1,3,5],:), ...
%     mua_trial_types(:,:,[1,3,5],:) -  mua_trial_types(:,:,[2,4,6],:));
mua_trial_types_pseudodiff = cat(3,mua_trial_types(:,:,[1,3,5],:), ...
    mua_trial_types(:,:,[2,4,6],:));

figure;
for curr_rxn = 1:length(rxn_time_bins)
    for curr_trial_type = 1:size(trial_types,2)/2
        
        subplot(size(trial_types,2)/2,length(rxn_time_bins)*2, ...
            curr_rxn + length(rxn_time_bins)*2*(curr_trial_type-1)); 
        hold on; set(gca,'ColorOrder',copper(n_depths));
        plot(t, ...
            mua_trial_types_pseudodiff(:,:,curr_trial_type,curr_rxn)', ...
            'linewidth',2);
        ylim([min(mua_trial_types(:)),max(mua_trial_types(:))]);
        line([0,0],ylim,'color','k')
        
        subplot(size(trial_types,2)/2,length(rxn_time_bins)*2, ...
            curr_rxn + length(rxn_time_bins) + length(rxn_time_bins)*2*(curr_trial_type-1)); 
        hold on; set(gca,'ColorOrder',copper(n_depths));
        plot(t, ...
            mua_trial_types_pseudodiff(:,:,curr_trial_type+3,curr_rxn)', ...
            'linewidth',2);
       ylim([min(mua_trial_types(:)),max(mua_trial_types(:))]);
        line([0,0],ylim,'color','k')
        
    end
end

%%% Get and plot wheel 
wheel_trial_types = nan(...
    size(wheel_velocity_allcat,2),size(trial_types,2),length(rxn_time_bins));
for curr_rxn = 1:length(rxn_time_bins)
    for curr_trial_type = 1:size(trial_types,2)
        
        curr_trials = trial_types(:,curr_trial_type) & ...
            move_t > rxn_time_bins{curr_rxn}(1) & ...
            move_t < rxn_time_bins{curr_rxn}(2);
        
        curr_data = wheel_velocity_allcat(curr_trials,:);
        if move_align
            % Re-align to movement onset
            t_leeway = 0.5;
            leeway_samples = round(t_leeway*sample_rate);
            curr_move_idx = move_idx(curr_trials);
            for i = 1:size(curr_data,1)
                curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
            end
        end
        wheel_trial_types(:,curr_trial_type,curr_rxn) = squeeze(nanmean(curr_data,1))';
               
    end
end

figure;
for curr_rxn = 1:length(rxn_time_bins)
    subplot(1,length(rxn_time_bins),curr_rxn); hold on;
    plot(t,wheel_trial_types(:,:,curr_rxn),'linewidth',2);
    ylim([-max(abs(wheel_trial_types(:))),max(abs(wheel_trial_types(:)))]);
    line([0,0],ylim,'color','k')
end

%% Fig 3b/c: regression from task events to cortex/striatum
% (clean this up, save everything necessary after regression)

%%% Load and prepare regression

% Load regression results
regression_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld\all_trial_activity_regressed.mat';
load(regression_fn);

% Get number of V's/depths
n_vs = size(fluor_allcat_predicted,3);
n_depths = size(mua_allcat_predicted,3);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% Get move onset index
[~,move_idx] = max(abs(wheel(:,:,1)) > 2,[],2);
move_t = t_downsample_diff(move_idx)';

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

% Set which components to use
use_svs = 200;

% Get time shifts in samples
sample_shifts = cellfun(@(x) round(x(1)*(sample_rate/downsample_factor)): ...
    round(x(2)*(sample_rate/downsample_factor)),t_shifts,'uni',false);

%%% Plot kernels

% Plot fluorescence kernels
for curr_regressor = 1:length(regressors)
    curr_k_v = permute(cell2mat(permute(fluor_kernel(curr_regressor,1:use_svs),[1,3,2])),[3,2,1]);
    curr_k_px = reshape(svdFrameReconstruct(U_master(:,:,1:use_svs), ...
        reshape(curr_k_v,use_svs,[])),size(U_master,1),size(U_master,2),[],size(curr_k_v,3));
    AP_image_scroll(curr_k_px,sample_shifts{curr_regressor}/(sample_rate/downsample_factor));
    axis image
    caxis([-prctile(abs(curr_k_px(:)),100),prctile(abs(curr_k_px(:)),100)]);
    colormap(brewermap([],'*RdBu'))
    AP_reference_outline('ccf_aligned','k');
    set(gcf,'Name',regressor_labels{curr_regressor});
end
% Plot fluorescence constant
curr_k_v = permute(cell2mat(permute(fluor_kernel(length(regressors)+1,1:use_svs),[1,3,2])),[3,2,1]);
curr_k_px = reshape(svdFrameReconstruct(U_master(:,:,1:use_svs), ...
    reshape(curr_k_v,use_svs,[])),size(U_master,1),size(U_master,2),[],size(curr_k_v,3));
figure;imagesc(curr_k_px);
axis image off
caxis([-prctile(abs(curr_k_px(:)),100),prctile(abs(curr_k_px(:)),100)]);
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k');
title('Constant');

% Plot fluorescence ROI kernels
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);
n_rois = numel(wf_roi);

roi_mask = cat(3,wf_roi.mask);

U_roi = bsxfun(@rdivide,transpose(reshape(U_master,[],size(U_master,3))'* ...
    reshape(roi_mask,[],size(roi_mask,3))), ...
    sum(reshape(roi_mask,[],size(roi_mask,3)),1)');

for curr_regressor = 1:length(regressors)
    curr_k_v = permute(cell2mat(permute(fluor_kernel(curr_regressor,1:use_svs),[1,3,2])),[3,2,1]);
    
    curr_k_roi = reshape(AP_svd_roi(U_master(:,:,1:n_vs), ...
        reshape(curr_k_v,n_vs,[]),[],[],roi_mask)', ...
        size(curr_k_v,2),size(curr_k_v,3),n_rois);
    
    if size(curr_k_roi,2) > 1
        curr_col = colormap_BlueWhiteRed(size(curr_k_roi,2)/2);
        curr_col(size(curr_k_roi,2)/2+1,:) = [];
    else
        curr_col = 'k'
    end
    
    figure; hold on;
    for curr_subk = 1:size(curr_k_roi,2)
        AP_stackplot(squeeze(curr_k_roi(:,curr_subk,:)), ...
            sample_shifts{curr_regressor}/(sample_rate/downsample_factor), ...
            1.5e-3,false,curr_col(curr_subk,:));
    end
    title(regressor_labels{curr_regressor});
end

% Plot MUA kernels
figure;
for curr_regressor = 1:length(regressors)
    
    curr_k_cat = permute(cat(3,mua_kernel{curr_regressor,:}),[2,1,3]);
    
    if size(curr_k_cat,2) > 1
        curr_col = colormap_BlueWhiteRed(size(curr_k_cat,2)/2);
        curr_col(size(curr_k_cat,2)/2+1,:) = [];
    else
        curr_col = 'k';
    end
    
    figure; hold on;
    for curr_subk = 1:size(curr_k_cat,2)
        AP_stackplot(squeeze(curr_k_cat(:,curr_subk,:)), ...
            sample_shifts{curr_regressor}/(sample_rate/downsample_factor), ...
            1,false,curr_col(curr_subk,:));
    end
    title(regressor_labels{curr_regressor});
end

%% Fig 4a: Passive stim responses in cortex/striatum
% (TO DO: normalize this the same as during task so they're directly
% comparable)

% Load data
n_aligned_depths = 4;
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\choiceworld';
% data_fn = ['all_trial_activity_Udf_kernel-str_passive_fullscreen_4_depths'];
data_fn = ['all_trial_activity_Udf_kernel-str_passive_choiceworld_4_depths_naive'];

load([data_path filesep data_fn]);
n_animals = length(D_all);

% Get time
framerate = 35;
raster_window = [-0.5,3];
upsample_factor = 3;
sample_rate = (framerate*upsample_factor);
t = raster_window(1):1/sample_rate:raster_window(2);

% % Load use experiments and cut out bad ones
% exclude_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\experiment_exclusion';
% exclude_fn{1} = 'bhv_use_experiments';
% % exclude_fn{2} = 'expl_var_use_experiments';
% use_experiments_all = {};
% for curr_exclude = 1:length(exclude_fn)
%     curr_use_experiments = load([exclude_path filesep exclude_fn{curr_exclude}]);
%     use_experiments_all = [use_experiments_all;curr_use_experiments.use_experiments];
% end
% use_experiments = arrayfun(@(x) all(vertcat(use_experiments_all{:,x}),1), ...
%     1:size(use_experiments_all,2),'uni',false);
% 
% D_all = cellfun(@(data,use_expts) data(use_expts),D_all,use_experiments','uni',false);
% fluor_all = cellfun(@(data,use_expts) data(use_expts),fluor_all,use_experiments','uni',false);
% mua_all = cellfun(@(data,use_expts) data(use_expts),mua_all,use_experiments','uni',false);
% wheel_all = cellfun(@(data,use_expts) data(use_expts),wheel_all,use_experiments','uni',false);

use_animals = cellfun(@(x) ~isempty(x),D_all);
D_all = D_all(use_animals);
fluor_all = fluor_all(use_animals);
mua_all = mua_all(use_animals);
wheel_all = wheel_all(use_animals);

% Get widefield ROIs
n_vs = 200;

% MUA depths
n_depths = size(mua_all{1}{1},3);

% Set up conditions
contrasts = [0,0.06,0.125,0.25,0.5,1];
sides = [-1,1];
choices = [-1,1];
timings = [1,2];
conditions = combvec(contrasts,sides,choices,timings)';
n_conditions = size(conditions,1);

% Normalize and concatenate activity/wheel
fluor_allcat = nan(0,length(t),n_vs,1);
mua_allcat = nan(0,length(t),n_depths,1);
wheel_allcat = nan(0,length(t),1);
reward_allcat = nan(0,length(t),1);

trial_contrast_allcat = [];
trial_side_allcat = [];
trial_choice_allcat = [];

D_cat = struct('stimulus',cell(1,1),'response',cell(1,1),'repeatNum',cell(1,1));

for curr_animal = 1:length(D_all)
    
    trial_day = cell2mat(cellfun(@(day,act) repmat(day,size(act,1),1), ...
        num2cell(1:length(fluor_all{curr_animal}))',fluor_all{curr_animal},'uni',false));
    
    % Concatenate fluorescence, get L-R, normalize (all together)
    fluor_cat = cat(1,fluor_all{curr_animal}{:});
    
    % Concatenate MUA and normalize (separately by day)
    mua_cat_raw = cat(1,mua_all{curr_animal}{:});
    
    n_align = size(fluor_cat,4);
    
    % NaN-out MUA trials with no spikes (no data collected - this should be
    % done when getting the activity I guess using some method besides
    % hist)
    filled_mua_trials = +cell2mat(cellfun(@(x) repmat(any(reshape(permute(x,[1,2,4,3]),[],n_depths),1), ...
        size(x,1),1),mua_all{curr_animal},'uni',false));    
    filled_mua_trials(~filled_mua_trials) = NaN;
    mua_cat_raw = bsxfun(@times,mua_cat_raw,permute(filled_mua_trials,[1,3,2,4]));
    
    % Smooth MUA with replicated edges
    smooth_size = 9; % MUST BE ODD 
    gw = gausswin(smooth_size,3)';
    smWin = gw./sum(gw);
    mua_cat_raw_smoothed = padarray(convn(mua_cat_raw,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both');
    
    t_baseline = t < 0;        
    mua_day_baseline = cell2mat(cellfun(@(x) ...
        repmat(permute(nanmean(reshape(permute(x(:,t_baseline,:,1),[1,2,4,3]),[],n_depths),1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    t_std = t > -0.5 & t < 0.5;
    mua_day_std = cell2mat(cellfun(@(x) ...
        repmat(permute(nanstd(reshape(permute(x(:,t_std,:,1),[1,2,4,3]),[],n_depths),[],1), ...
        [1,3,2,4]),[size(x,1),1]),  ...
        mat2cell(mua_cat_raw_smoothed,cellfun(@(x) size(x,1), mua_all{curr_animal}),length(t),n_depths,n_align),'uni',false));
    
    softnorm = 20;
    
    mua_cat_norm = bsxfun(@rdivide,bsxfun(@minus,mua_cat_raw_smoothed,mua_day_baseline),mua_day_std+softnorm);   
    
    % Concatenate wheel
    wheel_cat = cat(1,wheel_all{curr_animal}{:});
    
    % Concatenate stim
    D = struct;
    D.stimulus = cell2mat(cellfun(@(x) x.stimulus,D_all{curr_animal},'uni',false));
    D_cat.stimulus = vertcat(D_cat.stimulus,D.stimulus);
   
    % Concatenate everything
    fluor_allcat = [fluor_allcat;fluor_cat];
    mua_allcat = [mua_allcat;mua_cat_norm];
    wheel_allcat = [wheel_allcat;wheel_cat];
    
end

%%% REMOVE MOVE TRIALS FOR PASSIVE
move_trials = any(abs(wheel_allcat(:,t >= -0.5 & t <= 2) > 2),2);
fluor_allcat(move_trials,:,:) = [];
mua_allcat(move_trials,:,:) = [];
wheel_allcat(move_trials,:,:) = [];
D_cat.stimulus(move_trials) = [];
%%%

% Load the master U
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');

%%% Get and plot fluorescence
px_stim = nan(size(U_master,1),size(U_master,2), ...
    size(fluor_allcat,2)-1,length(unique(D_cat.stimulus)));
for curr_stim = unique(D_cat.stimulus)'
    
    curr_trials = D_cat.stimulus == curr_stim;
    
    %             % Straight fluorescence
    %             curr_data = fluor_allcat(curr_trials,:,:,plot_align);
    %             curr_baseline = nanmean(fluor_allcat(curr_trials,t > -0.2 & t < 0,:,1),2);
    %             curr_data = squeeze(nanmean(bsxfun(@minus,curr_data,curr_baseline),1))';
    %             curr_px = svdFrameReconstruct(U_master(:,:,1:n_rois),curr_data(:,1:end-1));
    
    %         % Fluorescence derivative
    %         curr_data = squeeze(nanmean(fluor_allcat(curr_trials,:,:),1))';
    %         curr_px = diff(svdFrameReconstruct(U_master(:,:,1:n_rois),curr_data),[],3);
    %         curr_px(curr_px < 0) = 0;
    
    % Smoothed Fluorescence derivative
    smooth_factor = 3;
    curr_data = convn(fluor_allcat(curr_trials,:,:), ...
        ones(1,smooth_factor)/smooth_factor,'same');
    
    curr_data_mean = squeeze(nanmean(curr_data,1))';   
    curr_px = diff(svdFrameReconstruct(U_master(:,:,1:n_vs),curr_data_mean),[],3);
%     curr_px(curr_px < 0) = 0;
    
    px_stim(:,:,:,curr_stim) = curr_px;
end

t_diff = conv(t,[1,1]/2,'valid');

% Choose stim to plot
plot_stim = unique(D_cat.stimulus);
% plot_stim = [2,5,8,10];

% Plot all together
AP_image_scroll(reshape(permute(px_stim(:,:,:,plot_stim),[1,2,4,3]), ...
    size(px_stim,1),[],size(px_stim,3)),t_diff);
axis image; caxis([-1.5e-3,1.5e-3]); 
colormap(brewermap([],'*RdBu'))
AP_reference_outline('ccf_aligned','k',[], ...
    [size(px_stim,1),size(px_stim,2),1,length(plot_stim)]);

% Plot striatum
plot_cols = lines(length(plot_stim));
figure; hold on;
p = gobjects(n_depths,length(plot_stim));
for curr_plot_condition = 1:length(plot_stim)  
    curr_trials = D_cat.stimulus == plot_stim(curr_plot_condition);
    curr_data_mean = squeeze(nanmean(mua_allcat(curr_trials,:,:),1));   
    curr_col = plot_cols(curr_plot_condition,:);   
    p(:,curr_plot_condition) = AP_stackplot(curr_data_mean,t,1,false,curr_col,1:n_depths);    
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');
legend(p(1,:),cellfun(@num2str,num2cell(unique(D_cat.stimulus)),'uni',false));

% Get fluorescence in widefield ROIs
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
wf_roi = wf_roi(:,1);
n_rois = numel(wf_roi);

roi_mask = cat(3,wf_roi.mask);

U_roi = bsxfun(@rdivide,transpose(reshape(U_master,[],size(U_master,3))'* ...
    reshape(roi_mask,[],size(roi_mask,3))), ...
    sum(reshape(roi_mask,[],size(roi_mask,3)),1)');

smooth_size = 1;
gw = gausswin(smooth_size,3)';
smWin = gw./sum(gw);

fluor_roi = permute(reshape(transpose(U_roi(:,1:n_vs)* ...
    reshape(permute(fluor_allcat(:,:,:,1),[3,2,1,4]),n_vs,[])), ...
    size(fluor_allcat,2),size(fluor_allcat,1),n_rois),[2,1,3]);

fluor_roi_diff = diff(padarray(convn(fluor_roi,smWin,'valid'),[0,floor(size(smWin,2)/2)],'replicate','both'),[],2);

% (zero negatives, normalize)
fluor_roi_diff(fluor_roi_diff < 0) = 0;
fluor_roi_diff = bsxfun(@rdivide,fluor_roi_diff,nanstd(reshape(fluor_roi_diff,[],1,size(wf_roi,1)),[],1));

t_diff =  conv(t,[1,1]/2,'valid');

figure; hold on;
plot_cols = lines(length(plot_stim));
for curr_plot_condition = 1:length(plot_stim)
    
    curr_trials = D_cat.stimulus == plot_stim(curr_plot_condition);
    curr_data = fluor_roi_diff(curr_trials,:,:);
   
    curr_data_mean = squeeze(nanmean(curr_data,1));

    curr_col = plot_cols(curr_plot_condition,:);   
    AP_stackplot(curr_data_mean,t_diff,3.5,false,curr_col,{wf_roi.area});
    
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');
title('Widefield ROI');












