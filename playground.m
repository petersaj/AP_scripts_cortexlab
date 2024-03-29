
%% Hemidiff in V (V-V_mirror) - keep for future reference
% (PUT THIS INTO A FUNCTION)

U_r = reshape(U_master(:,:,1:n_vs),[],n_vs);
U_mirror_r = reshape(AP_reflect_widefield(U_master(:,:,1:n_vs)),[],n_vs);
mirror_matrix = U_r'*U_mirror_r;
fluor_allcat_deconv_mirror = reshape(transpose( ...
    mirror_matrix*reshape(fluor_allcat_deconv,[],n_vs)'),size(fluor_allcat_deconv));

% sanity check: plot regular and mirrored
use_trials = trial_contrast_allcat > 0 & trial_side_allcat == 1 & ...
    trial_choice_allcat == -1 & move_t > 0.1 & move_t < 0.3;

V_mean = squeeze(nanmean(fluor_allcat(use_trials,:,:),1))';
V_mirror_mean = squeeze(nanmean(fluor_allcat_mirror(use_trials,:,:),1))';

px_mean = svdFrameReconstruct(U_master(:,:,1:n_vs),V_mean);
px_mirror_mean = svdFrameReconstruct(U_master(:,:,1:n_vs),V_mirror_mean);

AP_imscroll([px_mean,px_mirror_mean],t);
caxis([-prctile(abs(px_mean(:)),99),prctile(abs(px_mean(:)),99)]);
colormap(brewermap([],'*RdBu'));
axis image;

V_hemidiff = V_mean - V_mirror_mean;
px_hemidiff = svdFrameReconstruct(U_master(:,:,1:n_vs),V_hemidiff);
AP_imscroll(px_hemidiff,t);
caxis([-prctile(abs(px_hemidiff(:)),99),prctile(abs(px_hemidiff(:)),99)]);
colormap(brewermap([],'*RdBu'));
axis image;

% sanity check
use_trials = trial_contrast_allcat > 0 & trial_side_allcat == 1 & ...
    trial_choice_allcat == -1 & move_t > 0.1 & move_t < 0.3;

V_check = squeeze(nanmean(fluor_allcat_downsamp_filt(use_trials,:,:),1))';
px_check = svdFrameReconstruct(U_master(:,:,1:n_vs),V_check);

AP_imscroll(px_check,t);
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


%% Plot template waveform across channels

curr_template = 9;

figure; axes('YDir','reverse','visible','off'); hold on;
p = arrayfun(@(x) plot(0,0,'k','linewidth',2),1:size(templates,3));

template_xscale = 3;
template_yscale = 0.5;

template_y = permute(mean(templates(curr_template,:,:),1),[3,2,1]);
template_y = -template_y*template_yscale + channel_positions(:,2);
template_x = (1:size(templates,2)) + channel_positions(:,1)*template_xscale;

template_channel_amp = squeeze(range(templates(curr_template,:,:),2));
template_thresh = max(template_channel_amp)*0.2;
template_use_channels = template_channel_amp > template_thresh;
[~,max_channel] = max(max(abs(templates(curr_template,:,:)),[],2),[],3);

arrayfun(@(ch) set(p(ch),'XData',template_x(ch,:),'YData',template_y(ch,:)),1:size(templates,3));
arrayfun(@(ch) set(p(ch),'Color','r'),find(template_use_channels));
arrayfun(@(ch) set(p(ch),'Color','k'),find(~template_use_channels));
set(p(max_channel),'Color','b');

yrange = range(channel_positions(:,2))*0.03.*[-1,1];
ylim(channel_positions(max_channel,2) + yrange);


title(curr_template);



%% Waveform grab test

animal = 'AP024';
day = '2017-09-29';
experiment = 1;

[kilosort_path,kilsort_exists] = AP_cortexlab_filename(animal,day,[],'ephys');
[ephys_ap_filename,ephys_ap_exists] = AP_cortexlab_filename(animal,day,[],'ephys_ap');
[ephys_path,ephys_exists] = AP_cortexlab_filename(animal,day,[],'ephys_dir');

spike_times = readNPY([kilosort_path filesep 'spike_times.npy']);
spike_clusters = readNPY([kilosort_path filesep 'spike_clusters.npy']);

templates_whitened = readNPY([kilosort_path filesep 'templates.npy']);
channel_positions = readNPY([kilosort_path filesep 'channel_positions.npy']);
channel_map = readNPY([kilosort_path filesep 'channel_map.npy']);
winv = readNPY([kilosort_path filesep 'whitening_mat_inv.npy']);

% Default channel map/positions are from end: make from surface
channel_positions(:,2) = max(channel_positions(:,2)) - channel_positions(:,2);

% Unwhiten templates
templates = zeros(size(templates_whitened));
for t = 1:size(templates_whitened,1)
    templates(t,:,:) = squeeze(templates_whitened(t,:,:))*winv;
end

gwfparams.dataDir = kilosort_path;
gwfparams.fileName = ephys_ap_filename;
gwfparams.dataType = 'int16';
gwfparams.nCh = 384;
gwfparams.wfWin = [-40 80];
gwfparams.nWf = 200;
gwfparams.spikeTimes = spike_times;
gwfparams.spikeClusters = spike_clusters;

wf = getWaveForms(gwfparams);

cluster_mean_waveforms = permute(wf.waveFormsMean,[1,3,2]);

% wf_car_mean = permute(nanmean(wf.waveForms - nanmedian(wf.waveForms,2),2),[1,4,3,2]);
wf_car_mean = permute(wf.waveFormsMean,[1,3,2]);

curr_template = 1;



figure; hold on; axis off;
p = arrayfun(@(x) plot(0,0,'k','linewidth',2),1:size(wf_car_mean,3));

yscale = 40;
xscale = 7;

y = permute(wf_car_mean(curr_template,:,:),[3,2,1]);
y = y - channel_positions(:,2)*yscale;
x = (1:size(wf_car_mean,2)) + channel_positions(:,1)*xscale;

template_channel_amp = squeeze(range(wf_car_mean(curr_template,:,:),2));
template_thresh = max(template_channel_amp)*0.2;
template_use_channels = template_channel_amp > template_thresh;
[~,max_channel] = max(max(abs(wf_car_mean(curr_template,:,:)),[],2),[],3);

arrayfun(@(ch) set(p(ch),'XData',x(ch,:),'YData',y(ch,:)),1:size(wf_car_mean,3));
arrayfun(@(ch) set(p(ch),'Color','r'),find(template_use_channels));
arrayfun(@(ch) set(p(ch),'Color','k'),find(~template_use_channels));
set(p(max_channel),'Color','b');

yrange = range(channel_positions(:,2))*yscale*0.08.*[-1,1];
ylim([-channel_positions(max_channel,2)*yscale + yrange]);

title(curr_template);








%% Plot mean waveforms from save_phy

spike_clusters = readNPY(['C:\data_temp\phy\spike_clusters.npy']);
unique_clusters = unique(spike_clusters);
channel_positions = readNPY(['C:\data_temp\phy\channel_positions.npy']);

cluster_mean_waveforms_norm = ...
    (cluster_mean_waveforms-cluster_mean_waveforms(:,1,:))./ ...
    abs(max(max(cluster_mean_waveforms,[],2),[],3));

figure; hold on;
p = arrayfun(@(x) plot(0,0,'k'),1:size(cluster_mean_waveforms_norm,3));
for curr_template = 1:size(cluster_mean_waveforms_norm,1)
    
    y = permute(cluster_mean_waveforms_norm(curr_template,:,:),[3,2,1]);
    y = y + (max(channel_positions(:,2)) - channel_positions(:,2))/100;   
    x = (1:size(cluster_mean_waveforms_norm,2)) + channel_positions(:,1)*8;
    
    nonzero_channels = squeeze(any(cluster_mean_waveforms_norm(curr_template,:,:),2));
    [~,max_channel] = max(max(abs(cluster_mean_waveforms_norm(curr_template,:,:)),[],2),[],3);
    
    arrayfun(@(ch) set(p(ch),'XData',x(ch,:),'YData',y(ch,:)),1:size(cluster_mean_waveforms_norm,3));
    arrayfun(@(ch) set(p(ch),'Color','r'),find(nonzero_channels));
    arrayfun(@(ch) set(p(ch),'Color','k'),find(~nonzero_channels));
    set(p(max_channel),'Color','b');
    
%     ylim([min(reshape(y(nonzero_channels,:),[],1)), ...
%         max(reshape(y(nonzero_channels,:),[],1))]);
    
    curr_cluster = unique_clusters(curr_template);
    title(curr_cluster);
    waitforbuttonpress;
    
end

 
%% (to plot fluorescence by trial group)

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
    -ones(1,7),ones(1,5); ...
    ones(1,6),-ones(1,6)]';

% plot_conditions = ...
%     [contrasts(2:end),contrasts(2:end); ...
%     -ones(1,5),ones(1,5); ...
%     ones(1,5),-ones(1,5)]';

% plot_conditions = ...
%     [contrasts(2:end),contrasts(2:end); ...
%     -ones(1,5),ones(1,5); ...
%     -ones(1,5),ones(1,5)]';

% plot_conditions = ...
%     [1,0,1,1,0,1; ...
%     1,-1,-1,-1,-1,1; ...
%     -1,-1,-1,1,1,1]';

% plot_conditions = ...
%     [1,0,1,1; ...
%     -1,-1,1,-1; ...
%     -1,-1,-1,1]';

% plot_conditions = ...
%     [1; ...
%     1,; ...
%     -1]';

% plot_conditions = ...
%     [0.06,0.06,0.06; ...
%     1,1,-1; ...
%     -1,1,-1]';

use_rxn = move_t > 0 & move_t < 0.5;
    
[~,plot_id] = ismember( ...
    [trial_contrast_allcat,trial_side_allcat,trial_choice_allcat], ...
    plot_conditions,'rows');

spacing = nanstd(fluor_roi_deconv(:))*10;

% fluor_roi_deconv_hemidiff = ...
%     (fluor_roi_deconv(:,:,1:size(wf_roi,1))-fluor_roi_taskpred_reduced(:,:,1:size(wf_roi,1),2)) - ...
%     (fluor_roi_deconv(:,:,size(wf_roi,1)+1:end)-fluor_roi_taskpred_reduced(:,:,size(wf_roi,1)+1:end,2));

% fluor_roi_deconv_hemidiff = fluor_roi_deconv(:,:,1:size(wf_roi,1)) - fluor_roi_deconv(:,:,size(wf_roi,1)+1:end);

figure; hold on;
for curr_plot_condition = 1:size(plot_conditions,1)
    
    curr_trials = plot_id == curr_plot_condition & use_rxn;    
%     curr_data = fluor_roi_deconv_hemidiff(curr_trials,:,:); 
    curr_data = fluor_roi_deconv(curr_trials,:,:) - fluor_roi_taskpred_reduced(curr_trials,:,:,2);
%     curr_data = fluor_roi_deconv(curr_trials,:,:);
    
    % Re-align to movement onset
    t_leeway = 0.5;
    leeway_samples = round(t_leeway*sample_rate);
    curr_move_idx = move_idx(curr_trials);
    for i = 1:size(curr_data,1)
        curr_data(i,:,:) = circshift(curr_data(i,:,:),-curr_move_idx(i)+leeway_samples,2);
    end
    
    curr_data_mean = squeeze(nanmean(curr_data,1));

    curr_contrast_side = max(trial_contrast_allcat(curr_trials))*sign(max(trial_side_allcat(curr_trials)));
    contrast_side_idx = find(curr_contrast_side == contrast_side_val);
    curr_col = contrast_side_col(contrast_side_idx,:);   
    
    curr_linewidth = mean(trial_choice_allcat(curr_trials))+2;
    
    if curr_contrast_side == 0
        switch max(trial_choice_allcat(curr_trials))
            case -1
                curr_col = 'm';
            case 1
                curr_col = 'c';
        end
    end
    
    AP_stackplot(curr_data_mean,t,spacing,false,curr_col,{wf_roi.area});
    
end
line([0,0],ylim,'color','k');
xlabel('Time from stim');
title('Widefield ROI');



% (rank difference)
trial_groups = {'Stim','Move onset','Outcome'};
figure;
for curr_group = 1:length(trial_groups)
    
    curr_regressor_idx = strcmp(trial_groups{curr_group},regressor_labels);
    
    curr_fluor =  fluor_roi_deconv - fluor_roi_taskpred_reduced(:,:,:,curr_regressor_idx);
%     curr_fluor = curr_fluor(:,:,1:size(wf_roi,1)) - curr_fluor(:,:,size(wf_roi,1)+1:end);
    
    % (if movement, align to move onset)
    if strcmp(trial_groups{curr_group},'Move onset')
        t_leeway = -t(1);
        leeway_samples = round(t_leeway*(sample_rate));
        for i = 1:size(fluor_allcat_deconv,1)
            curr_fluor(i,:,:) = circshift(curr_fluor(i,:,:),-move_idx(i)+leeway_samples,2);
        end
    end
    
    % Set trials and grouping (-1,1) to use
    switch trial_groups{curr_group}
        case 'Stim'
            use_trials = move_t < 1;
            trial_group = trial_side_allcat;
            group1 = 1;
            group2 = -1;
        case 'Move onset'
            use_trials = move_t < 1 & trial_side_allcat == 1;
            trial_group = trial_choice_allcat;
            group1 = -1;
            group2 = 1;
        case 'Reward'
            use_trials = move_t < 1;
            trial_group = trial_choice_allcat == -trial_side_allcat;
            group1 = 1;
            group2 = 0;
    end
    
    act_rank = tiedrank(curr_fluor(use_trials,:,:));
    act_rank_difference = squeeze(( ...
        nanmean(act_rank(trial_group(use_trials) == group1,:,:),1) - ...
        nanmean(act_rank(trial_group(use_trials) == group2,:,:),1))./sum(use_trials));   
    
    p1 = subplot(1,length(trial_groups),curr_group); 
    hold on; set(gca,'ColorOrder',jet(size(curr_fluor,3)));
    plot(t,act_rank_difference,'linewidth',2);
    axis tight; line([0,0],ylim,'color','k');
    xlabel('Time');
    ylabel('Mean rank difference');
    title([trial_groups{curr_group}]);
   

end



% Get velocity and bins
% (to use max velocity regardless of final choice)
max_vel = AP_signed_max(wheel_velocity_allcat_move(:,t > 0 & t < 0.5),2);
% (to use summed velocity regardless of final choice)
% max_vel = sum(wheel_velocity_allcat_move(:,t > 0 & t < 0.5),2);
% (to use maximum cumulative velocity regardless of final choice)
% max_vel = AP_signed_max(cumsum(wheel_velocity_allcat_move(:,t > 0 & t < 0.5),2),2);

n_vel_bins = 3;
vel_amp_edges = prctile(abs(max_vel),linspace(0,100,n_vel_bins+1));
% vel_amp_edges = linspace(prctile(abs(max_vel),10),prctile(abs(max_vel),90),n_vel_bins);
vel_amp_edges = sort([vel_amp_edges,-vel_amp_edges]);

vel_amp_bins = discretize(max_vel,vel_amp_edges);
vel_amp_bins(vel_amp_bins == n_vel_bins+1) = NaN;



% Make move-aligned fluorescence
fluor_allcat_deconv_move = fluor_allcat_deconv;
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
for i = 1:size(fluor_allcat_deconv,1)
    fluor_allcat_deconv_move(i,:,:,:) = circshift(fluor_allcat_deconv_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
end

% Plot pixel fluorescence
plot_rxn_time = [0.5,Inf];

% plot_trials = ...
%     trial_contrast_allcat > 0 & ...
%     trial_side_allcat == 1 & ...
%     trial_choice_allcat == -1 & ...
%     move_t >= plot_rxn_time(1) & ...
%     move_t <= plot_rxn_time(2);

plot_trials = ...
    trial_stim_allcat > 0 & ...
    trial_choice_allcat == -1 & ...
    move_t >= plot_rxn_time(1) & ...
    move_t <= plot_rxn_time(2);

% plot_trials = ...
%     trial_contrast_allcat == 0 & ...
%     vel_amp_bins == 1 & ...
%     move_t >= plot_rxn_time(1) & ...
%     move_t <= plot_rxn_time(2);

% plot_trials = ...
%     trial_contrast_allcat > 0 & ...
%     trial_side_allcat == 1 & ...
%     vel_amp_bins == 1 & ...
%     move_t >= plot_rxn_time(1) & ...
%     move_t <= plot_rxn_time(2);

% plot_px = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
%     squeeze(nanmean( ...
%     fluor_allcat_deconv_move(plot_trials,:,:) - ...
%     fluor_taskpred_reduced_allcat_move(plot_trials,:,:,2),1))');

plot_px = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    squeeze(nanmean(fluor_allcat_deconv(plot_trials,:,:)))');

% plot_px = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
%     squeeze(nanmean(fluor_allcat_deconv_move(plot_trials,:,:)))');

AP_imscroll(plot_px,t);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
AP_reference_outline('ccf_aligned','k');




% Plot fluorescence mirroring and combining sides (side = ipsi/contra now)
mirror_matrix = reshape(U_master(:,:,1:n_vs),[],n_vs)'* ...
    reshape(AP_reflect_widefield(U_master(:,:,1:n_vs)),[],n_vs);
fluor_allcat_deconv_mirror = reshape(transpose( ...
    mirror_matrix*reshape(fluor_allcat_deconv,[],n_vs)'),size(fluor_allcat_deconv));
fluor_allcat_deconv_move_mirror = reshape(transpose( ...
    mirror_matrix*reshape(fluor_allcat_deconv_move,[],n_vs)'),size(fluor_allcat_deconv_move));
% fluor_allcat_deconv_move_mirror = reshape(transpose( ...
%     mirror_matrix*reshape(fluor_allcat_deconv_move - ...
%     fluor_taskpred_reduced_allcat_move(:,:,:,2), ...
%     [],n_vs)'),size(fluor_allcat_deconv_move));

plot_rxn_time = [0,0.5];

plot_side = 1;
plot_choice = -1;
plot_contrast = trial_contrast_allcat > 0;
plot_move_t = move_t >= plot_rxn_time(1) & move_t <= plot_rxn_time(2);

curr_trials_standard = trial_side_allcat == plot_side & ...
    trial_choice_allcat == plot_choice & ...
    plot_contrast & plot_move_t;
%     vel_amp_bins == 3;

curr_trials_mirror = trial_side_allcat == -plot_side & ...
    trial_choice_allcat == -plot_choice & ...
    plot_contrast & plot_move_t;
%     vel_amp_bins == 5;

curr_px = svdFrameReconstruct(U_master(:,:,1:n_vs),squeeze(nanmean( ...
    [fluor_allcat_deconv_move(curr_trials_standard,:,:); ...
    fluor_allcat_deconv_move_mirror(curr_trials_mirror,:,:)],1))');

% curr_px = svdFrameReconstruct(U_master(:,:,1:n_vs),squeeze(nanmean( ...
%     [fluor_allcat_deconv_move(curr_trials_standard,:,:); ...
%     fluor_allcat_deconv_move_mirror(curr_trials_mirror,:,:)],1))');

% curr_px = svdFrameReconstruct(U_master(:,:,1:n_vs),squeeze(nanmean( ...
%     [fluor_allcat_deconv_move(curr_trials_standard,:,:)-fluor_taskpred_reduced_allcat_move(curr_trials_standard,:,:,2); ...
%     fluor_allcat_deconv_move_mirror(curr_trials_mirror,:,:)],1))');

AP_imscroll(curr_px,t);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k');






% (for passive)

% Get trials with movement during stim to exclude
wheel_thresh = 0.025;
quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > wheel_thresh,2);

fluor_stim = nan(length(unique(D_allcat.stimulus)),length(t),n_vs);
for curr_v = 1:n_vs
    fluor_stim(:,:,curr_v) = grpstats(fluor_allcat_deconv(quiescent_trials,:,curr_v),D_allcat.stimulus(quiescent_trials),'nanmean');
end
fluor_stim_px = cell2mat(permute(arrayfun(@(x) svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    squeeze(fluor_stim(x,:,:))'),1:length(unique(D_allcat.stimulus)),'uni',false),[1,3,4,2]));

AP_imscroll(fluor_stim_px,t);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k');


wheel_thresh = 0.025;
quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > wheel_thresh,2);
plot_trials = quiescent_trials & trial_contrastside_allcat == 1;

plot_px = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    squeeze(nanmean(fluor_allcat_deconv(plot_trials,:,:)))');

AP_imscroll(plot_px,t);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
AP_reference_outline('ccf_aligned','k');





%% Cortex -> Striatum regression around select events

clear all
disp('Choiceworld trial activity')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);
mua_taskpred_k_all = cell(1,1);
mua_taskpred_all = cell(1,1);
mua_taskpred_reduced_all = cell(1,1);
mua_ctxpred_taskpred_k_all = cell(1,1);
mua_ctxpred_taskpred_all = cell(1,1);
mua_ctxpred_taskpred_reduced_all = cell(1,1);
fluor_taskpred_k_all = cell(1,1);
fluor_taskpred_all = cell(1,1);
fluor_taskpred_reduced_all = cell(1,1);
wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;    
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Align (only to stim)
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf,[],2)',time_bin_centers)';
        
        %%%% TESTING: DECONV
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        %%%%
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;        
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        %%%%% TEST
        warning('OVERRIDING LAMBDA')
        lambda = 5;
        %%%%%%%
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);      
        
%         %%%% TESTING DECONV
%         [ctx_str_k,predicted_spikes_std,explained_var] = ...
%             AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
%             binned_spikes_std,kernel_frames,lambda, ...
%             regression_params.zs,regression_params.cvfold, ...
%             true,regression_params.use_constant);
        
        %%%%%%%%%%%% TESTING USING SPECIFIC TIMES FOR REGRESSION
        % Only predict 0.5 seconds around leftward movement 
        
        %%% To use stim right times
        % Get right stim samples       
        stim_R_trace = histcounts(stimOn_times( ...
            signals_events.trialSideValues(1:length(stimOn_times)) == 1),time_bins);
        
        stim_R_surround_samples = round(0.5*framerate);
        stim_R_surround = conv(stim_R_trace,ones(1,stim_R_surround_samples),'same') > 0;
        
        % Set regression discontinuities
        sample_num = 1:length(stim_R_surround);
        discontinuities = [1,diff(sample_num(stim_R_surround)) > 1];
        
         % Do cortex->striatum regression around right stim
        [ctx_str_k,~,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,stim_R_surround), ...
            binned_spikes_std(:,stim_R_surround),kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant,discontinuities);
        
        
        
%         %%% To use movement L times
%         % Get trial-aligned wheel velocity
%         event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
%             wheel_velocity,t_peri_event);
%         
%         % Get reaction times
%         [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
%         move_idx(~move_trial) = NaN;
%         move_t = nan(size(move_idx));
%         move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
%         
%         % Get leftward movement samples
%         move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
%             find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == -1));
%         move_L_trace = histcounts(move_time_L_absolute,time_bins);
%         
%         move_L_surround_samples = round(0.5*framerate);
%         move_L_surround = conv(move_L_trace,ones(1,move_L_surround_samples),'same') > 0;
%         
%         % Set regression discontinuities
%         sample_num = 1:length(move_L_surround);
%         discontinuities = [1,diff(sample_num(move_L_surround)) > 1];
%         
%         % Do cortex->striatum regression around leftward movements
%         [ctx_str_k,~,explained_var] = ...
%             AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,move_L_surround), ...
%             binned_spikes_std(:,move_L_surround),kernel_frames,lambda, ...
%             regression_params.zs,regression_params.cvfold, ...
%             true,regression_params.use_constant,discontinuities);

        
        % Apply that kernel to the whole experiment
        predicted_spikes_std = nan(size(binned_spikes_std));
        for curr_depth = 1:n_depths
            curr_conv = zeros(1,size(binned_spikes_std,2));
            for curr_v = 1:length(regression_params.use_svs)
                curr_conv = curr_conv + convn(fVdf_deconv_resample(curr_v,:), ...
                    ctx_str_k{1}(curr_v,:,curr_depth),'same');
            end
            predicted_spikes_std(curr_depth,:) = ...
                curr_conv + ctx_str_k{2}(curr_depth);
        end
    
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Shove the k's into the master to save small V's 
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
                     
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        predicted_spikes = (predicted_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,predicted_spikes',t_peri_event)./raster_sample_rate;
        
        %%%%
        
        %%%% TESTING CTX->WHEEL VELOCITY
        % Resample wheel, predict both velocity and speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Shove the k's into the master to save small V's 
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%%%
        
        %%% Wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)        
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_bins(x,:)), ...
            find(trial_outcome == 1)','uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_bins(x,:)), ...
            find(trial_outcome == -1)','uni',false))) > 0;
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);    
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING: TASK REGRESSION (LONG TRACE)
                
        %%% Regression (separate regressors to get partial explained variance)
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times)).* ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end    
        
        % Stim move regressors (one for each stim when it starts to move)       
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));            
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
                        
        end      
        
        % Stim center regressors (one for each stim when it's stopped during reward)    
        unique_contrasts = unique(contrasts(contrasts > 0));
        stim_contrasts = ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ... 
                stim_contrasts == unique_contrasts(curr_contrast));
                        
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));                  
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
                        
        end    

        % Move onset regressors (L/R)          
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move ongoing regressors (L/R)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,wheel_velocity_interp < 0) = 1;
        move_ongoing_regressors(2,wheel_velocity_interp > 0) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        if length(signals_events.interactiveOnTimes) ~= length(move_t)
            error('Different number of interactive ons and move times')
        end
        
        go_cue_regressors = zeros(2,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        outcome_regressors = zeros(2,length(time_bin_centers));
        
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate regressors, set parameters
        regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        regressor_labels = {'Stim','Move onset','Go cue','Outcome'};     
        
        t_shifts = {[0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [-0.1,0.5]; ... % go cue
            [-0.5,1]}; % outcome
        
        sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;       
        
        % Regression task -> MUA
        mua_taskpred_k = cell(length(regressors)+return_constant,n_depths);
        mua_taskpred = nan(size(event_aligned_mua));
        mua_taskpred_reduced = ...
            repmat(nan(size(event_aligned_mua)),1,1,1,length(regressors));       
        for curr_depth = 1:n_depths
                        
            baseline = nanmean(reshape(event_aligned_mua(:,t < 0,curr_depth),[],1)*raster_sample_rate);
            activity = single(binned_spikes(curr_depth,:)) - baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
                   
            [task_kernel,activity_predicted,expl_var,activity_predicted_reduced] = ...
                AP_regresskernel(regressors,activity,sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            mua_taskpred_k(:,curr_depth) = task_kernel;
            
            mua_taskpred(:,:,curr_depth) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event)./raster_sample_rate;
            
            mua_taskpred_reduced(:,:,curr_depth,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event)./raster_sample_rate,permute(1:length(regressors),[1,3,4,2]),'uni',false));       
            
        end
        
        % Regression task -> MUA-ctxpred
        mua_ctxpred_taskpred_k = cell(length(regressors)+return_constant,n_depths);
        mua_ctxpred_taskpred = nan(size(event_aligned_mua));
        mua_ctxpred_taskpred_reduced = ...
            repmat(nan(size(event_aligned_mua)),1,1,1,length(regressors));       
        for curr_depth = 1:n_depths
                      
            baseline = nanmean(reshape(event_aligned_mua_ctxpred(:,t < 0,curr_depth),[],1)*raster_sample_rate);
            activity = single(predicted_spikes(curr_depth,:))-baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
            
            % (there are NaNs because of skipped edges)
            nan_samples = isnan(activity);
            activity_predicted = nan(size(activity));
            activity_predicted_reduced = nan(size(activity,1), ...
                size(activity,2),length(regressors));
            
            [task_kernel,activity_predicted(~nan_samples), ...
                expl_var,activity_predicted_reduced(:,~nan_samples,:)] = ...
                AP_regresskernel( ...
                cellfun(@(x) x(:,~nan_samples),regressors,'uni',false), ...
                activity(~nan_samples),sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            mua_ctxpred_taskpred_k(:,curr_depth) = task_kernel;
            
            mua_ctxpred_taskpred(:,:,curr_depth) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event)./raster_sample_rate;
            
            mua_ctxpred_taskpred_reduced(:,:,curr_depth,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event)./raster_sample_rate,permute(1:length(regressors),[1,3,4,2]),'uni',false));       
            
        end
        
%         % Regression task -> (master U, deconvolved) fluor
%         n_vs = length(use_components);
%         event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
%         fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);
%         
%         fluor_taskpred_k = cell(length(regressors)+return_constant,n_vs);
%         fluor_taskpred = nan(size(event_aligned_V));
%         fluor_taskpred_reduced = ...
%             repmat(nan(size(event_aligned_V)),1,1,1,length(regressors));       
%         for curr_v = 1:n_vs
%                       
%             baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,curr_v),[],1));
%             activity = single(fVdf_deconv_resample_recast(curr_v,:))-baseline;
%             
%             % Skip if nothing in this depth
%             if ~any(activity(:))
%                 continue
%             end
%             
%             % (there are NaNs because of skipped edges)
%             nan_samples = isnan(activity);
%             activity_predicted = nan(size(activity));
%             activity_predicted_reduced = nan(size(activity,1), ...
%                 size(activity,2),length(regressors));
%             
%             [task_kernel,activity_predicted(~nan_samples), ...
%                 expl_var,activity_predicted_reduced(:,~nan_samples,:)] = ...
%                 AP_regresskernel( ...
%                 cellfun(@(x) x(:,~nan_samples),regressors,'uni',false), ...
%                 activity(~nan_samples),sample_shifts, ...
%                 lambda,zs,cvfold,return_constant,use_constant);
%             
%             fluor_taskpred_k(:,curr_v) = task_kernel;
%             
%             fluor_taskpred(:,:,curr_v) = ...
%                 interp1(time_bin_centers,activity_predicted',t_peri_event);
%             
%             fluor_taskpred_reduced(:,:,curr_v,:) = cell2mat(arrayfun(@(x) ...
%                 interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
%                 t_peri_event),permute(1:length(regressors),[1,3,4,2]),'uni',false));  
%                         
%         end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        mua_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_taskpred_k;
        mua_taskpred_all{curr_animal,1}{curr_day,1} = mua_taskpred(use_trials,:,:,:);
        mua_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_taskpred_reduced(use_trials,:,:,:);
        
        mua_ctxpred_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_k;
        mua_ctxpred_taskpred_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred(use_trials,:,:,:);
        mua_ctxpred_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_reduced(use_trials,:,:,:);
        
%         fluor_taskpred_k_all{curr_animal,1}{curr_day,1} = fluor_taskpred_k;
%         fluor_taskpred_all{curr_animal,1}{curr_day,1} = fluor_taskpred(use_trials,:,:,:);
%         fluor_taskpred_reduced_all{curr_animal,1}{curr_day,1} = fluor_taskpred_reduced(use_trials,:,:,:);
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);   
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except curr_animal ...
        n_aligned_depths regression_params animals t ...
        fluor_all ...
        mua_all ...
        mua_ctxpred_all ...
        mua_taskpred_k_all ...
        mua_taskpred_all ...
        mua_taskpred_reduced_all ...
        mua_ctxpred_taskpred_k_all ...
        mua_ctxpred_taskpred_all ...
        mua_ctxpred_taskpred_reduced_all ...
        fluor_taskpred_k_all ...
        fluor_taskpred_all ...
        fluor_taskpred_reduced_all ...
        wheel_ctxpred_all ...
        ctx_str_k_all ...
        ctx_wheel_k_all ...
        wheel_all ...
        outcome_all ...
        D_all
    
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    fluor_taskpred_k_all ...
    fluor_taskpred_all ...
    fluor_taskpred_reduced_all ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    outcome_all ...
    D_all

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_STIMRTEST'];
save([save_path filesep save_fn],'-v7.3');



%% Cortex -> Striatum forward in time only

clear all
disp('Choiceworld trial activity')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [0,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);
mua_taskpred_k_all = cell(1,1);
mua_taskpred_all = cell(1,1);
mua_taskpred_reduced_all = cell(1,1);
mua_ctxpred_taskpred_k_all = cell(1,1);
mua_ctxpred_taskpred_all = cell(1,1);
mua_ctxpred_taskpred_reduced_all = cell(1,1);
fluor_taskpred_k_all = cell(1,1);
fluor_taskpred_all = cell(1,1);
fluor_taskpred_reduced_all = cell(1,1);
wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;    
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Align (only to stim)
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf,[],2)',time_bin_centers)';
        
        %%%% TESTING: DECONV
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        %%%%
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;        
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        %%%%% TEST
        warning('OVERRIDING LAMBDA')
        lambda = 5;
        %%%%%%%
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);      
        
        %%%% TESTING DECONV
        [ctx_str_k,predicted_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Shove the k's into the master to save small V's 
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
                     
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        predicted_spikes = (predicted_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,predicted_spikes',t_peri_event)./raster_sample_rate;
        
        %%%%
        
        %%%% TESTING CTX->WHEEL VELOCITY
        % Resample wheel, predict both velocity and speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Shove the k's into the master to save small V's 
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%%%
        
        %%% Wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)        
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_bins(x,:)), ...
            find(trial_outcome == 1)','uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_bins(x,:)), ...
            find(trial_outcome == -1)','uni',false))) > 0;
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);    
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING: TASK REGRESSION (LONG TRACE)
                
        %%% Regression (separate regressors to get partial explained variance)
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times)).* ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end    
        
        % Stim move regressors (one for each stim when it starts to move)       
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));            
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
                        
        end      
        
        % Stim center regressors (one for each stim when it's stopped during reward)    
        unique_contrasts = unique(contrasts(contrasts > 0));
        stim_contrasts = ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ... 
                stim_contrasts == unique_contrasts(curr_contrast));
                        
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));                  
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
                        
        end    

        % Move onset regressors (L/R)          
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move ongoing regressors (L/R)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,wheel_velocity_interp < 0) = 1;
        move_ongoing_regressors(2,wheel_velocity_interp > 0) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        if length(signals_events.interactiveOnTimes) ~= length(move_t)
            error('Different number of interactive ons and move times')
        end
        
        go_cue_regressors = zeros(2,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        outcome_regressors = zeros(2,length(time_bin_centers));
        
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate regressors, set parameters
        regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        regressor_labels = {'Stim','Move onset','Go cue','Outcome'};     
        
        t_shifts = {[0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [-0.1,0.5]; ... % go cue
            [-0.5,1]}; % outcome
        
        sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;       
        
        % Regression task -> MUA
        mua_taskpred_k = cell(length(regressors)+return_constant,n_depths);
        mua_taskpred = nan(size(event_aligned_mua));
        mua_taskpred_reduced = ...
            repmat(nan(size(event_aligned_mua)),1,1,1,length(regressors));       
        for curr_depth = 1:n_depths
                        
            baseline = nanmean(reshape(event_aligned_mua(:,t < 0,curr_depth),[],1)*raster_sample_rate);
            activity = single(binned_spikes(curr_depth,:)) - baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
                   
            [task_kernel,activity_predicted,expl_var,activity_predicted_reduced] = ...
                AP_regresskernel(regressors,activity,sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            mua_taskpred_k(:,curr_depth) = task_kernel;
            
            mua_taskpred(:,:,curr_depth) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event)./raster_sample_rate;
            
            mua_taskpred_reduced(:,:,curr_depth,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event)./raster_sample_rate,permute(1:length(regressors),[1,3,4,2]),'uni',false));       
            
        end
        
        % Regression task -> MUA-ctxpred
        mua_ctxpred_taskpred_k = cell(length(regressors)+return_constant,n_depths);
        mua_ctxpred_taskpred = nan(size(event_aligned_mua));
        mua_ctxpred_taskpred_reduced = ...
            repmat(nan(size(event_aligned_mua)),1,1,1,length(regressors));       
        for curr_depth = 1:n_depths
                      
            baseline = nanmean(reshape(event_aligned_mua_ctxpred(:,t < 0,curr_depth),[],1)*raster_sample_rate);
            activity = single(predicted_spikes(curr_depth,:))-baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
            
            % (there are NaNs because of skipped edges)
            nan_samples = isnan(activity);
            activity_predicted = nan(size(activity));
            activity_predicted_reduced = nan(size(activity,1), ...
                size(activity,2),length(regressors));
            
            [task_kernel,activity_predicted(~nan_samples), ...
                expl_var,activity_predicted_reduced(:,~nan_samples,:)] = ...
                AP_regresskernel( ...
                cellfun(@(x) x(:,~nan_samples),regressors,'uni',false), ...
                activity(~nan_samples),sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            mua_ctxpred_taskpred_k(:,curr_depth) = task_kernel;
            
            mua_ctxpred_taskpred(:,:,curr_depth) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event)./raster_sample_rate;
            
            mua_ctxpred_taskpred_reduced(:,:,curr_depth,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event)./raster_sample_rate,permute(1:length(regressors),[1,3,4,2]),'uni',false));       
            
        end
        
%         % Regression task -> (master U, deconvolved) fluor
%         n_vs = length(use_components);
%         event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
%         fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);
%         
%         fluor_taskpred_k = cell(length(regressors)+return_constant,n_vs);
%         fluor_taskpred = nan(size(event_aligned_V));
%         fluor_taskpred_reduced = ...
%             repmat(nan(size(event_aligned_V)),1,1,1,length(regressors));       
%         for curr_v = 1:n_vs
%                       
%             baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,curr_v),[],1));
%             activity = single(fVdf_deconv_resample_recast(curr_v,:))-baseline;
%             
%             % Skip if nothing in this depth
%             if ~any(activity(:))
%                 continue
%             end
%             
%             % (there are NaNs because of skipped edges)
%             nan_samples = isnan(activity);
%             activity_predicted = nan(size(activity));
%             activity_predicted_reduced = nan(size(activity,1), ...
%                 size(activity,2),length(regressors));
%             
%             [task_kernel,activity_predicted(~nan_samples), ...
%                 expl_var,activity_predicted_reduced(:,~nan_samples,:)] = ...
%                 AP_regresskernel( ...
%                 cellfun(@(x) x(:,~nan_samples),regressors,'uni',false), ...
%                 activity(~nan_samples),sample_shifts, ...
%                 lambda,zs,cvfold,return_constant,use_constant);
%             
%             fluor_taskpred_k(:,curr_v) = task_kernel;
%             
%             fluor_taskpred(:,:,curr_v) = ...
%                 interp1(time_bin_centers,activity_predicted',t_peri_event);
%             
%             fluor_taskpred_reduced(:,:,curr_v,:) = cell2mat(arrayfun(@(x) ...
%                 interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
%                 t_peri_event),permute(1:length(regressors),[1,3,4,2]),'uni',false));  
%                         
%         end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        mua_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_taskpred_k;
        mua_taskpred_all{curr_animal,1}{curr_day,1} = mua_taskpred(use_trials,:,:,:);
        mua_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_taskpred_reduced(use_trials,:,:,:);
        
        mua_ctxpred_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_k;
        mua_ctxpred_taskpred_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred(use_trials,:,:,:);
        mua_ctxpred_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_reduced(use_trials,:,:,:);
        
%         fluor_taskpred_k_all{curr_animal,1}{curr_day,1} = fluor_taskpred_k;
%         fluor_taskpred_all{curr_animal,1}{curr_day,1} = fluor_taskpred(use_trials,:,:,:);
%         fluor_taskpred_reduced_all{curr_animal,1}{curr_day,1} = fluor_taskpred_reduced(use_trials,:,:,:);
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);   
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except curr_animal ...
        n_aligned_depths regression_params animals t ...
        fluor_all ...
        mua_all ...
        mua_ctxpred_all ...
        mua_taskpred_k_all ...
        mua_taskpred_all ...
        mua_taskpred_reduced_all ...
        mua_ctxpred_taskpred_k_all ...
        mua_ctxpred_taskpred_all ...
        mua_ctxpred_taskpred_reduced_all ...
        fluor_taskpred_k_all ...
        fluor_taskpred_all ...
        fluor_taskpred_reduced_all ...
        wheel_ctxpred_all ...
        ctx_str_k_all ...
        ctx_wheel_k_all ...
        wheel_all ...
        outcome_all ...
        D_all
    
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    fluor_taskpred_k_all ...
    fluor_taskpred_all ...
    fluor_taskpred_reduced_all ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    outcome_all ...
    D_all

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_FWDCTXTEST'];
save([save_path filesep save_fn],'-v7.3');


%% Choiceworld trial activity (TRIAL TASK REGRESSORS - EARLY, NO GO CUE)

clear all
disp('Choiceworld trial activity')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);
mua_taskpred_k_all = cell(1,1);
mua_taskpred_all = cell(1,1);
mua_taskpred_reduced_all = cell(1,1);
mua_ctxpred_taskpred_k_all = cell(1,1);
mua_ctxpred_taskpred_all = cell(1,1);
mua_ctxpred_taskpred_reduced_all = cell(1,1);
fluor_taskpred_k_all = cell(1,1);
fluor_taskpred_all = cell(1,1);
fluor_taskpred_reduced_all = cell(1,1);
wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;    
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Align (only to stim)
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf,[],2)',time_bin_centers)';
        
        %%%% TESTING: DECONV
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        %%%%
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;        
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        %%%%% TEST
        warning('OVERRIDING LAMBDA')
        lambda = 5;
        %%%%%%%
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        %         [~,predicted_spikes_std,explained_var] = ...
        %             AP_regresskernel(dfVdf_resample, ...
        %             binned_spikes_std,kernel_frames,lambda, ...
        %             regression_params.zs,regression_params.cvfold, ...
        %             false,regression_params.use_constant);
        
        %%%% TESTING DECONV
        [ctx_str_k,predicted_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Shove the k's into the master to save small V's 
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
                     
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        predicted_spikes = (predicted_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,predicted_spikes',t_peri_event)./raster_sample_rate;
        
        %%%%
        
        %%%% TESTING CTX->WHEEL VELOCITY
        % Resample wheel, predict both velocity and speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Shove the k's into the master to save small V's 
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%%%
        
        %%% Wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)        
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_bins(x,:)), ...
            find(trial_outcome == 1)','uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_bins(x,:)), ...
            find(trial_outcome == -1)','uni',false))) > 0;
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);    
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING: TASK REGRESSION (LONG TRACE)
                
        %%% Regression (separate regressors to get partial explained variance)
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times)).* ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end    
        
        % Stim move regressors (one for each stim when it starts to move)       
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));            
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
                        
        end      
        
        % Stim center regressors (one for each stim when it's stopped during reward)    
        unique_contrasts = unique(contrasts(contrasts > 0));
        stim_contrasts = ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ... 
                stim_contrasts == unique_contrasts(curr_contrast));
                        
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));                  
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
                        
        end    

        % Move onset regressors (L/R)          
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move ongoing regressors (L/R choice for duration of movement)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        
        move_stopped_t = 0.5;
        move_stopped_samples = round(sample_rate*move_stopped_t);
        wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.02), ...
            ones(1,move_stopped_samples),'full') > 0;
        wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);
        
        move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_L_absolute','uni',false));
        move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_R_absolute','uni',false));
        
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
        move_ongoing_regressors(2,move_ongoing_R_samples) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        if length(signals_events.interactiveOnTimes) ~= length(move_t)
            error('Different number of interactive ons and move times')
        end
        
        go_cue_regressors = zeros(2,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        outcome_regressors = zeros(2,length(time_bin_centers));
        
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate regressors, set parameters
        regressors = {stim_regressors;move_onset_regressors;outcome_regressors};
        regressor_labels = {'Stim','Move onset','Outcome'};     
        
        t_shifts = {[0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [-0.5,1]}; % outcome
        
        sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;
        
        % Use trial structure (to not regress all trials)
        regress_trials = move_t < 0.5;
        trial_regressors = cellfun(@(regressor) ...
            reshape(interp1(time_bin_centers,regressor', ...
            t_peri_event(regress_trials,:)','nearest',0),[],size(regressor,1))', ...
            regressors,'uni',false);        
        discontinuities = reshape([ones(sum(regress_trials),1), ...
            zeros(sum(regress_trials),length(t)-1)]',[],1);
        
        % Regression task -> MUA
        mua_taskpred_k = cell(length(regressors)+return_constant,n_depths);
        mua_taskpred = nan(size(event_aligned_mua));
        mua_taskpred_reduced = ...
            repmat(nan(size(event_aligned_mua)),1,1,1,length(regressors));       
        for curr_depth = 1:n_depths
                        
            baseline = nanmean(reshape(event_aligned_mua(:,t < 0,curr_depth),[],1)*raster_sample_rate);
            %activity = single(binned_spikes(curr_depth,:)) - baseline;           
            activity = reshape(event_aligned_mua(regress_trials,:,curr_depth)'* ...
                raster_sample_rate,1,[]) - baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
                   
            [task_kernel,activity_predicted,expl_var,activity_predicted_reduced] = ...
                AP_regresskernel(trial_regressors,activity,sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant,discontinuities);
            
            mua_taskpred_k(:,curr_depth) = task_kernel;
            
            mua_taskpred(regress_trials,:,curr_depth) = ...
                reshape(activity_predicted,length(t),sum(regress_trials))'./raster_sample_rate;
            
            mua_taskpred_reduced(regress_trials,:,curr_depth,:) = ...
                cell2mat(arrayfun(@(x) reshape(activity_predicted_reduced(:,:,x), ...
                length(t),sum(regress_trials))'./raster_sample_rate, ...
                permute(1:length(regressors),[1,3,4,2]),'uni',false));                
            
        end
        
        % Regression task -> MUA-ctxpred
        mua_ctxpred_taskpred_k = cell(length(regressors)+return_constant,n_depths);
        mua_ctxpred_taskpred = nan(size(event_aligned_mua));
        mua_ctxpred_taskpred_reduced = ...
            repmat(nan(size(event_aligned_mua)),1,1,1,length(regressors));       
        for curr_depth = 1:n_depths
                      
            baseline = nanmean(reshape(event_aligned_mua_ctxpred(:,t < 0,curr_depth),[],1)*raster_sample_rate);
            %activity = single(predicted_spikes(curr_depth,:))-baseline;
            activity = reshape(event_aligned_mua_ctxpred(regress_trials,:,curr_depth)'* ...
                raster_sample_rate,1,[]) - baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
            
            % (there are NaNs because of skipped edges)
            nan_samples = isnan(activity);
            activity_predicted = nan(size(activity));
            activity_predicted_reduced = nan(size(activity,1), ...
                size(activity,2),length(regressors));
            
            [task_kernel,activity_predicted(~nan_samples), ...
                expl_var,activity_predicted_reduced(:,~nan_samples,:)] = ...
                AP_regresskernel( ...
                cellfun(@(x) x(:,~nan_samples),trial_regressors,'uni',false), ...
                activity(~nan_samples),sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            mua_ctxpred_taskpred_k(:,curr_depth) = task_kernel;
            
            mua_ctxpred_taskpred(regress_trials,:,curr_depth) = ...
                reshape(activity_predicted,length(t),sum(regress_trials))'./raster_sample_rate;
            
            mua_ctxpred_taskpred_reduced(regress_trials,:,curr_depth,:) = ...
                cell2mat(arrayfun(@(x) reshape(activity_predicted_reduced(:,:,x), ...
                length(t),sum(regress_trials))'./raster_sample_rate, ...
                permute(1:length(regressors),[1,3,4,2]),'uni',false));                
            
        end
        
        % Regression task -> (master U, deconvolved) fluor
        n_vs = length(use_components);
        event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
        fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);
        
        fluor_taskpred_k = cell(length(regressors)+return_constant,n_vs);
        fluor_taskpred = nan(size(event_aligned_V));
        fluor_taskpred_reduced = ...
            repmat(nan(size(event_aligned_V)),1,1,1,length(regressors));       
        for curr_v = 1:n_vs
                      
            baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,curr_v),[],1));
            %activity = single(fVdf_deconv_resample_recast(curr_v,:))-baseline;
            activity = reshape(event_aligned_V_deconv(regress_trials,:,curr_v)',1,[]) - baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
            
            % (there are NaNs because of skipped edges)
            nan_samples = isnan(activity);
            activity_predicted = nan(size(activity));
            activity_predicted_reduced = nan(size(activity,1), ...
                size(activity,2),length(regressors));
            
            [task_kernel,activity_predicted(~nan_samples), ...
                expl_var,activity_predicted_reduced(:,~nan_samples,:)] = ...
                AP_regresskernel( ...
                cellfun(@(x) x(:,~nan_samples),trial_regressors,'uni',false), ...
                activity(~nan_samples),sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            fluor_taskpred_k(:,curr_v) = task_kernel;
            
            fluor_taskpred(regress_trials,:,curr_v) = ...
                reshape(activity_predicted,length(t),sum(regress_trials))'./raster_sample_rate;
            
            fluor_taskpred_reduced(regress_trials,:,curr_v,:) = ...
                cell2mat(arrayfun(@(x) reshape(activity_predicted_reduced(:,:,x), ...
                length(t),sum(regress_trials))'./raster_sample_rate, ...
                permute(1:length(regressors),[1,3,4,2]),'uni',false));            
                        
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        mua_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_taskpred_k;
        mua_taskpred_all{curr_animal,1}{curr_day,1} = mua_taskpred(use_trials,:,:,:);
        mua_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_taskpred_reduced(use_trials,:,:,:);
        
        mua_ctxpred_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_k;
        mua_ctxpred_taskpred_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred(use_trials,:,:,:);
        mua_ctxpred_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_reduced(use_trials,:,:,:);
        
        fluor_taskpred_k_all{curr_animal,1}{curr_day,1} = fluor_taskpred_k;
        fluor_taskpred_all{curr_animal,1}{curr_day,1} = fluor_taskpred(use_trials,:,:,:);
        fluor_taskpred_reduced_all{curr_animal,1}{curr_day,1} = fluor_taskpred_reduced(use_trials,:,:,:);
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);   
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except curr_animal ...
        n_aligned_depths regression_params animals t ...
        fluor_all ...
        mua_all ...
        mua_ctxpred_all ...
        mua_taskpred_k_all ...
        mua_taskpred_all ...
        mua_taskpred_reduced_all ...
        mua_ctxpred_taskpred_k_all ...
        mua_ctxpred_taskpred_all ...
        mua_ctxpred_taskpred_reduced_all ...
        fluor_taskpred_k_all ...
        fluor_taskpred_all ...
        fluor_taskpred_reduced_all ...
        wheel_ctxpred_all ...
        ctx_str_k_all ...
        ctx_wheel_k_all ...
        wheel_all ...
        outcome_all ...
        D_all
    
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    fluor_taskpred_k_all ...
    fluor_taskpred_all ...
    fluor_taskpred_reduced_all ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    outcome_all ...
    D_all

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_EARLYMOVETEST'];
save([save_path filesep save_fn],'-v7.3');


%% Choiceworld trial activity (TESTING ONLY REGRESSING MOVEMENT-RELATED)

clear all
disp('Choiceworld trial activity')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);
mua_taskpred_k_all = cell(1,1);
mua_taskpred_all = cell(1,1);
mua_taskpred_reduced_all = cell(1,1);
mua_ctxpred_taskpred_k_all = cell(1,1);
mua_ctxpred_taskpred_all = cell(1,1);
mua_ctxpred_taskpred_reduced_all = cell(1,1);
fluor_taskpred_k_all = cell(1,1);
fluor_taskpred_all = cell(1,1);
fluor_taskpred_reduced_all = cell(1,1);
wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;    
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Align (only to stim)
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf,[],2)',time_bin_centers)';
        
        %%%% TESTING: DECONV
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        %%%%
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;        
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        %%%%% TEST
        warning('OVERRIDING LAMBDA')
        lambda = 5;
        %%%%%%%
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        %         [~,predicted_spikes_std,explained_var] = ...
        %             AP_regresskernel(dfVdf_resample, ...
        %             binned_spikes_std,kernel_frames,lambda, ...
        %             regression_params.zs,regression_params.cvfold, ...
        %             false,regression_params.use_constant);
        
        %%%% TESTING DECONV
        [ctx_str_k,predicted_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Shove the k's into the master to save small V's 
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
                     
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        predicted_spikes = (predicted_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,predicted_spikes',t_peri_event)./raster_sample_rate;
        
        %%%%
        
        %%%% TESTING CTX->WHEEL VELOCITY
        % Resample wheel, predict both velocity and speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Shove the k's into the master to save small V's 
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%%%
        
        %%% Wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)        
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_bins(x,:)), ...
            find(trial_outcome == 1)','uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_bins(x,:)), ...
            find(trial_outcome == -1)','uni',false))) > 0;
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);    
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING: TASK REGRESSION (LONG TRACE)
                
        %%% Regression (separate regressors to get partial explained variance)
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times)).* ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end    
        
        % Stim move regressors (one for each stim when it starts to move)       
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));            
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
                        
        end      
        
        % Stim center regressors (one for each stim when it's stopped during reward)    
        unique_contrasts = unique(contrasts(contrasts > 0));
        stim_contrasts = ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ... 
                stim_contrasts == unique_contrasts(curr_contrast));
                        
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));                  
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
                        
        end    

        % Move onset regressors (L/R)          
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move ongoing regressors (L/R choice for duration of movement)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        
        move_stopped_t = 0.5;
        move_stopped_samples = round(sample_rate*move_stopped_t);
        wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.02), ...
            ones(1,move_stopped_samples),'full') > 0;
        wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);
        
        move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_L_absolute','uni',false));
        move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_R_absolute','uni',false));
        
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
        move_ongoing_regressors(2,move_ongoing_R_samples) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        if length(signals_events.interactiveOnTimes) ~= length(move_t)
            error('Different number of interactive ons and move times')
        end
        
        go_cue_regressors = zeros(2,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        outcome_regressors = zeros(2,length(time_bin_centers));
        
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate regressors, set parameters
        regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        regressor_labels = {'Stim','Move onset','Go cue','Outcome'};     
        
        t_shifts = {[0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [-0.1,0.5]; ... % go cue
            [-0.5,1]}; % outcome
        
        sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;       
        
        % Regression task -> MUA
        mua_taskpred_k = cell(length(regressors)+return_constant,n_depths);
        mua_taskpred = nan(size(event_aligned_mua));
        mua_taskpred_reduced = ...
            repmat(nan(size(event_aligned_mua)),1,1,1,length(regressors));       
        binned_spikes_moveonly = nan(size(binned_spikes));
        for curr_depth = 1:n_depths
                        
            baseline = nanmean(reshape(event_aligned_mua(:,t < 0,curr_depth),[],1)*raster_sample_rate);
            activity = single(binned_spikes(curr_depth,:)) - baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
                   
            [task_kernel,activity_predicted,expl_var,activity_predicted_reduced] = ...
                AP_regresskernel(regressors,activity,sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            mua_taskpred_k(:,curr_depth) = task_kernel;
            
            mua_taskpred(:,:,curr_depth) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event)./raster_sample_rate;
            
            mua_taskpred_reduced(:,:,curr_depth,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event)./raster_sample_rate,permute(1:length(regressors),[1,3,4,2]),'uni',false));       
            
            binned_spikes_moveonly(curr_depth,:) = activity_predicted_reduced(:,:,2);
            
        end
        
        % Regression task -> MUA-ctxpred
        mua_ctxpred_taskpred_k = cell(length(regressors)+return_constant,n_depths);
        mua_ctxpred_taskpred = nan(size(event_aligned_mua));
        mua_ctxpred_taskpred_reduced = ...
            repmat(nan(size(event_aligned_mua)),1,1,1,length(regressors));       
        for curr_depth = 1:n_depths
                      
            baseline = nanmean(reshape(event_aligned_mua_ctxpred(:,t < 0,curr_depth),[],1)*raster_sample_rate);
            activity = single(predicted_spikes(curr_depth,:))-baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
            
            % (there are NaNs because of skipped edges)
            nan_samples = isnan(activity);
            activity_predicted = nan(size(activity));
            activity_predicted_reduced = nan(size(activity,1), ...
                size(activity,2),length(regressors));
            
            [task_kernel,activity_predicted(~nan_samples), ...
                expl_var,activity_predicted_reduced(:,~nan_samples,:)] = ...
                AP_regresskernel( ...
                cellfun(@(x) x(:,~nan_samples),regressors,'uni',false), ...
                activity(~nan_samples),sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            mua_ctxpred_taskpred_k(:,curr_depth) = task_kernel;
            
            mua_ctxpred_taskpred(:,:,curr_depth) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event)./raster_sample_rate;
            
            mua_ctxpred_taskpred_reduced(:,:,curr_depth,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event)./raster_sample_rate,permute(1:length(regressors),[1,3,4,2]),'uni',false));       
            
        end
        
        % Regression task -> (master U, deconvolved) fluor
        n_vs = length(use_components);
        event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
        fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);
        
        fluor_taskpred_k = cell(length(regressors)+return_constant,n_vs);
        fluor_taskpred = nan(size(event_aligned_V));
        fluor_taskpred_reduced = ...
            repmat(nan(size(event_aligned_V)),1,1,1,length(regressors));       
        
        fVdf_moveonly = nan(n_vs,length(time_bin_centers),'single');
        for curr_v = 1:n_vs
                      
            baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,curr_v),[],1));
            activity = single(fVdf_deconv_resample_recast(curr_v,:))-baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
            
            % (there are NaNs because of skipped edges)
            nan_samples = isnan(activity);
            activity_predicted = nan(size(activity));
            activity_predicted_reduced = nan(size(activity,1), ...
                size(activity,2),length(regressors));
            
            [task_kernel,activity_predicted(~nan_samples), ...
                expl_var,activity_predicted_reduced(:,~nan_samples,:)] = ...
                AP_regresskernel( ...
                cellfun(@(x) x(:,~nan_samples),regressors,'uni',false), ...
                activity(~nan_samples),sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            fluor_taskpred_k(:,curr_v) = task_kernel;
            
            fluor_taskpred(:,:,curr_v) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event);
            
            fluor_taskpred_reduced(:,:,curr_v,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event),permute(1:length(regressors),[1,3,4,2]),'uni',false));  
            
            fVdf_moveonly(curr_v,:) = activity_predicted_reduced(:,:,2);
                        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING: Do ctx->str move-related only
        binned_spikes_moveonly(isnan(binned_spikes_moveonly)) = 0;
        fVdf_moveonly(isnan(fVdf_moveonly)) = 0;
        
        binned_spikes_move_std = (binned_spikes-binned_spikes_moveonly)./ ...
            nanstd((binned_spikes-binned_spikes_moveonly),[],2);
        binned_spikes_move_std(isnan(binned_spikes_move_std)) = 0;        
        
        fVdf_recast_move = fVdf_deconv_resample_recast(1:n_vs,:) - fVdf_moveonly;
        
        lambda = 5;      
        [ctx_str_k,predicted_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_recast_move(regression_params.use_svs,:), ...
            binned_spikes_move_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        ctx_str_k_recast = ctx_str_k{1};
                             
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        predicted_spikes = (predicted_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,predicted_spikes',t_peri_event)./raster_sample_rate;      

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        mua_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_taskpred_k;
        mua_taskpred_all{curr_animal,1}{curr_day,1} = mua_taskpred(use_trials,:,:,:);
        mua_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_taskpred_reduced(use_trials,:,:,:);
        
        mua_ctxpred_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_k;
        mua_ctxpred_taskpred_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred(use_trials,:,:,:);
        mua_ctxpred_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_reduced(use_trials,:,:,:);
        
        fluor_taskpred_k_all{curr_animal,1}{curr_day,1} = fluor_taskpred_k;
        fluor_taskpred_all{curr_animal,1}{curr_day,1} = fluor_taskpred(use_trials,:,:,:);
        fluor_taskpred_reduced_all{curr_animal,1}{curr_day,1} = fluor_taskpred_reduced(use_trials,:,:,:);
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);   
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except curr_animal ...
        n_aligned_depths regression_params animals t ...
        fluor_all ...
        mua_all ...
        mua_ctxpred_all ...
        mua_taskpred_k_all ...
        mua_taskpred_all ...
        mua_taskpred_reduced_all ...
        mua_ctxpred_taskpred_k_all ...
        mua_ctxpred_taskpred_all ...
        mua_ctxpred_taskpred_reduced_all ...
        fluor_taskpred_k_all ...
        fluor_taskpred_all ...
        fluor_taskpred_reduced_all ...
        wheel_ctxpred_all ...
        ctx_str_k_all ...
        ctx_wheel_k_all ...
        wheel_all ...
        outcome_all ...
        D_all
    
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    fluor_taskpred_k_all ...
    fluor_taskpred_all ...
    fluor_taskpred_reduced_all ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    outcome_all ...
    D_all

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_CTXSTRMOVEONLY'];
save([save_path filesep save_fn],'-v7.3');


%% Choiceworld trial activity (CTX->STR RESIDUAL 2X)

clear all
disp('Choiceworld trial activity')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);
mua_taskpred_k_all = cell(1,1);
mua_taskpred_all = cell(1,1);
mua_taskpred_reduced_all = cell(1,1);
mua_ctxpred_taskpred_k_all = cell(1,1);
mua_ctxpred_taskpred_all = cell(1,1);
mua_ctxpred_taskpred_reduced_all = cell(1,1);
fluor_taskpred_k_all = cell(1,1);
fluor_taskpred_all = cell(1,1);
fluor_taskpred_reduced_all = cell(1,1);
wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;    
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Align (only to stim)
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf,[],2)',time_bin_centers)';
        
        %%%% TESTING: DECONV
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        %%%%
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;        
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        %%%%% TEST
        warning('OVERRIDING LAMBDA')
        lambda = 5;
        %%%%%%%
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);       
        
        [ctx_str_k,predicted_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        %%%%%% TESTING 2x regression to regress the residuals
        binned_spikes_std_residuals = binned_spikes_std - predicted_spikes_std;
        binned_spikes_std_residuals(isnan(binned_spikes_std_residuals)) = 0;
        
        [ctx_str_k,predicted_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std_residuals,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Shove the k's into the master to save small V's 
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
                     
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        predicted_spikes = (predicted_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,predicted_spikes',t_peri_event)./raster_sample_rate;
        
        %%%%
        
        %%%% TESTING CTX->WHEEL VELOCITY
        % Resample wheel, predict both velocity and speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Shove the k's into the master to save small V's 
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%%%
        
        %%% Wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)        
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_bins(x,:)), ...
            find(trial_outcome == 1)','uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_bins(x,:)), ...
            find(trial_outcome == -1)','uni',false))) > 0;
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);    
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING: TASK REGRESSION (LONG TRACE)
                
        %%% Regression (separate regressors to get partial explained variance)
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times)).* ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end    
        
        % Stim move regressors (one for each stim when it starts to move)       
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));            
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
                        
        end      
        
        % Stim center regressors (one for each stim when it's stopped during reward)    
        unique_contrasts = unique(contrasts(contrasts > 0));
        stim_contrasts = ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ... 
                stim_contrasts == unique_contrasts(curr_contrast));
                        
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));                  
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
                        
        end    

        % Move onset regressors (L/R)          
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move ongoing regressors (L/R choice for duration of movement)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        
        move_stopped_t = 0.5;
        move_stopped_samples = round(sample_rate*move_stopped_t);
        wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.02), ...
            ones(1,move_stopped_samples),'full') > 0;
        wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);
        
        move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_L_absolute','uni',false));
        move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_R_absolute','uni',false));
        
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
        move_ongoing_regressors(2,move_ongoing_R_samples) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        if length(signals_events.interactiveOnTimes) ~= length(move_t)
            error('Different number of interactive ons and move times')
        end
        
        go_cue_regressors = zeros(2,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        outcome_regressors = zeros(2,length(time_bin_centers));
        
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate regressors, set parameters
        regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        regressor_labels = {'Stim','Move onset','Go cue','Outcome'};     
        
        t_shifts = {[0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [-0.1,0.5]; ... % go cue
            [-0.5,1]}; % outcome
        
        sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;       
        
        % Regression task -> MUA
        mua_taskpred_k = cell(length(regressors)+return_constant,n_depths);
        mua_taskpred = nan(size(event_aligned_mua));
        mua_taskpred_reduced = ...
            repmat(nan(size(event_aligned_mua)),1,1,1,length(regressors));       
        for curr_depth = 1:n_depths
                        
            baseline = nanmean(reshape(event_aligned_mua(:,t < 0,curr_depth),[],1)*raster_sample_rate);
            activity = single(binned_spikes(curr_depth,:)) - baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
                   
            [task_kernel,activity_predicted,expl_var,activity_predicted_reduced] = ...
                AP_regresskernel(regressors,activity,sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            mua_taskpred_k(:,curr_depth) = task_kernel;
            
            mua_taskpred(:,:,curr_depth) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event)./raster_sample_rate;
            
            mua_taskpred_reduced(:,:,curr_depth,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event)./raster_sample_rate,permute(1:length(regressors),[1,3,4,2]),'uni',false));       
            
        end
        
        % Regression task -> MUA-ctxpred
        mua_ctxpred_taskpred_k = cell(length(regressors)+return_constant,n_depths);
        mua_ctxpred_taskpred = nan(size(event_aligned_mua));
        mua_ctxpred_taskpred_reduced = ...
            repmat(nan(size(event_aligned_mua)),1,1,1,length(regressors));       
        for curr_depth = 1:n_depths
                      
            baseline = nanmean(reshape(event_aligned_mua_ctxpred(:,t < 0,curr_depth),[],1)*raster_sample_rate);
            activity = single(predicted_spikes(curr_depth,:))-baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
            
            % (there are NaNs because of skipped edges)
            nan_samples = isnan(activity);
            activity_predicted = nan(size(activity));
            activity_predicted_reduced = nan(size(activity,1), ...
                size(activity,2),length(regressors));
            
            [task_kernel,activity_predicted(~nan_samples), ...
                expl_var,activity_predicted_reduced(:,~nan_samples,:)] = ...
                AP_regresskernel( ...
                cellfun(@(x) x(:,~nan_samples),regressors,'uni',false), ...
                activity(~nan_samples),sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            mua_ctxpred_taskpred_k(:,curr_depth) = task_kernel;
            
            mua_ctxpred_taskpred(:,:,curr_depth) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event)./raster_sample_rate;
            
            mua_ctxpred_taskpred_reduced(:,:,curr_depth,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event)./raster_sample_rate,permute(1:length(regressors),[1,3,4,2]),'uni',false));       
            
        end
        
        % Regression task -> (master U, deconvolved) fluor
        n_vs = length(use_components);
        event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
        fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);
        
        fluor_taskpred_k = cell(length(regressors)+return_constant,n_vs);
        fluor_taskpred = nan(size(event_aligned_V));
        fluor_taskpred_reduced = ...
            repmat(nan(size(event_aligned_V)),1,1,1,length(regressors));       
        for curr_v = 1:n_vs
                      
            baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,curr_v),[],1));
            activity = single(fVdf_deconv_resample_recast(curr_v,:))-baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
            
            % (there are NaNs because of skipped edges)
            nan_samples = isnan(activity);
            activity_predicted = nan(size(activity));
            activity_predicted_reduced = nan(size(activity,1), ...
                size(activity,2),length(regressors));
            
            [task_kernel,activity_predicted(~nan_samples), ...
                expl_var,activity_predicted_reduced(:,~nan_samples,:)] = ...
                AP_regresskernel( ...
                cellfun(@(x) x(:,~nan_samples),regressors,'uni',false), ...
                activity(~nan_samples),sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            fluor_taskpred_k(:,curr_v) = task_kernel;
            
            fluor_taskpred(:,:,curr_v) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event);
            
            fluor_taskpred_reduced(:,:,curr_v,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event),permute(1:length(regressors),[1,3,4,2]),'uni',false));  
                        
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        mua_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_taskpred_k;
        mua_taskpred_all{curr_animal,1}{curr_day,1} = mua_taskpred(use_trials,:,:,:);
        mua_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_taskpred_reduced(use_trials,:,:,:);
        
        mua_ctxpred_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_k;
        mua_ctxpred_taskpred_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred(use_trials,:,:,:);
        mua_ctxpred_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_reduced(use_trials,:,:,:);
        
        fluor_taskpred_k_all{curr_animal,1}{curr_day,1} = fluor_taskpred_k;
        fluor_taskpred_all{curr_animal,1}{curr_day,1} = fluor_taskpred(use_trials,:,:,:);
        fluor_taskpred_reduced_all{curr_animal,1}{curr_day,1} = fluor_taskpred_reduced(use_trials,:,:,:);
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);   
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except curr_animal ...
        n_aligned_depths regression_params animals t ...
        fluor_all ...
        mua_all ...
        mua_ctxpred_all ...
        mua_taskpred_k_all ...
        mua_taskpred_all ...
        mua_taskpred_reduced_all ...
        mua_ctxpred_taskpred_k_all ...
        mua_ctxpred_taskpred_all ...
        mua_ctxpred_taskpred_reduced_all ...
        fluor_taskpred_k_all ...
        fluor_taskpred_all ...
        fluor_taskpred_reduced_all ...
        wheel_ctxpred_all ...
        ctx_str_k_all ...
        ctx_wheel_k_all ...
        wheel_all ...
        outcome_all ...
        D_all
    
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    fluor_taskpred_k_all ...
    fluor_taskpred_all ...
    fluor_taskpred_reduced_all ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    outcome_all ...
    D_all

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_REGRESSRESIDUAL'];
save([save_path filesep save_fn],'-v7.3');


%% Choiceworld trial activity (TESTING MOVExSTIM INTERACTION REGRESSION)

clear all
disp('Choiceworld trial activity')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);
mua_taskpred_k_all = cell(1,1);
mua_taskpred_all = cell(1,1);
mua_taskpred_reduced_all = cell(1,1);
mua_ctxpred_taskpred_k_all = cell(1,1);
mua_ctxpred_taskpred_all = cell(1,1);
mua_ctxpred_taskpred_reduced_all = cell(1,1);
fluor_taskpred_k_all = cell(1,1);
fluor_taskpred_all = cell(1,1);
fluor_taskpred_reduced_all = cell(1,1);
wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;    
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Align (only to stim)
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf,[],2)',time_bin_centers)';
        
        %%%% TESTING: DECONV
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        %%%%
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;        
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        %%%%% TEST
        warning('OVERRIDING LAMBDA')
        lambda = 5;
        %%%%%%%
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        %         [~,predicted_spikes_std,explained_var] = ...
        %             AP_regresskernel(dfVdf_resample, ...
        %             binned_spikes_std,kernel_frames,lambda, ...
        %             regression_params.zs,regression_params.cvfold, ...
        %             false,regression_params.use_constant);
        
        %%%% TESTING DECONV
        [ctx_str_k,predicted_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Shove the k's into the master to save small V's 
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
                     
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        predicted_spikes = (predicted_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,predicted_spikes',t_peri_event)./raster_sample_rate;
        
        %%%%
        
        %%%% TESTING CTX->WHEEL VELOCITY
        % Resample wheel, predict both velocity and speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Shove the k's into the master to save small V's 
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%%%
        
        %%% Wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)        
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_bins(x,:)), ...
            find(trial_outcome == 1)','uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_bins(x,:)), ...
            find(trial_outcome == -1)','uni',false))) > 0;
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);    
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING: TASK REGRESSION (LONG TRACE)
                
        %%% Regression (separate regressors to get partial explained variance)
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times)).* ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end    
        
        % Stim move regressors (one for each stim when it starts to move)       
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));            
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
                        
        end      
        
        % Stim center regressors (one for each stim when it's stopped during reward)    
        unique_contrasts = unique(contrasts(contrasts > 0));
        stim_contrasts = ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ... 
                stim_contrasts == unique_contrasts(curr_contrast));
                        
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));                  
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
                        
        end    

        % Move onset regressors (L/R)          
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move onset x stim regressors (one for each contrast/side)
        move_onset_stim_time_absolute = arrayfun(@(curr_stim) ...
            arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & stim_contrastsides' == unique_stim(curr_stim))), ...
            1:length(unique_stim),'uni',false);
        
        move_onset_stim_regressors = zeros(10,length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
           move_onset_stim_regressors(curr_stim,:) = ...
               histcounts(move_onset_stim_time_absolute{curr_stim},time_bins); 
        end
        
        % Move ongoing regressors (L/R choice for duration of movement)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        
        move_stopped_t = 0.5;
        move_stopped_samples = round(sample_rate*move_stopped_t);
        wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.02), ...
            ones(1,move_stopped_samples),'full') > 0;
        wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);
        
        move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_L_absolute','uni',false));
        move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_R_absolute','uni',false));
        
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
        move_ongoing_regressors(2,move_ongoing_R_samples) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        if length(signals_events.interactiveOnTimes) ~= length(move_t)
            error('Different number of interactive ons and move times')
        end
        
        go_cue_regressors = zeros(2,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        outcome_regressors = zeros(2,length(time_bin_centers));
        
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate regressors, set parameters
        regressors = {stim_regressors;move_onset_regressors;move_onset_stim_regressors;go_cue_regressors;outcome_regressors};
        regressor_labels = {'Stim','Move onset','Move onset stim','Go cue','Outcome'};     
        
        t_shifts = {[0,0.5]; ... % stim
            [-0.5,1]; ... % move onset
            [-0.5,1]; ... % move onset stim
            [-0.1,0.5]; ... % go cue
            [-0.5,1]}; % outcome
        
        sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;       
        
        % Regression task -> MUA
        mua_taskpred_k = cell(length(regressors)+return_constant,n_depths);
        mua_taskpred = nan(size(event_aligned_mua));
        mua_taskpred_reduced = ...
            repmat(nan(size(event_aligned_mua)),1,1,1,length(regressors));       
        for curr_depth = 1:n_depths
                        
            baseline = nanmean(reshape(event_aligned_mua(:,t < 0,curr_depth),[],1)*raster_sample_rate);
            activity = single(binned_spikes(curr_depth,:)) - baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
                   
            [task_kernel,activity_predicted,expl_var,activity_predicted_reduced] = ...
                AP_regresskernel(regressors,activity,sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            mua_taskpred_k(:,curr_depth) = task_kernel;
            
            mua_taskpred(:,:,curr_depth) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event)./raster_sample_rate;
            
            mua_taskpred_reduced(:,:,curr_depth,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event)./raster_sample_rate,permute(1:length(regressors),[1,3,4,2]),'uni',false));       
            
        end
        
        % Regression task -> MUA-ctxpred
        mua_ctxpred_taskpred_k = cell(length(regressors)+return_constant,n_depths);
        mua_ctxpred_taskpred = nan(size(event_aligned_mua));
        mua_ctxpred_taskpred_reduced = ...
            repmat(nan(size(event_aligned_mua)),1,1,1,length(regressors));       
        for curr_depth = 1:n_depths
                      
            baseline = nanmean(reshape(event_aligned_mua_ctxpred(:,t < 0,curr_depth),[],1)*raster_sample_rate);
            activity = single(predicted_spikes(curr_depth,:))-baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
            
            % (there are NaNs because of skipped edges)
            nan_samples = isnan(activity);
            activity_predicted = nan(size(activity));
            activity_predicted_reduced = nan(size(activity,1), ...
                size(activity,2),length(regressors));
            
            [task_kernel,activity_predicted(~nan_samples), ...
                expl_var,activity_predicted_reduced(:,~nan_samples,:)] = ...
                AP_regresskernel( ...
                cellfun(@(x) x(:,~nan_samples),regressors,'uni',false), ...
                activity(~nan_samples),sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            mua_ctxpred_taskpred_k(:,curr_depth) = task_kernel;
            
            mua_ctxpred_taskpred(:,:,curr_depth) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event)./raster_sample_rate;
            
            mua_ctxpred_taskpred_reduced(:,:,curr_depth,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event)./raster_sample_rate,permute(1:length(regressors),[1,3,4,2]),'uni',false));       
            
        end
        
        % Regression task -> (master U, deconvolved) fluor
        n_vs = length(use_components);
        event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
        fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);
        
        fluor_taskpred_k = cell(length(regressors)+return_constant,n_vs);
        fluor_taskpred = nan(size(event_aligned_V));
        fluor_taskpred_reduced = ...
            repmat(nan(size(event_aligned_V)),1,1,1,length(regressors));       
        for curr_v = 1:n_vs
                      
            baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,curr_v),[],1));
            activity = single(fVdf_deconv_resample_recast(curr_v,:))-baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
            
            % (there are NaNs because of skipped edges)
            nan_samples = isnan(activity);
            activity_predicted = nan(size(activity));
            activity_predicted_reduced = nan(size(activity,1), ...
                size(activity,2),length(regressors));
            
            [task_kernel,activity_predicted(~nan_samples), ...
                expl_var,activity_predicted_reduced(:,~nan_samples,:)] = ...
                AP_regresskernel( ...
                cellfun(@(x) x(:,~nan_samples),regressors,'uni',false), ...
                activity(~nan_samples),sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            fluor_taskpred_k(:,curr_v) = task_kernel;
            
            fluor_taskpred(:,:,curr_v) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event);
            
            fluor_taskpred_reduced(:,:,curr_v,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event),permute(1:length(regressors),[1,3,4,2]),'uni',false));  
                        
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        mua_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_taskpred_k;
        mua_taskpred_all{curr_animal,1}{curr_day,1} = mua_taskpred(use_trials,:,:,:);
        mua_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_taskpred_reduced(use_trials,:,:,:);
        
        mua_ctxpred_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_k;
        mua_ctxpred_taskpred_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred(use_trials,:,:,:);
        mua_ctxpred_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_reduced(use_trials,:,:,:);
        
        fluor_taskpred_k_all{curr_animal,1}{curr_day,1} = fluor_taskpred_k;
        fluor_taskpred_all{curr_animal,1}{curr_day,1} = fluor_taskpred(use_trials,:,:,:);
        fluor_taskpred_reduced_all{curr_animal,1}{curr_day,1} = fluor_taskpred_reduced(use_trials,:,:,:);
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);   
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except curr_animal ...
        n_aligned_depths regression_params animals t ...
        fluor_all ...
        mua_all ...
        mua_ctxpred_all ...
        mua_taskpred_k_all ...
        mua_taskpred_all ...
        mua_taskpred_reduced_all ...
        mua_ctxpred_taskpred_k_all ...
        mua_ctxpred_taskpred_all ...
        mua_ctxpred_taskpred_reduced_all ...
        fluor_taskpred_k_all ...
        fluor_taskpred_all ...
        fluor_taskpred_reduced_all ...
        wheel_ctxpred_all ...
        ctx_str_k_all ...
        ctx_wheel_k_all ...
        wheel_all ...
        outcome_all ...
        D_all
    
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    fluor_taskpred_k_all ...
    fluor_taskpred_all ...
    fluor_taskpred_reduced_all ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    outcome_all ...
    D_all

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_MOVExSTIM'];
save([save_path filesep save_fn],'-v7.3');


%% Choiceworld trial activity (MSN ONLY)

clear all
disp('Choiceworld trial activity')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);
mua_taskpred_k_all = cell(1,1);
mua_taskpred_all = cell(1,1);
mua_taskpred_reduced_all = cell(1,1);
mua_ctxpred_taskpred_k_all = cell(1,1);
mua_ctxpred_taskpred_all = cell(1,1);
mua_ctxpred_taskpred_reduced_all = cell(1,1);
fluor_taskpred_k_all = cell(1,1);
fluor_taskpred_all = cell(1,1);
fluor_taskpred_reduced_all = cell(1,1);
wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;    
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;

        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Align (only to stim)
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);

        %%% MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            
%             curr_spikes = spike_times_timeline(depth_group == curr_depth);

            % (for only msns in depth group)
            curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
                ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf,[],2)',time_bin_centers)';
        
        %%%% TESTING: DECONV
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        %%%%
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            
%             curr_spikes = spike_times_timeline(depth_group == curr_depth);
            
            % (for only msns in depth group)
            curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
                ismember(spike_templates,find(msn)));
            
            binned_spikes(curr_depth,:) = histcounts(curr_spikes,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;        
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        %%%%% TEST
        warning('OVERRIDING LAMBDA')
        lambda = 5;
        %%%%%%%
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        %         [~,predicted_spikes_std,explained_var] = ...
        %             AP_regresskernel(dfVdf_resample, ...
        %             binned_spikes_std,kernel_frames,lambda, ...
        %             regression_params.zs,regression_params.cvfold, ...
        %             false,regression_params.use_constant);
        
        %%%% TESTING DECONV
        [ctx_str_k,predicted_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Shove the k's into the master to save small V's 
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
                     
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        predicted_spikes = (predicted_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,predicted_spikes',t_peri_event)./raster_sample_rate;
        
        %%%%
        
        %%%% TESTING CTX->WHEEL VELOCITY
        % Resample wheel, predict both velocity and speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Shove the k's into the master to save small V's 
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%%%
        
        %%% Wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)        
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_bins(x,:)), ...
            find(trial_outcome == 1)','uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_bins(x,:)), ...
            find(trial_outcome == -1)','uni',false))) > 0;
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);    
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING: TASK REGRESSION (LONG TRACE)
                
        %%% Regression (separate regressors to get partial explained variance)
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times)).* ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end    
        
        % Stim move regressors (one for each stim when it starts to move)       
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));            
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
                        
        end      
        
        % Stim center regressors (one for each stim when it's stopped during reward)    
        unique_contrasts = unique(contrasts(contrasts > 0));
        stim_contrasts = ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ... 
                stim_contrasts == unique_contrasts(curr_contrast));
                        
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));                  
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
                        
        end    

        % Move onset regressors (L/R)          
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move onset x stim regressors (one for each contrast/side)
        move_onset_stim_time_absolute = arrayfun(@(curr_stim) ...
            arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & stim_contrastsides' == unique_stim(curr_stim))), ...
            1:length(unique_stim),'uni',false);
        
        move_onset_stim_regressors = zeros(10,length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
           move_onset_stim_regressors(curr_stim,:) = ...
               histcounts(move_onset_stim_time_absolute{curr_stim},time_bins); 
        end
        
        % Move ongoing regressors (L/R choice for duration of movement)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        
        move_stopped_t = 0.5;
        move_stopped_samples = round(sample_rate*move_stopped_t);
        wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.02), ...
            ones(1,move_stopped_samples),'full') > 0;
        wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);
        
        move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_L_absolute','uni',false));
        move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_R_absolute','uni',false));
        
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
        move_ongoing_regressors(2,move_ongoing_R_samples) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        if length(signals_events.interactiveOnTimes) ~= length(move_t)
            error('Different number of interactive ons and move times')
        end
        
        go_cue_regressors = zeros(2,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        outcome_regressors = zeros(2,length(time_bin_centers));
        
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate regressors, set parameters
        regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        regressor_labels = {'Stim','Move onset','Go cue','Outcome'};     
        
        t_shifts = {[0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [-0.1,0.5]; ... % go cue
            [-0.5,1]}; % outcome
        
        sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;       
        
        % Regression task -> MUA
        mua_taskpred_k = cell(length(regressors)+return_constant,n_depths);
        mua_taskpred = nan(size(event_aligned_mua));
        mua_taskpred_reduced = ...
            repmat(nan(size(event_aligned_mua)),1,1,1,length(regressors));       
        for curr_depth = 1:n_depths
                        
            baseline = nanmean(reshape(event_aligned_mua(:,t < 0,curr_depth),[],1)*raster_sample_rate);
            activity = single(binned_spikes(curr_depth,:)) - baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
                   
            [task_kernel,activity_predicted,expl_var,activity_predicted_reduced] = ...
                AP_regresskernel(regressors,activity,sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            mua_taskpred_k(:,curr_depth) = task_kernel;
            
            mua_taskpred(:,:,curr_depth) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event)./raster_sample_rate;
            
            mua_taskpred_reduced(:,:,curr_depth,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event)./raster_sample_rate,permute(1:length(regressors),[1,3,4,2]),'uni',false));       
            
        end
        
        % Regression task -> MUA-ctxpred
        mua_ctxpred_taskpred_k = cell(length(regressors)+return_constant,n_depths);
        mua_ctxpred_taskpred = nan(size(event_aligned_mua));
        mua_ctxpred_taskpred_reduced = ...
            repmat(nan(size(event_aligned_mua)),1,1,1,length(regressors));       
        for curr_depth = 1:n_depths
                      
            baseline = nanmean(reshape(event_aligned_mua_ctxpred(:,t < 0,curr_depth),[],1)*raster_sample_rate);
            activity = single(predicted_spikes(curr_depth,:))-baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
            
            % (there are NaNs because of skipped edges)
            nan_samples = isnan(activity);
            activity_predicted = nan(size(activity));
            activity_predicted_reduced = nan(size(activity,1), ...
                size(activity,2),length(regressors));
            
            [task_kernel,activity_predicted(~nan_samples), ...
                expl_var,activity_predicted_reduced(:,~nan_samples,:)] = ...
                AP_regresskernel( ...
                cellfun(@(x) x(:,~nan_samples),regressors,'uni',false), ...
                activity(~nan_samples),sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            mua_ctxpred_taskpred_k(:,curr_depth) = task_kernel;
            
            mua_ctxpred_taskpred(:,:,curr_depth) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event)./raster_sample_rate;
            
            mua_ctxpred_taskpred_reduced(:,:,curr_depth,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event)./raster_sample_rate,permute(1:length(regressors),[1,3,4,2]),'uni',false));       
            
        end
        
        % Regression task -> (master U, deconvolved) fluor
        n_vs = length(use_components);
        event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
        fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);
        
        fluor_taskpred_k = cell(length(regressors)+return_constant,n_vs);
        fluor_taskpred = nan(size(event_aligned_V));
        fluor_taskpred_reduced = ...
            repmat(nan(size(event_aligned_V)),1,1,1,length(regressors));       
        for curr_v = 1:n_vs
                      
            baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,curr_v),[],1));
            activity = single(fVdf_deconv_resample_recast(curr_v,:))-baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
            
            % (there are NaNs because of skipped edges)
            nan_samples = isnan(activity);
            activity_predicted = nan(size(activity));
            activity_predicted_reduced = nan(size(activity,1), ...
                size(activity,2),length(regressors));
            
            [task_kernel,activity_predicted(~nan_samples), ...
                expl_var,activity_predicted_reduced(:,~nan_samples,:)] = ...
                AP_regresskernel( ...
                cellfun(@(x) x(:,~nan_samples),regressors,'uni',false), ...
                activity(~nan_samples),sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            fluor_taskpred_k(:,curr_v) = task_kernel;
            
            fluor_taskpred(:,:,curr_v) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event);
            
            fluor_taskpred_reduced(:,:,curr_v,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event),permute(1:length(regressors),[1,3,4,2]),'uni',false));  
                        
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        mua_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_taskpred_k;
        mua_taskpred_all{curr_animal,1}{curr_day,1} = mua_taskpred(use_trials,:,:,:);
        mua_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_taskpred_reduced(use_trials,:,:,:);
        
        mua_ctxpred_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_k;
        mua_ctxpred_taskpred_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred(use_trials,:,:,:);
        mua_ctxpred_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_reduced(use_trials,:,:,:);
        
        fluor_taskpred_k_all{curr_animal,1}{curr_day,1} = fluor_taskpred_k;
        fluor_taskpred_all{curr_animal,1}{curr_day,1} = fluor_taskpred(use_trials,:,:,:);
        fluor_taskpred_reduced_all{curr_animal,1}{curr_day,1} = fluor_taskpred_reduced(use_trials,:,:,:);
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);   
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except curr_animal ...
        n_aligned_depths regression_params animals t ...
        fluor_all ...
        mua_all ...
        mua_ctxpred_all ...
        mua_taskpred_k_all ...
        mua_taskpred_all ...
        mua_taskpred_reduced_all ...
        mua_ctxpred_taskpred_k_all ...
        mua_ctxpred_taskpred_all ...
        mua_ctxpred_taskpred_reduced_all ...
        fluor_taskpred_k_all ...
        fluor_taskpred_all ...
        fluor_taskpred_reduced_all ...
        wheel_ctxpred_all ...
        ctx_str_k_all ...
        ctx_wheel_k_all ...
        wheel_all ...
        outcome_all ...
        D_all
    
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    fluor_taskpred_k_all ...
    fluor_taskpred_all ...
    fluor_taskpred_reduced_all ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    outcome_all ...
    D_all

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_MSN'];
save([save_path filesep save_fn],'-v7.3');



%% Choiceworld trial activity (FSI ONLY)

clear all
disp('Choiceworld trial activity')

n_aligned_depths = 4; 

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);
mua_taskpred_k_all = cell(1,1);
mua_taskpred_all = cell(1,1);
mua_taskpred_reduced_all = cell(1,1);
mua_ctxpred_taskpred_k_all = cell(1,1);
mua_ctxpred_taskpred_all = cell(1,1);
mua_ctxpred_taskpred_reduced_all = cell(1,1);
fluor_taskpred_k_all = cell(1,1);
fluor_taskpred_all = cell(1,1);
fluor_taskpred_reduced_all = cell(1,1);
wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = false;
    load_parts.imaging = true;
    load_parts.ephys = true;    
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;

        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Align (only to stim)
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Fluorescence
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);

        %%% MUA
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            
%             curr_spikes = spike_times_timeline(depth_group == curr_depth);

            % (for only msns in depth group)
            curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
                ismember(spike_templates,find(fsi)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regressed fluorescence -> MUA
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Get upsampled dVdf's
        dfVdf_resample = interp1(conv2(frame_t,[1,1]/2,'valid'), ...
            diff(fVdf,[],2)',time_bin_centers)';
        
        %%%% TESTING: DECONV
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        %%%%
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            
%             curr_spikes = spike_times_timeline(depth_group == curr_depth);
            
            % (for only msns in depth group)
            curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
                ismember(spike_templates,find(fsi)));
            
            binned_spikes(curr_depth,:) = histcounts(curr_spikes,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;        
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        %%%%% TEST
        warning('OVERRIDING LAMBDA')
        lambda = 5;
        %%%%%%%
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        %         [~,predicted_spikes_std,explained_var] = ...
        %             AP_regresskernel(dfVdf_resample, ...
        %             binned_spikes_std,kernel_frames,lambda, ...
        %             regression_params.zs,regression_params.cvfold, ...
        %             false,regression_params.use_constant);
        
        %%%% TESTING DECONV
        [ctx_str_k,predicted_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Shove the k's into the master to save small V's 
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
                     
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        predicted_spikes = (predicted_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,predicted_spikes',t_peri_event)./raster_sample_rate;
        
        %%%%
        
        %%%% TESTING CTX->WHEEL VELOCITY
        % Resample wheel, predict both velocity and speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Shove the k's into the master to save small V's 
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%%%
        
        %%% Wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)        
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_bins(x,:)), ...
            find(trial_outcome == 1)','uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_bins(x,:)), ...
            find(trial_outcome == -1)','uni',false))) > 0;
        
        % Pick trials to use
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);    
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING: TASK REGRESSION (LONG TRACE)
                
        %%% Regression (separate regressors to get partial explained variance)
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times)).* ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end    
        
        % Stim move regressors (one for each stim when it starts to move)       
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));            
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
                        
        end      
        
        % Stim center regressors (one for each stim when it's stopped during reward)    
        unique_contrasts = unique(contrasts(contrasts > 0));
        stim_contrasts = ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ... 
                stim_contrasts == unique_contrasts(curr_contrast));
                        
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));                  
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
                        
        end    

        % Move onset regressors (L/R)          
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move onset x stim regressors (one for each contrast/side)
        move_onset_stim_time_absolute = arrayfun(@(curr_stim) ...
            arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & stim_contrastsides' == unique_stim(curr_stim))), ...
            1:length(unique_stim),'uni',false);
        
        move_onset_stim_regressors = zeros(10,length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
           move_onset_stim_regressors(curr_stim,:) = ...
               histcounts(move_onset_stim_time_absolute{curr_stim},time_bins); 
        end
        
        % Move ongoing regressors (L/R choice for duration of movement)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        
        move_stopped_t = 0.5;
        move_stopped_samples = round(sample_rate*move_stopped_t);
        wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.02), ...
            ones(1,move_stopped_samples),'full') > 0;
        wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);
        
        move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_L_absolute','uni',false));
        move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_R_absolute','uni',false));
        
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
        move_ongoing_regressors(2,move_ongoing_R_samples) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        if length(signals_events.interactiveOnTimes) ~= length(move_t)
            error('Different number of interactive ons and move times')
        end
        
        go_cue_regressors = zeros(2,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        outcome_regressors = zeros(2,length(time_bin_centers));
        
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate regressors, set parameters
        regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        regressor_labels = {'Stim','Move onset','Go cue','Outcome'};     
        
        t_shifts = {[0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [-0.1,0.5]; ... % go cue
            [-0.5,1]}; % outcome
        
        sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;       
        
        % Regression task -> MUA
        mua_taskpred_k = cell(length(regressors)+return_constant,n_depths);
        mua_taskpred = nan(size(event_aligned_mua));
        mua_taskpred_reduced = ...
            repmat(nan(size(event_aligned_mua)),1,1,1,length(regressors));       
        for curr_depth = 1:n_depths
                        
            baseline = nanmean(reshape(event_aligned_mua(:,t < 0,curr_depth),[],1)*raster_sample_rate);
            activity = single(binned_spikes(curr_depth,:)) - baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
                   
            [task_kernel,activity_predicted,expl_var,activity_predicted_reduced] = ...
                AP_regresskernel(regressors,activity,sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            mua_taskpred_k(:,curr_depth) = task_kernel;
            
            mua_taskpred(:,:,curr_depth) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event)./raster_sample_rate;
            
            mua_taskpred_reduced(:,:,curr_depth,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event)./raster_sample_rate,permute(1:length(regressors),[1,3,4,2]),'uni',false));       
            
        end
        
        % Regression task -> MUA-ctxpred
        mua_ctxpred_taskpred_k = cell(length(regressors)+return_constant,n_depths);
        mua_ctxpred_taskpred = nan(size(event_aligned_mua));
        mua_ctxpred_taskpred_reduced = ...
            repmat(nan(size(event_aligned_mua)),1,1,1,length(regressors));       
        for curr_depth = 1:n_depths
                      
            baseline = nanmean(reshape(event_aligned_mua_ctxpred(:,t < 0,curr_depth),[],1)*raster_sample_rate);
            activity = single(predicted_spikes(curr_depth,:))-baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
            
            % (there are NaNs because of skipped edges)
            nan_samples = isnan(activity);
            activity_predicted = nan(size(activity));
            activity_predicted_reduced = nan(size(activity,1), ...
                size(activity,2),length(regressors));
            
            [task_kernel,activity_predicted(~nan_samples), ...
                expl_var,activity_predicted_reduced(:,~nan_samples,:)] = ...
                AP_regresskernel( ...
                cellfun(@(x) x(:,~nan_samples),regressors,'uni',false), ...
                activity(~nan_samples),sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            mua_ctxpred_taskpred_k(:,curr_depth) = task_kernel;
            
            mua_ctxpred_taskpred(:,:,curr_depth) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event)./raster_sample_rate;
            
            mua_ctxpred_taskpred_reduced(:,:,curr_depth,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event)./raster_sample_rate,permute(1:length(regressors),[1,3,4,2]),'uni',false));       
            
        end
        
        % Regression task -> (master U, deconvolved) fluor
        n_vs = length(use_components);
        event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
        fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);
        
        fluor_taskpred_k = cell(length(regressors)+return_constant,n_vs);
        fluor_taskpred = nan(size(event_aligned_V));
        fluor_taskpred_reduced = ...
            repmat(nan(size(event_aligned_V)),1,1,1,length(regressors));       
        for curr_v = 1:n_vs
                      
            baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,curr_v),[],1));
            activity = single(fVdf_deconv_resample_recast(curr_v,:))-baseline;
            
            % Skip if nothing in this depth
            if ~any(activity(:))
                continue
            end
            
            % (there are NaNs because of skipped edges)
            nan_samples = isnan(activity);
            activity_predicted = nan(size(activity));
            activity_predicted_reduced = nan(size(activity,1), ...
                size(activity,2),length(regressors));
            
            [task_kernel,activity_predicted(~nan_samples), ...
                expl_var,activity_predicted_reduced(:,~nan_samples,:)] = ...
                AP_regresskernel( ...
                cellfun(@(x) x(:,~nan_samples),regressors,'uni',false), ...
                activity(~nan_samples),sample_shifts, ...
                lambda,zs,cvfold,return_constant,use_constant);
            
            fluor_taskpred_k(:,curr_v) = task_kernel;
            
            fluor_taskpred(:,:,curr_v) = ...
                interp1(time_bin_centers,activity_predicted',t_peri_event);
            
            fluor_taskpred_reduced(:,:,curr_v,:) = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,activity_predicted_reduced(:,:,x)', ...
                t_peri_event),permute(1:length(regressors),[1,3,4,2]),'uni',false));  
                        
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        mua_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_taskpred_k;
        mua_taskpred_all{curr_animal,1}{curr_day,1} = mua_taskpred(use_trials,:,:,:);
        mua_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_taskpred_reduced(use_trials,:,:,:);
        
        mua_ctxpred_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_k;
        mua_ctxpred_taskpred_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred(use_trials,:,:,:);
        mua_ctxpred_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_reduced(use_trials,:,:,:);
        
        fluor_taskpred_k_all{curr_animal,1}{curr_day,1} = fluor_taskpred_k;
        fluor_taskpred_all{curr_animal,1}{curr_day,1} = fluor_taskpred(use_trials,:,:,:);
        fluor_taskpred_reduced_all{curr_animal,1}{curr_day,1} = fluor_taskpred_reduced(use_trials,:,:,:);
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);   
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
    end
    
    clearvars -except curr_animal ...
        n_aligned_depths regression_params animals t ...
        fluor_all ...
        mua_all ...
        mua_ctxpred_all ...
        mua_taskpred_k_all ...
        mua_taskpred_all ...
        mua_taskpred_reduced_all ...
        mua_ctxpred_taskpred_k_all ...
        mua_ctxpred_taskpred_all ...
        mua_ctxpred_taskpred_reduced_all ...
        fluor_taskpred_k_all ...
        fluor_taskpred_all ...
        fluor_taskpred_reduced_all ...
        wheel_ctxpred_all ...
        ctx_str_k_all ...
        ctx_wheel_k_all ...
        wheel_all ...
        outcome_all ...
        D_all
    
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    fluor_taskpred_k_all ...
    fluor_taskpred_all ...
    fluor_taskpred_reduced_all ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    outcome_all ...
    D_all

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_FSI'];
save([save_path filesep save_fn],'-v7.3');




%% Choiceworld trial activity (TESTING 200 UM DEPTH ALIGN)

clear all
disp('Choiceworld trial activity (by 200 um depths)')

n_aligned_depths = 15; % (this gives 197 um segments)

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);

mua_taskpred_k_all = cell(1,1);
mua_taskpred_all = cell(1,1);
mua_taskpred_reduced_all = cell(1,1);
mua_taskpred_expl_var_total_all = cell(1,1);
mua_taskpred_expl_var_partial_all = cell(1,1);

mua_ctxpred_taskpred_k_all = cell(1,1);
mua_ctxpred_taskpred_all = cell(1,1);
mua_ctxpred_taskpred_reduced_all = cell(1,1);
mua_ctxpred_taskpred_expl_var_total_all = cell(1,1);
mua_ctxpred_taskpred_expl_var_partial_all = cell(1,1);

fluor_taskpred_k_all = cell(1,1);
fluor_taskpred_all = cell(1,1);
fluor_taskpred_reduced_all = cell(1,1);
fluor_taskpred_expl_var_total_all = cell(1,1);
fluor_taskpred_expl_var_partial_all = cell(1,1);

wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
movement_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    load_parts.cam = true;
    load_parts.imaging = true;
    load_parts.ephys = true;    
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'depth';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Get align times
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        
        %%% Trial-align cortex
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% Trial-align striatum
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        t_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regress cortex to striatum
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;

        % Deconvolve fluoresence
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        
        binned_spikes = zeros(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        binned_spikes_std(isnan(binned_spikes_std)) = 0;        
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        % Regress cortex to striatum
        [ctx_str_k,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Recast the k's into the master U
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
                     
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        ctxpred_spikes = (ctxpred_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,ctxpred_spikes',t_peri_event)./raster_sample_rate;
                       
        % Regress cortex to wheel velocity/speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Recast the k's into the master U
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
                
        %%% Trial-align wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Trial-align facecam movement
        event_aligned_movement = interp1(facecam_t(~isnan(facecam_t)), ...
            frame_movement(~isnan(facecam_t)),t_peri_event);
        
        %%% Trial-align outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)        
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_bins(x,:)), ...
            find(trial_outcome == 1)','uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_bins(x,:)), ...
            find(trial_outcome == -1)','uni',false))) > 0;
        
        % Pick trials to keep
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials) & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials) == -1;
        R_trials = signals_events.trialSideValues(1:n_trials) == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = 3-(abs((trial_choice(use_trials)'+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);    
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
                        
        %%% Regress task to cortex/striatum/cortex-predicted striatum
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times)).* ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end    
        
        % Stim move regressors (one for each stim when it starts to move)       
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));            
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
                        
        end      
        
        % Stim center regressors (one for each stim when it's stopped during reward)    
        unique_contrasts = unique(contrasts(contrasts > 0));
        stim_contrasts = ...
            signals_events.trialContrastValues(1:length(stimOn_times));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ... 
                stim_contrasts == unique_contrasts(curr_contrast));
                        
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));                  
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
                        
        end    

        % Move onset regressors (L/R)          
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times))' == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move onset x stim regressors (one for each contrast/side)
        move_onset_stim_time_absolute = arrayfun(@(curr_stim) ...
            arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & stim_contrastsides' == unique_stim(curr_stim))), ...
            1:length(unique_stim),'uni',false);
        
        move_onset_stim_regressors = zeros(10,length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
           move_onset_stim_regressors(curr_stim,:) = ...
               histcounts(move_onset_stim_time_absolute{curr_stim},time_bins); 
        end
        
        % Move ongoing regressors (L/R choice for duration of movement)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        
        move_stopped_t = 0.5;
        move_stopped_samples = round(sample_rate*move_stopped_t);
        wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.02), ...
            ones(1,move_stopped_samples),'full') > 0;
        wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);
        
        move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_L_absolute','uni',false));
        move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_R_absolute','uni',false));
        
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
        move_ongoing_regressors(2,move_ongoing_R_samples) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        if length(signals_events.interactiveOnTimes) ~= length(move_t)
            error('Different number of interactive ons and move times')
        end
        
        go_cue_regressors = zeros(2,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        outcome_regressors = zeros(2,length(time_bin_centers));
        
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate selected regressors, set parameters
        regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        regressor_labels = {'Stim','Move onset','Go cue','Outcome'};     
        
        t_shifts = {[0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [-0.1,0.5]; ... % go cue
            [-0.5,1]}; % outcome
        
        sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;  
        
        
        % Regression task -> MUA
        baseline = nanmean(reshape(event_aligned_mua(:,t < 0,:),[], ...
            size(event_aligned_mua,3))*raster_sample_rate,1)';
        activity = single(binned_spikes) - baseline;
        
        [mua_taskpred_k,mua_taskpred_long,mua_taskpred_expl_var,mua_taskpred_reduced_long] = ...
            AP_regresskernel(regressors,activity,sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        mua_taskpred = ...
            interp1(time_bin_centers,mua_taskpred_long',t_peri_event)./raster_sample_rate;
        
        mua_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,mua_taskpred_reduced_long(:,:,x)', ...
            t_peri_event)./raster_sample_rate,permute(1:length(regressors),[1,3,4,2]),'uni',false));
        
        % Regression task -> MUA-ctxpred
        baseline = nanmean(reshape(event_aligned_mua_ctxpred(:,t < 0,:),[], ...
            size(event_aligned_mua_ctxpred,3))*raster_sample_rate,1)';
        activity = single(ctxpred_spikes) - baseline;
        
        [mua_ctxpred_taskpred_k,mua_ctxpred_taskpred_long,mua_ctxpred_taskpred_expl_var,mua_ctxpred_taskpred_reduced_long] = ...
            AP_regresskernel(regressors,activity,sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        mua_ctxpred_taskpred = ...
            interp1(time_bin_centers,mua_ctxpred_taskpred_long',t_peri_event)./raster_sample_rate;
        
        mua_ctxpred_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,mua_ctxpred_taskpred_reduced_long(:,:,x)', ...
            t_peri_event)./raster_sample_rate,permute(1:length(regressors),[1,3,4,2]),'uni',false));
        
        % Regression task -> (master U, deconvolved) fluor
        event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
        fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);
        
        baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,:),[],size(event_aligned_V_deconv,3)))';
        activity = single(fVdf_deconv_resample_recast(use_components,:))-baseline;
        
        [fluor_taskpred_k,fluor_taskpred_long,fluor_taskpred_expl_var,fluor_taskpred_reduced_long] = ...
            AP_regresskernel(regressors,activity,sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        fluor_taskpred = ...
            interp1(time_bin_centers,fluor_taskpred_long',t_peri_event);
        
        fluor_taskpred_reduced = cell2mat(arrayfun(@(x) ...
                interp1(time_bin_centers,fluor_taskpred_reduced_long(:,:,x)', ...
                t_peri_event),permute(1:length(regressors),[1,3,4,2]),'uni',false));  
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        mua_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_taskpred_k;
        mua_taskpred_all{curr_animal,1}{curr_day,1} = mua_taskpred(use_trials,:,:,:);
        mua_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_taskpred_reduced(use_trials,:,:,:);
        mua_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = mua_taskpred_expl_var.total;
        mua_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = mua_taskpred_expl_var.partial;
        
        mua_ctxpred_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_k;
        mua_ctxpred_taskpred_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred(use_trials,:,:,:);
        mua_ctxpred_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_reduced(use_trials,:,:,:);
        mua_ctxpred_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_expl_var.total;
        mua_ctxpred_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_expl_var.partial;
        
        fluor_taskpred_k_all{curr_animal,1}{curr_day,1} = fluor_taskpred_k;
        fluor_taskpred_all{curr_animal,1}{curr_day,1} = fluor_taskpred(use_trials,:,:,:);
        fluor_taskpred_reduced_all{curr_animal,1}{curr_day,1} = fluor_taskpred_reduced(use_trials,:,:,:);
        fluor_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.total;
        fluor_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.partial;
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        movement_all{curr_animal,1}{curr_day,1} = event_aligned_movement(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);   
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars -except curr_animal animal protocol experiments load_parts curr_day ...
        n_aligned_depths regression_params animals t ...
        ...
        fluor_all ...
        mua_all ...
        mua_ctxpred_all ...
        ...
        mua_taskpred_k_all ...
        mua_taskpred_all ...
        mua_taskpred_reduced_all ...
        mua_taskpred_expl_var_total_all ...
        mua_taskpred_expl_var_partial_all ...
        ...
        mua_ctxpred_taskpred_k_all ...
        mua_ctxpred_taskpred_all ...
        mua_ctxpred_taskpred_reduced_all ...
        mua_ctxpred_taskpred_expl_var_total_all ...
        mua_ctxpred_taskpred_expl_var_partial_all ...
        ...
        fluor_taskpred_k_all ...
        fluor_taskpred_all ...
        fluor_taskpred_reduced_all ...
        fluor_taskpred_expl_var_total_all ...
        fluor_taskpred_expl_var_partial_all ...
        ...
        wheel_ctxpred_all ...
        ctx_str_k_all ...
        ctx_wheel_k_all ...
        wheel_all ...
        movement_all ...
        outcome_all ...
        D_all
    
    end  
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_taskpred_expl_var_total_all ...
    mua_taskpred_expl_var_partial_all ...
    ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    mua_ctxpred_taskpred_expl_var_total_all ...
    mua_ctxpred_taskpred_expl_var_partial_all ...
    ...
    fluor_taskpred_k_all ...
    fluor_taskpred_all ...
    fluor_taskpred_reduced_all ...
    fluor_taskpred_expl_var_total_all ...
    fluor_taskpred_expl_var_partial_all ...
    ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    movement_all ...
    outcome_all ...
    D_all

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_200umdepth'];
save([save_path filesep save_fn],'-v7.3');




%% temp plot for cosyne

figure; hold on

AP_stackplot(squeeze(nanmean(mua_allcat,1)),t,5,false,'k',1:4);

% AP_stackplot(squeeze(nanmean(mua_allcat(trial_contrastside_allcat > 0,:,:),1)),t,5,false,'r',1:4);
% AP_stackplot(squeeze(nanmean(mua_allcat(trial_contrastside_allcat < 0,:,:),1)),t,5,false,'b',1:4);
% 
% AP_stackplot(squeeze(nanmean(mua_allcat(trial_choice_allcat == -1,:,:),1)),t,5,false,[0.6,0,0.6],1:4);
% AP_stackplot(squeeze(nanmean(mua_allcat(trial_choice_allcat == 1,:,:),1)),t,5,false,[0,0.7,0],1:4);

% AP_stackplot(squeeze(nanmean(mua_allcat(move_t > 0 & move_t < 0.2,:,:),1)),t,5,false,[0,0.3,0.3],1:4);
% AP_stackplot(squeeze(nanmean(mua_allcat(move_t > 0 & move_t < 0.5,:,:),1)),t,5,false,[0,0.5,0.5],1:4);
% AP_stackplot(squeeze(nanmean(mua_allcat(move_t > 0.5 & move_t < 1,:,:),1)),t,5,false,[0,0.7,0.7],1:4);

% AP_stackplot(squeeze(nanmean(mua_allcat(trial_outcome_allcat == -1,:,:),1)),t,5,false,'m',1:4);
% AP_stackplot(squeeze(nanmean(mua_allcat(trial_outcome_allcat == 1,:,:),1)),t,5,false,'c',1:4);

spacing = 7;
figure; hold on;
plot_trials = move_t > 0 & move_t < 0.5 & trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & trial_choice_allcat == -1;
AP_stackplot(squeeze(nanmean(mua_allcat(plot_trials,:,:),1)),t,spacing,false,'k',1:4);
AP_stackplot(squeeze(nanmean(mua_ctxpred_allcat(plot_trials,:,:),1)),t,spacing,false,[0,0.7,0],1:4);
AP_stackplot(squeeze(nanmean(mua_allcat(plot_trials,:,:),1) - ...
    nanmean(mua_ctxpred_allcat(plot_trials,:,:),1)),t,spacing,false,[0.5,0.5,0.5],1:4);
line([0,0],ylim,'color','k');
line([0,0.5],[0,0],'linewidth',5);
line([0,0],[0,1],'linewidth',5);

figure; 
plot_trials = move_t > 0 & move_t < 0.5 & trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & trial_choice_allcat == -1;
AP_stackplot(squeeze(nanmean(mua_allcat(plot_trials,:,:),1)),t,5,false,'k',1:4);
line([0,0],ylim,'color','k');
line([0,0.5],[0,0],'linewidth',5);
line([0,0],[0,1],'linewidth',5);




figure;
plot_trials = move_t > 0 & move_t < 0.5 & trial_contrast_allcat > 0 & ...
    trial_side_allcat == 1 & trial_choice_allcat == -1;

subplot(1,2,1); hold on;
plot(t,nanmean(mua_allcat(plot_trials,:,2),1),'k','linewidth',2);
plot(t,nanmean(mua_allcat(plot_trials,:,2),1) - ...
    nanmean(mua_taskpred_reduced_allcat(plot_trials,:,2,1),1),'r','linewidth',2);

subplot(1,2,2); hold on;
plot(t,nanmean(mua_allcat_move(plot_trials,:,2),1),'k','linewidth',2);
plot(t,nanmean(mua_allcat_move(plot_trials,:,2),1) - ...
    nanmean(mua_taskpred_reduced_allcat_move(plot_trials,:,2,2),1),'color',[0.7,0,0.7],'linewidth',2);



%% (choiceworld trials: trying fit hemidiff)

clear all
disp('Choiceworld trial activity')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);

mua_taskpred_k_all = cell(1,1);
mua_taskpred_all = cell(1,1);
mua_taskpred_reduced_all = cell(1,1);
mua_taskpred_expl_var_total_all = cell(1,1);
mua_taskpred_expl_var_partial_all = cell(1,1);

mua_ctxpred_taskpred_k_all = cell(1,1);
mua_ctxpred_taskpred_all = cell(1,1);
mua_ctxpred_taskpred_reduced_all = cell(1,1);
mua_ctxpred_taskpred_expl_var_total_all = cell(1,1);
mua_ctxpred_taskpred_expl_var_partial_all = cell(1,1);

fluor_taskpred_k_all = cell(1,1);
fluor_taskpred_all = cell(1,1);
fluor_taskpred_reduced_all = cell(1,1);
fluor_taskpred_expl_var_total_all = cell(1,1);
fluor_taskpred_expl_var_partial_all = cell(1,1);

wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
movement_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Get align times
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        
        %%% Trial-align cortex
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% Trial-align striatum
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            if isempty(curr_spikes)
                continue
            end
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regress cortex to striatum
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Deconvolve fluoresence
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        
        % Bin spikes across the experiment
        binned_spikes = nan(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            % Skip if no spikes at this depth
            if isempty(curr_spike_times)
                continue
            end
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        % Regress cortex to striatum
        [ctx_str_k,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        %%%% TEST HEMIDIFF V'S HERE
        
        n_vs = size(U,3);
        mirror_matrix = reshape(Udf_aligned(:,:,1:n_vs),[],n_vs)'* ...
            reshape(AP_reflect_widefield(Udf_aligned(:,:,1:n_vs)),[],n_vs);
        fVdf_deconv_resample_mirror = reshape(transpose( ...
            mirror_matrix*reshape(fVdf_deconv_resample,[],n_vs)'),size(fVdf_deconv_resample));
        fVdf_deconv_resample_hemidiff = fVdf_deconv_resample - fVdf_deconv_resample_mirror;

        [ctx_str_k,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample_hemidiff(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        % Recast the k's into the master U
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
        
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        ctxpred_spikes = (ctxpred_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,ctxpred_spikes',t_peri_event)./raster_sample_rate;
        
        % Regress cortex to wheel velocity/speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Recast the k's into the master U
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%% Trial-align wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Trial-align facecam movement
        event_aligned_movement = interp1(facecam_t(~isnan(facecam_t)), ...
            frame_movement(~isnan(facecam_t)),t_peri_event);
        
        %%% Trial-align outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_peri_event_bins(x,:)), ...
            find(trial_outcome == 1),'uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_peri_event_bins(x,:)), ...
            find(trial_outcome == -1),'uni',false))) > 0;
        
        % Pick trials to keep
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials)' & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials)' == -1;
        R_trials = signals_events.trialSideValues(1:n_trials)' == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = 3-(abs((trial_choice(use_trials)+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%% Regress task to cortex/striatum/cortex-predicted striatum
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Build regressors (only a subset of these are used)
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times))'.* ...
            signals_events.trialContrastValues(1:length(stimOn_times))';
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end
        
        % Stim move regressors (one for each stim when it starts to move)
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
            
        end
        
        % Stim center regressors (one for each stim when it's stopped during reward)
        unique_contrasts = unique(contrasts(contrasts > 0));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ...
                abs(stim_contrastsides) == unique_contrasts(curr_contrast));
            
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
            
        end
        
        % Move onset regressors (L/R)
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move onset x stim regressors (one for each contrast/side)
        move_onset_stim_time_absolute = arrayfun(@(curr_stim) ...
            arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & stim_contrastsides == unique_stim(curr_stim))), ...
            1:length(unique_stim),'uni',false);
        
        move_onset_stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            move_onset_stim_regressors(curr_stim,:) = ...
                histcounts(move_onset_stim_time_absolute{curr_stim},time_bins);
        end
        
        % Move ongoing regressors (L/R choice for duration of movement)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        
        move_stopped_t = 0.5;
        move_stopped_samples = round(sample_rate*move_stopped_t);
        wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.02), ...
            ones(1,move_stopped_samples),'full') > 0;
        wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);
        
        move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_L_absolute','uni',false));
        move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_R_absolute','uni',false));
        
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
        move_ongoing_regressors(2,move_ongoing_R_samples) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        % (for go cue only on late move trials)
%         go_cue_regressors = histcounts( ...
%             signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        % (for go cue with early/late move trials)
        go_cue_regressors = zeros(1,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        % (regressors for hit only)
%         outcome_regressors = histcounts(reward_t_timeline,time_bins);
        % (regressors for both hit and miss)
        outcome_regressors = zeros(2,length(time_bin_centers));
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate selected regressors, set parameters
        
        task_regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        task_regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
                       
        task_t_shifts = { ...
            [0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [0,0.5]; ... % go cue
            [0,0.5]}; % outcome       
        
        % (old extended timings)
        %         task_t_shifts = {[0,0.5]; ... % stim
        %             [-0.5,1]; ... % move
        %             [-0.1,0.5]; ... % go cue
        %             [-0.5,1]}; % outcome

        % (to include stim x move)
        %         task_regressors = {stim_regressors;move_onset_regressors;move_onset_stim_regressors;go_cue_regressors;outcome_regressors};
        %         task_regressor_labels = {'Stim','Move','Stim x move','Go cue','Outcome'};
        %
        %         task_t_shifts = { ...
        %             [0,0.5]; ... % stim
        %             [-0.5,1]; ... % move
        %             [-0.5,1]; ... % stim x move
        %             [0,0.5]; ... % go cue
        %             [0,0.5]}; % outcome
        
        task_regressor_sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),task_t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;
        
        % Regression task -> MUA
        baseline = nanmean(reshape(event_aligned_mua(:,t < 0,:),[], ...
            size(event_aligned_mua,3))*raster_sample_rate,1)';
        activity = single(binned_spikes) - baseline;
        
        [mua_taskpred_k,mua_taskpred_long,mua_taskpred_expl_var,mua_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        mua_taskpred = ...
            interp1(time_bin_centers,mua_taskpred_long',t_peri_event)./raster_sample_rate;
        
        mua_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,mua_taskpred_reduced_long(:,:,x)', ...
            t_peri_event)./raster_sample_rate,permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        % Regression task -> MUA-ctxpred
        baseline = nanmean(reshape(event_aligned_mua_ctxpred(:,t < 0,:),[], ...
            size(event_aligned_mua_ctxpred,3))*raster_sample_rate,1)';
        activity = single(ctxpred_spikes) - baseline;
        
        [mua_ctxpred_taskpred_k,mua_ctxpred_taskpred_long,mua_ctxpred_taskpred_expl_var,mua_ctxpred_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        mua_ctxpred_taskpred = ...
            interp1(time_bin_centers,mua_ctxpred_taskpred_long',t_peri_event)./raster_sample_rate;
        
        mua_ctxpred_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,mua_ctxpred_taskpred_reduced_long(:,:,x)', ...
            t_peri_event)./raster_sample_rate,permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        % Regression task -> (master U, deconvolved) fluor
        event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
        fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);
        
        baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,:),[],size(event_aligned_V_deconv,3)))';
        activity = single(fVdf_deconv_resample_recast(use_components,:))-baseline;
        
        [fluor_taskpred_k,fluor_taskpred_long,fluor_taskpred_expl_var,fluor_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        fluor_taskpred = ...
            interp1(time_bin_centers,fluor_taskpred_long',t_peri_event);
        
        fluor_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,fluor_taskpred_reduced_long(:,:,x)', ...
            t_peri_event),permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        mua_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_taskpred_k;
        mua_taskpred_all{curr_animal,1}{curr_day,1} = mua_taskpred(use_trials,:,:,:);
        mua_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_taskpred_reduced(use_trials,:,:,:);
        mua_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = mua_taskpred_expl_var.total;
        mua_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = mua_taskpred_expl_var.partial;
        
        mua_ctxpred_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_k;
        mua_ctxpred_taskpred_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred(use_trials,:,:,:);
        mua_ctxpred_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_reduced(use_trials,:,:,:);
        mua_ctxpred_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_expl_var.total;
        mua_ctxpred_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_expl_var.partial;
        
        fluor_taskpred_k_all{curr_animal,1}{curr_day,1} = fluor_taskpred_k;
        fluor_taskpred_all{curr_animal,1}{curr_day,1} = fluor_taskpred(use_trials,:,:,:);
        fluor_taskpred_reduced_all{curr_animal,1}{curr_day,1} = fluor_taskpred_reduced(use_trials,:,:,:);
        fluor_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.total;
        fluor_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.partial;
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        movement_all{curr_animal,1}{curr_day,1} = event_aligned_movement(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars -except curr_animal animal protocol experiments load_parts curr_day ...
            n_aligned_depths regression_params animals t ...
            ...
            fluor_all ...
            mua_all ...
            mua_ctxpred_all ...
            ...
            mua_taskpred_k_all ...
            mua_taskpred_all ...
            mua_taskpred_reduced_all ...
            mua_taskpred_expl_var_total_all ...
            mua_taskpred_expl_var_partial_all ...
            ...
            mua_ctxpred_taskpred_k_all ...
            mua_ctxpred_taskpred_all ...
            mua_ctxpred_taskpred_reduced_all ...
            mua_ctxpred_taskpred_expl_var_total_all ...
            mua_ctxpred_taskpred_expl_var_partial_all ...
            ...
            fluor_taskpred_k_all ...
            fluor_taskpred_all ...
            fluor_taskpred_reduced_all ...
            fluor_taskpred_expl_var_total_all ...
            fluor_taskpred_expl_var_partial_all ...
            ...
            wheel_ctxpred_all ...
            ctx_str_k_all ...
            ctx_wheel_k_all ...
            wheel_all ...
            movement_all ...
            outcome_all ...
            D_all ...
            ...
            task_regressor_labels ...
            task_regressor_sample_shifts
        
    end
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_taskpred_expl_var_total_all ...
    mua_taskpred_expl_var_partial_all ...
    ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    mua_ctxpred_taskpred_expl_var_total_all ...
    mua_ctxpred_taskpred_expl_var_partial_all ...
    ...
    fluor_taskpred_k_all ...
    fluor_taskpred_all ...
    fluor_taskpred_reduced_all ...
    fluor_taskpred_expl_var_total_all ...
    fluor_taskpred_expl_var_partial_all ...
    ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    movement_all ...
    outcome_all ...
    D_all ...
    ...
    task_regressor_labels ...
    task_regressor_sample_shifts

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_hemidiff'];
save([save_path filesep save_fn],'-v7.3');



%% (choiceworld trials: trying str prediction from ctx L only)

clear all
disp('Choiceworld trial activity')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);

mua_taskpred_k_all = cell(1,1);
mua_taskpred_all = cell(1,1);
mua_taskpred_reduced_all = cell(1,1);
mua_taskpred_expl_var_total_all = cell(1,1);
mua_taskpred_expl_var_partial_all = cell(1,1);

mua_ctxpred_taskpred_k_all = cell(1,1);
mua_ctxpred_taskpred_all = cell(1,1);
mua_ctxpred_taskpred_reduced_all = cell(1,1);
mua_ctxpred_taskpred_expl_var_total_all = cell(1,1);
mua_ctxpred_taskpred_expl_var_partial_all = cell(1,1);

fluor_taskpred_k_all = cell(1,1);
fluor_taskpred_all = cell(1,1);
fluor_taskpred_reduced_all = cell(1,1);
fluor_taskpred_expl_var_total_all = cell(1,1);
fluor_taskpred_expl_var_partial_all = cell(1,1);

wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
movement_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Get align times
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        
        %%% Trial-align cortex
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% Trial-align striatum
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            if isempty(curr_spikes)
                continue
            end
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regress cortex to striatum
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Deconvolve fluoresence
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        
        % Bin spikes across the experiment
        binned_spikes = nan(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            % Skip if no spikes at this depth
            if isempty(curr_spike_times)
                continue
            end
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        % Regress cortex to striatum
        [ctx_str_k,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        %%%% TEST CtxL-only V'S HERE
        
        % Load bregma and master CCF tform
        bregma = allenCCFbregma;
        bregma(3) = bregma(3) + 0.5;
        ccf_tform_fn = ['C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF\ccf_tform'];
        load(ccf_tform_fn);       
        
        um2pixel = 20.6;
        bregma_resize = bregma*(10/um2pixel);
        bregma_align = [bregma_resize([3,1]),1]*ccf_tform.T;
        
        Udf_aligned_L = Udf_aligned;
        Udf_aligned_L(:,round(bregma_align(1)):end,:) = 0;
        
        n_vs = size(U,3);
        U_tform_matrix = reshape(Udf_aligned(:,:,1:n_vs),[],n_vs)'* ...
            reshape(Udf_aligned_L,[],n_vs);
        fVdf_deconv_resample_L = reshape(transpose( ...
            U_tform_matrix*reshape(fVdf_deconv_resample,[],n_vs)'),size(fVdf_deconv_resample));

        [ctx_str_k,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample_L(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Recast the k's into the master U
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned_L(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        
%         % Recast the k's into the master U
%         ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
%             reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
%             size(ctx_str_k{1}));
        
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        ctxpred_spikes = (ctxpred_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,ctxpred_spikes',t_peri_event)./raster_sample_rate;
        
        % Regress cortex to wheel velocity/speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Recast the k's into the master U
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%% Trial-align wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Trial-align facecam movement
        event_aligned_movement = interp1(facecam_t(~isnan(facecam_t)), ...
            frame_movement(~isnan(facecam_t)),t_peri_event);
        
        %%% Trial-align outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_peri_event_bins(x,:)), ...
            find(trial_outcome == 1),'uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_peri_event_bins(x,:)), ...
            find(trial_outcome == -1),'uni',false))) > 0;
        
        % Pick trials to keep
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials)' & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials)' == -1;
        R_trials = signals_events.trialSideValues(1:n_trials)' == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = 3-(abs((trial_choice(use_trials)+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%% Regress task to cortex/striatum/cortex-predicted striatum
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Build regressors (only a subset of these are used)
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times))'.* ...
            signals_events.trialContrastValues(1:length(stimOn_times))';
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end
        
        % Stim move regressors (one for each stim when it starts to move)
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
            
        end
        
        % Stim center regressors (one for each stim when it's stopped during reward)
        unique_contrasts = unique(contrasts(contrasts > 0));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ...
                abs(stim_contrastsides) == unique_contrasts(curr_contrast));
            
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
            
        end
        
        % Move onset regressors (L/R)
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move onset x stim regressors (one for each contrast/side)
        move_onset_stim_time_absolute = arrayfun(@(curr_stim) ...
            arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & stim_contrastsides == unique_stim(curr_stim))), ...
            1:length(unique_stim),'uni',false);
        
        move_onset_stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            move_onset_stim_regressors(curr_stim,:) = ...
                histcounts(move_onset_stim_time_absolute{curr_stim},time_bins);
        end
        
        % Move ongoing regressors (L/R choice for duration of movement)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        
        move_stopped_t = 0.5;
        move_stopped_samples = round(sample_rate*move_stopped_t);
        wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.02), ...
            ones(1,move_stopped_samples),'full') > 0;
        wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);
        
        move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_L_absolute','uni',false));
        move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_R_absolute','uni',false));
        
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
        move_ongoing_regressors(2,move_ongoing_R_samples) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        % (for go cue only on late move trials)
%         go_cue_regressors = histcounts( ...
%             signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        % (for go cue with early/late move trials)
        go_cue_regressors = zeros(1,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        % (regressors for hit only)
%         outcome_regressors = histcounts(reward_t_timeline,time_bins);
        % (regressors for both hit and miss)
        outcome_regressors = zeros(2,length(time_bin_centers));
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate selected regressors, set parameters
        
        task_regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        task_regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
                       
        task_t_shifts = { ...
            [0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [0,0.5]; ... % go cue
            [0,0.5]}; % outcome       
        
        % (old extended timings)
        %         task_t_shifts = {[0,0.5]; ... % stim
        %             [-0.5,1]; ... % move
        %             [-0.1,0.5]; ... % go cue
        %             [-0.5,1]}; % outcome

        % (to include stim x move)
        %         task_regressors = {stim_regressors;move_onset_regressors;move_onset_stim_regressors;go_cue_regressors;outcome_regressors};
        %         task_regressor_labels = {'Stim','Move','Stim x move','Go cue','Outcome'};
        %
        %         task_t_shifts = { ...
        %             [0,0.5]; ... % stim
        %             [-0.5,1]; ... % move
        %             [-0.5,1]; ... % stim x move
        %             [0,0.5]; ... % go cue
        %             [0,0.5]}; % outcome
        
        task_regressor_sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),task_t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;
        
        % Regression task -> MUA
        baseline = nanmean(reshape(event_aligned_mua(:,t < 0,:),[], ...
            size(event_aligned_mua,3))*raster_sample_rate,1)';
        activity = single(binned_spikes) - baseline;
        
        [mua_taskpred_k,mua_taskpred_long,mua_taskpred_expl_var,mua_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        mua_taskpred = ...
            interp1(time_bin_centers,mua_taskpred_long',t_peri_event)./raster_sample_rate;
        
        mua_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,mua_taskpred_reduced_long(:,:,x)', ...
            t_peri_event)./raster_sample_rate,permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        % Regression task -> MUA-ctxpred
        baseline = nanmean(reshape(event_aligned_mua_ctxpred(:,t < 0,:),[], ...
            size(event_aligned_mua_ctxpred,3))*raster_sample_rate,1)';
        activity = single(ctxpred_spikes) - baseline;
        
        [mua_ctxpred_taskpred_k,mua_ctxpred_taskpred_long,mua_ctxpred_taskpred_expl_var,mua_ctxpred_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        mua_ctxpred_taskpred = ...
            interp1(time_bin_centers,mua_ctxpred_taskpred_long',t_peri_event)./raster_sample_rate;
        
        mua_ctxpred_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,mua_ctxpred_taskpred_reduced_long(:,:,x)', ...
            t_peri_event)./raster_sample_rate,permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        % Regression task -> (master U, deconvolved) fluor
        event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
        fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);
        
        baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,:),[],size(event_aligned_V_deconv,3)))';
        activity = single(fVdf_deconv_resample_recast(use_components,:))-baseline;
        
        [fluor_taskpred_k,fluor_taskpred_long,fluor_taskpred_expl_var,fluor_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        fluor_taskpred = ...
            interp1(time_bin_centers,fluor_taskpred_long',t_peri_event);
        
        fluor_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,fluor_taskpred_reduced_long(:,:,x)', ...
            t_peri_event),permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        mua_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_taskpred_k;
        mua_taskpred_all{curr_animal,1}{curr_day,1} = mua_taskpred(use_trials,:,:,:);
        mua_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_taskpred_reduced(use_trials,:,:,:);
        mua_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = mua_taskpred_expl_var.total;
        mua_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = mua_taskpred_expl_var.partial;
        
        mua_ctxpred_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_k;
        mua_ctxpred_taskpred_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred(use_trials,:,:,:);
        mua_ctxpred_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_reduced(use_trials,:,:,:);
        mua_ctxpred_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_expl_var.total;
        mua_ctxpred_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_expl_var.partial;
        
        fluor_taskpred_k_all{curr_animal,1}{curr_day,1} = fluor_taskpred_k;
        fluor_taskpred_all{curr_animal,1}{curr_day,1} = fluor_taskpred(use_trials,:,:,:);
        fluor_taskpred_reduced_all{curr_animal,1}{curr_day,1} = fluor_taskpred_reduced(use_trials,:,:,:);
        fluor_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.total;
        fluor_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.partial;
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        movement_all{curr_animal,1}{curr_day,1} = event_aligned_movement(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars -except curr_animal animal protocol experiments load_parts curr_day ...
            n_aligned_depths regression_params animals t ...
            ...
            fluor_all ...
            mua_all ...
            mua_ctxpred_all ...
            ...
            mua_taskpred_k_all ...
            mua_taskpred_all ...
            mua_taskpred_reduced_all ...
            mua_taskpred_expl_var_total_all ...
            mua_taskpred_expl_var_partial_all ...
            ...
            mua_ctxpred_taskpred_k_all ...
            mua_ctxpred_taskpred_all ...
            mua_ctxpred_taskpred_reduced_all ...
            mua_ctxpred_taskpred_expl_var_total_all ...
            mua_ctxpred_taskpred_expl_var_partial_all ...
            ...
            fluor_taskpred_k_all ...
            fluor_taskpred_all ...
            fluor_taskpred_reduced_all ...
            fluor_taskpred_expl_var_total_all ...
            fluor_taskpred_expl_var_partial_all ...
            ...
            wheel_ctxpred_all ...
            ctx_str_k_all ...
            ctx_wheel_k_all ...
            wheel_all ...
            movement_all ...
            outcome_all ...
            D_all ...
            ...
            task_regressor_labels ...
            task_regressor_sample_shifts
        
    end
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_taskpred_expl_var_total_all ...
    mua_taskpred_expl_var_partial_all ...
    ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    mua_ctxpred_taskpred_expl_var_total_all ...
    mua_ctxpred_taskpred_expl_var_partial_all ...
    ...
    fluor_taskpred_k_all ...
    fluor_taskpred_all ...
    fluor_taskpred_reduced_all ...
    fluor_taskpred_expl_var_total_all ...
    fluor_taskpred_expl_var_partial_all ...
    ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    movement_all ...
    outcome_all ...
    D_all ...
    ...
    task_regressor_labels ...
    task_regressor_sample_shifts

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_ctxLonly'];
save([save_path filesep save_fn],'-v7.3');


%% (choiceworld trials: trying low-pass filtering str responses)
% (is the cortical estimation just because there's fast signals in str that
% aren't predictable from cortex?)

clear all
disp('Choiceworld trial activity')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);

mua_taskpred_k_all = cell(1,1);
mua_taskpred_all = cell(1,1);
mua_taskpred_reduced_all = cell(1,1);
mua_taskpred_expl_var_total_all = cell(1,1);
mua_taskpred_expl_var_partial_all = cell(1,1);

mua_ctxpred_taskpred_k_all = cell(1,1);
mua_ctxpred_taskpred_all = cell(1,1);
mua_ctxpred_taskpred_reduced_all = cell(1,1);
mua_ctxpred_taskpred_expl_var_total_all = cell(1,1);
mua_ctxpred_taskpred_expl_var_partial_all = cell(1,1);

fluor_taskpred_k_all = cell(1,1);
fluor_taskpred_all = cell(1,1);
fluor_taskpred_reduced_all = cell(1,1);
fluor_taskpred_expl_var_total_all = cell(1,1);
fluor_taskpred_expl_var_partial_all = cell(1,1);

wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
movement_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Get align times
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        
        %%% Trial-align cortex
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% Trial-align striatum
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            if isempty(curr_spikes)
                continue
            end
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regress cortex to striatum
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Deconvolve fluoresence
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        
        % Bin spikes across the experiment
        binned_spikes = nan(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            % Skip if no spikes at this depth
            if isempty(curr_spike_times)
                continue
            end
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        
        %%%%%%%%% TESTING: LOWPASS STR SPIKES
        
        lowpassCutoff = 10;
        [b100s, a100s] = butter(2, lowpassCutoff/((sample_rate)/2),'low');
        binned_spikes = filtfilt(b100s,a100s,binned_spikes')';
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        
        %%%%%%%%%%%%%        
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        % Regress cortex to striatum
        [ctx_str_k,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Recast the k's into the master U
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
        
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        ctxpred_spikes = (ctxpred_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,ctxpred_spikes',t_peri_event)./raster_sample_rate;
        
        % Regress cortex to wheel velocity/speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Recast the k's into the master U
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%% Trial-align wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Trial-align facecam movement
        event_aligned_movement = interp1(facecam_t(~isnan(facecam_t)), ...
            frame_movement(~isnan(facecam_t)),t_peri_event);
        
        %%% Trial-align outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_peri_event_bins(x,:)), ...
            find(trial_outcome == 1),'uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_peri_event_bins(x,:)), ...
            find(trial_outcome == -1),'uni',false))) > 0;
        
        % Pick trials to keep
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials)' & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials)' == -1;
        R_trials = signals_events.trialSideValues(1:n_trials)' == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = 3-(abs((trial_choice(use_trials)+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%% Regress task to cortex/striatum/cortex-predicted striatum
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Build regressors (only a subset of these are used)
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times))'.* ...
            signals_events.trialContrastValues(1:length(stimOn_times))';
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end
        
        % Stim move regressors (one for each stim when it starts to move)
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
            
        end
        
        % Stim center regressors (one for each stim when it's stopped during reward)
        unique_contrasts = unique(contrasts(contrasts > 0));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ...
                abs(stim_contrastsides) == unique_contrasts(curr_contrast));
            
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
            
        end
        
        % Move onset regressors (L/R)
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move onset x stim regressors (one for each contrast/side)
        move_onset_stim_time_absolute = arrayfun(@(curr_stim) ...
            arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & stim_contrastsides == unique_stim(curr_stim))), ...
            1:length(unique_stim),'uni',false);
        
        move_onset_stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            move_onset_stim_regressors(curr_stim,:) = ...
                histcounts(move_onset_stim_time_absolute{curr_stim},time_bins);
        end
        
        % Move ongoing regressors (L/R choice for duration of movement)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        
        move_stopped_t = 0.5;
        move_stopped_samples = round(sample_rate*move_stopped_t);
        wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.02), ...
            ones(1,move_stopped_samples),'full') > 0;
        wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);
        
        move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_L_absolute','uni',false));
        move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_R_absolute','uni',false));
        
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
        move_ongoing_regressors(2,move_ongoing_R_samples) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        % (for go cue only on late move trials)
%         go_cue_regressors = histcounts( ...
%             signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        % (for go cue with early/late move trials)
        go_cue_regressors = zeros(1,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        % (regressors for hit only)
%         outcome_regressors = histcounts(reward_t_timeline,time_bins);
        % (regressors for both hit and miss)
        outcome_regressors = zeros(2,length(time_bin_centers));
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate selected regressors, set parameters
        
        task_regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        task_regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
                       
        task_t_shifts = { ...
            [0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [0,0.5]; ... % go cue
            [0,0.5]}; % outcome       
        
        % (old extended timings)
        %         task_t_shifts = {[0,0.5]; ... % stim
        %             [-0.5,1]; ... % move
        %             [-0.1,0.5]; ... % go cue
        %             [-0.5,1]}; % outcome

        % (to include stim x move)
        %         task_regressors = {stim_regressors;move_onset_regressors;move_onset_stim_regressors;go_cue_regressors;outcome_regressors};
        %         task_regressor_labels = {'Stim','Move','Stim x move','Go cue','Outcome'};
        %
        %         task_t_shifts = { ...
        %             [0,0.5]; ... % stim
        %             [-0.5,1]; ... % move
        %             [-0.5,1]; ... % stim x move
        %             [0,0.5]; ... % go cue
        %             [0,0.5]}; % outcome
        
        task_regressor_sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),task_t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;
        
        % Regression task -> MUA
        baseline = nanmean(reshape(event_aligned_mua(:,t < 0,:),[], ...
            size(event_aligned_mua,3))*raster_sample_rate,1)';
        activity = single(binned_spikes) - baseline;
        
        [mua_taskpred_k,mua_taskpred_long,mua_taskpred_expl_var,mua_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        mua_taskpred = ...
            interp1(time_bin_centers,mua_taskpred_long',t_peri_event)./raster_sample_rate;
        
        mua_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,mua_taskpred_reduced_long(:,:,x)', ...
            t_peri_event)./raster_sample_rate,permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        % Regression task -> MUA-ctxpred
        baseline = nanmean(reshape(event_aligned_mua_ctxpred(:,t < 0,:),[], ...
            size(event_aligned_mua_ctxpred,3))*raster_sample_rate,1)';
        activity = single(ctxpred_spikes) - baseline;
        
        [mua_ctxpred_taskpred_k,mua_ctxpred_taskpred_long,mua_ctxpred_taskpred_expl_var,mua_ctxpred_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        mua_ctxpred_taskpred = ...
            interp1(time_bin_centers,mua_ctxpred_taskpred_long',t_peri_event)./raster_sample_rate;
        
        mua_ctxpred_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,mua_ctxpred_taskpred_reduced_long(:,:,x)', ...
            t_peri_event)./raster_sample_rate,permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        % Regression task -> (master U, deconvolved) fluor
        event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
        fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);
        
        baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,:),[],size(event_aligned_V_deconv,3)))';
        activity = single(fVdf_deconv_resample_recast(use_components,:))-baseline;
        
        [fluor_taskpred_k,fluor_taskpred_long,fluor_taskpred_expl_var,fluor_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        fluor_taskpred = ...
            interp1(time_bin_centers,fluor_taskpred_long',t_peri_event);
        
        fluor_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,fluor_taskpred_reduced_long(:,:,x)', ...
            t_peri_event),permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        mua_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_taskpred_k;
        mua_taskpred_all{curr_animal,1}{curr_day,1} = mua_taskpred(use_trials,:,:,:);
        mua_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_taskpred_reduced(use_trials,:,:,:);
        mua_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = mua_taskpred_expl_var.total;
        mua_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = mua_taskpred_expl_var.partial;
        
        mua_ctxpred_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_k;
        mua_ctxpred_taskpred_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred(use_trials,:,:,:);
        mua_ctxpred_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_reduced(use_trials,:,:,:);
        mua_ctxpred_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_expl_var.total;
        mua_ctxpred_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_expl_var.partial;
        
        fluor_taskpred_k_all{curr_animal,1}{curr_day,1} = fluor_taskpred_k;
        fluor_taskpred_all{curr_animal,1}{curr_day,1} = fluor_taskpred(use_trials,:,:,:);
        fluor_taskpred_reduced_all{curr_animal,1}{curr_day,1} = fluor_taskpred_reduced(use_trials,:,:,:);
        fluor_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.total;
        fluor_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.partial;
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        movement_all{curr_animal,1}{curr_day,1} = event_aligned_movement(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars -except curr_animal animal protocol experiments load_parts curr_day ...
            n_aligned_depths regression_params animals t ...
            ...
            fluor_all ...
            mua_all ...
            mua_ctxpred_all ...
            ...
            mua_taskpred_k_all ...
            mua_taskpred_all ...
            mua_taskpred_reduced_all ...
            mua_taskpred_expl_var_total_all ...
            mua_taskpred_expl_var_partial_all ...
            ...
            mua_ctxpred_taskpred_k_all ...
            mua_ctxpred_taskpred_all ...
            mua_ctxpred_taskpred_reduced_all ...
            mua_ctxpred_taskpred_expl_var_total_all ...
            mua_ctxpred_taskpred_expl_var_partial_all ...
            ...
            fluor_taskpred_k_all ...
            fluor_taskpred_all ...
            fluor_taskpred_reduced_all ...
            fluor_taskpred_expl_var_total_all ...
            fluor_taskpred_expl_var_partial_all ...
            ...
            wheel_ctxpred_all ...
            ctx_str_k_all ...
            ctx_wheel_k_all ...
            wheel_all ...
            movement_all ...
            outcome_all ...
            D_all ...
            ...
            task_regressor_labels ...
            task_regressor_sample_shifts
        
    end
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_taskpred_expl_var_total_all ...
    mua_taskpred_expl_var_partial_all ...
    ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    mua_ctxpred_taskpred_expl_var_total_all ...
    mua_ctxpred_taskpred_expl_var_partial_all ...
    ...
    fluor_taskpred_k_all ...
    fluor_taskpred_all ...
    fluor_taskpred_reduced_all ...
    fluor_taskpred_expl_var_total_all ...
    fluor_taskpred_expl_var_partial_all ...
    ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    movement_all ...
    outcome_all ...
    D_all ...
    ...
    task_regressor_labels ...
    task_regressor_sample_shifts

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_strfilt'];
save([save_path filesep save_fn],'-v7.3');


%% (choiceworld trials: trying regression from (ctx and task) -> str

clear all
disp('Choiceworld trial activity')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);

mua_taskpred_k_all = cell(1,1);
mua_taskpred_all = cell(1,1);
mua_taskpred_reduced_all = cell(1,1);
mua_taskpred_expl_var_total_all = cell(1,1);
mua_taskpred_expl_var_partial_all = cell(1,1);

mua_ctxpred_taskpred_k_all = cell(1,1);
mua_ctxpred_taskpred_all = cell(1,1);
mua_ctxpred_taskpred_reduced_all = cell(1,1);
mua_ctxpred_taskpred_expl_var_total_all = cell(1,1);
mua_ctxpred_taskpred_expl_var_partial_all = cell(1,1);

fluor_taskpred_k_all = cell(1,1);
fluor_taskpred_all = cell(1,1);
fluor_taskpred_reduced_all = cell(1,1);
fluor_taskpred_expl_var_total_all = cell(1,1);
fluor_taskpred_expl_var_partial_all = cell(1,1);

mua_ctxtaskpred_all = cell(1,1);

wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
movement_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Get align times
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        
        %%% Trial-align cortex
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% Trial-align striatum
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            if isempty(curr_spikes)
                continue
            end
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regress cortex to striatum
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Deconvolve fluoresence
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        
        % Bin spikes across the experiment
        binned_spikes = nan(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            % Skip if no spikes at this depth
            if isempty(curr_spike_times)
                continue
            end
            binned_spikes(curr_depth,:) = histcounts(curr_spike_times,time_bins);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                ctx_lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        % Regress cortex to striatum
        [ctx_str_k,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,ctx_lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Recast the k's into the master U
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
        
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        ctxpred_spikes = (ctxpred_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,ctxpred_spikes',t_peri_event)./raster_sample_rate;
        
        % Regress cortex to wheel velocity/speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,ctx_lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Recast the k's into the master U
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%% Trial-align wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Trial-align facecam movement
        event_aligned_movement = interp1(facecam_t(~isnan(facecam_t)), ...
            frame_movement(~isnan(facecam_t)),t_peri_event);
        
        %%% Trial-align outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_peri_event_bins(x,:)), ...
            find(trial_outcome == 1),'uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_peri_event_bins(x,:)), ...
            find(trial_outcome == -1),'uni',false))) > 0;
        
        % Pick trials to keep
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials)' & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials)' == -1;
        R_trials = signals_events.trialSideValues(1:n_trials)' == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = 3-(abs((trial_choice(use_trials)+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%% Regress task to cortex/striatum/cortex-predicted striatum
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Build regressors (only a subset of these are used)
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times))'.* ...
            signals_events.trialContrastValues(1:length(stimOn_times))';
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end
        
        % Stim move regressors (one for each stim when it starts to move)
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
            
        end
        
        % Stim center regressors (one for each stim when it's stopped during reward)
        unique_contrasts = unique(contrasts(contrasts > 0));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ...
                abs(stim_contrastsides) == unique_contrasts(curr_contrast));
            
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
            
        end
        
        % Move onset regressors (L/R)
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move onset x stim regressors (one for each contrast/side)
        move_onset_stim_time_absolute = arrayfun(@(curr_stim) ...
            arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & stim_contrastsides == unique_stim(curr_stim))), ...
            1:length(unique_stim),'uni',false);
        
        move_onset_stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            move_onset_stim_regressors(curr_stim,:) = ...
                histcounts(move_onset_stim_time_absolute{curr_stim},time_bins);
        end
        
        % Move ongoing regressors (L/R choice for duration of movement)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        
        move_stopped_t = 0.5;
        move_stopped_samples = round(sample_rate*move_stopped_t);
        wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.02), ...
            ones(1,move_stopped_samples),'full') > 0;
        wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);
        
        move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_L_absolute','uni',false));
        move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_R_absolute','uni',false));
        
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
        move_ongoing_regressors(2,move_ongoing_R_samples) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        % (for go cue only on late move trials)
%         go_cue_regressors = histcounts( ...
%             signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        % (for go cue with early/late move trials)
        go_cue_regressors = zeros(1,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        % (regressors for hit only)
%         outcome_regressors = histcounts(reward_t_timeline,time_bins);
        % (regressors for both hit and miss)
        outcome_regressors = zeros(2,length(time_bin_centers));
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate selected regressors, set parameters
        
        task_regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        task_regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
                       
        task_t_shifts = { ...
            [0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [0,0.5]; ... % go cue
            [0,0.5]}; % outcome       
        
        % (old extended timings)
        %         task_t_shifts = {[0,0.5]; ... % stim
        %             [-0.5,1]; ... % move
        %             [-0.1,0.5]; ... % go cue
        %             [-0.5,1]}; % outcome

        % (to include stim x move)
        %         task_regressors = {stim_regressors;move_onset_regressors;move_onset_stim_regressors;go_cue_regressors;outcome_regressors};
        %         task_regressor_labels = {'Stim','Move','Stim x move','Go cue','Outcome'};
        %
        %         task_t_shifts = { ...
        %             [0,0.5]; ... % stim
        %             [-0.5,1]; ... % move
        %             [-0.5,1]; ... % stim x move
        %             [0,0.5]; ... % go cue
        %             [0,0.5]}; % outcome
        
        task_regressor_sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),task_t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;
        
        % Regression task -> MUA
        baseline = nanmean(reshape(event_aligned_mua(:,t < 0,:),[], ...
            size(event_aligned_mua,3))*raster_sample_rate,1)';
        activity = single(binned_spikes) - baseline;
        
        [mua_taskpred_k,mua_taskpred_long,mua_taskpred_expl_var,mua_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        mua_taskpred = ...
            interp1(time_bin_centers,mua_taskpred_long',t_peri_event)./raster_sample_rate;
        
        mua_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,mua_taskpred_reduced_long(:,:,x)', ...
            t_peri_event)./raster_sample_rate,permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        % Regression task -> MUA-ctxpred
        baseline = nanmean(reshape(event_aligned_mua_ctxpred(:,t < 0,:),[], ...
            size(event_aligned_mua_ctxpred,3))*raster_sample_rate,1)';
        activity = single(ctxpred_spikes) - baseline;
        
        [mua_ctxpred_taskpred_k,mua_ctxpred_taskpred_long,mua_ctxpred_taskpred_expl_var,mua_ctxpred_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        mua_ctxpred_taskpred = ...
            interp1(time_bin_centers,mua_ctxpred_taskpred_long',t_peri_event)./raster_sample_rate;
        
        mua_ctxpred_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,mua_ctxpred_taskpred_reduced_long(:,:,x)', ...
            t_peri_event)./raster_sample_rate,permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        % Regression task -> (master U, deconvolved) fluor
        event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
        fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);
        
        baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,:),[],size(event_aligned_V_deconv,3)))';
        activity = single(fVdf_deconv_resample_recast(use_components,:))-baseline;
        
        [fluor_taskpred_k,fluor_taskpred_long,fluor_taskpred_expl_var,fluor_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        fluor_taskpred = ...
            interp1(time_bin_centers,fluor_taskpred_long',t_peri_event);
        
        fluor_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,fluor_taskpred_reduced_long(:,:,x)', ...
            t_peri_event),permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        %%%% TESTING: Regression (ctx and task) -> str
        fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master); 
        
        ctx_activity = single(fVdf_deconv_resample_recast(regression_params.use_svs,:));
        str_activity = single(binned_spikes);
        str_activity_std = str_activity./nanstd(str_activity,[],2);
        
        task_regressors_single = cellfun(@single,task_regressors,'uni',false);
        
        [ctx_task_k,mua_ctxtaskpred_std_long] = AP_regresskernel( ...
            [{ctx_activity};task_regressors_single], ... & Task and ctx regressors
            str_activity_std, ... % Str activity (std)
            [{kernel_frames};task_regressor_sample_shifts], ... & time shifts for ctx and task
            [ctx_lambda;zeros(size(task_regressors))], ... % lambdas for ctx and zeros for task
            [false,false], ... % don't auto zscore anything
            5,1,1); % cvfold, return and use constant
        
        % re-scale spikes 
        mua_ctxtaskpred_long = (mua_ctxtaskpred_std_long - squeeze(ctx_task_k{end})).* ...
            nanstd(str_activity,[],2) + ...
            nanstd(str_activity,[],2).*squeeze(ctx_task_k{end});
        
        mua_ctxtaskpred = ...
            interp1(time_bin_centers,mua_ctxtaskpred_long',t_peri_event)./raster_sample_rate;    
        
        %%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        mua_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_taskpred_k;
        mua_taskpred_all{curr_animal,1}{curr_day,1} = mua_taskpred(use_trials,:,:,:);
        mua_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_taskpred_reduced(use_trials,:,:,:);
        mua_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = mua_taskpred_expl_var.total;
        mua_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = mua_taskpred_expl_var.partial;
        
        mua_ctxpred_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_k;
        mua_ctxpred_taskpred_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred(use_trials,:,:,:);
        mua_ctxpred_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_reduced(use_trials,:,:,:);
        mua_ctxpred_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_expl_var.total;
        mua_ctxpred_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_expl_var.partial;
        
        fluor_taskpred_k_all{curr_animal,1}{curr_day,1} = fluor_taskpred_k;
        fluor_taskpred_all{curr_animal,1}{curr_day,1} = fluor_taskpred(use_trials,:,:,:);
        fluor_taskpred_reduced_all{curr_animal,1}{curr_day,1} = fluor_taskpred_reduced(use_trials,:,:,:);
        fluor_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.total;
        fluor_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.partial;
        
        mua_ctxtaskpred_all{curr_animal,1}{curr_day,1} = mua_ctxtaskpred(use_trials,:,:);
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        movement_all{curr_animal,1}{curr_day,1} = event_aligned_movement(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars -except curr_animal animal protocol experiments load_parts curr_day ...
            n_aligned_depths regression_params animals t ...
            ...
            fluor_all ...
            mua_all ...
            mua_ctxpred_all ...
            ...
            mua_taskpred_k_all ...
            mua_taskpred_all ...
            mua_taskpred_reduced_all ...
            mua_taskpred_expl_var_total_all ...
            mua_taskpred_expl_var_partial_all ...
            ...
            mua_ctxpred_taskpred_k_all ...
            mua_ctxpred_taskpred_all ...
            mua_ctxpred_taskpred_reduced_all ...
            mua_ctxpred_taskpred_expl_var_total_all ...
            mua_ctxpred_taskpred_expl_var_partial_all ...
            ...
            fluor_taskpred_k_all ...
            fluor_taskpred_all ...
            fluor_taskpred_reduced_all ...
            fluor_taskpred_expl_var_total_all ...
            fluor_taskpred_expl_var_partial_all ...
            ...
            mua_ctxtaskpred_all ...
            ...
            wheel_ctxpred_all ...
            ctx_str_k_all ...
            ctx_wheel_k_all ...
            wheel_all ...
            movement_all ...
            outcome_all ...
            D_all ...
            ...
            task_regressor_labels ...
            task_regressor_sample_shifts
        
    end
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_taskpred_expl_var_total_all ...
    mua_taskpred_expl_var_partial_all ...
    ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    mua_ctxpred_taskpred_expl_var_total_all ...
    mua_ctxpred_taskpred_expl_var_partial_all ...
    ...
    fluor_taskpred_k_all ...
    fluor_taskpred_all ...
    fluor_taskpred_reduced_all ...
    fluor_taskpred_expl_var_total_all ...
    fluor_taskpred_expl_var_partial_all ...
    ...
    mua_ctxtaskpred_all ...
    ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    movement_all ...
    outcome_all ...
    D_all ...
    ...
    task_regressor_labels ...
    task_regressor_sample_shifts

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_ctxtaskpred'];
save([save_path filesep save_fn],'-v7.3');


%% (choiceworld trial activity: smoothed striatum)

clear all
disp('Choiceworld trial activity')

n_aligned_depths = 4;

% Regression parameters
regression_params.use_svs = 1:50;
regression_params.skip_seconds = 20;
regression_params.upsample_factor = 1;
regression_params.kernel_t = [-0.5,0.5];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;

animals = {'AP024','AP025','AP026','AP027','AP028','AP029'};

% Initialize all saved variables for indexing
fluor_all = cell(1,1);
mua_all = cell(1,1);
mua_ctxpred_all = cell(1,1);

mua_taskpred_k_all = cell(1,1);
mua_taskpred_all = cell(1,1);
mua_taskpred_reduced_all = cell(1,1);
mua_taskpred_expl_var_total_all = cell(1,1);
mua_taskpred_expl_var_partial_all = cell(1,1);

mua_ctxpred_taskpred_k_all = cell(1,1);
mua_ctxpred_taskpred_all = cell(1,1);
mua_ctxpred_taskpred_reduced_all = cell(1,1);
mua_ctxpred_taskpred_expl_var_total_all = cell(1,1);
mua_ctxpred_taskpred_expl_var_partial_all = cell(1,1);

fluor_taskpred_k_all = cell(1,1);
fluor_taskpred_all = cell(1,1);
fluor_taskpred_reduced_all = cell(1,1);
fluor_taskpred_expl_var_total_all = cell(1,1);
fluor_taskpred_expl_var_partial_all = cell(1,1);

wheel_ctxpred_all = cell(1,1);
ctx_str_k_all = cell(1,1);
ctx_wheel_k_all = cell(1,1);
wheel_all = cell(1,1);
movement_all = cell(1,1);
outcome_all = cell(1,1);
D_all = cell(1,1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    
    experiments = experiments([experiments.imaging] & [experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        smooth_factor = 5;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment;
        
        % Load experiment
        str_align = 'kernel';
        AP_load_experiment;
        
        % Prepare fluorescence
        % Convert U to master U
        load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master.mat');
        Udf_aligned = single(AP_align_widefield(animal,day,Udf));
        fVdf_recast = ChangeU(Udf_aligned,fVdf,U_master);
        
        % Set components to keep
        use_components = 1:200;
        
        % (aligned striatum depths)
        n_depths = n_aligned_depths;
        depth_group = aligned_str_depth_group;
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        upsample_factor = 1;
        raster_sample_rate = 1/(framerate*upsample_factor);
        t = raster_window(1):raster_sample_rate:raster_window(2);
        
        % Get align times
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = bsxfun(@plus,use_align,t);
        t_peri_event_bins = [t_peri_event-raster_sample_rate/2,t_peri_event(:,end)+raster_sample_rate/2];
        
        %%% Trial-align cortex
        event_aligned_V = ...
            interp1(frame_t,fVdf_recast(use_components,:)',t_peri_event);
        
        %%% Trial-align striatum
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            % (for only msns in depth group)
            %                 curr_spikes = spike_times_timeline(depth_group == curr_depth & ...
            %                     ismember(spike_templates,find(msn)));
            
            if isempty(curr_spikes)
                continue
            end
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                smooth(histcounts(curr_spikes,t_peri_event_bins(x,:)),smooth_factor)', ...
                [1:size(t_peri_event,1)]','uni',false))./raster_sample_rate;
        end
        
        %%% Regress cortex to striatum
        
        % Get time points to bin
        sample_rate = framerate*regression_params.upsample_factor;
        time_bins = frame_t(find(frame_t > ...
            regression_params.skip_seconds,1)):1/sample_rate: ...
            frame_t(find(frame_t-frame_t(end) < ...
            -regression_params.skip_seconds,1,'last'));
        time_bin_centers = time_bins(1:end-1) + diff(time_bins)/2;
        
        % Deconvolve fluoresence
        fVdf_deconv = AP_deconv_wf(fVdf);
        fVdf_deconv(isnan(fVdf_deconv)) = 0;
        fVdf_deconv_resample = interp1(frame_t,fVdf_deconv',time_bin_centers)';
        
        % Bin spikes across the experiment
        binned_spikes = nan(n_depths,length(time_bin_centers));
        for curr_depth = 1:n_depths
            curr_spike_times = spike_times_timeline(depth_group == curr_depth);
            % Skip if no spikes at this depth
            if isempty(curr_spike_times)
                continue
            end
            binned_spikes(curr_depth,:) = smooth(histcounts(curr_spike_times,time_bins),smooth_factor);
        end
        binned_spikes_std = binned_spikes./nanstd(binned_spikes,[],2);
        
        % Load lambda from previously estimated and saved
        lambda_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\ephys_processing\ctx-str_lambda';
        load(lambda_fn);
        curr_animal_idx = strcmp(animal,{ctx_str_lambda.animal});
        if any(curr_animal_idx)
            curr_day_idx = strcmp(day,ctx_str_lambda(curr_animal_idx).day);
            if any(curr_day_idx)
                lambda = ctx_str_lambda(curr_animal_idx).best_lambda(curr_day_idx);
            end
        end
        
        kernel_frames = round(regression_params.kernel_t(1)*sample_rate): ...
            round(regression_params.kernel_t(2)*sample_rate);
        
        % Regress cortex to striatum
        [ctx_str_k,ctxpred_spikes_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            binned_spikes_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            true,regression_params.use_constant);
        
        % Recast the k's into the master U
        ctx_str_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_str_k{1},size(ctx_str_k{1},1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_str_k{1}));
        
        % Re-scale the prediction (subtract offset, multiply, add scaled offset)
        ctxpred_spikes = (ctxpred_spikes_std - squeeze(ctx_str_k{end})).* ...
            nanstd(binned_spikes,[],2) + ...
            nanstd(binned_spikes,[],2).*squeeze(ctx_str_k{end});
        
        event_aligned_mua_ctxpred = ...
            interp1(time_bin_centers,ctxpred_spikes',t_peri_event)./raster_sample_rate;
        
        % Regress cortex to wheel velocity/speed
        wheel_velocity_resample = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        wheel_velspeed_resample = [wheel_velocity_resample;abs(wheel_velocity_resample)];
        wheel_velspeed_resample_std = wheel_velspeed_resample./std(wheel_velocity_resample);
        
        [ctx_wheel_k,predicted_wheel_velspeed_std,explained_var] = ...
            AP_regresskernel(fVdf_deconv_resample(regression_params.use_svs,:), ...
            wheel_velspeed_resample_std,kernel_frames,lambda, ...
            regression_params.zs,regression_params.cvfold, ...
            false,false);
        
        predicted_wheel_velspeed = predicted_wheel_velspeed_std.* ...
            std(wheel_velocity_resample);
        
        % Recast the k's into the master U
        ctx_wheel_k_recast = reshape(ChangeU(Udf_aligned(:,:,regression_params.use_svs), ...
            reshape(ctx_wheel_k,size(ctx_wheel_k,1),[]),U_master(:,:,regression_params.use_svs)), ...
            size(ctx_wheel_k));
        
        event_aligned_wheel_ctxpred = ...
            interp1(time_bin_centers,predicted_wheel_velspeed(1,:)',t_peri_event);
        
        %%% Trial-align wheel velocity
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            wheel_velocity,t_peri_event);
        
        %%% Trial-align facecam movement
        event_aligned_movement = interp1(facecam_t(~isnan(facecam_t)), ...
            frame_movement(~isnan(facecam_t)),t_peri_event);
        
        %%% Trial-align outcome (reward page 1, punish page 2)
        % (note incorrect outcome imprecise from signals, but looks good)
        event_aligned_outcome = zeros(size(t_peri_event,1),size(t_peri_event,2),2);
        
        event_aligned_outcome(trial_outcome == 1,:,1) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(reward_t_timeline,t_peri_event_bins(x,:)), ...
            find(trial_outcome == 1),'uni',false))) > 0;
        
        event_aligned_outcome(trial_outcome == -1,:,2) = ...
            (cell2mat(arrayfun(@(x) ...
            histcounts(signals_events.responseTimes,t_peri_event_bins(x,:)), ...
            find(trial_outcome == -1),'uni',false))) > 0;
        
        % Pick trials to keep
        use_trials = ...
            trial_outcome ~= 0 & ...
            ~signals_events.repeatTrialValues(1:n_trials)' & ...
            stim_to_feedback < 1.5;
        
        % Get behavioural data
        D = struct;
        D.stimulus = zeros(sum(use_trials),2);
        
        L_trials = signals_events.trialSideValues(1:n_trials)' == -1;
        R_trials = signals_events.trialSideValues(1:n_trials)' == 1;
        
        D.stimulus(L_trials(use_trials),1) = signals_events.trialContrastValues(L_trials & use_trials);
        D.stimulus(R_trials(use_trials),2) = signals_events.trialContrastValues(R_trials & use_trials);
        
        D.response = 3-(abs((trial_choice(use_trials)+1)/2)+1);
        D.repeatNum = ones(sum(use_trials),1);
        
        D.outcome = reshape(trial_outcome(use_trials),[],1);
        
        %%% Regress task to cortex/striatum/cortex-predicted striatum
        
        % Get reaction time for building regressors
        [move_trial,move_idx] = max(abs(event_aligned_wheel) > 0.02,[],2);
        move_idx(~move_trial) = NaN;
        move_t = nan(size(move_idx));
        move_t(~isnan(move_idx) & move_trial) = t(move_idx(~isnan(move_idx) & move_trial))';
        
        % Build regressors (only a subset of these are used)
        
        % Stim regressors
        unique_stim = unique(contrasts(contrasts > 0).*sides');
        stim_contrastsides = ...
            signals_events.trialSideValues(1:length(stimOn_times))'.* ...
            signals_events.trialContrastValues(1:length(stimOn_times))';
        
        stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            curr_stim_times = stimOn_times(stim_contrastsides == unique_stim(curr_stim));
            stim_regressors(curr_stim,:) = histcounts(curr_stim_times,time_bins);
        end
        
        % Stim move regressors (one for each stim when it starts to move)
        stim_move_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            
            % (find the first photodiode flip after the stim azimuth has
            % moved past a threshold)
            
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) ~= 0 & ...
                stim_contrastsides == unique_stim(curr_stim));
            
            azimuth_move_threshold = 5; % degrees to consider stim moved
            stim_move_times_signals = ...
                signals_events.stimAzimuthTimes( ...
                abs(signals_events.stimAzimuthValues - 90) > azimuth_move_threshold);
            curr_stim_move_times_signals = arrayfun(@(x) ...
                stim_move_times_signals(find(stim_move_times_signals > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_stim_move_times_photodiode = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times > ...
                curr_stim_move_times_signals(x),1)),1:length(curr_stim_move_times_signals));
            
            stim_move_regressors(curr_stim,:) = histcounts(curr_stim_move_times_photodiode,time_bins);
            
        end
        
        % Stim center regressors (one for each stim when it's stopped during reward)
        unique_contrasts = unique(contrasts(contrasts > 0));
        
        stim_center_regressors = zeros(length(unique_contrasts),length(time_bin_centers));
        for curr_contrast = 1:length(unique_contrasts)
            
            % (find the last photodiode flip before the reward)
            curr_stimOn_times = stimOn_times(trial_outcome(1:length(stimOn_times)) == 1 & ...
                abs(stim_contrastsides) == unique_contrasts(curr_contrast));
            
            curr_reward_times = arrayfun(@(x) ...
                reward_t_timeline(find(reward_t_timeline > ...
                curr_stimOn_times(x),1)),1:length(curr_stimOn_times));
            
            curr_prereward_photodiode_times = arrayfun(@(x) ...
                photodiode_flip_times(find(photodiode_flip_times < ...
                curr_reward_times(x),1,'last')),1:length(curr_reward_times));
            
            stim_center_regressors(curr_contrast,:) = histcounts(curr_prereward_photodiode_times,time_bins);
            
        end
        
        % Move onset regressors (L/R)
        move_time_L_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == -1));
        move_time_R_absolute = arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & trial_choice(1:length(stimOn_times)) == 1));
        
        move_onset_regressors = zeros(2,length(time_bin_centers));
        move_onset_regressors(1,:) = histcounts(move_time_L_absolute,time_bins);
        move_onset_regressors(2,:) = histcounts(move_time_R_absolute,time_bins);
        
        % Move onset x stim regressors (one for each contrast/side)
        move_onset_stim_time_absolute = arrayfun(@(curr_stim) ...
            arrayfun(@(x) t_peri_event(x,move_idx(x)), ...
            find(~isnan(move_idx) & stim_contrastsides == unique_stim(curr_stim))), ...
            1:length(unique_stim),'uni',false);
        
        move_onset_stim_regressors = zeros(length(unique_stim),length(time_bin_centers));
        for curr_stim = 1:length(unique_stim)
            move_onset_stim_regressors(curr_stim,:) = ...
                histcounts(move_onset_stim_time_absolute{curr_stim},time_bins);
        end
        
        % Move ongoing regressors (L/R choice for duration of movement)
        wheel_velocity_interp = interp1(Timeline.rawDAQTimestamps,wheel_velocity,time_bin_centers);
        
        move_stopped_t = 0.5;
        move_stopped_samples = round(sample_rate*move_stopped_t);
        wheel_moving_conv = convn((abs(wheel_velocity_interp) > 0.02), ...
            ones(1,move_stopped_samples),'full') > 0;
        wheel_moving = wheel_moving_conv(end-length(time_bin_centers)+1:end);
        
        move_ongoing_L_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_L_absolute','uni',false));
        move_ongoing_R_samples = cell2mat(arrayfun(@(x) ...
            find(time_bin_centers > x): ...
            find(time_bin_centers > x & ~wheel_moving,1), ...
            move_time_R_absolute','uni',false));
        
        move_ongoing_regressors = zeros(2,length(time_bin_centers));
        move_ongoing_regressors(1,move_ongoing_L_samples) = 1;
        move_ongoing_regressors(2,move_ongoing_R_samples) = 1;
        
        % Go cue regressors - separate for early/late move
        % (using signals timing - not precise but looks good)
        % (for go cue only on late move trials)
%         go_cue_regressors = histcounts( ...
%             signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        % (for go cue with early/late move trials)
        go_cue_regressors = zeros(1,length(time_bin_centers));
        go_cue_regressors(1,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t <= 0.5),time_bins);
        go_cue_regressors(2,:) = histcounts( ...
            signals_events.interactiveOnTimes(move_t > 0.5),time_bins);
        
        % Outcome regressors
        % (using signals timing - not precise but looks good)
        % (regressors for hit only)
%         outcome_regressors = histcounts(reward_t_timeline,time_bins);
        % (regressors for both hit and miss)
        outcome_regressors = zeros(2,length(time_bin_centers));
        outcome_regressors(1,:) = histcounts( ...
            reward_t_timeline,time_bins);
        outcome_regressors(2,:) = histcounts( ...
            signals_events.responseTimes(trial_outcome == -1),time_bins);
        
        % Concatenate selected regressors, set parameters
        
        task_regressors = {stim_regressors;move_onset_regressors;go_cue_regressors;outcome_regressors};
        task_regressor_labels = {'Stim','Move onset','Go cue','Outcome'};
                       
        task_t_shifts = { ...
            [0,0.5]; ... % stim
            [-0.5,1]; ... % move
            [0,0.5]; ... % go cue
            [0,0.5]}; % outcome       
        
        % (old extended timings)
        %         task_t_shifts = {[0,0.5]; ... % stim
        %             [-0.5,1]; ... % move
        %             [-0.1,0.5]; ... % go cue
        %             [-0.5,1]}; % outcome

        % (to include stim x move)
        %         task_regressors = {stim_regressors;move_onset_regressors;move_onset_stim_regressors;go_cue_regressors;outcome_regressors};
        %         task_regressor_labels = {'Stim','Move','Stim x move','Go cue','Outcome'};
        %
        %         task_t_shifts = { ...
        %             [0,0.5]; ... % stim
        %             [-0.5,1]; ... % move
        %             [-0.5,1]; ... % stim x move
        %             [0,0.5]; ... % go cue
        %             [0,0.5]}; % outcome
        
        task_regressor_sample_shifts = cellfun(@(x) round(x(1)*(sample_rate)): ...
            round(x(2)*(sample_rate)),task_t_shifts,'uni',false);
        lambda = 0;
        zs = [false,false];
        cvfold = 5;
        use_constant = false;
        return_constant = false;
        
        % Regression task -> MUA
        baseline = nanmean(reshape(event_aligned_mua(:,t < 0,:),[], ...
            size(event_aligned_mua,3))*raster_sample_rate,1)';
        activity = single(binned_spikes) - baseline;
        
        [mua_taskpred_k,mua_taskpred_long,mua_taskpred_expl_var,mua_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        mua_taskpred = ...
            interp1(time_bin_centers,mua_taskpred_long',t_peri_event)./raster_sample_rate;
        
        mua_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,mua_taskpred_reduced_long(:,:,x)', ...
            t_peri_event)./raster_sample_rate,permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        % Regression task -> MUA-ctxpred
        baseline = nanmean(reshape(event_aligned_mua_ctxpred(:,t < 0,:),[], ...
            size(event_aligned_mua_ctxpred,3))*raster_sample_rate,1)';
        activity = single(ctxpred_spikes) - baseline;
        
        [mua_ctxpred_taskpred_k,mua_ctxpred_taskpred_long,mua_ctxpred_taskpred_expl_var,mua_ctxpred_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        mua_ctxpred_taskpred = ...
            interp1(time_bin_centers,mua_ctxpred_taskpred_long',t_peri_event)./raster_sample_rate;
        
        mua_ctxpred_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,mua_ctxpred_taskpred_reduced_long(:,:,x)', ...
            t_peri_event)./raster_sample_rate,permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        % Regression task -> (master U, deconvolved) fluor
        event_aligned_V_deconv = AP_deconv_wf(event_aligned_V);
        fVdf_deconv_resample_recast = ChangeU(Udf_aligned,fVdf_deconv_resample,U_master);
        
        baseline = nanmean(reshape(event_aligned_V_deconv(:,t < 0,:),[],size(event_aligned_V_deconv,3)))';
        activity = single(fVdf_deconv_resample_recast(use_components,:))-baseline;
        
        [fluor_taskpred_k,fluor_taskpred_long,fluor_taskpred_expl_var,fluor_taskpred_reduced_long] = ...
            AP_regresskernel(task_regressors,activity,task_regressor_sample_shifts, ...
            lambda,zs,cvfold,return_constant,use_constant);
        
        fluor_taskpred = ...
            interp1(time_bin_centers,fluor_taskpred_long',t_peri_event);
        
        fluor_taskpred_reduced = cell2mat(arrayfun(@(x) ...
            interp1(time_bin_centers,fluor_taskpred_reduced_long(:,:,x)', ...
            t_peri_event),permute(1:length(task_regressors),[1,3,4,2]),'uni',false));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Store everything
        fluor_all{curr_animal,1}{curr_day,1} = event_aligned_V(use_trials,:,:,:);
        mua_all{curr_animal,1}{curr_day,1} = event_aligned_mua(use_trials,:,:,:);
        
        ctx_str_k_all{curr_animal,1}{curr_day,1} = ctx_str_k_recast;
        mua_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_mua_ctxpred(use_trials,:,:,:);
        
        mua_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_taskpred_k;
        mua_taskpred_all{curr_animal,1}{curr_day,1} = mua_taskpred(use_trials,:,:,:);
        mua_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_taskpred_reduced(use_trials,:,:,:);
        mua_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = mua_taskpred_expl_var.total;
        mua_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = mua_taskpred_expl_var.partial;
        
        mua_ctxpred_taskpred_k_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_k;
        mua_ctxpred_taskpred_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred(use_trials,:,:,:);
        mua_ctxpred_taskpred_reduced_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_reduced(use_trials,:,:,:);
        mua_ctxpred_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_expl_var.total;
        mua_ctxpred_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = mua_ctxpred_taskpred_expl_var.partial;
        
        fluor_taskpred_k_all{curr_animal,1}{curr_day,1} = fluor_taskpred_k;
        fluor_taskpred_all{curr_animal,1}{curr_day,1} = fluor_taskpred(use_trials,:,:,:);
        fluor_taskpred_reduced_all{curr_animal,1}{curr_day,1} = fluor_taskpred_reduced(use_trials,:,:,:);
        fluor_taskpred_expl_var_total_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.total;
        fluor_taskpred_expl_var_partial_all{curr_animal,1}{curr_day,1} = fluor_taskpred_expl_var.partial;
        
        wheel_all{curr_animal,1}{curr_day,1} = event_aligned_wheel(use_trials,:,:);
        movement_all{curr_animal,1}{curr_day,1} = event_aligned_movement(use_trials,:,:);
        
        ctx_wheel_k_all{curr_animal,1}{curr_day,1} = ctx_wheel_k_recast;
        wheel_ctxpred_all{curr_animal,1}{curr_day,1} = event_aligned_wheel_ctxpred(use_trials,:,:);
        
        outcome_all{curr_animal,1}{curr_day,1} = event_aligned_outcome(use_trials,:,:);
        D_all{curr_animal,1}{curr_day,1} = D;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars -except curr_animal animal protocol experiments load_parts curr_day ...
            n_aligned_depths regression_params animals t ...
            ...
            fluor_all ...
            mua_all ...
            mua_ctxpred_all ...
            ...
            mua_taskpred_k_all ...
            mua_taskpred_all ...
            mua_taskpred_reduced_all ...
            mua_taskpred_expl_var_total_all ...
            mua_taskpred_expl_var_partial_all ...
            ...
            mua_ctxpred_taskpred_k_all ...
            mua_ctxpred_taskpred_all ...
            mua_ctxpred_taskpred_reduced_all ...
            mua_ctxpred_taskpred_expl_var_total_all ...
            mua_ctxpred_taskpred_expl_var_partial_all ...
            ...
            fluor_taskpred_k_all ...
            fluor_taskpred_all ...
            fluor_taskpred_reduced_all ...
            fluor_taskpred_expl_var_total_all ...
            fluor_taskpred_expl_var_partial_all ...
            ...
            wheel_ctxpred_all ...
            ctx_str_k_all ...
            ctx_wheel_k_all ...
            wheel_all ...
            movement_all ...
            outcome_all ...
            D_all ...
            ...
            task_regressor_labels ...
            task_regressor_sample_shifts
        
    end
end

clearvars -except ...
    n_aligned_depths regression_params animals t ...
    ...
    fluor_all ...
    mua_all ...
    mua_ctxpred_all ...
    ...
    mua_taskpred_k_all ...
    mua_taskpred_all ...
    mua_taskpred_reduced_all ...
    mua_taskpred_expl_var_total_all ...
    mua_taskpred_expl_var_partial_all ...
    ...
    mua_ctxpred_taskpred_k_all ...
    mua_ctxpred_taskpred_all ...
    mua_ctxpred_taskpred_reduced_all ...
    mua_ctxpred_taskpred_expl_var_total_all ...
    mua_ctxpred_taskpred_expl_var_partial_all ...
    ...
    fluor_taskpred_k_all ...
    fluor_taskpred_all ...
    fluor_taskpred_reduced_all ...
    fluor_taskpred_expl_var_total_all ...
    fluor_taskpred_expl_var_partial_all ...
    ...
    wheel_ctxpred_all ...
    ctx_str_k_all ...
    ctx_wheel_k_all ...
    wheel_all ...
    movement_all ...
    outcome_all ...
    D_all ...
    ...
    task_regressor_labels ...
    task_regressor_sample_shifts

disp('Finished loading all')

save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data';
save_fn = ['trial_activity_choiceworld_strsmooth'];
save([save_path filesep save_fn],'-v7.3');


%% [ctx L, ctxR] \ str

% Get flipped kernel ROIs and fluorescence (do once, not in load)
% kernel_roiR = AP_reflect_widefield(kernel_roi);
% 
% fluor_kernelroiR_deconv = permute(reshape( ...
%     AP_svd_roi(U_master(:,:,1:n_vs), ...
%     reshape(permute(fluor_allcat_deconv,[3,2,1]),n_vs,[]),[],[],kernel_roiR), ...
%     size(kernel_roiR,3),[],size(fluor_allcat_deconv,1)),[3,2,1]);
% 
% fluor_kernelroiR_deconv_std = ...
%     cell2mat(cellfun(@(x) repmat(nanstd(reshape(x,[],1,size(kernel_roi,3)),[],1),size(x,1),1), ...
%     mat2cell(fluor_kernelroiR_deconv,n_trials_day,length(t),size(kernel_roi,3)),'uni',false));
% 
% fluor_kernelroiR_deconv = fluor_kernelroiR_deconv./fluor_kernelroiR_deconv_std;
% 
% fluor_kernelroiR_taskpred_reduced = cell2mat(permute(arrayfun(@(x) ...
%     permute(reshape( ...
%     AP_svd_roi(U_master(:,:,1:n_vs), ...
%     reshape(permute(fluor_taskpred_reduced_allcat(:,:,:,x),[3,2,1]), ...
%     n_vs,[]),[],[],kernel_roiR), ...
%     size(kernel_roi,3),[],size(fluor_taskpred_reduced_allcat,1)),[3,2,1]), ...
%     1:size(fluor_taskpred_reduced_allcat,4),'uni',false),[1,3,4,2]))./fluor_kernelroiR_deconv_std;


figure; 

% Regress Str from Ctx L,R and get cross-validated expl var
use_trials = move_t < 0.5;

k = nan(length(t),3);
ev = nan(length(t),2);
for curr_t = 1:length(t)
    
    use_reduction = 1;
%     curr_data = [mua_allcat(:,curr_t,1) - mua_taskpred_reduced_allcat(:,curr_t,1,use_reduction), ...
%         fluor_kernelroi_deconv(:,curr_t,1) - fluor_kernelroi_taskpred_reduced(:,curr_t,1,use_reduction), ...
%         fluor_kernelroiR_deconv(:,curr_t,1) - fluor_kernelroiR_taskpred_reduced(:,curr_t,1,use_reduction)];
curr_data = [mua_allcat(:,curr_t,1) - mua_taskpred_reduced_allcat(:,curr_t,1,use_reduction), ...
        fluor_roi_deconv(:,curr_t,3) - fluor_roi_taskpred_reduced(:,curr_t,3,use_reduction), ...
        fluor_roi_deconv(:,curr_t,13) - fluor_roi_taskpred_reduced(:,curr_t,1,use_reduction)];
    
    nonan = ~any(isnan(curr_data),2);
    
    [kLR,predicted_signals,evLR] = ...
        AP_regresskernel(curr_data(use_trials & nonan,2:3)',curr_data(use_trials & nonan,1)',0,0,[false,false],10,1,1,[]);    
    k(curr_t,:) = cell2mat(kLR)';
    
    [kL,predicted_signals,evL] = ...
        AP_regresskernel(curr_data(use_trials & nonan,2)',curr_data(use_trials & nonan,1)',0,0,[false,false],10,1,1,[]);
    
    ev(curr_t,1) = evLR.total;
    ev(curr_t,2) = evL.total;
    
end
figure; 
subplot(2,1,1);
plot(t,k,'linewidth',2);
line([0,0],ylim,'color','k');
line(xlim,[0,0],'color','k');
legend({'Ctx_L','Ctx_R','Offset'});
xlabel('Time');
ylabel('Weight');
title('Str regression');

subplot(2,1,2);
plot(t,ev,'linewidth',2);
line([0,0],ylim,'color','k');
ylabel('Explained variance');
xlabel('Time');
legend({'Ctx L & Ctx R','Ctx L'});

% m = nan(length(t),2);
% for curr_t = 1:length(t)
%     
%     %     curr_data = [mua_allcat(:,curr_t,1), ...
%     %         fluor_kernelroi_deconv(:,curr_t,1), ...
%     %         fluor_kernelroiR_deconv(:,curr_t,1)];
%     
%     use_reduction = 1;
%     curr_data = [mua_allcat(:,curr_t,1) - mua_taskpred_reduced_allcat(:,curr_t,1,use_reduction), ...
%         fluor_kernelroi_deconv(:,curr_t,1) - fluor_kernelroi_taskpred_reduced(:,curr_t,1,use_reduction)];
%     
%     nonan = ~any(isnan(curr_data),2);
%     m(curr_t,:) = [curr_data(use_trials & nonan,2:end),ones(sum(use_trials & nonan),1)]\ ...
%         curr_data(use_trials & nonan,1);
%     
% end
% 
% figure;plot(t,m,'linewidth',2);
% legend({'Ctx_L','Offset'});
% xlabel('Time');
% ylabel('Weight');
%     
    




%% Fig 4g: Cortex-predicted striatum error (TESTING ERROR FIX)

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_contrast_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set windows to average activity
timeavg_labels = {'Stim','Move onset','Outcome'};
timeavg_t = {[0.05,0.15],[-0.05,0.05],[0.05,0.15]};
timeavg_align = {stim_align,move_align,outcome_align};
timeavg_trial_conditions = ...
    {[trial_contrastside_allcat > 0,trial_contrastside_allcat < 0], ...
    [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
    [trial_outcome_allcat == 1, trial_outcome_allcat == -1]};
timeavg_task_reduction = [1,2,4];

% timeavg_labels = {'Pre-stim','Stim'};
% timeavg_t = {[-0.2,-0.1],[0.05,0.15]};
% timeavg_align = {stim_align,stim_align};

% Set activity percentiles and bins
act_prctile = [10,90];
n_act_bins = 5;

% Set areas and conditions
plot_areas = [1,2,3,4];

% Loop across area pairs, plot binned predicted v measured activity
curr_act_allcat = mua_allcat;
curr_act_taskpred_reduced_allcat = mua_taskpred_reduced_allcat;

curr_act_pred_allcat = mua_ctxpred_allcat;
curr_act_pred_taskpred_reduced_allcat = mua_ctxpred_taskpred_reduced_allcat;

% curr_act_pred_allcat = mua_taskpred_allcat;
% curr_act_pred_taskpred_reduced_allcat = mua_taskpred_reduced_allcat;

% Get "fixing" matrix: difference between task predicted str/ctx-pred str
task_fix = mua_taskpred_allcat - mua_ctxpred_taskpred_allcat;
task_fix_reduced = mua_taskpred_reduced_allcat - mua_ctxpred_taskpred_reduced_allcat;

measured_v_pred_fig = figure('color','w');
for curr_area_idx = 1:length(plot_areas)
    
    plot_area = plot_areas(curr_area_idx);
    
    % Set up the summary values
    curr_act_pred_diff = nan(2,length(use_split),length(timeavg_labels));
    curr_act_pred_condition_diff = nan(length(use_split),length(timeavg_labels));
    n_shuff = 1000;
    curr_act_pred_condition_diff_shuff = nan(length(use_split),n_shuff,length(timeavg_labels));
    
    for curr_timeavg = 1:length(timeavg_labels)
        
        trial_conditions = timeavg_trial_conditions{curr_timeavg};
        curr_task_reduction = timeavg_task_reduction(curr_timeavg);
                
        curr_act_pred_fix = task_fix - task_fix_reduced(:,:,:,curr_task_reduction);
        
        % (re-align and split activity)
        curr_act = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_allcat(trial,:,plot_area) - curr_act_taskpred_reduced_allcat(trial,:,plot_area,curr_task_reduction), ..., ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        curr_act_pred = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_pred_allcat(trial,:,plot_area) - curr_act_pred_taskpred_reduced_allcat(trial,:,plot_area,curr_task_reduction), ..., ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
%         curr_act_pred = mat2cell(...
%             cell2mat(arrayfun(@(trial) circshift( ...
%             curr_act_pred_allcat(trial,:,plot_area) - curr_act_pred_taskpred_reduced_allcat(trial,:,plot_area,curr_task_reduction) ...
%             + curr_act_pred_fix(trial,:,plot_area), ...
%             timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
%             use_split,length(t));
                
        trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
        
        % (get average activity within window)
        curr_event_t = t > timeavg_t{curr_timeavg}(1) & t <= timeavg_t{curr_timeavg}(2);
        curr_act_avg = cellfun(@(x) double(squeeze(nanmean(x(:,curr_event_t,:),2))),curr_act,'uni',false);
        curr_act_pred_avg = cellfun(@(x) double(squeeze(nanmean(x(:,curr_event_t,:),2))),curr_act_pred,'uni',false);

        % (bin predicted data across percentile range)
        bin_range = prctile(cell2mat(curr_act_pred_avg),act_prctile);
        bin_edges = linspace(bin_range(1),bin_range(2),n_act_bins+1);
        bin_centers = bin_edges(1:end-1) + diff(bin_edges)./2;
        
        trial_bins = cellfun(@(x) discretize(x,bin_edges),curr_act_pred_avg,'uni',false);
        total_bins = cellfun(@(x) discretize(x,[bin_edges(1),bin_edges(end)]),curr_act_pred_avg,'uni',false);
        
        act_bin_range = prctile(cell2mat(curr_act_avg),act_prctile);
        act_bin_edges = linspace(act_bin_range(1),act_bin_range(2),n_act_bins+1);
        act_bin_centers = act_bin_edges(1:end-1) + diff(act_bin_edges)./2;
        act_trial_bins = cellfun(@(x) discretize(x,act_bin_edges),curr_act_avg,'uni',false);
        
        % (get trials to use: no NaNs in either time series or bins)
        use_trials = cellfun(@(act,act_pred,trial_bins) ...
            squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
            squeeze(~any(isnan(act_pred(:,curr_event_t)),2)) & ...
            ~isnan(trial_bins), ...
            curr_act,curr_act_pred,trial_bins,'uni',false);
        
        % (get average binned activity)
        act_binmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,trial_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_binmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,trial_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        
        act_use_trials = cellfun(@(act,act_pred,trial_bins) ...
            squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
            squeeze(~any(isnan(act_pred(:,curr_event_t)),2)) & ...
            ~isnan(trial_bins), ...
            curr_act,curr_act_pred,act_trial_bins,'uni',false);
        
        act_act_binmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,act_trial_bins,trial_conditions_exp,act_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        error_binmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,act_pred,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            abs(act(use_trials & trial_cond(:,condition)) - ...
            act_pred(use_trials & trial_cond(:,condition))), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,curr_act_pred_avg,act_trial_bins,trial_conditions_exp,act_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));

        % (get the average total activity)       
        act_totalmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [1,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,total_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_totalmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [1,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,total_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        % Get average act-pred difference (ALL TRIALS)
        curr_act_pred_diff(:,:,curr_timeavg) = ...
            cell2mat(arrayfun(@(cond) cellfun(@(act,pred,trial_cond,use_trials) ...
            nanmean(abs(act(trial_cond(:,cond)) - pred(trial_cond(:,cond)))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp,use_trials), ...
            1:size(trial_conditions,2),'uni',false))';
        
        % Get condition  difference (assume 2 conditions)
        curr_act_pred_condition_diff(:,curr_timeavg) = ...
            cellfun(@(act,pred,trial_cond,use_trials) ...       
            nanmean(abs(act(trial_cond(:,1)) - pred(trial_cond(:,1)))) - ...
            nanmean(abs(act(trial_cond(:,2)) - pred(trial_cond(:,2)))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp,use_trials);
        
        % Get condition-shuffled act-pred difference
        trial_conditions_exp_shuff = cellfun(@(x) AP_shake(repmat(x,1,1,n_shuff),1), ...
            trial_conditions_exp,'uni',false);   
                
        curr_act_pred_condition_diff_shuff(:,:,curr_timeavg) = ...
            cell2mat(arrayfun(@(shuff) ...
            cellfun(@(act,pred,trial_cond,use_trials) ...
            nanmean(abs(act(trial_cond(:,1,shuff)) - pred(trial_cond(:,1,shuff)))) - ...
            nanmean(abs(act(trial_cond(:,2,shuff)) - pred(trial_cond(:,2,shuff)))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp_shuff,use_trials), ...
            1:n_shuff,'uni',false));
        
        % Plot binned predicted v measured
        figure(measured_v_pred_fig);
        subplot(length(plot_areas),length(timeavg_labels)+2, ...
            sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),curr_timeavg,curr_area_idx));
        hold on;
        
%         errorbar( ...
%             squeeze(nanmean(act_pred_binmean,2)), ...
%             squeeze(nanmean(act_binmean,2)), ...
%             squeeze(AP_sem(act_binmean,2)), ...
%             'linewidth',2,'CapSize',0);
%         errorbar( ...
%             squeeze(nanmean(act_pred_totalmean,2)), ...
%             squeeze(nanmean(act_totalmean,2)), ...
%             squeeze(AP_sem(act_totalmean,2)), ...
%             squeeze(AP_sem(act_totalmean,2)), ...
%             squeeze(AP_sem(act_pred_totalmean,2)), ...
%             squeeze(AP_sem(act_pred_totalmean,2)), ...
%             '.','linewidth',3,'CapSize',0);
%         xlabel(['Predicted (' num2str(plot_area) ')']);
%         ylabel(['Measured (' num2str(plot_area) ')'])
%         ylim(xlim); axis square;
%         title([timeavg_labels{curr_timeavg} ' (' task_regressor_labels{curr_task_reduction} '-reduced)']);

        errorbar( ...
            squeeze(nanmean(error_binmean,2)), ...
            squeeze(AP_sem(act_act_binmean,2)),'linewidth',2);
        xlabel('Measured bin');
        ylabel('Ctx-predicted error');
        
    end
    
    % Plot measured - predicted
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+1,curr_area_idx));
     hold on;
    errorbar( ...
        permute(nanmean(curr_act_pred_diff,2),[3,1,2]), ...
        permute(AP_sem(curr_act_pred_diff,2),[3,1,2]),'linewidth',2,'CapSize',0);
    ylabel('Abs(Meas - Pred)');
    set(gca,'XTick',1:4,'XTickLabels',timeavg_labels,'XTickLabelRotation',45)
    xlim([0.5,length(timeavg_labels)+0.5]);
    
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+2,curr_area_idx));
     hold on;
    errorbar( ...
        nanmean(curr_act_pred_condition_diff,1), ...
        AP_sem(curr_act_pred_condition_diff,1),'linewidth',2,'CapSize',0);
    ylabel('Condition difference');
    set(gca,'XTick',1:4,'XTickLabels',timeavg_labels,'XTickLabelRotation',45)
    xlim([0.5,length(timeavg_labels)+0.5]);
    
    % Get and plot significance
    real_diff = nanmean(curr_act_pred_condition_diff,1); 
    alpha = [5/2/3,100-(5/2/3)];
    shuff_diff_ci = permute(nanmean(prctile(curr_act_pred_condition_diff_shuff,alpha,2),1),[2,3,1]);
    sig_diff = real_diff < shuff_diff_ci(1,:) | real_diff > shuff_diff_ci(2,:);
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+1,curr_area_idx));
    plot(find(sig_diff),repmat(max(ylim),sum(sig_diff),1),'*','color','k','MarkerSize',5);
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+2,curr_area_idx));
    plot(shuff_diff_ci','r');
    plot(find(sig_diff),repmat(max(ylim),sum(sig_diff),1),'*','color','k','MarkerSize',5);
    
end

% Link axes of all predicted v measured and all explained var
all_axes = get(measured_v_pred_fig,'Children');
meas_pred_axes = all_axes(setdiff(1:length(all_axes), ...
    [1:length(timeavg_labels)+2:length(all_axes), ...
    2:length(timeavg_labels)+2:length(all_axes)]));
linkaxes(meas_pred_axes)
linkaxes(all_axes(1:length(timeavg_labels)+2:length(all_axes)));
linkaxes(all_axes(2:length(timeavg_labels)+2:length(all_axes)));

for curr_axes = 1:length(meas_pred_axes)
   line(meas_pred_axes(curr_axes),xlim(meas_pred_axes(curr_axes)), ...
       xlim(meas_pred_axes(curr_axes)),'color','k'); 
end
for curr_axes = 1:length(timeavg_labels)+2:length(all_axes)
   line(all_axes(curr_axes),xlim(all_axes(curr_axes)),[0,0],'color','k');
end
for curr_axes = 2:length(timeavg_labels)+2:length(all_axes)
   line(all_axes(curr_axes),xlim(all_axes(curr_axes)),[0,0],'color','k');
end


%% Fig 4d,e,f: Striatum v Cortex by condition (TESTING ERROR FIX)
% I THINK THIS WAS JUST COPIED, NOT STARTED YET

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_contrast_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set windows to average activity
timeavg_labels = {'Stim','Move onset','Outcome'};
timeavg_t = {[0.05,0.15],[-0.05,0.05],[0.05,0.15]};
timeavg_align = {stim_align,move_align,outcome_align};
timeavg_trial_conditions = ...
    {[trial_contrastside_allcat > 0,trial_contrastside_allcat < 0], ...
    [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
    [trial_outcome_allcat == 1, trial_outcome_allcat == -1]};
timeavg_task_reduction = [1,2,4];

% Set activity percentiles and bins
act_prctile = [10,90];
n_act_bins = 5;

% Set areas and conditions
plot_str_ctx = {[1,1],[2,2],[3,3],[4,4]};

% Get "fixing" matrix: difference between task predicted str/ctx-pred str
task_fix = mua_taskpred_allcat - mua_ctxpred_taskpred_allcat;
task_fix_reduced = mua_taskpred_reduced_allcat - mua_ctxpred_taskpred_reduced_allcat;

% Loop across area pairs, plot binned predicted v measured activity
str_v_ctx_fig = figure('color','w');
for curr_str_ctx = 1:length(plot_str_ctx)
    
    plot_str = plot_str_ctx{curr_str_ctx}(1);
    plot_ctx = plot_str_ctx{curr_str_ctx}(2);

    for curr_mua = 1:2
        
        % (set striatum activity to use)
        switch curr_mua
            case 1
                curr_str_act_allcat = single(mua_allcat);
                curr_str_act_taskpred_reduced_allcat = single(mua_taskpred_reduced_allcat);
                mua_label = 'Measured';
            case 2
                curr_str_act_allcat = single(mua_ctxpred_allcat);
                curr_str_act_taskpred_reduced_allcat = single(mua_ctxpred_taskpred_reduced_allcat);
                
                curr_str_act_taskpred_reduced_allcat = ...
                    curr_str_act_taskpred_reduced_allcat + task_fix_reduced;
                
                
                mua_label = 'Ctx-pred';
        end
        
        for curr_timeavg = 1:length(timeavg_labels)
            
            trial_conditions = timeavg_trial_conditions{curr_timeavg};
            curr_task_reduction = timeavg_task_reduction(curr_timeavg);
            
            % (set cortex activity to use)
%             curr_ctx_act_allcat = fluor_roi_deconv;
            curr_ctx_act_allcat = fluor_kernelroi_deconv;   
            curr_ctx_act_taskpred_reduced_allcat = fluor_kernelroi_taskpred_reduced;
            
            % (re-align and split activity and conditions)
            curr_ctx_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_ctx_act_allcat(trial,:,:) - curr_ctx_act_taskpred_reduced_allcat(trial,:,:,curr_task_reduction), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_ctx_act_allcat,1)),'uni',false)), ...
                use_split,length(t),size(curr_ctx_act_allcat,3));
            
            curr_str_act = mat2cell(...
                cell2mat(arrayfun(@(trial) circshift( ...
                curr_str_act_allcat(trial,:,:) - curr_str_act_taskpred_reduced_allcat(trial,:,:,curr_task_reduction), ...
                timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_str_act_allcat,1)),'uni',false)), ...
                use_split,length(t),size(curr_str_act_allcat,3));
            
            trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
            
            % (get average activity within window)
            curr_event_t = t > timeavg_t{curr_timeavg}(1) & t < timeavg_t{curr_timeavg}(2);
            curr_ctx_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_ctx_act,'uni',false);
            curr_str_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_str_act,'uni',false);
            
            % (bin measured data across percentile range)
            bin_range = prctile(cell2mat(cellfun(@(x) x(:,plot_ctx),curr_ctx_act_avg,'uni',false)),act_prctile);
            bin_edges = linspace(bin_range(1),bin_range(2),n_act_bins+1);
            bin_centers = bin_edges(1:end-1) + diff(bin_edges)./2;
            
            trial_bins = cellfun(@(x) discretize(x,bin_edges),curr_ctx_act_avg,'uni',false);
            total_bins = cellfun(@(x) discretize(x,[bin_edges(1),bin_edges(end)]),curr_ctx_act_avg,'uni',false);
            
            % (get trials to use: no NaNs in either time series or bins)
            nonan_trials = cellfun(@(ctx_act,str_act,trial_bins) ...
                squeeze(~any(isnan(ctx_act(:,curr_event_t,plot_ctx)),2)) & ...
                squeeze(~any(isnan(str_act(:,curr_event_t,plot_str)),2)) & ...
                ~isnan(trial_bins), ...
                curr_ctx_act,curr_str_act,trial_bins,'uni',false);
            
            % (get average binned activity for measured/predicted by condition)
            ctx_act_binmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_ctx_act_avg,trial_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            str_act_binmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [n_act_bins,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_str_act_avg,trial_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            % (get the average total activity)
            ctx_act_totalmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [1,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_ctx_act_avg,total_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            str_act_totalmean = cell2mat(arrayfun(@(condition) ...
                cell2mat(permute(cellfun(@(act,bins,trial_cond,use_trials) cell2mat(arrayfun(@(area) ...
                accumarray(bins(use_trials(:,area) & trial_cond(:,condition),area), ...
                act(use_trials(:,area) & trial_cond(:,condition),area), ...
                [1,1],@nanmean,cast(NaN,class(act))), 1:size(act,2),'uni',false)), ...
                curr_str_act_avg,total_bins,trial_conditions_exp,nonan_trials,'uni',false),[2,3,1])), ...
                permute(1:size(trial_conditions,2),[1,3,4,2]),'uni',false));
            
            % Plot binned predicted v measured           
            subplot(length(plot_str_ctx),length(timeavg_labels), ...
                sub2ind(fliplr([length(plot_str_ctx),length(timeavg_labels)]),curr_timeavg,curr_str_ctx));
            hold on;
            col = lines(size(trial_conditions,2));
            switch curr_mua
                case 1
                    errorbar( ...
                        squeeze(nanmean(ctx_act_binmean(:,plot_ctx,:,:),3)), ...
                        squeeze(nanmean(str_act_binmean(:,plot_str,:,:),3)), ...
                        squeeze(AP_sem(str_act_binmean(:,plot_str,:,:),3)), ...
                        'linewidth',2,'CapSize',0);
%                     errorbar( ...
%                         squeeze(nanmean(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
%                         squeeze(nanmean(str_act_totalmean(:,plot_str,:,:),3)), ...
%                         squeeze(AP_sem(str_act_totalmean(:,plot_str,:,:),3)), ...
%                         squeeze(AP_sem(str_act_totalmean(:,plot_str,:,:),3)), ...
%                         squeeze(AP_sem(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
%                         squeeze(AP_sem(ctx_act_totalmean(:,plot_ctx,:,:),3)), ...
%                         '.','linewidth',3,'CapSize',0);
                case 2
                    for curr_cond = 1:size(trial_conditions,2)
                        AP_errorfill( ...
                            squeeze(nanmean(ctx_act_binmean(:,plot_ctx,:,curr_cond),3)), ...
                            squeeze(nanmean(str_act_binmean(:,plot_str,:,curr_cond),3)), ...
                            squeeze(AP_sem(str_act_binmean(:,plot_str,:,curr_cond),3)), ...
                            col(curr_cond,:),0.5,false);
                    end
            end
            xlabel(['Ctx (' num2str(plot_ctx) ')']);
            ylabel([' Str (' num2str(plot_str) ')'])
            title([timeavg_labels{curr_timeavg} '(' task_regressor_labels{curr_task_reduction} '-reduced)']);
            
        end
        
    end
end

% Link axes of all plots
linkaxes(get(str_v_ctx_fig,'Children'));

%% TESTING NLIN FIT

% Apply empirical static nonlinearity
figure;
mua_ctxpred_allcat_nlin = nan(size(mua_ctxpred_allcat));
% mua_ctxpred_taskpred_allcat_nlin = nan(size(mua_ctxpred_taskpred_allcat));

for curr_depth = 1:n_depths       
    measured_data = reshape(mua_allcat(:,:,curr_depth),[],1);
    predicted_data = double(reshape(mua_ctxpred_allcat(:,:,curr_depth),[],1));
%     predicted_data_task = double(reshape(mua_ctxpred_taskpred_allcat(:,:,curr_depth),[],1));
    
    n_bins = 1000;
    activity_bounds = linspace(-1,9,n_bins+1);
    activity_bin_centers = conv2(activity_bounds,[1,1]/2,'valid');
  
    predicted_bins = discretize(predicted_data,activity_bounds);
%     predicted_task_bins = discretize(predicted_data_task,activity_bounds);
    
    measured_data_binmean = accumarray(predicted_bins(~isnan(predicted_bins)), ...
        measured_data(~isnan(predicted_bins)),[n_bins,1],@mean,nan);
    predicted_data_binmean = accumarray(predicted_bins(~isnan(predicted_bins)),...
        predicted_data(~isnan(predicted_bins)),[n_bins,1],@mean,nan);
    
    % Smooth out the measured data binmean to get nonlinear transform
    measured_data_binmean_smooth = medfilt1(measured_data_binmean,30,'omitnan');
    
    predicted_data_nlin = nan(size(predicted_data));
    predicted_data_nlin(~isnan(predicted_bins)) = measured_data_binmean_smooth(predicted_bins(~isnan(predicted_bins)));
    
%     predicted_data_task_nlin = nan(size(predicted_data_task));
%     predicted_data_task_nlin(~isnan(predicted_task_bins)) = measured_data_binmean_smooth(predicted_task_bins(~isnan(predicted_task_bins)));
    
    predicted_data_nlin_binmean = accumarray( ...
        predicted_bins(~isnan(predicted_bins)), ...
        predicted_data_nlin(~isnan(predicted_bins)),[n_bins,1],@mean,nan);
    
    % Apply nonlinearity
    mua_ctxpred_allcat_nlin(:,:,curr_depth) = ...
        reshape(predicted_data_nlin, ...
        size(mua_ctxpred_allcat_nlin,1),size(mua_ctxpred_allcat_nlin,2));
    
%     mua_ctxpred_taskpred_allcat_nlin(:,:,curr_depth) = ...
%         reshape(predicted_data_task_nlin, ...
%         size(mua_ctxpred_taskpred_allcat_nlin,1),size(mua_ctxpred_taskpred_allcat_nlin,2));
    
    % Plot linear and nonlinear predictions 
    subplot(2,n_depths,curr_depth); hold on;
    AP_heatscatter(predicted_data,measured_data,100)
    plot(predicted_data_binmean,measured_data_binmean,'linewidth',2);
    plot(predicted_data_binmean,measured_data_binmean_smooth,'linewidth',2);
    xlim([-2,9]);ylim(xlim);
    line(xlim,ylim,'color','k');
    xlabel('Predicted')
    ylabel('Measured')
    axis square;
    
    subplot(2,n_depths,curr_depth + n_depths); hold on;
    AP_heatscatter(predicted_data_nlin,measured_data,100)
    plot(predicted_data_nlin_binmean,measured_data_binmean,'linewidth',2);
    xlim([-2,9]);ylim(xlim);
    line(xlim,ylim,'color','k');
    xlabel('Predicted (nonlinear)')
    ylabel('Measured')
    axis square;
    
end

%% TESTING NLIN FIT (by experiment)

mua_allcat_exp = mat2cell(mua_allcat,trials_recording,length(t),n_depths);
mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,trials_recording,length(t),n_depths);
mua_ctxpred_taskpred_allcat_exp = mat2cell(mua_ctxpred_taskpred_allcat,trials_recording,length(t),n_depths);

% Apply empirical static nonlinearity
use_split_cumsum = [0;cumsum(trials_recording)];

mua_ctxpred_allcat_nlin = nan(size(mua_ctxpred_allcat));
mua_ctxpred_taskpred_allcat_nlin = nan(size(mua_ctxpred_taskpred_allcat));

for curr_expt = 1:length(mua_allcat_exp)
    for curr_depth = 1:n_depths
        
        measured_data = reshape(mua_allcat_exp{curr_expt}(:,:,curr_depth),[],1);
        predicted_data = double(reshape(mua_ctxpred_allcat_exp{curr_expt}(:,:,curr_depth),[],1));
        predicted_data_task = double(reshape(mua_ctxpred_taskpred_allcat_exp{curr_expt}(:,:,curr_depth),[],1));
        
        if all(isnan(measured_data(:)))
            continue
        end
        
        n_bins = 1000;
        activity_bounds = linspace(-1,9,n_bins+1);
        activity_bin_centers = conv2(activity_bounds,[1,1]/2,'valid');
        
        predicted_bins = discretize(predicted_data,activity_bounds);
        predicted_task_bins = discretize(predicted_data_task,activity_bounds);
        
        measured_data_binmean = accumarray(predicted_bins(~isnan(predicted_bins)), ...
            measured_data(~isnan(predicted_bins)),[n_bins,1],@mean,nan);
        predicted_data_binmean = accumarray(predicted_bins(~isnan(predicted_bins)),...
            predicted_data(~isnan(predicted_bins)),[n_bins,1],@mean,nan);
        
        % Smooth out the measured data binmean to get nonlinear transform
        measured_data_binmean_smooth = medfilt1(measured_data_binmean,30,'omitnan');
        
        predicted_data_nlin = nan(size(predicted_data));
        predicted_data_nlin(~isnan(predicted_bins)) = measured_data_binmean_smooth(predicted_bins(~isnan(predicted_bins)));
        
        predicted_data_task_nlin = nan(size(predicted_data_task));
        predicted_data_task_nlin(~isnan(predicted_task_bins)) = measured_data_binmean_smooth(predicted_task_bins(~isnan(predicted_task_bins)));
        
        predicted_data_nlin_binmean = accumarray( ...
            predicted_bins(~isnan(predicted_bins)), ...
            predicted_data_nlin(~isnan(predicted_bins)),[n_bins,1],@mean,nan);
        
        % Apply nonlinearity
        mua_ctxpred_allcat_nlin(use_split_cumsum(curr_expt)+1:use_split_cumsum(curr_expt+1),:,curr_depth) = ...
            reshape(predicted_data_nlin,[],length(t));
        
        mua_ctxpred_taskpred_allcat_nlin(use_split_cumsum(curr_expt)+1:use_split_cumsum(curr_expt+1),:,curr_depth) = ...
            reshape(predicted_data_task_nlin,[],length(t));
              
    end
end


% Plot linear and nonlinear predictions
figure
for curr_depth = 1:n_depths
    subplot(2,n_depths,curr_depth); hold on;
    AP_heatscatter(reshape(mua_ctxpred_allcat(:,:,curr_depth),[],1), ...
        reshape(mua_allcat(:,:,curr_depth),[],1),200);
    xlim([-2,9]);ylim(xlim);
    line(xlim,ylim,'color','k');
    xlabel('Predicted')
    ylabel('Measured')
    axis square;
    
    subplot(2,n_depths,curr_depth + n_depths); hold on;
    AP_heatscatter(reshape(mua_ctxpred_allcat_nlin(:,:,curr_depth),[],1), ...
        reshape(mua_allcat(:,:,curr_depth),[],1),200);xlim([-2,9]);ylim(xlim);
    line(xlim,ylim,'color','k');
    xlabel('Predicted (nonlinear)')
    ylabel('Measured')
    axis square;
end

%% TESTING NLIN FIT (by experiment) (fit on "spontaneous"?)

mua_allcat_exp = mat2cell(mua_allcat - mua_taskpred_allcat,trials_recording,length(t),n_depths);
mua_ctxpred_allcat_minustask_exp = mat2cell(mua_ctxpred_allcat - mua_ctxpred_taskpred_allcat,trials_recording,length(t),n_depths);

mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,trials_recording,length(t),n_depths);
mua_ctxpred_taskpred_allcat_exp = mat2cell(mua_ctxpred_taskpred_allcat,trials_recording,length(t),n_depths);



% Apply empirical static nonlinearity
use_split_cumsum = [0;cumsum(trials_recording)];

mua_ctxpred_allcat_nlin = nan(size(mua_ctxpred_allcat));
mua_ctxpred_taskpred_allcat_nlin = nan(size(mua_ctxpred_taskpred_allcat));

for curr_expt = 1:length(mua_allcat_exp)
    for curr_depth = 1:n_depths
        
        measured_data = reshape(mua_allcat_exp{curr_expt}(:,:,curr_depth),[],1);
        predicted_data_minustask = double(reshape(mua_ctxpred_allcat_minustask_exp{curr_expt}(:,:,curr_depth),[],1));
        predicted_data = double(reshape(mua_ctxpred_allcat_exp{curr_expt}(:,:,curr_depth),[],1));
        predicted_data_task = double(reshape(mua_ctxpred_taskpred_allcat_exp{curr_expt}(:,:,curr_depth),[],1));
        
        if all(isnan(measured_data(:)))
            continue
        end
        
        n_bins = 1000;
        activity_bounds = linspace(-1,9,n_bins+1);
        activity_bin_centers = conv2(activity_bounds,[1,1]/2,'valid');
        
        predicted_bins = discretize(predicted_data_minustask,activity_bounds);
        predicted_task_bins = discretize(predicted_data_task,activity_bounds);
        
        measured_data_binmean = accumarray(predicted_bins(~isnan(predicted_bins)), ...
            measured_data(~isnan(predicted_bins)),[n_bins,1],@mean,nan);
        predicted_data_binmean = accumarray(predicted_bins(~isnan(predicted_bins)),...
            predicted_data_minustask(~isnan(predicted_bins)),[n_bins,1],@mean,nan);
        
        % Smooth out the measured data binmean to get nonlinear transform
        measured_data_binmean_smooth = medfilt1(measured_data_binmean,30,'omitnan');
        
        predicted_withtask_bins = discretize(predicted_data,activity_bounds);
        
        predicted_data_nlin = nan(size(predicted_data));
        predicted_data_nlin(~isnan(predicted_withtask_bins)) = measured_data_binmean_smooth(predicted_withtask_bins(~isnan(predicted_withtask_bins)));
        
        predicted_data_task_nlin = nan(size(predicted_data_task));
        predicted_data_task_nlin(~isnan(predicted_task_bins)) = measured_data_binmean_smooth(predicted_task_bins(~isnan(predicted_task_bins)));
                
        % Apply nonlinearity
        mua_ctxpred_allcat_nlin(use_split_cumsum(curr_expt)+1:use_split_cumsum(curr_expt+1),:,curr_depth) = ...
            reshape(predicted_data_nlin,[],length(t));
        
        mua_ctxpred_taskpred_allcat_nlin(use_split_cumsum(curr_expt)+1:use_split_cumsum(curr_expt+1),:,curr_depth) = ...
            reshape(predicted_data_task_nlin,[],length(t));
              
    end
end


% Plot linear and nonlinear predictions
figure
for curr_depth = 1:n_depths
    subplot(2,n_depths,curr_depth); hold on;
    AP_heatscatter(reshape(mua_ctxpred_allcat(:,:,curr_depth),[],1), ...
        reshape(mua_allcat(:,:,curr_depth),[],1),200);
    xlim([-2,9]);ylim(xlim);
    line(xlim,ylim,'color','k');
    xlabel('Predicted')
    ylabel('Measured')
    axis square;
    
    subplot(2,n_depths,curr_depth + n_depths); hold on;
    AP_heatscatter(reshape(mua_ctxpred_allcat_nlin(:,:,curr_depth),[],1), ...
        reshape(mua_allcat(:,:,curr_depth),[],1),200);xlim([-2,9]);ylim(xlim);
    line(xlim,ylim,'color','k');
    xlabel('Predicted (nonlinear)')
    ylabel('Measured')
    axis square;
end



%% Make toy data for sanity check


use_trials = move_t < 0.5;
toy_stim_response = grpstats(mua_allcat(use_trials,:,1),trial_contrastside_allcat(use_trials));
[~,cond_idx] = ismember(trial_contrastside_allcat,unique(trial_contrastside_allcat),'rows');

toy_stim_trial = toy_stim_response(cond_idx,:);
toy_noise_trial = rand(size(mua_allcat,1),size(mua_allcat,2));

toy_str = toy_stim_trial + toy_noise_trial + rand(size(mua_allcat,1),size(mua_allcat,2));
% toy_str_ctxpred = toy_stim_trial*0.2 + toy_noise_trial + rand(size(mua_allcat,1),size(mua_allcat,2));
toy_str_ctxpred = max(toy_stim_trial - 0.5,0) + toy_noise_trial + rand(size(mua_allcat,1),size(mua_allcat,2));

toy_stim_response_str = grpstats(toy_str(use_trials,:,1),trial_contrastside_allcat(use_trials));
toy_stim_response_str_ctxpred = grpstats(toy_str_ctxpred(use_trials,:,1),trial_contrastside_allcat(use_trials));

toy_fix_additive = toy_stim_response_str(cond_idx,:) - toy_stim_response_str_ctxpred(cond_idx,:);
toy_fix_multiplicative = toy_stim_response_str(cond_idx,:)./toy_stim_response_str_ctxpred(cond_idx,:);


%% (checking error from nlin)

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_contrast_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set windows to average activity
% timeavg_labels = {'Stim','Move onset','Outcome'};
% timeavg_t = {[0.05,0.15],[-0.05,0.05],[0.05,0.15]};
% timeavg_align = {stim_align,move_align,outcome_align};
% timeavg_trial_conditions = ...
%     {[sign(trial_contrastside_allcat) == 1,sign(trial_contrastside_allcat) == -1], ...
%     [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
%     [trial_outcome_allcat == 1, trial_outcome_allcat == -1]};
% timeavg_task_reduction = [1,2,4];

timeavg_labels = {'Pre-stim','Stim'};
timeavg_t = {[-0.2,-0.1],[0.05,0.15]};
timeavg_align = {stim_align,stim_align};
timeavg_trial_conditions = ...
    {[trial_contrastside_allcat > 0, trial_contrastside_allcat < 0], ...
    [trial_contrastside_allcat > 0, trial_contrastside_allcat < 0]};

% timeavg_labels = {'Pre-stim','Stim'};
% timeavg_t = {[-0.2,-0.1],[0.05,0.15]};
% timeavg_align = {stim_align,stim_align};

% Set activity percentiles and bins
act_prctile = [10,90];
n_act_bins = 5;

% Set areas and conditions
plot_areas = [1];

% Loop across area pairs, plot binned predicted v measured activity
curr_act_allcat = mua_allcat;
% curr_act_pred_allcat = mua_ctxpred_allcat;
curr_act_pred_allcat = mua_ctxpred_allcat_nlin;

% curr_act_allcat = toy_str;
% curr_act_pred_allcat = toy_str_ctxpred;

% Get "fixing" matrix: difference between task predicted str/ctx-pred str
% task_fix = (mua_taskpred_allcat - mua_ctxpred_taskpred_allcat);
task_fix = (mua_taskpred_allcat - mua_ctxpred_taskpred_allcat_nlin);

% task_fix = toy_fix_additive;

%%% TESTING: estimate a crappy version here?
% use_trials = move_t < 0.5;
% est_str_stim = grpstats(curr_act_allcat(use_trials,:,1),trial_contrastside_allcat(use_trials));
% est_str_ctx_stim = grpstats(curr_act_pred_allcat(use_trials,:,1),trial_contrastside_allcat(use_trials));
% est_str_ctx_stim_fix = max(est_str_stim - est_str_ctx_stim,0);
% 
% [~,cond_idx] = ismember(trial_contrastside_allcat,unique(trial_contrastside_allcat),'rows');
% task_fix = est_str_ctx_stim_fix(cond_idx,:);

measured_v_pred_fig = figure('color','w');
for curr_area_idx = 1:length(plot_areas)
    
    plot_area = plot_areas(curr_area_idx);
    
    % Set up the summary values
    curr_act_pred_diff = nan(2,length(use_split),length(timeavg_labels));
    curr_act_pred_condition_diff = nan(length(use_split),length(timeavg_labels));
    n_shuff = 1000;
    curr_act_pred_condition_diff_shuff = nan(length(use_split),n_shuff,length(timeavg_labels));
    
    for curr_timeavg = 1:length(timeavg_labels)
        
        trial_conditions = timeavg_trial_conditions{curr_timeavg};
        
        % (re-align and split activity)
        curr_act = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
%         curr_act_pred = mat2cell(...
%             cell2mat(arrayfun(@(trial) circshift( ...
%             curr_act_pred_allcat(trial,:,plot_area), ...
%             timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
%             use_split,length(t));
        
        curr_act_pred = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_pred_allcat(trial,:,plot_area) + task_fix(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
            use_split,length(t));
                
        trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
        
        % (get average activity within window)
        curr_event_t = t > timeavg_t{curr_timeavg}(1) & t < timeavg_t{curr_timeavg}(2);
        curr_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act,'uni',false);
        curr_act_pred_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_pred,'uni',false);

        % (bin predicted data across percentile range)
        bin_range = prctile(cell2mat(curr_act_pred_avg),act_prctile);
        bin_edges = linspace(bin_range(1),bin_range(2),n_act_bins+1);
        bin_centers = bin_edges(1:end-1) + diff(bin_edges)./2;
        
        trial_bins = cellfun(@(x) discretize(x,bin_edges),curr_act_pred_avg,'uni',false);
        total_bins = cellfun(@(x) discretize(x,[bin_edges(1),bin_edges(end)]),curr_act_pred_avg,'uni',false);
        
        % (get trials to use: no NaNs in either time series or bins)
        use_trials = cellfun(@(act,act_taskpred,trial_bins) ...
            squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
            squeeze(~any(isnan(act_taskpred(:,curr_event_t)),2)) & ...
            ~isnan(trial_bins), ...
            curr_act,curr_act_pred,trial_bins,'uni',false);
        
        % (get average binned activity)
        act_binmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,trial_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_binmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,trial_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));      

        % (get the average total activity)       
        act_totalmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [1,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,total_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_totalmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [1,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,total_bins,trial_conditions_exp,use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        % Get average act-pred difference
        curr_act_pred_diff(:,:,curr_timeavg) = ...
            cell2mat(arrayfun(@(cond) cellfun(@(act,pred,trial_cond,use_trials) ...
            nanmean(abs(act(trial_cond(:,cond) & use_trials) - pred(trial_cond(:,cond) & use_trials))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp,use_trials), ...
            1:size(trial_conditions,2),'uni',false))';
        
        % Get condition  difference (assume 2 conditions)
        curr_act_pred_condition_diff(:,curr_timeavg) = ...
            cellfun(@(act,pred,trial_cond,use_trials) ...       
            nanmean(abs(act(trial_cond(:,1) & use_trials) - pred(trial_cond(:,1) & use_trials))) - ...
            nanmean(abs(act(trial_cond(:,2) & use_trials) - pred(trial_cond(:,2) & use_trials))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp,use_trials);
        
        % Get condition-shuffled act-pred difference
        trial_conditions_exp_shuff = cellfun(@(x) AP_shake(repmat(x,1,1,n_shuff),1), ...
            trial_conditions_exp,'uni',false);   
                
        curr_act_pred_condition_diff_shuff(:,:,curr_timeavg) = ...
            cell2mat(arrayfun(@(shuff) ...
            cellfun(@(act,pred,trial_cond,use_trials) ...
            nanmean(abs(act(trial_cond(:,1,shuff)) - pred(trial_cond(:,1,shuff)))) - ...
            nanmean(abs(act(trial_cond(:,2,shuff)) - pred(trial_cond(:,2,shuff)))), ...
            curr_act_avg,curr_act_pred_avg,trial_conditions_exp_shuff,use_trials), ...
            1:n_shuff,'uni',false));
        
        % Plot binned predicted v measured
        figure(measured_v_pred_fig);
        subplot(length(plot_areas),length(timeavg_labels)+2, ...
            sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),curr_timeavg,curr_area_idx));
        hold on;
        
        errorbar( ...
            squeeze(nanmean(act_pred_binmean,2)), ...
            squeeze(nanmean(act_binmean,2)), ...
            squeeze(AP_sem(act_binmean,2)), ...
            'linewidth',2,'CapSize',0);
        errorbar( ...
            squeeze(nanmean(act_pred_totalmean,2)), ...
            squeeze(nanmean(act_totalmean,2)), ...
            squeeze(AP_sem(act_totalmean,2)), ...
            squeeze(AP_sem(act_totalmean,2)), ...
            squeeze(AP_sem(act_pred_totalmean,2)), ...
            squeeze(AP_sem(act_pred_totalmean,2)), ...
            '.','linewidth',3,'CapSize',0);
        xlabel(['Predicted (' num2str(plot_area) ')']);
        ylabel(['Measured (' num2str(plot_area) ')'])
        ylim(xlim); axis square;
        title([timeavg_labels{curr_timeavg} ' (' task_regressor_labels{curr_task_reduction} '-reduced)']);
        
    end
    
    % Plot measured - predicted
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+1,curr_area_idx));
     hold on;
    errorbar( ...
        permute(nanmean(curr_act_pred_diff,2),[3,1,2]), ...
        permute(AP_sem(curr_act_pred_diff,2),[3,1,2]),'linewidth',2,'CapSize',0);
    ylabel('Meas - Pred');
    set(gca,'XTick',1:4,'XTickLabels',timeavg_labels,'XTickLabelRotation',45)
    xlim([0.5,length(timeavg_labels)+0.5]);
    
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+2,curr_area_idx));
     hold on;
    errorbar( ...
        nanmean(curr_act_pred_condition_diff,1), ...
        AP_sem(curr_act_pred_condition_diff,1),'linewidth',2,'CapSize',0);
    ylabel('Condition difference');
    set(gca,'XTick',1:4,'XTickLabels',timeavg_labels,'XTickLabelRotation',45)
    xlim([0.5,length(timeavg_labels)+0.5]);
    
    % Get and plot significance
    real_diff = nanmean(curr_act_pred_condition_diff,1); 
    alpha = [5/2/3,100-(5/2/3)];
    shuff_diff_ci = permute(nanmean(prctile(curr_act_pred_condition_diff_shuff,alpha,2),1),[2,3,1]);
    sig_diff = real_diff < shuff_diff_ci(1,:) | real_diff > shuff_diff_ci(2,:);
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+1,curr_area_idx));
    plot(find(sig_diff),repmat(max(ylim),sum(sig_diff),1),'*','color','k','MarkerSize',5);
    subplot(length(plot_areas),length(timeavg_labels)+2, ...
        sub2ind(fliplr([length(plot_areas),length(timeavg_labels)+2]),length(timeavg_labels)+2,curr_area_idx));
    plot(shuff_diff_ci','r');
    plot(find(sig_diff),repmat(max(ylim),sum(sig_diff),1),'*','color','k','MarkerSize',5);
    
end

% Link axes of all predicted v measured and all explained var
all_axes = get(measured_v_pred_fig,'Children');
meas_pred_axes = all_axes(setdiff(1:length(all_axes), ...
    [1:length(timeavg_labels)+2:length(all_axes), ...
    2:length(timeavg_labels)+2:length(all_axes)]));
linkaxes(meas_pred_axes)
linkaxes(all_axes(1:length(timeavg_labels)+2:length(all_axes)));
linkaxes(all_axes(2:length(timeavg_labels)+2:length(all_axes)));

for curr_axes = 1:length(meas_pred_axes)
   line(meas_pred_axes(curr_axes),xlim(meas_pred_axes(curr_axes)), ...
       xlim(meas_pred_axes(curr_axes)),'color','k'); 
end
for curr_axes = 1:length(timeavg_labels)+2:length(all_axes)
   line(all_axes(curr_axes),xlim(all_axes(curr_axes)),[0,0],'color','k');
end
for curr_axes = 2:length(timeavg_labels)+2:length(all_axes)
   line(all_axes(curr_axes),xlim(all_axes(curr_axes)),[0,0],'color','k');
end




%% (checking if error is just amplitude-dependent)

% a = reshape(mua_allcat(:,:,1)',[],1);
% b = reshape(mua_ctxpred_allcat(:,:,1)',[],1);
% c = reshape(mua_ctxpred_allcat_nlin(:,:,1)',[],1);
% d = reshape(mua_taskpred_allcat(:,:,1)',[],1);

% see if error is better explained by contrast or by value?

% current error
m = convn(mua_allcat,ones(1,5)./5,'same');
r = m(:,15:30,1) - mua_ctxpred_allcat(:,15:30,1);

% additive offset by contrast
z = grpstats(r(:,:,1),trial_stim_allcat);
[~,cond_idx] = ismember(trial_stim_allcat,unique(trial_stim_allcat),'rows');
z_tr = z(cond_idx,:);
rz = r - z_tr;

% scaling?
m2 = m(:,15:30,1);
mc2 = mua_ctxpred_allcat(:,15:30,1);
nonan = ~isnan(m2) & ~isnan(mc2);
error_scale = [mc2(nonan),ones(sum(nonan(:)),1)]\reshape(m2(nonan),[],1);
rs = m(:,15:30,1) - (mua_ctxpred_allcat(:,15:30,1)*error_scale(1) + error_scale(2));

% rs = m(:,15:30,1) - mua_ctxpred_allcat(:,15:30,1)*1.7;

figure;
subplot(1,3,1);
AP_heatscatter(reshape(m(:,15:30,1),[],1),reshape(r,[],1),200)
line(xlim,xlim);
xlabel('measured');ylabel('error');

subplot(1,3,2);
AP_heatscatter(reshape(m(:,15:30,1),[],1),reshape(rz,[],1),200)
line(xlim,xlim);
xlabel('measured');ylabel('error');

subplot(1,3,3);
AP_heatscatter(reshape(m(:,15:30,1),[],1),reshape(rs,[],1),200)
line(xlim,xlim);
xlabel('measured');ylabel('error');

linkaxes(get(gcf,'Children'),'xy');

%% Testing Fig 4 bin plot examples

curr_act_cat = cell2mat(curr_act);
curr_act_pred_cat = cell2mat(curr_act_pred);

unique_stim = unique(trial_contrastside_allcat);

% Plot trials on top of each other, separated by contrast
plot_prctiles = linspace(10,90,10);

figure; 
col = colormap_BlueWhiteRed(5);
col(6,:) = [0.5,0.5,0.5];
for curr_stim = 1:length(unique_stim)
    subplot(1,length(unique_stim),curr_stim); hold on;
        
    curr_trials = find(~any(isnan(curr_act_pred_cat(:,curr_event_t)),2) & ...
        trial_contrastside_allcat == unique_stim(curr_stim));
    trial_act = nanmean(curr_act_pred_cat(curr_trials,curr_event_t),2);
    [~,sort_idx] = sort(trial_act);
    
    plot_trials = curr_trials(sort_idx(round(prctile(1:length(curr_trials),plot_prctiles))));
    
    plot(t,curr_act_pred_cat(plot_trials,:)','color',col(curr_stim,:));
    plot(t,nanmean(curr_act_pred_cat(curr_trials,:),1),'color',max(col(curr_stim,:)-0.2,0),'linewidth',2);
    title(['Stim: ' num2str(unique_stim(curr_stim))]);
end
linkaxes(get(gcf,'Children'),'xy');

% Plot average +/- percentile for each trial condition
figure;
for curr_cond = 1:size(trial_conditions,2)
    subplot(1,size(trial_conditions,2),curr_cond);
    curr_trials = trial_conditions(:,curr_cond) & ...
        ~any(isnan(curr_act_pred_cat(:,curr_event_t)),2);
    
    AP_errorfill(t,nanmean(curr_act_pred_cat(curr_trials,:),1), ...
        prctile(curr_act_pred_cat(curr_trials,:),[10,90],1) - ...
        nanmean(curr_act_pred_cat(curr_trials,:),1), ...
        'k',0.5,true);
    for curr_edge = 1:length(pred_bin_edges)
        line(xlim,repmat(pred_bin_edges(curr_edge),2,1),'color','k');
    end
end
linkaxes(get(gcf,'Children'),'xy');

% Plot distribution of activity by contrast
figure; hold on
distributionPlot(cell2mat(curr_act_pred_avg), ...
    'groups',trial_contrastside_allcat, ...
    'color',mat2cell(col,ones(11,1),3), ...
    'xvalues',1:11,'showMM',0,'globalNorm',1);
for curr_edge = 1:length(pred_bin_edges)
   line(xlim,repmat(pred_bin_edges(curr_edge),2,1),'color','k');
end
xlabel('Contrast*Side');
ylabel('Cortex');

% Split cortex into high/low, plot cortex and striatum
figure; 
for curr_stim = 1:length(unique_stim)
    subplot(1,length(unique_stim),curr_stim); hold on;
        
    curr_trials = find(~any(isnan(curr_act_pred_cat(:,curr_event_t)),2) & ...
        trial_contrastside_allcat == unique_stim(curr_stim));
    
    trial_act = nanmean(curr_act_pred_cat(curr_trials,curr_event_t),2);
    [~,sort_idx] = sort(trial_act);
    
    high_ctx_trials = curr_trials(sort_idx(round(prctile(1:length(curr_trials),50))+1:length(curr_trials)));
    low_ctx_trials = curr_trials(sort_idx(1:round(prctile(1:length(curr_trials),50))));
    
    subplot(2,length(unique_stim),curr_stim); hold on;
    plot(t,nanmean(curr_act_pred_cat(high_ctx_trials,:),1),'linewidth',2,'color',[0.3,0.8,0.3]);
    plot(t,nanmean(curr_act_pred_cat(low_ctx_trials,:),1),'linewidth',2,'color',[0,0.3,0]);
    title(['Stim: ' num2str(unique_stim(curr_stim))]);
    
    subplot(2,length(unique_stim),curr_stim + length(unique_stim)); hold on;
    plot(t,nanmean(curr_act_cat(high_ctx_trials,:),1),'linewidth',2,'color',[0.6,0.6,0.6]);
    plot(t,nanmean(curr_act_cat(low_ctx_trials,:),1),'linewidth',2,'color',[0.3,0.3,0.3]);
    
end
linkaxes(get(gcf,'Children'),'xy');


%% (make new alignment for all animals)

animals = {'AP024','AP025','AP026','AP027','AP028','AP029', ...
    'AP032','AP033','AP034','AP035','AP036'};

for curr_animal = 8:length(animals)
    
    animal = animals{curr_animal};
    % Find experiments
    % (use only behavior days because cortical recordings afterwards)
    protocol = 'vanillaChoiceworld';
    experiments = AP_find_experiments(animal,protocol);
    experiments(~([experiments.imaging])) = [];
    if isempty(experiments)
        % (if no behavior days then it was a naive mouse - use passive expt)
        protocol = 'AP_choiceWorldStimPassive';
        experiments = AP_find_experiments(animal,protocol);
        experiments(~([experiments.imaging] & [experiments.ephys])) = [];
    end
    
    avg_im_days = cell(size(experiments));
    for curr_day = 1:length(experiments)
        day = experiments(curr_day).day;
        [img_path,img_exists] = AP_cortexlab_filename(animal,day,[],'imaging');
        avg_im = readNPY([img_path filesep 'meanImage_purple.npy']);
        avg_im_days{curr_day} = avg_im;
    end
    
    AP_align_widefield(avg_im_days,animal,{experiments.day},'new_days');
        
end


%% (make new master retinotopy)

retinotopy_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment\retinotopy';
retinotopy_dir = dir(retinotopy_path);

animal_retinotopy_idx = cellfun(@(x) ~isempty(x), regexp({retinotopy_dir.name},'AP\d*_retinotopy'));
animals_tokens = cellfun(@(x) regexp({x},'(AP\d*)_retinotopy','tokens'),{retinotopy_dir.name});
animals = cellfun(@(x) cell2mat(x{:}),animals_tokens(cellfun(@(x) ~isempty(x),animals_tokens)),'uni',false);

vfs_animal = cell(length(animals),1);
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    load([retinotopy_path filesep animal '_retinotopy'])
    
    aligned_vfs = cell(length(retinotopy),1);
    for curr_day = 1:length(retinotopy)
        try
        aligned_vfs{curr_day} = AP_align_widefield(retinotopy(curr_day).vfs,animal,retinotopy(curr_day).day,'day_only');
        catch me
           disp(['Skipping ' animal ' ' retinotopy(curr_day).day]); 
        end
    end
    
    vfs_animal{curr_animal} = nanmean(cat(3,aligned_vfs{:}),3);
    
end

AP_align_widefield(vfs_animal,[],[],'create_master');

%% (re-align animals to master)

retinotopy_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment\retinotopy';
retinotopy_dir = dir(retinotopy_path);

animal_retinotopy_idx = cellfun(@(x) ~isempty(x), regexp({retinotopy_dir.name},'AP\d*_retinotopy'));
animals_tokens = cellfun(@(x) regexp({x},'(AP\d*)_retinotopy','tokens'),{retinotopy_dir.name});
animals = cellfun(@(x) cell2mat(x{:}),animals_tokens(cellfun(@(x) ~isempty(x),animals_tokens)),'uni',false);

vfs_animal = cell(length(animals),1);
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    load([retinotopy_path filesep animal '_retinotopy'])
    
    aligned_vfs = cell(length(retinotopy),1);
    for curr_day = 1:length(retinotopy)
        try
        aligned_vfs{curr_day} = AP_align_widefield(retinotopy(curr_day).vfs,animal,retinotopy(curr_day).day,'day_only');
        catch me
           disp(['Skipping ' animal ' ' retinotopy(curr_day).day]); 
        end
    end
    
    vfs_mean = nanmean(cat(3,aligned_vfs{:}),3);
    AP_align_widefield(vfs_mean,animal,[],'new_animal');

end

%% (check master)

retinotopy_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment\retinotopy';
retinotopy_dir = dir(retinotopy_path);

animal_retinotopy_idx = cellfun(@(x) ~isempty(x), regexp({retinotopy_dir.name},'AP\d*_retinotopy'));
animals_tokens = cellfun(@(x) regexp({x},'(AP\d*)_retinotopy','tokens'),{retinotopy_dir.name});
animals = cellfun(@(x) cell2mat(x{:}),animals_tokens(cellfun(@(x) ~isempty(x),animals_tokens)),'uni',false);

vfs_animal = cell(length(animals),1);
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    load([retinotopy_path filesep animal '_retinotopy'])
    
    aligned_vfs = cell(length(retinotopy),1);
    for curr_day = 1:length(retinotopy)
        try
        aligned_vfs{curr_day} = AP_align_widefield(retinotopy(curr_day).vfs,animal,retinotopy(curr_day).day);
%             aligned_vfs{curr_day} = AP_align_widefield_old(animal,retinotopy(curr_day).day,retinotopy(curr_day).vfs);
        catch me
           disp(['Skipping ' animal ' ' retinotopy(curr_day).day]); 
        end
    end
    
    vfs_animal{curr_animal} = nanmean(cat(3,aligned_vfs{:}),3);
    
end


%% Average full screen flicker

stim_px = nan(size(U_master,1),size(U_master,2),length(t),3);
for curr_stim = 1:3
   curr_v = permute(nanmean(fluor_allcat(trial_stim_allcat == curr_stim,:,:),1),[3,2,1]);
   curr_px = svdFrameReconstruct(U_master(:,:,1:n_vs),curr_v);
   stim_px(:,:,:,curr_stim) = curr_px;
end

AP_imscroll(stim_px,t);
axis image;
AP_reference_outline('ccf_aligned','k');
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));


%% Widefield: get px std for normalizing

fluor_allcat_deconv_exp = mat2cell(fluor_allcat_deconv,use_split,length(t),n_vs);

% this takes ages
fluor_day_baseline_px = cellfun(@(fluor) nanstd(svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    reshape(fluor(:,t_baseline,:),[],n_vs)'),[],3), ...
    fluor_allcat_deconv_exp,'uni',false);

% Put the std into V-space
U_master_reshape = reshape(U_master(:,:,1:n_vs),[],n_vs);
fluor_day_baseline = cellfun(@(fluor) ...
    U_master_reshape\reshape(fluor,[],1),fluor_day_baseline_px,'uni',false);

% Get the normalized fluor
% NOPE APPARENTLY CANT DO THIS IN V SPACE??
fluor_allcat_deconv_exp_norm = cellfun(@(fluor,norm) fluor./reshape(norm,1,1,[]), ...
    fluor_allcat_deconv_exp,fluor_day_baseline,'uni',false);


%% Passive stim average (ctx+str)

% %%%%%%%%% TESTING: AVERAGE BY EXPERIMENT to see if it's cleaner
% 
% stim_px = nan(size(U_master,1),size(U_master,2),length(t),length(unique_stim),length(use_split));
% for curr_exp = 1:length(use_split)
%     for curr_stim = 1:length(unique_stim)
%         curr_v = permute(nanmean(fluor_allcat_deconv(split_idx == curr_exp & ...
%             trial_stim_allcat == unique_stim(curr_stim),:,:),1),[3,2,1]);
%         curr_px = svdFrameReconstruct(U_master(:,:,1:n_vs),curr_v);
%         stim_px(:,:,:,curr_stim,curr_exp) = curr_px;
%     end
% end
% 
% trial_stim_allcat_exp = mat2cell(trial_stim_allcat,trials_recording,1);
% fluor_allcat_deconv_exp = mat2cell(fluor_allcat_deconv,use_split,length(t),n_vs);
% quiescent_trials_exp = mat2cell(quiescent_trials,trials_recording,1);
% 
% curr_v = cellfun(@(act,stim,quiesc) ...
%     permute(nanmean(act(stim == 1 & quiesc,:,:),1),[3,2,1]), ...
%     fluor_allcat_deconv_exp,trial_stim_allcat_exp,quiescent_trials_exp,'uni',false);
% 
% curr_v_mean = nanmean(cat(3,curr_v{:}),3);
% 
% curr_px = svdFrameReconstruct(U_master(:,:,1:n_vs),curr_v_mean);
% 
% AP_imscroll(curr_px,t)
% axis image;
% AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
% caxis([-max(abs(caxis)),max(abs(caxis))]);
% colormap(brewermap([],'*RdBu'));

%%%%%%%%%

% Get trials with movement during stim to exclude
quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > 0.025,2);
% quiescent_trials = move_t < 0.5;

unique_stim = unique(trial_stim_allcat);

% Plot widefield movie
stim_px = nan(size(U_master,1),size(U_master,2),length(t),length(unique_stim));
for curr_stim = 1:length(unique_stim)
   curr_v = permute(nanmean(fluor_allcat_deconv(quiescent_trials & ...
       trial_stim_allcat == unique_stim(curr_stim),:,:),1),[3,2,1]);
   curr_px = svdFrameReconstruct(U_master(:,:,1:n_vs),curr_v);
   stim_px(:,:,:,curr_stim) = curr_px;
end

AP_imscroll(stim_px,t);
axis image;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));

% Plot cortex roi / striatum over time
fluor_roi_stim_avg = nan(length(length(unique_stim)),length(t),n_rois,length(use_split));
fluor_kernelroi_stim_avg = nan(length(length(unique_stim)),length(t),n_depths,length(use_split));
mua_stim_avg = nan(length(length(unique_stim)),length(t),n_depths,length(use_split));
mua_ctxpred_stim_avg = nan(length(length(unique_stim)),length(t),n_depths,length(use_split));
for curr_exp = 1:length(use_split)
    for curr_stim = 1:length(unique_stim)
            
        fluor_roi_stim_avg(curr_stim,:,:,curr_exp) = ...
            nanmean(fluor_roi_deconv(split_idx == curr_exp & quiescent_trials ...
            & trial_stim_allcat == unique_stim(curr_stim),:,:),1);
        
        fluor_kernelroi_stim_avg(curr_stim,:,:,curr_exp) = ...
            nanmean(fluor_kernelroi_deconv(split_idx == curr_exp & quiescent_trials ...
            & trial_stim_allcat == unique_stim(curr_stim),:,:),1);
        
         mua_stim_avg(curr_stim,:,:,curr_exp) = ...
            nanmean(mua_allcat(split_idx == curr_exp & quiescent_trials ...
            & trial_stim_allcat == unique_stim(curr_stim),:,:),1);
        
         mua_ctxpred_stim_avg(curr_stim,:,:,curr_exp) = ...
            nanmean(mua_ctxpred_allcat(split_idx == curr_exp & quiescent_trials ...
            & trial_stim_allcat == unique_stim(curr_stim),:,:),1);
        
    end
end

stim_all = unique([0,0.06,0.125,0.25,0.5,1].*[-1;1]);
col_all = brewermap(11,'*RdBu');
col_all(6,:) = [0.5,0.5,0.5];
[~,stim_idx] = ismember(unique_stim,stim_all,'rows');
col = col_all(stim_idx,:);

figure;
for curr_depth = 1:n_depths
    
    subplot(n_depths,2,(curr_depth-1)*2+1); hold on;
    AP_errorfill(t,nanmean(fluor_kernelroi_stim_avg(:,:,curr_depth,:),4)', ...
        AP_sem(fluor_roi_stim_avg(:,:,curr_depth,:),4)',col);
    title('Fluor (kernel ROI)');
    
    subplot(n_depths,2,(curr_depth-1)*2+2); hold on;
    AP_errorfill(t,nanmean(mua_stim_avg(:,:,curr_depth,:),4)', ...
        AP_sem(mua_stim_avg(:,:,curr_depth,:),4)',col);
    title('MUA');
    
%     subplot(n_depths,3,(curr_depth-1)*3+3); hold on;
%     AP_errorfill(t,nanmean(mua_ctxpred_stim_avg(:,:,curr_depth,:),4)', ...
%         AP_sem(mua_ctxpred_stim_avg(:,:,curr_depth,:),4)',col);
%     title('MUA ctxpred');
    
end
linkaxes(get(gcf,'Children'),'xy');


% Plot cortex vs fluorescence in timepoint
plot_str = 1;

use_t = t > 0 & t < 0.2;

figure; 
subplot(1,2,1); hold on;
plot_x = squeeze(nanmean(mua_stim_avg(:,use_t,plot_str,:),2));
plot_y = squeeze(nanmean(fluor_kernelroi_stim_avg(:,use_t,plot_str,:),2));
errorbar(nanmean(plot_x,2),nanmean(plot_y,2), ...
    AP_sem(plot_y,2),AP_sem(plot_y,2), ...
    AP_sem(plot_x,2),AP_sem(plot_x,2),'linestyle','none','color','k');
scatter(nanmean(plot_x,2),nanmean(plot_y,2), ...
    30,col,'filled','MarkerEdgeColor','k');
xlabel(['Striatum (' num2str(plot_str) ')']);
ylabel(['Cortex ROI (' num2str(plot_str) ')']);
line(xlim,xlim,'color','k')

subplot(1,2,2); hold on;
plot_x = squeeze(nanmean(mua_stim_avg(:,use_t,plot_str,:),2));
plot_y = squeeze(nanmean(mua_ctxpred_stim_avg(:,use_t,plot_str,:),2));
errorbar(nanmean(plot_x,2),nanmean(plot_y,2), ...
    AP_sem(plot_y,2),AP_sem(plot_y,2), ...
    AP_sem(plot_x,2),AP_sem(plot_x,2),'linestyle','none','color','k');
scatter(nanmean(plot_x,2),nanmean(plot_y,2), ...
    30,col,'filled','MarkerEdgeColor','k');
xlabel(['Striatum (' num2str(plot_str) ')']);
ylabel(['Cortex kernel (' num2str(plot_str) ')']);
line(xlim,xlim,'color','k')




%% (temp data club) Plot measured/ctx-predicted spikes

% Set options
surround_window = [-0.5,3];
baseline_window = [-0.1,0];

surround_samplerate = 1/(framerate*1);
surround_time = surround_window(1):surround_samplerate:surround_window(2);
baseline_surround_time = baseline_window(1):surround_samplerate:baseline_window(2);

% (passive)
use_stims = find(stimIDs == 3);
% (choiceworld)
% stimIDs = trial_conditions(:,1).*trial_conditions(:,2);
% use_stims = find(stimIDs > 0);

use_stimOn_times = stimOn_times(use_stims);

stim_surround_times = bsxfun(@plus, use_stimOn_times(:), surround_time);
stim_baseline_surround_times = bsxfun(@plus, use_stimOn_times(:), baseline_surround_time);

psth_measured_stim = permute(interp1(time_bin_centers,binned_spikes_std',stim_surround_times),[3,2,1]);
psth_predicted_stim = permute(interp1(time_bin_centers,predicted_spikes',stim_surround_times),[3,2,1]);

psth_measured_baseline = permute(nanmean(interp1(time_bin_centers,binned_spikes_std',stim_baseline_surround_times),2),[3,2,1]);
psth_predicted_baseline = permute(nanmean(interp1(time_bin_centers,predicted_spikes',stim_baseline_surround_times),2),[3,2,1]);

psth_measured = nanmean(psth_measured_stim - psth_measured_baseline,3);
psth_predicted = nanmean(psth_predicted_stim - psth_predicted_baseline,3);

figure; hold on;
AP_stackplot(psth_measured',surround_time,2,[],'k');
AP_stackplot(psth_predicted',surround_time,2,[],'r');
% plot(surround_time,psth_measured,'k','linewidth',2);
% plot(surround_time,psth_predicted,'r','linewidth',2);

%% Plot lilrig water calibration

load //zserver.cortexlab.net/Code/Rigging/config/LILRIG-STIM/hardware.mat

figure; hold on;
n_cal = length(daqController.SignalGenerators(1, 1).Calibrations);
set(gca,'ColorOrder',copper(n_cal));
for plot_cal = 1:n_cal-1
    plot([daqController.SignalGenerators(1, 1).Calibrations(plot_cal).measuredDeliveries.durationSecs],...
        [daqController.SignalGenerators(1, 1).Calibrations(plot_cal).measuredDeliveries.volumeMicroLitres])
end

plot([daqController.SignalGenerators(1, 1).Calibrations(n_cal).measuredDeliveries.durationSecs],...
    [daqController.SignalGenerators(1, 1).Calibrations(n_cal).measuredDeliveries.volumeMicroLitres], ...
    'linewidth',3,'color','r')
title(datestr(daqController.SignalGenerators(1, 1).Calibrations(n_cal).dateTime))
xlabel('duration (sec)'); ylabel('volume (uL)')



%% Karolinska temp

use_t = find(t > -0.2 & t < 0.8);
use_im = ap_memmap(:,:,use_t,:);

framerate = 35/5;
color_map = colormap;
color_axis = caxis;
figure_position = get(gcf,'Position');
fn = 'avg_center_stim';
savefile = ['C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\presentations\\200302_data_club\' fn '.avi'];
t_annotation = arrayfun(@(x) sprintf('Time from stim: %.2f',t(x)),use_t,'uni',false);
AP_movie2avi(use_im,framerate,color_map,color_axis,figure_position,savefile,t_annotation)





plot_t = time_bin_centers > 455 & time_bin_centers < 4;
figure; hold on
plot(time_bin_centers(plot_t),binned_spikes_std(1,plot_t),'k')
plot(time_bin_centers(plot_t),ctxpred_spikes_std(1,plot_t),'g')
x = xlim;

line([x(1),x(1)],[0,2],'linewidth',2);
line([x(1),x(1)+1],[0,0],'linewidth',2);




% Plot str1/AM
use_trials = move_t < 0.5;
plot_t = t > -0.2 & t < 1;
ap_memmap = grpstats(mua_allcat(use_trials,:,1),trial_contrastside_allcat(use_trials));
b = grpstats(fluor_roi_deconv(use_trials,:,3),trial_contrastside_allcat(use_trials));

stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = 0.5;

figure;
subplot(1,2,1); hold on;
set(gca,'ColorOrder',stim_col);
plot(t(plot_t),ap_memmap(:,plot_t)','linewidth',2);
line([0,0],ylim,'color','k');
line(repmat(nanmean(move_t(use_trials)),2,1),ylim,'color','k');
line([0.5,0.5],ylim,'color','k');

line([0,0],[0,1],'linewidth',3,'color','k');
line([0,0.2],[0,0],'linewidth',3,'color','k');

subplot(1,2,2); hold on;
set(gca,'ColorOrder',stim_col);
plot(t(plot_t),b(:,plot_t)','linewidth',2);
line([0,0],ylim,'color','k');
line(repmat(nanmean(move_t(use_trials)),2,1),ylim,'color','k');
line([0.5,0.5],ylim,'color','k');

line([0,0],[0,0.004],'linewidth',3,'color','k');




% (get sum spikes in passive condition)

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(mua_allcat,1),1);

% Get trials with movement during stim to exclude
wheel_thresh = 0.025;
quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > wheel_thresh,2);

% Set windows to average activity
timeavg_labels = {'Stim'};
timeavg_t = {[0.05,0.15]};
timeavg_align = {stim_align};
timeavg_trial_conditions = ...
    {[trial_stim_allcat > 0 & quiescent_trials, ...
    trial_stim_allcat < 0 & quiescent_trials]};
timeavg_condition_colors = ...
    {[1,0,0;0,0,1]};

% Loop across area pairs, plot binned predicted v measured activity
curr_act_allcat = mua_allcat;

% (ctx-predicted)
curr_act_pred_allcat = mua_ctxpred_allcat;

% (fix by average stim response within experiment)
mua_allcat_exp = mat2cell(mua_allcat,use_split,length(t),n_depths);
mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,use_split,length(t),n_depths);
fluor_roi_deconv_exp = mat2cell(fluor_roi_deconv,use_split,length(t),numel(wf_roi));

trial_contrastside_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
trial_conditions_exp = mat2cell(timeavg_trial_conditions{1},use_split, ...
    size(timeavg_trial_conditions{1},2));

str_stim_exp = nan(length(unique(trial_stim_allcat)),length(t),n_depths,length(use_split));
str_ctxpred_stim_exp = nan(length(unique(trial_stim_allcat)),length(t),n_depths,length(use_split));
ctx_stim_exp = nan(length(unique(trial_stim_allcat)),length(t),n_rois,length(use_split));

curr_act_pred_fix_allcat_exp = mua_ctxpred_allcat_exp;
for curr_exp = 1:length(use_split)
    
    curr_stim = unique(trial_contrastside_allcat_exp{curr_exp});
    for curr_stim_idx = 1:length(curr_stim)
        curr_trials = any(trial_conditions_exp{curr_exp},2) & ...
            trial_contrastside_allcat_exp{curr_exp} == curr_stim(curr_stim_idx);
        curr_act_stim_avg = nanmean(mua_allcat_exp{curr_exp}(curr_trials,:,:),1);
        curr_act_pred_stim_avg = nanmean(mua_ctxpred_allcat_exp{curr_exp}(curr_trials,:,:),1);
        curr_act_ctx_stim_avg = nanmean(fluor_roi_deconv_exp{curr_exp}(curr_trials,:,:),1);
        
        str_stim_exp(curr_stim_idx,:,:,curr_exp) = curr_act_stim_avg;
        str_ctxpred_stim_exp(curr_stim_idx,:,:,curr_exp) = curr_act_pred_stim_avg;
        ctx_stim_exp(curr_stim_idx,:,:,curr_exp) = curr_act_ctx_stim_avg;   
    end
end

use_stim = unique(trial_stim_allcat);

use_t = t >= 0.05 & t <= 0.15;
ap_memmap = squeeze(nanmean(sum(str_stim_exp(:,use_t,1,:),2),4));
b = squeeze(nanmean(sum(str_ctxpred_stim_exp(:,use_t,1,:),2),4));
figure; hold on;
plot(use_stim,ap_memmap,'k','linewidth',2);
plot(use_stim,b,'g','linewidth',2);


ap_memmap = squeeze(sum(str_stim_exp(use_stim == 1,use_t,1,:),2));
b = squeeze(sum(str_ctxpred_stim_exp(use_stim == 1,use_t,1,:),2));
c = squeeze(sum(ctx_stim_exp(use_stim == 1,use_t,3,:),2));
d = [ap_memmap,b,c];
figure; errorbar(nanmean(d,1),AP_sem(d,1),'k','linewidth',2);
xlim([0.5,2.5]);




% Plot stim response
stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = 0.5;
stim_col_contrastside = unique([0,0.06,0.125,0.25,0.5,1].*[-1;1]);
[~,used_stim_idx] = ismember(unique(trial_contrastside_allcat), ...
    stim_col_contrastside,'rows');

str_stim_avg = nanmean(str_stim_exp,4);
figure('color','w');
for curr_area_idx = 1:length(plot_areas)
    curr_area = plot_areas(curr_area_idx);
    subplot(length(plot_areas),2,length(plot_areas)*(curr_area_idx-1)+1); 
    hold on;
    for curr_stim_idx = 1:size(str_stim_avg,1)
        AP_errorfill(t, ...
            nanmean(str_stim_exp(curr_stim_idx,:,curr_area,:),4), ...
            AP_sem(str_stim_exp(curr_stim_idx,:,curr_area,:),4), ...
            stim_col(used_stim_idx(curr_stim_idx),:),0.5,true);
    end
    line([0,0],ylim,'color','k');
    ylabel(['Measured (' num2str(curr_area) ')']);
    xlabel('Time from stim');
    
    subplot(length(plot_areas),2,length(plot_areas)*(curr_area_idx-1)+2); 
    hold on;
    for curr_stim_idx = 1:size(str_stim_avg,1)
        AP_errorfill(t, ...
            nanmean(ctx_stim_exp(curr_stim_idx,:,curr_area,:),4), ...
            AP_sem(ctx_stim_exp(curr_stim_idx,:,curr_area,:),4), ...
            stim_col(used_stim_idx(curr_stim_idx),:),0.5,true);
    end
    line([0,0],ylim,'color','k');
    ylabel(['Predicted (' num2str(curr_area) ')']);
    xlabel('Time from stim');
end


%% Testing fig 4 passive: ROI instead of prediction
% DON'T COPY THIS CODE FOR ANYTHING, PLAYED WITH LOTS OF THINGS

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(mua_allcat,1),1);

% Get trials with movement during stim to exclude
wheel_thresh = 0.025;
quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > wheel_thresh,2);

% Set windows to average activity
timeavg_labels = {'Stim'};
timeavg_t = {[-0.3,-0.2]};
% timeavg_t = {[0.05,0.15]};
timeavg_align = {stim_align};
timeavg_trial_conditions = ...
    {[trial_stim_allcat > 0 & quiescent_trials, ...
    trial_stim_allcat < 0 & quiescent_trials]};
timeavg_condition_colors = ...
    {[1,0,0;0,0,1]};

% Set activity percentiles and bins
act_prctile = [5,95];
n_act_bins = 5;

% Set areas and conditions
plot_areas = [1];

% Loop across area pairs, plot binned predicted v measured activity
curr_act_allcat = mua_allcat;

% (ctx-predicted)
curr_act_pred_allcat = mua_ctxpred_allcat;

% (fix by average stim response within experiment)
mua_allcat_exp = mat2cell(mua_allcat,use_split,length(t),n_depths);
mua_ctxpred_allcat_exp = mat2cell(mua_ctxpred_allcat,use_split,length(t),n_depths);
trial_stim_allcat_exp = mat2cell(trial_stim_allcat,use_split,1);
trial_conditions_exp = mat2cell(timeavg_trial_conditions{1},use_split, ...
    size(timeavg_trial_conditions{1},2));

str_stim_exp = nan(length(unique(trial_stim_allcat)),length(t),n_depths,length(use_split));
ctx_stim_exp = nan(length(unique(trial_stim_allcat)),length(t),n_depths,length(use_split));
curr_act_pred_fix_allcat_exp = mua_ctxpred_allcat_exp;
for curr_exp = 1:length(use_split)
    
    curr_stim = unique(trial_stim_allcat_exp{curr_exp});
    for curr_stim_idx = 1:length(curr_stim)
        curr_trials = any(trial_conditions_exp{curr_exp},2) & ...
            trial_stim_allcat_exp{curr_exp} == curr_stim(curr_stim_idx);
        curr_act_stim_avg = nanmean(mua_allcat_exp{curr_exp}(curr_trials,:,:),1);
        curr_act_pred_stim_avg = nanmean(mua_ctxpred_allcat_exp{curr_exp}(curr_trials,:,:),1);
        
        str_stim_exp(curr_stim_idx,:,:,curr_exp) = curr_act_stim_avg;
        ctx_stim_exp(curr_stim_idx,:,:,curr_exp) = curr_act_pred_stim_avg;
        
        curr_stim_fix = curr_act_stim_avg - curr_act_pred_stim_avg;        
        curr_act_pred_fix_allcat_exp{curr_exp}(curr_trials,:,:) = ...
            curr_act_pred_fix_allcat_exp{curr_exp}(curr_trials,:,:) + curr_stim_fix;       
    end
end
curr_act_pred_fix_allcat = cell2mat(curr_act_pred_fix_allcat_exp);

% Plot stim response
stim_col = colormap_BlueWhiteRed(5);
stim_col(6,:) = 0.5;
stim_col_contrastside = unique([0,0.06,0.125,0.25,0.5,1].*[-1;1]);
[~,used_stim_idx] = ismember(unique(trial_stim_allcat), ...
    stim_col_contrastside,'rows');

str_stim_avg = nanmean(str_stim_exp,4);
figure('color','w');
for curr_area_idx = 1:length(plot_areas)
    curr_area = plot_areas(curr_area_idx);
    subplot(length(plot_areas),2,length(plot_areas)*(curr_area_idx-1)+1); 
    hold on;
    for curr_stim_idx = 1:size(str_stim_avg,1)
        AP_errorfill(t, ...
            nanmean(str_stim_exp(curr_stim_idx,:,curr_area,:),4), ...
            AP_sem(str_stim_exp(curr_stim_idx,:,curr_area,:),4), ...
            stim_col(used_stim_idx(curr_stim_idx),:),0.5,true);
    end
    line([0,0],ylim,'color','k');
    ylabel(['Measured (' num2str(curr_area) ')']);
    xlabel('Time from stim');
    
    subplot(length(plot_areas),2,length(plot_areas)*(curr_area_idx-1)+2); 
    hold on;
    for curr_stim_idx = 1:size(str_stim_avg,1)
        AP_errorfill(t, ...
            nanmean(ctx_stim_exp(curr_stim_idx,:,curr_area,:),4), ...
            AP_sem(ctx_stim_exp(curr_stim_idx,:,curr_area,:),4), ...
            stim_col(used_stim_idx(curr_stim_idx),:),0.5,true);
    end
    line([0,0],ylim,'color','k');
    ylabel(['Predicted (' num2str(curr_area) ')']);
    xlabel('Time from stim');
end
linkaxes(get(gcf,'Children'),'xy');

measured_v_pred_fig = figure('color','w');
for curr_area_idx = 1:length(plot_areas)   
    
    plot_area = plot_areas(curr_area_idx);  
    
    for curr_timeavg = 1:length(timeavg_labels)
        
        trial_conditions = timeavg_trial_conditions{curr_timeavg};
        
        % (re-align and split activity)
        act_title = timeavg_labels{curr_timeavg};
        
        curr_act = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        curr_act_pred = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_pred_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        curr_act_predfix = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_pred_fix_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_fix_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
        
        % (get average activity within window)
        curr_event_t = t >= timeavg_t{curr_timeavg}(1) & t <= timeavg_t{curr_timeavg}(2);
        curr_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act,'uni',false);
        curr_act_pred_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_pred,'uni',false);
        curr_act_predfix_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_predfix,'uni',false);
        
        % (bin predicted data across percentile range)
        pred_bin_edges = prctile(cell2mat(curr_act_pred_avg),linspace(act_prctile(1),act_prctile(2),n_act_bins+1));
        pred_bin_centers = pred_bin_edges(1:end-1) + diff(pred_bin_edges)./2;        
        pred_trial_bins = cellfun(@(x) discretize(x,pred_bin_edges),curr_act_pred_avg,'uni',false);         
        
%         %TESTING:  (bin MEASURED data across percentile range)
%         pred_bin_edges = prctile(cell2mat(curr_act_avg),linspace(act_prctile(1),act_prctile(2),n_act_bins+1));
%         pred_bin_centers = pred_bin_edges(1:end-1) + diff(pred_bin_edges)./2;        
%         pred_trial_bins = cellfun(@(x) discretize(x,pred_bin_edges),curr_act_avg,'uni',false);         
        
        % (get activity binned by predicted)        
        pred_use_trials = cellfun(@(act,act_pred,trial_bins) ...
            squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
            squeeze(~any(isnan(act_pred(:,curr_event_t)),2)) & ...
            ~isnan(trial_bins), ...
            curr_act,curr_act_pred,pred_trial_bins,'uni',false);
        
        act_predbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_predbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        % (get "fixed" predicted activity binned by predicted)
        act_predfix_predbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_predfix_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        % Get condition difference significance from shuffle
        n_shuff = 1000;
        trial_conditions_shuff = cellfun(@(x) AP_shake(x,1), ...
            repmat(trial_conditions_exp,1,n_shuff),'uni',false);
        
        act_predbinmean_condshuff = cell2mat(arrayfun(@(curr_shuff) ...
            cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,pred_trial_bins,trial_conditions_shuff(:,curr_shuff), ...
            pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false)), ...
            permute(1:n_shuff,[1,3,4,2]),'uni',false));
        
        act_predfix_predbinmean_condshuff = cell2mat(arrayfun(@(curr_shuff) ...
            cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_predfix_avg,pred_trial_bins,trial_conditions_shuff(:,curr_shuff), ...
            pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false)), ...
            permute(1:n_shuff,[1,3,4,2]),'uni',false));
        
        % (measured data: null = no difference between conditions)
        cond_combos = nchoosek(1:size(trial_conditions,2),2);
        cond_sig_diff = false(n_act_bins,size(cond_combos,1));
        predfix_cond_sig_diff = false(n_act_bins,size(cond_combos,1));
        for curr_cond_combo = 1:size(cond_combos,1)
            curr_combo_diff = nanmean(act_predbinmean(:,:,cond_combos(curr_cond_combo,1)) - ...
                 act_predbinmean(:,:,cond_combos(curr_cond_combo,2)),2);
            shuff_prctile = squeeze(prctile(nanmean(act_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,1),:) - ...
                act_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,2),:),2),[2.5,97.5],4));
            cond_sig_diff(:,curr_cond_combo) = ...
                curr_combo_diff < shuff_prctile(:,1) | ...
                curr_combo_diff > shuff_prctile(:,2);
            
            predfix_curr_combo_diff = nanmean(act_predfix_predbinmean(:,:,cond_combos(curr_cond_combo,1)) - ...
                 act_predfix_predbinmean(:,:,cond_combos(curr_cond_combo,2)),2);
            predfix_shuff_prctile = squeeze(prctile(nanmean(act_predfix_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,1),:) - ...
                act_predfix_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,2),:),2),[2.5,97.5],4));
            predfix_cond_sig_diff(:,curr_cond_combo) = ...
                predfix_curr_combo_diff < predfix_shuff_prctile(:,1) | ...
                predfix_curr_combo_diff > predfix_shuff_prctile(:,2);
        end 
        
        % Plot measured v predicted in bins
        figure(measured_v_pred_fig)
        
        % (measured vs binned predicted)
        subplot(length(plot_areas),length(timeavg_labels), ...
            sub2ind(fliplr([length(plot_areas),length(timeavg_labels)]),curr_timeavg,curr_area_idx));
        hold on;
        set(gca,'ColorOrder',timeavg_condition_colors{curr_timeavg});

        fill_cols = min(timeavg_condition_colors{curr_timeavg} + 0.5,1);
%         for curr_cond = 1:size(trial_conditions,2)
%             AP_errorfill( ...
%                 squeeze(nanmean(act_pred_predbinmean(:,:,curr_cond),2)), ...
%                 squeeze(nanmean(act_predfix_predbinmean(:,:,curr_cond),2)), ...
%                 squeeze(AP_sem(act_predfix_predbinmean(:,:,curr_cond),2)), ...
%                 fill_cols(curr_cond,:),1,false);
%         end
        
        errorbar( ...
            squeeze(nanmean(act_pred_predbinmean,2)), ...
            squeeze(nanmean(act_predbinmean,2)), ...
            squeeze(AP_sem(act_predbinmean,2)), ...
            'linestyle','-','linewidth',2,'CapSize',10);
        
        n_col_bins = n_act_bins + 2;
        ctx_col = brewermap(n_col_bins,'*Greens');        
        scatter( ...
            reshape(squeeze(nanmean(act_pred_predbinmean,2)),[],1), ...
            reshape(squeeze(nanmean(act_predbinmean,2)),[],1), ...
            60,repmat(ctx_col(1:n_act_bins,:),size(trial_conditions,2),1),'filled');
        
        xlabel(['Predicted (' num2str(plot_area) ')']);
        ylabel(['Measured (' num2str(plot_area) ')'])
        axis square;
        title(act_title);
        
        % (plot significant measured condition differences)
        % (* and o = significant in both measured and "fixed" predicted)
        curr_ylim = max(ylim);
        for curr_cond_combo = 1:size(cond_combos,1)
            % (plot * for measured condition differences)
            sig_y = curr_ylim + 0.1*curr_cond_combo;
            sig_x = pred_bin_centers;
            if any(cond_sig_diff(:,curr_cond_combo))
                plot(sig_x(cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                    '*','MarkerSize',10,'color', ...
                    timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,1),:));
                plot(sig_x(cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                    '*','MarkerSize',5,'color', ...
                    timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,2),:));
            end
            % (plot o for [predicted condition differences)
            if any(predfix_cond_sig_diff(:,curr_cond_combo))
                plot(sig_x(predfix_cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                    'o','MarkerSize',15,'color', ...
                    timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,1),:));
                plot(sig_x(predfix_cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                    'o','MarkerSize',10,'color', ...
                    timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,2),:));
            end
        end
        
        drawnow;
        
    end
end

linkaxes(get(measured_v_pred_fig,'Children'),'xy');



% Plot long trace measured + predicted with stim marked
plot_trials = trial_stim_allcat == 1;

figure; hold on;
plot(reshape(mua_allcat(plot_trials,:,plot_areas)',[],1));
plot(reshape(mua_ctxpred_allcat(plot_trials,:,plot_areas)',[],1));

stim_matrix = zeros(size(mua_allcat(plot_trials,:,1)));
stim_matrix(:,find(t >= 0,1)) = 1;
stim_idx = find(stim_matrix');
for i = 1:length(stim_idx)
   line(repmat(stim_idx(i),2,1),ylim,'color','k');
end

% Plot distributions
curr_act_avg_cat = cell2mat(curr_act_avg);
curr_act_pred_avg_cat = cell2mat(curr_act_pred_avg);

figure;
subplot(1,2,1);
AP_heatscatter(reshape(mua_ctxpred_allcat(:,:,plot_areas),[],1),reshape(mua_allcat(:,:,plot_areas),[],1),50)
line(xlim,xlim)
xlabel('Fluorescence');
ylabel('MUA');
title('All points')

subplot(1,2,2); hold on;
use_bins = linspace(-2,10,100);
act_dist_1 = histcounts2(curr_act_pred_avg_cat(trial_stim_allcat == 1), ...
    curr_act_avg_cat(trial_stim_allcat == 1),use_bins,use_bins);
act_dist_2 = histcounts2(curr_act_pred_avg_cat(trial_stim_allcat == -1), ...
    curr_act_avg_cat(trial_stim_allcat == -1),use_bins,use_bins);
imagesc(imgaussfilt(act_dist_1',2)-imgaussfilt(act_dist_2',2));
line(xlim,xlim)
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(gca,brewermap([],'*RdBu'));
title('Stim act distr difference');

figure;
plot_trials = trial_stim_allcat ~= 0;
scatterhist(curr_act_pred_avg_cat(plot_trials), ...
    curr_act_avg_cat(plot_trials), ...
    'Group',trial_stim_allcat(plot_trials), ...
    'Kernel','on','Color','br','Marker','.');




figure; hold on;
use_bins = linspace(-2,10,100);
act_dist_1 = histcounts2( ...
    reshape(mua_ctxpred_allcat(trial_stim_allcat == 1,:,plot_areas)',[],1), ...
    reshape(mua_allcat(trial_stim_allcat == 1,:,plot_areas)',[],1), ...
    use_bins,use_bins);
act_dist_2 = histcounts2( ...
    reshape(mua_ctxpred_allcat(trial_stim_allcat == -1,:,plot_areas)',[],1), ...
    reshape(mua_allcat(trial_stim_allcat == -1,:,plot_areas)',[],1), ...
    use_bins,use_bins);
imagesc(imgaussfilt(act_dist_1',2)-imgaussfilt(act_dist_2',2));
line(xlim,xlim)
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(gca,brewermap([],'*RdBu'));
title('Stim act distr difference (whole trials)');


%% Testing fig 4: flip axes

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(mua_allcat,1),1);
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set windows to average activity
timeavg_labels = {'Stim','Move early','Move late','Outcome'};
timeavg_t = {[0.05,0.15],[-0.05,0.05],[0.05,0.15],[0.05,0.15]};
timeavg_align = {stim_align,move_align,move_align,outcome_align};
timeavg_trial_conditions = ...
    {[trial_stim_allcat > 0,trial_stim_allcat == 0,trial_stim_allcat < 0], ...
    [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
    [trial_choice_allcat == -1,trial_choice_allcat == 1], ...
    [trial_outcome_allcat == 1, trial_outcome_allcat == -1]};
timeavg_condition_colors = ...
    {[1,0,0;0,0,0;0,0,1], ...
    [0.6,0,0.6;0,0.6,0], ...
    [0.6,0,0.6;0,0.6,0], ...
    [0,0,0.7;0,0,0]};
timeavg_event_offset = {'Stim','Move onset','Move onset','Outcome'};

% Set activity percentiles and bins
act_prctile = [5,95];
n_act_bins = 5;

% Set areas and conditions
plot_areas = [1:n_depths];

% Loop across area pairs, plot binned predicted v measured activity
curr_act_allcat = mua_ctxpred_allcat;
curr_act_pred_allcat = mua_allcat;

% (old: to use all events)
% task_fix = mua_taskpred_allcat - mua_ctxpred_taskpred_allcat;
% curr_act_pred_fix_allcat = curr_act_pred_allcat + task_fix;

% (to use individual events)
task_fix = (mua_taskpred_allcat - mua_taskpred_reduced_allcat) - ...
    (mua_ctxpred_taskpred_allcat - mua_ctxpred_taskpred_reduced_allcat);
curr_act_pred_fix_allcat = curr_act_pred_allcat + task_fix;

measured_v_pred_fig = figure('color','w');
for curr_area_idx = 1:length(plot_areas)   
    
    plot_area = plot_areas(curr_area_idx);  
    
    for curr_timeavg = 1:length(timeavg_labels)
        
        trial_conditions = timeavg_trial_conditions{curr_timeavg};
        curr_event_offset_idx = find(strcmp(timeavg_event_offset{curr_timeavg},task_regressor_labels));
        curr_act_pred_fix_allcat_event = ...
            curr_act_pred_fix_allcat(:,:,:,curr_event_offset_idx);
        
        % (re-align and split activity)
        act_title = timeavg_labels{curr_timeavg};
        
        curr_act = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        curr_act_pred = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_pred_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        curr_act_predfix = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_pred_fix_allcat_event(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_fix_allcat_event,1)),'uni',false)), ...
            use_split,length(t));
        
        trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
        
        % (get average activity within window)
        curr_event_t = t >= timeavg_t{curr_timeavg}(1) & t <= timeavg_t{curr_timeavg}(2);
        curr_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act,'uni',false);
        curr_act_pred_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_pred,'uni',false);
        curr_act_predfix_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_predfix,'uni',false);
        
        % (bin predicted data across percentile range)
        pred_bin_edges = prctile(cell2mat(curr_act_pred_avg),linspace(act_prctile(1),act_prctile(2),n_act_bins+1));
        pred_bin_centers = pred_bin_edges(1:end-1) + diff(pred_bin_edges)./2;        
        pred_trial_bins = cellfun(@(x) discretize(x,pred_bin_edges),curr_act_pred_avg,'uni',false);         
        
        % (get activity binned by predicted)        
        pred_use_trials = cellfun(@(act,act_pred,trial_bins) ...
            squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
            squeeze(~any(isnan(act_pred(:,curr_event_t)),2)) & ...
            ~isnan(trial_bins), ...
            curr_act,curr_act_pred,pred_trial_bins,'uni',false);
        
        act_predbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_predbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        % (get "fixed" predicted activity binned by predicted)
        act_predfix_predbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_predfix_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        % Get condition difference significance from shuffle
        n_shuff = 1;
        trial_conditions_shuff = cellfun(@(x) AP_shake(x,1), ...
            repmat(trial_conditions_exp,1,n_shuff),'uni',false);
        
        act_predbinmean_condshuff = cell2mat(arrayfun(@(curr_shuff) ...
            cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,pred_trial_bins,trial_conditions_shuff(:,curr_shuff), ...
            pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false)), ...
            permute(1:n_shuff,[1,3,4,2]),'uni',false));
        
        act_predfix_predbinmean_condshuff = cell2mat(arrayfun(@(curr_shuff) ...
            cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_predfix_avg,pred_trial_bins,trial_conditions_shuff(:,curr_shuff), ...
            pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false)), ...
            permute(1:n_shuff,[1,3,4,2]),'uni',false));
        
        % (measured data: null = no difference between conditions)
        cond_combos = nchoosek(1:size(trial_conditions,2),2);
        cond_sig_diff = false(n_act_bins,size(cond_combos,1));
        predfix_cond_sig_diff = false(n_act_bins,size(cond_combos,1));
        for curr_cond_combo = 1:size(cond_combos,1)
            curr_combo_diff = nanmean(act_predbinmean(:,:,cond_combos(curr_cond_combo,1)) - ...
                 act_predbinmean(:,:,cond_combos(curr_cond_combo,2)),2);
            shuff_prctile = squeeze(prctile(nanmean(act_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,1),:) - ...
                act_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,2),:),2),[2.5,97.5],4));
            cond_sig_diff(:,curr_cond_combo) = ...
                curr_combo_diff < shuff_prctile(:,1) | ...
                curr_combo_diff > shuff_prctile(:,2);
            
            predfix_curr_combo_diff = nanmean(act_predfix_predbinmean(:,:,cond_combos(curr_cond_combo,1)) - ...
                 act_predfix_predbinmean(:,:,cond_combos(curr_cond_combo,2)),2);
            predfix_shuff_prctile = squeeze(prctile(nanmean(act_predfix_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,1),:) - ...
                act_predfix_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,2),:),2),[2.5,97.5],4));
            predfix_cond_sig_diff(:,curr_cond_combo) = ...
                predfix_curr_combo_diff < predfix_shuff_prctile(:,1) | ...
                predfix_curr_combo_diff > predfix_shuff_prctile(:,2);
        end 
      
        % Plot binned measured, predicted, and error (by predicted bins)
        measured_pred_fig = figure('color','w','Name', ...
            ['Str ' num2str(plot_area)' ', ' timeavg_labels{curr_timeavg}]);
        n_col_bins = n_act_bins + 2;
        
        [binned_act_pred_t,binned_act_pred_grp] = grpstats(cell2mat(curr_act_pred), ...
            [cell2mat(pred_trial_bins),trial_conditions],{'nanmean','gname'});
        binned_act_pred_grp = cellfun(@str2num,binned_act_pred_grp);        
        
        [binned_act_t,binned_act_grp] = grpstats(cell2mat(curr_act), ...
            [cell2mat(pred_trial_bins),trial_conditions],{'nanmean','gname'});
        binned_act_grp = cellfun(@str2num,binned_act_grp);
        
        % (plot predicted data timecourse)
        for curr_cond = 1:size(trial_conditions,2)           
            subplot(2,size(trial_conditions,2), ...
                sub2ind(fliplr([2,size(trial_conditions,2)]),curr_cond,1)); hold on;
                     
            curr_mean = nanmean(cell2mat(cellfun(@(act,use_trials,cond) ...
                act(use_trials & cond(:,curr_cond),:), ...
                curr_act_pred,pred_use_trials,trial_conditions_exp,'uni',false)),1);
            curr_std = nanstd(cell2mat(cellfun(@(act,use_trials,cond) ...
                act(use_trials & cond(:,curr_cond),:), ...
                curr_act_pred,pred_use_trials,trial_conditions_exp,'uni',false)),[],1);
            
%             AP_errorfill(t,curr_mean,curr_std,'k',0.5,true);
            
            set(gca,'ColorOrder',brewermap(n_col_bins,'*Greens'));
            plot(t,binned_act_pred_t(binned_act_pred_grp(:,curr_cond + 1) == 1,:)','linewidth',2);
            xlabel('Time'); ylabel('Fluorescence'); 
            title(['Condition ' num2str(curr_cond)]);
            line(repmat(timeavg_t{curr_timeavg}(1),2,1),ylim,'color','k');
            line(repmat(timeavg_t{curr_timeavg}(2),2,1),ylim,'color','k');            
        end
 
        % (plot measured data timecourse)        
        for curr_cond = 1:size(trial_conditions,2)           
            subplot(2,size(trial_conditions,2), ...
                sub2ind(fliplr([2,size(trial_conditions,2)]),curr_cond,2)); hold on;
            set(gca,'ColorOrder',[brewermap(n_col_bins,'*Greys')]);
            plot(t,binned_act_t(binned_act_grp(:,curr_cond + 1) == 1,:)','linewidth',2);
            xlabel('Time'); ylabel('MUA'); 
            title(['Condition ' num2str(curr_cond)]);
            line(repmat(timeavg_t{curr_timeavg}(1),2,1),ylim,'color','k');
            line(repmat(timeavg_t{curr_timeavg}(2),2,1),ylim,'color','k');            
        end        
        
        linkaxes(get(measured_pred_fig,'Children'),'xy');       
        
        % Plot measured v binned predicted
        figure(measured_v_pred_fig)
        
        % (measured vs binned predicted)
        subplot(length(plot_areas),length(timeavg_labels), ...
            sub2ind(fliplr([length(plot_areas),length(timeavg_labels)]),curr_timeavg,curr_area_idx));
        hold on;
        set(gca,'ColorOrder',timeavg_condition_colors{curr_timeavg});

%         fill_cols = min(timeavg_condition_colors{curr_timeavg} + 0.5,1);
%         for curr_cond = 1:size(trial_conditions,2)
%             AP_errorfill( ...
%                 squeeze(nanmean(act_pred_predbinmean(:,:,curr_cond),2)), ...
%                 squeeze(nanmean(act_predfix_predbinmean(:,:,curr_cond),2)), ...
%                 squeeze(AP_sem(act_predfix_predbinmean(:,:,curr_cond),2)), ...
%                 fill_cols(curr_cond,:),1,false);
%         end
        
        errorbar( ...
            squeeze(nanmean(act_pred_predbinmean,2)), ...
            squeeze(nanmean(act_predbinmean,2)), ...
            squeeze(AP_sem(act_predbinmean,2)), ...
            'linestyle','-','linewidth',2,'CapSize',10);
        
        ctx_col = brewermap(n_col_bins,'*Greens');        
        scatter( ...
            reshape(squeeze(nanmean(act_pred_predbinmean,2)),[],1), ...
            reshape(squeeze(nanmean(act_predbinmean,2)),[],1), ...
            60,repmat(ctx_col(1:n_act_bins,:),size(trial_conditions,2),1),'filled');
        
        xlabel(['Fluorescence (' num2str(plot_area) ')']);
        ylabel(['MUA (' num2str(plot_area) ')'])
        axis square tight;
        xlim(xlim + [-0.1,0.1]);
        title(act_title);
        
%         % (plot significant measured condition differences)
%         % (* and o = significant in both measured and "fixed" predicted)
%         curr_ylim = max(ylim);
%         for curr_cond_combo = 1:size(cond_combos,1)
%             % (plot * for measured condition differences)
%             sig_y = curr_ylim + 0.1*curr_cond_combo;
%             sig_x = pred_bin_centers;
%             if any(cond_sig_diff(:,curr_cond_combo))
%                 plot(sig_x(cond_sig_diff(:,curr_cond_combo)),sig_y, ...
%                     '*','MarkerSize',10,'color', ...
%                     timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,1),:));
%                 plot(sig_x(cond_sig_diff(:,curr_cond_combo)),sig_y, ...
%                     '*','MarkerSize',5,'color', ...
%                     timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,2),:));
%             end
%             % (plot o for [predicted condition differences)
%             if any(predfix_cond_sig_diff(:,curr_cond_combo))
%                 plot(sig_x(predfix_cond_sig_diff(:,curr_cond_combo)),sig_y, ...
%                     'o','MarkerSize',15,'color', ...
%                     timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,1),:));
%                 plot(sig_x(predfix_cond_sig_diff(:,curr_cond_combo)),sig_y, ...
%                     'o','MarkerSize',10,'color', ...
%                     timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,2),:));
%             end
%         end
                   
        drawnow;
    end
end

%% Testing fig 4: flip axes (passive)

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(mua_allcat,1),1);

% Get trials with movement during stim to exclude
wheel_thresh = 0.025;
quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > wheel_thresh,2);

% Set windows to average activity
timeavg_labels = {'Stim'};
timeavg_t = {[0,0.2]};
timeavg_align = {stim_align};
% timeavg_trial_conditions = ...
%     {[trial_stim_allcat > 0 & quiescent_trials, ...
%     trial_stim_allcat < 0 & quiescent_trials]};
timeavg_trial_conditions = ...
    {[trial_stim_allcat == 1, ...
    trial_stim_allcat == -1]};
% timeavg_trial_conditions = ...
%     {[trial_stim_allcat > 0, ...
%     trial_stim_allcat < 0]};
% timeavg_trial_conditions = ...
%     {[trial_choice_allcat == 1, ...
%     trial_choice_allcat == -1]};
% timeavg_trial_conditions = ...
%     {[trial_outcome_allcat == 1, ...
%     trial_outcome_allcat == -1]};
timeavg_condition_colors = ...
    {[1,0,0;0,0,1]};

% Set activity percentiles and bins
act_prctile = [5,95];
n_act_bins = 5;

% Set areas and conditions
plot_areas = [1:3];

% Loop across area pairs, plot binned predicted v measured activity
curr_act_allcat = fluor_kernelroi_deconv; % fluor_kernelroi_deconv, mua_ctxpred_allcat, mua_allcat

% (ctx-predicted)
curr_act_pred_allcat = mua_allcat; % fluor_kernelroi_deconv, mua_ctxpred_allcat, mua_allcat

figure('color','w');
for curr_area_idx = 1:length(plot_areas)
    curr_area = plot_areas(curr_area_idx);
  
    [grp_max,grp_idx] = max(timeavg_trial_conditions{1},[],2);
    grp_idx(grp_max == 0) = 0;
    curr_act_avg = nan(max(grp_idx),length(t),length(use_split));
    curr_act_pred_avg = nan(max(grp_idx),length(t),length(use_split));
    for curr_exp = 1:length(use_split)
        for curr_grp = 1:max(grp_idx)
            curr_act_avg(curr_grp,:,curr_exp) = ...
                nanmean(curr_act_allcat(grp_idx == curr_grp & ...
                split_idx == curr_exp & quiescent_trials,:,curr_area),1);
            curr_act_pred_avg(curr_grp,:,curr_exp) = ...
                nanmean(curr_act_pred_allcat(grp_idx == curr_grp & ...
                split_idx == curr_exp & quiescent_trials,:,curr_area),1);
        end
    end

    subplot(length(plot_areas),2,(curr_area_idx-1)*2+1); 
    AP_errorfill(t, ...
            nanmean(curr_act_avg,3)', ...
            AP_sem(curr_act_avg,3)', ...
            timeavg_condition_colors{1}); 
    line([0,0],ylim,'color','k');
    line(repmat(timeavg_t{1}(1),2,1),ylim,'color','g');
    line(repmat(timeavg_t{1}(2),2,1),ylim,'color','g');
    ylabel(['Measured (' num2str(curr_area) ')']);
    xlabel('Time from stim');
    
    subplot(length(plot_areas),2,(curr_area_idx-1)*2+2); 
    AP_errorfill(t, ...
            nanmean(curr_act_pred_avg,3)', ...
            AP_sem(curr_act_pred_avg,3)', ...
            timeavg_condition_colors{1}); 
    line([0,0],ylim,'color','k');
    line(repmat(timeavg_t{1}(1),2,1),ylim,'color','g');
    line(repmat(timeavg_t{1}(2),2,1),ylim,'color','g');
    ylabel(['Predicted (' num2str(curr_area) ')']);
    xlabel('Time from stim');
    
  
end
linkaxes(get(gcf,'Children'),'xy');


measured_v_pred_fig = figure('color','w');
for curr_area_idx = 1:length(plot_areas)   
    
    plot_area = plot_areas(curr_area_idx);  
    
    for curr_timeavg = 1:length(timeavg_labels)
        
        trial_conditions = timeavg_trial_conditions{curr_timeavg};
        
        % (re-align and split activity)
        act_title = timeavg_labels{curr_timeavg};
        
        curr_act = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        curr_act_pred = mat2cell(...
            cell2mat(arrayfun(@(trial) circshift( ...
            curr_act_pred_allcat(trial,:,plot_area), ...
            timeavg_align{curr_timeavg}(trial),2),transpose(1:size(curr_act_pred_allcat,1)),'uni',false)), ...
            use_split,length(t));
        
        trial_conditions_exp = mat2cell(trial_conditions,use_split,size(trial_conditions,2));
        
        % (get average activity within window)
        curr_event_t = t >= timeavg_t{curr_timeavg}(1) & t <= timeavg_t{curr_timeavg}(2);
        curr_act_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act,'uni',false);
        curr_act_pred_avg = cellfun(@(x) squeeze(nanmean(x(:,curr_event_t,:),2)),curr_act_pred,'uni',false);
        
        % (bin predicted data across percentile range)
        pred_bin_edges = prctile(cell2mat(curr_act_pred_avg),linspace(act_prctile(1),act_prctile(2),n_act_bins+1));
        pred_bin_centers = pred_bin_edges(1:end-1) + diff(pred_bin_edges)./2;        
        pred_trial_bins = cellfun(@(x) discretize(x,pred_bin_edges),curr_act_pred_avg,'uni',false);         
        
        % (get activity binned by predicted)        
        pred_use_trials = cellfun(@(act,act_pred,trial_bins) ...
            squeeze(~any(isnan(act(:,curr_event_t)),2)) & ...
            squeeze(~any(isnan(act_pred(:,curr_event_t)),2)) & ...
            ~isnan(trial_bins), ...
            curr_act,curr_act_pred,pred_trial_bins,'uni',false);
        
        act_predbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));
        
        act_pred_predbinmean = cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_pred_avg,pred_trial_bins,trial_conditions_exp,pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false));    
        
        % Get condition difference significance from shuffle
        n_shuff = 1;
        trial_conditions_shuff = cellfun(@(x) AP_shake(x,1), ...
            repmat(trial_conditions_exp,1,n_shuff),'uni',false);
        
        act_predbinmean_condshuff = cell2mat(arrayfun(@(curr_shuff) ...
            cell2mat(arrayfun(@(condition) ...
            cell2mat(cellfun(@(act,bins,trial_cond,use_trials) ...
            accumarray(bins(use_trials & trial_cond(:,condition)), ...
            act(use_trials & trial_cond(:,condition)), ...
            [n_act_bins,1],@nanmean,cast(NaN,class(act))), ...
            curr_act_avg,pred_trial_bins,trial_conditions_shuff(:,curr_shuff), ...
            pred_use_trials,'uni',false)'), ...
            permute(1:size(trial_conditions,2),[1,3,2]),'uni',false)), ...
            permute(1:n_shuff,[1,3,4,2]),'uni',false));
        
        % (measured data: null = no difference between conditions)
        cond_combos = nchoosek(1:size(trial_conditions,2),2);
        cond_sig_diff = false(n_act_bins,size(cond_combos,1));
        for curr_cond_combo = 1:size(cond_combos,1)
            curr_combo_diff = nanmean(act_predbinmean(:,:,cond_combos(curr_cond_combo,1)) - ...
                 act_predbinmean(:,:,cond_combos(curr_cond_combo,2)),2);
            shuff_prctile = squeeze(prctile(nanmean(act_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,1),:) - ...
                act_predbinmean_condshuff(:,:,cond_combos(curr_cond_combo,2),:),2),[2.5,97.5],4));
            cond_sig_diff(:,curr_cond_combo) = ...
                curr_combo_diff < shuff_prctile(:,1) | ...
                curr_combo_diff > shuff_prctile(:,2);
        end 
        
        % Plot measured v predicted in bins
        figure(measured_v_pred_fig)
        
        % (measured vs binned predicted)
        subplot(length(plot_areas),length(timeavg_labels), ...
            sub2ind(fliplr([length(plot_areas),length(timeavg_labels)]),curr_timeavg,curr_area_idx));
        hold on;
        set(gca,'ColorOrder',timeavg_condition_colors{curr_timeavg});

        fill_cols = min(timeavg_condition_colors{curr_timeavg} + 0.5,1);
%         for curr_cond = 1:size(trial_conditions,2)
%             AP_errorfill( ...
%                 squeeze(nanmean(act_pred_predbinmean(:,:,curr_cond),2)), ...
%                 squeeze(nanmean(act_predfix_predbinmean(:,:,curr_cond),2)), ...
%                 squeeze(AP_sem(act_predfix_predbinmean(:,:,curr_cond),2)), ...
%                 fill_cols(curr_cond,:),1,false);
%         end
        
        errorbar( ...
            squeeze(nanmean(act_pred_predbinmean,2)), ...
            squeeze(nanmean(act_predbinmean,2)), ...
            squeeze(AP_sem(act_predbinmean,2)), ...
            'linestyle','-','linewidth',2,'CapSize',10);
        
        n_col_bins = n_act_bins + 2;
        ctx_col = brewermap(n_col_bins,'*Greens');        
        scatter( ...
            reshape(squeeze(nanmean(act_pred_predbinmean,2)),[],1), ...
            reshape(squeeze(nanmean(act_predbinmean,2)),[],1), ...
            60,repmat(ctx_col(1:n_act_bins,:),size(trial_conditions,2),1),'filled');
        
        xlabel(['Striatum (' num2str(plot_area) ')']);
        ylabel(['Cortex (' num2str(plot_area) ')'])
        axis square;
        title(act_title);
        
        % (plot significant measured condition differences)
        % (* and o = significant in both measured and "fixed" predicted)
        curr_ylim = max(ylim);
        for curr_cond_combo = 1:size(cond_combos,1)
            % (plot * for measured condition differences)
            sig_y = curr_ylim + 0.1*curr_cond_combo;
            sig_x = pred_bin_centers;
            if any(cond_sig_diff(:,curr_cond_combo))
                plot(sig_x(cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                    '*','MarkerSize',10,'color', ...
                    timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,1),:));
                plot(sig_x(cond_sig_diff(:,curr_cond_combo)),sig_y, ...
                    '*','MarkerSize',5,'color', ...
                    timeavg_condition_colors{curr_timeavg}(cond_combos(curr_cond_combo,2),:));
            end
        end
        
        drawnow;
        
    end
end

linkaxes(get(measured_v_pred_fig,'Children'),'xy');




% Plot distributions
curr_act_avg_cat = cell2mat(curr_act_avg);
curr_act_pred_avg_cat = cell2mat(curr_act_pred_avg);

t_smooth_n = sum(t >= timeavg_t{1}(1) & t <= timeavg_t{1}(2));
t_smooth = ones(t_smooth_n,1)./t_smooth_n;

figure;
subplot(2,2,1);
AP_heatscatter(conv(reshape(curr_act_pred_allcat(:,t < 0,plot_areas),[],1),t_smooth,'same'), ...
    conv(reshape(curr_act_allcat(:,t < 0,plot_areas),[],1),t_smooth,'same'),100);
line(xlim,xlim)
xlabel('MUA');
ylabel('Fluorescence');
title('ITI points')

use_bins = linspace(-2,6,100);
use_bin_centers = conv(use_bins,[0.5,0.5],'valid');
act_dist_1 = histcounts2(curr_act_pred_avg_cat(timeavg_trial_conditions{1}(:,1)), ...
    curr_act_avg_cat(timeavg_trial_conditions{1}(:,1)),use_bins,use_bins);
act_dist_2 = histcounts2(curr_act_pred_avg_cat(timeavg_trial_conditions{1}(:,2)), ...
    curr_act_avg_cat(timeavg_trial_conditions{1}(:,2)),use_bins,use_bins);
act_dist_iti = histcounts2( ...
    conv(reshape(curr_act_pred_allcat(:,t < 0 | t > 1.5,plot_areas),[],1),t_smooth,'same'), ...
    conv(reshape(curr_act_allcat(:,t < 0 | t > 1.5,plot_areas),[],1),t_smooth,'same'), ...
    use_bins,use_bins);

subplot(2,2,2); hold on;
imagesc(use_bin_centers,use_bin_centers,imgaussfilt(act_dist_1',3)-imgaussfilt(act_dist_2',3));
line(xlim,xlim)
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(gca,brewermap([],'*RdBu'));
title('Stim act distr difference');

subplot(2,2,3); hold on;
imagesc(use_bin_centers,use_bin_centers,act_dist_1');
line(xlim,xlim)
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(gca,brewermap([],'*RdBu'));
title('Contra distribution');

subplot(2,2,4); hold on;
imagesc(use_bin_centers,use_bin_centers,act_dist_2');
line(xlim,xlim)
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(gca,brewermap([],'*RdBu'));
title('Ipsi distribution');

linkaxes(get(gcf,'Children'),'xy');


figure;

% stim_all = unique([0,0.06,0.125,0.25,0.5,1].*[-1;1]);
% unique_stim = unique(trial_stim_allcat);
% col_all = brewermap(11,'*RdBu');
% [~,stim_idx] = ismember(unique_stim,stim_all,'rows');
% col = col_all(stim_idx,:);

[~,grp_idx] = max(timeavg_trial_conditions{1},[],2);
scatterhist(curr_act_pred_avg_cat(grp_idx ~= 0), ...
    curr_act_avg_cat(grp_idx ~= 0), ...
    'Group',grp_idx(grp_idx ~= 0), ...
    'Kernel','on','Color',timeavg_condition_colors{1},'Marker','.');
xlabel('Striatum');
ylabel('Cortex-predicted');


figure;
scatterhist(conv(reshape(curr_act_pred_allcat(:,t < 0,plot_areas),[],1),t_smooth,'same'), ...
    conv(reshape(curr_act_allcat(:,t < 0,plot_areas),[],1),t_smooth,'same'), ...
    'Kernel','on','Color','k','Marker','.');
xlabel('Striatum');
ylabel('Cortex-predicted');
title('ITI activity');


figure;
imagesc(use_bin_centers,use_bin_centers,act_dist_iti');
line(xlim,xlim)
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(gca,brewermap([],'*RdBu'));
set(gca,'YDir','normal')
xlabel('Striatum');
ylabel('Cortex');
title('ITI distribution');

        

% Fit a line between cortex ROI and striatum
% (iti)
ap_memmap = conv(reshape(curr_act_pred_allcat(:,t < 0,plot_areas),[],1),t_smooth,'same');
b = conv(reshape(curr_act_allcat(:,t < 0,plot_areas),[],1),t_smooth,'same');


nanpoints = isnan(ap_memmap) | isnan(b);

% (fit each axis separately)
c1 = [ap_memmap(~nanpoints),ones(sum(~nanpoints),1)]\b(~nanpoints);
c2 = [b(~nanpoints),ones(sum(~nanpoints),1)]\ap_memmap(~nanpoints);
% (fit with PCA)
c3 = pca([ap_memmap,b]);

figure; 
subplot(1,2,1);
AP_heatscatter(ap_memmap,b,100);
line(xlim,xlim*c1(1)+c1(2),'color','b','linewidth',2);
line(ylim*c2(1)+c2(2),ylim,'color','b','linewidth',2);
line(xlim,xlim*(c3(2)/c3(1)),'color','g','linewidth',2)

title('ITI activity and fit lines');

% (stim)
ap_memmap = nanmean(curr_act_pred_allcat(trial_stim_allcat == 1,t > 0 & t < 0.15,plot_areas),2);
b = nanmean(curr_act_allcat(trial_stim_allcat == 1,t > 0 & t < 0.15,plot_areas),2);

nanpoints = isnan(ap_memmap) | isnan(b);

% (fit each axis separately)
c1 = [ap_memmap(~nanpoints),ones(sum(~nanpoints),1)]\b(~nanpoints);
c2 = [b(~nanpoints),ones(sum(~nanpoints),1)]\ap_memmap(~nanpoints);
% (fit with PCA)
c3 = pca([ap_memmap,b]);

subplot(1,2,2);
AP_heatscatter(ap_memmap,b,50);
line(xlim,xlim*c1(1)+c1(2),'color','b','linewidth',2);
line(ylim*c2(1)+c2(2),ylim,'color','b','linewidth',2);
line(xlim,xlim*(c3(2)/c3(1)),'color','g','linewidth',2)

title('Stim activity and fit lines');

linkaxes(get(gcf,'Children'),'xy');



% combining binning on each axis?
ap_memmap = reshape(nanmean(curr_act_pred_allcat(trial_stim_allcat == 1,t > 0 & t < 0.2,plot_areas),2),[],1);
b = reshape(nanmean(curr_act_allcat(trial_stim_allcat == 1,t > 0 & t < 0.2,plot_areas),2),[],1);

% c = reshape(nanmean(curr_act_pred_allcat(trial_stim_allcat == -1,t > 0 & t < 0.2,plot_areas),2),[],1);
% d = reshape(nanmean(curr_act_allcat(trial_stim_allcat == -1,t > 0 & t < 0.2,plot_areas),2),[],1);
c = convn(reshape(curr_act_pred_allcat(:,t < 0,plot_areas),[],1),t_smooth,'same');
d = convn(reshape(curr_act_allcat(:,t < 0,plot_areas),[],1),t_smooth,'same');

use_bins = linspace(-1,2,5);
bin_centers = (use_bins(1:end-1) + diff(use_bins)/2);

a_bins = discretize(ap_memmap,use_bins);
b_bins = discretize(b,use_bins);
a_b_bins = accumarray(b_bins(~isnan(b_bins)),ap_memmap(~isnan(b_bins)),[length(bin_centers),1],@nanmean);
b_a_bins = accumarray(a_bins(~isnan(a_bins)),b(~isnan(a_bins)),[length(bin_centers),1],@nanmean);


c_bins = discretize(c,use_bins);
d_bins = discretize(d,use_bins);
c_d_bins = accumarray(d_bins(~isnan(d_bins)),c(~isnan(d_bins)),[length(bin_centers),1],@nanmean);
d_c_bins = accumarray(c_bins(~isnan(c_bins)),d(~isnan(c_bins)),[length(bin_centers),1],@nanmean);

figure;
subplot(1,2,1); hold on
plot(c,d,'.k');
plot(ap_memmap,b,'.r');
plot(bin_centers,d_c_bins,'linewidth',2);
plot(c_d_bins,bin_centers,'linewidth',2);
plot(bin_centers,b_a_bins,'linewidth',2);
plot(a_b_bins,bin_centers,'linewidth',2);
subplot(1,2,2); hold on
plot(c_d_bins,d_c_bins,'linewidth',2);
plot(a_b_bins,b_a_bins,'linewidth',2);



%% Plot average image across recordings

animals = {'AP045','AP054','AP055','AP053','AP047','AP048'};

avg_im_cat = [];
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    protocol = 'AP_sparseNoise';
    experiments = AP_find_experiments(animal,protocol);
    experiments(~([experiments.imaging] & [experiments.ephys])) = [];
    
    for curr_day = 1:length(experiments)
        day = experiments(curr_day).day;
        [img_path,img_exists] = AP_cortexlab_filename(animal,day,[],'imaging');
        avg_im = readNPY([img_path filesep 'meanImage_blue.npy']);
        avg_im_cat = cat(3,avg_im_cat,AP_align_widefield(avg_im,animal,day));
    end
    
    AP_print_progress_fraction(curr_animal,length(animals));
end

AP_imscroll(avg_im_cat)
axis image;





%% Trying "VFS" of pixel correlation

wf_corr_borders_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_borders\wf_corr_borders.mat';
load(wf_corr_borders_fn);

% Spacing/downsampling (hardcoded - in preprocessing)
px_spacing = 20;
downsample_factor = 10;

% Get average correlation maps
wf_corr_map_recording = [wf_corr_borders(:).corr_map_downsamp];
wf_corr_map_cat = cat(3,wf_corr_map_recording{:});
wf_corr_map_mean = cell(size(wf_corr_map_cat,1),size(wf_corr_map_cat,2));
for curr_y = 1:size(wf_corr_map_cat,1)
    for curr_x = 1:size(wf_corr_map_cat,2)
        wf_corr_map_mean{curr_y,curr_x} = ...
            nanmean(cell2mat(wf_corr_map_cat(curr_y,curr_x,:)),3);
    end
end

[ap_memmap,b] = cellfun(@imgradient,wf_corr_map_mean,'uni',false);
a_upsample = imresize(nanmean(cat(3,ap_memmap{:}),3),downsample_factor,'bilinear');


[ap_memmap,b] = cellfun(@imgradientxy,wf_corr_map_mean,'uni',false);



[Uy,Ux] = size(wf_corr_map_mean{1});

[ny,nx] = size(wf_corr_map_mean);

use_method = 'com'; % max or com

stim_im_smoothed = reshape(permute(cat(3,wf_corr_map_mean{:}),[3,1,2]),ny,nx,[]);
stim_im_smoothed(isnan(stim_im_smoothed)) = 0;

% Conversely, do COM on original^2
[xx,yy] = meshgrid(1:size(stim_im_smoothed,2),1:size(stim_im_smoothed,1));
m_xr = reshape(sum(sum(bsxfun(@times,stim_im_smoothed.^2,xx),1),2)./sum(sum(stim_im_smoothed.^2,1),2),Uy,Ux);
m_yr = reshape(sum(sum(bsxfun(@times,stim_im_smoothed.^2,yy),1),2)./sum(sum(stim_im_smoothed.^2,1),2),Uy,Ux);

% Calculate and plot sign map (dot product between horz & vert gradient)

% 1) get gradient direction
[~,Vdir] = imgradient(m_yr);
[~,Hdir] = imgradient(m_xr);

% 3) get sin(difference in direction) if retinotopic, H/V should be
% orthogonal, so the closer the orthogonal the better (and get sign)
angle_diff = sind(Vdir-Hdir);
angle_diff(isnan(angle_diff)) = 0;

vfs = angle_diff;

vfs_upsample = imresize(vfs,downsample_factor,'nearest');
figure;imagesc(vfs_upsample);
colormap(brewermap([],'*RdBu'));
caxis([-1,1]);
axis image off
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

%% Matteo: bilaterality by area

wf_corr_borders_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_borders\wf_corr_borders.mat';
load(wf_corr_borders_fn);

% Spacing/downsampling (hardcoded - in preprocessing)
px_spacing = 20;
downsample_factor = 10;

% Get average correlation maps
wf_corr_map_recording = [wf_corr_borders(:).corr_map_downsamp];
wf_corr_map_cat = cat(3,wf_corr_map_recording{:});
wf_corr_map_mean = cell(size(wf_corr_map_cat,1),size(wf_corr_map_cat,2));
for curr_y = 1:size(wf_corr_map_cat,1)
    for curr_x = 1:size(wf_corr_map_cat,2)
        wf_corr_map_mean{curr_y,curr_x} = ...
            nanmean(cell2mat(wf_corr_map_cat(curr_y,curr_x,:)),3);
    end
end


px_y = (0:size(wf_corr_map_mean,1)-1)*px_spacing+1;
px_x = (0:size(wf_corr_map_mean,2)-1)*px_spacing+1;
ipsi_corr = nan(461,440);
contra_corr = nan(461,440);

for curr_y_idx = 1:length(px_y)
    for curr_x_idx = 1:length(px_x)
        curr_corr_upsample = ...
            imresize(wf_corr_map_mean{curr_y_idx,curr_x_idx}, ...
            downsample_factor,'nearest');
        curr_corr_upsample(isnan(curr_corr_upsample)) = 0;
        
        curr_mask = false(size(curr_corr_upsample));
        curr_mask(px_y(curr_y_idx),px_x(curr_x_idx)) = true;
        curr_mask_reflect = AP_reflect_widefield(+curr_mask) > 0;
        
        if ~any(curr_mask_reflect(:))
            continue
        end
        
        ipsi_corr(px_y(curr_y_idx),px_x(curr_x_idx)) = curr_corr_upsample(curr_mask);
        contra_corr(px_y(curr_y_idx),px_x(curr_x_idx)) = ...
            (curr_mask_reflect(:)'*curr_corr_upsample(:))./sum(curr_mask_reflect(:));
    end
end

ipsi_corr_inpaint = inpaint_nans(ipsi_corr);
contra_corr_inpaint = inpaint_nans(contra_corr);

figure
subplot(1,2,1);
imagesc(ipsi_corr_inpaint);
axis image;
caxis([0,1]);
AP_reference_outline('ccf_aligned','k');
subplot(1,2,2);
imagesc(contra_corr_inpaint);
axis image
caxis([0,1]);
AP_reference_outline('ccf_aligned','k');

figure;
imagesc(contra_corr_inpaint);
axis image off;
AP_reference_outline('ccf_aligned','k');
c = colorbar('FontSize',12);
ylabel(c,'Correlation with contra','FontSize',12);
colormap(brewermap([],'Greys'));
caxis([0.5,1]);

%% Bilateral correlation (like above but within experiment)

fVdf_deconv = AP_deconv_wf(fVdf);
aUdf = AP_align_widefield(Udf,animal,day);

% Get V covariance
Ur = reshape(aUdf, size(aUdf,1)*size(aUdf,2),[]); % P x S
covV = cov(fVdf_deconv'); % S x S % this is the only one that takes some time really
varP = dot((Ur*covV)', Ur'); % 1 x P

% Loop through pixels, correlate with mirror pixel
px_mirror_idx = round(AP_reflect_widefield(reshape(1:size(Ur,1), ...
    size(aUdf,1),size(aUdf,2))));

px_corr = NaN(size(aUdf,1),size(aUdf,2));

figure;
c = imagesc(px_corr);
colormap(brewermap([],'*RdBu'));
caxis([-1,1]);
for curr_px = 1:100:size(Ur,1)
    % Get mirror pixel (doesn't work just using mirrored indicies? why??)
    curr_px_mirror = px_mirror_idx(curr_px);
    
%     curr_mask = false(size(aUdf,1),size(aUdf,2));
%     curr_mask(curr_px) = true;
%     curr_px_mirror = find(AP_reflect_widefield(+curr_mask),1);
    
    % If mirrored pixel not imaged, skip
    if isempty(curr_px_mirror) || curr_px_mirror == 0
        continue
    end
    
    covP = Ur(curr_px,:)*covV*Ur(curr_px_mirror,:)';
    stdPxPy = varP(curr_px).^0.5 * varP(curr_px_mirror).^0.5;
    px_corr(curr_px) = covP./stdPxPy;
    
    set(c,'CData',px_corr);
    drawnow;
    
    AP_print_progress_fraction(curr_px,size(Ur,1))
end


% (old: gets full pixel correlation map, not just bilateral)

ySize = size(aUdf,1); xSize = size(aUdf,2);

px_spacing = 20;
use_y = 1:px_spacing:size(aUdf,1);
use_x = 1:px_spacing:size(aUdf,2);
corr_map = cell(length(use_y),length(use_x));
for curr_x_idx = 1:length(use_x)
    curr_x = use_x(curr_x_idx);
    for curr_y_idx = 1:length(use_y)
        curr_y = use_y(curr_y_idx);
        
        pixel = [curr_y,curr_x];
        pixelInd = sub2ind([ySize, xSize], pixel(1), pixel(2));
        
        covP = Ur(pixelInd,:)*covV*Ur'; % 1 x P
        stdPxPy = varP(pixelInd).^0.5 * varP.^0.5; % 1 x P
        corrMat = reshape(covP./stdPxPy,ySize,xSize); % 1 x P
        
        corr_map{curr_y_idx,curr_x_idx} = corrMat;
    end
    AP_print_progress_fraction(curr_x_idx,length(use_x))
end

% Get downsampled correlation maps
downsample_factor = px_spacing;
corr_map_downsamp = cellfun(@(x) ...
    imresize(x,1/downsample_factor,'bilinear'),corr_map,'uni',false);

wf_corr_map_mean = corr_map_downsamp;

%%% (FROM ABOVE CELL)
px_y = (0:size(wf_corr_map_mean,1)-1)*px_spacing+1;
px_x = (0:size(wf_corr_map_mean,2)-1)*px_spacing+1;
ipsi_corr = nan(461,440);
contra_corr = nan(461,440);

for curr_y_idx = 1:length(px_y)
    for curr_x_idx = 1:length(px_x)
        curr_corr_upsample = ...
            imresize(wf_corr_map_mean{curr_y_idx,curr_x_idx}, ...
            downsample_factor,'nearest');
        curr_corr_upsample(isnan(curr_corr_upsample)) = 0;
        
        curr_mask = false(size(curr_corr_upsample));
        curr_mask(px_y(curr_y_idx),px_x(curr_x_idx)) = true;
        curr_mask_reflect = AP_reflect_widefield(+curr_mask) > 0;
        
        if ~any(curr_mask_reflect(:))
            continue
        end
        
        ipsi_corr(px_y(curr_y_idx),px_x(curr_x_idx)) = curr_corr_upsample(curr_mask);
        contra_corr(px_y(curr_y_idx),px_x(curr_x_idx)) = ...
            (curr_mask_reflect(:)'*curr_corr_upsample(:))./sum(curr_mask_reflect(:));
    end
end

ipsi_corr_inpaint = inpaint_nans(ipsi_corr);
contra_corr_inpaint = inpaint_nans(contra_corr);

figure
subplot(1,2,1);
imagesc(ipsi_corr_inpaint);
axis image;
caxis([0,1]);
AP_reference_outline('ccf_aligned','k');
subplot(1,2,2);
imagesc(contra_corr_inpaint);
axis image
caxis([0,1]);
AP_reference_outline('ccf_aligned','k');

figure;
imagesc(contra_corr_inpaint);
axis image off;
AP_reference_outline('ccf_aligned','k');
c = colorbar('FontSize',12);
ylabel(c,'Correlation with contra','FontSize',12);
colormap(brewermap([],'Greys'));
caxis([0.5,1]);





%% Wheel move (any) aligned activity

% Align U and deconv fluor
aUdf = AP_align_widefield(Udf,animal,day);
fVdf_deconv = AP_deconv_wf(fVdf);

% Get fluorescence by ROI
load('C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_alignment\U_master');
wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);
fluor_roi_deconv = AP_svd_roi(aUdf,fVdf_deconv,[],[],cat(3,wf_roi.mask));

% Get movement onsets (0 velocity q time, > thresh for 80% of m time)
thresh_displacement = 0.025;
wheel_quiescence = 0.5; % quiescent for this long
wheel_move = 0.2; % moving for this long

sample_rate = Timeline.hw.daqSampleRate;
wheel_quiescence_samples = round(wheel_quiescence*sample_rate);
wheel_move_samples = round(wheel_move*sample_rate);

wheel_quiescent_thresh_full = conv(wheel_velocity == 0, ...
    ones(1,wheel_quiescence_samples)) == wheel_quiescence_samples;
wheel_quiescent_thresh = wheel_quiescent_thresh_full(1:size(wheel_velocity));

wheel_move_thresh_full = conv(abs(wheel_velocity) > thresh_displacement, ...
    ones(1,wheel_move_samples)) > wheel_move_samples*0.8;
wheel_move_thresh = wheel_move_thresh_full(end-size(wheel_velocity)+1:end);


% Make move-aligned MUA
mua_allcat_move = mua_allcat;
mua_ctxpred_allcat_move = mua_ctxpred_allcat;
mua_taskpred_allcat_move = mua_taskpred_allcat;
mua_taskpred_reduced_allcat_move = mua_taskpred_reduced_allcat;
mua_ctxpred_taskpred_allcat_move = mua_taskpred_allcat;
mua_ctxpred_taskpred_reduced_allcat_move = mua_taskpred_reduced_allcat;

t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
for i = 1:size(mua_allcat,1)
    mua_allcat_move(i,:,:) = circshift(mua_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
    mua_ctxpred_allcat_move(i,:,:) = circshift(mua_ctxpred_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
    mua_taskpred_allcat_move(i,:,:) = circshift(mua_taskpred_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
    mua_taskpred_reduced_allcat_move(i,:,:,:) = circshift(mua_taskpred_reduced_allcat_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
    mua_ctxpred_taskpred_allcat_move(i,:,:) = circshift(mua_ctxpred_taskpred_allcat_move(i,:,:),-move_idx(i)+leeway_samples,2);
    mua_ctxpred_taskpred_reduced_allcat_move(i,:,:,:) = circshift(mua_ctxpred_taskpred_reduced_allcat_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
end



% Make move-aligned fluorescence
fluor_allcat_deconv_move = fluor_allcat_deconv;
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
for i = 1:size(fluor_allcat_deconv,1)
    fluor_allcat_deconv_move(i,:,:,:) = circshift(fluor_allcat_deconv_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
end


%% Load natural images
im_dir = '\\zserver.cortexlab.net\Data\pregenerated_textures\JulieF\shapesAndNatImages';
n_im = 50;
nat_im = nan(375,500,n_im);
for curr_im_idx = 1:n_im
   curr_im = load([im_dir filesep 'img' num2str(curr_im_idx) '.mat']);
   nat_im(:,:,curr_im_idx) = curr_im.img;
   AP_print_progress_fraction(curr_im_idx,n_im);
end



%% Average stim-aligned cortex

% Get average stim-aligned fluorescence 
plot_trials = trial_stim_allcat > 0 & trial_choice_allcat == 1;
plot_trials_exp = mat2cell(plot_trials,use_split,1);

fluor_allcat_deconv_exp = mat2cell(fluor_allcat_deconv,use_split,length(t),n_vs);

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = zeros(size(trial_stim_allcat));
move_align = -move_idx + leeway_samples;
outcome_align = -outcome_idx + leeway_samples;

% Set windows to average activity
use_align_labels = {'Stim','Move onset','Outcome'};
use_align = {stim_align};
plot_t = [0.08,0,0.08];

for curr_align = 1:length(use_align)
    
    % (re-align activity)
    curr_ctx_act = cellfun(@(act,trials,shift) cell2mat(arrayfun(@(trial) ...
        circshift(act(trial,:,:),shift(trial),2), ...
        find(trials),'uni',false)), ...
        fluor_allcat_deconv_exp,plot_trials_exp, ...
        mat2cell(use_align{curr_align},use_split,1),'uni',false);
    
    curr_ctx_act_mean = ...
        permute(nanmean(cell2mat(cellfun(@(x) nanmean(x,1), ...
        curr_ctx_act,'uni',false)),1),[3,2,1]);
    
    curr_ctx_act_mean_px = svdFrameReconstruct(U_master(:,:,1:n_vs), ...
        curr_ctx_act_mean);
    
    AP_imscroll(curr_ctx_act_mean_px,t);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    axis image off;
    colormap(brewermap([],'PRGn'));
    caxis([-0.02,0.02]);
    title([use_align_labels{curr_align} ': ' num2str(plot_t(curr_align)) ' sec']);
    
end



%% Testing as a check: how well does V1 activity predict AM activity?

% Regress kernel ROI activity to striatum domain activity (per recording)
regression_params.kernel_t = [-0.1,0.1];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;
lambda = 0;
kernel_frames = floor(regression_params.kernel_t(1)*sample_rate): ...
    ceil(regression_params.kernel_t(2)*sample_rate);

% Set time to use
% (e.g. this type of code can be used to regress only ITI time)
use_t = true(size(t));

fluor_roi_deconv_exp = mat2cell(fluor_roi_deconv,trials_recording,length(t),numel(wf_roi));

fluorpred_k = nan(length(use_split),length(kernel_frames));
fluorpred_exp = cellfun(@(x) nan(size(x(:,:,1))),fluor_roi_deconv_exp,'uni',false);
for curr_exp = 1:length(trials_recording)

    regressor_roi = 1;
    predict_roi = 3;
    
    curr_mua = reshape(fluor_roi_deconv_exp{curr_exp}(:,use_t,predict_roi)',1,[]);
    curr_mua_std = curr_mua./nanstd(curr_mua);
    
    curr_fluor = reshape(permute( ...
        fluor_roi_deconv_exp{curr_exp}(:,use_t,regressor_roi),[2,1,3]),[],1)';
    curr_fluor_std = curr_fluor./nanstd(curr_fluor);

    % Set discontinuities in trial data
    trial_discontinuities = false(size(fluor_roi_deconv_exp{curr_exp}(:,use_t,1)));
    trial_discontinuities(:,1) = true;
    trial_discontinuities = reshape(trial_discontinuities',[],1)';
    
    % Full cortex regression
    [ctx_str_k,curr_mua_fluorpred_std,explained_var_trial] = ...
        AP_regresskernel(curr_fluor_std, ...
        curr_mua_std,kernel_frames,lambda, ...
        regression_params.zs,regression_params.cvfold, ...
        true,regression_params.use_constant,trial_discontinuities);
    
    % Re-scale the prediction (subtract offset, multiply, add scaled offset)
    ctxpred_spikes = (curr_mua_fluorpred_std - squeeze(ctx_str_k{end})).* ...
        nanstd(curr_mua,[],2) + ...
        nanstd(curr_mua,[],2).*squeeze(ctx_str_k{end});
        
    fluorpred_k(curr_exp,:) = ctx_str_k{1};
    fluorpred_exp{curr_exp}(:,use_t) = ...
        reshape(ctxpred_spikes,sum(use_t),[])';
    
    AP_print_progress_fraction(curr_exp,length(trials_recording));
      
end


% Get R^2 for task, cortex full, and cortex ROI predictions
fluor_roi_taskpred_exp = mat2cell(fluor_roi_taskpred,trials_recording,length(t),numel(wf_roi));

ctxpred_r2 = nan(max(split_idx),1);
ctxtaskpred_r2 = nan(max(split_idx),1);

for curr_exp = 1:max(split_idx)
    
    curr_data = reshape(permute(fluor_roi_deconv_exp{curr_exp}(:,:,predict_roi),[2,1,3]),[],1);
    curr_ctxpred_data = reshape(permute(fluorpred_exp{curr_exp},[2,1,3]),[],1);
    curr_ctxtaskpred_data = reshape(permute(fluor_roi_taskpred_exp{curr_exp}(:,:,predict_roi),[2,1,3]),[],1);

    % Set common NaNs
    nan_samples = isnan(curr_data) | isnan(curr_ctxpred_data) ...
        | isnan(curr_ctxtaskpred_data);
    curr_data(nan_samples) = NaN;
    curr_ctxpred_data(nan_samples) = NaN;
    curr_ctxtaskpred_data(nan_samples) = NaN;
    
    ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    ctxtaskpred_r2(curr_exp,:) = 1 - (nansum((curr_data-curr_ctxtaskpred_data).^2,1)./ ...
        nansum((curr_data-nanmean(curr_data,1)).^2,1));
    
end

figure; plot(ctxtaskpred_r2,ctxpred_r2,'.k')
xlabel('Task-predicted')
ylabel(['Cortex-predicted (from ' wf_roi(regressor_roi).area ')'])
title(wf_roi(predict_roi).area)
axis square
xlim([min([xlim,ylim]),max([xlim,ylim])]);
ylim(xlim);
line(xlim,xlim);

%% (from above: BUT PREDICT MUA FROM FLUOR ROI)

% Regress kernel ROI activity to striatum domain activity (per recording)
regression_params.kernel_t = [-0.1,0.1];
regression_params.zs = [false,false];
regression_params.cvfold = 5;
regression_params.use_constant = true;
lambda = 0;
kernel_frames = floor(regression_params.kernel_t(1)*sample_rate): ...
    ceil(regression_params.kernel_t(2)*sample_rate);

% Set time to use
% (e.g. this type of code can be used to regress only ITI time)
use_t = true(size(t));

mua_exp = vertcat(mua_all{:});

fluorpred_k = nan(length(use_split),length(kernel_frames));
fluorpred_exp = cellfun(@(x) nan(size(x(:,:,1))),fluor_roi_deconv_exp,'uni',false);
for curr_exp = 1:length(trials_recording)

    regressor_roi = 3;
    predict_roi = 1;
    
    curr_mua = reshape(mua_exp{curr_exp}(:,use_t,predict_roi)',1,[]);
    curr_mua_std = curr_mua./nanstd(curr_mua);
    
    curr_fluor = reshape(permute( ...
        fluor_roi_deconv_exp{curr_exp}(:,use_t,regressor_roi),[2,1,3]),[],1)';
    curr_fluor_std = curr_fluor./nanstd(curr_fluor);

    % Set discontinuities in trial data
    trial_discontinuities = false(size(fluor_roi_deconv_exp{curr_exp}(:,use_t,1)));
    trial_discontinuities(:,1) = true;
    trial_discontinuities = reshape(trial_discontinuities',[],1)';
    
    % Full cortex regression
    [ctx_str_k,curr_mua_fluorpred_std,explained_var_trial] = ...
        AP_regresskernel(curr_fluor_std, ...
        curr_mua_std,kernel_frames,lambda, ...
        regression_params.zs,regression_params.cvfold, ...
        true,regression_params.use_constant,trial_discontinuities);
    
    % Re-scale the prediction (subtract offset, multiply, add scaled offset)
    ctxpred_spikes = (curr_mua_fluorpred_std - squeeze(ctx_str_k{end})).* ...
        nanstd(curr_mua,[],2) + ...
        nanstd(curr_mua,[],2).*squeeze(ctx_str_k{end});
        
    fluorpred_k(curr_exp,:) = ctx_str_k{1};
    fluorpred_exp{curr_exp}(:,use_t) = ...
        reshape(ctxpred_spikes,sum(use_t),[])';
    
    AP_print_progress_fraction(curr_exp,length(trials_recording));
      
end

% Get R^2 for task and cortex
mua_taskpred_exp = vertcat(mua_taskpred_all{:});

taskpred_r2 = nan(length(mua_exp),1);
ctxpred_r2 = nan(length(mua_exp),1);
for curr_exp = 1:length(mua_exp)
    
    curr_str = reshape(permute(mua_exp{curr_exp}(:,:,predict_roi),[2,1,3]),[],1);
    curr_str_baselinesub = reshape(permute(mua_exp{curr_exp}(:,:,predict_roi),[2,1,3]),[],1) - ...
        (nanmean(reshape(mua_exp{curr_exp}(:,t < 0,predict_roi),[],1)));
    curr_taskpred = reshape(permute(mua_taskpred_exp{curr_exp}(:,:,predict_roi),[2,1,3]),[],1);
    curr_ctxpred = reshape(permute(fluorpred_exp{curr_exp},[2,1,3]),[],1);
    
    % Set common NaNs
    nan_samples = isnan(curr_str) | isnan(curr_str_baselinesub) | ...
        isnan(curr_taskpred) | isnan(curr_ctxpred);
    curr_str(nan_samples) = NaN;
    curr_str_baselinesub(nan_samples) = NaN;
    curr_taskpred(nan_samples) = NaN;
    curr_ctxpred(nan_samples) = NaN;
    
    % (task regressed from average baseline-subtracted data)
    taskpred_r2(curr_exp,:) = 1 - (nansum((curr_str_baselinesub-curr_taskpred).^2,1)./ ...
        nansum((curr_str_baselinesub-nanmean(curr_str_baselinesub,1)).^2,1));
    % (cortex regressed from raw data)
    ctxpred_r2(curr_exp,:) = 1 - (nansum((curr_str-curr_ctxpred).^2,1)./ ...
        nansum((curr_str-nanmean(curr_str,1)).^2,1));
    
end

figure; plot(taskpred_r2,ctxpred_r2,'.k')
line([0,1],[0,1])
xlabel('Task-predicted')
ylabel('Cortex-predicted')


%% Load widefield-only task data

trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\paper_unused';
data_fn = 'trial_activity_choiceworld_wfonly'; % Widefield-only days (no craniotomy, so cleaner)
exclude_data = false;

AP_load_concat_normalize_ctx_str;

% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_recording;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));


%% Check peaky activity (using trial data)

U_downsample_factor = 10;
U_master_d = imresize(U_master(:,:,1:n_vs),1/U_downsample_factor,'bilinear');
framerate = 1/mean(diff(t));

prominence_map_all = cell(size(animals));
for curr_animal = 1:length(animals)
    
    disp(sprintf('Prominence: %s',animals{curr_animal}));
    curr_animal_n_days = max(trial_day(trial_animal == curr_animal));
    prominence_map = nan(size(U_master_d,1),size(U_master_d,2),curr_animal_n_days);
    
    for curr_day = 1:curr_animal_n_days
        
        use_trials = trial_animal == curr_animal & trial_day == curr_day;
        
        curr_v = reshape(permute(fluor_allcat_deconv(use_trials,:,:), ...
            [2,1,3]),[],n_vs)';
        curr_px = reshape(AP_svdFrameReconstruct(U_master_d,curr_v),[],size(curr_v,2));
        
        peak_thresh = nanstd(reshape(curr_px,[],1));
        
        curr_prominence_map = nan(size(U_master_d,1),size(U_master_d,2));
        for curr_px_idx = 1:size(curr_px,1)
            [~,~,widths,prominences] = findpeaks( ...
                curr_px(curr_px_idx,:), ...
                framerate,'MinPeakProminence',peak_thresh);
            curr_prominence_map(curr_px_idx) = nanstd(prominences./widths);
        end
        
        prominence_map(:,:,curr_day) = curr_prominence_map;
        AP_print_progress_fraction(curr_day,curr_animal_n_days);
    end
    
    prominence_map_all{curr_animal} = prominence_map;
    
end

figure;
c = [0,prctile(reshape(cat(3,prominence_map_all{:}),[],1),[100])*0.5];
h = tight_subplot(length(animals),1);
for curr_animal = 1:length(animals)   
    subplot(h(curr_animal));
    imagesc(reshape(prominence_map_all{curr_animal}, ...
        size(prominence_map_all{curr_animal},1),[]));
    axis image off;
    colormap(hot);  
    caxis(c);
    title(animals{curr_animal});
end


%% Save widefield movie for presentations

% set frames as movie_px
movie_px = b(:,:,t>-0.2 & t<0.7,1);

movie_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\presentations\220801_GallegoLab';
movie_fn = fullfile(movie_path,'task_novice_hemidiff.avi');
movie_framerate = 35/2;
movie_cmap = colormap(AP_colormap('WP',[],1.5));
movie_caxis = [0,0.003];
movie_position = [43,350,415,415];
AP_movie2avi(movie_px,movie_framerate,movie_cmap,movie_caxis,movie_position,movie_fn);




%% Paper fig test: rxn vs activity task performance index plot: difference, day x v x+1

plot_rxn_diff = diff(plot_rxn,[],1);
plot_act_diff = diff(plot_act,[],1);

figure; hold on;
set(gca,'ColorOrder',animal_col);
plot(plot_rxn_diff,plot_act_diff,'.','MarkerSize',15);
xlabel({'Reaction time idx','(alt-meas/alt+meas)'})
ylabel(wf_roi(plot_roi).area);
axis square tight;
xlim(xlim + 0.1*range(xlim).*[-1,1]);
ylim(ylim + 0.1*range(ylim).*[-1,1]);

for curr_animal = 1:length(animals)
    plot(plot_rxn_diff(learned_day(curr_animal)-1,curr_animal), ...
        plot_act_diff(learned_day(curr_animal)-1,curr_animal), ...
        'ok','MarkerSize',5,'linewidth',1);
end


figure; 

subplot(1,2,1);hold on;
set(gca,'ColorOrder',animal_col);
plot(plot_rxn(1:end-1,:),plot_rxn(2:end,:),'.','MarkerSize',15);
axis square tight;
xlim(xlim + 0.1*range(xlim).*[-1,1]);
ylim(ylim + 0.1*range(ylim).*[-1,1]);
for curr_animal = 1:length(animals)
    plot(plot_rxn(learned_day(curr_animal)-1,curr_animal), ...
        plot_rxn(learned_day(curr_animal),curr_animal), ...
        'ok','MarkerSize',5,'linewidth',1);
end
h = refline(1);
h.Color = 'k';
xlabel({'Reaction time index (given day)'})
ylabel({'Reaction time index (next day)'})


subplot(1,2,2);hold on;
set(gca,'ColorOrder',animal_col);
plot(plot_act(1:end-1,:),plot_act(2:end,:),'.','MarkerSize',15);
axis square tight;
xlim(xlim + 0.1*range(xlim).*[-1,1]);
ylim(ylim + 0.1*range(ylim).*[-1,1]);
for curr_animal = 1:length(animals)
    plot(plot_act(learned_day(curr_animal)-1,curr_animal), ...
        plot_act(learned_day(curr_animal),curr_animal), ...
        'ok','MarkerSize',5,'linewidth',1);
end
h = refline(1);
h.Color = 'k';
xlabel(sprintf('%s (given day)',wf_roi(plot_roi).area));
ylabel(sprintf('%s (next day)',wf_roi(plot_roi).area));




