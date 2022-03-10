%% Figures for longitudinal widefield during operant/sensorimotor learning
% Preprocessing done in AP_operant_learning_preprocessing

%% Behavior
% (adapted from AP_operant_behavior)

% Load behavior
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
bhv_fn = [data_path filesep 'bhv_teto'];
load(bhv_fn);
animals = {bhv.animal};

% Load muscimol injection info
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
muscimol_fn = [data_path filesep 'muscimol.mat'];
load(muscimol_fn);

% Use days before muscimol
use_days = cell(size(bhv));
for curr_animal = 1:length(bhv)
    muscimol_animal_idx = ismember({muscimol.animal},bhv(curr_animal).animal);
    if ~any(muscimol_animal_idx)
        use_days{curr_animal} = true(length(bhv(curr_animal).day),1);
        continue
    end
    muscimol_day_idx = datenum(bhv(curr_animal).day) >= ...
        datenum(muscimol(muscimol_animal_idx).day(1));
    use_days{curr_animal} = ~muscimol_day_idx;
end

% Set max days (for padding) and plot days (minimum n)
max_days = max(cellfun(@sum,use_days));
min_n = 4;
plot_days = find(accumarray(cell2mat(cellfun(@find,use_days,'uni',false)'),1) >= min_n);

% Set bins for reaction time histograms
rxn_bins = [0:0.02:0.5];
rxn_bin_centers = rxn_bins(1:end-1) + diff(rxn_bins)./2;

% Reaction time histograms by day
animal_rxn = nan(max_days,length(rxn_bin_centers),length(animals));
for curr_animal = 1:length(bhv)
    for curr_day = find(use_days{curr_animal})'
        animal_rxn(curr_day,:,curr_animal,:) = ...
            histcounts(bhv(curr_animal).stim_move_t{curr_day}, ...
            rxn_bins,'normalization','probability');
    end
end
figure;
imagesc(rxn_bin_centers,[],nanmean(animal_rxn(plot_days,:,:),3));
colormap(brewermap([],'Greys'));
h = colorbar;ylabel(h,'Probability');
xlabel('Reaction time');
ylabel('Day');

% Get reaction/resampled alt reaction times for all regular days
% (exclude trials without alts: rare, from wheel click issues)
use_rxn = cellfun(@(x) cellfun(@(x) cellfun(@(x) ...
    ~isempty(x),x),x,'uni',false),{bhv.alt_stim_move_t},'uni',false);

rxn_measured = cellfun(@(rxn,use_days,use_trials) ...
    cellfun(@(rxn,use_trials) rxn(use_trials),rxn(use_days),use_trials(use_days),'uni',false), ...
    {bhv.stim_move_t},use_days,use_rxn,'uni',false)';

n_rxn_altsample = 1000;
rxn_alt = cellfun(@(rxn,use_days,use_trials) ...
    cellfun(@(rxn,use_trials) ...
    cell2mat(cellfun(@(x) datasample(x,n_rxn_altsample)',rxn(use_trials),'uni',false)), ...
    rxn(use_days),use_trials(use_days),'uni',false), ...
    {bhv.alt_stim_move_t},use_days,use_rxn,'uni',false)';

% Concatenate all
rxn_measured_cat = cell2mat(cellfun(@cell2mat,rxn_measured,'uni',false));
rxn_alt_cat = cell2mat(cellfun(@cell2mat,rxn_alt,'uni',false));

% Index each trial for daysplit/day/animal and concat
n_daysplit = 4;
trial_split_idx = cellfun(@(x,animal_num) ...
    cellfun(@(x,day) ...
    [min(floor(linspace(1,n_daysplit+1,length(x))),n_daysplit)', ...
    repmat(day,length(x),1),repmat(animal_num,length(x),1)], ...
    x,num2cell(1:length(x))','uni',false), ...
    rxn_measured,num2cell(1:length(bhv))','uni',false);
trial_split_idx_cat = ...
    cell2mat(cellfun(@(x) cat(1,x{:}),trial_split_idx,'uni',false));

% Get histogram of reaction times > day 7 for each animal
animal_rxn_measured_cathist = cell2mat(arrayfun(@(x) ...
    histcounts(rxn_measured_cat( ...
    trial_split_idx_cat(:,3) == x & ...
    trial_split_idx_cat(:,2) >= 7), ...
    rxn_bins,'normalization','probability')',1:length(bhv),'uni',false));

animal_rxn_alt_cathist = cell2mat(permute( ...
    arrayfun(@(rep) cell2mat(arrayfun(@(x) ...
    histcounts(rxn_alt_cat( ...
    trial_split_idx_cat(:,3) == x & ...
    trial_split_idx_cat(:,2) >= 7,rep), ...
    rxn_bins,'normalization','probability')',1:length(bhv),'uni',false)), ...
    1:n_rxn_altsample,'uni',false),[1,3,2]));

figure; hold on

animal_rxn_alt_cathist_ci = ...
    squeeze(prctile(nanmean(animal_rxn_alt_cathist,2),[5,95],3));
AP_errorfill(rxn_bin_centers,nanmean(nanmean(animal_rxn_alt_cathist,2),3), ...
    animal_rxn_alt_cathist_ci,'r',[],false);

AP_errorfill(rxn_bin_centers,nanmean(animal_rxn_measured_cathist,2), ...
    AP_sem(animal_rxn_measured_cathist,2),'k');
xlabel('Reaction time');
ylabel('Probability');
legend({'Null','Measured'});

% Get fraction of reaction times within window
rxn_frac_window = [0.1,0.2];

rxn_measured_frac = accumarray(trial_split_idx_cat, ...
    rxn_measured_cat > rxn_frac_window(1) & ...
    rxn_measured_cat < rxn_frac_window(2)', ...
    [n_daysplit,max_days,length(bhv)],@nanmean,NaN);
rxn_alt_frac = cell2mat(permute(arrayfun(@(x) ...
    accumarray(trial_split_idx_cat, ...
    rxn_alt_cat(:,x) > rxn_frac_window(1) & ...
    rxn_alt_cat(:,x) < rxn_frac_window(2)', ...
    [n_daysplit,max_days,length(bhv)],@nanmean,NaN), ...
    1:n_rxn_altsample,'uni',false),[1,3,4,2]));

% Put NaNs between days to plot with gaps
rxn_measured_frac_long = reshape(padarray(rxn_measured_frac,[1,0,0],NaN,'post'),[],length(animals));
rxn_alt_frac_long = reshape(padarray(rxn_alt_frac,[1,0,0],NaN,'post'),[],length(animals),n_rxn_altsample);

% Plot relative to training day
figure; 
daysplit_x = [1:size(rxn_measured_frac_long,1)]/(n_daysplit+1);
subplot(1,2,1);hold on;

rxn_alt_frac_long_ci = ...
    permute(prctile(nanmean(rxn_alt_frac_long,2),[5,95],3),[1,3,2]);
AP_errorfill(daysplit_x,nanmean(nanmean(rxn_alt_frac_long,2),3), ...
    rxn_alt_frac_long_ci,'r',[],false);

errorbar(daysplit_x,nanmean(rxn_measured_frac_long,2), ...
    AP_sem(rxn_measured_frac_long,2),'k','linewidth',2,'CapSize',0);

xlabel('Training day');
ylabel(sprintf('Frac reaction times %.2f-%.2f',rxn_frac_window(1),rxn_frac_window(2)));
axis tight;
xlim(xlim + [-0.5,0.5]);
ylim([0,1]);

% Plot relative to learned day
learned_day = cellfun(@(x) find(x,1),{bhv.learned_days});
learned_day_x = [1:max_days]'-learned_day;

learned_daysplit_x = cell2mat(cellfun(@(x) x+(0:n_daysplit)'/(n_daysplit+1), ...
    num2cell(learned_day_x),'uni',false));

[rxn_learn_mean,rxn_learn_sem,learned_day_grp,learned_day_n] = ...
    grpstats(rxn_measured_frac_long(:),learned_daysplit_x(:), ...
    {'nanmean','sem','gname','numel'});
rxn_alt_learn_mean = ...
    grpstats(reshape(rxn_alt_frac_long,[],n_rxn_altsample),learned_daysplit_x(:), ...
    {'nanmean'});

learned_day_grp = cellfun(@str2num,learned_day_grp);
plot_learned = learned_day_n >= min_n | isnan(rxn_learn_mean);

subplot(1,2,2); hold on;

rxn_alt_learn_ci = prctile(rxn_alt_learn_mean,[5,95],2);
AP_errorfill(learned_day_grp(plot_learned), ...
    nanmean(rxn_alt_learn_mean(plot_learned,:),2), ...
    rxn_alt_learn_ci(plot_learned,:),'r',[],false);

errorbar(learned_day_grp(plot_learned),rxn_learn_mean(plot_learned), ...
    rxn_learn_sem(plot_learned),'k','linewidth',2,'CapSize',0);
xlabel('Learned day');
axis tight;
xlim(xlim + [-0.5,0.5]);
ylim([0,1]);
line([0,0],ylim,'color','k','linestyle','--');


%% Behavior: muscimol

% Set animal to use
% (only use tetO mice with successful injections evidenced by retinotopic
% map loss - manually set)
muscimol_v1_animals = {'AP100','AP105','AP106','AP107','AP108'};

% Load behavior
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
bhv_fn = [data_path filesep 'bhv_teto'];
bhv_all = load(bhv_fn);
animals_all = {bhv_all.bhv.animal};
use_animals = ismember(animals_all,muscimol_v1_animals);
bhv = bhv_all.bhv(use_animals);

% Load muscimol injection info
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
muscimol_fn = [data_path filesep 'muscimol.mat'];
load(muscimol_fn);

% Grab day index of [V1 muscimol, washout] and data to plot
muscimol_v1_days = nan(length(bhv),2);
muscimol_v1_retinotopy = cell(size(muscimol_v1_days));
muscimol_v1_stim_surround_wheel = cell(size(muscimol_v1_days));
muscimol_v1_wheel_mm = nan(size(muscimol_v1_days));
for curr_animal = 1:length(bhv)

    muscimol_animal_idx = ismember({muscimol.animal},bhv(curr_animal).animal);
    if ~any(muscimol_animal_idx)
        continue
    end
    % (find last V1 muscimol day)
    curr_v1_muscimol = find(strcmp(lower( ...
        muscimol(muscimol_animal_idx).area),'v1'),1,'last');

    % (find first washout after V1 muscimol)
     curr_v1_washout = curr_v1_muscimol + ...
        find(strcmp(lower( ...
        muscimol(muscimol_animal_idx).area(curr_v1_muscimol+1:end)), ...
        'washout'),1,'first');

    curr_v1_muscimol_dayidx = find(strcmp(bhv(curr_animal).day, ...
        muscimol(muscimol_animal_idx).day(curr_v1_muscimol)));

    curr_v1_washout_dayidx = find(strcmp(bhv(curr_animal).day, ...
        muscimol(muscimol_animal_idx).day(curr_v1_washout)));

    % Grab days
     muscimol_v1_days(curr_animal,:) = ...
             [curr_v1_muscimol_dayidx,curr_v1_washout_dayidx];

     % Grab retinotopy
    muscimol_v1_retinotopy(curr_animal,:) = ...
        muscimol(muscimol_animal_idx).vfs([curr_v1_muscimol,curr_v1_washout]);

    % Grab stim-aligned movement
    % (NaN-out times during quiescence)
    stim_surround_t = bhv(1).stim_surround_t;
    curr_stim_surround_wheel = ...
        bhv(curr_animal).stim_surround_wheel([curr_v1_muscimol_dayidx,curr_v1_washout_dayidx]);
    curr_quiescence_t = ...
        bhv(curr_animal).quiescence_t([curr_v1_muscimol_dayidx,curr_v1_washout_dayidx]);
    for curr_cond = 1:2
        for curr_trial = 1:size(curr_stim_surround_wheel{curr_cond},1)
            q_time = stim_surround_t >= -curr_quiescence_t{curr_cond}(curr_trial) & ...
                stim_surround_t <= 0;
            curr_stim_surround_wheel{curr_cond}(curr_trial,q_time) = NaN;
        end
    end

    muscimol_v1_stim_surround_wheel(curr_animal,:) = curr_stim_surround_wheel;

    % Grab wheel travel/time
    muscimol_v1_wheel_mm(curr_animal,:) = ...
        bhv(curr_animal).wheel_mm([curr_v1_muscimol_dayidx,curr_v1_washout_dayidx])./ ...
        bhv(curr_animal).session_duration([curr_v1_muscimol_dayidx,curr_v1_washout_dayidx]);

end

% Plot retinotopy difference
figure;
for curr_cond = 1:2
    subplot(1,2,curr_cond);
    imagesc(nanmean(cat(3,muscimol_v1_retinotopy{:,curr_cond}),3));
    axis image off;
    caxis([-1,1]);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    switch curr_cond
        case 1
            title('V1 muscimol');
        case 2
            title('Washout');
    end
end
colormap(brewermap([],'*RdBu'));

% (set 'use_days' to copy code from above)
use_days = mat2cell(muscimol_v1_days,ones(length(muscimol_v1_animals),1),2)';

% Get reaction/resampled alt reaction times for all regular days
% (exclude trials without alts: rare, from wheel click issues)
use_rxn = cellfun(@(x) cellfun(@(x) cellfun(@(x) ...
    ~isempty(x),x),x,'uni',false),{bhv.alt_stim_move_t},'uni',false);

rxn_measured = cellfun(@(rxn,use_days,use_trials) ...
    cellfun(@(rxn,use_trials) rxn(use_trials),rxn(use_days),use_trials(use_days),'uni',false), ...
    {bhv.stim_move_t},use_days,use_rxn,'uni',false)';

n_rxn_altsample = 1000;
rxn_alt = cellfun(@(rxn,use_days,use_trials) ...
    cellfun(@(rxn,use_trials) ...
    cell2mat(cellfun(@(x) datasample(x,n_rxn_altsample)',rxn(use_trials),'uni',false)), ...
    rxn(use_days),use_trials(use_days),'uni',false), ...
    {bhv.alt_stim_move_t},use_days,use_rxn,'uni',false)';

% Concat muscimol/washout
rxn_measured_cat = cat(2,rxn_measured{:});
rxn_alt_cat = cat(2,rxn_alt{:});

% Set bins for reaction time histograms
rxn_bins = [0:0.02:0.5];
rxn_bin_centers = rxn_bins(1:end-1) + diff(rxn_bins)./2;

% Plot reaction time histogram for pre/post muscimol
rxn_measured_hist = cellfun(@(x) ...
    histcounts(x,rxn_bins,'normalization','probability'), ...
    rxn_measured_cat,'uni',false);

rxn_alt_hist = cellfun(@(x) cell2mat(arrayfun(@(rep)...
    histcounts(x(:,rep),rxn_bins,'normalization','probability'), ...
    1:n_rxn_altsample,'uni',false)'),rxn_alt_cat,'uni',false);

rxn_alt_hist_ci = ...
    arrayfun(@(cond) prctile(nanmean(cat(3,rxn_alt_hist{cond,:}),3), ...
    [5,95],1),1:2,'uni',false);

figure;
h = tiledlayout(2,2,'TileSpacing','compact','padding','compact');
for curr_cond = 1:2
    nexttile; hold on
    AP_errorfill(rxn_bin_centers,nanmean(nanmean(cat(3,rxn_alt_hist{1,:}),3),1), ...
        rxn_alt_hist_ci{curr_cond}','r',[],false);

    AP_errorfill(rxn_bin_centers, ...
        nanmean(cat(1,rxn_measured_hist{curr_cond,:}),1)', ...
        AP_sem(cat(1,rxn_measured_hist{curr_cond,:}),1)','k');

    legend({'Null','Measured'});
    xlabel('Reaction time');
    ylabel('Probability')
    switch curr_cond
        case 1
            title('V1 muscimol');
        case 2
            title('Washout');
    end
    
end
linkaxes(allchild(h),'xy');

% Plot stim-aligned wheel movement
muscimol_stim_wheel_avg = ...
    arrayfun(@(cond) cell2mat(cellfun(@(x) ...
    nanmean(x,1),muscimol_v1_stim_surround_wheel(:,cond),'uni',false)), ...
    1:2,'uni',false);
nexttile;
AP_errorfill(bhv(1).stim_surround_t', ...
    cell2mat(cellfun(@(x) nanmean(x,1)',muscimol_stim_wheel_avg,'uni',false)), ...
    cell2mat(cellfun(@(x) AP_sem(x,1)',muscimol_stim_wheel_avg,'uni',false)), ...
    [0,0,0;1,0,0]);
xline(0,'linestyle','--');
xlabel('Time from stim (s)');
ylabel('Probability of movement');

% Plot total velocity
nexttile;
errorbar(nanmean(muscimol_v1_wheel_mm,1), ...
    AP_sem(muscimol_v1_wheel_mm,1),'k','linewidth',2);
xlim(xlim+[-0.5,0.5]);
ylim([0,max(ylim)])
set(gca,'XTick',1:2,'XTickLabel',{'V1 muscimol','Washout'});
ylabel('Wheel mm/min');

%% Whisker movement (passive)

% Load facecam align
facecam_processing_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\facecam_processing';
facecam_align_fn = fullfile(facecam_processing_path,'facecam_align.mat');
load(facecam_align_fn);

% Align facecams
% (reference is first image)
im_unaligned_cat = cellfun(@double,[facecam_align.im],'uni',false);
im_tform_cat = cat(1,facecam_align.tform);
im_ref = im_unaligned_cat{1};

im_aligned = nan(size(im_ref,1),size(im_ref,2),length(im_unaligned_cat));
for curr_im =1:length(im_unaligned_cat)
    if isempty(im_unaligned_cat{curr_im})
        continue
    end
    tform_size = imref2d(size(im_ref));
    im_aligned(:,:,curr_im) = ...
        imwarp(im_unaligned_cat{curr_im},im_tform_cat{curr_im}, ...
        'OutputView',tform_size);
end

% Concat whisker masks (pad whisker masks to same size as im)
whisker_mask_pad = cellfun(@(x,y) [x,cell(1,length(y)-length(x))], ...
    {facecam_align.whisker_mask},{facecam_align.im},'uni',false);
whisker_mask_cat = horzcat(whisker_mask_pad{:});

% Plot sample frame/mask, average frame/mask
use_days = cellfun(@(x,y) ~isempty(x) & ~isempty(y), ...
    im_unaligned_cat,whisker_mask_cat);

figure;
subplot(1,2,1);
ref_plot = imoverlay(mat2gray(im_ref),whisker_mask_cat{1},'r');
imagesc(ref_plot);
axis image off;
title('Sample image and whisker ROI');
subplot(1,2,2)
avg_plot = imoverlay(mat2gray(nanmean(im_aligned(:,:,use_days),3)), ...
    whisker_mask_cat{1},'r');
imagesc(avg_plot);
axis image off;
title('Average image and whisker ROI')



%% Passive - widefield
% (fix the indexing in this to be like in task, don't average first)

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_passive_tetO';
AP_load_trials_operant;
n_naive = 3; % (number of naive passive-only days)
min_n = 4; % (minimum n to plot data)

% Load behavior and get learned day
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
bhv_fn = [data_path filesep 'bhv_teto'];
load(bhv_fn);
learned_day = cellfun(@(x) find(x,1),{bhv.learned_days})';

% Get animal and day index for each trial
trial_animal = cell2mat(arrayfun(@(x) ...
    x*ones(size(vertcat(wheel_all{x}{:}),1),1), ...
    [1:length(wheel_all)]','uni',false));
trial_day = cell2mat(cellfun(@(x) cell2mat(cellfun(@(curr_day,x) ...
    curr_day*ones(size(x,1),1),num2cell(1:length(x))',x,'uni',false)), ...
    wheel_all,'uni',false));

trial_learned_day = cell2mat(cellfun(@(x,ld) cell2mat(cellfun(@(curr_day,x) ...
    curr_day*ones(size(x,1),1),num2cell([1:length(x)]-ld)',x,'uni',false)), ...
    wheel_all,num2cell(learned_day + n_naive),'uni',false));

learned_day_unique = unique(trial_learned_day)';
[~,trial_learned_day_id] = ismember(trial_learned_day,learned_day_unique);

trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));

% Get trials with movement during stim to exclude
quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > 0,2);

% Get average fluorescence by animal, day, stim
stim_unique = unique(trial_stim_allcat);
stim_v_avg = cell(length(animals),1);
stim_roi_avg = cell(length(animals),1);
for curr_animal = 1:length(animals)
    for curr_day = 1:max(trial_day(trial_animal == curr_animal))
        for curr_stim_idx = 1:length(stim_unique)
            use_trials = quiescent_trials & ...
                trial_animal == curr_animal & ...
                trial_day == curr_day & ...
                trial_stim_allcat == stim_unique(curr_stim_idx);
            stim_v_avg{curr_animal}(:,:,curr_day,curr_stim_idx) = ...
                permute(nanmean(fluor_allcat_deconv(use_trials,:,:),1),[3,2,1]);
            stim_roi_avg{curr_animal}(:,:,curr_day,curr_stim_idx) = ...
                permute(nanmean(fluor_roi_deconv(use_trials,:,:),1),[3,2,1]);
        end
    end
end

% Average V's by naive/prelearn/postlearn, plot time average
stim_v_avg_stage = cell2mat(permute(cellfun(@(x,ld) ...
    cat(3, ...
    nanmean(x(:,:,1:n_naive,:),3), ...
    nanmean(x(:,:,n_naive+1:ld-1,:),3), ...
    nanmean(x(:,:,ld:end,:),3)),stim_v_avg,num2cell(learned_day), ...
    'uni',false),[2,3,4,5,1]));

stim_px_avg_stage_movie = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),squeeze(nanmean(stim_v_avg_stage(:,:,:,3,:),5)));
AP_image_scroll(reshape(permute(stim_px_avg_stage_movie,[1,2,4,3]),size(U_master,1),[],length(t)),t);
axis image off;
colormap(brewermap([],'PrGn'));
caxis([-max(abs(caxis)),max(abs(caxis))]);

use_t = t >= 0 & t <= 0.2;
stim_px_avg_stage_tavg = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    squeeze(nanmean(nanmean(stim_v_avg_stage(:,use_t,:,:,:),2),5)));

figure;
tiledlayout(3,3,'TileSpacing','compact','padding','compact');
c = (max(stim_px_avg_stage_tavg(:)).*[-1,1])*0.5;
for curr_stage = 1:3
    for curr_stim = 1:3
        nexttile;
        imagesc(stim_px_avg_stage_tavg(:,:,curr_stage,curr_stim));
        axis image off;
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        colormap(brewermap([],'PrGn'));
        caxis(c);
        title(sprintf('Stage %d, stim %d',curr_stage,curr_stim));
    end
end

% Average ROIs by naive/prelearn/postlearn
stim_roi_avg_stage = cell2mat(permute(cellfun(@(x,ld) ...
    cat(3, ...
    nanmean(x(:,:,1:n_naive,:),3), ...
    nanmean(x(:,:,n_naive+1:ld-1,:),3), ...
    nanmean(x(:,:,ld:end,:),3)),stim_roi_avg,num2cell(learned_day), ...
    'uni',false),[2,3,4,5,1]));

figure;
plot_rois = [1,5,6];
stage_col = [0,0,0;0.7,0,0;0,0.7,0];
tiledlayout(length(plot_rois),2,'TileSpacing','compact','padding','compact');
for curr_roi_idx = 1:length(plot_rois)
    curr_l_roi = plot_rois(curr_roi_idx);
    curr_r_roi = curr_l_roi + size(wf_roi,1);

    % (plot left ROI w/ right stim)
    hl = nexttile;
    AP_errorfill(t, ...
        squeeze(nanmean(stim_roi_avg_stage(curr_l_roi,:,:,stim_unique == 1,:),5)), ...
        squeeze(AP_sem(stim_roi_avg_stage(curr_l_roi,:,:,stim_unique == 1,:),5)),stage_col);
    xlabel('Time from stim (s)');
    ylabel('\DeltaF/F_0');
    title(wf_roi(curr_l_roi).area);
    xline(0,'linestyle','--');xline(0.5,'linestyle','--');

    % (plot right ROI with left stim)
    hr =nexttile;
    AP_errorfill(t, ...
        squeeze(nanmean(stim_roi_avg_stage(curr_r_roi,:,:,stim_unique == -1,:),5)), ...
        squeeze(AP_sem(stim_roi_avg_stage(curr_r_roi,:,:,stim_unique == -1,:),5)),stage_col);
    xlabel('Time from stim (s)');
    ylabel('\DeltaF/F_0');
    title(wf_roi(curr_r_roi).area);
    xline(0,'linestyle','--');xline(0.5,'linestyle','--');

    % (link axes with same ROI)
    linkaxes([hl,hr],'xy');
end

% Plot all trials (by stage)
plot_rois = [1,6];
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);
n_trial_smooth = 20;

figure;
tiledlayout(2,length(plot_rois),'TileSpacing','compact','padding','compact');
for curr_roi = plot_rois
    for curr_stage = 1:max(trial_learned_stage)
        use_trials = trial_stim_allcat == 1 & ...
            trial_learned_stage == curr_stage & quiescent_trials;
        
        curr_data = fluor_roi_deconv(use_trials,:,curr_roi);
        curr_data_smooth = convn(curr_data, ...
            ones(n_trial_smooth,1)./n_trial_smooth,'same');
        
        nexttile;
        imagesc(t,[],curr_data_smooth);hold on;
        colormap(brewermap([],'PrGn'));
        caxis(max(caxis)*[-1,1])
        xline(0,'color','r');
        xlabel('Time from stim (s)');
        ylabel('Trial');
        title(sprintf('%s, stage %d',wf_roi(curr_roi).area,curr_stage));
    end
end

% Plot time-max within ROIs across days
use_t = t >= 0 & t <= 0.5;
stim_roi_tmax = cellfun(@(x) ...
    permute(max(x(:,use_t,:,:),[],2),[1,3,4,2]), ...
    stim_roi_avg,'uni',false);
stim_roi_tmax_daycat = cat(2,stim_roi_tmax{:});

n_days_animal = accumarray(trial_animal,trial_day,[],@max);
trained_day_animal_x = cellfun(@(n) [1:n]',num2cell(n_days_animal),'uni',false);
learned_day_animal_x = cellfun(@(ld,n) [1:n]'-(ld+n_naive), ...
    num2cell(learned_day),num2cell(n_days_animal),'uni',false);
[~,learned_day_animal_id] = cellfun(@(x) ismember(x,learned_day_unique), ...
    learned_day_animal_x,'uni',false);

% (grab indicies)
[roi_idx,learned_day_idx,stim_idx] = ...
    ndgrid(1:n_rois,cat(1,learned_day_animal_id{:}),1:3);
[~,trained_day_idx,~] = ...
    ndgrid(1:n_rois,cat(1,trained_day_animal_x{:}),1:3);
training_data = trained_day_idx > n_naive; % (exclude naive passive-only)

% (get mean/sem by index during naive)
stim_roi_tmax_naive_mean = ...
    accumarray([roi_idx(~training_data),trained_day_idx(~training_data),stim_idx(~training_data)], ...
    stim_roi_tmax_daycat(~training_data), ...
    [n_rois,max(trained_day_idx(:)),length(stim_unique)],@nanmean,NaN('single'));
stim_roi_tmax_naive_sem = ...
    accumarray([roi_idx(~training_data),trained_day_idx(~training_data),stim_idx(~training_data)], ...
    stim_roi_tmax_daycat(~training_data), ...
    [n_rois,max(trained_day_idx(:)),length(stim_unique)],@AP_sem,NaN('single'));

% (get mean/sem by index during training)
stim_roi_tmax_learn_mean = ...
    accumarray([roi_idx(training_data),learned_day_idx(training_data),stim_idx(training_data)], ...
    stim_roi_tmax_daycat(training_data), ...
    [n_rois,max(learned_day_idx(:)),length(stim_unique)],@nanmean,NaN('single'));
stim_roi_tmax_learn_sem = ...
    accumarray([roi_idx(training_data),learned_day_idx(training_data),stim_idx(training_data)], ...
    stim_roi_tmax_daycat(training_data), ...
    [n_rois,max(learned_day_idx(:)),length(stim_unique)],@AP_sem,NaN('single'));
% (get number of elements)
stim_roi_tmax_learn_numel = ...
    accumarray([roi_idx(training_data),learned_day_idx(training_data),stim_idx(training_data)], ...
    stim_roi_tmax_daycat(training_data), ...
    [n_rois,max(learned_day_idx(:)),length(stim_unique)],@numel,NaN);
% (plot days with minimum n)
plot_learned_day = squeeze(stim_roi_tmax_learn_numel(1,:,1)) > min_n;

figure;
plot_rois = [1,5,6];
tiledlayout(length(plot_rois),2,'TileSpacing','compact','padding','compact');
for curr_roi_idx = 1:length(plot_rois)
    curr_l_roi = plot_rois(curr_roi_idx);
    curr_r_roi = curr_l_roi + size(wf_roi,1);

    stim_col = [0,0,1;0.5,0.5,0.5;1,0,0];

    % Plot L ROI
    hl = nexttile; hold on;
    % (naive)
    AP_errorfill(min(learned_day_unique(plot_learned_day))+[-n_naive:-1], ...
        squeeze(stim_roi_tmax_naive_mean(curr_l_roi,1:n_naive,:)), ...
        squeeze(stim_roi_tmax_naive_sem(curr_l_roi,1:n_naive,:)),stim_col);
    % (training)
    AP_errorfill(learned_day_unique(plot_learned_day), ...
        squeeze(stim_roi_tmax_learn_mean(curr_l_roi,plot_learned_day,:)), ...
        squeeze(stim_roi_tmax_learn_sem(curr_l_roi,plot_learned_day,:)),stim_col);

    title(wf_roi(curr_l_roi).area);
    xlabel('Learned day');
    ylabel('\DeltaF/F_0');
    xline(0,'linestyle','--');
    axis tight;
    xlim(xlim+[-0.5,0.5]);

    % Plot R ROI
    hr = nexttile;
    % (naive)
    AP_errorfill(min(learned_day_unique(plot_learned_day))+[-n_naive:-1], ...
        squeeze(stim_roi_tmax_naive_mean(curr_r_roi,1:n_naive,:)), ...
        squeeze(stim_roi_tmax_naive_sem(curr_r_roi,1:n_naive,:)),stim_col);
    % (training)
    AP_errorfill(learned_day_unique(plot_learned_day), ...
        squeeze(stim_roi_tmax_learn_mean(curr_r_roi,plot_learned_day,:)), ...
        squeeze(stim_roi_tmax_learn_sem(curr_r_roi,plot_learned_day,:)),stim_col);

    title(wf_roi(curr_r_roi).area);
    xlabel('Learned day');
    ylabel('\DeltaF/F_0');
    xline(0,'linestyle','--');
    axis tight;
    xlim(xlim+[-0.5,0.5]);

    % (link axes with same ROI)
    linkaxes([hl,hr],'xy');

end

% Get average pupil diameter and whisker movement to stimuli
[~,trial_stim_allcat_id] = ismember(trial_stim_allcat,stim_unique);

learned_day_unique = unique(trial_learned_day)';
[~,trial_learned_day_id] = ismember(trial_learned_day,learned_day_unique);

[stim_idx,t_idx] = ndgrid(trial_stim_allcat_id,1:length(t));
[animal_idx,~] = ndgrid(trial_animal,1:length(t));
[learned_day_id_idx,~] = ndgrid(trial_learned_day_id,1:length(t));

bhv_accum_idx = cat(3,stim_idx(quiescent_trials,:), ...
    t_idx(quiescent_trials,:), ...
    learned_day_id_idx(quiescent_trials,:), ... 
    animal_idx(quiescent_trials,:));

pupil_diff_stim_avg = accumarray(reshape(bhv_accum_idx,[],size(bhv_accum_idx,3)), ...
    reshape(pupil_diameter_allcat_diff(quiescent_trials,:),[],1),[length(stim_unique),length(t), ...
    max(trial_learned_day_id),length(animals)],@nanmean,NaN);

whisker_stim_avg = accumarray(reshape(bhv_accum_idx,[],size(bhv_accum_idx,3)), ...
    reshape(whisker_allcat(quiescent_trials,:),[],1),[length(stim_unique),length(t), ...
    max(trial_learned_day_id),length(animals)],@nanmean,NaN);

bhv_stim_avg = cat(5,pupil_diff_stim_avg,whisker_stim_avg);
bhv_label = {'Pupil diff','Whisker'};

[~,~,learning_stage] = unique(learned_day_unique >= 0);

figure;
h = tiledlayout(size(bhv_stim_avg,5),max(learning_stage), ...
    'TileSpacing','compact','padding','compact');
stim_col = [0,0,1;0,0,0;1,0,0];
for curr_bhv = 1:size(bhv_stim_avg,5)
    for curr_stage = 1:max(learning_stage)
        nexttile;
        curr_mean = nanmean(nanmean(bhv_stim_avg(:,:, ...
            learning_stage == curr_stage,:,curr_bhv),3),4)';
        curr_sem = AP_sem(nanmean(bhv_stim_avg(:,:, ...
            learning_stage == curr_stage,:,curr_bhv),3),4)';
        AP_errorfill(t,curr_mean,curr_sem,stim_col);
        xlabel('Time from stim (s)');
        ylabel(bhv_label{curr_bhv});
        title(sprintf('Stage %d',curr_stage));
        axis tight;
        xlim(xlim+[-0.1,0.1]);
        xline(0,'linestyle','--');xline(0.5,'linestyle','--');
    end
    curr_axes = allchild(h);
    linkaxes(curr_axes(1:max(learning_stage)),'xy');
end

% Plot whisker movement over time
use_t = t > 0 & t < 0.2;
whisker_stim_tmax = squeeze(max(whisker_stim_avg(:,use_t,:,:),[],2));
plot_days = min(sum(~isnan(whisker_stim_tmax),3)) > min_n;
figure; hold on
AP_errorfill(learned_day_unique(plot_days), ...
    nanmean(whisker_stim_tmax(:,plot_days,:),3)', ...
    AP_sem(whisker_stim_tmax(:,plot_days,:),3)',stim_col);
xline(0,'color','k','linestyle','--');
xlabel('Learned day');
ylabel('Whisker movement')

% Plot fluorescence/whisker by discretized whisker movement
plot_rois = [1,6];
min_trials = 10; % (minimum trials to use within recording)

trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);

% (discretize whisker movement within recording - use R stim trials)
use_t = t > 0 & t < 0.2;
whisker_allcat_tavg = nanmean(whisker_allcat(:,use_t),2);
whisker_grp_use_trials = quiescent_trials & ...
        ~all(isnan(whisker_allcat),2);
    
whisker_recording = mat2cell(whisker_allcat_tavg,trials_recording);
whisker_grp_use_trials_recording = mat2cell(whisker_grp_use_trials,trials_recording);
whisker_grp_use_recordings = cellfun(@(whisker,use_trials) ...
    ~all(isnan(whisker)) & sum(use_trials) > min_trials,...
    whisker_recording,whisker_grp_use_trials_recording);

n_whisker_grps = 4;
whisker_grp = cellfun(@(x) nan(size(x)),whisker_recording,'uni',false);
whisker_grp(whisker_grp_use_recordings) = cellfun(@(whisker,use_trials) ...
    discretize(whisker, ...
    prctile(whisker(use_trials),linspace(0,100,n_whisker_grps+1))), ...
    whisker_recording(whisker_grp_use_recordings), ...
    whisker_grp_use_trials_recording(whisker_grp_use_recordings),'uni',false);
whisker_grp = cell2mat(whisker_grp);

% (get indicies for grouping)
[animal_idx,t_idx,roi_idx] = ndgrid(trial_animal,1:length(t),1:n_rois);
[whisker_grp_idx,~,~] = ndgrid(whisker_grp,1:length(t),1:n_rois);
[learned_stage_idx,~,~] = ndgrid(trial_learned_stage,1:length(t),1:n_rois);
accum_whisker_idx = cat(4,whisker_grp_idx,t_idx,learned_stage_idx,roi_idx,animal_idx);

% (plot grouped whisker movements and ROI fluorescence)
line_fig = figure;
stim_col = [0,0,1;0,0,0;1,0,0];
whisker_grp_col = copper(n_whisker_grps);
ax = gobjects((length(plot_rois)+1)*max(trial_learned_stage),length(stim_unique));
for curr_stim_idx = 1:length(stim_unique)
    
    curr_stim = stim_unique(curr_stim_idx);
    use_trials = trial_stim_allcat == curr_stim & quiescent_trials & ...
        ~all(isnan(whisker_allcat),2) & ~isnan(whisker_grp);

    whisker_whiskergrp_avg = accumarray( ...
        reshape(accum_whisker_idx(use_trials,:,1,[1,2,3,5],:),[],size(accum_whisker_idx,4)-1), ...
        reshape(whisker_allcat(use_trials,:),[],1),[],@nanmean,NaN);
    
    fluor_roi_whiskergrp_avg = accumarray( ...
        reshape(accum_whisker_idx(use_trials,:,:,:),[],size(accum_whisker_idx,4)), ...
        reshape(fluor_roi_deconv(use_trials,:,:),[],1),[],@nanmean,NaN('single'));
    
    figure;
    h = tiledlayout(max(trial_learned_stage),length(plot_rois)+1, ...
        'TileSpacing','compact','padding','compact');
    for curr_stage = 1:max(trial_learned_stage)
        nexttile;
        AP_errorfill(t, ...
            nanmean(whisker_whiskergrp_avg(:,:,curr_stage,:),4)', ...
            AP_sem(whisker_whiskergrp_avg(:,:,curr_stage,:),4)', ...
            whisker_grp_col);
        title(sprintf('Whiskers, stage %d',curr_stage));
        axis tight;
        xlim(xlim+[-0.1,0.1]);
        xline(0,'--k');
        xlabel('Time from stim (s)');
        ylabel('Pixel intensity change');
        
        for curr_roi = plot_rois
            nexttile;
            AP_errorfill(t, ...
                squeeze(nanmean(fluor_roi_whiskergrp_avg(:,:,curr_stage,curr_roi,:),5))', ...
                squeeze(AP_sem(fluor_roi_whiskergrp_avg(:,:,curr_stage,curr_roi,:),5))', ...
                whisker_grp_col);
            title(sprintf('%s, stage %d',wf_roi(curr_roi).area,curr_stage));
            axis tight;
            xlim(xlim+[-0.1,0.1]);
            xline(0,'--k');
            xlabel('Time from stim (s)');
            ylabel('\DeltaF/F0');
        end
    end
    ax(:,curr_stim_idx) = allchild(h);
    title(h,sprintf('Whisker move groups (stim %d)',curr_stim));
    drawnow;
    
    % Plot line plot of fluorescence vs whisker
    stage_linestyle = {'--','-'};
    for curr_stage = 1:max(trial_learned_stage)
        for curr_roi_idx = 1:length(plot_rois)
            curr_roi = plot_rois(curr_roi_idx);
            curr_ax = subplot(1,length(plot_rois),curr_roi_idx,'Parent',line_fig);
            hold(curr_ax,'on');
            curr_whisker_whiskergrp_tavg = ...
                squeeze(nanmean(whisker_whiskergrp_avg(:,use_t,curr_stage,:),2));
            curr_fluor_roi_whiskergrp_tavg = ...
                squeeze(nanmean(fluor_roi_whiskergrp_avg(:,use_t,curr_stage,curr_roi,:),2));
            
            x = nanmean(curr_whisker_whiskergrp_tavg,2);
            x_err = AP_sem(curr_whisker_whiskergrp_tavg,2);
            y = nanmean(curr_fluor_roi_whiskergrp_tavg,2);
            y_err = AP_sem(curr_fluor_roi_whiskergrp_tavg,2);
            errorbar(curr_ax,x,y,-y_err,y_err,-x_err,x_err, ...
                'linewidth',2,'color',stim_col(curr_stim_idx,:), ...
                'linestyle',stage_linestyle{curr_stage});
            
            xlabel(curr_ax,'Whisker movement');
            ylabel(curr_ax,'\DeltaF/F_0');
            title(curr_ax,wf_roi(curr_roi).area);
        end
    end
end
% (link modality axes across figures)
for i = 1:length(plot_rois)+1
    linkaxes(ax(i:length(plot_rois)+1:end,:),'xy');
end
% (legend for line plots - yes, ridiculous code)
line_fig_ax = get(line_fig,'Children');
line_fig_ax_lines = get(line_fig_ax(1),'Children');
legend(line_fig_ax_lines, ...
    fliplr(flipud(cellfun(@(x,y) [x ', ' y], ...
    repmat(arrayfun(@(x) sprintf('Stage %d',x), ...
    1:max(trial_learned_stage),'uni',false),length(stim_unique),1), ...
    repmat(arrayfun(@(x) sprintf('Stim %d',x), ...
    stim_unique,'uni',false),1,max(trial_learned_stage)),'uni',false)')), ...
    'location','nw');


%% Task - widefield

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_task_teto';
AP_load_trials_operant;
min_n = 4; % (minimum n to plot data)

% Load behavior and get learned day
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
bhv_fn = [data_path filesep 'bhv_teto'];
load(bhv_fn);
learned_day = cellfun(@(x) find(x,1),{bhv.learned_days})';

% Get animal and day index for each trial
trial_animal = cell2mat(arrayfun(@(x) ...
    x*ones(size(vertcat(wheel_all{x}{:}),1),1), ...
    [1:length(wheel_all)]','uni',false));
trial_day = cell2mat(cellfun(@(x) cell2mat(cellfun(@(curr_day,x) ...
    curr_day*ones(size(x,1),1),num2cell(1:length(x))',x,'uni',false)), ...
    wheel_all,'uni',false));

trial_learned_day = trial_day - learned_day(trial_animal);

learned_day_unique = unique(trial_learned_day)';
[~,trial_learned_day_id] = ismember(trial_learned_day,learned_day_unique);

% Get number of trials by animal/recording
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));

% Plot average full + reduced stimulus pixel activity pre/post learning
stim_regressor = strcmp(task_regressor_labels,'Stim');
stim_v_act = fluor_allcat_deconv - fluor_taskpred_reduced_allcat(:,:,:,stim_regressor);

[learned_idx,t_idx,v_idx] = ...
    ndgrid((trial_learned_day >= 0)+1,1:length(t),1:n_vs);
[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_vs);

n_stages = max(learned_idx(:));
full_v_act_learn_avg = permute(accumarray( ...
    [learned_idx(:),t_idx(:),v_idx(:),animal_idx(:)], ...
    fluor_allcat_deconv(:), ...
    [n_stages,length(t),n_vs,length(animals)], ...
    @nanmean,NaN('single')),[3,2,1,4]);

stim_v_act_learn_avg = permute(accumarray( ...
    [learned_idx(:),t_idx(:),v_idx(:),animal_idx(:)], ...
    stim_v_act(:), ...
    [n_stages,length(t),n_vs,length(animals)], ...
    @nanmean,NaN('single')),[3,2,1,4]);

full_px_act_learn_avg = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    squeeze(nanmean(full_v_act_learn_avg,4)));

stim_px_act_learn_avg = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    squeeze(nanmean(stim_v_act_learn_avg,4)));

figure;
use_t = t >= 0 & t <= 0.2;
tiledlayout(n_stages,2,'TileSpacing','compact','padding','compact');
c_full = (max(full_px_act_learn_avg(:)).*[-1,1])*0.7;
c_stim = (max(stim_px_act_learn_avg(:)).*[-1,1])*0.7;
for curr_stage = 1:n_stages
    nexttile;
    imagesc(nanmean(full_px_act_learn_avg(:,:,use_t,curr_stage),3));
    caxis(c_full);
    axis image off;
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    colormap(brewermap([],'PrGn'));
    title(sprintf('Full, stage %d',curr_stage));
    
    nexttile;
    imagesc(nanmean(stim_px_act_learn_avg(:,:,use_t,curr_stage),3));
    caxis(c_stim);
    axis image off;
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    colormap(brewermap([],'PrGn'));
    title(sprintf('Stim, stage %d',curr_stage));
end

% Get ROI average/time-max raw/reduced activity by learned stage/daysplit
stim_regressor = strcmp(task_regressor_labels,'Stim');
stim_roi_act = fluor_roi_deconv - fluor_roi_taskpred_reduced(:,:,:,stim_regressor);

n_daysplit = 4;
trial_daysplit_idx = cell2mat(arrayfun(@(x) ...
    min(floor(linspace(1,n_daysplit+1,x)),n_daysplit)', ...
    trials_recording,'uni',false));

[learned_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);

[learned_stage_idx,~] = ndgrid(discretize(trial_learned_day,[-Inf,0,Inf]), ...
    1:length(t),1:n_rois);

[daysplit_idx,~,~] = ...
    ndgrid(trial_daysplit_idx,1:length(t),1:n_rois);

[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);

accum_learnedday_daysplit_idx = cat(4,learned_day_idx,daysplit_idx,t_idx,roi_idx,animal_idx);
accum_stage_idx = cat(4,learned_stage_idx,t_idx,roi_idx,animal_idx);

% (wheel avg: learned stage x t x animal)
wheel_learn_avg = accumarray( ...
    reshape(squeeze(accum_stage_idx(:,:,1,[1,2,4])),[],size(accum_stage_idx,4)-1), ...
    wheel_allcat(:), ...
    [max(learned_stage_idx(:)),length(t),length(animals)], ...
    @nanmean,NaN);

% (roi activity avg: learned stage x t x roi x animal)
raw_roi_act_learn_avg = accumarray( ...
    reshape(accum_stage_idx,[],size(accum_stage_idx,4)), ...
    fluor_roi_deconv(:), ...
    [max(learned_stage_idx(:)),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

stim_roi_act_learn_avg = accumarray( ...
    reshape(accum_stage_idx,[],size(accum_stage_idx,4)), ...
    stim_roi_act(:), ...
    [max(learned_stage_idx(:)),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

% (roi activity daysplit: learned day x (daysplit) x t x roi x animal)
stim_roi_act_learn_avg_daysplit = accumarray( ...
    reshape(accum_learnedday_daysplit_idx,[],size(accum_learnedday_daysplit_idx,4)), ...
    stim_roi_act(:), ...
    [max(learned_day_idx(:)),n_daysplit,length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

% Plot wheel by stage
stage_col = [0.7,0,0;0,0.7,0];
figure;
AP_errorfill(t, ...
    permute(nanmean(wheel_learn_avg,3),[2,1,3]), ...
    permute(AP_sem(wheel_learn_avg,3),[2,1,3]), ...
    stage_col);
xlabel('Time from stim (s)');
ylabel('Wheel velocity');
xline(0,'linestyle','--');

% Plot ROIs by stage (raw activity)
figure;
plot_rois = [1,6];
h = tiledlayout(length(plot_rois),2,'TileSpacing','compact','padding','compact');
for curr_roi_idx = 1:length(plot_rois)
    curr_l_roi = plot_rois(curr_roi_idx);
    curr_r_roi = curr_l_roi + size(wf_roi,1);

    % (plot left ROI w/ right stim)
    hl = nexttile;
    AP_errorfill(t, ...
        permute(nanmean(raw_roi_act_learn_avg(:,:,curr_l_roi,:),4),[2,1,3]), ...
        permute(AP_sem(raw_roi_act_learn_avg(:,:,curr_l_roi,:),4),[2,1,3]), ...
        stage_col);
    xlabel('Time from stim (s)');
    ylabel('\DeltaF/F_0');
    title(wf_roi(curr_l_roi).area);
    xline(0,'linestyle','--');

    % (plot right ROI with left stim)
    hr =nexttile;
    AP_errorfill(t, ...
        permute(nanmean(raw_roi_act_learn_avg(:,:,curr_r_roi,:),4),[2,1,3]), ...
        permute(AP_sem(raw_roi_act_learn_avg(:,:,curr_r_roi,:),4),[2,1,3]), ...
        stage_col);
    xlabel('Time from stim (s)');
    ylabel('\DeltaF/F_0');
    title(wf_roi(curr_r_roi).area);
    xline(0,'linestyle','--');

    % (link axes with same ROI)
    linkaxes([hl,hr],'xy');
end
title(h,'Raw activity');

% Plot ROIs by stage (reduced stim)
figure;
plot_rois = [1,6];
h = tiledlayout(length(plot_rois),2,'TileSpacing','compact','padding','compact');
for curr_roi_idx = 1:length(plot_rois)
    curr_l_roi = plot_rois(curr_roi_idx);
    curr_r_roi = curr_l_roi + size(wf_roi,1);

    % (plot left ROI w/ right stim)
    hl = nexttile;
    AP_errorfill(t, ...
        permute(nanmean(stim_roi_act_learn_avg(:,:,curr_l_roi,:),4),[2,1,3]), ...
        permute(AP_sem(stim_roi_act_learn_avg(:,:,curr_l_roi,:),4),[2,1,3]), ...
        stage_col);
    xlabel('Time from stim (s)');
    ylabel('\DeltaF/F_0');
    title(wf_roi(curr_l_roi).area);
    xline(0,'linestyle','--');

    % (plot right ROI with left stim)
    hr =nexttile;
    AP_errorfill(t, ...
        permute(nanmean(stim_roi_act_learn_avg(:,:,curr_r_roi,:),4),[2,1,3]), ...
        permute(AP_sem(stim_roi_act_learn_avg(:,:,curr_r_roi,:),4),[2,1,3]), ...
        stage_col);
    xlabel('Time from stim (s)');
    ylabel('\DeltaF/F_0');
    title(wf_roi(curr_r_roi).area);
    xline(0,'linestyle','--');

    % (link axes with same ROI)
    linkaxes([hl,hr],'xy');
end
title(h,'Reduced stim activity');

% Plot all trials (raw and reduced stim: by stage, sorted by reaction time)
plot_rois = [1,6];
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);
n_trial_smooth = 20;

figure;
h = tiledlayout(2,length(plot_rois),'TileSpacing','compact','padding','compact');
for curr_roi = plot_rois
    for curr_stage = 1:max(trial_learned_stage)
        use_trials = find(trial_learned_stage == curr_stage);
        [~,sort_idx] = sort(move_t(use_trials));
        curr_data_sort = fluor_roi_deconv(use_trials(sort_idx),:,curr_roi);
        curr_data_sort_smooth = convn(curr_data_sort, ...
            ones(n_trial_smooth,1)./n_trial_smooth,'same');
        
        nexttile;
        imagesc(t,[],curr_data_sort_smooth);hold on;
        colormap(brewermap([],'PrGn'));
        caxis(max(caxis)*[-1,1])
        xline(0,'color','r');
        plot(move_t(use_trials(sort_idx)),1:length(use_trials),'color',[0.6,0,0.6]);
        xlabel('Time from stim (s)');
        ylabel('Trial (rxn-time sorted)');
        title(sprintf('%s, stage %d',wf_roi(curr_roi).area,curr_stage));
    end
end
title(h,'Raw activity');

figure;
h = tiledlayout(2,length(plot_rois),'TileSpacing','compact','padding','compact');
for curr_roi = plot_rois
    for curr_stage = 1:max(trial_learned_stage)
        use_trials = find(trial_learned_stage == curr_stage);
        [~,sort_idx] = sort(move_t(use_trials));
        
        stim_regressor = strcmp(task_regressor_labels,'Stim');
        curr_data_sort = fluor_roi_deconv(use_trials(sort_idx),:,curr_roi) - ...
            fluor_roi_taskpred_reduced(use_trials(sort_idx),:,curr_roi,stim_regressor);        
        % (nanconv, and zero out NaNs for plotting)
        curr_data_sort_smooth = nanconv(curr_data_sort, ...
            ones(n_trial_smooth,1)./n_trial_smooth);
        curr_data_sort_smooth(isnan(curr_data_sort_smooth)) = 0;
        
        nexttile;
        imagesc(t,[],curr_data_sort_smooth);hold on;
        colormap(brewermap([],'PrGn'));
        caxis(max(caxis)*[-1,1])
        xline(0,'color','r');
        plot(move_t(use_trials(sort_idx)),1:length(use_trials),'color',[0.6,0,0.6]);
        xlabel('Time from stim (s)');
        ylabel('Trial (rxn-time sorted)');
        title(sprintf('%s, stage %d',wf_roi(curr_roi).area,curr_stage));
    end
end
title(h,'Reduced stim activity');


% (tmax activity: learned day x daysplit x roi x animal)
use_t = t >= 0 & t <= 0.5;
stim_roi_act_tmax_daysplit = ...
    squeeze(max(stim_roi_act_learn_avg_daysplit(:,:,use_t,:,:),[],3));

% Plot time-max by learned day
learned_day_x_range = minmax(learned_day_unique);
learned_day_x = [learned_day_x_range(1):learned_day_x_range(2)]';
learned_daysplit_x = learned_day_x + linspace(0,1,n_daysplit+1);

plot_learned_day = sum(~isnan(stim_roi_act_tmax_daysplit(:,1,1,:)),4) > min_n;

figure;
plot_rois = [1,5,6];
tiledlayout(length(plot_rois),1,'TileSpacing','compact','padding','compact');
for curr_roi_idx = 1:length(plot_rois)
    curr_l_roi = plot_rois(curr_roi_idx);

    nexttile;
   
    errorbar(reshape(learned_daysplit_x(plot_learned_day,:)',[],1), ...
        reshape(padarray(nanmean( ...
        stim_roi_act_tmax_daysplit(plot_learned_day,:,curr_l_roi,:),4), ...
        [0,1],NaN,'post')',[],1), ...
        reshape(padarray(AP_sem( ...
        stim_roi_act_tmax_daysplit(plot_learned_day,:,curr_l_roi,:),4), ...
        [0,1],NaN,'post')',[],1),'k','linewidth',2,'CapSize',0);

    title(wf_roi(curr_l_roi).area);
    xlabel('Learned day');
    ylabel('\DeltaF/F_0');
    xline(0,'linestyle','--');
    axis tight;
    xlim(xlim+[-0.5,0.5]);

end

% Plot ROIs overlaid, normalized by pre-learn average
norm_days = learned_day_x < 0;
stim_roi_act_tmax_daysplit_normval = ...
    nanmean(nanmean(stim_roi_act_tmax_daysplit(norm_days,:,:,:),2),1);
stim_roi_act_tmax_daysplit_norm = ...
    (stim_roi_act_tmax_daysplit-stim_roi_act_tmax_daysplit_normval)./ ...
    stim_roi_act_tmax_daysplit_normval;

figure;
plot_rois = [1,6];
errorbar( ...
    repmat(reshape(learned_daysplit_x(plot_learned_day,:)',[],1),1,length(plot_rois)), ...
    reshape(permute(padarray(nanmean( ...
    stim_roi_act_tmax_daysplit_norm(plot_learned_day,:,plot_rois,:),4), ...
    [0,1],NaN,'post'),[2,1,3]),[],length(plot_rois)), ...
    reshape(permute(padarray(AP_sem( ...
    stim_roi_act_tmax_daysplit_norm(plot_learned_day,:,plot_rois,:),4), ...
    [0,1],NaN,'post'),[2,1,3]),[],length(plot_rois)),'linewidth',2,'CapSize',0);
xlabel('Learned day');
ylabel('Fluorescence (normalized to pre-learn)');
xline(0,'linestyle','--');
legend({wf_roi(plot_rois).area},'location','nw')


%% Task - electrophysiology

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_task_ephys';

AP_load_trials_operant;

% Get animal and day index for each trial
trial_animal = cell2mat(arrayfun(@(x) ...
    x*ones(size(vertcat(wheel_all{x}{:}),1),1), ...
    [1:length(wheel_all)]','uni',false));
trial_day = cell2mat(cellfun(@(x) cell2mat(cellfun(@(curr_day,x) ...
    curr_day*ones(size(x,1),1),num2cell(1:length(x))',x,'uni',false)), ...
    wheel_all,'uni',false));

trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));

% Plot total average stim response
use_trials = true(size(trial_stim_allcat));
spacing = 0.5;

figure;
subplot(1,2,1);
mua_depth_centers = mua_depth_edges(1:end-1)+diff(mua_depth_edges)./2;
AP_stackplot(squeeze(nanmean(mua_depth_allcat(use_trials,:,:),1)), ...
    t,spacing,[],'k',mua_depth_centers);
line([0,0],ylim,'color','r');
line([0.5,0.5],ylim,'color','r');
xlabel('Time from stim');
ylabel('Cortical depth');

subplot(1,2,2);
[~,area_idx] = cellfun(@(x) ismember(x,mua_areas),mua_areas_cat,'uni',false);
area_recording_n = accumarray(cell2mat(area_idx),1);
plot_areas = area_recording_n == length(trials_recording);
AP_stackplot(squeeze(nanmean(mua_area_allcat(use_trials,:,plot_areas),1)), ...
    t,spacing,[],'k',mua_areas(plot_areas));
line([0,0],ylim,'color','r');
line([0.5,0.5],ylim,'color','r');
xlabel('Time from stim');
ylabel('Area');

% Get normalized task kernels
mua_taskpred_k_allcat_norm = arrayfun(@(regressor) ...
    squeeze(permute(cellfun(@(x) x{regressor}, ...
    cellfun(@(kernel_set,mua_norm) cellfun(@(kernel) ...
    kernel./(mua_norm/sample_rate),kernel_set,'uni',false), ...
    vertcat(mua_taskpred_k_all{:}),vertcat(mua_area_day_baseline{:}),'uni',false), ...
    'uni',false),[2,3,4,1])),1:length(task_regressor_labels),'uni',false);

mua_area_taskpred_k = cellfun(@(x) cellfun(@(x) ...
    nan(size(x,1),size(x,2),length(mua_areas)),x,'uni',false), ...
    mua_taskpred_k_allcat_norm,'uni',false);
for curr_recording = 1:length(mua_area_norm)
    [~,curr_area_idx] = ismember(mua_areas_cat{curr_recording},mua_areas);
    for curr_k = 1:length(task_regressor_labels)
        mua_area_taskpred_k{curr_k}{curr_recording}(:,:,curr_area_idx) = ...
            mua_taskpred_k_allcat_norm{curr_k}{curr_recording};
    end
end

mua_area_taskpred_k_avg = cellfun(@(x) nanmean(cat(4,x{:}),4), ...
    mua_area_taskpred_k,'uni',false);

% Plot kernels
plot_areas = area_recording_n == length(trials_recording);
figure;
h = tiledlayout('flow');
for curr_regressor = 1:length(mua_area_taskpred_k_avg)
    for curr_subregressor = 1:size(mua_area_taskpred_k_avg{curr_regressor},1)
        nexttile;
        plot(task_regressor_sample_shifts{curr_regressor}/sample_rate, ...
            squeeze(mua_area_taskpred_k_avg{curr_regressor}(curr_subregressor,:,plot_areas)));
        title(sprintf('%s (%d)',task_regressor_labels{curr_regressor},curr_subregressor));
        xline(0,'color','k','linestyle','--');
    end
end
linkaxes(allchild(h),'y');


%% Passive - electrophysiology

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_passive_ephys';

AP_load_trials_operant;

% Get animal and day index for each trial
trial_animal = cell2mat(arrayfun(@(x) ...
    x*ones(size(vertcat(wheel_all{x}{:}),1),1), ...
    [1:length(wheel_all)]','uni',false));
trial_day = cell2mat(cellfun(@(x) cell2mat(cellfun(@(curr_day,x) ...
    curr_day*ones(size(x,1),1),num2cell(1:length(x))',x,'uni',false)), ...
    wheel_all,'uni',false));

trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));

% Get trials with movement during stim to exclude
quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > 0,2);

% Plot total average stim response
stim_col = ['b','k','r'];
unique_stim = unique(trial_stim_allcat);

figure;
spacing = 0.8;
for curr_stim_idx = 1:length(unique_stim)
    use_trials = quiescent_trials &  ...
        trial_stim_allcat == unique_stim(curr_stim_idx);
    
    subplot(1,2,1); hold on;
    mua_depth_centers = mua_depth_edges(1:end-1)+diff(mua_depth_edges)./2;
    AP_stackplot(squeeze(nanmean(mua_depth_allcat(use_trials,:,:),1)), ...
        t,spacing,[],stim_col(curr_stim_idx),mua_depth_centers);
    line([0,0],ylim,'color',[0.5,0.5,0.5]);
    line([0.5,0.5],ylim,'color',[0.5,0.5,0.5]);
    xlabel('Time from stim');
    ylabel('Cortical depth');
    
    subplot(1,2,2); hold on;
    [~,area_idx] = cellfun(@(x) ismember(x,mua_areas),mua_areas_cat,'uni',false);
    area_recording_n = accumarray(cell2mat(area_idx),1);
    plot_areas = area_recording_n == length(trials_recording);
    AP_stackplot(squeeze(nanmean(mua_area_allcat(use_trials,:,plot_areas),1)), ...
        t,spacing,[],stim_col(curr_stim_idx),mua_areas(plot_areas));
    line([0,0],ylim,'color',[0.5,0.5,0.5]);
    line([0.5,0.5],ylim,'color',[0.5,0.5,0.5]);
    xlabel('Time from stim');
    ylabel('Area');
    
    % Scalebar
    line([-0.5,-0.5],[1,1.5],'color','m','linewidth',2);

end

% Plot overlaid animal mean with errorbars
use_trials = quiescent_trials & trial_stim_allcat == 3; 
[animal_idx,t_idx,area_idx] = ...
    ndgrid(trial_animal(use_trials),1:length(t),1:length(mua_areas));
mua_animal_avg = accumarray([animal_idx(:),t_idx(:),area_idx(:)], ...
    reshape(mua_area_allcat(use_trials,:,:),[],1), ...
    [length(animals),length(t),length(mua_areas)], ...
    @nanmean,NaN);

plot_areas = area_recording_n == length(trials_recording);
h = AP_errorfill(t',squeeze(nanmean(m(:,:,plot_areas),1)), ...
    squeeze(AP_sem(m(:,:,plot_areas),1)));
xline([0,0.5],'color','k');
xlabel('Time from stim (s)');
ylabel('\DeltaR/R_0');
legend(h,mua_areas(plot_areas));


%% Passive - V1 muscimol

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_passive_teto_muscimol';

AP_load_trials_operant;

% Get animal and day index for each trial
trial_animal = cell2mat(arrayfun(@(x) ...
    x*ones(size(vertcat(wheel_all{x}{:}),1),1), ...
    [1:length(wheel_all)]','uni',false));
trial_day = cell2mat(cellfun(@(x) cell2mat(cellfun(@(curr_day,x) ...
    curr_day*ones(size(x,1),1),num2cell(1:length(x))',x,'uni',false)), ...
    wheel_all,'uni',false));

trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));

% Get trials with movement during stim to exclude
quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > 0,2);

% Get average fluorescence by animal, day, stim
stim_unique = unique(trial_stim_allcat);
stim_v_avg = cell(length(animals),1);
stim_roi_avg = cell(length(animals),1);
stim_roi = cell(length(animals),1);
for curr_animal = 1:length(animals)
    for curr_day = 1:max(trial_day(trial_animal == curr_animal))
        for curr_stim_idx = 1:length(stim_unique)
            use_trials = quiescent_trials & ...
                trial_animal == curr_animal & ...
                trial_day == curr_day & ...
                trial_stim_allcat == stim_unique(curr_stim_idx);
            stim_v_avg{curr_animal}(:,:,curr_day,curr_stim_idx) = ...
                permute(nanmean(fluor_allcat_deconv(use_trials,:,:),1),[3,2,1]);
            stim_roi_avg{curr_animal}(:,:,curr_day,curr_stim_idx) = ...
                permute(nanmean(fluor_roi_deconv(use_trials,:,:),1),[3,2,1]);
            
            stim_roi{curr_animal}{curr_day}{curr_stim_idx} = ...
                fluor_roi_deconv(use_trials,:,:);
        end
    end
end

% Plot average pixels
stim_v_avg_stage = cat(5,stim_v_avg{:});

stim_px_avg_stage_movie = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),squeeze(nanmean(stim_v_avg_stage(:,:,:,3,:),5)));
AP_image_scroll(reshape(permute(stim_px_avg_stage_movie,[1,2,4,3]),size(U_master,1),[],length(t)),t);
axis image off;
colormap(brewermap([],'PrGn'));
caxis([-max(abs(caxis)),max(abs(caxis))]);

use_t = t >= 0.05 & t <= 0.2;
stim_px_avg_stage_tavg = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    squeeze(nanmean(nanmean(stim_v_avg_stage(:,use_t,:,:,:),2),5)));

figure;
tiledlayout(2,3,'TileSpacing','compact','padding','compact');
c = (max(stim_px_avg_stage_tavg(:)).*[-1,1])*0.5;
for curr_stage = 1:2
    for curr_stim = 1:3
        nexttile;
        imagesc(stim_px_avg_stage_tavg(:,:,curr_stage,curr_stim));
        axis image off;
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        colormap(brewermap([],'PrGn'));
        caxis(c);
        title(sprintf('Stage %d, stim %d',curr_stage,curr_stim));
    end
end

% Plot average ROIs
stim_roi_avg_stage = cat(5,stim_roi_avg{:});

figure;
plot_rois = [1,6];
stage_col = [0,0,0;1,0,0];
tiledlayout(length(plot_rois),2,'TileSpacing','compact','padding','compact');
for curr_roi_idx = 1:length(plot_rois)
    curr_l_roi = plot_rois(curr_roi_idx);
    curr_r_roi = curr_l_roi + size(wf_roi,1);

    % (plot left ROI w/ right stim)
    hl = nexttile;
    AP_errorfill(t, ...
        squeeze(nanmean(stim_roi_avg_stage(curr_l_roi,:,:,stim_unique == 1,:),5)), ...
        squeeze(AP_sem(stim_roi_avg_stage(curr_l_roi,:,:,stim_unique == 1,:),5)),stage_col);
    xlabel('Time from stim (s)');
    ylabel('\DeltaF/F_0');
    title(wf_roi(curr_l_roi).area);
    xline([0,0.5],'linestyle','--');

    % (plot right ROI with left stim)
    hr =nexttile;
    AP_errorfill(t, ...
        squeeze(nanmean(stim_roi_avg_stage(curr_r_roi,:,:,stim_unique == -1,:),5)), ...
        squeeze(AP_sem(stim_roi_avg_stage(curr_r_roi,:,:,stim_unique == -1,:),5)),stage_col);
    xlabel('Time from stim (s)');
    ylabel('\DeltaF/F_0');
    title(wf_roi(curr_r_roi).area);
    xline([0,0.5],'linestyle','--');

    % (link axes with same ROI)
    linkaxes([hl,hr],'xy');
end

% Plot ROI time-max pre/post muscimol
use_t = t >= 0 & t <= 0.2;
stim_roi_avg_stage_tmax = squeeze(max(stim_roi_avg_stage(:,use_t,:,:,:),[],2));

plot_rois = [1,6];
plot_stim = 3;
figure; hold on
curr_data = squeeze(stim_roi_avg_stage_tmax(plot_rois,:,plot_stim,:));
plot(squeeze(curr_data(1,2,:)),squeeze(curr_data(2,2,:)),'.k','MarkerSize',20);
plot(squeeze(curr_data(1,1,:)),squeeze(curr_data(2,1,:)),'.r','MarkerSize',20);
xlabel(wf_roi(plot_rois(1)).area);
ylabel(wf_roi(plot_rois(2)).area);
legend({'Washout','V1 Muscimol'},'location','nw');





