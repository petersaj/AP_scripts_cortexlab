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

% Grab day index of [V1 muscimol, washout] and retinotopy
muscimol_v1_days = nan(length(bhv),2);
muscimol_v1_retinotopy = cell(size(muscimol_v1_days));
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

    % Store days
     muscimol_v1_days(curr_animal,:) = ...
             [curr_v1_muscimol_dayidx,curr_v1_washout_dayidx];

     % Store retinotopy
    muscimol_v1_retinotopy(curr_animal,:) = ...
        muscimol(muscimol_animal_idx).vfs([curr_v1_muscimol,curr_v1_washout]);

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


% (set use_days to copy code from above)
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
for curr_cond = 1:2
    subplot(1,2,curr_cond); hold on
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
linkaxes(get(gcf,'children'),'xy');



%% Passive

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

% Average V's by naive/prelearn/postlearn, plot time average
stim_v_avg_stage = cell2mat(permute(cellfun(@(x,ld) ...
    cat(3,nanmean(x(:,:,1:n_naive,:),3), ...
    nanmean(x(:,:,n_naive+1:ld-1,:),3), ...
    nanmean(x(:,:,ld:end,:),3)),stim_v_avg,num2cell(learned_day), ...
    'uni',false),[2,3,4,5,1]));

stim_px_avg_stage_movie = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),squeeze(nanmean(stim_v_avg_stage(:,:,:,3,:),5)));
AP_image_scroll(reshape(permute(stim_px_avg_stage_movie,[1,2,4,3]),size(U_master,1),[],length(t)),t);
axis image off;
colormap(brewermap([],'PrGn'));
caxis([-max(abs(caxis)),max(abs(caxis))]);

use_t = t >= 0.05 & t <= 0.2;
stim_px_avg_stage_tavg = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    squeeze(nanmean(nanmean(stim_v_avg_stage(:,use_t,:,:,:),2),5)));

figure;
tiledlayout(3,3,'TileSpacing','tight','padding','tight');
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
    cat(3,nanmean(x(:,:,1:n_naive,:),3), ...
    nanmean(x(:,:,n_naive+1:ld-1,:),3), ...
    nanmean(x(:,:,ld:end,:),3)),stim_roi_avg,num2cell(learned_day), ...
    'uni',false),[2,3,4,5,1]));

figure;
plot_rois = [1,5,6];
stage_col = [0,0,0;0.7,0,0;0,0.7,0];
tiledlayout(length(plot_rois),2,'TileSpacing','tight','padding','compact');
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

% (grab indicies)
[roi_idx,learned_day_idx,stim_idx] = ...
    ndgrid(1:n_rois,cat(1,learned_day_animal_x{:}),1:3);
[~,trained_day_idx,~] = ...
    ndgrid(1:n_rois,cat(1,trained_day_animal_x{:}),1:3);
training_data = trained_day_idx > n_naive; % (exclude naive passive-only)
% (re-index learned day to be positive for accumarray)
learned_day_minidx = learned_day_idx-min(learned_day_idx(:))+1;

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
    accumarray([roi_idx(training_data),learned_day_minidx(training_data),stim_idx(training_data)], ...
    stim_roi_tmax_daycat(training_data), ...
    [n_rois,max(learned_day_minidx(:)),length(stim_unique)],@nanmean,NaN('single'));
stim_roi_tmax_learn_sem = ...
    accumarray([roi_idx(training_data),learned_day_minidx(training_data),stim_idx(training_data)], ...
    stim_roi_tmax_daycat(training_data), ...
    [n_rois,max(learned_day_minidx(:)),length(stim_unique)],@AP_sem,NaN('single'));
% (get number of elements)
stim_roi_tmax_learn_numel = ...
    accumarray([roi_idx(training_data),learned_day_minidx(training_data),stim_idx(training_data)], ...
    stim_roi_tmax_daycat(training_data), ...
    [n_rois,max(learned_day_minidx(:)),length(stim_unique)],@numel,NaN);
% (get corresponding learned day x values and days to plot)
learned_day_x_range = minmax(learned_day_idx);
learned_day_x = learned_day_x_range(1):learned_day_x_range(2);
plot_learned_day = squeeze(stim_roi_tmax_learn_numel(1,:,1)) > min_n;

figure;
plot_rois = [1,5,6];
tiledlayout(length(plot_rois),2,'TileSpacing','tight','padding','compact');
for curr_roi_idx = 1:length(plot_rois)
    curr_l_roi = plot_rois(curr_roi_idx);
    curr_r_roi = curr_l_roi + size(wf_roi,1);

    stim_col = [0,0,1;0.5,0.5,0.5;1,0,0];

    % Plot L ROI
    hl = nexttile; hold on;
    % (naive)
    AP_errorfill(min(learned_day_x(plot_learned_day))+[-n_naive:-1], ...
        squeeze(stim_roi_tmax_naive_mean(curr_l_roi,1:n_naive,:)), ...
        squeeze(stim_roi_tmax_naive_sem(curr_l_roi,1:n_naive,:)),stim_col);
    % (training)
    AP_errorfill(learned_day_x(plot_learned_day), ...
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
    AP_errorfill(min(learned_day_x(plot_learned_day))+[-n_naive:-1], ...
        squeeze(stim_roi_tmax_naive_mean(curr_r_roi,1:n_naive,:)), ...
        squeeze(stim_roi_tmax_naive_sem(curr_r_roi,1:n_naive,:)),stim_col);
    % (training)
    AP_errorfill(learned_day_x(plot_learned_day), ...
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


%% Task

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

% Get number of trials by animal/recording
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));

% Get ROI average/time-max stim reduced activity by learned day (daysplit)
stim_regressor = strcmp(task_regressor_labels,'Stim');
stim_roi_act = fluor_roi_deconv - fluor_roi_taskpred_reduced(:,:,:,stim_regressor);

n_daysplit = 4;
trial_daysplit_idx = cell2mat(arrayfun(@(x) ...
    min(floor(linspace(1,n_daysplit+1,x)),n_daysplit)', ...
    trials_recording,'uni',false));

[learned_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_learned_day,1:length(t),1:n_rois);
learned_day_minidx = learned_day_idx-min(learned_day_idx(:))+1;

[daysplit_idx,~,~] = ...
    ndgrid(trial_daysplit_idx,1:length(t),1:n_rois);

[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);

% (roi activity: learned day x daysplit x t x roi x animal)
stim_roi_act_learn_avg = accumarray( ...
    [learned_day_minidx(:),daysplit_idx(:),t_idx(:),roi_idx(:),animal_idx(:)], ...
    stim_roi_act(:), ...
    [max(learned_day_minidx(:)),n_daysplit,length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

% (tmax activity: learned day x daysplit x roi x animal)
use_t = t >= 0 & t <= 0.5;
stim_roi_act_tmax_daysplit = ...
    squeeze(max(stim_roi_act_learn_avg(:,:,use_t,:,:),[],3));

% Plot time-max by learned day
learned_day_x_range = minmax(learned_day_idx);
learned_day_x = [learned_day_x_range(1):learned_day_x_range(2)]';
learned_daysplit_x = learned_day_x + linspace(0,1,n_daysplit+1);

plot_learned_day = sum(~isnan(stim_roi_act_tmax_daysplit(:,1,1,:)),4) > min_n;

figure;
plot_rois = [1,5,6];
tiledlayout(length(plot_rois),1,'TileSpacing','tight','padding','compact');
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



