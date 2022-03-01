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
plot_col = copper(length(plot_days));
h = AP_errorfill(rxn_bin_centers,nanmean(animal_rxn(plot_days,:,:),3)', ...
    AP_sem(animal_rxn(plot_days,:,:),3)',plot_col);
xlabel('Reaction time (s)');
ylabel('Probability');
xline([0.1,0.2],'r')
legend(h,cellfun(@(x) sprintf('Day %d',x),num2cell(plot_days),'uni',false));

% Total (trial/alt number-matched)
alt_rxn_matched = ...
    cellfun(@(x) cellfun(@(x) cell2mat(cellfun(@(x) ...
    datasample(x,min(length(x),1)),x,'uni',false)),x,'uni',false), ...
    {bhv.alt_stim_move_t},'uni',false);

animal_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(vertcat(x{use}),rxn_bins,'normalization','pdf'), ...
    use_days,{bhv.stim_move_t},'uni',false)');

animal_alt_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(vertcat(x{use}),rxn_bins,'normalization','pdf'), ...
    use_days,alt_rxn_matched,'uni',false)');

figure;
subplot(2,1,1); hold on;
histogram('BinEdges',rxn_bins,'BinCounts',nanmean(animal_rxn,1),'EdgeColor','none')
histogram('BinEdges',rxn_bins,'BinCounts',nanmean(animal_alt_rxn,1),'EdgeColor','none')
xlabel('Reaction time');
ylabel('PDF');
legend({'Measured','Null (n trial matched)'});
subplot(2,1,2);
plot(rxn_bin_centers,nanmean(animal_rxn - animal_alt_rxn,1));
line(xlim,[0,0],'color','r');
ylabel('Distribution difference');
linkaxes(get(gcf,'children'),'xy');

% Fraction rxn within boundary
rxn_frac_window = [0.1,0.2];
rxn_frac_pad = cell2mat(cellfun(@(x) padarray(cellfun(@(x) ...
    nanmean(x > rxn_frac_window(1) & x < rxn_frac_window(2)),x), ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_move_t},use_days,'uni',false),'uni',false));

alt_rxn_frac_pad = cell2mat(cellfun(@(x) padarray(cellfun(@(x) ...
    nanmean(cell2mat(x) > rxn_frac_window(1) & cell2mat(x) < rxn_frac_window(2)),x), ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.alt_stim_move_t},use_days,'uni',false),'uni',false));

% (plot relative to training day and learned day)
learned_day = cellfun(@(x) find(x,1),{bhv.learned_days});

figure; 

subplot(1,2,1);hold on;
plot(rxn_frac_pad(plot_days,:),'color',[0.5,0.5,0.5]);
errorbar(nanmean(rxn_frac_pad(plot_days,:),2), ...
    AP_sem(rxn_frac_pad(plot_days,:),2),'k','linewidth',2,'CapSize',0);
errorbar(nanmean(alt_rxn_frac_pad(plot_days,:),2), ...
AP_sem(alt_rxn_frac_pad(plot_days,:),2),'r','linewidth',2','CapSize',0);

learned_day_x = [1:max_days]'-learned_day;
[rxn_grp_mean,learned_day_grp,learned_day_n] = ...
    grpstats(rxn_frac_pad(:),learned_day_x(:),{'nanmean','gname','numel'});
learned_day_grp = cellfun(@str2num,learned_day_grp);

subplot(1,2,2); hold on;
plot(learned_day_x,rxn_frac_pad,'color',[0.5,0.5,0.5]);
plot(learned_day_grp,rxn_grp_mean,'k','linewidth',2);
xlabel('Learned day');
ylabel('Rxn frac');
line([0,0],ylim,'color','k','linestyle','--');


% Day-split reaction time fraction
n_daysplit = 4;

% Index: within-day split, day, animal
trial_split_idx = cellfun(@(x,use_days,animal_num) ...
    cellfun(@(x,day) ...
    [min(floor(linspace(1,n_daysplit+1,length(x))),n_daysplit)', ...
    repmat(day,length(x),1),repmat(animal_num,length(x),1)], ...
    x(use_days),num2cell(1:sum(use_days))','uni',false), ...
    {bhv.stim_move_t},use_days,num2cell(1:length(bhv)),'uni',false);

stim_move_t_cat = cell2mat(cellfun(@(x,use_days) cell2mat(x(use_days)), ...
    {bhv.stim_move_t},use_days,'uni',false)');
alt_stim_move_t_cat = cell2mat(cellfun(@(x,use_days) cell2mat(x(use_days)), ...
    alt_rxn_matched,use_days,'uni',false)');
% alt_stim_move_t_cat = cell2mat(cellfun(@(x,use_days) ...
%     cell2mat(cellfun(@(x) cellfun(@nanmedian,x),x(use_days),'uni',false)), ...
%     {bhv.alt_stim_move_t},use_days,'uni',false)');

stim_move_t_stimtime = accumarray(cell2mat(cat(1,trial_split_idx{:})), ...
    stim_move_t_cat > rxn_frac_window(1) & ...
    stim_move_t_cat < rxn_frac_window(2)', ...
    [n_daysplit,max_days,length(bhv)],@nanmean,NaN);
alt_stim_move_t_stimtime = accumarray(cell2mat(cat(1,trial_split_idx{:})), ...
    alt_stim_move_t_cat > rxn_frac_window(1) & ...
    alt_stim_move_t_cat < rxn_frac_window(2)', ...
    [n_daysplit,max_days,length(bhv)],@nanmean,NaN);

% (plot relative to training day and learned day)
figure; 

subplot(1,2,1);hold on;
stim_move_t_stimtime_long = reshape(padarray(stim_move_t_stimtime,[1,0,0],NaN,'post'),[],length(animals));
alt_stim_move_t_stimtime_long = reshape(padarray(alt_stim_move_t_stimtime,[1,0,0],NaN,'post'),[],length(animals));

p1 = errorbar([1:size(stim_move_t_stimtime_long,1)]/(n_daysplit+1), ...
    nanmean(stim_move_t_stimtime_long,2),AP_sem(stim_move_t_stimtime_long,2),'k','linewidth',2,'CapSize',0);
p2 = errorbar([1:size(alt_stim_move_t_stimtime_long,1)]/(n_daysplit+1), ...
    nanmean(alt_stim_move_t_stimtime_long,2),AP_sem(alt_stim_move_t_stimtime_long,2),'r','linewidth',2,'CapSize',0);

xlabel('Day');
ylabel('Frac rxn 100-200 ms');
legend([p1,p2],{'Measured','Null'},'location','nw');

% (relative to learned day)
learned_day_x_daysplit = cell2mat(cellfun(@(x) x+(0:n_daysplit)'/(n_daysplit+1), ...
    num2cell(learned_day_x),'uni',false));

[rxn_grp_mean_daysplit,rxn_grp_sem_daysplit,learned_day_grp_daysplit,learned_day_n_daysplit] = ...
    grpstats(stim_move_t_stimtime_long(:),learned_day_x_daysplit(:), ...
    {'nanmean','sem','gname','numel'});
[alt_rxn_grp_mean_daysplit,alt_rxn_grp_sem_daysplit] = ...
    grpstats(alt_stim_move_t_stimtime_long(:),learned_day_x_daysplit(:), ...
    {'nanmean','sem'});

learned_day_grp_daysplit = cellfun(@str2num,learned_day_grp_daysplit);

plot_daysplit = learned_day_n_daysplit >= min_n | isnan(rxn_grp_mean_daysplit);

subplot(1,2,2); hold on;
errorbar(learned_day_grp_daysplit(plot_daysplit),rxn_grp_mean_daysplit(plot_daysplit), ...
    rxn_grp_sem_daysplit(plot_daysplit),'k','linewidth',2,'CapSize',0);
errorbar(learned_day_grp_daysplit(plot_daysplit),alt_rxn_grp_mean_daysplit(plot_daysplit), ...
    alt_rxn_grp_sem_daysplit(plot_daysplit),'r','linewidth',2,'CapSize',0);
line([0,0],ylim,'color','k','linestyle','--');
xlabel('Learned day');

%%%%% TESTING (get many-sampled alt rxn and comparable rxn)

% (exclude trials without alts: rare, from wheel click issues)
use_rxn = cellfun(@(x) cellfun(@(x) cellfun(@(x) ...
    ~isempty(x),x),x,'uni',false),{bhv.alt_stim_move_t},'uni',false);

rxn_measured = cellfun(@(rxn,use_days,use_trials) ...
    cellfun(@(rxn,use_trials) rxn(use_trials),rxn(use_days),use_trials(use_days),'uni',false), ...
    {bhv.stim_move_t},use_days,use_rxn,'uni',false)';

n_resample = 10;
rxn_alt = cellfun(@(rxn,use_days,use_trials) ...
    cellfun(@(rxn,use_trials) ...
    cell2mat(cellfun(@(x) datasample(x,n_resample)',rxn(use_trials),'uni',false)), ...
    rxn(use_days),use_trials(use_days),'uni',false), ...
    {bhv.alt_stim_move_t},use_days,use_rxn,'uni',false)';


%% ~~~~~ Load passive

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_passive_tetO';
AP_load_trials_wf;
n_naive = 3; % (number of naive passive-only days)

% Load behavior and get learned day
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
bhv_fn = [data_path filesep 'bhv_teto'];
load(bhv_fn);
learned_day = cellfun(@(x) find(x,1),{bhv.learned_days});

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
stim_v_avg = cell(size(animals));
stim_roi_avg = cell(size(animals));
stim_roi = cell(size(animals));
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

% Average V's by naive/prelearn/postlearn
stim_v_avg_stage = cell2mat(permute(cellfun(@(x,ld) ...
    cat(3,nanmean(x(:,:,1:n_naive,:),3), ...
    nanmean(x(:,:,n_naive+1:ld-1,:),3), ...
    nanmean(x(:,:,ld:end,:),3)),stim_v_avg,num2cell(learned_day), ...
    'uni',false),[1,3,4,5,2]));

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
    'uni',false),[1,3,4,5,2]));

figure;
plot_rois = [1,5,6];
stage_col = [0,0,0;0.7,0,0;0,0.7,0];
tiledlayout(length(plot_rois),2,'TileSpacing','tight','padding','tight');
for curr_roi_idx = 1:length(plot_rois)
    curr_l_roi = plot_rois(curr_roi_idx);
    curr_r_roi = curr_l_roi + size(wf_roi,1);

    % (plot left ROI w/ right stim)
    hl = nexttile;
    AP_errorfill(t, ...
        squeeze(nanmean(stim_roi_avg_stage(curr_l_roi,:,:,stim_unique == 1,:),5)), ...
        squeeze(AP_sem(stim_roi_avg_stage(curr_l_roi,:,:,stim_unique == 1,:),5)),stage_col);
    title(wf_roi(curr_l_roi).area);
    xline([0,0.5],'linestyle','--');

    % (plot right ROI with left stim)
    hr =nexttile;
    AP_errorfill(t, ...
        squeeze(nanmean(stim_roi_avg_stage(curr_r_roi,:,:,stim_unique == -1,:),5)), ...
        squeeze(AP_sem(stim_roi_avg_stage(curr_r_roi,:,:,stim_unique == -1,:),5)),stage_col);
    title(wf_roi(curr_r_roi).area);
    xline([0,0.5],'linestyle','--');

    % (link axes with same ROI)
    linkaxes([hl,hr],'xy');
end



% TODO: get time window activity aligned to learning



























