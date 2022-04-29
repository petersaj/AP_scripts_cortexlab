%% Operant learning widefield paper figures
% Preprocessing done in AP_operant_learning_preprocessing


%% Plot example performance

animal = 'AP106';
protocol = 'AP_stimWheelRight';
experiments = AP_find_experiments(animal,protocol);
plot_days = [1,3,7];

h = tiledlayout(length(plot_days),1,'TileSpacing','compact','padding','compact');
for curr_day = plot_days
 
    day = experiments(curr_day).day;
    experiment = experiments(curr_day).experiment;

    load_parts.imaging = false;
    AP_load_experiment;
    
    t = Timeline.rawDAQTimestamps;
    
    % (minimum ITI: new trial + trial quiescence)
    min_iti_t = signals_events.newTrialTimes + ...
        signals_events.trialQuiescenceValues;
    
    plot_t = [50,100];
    plot_t_idx = t > plot_t(1) & t < plot_t(2);
    plot_stim_idx = find(stimOn_times > plot_t(1) & stimOn_times < plot_t(2))';
    plot_min_iti_t_idx = find(min_iti_t > plot_t(1) & min_iti_t < plot_t(2));
    plot_reward_idx = find(reward_t_timeline > plot_t(1) & reward_t_timeline < plot_t(2));
    
    nexttile; hold on;
    line_height = 0.1;
    line_start = 0.2;
    plot(t(plot_t_idx),wheel_velocity(plot_t_idx),'k');
    for i = plot_stim_idx
        line(repmat(stimOn_times(i),2,1), ...
            line_start+line_height*1+[0,line_height],'color','r','linewidth',2);
    end
    for i = plot_reward_idx'
        line(repmat(reward_t_timeline(i),2,1), ...
            line_start+line_height*0+[0,line_height],'color','b','linewidth',2);
    end
    
    title(sprintf('Day %d',curr_day));
    
end

linkaxes(allchild(h),'y');


%% Behavior - [LOAD DATA]
% (adapted from AP_operant_behavior)

% Load behavior
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
bhv_fn = [data_path filesep 'bhv_teto'];
load(bhv_fn);
animals = {bhv.animal};

% Grab learned day
learned_day = cellfun(@(x) find(x,1),{bhv.learned_days})';

% Grab window for fast/learned reaction times
rxn_window = bhv(1).learned_days_rxn_window;

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


%% ^^ Behavior - plot all

% Set bins for reaction time histograms
rxn_bins = [0:0.01:0.5];
rxn_bin_centers = rxn_bins(1:end-1) + diff(rxn_bins)./2;

% Get reaction time histograms by day
animal_rxn = nan(max_days,length(rxn_bin_centers),length(animals));
for curr_animal = 1:length(bhv)
    for curr_day = find(use_days{curr_animal})'
        animal_rxn(curr_day,:,curr_animal,:) = ...
            histcounts(bhv(curr_animal).stim_move_t{curr_day}, ...
            rxn_bins,'normalization','probability');
    end
end

% Concatenate all
rxn_measured_cat = cell2mat(cellfun(@cell2mat,rxn_measured,'uni',false));
rxn_alt_cat = cell2mat(cellfun(@cell2mat,rxn_alt,'uni',false));

% Index each trial for daysplit/day/animal and concat
n_daysplit = 3;
trial_split_idx = cellfun(@(x,animal_num) ...
    cellfun(@(x,day) ...
    [min(floor(linspace(1,n_daysplit+1,length(x))),n_daysplit)', ...
    repmat(day,length(x),1),repmat(animal_num,length(x),1)], ...
    x,num2cell(1:length(x))','uni',false), ...
    rxn_measured,num2cell(1:length(bhv))','uni',false);
trial_split_idx_cat = ...
    cell2mat(cellfun(@(x) cat(1,x{:}),trial_split_idx,'uni',false));

% Get histogram of reaction times in day groups (across animals)
animal_idx = cell2mat(cellfun(@(x,animal) repmat(animal,length(cell2mat(x)),1), ...
    rxn_measured,num2cell(1:length(bhv))','uni',false));
day_idx = cell2mat(cellfun(@(x) cell2mat(cellfun(@(x,day) ...
    repmat(day,length(x),1),x,num2cell(1:length(x))','uni',false)), ...
    rxn_measured,'uni',false));

day_grps = [1,3,5,7,Inf];
day_grp = discretize(day_idx,day_grps);

figure;
for curr_day_grp = 1:length(day_grps)-1

    animal_rxn_measured_cathist = cell2mat(arrayfun(@(x) ...
        histcounts(rxn_measured_cat( ...
        animal_idx == x & ...
        day_grp == curr_day_grp), ...
        rxn_bins,'normalization','probability')',1:length(bhv),'uni',false));

    animal_rxn_alt_cathist = cell2mat(permute( ...
        arrayfun(@(rep) cell2mat(arrayfun(@(x) ...
        histcounts(rxn_alt_cat( ...
        animal_idx == x & ...
        day_grp == curr_day_grp,rep), ...
        rxn_bins,'normalization','probability')',1:length(bhv),'uni',false)), ...
        1:n_rxn_altsample,'uni',false),[1,3,2]));

    % Plot measured/null histogram
    subplot(length(day_grps),1,curr_day_grp); hold on;
    day_grp_col = copper(length(day_grps)-1);

    animal_rxn_alt_cathist_ci = ...
        squeeze(prctile(nanmean(animal_rxn_alt_cathist,2),[5,95],3));
    AP_errorfill(rxn_bin_centers,nanmean(nanmean(animal_rxn_alt_cathist,2),3), ...
        animal_rxn_alt_cathist_ci,[0.5,0.5,0.5],[],false);
    plot(rxn_bin_centers,nanmean(animal_rxn_measured_cathist,2),'k','linewidth',2);
    xlabel('Reaction time');
    ylabel('Frequency');
    title(sprintf('Day %d-%d',day_grps(curr_day_grp),day_grps(curr_day_grp+1)-1));

    % Plot cumulative measured-null 
    subplot(length(day_grps),1,length(day_grps)); hold on;

    rxn_hist_altdiff_cumsum = cumsum(animal_rxn_measured_cathist - ...
        nanmean(animal_rxn_alt_cathist,3),1);
    AP_errorfill(rxn_bin_centers,nanmean(rxn_hist_altdiff_cumsum,2), ...
        AP_sem(rxn_hist_altdiff_cumsum,2),day_grp_col(curr_day_grp,:));
    xlabel('Reaction time');
    ylabel('Cumulative meas-null');

end

% Plot histogram by day heatmap
figure;
imagesc([],rxn_bin_centers,nanmean(animal_rxn(plot_days,:,:),3)');
colormap(1-gray);
h = colorbar;ylabel(h,'Probability');
xlabel('Day');
ylabel('Reaction time');

% UNUSED: MAD
% Get reaction times MAD: whole day
% (exclude jumping-the-gun < 0.1s)
% rxn_measured_mad = accumarray(trial_split_idx_cat(:,2:end), ...
%     rxn_measured_cat.*AP_nanout(rxn_measured_cat < 0.1), ...
%     [max_days,length(bhv)],@(x) mad(x,1),NaN);
% rxn_alt_mad = cell2mat(permute(arrayfun(@(x) ...
%     accumarray(trial_split_idx_cat(:,2:end), ...
%     rxn_alt_cat(:,x).*AP_nanout(rxn_measured_cat < 0.1), ...
%     [max_days,length(bhv)],@(x) mad(x,1),NaN), ...
%     1:n_rxn_altsample,'uni',false),[1,3,2]));
% 
% figure; 
% 
% subplot(1,2,1);hold on
% rxn_alt_mad_ci = squeeze(prctile(nanmedian(rxn_alt_mad,2),[5,95],3));
% AP_errorfill([],nanmedian(rxn_alt_mad_ci(plot_days,:),2), ...
%     rxn_alt_mad_ci(plot_days,:),[0.5,0.5,0.5],[],false);
% errorbar(nanmedian(rxn_measured_mad(plot_days,:),2), ...
%     mad(rxn_measured_mad(plot_days,:),1,2),'k','CapSize',0,'linewidth',2)
% xlabel('Training day');
% ylabel('Reaction time MAD (s)');
% axis tight;
% xlim(xlim + [-0.5,0.5]);
% 
% subplot(1,2,2);hold on
% rxn_measured_mad_altdiff = rxn_measured_mad - nanmean(rxn_alt_mad,3);
% rxn_alt_mad_altdiff = rxn_alt_mad - nanmean(rxn_alt_mad,3);
% 
% rxn_alt_mad_altdiff_ci = squeeze(prctile(nanmedian(rxn_alt_mad_altdiff,2),[5,95],3));
% AP_errorfill([],nanmedian(rxn_alt_mad_altdiff_ci(plot_days,:),2), ...
%     rxn_alt_mad_altdiff_ci(plot_days,:),[0.5,0.5,0.5],[],false);
% errorbar(nanmedian(rxn_measured_mad_altdiff(plot_days,:),2), ...
%     mad(rxn_measured_mad_altdiff(plot_days,:),1,2),'k','CapSize',0,'linewidth',2)
% xlabel('Training day');
% ylabel('Reaction time MAD-null (s)');
% axis tight;
% xlim(xlim + [-0.5,0.5]);

% Reaction time median: whole day
% (exclude too-fast rxn < 0.1)
rxn_measured_med = accumarray(trial_split_idx_cat(:,2:end), ...
    rxn_measured_cat.*AP_nanout(rxn_measured_cat < 0.1), ...
    [max_days,length(bhv)],@(x) nanmedian(x),NaN);
rxn_alt_med = cell2mat(permute(arrayfun(@(x) ...
    accumarray(trial_split_idx_cat(:,2:end), ...
    rxn_alt_cat(:,x).*AP_nanout(rxn_alt_cat(:,x) < 0.1), ...
    [max_days,length(bhv)],@(x) nanmedian(x),NaN), ...
    1:n_rxn_altsample,'uni',false),[1,3,2]));

figure; 

subplot(1,2,1,'YScale','log');hold on
rxn_alt_mad_ci = squeeze(prctile(nanmean(rxn_alt_med,2),[5,95],3));
AP_errorfill([],nanmean(rxn_alt_mad_ci(plot_days,:),2), ...
    rxn_alt_mad_ci(plot_days,:),[0.5,0.5,0.5],[],false);
errorbar(nanmean(rxn_measured_med(plot_days,:),2), ...
    AP_sem(rxn_measured_med(plot_days,:),2),'k','CapSize',0,'linewidth',2)
xlabel('Training day');
ylabel('Median reaction time (s)');
axis tight;
xlim(xlim + [-0.5,0.5]);

subplot(1,2,2);hold on
rxn_measured_med_altdiff = rxn_measured_med - nanmean(rxn_alt_med,3);
rxn_alt_med_altdiff = rxn_alt_med - nanmean(rxn_alt_med,3);

rxn_alt_med_altdiff_ci = squeeze(prctile(nanmean(rxn_alt_med_altdiff,2),[5,95],3));
AP_errorfill([],nanmean(rxn_alt_med_altdiff_ci(plot_days,:),2), ...
    rxn_alt_med_altdiff_ci(plot_days,:),[0.5,0.5,0.5],[],false);
errorbar(nanmedian(rxn_measured_med_altdiff(plot_days,:),2), ...
    AP_sem(rxn_measured_med_altdiff(plot_days,:),2),'k','CapSize',0,'linewidth',2)
xlabel('Training day');
ylabel('Median reaction time (meas-null) (s)');
axis tight;
xlim(xlim + [-0.5,0.5]);


% Plot fraction of reaction times within window: whole day
rxn_measured_prct = accumarray(trial_split_idx_cat(:,2:end), ...
    rxn_measured_cat >= rxn_window(1) & ...
    rxn_measured_cat <= rxn_window(2), ...
    [max_days,length(bhv)],@(x) nanmean(x)*100,NaN);
rxn_alt_prct = cell2mat(permute(arrayfun(@(x) ...
    accumarray(trial_split_idx_cat(:,2:end), ...
    rxn_alt_cat(:,x) >= rxn_window(1) & ...
    rxn_alt_cat(:,x) <= rxn_window(2), ...
    [max_days,length(bhv)],@(x) nanmean(x)*100,NaN), ...
    1:n_rxn_altsample,'uni',false),[1,3,2]));

figure; hold on
rxn_alt_frac_ci = squeeze(prctile(nanmean(rxn_alt_prct,2),[5,95],3));
AP_errorfill([],nanmean(rxn_alt_frac_ci(plot_days,:),2), ...
    rxn_alt_frac_ci(plot_days,:),[0.5,0.5,0.5],[],false);
% plot(rxn_measured_prct(plot_days,:),'color',[0.5,0.5,0.5])
errorbar(nanmean(rxn_measured_prct(plot_days,:),2), ...
    AP_sem(rxn_measured_prct(plot_days,:),2),'k','CapSize',0,'linewidth',2)
% plot(nanmean(rxn_measured_prct(plot_days,:),2),'k','linewidth',2)
xlabel('Training day');
ylabel(sprintf('Fast reaction times: %d-%dms (%%)',1000*rxn_window(1),1000*rxn_window(2)));
axis tight;
xlim(xlim + [-0.5,0.5]);
ylim([0,100]);

% Plot all animals separately
figure('Name',sprintf('rxn: %.2g-%.2gms (%%)',rxn_window(1),rxn_window(2)))
h = tiledlayout('flow','TileSpacing','compact','padding','compact');
for curr_animal = 1:length(animals)    
    curr_rxn = rxn_measured_prct(:,curr_animal) - ...
        nanmean(rxn_alt_prct(:,curr_animal,:),3);

    curr_alt_rxn = rxn_alt_prct(:,curr_animal,:) - ...
        nanmean(rxn_alt_prct(:,curr_animal,:),3);

    curr_ci = permute(prctile(curr_alt_rxn,[5,95],3),[1,3,2]);

    nexttile; hold on;
    
    hx = xline(learned_day(curr_animal),'color','r','linewidth',2);
    warning off;shx = struct(hx);shx.Edge.Layer = 'back';warning on;

    AP_errorfill([],nanmean(curr_ci,2),curr_ci,[0.5,0.5,0.5],[],false);
    stairs(curr_rxn,'k','linewidth',2);
    xlabel('Training day');
    ylabel('Rxn-null (%)');
    title(animals{curr_animal});
end
linkaxes(allchild(h),'xy');

% Plot fast reaction times aligned to learning (whole day)
learned_day_x = [1:max_days]'-learned_day';

rxn_measured_prct_altdiff = rxn_measured_prct - nanmean(rxn_alt_prct,3);

[rxn_altdiff_learn_mean,rxn_altdiff_learn_sem,learned_day_grp,learned_day_n] = ...
    grpstats(rxn_measured_prct_altdiff(:),learned_day_x(:), ...
    {'nanmean','sem','gname','numel'});
learned_day_grp = cellfun(@str2num,learned_day_grp);
plot_learned = learned_day_n >= min_n;

rxn_alt_prct_altdiff_ci = nan(length(learned_day_grp),2);
for curr_learned_day = learned_day_grp'
    curr_alt_prct = ...
        reshape(rxn_alt_prct(repmat(learned_day_x == ...
        curr_learned_day,1,1,n_rxn_altsample)),[],n_rxn_altsample);
    curr_alt_prct_altdiff = curr_alt_prct-nanmean(curr_alt_prct,2);
    rxn_alt_prct_altdiff_ci(learned_day_grp==curr_learned_day,:) = ...
        prctile(nanmean(curr_alt_prct_altdiff,1),[1,99]);
end

% (stair plot: extend last step by 1/2 day to draw errorbar there)
figure; hold on
xline(0,'linestyle','--');
AP_errorfill(learned_day_grp(plot_learned),nanmean(rxn_alt_prct_altdiff_ci(plot_learned,:),2), ...
    rxn_alt_prct_altdiff_ci(plot_learned,:),[0.5,0.5,0.5],[],false);
stairs([learned_day_grp(plot_learned);learned_day_grp(find(plot_learned,1,'last'))+0.5], ...
    [rxn_altdiff_learn_mean(plot_learned);rxn_altdiff_learn_mean(find(plot_learned,1,'last'))], ...
    'k','linewidth',2);
errorbar(learned_day_grp(plot_learned)+0.5,rxn_altdiff_learn_mean(plot_learned), ...
    rxn_altdiff_learn_sem(plot_learned),'k','linewidth',2,'linestyle','none','CapSize',0);
% % (plot all mice: NaN-out day without minimum n)
% plot_rxn_stairs = rxn_measured_prct_altdiff;
% plot_rxn_stairs(~ismember(learned_day_x,learned_day_grp(plot_learned))) = NaN;
% stairs(learned_day_x,plot_rxn_stairs,'color','k');
xlabel('Learned day');
ylabel('Fast reaction times (% meas-null)')
axis tight;
xlim(xlim+[-0.5,0.5])

%%%%%%%%% TESTING (above with CIs)
[rxn_altdiff_learn_mean,rxn_altdiff_learn_sem,learned_day_grp,learned_day_n] = ...
    grpstats(rxn_measured_prct_altdiff(:),learned_day_x(:), ...
    {'nanmedian','sem','gname','numel'});
learned_day_grp = cellfun(@str2num,learned_day_grp);
plot_learned = learned_day_n >= min_n;

figure; hold on
xline(0,'linestyle','--');
AP_errorfill(learned_day_grp(plot_learned),nanmean(rxn_alt_prct_altdiff_ci(plot_learned,:),2), ...
    rxn_alt_prct_altdiff_ci(plot_learned,:),[0.5,0.5,0.5],[],false);
plot(learned_day_grp(plot_learned),rxn_altdiff_learn_mean(plot_learned),'k','linewidth',2);

[curr_cil,curr_ciu] = ...
    grpstats(rxn_measured_prct_altdiff(:),learned_day_x(:), ...
    {@(x) prctile(x,[25]),@(x) prctile(x,[75])});

plot(learned_day_grp(plot_learned),[curr_cil(plot_learned,:),curr_ciu(plot_learned,:)],'k','linewidth',1);

xlabel('Learned day');
ylabel('Fast reaction times (% meas-null)')
axis tight;
xlim(xlim+[-0.5,0.5])
%%%%%%%%% TESTING



% Plot fraction of reaction times within window: daysplit
rxn_measured_prct_daysplit = accumarray(trial_split_idx_cat, ...
    rxn_measured_cat >= rxn_window(1) & ...
    rxn_measured_cat <= rxn_window(2), ...
    [n_daysplit,max_days,length(bhv)],@(x) nanmean(x)*100,NaN);
rxn_alt_prct_daysplit = cell2mat(permute(arrayfun(@(x) ...
    accumarray(trial_split_idx_cat, ...
    rxn_alt_cat(:,x) >= rxn_window(1) & ...
    rxn_alt_cat(:,x) <= rxn_window(2), ...
    [n_daysplit,max_days,length(bhv)],@(x) nanmean(x)*100,NaN), ...
    1:n_rxn_altsample,'uni',false),[1,3,4,2]));

% Put NaNs between days to plot with gaps
rxn_measured_prct_long = reshape(padarray(rxn_measured_prct_daysplit,[1,0,0],NaN,'post'),[],length(animals));
rxn_alt_prct_long = reshape(padarray(rxn_alt_prct_daysplit,[1,0,0],NaN,'post'),[],length(animals),n_rxn_altsample);

% Plot relative to training day (daysplit)
figure; hold on;
daysplit_x = [1:size(rxn_measured_prct_long,1)]/(n_daysplit+1);

rxn_alt_prct_long_ci = ...
    permute(prctile(nanmean(rxn_alt_prct_long,2),[5,95],3),[1,3,2]);
AP_errorfill(daysplit_x,nanmean(nanmean(rxn_alt_prct_long,2),3), ...
    rxn_alt_prct_long_ci,[0.5,0.5,0.5],[],false);

errorbar(daysplit_x,nanmean(rxn_measured_prct_long,2), ...
    AP_sem(rxn_measured_prct_long,2),'k','linewidth',2,'CapSize',0);

xlabel('Training day');
ylabel(sprintf('Reaction times %.g-%.g (%%)',rxn_window(1),rxn_window(2)));
axis tight;
xlim(xlim + [-0.5,0.5]);
ylim([0,100]);

% Plot relative to learned day (daysplit)
learned_daysplit_x = cell2mat(cellfun(@(x) x+(0:n_daysplit)'/(n_daysplit+1), ...
    num2cell(learned_day_x),'uni',false));

[rxn_learn_mean_daysplit,rxn_learn_sem_daysplit,learned_day_grp_daysplit,learned_day_n_daysplit] = ...
    grpstats(rxn_measured_prct_long(:),learned_daysplit_x(:), ...
    {'nanmean','sem','gname','numel'});
rxn_alt_learn_mean_daysplit = ...
    grpstats(reshape(rxn_alt_prct_long,[],n_rxn_altsample),learned_daysplit_x(:), ...
    {'nanmean'});

learned_day_grp_daysplit = cellfun(@str2num,learned_day_grp_daysplit);
plot_learned = learned_day_n_daysplit >= min_n | isnan(rxn_learn_mean_daysplit);

figure; hold on;

rxn_alt_learn_ci = prctile(rxn_alt_learn_mean_daysplit,[5,95],2);
p1 = AP_errorfill(learned_day_grp_daysplit(plot_learned), ...
    nanmean(rxn_alt_learn_mean_daysplit(plot_learned,:),2), ...
    rxn_alt_learn_ci(plot_learned,:),[0.5,0.5,0.5],[],false);

% plot(learned_daysplit_x,rxn_measured_prct_long,'color',[0.5,0.5,0.5])
p2 = errorbar(learned_day_grp_daysplit(plot_learned),rxn_learn_mean_daysplit(plot_learned), ...
    rxn_learn_sem_daysplit(plot_learned),'k','linewidth',2,'CapSize',0);
ylabel(sprintf('Reaction times %.g-%.g (%%)',rxn_window(1),rxn_window(2)));
xlabel('Learned day');
axis tight;
xlim(xlim + [-0.5,0.5]);
ylim([0,100]);
line([0,0],ylim,'color','k','linestyle','--');
legend([p1(1),p2(1)],{'Null','Measured'});


% Plot relative to learned day (null difference: daysplit)
rxn_altdiff_learn_mean_daysplit = rxn_learn_mean_daysplit - ...
    nanmean(rxn_alt_learn_mean_daysplit,2);
rxn_alt_altdiff_learn_ci = prctile( ...
    rxn_alt_learn_mean_daysplit-nanmean(rxn_alt_learn_mean_daysplit,2), ...
    [5,95],2);

figure; hold on;

p1 = AP_errorfill(learned_day_grp_daysplit(plot_learned), ...
    nanmean(rxn_alt_learn_mean_daysplit(plot_learned,:),2), ...
    rxn_alt_altdiff_learn_ci(plot_learned,:),[0.5,0.5,0.5],[],false);

p2 = errorbar(learned_day_grp_daysplit(plot_learned), ...
    rxn_altdiff_learn_mean_daysplit(plot_learned), ...
    rxn_learn_sem_daysplit(plot_learned),'k','linewidth',2,'CapSize',0);
ylabel(sprintf('Reaction times %.g-%.g (meas-null%%)',rxn_window(1),rxn_window(2)));
xlabel('Learned day');
axis tight;
xlim(xlim + [-0.5,0.5]);
xline(0,'--');
legend([p1(1),p2(1)],{'Null','Measured'});


%%%%%% Median, daysplit

rxn_measured_med_daysplit = accumarray(trial_split_idx_cat, ...
    rxn_measured_cat.*AP_nanout(rxn_measured_cat < 0.1), ...
    [n_daysplit,max_days,length(bhv)],@(x) nanmedian(x),NaN);

rxn_alt_med_daysplit = cell2mat(permute(arrayfun(@(x) ...
    accumarray(trial_split_idx_cat, ...
    rxn_alt_cat(:,x).*AP_nanout(rxn_alt_cat(:,x) < 0.1), ...
    [n_daysplit,max_days,length(bhv)],@(x) nanmedian(x),NaN), ...
    1:n_rxn_altsample,'uni',false),[1,3,4,2]));


% Put NaNs between days to plot with gaps
rxn_measured_med_long = reshape(padarray(rxn_measured_med_daysplit,[1,0,0],NaN,'post'),[],length(animals));
rxn_alt_med_long = reshape(padarray(rxn_alt_med_daysplit,[1,0,0],NaN,'post'),[],length(animals),n_rxn_altsample);

% Plot relative to learned day (daysplit)
learned_daysplit_x = cell2mat(cellfun(@(x) x+(0:n_daysplit)'/(n_daysplit+1), ...
    num2cell(learned_day_x),'uni',false));

[rxn_learn_mean_daysplit,rxn_learn_sem_daysplit,learned_day_grp_daysplit,learned_day_n_daysplit] = ...
    grpstats(rxn_measured_med_long(:),learned_daysplit_x(:), ...
    {'nanmean','sem','gname','numel'});
rxn_alt_learn_mean_daysplit = ...
    grpstats(reshape(rxn_alt_med_long,[],n_rxn_altsample),learned_daysplit_x(:), ...
    {'nanmean'});

learned_day_grp_daysplit = cellfun(@str2num,learned_day_grp_daysplit);
plot_learned = learned_day_n_daysplit >= min_n | isnan(rxn_learn_mean_daysplit);

figure; hold on;set(gca,'YScale','log');

rxn_alt_learn_ci = prctile(rxn_alt_learn_mean_daysplit,[1,99],2);
p1 = AP_errorfill(learned_day_grp_daysplit(plot_learned), ...
    nanmean(rxn_alt_learn_mean_daysplit(plot_learned,:),2), ...
    rxn_alt_learn_ci(plot_learned,:),[0.5,0.5,0.5],[],false);

p2 = errorbar(learned_day_grp_daysplit(plot_learned),rxn_learn_mean_daysplit(plot_learned), ...
    rxn_learn_sem_daysplit(plot_learned),'k','linewidth',2,'CapSize',0);
ylabel('Median reaction time')
xlabel('Learned day');
axis tight;
xlim(xlim + [-0.5,0.5]);
line([0,0],ylim,'color','k','linestyle','--');
legend([p1(1),p2(1)],{'Null','Measured'});


% Plot relative to learned day (null difference: daysplit)
rxn_altdiff_learn_mean_daysplit = rxn_learn_mean_daysplit - ...
    nanmean(rxn_alt_learn_mean_daysplit,2);
rxn_alt_altdiff_learn_ci = prctile( ...
    rxn_alt_learn_mean_daysplit-nanmean(rxn_alt_learn_mean_daysplit,2), ...
    [5,95],2);

figure; hold on;

p1 = AP_errorfill(learned_day_grp_daysplit(plot_learned), ...
    nanmean(rxn_alt_learn_mean_daysplit(plot_learned,:),2), ...
    rxn_alt_altdiff_learn_ci(plot_learned,:),[0.5,0.5,0.5],[],false);

plot(learned_day_grp_daysplit(plot_learned), ...
    rxn_altdiff_learn_mean_daysplit(plot_learned),'k','linewidth',2);

p2 = errorbar(learned_day_grp_daysplit(plot_learned), ...
    rxn_altdiff_learn_mean_daysplit(plot_learned), ...
    rxn_learn_sem_daysplit(plot_learned),'k','linewidth',2,'CapSize',0);
ylabel(sprintf('Reaction times %.g-%.g (meas-null%%)',rxn_window(1),rxn_window(2)));
xlabel('Learned day');
axis tight;
xlim(xlim + [-0.5,0.5]);
xline(0,'--');
legend([p1(1),p2(1)],{'Null','Measured'});




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

% Grab day index of [last pre-muscimol, V1 muscimol, washout] and data to plot
n_conditions = 3;
condition_labels = {'Pre-muscimol','V1 muscimol','Washout'};

muscimol_v1_days = nan(length(bhv),n_conditions);
muscimol_v1_retinotopy = cell(size(muscimol_v1_days));
muscimol_v1_stim_surround_wheel = cell(size(muscimol_v1_days));
muscimol_v1_wheel_mm = nan(size(muscimol_v1_days));
for curr_animal = 1:length(bhv)
    
    muscimol_animal_idx = ismember({muscimol.animal},bhv(curr_animal).animal);
    if ~any(muscimol_animal_idx)
        continue
    end
    
    % (find last day before muscimol)
    curr_premuscimol_dayidx = find(datenum(bhv(curr_animal).day) < ...
        datenum(muscimol(muscimol_animal_idx).day{1}),1,'last');
    
    % (find last V1 muscimol day)
    curr_v1_muscimol = find(strcmp(lower( ...
        muscimol(muscimol_animal_idx).area),'v1'),1,'last');
    curr_v1_muscimol_dayidx = find(strcmp(bhv(curr_animal).day, ...
        muscimol(muscimol_animal_idx).day(curr_v1_muscimol)));
    
    % (find first washout after V1 muscimol)
    curr_v1_washout = curr_v1_muscimol + ...
        find(strcmp(lower( ...
        muscimol(muscimol_animal_idx).area(curr_v1_muscimol+1:end)), ...
        'washout'),1,'first');
    curr_v1_washout_dayidx = find(strcmp(bhv(curr_animal).day, ...
        muscimol(muscimol_animal_idx).day(curr_v1_washout)));
    
    % (combine conditions)
    condition_dayidx = [curr_premuscimol_dayidx,curr_v1_muscimol_dayidx,curr_v1_washout_dayidx];
    
    % Grab days
    muscimol_v1_days(curr_animal,:) = condition_dayidx;
    
    % Grab retinotopy
    % (pre-muscimol)
    animal = muscimol(muscimol_animal_idx).animal;
    retinotopy_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment\retinotopy';
    load([retinotopy_path filesep animal '_retinotopy'])
    curr_premuscimol_vfs_dayidx = find(datenum({retinotopy.day}) < ...
        datenum(muscimol(muscimol_animal_idx).day{1}),1,'last');
    curr_premuscimol_vfs = AP_align_widefield( ...
        retinotopy(curr_premuscimol_vfs_dayidx).vfs,animal, ...
        retinotopy(curr_premuscimol_vfs_dayidx).day);
    % (add post-muscimol)
    muscimol_v1_retinotopy(curr_animal,:) = ...
        [{curr_premuscimol_vfs} , ...
        muscimol(muscimol_animal_idx).vfs([curr_v1_muscimol,curr_v1_washout])];
    
    % Grab stim-aligned movement
    % (NaN-out times during quiescence)
    stim_surround_t = bhv(1).stim_surround_t;
    curr_stim_surround_wheel = ...
        bhv(curr_animal).stim_surround_wheel(condition_dayidx);
    curr_quiescence_t = ...
        bhv(curr_animal).quiescence_t(condition_dayidx);
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
        bhv(curr_animal).wheel_mm(condition_dayidx)./ ...
        bhv(curr_animal).session_duration(condition_dayidx);
    
end

% Plot retinotopy difference
figure;
h = tiledlayout(1,n_conditions,'TileSpacing','compact','padding','compact');
for curr_cond = 1:n_conditions
    nexttile;
    imagesc(nanmean(cat(3,muscimol_v1_retinotopy{:,curr_cond}),3));
    axis image off;
    caxis([-1,1]);
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title(condition_labels{curr_cond})
end
colormap(brewermap([],'*RdBu'));

% (set 'use_days' to copy code from above)
use_days = mat2cell(muscimol_v1_days,ones(length(muscimol_v1_animals),1),n_conditions)';

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
    [5,95],1),1:n_conditions,'uni',false);

figure;
h = tiledlayout(1,n_conditions,'TileSpacing','compact','padding','compact');
for curr_cond = 1:n_conditions
    nexttile; hold on
    AP_errorfill(rxn_bin_centers,nanmean(nanmean(cat(3,rxn_alt_hist{1,:}),3),1), ...
        rxn_alt_hist_ci{curr_cond}','r',[],false);

    AP_errorfill(rxn_bin_centers, ...
        nanmean(cat(1,rxn_measured_hist{curr_cond,:}),1)', ...
        AP_sem(cat(1,rxn_measured_hist{curr_cond,:}),1)','k');

    legend({'Null','Measured'});
    xlabel('Reaction time');
    ylabel('Probability')
    title(condition_labels{curr_cond});
end
linkaxes(allchild(h),'xy');

% Plot stim-aligned wheel movement and total velocity
figure;
h = tiledlayout(1,2);

nexttile;
errorbar(nanmean(muscimol_v1_wheel_mm,1), ...
    AP_sem(muscimol_v1_wheel_mm,1),'k','linewidth',2);
xlim(xlim+[-0.5,0.5]);
ylim([0,max(ylim)])
set(gca,'XTick',1:n_conditions,'XTickLabel',condition_labels);
ylabel('Wheel mm/min');

muscimol_stim_wheel_avg = ...
    arrayfun(@(cond) cell2mat(cellfun(@(x) ...
    nanmean(x,1),muscimol_v1_stim_surround_wheel(:,cond),'uni',false)), ...
    1:n_conditions,'uni',false);
nexttile;
p = AP_errorfill(bhv(1).stim_surround_t', ...
    cell2mat(cellfun(@(x) nanmean(x,1)',muscimol_stim_wheel_avg,'uni',false)), ...
    cell2mat(cellfun(@(x) AP_sem(x,1)',muscimol_stim_wheel_avg,'uni',false)));
xline(0,'linestyle','--');
xlabel('Time from stim (s)');
ylabel('Probability of movement');
legend(p,condition_labels);

% Plot fraction of reaction times within window 
rxn_window = bhv(1).learned_days_rxn_window;

rxn_measured_prct = cellfun(@(x) ...
    100*nanmean(x >= rxn_window(1) & x <= rxn_window(2)), ...
    rxn_measured_cat);
rxn_alt_prct = cellfun(@(x) ...
    100*nanmean(x >= rxn_window(1) & x <= rxn_window(2),1), ...
    rxn_alt_cat,'uni',false);

figure; hold on
rxn_alt_frac_ci = prctile(nanmean(cell2mat(permute(rxn_alt_prct,[1,3,2])),3),[5,95],2);
AP_errorfill([],nanmean(rxn_alt_frac_ci,2),rxn_alt_frac_ci,[0.5,0.5,0.5],[],false);
errorbar(mean(rxn_measured_prct,2),AP_sem(rxn_measured_prct,2),'k','linewidth',2)
xlim([0,n_conditions] + 0.5);
ylim([0,100]);
set(gca,'XTick',1:n_conditions,'XTickLabel',condition_labels);
ylabel(sprintf('Reaction times: %d-%dms (%%)',1000*rxn_window(1),1000*rxn_window(2)));

% Plot fraction of reaction times within window (null difference)
rxn_alt_prct_mean = cellfun(@nanmean,rxn_alt_prct);
rxn_measured_altdiff_prct = rxn_measured_prct - rxn_alt_prct_mean;

figure; hold on
rxn_alt_altdiff_frac_ci = prctile(nanmean(cell2mat(permute( ...
    cellfun(@(x) x-nanmean(x),rxn_alt_prct,'uni',false),[1,3,2])),3),[5,95],2);
AP_errorfill([],nanmean(rxn_alt_altdiff_frac_ci,2),rxn_alt_altdiff_frac_ci,'r',[],false);
plot(rxn_measured_altdiff_prct,'color',[0.5,0.5,0.5]);
plot(nanmean(rxn_measured_altdiff_prct,2),'k','linewidth',2);
xlim([0,n_conditions] + 0.5);
ylim([-10,100]);
set(gca,'XTick',1:n_conditions,'XTickLabel',condition_labels);
ylabel(sprintf('Fast reaction times: %d-%dms (%%meas-null)',1000*rxn_window(1),1000*rxn_window(2)));

% UNUSED: MAD
% % Get reaction times MAD
% rxn_measured_mad = cellfun(@(x) mad(x.*AP_nanout(x<0.1),1,1),rxn_measured_cat);
% rxn_alt_mad = cellfun(@(x) mad(x.*AP_nanout(x<0.1),1,1),rxn_alt_cat,'uni',false);
% 
% rxn_measured_mad_altdiff = rxn_measured_mad - cellfun(@nanmean,rxn_alt_mad);
% rxn_alt_mad_altdiff = cell2mat(cellfun(@(x) reshape(x-nanmean(x),1,1,[]),rxn_alt_mad,'uni',false));
% 
% figure; hold on
% rxn_alt_mad_altdiff_ci = permute(prctile(nanmean(rxn_alt_mad_altdiff,2),[5,95],3),[1,3,2]);
% AP_errorfill([],nanmean(rxn_alt_mad_altdiff_ci,2),rxn_alt_mad_altdiff_ci,[0.5,0.5,0.5],[],false);
% plot(rxn_measured_mad_altdiff,'k');
% xlim([0,n_conditions] + 0.5);
% set(gca,'XTick',1:n_conditions,'XTickLabel',condition_labels);
% ylabel('Reaction time MAD (s)');

% Get reaction times median (exclude too-fast <0.1)
rxn_measured_med = cellfun(@(x) nanmedian(x.*AP_nanout(x<0.1),1),rxn_measured_cat);
rxn_alt_med = cell2mat(cellfun(@(x) ...
    reshape(nanmedian(x.*AP_nanout(x<0.1),1),1,1,[]),rxn_alt_cat,'uni',false));

rxn_measured_med_altdiff = rxn_measured_med - nanmean(rxn_alt_med,3);
rxn_alt_med_altdiff = rxn_alt_med - nanmean(rxn_alt_med,3);

figure; 

subplot(1,2,1,'YScale','log'); hold on;
rxn_alt_med_ci = permute(prctile(nanmean(rxn_alt_med,2),[5,95],3),[1,3,2]);
AP_errorfill([],nanmean(rxn_alt_med_ci,2),rxn_alt_med_ci,[0.5,0.5,0.5],[],false);
errorbar(nanmean(rxn_measured_med,2),AP_sem(rxn_measured_med,2),'k','linewidth',2);
xlim([0,n_conditions] + 0.5);
set(gca,'XTick',1:n_conditions,'XTickLabel',condition_labels);
ylabel('Median reaction time (s)');

subplot(1,2,2); hold on
rxn_alt_med_altdiff_ci = permute(prctile(nanmean(rxn_alt_med_altdiff,2),[5,95],3),[1,3,2]);
AP_errorfill([],nanmean(rxn_alt_med_altdiff_ci,2),rxn_alt_med_altdiff_ci,[0.5,0.5,0.5],[],false);
errorbar(nanmean(rxn_measured_med_altdiff,2),AP_sem(rxn_measured_med_altdiff,2),'k','linewidth',2);
xlim([0,n_conditions] + 0.5);
set(gca,'XTick',1:n_conditions,'XTickLabel',condition_labels);
ylabel('Median reaction time (meas-null,s)');

% Get fraction of muscimol days with significant fast reaction times
muscimol_learned_days = ...
    cell2mat(cellfun(@(ld,use_days) ld(use_days), ...
    {bhv.learned_days},use_days,'uni',false))';
figure;
imagesc(muscimol_learned_days);
colormap(gray);
set(gca,'XTick',1:n_conditions,'XTickLabel',condition_labels, ...
    'YTick',1:length({bhv.animal}),'YTickLabel',{bhv.animal});
title('% Reaction times sig.');



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

%% Passive - [LOAD DATA]

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_passive_tetO';
AP_load_trials_operant;
n_naive = 3; % (number of naive passive-only days, just hard-coding)
min_n = 4; % (minimum n to plot data)

% Load behavior, exclude animals not in dataset, get learned day
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
bhv_fn = [data_path filesep 'bhv_teto'];
load(bhv_fn);
bhv = bhv(ismember({bhv.animal},animals)); 
learned_day = cellfun(@(x) find(x,1),{bhv.learned_days})';

% Get animal and day index for each trial
trial_animal = cell2mat(arrayfun(@(x) ...
    x*ones(size(vertcat(wheel_all{x}{:}),1),1), ...
    [1:length(wheel_all)]','uni',false));
trial_day = cell2mat(cellfun(@(x) cell2mat(cellfun(@(curr_day,x) ...
    curr_day*ones(size(x,1),1),num2cell(1:length(x))',x,'uni',false)), ...
    wheel_all,'uni',false));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));

% Get "learned day" for each trial
trial_learned_day = cell2mat(cellfun(@(x,ld) cell2mat(cellfun(@(curr_day,x) ...
    curr_day*ones(size(x,1),1),num2cell([1:length(x)]-ld)',x,'uni',false)), ...
    wheel_all,num2cell(learned_day + n_naive),'uni',false));

% Define learning stages
% (learning stages: day 1-3 passive-only, then pre/post learning)
trial_stage = 1*(trial_day <= 3) + ...
    2*(trial_day > 3 & trial_learned_day < 0) + ...
    3*(trial_learned_day >= 0);
n_stages = length(unique(trial_stage));

% Get trials with movement during stim to exclude
quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > 0,2);

% Turn values into IDs for grouping
stim_unique = unique(trial_stim_allcat);
[~,trial_stim_id] = ismember(trial_stim_allcat,stim_unique);

learned_day_unique = unique(trial_learned_day)';
[~,trial_learned_day_id] = ismember(trial_learned_day,learned_day_unique);

% Get average fluorescence by animal/day/stim
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



%% ^^ Passive - stim-aligned averages

% Average V/ROI by learning stage
% (combined naive and pre-learn)
stim_v_avg_stage = cell2mat(permute(cellfun(@(x,ld) ...
    cat(3, ...
    nanmean(x(:,:,1:ld-1,:),3), ...
    nanmean(x(:,:,ld:end,:),3)), ...
    stim_v_avg,num2cell(learned_day+n_naive), ...
    'uni',false),[2,3,4,5,1]));

stim_roi_avg_stage = cell2mat(permute(cellfun(@(x,ld) ...
    cat(3, ...
    nanmean(x(:,:,1:ld-1,:),3), ...
    nanmean(x(:,:,ld:end,:),3)), ...
    stim_roi_avg,num2cell(learned_day+n_naive), ...
    'uni',false),[2,3,4,5,1]));

n_stages = size(stim_v_avg_stage,3);

% Get pixels and pixel timemax by stage
stim_px_avg_stage = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),nanmean(stim_v_avg_stage,5));
AP_imscroll(reshape(permute(stim_px_avg_stage, ...
    [1,4,2,5,3]),size(U_master,1)*n_stages,size(U_master,2)*3,length(t)),t);
axis image off;
colormap(AP_colormap('KWG'));
caxis([-max(abs(caxis)),max(abs(caxis))]);

use_t = t >= 0 & t <= 0.5;
stim_px_avg_stage_tmax = ...
    squeeze(max(stim_px_avg_stage(:,:,use_t,:,:),[],3));

% Plot pixel timeavg
figure;
h = tiledlayout(n_stages,length(stim_unique));
c = (max(stim_px_avg_stage_tmax(:)).*[-1,1])*0.5;
for curr_stage = 1:n_stages
    for curr_stim = 1:length(stim_unique)

        curr_px = stim_px_avg_stage_tmax(:,:,curr_stage,curr_stim);
        curr_hemidiff =  curr_px - AP_reflect_widefield(curr_px);
        x_midpoint = size(U_master,2)/2;

        nexttile;
        imagesc(curr_px);
        axis image off;
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        colormap(gca,AP_colormap('KWG',[],1.5));
        caxis(c)
 
    end
end
linkaxes(allchild(h),'xy');

% Plot pixel timeavg (full + hemidiff)
figure;
h = tiledlayout(n_stages,6,'TileSpacing','compact','padding','compact');
c = (max(stim_px_avg_stage_tmax(:)).*[-1,1])*0.3;
for curr_stage = 1:n_stages
    for curr_stim = 1:length(stim_unique)

        curr_px = stim_px_avg_stage_tmax(:,:,curr_stage,curr_stim);
        curr_hemidiff =  curr_px - AP_reflect_widefield(curr_px);
        x_midpoint = size(U_master,2)/2;

        nexttile;
        imagesc(curr_px);
        axis image off;
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        colormap(gca,AP_colormap('KWG',[],1.5));
        caxis(c)
        
        nexttile;
        imagesc(curr_hemidiff(:,1:x_midpoint));
        axis image off;
        AP_reference_outline('ccf_aligned_lefthemi',[0.5,0.5,0.5]);
        colormap(gca,AP_colormap('BWR',[],1.5));
        caxis(c)
    end
end
linkaxes(allchild(h),'xy');

% (save the avg pixel image for plotting ephys below)
save_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
stim_px_fn = fullfile(save_data_path,'stim_px');
stim_px = squeeze(stim_px_avg_stage_tmax(:,:,end,3));
save(stim_px_fn,'stim_px');


% Plot ROIs by stage
figure;
plot_rois = [6] + [0,size(wf_roi,1)];
stage_col = repmat(linspace(0.3,0,n_stages)',[1,3]);
h = tiledlayout(length(plot_rois),3,'TileSpacing','compact','padding','compact');
for curr_roi = plot_rois
  for curr_stim = stim_unique'
    nexttile;
    AP_errorfill(t, ...
        squeeze(nanmean(stim_roi_avg_stage(curr_roi,:,:,stim_unique == curr_stim,:),5)), ...
        squeeze(AP_sem(stim_roi_avg_stage(curr_roi,:,:,stim_unique == curr_stim,:),5)), ...
        stage_col(:,:));
    xlabel('Time from stim (s)');
    ylabel(sprintf('%s \\DeltaF/F_0',wf_roi(curr_roi).area));
    axis tight;xlim(xlim+[-0.1,0.1])
  end
end
% Link all axes
linkaxes(allchild(h),'xy');
% (shade stim area and put in back)
arrayfun(@(x) patch(x,[0,0.5,0.5,0], ...
    reshape(repmat(ylim(x),2,1),[],1),[1,1,0.8], ...
    'linestyle','none'),allchild(h));
arrayfun(@(x) set(x,'children',circshift(get(x,'children'),-1)),allchild(h));



%% ^^ Passive - ROIs across days by animal

% Plot time-max within ROIs across days
use_t = t > 0 & t <= 0.2;

stim_roi_tmax = cellfun(@(x) ...
    permute(max(x(:,use_t,:,:),[],2),[1,3,4,2]), ...
    stim_roi_avg,'uni',false);

% roi_hemiflip = circshift(1:n_rois,n_rois/2);
% stim_roi_hemidiff_avg = cellfun(@(x) x-x(roi_hemiflip,:,:,:),stim_roi_avg,'uni',false);
% stim_roi_tmax = cellfun(@(x) ...
%     permute(max(x(:,use_t,:,:),[],2),[1,3,4,2]), ...
%     stim_roi_hemidiff_avg,'uni',false);

plot_roi = 6;

figure;
h = tiledlayout('flow');
for curr_animal = 1:length(animals)
    nexttile; hold on;
    stairs((1:size(stim_roi_tmax{curr_animal},2))-n_naive, ...
        stim_roi_tmax{curr_animal}(plot_roi,:,stim_unique == 1),'k','linewidth',2);
    xline(1);
    xline(learned_day(curr_animal),'color','r','linewidth',2);
    ylabel(wf_roi(plot_roi).area(1:end-2));
    title(animals{curr_animal});
end


figure;hold on;
for curr_animal = 1:length(animals)
    plot((1:size(stim_roi_tmax{curr_animal},2))-n_naive - learned_day(curr_animal), ...
        stim_roi_tmax{curr_animal}(plot_roi,:,stim_unique == 1));
    xline(0);
    ylabel(wf_roi(plot_roi).area(1:end-2));
end


% Plot stairs?

a = cellfun(@(x) squeeze(x(plot_roi,:,3))',stim_roi_tmax,'uni',false);
b = AP_padcatcell(a);

learned_day_x = ((1:size(b,1))-3)' - learned_day';

[m,e,n] = grpstats(b(:),learned_day_x(:),{'mean','sem','numel'});
plot_ld = n >= min_n;

% (stair plot: extend last step by 1/2 day to draw errorbar there)
figure; hold on
learned_day_x_unique = unique(learned_day_x(:));

xline(0,'linestyle','--');
stairs(learned_day_x_unique(plot_ld),m(plot_ld),'k','linewidth',2);
errorbar(learned_day_x_unique(plot_ld)+0.5,m(plot_ld),e(plot_ld), ...
    'k','linestyle','none','linewidth',2,'capsize',0);
xlabel('Learned day');
ylabel('\DeltaF/F');
axis tight;
xlim(xlim+[-0.5,0.5])

% (plot stack sorted by learned day)
figure;hold on;
[~,ld_sort] = sort(learned_day);
imagesc(b(:,ld_sort)');
colormap(flipud(gray))
plot(learned_day(ld_sort),1:length(animals),'r','linewidth',2);

%%% TEMP: put into grid to grab data for plotting over 

x = nan(length(learned_day_unique),length(animals));
for curr_animal = 1:length(animals)

    curr_act =   stim_roi_tmax{curr_animal}(plot_roi,:,stim_unique == 1);
    curr_ld = (1:length(curr_act))-(learned_day(curr_animal)+n_naive);

    [~,ld_idx] = ismember(curr_ld,learned_day_unique);
     x(ld_idx,curr_animal) = curr_act;

end

figure;imagesc([],learned_day_unique,x);
xlabel('Animal');
ylabel('Learned day');

%% ^^ Passive - ROI timecourse by learned day

% Get indicies for averaging 
[learned_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);
[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);
[stim_idx,~,~] = ...
    ndgrid(trial_stim_id,1:length(t),1:n_rois);

accum_learnedday_idx = cat(4,learned_day_idx,t_idx,roi_idx,stim_idx,animal_idx);

use_trials = quiescent_trials;

roi_learnday_avg = accumarray( ...
    reshape(accum_learnedday_idx(use_trials,:,:,:),[],size(accum_learnedday_idx,4)), ...
    reshape(fluor_roi_deconv(use_trials,:,:),[],1), ...
    [max(learned_day_idx(:)),length(t),n_rois,length(stim_unique),length(animals)], ...
    @nanmean,NaN('single'));

use_t = t > 0 & t <= 0.2;
stim_roi_learnday_avg_tmax = squeeze(max(roi_learnday_avg(:,use_t,:,:,:),[],2));


plot_roi = 6;
plot_stim = 1;
curr_tcourse = nanmean(squeeze(roi_learnday_avg(:,:,plot_roi,stim_unique == plot_stim,:)),3);
curr_tmax = squeeze(stim_roi_learnday_avg_tmax(:,plot_roi,stim_unique == plot_stim,:));

figure; 
subplot(1,2,1);
imagesc(t,learned_day_unique,curr_tcourse);
caxis(max(abs(caxis))*[-1,1]);
colormap(AP_colormap('KWG',[],1.5));
yline(-0.5);

subplot(1,2,2);
errorbar(learned_day_unique,nanmean(curr_tmax,2),AP_sem(curr_tmax,2),'k','linewidth',2);
xline(0);










%% ^^ Passive - ROIs tmax by learned day (maybe retire for above)

% Plot time-max within ROIs across days
use_t = t > 0 & t <= 0.2;

stim_roi_tmax = cellfun(@(x) ...
    permute(max(x(:,use_t,:,:),[],2),[1,3,4,2]), ...
    stim_roi_avg,'uni',false);

% roi_hemiflip = circshift(1:n_rois,n_rois/2);
% stim_roi_hemidiff_avg = cellfun(@(x) x-x(roi_hemiflip,:,:,:),stim_roi_avg,'uni',false);
% stim_roi_tmax = cellfun(@(x) ...
%     permute(max(x(:,use_t,:,:),[],2),[1,3,4,2]), ...
%     stim_roi_hemidiff_avg,'uni',false);

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
plot_learned_day = squeeze(stim_roi_tmax_learn_numel(1,:,1)) >= min_n;

figure;
plot_rois = [1,6]; % left, right ROIs will also be plotted
stim_col = [0,0,1;0.5,0.5,0.5;1,0,0];
h = tiledlayout(length(plot_rois),2, ...
    'TileSpacing','compact','padding','compact','TileIndexing','columnmajor');
for curr_hemi = 1:2
    for curr_roi = plot_rois
        
        curr_roi = curr_roi + (curr_hemi-1)*size(wf_roi,1);
        nexttile; hold on;

        % (naive)
        set(gca,'ColorOrder',stim_col);
        errorbar( ...
            repmat((min(learned_day_unique(plot_learned_day))+[-n_naive:-1])',1,length(stim_unique)), ...
            squeeze(stim_roi_tmax_naive_mean(curr_roi,1:n_naive,:)), ...
            squeeze(stim_roi_tmax_naive_sem(curr_roi,1:n_naive,:)), ...
            '.-','MarkerSize',20','linewidth',2,'capsize',0);
        % (training)
        errorbar( ...
            repmat((learned_day_unique(plot_learned_day))',1,length(stim_unique)), ...
            squeeze(stim_roi_tmax_learn_mean(curr_roi,plot_learned_day,:)), ...
            squeeze(stim_roi_tmax_learn_sem(curr_roi,plot_learned_day,:)), ...
            '.-','MarkerSize',20','linewidth',2,'capsize',0);

        title(wf_roi(curr_roi).area);
        xlabel('Learned day');
        ylabel('\DeltaF/F_0');
        xline(0,'linestyle','--');
        axis tight;
        xlim(xlim+[-0.5,0.5]);
    end
end
% (link ROI axes across hemispheres)
ax = reshape(allchild(h),length(plot_rois),2);
for curr_roi = 1:length(plot_rois)
   linkaxes(ax(curr_roi,:),'xy'); 
end

%% ^^ Passive - STATS: ROIs across days

% (to do)

% Get roi averages with shuffled -1/1 stim IDs

% Get average fluorescence by animal/day/stim
stim_roi_avg_lrshuff = cell(length(animals),1);
for curr_animal = 1:length(animals)
    for curr_day = 1:max(trial_day(trial_animal == curr_animal))
        for curr_stim_idx = 1:length(stim_unique)

            % shuffle left/right stim identity within day   
            use_trials_day = find(quiescent_trials & ...
                trial_animal == curr_animal & ...
                trial_day == curr_day);

            curr_trial_stim_shuff = trial_stim_allcat(use_trials_day);
            curr_trial_stim_shuff(ismember(curr_trial_stim_shuff,[-1,1])) = ...
                AP_shake(curr_trial_stim_shuff(ismember(curr_trial_stim_shuff,[-1,1])));

            use_trials = use_trials_day(curr_trial_stim_shuff == stim_unique(curr_stim_idx));

            stim_roi_avg_lrshuff{curr_animal}(:,:,curr_day,curr_stim_idx) = ...
                permute(nanmean(fluor_roi_deconv(use_trials,:,:),1),[3,2,1]);
            
        end
    end
end




% Plot time-max within ROIs across days
use_t = t > 0 & t <= 0.2;

stim_roi_lrshuff_tmax = cellfun(@(x) ...
    permute(max(x(:,use_t,:,:),[],2),[1,3,4,2]), ...
    stim_roi_avg_lrshuff,'uni',false);

stim_roi_lrshuff_tmax_daycat = cat(2,stim_roi_lrshuff_tmax{:});

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

% (get mean/sem by index during training)
stim_roi_lrshuff_tmax_learn_mean = ...
    accumarray([roi_idx(training_data),learned_day_idx(training_data),stim_idx(training_data)], ...
    stim_roi_lrshuff_tmax_daycat(training_data), ...
    [n_rois,max(learned_day_idx(:)),length(stim_unique)],@nanmean,NaN('single'));



roi_lrshuff_stimdiff = stim_roi_lrshuff_tmax_learn_mean(:,:,3) - ...
    stim_roi_lrshuff_tmax_learn_mean(:,:,1);

plot_roi = 6;
% figure; hold on
plot(learned_day_unique,roi_lrshuff_stimdiff(plot_roi,:));



%% ^^ Passive - ROI activity vs frac rxn times [requires behavior load]

% % Get hemidiff activity
% roi_hemiflip = circshift(1:n_rois,n_rois/2);
% stim_roi_hemidiff_avg = cellfun(@(x) x-x(roi_hemiflip,:,:,:),stim_roi_avg,'uni',false);

% Get time-max within ROIs across days
use_t = t >= 0 & t <= 0.2;
stim_roi_tmax = cellfun(@(x) ...
    permute(max(x(:,use_t,:,:),[],2),[1,3,4,2]), ...
    stim_roi_avg,'uni',false);

% Get fraction window reaction times minus null (load behavior above)
rxn_window = bhv(1).learned_days_rxn_window;
rxn_measured_prct_animal = cellfun(@(meas) cellfun(@(meas)...
    nanmean(meas >= rxn_window(1) & meas <= rxn_window(2)),...
    meas),rxn_measured,'uni',false);

rxn_measured_prct_altdiff_animal = cellfun(@(meas,alt) cellfun(@(meas,alt)...
    nanmean(meas >= rxn_window(1) & meas <= rxn_window(2)) - ...
    nanmean(nanmean(alt >= rxn_window(1) & alt <= rxn_window(2),1),2), ...
    meas,alt),rxn_measured,rxn_alt,'uni',false);

use_bhv = rxn_measured_prct_animal;

% Plot activity vs behavior for each animal (exclude naive activity)
plot_rois = [1,7,6];
animal_colors = max(brewermap(length(animals),'Set3')-0.2,0);
figure;
h = tiledlayout(length(plot_rois),1,'TileSpacing','compact','padding','compact');
for curr_roi = plot_rois
    nexttile; hold on;
    
    set(gca,'ColorOrder',animal_colors);
    % (plot all data as dots)
    p1 = cellfun(@(bhv,act) ...
        plot(bhv,act(curr_roi,n_naive+1:end,stim_unique == 1),'.','MarkerSize',20), ...
        use_bhv,stim_roi_tmax);
    
    % (plot average binned by reaction window percent)
    rxn_prct_cat = cell2mat(use_bhv);
    act_cat = cell2mat(cellfun(@(x) x(curr_roi,n_naive+1:end,stim_unique == 1)', ...
        stim_roi_tmax,'uni',false));

    rxn_prct_bin_edges = linspace(prctile(rxn_prct_cat,10),prctile(rxn_prct_cat,90),6);
    rxn_prct_bins = discretize(rxn_prct_cat,rxn_prct_bin_edges);

    rxn_rxnbin_mean = accumarray(rxn_prct_bins(~isnan(rxn_prct_bins)), ...
        rxn_prct_cat(~isnan(rxn_prct_bins)),[],@nanmean);
    act_rxnbin_mean = accumarray(rxn_prct_bins(~isnan(rxn_prct_bins)), ...
        act_cat(~isnan(rxn_prct_bins)),[],@nanmean);
    act_rxnbin_sem = accumarray(rxn_prct_bins(~isnan(rxn_prct_bins)), ...
        act_cat(~isnan(rxn_prct_bins)),[],@AP_sem);

    p2 = errorbar(rxn_rxnbin_mean,act_rxnbin_mean,act_rxnbin_sem,'k','linewidth',2);

    legend([p1(1),p2],{'Animal, Day','Average'},'location','nw')
    xlabel(sprintf('Reaction times %.g-%.g (%%)',rxn_window(1),rxn_window(2)));
    ylabel(sprintf('%s \\DeltaF/F_0',wf_roi(curr_roi).area))
end

%% ^^ Passive - ROI activity vs frac rxn times by learned day [requires behavior load]

% % Get hemidiff activity
% roi_hemiflip = circshift(1:n_rois,n_rois/2);
% stim_roi_hemidiff_avg = cellfun(@(x) x-x(roi_hemiflip,:,:,:),stim_roi_avg,'uni',false);

% Get time-max within ROIs across days
use_t = t >= 0 & t <= 0.2;
stim_roi_tmax = cellfun(@(x) ...
    permute(max(x(:,use_t,:,:),[],2),[1,3,4,2]), ...
    stim_roi_avg,'uni',false);

% Get fraction window reaction times minus null (load behavior above)
rxn_window = bhv(1).learned_days_rxn_window;
rxn_measured_prct_animal = cellfun(@(meas) cellfun(@(meas)...
    nanmean(meas >= rxn_window(1) & meas <= rxn_window(2)),...
    meas),rxn_measured,'uni',false);

rxn_measured_prct_altdiff_animal = cellfun(@(meas,alt) cellfun(@(meas,alt)...
    nanmean(meas >= rxn_window(1) & meas <= rxn_window(2)) - ...
    nanmean(nanmean(alt >= rxn_window(1) & alt <= rxn_window(2),1),2), ...
    meas,alt),rxn_measured,rxn_alt,'uni',false);

use_bhv = rxn_measured_prct_animal;

% Plot activity vs behavior for each animal (exclude naive activity)
plot_rois = [1,7,6];
animal_colors = max(brewermap(length(animals),'Set3')-0.2,0);
figure;
h = tiledlayout(length(plot_rois),1,'TileSpacing','compact','padding','compact');
for curr_roi = plot_rois
    nexttile; hold on;
    
    % (plot all data as dots)
    p1 = cellfun(@(bhv,act,ld) ...
        plot(bhv(1:ld-1),act(curr_roi,n_naive+1:n_naive+(ld-1),stim_unique == 1),'.k','MarkerSize',10), ...
        use_bhv,stim_roi_tmax,num2cell(learned_day));

    p2 = cellfun(@(bhv,act,ld) ...
        plot(bhv(ld:end),act(curr_roi,n_naive+ld:end,stim_unique == 1),'.r','MarkerSize',10), ...
        use_bhv,stim_roi_tmax,num2cell(learned_day));
 
end


%% ^^ Passive - ROI activity vs rxn mad [requires behavior load]

% Get hemidiff activity
roi_hemiflip = circshift(1:n_rois,n_rois/2);
stim_roi_hemidiff_avg = cellfun(@(x) x-x(roi_hemiflip,:,:,:),stim_roi_avg,'uni',false);

% Get time-max within ROIs across days
use_t = t >= 0 & t <= 0.2;
stim_roi_tmax = cellfun(@(x) ...
    permute(max(x(:,use_t,:,:),[],2),[1,3,4,2]), ...
    stim_roi_hemidiff_avg,'uni',false);

% Get reaction time mad minus null (load behavior above)
rxn_mad_altdiff = cellfun(@(meas,null) cellfun(@(meas,null)...
    mad(meas) - nanmean(mad(null,[],1),2),meas,null),rxn_measured,rxn_alt,'uni',false);

% Plot activity vs behavior for each animal (exclude naive activity)
plot_rois = [1,7,6];
animal_colors = max(brewermap(length(animals),'Set3')-0.2,0);
figure;
h = tiledlayout(length(plot_rois),1,'TileSpacing','compact','padding','compact');
for curr_roi = plot_rois
    nexttile; hold on;
    
    set(gca,'ColorOrder',animal_colors);
    % (plot all data as dots)
    p1 = cellfun(@(bhv,act) ...
        plot(bhv,act(curr_roi,n_naive+1:end,stim_unique == 1),'.','MarkerSize',20), ...
        rxn_mad_altdiff,stim_roi_tmax);
    
    % (plot average binned by reaction window percent)
    rxn_mad_altdiff_cat = cell2mat(rxn_mad_altdiff);
    act_cat = cell2mat(cellfun(@(x) ...
        x(curr_roi,n_naive+1:end,stim_unique == 1)', ...
        stim_roi_tmax,'uni',false));

    bhv_bin_edges = linspace(prctile(rxn_mad_altdiff_cat,10),prctile(rxn_mad_altdiff_cat,90),5);
    bhv_bins = discretize(rxn_mad_altdiff_cat,bhv_bin_edges);

    bhv_bhvbin_mean = accumarray(bhv_bins(~isnan(bhv_bins)), ...
        rxn_mad_altdiff_cat(~isnan(bhv_bins)),[],@nanmean);
    act_bhvbin_mean = accumarray(bhv_bins(~isnan(bhv_bins)), ...
        act_cat(~isnan(bhv_bins)),[],@nanmean);
    act_bhvbin_sem = accumarray(bhv_bins(~isnan(bhv_bins)), ...
        act_cat(~isnan(bhv_bins)),[],@AP_sem);

    p2 = errorbar(bhv_bhvbin_mean,act_bhvbin_mean,act_bhvbin_sem,'k','linewidth',2);

    legend([p1(1),p2],{'Animal, Day','Average'},'location','nw')
    xlabel(sprintf('Reaction time m.a.d. (s)',rxn_window(1),rxn_window(2)));
    ylabel(sprintf('%s \\DeltaF/F_0',wf_roi(curr_roi).area))
end


%% ^^ Passive - ROI activity vs rxn median [requires behavior load]

% Get time-max within ROIs across days
use_t = t >= 0 & t <= 0.2;
stim_roi_tmax = cellfun(@(x) ...
    permute(max(x(:,use_t,:,:),[],2),[1,3,4,2]), ...
    stim_roi_avg,'uni',false);

% Get reaction time mad minus null (load behavior above)
rxn_med_altdiff = cellfun(@(meas,null) cellfun(@(meas,null)...
    nanmedian(meas.*AP_nanout(meas < 0.1)) - ...
    nanmean(nanmedian(null.*AP_nanout(null < 0.1),1),2), ...
    meas,null),rxn_measured,rxn_alt,'uni',false);

% Plot activity vs behavior for each animal (exclude naive activity)
plot_rois = [1,7,6];
animal_colors = max(brewermap(length(animals),'Set3')-0.2,0);
figure;
h = tiledlayout(length(plot_rois),1,'TileSpacing','compact','padding','compact');
for curr_roi = plot_rois
    nexttile; hold on;
    
    set(gca,'ColorOrder',animal_colors);
    % (plot all data as dots)
    p1 = cellfun(@(bhv,act) ...
        plot(bhv,act(curr_roi,n_naive+1:end,stim_unique == 1),'.','MarkerSize',20), ...
        rxn_med_altdiff,stim_roi_tmax);
    
    % (plot average binned by reaction window percent)
    rxn_mad_altdiff_cat = cell2mat(rxn_med_altdiff);
    act_cat = cell2mat(cellfun(@(x) ...
        x(curr_roi,n_naive+1:end,stim_unique == 1)', ...
        stim_roi_tmax,'uni',false));

    bhv_bin_edges = linspace(prctile(rxn_mad_altdiff_cat,10),prctile(rxn_mad_altdiff_cat,90),5);
    bhv_bins = discretize(rxn_mad_altdiff_cat,bhv_bin_edges);

    bhv_bhvbin_mean = accumarray(bhv_bins(~isnan(bhv_bins)), ...
        rxn_mad_altdiff_cat(~isnan(bhv_bins)),[],@nanmean);
    act_bhvbin_mean = accumarray(bhv_bins(~isnan(bhv_bins)), ...
        act_cat(~isnan(bhv_bins)),[],@nanmean);
    act_bhvbin_sem = accumarray(bhv_bins(~isnan(bhv_bins)), ...
        act_cat(~isnan(bhv_bins)),[],@AP_sem);

    p2 = errorbar(bhv_bhvbin_mean,act_bhvbin_mean,act_bhvbin_sem,'k','linewidth',2);

    legend([p1(1),p2],{'Animal, Day','Average'},'location','nw')
    xlabel(sprintf('Reaction time median (s)',rxn_window(1),rxn_window(2)));
    ylabel(sprintf('%s \\DeltaF/F_0',wf_roi(curr_roi).area))
end



%% ^^ Passive - ROI activity on day X vs. X+1

% Get hemidiff activity
roi_hemiflip = circshift(1:n_rois,n_rois/2);
stim_roi_hemidiff_avg = cellfun(@(x) x-x(roi_hemiflip,:,:,:),stim_roi_avg,'uni',false);

% Get time-max within ROIs across days
use_t = t >= 0 & t <= 0.2;
stim_roi_tmax = cellfun(@(x) ...
    permute(max(x(:,use_t,:,:),[],2),[1,3,4,2]), ...
    stim_roi_avg,'uni',false);

% % (subtract pre-learning)
% stim_roi_tmax = cellfun(@(x,ld) x-nanmean(x(:,1:length(x) < ld,:),2), ...
%     stim_roi_tmax,num2cell(learned_day + n_naive),'uni',false);

plot_roi = 6;

figure; h = tiledlayout('flow');

nexttile; hold on;
animal_col = max(0,jet(length(animals))-0.2);
p2 = cellfun(@(x,ld,col) plot(x(plot_roi,ld-1,3),x(plot_roi,ld,3),'.', ...
    'color',col,'MarkerSize',30), ...
    stim_roi_tmax,num2cell(learned_day+n_naive), ...
    mat2cell(animal_col,ones(length(animals),1),3));
p1 = cellfun(@(x) plot(x(plot_roi,1:end-1,3),x(plot_roi,2:end,3), ...
    '.','color',[0.5,0.5,0.5],'MarkerSize',15),stim_roi_tmax);
line(xlim,xlim,'color','k');
xlabel('\DeltaF/F_0 day X');
ylabel('\DeltaF/F_0 day X+1');
axis tight;
axis equal;
legend([p1(1),p2(1)],{'All days','Learned day'},'location','nw');
title('All days');


nexttile; hold on;
animal_col = max(0,jet(length(animals))-0.2);
p2 = cellfun(@(x,ld,col) plot(x(plot_roi,ld-1,3),x(plot_roi,ld,3),'.', ...
    'color',col,'MarkerSize',30), ...
    stim_roi_tmax,num2cell(learned_day+n_naive), ...
    mat2cell(animal_col,ones(length(animals),1),3));
p1 = cellfun(@(x,ld) plot(x(plot_roi,1:ld-1,3),x(plot_roi,2:ld,3), ...
    '.','color',[0.5,0.5,0.5],'MarkerSize',15), ...
    stim_roi_tmax,num2cell(learned_day+n_naive));
line(xlim,xlim,'color','k');
xlabel('\DeltaF/F_0 day X');
ylabel('\DeltaF/F_0 day X+1');
axis tight;
axis equal;
legend([p1(1),p2(1)],{'All days','Learned day'},'location','nw');
title('To learned day');


nexttile; hold on;
animal_col = max(0,jet(length(animals))-0.2);
p2 = cellfun(@(x,ld) plot(x(plot_roi,1:ld-2,3),x(plot_roi,2:ld-1,3), ...
    '.','color',[0.5,0.5,0.5],'MarkerSize',15), ...
    stim_roi_tmax,num2cell(learned_day+n_naive));
p3 = cellfun(@(x,ld) plot(x(plot_roi,ld:end-1,3),x(plot_roi,ld+1:end,3), ...
    '.','color','k','MarkerSize',15), ...
    stim_roi_tmax,num2cell(learned_day+n_naive));
p1 = cellfun(@(x,ld,col) plot(x(plot_roi,ld-1,3),x(plot_roi,ld,3),'.', ...
    'color','r','MarkerSize',15), ...
    stim_roi_tmax,num2cell(learned_day+n_naive), ...
    mat2cell(animal_col,ones(length(animals),1),3));
line(xlim,xlim,'color','k');
xlabel('\DeltaF/F_0 day X');
ylabel('\DeltaF/F_0 day X+1');
axis tight;
axis equal;
legend([p1(1),p2(1),p3(1)],{'Learned','Pre-learned','Post-learned'},'location','nw');
title('All days grouped');

%% ^^ Passive - (STATS - Kenneth) test fluorescence is related to learning day

% Get hemidiff activity
roi_hemiflip = circshift(1:n_rois,n_rois/2);
stim_roi_hemidiff_avg = cellfun(@(x) x-x(roi_hemiflip,:,:,:),stim_roi_avg,'uni',false);

% Get time-max within ROIs across days
use_t = t >= 0 & t <= 0.2;
stim_roi_tmax = cellfun(@(x) ...
    permute(max(x(:,use_t,:,:),[],2),[1,3,4,2]), ...
    stim_roi_hemidiff_avg,'uni',false);


curr_roi = 6;
a = cellfun(@(x) squeeze(x(curr_roi,:,3))',stim_roi_tmax,'uni',false);
b = cell2mat(cellfun(@(x,ld) accumarray((1:length(x) >= ld)'+1,x,[],@nanmean), ...
    a,num2cell(learned_day + n_naive),'uni',false)');

figure;
errorbar(nanmean(b,2),AP_sem(b,2))
xlim(xlim + [-0.1,0.1]);


n_shuff = 10000;
c = nan(2,length(animals),n_shuff);
for curr_shuff = 1:n_shuff
    c(:,:,curr_shuff) = ...
        cell2mat(cellfun(@(x,ld) accumarray((1:length(x) >= ld)'+1,x,[],@nanmean), ...
        a,AP_shake(num2cell(learned_day + n_naive)),'uni',false)');
    AP_print_progress_fraction(curr_shuff,n_shuff);
end

m = reshape(cat(3,nanmean(diff(b,[],1),2),nanmean(diff(c,[],1),2)),[],1);
m_rank = tiedrank(m);
m_p = m_rank(1)./(n_shuff+1);

figure;
subplot(1,2,1);
histogram(learned_day)
xlabel('Learned day');
ylabel('# animals');

subplot(1,2,2);
histogram(m(2:end),'EdgeColor','none');
xline(m(1),'color','r','linewidth',2);
xlabel('Mean \DeltaF/F post-pre learning')
ylabel('Frequency (shuffle)');
title(sprintf('p = %.2f',m_p));
legend({'Shuffle','Measured'})



%% ^^ Passive - whisker/pupil

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
use_t = t > 0 & t <= 0.2;
whisker_stim_tmax = squeeze(max(whisker_stim_avg(:,use_t,:,:),[],2));
plot_days = min(sum(~isnan(whisker_stim_tmax),3)) >= min_n;
figure; hold on
set(gca,'ColorOrder',stim_col);
errorbar(repmat(learned_day_unique(plot_days)',1,length(stim_unique)), ...
    nanmean(whisker_stim_tmax(:,plot_days,:),3)', ...
    AP_sem(whisker_stim_tmax(:,plot_days,:),3)','linewidth',2,'CapSize',0);
axis tight;
xlim(xlim+[-0.5,0.5]);
xline(0,'color','k','linestyle','--');
xlabel('Learned day');
ylabel('Whisker movement');

% Plot fluorescence/whisker by discretized whisker movement
plot_rois = [1,6];
min_trials = 10; % (minimum trials to use within recording)

trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);

% (discretize whisker movement within recording - use R stim trials)
use_t = t > 0 & t <= 0.2;
whisker_allcat_tavg = max(whisker_allcat(:,use_t),[],2);
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

%% Task - [LOAD DATA]

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_task_teto';
AP_load_trials_operant;
min_n = 4; % (minimum n to plot data)

% Load behavior, exclude animals not in dataset, get learned day
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
bhv_fn = [data_path filesep 'bhv_teto'];
load(bhv_fn);
bhv = bhv(ismember({bhv.animal},animals)); 
learned_day = cellfun(@(x) find(x,1),{bhv.learned_days})';
learned_day_animal = cellfun(@(ld,n) [1:n]'-(ld), ...
    num2cell(learned_day),num2cell(cellfun(@length,trial_info_all)), ...
    'uni',false);

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

%% ^^ Task - pre/post-learning average movie

%%%%%% Get hemidiff V's

mirror_matrix = reshape(U_master(:,:,1:n_vs),[],n_vs)'* ...
    reshape(AP_reflect_widefield(U_master(:,:,1:n_vs)),[],n_vs);
fluor_allcat_deconv_mirror = reshape(transpose( ...
    mirror_matrix*reshape(fluor_allcat_deconv,[],n_vs)'),size(fluor_allcat_deconv));
fluor_allcat_deconv_hemidiff = fluor_allcat_deconv - fluor_allcat_deconv_mirror;


%%%%%%%


fluor_postlearn_v_avg = nan(n_vs,length(t),length(animals));
fluor_postlearn_v_hemidiff_avg = nan(n_vs,length(t),length(animals));
for curr_animal = 1:length(animals)
    use_trials = trial_animal == curr_animal & trial_learned_day >= 0;
    fluor_postlearn_v_avg(:,:,curr_animal) = permute(nanmean(...
        fluor_allcat_deconv(use_trials,:,:),1),[3,2,1]);
    fluor_postlearn_v_hemidiff_avg(:,:,curr_animal) = permute(nanmean(...
        fluor_allcat_deconv_hemidiff(use_trials,:,:),1),[3,2,1]);
end

fluor_postlearn_px_avg = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    nanmean(fluor_postlearn_v_avg,3));
fluor_postlearn_px_hemidiff_avg = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    nanmean(fluor_postlearn_v_hemidiff_avg,3));

AP_imscroll([fluor_postlearn_px_avg, ...
    fluor_postlearn_px_avg - ...
    AP_reflect_widefield(fluor_postlearn_px_avg), ...
    fluor_postlearn_px_hemidiff_avg],t); 
axis image;
caxis(max(abs(caxis))*[-1,1]);
colormap(brewermap([],'PrGn'));
% colormap(AP_colormap('KWG'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5],[],[size(U_master,1),size(U_master,2),1,3]);

% Plot time average
use_t = t >= 0.05 & t <= 0.2;
fluor_postlearn_px_tavg= nanmean(fluor_postlearn_px_avg(:,:,use_t),3);
fluor_postlearn_px_hemidiff_tavg= nanmean(fluor_postlearn_px_hemidiff_avg(:,:,use_t),3);
c = max(fluor_postlearn_px_tavg(:))*[-1,1];
figure; subplot(1,2,1);
imagesc(fluor_postlearn_px_tavg);
axis image off; caxis(c);
colormap(brewermap([],'PrGn'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
title('Full');

subplot(1,2,2);
ml_midpoint = ceil(size(U_master,2)/2);
imagesc(fluor_postlearn_px_hemidiff_tavg(:,1:ml_midpoint))
axis image off; caxis(c);
colormap(brewermap([],'PrGn'));
AP_reference_outline('ccf_aligned_lefthemi',[0.5,0.5,0.5]);
title('Hemidiff');

linkaxes(get(gcf,'children'),'xy');


%% ^^ Task - event-aligned pixels by pre/post learning

% Set trials to use
learned_day_grp_edges = [-Inf,0,Inf];
trial_learned_day_grp = discretize(trial_learned_day,learned_day_grp_edges);
n_learn_grps = length(learned_day_grp_edges)-1;

% Set alignment shifts
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));
stim_align = ones(size(trial_stim_allcat));
move_align = move_idx - leeway_samples;
outcome_align = outcome_idx - leeway_samples;

% Set windows to average activity
use_align_labels = {'Stim','Move','Outcome'};
use_align = {stim_align,move_align,outcome_align};
plot_t = {0.1,0,0.1};

% Loop through alignments and get pixels
align_px = cellfun(@(x) nan(size(U_master,1),size(U_master,2),length(cell2mat(plot_t))), ...
    use_align,'uni',false);
for curr_align = 1:length(use_align)
    
    % (re-align activity if not all aligned to first frame)
    if all(use_align{curr_align} == 1)
        curr_v_align = fluor_allcat_deconv;
    else
        curr_v_align = nan(size(fluor_allcat_deconv),'single');
        for curr_trial = find(~isnan(use_align{curr_align}))'
            curr_shift_frames = ...
                use_align{curr_align}(curr_trial) + [0:length(t)-1];
            curr_shift_frames_use = curr_shift_frames > 0 & curr_shift_frames <= length(t);

            curr_grab_frames = curr_shift_frames(curr_shift_frames_use);
            curr_fill_frames = find(curr_shift_frames_use,length(t));

            curr_v_align(curr_trial,curr_fill_frames,:) = ...
                fluor_allcat_deconv(curr_trial,curr_grab_frames,:);
        end
    end
    
    curr_v_align_avg = nan(n_vs,length(t),n_learn_grps, ...
        length(animals),class(curr_v_align));
    for curr_animal = 1:length(animals)
        for curr_learned_grp = 1:n_learn_grps
        curr_trials = trial_animal == curr_animal & ...
            trial_learned_day_grp == curr_learned_grp & ...
            trial_outcome_allcat == 1;
        curr_v_align_avg(:,:,curr_learned_grp,curr_animal) = ...
            permute(nanmean(curr_v_align(curr_trials,:,:),1),[3,2,1]);
        end
    end

    % (animal average the last training stage defined above)
    curr_v_align_avg_t = ...
        permute(interp1(t,permute(nanmean(curr_v_align_avg,4),[2,1,3]), ...
        plot_t{curr_align},'previous'),[2,1,3]);
      
    curr_px_align_avg_t = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
        curr_v_align_avg_t);
    
    align_px{curr_align} = curr_px_align_avg_t;

end

% Plot pixels
c = prctile(reshape(cat(3,align_px{:}),[],1),100)*[-1,1]*0.8;
c_hemi = c*0.5;
figure;
h = tiledlayout(n_learn_grps,2*length(cell2mat(plot_t)), ...
    'TileSpacing','compact','padding','compact');
for curr_learn_grp = 1:n_learn_grps
    for curr_align = 1:length(use_align)
        for curr_plot_t = 1:length(plot_t{curr_align})
            
            curr_px = align_px{curr_align}(:,:,curr_plot_t,curr_learn_grp);

            % (full brain)
            nexttile
            imagesc(curr_px);
            AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
            axis image off;
            colormap(gca,AP_colormap('KWG',[],1.5));
            caxis(c);
            title([use_align_labels{curr_align} ': ' num2str(plot_t{curr_align}(curr_plot_t)) ' sec']);

            % (hemi-difference)
            nexttile;
            x_midpoint = size(U_master,2)/2;
            curr_align_px_hemidiff = ...
                curr_px - AP_reflect_widefield(curr_px);
            imagesc(curr_align_px_hemidiff(:,1:x_midpoint));
            AP_reference_outline('ccf_aligned_lefthemi',[0.5,0.5,0.5]);
            axis image off;
            colormap(gca,AP_colormap('BWR',[],1.5));
            caxis(c_hemi);
        end
    end
end
linkaxes(allchild(h),'xy');


%% ^^ Task - event kernels

% Get task>cortex parameters
n_regressors = length(task_regressor_labels);
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

% Average kernels by learning stage
learn_stages = cellfun(@(x) discretize(x,[-Inf,0,Inf]),learned_day_animal,'uni',false);
 
n_learn_stages = max(cell2mat(learn_stages));
fluor_taskpred_k_stage_animals = cell(n_regressors,n_learn_stages,length(animals));
for curr_animal = 1:length(animals)
    curr_k = horzcat(fluor_taskpred_k_all{curr_animal}{:});
    for curr_stage = 1:n_learn_stages
        for curr_reg = 1:n_regressors
            curr_k_avg = nanmean(cat(4,curr_k{curr_reg,learn_stages{curr_animal} == curr_stage}),4);
            fluor_taskpred_k_stage_animals{curr_reg,curr_stage,curr_animal} = ...
                curr_k_avg;
        end
    end
end

fluor_taskpred_k_stage_avg = cell(n_regressors,1);
for curr_stage = 1:n_learn_stages
    for curr_reg = 1:n_regressors
        curr_k = nanmean(cat(4,fluor_taskpred_k_stage_animals{curr_reg,curr_stage,:}),4);
        fluor_taskpred_k_stage_avg{curr_reg}(:,:,:,curr_stage) = ...
            permute(curr_k,[3,2,1]);
    end
end

fluor_taskpred_px_stage_avg = cellfun(@(x) ...
    AP_svdFrameReconstruct(U_master(:,:,1:n_vs),x), ...
    fluor_taskpred_k_stage_avg,'uni',false);

% Plot stim kernel pre/post learning (within time)
stim_regressor = strcmp(task_regressor_labels,'Stim');
    
figure;
h = tiledlayout(1,n_learn_stages);
c = prctile(reshape(fluor_taskpred_px_stage_avg{stim_regressor},[],1),100).*[-1,1]*0.8;
for curr_stage = 1:n_learn_stages
    nexttile;
    imagesc(squeeze(mean( ...
        fluor_taskpred_px_stage_avg{stim_regressor}(:,:,:,:,curr_stage),[],3)));
    title(sprintf('Stage %d',curr_stage));
    axis image off;
    caxis(c);
    colormap(AP_colormap('KWG'));
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    title(h,sprintf('%s kernel', ...
        task_regressor_labels{stim_regressor}));
end


% Get and plot task kernels in ROIs
fluor_taskpred_k_roi_stage_animals = cellfun(@(x) ...
       AP_svd_roi(U_master(:,:,1:n_vs),permute(x,[3,2,1]),[],[],...
       cat(3,wf_roi.mask)),fluor_taskpred_k_stage_animals,'uni',false);

n_subregressors = cellfun(@(x) size(x,3),fluor_taskpred_k_stage_avg);

stim_col = [0.7,0,0];
move_col = [0.6,0,0.6;1,0.6,0];
outcome_col = [0.2,0.8,1;0,0,0];
task_regressor_cols = {stim_col,move_col,outcome_col};

plot_rois = [1,6,7];
figure;
h = tiledlayout(length(plot_rois),sum(n_subregressors),'TileSpacing','compact','padding','compact');
for curr_roi = plot_rois
    for curr_regressor = 1:n_regressors
        for curr_subregressor = 1:n_subregressors(curr_regressor)

            curr_k = cell2mat(permute(cellfun(@(x) x(curr_roi,:,curr_subregressor), ...
                fluor_taskpred_k_roi_stage_animals(curr_regressor,:,:), ...
                'uni',false),[3,1,2]));

            curr_col = min(1,linspace(0.5,0,n_learn_stages)' + ...
                task_regressor_cols{curr_regressor}(curr_subregressor,:));

            nexttile;
            AP_errorfill(task_regressor_t_shifts{curr_regressor}, ...
                squeeze(nanmean(curr_k,1)),squeeze(AP_sem(curr_k,1)),curr_col);
            title(sprintf('%s %d',task_regressor_labels{curr_regressor}, ...
                curr_subregressor));
            xlabel('Time from event');
            ylabel(sprintf('%s \\DeltaF/F_0',wf_roi(curr_roi).area));
            xline(0)
        end
    end
end
% (link axes with same ROI)
ax = reshape(allchild(h),[],length(plot_rois))';
for curr_roi = 1:length(plot_rois)
   linkaxes(ax(curr_roi,:),'y'); 
end


%% ^^ Task - trial activity

plot_rois = reshape([1,6,7] + [0;size(wf_roi,1)],1,[]);
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);
n_trial_smooth = 100;

% Get indicies for averaging 
[learned_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);
[learned_stage_idx,~] = ndgrid(trial_learned_stage,1:length(t),1:n_rois);

[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);

accum_stage_idx = cat(4,learned_stage_idx,t_idx,roi_idx,animal_idx);

% (roi activity avg: learned stage x t x roi x animal)
roi_act_learn_avg = accumarray( ...
    reshape(accum_stage_idx,[],size(accum_stage_idx,4)), ...
    fluor_roi_deconv(:), ...
    [max(learned_stage_idx(:)),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

% Plot trial activity
figure;
h = tiledlayout(2*max(trial_learned_stage),length(plot_rois), ...
    'TileSpacing','compact','padding','compact');
for curr_stage = 1:max(trial_learned_stage)
    for curr_roi = plot_rois
        use_trials = find(trial_learned_stage == curr_stage);
        [~,sort_idx] = sort(move_t(use_trials));
        curr_data_sort = fluor_roi_deconv(use_trials(sort_idx),:,curr_roi);
        curr_data_sort_smooth = convn(curr_data_sort, ...
            ones(n_trial_smooth,1)./n_trial_smooth,'same');
        
        nexttile([2,1]);
        imagesc(t,[],curr_data_sort_smooth);hold on;
        colormap(AP_colormap('KWG'));
        c = prctile(reshape(fluor_roi_deconv(:,:,curr_roi),[],1),99).*[-1,1];
        caxis(c);
        xline(0,'color','r');
        plot(move_t(use_trials(sort_idx)),1:length(use_trials),'color',[0.6,0,0.6]);
        xlabel('Time from stim (s)');
        ylabel('Trial (rxn-time sorted)');
        title(sprintf('%s, stage %d',wf_roi(curr_roi).area,curr_stage));
    end    
end


%% ^^ Task - [TESTING OUT KENNETH IDEA: L/R AS HSV]

% load in c1 = L and c2 = R
figure;
h = tiledlayout(2,2);

c_sum = (c1+c2)./max(c1(:));
c_diff = (c1-c2)./(max(c1(:)-c2(:)));

nexttile;imagesc(c1);caxis(max(c1(:))*[-1,1]);colormap(gca,AP_colormap('KWG'));title('L');
nexttile;imagesc(c_diff);caxis([-1,1]);colormap(gca,AP_colormap('BWR'));title('L-R');

nexttile;
c = zeros(size(c1,1),size(c1,2),3);
c(:,:,1) = c1./max(c1(:));
c(:,:,2) = c2./max(c1(:));
cr = c+0.8*(1-max(c,[],3));
image(cr)
title('RGB')

r_hue = 0;
g_hue = 1;

c_hue = max(0,g_hue-rescale(c_diff,r_hue,g_hue,'InputMin',-1,'InputMax',1));
c_sat = mat2gray(c1,double([0,max(c1(:))/2]));
c_val = 0.95*ones(size(c1));

r = reshape(hsv2rgb([c_hue(:),c_sat(:),c_val(:)]),size(c1,1),size(c1,2),3);
nexttile;
image(r);
title('HSV')



%% ^^ Task - trial activity (hemidiff)

plot_rois = [1,7,6];
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);
n_trial_smooth = 100;

% Get indicies for averaging 
[learned_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);
[learned_stage_idx,~] = ndgrid(trial_learned_stage,1:length(t),1:n_rois);

[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);

accum_stage_idx = cat(4,learned_stage_idx,t_idx,roi_idx,animal_idx);

% (roi activity avg: learned stage x t x roi x animal)
roi_hemiflip = circshift(1:n_rois,n_rois/2);
curr_act = fluor_roi_deconv - fluor_roi_deconv(:,:,roi_hemiflip);

% Plot trial activity
figure;
h = tiledlayout(length(plot_rois),2,'TileIndexing','ColumnMajor');
for curr_stage = 1:max(trial_learned_stage)
    for curr_roi = plot_rois
        use_trials = find(trial_learned_stage == curr_stage);
        [~,sort_idx] = sort(move_t(use_trials));
        curr_data_sort = curr_act(use_trials(sort_idx),:,curr_roi);
        curr_data_sort_smooth = convn(curr_data_sort, ...
            ones(n_trial_smooth,1)./n_trial_smooth,'same');
        
        nexttile;
        imagesc(t,[],curr_data_sort_smooth);hold on;
        colormap(AP_colormap('BWR',[],1.5));
        c = prctile(reshape(curr_act(:,:,curr_roi),[],1),95).*[-1,1];
        caxis(c);
        xline(0,'color','k');
        plot(move_t(use_trials(sort_idx)),1:length(use_trials),'color',[0.6,0,0.6]);
        xlabel('Time from stim (s)');
        ylabel('Trial (rxn-time sorted)');
        title(sprintf('%s, stage %d',wf_roi(curr_roi).area,curr_stage));
    end    
end

%% ^^ Task - trial activity (minus move hemiratio)

plot_rois = [1,7,6];
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);
n_trial_smooth = 100;

% Set activity
roi_hemiflip = circshift(1:n_rois,n_rois/2);
curr_act = fluor_roi_deconv - ...
    trial_move_hemiratio.*fluor_roi_deconv(:,:,roi_hemiflip);

% Plot trial activity
figure;
h = tiledlayout(length(plot_rois),2,'TileIndexing','ColumnMajor');
for curr_stage = 1:max(trial_learned_stage)
    for curr_roi = plot_rois
        use_trials = find(trial_learned_stage == curr_stage);
        [~,sort_idx] = sort(move_t(use_trials));
        curr_data_sort = curr_act(use_trials(sort_idx),:,curr_roi);
        curr_data_sort_smooth = convn(curr_data_sort, ...
            ones(n_trial_smooth,1)./n_trial_smooth,'same');
        
        nexttile;
        imagesc(t,[],curr_data_sort_smooth);hold on;
        colormap(AP_colormap('BWR',[],1.5));
        c = prctile(reshape(curr_act(:,:,curr_roi),[],1),95).*[-1,1];
        caxis(c);
        xline(0,'color','k');
        plot(move_t(use_trials(sort_idx)),1:length(use_trials),'color',[0.6,0,0.6]);
        xlabel('Time from stim (s)');
        ylabel('Trial (rxn-time sorted)');
        title(sprintf('%s, stage %d',wf_roi(curr_roi).area,curr_stage));
    end    
end

%% ^^ Task - trial activity (move-aligned minus move hemiratio)

plot_rois = [1,7,6];
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);
n_trial_smooth = 100;

% Set activity
curr_act_stim = fluor_roi_deconv(:,:,1:size(wf_roi,1)) - ...
    trial_move_hemiratio.*fluor_roi_deconv(:,:,size(wf_roi,1)+1:end);

% Move-align activity
curr_act = nan(size(curr_act_stim),class(curr_act_stim));
move_align = move_idx - sum(t<0);
for curr_trial = find(~isnan(move_align))'
    curr_shift_frames = ...
        move_align(curr_trial) + [0:length(t)-1];
    curr_shift_frames_use = curr_shift_frames > 0 & curr_shift_frames <= length(t);
    
    curr_grab_frames = curr_shift_frames(curr_shift_frames_use);
    curr_fill_frames = find(curr_shift_frames_use,length(t));
    
    curr_act(curr_trial,curr_fill_frames,:) = ...
        curr_act_stim(curr_trial,curr_grab_frames,:);
end
    

% Plot trial activity
figure;
h = tiledlayout(length(plot_rois),2,'TileIndexing','ColumnMajor');
for curr_stage = 1:max(trial_learned_stage)
    for curr_roi = plot_rois
        use_trials = find(trial_learned_stage == curr_stage);
        [~,sort_idx] = sort(move_t(use_trials));
        curr_data_sort = curr_act(use_trials(sort_idx),:,curr_roi);
        curr_data_sort_smooth = convn(curr_data_sort, ...
            ones(n_trial_smooth,1)./n_trial_smooth,'same');
        
        nexttile;
        imagesc(t,[],curr_data_sort_smooth);hold on;
        colormap(AP_colormap('BWR',[],1.5));
        c = prctile(reshape(curr_act(:,:,curr_roi),[],1),95).*[-1,1];
        caxis(c);
        xline(0,'color','k');
        plot(-move_t(use_trials(sort_idx)),1:length(use_trials),'color',[0.6,0,0.6]);
        xlabel('Time from stim (s)');
        ylabel('Trial (rxn-time sorted)');
        title(sprintf('%s, stage %d',wf_roi(curr_roi).area,curr_stage));
    end    
end


%% (was a test: stim (without move) and move broken and plotted together

%%%%%%%%%% TO DO HERE: split and concat stim/move align
%%% (just quick and dirty here to play with it first)

t_leeway = 0;

% (stim/move: nan-out overlapping times)
curr_act_tostim = curr_act;
curr_act_frommove = curr_act_move;
for curr_trial = 1:size(curr_act,1)
    if isnan(move_t(curr_trial))
        continue
    end
    curr_act_tostim(curr_trial,t >= move_t(curr_trial),:) = NaN;
    curr_act_frommove(curr_trial,t < -t_leeway,:) = NaN;
end

plot_rois = [1,7,6];
figure;
h = tiledlayout(2,length(plot_rois),'TileIndexing','columnmajor');
for curr_roi = plot_rois
    for curr_stage = 1:2
        nexttile; hold on;

        curr_trials = trial_learned_stage == curr_stage;

        curr_split = nanmedian(move_t(curr_trials));
        curr_stim_split = t < curr_split + t_leeway;
        curr_move_split = t > -t_leeway;

        plot_split = curr_split;

        xline(0);
        plot(t(curr_stim_split),nanmean(curr_act_tostim(curr_trials, ...
            curr_stim_split,curr_roi),1),'r','linewidth',1);
        plot(t(curr_move_split) + plot_split, ...
            nanmean(curr_act_frommove(curr_trials, ...
            curr_move_split,curr_roi),1),'color',[0.6,0,0.6],'linewidth',1);

        xlabel('Time from stim (s)');
        ylabel(sprintf('%s_{-%s} \\DeltaF/F_0', ...
            wf_roi(curr_roi).area, ...
            wf_roi(roi_hemiflip(curr_roi)).area(end)))
    end
end
% (link axes with same ROI)
ax = reshape(allchild(h),2,length(plot_rois));
for curr_roi = 1:length(plot_rois)
   linkaxes(ax(:,curr_roi),'xy'); 
end



%% ^^ Task - trial and average activity (hemidiff: stim/move aligned)

plot_rois = [6];
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);

n_trial_smooth = 20;

% Make hemidiff activity
roi_hemiflip = circshift(1:n_rois,n_rois/2);
fluor_roi_deconv_hemidiff = fluor_roi_deconv - ...
    fluor_roi_deconv(:,:,roi_hemiflip);

% Move-align activity
fluor_roi_deconv_hemidiff_move = nan(size(fluor_roi_deconv_hemidiff),class(fluor_roi_deconv_hemidiff));
move_align = move_idx - sum(t<0);
for curr_trial = find(~isnan(move_align))'
    curr_shift_frames = ...
        move_align(curr_trial) + [0:length(t)-1];
    curr_shift_frames_use = curr_shift_frames > 0 & curr_shift_frames <= length(t);
    
    curr_grab_frames = curr_shift_frames(curr_shift_frames_use);
    curr_fill_frames = find(curr_shift_frames_use,length(t));
    
    fluor_roi_deconv_hemidiff_move(curr_trial,curr_fill_frames,:) = ...
        fluor_roi_deconv_hemidiff(curr_trial,curr_grab_frames,:);
end
    
% Get indicies for averaging 
[learned_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);
[learned_stage_idx,~] = ndgrid(trial_learned_stage,1:length(t),1:n_rois);

[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);

accum_stage_idx = cat(4,learned_stage_idx,t_idx,roi_idx,animal_idx);

% (roi activity avg: learned stage x t x roi x animal) 
roi_stim_learn_avg = accumarray( ...
    reshape(accum_stage_idx,[],size(accum_stage_idx,4)), ...
    fluor_roi_deconv_hemidiff(:), ...
    [max(learned_stage_idx(:)),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

roi_move_learn_avg = accumarray( ...
    reshape(accum_stage_idx,[],size(accum_stage_idx,4)), ...
    fluor_roi_deconv_hemidiff_move(:), ...
    [max(learned_stage_idx(:)),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

% Plot reduced activity
for curr_roi = plot_rois
    figure;
    h = tiledlayout(2*max(trial_learned_stage)+1,2, ...
        'TileSpacing','compact','padding','compact','TileIndexing','columnmajor');
    for curr_event = 1:2
        switch curr_event
            case 1
                curr_trial_data = fluor_roi_deconv_hemidiff;
                curr_avg_data = roi_stim_learn_avg;
                curr_event_label = 'Stim';
                curr_align = zeros(size(fluor_roi_deconv,1),1);
                curr_cmap = AP_colormap('KWG');
            case 2
                curr_trial_data = fluor_roi_deconv_hemidiff_move;
                curr_avg_data = roi_move_learn_avg;
                curr_event_label = 'Move';
                curr_align = move_t;
                curr_cmap = AP_colormap('KWG');
        end
        
        % Plot trial activity
        for curr_stage = 1:max(trial_learned_stage)
            use_trials = find(trial_learned_stage == curr_stage);
            [~,sort_idx] = sort(move_t(use_trials));
            curr_data_sort = curr_trial_data(use_trials(sort_idx),:,curr_roi);
            curr_data_sort_smooth = nanconv(curr_data_sort, ...
                ones(n_trial_smooth,1)./n_trial_smooth);
            
            nexttile([2,1]);
            imagesc(t,[],curr_data_sort_smooth);hold on;
            colormap(gca,curr_cmap);
            c = prctile(reshape(curr_trial_data(:,:,curr_roi),[],1),99).*[-1,1];
            caxis(c);
            plot(zeros(length(use_trials),1) - curr_align(use_trials(sort_idx)), ...
                1:length(use_trials),'color','r');
            plot(move_t(use_trials(sort_idx)) - curr_align(use_trials(sort_idx)), ...
                1:length(use_trials),'color',[0.6,0,0.6]);
            xlabel(sprintf('Time from %s (s)',curr_event_label));
            ylabel('Trial (rxn-time sorted)');
            title(sprintf('%s, stage %d',curr_event_label,curr_stage));
        end

        % Plot average
        stage_col = curr_cmap(end-round(prctile(1:((size(curr_cmap,1)-1)/2),[75,25])),:);
        nexttile;
        AP_errorfill(t, ...
            permute(nanmean(curr_avg_data(:,:,curr_roi,:),4),[2,1,3]), ...
            permute(AP_sem(curr_avg_data(:,:,curr_roi,:),4),[2,1,3]), ...
            stage_col);
        axis tight; xlim(xlim+[-0.1,0.1]);
        xlabel(sprintf('Time from %s (s)',curr_event_label));
        ylabel('\DeltaF/F_0');
        title(curr_event_label);
        xline(0,'linestyle','--');
        
        title(h,wf_roi(curr_roi).area);
    end
end

%% ^^ Task - average activity (stim/move-aligned, raw/hemidiff)

plot_rois = [1,6,7];
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);

% Make hemidiff/sum activity
roi_hemiflip = circshift(1:n_rois,n_rois/2);
fluor_roi_deconv_hemidiff = fluor_roi_deconv - ...
    fluor_roi_deconv(:,:,roi_hemiflip);
fluor_roi_deconv_hemisum = fluor_roi_deconv + ...
    fluor_roi_deconv(:,:,roi_hemiflip);

% Move-align activity
fluor_roi_deconv_move = nan(size(fluor_roi_deconv),class(fluor_roi_deconv));
fluor_roi_deconv_hemidiff_move = nan(size(fluor_roi_deconv_hemidiff),class(fluor_roi_deconv_hemidiff));
fluor_roi_deconv_hemisum_move = nan(size(fluor_roi_deconv_hemidiff),class(fluor_roi_deconv_hemisum));
move_align = move_idx - sum(t<0);
for curr_trial = find(~isnan(move_align))'
    curr_shift_frames = ...
        move_align(curr_trial) + [0:length(t)-1];
    curr_shift_frames_use = curr_shift_frames > 0 & curr_shift_frames <= length(t);
    
    curr_grab_frames = curr_shift_frames(curr_shift_frames_use);
    curr_fill_frames = find(curr_shift_frames_use,length(t));
    
    fluor_roi_deconv_move(curr_trial,curr_fill_frames,:) = ...
        fluor_roi_deconv(curr_trial,curr_grab_frames,:);
    fluor_roi_deconv_hemidiff_move(curr_trial,curr_fill_frames,:) = ...
        fluor_roi_deconv_hemidiff(curr_trial,curr_grab_frames,:);
    fluor_roi_deconv_hemisum_move(curr_trial,curr_fill_frames,:) = ...
        fluor_roi_deconv_hemisum(curr_trial,curr_grab_frames,:);
end

% Get indicies for averaging (learned stage x t x roi x animal)
[learned_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);
[learned_stage_idx,~] = ndgrid(trial_learned_stage,1:length(t),1:n_rois);

[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);

accum_stage_idx = cat(4,learned_stage_idx,t_idx,roi_idx,animal_idx);

% Average activity and plot
stage_col = [0.5,0.5,0.5;0,0,0];

figure;
h = tiledlayout(4,length(plot_rois), ...
    'TileSpacing','compact','padding','compact','TileIndexing','columnmajor');
for curr_roi = plot_rois
    for curr_act_type = 1:4
        switch curr_act_type
            case 1
                curr_data = fluor_roi_deconv_hemisum;
                curr_event = 'stim';
                curr_roi_label = sprintf('%s{%s+%s}',wf_roi(curr_roi).area(1:end-1), ...
                    wf_roi(curr_roi).area(end),wf_roi(roi_hemiflip(curr_roi)).area(end));
            case 2
                curr_data = fluor_roi_deconv_hemidiff;
                curr_event = 'stim';
                curr_roi_label = sprintf('%s{%s-%s}',wf_roi(curr_roi).area(1:end-1), ...
                    wf_roi(curr_roi).area(end),wf_roi(roi_hemiflip(curr_roi)).area(end));
            case 3
                curr_data = fluor_roi_deconv_hemisum_move;
                curr_event = 'move';
                curr_roi_label = sprintf('%s{%s+%s}',wf_roi(curr_roi).area(1:end-1), ...
                    wf_roi(curr_roi).area(end),wf_roi(roi_hemiflip(curr_roi)).area(end));
            case 4
                curr_data = fluor_roi_deconv_hemidiff_move;
                curr_event = 'move';
                curr_roi_label = sprintf('%s{%s-%s}',wf_roi(curr_roi).area(1:end-1), ...
                    wf_roi(curr_roi).area(end),wf_roi(roi_hemiflip(curr_roi)).area(end));
        end

        curr_data_avg = accumarray( ...
            reshape(accum_stage_idx,[],size(accum_stage_idx,4)), ...
            curr_data(:), ...
            [max(learned_stage_idx(:)),length(t),n_rois,length(animals)], ...
            @nanmean,NaN('single'));

        hl = nexttile;
        AP_errorfill(t, ...
            permute(nanmean(curr_data_avg(:,:,curr_roi,:),4),[2,1,3]), ...
            permute(AP_sem(curr_data_avg(:,:,curr_roi,:),4),[2,1,3]), ...
            stage_col);
        axis tight; xlim(xlim+[-0.1,0.1]);
        xlabel(sprintf('Time from %s (s)',curr_event));
        ylabel(sprintf('%s %s \\DeltaF/F_0',curr_roi_label));
        xline(0,'linestyle','--');
        drawnow;

    end
end
% (link axes with same ROI)
ax = reshape(allchild(h),4,length(plot_rois));
for curr_comb = 1:4
   linkaxes(ax(curr_comb,:),'xy'); 
end

%% ^^ Task - average activity (learning change stim/move-align)

plot_rois = [1,6,7];
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);

% Make hemidiff/sum activity
roi_hemiflip = circshift(1:n_rois,n_rois/2);
fluor_roi_deconv_hemidiff = fluor_roi_deconv - fluor_roi_deconv(:,:,roi_hemiflip);

% Move-align activity
fluor_roi_deconv_hemidiff_move = nan(size(fluor_roi_deconv_hemidiff),class(fluor_roi_deconv_hemidiff));
move_align = move_idx - sum(t<0);
for curr_trial = find(~isnan(move_align))'
    curr_shift_frames = ...
        move_align(curr_trial) + [0:length(t)-1];
    curr_shift_frames_use = curr_shift_frames > 0 & curr_shift_frames <= length(t);
    
    curr_grab_frames = curr_shift_frames(curr_shift_frames_use);
    curr_fill_frames = find(curr_shift_frames_use,length(t));   

    fluor_roi_deconv_hemidiff_move(curr_trial,curr_fill_frames,:) = ...
        fluor_roi_deconv_hemidiff(curr_trial,curr_grab_frames,:);

end

% Get indicies for averaging (learned stage x t x roi x animal)
[learned_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);
[learned_stage_idx,~] = ndgrid(trial_learned_stage,1:length(t),1:n_rois);

[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);

accum_stage_idx = cat(4,learned_stage_idx,t_idx,roi_idx,animal_idx);

% (roi activity avg: learned stage x t x roi x animal) 
roi_hemidiff_stim_learn_avg = accumarray( ...
    reshape(accum_stage_idx,[],size(accum_stage_idx,4)), ...
    fluor_roi_deconv_hemidiff(:), ...
    [max(learned_stage_idx(:)),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

roi_hemidiff_move_learn_avg = accumarray( ...
    reshape(accum_stage_idx,[],size(accum_stage_idx,4)), ...
    fluor_roi_deconv_hemidiff_move(:), ...
    [max(learned_stage_idx(:)),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

% Average activity and plot
figure;
h = tiledlayout(1,length(plot_rois)*2);
for curr_roi = plot_rois

    nexttile;
    curr_data = permute(diff(roi_hemidiff_stim_learn_avg(:,:,curr_roi,:),[],1),[2,4,1,3]);
    AP_errorfill(t,nanmean(curr_data,2),AP_sem(curr_data,2),'r');
    axis tight; xlim(xlim+[-0.1,0.1]);
    xlabel('Time from stim (s)');
    ylabel(sprintf('%s_{L-R}',wf_roi(curr_roi).area(1:end-2)));
    xline(0,'linestyle','--');
%     xlim([-0.2,median(move_t)])

    nexttile;
    curr_data = permute(diff(roi_hemidiff_move_learn_avg(:,:,curr_roi,:),[],1),[2,4,1,3]);
    AP_errorfill(t,nanmean(curr_data,2),AP_sem(curr_data,2),'r');
    axis tight; xlim(xlim+[-0.1,0.1]);
    xlabel('Time from move (s)');
    ylabel(sprintf('%s_{L-R}',wf_roi(curr_roi).area(1:end-2)));
    xline(0,'linestyle','--');
%     xlim([0,median(move_t)+0.2])

end
% (link axes with same ROI)
ax = reshape(allchild(h),2,length(plot_rois));
for curr_roi = 1:length(plot_rois)
   linkaxes(ax(:,curr_roi),'y'); 
end


%% ^^ Task - trial and average activity (stim/move-reduced)

plot_rois = [6];
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);

n_trial_smooth = 20;

% Set task-predicted NaNs to 0 (no task prediction = subtract nothing)
fluor_roi_taskpred_reduced_nonan = fluor_roi_taskpred_reduced;
fluor_roi_taskpred_reduced_nonan(isnan(fluor_roi_taskpred_reduced_nonan)) = 0;

% Reduce activity: stim (set taskpred nan = 0)
stim_regressor = strcmp(task_regressor_labels,'Stim');
fluor_roi_deconv_stim = fluor_roi_deconv - fluor_roi_taskpred_reduced_nonan(:,:,:,stim_regressor);

% Reduce activity: move & move-align (set taskpred nan = 0)
move_regressor = strcmp(task_regressor_labels,'Move onset');
fluor_roi_deconv_move = nan(size(fluor_roi_deconv),class(fluor_roi_deconv));
move_align = move_idx - sum(t<0);
for curr_trial = find(~isnan(move_align))'
    curr_shift_frames = ...
        move_align(curr_trial) + [0:length(t)-1];
    curr_shift_frames_use = curr_shift_frames > 0 & curr_shift_frames <= length(t);
    
    curr_grab_frames = curr_shift_frames(curr_shift_frames_use);
    curr_fill_frames = find(curr_shift_frames_use,length(t));
    
    fluor_roi_deconv_move(curr_trial,curr_fill_frames,:) = ...
        fluor_roi_deconv(curr_trial,curr_grab_frames,:) - ...
        fluor_roi_taskpred_reduced_nonan(curr_trial,curr_grab_frames,:,move_regressor);

end
    
% Get indicies for averaging 
[learned_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);
[learned_stage_idx,~] = ndgrid(trial_learned_stage,1:length(t),1:n_rois);

[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);

accum_stage_idx = cat(4,learned_stage_idx,t_idx,roi_idx,animal_idx);

% (roi activity avg: learned stage x t x roi x animal) 
roi_stim_learn_avg = accumarray( ...
    reshape(accum_stage_idx,[],size(accum_stage_idx,4)), ...
    fluor_roi_deconv_stim(:), ...
    [max(learned_stage_idx(:)),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

roi_move_learn_avg = accumarray( ...
    reshape(accum_stage_idx,[],size(accum_stage_idx,4)), ...
    fluor_roi_deconv_move(:), ...
    [max(learned_stage_idx(:)),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

% Plot reduced activity
for curr_roi = plot_rois
    figure;
    h = tiledlayout(2*max(trial_learned_stage)+1,2, ...
        'TileSpacing','compact','padding','compact','TileIndexing','columnmajor');
    for curr_event = 1:2
        switch curr_event
            case 1
                curr_trial_data = fluor_roi_deconv_stim;
                curr_avg_data = roi_stim_learn_avg;
                curr_event_label = 'Stim';
                curr_align = zeros(size(fluor_roi_deconv,1),1);
                curr_cmap = AP_colormap('BWR');
            case 2
                curr_trial_data = fluor_roi_deconv_move;
                curr_avg_data = roi_move_learn_avg;
                curr_event_label = 'Move';
                curr_align = move_t;
                curr_cmap = AP_colormap('KWP');
        end
        
        % Plot trial activity
        for curr_stage = 1:max(trial_learned_stage)
            use_trials = find(trial_learned_stage == curr_stage);
            [~,sort_idx] = sort(move_t(use_trials));
            curr_data_sort = curr_trial_data(use_trials(sort_idx),:,curr_roi);
            curr_data_sort_smooth = nanconv(curr_data_sort, ...
                ones(n_trial_smooth,1)./n_trial_smooth);
            
            nexttile([2,1]);
            imagesc(t,[],curr_data_sort_smooth);hold on;
            colormap(gca,curr_cmap);
            c = prctile(reshape(curr_trial_data(:,:,curr_roi),[],1),99).*[-1,1];
            caxis(c);
            plot(zeros(length(use_trials),1) - curr_align(use_trials(sort_idx)), ...
                1:length(use_trials),'color','r');
            plot(move_t(use_trials(sort_idx)) - curr_align(use_trials(sort_idx)), ...
                1:length(use_trials),'color',[0.6,0,0.6]);
            xlabel(sprintf('Time from %s (s)',curr_event_label));
            ylabel('Trial (rxn-time sorted)');
            title(sprintf('%s, stage %d',curr_event_label,curr_stage));
        end

        % Plot average
        stage_col = curr_cmap(end-round(prctile(1:((size(curr_cmap,1)-1)/2),[75,25])),:);
        nexttile;
        AP_errorfill(t, ...
            permute(nanmean(curr_avg_data(:,:,curr_roi,:),4),[2,1,3]), ...
            permute(AP_sem(curr_avg_data(:,:,curr_roi,:),4),[2,1,3]), ...
            stage_col);
        axis tight; xlim(xlim+[-0.1,0.1]);
        xlabel(sprintf('Time from %s (s)',curr_event_label));
        ylabel('\DeltaF/F_0');
        title(curr_event_label);
        xline(0,'linestyle','--');
        
        title(h,wf_roi(curr_roi).area);
    end
end


%% ^^ Task - average activity (L/R/diff, stim/move align)

plot_rois = [1,7,6];
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);

% Move-align activity
fluor_roi_deconv_move = nan(size(fluor_roi_deconv),class(fluor_roi_deconv));
move_align = move_idx - sum(t<0);
for curr_trial = find(~isnan(move_align))'
    curr_shift_frames = ...
        move_align(curr_trial) + [0:length(t)-1];
    curr_shift_frames_use = curr_shift_frames > 0 & curr_shift_frames <= length(t);
    
    curr_grab_frames = curr_shift_frames(curr_shift_frames_use);
    curr_fill_frames = find(curr_shift_frames_use,length(t));
    
    fluor_roi_deconv_move(curr_trial,curr_fill_frames,:) = ...
        fluor_roi_deconv(curr_trial,curr_grab_frames,:);
end

% Get indicies for averaging (learned stage x t x roi x animal)
[learned_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);
[learned_stage_idx,~] = ndgrid(trial_learned_stage,1:length(t),1:n_rois);

[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);

accum_stage_idx = cat(4,learned_stage_idx,t_idx,roi_idx,animal_idx);

% (roi activity avg: learned stage x t x roi x animal) 
roi_stim_learn_avg = accumarray( ...
    reshape(accum_stage_idx,[],size(accum_stage_idx,4)), ...
    fluor_roi_deconv(:), ...
    [max(learned_stage_idx(:)),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

roi_move_learn_avg = accumarray( ...
    reshape(accum_stage_idx,[],size(accum_stage_idx,4)), ...
    fluor_roi_deconv_move(:), ...
    [max(learned_stage_idx(:)),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

% Average activity and plot (stim + move)
figure;
h = tiledlayout(length(plot_rois),4);
for curr_roi = plot_rois
    for curr_stage = 1:max(trial_learned_stage)

        curr_l_data_stim = permute(roi_stim_learn_avg(curr_stage,:,curr_roi,:),[2,4,1,3]);
        curr_r_data_stim = permute(roi_stim_learn_avg(curr_stage,:,curr_roi+size(wf_roi,1),:),[2,4,1,3]);
        curr_diff_data_stim = curr_l_data_stim - curr_r_data_stim;

        curr_l_data_move = permute(roi_move_learn_avg(curr_stage,:,curr_roi,:),[2,4,1,3]);
        curr_r_data_move = permute(roi_move_learn_avg(curr_stage,:,curr_roi+size(wf_roi,1),:),[2,4,1,3]);
        curr_diff_data_move = curr_l_data_move - curr_r_data_move;

        nexttile; hold on;
        AP_errorfill(t,nanmean(curr_l_data_stim,2),AP_sem(curr_l_data_stim,2),'r');
        AP_errorfill(t,nanmean(curr_r_data_stim,2),AP_sem(curr_r_data_stim,2),'b');
        AP_errorfill(t,nanmean(curr_diff_data_stim,2),AP_sem(curr_diff_data_stim,2),'k');
        title(sprintf('Stage %d',curr_stage));
        xlabel('Stim');
        ylabel(wf_roi(curr_roi).area(1:end-2));

        nexttile; hold on;
        AP_errorfill(t,nanmean(curr_l_data_move,2),AP_sem(curr_l_data_move,2),'r');
        AP_errorfill(t,nanmean(curr_r_data_move,2),AP_sem(curr_r_data_move,2),'b');
        AP_errorfill(t,nanmean(curr_diff_data_move,2),AP_sem(curr_diff_data_move,2),'k');
        title(sprintf('Stage %d',curr_stage));
        xlabel('Move');
        ylabel(wf_roi(curr_roi).area(1:end-2));

    end
end
% (link axes with same ROI)
ax = reshape(allchild(h),4,length(plot_rois))';
for curr_roi = 1:length(plot_rois)
   linkaxes(ax(curr_roi,:),'y'); 
end


% Average activity and plot (stim)
figure;
h = tiledlayout(2,length(plot_rois),'TileIndexing','columnmajor');
for curr_roi = plot_rois
    for curr_stage = 1:max(trial_learned_stage)

        curr_l_data_stim = permute(roi_stim_learn_avg(curr_stage,:,curr_roi,:),[2,4,1,3]);
        curr_r_data_stim = permute(roi_stim_learn_avg(curr_stage,:,curr_roi+size(wf_roi,1),:),[2,4,1,3]);
        curr_diff_data_stim = curr_l_data_stim - curr_r_data_stim;

        curr_l_data_move = permute(roi_move_learn_avg(curr_stage,:,curr_roi,:),[2,4,1,3]);
        curr_r_data_move = permute(roi_move_learn_avg(curr_stage,:,curr_roi+size(wf_roi,1),:),[2,4,1,3]);
        curr_diff_data_move = curr_l_data_move - curr_r_data_move;

        nexttile; hold on;
        AP_errorfill(t,nanmean(curr_l_data_stim,2),AP_sem(curr_l_data_stim,2),'r');
        AP_errorfill(t,nanmean(curr_r_data_stim,2),AP_sem(curr_r_data_stim,2),'b');
        AP_errorfill(t,nanmean(curr_diff_data_stim,2),AP_sem(curr_diff_data_stim,2),'k');
        title(sprintf('Stage %d',curr_stage));
        xlabel('Stim');
        ylabel(wf_roi(curr_roi).area(1:end-2));
        xline(0);set(gca,'children',circshift(get(gca,'children'),-1))

    end
end
% (link axes with same ROI)
ax = reshape(allchild(h),2,length(plot_rois))';
for curr_roi = 1:length(plot_rois)
   linkaxes(ax(curr_roi,:),'y'); 
end

%% ^^ Task - average activity (L/R stim/move stage, L-R stage/change)

plot_rois = [1,7,6];
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);

% Move-align activity
fluor_roi_deconv_move = nan(size(fluor_roi_deconv),class(fluor_roi_deconv));
wheel_allcat_move = nan(size(wheel_allcat),class(wheel_allcat));

move_align = move_idx - sum(t<0);
for curr_trial = find(~isnan(move_align))'
    curr_shift_frames = ...
        move_align(curr_trial) + [0:length(t)-1];
    curr_shift_frames_use = curr_shift_frames > 0 & curr_shift_frames <= length(t);
    
    curr_grab_frames = curr_shift_frames(curr_shift_frames_use);
    curr_fill_frames = find(curr_shift_frames_use,length(t));
    
    fluor_roi_deconv_move(curr_trial,curr_fill_frames,:) = ...
        fluor_roi_deconv(curr_trial,curr_grab_frames,:);
    wheel_allcat_move(curr_trial,curr_fill_frames,:) = ...
        wheel_allcat(curr_trial,curr_grab_frames,:);
end

% Get indicies for averaging (learned stage x t x roi x animal)
[learned_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);
[learned_stage_idx,~] = ndgrid(trial_learned_stage,1:length(t),1:n_rois);

[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);

accum_stage_idx = cat(4,learned_stage_idx,t_idx,roi_idx,animal_idx);

% (roi activity avg: learned stage x t x roi x animal) 
roi_stim_learn_avg = accumarray( ...
    reshape(accum_stage_idx,[],size(accum_stage_idx,4)), ...
    fluor_roi_deconv(:), ...
    [max(trial_learned_stage),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

roi_move_learn_avg = accumarray( ...
    reshape(accum_stage_idx,[],size(accum_stage_idx,4)), ...
    fluor_roi_deconv_move(:), ...
    [max(trial_learned_stage),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

% (wheel avg: learned stage x t x animal)
wheel_accum_idx = ...
    cat(3,repmat(trial_learned_stage,1,length(t)), ...
    repmat(1:length(t),length(trial_animal),1), ...
    repmat(trial_animal,1,length(t)));

wheel_stim_learn_avg = accumarray( ...
    reshape(wheel_accum_idx,[],size(wheel_accum_idx,3)), ...
    wheel_allcat(:), ...
    [max(trial_learned_stage),length(t),length(animals)], ...
    @nanmean,NaN);

wheel_move_learn_avg = accumarray( ...
    reshape(wheel_accum_idx,[],size(wheel_accum_idx,3)), ...
    wheel_allcat_move(:), ...
    [max(trial_learned_stage),length(t),length(animals)], ...
    @nanmean,NaN);

% Plot L/R stim/move activity, separately by learning
median_move_t = median(move_t(trial_learned_day>0));
x_limit = [-0.1,median_move_t];

figure;
h = tiledlayout(2,length(plot_rois)*2);
for curr_stage = 1:max(trial_learned_stage)
    for curr_roi = plot_rois

        curr_l_data_stim = permute(roi_stim_learn_avg(curr_stage,:,curr_roi,:),[2,4,1,3]);
        curr_r_data_stim = permute(roi_stim_learn_avg(curr_stage,:,curr_roi+size(wf_roi,1),:),[2,4,1,3]);
        curr_diff_data_stim = curr_l_data_stim - curr_r_data_stim;

        curr_l_data_move = permute(roi_move_learn_avg(curr_stage,:,curr_roi,:),[2,4,1,3]);
        curr_r_data_move = permute(roi_move_learn_avg(curr_stage,:,curr_roi+size(wf_roi,1),:),[2,4,1,3]);
        curr_diff_data_move = curr_l_data_move - curr_r_data_move;

        nexttile; hold on;
        AP_errorfill(t,nanmean(curr_l_data_stim,2),AP_sem(curr_l_data_stim,2),'r');
        AP_errorfill(t,nanmean(curr_r_data_stim,2),AP_sem(curr_r_data_stim,2),'b');
        title(sprintf('Stage %d',curr_stage));
        xlabel('Stim');
        ylabel(wf_roi(curr_roi).area(1:end-2));
        xline(0);set(gca,'children',circshift(get(gca,'children'),-1))
        xlim(x_limit);

        nexttile; hold on;
        AP_errorfill(t,nanmean(curr_l_data_move,2),AP_sem(curr_l_data_move,2),'r');
        AP_errorfill(t,nanmean(curr_r_data_move,2),AP_sem(curr_r_data_move,2),'b');
        title(sprintf('Stage %d',curr_stage));
        xlabel('Move');
        ylabel(wf_roi(curr_roi).area(1:end-2));
        xline(0);set(gca,'children',circshift(get(gca,'children'),-1))
        xlim(x_limit - x_limit(1));

    end
end
linkaxes(allchild(h),'y');

% Plot L-R stim/move, learning overlaid
figure;
h = tiledlayout(1,length(plot_rois)*2);
stage_col = [0.5,0.5,0.5;0,0,0];
for curr_roi = plot_rois

    curr_l_data_stim = permute(roi_stim_learn_avg(:,:,curr_roi,:),[2,1,4,3]);
    curr_r_data_stim = permute(roi_stim_learn_avg(:,:,curr_roi+size(wf_roi,1),:),[2,1,4,3]);
    curr_diff_data_stim = curr_l_data_stim - curr_r_data_stim;

    curr_l_data_move = permute(roi_move_learn_avg(:,:,curr_roi,:),[2,1,4,3]);
    curr_r_data_move = permute(roi_move_learn_avg(:,:,curr_roi+size(wf_roi,1),:),[2,1,4,3]);
    curr_diff_data_move = curr_l_data_move - curr_r_data_move;

    nexttile; hold on;
    AP_errorfill(t,nanmean(curr_diff_data_stim,3),AP_sem(curr_diff_data_stim,3),stage_col);
    xlabel('Stim');
    ylabel(wf_roi(curr_roi).area(1:end-2));
    xline(0);set(gca,'children',circshift(get(gca,'children'),-1))
    xlim(x_limit);

    nexttile; hold on;
    AP_errorfill(t,nanmean(curr_diff_data_move,3),AP_sem(curr_diff_data_move,3),stage_col);
    xlabel('Move');
    ylabel(wf_roi(curr_roi).area(1:end-2));
    xline(0);set(gca,'children',circshift(get(gca,'children'),-1))
    xlim(x_limit - x_limit(1));

end
% (link ROI y-axes)
ax = reshape(allchild(h),2,length(plot_rois));
for curr_roi = 1:length(plot_rois)
   linkaxes(ax(:,curr_roi),'y'); 
end


% Plot L-R learning change overlaid
figure;
h = tiledlayout(1,2);

curr_l_data_stim = permute(roi_stim_learn_avg(:,:,plot_rois,:),[2,3,1,4]);
curr_r_data_stim = permute(roi_stim_learn_avg(:,:,plot_rois+size(wf_roi,1),:),[2,3,1,4]);
curr_diff_learndiff_data_stim = diff(curr_l_data_stim,[],3) - ...
    diff(curr_r_data_stim,[],3);

curr_l_data_move = permute(roi_move_learn_avg(:,:,plot_rois,:),[2,3,1,4]);
curr_r_data_move = permute(roi_move_learn_avg(:,:,plot_rois+size(wf_roi,1),:),[2,3,1,4]);
curr_diff_learndiff_data_move = diff(curr_l_data_move,[],3) - ...
    diff(curr_r_data_move,[],3);

nexttile;
AP_errorfill(t,nanmean(curr_diff_learndiff_data_stim,4), ...
    AP_sem(curr_diff_learndiff_data_stim,4));
xlabel('Stim');
xline(0);set(gca,'children',circshift(get(gca,'children'),-1))
xlim(x_limit);

nexttile;
AP_errorfill(t,nanmean(curr_diff_learndiff_data_move,4), ...
    AP_sem(curr_diff_learndiff_data_move,4));
xlabel('Stim');
xline(0);set(gca,'children',circshift(get(gca,'children'),-1))
xlim(x_limit - x_limit(1));
legend(flipud(findobj(gca,'Type','line')),cellfun(@(x) [x(1:end-2) '_{L-R}'], ...
    {wf_roi(plot_rois).area},'uni',false),'location','ne');

linkaxes(allchild(h),'y');

% Plot wheel learning change
figure;
h = tiledlayout(1,2);

curr_wheel_learndiff_stim = squeeze(diff(wheel_stim_learn_avg,[],1));
curr_wheel_learndiff_move = squeeze(diff(wheel_move_learn_avg,[],1));

nexttile;
AP_errorfill(t,-nanmean(curr_wheel_learndiff_stim,2), ...
    AP_sem(curr_wheel_learndiff_stim,2),'k');
xlabel('Stim');
xline(0);set(gca,'children',circshift(get(gca,'children'),-1))
xlim(x_limit);
ylabel('Wheel velocity');

nexttile;
AP_errorfill(t,-nanmean(curr_wheel_learndiff_move,2), ...
    AP_sem(curr_wheel_learndiff_move,2),'k');
xlabel('Move');
xline(0);set(gca,'children',circshift(get(gca,'children'),-1))
xlim(x_limit - x_limit(1));
ylabel('Wheel velocity');


%% ^^ Task - ROI timecourse by learned day

roi_hemiflip = circshift(1:n_rois,n_rois/2);
curr_act = fluor_roi_deconv - ...
    trial_move_hemiratio.*fluor_roi_deconv(:,:,roi_hemiflip);

% Get indicies for averaging 
[learned_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);
[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);

accum_learnedday_idx = cat(4,learned_day_idx,t_idx,roi_idx,animal_idx);

use_trials = true(size(move_t));

roi_learnday_avg = accumarray( ...
    reshape(accum_learnedday_idx(use_trials,:,:,:),[],size(accum_learnedday_idx,4)), ...
    reshape(curr_act(use_trials,:,:),[],1), ...
    [max(learned_day_idx(:)),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

use_t = t > 0 & t <= 0.2;
stim_roi_learnday_avg_tmax = squeeze(max(roi_learnday_avg(:,use_t,:,:,:),[],2));


plot_roi = 6;
curr_tcourse = nanmean(squeeze(roi_learnday_avg(:,:,plot_roi,:)),3);
curr_tmax = squeeze(stim_roi_learnday_avg_tmax(:,plot_roi,:));

figure; 
subplot(1,2,1);
imagesc(t,learned_day_unique,curr_tcourse);
caxis(max(abs(caxis))*[-1,1]);
colormap(AP_colormap('BWR',[],1.5));
yline(-0.5);

subplot(1,2,2);
errorbar(learned_day_unique,nanmean(curr_tmax,2),AP_sem(curr_tmax,2),'k','linewidth',2);
xline(0);





%% ^^ Task - daysplit ROI fluorescence

% Activity to average 
% stim_roi_act = fluor_roi_deconv;

roi_hemiflip = circshift(1:n_rois,n_rois/2);
stim_roi_act = fluor_roi_deconv - fluor_roi_deconv(:,:,roi_hemiflip);

% % (testing: wheel)
% stim_roi_act = single(repmat(abs(wheel_allcat),[1,1,n_rois]));

% stim_regressor = strcmp(task_regressor_labels,'Stim');
% stim_roi_act = fluor_roi_deconv - fluor_roi_taskpred_reduced(:,:,:,stim_regressor);

% Get indicies for averaging 
[learned_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);

trial_learned_stage = discretize(trial_learned_day,[-Inf,-2,-1,0,1,2,Inf]);
[learned_stage_idx,~] = ndgrid(trial_learned_stage,1:length(t),1:n_rois);

[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);

n_daysplit = 1;
trial_daysplit_idx = cell2mat(arrayfun(@(x) ...
    min(floor(linspace(1,n_daysplit+1,x)),n_daysplit)', ...
    trials_recording,'uni',false));

[daysplit_idx,~,~] = ...
    ndgrid(trial_daysplit_idx,1:length(t),1:n_rois);
accum_learnedday_daysplit_idx = cat(4,learned_day_idx,daysplit_idx,t_idx,roi_idx,animal_idx);

% %%%%%%%%%%%% TEMP
% use_trials = move_t >= 0.1 & move_t <= 0.2;
% accum_learnedday_daysplit_idx = accum_learnedday_daysplit_idx(use_trials,:,:,:);
% stim_roi_act = stim_roi_act(use_trials,:,:);
% 
% % % (quickplot movies)
% % use_trials = move_t > 0.15 & trial_learned_day == 0;
% % curr_v = cat(1, ...
% %     nanmean(fluor_allcat_deconv(move_t > 0.15 & trial_learned_day == -2,:,:),1), ...
% %     nanmean(fluor_allcat_deconv(move_t > 0.15 & trial_learned_day == -1,:,:),1), ...
% %     nanmean(fluor_allcat_deconv(move_t > 0.15 & trial_learned_day == 0,:,:),1));
% % 
% % a = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),permute(curr_v,[3,2,1]));
% % AP_imscroll(a,t); axis image;
% % caxis(max(abs(caxis))*[-1,1]);
% % colormap(AP_colormap('KWG'));
% % AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
% 
% %%%%%%%%%%%%%%%%%

% (roi activity daysplit: learned day x (daysplit) x t x roi x animal)
stim_roi_act_learn_avg_daysplit = accumarray( ...
    reshape(accum_learnedday_daysplit_idx,[],size(accum_learnedday_daysplit_idx,4)), ...
    stim_roi_act(:), ...
    [max(learned_day_idx(:)),n_daysplit,length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

% (tmax activity: learned day x daysplit x roi x animal)
use_t = t > 0 & t <= 0.2;
stim_roi_act_tmax_daysplit = ...
    permute(max(stim_roi_act_learn_avg_daysplit(:,:,use_t,:,:),[],3),[1,2,4,5,3]);

% Plot time-max by learned day
x_day_spacing = 1; % (gaps between days)
learned_day_x_range = [min(learned_day_unique),max(learned_day_unique)];
learned_day_x = [learned_day_x_range(1):learned_day_x_range(2)]';
learned_daysplit_x = learned_day_x + linspace(0,1,n_daysplit+x_day_spacing);

plot_learned_day = sum(~isnan(stim_roi_act_tmax_daysplit(:,1,1,:)),4) >= min_n;

figure;
plot_rois = [6];
tiledlayout(length(plot_rois),1);
for curr_roi = plot_rois

    nexttile;
   
    errorbar(reshape(learned_daysplit_x(plot_learned_day,:)',[],1), ...
        reshape(padarray(nanmean( ...
        stim_roi_act_tmax_daysplit(plot_learned_day,:,curr_roi,:),4), ...
        [0,x_day_spacing],NaN,'post')',[],1), ...
        reshape(padarray(AP_sem( ...
        stim_roi_act_tmax_daysplit(plot_learned_day,:,curr_roi,:),4), ...
        [0,x_day_spacing],NaN,'post')',[],1),'k','linewidth',2,'CapSize',0);

    xlabel('Learned day');
    ylabel(sprintf('%s \\DeltaF/F_0',wf_roi(curr_roi).area));
    xline(0,'linestyle','--');
    axis tight;
    xlim(xlim+[-0.5,0.5]);

end

%% ^^ Task - daysplit ROI fluorescence (minus move hemiratio)

% Activity to average 
stim_roi_act = fluor_roi_deconv - ...
    trial_move_hemiratio.*fluor_roi_deconv(:,:,roi_hemiflip);

% Get indicies for averaging 
[learned_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);

trial_learned_stage = discretize(trial_learned_day,[-Inf,-2,-1,0,1,2,Inf]);
[learned_stage_idx,~] = ndgrid(trial_learned_stage,1:length(t),1:n_rois);

[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);

n_daysplit = 3;
x_day_spacing = 1; % (plot gaps between days)

trial_daysplit_idx = cell2mat(arrayfun(@(x) ...
    min(floor(linspace(1,n_daysplit+1,x)),n_daysplit)', ...
    trials_recording,'uni',false));

[daysplit_idx,~,~] = ...
    ndgrid(trial_daysplit_idx,1:length(t),1:n_rois);
accum_learnedday_daysplit_idx = cat(4,learned_day_idx,daysplit_idx,t_idx,roi_idx,animal_idx);

% %%%%%%%%%%%% TEMP
% use_trials = move_t >= 0.4 & move_t <= Inf;
% accum_learnedday_daysplit_idx = accum_learnedday_daysplit_idx(use_trials,:,:,:);
% stim_roi_act = stim_roi_act(use_trials,:,:);
% 
% % % (quickplot movies)
% % use_trials = move_t > 0.15 & trial_learned_day == 0;
% % curr_v = cat(1, ...
% %     nanmean(fluor_allcat_deconv(move_t > 0.15 & trial_learned_day == -2,:,:),1), ...
% %     nanmean(fluor_allcat_deconv(move_t > 0.15 & trial_learned_day == -1,:,:),1), ...
% %     nanmean(fluor_allcat_deconv(move_t > 0.15 & trial_learned_day == 0,:,:),1));
% % 
% % a = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),permute(curr_v,[3,2,1]));
% % AP_imscroll(a,t); axis image;
% % caxis(max(abs(caxis))*[-1,1]);
% % colormap(AP_colormap('KWG'));
% % AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
% 
% %%%%%%%%%%%%%%%%%

% (roi activity daysplit: learned day x (daysplit) x t x roi x animal)
stim_roi_act_learn_avg_daysplit = accumarray( ...
    reshape(accum_learnedday_daysplit_idx,[],size(accum_learnedday_daysplit_idx,4)), ...
    stim_roi_act(:), ...
    [max(learned_day_idx(:)),n_daysplit,length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

% (tmax activity: learned day x daysplit x roi x animal)
use_t = t > 0 & t <= 0.2;
stim_roi_act_tmax_daysplit = ...
    permute(max(stim_roi_act_learn_avg_daysplit(:,:,use_t,:,:),[],3),[1,2,4,5,3]);

% Plot time-max by learned day
learned_day_x_range = [min(learned_day_unique),max(learned_day_unique)];
learned_day_x = [learned_day_x_range(1):learned_day_x_range(2)]';
learned_daysplit_x = learned_day_x + ...
    (0:(n_daysplit+x_day_spacing-1))/(n_daysplit+x_day_spacing);

plot_learned_day = sum(~isnan(stim_roi_act_tmax_daysplit(:,1,1,:)),4) >= min_n;

figure;
plot_rois = [6];
tiledlayout(length(plot_rois),1);
for curr_roi = plot_rois

    nexttile;
   
    errorbar(reshape(learned_daysplit_x(plot_learned_day,:)',[],1), ...
        reshape(padarray(nanmean( ...
        stim_roi_act_tmax_daysplit(plot_learned_day,:,curr_roi,:),4), ...
        [0,x_day_spacing],NaN,'post')',[],1), ...
        reshape(padarray(AP_sem( ...
        stim_roi_act_tmax_daysplit(plot_learned_day,:,curr_roi,:),4), ...
        [0,x_day_spacing],NaN,'post')',[],1),'k','linewidth',2,'CapSize',0);

    xlabel('Learned day');
    ylabel(sprintf('%s \\DeltaF/F_0',wf_roi(curr_roi).area));
    xline(0,'linestyle','--');
    axis tight;
    xlim(xlim+[-0.5,0.5]);

end

% (plot image for combining later)
figure;
imagesc([],learned_day_unique,squeeze(stim_roi_act_tmax_daysplit(:,:,6,:)));
xlabel('Animal');
ylabel('Learned day');
title('Task');

%% ^^ Task - daysplit ROI fluorescence (L/R/diff overlaid)

warning('Daysplit: params still in progress here');

% Activity to average 
stim_roi_act = fluor_roi_deconv;

% Get indicies for averaging 
[learned_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);

trial_learned_stage = discretize(trial_learned_day,[-Inf,-2,-1,0,1,2,Inf]);
[learned_stage_idx,~] = ndgrid(trial_learned_stage,1:length(t),1:n_rois);

[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);

n_daysplit = 4;
trial_daysplit_idx = cell2mat(arrayfun(@(x) ...
    min(floor(linspace(1,n_daysplit+1,x)),n_daysplit)', ...
    trials_recording,'uni',false));

[daysplit_idx,~,~] = ...
    ndgrid(trial_daysplit_idx,1:length(t),1:n_rois);
accum_learnedday_daysplit_idx = cat(4,learned_day_idx,daysplit_idx,t_idx,roi_idx,animal_idx);

% %%%%%%%%%%%% TEMP
% use_trials = move_t > 0.15;
% accum_learnedday_daysplit_idx = accum_learnedday_daysplit_idx(use_trials,:,:,:);
% stim_roi_act = stim_roi_act(use_trials,:,:);
% 
% % % (quickplot movies)
% % use_trials = move_t > 0.15 & trial_learned_day == 0;
% % curr_v = cat(1, ...
% %     nanmean(fluor_allcat_deconv(move_t > 0.15 & trial_learned_day == -2,:,:),1), ...
% %     nanmean(fluor_allcat_deconv(move_t > 0.15 & trial_learned_day == -1,:,:),1), ...
% %     nanmean(fluor_allcat_deconv(move_t > 0.15 & trial_learned_day == 0,:,:),1));
% % 
% % a = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),permute(curr_v,[3,2,1]));
% % AP_imscroll(a,t); axis image;
% % caxis(max(abs(caxis))*[-1,1]);
% % colormap(AP_colormap('KWG'));
% % AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
% 
% %%%%%%%%%%%%%%%%%

% (roi activity daysplit: learned day x (daysplit) x t x roi x animal)
stim_roi_act_learn_avg_daysplit = accumarray( ...
    reshape(accum_learnedday_daysplit_idx,[],size(accum_learnedday_daysplit_idx,4)), ...
    stim_roi_act(:), ...
    [max(learned_day_idx(:)),n_daysplit,length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

% (tmax activity: learned day x daysplit x roi x animal)
use_t = t > 0 & t <= 0.2;
stim_roi_act_tmax_daysplit = ...
    permute(max(stim_roi_act_learn_avg_daysplit(:,:,use_t,:,:),[],3),[1,2,4,5,3]);

% (average and max of abs wheel velocity)
wheel_accum_idx = cat(3,learned_day_idx(:,:,1),daysplit_idx(:,:,1),t_idx(:,:,1),animal_idx(:,:,1));
wheel_learn_avg_daysplit = accumarray( ...
    reshape(wheel_accum_idx,[],size(wheel_accum_idx,3)), ...
    abs(wheel_allcat(:)), ...
    [max(learned_day_idx(:)),n_daysplit,length(t),length(animals)], ...
    @nanmean,NaN('double'));
wheel_tmax_daysplit = ...
    permute(max(wheel_learn_avg_daysplit(:,:,use_t,:,:),[],3),[1,2,4,3]);

% Plot time-max by learned day
x_day_spacing = 1; % (gaps between days)
learned_day_x_range = [min(learned_day_unique),max(learned_day_unique)];
learned_day_x = [learned_day_x_range(1):learned_day_x_range(2)]';
learned_daysplit_x = learned_day_x + linspace(0,1,n_daysplit+x_day_spacing);

plot_learned_day = sum(~isnan(stim_roi_act_tmax_daysplit(:,1,1,:)),4) >= min_n;
learned_daysplit_x_long = reshape(learned_daysplit_x(plot_learned_day,:)',[],1);

figure;
plot_rois = [1,7,6];
tiledlayout(length(plot_rois)+1,1);
for curr_roi = plot_rois

    nexttile; hold on;

    curr_l_act = reshape(permute(padarray( ...
        stim_roi_act_tmax_daysplit(plot_learned_day,:,curr_roi,:), ...
        [0,x_day_spacing],NaN,'post'),[2,1,4,3]),[],length(animals));

    curr_r_act = reshape(permute(padarray( ...
        stim_roi_act_tmax_daysplit(plot_learned_day,:,curr_roi+size(wf_roi,1),:), ...
        [0,x_day_spacing],NaN,'post'),[2,1,4,3]),[],length(animals));

    errorbar(learned_daysplit_x_long,nanmean(curr_l_act,2),...
        AP_sem(curr_l_act,2),'r','linewidth',2,'CapSize',0);
    errorbar(learned_daysplit_x_long,nanmean(curr_r_act,2),...
        AP_sem(curr_r_act,2),'b','linewidth',2,'CapSize',0);
    errorbar(learned_daysplit_x_long,nanmean(curr_l_act-curr_r_act,2),...
        AP_sem(curr_l_act-curr_r_act,2),'k','linewidth',2,'CapSize',0);

    xlabel('Learned day');
    ylabel(sprintf('%s \\DeltaF/F_0',wf_roi(curr_roi).area(1:end-1)));
    xline(0,'linestyle','--');
    axis tight;
    xlim(xlim+[-0.5,0.5]);
    legend({'L hemi','R hemi','Diff'},'location','nw')

end

nexttile;
wheel_act_long = reshape(permute(padarray( ...
    wheel_tmax_daysplit(plot_learned_day,:,:), ...
    [0,x_day_spacing],NaN,'post'),[2,1,4,3]),[],length(animals));
errorbar(learned_daysplit_x_long, ...
    nanmean(wheel_act_long,2),AP_sem(wheel_act_long,2),'k','linewidth',2,'CapSize',0);
xlabel('Learned day');
ylabel('Wheel velocity');
xline(0,'linestyle','--');
axis tight;
xlim(xlim+[-0.5,0.5]);



%% ^^ Task - daysplit ROI fluorescence (move-aligned: sanity check)

warning('Daysplit: params still in progress here');

% Move-align activity
wheel_allcat_move = nan(size(wheel_allcat));
fluor_roi_deconv_move = nan(size(fluor_roi_deconv),class(fluor_roi_deconv));
move_align = move_idx - sum(t<0);
for curr_trial = find(~isnan(move_align))'
    curr_shift_frames = ...
        move_align(curr_trial) + [0:length(t)-1];
    curr_shift_frames_use = curr_shift_frames > 0 & curr_shift_frames <= length(t);
    
    curr_grab_frames = curr_shift_frames(curr_shift_frames_use);
    curr_fill_frames = find(curr_shift_frames_use,length(t));
    
    wheel_allcat_move(curr_trial,curr_fill_frames) = ...
        wheel_allcat(curr_trial,curr_grab_frames);
    fluor_roi_deconv_move(curr_trial,curr_fill_frames,:) = ...
        fluor_roi_deconv(curr_trial,curr_grab_frames,:);
end

% Activity to average 
% stim_roi_act = fluor_roi_deconv_move;

% roi_hemiflip = circshift(1:n_rois,n_rois/2);
% stim_roi_act = fluor_roi_deconv_move - fluor_roi_deconv_move(:,:,roi_hemiflip);

roi_hemiflip = circshift(1:n_rois,n_rois/2);
stim_roi_act = fluor_roi_deconv_move - ...
    trial_move_hemiratio.*fluor_roi_deconv_move(:,:,roi_hemiflip);

% stim_roi_act = single(repmat(abs(wheel_allcat_move),[1,1,n_rois]));

% Get indicies for averaging 
[learned_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);

trial_learned_stage = discretize(trial_learned_day,[-Inf,-2,-1,0,1,2,Inf]);
[learned_stage_idx,~] = ndgrid(trial_learned_stage,1:length(t),1:n_rois);

[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);

n_daysplit = 4;
trial_daysplit_idx = cell2mat(arrayfun(@(x) ...
    min(floor(linspace(1,n_daysplit+1,x)),n_daysplit)', ...
    trials_recording,'uni',false));

[daysplit_idx,~,~] = ...
    ndgrid(trial_daysplit_idx,1:length(t),1:n_rois);
accum_learnedday_daysplit_idx = cat(4,learned_day_idx,daysplit_idx,t_idx,roi_idx,animal_idx);

% %%%%%%%%%%%% TEMP
% use_trials = move_t > 0.25;
% accum_learnedday_daysplit_idx = accum_learnedday_daysplit_idx(use_trials,:,:,:);
% stim_roi_act = stim_roi_act(use_trials,:,:);
% 
% % % (quickplot movies)
% % use_trials = move_t > 0.15 & trial_learned_day == 0;
% % curr_v = cat(1, ...
% %     nanmean(fluor_allcat_deconv(move_t > 0.15 & trial_learned_day == -2,:,:),1), ...
% %     nanmean(fluor_allcat_deconv(move_t > 0.15 & trial_learned_day == -1,:,:),1), ...
% %     nanmean(fluor_allcat_deconv(move_t > 0.15 & trial_learned_day == 0,:,:),1));
% % 
% % a = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),permute(curr_v,[3,2,1]));
% % AP_imscroll(a,t); axis image;
% % caxis(max(abs(caxis))*[-1,1]);
% % colormap(AP_colormap('KWG'));
% % AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
% 
% %%%%%%%%%%%%%%%%%

% (roi activity daysplit: learned day x (daysplit) x t x roi x animal)
stim_roi_act_learn_avg_daysplit = accumarray( ...
    reshape(accum_learnedday_daysplit_idx,[],size(accum_learnedday_daysplit_idx,4)), ...
    stim_roi_act(:), ...
    [max(learned_day_idx(:)),n_daysplit,length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

% (tmax activity: learned day x daysplit x roi x animal)
use_t = t >= -0.2 & t <= 0.2;
stim_roi_act_tmax_daysplit = ...
    permute(max(stim_roi_act_learn_avg_daysplit(:,:,use_t,:,:),[],3),[1,2,4,5,3]);

% Plot time-max by learned day
x_day_spacing = 1; % (gaps between days)
learned_day_x_range = [min(learned_day_unique),max(learned_day_unique)];
learned_day_x = [learned_day_x_range(1):learned_day_x_range(2)]';
learned_daysplit_x = learned_day_x + linspace(0,1,n_daysplit+x_day_spacing);

plot_learned_day = sum(~isnan(stim_roi_act_tmax_daysplit(:,1,1,:)),4) >= min_n;

figure;
plot_rois = [1,7,6];
tiledlayout(length(plot_rois),1,'TileSpacing','compact','padding','compact');
for curr_roi = plot_rois

    nexttile;
   
    errorbar(reshape(learned_daysplit_x(plot_learned_day,:)',[],1), ...
        reshape(padarray(nanmean( ...
        stim_roi_act_tmax_daysplit(plot_learned_day,:,curr_roi,:),4), ...
        [0,x_day_spacing],NaN,'post')',[],1), ...
        reshape(padarray(AP_sem( ...
        stim_roi_act_tmax_daysplit(plot_learned_day,:,curr_roi,:),4), ...
        [0,x_day_spacing],NaN,'post')',[],1),'k','linewidth',2,'CapSize',0);

    xlabel('Learned day');
    ylabel(sprintf('%s \\DeltaF/F_0',wf_roi(curr_roi).area));
    xline(0,'linestyle','--');
    axis tight;
    xlim(xlim+[-0.5,0.5]);

end

%% ^^ Task - activity by learned day binned/controlled by reaction time

plot_roi = 6;

% Set reaction time bins
rxn_bin_edges = [0.1,0.2,0.3,0.5,1];
rxn_bin_centers = rxn_bin_edges(1:end-1) + diff(rxn_bin_edges)./2;
n_rxn_bins = length(rxn_bin_edges)-1;
rxn_bins = discretize(move_t,rxn_bin_edges);

% Get indicies for averaging (learned stage x t x roi x animal)
[learned_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);
[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);
[rxn_bins_idx,~,~] = ...
    ndgrid(rxn_bins,1:length(t),1:n_rois);

act_accum_stage_idx = cat(4,learned_day_idx,t_idx,rxn_bins_idx,roi_idx,animal_idx);
wheel_accum_stage_idx = squeeze(act_accum_stage_idx(:,:,1,[1,2,3,5]));

% Average ROI activity and wheel movement
use_trials = ~any(isnan(reshape(act_accum_stage_idx,sum(trials_recording),[])),2);

act_learnday_rxnbin_avg = accumarray( ...
    reshape(act_accum_stage_idx(use_trials,:,:,:),[],size(act_accum_stage_idx,4)), ...
    reshape(fluor_roi_deconv(use_trials,:),[],1), ...
    [length(learned_day_unique),length(t),n_rxn_bins,n_rois,length(animals)], ...
    @nanmean,NaN(class(fluor_roi_deconv)));

wheel_learnday_rxnbin_avg = accumarray( ...
    reshape(wheel_accum_stage_idx(use_trials,:,:),[],size(wheel_accum_stage_idx,3)), ...
    reshape(wheel_allcat(use_trials,:),[],1), ...
    [length(learned_day_unique),length(t),n_rxn_bins,length(animals)], ...
    @nanmean,NaN(class(wheel_allcat)));

% Get hemisphere difference
roi_hemiflip = circshift(1:n_rois,n_rois/2);
act_hemidiff_learnday_rxnbin_avg = act_learnday_rxnbin_avg - ....
    act_learnday_rxnbin_avg(:,:,:,roi_hemiflip,:);

% Get baseline (pre-learning average)
act_hemidiff_learnday_rxnbin_avg_prelearnsub = act_hemidiff_learnday_rxnbin_avg - ...
    nanmean(act_hemidiff_learnday_rxnbin_avg(learned_day_unique < 0,:,:,:,:),1);

% Plot wheel and activity by reaction time bin
plot_learned_days = -2:1;

figure;
rxn_col = copper(length(rxn_bin_centers));
h = tiledlayout(3,length(plot_learned_days),'TileIndexing','ColumnMajor');
for curr_learned_day = plot_learned_days

    curr_learned_day_idx = learned_day_unique == curr_learned_day;

    % (flip wheel to plot leftwards up)
    curr_wheel = -squeeze(wheel_learnday_rxnbin_avg( ...
        curr_learned_day_idx,:,:,:));

    curr_act = squeeze(act_hemidiff_learnday_rxnbin_avg( ...
        curr_learned_day_idx,:,:,plot_roi,:));

    curr_act_prelearnsub = squeeze(act_hemidiff_learnday_rxnbin_avg_prelearnsub( ...
        curr_learned_day_idx,:,:,plot_roi,:));

    nexttile(h); hold on; set(gca,'ColorOrder',rxn_col);
    plot(t,nanmean(curr_wheel,3),'linewidth',2)
    arrayfun(@(x) line(repmat(rxn_bin_centers(x),2,1), ...
        ylim,'color',rxn_col(x,:)),1:length(rxn_bin_centers));
    line([0,0],ylim,'color','k','linestyle','--');
    ylabel('Wheel');
    title(sprintf('Learned day %d', curr_learned_day));

    nexttile(h); hold on; set(gca,'ColorOrder',rxn_col);
    plot(t,nanmean(curr_act,3)','linewidth',2)
    arrayfun(@(x) line(repmat(rxn_bin_centers(x),2,1), ...
        ylim,'color',rxn_col(x,:)),1:length(rxn_bin_centers));
    line([0,0],ylim,'color','k','linestyle','--');
    ylabel(sprintf('%s_{L-R}',wf_roi(plot_roi).area(1:end-2)));
    title(sprintf('Learned day %d', curr_learned_day));

    nexttile(h); hold on; set(gca,'ColorOrder',rxn_col);
    plot(t,nanmean(curr_act_prelearnsub,3)','linewidth',2)
    arrayfun(@(x) line(repmat(rxn_bin_centers(x),2,1), ...
        ylim,'color',rxn_col(x,:)),1:length(rxn_bin_centers));
    line([0,0],ylim,'color','k','linestyle','--');
    ylabel(sprintf('%s_{L-R} -prelearn avg',wf_roi(plot_roi).area(1:end-2)));
    title(sprintf('Learned day %d', curr_learned_day));

end
% (link wheel and activity y-axes)
ax = reshape(flipud(allchild(h)),3,[]);
linkaxes(ax(1,:),'y');
linkaxes(ax(2:3,:),'y');

% Get and plot average tmax activity across rxn bins
use_t = t >= 0 & t <= 0.2;

act_hemidiff_learnday_prelearnsub_tmax = ...
    squeeze(nanmean(max(act_hemidiff_learnday_rxnbin_avg_prelearnsub(:,use_t,:,:,:),[],2),3));

act_hemidiff_learnday_tmax = ...
    squeeze(nanmean(max(act_hemidiff_learnday_rxnbin_avg(:,use_t,:,:,:),[],2),3));
act_hemidiff_learnday_tmax_prelearnsub = ...
    act_hemidiff_learnday_tmax - ...
    nanmean(act_hemidiff_learnday_tmax(learned_day_unique < 0,:,:),1);

figure; hold on
errorbar(learned_day_unique, ...
    nanmean(act_hemidiff_learnday_prelearnsub_tmax(:,plot_roi,:),3), ...
    AP_sem(act_hemidiff_learnday_prelearnsub_tmax(:,plot_roi,:),3),'linewidth',2);

% (plot image to pull out data to combine with passive);
figure;
imagesc([],learned_day_unique,squeeze(act_hemidiff_learnday_prelearnsub_tmax(:,plot_roi,:)));
xlabel('Animal');
ylabel('Learned day');

%% ^^ Task - [as above, but minus move hemiratio]

plot_roi = 6;

% Set reaction time bins
rxn_bin_edges = [0.1,0.2,0.5,Inf];%[0.1,0.25,0.5,1];
rxn_bin_centers = rxn_bin_edges(1:end-1) + diff(rxn_bin_edges)./2;
n_rxn_bins = length(rxn_bin_edges)-1;
rxn_bins = discretize(move_t,rxn_bin_edges);

% Get indicies for averaging (learned stage x t x roi x animal)
[learned_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);
[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);
[rxn_bins_idx,~,~] = ...
    ndgrid(rxn_bins,1:length(t),1:n_rois);

act_accum_stage_idx = cat(4,learned_day_idx,t_idx,rxn_bins_idx,roi_idx,animal_idx);
wheel_accum_stage_idx = squeeze(act_accum_stage_idx(:,:,1,[1,2,3,5]));

% Average ROI activity and wheel movement
use_trials = ~any(isnan(reshape(act_accum_stage_idx,sum(trials_recording),[])),2);

% (MOVE HEMIRATIO)
curr_act = fluor_roi_deconv - ...
    trial_move_hemiratio.*fluor_roi_deconv(:,:,roi_hemiflip);

act_learnday_rxnbin_avg = accumarray( ...
    reshape(act_accum_stage_idx(use_trials,:,:,:),[],size(act_accum_stage_idx,4)), ...
    reshape(curr_act(use_trials,:),[],1), ...
    [length(learned_day_unique),length(t),n_rxn_bins,n_rois,length(animals)], ...
    @nanmean,NaN(class(fluor_roi_deconv)));

wheel_learnday_rxnbin_avg = accumarray( ...
    reshape(wheel_accum_stage_idx(use_trials,:,:),[],size(wheel_accum_stage_idx,3)), ...
    reshape(wheel_allcat(use_trials,:),[],1), ...
    [length(learned_day_unique),length(t),n_rxn_bins,length(animals)], ...
    @nanmean,NaN(class(wheel_allcat)));

% Get hemisphere difference
roi_hemiflip = circshift(1:n_rois,n_rois/2);
act_hemidiff_learnday_rxnbin_avg = act_learnday_rxnbin_avg;

% Get baseline (pre-learning average)
act_hemidiff_learnday_rxnbin_avg_prelearnsub = act_hemidiff_learnday_rxnbin_avg - ...
    nanmean(act_hemidiff_learnday_rxnbin_avg(learned_day_unique < 0,:,:,:,:),1);

% Plot wheel and activity by reaction time bin
plot_learned_days = -2:1;

figure;
rxn_col = copper(length(rxn_bin_centers));
h = tiledlayout(3,length(plot_learned_days),'TileIndexing','ColumnMajor');
for curr_learned_day = plot_learned_days

    curr_learned_day_idx = learned_day_unique == curr_learned_day;

    % (flip wheel to plot leftwards up)
    curr_wheel = -permute(wheel_learnday_rxnbin_avg( ...
        curr_learned_day_idx,:,:,:),[2,3,4,1]);

    curr_act = permute(act_hemidiff_learnday_rxnbin_avg( ...
        curr_learned_day_idx,:,:,plot_roi,:),[2,3,5,1,4]);

    curr_act_prelearnsub = permute(act_hemidiff_learnday_rxnbin_avg_prelearnsub( ...
        curr_learned_day_idx,:,:,plot_roi,:),[2,3,5,1,4]);

    nexttile(h); hold on; set(gca,'ColorOrder',rxn_col);
    plot(t,nanmean(curr_wheel,3),'linewidth',2)
    arrayfun(@(x) line(repmat(rxn_bin_centers(x),2,1), ...
        ylim,'color',rxn_col(x,:)),1:length(rxn_bin_centers));
    line([0,0],ylim,'color','k','linestyle','--');
    ylabel('Wheel');
    title(sprintf('Learned day %d', curr_learned_day));

    nexttile(h); hold on; set(gca,'ColorOrder',rxn_col);
    plot(t,nanmean(curr_act,3)','linewidth',2)
    arrayfun(@(x) line(repmat(rxn_bin_centers(x),2,1), ...
        ylim,'color',rxn_col(x,:)),1:length(rxn_bin_centers));
    line([0,0],ylim,'color','k','linestyle','--');
    ylabel(sprintf('%s_{L-R}',wf_roi(plot_roi).area(1:end-2)));
    title(sprintf('Learned day %d', curr_learned_day));

    nexttile(h); hold on; set(gca,'ColorOrder',rxn_col);
    plot(t,nanmean(curr_act_prelearnsub,3)','linewidth',2)
    arrayfun(@(x) line(repmat(rxn_bin_centers(x),2,1), ...
        ylim,'color',rxn_col(x,:)),1:length(rxn_bin_centers));
    line([0,0],ylim,'color','k','linestyle','--');
    ylabel(sprintf('%s_{L-R} -prelearn avg',wf_roi(plot_roi).area(1:end-2)));
    title(sprintf('Learned day %d', curr_learned_day));

end
% (link wheel and activity y-axes)
ax = reshape(flipud(allchild(h)),3,[]);
linkaxes(ax(1,:),'y');
linkaxes(ax(2:3,:),'y');

% Get and plot average tmax activity across rxn bins
use_t = t >= 0 & t <= 0.2;

act_hemidiff_learnday_prelearnsub_tmax = ...
    squeeze(nanmean(max(act_hemidiff_learnday_rxnbin_avg_prelearnsub(:,use_t,:,:,:),[],2),3));

figure; hold on
errorbar(learned_day_unique, ...
    nanmean(act_hemidiff_learnday_prelearnsub_tmax(:,plot_roi,:),3), ...
    AP_sem(act_hemidiff_learnday_prelearnsub_tmax(:,plot_roi,:),3),'linewidth',2);

% (plot image to pull out data to combine with passive);
figure;
imagesc([],learned_day_unique,squeeze(act_hemidiff_learnday_prelearnsub_tmax(:,plot_roi,:)));
xlabel('Animal');
ylabel('Learned day');


%% Task & Passive - combine daysplit [requires above plots made first]

% First make plots for passive and task daysplit mPFC activity

% (click task L-R)
xt = get(gco,'XData');yt = get(gco,'YData');et = get(gco,'UData');

% (click passive L hemi R stim)
xp = get(gco,'XData');ypl = get(gco,'YData');epl = get(gco,'UData');
% (click passive R hemi L stim)
xp = get(gco,'XData');ypr = get(gco,'YData');epr = get(gco,'UData');


% (click passive L hemi R stim - naive)
xpn = get(gco,'XData');ypln = get(gco,'YData');epln = get(gco,'UData');
% (click passive R hemi L stim - naive)
xpn = get(gco,'XData');yprn = get(gco,'YData');eprn = get(gco,'UData');


% Get task daysplit
n_daysplit = mode(diff(find(isnan(yt)))) - 1;
tp_daysplit_offset = linspace(0,1,n_daysplit+2);

% Re-distribute x-values to make room for passive
xtr = xt(1:n_daysplit+1:end) + tp_daysplit_offset';
ytr = padarray(reshape(yt,n_daysplit+1,[]),[1,0],NaN,'post');
etr = padarray(reshape(et,n_daysplit+1,[]),[1,0],NaN,'post');

xtp = xtr(end-2:end,:);

figure; hold on;

% (plot task redistributed with an extra space)
yyaxis left;
ht = errorbar(xtr(:),ytr(:),etr(:),'color','k','linewidth',2,'CapSize',0);

axis tight;
xlim(xlim + [-0.5,0.5]);
ylabel('\DeltaF/F_0');

% (plot passive)
yyaxis right;
hpl = errorbar(xtp(2,:),ypl,epl,'.','MarkerSize',20,'color','r','linestyle','none','linewidth',2,'CapSize',0);
hpr = errorbar(xtp(2,:),ypr,epr,'.','MarkerSize',20,'color','b','linestyle','none','linewidth',2,'CapSize',0);

% (plot passive naive)
hpl = errorbar(xpn,ypln,epln,'.','MarkerSize',20,'color','r','linestyle','none','linewidth',2,'CapSize',0);
hpr = errorbar(xpn,yprn,eprn,'.','MarkerSize',20,'color','b','linestyle','none','linewidth',2,'CapSize',0);

axis tight;
xlim(xlim + [-0.5,0.5]);
xline(0,'linestyle','--');
xlabel('Learned day');
ylabel('\DeltaF/F_0');

legend([ht,hpl,hpr],{'Task','Passive (L-R/R stim)','Passive (R-L/L stim)'},'location','nw');


%% Task & Passive - combine whole day [requires above plots made first]

% (click task L-R)
xt = get(gco,'XData');yt = get(gco,'YData');et = get(gco,'UData');

% (click passive L-R)
xp = get(gco,'XData');yp = get(gco,'YData');ep = get(gco,'UData');


% Subtract baseline (learning day < 0)
yt = yt - nanmean(yt(xt < 0));
yp = yp - nanmean(yp(xp < 0));

% Divide by post-learn
yt_div = nanmean(yt(xt > 0));
yp_div = nanmean(yp(xp > 0));

yt = yt./yt_div;
et = et./yt_div;
yp = yp./yp_div;
ep = ep./yp_div;

% Get shared days to connect
xtp_shared = intersect(xt,xp);
ytp = [yt(ismember(xt,xtp_shared));yp(ismember(xp,xtp_shared))];

% Plot task and passive together
x_offset_passive = 0.2;

figure; hold on;
plot( ...
    reshape(padarray(xtp_shared + [0;x_offset_passive],[1,0],NaN,'post'),[],1), ...
    reshape(padarray(ytp,[1,0],NaN,'post'),[],1),'k','linewidth',2);
p1 = errorbar(xt,yt,et,'.k','MarkerSize',20,'linewidth',2,'CapSize',0);
p2 = errorbar(xp+0.2,yp,ep,'.r','MarkerSize',20,'linewidth',2,'CapSize',0);
xlabel('Learned day');
ylabel('\DeltaF/F_0');
xline(0,'linestyle','--');
legend([p1,p2],{'Task','Passive'},'location','nw');


%% Task & Passive - combine whole day by animal [requires image plot]

% (click task)
xt = get(gco,'YData');
yt = get(gco,'CData');

% (click passive)
xp = get(gco,'YData');
yp = get(gco,'CData');


% (z-score?)
yt = (yt-nanmean(yt,1))./nanstd(yt);
yp = (yp-nanmean(yp,1))./nanstd(yp);


xtp = intersect(xt,xp);
ytp = cat(3,yt(ismember(xt,xtp),:),yp(ismember(xp,xtp),:));

xtp_long = reshape(xtp+linspace(0,1,3)',[],1);
ytp_long = reshape(permute(padarray(ytp,[0,0,1],NaN,'post'),[3,1,2]),[],length(animals));

figure; hold on;
animal_col = max(0,jet(length(animals))-0.2);
set(gca,'ColorOrder',animal_col);
plot(xtp_long,ytp_long);

plot(xtp_long,nanmean(ytp_long,2),'k','linewidth',2);




%% Task - ephys

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
    kernel./(mua_norm),kernel_set,'uni',false), ...
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
        xlabel('Time from event');
        ylabel('\DeltaR/R_0');
    end
end
legend(mua_areas(plot_areas));
linkaxes(allchild(h),'y');


%% Passive - ephys

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
use_trials = quiescent_trials & trial_stim_allcat == 1; 
[animal_idx,t_idx,area_idx] = ...
    ndgrid(trial_animal(use_trials),1:length(t),1:length(mua_areas));
mua_animal_avg = accumarray([animal_idx(:),t_idx(:),area_idx(:)], ...
    reshape(mua_area_allcat(use_trials,:,:),[],1), ...
    [length(animals),length(t),length(mua_areas)], ...
    @nanmean,NaN);

figure;
plot_areas = area_recording_n == length(trials_recording);
h = AP_errorfill(t', ...
    squeeze(nanmean(mua_animal_avg(:,:,plot_areas),1)), ...
    squeeze(AP_sem(mua_animal_avg(:,:,plot_areas),1)));
xline(0,'color','k');xline(0.5,'color','k');
xlabel('Time from stim (s)');
ylabel('\DeltaR/R_0');
legend(h,mua_areas(plot_areas));

% Plot animal mean with errorbars on stacked subplots
figure;
h = tiledlayout(sum(plot_areas),1,'TileSpacing','compact','padding','compact');
for curr_area = find(plot_areas)'
    nexttile;
    AP_errorfill(t', ...
        squeeze(nanmean(mua_animal_avg(:,:,curr_area),1)), ...
        squeeze(AP_sem(mua_animal_avg(:,:,curr_area),1)),'k');
    ylabel(mua_areas(curr_area));
end
linkaxes(allchild(h),'xy');
% (shade stim area and put in back)
arrayfun(@(x) patch(x,[0,0.5,0.5,0], ...
    reshape(repmat(ylim(x),2,1),[],1),[1,1,0.8], ...
    'linestyle','none'),allchild(h));
arrayfun(@(x) set(x,'children',circshift(get(x,'children'),-1)),allchild(h));


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

stim_px_avg_stage = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),squeeze(nanmean(stim_v_avg_stage(:,:,:,3,:),5)));
AP_imscroll(reshape(permute(stim_px_avg_stage,[1,2,4,3]),size(U_master,1),[],length(t)),t);
axis image off;
colormap(AP_colormap('KWG'));
caxis([-max(abs(caxis)),max(abs(caxis))]);

use_t = t >= 0.05 & t <= 0.2;
stim_px_avg_stage_tmax = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    squeeze(nanmean(nanmean(stim_v_avg_stage(:,use_t,:,:,:),2),5)));

figure;
tiledlayout(2,3,'TileSpacing','compact','padding','compact');
c = (max(stim_px_avg_stage_tmax(:)).*[-1,1])*0.5;
for curr_stage = 1:2
    for curr_stim = 1:3
        nexttile;
        imagesc(stim_px_avg_stage_tmax(:,:,curr_stage,curr_stim));
        axis image off;
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        colormap(AP_colormap('KWG'));
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


%% Ephys - plot probe position

animals = {'AP100','AP101','AP104','AP105','AP106'};

% Load CCF annotated volume
allen_atlas_path = fileparts(which('template_volume_10um.npy'));
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']);

figure;
animal_col = AP_colormap('KR',5);

% Set up 3D axes
ccf_3d_axes = subplot(1,4,1);
[~, brain_outline] = plotBrainGrid([],ccf_3d_axes);
set(ccf_3d_axes,'ZDir','reverse');
hold(ccf_3d_axes,'on');
axis vis3d equal off manual
view([-30,25]);
axis tight;
h = rotate3d(ccf_3d_axes);
h.Enable = 'on';

% Set up 2D axes
bregma_ccf = [540,44,570];
ccf_size = size(av);

ccf_axes = gobjects(3,1);
ccf_axes(1) = subplot(1,4,2,'YDir','reverse');
hold on; axis image off;
ccf_axes(2) = subplot(1,4,3,'YDir','reverse');
hold on; axis image off;
ccf_axes(3) = subplot(1,4,4,'YDir','reverse');
hold on; axis image off;
for curr_view = 1:3
    curr_outline = bwboundaries(squeeze((max(av,[],curr_view)) > 1));    
    cellfun(@(x) plot(ccf_axes(curr_view),x(:,2),x(:,1),'k','linewidth',2),curr_outline)
    
    curr_bregma = fliplr(bregma_ccf(setdiff(1:3,curr_view)));
    plot(ccf_axes(curr_view),curr_bregma(1),curr_bregma(2),'rx');
    
    curr_size = fliplr(ccf_size(setdiff(1:3,curr_view)));
    xlim(ccf_axes(curr_view),[0,curr_size(1)]);
    ylim(ccf_axes(curr_view),[0,curr_size(2)]);
end
linkaxes(ccf_axes);

probe_coords_mean_all = nan(length(animals),3);
for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    
    % Load animal probe histology
    [probe_ccf_fn,probe_ccf_fn_exists] = AP_cortexlab_filename(animal,[],[],'probe_ccf');
    load(probe_ccf_fn);
        
    % Get line of best fit through mean of marked points
    probe_coords_mean = mean(probe_ccf.points,1);
    % (store mean for plotting later)
    probe_coords_mean_all(curr_animal,:) = probe_coords_mean;
    xyz = bsxfun(@minus,probe_ccf.points,probe_coords_mean);
    [~,~,V] = svd(xyz,0);
    histology_probe_direction = V(:,1);
    
    % (make sure the direction goes down in DV - flip if it's going up)
    if histology_probe_direction(2) < 0
        histology_probe_direction = -histology_probe_direction;
    end
    
    % Evaluate line of best fit (length of probe to deepest point)
    [~,deepest_probe_idx] = max(probe_ccf.points(:,2));
    probe_deepest_point = probe_ccf.points(deepest_probe_idx,:);
    probe_deepest_point_com_dist = pdist2(probe_coords_mean,probe_deepest_point);
    probe_length_ccf = 3840/10; % mm / ccf voxel size
    
    probe_line_eval = probe_deepest_point_com_dist - [probe_length_ccf,0];
    probe_line = (probe_line_eval'.*histology_probe_direction') + probe_coords_mean;
    
    % Draw probe in 3D view
    line(ccf_3d_axes,probe_line(:,1),probe_line(:,3),probe_line(:,2), ...
        'linewidth',2,'color',animal_col(curr_animal,:))
    
    % Draw probes on 2D views
    line(ccf_axes(1),probe_line(:,3),probe_line(:,2),'linewidth',2,'color',animal_col(curr_animal,:));
    line(ccf_axes(2),probe_line(:,3),probe_line(:,1),'linewidth',2,'color',animal_col(curr_animal,:));
    line(ccf_axes(3),probe_line(:,2),probe_line(:,1),'linewidth',2,'color',animal_col(curr_animal,:));
    
    drawnow
end

% Plot scalebar
scalebar_length = 1000/10; % mm*voxel size
line(ccf_axes(3),[0,0],[0,scalebar_length],'color','m','linewidth',3);

% Make legend
h = get(ccf_axes(end),'Children');
legend(h(2:end-2),animals);


% Plot probe positions over widefield stim response

% Load CCF to widefield transform
um2pixel = 20.6;
alignment_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\widefield_alignment';
ccf_tform_fn = [alignment_path filesep 'ccf_tform.mat'];
load(ccf_tform_fn);

% Transform probe coordinates from CCF to widefield
probe_coords_mean_wf = ...
    [(10/um2pixel)*probe_coords_mean_all(:,[3,1]), ...
    ones(length(animals),1)]* ...
    ccf_tform.T;

% Plot probe coordinates over widefield image (passive, saved above)
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
stim_px_fn = fullfile(data_path,'stim_px');
load(stim_px_fn);

figure;
imagesc(stim_px);
axis image off
caxis(max(caxis)*0.5*[-1,1])
colormap(AP_colormap('KWG',[],1.5));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
AP_reference_outline('grid_aligned',[0.8,0.8,0.8]);

scatter(probe_coords_mean_wf(:,1),probe_coords_mean_wf(:,2),100,animal_col,'filled');


%% Plot widefield ROIs (hemisphere)

wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

plot_rois = [1,6,7];

roi_col = repmat([0.5,0.7,0.5],n_rois,1);
roi_cat = cat(3,wf_roi.mask);

figure; hold on
set(gca,'YDir','reverse');
AP_reference_outline('ccf_aligned_lefthemi',[0.7,0.7,0.7]);
for curr_roi = plot_rois
    curr_roi_boundary = cell2mat(bwboundaries(roi_cat(:,:,curr_roi)));
    patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),roi_col(curr_roi,:), ...
        'EdgeColor','none');
end
axis image off;




