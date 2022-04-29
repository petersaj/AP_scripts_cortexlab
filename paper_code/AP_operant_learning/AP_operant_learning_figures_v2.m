%% Operant learning widefield paper figures
% Preprocessing done in AP_operant_learning_preprocessing
%
% (mostly just subset of what's in v1) 


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


%% Behavior
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

% Concatenate and get indicies
rxn_measured_allcat = cell2mat(cellfun(@cell2mat,rxn_measured,'uni',false));
rxn_alt_allcat = cell2mat(cellfun(@cell2mat,rxn_alt,'uni',false));

animal_idx = cell2mat(cellfun(@(x,animal) repmat(animal,length(cell2mat(x)),1), ...
    rxn_measured,num2cell(1:length(bhv))','uni',false));
day_idx = cell2mat(cellfun(@(x) cell2mat(cellfun(@(x,day) ...
    repmat(day,length(x),1),x,num2cell(1:length(x))','uni',false)), ...
    rxn_measured,'uni',false));


% Set bins for reaction time histograms
rxn_bins = [0:0.01:0.5];
rxn_bin_centers = rxn_bins(1:end-1) + diff(rxn_bins)./2;

% Get rxn histograms by animal/day and plot heatmap
animal_rxn = nan(max_days,length(rxn_bin_centers),length(animals));
for curr_animal = 1:length(bhv)
    for curr_day = find(use_days{curr_animal})'
        animal_rxn(curr_day,:,curr_animal,:) = ...
            histcounts(bhv(curr_animal).stim_move_t{curr_day}, ...
            rxn_bins,'normalization','probability');
    end
end

figure;
imagesc([],rxn_bin_centers,nanmean(animal_rxn(plot_days,:,:),3)');
colormap(1-gray);
h = colorbar;ylabel(h,'Probability');
xlabel('Day');
ylabel('Reaction time');

% Get histogram of reaction times in day groups (across animals)
day_grps = [1,3,5,7,Inf];
day_grp = discretize(day_idx,day_grps);

figure;
tiledlayout(length(day_grps)-1,1);
for curr_day_grp = 1:length(day_grps)-1

    animal_rxn_measured_cathist = cell2mat(arrayfun(@(x) ...
        histcounts(rxn_measured_allcat( ...
        animal_idx == x & ...
        day_grp == curr_day_grp), ...
        rxn_bins,'normalization','probability')',1:length(bhv),'uni',false));

    animal_rxn_alt_cathist = cell2mat(permute( ...
        arrayfun(@(rep) cell2mat(arrayfun(@(x) ...
        histcounts(rxn_alt_allcat( ...
        animal_idx == x & ...
        day_grp == curr_day_grp,rep), ...
        rxn_bins,'normalization','probability')',1:length(bhv),'uni',false)), ...
        1:n_rxn_altsample,'uni',false),[1,3,2]));

    animal_rxn_alt_cathist_ci = ...
        squeeze(prctile(nanmean(animal_rxn_alt_cathist,2),[5,95],3));

    nexttile; hold on;
    
    AP_errorfill(rxn_bin_centers,nanmean(nanmean(animal_rxn_alt_cathist,2),3), ...
        animal_rxn_alt_cathist_ci,[0.5,0.5,0.5],[],false);
    plot(rxn_bin_centers,nanmean(animal_rxn_measured_cathist,2),'k','linewidth',2);
    xlabel('Reaction time');
    ylabel('Frequency');
    title(sprintf('Day %d-%d',day_grps(curr_day_grp),day_grps(curr_day_grp+1)-1));

end

% Reaction time median: whole day
% (exclude too-fast rxn < 0.1)
rxn_measured_med = accumarray([day_idx,animal_idx], ...
    rxn_measured_allcat.*AP_nanout(rxn_measured_allcat < 0.1), ...
    [max_days,length(bhv)],@(x) nanmedian(x),NaN);
rxn_alt_med = cell2mat(permute(arrayfun(@(x) ...
    accumarray([day_idx,animal_idx], ...
    rxn_alt_allcat(:,x).*AP_nanout(rxn_alt_allcat(:,x) < 0.1), ...
    [max_days,length(bhv)],@(x) nanmedian(x),NaN), ...
    1:n_rxn_altsample,'uni',false),[1,3,2]));

% (plot mice separately)
[~,learn_sort_idx] = sort(learned_day,'ascend');
figure; h = tiledlayout('flow');
for curr_animal = 1:length(animals)
    nexttile; hold on; set(gca,'YScale','log');
    plot(rxn_measured_med(:,curr_animal),'k','linewidth',2);
    AP_errorfill([],[],prctile(rxn_alt_med(:,curr_animal,:),[5,95],3),[0.5,0.5,0.5]);
    xline(learned_day(curr_animal),'color','r');
    set(gca,'children',circshift(get(gca,'children'),-1))
    title(animals{curr_animal});
    ylabel('Reaction time');
    xlabel('Training day');
end

% (plot average)
figure; 
subplot(1,2,1,'YScale','log');hold on
rxn_alt_med_ci = squeeze(prctile(nanmedian(rxn_alt_med,2),[5,95],3));
AP_errorfill([],[],rxn_alt_med_ci(plot_days,:),[0.5,0.5,0.5],[],false);
errorbar(nanmedian(rxn_measured_med(plot_days,:),2), ...
    mad(rxn_measured_med(plot_days,:),1,2),'k','CapSize',0,'linewidth',2)
xlabel('Training day');
ylabel('Median reaction time (s)');
axis tight;
xlim(xlim + [-0.5,0.5]);

subplot(1,2,2);hold on
rxn_measured_med_altdiff = rxn_measured_med - nanmedian(rxn_alt_med,3);
rxn_alt_med_altdiff = rxn_alt_med - nanmedian(rxn_alt_med,3);

rxn_alt_med_altdiff_ci = squeeze(prctile(nanmedian(rxn_alt_med_altdiff,2),[5,95],3));
AP_errorfill([],[],rxn_alt_med_altdiff_ci(plot_days,:),[0.5,0.5,0.5],[],false);
errorbar(nanmedian(rxn_measured_med_altdiff(plot_days,:),2), ...
    mad(rxn_measured_med_altdiff(plot_days,:),1,2),'k','CapSize',0,'linewidth',2)
xlabel('Training day');
ylabel('Median reaction time (meas-null) (s)');
axis tight;
xlim(xlim + [-0.5,0.5]);


% Reaction time median: daysplit
% (exclude too-fast rxn < 0.1)
n_daysplit = 3;
daysplit_idx = cell2mat(cellfun(@(x) ...
    min(floor(linspace(1,n_daysplit+1,length(x))),n_daysplit)', ...
    cat(1,rxn_measured{:}),'uni',false));

rxn_measured_med_daysplit = accumarray([daysplit_idx,day_idx,animal_idx], ...
    rxn_measured_allcat.*AP_nanout(rxn_measured_allcat < 0.1), ...
    [n_daysplit,max_days,length(bhv)],@(x) nanmedian(x),NaN);

rxn_alt_med_daysplit = cell2mat(permute(arrayfun(@(x) ...
    accumarray([daysplit_idx,day_idx,animal_idx], ...
    rxn_alt_allcat(:,x).*AP_nanout(rxn_alt_allcat(:,x) < 0.1), ...
    [n_daysplit,max_days,length(bhv)],@(x) nanmedian(x),NaN), ...
    1:n_rxn_altsample,'uni',false),[1,3,4,2]));

% Put NaNs between days to plot with gaps
rxn_measured_med_long = reshape(padarray(rxn_measured_med_daysplit,[1,0,0],NaN,'post'),[],length(animals));
rxn_alt_med_long = reshape(padarray(rxn_alt_med_daysplit,[1,0,0],NaN,'post'),[],length(animals),n_rxn_altsample);

% Plot relative to learned day (daysplit)
learned_day_x = [1:max_days]'-learned_day';
learned_daysplit_x = cell2mat(cellfun(@(x) x+(0:n_daysplit)'/(n_daysplit+1), ...
    num2cell(learned_day_x),'uni',false));

[rxn_learn_med_daysplit,rxn_learn_mad_daysplit,learned_day_grp_daysplit,learned_day_n_daysplit] = ...
    grpstats(rxn_measured_med_long(:),learned_daysplit_x(:), ...
    {'nanmedian',@(x) mad(x,1),'gname','numel'});
rxn_alt_learn_med_daysplit = ...
    grpstats(reshape(rxn_alt_med_long,[],n_rxn_altsample),learned_daysplit_x(:), ...
    {'nanmedian'});

learned_day_grp_daysplit = cellfun(@str2num,learned_day_grp_daysplit);
plot_learned = learned_day_n_daysplit >= min_n | isnan(rxn_learn_med_daysplit);

figure; hold on;set(gca,'YScale','log');

rxn_alt_learn_ci = prctile(rxn_alt_learn_med_daysplit,[1,99],2);
p1 = AP_errorfill(learned_day_grp_daysplit(plot_learned), ...
    nanmean(rxn_alt_learn_med_daysplit(plot_learned,:),2), ...
    rxn_alt_learn_ci(plot_learned,:),[0.5,0.5,0.5],[],false);

p2 = errorbar(learned_day_grp_daysplit(plot_learned),rxn_learn_med_daysplit(plot_learned), ...
    rxn_learn_mad_daysplit(plot_learned),'k','linewidth',2,'CapSize',0);
ylabel('Median reaction time')
xlabel('Learned day');
axis tight;
xlim(xlim + [-0.5,0.5]);
line([0,0],ylim,'color','k','linestyle','--');
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
    p1 = AP_errorfill(rxn_bin_centers,[], ...
        rxn_alt_hist_ci{curr_cond}',[0.5,0.5,0.5],[],false);

    p2 = plot(rxn_bin_centers, ...
        nanmean(cat(1,rxn_measured_hist{curr_cond,:}),1)','k','linewidth',2);

    legend({'Null','Measured'});
    xlabel('Reaction time');
    ylabel('Probability')
    title(condition_labels{curr_cond});
end
linkaxes(allchild(h),'xy');

% Plot total velocity
figure; hold on;
errorbar(nanmean(muscimol_v1_wheel_mm,1), ...
    AP_sem(muscimol_v1_wheel_mm,1),'k','linewidth',2,'capsize',0);
xlim(xlim+[-0.5,0.5]);
ylim([0,max(ylim)])
set(gca,'XTick',1:n_conditions,'XTickLabel',condition_labels);
ylabel('Wheel mm/min');

% Get reaction times median (exclude too-fast <0.1)
rxn_measured_med = cellfun(@(x) nanmedian(x.*AP_nanout(x < 0.1),1),rxn_measured_cat);
rxn_alt_med = cell2mat(cellfun(@(x) ...
    reshape(nanmedian(x.*AP_nanout(x < 0.1),1),1,1,[]),rxn_alt_cat,'uni',false));

rxn_measured_med_altdiff = rxn_measured_med - nanmean(rxn_alt_med,3);
rxn_alt_med_altdiff = rxn_alt_med - nanmean(rxn_alt_med,3);

figure; 

subplot(1,2,1,'YScale','log'); hold on;
rxn_alt_med_ci = permute(prctile(nanmean(rxn_alt_med,2),[5,95],3),[1,3,2]);
AP_errorfill([],nanmean(rxn_alt_med_ci,2),rxn_alt_med_ci,[0.5,0.5,0.5],[],false);
errorbar(nanmean(rxn_measured_med,2),AP_sem(rxn_measured_med,2),'k','linewidth',2,'capsize',0);
xlim([0,n_conditions] + 0.5);
set(gca,'XTick',1:n_conditions,'XTickLabel',condition_labels);
ylabel('Median reaction time (s)');

subplot(1,2,2); hold on
rxn_alt_med_altdiff_ci = permute(prctile(nanmean(rxn_alt_med_altdiff,2),[5,95],3),[1,3,2]);
AP_errorfill([],nanmean(rxn_alt_med_altdiff_ci,2),rxn_alt_med_altdiff_ci,[0.5,0.5,0.5],[],false);
errorbar(nanmean(rxn_measured_med_altdiff,2),AP_sem(rxn_measured_med_altdiff,2),'k','linewidth',2,'capsize',0);
xlim([0,n_conditions] + 0.5);
set(gca,'XTick',1:n_conditions,'XTickLabel',condition_labels);
ylabel('Median reaction time (meas-null,s)');



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


% Get movement hemisphere ratio (OLD: kernel-based)
% % (get ROI kernel for movements during iti)
% move_iti_k = cellfun(@(x) cellfun(@(x) permute(x{2}(2,:,:),[3,2,1]), ...
%     x,'uni',false),fluor_taskpred_k_all,'uni',false);
% 
% move_iti_k_roi = cellfun(@(x) cellfun(@(x) ...
%     AP_svd_roi(U_master(:,:,1:n_vs),x,[],[],...
%     cat(3,wf_roi.mask)),x,'uni',false),move_iti_k,'uni',false);
% 
% % (get ratio of L/R hemisphere)
% roi_hemiflip = circshift(1:n_rois,n_rois/2);
% move_iti_k_roi_hemiratio = cellfun(@(x) cell2mat(cellfun(@(x) ...
%     arrayfun(@(roi) x(roi_hemiflip(roi),:)'\x(roi,:)',1:n_rois),x,'uni',false)), ...
%     move_iti_k_roi,'uni',false);
% 
% % (replicate movement hemiratio for all trials)
% trial_move_hemiratio = ...
%     cell2mat(cellfun(@(ratio,trials) permute(repmat(ratio,trials,1),[1,3,2]), ...
%     mat2cell(cat(1,move_iti_k_roi_hemiratio{:}),ones(length(trials_recording),1),n_rois), ...
%     num2cell(trials_recording),'uni',false));

warning('Move hemi-ratio not finalized yet')

% Get movement hemisphere ratio (from average no-stim movement activity)

% (RATIO FROM AVERAGE RATIO)
% % (get ROI activity for movement)
% % (use either: fluor_move_nostim_rewardable/leftward_all)
% fluor_move_nostim_roi = cellfun(@(x) cellfun(@(x) ...
%     AP_svd_roi(U_master(:,:,1:n_vs), ...
%     permute(x,[3,2,1]),[],[],cat(3,wf_roi.mask)), ...
%     x,'uni',false),fluor_move_nostim_rewardable_all,'uni',false);
% 
% roi_hemiflip = circshift(1:n_rois,n_rois/2);
% fluor_move_nostim_roi_hemiratio = cellfun(@(x) cellfun(@(x) ...
%     arrayfun(@(roi) ...
%     x(roi_hemiflip(roi),:)'\x(roi,:)',1:n_rois), ...
%     x,'uni',false),fluor_move_nostim_roi,'uni',false);
% 
% % Replicate movement hemisphere ratio for each trial
% % (day-specific)
% trial_move_hemiratio = cell2mat(cellfun(@(x,tr) repmat(permute(x,[1,3,2]), ...
%     tr,1,1),vertcat(fluor_move_nostim_roi_hemiratio{:}),num2cell(trials_recording), ...
%     'uni',false));
% % % (animal median)
% % trial_move_hemiratio = cell2mat(cellfun(@(x,tr) ...
% %     permute(repmat(nanmedian(cell2mat(x),1),tr,1),[1,3,2]), ...
% %     fluor_move_nostim_roi_hemiratio,num2cell(trials_animal)','uni',false));

% (RATIO FROM AVERAGE ACTIVITY)
fluor_move_nostim_roi = cellfun(@(x) ...
    AP_svd_roi(U_master(:,:,1:n_vs), ...
    permute(nanmean(cell2mat(x),1),[3,2,1]),[],[],cat(3,wf_roi.mask)), ...
    fluor_move_nostim_rewardable_all,'uni',false);

roi_hemiflip = circshift(1:n_rois,n_rois/2);
fluor_move_nostim_roi_hemiratio = cellfun(@(x) ...
    arrayfun(@(roi) ...
    x(roi_hemiflip(roi),:)'\x(roi,:)',1:n_rois), ...
    fluor_move_nostim_roi,'uni',false);

trial_move_hemiratio = cell2mat(cellfun(@(x,tr) ...
    permute(repmat(x,tr,1),[1,3,2]), ...
    fluor_move_nostim_roi_hemiratio,num2cell(trials_animal)','uni',false));

% (RATIO BY PIXEL)
%  Cov(v1, v2)/var(v2)
fluor_move_nostim_rewardable_animalavg_px = ...
    cellfun(@(x) AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    permute(nanmean(cell2mat(x),1),[3,2,1])), ...
    fluor_move_nostim_rewardable_all,'uni',false);

% scaling as cov(v1,v2)/var(v2)
fluor_move_nostim_rewardable_animalavg_px_hemiratio = ...
    cellfun(@(x,x_flip) ...
    reshape(sum(reshape(x-nanmean(x,3),[],size(x,3)).* ...
    reshape(x_flip-nanmean(x_flip,3),[],size(x,3)),2) ./ ...
    sum(reshape(x_flip-nanmean(x_flip,3),[],size(x,3)).* ...
    reshape(x_flip-nanmean(x_flip,3),[],size(x,3)),2), ...
    size(x,1),size(x,2)), ...
    fluor_move_nostim_rewardable_animalavg_px, ...
    cellfun(@AP_reflect_widefield,fluor_move_nostim_rewardable_animalavg_px,'uni',false), ...
    'uni',false);


%% ^^ Task - stim-aligned movie/tmax by pre/post learning (raw and hemidiff)

% Set trials to use
learned_day_grp_edges = [-Inf,0,Inf];
trial_learned_day_grp = discretize(trial_learned_day,learned_day_grp_edges);
n_learn_grps = length(learned_day_grp_edges)-1;


fluor_learning_animal = nan(n_vs,length(t),n_learn_grps, ...
    length(animals),class(fluor_allcat_deconv));
for curr_animal = 1:length(animals)
    for curr_learned_grp = 1:n_learn_grps
        curr_trials = trial_animal == curr_animal & ...
            trial_learned_day_grp == curr_learned_grp & ...
            trial_outcome_allcat == 1;
        fluor_learning_animal(:,:,curr_learned_grp,curr_animal) = ...
            permute(nanmean(fluor_allcat_deconv(curr_trials,:,:),1),[3,2,1]);
    end
end

fluor_learning_animal_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),fluor_learning_animal);

fluor_learning_animal_px_avg = nanmean(fluor_learning_animal_px,5);
fluor_learning_animal_px_hemidiff_avg = nanmean(fluor_learning_animal_px - ...
    cat(5,fluor_move_nostim_rewardable_animalavg_px_hemiratio{:}).* ...
    AP_reflect_widefield(fluor_learning_animal_px),5);
% (white-out non-imaged pixels)
fluor_learning_animal_px_hemidiff_avg(isnan(fluor_learning_animal_px_hemidiff_avg)) = 0;

AP_imscroll(fluor_learning_animal_px_avg,t)
caxis(max(abs(caxis))*[-1,1]);
colormap(AP_colormap('KWG',[],1.5));
axis image off;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

AP_imscroll(fluor_learning_animal_px_hemidiff_avg(:,1:size(U_master,2)/2,:,:),t)
caxis(max(abs(caxis))*[-1,1]);
colormap(AP_colormap('BWR',[],1.5));
axis image off;
AP_reference_outline('ccf_aligned_lefthemi',[0.5,0.5,0.5]);

% Get time window fluorescence
use_t = t >= 0 & t <= 0.2;
figure; tiledlayout(2,2);
c = max(fluor_learning_animal_px_avg(:)).*[-1,1];
c_hemi = max(fluor_learning_animal_px_hemidiff_avg(:)).*[-1,1];
for curr_learn = 1:n_learn_grps

    nexttile;
    imagesc(max(fluor_learning_animal_px_avg(:,:,use_t,curr_learn),[],3));
    colormap(gca,AP_colormap('KWG',[],1.5)); caxis(c);
    axis image off;
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

    nexttile;
    imagesc(max(fluor_learning_animal_px_hemidiff_avg(:,1:size(U_master,2)/2,use_t,curr_learn),[],3));
    colormap(gca,AP_colormap('KWP',[],1.5)); caxis(c_hemi);
    axis image off;
    AP_reference_outline('ccf_aligned_lefthemi',[0.5,0.5,0.5]);

end



%% ^^ Task - event-aligned pixels by pre/post learning (raw and hemidiff)

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
align_v = cell(length(use_align),1);
for curr_align = 1:length(use_align)
    
    % (re-align activity if not all aligned to first frame)
    if all(use_align{curr_align} == 1)
        curr_v_align = fluor_allcat_deconv;
    else
        curr_v_align = nan(size(fluor_allcat_deconv),'single');
        for curr_trial = find(~isnan(use_align{curr_align}))'
            curr_shift_frames = ...
                use_align{curr_align}(curr_trial) + (0:length(t)-1);
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

    % (get time point and store)
    curr_v_align_t = permute( ...
        interp1(t,permute(curr_v_align_avg,[2,1,3,4]),plot_t{curr_align}), ...
        [2,1,3,4]);
    align_v{curr_align} = curr_v_align_t;

end

% Get average pixels (raw and hemidiff)
align_px_avg = cellfun(@(x) ...
    nanmean(AP_svdFrameReconstruct(U_master(:,:,1:n_vs),x),5), ...
    align_v,'uni',false);

align_px_hemidiff_avg = cellfun(@(x) ...
    nanmean(AP_svdFrameReconstruct(U_master(:,:,1:n_vs),x) - ...
    cat(5,fluor_move_nostim_rewardable_animalavg_px_hemiratio{:}).* ...
    AP_reflect_widefield(AP_svdFrameReconstruct(U_master(:,:,1:n_vs),x)),5), ...
    align_v,'uni',false);

% Plot pixels raw
c = prctile(reshape(cat(3,align_px_avg{:}),[],1),100)*[-1,1]*0.8;
figure;
h = tiledlayout(n_learn_grps,length(cell2mat(plot_t)), ...
    'TileSpacing','compact','padding','compact');
for curr_learn_grp = 1:n_learn_grps
    for curr_align = 1:length(use_align)
        for curr_plot_t = 1:length(plot_t{curr_align})
            
            curr_px = align_px_avg{curr_align}(:,:,curr_plot_t,curr_learn_grp);

            nexttile
            imagesc(curr_px);
            AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
            axis image off;
            colormap(gca,AP_colormap('KWG',[],1.5));
            caxis(c);
            title([use_align_labels{curr_align} ': ' num2str(plot_t{curr_align}(curr_plot_t)) ' sec']);

        end
    end
end
linkaxes(allchild(h),'xy');

% Plot pixels hemidiff
c_hemidiff = c*0.5;
figure;
h = tiledlayout(n_learn_grps,length(cell2mat(plot_t)), ...
    'TileSpacing','compact','padding','compact');
for curr_learn_grp = 1:n_learn_grps
    for curr_align = 1:length(use_align)
        for curr_plot_t = 1:length(plot_t{curr_align})
            
            curr_px = align_px_hemidiff_avg{curr_align}(:,1:size(U_master,2)/2,curr_plot_t,curr_learn_grp);
            curr_px(isnan(curr_px)) = 0; % make non-imaged pixels white

            nexttile
            imagesc(curr_px);
            AP_reference_outline('ccf_aligned_lefthemi',[0.5,0.5,0.5]);
            axis image off;
            colormap(gca,AP_colormap('BWR',[],1.5));
            caxis(c_hemidiff);
            title([use_align_labels{curr_align} ': ' num2str(plot_t{curr_align}(curr_plot_t)) ' sec']);

        end
    end
end
linkaxes(allchild(h),'xy');


%% ^^ Task - trial activity

plot_rois = [1,7,6];
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);
n_trial_smooth = 200;

% Set activity
roi_hemiflip = circshift(1:n_rois,n_rois/2);
curr_act = fluor_roi_deconv;

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
        colormap(AP_colormap('KWG',[],1.5));
        c = prctile(reshape(curr_act(:,:,curr_roi),[],1),95).*[-1,1];
        caxis(c);
        xline(0,'color','k');
        plot(move_t(use_trials(sort_idx)),1:length(use_trials),'color',[0.6,0,0.6]);
        xlabel('Time from stim (s)');
        ylabel('Trial (rxn-time sorted)');
        title(sprintf('%s, stage %d',wf_roi(curr_roi).area,curr_stage));
    end    
end


%% ^^ Task - trial activity (hemidiff)

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

%% ^^ Task - average activity (stim/move align pre/post learn)

plot_rois = [1,7,6];
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);

% Set activity
curr_act = fluor_roi_deconv;

% Move-align activity
curr_act_move = nan(size(curr_act),class(curr_act));

move_align = move_idx - sum(t<0);
for curr_trial = find(~isnan(move_align))'
    curr_shift_frames = ...
        move_align(curr_trial) + (0:length(t)-1);
    curr_shift_frames_use = curr_shift_frames > 0 & curr_shift_frames <= length(t);
    
    curr_grab_frames = curr_shift_frames(curr_shift_frames_use);
    curr_fill_frames = find(curr_shift_frames_use,length(t));
    
    curr_act_move(curr_trial,curr_fill_frames,:) = ...
        curr_act(curr_trial,curr_grab_frames,:);
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
    curr_act(:), ...
    [max(trial_learned_stage),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

roi_move_learn_avg = accumarray( ...
    reshape(accum_stage_idx,[],size(accum_stage_idx,4)), ...
    curr_act_move(:), ...
    [max(trial_learned_stage),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

% Plot L-R stim/move, learning overlaid
median_move_t = median(move_t(trial_learned_day>0));
x_limit = [-0.1,median_move_t];

figure;
h = tiledlayout(1,length(plot_rois)*2);
stage_col = [0.5,0.7,0.5;0,0.7,0];
for curr_roi = plot_rois

    curr_data_stim = permute(roi_stim_learn_avg(:,:,curr_roi,:),[2,1,4,3]);
    curr_data_move = permute(roi_move_learn_avg(:,:,curr_roi,:),[2,1,4,3]);

    nexttile; hold on;
    AP_errorfill(t,nanmean(curr_data_stim,3),AP_sem(curr_data_stim,3),stage_col);
    xlabel('Stim');
    ylabel(wf_roi(curr_roi).area(1:end-2));
    xline(0);set(gca,'children',circshift(get(gca,'children'),-1))
    xlim(x_limit);

    nexttile; hold on;
    AP_errorfill(t,nanmean(curr_data_move,3),AP_sem(curr_data_move,3),stage_col);
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

%% ^^ Task - average activity (stim/move align pre/post learn) (overlay L/R)

plot_rois = [1,7,6];
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);

% Set activity
curr_act = fluor_roi_deconv;

% Move-align activity
curr_act_move = nan(size(curr_act),class(curr_act));

move_align = move_idx - sum(t<0);
for curr_trial = find(~isnan(move_align))'
    curr_shift_frames = ...
        move_align(curr_trial) + (0:length(t)-1);
    curr_shift_frames_use = curr_shift_frames > 0 & curr_shift_frames <= length(t);
    
    curr_grab_frames = curr_shift_frames(curr_shift_frames_use);
    curr_fill_frames = find(curr_shift_frames_use,length(t));
    
    curr_act_move(curr_trial,curr_fill_frames,:) = ...
        curr_act(curr_trial,curr_grab_frames,:);
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
    curr_act(:), ...
    [max(trial_learned_stage),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

roi_move_learn_avg = accumarray( ...
    reshape(accum_stage_idx,[],size(accum_stage_idx,4)), ...
    curr_act_move(:), ...
    [max(trial_learned_stage),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

% Plot L-R stim/move, learning overlaid
median_move_t = median(move_t(trial_learned_day>0));
x_limit = [-0.1,median_move_t];

figure;
h = tiledlayout(2,length(plot_rois)*2);
hemi_col = [0.7,0,0;0,0,0.7];
for curr_stage = unique(trial_learned_stage)'
    for curr_roi = plot_rois

        curr_rois = curr_roi + [0,size(wf_roi,1)];

        curr_data_stim = permute(roi_stim_learn_avg(curr_stage,:,curr_rois,:),[2,3,4,1]);
        curr_data_move = permute(roi_move_learn_avg(curr_stage,:,curr_rois,:),[2,3,4,1]);

        nexttile; hold on;
        AP_errorfill(t,nanmean(curr_data_stim,3),AP_sem(curr_data_stim,3),hemi_col);
        xlabel('Stim');
        ylabel(wf_roi(curr_roi).area(1:end-2));
        xline(0);set(gca,'children',circshift(get(gca,'children'),-1))
        xlim(x_limit);

        nexttile; hold on;
        AP_errorfill(t,nanmean(curr_data_move,3),AP_sem(curr_data_move,3),hemi_col);
        xlabel('Move');
        ylabel(wf_roi(curr_roi).area(1:end-2));
        xline(0);set(gca,'children',circshift(get(gca,'children'),-1))
        xlim(x_limit - x_limit(1));

    end
end

% (link all y-axes);
linkaxes(allchild(h),'y');


%% ^^ Task - average activity (stim/move align pre/post learn) (HEMIDIFF)

plot_rois = [1,7,6];
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);

% Set activity
roi_hemiflip = circshift(1:n_rois,n_rois/2);
curr_act = fluor_roi_deconv - ...
    trial_move_hemiratio.*fluor_roi_deconv(:,:,roi_hemiflip);

% Move-align activity
curr_act_move = nan(size(curr_act),class(curr_act));

move_align = move_idx - sum(t<0);
for curr_trial = find(~isnan(move_align))'
    curr_shift_frames = ...
        move_align(curr_trial) + (0:length(t)-1);
    curr_shift_frames_use = curr_shift_frames > 0 & curr_shift_frames <= length(t);
    
    curr_grab_frames = curr_shift_frames(curr_shift_frames_use);
    curr_fill_frames = find(curr_shift_frames_use,length(t));
    
    curr_act_move(curr_trial,curr_fill_frames,:) = ...
        curr_act(curr_trial,curr_grab_frames,:);
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
    curr_act(:), ...
    [max(trial_learned_stage),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

roi_move_learn_avg = accumarray( ...
    reshape(accum_stage_idx,[],size(accum_stage_idx,4)), ...
    curr_act_move(:), ...
    [max(trial_learned_stage),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

% Plot L-R stim/move, learning overlaid
median_move_t = median(move_t(trial_learned_day>0));
x_limit = [-0.1,median_move_t];

figure;
h = tiledlayout(1,length(plot_rois)*2);
stage_col = min(1,[0.4,0,0.4] + [0.5;0]);
for curr_roi = plot_rois

    curr_data_stim = permute(roi_stim_learn_avg(:,:,curr_roi,:),[2,1,4,3]);
    curr_data_move = permute(roi_move_learn_avg(:,:,curr_roi,:),[2,1,4,3]);

    nexttile; hold on;
    AP_errorfill(t,nanmean(curr_data_stim,3),AP_sem(curr_data_stim,3),stage_col);
    xlabel('Stim');
    ylabel(wf_roi(curr_roi).area(1:end-2));
    xline(0);set(gca,'children',circshift(get(gca,'children'),-1))
    xlim(x_limit);

    nexttile; hold on;
    AP_errorfill(t,nanmean(curr_data_move,3),AP_sem(curr_data_move,3),stage_col);
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

% Get days to plot (by minimum n)
min_n = 4;
learned_days_n = sum(any(stim_roi_learnday_avg_tmax,2),3);
plot_learned_days = learned_days_n >= min_n;


plot_roi = 6;
curr_tcourse = nanmean(squeeze(roi_learnday_avg(plot_learned_days,:,plot_roi,:)),3);
curr_tmax = squeeze(stim_roi_learnday_avg_tmax(plot_learned_days,plot_roi,:));

figure; 
subplot(1,2,1);
imagesc(t,learned_day_unique(plot_learned_days),curr_tcourse);
caxis(max(abs(caxis))*[-1,1]);
colormap(AP_colormap('BWR',[],1.5));
xline(0);
yline(-0.5);
xlabel('Time from stim (s)');
ylabel('Learned day');

subplot(1,2,2);
errorbar(learned_day_unique(plot_learned_days),nanmean(curr_tmax,2),AP_sem(curr_tmax,2),'k','linewidth',2);
xline(0);
xlabel('Time from stim (s)')
ylabel(wf_roi(plot_roi).area);



figure; 
line_col = [0.4,0,0.4];
AP_stackplot(curr_tcourse',t,1.5*std(curr_tcourse(:)),[],line_col,learned_day_unique(plot_learned_days));



figure; 
line_col = [0.4,0,0.4];
AP_stackplot(flipud(curr_tcourse)',t,1.5*std(curr_tcourse(:)),[],line_col,fliplr(learned_day_unique(plot_learned_days)),true);



figure;
curr_roi = 6;
line_days = -3:2;
use_learned_days = ismember(learned_day_unique,line_days);
col = AP_colormap('KWP',1+sum(use_learned_days));
col(median(1:size(col,1)),:) = [];
p = AP_errorfill(t, ...
    squeeze(nanmean(roi_learnday_avg(use_learned_days,:,curr_roi,:),4))', ...
    squeeze(AP_sem(roi_learnday_avg(use_learned_days,:,curr_roi,:),4))', ...
    col);
xline(0);
xlim([-0.05,0.4]);
legend(p,cellfun(@(x) sprintf('LDday %d',x),num2cell(line_days),'uni',false));
ylabel(wf_roi(curr_roi).area);
xlabel('Time from stim (s)');




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

learned_day_unique = unique(trial_learned_day);
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
stim_px_avg_stage = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),stim_v_avg_stage);

use_t = t >= 0 & t <= 0.2;
stim_px_avg_stage_tmax = ...
    squeeze(max(stim_px_avg_stage(:,:,use_t,:,:,:),[],3));

% Plot pixel timeavg
figure;
h = tiledlayout(n_stages,length(stim_unique));
c = max(reshape(nanmean(stim_px_avg_stage_tmax,5),[],1)).*[-1,1]*0.5;
for curr_stage = 1:n_stages
    for curr_stim = 1:length(stim_unique)

        curr_px = nanmean(stim_px_avg_stage_tmax(:,:,curr_stage,curr_stim,:),5);

        nexttile;
        imagesc(curr_px);
        axis image off;
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        colormap(gca,AP_colormap('KWG',[],1.5));
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
plot_rois = reshape([6] + [0,size(wf_roi,1)]',1,[]);
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


% Plot ROIs by stage (L/R overlay)
figure;
plot_rois = 6;

hemi_stage_col = min(1,reshape([0.7,0,0;0,0,0.7]' + cat(3,0.5,0),3,[]))';

h = tiledlayout(length(plot_rois),3,'TileSpacing','compact','padding','compact');
for curr_roi = plot_rois
    for curr_stim = stim_unique'

        curr_rois = curr_roi + [0,size(wf_roi,1)];

        nexttile;
        AP_errorfill(t, ...
            reshape(permute(nanmean(stim_roi_avg_stage(curr_rois,:,:,stim_unique == curr_stim,:),5),[2,1,3]),length(t),[]), ...
            reshape(permute(AP_sem(stim_roi_avg_stage(curr_rois,:,:,stim_unique == curr_stim,:),5),[2,1,3]),length(t),[]), ...
            hemi_stage_col(:,:));
        xlabel('Time from stim (s)');
        ylabel(sprintf('%s \\DeltaF/F_0',wf_roi(curr_roi).area));
        axis tight;xlim([-0.2,1])

    end
end
% Link all axes
linkaxes(allchild(h),'xy');
% (shade stim area and put in back)
arrayfun(@(x) patch(x,[0,0.5,0.5,0], ...
    reshape(repmat(ylim(x),2,1),[],1),[1,1,0.8], ...
    'linestyle','none'),allchild(h));
arrayfun(@(x) set(x,'children',circshift(get(x,'children'),-1)),allchild(h));



%% ^^ Passive - ROI timecourse/tmax by learned day

% Set values for plotting
plot_rois = reshape([1,6]+[0,size(wf_roi,1)]',1,[]);
use_t = t > 0 & t <= 0.2;

% Get indicies for averaging 
[trained_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_day,1:length(t),1:n_rois);

[learned_day_idx,~,~] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);
[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);
[stim_idx,~,~] = ...
    ndgrid(trial_stim_id,1:length(t),1:n_rois);

accum_idx = cat(4,trained_day_idx,learned_day_idx,t_idx,roi_idx,stim_idx,animal_idx);

use_trials_naive = quiescent_trials & trial_day <= n_naive;
use_trials_training = quiescent_trials & trial_day > n_naive;

roi_naive_avg = accumarray( ...
    reshape(accum_idx(use_trials_naive,:,:,[1,3,4,5,6]),[],5), ...
    reshape(fluor_roi_deconv(use_trials_naive,:,:),[],1), ...
    [n_naive,length(t),n_rois,length(stim_unique),length(animals)], ...
    @nanmean,NaN('single'));

roi_learnday_avg = accumarray( ...
    reshape(accum_idx(use_trials_training,:,:,[2,3,4,5,6]),[],5), ...
    reshape(fluor_roi_deconv(use_trials_training,:,:),[],1), ...
    [length(learned_day_unique),length(t),n_rois,length(stim_unique),length(animals)], ...
    @nanmean,NaN('single'));

stim_roi_naive_avg_tmax = squeeze(max(roi_naive_avg(:,use_t,:,:,:),[],2));
stim_roi_learnday_avg_tmax = squeeze(max(roi_learnday_avg(:,use_t,:,:,:),[],2));

% Get days to plot (by minimum n)
min_n = 4;
learned_days_n = sum(any(any(stim_roi_learnday_avg_tmax,2),3),4);
plot_learned_days = learned_days_n >= min_n;

figure;
stim_color = [0,0,0.8;0.5,0.5,0.5;0.8,0,0];
h = tiledlayout(length(plot_rois),6);
for curr_roi = plot_rois   

    c = prctile(nanmean(roi_learnday_avg(plot_learned_days,:,curr_roi,:,:),5),100,'all').*[-1,1];
    for curr_stim = 1:length(stim_unique)
        nexttile;
        imagesc(t,learned_day_unique(plot_learned_days), ...
            nanmean(roi_learnday_avg(plot_learned_days,:,curr_roi,curr_stim,:),5));
        colormap(gca,AP_colormap('KWG',[],1.5));
        caxis(c)
        xline(0); yline(-0.5);
        title(sprintf('Stim %d',stim_unique(curr_stim)));
    end

    nexttile([1,3]); hold on;
    % (naive and training as split lines)
    curr_x = [learned_day_unique(1:3);NaN;learned_day_unique(plot_learned_days)];
    curr_data = squeeze([...
        padarray(stim_roi_naive_avg_tmax(:,curr_roi,:,:),1,NaN,'post'); ...
        stim_roi_learnday_avg_tmax(plot_learned_days,curr_roi,:,:)]);

    set(gca,'ColorOrder',stim_color);
    p = errorbar(repmat(curr_x,1,length(stim_unique)), ...
        nanmean(curr_data,3),AP_sem(curr_data,3),'linewidth',2,'capsize',0);
    xline(0);set(gca,'children',circshift(get(gca,'children'),-1))
    xlabel('Time from stim (s)');
    ylabel(wf_roi(curr_roi).area);
    legend(p,num2str(stim_unique));
    axis tight; xlim(xlim + [-0.5,0.5]);

end

figure;
h = tiledlayout(length(plot_rois),1);
for curr_roi = plot_rois   
    nexttile; hold on;
    plot(t,nanmean(roi_learnday_avg(plot_learned_days & learned_day_unique < 0,:,curr_roi,stim_unique == 1,:),5),'k','linewidth',2)
    plot(t,nanmean(roi_learnday_avg(plot_learned_days & learned_day_unique >= 0,:,curr_roi,stim_unique == 1,:),5),'color',[0,0.7,0],'linewidth',2)
    ylabel(wf_roi(curr_roi).area);
end


figure; 
curr_roi = 6;
curr_data = nanmean(roi_learnday_avg(plot_learned_days,:,curr_roi,stim_unique == 1,:),5);
AP_stackplot(curr_data',t,3*std(curr_data(:)),[],'k',learned_day_unique(plot_learned_days));


figure;
curr_roi = 6;
line_days = -3:2;
use_learned_days = ismember(learned_day_unique,line_days);
col = AP_colormap('KWR',1+sum(use_learned_days));
col(median(1:size(col,1)),:) = [];
p = AP_errorfill(t, ...
    squeeze(nanmean(roi_learnday_avg(use_learned_days,:,curr_roi,stim_unique == 1,:),5))', ...
    squeeze(AP_sem(roi_learnday_avg(use_learned_days,:,curr_roi,stim_unique == 1,:),5))', ...
    col);
xline(0);
xlim([-0.05,0.4]);
legend(p,cellfun(@(x) sprintf('LDday %d',x),num2cell(line_days),'uni',false));
ylabel(wf_roi(curr_roi).area);
xlabel('Time from stim (s)');


%% ^^ Passive - whisker/pupil

% Get average pupil diameter and whisker movement to stimuli
[~,trial_stim_allcat_id] = ismember(trial_stim_allcat,stim_unique);

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

%% Get example pupil video 

animal = 'AP106';
day = '2021-06-25';
experiment = 2;

load_parts.cam = true;
verbose = true;
AP_load_experiment;

AP_mousemovie(eyecam_fn,eyecam_t,eyecam_dlc);


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
stim_px_avg_stage = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),stim_v_avg_stage);

use_t = t >= 0 & t <= 0.2;
stim_px_avg_stage_tmax = ...
    squeeze(max(stim_px_avg_stage(:,:,use_t,:,:,:),[],3));

figure;
tiledlayout(2,3,'TileSpacing','compact','padding','compact');
c = max(reshape(nanmean(stim_px_avg_stage_tmax,5),[],1)).*[-1,1]*1;
for curr_stage = 1:2
    for curr_stim = 1:3
        nexttile;
        imagesc(nanmean(stim_px_avg_stage_tmax(:,:,curr_stage,curr_stim,3),5));
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
h = tiledlayout(length(plot_rois),2);
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
% (link ROI y-axes)
ax = reshape(flipud(allchild(h)),2,length(plot_rois))';
for curr_roi = 1:length(plot_rois)
   linkaxes(ax(:,curr_roi),'y'); 
end
linkaxes(ax,'x');


% Plot ROI time-max pre/post muscimol
use_t = t >= 0 & t <= 0.2;
stim_roi_avg_stage_tmax = squeeze(max(stim_roi_avg_stage(:,use_t,:,:,:),[],2));

plot_rois = [1,6];
plot_stim = 3;
figure; hold on
curr_data = squeeze(stim_roi_avg_stage_tmax(plot_rois,:,plot_stim,:));
plot(squeeze(curr_data(1,:,:)),squeeze(curr_data(2,:,:)),'k');
p1 = plot(squeeze(curr_data(1,2,:)),squeeze(curr_data(2,2,:)),'.k','MarkerSize',20);
p2 = plot(squeeze(curr_data(1,1,:)),squeeze(curr_data(2,1,:)),'.r','MarkerSize',20);
xlabel(wf_roi(plot_rois(1)).area);
ylabel(wf_roi(plot_rois(2)).area);
legend([p1,p2],{'Washout','V1 Muscimol'},'location','nw');


%% Task & Passive - combine daysplit [requires above plots made first]

% First make plots for passive and task daysplit mPFC activity

% (click task L-R)
xt = get(gco,'XData');yt = get(gco,'YData');et = get(gco,'UData');

% (click passive L hemi R stim)
xp = get(gco,'XData');ypl = get(gco,'YData');epl = get(gco,'UData');
% (click passive L hemi R stim)
xp = get(gco,'XData');ypr = get(gco,'YData');epr = get(gco,'UData');


% Get task daysplit
n_daysplit = mode(diff(find(isnan(yt)))) - 1;
tp_daysplit_offset = linspace(0,1,n_daysplit+2);

xp_offset = tp_daysplit_offset(end-1);

% Re-distribute x-values to make room for passive
xtr = xt(1:n_daysplit+1:end) + tp_daysplit_offset';
ytr = padarray(reshape(yt,n_daysplit+1,[]),[1,0],NaN,'post');
etr = padarray(reshape(et,n_daysplit+1,[]),[1,0],NaN,'post');

figure; hold on;

% (plot task redistributed with an extra space)
yyaxis left;
ht = errorbar(xtr(:),ytr(:),etr(:),'color','k','linewidth',2,'CapSize',0);

axis tight;
xlim(xlim + [-0.5,0.5]);
ylabel('\DeltaF/F_0');

% (plot passive)
yyaxis right;
hpl = errorbar(xp+xp_offset,ypl,epl,'.','MarkerSize',20,'color','r','linestyle','none','linewidth',2,'CapSize',0);
hpr = errorbar(xp+xp_offset,ypr,epr,'.','MarkerSize',20,'color','b','linestyle','none','linewidth',2,'CapSize',0);

axis tight;
xlim(xlim + [-0.5,0.5]);
xline(0,'linestyle','--');
xlabel('Learned day');
ylabel('\DeltaF/F_0');

legend([ht,hpl,hpr],{'Task','Passive (R stim)','Passive (L stim)'},'location','nw');


%% Ephys - Passive

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

% Get average response timecourses
stim_unique = unique(trial_stim_allcat);
[~,trial_stim_id] = ismember(trial_stim_allcat,stim_unique);

use_trials = quiescent_trials;

[animal_idx,t_idx,area_idx] = ...
    ndgrid(trial_animal(use_trials),1:length(t),1:length(mua_areas));
[stim_idx,~,~,] = ...
    ndgrid(trial_stim_id(use_trials),1:length(t),1:length(mua_areas));

mua_animal_avg = accumarray([animal_idx(:),t_idx(:),stim_idx(:),area_idx(:)], ...
    reshape(mua_area_allcat(use_trials,:,:),[],1), ...
    [length(animals),length(t),length(stim_unique),length(mua_areas)], ...
    @nanmean,NaN);

% Get DV position of each area for sorting
% (overkill: do it by loading in atlas and getting the volume)
allen_atlas_path = fileparts(which('template_volume_10um.npy'));
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

mua_areas_dvmin = nan(size(mua_areas));
for curr_area = mua_areas'
    curr_structure_idx = find(strcmp(st.safe_name,curr_area));
    curr_structure_id = st.structure_id_path{curr_structure_idx};
    curr_ccf_idx = find(cellfun(@(x) contains(x,curr_structure_id), ...
        st.structure_id_path));

    slice_spacing = 5;
    curr_ccf_volume = ...
        ismember(av(1:slice_spacing:end, ...
        1:slice_spacing:end,1:slice_spacing:end),curr_ccf_idx);

    curr_ccf_coronal_max = permute(max(curr_ccf_volume,[],1),[2,3,1]);
    curr_min_dv = find(any(curr_ccf_coronal_max,2),1);

    mua_areas_dvmin(strcmp(curr_area,mua_areas)) = curr_min_dv;
end

% Plot timecourses overlaid
% (get areas to plot: anything present in all recordings)
[~,area_idx] = cellfun(@(x) ismember(x,mua_areas),mua_areas_cat,'uni',false);
area_recording_n = accumarray(cell2mat(area_idx),1);
plot_areas = find(area_recording_n == length(trials_recording));
[~,plot_area_sort_idx] = sort(mua_areas_dvmin(plot_areas));

figure;
stim_color = [0,0,0.8;0.5,0.5,0.5;0.8,0,0];
h = tiledlayout(length(plot_areas),1);
for curr_area = plot_areas(plot_area_sort_idx)'
    nexttile;
    AP_errorfill(t', ...
        squeeze(nanmean(mua_animal_avg(:,:,:,curr_area),1)), ...
        squeeze(AP_sem(mua_animal_avg(:,:,:,curr_area),1)),stim_color);
    ylabel(mua_areas(curr_area));
end
linkaxes(allchild(h),'xy');
% (shade stim area and put in back)
arrayfun(@(x) patch(x,[0,0.5,0.5,0], ...
    reshape(repmat(ylim(x),2,1),[],1),[1,1,0.8], ...
    'linestyle','none'),allchild(h));
arrayfun(@(x) set(x,'children',circshift(get(x,'children'),-1)),allchild(h));


%% Ephys - plot probe position

animals = {'AP100','AP101','AP104','AP105','AP106'};

% Load CCF annotated volume
allen_atlas_path = fileparts(which('template_volume_10um.npy'));
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

figure;
animal_col = repmat([0,0,0],length(animals),1);

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

% Plot specific areas
plot_structure_names = {'Secondary motor area', ...
    'Anterior cingulate area','Prelimbic area','Infralimbic area'};
plot_structure_colors = lines(length(plot_structure_names));

for plot_structure_name = plot_structure_names
    plot_structure = find(strcmp(st.safe_name,plot_structure_name));

    % Get all areas within and below the selected hierarchy level
    plot_structure_id = st.structure_id_path{plot_structure};
    plot_ccf_idx = find(cellfun(@(x) contains(x,plot_structure_id), ...
        st.structure_id_path));

    % plot the structure
    slice_spacing = 5;
    plot_structure_color = plot_structure_colors( ...
        strcmp(plot_structure_name,plot_structure_names),:);

    % Get structure volume
    plot_ccf_volume = ismember(av(1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end),plot_ccf_idx);

    for curr_view = 1:3
        curr_outline = bwboundaries(squeeze((max(plot_ccf_volume,[],curr_view))));
        cellfun(@(x) plot(ccf_axes(curr_view),x(:,2)*slice_spacing, ...
            x(:,1)*slice_spacing,'color',plot_structure_color,'linewidth',2),curr_outline)
    end
end

% Plot probe locations
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
    
    % Draw probes on coronal + saggital
    line(ccf_axes(1),probe_line(:,3),probe_line(:,2),'linewidth',2,'color',animal_col(curr_animal,:));
    line(ccf_axes(3),probe_line(:,2),probe_line(:,1),'linewidth',2,'color',animal_col(curr_animal,:));

    % Draw probe mean on horizontal
    plot(ccf_axes(2), probe_coords_mean(:,3),probe_coords_mean(:,1), ...
        '.','MarkerSize',10,'color',animal_col(curr_animal,:));
    
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

scatter(probe_coords_mean_wf(:,1),probe_coords_mean_wf(:,2),20,'k');

%% Ephys - Single cell

% Open Julie's fraction plot (AP_fractionCellsStimMove_mod)
% click stim fraction
stim_frac = get(gco,'xdata');
% click stim+move frac
stim_move_frac = get(gco,'xdata');
% click move frac
move_frac = get(gco,'xdata');

event_frac = [stim_frac;move_frac;stim_move_frac];
event_frac(4,:) = 1-sum(event_frac,1);

figure;
h = tiledlayout(4,1);
for curr_area = 1:4
    nexttile;
    pie(event_frac(:,curr_area));
end
event_col = [0.8,0,0;0.7,0,0.7;0.8,0.5,0.5;0,0,0];
colormap(event_col);
set(findobj('type','text'),'FontSize',14);



%% Plot widefield ROIs (hemisphere)

wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

plot_rois = [1,6,7];

roi_col = copper(n_rois);
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







