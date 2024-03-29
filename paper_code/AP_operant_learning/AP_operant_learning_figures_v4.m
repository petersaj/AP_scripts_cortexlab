% Generate figures for operant learning paper
% (data is prepared in AP_operant_learning_preprocessing)
%
% Code blocks with symbols: 
% load data first: run code block at top with same symbol
%
% v4: for revision

%% ------- LOAD DATA ---------------------------

%% >> Load task data

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

% Get movement hemisphere ratio (from rewardable no-stim movements)
% (average movement activity across days for each mouse and get hemiratio)
% (by ROI)
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

% (by pixel)
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


%% ++ Load passive data

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_passive_teto';
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


%% ------- GENERATE FIGS -----------------------

%% [FIG 1B]: example performance

animal = 'AP106';
protocol = 'AP_stimWheelRight';
experiments = AP_find_experiments(animal,protocol);
plot_days = [1,3,7];

figure;
h = tiledlayout(length(plot_days),1);
for curr_day = plot_days
 
    day = experiments(curr_day).day;
    experiment = experiments(curr_day).experiment;

    load_parts.imaging = false;
    AP_load_experiment;
    
    t = Timeline.rawDAQTimestamps;
    
    % Convert wheel velocity from clicks/s to mm/s
    % (mm in clicks from +hw.DaqRotaryEncoder, Lilrig encoder = 100)
    wheel_click2mm = 0.4869;
    wheel_velocity_mm = wheel_velocity*wheel_click2mm;

    % (minimum ITI: new trial + trial quiescence)
    min_iti_t = signals_events.newTrialTimes + ...
        signals_events.trialQuiescenceValues;
    
    plot_t = [100,150];
    plot_t_idx = t > plot_t(1) & t < plot_t(2);
    plot_stim_idx = find(stimOn_times > plot_t(1) & stimOn_times < plot_t(2))';
    plot_min_iti_t_idx = find(min_iti_t > plot_t(1) & min_iti_t < plot_t(2));
    plot_reward_idx = find(reward_t_timeline > plot_t(1) & reward_t_timeline < plot_t(2));
    
    nexttile; hold on;
    
    % Plot stim and rewards
    yyaxis left; hold on;

    area(t(plot_t_idx),stimOn_epochs(plot_t_idx), ...
        'EdgeColor','none','FaceColor',[1,1,0.8])
    for i = plot_reward_idx
        line(repmat(reward_t_timeline(i),2,1), ...
            [0,1],'color','b','linewidth',2);
    end

    % Plot wheel velocity
    yyaxis right;
    plot(t(plot_t_idx),wheel_velocity_mm(plot_t_idx),'k');
    
    axis off;
    title(sprintf('Day %d',curr_day));
    
end

h_ax = flipud(allchild(h));
linkaxes(h_ax,'y');

% Draw scalebars
t_scale = 5;
vel_scale = 100;
AP_scalebar(t_scale,vel_scale);


%% [FIG 1C-F, FIG 3C, FIG S1]: behavior

% Load behavior
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
bhv_fn = [data_path filesep 'bhv_teto'];
load(bhv_fn);
animals = {bhv.animal};

% Grab learned day
learned_day = cellfun(@(x) find(x,1),{bhv.learned_days})';

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
caxis([0,0.12]);

% Get histogram of reaction times in day groups (across animals)
day_grps = [1,3,5,7,Inf];
day_grp = discretize(day_idx,day_grps);

figure;
h = tiledlayout(length(day_grps)-1,1);
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
    xlim([0,0.5]);
end
linkaxes(allchild(h),'xy');

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
figure; h = tiledlayout('flow','TileSpacing','tight','padding','tight');
for curr_animal = 1:length(animals)
    nexttile; hold on; set(gca,'YScale','log');
    plot(rxn_measured_med(:,curr_animal),'k','linewidth',2);
    AP_errorfill([],[],prctile(rxn_alt_med(:,curr_animal,:),[5,95],3),[0.5,0.5,0.5]);
    xline(learned_day(curr_animal),'color','r');
    set(gca,'children',circshift(get(gca,'children'),-1))
    title(sprintf('Mouse %d',curr_animal));
    ylim([0,2]);
    set(gca,'YTick',[0,0.2,2,20]);
    if curr_animal == 1
        ylabel('Reaction time');
        xlabel('Training session');
    end
end
linkaxes(allchild(h),'xy');

% (plot median across mice)
figure; 
subplot(1,2,1,'YScale','log');hold on
rxn_alt_med_ci = squeeze(prctile(nanmedian(rxn_alt_med,2),[5,95],3));
AP_errorfill([],[],rxn_alt_med_ci(plot_days,:),[0.5,0.5,0.5],[],false);
errorbar(nanmedian(rxn_measured_med(plot_days,:),2), ...
    mad(rxn_measured_med(plot_days,:),1,2),'k','CapSize',0,'linewidth',2)
xlabel('Training day');
ylabel('Median reaction time (s)');
axis tight;
ylim([0,3]);
set(gca,'YTick',[0,0.25,0.5,1,2]);
xlim(xlim + [-0.5,0.5]);
ytickformat('%.2f')

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
line([0,0],ylim,'color','k','linestyle','--');
legend([p1(1),p2(1)],{'Null','Measured'});
set(gca,'YTick',[0,0.25,0.5,1,2]);
xlim(xlim + [-0.5,0.5]);
ylim([0,5]);
ytickformat('%.2f')

% Plot histogram of learned days
figure;histogram(learned_day,[1;plot_days]-0.5,'EdgeColor','none','FaceColor','k')
xlabel('Learned day');
ylabel('Number of mice');
xlim([0.5,max(plot_days)+0.5]);

% Plot frequency of non-response/response movements
nonresponse_move_rate = AP_padcatcell(cellfun(@(x,use_days) x(use_days), ...
    {bhv.nonresponse_move_rate},use_days,'uni',false));

trial_rate = AP_padcatcell(cellfun(@(n_trials,duration,use_days) ...
    n_trials(use_days)./(duration(use_days)*60), ...
    {bhv.n_trials},{bhv.session_duration},use_days,'uni',false));

figure; hold on
errorbar(nanmean(nonresponse_move_rate(plot_days,:),2), ...
    AP_sem(nonresponse_move_rate(plot_days,:),2),'k','linewidth',2);
errorbar(nanmean(trial_rate(plot_days,:),2), ...
    AP_sem(trial_rate(plot_days,:),2),'b','linewidth',2);
legend({'Non-response move rate','Trial rate'})
xlabel('Training day');
ylabel('Rate (number/sec)')



%% [FIG 1F-G, FIG S2A-B]: muscimol behavior and retinotopy

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

p = anova1(muscimol_v1_wheel_mm,[],'off');
fprintf('Muscimol wheel travel, 1-way anova = %.2d\n',p);

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
set(gca,'YTick',[0,0.25,0.5,1,2]);
xlim(xlim + [-0.5,0.5]);
ylim([0,3]);
ytickformat('%.2f')

p = anova1(rxn_measured_med',[],'off');
fprintf('Muscimol reaction time, 1-way anova = %.2d\n',p);

subplot(1,2,2); hold on
rxn_alt_med_altdiff_ci = permute(prctile(nanmean(rxn_alt_med_altdiff,2),[5,95],3),[1,3,2]);
AP_errorfill([],nanmean(rxn_alt_med_altdiff_ci,2),rxn_alt_med_altdiff_ci,[0.5,0.5,0.5],[],false);
errorbar(nanmean(rxn_measured_med_altdiff,2),AP_sem(rxn_measured_med_altdiff,2),'k','linewidth',2,'capsize',0);
xlim([0,n_conditions] + 0.5);
set(gca,'XTick',1:n_conditions,'XTickLabel',condition_labels);
ylabel('Median reaction time (meas-null,s)');


%% >> [FIG 2A]: task pixels stim/move aligned by novice/learned

% Set trials to use
learned_day_grp_edges = [-Inf,0,Inf];
trial_learned_day_grp = discretize(trial_learned_day,learned_day_grp_edges);
n_learn_grps = length(learned_day_grp_edges)-1;

% Set timepoints post stim/move
plot_t = {0.1,0};
align_labels = {'Stim','Move'};

% Get interpolated time points aligned to stim and move
align_v = cell(2,1);
for curr_align = 1:2

    switch curr_align
        case 1
            curr_v_align = fluor_allcat_deconv;
        case 2
            curr_v_align = fluor_move_allcat_deconv;
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

% Get average pixels across animals
align_px_avg = cellfun(@(x) ...
    nanmean(AP_svdFrameReconstruct(U_master(:,:,1:n_vs),x),5), ...
    align_v,'uni',false);

% Plot pixels
c = [0,0.007];
figure;
h = tiledlayout(n_learn_grps,length(cell2mat(plot_t)), ...
    'TileSpacing','compact','padding','compact');
for curr_learn_grp = 1:n_learn_grps
    for curr_align = 1:2
        for curr_plot_t = 1:length(plot_t{curr_align})
            
            curr_px = align_px_avg{curr_align}(:,:,curr_plot_t,curr_learn_grp);

            nexttile
            imagesc(curr_px);
            AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
            axis image off;
            colormap(gca,AP_colormap('WG',[],1.5));
            caxis(c);
            title([align_labels{curr_align} ': ' num2str(plot_t{curr_align}(curr_plot_t)) ' sec']);

        end
    end
end
linkaxes(allchild(h),'xy');
colorbar;


%% >> [FIG 2B]: task pixels hemidiff stim window by novice/learned

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

% Plot movies (for presentations)
AP_imscroll(fluor_learning_animal_px_avg,t)
caxis(max(abs(caxis))*[-1,1]);
colormap(AP_colormap('KWG',[],1.5));
axis image off;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

AP_imscroll(fluor_learning_animal_px_hemidiff_avg(:,1:size(U_master,2)/2,:,:),t)
caxis(max(abs(caxis))*[-1,1]);
colormap(AP_colormap('KWP',[],1.5));
axis image off;
AP_reference_outline('ccf_aligned_lefthemi',[0.5,0.5,0.5]);

% Get time window fluorescence
use_t = t >= 0 & t <= 0.2;
figure; tiledlayout(2,2);
c = [0,0.007];
c_hemi = [0,0.002];
for curr_learn = 1:n_learn_grps
    nexttile;
    imagesc(max(fluor_learning_animal_px_avg(:,:,use_t,curr_learn),[],3));
    colormap(gca,AP_colormap('WG',[],1.5)); caxis(c);
    axis image off;
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    colorbar

    nexttile;
    imagesc(max(fluor_learning_animal_px_hemidiff_avg(:,1:size(U_master,2)/2,use_t,curr_learn),[],3));
    colormap(gca,AP_colormap('WP',[],1.5)); caxis(c_hemi);
    axis image off;
    AP_reference_outline('ccf_aligned_lefthemi',[0.5,0.5,0.5]);
    colorbar
end


%% >> [FIG 2C]: task average L/R ROI timecourse (stim/move align, novice/learned)

plot_rois = [1,7,6];
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);

% Set activity
curr_act = fluor_roi_deconv;
curr_act_move = fluor_move_roi_deconv;

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

% Draw scalebars
x_scale = 0.1;
y_scale = 0.002;
AP_scalebar(x_scale,y_scale);

%% >> [FIG 2C ALT]: ROI stim-aligned sorted trials

plot_rois = [1,7,6];
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);
n_trial_smooth = 200;

% Set activity
roi_hemiflip = circshift(1:n_rois,n_rois/2);
curr_act = fluor_roi_deconv;

% Plot trial activity
figure;
h = tiledlayout(2,length(plot_rois));
for curr_stage = 1:max(trial_learned_stage)
    for curr_roi = plot_rois
        use_trials = find(trial_learned_stage == curr_stage);
        [~,sort_idx] = sort(move_t(use_trials));
        curr_data_sort = curr_act(use_trials(sort_idx),:,curr_roi);
        curr_data_sort_smooth = convn(curr_data_sort, ...
            ones(n_trial_smooth,1)./n_trial_smooth,'same');

        nexttile;
        imagesc(t,[],curr_data_sort_smooth);hold on;
        caxis(0.007.*[-1,1])
        colormap(AP_colormap('KWG',[],1));
        xline(0,'color',[0.8,0,0],'linewidth',2);
        plot(move_t(use_trials(sort_idx)),1:length(use_trials),'color','k','linewidth',2);
        xlabel('Time from stim (s)');
        ylabel('Trial (rxn-time sorted)');
        title(sprintf('%s, stage %d',wf_roi(curr_roi).area,curr_stage));
    end

end

linkaxes(allchild(h),'xy');
xlim([-0.2,1]);

ax = reshape(flipud(allchild(h)),[],2)';
x_scale = 0.2;
y_scale = 2000;
axes(ax(1,end)); AP_scalebar(x_scale,y_scale);
axes(ax(2,end)); AP_scalebar(x_scale,y_scale);
colorbar;


%% >> [FIG 2D]: task average hemidiff ROI timecourse (stim/move align, novice/learned)

plot_rois = [1,7,6];
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);

% Set activity
roi_hemiflip = circshift(1:n_rois,n_rois/2);
curr_act = fluor_roi_deconv - ...
    trial_move_hemiratio.*fluor_roi_deconv(:,:,roi_hemiflip);
curr_act_move = fluor_move_roi_deconv - ...
    trial_move_hemiratio.*fluor_move_roi_deconv(:,:,roi_hemiflip);

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
stage_col = min(1,AP_colormap('WP',1) + [0.5;0]);
for curr_roi = plot_rois

    curr_data_stim = permute(roi_stim_learn_avg(:,:,curr_roi,:),[2,1,4,3]);
    curr_data_move = permute(roi_move_learn_avg(:,:,curr_roi,:),[2,1,4,3]);

    nexttile; hold on;
    AP_errorfill(t,nanmean(curr_data_stim,3),AP_sem(curr_data_stim,3),stage_col);
    xlabel('Stim');
    ylabel(wf_roi(curr_roi).area(1:end-2));
    xline(0);yline(0);set(gca,'children',circshift(get(gca,'children'),-2))
    xlim(x_limit);

    p_t = t >= 0 & t <= 0.2;
    p = anova1(permute(max(curr_data_stim(p_t,:,:),[],1),[3,2,1]));    
    fprintf('Stim-aligned hemidiff t-max %s, 1-way anova p = %.2g\n',wf_roi(curr_roi).area(1:end-2),p);

    % Draw scalebars
    x_scale = 0.1;
    y_scale = 2e-4;
    AP_scalebar(x_scale,y_scale);

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


%% ++ [FIG 2E-G] passive pixels and L/R ROI in stim window

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
c = [0,0.003];
for curr_stage = 1:n_stages
    for curr_stim = 1:length(stim_unique)

        curr_px = nanmean(stim_px_avg_stage_tmax(:,:,curr_stage,curr_stim,:),5);

        nexttile;
        imagesc(curr_px);
        axis image off;
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        colormap(gca,AP_colormap('WG',[],1.5));
        caxis(c)
 
    end
end
linkaxes(allchild(h),'xy');
colorbar;

% Plot ROIs by stage (L/R overlay)
figure;
plot_rois = [1,6];

hemi_stage_col = min(1,reshape([0.7,0,0;0,0,0.7]' + cat(3,0.5,0),3,[]))';

h = tiledlayout(length(plot_rois),3,'TileSpacing','compact','padding','compact');
for curr_roi = plot_rois

    curr_rois = curr_roi + [0,size(wf_roi,1)];

    for curr_stim = stim_unique'

        nexttile;
        AP_errorfill(t, ...
            reshape(permute(nanmean(stim_roi_avg_stage(curr_rois,:,:,stim_unique == curr_stim,:),5),[2,1,3]),length(t),[]), ...
            reshape(permute(AP_sem(stim_roi_avg_stage(curr_rois,:,:,stim_unique == curr_stim,:),5),[2,1,3]),length(t),[]), ...
            hemi_stage_col(:,:));
        xlabel('Time from stim (s)');
        ylabel(sprintf('%s \\DeltaF/F_0',wf_roi(curr_roi).area));
        axis tight;xlim([-0.2,1])
    end

    % Stats (by learning and stimuli)
    [t_idx,stage_idx,stim_idx,animal_idx] = ndgrid(1:length(t),1:2,1:length(stim_unique),1:length(animals));
    [p_l,~,stats] = anovan(reshape(stim_roi_avg_stage(curr_rois(1),:,:,:,:),[],1), ...
        [t_idx(:),stage_idx(:),stim_idx(:)],'continuous',1,'model','interaction','display','on');
    [p_r,~,stats] = anovan(reshape(stim_roi_avg_stage(curr_rois(2),:,:,:,:),[],1), ...
        [t_idx(:),stage_idx(:),stim_idx(:)],'continuous',1,'model','interaction','display','off');
    fprintf('%s t/stage/stim 3-way anova p(stage) = %.2g\n',wf_roi(curr_rois(1)).area,p_l(2));
    fprintf('%s t/stage/stim 3-way anova p(stage) = %.2g\n',wf_roi(curr_rois(2)).area,p_r(2));

    % Stats (L vs R ROI on right-hand stimuli)
    use_t = t > 0 & t <= 0.2;
    curr_Lstim_tmax_act = squeeze(max(stim_roi_avg_stage(curr_rois,use_t,:,1,:),[],2));
    [roi_idx,stage_idx] = ndgrid(1:2,1:2,1:length(animals));
    [p,~,stats] = anovan(reshape(curr_Lstim_tmax_act,[],1), ...
        [roi_idx(:),stage_idx(:)],'model','interaction','display','off');
    fprintf('%s hemi/stage 2-way anova p(hemi) = %.2g, p(stage) = %.2g\n', ...
        wf_roi(curr_rois(1)).area(1:end-2),p(1),p(2));

end
% Link ROI axes
ax = reshape(flipud(allchild(h)),[],length(plot_rois))';
for curr_roi = 1:length(plot_rois)
    linkaxes(ax(curr_roi,:),'xy');

    % Draw scalebars
    axes(ax(curr_roi,end));
    x_scale = 0.2;
    y_scale = 4e-4;
    AP_scalebar(x_scale,y_scale);
end
% (shade stim area and put in back)
arrayfun(@(x) patch(x,[0,0.5,0.5,0], ...
    reshape(repmat(ylim(x),2,1),[],1),[1,1,0.8], ...
    'linestyle','none'),allchild(h));
arrayfun(@(x) set(x,'children',circshift(get(x,'children'),-1)),allchild(h));


%% >> [FIG 3A1] task hemidiff ROI timecourse by learned day

roi_hemiflip = circshift(1:n_rois,n_rois/2);
curr_act = fluor_roi_deconv - ...
    trial_move_hemiratio.*fluor_roi_deconv(:,:,roi_hemiflip);

% Get indicies for averaging 
[learned_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);
[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);
[trained_day_idx,~,~] = ...
    ndgrid(trial_day,1:length(t),1:n_rois);

accum_learnedday_idx = cat(4,learned_day_idx,t_idx,roi_idx,animal_idx);
accum_trainedday_idx = cat(4,trained_day_idx,t_idx,roi_idx,animal_idx);

use_trials = true(size(move_t));

roi_trainday_avg = accumarray( ...
    reshape(accum_trainedday_idx(use_trials,:,:,:),[],size(accum_trainedday_idx,4)), ...
    reshape(curr_act(use_trials,:,:),[],1), ...
    [max(learned_day_idx(:)),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

roi_learnday_avg = accumarray( ...
    reshape(accum_learnedday_idx(use_trials,:,:,:),[],size(accum_learnedday_idx,4)), ...
    reshape(curr_act(use_trials,:,:),[],1), ...
    [max(learned_day_idx(:)),length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

use_t = t > 0 & t <= 0.2;
stim_roi_learnday_avg_tmax = squeeze(max(roi_learnday_avg(:,use_t,:,:,:),[],2));

% Plot ROI timecourse first trained day / peri-learned days
plot_rois = [6];
plot_trained_days = 1;
plot_learned_days = -2:1;

figure;
plot_col = AP_colormap('WP',1);
h = tiledlayout(length(plot_rois), ...
    length(plot_trained_days) + length(plot_learned_days));

for curr_roi = plot_rois

    for curr_td = plot_trained_days
        curr_data = squeeze(roi_trainday_avg(curr_td,:,curr_roi,:));
        nexttile; hold on;
        AP_errorfill(t,nanmean(curr_data,2),AP_sem(curr_data,2),plot_col);
        xline(0); yline(0);
        title(sprintf('Trained day %d',curr_td));
    end   

    for curr_ld = plot_learned_days
        curr_data = squeeze(roi_learnday_avg(learned_day_unique == curr_ld,:,curr_roi,:));

        nexttile; hold on;
        AP_errorfill(t,nanmean(curr_data,2),AP_sem(curr_data,2),plot_col);
        xline(0); yline(0);
        title(sprintf('Learned day %d',curr_ld));
    end   
    ylabel(wf_roi(curr_roi).area);

     % Stats: difference between ld -2/-1 vs -1/0
    p_t = t >= 0 & t <= 0.2;
    p_prelearn = signrank( ...
        squeeze(max(roi_learnday_avg(learned_day_unique == -2,p_t,curr_roi,:),[],2)), ...
        squeeze(max(roi_learnday_avg(learned_day_unique == -1,p_t,curr_roi,:),[],2)));
    p_postlearn = signrank( ...
        squeeze(max(roi_learnday_avg(learned_day_unique == -1,p_t,curr_roi,:),[],2)), ...
        squeeze(max(roi_learnday_avg(learned_day_unique == -0,p_t,curr_roi,:),[],2)));
    fprintf('%s t-max signrank: p(-2 v -1) = %.2g, p(-1 v 0) = %.2g\n', ...
        wf_roi(curr_roi).area,p_prelearn,p_postlearn);

end

% (link ROI y-axes)
ax = reshape(flipud(allchild(h)),[],length(plot_rois))';
for curr_roi_idx = 1:length(plot_rois)
   linkaxes(ax(curr_roi_idx,:),'xy'); 
   xlim([-0.1,0.3]);
   ylim(prctile(reshape(cell2mat(ylim(ax(curr_roi_idx,:))),[],1),[0,100]));
end

% Draw scalebars
x_scale = 0.2;
y_scale = 2e-4;
AP_scalebar(x_scale,y_scale);


%% >> [FIG 3D(1)] task hemidiff ROI stim window daysplit by learned day

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

% Get x-values for plotting
learned_day_x_range = [min(learned_day_unique),max(learned_day_unique)];
learned_day_x = [learned_day_x_range(1):learned_day_x_range(2)]';
learned_daysplit_x = learned_day_x + ...
    (0:(n_daysplit+x_day_spacing-1))/(n_daysplit+x_day_spacing);

% Get days to plot with minimum n
plot_learned_day_idx = sum(~isnan(stim_roi_act_tmax_daysplit(:,1,1,:)),4) >= min_n;

% Plot
figure('Name','mPFC task learned day');
plot_rois = [6];
tiledlayout(length(plot_rois),1);
for curr_roi = plot_rois

    nexttile;
   
    errorbar(reshape(learned_daysplit_x(plot_learned_day_idx,:)',[],1), ...
        reshape(padarray(nanmean( ...
        stim_roi_act_tmax_daysplit(plot_learned_day_idx,:,curr_roi,:),4), ...
        [0,x_day_spacing],NaN,'post')',[],1), ...
        reshape(padarray(AP_sem( ...
        stim_roi_act_tmax_daysplit(plot_learned_day_idx,:,curr_roi,:),4), ...
        [0,x_day_spacing],NaN,'post')',[],1),'k','linewidth',2,'CapSize',0);

    xlabel('Learned day');
    ylabel(sprintf('%s \\DeltaF/F_0',wf_roi(curr_roi).area));
    xline(0,'linestyle','--');
    axis tight;
    xlim(xlim+[-0.5,0.5]);

end



%% ++ [FIG 3B,D(2), FIG S3] passive - L/R ROI timecourse/stim window by learned day

% Set ROIs to plot
plot_rois = [6];

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

roi_trainday_avg = accumarray( ...
    reshape(accum_idx(use_trials_training,:,:,[1,3,4,5,6]),[],5), ...
    reshape(fluor_roi_deconv(use_trials_training,:,:),[],1), ...
    [max(trial_day),length(t),n_rois,length(stim_unique),length(animals)], ...
    @nanmean,NaN('single'));

roi_learnday_avg = accumarray( ...
    reshape(accum_idx(use_trials_training,:,:,[2,3,4,5,6]),[],5), ...
    reshape(fluor_roi_deconv(use_trials_training,:,:),[],1), ...
    [length(learned_day_unique),length(t),n_rois,length(stim_unique),length(animals)], ...
    @nanmean,NaN('single'));

% Plot ROI timecourse by learned day (L/R)
plot_trained_days = 1;
plot_learned_days = -2:1;

figure;
hemi_col = [0.7,0,0;0,0,0.7];
h = tiledlayout(length(plot_rois),length(plot_trained_days) + length(plot_learned_days));
for curr_roi = plot_rois

    curr_lr_roi = curr_roi + [0,size(wf_roi,1)];

    for curr_td = plot_trained_days

        curr_td_postnaive = curr_td + n_naive;
        curr_data = squeeze(roi_trainday_avg(curr_td_postnaive,:, ...
            curr_lr_roi,stim_unique == 1,:));

        nexttile; hold on;
        AP_errorfill(t,nanmean(curr_data,3),AP_sem(curr_data,3),hemi_col);
        xline(0); yline(0);
        title(sprintf('Trained day %d',curr_td));
    end
    for curr_ld = plot_learned_days
        curr_data = squeeze(roi_learnday_avg(learned_day_unique == curr_ld,:, ...
        curr_lr_roi,stim_unique == 1,:));

        nexttile; hold on;
        AP_errorfill(t,nanmean(curr_data,3),AP_sem(curr_data,3),hemi_col);
        xline(0); yline(0);
        title(sprintf('Learned day %d',curr_ld));
    end   
    ylabel(wf_roi(curr_roi).area(1:end-2));

    % Stats: difference between ld -2/-1 vs -1/0
    p_t = t >= 0 & t <= 0.2;
    p_prelearn = signrank( ...
        squeeze(max(roi_learnday_avg(learned_day_unique == -2,p_t,curr_roi,stim_unique == 1,:),[],2)), ...
        squeeze(max(roi_learnday_avg(learned_day_unique == -1,p_t,curr_roi,stim_unique == 1,:),[],2)));
    p_postlearn = signrank( ...
        squeeze(max(roi_learnday_avg(learned_day_unique == -1,p_t,curr_roi,stim_unique == 1,:),[],2)), ...
        squeeze(max(roi_learnday_avg(learned_day_unique == 0,p_t,curr_roi,stim_unique == 1,:),[],2)));
    fprintf('%s t-max signrank: p(-2 v -1) = %.2g, p(-1 v 0) = %.2g\n', ...
        wf_roi(curr_roi).area,p_prelearn,p_postlearn);

end

% (link ROI y-axes)
ax = reshape(flipud(allchild(h)),[],length(plot_rois))';
for curr_roi_idx = 1:length(plot_rois)
    linkaxes(ax(curr_roi_idx,:),'xy');
    xlim(ax(curr_roi_idx,:),[-0.1,0.3]);
    ylim(ax(curr_roi_idx,:), ...
        prctile(reshape(cell2mat(ylim(ax(curr_roi_idx,:))),[],1),[0,100]));
end

% Draw scalebars
x_scale = 0.2;
y_scale = 4e-4;
AP_scalebar(x_scale,y_scale);

% Get ROI activity within stim window
use_t = t > 0 & t <= 0.2;

stim_roi_naive_avg_tmax = squeeze(max(roi_naive_avg(:,use_t,:,:,:),[],2));
stim_roi_learnday_avg_tmax = squeeze(max(roi_learnday_avg(:,use_t,:,:,:),[],2));

% Get days to plot by minimum n
learned_days_n = sum(any(any(stim_roi_learnday_avg_tmax,2),3),4);
plot_learned_day_idx = learned_days_n >= min_n;

% Plot ROI stim window activity by learned day
figure('Name','mPFC passive learned day');
plot_stim = [1,-1];
plot_stim_linestyle = {'-',':'};
h = tiledlayout(length(plot_rois),1);
for curr_roi = plot_rois

    curr_lr_roi = curr_roi + [0,size(wf_roi,1)];

    nexttile; hold on;
    set(gca,'ColorOrder',hemi_col);
    for curr_stim = plot_stim
        curr_stim_idx = ismember(stim_unique,curr_stim);
        % (naive and training as split lines)
        curr_x = [learned_day_unique(1:3);NaN;learned_day_unique(plot_learned_day_idx)];
        curr_data = squeeze([...
            padarray(stim_roi_naive_avg_tmax(:,curr_lr_roi,curr_stim_idx,:),1,NaN,'post'); ...
            stim_roi_learnday_avg_tmax(plot_learned_day_idx,curr_lr_roi,curr_stim_idx,:)]);
        
        p = errorbar(repmat(curr_x,1,length(plot_stim)), ...
            nanmean(curr_data,3),AP_sem(curr_data,3),'linewidth',2,'capsize',0, ...
            'linestyle',plot_stim_linestyle{ismember(plot_stim,curr_stim)});
  
    end
    xlabel('Time from stim (s)');
    ylabel(wf_roi(curr_roi).area(1:end-2));
    axis tight; xlim(xlim + [-0.5,0.5]);
    xline(0);
    legend({wf_roi(curr_lr_roi).area})

end

% STATS - DONT KNOW WHAT TO DO HERE YET
%
%     %%%%% MAYBE JUST DO A SHUFFLE TEST?
%     [ld_idx,stim_idx,~] = ndgrid(1:size(curr_data,1),1:length(stim_unique),1:length(animals));
%     [p,~,stats] = anovan(curr_data(~isnan(curr_data)),[ld_idx(~isnan(curr_data)),stim_idx(~isnan(curr_data))],'model','interaction','display','on');
%     c = multcompare(stats,'Dimension',[1,2]);
%     fprintf('%s t-max, 2-way anova p(stim) = %.2g,p(stage) = %.2g\n',wf_roi(curr_rois(1)).area,p_l(1),p_l(2));
%
%     curr_lrstim_diff = nanmean(abs(diff(curr_data(:,ismember(stim_unique,[-1,1]),:),[],2)),3);
%     n_shuff = 1000;
%     curr_lrstim_diff_shuff = nan(length(curr_x),n_shuff);
%     for curr_shuff = 1:n_shuff
%         curr_lrstim_diff_shuff(:,curr_shuff) = ...
%             nanmean(abs(diff(AP_shake(curr_data(:,ismember(stim_unique,[-1,1]),:),2),[],2)),3);
%     end
%     curr_lrstim_diff_rank = cell2mat(arrayfun(@(x) ...
%         tiedrank([curr_lrstim_diff(x),curr_lrstim_diff_shuff(x,:)]), ...
%         (1:length(curr_x))','uni',false));
%     curr_lrstim_diff_p = curr_lrstim_diff_rank(:,1)./(n_shuff+1);


%% ((run FIG 3D1-2)) [FIG 3D combined]: task & passive daysplit

% First make plots for passive and task daysplit mPFC activity

% Grab figure handles
fig_task = findall(groot,'type','figure','name','mPFC task learned day');
fig_passive = findall(groot,'type','figure','name','mPFC passive learned day');
if isempty(fig_task) || isempty(fig_passive)
    error('Task and passive figures not open');
end

% Grab errorbar handles
% (flip passive to get order plotted)
task_errorbar = findall(fig_task,'type','ErrorBar');
passive_errorbar = flipud(findall(fig_passive,'type','ErrorBar'));

% Pull data from task errorbar
xt = get(task_errorbar,'XData');
yt = get(task_errorbar,'YData');
et = get(task_errorbar,'UData');

% Pull data from passive errorbar
% (plot order: stim-L hemi-L/R, stim-R hemi L/R)
xp = get(passive_errorbar,'XData');
ypl = get(passive_errorbar,'YData');
epl = get(passive_errorbar,'UData');

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
task_col = AP_colormap('WP',1);
errorbar(xtr(:),ytr(:),etr(:),'color',task_col,'linewidth',2,'CapSize',0);

ylim([0,max(ytr(:))*1.3])
ylabel('\DeltaF/F_0');
ax = gca; ax.YColor = task_col;

% Plot passive
yyaxis right;
hemi_col = [0.7,0,0;0,0,0.7];
% (plot L-hemi R-stim)
curr_passive_plot = 1;
errorbar(xp{curr_passive_plot}+xp_offset,ypl{curr_passive_plot},epl{curr_passive_plot},'.','MarkerSize',15, ...
    'color',hemi_col(1,:),'linestyle','none','linewidth',2,'CapSize',0);
% (plot R-hemi R-stim)
curr_passive_plot = 2;
errorbar(xp{curr_passive_plot}+xp_offset,ypl{curr_passive_plot},epl{curr_passive_plot},'.','MarkerSize',15, ...
    'color',hemi_col(2,:),'linestyle','none','linewidth',2,'CapSize',0);

% (TESTING DIFFERENT KINDS OF PLOTS)
% % (plot L-hemi L-stim)
% curr_passive_plot = 3;
% errorbar(xp{curr_passive_plot}+xp_offset,ypl{curr_passive_plot},epl{curr_passive_plot},'o','MarkerSize',5, ...
%     'color',hemi_col(1,:),'linestyle','none','linewidth',2,'CapSize',0);
% % (plot R-hemi L-stim)
% curr_passive_plot = 4;
% errorbar(xp{curr_passive_plot}+xp_offset,ypl{curr_passive_plot},epl{curr_passive_plot},'o','MarkerSize',5, ...
%     'color',hemi_col(2,:),'linestyle','none','linewidth',2,'CapSize',0);

% % (plot L-hemi L-stim)
% curr_passive_plot = 3;
% AP_errorfill(xp{curr_passive_plot}+xp_offset,ypl{curr_passive_plot},epl{curr_passive_plot}, ...
%     hemi_col(1,:),0.2,false);
% % (plot R-hemi L-stim)
% curr_passive_plot = 4;
% AP_errorfill(xp{curr_passive_plot}+xp_offset,ypl{curr_passive_plot},epl{curr_passive_plot}, ...
%     hemi_col(2,:),0.2,false);


% % (plot L-hemi L-stim)
% curr_passive_plot = 3;
% errorbar(xp{curr_passive_plot}+xp_offset,ypl{curr_passive_plot},epl{curr_passive_plot},'MarkerSize',5, ...
%     'color',hemi_col(1,:),'linestyle','-','linewidth',2,'CapSize',0);
% % (plot R-hemi L-stim)
% curr_passive_plot = 4;
% errorbar(xp{curr_passive_plot}+xp_offset,ypl{curr_passive_plot},epl{curr_passive_plot},'MarkerSize',5, ...
%     'color',hemi_col(2,:),'linestyle','-','linewidth',2,'CapSize',0);

% curr_passive_plot = 3;
% errorbar(xp{curr_passive_plot}+xp_offset,ypl{curr_passive_plot},epl{curr_passive_plot},'*','MarkerSize',8, ...
%     'color',hemi_col(1,:),'linestyle','none','linewidth',2,'CapSize',0);

xlim(prctile(horzcat(xp{:}),[0,100]) + [-0.5,1.5]);
ylim([0,max(horzcat(ypl{:}))*1.3])
xline(0);
xlabel('Learned day');
ylabel('\DeltaF/F_0');
ax = gca; ax.YColor = 'k';


%% (OLD NOW) [FIG 4A]: ephys - plot probe position

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
scalebar_length = 1000/10; % um/voxel size
line(ccf_axes(3),[0,0],[0,scalebar_length],'color','m','linewidth',3);


%% [FIG 4B]: Example single unit rasters

animal = 'AP100';
day = '2021-05-25';
plot_units = [292,252,260];

% stim: 282 290 292
% stim+move: 270 252 286 233
% move: 260 301 393

% Set up subplots
figure;
h = tiledlayout(length(plot_units)*2,2,'tileindexing','columnmajor');

% Set raster time bins
raster_window = [-0.2,0.7];
psth_bin_size = 0.001;
raster_t_bins = raster_window(1):psth_bin_size:raster_window(2);
raster_t = raster_t_bins(1:end-1) + diff(raster_t_bins)./2;

for experiment = 1:2

    % (load LFP to find cortex start)
    preload_vars = who;
    lfp_channel = 'all';
    ephys_align = 'cortex';
    ephys_quality_control = false;
    AP_load_experiment;

    %%% JF code for loading good single units from bombcell
    % unitType: 0 = noise, 1 = good, 2 = multiunit
    ephysDirPath = AP_cortexlab_filename(animal, day, experiment, 'ephys_dir');
    qMetric_fn = fullfile(ephysDirPath, 'qMetrics');
    load(fullfile(qMetric_fn, 'qMetric.mat'))
    load(fullfile(qMetric_fn, 'param.mat'))
    clearvars unitType;

    % DEFAULT CHANGE: eliminate amplitude cutoff
    % (for one recording it got rid of almost all cells, and
    % cells under amplitude cutoff still look good)
    param.minAmplitude = 0;

    % (classify good cells)
    unitType = nan(length(qMetric.percSpikesMissing), 1);
    unitType( ...
        qMetric.nPeaks > param.maxNPeaks | ...
        qMetric.nTroughs > param.maxNTroughs | ...
        qMetric.somatic ~= param.somatic | ...
        qMetric.spatialDecaySlope <=  param.minSpatialDecaySlope | ...
        qMetric.waveformDuration < param.minWvDuration |...
        qMetric.waveformDuration > param.maxWvDuration  | ...
        qMetric.waveformBaseline >= param.maxWvBaselineFraction) = 0;
    unitType( ...
        any(qMetric.percSpikesMissing <= param.maxPercSpikesMissing, 2)' & ...
        qMetric.nSpikes > param.minNumSpikes & ...
        any(qMetric.Fp <= param.maxRPVviolations, 2)' & ...
        qMetric.rawAmplitude > param.minAmplitude & isnan(unitType)') = 1;
    unitType(isnan(unitType)') = 2;

    % Templates already 1/re-indexed, grab good ones
    good_templates = unitType == 1;
    good_templates_idx = find(good_templates);

    % Throw out all non-good template data
    templates = templates(good_templates,:,:);
    template_depths = template_depths(good_templates);
    waveforms = waveforms(good_templates,:);
    templateDuration = templateDuration(good_templates);
    templateDuration_us = templateDuration_us(good_templates);

    % Throw out all non-good spike data
    good_spike_idx = ismember(spike_templates,good_templates_idx);
    spike_times = spike_times(good_spike_idx);
    spike_templates_0idx = spike_templates_0idx(good_spike_idx);
    template_amplitudes = template_amplitudes(good_spike_idx);
    spike_depths = spike_depths(good_spike_idx);
    spike_times_timeline = spike_times_timeline(good_spike_idx);

    % Rename the spike templates according to the remaining templates
    % (and make 1-indexed from 0-indexed)
    new_spike_idx = nan(max(spike_templates_0idx)+1,1);
    new_spike_idx(unique(spike_templates_0idx)+1) = 1:length(unique(spike_templates_0idx));
    spike_templates = new_spike_idx(spike_templates_0idx+1);

    % Get area of each unit
    probe_area_boundary_starts = cellfun(@(x) x(1),probe_area_boundaries);
    [~,area_sort] = sort(probe_area_boundary_starts);

    unit_area = probe_areas(area_sort(discretize(template_depths, ...
        [-Inf;probe_area_boundary_starts(area_sort(2:end));Inf])));

    % Get raster alignment times
    if contains(expDef,'stimWheel')
        % If task, get rewardable movements and raster all events
        wheel_moves_deg = arrayfun(@(x) wheel_position_deg( ...
            Timeline.rawDAQTimestamps >= wheel_starts(x) & ...
            Timeline.rawDAQTimestamps <= wheel_stops(x)) - ...
            wheel_position_deg(find(Timeline.rawDAQTimestamps >= wheel_starts(x),1)), ...
            1:length(wheel_starts),'uni',false);

        deg_reward = -90;
        deg_punish = 90;

        wheel_moves_deg_rewardlimit = find(cellfun(@(x) ...
            any(x <= deg_reward) && ...
            ~any(x >= deg_punish),wheel_moves_deg));

        wheel_move_nostim_rewardable_idx = ...
            intersect(wheel_move_nostim_idx, wheel_moves_deg_rewardlimit);

        % Raster events: rewardable delay movements
        use_align = wheel_starts(wheel_move_nostim_rewardable_idx);

        align_label = 'Move onset';

    elseif contains(expDef,'Passive')
        % If passive, get quiescent movement and raster right-hand stim trials
        wheel_window = [0,0.5];
        wheel_window_t = wheel_window(1):1/Timeline.hw.daqSampleRate:wheel_window(2);
        wheel_window_t_peri_event = bsxfun(@plus,stimOn_times,wheel_window_t);
        event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
            +wheel_move,wheel_window_t_peri_event,'previous');
        quiescent_trials = ~any(abs(event_aligned_wheel) > 0,2);

        % Raster event: quiescent right-hand stim onset
        plot_trials = quiescent_trials & stimIDs == 3;
        use_align = stimOn_times(plot_trials);

        align_label = 'Stim onset';
    end

    t_peri_event = use_align + raster_t_bins;

    for curr_unit = plot_units

        % Bin spikes (use only spikes within time range, big speed-up)
        curr_spikes_idx = ismember(spike_templates,curr_unit);
        curr_raster_spike_times = spike_times_timeline(curr_spikes_idx);
        curr_raster_spike_times(curr_raster_spike_times < min(t_peri_event(:)) | ...
            curr_raster_spike_times > max(t_peri_event(:))) = [];

        curr_raster = cell2mat(arrayfun(@(x) ...
            histcounts(curr_raster_spike_times,t_peri_event(x,:)), ...
            [1:size(t_peri_event,1)]','uni',false));

        % Get smoothed PSTH
        smooth_size = 51;
        gw = gausswin(smooth_size,3)';
        smWin = gw./sum(gw);
        bin_t = mean(diff(raster_t));

        curr_psth = nanmean(curr_raster,1);
        curr_smoothed_psth = conv2(padarray(curr_psth, ...
            [0,floor(length(smWin)/2)],'replicate','both'), ...
            smWin,'valid')./bin_t;

        % Plot PSTH and raster
        nexttile;
        plot(raster_t,curr_smoothed_psth,'k','linewidth',1);
        axis tight off
        xline(0);
        title(sprintf('%s %s %s',animal,day,expDef), ...
            sprintf('unit %d %s',curr_unit,unit_area{curr_unit}));

        nexttile;
        [raster_y,raster_x] = find(curr_raster);
        plot(raster_t(raster_x),raster_y,'.k');
        axis tight off
        set(gca,'ydir','reverse');
        xlabel(align_label);

        drawnow;
    end

    clearvars('-except',preload_vars{:})

end

ax = reshape(flipud(allchild(h)),2,[]);
% Link all psth/raster y-axes
linkaxes(ax(1,:),'y')
linkaxes(ax(2,:),'y')

% Link all x-axes
linkaxes(ax,'x')

% Scalebars
t_scale = 0.2;
rate_scale = 20;

axes(ax(end-1));
AP_scalebar(t_scale,rate_scale)



%% (OLD NOW) [FIG 4C]: ephys - passive stim response

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

trial_recording = cell2mat(cellfun(@(tr,day) ...
    day*ones(size(tr,1),1), ...
    cat(1,wheel_all{:}),num2cell(1:length(cat(1,wheel_all{:})))', ...
    'uni',false));

trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));

% Get trials with movement during stim to exclude
quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > 0,2);

% Get average response timecourses
stim_unique = unique(trial_stim_allcat);
[~,trial_stim_id] = ismember(trial_stim_allcat,stim_unique);

use_trials = quiescent_trials;

% (loop through cell types)
for curr_celltype = 1:size(mua_area_allcat,4)

    [recording_idx,t_idx,area_idx] = ...
        ndgrid(trial_recording(use_trials),1:length(t),1:length(mua_areas));
    [stim_idx,~,~,] = ...
        ndgrid(trial_stim_id(use_trials),1:length(t),1:length(mua_areas));

    mua_recording_avg = accumarray([recording_idx(:),t_idx(:),stim_idx(:),area_idx(:)], ...
        reshape(mua_area_allcat(use_trials,:,:,curr_celltype),[],1), ...
        [length(trials_recording),length(t),length(stim_unique),length(mua_areas)], ...
        @nanmean,NaN);

    % Get DV position of each area in CCF for sorting
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

    % (get areas to plot: anything present in all recordings)
    [~,area_idx] = cellfun(@(x) ismember(x,mua_areas),mua_areas_cat,'uni',false);
    area_recording_n = accumarray(cell2mat(area_idx),1);
    plot_areas = find(area_recording_n == length(trials_recording));
    [~,plot_area_sort_idx] = sort(mua_areas_dvmin(plot_areas));

    % Plot stim overlaid
    figure;
    stim_color = [0,0,0.8;0.5,0.5,0.5;0.8,0,0];
    h = tiledlayout(length(plot_areas),1);
    for curr_area = plot_areas(plot_area_sort_idx)'
        nexttile;
        AP_errorfill(t', ...
            squeeze(nanmean(mua_recording_avg(:,:,:,curr_area),1)), ...
            squeeze(AP_sem(mua_recording_avg(:,:,:,curr_area),1)),stim_color);
        yline(0);
        ylabel(mua_areas(curr_area));
    end
    linkaxes(allchild(h),'xy');
    % (shade stim area and put in back)
    arrayfun(@(x) patch(x,[0,0.5,0.5,0], ...
        reshape(repmat(ylim(x),2,1),[],1),[1,1,0.8], ...
        'linestyle','none'),allchild(h));
    arrayfun(@(x) set(x,'children',circshift(get(x,'children'),-1)),allchild(h));
    xlim([-0.2,0.7]);

    x_scale = 0.2;
    y_scale = 0.5;
    AP_scalebar(x_scale,y_scale);

    title(h,sprintf('Cell type %d',curr_celltype));
end


%% [FIG 4D] ephys - move (no stim) response

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

trial_recording = cell2mat(cellfun(@(tr,day) ...
    day*ones(size(tr,1),1), ...
    cat(1,wheel_all{:}),num2cell(1:length(cat(1,wheel_all{:})))', ...
    'uni',false));

trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));

% Plot average activity aligned to stim and move

% (loop through cell types)
celltype_label = {'Wide','Narrow'};
for curr_celltype = 1:length(celltype_label)

    % Get average response timecourses
    mua_move_nostim_recording_avg = ...
        mua_area_move_nostim_rewardable_allcat(:,:,:,curr_celltype);

    % Get DV position of each area in CCF for sorting
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

    % Plot move no-stim aligned
    figure;
    h = tiledlayout(length(plot_areas),1);
    for curr_area = plot_areas(plot_area_sort_idx)'
        nexttile;
        AP_errorfill(t', ...
            squeeze(nanmean(mua_move_nostim_recording_avg(:,:,curr_area),1)), ...
            squeeze(AP_sem(mua_move_nostim_recording_avg(:,:,curr_area),1)),'k');
        xlabel('Move');
        xline(0);yline(0);
        ylabel(mua_areas(curr_area));
    end
    linkaxes(allchild(h),'xy');
    xlim([-0.2,0.7]);

    x_scale = 0.2;
    y_scale = 0.5;
    AP_scalebar(x_scale,y_scale);

    title(h,celltype_label{curr_celltype});

end


%% [FIG 4E-F] ephys single cell classification and stim/move response

data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'single_unit_data_all';
load(fullfile(data_path,data_fn));

% Unpack (event rates are stored as event x unit x pre/post)
unit_area_cat = horzcat(single_unit_data_all.unit_area)';
n_units_exp = cellfun(@length,unit_area_cat);
waveform_duration_exp = horzcat(single_unit_data_all.waveform_duration_JF);

passive_stim_fr_cat = horzcat(single_unit_data_all.passive_stim_fr)';
task_move_fr_cat = horzcat(single_unit_data_all.task_move_fr)';

passive_stim_psth_allcat = cell2mat(horzcat(single_unit_data_all.passive_stim_psth)');
task_move_psth_allcat = cell2mat(horzcat(single_unit_data_all.task_move_psth)');
% (PSTH time - just re-make here)
raster_window = [-0.5,1];
raster_sample_rate = 50;
raster_sample_time = 1/raster_sample_rate;
t = raster_window(1):raster_sample_time:raster_window(2);

% Classify narrow vs. wide waveforms
% (cutoff from Bartho JNeurophys 2004)
% (wide = 1, narrow = 2);
waveform_duration_cutoff = 400;
celltype_exp = cellfun(@(x) (x < waveform_duration_cutoff)+1, ...
    waveform_duration_exp,'uni',false)';
n_celltypes = length(unique(cell2mat(celltype_exp)));
celltype_label = {'Wide','Narrow'};


% (sanity check: plot waveforms by cell type)
waveform_allcat = cell2mat(horzcat(single_unit_data_all.waveform)');
figure; hold on
celltype_col = lines(max(vertcat(celltype_exp{:})));
for curr_celltype = 1:max(vertcat(celltype_exp{:}))
    plot(waveform_allcat(cell2mat(celltype_exp) == curr_celltype,:)', ...
        'color',celltype_col(curr_celltype,:));
end
wp = gobjects(max(vertcat(celltype_exp{:})),1);
for curr_celltype = 1:max(vertcat(celltype_exp{:}))
    wp(curr_celltype) = ...
        plot(nanmean(waveform_allcat(cell2mat(celltype_exp) == curr_celltype,:),1), ...
        'color',max(0,celltype_col(curr_celltype,:)-0.4),'linewidth',3);
end
legend(wp,celltype_label);

% Get firing rate change for stim (passive) and movement (delay)
passive_stim_fr_diff = cellfun(@(x) nanmean(diff(x,[],3),1), ...
    passive_stim_fr_cat,'uni',false);
task_move_fr_diff = cellfun(@(x) nanmean(diff(x,[],3),1), ...
    task_move_fr_cat,'uni',false);

% Get significance from pre/post event shuffle
stim_sig = cell(length(unit_area_cat),1);
move_sig = cell(length(unit_area_cat),1);
for curr_exp = 1:length(unit_area_cat)

    curr_n_units = length(unit_area_cat{curr_exp});

    n_shuff = 1000;
    stim_fr_diff_shuff = nan(n_shuff,curr_n_units);
    move_fr_diff_shuff = nan(n_shuff,curr_n_units);
    for curr_shuff = 1:n_shuff
        stim_fr_diff_shuff(curr_shuff,:) = ...
            nanmean(diff(AP_shake(passive_stim_fr_cat{curr_exp},3),[],3),1);
        move_fr_diff_shuff(curr_shuff,:) = ...
            nanmean(diff(AP_shake(task_move_fr_cat{curr_exp},3),[],3),1);
    end

    % p-value of absolute difference
    stim_fr_diff_p = 1 - cell2mat(arrayfun(@(x) ...
        tiedrank(abs([passive_stim_fr_diff{curr_exp}(x);stim_fr_diff_shuff(:,x)])), ...
        1:curr_n_units,'uni',false))./(n_shuff+1);
    move_fr_diff_p = 1 - cell2mat(arrayfun(@(x) ...
        tiedrank(abs([task_move_fr_diff{curr_exp}(x);move_fr_diff_shuff(:,x)])), ...
        1:curr_n_units,'uni',false))./(n_shuff+1);

    % Define significance
    sig_thresh = 0.01;
    stim_fr_diff_sig = stim_fr_diff_p(1,:) <= sig_thresh;
    move_fr_diff_sig = move_fr_diff_p(1,:) <= sig_thresh;

    % Store 
    stim_sig{curr_exp} = stim_fr_diff_sig';
    move_sig{curr_exp} = move_fr_diff_sig';

end

% Hard-code areas to plot
plot_areas = ...
    {'Secondary motor area', ...
    'Anterior cingulate area dorsal part', ...
    'Prelimbic area', ...
    'Infralimbic area'};

% Get and plot fraction of responsive cells by region/experiment
[~,unit_area_idx] = cellfun(@(x) ismember(x,plot_areas),unit_area_cat,'uni',false);

stim_frac = cell2mat(permute(cellfun(@(area,celltype,stim_sig,move_sig) ...
    accumarray([area(area ~= 0),celltype(area ~= 0)], ...
    stim_sig(area ~= 0) & ~move_sig(area ~= 0), ...
    [length(plot_areas),2],@mean), ...
    unit_area_idx,celltype_exp,stim_sig,move_sig,'uni',false),[3,2,1]));

move_frac = cell2mat(permute(cellfun(@(area,celltype,stim_sig,move_sig) ...
    accumarray([area(area ~= 0),celltype(area ~= 0)], ...
    ~stim_sig(area ~= 0) & move_sig(area ~= 0), ...
    [length(plot_areas),2],@mean), ...
    unit_area_idx,celltype_exp,stim_sig,move_sig,'uni',false),[3,2,1]));

stimmove_frac = cell2mat(permute(cellfun(@(area,celltype,stim_sig,move_sig) ...
    accumarray([area(area ~= 0),celltype(area ~= 0)], ...
    stim_sig(area ~= 0) & move_sig(area ~= 0), ...
    [length(plot_areas),2],@mean), ...
    unit_area_idx,celltype_exp,stim_sig,move_sig,'uni',false),[3,2,1]));

class_frac = permute(cat(3,nanmean(stim_frac,3), ...
    nanmean(stimmove_frac,3),nanmean(move_frac,3)),[1,3,2]);
class_col = [0.8,0,0;0.8,0.5,0.5;0.5,0.5,0.5;1,1,1];

figure; h = tiledlayout(length(plot_areas),n_celltypes);
for curr_area = 1:length(plot_areas)
    for curr_celltype = 1:n_celltypes
        nexttile;
        curr_data = [class_frac(curr_area,:,curr_celltype), ...
            1-sum(class_frac(curr_area,:,curr_celltype))];
        pie(curr_data);
        title(sprintf('%s %s',plot_areas{curr_area},celltype_label{curr_celltype}));
    end
end
colormap(class_col);
title(h,'Average fraction across experiments')

% (plot as bar plot also)
class_frac_sem = permute(cat(3,AP_sem(stim_frac,3), ...
    AP_sem(stimmove_frac,3),AP_sem(move_frac,3)),[1,3,2]);

figure; h = tiledlayout(1,n_celltypes);
for curr_celltype = 1:n_celltypes
    nexttile; hold on; set(gca,'ColorOrder',class_col);
    curr_data = [class_frac(:,:,curr_celltype), ...
        1-sum(class_frac(:,:,curr_celltype),2)];
    bar(curr_data,'stacked');

    curr_sem = class_frac_sem(:,:,curr_celltype);
    errorbar(cumsum(curr_data(:,1:end-1),2),curr_sem,'k','linewidth',2,'linestyle','none');
    
    title(celltype_label{curr_celltype});
    set(gca,'XTick',1:length(plot_areas),'XTickLabel',plot_areas);
    ylabel('Fraction of cells');
end


% Get and plot fraction of responsive cells (all combined)
unit_area_idx_allcat = cell2mat(unit_area_idx);
celltype_allcat = cell2mat(celltype_exp);
stim_sig_allcat = cell2mat(stim_sig);
move_sig_allcat = cell2mat(move_sig);

use_cells = unit_area_idx_allcat ~= 0;
stim_frac_allcat = ...
    accumarray([unit_area_idx_allcat(use_cells),celltype_allcat(use_cells)], ...
    stim_sig_allcat(use_cells) & ~move_sig_allcat(use_cells),[],@mean);
move_frac_allcat = ...
    accumarray([unit_area_idx_allcat(use_cells),celltype_allcat(use_cells)], ...
    ~stim_sig_allcat(use_cells) & move_sig_allcat(use_cells),[],@mean);
stimmove_frac_allcat = ...
    accumarray([unit_area_idx_allcat(use_cells),celltype_allcat(use_cells)], ...
    stim_sig_allcat(use_cells) & move_sig_allcat(use_cells),[],@mean);

class_frac_allcat = permute(cat(3,stim_frac_allcat, ...
    stimmove_frac_allcat,move_frac_allcat),[1,3,2]);

figure; h = tiledlayout(length(plot_areas),n_celltypes);
for curr_area = 1:length(plot_areas)
    for curr_celltype = 1:n_celltypes
        nexttile;
        pie([class_frac_allcat(curr_area,:,curr_celltype), ...
            1-sum(class_frac_allcat(curr_area,:,curr_celltype))]);
        title(sprintf('%s %s',plot_areas{curr_area},celltype_label{curr_celltype}));
    end
end
colormap(class_col);
title(h,'All cells pooled')

% Print mean and sem of stim and move cells in MOs and ACA
stim_frac_dmpfc = cellfun(@(area,sig) ...
    nanmean(sig(ismember(area,[1,2]))), ... 
    unit_area_idx,stim_sig);

move_frac_dmpfc = cellfun(@(area,sig) ...
    nanmean(sig(ismember(area,[1,2]))), ... 
    unit_area_idx,move_sig);

fprintf('Stim dmPFC: %0.2f +- %0.2f\n',nanmean(stim_frac_dmpfc),std(stim_frac_dmpfc));
fprintf('Move dmPFC: %0.2f +- %0.2f\n',nanmean(move_frac_dmpfc),std(move_frac_dmpfc));

% Get chance of stim & move cells (combined areas of interest)
stimmove_frac_total = cellfun(@(area,stim,move) ...
    mean(stim(area~=0) & move(area~=0)),unit_area_idx,stim_sig,move_sig);

n_shuff = 1000;
stimmove_frac_total_shuff = nan(length(unit_area_cat),n_shuff);
for curr_shuff = 1:n_shuff

    % (shuffle movement significance by area)
    curr_move_sig_shuff = cellfun(@(sig,area) ...
        AP_shake(sig,1,area),move_sig,unit_area_cat,'uni',false);

    % (get fraction of sitim & move in shuffle)
    stimmove_frac_total_shuff(:,curr_shuff) = ...
        cellfun(@(area,stim,move) ...
        mean(stim(area~=0) & move(area~=0)), ...
        unit_area_idx,stim_sig,curr_move_sig_shuff);

end

stimmove_frac_rank = tiedrank([nanmean(stimmove_frac_total,1), ...
    nanmean(stimmove_frac_total_shuff,1)]);
stimmove_frac_p = 1 - stimmove_frac_rank(1)./(n_shuff+1);
fprintf('Stim & move shuffle p = %.2f\n',stimmove_frac_p);



%%%%%%%%%%%%%% IN PROGRESS ADDITION: PSTHS AND NARROW/WIDE DIFF STATS


% DO HERE: sort by stim FR diff, and then category?
passive_stim_fr_diff_cat = horzcat(passive_stim_fr_diff{:})';
task_move_fr_diff_cat = horzcat(task_move_fr_diff{:})';

baseline_t = t < -0.2;
softnorm = 1;
passive_stim_psth_allcat_norm = ...
    (passive_stim_psth_allcat - nanmean(passive_stim_psth_allcat(:,baseline_t),2))./ ...
    (nanmean(passive_stim_psth_allcat(:,baseline_t),2) + softnorm);
task_move_psth_allcat_norm = ...
    (task_move_psth_allcat - nanmean(task_move_psth_allcat(:,baseline_t),2))./ ...
    (nanmean(task_move_psth_allcat(:,baseline_t),2) + softnorm);

figure; h = tiledlayout(length(plot_areas),2*n_celltypes);
sig_idx = sum([stim_sig_allcat & ~move_sig_allcat, ...
    stim_sig_allcat & move_sig_allcat, ...
    ~stim_sig_allcat & move_sig_allcat].*[1:3],2);
for curr_area = 1:length(plot_areas)
    for curr_celltype = 1:n_celltypes

        curr_cells = find(unit_area_idx_allcat == curr_area & ...
            celltype_allcat == curr_celltype & stim_sig_allcat);
        
        % (sort by passive fr diff?)
        [~,sort_idx] = sortrows([passive_stim_fr_diff_cat(curr_cells), ...
            task_move_fr_diff_cat(curr_cells)],[1,2],'descend');

        plot_cells = curr_cells(sort_idx);

        smooth_filt = ones(1,3)/3;

        nexttile;
        imagesc(t,[],convn(passive_stim_psth_allcat_norm(plot_cells,:),smooth_filt,'same'));
        colormap(AP_colormap('BWR'));
        caxis([-1,1].*5);

        nexttile;
        imagesc(t,[],convn(task_move_psth_allcat_norm(plot_cells,:),smooth_filt,'same'));
        colormap(AP_colormap('BWR'));
        caxis([-1,1].*5);

    end
end
% linkaxes(allchild(h),'xy');

% Difference in stim responses in narrow/wide cells
stim_sig_meas = nanmean(cell2mat(permute(cellfun(@(area,celltype,stim_sig,move_sig) ...
    accumarray([area(area ~= 0),celltype(area ~= 0)], ...
    stim_sig(area ~= 0), ...
    [length(plot_areas),n_celltypes],@mean), ...
    unit_area_idx,celltype_exp,stim_sig,move_sig,'uni',false),[3,2,1])),3);

n_shuff = 1000;
stim_sig_shuff = nan(length(plot_areas),n_celltypes,n_shuff);
for curr_shuff = 1:n_shuff
    stim_sig_shuff(:,:,curr_shuff) = nanmean( ...
        cell2mat(permute(cellfun(@(area,celltype,stim_sig,move_sig) ...
        accumarray([area(area ~= 0),AP_shake(celltype(area ~= 0))], ...
        stim_sig(area ~= 0), ...
        [length(plot_areas),n_celltypes],@mean), ...
        unit_area_idx,celltype_exp,stim_sig,move_sig,'uni',false),[3,2,1])),3);
end

figure; hold on;
AP_errorfill([],[],prctile(diff(stim_sig_shuff,[],2),[5,95],3));
plot(diff(stim_sig_meas,[],2),'k','linewidth',2);
set(gca,'XTick',1:length(plot_areas),'XTickLabel',plot_areas);
ylabel('Narrow - Wide fraction stim-responsive');

%%%%% TRYING OUT: as above, but just combine dmPFC
%%%%%%%%%%%%%%%%%% --->>> I THINK THIS IS GOOD,= RUN WITH THIS

baseline_t = t < -0.3;
softnorm = 10;
passive_stim_psth_allcat_norm = ...
    (passive_stim_psth_allcat - nanmean(passive_stim_psth_allcat(:,baseline_t),2))./ ...
    (nanmean(passive_stim_psth_allcat(:,baseline_t),2) + softnorm);
task_move_psth_allcat_norm = ...
    (task_move_psth_allcat - nanmean(task_move_psth_allcat(:,baseline_t),2))./ ...
    (nanmean(task_move_psth_allcat(:,baseline_t),2) + softnorm);


% Plot PSTH for dmPFC
figure; h = tiledlayout(2,n_celltypes);
for curr_celltype = 1:n_celltypes

    curr_cells = find(ismember(unit_area_idx_allcat,[1,2]) & ...
        celltype_allcat == curr_celltype & stim_sig_allcat);

    % (sort by passive fr diff?)
    [~,sort_idx] = sortrows([move_sig_allcat(curr_cells), ...
        passive_stim_fr_diff_cat(curr_cells)],[1,2],'descend');
%     [~,sort_idx] = sort(passive_stim_fr_diff_cat(curr_cells),'descend');

    plot_cells = curr_cells(sort_idx);

    smooth_filt = ones(1,3)/3;

    nexttile;
    imagesc(t,[],convn(passive_stim_psth_allcat_norm(plot_cells,:),smooth_filt,'same'));
    colormap(AP_colormap('BWR'));
    caxis([-1,1].*2);
    xline(0,'color','k','linewidth',2);
    yline(sum(move_sig_allcat(curr_cells)),'color','k','linewidth',2);
    title(sprintf('dmPFC %s Stim',celltype_label{curr_celltype}));
    
    nexttile;
    imagesc(t,[],convn(task_move_psth_allcat_norm(plot_cells,:),smooth_filt,'same'));
    colormap(AP_colormap('BWR'));
    caxis([-1,1].*2);
    xline(0,'color','k','linewidth',2);
    yline(sum(move_sig_allcat(curr_cells)),'color','k','linewidth',2);
    title(sprintf('dmPFC %s Move',celltype_label{curr_celltype}));

end

dmpfc_units = cellfun(@(area) ismember(area,[1,2]),unit_area_idx,'uni',false);
class_idx = cellfun(@(stim_sig,move_sig) sum( ...
    [(stim_sig & ~move_sig), ...
    (stim_sig & move_sig), ...
    (~stim_sig & move_sig), ...
    (~stim_sig & ~move_sig)].*[1:4],2), ...
    stim_sig,move_sig,'uni',false); % (1:4 = stim,stim+move,move,none)

dmpfc_class_frac = cell2mat(permute(cellfun(@(dmpfc,celltype,class_idx) ...
    accumarray([celltype(dmpfc),class_idx(dmpfc)],1,[],@sum)./ ...
    accumarray(celltype(dmpfc),1,[],@sum), ...
    dmpfc_units,celltype_exp,class_idx,'uni',false),[2,3,1]));

figure; hold on;
class_col = [0.8,0,0;0.8,0.5,0.5;0.5,0.5,0.5;1,1,1];
set(gca,'ColorOrder',class_col);
bar(nanmean(dmpfc_class_frac,3),'stacked');
errorbar(cumsum(nanmean(dmpfc_class_frac,3),2), ...
    AP_sem(dmpfc_class_frac,3),'k','linewidth',2,'linestyle','none');
ylabel('Fraction units')
title('dmPFC')
set(gca,'XTick',1:n_celltypes,'XTickLabel',celltype_label);


%%% (shuffle for dmpfc)

n_shuff = 1000;
dmpfc_class_frac_shuff = nan(n_celltypes,max(cell2mat(class_idx)),n_shuff);
for curr_shuff = 1:n_shuff
    dmpfc_class_frac_shuff(:,:,curr_shuff) = nanmean( ...
        cell2mat(permute(cellfun(@(dmpfc,celltype,class_idx) ...
        accumarray([AP_shake(celltype(dmpfc)),class_idx(dmpfc)],1,[],@sum)./ ...
        accumarray(celltype(dmpfc),1,[],@sum), ...
        dmpfc_units,celltype_exp,class_idx,'uni',false),[2,3,1])),3);
end

dmpfc_celltype_diff_rank = ...
    tiedrank(vertcat(diff(nanmean(dmpfc_class_frac,3),[],1), ...
    permute(diff(dmpfc_class_frac_shuff,[],1),[3,2,1])));
dmpfc_celltype_diff_p = dmpfc_celltype_diff_rank(1,:)./(n_shuff+1);
fprintf('\n Stim p = %.2f\n Stim & Move p = %.2f\n Move p = %.2f\n', ...
    dmpfc_celltype_diff_p(1),dmpfc_celltype_diff_p(2),dmpfc_celltype_diff_p(3));

% dmPFC waveforms
dmpfc_units_allcat = cell2mat(dmpfc_units);
waveform_allcat = cell2mat(horzcat(single_unit_data_all.waveform)');
[dmpfc_waveform_mean,dmpfc_waveform_std] = ...
    grpstats(waveform_allcat(dmpfc_units_allcat,:), ...
    celltype_allcat(dmpfc_units_allcat),{'mean','std'});
figure; hold on
AP_errorfill([],dmpfc_waveform_mean',dmpfc_waveform_std');

celltype_col = lines(max(vertcat(celltype_exp{:})));
for curr_celltype = 1:max(vertcat(celltype_exp{:}))
    plot(waveform_allcat(cell2mat(celltype_exp) == curr_celltype,:)', ...
        'color',celltype_col(curr_celltype,:));
end
wp = gobjects(max(vertcat(celltype_exp{:})),1);
for curr_celltype = 1:max(vertcat(celltype_exp{:}))
    wp(curr_celltype) = ...
        plot(nanmean(waveform_allcat(cell2mat(celltype_exp) == curr_celltype,:),1), ...
        'color',max(0,celltype_col(curr_celltype,:)-0.4),'linewidth',3);
end
legend(wp,celltype_label);


%% [FIG S2C-E]: V1 muscimol ROI activity pre/post

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

% Plot average ROIs
stim_roi_avg_stage = cat(5,stim_roi_avg{:});

figure;
plot_rois = [1,6];
stage_col = [0.5,0.7,0.7;0,0,0];
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
    xline([0,0.5]);

    % (plot right ROI with left stim)
    hr =nexttile;
    AP_errorfill(t, ...
        squeeze(nanmean(stim_roi_avg_stage(curr_r_roi,:,:,stim_unique == -1,:),5)), ...
        squeeze(AP_sem(stim_roi_avg_stage(curr_r_roi,:,:,stim_unique == -1,:),5)),stage_col);
    xlabel('Time from stim (s)');
    ylabel('\DeltaF/F_0');
    title(wf_roi(curr_r_roi).area);
    xline([0,0.5]);

    % (link axes with same ROI)
    linkaxes([hl,hr],'xy');
end
% (link ROI y-axes)
ax = reshape(flipud(allchild(h)),2,length(plot_rois))';
for curr_roi = 1:length(plot_rois)
   linkaxes(ax(:,curr_roi),'y'); 
end
linkaxes(ax,'x');
xlim([-0.2,0.7]);
% (draw scalebars)
x_scale = 0.2;
y_scale = 3e-3;
AP_scalebar(x_scale,y_scale);

% Plot ROI time-max pre/post muscimol
use_t = t >= 0 & t <= 0.2;
stim_roi_avg_stage_tmax = squeeze(max(stim_roi_avg_stage(:,use_t,:,:,:),[],2));

plot_rois = [1,6];
plot_stim = 3;
figure; hold on
curr_data = squeeze(stim_roi_avg_stage_tmax(plot_rois,:,plot_stim,:));
plot(squeeze(curr_data(1,:,:)),squeeze(curr_data(2,:,:)),'k');
p1 = plot(squeeze(curr_data(1,2,:)),squeeze(curr_data(2,2,:)),'.','MarkerSize',30,'color',stage_col(2,:));
p2 = plot(squeeze(curr_data(1,1,:)),squeeze(curr_data(2,1,:)),'.','MarkerSize',30,'color',stage_col(1,:));
xlabel(wf_roi(plot_rois(1)).area);
ylabel(wf_roi(plot_rois(2)).area);
legend([p1,p2],{'Washout','V1 Muscimol'},'location','nw');
axis square;

ax = gca;
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;

% Stats (few values - do left-sided test)
p = signrank(squeeze(curr_data(2,1,:)),squeeze(curr_data(2,2,:)),'tail','left');
fprintf('%s muscimol left-sided signed-rank: p = %.2g\n',wf_roi(plot_rois(2)).area,p);


%% [FIG S4A]: example pupil and facecam frame

animal = 'AP106';
day = '2021-06-25';
experiment = 2;

load_parts.cam = true;
verbose = true;
AP_load_experiment;

% Pick a frame from the movie
AP_mousemovie(eyecam_fn,eyecam_t,eyecam_dlc);

% Get aligned whisker mask for animal/day
facecam_processing_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\facecam_processing';
facecam_align_fn = fullfile(facecam_processing_path,'facecam_align.mat');
load(facecam_align_fn);

facecam_align_animalidx = find(strcmp(animal,{facecam_align.animal}));
facecam_align_dayidx = find(strcmp(day,[facecam_align(facecam_align_animalidx).day]));

whisker_mask = facecam_align(facecam_align_animalidx).whisker_mask{facecam_align_dayidx};

% Grab sample facecam frame (last frame)
vr = VideoReader(facecam_fn);
grab_frame = vr.NumFrames;
facecam_sample_frame = read(vr,grab_frame);

% Plot facecam frame and whisker ROI overlaid
figure;
image(imoverlay(mat2gray(facecam_sample_frame),whisker_mask,'r'));
axis image off




%% ++ [FIG S4 C-G]: passive whisker/pupil and ROI activity

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
h = tiledlayout(size(bhv_stim_avg,5),max(learning_stage));
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
        xlim([-0.2,0.7]);
        xline(0);xline(0.5);
        yline(0);
    end
end
% Link bhv y-axes
ax = reshape(flipud(allchild(h)),2,2)';
for curr_bhv = 1:2
    linkaxes(ax(curr_bhv,:),'y');
    switch curr_bhv
        case 1
            y_scale = 0.04;
        case 2
            y_scale = 0.1;
    end
    axes(ax(curr_bhv,2));
    x_scale = 0.2;
    AP_scalebar(x_scale,y_scale);
end

% Plot whisker movement over time
use_t = t > 0 & t <= 0.2;
whisker_stim_tmax = squeeze(max(whisker_stim_avg(:,use_t,:,:),[],2));
plot_days = min(sum(~isnan(whisker_stim_tmax),3)) >= min_n;
figure; hold on
set(gca,'ColorOrder',stim_col);
errorbar(repmat(learned_day_unique(plot_days),1,length(stim_unique)), ...
    nanmean(whisker_stim_tmax(:,plot_days,:),3)', ...
    AP_sem(whisker_stim_tmax(:,plot_days,:),3)','linewidth',2,'CapSize',0);
axis tight;
xlim(xlim+[-0.5,0.5]);
xline(0,'color','k','linestyle','--');
xlabel('Learned day');
ylabel('Whisker movement');

% Plot fluorescence/whisker by discretized whisker movement
plot_rois = [6];
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
        xlim([-0.2,0.7]);
        xline(0);
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
            xlim([-0.2,0.7]);
            xline(0);
            xlabel('Time from stim (s)');
            ylabel('\DeltaF/F0');
        end
    end
    ax(:,curr_stim_idx) = flipud(allchild(h));
    title(h,sprintf('Whisker move groups (stim %d)',curr_stim));

    % Plot line plot of fluorescence vs whisker
    stage_col = [0.5,0.5,0.5;0,0,0];
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
                'linewidth',2,'capsize',0,'color', ...
                min(1,stim_col(curr_stim_idx,:) + stage_col(curr_stage,:)));
            
            xlabel(curr_ax,'Whisker movement');
            ylabel(curr_ax,'\DeltaF/F_0');
            title(curr_ax,wf_roi(curr_roi).area);
        end
    end
end
% (link modality axes across figures)
for i = 1:length(plot_rois)+1
    linkaxes(ax(i:length(plot_rois)+1:end,:),'xy');
    axes(ax(i,end));
    if i == 1
        y_scale = 0.2;
    else
        y_scale = 4e-4;
    end
    x_scale = 0.2;
    AP_scalebar(x_scale,y_scale);
    drawnow;
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


% Stats: ANOVA
for curr_stim = stim_unique'
    use_trials = trial_stim_allcat == curr_stim & quiescent_trials & ...
        ~all(isnan(whisker_allcat),2) & ~isnan(whisker_grp);

    fluor_roi_whiskergrp_avg = accumarray( ...
        reshape(accum_whisker_idx(use_trials,:,:,:),[],size(accum_whisker_idx,4)), ...
        reshape(fluor_roi_deconv(use_trials,:,:),[],1),[],@nanmean,NaN('single'));

    curr_data = squeeze(nanmean(fluor_roi_whiskergrp_avg(:,use_t,:,curr_roi,:),2));
    p = anova2(reshape(permute(curr_data,[3,1,2]),[],2),length(animals),'off');

    fprintf('%s stim %d 2-way anova: p(stage) = %.2g\n',wf_roi(curr_roi).area,curr_stim,p(1));
end



%% (diagram: plot widefield ROIs)

wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

plot_rois = reshape([1,6,7]' + [0,size(wf_roi,1)],[],1)';

roi_col = copper(n_rois);
roi_cat = cat(3,wf_roi.mask);

figure; 

subplot(1,2,1); hold on
set(gca,'YDir','reverse');
AP_reference_outline('ccf_aligned_lefthemi',[0.7,0.7,0.7]);
for curr_roi = plot_rois
    curr_roi_boundary = cell2mat(bwboundaries(roi_cat(:,:,curr_roi)));
    patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),roi_col(curr_roi,:), ...
        'EdgeColor','none');
end
axis image off;

subplot(1,2,2); hold on;
set(gca,'YDir','reverse');
AP_reference_outline('ccf_aligned',[0.7,0.7,0.7]);
for curr_roi = plot_rois
    curr_roi_boundary = cell2mat(bwboundaries(roi_cat(:,:,curr_roi)));
    patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),roi_col(curr_roi,:), ...
        'EdgeColor','none');
end
axis image off;



%% ------- REVISION: IN PROGRESS -----------------------

%% ++ Long-term post-learning widefield

% Get training data from subset of long-term animals
long_term_animals = {'AP113','AP114','AP115'};
use_animals = ismember(animals,long_term_animals);


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

% Get pixels and pixel timemax by stage
stim_px_avg_stage = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),stim_v_avg_stage);

use_t = t >= 0 & t <= 0.2;
stim_px_avg_stage_tmax = ...
    squeeze(max(stim_px_avg_stage(:,:,use_t,:,:,:),[],3));

% Set ROIs to plot

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

roi_learnday_avg = accumarray( ...
    reshape(accum_idx(use_trials_training,:,:,[2,3,4,5,6]),[],5), ...
    reshape(fluor_roi_deconv(use_trials_training,:,:),[],1), ...
    [length(learned_day_unique),length(t),n_rois,length(stim_unique),length(animals)], ...
    @nanmean,NaN('single'));

% Keep pixels and ROI data from training, then clear
plot_stim = 1;

% (pixels)
px_training_stage = nanmean(stim_px_avg_stage_tmax(:,:,:, ...
    stim_unique == plot_stim,use_animals),5);
% (ROI)
plot_roi = 6;
use_t = t > 0 & t <= 0.2;
mpfc_training_stage = squeeze(max(cat(1,...
    nanmean(roi_learnday_avg(learned_day_unique < 0,use_t, ...
    plot_roi,stim_unique == plot_stim,use_animals),1), ...
    nanmean(roi_learnday_avg(learned_day_unique >= 0,use_t, ...
    plot_roi,stim_unique == plot_stim,use_animals),1)),[],2));

clearvars -except px_training_stage mpfc_training_stage


%%% LOAD POST-TRAINING DATA

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_passive_teto_postlearn';
AP_load_trials_operant;

% Get animal and day index for each trial
trial_animal = cell2mat(arrayfun(@(x) ...
    x*ones(size(vertcat(wheel_all{x}{:}),1),1), ...
    [1:length(wheel_all)]','uni',false));
trial_day = cell2mat(cellfun(@(x) cell2mat(cellfun(@(curr_day,x) ...
    curr_day*ones(size(x,1),1),num2cell(1:length(x))',x,'uni',false)), ...
    wheel_all,'uni',false));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));

% Store day relative to last task day
task_relative_day = cellfun(@(day,taskday) day' - taskday, ...
    recording_day,last_task_right_day,'uni',false)';
retired_relative_day = cellfun(@(day,taskday) day' - taskday, ...
    recording_day,last_task_left_day,'uni',false)';
trial_postlearn_day = cell2mat(cellfun(@(x,day) cell2mat(cellfun(@(x,day) ...
    day*ones(size(x,1),1),x,num2cell(day),'uni',false)), ...
    wheel_all,task_relative_day,'uni',false));

% Get trials with movement during stim to exclude
quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > 0,2);

% Turn values into IDs for grouping
stim_unique = unique(trial_stim_allcat);
[~,trial_stim_id] = ismember(trial_stim_allcat,stim_unique);

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

%%% Get post-learning pixels and ROI

% Get relative post-learned days
% (all animals imaged in same intervals, so unique should be same as length)
task_relative_day_unique = unique(cat(1,task_relative_day{:}));
retired_relative_day_unique = unique(cat(1,retired_relative_day{:}));
if ~all(task_relative_day_unique == cat(2,task_relative_day{:}),'all')
    error('Different relative days across animals');
end

% Get average image by day 
use_stim = stim_unique == 1;
use_t = t >= 0 & t <= 0.2;

stim_v_dayavg = nanmean(cat(5,stim_v_avg{:}),5);
stim_px_dayavg = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),stim_v_dayavg(:,:,:,use_stim));

stim_px_dayavg_tmax = squeeze(max(stim_px_dayavg(:,:,use_t,:),[],3));

% Get average ROI by day
stim_roi_avg_cat = cat(5,stim_roi_avg{:});
stim_roi_day = squeeze(stim_roi_avg_cat(:,:,:,use_stim,:));

% Get ROI activity within stim window
use_t = t > 0 & t <= 0.2;
stim_roi_avg_tmax = squeeze(max(stim_roi_avg_cat(:,use_t,:,:,:),[],2));


% Grab pixels and ROI data from post-training, concatenate to pre-training
plot_stim = 1;

% (pixels)
px_posttraining_stage = stim_px_dayavg_tmax;
% (ROI)
plot_roi = 6;
mpfc_posttraining_stage = squeeze(stim_roi_avg_tmax(plot_roi,:,stim_unique == plot_stim,:));

px_prepost_training = cat(3,px_training_stage,px_posttraining_stage);
mpfc_prepost_training = cat(1,mpfc_training_stage,mpfc_posttraining_stage);

% Plot pixels and ROIs across pre/post-training
prepost_training_labels = [{'Novice','Trained', ...
    sprintf('%d days from R task',task_relative_day_unique(1))}, ...
    cellfun(@(x) sprintf('%d days devalued',x), ...
    num2cell(retired_relative_day_unique(2:end))','uni',false)];
figure;
h = tiledlayout(1,size(px_prepost_training,3));
c = [0,0.003];
for curr_stage = 1:size(px_prepost_training,3)
        nexttile;
        imagesc(px_prepost_training(:,:,curr_stage));
        axis image off;
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        colormap(gca,AP_colormap('WG',[],1.5));
        caxis(c);

        title(prepost_training_labels{curr_stage});
end
linkaxes(allchild(h),'xy');
colorbar;

figure; hold on;
plot(mpfc_prepost_training,'color',[0.5,0.5,0.5]);
plot(nanmean(mpfc_prepost_training,2),'color','k','linewidth',2);
ylabel(wf_roi(plot_roi).area);
set(gca,'XTick',1:length(prepost_training_labels), ...
    'XTickLabels',prepost_training_labels);



%% ++ Reaction time by activity plot
% (this is a bit messy - clean this up)

% Load behavior, exclude animals not in dataset, get learned day
data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
bhv_fn = [data_path filesep 'bhv_teto'];
load(bhv_fn);
bhv = bhv(ismember({bhv.animal},animals)); 
learned_day = cellfun(@(x) find(x,1),{bhv.learned_days})';
learned_day_animal = cellfun(@(ld,n) [1:n]'-(ld), ...
    num2cell(learned_day),num2cell(cellfun(@length,trial_info_all)), ...
    'uni',false);

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
max_days = max(cellfun(@sum,use_days));

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


% Get average ROI activity 

% Get indicies for averaging 
% (adjust trial day to only include training days)
trial_day_training = trial_day - n_naive;

[trained_day_idx,t_idx,roi_idx] = ...
    ndgrid(trial_day_training,1:length(t),1:n_rois);

[learned_day_idx,~,~] = ...
    ndgrid(trial_learned_day_id,1:length(t),1:n_rois);
[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_rois);
[stim_idx,~,~] = ...
    ndgrid(trial_stim_id,1:length(t),1:n_rois);

accum_idx = cat(4,trained_day_idx,learned_day_idx,t_idx,roi_idx,stim_idx,animal_idx);

use_trials_training = quiescent_trials & trial_day_training > 0;

roi_trainday_avg = accumarray( ...
    reshape(accum_idx(use_trials_training,:,:,[1,3,4,5,6]),[],5), ...
    reshape(fluor_roi_deconv(use_trials_training,:,:),[],1), ...
    [max(trial_day_training),length(t),n_rois,length(stim_unique),length(animals)], ...
    @nanmean,NaN('single'));

% Get ROI activity within stim window
use_t = t > 0 & t <= 0.2;
stim_roi_trainday_avg_tmax = squeeze(max(roi_trainday_avg(:,use_t,:,:,:),[],2));


% Plot activity against reaction times
plot_roi = 6;
plot_stim = 1;

plot_rxn = (nanmean(rxn_alt_med,3)-rxn_measured_med)./ ...
    (nanmean(rxn_alt_med,3)+rxn_measured_med);
plot_act = permute(stim_roi_trainday_avg_tmax(:,plot_roi, ...
    stim_unique == plot_stim,:),[1,4,2,3]);

animal_col = max(0,jet(length(animals))-0.2);

figure; 
% (plot single mouse example - lines and scatter)
example_animal = 1;

subplot(1,3,1); hold on; axis square;
yyaxis left;plot(plot_rxn(:,example_animal),'k','linewidth',2);
axis tight;
ylabel('Reaction time idx');
yyaxis right;plot(plot_act(:,example_animal),'r','linewidth',2);
axis tight;
ylabel(wf_roi(plot_roi).area);
title(animals(example_animal));
xlabel('Day');

subplot(1,3,2); hold on; axis square;
plot(plot_rxn(:,example_animal),plot_act(:,example_animal), ...
    '.-','MarkerSize',15,'color',animal_col(example_animal,:),'linewidth',2);
xlabel('Reaction time idx');
ylabel(wf_roi(plot_roi).area);
title(animals(example_animal));

% (plot all mice together)
subplot(1,3,3); hold on; axis square;
set(gca,'ColorOrder',animal_col);
plot(plot_rxn,plot_act,'.','MarkerSize',15);
xlabel({'Reaction time idx','(alt-meas/alt+meas)'})
ylabel(wf_roi(plot_roi).area);
axis tight;
xlim(xlim + 0.1*range(xlim).*[-1,1]);
ylim(ylim + 0.1*range(ylim).*[-1,1]);

% Plot all mice separately
figure; h = tiledlayout('flow');
for curr_animal = 1:length(animals)
    nexttile;
    yyaxis left;plot(plot_rxn(:,curr_animal),'k','linewidth',2);
    axis tight; ylim(ylim + 0.1*range(ylim).*[-1,1]);
    ylabel('Reaction time idx');
    yyaxis right;plot(plot_act(:,curr_animal),'r','linewidth',2);
    axis tight; ylim(ylim + 0.1*range(ylim).*[-1,1]);
    ylabel(wf_roi(plot_roi).area);
    title(animals(curr_animal));
    xlabel('Day'); 
    xlim(xlim+[-1,1]);
end


%% ITI turns: quantify frequency of turning with no stim










%% Naive ephys - plot probe position

animals = {'AP116','AP117','AP118','AP119'};

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
scalebar_length = 1000/10; % um/voxel size
line(ccf_axes(3),[0,0],[0,scalebar_length],'color','m','linewidth',3);


%% Naive ephys - naive stim response

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_passive_ephys_naive';

AP_load_trials_operant;

% Get animal and day index for each trial
trial_animal = cell2mat(arrayfun(@(x) ...
    x*ones(size(vertcat(wheel_all{x}{:}),1),1), ...
    [1:length(wheel_all)]','uni',false));
trial_day = cell2mat(cellfun(@(x) cell2mat(cellfun(@(curr_day,x) ...
    curr_day*ones(size(x,1),1),num2cell(1:length(x))',x,'uni',false)), ...
    wheel_all,'uni',false));

trial_recording = cell2mat(cellfun(@(tr,day) ...
    day*ones(size(tr,1),1), ...
    cat(1,wheel_all{:}),num2cell(1:length(cat(1,wheel_all{:})))', ...
    'uni',false));

trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));

% Get trials with movement during stim to exclude
quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > 0,2);

% Get average response timecourses
stim_unique = unique(trial_stim_allcat);
[~,trial_stim_id] = ismember(trial_stim_allcat,stim_unique);

use_trials = quiescent_trials;

% (loop through cell types)
for curr_celltype = 1:size(mua_area_allcat,4)

    [recording_idx,t_idx,area_idx] = ...
        ndgrid(trial_recording(use_trials),1:length(t),1:length(mua_areas));
    [stim_idx,~,~,] = ...
        ndgrid(trial_stim_id(use_trials),1:length(t),1:length(mua_areas));

    mua_recording_avg = accumarray([recording_idx(:),t_idx(:),stim_idx(:),area_idx(:)], ...
        reshape(mua_area_allcat(use_trials,:,:,curr_celltype),[],1), ...
        [length(trials_recording),length(t),length(stim_unique),length(mua_areas)], ...
        @nanmean,NaN);

    % Get DV position of each area in CCF for sorting
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

    % (get areas to plot: anything present in all recordings)
    [~,area_idx] = cellfun(@(x) ismember(x,mua_areas),mua_areas_cat,'uni',false);
    area_recording_n = accumarray(cell2mat(area_idx),1);
    plot_areas = find(area_recording_n == length(trials_recording));
    [~,plot_area_sort_idx] = sort(mua_areas_dvmin(plot_areas));

    % Plot stim overlaid
    figure;
    stim_color = [0,0,0.8;0.5,0.5,0.5;0.8,0,0];
    h = tiledlayout(length(plot_areas),1);
    for curr_area = plot_areas(plot_area_sort_idx)'
        nexttile;
        AP_errorfill(t', ...
            squeeze(nanmean(mua_recording_avg(:,:,:,curr_area),1)), ...
            squeeze(AP_sem(mua_recording_avg(:,:,:,curr_area),1)),stim_color);
        yline(0);
        ylabel(mua_areas(curr_area));
    end
    linkaxes(allchild(h),'xy');
    % (shade stim area and put in back)
    arrayfun(@(x) patch(x,[0,0.5,0.5,0], ...
        reshape(repmat(ylim(x),2,1),[],1),[1,1,0.8], ...
        'linestyle','none'),allchild(h));
    arrayfun(@(x) set(x,'children',circshift(get(x,'children'),-1)),allchild(h));
    xlim([-0.2,0.7]);

    x_scale = 0.2;
    y_scale = 0.5;
    AP_scalebar(x_scale,y_scale);

    title(h,sprintf('Cell type %d',curr_celltype));
end

%% --> (USE THIS ALT) Recording locations in naive and trained

% Set animal groups (naive, trained)
animal_groups = {{'AP116','AP117','AP118','AP119'}, ...
    {'AP100','AP101','AP104','AP105','AP106'}};

animal_group_col = [0.5,0.5,0.5;0,0,0];

% Load CCF annotated volume
allen_atlas_path = fileparts(which('template_volume_10um.npy'));
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

% Set up 3D axes
figure;
ccf_3d_axes = axes;
[~, brain_outline] = plotBrainGrid([],ccf_3d_axes);
set(ccf_3d_axes,'ZDir','reverse');
hold(ccf_3d_axes,'on');
axis vis3d equal off manual
view([-30,25]);
axis tight;
h = rotate3d(ccf_3d_axes);
h.Enable = 'on';

% Set up 2D axes
figure;
bregma_ccf = [540,44,570];
ccf_size = size(av);

ccf_axes = gobjects(3,1);
ccf_axes(1) = subplot(1,3,1,'YDir','reverse');
hold on; axis image off;
ccf_axes(2) = subplot(1,3,2,'YDir','reverse');
hold on; axis image off;
ccf_axes(3) = subplot(1,3,3,'YDir','reverse');
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

for animal_group = 1:length(animal_groups)

    animals = animal_groups{animal_group};
    animal_col = repmat(animal_group_col(animal_group,:),length(animals),1);

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
end

% Plot scalebar
scalebar_length = 1000/10; % um/voxel size
line(ccf_axes(3),[0,0],[0,scalebar_length],'color','m','linewidth',3);


%% --> (USE THIS ALT) Naive and trained passive stim ephys

mua_stim_stage = cell(2,1);

for curr_stage = 1:2 % (naive, trained)

    % Load data
    trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    switch curr_stage
        case 1
            data_fn = 'trial_activity_passive_ephys_naive';
        case 2
            data_fn = 'trial_activity_passive_ephys';
    end

    AP_load_trials_operant;

    % Get animal and day index for each trial
    trial_animal = cell2mat(arrayfun(@(x) ...
        x*ones(size(vertcat(wheel_all{x}{:}),1),1), ...
        [1:length(wheel_all)]','uni',false));
    trial_day = cell2mat(cellfun(@(x) cell2mat(cellfun(@(curr_day,x) ...
        curr_day*ones(size(x,1),1),num2cell(1:length(x))',x,'uni',false)), ...
        wheel_all,'uni',false));

    trial_recording = cell2mat(cellfun(@(tr,day) ...
        day*ones(size(tr,1),1), ...
        cat(1,wheel_all{:}),num2cell(1:length(cat(1,wheel_all{:})))', ...
        'uni',false));

    trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));

    % Get trials with movement during stim to exclude
    quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > 0,2);

    % Get average response timecourses
    stim_unique = unique(trial_stim_allcat);
    [~,trial_stim_id] = ismember(trial_stim_allcat,stim_unique);

    use_trials = quiescent_trials;

    [recording_idx,t_idx,area_idx] = ...
        ndgrid(trial_recording(use_trials),1:length(t),1:length(mua_areas));
    [stim_idx,~,~,] = ...
        ndgrid(trial_stim_id(use_trials),1:length(t),1:length(mua_areas));

    mua_recording_avg = accumarray([recording_idx(:),t_idx(:),stim_idx(:),area_idx(:)], ...
        reshape(mua_area_allcat(use_trials,:,:),[],1), ...
        [length(trials_recording),length(t),length(stim_unique),length(mua_areas)], ...
        @nanmean,NaN);

    % Set areas and order to plot
    % (used to do this programatically to be agnostic, but not necessary
    % now so just hard-coding)
    plot_areas = {'Secondary motor area', ...
        'Anterior cingulate area dorsal part', ...
        'Prelimbic area','Infralimbic area'};
    [~,plot_areas_idx] = ismember(plot_areas,mua_areas);

    % Store activity
    mua_stim_stage{curr_stage} = mua_recording_avg(:,:,:,plot_areas_idx);

    % Clear for next round
    clearvars -except t plot_areas mua_stim_stage

end

% Plot stim overlaid
plot_stim = 3;
stage_cols = [0.5,0.5,0.5;0,0,0];
figure;
for curr_stage = 1:2
    for curr_area = 1:length(plot_areas)
        subplot(length(plot_areas),1,curr_area); hold on
        AP_errorfill(t', ...
            squeeze(nanmean(mua_stim_stage{curr_stage}(:,:,plot_stim,curr_area),1)), ...
            squeeze(AP_sem(mua_stim_stage{curr_stage}(:,:,plot_stim,curr_area),1)), ...
            stage_cols(curr_stage,:));
        yline(0);xline(0);
        ylabel(plot_areas{curr_area});
    end
end
linkaxes(get(gcf,'Children'),'xy')



