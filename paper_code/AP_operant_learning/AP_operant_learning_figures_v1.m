%% Behavior
% Adapted from AP_operant_behavior

% (load data saved in AP_operant_behavior)
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


% %%%%% TESTING STATS rxn within window
% THIS IS GOOD, USE THIS
max_days = max(cellfun(@sum,use_days));
r_frac_p_cat = nan(max_days,length(animals));
for curr_animal = 1:length(animals)
    for curr_day = find(use_days{curr_animal})'
        
        n_sample = 10000;
        use_trials = ~cellfun(@isempty,bhv(curr_animal).alt_stim_move_t{curr_day});
%         use_trials = ~cellfun(@isempty,bhv(curr_animal).alt_stim_move_t{curr_day}) & ...
%             linspace(0,1,length(bhv(curr_animal).stim_move_t{curr_day}))' < 0.9;
%           use_trials = ~cellfun(@isempty,bhv(curr_animal).alt_stim_move_t{curr_day}) & ...
%             (length(bhv(curr_animal).stim_move_t{curr_day}):-1:1)' > 10;
        
        r = bhv(curr_animal).stim_move_t{curr_day}(use_trials);
        ar = cell2mat(cellfun(@(x) datasample(x,n_sample)', ...
            bhv(curr_animal).alt_stim_move_t{curr_day}(use_trials),'uni',false));
        
        rxn_window = [0.1,0.2];
        
        r_frac = nanmean(r >= rxn_window(1) & r <= rxn_window(2));
        ar_frac = nanmean(ar >= rxn_window(1) & ar <= rxn_window(2),1);
        
        r_frac_rank = tiedrank([r_frac,ar_frac]);
        r_frac_p = r_frac_rank(1)./(n_sample+1);
        
        r_frac_p_cat(curr_day,curr_animal) = r_frac_p;
        
    end
end
figure; 
imagesc(r_frac_p_cat > 0.99,'AlphaData',~isnan(r_frac_p_cat));
set(gca,'Color',[0.5,0.5,0.5])
title('rxn window')
xlabel('Animal');
ylabel('Day');
[~,learned_day] = max(r_frac_p_cat,[],1);


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
% alt_stim_move_t_cat = cell2mat(cellfun(@(x,use_days) cell2mat(x(use_days)), ...
%     alt_rxn_matched,use_days,'uni',false)');
alt_stim_move_t_cat = cell2mat(cellfun(@(x,use_days) ...
    cell2mat(cellfun(@(x) cellfun(@nanmedian,x),x(use_days),'uni',false)), ...
    {bhv.alt_stim_move_t},use_days,'uni',false)');

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
learned_day_grp_daysplit = cellfun(@str2num,learned_day_grp_daysplit);

subplot(1,2,2); hold on;
errorbar(learned_day_grp_daysplit,rxn_grp_mean_daysplit, ...
    rxn_grp_sem_daysplit,'k','linewidth',2,'CapSize',0);
line([0,0],ylim,'color','k','linestyle','--');
xlabel('Learned day');

%% ~~~~~ Load passive

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_passive_tetO';
AP_load_trials_wf;

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

n_naive = 3; % Number of days naive only (no behavior)
passive_learned_day_grp = cellfun(@(v,learned_day) ...
    discretize(1:size(v,3),[0,n_naive+1,learned_day,size(v,3)+1],'IncludedEdge','left'), ...
    stim_v_avg,num2cell(learned_day + n_naive),'uni',false);

% (TO DO HERE: AVERAGE V'S BY NAIVE/PRELEARN/POSTLEARN ABOVE)
stim_v_avg_prepostlearn
stim_v_avg










