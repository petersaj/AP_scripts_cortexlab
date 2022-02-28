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

max_days = max(cellfun(@sum,use_days));

% Set bins for reaction time histograms
rxn_bins = [0:0.01:1];
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


% Total (single day - 1:1 rxn)
plot_day = 1;

animal_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(vertcat(x{plot_day}),rxn_bins,'normalization','pdf'), ...
    use_days,{bhv.stim_move_t},'uni',false)');

animal_alt_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(cell2mat(vertcat(x{plot_day})),rxn_bins,'normalization','pdf'), ...
    use_days,{bhv.alt_stim_move_t},'uni',false)');

figure; hold on;
histogram('BinEdges',rxn_bins,'BinCounts',nanmean(animal_rxn,1),'EdgeColor','none')
histogram('BinEdges',rxn_bins,'BinCounts',nanmean(animal_alt_rxn,1),'EdgeColor','none')
xlabel('Reaction time');
ylabel('PDF');
legend({'Measured','Null'});
title(sprintf('Day %d',plot_day));

% Total (single animal)
plot_animal = 10;

animal_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(vertcat(x{:}),rxn_bins,'normalization','pdf'), ...
    use_days(plot_animal),{bhv(plot_animal).stim_move_t},'uni',false)');

animal_alt_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(cell2mat(vertcat(x{:})),rxn_bins,'normalization','pdf'), ...
    use_days(plot_animal),{bhv(plot_animal).alt_stim_move_t},'uni',false)');

figure; hold on;
histogram('BinEdges',rxn_bins,'BinCounts',nanmean(animal_rxn,1),'EdgeColor','none')
histogram('BinEdges',rxn_bins,'BinCounts',nanmean(animal_alt_rxn,1),'EdgeColor','none')
xlabel('Reaction time');
ylabel('PDF');
legend({'Measured','Null'});
title(sprintf('Animal %s',animals{plot_animal}));


% Total
animal_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(vertcat(x{use}),rxn_bins,'normalization','pdf'), ...
    use_days,{bhv.stim_move_t},'uni',false)');

animal_alt_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(cell2mat(vertcat(x{use})),rxn_bins,'normalization','pdf'), ...
    use_days,{bhv.alt_stim_move_t},'uni',false)');

figure; hold on;
histogram('BinEdges',rxn_bins,'BinCounts',nanmean(animal_rxn,1),'EdgeColor','none')
histogram('BinEdges',rxn_bins,'BinCounts',nanmean(animal_alt_rxn,1),'EdgeColor','none')
xlabel('Reaction time');
ylabel('PDF');
legend({'Measured','Null (all)'});

% Total (trial/alt number-matched)
alt_rxn_matched = ...
    cellfun(@(x) cellfun(@(x) cell2mat(cellfun(@(x) ...
    datasample(x,min(length(x),1)),x,'uni',false)),x,'uni',false), ...
    {bhv.alt_stim_move_t},'uni',false);

animal_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(vertcat(x{use}),rxn_bins,'normalization','pdf'), ...
    cellfun(@(x) x&([1:length(x)]' > 6),use_days,'uni',false),{bhv.stim_move_t},'uni',false)');

animal_alt_rxn = cell2mat(cellfun(@(use,x) ...
    histcounts(vertcat(x{use}),rxn_bins,'normalization','pdf'), ...
    cellfun(@(x) x&([1:length(x)]' > 6),use_days,'uni',false),alt_rxn_matched,'uni',false)');

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


% Median reaction time (padded)
rxn_median_pad = cell2mat(cellfun(@(x) padarray(cellfun(@(x) ...
    nanmedian(x),x), ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_move_t},use_days,'uni',false),'uni',false));

% Fraction rxn within boundary
max_days = max(cellfun(@sum,use_days));

rxn_frac_window = [0.1,0.26];
rxn_frac_pad = cell2mat(cellfun(@(x) padarray(cellfun(@(x) ...
    nanmean(x > rxn_frac_window(1) & x < rxn_frac_window(2)),x), ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_move_t},use_days,'uni',false),'uni',false));

alt_rxn_frac_pad = cell2mat(cellfun(@(x) padarray(cellfun(@(x) ...
    nanmean(cell2mat(x) > rxn_frac_window(1) & cell2mat(x) < rxn_frac_window(2)),x), ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.alt_stim_move_t},use_days,'uni',false),'uni',false));

figure; hold on;
plot(rxn_frac_pad);
errorbar(nanmean(rxn_frac_pad,2),AP_sem(rxn_frac_pad,2),'k','linewidth',2);
errorbar(nanmean(alt_rxn_frac_pad,2),AP_sem(alt_rxn_frac_pad,2),'r','linewidth',2);

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

figure; hold on;
stim_move_t_stimtime_long = reshape(padarray(stim_move_t_stimtime,[1,0,0],NaN,'post'),[],length(animals));
alt_stim_move_t_stimtime_long = reshape(padarray(alt_stim_move_t_stimtime,[1,0,0],NaN,'post'),[],length(animals));

plot([1:size(stim_move_t_stimtime_long,1)]/(n_daysplit+1),stim_move_t_stimtime_long,'color',[0.5,0.5,0.5]);
p1 = errorbar([1:size(stim_move_t_stimtime_long,1)]/(n_daysplit+1), ...
    nanmean(stim_move_t_stimtime_long,2),AP_sem(stim_move_t_stimtime_long,2),'k','linewidth',2);

plot([1:size(alt_stim_move_t_stimtime_long,1)]/(n_daysplit+1),alt_stim_move_t_stimtime_long,'color',[1,0.5,0.5]);
p2 = errorbar([1:size(alt_stim_move_t_stimtime_long,1)]/(n_daysplit+1), ...
    nanmean(alt_stim_move_t_stimtime_long,2),AP_sem(alt_stim_move_t_stimtime_long,2),'r','linewidth',2);

xlabel('Day');
ylabel('Frac rxn 100-250 ms');
legend([p1,p2],{'Measured','Null'},'location','nw');



