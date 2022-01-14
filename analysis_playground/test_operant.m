%% Test analysis for imaging with operant task (over learning)
% (note: this includes things brought from test_corticostriatal and
% test_learning)

%% ~~~~~~~~~ SINGLE-SESSION ~~~~~~~~~

%% Load example session

% animal = 'AP100';
% day = '2021-05-03';
% day = '2021-05-07';

animal = 'AP104';
% day = '2021-06-06';
% day = '2021-06-08';
% day = '2021-06-09';
% day = '2021-06-11';
day = '2021-06-13';

experiment = 1;
verbose = true;
load_parts.imaging = false;
AP_load_experiment;

t = Timeline.rawDAQTimestamps;

%% Get rewarded/unrewarded movement times/epochs

% Plot velocity/reward/stim
wheel_velocity_move = wheel_velocity;
wheel_velocity_move(~wheel_move) = NaN;

figure; hold on;
plot(t,wheel_velocity,'k');
plot(t,wheel_velocity_move,'r');
plot(reward_t_timeline,0,'.b','MarkerSize',10);
plot(t,stimOn_epochs*0.1,'g');

% Get real and null distribution of stim-to-move times
% (get gaps between iti movements and prior movements)
use_wheel_move_iti = wheel_move_iti_idx(wheel_move_iti_idx > 1); % no gap for first movement
iti_move_gaps = wheel_starts(use_wheel_move_iti) - wheel_stops(use_wheel_move_iti - 1);
% (get gaps between stim movements and prior movements)
use_wheel_move_stim = wheel_move_stim_idx(wheel_move_stim_idx > 1); % no gap for first movement
wheel_move_preresponse_stop = wheel_stops(use_wheel_move_stim-1);
% (make null distribution: positive prior movement + gap relative to stim)
% (from 2nd stim, in case no move before first)
n_shuff = 10000;
stim_to_move_shuff = ...
    stimOn_times(end-length(use_wheel_move_stim)+1:end) ...
    - (wheel_move_preresponse_stop + ...
    reshape(randsample(iti_move_gaps,length(wheel_move_preresponse_stop)*n_shuff,true),[],n_shuff));
stim_to_move_shuff(stim_to_move_shuff < 0) = NaN; % positives only, negatives would've been resets

% (TESTING: get move start time after and end time before reward)
iti_move_gaps = wheel_starts(wheel_move_iti_idx(2:length(wheel_move_iti_idx))) - ...
    wheel_starts(wheel_move_iti_idx(2:length(wheel_move_iti_idx))-1);

wheel_move_preresponse_stop = wheel_stops(wheel_move_stim_idx(2:length(wheel_move_stim_idx))-1);

r = wheel_starts(wheel_move_stim_idx(2:end)) - wheel_move_preresponse_stop;


a = wheel_starts(wheel_move_stim_idx(2:length(wheel_move_stim_idx))-1);
b = wheel_stops(wheel_move_stim_idx(2:length(wheel_move_stim_idx))-2);
c = a - b;
figure; hold on;
histogram(stim_to_move_shuff(:),0:0.02:1,'EdgeColor','none','normalization','probability')
histogram(stim_to_move,0:0.02:1,'EdgeColor','none','normalization','probability')


% (align movement to prestim wheel stop)
figure;subplot(1,2,1);
% use_align = wheel_move_preresponse_stop;
use_align = wheel_stops(wheel_move_stim_idx-2);
surround_t_centers = -10:0.1:10;
surround_times = use_align + surround_t_centers;
surround_move = interp1(t,wheel_move,surround_times,'previous');
imagesc(surround_move);
colormap(brewermap([],'Greys'));

% (align movement to second-last stop vs stim)
subplot(1,2,2);
use_align = stimOn_times;
surround_t_centers = -10:0.1:10;
surround_times = use_align + surround_t_centers;
surround_move = interp1(t,wheel_move,surround_times,'previous');
imagesc(surround_move);
colormap(brewermap([],'Greys'));


%% Passive PSTH

% Get quiescent trials
wheel_window = [0,0.5];
wheel_window_t = wheel_window(1):1/Timeline.hw.daqSampleRate:wheel_window(2);
wheel_window_t_peri_event = bsxfun(@plus,stimOn_times,wheel_window_t);
event_aligned_wheel = interp1(Timeline.rawDAQTimestamps, ...
    +wheel_move,wheel_window_t_peri_event,'previous');
quiescent_trials = ~any(abs(event_aligned_wheel) > 0,2);

% PSTH viewer
AP_cellraster(stimOn_times(quiescent_trials),stimIDs(quiescent_trials))


%% Operant PSTH

% Standard events
[~,rxn_sort_idx] = sort(stim_to_move);
AP_cellraster({stimOn_times,wheel_starts(wheel_move_response_idx),reward_t_timeline},rxn_sort_idx)

% PSTH with rewarded & ITI move starts
[~,rxn_sort_idx] = sort(stim_to_move);
[~,iti_move_sort_idx] = sort(cellfun(@(x) sum(abs(x)),t_movesplit(iti_move_epochs)));

AP_cellraster( ...
    {stimOn_times,wheel_starts(wheel_move_response_idx),wheel_starts(wheel_move_iti_idx),reward_t_timeline}, ...
    {rxn_sort_idx,rxn_sort_idx,iti_move_sort_idx,1:length(reward_t_timeline)})


%% Get average rewarded and other movements
% (this is old - better ways to get rewarded/other movements now)

wheel_split_n = cellfun(@length,wheel_velocity_split);
wheel_pos_split = mat2cell(wheel_position,wheel_split_n);

wheel_split_n = cellfun(@length,wheel_velocity_split)';
split_move = cellfun(@(x) x(1),mat2cell(wheel_move,wheel_split_n));

wheel_gain = 8; % (deg/mm - in signals protocol)

max_move = max(wheel_split_n(split_move == 1));
vel_pad = cell2mat(cellfun(@(x) padarray(x,max_move-length(x),NaN,'post'), ...
    wheel_velocity_split(split_move == 1),'uni',false)');
pos_pad_raw = cell2mat(cellfun(@(x) padarray(x,max_move-length(x),NaN,'post'), ...
    wheel_pos_split(split_move == 1),'uni',false)'); 
pos_pad = (pos_pad_raw - pos_pad_raw(1,:))*wheel_gain;

rewarded_epochs = cellfun(@(t) any(ismember(t,reward_t_timeline)), ...
    mat2cell(Timeline.rawDAQTimestamps',wheel_split_n));

rewarded_movements = rewarded_epochs(split_move == 1);

% Plot average rewarded and non-rewarded movements
figure; 

subplot(2,1,1);hold on
plot(nanmean(pos_pad(:,~rewarded_movements),2),'k');
plot(nanmean(pos_pad(:,rewarded_movements),2),'b');
legend({'Unrewarded','Rewarded'});
xlabel('Time from movement onset');
ylabel('Wheel position');
xlim([0,1000]);

subplot(2,1,2); hold on
plot(nanmean(vel_pad(:,~rewarded_movements),2),'k');
plot(nanmean(vel_pad(:,rewarded_movements),2),'b');
legend({'Unrewarded','Rewarded'});
xlabel('Time from movement onset');
ylabel('Wheel velocity');
xlim([0,1000]);

%% Stim-surround movement

t = Timeline.rawDAQTimestamps';

% (probability of any movement around stim)
stim_surround_t_centers = -10:0.1:10;
stim_surround_times = stimOn_times + stim_surround_t_centers;
stim_surround_move = interp1(t,+wheel_move,stim_surround_times,'previous');

figure;
subplot(1,3,1);
imagesc(stim_surround_t_centers,[],stim_surround_move);
line([0,0],ylim,'color','b')
ylabel('Trial order');

subplot(1,3,2);
trial_quiescence = signals_events.trialQuiescenceValues(1:n_trials);
[~,sort_idx] = sort(trial_quiescence);
imagesc(stim_surround_t_centers,[],stim_surround_move(sort_idx,:));
hold on;
plot(-trial_quiescence(sort_idx),1:length(sort_idx),'.r')
plot(stim_to_feedback(sort_idx),1:length(sort_idx),'.c')
line([0,0],ylim,'color','b')
ylabel('Sorted by quiescence')

subplot(1,3,3);
[~,sort_idx] = sort(stim_to_move);
imagesc(stim_surround_t_centers,[],stim_surround_move(sort_idx,:));
hold on;
plot(-trial_quiescence(sort_idx),1:length(sort_idx),'.r')
plot(stim_to_feedback(sort_idx),1:length(sort_idx),'.c')
line([0,0],ylim,'color','b')
ylabel('Sorted by rxn')

colormap(brewermap([],'Greys'));

%% Move-surround movement (and stim)

t = Timeline.rawDAQTimestamps;

% wheel start surround
surround_t_centers = -10:0.1:10;
surround_times = wheel_starts + surround_t_centers;
surround_move = interp1(t,wheel_move,surround_times,'previous');

figure;
subplot(1,2,1);
imagesc(surround_t_centers,[],surround_move(wheel_move_iti_idx,:));
title('ITI movements')

subplot(1,2,2);
imagesc(surround_t_centers,[],surround_move(wheel_move_stim_idx,:));
hold on;
plot(-stim_to_move,1:length(stim_to_move),'.r');
title('Stim movements')

colormap(brewermap([],'Greys'));

% wheel stop surround
surround_t_centers = 0:0.1:5;
surround_times = wheel_stops + surround_t_centers;
surround_move = interp1(t,wheel_move,surround_times,'previous');

figure;
subplot(1,2,1);
imagesc(surround_t_centers,[],surround_move(wheel_move_iti_idx(2:end)-1,:));
title('ITI movements')

subplot(1,2,2);
imagesc(surround_t_centers,[],surround_move(wheel_move_stim_idx-1,:));
hold on;
title('Stim movements')

colormap(brewermap([],'Greys'));

% Outcome surround
surround_t_centers = -2:0.1:60;
surround_times = signals_events.responseTimes' + surround_t_centers;
surround_move = interp1(t,+wheel_move,surround_times,'previous');

response_times = signals_events.responseTimes';
response_to_stim = stimOn_times(2:n_trials)-response_times(1:n_trials-1);
[~,sort_idx] = sort(response_to_stim);

figure;
imagesc(surround_t_centers,[],surround_move(sort_idx,:)); hold on;
plot(response_to_stim(sort_idx)+surround_t_centers(1),1:length(sort_idx),'r','linewidth',2)
line([0,0],ylim,'color','blue','linewidth',2)
xlabel('Time from outcome')
ylabel('Trial (sorted by stim time')
colormap(brewermap([],'Greys'));





%% Move time regression (Cox proportional hazards)

t = Timeline.rawDAQTimestamps';

% Plot velocity/reward/stim
wheel_velocity_move = wheel_velocity;
wheel_velocity_move(~wheel_move) = NaN;

figure; hold on;
plot(t,wheel_velocity,'k');
plot(t,wheel_velocity_move,'r');
plot(reward_t_timeline,0,'.b','MarkerSize',10);
plot(t,stimOn_epochs*0.1,'g');


% (carried over from AP_operant_behavior)
t_since_reward = t - ...
    interp1(reward_t_timeline,reward_t_timeline,t,'previous','extrap');

t_since_move = t - ...
    interp1(wheel_stops,wheel_stops,t,'previous','extrap');

t_since_stim = t - ...
    interp1(stimOn_times,stimOn_times,t,'previous','extrap');

% r = [t_since_reward;t_since_move;t_since_stim];
% s = diff([0;wheel_move(1:1000:end)]');s = s == 1;
% t_lags = [-500:500];
% [k,ps,ev] = AP_regresskernel(r(:,1:1000:end),s,t_lags,0);
% figure;plot(t_lags,k');

a = wheel_starts(2:end)-wheel_stops(1:end-1);
b = a(wheel_move_stim_idx-1);

r_idx = (wheel_move_stim_idx-1)+[-10:-1];
use_r = all(r_idx>0,2);

% [k,ps,ev] = AP_regresskernel(a(r_idx(use_r,:)),b(use_r),0,0,[],5,0,1);

x = [a(r_idx(use_r,:)),stim_to_move(use_r),ones(sum(use_r),1)]\b(use_r);


a = diff([0;wheel_move]) == 1;
ar = t_since_reward(a);
am = t_since_move(a);
as = t_since_stim(a);
figure; hold on;
plot3(ar,am,as,'.k')
plot3(ar(wheel_move_stim_idx),am(wheel_move_stim_idx),as(wheel_move_stim_idx),'or')

% Regress time from stim to move relative to previous gaps
use_last_gaps = 5;
wheel_gaps = wheel_starts(2:end)-wheel_stops(1:end-1);
wheel_gaps_regressor_idx = (wheel_move_stim_idx-1)+[-use_last_gaps:-1];
use_move_regressors = all(wheel_gaps_regressor_idx>0,2);
wheel_gaps_regressors = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors,:));

[b,logl,H,stats] = coxphfit(wheel_gaps_regressors,stim_to_move(use_move_regressors));
figure;errorbar(stats.beta,stats.se);
xlabel('Previous move gap');
ylabel('Cox \beta');


% Regress time from move to post-stim move relative to pre-stim gap/stim
use_last_gaps = 5;
wheel_gaps = wheel_starts(2:end)-wheel_stops(1:end-1);
wheel_gaps_regressor_idx = (wheel_move_stim_idx-1)+[-use_last_gaps:0];
use_move_regressors = all(wheel_gaps_regressor_idx>0,2);

wheel_gap_predict = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors,end));
wheel_gaps_regressors = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors,1:end-1));
wheel_stim_regressor = stimOn_times(use_move_regressors) - wheel_stops(wheel_move_stim_idx(use_move_regressors)-1);

[b,logl_stim,H,stats] = coxphfit([wheel_gaps_regressors,wheel_stim_regressor],wheel_gap_predict);
figure;errorbar(stats.beta,stats.se);
xlabel('Previous move gap');
ylabel('Cox \beta');

[b_nostim,logl_nostim,H_nostim,stats_nostim] = coxphfit([wheel_gaps_regressors],wheel_gap_predict);

disp(logl_stim-logl_nostim);

% Regress gap time to stim move with/without stim time
use_last_gaps = 5;
wheel_gaps = [NaN;wheel_starts(2:end)-wheel_stops(1:end-1)];
wheel_gaps_regressor_idx = (1:length(wheel_gaps))'+[-use_last_gaps:0];
use_move_regressors = find(all(wheel_gaps_regressor_idx>0,2));

use_move = intersect(use_move_regressors,wheel_move_stim_idx);

wheel_gap_predict = wheel_gaps(wheel_gaps_regressor_idx(use_move,end));
wheel_gaps_regressors = wheel_gaps(wheel_gaps_regressor_idx(use_move,1:end-1));

wheel_stop_to_stim = stimOn_times - wheel_stops(wheel_move_stim_idx-1);
use_stim_move_idx = ismember(wheel_move_stim_idx,use_move);

q = signals_events.trialQuiescenceValues(1:n_trials)';

% % [b,logl,H,stats] = coxphfit([wheel_gaps_regressors],wheel_gap_predict);
% [b,logl_stim,H,stats] = coxphfit([wheel_gaps_regressors,wheel_stop_to_stim(use_stim_move_idx)],wheel_gap_predict);
% n_shuff = 1000;
% logl_stim_shuff = nan(n_shuff,1);
% for curr_shuff = 1:n_shuff
%     [~,logl_stim_shuff(curr_shuff)] = ...
%         coxphfit([wheel_gaps_regressors,AP_shake(wheel_stop_to_stim(use_stim_move_idx))],wheel_gap_predict);
% end
% [b,logl,H,stats] = coxphfit([wheel_gaps_regressors],wheel_gap_predict);
[b,logl_stim,H,stats] = coxphfit([wheel_gaps_regressors,q(use_stim_move_idx)],wheel_gap_predict);
n_shuff = 1000;
logl_stim_shuff = nan(n_shuff,1);
for curr_shuff = 1:n_shuff
    [~,logl_stim_shuff(curr_shuff)] = ...
        coxphfit([wheel_gaps_regressors,AP_shake(q(use_stim_move_idx))],wheel_gap_predict);
end
figure;histogram(logl_stim_shuff)
line([logl_stim,logl_stim],ylim,'color','r')
figure;errorbar(stats.beta,stats.se);
xlabel('Previous move gap');
ylabel('Cox \beta');
% (does this measure make sense?)
logl_increase_frac = (nanmedian(logl_stim_shuff)-logl_stim)/logl_stim;
disp(logl_increase_frac);

% Regress time from move to ANY move (with/without stim info)
use_last_gaps = 10;
wheel_gaps = wheel_starts(2:end)-wheel_stops(1:end-1);
wheel_gaps_regressor_idx = (1:length(wheel_gaps))'+[-use_last_gaps:0];
use_move_regressors = all(wheel_gaps_regressor_idx>0,2);

wheel_gap_predict = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors,end));
wheel_gaps_regressors = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors,1:end-1));

[b,logl,H,stats] = coxphfit([wheel_gaps_regressors],wheel_gap_predict);
figure;errorbar(stats.beta,stats.se);
xlabel('Previous move gap');
ylabel('Cox \beta');

a = stimOn_times - wheel_stops(wheel_move_stim_idx-1);
idx = ismember(2:length(wheel_starts),wheel_move_stim_idx)';
m = zeros(size(idx));
m(idx) = a;

[b_stim,logl_stim,H,stats_stim] = coxphfit([wheel_gaps_regressors,m(use_move_regressors)],wheel_gap_predict);
figure;errorbar(stats_stim.beta,stats_stim.se);
xlabel('Previous move gap');
ylabel('Cox \beta');

disp(logl_stim-logl);


% Regress ITI gap time move only (sanity check?)
use_last_gaps = 10;
wheel_gaps = [NaN;wheel_starts(2:end)-wheel_stops(1:end-1)];
wheel_gaps_regressor_idx = (1:length(wheel_gaps))'+[-use_last_gaps:0];
use_move_regressors = find(all(wheel_gaps_regressor_idx>0,2));

use_move = intersect(use_move_regressors,wheel_move_iti_idx);

wheel_gap_predict = wheel_gaps(wheel_gaps_regressor_idx(use_move,end));
wheel_gaps_regressors = wheel_gaps(wheel_gaps_regressor_idx(use_move,1:end-1));

wheel_starts_t_since_reward = t_since_reward(ismember(t,wheel_starts));
wheel_starts_t_since_reward_regressors = wheel_starts_t_since_reward(use_move);

[b,logl,H,stats] = coxphfit([wheel_gaps_regressors],wheel_gap_predict);
% [b,logl,H,stats] = coxphfit([wheel_gaps_regressors,wheel_starts_t_since_reward_regressors],wheel_gap_predict);
figure;errorbar(stats.beta,stats.se);
xlabel('Previous move gap');
ylabel('Cox \beta');




% (predict ITI movements, test on stim movements)
use_last_gaps = 5;
wheel_gaps = wheel_starts(2:end)-wheel_stops(1:end-1);
wheel_gaps_regressor_idx = (1:length(wheel_gaps))'+[-use_last_gaps:0];
use_move_regressors = all(wheel_gaps_regressor_idx>0,2);

iti_idx = ismember(2:length(wheel_starts),wheel_move_iti_idx)';
stim_idx = ismember(2:length(wheel_starts),wheel_move_stim_idx)';

% -> fit on ITI
r = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors & iti_idx,1:end-1));
s = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors & iti_idx,end));
[k,ps,ev] = AP_regresskernel(r',s',0,0,[],1,1,1);
disp(ev)

% -> test on stim
r_test = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors & stim_idx,1:end-1));
s_test = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors & stim_idx,end));

s_test_pred = [r_test,ones(size(r_test,1),1)]*cell2mat(k);
r2 = 1 - (nansum((s_test-s_test_pred).^2,1)./ ...
        nansum((s_test-nanmean(s_test,1)).^2,1));
disp(r2)




% (predict ITI movements, test on stim movements)
use_last_gaps = 1;
wheel_gaps = wheel_starts(2:end)-wheel_stops(1:end-1);
wheel_gaps_regressor_idx = (1:length(wheel_gaps))'+[-use_last_gaps:0];
use_move_regressors = all(wheel_gaps_regressor_idx>0,2);

iti_idx = ismember(2:length(wheel_starts),wheel_move_iti_idx)';
stim_idx = ismember(2:length(wheel_starts),wheel_move_stim_idx)';

% -> fit on ITI
r = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors & stim_idx,1:end-1));
s = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors & stim_idx,end));
[k,ps,ev] = AP_regresskernel(r',s',0,0,[]);
disp(ev)

% -> test on stim
r_test = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors & stim_idx,1:end-1));
s_test = wheel_gaps(wheel_gaps_regressor_idx(use_move_regressors & stim_idx,end));

s_test_pred = [r_test,ones(size(r_test,1),1)]*cell2mat(k);
r2 = 1 - (nansum((s_test-s_test_pred).^2,1)./ ...
        nansum((s_test-nanmean(s_test,1)).^2,1));
disp(r2)

%% Reaction time if quiescence time was less than it really was
% (this is a bit weird but gives what looks like a decent measure)

t = Timeline.rawDAQTimestamps';

% Get trace giving time from last movement 
% (processed wheel movement)
t_from_move_trace = ...
    t - interp1(t(wheel_move),t(wheel_move),t,'previous','extrap');

% % (trying to re-create quiescence watch in signals)
% % (NOTE: signals doesn't use clicks, it uses cumulative 1mm)
% % (this is dirty but I don't have a better solution at the moment)
% % (seems like quiescencewatch is really unreliable??)
% % (this looked like it was working at some point, but now not at all?)
% t_wheel_block = interp1(block2timeline,timeline2block,block.inputs.wheelMMTimes);
% a = block.inputs.wheelMMValues;
% quiescence_reset = false(size(block.inputs.wheelMMValues));
% i = 1;
% while i < length(block.inputs.wheelMMValues)
%     a(i:end) = a(i:end)-a(i);
%     curr_diff = abs(block.inputs.wheelMMValues-block.inputs.wheelMMValues(i));
%     thresh = 0.95;
%     i = find(curr_diff(i:end) > thresh,1,'first') + (i-1);
%     quiescence_reset(i) = true;
%     AP_print_progress_fraction(i,length(block.inputs.wheelMMValues));
% end
% % (quiescence watch starts on trial onset)
% quiescence_reset_t = sort([t_wheel_block(quiescence_reset & ~isnan(t_wheel_block)),signals_events.trialNumTimes]);
% t_from_quiescence_reset = t - ...
%     interp1(quiescence_reset_t,quiescence_reset_t,t,'previous','extrap');
% t_from_move_trace = t_from_quiescence_reset;

% Get trace when out of ITI period
noniti_flag = logical(interp1( ...
    [0;signals_events.trialNumTimes(1:n_trials)';signals_events.responseTimes(1:n_trials)'], ...
    [0;ones(n_trials,1);zeros(n_trials,1)], ...
    Timeline.rawDAQTimestamps','previous','extrap'));

% % Get when the minimum quiescence would've been hit
% % % (if quiescence was always minimum)
% min_quiescence = min(signals_events.trialQuiescenceValues);
% min_quiescence_flag = t_from_move_trace >= min_quiescence & noniti_flag;

% (if quiescence was specific for each trial - can shorten artifically)
% quiescence_shorten = 0;
quiescenceMin = min(signals_events.trialQuiescenceValues);
quiescence_shorten = arrayfun(@(x) ...
    randsample(quiescenceMin:0.1:max(x-0.1,quiescenceMin),1), ...
    signals_events.trialQuiescenceValues);
trial_quiescence_trace = interp1(signals_events.trialQuiescenceTimes, ...
    signals_events.trialQuiescenceValues-quiescence_shorten,t,'previous','extrap');
min_quiescence_flag = t_from_move_trace >= trial_quiescence_trace & noniti_flag;

% Sanity check plot
wheel_velocity_move = wheel_velocity;
wheel_velocity_move(~wheel_move) = NaN;

figure; hold on;
plot(t,wheel_velocity,'k');
plot(t,wheel_velocity_move,'r');
plot(reward_t_timeline,0,'.b','MarkerSize',10);

plot(t,t_from_move_trace*0.1,'c');
plot(t,noniti_flag.*trial_quiescence_trace*0.1,'color',[0.5,0.5,0.5]);
plot(t,stimOn_epochs*0.1,'color',[0.7,0,0.7],'linewidth',2);
plot(t,min_quiescence_flag*0.1,'color',[0,0.7,0]);
% legend({'Wheel velocity','Wheel move','Reward','t from trace','non-iti','stimOn','min quiesc'});

% Get first time quiescence period is met
min_quiescence_onsets = t(diff([0;min_quiescence_flag]) == 1);
earliest_quiescence_t = arrayfun(@(x) ...
    min_quiescence_onsets(find(min_quiescence_onsets > x,1,'first')), ...
    signals_events.trialNumTimes(1:n_trials)');

% Get the first movement from min quiescence would have been
earliest_quiescence_wheelstarts = arrayfun(@(x) ...
    wheel_starts(find(wheel_starts > x,1,'first')), ...
    earliest_quiescence_t);

earliest_rxn = earliest_quiescence_wheelstarts - earliest_quiescence_t;

% Plot real and earliest reaction times
figure; hold on;
plot(stim_to_move,earliest_rxn,'.k')
line(xlim,xlim);
xlabel('Rxn');
ylabel('Alt-quiescence rxn');

% movement relative to earliest quiescence period
surround_t_centers = -10:0.1:10;
surround_times = earliest_quiescence_t + surround_t_centers;
surround_move = interp1(t,+wheel_move,surround_times,'previous');

stim_surround_times = stimOn_times + surround_t_centers;
stim_surround_move = interp1(t,+wheel_move,stim_surround_times,'previous');

figure;
subplot(4,2,2:2:6);
imagesc(surround_t_centers,[],stim_surround_move);
colormap(gca,brewermap([],'Greys'));
title('Stim')
subplot(4,2,1:2:6);
imagesc(surround_t_centers,[],surround_move);
colormap(gca,brewermap([],'Greys'));
title('Quiescence')
subplot(4,1,4); hold on;
plot(surround_t_centers,nanmean(stim_surround_move,1),'k');
plot(surround_t_centers,nanmean(surround_move,1),'r');
ylim([0,1]);
legend({'Stim','Quiescence'});

%% (testing above reaction-test in batch)

% tetO mice
animal_group = 'teto';
animals = {'AP100','AP101','AP103','AP104','AP105','AP106'};

% % corticostriatal mice
% animal_group = 'cstr';
% animals = {'AP092','AP093','AP094','AP095','AP096','AP097'};

protocol = 'AP_stimWheelRight';
flexible_name = false;
bhv = struct;

for curr_animal = 1:length(animals)
    
    preload_vars = who;
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    
    if isempty(experiments)
        disp(['No behavior data: ' animal]);
        continue
    end
    
    disp([animal ', day:'])
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment_num = experiments(curr_day).experiment;
        
        % If multiple experiments, only use the last one
        % (usually multiple happens if mess ups and final one is good)
        for curr_experiment = length(experiment_num)
            
            experiment = experiment_num(curr_experiment);
            
            % Load (behavior/timeline only)
            load_parts.imaging = false;
            AP_load_experiment;
            
            %%%%%%%% TESTING: reaction time test
            t = Timeline.rawDAQTimestamps';
            
            % Get trace giving time from last movement
            % (processed wheel movement)
            t_from_move_trace = ...
                t - interp1(t(wheel_move),t(wheel_move),t,'previous','extrap');
    
            % Get trace when out of ITI period
            noniti_flag = logical(interp1( ...
                [0;signals_events.trialNumTimes(1:n_trials)';signals_events.responseTimes(1:n_trials)'], ...
                [0;ones(n_trials,1);zeros(n_trials,1)], ...
                Timeline.rawDAQTimestamps','previous','extrap'));

%             % (if quiescence was specific for each trial - can shorten artifically)
%             % quiescence_shorten = 0;
%             quiescenceMin = min(signals_events.trialQuiescenceValues);
%             quiescence_shorten = arrayfun(@(x) ...
%                 randsample(quiescenceMin:0.1:max(x-0.1,quiescenceMin),1), ...
%                 signals_events.trialQuiescenceValues);
%             trial_quiescence_trace = interp1(signals_events.trialQuiescenceTimes, ...
%                 signals_events.trialQuiescenceValues-quiescence_shorten,t,'previous','extrap');
%             min_quiescence_flag = t_from_move_trace >= trial_quiescence_trace & noniti_flag;
            
            % (shorter quiescence time generate a lot)
            n_null = 10;
            quiescenceMin = min(signals_events.trialQuiescenceValues);
            quiescence_shorten = cell2mat(arrayfun(@(x) ...
                reshape(randsample(quiescenceMin:0.1:max(x-0.1,quiescenceMin),n_null,true),[],1), ...
                signals_events.trialQuiescenceValues,'uni',false));
            
            surround_move_max = nan(n_null,1);
            for curr_null = 1:n_null
                
                trial_quiescence_trace = interp1(signals_events.trialQuiescenceTimes, ...
                    signals_events.trialQuiescenceValues-quiescence_shorten(curr_null,:),t,'previous','extrap');
                min_quiescence_flag = t_from_move_trace >= trial_quiescence_trace & noniti_flag;
                
                % Get first time quiescence period is met
                min_quiescence_onsets = t(diff([0;min_quiescence_flag]) == 1);
                earliest_quiescence_t = arrayfun(@(x) ...
                    min_quiescence_onsets(find(min_quiescence_onsets > x,1,'first')), ...
                    signals_events.trialNumTimes(1:n_trials)');
                
                % Get the first movement from min quiescence would have been
                earliest_quiescence_wheelstarts = arrayfun(@(x) ...
                    wheel_starts(find(wheel_starts > x,1,'first')), ...
                    earliest_quiescence_t);
                
                earliest_rxn = earliest_quiescence_wheelstarts - earliest_quiescence_t;
                
                % movement relative to earliest quiescence period
                surround_t_centers = -10:0.1:10;
                surround_times = earliest_quiescence_t + surround_t_centers;
                surround_move = interp1(t,+wheel_move,surround_times,'previous');
                surround_move_max(curr_null) = max(nanmean(surround_move(:,surround_t_centers>0,1)));
                
                AP_print_progress_fraction(curr_null,n_null);
                
            end
            
            % movement relative to earliest quiescence period
            surround_t_centers = -10:0.1:10;
            stim_surround_times = stimOn_times + surround_t_centers;
            stim_surround_move = interp1(t,+wheel_move,stim_surround_times,'previous');
            stim_surround_move_max = max(nanmean(stim_surround_move(:,surround_t_centers>0,1)));
            
            quiescence_stim_move_ratio = ...
                (stim_surround_move_max - surround_move_max)./ ...
                (stim_surround_move_max + surround_move_max);
                       
            % Store in behavior structure
            bhv(curr_animal).animal = animal;
            bhv(curr_animal).day{curr_day} = day;
            
            bhv(curr_animal).quiescence_stim_move_ratio(curr_day) = quiescence_stim_move_ratio;

            
            AP_print_progress_fraction(curr_day,length(experiments));
        end
        
    end
end


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

% (get max days for padding)
max_days = max(cellfun(@sum,use_days));

% Pad move from quiescence
quiescence_stim_move_ratio_allpad = cell2mat(cellfun(@(x) ...
    padarray(x,[0,max_days-length(x)],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days), ...
    {bhv.quiescence_stim_move_ratio},use_days,'uni',false)','uni',false));

quiescence_stim_move_ratio_allpad_norm = ...
    quiescence_stim_move_ratio_allpad./max(quiescence_stim_move_ratio_allpad,[],2);
figure;plot(quiescence_stim_move_ratio_allpad_norm');

%% (testing above reaction-test in batch: daysplit)

% tetO mice
animal_group = 'teto';
animals = {'AP100','AP101','AP103','AP104','AP105','AP106'};

% % corticostriatal mice
% animal_group = 'cstr';
% animals = {'AP092','AP093','AP094','AP095','AP096','AP097'};

protocol = 'AP_stimWheelRight';
flexible_name = false;
bhv = struct;

for curr_animal = 1:length(animals)
    
    preload_vars = who;
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    
    if isempty(experiments)
        disp(['No behavior data: ' animal]);
        continue
    end
    
    disp([animal ', day:'])
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment_num = experiments(curr_day).experiment;
        
        % If multiple experiments, only use the last one
        % (usually multiple happens if mess ups and final one is good)
        for curr_experiment = length(experiment_num)
            
            experiment = experiment_num(curr_experiment);
            
            % Load (behavior/timeline only)
            load_parts.imaging = false;
            AP_load_experiment;
            
            %%%%%%%% TESTING: reaction time test
            t = Timeline.rawDAQTimestamps';
            
            % Get trace giving time from last movement
            % (processed wheel movement)
            t_from_move_trace = ...
                t - interp1(t(wheel_move),t(wheel_move),t,'previous','extrap');
            
            % Get trace when out of ITI period
            noniti_flag = logical(interp1( ...
                [0;signals_events.trialNumTimes(1:n_trials)';signals_events.responseTimes(1:n_trials)'], ...
                [0;ones(n_trials,1);zeros(n_trials,1)], ...
                Timeline.rawDAQTimestamps','previous','extrap'));
            
            %             % (if quiescence was specific for each trial - can shorten artifically)
            %             % quiescence_shorten = 0;
            %             quiescenceMin = min(signals_events.trialQuiescenceValues);
            %             quiescence_shorten = arrayfun(@(x) ...
            %                 randsample(quiescenceMin:0.1:max(x-0.1,quiescenceMin),1), ...
            %                 signals_events.trialQuiescenceValues);
            %             trial_quiescence_trace = interp1(signals_events.trialQuiescenceTimes, ...
            %                 signals_events.trialQuiescenceValues-quiescence_shorten,t,'previous','extrap');
            %             min_quiescence_flag = t_from_move_trace >= trial_quiescence_trace & noniti_flag;
            
            n_daysplit = 4;
            day_split_idx = min(floor(linspace(1,n_daysplit+1,n_trials)),n_daysplit);
            
            % (shorter quiescence time generate a lot)
            n_null = 10;
            quiescenceMin = min(signals_events.trialQuiescenceValues);
            quiescence_shorten = cell2mat(arrayfun(@(x) ...
                reshape(randsample(quiescenceMin:0.1:max(x-0.1,quiescenceMin),n_null,true),[],1), ...
                signals_events.trialQuiescenceValues,'uni',false));
            
            surround_move_max = nan(n_null,n_daysplit);
            for curr_null = 1:n_null
                
                trial_quiescence_trace = interp1(signals_events.trialQuiescenceTimes, ...
                    signals_events.trialQuiescenceValues-quiescence_shorten(curr_null,:),t,'previous','extrap');
                min_quiescence_flag = t_from_move_trace >= trial_quiescence_trace & noniti_flag;
                
                % Get first time quiescence period is met
                min_quiescence_onsets = t(diff([0;min_quiescence_flag]) == 1);
                earliest_quiescence_t = arrayfun(@(x) ...
                    min_quiescence_onsets(find(min_quiescence_onsets > x,1,'first')), ...
                    signals_events.trialNumTimes(1:n_trials)');
                
                % Get the first movement from min quiescence would have been
                earliest_quiescence_wheelstarts = arrayfun(@(x) ...
                    wheel_starts(find(wheel_starts > x,1,'first')), ...
                    earliest_quiescence_t);
                
                earliest_rxn = earliest_quiescence_wheelstarts - earliest_quiescence_t;
                
                % movement relative to earliest quiescence period
                surround_t_centers = -10:0.1:10;
                surround_times = earliest_quiescence_t + surround_t_centers;
                surround_move = interp1(t,+wheel_move,surround_times,'previous');
                
                for curr_split = 1:n_daysplit
                    surround_move_max(curr_null,curr_split) = ...
                        max(nanmean(surround_move(day_split_idx == curr_split,surround_t_centers>0,1)));
                end
                
%                 AP_print_progress_fraction(curr_null,n_null);
                
            end
            
            % movement relative to earliest quiescence period
            surround_t_centers = -10:0.1:10;
            stim_surround_times = stimOn_times + surround_t_centers;
            stim_surround_move = interp1(t,+wheel_move,stim_surround_times,'previous');
            stim_surround_move_max = max(nanmean(stim_surround_move(:,surround_t_centers>0,1)));
            
            stim_surround_move_max = nan(1,curr_split);
            for curr_split = 1:n_daysplit
                stim_surround_move_max(curr_split) = ...
                    max(nanmean(stim_surround_move(day_split_idx == curr_split,surround_t_centers>0,1)));
            end
            
            quiescence_stim_move_ratio = ...
                (stim_surround_move_max - surround_move_max)./ ...
                (stim_surround_move_max + surround_move_max);
            
            % Store in behavior structure
            bhv(curr_animal).animal = animal;
            bhv(curr_animal).day{curr_day} = day;
            
            bhv(curr_animal).stim_surround_move_max{curr_day} = stim_surround_move_max;
            bhv(curr_animal).surround_move_max{curr_day} = nanmedian(surround_move_max,1);
            bhv(curr_animal).quiescence_stim_move_ratio{curr_day} = nanmedian(quiescence_stim_move_ratio,1);
            
            AP_print_progress_fraction(curr_day,length(experiments));
            
        end
    end
end


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

% (get max days for padding)
max_days = max(cellfun(@sum,use_days));

% Pad values
stim_surround_move_max_cat = cell2mat(permute(cellfun(@(x,use_days) ...
    padarray(vertcat(x{use_days}),[max_days-sum(use_days),0],NaN,'post'), ...
    {bhv.stim_surround_move_max},use_days,'uni',false),[1,3,2]));
surround_move_max_cat = cell2mat(permute(cellfun(@(x,use_days) ...
    padarray(vertcat(x{use_days}),[max_days-sum(use_days),0],NaN,'post'), ...
    {bhv.surround_move_max},use_days,'uni',false),[1,3,2]));

stim_surround_move_max_mean = nanmean(stim_surround_move_max_cat,3);
stim_surround_move_max_sem = AP_sem(stim_surround_move_max_cat,3);
surround_move_max_mean = nanmean(surround_move_max_cat,3);
surround_move_max_sem = AP_sem(surround_move_max_cat,3);
figure; hold on;
errorbar( ...
    reshape(padarray(stim_surround_move_max_mean,[0,1],NaN,'post')',[],1), ...
    reshape(padarray(stim_surround_move_max_sem,[0,1],NaN,'post')',[],1),'k','linewidth',2);
errorbar( ...
    reshape(padarray(surround_move_max_mean,[0,1],NaN,'post')',[],1), ...
    reshape(padarray(surround_move_max_sem,[0,1],NaN,'post')',[],1),'b','linewidth',2);
ylabel(' probability');
xlabel('Day');
legend({'Stim move','Shuffle move'});

move_ratio_cat = cell2mat(permute(cellfun(@(x,use_days) ...
    padarray(vertcat(x{use_days}),[max_days-sum(use_days),0],NaN,'post'), ...
    {bhv.quiescence_stim_move_ratio},use_days,'uni',false),[1,3,2]));

move_ratio_mean = nanmean(move_ratio_cat,3);
move_ratio_sem = AP_sem(move_ratio_cat,3);
figure;errorbar( ...
    reshape(padarray(move_ratio_mean,[0,1],NaN,'post')',[],1), ...
    reshape(padarray(move_ratio_sem,[0,1],NaN,'post')',[],1),'k','linewidth',2);
ylabel('Stim/null move probability');
xlabel('Day');

%% Compare stim-aligned movement to "baseline" movement relative to outcome
% (Kenneth suggested this)

t = Timeline.rawDAQTimestamps';

% Non-stim move trace
wheel_move_nostim = +wheel_move;
wheel_move_nostim(stimOn_epochs) = NaN;

% Get day split index
n_daysplit = 4;
day_split_idx = min(floor(linspace(1,n_daysplit+1,n_trials)),n_daysplit);

figure; hold on;
for curr_daysplit = 1:n_daysplit
    
    use_trials = day_split_idx == curr_daysplit;
    
    % Get outcome-to-stim wheel move
    response_times = signals_events.responseTimes(use_trials)';
    
    iti_surround_t_centers = 1.1:0.1:30;
    iti_surround_times = response_times + iti_surround_t_centers;
    % (align move to ITI)
    iti_surround_move = interp1(t,wheel_move_nostim,iti_surround_times,'previous');
    % (align stim to ITI, index times after stim onset)
    iti_surround_stim = cumsum(interp1(t,+stimOn_epochs,iti_surround_times,'previous'),2) > 0;
    % (get ITI-surround move pre-stim)
    iti_surround_move_prestim = iti_surround_move;
    iti_surround_move_prestim(iti_surround_stim) = NaN;
    % (average to get probability of movement given no stim)
    outcome_to_stim_move = nanmean(iti_surround_move_prestim);
    
    % Stim-align movement and remove "baseline"
    % (interpolate time from outcome)
    t_from_outcome = t - ...
        interp1(response_times,response_times,t,'previous','extrap');
    % (discretize into baseline movement bins and index baseline movement)
    t_from_outcome_idx = discretize(t_from_outcome,[iti_surround_t_centers,Inf]);
    baseline_move = nan(size(wheel_move));
    baseline_move(~isnan(t_from_outcome_idx)) = outcome_to_stim_move(t_from_outcome_idx(~isnan(t_from_outcome_idx)));
    % (stim-align movement and baseline)
    stim_surround_t_centers = -10:0.1:10;
    stim_surround_times = stimOn_times(use_trials) + stim_surround_t_centers;
    stim_surround_move = interp1(t,+wheel_move,stim_surround_times,'previous');
    stim_surround_move_baseline = interp1(t,+baseline_move,stim_surround_times,'previous');
    % (get max post-stim movement-baseline)
    stim_surround_move_overbaseline = nanmean(stim_surround_move-stim_surround_move_baseline,1);
    stim_response_value = max(stim_surround_move_overbaseline(stim_surround_t_centers>0));
    
    subplot(n_daysplit,1,curr_daysplit); hold on;
    plot(stim_surround_t_centers,nanmean(stim_surround_move))
    plot(stim_surround_t_centers,nanmean(stim_surround_move_baseline))
    
end



%% (above in batch daysplit)

% tetO mice
animal_group = 'teto';
animals = {'AP100','AP101','AP103','AP104','AP105','AP106'};

% % corticostriatal mice
% animal_group = 'cstr';
% animals = {'AP092','AP093','AP094','AP095','AP096','AP097'};

protocol = 'AP_stimWheelRight';
flexible_name = false;
bhv = struct;

for curr_animal = 1:length(animals)
    
    preload_vars = who;
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    
    if isempty(experiments)
        disp(['No behavior data: ' animal]);
        continue
    end
    
    disp([animal ', day:'])
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment_num = experiments(curr_day).experiment;
        
        % If multiple experiments, only use the last one
        % (usually multiple happens if mess ups and final one is good)
        for curr_experiment = length(experiment_num)
            
            experiment = experiment_num(curr_experiment);
            
            % Load (behavior/timeline only)
            load_parts.imaging = false;
            AP_load_experiment;
            
            
            t = Timeline.rawDAQTimestamps';
            
            % Non-stim move trace
            wheel_move_nostim = +wheel_move;
            wheel_move_nostim(stimOn_epochs) = NaN;
            
            % Get day split index
            n_daysplit = 4;
            day_split_idx = min(floor(linspace(1,n_daysplit+1,n_trials)),n_daysplit);
            
            stim_response_value = nan(1,n_daysplit);
            for curr_daysplit = 1:n_daysplit
                
                use_trials = day_split_idx == curr_daysplit;
                
                % Get outcome-to-stim wheel move
                response_times = signals_events.responseTimes(use_trials)';
                
                iti_surround_t_centers = 1.1:0.1:30;
                iti_surround_times = response_times + iti_surround_t_centers;
                % (align move to ITI)
                iti_surround_move = interp1(t,wheel_move_nostim,iti_surround_times,'previous');
                % (align stim to ITI, index times after stim onset)
                iti_surround_stim = cumsum(interp1(t,+stimOn_epochs,iti_surround_times,'previous'),2) > 0;
                % (get ITI-surround move pre-stim)
                iti_surround_move_prestim = iti_surround_move;
                iti_surround_move_prestim(iti_surround_stim) = NaN;
                % (average to get probability of movement given no stim)
                outcome_to_stim_move = nanmean(iti_surround_move_prestim);
                
                % Stim-align movement and remove "baseline"
                % (interpolate time from outcome)
                t_from_outcome = t - ...
                    interp1(response_times,response_times,t,'previous','extrap');
                % (discretize into baseline movement bins and index baseline movement)
                t_from_outcome_idx = discretize(t_from_outcome,[iti_surround_t_centers,Inf]);
                baseline_move = nan(size(wheel_move));
                baseline_move(~isnan(t_from_outcome_idx)) = outcome_to_stim_move(t_from_outcome_idx(~isnan(t_from_outcome_idx)));
                % (stim-align movement and baseline)
                stim_surround_t_centers = -10:0.1:10;
                stim_surround_times = stimOn_times(use_trials) + stim_surround_t_centers;
                stim_surround_move = interp1(t,+wheel_move,stim_surround_times,'previous');
                stim_surround_move_baseline = interp1(t,+baseline_move,stim_surround_times,'previous');
                % (get max post-stim movement-baseline)
                stim_surround_move_overbaseline = nanmean(stim_surround_move-stim_surround_move_baseline,1);
                stim_response_value(curr_daysplit) = max(stim_surround_move_overbaseline(stim_surround_t_centers>0));               
                
            end
           
            % Store in behavior structure
            bhv(curr_animal).animal = animal;
            bhv(curr_animal).day{curr_day} = day;
            
            bhv(curr_animal).stim_response_value{curr_day,1} = stim_response_value;
            
            AP_print_progress_fraction(curr_day,length(experiments));
            
        end
    end
end


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

% (get max days for padding)
max_days = max(cellfun(@sum,use_days));

% Pad stim response
stim_response_value_cat = cell2mat(permute(cellfun(@(x,use_days) ...
    padarray(vertcat(x{use_days}),[max_days-sum(use_days),0],NaN,'post'), ...
    {bhv.stim_response_value},use_days,'uni',false),[1,3,2]));

stim_response_value_cat_mean = nanmean(stim_response_value_cat,3);
stim_response_value_cat_sem = AP_sem(stim_response_value_cat,3);
figure;errorbar( ...
    reshape(padarray(stim_response_value_cat_mean,[0,1],NaN,'post')',[],1), ...
    reshape(padarray(stim_response_value_cat_sem,[0,1],NaN,'post')',[],1),'k','linewidth',2);
ylabel('Stim move - baseline move prob.');
xlabel('Day');


%% Stim response vs alternate timing stim onset
% (Kenneth suggested this: like the quiescence-shortening, but using times
% from other trials and both ITI and quiescence)
% (LOOKS OK?)
%
% For each trial: 
% - find another trial params where stim would've appeared earlier
% - align movement to the "alternate stim on"
% - don't use if no alternate earlier trials or not outside analysis window
% 
% compare: reaction time histogram (shouldn't have peak but should have
% same pedestal), triggered-average wheel speed/binary

t = Timeline.rawDAQTimestamps';

% Get trace giving time from last movement 
% (processed wheel movement)
t_from_move_trace = ...
    t - interp1(t(wheel_move),t(wheel_move),t,'previous','extrap');

% (quiescence watch full re-creation, takes too much time)
% % (trying to re-create quiescence watch in signals)
% % (NOTE: signals doesn't use clicks, it uses cumulative 1mm)
% % (this is dirty but I don't have a better solution at the moment)
% t_wheel_block = interp1(block2timeline,timeline2block,block.inputs.wheelMMTimes);
% a = block.inputs.wheelMMValues;
% quiescence_reset = false(size(block.inputs.wheelMMValues));
% i = 1;
% while i < length(block.inputs.wheelMMValues)
%     a(i:end) = a(i:end)-a(i);
%     curr_diff = abs(block.inputs.wheelMMValues-block.inputs.wheelMMValues(i));
%     thresh = 0.95;
%     i = find(curr_diff(i:end) > thresh,1,'first') + (i-1);
%     quiescence_reset(i) = true;
%     AP_print_progress_fraction(i,length(block.inputs.wheelMMValues));
% end
% % (quiescence watch starts on trial onset)
% quiescence_reset_t = sort([t_wheel_block(quiescence_reset & ~isnan(t_wheel_block)),signals_events.trialNumTimes]);
% t_from_quiescence_reset = t - ...
%     interp1(quiescence_reset_t,quiescence_reset_t,t,'previous','extrap');

alt_stimOn_times = cell(n_trials,1);
real_stimOn_times = cell(n_trials,1);
% (skip first since no ITI)
for curr_trial = 2:n_trials
    
    curr_trial_t_idx = t >= signals_events.responseTimes(curr_trial-1) & ...
        t <= stimOn_times(curr_trial);
    
    curr_trial_t = t(curr_trial_t_idx);
    
    % Re-create quiescence watch: resets at >1mm cumulative
    % (this basically works, but not 100%: sometimes I see two clicks that
    % should've reset it but it doesn't catch it, maybe because of exact
    % processing times in block or something)
    t_wheel_block = interp1(block2timeline,timeline2block,block.inputs.wheelMMTimes,'linear','extrap');
    % (quiescence watch starts on new trial)
    curr_trial_t_block_idx = t_wheel_block >= signals_events.responseTimes(curr_trial-1) & ...
        t_wheel_block <= stimOn_times(curr_trial);
    % (use wheelMM in block, loop through to estimate resets)
    curr_wheel_mm_t = t_wheel_block(curr_trial_t_block_idx);
    curr_wheel_mm = block.inputs.wheelMMValues(curr_trial_t_block_idx);  
    curr_quiescence_reset = false(size(curr_wheel_mm));
    i = 1;
    while i < length(curr_wheel_mm)
        curr_diff = abs(curr_wheel_mm - curr_wheel_mm(i));
        thresh = 0.95; % (it's really 1, but sometimes finnicky?)
        i = find(curr_diff(i:end) > thresh,1,'first') + (i-1);
        curr_quiescence_reset(i) = true;
    end
    quiescence_reset_t = curr_wheel_mm_t(curr_quiescence_reset)';
    t_from_quiescence_reset_full = t - ...
        interp1(quiescence_reset_t,quiescence_reset_t,t,'previous','extrap');
    t_from_quiescence_reset_trial = t_from_quiescence_reset_full(curr_trial_t_idx);
    
%     % (sanity plot)
%     t_plot_scale = 0.1;
%     figure; hold on;
%     plot(t(curr_trial_t_idx),wheel_velocity(curr_trial_t_idx),'k')
%     plot(t(curr_trial_t_idx),t_from_move_trace(curr_trial_t_idx)*t_plot_scale,'r');
%     plot(t(curr_trial_t_idx),t_from_quiescence_reset_trial*t_plot_scale,'b');
%     plot(t(curr_trial_t_idx),[0;diff(wheel_position(curr_trial_t_idx))]*0.1,'g')    
%     line(repmat(curr_trial_t(1)+signals_events.trialITIValues(curr_trial-1),2,1),ylim);
%     line(xlim,repmat(signals_events.trialQuiescenceValues(curr_trial),2,1)*t_plot_scale,'color','m');
    
    % Find alternate stim times from all other trial timings
    other_trial_idx = setdiff(1:n_trials,curr_trial);
    alt_stim_timing_grid = ...
        ((t(curr_trial_t_idx) - curr_trial_t(1)) > signals_events.trialITIValues(other_trial_idx)) & ...
        (t_from_quiescence_reset_trial > signals_events.trialQuiescenceValues(other_trial_idx));
    [alt_stim_value,alt_stim_idx] = max(alt_stim_timing_grid,[],1);
        
    % Store possible stim on times (where timings were hit, and not within
    % an analysis window around the real stim on time)
    leeway_window = 0.5; % seconds from real stim onset to exclude
    leeway_idx = find(curr_trial_t-curr_trial_t(end) < -leeway_window,1,'last');
    alt_stimOn_times{curr_trial} = curr_trial_t( ...
        alt_stim_idx(alt_stim_value & alt_stim_idx <= leeway_idx));
    
    % (sanity check: store real stim on times for each trial)
    real_stimOn_times{curr_trial} = repmat(stimOn_times(curr_trial), ...
        length(alt_stimOn_times{curr_trial}),1);
    
%     AP_print_progress_fraction(curr_trial,n_trials)
    
end

% Align activity to stim/all alt-stim onsets
surround_t_centers = -10:0.01:10;

stim_surround_times = stimOn_times + surround_t_centers;
stim_surround_move = interp1(t,+wheel_move,stim_surround_times,'previous');

alt_stim_surround_times = cell2mat(alt_stimOn_times) + surround_t_centers;
alt_stim_surround_move = interp1(t,+wheel_move,alt_stim_surround_times,'previous');

% Approximate reaction times
[~,stim_surround_move_idx] = max(stim_surround_move(:,surround_t_centers > 0),[],2);
stim_rxn = surround_t_centers(find(surround_t_centers > 0,1) + stim_surround_move_idx);

[~,alt_stim_surround_move_idx] = max(alt_stim_surround_move(:,surround_t_centers > 0),[],2);
alt_stim_rxn = surround_t_centers(find(surround_t_centers > 0,1) + alt_stim_surround_move_idx);

rxn_bins = [-Inf,0:0.01:1,Inf];
figure; hold on;
histogram(alt_stim_rxn,rxn_bins,'EdgeColor','none','normalization','pdf');
histogram(stim_rxn,rxn_bins,'EdgeColor','none','normalization','pdf');
xlabel('Reaction time');
ylabel('Prob');
legend({'Null','Measured'});

% Plot stim-aligned move
figure;
subplot(1,3,1);
imagesc(surround_t_centers,[],stim_surround_move);
subplot(1,3,2);
imagesc(surround_t_centers,[],alt_stim_surround_move);
% (plot real stim on times for each alternate)
hold on
plot(cell2mat(real_stimOn_times)-cell2mat(alt_stimOn_times),1:length(cell2mat(real_stimOn_times)),'.r');
subplot(1,3,3); hold on;
plot(surround_t_centers,nanmean(stim_surround_move,1),'linewidth',2);
plot(surround_t_centers,nanmean(alt_stim_surround_move,1),'linewidth',2);
colormap(brewermap([],'Greys'));


% Compare same stim/alt-stim trials (sub-sample alt-stim)
n_daysplit = 4;
day_split_idx = min(floor(linspace(1,n_daysplit+1,n_trials)),n_daysplit)';

use_alt_trials = cellfun(@(x) ~isempty(x),alt_stimOn_times);

stim_move_max = nan(n_daysplit,1);
alt_stim_move_max = nan(n_daysplit,1);
for curr_daysplit = 1:n_daysplit   
    
    curr_trials = day_split_idx == curr_daysplit & use_alt_trials;
    
    surround_t_centers = 0:0.01:2;
    curr_stim_surround_times = stimOn_times(curr_trials) + surround_t_centers;
    curr_stim_surround_move = interp1(t,+wheel_move,curr_stim_surround_times,'previous');
    
    % Get random subset alt trials to match 1:1 to real trials
    n_alt = 1000;
    curr_alt_subset = cell2mat(cellfun(@(x) randsample(x,n_alt,true)', ...
        alt_stimOn_times(curr_trials),'uni',false));    
      
    curr_alt_stim_surround_times = permute(curr_alt_subset,[1,3,2]) + surround_t_centers;
    curr_alt_stim_surround_move = interp1(t,+wheel_move,curr_alt_stim_surround_times,'previous');    
    
    stim_move_max(curr_daysplit) = max(nanmean(curr_stim_surround_move,1));
    alt_stim_move_max(curr_daysplit) = nanmean(max(nanmean(curr_alt_stim_surround_move,1),[],2));
    
end
stim_move_response_idx = (stim_move_max-alt_stim_move_max)./(stim_move_max+alt_stim_move_max);



%% (above in batch)
%
% For each trial: 
% - find another trial params where stim would've appeared earlier
% - align movement to the "alternate stim on"
% - don't use if no alternate earlier trials or not outside analysis window
% 
% compare: reaction time histogram (shouldn't have peak but should have
% same pedestal), triggered-average wheel speed/binary


% tetO mice
animal_group = 'teto';
% animals = {'AP100','AP101','AP103','AP104','AP105','AP106'};
animals = {'AP107','AP108','AP109'};

% % corticostriatal mice
% animal_group = 'cstr';
% animals = {'AP092','AP093','AP094','AP095','AP096','AP097'};

protocol = 'AP_stimWheelRight';
flexible_name = false;
bhv = struct;

for curr_animal = 1:length(animals)
    
    preload_vars = who;
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    
    if isempty(experiments)
        disp(['No behavior data: ' animal]);
        continue
    end
    
    disp([animal ', day:'])
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment_num = experiments(curr_day).experiment;
        
        % If multiple experiments, only use the last one
        % (usually multiple happens if mess ups and final one is good)
        for curr_experiment = length(experiment_num)
            
            experiment = experiment_num(curr_experiment);
            
            % Load (behavior/timeline only)
            load_parts.imaging = false;
            AP_load_experiment;
            
            
            t = Timeline.rawDAQTimestamps';
            
            % Get trace giving time from last movement
            % (processed wheel movement)
            t_from_move_trace = ...
                t - interp1(t(wheel_move),t(wheel_move),t,'previous','extrap');
            
            alt_stimOn_times = cell(n_trials,1);
            % (skip first since no ITI)
            for curr_trial = 2:n_trials
                
                curr_trial_t_idx = t >= signals_events.responseTimes(curr_trial-1) & ...
                    t <= stimOn_times(curr_trial);
                
                curr_trial_t = t(curr_trial_t_idx);
                
                % Re-create quiescence watch: resets at >1mm cumulative
                % (this basically works, but not 100%: sometimes I see two clicks that
                % should've reset it but it doesn't catch it, maybe because of exact
                % processing times in block or something)
                t_wheel_block = interp1(block2timeline,timeline2block,block.inputs.wheelMMTimes,'linear','extrap');
                % (quiescence watch starts on new trial)
                curr_trial_t_block_idx = t_wheel_block >= signals_events.responseTimes(curr_trial-1) & ...
                    t_wheel_block <= stimOn_times(curr_trial);
                % (use wheelMM in block, loop through to estimate resets)
                curr_wheel_mm_t = t_wheel_block(curr_trial_t_block_idx);
                curr_wheel_mm = block.inputs.wheelMMValues(curr_trial_t_block_idx);
                curr_quiescence_reset = false(size(curr_wheel_mm));
                i = 1;
                while i < length(curr_wheel_mm)
                    curr_diff = abs(curr_wheel_mm - curr_wheel_mm(i));
                    thresh = 0.95; % (it's really 1, but sometimes finnicky?)
                    i = find(curr_diff(i:end) > thresh,1,'first') + (i-1);
                    curr_quiescence_reset(i) = true;
                end
                % (skip trial if < 2 quiescence resets - can't interpolate)
                if sum(curr_quiescence_reset) < 2
                    continue
                end
                quiescence_reset_t = curr_wheel_mm_t(curr_quiescence_reset)';
                t_from_quiescence_reset_full = t - ...
                    interp1(quiescence_reset_t,quiescence_reset_t,t,'previous','extrap');
                t_from_quiescence_reset_trial = t_from_quiescence_reset_full(curr_trial_t_idx);
                
                % Find alternate stim times from all other trial timings
                other_trial_idx = setdiff(1:n_trials,curr_trial);
                alt_stim_timing_grid = ...
                    ((t(curr_trial_t_idx) - curr_trial_t(1)) > signals_events.trialITIValues(other_trial_idx)) & ...
                    (t_from_quiescence_reset_trial > signals_events.trialQuiescenceValues(other_trial_idx));
                [alt_stim_value,alt_stim_idx] = max(alt_stim_timing_grid,[],1);
                
                % Store possible stim on times (where timings were hit, and not within
                % an analysis window around the real stim on time)
                leeway_window = 0.5; % seconds from real stim onset to exclude
                leeway_idx = find(curr_trial_t-curr_trial_t(end) < -leeway_window,1,'last');
                alt_stimOn_times{curr_trial} = curr_trial_t( ...
                    alt_stim_idx(alt_stim_value & alt_stim_idx <= leeway_idx));
                
            end
            
            % Get would-be reaction time after alt stim times
            stim_leeway = 0.1;
            wheel_move_alt_stim_idx = ...
                arrayfun(@(stim) find(wheel_starts > stim-stim_leeway,1,'first'), ...
                cell2mat(alt_stimOn_times));

            alt_stim_to_move = ...
                mat2cell(wheel_starts(wheel_move_alt_stim_idx) - cell2mat(alt_stimOn_times), ...
                cellfun(@length,alt_stimOn_times));
            
            % Compare same stim/alt-stim trials (sub-sample alt-stim)
            n_daysplit = 4;
            day_split_idx = min(floor(linspace(1,n_daysplit+1,n_trials)),n_daysplit)';
            
            use_alt_trials = cellfun(@(x) ~isempty(x),alt_stimOn_times);
            
            stim_move_max = nan(1,n_daysplit);
            alt_stim_move_max = nan(1,n_daysplit);
            for curr_daysplit = 1:n_daysplit
                
                curr_trials = day_split_idx == curr_daysplit & use_alt_trials;
                
                surround_t_centers = 0:0.01:2;
                curr_stim_surround_times = stimOn_times(curr_trials) + surround_t_centers;
                curr_stim_surround_move = interp1(t,+wheel_move,curr_stim_surround_times,'previous');
                
                % Get random subset alt trials to match 1:1 to real trials
                n_alt = 1000;
                curr_alt_subset = cell2mat(cellfun(@(x) randsample(x,n_alt,true)', ...
                    alt_stimOn_times(curr_trials),'uni',false));
                
                curr_alt_stim_surround_times = permute(curr_alt_subset,[1,3,2]) + surround_t_centers;
                curr_alt_stim_surround_move = interp1(t,+wheel_move,curr_alt_stim_surround_times,'previous');
                
                stim_move_max(curr_daysplit) = max(nanmean(curr_stim_surround_move,1));
                alt_stim_move_max(curr_daysplit) = nanmean(max(nanmean(curr_alt_stim_surround_move,1),[],2));
                
            end
            
            % Get "stim response index" (change to move after stim vs null)
            stim_response_idx = (stim_move_max-alt_stim_move_max)./(stim_move_max+alt_stim_move_max);
            
            
            % Store in behavior structure
            bhv(curr_animal).animal = animal;
            bhv(curr_animal).day{curr_day} = day;
            % (reaction times)
            bhv(curr_animal).stim_move_t{curr_day} = stim_to_move;
            bhv(curr_animal).alt_stim_move_t{curr_day} = alt_stim_to_move;
            % (stim response index)
            bhv(curr_animal).stim_move_max{curr_day,1} = stim_move_max;
            bhv(curr_animal).alt_stim_move_max{curr_day,1} = alt_stim_move_max;
            bhv(curr_animal).stim_response_idx{curr_day,1} = stim_response_idx;
            
            AP_print_progress_fraction(curr_day,length(experiments));
            
        end
    end
end


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

% (get max days for padding)
max_days = max(cellfun(@sum,use_days));

% Pad stim response
stim_move_max_plotcat = reshape(permute(padarray(cell2mat(permute(cellfun(@(x,use_days) ...
    padarray(vertcat(x{use_days}),[max_days-sum(use_days),0],NaN,'post'), ...
    {bhv.stim_move_max},use_days,'uni',false),[1,3,2])),[0,1],NaN,'post'),[2,1,3]),[],length(bhv));

alt_stim_move_max_plotcat = reshape(permute(padarray(cell2mat(permute(cellfun(@(x,use_days) ...
    padarray(vertcat(x{use_days}),[max_days-sum(use_days),0],NaN,'post'), ...
    {bhv.alt_stim_move_max},use_days,'uni',false),[1,3,2])),[0,1],NaN,'post'),[2,1,3]),[],length(bhv));

stim_response_idx_plotcat = reshape(permute(padarray(cell2mat(permute(cellfun(@(x,use_days) ...
    padarray(vertcat(x{use_days}),[max_days-sum(use_days),0],NaN,'post'), ...
    {bhv.stim_response_idx},use_days,'uni',false),[1,3,2])),[0,1],NaN,'post'),[2,1,3]),[],length(bhv));

figure; hold on;
errorbar(nanmean(stim_move_max_plotcat,2), ...
    AP_sem(stim_move_max_plotcat,2),'r','linewidth',2);
errorbar(nanmean(alt_stim_move_max_plotcat,2), ...
    AP_sem(alt_stim_move_max_plotcat,2),'b','linewidth',2);
errorbar(nanmean(stim_response_idx_plotcat,2), ...
    AP_sem(stim_response_idx_plotcat,2),'k','linewidth',2);
ylabel('Stim response index');
xlabel('Day');
legend({'Stim move max','Null move max','Stim response idx'});


%%% Reaction time distribution
% (already in other script, but could do something here with null distr)

curr_animal = 1;
curr_day = 8;

rxn_bins = 0:0.01:1;

curr_rxn = bhv(curr_animal).stim_move_t{curr_day};
curr_alt_rxn = bhv(curr_animal).alt_stim_move_t{curr_day};

use_alt_trials = cellfun(@(x) ~isempty(x),curr_alt_rxn);

n_alt = 1000;
curr_alt_rxn_subset = cell2mat(cellfun(@(x) randsample(x,n_alt,true)', ...
    curr_alt_rxn(use_alt_trials),'uni',false));

figure; hold on
histogram(curr_rxn,rxn_bins,'EdgeColor','none','normalization','pdf');
histogram(curr_alt_rxn_subset,rxn_bins,'EdgeColor','none','normalization','pdf');

% (not using atm)
rxn_window = [0.1,0.25];
a = mean(curr_rxn >= rxn_window(1) & curr_rxn <= rxn_window(2));
b = mean(curr_alt_rxn_subset >= rxn_window(1) & curr_alt_rxn_subset <= rxn_window(2),1);

%% ~~~~~~~~~ DRAW WIDEFIELD ROIs  ~~~~~~~~~
% NOTE: the reference was just average all passive

% Set ROIs to draw
roi_areas = {'V1p','V1c','AM','RSPa','RSPp','M2p','FRm','FRa','SMl'};

% Load reference image
wf_roi_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\wf_processing\wf_rois';
wf_ref_fn = 'wf_roi_ref.fig';
openfig([wf_roi_path filesep wf_ref_fn]);

wf_roi = struct('area',cell(length(roi_areas),2),'mask',cell(length(roi_areas),2));
for curr_area = 1:length(roi_areas)
    
    % Get ROI from left hemisphere
    title(['Draw ' roi_areas{curr_area} '_L']);
    curr_mask_L = roipoly;
    wf_roi(curr_area,1).area = [roi_areas{curr_area} '_L'];
    wf_roi(curr_area,1).mask = curr_mask_L;
    
    % Reflect ROI to right hemisphere
    curr_mask_R = AP_reflect_widefield(curr_mask_L) > 0;
    wf_roi(curr_area,2).area = [roi_areas{curr_area} '_R'];
    wf_roi(curr_area,2).mask = curr_mask_R;
    
    % Draw ROIs
    curr_roi_L = cell2mat(bwboundaries(curr_mask_L));
    curr_roi_R = cell2mat(bwboundaries(curr_mask_R));
    plot(curr_roi_L(:,2),curr_roi_L(:,1),'m','linewidth',2);
    plot(curr_roi_R(:,2),curr_roi_R(:,1),'m','linewidth',2);
    drawnow;
    
end

wf_roi_fn = [wf_roi_path filesep 'wf_roi'];
save(wf_roi_fn,'wf_roi');
disp('Saved new widefield ROIs');
    
%% Plot widefield ROIs

wf_roi_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\wf_processing\wf_rois\wf_roi';
load(wf_roi_fn);
n_rois = numel(wf_roi);

roi_cat = cat(3,wf_roi.mask);
roi_col = [autumn(size(wf_roi,1));winter(size(wf_roi,1))];

figure; hold on
set(gca,'YDir','reverse');
AP_reference_outline('ccf_aligned','k');%AP_reference_outline('retinotopy','m');
for curr_roi = 1:n_rois
    curr_roi_boundary = cell2mat(bwboundaries(roi_cat(:,:,curr_roi)));
    patch(curr_roi_boundary(:,2),curr_roi_boundary(:,1),roi_col(curr_roi,:));
    
    text(nanmean(curr_roi_boundary(:,2)),nanmean(curr_roi_boundary(:,1)), ...
        wf_roi(curr_roi).area,'FontSize',12,'HorizontalAlignment','center')
end
axis image off;



%% ~~~~~~~~~ PIP'S NAIVE EXPERIMENTS  ~~~~~~~~~

% PC052/3/4/5

% ephys path is something different
ephys_path = '\\128.40.224.65\Subjects\PC052\2021-09-13\ephys\kilosort\imec0';

% params has different name
header_path = [ephys_path filesep 'params.py'];

% sync file is a little different
sync_long = load('\\128.40.224.65\Subjects\PC052\2021-09-13\ephys\kilosort\imec0\sync.mat');
sync(4).timestamps = (find(diff(sync_long.sync) ~= 0)+1)'./ephys_sample_rate;

% (also looks like good/bad sort not done yet)

animal = 'PC052';day = '2021-09-13';experiment = 4;verbose = true;AP_load_experiment;




%% ~~~~~~~~~ GRAB & SAVE BATCH  ~~~~~~~~~

%% Passive - corticostriatal

clear all
disp('Passive trial activity (corticostriatal)')

% (AP089/90/91: not variable quiescence)
animals = {'AP092','AP093','AP094','AP095','AP096','AP097'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);
    
    % Get days with muscimol
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    muscimol_start_day = muscimol(muscimol_animal_idx).day{1};
    muscimol_experiments = datenum({experiments.day})' >= datenum(muscimol_start_day);
    
    % Set experiments to use (imaging, not muscimol)
    experiments = experiments([experiments.imaging] & ~muscimol_experiments);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        AP_load_experiment;
        
        % Pull out trial data
        operant_grab_trial_data;
        
        % Store trial data into master structure
        trial_data_fieldnames = fieldnames(trial_data);
        for curr_trial_data_field = trial_data_fieldnames'
            trial_data_all.(cell2mat(curr_trial_data_field)){curr_animal,1}{curr_day,1} = ...
                trial_data.(cell2mat(curr_trial_data_field));
        end
        
        % Store general info
        trial_data_all.animals = animals;
        trial_data_all.t = t;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
    
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = ['trial_activity_passive_cstr'];
save([save_path filesep save_fn],'-v7.3');
disp(['Saved: ' save_path filesep save_fn])

%% Task - corticostriatal

clear all
disp('Task trial activity (corticostriatal)')

% (AP089/90/91: not variable quiescence)
animals = {'AP092','AP093','AP094','AP095','AP096','AP097'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_stimWheelRight';
    experiments = AP_find_experiments(animal,protocol);
    
    % Get days with muscimol
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    muscimol_start_day = muscimol(muscimol_animal_idx).day{1};
    muscimol_experiments = datenum({experiments.day})' >= datenum(muscimol_start_day);
    
    % Set experiments to use (imaging, not muscimol)
    experiments = experiments([experiments.imaging] & ~muscimol_experiments);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        AP_load_experiment;
        
        % Pull out trial data
        operant_grab_trial_data;
        
        % Store trial data into master structure
        trial_data_fieldnames = fieldnames(trial_data);
        for curr_trial_data_field = trial_data_fieldnames'
            trial_data_all.(cell2mat(curr_trial_data_field)){curr_animal,1}{curr_day,1} = ...
                trial_data.(cell2mat(curr_trial_data_field));
        end
        
        % Store general info
        trial_data_all.animals = animals;
        trial_data_all.t = t;
        trial_data_all.task_regressor_labels = task_regressor_labels;
        trial_data_all.task_regressor_sample_shifts = task_regressor_sample_shifts;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = ['trial_activity_task_corticostriatal'];
save([save_path filesep save_fn],'-v7.3');


%% Passive - tetO

clear all
disp('Passive trial activity (tetO)')

animals = {'AP100','AP101','AP103','AP104','AP105','AP106','AP107','AP108','AP109'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);
    
    % Get days with muscimol
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    if any(muscimol_animal_idx)
        muscimol_start_day = muscimol(muscimol_animal_idx).day{1};
        muscimol_experiments = datenum({experiments.day})' >= datenum(muscimol_start_day);
    else
        muscimol_experiments = false(size({experiments.day}));
    end
    
    % Set experiments to use (imaging, not muscimol)
    experiments = experiments([experiments.imaging] & ~muscimol_experiments);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        AP_load_experiment;
        
        % Pull out trial data
        operant_grab_trial_data;
        
        % Store trial data into master structure
        trial_data_fieldnames = fieldnames(trial_data);
        for curr_trial_data_field = trial_data_fieldnames'
            trial_data_all.(cell2mat(curr_trial_data_field)){curr_animal,1}{curr_day,1} = ...
                trial_data.(cell2mat(curr_trial_data_field));
        end
        
        % Store general info
        trial_data_all.animals = animals;
        trial_data_all.t = t;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
    
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = ['trial_activity_passive_teto'];
save([save_path filesep save_fn],'-v7.3');
disp(['Saved: ' save_path filesep save_fn])


%% Task - tetO

clear all
disp('Task trial activity (tetO)')

animals = {'AP100','AP101','AP103','AP104','AP105','AP106','AP107','AP108','AP109'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_stimWheelRight';
    experiments = AP_find_experiments(animal,protocol);
       
    % Get days with muscimol
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    if any(muscimol_animal_idx)
        muscimol_start_day = muscimol(muscimol_animal_idx).day{1};
        muscimol_experiments = datenum({experiments.day})' >= datenum(muscimol_start_day);
    else
        muscimol_experiments = false(size({experiments.day}));
    end
    
    % Set experiments to use (imaging, not muscimol)
    experiments = experiments([experiments.imaging] & ~muscimol_experiments);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        AP_load_experiment;
        
        % Pull out trial data
        operant_grab_trial_data;
        
        % Store trial data into master structure
        trial_data_fieldnames = fieldnames(trial_data);
        for curr_trial_data_field = trial_data_fieldnames'
            trial_data_all.(cell2mat(curr_trial_data_field)){curr_animal,1}{curr_day,1} = ...
                trial_data.(cell2mat(curr_trial_data_field));
        end
        
        % Store general info
        trial_data_all.animals = animals;
        trial_data_all.t = t;
        trial_data_all.task_regressor_labels = task_regressor_labels;
        trial_data_all.task_regressor_sample_shifts = task_regressor_sample_shifts;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = ['trial_activity_task_teto'];
save([save_path filesep save_fn],'-v7.3');

%% Muscimol: Passive - tetO

clear all
disp('Muscimol: Passive trial activity (tetO)')

animals = {'AP100','AP101','AP103','AP104','AP105','AP106'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);
    
    % Get days with muscimol
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    muscimol_experiments = ismember({experiments.day},muscimol(muscimol_animal_idx).day);
    
    % Set experiments to use (imaging, not muscimol)
    experiments = experiments([experiments.imaging] & muscimol_experiments);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        AP_load_experiment;
        
        % Pull out trial data
        operant_grab_trial_data;
        
        % Store trial data into master structure
        trial_data_fieldnames = fieldnames(trial_data);
        for curr_trial_data_field = trial_data_fieldnames'
            trial_data_all.(cell2mat(curr_trial_data_field)){curr_animal,1}{curr_day,1} = ...
                trial_data.(cell2mat(curr_trial_data_field));
        end
        
        % Store general info
        trial_data_all.animals = animals;
        trial_data_all.t = t;
        
        % Store muscimol info
        muscimol_day_idx = ismember(muscimol(muscimol_animal_idx).day,day);
        trial_data_all.muscimol_area{curr_animal}{curr_day} = ...
            muscimol(muscimol_animal_idx).area{muscimol_day_idx};
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
    
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = ['trial_activity_passive_teto_muscimol'];
save([save_path filesep save_fn],'-v7.3');
disp(['Saved: ' save_path filesep save_fn])


%% Muscimol: Task - tetO

clear all
disp('Muscimol: Task trial activity (tetO)')

animals = {'AP100','AP101','AP103','AP104','AP105','AP106'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_stimWheelRight';
    experiments = AP_find_experiments(animal,protocol);
       
    % Get days with muscimol
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    muscimol_start_day = muscimol(muscimol_animal_idx).day{1};   
    muscimol_experiments = ismember({experiments.day},muscimol(muscimol_animal_idx).day);
    
    % Set experiments to use (imaging, not muscimol)
    experiments = experiments([experiments.imaging] & muscimol_experiments);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        AP_load_experiment;
        
        % Pull out trial data
        operant_grab_trial_data;
        
        % Store trial data into master structure
        trial_data_fieldnames = fieldnames(trial_data);
        for curr_trial_data_field = trial_data_fieldnames'
            trial_data_all.(cell2mat(curr_trial_data_field)){curr_animal,1}{curr_day,1} = ...
                trial_data.(cell2mat(curr_trial_data_field));
        end
        
        % Store general info
        trial_data_all.animals = animals;
        trial_data_all.t = t;
        trial_data_all.task_regressor_labels = task_regressor_labels;
        trial_data_all.task_regressor_sample_shifts = task_regressor_sample_shifts;
        
        % Store muscimol info
        muscimol_day_idx = ismember(muscimol(muscimol_animal_idx).day,day);
        trial_data_all.muscimol_area{curr_animal}{curr_day} = ...
            muscimol(muscimol_animal_idx).area{muscimol_day_idx};
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = ['trial_activity_task_teto_muscimol'];
save([save_path filesep save_fn],'-v7.3');

%% Muscimol: Passive - corticostriatal

clear all
disp('Muscimol: Passive trial activity (cstr)')

animals = {'AP092','AP093','AP094','AP095','AP096','AP097'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);
    
    % Get days with muscimol
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    muscimol_experiments = ismember({experiments.day},muscimol(muscimol_animal_idx).day);
    
    % Set experiments to use (imaging, not muscimol)
    experiments = experiments([experiments.imaging] & muscimol_experiments);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        AP_load_experiment;
        
        % Pull out trial data
        operant_grab_trial_data;
        
        % Store trial data into master structure
        trial_data_fieldnames = fieldnames(trial_data);
        for curr_trial_data_field = trial_data_fieldnames'
            trial_data_all.(cell2mat(curr_trial_data_field)){curr_animal,1}{curr_day,1} = ...
                trial_data.(cell2mat(curr_trial_data_field));
        end
        
        % Store general info
        trial_data_all.animals = animals;
        trial_data_all.t = t;
        
        % Store muscimol info
        muscimol_day_idx = ismember(muscimol(muscimol_animal_idx).day,day);
        trial_data_all.muscimol_area{curr_animal}{curr_day} = ...
            muscimol(muscimol_animal_idx).area{muscimol_day_idx};
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
    
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = ['trial_activity_passive_cstr_muscimol'];
save([save_path filesep save_fn],'-v7.3');
disp(['Saved: ' save_path filesep save_fn])



%% Passive: ephys

animals = {'AP100','AP101','AP104','AP105','AP106'};

stim_mua = cell(length(animals),1);

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
%     protocol = 'AP_lcrGratingPassive';
    protocol = 'AP_stimWheelRight';
    experiments = AP_find_experiments(animal,protocol);
    
    % Set experiments to use (ephys)
    experiments = experiments([experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        AP_load_experiment;
        
        % QUICK ANALYSIS        
        
        % Get event-aligned activity
        raster_window = [-0.5,2];
        raster_sample_rate = 50;
        
        raster_sample_time = 1/raster_sample_rate;
        t = raster_window(1):raster_sample_time:raster_window(2);
        
        % Get align times
        use_align = stimOn_times;
        use_align(isnan(use_align)) = 0;
        
        t_peri_event = use_align + t;
        t_peri_event_bins = [t_peri_event-raster_sample_time/2,t_peri_event(:,end)+raster_sample_time/2];

        % Align wheel move
        event_aligned_wheelmove = interp1(Timeline.rawDAQTimestamps, ...
            +wheel_move,t_peri_event,'previous');
        
        % Align multiunit by depth       
        n_depths = 4;
        depth_group_edges = linspace(1000,max(channel_positions(:,2)),n_depths+1);
        depth_group_centers = round(depth_group_edges(1:end-1)+diff(depth_group_edges)/2);
        depth_group = discretize(spike_depths,depth_group_edges);
               
        event_aligned_mua = nan(length(stimOn_times),length(t),n_depths);
        for curr_depth = 1:n_depths
            curr_spikes = spike_times_timeline(depth_group == curr_depth);
            
            % (skip if no spikes at this depth)
            if isempty(curr_spikes)
                continue
            end
            
            event_aligned_mua(:,:,curr_depth) = cell2mat(arrayfun(@(x) ...
                histcounts(curr_spikes,t_peri_event_bins(x,:)), ...
                [1:size(t_peri_event,1)]','uni',false))*raster_sample_rate;
        end
        
        % Store average across trials
        quiescent_trials = ~any(event_aligned_wheelmove(:,t > 0 & t < 0.5),2);
        use_trials = quiescent_trials & stimIDs == 3;
        
        curr_mua = nanmean(event_aligned_mua(use_trials,:,:),1);

        stim_mua{curr_animal}(curr_day,:,:) = curr_mua;
              
        figure;plot(t,squeeze(curr_mua));
        drawnow;
        
        % Prep next loop
        AP_print_progress_fraction(curr_day,length(experiments));     
        clearvars('-except',preload_vars{:});
        
    end
end

disp('Finished loading all')


raster_window = [-0.5,2];
raster_sample_rate = 50;

raster_sample_time = 1/raster_sample_rate;
t = raster_window(1):raster_sample_time:raster_window(2);

a = cell2mat(stim_mua);
a_baseline = nanmean(a(:,t<0,:),2);
a2 = permute(nanmean((a - a_baseline)./a_baseline,1),[3,2,1]);

figure;plot(t,a2')



%% ~~~~~~~~~ BATCH ANALYSIS ~~~~~~~~~

%% >> Passive 

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_passive_tetO';
% data_fn = 'trial_activity_passive_cstr';

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

% (average stim response for each day)
use_stim = 3;

n_days = cellfun(@length,trial_info_all);
max_days = max(n_days);
stim_v_avg_dayavg = nanmean(cell2mat(permute(cellfun(@(x) ...
    padarray(x,[0,0,max_days-size(x,3),0],NaN,'post'), ...
    stim_v_avg,'uni',false),[1,3,4,5,2])),5);
stim_px_avg_dayavg = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    stim_v_avg_dayavg(:,:,:,use_stim));
AP_image_scroll(stim_px_avg_dayavg,t);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

% (plot time average by day)
use_t = t >= 0.1 & t <= 0.2;
stim_px_avg_dayavg_tavg = squeeze(nanmean(stim_px_avg_dayavg(:,:,use_t,:),3));

min_days = min(n_days);
c = repmat(prctile(stim_px_avg_dayavg_tavg(:),95),1,2).*[-1,1];
figure;
for curr_day = 1:min_days
    subplot(1,min_days,curr_day);
    imagesc(stim_px_avg_dayavg_tavg(:,:,curr_day));
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    caxis(c)
    axis image off;
    colormap(brewermap([],'PrGn'));
    if curr_day <= 3
        title('Passive')
    else
       title(['Training ' num2str(curr_day-3)]) 
    end
end


% (average stim response "pre/post learning")
naive_days = 1:3;
unlearned_days = 4:6;
expert_days = 6:10;

stim_avg_naive = nanmean(stim_px_avg_dayavg(:,:,:,naive_days),4);
stim_avg_unlearned = nanmean(stim_px_avg_dayavg(:,:,:,unlearned_days),4);
stim_avg_learned = nanmean(stim_px_avg_dayavg(:,:,:,expert_days),4);

AP_image_scroll([stim_avg_naive,stim_avg_unlearned,stim_avg_learned]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;

use_t = t >= 0.05 & t <= 0.2;
stim_learning_tavg = cat(3,nanmean(stim_avg_naive(:,:,use_t),3), ...
    nanmean(stim_avg_unlearned(:,:,use_t),3),nanmean(stim_avg_learned(:,:,use_t),3));
c = repmat(prctile(stim_learning_tavg(:),95),1,2).*[-1,1];

figure; 
for i = 1:3
    subplot(1,3,i)
    imagesc(stim_learning_tavg(:,:,i));
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    caxis(c)
    axis image off;
    colormap(brewermap([],'PrGn'));
    switch i
        case 1
            title('Naive');
        case 2
            title('Early training');
        case 3
            title('Late training');
    end
end


% Plot stim response in single animal over days
use_animal = 8;
curr_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),stim_v_avg{use_animal}(:,:,:,3));
AP_image_scroll(curr_px,t);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
axis image;
set(gcf,'name',animals{use_animal});


%% Time-average ROI activity by day

curr_act = fluor_roi_deconv(:,:,7);
% curr_act = fluor_roi_deconv(:,:,7) - fluor_roi_deconv(:,:,17);

use_trials = trial_stim_allcat == 1 & quiescent_trials;

use_t = t >= 0.1 & t <= 0.2;
curr_act_timeavg = nanmean(double(curr_act(:,use_t)),2);

roi_animal_day_avg = accumarray([trial_animal(use_trials),trial_day(use_trials)], ...
    curr_act_timeavg(use_trials),[max(trial_animal),max(trial_day)], ...
    @nanmean,NaN);
figure; hold on;
plot([1:max(trial_day)]-3,roi_animal_day_avg');
plot([1:max(trial_day)]-3,nanmean(roi_animal_day_avg,1),'k','linewidth',2);

% Violin plot
figure; 
distributionPlot(curr_act_timeavg(trial_animal == 1),'groups',trial_day(trial_animal == 1),'showMM',0);

figure;
plotSpread(curr_act_timeavg,'categoryIdx',trial_animal,'distributionIdx',trial_day);

figure;
curr_animal = 1;
curr_act_idx = trial_animal == curr_animal;
plotSpread(curr_act_timeavg(curr_act_idx),'distributionIdx',trial_day(curr_act_idx));


%%%% split days
n_daysplit = 4;
day_split_idx = cell2mat(arrayfun(@(x) ...
    min(floor(linspace(1,n_daysplit+1,x)),n_daysplit)', ...
    trials_recording,'uni',false));
roi_animal_daysplit_avg = accumarray( ...
    [trial_animal(use_trials),trial_day(use_trials),day_split_idx(use_trials)], ...
    curr_act_timeavg(use_trials),[max(trial_animal),max(trial_day),n_daysplit], ...
    @nanmean,NaN);
figure; hold on;
roi_animal_day_avg_long = reshape(permute( ...
    padarray(roi_animal_daysplit_avg,[0,0,1],NaN,'post'),[3,2,1]),[],length(animals));
% plot([1:size(roi_animal_day_avg_long,1)]/(n_daysplit+1),roi_animal_day_avg_long);
errorbar([1:size(roi_animal_day_avg_long,1)]/(n_daysplit+1) - 3, ...
    nanmean(roi_animal_day_avg_long,2),AP_sem(roi_animal_day_avg_long,2),'k','linewidth',2);



%% ROI-ROI correlation

use_t = t >= 0.1 & t <= 0.2;

% % (single stim - noise correlations)
% use_stim = 3;
% stim_roi_avg_t = cellfun(@(x) cellfun(@(x) ...
%     squeeze(nanmean(x{use_stim}(:,use_t,:),2)),x,'uni',false),stim_roi,'uni',false);

% (to combine stim - signal correlations)
a = cellfun(@(x) cellfun(@(x) vertcat(x{:}),x,'uni',false),stim_roi,'uni',false);
stim_roi_avg_t = cellfun(@(x) cellfun(@(x) squeeze(nanmean(x(:,use_t,:),2)),x,'uni',false),a,'uni',false);

curr_roi = 1;

stim_roi_avg_t_corr = cellfun(@(x) cellfun(@corrcoef,x,'uni',false),stim_roi_avg_t,'uni',false);
curr_roi_avg_corr = cellfun(@(x) cell2mat(cellfun(@(x) x(curr_roi,:),x,'uni',false)'),stim_roi_avg_t_corr,'uni',false);

figure; 

subplot(2,1,1); hold on;
compare_roi = 3;
for i = 1:length(curr_roi_avg_corr);plot(curr_roi_avg_corr{i}(:,compare_roi));end
ylabel([wf_roi(curr_roi).area '-' wf_roi(compare_roi).area]);
ylim([-0.2,1]);

subplot(2,1,2); hold on;
compare_roi = 7;
for i = 1:length(curr_roi_avg_corr);plot(curr_roi_avg_corr{i}(:,compare_roi));end
ylabel([wf_roi(curr_roi).area '-' wf_roi(compare_roi).area]);
ylim([-0.2,1]);





%% Pixel-behavior correlation 

% Get average stim (skip first 3 days: no behavior)
stim_px = cellfun(@(x) ...
    squeeze(AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    x(:,t >= 0 & t <= 1,4:end,3))), ...
    stim_v_avg,'uni',false);

stim_t = [0.1,0.2];
stim_px_tavg = cellfun(@(x) ...
    squeeze(AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    nanmean(x(:,t >= stim_t(1) & t <= stim_t(2),4:end,3),2))), ...
    stim_v_avg,'uni',false);


% Correlate behavior measure with pixels
stim_move_t_median = cellfun(@(x) cellfun(@nanmedian,x),{bhv.stim_move_t},'uni',false);
stim_reward_t_median = cellfun(@(x) cellfun(@nanmedian,x),{bhv.stim_feedback_t},'uni',false);
stim_move_reward_t_median = cellfun(@(x,y) cellfun(@(x,y) nanmedian(y-x),x,y),{bhv.stim_move_t},{bhv.stim_feedback_t},'uni',false);

rxn_frac = cellfun(@(x) cellfun(@(x) nanmean(x > 0.1 & x < 0.25),x), ...
    {bhv.stim_move_t},'uni',false);

stim_response_idx = {bhv.stim_response_idx};

stim_surround_t = bhv(curr_animal).stim_surround_t;
poststim_move_frac_max = cellfun(@(x) cellfun(@(x) ...
    max(nanmean(x(:,stim_surround_t > 0),1)),x), ...
    {bhv.stim_surround_wheel},'uni',false);

% use_bhv = rxn_frac;
use_bhv = stim_response_idx;

r_long = cell2mat(permute(cellfun(@(px,bhv) reshape(corr(reshape(px,[],size(px,3))', ...
    bhv(1:size(px,3))),[size(px,1),size(px,2)]),stim_px_tavg,use_bhv,'uni',false),[1,3,2]));
figure;imagesc(nanmean(r_long,3));
axis image off
caxis([-1,1]);
colormap(brewermap([],'*RdBu'));
title('Fluor:bhv corr')
AP_reference_outline('ccf_aligned','k');


% For more data: correlate all trials instead of avg
use_bhv = stim_move_t_median';

trial_bhv = cell(size(animals));
for curr_animal = 1:length(animals)
    curr_bhv = [nan(1,3),use_bhv{curr_animal}(1:n_days(curr_animal)-3)];
    trial_bhv{curr_animal} = cellfun(@(x,y) repmat(x,y,1), num2cell(curr_bhv)', ...
        cellfun(@(x) size(x,1),wheel_all{curr_animal},'uni',false),'uni',false);   
end
trial_bhv_allcat = cell2mat(vertcat(trial_bhv{:}));

use_trials = quiescent_trials & trial_stim_allcat == 1;

U_downsample_factor = 10;
Ud = imresize(U_master,1/U_downsample_factor,'bilinear');
trial_px = AP_svdFrameReconstruct(Ud(:,:,1:n_vs), ...
    permute(fluor_allcat_deconv(use_trials,:,:),[3,2,1]));

corr_px = nan(size(Ud,1),size(Ud,2),length(t));
corr_px_p = nan(size(Ud,1),size(Ud,2),length(t));
for curr_t = 1:length(t)
    [r,p] = corr(reshape(trial_px(:,:,curr_t,:),[],sum(use_trials))', ...
        trial_bhv_allcat(use_trials),'rows','complete');
    corr_px(:,:,curr_t) = reshape(r,size(Ud,1),size(Ud,2));
end
AP_image_scroll(corr_px,t);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'*RdBu'));
axis image;


%% Average passive by behavior

% Load and prep behavior
data_fn_parse = regexp(data_fn,'_','split');
bhv_fn = [trial_data_path filesep 'bhv_' data_fn_parse{end}];
load(bhv_fn);

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
    
    pre_muscimol_day_idx = datenum(bhv(curr_animal).day) < ...
        datenum(muscimol(muscimol_animal_idx).day(1));
    muscimol_day_idx = ismember(datenum(bhv(curr_animal).day), ...
        datenum(muscimol(muscimol_animal_idx).day));
    
    use_days{curr_animal} = pre_muscimol_day_idx;    
end

% Stim response index
stim_response_idx_prepad = cellfun(@(x,use_days) padarray(x(use_days),[3,0],NaN,'pre'), ...
    {bhv.stim_response_idx},use_days,'uni',false)';
stim_response_idx_prepadcat = cell2mat(stim_response_idx_prepad);

% Group stim responses by behavioral performance bins
stim_v_avg_cat = cell2mat(permute(stim_v_avg,[1,3,2,4])); 

stim_response_idx_padcat_grp = ...
    discretize(stim_response_idx_prepadcat,[-Inf,0,0.2,Inf]);
stim_response_idx_padcat_grp(isnan(stim_response_idx_padcat_grp)) = 0;
stim_v_avg_cat_grp = reshape(grpstats( ...
    reshape(permute(stim_v_avg_cat,[1,2,4,3]),[],size(stim_v_avg_cat,3))', ...
    stim_response_idx_padcat_grp)',n_vs,length(t),length(stim_unique),[]);

curr_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),stim_v_avg_cat_grp);
AP_image_scroll(reshape(permute(curr_px,[1,2,4,3,5]), ...
    size(curr_px,1),[],size(curr_px,3),size(curr_px,5)),t);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5],[], ...
    [size(U_master,1),size(U_master,2),1,size(curr_px,4)]);

% Plot time-average
use_t = t >= 0.1 & t <= 0.2;
curr_px_avg = squeeze(nanmean(curr_px(:,:,use_t,:,:),3));
c = prctile(curr_px_avg(:),99.5).*[-1,1];
figure;
for curr_stim = 1:length(stim_unique)
    for curr_learn = 1:size(curr_px_avg,3)
        subplot(length(stim_unique),size(curr_px_avg,3), ...
            size(curr_px_avg,3)*(curr_stim-1) + curr_learn);
        imagesc(curr_px_avg(:,:,curr_stim,curr_learn));
        caxis(c);
        colormap(brewermap([],'PrGn'));
        axis image off
        AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    end
end

% Get fluorescence in ROIs
stim_v_avg_cat_grp_roi = ...
    AP_svd_roi(U_master(:,:,1:n_vs),stim_v_avg_cat_grp,[],[], ...
    cat(3,wf_roi.mask));
    

% ROI fluorescence by behavior for each animal
use_t = t >= 0.1 & t <= 0.2;
use_stim = 1;
roi_act_timeavg = cellfun(@(x) squeeze(nanmean(x(:,use_t,:,use_stim),2)),stim_roi_avg,'uni',false)';

animal_colors = max(brewermap(length(animals),'Set3')-0.1,0);
figure;
for curr_roi = 1:numel(wf_roi)
    subplot(2,length(wf_roi),curr_roi);
    hold on; set(gca,'ColorOrder',animal_colors);
    cellfun(@(x,y) plot(x,normalize(y(curr_roi,:),'range'),'.','MarkerSize',20), ...
        stim_response_idx_prepad,roi_act_timeavg);
    xlabel('Stim response idx');
    ylabel('Fluor');
    title(wf_roi(curr_roi).area);
    axis square;
end
legend(animals,'location','nw');


% Align days across animals by "learned day"
curr_roi = 7;
use_t = t >= 0.1 & t <= 0.2;

curr_act_timeavg = cellfun(@(x) squeeze(nanmean(x(curr_roi,use_t,:,use_stim),2)),stim_roi_avg,'uni',false)';

max_days = max(cellfun(@length,curr_act_timeavg));

curr_act_timeavg_pad = cell2mat(cellfun(@(x) ...
    padarray(x,[max_days-length(x),0],NaN,'post'), ...
    curr_act_timeavg,'uni',false)');


stim_response_idx_pad = cell2mat(cellfun(@(x) padarray(cat(1,nan(3,1),x), ...
    [max_days-length(x)-3,0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_response_idx},use_days,'uni',false),'uni',false));

rxn_frac_pad = cell2mat(cellfun(@(x) padarray(cat(1,nan(3,1),cellfun(@(x) nanmean(x > 0.1 & x < 0.25),x)), ...
    [max_days-length(x)-3,0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_move_t},use_days,'uni',false),'uni',false));

alt_rxn_frac_pad = cell2mat(cellfun(@(x) padarray(cat(1,nan(3,1),cellfun(@(x) ...
    nanmean(cell2mat(x) > 0.1 & cell2mat(x) < 0.25),x)), ...
    [max_days-length(x)-3,0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.alt_stim_move_t},use_days,'uni',false),'uni',false));

rxn_frac_pad_diff = rxn_frac_pad - alt_rxn_frac_pad;

[~,learned_day] = max(rxn_frac_pad>0.4,[],1);
% [~,learned_day] = max(rxn_frac_pad_diff>0.3,[],1);
% [~,learned_day] = max(stim_response_idx_pad>0.2,[],1);
learned_day_x = [1:max_days]'-learned_day;

[rxn_grp_mean,learned_day_grp] = ...
    grpstats(rxn_frac_pad(:),learned_day_x(:),{'nanmean','gname'});
rxndiff_grp_mean = ...
    grpstats(rxn_frac_pad_diff(:),learned_day_x(:),{'nanmean'});
stimresponse_grp_mean = ...
    grpstats(stim_response_idx_pad(:),learned_day_x(:),{'nanmean'});
act_grp_mean = grpstats(curr_act_timeavg_pad(:),learned_day_x(:),'nanmean');
learned_day_grp = cellfun(@str2num,learned_day_grp);

figure;
subplot(1,4,1); hold on;
plot(learned_day_x,rxn_frac_pad);
plot(learned_day_grp,rxn_grp_mean,'k','linewidth',2);
xlabel('Learned day');
ylabel('Rxn frac');
line([0,0],ylim,'color','k','linestyle','--');

subplot(1,4,2); hold on;
plot(learned_day_x,rxn_frac_pad_diff);
plot(learned_day_grp,rxndiff_grp_mean,'k','linewidth',2);
xlabel('Learned day');
ylabel('Rxn frac diff');
line([0,0],ylim,'color','k','linestyle','--');

subplot(1,4,3); hold on;
plot(learned_day_x,stim_response_idx_pad);
plot(learned_day_grp,stimresponse_grp_mean,'k','linewidth',2);
xlabel('Learned day');
ylabel('Stim response idx');
line([0,0],ylim,'color','k','linestyle','--');

subplot(1,4,4); hold on;
plot(learned_day_x,curr_act_timeavg_pad);
plot(learned_day_grp,act_grp_mean,'k','linewidth',2);
xlabel('Learned day');
ylabel(sprintf('%s fluorescence',wf_roi(curr_roi).area))
line([0,0],ylim,'color','k','linestyle','--');

linkaxes(get(gcf,'Children'),'x')


% (average stim response relative to "learned day")
use_stim = 3;

learned_days_animal = cellfun(@(x,y) (1:size(x,3))-y, ...
    stim_v_avg,num2cell(learned_day),'uni',false);

stim_v_avg_learnedday = reshape(grpstats( ...
    reshape(cell2mat(permute(cellfun(@(x) x(:,:,:,use_stim), ...
    stim_v_avg,'uni',false),[1,3,2])),n_vs*length(t),[])', ...
    cell2mat(learned_days_animal)')',n_vs,length(t),[]);

stim_px_avg_learnedday = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    stim_v_avg_learnedday);

AP_image_scroll(stim_px_avg_learnedday,t);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);






%% >> (Choiceworld task: post-learn only)

trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\paper_unused';
data_fn = 'trial_activity_choiceworld_wfonly';

AP_load_trials_wf;


%% >> Task

% Task
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_task_teto';
% data_fn = 'trial_activity_task_corticostriatal';

AP_load_trials_wf;

% Get animal and day index for each trial
trial_animal = cell2mat(arrayfun(@(x) ...
    x*ones(size(vertcat(wheel_all{x}{:}),1),1), ...
    [1:length(wheel_all)]','uni',false));
trial_day = cell2mat(cellfun(@(x) cell2mat(cellfun(@(curr_day,x) ...
    curr_day*ones(size(x,1),1),num2cell(1:length(x))',x,'uni',false)), ...
    wheel_all,'uni',false));

% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_animal;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));


% Load and prep behavior
data_fn_parse = regexp(data_fn,'_','split');
bhv_fn = [trial_data_path filesep 'bhv_' data_fn_parse{end}];
load(bhv_fn);

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
    
    pre_muscimol_day_idx = datenum(bhv(curr_animal).day) < ...
        datenum(muscimol(muscimol_animal_idx).day(1));
    muscimol_day_idx = ismember(datenum(bhv(curr_animal).day), ...
        datenum(muscimol(muscimol_animal_idx).day));
    
    use_days{curr_animal} = pre_muscimol_day_idx;    
end

% Stim to move time (NOT SAME NUM TRIALS - REWARDED ONLY)
stim_move_t = cell2mat(cellfun(@(x,y) cell2mat(x(y)), ...
    {bhv.stim_move_t},use_days,'uni',false)');

% Rxn stim fraction
rxn_frac = cell2mat(cellfun(@(x) cellfun(@(x) nanmean(x > 0.1 & x < 0.25),x), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_move_t},use_days,'uni',false),'uni',false)');
rxn_frac_allcat = cell2mat(cellfun(@(x,y) repmat(x,y,1), ...
    num2cell(rxn_frac),num2cell(trials_recording),'uni',false));

% Stim response index
stim_response_idx = cell2mat(cellfun(@(x,y) x(y), ...
    {bhv.stim_response_idx},use_days,'uni',false)');
stim_response_idx_allcat = cell2mat(cellfun(@(x,y) repmat(x,y,1), ...
    num2cell(stim_response_idx),num2cell(trials_recording),'uni',false));




%% Average trial activity by day

% % (move-align fluor?)
% fluor_allcat_deconv_move = fluor_allcat_deconv;
% fluor_roi_deconv_move = fluor_roi_deconv;
% t_leeway = -t(1);
% leeway_samples = round(t_leeway*(sample_rate));
% for i = 1:size(fluor_allcat_deconv,1)
%     fluor_allcat_deconv_move(i,:,:,:) = circshift(fluor_allcat_deconv_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
%     fluor_roi_deconv_move(i,:,:,:) = circshift(fluor_roi_deconv_move(i,:,:,:),-move_idx(i)+leeway_samples,2);
% end

% Average activity within animal/day
n_days = cellfun(@length,trial_info_all);
max_days = max(trial_day);
fluor_deconv_cat = nan(n_vs,length(t),max_days,length(animals));
for curr_animal = 1:length(animals)
    for curr_day = 1:n_days(curr_animal)
        % (all trials)
        curr_trials = trial_animal == curr_animal & trial_day == curr_day;
        
%         % (full activity)
%         fluor_deconv_cat(:,:,curr_day,curr_animal) = ...
%             permute(nanmean(fluor_allcat_deconv(curr_trials,:,:),1),[3,2,1]);
        
%         % (task-predicted activity)
%         fluor_deconv_cat(:,:,curr_day,curr_animal) = ...
%             permute(nanmean(fluor_taskpred_allcat(curr_trials,:,:),1),[3,2,1]);
        
        % (reduced activity)
        fluor_deconv_cat(:,:,curr_day,curr_animal) = ...
            permute(nanmean(fluor_allcat_deconv(curr_trials,:,:) - ...
            fluor_taskpred_reduced_allcat(curr_trials,:,:,1),1),[3,2,1]);
        
    end
end

% Plot average for each day
curr_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),nanmean(fluor_deconv_cat,4));
AP_image_scroll(curr_px,t);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'))
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);





%% Average trial activity by rxn/bhv

% Average all activity by reaction time depending on learned
% rxn_bins = [-0.05:0.02:1]';
% rxn_bins = [-Inf,Inf];
% rxn_bins = [0.1,0.5];
rxn_bins = [0:0.1:0.5];
rxn_bin_centers = rxn_bins(1:end-1) + diff(rxn_bins)./2;
move_t_discretize = discretize(move_t,rxn_bins);
fluor_rxn = nan(n_vs,length(t),length(rxn_bins)-1);
for curr_rxn = 1:length(rxn_bins)-1
        use_trials = move_t_discretize == curr_rxn & ...
            trial_outcome_allcat == 1;
        
%     use_trials = move_t_discretize == curr_rxn & ...
%         trial_outcome_allcat == 1 & ...
%         trial_stim_allcat == 1;
     
    % (raw activity)
    fluor_rxn(:,:,curr_rxn) = ...
        permute(nanmean(fluor_allcat_deconv(use_trials,:,:),1),[3,2,1]);    
    
%     % (raw activity - move aligned)
%     fluor_rxn(:,:,curr_rxn) = ...
%         permute(nanmean(fluor_allcat_deconv_move(use_trials,:,:),1),[3,2,1]);
    
%     % (activity - taskpred_reduced)
%         fluor_rxn(:,:,curr_rxn) = ...
%             permute(nanmean(fluor_allcat_deconv(use_trials,:,:) - ...
%             fluor_taskpred_reduced_allcat(use_trials,:,:,1),1),[3,2,1]);
        
end

curr_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),fluor_rxn);
AP_image_scroll(curr_px,t);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'))
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);


%% Reduced activity by behavior (animal day average)

% Get average fluorescence by animal, day, stim
stim_unique = unique(trial_stim_allcat);
stim_v_avg = cell(size(animals));
for curr_animal = 1:length(animals)        
    for curr_day = 1:max(trial_day(trial_animal == curr_animal))
        for curr_stim_idx = 1:length(stim_unique)
            use_trials = ...
                move_t >= 0.1 & ...
                trial_animal == curr_animal & ...
                trial_day == curr_day & ...
                trial_stim_allcat == stim_unique(curr_stim_idx);
            stim_v_avg{curr_animal}(:,:,curr_day,curr_stim_idx) = ...
                permute(nanmean( ...
                fluor_allcat_deconv(use_trials,:,:) - ...
                fluor_taskpred_reduced_allcat(use_trials,:,:,1),1),[3,2,1]);
        end       
    end
end

% Group stim responses by behavioral performance
stim_v_avg_cat = cell2mat(permute(stim_v_avg,[1,3,2,4])); 

move_prepost_max_ratio_grp = ...
    discretize(move_prepost_max_ratio,[-Inf,0.2,Inf]);
move_prepost_max_ratio_grp(isnan(move_prepost_max_ratio)) = 0;
stim_v_avg_cat_grp = reshape(grpstats( ...
    reshape(stim_v_avg_cat,[],size(stim_v_avg_cat,3))', ...
    move_prepost_max_ratio_grp)',n_vs,length(t),[]);

curr_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),stim_v_avg_cat_grp);
AP_image_scroll(curr_px,t);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

% Plot time-average
use_t = t >= 0.1 & t <= 0.2;
curr_px_avg = squeeze(nanmean(curr_px(:,:,use_t,:,:),3));
c = prctile(curr_px_avg(:),99.5).*[-1,1];
figure;
for curr_learn = 1:size(curr_px_avg,3)
    subplot(1,size(curr_px_avg,3),curr_learn);
    imagesc(curr_px_avg(:,:,curr_learn));
    caxis(c);
    colormap(brewermap([],'PrGn'));
    axis image off
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
end




%% Average trial activity by wheel velocity

use_trials = trial_outcome_allcat == 1 & move_t > 0.1 & move_t < 0.2;

wheel_vel_summed = sum(abs(wheel_allcat(:,t > 0 & t < 0.5)),2);

wheel_vel_bins = prctile(wheel_vel_summed(use_trials),[0:10:100]);
wheel_vel_discretize = discretize(wheel_vel_summed,wheel_vel_bins);

fluor_binned = nan(n_vs,length(t),length(wheel_vel_bins)-1);
for curr_bin = 1:length(wheel_vel_bins)-1
    curr_trials = use_trials & wheel_vel_discretize == curr_bin;
    % (raw activity)
    fluor_binned(:,:,curr_bin) = ...
        permute(nanmean(fluor_allcat_deconv(curr_trials,:,:),1),[3,2,1]);
    % (task reduced activity)
%     fluor_binned(:,:,curr_bin) = ...
%         permute(nanmean(fluor_allcat_deconv(curr_trials,:,:) - ...
%         fluor_taskpred_reduced_allcat(curr_trials,:,:,1),1),[3,2,1]);
end

curr_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),fluor_binned);
AP_image_scroll(curr_px,t);
axis image;
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'))
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);


%% Average ROI activity by reaction time

use_trials = trial_outcome_allcat == 1 & stim_response_idx_allcat > 0.2;

% curr_act = fluor_roi_deconv(use_trials,:,7);
% curr_act = fluor_roi_taskpred(use_trials,:,1);
% curr_act = fluor_roi_deconv(use_trials,:,7) - fluor_roi_deconv(use_trials,:,17);
% curr_act = fluor_roi_taskpred(use_trials,:,4) - fluor_roi_taskpred(use_trials,:,14);
curr_act = fluor_roi_deconv(use_trials,:,7) - fluor_roi_taskpred_reduced(use_trials,:,7,1);

rxn_bins = [0:0.01:0.5]';
rxn_bin_centers = rxn_bins(1:end-1) + diff(rxn_bins)./2;
move_t_discretize = discretize(move_t(use_trials),rxn_bins);

% By each time point
curr_act_rxn = accumarray( ...
    {reshape(repmat(1:length(t),sum(~isnan(move_t_discretize)),1),[],1), ...
    reshape(repmat(move_t_discretize(~isnan(move_t_discretize)),1,length(t)),[],1)}, ...
    reshape(curr_act(~isnan(move_t_discretize),:),[],1), ...
    [length(t),length(rxn_bin_centers)],@nanmean,single(NaN));

figure;imagesc(t,rxn_bin_centers,curr_act_rxn');
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
set(gca,'YDir','normal')
xlabel('Activity time from stim');
ylabel('Reaction time');

%% Time-average ROI activity by day/behavior

% curr_act = fluor_roi_deconv(:,:,7);
% curr_act = fluor_roi_deconv(:,:,9) - fluor_roi_deconv(:,:,18);
curr_act = fluor_roi_deconv(:,:,7) - fluor_roi_taskpred_reduced(:,:,7,1);

use_t = t >= 0.1 & t <= 0.2;
curr_act_timeavg = nanmean(double(curr_act(:,use_t)),2);

use_trials = true(size(move_t));
% use_trials = move_t >= 0.3;

roi_animal_day_avg = accumarray([trial_animal(use_trials),trial_day(use_trials)], ...
    curr_act_timeavg(use_trials),[max(trial_animal),max(trial_day)], ...
    @nanmean,NaN);
figure; hold on;
plot(roi_animal_day_avg');
plot(nanmean(roi_animal_day_avg,1),'k','linewidth',2);

%%%% split days
n_daysplit = 4;
day_split_idx = cell2mat(arrayfun(@(x) ...
    min(floor(linspace(1,n_daysplit+1,x)),n_daysplit)', ...
    trials_recording,'uni',false));
roi_animal_daysplit_avg = accumarray( ...
    [trial_animal(use_trials),trial_day(use_trials),day_split_idx(use_trials)], ...
    curr_act_timeavg(use_trials),[max(trial_animal),max(trial_day),n_daysplit], ...
    @nanmean,NaN);
figure; hold on;
roi_animal_daysplit_avg_long = reshape(permute( ...
    padarray(roi_animal_daysplit_avg,[0,0,1],NaN,'post'),[3,2,1]),[],length(animals));
% plot([1:size(roi_animal_day_avg_long,1)]/(n_daysplit+1),roi_animal_day_avg_long);
errorbar([1:size(roi_animal_daysplit_avg_long,1)]/(n_daysplit+1), ...
    nanmean(roi_animal_daysplit_avg_long,2),AP_sem(roi_animal_daysplit_avg_long,2),'k','linewidth',2);

% Split by pre/post move ratio
bhv_bins = [-1:0.1:1];
bhv_bin_centers = bhv_bins(1:end-1)+diff(bhv_bins)/2;
bhv_discretize = discretize(stim_response_idx_allcat,bhv_bins);

rxn_bins = [0:0.05:0.7];
rxn_bin_centers = rxn_bins(1:end-1) + diff(rxn_bins)./2;
rxn_discretize = discretize(move_t,rxn_bins);

use_trials = ~isnan(rxn_discretize);

curr_act_binned = accumarray( ...
    [bhv_discretize(use_trials),rxn_discretize(use_trials)], ...
    curr_act_timeavg(use_trials), ...
    [length(bhv_bin_centers),length(rxn_bin_centers)],@nanmean,NaN);


figure;imagesc(rxn_bin_centers,bhv_bin_centers,curr_act_binned);
xlabel('Reaction time');
ylabel('Move ratio');
colormap(hot);



%% Task kernels

% Get task>cortex parameters
n_regressors = length(task_regressor_labels);
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

% Get average task > cortex kernels separately for each day
n_days = cellfun(@length,trial_info_all);
max_days = max(n_days);

% Make grid of kernel/day/mouse
fluor_taskpred_k_cat = cell(n_regressors,max_days,length(animals));
for curr_animal = 1:length(animals)
    fluor_taskpred_k_cat(:,1:n_days(curr_animal),curr_animal) = ...
        horzcat(fluor_taskpred_k_all{curr_animal}{:});
end

% Get average kernel given behavior
rxn_frac_pad = cell2mat(cellfun(@(x) padarray(cellfun(@(x) nanmean(x > 0.1 & x < 0.25),x), ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_move_t},use_days,'uni',false),'uni',false));

use_kernel_days = rxn_frac_pad > 0.4;
curr_k = squeeze(fluor_taskpred_k_cat(1,:,:));
curr_k_avg_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    permute(nanmean(vertcat(curr_k{use_kernel_days}),1),[3,2,1]));
AP_image_scroll(curr_k_avg_px);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
axis image;

% Get kernel weights in ROIs
fluor_taskpred_k_cat_roi = cellfun(@(x) ...
    AP_svd_roi(U_master(:,:,1:n_vs),permute(x,[3,2,1]),[],[],cat(3,wf_roi.mask)), ...
    fluor_taskpred_k_cat,'uni',false);

% Average regressors by day
k_day_mean = cell(n_regressors,max_days);
for curr_regressor = 1:n_regressors
    for curr_day = 1:max_days
        k_day_mean{curr_regressor,curr_day} = permute(nanmean(cat(4, ...
            fluor_taskpred_k_cat{curr_regressor,curr_day,:}),4),[3,2,1]);
    end
end

% Get regressor pixels
k_day_mean_px = cellfun(@(x) AP_svdFrameReconstruct(U_master(:,:,1:n_vs),x), ...
    k_day_mean,'uni',false);

% Plot average kernels
% (by day)
curr_regressor = 1;
curr_subregressor = 1;
curr_regressor_plot = cell2mat(permute(cellfun(@(x) x(:,:,:,curr_subregressor), ...
    k_day_mean_px(curr_regressor,:),'uni',false),[1,3,4,2]));

AP_image_scroll(curr_regressor_plot);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
axis image;

figure;imagesc(reshape(squeeze(sum(curr_regressor_plot,3)),size(U_master,1),[]));
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image off;


% (by unlearned vs learned)
unlearned_days = 1:5;
learned_days = 6:10;
curr_regressor_unlearned = nanmean(curr_regressor_plot(:,:,:,unlearned_days),4);
curr_regressor_learned = nanmean(curr_regressor_plot(:,:,:,learned_days),4);
AP_image_scroll([curr_regressor_unlearned,curr_regressor_learned,curr_regressor_learned-curr_regressor_unlearned]);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;


% Plot kernel over days in single animal
use_animal = 4;
curr_k = cell2mat(cellfun(@(x) x{1}(1,:,:),fluor_taskpred_k_all{use_animal},'uni',false));
curr_k_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),permute(curr_k,[3,2,1]));
AP_image_scroll(curr_k_px);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;
set(gcf,'name',animals{use_animal})


%% Within-animal: correlate reaction time to trial activity?

% (package back into animal/day)
n_days = cellfun(@length,wheel_all);
fluor_deconv_animal = mat2cell(mat2cell(fluor_allcat_deconv,trials_recording,length(t),n_vs),n_days);
move_t_animal = mat2cell(mat2cell(move_t,trials_recording),n_days);

% Plot time to movement for each animal/day
figure('Name','Time to move'); p = tight_subplot(length(animals),1);
for curr_animal = 1:length(animals)
    subplot(p(curr_animal));
    
    plot(vertcat(move_t_animal{curr_animal}{:}),'.k');
    
    curr_n_trials = cumsum(cellfun(@length,move_t_animal{curr_animal}));
    for i = 1:length(curr_n_trials)
        line([curr_n_trials,curr_n_trials],ylim,'color','r','linewidth',2);
    end
    ylabel(animals(curr_animal));
    set(gca,'XTick',[]);
end

%% ROI activity aligned by "learned day"
% ROI fluorescence by behavior for each animal

curr_roi = 7;

% curr_act = fluor_roi_deconv(:,:,curr_roi);
curr_act = fluor_roi_deconv(:,:,curr_roi) - fluor_roi_taskpred_reduced(:,:,curr_roi,1);

use_t = t >= 0.1 & t <= 0.2;
curr_act_timeavg = nanmean(double(curr_act(:,use_t)),2);

use_trials = true(size(move_t));
% use_trials = move_t >= 0.3;

curr_act_timeavg_pad = accumarray([trial_day(use_trials),trial_animal(use_trials)], ...
    curr_act_timeavg(use_trials),[max(trial_day),max(trial_animal)], ...
    @nanmean,NaN);
figure; hold on;
plot(curr_act_timeavg_pad);
plot(nanmean(curr_act_timeavg_pad,2),'k','linewidth',2);

% Align days across animals by "learned day"
max_days = max(trial_day);

stim_response_idx_pad = cell2mat(cellfun(@(x) padarray(x, ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_response_idx},use_days,'uni',false),'uni',false));

rxn_frac_pad = cell2mat(cellfun(@(x) padarray(cellfun(@(x) nanmean(x > 0.1 & x < 0.25),x), ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_move_t},use_days,'uni',false),'uni',false));

alt_rxn_frac_pad = cell2mat(cellfun(@(x) padarray(cellfun(@(x) ...
    nanmean(cell2mat(x) > 0.1 & cell2mat(x) < 0.25),x), ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.alt_stim_move_t},use_days,'uni',false),'uni',false));

rxn_frac_pad_diff = rxn_frac_pad - alt_rxn_frac_pad;

[~,learned_day] = max(rxn_frac_pad>0.4,[],1);
% [~,learned_day] = max(rxn_frac_pad_diff>0.3,[],1);
% [~,learned_day] = max(stim_response_idx_pad>0.2,[],1);
learned_day_x = [1:max_days]'-learned_day;

[rxn_grp_mean,learned_day_grp] = ...
    grpstats(rxn_frac_pad(:),learned_day_x(:),{'nanmean','gname'});
rxndiff_grp_mean = ...
    grpstats(rxn_frac_pad_diff(:),learned_day_x(:),{'nanmean'});
stimresponse_grp_mean = ...
    grpstats(stim_response_idx_pad(:),learned_day_x(:),{'nanmean'});
act_grp_mean = grpstats(curr_act_timeavg_pad(:),learned_day_x(:),'nanmean');
learned_day_grp = cellfun(@str2num,learned_day_grp);

figure;
subplot(1,4,1); hold on;
plot(learned_day_x,rxn_frac_pad);
plot(learned_day_grp,rxn_grp_mean,'k','linewidth',2);
xlabel('Learned day');
ylabel('Rxn frac');
line([0,0],ylim,'color','k','linestyle','--');

subplot(1,4,2); hold on;
plot(learned_day_x,rxn_frac_pad_diff);
plot(learned_day_grp,rxndiff_grp_mean,'k','linewidth',2);
xlabel('Learned day');
ylabel('Rxn frac diff');
line([0,0],ylim,'color','k','linestyle','--');

subplot(1,4,3); hold on;
plot(learned_day_x,stim_response_idx_pad);
plot(learned_day_grp,stimresponse_grp_mean,'k','linewidth',2);
xlabel('Learned day');
ylabel('Stim response idx');
line([0,0],ylim,'color','k','linestyle','--');

subplot(1,4,4); hold on;
plot(learned_day_x,curr_act_timeavg_pad);
plot(learned_day_grp,act_grp_mean,'k','linewidth',2);
xlabel('Learned day');
ylabel(sprintf('%s fluorescence',wf_roi(curr_roi).area))
line([0,0],ylim,'color','k','linestyle','--');

linkaxes(get(gcf,'Children'),'x')


% Plot learning-aligned daysplit
n_daysplit = 4;
day_split_idx = cell2mat(arrayfun(@(x) ...
    min(floor(linspace(1,n_daysplit+1,x)),n_daysplit)', ...
    trials_recording,'uni',false));
roi_animal_daysplit_avg = accumarray( ...
    [trial_animal(use_trials),trial_day(use_trials),day_split_idx(use_trials)], ...
    curr_act_timeavg(use_trials),[max(trial_animal),max(trial_day),n_daysplit], ...
    @nanmean,NaN);

roi_animal_daysplit_avg_long = reshape(permute( ...
    padarray(roi_animal_daysplit_avg,[0,0,1],NaN,'post'),[3,2,1]),[],length(animals));

learned_day_x_daysplit = cell2mat(cellfun(@(x) x+(0:n_daysplit)'/(n_daysplit+1), ...
    num2cell(learned_day_x),'uni',false));

[roi_daysplit,roi_daysplit_grp] = ...
    grpstats(roi_animal_daysplit_avg_long(:), ...
    learned_day_x_daysplit(:),{'nanmean','gname'});
roi_daysplit_grp = cellfun(@str2num,roi_daysplit_grp);

figure; hold on;
plot(learned_day_x_daysplit,roi_animal_daysplit_avg_long);
plot(roi_daysplit_grp,roi_daysplit,'k','linewidth',2);
line([0,0],ylim,'color','k','linestyle','--');
xlabel('Learned day');
ylabel(sprintf('%s fluorescence',wf_roi(curr_roi).area))




%% >> Muscimol: task

trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_task_teto_muscimol';

AP_load_trials_wf;

% Choose split for data
trials_allcat = size(wheel_allcat,1);
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));
use_split = trials_animal;

split_idx = cell2mat(arrayfun(@(exp,trials) repmat(exp,trials,1), ...
    [1:length(use_split)]',reshape(use_split,[],1),'uni',false));

% Get task>cortex parameters
n_regressors = length(task_regressor_labels);
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

% Average kernels by muscimol area
muscimol_area_cat = horzcat(muscimol_area{:})';
unique_muscimol_area = unique(muscimol_area_cat);

muscimol_k = cell(n_regressors,1);
for curr_area = 1:length(unique_muscimol_area)
    curr_exps = cellfun(@(x) strcmp(x,unique_muscimol_area{curr_area}), ...
        muscimol_area,'uni',false)';
    curr_k = cellfun(@(x,y) x(y),fluor_taskpred_k_all,curr_exps,'uni',false);
    curr_k_cat = vertcat(curr_k{:});
    for curr_regressor = 1:n_regressors
        muscimol_k{curr_regressor}(:,:,:,curr_area) = ...
            nanmean(cell2mat(permute(cellfun(@(x) x{curr_regressor}, ...
            curr_k_cat,'uni',false),[2,3,4,1])),4);
    end
end

curr_regressor = 1;
curr_subregressor = 1;
curr_k_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ... 
    permute(muscimol_k{curr_regressor}(curr_subregressor,:,:,:),[3,2,4,1]));
AP_image_scroll(curr_k_px);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

% Plot kernel over days in single animal
use_animal = 3;
curr_regressor = 1;
curr_subregressor = 1;

curr_k = cell2mat(cellfun(@(x) x{curr_regressor}(curr_subregressor,:,:), ...
    fluor_taskpred_k_all{use_animal},'uni',false));
curr_k_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),permute(curr_k,[3,2,1]));
AP_image_scroll(curr_k_px);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;


%% >> Muscimol: passive

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_passive_teto_muscimol';
% data_fn = 'trial_activity_passive_cstr_muscimol';

AP_load_trials_wf;

% (package deconv back into animal/day)
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));

fluor_deconv_exp = mat2cell(fluor_allcat_deconv,trials_recording,length(t),n_vs);
fluor_roi_exp = mat2cell(fluor_roi_deconv,trials_recording,length(t),n_rois);
muscimol_area_cat = horzcat(muscimol_area{:})';

quiescent_trials = ~any(abs(wheel_allcat(:,t >= 0 & t <= 0.5)) > 0,2);
quiescent_trials_exp = mat2cell(quiescent_trials,trials_recording,1);
trial_stim_exp = mat2cell(trial_stim_allcat,trials_recording,1);

fluor_deconv_exp_stimavg = cell2mat(cellfun(@(act,quiescent,stim) ...
    nanmean(act(quiescent & stim == 1,:,:),1), ...
    fluor_deconv_exp,quiescent_trials_exp,trial_stim_exp,'uni',false));
fluor_roi_exp_stimavg = cell2mat(cellfun(@(act,quiescent,stim) ...
    nanmean(act(quiescent & stim == 1,:,:),1), ...
    fluor_roi_exp,quiescent_trials_exp,trial_stim_exp,'uni',false));

unique_muscimol_area = flipud(unique(muscimol_area_cat));
fluor_muscimol_avg = nan(n_vs,length(t),length(unique_muscimol_area));
for curr_area = 1:length(unique_muscimol_area)
    curr_exps = strcmp(muscimol_area_cat,unique_muscimol_area{curr_area});
    fluor_muscimol_avg(:,:,curr_area) = ...
        permute(nanmean(fluor_deconv_exp_stimavg(curr_exps,:,:),1),[3,2,1]);
end

fluor_muscimol_avg_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    fluor_muscimol_avg);

AP_image_scroll(fluor_muscimol_avg_px)
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

% Plot time averages and ROIs
use_t = t >= 0.1 & t <= 0.2;
fluor_roi_stimavg = squeeze(nanmean(fluor_roi_exp_stimavg(:,use_t,:),2));
fluor_muscimol_px_stimavg = squeeze(nanmean(fluor_muscimol_avg_px(:,:,use_t,:),3));

figure;
c = prctile(fluor_muscimol_px_stimavg(:),100).*[-1,1];
for curr_area = 1:length(unique_muscimol_area)
    subplot(1,length(unique_muscimol_area)+1,curr_area);
    imagesc(fluor_muscimol_px_stimavg(:,:,curr_area));
    caxis(c);
    colormap(brewermap([],'PrGn'));
    axis image off;
    AP_reference_outline('ccf_aligned','k');
    title(unique_muscimol_area{curr_area});
end

subplot(1,length(unique_muscimol_area)+1,length(unique_muscimol_area)+1);
plot_rois = [1,7];
area_col = lines(length(unique_muscimol_area));
gscatter(fluor_roi_stimavg(:,plot_rois(1)),fluor_roi_stimavg(:,plot_rois(2)), ...
    muscimol_area_cat,area_col);
axis square;
xlabel(wf_roi(plot_rois(1)).area);
ylabel(wf_roi(plot_rois(2)).area);



% Plot ROI pre/post inactivation
n_days = cellfun(@length,trial_info_all);

fluor_roi_stimavg_animal = reshape(mat2cell(permute(fluor_roi_exp_stimavg,[3,2,1]), ...
    n_rois,length(t),n_days),[],1);

use_t = t >= 0.1 & t <= 0.2;
plot_roi = 7;

fluor_roi_stimavg_animal_tavg = cellfun(@(x) squeeze(nanmean(x(plot_roi,use_t,:),2)), ...
    fluor_roi_stimavg_animal,'uni',false);

figure; hold on;
for curr_animal = 1:length(animals)
    curr_muscimol = find(strcmp(lower(muscimol_area{curr_animal}),'v1'),1,'last');
    % (get first washout after muscimol)
    curr_washout = curr_muscimol + ...
        find(strcmp(lower(muscimol_area{curr_animal}(curr_muscimol+1:end)),'washout'),1,'first');
    
    plot(fluor_roi_stimavg_animal_tavg{curr_animal}(curr_muscimol), ...
        fluor_roi_stimavg_animal_tavg{curr_animal}(curr_washout),'.','MarkerSize',20);
end
xlabel('Muscimol');
ylabel('Washout');
axis square;
line(ylim,ylim,'color','k');








