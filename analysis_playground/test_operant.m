%% Test analysis for imaging with operant task (over learning)
% (note: this includes things brought from test_corticostriatal and
% test_learning)

%% ~~~~~~~~~ SINGLE-SESSION ~~~~~~~~~

%% Load example session

animal = 'AP100';
% day = '2021-05-03';
day = '2021-05-07';

% animal = 'AP104';
% day = '2021-06-06';
% day = '2021-06-08';
% day = '2021-06-09';
% day = '2021-06-11';
% day = '2021-06-13';

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
AP_cellraster({stimOn_times,wheel_starts(wheel_move_response_idx),signals_events.responseTimes(1:n_trials)},rxn_sort_idx)

% PSTH with rewarded & ITI move starts
[~,rxn_sort_idx] = sort(stim_to_move);

AP_cellraster( ...
    {stimOn_times,wheel_starts(wheel_move_response_idx),wheel_starts(wheel_move_iti_idx),reward_t_timeline}, ...
    {rxn_sort_idx,rxn_sort_idx,1:length(wheel_move_iti_idx),1:length(reward_t_timeline)})


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


%% Stim response vs alternate timing stim onset (Kenneth alt)
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


%% Stim response vs alternate timing stim onset (only if same movement) (Kenneth alt 2)

% For each trial:
% - get all combinations of ITIs/quiescence periods that would have lead to
% the first movement being the same
% - this gives a slightly jittered stim on time

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
        t <= signals_events.responseTimes(curr_trial);
    
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
    
    % (sanity plot)
    t_plot_scale = 0.1;
    figure; hold on;
    plot(t(curr_trial_t_idx),wheel_velocity(curr_trial_t_idx),'k')
    plot(t(curr_trial_t_idx),t_from_move_trace(curr_trial_t_idx)*t_plot_scale,'r');
    plot(t(curr_trial_t_idx),t_from_quiescence_reset_trial*t_plot_scale,'b');
    plot(t(curr_trial_t_idx),[0;diff(wheel_position(curr_trial_t_idx))]*0.1,'g')
    line(repmat(curr_trial_t(1)+signals_events.trialITIValues(curr_trial-1),2,1),ylim);
    line(xlim,repmat(signals_events.trialQuiescenceValues(curr_trial),2,1)*t_plot_scale,'color','m');
    legend({'Wheel velocity','t from move','t from quiescence reset','wheel click','trial iti','trial quiescence'});
    
    % Find alternate stim times which would have given same first move
    param_timestep = 0.1; % (hardcoded in expDef)
    possible_iti = max([block.paramsValues.itiMin]):param_timestep:max([block.paramsValues.itiMax]);
    possible_quiescence = max([block.paramsValues.quiescenceMin]):param_timestep:max([block.paramsValues.quiescenceMax]);
    
    % (getting possible iti + quiescence = alternate stim times)
    alt_iti_reached = ((t(curr_trial_t_idx) - curr_trial_t(1)) > possible_iti);
    alt_quiescence_reached = (t_from_quiescence_reset_trial > possible_quiescence);
    [alt_stim_value,alt_stim_idx] = max(reshape(alt_iti_reached & ...
        permute(alt_quiescence_reached,[1,3,2]),length(curr_trial_t),[]),[],1);
    alt_stimOn_times_all = curr_trial_t(alt_stim_idx(alt_stim_value & alt_stim_idx));
    
    % (get alt stim times that would have happened within same quiescence)
    t_prestim_move = curr_trial_t(find(wheel_move(curr_trial_t_idx),1,'last'));
    alt_stimOn_times{curr_trial} = ...
        unique(alt_stimOn_times_all(alt_stimOn_times_all > t_prestim_move));
    
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

% wf_roi_fn = [wf_roi_path filesep 'wf_roi'];
% save(wf_roi_fn,'wf_roi');
% disp('Saved new widefield ROIs');

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

animals = {'AP100','AP101','AP103','AP104','AP105','AP106','AP107','AP108','AP109','AP111','AP112'};

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

animals = {'AP100','AP101','AP103','AP104','AP105','AP106','AP107','AP108','AP109','AP111','AP112'};

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

animals = {'AP100','AP105','AP106','AP107','AP108'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);
    
    % Get days with V1 muscimol and washout
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    
    v1_muscimol_idx = find(strcmpi(muscimol(muscimol_animal_idx).area,'v1'));
    washout_idx = find(strcmpi(muscimol(muscimol_animal_idx).area,'washout'));
    v1_washout_idx = washout_idx(find(washout_idx > v1_muscimol_idx,1,'first'));
    
    % Set experiments to use (V1 muscimol and subsequent washout)
    use_muscimol_experiments = ismember({experiments.day}, ...
        muscimol(muscimol_animal_idx).day([v1_muscimol_idx,v1_washout_idx]));
    experiments = experiments(use_muscimol_experiments);
    
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

animals = {'AP100','AP105','AP106','AP107','AP108'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_stimWheelRight';
    experiments = AP_find_experiments(animal,protocol);
    
    % Get days with V1 muscimol and washout
    data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
    muscimol_fn = [data_path filesep 'muscimol.mat'];
    load(muscimol_fn);
    muscimol_animal_idx = ismember({muscimol.animal},animal);
    
    v1_muscimol_idx = find(strcmpi(muscimol(muscimol_animal_idx).area,'v1'));
    washout_idx = find(strcmpi(muscimol(muscimol_animal_idx).area,'washout'));
    v1_washout_idx = washout_idx(find(washout_idx > v1_muscimol_idx,1,'first'));
    
    % Set experiments to use (V1 muscimol and subsequent washout)
    use_muscimol_experiments = ismember({experiments.day}, ...
        muscimol(muscimol_animal_idx).day([v1_muscimol_idx,v1_washout_idx]));
    experiments = experiments(use_muscimol_experiments);
    
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



%% Passive - ephys

clear all
disp('Passive trial activity (ephys)')

animals = {'AP100','AP101','AP104','AP105','AP106'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);
    
    % Set experiments to use (ephys)
    experiments = experiments([experiments.ephys]);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        preload_vars = who;
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        % (load LFP to find cortex start)
        lfp_channel = 'all';
        ephys_align = 'cortex';
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
        trial_data_all.mua_depth_edges = depth_group_edges;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
    
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = ['trial_activity_passive_ephys'];
save([save_path filesep save_fn],'-v7.3');
disp(['Saved: ' save_path filesep save_fn])


%% Task - ephys

clear all
disp('Task trial activity (ephys)')

animals = {'AP100','AP101','AP104','AP105','AP106'};

% Initialize save variable
trial_data_all = struct;

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
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
        % (load LFP to find cortex start)
        lfp_channel = 'all';
        ephys_align = 'cortex';
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
        trial_data_all.mua_depth_edges = depth_group_edges;
        
        AP_print_progress_fraction(curr_day,length(experiments));
        
        % Clear for next loop
        clearvars('-except',preload_vars{:});
        
    end
end

clearvars -except trial_data_all
disp('Finished loading all')

% Save
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
save_fn = ['trial_activity_task_ephys'];
save([save_path filesep save_fn],'-v7.3');

%% Passive reversal (tetO)

clear all
disp('Passive trial activity reversal (tetO)')

animals = {'AP107'};

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
        muscimol_experiments = ismember({experiments.day}, ...
            [muscimol(muscimol_animal_idx).day]);
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
save_fn = ['trial_activity_passive_teto_reversal'];
save([save_path filesep save_fn],'-v7.3');
disp(['Saved: ' save_path filesep save_fn])


%% Facecam troubleshooting

animals = {'AP100','AP101','AP103','AP104','AP105','AP106','AP107','AP108','AP109','AP111','AP112'};

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    protocol = 'AP_lcrGratingPassive';
%     protocol = 'AP_stimWheel';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    
    disp(['Loading ' animal]);
    
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment = experiments(curr_day).experiment(end);
        
        % Load experiment
        load_parts.cam = true;
        AP_load_experiment;
        
        
        %%%% Redo facecam strobes
        
        % FACECAM
        [facecam_dir,facecam_exists] = AP_cortexlab_filename(animal,day,experiment,'facecam');
        
        if facecam_exists
            % Get camera times
            facecam_fn = AP_cortexlab_filename(animal,day,experiment,'facecam');
            facecam_dir = fileparts(facecam_fn);
            facecam_t_savefile = [facecam_dir filesep 'facecam_t.mat'];
            
            % Get facecam strobes
            faceCamStrobe_idx = strcmp({Timeline.hw.inputs.name}, 'faceCamStrobe');
            faceCamStrobe_thresh = max(Timeline.rawDAQData(:,faceCamStrobe_idx))/5;
            faceCamStrobe = Timeline.rawDAQData(:,faceCamStrobe_idx) > faceCamStrobe_thresh;
            faceCamStrobe_up = find((~faceCamStrobe(1:end-1) & faceCamStrobe(2:end)))+1;
            faceCamStrobe_up_t = Timeline.rawDAQTimestamps(faceCamStrobe_up);
            
            % Get sync times for cameras
            [facecam_sync_frames,n_facecam_frames] = AP_get_cam_sync_frames(facecam_fn);            
            
            % Get the closest cam strobe to sync start, find offset and frame idx
            [~,facecam_strobe_sync] = min(abs(camSync_flip(1) - faceCamStrobe_up));
            facecam_frame_offset = facecam_sync_frames(1) - facecam_strobe_sync;
            facecam_frame_idx = [1:length(faceCamStrobe_up)] + facecam_frame_offset;
            
            % Check that the number of frames between syncs matches
            % video and timeline
            n_facecam_frames_syncd_movie = diff(facecam_sync_frames) + 1;
            [~,facecam_strobe_sync_end] = min(abs(camSync_flip(3) - faceCamStrobe_up));
            n_facecam_frames_syncd_timeline = facecam_strobe_sync_end - facecam_strobe_sync;
            if abs(n_facecam_frames_syncd_movie - n_facecam_frames_syncd_timeline) > 2
                warning('[%s %s %d] Facecam: different n frames video vs timeline',animal,day,experiment);
                [facecam_sync_frames,n_facecam_frames] = AP_get_cam_sync_frames(facecam_fn);
            else
                close(gcf)
            end
            
            % Get times of cam frames in timeline
            facecam_t = nan(n_facecam_frames,1);
            facecam_t(facecam_frame_idx(facecam_frame_idx > 0)) = faceCamStrobe_up_t(facecam_frame_idx > 0);
            
            try
                save(facecam_t_savefile,'facecam_t');
            catch me
                warning (['Can''t save: ' facecam_t_savefile])
            end
        end
        
    end
end


%% ~~~~~~~~~ BATCH ANALYSIS ~~~~~~~~~

%% >> Passive

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_passive_tetO';
% data_fn = 'trial_activity_passive_tetO_reversal';
% data_fn = 'trial_activity_passive_cstr';

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
use_animal = 1;
plot_stim = 1;
curr_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    stim_v_avg{use_animal}(:,:,:,plot_stim));
AP_image_scroll(curr_px,t);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
axis image;
set(gcf,'name',animals{use_animal});


%% Time-average ROI activity by day

plot_roi = 7;

curr_act = fluor_roi_deconv(:,:,plot_roi);

use_trials = trial_stim_allcat == 1 & quiescent_trials;

use_t = t >= 0.1 & t <= 0.2;
curr_act_timeavg = nanmean(double(curr_act(:,use_t)),2);

roi_animal_day_avg = accumarray([trial_animal(use_trials),trial_day(use_trials)], ...
    curr_act_timeavg(use_trials),[max(trial_animal),max(trial_day)], ...
    @nanmean,NaN);

figure; hold on;
passive_x = [1:max(trial_day)]-3;
plot(passive_x,roi_animal_day_avg','color',[0.5,0.5,0.5]);
errorbar(passive_x,nanmean(roi_animal_day_avg,1),AP_sem(roi_animal_day_avg,1),'k','linewidth',2);
xlabel('Day');
ylabel(wf_roi(plot_roi).area);

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
roi_animal_day_avg_long = reshape(permute( ...
    padarray(roi_animal_daysplit_avg,[0,0,1],NaN,'post'),[3,2,1]),[],length(animals));
% plot([1:size(roi_animal_day_avg_long,1)]/(n_daysplit+1),roi_animal_day_avg_long);

figure; hold on;
daysplit_x = [1:size(roi_animal_day_avg_long,1)]/(n_daysplit+1) - 3;
plot(daysplit_x,roi_animal_day_avg_long','color',[0.5,0.5,0.5],'linewidth',2);
errorbar(daysplit_x,nanmean(roi_animal_day_avg_long,2),AP_sem(roi_animal_day_avg_long,2),'k','linewidth',2);
xlabel('Day');
ylabel(wf_roi(plot_roi).area);


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

% Rxn fraction
rxn_frac_prepadcat = cell2mat(cellfun(@(x,use_days) padarray( ...
    cellfun(@(x) nanmean(x > 0.1 & x < 0.25),x(use_days)),[3,0],NaN,'pre'), ...
    {bhv.stim_move_t},use_days,'uni',false)');

% Group stim responses by behavioral performance bins
stim_v_avg_cat = cell2mat(permute(stim_v_avg,[1,3,2,4]));

bhv_padcat_grp = ...
    discretize(rxn_frac_prepadcat,[-Inf,0,0.5,Inf]);
bhv_padcat_grp(isnan(bhv_padcat_grp)) = 0;

stim_v_avg_cat_grp = reshape(grpstats( ...
    reshape(permute(stim_v_avg_cat,[1,2,4,3]),[],size(stim_v_avg_cat,3))', ...
    bhv_padcat_grp)',n_vs,length(t),length(stim_unique),[]);

curr_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),stim_v_avg_cat_grp);
AP_image_scroll(reshape(permute(curr_px,[1,2,4,3,5]), ...
    size(curr_px,1),[],size(curr_px,3),size(curr_px,5)),t);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5],[], ...
    [size(U_master,1),size(U_master,2),1,size(curr_px,4)]);

% % (make movie)
% movie_stim = 3;
% movie_px = squeeze(AP_svdFrameReconstruct(U_master(:,:,1:n_vs),stim_v_avg_cat_grp(:,:,movie_stim,:)));
% movie_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\passive_wf.avi';
% movie_framerate = 1/mean(diff(t));
% movie_cmap = colormap(brewermap([],'PrGn'));
% movie_caxis = prctile(movie_px(:),99.9).*[-1,1];
% movie_position = [48,500,800,270];
% movie_annotation = {'Naive','Pre-learned','Post-learned'};
% AP_movie2avi(movie_px,movie_framerate,movie_cmap,movie_caxis,movie_position,movie_fn,t,movie_annotation);

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
use_stim = 3;
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
use_stim = 3;

curr_act_timeavg = cellfun(@(x) squeeze(nanmean(x(curr_roi,use_t,:,use_stim),2)),stim_roi_avg,'uni',false)';

max_days = max(cellfun(@length,curr_act_timeavg));

curr_act_timeavg_pad = cell2mat(cellfun(@(x) ...
    padarray(x,[max_days-length(x),0],NaN,'post'), ...
    curr_act_timeavg,'uni',false)');

rxn_frac_pad = cell2mat(cellfun(@(x) padarray(cat(1,nan(3,1),cellfun(@(x) nanmean(x > 0.1 & x < 0.2),x)), ...
    [max_days-length(x)-3,0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_move_t},use_days,'uni',false),'uni',false));

alt_rxn_frac_pad = cell2mat(cellfun(@(x) padarray(cat(1,nan(3,1),cellfun(@(x) ...
    nanmean(cell2mat(x) > 0.1 & cell2mat(x) < 0.25),x)), ...
    [max_days-length(x)-3,0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.alt_stim_move_t},use_days,'uni',false),'uni',false));

rxn_frac_pad_diff = rxn_frac_pad - alt_rxn_frac_pad;

% (grab learned day from behavior script)
learned_day_x = [1:max_days]'-[learned_day+3];

[rxn_grp_mean,learned_day_grp] = ...
    grpstats(rxn_frac_pad(:),learned_day_x(:),{'nanmean','gname'});
rxndiff_grp_mean = ...
    grpstats(rxn_frac_pad_diff(:),learned_day_x(:),{'nanmean'});
[act_grp_mean,act_grp_sem,act_grp_n] = grpstats(curr_act_timeavg_pad(:),learned_day_x(:),{'nanmean','sem','numel'});
learned_day_grp = cellfun(@str2num,learned_day_grp);

plot_days = act_grp_n > 4;

figure;
subplot(1,3,1); hold on;
plot(learned_day_x,rxn_frac_pad);
plot(learned_day_grp(plot_days),rxn_grp_mean(plot_days),'k','linewidth',2);
xlabel('Learned day');
ylabel('Rxn frac');
line([0,0],ylim,'color','k','linestyle','--');

subplot(1,3,2); hold on;
plot(learned_day_x,rxn_frac_pad_diff);
plot(learned_day_grp(plot_days),rxndiff_grp_mean(plot_days),'k','linewidth',2);
xlabel('Learned day');
ylabel('Rxn frac diff');
line([0,0],ylim,'color','k','linestyle','--');

subplot(1,3,3); hold on;
plot(learned_day_x,curr_act_timeavg_pad);
errorbar(learned_day_grp(plot_days),act_grp_mean(plot_days),act_grp_sem(plot_days),'k','linewidth',2);
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




% Plot select ROIs with contra stim
plot_rois = [1,3,7,8];
figure;
for curr_roi_idx = 1:length(plot_rois)
    curr_roi = plot_rois(curr_roi_idx);
    curr_contra_roi = curr_roi + size(wf_roi,1);
    
    curr_roi_prelearned = cell2mat(cellfun(@(x,learned) ...
        nanmean(x(curr_roi,:,learned < 0,stim_unique == 1),3), ...
        stim_roi_avg,learned_days_animal,'uni',false)');
    curr_contra_roi_prelearned = cell2mat(cellfun(@(x,learned) ...
        nanmean(x(curr_contra_roi,:,learned < 0,stim_unique == -1),3), ...
        stim_roi_avg,learned_days_animal,'uni',false)');
    
    curr_roi_learned = cell2mat(cellfun(@(x,learned) ...
        nanmean(x(curr_roi,:,learned > 0,stim_unique == 1),3), ...
        stim_roi_avg,learned_days_animal,'uni',false)');
    curr_contra_roi_learned = cell2mat(cellfun(@(x,learned) ...
        nanmean(x(curr_contra_roi,:,learned > 0,stim_unique == -1),3), ...
        stim_roi_avg,learned_days_animal,'uni',false)');
    
    p1 = subplot(length(plot_rois),2,curr_roi_idx*2-1);
    title(wf_roi(curr_roi).area);
    AP_errorfill(t,nanmean(curr_roi_prelearned,1),std(curr_roi_prelearned,[],1),'k');
    AP_errorfill(t,nanmean(curr_roi_learned,1),std(curr_roi_learned,[],1),'r');
    line([0,0],ylim,'color','k','linestyle','--');
    line([0.5,0.5],ylim,'color','k','linestyle','--');
    
    p2 = subplot(length(plot_rois),2,curr_roi_idx*2);
    title(wf_roi(curr_contra_roi).area);
    h1 = AP_errorfill(t,nanmean(curr_contra_roi_prelearned,1),std(curr_contra_roi_prelearned,[],1),'k');
    h2 = AP_errorfill(t,nanmean(curr_contra_roi_learned,1),std(curr_contra_roi_learned,[],1),'r');
    line([0,0],ylim,'color','k','linestyle','--');
    line([0.5,0.5],ylim,'color','k','linestyle','--');
    
    linkaxes([p1,p2],'xy')
    
end
legend([h1,h2],{'Pre-learned','Learned'});




%% >> (Choiceworld task: post-learn only)

trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\wf_ephys_choiceworld\paper\data\paper_unused';
data_fn = 'trial_activity_choiceworld_wfonly';

AP_load_trials_operant;


%% >> Task

% Task
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_task_teto';
% data_fn = 'trial_activity_task_corticostriatal';

AP_load_trials_operant;

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

use_rxn = [0.1,0.25];

bhv_bins = [0,0.5,Inf];
bhv_discretize = discretize(rxn_frac_allcat,bhv_bins);

act_avg_animal_v = nan(n_vs,length(t),max(bhv_discretize),length(animals));
for curr_animal = 1:length(animals)
    for curr_bhv = 1:max(bhv_discretize)
        
        curr_trials = trial_outcome_allcat == 1 & ...
            curr_animal == trial_animal & ...
            move_t >= use_rxn(1) & move_t <= use_rxn(2) & ...
            bhv_discretize == curr_bhv;
        
        % (raw)
        %         act_avg_animal_v(:,:,curr_bhv,curr_animal) = ...
        %             permute(nanmean(fluor_allcat_deconv(curr_trials,:,:),1),[3,2,1]);
        
        % (reduced)
        act_avg_animal_v(:,:,curr_bhv,curr_animal) = ...
            permute(nanmean(fluor_allcat_deconv(curr_trials,:,:) - ...
            fluor_taskpred_reduced_allcat(curr_trials,:,:,1),1),[3,2,1]);
        
    end
end

act_avg_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),nanmean(act_avg_animal_v,4));

AP_image_scroll(act_avg_px,t);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);

% % (make movie)
% movie_px = act_avg_px;
% movie_fn = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\task_reduced_wf.avi';
% movie_framerate = 1/mean(diff(t));
% movie_cmap = colormap(brewermap([],'PrGn'));
% movie_caxis = prctile(movie_px(:),99.9).*[-1,1];
% movie_position = [48,500,520,270];
% movie_annotation = {'Pre-learned','Post-learned'};
% AP_movie2avi(movie_px,movie_framerate,movie_cmap,movie_caxis,movie_position,movie_fn,t,movie_annotation);




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

plot_roi = 7;

use_trials = trial_outcome_allcat == 1 & rxn_frac_allcat > 0.4;
% use_trials = trial_outcome_allcat == 1 & stim_response_idx_allcat > 0.2;

curr_act = fluor_roi_deconv(use_trials,:,plot_roi);
curr_act_taskpred = fluor_roi_taskpred(use_trials,:,plot_roi);
curr_act_reduced = fluor_roi_deconv(use_trials,:,plot_roi) - fluor_roi_taskpred_reduced(use_trials,:,plot_roi,1);

% Get activity by reaction time bin
rxn_bins = [0:0.05:0.4]';
% rxn_bins = [-inf,inf];
rxn_bin_centers = rxn_bins(1:end-1) + diff(rxn_bins)./2;
move_t_discretize = discretize(move_t(use_trials),rxn_bins);

curr_act_rxn = accumarray( ...
    {reshape(repmat(1:length(t),sum(~isnan(move_t_discretize)),1),[],1), ...
    reshape(repmat(move_t_discretize(~isnan(move_t_discretize)),1,length(t)),[],1)}, ...
    reshape(curr_act(~isnan(move_t_discretize),:),[],1), ...
    [length(t),length(rxn_bin_centers)],@nanmean,single(NaN));

curr_act_taskpred_rxn = accumarray( ...
    {reshape(repmat(1:length(t),sum(~isnan(move_t_discretize)),1),[],1), ...
    reshape(repmat(move_t_discretize(~isnan(move_t_discretize)),1,length(t)),[],1)}, ...
    reshape(curr_act_taskpred(~isnan(move_t_discretize),:),[],1), ...
    [length(t),length(rxn_bin_centers)],@nanmean,single(NaN));

curr_act_reduced_rxn = accumarray( ...
    {reshape(repmat(1:length(t),sum(~isnan(move_t_discretize)),1),[],1), ...
    reshape(repmat(move_t_discretize(~isnan(move_t_discretize)),1,length(t)),[],1)}, ...
    reshape(curr_act_reduced(~isnan(move_t_discretize),:),[],1), ...
    [length(t),length(rxn_bin_centers)],@nanmean,single(NaN));

% Plot activity as lines
figure; hold on;
rxn_col = copper(length(rxn_bin_centers));
set(gca,'ColorOrder',rxn_col);
plot(t,curr_act_rxn','linewidth',2)
arrayfun(@(x) line(repmat(rxn_bin_centers(x),2,1), ...
    ylim,'color',rxn_col(x,:)),1:length(rxn_bin_centers));
line([0,0],ylim,'color','k','linestyle','--');
ylabel(wf_roi(plot_roi).area);

% Stackplot activity/task-pred activity
figure; hold on;
line_spread = 2*nanstd(curr_act_rxn(:));
h1 = AP_stackplot(curr_act_rxn,t,line_spread,[],rxn_col);
h2 = AP_stackplot(curr_act_taskpred_rxn,t,line_spread,[],'b');
set(h2,'linewidth',1);
h3 = AP_stackplot(curr_act_reduced_rxn,t,line_spread,[],'r');
set(h3,'linewidth',1);
arrayfun(@(x) line(repmat(rxn_bin_centers(x),2,1), ...
    ylim,'color',rxn_col(x,:)),1:length(rxn_bin_centers));
line([0,0],ylim,'color','k','linestyle','--');
xlabel('Time');
ylabel(wf_roi(plot_roi).area);
legend([h1(1),h2(1),h3(1)],'Measured','Task-predicted','Stim (reduced)');


%% Time-average ROI activity by day/behavior

plot_roi = 7;

% curr_act = fluor_roi_deconv(:,:,plot_roi);
curr_act = fluor_roi_deconv(:,:,plot_roi) - fluor_roi_taskpred_reduced(:,:,plot_roi,1);

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
daysplit_x = [1:size(roi_animal_daysplit_avg_long,1)]/(n_daysplit+1);
plot(daysplit_x,roi_animal_daysplit_avg_long,'color',[0.5,0.5,0.5]);
errorbar(daysplit_x,nanmean(roi_animal_daysplit_avg_long,2), ...
    AP_sem(roi_animal_daysplit_avg_long,2),'k','linewidth',2);
xlabel('Day');
ylabel(wf_roi(plot_roi).area);

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

% Get learned day
learned_day = cellfun(@(x) find(x,1),{bhv.learned_days})';

n_days_animal = accumarray(trial_animal,trial_day,[],@max);
learned_day_animal_x = cellfun(@(ld,n) [1:n]'-(ld), ...
    num2cell(learned_day),num2cell(n_days_animal),'uni',false);

% Get task>cortex parameters
n_regressors = length(task_regressor_labels);
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

% Average kernels by learning stage
learn_stages = cellfun(@(x) discretize(x,[-Inf,0,Inf]),learned_day_animal_x,'uni',false);
 
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

% Plot max kernels
for curr_reg = 1:n_regressors
    n_subreg = size(fluor_taskpred_px_stage_avg{curr_reg},4);
    figure;
    h = tiledlayout(n_subreg,n_learn_stages, ...
        'TileSpacing','compact','padding','compact');
    c = max(fluor_taskpred_px_stage_avg{curr_reg}(:)).*[-1,1];
    for curr_subreg = 1:n_subreg
        for curr_stage = 1:n_learn_stages
            nexttile;
            imagesc(max(fluor_taskpred_px_stage_avg{curr_reg} ...
                (:,:,:,curr_subreg,curr_stage),[],3));
            title(sprintf('K %d, Stage %d',curr_subreg,curr_stage));
            axis image off;
            caxis(c);
            colormap(brewermap([],'PrGn'));
            AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
        end
    end
    title(h,task_regressor_labels{curr_reg});
end



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

curr_roi = 6;

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

rxn_frac_pad = cell2mat(cellfun(@(x) padarray(cellfun(@(x) nanmean(x > 0.1 & x < 0.25),x), ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.stim_move_t},use_days,'uni',false),'uni',false));

alt_rxn_frac_pad = cell2mat(cellfun(@(x) padarray(cellfun(@(x) ...
    nanmean(cell2mat(x) > 0.1 & cell2mat(x) < 0.25),x), ...
    [max_days-length(x),0],NaN,'post'), ...
    cellfun(@(x,use_days) x(use_days),{bhv.alt_stim_move_t},use_days,'uni',false),'uni',false));

rxn_frac_pad_diff = rxn_frac_pad - alt_rxn_frac_pad;

% (grab learned day from behavior script)
learned_day_x = [1:max_days]'-learned_day;

[rxn_grp_mean,learned_day_grp] = ...
    grpstats(rxn_frac_pad(:),learned_day_x(:),{'nanmean','gname'});
rxndiff_grp_mean = ...
    grpstats(rxn_frac_pad_diff(:),learned_day_x(:),{'nanmean'});
[act_grp_mean,act_grp_sem,act_grp_n] = grpstats(curr_act_timeavg_pad(:),learned_day_x(:),{'nanmean','sem','numel'});
learned_day_grp = cellfun(@str2num,learned_day_grp);

plot_days = act_grp_n > 4;

figure;
subplot(1,3,1); hold on;
plot(learned_day_x,rxn_frac_pad);
plot(learned_day_grp(plot_days),rxn_grp_mean(plot_days),'k','linewidth',2);
xlabel('Learned day');
ylabel('Rxn frac');
line([0,0],ylim,'color','k','linestyle','--');

subplot(1,3,2); hold on;
plot(learned_day_x,rxn_frac_pad_diff);
plot(learned_day_grp(plot_days),rxndiff_grp_mean(plot_days),'k','linewidth',2);
xlabel('Learned day');
ylabel('Rxn frac diff');
line([0,0],ylim,'color','k','linestyle','--');

subplot(1,3,3); hold on;
plot(learned_day_x,curr_act_timeavg_pad);
errorbar(learned_day_grp(plot_days),act_grp_mean(plot_days),act_grp_sem(plot_days),'k','linewidth',2);
xlabel('Learned day');
ylabel(sprintf('%s fluorescence',wf_roi(curr_roi).area))
line([0,0],ylim,'color','k','linestyle','--');

linkaxes(get(gcf,'Children'),'x')

figure; hold on;
animal_colors = max(brewermap(length(animals),'Set3')-0.2,0);
set(gca,'ColorOrder',animal_colors);
plot(rxn_frac_pad,curr_act_timeavg_pad,'MarkerSize',15);
plot(rxn_frac_pad,curr_act_timeavg_pad,'.','MarkerSize',15);
legend(animals,'location','nw')
xlabel('Rxn frac');
ylabel('FRm fluor');

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

[roi_daysplit,roi_daysplit_sem,roi_daysplit_grp,roi_daysplit_n] = ...
    grpstats(roi_animal_daysplit_avg_long(:), ...
    learned_day_x_daysplit(:),{'nanmean','sem','gname','numel'});
roi_daysplit_grp = cellfun(@str2num,roi_daysplit_grp);

plot_days_split = isnan(roi_daysplit) | roi_daysplit_n > 4;

figure; hold on;
plot(learned_day_x_daysplit,roi_animal_daysplit_avg_long,'color',[0.5,0.5,0.5]);
% plot(roi_daysplit_grp(plot_days_split),roi_daysplit(plot_days_split),'k','linewidth',2);
errorbar(roi_daysplit_grp(plot_days_split),roi_daysplit(plot_days_split), ...
    roi_daysplit_sem(plot_days_split),'k','linewidth',2);
line([0,0],ylim,'color','k','linestyle','--');
xlabel('Learned day');
ylabel(sprintf('%s fluorescence',wf_roi(curr_roi).area))

%% Movement activity (vs learning, vs stim/non-stim)

learned_day = cellfun(@(x) find(x,1),{bhv.learned_days})';
trial_learned_day = trial_day - learned_day(trial_animal);

% Plot average reduced stimulus pixel activity pre/post learning
stim_regressor = strcmp(task_regressor_labels,'Stim');
stim_v_act = fluor_allcat_deconv - fluor_taskpred_reduced_allcat(:,:,:,stim_regressor);

[learned_idx,t_idx,v_idx] = ...
    ndgrid((trial_learned_day >= 0)+1,1:length(t),1:n_vs);
[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_vs);

n_stages = max(learned_idx(:));
stim_v_act_learn_avg = permute(accumarray( ...
    [learned_idx(:),t_idx(:),v_idx(:),animal_idx(:)], ...
    stim_v_act(:), ...
    [n_stages,length(t),n_vs,length(animals)], ...
    @nanmean,NaN('single')),[3,2,1,4]);

stim_px_act_learn_avg = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    squeeze(nanmean(stim_v_act_learn_avg,4)));

figure;
use_t = t >= 0.05 & t <= 0.2;
tiledlayout(1,n_stages,'TileSpacing','tight','padding','compact');
c = (max(stim_px_act_learn_avg(:)).*[-1,1])*0.5;
for curr_stage = 1:n_stages
    nexttile;
    imagesc(nanmean(stim_px_act_learn_avg(:,:,use_t,curr_stage),3));
    caxis(c);
    axis image off;
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    colormap(brewermap([],'PrGn'));
    title(sprintf('Learned stage %d',curr_stage));
end

%% Move-align activity (and reduce)

% (minimum number n to plot)
min_n = 4;

% Align trial activity to movement
t_leeway = -t(1);
leeway_samples = round(t_leeway*(sample_rate));

fluor_roi_deconv_move = nan(size(fluor_roi_deconv),'single');
move_v_act = nan(size(fluor_allcat_deconv),'single');
move_roi_act = nan(size(fluor_roi_deconv),'single');

for curr_trial = find(~isnan(move_idx))'
    curr_shift_frames = (move_idx(curr_trial)-leeway_samples) + [0:length(t)-1];
    curr_shift_frames_use = curr_shift_frames > 0 & curr_shift_frames < length(t);
    
    curr_grab_frames = curr_shift_frames(curr_shift_frames_use);
    curr_fill_frames = find(curr_shift_frames_use,length(t));
   
    % (raw activity)
   fluor_roi_deconv_move(curr_trial,curr_fill_frames,:) = ...
       fluor_roi_deconv(curr_trial,curr_grab_frames,:);
   
   % (reduced movement activity)
   move_regressor = strcmp(task_regressor_labels,'Move onset');
   
   move_v_act(i,curr_fill_frames,:) = ...
       fluor_allcat_deconv(i,curr_grab_frames,:) - ...
       fluor_taskpred_reduced_allcat(i,curr_grab_frames,:,move_regressor);
   
   move_roi_act(i,curr_fill_frames,:) = ...
       fluor_roi_deconv(i,curr_grab_frames,:) - ...
       fluor_roi_taskpred_reduced(i,curr_grab_frames,:,move_regressor);
end

% Plot average pixel movement activity by learning
[learned_idx,t_idx,v_idx] = ...
    ndgrid((trial_learned_day >= 0)+1,1:length(t),1:n_vs);
[animal_idx,~,~] = ...
    ndgrid(trial_animal,1:length(t),1:n_vs);

n_stages = max(learned_idx(:));
move_v_act_learn_avg = permute(accumarray( ...
    [learned_idx(:),t_idx(:),v_idx(:),animal_idx(:)], ...
    move_v_act(:), ...
    [n_stages,length(t),n_vs,length(animals)], ...
    @nanmean,NaN('single')),[3,2,1,4]);

move_px_act_learn_avg = AP_svdFrameReconstruct(U_master(:,:,1:n_vs), ...
    squeeze(nanmean(move_v_act_learn_avg,4)));

figure;
use_t = t >= -0.1 & t <= 0.1;
tiledlayout(1,n_stages,'TileSpacing','compact','padding','compact');
c = (max(move_px_act_learn_avg(:)).*[-1,1])*0.5;
for curr_stage = 1:n_stages
    nexttile;
    imagesc(nanmean(move_px_act_learn_avg(:,:,use_t,curr_stage),3));
    caxis(c);
    axis image off;
    AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
    colormap(brewermap([],'PrGn'));
    title(sprintf('Learned stage %d',curr_stage));
end

% Plot ROIs by learning
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
move_roi_act_learn_avg = accumarray( ...
    [learned_day_minidx(:),daysplit_idx(:),t_idx(:),roi_idx(:),animal_idx(:)], ...
    move_roi_act(:), ...
    [max(learned_day_minidx(:)),n_daysplit,length(t),n_rois,length(animals)], ...
    @nanmean,NaN('single'));

% (tmax activity: learned day x daysplit x roi x animal)
use_t = t >= 0 & t <= 0.5;
move_roi_act_tmax_daysplit = ...
    squeeze(max(move_roi_act_learn_avg(:,:,use_t,:,:),[],3));

% Plot time-max by learned day
learned_day_x_range = minmax(learned_day_idx);
learned_day_x = [learned_day_x_range(1):learned_day_x_range(2)]';
learned_daysplit_x = learned_day_x + linspace(0,1,n_daysplit+1);

plot_learned_day = sum(~isnan(move_roi_act_tmax_daysplit(:,1,1,:)),4) > min_n;

figure;
plot_rois = [1:6];
tiledlayout(length(plot_rois),1,'TileSpacing','compact','padding','compact');
for curr_roi_idx = 1:length(plot_rois)
    curr_l_roi = plot_rois(curr_roi_idx);

    nexttile;
   
    errorbar(reshape(learned_daysplit_x(plot_learned_day,:)',[],1), ...
        reshape(padarray(nanmean( ...
        move_roi_act_tmax_daysplit(plot_learned_day,:,curr_l_roi,:),4), ...
        [0,1],NaN,'post')',[],1), ...
        reshape(padarray(AP_sem( ...
        move_roi_act_tmax_daysplit(plot_learned_day,:,curr_l_roi,:),4), ...
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
move_roi_act_tmax_daysplit_normval = ...
    nanmean(nanmean(move_roi_act_tmax_daysplit(norm_days,:,:,:),2),1);
move_roi_act_tmax_daysplit_norm = ...
    (move_roi_act_tmax_daysplit-move_roi_act_tmax_daysplit_normval)./ ...
    move_roi_act_tmax_daysplit_normval;

figure;
plot_rois = [1,6];
errorbar( ...
    repmat(reshape(learned_daysplit_x(plot_learned_day,:)',[],1),1,length(plot_rois)), ...
    reshape(permute(padarray(nanmean( ...
    move_roi_act_tmax_daysplit_norm(plot_learned_day,:,plot_rois,:),4), ...
    [0,1],NaN,'post'),[2,1,3]),[],length(plot_rois)), ...
    reshape(permute(padarray(AP_sem( ...
    move_roi_act_tmax_daysplit_norm(plot_learned_day,:,plot_rois,:),4), ...
    [0,1],NaN,'post'),[2,1,3]),[],length(plot_rois)),'linewidth',2,'CapSize',0);
xlabel('Learned day');
ylabel('Fluorescence (normalized to pre-learn)');
xline(0,'linestyle','--');
legend({wf_roi(plot_rois).area},'location','nw')


%% Describing mPFC as x*somatomotor + y*visual

roi_idx_mpfc = find(strcmpi({wf_roi.area},'mpfc_l'));
roi_idx_v1 = find(strcmpi({wf_roi.area},'v1_l'));
roi_idx_sm = find(strcmpi({wf_roi.area},'sm_l'));

% (timecourse of selected time together)
ev_all = nan(length(animals),length(learned_day_unique),3);
for curr_animal = 1:length(animals)
    for curr_ld_idx = 1:length(learned_day_unique)
        
        curr_ld = learned_day_unique(curr_ld_idx);
        use_trials = trial_animal == curr_animal & trial_learned_day == curr_ld;
        
        if ~any(use_trials)
            continue
        end
        
        use_t = t > 0 & t < 0.5;
        a = reshape(permute(fluor_roi_deconv(use_trials,use_t,:),[2,1,3]),[],n_rois)';
        discontinuities = reshape([true(sum(use_trials),1), ...
            false(sum(use_trials),sum(use_t)-1)]',[],1);

        [k,predicted_signals,ev_vis,predicted_signals_reduced] = ...
            AP_regresskernel(a(roi_idx_v1,:), ...
            a(roi_idx_mpfc,:),[-20:20], ...
            0,[],5,false,true,discontinuities);
        
        [k,predicted_signals,ev_mot,predicted_signals_reduced] = ...
            AP_regresskernel(a(roi_idx_sm,:), ...
            a(roi_idx_mpfc,:),[-20:20], ...
            0,[],5,false,true,discontinuities);
                
        [k,predicted_signals,ev_vismot,predicted_signals_reduced] = ...
            AP_regresskernel(a([roi_idx_v1,roi_idx_sm],:), ...
            a(roi_idx_mpfc,:),[-20:20], ...
            0,[],5,false,true,discontinuities);
       
        
        ev_all(curr_animal,curr_ld_idx,1) = ev_vis.total;
        ev_all(curr_animal,curr_ld_idx,2) = ev_mot.total;
        ev_all(curr_animal,curr_ld_idx,3) = ev_vismot.total;
        
    end
    AP_print_progress_fraction(curr_animal,length(animals));
end

figure;
errorbar(repmat(learned_day_unique,3,1)', ...
    squeeze(nanmean(ev_all,1)),squeeze(AP_sem(ev_all,1)))
xline(0);
xlabel('Learned day');
ylabel('mPCF frac expl var')
legend({'Vis','Mot','Vis+Mot'});


% (by time chunk)
t_chunk_size = [-0.1,0.1]; % size of sliding time window
t_chunk_centers = [-0.4:0.1:0.5]; % centers of sliding time window
t_chunk = mat2cell(t_chunk_centers' + ...
    t_chunk_size,ones(length(t_chunk_centers),1),2);

ev_all = nan(length(learned_day_unique),length(t_chunk),length(animals));
for curr_t_chunk = 1:length(t_chunk)
    for curr_animal = 1:length(animals)
        for curr_ld_idx = 1:length(learned_day_unique)
            
            curr_ld = learned_day_unique(curr_ld_idx);
            use_trials = trial_animal == curr_animal & trial_learned_day == curr_ld;
            
            if ~any(use_trials)
                continue
            end
            
            use_t = t > t_chunk{curr_t_chunk}(1) & t < t_chunk{curr_t_chunk}(2);
            a = reshape(permute(fluor_roi_deconv(use_trials,use_t,:),[2,1,3]),[],n_rois)';
            discontinuities = reshape([true(sum(use_trials),1), ...
                false(sum(use_trials),sum(use_t)-1)]',[],1);
            
            [k,predicted_signals,ev_vis,predicted_signals_reduced] = ...
                AP_regresskernel(a(roi_idx_v1,:), ...
                a(roi_idx_mpfc,:),[-5:5], ...
                0,[],5,false,true,discontinuities);
            
            ev_all(curr_ld_idx,curr_t_chunk,curr_animal) = ev_vis.total;
        
        end
    end
    AP_print_progress_fraction(curr_t_chunk,length(t_chunk));
end

plot_learned_day = sum(~all(isnan(ev_all),2),3) > min_n;
figure;
imagesc(t_chunk_centers,learned_day_unique(plot_learned_day), ...
    nanmean(ev_all(plot_learned_day,:,:),3));
xlabel('Time from stim onset');
ylabel('Learned day');
colormap(brewermap([],'Greys'));
c = colorbar;
ylabel(c,'V1->mPFC frac explained var');
yline(0-0.5,'r')


%%% SANITY CHECK: predict V1 from SM/mPFC
% (timecourse of selected time together)
ev_all = nan(length(animals),length(learned_day_unique),3);
for curr_animal = 1:length(animals)
    for curr_ld_idx = 1:length(learned_day_unique)
        
        curr_ld = learned_day_unique(curr_ld_idx);
        use_trials = trial_animal == curr_animal & trial_learned_day == curr_ld;
        
        if ~any(use_trials)
            continue
        end
        
        use_t = t > 0 & t < 0.5;
        a = reshape(permute(fluor_roi_deconv(use_trials,use_t,:),[2,1,3]),[],n_rois)';
        discontinuities = reshape([true(sum(use_trials),1), ...
            false(sum(use_trials),sum(use_t)-1)]',[],1);

        [k,predicted_signals,ev_vis,predicted_signals_reduced] = ...
            AP_regresskernel(a(roi_idx_mpfc,:), ...
            a(roi_idx_v1,:),[-20:20], ...
            0,[],5,false,true,discontinuities);
        
        [k,predicted_signals,ev_mot,predicted_signals_reduced] = ...
            AP_regresskernel(a(roi_idx_sm,:), ...
            a(roi_idx_v1,:),[-20:20], ...
            0,[],5,false,true,discontinuities);
                
        [k,predicted_signals,ev_vismot,predicted_signals_reduced] = ...
            AP_regresskernel(a([roi_idx_mpfc,roi_idx_sm],:), ...
            a(roi_idx_v1,:),[-20:20], ...
            0,[],5,false,true,discontinuities);
       
        
        ev_all(curr_animal,curr_ld_idx,1) = ev_vis.total;
        ev_all(curr_animal,curr_ld_idx,2) = ev_mot.total;
        ev_all(curr_animal,curr_ld_idx,3) = ev_vismot.total;
        
    end
    AP_print_progress_fraction(curr_animal,length(animals));
end

figure;
errorbar(repmat(learned_day_unique,3,1)', ...
    squeeze(nanmean(ev_all,1)),squeeze(AP_sem(ev_all,1)))
xline(0);
xlabel('Learned day');
ylabel('V1 frac expl var')
legend({'mPFC','Mot','mPFC+Mot'});


%% Trial data

% Plot all trials (raw and reduced stim: by stage, sorted by reaction time)
plot_rois = [1,6,7];
trial_learned_stage = discretize(trial_learned_day,[-Inf,0,Inf]);
n_trial_smooth = 20;

figure;
h = tiledlayout(2,length(plot_rois),'TileSpacing','compact','padding','compact');
for curr_stage = 1:max(trial_learned_stage)
    for curr_roi = plot_rois
        curr_contra_roi = curr_roi + size(wf_roi,1);
        
        use_trials = find(trial_learned_stage == curr_stage);
        [~,sort_idx] = sort(move_t(use_trials));
        
%         % (just ROI)
%         curr_data_sort = fluor_roi_deconv(use_trials(sort_idx),:,curr_roi);
        
        % (subtract R from L)
        curr_data_sort = (fluor_roi_deconv(use_trials(sort_idx),:,curr_roi) - ...
            fluor_roi_deconv(use_trials(sort_idx),:,curr_contra_roi));

        curr_data_sort_smooth = convn(curr_data_sort, ...
            ones(n_trial_smooth,1)./n_trial_smooth,'same');
        
        nexttile;
        imagesc(t,[],curr_data_sort_smooth);hold on;
        colormap(brewermap([],'PrGn'));
        c = prctile(reshape(fluor_roi_deconv(:,:,curr_roi),[],1),90).*[-1,1];
        caxis(c);
        xline(0,'color','r');
        plot(move_t(use_trials(sort_idx)),1:length(use_trials),'color',[0.6,0,0.6]);
        xlabel('Time from stim (s)');
        ylabel('Trial (rxn-time sorted)');
        title(sprintf('%s, stage %d',wf_roi(curr_roi).area,curr_stage));
    end
end
title(h,'Contra-subtracted activity');


% Scale somatomotor activity to ROI, subtract and plot

% Get scaling factor from somatomotor to other ROIs (trial max)
roi_idx_sm = find(strcmpi({wf_roi.area},'sm_l'));
% sm_roi_scale = squeeze(max(fluor_roi_deconv,[],2)./ ...
%     max(fluor_roi_deconv(:,:,roi_idx_sm),[],2));
     
% sm_roi_scale = squeeze(sqrt(sum(fluor_roi_deconv.^2,2))./ ...
%     sqrt(sum(fluor_roi_deconv(:,:,roi_idx_sm).^2,2)));

% sm_roi_scale = squeeze(std(fluor_roi_deconv,[],2)./ ...
%     std(fluor_roi_deconv(:,:,roi_idx_sm),[],2));

% (trial-by-trial scaling: must be efficient way to do this...)
sm_roi_scale = nan(size(fluor_roi_deconv,1),1,n_rois);
for curr_roi = 1:n_rois
   for curr_trial = 1:size(fluor_roi_deconv)
      sm_roi_scale(curr_trial,curr_roi) =  ...
          fluor_roi_deconv(curr_trial,:,roi_idx_sm)'\ ...
          fluor_roi_deconv(curr_trial,:,curr_roi)';
   end
end

fluor_roi_deconv_smsub = fluor_roi_deconv - ...
    fluor_roi_deconv(:,:,roi_idx_sm).*sm_roi_scale;

figure;
h = tiledlayout(2,length(plot_rois),'TileSpacing','compact','padding','compact');
for curr_stage = 1:max(trial_learned_stage)
    for curr_roi = plot_rois
        
        use_trials = find(trial_learned_stage == curr_stage);
        [~,sort_idx] = sort(move_t(use_trials));

        % (scale SM ROI and subtract)
        curr_data_sort = fluor_roi_deconv_smsub(use_trials(sort_idx),:,curr_roi);
        curr_data_sort_smooth = convn(curr_data_sort, ...
            ones(n_trial_smooth,1)./n_trial_smooth,'same');
        
        nexttile;
        imagesc(t,[],curr_data_sort_smooth);hold on;
        colormap(brewermap([],'*RdGy'));
        
        c = prctile(reshape(fluor_roi_deconv(:,:,curr_roi),[],1),90).*[-1,1];
        caxis(c);
        xline(0,'color','r');
        plot(move_t(use_trials(sort_idx)),1:length(use_trials),'color',[0.6,0,0.6]);
        xlabel('Time from stim (s)');
        ylabel('Trial (rxn-time sorted)');
        title(sprintf('%s, stage %d',wf_roi(curr_roi).area,curr_stage));
    end
end
title(h,'Activity - scaled SM');


        
        
        
        
%% >> Muscimol: task

trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_task_teto_muscimol';

AP_load_trials_operant;

% Get animal and day index for each trial
trial_animal = cell2mat(arrayfun(@(x) ...
    x*ones(size(vertcat(wheel_all{x}{:}),1),1), ...
    [1:length(wheel_all)]','uni',false));

% Choose split for data
trials_animal = arrayfun(@(x) size(vertcat(wheel_all{x}{:}),1),1:size(wheel_all));
trials_recording = cellfun(@(x) size(x,1),vertcat(wheel_all{:}));

% Get task>cortex parameters
n_regressors = length(task_regressor_labels);
task_regressor_t_shifts = cellfun(@(x) x/sample_rate,task_regressor_sample_shifts,'uni',false);

% Average kernels by muscimol area
muscimol_area_cat = horzcat(muscimol_area{:})';
muscimol_area_allcat = arrayfun(@(x) ...
    repmat(muscimol_area_cat(x),trials_recording(x),1), ...
    1:length(trials_recording),'uni',false)';
muscimol_area_allcat = vertcat(muscimol_area_allcat{:});

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

%% average task activity

use_rxn = [0.1,0.25];

act_avg_animal_v = nan(n_vs,length(t),length(unique_muscimol_area),length(animals));
act_avg_animal_roi = nan(n_rois,length(t),length(unique_muscimol_area),length(animals));

for curr_animal = 1:length(animals)
    for curr_muscimol = 1:length(unique_muscimol_area)
        
        curr_trials = trial_outcome_allcat == 1 & ...
            curr_animal == trial_animal & ...
            move_t >= use_rxn(1) & move_t <= use_rxn(2) & ...
            strcmp(muscimol_area_allcat,unique_muscimol_area{curr_muscimol});
        
        %         % (raw)
        %         act_avg_animal_v(:,:,curr_muscimol,curr_animal) = ...
        %             permute(nanmean(fluor_allcat_deconv(curr_trials,:,:),1),[3,2,1]);
        %         act_avg_animal_roi(:,:,curr_muscimol,curr_animal) = ...
        %             permute(nanmean(fluor_roi_deconv(curr_trials,:,:),1),[3,2,1]);
        %
        % (reduced)
        use_reduce = 2;
        act_avg_animal_v(:,:,curr_muscimol,curr_animal) = ...
            permute(nanmean(fluor_allcat_deconv(curr_trials,:,:) - ...
            fluor_taskpred_reduced_allcat(curr_trials,:,:,use_reduce),1),[3,2,1]);
        act_avg_animal_roi(:,:,curr_muscimol,curr_animal) = ...
            permute(nanmean(fluor_roi_deconv(curr_trials,:,:) - ...
            fluor_roi_taskpred_reduced(curr_trials,:,:,use_reduce),1),[3,2,1]);
        
    end
end

act_avg_px = AP_svdFrameReconstruct(U_master(:,:,1:n_vs),nanmean(act_avg_animal_v,4));

AP_image_scroll(act_avg_px,t);
caxis([-max(abs(caxis)),max(abs(caxis))]);
colormap(brewermap([],'PrGn'));
axis image;
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);


plot_roi = 7;
plot_muscimol = [3,4];
curr_roi_act = squeeze(act_avg_animal_roi(plot_roi,:,plot_muscimol,:));
figure;
h = AP_errorfill(t,nanmean(curr_roi_act,3),AP_sem(curr_roi_act,3));
xlabel('Time from stim');
ylabel(wf_roi(plot_roi).area);
legend(h,unique_muscimol_area(plot_muscimol));



%% >> Muscimol: passive

% Load data
trial_data_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\analysis\operant_learning\data';
data_fn = 'trial_activity_passive_teto_muscimol';
% data_fn = 'trial_activity_passive_cstr_muscimol';

AP_load_trials_operant;

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

AP_image_scroll(fluor_muscimol_avg_px,t)
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

%% >> Ephys passive

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
    
end




%% >> Ephys task

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
tiledlayout('flow');
figure;
for curr_regressor = 1:length(mua_area_taskpred_k_avg)
    for curr_subregressor = 1:size(mua_area_taskpred_k_avg{curr_regressor},1)
        nexttile;
        plot(task_regressor_sample_shifts{curr_regressor}/sample_rate, ...
            squeeze(mua_area_taskpred_k_avg{curr_regressor}(curr_subregressor,:,plot_areas)));
        title(sprintf('Regressor %d, sub %d',curr_regressor,curr_subregressor));
        xline(0,'color','k','linestyle','--');
    end
end









