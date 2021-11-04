%% Plot example performance

animal = 'AP100';
% day = '2021-05-02';
day = '2021-05-07';

experiment = 1;
verbose = true;
load_parts.imaging = false;
AP_load_experiment;

t = Timeline.rawDAQTimestamps;

% (minimum ITI: new trial + trial quiescence)
min_iti_t = signals_events.newTrialTimes + ...
    signals_events.trialQuiescenceValues;

plot_t = [0,100];
plot_t_idx = t > plot_t(1) & t < plot_t(2);
plot_stim_idx = find(stimOn_times > plot_t(1) & stimOn_times < plot_t(2))';
plot_min_iti_t_idx = find(min_iti_t > plot_t(1) & min_iti_t < plot_t(2));
plot_reward_idx = find(reward_t_timeline > plot_t(1) & reward_t_timeline < plot_t(2));

figure; hold on
plot(t(plot_t_idx),wheel_velocity(plot_t_idx),'k');
for i = plot_stim_idx
   line(repmat(stimOn_times(i),2,1),ylim,'color','r'); 
end
for i = plot_min_iti_t_idx'
   line(repmat(min_iti_t(i),2,1),ylim,'linestyle','--','color','g'); 
end
for i = plot_reward_idx'
   line(repmat(reward_t_timeline(i),2,1),ylim,'color','b'); 
end


%% Make movie after averaging activity

% use_movie_t = t > -0.2 & t <= 0.1;
% use_movie_t = t >= 0.1 & t <= 0.3;
% use_movie_t = t >= 0.3 & t <= 0.7;
use_movie_t = t > -0.2 & t <= 0.7;

movie_rate = sample_rate/5;
color_map = brewermap([],'PrGn');
% color_axis = [-0.008,0.008];
% color_axis = [-0.002,0.002];
color_axis = [-0.005,0.005];
% figure_position = [105 497 500 483];
figure_position = [25 600 1300 300];
save_path = 'C:\Users\Andrew\OneDrive for Business\Documents\CarandiniHarrisLab\presentations\211011_data_club';
savefile = [save_path filesep 'avg_fluor_teto_rxn_reduced'];
t_annotation = cellfun(@(x) sprintf('Time from stimulus: %0.2f sec',x),num2cell(t),'uni',false);
AP_movie2avi(curr_px(:,:,use_movie_t,:), ...
    movie_rate,color_map,color_axis,figure_position,savefile,t_annotation(use_movie_t));





















