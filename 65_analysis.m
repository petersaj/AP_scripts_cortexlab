
%% Run kilosort on raw data
% Commented so I don't accidentally do it again

% animal = '65';
% days = {'20151028','20151029','20151030','20151031','20151101','20151102','20151103'};
% input_board = 'whisper';
% 
% for curr_day = 1:length(days)
%     data_filename = ['\\zserver.cortexlab.net\Data\multichanspikes\' ...
%         animal filesep days{curr_day} filesep days{curr_day} '_1.dat'];
%     
%     AP_run_kilosort(data_filename,input_board);
%     
%     disp(['Finished ' num2str(curr_day) '/' num2str(length(days))]);
% end


%% Load

animal = '65';
days = {'20151028','20151029','20151030','20151031','20151101','20151102','20151103'};

curr_day = 6;
[spikes,xpr] = AP_load_ephys(animal,days{curr_day});




%% PSTH GUI

trial_condition_rl = (xpr.bhv.condition(1,:) > xpr.bhv.condition(2,:)) - ...
    (xpr.bhv.condition(2,:) > xpr.bhv.condition(1,:));
spike_times_sec = double(spikes.spike_times);
raster_window = [-2,2];
psthViewer(spike_times_sec,spikes.spike_clusters, ...
    xpr.bhv.stim_onset_ephys,raster_window,trial_condition_rl);















