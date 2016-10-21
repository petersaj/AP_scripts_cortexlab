function t_new = AP_clock_fix(t,source_sync,target_sync)
% t_new = AP_clock_fix(t,source_sync,target_sync)
% 
% Converts time from one clock to another
% t - times
% source_sync - synchronization corresponding to t
% target_sync - synchronization corresponding to desired t
% NOTE: source/target syncs can be n timepoints, only uses first and last

if length(source_sync) ~= length(target_sync)
    error('Different numbers of sync pulses')
end

t_offset = source_sync(1) - target_sync(1);
t_drift = 1 - (target_sync(end) - target_sync(1))/(source_sync(end) - source_sync(1));

t_new = t - t_offset - (t - t(1)).*t_drift;





