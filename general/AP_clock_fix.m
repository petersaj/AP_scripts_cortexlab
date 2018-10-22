function t_new = AP_clock_fix(t,source_sync,target_sync)
% t_new = AP_clock_fix(t,source_sync,target_sync)
% 
% NOTE: this is literally just interp1
%
% Converts time from one clock to another
% t - times
% source_sync - synchronization corresponding to t
% target_sync - synchronization corresponding to desired t

% The number of syncs must be the same
if length(source_sync) ~= length(target_sync)
    error('Different numbers of sync pulses - cannot align')
end

t_new = interp1(source_sync,target_sync,t);






