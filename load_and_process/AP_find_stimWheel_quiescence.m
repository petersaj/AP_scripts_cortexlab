function quiescence_reset_t = AP_find_stimWheel_quiescence
% Reconstruct sig.quiescenceWatch times for AP_stimWheel protocols
% (it isn't saved in the block because it's an internal process)

% Pull in expected variables from workspace
try
    n_trials = evalin('base','n_trials');
    signals_events = evalin('base','signals_events');
    wheel_move_stim_idx = evalin('base','wheel_move_stim_idx');
    wheel_starts = evalin('base','wheel_starts');
    t = evalin('base','Timeline.rawDAQTimestamps''');
    block2timeline = evalin('base','block2timeline');
    timeline2block = evalin('base','timeline2block');
    block = evalin('base','block');
    wheel_position = evalin('base','wheel_position');
    wheel_velocity = evalin('base','wheel_velocity');
    stimOn_times = evalin('base','stimOn_times');
catch me
    missing_var = textscan(me.message,'Undefined function or variable '' %s');
    error(['Variable missing from base workspace: ' cell2mat(missing_var{:})]);
end

% (to plot for debugging)
plot_trial = false;
if plot_trial
    plot_trial_fig = figure;
end

% (skip first trial: no ITI, also quiescenceWatch doesn't seem to be active?)
quiescence_reset_t_split = cell(n_trials,1);
for curr_trial = 2:n_trials    
    
    % Pull out current trial times (last response to first post-stim move)
    curr_trialpull_start = signals_events.responseTimes(curr_trial-1);
    curr_trialpull_end = signals_events.responseTimes(curr_trial);
    
    curr_trial_t_idx = t >= curr_trialpull_start & t <= curr_trialpull_end;
    curr_trial_t = t(curr_trial_t_idx);
    
    % Get trial params
    % (iti is defined in last trial)
    curr_iti = signals_events.trialITIValues(curr_trial-1);
    curr_trialstart = signals_events.newTrialTimes(curr_trial);
    curr_quiescence = signals_events.trialQuiescenceValues(curr_trial);
    
    % Re-create quiescence watch: resets at cumulative abs < 1mm
    
    thresh_mm = 1; % (hardcoded in protocol)
    
    %%% Quiescence watch in block (what signals really uses)
    t_wheel_block = interp1(block2timeline,timeline2block,block.inputs.wheelMMTimes,'linear','extrap');
    curr_trial_t_block_idx = t_wheel_block >= curr_trialpull_start & t_wheel_block <= curr_trialpull_end;
    curr_wheel_mm_t = t_wheel_block(curr_trial_t_block_idx);
    curr_wheel_mm = block.inputs.wheelMMValues(curr_trial_t_block_idx);
    curr_quiescence_reset_block = false(size(curr_wheel_mm));
    % (quiescence watch starts on new trial)
    i = find(curr_wheel_mm_t >= curr_trialstart,1);
    while i < length(curr_wheel_mm)
        next_thresh_cross = min([length(curr_wheel_mm),(i-1) + ...
            find(cumsum(abs(diff(curr_wheel_mm(i:end)))) > thresh_mm,1,'first')]);
       
        curr_quiescence_reset_block(next_thresh_cross) = true;
        
        i = next_thresh_cross + 1;
    end
    quiescence_reset_t_block = [curr_trialstart;curr_wheel_mm_t(curr_quiescence_reset_block)'];
    
    t_from_quiescence_reset_full_block = t - ...
        interp1(quiescence_reset_t_block,quiescence_reset_t_block,t,'previous','extrap');
    t_from_quiescence_reset_trial_block = t_from_quiescence_reset_full_block(curr_trial_t_idx);
    curr_reconstructed_stimOn_block = ...
        curr_trial_t(find(t_from_quiescence_reset_trial_block > curr_quiescence,1));
    
    %%% Quiescence watch in timeline (should closely match block)
    mm_per_wheelclick = min(abs(diff(block.inputs.wheelMMValues)));
    thresh_wheelclick = thresh_mm/mm_per_wheelclick;
    curr_wheel_timeline = wheel_position(curr_trial_t_idx);
    curr_wheel_timeline_diff = abs([0;diff(curr_wheel_timeline)]);
    curr_quiescence_reset_timeline = false(size(curr_wheel_timeline));
    % (quiescence watch starts on new trial)
    i = find(curr_trial_t >= curr_trialstart,1);
    while i < length(curr_wheel_timeline)       
        next_thresh_cross = min([length(curr_wheel_timeline),(i-1) + ...
            find(cumsum(abs(diff(curr_wheel_timeline(i:end)))) > thresh_wheelclick,1,'first')]);       
      
        curr_quiescence_reset_timeline(next_thresh_cross) = true;      
       
        i = next_thresh_cross;        
    end
    quiescence_reset_t_timeline = ...
        [curr_trialstart,curr_trial_t(curr_quiescence_reset_timeline)'];
    
    t_from_quiescence_reset_full_timeline = t - ...
        interp1(quiescence_reset_t_timeline,quiescence_reset_t_timeline,t,'previous','extrap');
    t_from_quiescence_reset_trial_timeline = t_from_quiescence_reset_full_timeline(curr_trial_t_idx);
    curr_reconstructed_stimOn_timeline = ...
        curr_trial_t(find(t_from_quiescence_reset_trial_timeline > curr_quiescence,1));
    
    if plot_trial
        figure(plot_trial_fig);
        clf;hold on;
        t_plot_scale = 0.1;
        plot(curr_trial_t,wheel_velocity(curr_trial_t_idx),'k')
        plot(curr_trial_t,[0;diff(wheel_position(curr_trial_t_idx))]*0.1,'r')
        plot(curr_wheel_mm_t,0,'.b');
        line(repmat(curr_trial_t(1)+curr_iti,2,1),ylim);
        line(xlim,repmat(curr_quiescence,2,1)*t_plot_scale,'color','m');
        plot(curr_trial_t,t_from_quiescence_reset_trial_block*t_plot_scale,'b');
        plot(curr_trial_t,t_from_quiescence_reset_trial_timeline*t_plot_scale,'r');
        
        line(repmat(stimOn_times(curr_trial),1,2),ylim,'color','k','linestyle','--');
        line(repmat(curr_reconstructed_stimOn_block,1,2),ylim,'color','b','linestyle','--');
        if ~isempty(curr_reconstructed_stimOn_timeline)
            line(repmat(curr_reconstructed_stimOn_timeline,1,2),ylim,'color','r','linestyle','--');
        end
        
        title(sprintf('Trial %d/%d',curr_trial,n_trials));
        
        drawnow;
    end
    
    % Sanity check: signals stimOn should match reconstructed timings
    % (extremely closely - within 5ms)
    reconstructed_stimOn_leeway = 0.005;   
    curr_reconstructed_stimOn_error = ...
        abs(curr_reconstructed_stimOn_block - signals_events.stimOnTimes(curr_trial));
    if curr_reconstructed_stimOn_error > reconstructed_stimOn_leeway
        warning('Timing-reconstructed stimOn time incorrect: debugging');
        keyboard
    end
    
    % Store quiescence reset times
    quiescence_reset_t_split{curr_trial} = quiescence_reset_t_block;
        
end

% Concatenate all reset times
quiescence_reset_t = cell2mat(quiescence_reset_t_split);






