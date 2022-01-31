function quiescence_reset_t = AP_find_stimWheel_quiescence
% Reconstruct sig.quiescenceWatch times for AP_stimWheel protocols
% (it isn't saved in the block because it's an internal process)

% Pull in expected variables from workspace
try
    % (from loading in experiment)
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
    
    % (read from expDef)
    quiescThreshold = evalin('base','quiescThreshold');
catch me
    missing_var = textscan(me.message,'Undefined function or variable '' %s');
    error(['Variable missing from base workspace: ' cell2mat(missing_var{:})]);
end

% (to plot for debugging)
plot_trial = false;
if plot_trial
    plot_trial_fig = figure;
end

% Loop through trials, get quiescence resets
% (skip first trial: quiescence watch doesn't work?)
quiescence_reset_t_split = cell(n_trials,1);
stimOn_estimation_error = nan(n_trials,1);
backwards_quiescence_estimation = false(n_trials,1);
for curr_trial = 2:n_trials
    
    % Pull out time quiescence watch time (new trial to stim)
    curr_qwatch_start = signals_events.newTrialTimes(curr_trial);
    curr_qwatch_end = signals_events.stimOnTimes(curr_trial) + 0.2; %(leeway)
    
    curr_trial_t_idx = t >= curr_qwatch_start & t <= curr_qwatch_end;
    curr_trial_t = t(curr_trial_t_idx);
    
    % Get trial params
    curr_trialstart = signals_events.newTrialTimes(curr_trial);
    curr_quiescence = signals_events.trialQuiescenceValues(curr_trial);
    
    % Re-create quiescence watch: resets at cumulative abs < quiescThreshold mm
        
    %%% Quiescence watch in block (what signals really uses)
    t_wheel_block = interp1(block2timeline,timeline2block,block.inputs.wheelMMTimes,'linear','extrap');
    curr_trial_t_block_idx = t_wheel_block >= curr_qwatch_start & t_wheel_block <= curr_qwatch_end;
    curr_wheel_mm_t = t_wheel_block(curr_trial_t_block_idx);
    curr_wheel_mm = block.inputs.wheelMMValues(curr_trial_t_block_idx);
    
    % (quiescence watch starts on new trial)
    q_start_idx = find(curr_wheel_mm_t >= curr_trialstart,1);
    
    curr_quiescence_reset_block = false(size(curr_wheel_mm));
    while q_start_idx < length(curr_wheel_mm)
        next_thresh_cross = (q_start_idx-1) + ...
            find(cumsum(abs(diff(curr_wheel_mm(q_start_idx:end)))) > ...
            quiescThreshold,1,'first');
        
        curr_quiescence_reset_block(next_thresh_cross) = true;
        
        q_start_idx = next_thresh_cross + 1;
    end
    quiescence_reset_t_block = [curr_trialstart;curr_wheel_mm_t(curr_quiescence_reset_block)'];
    
    t_from_quiescence_reset_full_block = t - ...
        interp1([quiescence_reset_t_block;curr_qwatch_end],...
        [quiescence_reset_t_block;curr_qwatch_end],t,'previous','extrap');
    t_from_quiescence_reset_trial_block = t_from_quiescence_reset_full_block(curr_trial_t_idx);
    
    curr_reconstructed_stimOn_block = ...
        curr_trial_t(find(t_from_quiescence_reset_trial_block > curr_quiescence,1));
    curr_reconstructed_stimOn_error = ...
        min([Inf,abs(curr_reconstructed_stimOn_block - signals_events.stimOnTimes(curr_trial))]);
    
    % NOTE: it looks like sometimes there was some quiescence watch
    % weirdness if a wheel click happened too close to the trial time.
    % So - if doing it the "right" way by going forward doesn't work,
    % attempt to reconstruct by going backwards.
    % This is risky: it will always find a right answer even if something
    % else (like the quiescence threshold) is wrong.
    
    % Sanity check: signals stimOn should match reconstructed timings
    % (extremely closely)
    reconstructed_stimOn_max_error = 0.01;
    
    if curr_reconstructed_stimOn_error > reconstructed_stimOn_max_error
        % Estimate quiescence watch backwards from stim
        last_quiescence_reset = signals_events.stimOnTimes(curr_trial) - curr_quiescence;
        % (get closest wheel click to last quiescence reset)
        [~,last_quiescence_reset_wheel_idx] = min(abs(curr_wheel_mm_t - last_quiescence_reset));
        
        q_end_idx = last_quiescence_reset_wheel_idx;
        
        curr_quiescence_reset_block = false(size(curr_wheel_mm));
        while q_end_idx > find(curr_wheel_mm_t >= curr_trialstart,1)
            curr_quiescence_reset_block(q_end_idx) = true;
            
            prior_thresh_cross = max([1,(1+q_end_idx) - ...
                find(cumsum(fliplr(abs(diff(curr_wheel_mm(1:q_end_idx))))) > quiescThreshold,1,'first')]);
            
            q_end_idx = prior_thresh_cross;
        end
        quiescence_reset_t_block = [curr_trialstart;curr_wheel_mm_t(curr_quiescence_reset_block)'];
        
        t_from_quiescence_reset_full_block = t - ...
            interp1([quiescence_reset_t_block;curr_qwatch_end],...
            [quiescence_reset_t_block;curr_qwatch_end],t,'previous','extrap');
        t_from_quiescence_reset_trial_block = t_from_quiescence_reset_full_block(curr_trial_t_idx);
        
        curr_reconstructed_stimOn_block = ...
            curr_trial_t(find(t_from_quiescence_reset_trial_block > curr_quiescence,1));
        curr_reconstructed_stimOn_error = ...
            min([Inf,abs(curr_reconstructed_stimOn_block - signals_events.stimOnTimes(curr_trial))]);
        
        backwards_quiescence_estimation(curr_trial) = true;
    end
    
    % If the estimated stim time still is over the error thresh, debug
    if curr_reconstructed_stimOn_error > reconstructed_stimOn_max_error
        warning('Quiescence watch estimation wrong, debugging')
        plot_trial = true;
    end
    
    %     %%% Quiescence watch in timeline (should closely match block)
    %     % (isn't actually used, but I did this for comparison)
    %     mm_per_wheelclick = min(abs(diff(block.inputs.wheelMMValues)));
    %     quiescThreshold_wheelclick = quiescThreshold/mm_per_wheelclick;
    %     curr_wheel_timeline = wheel_position(curr_trial_t_idx);
    %     curr_wheel_timeline_diff = abs([0;diff(curr_wheel_timeline)]);
    %     curr_quiescence_reset_timeline = false(size(curr_wheel_timeline));
    %     % (quiescence watch starts on new trial)
    %     q_start_idx = find(curr_trial_t >= curr_trialstart,1);
    %     while q_start_idx < length(curr_wheel_timeline)
    %         next_thresh_cross = (q_start_idx-1) + ...
    %             find(cumsum(abs(diff(curr_wheel_timeline(q_start_idx:end)))) > quiescThreshold_wheelclick,1,'first');
    %
    %         curr_quiescence_reset_timeline(next_thresh_cross) = true;
    %
    %         q_start_idx = next_thresh_cross;
    %     end
    %     quiescence_reset_t_timeline = ...
    %         [curr_trialstart;curr_trial_t(curr_quiescence_reset_timeline)];
    %
    %     t_from_quiescence_reset_full_timeline = t - ...
    %         interp1([quiescence_reset_t_timeline;curr_qwatch_end],...
    %         [quiescence_reset_t_timeline;curr_qwatch_end],t,'previous','extrap');
    %     t_from_quiescence_reset_trial_timeline = t_from_quiescence_reset_full_timeline(curr_trial_t_idx);
    %     curr_reconstructed_stimOn_timeline = ...
    %         curr_trial_t(find(t_from_quiescence_reset_trial_timeline > curr_quiescence,1));
    
    % Plot trial (and debug)
    if plot_trial
        if ~exist('plot_trial_fig','var') || ~ishandle(plot_trial_fig)
            plot_trial_fig = figure;
        end
        figure(plot_trial_fig);
        clf;hold on;
        t_plot_scale = 0.1;
        plot(curr_trial_t,wheel_velocity(curr_trial_t_idx),'k')
        plot(curr_trial_t,[0;diff(wheel_position(curr_trial_t_idx))]*0.1,'r')
        plot(curr_wheel_mm_t,0,'.b');
        line(repmat(curr_trialstart,2,1),ylim,'color','b');
        line(xlim,repmat(curr_quiescence,2,1)*t_plot_scale,'color','m');
        plot(curr_trial_t,t_from_quiescence_reset_trial_block*t_plot_scale,'b');
        
        line(repmat(signals_events.stimOnTimes(curr_trial),1,2),ylim,'color','k','linestyle','--');
        line(repmat(curr_reconstructed_stimOn_block,1,2),ylim,'color','b','linestyle','--');
        
        title(sprintf('Trial %d/%d',curr_trial,n_trials));
        
        drawnow;
        keyboard % (only plot on debug at the moment)
    end
    
    % Store quiescence reset times
    quiescence_reset_t_split{curr_trial} = quiescence_reset_t_block;
    stimOn_estimation_error(curr_trial) = curr_reconstructed_stimOn_error;
    
    %     AP_print_progress_fraction(curr_trial,n_trials);
    
end

% If more than x% backwards-estimated quiescence, error out
backwards_quiescence_estimation_frac = nanmean(backwards_quiescence_estimation);
if backwards_quiescence_estimation_frac > 0.5
   error('%d%% backwards-estimated quiescence trials', ...
       round(backwards_quiescence_estimation_frac*100));
end

% Add first quiescence reset as first trial time
quiescence_reset_t_split{1} = signals_events.newTrialTimes(1);

% Concatenate all reset times
quiescence_reset_t = cell2mat(quiescence_reset_t_split);






