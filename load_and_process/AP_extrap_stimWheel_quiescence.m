function quiescence_reset_t = AP_extrap_stimWheel_quiescence
% Extrapolate sig.quiescenceWatch times for AP_stimWheel protocols from the
% actual stimOn time (to find would-be quiescence periods given different
% trial parameters)

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
for curr_trial = 2:n_trials
    
    % Pull out time quiescence watch time (new trial to stim)
    curr_qwatch_start = signals_events.responseTimes(curr_trial-1);
    curr_qwatch_end = signals_events.responseTimes(curr_trial);
    
    curr_trial_t_idx = t >= curr_qwatch_start & t <= curr_qwatch_end;
    curr_trial_t = t(curr_trial_t_idx);
    
    % Get trial params
    curr_trialstart = signals_events.newTrialTimes(curr_trial);
    curr_quiescence = signals_events.trialQuiescenceValues(curr_trial);
    curr_iti = signals_events.trialITIValues(curr_trial-1);
    
    % Re-create quiescence watch: resets at cumulative abs <
    % quiescThreshold mm (all using block data)
    
    t_wheel_block = interp1(block2timeline,timeline2block,block.inputs.wheelMMTimes,'linear','extrap');
    curr_trial_t_block_idx = t_wheel_block >= curr_qwatch_start & t_wheel_block <= curr_qwatch_end;
    curr_wheel_mm_t = t_wheel_block(curr_trial_t_block_idx);
    curr_wheel_mm = block.inputs.wheelMMValues(curr_trial_t_block_idx);
    
    curr_quiescence_reset_block = false(size(curr_wheel_mm));
    
    %%% Extrapolate quiescence backwards from stim
    last_quiescence_reset = signals_events.stimOnTimes(curr_trial) - curr_quiescence;
    % (2 options: the ITI was the last reset, or there was post-ITI
    % movement that was the last reset)
    post_iti_q_reset_flag = min(abs(last_quiescence_reset - curr_trialstart)) >= ...
        min(abs(last_quiescence_reset-curr_wheel_mm_t));    
    
    if post_iti_q_reset_flag
        % (if there was a wheel click closest, use that)
        [~,last_quiescence_reset_wheel_idx] = min(abs(curr_wheel_mm_t - last_quiescence_reset));
    else
        % (if the trial start was closest, use the wheel click prior to that)
        last_quiescence_reset_wheel_idx = find(curr_wheel_mm_t < ...
            curr_trialstart,1,'last');
    end
    
    
    q_end_idx = last_quiescence_reset_wheel_idx;   
    while q_end_idx > 1
        curr_quiescence_reset_block(q_end_idx) = true;
        
        prior_thresh_cross = max([1,(1+q_end_idx) - ...
            find(cumsum(fliplr(abs(diff(curr_wheel_mm(1:q_end_idx))))) > quiescThreshold,1,'first')]);
        
        if q_end_idx == prior_thresh_cross
            break
        end
        
        q_end_idx = prior_thresh_cross;
    end

    %%% Extrapolate forwards from last quiecence reset    
    q_start_idx = max([1,find(curr_quiescence_reset_block,1,'last')+1]);
    while q_start_idx < length(curr_wheel_mm)
        next_thresh_cross = (q_start_idx-1) + ...
            find(cumsum(abs(diff(curr_wheel_mm(q_start_idx:end)))) > ...
            quiescThreshold,1,'first');
        
        curr_quiescence_reset_block(next_thresh_cross) = true;
        
        q_start_idx = next_thresh_cross + 1;
    end
    
    quiescence_reset_t_block = [curr_trial_t(1);curr_wheel_mm_t(curr_quiescence_reset_block)'];   
        
    t_from_quiescence_reset_full_block = t - ...
        interp1(quiescence_reset_t_block,quiescence_reset_t_block,t,'previous','extrap');
    t_from_quiescence_reset_trial_block = t_from_quiescence_reset_full_block(curr_trial_t_idx);
    
    % (get reconstructed stimOn time given actual trialstart/ITI
    quiescence_reset_t_block_curriti = sort([quiescence_reset_t_block;curr_trialstart]);
    
    t_from_quiescence_reset_full_block_curriti = t - ...
        interp1(quiescence_reset_t_block_curriti,quiescence_reset_t_block_curriti,t,'previous','extrap');
    t_from_quiescence_reset_trial_block_curriti = ...
        t_from_quiescence_reset_full_block_curriti(curr_trial_t_idx);
    
    curr_reconstructed_stimOn_block = ...
        curr_trial_t(find( ...
        t_from_quiescence_reset_trial_block_curriti > curr_quiescence & ...
        curr_trial_t > curr_trialstart,1));
    curr_reconstructed_stimOn_error = ...
        min([Inf,abs(curr_reconstructed_stimOn_block - signals_events.stimOnTimes(curr_trial))]);
    
    % Don't check for stimOn time error at the moment - return a valid set
    % of quiescence resets, but sometimes a skipped wheel click means it
    % can't replicate the real stimOn time
%     % If the estimated and matched stim times don't match, error out
%     reconstructed_stimOn_max_error = 0.05;
%     if curr_reconstructed_stimOn_error > reconstructed_stimOn_max_error
%         warning('Quiescence watch estimation wrong, debugging')
%         plot_trial = true;
%     end
    
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
        plot(curr_trial_t,t_from_quiescence_reset_trial_block_curriti*t_plot_scale,'m');
        plot(quiescence_reset_t_block,0,'og');
        
        line(repmat(signals_events.stimOnTimes(curr_trial),1,2),ylim,'color','k');
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

% Add first quiescence reset as first trial time
quiescence_reset_t_split{1} = signals_events.newTrialTimes(1);

% Concatenate all reset times
quiescence_reset_t = cell2mat(quiescence_reset_t_split);






