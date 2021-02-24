classdef choiceWorldJFExpPanel < eui.SignalsExpPanel
  %eui.SignalsExpPanel Basic UI control for monitoring an experiment
  %   TODO
  %
  % Part of Rigbox
  
  properties
    % Trials per minute window size
    WindowSize = 20
    RewardUnits = sprintf('%cl', char(956))
  end
  
  properties (Access = protected)
    PsychometricAxes % Handle to axes of psychometric plot
    ExperimentAxes % Handle to axes of wheel trace and threhold line plot
    ExperimentHands % handles to plot objects in the experiment axes
    ScreenAxes
    ScreenHands
    VelAxes
    VelHands
    InputSensorPosTime % Vector of timesstamps in seconds for plotting the wheel trace
    InputSensorPos % Vector of azimuth values for plotting the wheel trace
    InputSensorPosCount = 0 % Running total of azimuth samples recieved for axes plot
    ExtendThresholdLines = false % Flag for plotting dotted threshold lines during cue interactive delay.  Currently unused.
    lastEvtTime = now;
  end
  
  properties (Access = protected, SetObservable)
    contrastLeft = []
    contrastRight = []
  end
  
  methods
    function obj = choiceWorldJFExpPanel(parent, ref, params, logEntry)
      obj = obj@eui.SignalsExpPanel(parent, ref, params, logEntry);
      % Initialize InputSensor properties for speed
      obj.InputSensorPos = nan(1000*30, 1);
      obj.InputSensorPosTime = nan(1000*30, 1);
      obj.InputSensorPosCount = 0;
      obj.Block.numCompletedTrials = -1; % Store sum of events.newTrial
      obj.Block.totalReward = 0; % Store sum of outputs.reward event here
      obj.Block.trial = struct('contrastLeft', [], 'contrastRight', [], ...
        'response', [], 'repeatNum', [], 'feedback', [],...
        'wheelGain', [], 'newTrialTimes', [], ...
        'baselineRT', [], 'windowedRT', [], 'tMin', []);
      % Add our reward info field.  Unlike TrialCountLabel we keep this in
      % the LabelsMap so that it updates green when changing
      obj.LabelsMap('outputs.reward') = obj.addInfoField('Reward', ['0 ' obj.RewardUnits]);
      set(obj.LabelsMap('outputs.reward'), 'UserData', clock);
      % Set which events we wish to display as info fields
      obj.UpdatesFilter = [
        "events.repeatNum", ...
        "events.disengaged", ...
        "events.pctDecrease", ...
        "events.proportionLeft", ...
        "inputs.lick", ...
        "events.prestimQuiescentPeriod"];
      obj.Exclude = false; % Include only those in the UpdatesFilter list
      % Prettify the InfoLabels
      obj.FormatLabels = true;
      obj.Listeners = [obj.Listeners,...
        addlistener(obj,'contrastLeft','PostSet',@obj.setContrast), ...
        addlistener(obj,'contrastRight','PostSet',@obj.setContrast)];
    end
    
  end
  
  methods (Access = protected)
    function setContrast(obj, ~, ~)
      cL = obj.contrastLeft;
      cR = obj.contrastRight;
      obj.Parameters.Struct.stimulusContrast = [cL cR];
    end
    
    function processUpdates(obj)
      updates = obj.SignalUpdates(1:obj.NumSignalUpdates);
      obj.NumSignalUpdates = 0;
      
      if ~isempty(updates)
        %fprintf('processing %i signal updates\n', length(updates));
        
        % pull out wheel updates
        allNames = {updates.name};
        wheelUpdates = strcmp(allNames, 'inputs.wheelMM');
        
        if sum(wheelUpdates)>0
          x = -[updates(wheelUpdates).value];
          t = (24*3600*cellfun(@(x)datenum(x), {updates(wheelUpdates).timestamp}))-(24*3600*obj.StartedDateTime);
          
          nx = numel(x);
          obj.InputSensorPosCount = obj.InputSensorPosCount+nx;
          
          if obj.InputSensorPosCount>numel(obj.InputSensorPos)
            % full - drop the first half of the array and shift the
            % last half back
            halfidx = floor(numel(obj.InputSensorPos)/2);
            obj.InputSensorPos(1:halfidx) = obj.InputSensorPos(halfidx:2*halfidx-1);
            obj.InputSensorPos(halfidx+1:end) = NaN;
            obj.InputSensorPosTime(1:halfidx) = obj.InputSensorPosTime(halfidx:2*halfidx-1);
            obj.InputSensorPosTime(halfidx+1:end) = NaN;
            obj.InputSensorPosCount = obj.InputSensorPosCount-halfidx;
          end
          obj.InputSensorPos(obj.InputSensorPosCount-nx+1:obj.InputSensorPosCount) = x;
          obj.InputSensorPosTime(obj.InputSensorPosCount-nx+1:obj.InputSensorPosCount) = t;
          
        end
        
        % now plot the wheel
        plotwindow = [-5 1];
        lastidx = obj.InputSensorPosCount;
        
        if lastidx > 1
          
          firstidx = find(obj.InputSensorPosTime>obj.InputSensorPosTime(lastidx)+plotwindow(1),1);
          
          xx = obj.InputSensorPos(firstidx:lastidx);
          tt = obj.InputSensorPosTime(firstidx:lastidx);
          
          set(obj.ExperimentHands.wheelH, 'XData', xx, 'YData', tt);
          
          set(obj.ExperimentAxes.Handle, 'YLim', plotwindow + tt(end));
          % update the velocity tracker too
          [tt, idx] = unique(tt);
          if numel(tt) > 1
            recentX = interp1(tt, xx(idx), tt(end)+(-0.3:0.05:0));
            vel = mean(diff(recentX));
          else
            vel = 0;
          end
          set(obj.VelHands.Vel, 'XData', vel*[1 1]);
          obj.VelHands.MaxVel = max(abs([obj.VelHands.MaxVel vel]));
          set(obj.VelAxes, 'XLim', obj.VelHands.MaxVel*[-1 1]);
          
        end
        
        % now deal with other updates
        updates = updates(~wheelUpdates);
        allNames = allNames(~wheelUpdates);
        
        % first check if there is an events.newTrial
        if any(strcmp(allNames, 'events.newTrial'))
          
          obj.Block.numCompletedTrials = obj.Block.numCompletedTrials+1;
          
          % Step 1: finish up the last trial
          obj.PsychometricAxes.clear();
          if obj.Block.numCompletedTrials > 2
            psy.plot2AUFC(obj.PsychometricAxes.Handle, obj.Block);
            contrast = diff(obj.Parameters.Struct.stimulusContrast)*100;
            obj.PsychometricAxes.plot(repmat(contrast,1,2),[0 100],'k:')
          end
          
          % make sure we have all necessary data about new trial
          assert(all(ismember(...
            {'events.trialNum', 'events.repeatNum'}, allNames)), ...
            'exp panel did not find all the required data about the new trial!');
          
          % pull out the things we need to keep
          trNum = updates(strcmp(allNames, 'events.trialNum')).value;
          % save the trial time
          %  trT = cellfun(@datenum,{updates(strcmp(allNames, 'events.trialNum')).timestamp});
          obj.Block.trial(end+1).newTrialTimes = datenum(updates(strcmp(allNames, 'events.trialNum')).timestamp);
          if ~(trNum==obj.Block.numCompletedTrials+1)
            fprintf(1, 'trial number mismatch: %d, %d\n', trNum, obj.Block.numCompletedTrials+1);
            obj.Block.numCompletedTrials = trNum-1;
          end
          obj.Block.trial(trNum).repeatNum = ...
            updates(strcmp(allNames, 'events.repeatNum')).value;
        end
        
        
        for ui = 1:length(updates)
          signame = updates(ui).name;
          switch signame
            case 'events.contrastLeft'
              % Store the left contrast value for use in the wheel and
              % psychometric plots.
              cL = updates(ui).value;
              idx = find(datenum(updates(ui).timestamp) > [obj.Block.trial.newTrialTimes], 1, 'last')+1;
              obj.Block.trial(idx).contrastLeft = cL;
              obj.contrastLeft = cL;
            case 'events.contrastRight'
              % Store the right contrast value for use in the wheel and
              % psychometric plots.
              cR = updates(ui).value;
              idx = find(datenum(updates(ui).timestamp) > [obj.Block.trial.newTrialTimes], 1, 'last')+1;
              obj.Block.trial(idx).contrastRight = cR;
              obj.contrastRight = cR;
            case 'events.wheelGain'
              % Store the wheel gain value for use with the wheel plot
              idx = find(datenum(updates(ui).timestamp) > [obj.Block.trial.newTrialTimes], 1, 'last')+1;
              obj.Block.trial(idx).wheelGain = updates(ui).value;
            case 'events.interactiveOn'
              % Start of interactive period means we must update the wheel
              % axes to show the thresholds and to colour the correct side
              % green.
              
              % re-set the response window starting now
              ioTime = (24*3600*datenum(updates(ui).timestamp))-(24*3600*obj.StartedDateTime);
              
              p = obj.Parameters.Struct;
              respWin = Inf; if respWin>1000; respWin = 1000; end
              % Due to the order of the updates, we look for the gain of
              % the previous trial.
              gain = getOr(obj.Block.trial(end-1), 'wheelGain', p.normalGain);
              th = p.responseDisplacement/gain;
              startPos = obj.InputSensorPos(find(obj.InputSensorPosTime<ioTime,1,'last'));
              if isempty(startPos); startPos = obj.InputSensorPos(obj.InputSensorPosCount); end % for first trial
              tL = startPos-th;
              tR = startPos+th;
              
              set(obj.ExperimentHands.threshL, ...
                'XData', [tL tL], 'YData', ioTime+[0 respWin]);
              set(obj.ExperimentHands.threshR, ...
                'XData', [tR tR], 'YData', ioTime+[0 respWin]);
              
              yd = get(obj.ExperimentHands.threshLoff, 'YData');
              set(obj.ExperimentHands.threshLoff, 'XData', [tL tL], 'YData', [yd(1) ioTime]);
              set(obj.ExperimentHands.threshRoff, 'XData', [tR tR], 'YData', [yd(1) ioTime]);
              
              cL = obj.contrastLeft;
              cR = obj.contrastRight;
              if ~isempty(cL)&&~isempty(cR)
                if cL>0 && cL>cR
                  colorL = 'g'; colorR = 'r';
                elseif cL>0 && cL==cR
                  colorL = 'g'; colorR = 'g';
                elseif cR>0
                  colorL = 'r'; colorR = 'g';
                elseif isnan(cL)||isnan(cR)
                  colorL = 'k'; colorR = 'k';
                else
                  colorL = 'r'; colorR = 'r';
                end
                set(obj.ExperimentHands.threshL, 'Color', colorL);
                set(obj.ExperimentHands.threshR, 'Color', colorR);
              end
              obj.ExperimentAxes.XLim = startPos+1.5*th*[-1 1];
              
            case 'events.stimulusOn'
              % Update the wheel axes to reflect the stimulus appearing
              p = obj.Parameters.Struct;
              soTime = (24*3600*datenum(updates(ui).timestamp))-(24*3600*obj.StartedDateTime);
              gain = iff(isempty(obj.Block.trial(end).wheelGain), p.normalGain, obj.Block.trial(end).wheelGain);
              th = p.responseDisplacement/gain;
              startPos = obj.InputSensorPos(find(obj.InputSensorPosTime<soTime,1,'last'));
              if isempty(startPos); startPos = obj.InputSensorPos(obj.InputSensorPosCount); end % for first trial
              tL = startPos-th;
              tR = startPos+th;
              
              set(obj.ExperimentHands.threshLoff,  ...
                'XData', [tL tL], 'YData', soTime+[0 100]);
              set(obj.ExperimentHands.threshRoff, ...
                'XData', [tR tR], 'YData', soTime+[0 100]);
              set(obj.ExperimentHands.threshL, 'YData', [NaN NaN]);
              set(obj.ExperimentHands.threshR, 'YData', [NaN NaN]);
              
              set(obj.ExperimentHands.incorrIcon, 'XData', 0, 'YData', NaN);
              set(obj.ExperimentHands.corrIcon, 'XData', 0, 'YData', NaN);
              
              obj.ExperimentAxes.XLim = startPos+1.5*th*[-1 1];
            case 'events.stimulusOff'
              % re-set the response window starting now
              ioTime = (24*3600*datenum(updates(ui).timestamp))-(24*3600*obj.StartedDateTime);
              yd = get(obj.ExperimentHands.threshR, 'YData');
              set(obj.ExperimentHands.threshL, ...
                'YData', [yd(1) ioTime]);
              set(obj.ExperimentHands.threshR, ...
                'YData', [yd(1) ioTime]);
              
            case 'events.response'
              % Store the response value for updating the psychometic curve
              obj.Block.trial(obj.Block.numCompletedTrials+1).response = updates(ui).value;
              
            case 'events.feedback'
              % Store the feedback value for updating the psychometic curve
              % and plot the feedback icon on the wheel axes plot. 
              obj.Block.trial(obj.Block.numCompletedTrials+1).feedback = updates(ui).value;
              
              fbTime = (24*3600*datenum(updates(ui).timestamp))-(24*3600*obj.StartedDateTime);
              whIdx = find(obj.InputSensorPosTime<fbTime,1, 'last');
              
              if updates(ui).value > 0
                set(obj.ExperimentHands.corrIcon, ...
                  'XData', obj.InputSensorPos(whIdx), ...
                  'YData', obj.InputSensorPosTime(whIdx));
                set(obj.ExperimentHands.incorrIcon, ...
                  'XData', 0, ...
                  'YData', NaN);
              else
                set(obj.ExperimentHands.incorrIcon, ...
                  'XData', obj.InputSensorPos(whIdx), ...
                  'YData', obj.InputSensorPosTime(whIdx));
                set(obj.ExperimentHands.corrIcon, ...
                  'XData', 0, ...
                  'YData', NaN);
              end
              
            case 'events.trialNum'
              % Update the trail count label
              set(obj.TrialCountLabel, ...
                'String', num2str(updates(ui).value));
              
            case {'events.baselineRT', 'events.windowedRT'}
              % Store the RT values for use in the disengagement plot
              name = signame(8:end);
              % value = str2double(strrep(updates(ui).value, ' sec', ''));
              value = updates(ui).value;
              obj.Block.trial(obj.Block.numCompletedTrials+1).(name) = value;
              
            case 'events.newTrial'
              % Update RT plot at the start of each new trial
              trialTimes = [obj.Block.trial.newTrialTimes];
              N = obj.WindowSize;
              if length(trialTimes) > N
                tMin = 1/(mean(diff(trialTimes(end-N+1:end)))*24*60);
                p = obj.Parameters.Struct;
                rtCriterion = p.rtCriterion;
                baseRT = [obj.Block.trial.baselineRT];
                winRT = [obj.Block.trial.windowedRT];
                set(obj.ScreenHands.baseRT, 'XData', 1:length(baseRT), 'YData', baseRT*rtCriterion)
                startIdx = numel(baseRT)+1-numel(winRT);
                set(obj.ScreenHands.winRT, 'XData', startIdx:numel(baseRT), 'YData', winRT)
                tMin = iff(all(isnan(obj.ScreenHands.tMin.YData)), ...
                  tMin, [obj.ScreenHands.tMin.YData tMin]);
                set(obj.ScreenHands.tMin, ...
                  'XData', obj.WindowSize:obj.WindowSize-1+length(tMin), ...
                  'YData', tMin)
                obj.ScreenAxes.YLim = [0 max(baseRT(startIdx:end)*rtCriterion)+0.1];
              end
              
            case 'outputs.reward'
              % Keep track of total reward and display both
              obj.Block.totalReward = obj.Block.totalReward + updates(ui).value;
              str = sprintf('%.1f%s (%.1f%s)', ...
                obj.Block.totalReward, obj.RewardUnits, ...
                updates(ui).value, obj.RewardUnits);
              set(obj.LabelsMap(signame), 'String', str, 'UserData', clock,...
                'ForegroundColor', obj.RecentColour);
              
            otherwise
              % For any custom updates that pass the UpdatesFilter, simply
              % display their values
              onList = any(ismember(signame, obj.UpdatesFilter));
              if (obj.Exclude && ~onList) || (~obj.Exclude && onList)
                if ~isKey(obj.LabelsMap, signame) % If new update, add field
                  obj.LabelsMap(signame) = obj.addInfoField(signame, '');
                end
                str = toStr(updates(ui).value); % Convert the value to string
                set(obj.LabelsMap(signame), 'String', str, 'UserData', clock,...
                  'ForegroundColor', obj.RecentColour); % Update value
              end
          end
        end
        
      end
    end
    
    function build(obj, parent)
      build@eui.SignalsExpPanel(obj, parent);
      
      % Build the psychometric axes
      plotgrid = uiextras.VBox('Parent', obj.CustomPanel, 'Padding', 5, 'Spacing', 10);
      
      %       uiextras.Empty('Parent', plotgrid, 'Visible', 'off');
      
      obj.PsychometricAxes = bui.Axes(plotgrid);
      obj.PsychometricAxes.YLim = [-1 101];
      obj.PsychometricAxes.NextPlot = 'add';
      obj.PsychometricAxes.Title = 'Pyschometric plot';
      obj.PsychometricAxes.xLabel('Contrast');
%       yyaxis(obj.PsychometricAxes.Handle,'left')
      set(obj.PsychometricAxes.Handle, 'YColor', [0 .8 .8])
      obj.PsychometricAxes.yLabel('% Leftward');
%       yyaxis(obj.PsychometricAxes.Handle,'right')
%       set(obj.PsychometricAxes.Handle, 'YColor', 'm')
%       obj.PsychometricAxes.yLabel('% Rightward');
      
      %       uiextras.Empty('Parent', plotgrid, 'Visible', 'off');
      
      obj.ScreenAxes = bui.Axes(plotgrid);
      obj.ScreenAxes.NextPlot = 'add';
      obj.ScreenAxes.Title = 'Reaction time';
      obj.ScreenAxes.xLabel('# trials');
      obj.ScreenAxes.yLabel('RT (s)');
      x = obj.Parameters.Struct.minTrials;
      plot(obj.ScreenAxes, [x, x], [0, 60], 'k:');
      obj.ScreenAxes.YLim = [0 5];
      obj.ScreenHands.winRT = plot(obj.ScreenAxes, nan(1,2), nan(1,2), 'k');
      obj.ScreenHands.baseRT = plot(obj.ScreenAxes, nan(1,2), nan(1,2), 'r');
      yyaxis(obj.ScreenAxes.Handle,'right')
      obj.ScreenHands.tMin = plot(obj.ScreenAxes, nan(1,2), nan(1,2), 'Color', [.8,.8,.8]);
      obj.ScreenAxes.yLabel('Trials / min');
      set(obj.ScreenAxes.Handle, 'YColor', [.5,.5,.5])
      yyaxis(obj.ScreenAxes.Handle,'left')
      scH = -2;
      
      %       uiextras.Empty('Parent', plotgrid, 'Visible', 'off');
      
      obj.ExperimentAxes = bui.Axes(plotgrid);
      obj.ExperimentAxes.ActivePositionProperty = 'position';
      obj.ExperimentAxes.XTickLabel = [];
      obj.ExperimentAxes.NextPlot = 'add';
      obj.ExperimentHands.wheelH = plot(obj.ExperimentAxes,...
        [0 0],...
        [NaN NaN],...
        'Color', .75*[1 1 1]);
      obj.ExperimentHands.threshL = plot(obj.ExperimentAxes, ...
        [0 0],...
        [NaN NaN],...
        'Color', [1 1 1], 'LineWidth', 4);
      obj.ExperimentHands.threshR = plot(obj.ExperimentAxes, ...
        [0 0],...
        [NaN NaN],...
        'Color', [1 1 1], 'LineWidth', 4);
      obj.ExperimentHands.threshLoff = plot(obj.ExperimentAxes, ...
        [0 0],...
        [NaN NaN],...
        'Color', [0.5 0.5 0.5], 'LineWidth', 4);
      obj.ExperimentHands.threshRoff = plot(obj.ExperimentAxes, ...
        [0 0],...
        [NaN NaN],...
        'Color', [0.5 0.5 0.5], 'LineWidth', 4);
      obj.ExperimentHands.corrIcon = scatter(obj.ExperimentAxes, ...
        0, NaN, pi*10^2, 'b', 'filled');
      obj.ExperimentHands.incorrIcon = scatter(obj.ExperimentAxes, ...
        0, NaN, pi*10^2, 'rx', 'LineWidth', 4);
      
      uiextras.Empty('Parent', plotgrid, 'Visible', 'off');
      
      obj.VelAxes = axes('Parent', plotgrid);
      obj.VelHands.Zero = plot(obj.VelAxes, [0 0], [0 1], 'k--');
      hold(obj.VelAxes, 'on');
      obj.VelHands.Vel = plot(obj.VelAxes, [0 0], [0 1], 'r', 'LineWidth', 2.0);
      axis(obj.VelAxes, 'off');
      obj.VelHands.MaxVel = 1e-9;
      
      %       set(plotgrid, 'Sizes', [30 -2 30 scH 10 -4 5 -1]);
      set(plotgrid, 'Sizes', [-2 scH -4 5 -1])
    end
  end
end