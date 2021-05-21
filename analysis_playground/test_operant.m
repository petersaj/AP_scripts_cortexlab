%% Get average rewarded and other movements

[wheel_vel,wheel_move,wheel_vel_split] = AP_parse_wheel(wheel_position,Timeline.hw.daqSampleRate);
wheel_pos_split = mat2cell(wheel_position,wheel_split_n);

wheel_split_n = cellfun(@length,wheel_vel_split)';
split_move = cellfun(@(x) x(1),mat2cell(wheel_move,wheel_split_n));

max_move = max(wheel_split_n(split_move == 1));
vel_pad = cell2mat(cellfun(@(x) padarray(x,max_move-length(x),NaN,'post'), ...
    wheel_vel_split(split_move == 1),'uni',false)');
pos_pad_raw = cell2mat(cellfun(@(x) padarray(x,max_move-length(x),NaN,'post'), ...
    wheel_pos_split(split_move == 1),'uni',false)'); 
pos_pad = pos_pad_raw - pos_pad_raw(1,:);


rewarded_epochs = cellfun(@(t) any(ismember(t,reward_t_timeline)), ...
    mat2cell(Timeline.rawDAQTimestamps',wheel_split_n));

rewarded_movements = rewarded_epochs(split_move == 1);

% Plot average rewarded and non-rewarded movements
figure; hold on
plot(nanmean(vel_pad(:,~rewarded_movements),2),'k');
plot(nanmean(vel_pad(:,rewarded_movements),2),'b');



