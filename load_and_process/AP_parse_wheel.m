function [wheel_velocity,wheel_move,wheel_velocity_split] = ...
    AP_parse_wheel(wheel_position,sample_rate)
% [wheel_velocity,wheel_move,wheel_velocity_split] = AP_parse_wheel(wheel_position,sample_rate)
% 
% Get velocity and parse movements from wheel position trace 
%
% Inputs: 
% wheel_position: assumes rotary encoder +/- 1 steps, evenly sampled
% sample_rate: sample rate for wheel position
%
% Outputs: 
% wheel_velocity: velocity of the wheel
% wheel_move: binary vector of times with movement or quiescence
% wheel_velocity_split: velocity of the wheel split into move/quies blocks

% Ensure wheel_position is column vector
wheel_position = wheel_position(:);

% Get wheel velocity
% (pick smoothing timing)
wheel_smooth_t = 0.05; % seconds
wheel_smooth_samples = round(wheel_smooth_t*sample_rate);

% (get rid of single wheel clicks or balancing clicks within the window)
wheel_clicks = diff(wheel_position); % rotary encoder clicks
wheel_click_remove = ...
    conv(abs(wheel_clicks),ones(1,wheel_smooth_samples),'same') < 2  | ...
    abs(conv(wheel_clicks,ones(1,wheel_smooth_samples),'same')) < 1;
wheel_clicks(wheel_click_remove) = 0;

% (get velocity by smoothing and median filtering cleaned trace)
wheel_velocity = interp1(conv(1:length(wheel_position),[1,1]/2,'valid'), ...
    medfilt1(smooth(wheel_clicks,wheel_smooth_samples),wheel_smooth_samples), ...
    1:length(wheel_position),'linear','extrap')';

% Threshold wheel for movement, get start/stops
wheel_velocity_thresh = abs(wheel_velocity) > 0;
% (if constantly moving or not moving: set constant and return)
if length(unique(wheel_velocity_thresh)) ~= 2
    wheel_move = wheel_velocity_thresh;
    return
end

wheel_starts_all = find(diff([false;wheel_velocity_thresh;false]) == 1);
wheel_stops_all = find(diff([false;wheel_velocity_thresh;false]) == -1);

% Combine movements with small gaps between
combine_move_t = 0.3; % in s (empirical/arbitrary)
combine_move_samples = round(combine_move_t*sample_rate);
combine_move = find((wheel_starts_all(2:end) - wheel_stops_all(1:end-1)) > combine_move_samples);
wheel_starts_trim = wheel_starts_all([1;combine_move+1]);
wheel_stops_trim = wheel_stops_all([combine_move;end]);

% Make vector of movement/quiescence
wheel_move = logical(interp1([wheel_starts_trim;wheel_stops_trim;0], ...
    [ones(size(wheel_starts_trim));zeros(size(wheel_stops_trim));0], ...
    transpose(1:length(wheel_position)),'previous','extrap'));

% Split wheel velocity into quiescence/move segments
move_quiescence_blocks = ...
    diff(unique([0;find(diff(wheel_move) ~= 0);length(wheel_move)]));
wheel_velocity_split = mat2cell(wheel_velocity,move_quiescence_blocks);













