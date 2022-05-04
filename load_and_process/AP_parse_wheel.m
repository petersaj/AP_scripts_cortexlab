function [wheel_velocity,wheel_move] = ...
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
% wheel_velocity: velocity of the wheel (as clicks/s)
% wheel_move: binary vector of times with movement or quiescence

% Set timing for smoothing wheel velocity
wheel_smooth_t = 0.05; % seconds
wheel_smooth_samples = round(wheel_smooth_t*sample_rate);

% Ensure wheel_position is column vector
wheel_position = wheel_position(:);

% Turn position into single clicks (position diff with leading 0)
wheel_clicks = [0;diff(wheel_position)]; % rotary encoder clicks

% Get rid of single wheel clicks within smoothing window
% (otherwise it's below meaninful velocity detection threshold)
wheel_clicks_use = movsum(abs(wheel_clicks),wheel_smooth_samples) > 1;
wheel_clicks_clean = wheel_clicks.*wheel_clicks_use;

% Get velocity: sum in window, median filter, convert to clicks/s
wheel_velocity = ...
    medfilt1(movsum(wheel_clicks_clean,wheel_smooth_samples), ...
    wheel_smooth_samples)/wheel_smooth_t;

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













