function AP_print_progress_time(numerator,denominator)
% AP_print_progress_time(numerator,denominator)
% 
% Display the progress of a for loop by estimated time

% Get time from last loop
loop_t = round(toc);

if numerator == 1
    loop_t = 0;
end

estimated_total_time = denominator*loop_t;
estimated_remaining_time = (denominator-numerator)*loop_t;

% Print fraction (new line, can't delete because don't know last size)
fprintf('\n%d/%d',estimated_remaining_time,estimated_total_time);

% If at the end, print a new line
if numerator == denominator
   fprintf('\n'); 
end

% Start the clock for the next loop
tic

