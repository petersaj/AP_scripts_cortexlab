function AP_print_progress_fraction(numerator,denominator)
% AP_print_progress_fraction(numerator,denominator)
% 
% Display the progress of a for loop by fraction in the command line

% If previous print, delete line
if numerator > 1
    fprintf(repmat('\b',1,numel(num2str(numerator-1))+numel(num2str(denominator))+1));
end

% Print fraction
fprintf('%d/%d',numerator,denominator);

% If at the end, print a new line
if numerator == denominator
   fprintf('\n'); 
end