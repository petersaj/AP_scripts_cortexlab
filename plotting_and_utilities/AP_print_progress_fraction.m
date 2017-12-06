function AP_print_progress_fraction(numerator,denominator)
% AP_print_progress_fraction(numerator,denominator)
% 
% Display the progress of a for loop by fraction in the command line

if numerator > 1
    fprintf(repmat('\b',1,numel(num2str(numerator-1))+numel(num2str(denominator))+1));
end
fprintf('%d/%d',numerator,denominator);