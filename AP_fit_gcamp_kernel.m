function sse = AP_fit_gcamp_kernel(params,t,actual_trace)
% sse = AP_fit_gcamp_kernel(params,t,actual_trace,t0_est,tau_on_est)
%
% Fit GCaMP kernel: single exponential up, double exponential down

t0 = params(1);
tau_on = params(2);
A1 = params(3);
tau1 = params(4);
A2 = params(5);
tau2 = params(6);

fitted_curve = (1-exp(-(t-t0)./tau_on)).*(A1*exp(-(t-t0)./tau1) + A2*exp(-(t-t0)./tau2));
fitted_curve(fitted_curve < 0) = 0;

error_vector = fitted_curve(t > t0) - actual_trace(t > t0);
sse = sum(error_vector.^2);