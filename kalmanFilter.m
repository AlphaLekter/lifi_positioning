function [x_hat_k, P_k] = kalmanFilter(x_hat_k_prec, P_prec, shoot_noise,thermal_noise, measured_power)
% x_hat_k_prec, P_prec, noise_var, measured_power
%% Time Update
x_hat_minus_k = x_hat_k_prec;
P_minus_k = P_prec ;

% Measurement Update
K_k = (P_minus_k) /(P_minus_k + thermal_noise + shoot_noise);
x_hat_k = (x_hat_minus_k) + (K_k * (measured_power - x_hat_minus_k));
P_k = (1 - K_k) * (P_minus_k);


%% References
% [1] Greg Welch, Gary Bishop, "An Introduction to the Kalman Filter", 
%    University of North Carolina at Chapel Hill Department of Computer Science, 2001
% [2] M.S.Grewal, A.P. Andrews, "Kalman Filtering - 
%    Theory and Practice Using MATLAB", Wiley, 2001
end