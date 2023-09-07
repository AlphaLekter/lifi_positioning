function [n_shoot, n_thermal, var_shoot, var_thermal] = noiseEstimation(P_curr, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B)
% shoot noise
%S_pd = 0.54; %A/W

%P_rec = (A_pd/R_pd) * P_curr;
P_rec = P_curr;
% parameters from
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6880333

%var_shoot = 2*q_0*R_pd*P_rec*B + 2*q_0*I_bg*I_2*B;
var_shoot_first_term = 2*q_0*P_curr*B; % A * s * A * (1/s) % A^2
%var_shoot_first_term = 2*q_0*R_pd*P_rec*B; % A * s * A * (1/s) % A^2

%var_shoot_second_term = 2*q_0*I_bg*I_2*B; % I_BG
var_shoot_second_term = 2*q_0*(5e-12)*B; % caso I_DC = 5 pA

var_shoot = var_shoot_first_term + var_shoot_second_term;

% y = (Coulomb=A*s) * (A/W) * (W) * (1/s) + (Coulomb=A*s)* A * (1/s)
% y = [A^2] + [A^2] = [A^2]

% measurements_pickup = [randn(1) randn(1) randn(1) randn(1) randn(1)];
% var_vect = var_shoot_sqr*measurements_pickup;
% shoot = mean(var_vect);


shoot=randn(1)*sqrt(var_shoot);

% thermal noise
var_thermal = ((8*pi*k_B*T_k*eta*A_pd*I_2*B^2)/G_0 + ((4*pi)^2*k_B*T_k*Gamma*eta^2*A_pd^2*I_3*B^3)/g_m);
% y = y1+y2 = (Joule/Kelvin) * Kelvin * (F/m^2) * (m^2) * (Mb/s)^2 + (Joule/Kelvin) * Kelvin * (F^2/m^4) * m^4 * (Mb/s)^3 / Siemens
% y1 = [K] [F/m^2] [m^2] [1/s^2] = [K][s^4 A^2 / m^2 * kg] [1/s^2] [m^2 Kg/ s^2 K] = A^2
% y2 = [m^2 kg / s^2 K]*[K]*[1/mS] * [F^2/m^4] * [m^4] * [1/s^3] =
% = [m^2 kg / s^2 K]*[K]* [kg * m^2 / s^3 * A^2] * [s^8 * A^4 / m^4 *
% kg^2]*[1/s^3] = A^2


% measurements_pickup = [randn(1) randn(1) randn(1) randn(1) randn(1)];
% var_vect = var_thermal_sqr*measurements_pickup;
% thermal = mean(var_vect);

thermal =randn(1)*sqrt(var_thermal);


if P_curr > 0 
    n_shoot = shoot; 
    n_thermal = thermal;
    %noiseContribution = shoot+thermal;
else
    n_shoot = 0; 
    n_thermal = 0;
end

end
