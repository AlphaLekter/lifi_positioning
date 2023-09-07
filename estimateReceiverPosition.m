
function [x_eval, y_eval, z_eval] = estimateReceiverPosition(LED, RIS1, RIS2, RIS3, RIS4, PDect_pos, p, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, R_pd, q_0, k_B, T_k, eta, I_2, I_3, Gamma, g_m, I_bg, G_0, B, K_0, K_n)
% calcoli presi da [SAA+2022], %[c2022]
% parameters: LED, RIS1, RIS2, RIS3, RIS4, PDect_pos, p, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, R_pd, q_0, k_B, T_k, eta, I_2, I_3, Gamma, g_m, I_bg, G_0, B
debug_enabled = 0;

kalman_filter_enabled = 0;

m = -(log(2)/log(cosd(Psi)));
G = (a^2)/((sind(Phi_FoV)^2));

%% Received power from each entity

% kalman filter initialization
if kalman_filter_enabled
    x_LED = 0;
    x_RIS1 = 0;
    x_RIS2 = 0;
    x_RIS3 = 0;
    x_RIS4 = 0;
    
    
    P_LED = 1e-8;
    P_RIS1 = 1e-9;
    P_RIS2 = 1e-9;
    P_RIS3 = 1e-9;
    P_RIS4 = 1e-9;
end

received_current_LED_array = nan(1, K_0);
received_current_RIS1_array = nan(1, K_n);
received_current_RIS2_array = nan(1, K_n);
received_current_RIS3_array = nan(1, K_n);
received_current_RIS4_array = nan(1, K_n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: implementare somma rumore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%par
for sampling_idx = 1:K_0
    % LED contribution    
    received_current_LED_array(sampling_idx) =  R_pd*p*singleEntityContribution(LED, 0, PDect_pos, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 1);
    LED_no_noise_sig(sampling_idx) = received_current_LED_array(sampling_idx);
    [n_shoot_LED, n_thermal_LED, ~, ~] = noiseEstimation(received_current_LED_array(sampling_idx), q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
    received_current_LED_array(sampling_idx) = received_current_LED_array(sampling_idx) + n_shoot_LED + n_thermal_LED;
end

for sampling_idx = 1:K_n
    
    % RIS1 contribution    
    received_current_RIS1_array(sampling_idx) = R_pd*p*singleEntityContribution(LED, RIS1, PDect_pos, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 0);
    received_current_RIS1_LOS = R_pd*p*singleEntityContribution(LED, 0, PDect_pos, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 1); % Ampere
    RIS1_no_noise_sig(sampling_idx) = received_current_RIS1_array(sampling_idx);
    % LoS noise 
    [n_shoot_LoS, n_thermal_LoS, ~, ~] = ...
        noiseEstimation(received_current_RIS1_LOS, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
    % NLoS noise
    [n_shoot_RIS1, n_thermal_RIS1, ~, ~] = ...
        noiseEstimation(received_current_RIS1_LOS + received_current_RIS1_array(sampling_idx), q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
    
    received_current_RIS1_array(sampling_idx) = received_current_RIS1_array(sampling_idx) + n_shoot_RIS1 + n_thermal_RIS1 + n_shoot_LoS + n_thermal_LoS;
    
    
    % RIS2 contribution    
    received_current_RIS2_array(sampling_idx) = R_pd*p*singleEntityContribution(LED, RIS2, PDect_pos, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 0);
    received_current_RIS2_LOS = R_pd*p*singleEntityContribution(LED, 0, PDect_pos, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 1); % Ampere
    RIS2_no_noise_sig(sampling_idx) = received_current_RIS2_array(sampling_idx);
    % LoS noise 
    [n_shoot_LoS, n_thermal_LoS, ~, ~] = ...
        noiseEstimation(received_current_RIS2_LOS, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
    % NLoS noise
    [n_shoot_RIS2, n_thermal_RIS2, ~, ~] = ...
        noiseEstimation(received_current_RIS2_LOS + received_current_RIS2_array(sampling_idx), q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
    
    received_current_RIS2_array(sampling_idx) = received_current_RIS2_array(sampling_idx) + n_shoot_RIS2 + n_thermal_RIS2 + n_shoot_LoS + n_thermal_LoS;
    
    
    % RIS3 contribution    
    received_current_RIS3_array(sampling_idx) = R_pd*p*singleEntityContribution(LED, RIS3, PDect_pos, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 0);
    received_current_RIS3_LOS = R_pd*p*singleEntityContribution(LED, 0, PDect_pos, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 1); % Ampere
    RIS3_no_noise_sig(sampling_idx) = received_current_RIS3_array(sampling_idx);
    % LoS noise 
    [n_shoot_LoS, n_thermal_LoS, ~, ~] = ...
        noiseEstimation(received_current_RIS3_LOS, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
    % NLoS noise
    [n_shoot_RIS3, n_thermal_RIS3, ~, ~] = ...
        noiseEstimation(received_current_RIS3_LOS + received_current_RIS3_array(sampling_idx), q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
    
    received_current_RIS3_array(sampling_idx) = received_current_RIS3_array(sampling_idx) + n_shoot_RIS3 + n_thermal_RIS3 + n_shoot_LoS + n_thermal_LoS;
        
    % RIS4 contribution    
    received_current_RIS4_array(sampling_idx) = R_pd*p*singleEntityContribution(LED, RIS4, PDect_pos, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 0);
    received_current_RIS4_LOS = R_pd*p*singleEntityContribution(LED, 0, PDect_pos, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 1); % Ampere
    RIS4_no_noise_sig(sampling_idx) = received_current_RIS4_array(sampling_idx);
    % LoS noise 
    [n_shoot_LoS, n_thermal_LoS, ~, ~] = ...
        noiseEstimation(received_current_RIS4_LOS, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
    % NLoS noise
    [n_shoot_RIS4, n_thermal_RIS4, ~, ~] = ...
        noiseEstimation(received_current_RIS4_LOS + received_current_RIS4_array(sampling_idx), q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
    
    received_current_RIS4_array(sampling_idx) = received_current_RIS4_array(sampling_idx) + n_shoot_RIS4 + n_thermal_RIS4 + n_shoot_LoS + n_thermal_LoS;
end

%     %kalman filter
%     if kalman_filter_enabled
%         
%         LED_mean_sig(sampling_idx) = mean(received_current_LED_array(~isnan(received_current_LED_array)));
%         LED_kalman_approx(sampling_idx) = x_LED;
%         if sampling_idx == 1
%             P_LED = received_current_LED_array(sampling_idx);
%         end
%         [x_LED, P_LED] = kalmanFilter(x_LED, P_LED, shoot_var_LED,thermal_var_LED, received_current_LED_array(sampling_idx));
%         
%         
%         RIS1_mean_sig(sampling_idx) = mean(received_current_RIS1_array(~isnan(received_current_RIS1_array)));
%         RIS1_kalman_approx(sampling_idx) = x_RIS1;
%         if sampling_idx == 1
%             P_RIS1 = received_current_RIS1_array(sampling_idx);
%         end
%         [x_RIS1, P_RIS1] = kalmanFilter(x_RIS1, P_RIS1, shoot_var_RIS1,thermal_var_RIS1, received_current_RIS1_array(sampling_idx));
%         
%         
%         RIS2_mean_sig(sampling_idx) = mean(received_current_RIS2_array(~isnan(received_current_RIS2_array)));
%         RIS2_kalman_approx(sampling_idx) = x_RIS2;
%         if sampling_idx == 1
%             P_RIS2 = received_current_RIS2_array(sampling_idx);
%         end
%         [x_RIS2, P_RIS2] = kalmanFilter(x_RIS2, P_RIS2, shoot_var_RIS2,thermal_var_RIS2, received_current_RIS2_array(sampling_idx));
%         
%         
%         RIS3_mean_sig(sampling_idx) = mean(received_current_RIS3_array(~isnan(received_current_RIS3_array)));
%         RIS3_kalman_approx(sampling_idx) = x_RIS3;
%         if sampling_idx == 1
%             P_RIS3 = received_current_RIS3_array(sampling_idx);
%         end
%         [x_RIS3, P_RIS3] = kalmanFilter(x_RIS3, P_RIS3, shoot_var_RIS3,thermal_var_RIS3, received_current_RIS3_array(sampling_idx));
%         
%         
%         RIS4_mean_sig(sampling_idx) = mean(received_current_RIS4_array(~isnan(received_current_RIS4_array)));
%         RIS4_kalman_approx(sampling_idx) = x_RIS4;
%         if sampling_idx == 1
%             P_RIS4 = received_current_RIS4_array(sampling_idx);
%         end
%         [x_RIS4, P_RIS4] = kalmanFilter(x_RIS4, P_RIS4, shoot_var_RIS4,thermal_var_RIS4, received_current_RIS4_array(sampling_idx));
%         
%     end

received_current_LED = sum(received_current_LED_array)/K_0;
received_current_RIS1 = sum(received_current_RIS1_array)/K_n;
received_current_RIS2 = sum(received_current_RIS2_array)/K_n;
received_current_RIS3 = sum(received_current_RIS3_array)/K_n;
received_current_RIS4 = sum(received_current_RIS4_array)/K_n;


if debug_enabled
    figure();
    subplot(3,2,1);
    hold on;
    plot(received_current_LED_array);
    plot(LED_mean_sig,'LineWidth',3);
    plot(LED_kalman_approx,'LineWidth',3);
    plot(LED_no_noise_sig);
    legend("Live measurement", "Mean", "Kalman Filter", "Real");
    title("LED Measurement");
    xlabel('Number of measurement [#]');
    ylabel('Received power [A^2]');
    
    subplot(3,2,2);
    hold on;
    plot(received_current_RIS1_array);
    plot(RIS1_mean_sig,'LineWidth',3);
    plot(RIS1_kalman_approx,'LineWidth',3);
    plot(RIS1_no_noise_sig);
    legend("Live measurement", "Mean", "Kalman Filter", "Real");
    title("RIS1 Measurement");
    xlabel('Number of measurement [#]');
    ylabel('Received power [A^2]');
    
    subplot(3,2,3);
    hold on;
    plot(received_current_RIS2_array);
    plot(RIS2_mean_sig,'LineWidth',3);
    plot(RIS2_kalman_approx,'LineWidth',3);
    plot(RIS2_no_noise_sig);
    legend("Live measurement", "Mean", "Kalman Filter", "Real");
    title("RIS2 Measurement");
    xlabel('Number of measurement [#]');
    ylabel('Received power [A^2]');
    
    subplot(3,2,4);
    hold on;
    plot(received_current_RIS3_array);
    plot(RIS3_mean_sig,'LineWidth',3);
    plot(RIS3_kalman_approx,'LineWidth',3);
    plot(RIS3_no_noise_sig);
    legend("Live measurement", "Mean", "Kalman Filter", "Real");
    title("RIS3 Measurement");
    xlabel('Number of measurement [#]');
    ylabel('Received power [A^2]');
    
    
    subplot(3,2,5.5);
    hold on;
    plot(received_current_RIS4_array);
    plot(RIS4_mean_sig,'LineWidth',3);
    plot(RIS4_kalman_approx,'LineWidth',3);
    plot(RIS4_no_noise_sig);
    legend("Live measurement", "Mean", "Kalman Filter", "Real");
    title("RIS4 Measurement");
    xlabel('Number of measurement [#]');
    ylabel('Received power [A^2]');
end



%% testONLY => forzo la correttezza delle distanze
dLED_real = calculateDistance(LED, PDect_pos);
dRIS1_real = calculateDistance(RIS1, PDect_pos);
dRIS2_real = calculateDistance(RIS2, PDect_pos);
dRIS3_real = calculateDistance(RIS3, PDect_pos);
dRIS4_real = calculateDistance(RIS4, PDect_pos);
% d1 = dLED_real;
% d2 = dRIS1_real;
% d3 = dRIS2_real;
% d4 = dRIS3_real;
% %d4 = NaN;
% d5 = dRIS4_real;
% d5 = NaN;
% %
%% LSE method
d1 = lseDistanceEstimator(LED, '' , PDect_pos, R_pd, A_pd, p, m, T_of, G, received_current_LED_array, '', 1);
d2 = lseDistanceEstimator(LED, RIS1 ,PDect_pos, R_pd, A_pd, p, m, T_of, G, received_current_RIS1_array, rho,0);
d3 = lseDistanceEstimator(LED, RIS2 ,PDect_pos, R_pd, A_pd, p, m, T_of, G, received_current_RIS2_array, rho,0);
d4 = lseDistanceEstimator(LED, RIS3 ,PDect_pos, R_pd, A_pd, p, m, T_of, G, received_current_RIS3_array, rho,0);
d5 = lseDistanceEstimator(LED, RIS4 ,PDect_pos, R_pd, A_pd, p, m, T_of, G, received_current_RIS4_array, rho,0);

%% kalman filter distance evaluation
if kalman_filter_enabled
    d1_est = getLEDDistanceByEstimatedPower(abs(PDect_pos(3)-LED(3)) , p, A_pd, x_LED,  Psi, T_of, Phi_FoV,a,               R_pd);
    d2_est = getRISDistanceByEstimatedPower(abs(PDect_pos(3)-RIS1(3)), p, A_pd, x_RIS1, Psi, T_of, Phi_FoV,a, rho, cosd(calculateAngle(LED, RIS1, 0, 0)), R_pd, calculateDistance(LED, RIS1));
    d3_est = getRISDistanceByEstimatedPower(abs(PDect_pos(3)-RIS2(3)), p, A_pd, x_RIS2, Psi, T_of, Phi_FoV,a, rho, cosd(calculateAngle(LED, RIS2, 0, 0)), R_pd, calculateDistance(LED, RIS2));
    d4_est = getRISDistanceByEstimatedPower(abs(PDect_pos(3)-RIS3(3)), p, A_pd, x_RIS3, Psi, T_of, Phi_FoV,a, rho, cosd(calculateAngle(LED, RIS3, 0, 0)), R_pd, calculateDistance(LED, RIS3));
    d5_est = getRISDistanceByEstimatedPower(abs(PDect_pos(3)-RIS4(3)), p, A_pd, x_RIS4, Psi, T_of, Phi_FoV,a, rho, cosd(calculateAngle(LED, RIS4, 0, 0)), R_pd, calculateDistance(LED, RIS4));
    
    d1 = d1_est;
    d2 = d2_est;
    d3 = d3_est;
    d4 = d4_est;
    d5 = d5_est;
end

%%
rec_power_vect = [received_current_LED received_current_RIS1 received_current_RIS2 received_current_RIS3 received_current_RIS4];
distance_vect = [d1 d2 d3 d4 d5];
entity_vect(1).position = LED;
entity_vect(2).position = RIS1;
entity_vect(3).position = RIS2;
entity_vect(4).position = RIS3;
entity_vect(5).position = RIS4;

[~, orderedIdx] = sort(rec_power_vect, 'descend');
distance_vect = distance_vect(orderedIdx);
entity_vect = entity_vect(orderedIdx);

distance_vect = distance_vect(~isnan(distance_vect) & ~isinf(distance_vect));
entity_vect = entity_vect(~isnan(distance_vect) & ~isinf(distance_vect));

pos_estimated = leastSquareMethod(entity_vect, distance_vect, PDect_pos(3));

if isnan(pos_estimated)
    x_eval = NaN;
    y_eval = NaN;
    z_eval = NaN;
else
    x_eval = pos_estimated(1);
    y_eval = pos_estimated(2);
    z_eval = pos_estimated(3);
end
end
