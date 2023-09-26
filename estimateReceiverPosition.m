
function [x_eval, y_eval, z_eval] = estimateReceiverPosition(LED1, LED2, LED3, LED4, PDect_pos, p, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of, R_pd, q_0, k_B, T_k, eta, I_2, I_3, Gamma, g_m, I_bg, G_0, B, K_0)
% calcoli presi da [SAA+2022], %[c2022]
% parameters: LED1, LED2, LED3, LED4, PDect_pos, p, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of, R_pd, q_0, k_B, T_k, eta, I_2, I_3, Gamma, g_m, I_bg, G_0, B
    debug_enabled = 0;
    kalman_filter_enabled = 0;
    
    m = -(log(2)/log(cosd(Psi)));
    G = (a^2)/((sind(Phi_FoV)^2));
    
    %% Received power from each entity
    
    % kalman filter initialization
    if kalman_filter_enabled
        x_LED1 = 0;
        x_LED2 = 0;
        x_LED3 = 0;
        x_LED4 = 0;
        
        P_LED1 = 1e-9;
        P_LED2 = 1e-9;
        P_LED3 = 1e-9;
        P_LED4 = 1e-9;
    end
    
    received_current_LED1_array = nan(1, K_0);
    received_current_LED2_array = nan(1, K_0);
    received_current_LED3_array = nan(1, K_0);
    received_current_LED4_array = nan(1, K_0);
    
    %par
    for sampling_idx = 1:K_0
        % LED1 contribution    
        received_current_LED1_array(sampling_idx) = R_pd * p * singleEntityContribution(LED1, PDect_pos, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of);
        [n_shoot_LED1, n_thermal_LED1, ~, ~] = noiseEstimation(received_current_LED1_array(sampling_idx), q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
        received_current_LED1_array(sampling_idx) = received_current_LED1_array(sampling_idx) + n_shoot_LED1 + n_thermal_LED1;

        % LED2 contribution    
        received_current_LED2_array(sampling_idx) = R_pd * p * singleEntityContribution(LED2, PDect_pos, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of);
        [n_shoot_LED2, n_thermal_LED2, ~, ~] = noiseEstimation(received_current_LED2_array(sampling_idx), q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
        received_current_LED2_array(sampling_idx) = received_current_LED2_array(sampling_idx) + n_shoot_LED2 + n_thermal_LED2;
        
        % LED3 contribution    
        received_current_LED3_array(sampling_idx) = R_pd * p * singleEntityContribution(LED3, PDect_pos, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of);
        [n_shoot_LED3, n_thermal_LED3, ~, ~] = noiseEstimation(received_current_LED3_array(sampling_idx), q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
        received_current_LED3_array(sampling_idx) = received_current_LED3_array(sampling_idx) + n_shoot_LED3 + n_thermal_LED3;

        % LED4 contribution    
        received_current_LED4_array(sampling_idx) = R_pd * p * singleEntityContribution(LED4, PDect_pos, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of);
        [n_shoot_LED4, n_thermal_LED4, ~, ~] = noiseEstimation(received_current_LED4_array(sampling_idx), q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
        received_current_LED4_array(sampling_idx) = received_current_LED4_array(sampling_idx) + n_shoot_LED4 + n_thermal_LED4;
    end
    
    received_current_LED1 = sum(received_current_LED1_array)/K_0;
    received_current_LED2 = sum(received_current_LED2_array)/K_0;
    received_current_LED3 = sum(received_current_LED3_array)/K_0;
    received_current_LED4 = sum(received_current_LED4_array)/K_0;
    
    %% testONLY => forzo la correttezza delle distanze
    dLED1_real = calculateDistance(LED1, PDect_pos);
    dLED2_real = calculateDistance(LED2, PDect_pos);
    dLED3_real = calculateDistance(LED3, PDect_pos);
    dLED4_real = calculateDistance(LED4, PDect_pos);
    
    % d1 = dLED1_real;
    % d2 = dLED2_real;
    % d3 = dLED3_real;
    % d4 = dLED4_real;
    
    %% LSE method
    d1 = lseDistanceEstimator(LED1, PDect_pos, R_pd, A_pd, p, m, T_of, G, received_current_LED1_array);
    d2 = lseDistanceEstimator(LED2, PDect_pos, R_pd, A_pd, p, m, T_of, G, received_current_LED2_array);
    d3 = lseDistanceEstimator(LED3, PDect_pos, R_pd, A_pd, p, m, T_of, G, received_current_LED3_array);
    d4 = lseDistanceEstimator(LED4, PDect_pos, R_pd, A_pd, p, m, T_of, G, received_current_LED4_array);
    
    rec_power_vect = [
        received_current_LED1
        received_current_LED2
        received_current_LED3
        received_current_LED4
    ];
    distance_vect = [d1 d2 d3 d4];
    entity_vect(1).position = LED1;
    entity_vect(2).position = LED2;
    entity_vect(3).position = LED3;
    entity_vect(4).position = LED4;
    
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
