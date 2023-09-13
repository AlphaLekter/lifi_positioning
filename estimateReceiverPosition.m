
function [x_eval, y_eval, z_eval] = estimateReceiverPosition(LED1, LED2, LED3, LED4, PDect_pos, p, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of, R_pd, q_0, k_B, T_k, eta, I_2, I_3, Gamma, g_m, I_bg, G_0, B, K_0, K_n)
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TODO: implementare somma rumore
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %par
    for sampling_idx = 1:K_0
        % LED1 contribution    
        received_current_LED1_array(sampling_idx) = R_pd * p * singleEntityContribution(LED1, PDect_pos, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of);
        LED1_no_noise_sig(sampling_idx) = received_current_LED1_array(sampling_idx);
        [n_shoot_LED1, n_thermal_LED1, ~, ~] = noiseEstimation(received_current_LED1_array(sampling_idx), q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
        received_current_LED1_array(sampling_idx) = received_current_LED1_array(sampling_idx) + n_shoot_LED1 + n_thermal_LED1;

        % LED2 contribution    
        received_current_LED2_array(sampling_idx) = R_pd * p * singleEntityContribution(LED2, PDect_pos, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of);
        LED2_no_noise_sig(sampling_idx) = received_current_LED2_array(sampling_idx);
        [n_shoot_LED2, n_thermal_LED2, ~, ~] = noiseEstimation(received_current_LED2_array(sampling_idx), q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
        received_current_LED2_array(sampling_idx) = received_current_LED2_array(sampling_idx) + n_shoot_LED2 + n_thermal_LED2;
        
        % LED3 contribution    
        received_current_LED3_array(sampling_idx) = R_pd * p * singleEntityContribution(LED3, PDect_pos, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of);
        LED3_no_noise_sig(sampling_idx) = received_current_LED3_array(sampling_idx);
        [n_shoot_LED3, n_thermal_LED3, ~, ~] = noiseEstimation(received_current_LED3_array(sampling_idx), q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
        received_current_LED3_array(sampling_idx) = received_current_LED3_array(sampling_idx) + n_shoot_LED3 + n_thermal_LED3;

        % LED4 contribution    
        received_current_LED4_array(sampling_idx) = R_pd * p * singleEntityContribution(LED4, PDect_pos, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of);
        LED4_no_noise_sig(sampling_idx) = received_current_LED4_array(sampling_idx);
        [n_shoot_LED4, n_thermal_LED4, ~, ~] = noiseEstimation(received_current_LED4_array(sampling_idx), q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
        received_current_LED4_array(sampling_idx) = received_current_LED4_array(sampling_idx) + n_shoot_LED4 + n_thermal_LED4;
    
        %kalman filter
        if kalman_filter_enabled
            LED1_mean_sig(sampling_idx) = mean(received_current_LED1_array(~isnan(received_current_LED1_array)));
            LED1_kalman_approx(sampling_idx) = x_LED1;
            if sampling_idx == 1
                P_LED1 = received_current_LED1_array(sampling_idx);
            end
            [x_LED1, P_LED1] = kalmanFilter(x_LED1, P_LED1, shoot_var_LED1, thermal_var_LED1, received_current_LED1_array(sampling_idx));
    
            LED2_mean_sig(sampling_idx) = mean(received_current_LED2_array(~isnan(received_current_LED2_array)));
            LED2_kalman_approx(sampling_idx) = x_LED2;
            if sampling_idx == 1
                P_LED2 = received_current_LED2_array(sampling_idx);
            end
            [x_LED2, P_LED2] = kalmanFilter(x_LED2, P_LED2, shoot_var_LED2, thermal_var_LED2, received_current_LED2_array(sampling_idx));
    
            LED3_mean_sig(sampling_idx) = mean(received_current_LED3_array(~isnan(received_current_LED3_array)));
            LED3_kalman_approx(sampling_idx) = x_LED3;
            if sampling_idx == 1
                P_LED3 = received_current_LED3_array(sampling_idx);
            end
            [x_LED3, P_LED3] = kalmanFilter(x_LED3, P_LED3, shoot_var_LED3, thermal_var_LED3, received_current_LED3_array(sampling_idx));
    
            LED4_mean_sig(sampling_idx) = mean(received_current_LED4_array(~isnan(received_current_LED4_array)));
            LED4_kalman_approx(sampling_idx) = x_LED4;
            if sampling_idx == 1
                P_LED4 = received_current_LED4_array(sampling_idx);
            end
            [x_LED4, P_LED4] = kalmanFilter(x_LED4, P_LED4, shoot_var_LED4, thermal_var_LED4, received_current_LED4_array(sampling_idx));
        end
    end
    
    received_current_LED1 = sum(received_current_LED1_array)/K_0;
    received_current_LED2 = sum(received_current_LED2_array)/K_0;
    received_current_LED3 = sum(received_current_LED3_array)/K_0;
    received_current_LED4 = sum(received_current_LED4_array)/K_0;
    
    if debug_enabled
        figure();
    
        subplot(2,2,1);
        hold on;
        plot(received_current_LED1_array);
        plot(LED1_mean_sig,'LineWidth',3);
        plot(LED1_kalman_approx,'LineWidth',3);
        plot(LED1_no_noise_sig);
        legend("Live measurement", "Mean", "Kalman Filter", "Real");
        title("LED1 Measurement");
        xlabel('Number of measurement [#]');
        ylabel('Received power [A^2]');
        
        subplot(2,2,2);
        hold on;
        plot(received_current_LED2_array);
        plot(LED2_mean_sig,'LineWidth',3);
        plot(LED2_kalman_approx,'LineWidth',3);
        plot(LED2_no_noise_sig);
        legend("Live measurement", "Mean", "Kalman Filter", "Real");
        title("LED2 Measurement");
        xlabel('Number of measurement [#]');
        ylabel('Received power [A^2]');
    
        subplot(2,2,3);
        hold on;
        plot(received_current_LED3_array);
        plot(LED3_mean_sig,'LineWidth',3);
        plot(LED3_kalman_approx,'LineWidth',3);
        plot(LED3_no_noise_sig);
        legend("Live measurement", "Mean", "Kalman Filter", "Real");
        title("LED3 Measurement");
        xlabel('Number of measurement [#]');
        ylabel('Received power [A^2]');
    
        subplot(2,2,4);
        hold on;
        plot(received_current_LED4_array);
        plot(LED4_mean_sig,'LineWidth',3);
        plot(LED4_kalman_approx,'LineWidth',3);
        plot(LED4_no_noise_sig);
        legend("Live measurement", "Mean", "Kalman Filter", "Real");
        title("LED4 Measurement");
        xlabel('Number of measurement [#]');
        ylabel('Received power [A^2]');
    end
    
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
