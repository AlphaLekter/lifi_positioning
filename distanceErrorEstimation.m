
function distanceErrorEstimation(LED, number_of_plot_needed, ...
    p, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of, R_pd, q_0, k_B, T_k, eta, I_2, I_3, Gamma, g_m, I_bg, G_0, B, ...
    x_probe, y_probe, z_plot_index, number_of_samples, led_number)

% m = -(log(2)/log(cosd(Psi))); % Lambertian mode number
% G = (a^2)/((sind(Phi_FoV)^2));
% we sample the photodiode current at a rate of 1000 samples/s for a period of 40 s
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6868970
% number_of_samples = 1000000;
received_power_LED_array = zeros(1, number_of_samples);

figure();
for k=1:number_of_plot_needed
    
    entityError = zeros(length(x_probe), length(y_probe));
    % error_matrix = [];

    for i=1:length(x_probe)
        for j=1:length(y_probe)
            PDect_pos_LED = [x_probe(i) y_probe(j) z_plot_index(k)];
            
            for sampling_idx = 1:number_of_samples
                % Estimate LED contribution
                % received_power_LED_array(sampling_idx) = (R_pd/A_pd)*p*singleEntityContribution(LED, PDect_pos_LED, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of);
                received_power_LED_array(sampling_idx) = R_pd * p * ...
                    singleEntityContribution(LED, PDect_pos_LED, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of);
                [n_shoot_LED, n_thermal_LED, ~, ~] = noiseEstimation(received_power_LED_array(sampling_idx), q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
                received_power_LED_array(sampling_idx) = received_power_LED_array(sampling_idx) + n_shoot_LED + n_thermal_LED;
            end

            received_power_LED(i,j) = sum(received_power_LED_array) / number_of_samples;
            
            dLED_real = calculateDistance(LED, PDect_pos_LED);
            h_t = LED(3) - PDect_pos_LED(3);
            
            if received_power_LED(i,j) >= 0
                d1 = getLEDDistanceByEstimatedPower(h_t, p, A_pd, received_power_LED(i,j), Psi, T_of, Phi_FoV, a, R_pd);
                entityError(i,j) = abs(dLED_real - d1);
            else
                % d1 = NaN;
                entityError(i,j) = NaN;
            end
            
        end
    end

    subplot(number_of_plot_needed/3, 3, k);
    
    % if isempty(entityError(~isnan(entityError) | entityError ~= Inf))
    if sum(sum(isnan(entityError) + (entityError == Inf))) > 0
        mean_error_string = "...";
    else
        mean_error_string = string(mean(mean(entityError)));
    end
    
    colormap;
    surf(entityError);
    xlim([0 length(x_probe)]);
    ylim([0 length(y_probe)]);
    xlabel('Y');
    ylabel('X');
    title('Z = '+string(z_plot_index(k))+' [m] - mean ' + mean_error_string);
    
end

colormap;
sgtitle('[LED' + led_number + ' distance error]');

% RIS Contribution (da ignorare)
% figure();
% for k=1:number_of_plot_needed
% 
%     entityError = zeros(length(x_probe), length(y_probe));
%     %error_matrix = [];
%     for i=1:length(x_probe)
%         for j=1:length(y_probe)
%             PDect_pos_RIS1 = [x_probe(i) y_probe(j) z_plot_index(k)];
% 
%             for sampling_idx = 1:number_of_samples
%                 % 2. Estimate RIS1 contribution
%                 %received_power_RIS1_array(sampling_idx) = (R_pd/A_pd)*p*singleEntityContribution(LED, PDect_pos_RIS1, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of);
%                 received_power_RIS1_array(sampling_idx) = R_pd*p*singleEntityContribution(LED, PDect_pos_RIS1, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of);
%                 [n_shoot_RIS1, n_thermal_RIS1, ~, ~] = noiseEstimation(received_power_RIS1_array(sampling_idx), q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
%                 received_power_RIS1_array(sampling_idx) = received_power_RIS1_array(sampling_idx) + n_shoot_RIS1 + n_thermal_RIS1;
%             end
%             received_power_RIS1(i,j) = sum(received_power_RIS1_array)/number_of_samples;
% 
%             %dRIS1_diff = calculateDistance(LED, RIS1);
%             dRIS1_real = calculateDistance(RIS1, PDect_pos_RIS1);
%             %d2 = getDistanceByEstimatedPower(h_r, p, A_pd, received_power_RIS1 ,m, T_of, G) - dRIS1_diff;
%             hRIS = RIS1(3)-PDect_pos_RIS1(3);
%             theta_RS = calculateAngle(LED, RIS1, 0, 0);
%             d_RS = calculateDistance(LED, RIS1);
%             if received_power_RIS1(i,j) > 0
%                 d2 = getRISDistanceByEstimatedPower(hRIS, p, A_pd, received_power_RIS1(i,j), Psi, T_of, Phi_FoV, a, rho, cosd(theta_RS), R_pd, d_RS);
%                 entityError(i,j)                 = abs(dRIS1_real - d2);
%             else
%                 %d2 = NaN;
%                 entityError(i,j)                 = NaN;
%             end
% 
%         end
%     end
%     %subplot(3,3,k);
%     subplot(number_of_plot_needed/3, 3, k);
% 
%     if isempty(entityError(~isnan(entityError) | entityError ~= Inf))
%         mean_error_string = "...";
%     else
%         mean_error_string = string(mean(mean(entityError)));
%     end
% 
%     colormap;
%     surf(entityError);
%     xlim([0 length(x_probe)]);
%     ylim([0 length(y_probe)]);
%     xlabel('Y');
%     ylabel('X');
%     title('Z = '+string((z_plot_index(k)))+' [m] - mean ' + mean_error_string);
% 
% end
% colormap;
% sgtitle('[RIS1 distance error]');

end
