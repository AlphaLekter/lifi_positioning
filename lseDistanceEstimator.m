function distance = lseDistanceEstimator(LED, RIS, PDect_pos, R_pd, A_pd, p, m, T_of, G, received_current_array, rho,LED_or_RIS)
% Per il LED indicare LED_or_RIS == 1, 0 per le RIS.
if LED_or_RIS == 1
    %% LSE - LED distance
    h_LED = abs(LED(3)-PDect_pos(3)); % altezza
    %h0_LED = (R_pd/A_pd)*p *   (((m+1)*A_pd * (h_LED^(m+1) * T_of * G) )/ (2* pi));
    h0_LED = R_pd*p *   (((m+1)*A_pd * (h_LED^(m+1) * T_of * G) )/ (2* pi));
    degree = m+3;
    rad = h0_LED * (length(received_current_array) / sum(received_current_array));
    distance = nthroot(rad, degree);
elseif LED_or_RIS == 0
    %b = ((rho*((m+1)*A_pd* p)* cosd(calculateAngle(LED, RIS, 0, 0))^m * T_of * G * (RIS(3)-PDect_pos(3))* R_pd)/ (2*pi* A_pd));
    b = ((rho*((m+1)*A_pd* p)* cosd(calculateAngle(LED, RIS, 0, 0))^m * T_of * G * (RIS(3)-PDect_pos(3))* R_pd)/ (2*pi));
    b = b * (length(received_current_array) / sum(received_current_array));
    a = calculateDistance(LED, RIS);
    % soluzione REALE di (x+a)^2*x = b
    distance = 1/3 * ((3 * sqrt(3) * sqrt(4 * a^3 * b + 27 * b^2) + 2 * a^3 + 27 * b)^(1/3)/2^(1/3) + ...
        (2^(1/3) * a^2)/(3 * sqrt(3) * sqrt(4 * a^3 * b + 27 * b^2) + 2 * a^3 + 27 * b)^(1/3) - 2 * a);
else
    distance = NaN;
end

end
