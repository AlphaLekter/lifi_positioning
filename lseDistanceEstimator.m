function distance = lseDistanceEstimator(LED, PDect_pos, R_pd, A_pd, p, m, T_of, G, received_current_array)
    %% LSE - LED distance
    h_LED = abs(LED(3)-PDect_pos(3)); % altezza
    % h0_LED = (R_pd/A_pd) * p * (((m+1) * A_pd * (h_LED^(m+1) * T_of * G)) / (2*pi));
    h0_LED = R_pd * p * (((m+1) * A_pd * h_LED^(m+1) * T_of * G) / (2*pi));
    degree = m+3;
    rad = h0_LED * (length(received_current_array) / sum(received_current_array));
    % distance = nthroot(rad, degree);
    distance = nthroot(abs(rad), degree);
end
