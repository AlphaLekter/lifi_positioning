function [impulseResponse, number_of_contribution] = drawChannelResponse...
    (Psi, LED1, LED2, LED3, LED4, PDect_pos, Phi_FoV, A_pd, T_of, a, Entity_enabled, alpha, beta )
    % (Psi, LED1, LED2, LED3, LED4, PDect, Phi_FoV, A_pd, T_of, a, Entity_enabled, alpha, beta );
    
    % impulse response calc
    % m = -(log(2)/log(cosd(Psi))); % Lambertian mode number (??? da ignorare?)
    
    % 1. Estimate LED1 contribution
    hLED1_sd = singleEntityContribution(LED1, PDect_pos, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of);
    % 2. Estimate LED2 contribution
    hLED2_sd = singleEntityContribution(LED2, PDect_pos, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of);
    % 3. Estimate LED3 contribution
    hLED3_sd = singleEntityContribution(LED3, PDect_pos, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of);
    % 4. Estimate LED4 contribution
    hLED4_sd = singleEntityContribution(LED4, PDect_pos, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of);
    
    % Total channel response
    overall_contribution = [hLED1_sd hLED2_sd hLED3_sd hLED4_sd]; % LED1 LED2 LED3 LED4
    
    % find out if is possible to estimate position
    number_of_contribution = length(find(overall_contribution>0));
    
    if number_of_contribution < 3
        number_of_contribution = 0;
        % else
        %     number_of_contribution = 1;
    end
    
    impulseResponse = sum(overall_contribution .* Entity_enabled );

end