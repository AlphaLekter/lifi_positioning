function [impulseResponse, number_of_contribution] = drawChannelResponse...
    (Psi, LED, RIS, RIS2, RIS3, RIS4, PDect_pos, Phi_FoV, A_pd, T_of, rho, a, Entity_enabled, alpha, beta )

    %(Psi, LED, RIS1, RIS2, RIS3, RIS4, PDect, Phi_FoV, A_pd, T_of, rho, a, Entity_enabled, alpha, beta );
% impulse response calc
%m = -(log(2)/log(cosd(Psi)));% Lambertian mode number


% 1. Estimate LED contribution
hLoS_sd = singleEntityContribution(LED, RIS, PDect_pos, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 1);
% 2. Estimate RIS1 contribution
hRIS_sd = singleEntityContribution(LED, RIS, PDect_pos, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 0);
% 3. Estimate RIS2 contribution
hRIS2_sd = singleEntityContribution(LED, RIS2, PDect_pos, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 0);
% 4. Estimate RIS3 contribution
hRIS3_sd = singleEntityContribution(LED, RIS3, PDect_pos, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 0);
% 5. Estimate RIS4 contribution
hRIS4_sd = singleEntityContribution(LED, RIS4, PDect_pos, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 0);

% Total channel response
overall_contribution = [hLoS_sd hRIS_sd hRIS2_sd hRIS3_sd hRIS4_sd]; % LED RIS RIS2 RIS3

% find out if is possible to estimate position
number_of_contribution = length(find(overall_contribution>0));
if number_of_contribution < 3
    number_of_contribution = 0;
    % else
    %     number_of_contribution = 1;
end
impulseResponse = sum(overall_contribution .* Entity_enabled );

end