function impulseResponse = singleEntityContribution...
    (LED, PDect_pos, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of)   
% parameters: (LED, PDect_pos, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of)   
m = -(log(2)/log(cosd(Psi))); % Lambertian mode number
% G = (a^2)/((sind(Phi_FoV)^2)); (??? da rimuovere)

theta_SD = calculateAngle(LED, PDect_pos, alpha, beta);
% TODO:to check
d_SD = calculateDistance(LED, PDect_pos); % distance LED to PD
% optical concentrator gain
if theta_SD >= 0 && theta_SD <= Phi_FoV
    G = (a^2)/((sind(Phi_FoV)^2));
else
    G = 0;
end

impulseResponse = (((m+1)*A_pd * ((cosd(theta_SD)^m) * cosd(theta_SD) * T_of * G) )/ (2* pi * (d_SD)^2)) ;
    
% RIS contribution (da ignorare)
%     theta_RS = calculateAngle(LED, RIS, alpha, beta); % angle of irradiance % LED to RIS
%     Phi_RD = calculateAngle(RIS, PDect_pos, alpha, beta); % angle of incidence % RIS to PD
% 
%     d_SR = calculateDistance(LED, RIS); % distance LED to RIS
%     %d_SR = 1;
%     d_RD = calculateDistance(RIS, PDect_pos); % distance RIS to PD
% 
%     % optical concentrator gain
%     if Phi_RD >= 0 && Phi_RD <= Phi_FoV
%         G = (a^2)/((sind(Phi_FoV)^2));
%     else
%         G = 0;
%     end
% 
%     impulseResponse = ((rho*(m+1)*A_pd * (cosd(theta_RS)^m) * cosd(Phi_RD) * T_of * G)/ (2* pi * (d_RD+d_SR)^2)) ;

end