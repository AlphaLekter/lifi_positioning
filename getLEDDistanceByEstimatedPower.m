function distance = getLEDDistanceByEstimatedPower(h, p, A_pd, received_power, Psi, T_of, Phi_FoV,a, R_pd)

                                                 %(h_t, p, A_pd, received_power_LED  ,Psi, T_of, Phi_FoV,a, R_pd);
%abs(LED(3) - target(3)); % from [MCA+2022]
m = -(log(2)/log(cosd(Psi)));% Lambertian mode number
G = (a^2)/((sind(Phi_FoV)^2));
degree = m+3;
%rad = (((m+1)*h^(m+1)*A_pd*p)* T_of * G * R_pd)/ (2*pi*received_power*A_pd);
rad = (((m+1)*h^(m+1)*A_pd*p)* T_of * G * R_pd)/ (2*pi*received_power);

if rad > 0
    distance = nthroot(rad, degree);
elseif rad > 0 && mod(N, 2) == 1 % ??? Branch inutile?
    distance = nthroot(rad, degree);
else
    distance = NaN;
end

% if distance == 0
%     distance = NaN;
% end

end