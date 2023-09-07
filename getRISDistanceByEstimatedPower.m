function distance = getRISDistanceByEstimatedPower(h, p, A_pd, received_power, Psi, T_of, Phi_FoV,a, rho, cos_RIS, R_pd, d_RS)
m = -(log(2)/log(cosd(Psi)));% Lambertian mode number
G = (a^2)/((sind(Phi_FoV)^2));

%b = ((rho*((m+1)*A_pd* p)* cos_RIS^m * T_of * G * h* R_pd)/ (2*pi* received_power *  A_pd));
b = ((rho*((m+1)*A_pd* p)* cos_RIS^m * T_of * G * h* R_pd)/ (2*pi* received_power));
a = d_RS;
% soluzione REALE di (x+a)^2*x = b
distance = 1/3 * ((3 * sqrt(3) * sqrt(4 * a^3 * b + 27 * b^2) + 2 * a^3 + 27 * b)^(1/3) /...
    2^(1/3) + (2^(1/3) * a^2) / (3 * sqrt(3) * sqrt(4 * a^3 * b + 27 * b^2) + 2 * a^3 + 27 * b)^(1/3) - 2 * a);
% if distance == 0
%     distance = NaN;
% end
end