%% clear environment
clc; clear all; close all;
%% set entities position
% room size
x_max       = 5;            % room size x-axis         % [SCA+2022]
y_max       = 5;            % room size y-axis         % [SCA+2022]
z_max       = 3;            % room size z-axis         % 3 in [SCA+2022]
granularity = 0.1;          % plot accuracy

x1 = x_max/2; y1 = y_max/2; z1 = z_max;
x2 = x_max/2; y2 = 0; z2 = z_max/2+0.3*z_max;
x3 = 0; y3 = y_max/2; z3 = z_max/2+0.3*z_max;
x4 = x_max; y4 = y_max/2; z4 = z_max/2+0.3*z_max;
x5 = x_max/2; y5 = y_max; z5 = z_max/2+0.3*z_max;

% % position
% LED     = [x_max/2     , y_max/2      , z_max  ];
% RIS1    = [x_max/2      , 0            , z_max/2+0.3*z_max];
% RIS2    = [0            , y_max/2      , z_max/2+0.3*z_max];
% RIS3    = [x_max        , y_max/2      , z_max/2+0.3*z_max];
% RIS4    = [x_max/2      , y_max        , z_max/2+0.3*z_max];

% coordinate sonda
noise = 0;

x_probe = 0:granularity:x_max;
y_probe = 0:granularity:y_max;
z_probe = 0:9*granularity:z_max;
%fisso zeta
%xk = 3;
%yk = 3;
zk = 0;


error_matrix = zeros(length(x_probe), length(y_probe));

for i=1:length(x_probe)
    for j=1:length(y_probe)
        xk = x_probe(i);
        yk = y_probe(j);
        SONDA = [xk yk zk];
        
        %PD_to_findx = [x_probe(i) y_probe(j) z_probe(1)];
        %[Pos_estimated(1), Pos_estimated(2), Pos_estimated(3)] = estimateReceiverPosition(LED, RIS1, RIS2, RIS3, RIS4, PD_to_findx, ...
        %    p, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, R_pd, q_0, k_B, T_k, eta, I_2, I_3, Gamma, g_m, I_bg, G_0, B);
        
        
        d1 = sqrt((x1-xk)^2 + (y1-yk)^2 + (z1-zk)^2)+noise;
        d2 = sqrt((x2-xk)^2 + (y2-yk)^2 + (z2-zk)^2)+noise;
        d3 = sqrt((x3-xk)^2 + (y3-yk)^2 + (z3-zk)^2)+noise;
        d4 = sqrt((x4-xk)^2 + (y4-yk)^2 + (z4-zk)^2)+noise;
        d5 = sqrt((x5-xk)^2 + (y5-yk)^2 + (z5-zk)^2)+noise;
        
        %LEAST SQUARE
        A = 2* [...
            x2-x1, y2-y1 ,z2-z1;    ...
            x3-x1, y3-y1 ,z3-z1;    ...
            x4-x1, y4-y1 ,z4-z1;    ...
            x5-x1, y5-y1 ,z5-z1    ...
            ];
        
        
        b = [...
            d1^2-d2^2-(x1^2+y1^2+z1^2) + (x2^2+y2^2+z2^2);...
            d1^2-d3^2-(x1^2+y1^2+z1^2) + (x3^2+y3^2+z3^2);...
            d1^2-d4^2-(x1^2+y1^2+z1^2) + (x4^2+y4^2+z4^2);...
            d1^2-d5^2-(x1^2+y1^2+z1^2) + (x5^2+y5^2+z5^2)...
            ];
        
        Pos_estimated = pinv(A)*b;
        
        if isnan(Pos_estimated(1)) || isnan(Pos_estimated(2))
            error_matrix(i,j) = NaN;
        else
            error_matrix(i,j) = (sqrt( (Pos_estimated(1)-SONDA(1))^2 + (Pos_estimated(2)-SONDA(2))^2 ));
        end
        
    end
end
%toc
surf(error_matrix);
hold on;
colormap;


xlabel('Y');
ylabel('X');
title('Zoom su errore per terra con Z = '+string((z_probe(1)))+' m');




