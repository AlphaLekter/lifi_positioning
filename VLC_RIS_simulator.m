%%
% possibile parametro di ottimizzazione: distanza fra le RIS?

%% TODO
% 2. aggiustare yticks grafici boxplot

%% clear environment
clc; clear; close all;

%% parameters setting
Psi     = 70;               % LED half-power semiangle [degree]
% rho     = 0.95;             % Reflection coefficient (da ignorare)
% A_pd    = 1e-04;            % 1cm^2 - Physical area of the PD [m^2]  %  [C2022]
A_pd    = 0.2e-04;          % Physical area of the PD [m^2]
T_of    = 1;                % Optical Filter Gain
a       = 1.5;              % Refractive index
Phi_FoV = 70;               % Field of view [degree]
B       = 5e6;              % System bandwidth [Hz]
% B       = 100e6;            % System bandwidth [Hz]
% R_pd    = 2.2e-8;           % Responsivity [A][m^2] / [W] -> https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6868970
R_pd    = 0.54;             % (Sensitivity) [A/W]
q       = 3;                % Conversion ration of optical-to-electrical power
N       = 10e-21;           % Power spectral density [A^2/Hz]
lumen_level = 1000;
p = lumen_level / 683;      % transmission power 6000 [Lumens] -> [Watt]

% noise parameter
q_0     = 1.602e-19;        % electronic charge [Coulombs]
I_bg    = 84e-6;            % background light current [A] -> 84 [µA]
%I_bg    = 5e-12;           % 5pA
k_B     = 1.38064852e-23;   % Boltzmann constant [Joule/Kelvin]
T_k     = 295;              % absolute temperature [K]
G_0     = 10;               % open-loop voltage gain
eta     = 1.12e-6;          % fixed capacitance of photo detector per unit area [F/m^2]
Gamma   = 1.5;              % FET channel noise factor
g_m     = 0.030;            % FET transconductance [Siemens] [mS]
I_2     = 0.562;            % noise BW factor
I_3     = 0.0868;           % noise BW factor

% room size
x_max       = 5;            % room size x-axis         % [SCA+2022]
y_max       = 5;            % room size y-axis         % [SCA+2022]
z_max       = 3;            % room size z-axis         % 3 in [SCA+2022]

% room size
%  x_max       = 7;            % room size x-axis         % [SCA+2022]
%  y_max       = 7;            % room size y-axis         % [SCA+2022]
%  z_max       = 5;            % room size z-axis         % 3 in [SCA+2022]

granularity = .05;          % plot accuracy

number_of_samples = 100; % ??? Valori consigliati?
% plot info
Entity_enabled = [1 1 1 1]; % LED1 - LED2 - LED3 - LED4
PLOT3D_enabled = 0;
PLOT2D_enabled = 1;
plot_1         = 1;
plot_2         = 0;
plot_3         = 0;
plot_4         = 0;
plot_5         = 0;
plot_6         = 0;
plot_7         = 0;
number_of_plot_needed = 18;

% tilt info
alpha   = 0;
beta    = 0;

%% set entities position

LED1 = [x_max/4      , y_max/4       , z_max];
LED2 = [x_max/4      , y_max*3/4     , z_max];
LED3 = [x_max*3/4    , y_max/4       , z_max];
LED4 = [x_max*3/4    , y_max*3/4     , z_max];

%% Estimate light power
x_probe = 0:granularity:x_max;
y_probe = 0:granularity:y_max;
z_probe = 0:granularity:z_max;

x_len = length(x_probe);
y_len = length(y_probe);
z_len = length(z_probe);

matrix_size = [x_len, y_len, z_len];
distance_error_check = struct();

% ??? faccio un distance_error_check per ogni LED?
distance_error_check.LED1 = struct();

impulseMatrix  = zeros(matrix_size);
overlap_info   = zeros(matrix_size);
dataRateMatrix = zeros(matrix_size);

tic; % tempo di calcolo
% PDect = struct();

if plot_2
    for x_index = 1:x_len
        for y_index = 1:y_len
            for z_index = 1:z_len
                PDect = [x_probe(x_index), y_probe(y_index), z_probe(z_index)];

                % se attivato rallenta di molto il processo di calcolo
                if PLOT3D_enabled
                    plotCube(LED1, LED2, LED3, LED4, PDect, x_max, y_max, z_max);
                end

                % carica le matrici degli impulsi e delle sovrapposizioni con i dati
                [ impulseMatrix(x_index,y_index,z_index), overlap_info(x_index,y_index,z_index) ] = ...
                    drawChannelResponse( ...
                    Psi, LED1, LED2, LED3, LED4, PDect, Phi_FoV, ...
                    A_pd, T_of, a, Entity_enabled, alpha, beta ...
                    );

                lowerBoundDataRateMatrix(x_index,y_index,z_index) = ...
                    lowerBoundDataRate(impulseMatrix(x_index,y_index,z_index), B, p, R_pd, q, N);
            end
        end
    end
    toc;

    %% plot section

    pause(0.01);
    close all;
    if PLOT2D_enabled
        rel_min = 0;
        nan_idx = isnan(impulseMatrix);
        rel_max = mean(mean(mean(impulseMatrix(~nan_idx)))); % media di tutti i valori di impulso validi
        rel_max = rel_max * 2.5;

        z_index = 1:length(z_probe);
        plot_index = linspace(1, max(z_index)-2, 9);

        figure();
        colorbar;

        for i=1:9
            subplot(3, 3, i);
            matrix_slice = squeeze(impulseMatrix(:,:,floor(plot_index(i))));
            imagesc(matrix_slice);

            colorbar;
            clim([rel_min rel_max]);   % same scale
            xlabel('Y');
            ylabel('X');
            title('Z = '+string((z_probe(floor(plot_index(i)))))+' m');
        end
        colormap(jet(126));
        sgtitle('Z-sliced - Room size ' + string(x_max) +'x' + string(y_max)+'x'+string(z_max) +' meters');

        y_index = 1:length(y_probe);
        plot_index = linspace(1, max(y_index)-2, 9);

        figure();
        colorbar;

        for i=1:9
            subplot(3, 3, i);
            matrix_slice = squeeze(impulseMatrix(:,floor(plot_index(i)),:));
            imagesc(matrix_slice);
            camroll(90);
            colorbar;
            clim([rel_min rel_max]); % same scale
            xlabel('Z');
            ylabel('X');
            title('Y = '+string((y_probe(floor(plot_index(i)))))+' m el:'+plot_index(i));
        end
        colormap(jet(126));
        sgtitle('Y-sliced - Room size ' + string(x_max) +'x' + string(y_max)+'x'+string(z_max) +' meters');

        x_index = 1:length(x_probe);
        plot_index = linspace(1, max(x_index)-2, 9);

        figure();
        colorbar;

        for i=1:9
            subplot(3,3, i);
            matrix_slice = squeeze(impulseMatrix(floor(plot_index(i)),:,:));
            imagesc(matrix_slice);
            camroll(90)
            colorbar;
            clim([rel_min rel_max]);   % same scale
            xlabel('Z');
            ylabel('Y');
            title('X = '+string((x_probe(floor(plot_index(i)))))+' m el:'+plot_index(i));
        end
        colormap(jet(126));
        sgtitle('X-sliced - Room size ' + string(x_max) +'x' + string(y_max)+'x'+string(z_max) +' meters');

        % overlapping
        z_index = 1:length(z_probe);

        plot_index = linspace(1, max(z_index)-2, number_of_plot_needed);
        figure();
        colorbar;

        for i=1:number_of_plot_needed
            subplot(number_of_plot_needed/3, 3, i);
            matrix_slice = squeeze(overlap_info(:,:,floor(plot_index(i))));
            imagesc(matrix_slice);

            colorbar
            clim([0 5]);   % same scale
            xlabel('Y');
            ylabel('X');
            title('Z = '+string((z_probe(floor(plot_index(i)))))+' m');
        end

        sgtitle('Overlapping MAP Z-sliced | Room size: ' + ...
            string(x_max) + 'x' + string(y_max)+'x'+string(z_max) +' meters');

    end
end

%% single entity distance error estimation
if plot_1
    pause(0.01);
    z_plot_index = linspace(0, max(z_probe), number_of_plot_needed);
    distanceErrorEstimation(LED1, number_of_plot_needed, p, ...
        alpha, beta, Phi_FoV, a, Psi, A_pd, T_of, R_pd, q_0, k_B, T_k, eta, ...
        I_2, I_3, Gamma, g_m, I_bg, G_0, B, x_probe, y_probe, z_plot_index, number_of_samples, "1")
    distanceErrorEstimation(LED2, number_of_plot_needed, p, ...
        alpha, beta, Phi_FoV, a, Psi, A_pd, T_of, R_pd, q_0, k_B, T_k, eta, ...
        I_2, I_3, Gamma, g_m, I_bg, G_0, B, x_probe, y_probe, z_plot_index, number_of_samples, "2")
    distanceErrorEstimation(LED3, number_of_plot_needed, p, ...
        alpha, beta, Phi_FoV, a, Psi, A_pd, T_of, R_pd, q_0, k_B, T_k, eta, ...
        I_2, I_3, Gamma, g_m, I_bg, G_0, B, x_probe, y_probe, z_plot_index, number_of_samples, "3")
    distanceErrorEstimation(LED4, number_of_plot_needed, p, ...
        alpha, beta, Phi_FoV, a, Psi, A_pd, T_of, R_pd, q_0, k_B, T_k, eta, ...
        I_2, I_3, Gamma, g_m, I_bg, G_0, B, x_probe, y_probe, z_plot_index, number_of_samples, "4")
end

%% noise contribution
if plot_3
    shift_distance = 1; % meter
    granularity = 0.01;

    x_probe = 0:granularity:x_max;
    y_probe = 0:granularity:y_max;
    z_probe = 0:9*granularity:z_max;

    x_len = length(x_probe);
    y_len = length(y_probe);
    z_len = length(z_probe);

    B_array = [5e6, 20e6, 100e6, 400e6];
    marker_array = ["x", "o", "^", "p"];
    legend_array_led = ["\mu_0"];
    legend_array_led = [];
    legend_array_ris = ["\Delta_n"];
    legend_array_ris = [];

    for bandIdx = 1:length(B_array)
        B=B_array(bandIdx);

        pause(0.01);
        clc;
        disp("Starting to draw the noise contribution...");
        tic
        received_power_LED = zeros(1, z_len);
        n_shoot_LED_contribution = zeros(1, z_len);
        n_thermal_LED_contribution = zeros(1, z_len);
        distance_LED = zeros(1,z_len);

        received_power_NLoS = zeros(1, length(y_probe));
        n_shoot_RIS1_contribution = zeros(1, length(y_probe));
        n_thermal_RIS1_contribution = zeros(1, length(y_probe));
        distance_RIS1 = zeros(1, length(y_probe));

        max_rnd = number_of_samples;
        posLEDError = NaN(z_len, max_rnd);
        % par
        for z_index = 1:(z_len)

            % only LoS contribution
            PDect_pos_LED = [LED1(1), LED1(2), LED1(3) - z_probe(z_index)];
            received_power_LED(z_index) = R_pd*p*singleEntityContribution(LED1, 0, PDect_pos_LED, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 1); % Ampere
            distance_LED(z_index) = calculateDistance(LED1, PDect_pos_LED);

            for idx_rnd = 1:max_rnd
                [ n_shoot_LED, n_thermal_LED, nsh_var, nth_var] = noiseEstimation(received_power_LED(z_index), q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
                n_shoot_LED_contribution(z_index) = nsh_var;
                n_thermal_LED_contribution(z_index) = nth_var;
                received_power_LED_w_noise = received_power_LED(z_index) + n_shoot_LED + n_thermal_LED;
                dLED = getLEDDistanceByEstimatedPower(LED1(3) - PDect_pos_LED(3), p, A_pd, received_power_LED_w_noise, Psi, T_of, Phi_FoV, a, R_pd);

                if isnan(dLED)
                    posLEDError(z_index, idx_rnd) = NaN;
                else
                    posLEDError(z_index, idx_rnd) = abs(distance_LED(z_index) - dLED);
                end
            end

        end

        %posRISError = NaN(z_len, max_rnd);
        posRISError = NaN(z_len, max_rnd);

        for z_index = 1:z_len-2
            if LED1(3)-z_probe(z_index) > 0
                PDect_pos_RIS1 = [LED1(1), LED1(2)+y_probe(z_index), LED1(3)-z_probe(z_index)];
                % LoS Contribution
                received_power_LoS = R_pd*p*singleEntityContribution(LED1, 0, [LED1(1) LED1(2) LED1(3)-shift_distance], alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 1); % Ampere
                %distance_RIS1_LoS = calculateDistance(LED1, PDect_pos_RIS1);

                % NLoS Contribution
                received_power_NLoS(z_index) = R_pd*p*singleEntityContribution(LED1, LED1, PDect_pos_RIS1, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 0);
                distance_RIS1(z_index)=calculateDistance(LED1, PDect_pos_RIS1);

                for idx_rnd = 1:max_rnd
                    % noise due to LoS
                    [n_shoot_LoS, n_thermal_LoS, nsh_var_LoS, nth_var_LoS] = ...
                        noiseEstimation(received_power_LoS, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
                    %n_shoot_LoS_contribution(z_index) = nsh_var_LoS;
                    %n_thermal_LoS_contribution(z_index) = nth_var_LoS;

                    % noise due to LoS+NLoS
                    [n_shoot_RIS1, n_thermal_RIS1, nsh_var, nth_var] = ...
                        noiseEstimation(received_power_LoS + received_power_NLoS(z_index), q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
                    n_shoot_RIS1_contribution(z_index) = nsh_var_LoS + nsh_var;
                    n_thermal_RIS1_contribution(z_index) = nth_var_LoS + nth_var;

                    received_power_RIS1_w_noise  =  received_power_NLoS(z_index) + n_shoot_RIS1 + n_thermal_RIS1 + n_shoot_LoS + n_thermal_LoS;
                    dRIS1 = getRISDistanceByEstimatedPower(LED1(3) - PDect_pos_RIS1(3), p, A_pd, received_power_RIS1_w_noise, Psi, T_of, Phi_FoV, a, rho, cosd(calculateAngle(LED1, LED1, 0, 0)), R_pd, calculateDistance(LED1, LED1));

                    if isnan(dRIS1)
                        posRISError(z_index, idx_rnd) = NaN;
                    else
                        posRISError(z_index, idx_rnd) = abs(distance_RIS1(z_index) - dRIS1);
                    end
                end

            else
                break;
            end
        end


        %legend_array_led = [legend_array_led , "\sigma_s with B = "+string(B*1e-6)+" [MHz]", "\sigma_t with B = "+string(B*1e-6)+" [MHz]"];


        figure(100);
        %grid on;
        hold on;
        %title("Noise estimation LED1");
        ylabel('Signal Power [Ampere^2]');
        xlabel('Distance [m]');
        if bandIdx == 1
            %marker_array = ["x", "o", "^", "p"];
            plot(NaN, 'b-', 'LineWidth',2, 'MarkerSize', 15,'MarkerFaceColor','k');
            plot(NaN, 'r-.', 'LineWidth',2, 'MarkerSize', 15,'MarkerFaceColor','k');
            plot(NaN, 'k-', 'LineWidth',2, 'MarkerSize', 15,'MarkerFaceColor','k');
            plot(NaN, 'k--', 'LineWidth',2, 'MarkerSize', 15,'MarkerFaceColor','k');
            plot(NaN, 'kx', 'LineWidth',2, 'MarkerSize', 15,'MarkerFaceColor','k');
            plot(NaN, 'ko', 'LineWidth',2, 'MarkerSize', 15,'MarkerFaceColor','k');
            plot(NaN, 'k^', 'LineWidth',2, 'MarkerSize', 15,'MarkerFaceColor','k');
            plot(NaN, 'kp', 'LineWidth',2, 'MarkerSize', 15,'MarkerFaceColor','k');

            plot(distance_LED,received_power_LED.^2,"b-", 'LineWidth',2);
            legend_array_led = [legend_array_led , "LED", "IMR", "\mu_0^2 or \Delta_n^2", "\sigma^2", "B = 5 MHz", "B = 20 MHz", "B = 100 MHz", "B = 400 MHz"];
            legend(legend_array_led);
        end
        %plot(distance_LED,sqrt(n_shoot_LED_contribution), 'LineWidth',2);
        %plot(distance_LED,sqrt(n_thermal_LED_contribution), '-.', 'LineWidth',2);
        plot(distance_LED,(n_shoot_LED_contribution + n_thermal_LED_contribution), 'b--',...
            'LineWidth',2, 'Marker', marker_array(bandIdx), 'MarkerSize', 15,'MarkerFaceColor','b', 'MarkerIndices', [2 1:12:145]);

        %ylim([10e-25 10e-7]);
        %legend(legend_array_led);
        set(gca, 'YScale', 'log');
        set(gca, "FontName", "Times New Roman", "FontSize", 33);

        %legend_array_ris = [legend_array_ris , "\sigma_s with B = "+string(B*1e-6)+" [MHz]", "\sigma_t B = "+string(B*1e-6)+" [MHz]"];
        figure(100);
        %grid on;
        hold on;
        %title("Noise estimation RIS1");
        ylabel('Signal Power [Ampere^2]');
        xlabel('Distance [m]');
        if bandIdx == 1
            plot(distance_RIS1,received_power_NLoS.^2,"r-", 'LineWidth',2);
            %legend_array_ris = [legend_array_ris , "\sigma", "B = 5 MHz", "B = 20 MHz", "B = 100 MHz", "B = 400 MHz"];
            %legend(legend_array_ris);
        end
        %plot(distance_RIS1,sqrt(n_shoot_RIS1_contribution), 'LineWidth',2);
        %plot(distance_RIS1,sqrt(n_thermal_RIS1_contribution), '-.', 'LineWidth',2);
        plot(distance_RIS1,n_shoot_RIS1_contribution + n_thermal_RIS1_contribution, 'r--',...
            'LineWidth', 2, 'Marker', marker_array(bandIdx), 'MarkerSize', 15, 'MarkerFaceColor','r', 'MarkerIndices', [2 6:12:1301]);

        hold off;
        %ylim([10e-25 10e-7]);
        %legend(["SIGNAL", "Shot", "Thermal"]);
        %legend(legend_array_ris);
        set(gca, 'YScale', 'log');
        xlim([0 5]);
        ylim([1e-18 1e-6]);
        set(gca, "FontName", "Times New Roman", "FontSize", 33);
    end

    %% grafico SNR al variare della potenza - Fig 1 A e B
    % LED, RIS1
    figure(101);
    hold on;
    plot(NaN, 'b-', 'LineWidth',2, 'MarkerSize', 15,'MarkerFaceColor','k');
    plot(NaN, 'r-.', 'LineWidth',2, 'MarkerSize', 15,'MarkerFaceColor','k');
    plot(NaN, 'kx', 'LineWidth',2, 'MarkerSize', 15,'MarkerFaceColor','k');
    plot(NaN, 'ko', 'LineWidth',2, 'MarkerSize', 15,'MarkerFaceColor','k');
    plot(NaN, 'k^', 'LineWidth',2, 'MarkerSize', 15,'MarkerFaceColor','k');
    plot(NaN, 'kp', 'LineWidth',2, 'MarkerSize', 15,'MarkerFaceColor','k');
    legend(["LED", "IMR", "B = 5MHz", "B = 20MHz", "B = 100MHz", "B = 400MHz"]);

    for bandIdx = 1:length(B_array)

        B=B_array(bandIdx);

        received_power_LED = [];
        distance_LED = [];
        received_power_NLoS = [];
        distance_RIS1 = [];
        lumen_array = [1:10:10000];




        for lumenIdx = 1: length(lumen_array)
            lumen_level = lumen_array(lumenIdx);
            p = lumen_level / 683;      % transmission power 6000 [Lumens] -> [Watt]


            PDect_pos_LED = [LED(1), LED(2), LED(3) - shift_distance];
            received_power_LED(lumenIdx) = R_pd*p*singleEntityContribution(LED, 0, PDect_pos_LED, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 1); % Ampere
            distance_LED(lumenIdx) = calculateDistance(LED, PDect_pos_LED);
            [n_shoot_LED, n_thermal_LED, nsh_var_LED, nth_var_LED] = noiseEstimation(received_power_LED(lumenIdx), q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);

            SNR_LED(lumenIdx) = 10* log10(received_power_LED(lumenIdx).^2 / (nsh_var_LED+nth_var_LED));


            PDect_pos_RIS1 = [LED1(1),LED1(2), LED1(3) - shift_distance];
            received_power_LoS = R_pd*p*singleEntityContribution(LED, 0, PDect_pos_LED, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 1); % Ampere
            %distance_RIS1_LoS = calculateDistance(LED, PDect_pos_RIS1);
            % NLoS Contribution
            received_power_NLoS(lumenIdx) = R_pd*p*singleEntityContribution(LED, LED1, PDect_pos_RIS1, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 0);
            distance_RIS1(lumenIdx)=calculateDistance(LED1, PDect_pos_RIS1);
            [n_shoot_LoS, n_thermal_LoS, nsh_var_LoS, nth_var_LoS] = ...
                noiseEstimation(received_power_LoS, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
            %n_shoot_LoS_contribution(z_index) = nsh_var_LoS;
            %n_thermal_LoS_contribution(z_index) = nth_var_LoS;

            % noise due to LoS+NLoS
            [n_shoot_RIS1, n_thermal_RIS1, nsh_var, nth_var] = ...
                noiseEstimation(received_power_LoS + received_power_NLoS(lumenIdx), q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
            n_shoot_RIS1_contribution(lumenIdx) = nsh_var_LoS + nsh_var;
            n_thermal_RIS1_contribution(lumenIdx) = nth_var_LoS + nth_var;


            SNR_RIS(lumenIdx) = 10* log10(received_power_NLoS(lumenIdx).^2 / (n_shoot_RIS1_contribution(lumenIdx)  + n_thermal_RIS1_contribution(lumenIdx)));

        end

        plot(lumen_array, SNR_LED, 'b--', 'Marker', marker_array(bandIdx), 'MarkerIndices', 1:100:1000, 'MarkerSize', 15, 'MarkerFaceColor','b');

        plot(lumen_array, SNR_RIS, 'r-.', 'Marker', marker_array(bandIdx), 'MarkerIndices', 1:100:1000, 'MarkerSize', 15, 'MarkerFaceColor','r');
        ylabel('SNR [dB]');
        xlabel('Optical Power [Lumen]');
    end


    set(gca, "FontName", "Times New Roman", "FontSize", 39);

    %% boxplot per errore medio
    %     figure(102);
    %     grid on;
    %     hold on;
    %     title("BoxPlot Error LED");
    %     xlabel('Distanza [m]');
    %     ylabel('Errore distanza [m]');
    %
    %     boxchart(posLEDError', 'MarkerStyle','none');
    %     %ylim([0 max(max(posLEDError))]);
    %     ylim([0 y_max/6]);
    %     ax = gca;
    %     %ax.YAxis.Scale ="log";
    %
    %     figure();
    %     grid on;
    %     hold on;
    %     title("BoxPlot Error RIS1");
    %     xlabel('Distanza [m]');
    %     ylabel('Errore distanza [m]');
    %     boxchart(posRISError', 'MarkerStyle','none');
    %     ylim([0 y_max/2]);
    %     ax = gca;
    toc;

end

%% Fig 2 epsilon n al variare di qualcosa
if plot_5
    m = -(log(2)/log(cosd(Psi)));
    G = (a^2)/((sind(Phi_FoV)^2));

    % fissiamo
    K_0 = 5;
    K_n = 10;
    epsilon = 0.001;

    epsilon_0_representation = nan(x_len, y_len);
    epsilon_n_representation= nan(x_len, y_len);

    if plot_4
        for x_index = 2:x_len
            for y_index = 2:y_len
                %for z_index = 1:z_len
                PDect = [x_probe(x_index), y_probe(y_index),0];

                p_LoS = R_pd*p*singleEntityContribution(LED, 0, PDect, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 1); % Ampere
                p_NLoS = R_pd*p*singleEntityContribution(LED, LED1, PDect, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 0);

                theta_RS = calculateAngle(LED, LED1, alpha, beta);
                r_n = calculateDistance(LED, LED1);
                d_n = calculateDistance(LED1, PDect);

                %TODO

                [~, ~, nsh_var_LoS, nth_var_LoS] = ...
                    noiseEstimation(p_LoS, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);

                [~, ~, nsh_var_NLoS, nth_var_NLoS] = ...
                    noiseEstimation(p_LoS + p_NLoS, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);

                sigma_0 = sqrt(nsh_var_LoS + nth_var_LoS);
                sigma_n = sqrt(nsh_var_NLoS + nth_var_NLoS);

                beta_n  = ((rho * G * R_pd * A_pd * p*(m+1))/(2*pi)) * cosd(theta_RS)^m * LED1(3);

                delta_0 = sqrt((2/K_0) * sigma_0^2 * erfinv(1-epsilon)^2);

                delta_0_tresh = sqrt(2*sigma_0^2 * erfinv(1-epsilon));

                if delta_0 > delta_0_tresh
                    error('ERRORE: delta_0 > delta_0_tresh')
                end

                delta_n = sqrt((2/K_n)*(sigma_n * erfinv(1-epsilon))^2 + delta_0^2);

                %             delta_n_tresh;
                %
                %             if delta_n > delta_n_tresh
                %                 error('delta_n > delta_n_tresh')
                %             end


                if beta_n < (delta_n*(r_n + d_n)^2 * d_n)
                    error('ERRORE: beta_n < (delta_n*(r_n + d_n)^2 * d_n)')
                end

                gamma_n = ((r_n+d_n)^2 * d_n * beta_n)/(beta_n - (delta_n*(r_n + d_n)^2 * d_n));
                alpha_n = gamma_n;
                %
                %             a = r_n;
                %             A = gamma_n;
                %             b = d_n;
                %             epsilon_n = (1/3) * ((3 * sqrt(3) * sqrt(4 * a^3 * A + 27 * A^2) + 2 * a^3 + 27 * A)^(1/3) / (2^(1/3)) + (2^(1/3) * a^2) / ((3 * sqrt(3) * sqrt(4 * a^3 * A + 27 * A^2) + 2 * a^3 + 27 * A)^(1/3)) - 2 * a - 3 * b);

                rad = 1 / (sqrt((27 * r_n^3 * alpha_n) + (729/4)*alpha_n^2) + r_n^3 + (27/2)*alpha_n);
                tao_n = nthroot(rad,3);
                d_est_n = (tao_n/3) * (r_n - (1/tao_n))^2;
                epsilon_n = d_est_n - d_n;
                epsilon_n_representation(x_index, y_index) = epsilon_n;

                d_0 = calculateDistance(LED, PDect);
                rad_0 = (p_LoS)/(p_LoS-delta_0);
                sol_rad = nthroot(rad_0,m+3);
                epsilon_0 = d_0 * (sol_rad - 1);

                epsilon_0_representation(x_index, y_index) = epsilon_0;




                %end
            end
        end


    end
    figure(89);
    h = surf(x_probe, y_probe, epsilon_0_representation', 'EdgeColor', 'None');
    %rotate(h, [0 1 0], 90);
    ylabel('X');
    xlabel('Y');
    hold on;
    plot(LED(2),LED(1), 'o', 'MarkerSize', 25 );


    figure(90);
    surf(y_probe, x_probe, epsilon_n_representation);
    ylabel('X');
    xlabel('Y');
    hold on;
    plot(LED1(2),LED1(1), 'x', 'MarkerSize', 25 );
    plot(LED(2),LED(1), 'o', 'MarkerSize', 25 );

    figure(91);
    [C, h] = contour(x_probe, y_probe, epsilon_n_representation',[0.01 0.02 0.03 0.04 0.05 0.07 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.5],"ShowText",true, 'LineWidth', 3);
    clabel(C,h, 'FontSize', 30);
    ylabel('y [m]');
    xlabel('x [m]');
    hold on;
    % plot(RIS1(2),RIS1(1), 'x', 'MarkerSize', 25 );
    % plot(LED(2),LED(1), 'o', 'MarkerSize', 25 );


    set(gca, "FontName", "Times New Roman", "FontSize", 39);


    %%
    K0_array = 1:1000;
    Kn_array = 1:1000;

    epsilon_n_representation= nan(x_len, y_len);
    PDect = [x_max/2, y_max/2, 0];

    p_LoS = R_pd*p*singleEntityContribution(LED, 0, PDect, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 1); % Ampere

    p_NLoS = R_pd*p*singleEntityContribution(LED, LED1, PDect, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 0);

    theta_RS = calculateAngle(LED, LED1, alpha, beta);

    r_n = calculateDistance(LED, LED1);
    d_n = calculateDistance(LED1, PDect);

    %TODO

    [n_shoot_LoS, n_thermal_LoS, nsh_var_LoS, nth_var_LoS] = ...
        noiseEstimation(p_LoS, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);

    [n_shoot_RIS1, n_thermal_RIS1, nsh_var_NLoS, nth_var_NLoS] = ...
        noiseEstimation(p_LoS + p_NLoS, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);

    sigma_0 = sqrt(nsh_var_LoS + nth_var_LoS);
    sigma_n = sqrt(nsh_var_NLoS + nth_var_NLoS);

    beta_n  = ((rho * G * R_pd * A_pd * p*(m+1))/(2*pi)) * cosd(theta_RS)^m * LED1(3);

    if plot_4
        for x_index = 1:length(K0_array)
            for y_index = 1:length(Kn_array)
                %for z_index = 1:z_len

                K_0 = K0_array(x_index);
                K_n = Kn_array(y_index);


                delta_0 = sqrt((2/K_0) * sigma_0^2 * erfinv(1-epsilon)^2);

                delta_0_tresh = sqrt(2*sigma_0^2 * erfinv(1-epsilon));

                delta_n = sqrt((2/K_n)*(sigma_n * erfinv(1-epsilon))^2 + delta_0^2);

                gamma_n = ((r_n+d_n)^2 * d_n * beta_n)/(beta_n - (delta_n*(r_n + d_n)^2 * d_n));
                alpha_n = gamma_n;
                %
                %             a = r_n;
                %             A = gamma_n;
                %             b = d_n;
                %             epsilon_n = (1/3) * ((3 * sqrt(3) * sqrt(4 * a^3 * A + 27 * A^2) + 2 * a^3 + 27 * A)^(1/3) / (2^(1/3)) + (2^(1/3) * a^2) / ((3 * sqrt(3) * sqrt(4 * a^3 * A + 27 * A^2) + 2 * a^3 + 27 * A)^(1/3)) - 2 * a - 3 * b);

                rad = 1 / (sqrt((27 * r_n^3 * alpha_n) + (729/4)*alpha_n^2) + r_n^3 + (27/2)*alpha_n);
                tao_n = nthroot(rad,3);
                d_est_n = (tao_n/3) * (r_n - (1/tao_n))^2;
                epsilon_n = d_est_n - d_n;

                if epsilon_n > 0
                    epsilon_n_representation(x_index, y_index) = epsilon_n;
                else
                    epsilon_n_representation(x_index, y_index) = nan;
                end


                %end
            end
        end


    end
    figure(92);
    surf(Kn_array, K0_array, epsilon_n_representation);
    %surf(epsilon_n_representation, 'EdgeColor', 'None');
    figure(93);
    [C,h]  = contour(epsilon_n_representation, [0.01 0.15 0.02 0.25 0.030 0.04 0.05 0.07 0.1 0.12 0.15 0.2 0.25] ,"ShowText",true, 'LineWidth', 3);
    clabel(C,h, 'FontSize', 30);
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    %caxis([0 0.5]);
    ylabel('K_0');
    xlabel('K_n');
    hold on;

    set(gca, "FontName", "Times New Roman", "FontSize", 39);
end
%% estimate position

pause(0.01);
fprintf("Starting to minimizing position error...\n");
tic
%     lumen_level = 1000;
%     p = lumen_level / 683;      % transmission power 6000 [Lumens] -> [Watt]
%     B       = 100e6;             % System bandwidth [Hz]

number_of_plot_needed = 1;

if plot_6
    figure(300);
    x_probe = 0:granularity:x_max;
    y_probe = 0:granularity:y_max;
    z_probe = 0:9*granularity:z_max;

    %number_of_plot_needed = 18;
    z_plot_index = linspace(0, max(z_probe), number_of_plot_needed);
    z_plot_index = 0;
    colorbar;

    %for k = 1:length(z_probe)
    %fprintf("Sample n: %d", number_of_samples);


    for k=1:number_of_plot_needed
        error_matrix = zeros(length(x_probe), length(y_probe));

        for i=1:length(x_probe)
            %par
            for j=1:length(y_probe)
                PD_to_findx = [x_probe(i) y_probe(j) z_plot_index(k)];
                [p1, p2, p3] = estimateReceiverPosition(LED, LED1, LED2, LED3, LED4, PD_to_findx,...
                    p, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, R_pd, q_0, k_B, T_k, eta, I_2, ...
                    I_3, Gamma, g_m, I_bg, G_0, B, K_0, K_n);
                if isnan(p1) || isnan(p2) || isnan(p3)
                    error_matrix(i,j) = NaN;
                else
                    error_matrix(i,j) = real(sqrt((p1 - PD_to_findx(1))^2 + (p2 - PD_to_findx(2))^2 + (p3 - PD_to_findx(3))^2)); % RMSE
                end
            end
        end
        %subplot(3,3,k);
        %subplot(number_of_plot_needed/3, 3, k);
        %figure(200);
        hold on;
        surf(error_matrix);
        xlim([0 length(x_probe)]);
        ylim([0 length(y_probe)]);
        colormap;
        xlabel('Y');
        ylabel('X');

        error_matrix_vec = error_matrix(:);
        nanIdx = isnan(error_matrix_vec);
        infIdx = isinf(error_matrix_vec);
        error_matrix_vec = error_matrix_vec(~nanIdx & ~infIdx);


        if isempty(error_matrix_vec(~isnan(error_matrix_vec) | error_matrix_vec ~= Inf))
            mean_error_string = "...";
        else
            mean_error_string = string(mean(mean(error_matrix_vec)));
        end

        title('Z = '+string((z_plot_index(k)))+' [m] - mean ' + mean_error_string);

    end
    sgtitle('Z-sliced - Room size ' + string(x_max) +'x' + string(y_max)+'x'+string(z_max) +' [m]');

end



%% errore di posizionamento con 3 entità
if plot_7
    dynamic_epsilon_estimation = 1;
    epsilon = 0.001;
    K_0 = 50;
    K_n = 2*K_0;


    % scelta epsilon

    %     epsilon_0 = 0.3;
    %     epsilon_1 = 0.3;
    %     epsilon_2 = 0.3;
    %     epsilon_3 = 0.3;
    %     epsilon_4 = 0.3;
    %

    figure(998);
    x1 = LED(1);
    y1 = LED(2);
    x2 = LED1(1);
    y2 = LED1(2);
    x3 = LED2(1);
    y3 = LED2(2);

    A =[...
        x2-x1, y2-y1; ...
        x3-x1, y3-y1; ...
        ];

    for i=1:length(x_probe)
        %par
        for j=1:length(y_probe)
            PD_to_findx = [x_probe(i) y_probe(j) 0];

            d_0 = calculateDistance(LED, PD_to_findx);
            d_1 = calculateDistance(LED1, PD_to_findx);
            d_2 = calculateDistance(LED2, PD_to_findx);

            if dynamic_epsilon_estimation
                epsilon_0 = estimate_epsilon_0 (K_0, epsilon, LED, PD_to_findx, R_pd, p, q_0, k_B, T_k, eta, I_2, I_3, Gamma,g_m, I_bg, G_0, B, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of);
                epsilon_1 = estimate_epsilon_n(K_0, K_n, epsilon,LED, PD_to_findx, LED1, p, Phi_FoV, a,rho, Psi, T_of, alpha, beta,  q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
                epsilon_2 = estimate_epsilon_n(K_0, K_n, epsilon,LED, PD_to_findx, LED2, p, Phi_FoV, a, rho,Psi, T_of, alpha, beta, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
                %epsilon_3 = estimate_epsilon_n(K_0, K_n, epsilon,LED, PD_to_findx, RIS3, p, Phi_FoV, a, rho,Psi, T_of, alpha, beta, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
                %epsilon_4 = estimate_epsilon_n(K_0, K_n, epsilon,LED, PD_to_findx, RIS4, p, Phi_FoV, a, rho, Psi, T_of, alpha, beta, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
            end

            y_tilde = 0.5 * [...
                epsilon_0*(epsilon_0 + 2*d_0) - epsilon_1*(epsilon_1+2*d_1);...
                epsilon_0*(epsilon_0 + 2*d_0) - epsilon_2*(epsilon_2+2*d_2);...
                ];
            error_vector = pinv(A)* y_tilde;
            error_matrix(i,j) = sqrt(error_vector(1)^2 + error_vector(2)^2);

        end
    end

    % error_matrix(length(x_probe)-5, length(y_probe)-6) = 0.0;

    hold on;

    min_error = min(min(error_matrix));
    max_error = max(max(error_matrix));

    plot(NaN, 'ro', 'LineWidth',3, 'MarkerSize', 25);
    plot(NaN, 'rx', 'LineWidth',3, 'MarkerSize', 25);

    plot3(LED(2),LED(1),max_error+0.1, 'ro', 'LineWidth', 3,'MarkerSize', 25 );
    plot3(LED1(1),LED1(2),max_error+0.1, 'rx', 'LineWidth', 3, 'MarkerSize', 25 );
    plot3(LED2(1),LED2(2),max_error+0.1, 'rx','LineWidth', 3, 'MarkerSize', 25 );
    surf (x_probe, y_probe, error_matrix', 'EdgeColor', 'None' );


    cb = colorbar;
    cb.Label.String = "RMSE";
    cb.Label.Rotation = 270;
    cb.Label.VerticalAlignment = "bottom";

    xlabel('x [m]');
    ylabel('y [m]');
    caxis([0 0.25]);
    xlim([0 x_max]);
    ylim([0 y_max]);
    set(gca, "FontName", "Times New Roman", "FontSize", 39);

    %% errore di posizionamento con 5 entità

    figure(999);

    x1 = LED(1);    y1 = LED(2);
    x2 = LED1(1);   y2 = LED1(2);
    x3 = LED2(1);   y3 = LED2(2);
    x4 = LED3(1);   y4 = LED3(2);
    x5 = LED4(1);   y5 = LED4(2);


    A =[...
        x2-x1, y2-y1; ...
        x3-x1, y3-y1; ...
        x4-x1, y4-y1; ...
        x5-x1, y5-y1; ...
        ];

    for i=1:length(x_probe)
        %par
        for j=1:length(y_probe)
            PD_to_findx = [x_probe(i) y_probe(j) 0];

            d_0 = calculateDistance(LED, PD_to_findx);
            d_1 = calculateDistance(LED1, PD_to_findx);
            d_2 = calculateDistance(LED2, PD_to_findx);
            d_3 = calculateDistance(LED3, PD_to_findx);
            d_4 = calculateDistance(LED4, PD_to_findx);


            if dynamic_epsilon_estimation
                epsilon_0 = estimate_epsilon_0 (K_0, epsilon, LED, PD_to_findx, R_pd, p, q_0, k_B, T_k, eta, I_2, I_3, Gamma,g_m, I_bg, G_0, B, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of);
                epsilon_1 = estimate_epsilon_n(K_0, K_n, epsilon,LED, PD_to_findx, LED1, p, Phi_FoV, a,rho, Psi, T_of, alpha, beta, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
                epsilon_2 = estimate_epsilon_n(K_0, K_n, epsilon,LED, PD_to_findx, LED2, p, Phi_FoV, a, rho,Psi, T_of, alpha, beta, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
                epsilon_3 = estimate_epsilon_n(K_0, K_n, epsilon,LED, PD_to_findx, LED3, p, Phi_FoV, a, rho,Psi, T_of, alpha, beta, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
                epsilon_4 = estimate_epsilon_n(K_0, K_n, epsilon,LED, PD_to_findx, LED4, p, Phi_FoV, a, rho,Psi, T_of, alpha, beta, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);
            end

            y_tilde = 0.5 * [...
                epsilon_0*(epsilon_0 + 2*d_0) - epsilon_1*(epsilon_1+2*d_1);...
                epsilon_0*(epsilon_0 + 2*d_0) - epsilon_2*(epsilon_2+2*d_2);...
                epsilon_0*(epsilon_0 + 2*d_0) - epsilon_3*(epsilon_3+2*d_3);...
                epsilon_0*(epsilon_0 + 2*d_0) - epsilon_4*(epsilon_4+2*d_4);...
                ];
            error_vector = pinv(A)* y_tilde;
            error_matrix(i,j) = sqrt(error_vector(1)^2 + error_vector(2)^2);
        end
    end
    %[C,h] = contour(y_probe, x_probe, error_matrix ,"ShowText",true, 'LineWidth', 3);
    %clabel(C,h, 'FontSize', 30);
    hold on;
    plot(NaN, 'ro', 'LineWidth',3, 'MarkerSize', 25);
    plot(NaN, 'rx', 'LineWidth',3, 'MarkerSize', 25);
    surf(x_probe, y_probe, error_matrix', 'EdgeColor', 'None');

    xlabel('x [m]');
    ylabel('y [m]');

    plot3(LED(2),LED(1), max_error+0.1, 'ro', 'LineWidth', 3,'MarkerSize', 25 );
    plot3(LED1(2),LED1(1),max_error+0.1, 'rx', 'LineWidth', 3, 'MarkerSize', 25 );
    plot3(LED2(2),LED2(1),max_error+0.1, 'rx','LineWidth', 3, 'MarkerSize', 25 );
    plot3(LED3(2),LED3(1),max_error+0.1, 'rx','LineWidth', 3, 'MarkerSize', 25 );
    plot3(LED4(2),LED4(1),max_error+0.1, 'rx','LineWidth', 3, 'MarkerSize', 25 );

    min_error = min([min_error min(error_matrix)]);
    max_error = max([max_error max(error_matrix)]);

    %%%%%%%%%%%%%%%%%
    % min_error = max([min_error min(error_matrix)]);
    % max_error = min([max_error max(error_matrix)]);
    %%%%%%%%%%%%%%%%%

    caxis([min_error max_error-0.2]);
    legend('LED', 'IMR');

    cb = colorbar;
    cb.Label.String = "RMSE";
    cb.Label.Rotation = 270;
    cb.Label.VerticalAlignment = "bottom";
    set(gca, "FontName", "Times New Roman", "FontSize", 39);

    figure(998);
    xlim([0 x_max]);
    ylim([0 y_max]);
    caxis([min_error max_error-0.2]);
    legend('LED', 'IMR');
end


function epsilon_n = estimate_epsilon_n(K_0, K_n, epsilon, LED, PDect, RIS, p, Phi_FoV, a, rho, Psi, T_of, alpha, beta, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B)
m = -(log(2)/log(cosd(Psi)));

p_LoS = R_pd*p*singleEntityContribution(LED, 0, PDect, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 1); % Ampere

p_NLoS = R_pd*p*singleEntityContribution(LED, RIS, PDect, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 0);

theta_RS = calculateAngle(LED, RIS, alpha, beta);

r_n = calculateDistance(LED, RIS);
d_n = calculateDistance(RIS, PDect);

theta_SD = calculateAngle(RIS, PDect, alpha, beta);
if theta_SD >= 0 && theta_SD <= Phi_FoV
    G = (a^2)/((sind(Phi_FoV)^2));
else
    G = 0;
end


[n_shoot_LoS, n_thermal_LoS, nsh_var_LoS, nth_var_LoS] = ...
    noiseEstimation(p_LoS, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);

[n_shoot_RIS1, n_thermal_RIS1, nsh_var_NLoS, nth_var_NLoS] = ...
    noiseEstimation(p_LoS + p_NLoS, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);

sigma_0 = sqrt(nsh_var_LoS + nth_var_LoS);
sigma_n = sqrt(nsh_var_NLoS + nth_var_NLoS);

beta_n  = ((rho * G * R_pd * A_pd * p*(m+1))/(2*pi)) * cosd(theta_RS)^m * RIS(3);

delta_0 = sqrt((2/K_0) * sigma_0^2 * erfinv(1-epsilon)^2);

delta_0_tresh = sqrt(2*sigma_0^2 * erfinv(1-epsilon));

delta_n = sqrt((2/K_n)*(sigma_n * erfinv(1-epsilon))^2 + delta_0^2);

gamma_n = ((r_n+d_n)^2 * d_n * beta_n)/(beta_n - (delta_n*(r_n + d_n)^2 * d_n));
alpha_n = gamma_n;
%
%             a = r_n;
%             A = gamma_n;
%             b = d_n;
%             epsilon_n = (1/3) * ((3 * sqrt(3) * sqrt(4 * a^3 * A + 27 * A^2) + 2 * a^3 + 27 * A)^(1/3) / (2^(1/3)) + (2^(1/3) * a^2) / ((3 * sqrt(3) * sqrt(4 * a^3 * A + 27 * A^2) + 2 * a^3 + 27 * A)^(1/3)) - 2 * a - 3 * b);

rad = 1 / (sqrt((27 * r_n^3 * alpha_n) + (729/4)*alpha_n^2) + r_n^3 + (27/2)*alpha_n);
tao_n = nthroot(rad,3);
d_est_n = (tao_n/3) * (r_n - (1/tao_n))^2;

epsilon_n = d_est_n - d_n;

if epsilon_n > 0
    epsilon_n = epsilon_n;
else
    epsilon_n = nan;
end
end


function epsilon_0 = estimate_epsilon_0 (K_0, epsilon, LED, PDect, R_pd, p, q_0, k_B, T_k, eta, I_2, I_3, Gamma,g_m, I_bg, G_0, B, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of)
m = -(log(2)/log(cosd(Psi)));
p_LoS = R_pd*p*singleEntityContribution(LED, 0, PDect, alpha, beta, Phi_FoV, a, rho, Psi, A_pd, T_of, 1); % Ampere

[~, ~, nsh_var_LoS, nth_var_LoS] = ...
    noiseEstimation(p_LoS, q_0, R_pd, k_B, T_k, eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);

sigma_0 = sqrt(nsh_var_LoS + nth_var_LoS);

delta_0 = sqrt((2/K_0) * sigma_0^2 * erfinv(1-epsilon)^2);

d_0 = calculateDistance(LED, PDect);
rad_0 = (p_LoS)/(p_LoS-delta_0);
sol_rad = nthroot(rad_0,m+3);
epsilon_0 = d_0 * (sol_rad - 1);
end
