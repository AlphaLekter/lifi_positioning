% VLP Simulator
clc; clear; close all;

%% TODO
% aggiustare yticks grafici boxplot

%% LED/PD parameters
Psi     = 70;               % LED half-power semiangle [degree]
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
K_0 = 5;                    % Numero di campioni
epsilon = 0.001;

%% Noise Parameters
q_0     = 1.602e-19;        % electronic charge [Coulombs]
I_bg    = 84e-6;            % background light current [A] -> 84 [µA]
% I_bg    = 5e-12;           % 5pA
k_B     = 1.38064852e-23;   % Boltzmann constant [Joule/Kelvin]
T_k     = 295;              % absolute temperature [K]
G_0     = 10;               % open-loop voltage gain
eta     = 1.12e-6;          % fixed capacitance of photo detector per unit area [F/m^2]
Gamma   = 1.5;              % FET channel noise factor
g_m     = 0.030;            % FET transconductance [Siemens] [mS]
I_2     = 0.562;            % noise BW factor
I_3     = 0.0868;           % noise BW factor

%% Room Sizes
x_max       = 3;            % room size x-axis         % [SCA+2022]
y_max       = 3;            % room size y-axis         % [SCA+2022]
z_max       = 5;            % room size z-axis         % 3 in [SCA+2022]

%% Plot Settings

granularity = min([x_max, y_max, z_max]) / 20; % plot accuracy
number_of_samples = 30;
Entity_enabled = [1 1 1 1]; % LED1 - LED2 - LED3 - LED4

PLOT3D_enabled = 0;
PLOT2D_enabled = 1;
plot_1         = 1; % Light Power
plot_2         = 0; % Distance Error
plot_3         = 1; % Noise Contribution
plot_4         = 1; % Position Estimation
number_of_plot_needed = 9;

% tilt info
alpha   = 0;
beta    = 0;

%% LED positions
LED1 = [x_max/4   , y_max/4   , z_max];
LED2 = [x_max/4   , y_max*3/4 , z_max];
LED3 = [x_max*3/4 , y_max/4   , z_max];
LED4 = [x_max*3/4 , y_max*3/4 , z_max];

%% Estimate light power
x_probe = 0:granularity:x_max;
y_probe = 0:granularity:y_max;
z_probe = 0:granularity:z_max;

x_len = length(x_probe);
y_len = length(y_probe);
z_len = length(z_probe);

matrix_size = [x_len, y_len, z_len];

impulseMatrix  = zeros(matrix_size);
overlap_info   = zeros(matrix_size);
dataRateMatrix = zeros(matrix_size);
lowerBoundDataRateMatrix = zeros(x_len, y_len, z_len);

tic; % tempo di calcolo

if plot_1
    for x_index = 1:x_len
        for y_index = 1:y_len
            for z_index = 1:z_len
                PDect = [x_probe(x_index), y_probe(y_index), z_probe(z_index)];

                % se attivato rallenta di molto il processo di calcolo
                if PLOT3D_enabled
                    plotCube(LED1, LED2, LED3, LED4, PDect, x_max, y_max, z_max);
                end

                % costruisce le matrici degli impulsi e delle sovrapposizioni luminose
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

    pause(0.01);
    close all;
    if PLOT2D_enabled
        rel_min = 0;
        nan_idx = isnan(impulseMatrix);
        rel_max = mean(mean(mean(impulseMatrix(~nan_idx))));
        rel_max = rel_max * 2.5;
        
        % Z Sliced
        z_index = 1:length(z_probe);
        plot_index = linspace(1, max(z_index)-2, 9);
        
        figure('WindowState', 'maximized');
        colorbar;
        
        for i=1:9
            index = floor(plot_index(i));
            subplot(3, 3, i);
            matrix_slice = squeeze(impulseMatrix(:,:,index));
            imagesc(matrix_slice);
            colorbar;
            clim([rel_min rel_max]); % same scale
            xlabel('Y');
            ylabel('X');
            title('Z = '+string(z_probe(index))+' [m] - Element '+index);
        end

        colormap(jet(126));
        sgtitle('Z-sliced - Room size = ' + string(x_max) + ...
            'x' + string(y_max) + 'x' + string(z_max) + ' [m]');
        
        % Y Sliced
        y_index = 1:length(y_probe);
        plot_index = linspace(1, max(y_index)-2, 9);

        figure('WindowState', 'maximized');
        colorbar;

        for i=1:9
            index = floor(plot_index(i));
            subplot(3, 3, i);
            matrix_slice = squeeze(impulseMatrix(:,index,:));
            imagesc(matrix_slice);
            camroll(90);
            colorbar;
            clim([rel_min rel_max]); % same scale
            xlabel('Z');
            ylabel('X');
            title('Y = '+string(y_probe(index))+' [m] - Element '+index);
        end

        colormap(jet(126));
        sgtitle('Y-sliced - Room size = ' + string(x_max) + 'x' + ...
            string(y_max) + 'x' + string(z_max) + ' [m]');
        
        % X Sliced
        x_index = 1:length(x_probe);
        plot_index = linspace(1, max(x_index)-2, 9);

        figure('WindowState', 'maximized');
        colorbar;

        for i=1:9
            index = floor(plot_index(i));
            subplot(3, 3, i);
            matrix_slice = squeeze(impulseMatrix(index,:,:));
            imagesc(matrix_slice);
            camroll(90)
            colorbar;
            clim([rel_min rel_max]); % same scale
            xlabel('Z');
            ylabel('Y');
            title('X = '+string(x_probe(index))+' [m] - Element '+index);
        end

        colormap(jet(126));
        sgtitle('X-sliced - Room size = ' + string(x_max) + 'x' + ...
            string(y_max) + 'x' + string(z_max) + ' [m]');

        % overlapping
        z_index = 1:length(z_probe);

        plot_index = linspace(1, max(z_index)-2, number_of_plot_needed);
        figure('WindowState', 'maximized');
        colorbar;

        for i=1:number_of_plot_needed
            index = floor(plot_index(i));
            subplot(number_of_plot_needed/3, 3, i);
            matrix_slice = squeeze(overlap_info(:,:,index));
            imagesc(matrix_slice);
            colorbar
            clim([0 5]);   % same scale
            xlabel('Y');
            ylabel('X');
            title('Z = '+string(z_probe(index))+' [m]');
        end

        sgtitle('Overlapping MAP Z-sliced | Room size = ' + ...
            string(x_max) + 'x' + string(y_max) + 'x' + string(z_max) + ' [m]');
    end
end

%% Single Entity Distance Error Estimation
if plot_2
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

%% Noise Contribution
if plot_3
    shift_distance = 1; % meter

    x_probe = 0:granularity:x_max;
    y_probe = 0:granularity:y_max;
    z_probe = 0:4*granularity:z_max;

    x_len = length(x_probe);
    y_len = length(y_probe);
    z_len = length(z_probe);

    B_array = [5e6, 20e6, 100e6, 400e6];
    marker_array = ["x", "o", "^", "p"];

    for bandIdx = 1:length(B_array)
        B = B_array(bandIdx);

        pause(0.01);
        clc;
        disp("Starting to draw the noise contribution...");
        tic
        received_power_LED = NaN(1, z_len);
        n_shoot_LED_contribution = NaN(1, z_len);
        n_thermal_LED_contribution = NaN(1, z_len);
        distance_LED = zeros(1, z_len);

        max_rnd = number_of_samples;
        posLEDError = NaN(z_len, max_rnd);
        % par
        for z_index = 1:z_len
            PDect_pos_LED = [LED1(1), LED1(2), LED1(3) - z_probe(z_index)];
            received_power_LED(z_index) = R_pd * p * singleEntityContribution(...
                LED1, PDect_pos_LED, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of); % Ampere
            distance_LED(z_index) = calculateDistance(LED1, PDect_pos_LED);

            for idx_rnd = 1:max_rnd
                [n_shoot_LED, n_thermal_LED, nsh_var, nth_var] = noiseEstimation( ...
                    received_power_LED(z_index), q_0, R_pd, k_B, T_k, eta, I_2, I_3, ...
                    Gamma, A_pd, g_m, I_bg, G_0, B);
                n_shoot_LED_contribution(z_index) = nsh_var;
                n_thermal_LED_contribution(z_index) = nth_var;
                received_power_LED_w_noise = received_power_LED(z_index) ...
                    + n_shoot_LED + n_thermal_LED;
                dLED = getLEDDistanceByEstimatedPower(LED1(3) - PDect_pos_LED(3), ...
                    p, A_pd, received_power_LED_w_noise, Psi, T_of, Phi_FoV, a, R_pd);

                if isnan(dLED)
                    posLEDError(z_index, idx_rnd) = NaN;
                else
                    posLEDError(z_index, idx_rnd) = abs(distance_LED(z_index) - dLED);
                end
            end
        end

        %legend_array_led = [legend_array_led , "\sigma_s with B = "+string(B*1e-6)+" [MHz]", "\sigma_t with B = "+string(B*1e-6)+" [MHz]"];

        f = figure(100);
        f.WindowState = "maximized";
        %grid on;
        hold on;
        %title("Noise estimation LED1");
        ylabel('Signal Power [Ampere^2]');
        xlabel('Distance [m]');

        if bandIdx == 1
            %marker_array = ["x", "o", "^", "p"];
            plot(NaN, 'b-',  'LineWidth',2, 'MarkerSize', 15,'MarkerFaceColor','k');
            plot(NaN, 'kx',  'LineWidth',2, 'MarkerSize', 15,'MarkerFaceColor','k');
            plot(NaN, 'ko',  'LineWidth',2, 'MarkerSize', 15,'MarkerFaceColor','k');
            plot(NaN, 'k^',  'LineWidth',2, 'MarkerSize', 15,'MarkerFaceColor','k');
            plot(NaN, 'kp',  'LineWidth',2, 'MarkerSize', 15,'MarkerFaceColor','k');
            
            plot(distance_LED, received_power_LED.^2, "b-", 'LineWidth', 2);
            legend_array_led = ["LED", "B = 5 MHz", "B = 20 MHz", "B = 100 MHz", "B = 400 MHz"];
            legend(legend_array_led);
        end
        
        marks_indices = [2 1:12:length(distance_LED)];
        
        % plot(distance_LED,sqrt(n_shoot_LED_contribution), 'LineWidth',2);
        % plot(distance_LED,sqrt(n_thermal_LED_contribution), '-.', 'LineWidth',2);
        plot(distance_LED,(n_shoot_LED_contribution + n_thermal_LED_contribution), 'b--', ...
            'LineWidth', 2, 'Marker', marker_array(bandIdx), 'MarkerSize', 15, ...
            'MarkerFaceColor','b', 'MarkerIndices', marks_indices);
        
        % ylim([10e-25 10e-7]);
        % legend(legend_array_led);
        set(gca, 'YScale', 'log');
        set(gca, "FontName", "Arial", "FontSize", 25);
    end
    
    f = figure(101);
    f.WindowState = "maximized";
    hold on;
    plot(NaN, 'b-', 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', 'k');
    plot(NaN, 'kx', 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', 'k');
    plot(NaN, 'ko', 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', 'k');
    plot(NaN, 'k^', 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', 'k');
    plot(NaN, 'kp', 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor', 'k');
    legend(["LED", "B = 5MHz", "B = 20MHz", "B = 100MHz", "B = 400MHz"]);

    for bandIdx = 1:length(B_array)
        B=B_array(bandIdx);

        lumen_array = 1:10:10000;

        received_power_LED = zeros(1, length(lumen_array));
        distance_LED = zeros(1, length(lumen_array));
        SNR_LED = zeros(1, length(lumen_array));
        
        for lumenIdx = 1:length(lumen_array)
            lumen_level = lumen_array(lumenIdx);
            p = lumen_level / 683; % transmission power 6000 [Lumens] -> [Watt]

            PDect_pos_LED = [LED1(1), LED1(2), LED1(3) - shift_distance];
            received_power_LED(lumenIdx) = R_pd * p * singleEntityContribution( ...
                LED1, PDect_pos_LED, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of); % Ampere
            [n_shoot_LED, n_thermal_LED, nsh_var_LED, nth_var_LED] = noiseEstimation( ...
                received_power_LED(lumenIdx), q_0, R_pd, k_B, T_k, ...
                eta, I_2, I_3, Gamma, A_pd, g_m, I_bg, G_0, B);

            SNR_LED(lumenIdx) = 10 * log10(received_power_LED(lumenIdx).^2 ...
                / (nsh_var_LED+nth_var_LED));
        end

        marks_indices = 1:100:length(lumen_array);

        plot(lumen_array, SNR_LED, 'b--', 'Marker', marker_array(bandIdx), ...
            'MarkerIndices', marks_indices, 'MarkerSize', 15, 'MarkerFaceColor','b');
        ylabel('SNR [dB]');
        xlabel('Optical Power [Lumen]');
    end
    
    set(gca, "FontName", "Arial", "FontSize", 25);
    
    % f = figure(102);
    % f.WindowState = "maximized";
    % grid on;
    % hold on;
    % title("BoxPlot Error LED");
    % xlabel('Distanza [m]');
    % ylabel('Errore distanza [m]');
    % boxchart(posLEDError', 'MarkerStyle','none');
    % ylim([0 y_max/6]);
    % ax = gca;
    
    toc;
end

%% Estimate Position
pause(0.01);
fprintf("Starting to minimize position error...\n");
tic

%     lumen_level = 1000;
%     p = lumen_level / 683;  % transmission power 6000 [Lumens] -> [Watt]
%     B = 100e6;              % System bandwidth [Hz]

number_of_plot_needed = 1;

if plot_4
    x_probe = 0:granularity:x_max;
    y_probe = 0:granularity:y_max;
    z_probe = 0:9*granularity:z_max;

    % number_of_plot_needed = 18;
    % z_plot_index = linspace(0, max(z_probe), number_of_plot_needed);
    z_plot_index = [0]; %#ok<NBRAK2>

    f = figure(300);
    f.WindowState = "maximized";
    colorbar;

    % for k = 1:length(z_probe)
    % fprintf("Sample n: %d", number_of_samples);

    for k=1:number_of_plot_needed
        error_matrix = zeros(length(x_probe), length(y_probe));

        for i=1:length(x_probe)
            %par
            for j=1:length(y_probe)
                PD_to_findx = [x_probe(i) y_probe(j) z_plot_index(k)];
                [p1, p2, p3] = estimateReceiverPosition(LED1, LED2, LED3, LED4, PD_to_findx,...
                    p, alpha, beta, Phi_FoV, a, Psi, A_pd, T_of, R_pd, q_0, k_B, T_k, eta, I_2, ...
                    I_3, Gamma, g_m, I_bg, G_0, B, K_0);
                if isnan(p1) || isnan(p2) || isnan(p3)
                    error_matrix(i,j) = NaN;
                else
                    error_matrix(i,j) = real(sqrt((p1 - PD_to_findx(1))^2 + (p2 - PD_to_findx(2))^2 + (p3 - PD_to_findx(3))^2)); % RMSE
                end
            end
        end
        
        if number_of_plot_needed > 1
            subplot(number_of_plot_needed/3, 3, k);
        else
            subplot(1,1,1);
        end
        
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

        title('Z = '+string((z_plot_index(k)))+' [m] - mean = ' + mean_error_string);
    end

    t = sgtitle('Position error Z-sliced - Room size = ' + string(x_max) + ...
        'x' + string(y_max) + 'x' + string(z_max) + ' [m]');
    t.Margin = 5;
end
