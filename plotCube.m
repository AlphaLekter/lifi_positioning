
function plotCube(LED1, LED2, LED3, LED4, PD, x, y, z)
pause(eps);
% Set up the figure
figure(1);
clf;
hold on;
grid on;
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');

%set(figure(1), 'Position', [100 100 800 600]);
% Plot the LED as red points
scatter3(LED1(1), LED1(2), LED1(3), 100, 'r', 'filled'); % LED1
scatter3(LED2(1), LED2(2), LED2(3), 100, 'r', 'filled'); % LED2
scatter3(LED3(1), LED3(2), LED3(3), 100, 'r', 'filled'); % LED3
scatter3(LED4(1), LED4(2), LED4(3), 100, 'r', 'filled'); % LED3
scatter3(PD(1), PD(2), PD(3), 100, 'b', 'filled'); % PD

text(LED1(1), LED1(2), LED1(3)+0.3, 'LED1', 'FontSize', 14, 'HorizontalAlignment','center');
text(LED2(1), LED2(2), LED2(3)+0.3, 'LED2', 'FontSize', 14, 'HorizontalAlignment','center');
text(LED3(1), LED3(2), LED3(3)+0.3, 'LED3', 'FontSize', 14, 'HorizontalAlignment','center');
text(LED4(1), LED4(2), LED4(3)+0.3, 'LED3', 'FontSize', 14, 'HorizontalAlignment','center');
text(PD(1), PD(2), PD(3)+0.3, 'PD', 'FontSize', 14, 'HorizontalAlignment','center');

plot3([PD(1), LED1(1)], [PD(2), LED1(2)], [PD(3), LED1(3)],  'green','LineWidth', 2); % line between PD / LED1
plot3([PD(1), LED2(1)], [PD(2), LED2(2)], [PD(3), LED2(3)],  'k','LineWidth', 2); % line between PD / LED2
plot3([PD(1), LED3(1)], [PD(2), LED3(2)], [PD(3), LED3(3)],  'b','LineWidth', 2); % line between PD / LED2
plot3([PD(1), LED4(1)], [PD(2), LED4(2)], [PD(3), LED4(3)],  'y','LineWidth', 2); % line between PD / LED2

% Set the limits of the plot
xlim([0 x]);
ylim([0 y]);
zlim([0 z]);

view(27,18);
pause(eps);
end