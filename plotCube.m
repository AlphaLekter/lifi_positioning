
function plotCube(LED, RIS1, RIS2, RIS3, RIS4, PD, x, y, z)
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
% Plot the LED and RIS as red points
scatter3(LED(1), LED(2), LED(3), 100, 'r', 'filled'); % LED
scatter3(RIS1(1), RIS1(2), RIS1(3), 100, 'r', 'filled'); % RIS1
scatter3(RIS2(1), RIS2(2), RIS2(3), 100, 'r', 'filled'); % RIS2
scatter3(RIS3(1), RIS3(2), RIS3(3), 100, 'r', 'filled'); % RIS3
scatter3(RIS4(1), RIS4(2), RIS4(3), 100, 'r', 'filled'); % RIS3
scatter3(PD(1), PD(2), PD(3), 100, 'b', 'filled'); % PD

text(LED(1), LED(2), LED(3)+0.3, 'LED', 'FontSize', 14, 'HorizontalAlignment','center');
text(RIS1(1), RIS1(2), RIS1(3)+0.3, 'RIS1', 'FontSize', 14, 'HorizontalAlignment','center');
text(RIS2(1), RIS2(2), RIS2(3)+0.3, 'RIS2', 'FontSize', 14, 'HorizontalAlignment','center');
text(RIS3(1), RIS3(2), RIS3(3)+0.3, 'RIS3', 'FontSize', 14, 'HorizontalAlignment','center');
text(RIS4(1), RIS4(2), RIS4(3)+0.3, 'RIS3', 'FontSize', 14, 'HorizontalAlignment','center');
text(PD(1), PD(2), PD(3)+0.3, 'PD', 'FontSize', 14, 'HorizontalAlignment','center');

plot3([LED(1), RIS1(1)], [LED(2), RIS1(2)], [LED(3), RIS1(3)],'green', 'LineWidth', 2); % line between LED / RIS1
plot3([LED(1), RIS2(1)], [LED(2), RIS2(2)], [LED(3), RIS2(3)],'k', 'LineWidth', 2); % line between LED / RIS2
plot3([LED(1), RIS3(1)], [LED(2), RIS3(2)], [LED(3), RIS3(3)],'b', 'LineWidth', 2); % line between LED / RIS3
plot3([LED(1), RIS4(1)], [LED(2), RIS4(2)], [LED(3), RIS4(3)],'y', 'LineWidth', 2); % line between LED / RIS3

plot3([PD(1), RIS1(1)], [PD(2), RIS1(2)], [PD(3), RIS1(3)],  'green','LineWidth', 2); % line between PD / RIS1
plot3([PD(1), RIS2(1)], [PD(2), RIS2(2)], [PD(3), RIS2(3)],  'k','LineWidth', 2); % line between PD / RIS2
plot3([PD(1), RIS3(1)], [PD(2), RIS3(2)], [PD(3), RIS3(3)],  'b','LineWidth', 2); % line between PD / RIS2
plot3([PD(1), RIS4(1)], [PD(2), RIS4(2)], [PD(3), RIS4(3)],  'y','LineWidth', 2); % line between PD / RIS2

% Set the limits of the plot
xlim([0 x]);
ylim([0 y]);
zlim([0 z]);

view(27,18);
pause(eps);
end