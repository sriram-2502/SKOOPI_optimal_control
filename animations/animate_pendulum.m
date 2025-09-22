function animation = animate_pendulum(theta, sys_info)
    L = sys_info.L;
    % Gray line for upright position
    plot([0, 0], [0, 1.5], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.5); hold on;
    plot(0, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'black'); % Pivot point
    x_pendulum = -L * sin(theta);
    y_pendulum = L * cos(theta);
    plot([0, x_pendulum], [0, y_pendulum], 'k-', 'LineWidth', 2); % Pendulum rod
    plot(x_pendulum, y_pendulum, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'red'); % Pendulum bob
    xlim([-1.5, 1.5]);
    ylim([-1.5, 1.5]);
    axis square;       % Make the axis square
    grid on; 
    box on;
    xticks([]);
    yticks([]);
    title('Pendulum Animation');
    box on
    set(gca, 'LineWidth', 2);
    set(gca, 'FontSize', 20);