function animation = animate_quadrotor(xout, uout, sys_info)
% x: state vector [x; y; theta; vx; vy; omega]
% u: control vector [u1; u2]

% parse params
l = sys_info.l; 
ax_limt = 1;

LENGTH = l;
theta = xout(3);
pos = [xout(1); xout(2)];

% Rotation matrix
rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];

% Define the quadrotor parts
main_frame    = rot * [-LENGTH, LENGTH; 0, 0] + repmat(pos,1,2);
left_propeller  = rot * [-1.3*LENGTH, -0.7*LENGTH; 0.1, 0.1] + repmat(pos,1,2);
right_propeller = rot * [1.3*LENGTH, 0.7*LENGTH; 0.1, 0.1] + repmat(pos,1,2);
left_thrust   = rot * [LENGTH, LENGTH; 0.1, 0.1 + uout(1)*0.04] + repmat(pos,1,2);
right_thrust  = rot * [-LENGTH, -LENGTH; 0.1, 0.1 + uout(2)*0.04] + repmat(pos,1,2);

% Plot the quadrotor
clf; % clear current figure
hold on; grid on; box on;
axis([-ax_limt ax_limt -ax_limt ax_limt]);
plot(main_frame(1,:), main_frame(2,:), 'k', 'LineWidth', 6);
plot(left_propeller(1,:), left_propeller(2,:), 'b', 'LineWidth', 4);
plot(right_propeller(1,:), right_propeller(2,:), 'b', 'LineWidth', 4);
plot(left_thrust(1,:), left_thrust(2,:), 'r', 'LineWidth', 1);
plot(right_thrust(1,:), right_thrust(2,:), 'r', 'LineWidth', 1);
drawnow;

end