function [robot_params,euler_params,navigation_params,animate_params] = get_params_planarRR(x)

if(nargin<1)
    x = sym('x',2);
end
%% robotic arm parameters
% robot arm lengths
robot_params.l1 = 1; robot_params.l2 = 1;
robot_params.m1 = 1; robot_params.m2 = 1;
robot_params.g = 9.81;
robot_params.I1 = 1 ;robot_params.I2 = 1;
robot_params.D = diag([1 1]);

%% robot end effector paramters
% these params are used for normalizing the control input from the planner
% normalize u as max(u)/mu*mg
% has a saturation effect on the controller
% leads to more smooth trajectories
%navigation_params.mu = 0.1; % Friction coefficient [-]
%navigation_params.mg = 20; % Weight [N].

% params for motion plan only
navigation_params.saturation = 2; %(best value = 2 for default ID)
navigation_params.rad_from_goal = 0.01; %(best value = 0.01 for default ID)
navigation_params.ctrl_multiplier = 1; %(best value = 10 for default ID)

%% parameters for euler integration
% Target set radius | Radius when to stop using Density FB control

euler_params.n_steps = 500;
euler_params.step_size = 0.01;


%% animation paramters
animate_params.flag_movie = 1;
animate_params.skip_rate = 10;

% params for dynamic obs case
animate_params.dynamic_obs.traj = 0;
animate_params.dynamic_obs.density_gif = 0;

% params for static obs case
animate_params.static_obs = 1;

% params for plotting obstalce in joint space
animate_params.plot_joint_obs = 1;