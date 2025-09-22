function u_energy = get_energy_based(x_op)
    
    % get system info
    sys_info = cart_pole_info();

    % parse states
    x           = x_op(1);
    theta       = x_op(2);
    x_dot       = x_op(3);
    theta_dot   = x_op(4);

    % prase system params
    M   = sys_info.M;   % Mass of cart
    m   = sys_info.m;   % mass of pendulum
    I   = sys_info.I;   % MOI of Pendulum
    l   = sys_info.l;   % COM of Pendulum
    g   = sys_info.g;   % Gravity Constant
    b   = sys_info.b;   % viscous damping at pivot of Pendulum

    %% setup params for energy shaping
    % Total energy
    E = m*g*l*(1-cos(theta)) + (1/2)*(I + m*l^2)*(theta_dot^2);
    
    % Potential Energy
    Er  = 2*m*g*l;

    %% setup params for swing control
    n       = sys_info.n; 
    k_swing = sys_info.k_swing;

    % Energy based swing up control 
    accel = 2*(E-Er)*sign(theta_dot*cos(theta));
    accel = k_swing*g*(saturate_fun(accel, n*g, -n*g));

    % feedback Linearization
    u_energy = (M+m)*(accel)+ 0*x_dot-m*l*( (theta_dot)^2)*sin(theta)- m*l*(cos(theta))*( ( b*theta_dot + m*l*accel*cos(theta) + m*g*l*sin(theta) )/(I+m*l^2) );
end