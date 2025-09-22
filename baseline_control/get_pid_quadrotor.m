function [u_pid,u_ref] = get_pid_quadrotor(x_des, x_op)
    % PID controlf or tracking a 2D quadrotor
    % x_des = y_des, z_des, vy_des, vz_des, ay_des, az_des

    % check if x is symbolic (for ref control
    is_symbolic = isa(x_op, 'sym');


    % get system info
    sys_info = quadrotor_info();

    % parse info
    m = sys_info.m;
    l = sys_info.l; 
    g = sys_info.g;
    I = sys_info.I;

    % parse x_des
    y_des        = x_des(1);
    z_des        = x_des(2);
    y_dot_des    = x_des(3);
    z_dot_des    = x_des(4);
    y_ddot_des   = x_des(5);
    z_ddot_des   = x_des(6);

    % parse xop
    y           = x_op(1);
    z           = x_op(2);
    theta       = x_op(3);
    y_dot       = x_op(4);
    z_dot       = x_op(5);
    theta_dot   = x_op(6);


    % setup gains
    Kp_y     = 0.4;
    Kv_y     = 1.0;
    Kp_z     = 0.4;
    Kv_z     = 1.0;
    Kp_theta = 18;
    Kv_theta = 15;

    % get control 
    theta_des = -1/g * (y_ddot_des + Kv_y * (y_dot_des - y_dot) + Kp_y * (y_des - y));
    F = m * (g + z_ddot_des + Kv_z * (z_dot_des - z_dot) + Kp_z * (z_des - z));
    T = I * (Kv_theta * (-theta_dot) + Kp_theta * (theta_des - theta));
    u_ref = [F, T];
    
    if(~is_symbolic)
        % clamp control 
        [F_clamped, T_clamped] = clamp_ctrl_quadrotor(F, T, l);
    
        % parse control
        u_pid = [F_clamped,T_clamped];
    else
        u_pid = u_ref;
    end

end



