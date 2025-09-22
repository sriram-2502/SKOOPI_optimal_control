function x_next = rk4(dynamics, dt_sim, x_op, u, use_reverse, sys_info)
    % Runge-Kutta 4th Order (RK4) integration method
    % Inputs:
    %   dynamics - function handle for system dynamics, returning dxdt
    %   dt_sim - time step for integration
    %   x_op - current state vector
    %   u - control input vector
    %
    % Output:
    %   x_next - next state vector after time step dt_sim

    % Calculate the four increments (k1, k2, k3, k4)
    if(use_reverse)
        k1 = -dynamics(x_op, u, sys_info);
        k2 = -dynamics(x_op + 0.5 * dt_sim * k1, u, sys_info);
        k3 = -dynamics(x_op + 0.5 * dt_sim * k2, u, sys_info);
        k4 = -dynamics(x_op + dt_sim * k3, u, sys_info);
        % Compute the next state
        x_next = x_op + (dt_sim / 6) * (k1 + 2*k2 + 2*k3 + k4);
    else
        k1 = dynamics(x_op, u, sys_info);
        k2 = dynamics(x_op + 0.5 * dt_sim * k1, u, sys_info);
        k3 = dynamics(x_op + 0.5 * dt_sim * k2, u, sys_info);
        k4 = dynamics(x_op + dt_sim * k3, u, sys_info);
        % Compute the next state
        x_next = x_op + (dt_sim / 6) * (k1 + 2*k2 + 2*k3 + k4);
    end        

end