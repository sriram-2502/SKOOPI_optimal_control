function x_next = euler(dynamics, dt_sim, x_op, u, use_reverse, sys_info)
    % Perform Euler integration to compute the next state
    % 
    % Inputs:
    %   dynamics   - Function handle that returns dxdt given x_op and u
    %   dt_sim     - Time step for the simulation
    %   x_op       - Current state vector (n x 1)
    %   u          - Control input (could be a vector, depending on the dynamics)
    % 
    % Output:
    %   x_next     - Next state vector (n x 1)
    
    if(use_reverse)
        % Compute the rate of change of the state
        dxdt = -dynamics(x_op, u, sys_info);  % use -f(x) for reverse time simulation
    else
        % Compute the rate of change of the state
        dxdt = dynamics(x_op, u, sys_info);  % Call the dynamics function to get dxdt
    end
    
    % Update the state using Euler integration
    x_next = x_op + dt_sim * dxdt;  % Euler method: x_next = x_op + dt * dxdt
end