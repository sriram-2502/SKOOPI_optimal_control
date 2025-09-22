function [F_clamped, M_clamped] = clamp_ctrl_quadrotor(F, T, l)
    % Calculate initial motor inputs
    u1 = 0.5 * (F - T / l);
    u2 = 0.5 * (F + T / l);
    
    % Check for motor saturation
    if u1 < 0 || u1 > 1.7658 || u2 < 0 || u2 > 1.7658
        fprintf('Motor saturation: u1 = %.4f, u2 = %.4f\n', u1, u2);
    end
    
    % Clamp motor inputs
    u1_clamped = min(max(0, u1), 1.7658);
    u2_clamped = min(max(0, u2), 1.7658);
    
    % Calculate clamped force and moment
    F_clamped = u1_clamped + u2_clamped;
    M_clamped = (u2_clamped - u1_clamped) * l;
end