function u = compute_control(lqr_params_transformed, P_riccati_curr, phi_x_op, grad_phi_x_op)
    
    % check if system is linear or not
    if(isempty(grad_phi_x_op))
        use_linear = true;
    else
        use_linear = false;
    end

    % parse system info
    g = lqr_params_transformed.B; %lienarized version of g(x) in transformed coordinates
%     g = sys_info.dynamics_g(x_op');

    % parse lqr params
    R = lqr_params_transformed.R;

    % get control
    if(use_linear)
        % if linear, no need for grad_phi
        u = -inv(R)*g'*(P_riccati_curr*phi_x_op);
    else
        u = -inv(R)*g'*(phi_x_op*P_riccati_curr*grad_phi_x_op)';
    end

    
end