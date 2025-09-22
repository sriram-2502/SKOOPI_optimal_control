function lqr_params = get_lqr(A,B,Q,R)
    % inputs
    % A - linearization A matrix
    % B - linearization B matrix
    % Q - state cost
    % R - control cost
    %
    % output
    % lqr_params (struct)
    
    % init struct
    lqr_params = struct();
    
    % calculation of LQR gain
    [K_lqr,P_lqr,e] = lqr(A,B,Q,R);
    
    % parse outputs
    lqr_params.A            = A;
    lqr_params.B            = B;
    lqr_params.Q            = Q;
    lqr_params.Q_terminal   = P_lqr;
    lqr_params.R            = R;
    lqr_params.K_lqr        = K_lqr;
    lqr_params.P_lqr        = P_lqr;
end
