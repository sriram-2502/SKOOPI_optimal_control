function [t_riccati,P_riccati] = compute_riccati(lqr_params_transfomed, t_span)
    
    A_transformed    = lqr_params_transfomed.A;
    B_transformed    = lqr_params_transfomed.B;
    Q_transformed    = lqr_params_transfomed.Q;
    R_transformed    = lqr_params_transfomed.R;
    Q_T_transformed  = lqr_params_transfomed.Q_terminal;
    
    %% Solve (numerically) for Finite-horizon LQR
    t_span_rev = flipud(t_span');
    P_initial = Q_T_transformed;
    
    options = odeset('JConstant','on', 'RelTol',1e-6, 'AbsTol',1e-6);
    [t_riccati,P_riccati]=ode45(@(t,P)riccati_ode(t,P,A_transformed,B_transformed,Q_transformed,R_transformed),t_span_rev,P_initial,options);
    P_riccati = flipud(P_riccati);
    
end