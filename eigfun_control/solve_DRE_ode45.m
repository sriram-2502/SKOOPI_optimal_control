function [P_list, s_x_list, K_ff_list, K_fb_list, t_grid] = solve_DRE_ode45(x_ref, A, B_list, Q, R, QN, T, dt)
    % solve_DRE using ode45 integration for continuous-time finite-horizon LQR tracking

    n = size(A, 1);
    t_grid = 0:dt:T;     % Forward time grid
    t_span = [T, 0];         % Integrate backward

    % Initial conditions (terminal values)
    P0 = QN;
    s0 = -QN * x_ref(:, end);
    Z0 = [P0(:); s0];

    % Combine dynamics for ode45
    ode_rhs = @(t, Z) riccati_ode(t, Z, A, B_list, x_ref, Q, R, t_grid);
    opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
    [t_sol, Z_sol] = ode45(ode_rhs, t_span, Z0, opts);

    % Flip to forward time
    t_sol = flip(t_sol);
    Z_sol = flip(Z_sol, 1);

    % Extract solutions
    N = length(t_sol);
    P_list = cell(1, N);
    s_x_list = cell(1, N);
    K_fb_list = cell(1, N);
    K_ff_list = cell(1, N);

    for i = 1:N
        Z_i = Z_sol(i, :)';
        P = reshape(Z_i(1:n^2), n, n);
        s_x = Z_i(n^2+1:end);

        if iscell(B_list)
            B = B_list{i};
        else
            B = B_list;
        end

        P_list{i} = P;
        s_x_list{i} = s_x;
        K_fb_list{i} = R \ (B' * P);
        K_ff_list{i} = R \ (B' * s_x);
    end
end

function dZdt = riccati_ode(t, Z, A, B_list, x_ref, Q, R, t_grid)
    n = size(A,1);
    P = reshape(Z(1:n^2), n, n);
    s_x = Z(n^2+1:end);

    % Interpolate reference trajectory
    x_d = interp1(t_grid', x_ref', t, 'linear', 'extrap')';

    % Interpolate B (if varying with time)
    if iscell(B_list)
        % Simple linear index for now (could improve with interpolation)
        idx = round(interp1(t_grid, 1:length(t_grid), t, 'nearest', 'extrap'));
        B = B_list{max(1,min(length(B_list),idx))};
    else
        B = B_list;
    end

    % Riccati equations
    dP = -(Q - P*B*(R\B')*P + P*A + A'*P);
    ds_x = -(-Q*x_d + (A' - P*B*(R\B'))*s_x);

    % Return stacked derivatives
    dZdt = [dP(:); ds_x];
end
