function [P_list, s_x_list, K_ff_list, K_fb_list,  t_grid] = solve_DRE_tracking(x_ref, u_ref, A, B_list, Q, R, QN, T, dt)
% Finite-horizon continuous-time LQT via stiff ODE (ode23s), backward in time.
% Minimal changes from your original backward-Euler version.

    % ---- detect B input type ----
    if iscell(B_list), use_B_list = true; else, use_B_list = false; end

    % ---- sampling grid (matches your main script loop dt:dt:T) ----
    % N samples total (no enforced N+1); we’ll evaluate exactly on these points.
    t_grid = T-dt:-dt:0;         % descending
    N = numel(t_grid);

    % ---- sizes & terminal conditions ----
    n  = size(A,1);
    nP = n*n;

    P_T = QN;
    s_T = -QN * x_ref(:, end);   % consistent with your code
    y0  = [P_T(:); s_T];

    % ---- precompute stable solver for R (kept minimal: inv) ----
    Rinv = inv(R);

    % ---- index mapping t -> k for refs and B ----
    % Works whether refs/B have length N or N+1; uses left-nearest sample.
    Nt_ref = size(x_ref,2);
    idx_clamp = @(k) max(1, min(Nt_ref, k));
    t2k = @(t) idx_clamp(floor(t/dt) + 1);

    % ---- RHS of coupled ODEs (DRE + costate) ----
    function dy = dre_rhs(t, y)
        % unpack & symmetrize P
        P = reshape(y(1:nP), n, n);
        P = 0.5*(P + P.');
        s = y(nP+1:nP+n);

        % refs & B at current time
        k   = t2k(t);
        x_d = x_ref(:, k);
        u_d = u_ref(:, k);
        if use_B_list, B = B_list{k}; else, B = B_list; end

        % -dP/dt = A'P + PA - PBR^{-1}B'P + Q
        dP = -(Q - P*B*(Rinv*B')*P + P*A + A'*P);

        % -ds/dt = (A' - PBR^{-1}B') s + Q x_d - P B u_d
        ds = -((A' - P*B*Rinv*B')*s + Q*x_d - P*B*u_d);

        dy = [dP(:); ds];
    end

    % ---- stiff integration backward: T -> 0 ----
    opts = odeset('RelTol',1e-7,'AbsTol',1e-9);
    sol  = ode23s(@dre_rhs, [T 0], y0, opts);

    % ---- sample solution exactly on your t_grid (descending) ----
    Yg = deval(sol, t_grid);     % (n^2+n) x N

    % ---- unpack lists (still descending) ----
    P_list   = cell(1, N);
    s_x_list = cell(1, N);
    for k = 1:N
        Pk = reshape(Yg(1:nP, k), n, n);
        Pk = 0.5*(Pk + Pk.');               % enforce symmetry
        sk = Yg(nP+1:nP+n, k);

        P_list{k}   = Pk;
        s_x_list{k} = sk;
    end

    % ---- gains on grid (your original formulas) ----
    K_fb_list = cell(1, N-1);
    K_ff_list = cell(1, N-1);
    for k = 1:N-1
        kk = t2k(t_grid(k));
        if use_B_list, Bk = B_list{kk}; else, Bk = B_list; end
        Pk = P_list{k};
        sk = s_x_list{k};

        x_d = x_ref(:,kk);              
        u_d = u_ref(:,kk);               
        K_fb = Rinv * Bk' * Pk;
        K_ff = u_d + Rinv * Bk' * (sk - Pk * x_d);
        K_fb_list{k} = K_fb;
        K_ff_list{k} = K_ff;
    end

    % ---- flip to ascending time to match prior convention ----
    P_list    = flip(P_list);
    s_x_list  = flip(s_x_list);
    K_fb_list = flip(K_fb_list);
    K_ff_list = flip(K_ff_list);
    t_grid    = flip(t_grid);
end
