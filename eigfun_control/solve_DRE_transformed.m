function [L_list, s_x_list, K_ff_list, K_fb_list, t_grid] = ...
    solve_DRE_transformed(Lambda, u_ref, phi_list, D_list, D_list2, lqr_params_transformed, T, dt)
% Finite-horizon transformed DRE (tracking) via ode23s, backward in time.
% Returns L(t), s(t), and tracking gains K_fb^ψ(t), K_ff^ψ(t).

    %% parse LQR params
    Q   = lqr_params_transformed.Q;
    R   = lqr_params_transformed.R;
    QN  = lqr_params_transformed.QN;
    sxN = -QN * phi_list{end};

    % time grid (descending)
    t_grid = T-dt:-dt:0;       % length N
    N      = numel(t_grid);

    % sizes
    n  = size(QN,1);
    m  = size(R,1);
    nL = n*n;

    % terminal conditions (t = T)
    L_T = QN;
    s_T = sxN;
    y0  = [L_T(:); s_T];

    % handle complex Lambda (your helper if available)
    if ~isreal(Lambda)
        if exist('replace_complex','file') == 2
            Lam_vec  = replace_complex(diag(Lambda));
            Lambda_R = diag(Lam_vec);
        else
            Lambda_R = diag(real(diag(Lambda)));
        end
    else
        Lambda_R = Lambda;
    end

    % index mapping t -> k (left-nearest sample)
    Nt_phi = numel(phi_list);
    clamp  = @(k,Nmax) max(1, min(Nmax, k));
    t2k    = @(t) clamp(floor(t/dt) + 1, Nt_phi);

    % precompute R^{-1} (minimal change)
    Rinv = inv(R);

    %% RHS for ode23s: y = [vec(L); s]
    function dy = dre_rhs(t, y)
        L = reshape(y(1:nL), n, n);  L = 0.5*(L + L.');
        s = y(nL+1:nL+n);

        kk    = t2k(t);
        phi_d = phi_list{kk};
        D_t   = D_list{kk};
        D_t2  = D_list2{kk};
        u_d   = u_ref(:, kk);

        % costs per step
        Q_L  = Q;
        Q_sx = QN * phi_d;

        % DRE (your form): -dL/dt = LΛ + Λ' L - 2 L D L + Q_L
        dL = -( L*Lambda_R + Lambda_R'*L - 2*L*D_t*L + Q_L );

        % costate: -ds/dt = (Λ' - 2 L D) s - Q_sx + L D2 u_d
        ds = - ( Lambda_R'*s - 2*L*D_t*s - Q_sx ) + L*D_t2*u_d;

        dy = [dL(:); ds];
    end

    %% integrate backward: T -> 0
    opts = odeset('RelTol',1e-7,'AbsTol',1e-9);
    sol  = ode23s(@dre_rhs, [T 0], y0, opts);

    %% sample on t_grid (descending)
    Yg = deval(sol, t_grid);   % (n^2+n) x N

    % unpack
    L_list   = cell(1, N);
    s_x_list = cell(1, N);
    for k = 1:N
        Lk = reshape(Yg(1:nL, k), n, n);  Lk = 0.5*(Lk + Lk.');
        sk = Yg(nL+1:nL+n, k);
        L_list{k}   = Lk;
        s_x_list{k} = sk;
    end

    %% gains on the grid (tracking form)
    K_fb_list = cell(1, N-1);
    K_ff_list = cell(1, N-1);
    for k = 1:N-1
        kk    = t2k(t_grid(k));
        Lk    = L_list{k};
        sk    = s_x_list{k};
        D2k   = D_list2{kk};
        phi_d = phi_list{kk};
        u_d   = u_ref(:, kk);

        Kfb = Rinv * D2k' * Lk;                    % m x n
        Kff = u_d + Rinv * D2k' * (sk - Lk*phi_d); % m x 1

        K_fb_list{k} = Kfb;
        K_ff_list{k} = Kff;
    end

    %% flip to ascending time
    L_list    = flip(L_list);
    s_x_list  = flip(s_x_list);
    K_fb_list = flip(K_fb_list);
    K_ff_list = flip(K_ff_list);
    t_grid    = flip(t_grid);
end
