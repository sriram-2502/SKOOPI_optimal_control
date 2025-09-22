clc; clear; close all

% ---------- default plot options ----------
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultAxesFontSize',18)
set(0,'defaultfigurecolor',[1 1 1])

% ---------- add paths (expects your existing folders/functions) ----------
addpath('dynamics', 'baseline_control', 'eigfun_control', ...
        'compute_eigfuns', 'utils', 'animations');
rng(0)

% ---------- system setup ----------
sys_params.use_stable         = false;   % locally stable
sys_params.use_unstable       = true;    % locally unstable
sys_params.use_linear_riccati = false;

% Load dynamics & system info
sys_info   = nonlinear_sys_info(sys_params);
dynamics   = @dynamics_nonlinear;
n_states   = sys_info.state_dim;
n_ctrl     = sys_info.ctrl_dim;

A = sys_info.A;
B = sys_info.B;

% Analytical eigenvector matrix at origin (your helper)
W = sys_info.grad_eigen_fun_analytical([0;0]);
sys_info.eig_vectors = W;

% ---------- LQR (baseline, infinite-horizon) ----------
Q_baseline = 1e2*eye(n_states);
R_baseline = 1e1*eye(n_ctrl);
lqr_params_baseline = get_lqr(A, B, Q_baseline, R_baseline);

% ---------- SKOOPI / psi coordinates (your current approach uses A,B) ----------
Q  = Q_baseline;
R  = R_baseline;
QN = 1e3*Q_baseline;

A_transformed = A;
B_transformed = B;
Q_transformed = Q_baseline;
lqr_params_transformed = get_lqr(A_transformed, B_transformed, Q_transformed, R);

% ---------- simulation timing ----------
dt_sim = 0.01;
t_end  = 5.0;
t_grid = 0:dt_sim:t_end;         % 1 x N
Nsteps = numel(t_grid);

% ---------- Precompute finite-horizon gains (LQR-FT) ----------
x_desired = zeros(n_states,1);
g_list    = repmat({B}, 1, Nsteps-1);           % constant B each step
x_ref_mat = zeros(n_states, Nsteps-1);          % zero reference
u_ref_mat = zeros(n_ctrl,   Nsteps-1);          % zero reference control

% Solve DRE in original coordinates (feedback + feedforward)
[P_list, s_x_list2, K_ff_list, K_fb_list, t_grid_dre] = ...
    solve_DRE(x_ref_mat, u_ref_mat, A, g_list, Q, R, QN, t_end, dt_sim);

% ---------- Precompute Riccati solution for SKOOPI (time-varying P) ----------
[~, P_riccati] = compute_riccati(lqr_params_transformed, t_grid);
% P_riccati: Nsteps x (n_states^2); row i is vec(P_i)

% ---------- IC sampling (uniform in a disk around (1,2)) ----------
N_ic   = 100;
center = [-2; 3];
radius = 1;

theta = 2*pi*rand(N_ic,1);
rad   = radius*sqrt(rand(N_ic,1));
ICs   = center.' + [rad.*cos(theta), rad.*sin(theta)];  % N_ic x 2

% ---------- arrays to collect time traces for each method ----------
X1_all = cell(N_ic,1);   % LQR states
X2_all = cell(N_ic,1);   % SKOOPI states
X3_all = cell(N_ic,1);   % LQR-FT states

U1_all = cell(N_ic,1);   % LQR controls (Nsteps-1 x n_ctrl)
U2_all = cell(N_ic,1);   % SKOOPI controls
U3_all = cell(N_ic,1);   % LQR-FT controls

% ---------- figure & colors for phase portrait ----------
if ~ishandle(1)
    figure(1);  % create if it doesn't exist
end
figure(1);  % make fig.1 current

axp = subplot(4,4,[1 2 5 6]);  % phase panel
hold(axp,'on'); box(axp,'on'); grid(axp,'on');

xlabel(axp,'state, $x_1$','Interpreter','latex');
ylabel(axp,'state, $x_2$','Interpreter','latex');

% Colors & alpha
C = axp.ColorOrder;
alpha_line = 0.5;

% Re-plot disk & center
tt = linspace(0,2*pi,200);
plot(axp, center(1)+radius*cos(tt), center(2)+radius*sin(tt), ':k');
plot(axp, center(1), center(2), 'ok', 'MarkerFaceColor','k');

% ---------- per-IC simulation ----------
for k = 1:N_ic
    % initial states for each controller path
    x1 = ICs(k,:);   % LQR
    x2 = x1;         % SKOOPI
    x3 = x1;         % LQR-FT

    % logs for this IC (states)
    X1 = x1; X2 = x2; X3 = x3;

    % logs for this IC (controls)
    U1 = []; U2 = []; U3 = [];

    iter = 1;  % time index for time-varying gains
    % simulate forward in time
    for t_sim = dt_sim:dt_sim:t_end
        % --- Controller 1: infinite-horizon LQR ---
        u1 = -lqr_params_baseline.K_lqr * x1';        % (n_ctrl x 1)

        % --- Controller 3: finite-horizon LQR (feedback + feedforward) ---
        idx = min(iter, numel(K_fb_list));
        K_fb = K_fb_list{idx};
        K_ff = K_ff_list{idx};
        u3   = -K_fb * x3' - K_ff;

        % --- Controller 2: SKOOPI (time-varying quadratic form in lifted coords) ---
        psi_x_anal = x2' + inv(W') * sys_info.transform_fun_analytical(x2');
        idp   = min(iter, size(P_riccati,1));
        Pk    = reshape(P_riccati(idp,:), n_states, n_states);
        u2    = -inv(R)*B'*(Pk*psi_x_anal);

        % --- Controller 4: finite-horizon SKOOPI --- 
        psi_x_anal = x2' + inv(W') * sys_info.transform_fun_analytical(x2'); 
        idx = min(iter, numel(K_fb_list)); 
        K_fb = K_fb_list{idx}; 
        K_ff = K_ff_list{idx}; 
        u2 = -K_fb * (psi_x_anal) - K_ff;

        % switch to LQR when close / invalid
        if ~isfinite(u2) || any(abs(u2)>1e3) || norm(x2) <= 1e-3
            u2 = -lqr_params_baseline.K_lqr * x2';
        end

        % --- integrate one step (your RK4) ---
        use_reverse = false;
        x1n = rk4(dynamics, dt_sim, x1', u1, use_reverse, sys_info);
        x2n = rk4(dynamics, dt_sim, x2', u2, use_reverse, sys_info);
        x3n = rk4(dynamics, dt_sim, x3', u3, use_reverse, sys_info);

        % update states
        x1 = x1n'; x2 = x2n'; x3 = x3n';
        X1 = [X1; x1];
        X2 = [X2; x2]; 
        X3 = [X3; x3];

        % update controls (row-wise)
        U1 = [U1; u1']; 
        U2 = [U2; u2']; 
        U3 = [U3; u3']; 

        iter = iter + 1;
    end

    % store time traces
    X1_all{k} = X1;  U1_all{k} = U1;
    X2_all{k} = X2;  U2_all{k} = U2;
    X3_all{k} = X3;  U3_all{k} = U3;

    % plot phase-portrait trajectories with alpha (RGBA)
    % plot(axp, X1(:,1), X1(:,2), '-', 'Color', [C(3,:) alpha_line], 'HandleVisibility','off'); % LQR
    % plot(axp, X3(:,1), X3(:,2), '-', 'Color', [C(2,:) alpha_line], 'HandleVisibility','off'); % LQR-FT
    plot(axp, X2(:,1), X2(:,2), '-', 'Color', [C(1,:) alpha_line], 'HandleVisibility','off');   % SKOOPI
end

% Re-plot disk & center
tt = linspace(0,2*pi,200);
plot(axp, center(1)+radius*cos(tt), center(2)+radius*sin(tt), ':k');
plot(axp, center(1), center(2), 'ok', 'MarkerFaceColor','k');

% add legend (opaque reps for clarity)
plot(axp, NaN, NaN, '-', 'Color', C(1,:), 'LineWidth', 2.2, 'DisplayName','LQR');
plot(axp, NaN, NaN, '-', 'Color', C(2,:), 'LineWidth', 2.2, 'DisplayName','LQR-FT');
plot(axp, NaN, NaN, '-', 'Color', C(3,:), 'LineWidth', 2.2, 'DisplayName','SKOOPI');
xlim([-20,20]); ylim([-20,20]);

%% ---- New figures: state & control trajectories vs time for each method ----
t = t_grid(:);               % Nsteps x 1
tc = t(2:end);               % control timeline (Nsteps-1 x 1)

plot_method_time_and_controls(X1_all, U1_all, C(1,:), 'LQR',    alpha_line, t, tc, n_ctrl);
plot_method_time_and_controls(X3_all, U3_all, C(2,:), 'LQR-FT', alpha_line, t, tc, n_ctrl);
plot_method_time_and_controls(X2_all, U2_all, C(1,:), 'SKOOPI', alpha_line, t, tc, n_ctrl);

%% ---- RMSE and NRMSE (overall) ----
[rmse1, nrmse1] = overall_rmse_nrmse_states(X1_all);
[rmse2, nrmse2] = overall_rmse_nrmse_states(X2_all);
[rmse3, nrmse3] = overall_rmse_nrmse_states(X3_all);

fprintf('\n===== Overall tracking error (states vs zero) =====\n');
fprintf('LQR    : RMSE = %.4f, NRMSE = %.2f%%%%\n', rmse1, nrmse1);
fprintf('SKOOPI : RMSE = %.4f, NRMSE = %.2f%%%%\n', rmse2, nrmse2);
fprintf('LQR-FT : RMSE = %.4f, NRMSE = %.2f%%%%\n\n', rmse3, nrmse3);

%% ==================== local helpers ====================
function plot_method_time_and_controls(Xcell, Ucell, color_rgb, name_str, a, t, tc, n_ctrl)
    % Creates one figure per controller with states (top) and controls (bottom).
    figure(1)
    figure('Name',[name_str ' — states & controls vs time'],'Color','w');
    sgtitle([name_str ' — trajectories from ' num2str(numel(Xcell)) ' ICs'], 'Interpreter','latex');

    % ---- x1(t) in subplot(4,4,[1 2]) ----
    subplot(4,4,[1 2]); hold on; grid on; box on;
    for i = 1:numel(Xcell)
        plot(t, Xcell{i}(:,1), 'Color', [color_rgb a], 'HandleVisibility','off');
    end
    plot(t, Xcell{1}(:,1), 'Color', color_rgb, 'LineWidth', 2, 'DisplayName','$x_1$ (FT-LQR)');
    ylabel('$x_1(t)$','Interpreter','latex');
    legend('Interpreter','latex','Location','best');
    ylim auto; xlim([t(1) t(end)]);

    % ---- x2(t) in subplot(4,4,[5 6]) ----
    subplot(4,4,[5 6]); hold on; grid on; box on;
    for i = 1:numel(Xcell)
        plot(t, Xcell{i}(:,2), 'Color', [color_rgb a], 'HandleVisibility','off');
    end
    plot(t, Xcell{1}(:,2), 'Color', color_rgb, 'LineWidth', 2, 'DisplayName','$x_2$ (FT-LQR)');
    xlabel('time (s)','Interpreter','latex');
    ylabel('$x_2(t)$','Interpreter','latex');
    legend('Interpreter','latex','Location','best');
    ylim auto; xlim([t(1) t(end)]);

    % ---- controls: place u_1 and u_2 (if present) on lower rows ----
    % u_1 in subplot(4,4,[9 10])
    subplot(4,4,[9 10]); hold on; grid on; box on;
    for i = 1:numel(Ucell)
        % Each U is (Nsteps-1) x n_ctrl
        plot(tc, Ucell{i}(:,1), 'Color', [color_rgb a], 'HandleVisibility','off');
    end
    plot(tc, Ucell{1}(:,1), 'Color', color_rgb, 'LineWidth', 2, 'DisplayName','$u$ (FT-LQR)');
    xlabel('time (s)','Interpreter','latex');
    ylabel('$u_1(t)$','Interpreter','latex');
    legend('Interpreter','latex','Location','best');
    xlim([tc(1) tc(end)]);
end

function [rmse_overall, nrmse_pct] = overall_rmse_nrmse_states(Xcell)
% Compute overall RMSE and NRMSE (percent) across all states and all ICs
% RMSE: sqrt(mean(x^2)) w.r.t. zero reference, pooling all time/IC samples
% NRMSE: average over states of RMSE_state / range_state * 100
    % Concatenate all trajectories vertically
    Xcat = [];
    for i = 1:numel(Xcell)
        Xcat = [Xcat; Xcell{i}]; %#ok<AGROW>
    end
    % Overall RMSE (all states pooled)
    rmse_overall = sqrt(mean(Xcat(:).^2));

    % Per-state NRMSE, then average
    n_states_local = size(Xcat,2);
    nrmse_states = zeros(n_states_local,1);
    for d = 1:n_states_local
        xd = Xcat(:,d);
        rmse_d = sqrt(mean(xd.^2));
        range_d = max(xd) - min(xd);
        if range_d <= eps
            nrmse_states(d) = 0; % degenerate (already at/near zero)
        else
            nrmse_states(d) = 100 * rmse_d / range_d;
        end
    end
    nrmse_pct = mean(nrmse_states);
end
