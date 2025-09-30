clc; clear; close all

%% ---------- Default plot options ----------
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultAxesFontSize',18)
set(0,'defaultfigurecolor',[1 1 1])

%% ---------- Add paths (expects your existing folders/functions) ----------
addpath('dynamics', 'baseline_control', 'eigfun_control', ...
        'compute_eigfuns', 'utils', 'animations');
rng(0)

%% ---------- System setup ----------
show_diagnositcs = false;     
use_reverse       = false;

sys_params.use_stable         = true;  
sys_params.use_unstable       = false;

% Load dynamics & system info
sys_info  = nonlinear_sys_info(sys_params);
dynamics  = @dynamics_nonlinear;
n_states  = sys_info.state_dim;
n_ctrl    = sys_info.ctrl_dim;

A = sys_info.A;
B = sys_info.B;

% Analytical eigenvector matrix at origin (your helper)
W = sys_info.grad_eigen_fun_analytical([0;0]);
sys_info.eig_vectors = W;

%% ---------- Cost & solver params ----------
Q  = diag([1 1]);
R  = 1e2*eye(n_ctrl);
QN = 1e1*Q;

% (transformed for SKOOPI uses A, g_list1 below)

%% ---------- Time setup ----------
dt_sim   = 0.01;
t_start  = 0;
t_end    = 5.0;

%% ---------- Precompute reference & time-varying gains (tracking form) ----------
% Build an open-loop reference (same style as your script)
x_desired = [2;3];
psi_list  = {};               % for SKOOPI (psi-coords)
x_desired_list = {};
u_desired_list = {};
g_list1 = {}; g_list2 = {};   % g along ref for SKOOPI vs baseline (B)

for t_sim = dt_sim:dt_sim:t_end
    % open-loop input (feel free to change)
    u_o = 5*cos(0.5*t_sim) - 5;

    % g(x) along reference
    g_t = sys_info.dynamics_g(x_desired);

    % store psi_d for SKOOPI
    psi_list{end+1}       = sys_info.transform_fun_analytical(x_desired); 
    x_desired_list{end+1} = x_desired;                                     
    u_desired_list{end+1} = u_o;                                         
    g_list1{end+1}        = g_t;                                       
    g_list2{end+1}        = B;                                            

    % propagate ref
    x_desired = rk4(dynamics, dt_sim, x_desired, u_o, use_reverse, sys_info);
end

x_ref = cell2mat(x_desired_list);
u_ref = cell2mat(u_desired_list);

% Transformed-coords (SKOOPI) tracking gains
[~, ~, K_ff_list1, K_fb_list1, ~] = ...
    solve_DRE_tracking(x_ref, u_ref, A, g_list1, Q, R, QN, t_end, dt_sim);

% Baseline (original coords) tracking gains (LQR-FT)
[~, ~, K_ff_list2, K_fb_list2, ~] = ...
    solve_DRE_tracking(x_ref, u_ref, A, g_list2, Q, R, QN, t_end, dt_sim);

%% ---------- Single-IC quick run (optional / for context) ----------
start_x = [2;3];
x_op1   = start_x';   % LQR-FT
x_op2   = start_x';   % SKOOPI
Tout    = t_start;
Xout1   = x_op1;
Xout2   = x_op2;

idx = 1;
for t_sim = dt_sim:dt_sim:t_end-dt_sim
    % LQR-FT control
    Kfb1 = K_fb_list2{idx};  Kff1 = K_ff_list2{idx};
    xd   = x_ref(:, idx);
    u1   = -Kfb1 * (x_op1' - xd) + Kff1;
    
    % SKOOPI control (psi)
    Kfb2 = K_fb_list1{idx};  Kff2 = K_ff_list1{idx};
    psi_x = sys_info.transform_fun_analytical(x_op2');
    psi_d = psi_list{idx};
    u2    = -Kfb2 * (psi_x - psi_d) + Kff2;

    % clamp
    if ~isfinite(u1) || abs(u1)>1e3, u1=0; end
    if ~isfinite(u2) || abs(u2)>1e3, u2=0; end

    % step
    x1n = rk4(dynamics, dt_sim, x_op1', u1, use_reverse, sys_info);
    x2n = rk4(dynamics, dt_sim, x_op2', u2, use_reverse, sys_info);

    % logs
    Tout  = [Tout; t_sim];
    Xout1 = [Xout1; x1n'];
    Xout2 = [Xout2; x2n'];
    x_op1 = x1n'; x_op2 = x2n'; idx = idx + 1;
end

%% ---------- Phase portrait base (figure 1) ----------
if ~ishandle(1)
    figure(1);  % create if it doesn't exist
end
figure(1);  % make fig.1 current

axp = subplot(4,4,[1 2 5 6 9 10]);  % <- adjust if your phase panel is different
hold(axp,'on'); box(axp,'on'); grid(axp,'on');
axp = subplot(4,4,[1 2 5 6 9 10]); hold(axp,'on'); box(axp,'on'); grid(axp,'on')
xlabel(axp,'state, $x_1$','Interpreter','latex');
ylabel(axp,'state, $x_2$','Interpreter','latex');
axis(axp,'equal')

C = axp.ColorOrder;
colLQR   = C(1,:);
colSKOOP = C(2,:);
alphaLine = 0.18;

% plot single-IC trajectories (hidden from legends)
plot(axp, Xout1(:,1), Xout1(:,2), '-', 'Color', [colLQR 0.7], 'HandleVisibility','off');
plot(axp, Xout2(:,1), Xout2(:,2), '-', 'Color', [colSKOOP 0.7], 'HandleVisibility','off');
plot(axp, x_ref(1,:), x_ref(2,:), '--k', 'HandleVisibility','off');

% IC disk (will reuse for sweep as well)
center = [2; 3];
radius = 0.50;
tt = linspace(0,2*pi,200);
plot(axp, center(1)+radius*cos(tt), center(2)+radius*sin(tt), ':k', 'HandleVisibility','off');
plot(axp, center(1), center(2), 'ok', 'MarkerFaceColor','k', 'HandleVisibility','off');
xlim([-15,5])
ylim([-3,8])

%% ---------- IC SWEEP from a DISK (100 ICs) ----------
N_ic = 100;
theta = 2*pi*rand(N_ic,1);
rad   = radius * sqrt(rand(N_ic,1));
ICs   = center.' + [rad.*cos(theta), rad.*sin(theta)];  % N_ic x 2

% time vector for sweep (drop t=0 to align steps)
t_vec = t_start + dt_sim : dt_sim : t_end;
Tn = numel(t_vec);

% prealloc for error stats
E1_LQR   = nan(Tn, N_ic);  E1_SKOOP = nan(Tn, N_ic);
E2_LQR   = nan(Tn, N_ic);  E2_SKOOP = nan(Tn, N_ic);

% store time series if needed later
Xb_all = cell(N_ic,1);
Xb_all = cell(N_ic,1);
Ub_all = cell(N_ic,1);   % controls for LQR-FT (each cell: Tn x 1)
Us_all = cell(N_ic,1);   % controls for SKOOPI (each cell: Tn x 1)

for k = 1:N_ic
    x1 = ICs(k,:);    % LQR-FT
    x2 = x1;          % SKOOPI
    X_b = x1; X_s = x2;
    U_b = [];  
    U_s = [];

    idx_s = 1;
    for t_sim = dt_sim:dt_sim:t_end
        % LQR-FT
        if idx_s <= numel(K_fb_list2)
            Kfb1 = K_fb_list2{idx_s};  Kff1 = K_ff_list2{idx_s};
        else
            Kfb1 = K_fb_list2{end};    Kff1 = K_ff_list2{end};
        end
        xd = x_ref(:, min(idx_s, size(x_ref,2)));
        u1 = -Kfb1 * (x1' - xd) + Kff1;

        % SKOOPI (psi)
        if idx_s <= numel(K_fb_list1)
            Kfb2 = K_fb_list1{idx_s};  Kff2 = K_ff_list1{idx_s}; psi_d = psi_list{idx_s};
        else
            Kfb2 = K_fb_list1{end};    Kff2 = K_ff_list1{end};   psi_d = psi_list{end};
        end
        psi_x = sys_info.transform_fun_analytical(x2');
        u2    = -Kfb2 * (psi_x - psi_d) + Kff2;

        % clamp
        if ~isfinite(u1) || abs(u1)>1e3, u1=0; end
        if ~isfinite(u2) || abs(u2)>1e3, u2=0; end

        % step
        x1n = rk4(dynamics, dt_sim, x1', u1, use_reverse, sys_info);
        x2n = rk4(dynamics, dt_sim, x2', u2, use_reverse, sys_info);

        % log
        X_b = [X_b; x1n'];  X_s = [X_s; x2n'];
        x1 = x1n'; x2 = x2n'; idx_s = idx_s + 1;
        U_b = [U_b; u1];   
        U_s = [U_s; u2];  
    end

    Xb_all{k} = X_b;  Xb_all{k} = X_s;
    Ub_all{k} = U_b;      % Tn x 1
    Us_all{k} = U_s;      % Tn x 1

    % phase portrait overlay (transparent, no legend)
    plotTransparentLine(axp, X_b(:,1), X_b(:,2), colLQR,   alphaLine);
    plotTransparentLine(axp, X_s(:,1), X_s(:,2), colSKOOP, alphaLine);

    % tracking errors aligned with t_vec (drop initial row)
    Xb_use = X_b(2:end,:);    % Tn x 2
    Xs_use = X_s(2:end,:);
    xref_use = x_ref(:, 1:Tn).';
    E1_LQR(:,k)   = Xb_use(:,1) - xref_use(:,1);
    E1_SKOOP(:,k) = Xs_use(:,1) - xref_use(:,1);
    E2_LQR(:,k)   = Xb_use(:,2) - xref_use(:,2);
    E2_SKOOP(:,k) = Xs_use(:,2) - xref_use(:,2);
end

% ---------- Legend (phase portrait) — proxies only ----------
[hLQRp, hSKOOPp, hREFp] = legendProxies(axp, colLQR, colSKOOP);
legend(axp, [hLQRp hSKOOPp hREFp], {'LQR-FT','SKOOPI','reference'}, ...
       'Interpreter','latex','Location','northoutside','NumColumns',3, ...
       'AutoUpdate','off');

% IC disk (will reuse for sweep as well)
center = [2; 3];
radius = 0.50;
tt = linspace(0,2*pi,200);
plot(axp, center(1)+radius*cos(tt), center(2)+radius*sin(tt), ':k', 'HandleVisibility','off');
plot(axp, center(1), center(2), 'ok', 'MarkerFaceColor','k', 'HandleVisibility','off');

%% ---------- Time-series & errors (figure 2, 4×4 layout) ----------
figure(2); clf
ax_ts1 = subplot(4,4,[1 2]); hold(ax_ts1,'on'); box(ax_ts1,'on'); grid(ax_ts1,'on')
ax_ts2 = subplot(4,4,[5 6]); hold(ax_ts2,'on'); box(ax_ts2,'on'); grid(ax_ts2,'on')
ax_e1  = subplot(4,4,[3 4]); hold(ax_e1,'on');  box(ax_e1,'on');  grid(ax_e1,'on')
ax_e2  = subplot(4,4,[7 8]); hold(ax_e2,'on');  box(ax_e2,'on');  grid(ax_e2,'on')

% spaghetti for states (hidden from legends)
for k = 1:N_ic
    Xb_use = Xb_all{k}(2:end,:);   % Tn x 2
    Xs_use = Xb_all{k}(2:end,:);
    plotTransparentLine(ax_ts1, t_vec(:), Xb_use(:,1), colLQR,   alphaLine);
    plotTransparentLine(ax_ts1, t_vec(:), Xs_use(:,1), colSKOOP, alphaLine);
    plotTransparentLine(ax_ts2, t_vec(:), Xb_use(:,2), colLQR,   alphaLine);
    plotTransparentLine(ax_ts2, t_vec(:), Xs_use(:,2), colSKOOP, alphaLine);
end

% reference overlays (hidden from legends; we’ll use proxies)
plot(ax_ts1, t_vec, x_ref(1,1:Tn), '--k', 'HandleVisibility','off');
plot(ax_ts2, t_vec, x_ref(2,1:Tn), '--k', 'HandleVisibility','off');

ylabel(ax_ts1,'$x_1$','Interpreter','latex'); xlim(ax_ts1,[t_vec(1) t_vec(end)]);
ylabel(ax_ts2,'$x_2$','Interpreter','latex'); xlim(ax_ts2,[t_vec(1) t_vec(end)]);
xlabel(ax_ts2,'time (s)','Interpreter','latex');

% error spaghetti (hidden from legends)
for k = 1:N_ic
    plotTransparentLine(ax_e1, t_vec(:), E1_LQR(:,k),   colLQR,   alphaLine);
    plotTransparentLine(ax_e1, t_vec(:), E1_SKOOP(:,k), colSKOOP, alphaLine);
    plotTransparentLine(ax_e2, t_vec(:), E2_LQR(:,k),   colLQR,   alphaLine);
    plotTransparentLine(ax_e2, t_vec(:), E2_SKOOP(:,k), colSKOOP, alphaLine);
end

% mean ± std bands (also hidden from legends)
mu_e1_lqr   = mean(E1_LQR,   2);  sd_e1_lqr   = std(E1_LQR,   0, 2);
mu_e1_skoop = mean(E1_SKOOP, 2);  sd_e1_skoop = std(E1_SKOOP, 0, 2);
mu_e2_lqr   = mean(E2_LQR,   2);  sd_e2_lqr   = std(E2_LQR,   0, 2);
mu_e2_skoop = mean(E2_SKOOP, 2);  sd_e2_skoop = std(E2_SKOOP, 0, 2);

fillBetween(ax_e1, t_vec, mu_e1_lqr - sd_e1_lqr,   mu_e1_lqr + sd_e1_lqr,   colLQR,   0.10);
fillBetween(ax_e1, t_vec, mu_e1_skoop - sd_e1_skoop, mu_e1_skoop + sd_e1_skoop, colSKOOP, 0.10);
plot(ax_e1, t_vec, mu_e1_lqr,   '-','Color',colLQR,   'LineWidth',2,'HandleVisibility','off');
plot(ax_e1, t_vec, mu_e1_skoop, '-','Color',colSKOOP, 'LineWidth',2,'HandleVisibility','off');

fillBetween(ax_e2, t_vec, mu_e2_lqr - sd_e2_lqr,   mu_e2_lqr + sd_e2_lqr,   colLQR,   0.10);
fillBetween(ax_e2, t_vec, mu_e2_skoop - sd_e2_skoop, mu_e2_skoop + sd_e2_skoop, colSKOOP, 0.10);
plot(ax_e2, t_vec, mu_e2_lqr,   '-','Color',colLQR,   'LineWidth',2,'HandleVisibility','off');
plot(ax_e2, t_vec, mu_e2_skoop, '-','Color',colSKOOP, 'LineWidth',2,'HandleVisibility','off');

% zero error lines (hidden)
plot(ax_e1, t_vec, zeros(size(t_vec)), ':k', 'HandleVisibility','off');
plot(ax_e2, t_vec, zeros(size(t_vec)), ':k', 'HandleVisibility','off');

ylabel(ax_e1,'$e_1$','Interpreter','latex'); xlim(ax_e1,[t_vec(1) t_vec(end)]);
ylabel(ax_e2,'$e_2$','Interpreter','latex'); xlim(ax_e2,[t_vec(1) t_vec(end)]);
xlabel(ax_e2,'time (s)','Interpreter','latex');

% ---------- Legend (time-series top panel) — proxies only ----------
[hLQRt, hSKOOPt, hREFt] = legendProxies(ax_ts1, colLQR, colSKOOP);
legend(ax_ts1, [hLQRt hSKOOPt hREFt], {'LQR-FT','SKOOPI','reference'}, ...
       'Interpreter','latex','Location','northoutside','NumColumns',3, ...
       'AutoUpdate','off');

%% ---------- Controls: per-IC spaghetti + mean ± std (single axes) ----------
figure(3); clf
ax_u = subplot(4,4,[1 2 5 6]); hold(ax_u,'on'); box(ax_u,'on'); grid(ax_u,'on')
xlabel(ax_u,'time (s)','Interpreter','latex')
ylabel(ax_u,'$u$','Interpreter','latex')
xlim(ax_u, [t_vec(1) t_vec(end)])

% Ensure control matrices are Tn x N_ic
if exist('U1_mat','var')~=1 || isempty(U1_mat)
    % Build from cells: horizontal concat yields Tn x N_ic
    U1_mat = cell2mat(Ub_all.');   % LQR-FT
    U2_mat = cell2mat(Us_all.');   % SKOOPI
end

% 1) per-IC spaghetti with transparency (hidden from legend)
for k = 1:N_ic
    plotTransparentLine(ax_u, t_vec(:), U1_mat(:,k), colLQR,   alphaLine);  % LQR-FT
    plotTransparentLine(ax_u, t_vec(:), U2_mat(:,k), colSKOOP, alphaLine);  % SKOOPI
end

% 2) mean ± std overlays (hidden fills; solid means)
mu_u1  = mean(U1_mat, 2);   sd_u1 = std(U1_mat, 0, 2);
mu_u2  = mean(U2_mat, 2);   sd_u2 = std(U2_mat, 0, 2);

fillBetween(ax_u, t_vec, mu_u1 - sd_u1,   mu_u1 + sd_u1,   colLQR,   0.12); % hidden from legend
fillBetween(ax_u, t_vec, mu_u2 - sd_u2,   mu_u2 + sd_u2,   colSKOOP, 0.12); % hidden from legend
plot(ax_u, t_vec, mu_u1, '-', 'Color', colLQR,   'LineWidth', 2, 'HandleVisibility','off');
plot(ax_u, t_vec, mu_u2, '-', 'Color', colSKOOP, 'LineWidth', 2, 'HandleVisibility','off');

% 3) single clean legend via proxies (only two entries here)
[hLQRu, hSKOOPu, ~] = legendProxies(ax_u, colLQR, colSKOOP);
legend(ax_u, [hLQRu hSKOOPu], {'LQR-FT','SKOOPI'}, ...
       'Interpreter','latex','Location','northoutside','NumColumns',2, ...
       'AutoUpdate','off');

%% ---------- Overall RMSE (across time, ICs, and both states) ----------
% Per-state RMSEs (averaged over all time samples and ICs)
rmse_x1_lqr   = sqrt(mean(E1_LQR(:).^2));
rmse_x2_lqr   = sqrt(mean(E2_LQR(:).^2));
rmse_x1_skoop = sqrt(mean(E1_SKOOP(:).^2));
rmse_x2_skoop = sqrt(mean(E2_SKOOP(:).^2));

% Overall RMSE combining x1 and x2
rmse_lqr   = sqrt(mean([E1_LQR(:).^2;   E2_LQR(:).^2]));
rmse_skoop = sqrt(mean([E1_SKOOP(:).^2; E2_SKOOP(:).^2]));

% Print nicely
fprintf('\nRMSE (per state):\n');
fprintf('  LQR-FT   : x1 = %.4g, x2 = %.4g\n', rmse_x1_lqr,   rmse_x2_lqr);
fprintf('  SKOOPI   : x1 = %.4g, x2 = %.4g\n', rmse_x1_skoop, rmse_x2_skoop);

%%
% Assuming x_ref_mat: Nsteps x N_ic matrix of reference trajectory
x1_reff = x_ref(1,:);  % adjust indexing for your setup
x2_reff = x_ref(2,:);

% RMSE already computed: rmse_x1_lqr, rmse_x2_lqr, etc.

nrmse_x1_lqr   = rmse_x1_lqr   / rms(x1_reff);
nrmse_x1_skoop = rmse_x1_skoop / rms(x1_reff);
nrmse_x2_lqr   = rmse_x2_lqr   / rms(x2_reff);
nrmse_x2_skoop = rmse_x2_skoop / rms(x2_reff);

fprintf('NRMSE (x1): LQR-FT = %.3f, SKOOPI = %.3f\n', nrmse_x1_lqr, nrmse_x1_skoop);
fprintf('NRMSE (x2): LQR-FT = %.3f, SKOOPI = %.3f\n', nrmse_x2_lqr, nrmse_x2_skoop);

%% ---------- Overall NRMSE & Improvement ----------
% Reference RMS magnitude (combine x1 and x2)
xref_cat = [x_ref(1,:).'; x_ref(2,:).'];        % stack both states
ref_rms  = rms(xref_cat(:));                    % overall RMS of reference

% Raw RMSE (already combined across states & ICs)
rmse_lqr   = sqrt(mean([E1_LQR(:).^2;   E2_LQR(:).^2]));
rmse_skoop = sqrt(mean([E1_SKOOP(:).^2; E2_SKOOP(:).^2]));

% Normalized RMSE (dimensionless)
nrmse_lqr   = rmse_lqr   / ref_rms;
nrmse_skoop = rmse_skoop / ref_rms;

% Relative improvement (positive means SKOOPI is better)
improvement_pct = 100 * (1 - nrmse_skoop / nrmse_lqr);

% Display
fprintf('\nOverall RMSE   : LQR-FT = %.4f,  SKOOPI = %.4f\n', rmse_lqr, rmse_skoop);
fprintf('Overall NRMSE  : LQR-FT = %.4f,  SKOOPI = %.4f\n', nrmse_lqr, nrmse_skoop);
fprintf('Relative improvement (SKOOPI vs LQR-FT): %.2f%%\n', improvement_pct);


%% ==================== Helpers (keep at end of file) ====================

function plotTransparentLine(ax, x, y, rgb, a)
% Version-robust transparent polyline via PATCH, hidden from legends.
    x = x(:); y = y(:);
    if numel(x) ~= numel(y)
        error('plotTransparentLine: x and y must be same length.');
    end
    xx = [x; NaN]; yy = [y; NaN];
    patch('Parent',ax,'XData',xx,'YData',yy, ...
          'FaceColor','none','EdgeColor',rgb,'EdgeAlpha',a, ...
          'LineWidth',1.5,'LineJoin','round','LineStyle','-', ...
          'HandleVisibility','off');   % <<< hide from legend
end

function fillBetween(ax, t, lo, hi, rgb, a)
% Filled band (mean±std etc.), hidden from legends.
    t = t(:); lo = lo(:); hi = hi(:);
    [t, i] = sort(t); lo = lo(i); hi = hi(i);
    X = [t; flipud(t)];
    Y = [lo; flipud(hi)];
    patch('Parent',ax, 'XData',X, 'YData',Y, ...
          'FaceColor',rgb, 'FaceAlpha',a, 'EdgeColor','none', ...
          'HandleVisibility','off');  
end

function [hLQR, hSKOOP, hREF] = legendProxies(ax, colLQR, colSKOOP)
    hLQR   = plot(ax, NaN,NaN,'-','Color',colLQR,   'LineWidth',2,   'HandleVisibility','on');
    hSKOOP = plot(ax, NaN,NaN,'-','Color',colSKOOP, 'LineWidth',2,   'HandleVisibility','on');
    hREF   = plot(ax, NaN,NaN,'--k','LineWidth',1.5,'HandleVisibility','on');
end
