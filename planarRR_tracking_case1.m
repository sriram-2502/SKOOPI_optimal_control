%% SKOOPI for robotic arm — one-shot (no chunking), circular ref + ID feedforward
% Full script with 10 random initial-condition sweeps (joint/task space)
% Save as: skoopi_rr_multiIC.m

clc; clear; close all
rng(134);

% ---------- default plot options ----------
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultAxesFontSize',18)
set(0,'defaultfigurecolor',[1 1 1])

% ----- paths (adapt as needed) -----
addpath(genpath('robotic_arm'))
addpath('dynamics')
addpath('baseline_control')
addpath('eigfun_control')
addpath('compute_eigfuns')
addpath('utils')
addpath('animations')
addpath('planarRR')

% ----- flags & system toggles ----- 
show_diagnositcs = true;
sys_params.use_stable = false;
sys_params.use_unstable = false;

%% Problem setup (symbolics placeholders as your utils expect)
x = sym('x',[2 1],'real'); 
q = sym('q',[4 1],'real'); 
u = sym('u',[2 1],'real'); 
syms t;

% ----- get params / models -----
[robot_params,euler_params,navigation_params,animate_params] = get_params_planarRR(x);

% parse system information (provides sym M,C,G; numeric/handles A,B,g, etc.)
sys_info = planarRR_info(0, q, u, robot_params, sys_params);
dynamics = @dynamics_planarRR;

% unpack (M,C,G are symbolic in your setup)
M_sym = sys_info.M; 
C_sym = sys_info.C;  
G_sym = sys_info.G;  
A      = sys_info.A; 
B      = sys_info.B;
g_fun  = sys_info.dynamics_g;
W      = sys_info.eig_vectors;
x_eqb  = [0 0 0 0]; 

% ----- make W real (pair real/imag parts) -----
% assumes W has complex-conjugate structure in columns 1..4
W = [real(W(:,1)), imag(W(:,2)), real(W(:,3)), imag(W(:,4))];

% ===== Controller weights =====
Q  = 1e2*diag([1 1 0.01 0.01]);
R  = 1e3*diag([1 1]);
QN = 1e1*Q;

%% ===== Task-space circular reference (downward workspace) =====
xc = 0.0;  yc = 0.0;     % circle center (m)
r  = 1.5;                % radius (m)
omega = 0.2;             % rad/s
% phase-shifted so it starts at the bottom: (t -> t - pi/(2*omega))
time_varying_goal = [xc + r*cos(omega*t-pi/2), yc + r*sin(omega*t-pi/2)];
matlabFunction(time_varying_goal, 'File', 'planarRR/x_ref_f', 'Vars', {t}, 'Optimize', false);

% IK and time-derivatives symbols → functions
navigation_params.x_goal_f       = IK_planarRR(time_varying_goal);
navigation_params.x_vel_goal_f   = diff(navigation_params.x_goal_f);
navigation_params.x_accel_goal_f = diff(navigation_params.x_vel_goal_f);

matlabFunction(navigation_params.x_goal_f,       'File', 'planarRR/xd_f',       'Vars', {t}, 'Optimize', false);
matlabFunction(navigation_params.x_vel_goal_f,   'File', 'planarRR/x_dot_d_f',  'Vars', {t}, 'Optimize', false);
matlabFunction(navigation_params.x_accel_goal_f, 'File', 'planarRR/x_ddot_d_f', 'Vars', {t}, 'Optimize', false);
matlabFunction(jacobian(navigation_params.x_goal_f, t), 'File', 'planarRR/jac_xd_f', 'Vars', {t}, 'Optimize', false);
rehash path

% Convenience handle for EE reference (2x1 column)
x_ref_f = @x_ref_f;

%% ===== Discretize reference & build joint-space (q, qd, qdd) =====
t_vector = (0:euler_params.n_steps).*euler_params.step_size;
dt       = euler_params.step_size;
Tgrid    = numel(t_vector);
dt_sim   = euler_params.step_size;
t_end    = (Tgrid-1)*dt_sim;

% Evaluate IK along the circle
q_goal = zeros(Tgrid, 2);
for i = 1:Tgrid
    q_goal(i,:) = IK_planarRR(x_ref_f(t_vector(i)));
end

% Finite-difference velocities and accelerations (same length as q_goal)
q_dot_goal  = zeros(size(q_goal));
q_ddot_goal = zeros(size(q_goal));

% velocities (central where possible)
q_dot_goal(2:end-1,:) = (q_goal(3:end,:) - q_goal(1:end-2,:)) / (2*dt);
q_dot_goal(1,:)       = (q_goal(2,:)     - q_goal(1,:))       / dt;
q_dot_goal(end,:)     = (q_goal(end,:)   - q_goal(end-1,:))   / dt;

% accelerations
q_ddot_goal(2:end-1,:) = (q_goal(3:end,:) - 2*q_goal(2:end-1,:) + q_goal(1:end-2,:)) / (dt^2);
q_ddot_goal(1,:)        = (q_dot_goal(2,:)   - q_dot_goal(1,:)) / dt;
q_ddot_goal(end,:)      = (q_dot_goal(end,:) - q_dot_goal(end-1,:)) / dt;

% Store for potential downstream usage
navigation_params.q_goal      = q_goal;      
navigation_params.q_dot_goal  = q_dot_goal;   
navigation_params.q_ddot_goal = q_ddot_goal; 

if show_diagnositcs
    figure('Name','Joint references','Color','w');
    subplot(3,1,1); plot(t_vector, q_goal,'-o');      grid on; title('q_d')
    subplot(3,1,2); plot(t_vector, q_dot_goal,'-o');  grid on; title('\dot q_d')
    subplot(3,1,3); plot(t_vector, q_ddot_goal,'-o'); grid on; title('\ddot q_d')
end

%% ===== Feedforward torques u_ref(t) and eigenfunctions along reference =====
u_ref_series   = zeros(2, Tgrid);
psi_ref_series = zeros(4, Tgrid);

w_bar = waitbar(0, '1', 'Name', sprintf('getting reference'));
for i = 1:Tgrid
    waitbar(i / Tgrid, w_bar, sprintf('%d / %d', i, Tgrid));

    qd   = q_goal(i,:).';  
    dqd  = q_dot_goal(i,:).';
    ddq  = q_ddot_goal(i,:).';
    xref = [q_goal(i,:).'; q_dot_goal(i,:).' ];

    % compute phi and psi along ref traj
    phi = compute_path_integrals(xref, dynamics, sys_info);
    phi_x_op = phi.phi_x_op';
    phi_x_op = [ real(phi_x_op(1)); imag(phi_x_op(1)); ...
                 real(phi_x_op(3)); imag(phi_x_op(3)) ];
    psi_ref_series(:,i) = xref + inv(W') * phi_x_op;

    % feedforward torques
    params = planarRR_info(t_vector(i), xref, [0;0], robot_params, sys_params);
    M_t = params.M; C_t = params.C; G_t = params.G; D_t = params.D;
    u_ref_series(:,i) = M_t*ddq + C_t + G_t + D_t*dqd;
end
delete_valid_waitbars();

%% ===== Assemble reference state series & g-lists (one shot) =====
x_ref_series = zeros(4, Tgrid);
g_list_skoopi = cell(1, Tgrid);
g_list_base   = cell(1, Tgrid);

for i = 1:Tgrid
    x_ref_series(:,i) = [ q_goal(i,:).'; q_dot_goal(i,:).' ];
    g_list_skoopi{i}  = g_fun( x_ref_series(:,i) );
    g_list_base{i}    = B;    % baseline uses constant B
end

x_ref = x_ref_series;
u_ref = u_ref_series;

%% ===== Solve tracking DRE once (no chunking) =====
[L_b, P_b, K_ff_list1, K_fb_list1, sx_b] = solve_DRE_tracking( ...
    x_ref, u_ref, A, g_list_base,   Q, R, QN, t_end, dt_sim); 

[L_s, P_s, K_ff_list2, K_fb_list2, sx_s] = solve_DRE_tracking( ...
    x_ref, u_ref, A, g_list_skoopi, Q, R, QN, t_end, dt_sim); 

% diagnostics: show u_ref
if show_diagnositcs
    figure('Name','u_{ref}','Color','w');
    plot(t_vector, u_ref_series(1,:),'-', t_vector, u_ref_series(2,:),'--'); grid on
    xlabel('t [s]'); ylabel('\tau [N·m]'); legend('\tau_1^{ref}','\tau_2^{ref}')
end

% animate reference (optional)
animate_planarRR(t_vector,q_goal,[],[xc; yc]);

delete_valid_waitbars();

%% ======================= multi-IC simulation loop =======================

% ----- IC mode & perturbation sizes -----
ic_mode = "joint";            % "joint" or "task"
N_ic    = 10;                 % number of initial conditions
dq0_max = deg2rad(3);         % max |q0 offset| per joint
dqdot0_max = 0;%deg2rad(10);  % max |qdot0| per joint
task_radius = 0.5;            % meters (used if ic_mode=="task")

% convenience
dt       = euler_params.step_size;
T        = numel(t_vector);
q0_ref   = q_goal(1,:).';                     % starting joint reference
xEE0_ref = x_ref_f(t_vector(1)).';            % starting EE reference [x;y]

% precompute inv(W') once
WinvT = inv(W.');

% storage
X_skoopi_all = cell(1,N_ic);   U_skoopi_all = cell(1,N_ic);
X_base_all   = cell(1,N_ic);   U_base_all   = cell(1,N_ic);

w_bar = waitbar(0, '1', 'Name', sprintf('running %d ICs', N_ic));
for ii = 1:N_ic
    waitbar(ii / N_ic, w_bar, sprintf('IC %d / %d', ii, N_ic));

    % ----- initial state (q0, qdot0) near ref start -----
    switch ic_mode
        case "joint"
            q0 = q0_ref + [0;(rand(1,1)-0.5)*2*dq0_max];
        case "task"
            ang = 2*pi*rand;   rad = task_radius*sqrt(rand);
            xEE = xEE0_ref + [rad*cos(ang); rad*sin(ang)];
            q0  = IK_planarRR(xEE.');   % returns 1x2, transpose to column
            q0  = q0(:);
        otherwise
            error('ic_mode must be "joint" or "task".');
    end
    qdot0 = (rand(2,1)-0.5)*2*dqdot0_max;

    xk_s = [q0; qdot0];  % SKOOPI state
    xk_b = xk_s;         % Baseline state

    % logs for this IC
    X_skoopi = zeros(4,T);  U_skoopi = zeros(2,T);
    X_base   = zeros(4,T);  U_base   = zeros(2,T);
    X_skoopi(:,1) = xk_s;   X_base(:,1) = xk_b;

    % simulate
    for idx = 1:T-1
        t_now = t_vector(idx);

        % gains (clamp to end if idx exceeds lists)
        if idx <= numel(K_fb_list1), Kfb2 = K_fb_list1{idx}; else, Kfb2 = K_fb_list1{end}; end
        if idx <= numel(K_ff_list1), Kff2 = K_ff_list1{idx}; else, Kff2 = K_ff_list1{end}; end
        if idx <= numel(K_fb_list2), Kfb1 = K_fb_list2{idx}; else, Kfb1 = K_fb_list2{end}; end
        if idx <= numel(K_ff_list2), Kff1 = K_ff_list2{idx}; else, Kff1 = K_ff_list2{end}; end

        % -------- SKOOPI control (psi-space) --------
        phi = compute_path_integrals(xk_s, dynamics, sys_info);
        ph  = phi.phi_x_op(:);
        ph  = [real(ph(1)); imag(ph(1)); real(ph(3)); imag(ph(3))];
        psi_x = xk_s + WinvT * ph;

        % desired psi (precomputed) or x_ref fallback
        if exist('psi_ref_series','var') && size(psi_ref_series,1)==4
            psi_d = psi_ref_series(:,idx);
        else
            psi_d = x_ref_series(:,idx);
        end

        u2 = -Kfb2 * (psi_x - psi_d) + Kff2;
        U_skoopi(:,idx) = u2;

        dx_s = dynamics_planarRR(xk_s, u2, sys_info);
        xk_s = xk_s + dt * dx_s;
%         xk_s(1:2) = mod(xk_s(1:2), 2*pi);
        X_skoopi(:,idx+1) = xk_s;

        % -------- Baseline control (x-space) --------
        x_d = x_ref_series(:,idx);
        u1  = -Kfb1 * (xk_b - x_d) + Kff1;
        U_base(:,idx) = u1;

        sys_info_b = planarRR_info(t_now, xk_b, u1, robot_params, sys_params);
        dx_b = dynamics_planarRR(xk_b, u1, sys_info_b);
        xk_b = xk_b + dt * dx_b;
%         xk_b(1:2) = mod(xk_b(1:2), 2*pi);
        X_base(:,idx+1) = xk_b;
    end

    % store this IC
    X_skoopi_all{ii} = X_skoopi;
    U_skoopi_all{ii} = U_skoopi;
    X_base_all{ii}   = X_base;
    U_base_all{ii}   = U_base;
end
delete_valid_waitbars();

%% ============================== plots ==============================
alpha_line = 0.5;     % transparency for individual trials
lw_thick   = 2;       % thicker legend proxies

% color palette proxies
fig_tmp = figure('visible','off'); ax_tmp = axes(fig_tmp);
C = ax_tmp.ColorOrder; close(fig_tmp);
col_s = C(1,:);   % SKOOPI color
col_b = C(2,:);   % LQR-FT color

labs  = {'$q_1$','$q_2$','$\dot q_1$','$\dot q_2$'};
ulabs = {'$\tau_1$','$\tau_2$'};

fig = figure('Name','States & Controls: SKOOPI vs LQR-FT (10 ICs, 4x4 grid)','Color','w');
set(fig, ...
    'DefaultAxesFontSize',20, ...
    'DefaultTextFontSize',20, ...
    'DefaultLegendFontSize',20, ...
    'DefaultTextInterpreter','latex', ...
    'DefaultAxesTickLabelInterpreter','latex', ...
    'DefaultLegendInterpreter','latex', ...
    'DefaultAxesTitleFontSizeMultiplier',1, ...
    'DefaultAxesLabelFontSizeMultiplier',1);
% -------- States on subplots [1 2 5 6] --------
state_slots = [1 2 5 6];
for kplot = 1:4
    subplot(4,4,state_slots(kplot)); hold on; grid on; box on
    % all ICs (semi-transparent)
    for ii = 1:N_ic
        Xi_s = X_skoopi_all{ii};
        Xi_b = X_base_all{ii};
        plot(t_vector, Xi_s(kplot,:), 'LineWidth', 2, ...
            'Color', [col_s alpha_line], 'HandleVisibility','off');
        plot(t_vector, Xi_b(kplot,:), 'LineWidth', 2, ...
            'Color', [col_b alpha_line], 'HandleVisibility','off');
    end
    % reference
    if exist('x_ref_series','var') && ~isempty(x_ref_series)
        plot(t_vector, x_ref_series(kplot,:), 'k-.', 'LineWidth', 1, 'DisplayName','reference');
    end
    % legend proxies
    plot(nan, nan, '-', 'Color', col_s, 'LineWidth', lw_thick, 'DisplayName','SKOOPI (10 ICs)');
    plot(nan, nan, '-', 'Color', col_b, 'LineWidth', lw_thick, 'DisplayName','LQR-FT (10 ICs)');

    xlabel('time (s)'); ylabel(labs{kplot});
    legend('Location','best');
    xlim([0,5])
end

% -------- Controls on subplots [9 10] --------
ctrl_slots = [9 10];
for kplot = 1:2
    subplot(4,4,ctrl_slots(kplot)); hold on; grid on; box on
    for ii = 1:N_ic
        Ui_s = U_skoopi_all{ii};
        Ui_b = U_base_all{ii};
        plot(t_vector, Ui_s(kplot,:), 'LineWidth', 2, ...
            'Color', [col_s alpha_line], 'HandleVisibility','off');
        plot(t_vector, Ui_b(kplot,:), 'LineWidth', 2, ...
            'Color', [col_b alpha_line], 'HandleVisibility','off');
    end
    % legend proxies
    plot(nan, nan, '-', 'Color', col_s, 'LineWidth', lw_thick, 'DisplayName','SKOOPI (10 ICs)');
    plot(nan, nan, '-', 'Color', col_b, 'LineWidth', lw_thick, 'DisplayName','LQR-FT (10 ICs)');

    xlabel('time (s)'); ylabel(ulabs{kplot});
    legend('Location','best');
    xlim([0,5])
end

%% ============================== EE trajectories (task space) ==============================
% Plots: reference circle + 10 IC trajectories in (x,y) for SKOOPI and LQR-FT

joint_control.t = t_vector;
q_skoopi = [X_skoopi(1,:)', X_skoopi(2,:)'];
plot_figure(robot_params,navigation_params,euler_params,joint_control,q_skoopi);

q_lqr = [X_base(1,:)', X_base(2,:)'];
plot_figure(robot_params,navigation_params,euler_params,joint_control,q_lqr);

%%
alpha_line = 0.50;      % transparency for trajectories
lw_ref     = 2;

% colors (reuse earlier palette)
fig_tmp = figure('visible','off'); ax_tmp = axes(fig_tmp);
C = ax_tmp.ColorOrder; close(fig_tmp);

col_s = C(1,:);   % SKOOPI color
col_b = C(2,:);   % LQR-FT color

% reference EE trajectory (2 x T)
XY_ref = zeros(2, numel(t_vector));
for k = 1:numel(t_vector)
    xrk = x_ref_f(t_vector(k));           % may be 1x2 or 2x1
    XY_ref(:,k) = xrk(:);
end

figure('Name','End-Effector Trajectories (Task Space)','Color','w');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on; grid on

% all SKOOPI ICs
for ii = 1:N_ic
    Xi_s = X_skoopi_all{ii};            % 4 x T
    XYs  = fk_planarRR(Xi_s(1:2,:)-[pi/2; 0], robot_params);   % 2 x T
    plot(XYs(1,:), XYs(2,:), 'LineWidth', 2, ...
        'Color', [col_s alpha_line], 'HandleVisibility','off');
    % start markers (hidden from legend)
    plot(XYs(1,1), XYs(2,1), 'ko', 'MarkerSize', 5);
end

for ii = 1:N_ic
    Xi_b = X_base_all{ii};              % 4 x T
    XYb  = fk_planarRR(Xi_b(1:2,:)-[pi/2; 0], robot_params);   % 2 x T
    plot(XYb(1,:), XYb(2,:),'-', 'LineWidth', 2, ...
        'Color', [col_b alpha_line], 'HandleVisibility','off');
    plot(XYb(1,1), XYb(2,1), 'ko', 'MarkerSize', 5, 'HandleVisibility','off');
end

% reference
plot(XY_ref(1,:), XY_ref(2,:), 'k--', 'LineWidth', lw_ref, 'DisplayName','reference');
box on
axis square
% legend proxy
plot(nan, nan, '-', 'Color', col_s, 'LineWidth', 2, 'DisplayName','SKOOPI');
plot(nan, nan, '-', 'Color', col_b, 'LineWidth', 2, 'DisplayName','FT-LQR');
xlabel('x [m]'); ylabel('y [m]'); axis equal; %legend('Location','best');


%% ============================== metrics: RMSE & NRMSE ==============================
% Definitions:
%   RMSE_i = sqrt(mean_t ( (y_i(t) - y_i,ref(t))^2 ))
%   NRMSE_i = RMSE_i / range(y_i,ref), averaged across components
% Notes:
%   - State error uses angle-wrap for q1,q2 (dims 1:2) via atan2(sinΔ,cosΔ)
%   - Control error is measured vs u_ref_series

assert(N_ic >= 1, 'No ICs were simulated.');

% ---- per-IC aggregates (mean over state dims / control dims) ----
rmse_state_LQR   = zeros(N_ic,1);
rmse_state_SKO   = zeros(N_ic,1);

for ii = 1:N_ic
    % ----- states (x = [q1 q2 qd1 qd2]) -----
    Xi_b = X_base_all{ii};      % 4 x T
    Xi_s = X_skoopi_all{ii};    % 4 x T

    [rmse_dims_b, rmse_mean_b] = compute_rmse_nrmse_states(Xi_b, x_ref_series);
    [rmse_dims_s, rmse_mean_s] = compute_rmse_nrmse_states(Xi_s, x_ref_series);

    rmse_state_LQR(ii)  = rmse_mean_b;
    rmse_state_SKO(ii)  = rmse_mean_s;
end

% ---- aggregate across ICs (report mean ± std) ----
m_rmse_x_LQR  = mean(rmse_state_LQR);   s_rmse_x_LQR  = std(rmse_state_LQR);
m_rmse_x_SKO  = mean(rmse_state_SKO);   s_rmse_x_SKO  = std(rmse_state_SKO);

% ---- SKOOPI improvement over LQR-FT ----
imp_rmse_state  = 100 * (1 - m_rmse_x_SKO / max(m_rmse_x_LQR, eps));

% ---- pretty print ----
fprintf('\n================ Metrics (over %d ICs) ================\n', N_ic);

fprintf('State RMSE (mean ± std):   LQR-FT = %.4f ± %.4f | SKOOPI = %.4f ± %.4f\n', ...
        m_rmse_x_LQR, s_rmse_x_LQR, m_rmse_x_SKO, s_rmse_x_SKO);
% fprintf('   → SKOOPI improvement: RMSE = %.2f%%', imp_rmse_state);


%% animations
% skoopi
X_skoopi1 = X_skoopi_all{1};
q_skoopi = [X_skoopi1(1,:)', X_skoopi1(2,:)'];
animate_planarRR(t_vector,q_skoopi,[],[xc; yc])

%% ============================== helpers ==============================
function delete_valid_waitbars()
    F = findall(0,'type','figure','tag','TMWWaitbar');
    if ~isempty(F), delete(F); end
end

function [rmse_dims, rmse_mean] = compute_rmse_nrmse_states(X, Xref)
% State metrics with angle-wrap for q1,q2 (dims 1:2).
% Returns per-dim RMSE/NRMSE, plus simple averages across dims.
    assert(all(size(X)==size(Xref)), 'Dimension mismatch.');
    [n, ~] = size(X);
    rmse_dims  = zeros(n,1);
    for i = 1:n
        if i <= 2
            e = wrapdiff(X(i,:), Xref(i,:));   % joint angle error
        else
            e = X(i,:) - Xref(i,:);
        end
        rmse_dims(i) = sqrt(mean(e.^2));
        denom = max(range(Xref(i,:)), 1e-12);  % NRMSE normalized by reference range
    end
    rmse_mean  = mean(rmse_dims);
end


function r = range(v)
% Simple range helper (max - min)
    r = max(v) - min(v);
end

function d = wrapdiff(a, b)
% Angle difference wrapped to [-pi, pi]
    d = atan2(sin(a - b), cos(a - b));
end

function XY = fk_planarRR(Q, robot_params)
% Forward kinematics for planar 2R arm.
% Q: 2 x T array of joint angles [q1; q2] (radians).
% Returns XY: 2 x T of end-effector positions [x; y].
%
% Expects link lengths in robot_params as one of:
%   l1/l2  OR  L1/L2  OR  a1/a2
    L1 = get_link_len(robot_params, {'l1','L1','a1'});
    L2 = get_link_len(robot_params, {'l2','L2','a2'});

    q1 = Q(1,:);  q2 = Q(2,:);
    x = L1*cos(q1) + L2*cos(q1 + q2);
    y = L1*sin(q1) + L2*sin(q1 + q2);
    XY = [x; y];
end

function L = get_link_len(rp, fields)
% Try multiple candidate field names for a link length.
    for k = 1:numel(fields)
        f = fields{k};
        if isfield(rp, f)
            L = rp.(f);
            return;
        end
    end
    error('Link length not found. Please add robot_params.l1/l2 (or L1/L2, a1/a2).');
end

