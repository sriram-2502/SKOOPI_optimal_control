%% SKOOPI for robotic arm — one-shot (no chunking), circular ref + ID feedforward
clc; clear; close all
rng(0);

% ----- paths -----
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
[robot_params,euler_params,navigation_params,lqr_params,animate_params] = get_params_planarRR(x);

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

% make W real
W = [real(W(:,1)), imag(W(:,2)), real(W(:,3)), imag(W(:,4))];

% ===== Controller weights =====
Q  = 1e2*diag([1 1 0.01 0.01]);
R  = 1e3*diag([1 1]);
QN = 1e1*Q;

%% ===== Task-space circular reference (downward workspace) =====
xc = 0.0;  yc = 0.0;    % circle center (m)
r  = 1.5;               % radius (m)
omega = 0.2;               % rad/s
time_varying_goal = [xc + r*cos(omega*t), yc + r*sin(omega*t)];
matlabFunction(time_varying_goal, 'File', 'planarRR/x_ref_f', 'Vars', {t}, 'Optimize', false);

% IK and time-derivatives symbols → functions
navigation_params.x_goal_f       = IK_planarRR(time_varying_goal);
navigation_params.x_vel_goal_f   = diff(navigation_params.x_goal_f);
navigation_params.x_accel_goal_f = diff(navigation_params.x_vel_goal_f);

matlabFunction(navigation_params.x_goal_f,       'File', 'planarRR/xd_f',      'Vars', {t}, 'Optimize', false);
matlabFunction(navigation_params.x_vel_goal_f,   'File', 'planarRR/x_dot_d_f', 'Vars', {t}, 'Optimize', false);
matlabFunction(navigation_params.x_accel_goal_f, 'File', 'planarRR/x_ddot_d_f','Vars', {t}, 'Optimize', false);
matlabFunction(jacobian(navigation_params.x_goal_f, t), 'File', 'planarRR/jac_xd_f', 'Vars', {t}, 'Optimize', false);
rehash path

% ===== Discretize reference & build joint-space (q, qd, qdd) =====
t_vector = (0:euler_params.n_steps).*euler_params.step_size;
dt       = euler_params.step_size;
T        = numel(t_vector);
dt_sim   = euler_params.step_size;
t_end    = (T-1)*dt_sim;

% Evaluate IK along the circle
q_goal = zeros(T, 2);
for i = 1:T
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

% Store for later (if other code uses navigation_params.*)
navigation_params.q_goal      = q_goal;
navigation_params.q_dot_goal  = q_dot_goal;
navigation_params.q_ddot_goal = q_ddot_goal;

if show_diagnositcs
    figure('Name','Joint references','Color','w');
    subplot(3,1,1); plot(t_vector, q_goal,'-o');      grid on; title('q_d')
    subplot(3,1,2); plot(t_vector, q_dot_goal,'-o');  grid on; title('\dot q_d')
    subplot(3,1,3); plot(t_vector, q_ddot_goal,'-o'); grid on; title('\ddot q_d')
end

% ===== Feedforward torques u_ref(t) and eigenfunctions =====
u_ref_series   = zeros(2, T);
psi_ref_series = zeros(4, T);
w_bar = waitbar(0, '1', 'Name', sprintf('geting reference'), ...
    'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
for i = 1:T
    waitbar(i / T, w_bar, sprintf('%d / %d iter', i, T));
    qd  = q_goal(i,:).';
    dqd = q_dot_goal(i,:).';
    ddq = q_ddot_goal(i,:).';
    x_ref = [q_goal(i,:).'; q_dot_goal(i,:).' ];
    
    % compute phi and psi along ref traj
    phi = compute_path_integrals(x_ref, dynamics, sys_info);
    phi_x_op = phi.phi_x_op';
    phi_x_op = [real(phi_x_op(1)); imag(phi_x_op(1));real(phi_x_op(3)); imag(phi_x_op(3));];
    psi_ref_series(:,i) = x_ref + inv(W') * phi_x_op;

    grad_phi_x_op = compute_gradients(phi);
    grad_phi_x_op = [
        real(grad_phi_x_op(1,:));  % real part of row 1
        imag(grad_phi_x_op(1,:));  % imag part of row 1
        real(grad_phi_x_op(3,:));  % real part of row 3
        imag(grad_phi_x_op(3,:));  % imag part of row 3
    ];
    
    % feedforward torques
    params = planarRR_info(t, x_ref, [0;0], robot_params, sys_params);
    M_t = params.M; C_t = params.C; G_t = params.G; D_t = params.D;
    u_ref_series(:,i) = M_t*ddq + C_t + G_t + D_t*dqd;
    
end

% ===== Assemble reference state series & g-lists (one shot) =====
x_ref_series = zeros(4, T);
g_list_skoopi = cell(1, T);
g_list_base   = cell(1, T);

for i = 1:T
    x_ref_series(:,i) = [ q_goal(i,:).'; q_dot_goal(i,:).' ];
    g_list_skoopi{i}  = g_fun( x_ref_series(:,i) );
    g_list_base{i}    = B;    % baseline uses constant B
end

x_ref = x_ref_series;
u_ref = u_ref_series;

% ===== Solve tracking DRE once (no chunking) =====
[L_b, P_b, K_ff_list1, K_fb_list1, sx_b] = solve_DRE_tracking( ...
    x_ref, u_ref, A, g_list_base, Q, R, QN, t_end, dt_sim);

[L_s, P_s, K_ff_list2, K_fb_list2, sx_s] = solve_DRE_tracking( ...
    x_ref, u_ref, A, g_list_skoopi, Q, R, QN, t_end, dt_sim);

% diagnostics: show u_ref
if show_diagnositcs
    figure('Name','u_{ref}','Color','w');
    plot(t_vector, u_ref_series(1,:),'-', t_vector, u_ref_series(2,:),'--'); grid on
    xlabel('t [s]'); ylabel('\tau [N·m]'); legend('\tau_1^{ref}','\tau_2^{ref}')     
end

% animate reference
animate_planarRR(t_vector,q_goal,[],[xc; yc])

% delete progress bar
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% ======================= simulation loop =======================
% time grid
T  = numel(t_vector);
dt = euler_params.step_size;

% precompute inv(W') once
WinvT = inv(W.');

% init states from the first reference waypoint
xk_s = [q_goal(1,:).'; q_dot_goal(1,:).'];   % SKOOPI (psi-space) controller state
xk_b = xk_s;                                  % Baseline (x-space) controller state

% logs
X_skoopi = zeros(4,T);  U_skoopi = zeros(2,T);
X_base   = zeros(4,T);  U_base   = zeros(2,T);
X_skoopi(:,1) = xk_s;   X_base(:,1) = xk_b;

w_bar = waitbar(0, '1', 'Name', sprintf('running simulation'), ...
    'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
for idx = 1:T-1
    t_now = t_vector(idx);
    waitbar(idx / T, w_bar, sprintf('%d / %d iter', idx, T));
    % -------- SKOOPI control (psi-space) --------
    % gains (clamp index to last if lists shorter)
    if idx <= numel(K_fb_list1), Kfb2 = K_fb_list1{idx}; else, Kfb2 = K_fb_list1{end}; end
    if idx <= numel(K_ff_list1), Kff2 = K_ff_list1{idx}; else, Kff2 = K_ff_list1{end}; end

    % eigenfunctions & psi mapping
    phi = compute_path_integrals(xk_s, dynamics, sys_info);
    ph = phi.phi_x_op(:);
    ph = [real(ph(1)); imag(ph(1)); real(ph(3)); imag(ph(3))];
    psi_x = xk_s + WinvT * ph;

    % desired psi (use precomputed psi_ref_series if you built it; else fall back to x_ref)
    if exist('psi_ref_series','var') && size(psi_ref_series,1)==4
        psi_d = psi_ref_series(:,idx);
    else
        psi_d = x_ref_series(:,idx);
    end

    % control and simulate one Euler step
    u2 = -Kfb2 * (psi_x - psi_d) + Kff2;
    U_skoopi(:,idx) = u2;

    dx_s = dynamics_planarRR(xk_s, u2, sys_info);
    xk_s = xk_s + dt * dx_s;
%     xk_s(1:2) = mod(xk_s(1:2), 2*pi);  % wrap joints
    X_skoopi(:,idx+1) = xk_s;

    % -------- Baseline control (x-space) --------
    if ~isempty(K_fb_list2) && ~isempty(K_ff_list2)
        if idx <= numel(K_fb_list2), Kfb1 = K_fb_list2{idx}; else, Kfb1 = K_fb_list2{end}; end
        if idx <= numel(K_ff_list2), Kff1 = K_ff_list2{idx}; else, Kff1 = K_ff_list2{end}; end

        x_d = x_ref_series(:,idx);
        u1  = -Kfb1 * (xk_b - x_d) + Kff1;
        U_base(:,idx) = u1;

        sys_info_b = planarRR_info(t_now, xk_b, u1, robot_params, sys_params);
        dx_b = dynamics_planarRR(xk_b, u1, sys_info_b);
        xk_b = xk_b + dt * dx_b;
%         xk_b(1:2) = mod(xk_b(1:2), 2*pi);
        X_base(:,idx+1) = xk_b;
    end
end

% delete progress bar
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% plots

figure('Name','SKOOPI vs LQR-FT — States & Controls (4x4 slots)','Color','w');
labs  = {'$q_1$','$q_2$','$\dot q_1$','$\dot q_2$'};
ulabs = {'$\tau_1$','$\tau_2$'};

% --- States on (1, 2, 5, 6) ---
state_slots = [1 2 5 6];
for kplot = 1:4
    subplot(4,4,state_slots(kplot)); hold on; grid on; box on
    plot(t_vector, X_skoopi(kplot,:), 'LineWidth', 2, 'DisplayName','SKOOPI');
    plot(t_vector, X_base(kplot,:),   'LineWidth', 2, 'DisplayName','LQR-FT');
    if exist('x_ref_series','var') && ~isempty(x_ref_series)
        plot(t_vector, x_ref_series(kplot,:), 'k--', 'LineWidth', 1, 'DisplayName','reference');
    end
    ylabel(labs{kplot}, 'Interpreter','latex');
    if ismember(state_slots(kplot), [5 6])  % second row
        xlabel('t [s]','Interpreter','latex');
    end
%     legend('Location','northoutside', 'NumColumns',3);
end

% --- Controls on (9, 10) ---
ctrl_slots = [9 10];
for kplot = 1:2
    subplot(4,4,ctrl_slots(kplot)); hold on; grid on; box on
    plot(t_vector, U_skoopi(kplot,:), 'LineWidth', 2, 'DisplayName','SKOOPI');
    plot(t_vector, U_base(kplot,:),   'LineWidth', 2, 'DisplayName','LQR-FT');
    xlabel('t [s]','Interpreter','latex');
    ylabel(ulabs{kplot}, 'Interpreter','latex');
%     legend('Location','best');
end

%% find total costs
% ===== SKOOPI =====
T  = size(X_skoopi,2);
Ju_skoopi = 0;  Jx_skoopi = 0;

for k = 1:T
    x = X_skoopi(:,k);
    xref = x_ref(:,k);
    Jx_skoopi = Jx_skoopi + (x-xref)'*Q*(x-xref);          % + *dt if you want dt-weighted
    if k <= size(U_skoopi,2)
        u = U_skoopi(:,k);
        Ju_skoopi = Ju_skoopi + (u.'*R*u);      % + *dt if you want dt-weighted
    end
end

fprintf('SKOOPI: Jx = %.6g, Ju = %.6g, Total = %.6g\n', ...
        Jx_skoopi, Ju_skoopi, Jx_skoopi + Ju_skoopi);

% ===== LQR-FT =====
if exist('X_base','var') && ~isempty(X_base)
    T2 = size(X_base,2);
    Ju_base = 0;  Jx_base = 0;

    for k = 1:T2
        x = X_base(:,k);
        xref = x_ref(:,k);
        Jx_base = Jx_base + (x-xref)'*Q*(x-xref);          % + *dt if dt-weighted
        if k <= size(U_base,2)
            u = U_base(:,k);
            Ju_base = Ju_base + (u.'*R*u);      % + *dt if dt-weighted
        end
    end

    fprintf('LQR-FT: Jx = %.6g, Ju = %.6g, Total = %.6g\n', ...
            Jx_base, Ju_base, Jx_base + Ju_base);
end


%% animate skpoopi
q_skoopi = [X_skoopi(1,:)', X_skoopi(2,:)'];
animate_planarRR(t_vector,q_skoopi,[],[xc; yc])

%% animate LQR-FT
q_lqr = [X_base(1,:)', X_base(2,:)'];
animate_planarRR(t_vector,q_lqr,[],[xc; yc])
