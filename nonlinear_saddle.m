%% eigenfunctions for 2D analytical system with saddle-point origin
clc; clear; close all;
set(0,'DefaultLineLineWidth',2) %linewidh on plots
set(0,'defaultfigurecolor',[1 1 1])

%% system description
% nonlinear ode x_dot = f(x)
% ---------- add paths (expects your existing folders/functions) ----------
addpath('dynamics', 'baseline_control', 'eigfun_control', ...
        'compute_eigfuns', 'utils', 'animations');
rng(0)

% ---------- system setup (saddle) ----------
sys_params.use_stable         = false;   % locally stable
sys_params.use_unstable       = false;    % locally unstable
use_reverse       = false;

% Load dynamics & system info
sys_info   = nonlinear_sys_info(sys_params);
dynamics   = @dynamics_nonlinear;
n_states   = sys_info.state_dim;
n_ctrl     = sys_info.ctrl_dim;

x = sym('x',[n_states,1],'real');

% define eigenvalues
Lam = [1;-3];

phi_1_ana = @(x1,x2)sin(x1)-2*x2;
phi_2_ana = @(x1,x2)x1+sin(x2);

Phi_ana = @(x1,x2)[phi_1_ana(x1,x2); phi_2_ana(x1,x2)];
phi_x = Phi_ana(x(1),x(2));
dPhi_dx = simplify(jacobian(phi_x,x));
f_x = inv(dPhi_dx)*diag(Lam)*phi_x;
f = matlabFunction(f_x,'vars',{x(1), x(2)});
g_x = [2+cos(x(2)); -1+cos(x(1))];
g = matlabFunction(g_x,'vars',{x(1), x(2)});

%% Get the eigenvectors from the analytical eigenfunctions
xEq_struct = solve(f_x==0);
xEq = nan(n_states,1);
xEq_fieldNames = fieldnames(xEq_struct);
for i = 1:length(xEq_fieldNames)
    x_idx = xEq_fieldNames{i}; % Get the axis name from the cell array
    xEq(i,1) = xEq_struct.(x_idx); % Access the field data using dynamic field names
end

% Find the linear part of the eigenfunction
phi_x_linMat = double(subs(dPhi_dx,x,xEq));

% Separate the eigenvectors
w = {};
for i = 1:size(phi_x_linMat,2)
    w{i} = phi_x_linMat(i,:)';
end
W = [w{1},w{2}];

%% Get linearization and non-linear part
A = double(subs(jacobian(f_x),[x(1) x(2)]',xEq));
B = eval(subs((g_x),[x(1) x(2)]',xEq));

% define nonlinear part x_dot = Ax + fn(x)
fn = f_x - A*x;
dfn_dx = simplify(jacobian(fn,x));
% check if the nonlinear part is zero at origin
double(subs(dfn_dx,[x(1) x(2)]',xEq));

% define matlab functions for W'Fn
wFn = {};
for i = 1:length(w)
    wFn{i} = matlabFunction(w{i}'*fn,'vars',{x(1), x(2)});
end

%% get LQR gains for stabilization
% ---------- LQR (baseline, infinite-horizon) ----------
Q_baseline = 1e2*eye(n_states);
R_baseline = 1e2*eye(n_ctrl);
lqr_params_baseline = get_lqr(A, B, Q_baseline, R_baseline);

% ---------- SKOOPI / psi coordinates (your current approach uses A,B) ----------
Q  = Q_baseline;
R  = R_baseline;
QN = 1e1*Q_baseline;

A_transformed = A;
B_transformed = B;
Q_transformed = Q_baseline;
lqr_params_transformed = get_lqr(A_transformed, B_transformed, Q_transformed, R);

% ---------- simulation timing ----------
dt_sim = 0.01;
t_end  = 3.0;
t_grid = 0:dt_sim:t_end;         % 1 x N
Nsteps = numel(t_grid);

% ---------- Precompute finite-horizon gains (LQR-FT) ----------
x_desired = zeros(n_states,1);
g_list    = repmat({B}, 1, Nsteps-1);           % constant B each step
x_ref_mat = zeros(n_states, Nsteps-1);          % zero reference
u_ref_mat = zeros(n_ctrl,   Nsteps-1);          % zero reference control

% Solve DRE in original coordinates (feedback + feedforward)
[P_list, s_x_list2, K_ff_list_stabilization, K_fb_list_stabilization, t_grid_dre] = ...
    solve_DRE(x_ref_mat, u_ref_mat, A, g_list, Q, R, QN, t_end, dt_sim);

% ---------- Precompute Riccati solution for SKOOPI (time-varying P) ----------
[~, P_riccati] = compute_riccati(lqr_params_transformed, t_grid);
% P_riccati: Nsteps x (n_states^2); row i is vec(P_i)


%% for tracking
x_desired = [5;5];
psi_desired_list  = {}; h_desired_list = {};  % for SKOOPI (psi-coords)
x_desired_list = {};
u_desired_list = {};
g_list1 = {}; g_list2 = {};   % g along ref for SKOOPI vs baseline (B)

for t_sim = dt_sim:dt_sim:t_end
    % open-loop input (feel free to change)
    u_o = 1*cos(0.5*t_sim);

    % g(x) along reference
    g_t = sys_info.dynamics_g(x_desired);

    % store psi_d for SKOOPI
    phi_desired = eval(subs((phi_x),[x(1) x(2)]',x_desired));
    psi_desired_list{end+1}       = inv(W')*phi_desired; 
%     psi_desired_list{end+1}       = x_desired + inv(W') * sys_info.transform_fun_analytical(x_desired); 
%     h_desired_list{end+1}       = sys_info.transform_fun_analytical(x_desired); 
    h_desired_list{end+1}       = phi_desired-W'*x_desired; 
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
[~, ~, K_ff_list_tracking1, K_fb_list_tracking1, ~] = ...
    solve_DRE_tracking(x_ref, u_ref, A, g_list1, Q, R, QN, t_end, dt_sim);

% Baseline (original coords) tracking gains (LQR-FT)
[~, ~, K_ff_list_tracking2, K_fb_list_tracking2, ~] = ...
    solve_DRE_tracking(x_ref, u_ref, A, g_list2, Q, R, QN, t_end, dt_sim);

%% simulation loop

% ---------- IC sampling (uniform in a disk around (1,2)) ----------
N_ic   = 1;
center = [5; 5];
radius = 1;

theta = 2*pi*rand(N_ic,1);
rad   = radius*sqrt(rand(N_ic,1));
% ICs   = center.' + [rad.*cos(theta), rad.*sin(theta)];  % N_ic x 2
ICs = center';

% ---------- arrays to collect time traces for each method ----------
X1_all = cell(N_ic,1);   % SKOOPI states
U1_all = cell(N_ic,1);   % SKOOPI controls (Nsteps-1 x n_ctrl)

% logs for eigenfunctions (rows = time, cols = eigenfunction index)
phi_est_log   = zeros(Nsteps, n_states);
phi_anal_log  = zeros(Nsteps, n_states);
psi_est_log   = zeros(Nsteps, n_states);
psi_anal_log  = zeros(Nsteps, n_states);
h_est_log   = zeros(Nsteps, n_states);
h_anal_log  = zeros(Nsteps, n_states);

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

% eigfn setup
ballR = 2000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6, ...
    'events',@(t, x)offFrame(t, x, 2000));
tEnd = 20;
tspan = [0 tEnd];

for k = 1:N_ic
    k
    % initial states for each controller path
    x1 = ICs(k,:);   % SKOPPI
    x2 = ICs(k,:);   % SKOPPI

    % logs for this IC (states)
    X1 = x1;
    X2 = x2;

    % logs for this IC (controls)
    U1 = [];
    U2 = [];

    iter = 1;  % time index for time-varying gains

    % simulate forward in time
    for t_sim = dt_sim:dt_sim:t_end
        t_sim
        % compute eigenfunctions
        [phi_x_op,h_x_op] = compute_path_integrals(Lam,w,wFn,f,x1');
        phi_analytical    = eval(subs((phi_x),[x(1) x(2)]',x1'));

        % psi_x           = x1' + inv(W') * sys_info.transform_fun_analytical(x1');
        psi_x             = inv(W')*phi_analytical(:);
        psi_x_op          = inv(W')*phi_x_op';
        h_x               = phi_analytical-W'*x1';
        % h_x             = sys_info.transform_fun_analytical(x1');
       
        % logs
        phi_est_log(iter,:)   = phi_x_op(:); % estimated
        phi_anal_log(iter,:)  = phi_analytical(:)'; % analytical
        psi_est_log(iter,:)   = psi_x_op; % estimated
        psi_anal_log(iter,:)  = psi_x; % analytical 
        h_est_log(iter,:)     = h_x_op'; % estimated
        h_anal_log(iter,:)    = h_x; % analytical
 
        idx  = min(iter, numel(K_fb_list_stabilization)); 

        % --- Stabilization Controller: finite-horizon SKOOPI --- 
%         K_fb = K_fb_list_stabilization{idx}; 
%         K_ff = K_ff_list_stabilization{idx}; 
%         u1   = -K_fb * (psi_x_anal) - K_ff;
%         u2   = -K_fb * phi_x_op' - K_ff;

        % --- Tracking Controller: finite-horizon SKOOPI --- 
        Kfb1 = K_fb_list_tracking1{idx};  Kff1 = K_ff_list_tracking1{idx};
        psi_d = psi_desired_list{idx};  h_d = h_desired_list{idx};

        % control test
%         u1     = -Kfb1 * (psi_x - psi_d) + Kff1; % does not work
        u1     = -Kfb1 * (psi_x_op - psi_d) + Kff1;
        % h_x_anal = sys_info.transform_fun_analytical(x1');
%         u1       = -Kfb1 * (h_x - h_d) + Kff1; % works best        
        % u1     = -Kfb1 * (h_x_op - h_d) + Kff1;

        % --- Tracking Controller: finite-horizon LQR ---
        Kfb2 = K_fb_list_tracking2{idx};  Kff2 = K_ff_list_tracking2{idx};
        xd = x_ref(:, min(idx, size(x_ref,2)));
        u2 = -Kfb2 * (x1' - xd) + Kff2;

        if(u1>1e3)
            disp("path integral control blowing up")
        end
        % switch to LQR when close / invalid
        if ~isfinite(u1) || any(abs(u1)>1e3) || norm(x1) <= 1e-3
            u1 = -lqr_params_baseline.K_lqr * x1';
        end

        % --- integrate one step (your RK4) ---
        use_reverse = false;
        x1n = rk4(dynamics, dt_sim, x1', u1, use_reverse, sys_info);
        x2n = rk4(dynamics, dt_sim, x2', u2, use_reverse, sys_info);

        % update states
        x1 = x1n'; x2 = x2n';
        X1 = [X1; x1]; X2 = [X2; x2];

        % update controls (row-wise)
        U1 = [U1; u1']; U2 = [U2; u2'];
        iter = iter + 1;
    end

    % store time traces
    X1_all{k} = X1;  U1_all{k} = U1;
    X2_all{k} = X2;  U2_all{k} = U2;
    
    % store eigfns
    Phi_est_all{k}    = phi_est_log;
    Phi_anal_all{k}   = phi_anal_log;
    Psi_est_all{k}    = psi_est_log;
    Psi_anal_all{k}   = psi_anal_log;
    H_est_all{k}    = h_est_log;
    H_anal_all{k}   = h_anal_log;

    % plot phase-portrait trajectories with alpha (RGBA)
    plot(axp, X2(:,1), X2(:,2), '-', 'Color', [C(1,:) alpha_line], 'HandleVisibility','off');   % LQR
    plot(axp, X1(:,1), X1(:,2), '-', 'Color', [C(2,:) alpha_line], 'HandleVisibility','off');   % SKOOPI
    plot(axp, x_ref(1,:), x_ref(2,:), '--', 'Color', 'k', 'HandleVisibility','off');   % ref
end

% Re-plot disk & center
tt = linspace(0,2*pi,200);
plot(axp, center(1)+radius*cos(tt), center(2)+radius*sin(tt), ':k');
plot(axp, center(1), center(2), 'ok', 'MarkerFaceColor','k');

% ---------- Legend (phase portrait) — proxies only ----------
[hLQRp, hSKOOPp, hREFp] = legendProxies(axp, C(1,:), C(2,:));
legend(axp, [hLQRp hSKOOPp hREFp], {'LQR-FT','SKOOPI','reference'}, ...
       'Interpreter','latex','Location','northoutside','NumColumns',3, ...
       'AutoUpdate','off');

%% ---- RMSE computation between estimated and analytical eigenfunctions ----
rmse_vals = zeros(N_ic,n_states);
nrmse_vals = zeros(N_ic,n_states);

for k = 1:N_ic
    phi_est  = Phi_est_all{k};
    phi_anal = Phi_anal_all{k};

    for i = 1:n_states
        err = phi_est(:,i) - phi_anal(:,i);
        rmse_vals(k,i) = sqrt(mean(err.^2));

        % Normalized RMSE relative to signal range
        range_phi = max(phi_anal(:,i)) - min(phi_anal(:,i));
        if range_phi > 1e-8
            nrmse_vals(k,i) = rmse_vals(k,i) / range_phi;
        else
            nrmse_vals(k,i) = NaN;
        end
    end
end

% Compute mean and std across ICs
rmse_mean  = mean(rmse_vals,1,'omitnan');
rmse_std   = std(rmse_vals,0,1,'omitnan');
nrmse_mean = mean(nrmse_vals,1,'omitnan');
nrmse_std  = std(nrmse_vals,0,1,'omitnan');

% Display results
disp('--- RMSE statistics between estimated and analytical eigenfunctions ---');
for i = 1:n_states
    fprintf('Eigenfunction %d: RMSE = %.4f ± %.4f,  NRMSE = %.4f ± %.4f\n', ...
        i, rmse_mean(i), rmse_std(i), nrmse_mean(i), nrmse_std(i));
end

%% ---- plots ----
t = t_grid(:);               % Nsteps x 1
tc = t(2:end);               % control timeline (Nsteps-1 x 1)

% plot traj
plot_method_time_and_controls(X2_all, U2_all, C(1,:), 'FT-LQR',    alpha_line, t, tc, n_ctrl);
plot_method_time_and_controls(X1_all, U1_all, C(2,:), 'SKOOPI',    alpha_line, t, tc, n_ctrl);

% --- plot eigfun
figure('Name','Eigenfunction comparison','Color','w');
subplot(4,4,[1 2]); hold on; grid on; box on;
for k = 1:N_ic
    plot(t, Phi_est_all{k}(:,1), '-', 'LineWidth', 1.5, ...
         'DisplayName','Estimated ($\hat{\phi}_1$)');
    plot(t, Phi_anal_all{k}(:,1), 'k--', 'LineWidth', 2, ...
         'DisplayName','Analytical ($\phi_1$)'); 
end
ylabel('$\phi_1$', 'Interpreter','latex');
% xlabel('time (s)','Interpreter','latex');
xticklabels([])
lgd = legend();
lgd.Interpreter='latex';
% ylim([-1,7])
% xticks([1 2 3 4 5])
% yticks([0 3 6])

subplot(4,4,[5 6]); hold on; grid on; box on;
for k = 1:N_ic
    plot(t, Phi_est_all{k}(:,2), '-', 'LineWidth', 1.5, ...
         'DisplayName','Estimated ($\hat{\phi}_2$)');
    plot(t, Phi_anal_all{k}(:,2), 'k--', 'LineWidth', 2, ...
         'DisplayName','Analytical ($\phi_2$)'); 
end
ylabel('$\phi_2$', 'Interpreter','latex');
xlabel('time (s)','Interpreter','latex');
lgd = legend();
lgd.Interpreter='latex';
% ylim([-20,10])

% --- plot psifun
subplot(4,4,[3 4]); hold on; grid on; box on;
for k = 1:N_ic
    plot(t, Psi_est_all{k}(:,1), '-', 'LineWidth', 1.5, ...
         'DisplayName','Estimated ($\hat{\psi}_1$)');
    plot(t, Psi_anal_all{k}(:,1), 'k--', 'LineWidth', 2, ...
         'DisplayName','Analytical ($\psi_1$)'); 
end
ylabel('$\psi_1$', 'Interpreter','latex');
% xlabel('time (s)','Interpreter','latex');
xticklabels([])
lgd = legend();
lgd.Interpreter='latex';
% ylim([-1,7])
% xticks([1 2 3 4 5])
% yticks([0 3 6])

subplot(4,4,[7 8]); hold on; grid on; box on;
for k = 1:N_ic
    plot(t, Psi_est_all{k}(:,2), '-', 'LineWidth', 1.5, ...
         'DisplayName','Estimated ($\hat{\psi}_2$)');
    plot(t, Psi_anal_all{k}(:,2), 'k--', 'LineWidth', 2, ...
         'DisplayName','Analytical ($\psi_2$)'); 
end
ylabel('$\psi_2$', 'Interpreter','latex');
xlabel('time (s)','Interpreter','latex');
lgd = legend();
lgd.Interpreter='latex';
% ylim([-20,10])

% --- plot h(x) fun
subplot(4,4,[9 10]); hold on; grid on; box on;
for k = 1:N_ic
    plot(t, H_est_all{k}(:,1), '-', 'LineWidth', 1.5, ...
         'DisplayName','Estimated ($\hat{h}_1$)');
    plot(t, H_anal_all{k}(:,1), 'k--', 'LineWidth', 2, ...
         'DisplayName','Analytical ($h_1$)'); 
end
ylabel('$h_1$', 'Interpreter','latex');
% xlabel('time (s)','Interpreter','latex');
xticklabels([])
lgd = legend();
lgd.Interpreter='latex';
% ylim([-1,7])
% xticks([1 2 3 4 5])
% yticks([0 3 6])

subplot(4,4,[13 14]); hold on; grid on; box on;
for k = 1:N_ic
    plot(t, H_est_all{k}(:,2), '-', 'LineWidth', 1.5, ...
         'DisplayName','Estimated ($\hat{h}_2$)');
    plot(t, H_anal_all{k}(:,2), 'k--', 'LineWidth', 2, ...
         'DisplayName','Analytical ($h_2$)'); 
end
ylabel('$h_2$', 'Interpreter','latex');
xlabel('time (s)','Interpreter','latex');
lgd = legend();
lgd.Interpreter='latex';
% ylim([-20,10])

%% ==================== local helpers ====================
% compute path integerals for saddle
function  [phi_PI,nl_phi] = compute_path_integrals(Lam,w,wFn,f,x_i)
    ballR = 2000;
    options = odeset('RelTol',1e-6,'AbsTol',1e-6, ...
        'events',@(t, x)offFrame(t, x, 2000));
    tEnd = 20;
    tspan = [0 tEnd];

    % simulate the trajectories (forward and backward)
    % % Forward time simulation for un-stable eigenfunction
    [tf,yf] = ode45(@(t,y)f(y(1),y(2)),tspan,x_i,options);
    % % Backward time simulation for stable eigenfunction
    [tr,yr] = ode45(@(t,y)f(y(1),y(2)),-tspan,x_i,options);
    
    for i = 1:length(Lam)
        ev_i = Lam(i);
        wFn_i = wFn{i};
        w_i = w{i};
        if real(ev_i) >= 0
            % un-stable eigenfunction
            if norm(yf(end,:))<=ballR+1e-6
                tE(i) = tf(end);
                nl_phi(i) = getIntegral(tf,yf,ev_i,wFn_i);
                phi_PI(i) = w_i'*x_i + nl_phi(i);
            end
        else
            % stable eigenfunction
            if norm(yr(end,:))<=ballR+1e-6
                tE(i) = tr(end);
                nl_phi(i) = getIntegral(tr,yr,ev_i,wFn_i);
                phi_PI(i) = w_i'*x_i + nl_phi(i);
            end
        end
    end
end

% To stop the excution at boundaty
function [value,isterminal,direction]=offFrame(~, Y, Dom)
value = (norm(Y)>Dom) | (norm(Y)<1e-6);
isterminal=1;
direction=0;
end

function integralVal = getIntegral(t,y,ev,wfn)
if ev>0
    iTraj2 = exp(-t(:)*ev).*wfn(y(:,1),y(:,2));%plot(t,iTraj2)
    integralVal = trapz(t,iTraj2,1);
elseif ev<0
    iTraj2 = exp(abs(t(:))*ev).*wfn(y(:,1),y(:,2));
    integralVal = trapz(t,iTraj2,1);
end
end

% plots helper
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
    plot(t, Xcell{1}(:,1), 'Color', color_rgb, 'LineWidth', 2, 'DisplayName',['$x_1$ -' name_str] );
    ylabel('$x_1(t)$','Interpreter','latex');
    legend('Interpreter','latex','Location','best');
    ylim auto; xlim([t(1) t(end)]);

    % ---- x2(t) in subplot(4,4,[5 6]) ----
    subplot(4,4,[5 6]); hold on; grid on; box on;
    for i = 1:numel(Xcell)
        plot(t, Xcell{i}(:,2), 'Color', [color_rgb a], 'HandleVisibility','off');
    end
    plot(t, Xcell{1}(:,2), 'Color', color_rgb, 'LineWidth', 2, 'DisplayName',['$x_2$ - ' name_str] );
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
    plot(tc, Ucell{1}(:,1), 'Color', color_rgb, 'LineWidth', 2, 'DisplayName',['$u$ -' name_str] );
    xlabel('time (s)','Interpreter','latex');
    ylabel('$u_1(t)$','Interpreter','latex');
    legend('Interpreter','latex','Location','best');
    xlim([tc(1) tc(end)]);
end

function [hLQR, hSKOOP, hREF] = legendProxies(ax, colLQR, colSKOOP)
    hLQR   = plot(ax, NaN,NaN,'-','Color',colLQR,   'LineWidth',2,   'HandleVisibility','on');
    hSKOOP = plot(ax, NaN,NaN,'-','Color',colSKOOP, 'LineWidth',2,   'HandleVisibility','on');
    hREF   = plot(ax, NaN,NaN,'--k','LineWidth',1.5,'HandleVisibility','on');
end