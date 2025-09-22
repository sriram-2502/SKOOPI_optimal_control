function eigen_fn = compute_eigen_fn(x_local, x_eqb, dynamics, D, W, sys_info)
% parse inputs
n_dim = length(x_eqb);

% get modified linear system
if(sys_info.use_stable)
    A = sys_info.A_stable;
elseif(sys_info.use_unstable)
    A = sys_info.A_unstable;
else
    A = sys_info.A;
end
dynamics_linearized = @(x,u,sys_params) A*x;

% parse path integral setup
use_cart_pole=false;
use_quadrotor=false;
use_non_linear=false;
use_planarRR=false;
if(strcmp(sys_info.id,'cart_pole'))
    use_cart_pole=true;
    path_integral_params = cart_pole_PI_params;
elseif(strcmp(sys_info.id,'quadrotor'))
    path_integral_params = quadrotor_PI_params;
    use_quadrotor=true;
elseif(strcmp(sys_info.id,'planarRR'))
    path_integral_params = planarRR_PI_params;
    use_planarRR=true;
elseif(strcmp(sys_info.id,'non_linear'))
    path_integral_params = nonlinear_PI_params;
    use_non_linear=true;
elseif(strcmp(sys_info.id,'linear'))
    path_integral_params = linear_PI_params;
elseif(strcmp(sys_info.id,'pendulum'))
    path_integral_params = pendulum_PI_params;
elseif(strcmp(sys_info.id,'duffing'))
    path_integral_params = duffing_PI_params;
elseif(strcmp(sys_info.id,'nonlinear_saddle'))
    path_integral_params = nonlinear_saddle_PI_params;
else
    path_integral_params = PI_params;
end
% check for reverse time
use_reverse = path_integral_params.unstable_reverse;
if(use_reverse)
    eig_vals = -diag(D);
else
    eig_vals = diag(D);
end

%% open loop simualtion
t_start = 0;
dt_sim  = path_integral_params.dt_sim;
t_end   = path_integral_params.t_end;
x_op    = x_local;

% Predefine a linear dynamics handle that matches rk4's signature
dyn_lin = @(x,u,sys) sys.A*x;

Tout = []; Xout = [];
for t_sim = t_start:dt_sim:t_end
    % zero-control (your existing logic)
    if use_quadrotor
        u_zero = sys_info.u_zero - sys_info.k_poles*x_op;
    else
        u_zero = sys_info.u_zero;
    end

    % Full nonlinear step
    x_next_full = rk4(dynamics, dt_sim, x_op, u_zero, use_reverse, sys_info);

    % Linear-only step with same integrator
    x_next_lin  = rk4(dyn_lin, dt_sim, x_op, u_zero, use_reverse, sys_info);

    % Nonlinear contribution in the state (full minus linear)
    x_next_nl = x_next_full - x_next_lin;

    % Update the "operating" state with the full step
    x_op = x_next_full;

    % Logs
    Tout = [Tout; t_sim];
    Xout = [Xout; x_next_nl.'];   % store the nonlinear part of the state
end

% shift eqb point
Xout   = Xout - x_eqb';

%% compute nonlinear part of eigfun
integrand_convergence = cell(n_dim);
solution_convergence  = cell(n_dim);

for i = 1:n_dim
    % get eigval and eigvec
    lambda  = eig_vals(i);
    w       = W(:,i);

    % compute path integral
    integrand = exp(-Tout*lambda).*(w'*Xout')';
    phi_nonlinear{i} = trapz(Tout,integrand,1);
    phi_linear{i} = w'*x_local;
    phi{i} = phi_linear{i}  + phi_nonlinear{i};

    % check for convergence
    sol_conv = (w'*Xout')';
    integrand_convergence{i} = integrand(end);
    solution_convergence{i}  = sol_conv(end);
end

%% Loop through each element in phi and assign it to phi_forward.phi
for i = 1:n_dim
    eigen_fn.phi(i)           = phi{i};
    eigen_fn.phi_linear(i)    = phi_linear{i};
    eigen_fn.phi_nonlinear(i) = phi_nonlinear{i};
    eigen_fn.integrand(i)     = integrand_convergence{i};
    eigen_fn.sol_conv(i)      = solution_convergence{i};
end