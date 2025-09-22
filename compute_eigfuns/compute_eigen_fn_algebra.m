function eigen_fn_algebra = compute_eigen_fn_algebra(x_local, x_eqb, dynamics, D,D_stable, W, sys_info, k)
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

scale = 1;
threshhold = 1;
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

% setup for algebra
eig_vector_scale = path_integral_params.eigen_vector_scale;
[lambda_unstable, idx_unstable] = max(diag(D)); % find max positive eig val
idx_stable    = ismember(diag(D), D_stable);
lambda_stable = D_stable;   % negative eig val
w_stable      = W(:,idx_stable)./eig_vector_scale;    % eig vector ~ negative eig val
w_unstable    = W(:,idx_unstable)./eig_vector_scale;    % eig vector ~ positive eig val
Q             = w_stable*w_unstable';
lambda        = lambda_stable + k*lambda_unstable; % make sure lambda is always positive

%% open loop simualtion
t_start = 0;
dt_sim  = path_integral_params.dt_sim;
t_end   = path_integral_params.t_end;
Xout    = x_local';
Tout    = 0;
x_op    = x_local;
Xout_linear = x_local';
Xout_full   = x_local';

for t_sim = t_start:dt_sim:t_end
    % forward simulate using rk4 with no control
    if(use_quadrotor)
        u_zero = sys_info.u_zero - sys_info.k_poles*x_op;
    else
        u_zero = sys_info.u_zero;
    end
    x_next_full   = euler(dynamics,dt_sim,x_op,u_zero,use_reverse,sys_info);
    x_next_linear = euler(dynamics_linearized,dt_sim,x_op,u_zero,use_reverse,sys_info);

    % get nonlinear part only
    x_next = (x_next_full-x_op) / dt_sim - A*x_op;

    % update
    x_op = x_next_full;

    % logs
    Tout  = [Tout;t_sim];
    Xout  = [Xout;x_next'];
    Xout_full   = [Xout_full; x_next_full'];
    Xout_linear = [Xout_linear; x_next_linear'];

end

% shift eqb point
Xout        = Xout - x_eqb';
Xout_full   = Xout_full - x_eqb';

% scale
Xout = scale*Xout;
Xout_full = scale*Xout_full;
lambda = lambda/scale;

%% compute nonlinear part of eigfun
% compute path integral
integrand = exp(-Tout*lambda).*((w_unstable'*Xout_full').^(k-1)*(Xout_full)*(Q+k*Q')*(Xout'))';
phi_nonlinear = trapz(Tout,integrand,1);
phi_linear = (w_stable'*x_local)*(w_unstable'*x_local)^k;
phi = phi_linear + phi_nonlinear;

% check for convergence
sol_conv = (w_unstable'*Xout_full').^(k-1)*(Xout_full)*(Q+k*Q')*(Xout');
integrand_convergence = integrand(end);
solution_convergence  = sol_conv(end);


% Loop through each element in phi and assign it to phi_forward.phi
eigen_fn_algebra.phi           = phi./threshhold;
eigen_fn_algebra.phi_linear    = phi_linear;
eigen_fn_algebra.phi_nonlinear = phi_nonlinear;
eigen_fn_algebra.integrand     = integrand_convergence;
eigen_fn_algebra.sol_conv      = solution_convergence;
