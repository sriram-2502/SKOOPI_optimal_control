function eigen_fn = compute_eigen_fn_reverse(x_local, x_eqb, dynamics, D, W, sys_info)
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

% parse path integral setup
use_cart_pole=false;
use_quadrotor=false;
if(strcmp(sys_info.id,'cart_pole'))
    use_cart_pole=true;
    path_integral_params = cart_pole_PI_params;
elseif(strcmp(sys_info.id,'quadrotor'))
    path_integral_params = quadrotor_PI_params;
    use_quadrotor=true;
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
eig_vector_scale = path_integral_params.eigen_vector_scale;
W = W./eig_vector_scale;

% check for reverse time
use_reverse = path_integral_params.unstable_reverse;
if(use_reverse)
    eig_vals = -diag(D);
else
    eig_vals = diag(D);
end

%% open loop simualtion
if(use_cart_pole||use_quadrotor||use_non_linear)
    t_start = 0;
    dt_sim  = path_integral_params.dt_sim;
    t_end   = 2;
    Xout    = x_local';
    x_op    = x_local;
    Tout    = 0;

    % forward simulate using rk4 with no control
    for t_sim = t_start:dt_sim:t_end
        if(use_quadrotor)
            u_zero = sys_info.u_zero - sys_info.k_poles*x_op;
        else
            u_zero = sys_info.u_zero;
        end
        x_next_full = euler(dynamics,dt_sim,x_op,u_zero,use_reverse,sys_info);
        % get nonlinear part only
        x_next = (x_next_full-x_op) / dt_sim - A*x_op;
        % update
        x_op = x_next_full;
        % logs
        Tout  = [Tout;t_sim];
        Xout  = [Xout;x_next'];
    end
    % shift eqb point
    Xout   = Xout - x_eqb';
    
    % reverse simulate using rk4 with no control
    x_op_reverse = x_local;
    Tout_reverse = 0;
    Xout_reverse = x_local';
    t_end_reverse = 2; 
    for t_sim = t_start:dt_sim:t_end_reverse
        if(use_quadrotor)
            u_zero = sys_info.u_zero - sys_info.k_poles*x_op_reverse;
        else
            u_zero = sys_info.u_zero;
        end
        x_next_full_reverse = euler(dynamics,dt_sim,x_op_reverse,u_zero,true,sys_info);
        % get nonlinear part only
        x_next_reverse = -(x_next_full_reverse-x_op_reverse)/dt_sim; 
        x_next_reverse = x_next_reverse + A*x_op_reverse;
        x_next_reverse = x_next_full_reverse;
        % update
        x_op_reverse = x_next_full_reverse;
        % logs
        Tout_reverse  = [Tout_reverse;t_sim];
        Xout_reverse  = [Xout_reverse;x_next_reverse'];
    end
    % shift eqb point
    Xout_reverse   = Xout_reverse - x_eqb';

else
    t_start = 0;
    dt_sim  = path_integral_params.dt_sim;
    t_end   = path_integral_params.t_end_unstable;
    Xout    = x_local';
    x_op    = x_local;
    Tout    = t_start;

    for t_sim = t_start+dt_sim:dt_sim:t_end
        % forward simulate using rk4 with no control
        u_zero = zeros(1,sys_info.n_ctrl);
        x_next_full = euler(dynamics,dt_sim,x_op,u_zero,use_reverse,sys_info);

        % get nonlinear part only
        x_next = (x_next_full-x_op) / dt_sim - A*x_op;

        % update
        x_op = x_next_full;

        % logs
        Tout  = [Tout;t_sim];
        Xout  = [Xout;x_next'];
    end

    % shift eqb point
    Xout   = Xout - x_eqb';
end

%% compute nonlinear part of eigfun
integrand_convergence = cell(n_dim);
solution_convergence  = cell(n_dim);

for i = 1:n_dim
    % get eigval and eigvec
    lambda  = eig_vals(i);
    w       = W(:,i);

    if(lambda>0)
        lambda_val = lambda;
        Xout_val = Xout;
        Tout_val = Tout;
    else
        lambda_val = -lambda;
        Xout_val = Xout_reverse;
        Tout_val = Tout_reverse;
    end

    % compute path integral
    integrand = exp(-Tout_val*lambda_val).*(w'*Xout_val')';
    phi_nonlinear{i} = trapz(Tout_val,integrand,1);
    phi_linear{i} = w'*x_local;
    if(lambda>0)
        phi{i} = phi_linear{i}  + phi_nonlinear{i};
    else
        phi{i} = phi_linear{i}  - phi_nonlinear{i};
    end

    % check for convergence
    sol_conv = (w'*Xout_val')';
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