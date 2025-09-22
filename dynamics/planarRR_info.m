function sys_info = planarRR_info(t, x, u, robot_params, sys_params)

% parse robot info
m1 = robot_params.m1 ; m2 = robot_params.m2;
l1 = robot_params.l1 ; l2 = robot_params.l2;
g = robot_params.g;
I1 = robot_params.I1 ;I2 = robot_params.I2; 
lc1 = l1/2; lc2 = l2/2;
D = robot_params.D;

% parse states
q1 = x(1,:); q2 = x(2,:); dq1 = x(3,:); dq2 = x(4,:);

% define eqb point as upright position
x_eq = [0; 0; 0; 0];
u_eq = [0; 0];

% define parameters
m11 = I1 + I2 + m2*l1^2 + 2*m2*l1*lc2.*cos(q2);
m12 = I2 + m2*l1*lc2.*cos(q2);
m22 = I2*ones(1,size(x,2));
sys_info.M = [m11 m12; m12 m22];
sys_info.C = [-2*m2*l1*lc2.*sin(q2).*dq2.*dq1-m2*l1*lc2.*sin(q2).*dq2.*dq2;
     m2*l1*lc2.*sin(q2).*dq1.*dq1+zeros(1, size(x,2)).*dq2];
sys_info.D = D;

% G is negative for q1 = 0  about vertical
% G is positive for q1 = pi about vertical
sys_info.G = [m1*g*lc1.*sin(q1) + m2*g.*(l1.*sin(q1)+lc2.*sin(q1+q2));
     m2*g*lc2.*sin(q1+q2)];

% for acrobot
%B = [zeros(1,size(x,2)); ones(1,size(x,2))];
%ddq = M\(B.*u-C-G);

B = eye(2);
dq = [dq1; dq2]; 
ddq = sys_info.M\(B*u-sys_info.C-sys_info.G-D*dq);

% sys_info.f = [dq;ddq];
% f_x = [dq;sys_info.M\(-sys_info.C-sys_info.G-D*dq)];
sigma = 0.0; 
sys_info.tau_d = sigma * randn(2,1);
f_x = [dq;sys_info.M\(-sys_info.C-sys_info.G-D*dq -sys_info.tau_d)];
% sys_info.tau_d

g_x = [zeros(2,2);sys_info.M\(B)];

if isa(x, 'sym') || isa(x, 'symfun')
    f = matlabFunction(f_x, 'vars', {x});
    g = matlabFunction(g_x, 'vars', {x});
    
    % Linearize around x_eq = [pi; 0; 0; 0], u_eq = [0; 0]
    % Compute Jacobians
    A_sym = jacobian(f_x, x);
    
    % Substitute equilibrium point
    sys_info.A = double(subs(A_sym, [x; u], [x_eq; u_eq]));
    sys_info.B = double(subs(g_x, [x; u], [x_eq; u_eq]));

    [~,Lambda,W] = eig(sys_info.A);
    sys_info.eig_vals       = Lambda;
    sys_info.eig_vectors    = W;
    sys_info.dynamics_f     = f;
    sys_info.dynamics_g     = g;
else
    f = @(x) f_x;  % Returns the fixed numeric value of f_x
    g = @(x) g_x; 
    sys_info.dynamics_f = f;
    sys_info.dynamics_g = g;
end

sys_info.x_eqb          = x_eq;
sys_info.use_stable     = sys_params.use_stable;
sys_info.use_unstable   = sys_params.use_unstable; 
sys_info.id             = 'planarRR';
sys_info.n_states       = 4;
sys_info.n_ctrl         = 2;
sys_info.u_zero         = u_eq;

end