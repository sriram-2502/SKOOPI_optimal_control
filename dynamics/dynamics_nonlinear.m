function [dxdt,fx,gx] = dynamics_test_nonlinear(x,u, sys_info)

dxdt = sys_info.dynamics_f(x) + sys_info.dynamics_g(x)*u;
fx = sys_info.dynamics_f;
gx = sys_info.dynamics_g;

end