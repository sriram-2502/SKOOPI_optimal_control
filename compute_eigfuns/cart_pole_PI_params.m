function path_integral_params = cart_pole_PI_params()

path_integral_params.dt_sim = 0.01; % sim time step
path_integral_params.t_end = 2; % integration time (best is t=2 for algebra)
path_integral_params.unstable_reverse = false; % use reverse flow in unstable
path_integral_params.stable_reverse = false; % use reverse flow in unstable
path_integral_params.eigen_vector_scale = 1; % scaling for algebra methods
end