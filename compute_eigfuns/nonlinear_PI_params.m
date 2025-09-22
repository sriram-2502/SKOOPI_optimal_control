function path_integral_params = nonlinear_PI_params()

path_integral_params.dt_sim = 0.05; % sim time step
path_integral_params.t_end = 1; % end time for open loop simualtion and final time of path integral
path_integral_params.t_end_unstable = 1; % use for different integration time for unstable eig fun
path_integral_params.unstable_reverse = true; % use reverse flow in unstable
path_integral_params.stable_reverse = false; % use reverse flow in unstable
path_integral_params.eigen_vector_scale = 1; % scaling for eigfns

end