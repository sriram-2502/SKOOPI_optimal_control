function path_integral_params = pendulum_PI_params()

path_integral_params.dt_sim = 0.01; % sim time step
path_integral_params.t_end = 10; % end time for open loop simualtion and final time of path integral
path_integral_params.unstable_reverse = false; % use reverse flow in unstable
path_integral_params.stable_reverse = false; % use reverse flow in unstable
path_integral_params.eigen_vector_scale = 1000; % scaling for algebra methods
end