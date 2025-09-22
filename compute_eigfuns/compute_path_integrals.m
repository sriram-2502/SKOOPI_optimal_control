function phi = compute_path_integrals(x_op, dynamics, sys_info)
    
    show_wait_bar    = false;
    show_diagnostics = false;
    
    % parse info
    D     = sys_info.eig_vals;
    W     = sys_info.eig_vectors;
    x_eqb = sys_info.x_eqb;

    %% setup grid
    [local_grid, local_axes] = setup_local_grid(x_op);
    num_elements = numel(local_grid{1}); % Number of grid points
    n_dim = length(x_op); % Number of dimensions
    
    % Flatten each grid into a column and concatenate into a matrix of points
    grid_points = cellfun(@(grid) grid(:), local_grid, 'UniformOutput', false);
    grid_points = [grid_points{:}];  % Concatenate into a single matrix
    
    %% Local path integral computation
    % Initialize arrays to store computed values
    phi_complete  = nan(num_elements, n_dim);
    phi_linear    = nan(num_elements, n_dim);
    phi_nonlinear = nan(num_elements, n_dim); 
    phi_integrand = nan(num_elements, n_dim); 
    phi_x_op           = nan(1, n_dim);  
    phi_integrand_x_op = nan(1, n_dim); 
    
    % load a progress bar
    if(show_wait_bar)
        w_bar = waitbar(0,'1','Name','computing path integrals ...', 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    end

    for idx = 1:num_elements
        if(show_wait_bar)
            waitbar(idx/num_elements,w_bar,sprintf(string(idx)+'/'+string(num_elements)))
        end
        x_local = grid_points(idx,:)';

        % compute path integral around the operating point
        if(all((diag(D)) ~= 0))
            phi                          = compute_eigen_fn(x_local, x_eqb, dynamics, D, W, sys_info);
            phi_complete(idx, 1:n_dim)   = phi.phi(:)';
            phi_linear(idx, 1:n_dim)     = phi.phi_linear(:)';
            phi_nonlinear(idx, 1:n_dim)  = phi.phi_nonlinear(:)';
            phi_integrand(idx, 1:n_dim)  = phi.integrand(:)';

            % Extract values for current operating point
            if(norm(x_local - x_op) <= 1e-3)
                if(show_diagnostics)
                    disp('----- computing eig_fun at x_op -----')
                end
                phi_x_op(1:n_dim)           = phi.phi(:)';
                phi_integrand_x_op(1:n_dim) = phi_integrand(idx, 1:n_dim);
            end
        else
            disp('!!! system has zero eigen values. Cannot use path integrals !!!')
        end
    end

    if(show_wait_bar)
    % delete progress bar
        if ishandle(w_bar)  % Ensure the waitbar handle is still valid
            delete(w_bar);
        end
    end

    %% Reshape phi.phi into cells
    phi_complete_cell   = cell(n_dim, 1);
    phi_linear_cell     = cell(n_dim, 1);
    phi_nonlinear_cell  = cell(n_dim, 1);
    phi_integrand_cell  = cell(n_dim, 1);
    
    for i = 1:n_dim
        % Get the size of the grid in the first dimension (same for each dimension)
        grid_size = size(local_grid{1}); % Assuming all grids have the same size
        
        % Reshape phi values for each dimension into the grid shape
        phi_complete_cell{i}    = reshape(phi_complete(:, i), grid_size);
        phi_linear_cell {i}     = reshape(phi_linear(:, i), grid_size);
        phi_nonlinear_cell {i}  = reshape(phi_nonlinear(:, i), grid_size);
        phi_integrand_cell {i}  = reshape(phi_integrand(:, i), grid_size);
    end
    
    %% store everything in a struct
    phi = struct();

    % store the axis values
    phi.axis = local_axes;
    phi.grid = local_grid;

    % store the reshaped phi values as cells for each dimension
    phi.phi             = phi_complete_cell;
    phi.phi_linear      = phi_linear_cell;
    phi.phi_nonlinear   = phi_nonlinear_cell;
    phi.phi_integrand   = phi_integrand_cell;

    % store eigfun values at the operating point
    phi.phi_x_op           = phi_x_op;
    phi.phi_integrand_x_op = phi_integrand_x_op;
end
