function grad_phi_x_op = compute_gradients(phi)
% Compute the gradient of phi at the operating point x_op for each phi function.
%
% Inputs:
% x_op          : Operating point (1 x n_dim array)
% phi           : Struct containing phi with values for each function and axis grid
%
% Outputs:
% grad_phi_x_op : n_dim x n_dim matrix, each column containing the gradient of
%                 each phi function at x_op in each dimension
% NOTE: Transponse grad_phi_x_op to match with W' (linear eig vec
% transpose)

    n_dim = length(phi.phi_x_op);               % Number of dimensions
    grad_phi_x_op = nan(n_dim, n_dim);  % Matrix to store gradients for each phi function

    % Loop over each phi function
    for func_idx = 1:n_dim
        % Get the phi function for the current index and initialize gradient vector
        phi_func = phi.phi{func_idx};
        grad_phi = nan(n_dim, 1);  % Initialize gradient vector for this phi function

        % Find the center index in each dimension for the central point
        idx_center = zeros(1, n_dim);
        for dim = 1:n_dim
            idx_center(dim) = ceil(size(phi_func, dim) / 2);  % Midpoint index in each dimension
        end

        % Loop over each dimension to compute the central difference gradient
        for dim = 1:n_dim
            % Determine the step size in the current dimension
            grid = phi.axis{dim};       % Axis values for this dimension
            delta = grid(2) - grid(1);  % Assuming uniform grid spacing

            % Initialize index arrays for forward and backward points
            idx_forward = idx_center;
            idx_backward = idx_center;

            % Move forward and backward along the current dimension
            idx_forward(dim) = idx_center(dim) + 1;
            idx_backward(dim) = idx_center(dim) - 1;

            % Extract scalar values at the forward and backward points in the current dimension
            phi_forward = get_nd_value(phi_func, idx_forward);
            phi_backward = get_nd_value(phi_func, idx_backward);

            % Compute the central difference gradient for this dimension
            grad_phi(dim) = (phi_forward - phi_backward) / (2 * delta);
        end

        % Store the gradient vector for the current phi function
        grad_phi_x_op(:, func_idx) = grad_phi;
    end
    %Transponse grad_phi_x_op to match with W'
    grad_phi_x_op = grad_phi_x_op';
end

function val = get_nd_value(array, indices)
% Helper function to index into an n-dimensional array with a dynamic set of indices
    s.type = '()';
    s.subs = num2cell(indices);  % Convert indices to cell for dynamic indexing
    val = subsref(array, s);     % Use subsref to get the scalar value at specified indices
end
