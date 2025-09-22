function result = replace_complex(M)
% replace_complex Replaces complex conjugate pairs in a matrix.
%
%   result = replace_complex_matrix(M)
%   For each column of the matrix M, the function replaces complex conjugate pairs
%   a ± ib with (a, b). The size of the output matrix matches the input.
%
%   Input:
%       M - A matrix containing real numbers and complex conjugate pairs.
%
%   Output:
%       result - A matrix of the same size as M, with processed values.

    [rows, cols] = size(M);  % Get the size of the matrix
    result = zeros(rows, cols);  % Initialize output matrix with zeros
    
    for c = 1:cols
        result(:, c) = process_vector(M(:, c));  % Process each column
    end
end

function result = process_vector(v)
% PROCESS_VECTOR Processes a vector to replace complex conjugate pairs
    result = zeros(size(v));  % Initialize result with zeros
    i = 1;
    
    while i <= length(v)
        if isreal(v(i))
            % Keep real numbers as they are
            result(i) = v(i);
            i = i + 1;
        else
            % For complex conjugate pairs, replace first with real part, second with imaginary part
            a = real(v(i));
            b = imag(v(i));
            result(i) = a;      % Replace the first occurrence with the real part (a)
            result(i + 1) = b;  % Replace the second occurrence with the imaginary part (b)
            i = i + 2; % Skip the conjugate pair
        end
    end
end
