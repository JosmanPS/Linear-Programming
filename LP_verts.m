function verts = LP_verts(A, b)

% Input:
% This function gets the matrix A and the vector b that indicates the
% restrictions of a linear programming problem.

% Output:
% The function returns a matrix "verts". Each column of the matrix
% represents a vertex of the restriction area.



% Get the dimensions of A.
[l1, l2] = size(A);

% We need b to be a vector of (n * 1) not (1 * n).
if( size(b,2) > 1)
    b = b'
end

% Define the matrix verts.
verts = zeros(l2,1);

% A counter to indicate the column of verts
ind = 1;

% Get the different square sub-matrix of A.
for i = 1:(l2-2)
    for j = (i+1):(l2-1)
        for k = (j+1):l2
           
            % Square sub-matrix of A.
            matrix = A(:, [i,j,k]);
           
            % Solve the system
            x = inv(matrix) * b;
           
            % If X >= 0, then X is a vertex.
            if(sum(x >= 0) == length(x))
               
                % Write the vertex in the matrix.
                verts([i,j,k], ind) = x;
                ind = ind + 1;
               
            end       
        end
    end
end
