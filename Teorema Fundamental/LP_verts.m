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
if( size(b,2) > size(b,1))
    b = b';
end

% Define the matrix verts.
verts = [];

% Different combinations for square matrix.
l1 = nchoosek(1:l2, l1);

for i = 1:size(l1,1)
    
    % Square sub-matrix of A.
    matrix = A(:,l1(i,:));
    
    if(det(matrix) ~= 0)
        % Solve the system.
        x = zeros(l2, 1);
        x(l1(i,:),1) = linsolve(matrix, b);
        
        % If X >= 0, then X is a vertex.
        if(sum(x >= 0) == length(x))
            
            % Write the vertex in the matrix.
            verts = [verts, x];
            
        end
    end
end

       
