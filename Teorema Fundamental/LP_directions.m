function d = LP_directions(A)

% Input:
% This function gets the matrix A that indicates the restrictions of a 
% linear programming problem.

% Output:
% The function returns a matrix "d". Each column of the matrix
% represents a direction of the restriction area.

% Take a basis of the nullspace of A.
basis = null(A);
d = [];

% Select the non-negative vectors.
for j = 1:size(basis,2)
    if(sum(basis(:,j) >= 0) == size(basis,1))
        d = [d, basis(:,j)];
    end
end

if(norm(d) == 0)
    disp('There are not directions in A')
end

