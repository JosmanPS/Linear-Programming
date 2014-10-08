function [] = LP_TeoFund(A,b,c)

% *************************************************************************
% Input:
% This function gets the matrix A and the vector b that indicates the
% restrictions of a linear programming problem. The vector of costs c.

% Output:
% Return which case of the Fundamental Theorem of Linear Programming holds.
% *************************************************************************


% We need b to be a vector of (n * 1) not (1 * n).
if( size(b,2) >= size(b,1))
    b = b';
end

% First case
% F = { x in R^n : Ax = b, x >= 0} = {}
% If there are not vertex in F, then F = {}.
verts = LP_verts(A,b);

if(size(verts,1) == 0)
    disp('No solution: F = { x in R^n : Ax = b, x >= 0} = {}.');
    return
end

% Second case
% F ~= {}
% We look for the directions of F.

% If D = {directions of F} = {}, then third case holds
D = LP_directions(A);
if (size(D,1) == 0)
    disp('The problem has a solution in a vertex of F.');
    return
end

% If D ~= {}, then we look for d in D, such that c' * d < 0. If it holds,
% then there is no solution.

% We need c to be a vector of (1 * n) not (n * 1).
if( size(c,2) <= size(c,1))
    c = c';
end

values = c * D;
values = c < 0;

if(sum(values) == 0)
    disp('The problem has a solution in a vertex of F');
    return
else
    disp('No solution: Exists d in D, such that t(c) * D < 0, F is not bounded.');
    return
end





