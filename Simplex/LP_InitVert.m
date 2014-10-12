function [A, v] = LP_InitVert(A,b)
A = full(A);
b = full(b);

% Cambiamos los valores negativos
negs = b < 0;
b(negs) = -b(negs);
A(:,negs) = -A(:,negs);

% Construimos el problema artificial
[m,n] = size(A);
c = [zeros(n,1);ones(m,1)];
B = [A eye(m)];

% Resolvemos el problema
v = Aux_Simplex(c,B,b);
v((n+1):(n+m)) = [];


