%script file: vertices3.m
% se encuentran los vértices del poliedro definido por
%  A*x = b ,  x >= 0
%
% donde
%     [0 6 0 -1 3]       [ 18]
% A = [1 3 0  0 2]   b = [ 24]
%     [0 -1 1 0 1]       [ 4 ]
%
% Programación Lineal
%  ITAM
%  19 de septiembre de 2014
% Dr. Zeferino Parada
%------------------------------------------
% matriz y lado derecho
A =[ 0 6 0 -1 3;  1 3 0 0 2;  0 -1 1 0 1];
b = [18 24 4]';
%c = [ -5 1 8 -4 3]';

%-----------------------------------------------------------
% Casos de subconjuntos de tres elementos de un conjunto de 5 elementos
W =[1 2 3; 1 2 4; 1 2 5; 1 3 4; 1 3 5; 1 4 5;  2 3 4; 2 3 5; 2 4 5; 3 4 5];
[nr, nc] = size(W); % número de renglones (nr) y columnas(nc) en W
%----------------------------------------------------------
          
V = [];   % se recolectan los vértices del poliedro.

for j = 1:nr
    v = W(j,:);              
    A1 = A(:,v);         % se escogen tres columnas de A/ A1 es 3x3
    d = det(A1);         % el determinante de A1
    if(d ~= 0)           % el sistema lineal A1*xB = b tiene solución única
        xB = A1\b;       % solución del sistema
        if(xB >= 0)      % se verifica que la solución es no negativa
            x = zeros(5,1);
            x(v) = xB;   % vértice del poliedro
            V = [V x];   % se recolecta el nuevo vértice
        end
    end
end

%-------------------------------------
% Muestra los vértices
disp(' Vértices de A*x = b, x >= 0')
A
disp('-------------------------------------')
b
disp('--------------------------------------')
disp('Vértices')

V


%-------------------------------------
% Una iteracion de Simplex
c = [ -5 1 8 -4 3]';
x = V(:,2);
varB = [1 2 5];
B = A(:, varB);
cB = c(varB);

% Paso 2
vpi = B' \ cB;       % Resolvemos el sistema lineal.
% Costos reducidos
j = 3;
c3 = c(j) - vpi' * A(:,j)  % > 0
j = 4;
c4 = c(j) - vpi' * A(:,j)  % < 0

xi = B \ A(:,j)
xi = -xi;
eta4 = zeros(5,1);
eta4(varB) = xi; eta4(4) = 1;

theta = - (x(1) / eta4(1));

xn = x + theta * eta4

