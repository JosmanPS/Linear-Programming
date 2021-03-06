

% Definiremos un problema f�cil de ejemplo para probar el algoritmo
c = [-5; 1; 8; -4; 3];
A = [0 6 0 -1 3; 1 3 0 0 2; 0 -1 1 0 1];
b = [18; 24; 4];


% ************* SIMPLEX ********************************************

% ************* PASO 1 *********************************************

% Buscamos un v�rtice inicial en F y tomamos B la submatriz de A con 
% columas correspondientes a los valores del v�rtice distintos de 0
[x, I1] = LP_verts(A,b);
B = A(:,I1);
aux = zeros(1,length(x));
aux(I1) = 1;
I2 = 1:length(x);
I2 = I2(aux == 0);


while 1
% ************* PASO 2 *********************************************

% Resolvemos el sistema lineal B' * PI = C_B
PI = linsolve(B', c(I1,1));

% Calculamos los costos reducidos
Cj = c(I2) -  A(:,I2)' * PI;


% ************* PASO 3 *********************************************

if(norm(x((m+1):(m+n)) == 0)
    disp('La soluci�n es:')
    disp(x)
    disp(' El valor de la funci�n objetivo en x es:')
    f = c' * x;
    disp(f)
    return
end


% ************* PASO 4 *********************************************
jopt = I2(Cj < 0);
jopt = jopt(1);


% ************* PASO 5 *********************************************
z = linsolve(B,A(:,jopt));
if(sum(z <= 0) == length(z))
    disp('La funci�n objetivo no est� acotada inferiormente en F');
    return
end


% ************* PASO 6 *********************************************
Xj = x(I1) ./ z;
[Xjp, p] = min(Xj(Xj > 0));
aux = 1:length(Xj);
jp = aux(Xjp == Xj);
jp = I1(jp);
Xjp = Xjp(1);

% ************* PASO 7 *********************************************

% Actualizamos
theta = Xjp;
x(jopt) = theta;
x(I1) = x(I1) - theta * z; 
I1 = I1(I1 ~= jp);
I1 = [I1, jopt];
I1 = sort(I1);
B = A(:,I1);
aux = zeros(1,length(x));
aux(I1) = 1;
I2 = 1:length(x);
I2 = I2(aux == 0);

end




