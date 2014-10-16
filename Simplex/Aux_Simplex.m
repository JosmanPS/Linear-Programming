
% ************* SIMPLEX ********************************************
function [x, f, It] = Aux_Simplex(c,A,b)

% ************* PASO 1 *********************************************

% Buscamos un vértice inicial en F y tomamos B la submatriz de A con 
% columas correspondientes a los valores del vértice distintos de 0
[m,n] = size(A);
n = n - m;
x = zeros(n+m,1);
% Vértice inicial
x((n+1):(n+m)) = b;
I1 = (n+1):(n+m);
B = A(:,I1);
aux = zeros(1,length(x));
aux(I1) = 1;
I2 = 1:length(x);
I2 = I2(aux == 0);
It = 1; % Iteraciones

while 1
% ************* PASO 2 *********************************************

% Resolvemos el sistema lineal B' * PI = C_B
PI = linsolve(B',c(I1,1));

% Calculamos los costos reducidos
Cj = c(I2) -  A(:,I2)' * PI;


% ************* PASO 3 *********************************************

if(norm(x((m+1):(m+n))) == 0)
    return
end

stop = 0;
i = 1;
while stop == 0
% ************* PASO 4 *********************************************
jopt = I2(Cj < 0);
if i > numel(jopt)
    return
end
jopt = jopt(i);


% ************* PASO 5 *********************************************
z = linsolve(B, A(:,jopt));


% ************* PASO 6 *********************************************
Xj = x(I1) ./ z;
[Xjp, p] = min(Xj(z > 0)); %%%%%%%%%%%%%%%%%
aux = 1:length(Xj);
jp = aux(Xjp == Xj);
jp = I1(jp);
jp = jp(1);
Xjp = Xjp(1);

    if Xjp ~= 0
% ************* PASO 7 *********************************************

        % Actualizamos
        theta = Xjp;

        x(jopt) = theta; x(I1) = x(I1) - theta * z;

        I1 = I1(I1 ~= jp); I1 = [I1, jopt]; I1 = sort(I1);

        B = A(:,I1); %%%%%%%%%%%%

        aux = zeros(1,length(x)); aux(I1) = 1;

        I2 = 1:length(x); I2 = I2(aux == 0);

        It = It + 1; %%%%%%%

        stop = 1;

    else
        i = i + 1;

    end

end

end

end



