%
% Simplex (con m?todo de la Gran-M)
%
% Esta es una implementacion del algoritmo simplex utilizando tableus. Primero 
% se pide el numero de variables y la cantidad de restricciones. Despues se pide 
% si es un problema de maximizacion o un problema de minimizacion, cada coeficiente
% de acuerdo a lo ingresado anteriormente y los sentidos de las desigualdades.
%
% Una vez que se tiene toda la informacion, primero se crea el tableu inicial, 
% transformado a un problema de minimizacion si es necesario, y se busca 
% factibilidad primal, despues se busca una identidad inicial, y de ahi en adelante 
% se busca la solucion iterando normalmente. Si alguna de las etapas indica que el 
% problema no es factible se notifica al usuario, y si encuentra soluciones 
% multiples las caracteriza.
%
% Como resultados finales indica si el algoritmo termino exitosamente, el valor
% de M utilizado en el metodo, el numero de iteraciones realizadas, el numero de
% veces que se utilizo la Regla de Bland (se utiliza cada vez que hay presencia
% de ceros en B^(-1)*b), la estructura del vector x, los indices de la base
% final, el vector x, y el valor de la funcion objetivo evaluada en x.
%
% Omar Trejo Navarro
%
% 10 de abril de 2014
%

function simplex()
    format RAT;
    fprintf('\n=========================================');
    fprintf('\n==   Implementacion Simplex (Gran-M)   ==');
    fprintf('\n=========================================');
    [tableu, tipo, ind_vb, n, n_holguras, n_artificiales, M] = formarTableu();    
    ind_artificiales = (n + n_holguras + 1):(n + n_holguras + n_artificiales);
    fprintf('\n[+] Tableu 1:\n\n');
    disp(tableu);
    n_tableu = 2;
    terminado = 0;
    bland = 0;
    while any(tableu(end,1:(end - 1)) > 0) && not(terminado)
        if any(tableu(1:(end - 1), end) == 0)
            bland = bland + 1;
            [tableu, ind_vb, terminado] = iterarBland(tableu, ind_vb);
        else
            [tableu, ind_vb, terminado] = iterar(tableu, ind_vb);
        end
        fprintf('\n[+] Tableu %d:\n\n', n_tableu);
        disp(tableu);
        n_tableu = n_tableu + 1;
    end
    %
    % Revisar factibilidad
    %
    sigue_factible = revisarFactibilidad(...
        tableu, ind_vb, n_artificiales, ind_artificiales, M);
    if not(sigue_factible)
        fprintf('\n[!] El problema es no-factible porque los coeficientes con');
        fprintf('\n    M son todos no-positivos y las variables artificiales');
        fprintf('\n    siguen siendo positivas.\n');
        terminado = 1;
    else
        for i = 1:n_artificiales
            if not(isempty(ind_vb(ind_vb == ind_artificiales(i)))) && not(terminado)
                fprintf('\n[!] El problema es no-factible porque en la base optima');
                fprintf('\n    hay presencia de variables artificiales.\n');
                terminado = 1;
            end
        end
    end
    %
    % Revisar soluciones multiples
    %
    revisarMultiples(tableu, ind_vb);
    %
    % Imprimir solucion
    %
    solucion(tableu, tipo, ind_vb, n_tableu, n, ...
        n_holguras, n_artificiales, M, terminado, bland);
end

function [tableu, tipo, ind_vb, n, n_holguras, n_artificiales, M] = formarTableu()
    %
    % Datos de restricciones
    %
    fprintf('\n\n[+] Tipo de optimizacion:\n');
    fprintf('\n\tMaximizacion: max');
    fprintf('\n\tMinimizacion: min\n');
    tipo = input('\n[-] Tipo: ', 's');
    while not(strcmp(tipo, 'max') || strcmp(tipo, 'min'))
        fprintf('\n[!] Error: solo se puede ingresar max o min.\n');
        tipo = input('\n[-] Tipo: ', 's');
    end
    n = input('[-] Cantidad de variables: ');
    m = input('[-] Cantidad de restricciones: ');
    tableu = zeros((m+1), (n+1));
    n_holguras = 0;
    for i = 1:m
        fprintf('\n[+] Restriccion %d:\n\n', i)
        for j = 1:n
            fprintf('\t[-] Coeficiente de la variable %d: ', j);
            tableu(i,j) = input('');
        end
        fprintf('\t[-] Signo de la restriccion (=, <=, >=): ')
        s = input('', 's');
        while not(strcmp(s, '=') || strcmp(s, '>=') || strcmp(s, '<='))
            fprintf('\n[!] Error: solo se puede ingresar =, <= o >=.\n');
            fprintf('\n\t[-] Signo de la restriccion (=, <=, >=): ')
            s = input('', 's');
            fprintf('\n');
        end
        if strcmp(s, '<=')
            % Agregar holgura positiva.
            e = zeros(m+1, 1); e(i) = 1;
            tableu = [tableu(:, 1:(end-1)) e tableu(:, end)];
            n_holguras = n_holguras + 1;
        elseif strcmp(s, '>=')
            % Agregar holgura negativa.
            e = zeros(m+1, 1); e(i) = -1;
            tableu = [tableu(:, 1:(end-1)) e tableu(:, end)];
            n_holguras = n_holguras + 1;
        end
        fprintf('\t[-] Coeficiente del recurso (la constante): ');
        tableu(i, end) = input('');
    end
    %
    % Buscar factibilidad primal
    %
    if any(tableu(:,end) < 0)
        for i = 1:m
            if tableu(i, end) < 0
                if any(tableu(i, 1:(end-1)) < 0)
                    tableu(i,:) = -1 .* tableu(i, :);
                else
                    fprintf('\n\t [!] Problema infactible porque la fila %d es positiva', i);
                    fprintf('\n\t     y el valor del recurso es negativo, por lo que');
                    fprintf('\n\t     no se puede conseguir factibilidad primal.\n');
                end
            end
        end
    end
    %
    % Buscar identidad inicial
    %
    n_artificiales = 0;
    for i = 1:m
        e = zeros((m + 1),1); e(i) = 1;
        j = 1;
        while j <= (size(tableu, 2) - 1) && not(isequal(tableu(:, j), e))
            j = j + 1;
        end
        if j == size(tableu, 2)
            % Agregamos varible artificial
            tableu = [tableu(:, 1:(end-1)) e tableu(:, end)];
            n_artificiales = n_artificiales + 1;
        end
    end
    ind_vb = zeros(m,1);
    % Buscar indices de base inicial
    for i = 1:m
        e = zeros((m + 1),1); e(i) = 1;
        j = 1;
        while j <= (size(tableu, 2) - 1) && not(isequal(tableu(:, j), e))
            j = j + 1;
        end
        ind_vb(i) = j;
    end
    %
    % Datos de funcion objetivo
    %
    c = zeros(1, (size(tableu, 2) - 1));
    % Costo M de variables artificiales.
    M = 100;
    c((end-(n_artificiales-1)):end) = M;
    if strcmp(tipo, 'max')
        c((end-(n_artificiales-1)):end) = -1 * c((end-(n_artificiales-1)):end);
    end
    fprintf('\n[+] Funcion objetivo:\n\n')
    for j = 1:n
        fprintf('\t[-] Coeficiente de la variable %d: ', j);
        c(j) = input('');
    end
    % Transformar a problema de minimizacion.
    if strcmp(tipo, 'max')
        c = -1 * c;
    end
    %
    % Calcular costos reducidos
    %
    for j = 1:(size(tableu, 2) - 1)
        cr = 0;
        for i = 1:m
            cr = cr + tableu(i,j) * c(ind_vb(i));
        end
        tableu(end, j) = cr - c(j);
    end
    %
    % Calcular valor de la funcion objetivo
    %
    f = 0;
    for i = 1:(size(tableu, 1) -1)
        f = f + tableu(i,end) * c(ind_vb(i));
    end
    tableu(end,end) = f;
end

function [tableu, ind_vb, terminado] = iterar(tableu, ind_vb)
    %
    % Encontrar pivote
    %
    [~, cp] = max(tableu(end, 1:(end - 1)));
    cp = cp(1);
    if all(tableu(1:(end-1), cp) <= 0)
        terminado = 1;
        fprintf('\n[!] El valor de la funcion objetivo es no-acotado.\n');
        fprintf('\n    Existe un rayo optimo de la forma:\n');
        fprintf('\n    (');
        x = zeros(size(tableu, 2) - 1, 1);
        for i = 1:(size(tableu, 1) - 1)
            x(ind_vb(i)) = tableu(i,end);
        end
        for i = 1:size(x)
           fprintf(' %.2g ', x(i));  
        end
        y = zeros(size(x));
        for i = 1:size(ind_vb)
            y(ind_vb(i)) = -1 * tableu(i, cp);
        end
        y(cp) = 1;
        fprintf(') + x_%d * (', cp);
        for i = 1:size(y)
           fprintf(' %.2g ', y(i));  
        end
        fprintf(')\n\n    con x_%d >= 0 y f(x) = %.2g - %.2g * x_%d\n', ...
            cp, tableu(end,end), tableu(end, cp), cp);
    else
        b_a = tableu(1:(end - 1), end) ./ tableu(1:(end - 1), cp);
        if any(b_a) >= 0
            terminado = 0;
            min_b_a = min(b_a(b_a >= 0));
            fp = find(b_a == min_b_a);
            fp = fp(1);
            % El pivote esta en tableu(fp, cp)
            ind_vb(fp) = cp;
            % 
            % Ajustes con pivote
            %
            tableu(fp,:) = tableu(fp,:) ./ tableu(fp, cp);
            for i = 1:size(tableu, 1)
                if not(i == fp)
                    tableu(i,:) = tableu(i,:) - tableu(fp,:) .* tableu(i, cp);
                end
            end
        else
            fprintf('\n[+] Region factible no acotada.\n');
            terminado = 1;
        end 
    end
end

function [tableu, ind_vb, terminado] = iterarBland(tableu, ind_vb)
    %
    % Encontrar pivote
    %
    ops_col = tableu(end, 1:(end - 1));
    pos = ops_col(ops_col > 0);
    % Escogemos el minimo indice:
    pos = pos(1);
    cp = find(ops_col == pos);
    cp = cp(1);
    if all(tableu(1:(end-1), cp) <= 0)
        terminado = 1;
        fprintf('\n[!] El valor de la funciono objetivo es no-acotada.\n');
        fprintf('\n    Existe un rayo optimo de la forma:\n');
        fprintf('\n    (');
        x = zeros(size(tableu, 2) - 1, 1);
        for i = 1:(size(tableu, 1) - 1)
            x(ind_vb(i)) = tableu(i,end);
        end
        for i = 1:size(x)
           fprintf(' %.2g ', x(i));  
        end
        y = zeros(size(x));
        for i = 1:size(ind_vb)
            y(ind_vb(i)) = -1 * tableu(i, cp);
        end
        y(cp) = 1;
        fprintf(') + x_%d * (', cp);
        for i = 1:size(y)
           fprintf(' %.2g ', y(i));  
        end
        fprintf(')\n\n    con x_%d >= 0 y f(x) = %.2g - %.2g * x_%d\n', ...
            cp, tableu(end,end), tableu(end, cp), cp);
    else
        terminado = 0;
        ops_fil = tableu(1:(end - 1), cp) > 0;
        % Escogemos el minimo indice:
        fp = find(ind_vb == min(ind_vb(ops_fil)));
        fp = fp(1);
        % El pivote esta en tableu(fp, cp)
        ind_vb(fp) = cp;
        % 
        % Ajustes con pivote
        %
        tableu(fp,:) = tableu(fp,:) ./ tableu(fp, cp);
        for i = 1:size(tableu, 1)
            if not(i == fp)
                tableu(i,:) = tableu(i,:) - tableu(fp,:) .* tableu(i, cp);
            end
        end
    end
end

function sigue_factible = revisarFactibilidad(...
    tableu, ind_vb, n_artificiales, ind_artificiales, M)
    
    n_va_np = sum(tableu(end, 1:(end - 1)) <= -.75 * M);
    if n_va_np == n_artificiales && not(isempty(intersect(ind_vb, ind_artificiales)))
        sigue_factible = 0;
    else
        sigue_factible = 1;
    end
end

function revisarMultiples(tableu, ind_vb)
    ind_vnb = setdiff(1:(size(tableu, 2) - 1), ind_vb);
    if any(tableu(end, ind_vnb) == 0)
        cp = ind_vnb(tableu(end, ind_vnb) == 0);
        cp = cp(1);
        b_a = tableu(1:(end - 1), end) ./ tableu(1:(end - 1), cp);
        if any(b_a) >= 0
            min_b_a = min(b_a(b_a >= 0));
            fprintf('\n[+] Hay soluciones multiples de la forma:\n');
            fprintf('\n    (');            
            x = zeros(size(tableu, 2) - 1, 1);
            for i = 1:(size(tableu, 1) - 1)
                x(ind_vb(i)) = tableu(i,end);
            end
            for i = 1:size(x)
               fprintf(' %.2g ', x(i));  
            end
            y = zeros(size(x));
            for i = 1:size(ind_vb)
                y(ind_vb(i)) = -1 * tableu(i, cp);
            end
            y(cp) = 1;
            fprintf(') + x_%d * (', cp);
            for i = 1:size(y)
               fprintf(' %.2g ', y(i));  
            end
            fprintf(')\n\n    con 0 <= x_%d <= %.2g\n', cp, min_b_a);
        end
    end
end

function solucion(tableu, tipo, ind_vb, n_tableu, ...
    n, n_holguras, n_artificiales, M, terminado, bland)

    if not(terminado)
        fprintf('\n[+] El algoritmo termino exitosamente: SI');
        fprintf('\n[+] Resultados de la solucion optima.\n');
    else
        fprintf('\n[+] El algoritmo termino exitosamente: NO');
        fprintf('\n[+] Resultados de la ultima iteracion.\n');
    end
    fprintf('\n[+] Valor de M:\t\t\t%d', M);
    fprintf('\n[+] Iteraciones:\t\t%d', (n_tableu - 2));
    fprintf('\n[+] Regla de Bland:\t\t%d', bland);
    fprintf('\n[+] Estructura de x:\n');
    fprintf('\n\t%d variables de decision,', n);
    fprintf('\n\t%d variables de holgura,', n_holguras);
    fprintf('\n\t%d variables artificiales.\n', n_artificiales);
    fprintf('\n[+] Indices de la base:\n\n');
    disp(ind_vb);
    x = zeros(size(tableu, 2) - 1, 1);
    for i = 1:(size(tableu, 1) - 1)
        x(ind_vb(i)) = tableu(i,end);
    end
    fprintf('[+] x = \n\n');
    disp(x);
    if strcmp(tipo, 'max')
        fprintf('[+] Valor de la funcion objetivo evaluada en x: %.2g\n\n', -tableu(end,end));
    else
        fprintf('[+] Funcion objetivo evaluada en x: %.2g\n\n', tableu(end,end));
    end
end

