    
% Este programa nos regresa L (matriz triangular inferior) U (matriz
% triangular superior) y x (solución al sistema lineal Ax=b)
function  x = FactLU(A,b)
[m,n]=size(A);
if (m ~= n )
disp ( 'LR2 error: La matriz debe ser cuadrada' );
return
end
   %creamos nuestras matrices L y U de ceros para iniciar la factorización
  L=zeros(m,m);
  U=zeros(m,m);
  for i=1:m
  % Encontramos la matriz L
  for k=1:i-1
  L(i,k)=A(i,k);
  for j=1:k-1
  L(i,k)= L(i,k)-L(i,j)*U(j,k);
  end
  L(i,k) = L(i,k)/U(k,k);
  end
  % Encontramos la matriz U
  for k=i:m
  U(i,k) = A(i,k);
  for j=1:i-1
  U(i,k)= U(i,k)-L(i,j)*U(j,k);
  end
  end
  end
  for i=1:m
  L(i,i)=1;
  end
  % Borrar ";" si quieres ver los resultados que se muestran antes de
  % resolver el problema.
  U;
  L;
  % Usamos un vector "y" para resolver 'Ly=b'
  y=zeros(m,1);
  y(1)=b(1)/L(1,1);
  
  for i=2:m
    y(i)=b(i);
    for k=1:i-1
        y(i)=y(i)-L(i,k)*y(k);  
    end
    y(i)=y(i)/L(i,i);
  end  

% Ahora resolvemos con el vector "x" el problema  Ux = y para encontrar
% nuestra solución al problema Ax=b
x=zeros(m,1);
x(m)=y(m)/U(m,m);
for i=m-1:-1:1
    x(i)=y(i);
    for k=i+1:m
        x(i)=x(i)-U(i,k)*x(k);
    end
    x(i)=x(i)/U(i,i);
end



end
