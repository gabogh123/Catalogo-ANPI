function parte1_p2()
  clc; clear;

 %Se inicializan los valores para el método 
  q = [1:0.1:25];
  p = [1:0.1:25];
  m = 242;

  A = tridiagonal(p,q,m)
  b = ones(242,1);
  x = zeros(241,1);
  tol = 0.00001;
  iterMax = 1000;
  
  [xk k] = jacobimetodo(A,b,x,tol,iterMax)

endfunction



function [xk k] = jacobimetodo(A,b,x,tol,iterMax)
  
%Esta función calcula el método de Jacobi mediante la fórmula de sumatoria de los valores de la matriz A
  
  n = length(x);
  
  for k=1:iterMax %Primer criterio de parada
    
    
    xk=x; %Se inicializa el xk en el x0
    for i=1:n
      
      sum = A(i,1:i-1)*xk(1:i-1)+A(i,i+1:n)*xk(i+1:n); %Se realiza la sumatoria
      x(i) = (b(i)-sum) / A(i,i);
      
    endfor
    
    if norm(xk-x,inf)<tol %Segundo criterio de parada 
      
      break 
      
    endif
    
  endfor
  
endfunction

function res = tridiagonal(p,q,m)
  
%Función que calcula la amtriz de coeficientes a partir de los vectores p y q

  m = length(p);
  if (length(p) == m)
      
      A = zeros(m,m);
      
      A(1,1) = 2*q(1); %Se definen las esquinas de la matriz
      A(m,m) = 2*p(m);
      
      for k=2:m     
        A(k-1,k) = q(k-1); %Se rellena la matriz
        A(k,k-1) = p(k);  
      endfor
      
      for k=2:m-1
        
        A(k,k) = 2*(p(k)+q(k)); %Se rellena la matriz
       
      endfor
    
  else
    display('Error: Los vectores p y q no son del mismo tamaño');    
    
  endif 
  
  res = A
    
endfunction

