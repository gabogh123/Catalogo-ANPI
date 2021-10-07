function jacobi()
  clc; clear;
  
  q = [1:1:30];
  p = [31:1:60];
  m = 30;

  A = tridiagonal(p,q,m)
  b = [5 6 4]';
  x = [1 1 1]';
  tol = 0.0001;
  iterMax = 1000;
  
  [xk k] = jacobimetodo(A,b,x,tol,iterMax)

endfunction



function [xk k] = jacobimetodo(A,b,x,tol,iterMax)
  
  n = length(x);
  
  for k=1:iterMax
    
    xk=x;
    for i=1:n
      
      s = A(i,1:i-1)*xk(1:i-1)+A(i,i+1:n)*xk(i+1:n);
      x(i) = (b(i)-s) / A(i,i);
      
    endfor
    
    if norm(xk-x,inf)<tol
      
      break 
      
    endif
    
  endfor
  
endfunction

function res = tridiagonal(p,q,m)
  
  q = [1:1:30];
  p = [31:1:60];
  m = 30;  
   
  m = length(p);
  if (length(p) == m)
      
      A = zeros(m,m);
      
      A(1,1) = 2*q(1);
      A(m,m) = 2*p(m);
      
      for k=2:m     
        A(k-1,k) = q(k-1);
        A(k,k-1) = p(k);  
      endfor
      
      for k=2:m-1
        
        A(k,k) = 2*(p(k)+q(k));
       
      endfor
    
  else
    display('Error: Los vectores p y q no son del mismo tamaño');    
    
  endif 
  
  res = A
    
endfunction

