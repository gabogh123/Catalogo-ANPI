pkg load parallel;

function [xk,error] = parte1_p3(A,b,x,m,tol,iterMax)
  
  p = [1:0.1:25];
  q = [1:0.1:25];
  m = 242;
  
  A = tridiagonal(p,q,m);
  
  b = ones(242,1);
  x = zeros(241,1);
  tol = 0.0001;
  iterMax = 1000;
  
  error = 0;
  
  xk = x;
  
  for iter = 1:iterMax
    
    vector_i = 1:m;
    
    xk = pararrayfun(nproc, @(i) sumatoria(m,A,b,xk,i), vector_i);
    
    error = norm((A*xk')-b);
    
    if error < tol
      break
      
    endif
    
  endfor
endfunction 

function xk = sumatoria(A,b,m,xk,i)
  
  suma = 0;
  for j = 1:m
    
    if j != i
      
      suma = suma+(A(i,j)*xk(j));
      
    endif
    
  endfor
  
  xk = (1/A(i,i)) * (b(i) - suma);
  
  return;
endfunction 


function res = tridiagonal(p,q,m)
  
  q = [1:0.1:25];
  p = [1:0.1:25];
  m = 242;
   
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
  
  res = A;
    
endfunction