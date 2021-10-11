pkg load parallel;

function [xk,error] = parte1_p3 (A, b, x, m, tol, iterMax)
  
  q = [1:0.1:25];
  p = [1:0.1:25];
  m = 242;

  A = tridiagonal(p,q,m)
  b = ones(242,1);
  x = zeros(241,1);
  tol = 0.0001;
  iterMax = 1000;
  
  
  error = 0;
  xk = x;
  
  for iter = 1: iterMax
    
    vector_i = 1:m;
    
    sum = @(i)  (b(i)-(A(i,1:i-1)*xk(1:i-1)+A(i,i+1:n)*xk(i+1:n)))/ (A(i,i));
    
    
    xk = pararrayfun(nproc, sum, vector_i);
    
    if norm((A*xk)-b)<tol
      
      break
    endif
    
  endfor
 
endfunction   


function res = tridiagonal(p,q,m)
  

   
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