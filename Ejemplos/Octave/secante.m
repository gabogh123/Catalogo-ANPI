function secante()
  
  clc;clear;
  
  f = 'exp(-x**2)-x';
  x0 = 0;
  x1 = 1;
  tol = 0.00001;
  iterMax = 1000;
  
  [xk k error] = Secante(f,x0,x1,tol,iterMax)
  
  
endfunction

function [xk k error] = Secante(f,x0,x1,tol,iterMax)
  
  pkg load symbolic
  
  f1 = matlabFunction(sym(f));

  
  for k=1:iterMax
    
    n = f1(x1) * (x1-x0);
    d = f(x1)-f(x0);
    
    xk = x1 -(n/d);
    
    x0 = x1;
    x1 = xk;
    
    error = abs(f1(xk));
    
    if error<tol
          break
        endif
        
  endfor
  
endfunction
