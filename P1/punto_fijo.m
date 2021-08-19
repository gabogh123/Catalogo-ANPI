function punto_fijo()
  clc;clear;
  
  f='log(2*x+1)';
  x0 = 7;
  tol = 0.0000001;
  iterMax = 100;
  
  [xk k error] = puntoFijo(f,x0,tol,iterMax)
  
end


function [xk k error] = puntoFijo(f,x0,tol,iterMax)
  
  pkg load symbolic
  
  f1 = matlabFunction(sym(f));
 
  error = tol+1;
  xk = x0;
  k = 1;
  iteracion = 0;
  
  iter = [];
  er = [];
  
  for k=1:iterMax
    
    iteracion = 1.0*k;
    xk_n = f1(xk);
    error = abs(xk_n - xk);
    xk = xk_n;
    
    iter = [iter;k];
    er = [er;error];
    
    if error<tol
      break;
    endif
    
  endfor
  
  plot(iter,er)
  grid on;
  
end
  
  
  

