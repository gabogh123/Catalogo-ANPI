function biseccion()
  
  clc;clear;
  
  f = '(exp(x) -2*x - 10)';
  a = 2;
  b = 4;
  tol = 0.00001;
  iterMax = 100;
  
  [xk k error] = Biseccion(f,a,b,tol,iterMax)

 end 
  
 
  
function [xk k error] = Biseccion(f,a,b,tol,iterMax)
  
  pkg load symbolic
  
  f1 = matlabFunction(sym(f));
 
  
  for k = 1:iterMax
    
    if f1(a) * f1(b) < 0
      
      xk = (a+b)/2;
      
      if f1(a)*f1(xk)<0
        
        b = xk;
        
      else
        a = xk;
        
      endif
      
      error = abs(f1(xk));
    
      if error<tol
      
        break
      endif
    
      
    endif
   
  endfor
  
endfunction