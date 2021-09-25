function newton_raphson()
  clc; clear;
  
  
  f = 'cos(2*x)^2 -x^2' ;
  x0=3/4;
  tol=10^-6;
  iterMax=1000;
  [xk k error]= newtonRaphson(f,x0,tol,iterMax)
  
end

function [xk k error] = newtonRaphson(f,x0,tol,iterMax)
  
  pkg load symbolic
  
  f1 = matlabFunction(sym(f));
  
  %aproxi = [];
 % er = [];
  
  error=tol+1;
  
  for k=1:iterMax
       
    xk=x0;
    k=k+1;
    n=f1(xk);
    d=diff(f1(xk));
    xk=(xk-(n./d));
    
    error=abs(f1(xk));
    
    %aproxi = [aproxi;xk];
    %er = [er;error];
    
    if error<tol
      break
    endif
    
  endfor
 
  %plot(aproxi,er);
 
 end