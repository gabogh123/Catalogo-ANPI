function muller()
  
  f = 'sin(x) - x/2';
  x0 = 2;
  x1 = 2.2;
  x2 = 1.8;
  tol = 10**-5;
  iterMax = 1000;
  
  [xk k error] = Muller(f,x0,x1,x2,tol,iterMax)
  
end

function [xk k error] = Muller(f,x0,x1,x2,tol,iterMax)
  
  pkg load symbolic
  
  f1 = matlabFunction(sym(f));

  for k=1:iterMax
    
    n0 = x0-x1;
    n1 = x0-x2;
    n2 = x1-x2;
    
    deno = n0*n1*n2;
    
    h12 = f1(x1) - f1(x2);
    h02 = f1(x0) - f1(x2);
    
    a = (n2*h02-n1*h12)/deno;
    b =  (n1**2*h12 - n2**2*h02)/deno;
    c = f1(x2);
    
    d = b**2 - 4*a*c;
    
    if abs(b-d) < abs(b+d)
      aprox = (2*c)/(b+d);
      
    else
      aprox = (2*c)/(b-d);
      
    endif
    
    xk = x2-aprox;
    error = abs( f(xk));
    
    if error < tol
      break
    endif
    
    x0 = x1;
    x1 = x2;
    x2 = xk;
    
  endfor
  

endfunction
  
  