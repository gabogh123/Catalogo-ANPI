function falsa_posicion()
  
  f = 'cos(x)-x';
  a = 1/2;
  b = pi/4;
  tol = 10**-5;
  iterMax = 1000;
  
  [xk k error] = falsaPosicion(f,a,b,tol,iterMax)
  
end

function [xk k error] = falsaPosicion(f,a,b,tol,iterMax)
  
  pkg load symbolic
  f1 = matlabFunction(sym(f)); 
  

  
  for k=1:iterMax
  
    
    if f1(a)*f1(b) < 0
      
      xk = b - ((f1(b)*(b-a))/ (f1(b) - f1(a)));
      
      if f1(a)*f1(xk) < 0
        
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


