function conjugado()
  clc; clear; 
  %Calcular el cero de exp(x)-2*x-10=0
  %Paso 1: Conocer el intervalo donde 
  %        se encuentra el cero
  %        Sugerencia: graficar la funcion
  %        f(x)=exp(x)-2*x-10
  

  
  %Concluimos que la funcion tiene un cero 
  %en el intervalo [2,4]


  f = '(x-2)^4 + (x-2*y)^2';
  x0 = [0 ;3];
  tol = 0.0000001;
  iterMax = 1000;
  [xk k error] = gcln(f,x0,tol,iterMax)
end

function [xk k error] = gcln(f,x0,tol,iterMax)
  pkg load symbolic;
  f1 = matlabFunction(sym (f))
  error = tol + 1;
  k = 1;
  x = sym('x');
  y = sym('y');
  grad_f = gradient(f1, [x, y])
  g_k = grad_f(x0);
  d_k = -1.0*g_k;
  x_k = x0;
  %Creacion de vectores para la grafica
  iter = [];
  er = [];
  while k < iterMax
    
    alpha_k = 1; %delta en ]0, 1[
    sigma = 0.5;
    for i = 1:iterMax
      arr = x_k + alpha_k*d_k;
      izq = f(arr + alpha_k*d_k)-f(x_k); 
      der = sigma*alpha_k*sum(g_k*d_k);
      if izq < der
        break;
      else
        alpha_k = alpha_k/2;      
      endif
    endfor
     
    x_k = x_k + alpha_k*d_k;
    g_k1 = grad_f(x_k);
    error = norm(g_k1);
    iter = [iter;k]; %Se concatena la iteracion al vector de iteraciones
    er = [er;error]; %Se concatena el error al vector de errores
    if error < tol
      break;
    endif 
    
    beta = norm(g_k1).^2/norm(g_k).^2 
    g_k1 = g_k1;
    d_k = -1.0*g_k1 + d_k*beta
    k = k + 1
  endwhile
  plot(iter,er) %Se plotean las iteraciones vs error
  grid on; %Se muestra el grid
end
