function conjugado()
  clc; clear; 
  f = '(x-2)^4 + (x-2*y)^2'; %String que representa la funcion a evaluar
  x0 = [0 ;3]; %Vector con valores iniciales
  tol = 0.0000001; %Valor de tolerancia aceptada
  iterMax = 1000; %Valor máximo de iteraciones
  [x_k k error] = gradiente_conjugado(f,x0,tol,iterMax)
end

function [x_k k error] = gradiente_conjugado(f,x0,tol,iterMax)
  %El método de gradiente conjugado no lineal permite resolver 
  %problemas de optimización en variable variable, por ello, 
  %en la función gradiente_conjugado
  %se desea conocer un valor mínimo de la función.
  %Parámetros de entrada:
  %     f = función en estudio       
  %     x0 = vector inicial
  %     tol = Valor de tolerancia aceptada
  %     iterMax = Cantidad máxima de iteraciones
  
  pkg load symbolic;
  f1 = matlabFunction(sym (f))
  error = tol + 1;
  k = 1;
  x = sym('x');
  y = sym('y');
  grad_f = matlabFunction(gradient(f1,[x y])) %Se averigua la función gradiente de f(x,  y)
  g_k = grad_f(x0(1), x0(2)); %Se evalúa el vector incial en el gradiente
  d_k = -1.0*g_k;
  x_k = x0;
  %Creacion de vectores para la grafica
  iter = [];
  er = [];
  while k < iterMax
    
    alpha_k = 1; 
    sigma = 0.5; %sigma en ]0, 1[
    %Ciclo para calcular el tamaño de paso alpha_k
    for i = 1:iterMax
      arr = x_k + alpha_k*d_k;
      izq = f1(arr(1), arr(2)) - f1(x_k(1),x_k(2)); 
      der = sigma*alpha_k*sum(g_k'*d_k);
      if izq < der
        break;
      else
        alpha_k = alpha_k/2;      
      endif
    endfor
     
    x_k = x_k + alpha_k*d_k;
    g_k1 = grad_f(x_k(1), x_k(2));
    error = norm(g_k1);
    iter = [iter;k]; %Se concatena la iteracion al vector de iteraciones
    er = [er;error]; %Se concatena el error al vector de errores
    if error < tol
      break;
    endif 
    
    beta = norm(g_k1)^2/norm(g_k)^2; %Selección del parámetro Beta según la regla de Fletcher y Reeves
    g_k = g_k1;
    d_k = -1.0*g_k1 + d_k*beta;
    k = k + 1;
  endwhile
  plot(iter,er) %Se plotean las iteraciones vs error
  grid on; %Se muestra el grid
end
