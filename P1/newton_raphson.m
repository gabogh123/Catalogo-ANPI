function newton_raphson()
  pkg load symbolic;
  syms x;
 
  f = (cos(2*x)^2 - x^2); #Funcion, esto pueden cambiarlo a conveniencia
  iterMax = 100; %Cantidad máxima de iteraciones
  disp(x - ((f))/diff(f)); % formula
  _newton_raphson (0.75,0,f, iterMax);
  
endfunction

function [xk, error] = _newton_raphson(xi,xk,f, iterMax)
  
  % Esta funcion aproxima numericamente el resultado de la funcion mediante el metodo de newton raphson
% Parametros de entrada: f=funcion a resolver
                        %xi = valor inicial
                        %xk = aproximacion
                        %iterMax = Cantidad de iteraciones maximas    
% Parametros de salida: xk = valor aproximado
                        %error = error del metodo
  
  %Definicion de paramteros iniciales
  syms x;
  f_temp = x - ((f))/diff(f); 
  x = xi;
  xk = eval(f_temp);
  error = abs(eval(f));
  k = 1;
  iteracion = 0;

 
  
  %Creacion de vectores para la grafica
  iter = [];
  er = [];
  
  
  for k=1:iterMax %Condicion de parada
    
    %fprintf("\t\txk = %f\n", xk);
    %fprintf("\t\t|f(xk)| = %f\n", error);
    
    iteracion = 1.0*k;
    syms x;
    f_temp = x - ((f))/diff(f); %se crea una variable temporal equivalente a la formula de nr
    x = xi; %se sustituye la x por el valor de xi
    xk = eval(f_temp); %se define el valor de xk
    error = abs(eval(f)); %se calcula el error
    
    iter = [iter;k]; %Se concatena la iteracion al vector de iteraciones
    er = [er;error]; %Se concatena el error al vector de errores
    
    xi = xk;
    xk = x;
    
    if (error == 0) %Segunda condicion de parada 
      
      %fprintf('Resultado = %f\n',xk);
      %fprintf("\t\t|f(xk)| = %f\n", error);
      
      break;
    endif
  endfor
  plot(iter,er) %Se plotean las iteraciones vs error
  grid on; %Se muestra el grid
endfunction