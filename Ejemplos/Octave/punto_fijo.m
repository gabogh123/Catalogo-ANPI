function punto_fijo()
  clc;clear; %
  
  f='log(2*x+1)'; %String que representa la funcion a resolver
  x0 = 7; %Valor inicial
  tol = 0.0000001; %Tolerancia
  iterMax = 100; %Cantidad m?xima de iteraciones
  
  [xk k error] = puntoFijo(f,x0,tol,iterMax) %Se llama a la siguiente funcion y se definen los valores de retorno
  
end


function [xk k error] = puntoFijo(f,x0,tol,iterMax)
  
  pkg load symbolic %Se carga el paquete symbolic
  
  f1 = matlabFunction(sym(f)); %Se convierte de string a representacion matematica
 
 %Definicion de paramteros iniciales
  error = tol+1; 
  xk = x0;
  k = 1;
  iteracion = 0;
  
  %Creacion de vectores para la grafica
  iter = [];
  er = [];
  
  for k=1:iterMax %Condicion de parada
    
    iteracion = 1.0*k;
    xk_n = f1(xk); %Se asigna el valor a xk_n
    error = abs(xk_n - xk); %Se calcula el error
    xk = xk_n; %Se redefine el valor de xk
    
    iter = [iter;k]; %Se concatena la iteracion al vector de iteraciones
    er = [er;error]; %Se concatena el error al vector de errores
    
    if error<tol %Segunda condicion de parada 
      break;
    endif
    
  endfor
   
end
  
  
  

