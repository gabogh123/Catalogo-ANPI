function adam_bashford_4 ()
  
  warning('off','all')
    clc; clear;
    
  % Ejemplo presentacion 12, diapositiva 64

  fun = '1 + (x-y)^2';
  itervalo = [2 4];
  num_puntos = 11;
  ini = [1 1.191 1.5964 1.8883];
  [xv, yv, polinomio] = metodo_adam_bashford_4(fun, itervalo, num_puntos, ini)
  
  % Grafica de la solucion
  hold on
  x_g=itervalo(1):0.0001:itervalo(2);
  polinomio1 = matlabFunction(sym(polinomio));
  y_p=polinomio1(x_g);
  plot(x_g,y_p,'b')
  title('Metodo Adam Bashford')
  xlabel('x')
  ylabel('y(x)')
  grid on;
  

  
end


function pol=lagrange(xv,yv)
  
  % Metodo de Lagrange para el caclculo del polinomioinomio
  % de itervalopolinomioacion
  
  % Parametros de entrada:
  %   xv, yv: vetores de los pares ordenador a partir de
  %     los cuales se calculara el polinomioinomio.
  
  % p =  polinomioinomio de itervalopolinomioacion de forma simbolica.
  
  syms x
  n=length(xv)-1;
  pol=0;
  for k=0:n
    pol=pol+yv(k+1)*funcion_Lk(xv,k);
  end
  pol=expand(pol);
end

function funcion_Lk=funcion_Lk(xv,k)
  syms x
  n=length(xv)-1;
  funcion_Lk=1;
  for j=0:n
    if j~=k
      funcion_Lk=funcion_Lk*(x-xv(j+1))/(xv(k+1)-xv(j+1));
    end    
  end
  funcion_Lk=expand(funcion_Lk);
end


function [xv, yv, polinomio] = metodo_adam_bashford_4(fun, itervalo, num_puntos, ini)
  
%Funcion que calcula la solucion de una ecuacion diferencial 
%Entradas: fun: funcion a evaluar
%          itervalo: es el intervalo de la funcion a evaluar
%          num_puntos: numero de puntos a analizar
%          ini: corresponde al conjunto de valores iniciales
  
% Parametros de salida:
% xv, yv:  vectores de valores aproximados 
% polinomio: polinomio de interpolacion
  
  
  
%Inicializacion de valores

  pkg load symbolic;
  syms x y;
  f1 = matlabFunction(sym(fun));
  
%Se calcula el paso h
  a = itervalo(1);
  b = itervalo(2);
  
  h=(b-a)/(num_puntos-1);
  xv=a:h:b;
  
  % calculo de yv
  y0 = ini(1);
  y1 = ini(2);
  y2 = ini(3);
  y3 = ini(4);
  yv = [y0 y1 y2 y3];

  for n=4:num_puntos-1  
    fk = f1(xv(n), yv(n));
    fk_1 = f1(xv(n-1), yv(n-1));
    fk_2 = f1(xv(n-2), yv(n-2));
    fk_3 = f1(xv(n-3), yv(n-3));
    yv(n+1)= yv(n)+ (h/24) * (55*fk - 59*fk_1 + 37*fk_2 - 9*fk_3);
  end
  
%Se calcula el polinomio mediante el metodo de lagrange
  polinomio = lagrange(xv, yv);
  
end

