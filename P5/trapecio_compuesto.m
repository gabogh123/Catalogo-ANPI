function trapecio_compuesto()
  clc; clear;
  f = 'ln(x)';
  inter = [2; 5];
  n = 500;
  [aprox error] = trapecio_aux(f, inter, n)
end

function [aprox error] = trapecio_aux(f, inter, n)
  %La regla del trapecio compuesto es un m�todo utilizado para aproximar
  %la integral de f(x) para un intervalo [a, b]
  %Parametros de entrada: 
            %f : funcion string para aproximar integral
            % inter : intervalo de integraci�n [a, b]
            % n: cantidad de puntos a evaluar en el intervalo
  %Salida:
            %aprox: aproximaci�n de la integral
            %error: cota de error m�ximo para la aproximaci�n calculada 
  pkg load symbolic;
  syms x;
  fs = sym (f); %Funcion simb�lica
  f1 = matlabFunction(sym (f)); %funci�n num�rica
  a = inter(1);
  b = inter(2);
  xv = linspace(a, b, n);
  h = (b-a)/(n-1);
  aprox = 0;
  
  %Ciclo donde se calcula la aproximaci�n siguiendo la f�rmula
  % Atc = h   _n-1_
  %       - * \     (f(x[i]) + f(x[i+1]))
  %       2   /____ 
  %             i=0
  for k=1:n-1
    current = xv(k);
    next = xv(k+1);
    f_c = f1(current);
    f_n = f1(next);
    aprox += h*(f_c+f_n)/2;
  endfor
  %C�lculo de la cota de error para la regla del trapecio compuesto
  dxx = abs(diff(diff(fs))); %segunda derivada de |f(x)| simb�lica
  dxx_n = matlabFunction(dxx);
  dxx_aux = -dxx;            %segunda derivada de -|f(x)| simb�lica                 
  xmax = fminbnd(matlabFunction(dxx_aux),a,b);  % max |f"(x)| = min -|f"(x)| x en [a, b]
  alpha_max = dxx_n(xmax);
  error = (b-a)*(h^2/12)*alpha_max;
end 

  