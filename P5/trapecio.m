function trapecio ()
  clc; clear;
  f = 'ln(x)';
  inter = [2; 5];
  [aprox error] = trapecio_aux(f, inter)
end

function [aprox error]=trapecio_aux(f, inter)
  %La regla del trapecio es un m�todo utilizado para aproximar
  %la integral de f(x) para un intervalo [a, b]
  %Parametros de entrada: 
            %f : funcion string para aproximar integral
            % inter : intervalo de integraci�n [a, b]
  %Salida:
            %aprox: aproximaci�n de la integral
            %error: cota de error m�ximo para la aproximaci�n calculada 
  pkg load symbolic;
  syms x;
  fs = sym (f); %Funcion simb�lica
  f1 = matlabFunction(sym (f)); %funci�n num�rica
  a = inter(1);
  b = inter(2);
  aprox = (b-a)*(f1(a) + f1(b))/2;
  
  %Cota de error
  dxx = abs(diff(diff(fs))); %segunda derivada de |f(x)| simb�lica
  dxx_n = matlabFunction(dxx);
  dxx_aux = -dxx;            %segunda derivada de -|f(x)| simb�lica
  h = b-a;                  
  xmax = fminbnd(matlabFunction(dxx_aux),a,b);  % max |f"(x)| = min -|f"(x)| x en [a, b]
  alpha_max = dxx_n(xmax);
  error = (h^3/12)*alpha_max;
end

  