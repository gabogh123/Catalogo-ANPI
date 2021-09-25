
## Metodo iterativo de la posicion falsa para la solucion de ecuaciones
## no lineales.
## Entradas: rango inicial, tolerancia minima y funcion a evaluar.
## Salidas: aproximacion de la solucion y numero de iteraciones.
function [x_k, nIterations] = falsa_posicion (a, b, tolerance, functionStr)
  func = str2func(functionStr);
  nIterations = 0;
  x_k = b - func(b) * (b-a) / (func(b) - func(a));
  if func(a)*func(b) <= 0
    while abs(func(x_k)) >= tolerance
      if func(a)*func(x_k) < 0
        b = x_k;
      else
        a = x_k;
      end
      x_k = b - func(b) * (b-a) / (func(b) - func(a));
      nIterations += 1;
    endwhile
  endif
endfunction
