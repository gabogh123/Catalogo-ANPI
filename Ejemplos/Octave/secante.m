
## Metodo iterativo de la secante para la solucion de ecuaciones no lineales.
## Entradas: rango inicial, tolerancia minima y funcion a evaluar.
## Salidas: aproximacion de la solucion y numero de iteraciones.
function [x_k, nIterations] = secante(x0, x1, tolerance, functionStr)
  func = str2func(functionStr);
  x_k = x1;
  x_k_1 = x0;
  nIterations = 0;
  while abs(func(x_k)) >= tolerance
    temp = x_k;
    numerator = func(x_k) * (x_k - x_k_1);
    denominator = func(x_k) - func(x_k_1);
    x_k = x_k - numerator / denominator;
    x_k_1 = temp;
    nIterations += 1;
  endwhile
endfunction
