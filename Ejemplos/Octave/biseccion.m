function [x_k, nIterations] = biseccion (a, b, tolerance, functionStr)
  func = str2func(functionStr);
  x_k = (a + b) / 2;
  nIterations = 0;
  if (func(a) * func(b) < 0)
    while abs(func(x_k)) >= tolerance
      if func(a) * func(x_k) < 0
        b = x_k;
      else
        a = x_k;
      end
      nIterations += 1;
      x_k = (a + b) / 2;
    endwhile
  endif
endfunction
