## Metodo iterativo de gauss hacia atras para solucion de sistemas de 
## ecuaciones.
## Entradas: Matriz de Coeficientes, matriz de terminos independientes,
## Valor inicial y tolerancia.
## Salidas: aproximacion y cantidad de iteraciones.

function [x, i] = gauss_atras(A, b, x0, tol)
  n = length(A);
  L = D = U = zeros(n);
  
  for i = 1:n
    for j = 1:n
      if i > j
        L(i, j) = A(i, j);
      elseif i < j
        U(i, j) = A(i, j);
      else
        D(i, j) = A(i, j);
      endif
    endfor
  endfor
  
  x = x0;
  error = norm(A * x - b);
  i = 0;
  
  m_inv = inv(L + D);
  m_times_n = m_inv * (-U);
  m_times_b = m_inv * b;
  
  while error >= tol
    x = m_times_n * x + m_times_b;
    error = norm(A * x - b);
    i += 1;
  endwhile
  
endfunction
