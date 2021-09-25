
## Metodo directo de factorizacion LU para solucion de sistemas de 
## ecuaciones.
## Entradas: Matriz de Coeficientes y matriz de términos independientes.
# Salidas: Solucion .
function [result] = factorizacion_LU(a, b)
    n = length(a);
    
    function [L, a] =  get_LU(a, n)
        L = zeros(n);
        for j = 1:n
            for i = 1:n
                if i > j && a(i, j) != 0
                    multiplier = -1 * a(i, j) / a(j, j);
                    L(i, j) = -1 * multiplier;
                    for k = 1:n
                        a(i, k) = a(i, k) + a(j, k) * multiplier;
                    endfor
                elseif i == j
                    L(i, j) = 1;
                end
            endfor
        endfor
    endfunction 
    
    function [result] = sustitucion_hacia_atras(a, b, n)
        result = zeros(1, n);
        i = n;
        while i > 0
            resultTemp = b(i);
            j = n;
            while j > 0
                if i == j
                    result(j) = resultTemp/a(i,j);
                    break
                else
                    resultTemp -= a(i,j)*result(j);
                end
                j -= 1;
            endwhile    
            i -= 1;
        endwhile
    endfunction
    
    function [result] = sustitucion_hacia_adelante(a, b, n)
        result = zeros(1, n);
        for i = 1:n
            resultTemp = b(i);
            for j = 1:n
                if i == j
                    result(j) = resultTemp/a(i,j);
                    break
                else
                    resultTemp -= a(i,j)*result(j);
                end
            endfor
        endfor    
    endfunction
    
    [L, U] = get_LU(a, n);
    y = sustitucion_hacia_adelante(L, b, n);
    result = sustitucion_hacia_atras(U, y, n);
endfunction
   