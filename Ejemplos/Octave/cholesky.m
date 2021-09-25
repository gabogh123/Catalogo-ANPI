

## Metodo directo de factorizacion de Cholesky para solucion de sistemas de 
## ecuaciones.
## Entradas: Matriz de Coeficientes y matriz de términos independientes.
# Salidas: Solucion .
function [result] = cholesky(a, b)
    n = length(a);
    
    function [L] = get_L(a, n)
        L = zeros(n);
        for i = 1:n
            for j = 1:n
                if i == j
                    sum1 = 0;
                    for k = 1:j
                        sum1 += L(j, k)**2;
                    endfor
                    L(i, j) = sqrt(a(i, j) - sum1);
                else 
                    if i > j
                        sum2 = 0;
                        for k = 1:j
                            sum2 += L(i, k)*L(j, k);
                        endfor    
                        L(i, j) = (a(i, j) - sum2) / L(j, j);
                    end
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
    L = get_L(a, n);
    y = sustitucion_hacia_adelante(L, b, n);
    result = sustitucion_hacia_atras(transpose(L), y, n);
endfunction    



