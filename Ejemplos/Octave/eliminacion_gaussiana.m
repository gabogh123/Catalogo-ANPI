## Metodo directo de eliminacion_gaussiana para solucion de sistemas de 
## ecuaciones.
## Entradas: Matriz de Coeficientes y matriz de términos independientes.
# Salidas: Solucion.
function [result] = eliminacion_gaussiana(a, b)
    n = length(a);
    function [a, b] =  triangular_superior(a, b, n)
        for j = 1:n
            for i = 1:n
                if i > j && a(i,j) != 0
                    multiplier = -1*a(i,j)/a(j,j);
                    for k = 1:n
                        a(i,k) = a(i,k) + a(j,k) * multiplier;
                    endfor
                    b(i) = b(i) + b(j) * multiplier;
                endif
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
    
    [a, b] = triangular_superior(a, b, n);
    result = sustitucion_hacia_atras(a, b, n);
endfunction
