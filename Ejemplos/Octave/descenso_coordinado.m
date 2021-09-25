
## Metodo iterativo del Descenso Coordinado para optimizacion en varias
## variables.
## Entradas: valores iniciales, variables de la funcion, tolerancia minima
## del resultado y funcion a evaluar en string.
## Salidas:Valor aproximado y numero de iteraciones
function [x0, iteration] = descenso_coordinado(x0, str_variables, tol, funcion)
    warning('off');
    n = length(str_variables);
    f = str2func(funcion);
    error = tol;
    iteration = 0;
    for i = 1:n
        variables(i) = sym(str_variables(i));
    endfor
    
    function x0 = replace(x0, indice, variables, f)
        temp = num2cell(x0);
        temp(indice) = variables(indice);
        df = diff(f(temp{:}));
        soluciones = solve(df);
        m = length(soluciones);
        for i = 1:m
            temp(indice) = double(soluciones(i));
            if f(temp{:}) < f(num2cell(x0){:})
            x0 = cell2mat(temp);
            endif
        endfor
    endfunction
    
    
    while error >= tol
        prev = x0;
        for j = 1:n
            x0 = replace(x0, j, variables, f);
        endfor
        iteration += 1;
        error = norm(x0 - prev);
    endwhile
endfunction