
## Metodo iterativo del Gradiente Conjugado para optimizacion en varias
## variables.
## Entradas: valores iniciales, variables de la funcion, tolerancia minima
## del resultado, funcion a evaluar en string y regla a utilizar.
## Salidas:Valor aproximado y numero de iteraciones
function [x_k, iteration] = gradiente_conjugado(x0, str_variables, tol, funcion, regla_bk = 'FR')
    n = length(str_variables);
    f = str2func(funcion);
    derivatives = sym(zeros(1, n));
    for i = 1:n
        variables(i) = sym(str_variables(i));
    endfor
    
    for i = 1:n
        derivatives(i) = diff(f, variables(i));
        i += 1;
    endfor
    x_k = x0;
    error = tol;
    iteration = 0;
    
    function [grad] = gradiente(derivatives, variables, values, n)
        grad = zeros(1, n);
        
        for i = 1:n
            grad(i) = matlabFunction(derivatives(i))(num2cell(values){:});
            i += 1;
        endfor
    endfunction
    
    g_k_1 = g_k = gradiente(derivatives, variables, x_k, n);
    d_k = -1 * g_k;
    
    function [alpha] = get_alpha(f, variables, x_k, d_k, g_k)
        alpha = 1;
        delta = rand(1);
        f_arg = x_k + alpha * d_k;
        while double(f(num2cell(f_arg){:})) - double(f(num2cell(x_k){:})) > delta * alpha * g_k * transpose(d_k)
            alpha /= 2;
            f_arg = x_k + alpha * d_k;
        endwhile    
    endfunction
   
    function [beta] = get_beta(regla, g_k_1, g_k, d_k)
        if lower(regla) == 'cd'
            beta = norm(g_k_1)**2 / ((-1 * d_k)*transpose(g_k));
        elseif lower(regla) == 'dy'
            beta = norm(g_k_1)**2 / d_k*transpose(g_k_1 - g_k);
        else
            beta = norm(g_k_1)**2 / norm(g_k)**2;
        end    
    endfunction
    
    while error >= tol
        alpha = get_alpha(f, variables, x_k, d_k, g_k_1);
        x_k = x_k + alpha * d_k;
        g_k = g_k_1;
        g_k_1 = gradiente(derivatives, variables, x_k, n);
        d_k = -1 * g_k_1 + get_beta(regla_bk, g_k_1, g_k, d_k) * d_k;
        error = norm(g_k_1);
        iteration += 1;
    endwhile
endfunction


