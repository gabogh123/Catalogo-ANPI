pkg load symbolic %carga paquete simbolico
clc; 

function [aproximacion] = potencia(A, x)
  % funcion que aproxima el vector propio 
  % entradas: una matriz y un vector inicial
  % salidas: la aproximacion de los valores y vectores propios
    tol = 10^-10
    norma_inf_ant = 0;
    
    while 1     
        % matriz multiplicada por el vector inicial
        y = A * x;
        % calcilo de la norma inferior
        normaInf = norm(y, inf);
        y = y/normaInf;
        % calculo de ck
        ck = abs(norma_inf_ant - normaInf);
        % calculo del error
        v = norm(x - y);
        error = max(ck, v);

        % condicion de parada
        if error < tol
            break;
        endif
        
        % se asigna nuevo valor de matriz inferior
        norma_inf_ant = normaInf;
        x = y;     
    endwhile

    % resultado
    aproximacion = normaInf;
    
    % grafico de error vs iteraciones
    plot(0:length(error)-1, error)
    xlabel('Iteracion')
    ylabel('Error')
endfunction

% prueba del metodo para la funcion
A = [3 -1 0; -1 2 -1; 0 -1 3]
x = [1 1 1]'
[aproximacion] = potencia(A, x)