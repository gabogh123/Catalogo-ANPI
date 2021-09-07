function p3_solucion_aplicacion()
  pkg load symbolic %Se carga el paquete symbolic
 
  f = '(5512026*x^2)+(28471210*y^2)-(5000*x)-(8660*y)'; %String que representa la funcion a evaluar
  xk = [0.000786 ; 0.0000878] % se escogen estos valores de x y y que representan el minimo de la funcion
  tol = 0.00001; %Valor de tolerancia
  iterMax = 1000;
  tolerancia = 0.00001;
  
 bfgs(f,xk,tol,iterMax, tolerancia);
  
  
end

function bfgs(f,xk,tol,iterMax, tolerancia) 
  
  % se define la funcion de prueba
 
  f1 = matlabFunction(sym (f))
  x = sym('x');
  y = sym('y');

  % se define el gradiente de la funcion de prueba
  grad_f = matlabFunction(gradient(f1,[x,y])) 

  %calculo del gradiente
  gk = grad_f(xk(1), xk(2));
  
  n = length(xk);
  Bk = eye(n);   % matriz identidad nxn donde n=cantidad de variables de entrada
  
  iteraciones = [];
  errores = []; 
  
  for i=1:iterMax
    
   % Obtener pk de la ecuacion: Bk*pk = -g(xk), donde g(xk) es el gradiente de la funcion, entonces pk = Bk^-1 * -gk
   Bk_inv = inv(Bk);
   pk = Bk_inv * -gk;
   
   % Determinar el stepsize con λk > 0 a partir de: f(xk + λk*pk) <= f(xk) + delta*λk*g(xk)^T * pk
   lambda = stepsize(xk, pk, gk, f1);
   
   % A partir de estos valores calcular xk+1 = xk_pasado + λk*pk
   
   % primero guardamos los xk y gk anteriores para usarlos al calcular xk+1, sk y yk
   xk_pasado = xk;
   gk_pasado = gk;
   
   % se calcula el nuevo xk y gk
   xk = xk + (lambda*pk);
   gk = grad_f(xk(1), xk(2));
   
   
   % Luego calcular Bk+1 con su respectiva funcion (2.10 en el documento) donde sk = xk - xk_pasado y yk = gk - gk_pasado
   
   % calculando sk y yk
   sk = xk - xk_pasado;
   yk = gk - gk_pasado;
   
   sk_t = transpose(sk);
   yk_t = transpose(yk);
   
    % Calculando el valor de Bk donde si se cumple cierta condicion cambia y si no la cumple se mantiene igual.
   if (((yk_t * sk) / norm(sk)^2) >= norm(gk))
      der = ((yk*yk_t) / (yk_t*sk));
      izq = ((Bk*sk*sk_t*Bk) / (sk_t*Bk*sk));   
      Bk = Bk - izq + der;
   else
      Bk = Bk;
   endif
   
   % Agregando los valores para graficar
   iteraciones = [iteraciones ; i];
   errores = [errores ; norm(gk)];
   
   % condicion de parada 
   if (norm(gk) < tolerancia)
     break;
   endif
   
  endfor
  disp("Valor final de xk:");
  xk
  plot(iteraciones, errores);  
  title("Metodo BFGS (Iteraciones vs Error)");
  xlabel("Iteraciones");
  ylabel("Error");
  set(gca, "fontsize", 16);
  grid on; %Se muestra el grid   
  
end

% esta función obtiene el gradiente de la funcion de prueba y lo retorna como un vector columna
function res = get_gradient(xk, grad_f)
  
  g_vec = [];
  for i = 1:length(xk)    
    g = grad_f(xk(i));
    g_vec = [g_vec ; g];
  endfor
  res = g_vec;
  
end

% esta funcion se utiliza para escoger el paso de busqueda λk
function res = stepsize(xk, pk, gk, fcn)
  delta = 0.5;
  lambda = 0.8;
  xk_lambda = xk + lambda * pk;
  f1 = fcn(xk_lambda(1), xk_lambda(2));
  f2 = fcn(xk(1), xk(2));
  
  while (f1 > f2 + delta*lambda*transpose(gk)*pk)
    lambda = lambda*0.8;
    xk_lambda = xk + lambda*pk;
    f1 = fcn(xk_lambda(1), xk_lambda(2));
    endwhile
    
    res = lambda;
    
end

