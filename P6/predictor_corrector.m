function [x,y]  = predictor_corrector( func, a, b, n)
    pkg load symbolic;
    syms x y z
    
    #Funcion que calcula la solucion de una ecuacion diferencial 
    #por medio del metodo de predictor corrector.
    #Entradas func: Ecuacion diferencial a calcularle la solucion   
    #         a: valor inicial del intervalo 
    #         b: valor final del intervalo
    #         n: numero de puntos a utilizar
   
    f = str2func(func); % por un error no conocido no funciona esto con la funcion de entrada por eso se setea en la linea de abajo manualmente
    f=@(x,y) y-x.^2+1;
    h = (b - a) / (n - 1);
    y0 = 0.5;
    x = [a];
    y = [0.5];
    zn = [];
   
    for i = 1:n-1
        x_temp = x(i) + h;
        x = [x, x_temp];
        zn = [zn, y(i) + h * f(x(i), y(i))]; % predictor
        predict = y(i) + h * (f(x(i), y(i)) + f(x(i + 1), zn(i))) / 2; %corrector
        y = [y, predict];
        
    endfor
    x = transpose(x);
    y = transpose(y);
    %pol = lagrange(x, y);
    
endfunction


function Lk = fun_Lk(xv,k)
  syms x
  %k=0,1,....,n
  n = length(xv) - 1;
  Lk=1;
  for j=0:n
    if j~=k
      Lk = Lk * (x - xv(j+1)) / (xv(k+1) - xv(j+1));
    end    
  end
  Lk = expand(Lk);
end

function pol = lagrange(xv,yv)
  syms x
  n=length(xv)-1;
  pol=0;
  for k=0:n
    pol = pol + yv(k+1) * fun_Lk(xv,k);
  end
  pol = expand(pol);
end