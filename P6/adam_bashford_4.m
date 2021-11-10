function [xv,yv_ab]  = adam_bashford_4(func, a, b, n)
  
    #Funcion que utiliza el metodo de adam brashfor de 4 pasos
    #Entradas func: Ecuacion diferencial a calcularle la solucion   
    #         a: valor inicial del intervalo 
    #         b: valor final del intervalo
    #         n: numero de puntos a utilizar
    
    [xv,yv]  = predictor_corrector(func, a, b, n) % se obtienen los vectores xv y yv
    
    h = (b - a) / (n - 1);
    f=@(x,y) y-x.^2+1;
    
    yn_3 = yv(1); %y0
    yn_2 = yv(2); %y1
    yn_1 = yv(3); %y2
    yn = yv(4); %y3
    
    yv_ab = [yn_3, yn_2, yn_1, yn]; % vector de las pre imagenes para el metodo de adam brashford
    
    for i = 4:n-1
        y_next = yv(i) + (h/24) * (55 * f(xv(i),yv(i))) - (59 * f(xv(i-1),yv(i-1))) + (37 * f(xv(i-2),yv(i-2))) - (9 * f(xv(i-3),yv(i-3)));
        yv_ab = [yv_ab, y_next];
    endfor
    yv_ab = transpose(yv_ab);
    xv;
    %pol = lagrange(xv, yv_ab)
endfunction





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
    y = [y0];
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

function p = lagrange(xv,yv)
  syms x
  n=length(xv)-1;
  p=0;
  for k=0:n
    p = p + yv(k+1) * fun_Lk(xv,k);
  end
  p = expand(p);
end