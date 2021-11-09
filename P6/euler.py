def euler(f, interv, y0, n):

    import sympy as sp
    import numpy as np
    import matplotlib.pyplot as plt
    import lagrange as lg
    
  
#Funcion que calcula la solucion de la ecuacion diferencial
#y'(x) = f(x, y(x))
#por medio del metodo de euler.
#Entradas f: Ecuacion diferencial a calcularle la solucion   
#         interv: array con valor inicial y final de intervalo 
#         y0: imagen del valor inicial
#         n: numero de puntos a utilizar


#Inicialización de variables
    x = sp.Symbol('x')
    y = sp.Symbol('y')
    f1 = sp.sympify(f)

#Calculo del paso h
    a = interv[0]
    b = interv[1]
    h = (b-a)/n

#Calculo las preimaegnes
    xv = []
    for j in range(n+1):
        aux = a + j*h
        xv.append(aux)

#Calculo de las imagenes
    yv = [y0]
    for i in range(n):
        aux = yv[i] + h * sp.N(f1.subs({x:xv[i], y:yv[i]}))
        yv.append(aux)

#Calculo polinomio de interpolacion
    pol = lg.lagrange(xv, yv)

    # grafica
    plt.rcParams.update({'font.size': 14})
    fig, graf = plt.subplots()
    graf.plot(xv, yv, 'b', marker='o', markerfacecolor='red', markersize=10)
    graf.set_xlabel('Solucion  ($y(x)$)')
    graf.set_ylabel('$x$')
    graf.set_title('Metodo de Euler')
    graf.grid(True)
    plt.show()

    print( "Preimagenes: " + str(xv) + "\n\n"

            + " Imagenes: " + str(yv) + "\n\n"

           + "Polinomio: " + str(pol))

# Ejemplo diapositiva 13, presetnación 12

f = 'y-x^2+1'
interv = [0, 2]
n = 11
y0 = 0.5

z = euler(f, interv, y0, n)


