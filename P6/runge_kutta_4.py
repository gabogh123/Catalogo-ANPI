import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import lagrange as lg

def runge_kutta_4(f, ab, n):
    """
    Método de runge kutta de orden 4 para aproximar la solución a la ecuación
    diferencial de orden 1 de la forma dy/dx = f(x, y)
    Parámetros de entrada:
       f: función a solucionar f(x, y)
       ab: vector con el rango en que se pretende aproximar la función
       n: cantidad de puntos con los que se aproxima la solución
    Valores de salida, pares ordenados:
        xv: vector de puntos en x
        yv: vector de puntos en y
        polinom: polinomio que describe la aproximación de la
        solución
    """

    x = sp.Symbol('x')
    y = sp.Symbol('y')
    f1 = sp.sympify(f)
    a = ab[0]
    b = ab[1]
    h = (b-a)/(n-1)
    xv = [] #vector de preimágenes a usar de acuerdo con una distancia h
    yv = [1] #y0 = 1, condición inicial para la derivada
    for j in range(n):
        aux = a + j*h
        xv.append(aux)
        
    #calculo de las 4 etapas
    for i in range(n-1): 
        k1 = sp.N(f1.subs({x:xv[i], y:yv[i]}))
        k2 = sp.N(f1.subs({x:xv[i] + h/2, y:yv[i]+ h*k1/2}))
        k3 = sp.N(f1.subs({x:xv[i] + h/2, y:yv[i]+ h*k2/2}))
        k4 = sp.N(f1.subs({x:xv[i] + h, y:yv[i] + h*k3}))
        yv_siguiente = yv[i] + h*(k1 + 2*k2 + 2*k3 + k4)/6
        yv.append(yv_siguiente)
        
    #calculo del polinomio a partir de los puntos
    polinom = lg.lagrange(xv, yv) 

    print("Preimágenes \n", xv)
    print("Imágenes \n", yv)
    print("Polinomio graficado: \n", polinom)
    graficar(xv,yv)

def graficar(xv, yv):
    """
    Función para hacer la gráfica de la solución a la ecuación
    diferencial aproximada.
    
    """
    plt.rcParams.update({'font.size': 14})
    fig, graf = plt.subplots()
    graf.plot(xv, yv, 'b', marker='o', markerfacecolor='red', markersize=5)
    graf.set_xlabel('$x$')
    graf.set_ylabel('Solución  $f(x)$')
    graf.set_title('Método de Runge-Kutta 4')
    graf.grid(True)
    plt.show()

#Ejemplo presentación 12, diap 48
runge_kutta_4('-x*y + 4*x/y', [0, 1], 11)
