import sympy as sp
import numpy as np
import matplotlib.pyplot as plt


def dd_newton(xv, yv):
    # Funcion que aproxima el polinomio de interpolación
    # Parametro de entrada:
    #    xv: arreglo de preimagenes
    #    yv: arreglo de imagenes de los puntos
    # Parametros de salida:
    #    polinomio: polinomio simbolico generado

    # Definición de constantes
    x = sp.Symbol('x')
    n = len(xv)
    pol = 0  # Inicialización del polinomio

    for j in range(n):

        # Se llama a la función auxiliar de diferencias divididas
        val = dd_newton_aux(0, j, xv, yv)

        pol_var = 1
        for i in range(j):
            pol_var *= (x - xv[i])

        pol = pol + val * pol_var

    pol_expand = sp.expand(pol)
    return pol_expand


def dd_newton_aux(i, j, xv, yv):
    # Función que calcula los coeficientes  de los sumandos
    #       de diferencias divididas

    if i == j:
        return yv[i]

    else:
        n = dd_newton_aux(i + 1, j, xv, yv) - dd_newton_aux(i, j - 1, xv, yv)
        d = xv[j] - xv[i]
        return n / d

def calculo_el(polinomio_ventas, polinomio_costos, xv):
    """
    Funcion que realiza el calculo de la elasticidad
    se calcula la elasticidad según la formula e =  C*dV
                                                   -----
                                                    V*dC
    Entradas:
        polinomio_ventas, polinomio_costos: funcion aproximada con Interpolación de Newton
        xv: vector de preimagenes
    Salida:
        elasticidad: vector de imagenes para la función de elasticidad
        
    """
    
    n = len(xv)
    x = sp.Symbol('x')

    elasticidad = []
    costos = eval_fun(xv, polinomio_costos)
    ventas = eval_fun(xv, polinomio_ventas)
    #Primero se deriva cada uno de los polinomios
    dv = sp.diff(polinomio_ventas, x)
    dc = sp.diff(polinomio_costos, x)
    dx_ventas = eval_fun(xv, dv)
    dx_costos = eval_fun(xv, dc)


    for i in range(n):
        actual = costos[i]*dx_ventas[i]/(ventas[i]*dx_costos[i])
        elasticidad.append(actual)
    return elasticidad

def eval_fun(xv, fun):
    """
    Funcion que evalua las preimágenes de una funcion dada
    entradas
        xv: vector de preimágenes
        fun: funcion simbolica a evaluar
    salida:
        yv = vector de imagenes
    
    """
    x = sp.Symbol('x')
    yv = []
    n = len(xv)
    for i in range(n):
        y = sp.N(fun.subs(x, xv[i]))
        yv.append(y)
    return yv
    
def plot_datos(xv, yv, nombre):
    etiqueta_y = "$" + nombre + "$"
    """
    Funcion que se encarga de graficar el comportamiento de
    y = f(x)
    :param xv: vector de puntos independientes
    :param yv: vector de imágenes de xv
    """
    plt.rcParams.update({'font.size': 10})
    ejex = xv
    fig, graf = plt.subplots()
    #Se plotean la iteracion y el error
    graf.plot(ejex,yv,'b',marker='o',markerfacecolor='red',markersize=5)
    graf.set_xlabel('($x$)')#Nombre del eje x
    graf.set_ylabel(etiqueta_y) #Nombre eje y
    #Titulo grafica
    graf.set_title("Grafica de " + nombre);
    graf.grid(True) #Mostrar grid
    plt.show() #Mostrar grafica

#Calculando los polinomios de ventas y costos a partir de valores muestra  
xv = [1, 2, 3, 4, 5, 6 ,7, 8]
yv_ventas = [19819634, 24411826, 32429311,38076878,43133089,51273750,65564512,67917477]
yv_costos = [12087296, 14802522, 19872418,22015823,24031723,28651153,36878941,37944221]

polinomio_ventas = dd_newton(xv, yv_ventas)
polinomio_costos = dd_newton(xv,yv_costos)

print("Polinomio de Ventas de Mercancía\n")
print(polinomio_ventas)

print("\nPolinomio de Costos de Venta")
print(polinomio_costos)

#Calculando el polinomio de elasticidad para un rango de muestras
xv = [0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
elasticidad = calculo_el(polinomio_ventas,polinomio_costos, xv)
polinomio_elasticidad = dd_newton(xv,elasticidad)
print("\nPolinomio de elasticidad")
print(polinomio_elasticidad)

#Graficar funciones de ventas, costo y elasticidad para un rango de puntos
xt = np.linspace(1, 8, 35) #usando 35 puntos para ventas y costos como en el documento de la referencia
xv = np.linspace(0, 10, 50) #usando 50 puntos para elasticidad
yventas = eval_fun(xt, polinomio_ventas)
ycostos = eval_fun(xt, polinomio_costos)
yelasticidad = eval_fun(xv, polinomio_elasticidad)
plot_datos(xv, yelasticidad, "elasticidad")
plot_datos(xt, yventas, "ventas")
plot_datos(xt, ycostos, "costos")

