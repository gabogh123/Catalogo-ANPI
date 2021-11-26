import sympy as sp

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

def calculo_el(polinomio_ventas, polinomio_costos):
    #Funcion que realiza el calculo de la elasticidad
    x = sp.Symbol('x')
    #Primero se deriva cada uno de los polinomios
    dpv = sp.diff(polinomio_ventas, x, 1)
    dpc = sp.diff(polinomio_costos, x, 1)
    print("Derivada del Polinomio de Ventas de Mercancía")
    print(dpv)
    print("Derivada del Polinomio de Costos de Venta")
    print(dpc)
    #Se retorna la derivada de ventas/costos
    return

# Ejercicio polinomio grado 3 diap 33, presentacion 8
xv = [1, 2, 3, 4, 5, 6 ,7, 8]
yv_ventas = [19819634, 24411826, 32429311,38076878,43133089,51273750,65564512,67917477]
yv_costos = [12087296, 14802522, 19872418,22015823,24031723,28651153,36878941,37944221]

polinomio_ventas = dd_newton(xv, yv_ventas)
polinomio_costos = dd_newton(xv,yv_costos)

print("Polinomio de Ventas de Mercancía")
print(polinomio_ventas)

print("Polinomio de Costos de Venta")
print(polinomio_costos)

elasticidad_ventas_costos = calculo_el(polinomio_ventas,polinomio_costos)

print(elasticidad_ventas_costos)
