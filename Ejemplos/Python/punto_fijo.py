import math
import sympy


"""
Metodo iterativo del punto fijo, para solucion de ecuaciones no lineales.
Entradas: valor incial de iteracion, tolerancia minima del resultado y
funcion a evaluar en formato lambda.
Salidas:Valor aproximado y numero de iteraciones
"""
def punto_fijo(x0, tolerance, phi):
    nIterations = 0;
    x_k = x0;
    while abs(eval(phi)(x_k) - x_k) >= tolerance:
        x_k = eval(phi)(x_k)
        nIterations += 1
    return x_k, nIterations
