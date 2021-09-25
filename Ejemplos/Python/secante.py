import math
import sympy


"""
Metodo iterativo de la secante, para solucion de ecuaciones no lineales.
Entradas: valores iniciales de iteracion, tolerancia minima del resultado y
funcion a evaluar en formato lambda.
Salidas:Valor aproximado y numero de iteraciones
"""
def secante(x0, x1, tolerance, function):
    x_k = x1
    x_k_1 = x0
    nIterations = 0
    while abs(eval(function)(x_k))  >= tolerance:
        temp = x_k
        numerator = (eval(function)(x_k) * (x_k - x_k_1))
        denominator = (eval(function)(x_k) - eval(function)(x_k_1))
        x_k = x_k - numerator / denominator
        x_k_1 = temp
        nIterations += 1
    return x_k, nIterations
