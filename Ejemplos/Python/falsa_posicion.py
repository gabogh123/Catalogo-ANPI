import math
import sympy


"""
Metodo iterativo de falsa posicion, para solucion de ecuaciones no lineales.
Entradas: limites de rango incial de iteracion, tolerancia minima del
resultado y funcion a evaluar en formato lambda.
Salidas:Valor aproximado y numero de iteraciones
"""
def falsa_posicion(a, b, tolerance, function):
    if eval(function)(a) * eval(function)(b) > 0:
        return 0, 0
    nIterations = 0
    while True:
        denominator = (eval(function)(b) - eval(function)(a))
        x_k = b - eval(function)(b) * (b - a) / denominator
        if abs(eval(function)(x_k)) < tolerance:
            break
        elif eval(function)(a) * eval(function)(x_k) < 0:
            b = x_k
        else:
            a = x_k
        nIterations += 1
    return x_k, nIterations
