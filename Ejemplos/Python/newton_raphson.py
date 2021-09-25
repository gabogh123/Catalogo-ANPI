import math
import sympy


"""
Metodo iterativo de Newton-Raphson, para solucion de ecuaciones no lineales.
Entradas: valor inicial de iteracion, tolerancia minima del resultado y
funcion a evaluar en formato lambda.
Salidas:Valor aproximado y numero de iteraciones
"""
def newton_raphson(x0, tolerance, function):
    x = sympy.symbols('x')
    nIterations = 0
    x_k = x0
    while abs(eval(function)(x_k))  >= tolerance:
        deriv = sympy.diff(eval(function)(x), x)
        x_k = x_k - eval(function)(x_k) / deriv.evalf(subs={x: x_k})
        nIterations += 1
    return x_k, nIterations
