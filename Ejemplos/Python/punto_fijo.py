import sympy as sp
import numpy as np
import math

def punto_fijo(f,x0,tol,iterMax):

    x = sp.Symbol('x')
    f1 = sp.sympify(f)

    k = 0
    error = tol + 1

    xk = x0

    while k<iterMax:

        xk2 = sp.N(f1.subs(x,xk))

        xk = xk2

        error = abs(xk2 - xk)

        k += 1

        if error<tol:
            break

    return [xk2,k,error]



f = 'log(2*x+1)'
x0 = 7
tol = 10**-5
iterMax = 1000


y = punto_fijo(f,x0,tol,iterMax)
print(y)
    
