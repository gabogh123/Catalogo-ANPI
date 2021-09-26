import numpy as np
import sympy as sp

def secante(f,x0,x1,tol,iterMax):

    x = sp.Symbol('x')
    f1 = sp.sympify(f)

    k = 0
    error = tol + 1

    while k < iterMax:

        fx0 = sp.N(f1.subs(x,x0))
        fx1 = sp.N(f1.subs(x,x1))

        xk = x1 - fx1*((x1-x0) / (fx1 - fx0))

        x0 = x1
        x1 = xk

        error = abs(f1.subs(x,xk))

        if error<tol:
            break

        k = k+1

    return [xk , k , error]


f = 'exp(-x**2)-x'
x0 = 0
x1 = 1
tol = 10**-5
iterMax = 1000

y = secante(f,x0,x1,tol,iterMax)
print(y)
