import math
import sympy as sp

def newton_raphson(f,x0,tol,iterMax):

    x = sp.Symbol('x')

    f1 = sp.sympify(f)

    df1 = sp.diff(f1,x)

    error = tol+1

    k = 0

    xk = x0

    while k < iterMax:

        n = sp.N(f1.subs(x,xk))
        d = sp.N(df1.subs(x,xk))

        k = k+1

        xk = xk - (n/d)

        error = abs(f1.subs(x,xk))

        if error<tol:
            break

    return [xk,k,error]



f = 'cos(2*x)**2 - x**2'
x0 = 3/4
tol = 10**-5
iterMax = 1000

y = newton_raphson(f,x0,tol,iterMax)
print(y)

        
        
