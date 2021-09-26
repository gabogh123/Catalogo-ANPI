import numpy as np
import sympy as sp

def biseccion(f,a,b,tol,iterMax):

    x = sp.Symbol('x')
    f1=sp.sympify(f)
    error = tol + 1
    k = 0

    if sp.N(f1.subs('x',a)) * sp.N(f1.subs('x',b)) < 0:

        for k in range (1,iterMax):

            xk = (a+b)/2


            if sp.N(f1.subs('x',a)) * sp.N(f1.subs('x',xk)) < 0:

                b = xk

            else:

                a = xk

            error = abs(sp.N(f1.subs('x',xk)))

            if error < tol:
                break

        

    else:
        print('El resutlado inicial no cumple el teorema de bolzano')

    return [xk,k,error]    


f = 'exp(x) -2*x - 10'
a = 2
b = 4
tol = 10**-5
iterMax = 1000

y = biseccion(f,a,b,tol,iterMax)
print(y)
 
        

        
