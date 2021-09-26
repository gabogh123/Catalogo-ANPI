import numpy as np
import sympy as sp
import math

def falsa_posicion(f,a,b,iterMax,tol):

    x = sp.Symbol('x')

    f1 = sp.sympify(f)

    fx0 = sp.N(f1.subs(x,a))
    fx1 = sp.N(f1.subs(x,b))

    error = tol+1
    k = 0

    while k < iterMax:

        if sp.N(f1.subs(x,a))*sp.N(f1.subs(x,b)) < 0: #Bolzano

            xk = b - (fx1*(b-a) / (fx1-a))

            if sp.N(f1.subs(x,a))*sp.N(f1.subs(x,xk))<0:

                b = xk

            else:

                a = xk

            error=abs(sp.N(f1.subs('x',xk)))

            if error < tol:
                break
        k+=1

    return [xk,k,error]


f = 'cos(x)-x'
a = 1/2
b = math.pi/4
tol = 10**-5
iterMax = 1000

y = falsa_posicion(f,a,b,iterMax,tol)
print(y)

            

    
