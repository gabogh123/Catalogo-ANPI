import sympy as sp
import numpy as np
from sympy import solve


def variables_simbolicas(variables):
    n = len(variables)
    tam = np.arange(0,n,2)
    v_var = []
    for i in tam:
        v_var.append(sp.Symbol(variables[i]))
    return v_var

def gradiente(f,v_var):
    n=len(v_var)
    g=[]
    for i in np.arange(0,n,1):
        g.append(sp.diff(f,v_var[i]))
    return g




def descenso():

    x=sp.Symbol('x')
    y = sp.Symbol('y')
    f1=sp.sympify(f)

    grad = gradiente(f1,v_var)
    
    error = tol + 1
    vk = v0
    xk = x0
    yk = y0
    

    #while k<iterMax and error>tol:

    f2 = f1.subs('y',v0[1])
           
    xk = solve(sp.diff(f2,x))[0]

    f3 = f1.subs('x',xk)
        
    yk = solve(sp.diff(f3,y))[0]

    vk = [xk,yk]

    
 
    
   
    print(error)
                
                
    

            

variables = 'x y'
v_var = variables_simbolicas(variables)
x0=1
y0=1

f = '(x-2)**2+(y+3)**2+x*y'
v0 = [x0,y0]
iterMax = 1000
tol = 10**-9

k = 0

