from metodos_p2 import *

f='cos(2*x)-x'
x0=0.739085
tol=10**-5
iterMax=500


def prueba_newton_H_m1():

    print("Método 1: H(x) = 1")
    
    print("f = cos(2x)-x")

    y = newton_H_m1(f,x0,tol,iterMax)
    return y
        

def prueba_newton_H_m2():

    print("Método 2: H(x) = 1/(1+Bu(x))")
    
    print("f = cos(2x)-x")

    y = newton_H_m2(f,x0,tol,iterMax)
    print(y)



def prueba_newton_G_m1():

    
    
    print("Método 3: G(x) = 2/(2-w(x))")
    
    print("f = cos(2x)-x")

    y = newton_G_m1(f,x0,tol,iterMax)
    print(y)


def prueba_newton_G_m2():

    
    print("Método 4: G(x) = 1+ (w(x)/2")
    
    print("f = cos(2x)-x")

    y = newton_G_m2(f,x0,tol,iterMax)
    print(y)
