import sympy as sp

def funcion_Lk(xv,k):


    x = sp.Symbol('x')
    m = len(xv)
    Lk = 1

    for i in range(m):

        if i != k:

            n = x-xv[i]
            d = xv[k] - xv[i]

            Lk = Lk* (n/d)

    exp = sp.expand(Lk)

    print("\nPolinomio iteración  " + str(k))
    print(exp)

    return exp


def lagrange(xv,yv):


    x = sp.Symbol('x')
    m = len(xv)
    p = 0

    for k in range(m):
        p = p+yv[k]* funcion_Lk(xv,k)

    polinomio = sp.expand(p)
    return polinomio


xv = [-2,0,1]
yv = [0,1,-1]
prueba = lagrange(xv,yv)
print(prueba)
