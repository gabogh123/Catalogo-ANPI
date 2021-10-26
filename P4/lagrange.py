import sympy as sp

def lagrange_aux(xv,k):

#Calculo de factores de variable de Lagrange
#Entradas: xv: Vector de preimagines
#           k: iteracion actual

    x = sp.Symbol('x')
    m = len(xv)
    Lk = 1

    for i in range(m):

        if i != k:

            n = x-xv[i]
            d = xv[k] - xv[i]

            Lk = Lk* (n/d)

    exp = sp.expand(Lk)



    return exp


def lagrange(xv,yv):
#Función que aproxima el polinomio de interpolacion
#Entradas: xv: Arreglo de preimagenes
#          yv: Arreglo de puntos considerados en el polinomio

#Salida: pol: polinomio simbolico a partir de diferencias divididas de Newton

    x = sp.Symbol('x')
    m = len(xv)
    p = 0


#Se forma el polinomio con la funcion auxiliar
    for k in range(m):
        p = p+yv[k]* lagrange_aux(xv,k)

#Distribución de factores
    polinomio = sp.expand(p)

    print("Polinomio de interpolación: ")
    print(str(polinomio))
5

xv = [-2,0,1]
yv = [0,1,-1]
prueba = lagrange(xv,yv)

