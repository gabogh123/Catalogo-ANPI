import sympy as sp


def dd_newton(xv, yv):

#Funcion que aproxima el polinomio de interpolación
#Parametro de entrada:
#    xv: arreglo de preimagenes
#    yv: arreglo de imagenes de los puntos
#Parametros de salida:
#    polinomio: polinomio simbolico generado

    #Definición de constantes
    x = sp.Symbol('x')
    n = len(xv)
    pol = 0 #Inicialización del polinomio

    for j in range(n):
        
        #Se llama a la función auxiliar de diferencias divididas
        val = dd_newton_aux(0, j, xv, yv) 


        pol_var = 1
        for i in range (j):
            pol_var *= (x-xv[i])

        pol = pol + val * pol_var

    pol_expand = sp.expand(pol)
    return pol_expand

def dd_newton_aux(i, j, xv, yv):

#Función que calcula los coeficientes  de los sumandos
#       de diferencias divididas

    if i == j:
        return yv[i]

    else:
        n = dd_newton_aux(i+1, j, xv, yv) - dd_newton_aux(i, j-1, xv, yv)
        d = xv[j]-xv[i]
        return n/d


#Ejercicio polinomio grado 3 diap 33, presentacion 8
xv = [2, 3, 5, 6]
yv = [2/3, 1, -1, 0]
polinimio = dd_newton(xv, yv)
print(polinimio)
