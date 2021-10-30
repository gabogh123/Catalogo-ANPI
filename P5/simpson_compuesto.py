# regla compuesta de Simpson y cota de error
import numpy
import numpy as np
import sys
from sympy import *
import sympy as sp
from numpy import log as ln
from scipy import optimize

# funcion que aproxima la integral definida utilizando la regla de Simpson
# entradas: una funcion, numero de puntos y un intervalo
# salidas: aproximacion y cota de error

def simpson_compuesto(fEntrada, mPuntos, ab):
    # se cambia la funcion de entrada a simbolico
    fs = sp.sympify(fEntrada)
    print("Funcion a utilizar: " , end ="")
    print(fs)
    print("Numero de puntos: ", end ="")
    print(mPuntos)
    print("Dentro del intervalo: ", end ="")
    print(ab)

    # se verifica que m sea impar
    # si m es par se termina la funcion
    if mPuntos % 2 == 0:
        return print('Error, la cantidad de puntos tiene que ser impar')

    # si m es impar, se procede con el calculo
    else:
        # se saca el largo solo para obtener b
        n = len(ab)
        # se obtienen los valores extremos del intervalo 
        a = ab[0] # extremo minimo
        b = ab[n-1] # extremo maximo
        # se averigua h
        h = (b-a)/(m-1)
        # se conoce h, entonces calculamos h/3 -> c
        c = h/3

        # para averiguar los puntos xi
        lista = [a] # inicializo una lista con el elem a ya incluido
        iAux = a # punto xi
        # se averiguan los puntos segun h
        for i in range(0,mPuntos-1):
            iAux = iAux + h
            lista.append(iAux)
        lista # lista con todos los puntos dentro del intervalo

        # para calcular la suma
        fPar = 0
        fImpar = 0
        # se averigua el largo
        n = len(lista)
        # se evaluan las funciones si son par o impar
        # excluyendo los extremos
        for i in range(1, n-1):
            if i%2 == 0:
                fPar = fPar + fs.subs('x',float(lista[i]))
            else:
                fImpar = fImpar + fs.subs('x',float(lista[i]))
        # evaluamos la funcion en a y b
        f_a = fs.subs('x',float(a))
        f_b = fs.subs('x',float(b))

        # calculamos la aproximacion
        aproxSimpsonComp = c * (f_a + (2 * fPar) + (4 * fImpar) + f_b)
        print("Aproximacion con Regla de Simpson Compuesto: " , end ="")
        print(aproxSimpsonComp)

        ######### para el calculo de la cota de error ############

        # se calcula la cuarta derivada de la funcion de entrada
        fs4 = sp.diff(fs,'x', 4)
        # se crea una funcion auxiliar negativa de la derivada
        faux = -(abs(fs4))
        # se pasa la funcion que entienda Scipy
        fs4numeric = sp.lambdify('x',faux)
        # se calcula el valor en x donde la funcion derivada es max
        x_max = optimize.fminbound(fs4numeric, a, b)

        # se pasa la funcion derivada a numeral
        fs4aux = sp.lambdify('x',abs(fs4))

        # se calcula la cota de error
        cota_de_error = (((b-a)*(h**4))/180)*fs4aux(x_max)
        print("Cota del error: " , end ="")
        print(cota_de_error)

funcion = 'ln(x)' # funcion a utilizar
m = 7 # numero de puntos
ab = [2,5] # intervalo de [2,3]
simpson_compuesto(funcion,m,ab)
