# regla de Simpson y cota de error
import numpy
import numpy as np
import sys
from sympy import *
import sympy as sp
from numpy import log as ln
from scipy import optimize

# funcion que aproxima la integral definida utilizando la regla de Simpson
# entradas: una funcion, un intervalo
# salidas: aproximacion y cota de error

def simpson(fEntrada, ab):
    # se cambia la funcion de entrada a simbolico
    fs = sp.sympify(fEntrada)
    print("Funcion a utilizar: " , end ="")
    print(fs)
    print("Dentro del intervalo: ", end ="")
    print(ab)

    # se saca el largo solo para obtener b
    n = len(ab)
    # se obtienen los valores extremos del intervalo 
    a = ab[0]
    b = ab[n-1]
    print('a',a)
    print('b',b)
    # se conoce a y b, entonces calculamos (a+b)/2 -> c
    c = (a+b)/2

    #calculando (b-a)/6 -> d
    d = (b-a)/6
    
    # evaluamos la funcion en a, b y c
    f_a = fs.subs('x',a)
    f_b = fs.subs('x',b)
    f_c = fs.subs('x',c)
    print('fa',f_a)
    print('fb',f_b)
    print('fc',f_c)

    # calculamos la aproximacion
    aproxSimpson = d * (f_a + (4 * f_c) + f_b)
    print("Aproximacion con Regla de Simpson: " , end ="")
    print(aproxSimpson)

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
    cota_de_error = (((b-a)**5)/2880)*fs4aux(x_max)
    print("Cota del error: " , end ="")
    print(cota_de_error)
    
funcion = 'ln(x)' # funcion a utilizar
ab = [2,5] # intervalo de [2,3]
simpson(funcion, ab)
