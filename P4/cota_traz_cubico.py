from sympy import *
import numpy as np

def cota_traz_cubico(f, x_n):
    """
        Función que calcula la cota de error para el trazador
        cúbico natural
        input:
            f : función a la que se aproxima el trazador cúbico natural
            x_n : vector de preimágenes para un intervalo en f
        output:
            cota: valor de error máximo para una imagen del trazador aproximado
    """
    n = len(x_n)
    x = Symbol("x")
    f1 = simplify(f)
    dx = f1.diff(x)
    dxx = dx.diff(x)
    dxxx = dxx.diff(x)
    dxxxx = dxxx.diff(x)
    m = len(x_n)
    imagenes = []
    for i in range(0, m):
        xi = x_n[i]
        yi = float(dxxxx.subs(x,xi))
        imagenes.append(abs(yi))
    h_i = calc_hi(x_n)
    h = max(h_i)
    alpha = max(imagenes)
    cota = (5*h**4)*alpha/384
    return cota
    
def calc_hi(x_n):
    """
        Función para calcular los diferenciales entre las preimágenes
        input:
            x_n: vector de preimágenes
        output:
            h_i: vector de diferenciales
    """
    h_i = []
    n = len(x_n)
    for i in range(n-1):
        h = x_n[i+1] - x_n[i]
        h_i.append(h)
    return h_i
#--------------------Ejemplo presentación 9, diap 25-----------------------------#
x_n = [1, 1.05, 1.07, 1.1]
cota = cota_traz_cubico('3*x*exp(x) - 2*exp(x)', x_n)
print(cota)
