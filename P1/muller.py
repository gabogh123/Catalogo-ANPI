#Se importan las bibliotecas necesarias para los cálculos
#Python v 3.6
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

def muller(f, x0, x1, x2, tol, iterMax):
    """
    Esta función aproxima la solución de la ecuación f(x) = 0
    Sintaxis: y = muller(f, x0, x1, x2, tol, iterMax)
    Parámetros de entrada:
        f = string que representa la función f
        x0, x1, x2 = numeros enteros, valores iniciales establecidos para 
        obtener la ecuación cuadrática según Müller
        tol = número entero que representa la tolerancia para la condición de parada
        iterMax = numero entero, es la cantidad de iteraciones permitida

    Parametrós de Salida:
        x3 = aproximación para la solución de la ecuación f(x) = 0
        k = entero, cantidad de iteraciones realizadas
        error = grado de error en el cálculo de la solución para k iteraciones    
    """
    f1 = sp.sympify(f) #conversión del string f a simbólico
    er = []            #arreglo para almacenar el error por cada iteración
    error = tol + 1    #grado de error para comparar con tolerancia
    k = 0
    x3 = 0
    while k <= iterMax:
        """
        Según el método de Müller a partir de los valores iniciales
        se obtiene un sistema de ecuaciones en donde se resuelve para los
        coeficientes a, b y c, siendo que:
        a =  (x1 - x2)[f(x0) - f(x2)] - (x0 - x2)[f(x1) - f(x2)]
             ___________________________________________________
                        (x0 - x1)(x0 - x2)(x1 - x2)

        b =  (x0 - x2)^2 [f(x1) - f(x2)] - (x1 - x2)^2 [f(x0) - f(x2)]
             _________________________________________________________
                        (x0 - x1)(x0 - x2)(x1 - x2)

        c = f(x2)
        """
        k += 1
        #Variables de las diferencias entre los valores iniciales
        n0 = x0 - x1
        n1 = x0 - x2
        n2 = x1 - x2
        #Denominador para las soluciones de a y b
        deno = n0*n1*n2
        #Valores de restas en los numeradores de las soluciones para a y b
        h12 = sp.N(f1.subs('x', x1)) - sp.N(f1.subs('x', x2))
        h02 = sp.N(f1.subs('x', x0)) - sp.N(f1.subs('x', x2))
        #Coeficientes de la ecuación p(x) = a(x - x2)^2 + b(x - x2) + c 
        a  = (n2*h02 - n1*h12)/deno
        b = (n1**2*h12 - n2**2*h02)/deno
        c = sp.N(f1.subs('x', x2))   
        """
        Cálculo del discriminante  para la iteración actual 
                     r = x2 - 2c/[b + sgn(b)sqrt(b^2 -4ac)]                       
        """       
        d = (complex(b**2-4*a*c))**0.5  
        #validación del signo en el discriminante y su proximidad al valor x2
        if abs(b-d) < abs(b+d):
            aprox = 2*c/(b+d)
        else:
            aprox = 2*c/(b-d)
        x3 = x2 - aprox
        error = abs(sp.N(f1.subs('x', x3)) - sp.N(f1.subs('x', x2)))
        er.append(sp.N(error))
        if error < tol:
            break
        x0 = x1
        x1 = x2
        x2 = x3
    #Llamadas a métodos para el gráfico de iteraciones vs error     
    plt.rcParams.update({'font.size': 14})
    ejex=np.arange(1,k+1,1)
    fig, graf=plt.subplots()
    graf.plot(ejex,er,'b',marker='o',markerfacecolor='red',markersize=6)
    graf.set_xlabel('Iteraciones ($k$)')
    graf.set_ylabel('$|f(x3)-f(x2)|$')
    graf.set_title('Método de Müller (Iteraciones vs Error)')
    graf.grid(True)
    plt.show()

    return [x3, k, error]

#Definición de parámetros de salida. Ejemplo presentación 1, pág 107
f = 'sin(x) - x/2'
x0 = 2
x1 = 2.2
x2 = 1.8
tol = 10**-5
iterMax = 1000
y = muller(f, x0, x1, x2, tol, iterMax)
print(y)    


