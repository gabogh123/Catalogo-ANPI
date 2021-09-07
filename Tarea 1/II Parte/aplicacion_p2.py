#Importacion de las librerias
import numpy as np
import sympy as sp
from sympy import *
import math


def aplicacion_p2():
#Esta funcion aproxima la distancia en redes de sensores inalámbricos
    
#Parametros de salida:  dk2 = aproximacion de la distancia de los sensores                   
#                       error = grado de error del calculo --> |f(xk)|

    d = symbols('d') #Se hace simbolica la variable d

# Inicialización de parámetros
    d0 = 1 
    tol = 10**-5
    error = tol+1
    iterMax = 100

    dk = d0

    k = 0

#Función g dependiente de d

    g = ( 200*acos(d/20)) - (d*sqrt(100-((d**2)/4)))

# Función f dependiendo de d y de la función g 

    f = ((log(7/d,10))/( 0.01*log(10)))  + ((d*(6-d)) / (((g**2) / (2*1*(40/log(10))**2)) * ( (1/g) + (1/(pi*100)) ))) #Listo   

    df = sp.diff(f,d) #Derivada de la función f

    while  k<iterMax: #Condición de parada 

        dk2 = dk -((f.subs(d,dk))/(df.subs(d,dk))) #Redefiición de dk newton-raphson

        error = abs((dk2-dk)/dk2) #Calculo del error con el valor absoluto de las iteraciones

        dk = float(dk2) # Redefinición del dk 

        k=k+1 #Se avanza en la iteración

    
#Se retornar la aproximación y el error
    return [float(dk2),float(error)] 

        
#Se imprime el resultado de la función
print(aplicacion_p2())   



    
    

    
