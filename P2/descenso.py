#Se importan las librerias necesarias para los calculos
import sympy as sp
import numpy as np
from sympy import solve
import math
import matplotlib.pyplot as plt


def variables_simbolicas(variables):
#Esta funcion aproxima retorna las variables simbolicas del programa
#Parametros de entrada: variables: string con las variables simbolicas a incluir en la funcion
    
#Parametros de salida: v_var: lista con las variables simbolicas del sistema
 
    n = len(variables) #Largo del string de entrada para la definicion del tamaño
    tam = np.arange(0,n,2) #Se toman los valores del string
    v_var = [] #Creacion de la lista de variables
    for i in tam:
        v_var.append(sp.Symbol(variables[i])) #Se agregan a la lista las variables simbolicas
    return v_var #Se retorna la lista de variables simbolicas



def gradiente(f,v_var):
    
#Esta funcion calcula el gradiente de una funcion
#Parametros de entrada: f: funcion a la cual se le calcula el gradiente
#                       v_var: lista de variables simbolicas                         
#Parametros de salida:  g : gradiente de la funcion

    
    n=len(v_var)
    g=[] #Se define el vector gradiente
    for i in np.arange(0,n,1): #Criterio de parada cuando llega al largo de la lista de variables simbolicas
        g.append(sp.diff(f,v_var[i])) # Se aplica la derivada de acuerdo a la variable en la que se encuentre el contador
    return g # Se retorna el gradiente




def descenso(f,v0,iterMax,tol):

#Esta funcion aproxima el vector de convergencia del problema de optimizacion mediante el método de descenso coordinado 
#Parametros de entrada: f: funcion a la cual se le calcula el gradiente
#                       v0: vector inicial
#                       iterMax: cantidad de iteraciones maximas del programa
#                       tol: tolerancia del programa
#Parametros de salida:  vk: vector aproximado
#                       error: cuota de error del programa



    x=sp.Symbol('x')
    y = sp.Symbol('y') #Se definen las variables simbolicas
    f1=sp.sympify(f) #Se pasa de string a funcion simbolica

    grad = gradiente(f1,v_var)   #Se calcula el gradiente de la funcion estudiada

    x0=v0[0] #Se obtiene el primer elemento del vector inicial (x inicial)
    y0=v0[1] #Se obtiene el segundo elemento del vector inicial (y inicial)

    k = 0 #Se inicializa el k de iteraciones
    
    error = tol + 1 #Se define el error del programa

#Se inicializan las variables de vector, x y y a utilizar
    vk = v0 
    xk = x0
    yk = y0
    

    while k<iterMax: #Primera condicion de parada: k debe ser menor a la cantidad maxima de iteraciones

        f2 = f1.subs('y',yk) # Se sustituye el yk en la funcion
               
        xk = float(solve(sp.diff(f2,x))[0]) #Se calcula el xk como la solucion del gradiente evaluado en x

        f3 = f1.subs('x',xk) # Se sustituye el xk en la funcion
            
        yk = float(solve(sp.diff(f3,y))[0]) #Se calcula el yk como la solucion del gradiente evaluado en y

        vk = [xk,yk] #Se redefine el vector de aproximacion


#       Calculo del error como la norma del vector

        error = float(math.sqrt( ((grad[0].subs('x',xk)).subs('y',yk))**2 + ((grad[1].subs('x',xk)).subs('y',yk))**2  ))

        k=k+1 # Se avanza a la siguiente iteracion

        if error<tol: #Segunda condicion de parada: el error debe ser mayor a la tolerancia
            break
        
         
    
   
    return [vk, error] #Se retorna el vector de aproximacion y el error
                
                
    

            
#Definicion de variables simbolicas
variables = 'x y'
v_var = variables_simbolicas(variables)

#Definicion de parametros de entrada

f = '(x-2)**2+(y+3)**2+x*y'
v0 = [1,1]
iterMax = 1000
tol = 10**-9


y = descenso(f,v0,iterMax,tol)
print(y)


