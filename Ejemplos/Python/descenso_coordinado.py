#Se importan las librerias necesarias para los calculos
import sympy as sp
import numpy as np
from sympy import solve
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


def variables_simbolicas(variables):
#Esta funcion aproxima retorna las variables simbolicas del programa
#Parametros de entrada: variables: string con las variables
#simbolicas a incluir en la funcion
    
#Parametros de salida: v_var: lista con las variables simbolicas del sistema
 
    n = len(variables) #Largo del string de entrada 
    tam = np.arange(0,n,2) #Se toman los valores del string
    v_var = [] #Creacion de la lista de variables
    for i in tam:
        #Se agregan a la lista las variables simbolicas
        v_var.append(sp.Symbol(variables[i])) 
    return v_var #Se retorna la lista de variables simbolicas



def gradiente(f,v_var):
    
#Esta funcion calcula el gradiente de una funcion
#Parametros de entrada: f: funcion a la cual se le calcula el gradiente
#                       v_var: lista de variables simbolicas                         
#Parametros de salida:  g : gradiente de la funcion

    
    n=len(v_var)
    g=[] #Se define el vector gradiente

    #Criterio de parada cuando llega al largo de la lista de variables simbolicas
    for i in np.arange(0,n,1):
    # Se aplica la derivada de acuerdo a la variable en la que se encuentre
        g.append(sp.diff(f,v_var[i])) 
    return g # Se retorna el gradiente




def coordinado(f,v0,iterMax,tol):

#Esta funcion aproxima el vector de convergencia del problema de optimizacion
#mediante el metodo de descenso coordinado 
#Parametros de entrada: f: funcion a la cual se le calcula el gradiente
#                       v0: vector inicial
#                       iterMax: cantidad de iteraciones maximas del programa
#                       tol: tolerancia del programa
#Parametros de salida:  vk: vector aproximado
#                       error: cuota de error del programa



    x=sp.Symbol('x')
    y = sp.Symbol('y') #Se definen las variables simbolicas
    f1=sp.sympify(f) #Se pasa de string a funcion simbolica
    er = [] #Vector de errores

    grad = gradiente(f1,v_var)   #Se calcula el gradiente de la funcion estudiada

    x0=v0[0] #Se obtiene el primer elemento del vector inicial (x inicial)
    y0=v0[1] #Se obtiene el segundo elemento del vector inicial (y inicial)

    k = 0 #Se inicializa el k de iteraciones
    
    error = tol + 1 #Se define el error del programa

#Se inicializan las variables de vector, x y y a utilizar
    vk = v0 
    xk = x0
    yk = y0
    
#Primera condicion de parada: k debe ser menor a la cantidad maxima de iteraciones
    while k<iterMax: 

        f2 = f1.subs('y',yk) # Se sustituye el yk en la funcion

        #Se calcula el xk como la solucion del gradiente evaluado en x
        xk = float(solve(sp.diff(f2,x))[0]) 

        f3 = f1.subs('x',xk) # Se sustituye el xk en la funcion

        #Se calcula el yk como la solucion del gradiente evaluado en y    
        yk = float(solve(sp.diff(f3,y))[0]) 

        vk = [xk,yk] #Se redefine el vector de aproximacion


#       Calculo del error como la norma del vector

        error = float(math.sqrt( ((grad[0].subs('x',xk)).subs('y',yk))**2 +
                                 ((grad[1].subs('x',xk)).subs('y',yk))**2  ))
        er.append(sp.N(error))
    
        k=k+1 # Se avanza a la siguiente iteracion

    #Segunda condicion de parada: el error debe ser mayor a la tolerancia
        if error<tol: 
            break
        
    plt.rcParams.update({'font.size': 14})
    ejex=np.arange(1,k+1,1)
    fig, graf=plt.subplots()     

    graf.plot(ejex,er,'b',marker='o',markerfacecolor='red',markersize=10)
    graf.set_xlabel('Iteraciones (k)')#Nombre del eje x
    graf.set_ylabel('$|f(x_k)|$') #Nombre eje y
    #Titulo grafica
    graf.set_title('Método del Descenso Coordinado');  
    graf.grid(True) #Mostrar grid
    plt.show() #Mostrar grafica
   
    return [vk,error] #Se retorna el vector de aproximacion y el error
                
                
    

            
#Definicion de variables simbolicas
variables = 'x y'
v_var = variables_simbolicas(variables)

#Definicion de parametros de entrada

f = '(x-2)**2+(y+3)**2+x*y'
v0 = [1,1]
iterMax = 1000
tol = 10**-9


y = coordinado(f,v0,iterMax,tol)
print(y)


