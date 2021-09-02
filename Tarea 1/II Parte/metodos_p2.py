#Importacion de librerias necesarias
import numpy as np
import sympy as sp

def newton_H_m1(f,x0,tol,iterMax):

#Esta funcion aproxima la solucion de la ecuacion f(x)=0
#Parametros de entrada: f = string que representa la funcion f,
#                       x0 = valor inicial
#                       tol = numero que representa la tolerancia
#                       iterMax = cantidad de iteraciones maximas
#Parametros de salida:  xk = aproximacion del cero o solucion de la funcion f                    
#                       error = grado de error del calculo --> |f(xk)|


#Se inicializan los parametros y se convierte la expresion de string a funcion
    x=sp.Symbol('x')

    
    f1=sp.sympify(f)
    df1 = sp.diff(f1,x)

    er=[]
    err=tol+1
    k=0
    xk=x0

    
#Funcion h de la derivacion del metodo de newton-raphson     
    h = 1
    h1 = sp.sympify(h)


    while err>tol and k<iterMax: #Criterios de parada del programa
        k=k+1
        
        n=sp.N(f1.subs(x,xk)) #Sustitucion de xk en la funcion f
        d=sp.N(df1.subs(x,xk)) #Sustitucion de xk en la derivada de f
        
        xk=xk-(h*n/d) #Recalculo de xk

        
        err=abs(f1.subs(x,xk)) #Calculo del error
        

#Se imprimen la aproximacion y el error
    print( "Aproximación: " + str(xk))
    print ("Error: " + str(err))


        
    

####################################################################################################################################

def newton_H_m2(f,x0,tol,iterMax):

#Esta funcion aproxima la solucion de la ecuacion f(x)=0
#Parametros de entrada: f = string que representa la funcion f,
#                       x0 = valor inicial
#                       tol = numero que representa la tolerancia
#                       iterMax = cantidad de iteraciones maximas
#Parametros de salida:  xk = aproximacion del cero o solucion de la funcion f                    
#                       error = grado de error del calculo --> |f(xk)|


#Se inicializan los parametros y se convierte la expresion de string a funcion

    x = sp.Symbol('x')
    B=2 #Constante beta

    f1 = sp.sympify(f) 
    df1 = sp.diff(f1,x)

    err = tol+1
    k=0
    xk = x0

    while err>tol and k<iterMax: #Criterios de parada del programa


        h = 1/(1+B*(sp.N(f1.subs(x,xk))/sp.N(df1.subs(x,xk)))) #Funcion h de la derivacion del metodo de newton-raphson   

        n = sp.N(f1.subs(x,xk))#Sustitucion de xk en la funcion f
        d = sp.N(df1.subs(x,xk))#Sustitucion de xk en la derivada de f

        xk = xk - (h*(n/d))#Recalculo de xk

        err = abs(f1.subs(x,xk))#Calculo del error

        k=k+1
        
    #Se imprimen la aproximacion y el error
    return ['Aproximación: ']+[xk]+ [   'Error: ']+ [err]




###########################################################################################################################


def newton_G_m1(f,x0,tol,iterMax):

#Esta funcion aproxima la solucion de la ecuacion f(x)=0
#Parametros de entrada: f = string que representa la funcion f,
#                       x0 = valor inicial
#                       tol = numero que representa la tolerancia
#                       iterMax = cantidad de iteraciones maximas
#Parametros de salida:  xk = aproximacion del cero o solucion de la funcion f                    
#                       error = grado de error del calculo --> |f(xk)|


#Se inicializan los parametros y se convierte la expresion de string a funcion

    x = sp.Symbol('x')

    f1 = sp.sympify(f)
    df1 = sp.diff(f1,x)
    df2 = sp.diff(df1,x)

    err = tol+1
    k = 0
    xk = x0

    while err>tol and k<iterMax: #Criterios de parada del programa

        g = 2/( 2 - ((sp.N(f1.subs(x,xk)))* (sp.N(df2.subs(x,xk))) / (sp.N(df1.subs(x,xk)))**2 )) #Funcion g de la derivacion del metodo de newton-raphson

        n = sp.N(f1.subs(x,xk))#Sustitucion de xk en la funcion f
        d = sp.N(df1.subs(x,xk))#Sustitucion de xk en la derivada de f

        
        xk = xk - (g*(n/d))#Recalculo de xk

        err = abs(f1.subs(x,xk))#Calculo del error

        
        k=k+1
    #Se imprimen la aproximacion y el error
        
    return ['Aproximación: ']+[xk]+ [   'Error: ']+ [err]
    


#############################################################################################################################


def newton_G_m2(f,x0,tol,iterMax):

#Esta funcion aproxima la solucion de la ecuacion f(x)=0
#Parametros de entrada: f = string que representa la funcion f,
#                       x0 = valor inicial
#                       tol = numero que representa la tolerancia
#                       iterMax = cantidad de iteraciones maximas
#Parametros de salida:  xk = aproximacion del cero o solucion de la funcion f                    
#                       error = grado de error del calculo --> |f(xk)|


#Se inicializan los parametros y se convierte la expresion de string a funcion

    
    x = sp.Symbol('x')

    f1 = sp.sympify(f)
    df1 = sp.diff(f1,x)
    df2 = sp.diff(df1,x)

    err = tol+1
    k = 0
    xk = x0

    while err>tol and k<iterMax: #Criterios de parada del programa
    

        g = 1 + ( ((sp.N(f1.subs(x,xk)))* (sp.N(df2.subs(x,xk))) / (sp.N(df1.subs(x,xk)))**2 ) /2 ) #Funcion g de la derivacion del metodo de newton-raphson

        n = sp.N(f1.subs(x,xk))#Sustitucion de xk en la funcion f
        d = sp.N(df1.subs(x,xk))#Sustitucion de xk en la derivada de f

        xk = xk - (g*(n/d))#Recalculo de xk

        err = abs(f1.subs(x,xk))#Calculo del error
        
        k=k+1

        #Se imprimen la aproximacion y el error

    return ['Aproximación: ']+[xk]+ [   'Error: ']+ [err]



