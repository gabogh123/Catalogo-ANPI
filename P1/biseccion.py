def biseccion(f,a,b,tol,iterMax):

#Esta funcion aproxima la solucion de la ecuacion f(x)=0
#Parametros de entrada: f = string que representa la funcion f,
#                       a,b = extremos del intervalo [a,b]
#                       tol = numero que representa la tolerancia
#                       iterMax = cantidad de iteraciones maximas
#Parametros de salida:  xk = aproximacion del cero o solucion de la funcion f
#                       k = numero de iteraciones realizadas
#                       error = grado de error del calculo --> |f(xk)|


#Se importan las librerias necesarias para los calculos
    import numpy as np
    import sympy as sp
    import matplotlib.pyplot as plt
    

#Se inicializan los parametros y se convierte la expresion de string a funcion
    x=sp.Symbol('x')
    f1=sp.sympify(f)
    er=[]
    error=tol+1
    k=0

    if sp.N(f1.subs('x',a)) * sp.N(f1.subs('x',b))<0:

        for k in range(1,iterMax):

            xk=(a+b)/2 #Se crea el xk a la mitad del intervalo

            #Se cumple el teorema de Bolzano para el  subintervalo [a,xk]

            if sp.N(f1.subs('x',a))*sp.N(f1.subs('x',xk))<0: 
                b=xk #Se redefine el b del intervalo

            else:
                a=xk #Se redefine el a del intervalo

            error=abs(sp.N(f1.subs('x',xk))) #Se cumple el teorema de Bolzano para el  subintervalo [xk,b]
            er.append(sp.N(error))

        #Condicion de parada cuando el error sea menor a la tolerancia establecida
            if error<tol: 
                break

    else:
        #Si el teorema de Bolzano no se cumple muestra el mensaje
        print("El intervalo seleccionado no cumple con el teorema de Bolzano") 

    plt.rcParams.update({'font.size': 14})
    ejex=np.arange(1,k+1,1)
    fig, graf=plt.subplots()
    #Se plotean la iteracion y el error
    graf.plot(ejex,er,'b',marker='o',markerfacecolor='red',markersize=10)
    graf.set_xlabel('Iteraciones ($k)')#Nombre del eje x
    graf.set_ylabel('$|f(x_k)|$') #Nombre eje y
    #Titulo grafica
    graf.set_title('Método de Newton-Raphson (Iteraciones vs Error)');  
    graf.grid(True) #Mostrar grid
    plt.show() #Mostrar grafica
    
    return [xk,k,error] #Luego de las iteraciones retorna la aproximacion y el error
        

  
#Definicion de parametros de entrada
    
f='exp(x)-2*x-10' 
a=2
b=4
tol=10**-5;
iterMax=1000;

y= biseccion(f,a,b,tol,iterMax)
print(y)


