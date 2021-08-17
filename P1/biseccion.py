def biseccion(f,a,b,tol,iterMax):

#Esta función aproxima la solución de la ecuación f(x)=0
#Parámetros de entrada: f = string que representa la función f,
#                       a,b = extremos del intervalo [a,b]
#                       tol = número que representa la tolerancia para la condición de parada
#                       iterMax = cantidad de iteraciones máximas
#Parámetros de salida:  xk = aproximación del cero o solución de la función f
#                       k = número de iteraciones realizadas
#                       error = grado de error del cálculo --> |f(xk)|


#Se importan las librerías necesarias para los cálculos
    import numpy as np
    import sympy as sp
    import matplotlib.pyplot as plt
    

#Se inicializan los parámetros y se convierte la expresión de string a función
    x=sp.Symbol('x')
    f1=sp.sympify(f)
    er=[]
    error=tol+1
    k=0

    if sp.N(f1.subs('x',a)) * sp.N(f1.subs('x',b))<0:

        for k in range(1,iterMax):

            xk=(a+b)/2 #Se crea el xk a la mitad del intervalo

            if sp.N(f1.subs('x',a))*sp.N(f1.subs('x',xk))<0: #Se cumple el teorema de Bolzano para el  subintervalo [a,xk]

                b=xk #Se redefine el b del intervalo

            else:
                a=xk #Se redefine el a del intervalo

            error=abs(sp.N(f1.subs('x',xk))) #Se cumple el teorema de Bolzano para el  subintervalo [xk,b]
            er.append(sp.N(error))

            if error<tol: #Condición de parada cuando el error sea menor a la tolerancia establecida
                break

    else:

        print("El intervalo seleccionado no cumple con el teorema de Bolzano") #Si el teorema de Bolzano no se cumple muestra el mensaje

    plt.rcParams.update({'font.size': 14})
    ejex=np.arange(1,k+1,1)
    fig, graf=plt.subplots()
    graf.plot(ejex,er,'b',marker='o',markerfacecolor='red',markersize=10)
    graf.set_xlabel('Iteraciones ($k)')
    graf.set_ylabel('$|f(x_k)|$')
    graf.set_title('Método de Newton-Raphson (Iteraciones vs Error)');
    graf.grid(True)
    plt.show()
    
    return [xk,k,error]
        

  
#Definición de parámetros de salida (pueden ser cambiados de acuerdo a lo que se desee calcular)
    
f='exp(x)-2*x-10' 
a=2
b=4
tol=10**-5;
iterMax=1000;

y= biseccion(f,a,b,tol,iterMax)
print(y)


