#Se importan la libreria numpy
import numpy as np

def jacobi(A,b,x0,tol,iterMax):

#Esta funcion aproxima la solcuión de un sistema de ecuacuiones
#con matrices mediante el método de Jacobi
#Parametros de entrada: A: matriz de coeficientes
#                       b: vector de términos independientes
#                       x0: valor inicial
#                       tol: tolerancia del metodo
#                       iterMax: iteraciones máximas del método
#Parametros de salida: xk: matriz solucion aproximada
#                      error: error del método de jacobi

#Se inicializa el k y el error  
    k = 0
    error = tol+1

#Se define el n como el tamaño de la matriz
    n = len(A)

#Se crean las matrices L, D y U rellenadas con ceros
    L = np.zeros((n,n))
    D = np.zeros((n,n))
    U = np.zeros((n,n))

    
#Se rellenan las matrices L,D,U, de acuerdo a las posiciones de la matriz A
    for i in range(n):
        for j in range(n):

            if i>j:
                L[i][j] = A[i][j]
            elif i==j:
                D[i][j] = A[i][j]
            else:
                U[i][j] = A[i][j]


#Se hace el xk igual a l x0 para la primera iteracion
    xk = np.array(x0)
    
#Se realiza la inversa de D
    D_inversa = np.linalg.inv(D)
    
#Se define la expresion D inversa * -(L+U)
    
    m_punto_n = np.dot(D_inversa, (-1 * (L + U)))


#Se define la expresion D inversa * b
    D_inversa_punto_b = np.dot(D_inversa, b)
    
#Se valida que la norma de M*N sea menor a 1, condición para realizar jacobi
    if np.linalg.norm(m_punto_n) < 1: 

#Se realiza mientras el k de iteraciones sean menor a la iteraciones máximas
        while k < iterMax: 

#Se redefine el xk mediante la función punto de M*N y xk sumado con D_inv*b
            xk = np.dot(m_punto_n, xk) + D_inversa_punto_b

#Se calcula el error como la norma de A*x - b            
            error = np.linalg.norm(np.dot(A, xk) - b)

            k += 1 #Se avanza en la iteracióm

            if error < tol: #Segunda condición de parada, error menor a la tol
                break
#Si la matriz no cumple la condición para ser resuleta por jacpbi retorna el mensaje
    else:
        print("La matriz no cumple la condición de que la norma infinita de -D_inv*(L+U) sea menor a cero")
   
#Se retorna la matriz final y el error obtenido    
    print( 'Aproximación: \n',xk, '\n Error:' , error)


A = [[5,1,1],[1,5,1],[1,1,5]]
b = [[7],[7],[7]]
x0 = [[0],[0],[0]]
tol = 0.0000001
iterMax = 1000

y = jacobi(A,b,x0,tol,iterMax)



 
    
    
    
    
