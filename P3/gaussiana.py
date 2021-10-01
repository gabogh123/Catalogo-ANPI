#Importacion de las librerias
import numpy as np
import math
import sys
eps = float(sys.float_info.epsilon)

def gaussiana(A, b):
    print("Sistema lineal a resolver: \n ", A)
    print("Vector b: ", b)
    
    _A, _b = triang_sup(A, b)
    X = sust_atras(_A, _b)
    print("El vector solución al sistema es: ", X)

def triang_sup(A, b):
    """
    Función que realiza las operaciones de fila para convertir
    al sistema de matriz densa a uno de tipo trangular superior.
    Parámetros de entrada:
                A: matriz de dimensiones nxn invertible
                b: vector xi elementos, donde i es la cantidad de incógnitas
    Salida: matriz triangular superior _A, vector de matriz aumentada _b        
    
    """
    n = len(b)
      
    for k in range(0, n-1):
        if A[k, k] <= eps:
            print("Error division por 0")
            return None
        for i in range(k + 1, n):
            pivot = A[i, k]/A[k, k]
            for j in range(k, n):
                A[i, j] -= pivot*A[k, j]
            b[i] -= pivot*b[k]
    _A = A
    _b = b
    return _A, _b
            

def sust_atras(_A, _b):
    """
    Función que realiza las operaciones para resolver el
    el sistema triangular superior mediante el método
    de sustitución hacia atrás.
    Parametros de entrada:
                _A: matriz triangular superior de dimensiones nxn
                _b: vector modificado según las operaciones de fila
    Salida: vector solución al sistema
    """
    m = _b.size
    x = np.zeros_like(_b)
    if _A[m - 1, m - 1] <= eps:
        error = "Error: denominador no puede ser cero"
        return (error)
    x[m-1] = _b[m-1]/_A[m-1, m-1]
    colAux = np.zeros((m, m))
    
    for i in range(m - 2, -1, -1):
        aux = 0
        for j in range(i + 1, m):
            aux += _A[i, j]*x[j]
        colAux[i, i] = _b[i] - aux
        x[i] = colAux[i, i]/_A[i, i]
    
    return x

    
if __name__ == "__main__":
    """
    Ejemplo presentación 4 diap 23
    """
    A = np.array([[2,-6,12,16],
               [1, -2, 6,6],
               [-1,3,-3,-7],
               [0,4,3,-6]])
    b1 = np.array([70, 26,-30,-26])
    b = b1.T
    gaussiana(A, b)

