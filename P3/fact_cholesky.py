#Importacion de las librerias
import numpy as np
from math import *
import sys
eps = float(sys.float_info.epsilon)


def fact_cholesky(A, b):
    """
    Función que calcula la solución al sistema lineal Ax = b mediante el
    método de factorización de Cholesky
    Parámetros de entrada:
                A: Matriz de dimensiones nxn simétrica positiva definida
                b: vector xi elementos, donde i es la cantidad de incógnitas
    Salida: vector solución x del sistema lineal
       
    """
    n = len(A)
    L = [[0 for i in range(n)] for j in range(n)]
    if esSimetrica(A) and esPosDefinida(A):

        for i in range(n):
            for j in range(n):
                if i == j:
                    ljj = 0
                    for k in range(j):
                        ljj += L[j][k]**2
                    L[i][j] = sqrt(A[i][j] - ljj)
                elif i > j:
                    lij = 0
                    for k in range(j):
                        lij += L[i][k]*L[j][k]
                    L[i][j] = (A[i][j] - lij) / L[j][j]

                    
        Lt =[[L[j][i] for j in range(n)] for i in range(n)]
        y = sust_adelan(L, b, n)
        return sust_atras(Lt, y, n)
    else:
        print("La matriz no es simétrica positiva definida")
    

def sust_atras(A, b, n):
    """
    Función que realiza las operaciones para resolver un
    el sistema triangular superior mediante el método
    de sustitución hacia atrás.
    Parametros de entrada:
                A: matriz triangular superior de dimensiones nxn
                b: vector modificado según las operaciones de fila
                n: dimensión de la matriz
    Salida: vector solución del sistema
    """
    result = [None]*n
    i = n - 1
    while i >= 0:
        aux = b[i]
        j = n - 1
        while j >= 0:
            if i == j:
                result[j] = aux/A[i][j]
                break
            else:
                aux -= A[i][j]*result[j]
            j -= 1
        i -= 1
    return result

def sust_adelan(A, b, n):
    """
    Función que realiza las operaciones para resolver un
    el sistema triangular inferior mediante el método
    de sustitución hacia atrás.
    Parámetros de entrada:
                A: matriz triangular inferior de dimensiones nxn
                b: vector modificado según las operaciones de fila
                n: dimensión de la matriz
    Salida: vector solución
    """
    result = [None]*n
    for i in range (n):
        resultTemp = b[i]
        for j in range(n):
            if i == j:
                result[j] = resultTemp/A[i][j]
                break
            else:
                resultTemp -= A[i][j]*result[j]
    return result

def esSimetrica(A):
    """
    Función que valida si la matriz en cuestión es
    simétrica.
    Parámetros de entrada:
                A: matriz nxn
    Salida: valor booleano
    """
    n = len(A)
    At = [[A[j][i] for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] != At[i][j]:
                return False
    return True

def esPosDefinida(A):
    """
    Función que valida si la matriz en cuestión es
    positiva definida.
    Parámetros de entrada:
                A: matriz nxn
    Salida: valor booleano
    """
    A = np.array(A)
    n = A.shape[0]
    j = 0
    for i in range(n+1): 
        step = A[:i,:j]
        j += 1
        det = np.linalg.det(step)
        if det < 0 : 
            return False

    return True    
 
if __name__ == "__main__":
       
    #Ejemplo presentación 4 diap 69
    A , b = [[25, 15, -5, -10], [15, 10, 1, -7], [-5, 1, 21, 4], [-10, -7, 4, 18]], [-25, -19, -21, -5]
    print("El sistema de lineal a resolver Ax = b, donde A:\n ", A)
    print("Con vector b: ", b)
    print("Tiene como solución x:")
    print(fact_cholesky(A, b))
