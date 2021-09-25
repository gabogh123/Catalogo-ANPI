import numpy as np


"""
Metodo iterativo de gauss hacia adelante para solucion de sistemas de 
ecuaciones.
Entradas: Matriz de Coeficientes, matriz de terminos independientes,
valor inicial y tolerancia.
Salidas: aproximacion y cantidad de iteraciones.
"""
def gauss_adelante(A, b, x0, tol):
    n = len(A)
    A = np.array(A)
    b = np.array(b)
    L = np.zeros((n, n))
    D = np.zeros((n, n))
    U = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if i > j:
                L[i][j] = A[i][j]
            elif i < j:
                U[i][j] = A[i][j]
            else:
                D[i][j] = A[i][j]

    x = np.array(x0)
    error = np.linalg.norm(np.dot(A, x) - b)
    i = 0

    m_inv = np.linalg.inv(D + U)
    m_times_n = np.dot(m_inv, (-1 * L))
    m_times_b = np.dot(m_inv, b)

    while error >= tol:
        x = np.dot(m_times_n, x) + m_times_b
        error = np.linalg.norm(np.dot(A, x) - b)
        i += 1

    return x, i


print(gauss_adelante([[3,1,1],[1,3,1],[1,1,3]], [[5],[5],[5]], [[-1],[0],[-1]], 0.001))