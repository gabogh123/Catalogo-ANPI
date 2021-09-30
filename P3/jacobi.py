import numpy as np

def jacobi(A,b,x0,tol,iterMax):

    k = 0
    
    n = len(A)




    L = np.zeros((n,n))
    D = np.zeros((n,n))
    U = np.zeros((n,n))

    for i in range(n):
        for j in range(n):

            if i>j:
                L[i][j] = A[i][j]
            elif i==j:
                D[i][j] = A[i][j]
            else:
                U[i][j] = A[i][j]

    xk = np.array(x0)
    error = np.linalg.norm(np.dot(A, xk) - b)
    

    D_inversa = np.linalg.inv(D)

    m_punto_n = np.dot(D_inversa, (-1 * (L + U)))

    D_inversa_punto_b = np.dot(D_inversa, b)

    if np.linalg.norm(m_punto_n) < 1:

        while k < iterMax:

            xk = np.dot(m_punto_n, xk) + D_inversa_punto_b
            error = np.linalg.norm(np.dot(A, xk) - b)

            k += 1

            if error < tol:
                break
    else:
        print("La matriz no cumple la condición de que la norma infinita de -D_inv*(L+U) sea menor a cero")
   
    
    print(xk , error)


A = [[5,1,1],[1,5,1],[1,1,5]]

b = [[7],[7],[7]]

x0 = [[0],[0],[0]]

tol = 0.0000001

iterMax = 1000

y = jacobi(A,b,x0,tol,iterMax)



 
    
    
    
    
