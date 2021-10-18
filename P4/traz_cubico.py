from sympy import *
import numpy as np

def traz_cubico(x_n, y_n):
    """
        Métodos de interpolación polinomial
        -Trazador Cúbico Natural
        Parámetros de entrada:
            x_n: vector de preimágenes
            y_n: vector de imágenes
        Parámetros de salida:
            S_i: vector simbólico de trazadores cúbicos
    """
    x = Symbol('x')
    h_i = calc_hi(x_n)          #Paso 1: calcular h_i
    u_i = calc_ui(y_n, h_i)     #Paso 2: calcular u_i
    A = crear_A(h_i)            #Paso 3: crear matriz tridiagonal A
    m = thomas(A, u_i)          #Paso 4: resolver el sistema Am = u_i
    n = len(x_n)        
    m_i = np.zeros(n)           #Paso 5: inicializar m_0 y m_n en 0

    for i in range(1, n-1):     #Ubica los valores obtenidos en la solución del sistema Am = u_i 
        m_i[i] = m[i-1]         #en la matriz m_i
    a_i = np.zeros(n-1)
    b_i = np.zeros(n-1)
    c_i = np.zeros(n-1)
    d_i = np.zeros(n-1)

    #-----------------Paso 6: calcular coeficientes--------------------#
    for i in range(n-1): 
        a_i[i] = (m_i[i+1]-m_i[i])/(6*h_i[i])
        b_i[i] = m_i[i]/2
        d_i[i] = y_n[i]
        c_l    = (y_n[i+1]-y_n[i])/h_i[i]
        c_r    = h_i[i]*(m_i[i+1] + 2*m_i[i])/6
        c_i[i] = c_l - c_r

    S_i = []
    #-----------------Paso 7: calcular los trazadores y simplificarlos--------------------#
    for i in range(n-1):
        s = a_i[i]*(x - x_n[i])**3  + b_i[i]*(x - x_n[i])**2 + c_i[i]*(x - x_n[i]) + d_i[i]
        s = simplify(s)
        S_i.append(s)

    return S_i

def calc_hi(x_n):
    """
        Función para calcular los diferenciales entre las preimágenes
        input:
            x_n: vector de preimágenes
        output:
            h_i: vector de diferenciales
    """
    h_i = []
    n = len(x_n)
    for i in range(n-1):
        h = x_n[i+1] - x_n[i]
        h_i.append(h)
    return h_i

def calc_ui(y_n, h_i):
    """
        Función para calcular el vector u_i
        input:
            y_n: vector de imágenes
            h_i: vector de diferenciales
        output:
            u_i: vector del sistema Am = u_i
    """
    u_i = []
    n = len(y_n)
    for j in range(n-2):
        u_l = (y_n[j+2]- y_n[j+1])/h_i[j+1]
        u_r = (y_n[j+1]- y_n[j])/h_i[j]
        u = 6*(u_l - u_r)
        u_i.append(u)
    return u_i

def crear_A(h_i):
    """
        Función para inicializar la matriz A
        input:
            h_i: vector de diferenciales
        output:
            A: matriz tridiagonal
    """
    n = len(h_i) - 1
    A = np.zeros((n, n))
    for i in range(1, n):
        A[i-1][i] = h_i[i]
        A[i][i-1] = h_i[i]
    for i in range(n):
        A[i][i] = 2*(h_i[i] + h_i[i+1])
    return A

#-------------------------Método de Thomas-----------------------------#
def thomas(A, u_i):
    n = len(u_i)
    a, b, c = vectores(A, n)
    p_i, q_i = get_pq(a, b, c, u_i)
    x_i = np.zeros(n)
    x_i[n - 1] = q_i[n-1]
    for i in np.arange(n-2, -1, -1):
        x_i[i] = q_i[i] - p_i[i]*x_i[i+1]

    return x_i
    
def vectores(A, n):
    a = [0] 
    b = []
    c = []
    for i in range(1, n):
        c_i = A[i-1][i]
        a_i = A[i][i-1]
        c.append(c_i)
        a.append(a_i)
        if i == n-1:
            c.append(0)
            
    for i in range(n):
        b_i = A[i][i]
        b.append(b_i)
    return a, b, c

def get_pq(a, b, c, u_i):
    n = len(u_i)
    p_i = []
    q_i = []
    for i in range(n-1):
        if i == 0:
            p = c[i]/b[i]
            p_i.append(p)
        else:
            p = c[i]/(b[i]-p_i[i-1]*a[i])
            p_i.append(p)
    for i in range(n):
        if i == 0:
            q = u_i[i]/b[i]
            q_i.append(q)
        else:
            q = (u_i[i] - q_i[i-1]*a[i])/(b[i]-p_i[i-1]*a[i])
            q_i.append(q)
            
    return p_i, q_i

#--------------------Ejemplo presentación 9, diap 25-----------------------------#
"""
Donde cada par (x, y) pertenece a la función f(x) = 3*x*exp(x) - 2*exp(x)
"""
x_n = [1, 1.05, 1.07, 1.1]
y_n = [2.718282, 3.286299, 3.527609, 3.905416]
y = traz_cubico(x_n, y_n)
for i in range(len(y)):
    print(y[i])
