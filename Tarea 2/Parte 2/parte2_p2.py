# Se importan las bibliotecas necesarias para los cálculos
# Python v 3.6
import numpy as np
import sympy as sp
import math
import matplotlib.pyplot as plt


def newton_raphson(x_0, x, f, tol, iterMax):
    """
    Funcion que aplica el método iterativo de Newton-Raphson para
    sistemas de ecuaciones no lineales como estrategia de solución al sistema
    con una aproximación dada de acuerdo a la tolerancia ingresada
    :param x_0: vector de valores iniciales
    :param x: vector string de variables del sistema
    :param f: vector string de funciones del sistema de la forma F(x_i) = 0
    :param tol: tolerancia aceptada del sistema
    :param iterMax: cantidad máxima de iteraciones
    :return: x_k como solución del sistema, k como cantidad de iteraciones,
            error como error de la última iteración donde se alcanza la tolerancia aceptada,
            er como vector de errores registrados para hacer su gráfica
    """
    k = 0
    f_ = simplify(f)
    my_vars = symbols(x)
    x_k = np.array(x_0)
    jacobian = get_jacobian(x, f_)
    er = []
    while k < iterMax:
        f_xk = eval_f(x_k, my_vars, f_)
        j_k = eval_jacobian(my_vars, x_k, jacobian)
        y = np.linalg.solve(j_k, f_xk)
        x_k = x_k - y
        error = np.linalg.norm(eval_f(x_k, my_vars, f_), ord=2)
        er.append(error)
        if error < tol:
            break
        k += 1

    return [x_k, k, error, er]


def simplify(f):
    """
    Funcion para simplificar la expresion de cada
    funcion del sistema
    :param f: sistema de ecuaciones
    :return: vector de funciones simplificadas
    """
    system = len(f)
    length = np.arange(0, system, 1)
    simp_f = []
    for i in length:
        simp_f.append(sp.simplify(f[i]))
    return simp_f


def symbols(vars):
    """
    Funcion que convierte un vector string de variables a
    vector simbolico de variables
    :param my_vars: vector string
    :return: vector simbolico
    """
    n = len(vars)
    length = np.arange(0, n, 1)
    my_vars = []
    for i in length:
        my_vars.append(sp.Symbol(vars[i]))
    return my_vars

def eval_f(x_k, my_vars, f_):
    """
    Funcion que se encarga de evaluar por sustitucion numerica
    la cada funcion del sistema de ecuaciones no lineales f_
    :param my_vars: variables de las que depende el sistema
    :param x_k: vector actual en el que se va a evaluar la función
    :param f_: vector de funciones simbólicas
    :return: b: vector con los valores obtenidos de evaluar x_k en el sistema
    """
    n = len(f_)
    length = np.arange(0, n, 1)
    f_ns = []
    b_k = []
    for i in length:
        f_n = sp.lambdify(my_vars, f_[i])
        f_ns.append(f_n)
    for i in length:
        if n == 2:
            current_f = f_ns[i]
            current_eval = current_f(x_k[0], x_k[1])
            b_k.append(current_eval)
        elif n == 3:
            current_f = f_ns[i]
            current_eval = current_f(x_k[0], x_k[1], x_k[2])
            b_k.append(current_eval)
        elif n == 4:
            current_f = f_ns[i]
            current_eval = current_f(x_k[0], x_k[1], x_k[2], x_k[3])
            b_k.append(current_eval)
    b = np.array(b_k)
    return b


def get_jacobian(my_vars, f_):
    """
    Función que calcula la matriz jacobiana del sistema
    :param my_vars: variables de las que depende el sistema
    :param f_: sistema de ecuaciones en forma simbolica
    :return: matriz jacobiana
    """
    n = len(f_)
    length = np.arange(0, n, 1)
    jacobian = []
    for i in length:
        current_jaco = derivative(f_[i], my_vars)
        jacobian.append(current_jaco)
    jacobian_i = np.array(jacobian)
    return jacobian_i


def eval_jacobian(my_vars, x_k, g_):
    """
    Funcion que se encarga de evaluar por sustitucion numerica
    la matriz jacobiana del sistema de ecuaciones
    no lineales f_
    :param my_vars: variables de las que depende el sistema
    :param x_k: vector actual en el que se va a evaluar la matriz
    :param g_: matriz jacobiana simbolica
    :return: matriz de numpy con los valores evaluados para la matriz jacobiana en x_k
    """
    n = len(g_)
    length = np.arange(0, n, 1)
    g_ns = []
    j_k = []
    for i in length:
        g_n = sp.lambdify(my_vars, g_[i]) #lambdify para habilitar la sustitucion numerica de cada función derivada parcialmente
        g_ns.append(g_n)
    for i in length:
        if n == 2:  #si el sistema depende de dos variables
            current_g = g_ns[i]
            current_eval = current_g(x_k[0], x_k[1])
            j_k.append(current_eval)
        elif n == 3: #si el sistema depende de tres variables
            current_g = g_ns[i]
            current_eval = current_g(x_k[0], x_k[1], x_k[2])
            j_k.append(current_eval)
        elif n == 4: #si el sistema depende de cuatro variables
            current_g = g_ns[i]
            current_eval = current_g(x_k[0], x_k[1], x_k[2], x_k[3])
            j_k.append(current_eval)
    j = np.array(j_k)
    return j


def derivative(f, my_vars):
    """
    Funcion que se encarga de calcular la primera derivada parcial de una función dada
    según sus variables
    :param f: funcion simbolica a derivar
    :param my_vars: variables simbolicas de las que depende f
    :return: vector con la derivada parcial en cada variable
    """
    n = len(my_vars)
    length = np.arange(0, n, 1)
    df_i = []
    for i in length:
        df_i.append(sp.diff(f, my_vars[i]))
    return df_i


def errors_plot(k, er):
    """
    Funcion que se encarga de graficar el comportamiento de
    iteraciones vs error
    :param k: cantidad de iteraciones
    :param er: vector con el registro de errores por cada iteracion
    """
    plt.rcParams.update({'font.size': 10})
    ejex=np.arange(0,k+1,1)
    fig, graf=plt.subplots()
    #Se plotean la iteracion y el error
    graf.plot(ejex,er,'b',marker='o',markerfacecolor='red',markersize=10)
    graf.set_xlabel('Iteraciones ($k$)')#Nombre del eje x
    graf.set_ylabel('$||f(x^{(k)})||_{2}$') #Nombre eje y
    #Titulo grafica
    graf.set_title('N-R Sistemas de Ecuaciones No Lineales (Iteraciones vs Error)');
    graf.grid(True) #Mostrar grid
    plt.show() #Mostrar grafica