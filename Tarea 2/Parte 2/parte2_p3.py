from parte2_p2 import *

'''En el presente script se presenta la aproximación, por medio del método de Newton-Raphson, de la lista de 
funciones presentada en: A. Cordero, E. Martínez and J. R. Torregrosa, "Iterative methods of order four and five for 
systems of nonlinear equations," J. Comput. Appl. Math., vol. 231, (2), pp. 541-551, 2009. '''


def ejemplo_a():
    x = ['x1', 'x2']
    f = ['exp(x1^2) - exp(2^(0.5)*x1)', 'x1 - x2']
    x0 = [2.3, 2.3]  # solucion de artículo: x = (sqrt(2), sqrt(2))
    caller(f, x, x0)


def ejemplo_b():
    x = ['x1', 'x2']
    f = ['x1 + exp(x2) - cos(x2)', '3*x1 - x2 - sin(x2)']
    x0 = [1.5, 2]  # solucion de artículo: x = (0, 0)
    caller(f, x, x0)


def ejemplo_c():
    x = ['x1', 'x2']
    f = ['x1^2 - 2*x1 - x2 + 0.5', 'x1^2 + 4*x2^2 - 4']
    x0 = [3, 2]  # solucion de artículo: x = (1.9007, 0.3112)
    caller(f, x, x0)


def ejemplo_d():
    x = ['x1', 'x2']
    f = ['x1^2 + x2^2 - 1', 'x1^2 - x2^2 + 0.5']
    x0 = [0.7, 1.2]  # solucion de artículo: x = (0.5, sqrt(3)/2)
    caller(f, x, x0)


def ejemplo_e():
    x = ['x1', 'x2']
    f = ['sin(x1) + x2*cos(x1)', 'x1 - x2']
    x0 = [1.2, -1.5]  # solucion de artículo: x = (0, 0)
    caller(f, x, x0)


def ejemplo_g():
    x = ['x1', 'x2', 'x3', 'x4']
    f = ['x2*x3 + x4*(x2 + x3)',
         'x1*x3 + x4*(x1 + x3)',
         'x1*x2 + x4*(x1 + x2)',
         'x1*x2 + x1*x3 + x2*x3 - 1']
    x0 = [-1, -1, -1, -1]  # solucion de artículo: x = (-1/sqrt(3), -1/sqrt(3), -1/sqrt(3), 1/(2*sqrt(3)))
    caller(f, x, x0)


def caller(f, x, x0):
    """
    Método que hace la llamada de la función newton_raphson
    para los ejemplos del artículo.
    :param f:vector de funciones
    :param x: vector de variables
    :param x0: vector de valores iniciales
    """
    tol = 0.00000001
    iterMax = 50
    y = newton_raphson(x0, x, f, tol, iterMax)
    solucion = y[0]
    iteraciones = y[1]
    error = y[2]
    er = y[3]
    print(solucion)
    print(iteraciones)
    print(error)
    errors_plot(iteraciones, er)

ejemplo_a()
# ejemplo_b()
# ejemplo_c()
# ejemplo_d()
# ejemplo_e()
# ejemplo_g()
