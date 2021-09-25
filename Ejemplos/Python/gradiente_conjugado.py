from sympy import *
import numpy as np
from numpy import linalg
import random


def gradiente_conjugado(x0, variables, tol, funcion, regla_bk='FR'):
    variables = symbols(variables)
    funcion = sympify(funcion)
    f = lambdify(variables, funcion, "numpy")
    derivatives = []
    for variable in variables:
        derivatives.append(lambdify(variables, diff(funcion, variable)))
    x_k = np.array([x0])
    g_k_1 = g_k = gradiente(derivatives, x_k)
    d_k = -1 * g_k
    iteracion = 0
    error = tol
    while error >= tol:
        alpha = get_alpha(f, x_k, d_k, g_k_1, 1)
        x_k = x_k + alpha * d_k
        g_k = g_k_1
        g_k_1 = gradiente(derivatives, x_k)
        d_k = -1 * g_k_1 + beta(regla_bk, g_k_1, g_k, d_k) * d_k
        error = linalg.norm(g_k_1)
        iteracion += 1
    return x_k[0].tolist(), iteracion


def gradiente(derivatives, values):
    result = []
    for derivative in derivatives:
        result.append(derivative(*values[0]))
    return np.array([result])


def get_alpha(f, x_k, d_k, g_k, alpha):
    delta = random.uniform(0, 1) #0.5 #random.randint(1, 9) * 0.1
    f_arg = x_k + alpha * d_k
    while f(*f_arg[0]) - f(*x_k[0]) > delta * alpha * g_k.dot(d_k.T):
        alpha /= 2
        f_arg = x_k + alpha * d_k
    return alpha


def beta(regla, g_k_1, g_k, d_k):
    if regla.lower() == 'cd':
        b = linalg.norm(g_k_1)**2 / ((-1 * d_k).dot(g_k.T))
    elif regla.lower() == 'dy':
        b = linalg.norm(g_k_1)**2 / d_k.dot((g_k_1 - g_k).T)
    else:
        b = linalg.norm(g_k_1)**2 / linalg.norm(g_k)**2
    return b


print(gradiente_conjugado([0, 3], ['x', 'y'], 0.001, '(x-2)**4+(x-2*y)**2', 'fr'))
print(gradiente_conjugado([0, 3], ['x', 'y'], 0.001, '(x-2)**4+(x-2*y)**2', 'cd'))
print(gradiente_conjugado([0, 3], ['x', 'y'], 0.001, '(x-2)**4+(x-2*y)**2', 'dy'))

print(gradiente_conjugado([1, 1, 1], ['x', 'y', 'z'], 0.001, '(x-2)**2+ (y+3)**2 + (x + y + z)**2', 'fr'))
print(gradiente_conjugado([1, 1, 1], ['x', 'y', 'z'], 0.001, '(x-2)**2+ (y+3)**2 + (x + y + z)**2', 'cd'))
print(gradiente_conjugado([1, 1, 1], ['x', 'y', 'z'], 0.001, '(x-2)**2+ (y+3)**2 + (x + y + z)**2', 'dy'))
