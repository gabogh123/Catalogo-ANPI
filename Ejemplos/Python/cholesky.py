from math import *
import numpy as np

def cholesky(a, b):
    n = len(a)
    L = get_L(a, b, n)
    Lt =[[L[j][i] for j in range(n)] for i in range(n)]
    y = sustitucion_hacia_adelante(L, b, n)
    return sustitucion_hacia_atras(Lt, y, n)

def get_L(a, b, n):
    L = [[0 for i in range(n)] for j in range(n)]

    for i in range(n):
        for j in range(n):
            if i == j:
                sum1 = 0
                for k in range(j):
                    sum1 += L[j][k]**2
                L[i][j] = sqrt(a[i][j] - sum1)
            elif i > j:
                sum2 = 0
                for k in range(j):
                    sum2 += L[i][k]*L[j][k]
                L[i][j] = (a[i][j] - sum2) / L[j][j]
    return L

def sustitucion_hacia_atras(a, b, n):
    result = [None]*n
    i = n - 1
    while i >= 0:
        resultTemp = b[i]
        j = n - 1
        while j >= 0:
            if i == j:
                result[j] = resultTemp/a[i][j]
                break
            else:
                resultTemp -= a[i][j]*result[j]
            j -= 1
        i -= 1
    return result


def sustitucion_hacia_adelante(a, b, n):
    result = [None]*n
    for i in range (n):
        resultTemp = b[i]
        for j in range(n):
            if i == j:
                result[j] = resultTemp/a[i][j]
                break
            else:
                resultTemp -= a[i][j]*result[j]
    return result


print(cholesky([[25, 15, -5, -10], [15, 10, 1, -7], [-5, 1, 21, 4], [-10, -7, 4, 18]], [-25, -19, -21, -5]))

