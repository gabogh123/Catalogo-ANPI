def factorizacion_LU(a, b):
    n = len(a)
    L, U = get_LU(a, n)
    y = sustitucion_hacia_adelante(L, b, n)
    return sustitucion_hacia_atras(U, y, n)




def get_LU(a, n):
    L = [[0 for i in range(n)] for j in range(n)]
    for j in range(n):
        for i in range(n):
            if i > j and a[i][j] != 0:
                multiplier = -1 * a[i][j] / a[j][j]
                L[i][j] = -1 * multiplier
                for k in range(n):
                    a[i][k] = a[i][k] + a[j][k] * multiplier
            elif i == j:
                L[i][j] = 1
    return L, a

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

print(factorizacion_LU([[4, -2, 1], [20, -7, 12], [-8, 13, 17]], [11, 70, 17]))
