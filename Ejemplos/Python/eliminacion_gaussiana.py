def eliminacion_gaussiana(a, b):
    n = len(a)
    (a, b) = triangular_superior(a, b, n)
    result = sustitucion_hacia_atras(a, b, n)
    print(result)


def triangular_superior(a, b, n):
    for j in range(n):
        for i in range(n):
            if i > j and a[i][j] != 0:
                multiplier = -1*a[i][j]/a[j][j]
                for k in range(n):
                    a[i][k] = a[i][k] + a[j][k] * multiplier
                b[i] = b[i] + b[j]*multiplier
    return a, b


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
                result[j] = [resultTemp/a[i][j]]
                break
            else:
                resultTemp -= a[i][j]*result[j]
    return result


eliminacion_gaussiana([[2, -6, 12, 16], [1, -2, 6, 6], [-1, 3, -3, -7], [0, 4, 3, -6]], [70, 26, -30, -26])