//
// Created by juanpr on 23/8/21.
#include <boost/multiprecision/cpp_int.hpp>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <boost/multiprecision/cpp_bin_float.hpp>
//
using namespace boost::multiprecision;
using namespace std;

cpp_bin_float_100 factor(cpp_bin_float_100 x){
    /**
     * @Brief Función para calcular el factorial del dato ingresado
     * @param x es el valor al que se desea calcular su factorial
     * @return x!
     */
    cpp_bin_float_100 num = 1.0;
    cpp_bin_float_100 i = 1.0;
    while (i <= x){
        num = i*num;
        i = i + 1.0;
    }
    return num;
}

cpp_bin_float_100 div_tt(cpp_bin_float_100 x){
    /**
     * @Brief Función para calcular el inverso multiplicativo del dato ingresado
     * @param x es el valor al que se desea calcular su inverso multiplicativo
     * @return 1 / x
     */
     /// Definicion de parametros condicionales
    int k = 0;
    int iterMax = 5000;
    cpp_bin_float_100 tol = 0.00000001;
    cpp_bin_float_100 error = tol+1;

    cpp_bin_float_100 flag = 0.0; ///Comparador del condicional para validar que la entrada no sea incorrecta
    if ( flag < x) {
        cpp_bin_float_100 x0 = 0.0;
        ///Variable que inicializa el valor de la precision relativa de punto flotante
        cpp_bin_float_100 eps = 2.220446049250313 * pow(10, -16);
        /// condicionales para validar el valor inicial segun la magnitud de x
        if (factor(0.00) < x && x <= factor(20.00)) {
            x0 = pow(eps, 2);
        } else if (factor(20.00) < x && x <= factor(40.00)) {
            x0 = pow(eps, 4);
        } else if (factor(40.00) < x && x <= factor(60.00)) {
            x0 = pow(eps, 8);
        } else if (factor(60.00) < x && x <= factor(80.00)) {
            x0 = pow(eps, 11);
        } else {
            x0 = pow(eps, 15);
        }
        cpp_bin_float_100 x1  = 0.0;
        while (k <= iterMax) {
            ///Aproximaciòn de la soluciòn por cada iteraciòn
            x1  = x0 * (2 -  (cpp_bin_float_100)x * x0) ;
            ///Error asociado al calculo de la aproximaciòn
            error = abs((x1 - x0) / x1);
            if (error < tol) break;
            x0 = x1;
            k++;
        }
        return x1;
    } else printf("el valor de entrada debe ser entero mayor que 0");
    return flag;
}

cpp_bin_float_100 sin_t(cpp_bin_float_100 x){
    /**
     * @Brief Función para calcular el valor de la funcion matematica sen(x) del dato ingresado
     * de acuerdo con la aproximacion por series de Taylor
     * @param x es el argumento de la funcion
     * @return sen(X)
     */
    /// Definicion de parametros condicionales
    int k = 0;
    int iterMax = 5000;
    cpp_bin_float_100 tol = 0.00000001;
    cpp_bin_float_100 error = tol+1;

    cpp_bin_float_100 sk = 0.0;
    cpp_bin_float_100 sk1 = 0.0;
    cpp_bin_float_100 n1 = 0.0;
    cpp_bin_float_100 n2 = 0.0;
    while (k <= iterMax){
        n1 = pow(-1, k);
        n2 = pow(x, 2*k +1);
        sk1 = n1*n2*div_tt(factor(2.0*k+1)) + sk1;
        error = abs(sk1-sk);
        if (error < tol) break;
        sk = sk1;
        k++;
    }
    return sk1;
}

cpp_bin_float_100 cos_t(cpp_bin_float_100 x){
    /**
     * @Brief Función para calcular el valor de la funcion matematica cos(x) del dato ingresado
     * de acuerdo con la aproximacion por series de Taylor
     * @param x es el argumento de la funcion
     * @return cos(X)
     */
    /// Definicion de parametros condicionales
    int k = 0;
    int iterMax = 5000;
    cpp_bin_float_100 tol = 0.00000001;
    cpp_bin_float_100 error = tol+1;

    cpp_bin_float_100 sk = 0.0;
    cpp_bin_float_100 sk1 = 0.0;
    cpp_bin_float_100 n1 = 0.0;
    cpp_bin_float_100 n2 = 0.0;
    while (k <= iterMax){
        n1 = pow(-1, k);
        n2 = pow(x, 2*k);
        sk1 = n1*n2*div_tt(factor(2.0*k)) + sk1;
        error = abs(sk1-sk);
        if (error < tol) break;
        sk = sk1;
        k ++;
    }
    return sk1;
}

cpp_bin_float_100 tan_t(cpp_bin_float_100 x){
/**
     * @Brief Función para calcular el valor de la funcion matematica tan(x) del dato ingresado
     * de acuerdo con la propiedad matematica tan(x) = sen(x)/cos(x)
     * @param x es el argumento de la funcion
     * @return tan(X)
     */
    cpp_bin_float_100 a = sin_t(x);
    cpp_bin_float_100 b = div_tt(cos_t(x));
    cpp_bin_float_100 c = a*b;
    return c;
}

cpp_bin_float_100 ln_t(cpp_bin_float_100 x){
    /**
     * @Brief Función para calcular el valor de la funcion matematica Ln(x) del dato ingresado
     * de acuerdo con la aproximacion por series de Taylor
     * @param x es el argumento de la funcion
     * @return ln(x)
     */
    /// Definicion de parametros condicionales
    int k = 0;
    int iterMax = 5000;
    cpp_bin_float_100 tol = 0.00000001;
    cpp_bin_float_100 error = tol+1;

    cpp_bin_float_100 p1 = (x - 1.0);
    cpp_bin_float_100 a = 2.0*p1*div_tt(x + 1.0);
    cpp_bin_float_100 sk = 0.0;
    cpp_bin_float_100 sk1 = sk;
    if (x != 0.0) {
        while (k <= iterMax) {
            cpp_bin_float_100 expo= pow(p1 * div_tt(x + 1.0), 2 * k);
            sk1 = a * expo * div_tt(2.0 * k + 1.0) + sk1;
            error = abs(sk1 - sk);
            if (error < tol) break;
            sk = sk1;
            k++;
        }
        return sk1;
    }
    else {
        cout << "Error. Ingrese un numero diferente de 0."<<endl;
        return 0.0;
    }
}

cpp_bin_float_100 exp_t(cpp_bin_float_100 x){
    /**
     * @Brief Función para calcular el valor de la funcion matematica e⋀x del dato ingresado
     * de acuerdo con la aproximacion por series de Taylor
     * @param x es el argumento de la funcion
     * @return e⋀x
     */
    /// Definicion de parametros condicionales
    int k = 0;
    int iterMax = 5000;
    cpp_bin_float_100 tol = 0.00000001;
    cpp_bin_float_100 error = tol+1;

    cpp_bin_float_100 sk = 0.0;
    cpp_bin_float_100 sk1 = 0.0;
    cpp_bin_float_100 n1 = 0.0;
    if (x != 0.0) {
        while (k <= iterMax) {
            n1 = (cpp_bin_float_100)pow(x, k);
            sk1 = n1 * div_tt(factor(1.0*k)) + sk1;
            error = abs(sk1 - sk);
            if (error < tol) break;
            sk = sk1;
            k++;
        }
        return sk1;
    }
    else return 1.0;
}

cpp_bin_float_100 f(cpp_bin_float_100 x, int p, cpp_bin_float_100 a){
    /**
     * @Brief Función auxiliar para evaluar valores en la funcion matematica
     * f(x) = x⋀p - a
     * @param x es la variable a aproximar
     * @param p es índice de la raíz original
     * @param a es el subradical de la función original
     * @return f(x)
     */
    return pow(x, p) - a;
}

cpp_bin_float_100 dfx(cpp_bin_float_100 x, int p){
    /**
     * @Brief Función derivada de f(x)
     * @param x es la variable independiente
     * @param p potencia de la función f(x)
     * @return f'(x)
     */
    return p*pow(x, p-1);
}

cpp_bin_float_100 root_t(cpp_bin_float_100 x, int a){
    /**
     * @Brief Función para calcular la raíz enésima de cualquier valor real a√x
     * mediante el cálculo de la raíz f(x) = 0 de f(x) = x⋀p - a con el método de Newton-Raphson
     * utizando las funciones auxiliares f(x, p, a) y dfx(x , p)
     * @param x es el subíndice del radical
     * @param a es el ìndice
     * @return f(x)
     */
    int k = 0;
    int iterMax = 5000;
    cpp_bin_float_100 tol = 0.00000001;
    cpp_bin_float_100 error = tol+1;

    cpp_bin_float_100 x0 = x*div_tt(2.0);
    cpp_bin_float_100 h = f(x0, a, x)/dfx(x0, a);
    cpp_bin_float_100 x1 = 0.0;
    while (k <= iterMax){
        h = f(x0, a, x)/dfx(x0, a);
        x0 = x0 - h;
        error = abs((x1-x0)/x1);
        if (error < tol) break;
        x1 = x0;
        k++;
    }
    return x1;
}
cpp_bin_float_100 power_t(cpp_bin_float_100 a, cpp_bin_float_100 x){
    /**
     * @Brief Funcion para calcular el valor de funciones de la forma f(x) = a^x
     * utilizando las propiedades de logaritmos
     * 1. e^ln[f(x)] = f(x)
     * 2. ln(a^x) = x * ln(a)
     */
    cpp_bin_float_100 logNatural  = x*ln_t(a);
    cpp_bin_float_100 result = exp_t(logNatural);
    return result;
}
cpp_bin_float_100 pi(){
    /**
     * @Brief Función para  aproximar la constante matematica Pi
     * Mediante la serie de Leibniz en https://en.wikipedia.org/wiki/Leibniz_formula_for_%CF%80
     * @return pi
     */
    int k = 0;
    int iterMax = 5000;
    cpp_bin_float_100 tol = 0.00000001;
    cpp_bin_float_100 error = tol+1;

    cpp_bin_float_100 sk1 = 0.0;
    cpp_bin_float_100 num = 0.0;
    cpp_bin_float_100 den = 0.0;
    while (k <= iterMax){
        num = 4.0*pow(-1, k);
        den = div_tt(2*k + 1);
        sk1 = num*den + sk1;
        k++;
    }
    return sk1;
}

cpp_bin_float_100 log_t(cpp_bin_float_100 x, cpp_bin_float_100 a){
    /**
     * @Brief Función para  aproximar el valor de un logarítmo utilizando
     * la propiedad de cambio de base log_a(x) = ln(x)/ln(a)
     * @param x es el argumento de la función
     * @param x es la base del la función
     * @return log_a(x)
     */
    int k = 0;
    int iterMax = 5000;
    cpp_bin_float_100 tol = 0.00000001;
    cpp_bin_float_100 error = tol+1;

    cpp_bin_float_100 arg = ln_t(x);
    cpp_bin_float_100 base = div_tt(ln_t(a));
    cpp_bin_float_100 result = arg*base;
    return result;
}

cpp_bin_float_100 asin_t(cpp_bin_float_100 x){
    /**
     * @Brief Función para calcular el valor de la funcion matematica arcsin(x) del dato ingresado
     * de acuerdo con la aproximacion por series de Taylor que converge cuando |x| <= 1
     * @param x es el argumento de la funcion
     * @return sin^⁻1 (x)
     */
    /// Definicion de parametros condicionales
    int k = 0;
    int iterMax = 5000;
    cpp_bin_float_100 tol = 0.00000001;
    cpp_bin_float_100 error = tol+1;

    cpp_bin_float_100 sk = 0.0;
    cpp_bin_float_100 sk1 = 0.0;
    cpp_bin_float_100 num = 0.0;
    cpp_bin_float_100 den = 0.0;
    while (k <= iterMax){
        num = factor(2.0*k)*pow(x, 2*k+1);
        den = (2.0*k + 1)*pow(2, 2*k)*pow(factor(1.0*k), 2);
        sk1 = num*div_tt(den) + sk;
        error = abs((sk1-sk)/sk1);
        if (error < tol) break;
        sk = sk1;
        k++;
    }
    return sk1;
}

cpp_bin_float_100 sqrt_t(cpp_bin_float_100 x){
    /**
     * @Brief Función para  aproximar la raíz cuadrada de un x dado utilizando
     * el método Babilonio
     * @param x es el subíndice del radical
     * @return √x
     */
    int k = 0;
    int iterMax = 5000;
    cpp_bin_float_100 tol = 0.00000001;
    cpp_bin_float_100 error = tol+1;

    cpp_bin_float_100 sk = 600.00;
    cpp_bin_float_100 sk1 = 0.0;
    cpp_bin_float_100 a = 0.0;
    cpp_bin_float_100 b = 0.0;
    if (0.0 < x){
        while (k <= iterMax){
            a = 0.5*sk;
            b = 0.5*x*div_tt(sk);
            sk1 = a + b;
            error = abs((sk1-sk)/sk1);
            if (error < tol) break;
            sk = sk1;
            k++;
        }
        return sk1;
    }
    else cout<< "Error. Ingrese un numero mayor a cero."<<endl;
    return 0.0;
}

void example(){
    /**
     * @Brief Función de prueba para aproximar el valor de la operación
     * ejemplo de la especificación de la parte 1 en este trabajo
     *              _________________
     * total =  \3/ sen(3/7) + ln(2)    + log_pi(e^⁻1)
     *             ------------------
     *                   cos(√2)
     *
     */
    double sen =  (double)sin_t(3.0*div_tt(7.0));
    double ln = (double)ln_t(2.00);
    double sqrt = (double)sqrt_t(2.00);
    double cos = (double)cos_t(sqrt);
    double euler = (double)exp_t(-1.0);
    double pii = (double)pi();
    double total = (double)(root_t(sen + ln, 3)*div_tt(cos) + log_t(euler, pii));
    cout <<"Resultado: " << total << endl;
}