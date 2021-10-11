#include <iostream>
#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;



vector<mat> get_X0X1(mat A){
    vector<mat> initial_mats;
    mat A_t = A.t();
    double a1 = 5 * 1e-10;
    double a2 = 2 * 1e-11;
    mat Xk_anterior = a1 * A_t; // X0
    mat Xk = a2 * A_t; // X1
    initial_mats.push_back(Xk_anterior);
    initial_mats.push_back(Xk);
    
    return initial_mats;
}

mat pseudo_inv(mat A, mat b) {
    /**
        Método iterativo para aproximar la pseudoinversa de una matriz A, el cual se deduce del método de la secante para aproximar un cero de la función no lineal.
        Parámetros de entrada:
            A: Matriz de dimensiones 45x30
            b: vector xi elementos, donde i es la cantidad de incógnitas
        Salida: vector xk
        **/
    
    int iterMax = 1000;
    double tol = 1e-5;
    
    
    vector<mat> initial_mats = get_X0X1(A);
    mat Xk = initial_mats.at(1); // Xk
    mat Xk_anterior = initial_mats.at(0); // Xk - 1
    mat Xk_siguiente; // Xk + 1

    
    mat to_norm;
    
    double error = tol + 1.0;
    
    for (int k=1; k<iterMax; k++){
        
        Xk_siguiente = Xk_anterior + Xk - (Xk_anterior * A * Xk); // Calculo de Xk + 1
        
        to_norm = (A * Xk * A) - A; // Se obtiene la matriz por normalizar
        
        error = norm(to_norm, "fro"); //Se calcula la norma de frobenius
        
        if (error < tol){ // Condicion de parada
            break;
        }
        
        // Se redefinen las variables para la siguiente iteración
        Xk_anterior = Xk;
        Xk = Xk_siguiente;
    }
  
    return Xk;
}

/**
----------------------------------------------------APLICACIÒN----------------------------------------------
*/





mat get_c(mat delta) {
    /**
        Esta función obtiene el vector c el cual permite conocer si los valores del delta se salen del rango
        Parámetros de entrada:
            delta: vector delta
            
        Salida: vector c
        **/
    mat c = delta;
    for (int i=0; i<delta.size(); i++) {
        if (-0.75 > delta(i) > 0.75){ // Revisa si esta fuera del rango
            c(i) = 0.75;
        } else {
            c(i) = 0;
        }
    }
    return c;
}



mat parte3_p2(mat A, mat b) {
    
    mat pinv_A = pseudo_inv(A,b); // Calcula la pinv
    mat delta = pinv_A * b; // Calcula el delta inicial
    mat c = get_c(delta); // Obtiene el vector c
    mat A_mod = A;
    
    for (int i=0; i<A_mod.n_rows; i++) { // Se crea la matriz A modificada con la ultima columna con valores iguales a cero
            for (int j=0; j<A_mod.n_cols; j++) {
                if (j == A_mod.n_cols-1) {
                    A_mod(i,j) = 0;
                }
            }
        }
    
    mat pinv_A_mod = pinv_A;
    
    for (int k=1; k<1000; k++){
        
        if (c.is_zero()) { // Se revisa si el vector c no tiene ningun valor diferente a cero
            break;
        } else { // Se calcula de nuevo el delta y el c
            pinv_A_mod = pseudo_inv(A_mod,b);
            delta = -c + pinv_A_mod * (b + (A * c));
            get_c(delta);
        }
        
    }
    cout << delta << endl;
    return delta;
}

int main() {
    mat A={{2,-2,-2,-1},{1,1,-3,-2},{2,-2,-1,-1}};
    mat b={0.5,1,1};
    b = b.t();
    parte3_p2(A, b);
   
    return 0;
}
