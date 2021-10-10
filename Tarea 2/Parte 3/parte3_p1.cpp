#include <iostream>
#include <armadillo>


using namespace std;
using namespace arma;




mat parte3_p1() {
    
    int iterMax = 1000;
    int tol = 1e-15;
    
    mat A(45,30);
    
    for (int i=0; i<45; i++) {
        for (int j=0; j<30; j++) {
          
            
            A(i,j) = (pow(i+1,2) + pow(j+1,2));
        }
    }
    mat A_t = A.t();
    cout << "Matriz A: " << endl;
    cout << A << endl;
    
    
    
    
    double a1 = 5 * 1e-10;
    double a2 = 2 * 1e-11;
    
    //cout << "a1: " << a1 << endl;
    //cout << "a2: " << a2 << endl;
    
 
    for (int i=0; i<45; i++) {
        for (int j=0; j<30; j++) {
          
            
            A(i,j) = (pow(i+1,2) + pow(j+1,2));
        }
    }
    
    mat Xk = a1 * A_t; // Xk
    mat Xk_anterior = a2 * A_t; // Xk - 1
    mat Xk_siguiente; // Xk + 1
    
    
    //cout << "Xk - 1 inicial: " << Xk << endl;
    //cout << "Xk inicial: " << Xk_anterior << endl;
    
    
    
    mat to_norm;
    
    int error;
    
    for (int k=1; k<iterMax; k++){
        
        Xk_siguiente = Xk_anterior + Xk - (Xk_anterior * A * Xk); // Calculo de Xk + 1
        
        to_norm = (A * Xk * A) - A; // Se obtiene la matriz por normalizar
        
        error = norm(to_norm, "fro"); //Se calcula la norma de frobenius
        
        if (error < tol){ // Condicion de parada
            cout << "break" << endl;
            cout << "Error: " << error;
            break;
        }
        
        // Se redefinen las variables para la siguiente iteración
        Xk_anterior = Xk;
        Xk = Xk_siguiente;
    }
    
    
    
    
    cout << "-------------------- . -------------------- . -------------------- . -------------------- . --------------------" << endl;
    cout << "El Xk final es el siguiente: " << endl;
    cout << "-------------------- . -------------------- . -------------------- . -------------------- . --------------------" << endl;
    cout << Xk << endl;
    cout << "-------------------- . -------------------- . -------------------- . -------------------- . --------------------" << endl;
    
    return A;
}

int main() {
    parte3_p1();
    return 0;
}
