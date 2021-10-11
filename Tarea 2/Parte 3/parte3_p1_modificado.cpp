#include <iostream>
#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

mat get_A(){
	int iterMax = 1000;
    int tol = 1e-15;
    
    mat A(45,30);
    
    for (int i=0; i<45; i++) {
        for (int j=0; j<30; j++) {
          
            
            A(i,j) = (pow(i+1,2) + pow(j+1,2));
        }
    }
    mat A_t = A.t();

    return A;
}

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

mat parte3_p1(mat A) {
    
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
            cout << "break" << endl;
            cout << "Error: " << error;
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



int main() {
    mat A = get_A();
    mat pinv_A = parte3_p1(A);
    cout<<"Xk final: \n"<<pinv_A<<endl;
    return 0;
}
