#include <iostream>
#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;



mat potencia_inversa(mat A,mat vec_inicial,int iter_max) {
    
    
    mat iteraciones;
    mat errores;
    double tol = 1e-5;
    
    
    mat yk1;
    float ck1;
    mat xk1;
    mat xk = vec_inicial;
    double error = 1.0 + tol;
    
    for (int k=1; k<iter_max; k++){
        
        yk1 = solve(A, xk);
        
        ck1 = norm(yk1, "inf");
        
        xk1 = (1/ck1) * yk1;
        
        error = norm(xk1-xk);
        
        xk = xk1;
        
        if (error < tol) {
            cout << "Final error: "<< error << endl;
            cout << "Final k: "<< k << endl;
            cout << "Final ck: "<< ck1 << endl;
            cout << "Final xk: "<< xk << endl;
            break;
        }
        cout << "Error: "<< error << endl;
        cout << "k: "<< k << endl;
        cout << "ck: "<< ck1 << endl;
        cout << "xk: "<< xk << endl;
        
    }
    
    return A;
}

int main() {
    mat A={{3,-1,0,},
           {-1,2,-1},
           {0,-1,3}};
    
    mat b={1,1,1};
    b = b.t();
    potencia_inversa(A,b,20);
   
    return 0;
}

