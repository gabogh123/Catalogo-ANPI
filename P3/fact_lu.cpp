#include <iostream>
#include <armadillo>
#include <limits>
#include <vector>

using namespace std;
using namespace arma;


vector<mat> fact_lu(mat A, mat b){
    /**
    Función que calcula la solución al sistema lineal Ax = b mediante el
    método de factorización LU
    Parámetros de entrada:
        A: Matriz de dimensiones nxn
        b: vector xi elementos, donde i es la cantidad de incógnitas
    Salida: vector con las matrices triangular inferior L y triangular superior u
    **/
    vector<mat> LU;		
    int n = A.n_rows;
    mat L = zeros<mat>(n,n);
    mat U = zeros<mat>(n,n);
    double aux = 0.00;
    for (int j = 0; j < n; j++){  
        for (int i = 0; i < n; i++){  	
            if (i > j && A(i,j) != numeric_limits<double>::epsilon()){
                aux = -1*A(i, j)/A(j, j);
                L(i,j) = -1*aux;
                for(int k = 0; k < n; k++){
                     A(i,k) = A(i,k) + A(j,k) * aux;
                }
            }
            else if(i==j){
                    L(i,j)  = 1;
            }  		
        }
    }
    LU.push_back(L);
    LU.push_back(A); 
    return LU;
}

mat sust_atras(mat A, mat b){
    /**
    Función que realiza las operaciones para resolver un
    el sistema triangular superior mediante el método
    de sustitución hacia atrás.
    Parametros de entrada:
        A: matriz triangular superior de dimensiones nxn
        b: vector modificado según las operaciones de fila
    Salida: vector solución del sistema    
    **/
    int m = A.n_rows;
    mat X = zeros<mat>(m ,1);
    int i = m - 1;

    double aux = 0.00;
    while (i >= 0){
        aux = b(i);
        int j = m - 1;
        while (j >= 0){
            if (i == j){
                X(j,0) = aux/A(i, j);
                break;
            }
            else{
                aux -= A(i, j)*X(j, 0);
            }
            j -= 1;
        }
    i -= 1;
    }	
    return X;
}
mat sust_adelan(mat A, mat b){
    /**
    
    Función que realiza las operaciones para resolver un
    el sistema triangular inferior mediante el método
    de sustitución hacia atrás.
    Parámetros de entrada:
        A: matriz triangular inferior de dimensiones nxn
        b: vector modificado según las operaciones de fila  
    Salida: vector solución
    
    **/
    int m = A.n_rows;
    mat X = zeros<mat>(m,1);
    double aux = 0.00;
    for (int i = 0; i < m; i++){
        aux = b(i);
        for (int j = 0; j < m; j++){
            if (i == j){
                X(j, 0) = aux/A(i, j); 
                break; 
            } 		    
            else{
                aux -= A(i, j)*X(j);
                }		        
        }	
    }
    return X;
}



bool esValida(mat A){
	int n = A.n_rows;
	int j = 0;
	for (int i = 0; i < n; i++){
		mat cols = A.cols(0, j);
		mat rows = A.rows(0,i);
		mat submat = A( span(0, i), span(0, j)); 
		cout << submat << endl; 
		j += 1;
		double deter = (double)det(submat);
		
		if (deter == numeric_limits<double>::epsilon()){
			return false;
			} 
	}
	return true;
}

int main(){
    mat A = {{4, -2, 1}, {20, -7, 12 }, {-8, 13, 17}};
    mat b1 = {11, 70, 17};
    mat b = b1.t();
    cout << "Sistema Ax = b a resolver:  \n" << endl;
    cout << "Matriz A: \n "<< A << "\n" << endl;
    cout << "\n Vector b: "<< b1 << "\n" <<  endl;
    vector<mat> LU = fact_lu(A, b);
    mat L = LU.at(0);
    mat U = LU.at(1);
    mat y = sust_adelan(L, b);
    cout << "\n Vector solución X \n "<< endl;
    mat X = sust_atras(U, y);
    cout << X << endl;  

    return 0;
}
