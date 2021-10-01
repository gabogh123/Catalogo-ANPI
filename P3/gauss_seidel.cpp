#include <iostream>
#include <armadillo>
#include "mgl2/mgl.h"

using namespace std;
using namespace arma;


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

void gauss_seidel(mat A, mat b, mat xk, double tol, double iterMax){
    /**
        Función que calcula la solución al sistema lineal Ax = b mediante el
        método de factorización Gauss Seidel
        Parámetros de entrada:
            A: Matriz de dimensiones nxn
            b: vector xi elementos, donde i es la cantidad de incógnitas
            xk: vector xk de valor inicial
            tol: tolerancia
            iterMax: iteraciones maximas
        Salida: vector con el xk final, el numero de iteraciones y el error
        **/
    double det_A = det(A);
    
    if (det_A != 0 && A.is_square()){
        double error;
        vector<mat> res1;
        vector<double> res2;
        
        // Paso 1
        mat L = trimatl(A,-1);
        mat D = diagmat(A);
        mat U = trimatu(A,1);
        
        // Paso 2
        mat LD = L + D;
        mat y = sust_adelan(LD, b);
        
        
        //Variables para graficar
        mglData dataX(40);
        mglData dataY(40);
        vector<int> iteraciones;
        vector<double> errores;
        
        mat zk;
        mat xk_n;
        int k = 0;
        
        // Paso 3
        while (k < iterMax) {
            // Paso 3.1
            zk = sust_adelan(LD, xk);
            // Paso 3.2
            xk_n = -zk + y;
            // Paso 3.3
            error = norm(A * xk_n - b, 2);
            // Paso 3.4
            if (error < tol) {
                break;
            }
            dataX.a[k] = k;
            dataY.a[k] = error;
            k++;
            
        }
        cout << "xk: " << xk_n << endl;
        cout << "Iteraciones: " << k << endl;
        cout << "Error: " << error << endl;
        
        mglGraph gr;
        gr.Title("Pseudoinversa");
        gr.SetOrigin(0,0);
        gr.SetRanges(0,iter_max+1,dataY[sizeof dataY],dataY[0]); // rangos de los ejes
        // etiquetas de los ejes
        gr.Label('x',"Iteraciones");
        gr.Label('y',"Error");
        gr.Plot(dataX,dataY ,"o!r"); // plot de la grafica con los valores de los vectores
        gr.Axis();
        gr.Grid();
        
        gr.WriteFrame("Pseudoinversa.png");
        res1.push_back(xk_n);
        res2.push_back(k);
        res2.push_back(error);
        
    } else {
        cout << "La matriz no es invertible" << endl;
    }
    
   
    
    
   
}

int main() {
    mat A={{5,1,1},{1,5,1},{1,1,5}};
    mat b = {7,7,7};
    b = b.t();
    
    double tol = 1e-15;
    double iterMax = 1000;
    mat xk = {1,1,1};
    xk = xk.t();
    
    gauss_seidel(A, b, xk, tol, iterMax);

  return 0;
}
