#include <iostream>
#include <armadillo>
#include "mgl2/mgl.h"

using namespace std;
using namespace arma;


mat newton_schulz(mat A, mat b, double tol, double iter_max){
    /**
        Función que calcula la solución de pseudoinversa  mediante el
        método de Newton Schultz
        Parámetros de entrada:
            A: Matriz de dimensiones nxn
            b: vector xi elementos, donde i es la cantidad de incógnitas
            tol: tolerancia
            iterMax: iteraciones maximas
            Salida: vector con el xk final, el numero de iteraciones y el error
        **/
    double error;

    // Calculos necesarios para obtener normsq (||A||2)
    mat mat_aux= A * A.t();
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, mat_aux);
    int normsq= max(eigval); // Calculo de normsq
    
    
    mat Xk = ((1+0.0)/normsq) * A.t(); // Calculo de Xk
    mat xk = Xk*b;
  
    
    int n = A.n_rows;
    mat I = eye(n,n); //matriz Identidad
    
    mat multiplicador;
    mat xk_n;
    
    //Variables para graficar
    mglData dataX(40);
    mglData dataY(40);
    vector<int> iteraciones;
    vector<double> errores;

    // Iteracion de Newton-Schultz
    for (int k=1; k<iter_max; k++){
        multiplicador = 2*I - (A*Xk);
        Xk = Xk * multiplicador; // Calculo de la aproximacion
        xk_n = Xk*b;
        error = norm(xk_n - xk,2) / norm(xk_n,2); // Calculo del error
        if (error < tol){ // Condicion de parada
            cout << "Numero de iteraciones: "<< k << endl;
            break;
        }
        xk = xk_n; // Se actualiza xk
        dataX.a[k] = k;
        dataY.a[k] = error;
    }
    cout << "Error: "<< error << endl;
    cout << "xk: "<< xk << endl;
    
    
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
    return xk;
}


mat pseudoinversa(mat A, mat b, double tol, double iter_max){
    /**
        Función que calcula la solución al sistema lineal Ax = b mediante el
        método de factorización Pseudoinversa
        Parámetros de entrada:
            A: Matriz de dimensiones nxn
            b: vector xi elementos, donde i es la cantidad de incógnitas
            tol: tolerancia
            iterMax: iteraciones maximas
        Salida: vector con el xk final, el numero de iteraciones y el error
        **/
    mat At = newton_schulz(A, b, tol, iter_max);
    
    cout << "Resultado de la aproximacion de la pseudoinversa: " << At << endl;
    return At;
}

int main() {
  mat A={{1,2,3},{1,8,9},{1,4,1},{44,0,1},{1,4,5}};  // set de la matriz a realizar el calculo
  mat b={1,1,1,1,1}; // set del vector determinos independientes

  b=b.t(); // traspuesta de b

  pseudoinversa(A, b, 1e-15, 1000);
  return 0;
}

