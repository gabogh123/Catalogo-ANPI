//Se importan las librerias a utilizar.
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
//Fin de la importacion de librerias.


using namespace std;
//Se define la funcion a utilizar 
float f(double x){
    return pow(x,6)+pow(x,2)-1;
}

//Definicion y asignacion de valores a parametros iniciales 
double x2;
int k=1;
int iterMax = 1000;
double x0=0;
double x1=1;
double tol=0.0001;
double error = tol+1;

void secante(float f(double x), double x0, double x1, double tol, int iterMax){

// Esta funcion aproxima numericamente el valor de una funcion
//mediante el metodo de la secante
// Parametros de entrada: f=funcion a aproximar
                        //x0,x1= valores iniciales del intervalo
                        //tol= tolerancia
                        //iterMax= cantidad maxima de iteraciones
// Parametros de salida:  x2=aproximacion numerica de la funcion
                        //error=error del metodo iterativo

    while(error>tol && k<iterMax){ //Se definen las condiciones de parada

        x2=x1-(f(x1)*(x1-x0)/(f(x1)-f(x0)));
        x0=x1; //Redefinicion de x0
        x1=x2; //Redefinicion de x1
        
        error = abs(f(x2)); //Calculo del error
        k++; //Se aumenta el contador para la proxima iteracion

       // aprox = {x2}+aprox;
       // er = {error}+er;

    }

    cout<<"Aprox: "<<x2<<"\t"<<"  Error: "<<error<<endl; //Se imprimen los valores


}


//Funcion main que llama a la funcion secante
int main(){

	secante(f, x0, x1, tol, iterMax);

	return 0;
}