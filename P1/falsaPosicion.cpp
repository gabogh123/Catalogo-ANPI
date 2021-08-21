//Se importan las librerias a utilizar.

#include <iostream>
#include <cmath>
#include <iomanip>
//Fin de la importacion de librerias.



using namespace std;

//Se define la funcion a utilizar 

float f(double x){
    return cos(x)-x;

}

//Definicion y asignacion de valores a parametros iniciales 

double x2,x3;
int k=1;
double a = 0;
double b = 2;
int iterMax = 1000;
double tol = 0.0001;
double error = tol+1;

void falsa_posicion(float f(double x), double a, double b, double tol, int iterMax){

// Esta funcion aproxima numericamente el valor de una funcion mediante 
//el metodo de la falsa posicion
// Parametros de entrada: f=funcion a aproximar
                        //a,b= valores iniciales del intervalo
                        //tol= tolerancia
                        //iterMax= cantidad maxima de iteraciones
// Parametros de salida:  x3=aproximacion numerica de la funcion
                        //error=error del metodo iterativo


    while (error>tol &&  k<iterMax){ //Se definen las condiciones de parada

        if (f(a)*f(b)<0){ //Se verifica el teorema de Bolzano

            //Se define el valor x2 y se le asigna el valor
            x2 = b -( (f(b)*(b-a)) / (f(b)-f(a))); 

            if(f(x2)*f(a)<0){ //Primer subintervalo del teorema de Bolzano

                b=x2; //Se asigna el valor del x2
                
            //Se realiza la siguiente iteracion y se le asigna el valor a x3
                x3 = x2-( (f(x2)*(x2-a)) / (f(x2)-f(a))); 


            }
            else{

                a=x2; //Segundo subintervalo del teorema de Bolzano


        //Se realiza la siguiente iteracion y se le asigna el valor a x3
                x3 = x2-( (f(x2)*(x2-b)) / (f(x2)-f(b))); 

            }
  

        }

        error = abs(f(x3)); //Calculo del error
        k++; //Se aumenta el contador para la siguiente iteracion

    }
    
    //Se retornan la aproximacion y el error
    cout<<"Aprox: "<<x3<<"\t"<<"  Error: "<<error<<endl; 

}

//Funcion main que llama a la funcion secante

int main(){

    falsa_posicion(f, a, b, tol, iterMax);

    return 0;

}
