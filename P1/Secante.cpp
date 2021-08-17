#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>


using namespace std;

float f(double x){
    return pow(x,6)+pow(x,2)-1;
}

double x2;
int k=1;
int iterMax = 1000;
double x0=0;
double x1=1;
double tol=0.0001;
double error = tol+1;

void secante(float f(double x), double x0, double x1, double tol, int iterMax){


    while(error>tol && k<iterMax){

        x2=x1-(f(x1)*(x1-x0)/(f(x1)-f(x0)));
        x0=x1;
        x1=x2;
        
        error = abs(f(x2));
        k++;

       // aprox = {x2}+aprox;
       // er = {error}+er;

    }

    cout<<"Aprox: "<<x2<<"\t"<<"  Error: "<<error<<endl;


}

int main(){

	secante(f, x0, x1, tol, iterMax);

	return 0;
}
