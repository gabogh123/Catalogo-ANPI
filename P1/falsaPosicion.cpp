#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double f(double x){
    return cos(x)-x;

}

double x2,x3;
int k=1;
double a = 0;
double b = 2;
int iterMax = 1000;
double tol = 0.0001;
double error = tol+1;

void falsa_posicion(){

    while (error>tol &&  k<iterMax){

        if (f(a)*f(b)<0){

            x2 = b -( (f(b)*(b-a)) / (f(b)-f(a)));

            if(f(x2)*f(a)<0){

                b=x2;
                
                x3 = x2-( (f(x2)*(x2-a)) / (f(x2)-f(a)));


            }
            else{

                a=x2;

                x3 = x2-( (f(x2)*(x2-b)) / (f(x2)-f(b))); 

            }
  

        }

        error = abs(f(x3));
        k++;

    }
    
    cout<<"Aprox: "<<x3<<"\t"<<"  Error: "<<error<<endl;

}

int main(){

    falsa_posicion();

    return 0;

}
