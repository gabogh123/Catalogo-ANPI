import numpy as np
import sympy as sp
from sympy import *
import math


def aplicacion_p2():

    d = symbols('d')

    d0 = 1
    tol = 10**-5
    error = tol+1
    iterMax = 100

    dk = d0

    k = 0

    g = ( 200*acos(d/20)) - (d*sqrt(100-((d**2)/4))) #Listo

    f = ((log(7/d,10))/( 0.01*log(10)))  + ((d*(6-d)) / (((g**2) / (2*1*(40/log(10))**2)) * ( (1/g) + (1/(pi*100)) ))) #Listo   

    df = sp.diff(f,d)

    while error>tol and k<iterMax:

        dk2 = dk -((f.subs(d,dk))/(df.subs(d,dk)))

        error = abs((dk2-dk)/dk2)

        dk = dk2

        k=k+1

        return [float(dk2),float(error)]

        

    



    
    

    
