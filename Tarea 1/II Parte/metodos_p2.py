import numpy as np
import sympy as sp

def newton_H_m1(f,x0,tol,iterMax):

    x=sp.Symbol('x')

    
    f1=sp.sympify(f)
    df1 = sp.diff(f1,x)
    
    h = 1
    h1 = sp.sympify(h)

    

    er=[]
    err=tol+1
    k=0
    xk=x0

    while err>tol and k<iterMax:
        k=k+1
        
        n=sp.N(f1.subs(x,xk))
        d=sp.N(df1.subs(x,xk))
        
        xk=xk-(h*n/d)

        
        err=abs(f1.subs(x,xk))
        


    print( "Aproximación: " + str(xk))
    print ("Error: " + str(err))


        
    

####################################################################################################################################

def newton_H_m2(f,x0,tol,iterMax):

    x = sp.Symbol('x')
    B=2

    f1 = sp.sympify(f)
    df1 = sp.diff(f1,x)

    err = tol+1
    k=0
    xk = x0

    while err>tol and k<iterMax:


        h = 1/(1+B*(sp.N(f1.subs(x,xk))/sp.N(df1.subs(x,xk))))

        n = sp.N(f1.subs(x,xk))
        d = sp.N(df1.subs(x,xk))

        xk = xk - (h*(n/d))

        err = abs(f1.subs(x,xk))

        k=k+1

    return ['Aproximación: ']+[xk]+ [   'Error: ']+ [err]




###########################################################################################################################


def newton_G_m1(f,x0,tol,iterMax):

    x = sp.Symbol('x')

    f1 = sp.sympify(f)
    df1 = sp.diff(f1,x)
    df2 = sp.diff(df1,x)

    err = tol+1
    k = 0
    xk = x0

    while err>tol and k<iterMax:

        g = 2/( 2 - ((sp.N(f1.subs(x,xk)))* (sp.N(df2.subs(x,xk))) / (sp.N(df1.subs(x,xk)))**2 ))

        n = sp.N(f1.subs(x,xk))
        d = sp.N(df1.subs(x,xk))

        
        xk = xk - (g*(n/d))

        err = abs(f1.subs(x,xk))

        
        k=k+1

    return ['Aproximación: ']+[xk]+ [   'Error: ']+ [err]
    


#############################################################################################################################


def newton_G_m2(f,x0,tol,iterMax):

    
    x = sp.Symbol('x')

    f1 = sp.sympify(f)
    df1 = sp.diff(f1,x)
    df2 = sp.diff(df1,x)

    err = tol+1
    k = 0
    xk = x0

    while err>tol and k<iterMax:
    

        g = 1 + ( ((sp.N(f1.subs(x,xk)))* (sp.N(df2.subs(x,xk))) / (sp.N(df1.subs(x,xk)))**2 ) /2 )

        n = sp.N(f1.subs(x,xk))
        d = sp.N(df1.subs(x,xk))

        xk = xk - (g*(n/d))

        err = abs(f1.subs(x,xk))
        
        k=k+1

    return ['Aproximación: ']+[xk]+ [   'Error: ']+ [err]



