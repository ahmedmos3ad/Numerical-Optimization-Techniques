import sympy
import numpy

#defining x,y,l as symbols
x,y,Lamda = sympy.symbols('x y Lamda')
#function given in lecture
function=x-y+2*x**2 + 2*x*y + y**2

#Powell's Method
def Powell(X1,e,NoOfIterations):
    print(f"POWELL'S METHOD")
    s = numpy.array([1,0])
    X = numpy.zeros((NoOfIterations+1,2))
    X[0]=X1
    LAMBDA = numpy.zeros(NoOfIterations)
    for i in range(NoOfIterations):
        if i != 3:
            # first cycle
            #s1
            if s[1] == 0:
                s = numpy.array([0,1])
            #s2
            else:
                s = numpy.array([1,0])
        # second cycle
        else:
            #generate first pattern direction
            s = X[3] - X[1]
        #printing values of s to console
        print("S",s)
        # substituting in function to get f+ and f-
        # we calculate f to compare with fpos and fneg later on to determine which pole it decreases with
        # substitute with values obtained from X
        f = function.subs({x:X[i,0],y:X[i,1]})
        fpos = function.subs({x:X[i,0] +  e * s[0],y:X[i,1] +  e * s[1]})
        fneg = function.subs({x:X[i,0] -  e * s[0],y:X[i,1] -  e * s[1]})
        print(f'f = {round(f,4)}, fpos = {round(fpos,4)}, fneg = {round(fneg,4)}')
        # if it decreases along -s multiply s by -1
        if fneg < f:
            s*=-1
            print("f decreases along -S")
        else:
            print("f decreases along +S")
        function_lamda = function.subs({x:X[i,0] + s[0]* Lamda, y: X[i,1] + s[1]*Lamda})
        print("f(",round(X[i,0],4) + s[0]* Lamda,",", round(X[i,1],4) + round(s[1],4)*Lamda,")=",function_lamda)
        dfdLamda = sympy.Derivative(function_lamda,Lamda).doit()
        print("dF/dLamda =",dfdLamda)
        #solve derv lamda
        LAMBDA[i] = sympy.solve(dfdLamda, Lamda)[0]
        print("Lamda =",LAMBDA[i])
        X[i+1] = X[i] + LAMBDA[i]*s
        #X2+++
        print('X',i+2,"=",X[i+1],"\n")
    print(f"==============================================================================================================")

#Steepest Descent (Cauchy’s) Method
def Steepest(X1,NoOfIterations):
    print(f"STEEPEST DESCENT METHOD")
    #def with respect to x
    dfdx = sympy.Derivative(function,x).doit()
    #def with respect to y
    dfdy = sympy.Derivative(function,y).doit()
    
    print("deltaF = {",dfdx," , ",dfdy,"}\n")
    deltaF = numpy.zeros((NoOfIterations,2))
    X = X1
    LAMBDA = numpy.zeros(NoOfIterations)
    for i in range(NoOfIterations):
        #substitute in the derived functions
        deltaF[i,0] = dfdx.subs({x:X[0],y:X[1]})
        deltaF[i,1] = dfdy.subs({x:X[0],y:X[1]})
        print("deltaF",i+1," = ",deltaF[i])
        function_lamda = function.subs({x:X[0] - deltaF[i,0]* Lamda, y: X[1] - deltaF[i,1]* Lamda})
        print("f(",X[0] - deltaF[i,0]* Lamda,",", X[1] - deltaF[i,1]* Lamda,")=",function_lamda)
        dfdLamda = sympy.Derivative(function_lamda,Lamda).doit()
        print("df/dLamda =",dfdLamda)
        LAMBDA[i] = sympy.solve(dfdLamda, Lamda)[0]
        print("Lamda =",LAMBDA[i])
        X = X - LAMBDA[i]*deltaF[i]
        #X1++
        print("X",i+1,"=",X,"\n")


#epsilon=0.01
#starting from x1 = [0,0]
#number of iterations 5
Powell([0,0], 0.01, 5)
#number of iterations for steepest descent is 4
#starting x is the same value
Steepest([0,0],4)