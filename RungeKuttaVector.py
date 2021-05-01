#runge kutta vector

import numpy as np
import matplotlib.pyplot as plt

def rungekevect(F,ab,n,ya, orden):      #F funcion, ab intervalo, n numero de puntos a evaluar, ya condiciones iniciales, orden si queremos orden 2 o 4
    wi=ya.reshape((len(ya),1))           #vector inicial, lo transponemos por que vamos a hacer una matriz concatenando todos los resultados 
    a=ab[0]
    b=ab[1]
    h=(b-a)/float(n)            #espaciado
    ti=a   #primer elemento de la variable t. en el se define la condicion inicial
    T=[ti]
    resu=wi
    if orden==2: 
        for i in range(0,n):
            k1=F(ti,wi)         #coeficiente k1
            k2=F(ti+0.5*h,wi+0.5*h*k1)      #coeficiente k2
            wi1=wi+h*k2         #calculamos asi 
            resu=np.concatenate((resu,wi1),axis=1)              #concatenamos la matriz, cada resultado wi1 sera una columna nueva 
            ti+=h
            T.append(ti)            #valores de la variable t
            wi=wi1
    elif orden==4:
        for i in range(0,n):
            k1=F(ti,wi)             #coeficientes de orden 4
            k2=F(ti+0.5*h,wi+0.5*h*k1)
            k3=F(ti+0.5*h,wi+0.5*h*k2)
            k4=F(ti+h,wi+h*k3)
            wi1=wi+1/6.*h*(k1+2*k2+2*k3+k4) 
            resu=np.concatenate((resu,wi1),axis=1)
            ti+=h
            T.append(ti)
            wi=wi1
    return T,resu
    
f=lambda t,Y:np.array([Y[1],2/(t**3+4)-Y[0]*Y[1]-9*(Y[0])**3])

ab=np.array([1,2])
ya=np.array([1/5.,7/25.])

T,resu=rungekevect(f,ab,10,ya,4)
print(T)
print(len(T))

X=resu[0,:]         #la primera fila de la matriz son los resultados en X(t), la segunda los de Y(t)
Y=resu[1,:]
print('valores de la funcion')
print(X)
print('valores de la derivada')
print(Y)
plt.figure()
plt.plot(T,X,'r',T,Y,'b')
plt.grid()
plt.show()
