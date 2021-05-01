# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 16:05:10 2018

@author: UO260284
"""

#aplicacion del metodo disparo a las edos y sedos con condiciones de contorno
#tratamos de reescribir el problema de contorno como uno de valor inicial
#nos depende la exactitud del metodo de valor inicial que usemos

import numpy as np
import matplotlib.pyplot as plt

def disparo(F,ab,cc,iter=50,tol=1E-6,mi=[0,1]):         #funcion, intervalo, condiciones de contorno, tolerancia del error, pendiente inicial 
    n=30
    yo=cc[0]
    y1=cc[1]
    a=ab[0]
    b=ab[1]
    T=np.linspace(0,1,n+1)
    h=(b-a)/float(n)         #paso a usar en el metodo de RK
    W=[]            #lista de valores finales anteriores 
    for j in range(0,iter):         # numero de iteraciones que vamos a hacer hasta obtener la gamma
        ya=np.array([yo,j])         #para las pendientes 0 y 1 lo definimos asi, tal cual empieza el bucle
        if j>1:
            g=mi[j-1]+(y1-W[j-1])*(mi[j-1]-mi[j-2])/(W[j-1]-W[j-2])            #para los valores de j mayores que 1, el valor de la derivada sera g
            ya=np.array([yo,g])
            mi.append(g)            #engadimos en la lista los valores que obtengamos de gamma, para usarlos despues 
        wi=ya.reshape((len(ya),1))
        resu=wi         #matriz de resultados 
        for i in range(0,n):
            k1=F(T[i],wi)             #coeficientes de orden 4
            k2=F(T[i]+0.5*h,wi+0.5*h*k1)
            k3=F(T[i]+0.5*h,wi+0.5*h*k2)
            k4=F(T[i]+h,wi+h*k3)
            wi1=wi+1/6.*h*(k1+2*k2+2*k3+k4) 
            resu=np.concatenate((resu,wi1),axis=1)
            wi=wi1
        W.append(resu[0,n])         #engadimos el valor final en el extremo
        if abs(resu[0,n]-y1)<=tol:
            return T,resu,g         #nos devuelve el soporte de puntos, de soluciones y el valor mas cercano de la derivada en el punto inicial 


ab=np.array([0,1])
cc=np.array([0,np.exp(2)])
f=lambda t,Y:np.array([Y[1],2*Y[1]-Y[0]]) 
T,resu,g=disparo(f,ab,cc)
print(g)
print(T)
X=resu[0,:]
Y=resu[1,:]
print(X)

#representamos los resultados 
plt.figure()
plt.subplot(211)
plt.plot(T,X)
plt.grid()
plt.title('Aproximacion por el metodo de disparo\n Funcion')
plt.subplot(212)
plt.plot(T,Y)
plt.title('Derivada de la Funcion')
plt.grid()
plt.show()

#probar para la siguiente funcion
