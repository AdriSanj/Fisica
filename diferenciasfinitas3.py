# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 12:13:17 2018

@author: UO260284
"""

#metodo de diferencias finitas con 3 puntos
#condiciones de contorno de tipo dirichlet


import numpy as np
import matplotlib.pyplot as plt
from math import *


def bvpdfm(F,ab,n,cc):
    a=ab[0]             #punto inicial
    b=ab[1]             #punto final
    x=np.linspace(a,b,n+1)               #array del soporte, metemos el numero de subintervalos
    h=(b-a)/n           #interespaciado
    fo=cc[0]            
    ff=cc[1]             
    B=F(x)[2]
    B[0]=fo         #primer elemento del array lo sustituimos por su valor inicial 
    B[n]=ff         #engadimos el ultimo coeficiente, longitud 5 del vector,lo mismo pero con el final
    P=F(x)[0]       #vector de valores de P
    Q=F(x)[1]       #vector de valores de Q
    SOL=np.zeros((n+1,n+1))         #matriz que tendremos que invertir
    SOL[0,0]=1
    SOL[n,n]=1
    for j in range(1,n):
        SOL[j,j-1]=1/(h**2)-1/(2*h)*P[j]
        SOL[j,j]=-2/(h**2)+Q[j]
        SOL[j,j+1]=1/h**2-1/(2*h)*P[j]
    SOL=np.linalg.solve(SOL,B)
    return B,x,P,Q,SOL
    
    
    
    
    
    
    
F= lambda x: [2.*x,x+np.ones(len(x)),3*np.ones(len(x))]

ab=np.array([-1,1])
cc=np.array([4,6])
n=4
B,x,P,Q,s=bvpdfm(F,ab,n,cc)


print(x)
#print(P)
#print(Q)
print(s)
plt.plot(x,s,'-r',x,s)
plt.grid()
plt.title('Comparacion')
plt.show()
