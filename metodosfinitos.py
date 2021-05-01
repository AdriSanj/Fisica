# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 16:05:36 2018

@author: UO260284
"""

#metodo de elementos finitos
#funciones lineales ek=[x(k-1),xk]    hk=xk-(x(k-1))
#v(k-1)=1/hk*(xk-x)    /   vk=1/hk*(x-x(k-1))

import numpy as np
import matplotlib.pyplot as plt
from math import *
from scipy.integrate import quad
#p,q,f funciones numericas
def bvpfeme(p,q,f,a,b):         #matriz caja 2x2, la integral usar quad
    h=b-a  #h,b y a seran variables globales para las funciones v0 y v1
    def v0(x):           #xka es x(k-1)
        return 1/h*(b-x)        #definimos asi las funciones triangulares para hacer los calculos    EN EL EXAMEN PONEMOS LAS QUE MANDEN
    def v1(x):
        return 1/h*(x-a)
    A=np.zeros((2,2))       #matriz de coeficientes
    B=np.zeros(2)
    def k1(x):       #no podemos multiplicar funciones numericas dentro de los quad, las definimos anteriormente 
        return f(x)*v0(x)
    def k2(x):
        return f(x)*v1(x)
    def m11(x):
        return p(x)/(h**2)+q(x)*(v0(x)**2)
    def m12(x):
        return q(x)*(v0(x)*v1(x))-p(x)/(h**2)
    def m22(x):
        return q(x)*(v1(x)**2)+p(x)/(h**2)
    A[0,0]=quad(m11,a,b)[0]     #evaluamos la integral de cada elemento
    A[0,1]=quad(m12,a,b)[0]
    A[1,0]=quad(m12,a,b)[0]
    A[1,1]=quad(m22,a,b)[0]
    B[0]=quad(k1,a,b)[0]        #vector 
    B[1]=quad(k2,a,b)[0]
    return A,B
    

def bvpfem(p,q,f,x,cc):
    SOL=np.zeros((len(x),len(x)))
    EFE=np.zeros(len(x))        #hacemos un bucle, SOL es la matriz donde incluimos los resultados, y EFE el vector para resolver el sistema
    M=1e10      #valor de aproximacion dado porque si
    for i in range(0,len(x)-1):
        A,B=bvpfeme(p,q,f,x[i],x[i+1])      #recorremos cada subintervalo
        SOL[i,i]=SOL[i,i]+A[0,0]            #sumamos cada elemento mas el siguiente que se le anhade. en los extremos solo hay un numero, en el resto de la diagonal sumas
        SOL[i,i+1]=SOL[i,i+1]+A[0,1]
        SOL[i+1,i]=SOL[i+1,i]+A[1,0]
        SOL[i+1,i+1]=SOL[i+1,i+1]+A[1,1]
        EFE[i]=EFE[i]+B[0]
        EFE[i+1]=EFE[i+1]+B[1]
    SOL[0,0]=M          #multiplicamos por la m en los extremos de la matriz y el vector
    SOL[len(x)-1,len(x)-1]=M
    EFE[0]=M*cc[0]
    EFE[len(x)-1]=M*cc[1]
    resultado=np.linalg.solve(SOL,EFE)          #resolvemos el sistema
    return resultado
    
def f(x):
    return np.exp(x)+np.e/16.*(2+3*x-2*x**2)
def p(x):
    return 9/2.-x
def q(x):
    return 3/2.-x

soporte=np.array([0.0,0.2,0.4,0.6,0.8,1.0])
cc=np.array([0.0,0.0])


resu=bvpfem(p,q,f,soporte,cc)
print resu


