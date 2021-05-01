# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 16:04:33 2018

@author: UO260284
"""

#metodos de runge-kutta
#Yo= alfa, Ym=Y(m-1)+h*suma(i,n,Cj*Kj)
#Kj=f(t(i-1)+p*h,Y(i-1)+h*suma(1,s,Qj.l*Kl))


#programar una funcion de runge-kutta para resolver una EDO de PVI

import numpy as np
import matplotlib.pyplot as plt

def rungeke(F,ab,n,ya,orden): #F=y', ab intervalo, n numero de subintervalos
    wi=ya      #valor inicial de las y
    a=ab[0]         #punto inicial y final
    b=ab[1]     
    h=(b-a)/float(n)            #definimos el paso
    ti=a       #valor inicial de las t
    T=[ti]
    Y=[wi]              #metemos el valor inicial en la lista porque en el for no lo coje 
    if orden==2:
        for i in range(0,n):
            k1=F(ti,wi)         #coeficiente k1
            k2=F(ti+0.5*h,wi+0.5*h*k1)      #coeficiente k2
            wi1=wi+h*k2          #valor de la Y en el punto ti+h
            Y.append(wi1)
            ti+=h
            T.append(ti)
            wi=wi1       #el valor siguiente de wi sera el calculado
    elif orden==4:
        for i in range(0,n):
            k1=F(ti,wi)             #coeficientes de orden 4
            k2=F(ti+0.5*h,wi+0.5*h*k1)
            k3=F(ti+0.5*h,wi+0.5*h*k2)
            k4=F(ti+h,wi+h*k3)
            wi1=wi+1/6.*h*(k1+2*k2+2*k3+k4)     #siguiente valor de la y
            Y.append(wi1)
            ti+=h
            T.append(ti) 
            wi=wi1 
    return T,Y

def edo(t,y):       #funÇao na que definimos a edo
    return -y
ab=np.array([0,5])

T1,Y1=rungeke(edo,ab,10,1,2)
T2,Y2=rungeke(edo,ab,10,1,4)
plt.plot(T1,Y1,'r.',T2,Y2,'b*')
plt.grid()
plt.show()


















