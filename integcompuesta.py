# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 12:13:20 2018

@author: UO260284
"""

#division en subintervalos. N es el numero de subsubintervalos
#EJERCICIO 5 DE INTEGRACION
from math import *

def simpson(a,b,f,N):           #formula de simpson en subintervalos 
    #numero de subintervalos es N
    h=(b-a)/(2.*N)
    suma0=f(a)+f(b)     #suma del primer y ultimo punto
    suma1=0                 #suma de las xk impares
    suma2=0         #suma de las xk pares
    #sumatoria de k=1 hasta N
    for i in range(1,N+1):          #integramos cada subintervalo
        xk1=a+(2*i-1)*h        #termino x(2k-1) (impar)
        xk2=a+(2*i)*h          #termino x(2k)   (par)
        R1=4*f(xk1)         
        suma1+=R1           #sumamos cada punto
    for m in range(1,N):
        xk1=a+(2*m-1)*h        #termino x(2k-1) (impar)
        xk2=a+(2*m)*h          #termino x(2k)   (par)
        R2=2*f(xk2)         
        suma2+=R2
    sumatotal=h/3.*(suma1+suma2+suma0)              #esta es la formula de simpson
    return sumatotal

def trapcomp(a,b,f,N):
    h=(b-a)/N
    suma0=f(a)+f(b)     #suma de los puntos extremos del intervalo a b
    suma1=0                 #suma de los puntos 
    for i in range(1,N):
        xk=a+i*h
        suma1+=f(xk)
    sumatotal=h/2.*(suma0+2*suma1)
    return sumatotal

def f1(x):          #[0,5]
    f=1/(1+x**2)
    return f

def f2(x):          #[0,4]
    f=e**(-x**2)
    return f

def f3(x):          #[-1,1]
    f=sqrt(1+x**2)
    return f

S=simpson(-1,1,f3,5)
print(S)

T=trapcomp(-1,1,f3,5)
print(T)

