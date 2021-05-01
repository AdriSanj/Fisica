#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 16:04:43 2018

@author: uo260284
"""

from math import *
#FORMULA DE REGLAS CERRADAS
def newton_cotes_cerrada(a,b,f):           #antes de poner el error
    h1=(b-a)                #a y b los extremos del intervalo, las hi son los espacios dependiendo de cuantos puntos intermedios tratemos
    h2=(b-a)/2.
    h3=(b-a)/3.
    h4=(b-a)/4
    R1=h1/2.*(f(a)+f(b))           #trapecio
    R2=h2/3.*(f(a)+4*f(a+h2)+f(b))     #simpson
    R3=3*h3/8.*(f(a)+3*f(a+h3)+3*f(b-h3)+f(b))        #simpson 3/8
    R4=2*h4/45.*(7*f(a)+32*f(a+h4)+12*f(a+2*h4)+32*f(b-h4)+f(b))            #boole
    return R1,R2,R3,R4

def f1(x):          #[0,5]
    f=1/(1+x**2)
    return f

def f2(x):          #[0,4]
    f=e**(-x**2)
    return f

def f3(x):          #[-1,1]
    f=sqrt(1+x**2)
    return f
    
s1,s2,s3,s4=newton_cotes_cerrada(0,5,f1)

print('Aproximaciones cerrada')
print(s1)
print(s2)
print(s3)
print(s4)

print(' ')

#multiplicadores de la derivada del error, pero no juega parte a la hora de hacer la aproximacion

#-((h1)**3)/12. 
#-((h2)**5)/90. 
#-3*((h3)**5)/80.  
#-8*((h4)**7)/945.

#FORMULA DE REGLAS ABIERTAS

def newton_cotes_abierta(a,b,f):
    h0=(b-a)/2          #n=0
    h1=(b-a)/3          #n=1
    h2=(b-a)/4          #n=2
    h3=(b-a)/5          #n=3
    R1=2*h0*f(a+h0)         #punto medio
    R2=3*h1/2.*(f(a+h1)+f(a+2*h1))          #trapecio abierta
    R3=4*h2/3.*(2*f(a+h2)-f(a+2*h2)+2*f(a+3*h2))        #Milne
    R4=5*h3/24.*(11*f(a+h3)+f(a+2*h3)+f(a+3*h3)+11*f(a+4*h3))           #wilms
    return R1,R2,R3,R4

w1,w2,w3,w4=newton_cotes_abierta(0,5,f1)
print('Aproximaciones abiertas')
print(w1)
print(w2)
print(w3)
print(w4)

print(' ')
#diferencia de metodos
print('Diferencias entre ambos metodos')
D1=abs(s1-w1)
D2=abs(s2-w2)   
D3=abs(s3-w3)  
D4=abs(s4-w4) 

print(D1)    
print(D2)     
print(D3)     
print(D4)     
    