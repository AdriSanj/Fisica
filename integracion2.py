#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 16:04:43 2018

@author: uo260284
"""

from math import *

def newton_cotes(a,b,f):           #antes de poner el error
    h1=(a-b)
    h2=(a-b)/2.
    h3=(a-b)/3.
    h4=(a-b)/4
    R1=h1/2.*(f(a)+f(b))           #trapecio
    R2=h2/3.*(f(a)+4*f(a+h2)+f(b))     #simpson
    R3=3*h3/8.*(f(a)+3*f(a+h3)+3*f(b-h3)+f(b))        #simpson 3/8
    R4=2*h4/45.*(7*f(a)+32*f(a+h4)+12*f(a+2*h4)+32*f(b-h4)+f(b))            #boole
    return R1,R2,R3,R4

def f1(x):
    f=1/(1+x**2)
    return f

def f2(x):
    f=e**(-x**2)
    return f

def f3(x):
    f=sqrt(1+x**2)
    return f
    
s1,s2,s3,s4=newton_cotes(0,5,f2)

print(s1)
print(s2)
print(s3)
print(s4)

#multiplicadores de la derivada del error, pero no juega parte a la hora de hacer la aproximacion

#-((h1)**3)/12. 
#-((h2)**5)/90. 
#-3*((h3)**5)/80.  
#-8*((h4)**7)/945.
