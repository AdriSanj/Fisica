# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 12:23:26 2018

@author: UO260284
"""
import sympy as sp
n,x,y=sp.symbols('n x y')
I=sp.integrate(sp.cos(n*x),x)
Z=sp.integrate(sp.cos(x**2),(x,0,sp.pi/2))
print(I)
print(Z)
W=sp.integrate(sp.exp(-x),(x,0,sp.oo))          #infinito sp.oo
print(W)
#integral doble
Doble=sp.integrate(sp.exp(-x**2-y**2),(x,-sp.oo,sp.oo),(y,0,sp.oo))
print(Doble)
P=sp.integrate(x**x,x)
print(P)    #no hay primitiva 

K=sp.Integral(sp.log(x)**2,x)           #integral diferida, no la calcula, es la representacion 
K.doit()            #esto calcula la primitiva
print(K)
#entonces Integral es la representacion simbolica, con el doit para resolver simbolicamente
#para resolver la definida es integral
#lo usamos para hacer comparaciones (el Integral)

