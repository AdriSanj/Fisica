# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
"""

#derivacion numerica

from sympy import Symbol , diff , sin, sqrt		#hace derivadas de funciones, no numericas
x=Symbol('x')		#asi x sera x como tal,simbolica,no un numero
y=Symbol('y')
#alternativa: from sympy.abc import x , y
#definimos funciones simbolicas
f=sqrt(x**2+y**2)
g=sin(x*y)
#derivamos
dfx=diff(f,x)
dgx=diff(g,x,1,y,2)
dfx
dgx
#evaluacion en un punto

valorf=dfx.subs(x,1)    #sustituye la x en el valor 1
valorf
valorf=dfx.subs([(x,1),(y,2)])
valorf
valorg=dgx.subs({x:1,y:2})
valorg

del x,y



#para las derivadas numericas, el modulo a usar es scipy.misc.derivative

from scipy.misc import derivative    #importamos derivative, generamos una funcion, y hacemos su derivada en dicho punto
funcion=lambda z: z**2+1/z
derivative(funcion, 1.0, dx=1e-6)
#derivative(func, x0, dx=1.0, n=1, args=(), order=3)  funcion, punto, diferencial, n, order numero de puntos
#si no decimos nada, toma una formula de derivacion con 3 puntos


#definicion de la formula numerica 
def deriva1(funcion,x,h=0.001):
    return(funcion(x+h)-funcion(x))/h
    
#definimos simbolicamente para comparar resultados
from sympy import *
from sympy.abc import x
f_sym=sin(pi*x)
df_sym=diff(f_sym,x)
f_num=lambdify(x,f_sym,'numpy')
df_num=lambdify(x,df_sym,'numpy')
#las funciones numericas son para ser evaluar, lamdify la convierte en numerica
#asi definimos funciones numericas a traves de simbolicas, podemos comparar objetos asi

#usamos un punto para comparar
pto=0
print(deriva1(f_num,pto))
print(df_num(pto))



#verificamos en varios puntos

import numpy as np
pts=np.arange(5)/4.0
print(deriva1(f_num,pts))
print(df_num(pts))
print(df_num(pts)-deriva1(f_num,pts))           #diferencia entre valores






#ejemplo, representar error cometido con h=0.1 en el intervalo de extremos cero y pi
#EJERCICO 1
pts=np.linspace(0,np.pi)
import matplotlib.pyplot as plt
fig=plt.figure(1)
plt.subplot(211)
plt.plot(pts,df_num(pts),pts,deriva1(f_num,pts,0.1))
plt.title(r'Derivación numérica (2 pts)', fontsize=20)
plt.ylabel(r'$f\left(x\right)$', fontsize=20)
plt.legend(('Exacta','Numérica'))
plt.grid()
plt.subplot(212)
plt.plot(pts,df_num(pts)-deriva1(f_num,pts,0.1))
plt.xlabel(r'$x$', fontsize=20)
plt.ylabel(r'$e\left(x\right)$', fontsize=20)
plt.grid()
plt.show()











