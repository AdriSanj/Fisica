# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 16:05:59 2018

@author: UO260284
"""

#ejemplo 1
from sympy import *
from sympy.abc import x, a
y = Function('y')           #definimos el objeto funcion simbolica
ecu= Derivative(y(x),x)+a*y(x)      #generamos la ecuacion y la iguala a cero
s=dsolve(ecu,y(x))      #resolvemos y de x
Y=s.rhs        #el lado derecho de la ecuacion lo llamamos Y
checkodesol(edo,Y)      #comprobamos que es la solucion
print(Y)

%reset -f



#ejemplo 2
import sympy as sym 
x, a=sym.symbols('x a')         #importamos simbolos
f=sym.symbols('f',cls=sym.Function)             #define f como una funcion simbolica
edo=sym.Eq(f(x).diff(x, x) - 2*f(x).diff(x) + f(x), sym.sin(x))  #definimos la edo. la parte de izquierda es una parte de la ecuacion, la de la derecha es la parte derecha de la ecu izq=dere 
s=sym.dsolve(edo,f(x))      #resolvemos. Eq es un tipo de objeto "ecuacion". Contiene info de la ecuacion
Y=s.rhs         #definimos la solucion
sym.checkodesol(edo,Y)      #comprueba
print(Y)



#ejemplo 3
#sympy no esta programado para resolver problemas de valor inicial. solo resuelve ecuaciones diferenciales , no problemas
# si queremos resolver un problema de valor inicial, hasta la s en este ejemplo lo hacemos igual, pero luego resolvemos
#las condiciones que yo pongo para obtener los valores de dichas constantes 
import sympy as sym 
from sympy.abc import x, a 
y = sym.Function('y') 
edo = sym.Derivative(y(x),x)+a*y(x) 
s=sym.dsolve(edo,y(x))          # hasta aqui la solucion general 
Yg=s.rhs        #Yg=C1*exp(-a*x) , sol general de la edo
ctes=sym.solve(Yg.subs(x,0)-1,dict=True) #tenemos la condicion inicial y(0)=1
#sustituimos el valor 0 en x y tiene que satisfacer y(0)-1=0 segun el comando anterior. dict = true  genera el diccionario (es un array) tal que {c1:1} asigna los valores con las variables que tengo libres
#entonces Yp=exp(-a*x)
#si no pones dict = true, pone [1].
#en SEDOS genera un array de diccionarios.  hay que tomar despues las soluciones por separadas
Yp=Yg.subs(ctes[0])         #sustituye el primer elemento
print(Yp)
#matlab deja que esten las condiciones de contorno de la edo. python no



#ejemplo 4
#EDO DE GRADO 2
# las soluciones no suelen ser unicas. pero esta funcion solo devuelve una solucion
import sympy as sym 
from sympy.abc import x 
y = sym.Function('y') 
edo = sym.Eq((y(x).diff(x))**2+(y(x))**2,1) 
s=sym.dsolve(edo,y(x)) 
Yg=s.rhs
print(Yg)



#ejemplo 5
#SEDOS. haremos una lista de objetos ecuaciones con Eq
import sympy as sym 
from sympy.abc import t 
x,y=sym.symbols('x y',cls=sym.Function) 
edo = (sym.Eq((x(t).diff(t))-y(t)),sym.Eq((y(t).diff(t))+x(t))) 
s=sym.dsolve(edo) 
Xg=s[0].rhs 
Yg=s[1].rhs
print(Xg)
print(Yg)
#se pueden imponer condiciones iniciales a estos sistemas

#EJERCICIO 1
#y'=1+y/t  en [1,2] y ademas y(1)=2

import sympy as sym 
x=sym.symbols('x')         #importamos simbolos
y=sym.symbols('y',cls=sym.Function)             #define f como una funcion simbolica
edo=sym.Eq(y(x).diff(x) , 1+y(x)/x)     #definimos la edo
s=sym.dsolve(edo,y(x))      #resolvemos
Yg=s.rhs
print(Yg)
ctes=sym.solve(Yg.subs(x,1)-2,dict=True)        #sabemos que y(1)=-2. entonces y(1)-2=0
Yp=Yg.subs(ctes[0])  #sustituimos
print(Yp)           #booooom

#ejercicio 2
#y'=1+(x-y)**2  en [2,3]   e y(2)=1
import sympy as sym 
x=sym.symbols('x')         #importamos simbolos
y=sym.symbols('y',cls=sym.Function)             #define f como una funcion simbolica
edo=sym.Eq(y(x).diff(x) , 1+(x-y(x))**2)     #definimos la edo
s=sym.dsolve(edo,y(x))      #resolvemos
Yg=s.rhs
sym.checkodesol(edo,Yg) #nos da una solucion aproximada, no es capaz de resolverlo 
print(Yg)
ctes=sym.solve(Yg.subs(x,2)-1,dict=True)        #sabemos que y(1)=-2. entonces y(1)-2=0
Yp=Yg.subs(ctes[0])  #sustituimos
print(Yp)   
