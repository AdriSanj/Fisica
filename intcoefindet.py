# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 12:51:16 2018

@author: UO260284
"""

import numpy as np
import sympy as sp

x=sp.symbols('x')
def matriz(s):          #funcion que hace la matriz de los valores 
    M=np.ones((len(s),len(s)))
    for i in range(len(s)):
        M[i,]=s**i      #cada fila de la matriz es un x elevado al numero de la fila que sea
    return M

def soluciones(r,a,b):      #nos da las soluciones a la integral normal. a y b es el intervalo de integracion
    N=[]
    for i in range(0,len(r)):         #anhadimos a la lista cada valor que tendra la integral dependiendo del grado de la x(el grado maximo es la longitud del vector del soporte)
        z=sp.integrate(x**i,(x,a,b))
        N.append(z)
    N=np.asarray(N)
    return N


soporte=np.array([0.25,0.5,0.75])           #soporte. podemos poner cualquier numero
#print(soporte)

K=matriz(soporte)       #calculamos la matriz 
#print(K)
                    #invertimos 
K=np.linalg.inv(K)
#print(K)

                #vector de soluciones
f=soluciones(soporte,0,1)
#print(f)

                #multiplicamos la inversa de K por f
C=K.dot(f)      
print(C)             #estos son los coeficientes

#orden

        
