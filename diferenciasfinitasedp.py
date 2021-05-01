# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 12:09:27 2018

@author: UO260284
"""

#EDP
import numpy as np
import matplotlib.pyplot as plt
from math import *

def mdf2elip(X,Y,edpfun,ccfun):   #numero de interespaciados. Esto lo cambiamos como queramos
    nx=X[2]
    ny=Y[2]
    Xt=np.linspace(X[0],X[1],nx+1)
    Yt=np.linspace(Y[0],Y[1],ny+1)
    m=len(Xt)*len(Yt)   
    sol=np.zeros((m,m))            #la matriz sera del tamaño MxN,MxN
    F=[]      #lista donde anhadiremos los valores de la funcion f vectorial
    for i in range(0,len(Xt)):
         for j in range(0,len(Yt)):
             F.append(edpfun(Xt[i],Yt[j])[1])   
    F=np.asarray(F)     convertimos la lista en array
    return Xt,Yt,F
    
    
    
    
def edpfun(x,y):    #funciones que aparecen en el programa
    return np.array([0,-4])     #alfa es la primera, beta la segunda

def ccfun(x,y):         #condiciones de dirichlet
    return np.array([np.pow(x,2),pow(y,2),pow(x-2,2),pow(y-1,2)]) #en orden: u(x,0);u(0,y);u(x,2);u(1,y)

X=np.array([0,1,4])
Y=np.array([0,2,5])

a,b,F=mdf2elip(X,Y,edpfun,ccfun)
print(F)
print(b)
