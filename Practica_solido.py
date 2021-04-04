# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 11:11:50 2019

@author: adria
"""

import numpy as np
from math import *
import matplotlib.pyplot as plt

A2theta=[43.524071,50.676231,74.364532,90.314720,95.534195,117.21842,117.935143,136.964066,137.670547,145.336716,146.183456]
Aintensidad=[35.45,14.67,5.23,7.17,2.73,1.16,0.36,3.06,1.07,5.31,2.49]
A2err=[0.006963,0.010604,0.009818,0.020884,0.048211,0.044181,0.08654,0.017479,0.035951,0.023749,0.035953]

B2theta=[22.590748,32.12957,39.589741,46.017071,51.820309,57.114616,62.219395,66.992477,67.183311,71.749184,76.220078]
Bintensidad=[68.33,13.61,30.85,10.74,13.82,3.01,16.4,0.64,0.53,7.73,4.8]
B2err=[0.005056,0.009898,0.0069,0.015193,0.011812,0.016609,0.010521,0.019543,0.021384,0.01773,0.027701]
#calculo de los indices de Miller

hkl=[]
for i in range(0,6):
    for j in range(0,6):
        for k in range(0,6):            
            c=[i,j,k]
            hkl.append(c)
hkl.pop(0)      #eliminamos asi el 0,0,0       
#print(hkl)

print("--------------------------------------------------------------------------")
print("--------------------------------------------------------------------------")
#calculamos las reglas de seleccion para quedarnos con los hkl correspondientes a cada cristal

def regselbcc(l):   #nos interesan o todos pares o dos impares mas uno par
    bcc=[]
    for e in l:
        a=e[0]+e[1]+e[2]
        if a%2==0:
            bcc.append(e)
    return bcc

bcc=regselbcc(hkl)
print("HKL permitidos por una BCC")
print(" ")
print(bcc)

print("--------------------------------------------------------------------------")
print("--------------------------------------------------------------------------")

def regselfcc(l):
    fcc=[]
    for e in l:
        a=e[0]+e[1]
        b=e[0]+e[2]
        c=e[1]+e[2]
        if a%2==0 and b%2==0 and c%2==0:
            fcc.append(e)
    return fcc

fcc=regselfcc(hkl)
print("HKL permitidos por una FCC")
print(" ")
print(fcc)

print("--------------------------------------------------------------------------")
print("--------------------------------------------------------------------------")

#calculo de los cocientes experimentales y teoricos

#muestra 2
def Rexp(theta2):
    theta=[]
    for e in theta2:
        theta.append(e/2.)
        #print(Atheta)
    Rexp=[]
    for k in theta:
        a=(sin(k*pi/180.)/sin(theta[0]*pi/180.))**2
        Rexp.append(a)
    return Rexp
print("Las R experimentales para la muestra 2 son:")
RexpA=Rexp(A2theta)
for i in RexpA:
    print(i)

#es una fcc, asi que vamos a ver los criterios de fcc. Usamos num1 y num 2 como listas auxiliares
#para asi ordenar mejor los datos. 
def Rteoricocubica(l):
    num1=[]
    Rteofcc=[]
    for e in l:
        a=e[0]**2+e[1]**2+e[2]**2
        num1.append(a)
    num1=sorted(num1)
    num2 = []
    for i in num1:
        if i not in num2:
            num2.append(i)
    print(num2)         #si hacemos el print vemos la suma de los hkl. en fcc es 1+1+1, en bcc 1+1
    for k in num2:
        Rteofcc.append(k/num2[0])        
    return Rteofcc

print("Las R teoricas de la muestra 2 son:")
Rfccteo=Rteoricocubica(fcc)

for i in Rfccteo:
    print(i)
    
print("--------------------------------------------------------------------------")
print("--------------------------------------------------------------------------")
RexpB=Rexp(B2theta)
print("Las R experimentales de la muestra 4 son:")
for k in RexpB:
    print(k)
Rbccteo=Rteoricocubica(bcc)
print("Las R teoricas de la muestra 4 son:")
for k in Rbccteo:
    print(k)
    
print("--------------------------------------------------------------------------")
print("--------------------------------------------------------------------------")

#calculo del parametro D
#modificamos las listas para que queden los datos en orden
#anhadimos los terminos teoricos de los repetidos
#en la fcc se repiten el 5.33333333,6.33333333 y 6.6666666
#aplicamos sorted para hacer que las listas 
Rfccteo.append(5.333333333333333)
Rfccteo.append(6.333333333333333)
Rfccteo.append(6.666666666666667)
Rfccteo=sorted(Rfccteo)

#en la bcc se repite el 8
Rbccteo.append(8.0)
Rbccteo=sorted(Rbccteo)

#calculamos primero el factor de forma a

def dexp(R,theta):      #lista de R experimentales, angulo primero
    lamda=1.5418 #en amstrongs
    d=[]
    for i in R:
        d.append(lamda/(2*sqrt(i)*sin(theta*pi/180)))       #sale de despejar la ley de Laue en funcion de Rexp 
    return d

#calculamos las d de ambas estructuras 

dbcc=dexp(RexpB,B2theta[0]/2.)
#print(dbcc)   #si descomentamos este print y el siguiente, obtenemos los dhkl experimentales

dfcc=dexp(RexpA,A2theta[0]/2.)
#print(dfcc)

def deRD(R,d,b):  #las R, las d, el primer termino que hace que R sea 1
    sum1=0
    sum2=0
    for i in range(len(d)):
        sum1+=(b*R[i])**2
        sum2+=b*R[i]/(d[i]**2)
    return sqrt(sum1/sum2)

Abcc=deRD(Rbccteo,dbcc,2)
print("La constante de red de la bcc es:")
print(Abcc)

Afcc=deRD(Rfccteo,dfcc,3)
print("La de la fcc es:")
print(Afcc)

print("--------------------------------------------------------------------------")
print("--------------------------------------------------------------------------")

#como d teorico es a/sqrt(h**2+k**2+l**2), procedemos a calcularlo

def dteo(R,a,b):
    d=[]
    for i in R:
        d.append(a/(sqrt(i*b)))
    return d

dteobcc=dteo(Rbccteo,Abcc,2)
#print(dteobcc)        #si descomentamos tanto este print como el siguiente, obtenemos los valores teoricos de dhkl
dteofcc=dteo(Rfccteo,Afcc,3)
#print(dteofcc)

#calculamos la D

def D(dteo,dexp):
    resu=0
    for i in range(len(dexp)):
        resu+=(1/(dteo[i])**2-1/(dexp[i])**2)**2
    return resu

Dbcc=D(dteobcc,dbcc)
print("La D de la bcc es:")
print(Dbcc)

Dfcc=D(dteofcc,dfcc)
print("La D de la fcc es:")
print(Dfcc)

print("--------------------------------------------------------------------------")
print("--------------------------------------------------------------------------")
#PROPAGACION DE ERRORES

#primero de la R. Aplicamos logaritmo diferencial para propagar el error en los valores experimentales

def errRexp(theta,delta,Rexp):
    err=[]
    for j in range(len(theta)):
        theta[j]=theta[j]/2.
        delta[j]=delta[j]/2.        #el angulo es la mitad, el error tambien
    for i in range(len(Rexp)):
        err.append(2*Rexp[i]*((delta[i]*pi/180.)/tan(theta[j]*pi/180.)+(delta[0]*pi/180.)/tan(theta[0]*pi/180.)))
    return err

deltaRbcc=errRexp(B2theta,B2err,RexpB)
print("el error en la medida de Rexp de la bcc es:")
print(deltaRbcc)

deltaRfcc=errRexp(A2theta,A2err,RexpA)
print("El error en la medida de Rexp de la fcc es:")
print(deltaRfcc)
print("--------------------------------------------------------------------------")
print("--------------------------------------------------------------------------")

#ahora el de d

def errdexp(theta,delta,dexp):
    err=[]
    for j in range(len(theta)):
        theta[j]=theta[j]/2.
        delta[j]=delta[j]/2.
    for i in range(len(theta)):
        err.append(dexp[i]*(delta[i]*pi/180.)/tan(theta[i]*pi/180.))
    return err

errdbcc=errdexp(B2theta,B2err,dbcc)
print("el error en la d experimental de la bcc es:")
print(errdbcc)
    
errdfcc=errdexp(A2theta,A2err,dfcc)
print("el error en la d experimental de la fcc es:")
print(errdfcc)

print("--------------------------------------------------------------------------")
print("--------------------------------------------------------------------------")
#error en a
def erra(Rteo,dexp,errd,a,b):
    suma1=0
    suma2=0
    for i in range(len(dexp)):
        suma1+=Rteo[i]*b*errd[i]/(dexp[i])**3
        suma2+=Rteo[i]*b/(dexp[i])**2
    return a*suma1/suma2

errabcc=erra(Rbccteo,dbcc,errdbcc,Abcc,2)
print("el error en a de la bcc es:")
print(errabcc)

errafcc=erra(Rfccteo,dfcc,errdfcc,Afcc,3)
print("el error en a de la fcc es:")
print(errafcc)