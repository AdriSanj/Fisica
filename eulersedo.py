#resolucion de SEDOS  por el metodo de euler 

import numpy as np
import matplotlib.pyplot as plt

def eulervect(F,ab,n,ya):           #funcion, intervalo, numero de puntos, valores en el momento inicial
    a=ab[0]         #punto inicial y final
    b=ab[1]
    h=(b-a)/float(n)            #paso
    ti=a
    yi=ya.reshape((len(ya),1))          #la transponemos, queremos que quede en columnas
    T=[ti]        #variable t
    resu=yi
    for i in range(0,n):
        ti1=ti+h       #o valor de t sera o da suma do seu valor anterior e mailo paso
        yi1=yi+h*F(ti,yi)
        T.append(ti1)
        resu=np.concatenate((resu,yi1),axis=1)          #concatenamos por columna hasta obtener la matriz de resultados, tantas columnas como puntos queramos, tantas filas como funciones con la variable dependiente t
        ti=ti1
        yi=yi1
    return T,resu

f=lambda t,Y:np.array([Y[1],2./(t**3+3.5)-9*(((Y[0])**3)+Y[0]*Y[1])])

ab=np.array([1,2])
ya=np.array([2/9.,8/27.])


T,resu=eulervect(f,ab,14,ya)
X=resu[0,:]         #la primera fila de la matriz son los resultados en X(t), la segunda los de Y(t)
Y=resu[1,:]
print(T)
print len(T)
print('funcion')
print(len(X))
print(X)
print('derivada')
print(Y)
plt.figure()
plt.plot(T,X,'r',T,Y,'b')
plt.grid()
plt.show()
