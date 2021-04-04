# -*- coding: utf-8 -*-


from numpy import *
import matplotlib.pyplot as plt
from numpy.random import rand,randint
    
def choque(x,v,l):
    for i in range(len(x)-1):       #analizamos el array de velocidades consigo mismo, tomamos las componentes x e y de cada particula y las comparamos con la de cualquier otra
        for j in range(len(x)-1):   #en caso de que su distancia en modulo sea menor que un epsilon (en este caso, menor que el doble del radio) y no sea consigo misma, calcula la nueva velocidad debido al choque
            M=sqrt((x[i][0]-x[j][0])**2+(x[i][1]-x[j][1])**2)
            if M!=0 and M<l*1e-2:           #si no esta evaluando la particula a si misma y la distancia entre ellas es menor que dos veces el radio, se ejecuta el if
                M=sqrt((x[i][0]-x[j][0])**2+(x[i][1]-x[j][1])**2)
                r12=x[i][:]-x[j][:]     #formuals derivadas de la mecanica clasica y conservacion de momentos y energia cinetica
                uparalelo=r12/M     #vector unitario del choque
                uperpendicular=array([-abs(uparalelo[1]),uparalelo[0]])
                V1par=dot(v[i][:],uparalelo)
                V1per=dot(v[i][:],uperpendicular)
                V2par=dot(v[j][:],uparalelo)
                V2per=dot(v[j][:],uperpendicular)
                v[i][0]=V2par
                v[i][1]=V1per
                v[j][0]=V1par  
                v[j][1]=V2per                                               
    return v                  
N=50       #numero de particulas
Kb=0.01     #constante de boltzmann
l=15   #lado de la caja
T=100.           #temperatura
eps=0.25
ptosx=array([0,0,l,l,0])
ptosy=array([0,l,l,0,0])
E=random.exponential(Kb*T,size=(N,2))       #funcion exponencial que nos genera energias aleatorias 
v=sqrt(2*E)  #velocidad sacada de la anterior energia 
L=15           #longitud de la caja
Lx=L   #la Lx corresponde a la longitud del eje x, por eso que las particulas del eje y suman la presion sobre ella. Esta arriba y abajo
Ly=L   #lo mismo para la Ly, son las del eje x las que chocan 
t=0
dt=0.1
x=abs(randint(1,l,(N,2)))      #posiciones iniciales
#v=2*l*rand(N,2)
Tmax=100          #Tmax= 1000 y dt=0.1 asegura 10000 iteraciones como nos pide el ejercicio
plt.figure('part',figsize=(10,8))
Presiones=[]        #hacemos una lista que va a estar vacia cada diez intervalos. con ella sumamos la presion para obtener la media en ese rango de tiempo 
Ppause=0            #lo usamos para definir el momento donde lo medimos
Presionesmedidas=[]
Tiempopresiones=[]
while t<Tmax:
    v=choque(x,v,l)     #comprobamos si choca alguna 
    ixl=where(x[::,0]>(l-eps))[0]           #comprobamos si se salen de la caja
    iyl=where(x[::,1]>(l-eps))[0]           #para aquellas componentes que se salga, primero hacemos el cla para que se queden en su sitio
    ix0=where(x[::,0]<eps)[0]               #luego, esa componente de la velocidad es cambiada de direccion
    iy0=where(x[::,1]<eps)[0]
    v[ix0,0]=abs(v[ix0,0])
    npartx0=len(v[ix0,0])
    Fx0=2*abs(v[ix0,0])/dt              #calculamos la fuerza de cada particula en una rray apra cada pared
    Px0=sum(Fx0)/Ly                     #sumamos todas las fuerzas y obtenemos lapresion si la dividimos entre la longitud de la pared
    v[ixl,0]=-abs(v[ixl,0])
    npartxl=len(v[ixl,0])
    FxL=2*abs(v[ixl,0])/dt
    PxL=sum(FxL)/Ly
    v[iy0,1]=abs(v[iy0,1])
    nparty0=len(v[iy0,1])
    Fy0=2*abs(v[iy0,1])/dt
    Py0=sum(Fy0)/Lx
    v[iyl,1]=-abs(v[iyl,1])
    npartyl=len(v[iyl,1])
    FyL=2*abs(v[iyl,1])/dt
    PyL=sum(FyL)/Lx
    x=x+v*dt
    t=t+dt
    #print('-')
    plt.subplot(211)        #hacemos  un subplot  y representamos
    plt.cla()
    plt.plot(ptosx,ptosy,'b')
    plt.plot(x[::,0],x[::,1],'r.')
    plt.draw()
    #print('Numero total de particulas que chocan contra las paredes en este instante')
    #print(npartx0+npartxl+nparty0+npartyl)
    Presiones.append(Px0+PxL+Py0+PyL)
    Ppause+=1
    print(' ')
    if Ppause%10==0:
        #print('Presion media')
        total=0
        for i in Presiones:
            total+=i
        #print(total/10)
        Presionesmedidas.append(total/10)
        Tiempopresiones.append(t)
        Presiones.clear()
    E=ones(len(x))      #histograma de la energia, en cada momento. Hacemos el producto escalar del vector velocidad de cada particula 
    for k in range(len(E)):
        E[k]=dot(v[k][:],v[k][:])
    E=0.5*E
    #print(E)
    plt.subplot(212)
    plt.cla()
    plt.title('Distribucion de energia')
    plt.ylabel('Numero de particulas')
    plt.xlabel('Energias')
    plt.hist(E,bins=20)
    plt.pause(0.001)

#representamos la presion en intervalos de tiempo al acabar de ejecutar el programa
plt.figure()
plt.title('Presion vs tiempo')
plt.plot(Tiempopresiones,Presionesmedidas)
plt.xlabel('Tiempo')
plt.ylabel('Presion')
plt.grid()
plt.show()

print(mean(asarray(Presionesmedidas)))


#PARA EL CASO DE DOS TIPOS DE PARTICULAS CON DIFERENTE MASA
#v1'paralelo=(v1paralelo*(m1-m2)+2*v2paralelo*m2)/(m1+m2)
#v2'paralelo=(v2paralelo*(m2-m1)+2*v1paralelo*m1)/(m1+m2)
#N1=50
#N2=50
#T1= 0 (temperaturas)
#T2=50
#m1=1
#m2=50
