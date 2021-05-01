#trabajo galaxias
from numpy import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from astropy.table import Table
from visual import *
import Image

#PARTE 1, DETERMINACION DE LA CONSTANTE FRENTE A SUS VALORES
def hubble(r,K):
    '''
    Funcion de la obtencion de las velocidades de alejamiento
    o acercamiento del objeto respecto los valores calculados
    por Hubble. La constante es un parametro de entrada que se
    calcula por la anterior relacion lineal
    '''
    return K*r
    
def f (x, a): return a*x			#funcion lineal de recurso para aproximar la constante de hubble
	
def incertidumbre(b):
	''' funcion para el calculo
		de la desviacion sigma
	'''
	N=float(len(b))
	media=sum(b)/N
	desv=abs(b-media)
	d=sqrt((sum(desv**2))/N)
	return desv,d
	
def residuo(v,u,desv):
	'''
	funcion que calcula chi cuadrado
	'''	
	N=float(len(u))		# u corresponde a los valores de velocidad observados, v a los del ajuste
	m=sum(u)/N			#hallamos su media
	chi2=sum((abs(u-v)**2)/(m-u)**2)			#resultado de chi
	return chi2

datos=loadtxt('datoshubble.txt',dtype='float',comments='#',delimiter=',')		#cargamos datos antiguos
obj=datos[:,0]         			#nombres de los objetos celestes
r=datos[:,1]		#leemos distancias en parsecs *10**6 (megaparsec) y las velocidades en km/s
v=datos[:,2]

#obtencion de la constante con curvefit y su error
desv,d=incertidumbre(v)			#valor necesario para calcular la pendiente por curvefit
hcurvefit,cov=curve_fit(f, r, v, sigma=desv) 
errh=sqrt(diag(cov))				#la bondad del ajuste de curvefit

#obtencion mediante linealsqrt
desv.shape=(21,1)
r.shape=(21,1)				#cambiamos formatos para cuadrar con leastsquare
v.shape=(21,1)
hlstsq, res, rango, singa = linalg.lstsq(r,v)       #ajuste por minimos

#tablas de resultados
resul11=hubble(r,hlstsq)  #linsqrt
resul12=hubble(r,hcurvefit[0])			#resultados con curve_fit
residuos1=residuo(resul11,v,desv)			#bondad del ajuste de leastsquare mediante una funcion que nos devuelve chi cuadrado
Data=Table([obj,r,v,desv,resul11,resul12],names=('N.G.C.','distance (Mpc)','velocities(km/s)','desv resp mean','estimated vel(km/s) lins','estimated vel(km/s) curve'))		#funcion de astropy para generar tablas de datos, el numero 1 nos indica en que dimension esta el array
print Data
print ('Old Hubble constant value with least_square= %4.3f ; residual= %4.3f '%(hlstsq[0],residuos1))
print ('Old Hubble constant value with curve_fit= %4.3f ; residual= %4.3f '%(hcurvefit[0],errh[0]))

#- - - - - - - - - - - - - - - - - - - - -
p1=raw_input('Press enter to see the taken 90s values...')			#continuacion del programa
#- - - - - - - - - - - - - - - - - - - - -
datos2=loadtxt('hubbledatos1991.csv',dtype='float',comments='#',delimiter=';')			#volvemos a cargar datos, los nuevos
obj2=datos2[:,0]
r2=datos2[:,2]
v2=datos2[:,1]

#obtencion mediante curve_fit
desv2,d2=incertidumbre(v2)
hcurvefit2,cov2=curve_fit(f, r2, v2, sigma=desv2)
errh2=sqrt(diag(cov2))			#bondad del segundo ajuste por curvefit

#mediante linsquare
r2.shape=(24,1)			#convertimos para que el formato se adapte al lstsq
v2.shape=(24,1)
hlstsq2, res2, rango2, singa2 = linalg.lstsq(r2,v2)			#calculo con least square

#segunda tabla de resultados
resul21=hubble(r2,hlstsq2)			#constante con lsrsq	
resul22=hubble(r2,hcurvefit2[0])		#constante con curve_fit
residuos2=residuo(resul21,v2,desv2)		#bondad del segundo ajuste por leastsquare
Data2=Table([obj2,r2,v2,desv2,resul21,resul22],names=('N.G.C.','distance (Mp)','velocities(km/s)','desv resp mean','estimated vel(km/s) lins','estimated vel(km/s) curve'))
print Data2
print ('Actual Hubble constant value with least_square= %4.3f ; residual= %4.3f '%(hlstsq2[0],residuos2))
print ('Actual Hubble constant value with curve_fit= %4.3f ; residual= %4.3f '%(hcurvefit2[0],errh2[0]))

#- - - - - - - - - - - - - - - - - - - - -
p2=raw_input('Press enter to see de graphic with the comparison of values...')
#- - - - - - - - - - - - - - - - - - - - -

dandr=0.775		#distancia de andromeda, en megaparsecs
andromeda1=hubble(dandr,hlstsq)			#valor de la velocidad de andromeda segun el ajuste dado
andromeda2=hubble(dandr,hcurvefit[0])
andromeda3=hubble(dandr,hlstsq2)
amdromeda4=hubble(dandr,hcurvefit2[0])

#comparacion entre resultados obtenidos

#datos de hubble
plt.subplot(2,1,1)
plt.plot(r,v,'ro',label='Hubble Data')
plt.plot(0.775,andromeda1,'m*',label='Andromeda Lstsq')         #ejemplo concreto de donde se situaria mediante un ajuste u otro,con su distancia
plt.plot(0.775,andromeda2,'b*',label='Andromeda CurveFit')
plt.plot(r,resul11,'m-',label='Lstsq results')
plt.plot(r, resul12, 'b-', label='Curvefit results')
plt.grid()
plt.legend(loc=4)				#localizacion de la leyenda
plt.title('Distances vs Velocities')
plt.ylabel('V (Km/s)')
plt.xlabel('r (Megaparsecs)')

#datos de los 90
plt.subplot(2,1,2)
plt.plot(r2,v2,'ro',label='1991 Data')
plt.plot(0.775,andromeda1,'m*',label='Andromeda Lstsq')        
plt.plot(0.775,andromeda2,'b*',label='Andromeda CurveFit')
plt.plot(r2,resul21,'m-',label='Lstsq results')
plt.plot(r2, resul22, 'b-', label='Curvefit results')
plt.grid()
plt.legend(loc=4)			#localizacion de la leyenda
plt.title('Distances vs Velocities')
plt.ylabel('V (Km/s)')
plt.xlabel('r (Megaparsecs)')
plt.show(block=False)

#- - - - - - - - - - - - - - - - - - - - -

#- - - - - - - - - - - - - - - - - - - - -
p3=raw_input('Press enter to watch a visual representation of two galaxies getting away with the velocities determined with the Hubble constant...')
#- - - - - - - - - - - - - - - - - - - - -

im = Image.open('galaxy1.jpg')                                  #anadimos imagen a las galaxias  
tex1 = materials.texture(data=im, mapping='rectangular')
im = Image.open('galaxy2.jpg')  
tex2 = materials.texture(data=im, mapping='rectangular')

#tierra
tierra=sphere(pos=(0,0,0),radius=0.4,material=materials.earth)              #sistema de referencia la tierra

#galaxia 1
g1=ellipsoid(pos=(3.13,0,0),length=2,height=1,width=1,color=color.white, material=tex2)
g1.rotate(angle=pi)

#galaxia2
g2=ellipsoid(pos=(-14.79/5.,0,0),length=2,height=1,width=1,color=color.white, material=tex1) #distancias en megaparsecs
g2.rotate(angle=pi)

g1.vel=vector(343,0,0) #km/s                    #vectores de veocidades iniciales de los datos
g2.vel=vector(-1350,0,0) #km/s 
dt=0.00003

label1 = label(pos=g1.pos, xoffset=20, yoffset=20, height=10, border=6, font='sans')            #etiquetas para las galaxias
label2 = label(pos=g2.pos, xoffset=20, yoffset=20, height=10, border=6, font='sans')

while 1:                #bucle indefinido
	rate(50)
	g1.vel.x=g1.pos.x*72.362
	g2.vel.x=g2.pos.x*72.362
	g1.pos=g1.pos+g1.vel*dt
	g2.pos=g2.pos+g2.vel*dt
	if abs(g1.pos.x)>15:               #cuando superen dicha distancia, veremos el efecto de corrimiento al rojo
		g1.color=color.red
	if abs(g2.pos.x)>15:	
		g2.color=color.red	
	label1.pos=g1.pos
	label1.text='v=%2.6f \n x=%2.6f'%(g1.vel.x,g1.pos.x)
	label2.pos=g2.pos
	label2.text='v=%2.6f \n x=%2.6f'%(g2.vel.x,g2.pos.x)

p4=raw_input('Press enter to end the program...')

