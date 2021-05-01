#derivacion por python

import numpy as np
import matplotlib.pyplot as plt
from math import *
from scipy.misc import derivative

#derivacion por dos puntos
def derivacion5(f,a,h):
    d1=-1/(12*h)*(-25*f(a)+48*f(a-h)-36*f(a-2*h)+16*f(a-3*h)-3*f(a-4*h)) #lateral atrasada
    d2=1/(12*h)*(-3*f(a-h)-10*f(a)+18*f(a+h)-6*f(a+2*h)+f(a+3*h)) #semilateral adelantada
    d3=1/(12*h)*(f(a-2*h)-8*f(a-h)+8*f(a+h)-f(a+2*h))           #centrada
    return d1,d2,d3

def derivacion3(f,a,h):
    d1=1/(2*h)*(-3*f(a)+4*f(a+h)-f(a+2*h))
    d2=1/(2*h)*(-f(a-h)+f(a+h))
    d3=-1/(2*h)*(-3*f(a)+4*f(a-h)-f(a-2*h))
    return d1,d2,d3


def f(x):
    return sqrt(log(x**2+1))

D1,D2,D3=derivacion5(f,0.46,0.2019)
print D1
print D2
print D3

print(derivative(f,0.46,dx=1e-6)-D1)
