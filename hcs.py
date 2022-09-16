#instability_dens.py
# coding=utf-8
# from matplotlib.lines import _LineStyle
import numpy as np
import pandas as pd
import scipy.special as special
from numba import jit
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True
from matplotlib.transforms import (
    Bbox, TransformedBbox, blended_transform_factory)
import random
# import scipy.special as special

# Import math Library
import math 

from sympy import *

# PARA USAR A LA HORA DE GUARDAR DOS COLUMNAS EN UN ARCHIVO
import csv

import cmath

import functions_d as f

from sklearn.linear_model import LinearRegression

from sklearn.metrics import mean_absolute_error


# !! Definir parametros 
alfa = 0.8
d= 3


z = f.Panel(d).zeta(alfa)
eta = f.Panel(d).eta(alfa)
c = f.Panel(d).kapa(alfa)
m =  f.Panel(d).mu(alfa)
# print("zeta")
# print(z)
# print("kapa")
# print(kapa)
# c= Symbol("kapa",positive=True)
# k = Symbol("k")
# m=Symbol("mu",positive=True)
# eta= Symbol("eta")
# z= Symbol("zeta")
# d = Symbol("d")
k= Symbol("k")
x = Symbol('x')


poly_d  =  x**3 + ( ( d + 2 ) * c / ( 2 * ( d -1 ) ) 
         
         + ( d - 1 ) * eta / d) *  k**2 * x**2 + ( ( d + 2 ) * k**4 * z * c / (2 * d ) 
         
         + ( k**2 / d )* ( d + 2 + ( d - 1 ) * z * eta - d * ( d + 2 ) * z * c / ( 2 * ( d - 1 ) ) )
         
         - z**2 ) * x + ( ( d + 2 ) * ( c - m )* k**2 / (2 * ( d - 1 ) ) - z ) * k**2 

print(poly_d)


quadratic_equation = Eq(poly_d, 0)

solucion=solve(quadratic_equation, x)
 
print(solucion[0])

X = np.linspace(0, 0.6,1000)
# Won't work -- f is not a function.
# Y = f(X)
# Works.
stangent = z - 0.5*eta*k**2  
Y = lambdify(k, (solucion[0]), "numpy")(X)
Z = lambdify(k, (solucion[1]), "numpy")(X)
W = lambdify(k, (solucion[2]), "numpy")(X)
U = lambdify(k,stangent , "numpy")(X)


fig = plt.figure()

plt.plot(X, (Y),color='C0',label=r'$\lambda_1$ ')
plt.plot(X, (Z),color='C1',label=r'$\lambda_2$ ')
plt.plot(X, (W),color='C2',label=r'$\lambda_3$ ')
plt.plot(X, (U),color='C3',label=r'$\lambda_3$ ')
plt.legend(loc=0,fontsize=30)
plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)

plt.xlabel( r' $k$ ', fontsize=30)




# ! plot hydrodynamics coefficients
alpha = np.linspace(0.0,0.88,1000)
fig2 = plt.figure()
plt.plot(alpha, f.Panel(d).kapa(alpha),color='C0',label=r'$\kappa(\alpha), \, d=1$ ')
plt.plot(alpha, f.Panel(d).mu(alpha),color='C1',label=r'$\mu(\alpha),\, d=1$ ')

plt.legend(loc=0,fontsize=30)
plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)

plt.xlabel( r' $\alpha$ ', fontsize=30)

plt.show()