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
alfa = 0.98
d = 2


ze = f.Panel(d).zet(alfa)
ze1d = f.Panel(1).zet(alfa)

nu1= f.Panel(d).nu1(alfa)

eta = f.Panel(d).eta(alfa)
eta1d = f.Panel(1).eta(alfa)

kapa = f.Panel(d).kapa(alfa)
kapa1d = f.Panel(1).kapa1d(alfa)

mu =  f.Panel(d).mu(alfa)
mu1d =  f.Panel(1).mu1d(alfa)

a2= f.Panel(d).a2d(alfa)


# print("kapa")
# print(kapa)
print("a2")
print(a2)

# c= Symbol("kapa",positive=True)
# k = Symbol("k")
# m=Symbol("mu",positive=True)
# eta= Symbol("eta")
# z= Symbol("zeta")
# d = Symbol("d")
k= Symbol("k")
x = Symbol('x')
kc = np.sqrt(2*(d-1)*ze/((d+2)*(kapa-mu)))
# ! Polinomio característico de las ecuaciones hydrodinamicas en d dimensiones
poly_d  =  x**3 + ( ( d + 2 ) * kapa / ( 2 * ( d -1 ) ) 
         
         + ( d - 1 ) * eta / d) *  k**2 * x**2 + ( ( d + 2 ) * k**4 * ze * kapa / (2 * d ) 
         
         + ( k**2 / d )* ( d + 2 + ( d - 1 ) * ze * eta - d * ( d + 2 ) * ze * kapa / ( 2 * ( d - 1 ) ) )
         
         - ze**2 ) * x + ( ( d + 2 ) * ( kapa - mu )* k**2 / (2 * ( d - 1 ) ) - ze ) * k**2 
         
# print(poly_d)

# ! Establecer el polinomio igualado a 0 para resolverlo posteriormente 
# quadratic_equation = Eq(poly_d, 0)
# ! Resolver polinomio en función de x 
# solucion=solve(quadratic_equation, x)

# ! printear las raices del polinomio 
# print("solucion 1")
# print(solucion[0])
# print("solucion 2")
# print(solucion[1])
# print("solucion 3")
# print(solucion[2])
# ! Raiz para el modo transversal 
stangent = ze - 0.5*eta*k**2  

poly_1d  =  x**3 + ( ( 3 ) * kapa1d / ( 2  ) ) *  k**2 * x**2 + ( ( 3 ) * k**4 * ze1d * kapa1d / (2 ) 
         
         + ( k**2  )* ( 3  - 3  * ze1d * kapa1d / ( 2  ) )
         
         - ze1d**2 ) * x + ( ( 3) * ( kapa1d - mu1d )* k**2 / (2  ) - ze1d ) * k**2 
print(poly_1d)

# ! Establecer el polinomio igualado a 0 para resolverlo posteriormente 
# quadratic_equation = Eq(poly_1d, 0)
# ! Resolver polinomio en función de x 
# solucion1d=solve(quadratic_equation, x)
# ! printear las raices del polinomio 
# print("solucion 1")
# print(solucion1d[0])
# print("solucion 2")
# print(solucion1d[1])
# print("solucion 3")
# print(solucion1d[2])


X = np.linspace(0, 0.6,1000)
# Won't work -- f is not a function.
# Y = f(X)
# Works.

# ! Establecer las soluciones del polinomio como funciones de k para representarlas gráficamente 
# np.seterr(divide='ignore', invalid='ignore')

# Y = lambdify(k, (solucion[0]), "numpy")(X)
# Z = lambdify(k, (solucion[1]), "numpy")(X)
# W = lambdify(k, (solucion[2]), "numpy")(X)
# U = lambdify(k,stangent , "numpy")(X)

# Y1d = lambdify(k, (solucion1d[0]), "numpy")(X)
# Z1d = lambdify(k, (solucion1d[1]), "numpy")(X)
# W1d = lambdify(k, (solucion1d[2]), "numpy")(X)




fig = plt.figure()
# plt.plot(X, (Y),color='C0',label=r'$\lambda_1$ ')
plt.plot(X, f.Panel(3).sol1(X,alfa),color='C0',label=r'$\lambda_1(d=3)$ ')
plt.plot(X, f.Panel(d).sol1(X,alfa),color='C0',linestyle = "--",label=r'$\lambda_1(d=2)$ ')
# plt.plot(X, f.Panel(1).sol11d(X,alfa),color='C0',linestyle = ":" ,label=r'$\lambda_1(d=1)$ ')
plt.plot( f.Panel(2).kc(alfa),0.0, marker="o", markersize=20, markeredgecolor="red", markerfacecolor="green")
plt.plot( f.Panel(1).kc1d(alfa),0.0, marker="o", markersize=20, markeredgecolor="red", markerfacecolor="green")

plt.plot(X, f.Panel(3).sol2(X,alfa),color='C1',label=r'$\lambda_2(d=3)$ ')
plt.plot(X, f.Panel(d).sol2(X,alfa),color='C1',linestyle = "--",label=r'$\lambda_2(d=2)$ ')
plt.plot(X,f.Panel(1).sol21d(X,alfa),color='C1',linestyle = ":" ,label=r'$\lambda_2(d=1)$ ')


plt.plot(X, f.Panel(3).sol3(X,alfa),color='C2',label=r'$\lambda_3(d=3)$ ')
plt.plot(X, f.Panel(d).sol3(X,alfa),color='C2',linestyle = "--",label=r'$\lambda_3(d=2)$ ')
plt.plot(X, f.Panel(1).sol31d(X,alfa),color='C2',linestyle = ":" ,label=r'$\lambda_3(d=1)$ ')

# plt.plot(X, f.Panel(1).sol_t(X,alfa),color='C3',linestyle = ":",label=r'$\lambda_{||}(d=1)$ ')


plt.plot(X, f.Panel(3).sol_t(X,alfa),color='C3',label=r'$\lambda_{||}(d=3)$ ')
plt.plot(X, f.Panel(d).sol_t(X,alfa),color='C3',linestyle = "--",label=r'$\lambda_{||}(d=2)$ ')

# plt.plot(X, (U),color='C3',label=r'$\lambda_{||}$ ')
plt.legend(loc=0,fontsize=30)
plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)

plt.xlabel( r' $k$ ', fontsize=30)










# ! plot hydrodynamics coefficients
alpha = np.linspace(0.0,0.99,1000)
fig2 = plt.figure()
plt.plot(alpha, f.Panel(3).kapa(alpha),color='C0',label=r'$\kappa(\alpha), \, d=3$ ')
plt.plot(alpha, f.Panel(3).mu(alpha),color='C1',label=r'$\mu(\alpha),\, d=3$ ')
plt.plot(alpha, f.Panel(d).kapa(alpha),color='C0',linestyle ="--",label=r'$\kappa(\alpha), \, d=2$ ')
plt.plot(alpha, f.Panel(d).mu(alpha),color='C1',linestyle ="--",label=r'$\mu(\alpha),\, d=2$ ')
plt.plot(alpha, f.Panel(1).kapa1d(alpha),color='C0',linestyle =":",label=r'$\kappa(\alpha), \, d=1$ ')
plt.plot(alpha, f.Panel(1).mu1d(alpha),color='C1',linestyle =":",label=r'$\mu(\alpha),\, d=1$ ')

plt.legend(loc=0,fontsize=30)
plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)

plt.xlabel( r' $\alpha$ ', fontsize=30)

plt.show()