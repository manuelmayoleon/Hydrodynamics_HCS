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
alfa = 0.999
d= 1


z=f.Panel(d).zeta(alfa)
nu1=f.Panel(d).nu1(alfa)
kapa=f.Panel(d).kapa(alfa)

print("zeta")
print(z)
print("kapa")
print(kapa)



alpha = np.linspace(0.0,0.88,1000)

plt.plot(alpha, f.Panel(d).kapa(alpha),color='C0',label=r'$\kappa(\alpha), \, d=1$ ')
plt.plot(alpha, f.Panel(d).mu(alpha),color='C1',label=r'$\mu(\alpha),\, d=1$ ')
plt.plot(alpha,2*f.Panel(d).zeta(alpha)*(f.Panel(d).kapa(alpha)-(d-1)*f.Panel(d).a2d(alpha)/(d*f.Panel(d).zeta(alpha)))/(2*f.Panel(d).nu2(alpha)-3*f.Panel(d).zeta(alpha)),color='C2',label=r'$2\nu_2(\alpha)-3\zeta(\alpha),\, d=1$ ')

# plt.plot(alpha,Panel(1).kc(alpha),color='C0',label=r'$k_c(\alpha),\, d=1$ ')


# plt.plot(alpha, a2(alpha),color='C0',label=r'$a_2( \alpha )\, d=1$ ')

# plt.plot(a, Panel.a2d(a,2),color='C2',label=r'$a_2( \alpha )\, d=2$ ')

# plt.plot(a, Panel.a2d(a,3),color='C3',label=r'$a_2( \alpha )\, d=3$ ')
plt.legend(loc=0,fontsize=30)
plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)

plt.xlabel( r' $\alpha$ ', fontsize=30)

plt.show()