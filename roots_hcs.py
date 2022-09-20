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





k= Symbol("k")
x = Symbol('x')
kapa= Symbol("self.kapa(alpha)",positive=True)
kapa1d= Symbol("self.kapa1d(alpha)",positive=True)
k = Symbol("k")
mu=Symbol("self.mu(alpha)",positive=True)
mu1d=Symbol("self.mu1d(alpha)",positive=True)
eta= Symbol("self.eta(alpha)")
zet= Symbol("self.zet(alpha)")
d = Symbol("self.d")


# ! Polinomio característico de las ecuaciones hydrodinamicas en d dimensiones
poly_d  =  x**3 + ( ( d + 2 ) * kapa / ( 2 * ( d -1 ) ) 
         
         + ( d - 1 ) * eta / d) *  k**2 * x**2 + ( ( d + 2 ) * k**4 * zet * kapa / (2 * d ) 
         
         + ( k**2 / d )* ( d + 2 + ( d - 1 ) * zet * eta - d * ( d + 2 ) * zet * kapa / ( 2 * ( d - 1 ) ) )
         
         - zet**2 ) * x + ( ( d + 2 ) * ( kapa - mu )* k**2 / (2 * ( d - 1 ) ) - zet ) * k**2  
         
# print(poly_analitic)

# ! Establecer el polinomio igualado a 0 para resolverlo posteriormente 
quadratic_equation = Eq(poly_d, 0)
# ! Resolver polinomio en función de x 
solucion=solve(quadratic_equation, x)

# ! printear las raices del polinomio 
# print("solucion 1")
# print(solucion[0])
# print("solucion 2")
# print(solucion[1])
# print("solucion 3")
# print(solucion[2])
# ! Raiz para el modo transversal 
# stangent = z - 0.5*eta*k**2


poly_1d  =  x**3 + ( ( 3 ) * kapa1d / ( 2  ) ) *  k**2 * x**2 + ( ( 3 ) * k**4 * zet * kapa1d / (2 ) 
         
         + ( k**2  )* ( 3  - 3  * zet * kapa1d / ( 2  ) )
         
         - zet**2 ) * x + ( ( 3 ) * ( kapa1d - mu1d )* k**2 / (2  ) - zet ) * k**2 
# print(poly_1d)

# ! Establecer el polinomio igualado a 0 para resolverlo posteriormente 
quadratic_equation = Eq(poly_1d, 0)
# ! Resolver polinomio en función de x 
solucion1d=solve(quadratic_equation, x)
# ! printear las raices del polinomio 
print("solucion 1")
print(solucion1d[0])
print("solucion 2")
print(solucion1d[1])
print("solucion 3")
print(solucion1d[2])
