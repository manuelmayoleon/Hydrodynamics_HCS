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


#!/usr/bin/python
# -*- coding: ascii -*-
import os, sys
##!! AVISOS URGENTES 
##** CONSEJOS 
##???? EXPLICACIONES 
##///  tachado

density= pd.read_csv("densityprom_0.990.txt",header=None,sep='\s+', names=['re','im'])

t= pd.read_csv("tiemposdecol_0.990colpp_500.0.txt",header=None,sep='\s+', names=['t'])

# regd= pd.read_csv("reg_dens.txt",header=None,sep='\s+', names=['dens','er','alfa'])
info= pd.read_csv("data.txt",header=None,sep='\s+')


info.columns = info.iloc[0]

colp= float(info['colisiones.p.p.'].iloc[-1])


print("colisiones por partícula")

print(colp)


 


# cols= 2473
# alfa= 0.990


# print(data[12])

# ss= pd.read_csv("sparam_0.95.txt",header=None,sep='\s+', names=['s'])





##!!Ajuste lineal del perfil de densidad obtenido con MD ###

min_index = np.argmin(density['re'])
max_index = np.argmax(density['re']) 
# print(min_index) 
# print(max_index) 


##* DETERMINACION DE DENSIDAD, SU ESCALA LOG
density=density['re']+1j*density['im']
denslog=np.log(density)

x=t['t']

##* Eliminamos el régimen anterior y posterior al incremento exponencial 

colpp= np.linspace(0,colp,len(denslog)) ##! n. de col. para la última simulacion 
logcolpp= np.linspace(0,colp,len(denslog)) ##! n. de col. para la última simulacion 

le= len(np.where(colpp<=200)[-1])
lh=0
lh= len(np.where(colpp<=300)[-1])
if lh != 0:  
    logcolpp=colpp[le-1:lh-1]

    denslog=denslog[le-1:lh-1]
else: 
    logcolpp=colpp[le-1:]
  
    denslog=denslog[le-1:]
#** Es la regresion en escala temporal
#!! reg = LinearRegression().fit(x.values.reshape((-1, 1)),np.real(denslog))

##?? Regresion lineal de los datos obtenidos para n_k a traves de MD

reg = LinearRegression().fit(logcolpp.reshape(-1,1),np.real(denslog))

r_sq = reg.score(logcolpp.reshape((-1, 1)), np.real( denslog))
inferred= reg.predict(logcolpp.reshape((-1, 1)))
model_error = mean_absolute_error(np.real(denslog), inferred)
print('coefficient of determination:', r_sq)
print('intercept:', reg.intercept_)
print('slope:', reg.coef_)
print('model error:',model_error)
linear_reg=reg.coef_*logcolpp+reg.intercept_

coef=reg.coef_

#?? Escribir en un archivo .txt la pendiente con su error y el coef de restitucion
# with open('reg_dens.txt', 'a',newline='\n') as f:
#     writer = csv.writer(f, delimiter='\t')
#     writer.writerows(zip( coef,[model_error],[alfa]))
  

##?? Prediccion teorica de la pendiente 
# tiemp=np.linspace(0.0,100.0,len(t['t']))
# densteo=np.exp(functions.eigenvalue1(functions.Panel(1).mu(alfa,lin_dens),functions.Panel(1).kapa(alfa,lin_dens),functions.lamda1(alfa,epsilon),k)*tt*dens/20)




##!!Representación del perfil de densidad ###

fig = plt.figure()
# plt.plot(tt,densteo,color='C2',label="$n_y$ ")
plt.plot(colpp,density,color='C1',label="$n_{\frac{\pi}{L}}$ (MD)")   
# plt.plot(colpp,denslog,color='C1',label="$n_{\frac{\pi}{L}}$ (MD)")   
#plt.plot(colpp,linear_reg,color='C2',label="$n_{\frac{\pi}{L}}$ (MD)")   
# print((functions.eigenvalue3(functions.Panel(1).mu(alfa),functions.Panel(1).kapa(alfa),functions.lamda1(alfa,epsilon),kk)).real)
# plt.plot(tt,expdens,color='C2',label="$n_y \;(k=2\pi/L)$ ")

# print(functions.Panel(1,sigma,epsilon,vp).T_s(alfa,lin_dens))
plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
plt.xlabel ( r'$s$', fontsize=30)
plt.ylabel ( r' $n_2$ ',rotation=0.0,fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

# plt.title ( r' \textbf {Autovalores en función de $\alpha$}' ,fontsize=40)
plt.legend(loc=0,fontsize=30)

fig11 = plt.figure()
plt.plot(logcolpp,denslog,color='C1',label="$n_{\frac{\pi}{L}}$ (MD)")   
plt.plot(logcolpp,linear_reg,color='C2',label="$n_{\frac{\pi}{L}}$ (MD)")   
# print((functions.eigenvalue3(functions.Panel(1).mu(alfa),functions.Panel(1).kapa(alfa),functions.lamda1(alfa,epsilon),kk)).real)
# plt.plot(tt,expdens,color='C2',label="$n_y \;(k=2\pi/L)$ ")

# print(functions.Panel(1,sigma,epsilon,vp).T_s(alfa,lin_dens))
plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
plt.xlabel ( r'$s$', fontsize=30)
plt.ylabel ( r' $n_2$ ',rotation=0.0,fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)



plt.legend(loc=0,fontsize=30)





plt.show()

