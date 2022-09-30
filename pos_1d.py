from numpy.core.function_base import linspace
import pandas as pd
import numpy as np
import scipy . stats as ss
import math
import matplotlib . mlab as mlab
import matplotlib . pyplot as plt
import matplotlib.animation as manimation
from scipy import optimize
import matplotlib
matplotlib.rcParams['text.usetex'] = True
from matplotlib.transforms import (
    Bbox, TransformedBbox, blended_transform_factory)
import os 





# r=np.genfromtxt("pos.txt",names=["ry","rz"])
# v = np.genfromtxt("velocidad.txt", names=["vy","vz"])
# rinicial=np.genfromtxt("posiciones_init.txt",names=["ry","rz"])
# vinit=np.genfromtxt("velocidad_init.txt", names=["vy","vz"])
# temp=np.genfromtxt("temperaturas.txt" ,names= ["y","z"])
# tiempo=np.genfromtxt("tiemposdecol.txt",names=["t"])
path = os.getcwd()


print(path)
# for filename in os.listdir(path):
#     if filename.startswith('pos_0.995_tiempo'):
#         ri[i] = pd.read_csv(filename,header=None,sep='\s+',names=["ryy","rzz"])
#     i=i+1   

# r= pd.read_csv("pos_0.995_tiempo_50000000.txt",header=None,sep='\s+',names=["ry","rz"])
# r= pd.read_csv("pos_0.995_tiempo_50000000.txt",header=None,sep='\s+',names=["ry","rz"])
r= pd.read_csv("pos_0.990colpp_100.000.txt",header=None,sep='\s+',names=["ry"])
v =  pd.read_csv("vel_0.990colpp_100.000.txt" ,header=None,sep='\s+' , names=["vy"])
rinicial= pd.read_csv("posiciones_init.txt",header=None,sep='\s+',names=["ry"])
vinit= pd.read_csv("velocidad_init.txt",header=None,sep='\s+', names=["vy"])
temp= pd.read_csv("temperaturas_0.990colpp_100.000.txt",header=None,sep='\s+' ,names= ["y"])
tiempo= pd.read_csv("tiemposdecol_0.990colpp_100.000.txt",header=None,sep='\s+',names=["t"])




# print(rho)



# cols = 2795
# alfa = 0.950


#*Esta parte simplemente sirve para poder plotear a la vez
#* el fit y el histograma
num_bins =100
fig22,ax = plt.subplots (figsize=(10,10))

plt.xlabel ( r' $v_i$ ', fontsize=20)
plt.ylabel ( r' Frecuencia ',fontsize=20)

plt.title ( r' \textbf {Histograma de la velocidad en el eje y}  ',fontsize=30)
# plt.xlim (1 ,9)

#! Hacer el histograma 
n,bins,patches = ax.hist(v['vy'],num_bins,density ='true',facecolor ='C0',edgecolor='white',label='$v$ ')
n1,bins1,patches1 = ax.hist(vinit['vy'],num_bins,density ='true',facecolor ='C1',edgecolor='white',alpha=0.5,label='$v_{init}$ ')

# ,edgecolor='yellow'
plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
#! Ajuste gaussiano 
N=100

def gaussian(x,mu,sig):
    return np.exp(-np.power(x-mu,2.)/(2*sig))/np.sqrt(2*np.pi*sig)
x1 =np.linspace(min(v['vy']),max(v['vy']),N)
x2 =np.linspace(min(vinit['vy']),max(vinit['vy']),N)

#Hacer el fiting de la ley de potencias

# # def fit_func (x ,a , b ) :
# #     return a * x **( b )

params , params_covariance = optimize.curve_fit(gaussian,bins[1:],n,method='dogbox')

params1 , params_covariance1 = optimize.curve_fit(gaussian,bins1[1:],n1,method='dogbox')

# # print('param')
# # print (params)
plt.plot(x1,gaussian(x1,params[0],params[1]),color='C3',label='Ajuste Gaussiano de $v_y$')
plt.plot(x2,gaussian(x2,params1[0],params1[1]),color='C2',label='Ajuste Gaussiano de $v_{init}$')

plt.legend(loc=0,fontsize=20)
# plt.xlim(-5,5)    

# plt.savefig (path+'/histograms/histogram_dist.pdf',format ='pdf')

# def fit_func (x ,a , b ) :
#     return a * x+ b 
fig1,ax1 = plt.subplots(1,1)

plt.plot(tiempo["t"],temp["y"],label="$T_y$")

plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
plt.xlabel ( r' $t$ ', fontsize=20)
plt.ylabel ( r' $T$ ',fontsize=20)
plt.legend(loc=0,fontsize=20)
# plt.xlim(0,1500)
# plt.ylim(0.0,1.0)
# plt.savefig (path+'/temps/temp%i'%cols+'alfa%1.3f'%alfa+'epsilon%1.3f'%epsilon+'.pdf',format ='pdf')


fig3, ax3 = plt.subplots(1,1)
    # plt.plot(r["ry"][0:len(r["ry"])-1],r["rz"][0:len(r["ry"])-1], "o")
plt.plot(r["ry"][0:len(r["ry"])-1],v["vy"][0:len(r["ry"])-1], "o",c='C0')
# plt.plot(rinicial["ry"],vinit["vy"], "o")
# # params , params_covariance = optimize.curve_fit(fit_func,corr['tau'],np.log(corr['corr']),method='dogbox')
# # print('param')
# # print (params)
# # print(params_covariance)
# # print(np.exp(params[0]))
# ax1 = plt.subplots(1,1)
# # plt.yscale("log")
# tau=np.linspace(0,2.5,100)
# plt.plot(corr['tau'],np.log(corr['corr']),'o')
# # plt.plot(tau,fit_func(tau,params[0],params[1]))
# plt.xlim(-7500,7500)
# plt.ylim(0.5,1.0)
plt.xlabel ( r' $y$ ', fontsize=20)
plt.ylabel ( r' $v_y$ ',fontsize=20)
# plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
plt.title ( r' \textbf {Espacio de fases }  ',fontsize=30)

plt.show()