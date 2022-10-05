import numpy as np
import pandas as pd
from numba import jit
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True
from matplotlib.transforms import (
    Bbox, TransformedBbox, blended_transform_factory)
from sklearn.linear_model import LinearRegression

from sklearn.metrics import mean_absolute_error



temp= pd.read_csv("temperaturas_0.990colpp_500.0.txt" ,header=None,sep='\s+' ,names= ["y"])
tiempo= pd.read_csv("tiemposdecol_0.990colpp_500.0.txt",names=["t"])
info= pd.read_csv("data.txt",header=None,sep='\s+')


alfa = 0.990




info.columns = info.iloc[0]

colp= float(info['colisiones.p.p.'].iloc[-1])
tini= float(info['T'].iloc[-1])
# colp = 500.0

print("colisiones por partícula")

print(colp)



# n=500
npar = 1000
sigma=1.0e-06

# alfa=0.950
# epsilon=0.5

# lin_dens=n/l            
l = 10 

rho=npar*sigma/l

k=2*np.pi/l



# print(datos)                   
#?? Datos para el método RK


# y0=5.0
# x0=ts+k
# y0=gamma*x0

#?? resolucion de las ecuaciones en funcion del tiempo
x0= temp['y'][0] /temp['y'][0] 
a=0
# b=int(tiempo["t"][len(tiempo)-1])
b= colp
# b=5

#??Resolucion de las ecuaciones en colpp
# x0= temp['y'][0]/vp**2
# y0=temp['z'][0]/vp**2
# a=0
# b=int(cols)

# b=30050
h=0.002

# print(temp)
# print(tiempo)

class Panel:
    def __init__(self,alfa=0.95,rho=0.03,sigma=1.0, npar = 500, l = 10 ):
        self.rho = rho
        self.sigma = sigma
        self.alfa = alfa
        self.npar = npar 
        self.l = l
    def f(self,t1,t):
        
        # return np.sqrt(np.pi)*(1+alfa)*epsilon*rho*np.sqrt(t1)*( -(1-alfa)*t1+epsilon**2.0*(-(5*alfa-1)*t1 +(3*alfa+1)*t2 )/12)
        # return 2.0*(1+self.alfa)*self.rho*np.sqrt(t1)* (-(1-self.alfa)*t1)/self.sigma
        # return -2.0*(1-self.alfa**2)*self.npar*np.sqrt(t1)*(t1)/self.l
        return -(1-self.alfa)*(t1)
       
        # return 4.0*epsilon**3.0*rho*np.sqrt(tx)*( ty -tx )/(3.0*np.sqrt(np.pi))

 
    
    
        
    # return -4.0*epsilon**3.0*rho*np.sqrt(tx)*( ty -tx )/(3.0*np.sqrt(np.pi))
#?? Representacion en funcion de las colisiones por particula
colpp=np.linspace(0,colp,len(temp["y"]))
plt.plot(colpp[0:500*npar],temp["y"][0:500*npar]/temp["y"][0],color='C1',label="$T$ (MD)")      
plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
plt.xlabel ( r' $s$ ', fontsize=30)
plt.ylabel ( r' $\tilde{T}$ ',rotation=0.0,fontsize=30)

# plt.show()

# #?? Representacion en funcion del tiempo de colision
# time = np.zeros(len(tiempo["t"]))
# time = tiempo["t"] - tiempo["t"][0]
# plt.plot(time,temp["y"]/temp["y"][0],color='C1',label="$T_y$ (MD)")      
# plt.grid(color='k', linestyle='--', linewidth=0.5,alpha=0.2)
# plt.xlabel ( r' $t(T_0/m\sigma^2)^{1/2}$ ', fontsize=30)
# plt.ylabel ( r' $T$ ',rotation=0.0,fontsize=30)



plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

plt.title ( r' \textbf {Solución de las ecuaciones para las temperaturas}  ',fontsize=40)


# @jit
@jit(forceobj=True)
def runge_kutta_system(f, x0, a, b, h):
    t = np.arange(a, b + h, h)
   
    n = len(t)
    x = np.zeros(n)
    # print(x)
    x[0] = x0
    
    for i in range(n - 1):
        k1 = h * f(x[i], t[i])
       
        k2 = h * f(x[i] + k1 / 2,  t[i] + h / 2)
       
        k3 = h * f(x[i] + k2 / 2,  t[i] + h / 2)
        
        k4 = h * f(x[i] + k3, t[i] + h)
       
        x[i + 1] = x[i] + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + 2 * k4)
    
    plt.plot(t, x ,color='C2',label='$T_y$ ')
    
    
    # reg = LinearRegression().fit(t.reshape(-1,1),np.real(np.log(x)))

    # r_sq = reg.score(t.reshape((-1, 1)), np.real( np.log(x)))
    # inferred= reg.predict(t.reshape((-1, 1)))
    # model_error = mean_absolute_error(np.real(np.log(x)), inferred)
    # print('coefficient of determination:', r_sq)
    # print('intercept:', reg.intercept_)
    # print('slope:', reg.coef_)
    # print('model error:',model_error)
    # linear_reg=reg.coef_*t+reg.intercept_

    
    print(min(x))
    plt.legend(loc=0,fontsize=30)
    plt.show()
# np.seterr('raise')
# #?? Aplicación RK  a ecuacion en escala de colpp
# runge_kutta_system(Panel(alfa,epsilon,rho,vp).f_colpp,Panel(alfa,epsilon,rho,vp).g_colpp,x0,y0,a,b,h)
#?? Aplicación RK  a ecuacion en escala temporal 
runge_kutta_system(Panel(alfa,rho,npar,l).f,x0,a,b,h)



# reg = LinearRegression().fit(colpp.reshape(-1,1),np.real(np.log(temp["y"]/temp["y"][0])))

# r_sq = reg.score(colpp.reshape((-1, 1)), np.real( np.log(temp["y"]/temp["y"][0])))
# inferred= reg.predict(colpp.reshape((-1, 1)))
# model_error = mean_absolute_error(np.real(np.log(temp["y"]/temp["y"][0])), inferred)
# print('coefficient of determination:', r_sq)
# print('intercept:', reg.intercept_)
# print('slope:', reg.coef_)
# print('model error:',model_error)
# linear_reg=reg.coef_*colpp+reg.intercept_

# coef=reg.coef_




# plt.show()