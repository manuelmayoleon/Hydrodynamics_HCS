#functions_d.py

import numpy as np
import pandas as pd
import scipy.special as special
from numba import jit
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True
from matplotlib.transforms import (
    Bbox, TransformedBbox, blended_transform_factory)
# Import math Library
import math 

from sympy import *

# PARA USAR A LA HORA DE GUARDAR DOS COLUMNAS EN UN ARCHIVO
import csv

import cmath

class Panel:
    def __init__(self,d=1,sigma=1.0,epsilon= 0.5 ,vp=0.0001):
        self.d = d
        self.sigma = sigma
        self.vp=  vp
        self.epsilon = epsilon 
        self.eta0 = (2+self.d)*np.sqrt(np.pi/2)/8
        self.kapa0= (2+self.d)*self.eta0/2
    def a2d(self,alpha):
            return  16*(1-2*alpha**2.00)*(1-alpha)/(9+24*self.d+(8*self.d-41)*alpha+30*alpha**2.00*(1-alpha))
    def zeta(self, alpha):
            return ((2+self.d)/(4*self.d))*(1-alpha**2.00)*(1+3*self.a2d(alpha)/16)
    
    def nu1(self, alpha):
            return (3-3*alpha+2*self.d)*(1+alpha)*(1-self.a2d(alpha)/32)/(4*self.d)
    def nu2(self,alpha):
            return (1+alpha)*((self.d-1)/2+3*(self.d+8)*(1-alpha)/16+(4+5*self.d-3*(4-self.d)*alpha)*self.a2d(alpha)/512)    
    
    def eta(self,alpha):
                return 1/(self.nu1(alpha)-self.zeta(alpha)/2)       
    def kapa(self,alpha):
        return (1+2*self.a2d(alpha))/(self.nu2(alpha)-2*self.d*self.zeta(alpha))*(self.d-1)
    def mu(self,alpha):
        return 2*self.zeta(alpha)* (self.kapa(alpha)+(self.d-1)*self.a2d(alpha)/(self.d*self.zeta(alpha)))/(2*self.nu2(alpha)/self.d-3*self.zeta(alpha))
      
            
    def zeta_max(self, alpha):
            return ((2+self.d)/4*self.d)*(1-alpha**2.00)
    def gamma(self,alpha):
            return (4*alpha*self.epsilon**2.00+12*(1-alpha))/((1+3*alpha)*self.epsilon**2.00)
    def T_s(self,alpha,density):
            return ((3*np.sqrt(np.pi)*self.gamma(alpha))/((1+alpha)*(self.gamma(alpha)-(1+alpha)/2)*self.epsilon**3.00*density*self.sigma))**2.00*self.vp**2
    
    def nu_max(self,alpha):
            return (1+alpha)*((self.d-1)/2+3*(self.d+8)*(1-alpha)/16)
    def eta(self,alpha):
            return 1/(self.nu2(alpha)-self.zeta(alpha)/2) 
   
    def kapa_max(self,alpha):
            return (1)/(self.nu_max(alpha)-2*self.d*self.zeta_max(alpha))*self.kapa0
    def mu_max(self,alpha):
            return 2*self.zeta_max(alpha)*(self.kapa_max(alpha))/(2*self.nu_max(alpha)/self.d-3*self.zeta_max(alpha))
    def kc(self,alpha):
            return np.sqrt(2/(self.d+2)*(self.kapa(alpha)-self.mu(alpha)))
    def factor(self,lin_dens,rho,alpha):
         return (lin_dens*np.sqrt(np.pi/2)/(rho*self.epsilon*(1+alpha)*self.sigma))*(2*np.sqrt( self.T_s(alpha,lin_dens)))/(2+self.T_s(alpha,lin_dens))
