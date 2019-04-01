# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 21:15:37 2018

@author: owner
"""

import os, subprocess, sys, time
from numpy import pi, array, zeros,linspace,sqrt,arctan,sin,cos,tanh
from matplotlib.pyplot import plot,show,subplots,savefig,xlabel,ylabel,clf,close,xlim,ylim
from scipy.integrate import odeint



# v points
v = linspace(0,1)

def K(v):
    k0 = 12
    p = 10
    return k0*(1 + v**p/(1.01-v**p))
    
def invK(v):
    invk0 = 1./12 
    w = 0.13
    return invk0*tanh((1-v)/w)

# solve ODE

#print(y[0])
# plot results
clf;close()
#plot(v,K(v))
plot(v,1/K(v))
plot(v,invK(v))
xlim([0,1])
#ylim([0,4])
show()