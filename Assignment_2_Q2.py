# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:23:57 2019

"""
import numpy as np
import matplotlib.pyplot as plt
import cmath
from scipy.optimize import fsolve


Vnom = 690
Pnom = (4.2*10**6)
nPoles = 40
Ke = 7.8019
Ls = 1.2*10**-3
Rs= 4*10**-3
Kmech = 373963.3497

RPM = np.linspace(3, 15, 50)
#RPM = np.array([4])
omegaMech = 2*np.pi*RPM/60
omegaEl = omegaMech*nPoles
Xs = 1j*omegaEl*Ls
Z = Rs + Xs

tol = 1E-15
#dVals = 0*RPM

#dVals = 0*RPM
#dValstest = 0*RPM
#d0 = np.pi/10
#for i in range(len(omegaEl)):
#    minFunc = lambda d: Kmech*omegaMech[i]**3-\
#                    Ke**2*omegaEl[i]*np.sin(2*d)/Ls
#    di = fsolve(minFunc, d0)
#    dValstest[i] = di

Pmech = Kmech*omegaMech**3

# From solving the SoE: (see mathematica)
dVals = np.arcsin(2*Pmech*Ls/(Ke**2*omegaEl))/2
Ig = Ke*np.sin(dVals)/Ls

Vg = Ke*omegaEl*np.cos(dVals)-Ig*Rs
PlossPhase = Rs*Ig**2
Ploss = 3*PlossPhase

PlossCheck = 3*Kmech*(2*np.pi*15/60)**3 - Pnom


plt.figure()
plt.plot(RPM, Ploss)
plt.xlabel('RPM')
plt.ylabel('$P_{loss} [W]$')
plt.grid()

plt.figure()
plt.plot(RPM, Vg)
plt.xlabel('RPM')
plt.ylabel('$Vg [V]$')
plt.grid()

plt.figure()
plt.plot(RPM, Ig)
plt.xlabel('RPM')
plt.ylabel('Ig [A]')
plt.grid()

