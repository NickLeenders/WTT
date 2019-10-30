#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 14:32:55 2019

@author: nickleenders
"""

import numpy as np
import cmath
import math
import matplotlib.pyplot as plt

"""From Task 2 (Generator):"""
Vnom = 690
Pnom = (4.2*10**6)
nPoles = 40
Ke = 7.8019
Ls = 1.2*10**-3
Rs= 4*10**-3
Kmech = 373963.3497

RPM = np.linspace(3, 15, 50)
omegaMech = 2*np.pi*RPM/60
omegaEl = omegaMech*nPoles
Xs = 1j*omegaEl*Ls
Z = Rs + Xs

tol = 1E-15

Pmech = Kmech*omegaMech**3
dVals = np.arcsin(2*Pmech*Ls/(Ke**2*omegaEl))/2
Ig = Ke*np.sin(dVals)/Ls

Vg = Ke*omegaEl*np.sin(dVals)-Ig*Rs
"""End Task 2"""

"""Transformer and LOC"""
L1 = 2*10**-6 #unit H
R1 = 2*10**-3 #unit Ohm
L2 = 2*10**-6 #unit H
R2 = 2*10**-3 #unit Ohm


Z1 = R1 + L1*1j*omegaEl
Z2 = R1 + L1*1j*omegaEl

R_cable = 1 #Unit Ohm
L_cable = 5*10**-3 #Unit H
C_cable = 1*10**-6 #Unit F

Xl_cable = L_cable*omegaEl
Xc_cable = -1/(omegaEl*C_cable)

#Z_cable = R_cable*Xl_cable
Z_cable = (R_cable*Xc_cable**2)/(R_cable**2+(Xl_cable+Xc_cable)**2)+1j*((R_cable**2*Xc_cable)+(Xc_cable*Xl_cable)*(Xc_cable-Xl_cable))/(R_cable**2+(Xl_cable+Xc_cable)**2)

Ztotal = Z1+Z2+Z_cable

Vlow = 690 #unit V
Vhigh = 33*10**3 #unit V
N_ratio = Vlow/Vhigh # (N1/N2) no unit
I1 = Ig # unit W from Task 2
I2 = I1*N_ratio
P2 = I2**2*Ztotal.real
S2 = I2**2*(Ztotal.real**2+Ztotal.imag**2)
Q2 = I2**2*Ztotal.imag


PF = np.divide(P2,S2)



plt.plot(RPM,abs(P2))
plt.plot(RPM,abs(Q2))
plt.xlabel('n (RPM)')
plt.ylabel('Power (W)')
plt.legend(['Active power','Reactive power'])
