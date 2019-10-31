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
from Assignment_2_Q2 import Vg,Ig,omegaEl,RPM


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

Xl_cable = L_cable*omegaEl #Xl for 
Xc_cable = -1/(omegaEl*C_cable) #

Z_cable = (R_cable*Xc_cable**2)/(R_cable**2+(Xl_cable+Xc_cable)**2)+1j*((R_cable**2*Xc_cable)+(Xc_cable*Xl_cable)*(Xc_cable-Xl_cable))/(R_cable**2+(Xl_cable+Xc_cable)**2) #Cable connection impedance

Ztotal = Z1+Z2+Z_cable # Total impedance

Vlow = 690 #unit V
Vhigh = 33*10**3 #unit V
N_ratio = Vlow/Vhigh # (N1/N2) no unit
I1 = Ig # unit W from Task 2
I2 = I1*(Ztotal.real**2+Ztotal.imag**2) # Current after transformer and LOC
P2 = I2**2*Ztotal.real # True power
S2 = I2**2*(Ztotal.real**2+Ztotal.imag**2) #Reactive power
Q2 = I2**2*Ztotal.imag # Appartent power


PF = np.divide(P2,S2) # Power factor

plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window

plt.plot(RPM,abs(P2))
plt.plot(RPM,abs(Q2))
plt.xlabel('n (RPM)')
plt.ylabel('Power (W)')
plt.legend(['Active power','Reactive power'])
