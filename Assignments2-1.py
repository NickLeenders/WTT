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
n_poles = 80                            #no. of generator poles
k_e = 5.34                              #From previous calculations
n_rpm = np.linspace(3,15,13)            #rpm range
omega_e = n_rpm*2*np.pi/60*n_poles      #Electrical frenquency
R_s = 4*1e-3
L_s = 1.2*1e-3

#Impedance
Z_s = R_s + 1j*omega_e*L_s
#EMF voltage (rms)
V_EMF = k_e*omega_e
I_gen = V_EMF/Z_s
"""End Task 2"""

"""Transformer and LOC"""
L1 = 2*10**-6 #unit H
R1 = 2*10**-3 #unit Ohm
L2 = 2*10**-6 #unit H
R2 = 2*10**-3 #unit Ohm


Z1 = R1 + L1*1j*omega_e
Z2 = R1 + L1*1j*omega_e

R_cable = 1 #Unit Ohm
L_cable = 5*10**-3 #Unit H
C_cable = 1*10**-6 #Unit F

Xl_cable = math.sqrt(3)*L_cable*omega_e

Z_cable = R_cable*Xl_cable

Xc_cable = math.sqrt(3)*-1*(omega_e*C_cable)

Ztotal = Z1+Z2+Z_cable

Vlow = 690 #unit V
Vhigh = 33*10**3 #unit V
N_ratio = Vlow/Vhigh # (N1/N2) no unit
I1 = 4.2*10**6/3/(Vlow/math.sqrt(3)) # unit W from Task 2
I2 = I1*Ztotal
V2 = Vlow*Ztotal
P2 = I2*V2

R_cable = 1 #Unit Ohm
L_cable = 5*10**-3 #Unit H
C_cable = 1*10**-6 #Unit F
Xl_cable = math.sqrt(3)*L_cable*omega_e
Xc_cable = math.sqrt(3)*-1*(omega_e*C_cable)

plt.plot(n_rpm,abs(P2.real))
plt.plot(n_rpm,abs(P2.imag))
plt.xlabel('n (RPM)')
plt.ylabel('Power (W)')
plt.legend(['Active power','Reactive power'])
