#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 11:36:30 2019

@author: nickleenders
"""

import numpy as np
import cmath
import math
import matplotlib.pyplot as plt
from Assignment_2_Q2 import Vg,Ig,omegaEl,RPM
plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window

L1 = 2*10**-6 #unit H (Transformer)
R1 = 2*10**-3 #unit Ohm (Transformer)
L2 = 2*10**-6 #unit H (Transformer)
R2 = 2*10**-3 #unit Ohm (Transformer)

R_cable = 1 #Unit Ohm
L_cable = 5*10**-3 #Unit H
C_cable = 1*10**-6 #Unit F

Xl_cable = L_cable*omegaGrid #Xl for cable connection
Xc_cable = -1/(omegaGrid*C_cable) # Xc for cable connection

# Circuit 1
Vpoc = 33*10**3*np.sqrt(3)
omegaGrid = 2*np.pi*50
Z1 = R1 + L1*1j*omegaGrid
Z2 = R1 + L1*1j*omegaGrid
Z_transformer = Z1+Z2
Z_cable = omegaGrid*C_cable
Z_circuit1 = Z1+Z2+R_cable+Xl_cable*1j
Ipoc = Ig

# Circuit 2
Z_circuit2 = Z1+Z2+R_cable+(Xl_cable+Xc_cable)*1j

# Cable
Ic = Vpoc/(1j/(omegaGrid*C_cable))+Ipoc
Vc = Ic*(R_cable+1j*omegaGrid*L_cable)+Vpoc

P2 = Ic.real*Vc.real
Q2 = abs(Ic.imag*Vc.imag)

# Transformer
Vlow = 690 #unit V
Vhigh = 33*10**3 #unit V
N_ratio = Vlow/Vhigh # (N1/N2) no unit
Vt = Vc*N_ratio
It = Ic*N_ratio


plt.plot(RPM,abs(P2))
plt.plot(RPM,abs(Q2))
plt.xlabel('n (RPM)')
plt.ylabel('Power (W)')
plt.legend(['Active power','Reactive power'])