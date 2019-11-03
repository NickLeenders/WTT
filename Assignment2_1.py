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
from Assignment_2_Q2 import Vg,Ig,omegaEl,RPM,Ploss
from assignment2 import phi

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
omegaGrid = 2*np.pi*50


Xl_cable = L_cable*omegaGrid #Xl for cable connection
Xc_cable = -1/(omegaGrid*C_cable) # Xc for cable connection

# Transformer
Vlow = 690 #unit V
Vhigh = 33*10**3 #unit V
N_ratio = Vlow/Vhigh # (N1/N2) no unit

Ic = Vhigh/np.sqrt(3)
P_poc = 3*Vg*Ig - 3*(R1+R2+R_cable)*Ig**2
Q_poc = -omegaEl*(L1+L2+L_cable)*Ig**2+Ic**2/(1/(omegaEl*C_cable))

PF = np.divide(P2,S2) # Power factor

plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window
plt.plot(RPM,abs(P_poc))
plt.plot(RPM,abs(Q_poc))
plt.xlabel('n (RPM)')
plt.ylabel('Power (W)')
plt.legend(['Active power','Reactive power'])