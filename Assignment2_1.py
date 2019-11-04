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
from Assignment_2_Q2 import Vg,Ig,omegaEl,RPM,PlossPhase,Pmech
from assignment2 import phi

plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window

L1 = 2*10**-6 #unit H (Transformer)
R1 = 2*10**-3 #unit Ohm (Transformer)
L2 = 2*10**-6 #unit H (Transformer)
R2 = 2*10**-3 #unit Ohm (Transformer)
# Transformer
Vlow = 690 #unit V
Vhigh = 33*10**3 #unit V
N_ratio = Vlow/Vhigh # (N1/N2) no unit

R_cable = 1*N_ratio**2 #Unit Ohm
L_cable = 5*10**-3*N_ratio**2 #Unit H
C_cable = 1*10**-6 #Unit F

V_rms = Vhigh/np.sqrt(3)
V_poc = V_rms*np.cos(phi)+1j*V_rms*np.sin(phi)


Ic = (V_poc)/(-1/(omegaEl*C_cable))
#Vt = (-1j*omegaEl*(L1+L2+L_cable)+(R1+R2+R_cable))*Ig
#Pg = 3*(Pmech - PlossPhase)
S_poc = 3*(V_poc)*np.conj(Ig*N_ratio + 1j*Ic)
#P_poc = Pg - 3*(R1+R2+R_cable)*Ig**2
#Q_poc = 3*(-omegaEl*(L1+L2+L_cable)*Ig**2 + Ic**2*(1/(omegaEl*C_cable)))
P_poc = S_poc.real
Q_poc = S_poc.imag
PhaseShift = np.arctan(Q_poc/P_poc)

PF = P_poc/np.sqrt(P_poc**2+Q_poc**2) # Power factor

Efficiency = P_poc/(3*Pmech)

plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window

plt.figure()
plt.plot(RPM,abs(P_poc))
plt.plot(RPM,abs(Q_poc))
plt.xlabel('n (RPM)')
plt.ylabel('Power (W)')
plt.legend(['Active power','Reactive power'])