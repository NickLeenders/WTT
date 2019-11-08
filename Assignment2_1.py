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

Pg = Vg*Ig

plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window

L1 = 2*10**-6 #unit H (Transformer)
R1 = 2*10**-3 #unit Ohm (Transformer)
L2 = 2*10**-6 #unit H (Transformer)
R2 = 2*10**-3 #unit Ohm (Transformer)

omegaGrid = 2*np.pi*50
# Transformer
Vlow = 690 #unit V
Vhigh = 33*10**3 #unit V
N_ratio = Vlow/Vhigh # (N1/N2) no unit
Lcable = 10 #Unit km
R_cable = 1 #Unit Ohm
L_cable = 5*10**-3 #Unit H
C_cable = 1*10**-6 #Unit F

V_rms = Vhigh/np.sqrt(3)

Z1 = R1 + L1*1j*omegaGrid
Z2 = R1 + L1*1j*omegaGrid
Zc = R_cable + 1j*omegaGrid*L_cable
Ztotal = Z1/N_ratio**2+Z2/N_ratio**2+Zc
Zcap = -1j/(omegaGrid*C_cable)

I1 = Ig*N_ratio
V2 = Vg/N_ratio
V_poc = -I1*Ztotal+V2



Ic = -(I1*Ztotal-V2)/(Zcap)
Igrid = (I1*Zcap+I1*Ztotal-V2)/(Zcap)

#Vt = (-1j*omegaEl*(L1+L2+L_cable)+(R1+R2+R_cable))*Ig
#Pg = 3*(Pmech - PlossPhase)
S_poc = 3*(V_poc)*np.conj(Igrid + Ic)
#P_poc = Pg - 3*(R1+R2+R_cable)*Ig**2
#Q_poc = 3*(-omegaEl*(L1+L2+L_cable)*Ig**2 + Ic**2*(1/(omegaEl*C_cable)))
P_poc = S_poc.real
Q_poc = S_poc.imag
PhaseShift = np.arctan(Q_poc/P_poc)

PF = P_poc/np.sqrt(P_poc**2+Q_poc**2) # Power factor

Efficiency = P_poc/(3*Pmech)*100

plt.cla()   # Clear axis
plt.clf()   # Clear figure
plt.close() # Close a figure window

plt.figure()
P_poc_plot = P_poc/10**6
Q_poc_plot = Q_poc/10**6
plt.plot(RPM,abs(P_poc_plot))
plt.plot(RPM,abs(Q_poc_plot))
plt.grid()
plt.xlabel('n (RPM)')
plt.ylabel('Power (MW or MVAr)')
plt.legend(['Active power (MW)','Reactive power (MVAr)'])
plt.savefig('S_figure.pdf')

plt.figure()
plt.plot(RPM,Efficiency)
plt.grid()
plt.xlabel('n (RPM)')
plt.ylabel('Efficiency (%)')
plt.savefig('Efficiency.pdf')