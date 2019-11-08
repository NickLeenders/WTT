
"""
Created on Sun Nov  3 11:36:30 2019

@author: nickleenders
"""
from scipy.optimize import fsolve
import numpy as np
import cmath
import math
import matplotlib.pyplot as plt
from sympy import *
from sympy.solvers.solveset import linsolve
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
Vgrid = 33*10**3
Vgridlow = 690/np.sqrt(3)

V_rms = Vhigh/np.sqrt(3)
A = []
b = []
sol = []
I2 = Ig*N_ratio
V2 = Vg/N_ratio
Z1 = R1 + L1*1j*omegaGrid
Z2 = R1 + L1*1j*omegaGrid
Zc = R_cable + 1j*omegaGrid*L_cable
Ztotal = Z1/N_ratio**2+Z2/N_ratio**2+Zc
Zcap = 1j/(omegaGrid*C_cable)
#for i,text in enumerate(I2):
#    A.append(np.array([[1, -1, 1],
#        [0, -Zc, 0],
#        [0, 0, -1]]))
#    b.append([[I2[i]],[I2[i]*Ztotal],[I2[i]*Ztotal+V2[i]]])
#    sol.append(np.linalg.solve(A,b).squeeze())

Ic = Vgrid /(Zcap)
Vc = Vgrid

# knowns Ztotal, Zc, Vgrid, Pg, 
PL, Pgr, Qc1, I11, I211 = symbols('PL Pgr Qc1 I11 I211')
PL      = 3*(R1+R2+R_cable)*(I11)**2
Pgr     = 3*Vlow*I211
Qc1     = 3*omegaGrid*(L1+L2+L_cable)*(I11)**2-3*omegaGrid*C_cable/N_ratio**2*Vgridlow**2

eqn2    = -Pmech[-1] + (PL**2+Pgr**2+Qc1**2+2*PL*Pgr)**0.5
eqn1 = I211 - I11 + 1j*omegaGrid*C_cable/N_ratio**2 *Vgridlow

sol = solve([eqn1,eqn2],[I11,I211])

#sol = solve([eq1,eq2,eq3,eq4,eq5],[I1, Igrid, V1])
# IC and VC stay the same so capacitor delivers a constant reactive power, the reactive power is going to differ because of the variable reactive power by 
