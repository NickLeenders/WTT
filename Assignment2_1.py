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

# Circuit 1
Vpoc = 33*10**3*np.sqrt(3)
omegaGrid = 2*np.pi*50
Z1 = R1 + L1*1j*omegaGrid
Z2 = R1 + L1*1j*omegaGrid
Z_circuit1 = Z1+Z2+R_cable+L_cable*omegaGrid*1j
Ipoc = Ig

# Cable
Ic = Vpoc/(1j/(omegaGrid*C_cable))+Ipoc
Vc = Ic*(R_cable+1j*omegaGrid*L_cable)+Vpoc