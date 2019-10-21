# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 11:57:25 2019

@author: dreth
"""

import numpy as np
from matplotlib import pyplot as plt
import cmath

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
P_loss = 3*R_s*abs(I_gen)**2
