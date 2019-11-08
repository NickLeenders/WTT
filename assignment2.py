# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 11:57:25 2019

@author: dreth
"""

import numpy as np
from matplotlib import pyplot as plt
import cmath

saveFig = 0

#%%

n_poles = 80                            #no. of generator poles
k_e = 5.34                              #From previous calculations
k_mech = 1.333*1e6                      #From previous calculations
n_rpm = np.linspace(3,15,13)            #rpm range
omega_mech = n_rpm*2*np.pi/60           #Mechanical angular velocity
omega_e = omega_mech*n_poles/2          #Electrical angular velocity
R_s = 4*1e-3
L_s = 1.2*1e-3

##Impedance
#Z_s = R_s + 1j*omega_e*L_s
##EMF voltage (rms)
#V_EMF = k_e*omega_e
#P_EMF = k_mech*omega_mech**3/3
#I_gen = P_EMF/V_EMF
#P_loss = 3*R_s*abs(I_gen)**2        #For all three phases combined

grid_omega = 50*2*np.pi
ratio = 690/(33*1e3)
R1 = 2*1e-3
R2t = R1
L1 = 2*1e-6
L2t = L1
Rc = 10*0.1
Lc = 10*0.5*1e-3
Cc = 10*0.1*1e-3
Zc = Rc #+ 1j*grid_omega*Lc
Z1 = R1 + 1j*grid_omega*L1
Z2t = R2t + 1j*grid_omega*L2t
Zct = Zc*ratio**2
Zt = Z1 + Z2t + Zct
Zt2 = Zt/ratio**2
phi = np.arctan(np.imag(Zt)/np.real(Zt))
phid = np.rad2deg(phi)

 #%%  plot acceleration at tip of blade 1
plt.figure('Generator power loss',figsize=(5,4))
#plt.plot(n_rpm,P_loss,'xkcd:amber',label = 'Generator power loss')
plt.grid(c='k', alpha=.3)
plt.xlabel('RPM [-]', fontsize=14)
plt.ylabel('Power loss [W]', fontsize=14)
plt.tick_params(labelsize=12)
plt.legend(fontsize = 12)
if saveFig:
    plt.savefig('accelerationBlade1Tip_y.pdf',bbox_inches='tight')
