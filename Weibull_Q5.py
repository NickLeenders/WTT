# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 15:40:03 2019

@author: Jan Koene
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats as stats


# Weibull distribution parameters
A = 9
k = 1.9
Vin = 4
Vout = 25

V0 = np.linspace(Vin,Vout,50)

# Supply array of Pi values for each Vi in V0
P = np.zeros(len(V0))+5000

# =============================================================================
# Comparing given formula to scipy's weibull_min
# f01 = np.exp(-(Vin/A)**k)-np.exp(-(Vout/A)**k)
# 
# # weibull pdf and cdf
# #f_weib = stats.weibull_min.pdf(V0, c=k, scale=A)
# #F_weib = stats.weibull_min.cdf(V0, c=k, scale=A)
# #plt.plot(V0, F_weib)
# 
# # Probability Vin<V0<Vout
# Po = 1 - stats.weibull_min.cdf(Vout, c=k, scale=A)
# Pu = stats.weibull_min.cdf(Vin, c=k, scale=A)
# Ptotal = 1 -(Po + Pu)
# =============================================================================

fi = np.zeros(len(V0)-1)

for i in range(len(V0)-1):
    fi[i] = 0.5*(P[i]+P[i+1])*8760*(np.exp(-(V0[i]/A)**k)-np.exp(-(V0[i+1]/A)**k))
    
AEO = sum(fi)/10**6
print('GWh = ' + str(AEO))

# Changing Vout to 20m/s
fi_2 = np.copy(fi)
Vout_2 = 20
V0_2 = np.linspace(Vin,Vout_2,50)
for i in range(len(V0_2)-1):
    fi[i] = 0.5*(P[i]+P[i+1])*8760*(np.exp(-(V0_2[i]/A)**k)-np.exp(-(V0_2[i+1]/A)**k))
    
AEO_2 = sum(fi)/10**6
print('Energy lost in GWh = ' + str(AEO - AEO_2))