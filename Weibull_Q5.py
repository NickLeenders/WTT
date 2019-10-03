# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 15:40:03 2019

@author: I. Koune
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats as stats
from scipy.interpolate import interp2d

#%%Parameters
R = 89.17
A = np.pi*R**2
rho = 1.225


#V0 = 9.0
V0 = np.linspace(4,11.19,100)
#omega = np.array([1.01])
#omega = np.linspace(0.5,1.01,20)
#TSR = omega*R/V0
omega = 8.0/R*V0
theta_p = np.deg2rad(np.array([0]))
#theta_p = np.deg2rad(np.linspace(-2,5,15)) 
B =  3

# Export figures as pdf
saveFig = 1

#%% Reference wind turbine blade description
r = np.zeros([18])
r[:] = np.array([2.80, 11.00, 16.87, 22.96, 32.31, 41.57, 50.41, 
                     58.53, 65.75, 71.97, 77.19, 78.71, 80.14, 82.71, 
                     84.93, 86.83, 88.45, 89.17-(0.2)])
r_int = np.append(r[:],R) #For tip integration
c = np.array([5.38, 5.45, 5.87, 6.18, 6.02, 5.42, 4.70, 4.00, 3.40, 
                2.91, 2.54, 2.43, 2.33, 2.13, 1.90, 1.63, 1.18, 0.60])
beta = np.deg2rad(np.array([14.50, 14.43, 12.55, 8.89, 6.38, 4.67, 
                              2.89, 1.21, -0.13, -1.11, -1.86, -2.08, 
                              -2.28, -2.64, -2.95, -3.18, -3.36, 
                              -3.43]))
ratio = np.array([100.00, 86.05, 61.10, 43.04, 32.42, 27.81,
                    25.32, 24.26, 24.10, 24.10, 24.10, 24.10, 24.10,
                    24.10, 24.10, 24.10, 24.10, 24.10])

#Import of data
airFoilThickness = np.array([24.1, 30.1, 36.0, 48.0, 60.0, 100])
suffixData = np.array(['FFA-W3-241.txt', 'FFA-W3-301.txt', 
                       'FFA-W3-360.txt', 'FFA-W3-480.txt', 
                       'FFA-W3-600.txt', 'cylinder.txt'])
sizeTest = np.genfromtxt('FFA-W3-241.txt',
                     dtype=None,
                     delimiter=None)
FFAdata = np.zeros((len(sizeTest),len(sizeTest[0]),len(suffixData)), 
                   dtype='float64')
for i in range(0,len(suffixData)):
    name = suffixData[i]
    FFAdata[:,:,i] = np.loadtxt(name,
                     dtype=None,
                     delimiter=None)

del sizeTest
del name

#Interpolation of data
intCl = interp2d(airFoilThickness,FFAdata[:,0,0],FFAdata[:,1,:], 
                 kind='linear') # Thickness, angle OUTPUT: Is in array
intCd = interp2d(airFoilThickness,FFAdata[:,0,0],FFAdata[:,2,:], 
                 kind='linear') # Thickness, angle
intCm = interp2d(airFoilThickness,FFAdata[:,0,0],FFAdata[:,3,:], 
                 kind='linear')

#%%Initialize variables
a_crit = 1/3
Pn = np.zeros([len(r)+1])
Pt = np.copy(Pn)
tol = 1E-5
CP = np.zeros([len(theta_p),len(omega)])
CT = np.copy(CP)
Power = np.zeros(len(omega))

"""BEGIN PITCH LOOP"""
for k in range(len(theta_p)):

    """BEGIN OMEGA LOOP"""
    for j in range(len(omega)):
    
        """BEGIN BEM LOOP"""
        for i in range(len(r)):
            a = 0; a_prime = 0;
            count = 0
        
            while True:
                count +=1
                phi = np.arctan(((1-a)*V0[j])/((1+a_prime)*omega[j]*r[i]))
                alpha = phi-(beta[i]+theta_p[k])
                #Table lookup for Cl, Cd
                Cl = intCl(ratio[i],np.rad2deg(alpha))
                Cd = intCd(ratio[i],np.rad2deg(alpha))
                Cn = Cl*np.cos(phi) + Cd*np.sin(phi)
                Ct = Cl*np.sin(phi) - Cd*np.cos(phi)
                #Prandt's Tip Loss
                F = 2/np.pi*np.arccos(np.exp((-B/2)*(R-r[i])/(r[i]*np.sin(abs(phi)))))
                
                sigma = c[i]*B/(2*np.pi*r[i])       #Solidity
                #Glauert Correction
                if a <= a_crit:
                    anew = 1/(4*F*(np.sin(phi)**2)/(sigma*Cn)+1)
                else:
                    CT_glauert = (1-a)**2*Cn*sigma/np.sin(phi)**2
                    astar = CT_glauert/(4*F*(1-0.25*(5-3*a)*a))
                    anew = 0.1*astar+0.9*a
                a_prime = 1/(4*F*(np.sin(phi)*np.cos(phi))/(sigma*Ct)-1)
                if abs(anew-a) <= tol:
                    a = anew
                    break
                elif count == 200:
                    #print('count > 200')
                    break
                else:
                    a = anew
            #print('Count = ' + str(count))
            Vrel = V0[j]*(1-a)/np.sin(phi)
            Fl = 1/2*rho*Vrel**2*c[i]*Cl
            Fd = 1/2*rho*Vrel**2*c[i]*Cd
            Pn[i] = Fl*np.cos(phi)+Fd*np.sin(phi)
            Pt[i] = Fl*np.sin(phi)-Fd*np.cos(phi)
        """END BEM LOOP"""
        
        #Power, Thrust and coefficients
        Power[j] = omega[j]*B*np.trapz(Pt*r_int, r_int)        #Rotor mechanical power
        CP[k,j] = Power[j]/(0.5*rho*V0[j]**3*A)
        Thrust = B*np.trapz(Pn,r_int)                       #Rotor thrust
        CT[k,j] = Thrust/(0.5*rho*V0[j]**2*A)
    """END OMEGA LOOP"""
"""END PITCH LOOP"""

plt.figure('Test', figsize = (5,4))
plt.plot(V0, Power)
#%% Start Question 5

# Create an array of power values corresponding to V0
Vin = 4
Vout =  25
Vout_2 = 20

nAppend = 500
Append_list = np.linspace(11.19, Vout, nAppend)
for item in Append_list:
    print(item)
    Power = np.append(Power, Power[-1])
    V0 = np.append(V0, item)


#%% Weibull distribution parameters
A = 9
k = 1.9

idx = (np.abs(V0-Vout_2)).argmin()
V0_2 = V0[0:idx+1]
#Power_2 = Power[0:idx+1]
    

# Supply array of Pi values for each Vi in V0

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
#
# Same results, either one can be used
# =============================================================================


fi = np.zeros(len(V0)-1)
#fi_2 = np.zeros(len(V0_2)-1)
fi_2 = np.zeros(idx)
fi_check = np.zeros(len(V0)-1)

for i in range(len(V0)-1):
    fi[i] = 0.5*(Power[i]+Power[i+1])*8760*(np.exp(-(V0[i]/A)**k)-np.exp(-(V0[i+1]/A)**k))

# Total annual estimated production   
AEO = sum(fi)/10**9

# When Vout = 20 m/s
for i in range(len(V0_2)-1):
    fi_2[i] = 0.5*(Power[i]+Power[i+1])*8760*(np.exp(-(V0_2[i]/A)**k)-np.exp(-(V0_2[i+1]/A)**k))    
AEO_2 = sum(fi_2)/10**9

# AEO - AEO2 should converge to this
test_diff = Power[-1]*8760*(np.exp(-(20/A)**k)-np.exp(-(25/A)**k))/10**9

# Plot pdf for given distribution
f_weib_V0 = stats.weibull_min.pdf(V0, c=k, scale=A)
f_weib_V2 = stats.weibull_min.pdf(V0_2, c=k, scale=A)

# PDF plot
plt.figure('Weibull PDF', figsize = (5,4))
plt.plot(V0, f_weib_V0, 'red', label = 'Weibull PDF, $V_{out}$ = 25 m/s')
plt.plot(V0_2, f_weib_V2, 'blue', label = 'Weibull PDF, $V_{out}$ = 20 m/s')
plt.grid(c='k', alpha = .3)
plt.xlabel('Wind speed $V_{0}$ [m/s]', fontsize = 14)
plt.ylabel('Probability density', fontsize = 14)
plt.fill_between(V0, f_weib_V0, color = 'red')
plt.fill_between(V0_2, f_weib_V2, color = 'blue')
if saveFig:
    plt.savefig('Weibull_PDF.pdf',bbox_inches='tight')
    
# CDF Plot
F_weib = stats.weibull_min.cdf(V0, c=k, scale=A)
F_weib_2 = stats.weibull_min.cdf(V0_2, c=k, scale=A)
plt.figure('Weibull CDF', figsize = (5,4))
plt.plot(V0, F_weib, 'red', label = 'Weibull CDF, $V_{out}$ = 25 m/s')
plt.plot(V0_2, F_weib_2, 'blue', label = 'Weibull CDF, $V_{out}$ = 20 m/s')
plt.grid(c='k', alpha = .3)
plt.xlabel('Wind speed $V_{0}$ [m/s]', fontsize = 14)
plt.ylabel('Cumulative probability', fontsize = 14)
if saveFig:
    plt.savefig('Weibull_CDF.pdf',bbox_inches='tight')

# Probability of exceedance of 20 m/s windspeed
F_20 = stats.weibull_min.cdf(20, c=k, scale=A)

# Percent decrease in AEP
Perc_decrease = (AEO-AEO_2)/AEO*100
    
# Money loss assuming 0.77 kr./kWh
0.77*(10**6)*(AEO-AEO_2)