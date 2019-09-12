import numpy as np
import matplotlib.pyplot as plt
from main_assignment1 import Values1

r = Values1.blade_table['r']
c = Values1.blade_table['c']

Fd = np.copy(r)
Fl = np.copy(r)
Pn = np.copy(r)
Pt = np.copy(r)
R = 89.17
rho = 1.225
V0 = 8.0
omega = 2.61
theta_p = np.deg2rad(-3.0) 
beta = Values1.blade_table['beta']/(180/np.pi)
B =  3
tc = Values1.blade_table['t/c']
# Export figures as pdf
saveFig = 0
cl_list = []
cd_list = []
alpha_list = []
phi_list = []
a_list = []
#Initialize alpha
a_crit = 1/3

tol = 1E-6
for i in range(len(r)):
    a = 0; a_prime = 0;
    count = 0

    while True:
        count +=1
        phi = np.arctan(((1-a)*V0)/((1+a_prime)*omega*r.values[i]))
        alpha = phi-(beta.values[i]+theta_p)
        Cl, Cd = Values1.interpolate_tc(tc.values[i], alpha)
        #Table lookup for Cl, Cd
        Cn = Cl*np.cos(phi) + Cd*np.sin(phi)
        Ct = Cl*np.sin(phi) - Cd*np.cos(phi)
        
        F = 2/np.pi*np.arccos(np.exp((-B/2)*(R-r.values[i])/(r.values[i]*np.sin(abs(phi)))))
        #F = 0.981
        
        sigma = c.values[i]*B/(2*np.pi*r.values[i])
        #a = 1/(4*F*(np.sin(phi)**2)/(sigma*Cn)-1)
           
        if a <= a_crit:
            anew = 1/(4*F*(np.sin(phi)**2)/(sigma*Cn)+1)
        else:
            CT = (1-a)**2*Cn*sigma/np.sin(phi)**2
            astar = CT/(4*F*(1-0.25*(5-3*a)*a))
            anew = 0.1*astar+0.9*a
        a_prime = 1/(4*F*(np.sin(phi)*np.cos(phi))/(sigma*Ct)-1)
        if abs(anew-a) <= tol:
            a = anew
            break
        elif count == 200:
            print('count > 200')
            break
        else:
            a = anew
    print('Count = ' + str(count))
    cl_list.append(Cl)
    cd_list.append(Cd)
    alpha_list.append(alpha*(180/np.pi))
    phi_list.append(phi*(180/np.pi))
    a_list.append(a)
    Vrel = V0*(1-a)/np.sin(phi)
    Fl = 1/2*rho*Vrel**2*c.values[i]*Cl
    Fd = 1/2*rho*Vrel**2*c.values[i]*Cd
    Pn[i] = Fl*np.cos(phi)+Fd*np.sin(phi)
    Pt[i] = Fl*np.sin(phi)-Fd*np.cos(phi)

# Rotor torque along the shaft
MT_blade = np.trapz(Pt*r, r)
MT_tot = B*MT_blade


plt.figure('Tip speed ratio',figsize=(5,4))
plt.plot(r, Pt, 'xkcd:amber',
         label = 'Tangential force distribution')
plt.grid(c='k', alpha=.3)
plt.xlabel('Radius [m]', fontsize=14)
plt.ylabel('Pt [N]', fontsize=14)
plt.tick_params(labelsize=12)
plt.legend(fontsize = 12)
if saveFig:
    plt.savefig('tipSpeedRatio.pdf',bbox_inches='tight')
   
# Total power and power coefficient
Ptot = MT_tot*omega
CP = Ptot/(0.5*rho*V0**3*np.pi*R**2)
pass
