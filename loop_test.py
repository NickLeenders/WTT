import numpy as np
import matplotlib.pyplot as plt
from .user import User
from .dir import Dir

r = np.linspace(1, 25, 40)
Fd = np.copy(r)
Fl = np.copy(r)
Pn = np.copy(r)
Pt = np.copy(r)
R = 31
rho = 1.225
V0 = 8.0
omega = 2.61
theta_p = np.deg2rad(-3.0) 
beta = np.deg2rad(2.0)
c = 1.5
Cl = 0.5
Cd = 0.01
B =  3

# Export figures as pdf
saveFig = 0

#Initialize alpha
a_crit = 1/3

tol = 1E-6
for i in range(len(r)):
    a = 0; a_prime = 0;
    count = 0
    while True:
        count +=1
        phi = np.arctan(((1-a)*V0)/((1+a_prime)*omega*r[i]))
        alpha = phi-(beta+theta_p)
        
        #Table lookup for Cl, Cd
        Cn = Cl*np.cos(phi) + Cd*np.sin(phi)
        Ct = Cl*np.sin(phi) - Cd*np.cos(phi)
        
        F = 2/np.pi*np.arccos(np.exp((-B/2)*(R-r[i])/(r[i]*np.sin(abs(phi)))))
        #F = 0.981
        
        sigma = c*B/(2*np.pi*r[i])
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

    Vrel = V0*(1-a)/np.sin(phi)
    Fl = 1/2*rho*Vrel**2*c*Cl
    Fd = 1/2*rho*Vrel**2*c*Cd
    Pn[i] = Fl*np.cos(phi)+Fd*np.sin(phi)
    Pt[i] = Fl*np.sin(phi)-Fd*np.cos(phi)

# Rotor torque along the shaft
MT_blade = np.trapz(Pt*r, r)
MT_tot = B*MT_blade


plt.figure('Tip speed ratio',figsize=(5,4))
plt.plot(r, Pt, 'xkcd:amber',
         label = 'Tip speed ratio, $\lambda$')
plt.grid(c='k', alpha=.3)
plt.xlabel('Time [$s$]', fontsize=14)
plt.ylabel('$\lambda$ [-]', fontsize=14)
plt.tick_params(labelsize=12)
plt.legend(fontsize = 12)
if saveFig:
    plt.savefig('tipSpeedRatio.pdf',bbox_inches='tight')
   
# Total power and power coefficient
Ptot = MT_tot*omega
CP = Ptot/(0.5*rho*V0**3*np.pi*R**2)
