import numpy as np

r = np.array([24.5])
R = 31
rho = 1.225
V0 = 8.0
omega = 2.61
theta_p = np.deg2rad(-3.0) #in degrees
beta = np.deg2rad(2.0)
c = 0.5
Cl = 0.5
Cd = 0.01
B =  3

#Initialize alpha
a = 0; a_prime = 0;

tol = 1E-6
for i in range(len(r)):
    while (abs(an-a)>tol):
        phi = np.arctan(((1-a)*V0)/((1+a_prime)*omega*r[i]))
        alpha = phi-(beta+theta_p)
        
        #Table lookup for Cl, Cd
        
        Cn = Cl*np.cos(phi) + Cd*np.sin(phi)
        #F = (2/np.pi*np.arccos((np.exp((-B/2)*(R-r[i])/(r[i]*np.sin(abs(phi))))))
        
        
        F = 0.981
        sigma = c*B/(2*np.pi*r)
        a = 1/(4*F*(np.sin(phi)**2)/(sigma*Cn)-1)
        
        
        
        if a <= 1/3:
            anew = 1/(4*F*(np.sin(phi)**2)/(sigma*Cn)-1)
        else:
            CT = ((1-a)**2)*sigma*Cn/(np.sin(phi)**2)
        
        
        


