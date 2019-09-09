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
a = 0; a_prime = 0; a_crit = 1/3

tol = 1E-6
count = 0
for i in range(len(r)):
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
            anew = 1/(4*F*(np.sin(phi)**2)/(sigma*Cn)-1)
        else:
            K = 4*F*(np.sin(phi)**2)/(sigma*Cn)
            anew = 1/2*(2+K*(1-2*a_crit)-np.sqrt((K*(1-2*a_crit)+2)**2+4*(K*a_crit**2-1)))
        a_prime = 1/(4*F*(np.sin(phi)*np.cos(phi))/(sigma*Ct)-1)
        if abs(anew-a) <= tol:
            a = anew
            break
        else:
            a = anew 

        


