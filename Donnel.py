import numpy as np
import matplotlib.pyplot as plt

# define parameters of the shell
# density
rho_s = 7850
# elasticity modulus
E = 200e9
# Poission's ration
sigma = 0.3
# radius of the neutral surface
a = 0.5
# thickness
h= 0.05
# longitudinal wave velocity
cp = np.sqrt(E/(rho_s*(1-sigma**2)))
# thickness factor
beta = np.sqrt(h**2/(12*a**2))

# define parameters of EM excitations
# amplitude
F0 = 1000
# order
n = 0
# frequency
f = 50

# plot mode
mode = 1 # Fix order and change frequency
# mode = 2 # Fix frequency and change order

if (mode == 1):
  f = np.arange(500, 5000, 100)
  omega = 2. * np.pi * f
  
elif (mode == 2):
  n = np.arange(0, 11, 2) 
  
omega = 2. * np.pi * f
Omega = omega*a/cp
Omega1 = np.sqrt( 0.5*(1+n**2+beta**2*n**4+np.sqrt((1+n**2+beta**2*n**4)**2-4*beta**2*n**6)) )
Omega2 = np.sqrt( 0.5*(1+n**2+beta**2*n**4-np.sqrt((1+n**2+beta**2*n**4)**2-4*beta**2*n**6)) )
w = -(a**2*(1-sigma**2)) / (E*h) * (F0*(Omega**2 - n**2)) / ( (Omega**2 - Omega1**2)*(Omega**2 - Omega2**2) )


fig, ax = plt.subplots()
if (mode==1):
  ax.plot(f, np.abs(w))
elif (mode==2):
  ax.plot(n, np.abs(w))
ax.set_yscale('log')
plt.show()


#print(w)
