import numpy as np
import matplotlib.pyplot as plt


# Define parameters
# Density
rho = 1.29 # kg/m^3
# Sound speed of air
c = 331.4 # m/s
# Angular frequency
omega = 2*np.pi*1000
# Wave number
k = omega/c
# Volumn velocity
Q = 0.005

azimuths = np.radians(np.linspace(0, 360, 72))
zeniths = np.linspace(0.001, 1, 500)
r, theta = np.meshgrid(zeniths, azimuths)

# Acoustic pressure
P = -1j * omega * rho * Q * np.exp(1j * k * r)/(4 * np.pi * r)
v = np.amax([np.amax(P.real), np.amin(P.real)])
print(v)
# Acoustic velocity
Vr = (Q * k)/(1j * 4 * np.pi * r) * (1- 1/(1j *k * r)) * np.exp(1j * k * r)

fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
Pre = plt.contourf(theta, r, P.real, levels=100, vmax=v, vmin=-v, cmap='bwr')
#Vel = plt.contourf(theta, r, Vr.real, levels=100, cmap='hsv')
plt.colorbar()



plt.show()